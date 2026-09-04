# =============================================================================
# CAMEROON COMPUTABLE GENERAL EQUILIBRIUM (CGE) MODEL
# =============================================================================
# A static, single-country CGE model calibrated to Cameroon's 1979-80 base year.
# Follows the Davis-De Melo-Robinson framework (late 1970s), originally coded
# in GAMS by Condon, Dahl, and Devarajan (World Bank, 1987).
#
# Author : John-Mary Matovu
#          Research Fellow, Research for Transformation and Development (RTD)
#          Kampala, Uganda  |  www.rtdug.org
#
# Dependencies (all free/open-source):
#   XLSX.jl   — reads Excel calibration data
#   JuMP.jl   — algebraic modelling language for optimisation
#   Ipopt.jl  — interior-point NLP solver (COIN-OR)
#
# Usage (see README.md):
#   using Pkg; Pkg.activate(@__DIR__); using CGECameroon
#   baseline, sim1 = CGECameroon.example()
#
# This file was refactored (eps-integration branch) from the original top-to-bottom
# script into a module with no include-time solves and no module-level mutable
# parameter state, so it can be embedded safely in a long-lived, multi-request host
# process (e.g. the Economic Policy Simulator web app). Every equation is
# mathematically identical to the original script -- only *where the numbers live*
# changed: `load_data` reads the workbook into a `RawData`, `calibrate` turns that
# into a `Params` (a plain, deep-copyable struct holding every calibrated parameter
# and closure value that used to be a module-level `Dict` or scalar -- including the
# household saving rate, promoted out of a bare literal into `Params.mps0`, and the
# exchange rate `Params.er`), and `solve` builds a fresh JuMP model from a `Params`
# argument (never from a global). `with_shocks` returns a shocked deep copy of a
# `Params`, so running two scenarios back to back never leaks state between them.
# =============================================================================

module CGECameroon

using XLSX
using JuMP
using Ipopt

export RawData, Params, Simulation, load_data, calibrate, solve, with_shocks, example

# =============================================================================
# SETS
# =============================================================================
# SEC  — all 11 production sectors
# IT   — 9 traded sectors (subject to Armington imports and CET exports)
# ITN  — 2 non-traded sectors (domestic supply = domestic demand, no trade)
# LC   — 3 labour categories
#
# These are structural (today: Cameroon-only, hardcoded) rather than scenario
# parameters, so they stay as plain constants -- but every function below reads them
# off `RawData`/`Params`, never off these module globals directly, so a future
# multi-country loader only has to change `load_data`.
# =============================================================================

const SEC = [
    :agsubsist,   # Food crops (subsistence agriculture)
    :agexpind,    # Cash crops (export-oriented agriculture)
    :sylvicult,   # Forestry and logging
    :indalim,     # Food processing industry
    :bienscons,   # Consumer goods manufacturing
    :biensint,    # Intermediate goods manufacturing
    :cimint,      # Construction materials
    :bienscap,    # Capital goods manufacturing
    :construct,   # Construction (non-traded)
    :services,    # Private services
    :publiques    # Public services (non-traded)
]

const IT  = [:agsubsist, :agexpind, :sylvicult, :indalim, :bienscons,
             :biensint, :cimint, :bienscap, :services]
const ITN = [:construct, :publiques]

const LC  = [:rural, :urbanunsk, :urbanskil]  # Rural, Urban unskilled, Urban skilled

# =============================================================================
# RAW DATA (returned by `load_data`)
# =============================================================================

"""
    RawData

The five `camdata.xlsx` sheets, read exactly as the original script did (same
hardcoded cell ranges and row order -- see `load_data`), plus the sector/labour sets.
No calibration algebra has run yet; this is pass-through data.
"""
struct RawData
    SEC::Vector{Symbol}
    IT::Vector{Symbol}
    ITN::Vector{Symbol}
    LC::Vector{Symbol}
    io::Dict{Tuple{Symbol,Symbol},Float64}     # iotable sheet
    imat::Dict{Tuple{Symbol,Symbol},Float64}   # imat sheet
    wdist::Dict{Tuple{Symbol,Symbol},Float64}  # wagedist sheet
    xle::Dict{Tuple{Symbol,Symbol},Float64}    # employment sheet
    zz::Dict{Tuple{Symbol,Symbol},Float64}     # miscellaneous sheet (17 rows x 11 sectors)
end

"""
    load_data(path::AbstractString) -> RawData

Reads `path` (a `camdata.xlsx`-layout workbook: sheets `iotable`, `imat`, `wagedist`,
`employment`, `miscellaneous`) into a `RawData`. Cell ranges and the `miscellaneous`
row order are hardcoded, exactly as in the original script -- there is still no
header-based lookup or shape validation.
"""
function load_data(path::AbstractString)
    dataio = XLSX.readdata(path, "iotable", "B3:L13")
    io = Dict{Tuple{Symbol,Symbol},Float64}()
    for i in 1:length(SEC), j in 1:length(SEC)
        io[SEC[i], SEC[j]] = dataio[i, j]
    end

    datacap = XLSX.readdata(path, "imat", "B3:L13")
    imat = Dict{Tuple{Symbol,Symbol},Float64}()
    for i in 1:length(SEC), j in 1:length(SEC)
        imat[SEC[i], SEC[j]] = datacap[i, j]
    end

    dataw = XLSX.readdata(path, "wagedist", "B3:D13")
    wdist = Dict{Tuple{Symbol,Symbol},Float64}()
    for i in 1:length(SEC), j in 1:length(LC)
        wdist[SEC[i], LC[j]] = dataw[i, j]
    end

    datae = XLSX.readdata(path, "employment", "B3:D13")
    xle = Dict{Tuple{Symbol,Symbol},Float64}()
    for i in 1:length(SEC), j in 1:length(LC)
        xle[SEC[i], LC[j]] = datae[i, j]
    end

    datazz = XLSX.readdata(path, "miscellaneous", "B3:L19")
    rowz = [:m0, :e0, :xd0, :k, :depr, :rhoc, :rhot, :eta, :pd0, :tm0,
            :itax, :cles, :gles, :kio, :dstr, :dst, :id]
    zz = Dict{Tuple{Symbol,Symbol},Float64}()
    for i in 1:length(rowz), j in 1:length(SEC)
        zz[rowz[i], SEC[j]] = datazz[i, j]
    end

    return RawData(SEC, IT, ITN, LC, io, imat, wdist, xle, zz)
end

# =============================================================================
# CALIBRATED PARAMETERS (returned by `calibrate`)
# =============================================================================

"""
    Params

Every calibrated parameter and closure value the model needs to build and solve --
what used to be ~30 module-level `Dict`s and scalars. Mutable so `with_shocks` can
`deepcopy` and edit a fresh instance per scenario; callers should otherwise treat a
`Params` as immutable (never mutate one that another in-flight `solve` might read).
"""
mutable struct Params
    SEC::Vector{Symbol}
    IT::Vector{Symbol}
    ITN::Vector{Symbol}
    LC::Vector{Symbol}

    # economy-wide scalars --------------------------------------------------------
    er::Float64        # exchange rate, CFAF/USD (was a bare module-level scalar)
    gr0::Float64        # base government revenue (start value only)
    gdtot0::Float64      # government real spending, closure-fixed
    cdtot0::Float64      # total private consumption at base (calibration only)
    fsav0::Float64       # foreign savings, closure-fixed
    mps0::Float64        # household saving rate, closure-fixed (was the literal 0.09305)
    y0::Float64          # base private GDP (start value only)

    # per-sector (Symbol => Float64) -----------------------------------------------
    depr::Dict{Symbol,Float64}
    rhoc::Dict{Symbol,Float64}
    rhot::Dict{Symbol,Float64}
    eta::Dict{Symbol,Float64}
    tm0::Dict{Symbol,Float64}
    te::Dict{Symbol,Float64}
    itax::Dict{Symbol,Float64}
    cles::Dict{Symbol,Float64}
    gles::Dict{Symbol,Float64}
    kio::Dict{Symbol,Float64}
    dstr::Dict{Symbol,Float64}
    m0::Dict{Symbol,Float64}
    e0::Dict{Symbol,Float64}
    xd0::Dict{Symbol,Float64}
    k0::Dict{Symbol,Float64}
    pd0::Dict{Symbol,Float64}
    pm0::Dict{Symbol,Float64}
    pe0::Dict{Symbol,Float64}
    pwm0::Dict{Symbol,Float64}
    pwe0::Dict{Symbol,Float64}
    pva0::Dict{Symbol,Float64}
    xxd0::Dict{Symbol,Float64}
    dst0::Dict{Symbol,Float64}
    id0::Dict{Symbol,Float64}
    cd0::Dict{Symbol,Float64}
    int0::Dict{Symbol,Float64}
    x0::Dict{Symbol,Float64}
    delta::Dict{Symbol,Float64}   # Armington share (IT only)
    ac::Dict{Symbol,Float64}      # Armington scale (IT only)
    gamma::Dict{Symbol,Float64}   # CET share (IT, plus 0 for any zero-import sector)
    at::Dict{Symbol,Float64}      # CET scale (IT only)
    ad::Dict{Symbol,Float64}      # Cobb-Douglas TFP

    # per-labour (Symbol => Float64) ------------------------------------------------
    wa0::Dict{Symbol,Float64}
    ls0::Dict{Symbol,Float64}

    # per (sector, sector) ------------------------------------------------------------
    io::Dict{Tuple{Symbol,Symbol},Float64}
    imat::Dict{Tuple{Symbol,Symbol},Float64}

    # per (sector, labour) or (labour, sector) -----------------------------------------
    wdist::Dict{Tuple{Symbol,Symbol},Float64}   # keyed (sector, labour)
    xle::Dict{Tuple{Symbol,Symbol},Float64}     # keyed (sector, labour)
    alphl::Dict{Tuple{Symbol,Symbol},Float64}   # keyed (labour, sector)
end

"""
    calibrate(raw::RawData) -> Params

Closed-form calibration algebra, identical to the original script's stage 3
(cge.jl:182-273 in the pre-refactor file): backs out share/shift parameters from the
base-year data so that every model equation holds exactly at the base-year point. No
optimisation happens here -- only the counterfactual `solve` calls use Ipopt. The
order of the blocks below matters (later blocks read earlier ones) and mirrors the
original exactly.
"""
function calibrate(raw::RawData)
    SEC, IT, ITN, LC = raw.SEC, raw.IT, raw.ITN, raw.LC
    io, imat, wdist, xle, zz = raw.io, raw.imat, raw.wdist, raw.xle, raw.zz

    # Base-year scalars that were set directly in code (not read from the workbook).
    wa0 = Dict{Symbol,Float64}(:rural => 0.11, :urbanunsk => 0.15678, :urbanskil => 1.8657)
    er     = 0.21     # Real exchange rate: CFAF per US dollar
    gr0    = 179.0    # Government revenue (billion CFAF)
    gdtot0 = 135.03   # Total government consumption (billion CFAF)
    cdtot0 = 947.98   # Total private consumption (billion CFAF)
    fsav0  = 36.841   # Foreign savings / current-account deficit (billion USD)
    mps0   = 0.09305  # Household saving rate (was a bare literal in `closuremp`)

    depr  = Dict{Symbol,Float64}(i => zz[:depr, i]      for i in SEC)
    rhoc  = Dict{Symbol,Float64}(i => 1/zz[:rhoc, i] - 1 for i in SEC)  # Armington exponent
    rhot  = Dict{Symbol,Float64}(i => 1/zz[:rhot, i] + 1 for i in SEC)  # CET exponent
    eta   = Dict{Symbol,Float64}(i => zz[:eta, i]        for i in SEC)
    tm0   = Dict{Symbol,Float64}(i => zz[:tm0, i]        for i in SEC)
    te    = Dict{Symbol,Float64}(i => 0.0                for i in SEC)  # No export duties at base year
    itax  = Dict{Symbol,Float64}(i => zz[:itax, i]       for i in SEC)
    cles  = Dict{Symbol,Float64}(i => zz[:cles, i]       for i in SEC)
    gles  = Dict{Symbol,Float64}(i => zz[:gles, i]       for i in SEC)
    kio   = Dict{Symbol,Float64}(i => zz[:kio, i]        for i in SEC)
    dstr  = Dict{Symbol,Float64}(i => zz[:dstr, i]       for i in SEC)

    # Replace zero-employment cells with 1 so 0^alpha doesn't produce NaN below.
    xllb = Dict{Tuple{Symbol,Symbol},Float64}((i, l) => xle[i, l] + (1 - sign(xle[i, l]))
                                              for i in SEC, l in LC)

    m0  = Dict{Symbol,Float64}(i => zz[:m0, i]  for i in SEC)
    e0  = Dict{Symbol,Float64}(i => zz[:e0, i]  for i in SEC)
    xd0 = Dict{Symbol,Float64}(i => zz[:xd0, i] for i in SEC)
    k0  = Dict{Symbol,Float64}(i => zz[:k, i]   for i in SEC)
    pd0 = Dict{Symbol,Float64}(i => zz[:pd0, i] for i in SEC)

    # At base year all domestic, import, and export prices are equal (normalisation).
    pm0 = copy(pd0)
    pe0 = copy(pd0)

    pwm0 = Dict{Symbol,Float64}(i => pm0[i] / ((1 + tm0[i]) * er) for i in SEC)
    pwe0 = Dict{Symbol,Float64}(i => pe0[i] / ((1 + te[i])  * er) for i in SEC)

    pva0 = Dict{Symbol,Float64}(i => pd0[i] - sum(io[j, i]*pd0[j] for j in SEC) - itax[i]
                                for i in SEC)

    xxd0 = Dict{Symbol,Float64}(i => xd0[i] - e0[i] for i in SEC)  # Domestic sales = output - exports
    dst0 = Dict{Symbol,Float64}(i => zz[:dst, i]    for i in SEC)
    id0  = Dict{Symbol,Float64}(i => zz[:id, i]     for i in SEC)
    ls0  = Dict{Symbol,Float64}(l => sum(xle[i, l] for i in SEC) for l in LC)
    y0   = sum(pva0[i]*xd0[i] - depr[i]*k0[i] for i in SEC)
    cd0  = Dict{Symbol,Float64}(i => cles[i]*cdtot0 for i in SEC)

    # Armington share parameter (delta): from first-order cost-minimisation condition.
    delta = Dict{Symbol,Float64}()
    for i in IT
        d = pm0[i]/pd0[i]*(m0[i]/xxd0[i])^(1 + rhoc[i])
        delta[i] = d / (1 + d)
    end

    # Composite good volumes at base (value terms for traded/non-traded sectors).
    x0 = Dict{Symbol,Float64}()
    for i in IT
        x0[i] = pd0[i]*xxd0[i] + pm0[i]*m0[i]
    end
    for i in ITN
        x0[i] = pd0[i]*xxd0[i]
    end

    # Armington shift parameter (ac): ensures armington() holds exactly at base data.
    ac = Dict{Symbol,Float64}(
        i => x0[i] / (delta[i]*m0[i]^(-rhoc[i]) + (1 - delta[i])*xxd0[i]^(-rhoc[i]))^(-1/rhoc[i])
        for i in IT)

    # Leontief intermediate input demands at base year.
    int0 = Dict{Symbol,Float64}(i => sum(io[i, j]*xd0[j] for j in SEC) for i in SEC)

    # CET share parameter (gamma): from export-supply first-order condition at base.
    gamma = Dict{Symbol,Float64}(
        i => 1/(1 + pd0[i]/pe0[i]*(e0[i]/xxd0[i])^(rhot[i]-1))
        for i in IT)
    for i in SEC  # Sectors with zero imports cannot export
        m0[i] == 0 && (gamma[i] = 0.0)
    end

    # Labour shares in Cobb-Douglas value added: alpha_l = wage bill / value added.
    alphl = Dict{Tuple{Symbol,Symbol},Float64}(
        (l, i) => (wdist[i, l]*wa0[l]*xle[i, l]) / (pva0[i]*xd0[i])
        for l in LC, i in SEC)

    # Production function TFP parameter (ad): calibrated from base output quantities.
    qd = Dict{Symbol,Float64}()
    for i in SEC
        qd[i] = (xllb[i, :rural]^alphl[:rural, i]) *
                (xllb[i, :urbanunsk]^alphl[:urbanunsk, i]) *
                (xllb[i, :urbanskil]^alphl[:urbanskil, i]) *
                (k0[i]^(1 - sum(alphl[l, i] for l in LC)))
    end
    ad = Dict{Symbol,Float64}(i => xd0[i] / qd[i] for i in SEC)

    # Recompute ac after ad calibration (ensures ordering consistency, as in the original).
    for i in IT
        x0[i] = pd0[i]*xxd0[i] + pm0[i]*m0[i]
    end
    for i in IT
        ac[i] = x0[i] / (delta[i]*m0[i]^(-rhoc[i]) + (1 - delta[i])*xxd0[i]^(-rhoc[i]))^(-1/rhoc[i])
    end

    # CET shift parameter (at): calibrated so cet() holds at base output/exports.
    at = Dict{Symbol,Float64}(
        i => xd0[i] / (gamma[i]*e0[i]^rhot[i] + (1 - gamma[i])*xxd0[i]^rhot[i])^(1/rhot[i])
        for i in IT)

    return Params(SEC, IT, ITN, LC,
                  er, gr0, gdtot0, cdtot0, fsav0, mps0, y0,
                  depr, rhoc, rhot, eta, tm0, te, itax, cles, gles, kio, dstr,
                  m0, e0, xd0, k0, pd0, pm0, pe0, pwm0, pwe0, pva0, xxd0, dst0, id0, cd0, int0, x0,
                  delta, ac, gamma, at, ad,
                  wa0, ls0, io, imat, wdist, xle, alphl)
end

# =============================================================================
# RESULTS CONTAINER
# Bundles all 38 endogenous variable arrays for one model run.
# Any is used because JuMP.value.() returns non-standard DenseAxisArray types.
# =============================================================================

struct Simulation
    pd::Any; pm::Any; pe::Any; pk::Any; px::Any; p::Any; pva::Any
    pwm::Any; pwe::Any; tm::Any; x::Any; xd::Any; xxd::Any; e::Any; m::Any
    k::Any; wa::Any; ls::Any; labd::Any; int::Any; cd::Any; gd::Any; id::Any; dst::Any
    y::Any; gr::Any; tariff::Any; indtax::Any; duty::Any; gdtot::Any; mps::Any
    hhsav::Any; govsav::Any; deprecia::Any; savings::Any; fsav::Any; dk::Any; omega::Any
end

# =============================================================================
# SOLVE
# Builds a fresh JuMP model from a Params, solves it, and extracts a Simulation.
# =============================================================================

"""
    solve(p::Params; silent=true, tol=1e-8, max_iter=3000)
        -> (sim::Simulation, status::String, converged::Bool, iterations, elapsed)

Builds a fresh JuMP model reading every parameter from `p` (never from a global),
solves it with Ipopt, and returns the solution plus solver diagnostics. `status` is
the raw JuMP `termination_status` as a string (e.g. `"LOCALLY_SOLVED"`); `converged`
is `true` for `OPTIMAL`, `LOCALLY_SOLVED`, or `ALMOST_LOCALLY_SOLVED` (Ipopt reports
"Solved To Acceptable Level" for both the baseline and `sim1` at default tolerances --
see the exploration report). `p` itself is never mutated.
"""
function solve(P::Params; silent::Bool=true, tol::Float64=1e-8, max_iter::Int=3000)
    SEC, IT, ITN, LC = P.SEC, P.IT, P.ITN, P.LC
    io, imat, wdist, xle, alphl = P.io, P.imat, P.wdist, P.xle, P.alphl
    depr, rhoc, rhot, eta   = P.depr, P.rhoc, P.rhot, P.eta
    tm0, te, itax           = P.tm0, P.te, P.itax
    cles, gles, kio, dstr   = P.cles, P.gles, P.kio, P.dstr
    m0, e0, xd0, k0         = P.m0, P.e0, P.xd0, P.k0
    pd0, pm0, pe0, pva0     = P.pd0, P.pm0, P.pe0, P.pva0
    pwm0, pwe0              = P.pwm0, P.pwe0
    x0, xxd0, dst0, id0, cd0, int0 = P.x0, P.xxd0, P.dst0, P.id0, P.cd0, P.int0
    delta, ac, gamma, at, ad = P.delta, P.ac, P.gamma, P.at, P.ad
    wa0, ls0                = P.wa0, P.ls0
    er, gr0, gdtot0, fsav0, mps0, y0 = P.er, P.gr0, P.gdtot0, P.fsav0, P.mps0, P.y0

    cgecam = Model(Ipopt.Optimizer)
    silent && set_silent(cgecam)
    set_optimizer_attribute(cgecam, "tol", tol)
    set_optimizer_attribute(cgecam, "max_iter", max_iter)

    # -------------------------------------------------------------------------
    # DECISION VARIABLES  (all strictly positive via lower bound 1e-6)
    # -------------------------------------------------------------------------
    @variables cgecam begin
        # Prices
        pd[i in SEC]  >= 1e-6, (start = pd0[i])   # Domestic goods price
        pm[i in SEC]  >= 1e-6, (start = pm0[i])   # Domestic price of imports
        pe[i in SEC]  >= 1e-6, (start = pe0[i])   # Domestic price of exports
        pk[i in SEC]  >= 1e-6, (start = pd0[i])   # Capital rental rate by sector
        px[i in SEC]  >= 1e-6, (start = pd0[i])   # Average output price (net of indirect tax)
        p[i in SEC]   >= 1e-6, (start = pd0[i])   # Composite (Armington) goods price
        pva[i in SEC] >= 1e-6, (start = pva0[i])  # Value-added price
        pwm[i in SEC] >= 1e-6, (start = pwm0[i])  # World import price (foreign currency)
        pwe[i in SEC] >= 1e-6, (start = pwe0[i])  # World export price (foreign currency)
        tm[i in SEC]  >= 1e-6, (start = tm0[i])   # Tariff rates

        # Production quantities
        x[i in SEC]   >= 1e-6, (start = x0[i])    # Composite good supply (Armington)
        xd[i in SEC]  >= 1e-6, (start = xd0[i])   # Domestic output by sector
        xxd[i in SEC] >= 1e-6, (start = xxd0[i])  # Domestic sales (output minus exports)
        e[i in SEC]   >= 1e-6, (start = e0[i])    # Exports by sector
        m[i in SEC]   >= 1e-6, (start = m0[i])    # Imports by sector

        # Factors
        k[i in SEC]            >= 1e-6, (start = k0[i])    # Capital stock by sector
        wa[l in LC]            >= 1e-6, (start = wa0[l])   # Economy-wide wage by labour type
        ls[l in LC]            >= 1e-6, (start = ls0[l])   # Labour supply by category
        labd[i in SEC, l in LC] >= 1e-6, (start = xle[i,l]) # Employment (sector × labour type)

        # Demand aggregates
        int[i in SEC]  >= 1e-6, (start = int0[i])   # Intermediate input demand
        cd[i in SEC]   >= 1e-6, (start = cd0[i])    # Private consumption demand
        gd[i in SEC]   >= 1e-6                       # Government consumption demand
        id[i in SEC]   >= 1e-6, (start = id0[i])    # Investment demand by sector of origin
        dst[i in SEC]  >= 1e-6, (start = dst0[i])   # Inventory investment
        y              >= 1e-6, (start = y0)          # Private GDP (value added minus depreciation)
        gr             >= 1e-6, (start = gr0)         # Government revenue
        tariff         >= 1e-6, (start = 76.548)      # Tariff revenue
        indtax         >= 1e-6                         # Indirect tax revenue
        duty           >= 1e-6                         # Export duty revenue
        gdtot          >= 1e-6                         # Total government consumption volume
        mps            >= 1e-6                         # Marginal propensity to save (households)
        hhsav          >= 1e-6                         # Household savings
        govsav         >= 1e-6                         # Government savings (budget surplus/deficit)
        deprecia       >= 1e-6                         # Economy-wide depreciation expenditure
        savings        >= 1e-6                         # Total savings (= total investment)
        fsav           >= 1e-6, (start = fsav0)       # Foreign savings (current-account deficit)
        dk[i in SEC]   >= 1e-6                         # Investment by sector of destination

        # Welfare
        omega                                           # Cobb-Douglas utility (welfare indicator)
    end

    # -------------------------------------------------------------------------
    # EQUILIBRIUM CONDITIONS
    # -------------------------------------------------------------------------
    @NLconstraints cgecam begin

        # --- Price Block ---

        # Import price: domestic = world price × exchange rate × (1 + tariff)
        pmdef[it in IT],          pm[it] == pwm[it]*er*(1 + tm[it])

        # Export price: world price = domestic export price ÷ (er × (1 + export duty))
        pedef[it in IT],          pe[it]*(1 + te[it]) == pwe[it]*er

        # Armington price identity: value of composite supply = domestic sales + imports
        absorption1[i in IT],     p[i]*x[i] == pd[i]*xxd[i] + pm[i]*m[i]

        # Non-traded: composite supply equals domestic sales only (no imports)
        absorption2[itn in ITN],  p[itn]*x[itn] == pd[itn]*xxd[itn]

        # Revenue identity: output revenue = domestic sales revenue + export revenue
        sales[i in SEC],          px[i]*xd[i] == pd[i]*xxd[i] + pe[i]*e[i]

        # Activity price: net-of-tax output price = value added + intermediate cost
        actp[i in SEC],           px[i]*(1 - itax[i]) == pva[i] + sum(io[j,i]*p[j] for j in SEC)

        # Capital goods price: weighted average of composite goods prices (via imat)
        pkdef[i in SEC],          pk[i] == sum(p[j]*imat[j,i] for j in SEC)

        # --- Production and Factor Block ---

        # Cobb-Douglas production with heterogeneous labour categories and fixed capital:
        #   xd[i] = ad[i] * Π_l( labd[i,l]^alphl[l,i] ) * k[i]^(1 - Σ_l alphl[l,i])
        activity[i in SEC],       xd[i] == ad[i] * prod(labd[i,l]^alphl[l,i] for l in LC) *
                                            k[i]^(1 - sum(alphl[l,i] for l in LC))

        # Profit maximisation (Shephard's lemma / FOC for labour):
        #   wage × wdist × employment = value-added price × output × labour share
        profitmax[i in SEC, l in LC], wa[l]*wdist[i,l]*labd[i,l] == xd[i]*pva[i]*alphl[l,i]

        # Labour market clearing: total sectoral demand equals exogenous supply
        lmequil[l in LC],         sum(labd[i,l] for i in SEC) == ls[l]

        # CET function: allocates output between domestic market and exports.
        # Substitution elasticity = 1/(rhot-1).
        cet[it in IT],            xd[it] == at[it] * (gamma[it]*e[it]^rhot[it] +
                                             (1 - gamma[it])*xxd[it]^rhot[it])^(1/rhot[it])

        # Export demand curve: foreign demand is downward-sloping in world price
        edemand[it in IT],        e[it]/e0[it] == (pwe0[it]/pwe[it])^eta[it]

        # Export supply ratio: export/domestic split from CET profit maximisation
        esupply[it in IT],        e[it]/xxd[it] == (pe[it]/pd[it] * (1-gamma[it])/gamma[it])^(1/(rhot[it]-1))

        # Armington function: CES aggregation of imports and domestically-produced goods
        armington[it in IT],      x[it] == ac[it] * (delta[it]*m[it]^(-rhoc[it]) +
                                            (1-delta[it])*xxd[it]^(-rhoc[it]))^(-1/rhoc[it])

        # Cost minimisation FOC: import/domestic ratio from Armington expenditure minimisation
        costmin[it in IT],        m[it]/xxd[it] == (pd[it]/pm[it] * delta[it]/(1-delta[it]))^(1/(1+rhoc[it]))

        # Non-traded sector: domestic sales equal total output (no trade possible)
        xxdsn[itr in ITN],        xxd[itr] == xd[itr]
        xsn[itr in ITN],          x[itr]   == xxd[itr]

        # --- Demand Block ---

        # Intermediate inputs: Leontief (fixed io coefficients), summed over using sectors
        inteq[j in SEC],          int[j] == sum(io[j,i]*xd[i] for i in SEC)

        # Inventory investment: fixed proportion of sectoral gross output
        dsteq[i in SEC],          dst[i] == dstr[i]*xd[i]

        # Private consumption: linear Engel curves, share cles[i] of disposable income
        cdeq[i in SEC],           p[i]*cd[i] == cles[i]*(1 - mps)*y

        # Private GDP: aggregate value added minus economy-wide depreciation
        gdp,                      y == sum(pva[i]*xd[i] for i in SEC) - deprecia

        # Household savings: fixed marginal propensity to save applied to income
        hhsaveq,                  hhsav == mps*y

        # Government revenue: sum of all tax instruments
        greq,                     gr == tariff + duty + indtax

        # Government budget: revenue finances consumption plus budget surplus/deficit
        gruse,                    gr == sum(p[i]*gd[i] for i in SEC) + govsav

        # Government consumption: each commodity gets a fixed share of total volume
        gdeq[i in SEC],           gd[i] == gles[i]*gdtot

        # Tariff revenue: ad-valorem tariffs applied to import values at world prices
        tariffdef,                tariff == sum(tm[it]*m[it]*pwm[it] for it in IT)*er

        # Indirect tax revenue: ad-valorem production taxes on gross output values
        indtaxdef,                indtax == sum(itax[i]*px[i]*xd[i] for i in SEC)

        # Export duty revenue (zero at base since te = 0 everywhere)
        dutydef,                  duty == sum(te[it]*e[it]*pe[it] for it in IT)

        # Aggregate depreciation expenditure (capital consumption)
        depreq,                   deprecia == sum(depr[i]*pk[i]*k[i] for i in SEC)

        # Savings-investment balance: total savings = total investment
        totsav,                   savings == hhsav + govsav + deprecia + fsav*er

        # Investment by sector of destination: savings allocated by fixed shares kio,
        # net of inventory investment (inventories are not part of capital formation)
        prodinv[i in SEC],        pk[i]*dk[i] == kio[i]*savings - kio[i]*sum(dst[j]*p[j] for j in SEC)

        # Investment demand by sector of origin: derived from capital composition matrix
        ieq[i in SEC],            id[i] == sum(imat[i,j]*dk[j] for j in SEC)

        # Current account balance: import bill = export revenues + net foreign savings
        caeq,                     sum(pwm[it]*m[it] for it in IT) == sum(pwe[it]*e[it] for it in IT) + fsav

        # Market clearing: composite supply = all categories of final and intermediate demand
        equil[i in SEC],          x[i] == int[i] + cd[i] + gd[i] + id[i] + dst[i]

        # Welfare indicator: Cobb-Douglas utility over private consumption quantities
        obj,                      omega == prod(cd[i]^cles[i] for i in SEC if cles[i] > 0.0)

        # --- Closure Rules ---
        # Fix exogenous variables; determines which variables adjust to clear markets.
        closurek[i in SEC],   k[i]    == k0[i]      # Capital stocks fixed (static model)
        closurep[i in SEC],   pwm[i]  == pwm0[i]    # World prices fixed (small open economy)
        closurel[l in LC],    ls[l]   == ls0[l]     # Labour supply fixed (no labour-leisure)
        closuret[it in IT],   tm[it]  == tm0[it]    # Tariffs fixed (policy instrument)
        closuref,             fsav    == fsav0       # Foreign savings exogenous
        closuremp,             mps     == mps0        # Household saving rate fixed
        closureg,             gdtot   == gdtot0      # Government spending volume fixed
        closurem[itn in ITN], m[itn]  == 0           # Non-traded: no imports
        closurlp,             labd[:publiques, :rural]    == 0  # No rural labour in public sector
        closurla,             labd[:agsubsist, :urbanskil] == 0  # No skilled labour in subsistence farming
        closuree[itn in ITN], e[itn]  == 0           # Non-traded: no exports
    end

    # Feasibility problem: solve for any point satisfying all equilibrium conditions
    @NLobjective(cgecam, Min, 1)

    elapsed = @elapsed JuMP.optimize!(cgecam)

    raw_ts = termination_status(cgecam)
    status = string(raw_ts)
    converged = raw_ts in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED)
    iterations = try
        MOI.get(cgecam, MOI.BarrierIterations())
    catch
        -1
    end

    # Extract solution values from JuMP variable containers
    pd  = JuMP.value.(pd);   pm  = JuMP.value.(pm);   pe  = JuMP.value.(pe)
    pk  = JuMP.value.(pk);   px  = JuMP.value.(px);   p   = JuMP.value.(p)
    pva = JuMP.value.(pva);  pwm = JuMP.value.(pwm);  pwe = JuMP.value.(pwe)
    tm  = JuMP.value.(tm)

    x   = JuMP.value.(x);   xd  = JuMP.value.(xd);   xxd = JuMP.value.(xxd)
    e   = JuMP.value.(e);   m   = JuMP.value.(m)

    k    = JuMP.value.(k);   wa   = JuMP.value.(wa);   ls   = JuMP.value.(ls)
    labd = JuMP.value.(labd)

    int  = JuMP.value.(int);  cd   = JuMP.value.(cd);   gd   = JuMP.value.(gd)
    id   = JuMP.value.(id);   dst  = JuMP.value.(dst);  y    = JuMP.value.(y)
    gr   = JuMP.value.(gr)

    tariff   = JuMP.value.(tariff);   indtax   = JuMP.value.(indtax)
    duty     = JuMP.value.(duty);     gdtot    = JuMP.value.(gdtot)
    mps      = JuMP.value.(mps);      hhsav    = JuMP.value.(hhsav)
    govsav   = JuMP.value.(govsav);   deprecia = JuMP.value.(deprecia)
    savings  = JuMP.value.(savings);  fsav     = JuMP.value.(fsav)
    dk       = JuMP.value.(dk);       omega    = JuMP.value.(omega)

    sim = Simulation(pd, pm, pe, pk, px, p, pva, pwm, pwe, tm,
                     x, xd, xxd, e, m, k, wa, ls, labd,
                     int, cd, gd, id, dst, y, gr, tariff, indtax, duty, gdtot,
                     mps, hhsav, govsav, deprecia, savings, fsav, dk, omega)

    return sim, status, converged, iterations, elapsed
end

# =============================================================================
# SHOCKS
# =============================================================================

"""
    with_shocks(p::Params, overrides::Dict{Symbol,Any}) -> Params

Returns a `deepcopy` of `p` with `overrides` applied; `p` itself is never mutated, so
running scenarios back to back (even concurrently) cannot leak state between them.
Each key of `overrides` must name a field of `Params`:

- a scalar value (e.g. `:gdtot0 => 145.0`, `:er => 0.25`, `:fsav0 => 55.0`) replaces a
  scalar field outright;
- a `Dict{Symbol,<:Real}` (e.g. `:k0 => Dict(:agsubsist => 363.5)`) merges new LEVELS
  into the given sector/labour keys of that field, leaving every other key untouched.

Values are always the new LEVEL of the parameter, never a percent change -- callers
compute the level from whatever shock convention they use (see `apply_mode` in the
app's `Common.jl`).
"""
function with_shocks(p::Params, overrides::Dict{Symbol,Any})
    p2 = deepcopy(p)
    for (field, val) in overrides
        hasfield(Params, field) || error("with_shocks: Params has no field :$field")
        if val isa AbstractDict
            cur = getfield(p2, field)
            for (k, v) in val
                cur[Symbol(k)] = Float64(v)
            end
        else
            setfield!(p2, field, Float64(val))
        end
    end
    return p2
end

# =============================================================================
# EXAMPLE
# Reproduces the baseline + sim1 solves the original script ran at include-time.
# =============================================================================

"""
    example()

Loads `data/camdata.xlsx`, calibrates it, solves the baseline, then solves `sim1`
("public investment in subsistence agriculture": a +10% capital-stock shock to
`agsubsist`), prints a short summary of each, and returns
`(baseline = ..., sim1 = ...)` as `Simulation`s. This is the moral equivalent of the
original script's include-time baseline/sim1 solves, just callable on demand instead
of running automatically.
"""
function example()
    raw = load_data(joinpath(@__DIR__, "data", "camdata.xlsx"))
    params = calibrate(raw)

    baseline, status0, converged0, iters0, t0 = solve(params)
    println("baseline: status=$status0 converged=$converged0 iterations=$iters0 ($(round(t0, digits=3))s)")
    println("  y=$(baseline.y)  omega=$(baseline.omega)  gr=$(baseline.gr)")

    shocked = with_shocks(params, Dict{Symbol,Any}(:k0 => Dict(:agsubsist => params.k0[:agsubsist] * 1.10)))
    sim1, status1, converged1, iters1, t1 = solve(shocked)
    println("sim1 (+10% capital stock, agsubsist): status=$status1 converged=$converged1 iterations=$iters1 ($(round(t1, digits=3))s)")
    println("  y=$(sim1.y)  omega=$(sim1.omega)  gr=$(sim1.gr)")

    return (baseline = baseline, sim1 = sim1)
end

end # module CGECameroon
