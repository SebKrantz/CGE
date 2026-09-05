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
# CAMEROON DEFAULTS (legacy loader only)
# =============================================================================
# The model itself is N-sector / L-labour-category: `calibrate` and `_solve_once` read
# every set off the `RawData`/`Params` they are handed, never off these constants. The
# constants survive only as the defaults `_load_legacy` stamps onto the legacy Cameroon
# workbook, which carries no `sectors`/`labour`/`scalars` sheets of its own. A generic
# workbook (see `load_data`) supplies all of them as data, for any sector count, any
# traded/non-traded partition and any number of labour categories.
#
# SEC  — the 11 Cameroon production sectors
# IT   — its 9 traded sectors (subject to Armington imports and CET exports)
# ITN  — its 2 non-traded sectors (domestic supply = domestic demand, no trade)
# LC   — its 3 labour categories
# WA0  — its base wage per labour category (million CFAF per 1000 workers)
# SCALARS — its economy-wide scalars, also the fallback for any name a generic
#           workbook's `scalars` sheet omits
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

const WA0 = Dict{Symbol,Float64}(:rural => 0.11, :urbanunsk => 0.15678, :urbanskil => 1.8657)

const SCALARS = Dict{Symbol,Float64}(
    :er     => 0.21,     # Real exchange rate: CFAF per US dollar
    :gr0    => 179.0,    # Government revenue (billion CFAF) -- start value only
    :gdtot0 => 135.03,   # Total government consumption (billion CFAF)
    :cdtot0 => 947.98,   # Total private consumption (billion CFAF)
    :fsav0  => 36.841,   # Foreign savings / current-account deficit (billion USD)
    :mps0   => 0.09305,  # Household saving rate (was a bare literal in `closuremp`)
    :td0    => 0.0,      # Household direct-tax rate (Cameroon 1979-80 has no direct tax)
)

"Row names of the `miscellaneous` sheet, in the legacy workbook's fixed row order."
const MISC_ROWS = [:m0, :e0, :xd0, :k, :depr, :rhoc, :rhot, :eta, :pd0, :tm0,
                   :itax, :cles, :gles, :kio, :dstr, :dst, :id]

"Names accepted in the `value` column of a generic workbook's `scalars` sheet, on top of
one `wa0_<labour code>` per labour category. `tariff0` is the only one with no entry in
`SCALARS`: omitted, `calibrate` derives it from the data (`Σ_IT m0·tm0/(1+tm0)`)."
const SCALAR_NAMES = (:er, :gr0, :gdtot0, :cdtot0, :fsav0, :mps0, :td0, :tariff0)

# =============================================================================
# RAW DATA (returned by `load_data`)
# =============================================================================

"""
    RawData

One base-year dataset, straight off the workbook: the sector set, the traded/non-traded
partition, the labour categories, the five data matrices, the base wage per labour
category and the economy-wide scalars. No calibration algebra has run yet; this is
pass-through data, and it -- not any module constant -- is the source of truth for every
set the rest of the model indexes over.
"""
struct RawData
    SEC::Vector{Symbol}
    IT::Vector{Symbol}
    ITN::Vector{Symbol}
    LC::Vector{Symbol}
    io::Dict{Tuple{Symbol,Symbol},Float64}     # iotable sheet, keyed (supplier, user)
    imat::Dict{Tuple{Symbol,Symbol},Float64}   # imat sheet, keyed (origin, destination)
    wdist::Dict{Tuple{Symbol,Symbol},Float64}  # wagedist sheet, keyed (sector, labour)
    xle::Dict{Tuple{Symbol,Symbol},Float64}    # employment sheet, keyed (sector, labour)
    zz::Dict{Tuple{Symbol,Symbol},Float64}     # miscellaneous sheet, keyed (MISC_ROWS name, sector)
    wa0::Dict{Symbol,Float64}                  # base wage by labour category
    scalars::Dict{Symbol,Float64}              # economy-wide scalars (SCALAR_NAMES)
end

"""
    load_data(path::AbstractString) -> RawData

Reads a base-year workbook into a `RawData`, in either of two layouts, chosen by whether
the workbook has a `sectors` sheet:

**Generic (has a `sectors` sheet)** — any number of sectors and labour categories, in any
order, all read by *header lookup* rather than by cell position:

| Sheet | Shape | Contents |
|---|---|---|
| `sectors` | header row + one row per sector | `code`, `label`, `traded` (TRUE/FALSE) — defines `SEC`, and `IT`/`ITN` from `traded` |
| `labour` | header row + one row per labour category | `code`, `label` — defines `LC` |
| `iotable` | labelled N×N | `io[i,j]`, rows = supplying sector, columns = using sector |
| `imat` | labelled N×N | `imat[i,j]`, rows = origin sector, columns = destination sector |
| `employment` | labelled N×L | `xle[i,l]` |
| `wagedist` | labelled N×L | `wdist[i,l]` |
| `miscellaneous` | labelled 17×N | rows named `m0 e0 xd0 k depr rhoc rhot eta pd0 tm0 itax cles gles kio dstr dst id`, columns = sector codes |
| `scalars` | header row + `name`,`value` rows | `er gr0 gdtot0 cdtot0 fsav0 mps0 td0 tariff0` and `wa0_<labour code>`; any omitted name falls back to the Cameroon literal in `SCALARS`/`WA0` (`wa0_*` to 1.0 for a labour code Cameroon does not have, `tariff0` to the value `calibrate` derives from the data) |

On a labelled sheet, cell A1 is a free corner, row 1 holds the column headers and column A
the row labels; order is irrelevant and a missing header, duplicate header or blank /
non-numeric data cell is an error naming the sheet and the offending key.

**Legacy (no `sectors` sheet)** — the original Cameroon `camdata.xlsx`: sheets `iotable`,
`imat`, `wagedist`, `employment`, `miscellaneous` read at the original hardcoded cell
ranges and row order, with `SEC`/`IT`/`ITN`/`LC`/`WA0`/`SCALARS` supplied as the module's
Cameroon defaults. Unchanged, byte for byte, from the pre-`n-sector` loader.
"""
function load_data(path::AbstractString)
    xf = XLSX.readxlsx(path)
    return "sectors" in XLSX.sheetnames(xf) ? _load_generic(xf) : _load_legacy(path)
end

function _load_legacy(path::AbstractString)
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
    zz = Dict{Tuple{Symbol,Symbol},Float64}()
    for i in 1:length(MISC_ROWS), j in 1:length(SEC)
        zz[MISC_ROWS[i], SEC[j]] = datazz[i, j]
    end

    return RawData(SEC, IT, ITN, LC, io, imat, wdist, xle, zz, copy(WA0), copy(SCALARS))
end

function _load_generic(xf)
    SECg, ITg, ITNg = _read_sectors(xf)
    LCg = _read_labour(xf)
    io    = _labelled(xf, "iotable",       SECg,      SECg)
    imat  = _labelled(xf, "imat",          SECg,      SECg)
    wdist = _labelled(xf, "wagedist",      SECg,      LCg)
    xle   = _labelled(xf, "employment",    SECg,      LCg)
    zz    = _labelled(xf, "miscellaneous", MISC_ROWS, SECg)
    wa0, scalars = _read_scalars(xf, LCg)
    return RawData(SECg, ITg, ITNg, LCg, io, imat, wdist, xle, zz, wa0, scalars)
end

# --- generic-workbook cell helpers -------------------------------------------------

"A header/label cell as a `Symbol` (whitespace trimmed), or `nothing` if it is blank."
function _key(v)
    v === missing && return nothing
    s = strip(string(v))
    return isempty(s) ? nothing : Symbol(s)
end

"Row-1 headers of `t` as lowercase `Symbol => column index`."
function _headers(t::AbstractMatrix)
    d = Dict{Symbol,Int}()
    for j in axes(t, 2)
        k = _key(t[1, j])
        k === nothing && continue
        d[Symbol(lowercase(String(k)))] = j
    end
    return d
end

"A data cell as a `Float64`, erroring with the sheet and cell named if it is not one."
function _num(v, sheet::AbstractString, what::AbstractString)
    v isa Real && return Float64(v)
    if v isa AbstractString
        x = tryparse(Float64, strip(v))
        x === nothing || return x
    end
    error("load_data: sheet `$sheet`, $what: expected a number, got $(repr(v))")
end

"A `traded` cell as a `Bool`, accepting Excel booleans, 0/1 and TRUE/FALSE-style strings."
function _bool(v, sheet::AbstractString, what::AbstractString)
    v isa Bool && return v
    v isa Real && return v != 0
    if v isa AbstractString
        s = lowercase(strip(v))
        s in ("true", "t", "yes", "y", "1")  && return true
        s in ("false", "f", "no", "n", "0")  && return false
    end
    error("load_data: sheet `$sheet`, $what: expected TRUE or FALSE, got $(repr(v))")
end

function _sheet(xf, sheet::AbstractString)
    sheet in XLSX.sheetnames(xf) ||
        error("load_data: workbook has no `$sheet` sheet (its sheets are " *
              join(XLSX.sheetnames(xf), ", ") * ")")
    return xf[sheet][:]
end

"""
    _labelled(xf, sheet, rowkeys, colkeys) -> Dict{Tuple{Symbol,Symbol},Float64}

Reads a labelled matrix sheet -- cell A1 a free corner, row 1 the column headers, column A
the row labels -- BY HEADER LOOKUP, so the order of rows and columns in the file is
irrelevant, and returns every `(rowkey, colkey)` cell of `rowkeys × colkeys`. Rows and
columns beyond those keys are ignored; a missing or duplicated header, or a blank /
non-numeric data cell, is an error naming the sheet and the offending key.
"""
function _labelled(xf, sheet::AbstractString, rowkeys::Vector{Symbol}, colkeys::Vector{Symbol})
    t = _sheet(xf, sheet)
    index(keys_, cells, what) = begin
        at = Dict{Symbol,Int}()
        for (n, v) in cells
            k = _key(v)
            k === nothing && continue
            haskey(at, k) && error("load_data: sheet `$sheet` has two $what labelled `$k`")
            at[k] = n
        end
        for k in keys_
            haskey(at, k) ||
                error("load_data: sheet `$sheet` has no $(what[1:end-1]) labelled `$k`" *
                      " (it has " * join(sort!(string.(collect(keys(at)))), ", ") * ")")
        end
        at
    end
    colat = index(colkeys, ((j, t[1, j]) for j in 2:size(t, 2)), "columns")
    rowat = index(rowkeys, ((i, t[i, 1]) for i in 2:size(t, 1)), "rows")

    d = Dict{Tuple{Symbol,Symbol},Float64}()
    for r in rowkeys, c in colkeys
        d[r, c] = _num(t[rowat[r], colat[c]], sheet, "row `$r` × column `$c`")
    end
    return d
end

"Reads the `sectors` sheet into (SEC, IT, ITN); `IT`/`ITN` come from the `traded` column."
function _read_sectors(xf)
    t = _sheet(xf, "sectors")
    h = _headers(t)
    for k in (:code, :traded)
        haskey(h, k) || error("load_data: sheet `sectors` needs a `$k` column")
    end
    SECg, ITg, ITNg = Symbol[], Symbol[], Symbol[]
    for i in 2:size(t, 1)
        code = _key(t[i, h[:code]])
        code === nothing && continue
        code in SECg && error("load_data: sheet `sectors` lists sector `$code` twice")
        push!(SECg, code)
        traded = _bool(t[i, h[:traded]], "sectors", "row $i (`$code`), column `traded`")
        push!(traded ? ITg : ITNg, code)
    end
    isempty(SECg) && error("load_data: sheet `sectors` lists no sectors")
    isempty(ITg)  && error("load_data: sheet `sectors` marks no sector as traded")
    return SECg, ITg, ITNg
end

"Reads the `labour` sheet into LC."
function _read_labour(xf)
    t = _sheet(xf, "labour")
    h = _headers(t)
    haskey(h, :code) || error("load_data: sheet `labour` needs a `code` column")
    LCg = Symbol[]
    for i in 2:size(t, 1)
        code = _key(t[i, h[:code]])
        code === nothing && continue
        code in LCg && error("load_data: sheet `labour` lists labour category `$code` twice")
        push!(LCg, code)
    end
    isempty(LCg) && error("load_data: sheet `labour` lists no labour categories")
    return LCg
end

"""
Reads the optional `scalars` sheet (`name`/`value` columns) into `(wa0, scalars)`. Every
name it omits falls back to the Cameroon literal: `SCALARS[name]` for an economy-wide
scalar, `WA0[l]` (or 1.0, for a labour code Cameroon does not have) for `wa0_<l>`.
`tariff0` has no fallback here -- left out of `scalars`, `calibrate` derives it.
"""
function _read_scalars(xf, LCg::Vector{Symbol})
    scalars = copy(SCALARS)
    wa0 = Dict{Symbol,Float64}(l => get(WA0, l, 1.0) for l in LCg)
    "scalars" in XLSX.sheetnames(xf) || return wa0, scalars

    t = xf["scalars"][:]
    h = _headers(t)
    for k in (:name, :value)
        haskey(h, k) || error("load_data: sheet `scalars` needs a `$k` column")
    end
    for i in 2:size(t, 1)
        name = _key(t[i, h[:name]])
        name === nothing && continue
        v = _num(t[i, h[:value]], "scalars", "row $i (`$name`), column `value`")
        s = String(name)
        if startswith(s, "wa0_")
            l = Symbol(s[5:end])
            l in LCg || error("load_data: sheet `scalars` sets `$name`, but `$l` is not a " *
                              "labour category in the `labour` sheet")
            wa0[l] = v
        else
            name in SCALAR_NAMES ||
                error("load_data: sheet `scalars` has unknown name `$name` (expected one of " *
                      join(SCALAR_NAMES, ", ") * ", or wa0_<labour code>)")
            scalars[name] = v
        end
    end
    return wa0, scalars
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
    fsav0::Float64       # foreign savings, closure-fixed (may be negative: trade surplus)
    mps0::Float64        # household saving rate, closure-fixed (was the literal 0.09305)
    td0::Float64         # household direct-tax rate, closure-fixed (0 => no direct tax)
    tariff0::Float64     # base tariff revenue (start value only), Σ_IT m0·tm0/(1+tm0)
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

    # Base-year economy-wide scalars: from the workbook (a generic workbook's `scalars`
    # sheet, or the Cameroon defaults the legacy loader stamps on), never from a literal
    # here. `SCALARS[name]` is the fallback for anything a hand-built `RawData` omits.
    sc = raw.scalars
    scalar(name::Symbol) = get(sc, name, SCALARS[name])
    wa0    = Dict{Symbol,Float64}(l => raw.wa0[l] for l in LC)
    er     = scalar(:er)      # Real exchange rate (currency per unit foreign currency)
    gr0    = scalar(:gr0)     # Government revenue (start value only)
    gdtot0 = scalar(:gdtot0)  # Total government consumption
    cdtot0 = scalar(:cdtot0)  # Total private consumption
    fsav0  = scalar(:fsav0)   # Foreign savings (negative = base-year trade surplus)
    mps0   = scalar(:mps0)    # Household saving rate
    td0    = scalar(:td0)     # Household direct-tax rate

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

    # Every traded sector needs strictly positive base-year imports AND exports, or the
    # calibration algebra below silently produces NaN/Inf and the solve returns
    # INVALID_MODEL: m0[i] == 0 forces gamma[i] = 0 (line ~"Sectors with zero imports"),
    # which makes `ac[i]` NaN and `esupply`'s (1-gamma)/gamma divide by zero, while
    # e0[i] == 0 gives gamma[i] = 1, at[i] = Inf and an `edemand` that divides by e0.
    # Fail here instead, naming the sector, so the data builder can fix the file.
    for i in IT
        m0[i] > 0 && e0[i] > 0 && continue
        error("calibrate: sector `$i` is marked traded but has m0 = $(m0[i]), e0 = $(e0[i]); " *
              "every traded sector needs m0 > 0 AND e0 > 0. Fix the data: aggregate the " *
              "sector into a neighbour, mark it non-traded (`traded = FALSE` in the " *
              "`sectors` sheet, with m0 = e0 = 0), or give it a token trade flow and " *
              "rebalance its commodity row.")
    end

    # Base tariff revenue -- a JuMP start value only (`tariff` is determined by
    # `tariffdef`), so a workbook may just state it; otherwise derive it from the data.
    # At base prices this is Σ_IT tm0·m0·pwm0·er with pwm0·er = pd0/(1+tm0).
    tariff0 = haskey(sc, :tariff0) ? sc[:tariff0] :
              sum(m0[i]*tm0[i]/(1 + tm0[i]) for i in IT)

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
        q = 1.0
        for l in LC          # was three hardcoded Cameroon labour categories, same product
            q *= xllb[i, l]^alphl[l, i]
        end
        qd[i] = q * (k0[i]^(1 - sum(alphl[l, i] for l in LC)))
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
                  er, gr0, gdtot0, cdtot0, fsav0, mps0, td0, tariff0, y0,
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
    solve(p::Params; silent=true, tol=1e-8, max_iter=3000,
          start::Union{Nothing,Simulation}=nothing, steps::Int=1,
          base::Union{Nothing,Params}=nothing)
        -> (sim::Simulation, status::String, converged::Bool, iterations, elapsed)

Builds a fresh JuMP model reading every parameter from `p` (never from a global),
solves it with Ipopt, and returns the solution plus solver diagnostics. `status` is
the raw JuMP `termination_status` as a string (e.g. `"LOCALLY_SOLVED"`); `converged`
is `true` for `OPTIMAL`, `LOCALLY_SOLVED`, or `ALMOST_LOCALLY_SOLVED` (Ipopt reports
"Solved To Acceptable Level" for both the baseline and `sim1` at default tolerances --
see the exploration report). `p` itself is never mutated.

Every solve of a *shocked* `Params` still starts every variable from the raw base-year
`*0` data by default (`start = nothing`), exactly like the original script. Two knobs
help a modest policy shock converge reliably instead of landing on a spurious
`LOCALLY_INFEASIBLE` (see `test/robustness_grid.jl` and the exploration report for why
that status shows up on a cold start even though the shocked equations are, in fact,
satisfiable a hair away from where Ipopt gives up):

- `start`: a previously-solved `Simulation` (typically the cached baseline) whose values
  seed every variable's start value via `set_start_value`, instead of the base-year `*0`
  data. This is the preferred fix -- it costs nothing extra and does not change what
  equations are solved.
- `steps`: solve a homotopy path in `steps` equal linear increments from `base` (the
  *pre-shock* `Params`) to `p` (the shocked `Params`), warm-starting each increment from
  the previous one's solution. Only use this if `start` alone still fails -- it requires
  `base` and costs `steps` solves instead of one.
"""
function solve(P::Params; silent::Bool=true, tol::Float64=1e-8, max_iter::Int=3000,
               start::Union{Nothing,Simulation}=nothing, steps::Int=1,
               base::Union{Nothing,Params}=nothing)
    steps >= 1 || error("solve: steps must be >= 1")
    if steps == 1
        return _solve_once(P; silent, tol, max_iter, start)
    end
    base === nothing && error("solve: steps > 1 requires `base` (the pre-shock Params) so " *
                              "the shock can be applied gradually via linear interpolation")
    cur_start = start
    local sim, status, converged, iterations, elapsed
    total_elapsed = 0.0
    for k in 1:steps
        frac = k / steps
        Pk = _interp_params(base, P, frac)
        sim, status, converged, iterations, elapsed =
            _solve_once(Pk; silent, tol, max_iter, start = cur_start)
        total_elapsed += elapsed
        converged || break
        cur_start = sim
    end
    return sim, status, converged, iterations, total_elapsed
end

"""
    _interp_params(base::Params, target::Params, frac::Float64) -> Params

Linear homotopy step: a deep copy of `base` with every scalar `Float64` field and every
`Dict` field's values moved a `frac` fraction of the way from `base`'s value to
`target`'s value (fields that don't differ between `base` and `target` -- i.e. everything
`with_shocks` didn't touch -- are unaffected since the interpolation is a no-op). Used
only by `solve`'s `steps` homotopy.
"""
function _interp_params(base::Params, target::Params, frac::Float64)
    p = deepcopy(base)
    for f in fieldnames(Params)
        bv = getfield(base, f)
        tv = getfield(target, f)
        if bv isa AbstractDict
            d = getfield(p, f)
            for k in keys(tv)
                b = get(bv, k, tv[k])
                d[k] = b + frac * (tv[k] - b)
            end
        elseif bv isa Float64
            setfield!(p, f, bv + frac * (tv - bv))
        end
    end
    return p
end

"""
    _solve_once(P::Params; silent, tol, max_iter, start) -> (sim, status, converged, iterations, elapsed)

Builds and solves one JuMP model from `P` -- the actual model definition (was the whole
body of `solve` before warm-starting/homotopy were added). `start`, if given, seeds every
variable's start value from that `Simulation` instead of `P`'s base-year `*0` data.
"""
function _solve_once(P::Params; silent::Bool=true, tol::Float64=1e-8, max_iter::Int=3000,
                     start::Union{Nothing,Simulation}=nothing)
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
    td0, tariff0            = P.td0, P.tariff0

    cgecam = Model(Ipopt.Optimizer)
    silent && set_silent(cgecam)
    set_optimizer_attribute(cgecam, "tol", tol)
    set_optimizer_attribute(cgecam, "max_iter", max_iter)
    # `bound_relax_factor` (default 1e-8): raised to 1e-5 because this model now has
    # ~20 variables genuinely FIXED at exactly 0 (see the "STRUCTURALLY ZERO" block
    # below) -- Ipopt's own doc for this option describes exactly that situation ("if
    # too many bounds are exactly satisfied, this can cause problems for the
    # algorithm"). Verified necessary: with the structurally-zero cells fixed but this
    # left at Ipopt's default (or even at 1e-6), the baseline solve itself
    # intermittently reports a spurious LOCALLY_INFEASIBLE even though the returned
    # point is a perfectly good solution (every real equation satisfied to high
    # precision) -- 1e-5 was the smallest value that converged reliably across repeated
    # checks. Still tiny: every fixed cell lands within ~1e-8 of its true value (see
    # test/runtests.jl's `skip_cell`, which excludes exactly those cells from the
    # reference.json check), and every other variable sits comfortably off its bound in
    # equilibrium, so this is a no-op for them -- confirmed 0 relative difference on
    # every non-structurally-zero cell.
    set_optimizer_attribute(cgecam, "bound_relax_factor", 1e-5)
    # NOTE on Ipopt tuning (see the exploration report and test/robustness_grid.jl): the
    # obvious feasibility-problem knob, `mu_strategy = "adaptive"`, was tried and rejected --
    # it actually breaks the *baseline* solve (LOCALLY_INFEASIBLE instead of
    # ALMOST_LOCALLY_SOLVED) with Ipopt's default monotone barrier schedule working fine, so
    # this model's ill-conditioning is specific enough that a "generally more robust" solver
    # setting is not a safe global default. Warm-starting from a previously-solved
    # `Simulation` (the `start` keyword above) is what actually fixes the grid's spurious
    # LOCALLY_INFEASIBLE cases -- see `solve`'s docstring.

    # -------------------------------------------------------------------------
    # Start-value helpers: read from `start` (a previously-solved Simulation, typically the
    # cached baseline) when given, otherwise fall back to `default` (the base-year `*0` data,
    # or `nothing` for the handful of variables the original script never seeded).
    # -------------------------------------------------------------------------
    sv0(field::Symbol, default) = start === nothing ? default : getfield(start, field)
    sv1(field::Symbol, i, default) = start === nothing ? default : getfield(start, field)[i]
    sv2(field::Symbol, i, l, default) = start === nothing ? default : getfield(start, field)[i, l]

    # -------------------------------------------------------------------------
    # DECISION VARIABLES  (strictly positive via lower bound 1e-6, except where noted)
    #
    # `hhsav`, `govsav` and `fsav` are declared FREE. All three are accounting residuals
    # that a perfectly ordinary base year can show as negative -- a government running a
    # deficit, an economy running a trade surplus, a household dissaving -- and a `>= 1e-6`
    # bound on them made every such dataset infeasible rather than merely unusual (the
    # identities analysis' §9 pitfalls 1a-1c). `td` is free for the same reason plus one
    # more: `closuretd` pins it to `td0`, which is 0 whenever there is no direct tax, and
    # a bound of 1e-6 against an equation demanding exactly 0 is the bound-vs-equation
    # conflict the "STRUCTURALLY ZERO" block below exists to avoid. Ipopt's
    # `bound_relax_factor` (set above) is unaffected: it relaxes bounds that exist, and
    # these four now have none.
    # -------------------------------------------------------------------------
    @variables cgecam begin
        # Prices
        pd[i in SEC]  >= 1e-6, (start = sv1(:pd, i, pd0[i]))   # Domestic goods price
        pm[i in SEC]  >= 1e-6, (start = sv1(:pm, i, pm0[i]))   # Domestic price of imports
        pe[i in SEC]  >= 1e-6, (start = sv1(:pe, i, pe0[i]))   # Domestic price of exports
        pk[i in SEC]  >= 1e-6, (start = sv1(:pk, i, pd0[i]))   # Capital rental rate by sector
        px[i in SEC]  >= 1e-6, (start = sv1(:px, i, pd0[i]))   # Average output price (net of indirect tax)
        p[i in SEC]   >= 1e-6, (start = sv1(:p, i, pd0[i]))    # Composite (Armington) goods price
        pva[i in SEC] >= 1e-6, (start = sv1(:pva, i, pva0[i])) # Value-added price
        pwm[i in SEC] >= 1e-6, (start = sv1(:pwm, i, pwm0[i])) # World import price (foreign currency)
        pwe[i in SEC] >= 1e-6, (start = sv1(:pwe, i, pwe0[i])) # World export price (foreign currency)
        tm[i in SEC]  >= 1e-6, (start = sv1(:tm, i, tm0[i]))   # Tariff rates

        # Production quantities
        x[i in SEC]   >= 1e-6, (start = sv1(:x, i, x0[i]))     # Composite good supply (Armington)
        xd[i in SEC]  >= 1e-6, (start = sv1(:xd, i, xd0[i]))   # Domestic output by sector
        xxd[i in SEC] >= 1e-6, (start = sv1(:xxd, i, xxd0[i])) # Domestic sales (output minus exports)
        e[i in SEC]   >= 1e-6, (start = sv1(:e, i, e0[i]))     # Exports by sector
        m[i in SEC]   >= 1e-6, (start = sv1(:m, i, m0[i]))     # Imports by sector

        # Factors
        k[i in SEC]            >= 1e-6, (start = sv1(:k, i, k0[i]))     # Capital stock by sector
        wa[l in LC]            >= 1e-6, (start = sv1(:wa, l, wa0[l]))   # Economy-wide wage by labour type
        ls[l in LC]            >= 1e-6, (start = sv1(:ls, l, ls0[l]))   # Labour supply by category
        labd[i in SEC, l in LC] >= 1e-6, (start = sv2(:labd, i, l, xle[i,l])) # Employment (sector × labour type)

        # Demand aggregates
        int[i in SEC]  >= 1e-6, (start = sv1(:int, i, int0[i]))  # Intermediate input demand
        cd[i in SEC]   >= 1e-6, (start = sv1(:cd, i, cd0[i]))    # Private consumption demand
        gd[i in SEC]   >= 1e-6, (start = sv1(:gd, i, nothing))    # Government consumption demand
        id[i in SEC]   >= 1e-6, (start = sv1(:id, i, id0[i]))    # Investment demand by sector of origin
        dst[i in SEC]  >= 1e-6, (start = sv1(:dst, i, dst0[i]))  # Inventory investment
        y              >= 1e-6, (start = sv0(:y, y0))             # Private GDP (value added minus depreciation)
        gr             >= 1e-6, (start = sv0(:gr, gr0))           # Government revenue
        tariff         >= 1e-6, (start = sv0(:tariff, tariff0))    # Tariff revenue
        indtax         >= 1e-6, (start = sv0(:indtax, nothing))    # Indirect tax revenue
        duty           >= 1e-6, (start = sv0(:duty, nothing))      # Export duty revenue
        gdtot          >= 1e-6, (start = sv0(:gdtot, nothing))     # Total government consumption volume
        mps            >= 1e-6, (start = sv0(:mps, nothing))       # Marginal propensity to save (households)
        td,             (start = td0)                              # Household direct-tax rate (free; see above)
        hhsav,          (start = sv0(:hhsav, nothing))             # Household savings (free)
        govsav,         (start = sv0(:govsav, nothing))            # Government savings (budget surplus/deficit; free)
        deprecia       >= 1e-6, (start = sv0(:deprecia, nothing))  # Economy-wide depreciation expenditure
        savings        >= 1e-6, (start = sv0(:savings, nothing))   # Total savings (= total investment)
        fsav,           (start = sv0(:fsav, fsav0))                # Foreign savings (current-account deficit; free)
        dk[i in SEC]   >= 1e-6, (start = sv1(:dk, i, nothing))     # Investment by sector of destination

        # Welfare
        omega, (start = sv0(:omega, nothing))                       # Cobb-Douglas utility (welfare indicator)
    end

    # -------------------------------------------------------------------------
    # STRUCTURALLY ZERO / INDETERMINATE CELLS
    #
    # A number of cells are forced to exactly zero by an equation for ANY value of
    # every other parameter or variable -- not merely small, or zero only at this
    # particular shock -- yet were given the same >= 1e-6 lower bound as every other
    # variable above. That makes the feasibility problem inconsistent by construction
    # (the equation demands 0, the bound forbids anything below 1e-6): Ipopt tolerates
    # the resulting ~1e-6 constraint violation near the base-year start point (see
    # test/reference.json, where these cells sit at ~9.9e-7, just under the bound), but
    # it is the root cause of the spurious LOCALLY_INFEASIBLE reports on shocked solves
    # documented in `solve`'s docstring and test/robustness_grid.jl (the concentration
    # of residual failures in cimint/bienscap/biensint/services was the tell -- those
    # sectors carry several of the cells below). Each such cell is fixed directly below
    # with `fix(...; force = true)` instead of a lower bound, and the now-redundant
    # equation that used to "force" it is dropped, or its index set narrowed, for that
    # cell only (never for a non-zero cell -- see the comments at `gdeq`, `dutydef`,
    # `closurem`, `closuree`, `cdeq`, `dsteq`, `ieq`, `profitmax`/`closurlp`/`closurla`
    # below, where these used to live):
    #
    #   gd[i] for i with gles[i] == 0 -- gdeq: gd[i] == gles[i]*gdtot == 0 regardless
    #       of gdtot (10 of 11 sectors in this data; only `publiques` has gles > 0)
    #   cd[i] for i with cles[i] == 0 -- cdeq: p[i]*cd[i] == cles[i]*(1-mps)*y == 0
    #       regardless of p/mps/y (sylvicult, cimint, bienscap in this data; `obj`
    #       already excludes these via its own `cles[i] > 0.0` filter -- cdeq did not)
    #   dst[i] for i with dstr[i] == 0 -- dsteq: dst[i] == dstr[i]*xd[i] == 0
    #       regardless of xd (cimint, construct, services, publiques in this data)
    #   id[i] for i with row i of `imat` entirely zero -- ieq: id[i] ==
    #       sum(imat[i,j]*dk[j] for j in SEC) == 0 regardless of dk, because every
    #       coefficient in the sum is 0 (only agsubsist/bienscap/construct actually
    #       supply capital goods in this economy's `imat`; the other 8 sectors' rows
    #       are all-zero)
    #   labd[i,l] for (i,l) with alphl[l,i] == 0 -- profitmax[i,l]:
    #       wa[l]*wdist[i,l]*labd[i,l] == xd[i]*pva[i]*alphl[l,i] has wdist[i,l] == 0
    #       too at every such cell in this data, so BOTH sides are identically zero
    #       regardless of wa/xd/pva/labd -- a zero-Jacobian-row equation carrying no
    #       information, worse than merely redundant (this is exactly what closurlp/
    #       closurla used to patch for their one cell each; profitmax's own zero row
    #       means they were the only real source of the closure, just as
    #       inconsistently bounded as everywhere else here). Two cells in this data:
    #       labd[:publiques,:rural], labd[:agsubsist,:urbanskil] (zero base employment)
    #   duty -- dutydef: duty == sum(te[it]*e[it]*pe[it] for it in IT) == 0 regardless
    #       of e/pe, since te[it] == 0 for every it (no export duties in this model)
    #   e[itn], m[itn] for itn in ITN -- closuree/closurem: non-traded sectors have no
    #       exports or imports by construction, for any parameter value
    #
    # Fixing e[itn]/m[itn] leaves four cells indeterminate rather than zero (they now
    # multiply an identically-zero quantity in the one equation that mentions them, or
    # appear in no equation at all) -- reference.json shows exactly this: pe[itn] lands
    # on an arbitrary small value that differs between baseline and sim1, and
    # pm[itn]/pwe[itn]/tm[itn] (which appear in NO equation for itn -- pmdef/
    # absorption1/costmin/pedef/tariffdef/closuret/caeq all index or sum over IT only)
    # run away to ~1e5. Per cell, pin it to its base value (dropping the cell from the
    # variable set is the alternative -- rejected here so every Simulation field keeps
    # its full SEC-indexed shape, which callers rely on):
    #
    #   pe[itn]  -- appears only in sales[itn]: px*xd == pd*xxd + pe[itn]*e[itn]
    #   pm[itn], pwe[itn], tm[itn] -- appear in no equation at all
    #
    # NOTE: tm[services] is a different case, deliberately NOT touched here: tm0
    # happens to be 0 in this data's base year, but tm[services] is a live,
    # shockable policy instrument (tariff scenarios on `services` set tm0[services] >
    # 0), not a structural zero for any parameter value -- fixing it would silently
    # break every such shock. Likewise id0/dst0/cd0 having zero cells beyond the ones
    # above (e.g. id0[:construct], id0[:bienscap]) is a base-year coincidence, not a
    # structural forcing -- imat's row for those sectors is NOT all-zero, so id there
    # can legitimately move away from 0 under a shock, and is left untouched.
    # -------------------------------------------------------------------------
    zero_gd    = [i for i in SEC if gles[i] == 0.0]
    nonzero_gd = [i for i in SEC if gles[i] != 0.0]
    zero_cd    = [i for i in SEC if cles[i] == 0.0]
    nonzero_cd = [i for i in SEC if cles[i] != 0.0]
    zero_dst    = [i for i in SEC if dstr[i] == 0.0]
    nonzero_dst = [i for i in SEC if dstr[i] != 0.0]
    zero_id    = [i for i in SEC if all(imat[i, j] == 0.0 for j in SEC)]
    nonzero_id = [i for i in SEC if !(i in zero_id)]
    zero_labd  = [(i, l) for i in SEC, l in LC if alphl[l, i] == 0.0]

    for i in zero_gd;  fix(gd[i],  0.0; force = true); end
    for i in zero_cd;  fix(cd[i],  0.0; force = true); end
    for i in zero_dst; fix(dst[i], 0.0; force = true); end
    for i in zero_id;  fix(id[i],  0.0; force = true); end
    for (i, l) in zero_labd
        fix(labd[i, l], 0.0; force = true)
    end
    fix(duty, 0.0; force = true)
    for itn in ITN
        fix(e[itn],   0.0;       force = true)
        fix(m[itn],   0.0;       force = true)
        fix(pe[itn],  pd0[itn];  force = true)
        fix(pm[itn],  pm0[itn];  force = true)
        fix(pwe[itn], pwe0[itn]; force = true)
        fix(tm[itn],  0.0;       force = true)
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
        # The product skips l with alphl[l,i] == 0 -- mathematically a no-op (x^0 == 1
        # for any x, so the value is unchanged), but symbolically necessary: at such a
        # cell labd[i,l] is fixed to exactly 0 above, and automatic differentiation of
        # x^0 is p*x^(p-1) = 0*x^(-1), which evaluates to 0*Inf = NaN at x == 0 even
        # though the function value itself (1) and true derivative (0) are fine.
        activity[i in SEC],       xd[i] == ad[i] * prod(labd[i,l]^alphl[l,i] for l in LC if alphl[l,i] != 0.0) *
                                            k[i]^(1 - sum(alphl[l,i] for l in LC))

        # Profit maximisation (Shephard's lemma / FOC for labour):
        #   wage × wdist × employment = value-added price × output × labour share
        # Restricted to (i,l) with alphl[l,i] != 0 -- at every zero-employment cell in
        # this data wdist[i,l] == 0 too, so both sides are identically zero regardless
        # of wa/xd/pva/labd (a zero-Jacobian-row equation, not merely redundant); those
        # cells' labd are fixed directly above instead.
        profitmax[i in SEC, l in LC; alphl[l,i] != 0.0], wa[l]*wdist[i,l]*labd[i,l] == xd[i]*pva[i]*alphl[l,i]

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

        # Inventory investment: fixed proportion of sectoral gross output. Restricted
        # to sectors with dstr[i] != 0 -- for dstr[i] == 0 this equation would force
        # dst[i] == 0 for any xd, so those cells are fixed directly above instead.
        dsteq[i in nonzero_dst],  dst[i] == dstr[i]*xd[i]

        # Private consumption: linear Engel curves, share cles[i] of disposable
        # income -- income NET of the direct tax, `(1 - td)*y`. Restricted to sectors
        # with cles[i] != 0 -- for cles[i] == 0 this equation would force cd[i] == 0
        # for any p/mps/td/y, so those cells are fixed directly above instead (mirrors
        # the `cles[i] > 0.0` filter `obj` already uses for the same reason).
        cdeq[i in nonzero_cd],    p[i]*cd[i] == cles[i]*(1 - mps)*(1 - td)*y

        # Private GDP: aggregate value added minus economy-wide depreciation
        gdp,                      y == sum(pva[i]*xd[i] for i in SEC) - deprecia

        # Household savings: fixed marginal propensity to save applied to DISPOSABLE
        # income (income net of the direct tax). Free to be negative (dissaving).
        hhsaveq,                  hhsav == mps*(1 - td)*y

        # Government revenue: sum of all tax instruments, including the direct tax on
        # household income (`td*y`) -- the counterpart of the `(1 - td)` the household
        # keeps in `cdeq`/`hhsaveq`, so the two sides of the transfer always net out.
        greq,                     gr == tariff + duty + indtax + td*y

        # Government budget: revenue finances consumption plus budget surplus/deficit
        gruse,                    gr == sum(p[i]*gd[i] for i in SEC) + govsav

        # Government consumption: each commodity gets a fixed share of total volume.
        # Restricted to sectors with gles[i] != 0 -- for gles[i] == 0 this equation
        # would force gd[i] == 0 for any gdtot, so those cells are fixed directly
        # above instead (see the "STRUCTURALLY ZERO" block).
        gdeq[i in nonzero_gd],    gd[i] == gles[i]*gdtot

        # Tariff revenue: ad-valorem tariffs applied to import values at world prices
        tariffdef,                tariff == sum(tm[it]*m[it]*pwm[it] for it in IT)*er

        # Indirect tax revenue: ad-valorem production taxes on gross output values
        indtaxdef,                indtax == sum(itax[i]*px[i]*xd[i] for i in SEC)

        # Export duty revenue: always 0, since te[it] == 0 for every it (no export
        # duties in this model version) -- `duty` is fixed directly above instead of
        # via this equation (see the "STRUCTURALLY ZERO" block); keeping both would
        # make `duty` doubly-determined (bound and equation both forcing the same 0).

        # Aggregate depreciation expenditure (capital consumption)
        depreq,                   deprecia == sum(depr[i]*pk[i]*k[i] for i in SEC)

        # Savings-investment balance: total savings = total investment
        totsav,                   savings == hhsav + govsav + deprecia + fsav*er

        # Investment by sector of destination: savings allocated by fixed shares kio,
        # net of inventory investment (inventories are not part of capital formation)
        prodinv[i in SEC],        pk[i]*dk[i] == kio[i]*savings - kio[i]*sum(dst[j]*p[j] for j in SEC)

        # Investment demand by sector of origin: derived from capital composition
        # matrix. Restricted to sectors whose `imat` row is not entirely zero -- for
        # an all-zero row this equation would force id[i] == 0 for any dk (only
        # agsubsist/bienscap/construct actually supply capital goods in this data),
        # so those cells are fixed directly above instead.
        ieq[i in nonzero_id],     id[i] == sum(imat[i,j]*dk[j] for j in SEC)

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
        closuretd,            td      == td0         # Household direct-tax rate fixed
        closureg,             gdtot   == gdtot0      # Government spending volume fixed
        # Non-traded: no imports/no exports -- m[itn]/e[itn] are fixed directly above
        # instead of via an equation here (see the "STRUCTURALLY ZERO" block); keeping
        # both would make them doubly-determined (bound and equation both forcing 0).
        # No rural labour in public sector / no skilled labour in subsistence farming
        # (closurlp/closurla) -- labd[:publiques,:rural]/labd[:agsubsist,:urbanskil]
        # are fixed directly above instead (see the "STRUCTURALLY ZERO" block: this is
        # the (i,l) pair where profitmax itself has a zero Jacobian row).
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
