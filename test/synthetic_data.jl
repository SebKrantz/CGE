# A small synthetic base-year dataset for the N-sector loader/model: 3 sectors (2 traded,
# 1 non-traded) and 2 labour categories, with a base-year TRADE SURPLUS (fsav0 < 0), a
# GOVERNMENT DEFICIT (govsav0 < 0) and a positive household DIRECT TAX (td0 = 5%) --
# i.e. exactly the three things the pre-`n-sector` model could not represent (`fsav`,
# `govsav` bounded below at 1e-6; no direct tax at all).
#
# Everything is derived from a handful of primitives so the SAM balances by construction:
# household consumption `cd0` is the residual of each commodity row, and Walras' law then
# makes `Σ cd0` come out equal to `(1-mps0)(1-td0)y0` on its own (asserted in runtests.jl).
# See /Users/…/docs/cge_static_identities.md §§0-10 for every identity this satisfies.
#
# Used by test/runtests.jl; not loaded by cge.jl.

"Builds the synthetic dataset as a `CGECameroon.RawData` (no file involved)."
function synthetic_rawdata()
    SEC = [:agri, :manuf, :serv]
    IT  = [:agri, :manuf]          # traded
    ITN = [:serv]                  # non-traded
    LC  = [:unskilled, :skilled]

    er = 1.0                       # SAM in "million USD": er is a pure units choice

    # --- technology and taxes ------------------------------------------------------
    xd0  = Dict(:agri => 200.0, :manuf => 300.0, :serv => 500.0)
    io   = Dict{Tuple{Symbol,Symbol},Float64}(          # io[supplier, user]
        (:agri, :agri) => 0.10, (:agri, :manuf) => 0.15, (:agri, :serv) => 0.02,
        (:manuf,:agri) => 0.08, (:manuf,:manuf) => 0.20, (:manuf,:serv) => 0.10,
        (:serv, :agri) => 0.05, (:serv, :manuf) => 0.10, (:serv, :serv) => 0.15)
    itax = Dict(:agri => 0.02, :manuf => 0.05, :serv => 0.03)
    pva0 = Dict(i => 1.0 - sum(io[j, i] for j in SEC) - itax[i] for i in SEC)  # 0.75/0.50/0.70
    va   = Dict(i => pva0[i]*xd0[i] for i in SEC)                              # 150/150/350

    k0   = Dict(:agri => 400.0, :manuf => 500.0, :serv => 600.0)
    depr = Dict(:agri => 0.02,  :manuf => 0.03,  :serv => 0.025)
    deprecia0 = sum(depr[i]*k0[i] for i in SEC)                                # 38
    y0   = sum(va[i] for i in SEC) - deprecia0                                 # 612

    # --- external account: a base-year TRADE SURPLUS -------------------------------
    tm0 = Dict(:agri => 0.10, :manuf => 0.20, :serv => 0.0)
    cif = Dict(:agri => 20.0, :manuf => 50.0, :serv => 0.0)                    # cif imports
    m0  = Dict(i => cif[i]*(1 + tm0[i]) for i in SEC)                          # duty-paid
    e0  = Dict(:agri => 40.0, :manuf => 70.0, :serv => 0.0)
    fsav0 = sum(cif[i] for i in SEC) - sum(e0[i] for i in SEC)                 # -40

    # --- government: a direct tax, and a base-year DEFICIT --------------------------
    td0     = 0.05
    tariff0 = sum(m0[i]*tm0[i]/(1 + tm0[i]) for i in IT)                       # 12
    indtax0 = sum(itax[i]*xd0[i] for i in SEC)                                 # 34
    gr0     = tariff0 + indtax0 + td0*y0                                       # 76.6
    gdtot0  = 100.0
    gles    = Dict(:agri => 0.05, :manuf => 0.15, :serv => 0.80)               # sums to 1
    gd0     = Dict(i => gles[i]*gdtot0 for i in SEC)
    govsav0 = gr0 - sum(gd0[i] for i in SEC)                                   # -23.4

    # --- household and savings-investment ------------------------------------------
    mps0     = 0.10
    cdtot0   = (1 - mps0)*(1 - td0)*y0                                         # 523.26
    hhsav0   = mps0*(1 - td0)*y0                                               # 58.14
    savings0 = hhsav0 + govsav0 + deprecia0 + fsav0*er                         # 32.74

    dstr = Dict(:agri => 0.01, :manuf => 0.005, :serv => 0.0)
    dst0 = Dict(i => dstr[i]*xd0[i] for i in SEC)
    netinv = savings0 - sum(dst0[i] for i in SEC)                              # 29.24
    kio  = Dict(:agri => 0.2, :manuf => 0.3, :serv => 0.5)                     # all > 0, sums to 1
    dk0  = Dict(i => kio[i]*netinv for i in SEC)
    # Capital goods come from manufacturing and services only: agri's row is all-zero, so
    # the model fixes id[:agri] == 0 (exercises the structurally-zero `ieq` path). Columns
    # sum to 1, hence pk == 1 at base.
    imat = Dict{Tuple{Symbol,Symbol},Float64}()
    for j in SEC
        imat[:agri,  j] = 0.0
        imat[:manuf, j] = 0.6
        imat[:serv,  j] = 0.4
    end
    id0 = Dict(i => sum(imat[i, j]*dk0[j] for j in SEC) for i in SEC)

    # --- commodity balance closes on household consumption -------------------------
    int0 = Dict(i => sum(io[i, j]*xd0[j] for j in SEC) for i in SEC)
    x0   = Dict(i => xd0[i] - e0[i] + m0[i] for i in SEC)
    cd0  = Dict(i => x0[i] - int0[i] - gd0[i] - id0[i] - dst0[i] for i in SEC)
    cles = Dict(i => cd0[i]/cdtot0 for i in SEC)
    cles[SEC[end]] = 1.0 - sum(cles[i] for i in SEC[1:end-1])   # sums to exactly 1 in Float64

    # --- factor markets (the DATA.md "wa0 = 1, wdist = 1, xle = wage bill" convention) --
    xle = Dict{Tuple{Symbol,Symbol},Float64}(
        (:agri,  :unskilled) =>  70.0, (:agri,  :skilled) =>  20.0,
        (:manuf, :unskilled) =>  40.0, (:manuf, :skilled) =>  35.0,
        (:serv,  :unskilled) => 100.0, (:serv,  :skilled) => 145.0)
    wdist = Dict{Tuple{Symbol,Symbol},Float64}((i, l) => 1.0 for i in SEC, l in LC)
    wa0   = Dict{Symbol,Float64}(l => 1.0 for l in LC)

    # --- the `miscellaneous` rows ---------------------------------------------------
    rows = Dict(
        :m0 => m0, :e0 => e0, :xd0 => xd0, :k => k0, :depr => depr,
        :rhoc => Dict(:agri => 2.0, :manuf => 1.5, :serv => 0.4),   # Armington sigma (!= 1)
        :rhot => Dict(:agri => 3.0, :manuf => 2.0, :serv => 0.4),   # CET sigma
        :eta  => Dict(i => 4.0 for i in SEC),                       # export-demand elasticity
        :pd0  => Dict(i => 1.0 for i in SEC),                       # base prices normalised to 1
        :tm0 => tm0, :itax => itax, :cles => cles, :gles => gles, :kio => kio,
        :dstr => dstr, :dst => dst0, :id => id0)
    zz = Dict{Tuple{Symbol,Symbol},Float64}()
    for r in CGECameroon.MISC_ROWS, i in SEC
        zz[r, i] = rows[r][i]
    end

    # `tariff0` is deliberately left out of `scalars`, so `calibrate` derives it.
    scalars = Dict{Symbol,Float64}(:er => er, :gr0 => gr0, :gdtot0 => gdtot0,
                                   :cdtot0 => cdtot0, :fsav0 => fsav0,
                                   :mps0 => mps0, :td0 => td0)

    raw = CGECameroon.RawData(SEC, IT, ITN, LC, io, imat, wdist, xle, zz, wa0, scalars)
    expected = (; y0, deprecia0, tariff0, indtax0, gr0, govsav0, hhsav0, cdtot0,
                  savings0, fsav0, td0, mps0, gdtot0)
    return raw, expected
end
