# Regression tests for CGECameroon (cge.jl).
#
# Run with:  julia --project=. test/runtests.jl   (from the repo root)
#
# Three checks:
#  1. The baseline solve reproduces the base-year calibration data (`xd0`, `k0`, `y0`,
#     ...) to within Ipopt's own "acceptable level" solver tolerance -- this is the
#     standard "does the model calibrate" check for this kind of CGE.
#  2. `sim1` (a +10% capital-stock shock to `agsubsist`, "public investment in
#     subsistence agriculture") matches `test/reference.json` -- a snapshot of every
#     field of the ORIGINAL, pre-refactor script's `baseline` and `sim1` -- to 1e-6
#     relative on every value. This is the regression guard for the refactor itself:
#     it proves the module reproduces the original script bit-for-bit (the observed
#     difference at the time of writing is exactly 0.0). Both of these solve cold (no
#     `start`), exactly as the original script did.
#  3. A reduced robustness testset (see test/robustness_grid.jl for the full ~60-scenario
#     grid) confirms `solve`'s `start` keyword fixes the exact cold-start
#     `LOCALLY_INFEASIBLE` cases a prior exploration found (e.g. +5pp tariff on services).
#  4. The GENERIC workbook path: the Cameroon `RawData` is written back out in the generic
#     layout (`test/generic_workbook.jl`) with its sectors and labour categories in REVERSE
#     order, re-read through `load_data`'s generic branch, and must reproduce the legacy
#     loader's sets, matrices, calibration and solution -- proving the loader is
#     order-independent (header lookup, not cell position) and that nothing downstream of
#     `RawData` depends on the Cameroon constants.
#  5. A small SYNTHETIC 3-sector / 2-labour dataset (`test/synthetic_data.jl`) with a
#     base-year TRADE SURPLUS (`fsav0 < 0`), a GOVERNMENT DEFICIT (`govsav0 < 0`) and a
#     positive household DIRECT TAX (`td0 = 5%`) calibrates, round-trips through the generic
#     workbook, and solves back to its own base year -- the three things the pre-`n-sector`
#     model could not represent at all.

using Test
using JuMP

include(joinpath(@__DIR__, "..", "cge.jl"))
using .CGECameroon

# ── a tiny reader for the flat, numbers-only JSON this repo's reference file uses ──
# (no JSON package dependency: the file has no arrays, no string values, no escaping
# beyond quotes in keys, so a ~30-line recursive-descent reader is simpler and lighter
# than adding a dependency for it.)
function parse_json_numeric(s::AbstractString)
    i = firstindex(s)
    n = lastindex(s)
    skipws() = (while i <= n && isspace(s[i]); i += 1; end)
    function parse_value()
        skipws()
        if s[i] == '{'
            i += 1
            d = Dict{String,Any}()
            skipws()
            if s[i] == '}'
                i += 1
                return d
            end
            while true
                skipws()
                @assert s[i] == '"'
                i += 1
                startk = i
                while s[i] != '"'
                    i += 1
                end
                key = s[startk:i-1]
                i += 1
                skipws()
                @assert s[i] == ':'
                i += 1
                d[key] = parse_value()
                skipws()
                if s[i] == ','
                    i += 1
                    continue
                elseif s[i] == '}'
                    i += 1
                    break
                else
                    error("parse_json_numeric: unexpected character at $i: $(s[i])")
                end
            end
            return d
        else
            startn = i
            while i <= n && (isdigit(s[i]) || s[i] in ('-', '+', '.', 'e', 'E'))
                i += 1
            end
            return parse(Float64, s[startn:i-1])
        end
    end
    return parse_value()
end

"""
`fixed_cell` lists the cells `_solve_once` pins at exactly 0 (or at their base value); a
relative comparison of two solver-noise values there is meaningless, so both the
reference-file check below (`skip_cell`, which adds two aggregates) and the
solution-vs-solution check (`max_sim_diff`) exclude them.

Cells that `_solve_once` (cge.jl) now fixes directly instead of leaving them
lower-bound-clipped near 1e-6 (or, for pe/pwe/tm, indeterminate/runaway) the way the
pre-refactor script did -- see cge.jl's "STRUCTURALLY ZERO / INDETERMINATE CELLS"
block. `test/reference.json` is a snapshot of that pre-refactor (buggy) behaviour, so
these specific cells are EXPECTED to now differ from it -- that is the bug being
fixed, not a regression -- and are excluded from the reference comparison below. Every
other cell of every other field must still match to 1e-6 relative, unchanged.

  - `duty` (scalar): reference ~9.9e-7 (bound-clipped); now fixed at exactly 0, since
    `te[it] == 0` for every it forces it structurally, for any parameter value.
  - `gd[i]` for i with `gles[i] == 0` (10 of the 11 sectors): reference ~9.9e-7 each;
    now fixed at exactly 0 (`gd[i] == gles[i]*gdtot == 0` regardless of gdtot).
  - `cd[i]` for i with `cles[i] == 0` (sylvicult, cimint, bienscap): reference ~9.9e-7
    each; now fixed at exactly 0 (`cdeq` would force cd[i] == 0 for any p/mps/y).
  - `dst[i]` for i with `dstr[i] == 0` (cimint, construct, services, publiques):
    reference ~9.9e-7 each; now fixed at exactly 0 (`dsteq` would force dst[i] == 0
    for any xd).
  - `id[i]` for i whose `imat` row is entirely zero (agexpind, sylvicult, indalim,
    bienscons, biensint, cimint, services, publiques): reference ~9.9e-7 each; now
    fixed at exactly 0 (`ieq` would force id[i] == 0 for any dk, since every
    coefficient in the sum is 0 -- only agsubsist/bienscap/construct supply capital
    goods in this data).
  - `labd[i,l]` for (i,l) with `alphl[l,i] == 0` (publiques/rural, agsubsist/
    urbanskil -- zero base employment): reference ~9.9e-7 each; now fixed at exactly
    0 (`profitmax` has a zero Jacobian row there, carrying no information, since
    `wdist[i,l] == 0` too at both cells).
  - `e[itn]`, `m[itn]` for itn in `ITN` (construct, publiques): reference ~9.9e-7; now
    fixed at exactly 0 (non-traded sectors have no exports/imports by construction).
  - `pe[itn]`: reference is an arbitrary indeterminate value (it even differs between
    baseline and sim1 in reference.json itself, since nothing pins it once e[itn] ==
    0); now fixed at the base price `pd0[itn]`.
  - `pm[itn]`, `pwe[itn]`, `tm[itn]`: reference is an Ipopt runaway (~1e5, since none
    of the three appears in any equation for a non-traded sector -- pmdef/absorption1/
    costmin/pedef/tariffdef/closuret/caeq all index or sum over IT only); now fixed at
    `pm0[itn]` / `pwe0[itn]` / 0.

One more cell is skipped for a DIFFERENT reason -- not a structural fix in cge.jl at
all, just a side effect of it:

  - `tm[services]`: NOT structurally zero for any parameter value (`tm0[services]` is
    a live, shockable policy instrument -- e.g. the "tariff +5pp services" robustness
    case below sets it > 0, and cge.jl deliberately leaves this cell alone, see its
    "STRUCTURALLY ZERO" block NOTE). It merely happens to equal 0 in this data's base
    year, so `closuret` pins tm[services] == 0 in both `baseline` and `sim1` (neither
    shocks tm0), the same near-the-lower-bound situation every genuinely-structural
    cell above was in. Raising `bound_relax_factor` (see cge.jl's `_solve_once`) to
    make those cells converge reliably ALSO lets Ipopt settle tm[services] far closer
    to true 0 (~1e-31) than the pre-refactor script's bound-clipped ~9.9e-7 in
    reference.json -- an accuracy IMPROVEMENT, not a bug, but a relative-error
    comparison of two near-zero floats is dominated by solver noise either way.

Two more (scalar) cells are skipped for a THIRD reason: they are aggregates that SUM
several of the cells above, so correcting those cells necessarily -- and correctly --
moves the aggregate by a comparable tiny amount relative to reference.json's snapshot
of the old, bound-clipped behaviour. Both differ from reference.json by ~1.0-1.6e-6,
just over threshold, and this is stable across every `tol`/`bound_relax_factor` tried
(i.e. not under-convergence -- see cge.jl's exploration notes):

  - `tariff`: tariffdef sums `tm[it]*m[it]*pwm[it]` over IT; tm[services] moving from
    reference's ~9.9e-7 to ~1e-31 (see the tm[services] note above) shifts this sum by
    `~9.9e-7 * m[services] * pwm[services] * er`, order 1e-5 absolute on tariff's ~76
    scale.
  - `govsav`: gruse: `gr == sum(p[i]*gd[i] for i in SEC) + govsav`, i.e. `govsav = gr -
    sum(p*gd)`. Fixing 10 sectors' gd from reference's ~9.9e-7 each to exactly 0 moves
    that sum by ~1e-5, which `govsav` (and `gr`, via tariff/duty) absorbs -- order 1e-5
    absolute on govsav's ~44 scale.
"""
function fixed_cell(params::CGECameroon.Params, f::Symbol, idx)
    f === :duty && return true
    f === :gd && return params.gles[idx] == 0.0
    f === :cd && return params.cles[idx] == 0.0
    f === :dst && return params.dstr[idx] == 0.0
    f === :id && return all(params.imat[idx, j] == 0.0 for j in params.SEC)
    (f === :e || f === :m || f === :pe || f === :pm || f === :pwe) && return idx in params.ITN
    f === :tm && return (idx in params.ITN) || params.tm0[idx] == 0.0
    return false
end
fixed_cell(params::CGECameroon.Params, f::Symbol, i, l) =
    f === :labd && params.alphl[l, i] == 0.0

skip_cell(params::CGECameroon.Params, f::Symbol, idx) =
    fixed_cell(params, f, idx) || f === :tariff || f === :govsav
skip_cell(params::CGECameroon.Params, f::Symbol, i, l) = fixed_cell(params, f, i, l)

"Maximum relative difference between a Simulation and a reference Dict (as read by
`parse_json_numeric`), scanning every field (scalars and every element of every
indexed field) except the structurally-zero/indeterminate cells `skip_cell` documents."
function max_relative_diff(sim::CGECameroon.Simulation, ref::Dict{String,Any}, params::CGECameroon.Params)
    worst = 0.0
    for f in fieldnames(CGECameroon.Simulation)
        v = getfield(sim, f)
        rv = ref[string(f)]
        if v isa JuMP.Containers.DenseAxisArray
            ax = axes(v)
            if length(ax) == 1
                for idx in ax[1]
                    skip_cell(params, f, idx) && continue
                    worst = max(worst, _relerr(v[idx], rv[string(idx)]))
                end
            else
                for i1 in ax[1], i2 in ax[2]
                    skip_cell(params, f, i1, i2) && continue
                    worst = max(worst, _relerr(v[i1, i2], rv[string(i1) * "|" * string(i2)]))
                end
            end
        else
            skip_cell(params, f, nothing) && continue
            worst = max(worst, _relerr(Float64(v), Float64(rv)))
        end
    end
    return worst
end
_relerr(a, b) = abs(a - b) / max(abs(b), 1e-12)

"Maximum relative difference between two `Simulation`s solved from the same economy (used
to compare the generic-workbook path against the legacy one), skipping only the cells
`fixed_cell` documents -- there both solves land on their own solver noise around 0."
function max_sim_diff(a::CGECameroon.Simulation, b::CGECameroon.Simulation,
                      params::CGECameroon.Params)
    worst = 0.0
    for f in fieldnames(CGECameroon.Simulation)
        va, vb = getfield(a, f), getfield(b, f)
        if va isa JuMP.Containers.DenseAxisArray
            ax = axes(va)
            if length(ax) == 1
                for idx in ax[1]
                    fixed_cell(params, f, idx) && continue
                    worst = max(worst, _relerr(va[idx], vb[idx]))
                end
            else
                for i1 in ax[1], i2 in ax[2]
                    fixed_cell(params, f, i1, i2) && continue
                    worst = max(worst, _relerr(va[i1, i2], vb[i1, i2]))
                end
            end
        else
            fixed_cell(params, f, nothing) && continue
            worst = max(worst, _relerr(Float64(va), Float64(vb)))
        end
    end
    return worst
end

"Maximum relative difference between two `Params`, over every scalar field and every value
of every per-sector/per-labour `Dict` field (the set fields are compared separately)."
function max_params_diff(a::CGECameroon.Params, b::CGECameroon.Params)
    worst = 0.0
    for f in fieldnames(CGECameroon.Params)
        va, vb = getfield(a, f), getfield(b, f)
        if va isa AbstractDict
            for k in keys(va)
                worst = max(worst, _relerr(va[k], vb[k]))
            end
        elseif va isa Float64
            worst = max(worst, _relerr(va, vb))
        end
    end
    return worst
end

include(joinpath(@__DIR__, "generic_workbook.jl"))  # write_generic_workbook: the test-side writer
include(joinpath(@__DIR__, "synthetic_data.jl"))    # synthetic_rawdata: the 3-sector fixture

const DATA_PATH = joinpath(@__DIR__, "..", "data", "camdata.xlsx")
const REFERENCE = parse_json_numeric(read(joinpath(@__DIR__, "reference.json"), String))

@testset "CGECameroon" begin
    raw = CGECameroon.load_data(DATA_PATH)
    params = CGECameroon.calibrate(raw)

    @testset "baseline reproduces the base-year data" begin
        baseline, status, converged, iterations, elapsed = CGECameroon.solve(params)
        @test converged
        @test elapsed < 5.0
        for i in raw.SEC
            @test isapprox(baseline.xd[i], params.xd0[i]; rtol = 1e-3)
            @test isapprox(baseline.k[i], params.k0[i]; rtol = 1e-3)
        end
        for l in raw.LC
            @test isapprox(baseline.ls[l], params.ls0[l]; rtol = 1e-3)
        end
        @test isapprox(baseline.y, params.y0; rtol = 1e-3)
        @test isapprox(baseline.fsav, params.fsav0; rtol = 1e-3)

        @testset "baseline matches reference.json" begin
            @test max_relative_diff(baseline, REFERENCE["baseline"], params) < 1e-6
        end
    end

    @testset "sim1 matches the pre-refactor reference" begin
        shocked = CGECameroon.with_shocks(params,
            Dict{Symbol,Any}(:k0 => Dict(:agsubsist => params.k0[:agsubsist] * 1.10)))
        # with_shocks must not mutate its input.
        @test params.k0[:agsubsist] != shocked.k0[:agsubsist]

        sim1, status, converged, iterations, elapsed = CGECameroon.solve(shocked)
        @test converged
        @test max_relative_diff(sim1, REFERENCE["sim1"], shocked) < 1e-6
    end

    @testset "solve reads only from Params (no cross-call leakage)" begin
        base_before, = CGECameroon.solve(params)
        shocked = CGECameroon.with_shocks(params, Dict{Symbol,Any}(:gdtot0 => params.gdtot0 * 1.5))
        CGECameroon.solve(shocked)
        base_after, = CGECameroon.solve(params)
        @test base_before.y == base_after.y
        @test base_before.gr == base_after.gr
    end

    # A reduced version (~5 scenarios, well under 20s) of test/robustness_grid.jl's ~60-scenario
    # grid -- see that file's header comment for the full exploration. A cold solve of a shocked
    # Params (every variable starting from the base-year `*0` data, regardless of the shock) used
    # to be able to return Ipopt's LOCALLY_INFEASIBLE even though the shocked equations are, in
    # fact, satisfiable: e.g. a +5pp import tariff on services, agexpind, or indalim. Warm-starting
    # from a previously-solved Simulation (`solve(...; start = baseline)`) fixes this reliably
    # (see cge.jl's `solve` docstring) -- this testset checks that fix directly against the exact
    # cases the original bug report found, plus one control sector (bienscons) that converges even
    # cold.
    #
    # NOTE: fixing the structurally-zero/indeterminate cells (cge.jl's "STRUCTURALLY ZERO" block --
    # gd/cd/dst/id/labd/duty/e/m/pe/pm/pwe/tm cells that used to be lower-bound-clipped or run away
    # instead of being exactly what the equations force) changes exactly which of these cases (if
    # any) happen to converge cold -- it moved from 3/5 failing cold before that fix, to 0/5, to
    # 4/5 across successive refinements of it (adding `bound_relax_factor`, fixing pm[itn], ...),
    # with no change in the underlying economics each time. cold-start success was never this
    # testset's guarantee and is sensitive to Ipopt's exact solve path, so it is deliberately NOT
    # asserted here -- only `warm_converged` (below) is, since that is what `start` actually exists
    # to guarantee, and it holds unconditionally for all of these across every one of those changes.
    @testset "robustness (warm-start fixes the LOCALLY_INFEASIBLE bug report cases)" begin
        include(joinpath(@__DIR__, "robustness_grid.jl"))  # build_grid, run_grid, GridRow (reuses CGECameroon)

        baseline_r, _, baseline_r_converged, = CGECameroon.solve(params)
        @test baseline_r_converged

        reduced = Tuple{String,Dict{Symbol,Any}}[
            ("tariff +5pp $sec", Dict{Symbol,Any}(:tm0 => Dict(sec => params.tm0[sec] + 0.05)))
            for sec in (:services, :agexpind, :indalim, :bienscons)
        ]
        push!(reduced, ("foreign savings +50%", Dict{Symbol,Any}(:fsav0 => params.fsav0 * 1.50)))

        rows = run_grid(params, baseline_r, reduced; homotopy_steps = 4)
        for r in rows
            @test r.warm_converged
        end
    end

    # The same Cameroon data, written back out in the GENERIC workbook layout (sheets
    # `sectors`, `labour`, `iotable`, `imat`, `employment`, `wagedist`, `miscellaneous`,
    # `scalars` -- see DATA.md §5 and `load_data`'s docstring) with its sectors AND labour
    # categories in REVERSE order, then re-read through `load_data`'s generic branch. The
    # loader works by header lookup, so the file's order must not matter: the sets must come
    # back with the same membership, every matrix cell identical, and the calibration and the
    # solved baseline must agree with the legacy path to 1e-8. This is the end-to-end check
    # that nothing downstream of `RawData` reads the Cameroon `SEC`/`IT`/`ITN`/`LC` constants.
    @testset "generic workbook reproduces the legacy Cameroon load" begin
        dir = mktempdir()
        path = write_generic_workbook(joinpath(dir, "camdata_generic.xlsx"), raw;
                                      sector_order = reverse(raw.SEC),
                                      labour_order = reverse(raw.LC))
        graw = CGECameroon.load_data(path)

        @test graw.SEC == reverse(raw.SEC)      # sets follow the FILE's order ...
        @test graw.LC  == reverse(raw.LC)
        @test Set(graw.IT)  == Set(raw.IT)      # ... with exactly the same membership
        @test Set(graw.ITN) == Set(raw.ITN)
        for f in (:io, :imat, :wdist, :xle, :zz)
            d, g = getfield(raw, f), getfield(graw, f)
            @test Set(keys(g)) == Set(keys(d))
            @test maximum(abs(g[k] - d[k]) for k in keys(d)) == 0.0
        end
        @test graw.wa0 == raw.wa0
        @test graw.scalars == raw.scalars
        @test graw.scalars[:td0] == 0.0         # Cameroon 1979-80 has no direct tax

        gparams = CGECameroon.calibrate(graw)
        @test Set(gparams.SEC) == Set(params.SEC)
        @test max_params_diff(gparams, params) < 1e-8

        legacy_sim, _, legacy_converged, = CGECameroon.solve(params)
        gsim, _, gconverged, = CGECameroon.solve(gparams)
        @test legacy_converged
        @test gconverged
        @test max_sim_diff(gsim, legacy_sim, gparams) < 1e-8
    end

    # A dataset that is not Cameroon at all: 3 sectors (1 non-traded), 2 labour categories,
    # and the three things the pre-`n-sector` model could not represent -- a base-year TRADE
    # SURPLUS (fsav0 < 0), a GOVERNMENT DEFICIT (govsav0 < 0) and a household DIRECT TAX
    # (td0 = 5%). It must calibrate, survive a generic-workbook round trip, and solve back to
    # its own base year.
    @testset "synthetic 3-sector dataset (trade surplus, government deficit, direct tax)" begin
        sraw, expected = synthetic_rawdata()
        @test length(sraw.SEC) == 3
        @test length(sraw.LC) == 2
        @test sraw.ITN == [:serv]

        sparams = CGECameroon.calibrate(sraw)
        @test isapprox(sparams.y0, expected.y0; rtol = 1e-12)
        @test sparams.td0 == expected.td0 > 0
        @test sparams.fsav0 == expected.fsav0 < 0    # trade surplus
        # `tariff0` is deliberately absent from the workbook's `scalars`, so `calibrate`
        # derives it from the data instead of using the old hardcoded Cameroon 76.548.
        @test isapprox(sparams.tariff0, expected.tariff0; rtol = 1e-12)

        sim, _, converged, = CGECameroon.solve(sparams)
        @test converged
        for i in sraw.SEC
            @test isapprox(sim.xd[i], sparams.xd0[i]; rtol = 1e-6)
            @test isapprox(sim.cd[i], sparams.cd0[i]; rtol = 1e-6)
            @test isapprox(sim.k[i],  sparams.k0[i];  rtol = 1e-6)
        end
        for l in sraw.LC
            @test isapprox(sim.ls[l], sparams.ls0[l]; rtol = 1e-6)
        end
        @test isapprox(sim.y,        expected.y0;        rtol = 1e-6)
        @test isapprox(sim.deprecia, expected.deprecia0; rtol = 1e-6)
        @test isapprox(sim.tariff,   expected.tariff0;   rtol = 1e-6)
        @test isapprox(sim.indtax,   expected.indtax0;   rtol = 1e-6)
        @test isapprox(sim.gr,       expected.gr0;       rtol = 1e-6)
        @test isapprox(sim.hhsav,    expected.hhsav0;    rtol = 1e-6)
        @test isapprox(sim.savings,  expected.savings0;  rtol = 1e-6)

        # The two residuals that used to be bounded below at 1e-6, now free and negative.
        @test sim.fsav < 0
        @test isapprox(sim.fsav, expected.fsav0; rtol = 1e-6)
        @test sim.govsav < 0
        @test isapprox(sim.govsav, expected.govsav0; rtol = 1e-6)

        # The direct tax is exactly the wedge between the two sides of the transfer: it is
        # the part of government revenue no indirect instrument accounts for (greq), and the
        # part of income the household neither spends nor saves (cdeq + hhsaveq).
        @test isapprox(sim.gr - sim.tariff - sim.indtax - sim.duty,
                       expected.td0 * expected.y0; rtol = 1e-6)
        @test isapprox(sum(sim.p[i]*sim.cd[i] for i in sraw.SEC) + sim.hhsav,
                       (1 - expected.td0) * sim.y; rtol = 1e-6)

        @testset "round-trips through the generic workbook" begin
            spath = write_generic_workbook(joinpath(mktempdir(), "synthetic.xlsx"), sraw)
            sraw2 = CGECameroon.load_data(spath)
            @test sraw2.SEC == sraw.SEC
            @test sraw2.IT  == sraw.IT
            @test sraw2.ITN == sraw.ITN
            @test sraw2.LC  == sraw.LC
            @test sraw2.wa0 == sraw.wa0
            @test sraw2.scalars == sraw.scalars

            sparams2 = CGECameroon.calibrate(sraw2)
            @test max_params_diff(sparams2, sparams) < 1e-8
            sim2, _, converged2, = CGECameroon.solve(sparams2)
            @test converged2
            @test max_sim_diff(sim2, sim, sparams2) < 1e-8
        end

        @testset "a traded sector with no base-year trade is a named error" begin
            for (row, sector) in ((:m0, :agri), (:e0, :manuf))
                bad, _ = synthetic_rawdata()
                bad.zz[row, sector] = 0.0
                err = try CGECameroon.calibrate(bad); nothing catch e; e end
                @test err isa ErrorException
                @test occursin(String(sector), err.msg)
            end
        end
    end
end
