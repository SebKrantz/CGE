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

"Maximum relative difference between a Simulation and a reference Dict (as read by
`parse_json_numeric`), scanning every field (scalars and every element of every
indexed field)."
function max_relative_diff(sim::CGECameroon.Simulation, ref::Dict{String,Any})
    worst = 0.0
    for f in fieldnames(CGECameroon.Simulation)
        v = getfield(sim, f)
        rv = ref[string(f)]
        if v isa JuMP.Containers.DenseAxisArray
            ax = axes(v)
            if length(ax) == 1
                for idx in ax[1]
                    worst = max(worst, _relerr(v[idx], rv[string(idx)]))
                end
            else
                for i1 in ax[1], i2 in ax[2]
                    worst = max(worst, _relerr(v[i1, i2], rv[string(i1) * "|" * string(i2)]))
                end
            end
        else
            worst = max(worst, _relerr(Float64(v), Float64(rv)))
        end
    end
    return worst
end
_relerr(a, b) = abs(a - b) / max(abs(b), 1e-12)

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
            @test max_relative_diff(baseline, REFERENCE["baseline"]) < 1e-6
        end
    end

    @testset "sim1 matches the pre-refactor reference" begin
        shocked = CGECameroon.with_shocks(params,
            Dict{Symbol,Any}(:k0 => Dict(:agsubsist => params.k0[:agsubsist] * 1.10)))
        # with_shocks must not mutate its input.
        @test params.k0[:agsubsist] != shocked.k0[:agsubsist]

        sim1, status, converged, iterations, elapsed = CGECameroon.solve(shocked)
        @test converged
        @test max_relative_diff(sim1, REFERENCE["sim1"]) < 1e-6
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
    # Params (every variable starting from the base-year `*0` data, regardless of the shock) can
    # return Ipopt's LOCALLY_INFEASIBLE even though the shocked equations are, in fact,
    # satisfiable: e.g. a +5pp import tariff on services, agexpind, or indalim. Warm-starting
    # from a previously-solved Simulation (`solve(...; start = baseline)`) fixes this reliably
    # (see cge.jl's `solve` docstring) -- this testset checks that fix directly against the
    # exact cases the original bug report found, plus one control sector (bienscons) that
    # converges even cold.
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
        # services/agexpind/indalim must reproduce the original cold-start bug (bienscons and
        # foreign savings converge even cold, so this counts only the known-broken cases).
        @test count(r -> !r.cold_converged, rows) == 3
    end
end
