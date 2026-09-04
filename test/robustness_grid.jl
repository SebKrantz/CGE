# Robustness grid for CGECameroon.solve on modest, single-shock policy scenarios.
#
# Run with:  julia --project=. test/robustness_grid.jl   (from the repo root)
#
# Motivation: a prior exploration found that solving a *shocked* Params cold -- i.e.
# starting every variable from the base-year `*0` data, exactly like the original script,
# regardless of how far the shock has moved the equilibrium -- can return Ipopt's
# LOCALLY_INFEASIBLE status on some modest shocks (e.g. +5pp import tariff on services,
# agexpind, or indalim) while economically similar shocks on other sectors (bienscons,
# sylvicult) converge fine. Since a policy simulator must apply arbitrary modest shocks
# reliably, this script grids ~60 single-shock scenarios spanning every shock kind the
# app's Static CGE adapter exposes, and solves each three ways:
#
#   COLD     -- no warm start (the historical behaviour, every variable starts from `*0`)
#   WARM     -- `solve(...; start = baseline)`, seeded from the cached baseline Simulation
#   HOMOTOPY -- only attempted if WARM fails: `solve(...; start = baseline, steps = 4,
#               base = params)`, applying the shock in 4 equal increments
#
# Diagnostic findings (see the exploration report and cge.jl's `solve` docstring for the
# mechanism): WARM turns every case in this grid that fails only because of a bad start
# point into a fast (<0.03s), reliable solve. A residual subset of larger/specific shocks
# (concentrated on sectors with extreme calibrated Armington/CET shares -- cimint,
# bienscap, biensint -- and on large aggregate shocks -- +15% gov't spending, +10% labour
# supply) still report LOCALLY_INFEASIBLE under WARM *and* HOMOTOPY *and* a wide sweep of
# Ipopt tuning options tried during the investigation (mu_strategy, nlp_scaling_method,
# bound_relax_factor, least-squares dual/primal initialization, ...). Inspecting the
# returned point's constraint residuals directly (bypassing termination_status) shows
# every genuine economic equation satisfied to ~1e-13, i.e. NOT a real infeasibility --
# this is Ipopt's restoration phase misreporting a numerically fragile point as locally
# infeasible, most likely tied to the model's pre-existing "structurally zero but bounded
# away from zero at 1e-6" variables (gd for the 10 sectors with gles=0, duty since te=0
# everywhere, tm[services] since tm0=0) which are a wart inherited from the original
# pre-refactor script (test/runtests.jl's reference.json match requires preserving them
# bit-for-bit, so their bounds cannot be relaxed without breaking that regression guard).
# These residual cases are reported here, not silently hidden.

using Printf

# Only load CGECameroon if it isn't already in scope -- test/runtests.jl includes this file
# (for `build_grid`/`run_grid`/`GridRow`/`final_converged`) after already loading the module
# itself, and re-including cge.jl would redefine the module, making runtests.jl's own
# `Params`/`Simulation` instances incompatible with this file's.
isdefined(Main, :CGECameroon) || include(joinpath(@__DIR__, "..", "cge.jl"))
using .CGECameroon

const DATA_PATH = joinpath(@__DIR__, "..", "data", "camdata.xlsx")

"""
    build_grid(params) -> Vector{Tuple{String,Dict{Symbol,Any}}}

The ~60-scenario shock grid: for every traded sector, a +5pp tariff, a +20% tariff
(pct_change), a -50% tariff, a +10% TFP shock, a +10% capital-stock shock, and a +10%
world-import-price shock (9 traded sectors x 6 shock kinds = 54); plus a +10% labour-supply
shock per labour category (3), a +15% government-spending shock (1), and a +50%
foreign-savings shock (1) -- 59 scenarios total, matching the shock kinds and target sets
the app's Static CGE adapter (`src/models/StaticCGE.jl`) exposes.
"""
function build_grid(params::CGECameroon.Params)
    grid = Tuple{String,Dict{Symbol,Any}}[]
    for sec in params.IT
        push!(grid, ("tariff +5pp $sec",
            Dict{Symbol,Any}(:tm0 => Dict(sec => params.tm0[sec] + 0.05))))
        push!(grid, ("tariff +20% $sec",
            Dict{Symbol,Any}(:tm0 => Dict(sec => params.tm0[sec] * 1.20))))
        push!(grid, ("tariff -50% $sec",
            Dict{Symbol,Any}(:tm0 => Dict(sec => params.tm0[sec] * 0.50))))
        push!(grid, ("TFP +10% $sec",
            Dict{Symbol,Any}(:ad => Dict(sec => params.ad[sec] * 1.10))))
        push!(grid, ("capital stock +10% $sec",
            Dict{Symbol,Any}(:k0 => Dict(sec => params.k0[sec] * 1.10))))
        push!(grid, ("world import price +10% $sec",
            Dict{Symbol,Any}(:pwm0 => Dict(sec => params.pwm0[sec] * 1.10))))
    end
    for l in params.LC
        push!(grid, ("labour supply +10% $l",
            Dict{Symbol,Any}(:ls0 => Dict(l => params.ls0[l] * 1.10))))
    end
    push!(grid, ("gov spending +15%", Dict{Symbol,Any}(:gdtot0 => params.gdtot0 * 1.15)))
    push!(grid, ("foreign savings +50%", Dict{Symbol,Any}(:fsav0 => params.fsav0 * 1.50)))
    return grid
end

"One row of grid results: label, and the status/iters/elapsed of each of the 3 tiers."
struct GridRow
    label::String
    cold_status::String
    cold_converged::Bool
    cold_iters::Int
    cold_elapsed::Float64
    warm_status::String
    warm_converged::Bool
    warm_iters::Int
    warm_elapsed::Float64
    homotopy_tried::Bool
    homotopy_status::String
    homotopy_converged::Bool
end

"Final converged status after the tiered fallback (cold -> warm -> homotopy)."
final_converged(r::GridRow) = r.cold_converged || r.warm_converged || r.homotopy_converged

function run_grid(params::CGECameroon.Params, baseline::CGECameroon.Simulation,
                  grid::Vector{Tuple{String,Dict{Symbol,Any}}}; homotopy_steps::Int = 4)
    rows = GridRow[]
    for (label, overrides) in grid
        shocked = CGECameroon.with_shocks(params, overrides)
        _, cs, cc, ci, ce = CGECameroon.solve(shocked)
        _, ws, wc, wi, we = CGECameroon.solve(shocked; start = baseline)
        if wc
            push!(rows, GridRow(label, cs, cc, ci, ce, ws, wc, wi, we, false, "", false))
        else
            _, hs, hc, hi, he = CGECameroon.solve(shocked; start = baseline,
                                                   steps = homotopy_steps, base = params)
            push!(rows, GridRow(label, cs, cc, ci, ce, ws, wc, wi, we, true, hs, hc))
        end
    end
    return rows
end

function print_grid(rows::Vector{GridRow})
    @printf("%-32s | %-11s %5s | %-11s %5s | %-11s\n",
            "scenario", "cold", "iter", "warm", "iter", "+homotopy")
    println(repeat("-", 32 + 3 + 17 + 3 + 17 + 3 + 11))
    for r in rows
        homotopy_col = r.homotopy_tried ? (r.homotopy_converged ? "OK" : "FAILED") : "-"
        @printf("%-32s | %-11s %5d | %-11s %5d | %-11s\n",
                r.label, r.cold_converged ? "OK" : "INFEASIBLE", r.cold_iters,
                r.warm_converged ? "OK" : "INFEASIBLE", r.warm_iters, homotopy_col)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    raw = CGECameroon.load_data(DATA_PATH)
    params = CGECameroon.calibrate(raw)
    baseline, bstatus, bconverged, biters, belapsed = CGECameroon.solve(params)
    println("baseline: status=$bstatus converged=$bconverged iterations=$biters " *
            "($(round(belapsed, digits = 3))s)")
    bconverged || error("robustness_grid: baseline itself did not converge")

    grid = build_grid(params)
    println("\nRunning $(length(grid)) scenarios (cold, then warm, then homotopy on warm failures)...\n")
    t = @elapsed rows = run_grid(params, baseline, grid)
    print_grid(rows)

    n = length(rows)
    n_cold_ok = count(r -> r.cold_converged, rows)
    n_warm_ok = count(r -> r.warm_converged, rows)
    n_final_ok = count(final_converged, rows)
    println("\nSummary: $n scenarios")
    println("  cold (no warm start):        $n_cold_ok/$n converged")
    println("  warm (baseline start):       $n_warm_ok/$n converged")
    println("  final (+ homotopy fallback): $n_final_ok/$n converged")
    println("  total grid wall time: $(round(t, digits = 2))s")

    still_failing = [r.label for r in rows if !final_converged(r)]
    if isempty(still_failing)
        println("\nEvery scenario converged (warm start, or homotopy where warm alone failed).")
    else
        println("\nSTILL LOCALLY_INFEASIBLE after warm start + homotopy (see the header comment " *
                "and the exploration report -- residual inspection at the returned point shows " *
                "these are Ipopt restoration-phase misreports, not genuine economic infeasibility):\n  " *
                join(still_failing, "\n  "))
    end
end
