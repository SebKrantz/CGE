# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A static, single-country Computable General Equilibrium (CGE) model of Cameroon, calibrated to the
1979-80 national accounts and input-output table (Davis-De Melo-Robinson framework). Implemented in
Julia with JuMP (algebraic modelling) and Ipopt (interior-point NLP solver). The entire model —
sets, calibration, the model function, and two simulation runs — lives in the single file `cge.jl`.
See `README.md` for the full equation-by-equation description of the economics and for the current
function-level API.

**`n-sector` branch:** the model is no longer Cameroon-shaped. Sets and economy-wide scalars come
from the data (see Architecture stages 1–2 and DATA.md §5), a household direct tax `td0` was added,
and `hhsav`/`govsav`/`fsav`/`td`/`mps`/`tm`/`tariff`/`indtax`/`gr` are free variables so a
government deficit, a trade surplus, a dissaving household, a zero tariff line or net production
subsidies are representable. Variable lower bounds and Ipopt's bound/barrier constants are
**scale-relative**, not the absolute `1e-6`/`1e-5`/`0.01`/`0.1` that were implicitly calibrated to
Cameroon's billion-CFAF magnitudes (see stage 4 below). With `td0 = 0` the Cameroon results are
unchanged — `test/runtests.jl` still checks them bit-for-bit against `test/reference.json`.

**`eps-integration` branch (do not merge into `main`'s history without review):** `cge.jl` was
refactored from a top-to-bottom script with module-level mutable globals into a module
(`CGECameroon`) with **no include-time solves and no module-level mutable parameter state**, so it
can be embedded safely in a long-lived, multi-request host process (the Economic Policy Simulator
app includes this file inside its own wrapper module). Every equation is mathematically identical to
the pre-refactor script — see `test/runtests.jl`, which checks the refactor bit-for-bit against a
snapshot of the original script's output (`test/reference.json`).

## Running the model

`Project.toml`/`Manifest.toml` now pin `XLSX`, `JuMP`, `Ipopt` (added via `Pkg.add` inside this
repo's own environment — see `git log` on `eps-integration`). From the repo root:

```julia
import Pkg; Pkg.activate(@__DIR__)
include("cge.jl")
using .CGECameroon
baseline, sim1 = CGECameroon.example()   # loads data, calibrates, solves both runs, prints a summary
```

Or drive it step by step — see README.md's "API" table for `load_data`, `calibrate`, `solve`,
`with_shocks`, `example`.

**Compatibility note (verified against Julia 1.12 / current JuMP & Ipopt):** `solve()` builds the
model with `Model(Ipopt.Optimizer)`. An earlier version of this script used the pre-1.0 JuMP API
`Model(with_optimizer(Ipopt.Optimizer))`, which raises `UndefVarError: with_optimizer not defined` on
any current JuMP version — that's been fixed. The rest of the model (including `@NLconstraints`/
`@NLobjective`, still legacy-supported) runs unmodified and both `baseline` and `sim1` solve
successfully ("Solved To Acceptable Level", i.e. JuMP termination status `ALMOST_LOCALLY_SOLVED`,
which `solve()` treats as converged alongside `OPTIMAL`/`LOCALLY_SOLVED`).

The README's "Running the model" section used to refer to the file as `camcge.jl` — the actual
filename in this repo is, and always was, `cge.jl`.

`test/runtests.jl` is the test suite (`julia --project=. test/runtests.jl`); there is no linter or
build process.

## Architecture

`cge.jl` (module `CGECameroon`) is organised as four functions, in the same four stages the
pre-refactor script ran top-to-bottom (the order is preserved because later stages depend on earlier
ones — don't reorder without checking what each later block reads):

1. **Sets** — **data, not constants** (`n-sector`): the sector list, the traded/non-traded partition
   and the labour categories are read off the workbook into `RawData`, carried through `Params`, and
   used as the index sets of every `@variable`/`@NLconstraint`, so the model runs on any N sectors ×
   L labour categories. The module `const`s at the top (`SEC`, `IT`, `ITN`, `LC`, `WA0`, `SCALARS`)
   are now only the **defaults the legacy loader stamps on `data/camdata.xlsx`**, which has no sheets
   for them; never index off them anywhere else.

2. **`load_data(path) -> RawData`**: two layouts, chosen by whether the workbook has a `sectors`
   sheet. **Generic** (`_load_generic`): sheets `sectors`, `labour`, labelled `iotable`, `imat`,
   `employment`, `wagedist`, `miscellaneous` and `scalars`, all read by **header lookup** so any
   sector/labour order works — this is the format for a new country, fully specified in DATA.md §5
   (`test/generic_workbook.jl` is a reference writer for it). **Legacy** (`_load_legacy`): the
   original hardcoded `camdata.xlsx` cell ranges and `MISC_ROWS` row order, byte-for-byte unchanged.
   Either way the result is `Dict`s keyed by `(sector, sector)` / `(sector, labour)` tuples plus the
   sets, the per-labour base wages `wa0` and the economy-wide `scalars`. `_sheet` (generic path
   only) recomputes a worksheet's true extent from its cells when the file's own `<dimension>` tag
   under-reports it — openxlsx stamps a placeholder `<dimension ref="A1"/>` on every sheet, so every
   R-built country workbook otherwise read as one cell per sheet and failed on the first header
   lookup.

3. **`calibrate(raw::RawData) -> Params`**: closed-form algebra that backs out share/shift parameters
   (`delta`, `ac`, `gamma`, `at`, `ad`, `alphl`, ...) so the model reproduces the base-year SAM
   exactly, returned as a `Params` (a plain mutable struct holding every calibrated parameter and
   closure value — including `mps0`, the household saving rate promoted out of a bare literal, `er`,
   the exchange rate, `td0`, the household direct-tax rate, and `tariff0`, the base tariff revenue
   start value, derived from the data instead of the old hardcoded `76.548`). It raises a named error
   for any traded sector with `m0 <= 0` or `e0 <= 0`, which used to calibrate to NaN/Inf. Several
   quantities are computed twice deliberately (e.g. `x0`/`ac` before and after `ad` is derived) to
   maintain internal consistency, exactly as in the original.

4. **`solve(p::Params; silent=true, tol=1e-8, max_iter=3000, start=nothing, steps=1, base=nothing) ->
   (sim, status, converged, iterations, elapsed)`**: `solve` itself is now a thin wrapper — when
   `steps == 1` (default) it calls `_solve_once(p; ..., start)` directly; when `steps > 1` it walks a
   linear homotopy from `base` (the pre-shock `Params`, via `_interp_params`) to `p` in `steps` equal
   increments, warm-starting each from the previous one's solution, and returns the last step's
   result. `_solve_once` walks `MU_INIT_LADDER` — it calls `_solve_attempt` at the small barrier
   start that suits a good start point and, only if that fails to converge, retries once at a larger
   one for the inconsistent-data case (see the `MU_INIT_LADDER` docstring). `_solve_attempt` is the
   original model-building body: builds a fresh JuMP `Model` reading
   every parameter from `p` (never from a global), declares ~36 endogenous variable groups, ~30
   `@NLconstraint` equation blocks (price, production/factor, trade (Armington/CET), demand, and
   closure-rule blocks — see the README's "Key equations" table), sets a dummy objective (`Min 1`) so
   the solve is a pure feasibility problem, calls `set_silent`/sets Ipopt `tol`/`max_iter`, calls
   `JuMP.optimize!`, and extracts all `JuMP.value.(...)` results into a `Simulation`. Every
   `@variable`'s `start` clause now runs through a small `sv0`/`sv1`/`sv2` closure that reads from
   `start::Union{Nothing,Simulation}` when given, falling back to exactly the old base-year `*0`
   default (or `nothing`, for the handful of variables that never had one) when `start === nothing`.
   See the README's "Solving shocked scenarios reliably" section and `test/robustness_grid.jl`'s
   header comment for why `start`/`steps` exist (a cold solve of a shocked `Params` can report a
   spurious Ipopt `LOCALLY_INFEASIBLE`). `mu_strategy = "adaptive"` was tried and rejected as a global
   default — it actually breaks the *baseline* solve.

   **Nothing in the model may be an absolute number in the data's units.** `_solve_once` used to
   carry four, all silently calibrated to Cameroon (billion CFAF; sector outputs 30–400, smallest
   positive base-year cell 0.023): the `>= 1e-6` lower bound on every variable, Ipopt's
   `bound_relax_factor = 1e-5`, and Ipopt's own defaults `bound_push`/`bound_frac = 0.01` and
   `mu_init = 0.1`. A generic workbook in billion USD (DATA.md §5) spans 1e2 down to 1e-11 in the
   same columns, and there each of the four breaks the solve outright — the bound excludes the base
   year, the relaxation lets quantities go negative into `x^(-rhoc)` (NaN → restoration failure), the
   push throws away every start value below 0.01, and the barrier gives 1e-11-scale cells 1e10-scale
   duals. All four are now relative: a `qlb(v) = 1e-6*v` bound (a millionth of the variable's *own*
   base level, computed in the new "BASE-YEAR LEVELS" block, which also supplies the eleven start
   values that used to be missing), `bound_relax_factor`/`bound_push`/`bound_frac = 1e-12` (only
   large enough to keep the `fix`ed cells from being eliminated as parameters, which would leave one
   more equation than free variable — Walras' law makes one redundant) and `mu_init = 1e-9` (this is
   a pure feasibility problem started at a point that is already a solution; the barrier buys
   nothing). Cameroon is untouched by all of it; on the 188 30-sector country databases under
   `data/`, cold convergence went 88/188 → 188/188 and exact base-year replication 8 → 104.

   **Structurally-zero / indeterminate cells (fixed, not just worked around):** right after the
   `@variables` block, `_solve_once` `fix(...; force = true)`s ~20 cells that an equation forces to
   exactly 0 (or leaves fully indeterminate) *for any parameter value* — `gd[i]`/`cd[i]`/`dst[i]` for
   sectors with `gles[i]`/`cles[i]`/`dstr[i] == 0`, `id[i]` for sectors whose `imat` row is entirely
   zero, `labd[i,l]` for zero-base-employment cells, `duty` (`te ≡ 0`), and `e`/`m`/`pe`/`pm`/`pwe`/`tm`
   for the two non-traded sectors — and narrows the corresponding equation's index set to match (see
   the block's own comment for the equation that forces each one). These previously carried the same
   `>= 1e-6` lower bound as every other variable, an inherent bound-vs-equation conflict (not a bad
   start point) that was the actual root cause of most of `test/robustness_grid.jl`'s residual
   `LOCALLY_INFEASIBLE` cases. Fixing it took the grid from 44/59 to 59/59 converged warm-started
   (and, since the scale-relative rework above, 59/59 cold as well). This makes `start = nothing`
   (cold, the default)
   **no longer** byte-identical to the pre-fix script on these specific cells — expected, since they
   used to be wrong (lower-bound-clipped to ~9.9e-7, or, where nothing pinned them, an Ipopt runaway
   to ~1e5); `test/runtests.jl`'s reference.json check now excludes exactly these cells (plus
   `tm[services]` and the `tariff`/`govsav` aggregates that ripple from them by a comparable tiny
   amount) via `skip_cell`, and still requires 1e-6 relative agreement on every other cell.

`with_shocks(p::Params, overrides::Dict{Symbol,Any}) -> Params` replaces the old "mutate a global,
call `cammodel()` again" pattern: it returns a `deepcopy` of `p` with `overrides` applied (a scalar
replaces a scalar field; a `Dict{Symbol,<:Real}` merges new levels into a per-sector/per-labour
field), so `p` itself is never mutated and two scenarios (even solved back to back, or from two
concurrent callers) cannot leak into each other. `example()` reproduces the original script's
include-time `baseline`/`sim1` runs on demand.

## Data

All calibration data live in `data/camdata.xlsx`. See the README's "Data" table for sheet names,
contents, and dimensions (11×11 IO/capital matrices, 11×3 wage/employment matrices, 17×11 scalar
parameter sheet). Its cell ranges are hardcoded in `_load_legacy` (stage 2 above); if that
spreadsheet's layout changes, the `XLSX.readdata` ranges and the `MISC_ROWS` ordering must be updated
to match. **A new dataset should not use that layout** — write the generic workbook of DATA.md §5
instead, which is read by header lookup and carries its own sets and scalars.
