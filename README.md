# Cameroon CGE Model

A static, single-country Computable General Equilibrium model calibrated to Cameroon's **1979-80** national accounts and input-output table.

---

## Background

This model follows the **Davis-De Melo-Robinson** framework described in Dervis, De Melo and Robinson (1982) and was originally implemented in GAMS by Condon, Dahl and Devarajan (1987) for the World Bank.  This version is a Julia reimplementation using **JuMP** for model specification and **Ipopt** for solution.

---

## Model structure

### Sets

| Set | Name | Members |
|-----|------|---------|
| `SEC` | All sectors | 11 |
| `IT` | Traded sectors | 9 (agsubsist, agexpind, sylvicult, indalim, bienscons, biensint, cimint, bienscap, services) |
| `ITN` | Non-traded sectors | 2 (construct, publiques) |
| `lc` | Labour categories | rural, urbanunsk, urbanskil |

The sets are **data, not constants**: sector codes, the traded/non-traded partition and the
labour categories are read from the workbook into `RawData` and carried through `Params` into
every `@variable`/`@NLconstraint` index set, so the model runs on any N sectors × L labour
categories. The members above are Cameroon's, supplied as defaults for the legacy
`data/camdata.xlsx` (which carries no sets of its own) — see DATA.md §5 for the generic
workbook format that states them as data.

### Production

- Cobb-Douglas production function with heterogeneous labour types and fixed capital
- Sector-specific total factor productivity (TFP) parameter `ad` calibrated to base year
- Leontief intermediate input demands (fixed IO coefficients `io[i,j]`)
- Depreciation reduces capital services available to each sector

### Trade (Armington and CET)

| Direction | Mechanism | Elasticity |
|-----------|-----------|------------|
| Imports | Armington CES: composite supply = domestic + imports | `1/(1+rhoc)` |
| Exports | CET: output allocated between domestic sales and exports | `1/(rhot-1)` |
| Foreign demand | Downward-sloping export demand curve | `eta` |

Non-traded sectors (`construct`, `publiques`) have neither imports nor exports.

### Government

- Revenue: tariffs + export duties + indirect production taxes + a direct tax on household income (`td·y`)
- Expenditure: fixed-share consumption (`gles`) of fixed volume (`gdtot`) + budget surplus/deficit (`govsav`, free to be negative)

### Savings and Investment

- Household savings: fixed marginal propensity to save out of **disposable** income, `hhsav = mps·(1−td)·y` (Cameroon: `mps0 = 0.09305`, `td0 = 0`)
- Total savings = household + government + depreciation + foreign savings
- Investment allocated across sectors by fixed shares `kio`; capital goods sourced via `imat`

### Closure rules

| Variable | Exogenous value |
|----------|----------------|
| Capital stocks `k` | Fixed at base year `k0` |
| World prices `pwm` | Fixed (small open economy) |
| Labour supplies `ls` | Fixed (inelastic labour supply) |
| Tariff rates `tm` | Fixed at base `tm0` (policy instrument) |
| Foreign savings `fsav` | Fixed at base `fsav0` (negative = trade surplus) |
| Household saving rate `mps` | Fixed at base `mps0` |
| Household direct-tax rate `td` | Fixed at base `td0` (`0` = no direct tax) |
| Government consumption volume `gdtot` | Fixed at base `gdtot0` |

`hhsav`, `govsav`, `fsav` and `td` are declared **free** (no lower bound), so a base year with a
government deficit, a trade surplus or household dissaving is representable; every price and
quantity variable keeps its `>= 1e-6` bound.

---

## Data

All calibration data are in [`data/camdata.xlsx`](data/camdata.xlsx):

| Sheet | Contents | Dimensions |
|-------|----------|------------|
| `iotable` | Leontief input-output technical coefficients | 11 × 11 |
| `imat` | Capital composition matrix (investment demand by sector of origin) | 11 × 11 |
| `wagedist` | Wage proportionality factors by sector and labour type | 11 × 3 |
| `employment` | Base-year employment in 1000 persons | 11 × 3 |
| `miscellaneous` | Scalar parameters: volumes, prices, elasticities, shares | 17 × 11 |

`load_data` reads this **legacy** layout at fixed cell ranges. A workbook with a `sectors` sheet is
instead read as a **generic** workbook — `sectors`, `labour`, labelled `iotable`, `imat`,
`employment`, `wagedist`, `miscellaneous`, and a `scalars` sheet (`er`, `gr0`, `gdtot0`, `cdtot0`,
`fsav0`, `mps0`, `td0`, `tariff0`, `wa0_<labour code>`) — with any sector count, any
traded/non-traded partition, any labour categories, in any order, all by header lookup. **DATA.md
§5 is the full specification**; `test/generic_workbook.jl` is a reference writer for it.

---

## Running the model

As of the `eps-integration` branch, `cge.jl` is a module (`CGECameroon`) with **no
include-time solves and no module-level mutable parameter state** — every calibrated
parameter lives in a `Params` value you pass around explicitly, so two scenarios never
leak into each other. The file is still named `cge.jl` (not `camcge.jl`).

```julia
import Pkg; Pkg.activate(@__DIR__)   # uses this repo's own Project.toml
include("cge.jl")
using .CGECameroon

raw    = load_data(joinpath(@__DIR__, "data", "camdata.xlsx"))  # read the workbook
params = calibrate(raw)                                          # closed-form calibration
baseline, status, converged, iterations, elapsed = solve(params) # fresh JuMP model, Ipopt
```

Or, to reproduce the two runs below in one call:

```julia
baseline, sim1 = CGECameroon.example()
```

Access results (unchanged — `Simulation` still has the same 38 fields):

```julia
baseline.xd   # domestic output by sector (baseline)
sim1.xd       # domestic output after +10% agricultural capital shock
sim1.y        # private GDP in Simulation 1
```

Run a shock with `with_shocks`, which returns a shocked **deep copy** of `params` — the
original is never mutated:

```julia
shocked = with_shocks(params, Dict{Symbol,Any}(:k0 => Dict(:agsubsist => params.k0[:agsubsist] * 1.10)))
sim1, status, converged, iterations, elapsed = solve(shocked)
```

A shock value is a dict key naming a `Params` field and either a new scalar level
(e.g. `:gdtot0 => 145.0`, `:fsav0 => 55.0`, `:er => 0.25`) or a
`Dict{Symbol,<:Real}` of new per-sector/per-labour levels (e.g.
`:tm0 => Dict(:services => 0.30)`). Values are always the new **level**, never a
percent change.

### Solving shocked scenarios reliably: `start` and `steps`

Solving a shocked `Params` cold — every variable starting from the base-year `*0` data,
exactly like the original script, regardless of how far the shock has moved the
equilibrium — can still return Ipopt's `LOCALLY_INFEASIBLE` even though the shocked
equations are, in fact, satisfiable (e.g. a +5pp import tariff on `services`). See
`test/robustness_grid.jl`'s header comment for the full exploration; in short, this is a
solver artifact of starting far from the new equilibrium with a pure feasibility objective
(`Min 1`), not a real economic infeasibility.

Fix it by warm-starting from a previously-solved `Simulation` — typically the baseline:

```julia
baseline, = solve(params)
shocked = with_shocks(params, Dict{Symbol,Any}(:tm0 => Dict(:services => params.tm0[:services] + 0.05)))
sim, status, converged, iterations, elapsed = solve(shocked; start = baseline)
```

Warm-starting fixes every scenario in `test/robustness_grid.jl`'s ~60-scenario grid
(59/59, up from 44/59 before the structural-zero fix below) and costs nothing extra (a
warm solve takes well under 0.1s, vs. ~2s cold). `steps`-increment homotopy from the
pre-shock `Params` (`base`) remains available for anything that somehow still fails warm:

```julia
sim, status, converged, iterations, elapsed =
    solve(shocked; start = baseline, steps = 4, base = params)
```

**Root cause fixed, not just worked around.** The ~15/59 scenarios that used to resist
even warm-starting were traced to a real bug, not just solver fragility: a number of
cells (`gd`, `cd`, `dst`, `id`, `labd` for particular sector/labour combinations, `duty`,
and `e`/`m`/`pe`/`pm`/`pwe`/`tm` for the two non-traded sectors) are forced to exactly 0
(or otherwise indeterminate) by an equation for *any* parameter value, yet were declared
with the same `>= 1e-6` lower bound as every other variable — a hard bound-vs-equation
conflict Ipopt could only paper over near the base-year start point (see
`test/reference.json`, where these cells sit at ~9.9e-7, or, where nothing pins them at
all, run away to ~1e5). `_solve_once` now fixes each such cell directly with
`fix(...; force = true)` and narrows the equation that used to (mis)determine it — see
cge.jl's "STRUCTURALLY ZERO / INDETERMINATE CELLS" block for the full list and the
equation that forces each one.

### API

| Function | Signature | Purpose |
|---|---|---|
| `load_data` | `(path) -> RawData` | Reads a base-year workbook: the **generic** layout (8 sheets, header lookup, sets and scalars as data — DATA.md §5) when it has a `sectors` sheet, the **legacy** `camdata.xlsx` layout (5 sheets, hardcoded ranges, Cameroon sets as defaults) otherwise |
| `calibrate` | `(raw::RawData) -> Params` | Closed-form calibration algebra (unchanged economics) |
| `solve` | `(p::Params; silent=true, tol=1e-8, max_iter=3000, start=nothing, steps=1, base=nothing) -> (sim, status, converged, iterations, elapsed)` | Builds a fresh JuMP model from `p`, solves with Ipopt; `start` warm-starts from a previous `Simulation`, `steps`>1 homotopies from `base` (see above) |
| `with_shocks` | `(p::Params, overrides::Dict{Symbol,Any}) -> Params` | Deep-copies `p` and applies level overrides |
| `example` | `() -> (baseline=Simulation, sim1=Simulation)` | Reproduces the two runs below and prints a summary |

`solve`'s `converged` is `true` for JuMP termination statuses `OPTIMAL`,
`LOCALLY_SOLVED`, or `ALMOST_LOCALLY_SOLVED` (Ipopt reports "Solved To Acceptable
Level" — `ALMOST_LOCALLY_SOLVED` — for both runs below at default tolerances); `status`
is the raw termination status as a string.

The household saving rate, previously a bare literal (`0.09305`) inside the
`closuremp` constraint, is now `Params.mps0`; the exchange rate, previously a
module-level scalar, is now `Params.er`; the direct-tax rate is `Params.td0` and the
base tariff revenue `Params.tariff0` (a start value, derived from the data instead of
the old hardcoded `76.548`). All are shockable via `with_shocks` like any other
parameter.

---

## Simulations

| # | Description | Change |
|---|-------------|--------|
| Baseline | Replicates 1979-80 Cameroon equilibrium | — |
| Simulation 1 | Public investment in subsistence agriculture | `k0[:agsubsist] × 1.10` |

### Tests

`test/runtests.jl` (`julia --project=. test/runtests.jl`) checks that the baseline
reproduces the base-year data and that `sim1` matches `test/reference.json` — a
snapshot of the original, pre-refactor script's output — to 1e-6 relative on every
value (both solved cold, exactly as the original script did, so this remains a
bit-for-bit regression guard on the refactor), **except** the structurally-zero/
indeterminate cells listed above plus two aggregates that legitimately shift by a
comparable tiny amount because they sum those cells (`tariff`, `govsav`) — see
`test/runtests.jl`'s `skip_cell` for the exact list and the reasoning for each. It also
runs a reduced (~5-scenario) robustness testset confirming `start`'s warm-start fix on
the exact cases the original bug report found.

Two further testsets cover the N-sector work: a **generic-workbook round trip** (the Cameroon
`RawData` written out in the generic layout with its sectors and labour categories reversed, re-read
through `load_data`'s generic branch, and required to reproduce the legacy calibration and solution
to 1e-8) and a **synthetic 3-sector / 2-labour dataset** (`test/synthetic_data.jl`) with a base-year
trade surplus, a government deficit and a 5% direct tax, which calibrates and solves back to its own
base year.

`test/robustness_grid.jl` (`julia --project=. test/robustness_grid.jl`) is the full
~60-scenario grid behind that testset — see its header comment for the shock list and
the cold/warm/homotopy convergence counts (59/59 converge warm-started).

---

## Key equations

| Equation | Description |
|----------|-------------|
| `pmdef` | Import price = world price × ER × (1 + tariff) |
| `absorption1/2` | Armington composite price identity |
| `actp` | Activity price = value added + intermediate costs − indirect tax |
| `activity` | Cobb-Douglas production function |
| `profitmax` | Labour demand from profit maximisation FOC |
| `cet` | CET output transformation (domestic vs. export) |
| `armington` | CES import-domestic aggregation |
| `equil` | Market clearing: supply = intermediate + consumption + government + investment + inventory |
| `totsav` | Savings-investment balance |
| `caeq` | Current account balance |

---

## References

- Dervis, K., De Melo, J., and Robinson, S. (1982). *General Equilibrium Models for Development Policy*. Cambridge University Press.
- Condon, T., Dahl, H., and Devarajan, S. (1987). "Implementing a computable general equilibrium model on GAMS: The Cameroon model." World Bank Discussion Paper DRD290.
