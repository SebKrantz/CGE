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

### Production

- Cobb-Douglas production function with three heterogeneous labour types and fixed capital
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

- Revenue: tariffs + export duties + indirect production taxes
- Expenditure: fixed-share consumption (`gles`) of fixed volume (`gdtot`) + budget surplus/deficit (`govsav`)

### Savings and Investment

- Household savings: fixed marginal propensity to save (`mps = 0.09305`)
- Total savings = household + government + depreciation + foreign savings
- Investment allocated across sectors by fixed shares `kio`; capital goods sourced via `imat`

### Closure rules

| Variable | Exogenous value |
|----------|----------------|
| Capital stocks `k` | Fixed at base year `k0` |
| World prices `pwm` | Fixed (small open economy) |
| Labour supplies `ls` | Fixed (inelastic labour supply) |
| Tariff rates `tm` | Fixed at base `tm0` (policy instrument) |
| Foreign savings `fsav` | Fixed at base `fsav0` |
| Household saving rate `mps` | Fixed at 0.09305 |
| Government consumption volume `gdtot` | Fixed at base `gdtot0` |

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
equilibrium — can return Ipopt's `LOCALLY_INFEASIBLE` even though the shocked equations
are, in fact, satisfiable (e.g. a +5pp import tariff on `services`). See
`test/robustness_grid.jl`'s header comment for the full exploration; in short, this is a
solver artifact of starting far from the new equilibrium with a pure feasibility objective
(`Min 1`), not a real economic infeasibility.

Fix it by warm-starting from a previously-solved `Simulation` — typically the baseline:

```julia
baseline, = solve(params)
shocked = with_shocks(params, Dict{Symbol,Any}(:tm0 => Dict(:services => params.tm0[:services] + 0.05)))
sim, status, converged, iterations, elapsed = solve(shocked; start = baseline)
```

This alone fixes most modest single-shock scenarios (in the grid: 27/59 cold → 44/59 warm)
and costs nothing extra (a warm solve typically takes <0.03s, vs. ~2s cold). For a
residual set of larger or sector-specific shocks that still fail even warm, retry with a
`steps`-increment homotopy from the pre-shock `Params` (`base`):

```julia
sim, status, converged, iterations, elapsed =
    solve(shocked; start = baseline, steps = 4, base = params)
```

This applies the shock in `steps` equal linear increments, warm-starting each from the
previous one's solution. (In the grid's residual ~15/59 failures, homotopy up to 100 steps
did not help either — residual inspection at the returned point shows every real equation
already satisfied to ~1e-13, so these are Ipopt restoration-phase misreports rather than
genuine infeasibilities, most likely tied to the model's pre-existing "structurally zero
but bounded away from zero at 1e-6" variables, e.g. `gd` for the 10 sectors with
`gles = 0`. Relaxing those bounds isn't an option here — it would break the bit-for-bit
match to `test/reference.json`.)

### API

| Function | Signature | Purpose |
|---|---|---|
| `load_data` | `(path) -> RawData` | Reads the 5 `camdata.xlsx` sheets (same hardcoded ranges as before) |
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
module-level scalar, is now `Params.er`. Both are shockable via `with_shocks` like any
other parameter.

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
bit-for-bit regression guard on the refactor). It also runs a reduced (~5-scenario)
robustness testset confirming `start`'s warm-start fix on the exact cases the original
bug report found.

`test/robustness_grid.jl` (`julia --project=. test/robustness_grid.jl`) is the full
~60-scenario grid behind that testset — see its header comment for the shock list and
the cold/warm/homotopy convergence counts.

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
