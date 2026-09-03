# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A static, single-country Computable General Equilibrium (CGE) model of Cameroon, calibrated to the
1979-80 national accounts and input-output table (Davis-De Melo-Robinson framework). Implemented in
Julia with JuMP (algebraic modelling) and Ipopt (interior-point NLP solver). The entire model —
sets, calibration, the model function, and two simulation runs — lives in the single file `cge.jl`.
See `README.md` for the full equation-by-equation description of the economics.

## Running the model

There is no `Project.toml`/`Manifest.toml` in the repo, so dependencies are not pinned. Install once
per environment:

```julia
import Pkg; Pkg.add(["XLSX", "JuMP", "Ipopt"])
```

Then, from the repo root (data is loaded via a path relative to the script's own location, so the
working directory doesn't matter, but the include path does):

```julia
include("cge.jl")
```

**Compatibility note (verified against Julia 1.12 / current JuMP & Ipopt):** `cammodel()` builds the
model with `Model(Ipopt.Optimizer)`. An earlier version of this script used the pre-1.0 JuMP API
`Model(with_optimizer(Ipopt.Optimizer))`, which raises `UndefVarError: with_optimizer not defined` on
any current JuMP version — that's been fixed. The rest of the script (including `@NLconstraints`/
`@NLobjective`, still legacy-supported) runs unmodified and both `baseline` and `sim1` solve
successfully ("Solved To Acceptable Level").

The README's "Running the model" section refers to the file as `camcge.jl` — the actual filename in
this repo is `cge.jl`.

No test suite, linter, or build process exists in this repo.

## Architecture

`cge.jl` executes top-to-bottom as a script, in four stages that must stay in this order:

1. **Sets and empty parameter containers** (~L27-114): `SEC` (11 sectors), `IT`/`ITN` (traded /
   non-traded split), `lc` (3 labour categories). All CGE parameters (`delta`, `ac`, `rhoc`, `rhot`,
   `gamma`, `ad`, `alphl`, `cles`, etc.) are declared as empty `Dict()`s at module scope here and
   populated later — they are genuine globals, not passed as function arguments.

2. **Data loading** (~L139-180): reads fixed cell ranges from `data/camdata.xlsx` (sheets `iotable`,
   `imat`, `wagedist`, `employment`, `miscellaneous`) into `Dict`s keyed by `(sector, sector)` or
   `(sector, labour)` tuples. Row order in the `miscellaneous` sheet (`rowz`) must match the sheet's
   actual layout — there's no header-based lookup.

3. **Calibration** (~L188-273): closed-form algebra that backs out share/shift parameters (`delta`,
   `ac`, `gamma`, `at`, `ad`, `alphl`, ...) so the model reproduces the base-year SAM exactly. Several
   quantities are computed twice deliberately (e.g. `x0`/`ac` before and after `ad` is derived) to
   maintain internal consistency — this is order-dependent; don't reorder blocks without checking
   what each later block depends on.

4. **`cammodel()`** (~L283-507): builds a fresh JuMP `Model`, declares ~36 endogenous variable
   groups, ~30 `@NLconstraint` equation blocks (price, production/factor, trade (Armington/CET),
   demand, and closure-rule blocks — see the README's "Key equations" table), sets a dummy objective
   (`Min 1`) so the solve is a pure feasibility problem, calls `JuMP.optimize!`, extracts all
   `JuMP.value.(...)` results, and returns them as a large tuple. **It reads the calibrated globals
   from stage 3 as constants** (`ad`, `alphl`, `delta`, `ac`, `rhoc`, `rhot`, `gamma`, `at`, `io`,
   `imat`, `kio`, `cles`, `gles`, `dstr`, `itax`, `tm0`, `k0`, `ls0`, `fsav0`, `gdtot0`, ...) rather
   than taking them as arguments.

Results from each call to `cammodel()` are wrapped in a `Simulation` struct (all fields `::Any`,
since `JuMP.value.()` returns `DenseAxisArray`s). `baseline` is solved directly from calibrated
parameters. `sim1` is produced by **mutating a global parameter** (`k0[:agsubsist] = 1.10 * k0[:agsubsist]`)
and calling `cammodel()` again — this is the pattern for adding any new counterfactual: mutate the
relevant global parameter/closure value, call `cammodel()`, wrap the result in a new `Simulation`.
Because closure rules (`closurek`, `closurep`, `closurel`, `closuret`, `closuref`, `closuremp`,
`closureg`, ...) fix most parameters to their `*0` base-year values by default, a new scenario
typically means changing one of those `*0` globals (or editing a `closure*` constraint directly)
before re-solving.

## Data

All calibration data live in `data/camdata.xlsx`. See the README's "Data" table for sheet names,
contents, and dimensions (11×11 IO/capital matrices, 11×3 wage/employment matrices, 17×11 scalar
parameter sheet). Cell ranges are hardcoded in the loading code (stage 2 above); if the spreadsheet
layout changes, the `XLSX.readdata` ranges and `rowz` ordering in `cge.jl` must be updated to match.
