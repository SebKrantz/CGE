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

```julia
include("camcge.jl")
# Baseline and Simulation 1 run automatically and are stored as Simulation structs.
```

Access results:

```julia
baseline.xd   # domestic output by sector (baseline)
sim1.xd       # domestic output after +10% agricultural capital shock
sim1.y        # private GDP in Simulation 1
```

---

## Simulations

| # | Description | Change |
|---|-------------|--------|
| Baseline | Replicates 1979-80 Cameroon equilibrium | — |
| Simulation 1 | Public investment in subsistence agriculture | `k0[:agsubsist] × 1.10` |

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
