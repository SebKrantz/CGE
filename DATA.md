# DATA.md — Input data contract for the Cameroon CGE model

This document specifies what a base-year dataset for this model must contain, in what exact
shape, so that someone can build a new country's input workbook **without reading `cge.jl`**.
Section 3 documents the current workbook (`data/camdata.xlsx`) cell-by-cell, as actually read by
`load_data`/`calibrate`. Section 5 proposes a generic long-format schema for future countries (not
implemented). Section 6 only suggests where raw numbers might come from — not part of the data
contract, since `load_data(path)` takes a path argument and knows nothing about the file's origin.

All `cge.jl:N` references are to this repo's current `eps-integration` branch file (930 lines).
All sheet contents quoted below were read directly from `data/camdata.xlsx` with `XLSX.jl`.

---

## 1. What the model needs, conceptually

The model is a **base-year social accounting picture for one country**: a single period, single
representative household, single "rest of world" account, **11 production sectors** (2 of which
are non-traded) and **3 labour categories**. Concretely, a new country's data must supply:

- An **input-output table**: intermediate use of each sector's output by every other sector
  (sector × sector, "who buys from whom" as a share of the buying sector's gross output).
- **Imports and exports by sector** (for the traded sectors only — the two non-traded sectors have
  neither).
- **Final demands**: household consumption by sector, government consumption by sector,
  investment by sector of *destination* (via fixed shares) translated into investment by sector of
  *origin* via a **capital composition matrix** (`imat`) — i.e. "a unit of investment destined for
  sector X is made up of these shares of goods from sectors A, B, C…".
- **Value added**, split into **labour income by 3 labour categories** and a **capital** residual,
  plus a **depreciation rate** on capital.
- **Taxes**: import tariffs, export duties, and indirect/production taxes — all as sector-specific
  ad-valorem rates.
- The **external account**: foreign savings (the current-account deficit, in foreign currency).
- **Employment by sector × labour category** (levels) and a **wage-distribution matrix** (relative
  wages of each labour category within each sector).
- A handful of **economy-wide scalars**: the exchange rate, total government spending, total
  household consumption, the household saving rate, government revenue (a start value only).

**Units.** Almost everything is in **billion base-year CFAF** (flows and capital stocks alike).
Employment is in **1000 persons**. Wages are in **million CFAF per 1000 workers**. The exchange
rate and foreign savings are the only quantities in **USD** (billion), because the current-account
identity is expressed at world prices (`caeq`, cge.jl:792). All base-year prices are normalised to
**1** (confirmed: every `pd0` cell in the workbook is exactly `1.0`).

**Base year convention.** The model calibrates to exactly one year (here, Cameroon's 1979–80).
There is no time dimension in the data or the equations — "base year" means the single
cross-section the SAM-like data describes, and every calibrated parameter is derived so the model
reproduces that cross-section exactly (`calibrate`, cge.jl:229–347; see §4).

**Scope note.** This is a **single-country, single-household** model: one representative
household (no income deciles, no rural/urban split), one government, and the rest of the world
appears only through the trade/current-account equations — not as a modelled region. The set
`SEC` fixes **11 sectors** and `LC` fixes **3 labour categories**; of the 11 sectors, **2 are
non-traded** (`construct`, `publiques`) and have no import/export data at all (cge.jl:57–75).

---

## 2. Sets: sectors and labour categories

| # | Code | Label (from `cge.jl` comments) | Traded? |
|---|---|---|---|
| 1 | `agsubsist` | Food crops (subsistence agriculture) | traded |
| 2 | `agexpind` | Cash crops (export-oriented agriculture) | traded |
| 3 | `sylvicult` | Forestry and logging | traded |
| 4 | `indalim` | Food processing industry | traded |
| 5 | `bienscons` | Consumer goods manufacturing | traded |
| 6 | `biensint` | Intermediate goods manufacturing | traded |
| 7 | `cimint` | Construction materials | traded |
| 8 | `bienscap` | Capital goods manufacturing | traded |
| 9 | `construct` | Construction | **non-traded** |
| 10 | `services` | Private services | traded |
| 11 | `publiques` | Public services | **non-traded** |

Source: `SEC` cge.jl:57–69, `IT` (the 9 traded sectors) cge.jl:71–72, `ITN` (the 2 non-traded
sectors) cge.jl:73. **Sector position and count are structural** — they are Julia `const`s, not
data read from any sheet. The 11-way split and the traded/non-traded partition (positions 9 and 11
are non-traded) cannot change without editing `cge.jl` itself.

| Labour category (`LC`) | Meaning |
|---|---|
| `rural` | Rural labour |
| `urbanunsk` | Urban unskilled labour |
| `urbanskil` | Urban skilled labour |

Source: cge.jl:75.

---

## 3. The exact current workbook contract (`data/camdata.xlsx`)

Five sheets, read by `load_data(path)` (cge.jl:108–142) with **hardcoded cell ranges and no
header-based lookup** — column/row headers in the workbook are free text for a human reader only;
the code never looks at them. Data is read positionally: **column order and row order must match
the 11-sector / 3-labour ordering in §2 exactly.**

### 3.1 `iotable` — input-output technical coefficients

Read range: `B3:L13` → **11×11** (cge.jl:109). Actual sheet dimension: `A1:L13` (no extra rows).

| Cell | Meaning |
|---|---|
| Row 1 | Free-text title (`"Table io(i,j) Input-Output Coefficients (Unity)'"`), not read |
| Row 2 | Sector-name header, not read (some cells have trailing spaces, e.g. `"agexpind "`) |
| Rows 3–13, Col A | Sector row label, not read (row *position* is what's used) |
| Rows 3–13 × Cols B–L | `io[i,j]` = value of good *i* (row) used as an intermediate input per unit **gross output value** of sector *j* (column). Unitless share. |

Row = supplying sector *i*, column = using sector *j*, in §2 order; all cells non-negative shares.
Sample cells (exact values, for orientation): `io[agsubsist,indalim]=0.30266` (subsistence-agriculture
input into food processing), `io[cimint,cimint]=0.27608` (construction materials' own use of
itself), `io[services,agexpind]=0.30649` (services input into cash-crop agriculture). `io[i,j]`
columns are generally sparse — e.g. `sylvicult` (forestry) supplies non-zero inputs to only 2 of
the 11 sectors (`indalim` at `0.00243`, `biensint` at `0.02106`).

### 3.2 `imat` — capital composition matrix

Read range: `B3:L13` → **11×11** (cge.jl:115). **Actual sheet dimension is `A1:L19`** — 6 extra
blank rows (14–19) exist below the data but are not read; this is a real quirk of the file, not a
guess.

| Cell | Meaning |
|---|---|
| Rows 3–13 × Cols B–L | `imat[i,j]` = share of good *i* (row) in **one unit of investment good** destined for sector *j* (column) — the "recipe" that turns investment-by-destination into investment-by-origin (`pkdef`, cge.jl:681; `ieq`, cge.jl:789). |

In this workbook only **3 of 11 rows are non-zero**: `agsubsist` (only supplies itself, share
`0.23637` on its own diagonal, 0 elsewhere), `bienscap` (shares `0.5953`–`0.78723` across all 11
destination columns), and `construct` (shares `0.16833`–`0.8239`). The other 8 rows
(`agexpind, sylvicult, indalim, bienscons, biensint, cimint, services, publiques`) are **entirely
zero** — i.e. this economy's capital goods are supplied only by subsistence agriculture, capital
goods manufacturing, and construction. **Each column sums to 1.0** (e.g. the `agsubsist` column:
`0.23637 + 0.5953 + 0.16833 = 1.00000`; the `agexpind` column: `0 + 0.60608 + 0.39392 = 1.0`) —
this is the "capital composition columns sum to 1" consistency rule referenced in §4.

### 3.3–3.4 `employment` and `wagedist` — paired sector × labour-category sheets

`employment`: read range `B3:D13` → **11×3** (cge.jl:127), sheet dimension `A1:D13` (exact match).
Units **1000 persons**; cell meaning `xle[i,l]` = employment of labour category *l* in sector *i*.
`wagedist`: read range `B3:D13` → **11×3** (cge.jl:121), sheet dimension `A1:L13` (12 columns wide
in the file, but only B:D are populated — E:L are blank/`missing` and not read). Unitless
multiplier on the economy-wide wage `wa[l]` (a per-sector "wage premium/discount"); cell meaning
`wdist[i,l]`. Shown together because the two sheets must be **zero-paired** (see below):

| Sector | employment: rural | urbanunsk | urbanskil | wagedist: rural | urbanunsk | urbanskil |
|---|---:|---:|---:|---:|---:|---:|
| agsubsist | 1654.430 | 162.890 | **0** | 1.0189 | 0.71491 | **0** |
| agexpind | 399.930 | 45.508 | 5.057 | 0.49556 | 0.34774 | 0.29222 |
| sylvicult | 7.662 | 1.789 | 0.597 | 3.2628 | 2.2890 | 1.9232 |
| indalim | 12.989 | 9.434 | 2.358 | 1.4571 | 1.0223 | 0.85902 |
| bienscons | 28.344 | 37.462 | 12.488 | 1.1335 | 0.79531 | 0.66829 |
| biensint | 18.331 | 16.553 | 8.300 | 3.1074 | 2.1806 | 1.8323 |
| cimint | 1.458 | 1.317 | 0.660 | 6.3224 | 4.4364 | 3.7277 |
| bienscap | 3.112 | 2.820 | 1.208 | 2.5035 | 1.7552 | 1.4758 |
| construct | 22.584 | 28.462 | 7.116 | 2.9204 | 2.0492 | 1.7220 |
| services | 121.200 | 125.800 | 61.960 | 1.4039 | 0.98502 | 0.82776 |
| publiques | **0** | 83.029 | 32.771 | **0** | 1.3263 | 1.1146 |

The zero cells (`agsubsist×urbanskil`, `publiques×rural`) **coincide exactly** between the two
sheets — a required data-consistency pairing (see §4): wherever employment for a (sector, labour)
cell is zero, the calibrated labour share `alphl[l,i]` is zero regardless of `wdist`, and the model
fixes that `labd` cell to exactly zero rather than solving for it (cge.jl:592–600, 636, 642–644).

### 3.5 `miscellaneous` — 17 scalar-parameter rows × 11 sectors

Read range: `B3:L19` → **17×11** (cge.jl:133). Actual sheet dimension: `A1:L19` (exact match). Row
order is significant and is **not** looked up by row label — it is the hardcoded `rowz` array
(cge.jl:134–135):

`m0, e0, xd0, k, depr, rhoc, rhot, eta, pd0, tm0, itax, cles, gles, kio, dstr, dst, id`

Full sheet (rounded for display; all values non-negative):

| Row | agsubsist | agexpind | sylvicult | indalim | bienscons | biensint | cimint | bienscap | construct | services | publiques |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `m0` (imports, CFAF bn) | 2.461 | 8.039 | 0.023 | 17.961 | 37.062 | 138.57 | 49.616 | 134.72 | **0** | 74.439 | **0** |
| `e0` (exports, CFAF bn) | 4.594 | 125.07 | 22.337 | 23.451 | 5.864 | 101.33 | 10.501 | 3.838 | **0** | 81.626 | **0** |
| `xd0` (gross output, CFAF bn) | 330.48 | 131.45 | 29.503 | 72.024 | 118.43 | 284.38 | 34.169 | 10.298 | 174.12 | 615.79 | 163.98 |
| `k` (capital stock, CFAF bn) | 495.73 | 170.89 | 73.76 | 140.0 | 236.87 | 853.13 | 102.51 | 20.6 | 435.29 | 769.73 | 180.36 |
| `depr` (depreciation rate) | 0.0246 | 0.0472 | 0.0244 | 0.0144 | 0.0212 | 0.0335 | 0.0335 | 0.0111 | 0.0232 | 0.0637 | 0.0637 |
| `rhoc` (Armington σ) | 1.5 | 0.9 | 0.4 | 1.25 | 1.25 | 0.5 | 0.75 | 0.4 | 0.4 | 0.4 | 0.4 |
| `rhot` (CET σ) | 1.5 | 0.9 | 0.4 | 1.25 | 1.25 | 0.5 | 0.75 | 0.4 | 0.4 | 0.4 | 0.4 |
| `eta` (export-demand elasticity) | 1.0 | 1.0 | 1.0 | 4.0 | 4.0 | 4.0 | 4.0 | 4.0 | 4.0 | 4.0 | 4.0 |
| `pd0` (base domestic price) | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| `tm0` (tariff rate) | 0.2205 | 0.233 | 0.278 | 0.3534 | 0.3826 | 0.1768 | 0.2633 | 0.268 | **0** | **0** | **0** |
| `itax` (indirect tax rate) | 0.002 | 0.191 | 0.057 | 0.038 | 0.096 | 0.026 | 0.014 | 0.029 | 0.034 | 0.076 | **0** |
| `cles` (consumption share) | 0.2744 | 0.00445 | **0** | 0.05599 | 0.14099 | 0.17738 | **0** | **0** | 0.004 | 0.31921 | 0.02358 |
| `gles` (govt. spending share) | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **1** |
| `kio` (investment-destination share) | 0.11 | 0.09 | 0.06 | 0.01 | 0.04 | 0.14 | 0.02 | 0.01 | 0.08 | 0.34 | 0.1 |
| `dstr` (inventory/output ratio) | 0.012203 | 0.026694 | 0.034742 | 0.044291 | 0.059958 | 0.012287 | **0** | 0.042047 | **0** | **0** | **0** |
| `dst` (inventory investment, CFAF bn) | 4.033 | 3.509 | 1.025 | 3.19 | 7.101 | 3.494 | 0 | 0.433 | 0 | 0 | 0 |
| `id` (investment by origin, CFAF bn) | 6.71 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Notes tying this table back to §1's units: `m0, e0, xd0, k, dst, id` are CFAF bn; `depr, tm0, itax,
cles, gles, kio, dstr` are unitless rates/shares; `rhoc, rhot` are elasticities of substitution
*as reported* (σ, **not** the CES/CET curvature exponent ρ the equations actually use — `calibrate`
converts `rhoc_exponent = 1/σ−1` and `rhot_exponent = 1/σ+1`, cge.jl:243–244); `eta` is a demand
elasticity; `pd0` is a normalised price (always 1 at base). `dst` and `id` are **not** fed back
into calibration algebra — they are used only as JuMP solver start values (cge.jl:540–541), so a
new country's data for these two rows only needs to be a reasonable base-year cross-check, not an
input the calibration depends on.

Several economy-wide scalars are **not** on any sheet — hardcoded in `calibrate` (cge.jl:234–240)
and needed for a new country outside the workbook: base wages by labour category `wa0` (million
CFAF/1000 workers: rural `0.11`, urban-unskilled `0.15678`, urban-skilled `1.8657`) and the scalars
`er=0.21` (CFAF/USD), `gr0=179.0`, `gdtot0=135.03`, `cdtot0=947.98`, `fsav0=36.841` (billion USD),
`mps0=0.09305`. `te` (export duty rate) is **not read from any sheet** — hardcoded to `0.0` for
every sector (cge.jl:247): today's contract has no export-duty row.

---

## 4. Derived parameters — what `calibrate` computes, and the consistency rules it assumes

`calibrate(raw::RawData) -> Params` (cge.jl:229–347) is closed-form, non-econometric algebra: it
backs out share/shift parameters so every model equation holds **exactly** at the base-year data
point. In order of computation:

| Parameter | Formula (in words) | Reference |
|---|---|---|
| `pm0[i]`, `pe0[i]` | = `pd0[i]` — at base year, domestic, import, and export prices are assumed equal (all normalised to 1) | cge.jl:265–266 |
| `pwm0[i]` | implied world import price = `pm0[i] / ((1+tm0[i])·er)` | cge.jl:268 |
| `pwe0[i]` | implied world export price = `pe0[i] / ((1+te[i])·er)` | cge.jl:269 |
| `pva0[i]` | value-added price = `pd0[i] − Σ_j io[j,i]·pd0[j] − itax[i]` (output price net of intermediate cost and indirect tax); **must be positive**, i.e. `Σ_j io[j,i] + itax[i] < 1` for every sector | cge.jl:271–272 |
| `xxd0[i]` | domestic sales = `xd0[i] − e0[i]` | cge.jl:274 |
| `ls0[l]` | economy-wide labour supply = `Σ_i xle[i,l]` | cge.jl:277 |
| `y0` | base private GDP = `Σ_i (pva0[i]·xd0[i] − depr[i]·k0[i])` | cge.jl:278 |
| `cd0[i]` | base consumption = `cles[i]·cdtot0` | cge.jl:279 |
| `delta[i]` (Armington share, traded only) | from the base-year cost-minimisation FOC given `pm0/pd0` and the import/domestic ratio `m0/xxd0` | cge.jl:282–286 |
| `x0[i]` | composite-supply value: traded = `pd0·xxd0 + pm0·m0`; non-traded = `pd0·xxd0` | cge.jl:289–295 |
| `ac[i]` (Armington scale, traded only) | chosen so the CES aggregator reproduces `x0[i]` exactly given `delta[i]`, `m0[i]`, `xxd0[i]` | cge.jl:298–300, 329–334 |
| `int0[i]` | base intermediate demand = `Σ_j io[i,j]·xd0[j]` | cge.jl:303 |
| `gamma[i]` (CET share, traded only) | from the base-year export-supply FOC given `pd0/pe0` and `e0/xxd0`; **forced to 0 for any sector with `m0[i]==0`** (a sector with no imports is assumed to have no independent export-supply curve) | cge.jl:306–311 |
| `alphl[l,i]` | labour value share = `(wdist[i,l]·wa0[l]·xle[i,l]) / (pva0[i]·xd0[i])` | cge.jl:314–316 |
| `ad[i]` (Cobb-Douglas TFP) | `xd0[i]` ÷ (the Cobb-Douglas quantity index built from `xllb` — employment with zero cells replaced by 1 to avoid `0^α = NaN` — and `k0`) | cge.jl:318–326 |
| `at[i]` (CET scale, traded only) | chosen so the CET transformation function reproduces `xd0[i]` exactly given `gamma[i]`, `e0[i]`, `xxd0[i]` | cge.jl:337–339 |

### Consistency requirements a new dataset must satisfy

1. **Column sums.** `Σ_i io[i,j] + itax[j] < 1` for every sector `j`, so `pva0[j] > 0` (verified,
   e.g. `cimint`'s `io` column sums to ≈0.604, `itax[cimint]=0.014`). `Σ_i imat[i,j] = 1` for
   every `j` (verified exactly in §3.2).
2. **Budget shares sum to 1.** `Σ_i cles[i] = 1` over sectors with `cles[i] > 0` (verified: the 8
   non-zero `cles` cells in §3.5 sum to exactly `1.00000`); `Σ_i gles[i] = 1` (only `publiques=1`
   here); `Σ_i kio[i] = 1` (verified: sums to exactly `1.00`).
3. **Zero-import ⇒ zero CET share.** `gamma[i]` is forced to 0 wherever `m0[i]==0` (cge.jl:309–311).
   In this dataset that rule is only ever triggered by the two **non-traded** sectors (`construct`,
   `publiques`, both `m0=0`), which are outside `IT` anyway and never use `gamma`; **no traded
   sector in this Cameroon data has `m0[i]==0`**, so the rule is defensive, not actually load-
   bearing here — a future country whose data does zero out imports for a genuinely traded sector
   would exercise it for real.
4. **Non-traded sectors must have zero trade data.** For `i in ITN` (`construct`, `publiques`),
   `m0[i]` and `e0[i]` must be exactly `0` in the sheet (verified) — the model additionally *forces*
   `m[i]=e[i]=0` at solve time regardless of what the sheet says (`fix(...)`, cge.jl:646–652), so a
   non-zero value here would still corrupt `xxd0 = xd0 − e0` in calibration even though the solved
   model discards it.
5. **Non-negative taxes.** `tm0[i]`, `itax[i]` (and the hardcoded `te[i]=0`) are all `≥ 0` throughout.
6. **Employment/wage-distribution zero pairing.** Wherever `xle[i,l]=0`, `alphl[l,i]` calibrates to
   0 regardless of `wdist[i,l]`, and the solved model fixes that `labd[i,l]` cell to exactly 0
   (cge.jl:592–600, 628–653). This fixing logic — like the analogous zero handling for `gd`, `cd`,
   `dst`, `id` (cge.jl:628–635, from `gles[i]==0`, `cles[i]==0`, `dstr[i]==0`, an all-zero `imat`
   row) — is computed dynamically from whichever `Params` are passed in, **not** from a hardcoded
   list of sector names, so it already generalises to any dataset of the same 11×3 shape. Only the
   two non-traded sectors' forced-zero trade cells and `duty≡0` (`te≡0` always) are structural
   regardless of the data.
7. **Base-year replication check.** Solving `calibrate(load_data(path))` cold must reproduce
   `xd0`, `k0`, `ls0`, `y0`, `fsav0` to high precision — this is exactly what
   `test/runtests.jl`'s first testset checks (`rtol=1e-3` on every sector/labour cell,
   `test/runtests.jl:200–217`). If a new country's data fails this check, the SAM itself is
   internally inconsistent (rows/columns don't balance) before any counterfactual is even run.

### Validation steps to run against a candidate dataset

```julia
import Pkg; Pkg.activate(@__DIR__)
include("cge.jl"); using .CGECameroon

raw = load_data(path)            # will throw on a dimension mismatch, nothing else is checked
params = calibrate(raw)          # closed-form algebra; NaN/Inf here usually means a zero divided
                                  # by zero (e.g. m0==0 and xxd0==0 in the same traded sector)
baseline, status, converged, iters, elapsed = solve(params)
converged || error("does not reproduce the base year: $status")
# compare baseline.xd/.k/.ls/.y/.fsav against params.xd0/.k0/.ls0/.y0/.fsav0
```
Or simply run `CGECameroon.example()` (cge.jl:914–928, reproduces baseline + the `sim1` shock and
prints a summary), and the full regression suite: `julia --project=. test/runtests.jl`
(`test/runtests.jl`, 3 testsets: base-year replication, a bit-for-bit regression against
`test/reference.json`, and a warm-start robustness check) plus the optional ~60-scenario grid in
`test/robustness_grid.jl`.

---

## 5. Proposed generic long-format schema (NOT YET IMPLEMENTED)

Today's loader is Cameroon-only: `SEC`/`IT`/`ITN`/`LC` are hardcoded constants and every cell range
in `load_data` assumes exactly this 11-sector/3-labour/2-non-traded layout (cge.jl:57–142). Onboarding
a genuinely different country under the *current* code means literally reusing these 11 symbols
positionally (a new country's "sector 9" must still be its non-traded construction-like sector) or
editing `cge.jl` itself. The following is a **proposal**, not built: a single generic long-format
schema every future country's data would be transformed into, so one generic loader could serve all
of them regardless of sector count or names.

| File | Grain | Columns | Maps onto today's workbook |
|---|---|---|---|
| `sectors.csv` | one row per sector | `code, label, traded, non_traded` | §2 table (`SEC`/`IT`/`ITN`, cge.jl:57–75) |
| `io.csv` | one row per (from,to) pair with a non-zero flow | `from_sector, to_sector, value` | `iotable` sheet, §3.1 (`io[i,j]`) |
| `final_demand.csv` | one row per (sector, component) | `sector, component, value` (`component` ∈ `consumption, government, investment_origin, investment_destination, inventory`) | `cd0`/`cles`, `gd`/`gles`, `id`/`imat`, `dst`/`dstr` in the `miscellaneous` and `imat` sheets, §3.2/§3.5 |
| `factors.csv` | one row per (sector, labour_category) | `sector, labour_category, employment, wage_share` | `employment` (§3.3) + `wagedist` (§3.4) sheets |
| `trade.csv` | one row per traded sector | `sector, imports, exports, world_import_price, tariff_rate, export_tax_rate` | `m0, e0, pwm0, tm0` (`te` today is a hardcoded scalar, not sheet data), §3.5 |
| `scalars.csv` | one row per named scalar | `name, value` (`government_spending, foreign_savings, exchange_rate, saving_rate, capital_stock_by_sector, depreciation_rate_by_sector, elasticities…`) | the hardcoded scalars of cge.jl:234–240 plus the `k`, `depr`, `rhoc`, `rhot`, `eta`, `kio` rows of `miscellaneous` |
| `dataset.toml` | one file | `country, year, currency, source` | not present today at all — see §6, `registry.toml`'s per-dataset metadata is the closest existing analogue |

Every current workbook cell has a home in this schema (see the "Maps onto" column above); nothing
in §3 is dropped, it is only re-shaped from fixed positional ranges into keyed long rows so a
loader can validate shape/keys and support an arbitrary sector list and count. This schema is
**aspirational** — implementing it means writing a new `load_data`-equivalent, generalising
`SEC`/`IT`/`ITN`/`LC` into data read from `sectors.csv`/`factors.csv` rather than module `const`s,
and parameterising every set-indexed `@variable`/`@NLconstraint` in `_solve_once` (cge.jl:510–816)
by those data-driven sets instead of the hardcoded ones — real refactoring work, not yet started
on this branch.

---

## 6. Possible sources for a new country's data (suggestions only — not part of the contract)

`load_data` takes a **file path** (cge.jl:108); it has no notion of where that file's numbers
originally came from. The pointers below are suggestions for where the raw numbers behind §3 or §5
might come from — they are not a required directory layout, not machine-specific paths, and not
something the model or loader depends on.

- **National accounts** (e.g. UN National Accounts official country data, or an analogous NSO
  source) are a natural source for GDP by activity/institutional sector, final-demand components
  (household/government consumption, investment), and taxes-less-subsidies at an aggregate level —
  useful for the `scalars.csv`/economy-wide totals and as a cross-check on `gles`/`cles`
  aggregates, but not sector-by-sector detail.
- **A multi-regional input-output database** (e.g. EMERGING/ICIO-style tables, or GTAP once
  available) is the natural source for the sector × sector **intermediate-use matrix** (`io.csv` /
  the `iotable` sheet), **imports and exports by sector** (`trade.csv`), and a **value-added**
  split by sector — these databases typically classify sectors on ISIC or a proprietary code list
  with many more sectors than this model's 11, so an **aggregation/concordance step** is needed
  (example, illustrative only — the actual mapping depends on the source classification):

  | Model sector | Example source-sector aggregation |
  |---|---|
  | `agsubsist` | ISIC food-crop agriculture codes |
  | `agexpind` | ISIC export/cash-crop agriculture codes |
  | `sylvicult` | ISIC forestry & logging |
  | `indalim` | ISIC food, beverages & tobacco manufacturing |
  | `bienscons` | ISIC light consumer-goods manufacturing (textiles, furniture, etc.) |
  | `biensint` | ISIC intermediate-goods manufacturing (chemicals, basic metals, etc.) |
  | `cimint` | ISIC non-metallic mineral products (cement, construction materials) |
  | `bienscap` | ISIC machinery & equipment manufacturing |
  | `construct` | ISIC construction |
  | `services` | ISIC private services (trade, transport, finance, etc.) |
  | `publiques` | ISIC public administration & other government services |

- **A household survey or labour force survey** is generally needed for **employment by sector ×
  labour category** and the **wage-distribution matrix** (`factors.csv` / the `employment` and
  `wagedist` sheets) — this level of factor-market detail is not typically present in national
  accounts or MRIO databases and needs a dedicated labour-statistics source.
- **Balancing.** Raw data assembled from several sources at different levels of aggregation will
  generally not satisfy §4's consistency rules (columns summing to 1, rows/columns of the SAM
  balancing) out of the box; a standard **RAS** (or similar bi-proportional) balancing step against
  known row/column totals is the usual way to reconcile the assembled matrices before they are
  written into the schema in §5, followed by the same checks listed in §4 (column sums, budget
  shares summing to 1, non-negative taxes, `pva0>0`, and the base-year replication solve).

---

## 7. The web-app registry entry

The "Economic Policy Simulator" web app (`CGE-intelligence`, a separate repo) selects a dataset via
`data/registry.toml`, one `[[dataset]]` table per (model, dataset) pair, parsed by
`src/Registry.jl`. The current Cameroon entry:

```toml
[[dataset]]
id          = "cmr_1979"
model       = "static"
label       = "Cameroon 1979-80 (Condon, Dahl & Devarajan 1987)"
country     = "CMR"
base_year   = 1979
currency    = "CFAF bn"
path        = "data/static/cmr_1979/camdata.xlsx"
loader      = "camdata_xlsx"
n_sectors   = 11
caveats     = ["Historic illustrative SAM", "Single household; fixed exchange rate closure"]
```

`Registry.jl` requires `id, model, label, path, loader, base_year`; every other key
(`n_sectors`, etc.) is carried through as free-form `options`. The app's `StaticCGEAdapter`
(`src/models/StaticCGE.jl:273–280`) dispatches on `loader == "camdata_xlsx"` and calls
`CGECameroon.load_data(entry.path)` then `calibrate(raw)` — i.e. **today's registry loader is this
repo's exact `load_data`/`calibrate` pair**, nothing more generic.

**How a new dataset is registered today:** add a new `[[dataset]]` entry with `model = "static"`,
`loader = "camdata_xlsx"`, and a `path` to a workbook that matches §3's contract **exactly**
(same 5 sheet names, same cell ranges, same 11-sector/3-labour positional ordering, same
`rowz` row order) — i.e. a country's data must literally fit inside Cameroon's positional shape,
not just have "11 sectors and 3 labour categories" in some other order or naming. Registering a
dataset that does not fit this shape requires either (a) hand-editing that country's numbers into
the exact Cameroon-shaped workbook layout, or (b) implementing §5's generic schema plus a new
`loader` value and a matching branch in `StaticCGEAdapter.load_dataset`, at which point the
registry entry's `loader` field would name the new loader instead of `"camdata_xlsx"`.
