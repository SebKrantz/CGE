# DATA.md — Input data contract for the Cameroon CGE model

This document specifies what a base-year dataset for this model must contain, in what exact
shape, so that someone can build a new country's input workbook **without reading `cge.jl`**.

`load_data` reads **two** layouts and picks between them automatically, on whether the workbook has
a `sectors` sheet: the **legacy** Cameroon workbook (`data/camdata.xlsx`), documented cell-by-cell
in section 3, and the **generic** workbook, documented in section 5 — the one a new country's data
should be written in, since it carries its own sectors, traded/non-traded partition, labour
categories and economy-wide scalars, for any sector and labour-category count. Section 6 only
suggests where raw numbers might come from — not part of the data contract, since `load_data(path)`
takes a path argument and knows nothing about the file's origin.

All `cge.jl:N` references were written against the pre-`n-sector` file (930 lines) and are now
roughly 200–300 lines low (the generic loader and the data-driven sets were inserted ahead of most
of them); every function, variable and equation name they cite is unchanged. All sheet contents
quoted below were read directly from `data/camdata.xlsx` with `XLSX.jl`.

---

## 1. What the model needs, conceptually

The model is a **base-year social accounting picture for one country**: a single period, single
representative household, single "rest of world" account, **N production sectors** (any number of
which may be non-traded) and **L labour categories**. N, L and the traded/non-traded partition are
read from the workbook (§5); the Cameroon dataset of §2–§3 is one instance of it, with N = 11
(2 non-traded) and L = 3. Concretely, a new country's data must supply:

- An **input-output table**: intermediate use of each sector's output by every other sector
  (sector × sector, "who buys from whom" as a share of the buying sector's gross output).
- **Imports and exports by sector** (for the traded sectors only — the two non-traded sectors have
  neither).
- **Final demands**: household consumption by sector, government consumption by sector,
  investment by sector of *destination* (via fixed shares) translated into investment by sector of
  *origin* via a **capital composition matrix** (`imat`) — i.e. "a unit of investment destined for
  sector X is made up of these shares of goods from sectors A, B, C…".
- **Value added**, split into **labour income by labour category** and a **capital** residual,
  plus a **depreciation rate** on capital.
- **Taxes**: import tariffs, export duties, and indirect/production taxes — all as sector-specific
  ad-valorem rates — plus one economy-wide **household direct-tax rate** `td0` (`0` if the dataset
  has no direct tax, as Cameroon 1979–80 does not).
- The **external account**: foreign savings (the current-account deficit, in foreign currency).
- **Employment by sector × labour category** (levels) and a **wage-distribution matrix** (relative
  wages of each labour category within each sector).
- A handful of **economy-wide scalars**: the exchange rate, total government spending, total
  household consumption, the household saving rate, the direct-tax rate, government revenue and
  base tariff revenue (the last two are solver start values only). See §5's `scalars` sheet.

**Units.** Almost everything is in **billion base-year CFAF** (flows and capital stocks alike).
Employment is in **1000 persons**. Wages are in **million CFAF per 1000 workers**. The exchange
rate and foreign savings are the only quantities in **USD** (billion), because the current-account
identity is expressed at world prices (`caeq`, cge.jl:792). All base-year prices are normalised to
**1** (confirmed: every `pd0` cell in the workbook is exactly `1.0`). None of this is enforced: the
model is unit-agnostic, so a new dataset may state everything in its own currency with `er = 1`
(that is what the synthetic fixture in `test/synthetic_data.jl` does).

**Base year convention.** The model calibrates to exactly one year (here, Cameroon's 1979–80).
There is no time dimension in the data or the equations — "base year" means the single
cross-section the SAM-like data describes, and every calibrated parameter is derived so the model
reproduces that cross-section exactly (`calibrate`, cge.jl:229–347; see §4).

**Scope note.** This is a **single-country, single-household** model: one representative
household (no income deciles, no rural/urban split), one government, and the rest of the world
appears only through the trade/current-account equations — not as a modelled region. Sector count,
sector names, the traded/non-traded partition and the labour categories are **data**, read off the
workbook into `RawData` and carried through `Params` into every index set of the model
(`_solve_once`); the module `const`s `SEC`/`IT`/`ITN`/`LC`/`WA0`/`SCALARS` survive only as the
defaults the legacy loader stamps on `data/camdata.xlsx`, which has no sheets of its own for them.

**Signs.** `fsav0` (foreign savings) is negative for a base-year **trade surplus**; `mps0` is
negative for a **dissaving household**; and the solved `govsav` (government budget balance),
`hhsav` (household saving) and `indtax`/`gr` (net indirect taxes, government revenue) may likewise
come out negative — all of them are declared free in `_solve_attempt`. See §5.5 for the full list
and why each occurs in real data.

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
sectors) cge.jl:73. These `const`s are **the legacy workbook's defaults, not a structural limit**:
`data/camdata.xlsx` carries no `sectors`/`labour` sheets, so `load_data` stamps them onto the
`RawData` it builds from it. A generic workbook (§5) states its own sector list, traded flags and
labour categories, in any order and any count, and every set-indexed variable and equation in the
model follows that data.

| Labour category (`LC`) | Meaning |
|---|---|
| `rural` | Rural labour |
| `urbanunsk` | Urban unskilled labour |
| `urbanskil` | Urban skilled labour |

Source: cge.jl:75. Their base wages `wa0` (million CFAF per 1000 workers: `rural 0.11`,
`urbanunsk 0.15678`, `urbanskil 1.8657`) are the module `const` `WA0`, likewise a legacy-workbook
default; a generic workbook states one `wa0_<code>` per labour category on its `scalars` sheet.

---

## 3. The legacy workbook contract (`data/camdata.xlsx`)

This is the **legacy** layout, taken whenever the workbook has **no `sectors` sheet**; a new
country's dataset should use the generic layout of §5 instead. Five sheets, read by `_load_legacy`
(cge.jl:164–196) with **hardcoded cell ranges and no header-based lookup** — column/row headers in the workbook are free text for a human reader only;
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

Several economy-wide scalars are **not** on any legacy sheet: the module `const`s `WA0` (base wages
by labour category, million CFAF/1000 workers: rural `0.11`, urban-unskilled `0.15678`,
urban-skilled `1.8657`) and `SCALARS` (`er=0.21` CFAF/USD, `gr0=179.0`, `gdtot0=135.03`,
`cdtot0=947.98`, `fsav0=36.841` billion USD, `mps0=0.09305`, `td0=0.0`) supply them, and
`_load_legacy` copies both onto the `RawData` — `calibrate` reads them from there, never from a
literal. A generic workbook states them on its own `scalars` sheet (§5) and needs nothing from these
constants. `te` (export duty rate) is **not read from any sheet** in either layout — hardcoded to
`0.0` for every sector: today's contract has no export-duty row.

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
3. **Every traded sector must actually trade.** `m0[i] > 0` **and** `e0[i] > 0` for every `i in IT`.
   Otherwise the calibration algebra degenerates silently — `m0[i]==0` forces `gamma[i]=0`, which
   makes `ac[i]` NaN and divides by zero in `esupply`; `e0[i]==0` gives `gamma[i]=1`, `at[i]=Inf`
   and an `edemand` that divides by `e0` — so `calibrate` now **raises an error naming the sector**
   instead of returning NaN/Inf and letting the solve fail with `INVALID_MODEL`. Fixes: aggregate
   the sector into a neighbour, mark it `traded = FALSE` (with `m0 = e0 = 0`), or give it a token
   trade flow and rebalance its commodity row.
4. **Non-traded sectors must have zero trade data.** For `i in ITN` (`construct`, `publiques`),
   `m0[i]` and `e0[i]` must be exactly `0` in the sheet (verified) — the model additionally *forces*
   `m[i]=e[i]=0` at solve time regardless of what the sheet says (`fix(...)`, cge.jl:646–652), so a
   non-zero value here would still corrupt `xxd0 = xd0 − e0` in calibration even though the solved
   model discards it.
5. **Non-negative taxes.** `tm0[i]`, `itax[i]` (and the hardcoded `te[i]=0`) are all `≥ 0` throughout
   *in this Cameroon dataset*. Neither is required by the model: `tm` and `indtax`/`gr` are free
   variables, so a zero-rated tariff line or a sector with net production subsidies (`itax[i] < 0`)
   is representable — see §5.5.
6. **Employment/wage-distribution zero pairing.** Wherever `xle[i,l]=0`, `alphl[l,i]` calibrates to
   0 regardless of `wdist[i,l]`, and the solved model fixes that `labd[i,l]` cell to exactly 0
   (cge.jl:592–600, 628–653). This fixing logic — like the analogous zero handling for `gd`, `cd`,
   `dst`, `id` (cge.jl:628–635, from `gles[i]==0`, `cles[i]==0`, `dstr[i]==0`, an all-zero `imat`
   row) — is computed dynamically from whichever `Params` are passed in, **not** from a hardcoded
   list of sector names, so it generalises to any N-sector × L-labour dataset. Only the
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

raw = load_data(path)            # generic layout: errors naming the sheet and the offending
                                  # sector/label on a missing header or a non-numeric cell
params = calibrate(raw)          # closed-form algebra; errors naming any traded sector with
                                  # m0 <= 0 or e0 <= 0 (rule 3 above)
baseline, status, converged, iters, elapsed = solve(params)
converged || error("does not reproduce the base year: $status")
# compare baseline.xd/.k/.ls/.y/.fsav against params.xd0/.k0/.ls0/.y0/.fsav0
```
Or simply run `CGECameroon.example()` (cge.jl:914–928, reproduces baseline + the `sim1` shock and
prints a summary), and the full regression suite: `julia --project=. test/runtests.jl`
(`test/runtests.jl`, 5 testsets: base-year replication, a bit-for-bit regression against
`test/reference.json`, a warm-start robustness check, a generic-workbook round trip of the Cameroon
data, and the synthetic 3-sector dataset of `test/synthetic_data.jl`) plus the optional
~60-scenario grid in `test/robustness_grid.jl`.

---

## 5. The generic workbook format (implemented — use this for a new country)

`load_data(path)` chooses its loader by **whether the workbook has a sheet named `sectors`**: if it
does, the file is read as a *generic* workbook — any number of sectors, any traded/non-traded
partition, any number of labour categories, in **any order** — and if it does not, as the legacy
Cameroon workbook of §3. Nothing else distinguishes the two, and nothing about the model is
Cameroon-specific once the file is loaded: `RawData` carries the sets, and `calibrate`/`_solve_once`
index every parameter, variable and equation by them.

`test/generic_workbook.jl` (`write_generic_workbook`) is the reference *writer* for this format —
the round-trip test writes `data/camdata.xlsx` out in it, with sectors and labour categories in
reverse order, and requires the result to reproduce the legacy load's calibration and solution to
1e-8.

### 5.1 Sheets

Eight sheets, all **lower-case names, matched exactly**; sheet order in the file is irrelevant, and
any additional sheet is ignored.

| Sheet | Shape | Contents |
|---|---|---|
| `sectors` | header row + one row per sector | `code`, `label`, `traded` — defines the sector set, and the traded/non-traded partition from `traded` |
| `labour` | header row + one row per labour category | `code`, `label` — defines the labour-category set |
| `iotable` | labelled N×N | `io[i,j]`: **rows = supplying** sector, **columns = using** sector (§3.1) |
| `imat` | labelled N×N | `imat[i,j]`: **rows = origin** sector, **columns = destination** sector (§3.2) |
| `employment` | labelled N×L | `xle[i,l]`: rows = sector, columns = labour category (§3.3) |
| `wagedist` | labelled N×L | `wdist[i,l]`: rows = sector, columns = labour category (§3.4) |
| `miscellaneous` | labelled 17×N | rows = the 17 parameter names below, columns = sector codes (§3.5) |
| `scalars` | header row + `name`,`value` rows | the economy-wide scalars (§5.4); **optional**, but a new country should always write it |

### 5.2 Labelled matrix sheets (`iotable`, `imat`, `employment`, `wagedist`, `miscellaneous`)

```
     A          B          C          D
1    <corner>   agri       manuf      serv        ← row 1: column headers
2    agri       0.10       0.15       0.02        ← column A: row labels
3    manuf      0.08       0.20       0.10
4    serv       0.05       0.10       0.15
```

- **Cell A1 is a free corner** — anything (a title, or blank); the loader never reads it.
- **Row 1** holds the column labels, **column A** the row labels. Both are looked up **by name**, so
  row and column order are irrelevant and need not agree between sheets or with the `sectors` sheet.
- Labels are **case-sensitive** and whitespace-trimmed, and must match the `code` values of the
  `sectors`/`labour` sheets exactly (unlike the *header* names of §5.3/§5.4, which are matched
  case-insensitively).
- Rows and columns beyond the expected labels are **ignored** (so a total row/column is harmless);
  a **missing** label, a **duplicate** label, or a **blank or non-numeric** data cell in the
  expected block is an error naming the sheet and the offending key.
- Every cell of the N×N / N×L / 17×N block must be present — write explicit `0`s, not blanks.

The `miscellaneous` sheet's 17 row labels, exactly (order irrelevant, all required, one column per
sector — see §3.5 for the meaning and units of each):

`m0`, `e0`, `xd0`, `k`, `depr`, `rhoc`, `rhot`, `eta`, `pd0`, `tm0`, `itax`, `cles`, `gles`, `kio`,
`dstr`, `dst`, `id`

### 5.3 `sectors` and `labour`

```
sectors:                              labour:
     A          B              C           A            B
1    code       label          traded  1   code         label
2    agri       Agriculture    TRUE    2   unskilled    Unskilled labour
3    manuf      Manufacturing  TRUE    3   skilled      Skilled labour
4    serv       Services       FALSE
```

- Row 1 holds the column **headers**, looked up by name and **case-insensitively**; column order is
  irrelevant and extra columns are ignored.
- `sectors` requires `code` and `traded`; `labour` requires `code`. `label` is optional and is read
  by nobody — it is there for the human reading the file.
- `code` becomes a Julia `Symbol` and is the key used everywhere else (matrix labels,
  `miscellaneous` columns, `with_shocks` overrides), so keep it short, ASCII and space-free.
  Duplicate codes are an error; a row with a blank `code` is skipped, so trailing blank rows are
  harmless.
- `traded` accepts an Excel boolean (`TRUE`/`FALSE`), a number (`0`/`1`, non-zero = traded), or a
  string (`true/false`, `t/f`, `yes/no`, `y/n`, `1/0`, any case). At least one sector must be
  traded; every traded sector must have `m0 > 0` and `e0 > 0` (§4 rule 3), and every non-traded one
  must have `m0 = e0 = 0`.

### 5.4 `scalars`

```
     A             B
1    name          value
2    er            1.0
3    gr0           76.6
4    gdtot0        100.0
5    cdtot0        523.26
6    fsav0         -40.0
7    mps0          0.10
8    td0           0.05
9    wa0_unskilled 1.0
10   wa0_skilled   1.0
```

Two columns, `name` and `value` (headers matched case-insensitively, order irrelevant). One row per
scalar, in any order; a blank `name` skips the row. An **unrecognised** name is an error.

| Name | Meaning | Omitted ⇒ |
|---|---|---|
| `er` | exchange rate, domestic currency per unit foreign currency (`1.0` if the whole dataset is in one currency) | `0.21` |
| `gr0` | base government revenue — **solver start value only** | `179.0` |
| `gdtot0` | total government consumption; closure-fixed (`closureg`) | `135.03` |
| `cdtot0` | total private consumption; used by calibration (`cd0 = cles·cdtot0`) | `947.98` |
| `fsav0` | foreign savings = current-account deficit, in **foreign** currency; closure-fixed (`closuref`). **Negative for a trade surplus** | `36.841` |
| `mps0` | household saving rate on disposable income; closure-fixed (`closuremp`) | `0.09305` |
| `td0` | household **direct-tax** rate; closure-fixed (`closuretd`). `0` = no direct tax | `0.0` |
| `tariff0` | base tariff revenue — **solver start value only** | derived: `Σ_traded m0·tm0/(1+tm0)` |
| `wa0_<labour code>` | base wage of that labour category, one row each (`wa0_skilled`, …) | the Cameroon `WA0` value for a code Cameroon has, else `1.0` |

The "omitted" column is the fallback: the whole sheet may be left out, and any single name may be —
the module's Cameroon `const`s (`SCALARS`/`WA0`) fill the gap. **A new country should write every
name explicitly**; silently inheriting Cameroon's exchange rate or saving rate is a much worse
failure than a missing-sheet error.

`wa0` and `wdist` are only ever used as the product `wa0[l]·wdist[i,l]·xle[i,l]`, the base-year wage
bill of cell `(i,l)`, so a dataset that knows wage bills rather than wages and employment can set
`wa0 = 1`, `wdist = 1` and put the **wage bill** in `employment` — that is what
`test/synthetic_data.jl` does.

### 5.5 What the direct tax does, and which balances may be negative

`td0` is the one genuinely new parameter relative to §3's Cameroon contract. It enters three
equations: the household spends and saves out of **disposable** income `(1 − td)·y`
(`cdeq`, `hhsaveq`), and government revenue gains `td·y` (`greq`), so the transfer nets out exactly;
`closuretd` fixes `td = td0`. With `td0 = 0` every equation reduces to the pre-`n-sector` model —
the Cameroon results are unchanged, which is what `test/reference.json` still checks bit-for-bit.

Nine variables are declared **free** in `_solve_attempt` (they used to carry the same `>= 1e-6`
lower bound as every quantity variable): `hhsav`, `govsav`, `fsav`, `td`, plus `mps`, `tm`,
`tariff`, `indtax` and `gr`. Each is a rate or an accounting residual that one equation pins or
defines outright, so a bound on it can only make a satisfiable base year infeasible. All of the
following are therefore ordinary data now, rather than an infeasible model:

| Data feature | Base-year value | Why it happens |
|---|---|---|
| government deficit | `govsav0 < 0` | spending exceeds tax revenue |
| trade surplus | `fsav0 < 0` | exports exceed imports |
| **dissaving household** | **`mps0 < 0`** | consumption exceeds net value added, because remittances / aid / transfers finance the gap and this model carries none of them (86 of the 188 30-sector country databases under `data/`) |
| **zero tariff line** | **`tm0[i] = 0` for a traded `i`** | most tariff schedules have zero-rated lines (186 of 188 databases) |
| **net production subsidies** | **`indtax0 < 0`, hence possibly `gr0 < 0`** | subsidies exceed indirect taxes (12 of 188 databases) |

The `savings = hhsav + govsav + deprecia + fsav·er` total is still bounded positive, as are all
prices and quantities — but their lower bound is now **relative**, a millionth of each variable's
own base-year level, so a legitimately tiny cell (a 1e-11 government-consumption or employment
share of a small economy in billion USD) no longer collides with it. Nothing in the model is an
absolute quantity in the data's units any more; see the README's "Units and scale".

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
