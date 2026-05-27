# =============================================================================
# CAMEROON COMPUTABLE GENERAL EQUILIBRIUM (CGE) MODEL
# =============================================================================
# A static, single-country CGE model calibrated to Cameroon's 1979-80 base year.
# Follows the Davis-De Melo-Robinson framework (late 1970s), originally coded
# in GAMS by Condon, Dahl, and Devarajan (World Bank, 1987).
#
# Author : John-Mary Matovu
#          Research Fellow, Research for Transformation and Development (RTD)
#          Kampala, Uganda  |  www.rtdug.org
#
# Dependencies (all free/open-source):
#   XLSX.jl   — reads Excel calibration data
#   JuMP.jl   — algebraic modelling language for optimisation
#   Ipopt.jl  — interior-point NLP solver (COIN-OR)
#
# Usage:
#   include("camcge.jl")
#   # Results stored in `baseline` and `sim1` (Simulation structs)
# =============================================================================

using XLSX
using JuMP
using Ipopt

# =============================================================================
# SETS
# =============================================================================
# SEC  — all 11 production sectors
# IT   — 9 traded sectors (subject to Armington imports and CET exports)
# ITN  — 2 non-traded sectors (domestic supply = domestic demand, no trade)
# lc   — 3 labour categories
# =============================================================================

SEC = [
    :agsubsist,   # Food crops (subsistence agriculture)
    :agexpind,    # Cash crops (export-oriented agriculture)
    :sylvicult,   # Forestry and logging
    :indalim,     # Food processing industry
    :bienscons,   # Consumer goods manufacturing
    :biensint,    # Intermediate goods manufacturing
    :cimint,      # Construction materials
    :bienscap,    # Capital goods manufacturing
    :construct,   # Construction (non-traded)
    :services,    # Private services
    :publiques    # Public services (non-traded)
]

IT  = [:agsubsist, :agexpind, :sylvicult, :indalim, :bienscons,
       :biensint, :cimint, :bienscap, :services]
ITN = [:construct, :publiques]

lc  = [:rural, :urbanunsk, :urbanskil]  # Rural, Urban unskilled, Urban skilled

# =============================================================================
# PARAMETER DICTIONARIES
# Declared empty here; calibrated from base-year data in the section below.
# =============================================================================

# Armington (import aggregation) CES function parameters
delta = Dict()  # Share parameter:  higher delta → more import-intensive sector
ac    = Dict()  # Shift (TFP-like scale) parameter
rhoc  = Dict()  # Substitution exponent: rhoc = 1/sigma_armington - 1

# CET (export allocation) function parameters
rhot  = Dict()  # Transformation exponent: rhot = 1/sigma_cet + 1
at    = Dict()  # Shift parameter
gamma = Dict()  # Share parameter: higher gamma → more export-oriented sector

# Trade and tax parameters
eta   = Dict()  # Export demand elasticity (foreign demand curve slope)
tm0   = Dict()  # Base-year ad-valorem tariff rates (policy instrument)
te    = Dict()  # Export duty rates (zero at base year)
itax  = Dict()  # Indirect (production) tax rates

# Production function parameters
ad    = Dict()  # Cobb-Douglas total-factor-productivity (shift) parameter
alphl = Dict()  # Labour share in value added, by labour type l and sector i
depr  = Dict()  # Depreciation rates by sector

# Consumption and investment parameters
cles  = Dict()  # Private consumption expenditure shares (Engel curve)
gles  = Dict()  # Government consumption expenditure shares
kio   = Dict()  # Investment demand shares by sector of destination
dstr  = Dict()  # Ratio of inventory investment to gross output

# Base-year quantity variables (calibration benchmarks)
m0    = Dict()  # Import volumes           (billion 1979-80 CFAF)
e0    = Dict()  # Export volumes           (billion 1979-80 CFAF)
xd0   = Dict()  # Domestic output          (billion 1979-80 CFAF)
k0    = Dict()  # Capital stocks           (billion 1979-80 CFAF)
id0   = Dict()  # Investment by sector of origin
dst0  = Dict()  # Inventory investment
int0  = Dict()  # Intermediate input demands
xxd0  = Dict()  # Domestic sales = output minus exports
x0    = Dict()  # Composite good supply (Armington aggregate)

# Base-year prices (all normalised near unity at base year)
pwe0  = Dict()  # World export prices (foreign currency)
pwm0  = Dict()  # World import prices (foreign currency)
pd0   = Dict()  # Domestic goods prices
pe0   = Dict()  # Domestic price of exports (= pd0 at base)
pm0   = Dict()  # Domestic price of imports (= pd0 at base)
pva0  = Dict()  # Value-added prices by sector

# Labour variables
wa0   = Dict()  # Base-year wage by labour category (million CFAF / 1000 workers)
ld    = Dict()  # Labour demand derived from profit-maximisation (calibration check)
ls0   = Dict()  # Aggregate labour supply by category (1000 persons)

# Auxiliary calibration helpers
qd    = Dict()  # Intermediate Cobb-Douglas output (used to back out ad)
xllb  = Dict()  # Employment with zeros replaced by 1 (avoids 0^alpha = NaN)
cd0   = Dict()  # Base-year private consumption by sector

# =============================================================================
# BASE-YEAR SCALARS AND WAGES
# =============================================================================

# Average wages (million 1979-80 CFAF per 1000 workers)
[wa0[l] = 0.1 for l in lc]       # default initialisation (overwritten below)
wa0[:rural]     =  0.11
wa0[:urbanunsk] =  0.15678
wa0[:urbanskil] =  1.8657

# Economy-wide scalars
er      =   0.21    # Real exchange rate: CFAF per US dollar
gr0     = 179.0     # Government revenue (billion CFAF)
gdtot0  = 135.03    # Total government consumption (billion CFAF)
cdtot0  = 947.98    # Total private consumption (billion CFAF)
fsav0   =  36.841   # Foreign savings / current-account deficit (billion USD)
y0      =   0.0     # Placeholder; overwritten after value-added prices are computed

# =============================================================================
# DATA LOADING FROM EXCEL
# Path is relative to this file's directory so the model runs from any machine.
# =============================================================================

data_path = joinpath(@__DIR__, "data", "camdata.xlsx")

# Input-output table: io[i, j] = intermediate use of good i per unit output of j
dataio = XLSX.readdata(data_path, "iotable", "B3:L13")
io = Dict()
for i = 1:length(SEC), j = 1:length(SEC)
    io[SEC[i], SEC[j]] = dataio[i, j]
end

# Capital composition matrix: imat[i, j] = share of sector i goods in one unit
# of capital investment destined for sector j  (maps savings → investment demands)
datacap = XLSX.readdata(data_path, "imat", "B3:L13")
imat = Dict()
for i = 1:length(SEC), j = 1:length(SEC)
    imat[SEC[i], SEC[j]] = datacap[i, j]
end

# Wage proportionality factors: wdist[i, l] = relative wage of labour type l in
# sector i.  Integrates heterogeneous labour quality within each labour category.
dataw = XLSX.readdata(data_path, "wagedist", "B3:D13")
wdist = Dict()
for i = 1:length(SEC), j = 1:length(lc)
    wdist[SEC[i], lc[j]] = dataw[i, j]
end

# Employment by sector and labour category (1000 persons)
datae = XLSX.readdata(data_path, "employment", "B3:D13")
xle = Dict()
for i = 1:length(SEC), j = 1:length(lc)
    xle[SEC[i], lc[j]] = datae[i, j]
end

# Miscellaneous parameters (17 rows × 11 sectors).
# Row order: m0, e0, xd0, k, depr, rhoc, rhot, eta, pd0, tm0, itax, cles,
#            gles, kio, dstr, dst, id
datazz = XLSX.readdata(data_path, "miscellaneous", "B3:L19")
rowz = [:m0, :e0, :xd0, :k, :depr, :rhoc, :rhot, :eta, :pd0, :tm0,
        :itax, :cles, :gles, :kio, :dstr, :dst, :id]
zz = Dict()
for i = 1:length(rowz), j = 1:length(SEC)
    zz[rowz[i], SEC[j]] = datazz[i, j]
end

# =============================================================================
# PARAMETER CALIBRATION
# Parameters are chosen so that all model equations hold exactly at base-year
# data.  This guarantees the model replicates the Social Accounting Matrix.
# =============================================================================

[depr[i]  = zz[:depr, i]              for i in SEC]
[rhoc[i]  = 1/zz[:rhoc, i] - 1       for i in SEC]   # Armington exponent (elasticity transform)
[rhot[i]  = 1/zz[:rhot, i] + 1       for i in SEC]   # CET exponent (elasticity transform)
[eta[i]   = zz[:eta,  i]             for i in SEC]
[tm0[i]   = zz[:tm0,  i]             for i in SEC]
[te[i]    = 0                        for i in SEC]    # No export duties at base year
[itax[i]  = zz[:itax, i]             for i in SEC]
[cles[i]  = zz[:cles, i]             for i in SEC]
[gles[i]  = zz[:gles, i]             for i in SEC]
[kio[i]   = zz[:kio,  i]             for i in SEC]
[dstr[i]  = zz[:dstr, i]             for i in SEC]

# Replace zeros in employment matrix with 1 so that 0^alpha does not produce
# NaN in the Cobb-Douglas function.  Sectors where a labour type is absent have
# alphl = 0, making the replacement harmless.
[xllb[i, l] = xle[i, l] + (1 - sign(xle[i, l])) for i in SEC, l in lc]

# Base-year quantities from the miscellaneous sheet
[m0[i]    = zz[:m0,  i]              for i in SEC]
[e0[i]    = zz[:e0,  i]              for i in SEC]
[xd0[i]   = zz[:xd0, i]             for i in SEC]
[k0[i]    = zz[:k,   i]             for i in SEC]
[pd0[i]   = zz[:pd0, i]             for i in SEC]

# At base year all domestic, import, and export prices are equal (normalisation)
pm0 = pd0
pe0 = pd0

# World prices back-calculated from domestic prices and the exchange rate
[pwm0[i] = pm0[i] / ((1 + tm0[i]) * er)  for i in SEC]
[pwe0[i] = pe0[i] / ((1 + te[i])  * er)  for i in SEC]

# Value-added price: residual after subtracting indirect taxes and intermediate costs
[pva0[i] = pd0[i] - sum(io[j,i]*pd0[j] for j in SEC) - itax[i]  for i in SEC]

[xxd0[i] = xd0[i] - e0[i]                for i in SEC]   # Domestic sales = output - exports
[dst0[i] = zz[:dst, i]                   for i in SEC]
[id0[i]  = zz[:id,  i]                   for i in SEC]
[ls0[l]  = sum(xle[i,l] for i in SEC)    for l in lc]    # Aggregate labour supplies
y0       = sum(pva0[i]*xd0[i] - depr[i]*k0[i] for i in SEC)  # Private GDP at base
[cd0[i]  = cles[i]*cdtot0                for i in SEC]

# Armington share parameter (delta): from first-order cost-minimisation condition
[delta[i] = pm0[i]/pd0[i]*(m0[i]/xxd0[i])^(1+rhoc[i])  for i in IT]
[delta[i] = delta[i]/(1 + delta[i])                      for i in IT]

# Composite good volumes at base (value terms for traded/non-traded sectors)
[x0[i] = pd0[i]*xxd0[i] + pm0[i]*m0[i]  for i in IT]
[x0[i] = pd0[i]*xxd0[i]                  for i in ITN]

# Armington shift parameter (ac): ensures armington() holds exactly at base data
[ac[i] = x0[i] / (delta[i]*m0[i]^(-rhoc[i]) + (1-delta[i])*xxd0[i]^(-rhoc[i]))^(-1/rhoc[i])
         for i in IT]

# Leontief intermediate input demands at base year
[int0[i] = sum(io[i,j]*xd0[j] for j in SEC)  for i in SEC]

# CET share parameter (gamma): from export-supply first-order condition at base
[gamma[i] = 1/(1 + pd0[i]/pe0[i]*(e0[i]/xxd0[i])^(rhot[i]-1))  for i in IT]
[gamma[i] = 0  for i in SEC if m0[i] == 0]   # Sectors with zero imports cannot export

# Labour shares in Cobb-Douglas value added: alpha_l = wage bill / value added
[alphl[l, i] = (wdist[i,l]*wa0[l]*xle[i,l]) / (pva0[i]*xd0[i])
               for l in lc, i in SEC]

# Production function TFP parameter (ad): calibrated from base output quantities
[qd[i] = (xllb[i,:rural]^alphl[:rural,i]) *
          (xllb[i,:urbanunsk]^alphl[:urbanunsk,i]) *
          (xllb[i,:urbanskil]^alphl[:urbanskil,i]) *
          (k0[i]^(1 - sum(alphl[l,i] for l in lc)))
          for i in SEC]
[ad[i] = xd0[i] / qd[i]  for i in SEC]

# Recompute ac after ad calibration (ensures ordering consistency)
[x0[i] = pd0[i]*xxd0[i] + pm0[i]*m0[i]  for i in IT]
[ac[i] = x0[i] / (delta[i]*m0[i]^(-rhoc[i]) + (1-delta[i])*xxd0[i]^(-rhoc[i]))^(-1/rhoc[i])
         for i in IT]

# Labour demand at base (derived from profit maximisation; should match xle)
[ld[l] = sum((xd0[i]*pva0[i]*alphl[l,i]) / (wdist[i,l]*wa0[l])
             for i in SEC if wdist[i,l] > 0)
         for l in lc]

# CET shift parameter (at): calibrated so cet() holds at base output/exports
[at[i] = xd0[i] / (gamma[i]*e0[i]^rhot[i] + (1-gamma[i])*xxd0[i]^rhot[i])^(1/rhot[i])
         for i in IT]

# =============================================================================
# CGE MODEL FUNCTION
# Solves for the general equilibrium given current parameter values.
# Returns all 36 endogenous variable arrays.
# The model is posed as a feasibility NLP (minimise dummy constant objective)
# subject to all equilibrium conditions; any feasible point is an equilibrium.
# =============================================================================

function cammodel()
    cgecam = Model(with_optimizer(Ipopt.Optimizer))

    # -------------------------------------------------------------------------
    # DECISION VARIABLES  (all strictly positive via lower bound 1e-6)
    # -------------------------------------------------------------------------
    @variables cgecam begin
        # Prices
        pd[i in SEC]  >= 1e-6, (start = pd0[i])   # Domestic goods price
        pm[i in SEC]  >= 1e-6, (start = pm0[i])   # Domestic price of imports
        pe[i in SEC]  >= 1e-6, (start = pe0[i])   # Domestic price of exports
        pk[i in SEC]  >= 1e-6, (start = pd0[i])   # Capital rental rate by sector
        px[i in SEC]  >= 1e-6, (start = pd0[i])   # Average output price (net of indirect tax)
        p[i in SEC]   >= 1e-6, (start = pd0[i])   # Composite (Armington) goods price
        pva[i in SEC] >= 1e-6, (start = pva0[i])  # Value-added price
        pwm[i in SEC] >= 1e-6, (start = pwm0[i])  # World import price (foreign currency)
        pwe[i in SEC] >= 1e-6, (start = pwe0[i])  # World export price (foreign currency)
        tm[i in SEC]  >= 1e-6, (start = tm0[i])   # Tariff rates

        # Production quantities
        x[i in SEC]   >= 1e-6, (start = x0[i])    # Composite good supply (Armington)
        xd[i in SEC]  >= 1e-6, (start = xd0[i])   # Domestic output by sector
        xxd[i in SEC] >= 1e-6, (start = xxd0[i])  # Domestic sales (output minus exports)
        e[i in SEC]   >= 1e-6, (start = e0[i])    # Exports by sector
        m[i in SEC]   >= 1e-6, (start = m0[i])    # Imports by sector

        # Factors
        k[i in SEC]            >= 1e-6, (start = k0[i])    # Capital stock by sector
        wa[l in lc]            >= 1e-6, (start = wa0[l])   # Economy-wide wage by labour type
        ls[l in lc]            >= 1e-6, (start = ls0[l])   # Labour supply by category
        labd[i in SEC, l in lc] >= 1e-6, (start = xle[i,l]) # Employment (sector × labour type)

        # Demand aggregates
        int[i in SEC]  >= 1e-6, (start = int0[i])   # Intermediate input demand
        cd[i in SEC]   >= 1e-6, (start = cd0[i])    # Private consumption demand
        gd[i in SEC]   >= 1e-6                       # Government consumption demand
        id[i in SEC]   >= 1e-6, (start = id0[i])    # Investment demand by sector of origin
        dst[i in SEC]  >= 1e-6, (start = dst0[i])   # Inventory investment
        y              >= 1e-6, (start = y0)          # Private GDP (value added minus depreciation)
        gr             >= 1e-6, (start = gr0)         # Government revenue
        tariff         >= 1e-6, (start = 76.548)      # Tariff revenue
        indtax         >= 1e-6                         # Indirect tax revenue
        duty           >= 1e-6                         # Export duty revenue
        gdtot          >= 1e-6                         # Total government consumption volume
        mps            >= 1e-6                         # Marginal propensity to save (households)
        hhsav          >= 1e-6                         # Household savings
        govsav         >= 1e-6                         # Government savings (budget surplus/deficit)
        deprecia       >= 1e-6                         # Economy-wide depreciation expenditure
        savings        >= 1e-6                         # Total savings (= total investment)
        fsav           >= 1e-6, (start = fsav0)       # Foreign savings (current-account deficit)
        dk[i in SEC]   >= 1e-6                         # Investment by sector of destination

        # Welfare
        omega                                           # Cobb-Douglas utility (welfare indicator)
    end

    # -------------------------------------------------------------------------
    # EQUILIBRIUM CONDITIONS
    # -------------------------------------------------------------------------
    @NLconstraints cgecam begin

        # --- Price Block ---

        # Import price: domestic = world price × exchange rate × (1 + tariff)
        pmdef[it in IT],          pm[it] == pwm[it]*er*(1 + tm[it])

        # Export price: world price = domestic export price ÷ (er × (1 + export duty))
        pedef[it in IT],          pe[it]*(1 + te[it]) == pwe[it]*er

        # Armington price identity: value of composite supply = domestic sales + imports
        absorption1[i in IT],     p[i]*x[i] == pd[i]*xxd[i] + pm[i]*m[i]

        # Non-traded: composite supply equals domestic sales only (no imports)
        absorption2[itn in ITN],  p[itn]*x[itn] == pd[itn]*xxd[itn]

        # Revenue identity: output revenue = domestic sales revenue + export revenue
        sales[i in SEC],          px[i]*xd[i] == pd[i]*xxd[i] + pe[i]*e[i]

        # Activity price: net-of-tax output price = value added + intermediate cost
        actp[i in SEC],           px[i]*(1 - itax[i]) == pva[i] + sum(io[j,i]*p[j] for j in SEC)

        # Capital goods price: weighted average of composite goods prices (via imat)
        pkdef[i in SEC],          pk[i] == sum(p[j]*imat[j,i] for j in SEC)

        # --- Production and Factor Block ---

        # Cobb-Douglas production with heterogeneous labour categories and fixed capital:
        #   xd[i] = ad[i] * Π_l( labd[i,l]^alphl[l,i] ) * k[i]^(1 - Σ_l alphl[l,i])
        activity[i in SEC],       xd[i] == ad[i] * prod(labd[i,l]^alphl[l,i] for l in lc) *
                                            k[i]^(1 - sum(alphl[l,i] for l in lc))

        # Profit maximisation (Shephard's lemma / FOC for labour):
        #   wage × wdist × employment = value-added price × output × labour share
        profitmax[i in SEC, l in lc], wa[l]*wdist[i,l]*labd[i,l] == xd[i]*pva[i]*alphl[l,i]

        # Labour market clearing: total sectoral demand equals exogenous supply
        lmequil[l in lc],         sum(labd[i,l] for i in SEC) == ls[l]

        # CET function: allocates output between domestic market and exports.
        # Substitution elasticity = 1/(rhot-1).
        cet[it in IT],            xd[it] == at[it] * (gamma[it]*e[it]^rhot[it] +
                                             (1 - gamma[it])*xxd[it]^rhot[it])^(1/rhot[it])

        # Export demand curve: foreign demand is downward-sloping in world price
        edemand[it in IT],        e[it]/e0[it] == (pwe0[it]/pwe[it])^eta[it]

        # Export supply ratio: export/domestic split from CET profit maximisation
        esupply[it in IT],        e[it]/xxd[it] == (pe[it]/pd[it] * (1-gamma[it])/gamma[it])^(1/(rhot[it]-1))

        # Armington function: CES aggregation of imports and domestically-produced goods
        armington[it in IT],      x[it] == ac[it] * (delta[it]*m[it]^(-rhoc[it]) +
                                            (1-delta[it])*xxd[it]^(-rhoc[it]))^(-1/rhoc[it])

        # Cost minimisation FOC: import/domestic ratio from Armington expenditure minimisation
        costmin[it in IT],        m[it]/xxd[it] == (pd[it]/pm[it] * delta[it]/(1-delta[it]))^(1/(1+rhoc[it]))

        # Non-traded sector: domestic sales equal total output (no trade possible)
        xxdsn[itr in ITN],        xxd[itr] == xd[itr]
        xsn[itr in ITN],          x[itr]   == xxd[itr]

        # --- Demand Block ---

        # Intermediate inputs: Leontief (fixed io coefficients), summed over using sectors
        inteq[j in SEC],          int[j] == sum(io[j,i]*xd[i] for i in SEC)

        # Inventory investment: fixed proportion of sectoral gross output
        dsteq[i in SEC],          dst[i] == dstr[i]*xd[i]

        # Private consumption: linear Engel curves, share cles[i] of disposable income
        cdeq[i in SEC],           p[i]*cd[i] == cles[i]*(1 - mps)*y

        # Private GDP: aggregate value added minus economy-wide depreciation
        gdp,                      y == sum(pva[i]*xd[i] for i in SEC) - deprecia

        # Household savings: fixed marginal propensity to save applied to income
        hhsaveq,                  hhsav == mps*y

        # Government revenue: sum of all tax instruments
        greq,                     gr == tariff + duty + indtax

        # Government budget: revenue finances consumption plus budget surplus/deficit
        gruse,                    gr == sum(p[i]*gd[i] for i in SEC) + govsav

        # Government consumption: each commodity gets a fixed share of total volume
        gdeq[i in SEC],           gd[i] == gles[i]*gdtot

        # Tariff revenue: ad-valorem tariffs applied to import values at world prices
        tariffdef,                tariff == sum(tm[it]*m[it]*pwm[it] for it in IT)*er

        # Indirect tax revenue: ad-valorem production taxes on gross output values
        indtaxdef,                indtax == sum(itax[i]*px[i]*xd[i] for i in SEC)

        # Export duty revenue (zero at base since te = 0 everywhere)
        dutydef,                  duty == sum(te[it]*e[it]*pe[it] for it in IT)

        # Aggregate depreciation expenditure (capital consumption)
        depreq,                   deprecia == sum(depr[i]*pk[i]*k[i] for i in SEC)

        # Savings-investment balance: total savings = total investment
        totsav,                   savings == hhsav + govsav + deprecia + fsav*er

        # Investment by sector of destination: savings allocated by fixed shares kio,
        # net of inventory investment (inventories are not part of capital formation)
        prodinv[i in SEC],        pk[i]*dk[i] == kio[i]*savings - kio[i]*sum(dst[j]*p[j] for j in SEC)

        # Investment demand by sector of origin: derived from capital composition matrix
        ieq[i in SEC],            id[i] == sum(imat[i,j]*dk[j] for j in SEC)

        # Current account balance: import bill = export revenues + net foreign savings
        caeq,                     sum(pwm[it]*m[it] for it in IT) == sum(pwe[it]*e[it] for it in IT) + fsav

        # Market clearing: composite supply = all categories of final and intermediate demand
        equil[i in SEC],          x[i] == int[i] + cd[i] + gd[i] + id[i] + dst[i]

        # Welfare indicator: Cobb-Douglas utility over private consumption quantities
        obj,                      omega == prod(cd[i]^cles[i] for i in SEC if cles[i] > 0.0)

        # --- Closure Rules ---
        # Fix exogenous variables; determines which variables adjust to clear markets.
        closurek[i in SEC],   k[i]    == k0[i]      # Capital stocks fixed (static model)
        closurep[i in SEC],   pwm[i]  == pwm0[i]    # World prices fixed (small open economy)
        closurel[l in lc],    ls[l]   == ls0[l]     # Labour supply fixed (no labour-leisure)
        closuret[it in IT],   tm[it]  == tm0[it]    # Tariffs fixed (policy instrument)
        closuref,             fsav    == fsav0       # Foreign savings exogenous
        closuremp,            mps     == 0.09305     # Household saving rate fixed
        closureg,             gdtot   == gdtot0      # Government spending volume fixed
        closurem[itn in ITN], m[itn]  == 0           # Non-traded: no imports
        closurlp,             labd[:publiques, :rural]    == 0  # No rural labour in public sector
        closurla,             labd[:agsubsist, :urbanskil] == 0  # No skilled labour in subsistence farming
        closuree[itn in ITN], e[itn]  == 0           # Non-traded: no exports
    end

    # Feasibility problem: solve for any point satisfying all equilibrium conditions
    @NLobjective(cgecam, Min, 1)

    JuMP.optimize!(cgecam)

    # Extract solution values from JuMP variable containers
    pd  = JuMP.value.(pd);   pm  = JuMP.value.(pm);   pe  = JuMP.value.(pe)
    pk  = JuMP.value.(pk);   px  = JuMP.value.(px);   p   = JuMP.value.(p)
    pva = JuMP.value.(pva);  pwm = JuMP.value.(pwm);  pwe = JuMP.value.(pwe)
    tm  = JuMP.value.(tm)

    x   = JuMP.value.(x);   xd  = JuMP.value.(xd);   xxd = JuMP.value.(xxd)
    e   = JuMP.value.(e);   m   = JuMP.value.(m)

    k    = JuMP.value.(k);   wa   = JuMP.value.(wa);   ls   = JuMP.value.(ls)
    labd = JuMP.value.(labd)

    int  = JuMP.value.(int);  cd   = JuMP.value.(cd);   gd   = JuMP.value.(gd)
    id   = JuMP.value.(id);   dst  = JuMP.value.(dst);  y    = JuMP.value.(y)
    gr   = JuMP.value.(gr)

    tariff   = JuMP.value.(tariff);   indtax   = JuMP.value.(indtax)
    duty     = JuMP.value.(duty);     gdtot    = JuMP.value.(gdtot)
    mps      = JuMP.value.(mps);      hhsav    = JuMP.value.(hhsav)
    govsav   = JuMP.value.(govsav);   deprecia = JuMP.value.(deprecia)
    savings  = JuMP.value.(savings);  fsav     = JuMP.value.(fsav)
    dk       = JuMP.value.(dk)

    return pd, pm, pe, pk, px, p, pva, pwm, pwe, tm,
           x, xd, xxd, e, m, k, wa, ls, labd,
           int, gd, id, dst, y, gr, tariff, indtax, duty, gdtot,
           mps, hhsav, govsav, deprecia, savings, fsav, dk
end

# =============================================================================
# RESULTS CONTAINER
# Bundles all 36 endogenous variable arrays for one model run.
# Any is used because JuMP.value.() returns non-standard DenseAxisArray types.
# =============================================================================

struct Simulation
    pd::Any; pm::Any; pe::Any; pk::Any; px::Any; p::Any; pva::Any
    pwm::Any; pwe::Any; tm::Any; x::Any; xd::Any; xxd::Any; e::Any; m::Any
    k::Any; wa::Any; ls::Any; labd::Any; int::Any; gd::Any; id::Any; dst::Any
    y::Any; gr::Any; tariff::Any; indtax::Any; duty::Any; gdtot::Any; mps::Any
    hhsav::Any; govsav::Any; deprecia::Any; savings::Any; fsav::Any; dk::Any
end

# =============================================================================
# BASELINE SOLUTION
# Solve the model at base-year parameters to reproduce the calibration data.
# =============================================================================

pd, pm, pe, pk, px, p, pva, pwm, pwe, tm,
x, xd, xxd, e, m, k, wa, ls, labd,
int, gd, id, dst, y, gr, tariff, indtax, duty, gdtot,
mps, hhsav, govsav, deprecia, savings, fsav, dk = cammodel()

baseline = Simulation(pd, pm, pe, pk, px, p, pva, pwm, pwe, tm,
                      x, xd, xxd, e, m, k, wa, ls, labd,
                      int, gd, id, dst, y, gr, tariff, indtax, duty, gdtot,
                      mps, hhsav, govsav, deprecia, savings, fsav, dk)

# =============================================================================
# SIMULATION 1 — Agricultural Capital Shock
# Policy: Increase capital stock in the food-crop (subsistence agriculture) sector
# by 10 percent.  Represents a public investment programme or improved access to
# agricultural capital (machinery, irrigation infrastructure, etc.).
# Compare sim1 with baseline to measure the general-equilibrium effects.
# =============================================================================

k0[:agsubsist] = 1.10 * k0[:agsubsist]

pd, pm, pe, pk, px, p, pva, pwm, pwe, tm,
x, xd, xxd, e, m, k, wa, ls, labd,
int, gd, id, dst, y, gr, tariff, indtax, duty, gdtot,
mps, hhsav, govsav, deprecia, savings, fsav, dk = cammodel()

sim1 = Simulation(pd, pm, pe, pk, px, p, pva, pwm, pwe, tm,
                  x, xd, xxd, e, m, k, wa, ls, labd,
                  int, gd, id, dst, y, gr, tariff, indtax, duty, gdtot,
                  mps, hhsav, govsav, deprecia, savings, fsav, dk)
