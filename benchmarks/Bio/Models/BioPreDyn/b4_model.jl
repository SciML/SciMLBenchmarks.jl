#=
BioPreDyn Benchmark B4: CHO Cell Metabolism Model
Translated from Matlab (BioPreDynBenchFiles v15102014B)

35 states, 32 reactions, 117 parameters
Log-linear elasticity kinetics with Michaelis-Menten for protein synthesis (r19)
=#

# ---------------------------------------------------------------------------- #
#  State variable names (indices 1..35)
# ---------------------------------------------------------------------------- #
const B4_STATE_NAMES = [
    "beta_D_Glucose_f",             # 1
    "L_Lactate_f",                  # 2
    "L_Leucine_f",                  # 3
    "L_Methionine_f",               # 4
    "Productprotein_f",             # 5
    "L_Glutamate_m",                # 6
    "NAD_m",                        # 7
    "2_Oxoglutarate_m",             # 8
    "NADH_m",                       # 9
    "L_Glutamine_c",                # 10
    "ADP_c",                        # 11
    "L_Glutamate_c",                # 12
    "ATP_c",                        # 13
    "L_Aspartate_m",                # 14
    "Oxaloacetate_m",               # 15
    "L_Malate_m",                   # 16
    "CoQH_radical_m",               # 17
    "H_out_m",                      # 18
    "CoQ_m",                        # 19
    "H_in_m",                       # 20
    "Pyruvate_c",                   # 21
    "Phosphoenolpyruvate_c",        # 22
    "NADH_c",                       # 23
    "D_Glycerate_3_phosphate_c",    # 24
    "NAD_c",                        # 25
    "beta_D_Glucose_c",             # 26
    "L_Malate_c",                   # 27
    "2_Oxoglutarate_c",             # 28
    "L_Aspartate_c",                # 29
    "ATP_m",                        # 30
    "Orthophosphate_m",             # 31
    "ADP_m",                        # 32
    "L_Leucine_c",                  # 33
    "L_Methionine_c",               # 34
    "Fermenter_volume",             # 35
]

# ---------------------------------------------------------------------------- #
#  Reaction names (indices 1..32)
# ---------------------------------------------------------------------------- #
const B4_REACTION_NAMES = [
    "EXCHANGE_glc",   "EXCHANGE_lac",   "EXCHANGE_leu",   "EXCHANGE_met",
    "EXCHANGE_product", "Glud1",        "Glul",           "Got2",
    "Mdh2",           "ND1_a",          "Pklr",           "Subset0",
    "Subset1",        "Subset2",        "Subset26",       "Subset3",
    "Subset35",       "Subset37",       "Subset4",        "Subset5",
    "Ucp",            "adencarr",       "akgcarr",        "aspglucarr",
    "atpase",         "atpase1",        "dicarr",         "feed_glc",
    "feed_leu",       "feed_met",       "glucarr",        "mitphocarr",
]

# ---------------------------------------------------------------------------- #
#  Nominal parameter vector (117 elements, all stored as |value|)
# ---------------------------------------------------------------------------- #
const B4_NOMINAL_PARAMETERS = Float64[
    # Km parameters for Subset4 (Michaelis-Menten), par[1:8]
    1000.0,   # 1  Km_Subset4_D-Glycerate_3-phosphate_c
    1000.0,   # 2  Km_Subset4_NAD_c
    1000.0,   # 3  Km_Subset4_L-Glutamate_c
    1000.0,   # 4  Km_Subset4_L-Leucine_c
    1000.0,   # 5  Km_Subset4_L-Methionine_c
    1000.0,   # 6  Km_Subset4_L-Aspartate_c
    1000.0,   # 7  Km_Subset4_L-Glutamine_c
    1000.0,   # 8  Km_Subset4_ATP_c
    # Elasticity parameters, par[9:117]
    1.0,      # 9   e_substrate_Glud1_2-Oxoglutarate_m
    1.0,      # 10  e_substrate_Glud1_NADH_m
    0.7,      # 11  e_product_Glud1_L-Glutamate_m        (negated in model)
    0.7,      # 12  e_product_Glud1_NAD_m                 (negated)
    0.7,      # 13  e_substrate_Glul_L-Glutamate_c
    0.7,      # 14  e_substrate_Glul_ATP_c
    0.2,      # 15  e_product_Glul_L-Glutamine_c          (negated)
    0.2,      # 16  e_product_Glul_ADP_c                  (negated)
    1.0,      # 17  e_substrate_Got2_Oxaloacetate_m
    1.0,      # 18  e_substrate_Got2_L-Glutamate_m
    0.7,      # 19  e_product_Got2_2-Oxoglutarate_m       (negated)
    0.7,      # 20  e_product_Got2_L-Aspartate_m          (negated)
    0.05,     # 21  e_activator_Got2_L-Malate_m
    1.0,      # 22  e_substrate_Mdh2_L-Malate_m
    0.5,      # 23  e_substrate_Mdh2_NAD_m
    0.7,      # 24  e_product_Mdh2_Oxaloacetate_m         (negated)
    0.5,      # 25  e_product_Mdh2_NADH_m                 (negated)
    2.0,      # 26  e_substrate_ND1_a_NADH_m
    0.7,      # 27  e_substrate_ND1_a_CoQ_m
    0.7,      # 28  e_substrate_ND1_a_H_in_m
    2.0,      # 29  e_product_ND1_a_NAD_m                 (negated)
    0.2,      # 30  e_product_ND1_a_CoQH_radical_m        (negated)
    0.2,      # 31  e_product_ND1_a_H_out_m               (negated)
    1.0,      # 32  e_substrate_Pklr_Phosphoenolpyruvate_c
    0.7,      # 33  e_substrate_Pklr_ADP_c
    0.2,      # 34  e_product_Pklr_Pyruvate_c             (negated)
    0.21,     # 35  e_product_Pklr_ATP_c                  (negated)
    1.0,      # 36  e_substrate_Subset0_Oxaloacetate_m
    0.7,      # 37  e_substrate_Subset0_NAD_m
    1.0,      # 38  e_substrate_Subset0_Pyruvate_c
    0.2,      # 39  e_product_Subset0_2-Oxoglutarate_m    (negated)
    0.2,      # 40  e_product_Subset0_NADH_m              (negated)
    0.05,     # 41  e_activator_Subset0_ADP_m
    0.01,     # 42  e_inhibitor_Subset0_ATP_m             (negated)
    0.7,      # 43  e_substrate_Subset1_NAD_c
    1.0,      # 44  e_substrate_Subset1_beta-D-Glucose_c
    0.5,      # 45  e_product_Subset1_NADH_c              (negated)
    1.0,      # 46  e_product_Subset1_D-Glycerate_3-phosphate_c (negated)
    0.04,     # 47  e_activator_Subset1_ADP_c
    0.01,     # 48  e_inhibitor_Subset1_ADP_c             (negated)
    0.01,     # 49  e_inhibitor_Subset1_ATP_c             (negated)
    0.01,     # 50  e_inhibitor_Subset1_Phosphoenolpyruvate_c (negated)
    1.0,      # 51  e_substrate_Subset2_2-Oxoglutarate_c
    1.0,      # 52  e_substrate_Subset2_L-Aspartate_c
    0.5,      # 53  e_substrate_Subset2_NADH_c
    1.0,      # 54  e_product_Subset2_L-Glutamate_c       (negated)
    1.0,      # 55  e_product_Subset2_L-Malate_c          (negated)
    0.5,      # 56  e_product_Subset2_NAD_c               (negated)
    0.1,      # 57  e_activator_Subset2_L-Malate_c
    0.1,      # 58  e_inhibitor_Subset2_L-Glutamine_c     (negated)
    0.5,      # 59  e_substrate_Subset26_ADP_m
    1.0,      # 60  e_substrate_Subset26_Phosphoenolpyruvate_c
    1.0,      # 61  e_substrate_Subset26_L-Malate_m
    1.0,      # 62  e_product_Subset26_Oxaloacetate_m     (negated)
    0.3,      # 63  e_product_Subset26_ATP_m              (negated)
    1.0,      # 64  e_product_Subset26_L-Malate_c         (negated)
    1.0,      # 65  e_substrate_Subset3_Pyruvate_c
    0.7,      # 66  e_substrate_Subset3_NADH_c
    0.4,      # 67  e_product_Subset3_L-Lactate_f         (negated)
    0.7,      # 68  e_product_Subset3_NAD_c               (negated)
    1.0,      # 69  e_substrate_Subset35_D-Glycerate_3-phosphate_c
    0.7,      # 70  e_product_Subset35_Phosphoenolpyruvate_c (negated)
    0.7,      # 71  e_substrate_Subset37_H_in_m
    2.0,      # 72  e_substrate_Subset37_CoQH_radical_m
    0.2,      # 73  e_product_Subset37_H_out_m            (negated)
    2.0,      # 74  e_product_Subset37_CoQ_m              (negated)
    1.0,      # 75  e_substrate_Subset5_2-Oxoglutarate_m
    0.7,      # 76  e_substrate_Subset5_NAD_m
    0.7,      # 77  e_substrate_Subset5_CoQ_m
    0.2,      # 78  e_substrate_Subset5_Orthophosphate_m
    0.7,      # 79  e_substrate_Subset5_ADP_m
    0.5,      # 80  e_product_Subset5_L-Malate_m          (negated)
    0.21,     # 81  e_product_Subset5_NADH_m              (negated)
    0.2,      # 82  e_product_Subset5_CoQH_radical_m      (negated)
    0.2,      # 83  e_product_Subset5_ATP_m               (negated)
    0.01,     # 84  e_inhibitor_Subset5_Oxaloacetate_m    (negated)
    2.0,      # 85  e_substrate_adencarr_ADP_c
    2.0,      # 86  e_substrate_adencarr_ATP_m
    2.0,      # 87  e_product_adencarr_ADP_m              (negated)
    2.0,      # 88  e_product_adencarr_ATP_c              (negated)
    2.0,      # 89  e_substrate_akgcarr_2-Oxoglutarate_m
    2.0,      # 90  e_substrate_akgcarr_L-Malate_c
    2.0,      # 91  e_product_akgcarr_2-Oxoglutarate_c    (negated)
    2.0,      # 92  e_product_akgcarr_L-Malate_m          (negated)
    1.0,      # 93  e_substrate_aspglucarr_L-Aspartate_m
    1.0,      # 94  e_substrate_aspglucarr_L-Glutamate_c
    1.0,      # 95  e_product_aspglucarr_L-Aspartate_c    (negated)
    1.0,      # 96  e_product_aspglucarr_L-Glutamate_m    (negated)
    1.0,      # 97  e_substrate_atpase_ATP_c
    0.5,      # 98  e_product_atpase_ADP_c                (negated)
    0.7,      # 99  e_substrate_atpase1_ADP_m
    0.7,      # 100 e_substrate_atpase1_Orthophosphate_m
    2.0,      # 101 e_substrate_atpase1_H_out_m
    0.2,      # 102 e_product_atpase1_ATP_m               (negated)
    2.0,      # 103 e_product_atpase1_H_in_m              (negated)
    1.0,      # 104 e_substrate_dicarr_L-Malate_c
    1.0,      # 105 e_substrate_dicarr_Orthophosphate_m
    0.7,      # 106 e_product_dicarr_L-Malate_m           (negated)
    1.0,      # 107 e_substrate_feed_glc_beta-D-Glucose_f
    0.7,      # 108 e_product_feed_glc_beta-D-Glucose_c   (negated)
    1.0,      # 109 e_substrate_feed_leu_L-Leucine_f
    0.7,      # 110 e_product_feed_leu_L-Leucine_c        (negated)
    1.0,      # 111 e_substrate_feed_met_L-Methionine_f
    0.7,      # 112 e_product_feed_met_L-Methionine_c     (negated)
    2.0,      # 113 e_substrate_glucarr_L-Glutamate_m
    2.0,      # 114 e_product_glucarr_L-Glutamate_c       (negated)
    0.7,      # 115 e_substrate_mitphocarr_H_out_m
    0.2,      # 116 e_product_mitphocarr_Orthophosphate_m (negated)
    0.2,      # 117 e_product_mitphocarr_H_in_m           (negated)
]

# ---------------------------------------------------------------------------- #
#  Initial conditions (35 states)
# ---------------------------------------------------------------------------- #
const B4_U0 = Float64[5.0, 1.0, 5.0, 5.0, ones(30)..., 6.0]

# ---------------------------------------------------------------------------- #
#  Time span
# ---------------------------------------------------------------------------- #
const B4_TSPAN = (0.0, 300.0)

# ---------------------------------------------------------------------------- #
#  Pre-computed model constants (allocated once, used inside the ODE)
# ---------------------------------------------------------------------------- #
struct B4Constants{T}
    rho::T
    vol_1::T
    vol_2::T
    vol_ratio::T        # vol_1 / vol_2
    inv_rho::T          # 1 / rho
    rss::Vector{T}      # 32 reference steady-state rates
    css::Vector{T}      # 34 steady-state concentrations
    fp::T               # feed pulse magnitude
end

function build_b4_constants(::Type{T} = Float64) where {T}
    rho   = T(141.4710605)
    vol_1 = T(0.9)
    vol_2 = T(0.1)
    vr    = vol_1 / vol_2   # 9.0

    rss = zeros(T, 32)
    rss[1]  = T(83401.9654177705)
    rss[2]  = T(127877.37232994902)
    rss[3]  = T(603.4140712332287)
    rss[4]  = T(603.4140712332287)
    rss[5]  = T(5.028450593613567)
    rss[6]  = T(2413.6562849217303) * vr
    rss[7]  = T(603.414071180172)
    rss[8]  = T(40133.38664791273) * vr
    rss[9]  = T(74836.04665510054) * vr
    rss[10] = T(180754.26889005644) * vr
    rss[11] = T(164390.2745505474)
    rss[12] = T(36512.902220891134)
    rss[13] = T(166803.93083553354)
    rss[14] = T(39529.972576829896)
    rss[15] = T(1810.242213701644)
    rss[16] = T(127877.37232994902)
    rss[17] = T(166200.51676427448)
    rss[18] = T(216060.34296802033) * vr
    rss[19] = T(5.028450593613567)
    rss[20] = T(35306.0740779408) * vr
    rss[21] = T(0.0) * vr    # Ucp: zero rate
    rss[22] = T(532531.1670427525)
    rss[23] = T(38926.55850561767)
    rss[24] = T(40133.38664808034)
    rss[25] = T(689982.1797742102)
    rss[26] = T(495414.85075147304) * vr
    rss[27] = T(2413.6562849329107)
    rss[28] = T(83401.9654177705)
    rss[29] = T(603.4140712332287)
    rss[30] = T(603.4140712332287)
    rss[31] = T(2413.656284931604)
    rss[32] = T(533134.5811139438)

    css = zeros(T, 34)
    css[1]  = T(100000); css[2]  = T(1);     css[3]  = T(1000);  css[4]  = T(1000)
    css[5]  = T(30);     css[6]  = T(1000);  css[7]  = T(1000);  css[8]  = T(1000)
    css[9]  = T(1000);   css[10] = T(1000);  css[11] = T(1000);  css[12] = T(1000)
    css[13] = T(3000);   css[14] = T(1000);  css[15] = T(1000);  css[16] = T(1000)
    css[17] = T(100);    css[18] = T(100);   css[19] = T(100);   css[20] = T(100)
    css[21] = T(1000);   css[22] = T(1000);  css[23] = T(1000);  css[24] = T(1000)
    css[25] = T(1000);   css[26] = T(1000);  css[27] = T(1000);  css[28] = T(1000)
    css[29] = T(1000);   css[30] = T(3000);  css[31] = T(100000); css[32] = T(1000)
    css[33] = T(1000);   css[34] = T(1000)

    return B4Constants{T}(rho, vol_1, vol_2, vr, one(T) / rho, rss, css, zero(T))
end

const B4_CONSTANTS = build_b4_constants()

# ---------------------------------------------------------------------------- #
#  Map from the 117-element parameter vector to the 128-element elasticity
#  array used internally.  This exactly mirrors the Matlab index juggling.
# ---------------------------------------------------------------------------- #
@inline function _build_elasticities!(e::AbstractVector, par::AbstractVector)
    # par[1:8] -> e[1:8]  (Km values for Subset4)
    @inbounds for i in 1:8
        e[i] = par[i]
    end
    # Fixed zeros / ones
    @inbounds begin
        e[9]  = 0.0
        e[10] = 1.0
        e[11] = 0.0
        e[12] = 0.0
        e[13] = 1.0
    end
    # par[9:35] -> e[14:40]
    @inbounds for i in 9:35
        e[i + 5] = par[i]
    end
    @inbounds e[41] = 0.0
    # par[36:41] -> e[42:47]
    @inbounds for i in 36:41
        e[i + 6] = par[i]
    end
    @inbounds e[48] = 0.0
    # par[42:83] -> e[49:90]
    @inbounds for i in 42:83
        e[i + 7] = par[i]
    end
    @inbounds begin
        e[91] = 0.0
        e[92] = 0.0
        e[93] = par[84]
        e[94] = 0.7
        e[95] = -0.2
    end
    # par[85:117] -> e[96:128]
    @inbounds for i in 85:117
        e[i + 11] = par[i]
    end

    # Negate the product / inhibitor elasticities
    @inbounds begin
        e[16]  = -par[11];  e[17]  = -par[12]
        e[20]  = -par[15];  e[21]  = -par[16]
        e[24]  = -par[19];  e[25]  = -par[20]
        e[29]  = -par[24];  e[30]  = -par[25]
        e[34]  = -par[29];  e[35]  = -par[30];  e[36] = -par[31]
        e[39]  = -par[34];  e[40]  = -par[35]
        e[45]  = -par[39];  e[46]  = -par[40]
        e[49]  = -par[42]
        e[52]  = -par[45];  e[53]  = -par[46]
        e[55]  = -par[48];  e[56]  = -par[49];  e[57] = -par[50]
        e[61]  = -par[54];  e[62]  = -par[55];  e[63] = -par[56]
        e[65]  = -par[58]
        e[69]  = -par[62];  e[70]  = -par[63];  e[71] = -par[64]
        e[74]  = -par[67];  e[75]  = -par[68]
        e[77]  = -par[70]
        e[80]  = -par[73];  e[81]  = -par[74]
        e[87]  = -par[80];  e[88]  = -par[81];  e[89] = -par[82];  e[90] = -par[83]
        e[93]  = -par[84]
        e[98]  = -par[87];  e[99]  = -par[88]
        e[102] = -par[91];  e[103] = -par[92]
        e[106] = -par[95];  e[107] = -par[96]
        e[109] = -par[98]
        e[113] = -par[102]; e[114] = -par[103]
        e[117] = -par[106]
        e[119] = -par[108]
        e[121] = -par[110]
        e[123] = -par[112]
        e[125] = -par[114]
        e[127] = -par[116]; e[128] = -par[117]
    end
    return nothing
end

# ---------------------------------------------------------------------------- #
#  ODE right-hand side
# ---------------------------------------------------------------------------- #
"""
    biopredyn_b4!(du, u, p, t)

BioPreDyn B4 CHO cell metabolism ODE.

- `u`: 35-element state vector (dimensionless, scaled by css)
- `p`: 117-element parameter vector (absolute values; signs applied internally)
- `t`: time
- `du`: derivative vector (mutated in-place)

Uses the pre-built `B4_CONSTANTS` for rss, css, rho, and volume fractions.
"""
function biopredyn_b4!(du, u, p, t)
    # Unpack constants
    mc = B4_CONSTANTS
    rss = mc.rss
    css = mc.css
    inv_rho = mc.inv_rho
    vr  = mc.vol_ratio      # vol_1/vol_2 = 9.0
    fp  = mc.fp              # feed pulse magnitude (0.0)

    # Build the 128-element elasticity array from the 117 parameters
    e = zeros(eltype(u), 128)
    _build_elasticities!(e, p)

    # Shorthand for states
    c = u   # c[i] = u[i], 1-indexed

    # Log of concentrations (clamped: if c[i] < 0, use 0.0 instead of log)
    clogk = ntuple(i -> (c[i] < zero(eltype(u)) ? zero(eltype(u)) : log(c[i])), Val(35))

    # Feed pulse: active only during t in (5, 10)
    feed_val = (t > 5.0 && t < 10.0) ? fp : zero(eltype(u))

    # ---- Reactions ----
    # r[1..5] = 0 (exchange reactions, set to zero)
    r1  = zero(eltype(u))
    r2  = zero(eltype(u))
    r3  = zero(eltype(u))
    r4  = zero(eltype(u))
    r5  = zero(eltype(u))

    # Glud1 (r6)
    r6  = rss[6]  * (1.0 + (e[14]*clogk[8] + e[15]*clogk[9]) +
                             (e[16]*clogk[6] + e[17]*clogk[7]))

    # Glul (r7)
    r7  = rss[7]  * (1.0 + (e[18]*clogk[12] + e[19]*clogk[13]) +
                             (e[20]*clogk[10] + e[21]*clogk[11]))

    # Got2 (r8)
    r8  = rss[8]  * (1.0 + ((e[22]*clogk[15] + e[23]*clogk[6]) +
                             (e[24]*clogk[8]  + e[25]*clogk[14])) +
                             e[26]*clogk[16])

    # Mdh2 (r9)
    r9  = rss[9]  * (1.0 + (e[27]*clogk[16] + e[28]*clogk[7]) +
                             (e[29]*clogk[15] + e[30]*clogk[9]))

    # ND1_a (r10)
    r10 = rss[10] * (1.0 + ((e[31]*clogk[9]  + e[32]*clogk[19]) + e[33]*clogk[20]) +
                            ((e[34]*clogk[7]  + e[35]*clogk[17]) + e[36]*clogk[18]))

    # Pklr (r11)
    r11 = rss[11] * (1.0 + ((e[37]*clogk[22] + e[38]*clogk[11]) +
                             (e[39]*clogk[21] + e[40]*clogk[13])) +
                             e[41]*clogk[13])

    # Subset0 (r12)
    r12 = rss[12] * (1.0 + ((((e[42]*clogk[15] + e[43]*clogk[7]) + e[44]*clogk[21]) +
                              (e[45]*clogk[8]  + e[46]*clogk[9])) + e[47]*clogk[32]) +
                             e[48]*clogk[9] + e[49]*clogk[30])

    # Subset1 (r13)
    r13 = rss[13] * (1.0 + ((((e[50]*clogk[25] + e[51]*clogk[26]) +
                              (e[52]*clogk[23] + e[53]*clogk[24])) + e[54]*clogk[11]) +
                             e[55]*clogk[11]) + e[56]*clogk[13] + e[57]*clogk[22])

    # Subset2 (r14)
    r14 = rss[14] * (1.0 + (((e[58]*clogk[28] + e[59]*clogk[29]) + e[60]*clogk[23]) +
                             ((e[61]*clogk[12] + e[62]*clogk[27]) + e[63]*clogk[25])) +
                             e[64]*clogk[27] + e[65]*clogk[10])

    # Subset26 (r15)
    r15 = rss[15] * (1.0 + ((e[66]*clogk[32] + e[67]*clogk[22]) + e[68]*clogk[16]) +
                            ((e[69]*clogk[15] + e[70]*clogk[30]) + e[71]*clogk[27]))

    # Subset3 (r16)
    r16 = rss[16] * (1.0 + (e[72]*clogk[21] + e[73]*clogk[23]) +
                             (e[74]*clogk[2]  + e[75]*clogk[25]))

    # Subset35 (r17)
    r17 = rss[17] * (1.0 + e[76]*clogk[24] + e[77]*clogk[22])

    # Subset37 (r18)
    r18 = rss[18] * (1.0 + (e[78]*clogk[20] + e[79]*clogk[17]) +
                             (e[80]*clogk[18] + e[81]*clogk[19]))

    # Subset4 (r19) -- Michaelis-Menten multi-substrate kinetics
    # Reference saturation (at steady state concentrations css)
    Km = @view e[1:8]
    # The 8 substrates of Subset4 map to states: 24, 25, 12, 33, 34, 29, 10, 13
    sub_idx = (24, 25, 12, 33, 34, 29, 10, 13)

    # Numerator_ss = prod(css[sub_idx[k]] / Km[k])
    # Denom_ss     = prod(1 + css[sub_idx[k]] / Km[k])
    # Numerator     = prod(c[sub_idx[k]] * css[sub_idx[k]] / Km[k])
    # Denom         = prod(1 + c[sub_idx[k]] * css[sub_idx[k]] / Km[k])
    # r19 = rss[19] / (Numerator_ss / Denom_ss) * (Numerator / Denom)
    num_ss  = one(eltype(u))
    den_ss  = one(eltype(u))
    num_cur = one(eltype(u))
    den_cur = one(eltype(u))
    @inbounds for k in 1:8
        si = sub_idx[k]
        ratio_ss  = css[si] / Km[k]
        ratio_cur = c[si] * css[si] / Km[k]
        num_ss  *= ratio_ss
        den_ss  *= (1.0 + ratio_ss)
        num_cur *= ratio_cur
        den_cur *= (1.0 + ratio_cur)
    end
    r19 = rss[19] / (num_ss / den_ss) * (num_cur / den_cur)

    # Subset5 (r20)
    r20 = rss[20] * (1.0 + (((((((e[82]*clogk[8]  + e[83]*clogk[7]) + e[84]*clogk[19]) +
                                  e[85]*clogk[31]) + e[86]*clogk[32]) +
                                (((e[87]*clogk[16] + e[88]*clogk[9]) + e[89]*clogk[17]) +
                                  e[90]*clogk[30])) + e[91]*clogk[31]) +
                             e[92]*clogk[9]) + e[93]*clogk[15])

    # Ucp (r21) -- rss[21] = 0, so r21 = 0 always
    r21 = rss[21] * (1.0 + e[94]*clogk[18] + e[95]*clogk[20])

    # adencarr (r22)
    r22 = rss[22] * (1.0 + (e[96]*clogk[11]  + e[97]*clogk[30]) +
                             (e[98]*clogk[32]  + e[99]*clogk[13]))

    # akgcarr (r23)
    r23 = rss[23] * (1.0 + (e[100]*clogk[8]  + e[101]*clogk[27]) +
                             (e[102]*clogk[28] + e[103]*clogk[16]))

    # aspglucarr (r24)
    r24 = rss[24] * (1.0 + (e[104]*clogk[14] + e[105]*clogk[12]) +
                             (e[106]*clogk[29] + e[107]*clogk[6]))

    # atpase (r25)
    r25 = rss[25] * (1.0 + e[108]*clogk[13] + e[109]*clogk[11])

    # atpase1 (r26)
    r26 = rss[26] * (1.0 + ((e[110]*clogk[32] + e[111]*clogk[31]) + e[112]*clogk[18]) +
                             (e[113]*clogk[30] + e[114]*clogk[20]))

    # dicarr (r27)
    r27 = rss[27] * (1.0 + (e[115]*clogk[27] + e[116]*clogk[31]) + e[117]*clogk[16])

    # feed_glc (r28)
    r28 = rss[28] * (1.0 + e[118]*clogk[1]  + e[119]*clogk[26])

    # feed_leu (r29)
    r29 = rss[29] * (1.0 + e[120]*clogk[3]  + e[121]*clogk[33])

    # feed_met (r30)
    r30 = rss[30] * (1.0 + e[122]*clogk[4]  + e[123]*clogk[34])

    # glucarr (r31)
    r31 = rss[31] * (1.0 + e[124]*clogk[6]  + e[125]*clogk[12])

    # mitphocarr (r32)
    r32 = rss[32] * (1.0 + e[126]*clogk[18] + (e[127]*clogk[31] + e[128]*clogk[20]))

    # ---- ODEs ----
    # Fermenter volume (state 35)
    du[35] = feed_val
    dV = du[35]

    # beta-D-Glucose_f
    du[1]  = (r1  - r28 * inv_rho - c[1]  * css[1]  / c[35] * dV) / css[1]
    # L-Lactate_f
    du[2]  = (-r2  + r16 * inv_rho - c[2]  * css[2]  / c[35] * dV) / css[2]
    # L-Leucine_f
    du[3]  = (r3  - r29 * inv_rho - c[3]  * css[3]  / c[35] * dV) / css[3]
    # L-Methionine_f
    du[4]  = (r4  - r30 * inv_rho - c[4]  * css[4]  / c[35] * dV) / css[4]
    # Productprotein_f
    du[5]  = (-r5  + r19 * inv_rho - c[5]  * css[5]  / c[35] * dV) / css[5]
    # L-Glutamate_m
    du[6]  = (r6  - r8  + r24 * vr - r31 * vr) / css[6]
    # NAD_m
    du[7]  = (r6  - r9  + r10 - 2.0 * r12 * vr - r20) / css[7]
    # 2-Oxoglutarate_m
    du[8]  = (-r6  + r8  + r12 * vr - r20 - r23 * vr) / css[8]
    # NADH_m
    du[9]  = (-r6  + r9  - r10 + 2.0 * r12 * vr + r20) / css[9]
    # L-Glutamine_c
    du[10] = (r7  - 120.0 * r19) / css[10]
    # ADP_c
    du[11] = (r7  - r11 + 1260.0 * r19 - r22 + r25) / css[11]
    # L-Glutamate_c
    du[12] = (-r7  + r14 - 240.0 * r19 - r24 + r31) / css[12]
    # ATP_c
    du[13] = (-r7  + r11 - 1260.0 * r19 + r22 - r25) / css[13]
    # L-Aspartate_m
    du[14] = (r8  - r24 * vr) / css[14]
    # Oxaloacetate_m
    du[15] = (-r8  + r9  - r12 * vr + r15 * vr) / css[15]
    # L-Malate_m
    du[16] = (-r9  - r15 * vr + r20 + r23 * vr + r27 * vr) / css[16]
    # CoQH_radical_m
    du[17] = (2.0 * r10 - 2.0 * r18 + 2.0 * r20) / css[17]
    # H_out_m
    du[18] = (4.0 * r10 + 6.0 * r18 - r21 - 3.0 * r26 - r32 * vr) / css[18]
    # CoQ_m
    du[19] = (-2.0 * r10 + 2.0 * r18 - 2.0 * r20) / css[19]
    # H_in_m
    du[20] = (-4.0 * r10 - 6.0 * r18 + r21 + 3.0 * r26 + r32 * vr) / css[20]
    # Pyruvate_c
    du[21] = (r11 - r12 - r16) / css[21]
    # Phosphoenolpyruvate_c
    du[22] = (-r11 - r15 + r17) / css[22]
    # NADH_c
    du[23] = (r13 - r14 - r16 + 120.0 * r19) / css[23]
    # D-Glycerate_3-phosphate_c
    du[24] = (r13 - r17 - 120.0 * r19) / css[24]
    # NAD_c
    du[25] = (-r13 + r14 + r16 - 120.0 * r19) / css[25]
    # beta-D-Glucose_c
    du[26] = (-0.5 * r13 + r28) / css[26]
    # L-Malate_c
    du[27] = (r14 + r15 - r23 - r27) / css[27]
    # 2-Oxoglutarate_c
    du[28] = (-r14 + 120.0 * r19 + r23) / css[28]
    # L-Aspartate_c
    du[29] = (-r14 - 120.0 * r19 + r24) / css[29]
    # ATP_m
    du[30] = (r15 * vr + r20 - r22 * vr + r26) / css[30]
    # Orthophosphate_m
    du[31] = (-r20 - r26 - r27 * vr + r32 * vr) / css[31]
    # ADP_m
    du[32] = (-r15 * vr - r20 + r22 * vr - r26) / css[32]
    # L-Leucine_c
    du[33] = (-120.0 * r19 + r29) / css[33]
    # L-Methionine_c
    du[34] = (-120.0 * r19 + r30) / css[34]

    return nothing
end
