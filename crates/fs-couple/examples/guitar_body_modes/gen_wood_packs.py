"""Write the three shell-ready wood packs (spruce top/braces, rosewood back/ribs, cedar back bars) of the Esteso
guitar body from the committed FPL-GTR-282 seed packs, Bremaud, Gril and Thibaut (2011) Table 4 and the committed
Sitka-spruce loss pack (Qiu 2026).

Stiffness: handbook values (static, about 12 percent moisture) used at audio frequency.
Loss: directional loss factors eta_c (measured on other specimens) are turned into the nine Kelvin-Voigt viscosities the shell owner
consumes by V = D^(1/2) Q D^(1/2) / omega_ref (Q the reduced stiffness, D = diag(eta_L, eta_R, eta_LR)), omega_ref = 2 pi 200 rad/s. Each pack
states omega_ref as constant_loss_reference_omega, so the shell reduction applies the loss factors at every mode's own frequency
(fs_couple::shear_shell::constant_loss). Read without it, the same viscosities are Kelvin-Voigt, and a mode's loss factor
grows as omega / omega_ref.

usage: python3 -I gen_wood_packs.py REPO_ROOT
Prints every derived number with its arithmetic; writes data/matdb/seed-v1/<pack>/{manifest.tsv,properties.tsv}.
"""
import math
import os
import sys

REPO = sys.argv[1].rstrip("/")
SEED = REPO + "/data/matdb/seed-v1/"
OMEGA_MAX = 12566.37  # 2 pi 2000 rad/s: declared use band of static-derived moduli, not a measured band
OMEGA_REF = 2.0 * math.pi * 200.0


def read_scalars(pack):
    out = {}
    for line in open(SEED + pack + "/properties.tsv", encoding="utf-8").read().splitlines():
        f = line.split("\t")
        if f[0] == "scalar":
            out[f[3]] = float(f[4])
    return out


sitka = read_scalars("spruce-sitka-fpl-gtr282")
rose = read_scalars("rosewood-brazilian-fpl-gtr282")
cedar = read_scalars("cedar-western-red-fpl-gtr282")
qiu = read_scalars("spruce-sitka-loss-qiu2026")
hardwoods = ["birch-yellow", "mahogany-african", "mahogany-honduras", "sweetgum", "walnut-black", "yellow-poplar"]
hw = [read_scalars(h + "-fpl-gtr282") for h in hardwoods]
hw_grt = sum(h["grt_over_el"] for h in hw) / len(hw)
hw_nulr = sum(h["nu_lr"] for h in hw) / len(hw)
print("hardwood mean G_RT/E_L = %.5f over %s" % (hw_grt, hardwoods))
print("hardwood mean nu_LR    = %.5f" % hw_nulr)

# Bremaud, Gril, Thibaut (2011) Wood Sci Technol 45:735-754, Table 4, [1,2,3] = [R,T,L]. Stiffness ratios (E_L over ...):
PADAUK = dict(EL_over_G23=9.7, EL_over_E2=11.5, EL_over_G13=9.0, EL_over_E1=7.0)
# Loss rows of the same table (the first column is tan delta_L as printed, in units of 1e-3; the ratios are to tan delta_L):
#   column pairs: tan dG23/tan d3 (shear TL), tan d2/tan d3 (T), tan dG13/tan d3 (shear RL), tan d1/tan d3 (R)
LOSS_ROWS = {
    "resonance spruce": dict(eta_L=6.8e-3, G23=2.32, E2=3.36, G13=2.22, E1=2.92),
    "softwoods": dict(eta_L=7.1e-3, G23=2.14, E2=3.08, G13=2.16, E1=2.90),
    "padauk": dict(eta_L=4.8e-3, G23=2.13, E2=3.5, G13=2.1, E1=3.1),
}

# FPL-GTR-282 Table 4-3b shrinkage green to ovendry, Brazilian rosewood: radial 2.9, tangential 4.6, volumetric 8.5 percent
SV = 8.5
SG_GREEN = rose["specific_gravity"]
rho_rose = 1000.0 * SG_GREEN * 1.12 / (1.0 - SV / 100.0 * (30.0 - 12.0) / 30.0)
print("rosewood rho(12%%) = 1000*%.2f*1.12/(1-%.3f*(30-12)/30) = %.2f kg/m3" % (SG_GREEN, SV / 100.0, rho_rose))

# Spruce loss: the committed Qiu et al. (2026) Sitka pack gives eta_L 0.01424 at 62.36 Hz (beam mode 1) and 0.01134 (beam mode 5, about 830 Hz
# for a 500 x 25 x 3 mm free-free beam); the guitar modes lie between: use the mean, 0.01279, at the geometric midpoint (about 230 Hz).
eta_L_spruce = 0.5 * (qiu["loss_factor_longitudinal"] + qiu["loss_factor_longitudinal_mode5"])
print("spruce eta_L = mean(%.5f, %.5f) = %.5f" % (qiu["loss_factor_longitudinal"], qiu["loss_factor_longitudinal_mode5"], eta_L_spruce))


def etas(row, eta_L):
    r = LOSS_ROWS[row]
    return dict(L=eta_L, R=eta_L * r["E1"], LR=eta_L * r["G13"], LT=eta_L * r["G23"], RT=eta_L * r["G23"])  # RT: G_RT loss unknown, taken as G_LT's


def viscosities(E1, E2, nu12, G12, G13, G23, eta):
    nu21 = nu12 * E2 / E1
    den = 1.0 - nu12 * nu21
    q11, q22, q12, q66 = E1 / den, E2 / den, nu12 * E2 / den, G12
    w = OMEGA_REF
    return {
        "plane_stress_viscosity_11": eta["L"] * q11 / w,
        "plane_stress_viscosity_22": eta["R"] * q22 / w,
        "plane_stress_viscosity_66": eta["LR"] * q66 / w,
        # arithmetic mean of the two in-plane loss factors on the Poisson cross term: the Adams-Bacon split of the stiffness by loss channel
        # (the constant-loss-factor shell path), so the Kelvin-Voigt and constant-eta realisations of one pack agree on every entry
        "plane_stress_viscosity_12": 0.5 * (eta["L"] + eta["R"]) * q12 / w,
        "plane_stress_viscosity_16": 0.0,
        "plane_stress_viscosity_26": 0.0,
        "transverse_shear_viscosity_13": eta["LT"] * G13 / w,
        "transverse_shear_viscosity_23": eta["RT"] * G23 / w,
        "transverse_shear_viscosity_13_23": 0.0,
    }


LOSS_NOTE = (
    "Loss: five dimensionless loss factors loss_factor_1 (L), loss_factor_2 (R), loss_factor_12 (L-R shear), loss_factor_13 (L-T shear) and loss_factor_23 (R-T shear) in the material axes, valid for 2 pi 50 to 2 pi 1000 rad/s (a band declared by the pack author, not measured), for the constant-loss-factor shell path; "
    "and the nine Kelvin-Voigt viscosities built from the same loss factors as eta x stiffness / omega_ref with omega_ref = 2 pi 200 rad/s (1256.64 rad/s): plane_stress_viscosity_11, _22, _66 and the through-thickness shears are eta_L, eta_R, eta_LR, eta_LT, eta_RT times the corresponding stiffness over omega_ref, "
    "the Poisson cross term _12 uses the ARITHMETIC mean of eta_L and eta_R times Q12 / omega_ref (the Adams-Bacon split of the stiffness by loss channel, so both realisations agree on every entry), and the 16 and 26 entries are zero by the orthotropic axes. "
    "Wood's loss factor is nearly constant, so the pack states omega_ref as constant_loss_reference_omega, and the shell reduction applies each loss factor at every mode's own frequency; read as Kelvin-Voigt, the same viscosities would make the loss factor proportional to frequency. "
    "The loss factors are measured on other specimens and are NOT a measurement of the Esteso guitar."
)

COMMON = (
    "Axis convention: material axis 1 = L (grain), axis 2 = R (across the growth rings, in the plane of a quartersawn plate), "
    "thickness = T, so shear_modulus_12 = G_LR, shear_modulus_13 = G_LT and shear_modulus_23 = G_RT, poisson_ratio_12 = nu_LR "
    "(strain along R per strain along L under stress along L). Moduli are static or quasi-static handbook values at about 12 percent "
    "moisture content, used here at audio frequency; the validity band omega 0 to 2 pi 2000 rad/s and T 280 to 320 K is the band in which "
    "this use is DECLARED by the pack author, not a measured band (wood storage moduli typically rise a few percent per decade of frequency). "
    + LOSS_NOTE
)

VISC_ORDER = [
    "plane_stress_viscosity_11", "plane_stress_viscosity_22", "plane_stress_viscosity_66",
    "plane_stress_viscosity_12", "plane_stress_viscosity_16", "plane_stress_viscosity_26",
    "transverse_shear_viscosity_13", "transverse_shear_viscosity_23", "transverse_shear_viscosity_13_23",
]


def fmt(v):
    return "%.9g" % v


def write_pack(pack_id, obs_main, rows, citation, extra_obs=(), fpl_only=True):
    d = SEED + pack_id
    os.makedirs(d, exist_ok=True)
    lines = ["frankensim.matdb-source.v1"]
    lines.append("\t".join(["observation"] + list(obs_main)))
    for o in extra_obs:
        lines.append("\t".join(["observation"] + list(o)))
    for (pid, obs, name, value, unit) in rows:
        lines.append("\t".join(["scalar", pid, obs, name, fmt(value), unit, "constant"]))
        lines.append("\t".join(["uncertainty", pid, "unstated", "-", "-", "-", "-"]))
        lines.append("\t".join(["validity", pid, "T", "280", "320", "K"]))
        if name.startswith("loss_factor_"):
            lines.append("\t".join(["validity", pid, "omega", fmt(LOSS_BAND[0]), fmt(LOSS_BAND[1]), "rad/s"]))
        else:
            lines.append("\t".join(["validity", pid, "omega", "0", fmt(OMEGA_MAX), "rad/s"]))
    open(d + "/properties.tsv", "w", encoding="utf-8").write("\n".join(lines) + "\n")
    man = [
        "frankensim.matdb-manifest.v1",
        "pack_id\t" + pack_id,
        ("redistribution\tpermitted\tUnited States Government work (USDA Forest Service Forest Products Laboratory) and published scientific facts with attribution retained (Bremaud, Gril and Thibaut 2011 ratios; Qiu et al. 2026 loss factors, CC BY 4.0); retain the report number, chapter, table numbers, and URL" if fpl_only else "redistribution\tpermitted\tPublished scientific facts (numbers) with attribution retained: Wood Handbook FPL-GTR-282 is a United States Government work, public domain in the United States; Bremaud, Gril and Thibaut (2011) ratios are transcribed facts; retain the citations"),
        "citation\t" + citation,
        ("license\tPublished-Scientific-Facts-Attribution-Retained"),
        "source\tprimary\tproperties.tsv\tmaterial-tsv-v1",
    ]
    open(d + "/manifest.tsv", "w", encoding="utf-8").write("\n".join(man) + "\n")
    print("wrote", pack_id)
    for r in rows:
        print("   %-34s %s %s" % (r[2], fmt(r[3]), r[4]))


LOSS_FACTOR_ROWS = [("loss_factor_1", "L"), ("loss_factor_2", "R"), ("loss_factor_12", "LR"), ("loss_factor_13", "LT"), ("loss_factor_23", "RT")]
LOSS_BAND = (2.0 * math.pi * 50.0, 2.0 * math.pi * 1000.0)  # rad/s: declared by the pack author, not a measured band


def rows_for(obs, vals, visc, eta):
    base = [("density", "kg/m^3"), ("young_modulus_1", "Pa"), ("young_modulus_2", "Pa"), ("poisson_ratio_12", "1"),
            ("shear_modulus_12", "Pa"), ("shear_modulus_13", "Pa"), ("shear_modulus_23", "Pa")]
    out = [(k, obs[k] if isinstance(obs, dict) else obs, k, vals[k], u) for k, u in base]
    out += [(k, obs[k] if isinstance(obs, dict) else obs, k, visc[k], "Pa*s") for k in VISC_ORDER]
    # dimensionless loss factors in the material axes for the constant-loss-factor shell path (same source observation as the viscosities)
    loss_obs = obs[VISC_ORDER[0]] if isinstance(obs, dict) else obs
    out += [(k, loss_obs, k, eta[axis], "1") for k, axis in LOSS_FACTOR_ROWS]
    # The angular frequency at which the nine viscosities above were written: a reader applies each loss factor at every
    # mode's own frequency with it (fs_couple::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA); without it they are Kelvin-Voigt.
    out.append(("constant_loss_reference_omega", loss_obs, "constant_loss_reference_omega", OMEGA_REF, "rad/s"))
    return out


# ------------------------------------------------------------------ spruce (Sitka stiffness, Qiu loss, Bremaud resonance-spruce loss ratios)
EL = sitka["young_modulus_longitudinal"] * 1e6
sv = dict(density=sitka["density"], young_modulus_1=EL, young_modulus_2=sitka["er_over_el"] * EL, poisson_ratio_12=sitka["nu_lr"],
          shear_modulus_12=sitka["glr_over_el"] * EL, shear_modulus_13=sitka["glt_over_el"] * EL, shear_modulus_23=sitka["grt_over_el"] * EL)
e_sp = etas("resonance spruce", eta_L_spruce)
vis_sp = viscosities(sv["young_modulus_1"], sv["young_modulus_2"], sv["poisson_ratio_12"], sv["shear_modulus_12"], sv["shear_modulus_13"], sv["shear_modulus_23"], e_sp)
print("spruce eta", e_sp)
obs_sp = {k: "wh-sitka-shell" for k in sv}
obs_sp.update({k: "qiu-bremaud-loss" for k in VISC_ORDER})
write_pack(
    "guitar-spruce-sitka-fpl-gtr282-shell",
    ("wh-sitka-shell", "Sitka spruce (Picea sitchensis) clear straight-grained wood, soundboard stand-in",
     "FPL-GTR-282 chapter 5: Table 5-3a bending MOE 10800 MPa times 1.10 (the chapter's shear-deflection correction) gives E_L = 11880 MPa; Table 5-1 ratios E_R/E_L 0.078, G_LR/E_L 0.064, G_LT/E_L 0.061, G_RT/E_L 0.003; Table 5-2 nu_LR 0.372; Table 5-3a SG 0.40 gives density = 1000 x 0.40 x 1.12 = 448 kg/m3. The absolute E_R, G_LR, G_LT and G_RT are those ratios times E_L (arithmetic done by the pack author from the committed spruce-sitka-fpl-gtr282 pack)",
     "The species of the Esteso 1964 top is UNKNOWN (the plan says only spruce); Sitka is the FPL species with a full orthotropic set, not a claim about the specimen. " + COMMON),
    rows_for(obs_sp, sv, vis_sp, e_sp),
    "Forest Products Laboratory (2021), Wood Handbook, FPL-GTR-282, Chapter 5 Mechanical Properties of Wood, Tables 5-1, 5-2 and 5-3a; https://research.fs.usda.gov/treesearch/62200; Qiu et al. (2026) doi:10.3390/f17060669 (committed pack spruce-sitka-loss-qiu2026); Bremaud, Gril and Thibaut (2011) doi:10.1007/s00226-010-0393-8 Table 4; derived by the pack author",
    extra_obs=[("qiu-bremaud-loss", "Sitka spruce loss factors (Qiu et al. 2026) with the directional loss ratios of resonance spruce (Bremaud, Gril and Thibaut 2011 Table 4)",
                "eta_L = mean of the committed Qiu pack's loss_factor_longitudinal %.5f (beam mode 1, 62.36 Hz) and loss_factor_longitudinal_mode5 %.5f (beam mode 5, about 830 Hz for a 500 x 25 x 3 mm free-free beam) = %.5f; directional ratios from the resonance-spruce median row: eta_R/eta_L = 2.92, eta_LR/eta_L = 2.22, eta_LT/eta_L = 2.32; eta_RT is not given by either source and is taken equal to eta_LT" % (qiu["loss_factor_longitudinal"], qiu["loss_factor_longitudinal_mode5"], eta_L_spruce),
                "SURROGATE for the Esteso top: the loss factors are Sitka beams cut along the grain at 12 percent moisture, the ratios belong to a different spruce set. Bremaud's median resonance-spruce eta_L is 0.0068, half of Qiu's: report Q with both.")],
    fpl_only=True,
)

# ------------------------------------------------------------------ rosewood (Brazilian: FPL E_L and density, Bremaud padauk ratios, FPL hardwood means; padauk loss)
ELr = rose["young_modulus_longitudinal"] * 1e6
rv = dict(density=rho_rose, young_modulus_1=ELr, young_modulus_2=ELr / PADAUK["EL_over_E1"], poisson_ratio_12=hw_nulr,
          shear_modulus_12=ELr / PADAUK["EL_over_G13"], shear_modulus_13=ELr / PADAUK["EL_over_G23"], shear_modulus_23=hw_grt * ELr)
e_ro = etas("padauk", LOSS_ROWS["padauk"]["eta_L"])
vis_ro = viscosities(rv["young_modulus_1"], rv["young_modulus_2"], rv["poisson_ratio_12"], rv["shear_modulus_12"], rv["shear_modulus_13"], rv["shear_modulus_23"], e_ro)
print("rosewood eta", e_ro)
obs_ro = {"density": "wh-rosewood-density", "young_modulus_1": "wh-rosewood-el", "young_modulus_2": "bremaud-padauk-ratios",
          "poisson_ratio_12": "wh-hardwood-mean", "shear_modulus_12": "bremaud-padauk-ratios", "shear_modulus_13": "bremaud-padauk-ratios",
          "shear_modulus_23": "wh-hardwood-mean"}
obs_ro.update({k: "bremaud-padauk-loss" for k in VISC_ORDER})
write_pack(
    "guitar-rosewood-brazilian-fpl-gtr282-bremaud2011-shell",
    ("wh-rosewood-el", "Brazilian rosewood (Dalbergia nigra), the Rio rosewood of the plan's back and sides",
     "FPL-GTR-282 Table 5-5a, 12 percent moisture: bending MOE 13000 MPa times 1.10 gives E_L = 14300 MPa",
     "E_L only. The Wood Handbook prints no elastic ratios for any rosewood. " + COMMON),
    rows_for(obs_ro, rv, vis_ro, e_ro),
    "Forest Products Laboratory (2021) Wood Handbook FPL-GTR-282 Tables 4-3b and 5-5a, and Bremaud I, Gril J, Thibaut B (2011) Anisotropy of wood vibrational properties: dependence on grain angle and review of literature data, Wood Science and Technology 45(4):735-754, doi:10.1007/s00226-010-0393-8, Table 4 (African padauk row); https://hal.science/hal-00804242",
    extra_obs=[
        ("wh-rosewood-density", "Brazilian rosewood (Dalbergia nigra)",
         "FPL-GTR-282 Table 5-5a specific gravity 0.80 on the OVENDRY-WEIGHT OVER GREEN-VOLUME basis, and Table 4-3b shrinkage green to ovendry radial 2.9, tangential 4.6, volumetric 8.5 percent. Density at 12 percent moisture = 1000 x SG x 1.12 / (1 - S_v x (30 - 12)/30) with S_v = 0.085 (linear shrinkage between the fibre saturation point at 30 percent and ovendry, Wood Handbook chapter 4): 896 / 0.949 = %.1f kg/m3" % rho_rose,
         "Derived by the pack author with the Wood Handbook's own shrinkage rule; the frankensim FPL seed pack prints no density for this species because its SG basis differs from Table 5-3a."),
        ("bremaud-padauk-ratios", "African padauk (Pterocarpus soyauxii), a Leguminosae of reduced anisotropy used as the surrogate for Rio rosewood's anisotropy",
         "Bremaud, Gril and Thibaut (2011) Table 4, padauk row, [1,2,3] = [R,T,L]: E_L/G_TL = 9.7, E_L/E_T = 11.5, E_L/G_RL = 9, E_L/E_R = 7 (the table marks part of the row as values adjusted by comparing two species of one genus); the same paper compares measured Rio rosewood (Dalbergia nigra, Yano et al. 1995) grain-angle dependence of storage modulus and loss factor with calculation from the PADAUK parameters and reports good agreement. Used here: E_R = E_L/7, G_LR = E_L/9, G_LT = E_L/9.7 with E_L from the Wood Handbook",
         "SURROGATE: ratios are padauk's, not measured on Rio rosewood. Sensitivity to this choice must be reported with any result."),
        ("wh-hardwood-mean", "mean of six hardwoods with full Table 5-1/5-2 sets in the frankensim seed packs: yellow birch, African mahogany, Honduras mahogany, sweetgum, black walnut, yellow-poplar",
         "FPL-GTR-282 Tables 5-1 and 5-2: mean G_RT/E_L = %.5f and mean nu_LR = %.5f over those six species; used for the two constants (G_RT, nu_LR) that neither the Wood Handbook nor Bremaud et al. give for rosewood" % (hw_grt, hw_nulr),
         "SURROGATE: a mean over other hardwoods, not a rosewood measurement. G_RT enters only through-thickness shear across the grain of a thin plate; nu_LR enters plate stiffness weakly."),
        ("bremaud-padauk-loss", "African padauk loss factors as the surrogate for Rio rosewood",
         "Bremaud, Gril and Thibaut (2011) Table 4, padauk row: tan delta_L printed 4.8 (units of 1e-3, consistent with E''/rho 86 MPa over E'/rho 17.3 GPa), tan dG23/tan d3 = 2.13, tan d2/tan d3 = 3.5, tan dG13/tan d3 = 2.1, tan d1/tan d3 = 3.1: eta_L = 0.0048, eta_R = 0.01488, eta_LR = 0.01008, eta_LT = 0.010224; eta_RT taken equal to eta_LT",
         "SURROGATE and a low value: the same paper states that Rio rosewood and padauk both have abnormally low damping; a dynamic-mechanical test of Dalbergia nigra elsewhere gives tan delta 0.022 under different conditions (frequency and temperature not comparable). Report Q with both."),
    ],
    fpl_only=False,
)

# ------------------------------------------------------------------ cedar (western red cedar, FPL; softwood loss row)
ELc = cedar["young_modulus_longitudinal"] * 1e6
cv = dict(density=cedar["density"], young_modulus_1=ELc, young_modulus_2=cedar["er_over_el"] * ELc, poisson_ratio_12=cedar["nu_lr"],
          shear_modulus_12=cedar["glr_over_el"] * ELc, shear_modulus_13=cedar["glt_over_el"] * ELc, shear_modulus_23=cedar["grt_over_el"] * ELc)
e_ce = etas("softwoods", LOSS_ROWS["softwoods"]["eta_L"])
vis_ce = viscosities(cv["young_modulus_1"], cv["young_modulus_2"], cv["poisson_ratio_12"], cv["shear_modulus_12"], cv["shear_modulus_13"], cv["shear_modulus_23"], e_ce)
print("cedar eta", e_ce)
obs_ce = {k: "wh-cedar-shell" for k in cv}
obs_ce.update({k: "bremaud-softwood-loss" for k in VISC_ORDER})
write_pack(
    "guitar-cedar-western-red-fpl-gtr282-shell",
    ("wh-cedar-shell", "Western redcedar (Thuja plicata) clear straight-grained wood, stand-in for the plan's unnamed cedar",
     "FPL-GTR-282 chapter 5: Table 5-3a bending MOE 7700 MPa times 1.10 gives E_L = 8470 MPa; Table 5-1 ratios E_R/E_L 0.081, G_LR/E_L 0.087, G_LT/E_L 0.086, G_RT/E_L 0.005; Table 5-2 nu_LR 0.378; Table 5-3a SG 0.32 gives density = 1000 x 0.32 x 1.12 = 358.4 kg/m3. Absolute values are those ratios times E_L",
     "The plan names cedar for the back bars and the neck without a species; Spanish cedar (Cedrela odorata) is the usual guitar cedar and is not in the Wood Handbook's elastic tables, so western redcedar is a SURROGATE of lower density. " + COMMON),
    rows_for(obs_ce, cv, vis_ce, e_ce),
    "Forest Products Laboratory (2021), Wood Handbook, FPL-GTR-282, Chapter 5, Tables 5-1, 5-2 and 5-3a; https://research.fs.usda.gov/treesearch/62200; Bremaud, Gril and Thibaut (2011) Table 4 softwoods row; derived by the pack author",
    extra_obs=[("bremaud-softwood-loss", "mean softwoods (Bremaud, Gril and Thibaut 2011 Table 4) as the loss surrogate for cedar",
                "median softwood row: tan delta_L 7.1 (units of 1e-3), tan dG23/tan d3 = 2.14, tan d2/tan d3 = 3.08, tan dG13/tan d3 = 2.16, tan d1/tan d3 = 2.90: eta_L = 0.0071, eta_R = 0.02059, eta_LR = 0.015336, eta_LT = 0.015194; eta_RT taken equal to eta_LT",
                "SURROGATE: no cedar loss measurement is used.")],
    fpl_only=True,
)
