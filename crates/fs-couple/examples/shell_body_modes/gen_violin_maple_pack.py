"""Write the shell-ready sycamore maple pack of the il Cannone violin body (back, ribs, bridge) from the committed Sonderegger et al. (2013) and Laemmlein et al. (2020) packs.

Stiffness (one tree, 626 kg/m3 at about 12 percent moisture, 20 C / 65 percent RH): E_L = the static bending modulus, E_R = the mean of the static radial tension and compression moduli,
nu_LR = the mean Poisson ratio, G_LR, G_LT and G_RT = the ULTRASONIC shear moduli (1 MHz: NOT audio band, declared in the pack text).
Loss: eta_L and eta_LR are Laemmlein's violin-back maple plates; eta_R, eta_LT and eta_RT are UNKNOWN for maple and are set equal to eta_LR as a declared surrogate (mail 5280 of the
orchestrator: T and RT/LT shear loss are unknown in every open source). The nine Kelvin-Voigt viscosities the shell owner consumes are V = D^(1/2) Q D^(1/2) / omega_ref with
omega_ref = 2 pi 200 rad/s, exactly as in the guitar packs (gen_wood_packs.py); they do not enter the frequencies. A constant loss factor per mode is recovered by the example as
C_hyst(omega) = (omega_ref / omega) C_KV, and the packs of the wood-loss lane replace these numbers.

usage: python3 -I gen_violin_maple_pack.py REPO_ROOT OUT_ROOT
Reads REPO_ROOT/data/matdb/seed-v1/{maple-sycamore-sonderegger2013,maple-sycamore-loss-lammlein2020}/properties.tsv and writes OUT_ROOT/violin-maple-sycamore-sonderegger2013-shell/{manifest.tsv,properties.tsv}.
"""
import math
import os
import sys

REPO, OUT = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
SEED = REPO + "/data/matdb/seed-v1/"
OMEGA_MAX = 12566.37
OMEGA_REF = 2.0 * math.pi * 200.0
PACK = "violin-maple-sycamore-sonderegger2013-shell"


def scalars(pack):
    # Keyed by property name AND by row id: the Sonderegger pack has two density rows (density-normal-climate and density-ovendry), and the
    # property name alone kept whichever came last (the oven-dry 563 kg/m3, against moduli measured at about 12 percent moisture).
    out = {}
    for line in open(SEED + pack + "/properties.tsv", encoding="utf-8").read().splitlines():
        f = line.split("\t")
        if f[0] == "scalar":
            out[f[3]] = float(f[4])
            out["row:" + f[1]] = float(f[4])
    return out


so = scalars("maple-sycamore-sonderegger2013")
la = scalars("maple-sycamore-loss-lammlein2020")
E_L = so["young_modulus_longitudinal_bending"] * 1e6
E_R = 0.5 * (so["young_modulus_radial_tension"] + so["young_modulus_radial_compression"]) * 1e6
vals = dict(density=so["row:density-normal-climate"], young_modulus_1=E_L, young_modulus_2=E_R, poisson_ratio_12=so["nu_lr"], shear_modulus_12=so["shear_modulus_lr_ultrasonic"] * 1e6,
            shear_modulus_13=so["shear_modulus_lt_ultrasonic"] * 1e6, shear_modulus_23=so["shear_modulus_rt_ultrasonic"] * 1e6)
eta = dict(L=la["loss_factor_longitudinal"], LR=la["loss_factor_shear_lr"])
eta.update(R=eta["LR"], LT=eta["LR"], RT=eta["LR"])
print("E_L %.4g Pa (static bending), E_R %.4g Pa (mean of radial tension and compression), nu_LR %.3f, G_LR %.4g, G_LT %.4g, G_RT %.4g Pa (ultrasonic)" % (
    E_L, E_R, vals["poisson_ratio_12"], vals["shear_modulus_12"], vals["shear_modulus_13"], vals["shear_modulus_23"]))
print("eta", eta)

nu12 = vals["poisson_ratio_12"]
nu21 = nu12 * E_R / E_L
den = 1.0 - nu12 * nu21
q11, q22, q12, q66 = E_L / den, E_R / den, nu12 * E_R / den, vals["shear_modulus_12"]
w = OMEGA_REF
visc = {
    "plane_stress_viscosity_11": eta["L"] * q11 / w,
    "plane_stress_viscosity_22": eta["R"] * q22 / w,
    "plane_stress_viscosity_66": eta["LR"] * q66 / w,
    "plane_stress_viscosity_12": math.sqrt(eta["L"] * eta["R"]) * q12 / w,
    "plane_stress_viscosity_16": 0.0,
    "plane_stress_viscosity_26": 0.0,
    "transverse_shear_viscosity_13": eta["LT"] * vals["shear_modulus_13"] / w,
    "transverse_shear_viscosity_23": eta["RT"] * vals["shear_modulus_23"] / w,
    "transverse_shear_viscosity_13_23": 0.0,
}
VISC_ORDER = ["plane_stress_viscosity_11", "plane_stress_viscosity_22", "plane_stress_viscosity_66", "plane_stress_viscosity_12", "plane_stress_viscosity_16", "plane_stress_viscosity_26",
              "transverse_shear_viscosity_13", "transverse_shear_viscosity_23", "transverse_shear_viscosity_13_23"]
BASE = [("density", "kg/m^3"), ("young_modulus_1", "Pa"), ("young_modulus_2", "Pa"), ("poisson_ratio_12", "1"), ("shear_modulus_12", "Pa"), ("shear_modulus_13", "Pa"), ("shear_modulus_23", "Pa")]

AXES = ("Axis convention: material axis 1 = L (grain), axis 2 = R (across the growth rings, in the plane of a quartersawn plate), thickness = T, so shear_modulus_12 = G_LR, "
        "shear_modulus_13 = G_LT and shear_modulus_23 = G_RT, poisson_ratio_12 = nu_LR (strain along R per strain along L under stress along L; the source's subscripts were converted to this "
        "order in the committed Sonderegger pack). The validity band omega 0 to 2 pi 2000 rad/s and T 280 to 320 K is the band in which this use is DECLARED by the pack author, not a measured band.")
OBS = [
    ("sond-shell", "Sycamore maple (Acer pseudoplatanus L.), clear wood of ONE tree of Eastern Switzerland, 626 kg/m3 at about 12 percent moisture, 20 C / 65 percent RH",
     "Sonderegger et al. (2013) Eur J Wood Prod 71:91-99 as committed in the pack maple-sycamore-sonderegger2013: E_L = the static bending modulus %.0f MPa; E_R = the mean of the static radial tension %.0f and "
     "compression %.0f MPa = %.0f MPa; nu_LR = the mean over tension and compression %.2f; G_LR %.0f, G_LT %.0f and G_RT %.0f MPa are ULTRASONIC values (1 MHz transverse waves, 10 mm cubes)" % (
         so["young_modulus_longitudinal_bending"], so["young_modulus_radial_tension"], so["young_modulus_radial_compression"], E_R / 1e6, so["nu_lr"], so["shear_modulus_lr_ultrasonic"],
         so["shear_modulus_lt_ultrasonic"], so["shear_modulus_rt_ultrasonic"]),
     "The shear moduli are not audio-band (ultrasonic moduli exceed static ones) and the single tree is not the violin's wood. A violin-grade (curly, resonance) maple of the committed Spycher pack has "
     "E_L about 9700 to 10000 MPa and E_R about 1650 to 1700 MPa dynamically: carry that as the sensitivity. " + AXES),
    ("laemmlein-maple-loss", "Sycamore maple violin-back wood, 30 quarter-sawn plates 95 x 60 x 2.6 mm, unvarnished",
     "Laemmlein et al. (2020) Holzforschung 74:765-776 as committed in the pack maple-sycamore-loss-lammlein2020: loss_factor_longitudinal %.4f and loss_factor_shear_lr %.4f; eta_L = %.4f, eta_LR = %.4f; "
     "eta_R, eta_LT and eta_RT are UNKNOWN for maple and set equal to eta_LR" % (la["loss_factor_longitudinal"], la["loss_factor_shear_lr"], eta["L"], eta["LR"]),
     "SURROGATE for three of the five directional loss factors. The nine viscosities are Kelvin-Voigt values V = eta x stiffness / omega_ref with omega_ref = 2 pi 200 rad/s (1256.64 rad/s); the example "
     "recovers a constant loss factor per mode by omega_ref / omega. They do not enter any frequency."),
]


def fmt(v):
    return "%.9g" % v


os.makedirs(OUT + "/" + PACK, exist_ok=True)
lines = ["frankensim.matdb-source.v1"]
for o in OBS:
    lines.append("\t".join(("observation",) + o))
for k, u in BASE:
    for row in (("scalar", k, "sond-shell", k, fmt(vals[k]), u, "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", "0", fmt(OMEGA_MAX), "rad/s")):
        lines.append("\t".join(row))
for k in VISC_ORDER:
    for row in (("scalar", k, "laemmlein-maple-loss", k, fmt(visc[k]), "Pa*s", "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", "0", fmt(OMEGA_MAX), "rad/s")):
        lines.append("\t".join(row))
open(OUT + "/" + PACK + "/properties.tsv", "w", encoding="utf-8").write("\n".join(lines) + "\n")
man = ["frankensim.matdb-manifest.v1", "pack_id\t" + PACK,
       "redistribution\tpermitted\tNumerical measurement values retained as scientific facts with attribution; no article text, figure or scanned table is reproduced; derived by the pack author from the committed packs maple-sycamore-sonderegger2013 and maple-sycamore-loss-lammlein2020 (whose own terms apply to the loss values: publisher copyright on the article, ETH Research Collection copy In Copyright - Non-Commercial Use Permitted)",
       "citation\tSonderegger W., Martienssen A., Nitsche C., Ozyhar T., Kaliske M., Niemz P. (2013), Eur J Wood Prod 71:91-99, doi:10.1007/s00107-012-0641-8; Laemmlein S.L., Van Damme B., Mannes D., Schwarze F.W.M.R., Burgert I. (2020), Holzforschung 74(8):765-776, doi:10.1515/hf-2019-0182; derived by the pack author",
       "license\tPublished-Scientific-Facts-Attribution-Retained", "source\tprimary\tproperties.tsv\tmaterial-tsv-v1"]
open(OUT + "/" + PACK + "/manifest.tsv", "w", encoding="utf-8").write("\n".join(man) + "\n")
print("wrote", OUT + "/" + PACK)
for k in VISC_ORDER:
    print("   %-36s %s Pa*s" % (k, fmt(visc[k])))
