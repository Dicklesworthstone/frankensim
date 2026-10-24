"""Write the shell-ready Norway spruce pack of the il Cannone violin body (top, bass bar, soundpost) from the committed Viala et al. (2020) violin-quarter pack.

Stiffness and density: the averages over 10 Norway spruce (Picea abies) violin quarters identified at audio frequency by finite element model updating
(spruce-norway-viala2020-quarters): E_L, E_R, G_LR, G_LT (the source's G_TL) and G_RT. nu_LR is not identified by the source; the Sitka value of the
committed Wood Handbook pack (spruce-sitka-fpl-gtr282, nu-lr) is a declared SURROGATE. Loss: eta_L, eta_R and eta_LR of the same quarters; eta_LT and
eta_RT are UNKNOWN and set equal to eta_LR, as in the maple pack. The nine Kelvin-Voigt viscosities are V = eta x stiffness / omega_ref with
omega_ref = 2 pi 200 rad/s, exactly as in gen_violin_maple_pack.py; they do not enter the frequencies. The pack states omega_ref as
constant_loss_reference_omega, so the shell reduction reads each loss factor at every mode's own frequency (fs_couple::shear_shell::constant_loss).

usage: python3 -I gen_violin_spruce_pack.py REPO_ROOT OUT_ROOT
Reads REPO_ROOT/data/matdb/seed-v1/{spruce-norway-viala2020-quarters,spruce-sitka-fpl-gtr282}/properties.tsv and writes OUT_ROOT/violin-spruce-norway-viala2020-shell/{manifest.tsv,properties.tsv}.
"""
import math
import os
import sys

REPO, OUT = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
SEED = REPO + "/data/matdb/seed-v1/"
OMEGA_MAX = 12566.37
OMEGA_REF = 2.0 * math.pi * 200.0
PACK = "violin-spruce-norway-viala2020-shell"


def scalars(pack):
    out = {}
    for line in open(SEED + pack + "/properties.tsv", encoding="utf-8").read().splitlines():
        f = line.split("\t")
        if f[0] == "scalar":
            out[f[3]] = float(f[4])
    return out


vi = scalars("spruce-norway-viala2020-quarters")
fpl = scalars("spruce-sitka-fpl-gtr282")
E_L = vi["young_modulus_longitudinal"] * 1e6
E_R = vi["young_modulus_radial"] * 1e6
vals = dict(density=vi["density"], young_modulus_1=E_L, young_modulus_2=E_R, poisson_ratio_12=fpl["nu_lr"], shear_modulus_12=vi["shear_modulus_lr"] * 1e6,
            shear_modulus_13=vi["shear_modulus_tl"] * 1e6, shear_modulus_23=vi["shear_modulus_rt"] * 1e6)
eta = dict(L=vi["loss_factor_longitudinal"], R=vi["loss_factor_radial"], LR=vi["loss_factor_shear_lr"])
eta.update(LT=eta["LR"], RT=eta["LR"])
print("E_L %.4g, E_R %.4g Pa, nu_LR %.3f (Sitka surrogate), G_LR %.4g, G_LT %.4g, G_RT %.4g Pa, density %.0f kg/m3" % (
    E_L, E_R, vals["poisson_ratio_12"], vals["shear_modulus_12"], vals["shear_modulus_13"], vals["shear_modulus_23"], vals["density"]))
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
        "shear_modulus_13 = G_LT and shear_modulus_23 = G_RT, poisson_ratio_12 = nu_LR (strain along R per strain along L under stress along L). The validity band omega 0 to 2 pi 2000 rad/s "
        "and T 280 to 320 K is the band in which this use is DECLARED by the pack author, not a measured band.")
OBS = [
    ("viala-shell", "Norway spruce (Picea abies) violin-top tonewood, 10 split quarters from one seller, 9-11 percent moisture",
     "Viala et al. (2020) J. Cult. Herit. 42:108-116 as committed in the pack spruce-norway-viala2020-quarters (audio-band finite element model updating, 100-5000 Hz): density %.0f kg/m3, "
     "E_L %.0f MPa, E_R %.0f MPa, G_LR %.0f MPa, G_LT (the source's G_TL) %.0f MPa, G_RT %.0f MPa; nu_LR %.3f from the committed pack spruce-sitka-fpl-gtr282" % (
         vi["density"], vi["young_modulus_longitudinal"], vi["young_modulus_radial"], vi["shear_modulus_lr"], vi["shear_modulus_tl"], vi["shear_modulus_rt"], fpl["nu_lr"]),
     "nu_LR is a SURROGATE (Sitka, Wood Handbook): the source does not identify it. The quarters are violin tonewood, not the il Cannone top, whose species the ESRF record does not state. " + AXES),
    ("viala-shell-loss", "the same quarters",
     "Viala et al. (2020) Table 4 average row as committed in spruce-norway-viala2020-quarters: eta_L %.4f, eta_R %.4f, eta_LR %.4f; eta_LT and eta_RT are UNKNOWN and set equal to eta_LR" % (
         eta["L"], eta["R"], eta["LR"]),
     "SURROGATE for two of the five directional loss factors. The nine viscosities are Kelvin-Voigt values V = eta x stiffness / omega_ref with omega_ref = 2 pi 200 rad/s (1256.64 rad/s); the pack "
     "states omega_ref as constant_loss_reference_omega, and the shell reduction applies each loss factor at every mode's own frequency. They do not enter any frequency."),
]


def fmt(v):
    return "%.9g" % v


os.makedirs(OUT + "/" + PACK, exist_ok=True)
lines = ["frankensim.matdb-source.v1"]
for o in OBS:
    lines.append("\t".join(("observation",) + o))
for k, u in BASE:
    for row in (("scalar", k, "viala-shell", k, fmt(vals[k]), u, "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", "0", fmt(OMEGA_MAX), "rad/s")):
        lines.append("\t".join(row))
for k in VISC_ORDER:
    for row in (("scalar", k, "viala-shell-loss", k, fmt(visc[k]), "Pa*s", "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", "0", fmt(OMEGA_MAX), "rad/s")):
        lines.append("\t".join(row))
# The angular frequency at which the nine viscosities above were written: a reader applies each loss factor at every
# mode's own frequency with it (fs_couple::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA); without it they are Kelvin-Voigt.
k = "constant_loss_reference_omega"
for row in (("scalar", k, "viala-shell-loss", k, fmt(OMEGA_REF), "rad/s", "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", "0", fmt(OMEGA_MAX), "rad/s")):
    lines.append("\t".join(row))
open(OUT + "/" + PACK + "/properties.tsv", "w", encoding="utf-8").write("\n".join(lines) + "\n")
man = ["frankensim.matdb-manifest.v1", "pack_id\t" + PACK,
       "redistribution\tpermitted\tNumerical measurement values retained as scientific facts with attribution; no article text, figure or scanned table is reproduced; derived by the pack author from the committed packs spruce-norway-viala2020-quarters and spruce-sitka-fpl-gtr282 (whose own terms apply)",
       "citation\tViala R., Placet V., Cogan S. (2020), J. Cult. Herit. 42:108-116, doi:10.1016/j.culher.2019.09.004; USDA Forest Products Laboratory (2021), Wood Handbook, FPL-GTR-282, Table 5-2 (nu_LR); derived by the pack author",
       "license\tPublished-Scientific-Facts-Attribution-Retained", "source\tprimary\tproperties.tsv\tmaterial-tsv-v1"]
open(OUT + "/" + PACK + "/manifest.tsv", "w", encoding="utf-8").write("\n".join(man) + "\n")
print("wrote", OUT + "/" + PACK)
