"""Write the shell-ready sycamore maple pack of the il Cannone violin body (back, ribs, bridge) from the committed Sonderegger et al. (2013) and Laemmlein et al. (2020) packs.

Stiffness (one tree, 626 kg/m3 at about 12 percent moisture, 20 C / 65 percent RH): E_L = the static bending modulus, E_R = the mean of the static radial tension and compression moduli,
nu_LR = the mean Poisson ratio, G_LT and G_RT = the ULTRASONIC shear moduli (1 MHz: NOT audio band; no audio-band source was found, declared in the pack text). G_LR = the AUDIO-BAND in-plane
shear modulus of Laemmlein's 30 violin-back maple plates (maple-sycamore-lammlein2020-plates; plate modes near 0.7 kHz), replacing the ultrasonic value (orchestrator ruling, mail 6065).
Loss: eta_L and eta_LR are Laemmlein's violin-back maple plates; eta_R, eta_LT and eta_RT are UNKNOWN for maple and are set equal to eta_LR as a declared surrogate (mail 5280 of the
orchestrator: T and RT/LT shear loss are unknown in every open source). The nine Kelvin-Voigt viscosities the shell owner consumes are V = D^(1/2) Q D^(1/2) / omega_ref with
omega_ref = 2 pi 200 rad/s, exactly as in the guitar packs (gen_wood_packs.py); they do not enter the frequencies. A constant loss factor per mode is recovered by the example as
C_hyst(omega) = (omega_ref / omega) C_KV, and the packs of the wood-loss lane replace these numbers.

usage: python3 -I gen_violin_maple_pack.py REPO_ROOT OUT_ROOT [shear=ultrasonic|shear=fpl] [inplane=laemmlein]
Reads REPO_ROOT/data/matdb/seed-v1/{maple-sycamore-sonderegger2013,maple-sycamore-loss-lammlein2020,maple-sycamore-lammlein2020-plates}/properties.tsv and writes
OUT_ROOT/violin-maple-sycamore-sonderegger2013-shell/{manifest.tsv,properties.tsv}.
shear=ultrasonic writes the earlier pack (G_LR ultrasonic as well, the scalars of 4f9abefe1) under the name ...-shell-ultrasonicshear.
inplane=laemmlein (a declared ALTERNATIVE, not the registered model) takes the density, E_L, E_R and G_LR all from Laemmlein's plates (audio band, one consistent source) with G_LT, G_RT and
nu_LR as above; the pack is then named ...-shell-laemmlein.
With shear=fpl (a declared SENSITIVITY, not the registered model) G_LR and G_LT are the Wood Handbook sugar maple ratios times the same E_L (glr_over_el 0.111, glt_over_el 0.063 of the committed pack
maple-sugar-fpl-gtr282) instead of the ultrasonic values; G_RT stays ultrasonic (the Wood Handbook set has no G_RT ratio); the pack is then named ...-fplshear-shell.
"""
import math
import os
import sys

REPO, OUT = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
FPL_SHEAR = "shear=fpl" in sys.argv[3:]
ULTRASONIC_SHEAR = "shear=ultrasonic" in sys.argv[3:]
INPLANE_LAEMMLEIN = "inplane=laemmlein" in sys.argv[3:]
SEED = REPO + "/data/matdb/seed-v1/"
OMEGA_MAX = 12566.37
OMEGA_REF = 2.0 * math.pi * 200.0
PACK = ("violin-maple-sycamore-sonderegger2013-shell" + ("-fplshear" if FPL_SHEAR else "") + ("-ultrasonicshear" if ULTRASONIC_SHEAR else "")
        + ("-laemmlein" if INPLANE_LAEMMLEIN else ""))


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
lp = scalars("maple-sycamore-lammlein2020-plates")
E_L = so["young_modulus_longitudinal_bending"] * 1e6
E_R = 0.5 * (so["young_modulus_radial_tension"] + so["young_modulus_radial_compression"]) * 1e6
density = so["row:density-normal-climate"]
stiffness_text = ("E_L = the static bending modulus %.0f MPa; E_R = the mean of the static radial tension %.0f and compression %.0f MPa = %.0f MPa" % (
    so["young_modulus_longitudinal_bending"], so["young_modulus_radial_tension"], so["young_modulus_radial_compression"], E_R / 1e6))
if INPLANE_LAEMMLEIN:
    E_L, E_R, density = lp["young_modulus_longitudinal"] * 1e6, lp["young_modulus_radial"] * 1e6, lp["density"]
    stiffness_text = ("density %.1f kg/m3, E_L %.0f MPa and E_R %.0f MPa of Laemmlein's violin-back plates (maple-sycamore-lammlein2020-plates, audio band)" % (density, E_L / 1e6, E_R / 1e6))
g_lr, g_lt, g_rt = so["shear_modulus_lr_ultrasonic"], so["shear_modulus_lt_ultrasonic"], so["shear_modulus_rt_ultrasonic"]  # MPa
shear_text = "ULTRASONIC values (1 MHz transverse waves, 10 mm cubes)"
if FPL_SHEAR:
    fpl = scalars("maple-sugar-fpl-gtr282")
    g_lr, g_lt = fpl["glr_over_el"] * E_L / 1e6, fpl["glt_over_el"] * E_L / 1e6
    shear_text = "SENSITIVITY values: G_LR and G_LT from the Wood Handbook sugar maple ratios times E_L, G_RT ultrasonic"
    print("shear=fpl: G_LR = %.3f x E_L = %.0f MPa, G_LT = %.3f x E_L = %.0f MPa (sugar maple ratios, Wood Handbook GTR-282)" % (fpl["glr_over_el"], g_lr, fpl["glt_over_el"], g_lt))
elif not ULTRASONIC_SHEAR:
    g_lr = lp["shear_modulus_lr"]
    shear_text = ("G_LR from Laemmlein's violin-back plates (AUDIO band, identified from the (1,1) and (2,2) plate modes); G_LT and G_RT are ULTRASONIC values (1 MHz transverse "
                  "waves, 10 mm cubes; no audio-band source was found)")
vals = dict(density=density, young_modulus_1=E_L, young_modulus_2=E_R, poisson_ratio_12=so["nu_lr"], shear_modulus_12=g_lr * 1e6, shear_modulus_13=g_lt * 1e6,
            shear_modulus_23=g_rt * 1e6)
eta = dict(L=la["loss_factor_longitudinal"], LR=la["loss_factor_shear_lr"])
eta.update(R=eta["LR"], LT=eta["LR"], RT=eta["LR"])
print("density %.1f kg/m3, E_L %.4g Pa, E_R %.4g Pa, nu_LR %.3f, G_LR %.4g, G_LT %.4g, G_RT %.4g Pa" % (
    vals["density"], E_L, E_R, vals["poisson_ratio_12"], vals["shear_modulus_12"], vals["shear_modulus_13"], vals["shear_modulus_23"]))
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
    ("sond-shell", "Sycamore maple (Acer pseudoplatanus L.): ONE tree of Eastern Switzerland (Sonderegger, 626 kg/m3 at about 12 percent moisture, 20 C / 65 percent RH) and 30 violin-back plates "
     "(Laemmlein, 665 kg/m3, 20 C / 65 percent RH)",
     "Sonderegger et al. (2013) Eur J Wood Prod 71:91-99 as committed in the pack maple-sycamore-sonderegger2013 and Laemmlein et al. (2020) Holzforschung 74:765-776 as committed in "
     "maple-sycamore-lammlein2020-plates: %s; nu_LR = Sonderegger's mean over tension and compression %.2f; G_LR %.0f, G_LT %.0f and G_RT %.0f MPa: %s" % (
         stiffness_text, so["nu_lr"], g_lr, g_lt, g_rt, shear_text),
     "Ultrasonic moduli are not audio-band. Laemmlein's audio-band plates give G_LR 1580 MPa against Sonderegger's ultrasonic 1680, and E_L 12500 / E_R 1550 MPa against Sonderegger's static "
     "10900 / 1175 (inplane=laemmlein). Neither source is the il Cannone's wood. A violin-grade (curly, resonance) maple of the committed Spycher pack has E_L about 9700 to 10000 MPa and "
     "E_R about 1650 to 1700 MPa dynamically. " + AXES),
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
       "redistribution\tpermitted\tNumerical measurement values retained as scientific facts with attribution; no article text, figure or scanned table is reproduced; derived by the pack author from the committed packs maple-sycamore-sonderegger2013, maple-sycamore-lammlein2020-plates and maple-sycamore-loss-lammlein2020 (whose own terms apply to the Laemmlein values: publisher copyright on the article, ETH Research Collection copy In Copyright - Non-Commercial Use Permitted)",
       "citation\tSonderegger W., Martienssen A., Nitsche C., Ozyhar T., Kaliske M., Niemz P. (2013), Eur J Wood Prod 71:91-99, doi:10.1007/s00107-012-0641-8; Laemmlein S.L., Van Damme B., Mannes D., Schwarze F.W.M.R., Burgert I. (2020), Holzforschung 74(8):765-776, doi:10.1515/hf-2019-0182; derived by the pack author",
       "license\tPublished-Scientific-Facts-Attribution-Retained", "source\tprimary\tproperties.tsv\tmaterial-tsv-v1"]
open(OUT + "/" + PACK + "/manifest.tsv", "w", encoding="utf-8").write("\n".join(man) + "\n")
print("wrote", OUT + "/" + PACK)
for k in VISC_ORDER:
    print("   %-36s %s Pa*s" % (k, fmt(visc[k])))
