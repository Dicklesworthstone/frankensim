"""Write the shell-ready sycamore maple pack of the il Cannone violin body (back, ribs, bridge) from the committed Sonderegger et al. (2013) and Laemmlein et al. (2020) packs.

Stiffness (one tree, 626 kg/m3 at about 12 percent moisture, 20 C / 65 percent RH): E_L = the static bending modulus, E_R = the mean of the static radial tension and compression moduli,
nu_LR = the mean Poisson ratio, G_LT and G_RT = the ULTRASONIC shear moduli (1 MHz: NOT audio band; no audio-band source was found, declared in the pack text). G_LR = the AUDIO-BAND in-plane
shear modulus of Laemmlein's 30 violin-back maple plates (maple-sycamore-lammlein2020-plates; plate modes near 0.7 kHz), replacing the ultrasonic value (orchestrator ruling, mail 6065).
Loss: eta_L and eta_LR are Laemmlein's violin-back maple plates. eta_R is the MEASURED radial loss of Alkadri et al. (2018)'s violin-back plate C (maple-sycamore-wavy-loss-alkadri2018,
loss_factor_radial_plate_c), the one plate of the twelve whose printed longitudinal loss holds Laemmlein's eta_L inside its printed sigma (orchestrator ruling, mail 13541; asserted below),
in place of the eta_LR surrogate. eta_LT and eta_RT are UNKNOWN for maple and are set equal to eta_LR as a declared surrogate (mail 5280 of the
orchestrator: T and RT/LT shear loss are unknown in every open source). The nine Kelvin-Voigt viscosities the shell owner consumes are V = D^(1/2) Q D^(1/2) / omega_ref with
omega_ref = 2 pi 200 rad/s, exactly as in the guitar packs (gen_wood_packs.py); they do not enter the frequencies. The pack states omega_ref as
constant_loss_reference_omega, so the shell reduction reads each loss factor at every mode's own frequency (fs_couple::shear_shell::constant_loss), and the
packs of the wood-loss lane replace these numbers.

usage: python3 -I gen_violin_maple_pack.py REPO_ROOT OUT_ROOT [shear=ultrasonic|shear=fpl] [inplane=laemmlein] [bridge=merhar2020]
Reads REPO_ROOT/data/matdb/seed-v1/{maple-sycamore-sonderegger2013,maple-sycamore-loss-lammlein2020,maple-sycamore-lammlein2020-plates,maple-sycamore-wavy-loss-alkadri2018}/properties.tsv and writes
OUT_ROOT/violin-maple-sycamore-sonderegger2013-shell/{manifest.tsv,properties.tsv}.
shear=ultrasonic writes the earlier pack (G_LR ultrasonic as well, the scalars of 4f9abefe1) under the name ...-shell-ultrasonicshear.
inplane=laemmlein (the audio-band pack, registered as violin-maple-sycamore-sonderegger2013-shell-laemmlein; the violin, viola and cello maple parts use it) takes the density, E_L, E_R and G_LR all from Laemmlein's plates (audio band, one consistent source) with G_LT, G_RT and
nu_LR as above; the pack is then named ...-shell-laemmlein.
With shear=fpl (a declared SENSITIVITY, not the registered model) G_LR and G_LT are the Wood Handbook sugar maple ratios times the same E_L (glr_over_el 0.111, glt_over_el 0.063 of the committed pack
maple-sugar-fpl-gtr282) instead of the ultrasonic values; G_RT stays ultrasonic (the Wood Handbook set has no G_RT ratio); the pack is then named ...-fplshear-shell.
bridge=merhar2020 (with inplane=laemmlein only) writes the violin BRIDGE's pack: the audio-band pack with its density replaced by the measured density of a maple violin bridge
(maple-bridge-aubert-merhar2020: Merhar and Humar 2020, Aubert bridge No. 1, from the bridge's own mass and volume) and every other row unchanged, so the stiffness and the loss
stay the back maple's, a declared borrowing (orchestrator ruling, mail 14630); the nine viscosities hold no density, so they are the parent pack's byte for byte. The pack is
then named ...-shell-laemmlein-bridge-merhar2020.
"""
import math
import os
import sys

REPO, OUT = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
FPL_SHEAR = "shear=fpl" in sys.argv[3:]
ULTRASONIC_SHEAR = "shear=ultrasonic" in sys.argv[3:]
INPLANE_LAEMMLEIN = "inplane=laemmlein" in sys.argv[3:]
BRIDGE_MERHAR = "bridge=merhar2020" in sys.argv[3:]
assert not BRIDGE_MERHAR or (INPLANE_LAEMMLEIN and not FPL_SHEAR and not ULTRASONIC_SHEAR), "bridge=merhar2020 is built on the audio-band pack only (inplane=laemmlein)"
SEED = REPO + "/data/matdb/seed-v1/"
OMEGA_MIN, OMEGA_MAX = 2.0 * math.pi * 150.0, 2.0 * math.pi * 10000.0  # the narrower of the sources' bands: Laemmlein's sweep, 150-10000 Hz
OMEGA_REF = 2.0 * math.pi * 200.0
PACK = ("violin-maple-sycamore-sonderegger2013-shell" + ("-fplshear" if FPL_SHEAR else "") + ("-ultrasonicshear" if ULTRASONIC_SHEAR else "")
        + ("-laemmlein" if INPLANE_LAEMMLEIN else "") + ("-bridge-merhar2020" if BRIDGE_MERHAR else ""))


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
DENSITY_OBS = "sond-shell"
if BRIDGE_MERHAR:
    # The violin bridge's one measured property: its density, from the mass and volume of Aubert bridge No. 1 (Merhar and Humar 2020, p.164 and Table 2,
    # unmodified). It replaces the back plates' 665.1 kg/m3 on the bridge only. The stiffness stays the back maple's, a declared borrowing: the paper's
    # elastic set (Table 2) is literature ratios scaled to one clamped-bridge mode, a fit, and is not used.
    mh = scalars("maple-bridge-aubert-merhar2020")
    stiffness_text = ("E_L %.0f MPa and E_R %.0f MPa of Laemmlein's violin-back plates (maple-sycamore-lammlein2020-plates, audio band); the plates' density %.1f kg/m3 is NOT "
                      "used: this bridge pack carries the measured density of a violin bridge (observation merhar-bridge)" % (E_L / 1e6, E_R / 1e6, density))
    density, DENSITY_OBS = mh["row:density-bridge-unmodified"], "merhar-bridge"
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
# eta_R: Alkadri's plate C, chosen by the pack's own eta_L. Alkadri measured tan delta_L and tan delta_R on the same twelve violin-back plates (Table 1);
# plate C (L 13.8 (2.2), R 24.0 (0.6) x 1e-3) is the only plate whose printed L mean +/- sigma holds Laemmlein's eta_L (every other plate's mean + sigma
# lies below it, at most 12.1e-3: plate J). Neither the range midpoint nor an R/L ratio carried onto eta_L is used.
ak = scalars("maple-sycamore-wavy-loss-alkadri2018")
ak_l, ak_l_sigma = ak["loss_factor_longitudinal_plate_c"], ak["loss_factor_longitudinal_plate_c_printed_sigma"]
assert abs(eta["L"] - ak_l) <= ak_l_sigma, "plate C's printed L loss no longer holds eta_L: the selection rule must be revisited"
assert ak["loss_factor_radial_printed_range_low"] <= ak["loss_factor_radial_plate_c"] <= ak["loss_factor_radial_printed_range_high"]
eta.update(R=ak["loss_factor_radial_plate_c"], LT=eta["LR"], RT=eta["LR"])
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
        "order in the committed Sonderegger pack). The validity band omega 2 pi 150 to 2 pi 10000 rad/s is the narrower of the sources' bands: Laemmlein's loss factors and audio-band G_LR "
        "come from a loudspeaker sweep over 150-10000 Hz (four plate modes of it, whose frequencies the source does not state), and the static and 1 MHz ultrasonic moduli bracket the "
        "audio band. Using a static or ultrasonic modulus inside the band is a STATED ASSUMPTION that the audio-band value lies between the two, with no dispersion applied; it is not a "
        "fit. T 280 to 320 K is the band in which this use is DECLARED by the pack author, not a measured band.")
OBS = [
    ("sond-shell", "Sycamore maple (Acer pseudoplatanus L.): ONE tree of Eastern Switzerland (Sonderegger, 626 kg/m3 at about 12 percent moisture, 20 C / 65 percent RH) and 30 violin-back plates "
     "(Laemmlein, 665 kg/m3, 20 C / 65 percent RH)",
     "Sonderegger et al. (2013) Eur J Wood Prod 71:91-99 as committed in the pack maple-sycamore-sonderegger2013 and Laemmlein et al. (2020) Holzforschung 74:765-776 as committed in "
     "maple-sycamore-lammlein2020-plates: %s; nu_LR = Sonderegger's mean over tension and compression %.2f; G_LR %.0f, G_LT %.0f and G_RT %.0f MPa: %s" % (
         stiffness_text, so["nu_lr"], g_lr, g_lt, g_rt, shear_text),
     "Ultrasonic moduli are not audio-band. Laemmlein's audio-band plates give G_LR 1580 MPa against Sonderegger's ultrasonic 1680, and E_L 12500 / E_R 1550 MPa against Sonderegger's static "
     "10900 / 1175 (inplane=laemmlein). Neither source is the il Cannone's wood. A violin-grade (curly, resonance) maple of the committed Spycher pack has E_L about 9700 to 10000 MPa and "
     "E_R about 1650 to 1700 MPa dynamically. " + AXES),
    ("laemmlein-maple-loss", "Sycamore maple violin-back wood: 30 quarter-sawn plates 95 x 60 x 2.6 mm, unvarnished (Laemmlein), and violin-back plate C of 12 mostly wavy sycamore plates, "
     "120 x 12 x 2 mm radial strips (Alkadri)",
     "Laemmlein et al. (2020) Holzforschung 74:765-776 as committed in the pack maple-sycamore-loss-lammlein2020: loss_factor_longitudinal %.4f and loss_factor_shear_lr %.4f; eta_L = %.4f, eta_LR = %.4f; "
     "Alkadri et al. (2018) IAWA J 39:63-86 as committed in the pack maple-sycamore-wavy-loss-alkadri2018: loss_factor_radial_plate_c %.4f (Table 1, tan delta_R of plate C, n 4, free-free "
     "bending at 80-500 Hz); eta_R = %.4f, taken from the one plate whose printed tan delta_L %.4f +/- %.4f holds Laemmlein's eta_L; eta_LT and eta_RT are UNKNOWN for maple and set equal to "
     "eta_LR" % (la["loss_factor_longitudinal"], la["loss_factor_shear_lr"], eta["L"], eta["LR"], ak["loss_factor_radial_plate_c"], eta["R"], ak_l, ak_l_sigma),
     "SURROGATE for two of the five directional loss factors (eta_LT and eta_RT). eta_R and eta_L come from different specimen sets matched by eta_L, not from one specimen; Alkadri's text "
     "(HAL PDF p.16) prints 0.022-0.026 radial, which disagrees with its Table 1, and the pack follows the table. The nine viscosities are Kelvin-Voigt values V = eta x stiffness / omega_ref "
     "with omega_ref = 2 pi 200 rad/s (1256.64 rad/s); the pack states omega_ref as constant_loss_reference_omega, and the shell reduction applies each loss factor at every mode's own "
     "frequency. They do not enter any frequency."),
]
if BRIDGE_MERHAR:
    OBS.append(("merhar-bridge", "Maple violin bridge No. 1 made by Aubert (one of three bridges: Nos. 1 and 2 by Aubert, No. 3 by Teller), as supplied, before the authors' thermal "
                "modification; the maple species is not named",
                "Merhar and Humar (2020) Drvna Industrija 71(2):163-169 as committed in the pack maple-bridge-aubert-merhar2020: density %.0f kg/m3 of the unmodified bridge, determined "
                "from measurements of the mass and volume of the bridge (p.164; Table 2)" % density,
                "ONE bridge, density only. The moisture content and the climate of the unmodified bridge's measurement are not printed (the modified pieces were stored a week at 23 C and "
                "65 percent RH). The bridge's stiffness and loss are not measured and stay the violin-back plates' values of this pack, a DECLARED BORROWING: the paper's elastic set is "
                "literature ratios (Kollmann 1975) adjusted so that the computed main mode of the clamped bridge matches its measured 4720 Hz, a fit, and its damping ratio 0.029 is the "
                "clamped bridge's modal bandwidth, not a material loss; neither is used."))


def fmt(v):
    return "%.9g" % v


os.makedirs(OUT + "/" + PACK, exist_ok=True)
lines = ["frankensim.matdb-source.v1"]
for o in OBS:
    lines.append("\t".join(("observation",) + o))
for k, u in BASE:
    for row in (("scalar", k, DENSITY_OBS if k == "density" else "sond-shell", k, fmt(vals[k]), u, "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", fmt(OMEGA_MIN), fmt(OMEGA_MAX), "rad/s")):
        lines.append("\t".join(row))
for k in VISC_ORDER:
    for row in (("scalar", k, "laemmlein-maple-loss", k, fmt(visc[k]), "Pa*s", "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", fmt(OMEGA_MIN), fmt(OMEGA_MAX), "rad/s")):
        lines.append("\t".join(row))
# The angular frequency at which the nine viscosities above were written: a reader applies each loss factor at every
# mode's own frequency with it (fs_couple::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA); without it they are Kelvin-Voigt.
k = "constant_loss_reference_omega"
for row in (("scalar", k, "laemmlein-maple-loss", k, fmt(OMEGA_REF), "rad/s", "constant"), ("uncertainty", k, "unstated", "-", "-", "-", "-"), ("validity", k, "T", "280", "320", "K"), ("validity", k, "omega", fmt(OMEGA_MIN), fmt(OMEGA_MAX), "rad/s")):
    lines.append("\t".join(row))
open(OUT + "/" + PACK + "/properties.tsv", "w", encoding="utf-8").write("\n".join(lines) + "\n")
man = ["frankensim.matdb-manifest.v1", "pack_id\t" + PACK,
       "redistribution\tpermitted\tNumerical measurement values retained as scientific facts with attribution; no article text, figure or scanned table is reproduced; derived by the pack author from the committed packs maple-sycamore-sonderegger2013, maple-sycamore-lammlein2020-plates, maple-sycamore-loss-lammlein2020 and maple-sycamore-wavy-loss-alkadri2018 (whose own terms apply to the Laemmlein values: publisher copyright on the article, ETH Research Collection copy In Copyright - Non-Commercial Use Permitted; and to the Alkadri value: HAL author file under the HAL deposit authorisation only, no Creative Commons licence, publisher copyright on the article)",
       "citation\tSonderegger W., Martienssen A., Nitsche C., Ozyhar T., Kaliske M., Niemz P. (2013), Eur J Wood Prod 71:91-99, doi:10.1007/s00107-012-0641-8; Laemmlein S.L., Van Damme B., Mannes D., Schwarze F.W.M.R., Burgert I. (2020), Holzforschung 74(8):765-776, doi:10.1515/hf-2019-0182; Alkadri A., Carlier C., Wahyudi I., Gril J., Langbour P., Bremaud I. (2018), IAWA Journal 39(1):63-86, doi:10.1163/22941932-20170185, Table 1; derived by the pack author",
       "license\tPublished-Scientific-Facts-Attribution-Retained", "source\tprimary\tproperties.tsv\tmaterial-tsv-v1"]
if BRIDGE_MERHAR:
    man[2] += ("; and the bridge density of the committed pack maple-bridge-aubert-merhar2020 (Merhar and Humar 2020, Creative Commons Attribution 4.0 (CC BY) open-access "
               "article)")
    man[3] = man[3].replace("; derived by the pack author", "; Merhar M., Humar M. (2020), Drvna Industrija 71(2):163-169, doi:10.5552/drvind.2020.1966, p.164 and Table 2; "
                            "derived by the pack author")
open(OUT + "/" + PACK + "/manifest.tsv", "w", encoding="utf-8").write("\n".join(man) + "\n")
print("wrote", OUT + "/" + PACK)
for k in VISC_ORDER:
    print("   %-36s %s Pa*s" % (k, fmt(visc[k])))
