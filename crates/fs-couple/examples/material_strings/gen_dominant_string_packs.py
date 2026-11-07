"""Write the four Thomastik Dominant violin strings of Pickering (1985) as equivalent homogeneous material sections for
material_strings: one matdb source per string (E s01, A s12, D s22, G s32 of the committed seed
violin-strings-whole-pickering-jcas1985).

Per string, from the seed alone (Pickering's 329 mm test bed, the printed tension):
- radius r: half the mean of the printed maximum and minimum relaxed diameters;
- density = linear density / (pi r^2), so the section carries the measured mass per length;
- young_modulus = EI / (pi r^4 / 4) with EI = B T L^2 / pi^2. B is the least-squares inharmonicity of the printed partial
  deviations (partials 2-8, integer cents re n f1): 1200 log2(sqrt(1 + B n^2) / sqrt(1 + B));
- shear_modulus = GJ / (pi r^4 / 2) with GJ = k L / 4: the pendulum clamps the string's centre and both halves twist;
  k is the printed string contribution [g cm/rad].

So E and G are the whole wound string's bending and twisting stiffness, assigned to a cylinder of its own radius, and are
not the moduli of any one of its materials. A wound string is stiffer in twist than in bending, and G > E here.

Torsional loss: Pickering's damping coefficient is 100x the logarithmic decrement of the pendulum and the string swinging
TOGETHER, near 7 Hz. His method text (p. 6) gives the pendulum alone: 4.12 g, 5.4016 Hz, effective torsional stiffness
30.355 g cm/rad. The string holds k / (k + 30.355) of the stored energy, so its own loss factor is
eta_t = (eta_sys (k + k_p) - eta_p k_p) / k with eta_sys = decrement / pi. eta_p, the pendulum's own loss, is unstated.
It lies between 0 and the least, over every string of the table measured on that pendulum, of the value that gives the
string all of its system loss, eta_sys (k + k_p) / k_p. The pack writes the midpoint of that bracket and states both ends.

The printed k is used as printed. By the same method Pickering's solid-steel E strings (s02, s03, s05) read a shear
modulus 4 to 5 times steel's, so k may hold more than the string's own twisting stiffness; each section says so.

eta_t is written as a constant loss factor (fs_couple::bowed_string::runtime::reciprocal::bank::torsion): the Kelvin-Voigt
shear viscosity eta_t G / omega_ref with omega_ref = 2 pi 200 rad/s, and the anchor
shear_viscosity_constant_loss_reference_omega = omega_ref. Measured torsional Q of wound strings is nearly the same on
the first few modes (Wollman, Smith & Wolfe 2010, ISMA, Table 1; Woodhouse & Loach 1999, Acustica 85:734-740). A
dashpot would halve it from mode 1 to mode 2.

Bending loss is UNKNOWN: Pickering measures none. kelvin_voigt_bending_viscosity is written as 0, the lower end. The
string's transverse loss is then its air drag (from the performance's gas) and its supports only.

usage: python3 -I gen_dominant_string_packs.py REPO_ROOT OUT_ROOT
Reads REPO_ROOT/data/matdb/seed-v1/violin-strings-whole-pickering-jcas1985/properties.tsv and writes
OUT_ROOT/violin-string-dominant-{e,a,d,g}-pickering1985-section/{manifest.tsv,properties.tsv}.
"""
import math
import os
import sys

REPO, OUT = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
SEED = "violin-strings-whole-pickering-jcas1985"
G_CM = 9.80665e-5  # one gram-force centimetre [N m]
BED = 0.329  # Pickering's bridges [m]
PENDULUM_STIFFNESS = 30.355  # [g cm/rad], Pickering (1985) p. 6
PENDULUM_HZ = 5.4016  # Pickering (1985) p. 6
OMEGA_REF = 2.0 * math.pi * 200.0
OMEGA_MAX = 2.0 * math.pi * 40_000.0
STRINGS = [("e", "s01"), ("a", "s12"), ("d", "s22"), ("g", "s32")]
SOLID_STEEL = ("s02", "s03", "s05")

rows = {}
notes = {}
for line in open(REPO + "/data/matdb/seed-v1/" + SEED + "/properties.tsv", encoding="utf-8").read().splitlines():
    f = line.split("\t")
    if f[0] == "scalar":
        rows[f[3]] = float(f[4])
    elif f[0] == "observation":
        notes[f[1]] = f[2]


def inharmonicity(s):
    """Least-squares B over the printed partial deviations; golden-section search on [0, 2e-3]."""
    ns = [n for n in range(2, 9) if "partial_%d_deviation_cents_%s" % (n, s) in rows]
    cents = [rows["partial_%d_deviation_cents_%s" % (n, s)] for n in ns]

    def err(b):
        return sum((1200.0 * math.log2(math.sqrt(1.0 + b * n * n) / math.sqrt(1.0 + b)) - c) ** 2 for n, c in zip(ns, cents))

    lo, hi = 0.0, 2e-3
    g = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(200):
        a, b = hi - g * (hi - lo), lo + g * (hi - lo)
        if err(a) < err(b):
            hi = b
        else:
            lo = a
    b = 0.5 * (lo + hi)
    return b, math.sqrt(err(b) / len(ns)), ns


def system_loss(s):
    return rows["damping_coefficient_printed_x100_" + s] / 100.0 / math.pi


def pendulum_bound(s):
    """The pendulum's own loss factor that would give string s all of its system loss."""
    return system_loss(s) * (rows["torsion_stiffness_printed_g_cm_per_rad_" + s] + PENDULUM_STIFFNESS) / PENDULUM_STIFFNESS


TABLE = sorted(key[len("damping_coefficient_printed_x100_"):] for key in rows
               if key.startswith("damping_coefficient_printed_x100_")
               and "torsion_stiffness_printed_g_cm_per_rad_" + key[len("damping_coefficient_printed_x100_"):] in rows)
BOUNDING = min(TABLE, key=pendulum_bound)
pendulum_loss_max = pendulum_bound(BOUNDING)


def solid_steel_shear_ratio():
    """Shear modulus by this method on Pickering's solid-steel E strings, over steel's 80 GPa: (least, greatest)."""
    ratios = []
    for s in SOLID_STEEL:
        r = 0.25 * (rows["diameter_max_" + s] + rows["diameter_min_" + s])
        gj = rows["torsion_stiffness_printed_g_cm_per_rad_" + s] * G_CM * BED / 4.0
        ratios.append(gj / (math.pi * r ** 4 / 2.0) / 80e9)
    return min(ratios), max(ratios)


SOLID_LO, SOLID_HI = solid_steel_shear_ratio()


def fmt(v):
    return "%.9g" % v


for name, s in STRINGS:
    pack = "violin-string-dominant-%s-pickering1985-section" % name
    r = 0.25 * (rows["diameter_max_" + s] + rows["diameter_min_" + s])
    mu, tension = rows["linear_density_" + s], rows["tension_" + s]
    b, rms, ns = inharmonicity(s)
    ei = b * tension * BED ** 2 / math.pi ** 2
    k = rows["torsion_stiffness_printed_g_cm_per_rad_" + s]
    gj = k * G_CM * BED / 4.0
    density, young, shear = mu / (math.pi * r * r), ei / (math.pi * r ** 4 / 4.0), gj / (math.pi * r ** 4 / 2.0)
    eta_sys = system_loss(s)
    eta_hi = eta_sys * (k + PENDULUM_STIFFNESS) / k
    eta_lo = max(0.0, (eta_sys * (k + PENDULUM_STIFFNESS) - pendulum_loss_max * PENDULUM_STIFFNESS) / k)
    eta_t = 0.5 * (eta_lo + eta_hi)
    swing_hz = PENDULUM_HZ * math.sqrt(1.0 + k / PENDULUM_STIFFNESS)
    vals = [("density", density, "kg/m^3", "sec"), ("young_modulus", young, "Pa", "sec"), ("shear_modulus", shear, "Pa", "sec"),
            ("kelvin_voigt_bending_viscosity", 0.0, "Pa*s", "bend"),
            ("kelvin_voigt_shear_viscosity", eta_t * shear / OMEGA_REF, "Pa*s", "twist"),
            ("shear_viscosity_constant_loss_reference_omega", OMEGA_REF, "rad/s", "twist")]
    print("%s (%s): r %.5g m, density %.6g kg/m3, B %.4e (rms %.2f cents over partials %d-%d), E %.5g Pa, G %.5g Pa; "
          "system loss %.5f at %.3f Hz, string eta_t %.5f-%.5f, written %.5f" % (
              name.upper(), s, r, density, b, rms, ns[0], ns[-1], young, shear, eta_sys, swing_hz, eta_lo, eta_hi, eta_t))
    obs = [
        ("dominant-%s-section" % name, notes["pick85-" + s],
         "Equivalent homogeneous section from Pickering (1985) as committed in %s, row %s: radius %s m (mean relaxed "
         "diameter / 2); density = linear density %s kg/m / (pi r^2); young_modulus = B T L^2 / pi^2 / (pi r^4 / 4) with "
         "B %s, the least-squares inharmonicity of the printed deviations of partials %d-%d (rms residual %.2f cents), "
         "the printed tension T %s N and the test bed L %s m; shear_modulus = (k L / 4) / (pi r^4 / 2) with the printed "
         "string torsional stiffness k %s g cm/rad" % (SEED, s, fmt(r), fmt(mu), fmt(b), ns[0], ns[-1], rms, fmt(tension), fmt(BED), fmt(k)),
         "The stiffnesses are the whole wound string's, assigned to a cylinder of its own radius; they are not the moduli "
         "of any one of its materials. The printed torsional stiffness is used as printed, but by this same method "
         "Pickering's solid-steel E strings (%s) read a shear modulus %.1f to %.1f times steel's 80 GPa, so it may hold "
         "more than the string's own twisting stiffness. The deviations are printed to whole cents, which bounds the "
         "precision of B. " % (", ".join(SOLID_STEEL), SOLID_LO, SOLID_HI) +
         "Temperature and humidity are not stated by the source. The validity band omega 0 to 2 pi 40 kHz and T 280 to "
         "320 K is the band in which this use is DECLARED by the pack author, not a measured band."),
        ("dominant-%s-twist-loss" % name, "the same string",
         "Pickering (1985) printed damping coefficient %s = 100 x the logarithmic decrement of his torsion pendulum (p. 6: "
         "4.12 g, %s Hz, effective torsional stiffness %s g cm/rad) clamped at the string's centre, swinging at %.3f Hz: "
         "system loss factor %.6f; the string's own torsional loss factor %.6f to %.6f (the pendulum's own loss from 0 to "
         "%.6f, the least over the table's %d strings of the value that gives a string all of its system loss, reached by "
         "%s); written %.6f, the midpoint" % (
             fmt(rows["damping_coefficient_printed_x100_" + s]), fmt(PENDULUM_HZ), fmt(PENDULUM_STIFFNESS), swing_hz, eta_sys,
             eta_lo, eta_hi, pendulum_loss_max, len(TABLE), BOUNDING, eta_t),
         "Measured near 7 Hz and applied as a CONSTANT loss factor at every twist mode's own frequency: "
         "kelvin_voigt_shear_viscosity = eta_t G / omega_ref with shear_viscosity_constant_loss_reference_omega = "
         "omega_ref = 2 pi 200 rad/s, an encoding frequency only. Wound strings' measured torsional Q is nearly the same on "
         "the first few modes (Wollman, Smith & Wolfe 2010, ISMA, Table 1; Woodhouse & Loach 1999, Acustica 85:734-740)."),
        ("dominant-%s-bending-loss" % name, "the same string",
         "No measurement: Pickering (1985) reports no transverse damping",
         "UNKNOWN. kelvin_voigt_bending_viscosity is written as 0, the lower end, so the string's transverse loss is its air "
         "drag and its supports only. Measured bending loss factors of other wound strings: 0.047 to 0.13 (D'Addario Kaplan "
         "cello strings; Mansour, Woodhouse & Scavone 2016, Acta Acust. 102:1082-1093, Table 1) and 0.002 (a steel violin "
         "E; Galluzzo, Woodhouse & Mansour 2017, Acta Acust. 103:1080-1099, Table 1)."),
    ]
    lines = ["frankensim.matdb-source.v1"]
    for o in obs:
        lines.append("\t".join(("observation",) + o))
    for key, value, unit, part in vals:
        observation = {"sec": obs[0][0], "twist": obs[1][0], "bend": obs[2][0]}[part]
        for row in (("scalar", key, observation, key, fmt(value), unit, "constant"), ("uncertainty", key, "unstated", "-", "-", "-", "-"),
                    ("validity", key, "T", "280", "320", "K"), ("validity", key, "omega", "0", fmt(OMEGA_MAX), "rad/s")):
            lines.append("\t".join(row))
    os.makedirs(OUT + "/" + pack, exist_ok=True)
    open(OUT + "/" + pack + "/properties.tsv", "w", encoding="utf-8").write("\n".join(lines) + "\n")
    man = ["frankensim.matdb-manifest.v1", "pack_id\t" + pack,
           "redistribution\tpermitted\tNumerical measurement values retained as scientific facts with attribution; no article text, figure "
           "or scanned table is reproduced; derived by the pack author from the committed pack " + SEED + " (whose own terms apply) "
           "and the pendulum constants of the same source's method text",
           "citation\tPickering N.C. (1985), Physical Properties of Violin Strings, Journal of the Catgut Acoustical Society 44:6-8 "
           "(November 1985), the table on p. 7 and the method text on p. 6; derived by the pack author",
           "license\tPublished-Scientific-Facts-Attribution-Retained", "source\tprimary\tproperties.tsv\tmaterial-tsv-v1"]
    open(OUT + "/" + pack + "/manifest.tsv", "w", encoding="utf-8").write("\n".join(man) + "\n")
print("pendulum's own loss factor bound %.6f" % pendulum_loss_max)
