"""Write Pickering's (1985) Super Sensitive E string, solid steel, medium gauge (row s02 of the committed seed
violin-strings-whole-pickering-jcas1985), as an equivalent homogeneous material section for material_strings, with the
three thermal values Heptonstall, Waller & Robertson print for high-carbon music wire (LIGO P1400229, Table 3), so the
string can select `bending-loss thermoelastic`.

Assumption, named: Pickering's solid-steel E is high-carbon music wire. The section therefore takes:
- from Pickering's row s02: radius r = half the mean printed relaxed diameter, density = linear density / (pi r^2);
- music wire's moduli, NOT the stiffnesses Pickering's inharmonicity and pendulum imply (by the method of
  gen_dominant_string_packs.py they read E 2.3 and G 5.3 times the values written here on this string): young_modulus 212 GPa
  (Heptonstall, Waller & Robertson, Section III.A) and shear_modulus 11.5 Msi = 79.29 GPa (Fuchs 1968, Table 2,
  quoting the SAE Manual on Helical Springs; a secondary spring-industry value);
- Heptonstall's Table 3 expansion, heat capacity and conductivity (the authors' inputs, not measured on their
  specimens), at 293.15 K only, keeping that validity.
So the inharmonicity follows from steel's E, not from Pickering's printed cents.

Twist loss: Pickering's s02 decrement by the method of gen_dominant_string_packs.py, the string's share of the swing
energy from his printed k (the measured pendulum frequency), the pendulum's own loss bounded by every string of his
table; the midpoint is written as a constant loss factor on the music-wire G.

usage: python3 -I gen_super_sensitive_e_pack.py REPO_ROOT OUT_ROOT
"""
import math
import os
import sys

REPO, OUT = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
SEED = "violin-strings-whole-pickering-jcas1985"
THERMAL = "steel-music-wire-string-heptonstall-ligo-p1400229"
SHEAR = "music-wire-spring-moduli-fuchs-1968"
PACK = "violin-string-super-sensitive-e-pickering1985-music-wire-section"
ROW = "s02"
G_CM = 9.80665e-5
BED = 0.329
PENDULUM_STIFFNESS = 30.355
PENDULUM_HZ = 5.4016
OMEGA_REF = 2.0 * math.pi * 200.0
OMEGA_MAX = 2.0 * math.pi * 40_000.0


def read(pack):
    rows, notes, keys = {}, {}, {}
    for line in open(REPO + "/data/matdb/seed-v1/" + pack + "/properties.tsv", encoding="utf-8").read().splitlines():
        f = line.split("\t")
        if f[0] == "scalar":
            rows[f[3]] = float(f[4])
            keys[f[3]] = f
        elif f[0] == "observation":
            notes[f[1]] = f[2:]
    return rows, notes, keys


rows, notes, _ = read(SEED)
hrows, hnotes, hkeys = read(THERMAL)
frows, fnotes, fkeys = read(SHEAR)
assert hkeys["young_modulus"][5] == "MPa" and fkeys["shear_modulus"][5] == "MPa", (hkeys["young_modulus"], fkeys["shear_modulus"])
YOUNG = hrows["young_modulus"] * 1e6
SHEAR_MODULUS = frows["shear_modulus"] * 1e6


def inharmonicity(s):
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


def fmt(v):
    return "%.9g" % v


table = sorted(k[len("damping_coefficient_printed_x100_"):] for k in rows if k.startswith("damping_coefficient_printed_x100_"))
bounds = {s: system_loss(s) * (rows["torsion_stiffness_printed_g_cm_per_rad_" + s] + PENDULUM_STIFFNESS) / PENDULUM_STIFFNESS
          for s in table if "torsion_stiffness_printed_g_cm_per_rad_" + s in rows}
bounding = min(bounds, key=bounds.get)
pendulum_loss_max = bounds[bounding]

s = ROW
r = 0.25 * (rows["diameter_max_" + s] + rows["diameter_min_" + s])
mu, tension = rows["linear_density_" + s], rows["tension_" + s]
b, rms, ns = inharmonicity(s)
ei = b * tension * BED ** 2 / math.pi ** 2
k = rows["torsion_stiffness_printed_g_cm_per_rad_" + s]
gj = k * G_CM * BED / 4.0
density, young_p, shear_p = mu / (math.pi * r * r), ei / (math.pi * r ** 4 / 4.0), gj / (math.pi * r ** 4 / 2.0)
young, shear = YOUNG, SHEAR_MODULUS
b_steel = math.pi ** 2 * young * (math.pi * r ** 4 / 4.0) / (tension * BED ** 2)
eta_sys = system_loss(s)
eta_hi = eta_sys * (k + PENDULUM_STIFFNESS) / k
eta_lo = max(0.0, (eta_sys * (k + PENDULUM_STIFFNESS) - pendulum_loss_max * PENDULUM_STIFFNESS) / k)
eta_t = 0.5 * (eta_lo + eta_hi)
swing_hz = PENDULUM_HZ * math.sqrt(1.0 + k / PENDULUM_STIFFNESS)
print("%s: r %.6g m, density %.6g kg/m3; written E %.6g Pa, G %.6g Pa (Pickering's method gives E %.6g = %.2fx, G %.6g "
      "= %.2fx; his cents fit B %.4e, steel's E gives B %.4e on his bed); system loss %.6f at %.3f Hz, string eta_t "
      "%.6f-%.6f, written %.6f; pendulum bound %.6f from %s (of %d rows)" % (
          s, r, density, young, shear, young_p, young_p / young, shear_p, shear_p / shear, b, b_steel, eta_sys, swing_hz,
          eta_lo, eta_hi, eta_t, pendulum_loss_max, bounding, len(bounds)))

thermal_keys = [("thermal-expansion", "linear_thermal_expansion_coefficient", "K^-1"),
                ("specific-heat", "specific_heat_capacity", "J/kg/K"),
                ("thermal-conductivity", "thermal_conductivity", "W/m/K")]
for _, key, unit in thermal_keys:
    f = hkeys[key]
    assert f[5] == unit and f[6] == "constant", f
assert hnotes["hl-table3-thermal"][0].startswith("Music wire thermal parameters used by the authors"), hnotes["hl-table3-thermal"]
thermal_text = hnotes["hl-table3-thermal"]

obs = [
    ("super-sensitive-e-section", notes["pick85-" + s][0],
     "Solid circular section of Pickering (1985) as committed in %s, row %s: radius %s m (mean relaxed diameter / 2); "
     "density = linear density %s kg/m / (pi r^2). The printed tension %s N tuned it to E on his %s m bed." % (
         SEED, s, fmt(r), fmt(mu), fmt(tension), fmt(BED)),
     "A solid steel string (the table's note 2). Temperature and humidity are not stated by the source. The validity "
     "band omega 0 to 2 pi 40 kHz and T 280 to 320 K is the band in which this use is DECLARED by the pack author, not a "
     "measured band."),
    ("music-wire-moduli", "High-carbon steel music wire",
     "young_modulus %s Pa from the committed pack %s (Heptonstall, Waller & Robertson, Section III.A, tensile tests on "
     "untreated music wire); shear_modulus %s Pa from the committed pack %s (Fuchs 1968, Table 2: music wire, cold wound, "
     "G = 11.5 Msi, quoting the SAE Manual on Helical Springs)" % (fmt(young), THERMAL, fmt(shear), SHEAR),
     "ASSUMPTION named by the pack author: Pickering's solid-steel E is high-carbon music wire. Pickering's own "
     "measurements do not give these moduli: by the method of gen_dominant_string_packs.py his partials (B %s, rms %.2f "
     "cents over %d-%d) and torsional stiffness (k %s g cm/rad) read E %s Pa and G %s Pa on this string, %.1f and %.1f "
     "times these; steel's E gives B %s on his bed. So this string's inharmonicity follows from steel's E, not from his "
     "printed cents. The G value is a secondary spring-industry figure. Both are written over the declared use band of "
     "this pack (T 280 to 320 K), not their sources' single temperature." % (
         fmt(b), rms, ns[0], ns[-1], fmt(k), fmt(young_p), fmt(shear_p), young_p / young, shear_p / shear, fmt(b_steel))),
    ("super-sensitive-e-twist-loss", "the same string",
     "Pickering (1985) printed damping coefficient %s = 100 x the logarithmic decrement of his torsion pendulum (p. 6: "
     "4.12 g, %s Hz, effective torsional stiffness %s g cm/rad) clamped at the string's centre, swinging at %.3f Hz: "
     "system loss factor %.6f; the string's own torsional loss factor %.6f to %.6f (the pendulum's own loss from 0 to "
     "%.6f, the least over the table's %d strings of the value that gives a string all of its system loss, reached by "
     "%s); written %.6f, the midpoint" % (
         fmt(rows["damping_coefficient_printed_x100_" + s]), fmt(PENDULUM_HZ), fmt(PENDULUM_STIFFNESS), swing_hz, eta_sys,
         eta_lo, eta_hi, pendulum_loss_max, len(bounds), bounding, eta_t),
     "The string's share of the swing energy uses Pickering's printed k, which his measured swing frequency gives; that "
     "k reads more than music wire's G implies (see music-wire-moduli), so the share, and with it this loss factor, "
     "belong to his clamped specimen. Measured near 7 Hz and applied as a CONSTANT loss factor at every twist mode's own "
     "frequency, as in the Dominant packs: kelvin_voigt_shear_viscosity = eta_t G / omega_ref with the music-wire G and "
     "shear_viscosity_constant_loss_reference_omega = omega_ref = 2 pi 200 rad/s, an encoding frequency only. The "
     "constant-Q form is measured for wound strings (Wollman, Smith & Wolfe 2010; Woodhouse & Loach 1999); for a solid "
     "string it is an assumption."),
    ("super-sensitive-e-bending-loss", "the same string",
     "No measurement: Pickering (1985) reports no transverse damping",
     "kelvin_voigt_bending_viscosity is written as 0. The thermal values below let the card select the circular "
     "thermoelastic bending loss of a homogeneous solid section; any further anelastic loss stays UNKNOWN (a steel violin "
     "E's measured bending loss factor 0.002: Galluzzo, Woodhouse & Mansour 2017, Acta Acust. 103:1080-1099, Table 1)."),
    ("music-wire-thermal", thermal_text[0],
     thermal_text[1] + " Copied from the committed pack " + THERMAL + " (observation hl-table3-thermal).",
     "ASSUMPTION named by the pack author: Pickering's solid-steel E is high-carbon music wire. Pickering gives no steel "
     "grade or supplier specification for the string. The values are held at 293.15 K, the only temperature they are "
     "given for; a card held at another temperature refuses them."),
]
vals = [("density", density, "kg/m^3", 0), ("young_modulus", young, "Pa", 1), ("shear_modulus", shear, "Pa", 1),
        ("kelvin_voigt_bending_viscosity", 0.0, "Pa*s", 3),
        ("kelvin_voigt_shear_viscosity", eta_t * shear / OMEGA_REF, "Pa*s", 2),
        ("shear_viscosity_constant_loss_reference_omega", OMEGA_REF, "rad/s", 2)]
assert [o[0] for o in obs] == ["super-sensitive-e-section", "music-wire-moduli", "super-sensitive-e-twist-loss",
                               "super-sensitive-e-bending-loss", "music-wire-thermal"], [o[0] for o in obs]
lines = ["frankensim.matdb-source.v1"]
for o in obs:
    lines.append("\t".join(("observation",) + o))
for key, value, unit, part in vals:
    for row in (("scalar", key, obs[part][0], key, fmt(value), unit, "constant"), ("uncertainty", key, "unstated", "-", "-", "-", "-"),
                ("validity", key, "T", "280", "320", "K"), ("validity", key, "omega", "0", fmt(OMEGA_MAX), "rad/s")):
        lines.append("\t".join(row))
for _, key, unit in thermal_keys:
    for row in (("scalar", key, obs[4][0], key, hkeys[key][4], unit, "constant"), ("uncertainty", key, "unstated", "-", "-", "-", "-"),
                ("validity", key, "T", "293.15", "293.15", "K")):
        lines.append("\t".join(row))
os.makedirs(OUT + "/" + PACK, exist_ok=True)
open(OUT + "/" + PACK + "/properties.tsv", "w", encoding="utf-8").write("\n".join(lines) + "\n")
# Manifest v3: thermoelastic bending needs exactly one semantically typed absolute-temperature axis, so T is declared as
# one (as in the Heptonstall pack); omega keeps the untyped meaning it has in a v1 manifest.
man = ["frankensim.matdb-manifest.v3", "pack_id\t" + PACK,
       "redistribution\tpermitted\tNumerical measurement values retained as scientific facts with attribution; no article text, figure "
       "or scanned table is reproduced; derived by the pack author from the committed packs " + SEED + ", " + THERMAL +
       " and " + SHEAR + " (whose own terms apply) and the pendulum constants of Pickering's method text",
       "citation\tPickering N.C. (1985), Physical Properties of Violin Strings, Journal of the Catgut Acoustical Society 44:6-8 "
       "(November 1985), the table on p. 7 and the method text on p. 6; Heptonstall A., Waller M., Robertson N.A., Investigation "
       "of mechanical properties of cryogenically treated music wire, LIGO document P1400229 / arXiv:1506.03778, Section III.A "
       "and Table 3; Fuchs H.O. (1968), Minimum-Weight Springs, 3rd Aerospace Mechanisms Symposium (JPL), NASA, Table 2 "
       "(quoting the SAE Manual on Helical Springs); derived by the pack author",
       "license\tPublished-Scientific-Facts-Attribution-Retained", "source\tprimary\tproperties.tsv\tmaterial-tsv-v1",
       "axis\tT\tT\tabsolute-temperature", "axis\tomega\tomega\tdimensional"]
open(OUT + "/" + PACK + "/manifest.tsv", "w", encoding="utf-8").write("\n".join(man) + "\n")
