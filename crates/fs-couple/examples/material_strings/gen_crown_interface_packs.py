"""Write the lateral stick joint of each Dominant string on the il Cannone bridge crown as a matdb interface source, from
contact mechanics, for the string rows of cannone-dominant-g3.material-strings (MagentaBadger 7145/7213, H1).

Every input is committed data:
- the bridge's maple shell pack (violin-maple-sycamore-sonderegger2013-shell);
- the string surfaces' metal packs;
- the CT contact table (break angles) and the CT dimensions (crown thickness);
- the performance's own string rows: each string's pack, radius and tension.

Contact. Crossed cylinders: the string (radius r_s, axis along the string) on the crown edge (radius r_b, axis across
the bridge top).
- Equivalent sphere R = sqrt(r_s r_b), normal load N = 2 T sin(theta/2), theta the break angle.
- Hertz radius a = (3 N R / (4 E*))^(1/3), 1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2 (Johnson 1985, eq. 4.53).
- Mindlin's no-slip tangential stiffness k_t = 8 G* a, 1/G* = (2-nu1)/G1 + (2-nu2)/G2 (Johnson 1985, eq. 7.42).

STATED UNKNOWNs, each enveloped at both ends. k_t grows with G and falls as E^(-1/3), so every compliant choice pairs
the smaller shear modulus with the larger Young modulus, and the stiff end the reverse:
- Maple end grain under a load along R: the orthotropic half-space lies between isotropic ones over the pack's moduli.
  The compliant end is G_RT with E_L, the stiff end G_LR with E_R.
- The crown edge radius runs from 0.1 mm to half the CT's minimum bridge-top thickness; k_t ~ r_b^(1/6).
- The wrap metal: Pickering's table names a wound steel E and three strings wound on a Perlon core, but not the
  winding metal. Wound strings are enveloped by the committed aluminium and silver moduli, and the E by aluminium and
  music wire.
- The winding's turn-by-turn contact is not resolved: the string is a smooth cylinder of its own radius.

The elastic envelope's mean pressures are far beyond what end-grain maple bears, so a played crown is bedded into a
conforming groove. Bedding only enlarges the contact: the compliant end is a lower bound, and the stiff end can be
exceeded.

Joint loss: the maple's Kelvin-Voigt retardation time in shear (pack viscosity / modulus) times the maple's share of the
tangential compliance, D = share tau k_t. The maple pack writes its viscosities at its constant-loss anchor omega_ref
(tau = eta / omega_ref), and each source states the same anchor as constant_loss_reference_omega: the joint's loss factor
is share x eta at every frequency, not the Kelvin-Voigt share x eta x omega / omega_ref. The metal's own loss is not in
its packs; its compliance share (printed) bounds what that omission contributes.

Applicability (MagentaBadger 7268 (a)): every row claims 0 to 2 pi 200 kHz. The contact is quasi-static there (the
maple's shear wavelength stays far above the contact radius), so the network may ask the law at a free coordinate and
at its coupled stiffness bound. The loss outside Laemmlein's measured 150 Hz-10 kHz (the maple pack's own band) is a
STATED UNKNOWN: there it is the network's five-branch constant-loss shape, extrapolated.

Each source carries both ends per area of that end's own contact: crown_lateral_{stiffness,viscosity}_per_area_{compliant,
stiff}. The string row's AREA is that end's pi a^2, printed here.

usage: python3 -I gen_crown_interface_packs.py REPO_ROOT PERFORMANCE OUT_ROOT
Reads the committed data under REPO_ROOT and the string rows of PERFORMANCE, and writes
OUT_ROOT/violin-cannone-crown-dominant-{e,a,d,g}-interface/{manifest.tsv,interface.tsv}.
"""
import math
import os
import sys

REPO, PERFORMANCE, OUT = sys.argv[1].rstrip("/"), sys.argv[2], sys.argv[3].rstrip("/")
SEED = os.path.join(REPO, "data/matdb/seed-v1")
GEOMETRY = os.path.join(REPO, "data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026")
MAPLE = "violin-maple-sycamore-sonderegger2013-shell"
ALUMINIUM = "aluminum-2024-t3-sheet-mil-hdbk-5j"
SILVER = ("silver-pure-smith-fickett-nist1995", "silver-pure-wrap-nist1995")
STEEL = "music-wire-spring-moduli-fuchs-1968"
PREFIXED = {"MPa": 1e6, "GPa": 1e9}


def scalars(pack):
    """property -> value of every scalar row of one seed source, with MPa and GPa moduli in Pa."""
    out = {}
    for line in open(os.path.join(SEED, pack, "properties.tsv"), encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "scalar":
            out[f[3]] = float(f[4]) * PREFIXED.get(f[5], 1.0)
    return out


maple = scalars(MAPLE)
tau = max(maple[v] / maple[g] for v, g in (("plane_stress_viscosity_66", "shear_modulus_12"),
                                            ("transverse_shear_viscosity_13", "shear_modulus_13"),
                                            ("transverse_shear_viscosity_23", "shear_modulus_23")))
nu_maple = maple["poisson_ratio_12"]
# (E, G) of the crown at the compliant and the stiff end.
crown = ((maple["young_modulus_1"], maple["shear_modulus_23"]), (maple["young_modulus_2"], maple["shear_modulus_12"]))
omega = [f for f in (line.rstrip("\n").split("\t") for line in open(os.path.join(SEED, MAPLE, "properties.tsv")))
         if f[0] == "validity" and f[2] == "omega"]
omega_lo, omega_hi = omega[0][3], omega[0][4]
anchor = maple["constant_loss_reference_omega"]
CLAIM = ("0", "%.9g" % (2 * math.pi * 200e3))
temperature = [f for f in (line.rstrip("\n").split("\t") for line in open(os.path.join(SEED, MAPLE, "properties.tsv")))
               if f[0] == "validity" and f[2] == "T"]
t_lo, t_hi = temperature[0][3], temperature[0][4]

al = scalars(ALUMINIUM)
ag, ag_shear = scalars(SILVER[0]), scalars(SILVER[1])
wire = scalars(STEEL)
# (name, E, nu, G) of each candidate surface metal.
aluminium = ("aluminium (" + ALUMINIUM + ")", al["young_modulus"], al["poisson_ratio"], al["shear_modulus"])
silver = ("silver (" + ", ".join(SILVER) + ")", ag["young_modulus"], ag["poisson_ratio"], ag_shear["shear_modulus"])
steel = ("steel (" + STEEL + ")", wire["young_modulus"], wire["young_modulus"] / (2 * wire["shear_modulus"]) - 1,
         wire["shear_modulus"])
SURFACES = {"e": (aluminium, steel), "a": (aluminium, silver), "d": (aluminium, silver), "g": (aluminium, silver)}

contacts = {}
for line in open(os.path.join(GEOMETRY, "string-contacts.tsv"), encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    if len(f) == 7 and f[1] == "bridge":
        contacts[f[0]] = math.atan(-float(f[5])) + math.atan(float(f[6]))
crown_m = None
for line in open(os.path.join(GEOMETRY, "dimensions.tsv"), encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    if len(f) > 4 and f[0] == "bridge" and f[2] == "thickness-z-min-above-y-minus-52" and f[4] == "m":
        crown_m = float(f[3])
radius = (1.0e-4, 0.5 * crown_m)

# The performance's string rows: pack, length, radius, tension.
rows = {}
for line in open(PERFORMANCE, encoding="utf-8"):
    f = line.split()
    if f and f[0] == "string" and f[4] == "tension":
        name = f[1].split("violin-string-dominant-")[1][0]
        rows[name] = (f[1].removesuffix(".fsmatpk"), float(f[3]), float(f[5]))
assert sorted(rows) == ["a", "d", "e", "g"], "the performance must string the four Dominant packs"


def ends(theta, r_s, tension, surfaces):
    """Both envelope ends: dicts of the contact's numbers."""
    load = 2 * tension * math.sin(theta / 2)
    out = []
    for r_b, (e_b, g_b), (metal, e_s, nu_s, g_s) in zip(radius, crown, surfaces):
        e_star = 1 / ((1 - nu_s**2) / e_s + (1 - nu_maple**2) / e_b)
        maple_shear, metal_shear = (2 - nu_maple) / g_b, (2 - nu_s) / g_s
        g_star = 1 / (maple_shear + metal_shear)
        a = (3 * load * math.sqrt(r_s * r_b) / (4 * e_star)) ** (1 / 3)
        k = 8 * g_star * a
        share = maple_shear / (maple_shear + metal_shear)
        area = math.pi * a * a
        out.append(dict(r_b=r_b, metal=metal, a=a, area=area, k=k, k_n=2 * e_star * a, pressure=load / area,
                        share=share, damping=share * tau * k))
    return load, out


areas = {}
print(f"maple: E_L {crown[0][0]:.4g} Pa, E_R {crown[1][0]:.4g} Pa, G_RT {crown[0][1]:.4g} Pa, G_LR {crown[1][1]:.4g} Pa, "
      f"nu {nu_maple}; shear retardation time {tau:.6g} s (loss factor {tau * anchor:.4g} at the anchor {anchor} rad/s); "
      f"maple validity T {t_lo}-{t_hi} K, omega {omega_lo}-{omega_hi} rad/s; joint claims omega {CLAIM[0]}-{CLAIM[1]} rad/s")
print(f"crown edge radius STATED UNKNOWN: {radius[0] * 1e3:.2f}-{radius[1] * 1e3:.2f} mm "
      f"(half the CT bridge-top thickness {crown_m * 1e3:.2f} mm)")
for name in ("e", "a", "d", "g"):
    pack, r_s, tension = rows[name]
    theta = contacts[name]
    load, both = ends(theta, r_s, tension, SURFACES[name])
    length = next(float(line.split()[2]) for line in open(PERFORMANCE) if line.startswith("string " + pack))
    own = tension / length
    print(f"{name}: {pack}, r {r_s * 1e3:.5g} mm, T {tension:.4f} N, break {math.degrees(theta):.2f} deg, "
          f"N {load:.4f} N; the string's own T/L {own:.0f} N/m")
    for label, e in zip(("compliant", "stiff"), both):
        print(f"  {label:9s} k_t {e['k']:.5g} N/m ({e['k'] / own:.0f}x T/L), k_n {e['k_n']:.4g} N/m, "
              f"a {e['a'] * 1e6:.1f} um, mean pressure {e['pressure'] / 1e6:.0f} MPa, {e['metal']}, "
              f"maple share {e['share']:.4f}, D {e['damping']:.5g} N s/m; string row AREA {e['area']:.9g} m^2")
    pack_id = f"violin-cannone-crown-dominant-{name}-interface"
    areas[pack] = (pack_id, [e["area"] for e in both])
    out = os.path.join(OUT, pack_id)
    os.makedirs(out, exist_ok=True)
    lines = [
        "frankensim.matdb-source.v1",
        f"surface_a\tsource-pack:{pack}\tdeclared solid model\tfrozen source-pack realization\t0\t"
        "outer-string-surface-at-the-bridge-crown",
        f"surface_b\tsource-pack:{MAPLE}\tdeclared solid model\tfrozen source-pack realization\t0\t"
        "bridge-crown-end-grain",
        "context\tdry-preloaded-contact\t-\tstring-tension-downbearing-at-the-cannone-break-angle\t"
        "stress-free-reference-about-the-static-preload",
        f"observation\tcrown-contact\tPickering 1985 Dominant {name.upper()} on the il Cannone maple bridge crown; "
        f"tension {tension:.4f} N (cannone-dominant-g3.material-strings), break angle {math.degrees(theta):.2f} deg "
        f"(CT string-contacts.tsv), normal load {load:.4f} N\tderived by gen_crown_interface_packs.py: Hertz crossed "
        "cylinders and Mindlin no-slip tangential stiffness (Johnson 1985 eqs. 4.53 and 7.42) from the cited packs; "
        "STATED UNKNOWN and enveloped at both ends: the maple orthotropy (G_RT with E_L, G_LR with E_R), the crown "
        "edge radius (0.1 mm to half the CT bridge-top thickness) and the unnamed wrap metal "
        f"({both[0]['metal']} to {both[1]['metal']}); joint loss the maple shear retardation time times its "
        f"compliance share, written at the maple's constant-loss anchor {anchor} rad/s (constant loss factor share x "
        "eta)\tthe elastic envelope's mean pressures exceed what end-grain maple bears, so a played "
        "crown is bedded: a conforming groove only enlarges the contact, the compliant end is a lower bound and "
        "the stiff end can be exceeded; the winding's turn-by-turn contact is not resolved and the metal's own "
        "loss is not in its packs; claimed 0 to 2 pi 200 kHz (quasi-static contact), the loss outside the maple's "
        f"measured {omega_lo}-{omega_hi} rad/s is a STATED UNKNOWN",
    ]
    for label, e in zip(("compliant", "stiff"), both):
        for kind, value, unit in (("stiffness", e["k"] / e["area"], "Pa/m"),
                                  ("viscosity", e["damping"] / e["area"], "Pa*s/m")):
            key = f"crown-{kind}-{label}"
            lines += [f"scalar\t{key}\tcrown-contact\tcrown_lateral_{kind}_per_area_{label}\t{value:.9g}\t{unit}\t"
                      "constant",
                      f"uncertainty\t{key}\tunstated\t-\t-\t-\t-",
                      f"validity\t{key}\tT\t{t_lo}\t{t_hi}\tK",
                      f"validity\t{key}\tomega\t{CLAIM[0]}\t{CLAIM[1]}\trad/s"]
    key = "crown-loss-anchor"
    lines += [f"scalar\t{key}\tcrown-contact\tconstant_loss_reference_omega\t{anchor:.9g}\trad/s\tconstant",
              f"uncertainty\t{key}\tunstated\t-\t-\t-\t-",
              f"validity\t{key}\tT\t{t_lo}\t{t_hi}\tK",
              f"validity\t{key}\tomega\t{CLAIM[0]}\t{CLAIM[1]}\trad/s"]
    open(os.path.join(out, "interface.tsv"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    metals = sorted({p for s in SURFACES[name] for p in s[0].split("(")[1].rstrip(")").split(", ")})
    open(os.path.join(out, "manifest.tsv"), "w", encoding="utf-8").write("\n".join([
        "frankensim.matdb-manifest.v1",
        f"pack_id\t{pack_id}",
        "redistribution\tpermitted\tFrankenSim derivation from committed packs and CT tables; no external rights",
        f"citation\tDerived from data/matdb/seed-v1 packs {MAPLE}, {', '.join(metals)} and the il Cannone CT tables "
        "string-contacts.tsv and dimensions.tsv; Johnson, Contact Mechanics (1985) eqs. 4.53 and 7.42",
        "license\tFrankenSim-Original-Work",
        "source\tprimary\tinterface.tsv\tinterface-tsv-v1"]) + "\n")

# The performance at each envelope end: every string row's scalar joint (the last two fields) becomes that end's interface.
# The stiff joints' coupled stiffness bound exceeds the network's 0.9 pi f_s guard at 144 kHz, so the stiff end renders
# at 384 kHz over the same duration (MagentaBadger 7268 (b)).
STIFF_RATE = 384000
stem = os.path.basename(PERFORMANCE).removesuffix(".material-strings")
for end, label in enumerate(("compliant", "stiff")):
    text = []
    for line in open(PERFORMANCE, encoding="utf-8").read().splitlines():
        f = line.split(" ")
        if f[0] == "audio" and label == "stiff":
            rate, samples = int(f[1]), int(f[2])
            if samples * STIFF_RATE % rate:
                sys.exit(f"the note's {samples} samples at {rate} Hz are not a whole number at {STIFF_RATE} Hz")
            f[1:3] = [str(STIFF_RATE), str(samples * STIFF_RATE // rate)]
        if f[0] == "string" and len(f) == 11 and f[1].removesuffix(".fsmatpk") in areas:
            pack_id, area = areas[f[1].removesuffix(".fsmatpk")]
            f[9:11] = ["interface", pack_id + ".fsintpk", f"{area[end]:.9g}",
                       f"crown_lateral_stiffness_per_area_{label}", f"crown_lateral_viscosity_per_area_{label}"]
        text.append(" ".join(f))
    open(os.path.join(OUT, f"{stem}-crown-{label}.material-strings"), "w", encoding="utf-8").write("\n".join(text) + "\n")
