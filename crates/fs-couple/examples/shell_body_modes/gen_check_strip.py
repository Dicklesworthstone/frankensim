"""Tables for an independent check of the shell and tie machinery of shell_body_modes against closed-form beam frequencies (no specimen data involved).

usage: python3 -I gen_check_strip.py REPO_ROOT KIND OUT.tsv NX [key=value ...]
  KIND = plate : one free orthotropic strip (grain along its length), L x b x h
         split : the same strip as two components tied end to end at the midsurface (rigid, full rotations)
         tbeam : the strip with a thin wall (rib) standing on its upper skin along the whole length, tied on the skin (side +1), grain along the length: a T section
  NX = number of elements along the length (width has 2 elements; the rib has K layers over its height).
Options: L_mm (500), b_mm (25), h_mm (3), H_mm (12 rib height), w_mm (3 rib thickness), K (4 element layers over the rib height), pack (spruce).

The material is the guitar spruce shell pack (REPO_ROOT/data/matdb/seed-v1/guitar-spruce-sitka-fpl-gtr282-shell): density, E_L, E_R, nu_LR as the pack states them. Prints the closed-form free-free
Euler-Bernoulli frequencies of the section (plate: E_L I / (rho A); tbeam: the T section by the parallel-axis theorem, flange and web both E_L), the plane-strain-like correction of the orthotropic
flange E_L / (1 - nu_LR nu_RL) as a bracket, and a Timoshenko shear estimate; the table is written in the format frankensim.shell-assembly.v1 that the example reads. Standard library only.
"""
import math
import sys

repo, kind, out, nx = sys.argv[1].rstrip("/"), sys.argv[2], sys.argv[3], int(sys.argv[4])
opt = dict(a.split("=", 1) for a in sys.argv[5:])
L = float(opt.get("L_mm", "500")) / 1000.0
b = float(opt.get("b_mm", "25")) / 1000.0
h = float(opt.get("h_mm", "3")) / 1000.0
H = float(opt.get("H_mm", "12")) / 1000.0
w = float(opt.get("w_mm", "3")) / 1000.0
pack = opt.get("pack", "spruce")
ny = 2
K = int(opt.get("K", "4"))  # element layers over the rib height

props = {}
for line in open(repo + "/data/matdb/seed-v1/guitar-spruce-sitka-fpl-gtr282-shell/properties.tsv", encoding="utf-8").read().splitlines():
    f = line.split("\t")
    if f[0] == "scalar":
        props[f[3]] = float(f[4])
rho, EL, ER, nu = props["density"], props["young_modulus_1"], props["young_modulus_2"], props["poisson_ratio_12"]
GLT = props["shear_modulus_13"]
nuRL = nu * ER / EL
BL = [4.730040745, 7.853204624, 10.995607838, 14.137165491]  # roots of cos(x) cosh(x) = 1


def freqs(EI, rhoA):
    return [(x * x) / (2 * math.pi * L * L) * math.sqrt(EI / rhoA) for x in BL]


A1 = b * h
if kind == "tbeam":
    A2 = w * H
    zc = (A2 * (h / 2 + H / 2)) / (A1 + A2)
    I = b * h ** 3 / 12 + A1 * zc ** 2 + w * H ** 3 / 12 + A2 * (h / 2 + H / 2 - zc) ** 2
    rhoA = rho * (A1 + A2)
    ref = freqs(EL * I, rhoA)
    ref_flange = freqs(EL * (I - 0) / 1.0, rhoA)
    I_flange_corr = b * h ** 3 / 12 / (1 - nu * nuRL) + A1 * zc ** 2 + w * H ** 3 / 12 + A2 * (h / 2 + H / 2 - zc) ** 2
    ref_corr = freqs(EL * I_flange_corr, rhoA)
else:
    I = b * h ** 3 / 12
    rhoA = rho * A1
    ref = freqs(EL * I, rhoA)
    ref_corr = freqs(EL * I / (1 - nu * nuRL), rhoA)
print("material: rho %.1f kg/m3, E_L %.4g Pa, E_R %.4g, nu_LR %.3f (nu_RL %.4f), G_LT %.4g" % (rho, EL, ER, nu, nuRL, GLT))
print("section %s: L %.3f m, b %.1f mm, h %.1f mm%s" % (kind, L, b * 1e3, h * 1e3, (", rib %.1f x %.1f mm (height x thickness)" % (H * 1e3, w * 1e3)) if kind == "tbeam" else ""))
print("closed-form free-free bending frequencies (Euler-Bernoulli, E_L I / rho A): " + ", ".join("%.3f" % f for f in ref) + " Hz")
print("with the flange stiffness E_L / (1 - nu_LR nu_RL) (plate-like bracket):        " + ", ".join("%.3f" % f for f in ref_corr) + " Hz")
beta = [x / L for x in BL]
tim = []
for f0, bt in zip(ref, beta):
    r2 = (I / (A1 if kind != "tbeam" else A1 + w * H))
    tim.append(f0 / (1 + 0.5 * bt * bt * r2 * (1 + EL / (5.0 / 6.0 * GLT))))
print("with a Timoshenko estimate (kappa 5/6, G_LT, rotary inertia):                   " + ", ".join("%.3f" % f for f in tim) + " Hz")


# ---- the table
comps, ties = [], []


def plate(name, x0, x1, nxe):
    nodes, tris = [], []
    for j in range(ny + 1):
        for i in range(nxe + 1):
            nodes.append((x0 + (x1 - x0) * i / nxe, -b / 2 + b * j / ny, 0.0, 0.0, 0.0, 1.0))
    idx = lambda i, j: j * (nxe + 1) + i
    for j in range(ny):
        for i in range(nxe):
            a, c, d, e = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            if (i + j) % 2 == 0:
                tris += [(a, c, d), (a, d, e)]
            else:
                tris += [(a, c, e), (c, d, e)]
    comps.append(dict(name=name, t=h, nodes=nodes, tris=tris, axis=(1.0, 0.0, 0.0)))
    return idx


if kind == "split":
    n1 = nx // 2
    idx_a = plate("strip-a", -L / 2, 0.0, n1)
    idx_b = plate("strip-b", 0.0, L / 2, nx - n1)
    ties.append(("strip-a", "strip-b", 0, [idx_a(n1, j) for j in range(ny + 1)], [idx_b(0, j) for j in range(ny + 1)]))
else:
    idx = plate("plate", -L / 2, L / 2, nx)
    if kind == "tbeam":
        nodes, tris = [], []
        for k in range(K + 1):
            for i in range(nx + 1):
                nodes.append((-L / 2 + L * i / nx, 0.0, h / 2 + H * k / K, 0.0, 1.0, 0.0))
        wi = lambda i, k: k * (nx + 1) + i
        for k in range(K):
            for i in range(nx):
                a, c, d, e = wi(i, k), wi(i + 1, k), wi(i + 1, k + 1), wi(i, k + 1)
                # normal +y: with x right and z up cross(x, z) = -y, so the winding is reversed
                if (i + k) % 2 == 0:
                    tris += [(a, d, c), (a, e, d)]
                else:
                    tris += [(a, e, c), (c, e, d)]
        comps.append(dict(name="rib", t=w, nodes=nodes, tris=tris, axis=(1.0, 0.0, 0.0)))
        ties.append(("plate", "rib", 1, [idx(i, 1) for i in range(nx + 1)], [wi(i, 0) for i in range(nx + 1)]))

with open(out, "w", encoding="utf-8") as f:
    f.write("frankensim.shell-assembly.v1\n")
    f.write("meta\tcheck\tindependent closed-form check of the shell and tie machinery (%s, %d elements along the length); closed-form free-free bending frequencies %s Hz; no specimen data\n" % (kind, nx, ", ".join("%.3f" % x for x in ref)))
    for c in comps:
        f.write("component\t%s\t%s\t%.9f\t1\t%d\t%d\t0\t0\n" % (c["name"], pack, c["t"], len(c["nodes"]), len(c["tris"])))
        f.write("region\t%.9f\t%.9f\t%.9f\n" % c["axis"])
        for n in c["nodes"]:
            f.write("node\t%.9f\t%.9f\t%.9f\t%.12f\t%.12f\t%.12f\n" % n)
        for t in c["tris"]:
            f.write("tri\t%d\t%d\t%d\t0\n" % t)
    for m, s, side, mn, sn in ties:
        f.write("tie\t%s\t%s\t%d\t%d\t%d\n" % (m, s, side, len(mn), len(sn)))
        f.write("master\t%s\n" % " ".join(map(str, mn)))
        f.write("slave\t%s\n" % " ".join(map(str, sn)))
    f.write("cavity\t0.001\t0.010\t0.003\n")
    f.write("drive\tmid-length, normal to the strip\t%s\t0.0\t0.0\t0.0\t0\t0\t1\n" % comps[0]["name"] if kind != "split" else "drive\tmid-length, normal to the strip\tstrip-a\t-0.001\t0.0\t0.0\t0\t0\t1\n")
print("wrote %s: %d components, %d nodes, %d triangles, %d ties" % (out, len(comps), sum(len(c["nodes"]) for c in comps), sum(len(c["tris"]) for c in comps), len(ties)))
