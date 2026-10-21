"""scale_table.py IN.tsv OUT.tsv S [K0_1 K0_2 A_1 A_2] : the declared plan-view scale sensitivity of a shell assembly table WITHOUT remeshing: every node keeps its topology and its height, x and z are
scaled by S about the body centre (x 0, z -0.087 m), the nodal directors follow the surface (n -> (n_x / S, n_y, n_z / S) normalised), the slave nodes of the skin ties are put back exactly on the
master skin (the owner needs them there), the cavity volume goes with S^2 and the neck is recomputed from the f-hole conductance (thin-plate conductance K0 of each hole times S in series with
t / (A S^2); the radius from K = pi r^2 / (t + 1.6432 r) with the length t), the drive points are moved with the plan. Default K0 and A are those printed in the table's neck and f_holes rows
of the base table. A remeshed scaled table changes the zones, the ties and the mesh together and so confounds the scale with the discretisation. Standard library only."""
import math
import sys

src, dst, S = sys.argv[1], sys.argv[2], float(sys.argv[3])
K0 = [float(sys.argv[4]), float(sys.argv[5])] if len(sys.argv) > 5 else [32.7, 28.8]
AH = [float(sys.argv[6]), float(sys.argv[7])] if len(sys.argv) > 7 else [403.0, 309.0]
CZ = -0.087
lines = open(src, encoding="utf-8").read().splitlines()
out = []
comps = {}
i = 1
head = lines[0]
records = []  # (kind, payload)
while i < len(lines):
    f = lines[i].split("\t")
    if f[0] == "component":
        name, nr, nn, nt = f[1], int(f[4]), int(f[5]), int(f[6])
        reg = lines[i + 1:i + 1 + nr]
        nodes = [[float(x) for x in l.split("\t")[1:]] for l in lines[i + 1 + nr:i + 1 + nr + nn]]
        tris = lines[i + 1 + nr + nn:i + 1 + nr + nn + nt]
        comps[name] = dict(head=lines[i], t=float(f[3]), reg=reg, nodes=nodes, tris=tris)
        records.append(("component", name))
        i += 1 + nr + nn + nt
    elif f[0] == "tie":
        records.append(("tie", (f, lines[i + 1], lines[i + 2])))
        i += 3
    else:
        records.append((f[0], lines[i]))
        i += 1


def unit(v):
    n = math.sqrt(sum(a * a for a in v))
    return [a / n for a in v]


# 1. scale the plan view and the directors
for c in comps.values():
    for nd in c["nodes"]:
        nd[0] = nd[0] * S
        nd[2] = CZ + (nd[2] - CZ) * S
        d = unit([nd[3] / S, nd[4], nd[5] / S])
        nd[3], nd[4], nd[5] = d
    c["reg"] = [("region\t%.9f\t%.9f\t%.9f" % tuple(unit([float(x) for x in r.split("\t")[1:4]]))) for r in c["reg"]]
# 2. skin ties: put the slave nodes exactly on the master skin
moved = 0.0
for kind, payload in records:
    if kind != "tie":
        continue
    f, mline, sline = payload
    master, slave, side = f[1], f[2], int(f[3])
    mn = [int(x) for x in mline.split("\t")[1].split()]
    sn = [int(x) for x in sline.split("\t")[1].split()]
    assert len(mn) == len(sn)
    if side == 0:
        continue
    cm, cs = comps[master], comps[slave]
    for a, b in zip(mn, sn):
        pm, dm = cm["nodes"][a], cm["nodes"][a][3:]
        new = [pm[k] + side * 0.5 * cm["t"] * dm[k] for k in range(3)]
        moved = max(moved, math.dist(new, cs["nodes"][b][:3]))
        cs["nodes"][b][:3] = new
# 3. cavity and neck, drives
for k, (kind, payload) in enumerate(records):
    if kind == "cavity":
        f = payload.split("\t")
        v = float(f[1]) * S * S
        t_mm = float(f[3]) * 1000.0
        kk = sum(1.0 / (1.0 / (K0[h] * S) + t_mm / (AH[h] * S * S)) for h in range(2))
        r = (1.6432 * kk + math.sqrt((1.6432 * kk) ** 2 + 4 * math.pi * kk * t_mm)) / (2 * math.pi)
        records[k] = ("cavity", "cavity\t%.9f\t%.9f\t%.9f" % (v, r / 1000.0, t_mm / 1000.0))
    elif kind == "drive":
        f = payload.split("\t")
        x, y, z = float(f[3]) * S, float(f[4]), CZ + (float(f[5]) - CZ) * S
        records[k] = ("drive", "\t".join(f[:3] + ["%.6f" % x, "%.6f" % y, "%.6f" % z] + f[6:]))
    elif kind == "meta":
        pass
print("largest move of a slave node onto its master skin: %.3e m; scale %.4f" % (moved, S))
with open(dst, "w", encoding="utf-8") as fo:
    fo.write(head + "\n")
    fo.write("meta\tscale_table\tthe base table scaled in plan by %.4f about the body centre by scale_table.py with its topology, heights and ties unchanged; directors, skin-tie slave nodes, cavity volume, neck and drive points follow\n" % S)
    for kind, payload in records:
        if kind == "component":
            c = comps[payload]
            fo.write(c["head"] + "\n")
            for r in c["reg"]:
                fo.write(r + "\n")
            for nd in c["nodes"]:
                fo.write("node\t%.9f\t%.9f\t%.9f\t%.12f\t%.12f\t%.12f\n" % tuple(nd))
            for t in c["tris"]:
                fo.write(t + "\n")
        elif kind == "tie":
            f, m, s_ = payload
            fo.write("\t".join(f) + "\n" + m + "\n" + s_ + "\n")
        else:
            fo.write(payload + "\n")
