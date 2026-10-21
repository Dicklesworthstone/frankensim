"""audit_boundaries.py TABLE : the complement of audit_welds.py for tables whose ties are skin-offset (slave nodes are not coincident with
their master): (1) coincident nodes INSIDE one component (a crack in a component), (2) every boundary node of every component (an edge used by
one triangle) that is in no tie, grouped by component, with the lowest and highest z of those nodes, so a boundary that should be welded and is
not shows up as a component whose free nodes are not the documented free edge (a wall's top row and ends, the sound-hole rim). Numpy-free.
"""
import sys
from collections import defaultdict

lines = open(sys.argv[1], encoding="utf-8").read().splitlines()
comps = {}
ties = []
i = 1
while i < len(lines):
    f = lines[i].split("\t")
    if f[0] == "component":
        name, nr, nn, nt = f[1], int(f[4]), int(f[5]), int(f[6])
        nodes = [tuple(float(x) for x in l.split("\t")[1:4]) for l in lines[i + 1 + nr : i + 1 + nr + nn]]
        tris = [tuple(int(x) for x in l.split("\t")[1:4]) for l in lines[i + 1 + nr + nn : i + 1 + nr + nn + nt]]
        comps[name] = (nodes, tris)
        i += 1 + nr + nn + nt
    elif f[0] == "tie":
        ties.append((f[1], f[2], int(f[3]), [int(x) for x in lines[i + 1].split("\t")[1].split()], [int(x) for x in lines[i + 2].split("\t")[1].split()]))
        i += 3
    else:
        i += 1

tied = defaultdict(set)
for m, s, side, mn, sn in ties:
    tied[m].update(mn)
    tied[s].update(sn)

total_dup = 0
total_free = 0
print("%-34s %7s %7s %9s %9s %s" % ("component", "dups", "bnd", "bnd_free", "z_lo_mm", "z_hi_mm"))
for name, (nodes, tris) in comps.items():
    # coincident nodes inside one component
    seen = defaultdict(list)
    for k, p in enumerate(nodes):
        seen[tuple(round(c, 8) for c in p)].append(k)
    dups = sum(len(v) - 1 for v in seen.values() if len(v) > 1)
    # boundary nodes
    edges = defaultdict(int)
    for (a, b, c) in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            edges[(min(u, v), max(u, v))] += 1
    boundary = set()
    for (u, v), n in edges.items():
        if n == 1:
            boundary.update((u, v))
    free = sorted(boundary - tied[name])
    total_dup += dups
    total_free += len(free)
    zs = [nodes[k][2] * 1000 for k in free]
    print("%-34s %7d %7d %9d %9s %s" % (name, dups, len(boundary), len(free), "%.1f" % min(zs) if zs else "-", "%.1f" % max(zs) if zs else "-"))
print("coincident nodes inside components: %d; boundary nodes in no tie: %d" % (total_dup, total_free))
