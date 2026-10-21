"""audit_welds.py TABLE : coincident nodes of different components (same position to 1e-9 m) and whether the table welds them: every such group must have all its nodes but one as slave nodes of ties
(any side-0 tie between the components). Prints the number of coincident groups, of unwelded nodes, and by component pair. Numpy-free."""
import sys
from collections import defaultdict

lines = open(sys.argv[1], encoding="utf-8").read().splitlines()
comps, ties = {}, []
i = 1
while i < len(lines):
    f = lines[i].split("\t")
    if f[0] == "component":
        name, nr, nn, nt = f[1], int(f[4]), int(f[5]), int(f[6])
        nodes = [tuple(float(x) for x in l.split("\t")[1:4]) for l in lines[i + 1 + nr:i + 1 + nr + nn]]
        comps[name] = nodes
        i += 1 + nr + nn + nt
    elif f[0] == "tie":
        ties.append((f[1], f[2], int(f[3]), [int(x) for x in lines[i + 1].split("\t")[1].split()], [int(x) for x in lines[i + 2].split("\t")[1].split()]))
        i += 3
    else:
        i += 1
slave = set()
master = set()
for m, s, side, mn, sn in ties:
    for a in mn:
        master.add((m, a))
    for b in sn:
        slave.add((s, b))
groups = defaultdict(list)
for name, nodes in comps.items():
    for k, p in enumerate(nodes):
        groups[tuple(round(c, 8) for c in p)].append((name, k))
coincident = [g for g in groups.values() if len({n for n, _ in g}) > 1]
unwelded = 0
pairs = defaultdict(int)
for g in coincident:
    free = [x for x in g if x not in slave]
    if len(free) > 1:
        unwelded += len(free) - 1
        for x in free[1:]:
            pairs[(free[0][0], x[0])] += 1
print("coincident groups %d; nodes left unwelded (not a slave of any tie, beyond one per group) %d of %d duplicate nodes" % (len(coincident), unwelded, sum(len(g) - 1 for g in coincident)))
for k, v in sorted(pairs.items(), key=lambda kv: -kv[1])[:15]:
    print("  %s ~ %s: %d" % (k[0], k[1], v))
