"""check_shell_assembly.py TABLE : audit of a frankensim.shell-assembly.v1 table before it is run: per component the facet-against-director alignment (the example refuses below 0.9) and the
triangle quality, per tie the gap between each slave node and its master partner (master node moved by side x thickness / 2 along the master director; exact pairs only), the enclosed
cavity volume, and the mass of every component from the pack densities given on the command line (name=kg/m3). Standard library and numpy."""
import sys

import numpy as np

path = sys.argv[1]
dens = {a.split("=")[0]: float(a.split("=")[1]) for a in sys.argv[2:]}
lines = open(path, encoding="utf-8").read().splitlines()
assert lines[0] == "frankensim.shell-assembly.v1", lines[0]
comps, ties_, cavity = {}, [], None
i = 1
while i < len(lines):
    f = lines[i].split("\t")
    if f[0] == "meta":
        i += 1
    elif f[0] == "component":
        name, pack, t, nr, nn, nt = f[1], f[2], float(f[3]), int(f[4]), int(f[5]), int(f[6])
        nodes = np.array([[float(x) for x in lines[i + 1 + nr + k].split("\t")[1:]] for k in range(nn)])
        tris = np.array([[int(x) for x in lines[i + 1 + nr + nn + k].split("\t")[1:4]] for k in range(nt)])
        comps[name] = dict(pack=pack, t=t, P=nodes[:, :3], D=nodes[:, 3:], T=tris, cav=f[8] == "1", ext=f[7] == "1")
        i += 1 + nr + nn + nt
    elif f[0] == "tie":
        ties_.append((f[1], f[2], int(f[3]), [int(x) for x in lines[i + 1].split("\t")[1].split()], [int(x) for x in lines[i + 2].split("\t")[1].split()]))
        i += 3
    elif f[0] == "drive":
        i += 1
    elif f[0] == "cavity":
        cavity = [float(x) for x in f[1:]]
        i += 1
    else:
        raise SystemExit("unknown record " + f[0])

print("components %d, ties %d, nodes %d, triangles %d" % (len(comps), len(ties_), sum(len(c["P"]) for c in comps.values()), sum(len(c["T"]) for c in comps.values())))
worst_dot, total_mass, vol, worst_pair_all = 1.0, 0.0, 0.0, 1.0
for name, c in sorted(comps.items()):
    P, T, Dn = c["P"], c["T"], c["D"]
    n = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    area = 0.5 * np.linalg.norm(n, axis=1)
    fn = n / np.linalg.norm(n, axis=1, keepdims=True)
    dots = np.array([[fn[t] @ Dn[v] for v in T[t]] for t in range(len(T))])
    e = [np.linalg.norm(P[T[:, (k + 1) % 3]] - P[T[:, k]], axis=1) for k in range(3)]
    ang = []
    for k in range(3):
        e1, e2 = P[T[:, (k + 1) % 3]] - P[T[:, k]], P[T[:, (k + 2) % 3]] - P[T[:, k]]
        ang.append(np.degrees(np.arccos(np.clip((e1 * e2).sum(1) / (np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1)), -1, 1))))
    # the owner also refuses two facets that share a node when their normals differ by more than the same 0.9 in cosine
    inc = [[] for _ in range(len(P))]
    for t_i, tri in enumerate(T):
        for v in tri:
            inc[v].append(fn[t_i])
    worst_pair = 1.0
    for lst in inc:
        for a_ in range(len(lst)):
            for b_ in range(a_ + 1, len(lst)):
                worst_pair = min(worst_pair, float(lst[a_] @ lst[b_]))
    worst_pair_all = min(worst_pair_all, worst_pair)
    # the owner refuses a component that is not one edge-connected body
    parent = list(range(len(T)))

    def find(a_):
        while parent[a_] != a_:
            parent[a_] = parent[parent[a_]]
            a_ = parent[a_]
        return a_

    owner = {}
    for t_i, tri in enumerate(T):
        for k in range(3):
            e = (min(tri[k], tri[(k + 1) % 3]), max(tri[k], tri[(k + 1) % 3]))
            if e in owner:
                parent[find(t_i)] = find(owner[e])
            else:
                owner[e] = t_i
    # one fan at every vertex: the triangles round a vertex must be connected through edges that contain it
    vt = [[] for _ in range(len(P))]
    for t_i, tri in enumerate(T):
        for v in tri:
            vt[v].append(t_i)
    multi = 0
    for v, lst in enumerate(vt):
        seen_, comps_n = set(), 0
        for a_ in lst:
            if a_ in seen_:
                continue
            comps_n += 1
            stack = [a_]
            seen_.add(a_)
            while stack:
                x_ = stack.pop()
                for y_ in lst:
                    if y_ not in seen_ and len(set(T[x_].tolist()) & set(T[y_].tolist())) >= 2:
                        seen_.add(y_)
                        stack.append(y_)
        if comps_n > 1:
            multi += 1
    if multi:
        print("NOT ONE FAN AT EVERY VERTEX: %s has %d vertices with several fans" % (name, multi))
    pieces = len({find(t_i) for t_i in range(len(T))})
    if pieces != 1:
        print("NOT EDGE-CONNECTED: %s has %d pieces" % (name, pieces))
    mass = area.sum() * c["t"] * dens.get(c["pack"], float("nan"))
    total_mass += mass if mass == mass else 0.0
    worst_dot = min(worst_dot, dots.min())
    if c["cav"]:
        vol += float(np.sum(((P[T[:, 0]] + P[T[:, 1]] + P[T[:, 2]]) / 3.0 * (n / 2.0)).sum(1)) / 3.0)
    print("%-14s %-7s t=%.3f mm tris %4d nodes %4d area %9.1f mm2 min facet.director %.3f min angle %5.1f deg min facet.facet %.3f mass %.2f g" % (name, c["pack"], c["t"] * 1e3, len(T), len(P), area.sum() * 1e6, dots.min(), np.min(ang), worst_pair, mass * 1e3))
print("worst facet.director %.3f, worst facet.facet at a shared node %.3f (the example refuses below 0.9 for both); total mass of the modelled components %.1f g (pack densities given: %s)" % (worst_dot, worst_pair_all, total_mass * 1e3, ", ".join("%s=%g" % kv for kv in dens.items())))
print("midsurface cavity volume from the table: %.4f L (cavity record %.4f L, neck radius %.3f mm, length %.3f mm)" % (vol * 1e3, cavity[0] * 1e3, cavity[1] * 1e3, cavity[2] * 1e3))
gaps = []
for m, s, side, mn, sn in ties_:
    cm, cs = comps[m], comps[s]
    if len(mn) != len(sn):
        print("tie %s -> %s: %d master vs %d slave nodes (not paired; skipped)" % (m, s, len(mn), len(sn)))
        continue
    expect = cm["P"][mn] + side * 0.5 * cm["t"] * cm["D"][mn]
    g = np.linalg.norm(cs["P"][sn] - expect, axis=1)
    gaps.append((g.max(), m, s, side, len(mn)))
gaps.sort(reverse=True)
print("tie gaps (largest 8 of %d): " % len(gaps) + "; ".join("%s -> %s side %+d n=%d: %.2e m" % (m, s, side, n_, g) for g, m, s, side, n_ in gaps[:8]))
print("largest gap over all paired ties: %.3e m (the example's tolerance is 2e-4 m)" % (gaps[0][0] if gaps else 0.0))
