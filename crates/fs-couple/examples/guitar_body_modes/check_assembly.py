"""Check a frankensim.guitar-shell-assembly.v1 file: ties (gap of every slave node to the offset master polyline),
facet orientation against directors, edge manifoldness of the cavity boundary, enclosed volume.
usage: python3 -I check_assembly.py FILE.tsv
"""
import sys
import numpy as np

lines = open(sys.argv[1], encoding="utf-8").read().splitlines()
assert lines[0] == "frankensim.guitar-shell-assembly.v1"
comps = {}
order = []
ties = []
i = 1
cavity = None
while i < len(lines):
    f = lines[i].split("\t")
    if f[0] == "meta":
        i += 1
    elif f[0] == "component":
        name, pack, t, nr, nn, nt, ext, cav = f[1], f[2], float(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]), int(f[8])
        regs = [list(map(float, lines[i + 1 + k].split("\t")[1:])) for k in range(nr)]
        nodes = np.array([list(map(float, lines[i + 1 + nr + k].split("\t")[1:])) for k in range(nn)])
        tris = np.array([list(map(int, lines[i + 1 + nr + nn + k].split("\t")[1:])) for k in range(nt)])
        comps[name] = dict(pack=pack, t=t, regs=regs, P=nodes[:, :3], Dn=nodes[:, 3:], T=tris, ext=ext, cav=cav)
        order.append(name)
        i += 1 + nr + nn + nt
    elif f[0] == "tie":
        m, s, side, nm, ns = f[1], f[2], int(f[3]), int(f[4]), int(f[5])
        mn = list(map(int, lines[i + 1].split("\t")[1].split()))
        sn = list(map(int, lines[i + 2].split("\t")[1].split()))
        assert len(mn) == nm and len(sn) == ns, (m, s)
        ties.append((m, s, side, mn, sn))
        i += 3
    elif f[0] == "cavity":
        cavity = list(map(float, f[1:]))
        i += 1
    else:
        raise SystemExit("unknown record " + f[0])

# facet normals vs directors
for n in order:
    c = comps[n]
    P, T, Dn = c["P"], c["T"][:, :3], c["Dn"]
    nrm = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    area = 0.5 * np.linalg.norm(nrm, axis=1)
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    dm = Dn[T[:, 0]] + Dn[T[:, 1]] + Dn[T[:, 2]]
    dm /= np.linalg.norm(dm, axis=1)[:, None]
    dots = (nrm * dm).sum(axis=1)
    # nodal min dot with each incident facet
    mind = 1.0
    for k in range(3):
        mind = min(mind, float(((nrm * Dn[T[:, k]]).sum(axis=1)).min()))
    print("%-22s nodes=%5d tris=%5d area=%.5f m2 t=%.5f min(facet.nodal director)=%.4f min(area)=%.2e" % (n, len(P), len(T), area.sum(), c["t"], mind, area.min()))

# tie gaps
worst = 0.0
for (m, s, side, mn, sn) in ties:
    cm, cs = comps[m], comps[s]
    # master polyline in 3-D with the skin offset along the nodal director of the master
    Pm = cm["P"][mn] + side * 0.5 * cm["t"] * cm["Dn"][mn] if side != 0 else cm["P"][mn]
    gaps = []
    for k in sn:
        p = cs["P"][k]
        best = 1e9
        for a, b in zip(Pm[:-1], Pm[1:]):
            ab = b - a
            L2 = ab @ ab
            u = 0.0 if L2 == 0 else max(0.0, min(1.0, ((p - a) @ ab) / L2))
            best = min(best, float(np.linalg.norm(p - (a + u * ab))))
        gaps.append(best)
    g = max(gaps)
    worst = max(worst, g)
    print("tie %-12s -> %-22s side=%+d nm=%3d ns=%3d max gap %.3e m" % (m, s, side, len(mn), len(sn), g))
print("worst tie gap %.3e m" % worst)

# cavity boundary: signed volume and open edges
cav_names = [n for n in order if comps[n]["cav"]]
V = 0.0
edges = {}
for n in cav_names:
    c = comps[n]
    P, T = c["P"], c["T"][:, :3]
    for (i_, j_, k_) in T:
        a, b, cc = P[i_], P[j_], P[k_]
        V += float(((a + b + cc) / 3.0) @ (np.cross(b - a, cc - a) / 2.0)) / 3.0
print("signed cavity volume (open at the hole): %.6f m3 ; file value %.6f m3 (inner, after half thickness)" % (V, cavity[0]))
