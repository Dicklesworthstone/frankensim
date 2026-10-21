"""Plot the vertical displacement of the top and the back for the first air-coupled modes written by guitar_body_modes --shapes.
usage: python3 -I plot_modes.py SHAPES.tsv OUT.png [N_MODES]
(uv run --with matplotlib --with numpy)
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

src, out = sys.argv[1], sys.argv[2]
nm_max = int(sys.argv[3]) if len(sys.argv) > 3 else 9
rows = [l.split("\t") for l in open(src, encoding="utf-8").read().splitlines()]
nm = min(nm_max, len(rows[0]) - 4)
data = {"top": [], "back": []}
for r in rows[1:]:
    data[r[0]].append([float(v) for v in r[2:]])
top = np.array(data["top"])
back = np.array(data["back"])
fig, axes = plt.subplots(nm, 2, figsize=(9, 2.2 * nm))
for m in range(nm):
    for j, (name, arr) in enumerate((("top", top), ("back", back))):
        ax = axes[m][j]
        x, y, w = arr[:, 0], arr[:, 1], arr[:, 2 + m]
        tri = mtri.Triangulation(x, y)
        lim = np.max(np.abs(w)) or 1.0
        ax.tricontourf(tri, w, levels=np.linspace(-lim, lim, 21), cmap="RdBu_r")
        ax.tricontour(tri, w, levels=[0.0], colors="k", linewidths=0.6)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("mode %d %s (max |w| %.2e)" % (m + 1, name, lim), fontsize=7)
fig.tight_layout()
fig.savefig(out, dpi=70)
print("wrote", out)
