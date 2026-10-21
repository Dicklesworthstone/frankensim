"""mesh_resolution.py : how coarse are the published meshes against the 35.66 um voxel, and does decimation change what we measure? For the whole-violin STLs (the file named 200K and
the one named 10M) and the plate and rib meshes: triangle count, edge-length statistics (median, 10th, 90th percentile, maximum, in mm and in voxels), surface area and enclosed
volume, and the area of the cross-section by a few planes (z = -100, z = -190, y = 5, x = 30) for the two STLs. Read-only; numpy only."""
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cannone_io as io  # noqa: E402
import cannone_sections as cs  # noqa: E402

VOX = 0.03566
R = "/data/tmp/wa-inst/geometry/violin-cannone-esrf/raw"
FILES = [
    ("STL named 200K", R + "/stage1/2416550404_Cannone_FS2-accorded_surface_200K.stl"),
    ("STL named 10M", R + "/stage3a/2416550458_Cannone_FS2-accorded_surface_10M.stl"),
    ("top plate mesh", R + "/stage2/2416550405_table_all.glb"),
    ("back plate mesh", R + "/stage2/2416550427_dos_all.glb"),
    ("ribs mesh", R + "/stage2/2416550419_eclisses_all.glb"),
    ("soundpost mesh", R + "/stage1/2416550418_ame.glb"),
]
sec_planes = ((2, -100.0), (2, -190.0), (1, 5.0), (0, 30.0))
out = {}
for name, path in FILES:
    if not os.path.exists(path):
        print("==", name, "not fetched")
        continue
    v, t, info = io.stl_triangles(path) if path.endswith(".stl") else io.glb_triangles(path)
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    e = np.concatenate([np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1), np.linalg.norm(a - c, axis=1)])
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
    vol = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0
    p10, p50, p90 = np.percentile(e, [10, 50, 90])
    print("== %s: %d triangles; edge mm median %.3f (%.1f voxels), p10 %.3f, p90 %.3f, max %.2f; area %.0f mm2; volume %.0f mm3" % (name, len(t), p50, p50 / VOX, p10, p90, e.max(), area, vol), flush=True)
    if "STL" in name:
        secs = []
        for ax, s in sec_planes:
            rows = cs.sections(v, t, ax, 1.0, lo=s, hi=s + 0.5)
            secs.append(rows[0][1] if rows else float("nan"))
        out[name] = secs
        print("   section areas (planes %s): %s" % (", ".join("%s=%g" % ("xyz"[ax], s) for ax, s in sec_planes), ", ".join("%.1f" % q for q in secs)), flush=True)
if len(out) == 2:
    k = list(out)
    print("section area ratio 10M / 200K:", ", ".join("%.5f" % (y / x) for x, y in zip(out[k[0]], out[k[1]])))
