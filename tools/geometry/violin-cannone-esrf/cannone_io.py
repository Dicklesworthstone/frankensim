"""cannone_io.py : readers (numpy only) for the il Cannone GLB and STL files;
the original probe description was: read-only probe of GLB (glTF binary) and STL meshes with numpy only. For each file: format facts (glTF asset/extensions/node transforms, or STL kind),
vertex/triangle counts, bounding box, surface area, signed enclosed volume, and edge-manifold statistics (edges used once / twice / more). Writes nothing."""
import json
import struct
import sys

import numpy as np

COMP = {5120: np.int8, 5121: np.uint8, 5122: np.int16, 5123: np.uint16, 5125: np.uint32, 5126: np.float32}
NCOMP = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4, "MAT4": 16}


def read_glb(path):
    b = open(path, "rb").read()
    magic, ver, total = struct.unpack("<4sII", b[:12])
    assert magic == b"glTF", magic
    pos, js, binc = 12, None, None
    while pos + 8 <= len(b):
        ln, ty = struct.unpack("<I4s", b[pos:pos + 8])
        body = b[pos + 8:pos + 8 + ln]
        if ty == b"JSON":
            js = json.loads(body.decode("utf-8"))
        elif ty == b"BIN\x00":
            binc = body
        pos += 8 + ln
    return js, binc


def accessor(js, binc, idx):
    a = js["accessors"][idx]
    bv = js["bufferViews"][a["bufferView"]]
    dt = np.dtype(COMP[a["componentType"]])
    n = NCOMP[a["type"]]
    off = bv.get("byteOffset", 0) + a.get("byteOffset", 0)
    stride = bv.get("byteStride")
    cnt = a["count"]
    if stride and stride != dt.itemsize * n:
        raw = np.frombuffer(binc, dtype=np.uint8, count=stride * cnt, offset=off).reshape(cnt, stride)
        arr = np.ascontiguousarray(raw[:, :dt.itemsize * n]).view(dt).reshape(cnt, n)
    else:
        arr = np.frombuffer(binc, dtype=dt, count=cnt * n, offset=off).reshape(cnt, n)
    return arr if n > 1 else arr[:, 0]


def node_matrix(node):
    m = np.eye(4)
    if "matrix" in node:
        return np.array(node["matrix"], dtype=np.float64).reshape(4, 4).T
    t = np.array(node.get("translation", [0, 0, 0]), dtype=np.float64)
    q = np.array(node.get("rotation", [0, 0, 0, 1]), dtype=np.float64)
    s = np.array(node.get("scale", [1, 1, 1]), dtype=np.float64)
    x, y, z, w = q
    r = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    m[:3, :3] = r * s
    m[:3, 3] = t
    return m


def glb_triangles(path):
    js, binc = read_glb(path)
    info = {"asset": js.get("asset"), "extensionsUsed": js.get("extensionsUsed"), "extensionsRequired": js.get("extensionsRequired"),
            "nodes": len(js.get("nodes", [])), "meshes": len(js.get("meshes", [])), "materials": [m.get("name") for m in js.get("materials", [])][:8]}
    verts, tris, base = [], [], 0
    # world transform of each mesh instance through the scene graph
    parent = {}
    for i, nd in enumerate(js.get("nodes", [])):
        for c in nd.get("children", []):
            parent[c] = i

    def world(i):
        m = node_matrix(js["nodes"][i])
        while i in parent:
            i = parent[i]
            m = node_matrix(js["nodes"][i]) @ m
        return m

    tf = []
    for i, nd in enumerate(js.get("nodes", [])):
        if "mesh" in nd:
            M = world(i)
            tf.append((nd.get("name"), nd["mesh"], M))
            for prim in js["meshes"][nd["mesh"]]["primitives"]:
                if prim.get("mode", 4) != 4:
                    continue
                if "extensions" in prim and "KHR_draco_mesh_compression" in prim["extensions"]:
                    raise RuntimeError("Draco-compressed primitive: not decodable with numpy only")
                p = accessor(js, binc, prim["attributes"]["POSITION"]).astype(np.float64)
                p = p @ M[:3, :3].T + M[:3, 3]
                if "indices" in prim:
                    ix = accessor(js, binc, prim["indices"]).astype(np.int64).reshape(-1, 3)
                else:
                    ix = np.arange(len(p)).reshape(-1, 3)
                verts.append(p)
                tris.append(ix + base)
                base += len(p)
    info["mesh_nodes"] = [(n, mi, np.round(M[:3, :3].diagonal(), 6).tolist(), np.round(M[:3, 3], 4).tolist()) for n, mi, M in tf][:6]
    return np.concatenate(verts), np.concatenate(tris), info


def stl_triangles(path):
    b = open(path, "rb").read()
    if len(b) >= 84 and 84 + 50 * struct.unpack("<I", b[80:84])[0] == len(b):
        n = struct.unpack("<I", b[80:84])[0]
        raw = np.frombuffer(b, dtype=np.uint8, count=50 * n, offset=84).reshape(n, 50)
        v = np.ascontiguousarray(raw[:, 12:48]).view("<f4").reshape(n, 3, 3).astype(np.float64).reshape(-1, 3)
        return v, np.arange(len(v)).reshape(-1, 3), {"kind": "binary STL", "header": b[:80].rstrip(b"\x00 ").decode("latin1")[:60], "triangles": n}
    txt = b.decode("latin1")
    vs = np.array([[float(x) for x in ln.split()[1:4]] for ln in txt.splitlines() if ln.strip().startswith("vertex")])
    return vs, np.arange(len(vs)).reshape(-1, 3), {"kind": "ASCII STL", "triangles": len(vs) // 3}


def stats(v, t):
    a, b_, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    cr = np.cross(b_ - a, c - a)
    area = 0.5 * np.linalg.norm(cr, axis=1).sum()
    vol = np.einsum("ij,ij->i", a, np.cross(b_, c)).sum() / 6.0
    # weld vertices exactly (STL has none; GLB may split) before the edge count
    q = np.round(v * 1e3).astype(np.int64) + (1 << 20)  # micrometre grid, coordinates in mm, |x| < 1000 mm
    assert q.min() >= 0 and q.max() < (1 << 21), "coordinates outside the weld grid"
    key = (q[:, 0] << 42) | (q[:, 1] << 21) | q[:, 2]
    _, inv = np.unique(key, return_inverse=True)
    inv = inv.reshape(-1)
    tt = inv[t]
    e = np.concatenate([tt[:, [0, 1]], tt[:, [1, 2]], tt[:, [2, 0]]])
    e.sort(axis=1)
    _, cnt = np.unique(e[:, 0] * (inv.max() + 1) + e[:, 1], return_counts=True)
    return area, vol, {int(k): int(np.sum(cnt == k)) for k in (1, 2)} | {"3+": int(np.sum(cnt >= 3))}, int(inv.max() + 1)
