"""Out-of-process gmsh/OpenCASCADE adapter for the IRP M1.3-W4 import fixtures.

Run offline (never a build or link dependency, ADPT-2026-07):

    uv run --no-project --with gmsh==4.15.2 python -I make_holed_block.py <out-dir>

Writes, in millimetres (Gmsh files carry no unit; the importer declares it):
- holed_block_tet10.msh     second-order (Gmsh type 11) volume mesh, MSH 2.2 ASCII
- holed_block.step          the OpenCASCADE solid
- holed_block_surface.msh   a 3-node triangle tessellation of that STEP, MSH 2.2 ASCII
- holed_block.receipt.txt   key=value receipt: tool identity, settings, census, volumes

Everything here is Estimate authority. The geometry is a declared synthetic
design: a 24 x 12 x 8 mm block with a 3 mm radius through-hole along z.
"""

import hashlib
import math
import os
import sys

import gmsh

BOX_MM = (24.0, 12.0, 8.0)
HOLE_CENTRE_MM = (12.0, 6.0)
HOLE_RADIUS_MM = 3.0
VOLUME_SETTINGS = {
    "General.NumThreads": 1,
    "Mesh.MaxNumThreads3D": 1,
    "Mesh.Algorithm": 6,
    "Mesh.Algorithm3D": 1,
    "Mesh.MeshSizeMax": 4.0,
    "Mesh.ElementOrder": 2,
    "Mesh.SecondOrderLinear": 0,
    "Mesh.HighOrderOptimize": 0,
    "Mesh.Optimize": 1,
    "Mesh.MshFileVersion": 2.2,
    "Mesh.Binary": 0,
}
SURFACE_SETTINGS = {
    "General.NumThreads": 1,
    "Mesh.Algorithm": 6,
    "Mesh.MeshSizeMax": 3.0,
    "Mesh.ElementOrder": 1,
    "Mesh.MshFileVersion": 2.2,
    "Mesh.Binary": 0,
}


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def settings_fingerprint(settings):
    text = ";".join(f"{k}={settings[k]!r}" for k in sorted(settings))
    return int(hashlib.sha256(text.encode()).hexdigest()[:16], 16)


def sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def min_dihedral_deg(p):
    faces = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
    normals = []
    for k, (a, b, c) in enumerate(faces):
        n = cross(sub(p[b], p[a]), sub(p[c], p[a]))
        if dot(n, sub(p[k], p[a])) > 0.0:
            n = [-v for v in n]
        normals.append(n)
    best = 180.0
    for k in range(4):
        for m in range(k + 1, 4):
            c = dot(normals[k], normals[m]) / math.sqrt(
                dot(normals[k], normals[k]) * dot(normals[m], normals[m])
            )
            best = min(best, 180.0 - math.degrees(math.acos(max(-1.0, min(1.0, c)))))
    return best


def normalize_step_header(path):
    """Put the Part-21 header entities in the order ISO 10303-21 requires
    (FILE_DESCRIPTION, FILE_NAME, FILE_SCHEMA). OpenCASCADE 7.8 writes
    FILE_NAME first, which a conforming reader refuses. Nothing else moves."""
    with open(path) as f:
        text = f.read()
    start = text.index("HEADER;") + len("HEADER;")
    end = text.index("ENDSEC;", start)
    # Split at semicolons outside Part-21 strings ('' escapes a quote).
    statements, current, quoted = [], "", False
    for ch in text[start:end]:
        if ch == "'":
            quoted = not quoted
        if ch == ";" and not quoted:
            statements.append(current.strip())
            current = ""
        else:
            current += ch
    assert not quoted and not current.strip(), current
    statements = [s for s in statements if s]
    order = ["FILE_DESCRIPTION", "FILE_NAME", "FILE_SCHEMA"]
    by_name = {s.split("(", 1)[0].strip(): s for s in statements}
    assert sorted(by_name) == sorted(order), by_name.keys()
    original = [s.split("(", 1)[0].strip() for s in statements]
    header = "".join(f"\n{by_name[name]};" for name in order) + "\n"
    with open(path, "w") as f:
        f.write(text[:start] + header + text[end:])
    return ",".join(original)


def duffy_rule():
    """Collapsed 3-point Gauss-Legendre rule on the unit tetrahedron: exact
    for the cubic det J of an affine-parametrised quadratic tet."""
    gl = [(0.5 - math.sqrt(0.15), 5.0 / 18.0), (0.5, 8.0 / 18.0), (0.5 + math.sqrt(0.15), 5.0 / 18.0)]
    points, weights = [], []
    for a, wa in gl:
        for b, wb in gl:
            for c, wc in gl:
                points += [a, b * (1.0 - a), c * (1.0 - a) * (1.0 - b)]
                weights.append(wa * wb * wc * (1.0 - a) ** 2 * (1.0 - b))
    return points, weights


def apply(settings):
    for key, value in settings.items():
        gmsh.option.setNumber(key, value)


def main(out):
    receipt = {}
    gmsh.initialize(["-noenv"])
    gmsh.option.setNumber("General.Terminal", 0)
    receipt["adapter.tool"] = "gmsh"
    receipt["adapter.version"] = gmsh.option.getString("General.Version")
    receipt["adapter.kernel"] = "OpenCASCADE"
    receipt["adapter.python"] = sys.version.split()[0]
    receipt["adapter.script_sha256"] = sha256(__file__)
    receipt["authority"] = "Estimate (ADPT-2026-07)"
    receipt["unit"] = "mm"

    gmsh.model.add("holed_block")
    box = gmsh.model.occ.addBox(0.0, 0.0, 0.0, *BOX_MM)
    hole = gmsh.model.occ.addCylinder(*HOLE_CENTRE_MM, 0.0, 0.0, 0.0, BOX_MM[2], HOLE_RADIUS_MM)
    solid, _ = gmsh.model.occ.cut([(3, box)], [(3, hole)])
    gmsh.model.occ.synchronize()
    volume_tag = solid[0][1]
    receipt["cad.volume_mm3"] = repr(gmsh.model.occ.getMass(3, volume_tag))
    receipt["cad.exact_volume_mm3"] = repr(
        BOX_MM[0] * BOX_MM[1] * BOX_MM[2] - math.pi * HOLE_RADIUS_MM**2 * BOX_MM[2]
    )
    step = f"{out}/holed_block.step"
    gmsh.write(step)
    receipt["step.raw_sha256"] = sha256(step)
    receipt["step.raw_header_order"] = normalize_step_header(step)
    receipt["step.header_order"] = "FILE_DESCRIPTION,FILE_NAME,FILE_SCHEMA"
    receipt["step.sha256"] = sha256(step)
    receipt["step.bytes"] = str(os.path.getsize(step))

    # Second-order volume mesh: only the volume is a physical group, so the
    # file holds TET10 elements alone.
    gmsh.model.addPhysicalGroup(3, [volume_tag], 1, "block")
    apply(VOLUME_SETTINGS)
    gmsh.model.mesh.generate(3)
    msh = f"{out}/holed_block_tet10.msh"
    gmsh.write(msh)
    receipt["tet10.sha256"] = sha256(msh)
    receipt["tet10.bytes"] = str(os.path.getsize(msh))
    receipt["tet10.settings_fingerprint"] = f"{settings_fingerprint(VOLUME_SETTINGS):016x}"
    types, tags, nodes = gmsh.model.mesh.getElements(3, volume_tag)
    assert list(types) == [11], types
    tags = list(tags[0])
    corner = {}
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    for i, t in enumerate(node_tags):
        corner[int(t)] = [coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]]
    connectivity = list(nodes[0])
    # Curved element volumes from gmsh's own Jacobians; its "volume" quality
    # uses the corners only and is recorded as the straight-sided volume.
    points, weights = duffy_rule()
    _, det, _ = gmsh.model.mesh.getJacobians(11, points, volume_tag)
    n = len(weights)
    volumes = [math.fsum(det[e * n + q] * weights[q] for q in range(n)) for e in range(len(tags))]
    assert min(volumes) > 0.0
    straight = list(gmsh.model.mesh.getElementQualities(tags, "volume"))
    sicn = list(gmsh.model.mesh.getElementQualities(tags, "minSICN"))
    dihedral = min(
        min_dihedral_deg([corner[int(connectivity[10 * e + k])] for k in range(4)])
        for e in range(len(tags))
    )
    receipt["tet10.elements"] = str(len(tags))
    receipt["tet10.nodes"] = str(len(set(int(n) for n in connectivity)))
    receipt["tet10.volume_mm3"] = repr(math.fsum(volumes))
    receipt["tet10.straight_sided_volume_mm3"] = repr(math.fsum(straight))
    receipt["tet10.element_volume_min_mm3"] = repr(min(volumes))
    receipt["tet10.element_volume_max_mm3"] = repr(max(volumes))
    receipt["tet10.min_corner_dihedral_deg"] = repr(dihedral)
    receipt["tet10.min_sicn"] = repr(min(sicn))

    # Identified tessellation of the STEP solid.
    gmsh.model.remove()
    gmsh.model.add("holed_block_surface")
    gmsh.model.occ.importShapes(step)
    gmsh.model.occ.synchronize()
    surfaces = [t for _, t in gmsh.model.getEntities(2)]
    gmsh.model.addPhysicalGroup(2, surfaces, 1, "boundary")
    apply(SURFACE_SETTINGS)
    gmsh.model.mesh.generate(2)
    surface = f"{out}/holed_block_surface.msh"
    gmsh.write(surface)
    receipt["tessellation.sha256"] = sha256(surface)
    receipt["tessellation.bytes"] = str(os.path.getsize(surface))
    receipt["tessellation.settings_fingerprint"] = (
        f"{settings_fingerprint(SURFACE_SETTINGS):016x}"
    )
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    position = {int(t): coords[3 * i : 3 * i + 3] for i, t in enumerate(node_tags)}
    triangles = []
    for s in surfaces:
        types, _, nodes = gmsh.model.mesh.getElements(2, s)
        assert list(types) == [2], types
        triangles += [
            [int(nodes[0][3 * i + k]) for k in range(3)] for i in range(len(nodes[0]) // 3)
        ]
    # Divergence-theorem volume, orientation-agnostic, and the largest gap
    # between a hole-wall triangle's edge midpoints and the cylinder.
    signed = math.fsum(
        dot(position[a], cross(position[b], position[c])) / 6.0 for a, b, c in triangles
    )
    deviation = 0.0
    for tri in triangles:
        ps = [position[v] for v in tri]
        radii = [math.hypot(p[0] - HOLE_CENTRE_MM[0], p[1] - HOLE_CENTRE_MM[1]) for p in ps]
        if all(abs(r - HOLE_RADIUS_MM) < 1e-9 for r in radii):
            for i in range(3):
                m = [(ps[i][k] + ps[(i + 1) % 3][k]) / 2.0 for k in range(3)]
                rho = math.hypot(m[0] - HOLE_CENTRE_MM[0], m[1] - HOLE_CENTRE_MM[1])
                deviation = max(deviation, HOLE_RADIUS_MM - rho)
    receipt["tessellation.triangles"] = str(len(triangles))
    receipt["tessellation.nodes"] = str(len({v for t in triangles for v in t}))
    receipt["tessellation.enclosed_volume_mm3"] = repr(abs(signed))
    receipt["tessellation.max_hole_chord_deviation_mm"] = repr(deviation)
    gmsh.finalize()

    with open(f"{out}/holed_block.receipt.txt", "w") as f:
        for key in sorted(receipt):
            f.write(f"{key}={receipt[key]}\n")


if __name__ == "__main__":
    main(sys.argv[1])
