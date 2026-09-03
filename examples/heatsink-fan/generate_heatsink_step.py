#!/usr/bin/env python3
"""Generate heatsink.step for the heatsink-fan worked example: ONE closed manifold shell
represented as an ISO 10303-21 Part 21 triangular FACETED_BREP.

This provides the CAD / B-rep import path parity with `heatsink.stl` (bead frankensim-rc-root-q61wp.51).
The geometry is identical to `heatsink.stl`: base plate 80 x 60 x 5 mm at the origin with 4 fins
6 x 60 x 20 mm, forming a 108-facet closed 2-manifold shell.

Usage:
    python3 examples/heatsink-fan/generate_heatsink_step.py [OUT] [NFINS]
Default OUT is examples/heatsink-fan/heatsink.step with 4 fins.
"""
import json
import os
import sys

from generate_heatsink_stl import build


def generate_step(facets, vol, nfins, step_path, prov_path=None):
    # Collect unique vertices preserving determinism
    unique_verts = []
    vert_map = {}
    for _, a, b, c in facets:
        for p in (a, b, c):
            key = (round(p[0], 9), round(p[1], 9), round(p[2], 9))
            if key not in vert_map:
                vert_map[key] = len(unique_verts) + 1  # 1-based index
                unique_verts.append(p)

    lines = [
        "ISO-10303-21;",
        "HEADER;",
        "FILE_DESCRIPTION(('heatsink triangular faceted B-rep CAD model'),'2;1');",
        "FILE_NAME('heatsink.step','2026-09-03T00:00:00',('fs-io'),('FrankenSim'),'fs-io','FrankenSim','');",
        "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));",
        "ENDSEC;",
        "DATA;",
    ]

    # Entity numbering layout:
    # Points: 1 .. N_verts
    # Loops:  1000 .. 1000 + N_facets - 1
    # Bounds: 2000 .. 2000 + N_facets - 1
    # Faces:  3000 .. 3000 + N_facets - 1
    # Shell:  5000
    # Brep:   5001 (root_id = 5001)

    for i, p in enumerate(unique_verts, start=1):
        lines.append(f"#{i}=CARTESIAN_POINT('',({p[0]:.6f},{p[1]:.6f},{p[2]:.6f}));")

    face_ids = []
    for i, (_, a, b, c) in enumerate(facets):
        p1 = vert_map[(round(a[0], 9), round(a[1], 9), round(a[2], 9))]
        p2 = vert_map[(round(b[0], 9), round(b[1], 9), round(b[2], 9))]
        p3 = vert_map[(round(c[0], 9), round(c[1], 9), round(c[2], 9))]

        loop_id = 1000 + i
        bound_id = 2000 + i
        face_id = 3000 + i
        face_ids.append(face_id)

        lines.append(f"#{loop_id}=POLY_LOOP('',(#{p1},#{p2},#{p3}));")
        lines.append(f"#{bound_id}=FACE_OUTER_BOUND('',#{loop_id},.T.);")
        lines.append(f"#{face_id}=FACE('',(#{bound_id}));")

    shell_faces_str = ",".join(f"#{fid}" for fid in face_ids)
    shell_id = 5000
    root_id = 5001
    lines.append(f"#{shell_id}=CLOSED_SHELL('',({shell_faces_str}));")
    lines.append(f"#{root_id}=FACETED_BREP('',#{shell_id});")
    lines.append("ENDSEC;")
    lines.append("END-ISO-10303-21;")
    lines.append("")

    content = "\n".join(lines)
    with open(step_path, "w") as f:
        f.write(content)

    provenance = {
        "generator": "examples/heatsink-fan/generate_heatsink_step.py",
        "generator_version": "1.0.0",
        "command": f"python3 examples/heatsink-fan/generate_heatsink_step.py {step_path} {nfins}",
        "created_utc": "2026-09-03T00:00:00Z",
        "entity_type": "FACETED_BREP",
        "root_id": root_id,
        "shell_id": shell_id,
        "vertex_count": len(unique_verts),
        "face_count": len(facets),
        "volume_m3": vol,
        "format": "step",
        "schema": "CONFIG_CONTROL_DESIGN",
        "unit": "m",
        "manifold": True,
    }

    if prov_path:
        with open(prov_path, "w") as f:
            json.dump(provenance, f, indent=2)
            f.write("\n")

    print(
        f"Generated {step_path}: vertices={len(unique_verts)} faces={len(facets)} "
        f"root_id={root_id} volume={vol:.6e} m3"
    )
    return root_id, provenance


def main():
    default_step = os.path.join(os.path.dirname(__file__), "heatsink.step")
    default_prov = os.path.join(os.path.dirname(__file__), "heatsink.step.provenance.json")

    step_out = sys.argv[1] if len(sys.argv) > 1 else default_step
    nfins = int(sys.argv[2]) if len(sys.argv) > 2 else 4

    facets, vol = build(nfins)
    prov_out = step_out + ".provenance.json" if step_out != default_step else default_prov
    generate_step(facets, vol, nfins, step_out, prov_out)


if __name__ == "__main__":
    main()
