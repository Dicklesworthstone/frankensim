#!/usr/bin/env python3
"""Level-B cross-code reference deck runner (scikit-fem side).

This script is the EXTERNAL half of the FrankenSim Level-B thermal
cross-code corpus (bead frankensim-extreal-program-f85xj.4.3). It solves
steady heat conduction cases described by committed case-spec JSON decks
using scikit-fem + scipy (SuperLU) — an implementation that shares no
code with the FrankenSim workspace — and freezes probe temperatures into
a TSV reference manifest that `fs-vvreg` parses fail-closed.

Reproducibility contract
------------------------
* The environment is pinned by `tools/vvref/pyproject.toml` + `uv.lock`
  (committed). Re-derive with:
      uv run --project tools/vvref python tools/vvref/solve_skfem.py \
          --cases <case.json>... --out <references.tsv>
* The mesh is the exact Kuhn/Freudenthal structured subdivision used by
  `fs_conduction::fixtures::structured_tets` (same vertex indexing, same
  6 permutations, same odd-permutation orientation swap, same IEEE-754
  position arithmetic `i * (extent/n)`), witnessed by a BLAKE3 hash over
  canonical mesh bytes that the Rust comparison lane recomputes.
* The discrete system matches fs-conduction's documented quadrature
  claim (crates/fs-conduction/src/assemble.rs): P1 stiffness with K
  evaluated once per element at the element-mean temperature, consistent
  tet source mass V(1+delta)/20, consistent triangle Robin mass
  A(1+delta)/12. Agreement is therefore same-discretization cross-code
  parity: it checks independent ASSEMBLY/BC/SOLVER implementations, not
  discretization error.
* Two codes agreeing is NOT physical truth: the corpus registers these
  rows as Estimated-colour Level-B (CrossCode) references.

Self-checks (known-answer probes proving the harness plumbing) run
before any case is solved; a failure aborts the freeze.
"""

from __future__ import annotations

import argparse
import json
import platform
import struct
import sys

import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import skfem
from blake3 import blake3
from skfem import Basis, ElementTetP1, FacetBasis, MeshTet
from skfem.assembly import BilinearForm, LinearForm
from skfem.helpers import dot, grad

# The 6 axis orderings of fs_conduction::fixtures::PERMS, in order.
PERMS = [
    (0, 1, 2),
    (0, 2, 1),
    (1, 0, 2),
    (1, 2, 0),
    (2, 0, 1),
    (2, 1, 0),
]
# Permutations whose Kuhn tet needs the orientation swap in the Rust
# fixture (`odd` in structured_tets).
ODD_PERMS = {(0, 2, 1), (1, 0, 2), (2, 1, 0)}

FACE_AXES = {"x": 0, "y": 1, "z": 2}
FACE_BAND = 1e-9  # fs_conduction::fixtures::on_box_face band.


def kuhn_box_grid(counts, extent):
    """Port of fs_conduction::fixtures::box_grid, bit-identical positions.

    Returns (positions[(n,3) float64], tets[(m,4) int64]).
    """
    nx, ny, nz = counts
    ex, ey, ez = extent
    hx = ex / nx
    hy = ey / ny
    hz = ez / nz
    px, py, pz = nx + 1, ny + 1, nz + 1
    positions = np.empty((px * py * pz, 3), dtype=np.float64)
    row = 0
    for i in range(px):
        for j in range(py):
            for k in range(pz):
                # Scalar IEEE ops in the same order as the Rust fixture.
                positions[row, 0] = i * hx
                positions[row, 1] = j * hy
                positions[row, 2] = k * hz
                row += 1

    def idx(i, j, k):
        return i * py * pz + j * pz + k

    tets = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for perm in PERMS:
                    corners = [(i, j, k)]
                    for axis in perm:
                        prev = corners[-1]
                        nxt = list(prev)
                        nxt[axis] += 1
                        corners.append(tuple(nxt))
                    v = [idx(*c) for c in corners]
                    if perm in ODD_PERMS:
                        tets.append([v[0], v[2], v[1], v[3]])
                    else:
                        tets.append([v[0], v[1], v[2], v[3]])
    return positions, np.asarray(tets, dtype=np.int64)


def mesh_blake3(counts, positions, tets):
    """Canonical mesh identity: recomputed by the Rust comparison lane."""
    h = blake3()
    h.update(b"frankensim-vvref-mesh-v1\n")
    for c in counts:
        h.update(struct.pack("<Q", c))
    h.update(struct.pack("<Q", positions.shape[0]))
    for row in positions:
        for value in row:
            h.update(struct.pack("<d", value))
    h.update(struct.pack("<Q", tets.shape[0]))
    for tet in tets:
        for v in tet:
            h.update(struct.pack("<I", v))
    return h.hexdigest()


def face_value(spec, extent):
    axis = FACE_AXES[spec[0]]
    target = extent[axis] if spec[1] == "max" else float(spec[1])
    return axis, target


def nodes_on_face(positions, axis, target):
    return np.where(np.abs(positions[:, axis] - target) < FACE_BAND)[0]


def conductivity_tensor(material, tbar=None):
    """Per-element 3x3 conductivity for the declared material model."""
    kind = material["kind"]
    if kind == "isotropic":
        k = float(material["k"])
        return np.eye(3) * k
    if kind == "tensor":
        return np.asarray(material["k"], dtype=np.float64)
    raise ValueError(f"tensor form unsupported for material kind {kind}")


def kt_eval(knots, t):
    """Piecewise-linear k(T); refuses outside the sampled span."""
    ts = np.asarray([k[0] for k in knots], dtype=np.float64)
    ks = np.asarray([k[1] for k in knots], dtype=np.float64)
    if np.any(t < ts[0]) or np.any(t > ts[-1]):
        raise ValueError(
            f"temperature iterate leaves sampled span [{ts[0]}, {ts[-1]}]: "
            f"[{t.min()}, {t.max()}]"
        )
    return np.interp(t, ts, ks)


def source_nodal(source, positions, extent):
    """Nodal volumetric source, W/m^3, with a documented op order."""
    kind = source["kind"]
    if kind == "none":
        return np.zeros(positions.shape[0], dtype=np.float64)
    if kind == "poly-xy":
        # Exact op order (both codes): q0 * (x / ex) * (1.0 - y / ey)
        q0 = float(source["q0"])
        ex, ey = extent[0], extent[1]
        out = np.empty(positions.shape[0], dtype=np.float64)
        for n in range(positions.shape[0]):
            x = positions[n, 0]
            y = positions[n, 1]
            out[n] = q0 * (x / ex) * (1.0 - y / ey)
        return out
    raise ValueError(f"unknown source kind {kind}")


def solve_case(case):
    counts = [int(c) for c in case["mesh"]["counts"]]
    extent = [float(e) for e in case["mesh"]["extent"]]
    positions, tets = kuhn_box_grid(counts, extent)
    mhash = mesh_blake3(counts, positions, tets)

    mesh = MeshTet(positions.T.copy(), tets.T.copy())
    basis = Basis(mesh, ElementTetP1())
    nverts = positions.shape[0]

    material = case["material"]
    nonlinear = material["kind"] == "linear_kt"
    if nonlinear and material.get("evaluation") != "element-mean-temperature":
        raise ValueError("k(T) cases must declare element-mean-temperature evaluation")

    q_nodal = source_nodal(case.get("source", {"kind": "none"}), positions, extent)

    # Boundary partition.
    dirichlet_nodes = {}
    robin_faces = []  # (axis, target, h, t_inf)
    for bc in case["bcs"]:
        axis, target = face_value(bc["face"], extent)
        if bc["kind"] == "dirichlet":
            for n in nodes_on_face(positions, axis, target):
                dirichlet_nodes[int(n)] = float(bc["t_k"])
        elif bc["kind"] == "robin":
            robin_faces.append((axis, target, float(bc["h"]), float(bc["t_inf_k"])))
        elif bc["kind"] == "adiabatic":
            pass
        else:
            raise ValueError(f"unknown bc kind {bc['kind']}")

    fixed = np.asarray(sorted(dirichlet_nodes), dtype=np.int64)
    free = np.asarray(
        [n for n in range(nverts) if n not in dirichlet_nodes], dtype=np.int64
    )
    xp = np.zeros(nverts, dtype=np.float64)
    for n, value in dirichlet_nodes.items():
        xp[n] = value

    @LinearForm
    def source_form(v, w):
        return w["q"] * v

    b_source = source_form.assemble(basis, q=basis.interpolate(q_nodal))

    # Robin blocks: consistent facet mass and load per declared face set.
    robin_A = scipy.sparse.csr_matrix((nverts, nverts))
    robin_b = np.zeros(nverts, dtype=np.float64)
    for axis, target, h, t_inf in robin_faces:
        facets = mesh.facets_satisfying(lambda x: np.abs(x[axis] - target) < FACE_BAND)
        expected = 2 * np.prod([counts[a] for a in range(3) if a != axis])
        if facets.shape[0] != int(expected):
            raise ValueError(
                f"face selection ({axis},{target}) matched {facets.shape[0]} "
                f"facets, expected {expected}"
            )
        fbasis = FacetBasis(mesh, ElementTetP1(), facets=facets)

        @BilinearForm
        def robin_mass(u, v, w):
            return u * v

        @LinearForm
        def robin_load(v, w):
            return v

        robin_A = robin_A + h * robin_mass.assemble(fbasis)
        robin_b = robin_b + (h * t_inf) * robin_load.assemble(fbasis)

    def assemble_stiffness(temperature):
        if nonlinear:
            tbar = temperature[tets].mean(axis=1)
            kbar = kt_eval(material["knots"], tbar)
            nqp = basis.X.shape[1]
            kk = np.broadcast_to(kbar[:, None], (kbar.shape[0], nqp))

            @BilinearForm
            def conduction(u, v, w):
                return w["kk"] * dot(grad(u), grad(v))

            return conduction.assemble(basis, kk=kk)
        ktensor = conductivity_tensor(material)

        @BilinearForm
        def conduction(u, v, w):
            gu = grad(u)
            gv = grad(v)
            out = 0.0
            for i in range(3):
                for j in range(3):
                    if ktensor[i, j] != 0.0:
                        out = out + ktensor[i, j] * gu[j] * gv[i]
            return out

        return conduction.assemble(basis)

    def linear_solve(temperature):
        A = (assemble_stiffness(temperature) + robin_A).tocsr()
        b = b_source + robin_b
        A_ff = A[free][:, free]
        A_fp = A[free][:, fixed]
        rhs = b[free] - A_fp @ xp[fixed]
        x_f = scipy.sparse.linalg.spsolve(A_ff.tocsc(), rhs)
        full = xp.copy()
        full[free] = x_f
        return full

    picard = case.get("picard", {"tol": 1e-13, "max_iterations": 60})
    t = xp.copy()
    if fixed.shape[0] > 0:
        t[free] = np.mean([dirichlet_nodes[int(n)] for n in fixed])
    else:
        # Pure-Robin case: start from the mean of the declared T_inf values.
        t[:] = np.mean([r[3] for r in robin_faces])

    iterations = 0
    if nonlinear:
        for iterations in range(1, int(picard["max_iterations"]) + 1):
            t_new = linear_solve(t)
            step = np.max(np.abs(t_new - t))
            t = t_new
            if step <= float(picard["tol"]) * max(1.0, np.max(np.abs(t))):
                break
        else:
            raise ValueError("Picard iteration budget exhausted")
    else:
        t = linear_solve(t)
        iterations = 1

    # Probe extraction by GRID INDEX (no floating-point matching).
    py, pz = counts[1] + 1, counts[2] + 1
    probes = []
    for gi, gj, gk in case["probes"]:
        n = int(gi) * py * pz + int(gj) * pz + int(gk)
        probes.append(
            (
                (int(gi), int(gj), int(gk)),
                (positions[n, 0], positions[n, 1], positions[n, 2]),
                t[n],
            )
        )
    return {
        "mesh_blake3": mhash,
        "picard_iterations": iterations,
        "probes": probes,
        "temperature_range": (float(t.min()), float(t.max())),
    }


def spec_echo_rows(case):
    """Flat numeric echo of every load-bearing spec value.

    The Rust catalog compares these against its own typed constants
    bit-exactly, binding the frozen reference to the case it solved.
    """
    rows = []
    counts = case["mesh"]["counts"]
    extent = case["mesh"]["extent"]
    for axis, name in enumerate(["x", "y", "z"]):
        rows.append((f"mesh.counts.{name}", float(counts[axis])))
        rows.append((f"mesh.extent.{name}", float(extent[axis])))
    material = case["material"]
    if material["kind"] == "isotropic":
        rows.append(("material.isotropic.k", float(material["k"])))
    elif material["kind"] == "tensor":
        k = material["k"]
        for i in range(3):
            for j in range(3):
                rows.append((f"material.tensor.k{i}{j}", float(k[i][j])))
    elif material["kind"] == "linear_kt":
        for n, (tk, kk) in enumerate(material["knots"]):
            rows.append((f"material.linear_kt.knot{n}.t", float(tk)))
            rows.append((f"material.linear_kt.knot{n}.k", float(kk)))
        picard = case.get("picard", {"tol": 1e-13, "max_iterations": 60})
        rows.append(("picard.tol", float(picard["tol"])))
        rows.append(("picard.max_iterations", float(picard["max_iterations"])))
    source = case.get("source", {"kind": "none"})
    if source["kind"] == "poly-xy":
        rows.append(("source.poly_xy.q0", float(source["q0"])))
    for bc in case["bcs"]:
        prefix = f"bc.{bc['name']}"
        axis, target = face_value(bc["face"], case["mesh"]["extent"])
        rows.append((f"{prefix}.axis", float(axis)))
        rows.append((f"{prefix}.target", float(target)))
        if bc["kind"] == "dirichlet":
            rows.append((f"{prefix}.dirichlet.t", float(bc["t_k"])))
        elif bc["kind"] == "robin":
            rows.append((f"{prefix}.robin.h", float(bc["h"])))
            rows.append((f"{prefix}.robin.t_inf", float(bc["t_inf_k"])))
    return rows


def self_check():
    """Known-answer probes proving the harness plumbing.

    1. Dirichlet-Dirichlet slab: exact linear profile is IN the P1
       space, so nodal agreement must be round-off level.
    2. Dirichlet-Robin slab: exact solution is again linear; checks the
       Robin matrix/load plumbing against the closed form.
    """
    case = {
        "id": "self-check-slab-dd",
        "mesh": {"counts": [6, 3, 3], "extent": [0.2, 1.0, 1.0]},
        "material": {"kind": "isotropic", "k": 20.0},
        "source": {"kind": "none"},
        "bcs": [
            {"name": "hot", "face": ["x", 0.0], "kind": "dirichlet", "t_k": 340.0},
            {"name": "cold", "face": ["x", "max"], "kind": "dirichlet", "t_k": 300.0},
        ],
        "probes": [[3, 1, 1]],
    }
    result = solve_case(case)
    (_, (x, _, _), t) = result["probes"][0]
    exact = 340.0 + (300.0 - 340.0) * (x / 0.2)
    if abs(t - exact) > 1e-9:
        raise AssertionError(f"self-check DD slab: |{t} - {exact}| > 1e-9")

    k, length, h, t_hot, t_inf = 20.0, 0.2, 50.0, 340.0, 300.0
    case = {
        "id": "self-check-slab-dr",
        "mesh": {"counts": [6, 3, 3], "extent": [length, 1.0, 1.0]},
        "material": {"kind": "isotropic", "k": k},
        "source": {"kind": "none"},
        "bcs": [
            {"name": "hot", "face": ["x", 0.0], "kind": "dirichlet", "t_k": t_hot},
            {
                "name": "cooled",
                "face": ["x", "max"],
                "kind": "robin",
                "h": h,
                "t_inf_k": t_inf,
            },
        ],
        "probes": [[6, 1, 1]],
    }
    result = solve_case(case)
    (_, _, t_l) = result["probes"][0]
    # Exact wall temperature: T(L) = T_inf + (T_hot - T_inf) / (1 + hL/k).
    exact_l = t_inf + (t_hot - t_inf) / (1.0 + h * length / k)
    if abs(t_l - exact_l) > 1e-9:
        raise AssertionError(f"self-check DR slab: |{t_l} - {exact_l}| > 1e-9")
    return "pass"


def fmt(value):
    """Shortest round-trip decimal for an f64 (Python repr)."""
    return repr(float(value))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    check = self_check()

    lines = [
        "# frankensim Level-B cross-code reference manifest v1",
        "# generated by tools/vvref/solve_skfem.py (pinned by tools/vvref/uv.lock)",
        "# columns: kind\tcase_id\tfield...\tvalue",
    ]
    for path in args.cases:
        with open(path, "rb") as f:
            deck_bytes = f.read()
        case = json.loads(deck_bytes)
        result = solve_case(case)
        cid = case["id"]
        meta = [
            ("deck_blake3", blake3(deck_bytes).hexdigest()),
            ("external_code", f"scikit-fem {skfem.__version__}"),
            ("linear_solver", f"scipy {scipy.__version__} spsolve/SuperLU"),
            ("mesh_blake3", result["mesh_blake3"]),
            ("numpy", np.__version__),
            ("picard_iterations", str(result["picard_iterations"])),
            ("platform", platform.machine()),
            ("python", platform.python_version()),
            ("self_check", check),
            ("t_max_k", fmt(result["temperature_range"][1])),
            ("t_min_k", fmt(result["temperature_range"][0])),
        ]
        for key, value in meta:
            lines.append(f"case_meta\t{cid}\t{key}\t{value}")
        for key, value in spec_echo_rows(case):
            lines.append(f"spec_echo\t{cid}\t{key}\t{fmt(value)}")
        for index, ((gi, gj, gk), (x, y, z), t) in enumerate(result["probes"]):
            lines.append(
                f"probe\t{cid}\t{index}\t{gi}\t{gj}\t{gk}\t"
                f"{fmt(x)}\t{fmt(y)}\t{fmt(z)}\t{fmt(t)}"
            )
    payload = "\n".join(lines) + "\n"
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(payload)
    print(f"wrote {args.out} ({len(payload)} bytes), self_check={check}")


if __name__ == "__main__":
    sys.exit(main())
