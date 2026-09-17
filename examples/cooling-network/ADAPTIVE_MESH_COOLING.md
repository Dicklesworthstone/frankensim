# Goal-weighted local refinement for coupled cooling

```bash
frankensim --json cooling-network examples/cooling-network/adaptive-contact-hotspot.json
```

This is a real solve/adjoint/mark/refine loop in the existing steady cooling
command. It retains the original coupled conduction, contact, nonlinear-material,
air-mixing and energy checks. The existing uniform mesh study remains the default.

## Request

Use `objective.gradient=false` and declare the marker explicitly:

```json
"mesh_convergence": {
  "strategy": "goal-recovery",
  "marking_fraction": 0.5,
  "max_refinements": 12,
  "consecutive_passes": 2,
  "temperature_tolerance_k": 0.02,
  "max_vertices": 20000,
  "max_tetrahedra": 100000
}
```

`marking_fraction` is in (0,1] and is required for `goal-recovery`. The local
strategy permits up to 32 refinement rounds; the uniform strategy retains its
six-round cap. The round limit includes global confirmation refinements.
Other mesh-study controls retain their existing meaning. Transient and nested
design studies remain excluded. This does not wire the native `.fsim` fidelity
controller or complete the separate rigorous DWR/error-bound work.

## Marking and conformity

Each solved mesh supplies its temperature field and the ACTUAL total coupled
nodal-load adjoint for the selected objective. The adjoint includes the existing
mixed-air feedback and material Jacobian. Vertex/material patches recover the
piecewise-constant P1 gradients. A cell's volume-weighted product of primal and
adjoint recovery defects supplies a refinement SCORE. Global normalization avoids
arbitrary changes from multiplying the adjoint or adding a uniform temperature
offset. Recovery does not average across material labels or independent contact
trace IDs.

Largest scores are marked until their sum reaches the requested fraction of the
complete score. Equal scores break ties by cell index. Zero total score triggers
a uniform probe; it never means zero discretization error or permission to stop.
These are heuristic recovery scores, NOT dual-weighted residual estimates,
Kelvin error bars, equilibrated fluxes, or outward-rounded certificates.

Each marked cell nominates its longest edge. All tetrahedra incident to that
edge are bisected, including unmarked neighbours required for conformity. Contact
edge correspondences close the selection on both traces. Common geometric edge
order and geometric child-face matching preserve the paired triangulations even
when opposite sides use different vertex numbering/order. Midpoints are not
welded across contact gaps; resistance is not rescaled.

P1 heating fields are prolonged WITHOUT re-normalizing or shrinking their
original support. Every child inherits its actual parent material, not an assumed
eight-child index. Cooling face names and constitutive declarations are preserved.
The resulting request is re-admitted by the ordinary cooling parser. Orientation,
source-watt, geometry and solver failures are refusals; there is no implicit
repair or clipping. Local bisection carries no universal shape-regularity or
Delaunay guarantee. Existing downstream geometry admission still applies.

## Global check before success

Two or more consecutive small objective changes are required. When local
refinement reaches that streak, the next round refines the ENTIRE mesh with the
existing uniform red-refinement kernel. An adaptive result can be returned only
when a complete uniform comparison also meets the original temperature tolerance.
If it fails, the streak resets and local refinement continues. If its mesh or
round count exceeds a budget, the command refuses rather than accepting the local
streak. This safeguard can cost substantially more than uniform-only refinement.

Agreement still does not prove a continuum error bound or physical accuracy.
A common bias or unresolved asymptotic regime can survive both local and uniform
comparisons. Maxima concern the selected discrete field; ties retain the existing
active-vertex rule. No globally optimal mesh, universal cell reduction, physical
compliance or measured hardware validation is inferred.

## Results and budgets

The result keeps `status=successive-mesh-tolerance-met`, with
`method=goal-recovery-edge-bisection` and `global_confirmation=true`.
Each history row records the refinement used to arrive at that mesh and the
marker score sum, marked count, captured fraction and adjoint sweeps. The count
is the nominated cells, not the additional cells required by edge/contact closure.
`score_is_error_bound` is always false.

`total_solid_solves` counts primal coupling solves, as before;
`total_adjoint_sweeps` separately reports the marker's additional interface work,
not every inner Krylov iteration. Existing derivative and linear limits apply on
each solved mesh. The one original wall deadline covers the entire loop,
including adjoints and refinement. The output-count limits are not total allocator
memory guarantees.

The final `resolved_request` owns the published temperature vector and includes
its actual refined nodal source, material assignment and contact faces. It keeps
`gradient=false`, has no recursive mesh-study option, and can be independently
solved by the ordinary command. Raw adjoint fields used only for marking are not
misrepresented as requested output gradients. Unrequested uniform studies do no
new marker-adjoint work.

## Focused tests and independent reference

```bash
cargo test -p fs-mesh --test marked_tet
cargo test -p fs-cli --test cooling_adaptive_mesh
cargo test -p fs-cli --test cooling_mesh_convergence
cargo test -p fs-cli --bin frankensim mesh_convergence
```

Eleven new Rust tests cover edge-star conformity, remote-cell preservation,
signed P1/source-volume conservation, material parentage, shuffled contact traces,
mark-order replay, actual nonlinear cooling/refined-request replay, global
confirmation, and cancellation/resource/adjoint failures. These tests were NOT
executed in the authoring environment, which had no Rust toolchain.

An independent sparse NumPy/SciPy FEM calculation, with analytic air elimination
and a direct transpose solve, exercised 72 random-mark/contact/source-transfer
cases. Its nonlinear one-watt contact example used 12, 24, 36 and 288 cells,
with peaks 301.955112, 301.961392, 301.974972 and 301.978296 K. The final global
comparison changed the peak by 0.003324 K. A separate 6144-cell uniform reference
was 301.971581 K; this observed discrepancy is not a rigorous bound.

A second, non-contact hotspot exposed why the global check matters: two apparently
settled local sequences were followed by global changes of 0.108092 and 0.065080 K,
both above the 0.02 K tolerance. It needed 25792 final cells before the global
check passed. There is therefore NO universal speedup/cell-saving claim. These
are independent numerical reference results, not executions or benchmarks of
the committed Rust implementation. The example inputs are illustrative.
