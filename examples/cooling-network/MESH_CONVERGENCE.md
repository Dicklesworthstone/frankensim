# Mesh-refinement studies for coupled cooling

The experimental `cooling-network` command can execute the same declared steady
cooling problem on successively refined tetrahedral meshes. This is a real
solve/refine/solve loop, not a mesh-quality label or a single-solve claim.

```bash
frankensim --json cooling-network \
  examples/cooling-network/mesh-convergence-contact-hotspot.json
```

The example combines a one-watt localized P1 source, two temperature-dependent
conductivity laws, a finite-resistance matching contact, and fan-driven bypass
mixing. All physical inputs are illustrative declarations, not measured data.

## Explicit stopping policy

Add this object at the request root:

```json
"mesh_convergence": {
  "max_refinements": 3,
  "consecutive_passes": 2,
  "temperature_tolerance_k": 0.02,
  "max_vertices": 20000,
  "max_tetrahedra": 100000
}
```

The base mesh is solved first. Each refinement splits every tetrahedron into
eight children and every selected triangle into four, using the existing
`fs-mesh` red-refinement kernel. The next mesh is re-admitted by the ordinary
cooling parser and solved by the same coupled conduction/air producer.

A comparison passes when the absolute change in the selected objective from the
previous solved mesh is at most `temperature_tolerance_k`. Success requires the
declared number of CONSECUTIVE passing comparisons. At least two comparisons
are mandatory, so at least three meshes must actually be solved. A failed
comparison resets the streak. `max_refinements` excludes the base solve and
must be between two and six; `consecutive_passes` is between two and that cap.

The result's `mesh_convergence.status` is `successive-mesh-tolerance-met`, NOT a
certified continuum convergence verdict. Its history records every solved
mesh's vertex/cell counts, objective, change, source watts, solid-solve count and
pass streak. `total_solid_solves` sums real coupling solves across the ladder.

This is an observed UNIFORM h-ladder, not goal-oriented/DWR marking, an adaptive
mesh-quality improver, Richardson extrapolation, or a guaranteed error bound.
No monotonicity or asymptotic rate is assumed. Agreement can occur before the
asymptotic regime or miss a common bias. A surface-mean agreement does not bound
a maximum elsewhere. A small difference does not produce a physical compliance
margin or change the result's nominal-estimate authority.

## Preserve the actual problem under refinement

The source transfer is essential. A component's original `vertices` select a
P1 basis-function footprint, not a mesh-independent sharp subvolume. Selecting
those SAME vertex IDs on a finer mesh shrinks the footprint. Renormalizing it
to the same total watts still changes the problem and can strongly change its
peak temperature.

The study instead resolves `component_power` ONCE on the base mesh, then
prolongates that piecewise-linear volumetric density: original values are
retained and each edge midpoint receives the mean of its endpoint densities.
No renormalization or new component-footprint interpretation occurs. The
physical P1 source and its integral are preserved in exact arithmetic. Every
solved rung also checks its assembled watts against the base result under the
original `heat_w` tolerance. A source-total change is a refusal, not convergence.

This field is represented by the new explicit `solid.nodal_source_w_m3` mode:
exactly one finite density per node, mutually exclusive with `source_w_m3` and
`component_power`. The values are W/m3, not nodal watts. Signed source densities
remain allowed just as signed uniform sources were. Ordinary uniform sources
remain uniform. Results disclose the actual nodal density when this mode is used.

Each child inherits its parent's material assignment and full constitutive law,
including anisotropic tensors or bounded scalar k(T). Cooling faces subdivide
while retaining their surface name and coefficient/correlation declaration.
Matching contact triangles subdivide on BOTH sides under their geometric vertex
correspondence. Coincident but separately numbered traces are never welded,
and contact resistance is not rescaled. This refinement path requires exactly
coincident corresponding coordinates; a tolerance-only match is refused rather
than silently moving geometry.

Mean-wall and wall/whole-solid maximum selectors remain the same geometric
functionals; refined maxima include newly introduced vertices. `max_vertices`
continues to observe the specified ORIGINAL vertices, which keep their IDs. It
does not silently expand a point-set query into a whole-solid maximum.

## Retain a usable final mesh and source

`mesh_convergence.resolved_request` is the complete final refined input without
the mesh-study option. Its node numbering owns the published temperature vector
and active-vertex ID. Save it as JSON and pass it to `cooling-network` to evaluate
the same refined problem independently. It retains the actual P1 density,
material assignments, refined cooling/contact faces, and unchanged hydraulics.
For example, with `jq` installed:

```bash
frankensim --json cooling-network \
  examples/cooling-network/mesh-convergence-contact-hotspot.json > study.json
jq '.mesh_convergence.resolved_request' study.json > resolved-cooling.json
frankensim --json cooling-network resolved-cooling.json
```

Exact same-profile replay of temperatures and objective is asserted by the
committed actual-command regression; it has not been executed in the authoring
environment. The resolved field source no longer has named component controls;
those base declarations must not be reattached as a new fine-mesh footprint.

## Resource and feature boundaries

`max_vertices` and `max_tetrahedra` bound each complete refined mesh, at most
20,000 vertices and 100,000 cells. Output counts are checked before splitting;
they are not a total allocator/workspace guarantee. Refined inputs also obey
the original 16 MiB input cap. Cancellation is polled during admission/transfer,
around the bounded split-kernel call, and through the existing numerical solves.
The original wall deadline covers the complete study; each rung retains the
existing graph, coupling and linear-iteration budgets without increasing them.

Count limits, refinement exhaustion or wall exhaustion return budget exit 6,
not an accepted temperature result. Solver, material-range, contact and source
failures propagate. No partial mesh-study field or resumable checkpoint is
published on failure. This does not add a second checkpoint framework.

The first product path requires steady input, `objective.gradient=false`, and
no nested effective-h/fan design or transient schedule. Ordinary solves without
`mesh_convergence` keep their existing path. This is not native `.fsim` fidelity
adaptation and does not close the remaining DWR/guaranteed-bound plan obligations.

## Numerical references and tests

Independent sparse P1 FEM calculations eliminate the two air-reference equations
analytically rather than repeating the production staggered loop. For the
nonlinear contact example, they give:

| Mesh | Vertices | Tetrahedra | Peak, K | Successive change, K |
|---|---:|---:|---:|---:|
| Base | 16 | 12 | 301.95511185 | — |
| Refined once | 54 | 96 | 301.93619665 | 0.01891520 |
| Refined twice | 250 | 768 | 301.95254839 | 0.01635174 |

The peak is non-monotone. Both changes meet the example's 0.02 K observation
criterion, but that is NOT a 0.02 K continuum error certificate. Rebinding the
original hot vertex on the final mesh while preserving one total watt instead
gives 304.22214356 K: a different spatial source, not a refinement estimate.

A separate uniformly heated linear slab has an independently derived continuum
maximum of 301.21693753 K. Its reference P1 peak error decreases from 0.00636474 K
on 12 cells to 0.00077914 K on 6,144 cells, again without monotone peak values.
There is no universal rate or efficiency claim from this fixture.

```bash
cargo test -p fs-mesh --test tet_refinement
cargo test -p fs-cli --bin frankensim network_command::mesh_convergence
cargo test -p fs-cli --test cooling_mesh_convergence
```

Fifteen added Rust tests include four actual-command regressions. They cover
source interpolation/integration, independent coincident traces, material and
contact inheritance, complete refined-request replay, observable selection,
output limits, cancellation and refusal. These Rust tests were NOT executed in
the authoring environment, which has no Rust toolchain. The independent Python
reference ran seven mesh/source cases with seven deliberately altered-footprint
controls, four analytic-slab comparisons, and 36 signed-source transfer checks.
Those calculations are not Rust compilation, CLI execution or physical validation.
