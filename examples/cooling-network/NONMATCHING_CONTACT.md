# Cooling across independently meshed planar contacts

```bash
frankensim --json cooling-network \
  examples/cooling-network/nonmatching-contact-hotspot.json
```

This example has 32 vertices and 42 tetrahedra: six cells in one solid and
36 in the other. The common 0.01 m2 contact has two triangles on side A and
twelve on side B. Neither volume mesh is remeshed or welded for contact binding.
A 1 W localized P1 heat source, two materials and fan-driven bypass/mixing
complete the actual cooling problem. Hardware, material and resistance inputs
are illustrative, not measured or experimentally validated data.

## Declare both complete traces

Use `nonmatching` INSTEAD OF `face_pairs` in one `solid.contacts` entry. Keep
its existing name, source, side-material labels and `resistance_m2_k_w`:

```json
"nonmatching": {
  "side_a_faces": [[1,3,7],[1,5,7]],
  "side_b_faces": [[8,10,16],[8,14,16],[10,12,18],[10,16,18],
                   [14,16,22],[14,20,22],[16,18,24],[16,22,24],
                   [20,22,28],[20,26,28],[22,24,30],[22,28,30]],
  "plane_tolerance_m": 1e-12,
  "coverage_relative_tolerance": 1e-10,
  "max_pair_tests": 10000,
  "max_overlap_triangles": 10000
}
```

These indices refer to the complete example mesh, not arbitrary geometry.
The faces must be exterior trace triangles, separately numbered on the two
solids, with opposing outward normals. The complete declared surface is planar;
within-side overlapping facets, incomplete coverage, duplicate ownership and
external boundary conditions on contact faces refuse. Declare separate contacts
for separate planes. A face cannot belong to two contacts or to contact and air.

Plane tolerance is a numerical admission allowance, not a gap constitutive
law: vertices are projected to the admitted reference plane. It must be finite,
nonnegative and at most one millionth of the trace span. Coverage tolerance
must be in (0,1e-6] and applies to EACH face's overlap area, not only the total.
Neither test is an exact/certified geometry predicate. A larger tolerance is
not a license to model a physical open gap or misaligned parts.

The resistance remains one strictly positive, constant area-specific value
R'' in m2 K/W. It is neither a penalty parameter nor inferred perfect contact.
Side-material names and the inline material card retain their existing
caller-declaration authority; this feature does not upgrade those labels to
independent material evidence.

## What the operator does

The `fs_conduction::interface::NonmatchingSurface` and
`ThermalInterfaces::with_nonmatching` producers construct the common refinement
of the selected trace triangles. On each overlap, the contribution is

```
integral (T_A - T_B) * (v_A - v_B) / R'' dA.
```

Both original P1 bases are evaluated at the SAME quadrature points. Three
positive degree-two weights integrate their product on each overlap triangle.
This is not nearest-node matching, averaging each surface to one temperature,
or replacing the layer with a convection coefficient. The real-arithmetic
bilinear form is symmetric and nonnegative, annihilates common constants and
transfers equal/opposite heat. Floating-point assembly remains subject to the
existing residual and energy gates, not an exact-conservation certificate.

An exactly matching subpatch is delegated to the unchanged original matching
producer once; it is not double counted. Matching-only `face_pairs` inputs
retain that original numerical path. Named flux reports merge any exact and
nonmatching pieces into one physical interface. For a nonmatching report,
`face_pairs` is null and `discretization` is `planar-common-refinement-P1`.
It does not invent a one-to-one pairing where none exists.

## Solves, sensitivities and uncertainty

The operator enters the existing steady, temperature-dependent, transient and
radiative conduction assemblers. The same matrix enters their tangents and
transposes, so earlier solid-storage and coupled-air derivatives include its
thermal influence. The example requests a steady contact-resistance derivative.
For p=ln(R''), dK/dp=-K, and the total coupled load multiplier contracts as
`lambda^T K_contact T`. This uses the full overlap integrals, not the product
of two mean jumps. Geometry and the overlap topology stay fixed.

This does not add a transient contact-resistance control: existing transient
power/fan/history/capacity adjoints can use the contact operator, while their
resistance remains fixed. It does not add shape or contact-search derivatives.

For UQ, set the base request's `objective.gradient=false` and use the existing
parameter target, for example:

```json
{
  "target": {"kind":"contact-resistance","contact":"bondline"},
  "distribution": {"kind":"uniform","lo":0.005,"hi":0.025}
}
```

Each draw changes the real interface law before the cooling child runs. Both
meshes, source footprint, contact sides and geometry policies are unchanged.
A transient holds one draw throughout its trajectory. Existing checkpoint,
exact-ordinal resume and confidence-decision semantics apply. Invalid physical
samples refuse; they are not clipped, filtered or replaced by redraws.

## Uniform and locally adaptive mesh studies

```bash
frankensim --json cooling-network \
  examples/cooling-network/adaptive-nonmatching-radiative-hotspot.json
```

Steady `mesh_convergence` now supports these independent side declarations,
including nonlinear materials, radiation and total-adjoint goal-recovery marking.
Each contact side follows its OWN volume refinement. A mark in one solid does
not create linked edges in its nonmatching partner; the partner may remain
coarse. Matching-only contacts still synchronize their paired edges as before.
The same contact binder reconstructs overlap quadrature on every refined mesh,
including any subpatch that changes between exact matching and nonmatching.
Cached quadrature points are not prolonged or reused on a changed mesh.

The original source field is prolonged in P1 without shrinking its footprint or
renormalizing it. Materials inherit by parent cell. Contact resistance, side
labels, patch areas, radiation declarations and geometry tolerances are not
rescaled. Independently numbered nodes remain separate even when their physical
coordinates coincide. `resolved_request` retains the final side faces, mesh and
source for ordinary `cooling-network` replay, with no recursive mesh study.

The example requests two consecutive temperature changes below 0.05 K and a
complete uniform-refinement confirmation after local agreement. It explicitly
allocates 1,000,000 pair tests and 100,000 overlap triangles per binding. These
are example INPUT budgets, not automatic solver increases. A refined interface
with n total side faces requires n*(n-1)/2 pair tests, including same-side
checks. That count is admitted before rebinding; exhaustion returns a mesh-budget
refusal rather than accepting an unconfirmed result. Overlap-output, derivative,
mesh-size and original wall budgets also retain their existing refusals.

Both strategies use the full same-model cooling producer on every mesh. Internal
marker derivatives include contact and radiative feedback but are not published
as user-requested sensitivities when `objective.gradient=false`. Repeated mesh
agreement and the global probe are observations, not a continuum error bound,
point-maximum certificate, monotonicity guarantee or proof of physical accuracy.

## Resource and scope boundaries

Pair tests include both same-side disjointness checks and cross-side overlaps.
The complete pair count is admitted before clipping. Overlap output has its own
explicit cap. The CLI caps these at 1,000,000 pair tests and 100,000 overlap
triangles per interface. Core binding accepts a cancellation context and polls
it. Initial CLI geometry binding occurs during request parsing, BEFORE the
thermal solve wall watchdog. Geometry rebinding during mesh studies also uses
the parser's bounded preprocessing context rather than an independently timed
clipping deadline. Its pair/output caps bound that work; thermal solves and
refinement polls still use the original study deadline. No new per-level time
allowance is created.

No `.fsim` nonmatching contact lowering, curved contact, pressure-dependent
resistance, general mortar constraint, contact detection or automatic volume-mesh
repair is added. Exact coincident matching triangles still require ownership
everywhere; this explicit path does NOT search the whole mesh for undeclared
nonmatching partners. Transient mesh studies and space-time adaptation remain
unsupported.

## Focused verification

```bash
cargo test -p fs-conduction --test nonmatching_contact
cargo test -p fs-cli --test cooling_nonmatching_contact
cargo test -p fs-cli --bin frankensim contact_transfer
cargo test -p fs-cli --test cooling_adaptive_mesh
cargo test -p fs-cli --test cooling_radiation_mesh
```

The original thirteen nonmatching-contact Rust regressions remain. Six additional
tests cover one-sided local refinement, nonuniform contact bilinear forms,
source conservation, complete uniform/adaptive radiative command runs, exact
refined-field replay, face-order determinism, global confirmation and unchanged
contact work caps. An older mesh-refusal test now checks actual pair-budget
exhaustion rather than missing functionality. These Rust tests have NOT been
executed in the authoring environment, which lacks Rust/RCH. Compilation and
actual command behavior remain unverified.

Independent Python calculations use exact-rational contact intersections and
analytic polynomial moments with sparse volume FEM and direct transpose solves.
Sixteen signed-field refinement checks preserved the contact bilinear form to
within 1.67e-16. Refining only the first volume produced 2, 4 and 8 contact
triangles while the other side stayed at 12; its zero-mean nonuniform jump still
integrated to the analytic value 1/12, not zero.

For the new nonlinear radiative example, the independent local ladder used
42, 48, 60 and 480 cells with peaks 301.699226, 301.691196, 301.709776 and
301.690594 K. The last mesh is the complete global probe; its change was
0.019182 K. A separate uniform ladder used 42, 336 and 2688 cells, ending at
301.686480 K. These are observed comparisons, NOT Rust execution, certified
error bounds or measured speedups. Removing radiation on the final local mesh
raised the reference peak to 301.792381 K. The retained source is one watt,
with separate air and radiation accounting throughout.
