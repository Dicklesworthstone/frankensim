# Cooling across independently meshed planar contacts

```bash
frankensim --json cooling-network \
  examples/cooling-network/nonmatching-contact-hotspot.json
```

This example has 32 vertices and 42 tetrahedra: six cells in one solid and
36 in the other. The common 0.01 m2 contact has two triangles on side A and
twelve on side B. Neither volume mesh is remeshed or welded. A 1 W localized
P1 heat source, two materials and fan-driven bypass/mixing complete the actual
cooling problem. Hardware, material and resistance inputs are illustrative,
not measured or experimentally validated data.

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

The new `fs_conduction::interface::NonmatchingSurface` and
`ThermalInterfaces::with_nonmatching` construct the common refinement of the
selected trace triangles. On each overlap, the contribution is

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

## Resource and scope boundaries

Pair tests include both same-side disjointness checks and cross-side overlaps.
The complete pair count is admitted before clipping. Overlap output has its own
explicit cap. The CLI caps these at 1,000,000 pair tests and 100,000 overlap
triangles per interface. Core binding accepts a cancellation context and polls
it. CLI geometry binding occurs during existing request parsing, BEFORE the
thermal solve wall watchdog; its pair/output caps bound this preprocessing,
not a newly claimed parse-time deadline. Solve assembly then uses the existing
solve cancellation, residual, energy and iteration budgets.

`mesh_convergence` currently refuses this new side-declaration form before its
first physical solve: its transfer code only knows matching face pairs.
Explicit nonmatching meshes may be solved separately. No `.fsim` nonmatching
contact lowering, curved contact, pressure-dependent resistance, general mortar
constraint, contact detection or automatic volume-mesh repair is added.
Exact coincident matching triangles still require ownership everywhere; the
new path does NOT search the whole mesh for undeclared nonmatching partners.
This narrowly scoped explicit path extends the older matching-only summaries.

## Focused verification

```bash
cargo test -p fs-conduction --test nonmatching_contact
cargo test -p fs-cli --test cooling_nonmatching_contact
cargo test -p fs-cli --bin frankensim uq_command::model::contact
```

Thirteen new Rust tests include six numerical-library regressions, six
actual-command regressions and one target-mutation regression. They cover
unequal meshes and analytic slab fields, mean-zero contact modes, exact-pair
delegation, signed heat under side reversal, geometry/budget refusals,
nonlinear/radiative steady gradients, repeated radiating storage adjoints,
unchanged forward fields, actual uncertain samples and exact chunked replay.
These Rust tests have NOT been executed in the authoring environment, which
has no Rust toolchain. Compilation and actual CLI behavior remain unverified.

Independent Python calculations use exact-rational intersection enumeration
and analytic polynomial moments as an oracle for the clipping/quadrature
mathematical mirror. Five slab/scale cases had at most 6.82e-13 K field error
against the independent analytic solution. The largest scaled matrix
discrepancy against the exact-moment oracle was 6.66e-16. Four complete cooling
finite-difference comparisons, including nonlinear material and radiation
feedback, differed by at most 2.70e-9 K per unit log-resistance change.

The example's independent reference peak is 301.810057 K; its contact carries
0.359936 W and d(peak)/dln(R'') is 0.168944 K. These calculations are NOT
executions of FrankenSim, mesh-error certificates or physical validation.
