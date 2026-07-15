# fs-rep-nurbs — CONTRACT

Rational B-spline charts (plan §7.2): bounded EXACT spline algebra, trimmed
patches with certified classification, measured closest-point brackets,
and the HONEST Boolean position.

Ambition tags: exact algebra / trims / closest-point [F per bead label;
what is shipped is tested to [S] discipline]; direct B-rep Booleans
deliberately NOT shipped (see the Boolean position).

## Purpose and layer

Layer **L2** (MORPH). Runtime deps: `std`, fs-evidence, fs-exec, fs-geom,
fs-ivl, fs-math. Consumers:
fs-iga (geometry basis = analysis basis), fs-render NURBS tracing
(shares the clipping/Newton machinery), the NURBS↔SDF converter beads
(wqd.11/12).

## Public types and semantics

- `Rat` — exact i128 rational scalar (gcd-reduced, overflow-CHECKED:
  leaving the exactness domain is never silent wraparound). Fallible exact
  helpers return `NurbsError::Exactness`; the present operator traits cannot
  transport it and retain a named panic at that boundary.
- `Scalar` — the field abstraction; the SAME generic basis/curve/surface
  code runs at `f64` (fast) and `Rat` (exact). The conformance suite
  also instantiates it at a test dual, so derivative checks flow through
  the identical code path. Scalar domains explicitly define weight
  admissibility; the f64 path rejects subnormal weights whose multiplication
  by a basis value can underflow an otherwise positive denominator.
- `KnotVector` (finite, nondecreasing, exact clamped-end multiplicity
  `degree+1`, no knot run above `degree+1`), Cox–de Boor basis
  (Piegl–Tiller A2.2 shape). Knot storage and degree are sealed behind checked
  construction and read-only accessors. `AdmittedKnotVector` borrows one exact
  immutable snapshot, so its domain/span/basis path reuses structural validation
  and safe mutation while the view is live is unrepresentable. Owning
  convenience calls remain fallible; basis allocation is fallible and
  triangular degree work shares the defensive legacy ceiling.
  `KnotVector::admit_with_cx` returns `KnotAdmissionRun` and publishes the
  lifetime-bound admitted authority only after every validation pass and its
  final checkpoint complete.
  `AdmittedKnotVector::span_with_cx` returns transactional `KnotSpanRun` state:
  the complete scalar span index or cancellation with no published index.
  `AdmittedKnotVector::basis_with_cx` returns transactional `BasisRun` state:
  complete span/values or cancellation with no partial row.
- `NurbsCurve<S, DIM>` — homogeneous de Boor evaluation in dimensions 0–3;
  construction rejects non-finite homogeneous products rather than accepting
  finite source values whose multiplication overflowed. f64 derivatives
  pass both the explicit legacy order ceiling 64 and structure-sensitive
  retained-net work plus 64 MiB payload ceilings, with fallible reservation
  (homogeneous differencing +
  rational Leibniz, including nonzero rational derivatives above polynomial
  degree). At an interior repeated knot, this API returns only ordinary
  derivatives through the knot's continuity order and refuses requests that
  would silently mean a one-sided jet;
  checked Cartesian or homogeneous construction seals the knot/control
  representation. Owning `try_clone_with_cx` returns transactional
  `CurveCloneRun` state and publishes only a complete sealed copy after a
  checked work and 64 MiB retained-output gate. `AdmittedNurbsCurve` binds a
  validated immutable snapshot and evaluation consumes its admitted knot view
  without a second structural scan. Owning `admit_with_cx` returns
  `CurveAdmissionRun` and gates the
  lifetime-bound authority after knot and control validation. Its admitted-only
  `eval_with_cx` returns transactional
  `CurveEvaluationRun` state and never publishes a partial Cartesian point.
  Its f64-only admitted `derivatives_with_cx` returns transactional
  `CurveDerivativesRun` state and publishes only a complete Cartesian jet.
  Its admitted-only `insert_knot_with_cx` returns transactional
  `CurveInsertionRun` state and publishes only a complete validated exact
  refinement generation.
  Its admitted-only `remove_knot_with_cx` returns transactional
  `CurveRemovalRun` state and publishes only a complete validated exact
  coarsening after a cancellable full-representation reinsertion verifier.
  The open-domain check and count-derived aggregate work plus 64 MiB
  simultaneously-live derived-storage envelope precede cancellation; one
  gate then spans occurrence discovery, fallible knot/control allocation,
  exact reconstruction, derived validation, reinsertion, comparison, and final
  publication; restored verifier storage is dropped before the final
  checkpoint.
  Its admitted `span_boxes_with_cx` returns transactional
  `CurveSpanBoxesRun` state and publishes only the complete ordered box table.
  Its admitted-only `to_bezier_form_with_cx` returns transactional
  `CurveBezierRun` state and publishes only a fully validated exact derived
  generation.
  EXACT Boehm `insert_knot` under checked aggregate work and a 64 MiB derived
  output envelope; EXACT `remove_knot` under checked aggregate work and a
  64 MiB simultaneously-live derived-storage envelope (reconstruction checked
  with scalar EQUALITY — in `Rat` a proof); Bézier decomposition; EXACT
  `elevate_degree` (per-segment binomial elevation reassembled on a
  full-multiplicity knot vector — valid and evaluation-identical, including
  legal discontinuous full breaks whose independent endpoints and raised
  multiplicity are preserved; minimal-multiplicity elevation is a follow-up).
- `NurbsSurface<S>` — sealed tensor-product representation with checked
  Cartesian/homogeneous construction. Owning `try_clone_with_cx` returns
  transactional `SurfaceCloneRun` state and publishes only a complete sealed
  copy after checked work and a 64 MiB retained-output gate.
  `AdmittedNurbsSurface` reuses one source validation through tensor evaluation,
  first partials, and per-span control boxes. Owning `admit_with_cx` returns
  `SurfaceAdmissionRun` and gates the
  lifetime-bound structural authority after U-knot, V-knot, and row-major
  control-net validation. Its admitted-only `eval_with_cx` returns transactional
  `SurfaceEvaluationRun` state and never publishes a partial Cartesian point;
  its f64-only admitted `partials_with_cx` returns transactional
  `SurfacePartialsRun` state and publishes value plus both first partials only
  as one complete result;
  its admitted `span_boxes_with_cx` returns transactional
  `SurfaceSpanBoxesRun` state and publishes only the complete U-major, V-minor
  box table;
  its admitted directional insertions reuse the immutable source authority,
  while `insert_knot_u_with_cx` / `insert_knot_v_with_cx` return transactional
  `SurfaceInsertionRun` state and publish only a complete validated exact
  derived surface.
- `TrimLoop`/`TrimmedPatch` — trim curves in EXACT RATIONAL form
  (closure and continuity validated by rational equality, including exact
  left/right-limit agreement at any full knot break). Loop, curve, and
  subdivision storage is sealed; borrowed admitted views bind closure and
  structure validation to the immutable source. Admitted
  `reversed_for_hole_with_cx` returns transactional `TrimLoopReversalRun`
  state and publishes only a complete opposite-orientation loop after checked
  aggregate work/retained-storage admission and full derived validation.
  Aggregate loop-count/structure validation is admitted before the first deep
  scan and spends the same defensive classification budget;
  shell ingestion shares a bounded construction-validation ledger across all
  trims. `classify` is CERTIFIED:
  outside every Bézier span hull ⇒ the exactly-computed control-polygon
  winding equals the curve winding (hull-homotopy argument); ambiguous
  points subdivide EXACTLY (rational midpoints) and become `Boundary`
  after the budget — never a guessed in/out.
- `closest_point_curve`/`closest_point_surface` — best-first
  branch-and-bound over rational-Bézier hull boxes (positive weights
  survive de Casteljau, so every sub-segment hull bounds it), one local
  split per iteration, split junctions as free upper-bound samples,
  Newton polish on curves. Returns
  `DistanceBracketEstimate { lower, upper }` with hulls heuristically
  expanded by one ULP. Cartesian division, evaluation and norm
  arithmetic are ordinary f64, so dense-oracle containment is measured
  evidence rather than a rigorous enclosure.
  An admitted curve or surface also exposes `closest_point_with_cx`, whose
  `ClosestPointRun` publishes either that complete estimate or `Cancelled`;
  cancellation never publishes a partial frontier-derived bracket.
- `boolean(op, policy)` — THE BOOLEAN POSITION: always a structured
  `BooleanRefusal` in v0. Default policy routes through SDF (convert →
  implicit CSG → re-fit); `DirectCertificateGated` refuses pending a
  coverage-complete continuum watertightness certificate. The current
  sampled sheaf interface-agreement evidence is insufficient. An attempt,
  never a promise.

- `sdf` module (plan §7.3 edge 3, bead wqd.11; [F], behind the
  `nurbs-sdf` feature until its Gauntlet tier is green): the current measured
  NURBS → distance-field approximation. `ShellSdf::distance` uses
  convex-hull branch-and-bound plus damped Gauss–Newton, but its bracket is
  not outward-rounded. When the found point is trimmed away or in the boundary
  band, the query carries an infinite upper value and the chart returns
  `NoClaim`; distance-to-kept-region search remains successor work.
  `ShellSdfChart` construction validates the configured aggregate query work
  and support padding and returns a structured error if it could only refuse at
  evaluation time. `ShellSdf::control_aabb` likewise rejects malformed padding
  and any outward expansion that leaves the finite AABB domain through `Result`
  rather than admitting unbounded support or panicking on public configuration.
  It presents `nurbs-sdf/estimated-signed` under declared local
  orientation or `nurbs-sdf/estimated-unsigned` otherwise; it emits
  `NumericalKind::Estimate`, no Lipschitz authority, and no continuity claim:
  finite-budget best-first selection can switch sampled witnesses.
  `generate_tile` is effort-adaptive under defensive static sample/split
  ceilings, not caller-budgeted (P4 remains successor work): refinement fires
  within two cell diagonals of the surface, and achieved pre-storage measured
  widths, outward-expanded f64→f32 quantization errors, and branch-and-bound
  splits are reported per tile. Finite underflow may produce a signed-zero f32
  sample, but its absolute loss remains visible in the storage-error field. A
  trim-downgraded cell emits
  an infinite sentinel/no usable distance rather than the finite distance to a
  point that is not on the kept surface.

- `refit` module (plan §7.3 edge 4, bead wqd.12; [F], behind the
  `nurbs-refit` feature — which enables `nurbs-sdf`, since the CSG
  sampled forward-converter and interface-evidence loop): the SDF → NURBS
  RE-FIT, the edge that makes the
  honest Boolean policy work (§7.2: Booleans route through F-rep, then
  re-fit when a spline chart is required). v1 pipeline for star-shaped
  domains: radial-bisection sampling ON THE FIELD ITSELF, tensor-product
  B-spline least squares with discrete thin-plate (control-lattice
  Laplacian) regularization, exact G⁰ seam closure by control-column
  tying (G¹ measured). ERROR HONESTY: the spline→field direction separately
  reports `max sampled |f(S)|` and the geometric probe-spacing estimate
  `(L_u+L_v)·Δ/2`, where `Δ = 1/probe` is the per-axis probe spacing, with
  each hodograph coefficient evaluated using its actual clamped knot span,
  `L_u = max_i p‖ΔC_i‖/(u_{i+p+1}−u_{i+1})` (and analogously for v).
  The generic closure API supplies no units, admitted Lipschitz witness, or
  metric-error-bound witness, so it does not add these potentially
  dimensionally incompatible quantities. All arithmetic is ordinary f64 rather
  than outward-rounded. Neither quantity is a continuum field enclosure or a
  Hausdorff certificate. The other direction is the worst retained radial
  sign-bracket-target fit residual; a generic closure supplies no continuity
  or root-existence witness, so this is not named source-surface coverage.
  The report records all three estimates. Features encountered by retained
  projection rays can produce STRUCTURED WARNINGS with parameter and world
  locations; a feature missed by every ray is invisible to this API. The
  patch-density knobs are the ErrBudget trade, ledgered.
  The G¹ seam-angle diagnostic evaluates the open `v` interior at the exact
  two `u` seam endpoints. Its report explicitly marks `v=0,1` as excluded,
  because pole-chart tangent directions may be undefined; those endpoints gain
  no implied G¹ authority pending the chart-aware pole audit. After constructing
  the fitted surface, the report binds one admitted immutable snapshot across
  every dense probe, hodograph estimate, and seam-partial evaluation; construction
  validation remains a separate stage, but the report does not repeatedly scan
  the whole control net and knot vectors.

## Invariants

1. **Refinement exactness (the definitive test)**: knot insertion,
   degree elevation, and surface insertion leave evaluation EXACTLY
   equal (rational `==`) at common parameters; insert+remove recovers
   the control net exactly.
2. **G0 laws in exact arithmetic**: partition of unity sums to exactly
   1; clamped endpoints interpolate exactly.
3. **Derivatives match duals**: the generic evaluator instantiated at a
   dual scalar equals the analytic derivative pipeline.
4. **Trim classification never lies**: adversarial fixtures (nested
   square/diamond/island, near-tangent vertices, sliver holes) classify
   correctly; on-boundary points report `Boundary`.
5. **Closest-point estimates cross-check sampled oracles**: dense-sampling
   oracles (1e5 curve samples, 300² surface grid) fell inside `[lower, upper]`
   in the retained fixtures. This is measured evidence, not a universal
   containment theorem.

## Error model

`NurbsError`: `Structure`, `Domain`, `Exactness`. Fallible exact helpers use
`Exactness` for an unrepresentable reduced result. Current `Rat` operator-trait
overflow still panics with a named message because those traits cannot return a
typed refusal (exactness-domain exit — a documented boundary, never wraparound).

## Determinism class

**D0**: exact arithmetic on the `Rat` path; fixed iteration orders
everywhere; best-first ties broken by a monotone logical insertion identity.

## Cancellation behavior

Trim and closest-point subdivision loops carry explicit static iteration
limits. Direct `KnotVector` domain/span/basis admission charges live-knot
validation, worst-case span search, and Cox–de Boor triangular degree work
before allocating or iterating. Borrowed admitted knot, curve, and surface
views make repeated source validation unnecessary for basis, evaluation,
partials, and span boxes. The admitted-only `span_with_cx` path preserves
checked work-refusal then parameter-refusal precedence, polls the directional
linear search after at most 64 span steps, and gates final index publication.
It allocates no payload and does not claim owning admission, wall-time
preemption, caller-budget consumption, or executor drain/finalize authority;
individual generic-scalar comparisons remain non-preemptible. The admitted-only
`basis_with_cx` path preserves
constant-time parameter/work-refusal precedence, polls before allocation and
publication plus every 64 logical span/initialization/triangle/finite-check
operations, and drops all scratch on `BasisRun::Cancelled`. It deliberately
does not claim caller-budget consumption or executor drain/finalize authority.
`AdmittedNurbsCurve::eval_with_cx` carries the same cancellation gate through
basis construction, polls homogeneous accumulation every 64 logical scalar
updates, and gates final Cartesian publication. `CurveEvaluationRun::Cancelled`
contains no partial point; owning curve admission and caller-owned affine
budget/drain/finalize semantics remain outside this primitive claim.
`NurbsCurve::try_clone_with_cx` preserves count-derived copy-work and 64 MiB
retained-output refusal precedence, then carries one fixed-stride gate through
fallible knot allocation/copy, fallible control allocation/copy, and final
publication. `CurveCloneRun::Cancelled` exposes no partial copy and drops all
partial output storage. The borrowed source is excluded from the output
envelope; allocator calls, scalar/array copies, and destructors are
non-preemptible, and the primitive adds no source revalidation, wall-time,
exact caller-budget, drain/finalize, resumability, or geometric-certificate
claim.
`NurbsSurface::try_clone_with_cx` likewise preserves count-derived copy-work
and 64 MiB retained-output refusal precedence before carrying one gate through
the ordered U-knot copy, V-knot copy, outer row table, every fallible inner-row
allocation and row-major control copy, and final publication.
`SurfaceCloneRun::Cancelled` exposes no partial surface and drops all nested
output storage. The retained envelope includes both knot payloads, every inner
row header, and every homogeneous control while excluding the borrowed source.
Allocator calls, scalar/array copies, and nested-vector destruction are
non-preemptible; the primitive performs no source revalidation and claims no
wall-time bound, exact caller-budget consumption, executor drain/finalize,
resumability, surface regularity, topology, or geometric certificate.
`AdmittedNurbsCurve::derivatives_with_cx` preserves request, ordinary-continuity,
checked work, and retained-memory refusal precedence, then polls before every
allocation and through derivative-net copies and differences, reduced knot and
basis scans, homogeneous accumulation, the rational quotient recurrence, and
final jet publication. `CurveDerivativesRun::Cancelled` carries no partial jet
and drops all scratch. The f64 scalar operations and allocator calls themselves
are not preemptible, the conservative legacy envelope is not an exact caller
work-budget ledger, and owning admission, one-sided jets, executor
drain/finalize, and resumability remain outside the claim.
`AdmittedNurbsCurve::insert_knot_with_cx` preserves open-domain, checked
aggregate-work, and 64 MiB derived-output refusal precedence, then carries one
gate through the admitted span search, output allocations and copies,
homogeneous blends, both derived structural-validation passes, and final
generation publication. `CurveInsertionRun::Cancelled` carries no partial
curve and drops all derived storage. The borrowed source is excluded from the
output envelope; allocator calls and individual generic-scalar operations are
not preemptible, and this primitive adds no wall-time, exact caller-budget,
owning-admission, drain/finalize, resumability, or geometric-certificate claim.
`AdmittedNurbsCurve::remove_knot_with_cx` preserves the open-domain check and
count-derived aggregate-work plus 64 MiB simultaneously-live derived-storage
refusal before observing cancellation. Its fixed-stride gate covers the
ordered knot-occurrence scan, fallible knot and reconstruction/control
allocations, the exact forward recurrence and meet check, both candidate
validation passes, restoring insertion, full representation comparison, and
final candidate publication. Restored-generation storage is dropped before the
final publication checkpoint, but the individual destructor is non-preemptible.
`CurveRemovalRun::Cancelled` contains no partial curve. Interior-knot absence
is discovered by the cancellable source scan, so an already-observed
cancellation may win before that data-dependent refusal. The borrowed source
is excluded from the derived-storage envelope; allocator, destructor, and
individual generic-scalar operations are not preemptible, and this primitive
adds no wall-time, exact caller-budget, owning-admission, drain/finalize,
resumability, or geometric-certificate claim.
Admitted curve and surface `span_boxes_with_cx` preserve their existing checked
traversal-work and 64 MiB retained-output envelopes before carrying bounded
polling through output allocation, candidate spans (including skipped
zero-width spans), Cartesian projection/bound updates, and final table
publication. Their Cancelled variants expose no partial box vector. Individual
allocator calls and generic scalar division/comparison remain non-preemptible;
the APIs add no wall-time, Cx-budget, drain/finalize, tight-box, topology, or
certificate claim, and the documented exact-scalar overflow boundary is
unchanged.
`NurbsCurve::admit_with_cx` preserves dimension and static work-refusal
precedence, then carries one cancellation gate through both knot validation and
the homogeneous-control weight, finite, quotient, and inactive-lane scans. A
final checkpoint gates `CurveAdmissionRun::Complete`; construction, `Cx` budget
consumption, and executor drain/finalize remain outside the admission claim.
`AdmittedNurbsCurve::to_bezier_form_with_cx` preserves the existing checked
work and two-generation retained-byte envelope, then polls exact planning,
source copies, repeated target scans, every insertion copy/blend, both derived
knot-validation passes, control validation, and final publication at phase
boundaries and after at most 64 logical operations. Cancellation drops all
partial owned generations and returns `CurveBezierRun::Cancelled`. Individual
allocator calls and scalar operations are not preemptible; the API claims a
logical-operation bound, not a wall-time bound, and does not claim owning curve
admission, `Cx` budget consumption, or executor drain/finalize authority.
`AdmittedNurbsSurface::eval_with_cx` evaluates the U basis and then the V basis
under one cancellation gate, preserving the synchronous U-major/V-minor tensor
arithmetic order. It polls the pair product and homogeneous lane updates after
at most 64 logical operations and gates final Cartesian publication.
`SurfaceEvaluationRun::Cancelled` carries no partial accumulator or point and
drops any allocated basis workspaces. `NurbsSurface::admit_with_cx` preserves
checked static work-refusal precedence, then carries one cancellation gate
through U-knot validation, V-knot validation, grid-shape and row-major control
validation, and final authority publication. It allocates no payload and does
not certify regularity, trimming, topology, closest-point bounds, or rendering.
Construction, aggregate affine-budget consumption, and executor drain/finalize
remain caller responsibilities; allocator calls and individual scalar
operations are not wall-time preemptible. The owning
`KnotVector::admit_with_cx` path applies the
same fixed stride across finite, ordering, multiplicity, and clamping validation
and gates admitted authority at publication; `KnotVector::new` construction
remains outside that cancellation claim. The production `fs-render` ray path
preflights sealed metadata and cancellation, then binds one admitted surface
across domain lookup, seed evaluation, and Newton partials.
Admitted U/V surface insertion preserves open-domain, conservative aggregate
work, and 64 MiB requested derived-payload refusal precedence, then uses one
gate for directional span lookup, direct tensor Boehm row/control assembly,
both knot copies, inserted-knot and complete surface validation, and final
publication. `SurfaceInsertionRun::Cancelled` carries no partial surface and
drops all derived storage. The payload model includes both knot arrays, every
row-vector header, and all homogeneous controls while excluding the borrowed
source; allocator rounding and individual generic-scalar operations are not
preemptible. These methods claim no wall-time, exact caller-budget,
owning-admission, drain/finalize, resumability, or geometric certificate.
`AdmittedNurbsSurface::partials_with_cx` preserves U-then-V parameter,
aggregate-envelope, and ordinary-derivative refusal order before carrying one
gate through both basis rows, the sequential U and V isocurve contractions,
both shared curve-derivative engines, and final tuple publication.
`SurfacePartialsRun::Cancelled` exposes neither the value nor one directional
jet and drops all temporary nets and derivative scratch. It proves only one f64
value/first-partial request on an admitted surface: construction, owning
admission, one-sided jets, regularity, normals/orientation, geometric
certificates, exact caller-budget consumption, wall time, drain/finalize, and
resumability remain outside the claim.
`TrimLoop::admit_with_cx` carries one gate through exact curve/knot admission,
both endpoint evaluations, full-break continuity traversal, and final
lifetime-bound authority publication. `TrimLoopAdmissionRun::Cancelled`
publishes no admitted view. Polling is fixed-stride between logical knot-run
and continuity checks; an individual exact-rational operation is not
preemptible. The primitive does not consume caller budget or own surrounding
request drain/finalize, wall-time, or resumability semantics.
`AdmittedTrimLoop::reversed_for_hole_with_cx` preserves count-derived aggregate
work and simultaneously-live retained-storage refusal precedence, then carries
one `Cx` through fallible reversed-knot allocation, exact same-sign-safe knot
mirroring, fallible reversed-control allocation, ordered fixed-stride copies,
full derived loop admission, and final owned publication.
`TrimLoopReversalRun::Cancelled` exposes no partial loop and drops all derived
storage. The borrowed source is excluded from the retained envelope; allocator
calls, individual exact-rational operations, and destructors are
non-preemptible. The primitive adds no `Cx` budget consumption, wall-time,
drain/finalize, resumability, topology, or geometric-certificate claim.
`TrimmedPatch::admit_with_cx` retains the constant-time minimum loop-count work
refusal ahead of cancellation, then carries the same caller gate through the
exact aggregate validation-work scan, every nested loop/curve admission, and
final patch-authority publication. Plan and loop-table traversal polls at most
every 64 logical entries in addition to the nested curve/loop gates.
`TrimmedPatchAdmissionRun::Cancelled` exposes no partially admitted loop table.
Individual exact-rational operations remain non-preemptible, and the primitive
does not consume the `Cx` budget or own caller drain/finalize, wall-time, or
resumability semantics.
Owning and admitted trim `classify_with_cx`/`classify_box_with_cx` preserve the
componentwise box-order and constant count-only work refusals ahead of
cancellation. One caller gate then spans owning patch admission when needed,
persistent-source accounting, exact U-then-V witness construction, Bezier plan
and conversion scans, span-box assembly, ordered offending-interval scans,
midpoint insertion and reconversion, exact control-polygon winding, cleanup,
and terminal verdict publication. `TrimClassificationRun::Cancelled` carries
neither a partial winding nor `Inside`, `Outside`, or `Boundary`. Linear loop,
box, interval, and polygon traversal polls at most every 64 logical entries in
addition to the nested curve gates. Individual allocations, exact-rational
operations, and derived-curve/vector destructors remain non-preemptible. These
paths retain the static defensive work ceiling and do not consume an affine
`Cx` budget, own request drain/finalize, promise wall-time preemption, or expose
resumable subdivision state.
Owning trim classification now binds one admitted patch/loop/curve generation
through exact Bezier conversion, span boxes, and winding. Its checked conversion
plan charges scan/insertion work and old-plus-new curve storage before the first
clone; all patch-owned source curves, each subdivision insertion, and the
simultaneously retained converted curve, span boxes, and offending intervals are
admitted under a 64 MiB defensive peak before allocation. Offending intervals
use fallible reservation and winding projects controls on demand without
retaining a second polygon. An admitted patch can classify repeatedly without
revalidating its immutable source.
Curve closest-point execution now validates its immutable source once, consumes
the exact Bezier conversion plan, pre-reserves the admitted heap extent, and
reuses the admitted source for optional derivative polish. Its 256 MiB payload
ceiling composes the borrowed source with the maximum of conversion,
converted-curve plus queue/frontier/scratch, and post-release derivative-polish
phases; the converted curve and queue are dropped before polish. The owning
wrapper rejects malformed requests before scanning the source, while an
`AdmittedNurbsCurve` can repeat closest-point calls without rescanning source
structure. Its admitted-only `closest_point_with_cx` path keeps
malformed-request and count-only pre-scan work refusals ahead of cancellation,
then uses the same `Cx` through Bezier planning/conversion, seed copies and
hulls, heap traversal, de Casteljau splits, Newton polish, final evaluation,
and publication. Linear and triangular logical work polls at most every 64
operations, with explicit gates around owned allocations and bounded heap
operations. `Cancelled` publishes no partial estimate; cancellation inside
optional derivative/evaluation kernels is not downgraded to a polish miss.
Individual allocator, heap, scalar, and destructor calls remain
non-preemptible. This measured primitive does not consume the caller budget,
own request-drain-finalize, promise wall-time preemption, or provide resumable
frontier state. Surface closest-point execution follows the same
malformed-request before source-scan order and exposes an admitted repeat path.
It reuses that source admission through conversion planning and final
evaluation, pre-reserves the exact `seed_leaves + max_splits` heap extent,
bounds every push, and prices
log-height heap operations. Its stage-faithful knot-insertion, expanded-grid,
run-scan, and queue-seeding estimate is charged before conversion, including
when the requested split budget is zero; degree-scaled de Casteljau split work
and the worst retained queue/scratch frontier are admitted before allocation.
The conversion preflight reuses the directional insertion engine's own
per-generation work and 64 MiB derived-output envelope at the largest planned
generation, so an admitted closest-point request cannot first discover that
nested refusal after conversion begins.
Its 256 MiB requested-payload ceiling composes the borrowed source with the
maximum of exact-conversion allocations, the converted surface plus traversal
frontier, and that retained search state plus final basis-evaluation workspace.
The conversion bound covers surface row tables and the exact overlap of the
largest old/new direct tensor insertion generations; it does not claim
allocator metadata, rounding, or pre-existing spare source capacity.
The admitted surface `closest_point_with_cx` path keeps request and count-only
pre-scan refusals ahead of cancellation, then polls the ordered U/V planning
scans, fallible source clone, alternating exact U-then-V conversion, U-major
and V-minor seed copies/hulls, both patch-split axes, heap traversal, optional
center evaluation, cleanup, and final publication. Fixed-stride loops observe
the gate at most every 64 logical operations; `Cancelled` never carries a
partial frontier or bracket, and cancellation from optional evaluation is not
treated as an ordinary evaluation miss. Converted storage and the frontier
remain live through final evaluation as admitted by the aggregate envelope.
Individual allocations, heap calls, scalar operations, and nested-`Vec`
destructors remain non-preemptible; cleanup latency, caller-budget consumption,
wall-time preemption, resumability, and request-drain-finalize ownership remain
explicit no-claims.
Owning derivative and refit construction paths are not all
migrated yet; they make no claim of caller-budgeted preflight or end-to-end
validate-once execution. The SDF shell rejects malformed point/tolerance input
before surface planning, admits each immutable surface once per distance query,
reuses that admission through closest and Gauss-Newton polish, and carries the
winning admission into gradient, orientation, and regular-witness sign work.
It reuses the same split/frontier model, exports worst-case seed and split heap
coefficients into its aggregate work gate, and additionally charges
structure-sensitive polishing, sign-repair, and trim coefficients; a shell for
which even a zero-split query exceeds the ceiling is not constructible. The SDF
and legacy refit paths retain static process ceilings with no effective `Cx`.
Trim classification now has effective primitive checkpoints but still uses its
static ceiling rather than a caller-owned affine budget. None of these APIs yet
owns end-to-end request→drain→finalize semantics, so they are not P7
cancellation-correct. Successor budgeted interfaces must close those remaining
orchestration gaps before promotion.

Structural admission proves representation well-formedness only. It does not
grant enclosure, geometric certification, cancellation, caller-budget,
allocation-completeness, replay, or downstream publication authority.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

All OFF by default per the Ambition-Tag rule:
- `nurbs-sdf` [F] — the NURBS→SDF converter; disabled until its
  Gauntlet tier is green. Gates the `sdf` integration target.
- `nurbs-refit` [F] — the SDF→NURBS refit (implies `nurbs-sdf`; its current
  acceptance evidence is sampled and does not close a continuum certificate).
  Gates the `refit` integration target.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-rep-nurbs/conformance`):
nb-001 exact G0 laws + derivative-vs-dual through the generic evaluator;
nb-002 six random rational curves: insertion/elevation evaluation-exact +
insert/remove lossless round trip; nb-003 adversarial trim battery
(nesting, tangency, slivers, boundary honesty); nb-004 measured curve +
surface brackets vs dense oracles plus non-finite, scaled-norm,
adjacent-float termination, multiplicity, and rational-high-derivative
regressions (iteration stats logged, no enclosure claim); nb-005 the
Boolean policy refusals with teaching routes; nb-006 exact surface
refinement + partials vs central differences.

## No-claim boundaries

- **No direct B-rep Booleans** — by design, not omission. The SDF route
  is the product; the direct path requires the successor coverage-complete
  continuum watertightness certificate. The sampled interface evidence from
  wqd.13 is intentionally insufficient for that authority.
- **Degree elevation emits full-multiplicity knot vectors** (valid,
  evaluation-identical, with discontinuous full breaks preserved rather than
  silently healed); minimal-multiplicity reassembly is follow-up.
- **Closest-point and NURBS→distance brackets are measured estimates.** The
  one-ULP hull expansion is heuristic and cannot authorize `Enclosure`, exact sign,
  a 1-Lipschitz field, or no-tunneling. The fs-ivl/Taylor path with outward
  rational projection, norm bounds and interval Newton is the certified upgrade.
- **Greville/Gauss quadrature tables for IGA** land with the fs-iga
  consumer (tfz.9), which owns the quadrature accuracy claims.
- **`Rat` is i128-bounded**: deep repeated refinement can exceed the
  exactness domain; the current failure is loud and named. Converting that
  boundary to a fully fallible exact-arithmetic API remains required before
  hostile inputs can be treated as an ordinary data path.
- Ray intersection shares this machinery but ships with the LUMEN
  chart-backend bead (qfx.2).

## No-claim boundaries (sdf converter)

- Lower estimates come from f64 control-hull boxes; near MEDIAL AXES (many
  equidistant patches) and pole-degenerate parameterizations, bracket
  widths converge slowly with splits — the documented budget trade
  (~1e-3 at 2000 splits/cell, 2.6e-4 at 8000 on the unit-sphere tile).
  Interval-Newton (Krawczyk) contraction of the projection equations is
  the upgrade path when fs-ivl grows 2-D machinery.
- Sign is estimated from the declared local B-rep orientation; the winding-style
  fallback for imperfect shells lands with the quarantine/census beads
  (fs-io owns mesh-side honesty).
- Trim downgrades widen the estimate; distance-to-kept-region
  (excluding trimmed area from the B&B itself) is future work — the
  current f64 lower value has no enclosure authority even for the untrimmed
  surface.

## No-claim boundaries (refit)

- v1 parameterization is RADIAL (star-shaped domains around the given
  center; the bracket failure is a structured teaching error).
  General-topology segmentation over the dual-contoured mesh (wqd.10)
  is the upgrade path.
- The generic refit accepts an arbitrary scalar closure and therefore grants
  no field, zero-set, or distance authority. Even a proven 1-Lipschitz field
  satisfies `|f(x)| ≤ dist(x,Z)` (the opposite direction from the upper
  geometric bound a Hausdorff claim needs); a geometric promotion requires a
  separate metric-error-bound witness such as `dist(x,Z) ≤ κ|f(x)|`, plus
  outward-rounded evaluation and coverage evidence. Min/max CSG fields do not
  acquire that authority automatically.
- Retained radial sign-bracket-target → spline evidence is sampled at the
  projection grid: a feature no ray hits is invisible, and without admitted
  continuity/root-existence evidence the target is not authorized as a point
  on a source surface. Density is the caller's knob. The warning channel
  reports only a retained paired-parameter residual above threshold; it does
  not distinguish thin geometry from smoothing, inadequate global density,
  conditioning, noise, or stateful caller-field behavior.
- The legacy closure API uses validated dimensions plus a conservative live-
  payload allocation estimate, fallible checked reservations before initializing
  every stage-owned numerical `Vec`/matrix buffer, and fixed
  probe/algorithmic-work ceilings so malformed configurations refuse before
  evaluating the field. Allocator metadata, diagnostic strings, arbitrary
  closure cost, and every small transient are not a complete memory or time
  budget. These caps are process constants, not caller budgets;
  the dense fit has no `Cx` cancellation points. Typed admitted-field authority,
  ledgered caller budgets, bounded cancellation latency, and outward-rounded
  geometric certification remain explicit successor work.
- G¹ cross-seam continuity is measured, not enforced; pole rows are
  least-squares-collapsed, not pinned.
