# CONTRACT: fs-geom

## Purpose and layer
The Region/Chart abstraction (plan §7.1): abstract regions presented
through charts; agreement between presentations as a checkable, localized
proposition; certified conversion receipts (fs-evidence) as the Error
Ledger's geometry feed; no privileged representation, ever. Layer: L2.
Depends on fs-exec (Cx), fs-blake3 (typed canonical identities), fs-evidence,
fs-ivl, fs-alloc, fs-obs, fs-sparse.

## Public types and semantics
- `Point3`/`Vec3`/`Aabb` — minimal geometry-local types (fs-la owns real
  linear algebra); `Aabb` normalizes ordered numeric corners without
  laundering NaNs, offers containment/union/inflation/intersection, and uses
  `WHOLE_SPACE` for honest unbounded support. Set union preserves a malformed
  public operand for structured admission rather than laundering it into a
  plausible finite support.
- `SamplingDomain` is the mandatory finite-domain admission boundary before
  midpoint, span, diagonal, count, allocation, or sampling arithmetic. It
  validates raw extended supports before set operations, resolves unbounded
  supports only through an explicit finite positive-volume clip, and returns
  axis-attributed `SamplingDomainError` refusals for malformed, unresolved,
  degenerate, or non-representable domains.
- `ClippedChart` represents the geometric intersection of a source chart and
  a finite clip with `max(source_field, exact_box_sdf)`. Its support and sign
  are honest; its composite magnitude, gradient ties, abstract-distance error,
  and ray-step theorem retain conservative C0/`NoClaim` semantics.
- `BettiBounds` — per-dimension (lower, upper) topology hints;
  `unknown()` is the honest default (certificates are the wqd.7 bead's).
- `ChartSample { signed_distance, gradient: Option, lipschitz: Option,
  error: NumericalCertificate }` — plan Appendix B's value + gradient +
  certified Lipschitz + DECLARED error model relative to the abstract
  region. SDF sign convention: negative inside.
- `TraceStepClaim::{NoClaim, ExactDistance, LipschitzImplicit}` — the typed
  theorem available to a ray stepper. The default is `NoClaim`: a sample-level
  `Some(lipschitz)` alone cannot mint a no-tunneling certificate.
  `ExactDistance` states that the represented real field is the exact signed
  distance and requires either a genuinely exact singleton or a rigorous
  outward enclosure of each rounded evaluation. A stepper uses the enclosure
  endpoint closest to zero for its no-tunneling radius and the farthest endpoint
  for its hit residual. `LipschitzImplicit` states that the field has the
  represented region's exact sign and zero set, is continuous on every finite
  line segment, and that each sample's bound is valid over the entire closed
  `|f|/L` step ball. Its separate
  `trace_value_enclosure` rigorously encloses `f(p)`; `ChartSample.error` remains
  relative to abstract Euclidean signed distance and may honestly be only an
  `Estimate`. The resulting radius is safe but not a geometric-distance upper
  bound. Finite-segment continuity does let consumers use rigorously opposite
  endpoint signs as existence evidence for a zero inside a short segment.
- `Chart` (object-safe): `eval(x, &Cx)`, `support()`, `trace_step_claim()`,
  `trace_value_enclosure(x, sample, &Cx)`, `topology_hint()`, `name()`,
  `differentiability()`, provided `inside()`.
  Implementations
  poll `cx.checkpoint()` at bounded strides. The plan's `type Param`
  lives on the `DesignChart: Chart` subtrait so `Region` can hold
  heterogeneous `Arc<dyn Chart>` (same contract, object-safe core;
  fs-xform builds on `DesignChart`).
- `Region` — charts + per-chart `ProvenanceHash`; deterministic
  `primary()` (first); `check_agreement(&AgreementConfig, &Cx) ->
  Result<AgreementReport, Cancelled>` — seeded sampling over the inflated
  union support, or over `AgreementConfig::sampling_clip`, with an explicit
  `Agreed | Disagreed | Unknown` verdict.
  Two valid chart samples agree at x iff their declared signed-distance
  intervals overlap after configured slack. Zero samples, fewer than two
  charts, invalid configuration/support, non-finite outputs, malformed
  certificates, and `NoClaim` all produce `Unknown`, never vacuous agreement.
  Reports retain the weakest certificate class used, the strongest class
  supporting a counterexample, exact total counts, signed worst excess (so a
  negative agreement margin is not rounded away), and first-K localized
  disagreement/unknown diagnostics in canonical JSON. `AgreementReport::scope`
  distinguishes global-support evidence from explicitly clipped local evidence,
  and `sampling_domain` records the exact admitted finite box.
- `Convert<Dst>: Chart` — `convert(ErrBudget, &Cx) ->
  Result<Certified<Dst>, ConvertDiag>` promotes only rigorous global
  abstract-distance evidence, which requires the source's global
  `TraceStepClaim::ExactDistance` theorem in addition to rigorous sample
  certificates. `convert_with_domain` returns plain `Evidence` so weaker
  source fields remain usable as estimates without laundering sampled LOCAL
  Lipschitz values into global authority;
  `convert_clipped` converts the actual `ClippedChart` composite and returns
  plain `Evidence`, because generic `max(source_field, exact_box_sdf)` has no
  abstract signed-distance theorem. The receipt QoI is the total
  reconstruction-plus-source bound when available, and its numerical kind is
  the weakest sampled source authority after interpolation demotes `Exact` to
  `Enclosure`. Exact clip endpoint bits participate in provenance.
  `ConvertDiag` refuses EARLY with ranked fixes (`BudgetInfeasible`
  computes the needed resolution vs the cap; `NoLipschitzBound` when the
  source certifies none; `UnrepresentableGrid` when a translated/tiny interval
  cannot hold the required number of distinct f64 nodes), rejects every
  non-finite source signed-distance sample before publishing a chart or
  `Certified` receipt, and returns a stage/count-attributed `Cancelled` refusal
  when its own `Cx` checkpoints observe cancellation.
- `SampledSdf` — dense trilinear sampled-field target with a finite nominal
  reconstruction bound kept distinct from authority relative to abstract
  region signed distance. It stores strictly increasing representable nodes on
  each axis and locates/interpolates against those actual cells. Its rigorous
  reconstruction radius is the outward product of `L` and the largest actual
  full-cell diagonal (valid for arbitrary convex trilinear weights), plus a
  finite outward interpolation-roundoff allowance. It validates and outwardly
  composes every sampled `ChartSample::error`: malformed claims and `NoClaim`
  fail closed, estimates stay estimates, rigorous source radii are added to
  reconstruction error, and interpolation never claims `Exact`. Even rigorous
  point samples are demoted to `Estimate` when the chart lacks `ExactDistance`,
  because a finite grid cannot establish a global interpolation theorem.
  In-box certificate endpoints and the Lipschitz outside-box extension are
  outward-rounded and contain the published nominal. Outside-box evaluations
  retain the same authority grade; blanket
  `impl<C: Chart> Convert<SampledSdf> for C` (specialized edges arrive
  with rep-* beads). Resolution cap `SAMPLED_SDF_MAX_RESOLUTION = 96`.
- `fixtures` (PUBLIC on purpose — the shared MORPH test vocabulary):
  valid positive-radius `SphereChart`, finite strictly three-dimensional
  `BoxChart`, and valid ring `TorusChart` instances (unit
  Lipschitz, outward-rounded evaluation enclosures, `ExactDistance` trace claim,
  known Betti numbers); degenerate/invalid boxes downgrade to `NoClaim` and unknown
  topology, while horn/spindle torus parameters downgrade to
  `LipschitzImplicit`, retain an abstract-distance `Estimate`, publish a
  separate outward trace-value enclosure, and report unknown topology.
  `LyingSphereChart` is
  deliberately biased with a lying error model and the default `NoClaim` for
  detection tests.

- `derived` module (RD.1a, `[M]`, behind `derived-geometry`) defines the finite
  admitted *object language* for stratified derived machine geometry. A
  `DerivedGeometryIrV1` binds an immutable subject/model version to an explicit
  algebraic or restricted-analytic category, coefficient field/ring and
  real-versus-complex semantics, finite configuration charts, frame and unit
  conventions, locality and compactness, equality germs, ordered inequality
  germs and active sets, relative boundaries, unilateral contacts, constitutive
  metadata, tangent/cotangent/deformation-obstruction complexes, finite
  resolutions, local singular-model classes, a finite stratum poset, compact
  local links, and proof-state metadata. Equality, inequality, contact, and
  constitutive IDs are nominally different Rust types and live in different
  collections; matching digest bytes cannot interchange their roles.
  `admit_derived_geometry_v1(ir, budget, &Cx)` checks collection and aggregate
  rank ceilings before canonical sorting, rejects declared or observed basis
  dimensions above the absolute finite-dimension ceiling, and rejects unbounded
  or infinite-dimensional locality, opaque/infinite computation, unsupported
  analytic callbacks,
  unordered complex inequality/contact semantics, mixed frames/units without
  an explicit successor morphism, malformed complexes/references/incidences,
  and invalid proof scope. It polls cancellation at a fixed nested-rank stride,
  between every canonical sort, per finite object, and while streaming identity
  bytes. Success returns an opaque
  `AdmittedDerivedGeometryV1` with canonically ordered IR and an fs-blake3
  `IdentityReceipt<DerivedGeometryIdV1>`; the receipt is content-addressed
  structural admission, not a theorem or physical-validity certificate.
  Literal/fixed-resolution/external presentation scopes make the RD.1b
  equivalence boundary explicit rather than assuming generator or coordinate
  invariance.

- `derived_morphism` module (RD.1b foundation, `[M]`, behind
  `derived-geometry`) admits one rank-neutral identity arrow per exact geometry
  and structurally strict maps between exact `AdmittedDerivedGeometryV1`
  endpoints. Strict maps require
  equal physical subject, immutable model version, mathematical category,
  coefficient semantics, frame, and unit system; version/coordinate/unit
  transitions are not silently inferred.
  `DeclaredChartMap` is a distinct primitive family: it binds source/target
  chart IDs, one nominal overlap-scope ID, and one nominal forward-map ID.
  Admission resolves both charts in their exact endpoint objects and requires
  equal coordinate/ambient dimensions, endpoint frames, coordinate unit
  systems, quantities, and exact scale bits. These checks establish a typed
  structural map declaration only; they do not prove that the overlap exists,
  that the map implements the declaration, or that it is invertible. Homogeneous
  chart-map paths compose only across an exact middle chart and retain their
  first/last chart IDs. The sealed token exposes the ordered primitive geometry,
  chart, overlap, and map IDs so downstream artifact resolution does not depend
  on the caller's raw request sidecar. It also exposes one ordered typed
  primitive path. Mixed generic-strict/chart-map composition produces a
  `HeterogeneousPath` without erasing primitive family or endpoints; any two
  adjacent chart-map primitives still require an exact chart seam.
  `DerivedChartTransitionInverseLawCandidateIrV1` separately binds two direct
  declared chart-map children with exactly reversed geometry and chart
  endpoints, one exact common nominal overlap, two nominal round-trip
  declarations, and a no-authority artifact. Admission additionally requires
  the declared evidence transports to have one variance and to compose across
  both ordered artifact/rank seams. The packet does not execute either map or
  turn either round-trip declaration into equality, inverse, or equivalence
  authority.
  `DeclaredComplexRefinement` is a coarse-source to refined-target finite-
  complex rank-envelope declaration. It binds exact source/target complex and
  resolution IDs plus nominal aggregate prolongation and differential-
  commutation artifacts. Admission resolves both complexes in their exact
  endpoint objects, requires equal roles and an exactly unchanged chart
  presentation, covers every source degree with the same quantity and no rank
  decrease, and requires strict progress through positive-rank degree growth,
  degreewise rank growth, higher positive truncation order, or removal of
  truncation. An untruncated source cannot regress to a truncated target.
  Homogeneous refinement paths compose only across exact middle complex and
  resolution IDs; the same seam remains enforced when a refinement is adjacent
  to another refinement inside a heterogeneous path.
  `DeclaredInclusion` is a distinct whole-object directed-map family. It binds
  one nominal declared inclusion-map ID and one nominal claimed-containment
  artifact while reusing the exact subject, immutable model version, category,
  coefficient, frame, unit, evidence-orientation/rank, and no-equivalence
  gates. Admission deliberately does not infer containment from chart,
  constraint, complex, local-model, or stratum collection syntax. Homogeneous
  inclusion paths compose through the exact middle geometry; mixed paths retain
  each inclusion primitive and its factor-local artifacts.
  `DerivedStratumMorphismIrV1` belongs to a separate stratum-scoped category,
  not the whole-geometry directed-map algebra. Each object is the exact triple
  `(geometry, stratification, stratum)`. A declared component binds exact source
  and target triples plus nominal component-map and constructibility artifacts;
  admission resolves both selectors in their exact admitted geometries but does
  not require equal charts or dimensions. Component paths compose only across
  an exactly equal middle triple and flatten ordered component/no-authority
  lineage. There is no conversion to `AdmittedDerivedMorphismV1`, heterogeneous
  whole-object composition, or geometry-wide evidence transport. One component
  deliberately says nothing about unlisted strata or a whole stratified map.
  `DerivedExhaustiveStratifiedMapCandidateIrV1` is a standalone finite assembly
  candidate over those sealed components. In canonical source-stratum order it
  requires exactly one direct component binding for every source stratum, with
  exact source and target geometry/stratification selectors and explicit raw-to-
  sealed receipt identity. Target strata may repeat and need not be exhaustive.
  The token retains nominal whole-family and global-constructibility declaration
  IDs plus a no-authority artifact. It proves finite source-selector coverage,
  not that the components execute or glue into a continuous, constructible, or
  incidence-preserving whole map. The base token has no conversion, execution,
  evidence-transport, inverse, or equivalence API.
  `DerivedExhaustiveStratifiedMapCompositionCandidateIrV1` separately retains
  two exact exhaustive-map candidates in caller-significant order. Source,
  middle, and target geometry/stratification selectors are derived from the
  sealed children, and admission requires exact equality of both middle
  selectors. This packet does not compose component paths, synthesize a direct
  exhaustive or whole-geometry map, authenticate assembly/constructibility, or
  grant gluing, continuity, preservation, categorical-law, theorem, inverse,
  equivalence, evidence, or physical authority.
  `DerivedStratificationRefinementCandidateIrV1` binds one exact exhaustive
  fine-to-coarse assembly to exact refined/coarse geometry and stratification
  selectors. Admission requires canonical coverage of every refined stratum,
  at least one selected preimage for every coarse stratum, and no increase from
  a refined stratum's dimension to its selected coarse stratum. A separate
  `DerivedStratificationRefinementCompositionCandidateIrV1` retains two exact
  candidates only when both middle geometry and stratification selectors match.
  Neither token mints a direct refinement map, containment/incidence/link
  theorem, transitivity, evidence transport, or equivalence authority.
  `DerivedParallelMorphismComparisonCandidateIrV1` binds two exact ordered
  structural paths with one exact common source and target plus nominal scope,
  relation, and no-authority artifacts. It is a comparison input only: direct
  and composite paths, identities and cycles, may be retained without claiming
  path equality, commutativity, homotopy, coherence, execution, or equivalence.
  `DerivedSpanCorrespondenceIrV1` is deliberately separate from the directed
  morphism-kind algebra. It declares `source <- apex -> target` by binding exact
  source/apex/target geometry IDs and two already-admitted morphism receipts
  oriented `apex -> source` and `apex -> target`, plus an explicit no-authority
  artifact. The sealed span exposes no direct source-to-target evidence
  transport and has no composition method.
  `DerivedSpanEvidenceTransportCandidateIrV1` separately binds one exact sealed
  span to one exact left-leg restriction-polarity child and one exact right-leg
  balance-corestriction-polarity child. The restriction child must wrap the
  span's `apex -> source` leg, the corestriction child must wrap its
  `apex -> target` leg, and admission requires exact evidence-artifact, rank,
  and polarity equality where the two declared lanes meet at the apex. The
  token derives its three geometry selectors and all outer/apex evidence fields
  from the sealed children. It retains one nominal aggregate declaration and a
  mandatory no-authority artifact but exposes no executable transport, direct
  source-to-target map, functionality, pull-push, base-change, preservation,
  theorem, inverse, equivalence, or physical authority.
  `DerivedSpanMorphismCandidateIrV1` binds two fixed-foot spans, one exact apex
  morphism, and exact ordered comparison children between each source leg and
  the recomposed apex-map-then-target-leg route. It binds the exact structural
  wiring, not either commuting triangle or a span 2-cell.
  `DerivedSpanPullbackSquareCandidateIrV1` binds two structurally composable
  spans, exact projections from one proposed apex, and one exact ordered
  comparison of the recomposed routes to the common middle geometry.
  `DerivedSpanCompositionCandidateIrV1` then binds that exact square, one exact
  proposed outer span, and the exact recomposed routes to both outer feet.
  These packets retain a proposed pullback and outer composite shape only; they
  grant no square commutativity, pullback existence/universality, categorical
  or bicategorical composition, associativity, base-change, pull-push, or
  evidence authority.
  `DerivedMorphismInverseLawCandidateIrV1` binds exact oppositely oriented
  morphisms, derives both endpoint identities and both round-trip composites,
  and requires exact ordered identity-versus-round-trip comparison children.
  Caller-authored identity slots are intentionally absent. The packet grants no
  path equality, inverse, isomorphism, homotopy, simplification, or equivalence.
  `DerivedSpanEquivalenceCandidateIrV1` assembles exact source/target spans,
  ordered forward/reverse fixed-foot span-morphism candidates, and the exact
  inverse-law candidate for their ordered apex morphisms. Admission rechecks
  every outer-foot, span, apex, direction, and apex-map seam transitively. The
  resulting structural packet deliberately cannot convert to
  `DerivedEquivalenceBoundaryV1` and grants no span equivalence, invertible
  2-cell, unit/counit, triangle, coherence, interchange, semantic, or physical
  authority.
  `DerivedFixedResolutionQuasiIsomorphismCandidateIrV1` is another standalone
  declaration outside the morphism-kind and equivalence algebras. It binds one
  exact sealed homogeneous refinement path to exact source/target local-model,
  complex, resolution, and complex-role selectors. Both local models must own
  the selected role complex, share one exact locality, and explicitly carry a
  `FixedResolution` presentation scope matching their own selected resolution.
  Literal and merely `ExternallyChecked` presentation scopes fail closed. The
  sealed candidate retains both endpoint scope-witness IDs plus nominal theorem,
  checker, and external check-receipt IDs and an explicit no-authority artifact.
  These are structured inputs for RD.1c, not a quasi-isomorphism or equivalence
  capability; candidate composition and inversion are intentionally absent.
  `DerivedLocalPresentationCorrespondenceCandidateIrV1` is a separate
  exhaustive finite relation over exact source/target local models. Equality,
  explicitly active inequality, explicitly active contact, and constitutive
  families use distinct typed edge collections; constitutive metadata is never
  reclassified as a geometric constraint. Each edge binds exact source/target
  members plus one nominal relation artifact. Admission requires exact endpoint
  conventions, chart, locality, and two-sided coverage of every member in every
  selected family. Edges are canonicalized by `(source, target, relation)`;
  duplicate pairs refuse even when their nominal relation artifacts differ.
  Repeated sources and targets otherwise remain valid, so the token is a finite
  relation rather than a hidden function or bijection. One nominal aggregate
  declaration and a no-authority artifact are retained. Payload semantics,
  functionality, composition, evidence transport, inverse, and equivalence are
  intentionally absent and require independent RD.1c promotion.
  Across every RD.1b candidate, nominal relation/theorem/comparison/round-trip/
  composition/equivalence declarations and no-authority artifacts are opaque
  nonzero identities. Structural admission retains and domain-separates those
  bytes; it never dereferences, executes, authenticates, or upgrades them.
  Evidence transport is explicitly contravariant restriction or covariant
  balance corestriction and binds exact input/output geometry identities plus
  nominal caller-declared input/output evidence-artifact identities and ranks.
  Admission rejects a declared output rank above the declared input rank.
  `DerivedEvidencePolarityTransportCandidateIrV1` separately binds one exact
  sealed nonidentity morphism child and derives its map endpoints, evidence
  variance, artifact IDs, and ranks from that child. It retains explicit
  `Affirmative`, `Refutation`, or `NoClaim` input/output polarity tags and
  admits only an exactly equal pair. Identity children refuse because they
  carry no nominal evidence artifacts or ranks. Its 12-field parent receipt
  recursively validates the six-field morphism child for an exact 18-field
  limit. The packet prevents unequal polarity tags within one candidate; it
  does not authenticate the chosen tag or stop a caller from proposing a
  distinct nominal packet for the same child.
  Composition requires an exact middle endpoint and matching nominal evidence
  identity/rank seam, flattens primitive factors and no-equivalence artifacts in
  semantic order, and produces an associative content-addressed receipt with
  exact endpoint identities and ordered typed primitive lineage. Homogeneous
  strict, chart, complex-refinement, and inclusion composites retain family-
  specific class encoding; a distinct heterogeneous composite tag plus the
  flattened primitive-factor identities binds mixed-family order. Immediate composite
  operand IDs are intentionally
  not encoded because that would make the receipt depend on parenthesization.
  Only an identity arrow carries identity equivalence; every nonidentity
  primitive/composite retains explicit no-equivalence artifacts. The evidence
  artifact IDs and `ColorRank` values are structural
  declarations with zero payload authority: they do not authenticate evidence,
  establish payload preservation or validity-domain inclusion, prove theorem
  truth, or establish physical equivalence. RD.1b now provides separate
  structural no-authority candidates for direct chart-transition inverse laws,
  ordered exhaustive-map composition, finite stratification refinement and its
  ordered composition, evidence polarity, span restriction/corestriction
  routes, parallel paths, fixed-foot span morphisms,
  pullback/composition-shaped spans, generic inverse laws, and span equivalence;
  none is represented as a strict map or authority-bearing equivalence. The L6
  Machine-IR crosswalk likewise lives at the `fs-ir` feature seam rather than in
  this L2 morphism algebra. Authenticated subobject inclusion, geometric and
  stratified refinement theorems, authoritative constructible/global maps,
  actual span pullbacks/composites/equivalences, quasi-isomorphisms, and physical
  crosswalks require independent RD.1c or higher-layer checking of the retained
  structural candidates and artifacts.

- `exit_path` module (RD.X1, `[M]`, behind `derived-geometry`) admits a bounded,
  executable statement language for maximal exit/entrance-path approximation
  and constructibility theorem families over an admitted RD.1a geometry. Each
  `ExitPathFamilyIrV1` fixes direction, sheaf/cosheaf variance, constructible
  coefficient system, stratified-path equivalence, link fidelity, monodromy/
  local-system fidelity, constructibility, properness, refinement, homotopy
  truncation, theorem lifecycle, TCB, budgets, and a canonical set containing
  all four required adversarial classes. Admission rejects zero identity sentinels,
  any geometry/model/stratification/frame/unit tuple that does not match the
  supplied `AdmittedDerivedGeometryV1`, unsupported/infinite presentations,
  invalid truncations, insufficient referenced-artifact or implication
  budgets, malformed falsifiers, mismatched theorem cards, and dangling
  refutations. It derives every requested node from incidence poset through a
  directed exit/entrance one-category enriched by within-stratum groupoids,
  simplicial category, finite higher truncations, and the full higher category.
  A node is `SufficientStatement` only when every fidelity axis reaches that node's
  degree; otherwise it retains a precise `Unknown` reason. Trivial fundamental
  group/monodromy data reaches only degree one. An unauthenticated refutation
  record names one exact lattice node and does not erase independently eligible
  richer fallbacks or mint a scientific verdict.
  `SufficientStatement` means that the implication is well formed under the
  declared hypotheses, never that the implication has been proved. The
  `ExitPathFamilySnapshotIdV1` intentionally binds lifecycle, TCB, falsifier,
  and budget metadata; it identifies a complete statement/evidence/operation
  snapshot, not a stable theorem statement across checker runs.
  `constructible_coefficients` names the target category for constructible
  (co)sheaves; it is intentionally independent of the admitted geometry's
  coordinate algebra and is therefore identity-bound rather than equality-
  checked against `DerivedGeometryIrV1::coefficients`.

- `router` (the Rep Router, Bet 1): converter-edge registry
  (`ConverterSpec`: cost model, error model with declared composition rule,
  declared certificate availability), bounded Pareto label-correcting planner over
  (cost, composed absolute error, uncertified-edge count) with the
  deterministic winner rule certified-preferred → min cost → min error →
  lexicographic path; `explain()` returns every Pareto candidate and why
  the winner won; refusals name the binding constraint (error/cost/
  no-path) with ranked relaxations; `execute()` runs chains through
  `EdgeRunner`s, composes local-error Evidence receipts with one shared
  directed-rounding algebra (additive sum; relative upstream amplification plus
  local receipt; exact requires zero), and records actuals through `CostOracle`
  (an L2-clean abstraction — HELM
  wires the ledger tune table behind it; `MemoryCostOracle` in-process).
  `ChainOutcome` is opaque outside `fs-geom`; read-only receipt and measured-cost
  accessors prevent callers from fabricating route authority.
  Oracle reads are fallible and scoped to the exact `ConverterSpec`; one-pass
  read snapshots are identity-bound into opaque `RoutePlan`s and rechecked
  before and after execution. `CostOracle::record_batch` is fallible:
  invalid/nonfinite/overflowing/capacity-
  exceeding evidence returns `CostOracleError`, and `execute()` propagates it
  as edge-attributed `ExecuteErrorKind::OracleRecord` instead of reporting a
  successful chain whose actuals were silently dropped. `MemoryCostOracle`
  updates cost sums/counts and observed error maxima atomically and bounds
  distinct specs. Learned observations can only increase an uncertified
  additive declaration; retrospective means/quantiles never tighten hard error
  authority. Router edges/nodes, path length, total/per-node labels, and
  candidate expansions have deterministic caps with typed refusals.
  Identity routes skip empty oracle writes; execution polls cancellation before
  each edge and before evidence persistence. Optional sheaf rerouting retains a
  structured `RoutePlanError` instead of silently dropping malformed authority.

- `sheaf` module (bead wqd.13, Bet 11 [F/M]): a constant-scalar,
  graph-gauge base for cellular-sheaf sampled
  interface-agreement evidence, with continuum watertightness retained as the
  ambitious successor theorem. General function-valued stalks and admitted
  trace/conversion restriction maps remain the target architecture rather than
  a property of this scalar base. Fallible `SheafComplex::from_charts` discovers
  interfaces by support overlap + shared zero-band sampling and returns an
  immutable `AdmittedSheafComplex`. Raw public `SheafComplex` parts remain
  available for fixtures and incidence diagnostics but can emit only
  `Unknown`/`NoClaim`; callers cannot construct or mutate admitted evidence.
  Geometry-derived
  seeds make the retained two-chart swapped-pair mismatch bound bitwise stable,
  but full evidence/provenance binds patch labels and the finite-iteration
  multi-patch gauge diagnostic is not yet permutation-invariant. Plus
  pairwise-interface 3-cliques as candidate 2-cells. Their three pairwise
  rejection-sample sets are independent; a clique and the retained minimum
  pairwise count do not prove a common triple zero-band point or aligned
  restriction samples and carry no Čech/topology authority. The fallible
  `SheafSkeleton::of` structurally validates caller-supplied public adjacency
  for diagnostic algebra and then omits those candidate triples; it does not
  authenticate builder origin or confer topology authority. Verified common
  intersections and an admitted extractor remain successor work. The edge-level
  `delta0_edges`/`delta1` maps
  assemble as fs-sparse matrices with entries in {−1, 0, +1}; their
  `δ¹·δ⁰ = 0` identity holds BITWISE (small-integer f64). The public `delta0`
  sample-row restriction incidence is deliberately a different cochain space.
  All three raw-complex incidence entry points and `section_solve` return
  `Result`, validate ordering/indices before sparse pushes, and enforce the same
  static patch/interface/sample/triple ceilings before sparse construction.
  `section_solve` fallibly reserves its graph workspaces and returns typed
  indeterminate-sample or numerical-overflow refusals rather than successful
  non-finite diagnostics. The current `fs-sparse::Coo` staging/assembly API
  still allocates internally through infallible vectors after these caps; a
  fully fallible sparse-builder successor remains required before
  `ResourceExhausted` can cover the entire incidence path. Triple
  discovery builds fallibly reserved vector adjacency and probes the smaller deterministic
  adjacency set explicitly, counts every membership probe (not only emitted
  triangles), polls cancellation at bounded strides and before publication, and
  enforces `SHEAF_MAX_TRIPLE_CANDIDATES` with a structured work refusal. The
  builder also preflights chart count, O(P²) pair candidates, worst-case chart
  evaluations, and retained interface samples before evaluating a chart;
  bounded support/domain/interface/sample/triple allocations return named
  resource refusals. Cancellation reports pair context plus an explicit typed
  progress unit (charts, pairs, draws, edges, probes, or retained triples).
  Every producer sample is checked immediately after evaluation; a non-finite
  signed distance is a pair-, chart-, and point-attributed build refusal rather
  than an implicit rejection-sampling miss. The legacy-named
  `AdmittedSheafComplex::watertightness(tol)` returns
  `Evidence<SheafVerdict>`: PASS requires a valid
  retained sampling scope and a uniquely ordered public complex with at least
  one nonempty interface and
  every sample's |mismatch| enclosure INSIDE `[0, tol]` via fs-ivl's
  sound predicates (no bound extraction); FAIL requires an enclosure
  ENTIRELY above tol — an interval-proven interface violation with the
  offending interface cells and magnitudes attached; anything else is an
  honest Unknown. A proven leak remains FAIL even when unrelated intervals are
  indeterminate; the optional gauge diagnostic is then unavailable. Positive
  sampled agreement does not establish between-sample coverage, continuum
  watertightness, cocycle membership, or non-exactness and therefore emits no
  H¹ or topological-obstruction claim.
  `AdmittedSheafComplex::mismatch_bounds()` exposes only immutable,
  context-free numeric enclosures. Tolerance predicates take the tolerance at
  use time, so a loose-tolerance boolean cannot be detached and reused as a
  strict-tolerance result; the bounds still are not a replay-complete source
  receipt.
  Unresolved unbounded overlaps refuse with pair attribution;
  `from_charts_clipped` admits an explicit finite local scope and preflights
  that clip even for empty/disjoint inputs. `SheafComplex::sampling_clip`
  retains the exact caller scope (`None` means admitted global supports), and
  sampled-agreement fingerprint input binds the global/local discriminator, all six
  clip endpoint bit patterns, complete interface sample points/enclosures,
  triples, assessed bounds, and the verdict payload. Those legacy canonical bytes
  are fed incrementally into FNV rather than duplicated in an unbounded transcript
  `String`. The v1 complex does not
  retain source chart identities, and its live `ProvenanceHash` is a 64-bit FNV
  fingerprint. It is neither collision-resistant content identity nor a
  replay-complete source-model receipt; migration to the shared strong identity
  substrate is tracked by the geometry-specific
  `frankensim-sj31i.52.4`. `section_solve` computes
  per-patch gauge offsets, pinning the smallest patch in every connected
  component (including isolates),
  over the adjacency Laplacian and reports the fractional reduction in
  uncentered sample-level midpoint-mismatch mean-square energy from that graph gauge. This
  diagnostic is not an
  exact/coexact/harmonic classifier; the feature-gated `sheaf_repair` module
  owns those claims. `validate_outside_ray_samples` is not an independent
  topology falsifier: with both endpoints proven outside each chart by an
  excluding support or rigorous positive distance enclosure,
  a boolean sign sequence necessarily has an even number of toggles. It is only
  a fallible, work-capped input/sampling diagnostic whose report retains work
  and transition counts: empty inputs, zero steps, non-finite
  endpoints or chart values, endpoints not nominally outside, endpoints whose
  outside status is unproven, and
  unrepresentable interpolation are structured refusals. Segment points use
  convex interpolation rather than a potentially overflowing endpoint
  subtraction, and producer-side cancellation is checked after every chart
  evaluation. It carries no promotion authority. Authentic independent
  cross-examination requires certified oriented surface intersections or
  winding/degree evidence and remains tracked work.

- `ident` module (the R3 AMENDMENT, bead lmp4.10): STABLE PERSISTENT
  ENTITY IDENTITY is a hard core requirement — `EntityId`s are assigned
  at creation and transformed EXPLICITLY by ledgered edits
  (`IdTransform`: Preserved/Replaced/Split/Merged/Created/Deleted;
  `IdentityMap::ops_touching` walks the full replace/split/merge
  ancestry). Identity is a kernel invariant, never a heuristic
  reconstruction. UNGATED: every new chart-producing operation must
  record its transforms.
- `diff` module ([F], behind the `semantic-diff` feature until its
  Gauntlet tier + kill metric are green): the PHYSICS diff.
  Fallible `semantic_diff` aligns worlds by `EntityId`, measures field
  differences on shared support (the sheaf band-sampling machinery),
  and attributes each finding to a RANKED list of contributing causal
  edits with per-edit contributions MEASURED across generation
  snapshots when supplied (unpartitioned-but-flagged otherwise).
  Unidentified entities degrade to a geometric fallback FLAGGED
  `attributed: false`, and the fallback fraction is the R3
  early-warning metric. `semantic_diff_clipped` provides explicit finite local
  scope and retains it in `DiffReport::sampling_clip`; unresolved unbounded
  comparisons otherwise return typed refusal. Invalid tolerances,
  non-representable sample coordinates, non-finite chart values, and overflow
  of a finite pair's difference are also typed refusals rather than false-clean
  reports; the consumer polls cancellation immediately after each producer
  evaluation and once more before publishing the final report.
  `DiffReport::filter` triages by
  region/quantity/magnitude.

- `sheaf_repair` module (patch Rev L, bead wqd.14; [M], behind the
  `sheaf-repair` feature until certifier trials pass): DIAGNOSIS → ALGEBRAIC
  GAUGE-CORRECTION PLANNING. The legacy `hodge_decompose` performs a deterministic,
  fixed-iteration sequence: least-squares fit to the coboundary image, then a
  least-squares fit of the residual to the retained triangle-coboundary image,
  followed by a remainder. The retained fixtures compare the first fit with an
  independent dense reference, but a generic result is not yet a certified
  orthogonal Hodge decomposition. `hodge_decompose_bounded` runs that diagnostic
  over an opaque admitted skeleton with explicit sweep/operator/memory/poll
  budgets and `Cx` cancellation, returning retained usage or a typed refusal.
  The parallel `assess_hodge_decomposition_bounded` path is the first
  numerical-authority slice: it pins the smallest patch in every connected
  component (including isolates), retains candidate-named projection vectors,
  and returns exactly `Converged`, `Indeterminate`, or `Refused`. Promotion to
  the opaque converged view requires every outward `fs-ivl` enclosure for the
  two normal equations, `delta0^T` remainder, pairwise orthogonality, and
  reconstruction to meet the caller's finite non-negative relative tolerance
  under the versioned Frobenius/L2 normalization. The operators and
  reconstruction are interval-evaluated from exact point intervals around the
  retained candidates; a norm of an already-rounded residual cannot mint the
  token. Absolute residual enclosures, normalized enclosures, interval replay
  vectors, exact admitted incidence and mismatch values, stopping reason,
  budget, and measured work remain attached. Its spectral field is explicitly
  `Unknown`: for the admitted unweighted `delta0^T delta0` graph Laplacian it
  records exact connected-component nullity bounds, one set-valued structural
  `[0,0]` cluster with that multiplicity, deterministic component gauge
  representatives (not eigenvectors), and the versioned positive
  `max(1, 2 * maximum_degree)` Gershgorin scale. No eigensolver output is
  consumed, so candidate clusters and the requested/covered ranges remain
  empty and every numerical mode remains unresolved until a source-bound
  `fs-spectral` coverage receipt is consumed. The exact-nullity statement does
  not extend to weighted, signed, or quotiented operators. `Converged`
  therefore proves only
  the named finite-dimensional residual obligations; it does not prove
  continuum coverage, chart realizability, topology, H1, non-exactness, repair
  feasibility, or merge/publication authority. The compatibility assessment
  still refuses cancellation transactionally without publishing an in-flight
  candidate. The additive
  `assess_hodge_decomposition_retaining_sweeps_bounded` lane pre-binds the
  admitted source and uses an accounted two-buffer projection schedule: its
  committed iterate is changed only at a complete sweep boundary, while work
  performed in a discarded partial sweep remains charged. After at least one
  exact sweep, runtime cancellation or an ambient deadline, poll, or cost
  refusal may return `Interrupted` with an opaque in-memory
  `SheafSweepCheckpoint`. Public checkpoint metadata is limited to source,
  tolerance, typed phase, component count, original budget, measured usage,
  ordered ambient-attempt history, and committed exact/coexact sweep counts.
  Solver vectors remain private; the token grants no candidate, residual,
  convergence, spectrum, topology, repair, or publication authority. Invalid
  input, admission failure, local cap failure, allocation failure, and
  numerical/arithmetic failure remain ordinary `Finished(Refused)` outcomes.
  The additive
  `resume_hodge_decomposition_retaining_sweeps_bounded` entry point linearly
  consumes that opaque token under a fresh `Cx` and accepts no replacement
  source, tolerance, or local budget. It reconstructs only deterministic
  scratch, continues unsealed exact/coexact sweeps from the last committed
  vectors, and restarts verification from a private completed split when one
  exists. Local sweep/operator/work counters remain cumulative across attempts,
  including discarded partial-sweep and reconstructed-scratch work, while each
  ambient `BudgetConsumption` remains a distinct ordered entry in a fixed
  allocation-free attempt history. Pre-cancellation, fresh ambient admission
  refusal, cumulative local-cap refusal, or attempt-cap-plus-one returns
  `NotAdmitted` with the exact owned checkpoint logically unchanged; only a
  successfully admitted attempt may move its source or solver state. The
  resumed terminal report must reproduce uninterrupted candidate/receipt bits
  while committing each projection sweep exactly once. This is a same-process
  rollback/resume slice only: checkpoints still have no clone, serialization,
  migration, fork, loaded-image, durable-ledger, or cross-ISA authority, so the
  full P7/G4 checkpoint program remains open.
  Its INTERPRETATION CONTRACT is: fitted exact
  component → a sampled-mismatch 0-cochain correction candidate
  ONLY when every per-patch offset fits that chart's declared error
  budget — a repair never silently distorts geometry); coexact → a circulation
  candidate over retained triangle cells, with converter, chart/model,
  junction, and sampling causes left as ranked hypotheses until provenance or
  intervention distinguishes them; harmonic → a retained candidate remainder and its interface
  support after the current deterministic fits. The three reported squared-norm
  ratios use the original mismatch norm as a common denominator; without a
  per-result orthogonality proof they are diagnostics and need not sum to one.
  The fixed-iteration numerical solves carry no per-result convergence
  certificate, so
  `plan_repair` does not generically certify non-exactness, impossibility of a
  gauge repair, or a required geometry-topology change. The ring conformance
  fixture separately checks a retained nonzero remainder for closure and
  non-exactness before labeling that fixture H¹.
  `plan_repair -> Result<RepairPlan, SheafRepairError>` validates raw incidence,
  exact mismatch/budget cardinalities, finite mismatch values, and finite
  non-negative per-patch budgets before decomposition or proposal allocation;
  refusal publishes no partial plan. It emits ranked agent-facing proposals.
  A finite post-repair seam norm is attached only to
  the constructive gauge proposal; diagnostic and reroute proposals use `+∞`
  to mean that no comparable post-state seam norm has been established.
  Optional Rep-Router reroutes retain modeled cost. `try_apply_gauge` validates
  raw incidence plus finite, exact-length mismatch/gauge cochains and returns a
  typed allocation/arithmetic refusal. The compatibility `apply_gauge` keeps
  the original `Vec` result for downstream merge callers: valid inputs preserve
  historical behavior, while malformed inputs return a fixed non-finite refusal
  sentinel instead of panicking; existing merge callers fail closed on that
  sentinel. New authority paths must use the typed API.
  Gauge application subtracts one explicit coboundary from the retained mismatch
  vector; it does not edit or re-evaluate a chart, commit geometry, or prove
  realization of that correction.
  Applying the same nonzero correction twice is not idempotent, while re-planning
  a converged algebraic mismatch can produce a zero follow-up gauge. Transactional
  chart mutation plus revalidation remains explicit successor work.

- `sheaf_merge` module (Proposal 10's CROWN JEWEL, bead lmp4.12; [M],
  behind the `sheaf-merge` feature — which enables `sheaf-repair` —
  until its Gauntlet tier + kill metric are green): the sheaf operators
  support a guarded merge heuristic, not an H¹ certifier.
  `three_way_merge` forms the union of edits (X + Y − B at the cochain
  level), runs deterministic fixed-iteration gauge reconciliation, and
  reports `Resolved` only when the resulting nominal `f64` infinity norm
  passes the requested tolerance. The v1 `MergeResidualReceipt` is only that
  local post-state numerical observation; it is not yet content-bound to the
  skeleton, parents, union, gauge, or algorithm and therefore is neither a
  portable replay receipt nor an interval watertightness certificate. When the
  check fails, a dominant decomposition remainder is retained as a
  `CandidateRemainderConflict` localized to its supra-tolerance cells with both
  caller-supplied parent labels (not authenticated identities); it proves neither closure,
  non-exactness, H¹, nor the absence of another local repair. Other
  numerical residues ESCALATE unresolved. A separately retained valid
  skeleton cochain may receive an H¹ interpretation only after executable
  closure and non-exactness evidence, as in the ring conformance fixture.
  Because the v1 call has no base assignment map and its success variants
  cannot return merged assignments, every assignment payload currently refuses
  before decomposition instead of being silently discarded. Even same-key,
  different branch values are only pairwise conflict candidates: either value
  may equal the unknown base. True base-aware assignment merge is tracked
  successor work. Trust is
  conditioned on `spectral_gap` (the thresholded second sorted
  weighted vertex-Laplacian eigenvalue) with LowGap flagging on `Resolved`,
  `Conflicted`, and `EscalatedUnresolved` numerical outcomes (R5).
  `candidate_remainder_conflict_rate` returns a typed error for invalid or
  refused trial sets and measures only the candidate-remainder numerator. It
  is diagnostic-only: without the full
  Resolved/Trivial/Conflicted/Escalated/Refused outcome histogram it cannot
  satisfy the 25% kill criterion, and it is not a certified cohomology rate.

## Invariants
1. Trait laws (G0, geo-001, 12k seeded queries): `inside(x)` ⇔ `sd(x) <
   0`; `support()` bounds the region (no negative sd outside, to
   tolerance); certified Lipschitz bounds hold along random steps;
   claimed gradients are unit-norm and match central differences on
   smooth fixtures.
2. Agreement soundness: identical valid presentations always agree; a
   disagreement implies at least one chart's geometry or error declaration is
   wrong — proven by detecting an undeclared 0.03 bias with exact-strength,
   gap-localized diagnostics naming the lying chart (geo-003). Missing or
   malformed evidence is structurally `Unknown`, including when diagnostic
   retention is disabled.
3. Conversion receipts are conservative: empirical |sampled − exact| over
   10k seeded points never exceeds the receipt's QoI bound (geo-004);
   receipts satisfy the `Certified` discipline (enclosure-grade, chained
   provenance).
4. Budget infeasibility refuses BEFORE any sampling runs, with ranked
   fixes (geo-004).
5. Agreement checks are seeded-deterministic (same config ⇒ identical
   JSON, G5) and poll cancellation every sample (geo-002/005).
6. Every sampling consumer validates raw supports before intersection/union
   and admits a finite `SamplingDomain` before evaluating charts. Unbounded
   geometry requires explicit finite scope; bounded disjoint pairs remain
   ordinary no-overlap rather than errors.
7. RD.1a admission is finite and type preserving: collection counts and
   aggregate rank are checked before sorting; dimension and canonical-byte
   ceilings are enforced before publication. Canonical order is independent of
   input collection order. Duplicate identity families stop at a fixed-order
   ambiguity barrier before per-record validation; every reference resolves
   within its nominal family, active local-model references resolve only to
   `Active` inequalities/contacts, local links bind real incidence edges with
   `dim(link) + 1 = codim`, and ordered
   inequality/contact semantics never admit complex coefficients. The typed
   BLAKE3 receipt binds all semantic fields and the immutable model version.
   Graded spaces are contiguous (including explicit zero-dimensional degrees),
   differentials cover every adjacent degree, and boundary/incidence/
   constitutive references cannot cross chart ownership.
8. RD.1b structural morphisms have one rank-neutral identity per exact admitted
   geometry. Strict composition is defined only across exact geometry and
   declared evidence-artifact/rank seams. Primitive morphism identities and
   no-equivalence artifacts remain in semantic order, making receipt identity
   independent of composition parenthesization but sensitive to factor order.
   Primitive chart-map receipts domain-separate their family and bind exact
   source chart, target chart, overlap, and forward-map IDs. Homogeneous chart
   paths additionally require an exact chart seam. Mixed strict/chart paths
   retain one ordered typed primitive sequence, and adjacent chart primitives
   preserve the same exact seam obligation inside a heterogeneous path.
   Primitive finite-complex refinement receipts additionally bind exact
   source/target complex and resolution IDs plus nominal prolongation and
   commutation IDs. Homogeneous refinement paths require exact complex and
   resolution seams; heterogeneous paths preserve that obligation for adjacent
   refinement primitives.
   Primitive whole-object inclusion receipts bind nominal declared-map and
   claimed-containment IDs under a separate family tag. Homogeneous inclusion
   paths require the exact geometry seam already enforced for every arrow;
   composition retains ordered factor declarations but proves no containment
   theorem.
9. Standalone stratum morphisms use a separate schema and exact
   `(geometry, stratification, stratum)` objects. Identity is available only on
   one exact triple. A component binds nominal map and constructibility IDs plus
   a nonzero no-authority artifact. Composition requires equality of the whole
   middle triple and flattens ordered primitive and no-authority lineage. These
   tokens cannot enter whole-geometry heterogeneous composition or carry
   geometry-wide evidence transport.
10. Exhaustive stratified-map assembly candidates bind exact endpoint geometry
   and stratification IDs, exactly one direct sealed component per source
   stratum in canonical source order, explicit source/target stratum selectors,
   nominal assembly/global-constructibility IDs, and a nonzero no-authority
   artifact. Repeated target strata are allowed; target coverage and every
   global map theorem remain unclaimed. Ordered two-step exhaustive-map
   candidates bind two exact recursively typed children, derive all six
   source/middle/target selectors, require both exact middle seams, and retain
   one nominal composition declaration plus a nonzero no-authority artifact.
   They neither rescan nor synthesize a composed component family or global map.
11. Standalone declared spans bind exact `source <- apex -> target` geometry
   identities, exact ordered admitted leg-morphism IDs, the common-apex/outer-
   endpoint orientation of both sealed legs, and one nonzero no-authority
   artifact. They do not enter directed-morphism composition or expose a direct
   evidence transport.
12. Fixed-resolution quasi-isomorphism candidates bind an exact nonempty
   homogeneous refinement path; exact endpoint geometry, local-model, role,
   complex, and resolution selectors; each endpoint's exact fixed-resolution
   scope witness; nominal theorem/checker/check-receipt IDs; and an explicit
   no-authority artifact. They do not enter morphism composition or the
   equivalence boundary, and no field can be omitted without changing or
   refusing the candidate receipt.
13. Local-presentation correspondence candidates bind exact endpoint geometry
   and local-model selectors; separate canonical sets for equality, active-
   inequality, active-contact, and constitutive relation edges; one nominal
   aggregate declaration; and one nonzero no-authority artifact. Every selected
   source and target family member has at least one edge, exact duplicate pairs
   refuse, and caller ordering is nonsemantic. Repeated sources and targets are
   admitted explicitly, so coverage cannot be mistaken for functionality or a
   bijection.
14. Scoped presentation-equivalence candidate assemblies bind exact sealed
   child IDs for the tangent, cotangent, and deformation-obstruction fixed-
   resolution quasi-isomorphism candidate roles plus one exhaustive local-
   presentation correspondence candidate. Every child must retain the exact
   declared geometry and local-model endpoints. The three role candidates must
   retain identical source/target resolution IDs and fixed-scope witnesses;
   those four selectors are derived from sealed children, not supplied by the
   caller. The assembly requires a nonzero no-authority artifact and grants no
   equivalence, inverse, homotopy, composition, or evidence-transport API.
15. Finite stratification-refinement candidates bind exact refined/coarse
   selectors to one exact exhaustive child, cover both finite stratum families
   in their declared directions, permit repeated coarse targets, and enforce
   dimension monotonicity. Ordered two-step candidates require exact equality
   of both middle selectors and do not synthesize a direct refinement.
16. Direct chart-transition inverse-law candidates accept only two direct
   declared chart maps with exact reversed geometry/chart endpoints, one common
   overlap, one evidence variance, and both exact ordered evidence seams.
   Parallel-path candidates independently require exact common endpoints and
   preserve caller-significant left/right order; neither candidate proves its
   nominal relation.
17. Fixed-foot span-morphism candidates bind exact source/target spans, one
   exact apex morphism with the required orientation, all four exact parent
   legs, both recomposed apex-to-foot routes, and exact ordered comparison
   children. No commuting triangle or 2-cell is inferred.
18. Span pullback-square candidates bind exact composable parent spans,
   projections from one exact proposed apex, and exact ordered routes to the
   common middle geometry. Span-composition candidates bind one exact square,
   one exact proposed outer span, and both exact recomposed outer legs. No
   pullback or composition authority is inferred.
19. Generic inverse-law candidates bind exact forward/reverse morphisms to
   reversed endpoints, derive both canonical identities and round trips, and
   require exact ordered identity-versus-round-trip comparison children. They
   expose no caller-authored identity slot and mint no inverse.
20. Span-equivalence candidates bind exact source/target spans, ordered
   forward/reverse span-morphism candidates, and the exact apex inverse-law
   candidate. All direct IDs and transitive outer-foot/span/apex/morphism seams
   are rechecked; the mandatory no-authority artifact cannot be converted into
   an equivalence boundary.
21. Evidence-polarity transport candidates bind one exact recursively typed
   nonidentity morphism child, derive every endpoint/variance/artifact/rank
   selector from it, and require exact equality of the explicit input/output
   polarity tags. Identity children and all six unequal pairs refuse. The
   nominal declaration and mandatory no-authority artifact grant no canonical
   classification, payload authenticity, preservation theorem, inverse,
   equivalence, or physical authority.
22. Span evidence-transport candidates bind one exact recursively typed span,
   one exact recursively typed left restriction-polarity child, and one exact
   recursively typed right balance-corestriction-polarity child. Admission
   rechecks both span-leg morphism IDs, all four leg endpoints, both required
   variances, and exact evidence-artifact/rank/polarity continuity at the apex.
   The nominal aggregate declaration and mandatory no-authority artifact do not
   make the span functional or grant executable pull-push, base-change,
   preservation, theorem, inverse, equivalence, or physical authority.

## Error model
Structured teaching values throughout: `ConvertDiag` (ranked fixes),
`fs_exec::Cancelled` for interrupted checks, `AgreementStatus::Unknown` plus
structured `AgreementUnknownReason` for non-evaluable checks, and
`fs_evidence::CertifyError` through receipts. Router execution additionally
returns `ExecuteError` with a stable missing-runner/runner/oracle-record class;
oracle evidence refusals use `CostOracleError`. Fallible sampling and conversion
entry points return structured refusals, and `Aabb::new` is total and normalizes.
The feature-gated legacy raw `hodge_decompose` and merge diagnostic APIs still
assert some shape/index preconditions on caller-built skeletons; their total,
budgeted `Result` replacements remain required before promotion. `plan_repair`
and `try_apply_gauge` are structured refusals; the source-compatible
`apply_gauge` adapter is total but carries only a non-finite refusal sentinel,
not the structured cause. The base
`SheafComplex` incidence and section APIs are now structured refusals for
malformed/oversized raw parts. RD.1a uses `DerivedAdmissionReportV1`:
unsupported schema/category/scope/encoding, typed-reference defects, explicit
cross-chart references, mixed units/frames, finite-complex and stratification
defects, resource exhaustion, cancellation, and canonical identity failures
publish no admitted token.
RD.1b uses `DerivedMorphismErrorV1` for schema/endpoint/convention defects,
model-version drift, chart and complex ownership/shape defects, structural rank-
envelope or truncation regression, evidence orientation, missing nominal
artifact IDs (including inclusion declarations), declared-rank strengthening,
equivalence laundering, typed chart/refinement/evidence seams, bounded lineage,
allocation, cancellation, and canonical identity failures.
Refusal publishes no admitted morphism.
Evidence-polarity packet admission uses
`DerivedEvidencePolarityTransportCandidateErrorV1` for schema, zero nominal
identity, raw/sealed morphism-child mismatch, unsupported identity evidence,
unequal input/output polarity, cancellation, and canonical-identity defects.
Refusal publishes no polarity packet and cannot mutate or upgrade the supplied
sealed morphism.
Standalone stratum-scoped admission and composition use
`DerivedStratumMorphismErrorV1` for exact geometry/stratification/stratum
ownership, model-version and convention compatibility, nominal identity and
authority laundering, exact middle-object seams, bounded lineage, allocation,
cancellation, and canonical identity failures. Refusal publishes no stratum
morphism and cannot affect a whole-geometry morphism.
Exhaustive finite assembly admission uses
`DerivedExhaustiveStratifiedMapCandidateErrorV1` for exact endpoint and
stratification binding, missing nominal identities, full canonical source-
stratum coverage, raw/sealed component identity and endpoint seams, direct-
component shape, bounded retention, allocation, cancellation, and canonical
identity failures. Refusal publishes no assembly candidate and grants no global
map authority.
Ordered exhaustive-map composition admission uses
`DerivedExhaustiveStratifiedMapCompositionCandidateErrorV1` for unsupported
schema, zero child/declaration/no-authority identities, raw/sealed child drift,
both exact middle seams, cancellation, and canonical-identity defects. Refusal
publishes no partial candidate, component composition, or direct map.
Finite-refinement candidate admission uses
`DerivedStratificationRefinementCandidateErrorV1` for raw/sealed selector and
child mismatches, two-sided coverage, dimension increase, resource/allocation,
cancellation, and identity defects. Ordered refinement composition uses
`DerivedStratificationRefinementCompositionCandidateErrorV1` for raw/sealed
children and both exact middle seams. Neither refusal publishes a partial
candidate or direct refinement.
Standalone span admission uses `DerivedSpanCorrespondenceErrorV1` for schema,
zero no-authority identity, raw-leg/sealed-leg mismatch, leg orientation,
cancellation, and canonical-identity defects. Refusal publishes no span token.
Span evidence-transport admission uses
`DerivedSpanEvidenceTransportCandidateErrorV1` for unsupported schema, zero
child/declaration/no-authority identities, raw/sealed child drift, exact
span-leg binding, all four endpoint relations, required left restriction and
right balance-corestriction variance, exact apex evidence/rank/polarity seams,
cancellation, and canonical-identity defects. Refusal publishes no partial
route and cannot upgrade the span or either polarity child.
Parallel-path, span-morphism, pullback-square, and span-composition admission use
`DerivedParallelMorphismComparisonCandidateErrorV1`,
`DerivedSpanMorphismCandidateErrorV1`,
`DerivedSpanPullbackSquareCandidateErrorV1`, and
`DerivedSpanCompositionCandidateErrorV1`, respectively. They distinguish
raw/sealed child mismatch from endpoint/orientation, exact parent-leg or
projection binding, recomposition refusal, ordered comparison-route drift, and
proposed outer-span drift; no refusal publishes partial categorical authority.
Direct chart and generic inverse-law admission use
`DerivedChartTransitionInverseLawCandidateErrorV1` and
`DerivedMorphismInverseLawCandidateErrorV1` for direct-child shape,
reversed-endpoint, overlap/variance/evidence-seam, derived identity/round-trip,
and ordered comparison defects. Span-equivalence assembly uses
`DerivedSpanEquivalenceCandidateErrorV1` for raw/sealed child mismatch,
ordered outer-foot/span-morphism binding, exact apex inverse-law binding,
cancellation, and canonical-identity defects. All three refuse without
upgrading any supplied child.
Fixed-resolution quasi-isomorphism candidate admission uses
`DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1` for schema, exact
endpoint/path/selector, local-model ownership, role/complex/chart/locality,
fixed-resolution scope, nominal identity, cancellation, and canonical-identity
defects. Refusal publishes no candidate token and cannot affect the supplied
sealed morphism.
Local-presentation correspondence candidate admission uses
`DerivedLocalPresentationCorrespondenceCandidateErrorV1` for schema, endpoint
and local-model ownership, geometry conventions, exact chart/locality,
relation-member ownership, two-sided family coverage, duplicate source/target
pairs, nominal identities, aggregate resource/allocation limits, cancellation,
and canonical-identity defects. Refusal publishes no relation candidate.
Scoped presentation-equivalence candidate assembly uses
`DerivedScopedPresentationEquivalenceCandidateAssemblyErrorV1` for schema,
nonzero local-model/no-authority identities, raw/sealed child identity, exact
role slots, endpoint and local-model agreement, common resolution/scope-
witness selectors, cancellation, and canonical-identity defects. Refusal
publishes no assembly and cannot promote or mutate any supplied sealed child.

## Determinism class
Deterministic: seeded sampling, insertion-ordered charts, canonical JSON
renderings; no clocks, no addresses. Float behavior inherits
fs-math-class scalar arithmetic. RD.1a sorts every set-valued object collection
by nominal identity, uses an unambiguous length-framed nested encoding and the
shared schema-typed BLAKE3 canonical encoder, stops duplicate identity families
at a deterministic ambiguity barrier, and retains exact f64 unit/contact
coefficient bits. Negative-zero Coulomb coefficients are refused as
noncanonical input. Reordering input collections does not move identity;
changing category, coefficients, units, frame, model version, local singular class,
stratification, or proof metadata does.
RD.1b encodes exact endpoints, map witness/class, structural evidence transport,
ordered no-equivalence artifacts, and ordered primitive identities. Flattened
composition makes equal ordered factor sequences replay to the same receipt
regardless of parenthesization; reversing factors changes identity.
Evidence-polarity packets use a separate schema/domain and recursively bind the
unchanged six-field morphism-child schema. Changing either polarity tag,
nominal declaration, no-authority ID, child identity, or child-derived evidence
metadata changes the 12-field parent receipt; no existing morphism identity
material is reinterpreted.
Declared chart-map primitive class bytes additionally encode the source/target
chart, nominal overlap, and nominal map IDs under a distinct family tag.
Heterogeneous composites use a distinct class tag while the ordered primitive
factor IDs bind every family-specific nominal ID and structural field without
encoding parenthesized intermediate operands or authenticating artifact payloads.
Finite-complex refinement primitives use a distinct tag and bind source/target
complex, resolution, prolongation, and commutation IDs; homogeneous composites
use another distinct tag while flattened factors preserve parenthesization-
independent identity.
Whole-object inclusion primitives likewise use a distinct tag and bind exact
nominal map and containment-artifact IDs; their homogeneous composite tag plus
flattened primitive factors preserves family, order, and parenthesization
without authenticating either artifact.
Standalone stratum-morphism receipts use a separate schema/domain. They encode
the exact source and target geometry/stratification/stratum triples, identity or
primitive/composite class tag, ordered no-authority artifacts, and flattened
primitive receipt IDs. Primitive class bytes bind the nominal component-map and
constructibility IDs. Equal ordered paths replay independently of
parenthesization; moving any selector or reversing factors changes identity.
Exhaustive assembly-candidate receipts use another separate schema/domain and
encode exact endpoint geometry/stratification IDs, canonically source-ordered
sealed component IDs (which transitively bind each explicit stratum selector),
nominal assembly/global-constructibility IDs, and the no-authority artifact.
Replaying the same full coverage is stable; moving any retained field changes
identity without authenticating components or nominal payloads.
Ordered exhaustive-map composition receipts use a further separate schema/domain
and bind six selectors derived from the sealed children, first then second typed
child IDs, one nominal composition declaration, and one no-authority artifact.
The ten-field parent plus both complete eight-field child schemas yields an exact
recursive field limit of 26. Replay is stable; moving any derived selector,
child, order, declaration, or no-authority identity changes the receipt.
Standalone span receipts use a separate schema/domain and encode exact
source/apex/target geometry IDs, left then right admitted-leg IDs, and the
no-authority artifact. Replaying the same ordered legs is stable; swapping valid
distinct equal-endpoint legs changes identity.
Span evidence-transport receipts use another separate schema/domain. Their 17
parent fields encode the derived source/apex/target selectors, outer/apex
evidence artifacts, ranks and polarities, the typed span child, typed left then
right polarity children, nominal aggregate declaration, and no-authority
artifact. Complete recursive validation of the six-field span child and both
18-field polarity-child trees yields the exact 59-field limit. Replay is stable;
moving any retained field, child, child domain, or child order changes or
refuses the receipt without authenticating the declared route.
Recent structural-candidate receipts each use a separate schema/domain and
recursively bind complete typed child descriptors. Their exact parent/recursive
field counts are: ordered exhaustive-map composition 10/26, stratification
refinement 7/15, ordered refinement composition 10/40, direct chart-transition
inverse law 10/22, parallel-path comparison 7/19,
span evidence transport 17/59, span morphism 11/67, pullback square 13/56,
span composition 10/84, generic inverse law 8/58, and span equivalence 11/215.
Each receipt encodes selectors
derived from or exact-validated against sealed children and admitted endpoints,
plus ordered direct child IDs and nominal/no-authority fields; changing
child order, a derived selector, or any retained declaration changes identity.
Shared child IDs are counted again at each typed edge, and no parent receipt
authenticates a nominal relation merely by retaining it.
Fixed-resolution quasi-isomorphism candidate receipts use another separate
schema/domain and encode the exact refinement path; ordered endpoint geometry,
local-model, complex, and resolution IDs; the selected complex role; both
endpoint scope witnesses; nominal theorem/checker/check-receipt IDs; and the
no-authority artifact. Exact replay is stable, while moving any retained field
changes identity without authenticating its payload.
Local-presentation correspondence candidate receipts use a separate domain and
encode exact endpoint geometry/local-model IDs; four independently sorted
canonical sets of complete typed relation edges; the nominal aggregate ID; and
the no-authority artifact. Permuting a relation collection is identity-neutral;
moving a member, edge artifact, aggregate, endpoint, or authority field changes
identity without authenticating relation semantics.
Scoped presentation-equivalence candidate assembly receipts use another
separate domain and exactly 13 ordered fields: geometry and local-model
endpoints; source/target resolution IDs and scope witnesses derived from sealed
children; the tangent, cotangent, deformation-obstruction, and correspondence
candidate receipts as recursively schema-bound typed children; and the
no-authority artifact. Exact replay is stable, child schema/version/context is
part of the parent schema identity, and moving any one field changes identity
without promoting any child claim.

## Cancellation behavior
Chart evaluation and production sampling paths take `&Cx`.
`check_agreement` polls per sample point; conversion grids, sheaf interface
draws, triple discovery, the ray-sequence validator, and semantic-diff field
draws poll `Cx` directly at deterministic bounded strides and return typed
cancellation diagnostics without publishing partial authoritative output.
Incidence assembly and mismatch assessment/section solve do not yet accept
`Cx`. `hodge_decompose_bounded` does poll one caller context under explicit
work/memory limits, but high-level proposal construction and merge diagnostics
remain without `Cx` or complete fallible-allocation accounting, so these APIs
are not P7-complete. Chart-local polls remain an additional inner-kernel
obligation. RD.1a admission takes `&Cx`, checks it
before sorting, at a fixed stride while preflighting nested graded spaces, after
every bounded canonical sort, once per duplicate identity family, once per
chart/constraint/complex/local-model/stratum/incidence/link, during nested
canonical item encoding, inside the
streaming identity encoder, and immediately before publication. Cancellation
returns a stage and completed-work count and cannot expose a partial admitted
object.
RD.1b admission polls at entry, before canonical identity construction, inside
finite-complex refinement rank-envelope scans, inside the streaming encoder, and
immediately before publication. Composition polls at entry, at a fixed stride
while copying bounded typed-primitive, factor, and no-claim lineage, inside
identity construction, and before publication. Cancellation exposes no partial
admitted morphism.
Standalone stratum admission polls at entry, at a fixed stride while resolving
finite stratum selectors, before and inside canonical identity construction,
and immediately before publication. Composition additionally polls while
validating and copying bounded scoped primitive/factor/no-authority lineage.
Cancellation exposes no partial stratum morphism.
Exhaustive assembly-candidate admission polls at entry, at a fixed stride while
checking and retaining exact component coverage, before and inside canonical
identity construction, and immediately before publication. Cancellation exposes
no partial assembly candidate.
Ordered exhaustive-map composition polls at entry, inside canonical identity
encoding, and immediately before publication. It scans and allocates no child
collection, and cancellation exposes no partial composition candidate.
Standalone span admission polls at entry, before and inside identity encoding,
and immediately before publication. Cancellation exposes no partial span token.
Span evidence-transport admission polls at entry, inside its bounded canonical
encoder, and immediately before publication. It scans and allocates no child
collection, and cancellation exposes no partial route candidate.
Finite-refinement admission additionally polls at a fixed stride while scanning
canonical two-sided stratum coverage and uses fallible bounded coarse-coverage
storage; ordered refinement composition polls at entry, inside identity
encoding, and before publication. Span morphism propagates cancellation from
both target-route compositions, pullback-square from both middle-route
compositions, span composition from both outer-route compositions, and generic
inverse law from both canonical identity constructions and both round-trip
compositions. Direct chart inverse-law checks exactly two ordered evidence-cycle
compositions rather than an unbounded loop. Parallel and span-equivalence
assembly allocate no child collection. Every candidate polls at entry, inside
its bounded canonical encoder, and before publication; cancellation exposes no
partial comparison, refinement, polarity, inverse-law, pullback/composition, or
equivalence token.
Fixed-resolution quasi-isomorphism candidate admission additionally polls while
scanning the bounded typed refinement path, before and inside identity encoding,
and immediately before publication. Cancellation exposes no partial candidate
token.
Local-presentation correspondence candidate admission polls at entry; before
and after each canonical relation sort; at fixed strides during membership,
two-sided coverage, and relation encoding; inside canonical identity encoding;
and immediately before publication. Cancellation exposes no partial relation
candidate.
Scoped presentation-equivalence candidate assembly polls at entry, inside its
bounded 13-field identity encoder, and immediately before publication. It
allocates no child collection and cancellation exposes no partial assembly.
RD.X1 statement admission polls before validation, once for every derived
theorem-lattice node, before identity construction, and inside the streaming
encoder. Its falsifier set and truncation lattice have hard versioned caps;
cancellation cannot expose a partial `ValidatedExitPathFamilyV1`.

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
All OFF by default per the Ambition-Tag rule (the default-path chart
abstractions remain unflagged `[S]`):
- `derived-geometry` [M] — RD.1a finite admitted derived/stratified object
  language, the RD.1b structural morphism spine, and the RD.X1 exit-path
  theorem-statement lattice; disabled until
  RD.1b/RD.1c/RD.X2 theorem, equivalence, artifact, and independent-checker
  lanes establish promotion evidence.
- `semantic-diff` [F] — semantic design diff; disabled until its
  Gauntlet tier + kill metric (R3 fallback fraction) are green.
- `sheaf-repair` [M] — sheaf-adjudicated repair; disabled until
  certifier trials pass (milestone P6).
- `sheaf-merge` [M] — the sheaf-adjudicated merge (crown jewel);
  implies `sheaf-repair`.
Each gates its own integration target (required-features declared).

## Conformance tests
tests/conformance.rs, cases geo-001..geo-005 (JSON-line verdicts; seeded
cases carry seeds): the fixture trait-law battery, multi-chart agreement
within composed bounds + G5 replay, lying-chart detection with localized
diagnostics, rigorous conversion receipts + teaching refusals, and
cancellation. In-module suites cover Aabb/vector laws, fixture known values,
agreement determinism, cancellation, zero-evidence/one-chart refusal,
non-finite configuration and chart output, `NoClaim`, malformed certificates
and support, and exact disagreement with zero diagnostic retention. The 30-case
library battery additionally covers router/brute-force Pareto agreement,
directed additive/relative exact-bound composition, sealed plan/oracle/spec
identity, read/write failure, retrospective error maxima, conditional
refusals, winner policy, cancellation, identity routes, and registry/path/label
limit+1 refusals.
`tests/sheaf.rs` additionally locks admitted-vs-raw authority, tolerance-at-use
mismatch predicates, malformed/oversized raw-algebra refusal, bitwise δδ,
component-rooted graph-gauge fitting, exact dense-overlap preflight with zero
chart evaluations, endpoint outside-enclosure proof, finite/non-finite producer
behavior, and cancellation progress units. Allocation-fault injection and
small-limit exhaustive cap testing remain explicit successor coverage.
`tests/metamorphic.rs` applies the shared G3 refinement-monotonicity engine to
the production `SphereChart::convert::<SampledSdf>` edge. Across a bounded
catalog of feasible exact-distance spheres, it tightens the requested error by
a nonidentity power of two and requires the certified receipt QoI not to
increase and the dense grid resolution not to decrease. The fixed conversion,
cancellation, and refusal pins in `tests/conformance.rs` remain authoritative
for their individual semantics.
`tests/derived.rs` is the RD.1a G0/G3/G5 battery: regular linkage admission and
schema replay, canonical redundant-equation ordering, distinct cusp/node local
models, relative boundary/contact-corner/local-link typing, active-set refusal,
mixed-unit and real/complex category mutations, admitted restricted-analytic
programs versus opaque callback refusal, unbounded/infinite model refusal,
complex-role and incidence/link mutations, resource caps, proof-scope refusal,
cross-chart reference refusal, permutation-stable ambiguous-duplicate reports,
permutation-stable fail-fast rank refusal, negative-zero contact and unit-scale
refusal, declared basis-dimension and local-link overflow refusal, model-version
identity movement, and pre-publication cancellation. The `derived` module's
private unit suite injects
mid-canonicalization cancellation deterministically without timing races.
The `derived_morphism` private G0/G3/G5 suite covers unique neutral identities,
associative ordered receipt replay, factor/no-claim order, exact middle-object
refusal, declared evidence-seam and mixed-variance refusal, declared-rank
strengthening refusal, both variance directions, equivalence laundering,
model-version/convention mismatch, the lineage cap/cap+1 boundary, deterministic
replay, and already-requested entry cancellation. Public-wrapper coverage now
uses fully admitted RD.1a fixtures for strict admission, identity, composition
neutrality, replay, ordered raw endpoint refusal, and immutable model-version
refusal. A private deterministic probe seam now replays cancellation twice at
every base-morphism identity-entry/streaming-encoder poll and exactly at the
pre-publication linearization point. A non-cancelling probe must reproduce the
public admission byte-for-byte, and every injected cancellation returns no
admitted morphism; the seam is not public API and changes no receipt schema.
The evidence-polarity candidate suite covers its separate domain and exact
12-parent/18-recursive-field contract, both variance lanes, all three equal
polarity pairs with identity separation, recursive child-domain refusal, raw/
sealed child drift, zero nominal identities, identity-evidence refusal, all six
unequal-polarity pairs, replay/accessors, retained-field identity movement, and
already-requested entry/identity-entry cancellation. The same private probe
pattern replays every evidence-polarity identity checkpoint twice and injects
exactly before publication; the non-cancelling path must equal public admission
and every injected path returns no candidate.
The standalone stratum suite covers domain-separated class encoding, exact
geometry/stratification/stratum ownership, identity neutrality, associative
flattening, exact middle-triple refusal, nominal-field identity movement,
authority-laundering refusal, deterministic replay, and entry cancellation.
The exhaustive assembly-candidate suite covers separate schema identity,
replay/accessors, nominal-field identity movement, exact full source coverage,
canonical order, raw/sealed ID binding, repeated target acceptance, identity-
component acceptance, composite-path refusal, missing authority fields, the
component cap, and entry cancellation.
The ordered exhaustive-map composition suite covers its separate ten-parent/
26-recursive-field schema, deterministic replay and every accessor, movement of
all four caller fields and all six child-derived selectors, typed child-domain
binding, child order, zero and raw/sealed identity refusals, both middle seams,
reversed order, and entry cancellation without claiming or executing a composed
map. The fixture uses disjoint complex, resolution, and local-model seed ranges.
The local-presentation correspondence-candidate suite covers domain-separated
replay and accessors for all four typed families, caller-order canonicalization,
many-to-many source/target acceptance, two-sided missing-coverage refusal,
top-level-but-out-of-model member refusal, duplicate-pair and zero-artifact
refusal, exact endpoint/model/version/locality binding, nominal semantic and
physical non-authority, aggregate capacity, every receipt field, and entry
cancellation. Exact-scope cases also cover every geometry convention and full
admitted-chart semantics rather than trusting a reused nominal chart ID.
The scoped presentation-equivalence candidate assembly suite covers its
separate 13-field schema/domain, deterministic replay and every accessor,
identity movement for every retained field, exact raw/sealed binding for all
four children, duplicate/wrong role refusal, raw and child endpoint/model
seams, common resolution/scope-witness refusal, the nonzero no-authority
selector, and already-requested entry cancellation. Its fixture supplies all
three roles
under one fixed-resolution presentation scope without asserting any promoted
equivalence theorem.
It also covers declared chart-map ownership, missing IDs, dimension/frame/unit/
quantity/scale mutations, typed-ID-bound receipt replay and public primitive
retention, homogeneous associativity, exact chart seams, identity neutrality,
mixed-family associativity, typed primitive order, and chart-seam retention
through heterogeneous paths.
Finite-complex refinement coverage adds exact selector/resolution ownership,
role/chart/rank/quantity/truncation/progress refusals, nominal-artifact receipt
movement, homogeneous associativity and identity, exact refinement seams, and
seam retention across heterogeneous parenthesization.
Whole-object inclusion coverage adds frozen family tags, nominal map and
containment receipt movement, strict-family domain separation, missing-ID and
authority-laundering refusals, inherited model-version compatibility, exact
identity/endpoint laws, homogeneous associativity, and mixed-family typed order.
Standalone span coverage adds deterministic replay and exact accessors, no-claim
identity movement, left/right order sensitivity, raw-leg ID mismatch, all four
apex/outer-endpoint orientation refusals, identity-left graph shape, and already-
requested entry cancellation. No test claims that a graph-shaped span is
functional or that arbitrary spans compose.
Span evidence-transport coverage adds its separate 17-parent/59-recursive-field
contract, exact typed-child domains and order, deterministic replay and every
accessor, identity movement for every retained and derived field, all three
polarity classes, zero-ID and raw/sealed child drift, both span-leg bindings,
all four endpoint relations, both reversed variances, exact apex evidence/rank/
polarity seams, and already-requested entry cancellation. These tests bind only
the structural restriction-to-corestriction route; they neither execute nor
authenticate an outer evidence transport.
Finite stratification-refinement coverage adds exact recursive child schemas,
two-sided coverage, repeated coarse targets, dimension monotonicity, every raw
and sealed selector seam, reversed orientation, deterministic replay, and entry
cancellation. Its ordered composition coverage adds both exact middle seams,
isolated child order, raw/sealed mismatch, every receipt field, and explicit
no-transitivity authority.
Direct chart-transition inverse-law coverage exercises both covariant and
contravariant closed evidence cycles, direct-child and reversed endpoint/chart
binding, common overlap, mixed variance, both ordered artifact/rank seam
failures, deterministic receipt replay, and cancellation. Parallel-path
coverage exercises identity-versus-cycle and direct-versus-composite pairs,
isolated left/right order, both endpoint mismatches and precedence, raw/sealed
children, every receipt field, and cancellation without claiming equality.
The layered span suite covers fixed-foot span morphisms, proposed pullback
squares, proposed outer composition spans, generic inverse-law packets, and the
final structural span-equivalence assembly. Genuine fixtures recompose every
retained route. Adversarial cases cover typed-child domains, zero and nonzero
raw/sealed drift, representative direct outer-foot/span/apex/projection/leg
bindings, transitive aggregate seams, evidence seams, reversed route and child
order, stale or forged transitive tokens,
derived identity/round-trip mismatch, exact recursive receipt fields, refusal
precedence, replay, and entry cancellation. These tests establish fail-closed
structural wiring only; they do not test or mint categorical laws.
`tests/exit_path.rs` supplies RD.X1 G0/G3 examples and a bounded-cancellation
regression: regular-cell poset sufficiency, cone/cusp groupoid-enriched
one-category fallback, circular-stratum local systems, finite-versus-full
higher fidelity, weak path-equivalence and monodromy
boundaries, hypothesis deletion, node-scoped refutation recording with richer fallback,
admitted-subject cross-binding, identity movement across every encoded family,
canonical and mandatory falsifiers, zero-identity refusal, exact referenced-
artifact budget boundaries, schema/truncation/implication/falsifier caps, and
already-requested cancellation. Mid-lattice cancellation storms and cross-
thread replay remain RD.X2/batch-verification work; this test target does not
claim those stronger G4/G5 results.

## No-claim boundaries
- NO watertightness/manifoldness/self-intersection certificates here —
  those are wqd.7 (validity certificates) and the sheaf bead; agreement
  checking is SAMPLED evidence, not a proof (`Agreed` means "no
  counterexample found at these seeded points under the reported certificate
  strength + declared intervals + configured slack").
- Two registered presentations are sufficient to run the pairwise check; this
  module does not certify implementation or provenance independence between
  them. Independence must be established by the campaign that cites the result.
- Explicitly clipped agreement, sheaf, semantic-diff, and conversion evidence
  is local to the stated clip. A clipped conversion receipt is plain
  `Evidence<SampledSdf>` with `NoClaim` abstract-distance authority: its finite
  `nominal_field_bound` describes reconstruction of the sampled
  `max(source_field, exact_box_sdf)` composite only and cannot be promoted to
  `Certified`.
- `SampledSdf` claims no Lipschitz bound for its interpolant and no
  gradients (rep-sdf's job); its outside-box enclosure relies on the
  source theorem recorded during conversion. For weak sources the same formula
  is only estimate/no-claim evidence, never a rigorous enclosure.
- The G3 refinement relation covers only the declared finite sphere/budget
  catalog and receipt monotonicity for this dense converter. It does not prove
  strict improvement, a convergence order, empirical field accuracy, adaptive
  representation behavior, a common sampling domain (converter padding depends
  on the requested budget), or monotonicity for arbitrary charts and budgets.
- `TraceStepClaim::LipschitzImplicit` certifies no-tunneling step radii, not
  Euclidean proximity from a small normalized residual. Consumers must retain
  that distinction in hit/error language. A short opposite-sign bracket can
  separately prove boundary existence; same-sign or indeterminate endpoint
  evidence, including a generic tangency, cannot.
- `ConverterSpec::certified` is a declaration, not an authenticated admission
  receipt. Runtime certified runners must return `Certified<f64>` local-error
  evidence and fs-ir rejects routes containing an estimated declaration, but a
  malicious caller can still lie in the Boolean until the opaque admitted-color
  / converter-authority migration lands.
- The dense converter deliberately uses the conservative largest actual full
  cell diagonal instead of a sharper trilinear-weight distance. Tighter
  geometry-dependent constants and interval-verified source sampling arrive
  with rep-sdf.
- Curvature, closest-point, ray-intersection, and integral queries are
  declared in the plan but NOT in this trait yet — added as capability
  traits when their first consumers land (router capability negotiation).
- `topology_hint` is a HINT; nothing verifies it (persistence
  certificates are wqd.19's).
- Cost models, chart selection, and the Pareto routing plane are the Rep
  Router bead's; `Region::primary` is insertion-order only.

## No-claim boundaries (derived geometry)

- `AdmittedDerivedGeometryV1` proves only bounded structural well-formedness,
  nominal typing, deterministic canonical identity, and the stated finite
  admission predicates. It does not prove a derived intersection, exactness or
  quasi-isomorphism of a complex, virtual-dimension theorem, smoothability,
  Whitney A/B or Thom conditions, link topology, transversality, rigidity,
  contact well-posedness, constitutive admissibility, or physical validity.
- `DerivedProofStateV1::ExternallyChecked` retains theorem/checker/receipt/scope
  identities but does not authenticate or independently execute them. RD.1c
  owns that promotion. Likewise a `StratificationClassV1::Whitney*` or `Thom`
  payload is exact metadata attached to a witness, not authority minted by this
  structural checker.
- `PresentationScopeV1` is deliberately scoped. Literal presentations receive
  no coordinate/generator equivalence; fixed-resolution and externally checked
  variants identify their asserted scope but do not supply the morphisms,
  composition laws, quasi-isomorphism receipts, refinement variance, or
  physical crosswalks owned by RD.1b.
- `AdmittedDerivedMorphismV1` proves only bounded structural compatibility,
  deterministic identity, declared evidence orientation/rank shape, and an
  explicit no-equivalence boundary. Nominal `DerivedEvidenceArtifactIdV1`
  values are not authenticated evidence receipts; callers can declare their
  bytes and ranks. RD.1c or a successor admitted-evidence type must bind and
  independently validate payload authority before any evidence-preservation
  claim is promotable.
- `AdmittedDerivedEvidencePolarityTransportCandidateV1` proves only that one
  exact sealed nonidentity morphism child retained the derived variance,
  evidence artifacts, and ranks encoded recursively in the packet and that its
  two caller-declared polarity tags are equal. It does not authenticate a
  proposition or payload, make one polarity canonical, prevent another nominal
  packet for the same child, execute the transport, or prove evidence
  preservation, functoriality, naturality, theorem truth, inverse, equivalence,
  or physical balance. It has no conversion to an authority-bearing evidence
  receipt or `DerivedEquivalenceBoundaryV1`.
- `DeclaredChartMap` proves neither overlap coverage nor any property of the
  nominal forward-map artifact. It does not supply an inverse, round-trip law,
  atlas compatibility, coordinate equivalence, or physical correspondence.
  Those require a separate scoped-equivalence receipt with independently
  checked inverse laws; `IdentityOnly` is refused for every declared chart map.
- `AdmittedDerivedChartTransitionInverseLawCandidateV1` proves only that two
  direct declared chart-map children have exact reversed selectors, one nominal
  overlap, one variance, and structurally composable evidence seams in both
  orders. It does not execute the maps, authenticate either nominal round trip,
  prove composite identity, atlas compatibility, inverse, coordinate
  equivalence, evidence preservation, or physical correspondence.
- `DeclaredComplexRefinement` proves only a same-chart, structurally monotone
  finite graded-rank/truncation envelope. It does not prove that the target is
  geometrically or numerically finer; that prolongation exists, is linear,
  injective, degree-preserving, or unit-preserving; that differentials commute;
  or that exactness, cohomology, quasi-isomorphism, remainder inclusion, error
  reduction, convergence, constraints, strata, physics, or evidence authority
  are preserved. Nonzero prolongation and commutation IDs name nominal artifacts
  only, and `IdentityOnly` is refused.
- `AdmittedDerivedStratificationRefinementCandidateV1` proves only finite
  selector coverage and dimension monotonicity for one exact exhaustive
  fine-to-coarse component family. It proves no containment, incidence/frontier
  or link compatibility, Whitney/Thom condition, component execution/gluing,
  refinement theorem, evidence preservation, or equivalence. Its ordered
  composition candidate proves only exact middle-selector alignment; it does
  not mint a direct map or transitivity.
- `AdmittedDerivedFixedResolutionQuasiIsomorphismCandidateV1` proves only that
  exact endpoint local models select the declared role complexes under matching
  endpoint-local `FixedResolution` scopes, and that an exact homogeneous
  refinement path connects those selectors. It does not validate differential
  payloads or square-zero laws; prolongation existence, linearity, degree/unit/
  quantity preservation, or commutation; induced cohomology maps or their
  injectivity/surjectivity; Betti, Euler, or torsion agreement; an inverse,
  zigzag, or chain homotopy; presentation, generator, constraint, stratum,
  geometry, evidence, or physical equivalence; refinement naturality or
  invariance; or error reduction and convergence. Its theorem, checker, check-
  receipt, scope-witness, prolongation, and commutation IDs remain nominal.
  RD.1c must independently resolve the exact candidate and path, verify every
  differential and factor map (including truncation boundaries), establish the
  induced degreewise cohomology isomorphism under the declared coefficient
  semantics, and bind checker/version/TCB/budget evidence in a new authority
  receipt. Stronger presentation equivalence additionally requires semantic
  generator/constraint correspondence plus inverse or zigzag and homotopy
  evidence; extension to later resolutions requires a separate naturality or
  coherent-zigzag theorem.
- `AdmittedDerivedLocalPresentationCorrespondenceCandidateV1` proves only exact
  endpoint/model ownership, shared structural conventions and local scope,
  canonical duplicate-free relation retention, and finite two-sided coverage of
  each selected equality, active-inequality, active-contact, and constitutive
  family. It does not compare or authenticate equation/function payloads,
  inequality senses, active-set witnesses, normal cones, contact sides/gaps/
  laws, constitutive roles/state dimensions, regularity, units, computability,
  locality maps, or the nominal edge/aggregate artifacts. Coverage does not
  establish functionality, injectivity, surjectivity, a bijection, generator
  equivalence, coordinate equivalence, physical correspondence, inverse,
  homotopy, or presentation equivalence. Repeated source and target members are
  therefore permitted; unrepresented added members refuse. RD.1c must resolve
  and semantically check every relation edge and combine any promoted result
  with the required role-specific quasi-isomorphisms and inverse/zigzag/homotopy
  evidence in a new authority receipt.
- `AdmittedDerivedScopedPresentationEquivalenceCandidateAssemblyV1` proves only
  that exact sealed child IDs were supplied for the tangent, cotangent, and
  deformation-obstruction fixed-resolution quasi-isomorphism candidate roles
  and for one exhaustive local-presentation correspondence candidate; that all
  four children retain the exact declared geometry/local-model endpoints; and
  that the three role candidates retain identical source/target resolution-ID
  and fixed-scope-witness selectors. Those selectors are derived from sealed
  children rather than accepted from the caller. The assembly does not compare
  the per-complex `FiniteResolutionV1` envelopes behind equal nominal selectors;
  equal selector IDs therefore do not prove equal rank, degree, truncation, or
  remainder envelopes. This assembly also does not
  execute or authenticate any child theorem, checker, path, relation, or
  artifact and does not establish chain-map laws, cohomology isomorphism,
  generator or constraint semantics, an inverse, a zigzag, a homotopy,
  coherence, naturality across refinements, evidence preservation, coordinate
  equivalence, or physical equivalence. It intentionally has no conversion to
  `DerivedEquivalenceBoundaryV1`, composition, inverse, or transport API. RD.1c
  must independently resolve and validate every child and then bind explicit
  inverse-or-zigzag plus homotopy/coherence evidence in a new authority-bearing
  receipt; the RD.1b assembly identity itself is only a deterministic structural
  packet and its `no_authority` artifact remains mandatory.
- `DeclaredInclusion` is a typed whole-object inclusion declaration, not a
  containment certificate. It proves no actual subset or validity-domain
  inclusion; map execution or payload truth; injectivity, monicity, embedding,
  image characterization, dimension/codimension relation, or preservation of
  charts, constraints, strata, local models, complexes, topology, physics, or
  evidence. Composition retains an ordered path of nominal factor declarations
  but does not prove transitive containment. Nonzero map and containment IDs are
  nominal only, and `IdentityOnly` is refused.
- `AdmittedDerivedStratumMorphismV1` proves only exact ownership of its endpoint
  triples and retains nominal map/constructibility artifacts for each component.
  It proves neither that a component executes nor that its image lies in the
  target stratum; it supplies no exhaustive source-stratum coverage or globally
  stratified/constructible map. Continuity, smoothness, properness, incidence or
  frontier preservation, Whitney/Thom compatibility, local-link preservation,
  constructible sheaf/cosheaf pullback or pushforward, naturality, base change,
  evidence authority, inverse, equivalence, theorem truth, and physical
  correspondence are all outside this token. Composition retains exact adjacent
  triples and nominal factor declarations only, and no API promotes the path to
  a whole-geometry morphism.
- `AdmittedDerivedExhaustiveStratifiedMapCandidateV1` proves only that every
  exact source stratum has one canonically ordered direct sealed component into
  the exact target stratification. It does not prove target coverage or
  uniqueness; map/component execution; image landing; component gluing;
  continuity, smoothness, properness, or totality; incidence/frontier,
  local-link, Whitney, or Thom compatibility; global constructibility;
  sheaf/cosheaf functoriality, naturality, base change, or evidence transport;
  inverse, equivalence, theorem truth, or physical correspondence. Its assembly
  and global-constructibility IDs remain nominal. RD.1c must resolve and execute
  every component, check the selected global policy and retained falsifiers,
  and mint a new authority receipt rather than upgrading this candidate in
  place.
- `AdmittedDerivedExhaustiveStratifiedMapCompositionCandidateV1` proves only
  that two exact sealed exhaustive-map candidates are retained in order and
  that the first target geometry/stratification exactly equals the second
  source geometry/stratification. It does not compose or execute components,
  synthesize a direct exhaustive or whole-geometry map, authenticate either
  child assembly/constructibility declaration or its own nominal composition,
  or prove total/function semantics, gluing, continuity, constructibility,
  incidence/frontier/link/Whitney/Thom preservation, categorical units or
  associativity, evidence transport, theorem truth, inverse, equivalence, or
  physical correspondence. RD.1c must independently execute and check any
  promoted composition and mint a distinct authority receipt.
- `AdmittedDerivedParallelMorphismComparisonCandidateV1` proves only that two
  exact ordered structural paths share exact endpoints and are retained under
  nominal comparison metadata. It proves no equality, commutativity, homotopy,
  naturality, coherence, execution agreement, evidence transport, inverse, or
  equivalence or physical agreement. Identity-versus-cycle and direct-versus-
  composite packets remain merely well-typed comparison inputs.
- `AdmittedDerivedSpanCorrespondenceV1` proves only that two exact sealed legs
  share the declared apex and land at the declared outer endpoints. It proves no
  totality, single-valuedness, functionality, injectivity, surjectivity,
  nonemptiness, properness, closedness, graph/image truth, physical
  correspondence, direct outer-endpoint evidence transport, inverse, or
  equivalence. No pullback, base-change, pull-push, Beck-Chevalley, projection-
  formula, or composition authority is available. Even an identity-left
  graph-shaped span remains only a structural declaration.
- `AdmittedDerivedSpanEvidenceTransportCandidateV1` proves only that one exact
  sealed span is cross-bound to the exact typed polarity packets wrapping its
  left restriction and right balance-corestriction legs, with exact artifact,
  rank, and polarity continuity at their apex seam. It does not authenticate a
  proposition or evidence payload; execute or construct a direct outer
  transport; make the correspondence total, single-valued, or functional; or
  prove functoriality, naturality, preservation, pull-push, base-change,
  Beck-Chevalley, projection-formula, theorem, inverse, equivalence, or physical
  authority. It exposes no conversion to an authority-bearing evidence receipt,
  `AdmittedDerivedMorphismV1`, or `DerivedEquivalenceBoundaryV1`.
- `AdmittedDerivedSpanMorphismCandidateV1` proves only exact fixed-foot/apex-map
  wiring and exact ordered comparison routes. Neither nominal comparison proves
  a commuting triangle or any path equality/homotopy, so the token is not a
  natural transformation, invertible 2-cell, inverse, span identity/composite,
  or equivalence. It grants no execution, evidence preservation, functionality,
  physical correspondence, associativity/coherence/interchange, pullback,
  base-change, Beck-Chevalley, or projection-formula authority.
- `AdmittedDerivedSpanPullbackSquareCandidateV1` and
  `AdmittedDerivedSpanCompositionCandidateV1` prove only exact proposed-apex,
  projection, middle-route, outer-span, and outer-route wiring. They prove no
  square commutativity, pullback existence/universality/uniqueness/nonemptiness,
  outer-correspondence functionality, categorical or bicategorical composition,
  unit/associativity/coherence, base-change, Beck-Chevalley,
  projection-formula, pull-push, evidence, equivalence, or physical-
  correspondence law.
- `AdmittedDerivedMorphismInverseLawCandidateV1` proves only exact reversed
  endpoints and exact ordered identity-versus-round-trip comparison wiring.
  Its comparison relations remain nominal; it proves no path equality,
  left/right inverse, isomorphism, homotopy or quasi-isomorphism, evidence-cycle
  identity, execution, evidence preservation, naturality, simplification,
  coherence, or physical correspondence.
- `AdmittedDerivedSpanEquivalenceCandidateV1` proves only that exact source and
  target spans, ordered forward/reverse span-morphism packets, and the exact
  apex inverse-law packet cross-bind transitively. Because none of those nested
  relations is authenticated, the assembly proves no span equivalence,
  commuting triangle, path equality, invertible 2-cell, apex inverse,
  Morita-style relation, unit/counit or triangle identity, naturality,
  coherence/interchange, composition/pullback/base-change law, execution,
  functionality, semantic agreement, evidence preservation, or physical
  correspondence. It intentionally has no conversion to
  `DerivedEquivalenceBoundaryV1`.
- The v1 sublanguage containing identities, generic strict maps, and declared
  chart maps, finite-complex refinements, and whole-object inclusion declarations
  closes mixed map-family composition through ordered typed heterogeneous paths
  whenever its exact geometry, evidence, adjacent-chart, and adjacent-refinement
  seams pass. This removes the implemented-family closure gap but does not
  promote the full RD.1b category claim. Separate structural candidates now
  retain ordered exhaustive-map compositions, finite stratification refinements,
  evidence polarity, parallel comparisons, direct chart and generic inverse-law
  shapes,
  pullback/composition-shaped spans, span equivalence, scoped local-presentation
  equivalence, and the L6 Machine-IR crosswalk. Those packets close typed seam-
  shape gaps only. Authenticated
  subobject inclusion, geometric/stratified refinement, map execution,
  categorical pullback/composition/equivalence, quasi-isomorphism, evidence
  preservation, and physical crosswalk authority still require independent
  checker receipts before the broader claim is promotable.
- V1 refuses unbounded and infinite-dimensional local models, opaque external
  analytic functions, unknown compactness/regularity, and infinite computation.
  These are admitted-class limits, not claims that the excluded mathematics is
  invalid. A later version may expand the finite theorem/checker envelope
  without reinterpreting v1 receipts.

## No-claim boundaries (exit-path theorem families)

- `ValidatedExitPathFamilyV1` establishes only bounded structural admission,
  canonical identity, and deterministic evaluation of declared sufficient-
  hypothesis predicates. It does not prove an exit-path equivalence,
  constructible descent, link contractibility, monodromy triviality, refinement
  invariance, or any theorem-card claim. V1 deliberately exposes only
  `ScientificCorrectnessNotProven`; RD.X2 owns retained finite artifacts and an
  independent checker lane owns any future authority promotion.
- Opaque witness, countermodel, theorem, checker, and no-claim identities are
  content references, not authenticated evidence. Admission checks that they
  are nonzero and internally referenced consistently, but does not dereference,
  execute, or authenticate them. The preregistered same-incidence countermodels
  therefore lock the required adversarial *classes* and artifacts, not the
  truth of their topological premises.
- `declared_wall_seconds` is a positive, identity-bound semantic declaration;
  it is not an execution receipt. Cooperative `Cx` cancellation and hard
  referenced-artifact/truncation caps bound the local work, while a deadline-
  enforcing executor remains a successor responsibility.

## No-claim boundaries (sheaf)

- Restriction maps are POINT SAMPLERS on the shared zero band;
  spline-trace and mesh-edge restriction assemblies land with their
  consuming beads (fs-iga mortar, MORPH conformance).
- Reported margins are aggregated directly from fs-ivl's outward interval
  endpoints. A non-finite/whole mismatch interval keeps an infinite upper
  report; an `Unknown` verdict or any indeterminate interface publishes
  `NumericalCertificate::NoClaim`, never an enclosure reconstructed from a
  human-facing approximation.
- The gauge-fit-share field is optional: finite least-squares diagnostics
  report fractional reduction in uncentered weighted mean-square mismatch
  energy in `[0, 1]`; if their unscaled public diagnostic arithmetic is not
  representable, the field is `None` rather than NaN or a fabricated split.
  The Laplacian optimizer is assembled from per-interface sample sums and
  counts, while the reported before/after energy is evaluated over every
  retained sample midpoint mismatch. It is a graph-gauge diagnostic, not
  evidence of cocycle membership, exactness, or an H¹ class.
- BDDC-style coarse spaces from harmonic sections (the second consumer)
  belong to the solver-dd bead; the spectral-gap confidence signal to
  Proposal 5.
- Scaling target (hundreds of patches) is structural (O(P²) overlap
  discovery + linear sampling); measured perf gates land with MORPH
  conformance.
- Chart samples with `Estimate`, `NoClaim`, or malformed rigorous error
  certificates poison their own interfaces to infinite enclosures. They can
  never contribute positive authority; an independent, determinate interface
  may still establish a localized `Fail`, while the aggregate numerical
  certificate remains `NoClaim`.

## No-claim boundaries (identity + diff)

- The diff compares CHART FIELDS (signed distance); solver-field diffs
  (stress, velocity) join when field charts land — the quantity tag is
  already plural-ready.
- Contribution measurement requires generation snapshots (one world per
  divergent op); without them, causes carry the total on the first
  touching op — explicitly unpartitioned, never silently split.
- fs-ledger `explain()` integration (walking real provenance trees
  instead of caller-supplied divergent-op lists) lands with the bisect
  bead, which owns deterministic replay.
- The R3 kill-metric wiring (quarterly fallback-fraction review) is
  governance (xpck.6); this module measures and reports the number.

## No-claim boundaries (repair)

- Solves are Gauss–Seidel over small complexes; fs-la eigensolver
  integration (spectral gap → merge confidence, Proposal 5) lands with
  the knh1.3 bead.
- Gauge repair adjusts patch potentials (constant offsets); non-constant
  boundary control-point projection (the NURBS example) lands with the
  converter beads that own those charts.
- The auto-apply POLICY (when to apply without human acceptance) is
  session governance. For the current constant-gauge model, this module first
  intersects the admissible constant-shift interval independently on each
  connected component, chooses its deterministic maximum-slack midpoint (or a
  finite edge for half-unbounded intervals), and computes
  `gauge_step_eligible` against that gauge. This field authorizes only the
  budget-feasible exact-component step, never complete repair. It refuses when no
  representative fits every per-patch budget.
- Reported harmonic support first requires the entire remainder's squared-norm
  ratio to exceed `COMPONENT_FLOOR`, then retains entries above the
  within-component relative amplitude floor
  `sqrt(COMPONENT_FLOOR) * max(|h|)`. The raw split still retains sub-floor
  residue for diagnosis. This two-stage rule is unit-rescaling stable but is
  not a graph-theoretic minimal cut-set; graph-min-cut refinement over weighted
  topologies is future work if fixtures demand it.

## No-claim boundaries (merge)

- `Trivial` proves only exact branch/base equality and carries no post-state
  residual receipt. `MergeResidualReceipt` checks nominal `f64` edge values;
  it is not the interval seam certificate, and the outside-ray sample validator
  is itself only a bounded sampling/input diagnostic rather than an independent
  falsifier.
- A `CandidateRemainderConflict` is an operational auto-merge refusal from a
  fixed-iteration heuristic. It is not an H¹ class, a graph-minimal cut, or
  evidence that no local repair exists. Diagnostic component ratios are not a
  certified orthogonal energy partition.
- Malformed skeletons (including duplicate cells), lengths, tolerances,
  weights, and non-finite cochains return `Refused` before fast-path or
  resolution authority. Non-finite weighted-Laplacian, union, decomposition,
  or post-gauge arithmetic also refuses.
- The gap proxy is the weighted VERTEX-Laplacian algebraic
  connectivity; full sheaf-Laplacian edge spectra land with the
  Proposal-5/eigensolver integration (knh1.3).
- Coupling-graph LEGALITY of merged assignments is fs-iface's contract at its
  own layer. The current merge API lacks base assignments and refuses every
  assignment payload; its pairwise difference helper cannot by itself prove a
  three-way collision.
- Merge operates on interface cochains (gauge states); merging chart
  GEOMETRY payloads routes through the converters + semantic diff.
- The seeded harness measures candidate-remainder conflicts only. A full kill
  measurement must retain candidate conflicts, escalations, refusals, and type
  conflicts from the same realistic swarm-concurrency trials; that quarterly
  trial and any fallback to ownership partitioning are governance decisions
  (xpck.6).
