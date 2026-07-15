# CONTRACT: fs-geom

## Purpose and layer
The Region/Chart abstraction (plan §7.1): abstract regions presented
through charts; agreement between presentations as a checkable, localized
proposition; certified conversion receipts (fs-evidence) as the Error
Ledger's geometry feed; no privileged representation, ever. Layer: L2.
Depends on fs-exec (Cx), fs-evidence, fs-ivl, fs-alloc, fs-obs, fs-sparse.

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
  GAUGE-CORRECTION PLANNING. `hodge_decompose` currently performs a deterministic,
  fixed-iteration sequence: least-squares fit to the coboundary image, then a
  least-squares fit of the residual to the retained triangle-coboundary image,
  followed by a remainder. The retained fixtures compare the first fit with an
  independent dense reference, but a generic result is not yet a certified
  orthogonal Hodge decomposition. Its INTERPRETATION CONTRACT is: fitted exact
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
  non-exactness before labeling that fixture H¹. `plan_repair` emits ranked
  agent-facing proposals. A finite post-repair seam norm is attached only to
  the constructive gauge proposal; diagnostic and reroute proposals use `+∞`
  to mean that no comparable post-state seam norm has been established.
  Optional Rep-Router reroutes retain modeled cost. `apply_gauge` subtracts one
  explicit coboundary from the retained mismatch vector; it does not edit or
  re-evaluate a chart, commit geometry, or prove realization of that correction.
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

## Error model
Structured teaching values throughout: `ConvertDiag` (ranked fixes),
`fs_exec::Cancelled` for interrupted checks, `AgreementStatus::Unknown` plus
structured `AgreementUnknownReason` for non-evaluable checks, and
`fs_evidence::CertifyError` through receipts. Router execution additionally
returns `ExecuteError` with a stable missing-runner/runner/oracle-record class;
oracle evidence refusals use `CostOracleError`. Fallible sampling and conversion
entry points return structured refusals, and `Aabb::new` is total and normalizes.
Feature-gated public Hodge/repair and merge diagnostic APIs still assert some
shape/index preconditions on caller-built skeletons; their total, budgeted
`Result` replacements remain required before promotion. The base
`SheafComplex` incidence and section APIs are now structured refusals for
malformed/oversized raw parts.

## Determinism class
Deterministic: seeded sampling, insertion-ordered charts, canonical JSON
renderings; no clocks, no addresses. Float behavior inherits
fs-math-class scalar arithmetic.

## Cancellation behavior
Chart evaluation and production sampling paths take `&Cx`.
`check_agreement` polls per sample point; conversion grids, sheaf interface
draws, triple discovery, the ray-sequence validator, and semantic-diff field
draws poll `Cx` directly at deterministic bounded strides and return typed
cancellation diagnostics without publishing partial authoritative output.
Incidence assembly, mismatch assessment/section solve, Hodge/repair, and merge
diagnostics do not yet accept `Cx`; base-complex work is statically bounded and
fallibly allocated, but these APIs are not P7-complete. Chart-local polls remain
an additional inner-kernel obligation.

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
All OFF by default per the Ambition-Tag rule (the default-path chart
abstractions remain unflagged `[S]`):
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
