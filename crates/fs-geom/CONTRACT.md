# CONTRACT: fs-geom

## Purpose and layer
The Region/Chart abstraction (plan §7.1): abstract regions presented
through charts; agreement between presentations as a checkable, localized
proposition; certified conversion receipts (fs-evidence) as the Error
Ledger's geometry feed; no privileged representation, ever. Layer: L2.
Depends on fs-exec (Cx), fs-evidence, fs-alloc, fs-obs.

## Public types and semantics
- `Point3`/`Vec3`/`Aabb` — minimal geometry-local types (fs-la owns real
  linear algebra); `Aabb` normalizes corners, offers containment, union,
  inflation.
- `BettiBounds` — per-dimension (lower, upper) topology hints;
  `unknown()` is the honest default (certificates are the wqd.7 bead's).
- `ChartSample { signed_distance, gradient: Option, lipschitz: Option,
  error: NumericalCertificate }` — plan Appendix B's value + gradient +
  certified Lipschitz + DECLARED error model relative to the abstract
  region. SDF sign convention: negative inside.
- `TraceStepClaim::{NoClaim, ExactDistance, LipschitzImplicit}` — the typed
  theorem available to a ray stepper. The default is `NoClaim`: a sample-level
  `Some(lipschitz)` alone cannot mint a no-tunneling certificate.
  `ExactDistance` additionally requires an exact singleton numerical
  certificate; `LipschitzImplicit` states that the field has the represented
  region's exact sign and zero set and that each sample's bound is valid over
  the entire closed `|f|/L` step ball, so that radius is safe but not a
  geometric-distance upper bound.
- `Chart` (object-safe): `eval(x, &Cx)`, `support()`, `trace_step_claim()`,
  `topology_hint()`, `name()`, `differentiability()`, provided `inside()`.
  Implementations
  poll `cx.checkpoint()` at bounded strides. The plan's `type Param`
  lives on the `DesignChart: Chart` subtrait so `Region` can hold
  heterogeneous `Arc<dyn Chart>` (same contract, object-safe core;
  fs-xform builds on `DesignChart`).
- `Region` — charts + per-chart `ProvenanceHash`; deterministic
  `primary()` (first); `check_agreement(&AgreementConfig, &Cx) ->
  Result<AgreementReport, Cancelled>` — seeded sampling over the inflated
  union support with an explicit `Agreed | Disagreed | Unknown` verdict.
  Two valid chart samples agree at x iff their declared signed-distance
  intervals overlap after configured slack. Zero samples, fewer than two
  charts, invalid configuration/support, non-finite outputs, malformed
  certificates, and `NoClaim` all produce `Unknown`, never vacuous agreement.
  Reports retain the weakest certificate class used, the strongest class
  supporting a counterexample, exact total counts, signed worst excess (so a
  negative agreement margin is not rounded away), and first-K localized
  disagreement/unknown diagnostics in canonical JSON.
- `Convert<Dst>: Chart` — `convert(ErrBudget, &Cx) ->
  Result<Certified<Dst>, ConvertDiag>`: the receipt's QoI is the achieved
  absolute sd-error bound (enclosed `[0, achieved]`), provenance chained.
  `ConvertDiag` refuses EARLY with ranked fixes (`BudgetInfeasible`
  computes the needed resolution vs the cap; `NoLipschitzBound` when the
  source certifies none).
- `SampledSdf` — dense trilinear sampled-SDF target with a rigorous
  in-box bound (`L·h`, conservative) and outside-box enclosures leaning on
  the source's certified Lipschitz constant; blanket
  `impl<C: Chart> Convert<SampledSdf> for C` (specialized edges arrive
  with rep-* beads). Resolution cap `SAMPLED_SDF_MAX_RESOLUTION = 96`.
- `fixtures` (PUBLIC on purpose — the shared MORPH test vocabulary):
  valid positive-radius `SphereChart`, finite strictly three-dimensional
  `BoxChart`, and valid ring `TorusChart` instances (unit
  Lipschitz, Exact error models, `ExactDistance` trace claim, known Betti
  numbers); degenerate/invalid boxes downgrade to `NoClaim` and unknown
  topology, while horn/spindle torus parameters downgrade to
  `LipschitzImplicit`, `Estimate`, and unknown topology. `LyingSphereChart` is
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

- `sheaf` module (bead wqd.13, Bet 11 [F/M]): cellular-sheaf
  WATERTIGHTNESS certificates. `SheafComplex::from_charts` discovers
  interfaces by support overlap + shared zero-band sampling
  (geometry-seeded, index-free — re-index invariance is exact), plus
  triple junctions as 2-cells. δ⁰/δ¹ assemble as fs-sparse matrices with
  entries in {−1, 0, +1}; `δ¹·δ⁰ = 0` holds BITWISE (small-integer f64).
  `watertightness(tol)` returns `Evidence<SheafVerdict>`: PASS requires
  every sample's |mismatch| enclosure INSIDE `[0, tol]` via fs-ivl's
  sound predicates (no bound extraction); FAIL requires an enclosure
  ENTIRELY above tol — the H¹ obstruction with the offending interface
  cells and magnitudes attached; anything else is an honest Unknown.
  `section_solve` computes per-patch gauge offsets over the adjacency
  Laplacian, splitting mismatch into a reconcilable coboundary share and
  the structural residual — the exact split Proposal 10's merge
  semantics reuses. `ray_parity_falsifier` is the independent
  cross-examination (registry pairing: watertightness → ray-parity).

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
  `semantic_diff` aligns worlds by `EntityId`, measures field
  differences on shared support (the sheaf band-sampling machinery),
  and attributes each finding to a RANKED list of contributing causal
  edits with per-edit contributions MEASURED across generation
  snapshots when supplied (unpartitioned-but-flagged otherwise).
  Unidentified entities degrade to a geometric fallback FLAGGED
  `attributed: false`, and the fallback fraction is the R3
  early-warning metric. `DiffReport::filter` triages by
  region/quantity/magnitude.

- `sheaf_repair` module (patch Rev L, bead wqd.14; [M], behind the
  `sheaf-repair` feature until certifier trials pass): DIAGNOSIS →
  CONSTRUCTIVE REPAIR. `hodge_decompose` splits the interface mismatch
  cochain into exact ⊕ coexact ⊕ harmonic over the skeleton
  (Gauss–Seidel normal equations; dense-oracle-verified), with the
  INTERPRETATION CONTRACT: exact → local gauge repair (auto-appliable
  ONLY when every per-patch offset fits that chart's declared error
  budget — a repair never silently distorts geometry); coexact →
  circulation around triple junctions, diagnosed CONVERTER-side (not a
  geometry edit); harmonic → true topological obstruction, declared
  unrepairable-locally with the interface cut-set. `plan_repair` emits
  ranked agent-facing proposals (expected post-repair norms; optional
  Rep-Router reroute with modeled cost); `apply_gauge` is the
  constructive, idempotent step.

- `sheaf_merge` module (Proposal 10's CROWN JEWEL, bead lmp4.12; [M],
  behind the `sheaf-merge` feature — which enables `sheaf-repair` —
  until its Gauntlet tier + kill metric are green): the sheaf machinery
  as a merge-conflict classifier. `three_way_merge` forms the union of
  edits (X + Y − B at the cochain level), Hodge-decomposes, applies the
  canonical least-squares gauge reconciliation, and RE-VERIFIES the
  reconciled state's own certificate before reporting resolved (Sev-0:
  a passing certificate is never attached over a failing state).
  Verification failures classify: a dominant harmonic residue above
  tolerance is a STRUCTURAL CONFLICT localized to its supporting cells
  with both parents' provenance; anything else (coexact circulation)
  ESCALATES unresolved. Auto-resolution is licensed exactly when the
  reconciled state passes — a harmonic remnant below the watertightness
  tolerance is not an obstruction (a lesson the kill-criterion harness
  taught: machine-floor triggers made every noisy merge conflict).
  Type-level collisions (same key, different values) are caught BEFORE
  decomposition; trust is conditioned on `spectral_gap` (weighted
  algebraic connectivity, Jacobi eigenvalues) with LowGap flagging
  (R5). `harmonic_conflict_rate` is the kill-criterion measurement
  (25% line).

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

## Error model
Structured teaching values throughout: `ConvertDiag` (ranked fixes),
`fs_exec::Cancelled` for interrupted checks, `AgreementStatus::Unknown` plus
structured `AgreementUnknownReason` for non-evaluable checks, and
`fs_evidence::CertifyError` through receipts. Router execution additionally
returns `ExecuteError` with a stable missing-runner/runner/oracle-record class;
oracle evidence refusals use `CostOracleError`. Constructors are total
(`Aabb::new` normalizes); no panics cross the boundary.

## Determinism class
Deterministic: seeded sampling, insertion-ordered charts, canonical JSON
renderings; no clocks, no addresses. Float behavior inherits
fs-math-class scalar arithmetic.

## Cancellation behavior
Every query path takes `&Cx`. `check_agreement` polls per sample point;
`convert` grids poll through the source chart's `eval`. Long geometry is
interruptible like any kernel (P7); fixtures are O(1) per query.

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

## No-claim boundaries
- NO watertightness/manifoldness/self-intersection certificates here —
  those are wqd.7 (validity certificates) and the sheaf bead; agreement
  checking is SAMPLED evidence, not a proof (`Agreed` means "no
  counterexample found at these seeded points under the reported certificate
  strength + declared intervals + configured slack").
- Two registered presentations are sufficient to run the pairwise check; this
  module does not certify implementation or provenance independence between
  them. Independence must be established by the campaign that cites the result.
- `SampledSdf` claims no Lipschitz bound for its interpolant and no
  gradients (rep-sdf's job); its outside-box enclosure relies on the
  SOURCE's certified Lipschitz constant being truthful.
- `TraceStepClaim::LipschitzImplicit` certifies no-tunneling step radii, not
  Euclidean proximity from a small normalized residual. Consumers must retain
  that distinction in hit/error language.
- `ConverterSpec::certified` is a declaration, not an authenticated admission
  receipt. Runtime certified runners must return `Certified<f64>` local-error
  evidence and fs-ir rejects routes containing an estimated declaration, but a
  malicious caller can still lie in the Boolean until the opaque admitted-color
  / converter-authority migration lands.
- The trilinear bound uses the conservative `L·h` constant; sharper
  (√3/2) constants and fs-ivl interval-verified sampling arrive with
  rep-sdf.
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
- Reported margins are midpoint±width reconstructions (≤1 ulp); the
  VERDICT itself is decided only by fs-ivl's sound interval predicates.
- BDDC-style coarse spaces from harmonic sections (the second consumer)
  belong to the solver-dd bead; the spectral-gap confidence signal to
  Proposal 5.
- Scaling target (hundreds of patches) is structural (O(P²) overlap
  discovery + linear sampling); measured perf gates land with MORPH
  conformance.
- Charts with NoClaim error certificates poison their interfaces to
  infinite enclosures — such models can only ever be Unknown (honest).

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
  session governance; this module computes `auto_repairable` and
  refuses over-budget repairs in the proposal text.
- The harmonic cut-set is the harmonic component's support, minimal for
  the reported cochain; graph-min-cut refinement over weighted
  topologies is future work if fixtures demand it.

## No-claim boundaries (merge)

- The gap proxy is the weighted VERTEX-Laplacian algebraic
  connectivity; full sheaf-Laplacian edge spectra land with the
  Proposal-5/eigensolver integration (knh1.3).
- Coupling-graph LEGALITY of merged assignments is fs-iface's contract
  at its own layer; this module catches keyed collisions.
- Merge operates on interface cochains (gauge states); merging chart
  GEOMETRY payloads routes through the converters + semantic diff.
- The kill measurement here is the harness; the quarterly swarm
  concurrency TRIAL and any fallback to ownership partitioning are
  governance decisions (xpck.6).
