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
- `Chart` (object-safe): `eval(x, &Cx)`, `support()`, `topology_hint()`,
  `name()`, `differentiability()`, provided `inside()`. Implementations
  poll `cx.checkpoint()` at bounded strides. The plan's `type Param`
  lives on the `DesignChart: Chart` subtrait so `Region` can hold
  heterogeneous `Arc<dyn Chart>` (same contract, object-safe core;
  fs-xform builds on `DesignChart`).
- `Region` — charts + per-chart `ProvenanceHash`; deterministic
  `primary()` (first); `check_agreement(&AgreementConfig, &Cx) ->
  Result<AgreementReport, Cancelled>` — seeded sampling over the inflated
  union support; two charts agree at x iff their signed distances differ
  by at most the SUM of their declared error half-widths plus the
  configured slack; failures localize (point, chart names, gap vs
  allowed, first-K cap) and reports render canonical JSON.
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
  exact `SphereChart`/`BoxChart`/`TorusChart` (unit Lipschitz, Exact error
  models, known Betti numbers) and `LyingSphereChart` (deliberately biased
  with a lying error model) for detection tests.

- `router` (the Rep Router, Bet 1): converter-edge registry
  (`ConverterSpec`: cost model, error model with declared composition rule,
  certificate availability), exact Pareto label-correcting planner over
  (cost, composed absolute error, uncertified-edge count) with the
  deterministic winner rule certified-preferred → min cost → min error →
  lexicographic path; `explain()` returns every Pareto candidate and why
  the winner won; refusals name the binding constraint (error/cost/
  no-path) with ranked relaxations; `execute()` runs chains through
  `EdgeRunner`s, composes per-edge Evidence receipts additively, and
  records actuals through `CostOracle` (an L2-clean abstraction — HELM
  wires the ledger tune table behind it; `MemoryCostOracle` in-process).
  Learned measurements replace declared error magnitudes ONLY on
  uncertified additive edges; certificates are never learned away.

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

## Invariants
1. Trait laws (G0, geo-001, 12k seeded queries): `inside(x)` ⇔ `sd(x) <
   0`; `support()` bounds the region (no negative sd outside, to
   tolerance); certified Lipschitz bounds hold along random steps;
   claimed gradients are unit-norm and match central differences on
   smooth fixtures.
2. Agreement soundness: identical presentations always agree; a
   disagreement implies at least one chart's geometry or error
   declaration is wrong — proven by detecting an undeclared 0.03 bias
   with gap-localized diagnostics naming the lying chart (geo-003).
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
`fs_exec::Cancelled` for interrupted checks, `fs_evidence::CertifyError`
through receipts. Constructors are total (`Aabb::new` normalizes); no
panics cross the boundary.

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
None yet. `[F]`-tagged bead: the abstractions are frontier by AMBITION
(category-of-charts doctrine), but everything shipped here is exercised
by the conformance suite — no flag needed until moonshot certificates
(sheaf gluing) arrive.

## Conformance tests
tests/conformance.rs, cases geo-001..geo-005 (JSON-line verdicts; seeded
cases carry seeds): the fixture trait-law battery, multi-chart agreement
within composed bounds + G5 replay, lying-chart detection with localized
diagnostics, rigorous conversion receipts + teaching refusals, and
cancellation. In-module suites cover Aabb/vector laws, fixture known
values, agreement determinism, and cancellation.

## No-claim boundaries
- NO watertightness/manifoldness/self-intersection certificates here —
  those are wqd.7 (validity certificates) and the sheaf bead; agreement
  checking is SAMPLED evidence, not a proof (its verdict is "no
  counterexample found at these seeded points + declared errors").
- `SampledSdf` claims no Lipschitz bound for its interpolant and no
  gradients (rep-sdf's job); its outside-box enclosure relies on the
  SOURCE's certified Lipschitz constant being truthful.
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
