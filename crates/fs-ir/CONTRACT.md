# CONTRACT: fs-ir

> Status: ACTIVE (FrankenScript core, IR language v3). Owns the typed AST,
> both concrete syntaxes, study recognition, and verb lowering. Admission
> (dimensional/chart/budget/capability checks) is the gp3.5 bead's;
> the operator catalog is gp3.6's.

## Purpose and layer

FrankenScript — the system's one true interface (plan §11.1, Decalogue
P10): a typed, versioned IR with two isomorphic concrete syntaxes
(canonical s-expressions; lossless JSON mapping), both parsing to the same
typed AST. Layer: L6 (HELM). Production dependencies are declared in
`Cargo.toml`; fs-blake3 supplies the default [S] identity kernel, while the
default fs-exec/fs-scenario seam supplies cancellation-correct PR-5 scenario
projection, and the moonshot `derived-crosswalk` feature enables fs-geom's
admitted derived-geometry boundary.

## Public types and semantics

- `Node`/`NodeKind`/`Span` — every node carries a byte span. Atoms are the
  real nouns: `Int`, `Float` (finite only), `Qty` (fs-qty SI value + dims +
  validated source literal text — fs-qty normalizes 65deg → rad, so the
  source spelling is retained in the parsed AST while canonical identities
  use one checked SI-base encoder; equality is value+dims), `Count`
  (information/core grants: B/KiB/MiB/GiB/cores —
  deliberately outside fs-qty's SI domain), `Seed` (0x… u64), `Str`,
  `Symbol`, `Keyword`, `List`.
- `CountValue` preserves bare integers as exact `u128` and decimal/exponent
  spellings as a bounded exact decimal (`u128` significand + base-10
  exponent). Whole-byte/core enforcement uses checked integer arithmetic;
  binary floating point is reporting-only. Mixed syntax classes remain
  distinct identities (`2B` differs from `2.0B`), while each class has one
  deterministic canonical spelling.
- `Node::same_shape` — semantic equality ignoring spans and Qty
  presentation; the isomorphism property is stated in terms of it.
- `VersionedProgram` — the persisted/replayed artifact boundary shared by both
  syntaxes. It canonically wraps a program as
  `(frankensim-ir :version 3 :program <node>)`, binds the language version into
  serialized identity, requires its persisted s-expression or JSON input to be
  in the one canonical byte encoding, and refuses older/newer semantics unless
  a caller first performs an explicit audited migration. Bare parsers remain
  syntax-only and may be used as the explicit normalization/migration entry.
- `sexpr::parse/print` — total reader with spans, comments (`;`),
  string escapes, deterministic atom classification (numeric-leading
  tokens MUST fully parse — a number with a garbage suffix is a structured
  error, never silently a symbol), depth cap (adversarial nesting refuses
  structurally). Printer output reparses to the same shape.
- `json::parse/print` — the lossless mapping (single-key tagged objects:
  i/f/q/c/seed/s/sym/kw; arrays = lists). Qty/Count/Seed reuse the s-expr
  literal grammar inside strings so ONE classifier owns numeric semantics
  for both syntaxes. Numbers follow RFC 8259 exactly, JSON strings reject raw
  controls and unpaired UTF-16 surrogates while decoding valid surrogate pairs,
  and unknown tags/tag-literal mismatches refuse with spans.
- `Study::from_node` — recognizes Appendix C study forms: name, seed,
  versions/budget/capability clauses, `(let …)` bindings, body;
  `constellation_lock()` extracts the versions pin. Duplicate Five-Explicit
  pillars and duplicate let names refuse as ambiguous instead of replacing an
  earlier declaration. Extraction only — validity POLICY lives in `admission`
  (below).
- `lower::lower` — high-level verbs (`optimize-shape`, `simulate-pour`)
  expand to explicit IR with an inspectable trace naming every injected
  default (progressive disclosure with nothing hidden); idempotent;
  malformed verb usage refuses with the verb's span. The public boundary first
  validates the complete caller-provided AST, including its depth and exact
  invalid-atom path, before recursive descent.
- `IrError` — span + stable kind code + detail + fix hint (refusals
  teach). `IR_VERSION` — the language version this build reads/writes.
- `derived_crosswalk` (`derived-crosswalk`, [M]) —
  `admit_derived_machine_model_crosswalk_candidate_v1` binds one exact sealed
  `AdmittedDerivedGeometryV1` to nominal Machine-IR model, subject, immutable
  version, frame, and unit selectors. Its domain-separated evidence-node
  receipt retains the current `IR_VERSION`, exact typed derived-geometry child plus
  redundant subject/version/frame/unit selectors, four endpoint-specific
  nominal mapping artifacts, one aggregate declaration, and a mandatory
  no-authority artifact. Admission refuses zero identities, stale schema or IR
  versions, a raw geometry ID that does not name the supplied sealed object,
  and every redundant selector mismatch. The token is structural lineage only:
  it does not bind or inspect an admitted Machine-IR graph.
- `machine` (Machine-IR E0 PR-1 through PR-5, [S]) — six nominally distinct durable entity
  types (`BodyId`, `SurfacePatchId`, `ContactFeatureId`, `TerminalId`,
  `PortId`, `StateSlotId`) use `fs-blake3::identity::EntityId` under six
  different static schemas. Their bounded, human-auditable hierarchical keys
  use the exact `[a-z][a-z0-9-]*` segment grammar; numeric-only path segments,
  implicit normalization, array indices, uppercase aliases, and unbounded text
  refuse before an ID is published. Each value retains its complete canonical
  preimage receipt for collision adjudication.
- `machine::LineageRecord::admit` canonicalizes bounded split, merge, remesh,
  wear, and fracture source-target relations by the durable IDs rather than
  caller order. Downstream cache/contact/winding/adjoint attachments are
  rebound only across a unique one-to-one source relation. A live attachment
  on a one-to-many source returns `LineageRefusal::Ambiguous` carrying a
  domain-separated `LineageInvalidation`: the complete refused relation set,
  every considered dependent, its ambiguous relation subset, every invalidated
  dependent, and its canonical identity receipt. No success record or guessed
  target is published on that path.
- `machine::LineageRecord::admit_with_decision` retains every attempt as a
  bounded `LineageAdmissionDecision`: attempted event and input counts, a stable
  decision/rule code, and either the admitted receipt or complete typed refusal.
  The core prints nothing. This is an outcome summary for structured tracing,
  not a digest of early-refused inputs; replay ledgers must retain those inputs
  separately.
- `machine::MachineGraphDraft::admit` is the dependency-neutral PR-2 graph
  boundary. Four additional graph-local durable IDs name subsystems,
  relations, clocks, and interfaces without widening the closed six-role
  topology/lineage enum. The draft contains bounded, canonically ordered
  subsystem/model declarations, typed terminals and effort/flow ports,
  directed algebraic or stateful relations, logical clocks, material-card
  bindings, and role-oriented interface bindings. External model, material,
  interface, and solve-policy references are nominal, versioned, nonzero
  semantic digests under canonical namespaces; `fs-ir` neither imports nor
  impersonates their runtime owners.
- `machine::codec` is the additive FrankenScript-to-Machine graph/behavior-v1
  boundary. One current-version `VersionedProgram` body decodes through either
  closed positional grammar into an authority-free draft; only the existing
  graph or behavior admission may publish semantic identity. Behavior syntax
  carries the exact `MachineGraphIdV1` it extends, and admission refuses an
  attempted rebind before inspecting behavior semantics. Reverse writers
  accept only admitted values, emit their canonical row order, and re-enter the
  global version envelope. Both concrete syntaxes therefore share one AST
  grammar rather than independent Machine parsers. Role-specific IDs are
  reconstructed from canonical keys; external references require a canonical
  namespace, nonzero full-range u64 version encoded as a decimal string, and
  exact lowercase 64-hex digest. Syntax/resource refusals retain a stable code,
  source span, structural path, detail, and repair hint without flattening
  semantic `MachineGraphFinding`s or `MachineBehaviorFinding`s. Graph decode
  preflights at most 262,144 AST nodes, 16 MiB aggregate atom text, all public
  graph collection caps, and the 16,384-element aggregate ownership cap.
  Behavior decode separately preflights a derived 2,359,360-node canonical
  envelope, 192 MiB atom text, all six public collection caps, and 32,768
  aggregate history-event, guard-read, reset-write, and dependence-member
  references. The admitted-value writer runs those same preflights. All
  resource passes precede recursive AST validation and decoded-vector
  allocation.
- Graph admission accumulates sorted, duplicate-free `MachineGraphFinding`s
  before publishing any `MachineGraphIdV1`. It rejects dangling or duplicate
  ownership, unsupported semantic scalar forms, quantity/shape/frame/
  orientation/clock/causality gaps, effort-flow pairs without a shared shape or
  whose checked six-base dimensions are not power, missing or multiply
  produced inputs, unowned or
  multiply/unwritten state, instantaneous cycles without an explicit
  solve-policy boundary, incomplete body-material closure, and invalid or
  conflicting interfaces. `MachineGraphAdmissionDecision` retains submitted
  collection counts plus the admitted graph or complete refusal and prints
  nothing.
- `machine::semantics::MachineBehaviorDraft::admit_against` is the PR-3
  behavior overlay and binds one exact admitted `MachineGraphIdV1` without
  modifying the already-published graph-v1 identity grammar. It gives every
  owned state slot an exact quantity/shape/clock/frame contract and exactly one
  initial source, gives every `ExternalInput` terminal exactly one matching
  boundary source, and gives every owned body an explicit static or prescribed
  motion. Conditions may name fixed values, distributions, or histories whose
  discontinuities are tied to durable event IDs; lower-layer PDE boundary kinds
  are deliberately not part of this graph language.
- Events use durable IDs, event-driven clocks, explicit state/terminal guard
  dependencies, oriented crossing semantics, and deterministic, set-valued,
  terminal, or honest-unknown reset relations. Reset writes are closed over
  owned states. Potential simultaneity is never resolved by collection order:
  every event on one clock participates in one coherent policy consisting of
  either unique superdense microsteps or one externally identified set-valued
  group with at least two members.
- Tolerances bind durable graph targets and nominal parameter artifacts to exact
  quantity/shape declarations. Bounded envelopes retain canonical asymmetric
  finite widths; random tolerances retain a positive scale plus exact law and
  marginal references. V1 admits at most one tolerance ID/law for each exact
  target-parameter pair; composition requires a future explicit relation.
  Every distribution-valued condition and random
  tolerance appears exactly once in one global mutual-independence or joint
  correlation declaration, while fixed/history conditions and bounded-only
  tolerances appear in none. The canonical member list is the correlation
  artifact's declared axis order. Correlation models are nominal, versioned
  external references: behavior admission proves membership closure, not
  covariance validity, PSD, or manufacturing-population representativeness.
- Behavior admission checks public collection limits before deeper traversal,
  canonicalizes every outer collection and nested event/reset/dependence set,
  accumulates sorted duplicate-free `MachineBehaviorFinding`s, and publishes no
  `MachineBehaviorIdV1` on refusal. `MachineBehaviorAdmissionDecision` retains
  submitted counts plus the outcome for structured tracing; it is not a replay
  digest of refused input.
- `machine::assurance::MachineAssuranceDraft::admit_against` is the PR-4
  operational overlay. It binds the exact graph and behavior identities plus
  every selected `fs-evidence::vv::SchemaAdmissionReceipt` (schema/ruleset,
  case hash, receipt hash, Context-of-Use, validation plan, experiment hashes,
  and exact QoI set). Acceptance criteria remain owned by the admitted V&V
  case; Machine IR cannot introduce a second, contradictory acceptance row.
- Sensors bind a durable sensor ID and owner to an exact terminal/state
  contract. Direct observations share the target clock; modeled resampling
  names an explicit clock-bridge artifact. A plant sensor exposes a matching
  owner/output-terminal contract, while an experiment-only sensor says so
  explicitly. Each experiment maps machine sensors one-to-one to the exact
  calibrated instrument IDs inside its admitted `ExperimentArtifact`; all
  case experiments and Context-qualified QoIs are closed.
- Hazards use durable `HazardId`s, whole-machine/subsystem/element/relation/
  interface scopes, exact Context-of-Use, requirement, operating-envelope and
  safety-case references, and nonempty links to the context case's admitted
  assumptions-ledger rows. Thus predicate/evidence/monitor/violation/owner/
  review semantics have one authority. Modeled hazards require a covering
  fault; honestly unmodeled hazards carry an explicit no-claim reference and
  cannot simultaneously receive a fault edge.
- Accounting windows name an exact context, clock, interval, signed boundary,
  extensive balance kind, quantity, contribution targets, role/orientation,
  audit policy, and unique loss owner for each dissipative contribution. Known
  energy/enthalpy/momentum/mass/charge/amount/entropy/exergy balances check
  six-base dimensions; species, element and custom balances bind an external
  law. A target has one role per window, preventing structural double count.
- Fidelity policy contains exactly one baseline per graph subsystem, and that
  baseline uses the graph's exact `ModelRef`. Every bounded rung names its
  validity domain, cost/error model, graph-model crosswalk, Context-qualified
  QoIs and falsifiers. Every rung is reachable from its subsystem baseline and
  has exactly one same-subsystem outgoing trigger: either an explicit
  state-transfer/model-crosswalk edge or refusal. The resulting graph is
  acyclic and terminating; no collection order or implicit fallback chooses a
  model. The fixed-fidelity replay oracle remains identity-semantic.
- Assurance admission checks public and nested limits before deeper traversal,
  canonicalizes every outer and nested collection, returns sorted duplicate-
  free `MachineAssuranceFinding`s, and publishes no `MachineAssuranceIdV1` on
  any graph, evidence, scope, accounting, or fidelity refusal.
- `machine::lowering::admit_machine_domain_lowering_v1` is the PR-5 one-way
  projection boundary. It accepts the exact admitted graph/behavior/assurance
  chain, one caller-materialized concrete `fs_scenario::Scenario`, exact
  versioned external domain-artifact references, an exact versioned projection
  policy, and one mapping-law-bound crosswalk row for every required Machine
  source. Required sources comprise the three aggregate identities plus every
  separately durable clock, subsystem, Machine element, relation, interface,
  event, tolerance, sensor, experiment, hazard, fault, accounting window, and
  fidelity rung. Aggregate graph/behavior/assurance identities transitively
  bind declarations that deliberately have no separate local ID; collection
  positions never become identity.
- Projection admission verifies the source chain before work, preflights
  public and aggregate resource bounds (including a hard fixed-record bound
  before the legacy scenario-plan scan), canonicalizes artifact/crosswalk
  order, rejects missing/duplicate/foreign sources and missing/ambiguous
  scenario locators, rejects unknown or orphan external artifacts, and runs
  `Scenario::validate_with_budget` under the supplied `Cx`. The
  `admit_machine_domain_lowering_with_decision_v1` variant retains bounded
  submitted artifact/crosswalk/target counts plus stable admitted/refused and
  detailed-refusal codes for structured tracing. Admission then requires a
  current-version, migration-free, value-equal, byte-identical
  `write_ir -> parse_ir -> write_ir` round trip before publishing
  `MachineDomainLoweringIdV1`.
- The projection identity binds PR-5 and FrankenScript IR versions; exact
  graph, behavior, and assurance IDs; scenario IR version and canonical byte
  hash; the explicit validation budget and preflight plan; projection policy;
  canonical external-artifact rows; every crosswalk source, target, and law;
  and the canonical manifest hash. The admitted value retains that manifest
  plus a framed portable payload containing the manifest and scenario bytes.
  `verify_replay` independently reparses and re-derives scenario, manifest,
  payload, typed identity, and canonical-preimage receipt. The standalone
  `parse_machine_domain_portable_payload_v1` boundary independently checks
  payload/manifest framing and bounds, schema/build versions, typed aggregate
  IDs, complete opaque row framing, canonical scenario parse/reprint, and
  scenario hash binding without reconstructing authority-bearing domain
  artifacts or claiming that row bodies are canonical.
- `query` (addendum Proposal 8 — declarative query language v0): a query is
  `(QoI, Target, budget_usd, deadline_s)` where `Qoi` is a fixed MENU —
  `MaxOverRegion`, `Integral` (linear), `Exceedance` (probabilistic, needs a
  named environment) — each advertising `QoiMeta { linear, adjoint_available,
  ladder_applicable }` for the planner and a `value_dims(field_dims)`
  (max→field dims; integral→field·m³; exceedance→dimensionless). `Target` is
  `Tolerance{value,dims}`, `Confidence(f64)`, or `ToleranceAndConfidence`.
  `Query::admit(&FieldRegistry) -> QueryAdmission` type-checks a query in
  constant time (no solves) over six fixed-order checks — `query.field`
  (the QoI's field must exist), `query.budget` (finite positive $),
  `query.deadline` (finite positive s), `query.confidence` (strictly in
  `(0,1)` — 100% is uncertifiable), `query.target` (finite positive
  tolerance), `query.dimensions` (tolerance dims == QoI value dims; exceedance
  threshold dims == field dims) — emitting the admission bead's `Finding`s
  with ranked teaching `RankedFix`es. `Query::from_node`/`to_node` give the
  `(query …)` IR surface, round-tripping under `same_shape`. This is the
  addendum's declarative surface; imperative solver access is the
  internal/expert path.

- `admission` (the gp3.5 bead): `admit(node, &AdmissionContext) ->
  AdmissionReport` first binds raw syntax to the current `VersionedProgram`,
  fully lowers shorthand, and binds the lowered versioned identity before any
  authority check. The receipt retains exact canonical raw and lowered
  envelopes. A caller-forged raw AST that cannot be canonicalized has no raw
  identity; a lowering/output-envelope refusal retains the raw identity but has
  no lowered identity. Both paths return one deterministic structured refusal
  with no authority timings, never a partial receipt or panic. Admission then
  runs six timed dimensions over only the lowered semantics — Five Explicits structure,
  dimensional analysis (fs-qty dims inferred bottom-up through `+ - * /
  min max` and comparisons; unknown verbs never false-reject), budget
  feasibility (fs-plan SEALED cost models, bead 2pmb: receipt-backed
  authority mintable only by the exact roofline loader, and costing any
  verb through a `ProvisionalUnaudited` model adds a once-per-verb
  `Warn` finding naming the evidence class and scope, so an admitted
  study's record states what its wall-cost evidence was; bead l2k92:
  when the caller supplies a `CostFreshnessContext` (now, machine
  fingerprint, `FreshnessPolicy`), each exact model is re-assessed via
  `SealedCostModel::assess` at the decision point, and an aged-out,
  machine-drifted, or future-recorded receipt admits with a
  once-per-verb `CostModelStale` warning naming the exact
  `StalenessVerdict` cause — downgraded evidence never silently steers
  admission as if sealed-fresh; without the context, behavior is
  unchanged (`NotApplicable`); models range
  over exact numeric `:dof`/`:size`/`:modes`
  features, with malformed or duplicate explicit features refused rather than
  priced as unit size; p90
  totals vs the `(budget (wall …))` bound, with RANKED cost-model-derived
  fixes: coarsen / surrogate-screen / relax), capability sufficiency
  (finite non-negative session grants, session-token and self-contained
  explicit globs vs namespaced verbs, and finite declared asks). Capability
  fields are exact keyword/value pairs; operator grants are exact names or
  namespace wildcards of the form `foo.*`. Wall/memory budget clauses have
  exact arity; structured operator-specific budget clauses remain extensible
  until the catalog lands. Chart
  routability (fs-geom Router as an admission predicate with the
  `RouteRefusal`'s own fixes attached; malformed/spec-mismatched oracle
  authority and bounded-search exhaustion reject distinctly; a route containing
  any converter not declared certificate-backed is estimated and cannot
  authorize admission), and regime gating (explicit
  `(assert (regime.allows …))` plus `flux.*` verbs checked against an
  fs-regime report; policy-graded Reject/Warn). Findings carry spans,
  diagnoses, and `RankedFix { action, predicted_wall_s, qoi_impact }`.

- `planner` module (addendum Proposal 8, bead lmp4.16; [F], behind
  `ladder-planner` → optional fs-verify dep): the GREEDY LADDER WALK —
  not a general planner (Governance Rule 1). The operator menu
  `{cache, speculate, solve-rung, refine-to-target, climb}` runs
  greedily with costs LEARNED from telemetry (`CostTable`; cold
  entries fall back to the conservative default). Speculation
  verifies a prolongated coarse answer WITHOUT solving; refinement is
  the textbook equidistribution criterion (split every element above
  the per-element target `tol²/n` with per-element depth from its own
  gap). A discharged answer's bound is a REAL equilibrated enclosure
  (VERIFIED color) and never violates the query's certified-accuracy
  contract. `ProblemFamily`, `CachedAnswer`, and `CostTable` have checked
  constructors and private state; `plan -> Result<PlanOutcome, PlanError>`
  validates finite theta/tolerance/budget, a non-empty strictly increasing
  non-zero rung ladder, family coefficients/boundaries, meshes, candidates,
  telemetry, verifier enclosures/refusals, and every arithmetic result. Bounded
  FEM refusals retain their structured `Fem1dError` source in `PlanError`
  instead of being flattened into a plausible result. Cache-declared
  bounds are never trusted: hits are independently re-verified against a
  lower-layer canonical MMS class identity before discharge. The family owns
  an admitted immutable `fs-verify::MmsClass`; theta scaling constructs a new
  admitted class, and cache keys bind that class's versioned canonical bytes
  instead of independently serializing planner coefficients. The retained key
  grammar is explicitly `PLANNER_CACHE_KEY_VERSION = 3` with the already-shipped
  exact `fs-ir-ladder:v3:` domain/prefix; this declaration does not re-key
  existing entries. The lowercase-hex payload is the exact canonical byte
  stream of the theta-scaled lower-layer class, so one-ULP theta changes that
  survive lower-layer admission remain distinct even when display formatting is
  identical. Signed-zero theta and signed/trailing-zero polynomial spelling are
  intentional admission normalizations. Retained adapters must call
  `admit_planner_cache_key`: changed domains, stale versions, version aliases,
  uppercase hex, and malformed payloads fail closed rather than being guessed.
  Tolerance, budget, ladder, telemetry, and observer state are deliberately not
  key fields: tolerance is enforced at lookup and every hit is independently
  re-verified, while the others affect execution policy rather than answer
  identity. Cache records
  normalize signed-zero mesh and nodal values before retention. Ladder transfers
  use deterministic P1 interpolation over the actual coarse coordinates,
  including non-dyadic ladders and adaptive-to-uniform moves. Actual solve/speculation
  costs are admitted before execution, so spend never exceeds budget. Learned
  costs rank the zero-cost refine/climb transition only; they never veto exact
  affordable work. Transition telemetry is pending until its first downstream
  verification/solve executes, records that actual compute cost before the
  resulting certificate checkpoint, and is dropped when admission aborts before
  downstream work. A rejected climb speculation completes that transition;
  any subsequent fine-rung solve is separately charged as `SolveRung`. The
  family degree is capped at five (six coefficients), matching the exactness
  envelope of the verifier's five-point squared-residual quadrature; its work
  preflight counts all five quadrature evaluations per coefficient/cell. The
  homogeneous trace, six-coefficient cap, signed/trailing-zero
  normalization, derived forcing, and stable identity all come from the single
  lower-layer fs-verify admission type. Parameter scaling is refused if
  independently rounded products no longer sum exactly to zero at `x=1`. The
  cannot-discharge boundary refuses with the best achieved certified interval;
  if no solve is affordable it returns `RefusedWithoutAnswer` with no interval
  or color rather than fabricating evidence. Operator choice tie-breaks
  deterministically (G5 replay).
  Every finite audit-log bound retains a private-constructor
  `VerifierCertificate`: the guarded `fs-evidence` color, stable verifier-family
  identity, and reconstructed-flux hash travel together instead of reminting a
  bare `Color::Verified` from a discarded scalar. The cell cap is derived as
  one less than the lower-layer `fs-verify` node cap. Resource admission caps cells,
  family coefficients, fidelity-rung count, and their coefficient-by-cell-by-
  quadrature-point work product. Uniform/adaptive mesh, indicator,
  prolongation, trajectory, and
  family-scaling allocations use fallible reservation; exact downstream solve
  cost and the combined compute envelope are admitted before mesh allocation.
  `baseline_uniform` is the fixed control the kill criterion measures
  against and is fallible under the same numerical/family boundary.

- `anytime` module (addendum Proposal 8, bead lmp4.17; ships behind
  `ladder-planner` but its CONTRACT survives even a frozen planner —
  the product win): `run_anytime_observed` drives one cumulative planner
  execution through operational certification/pre-work checkpoints. It emits
  each affordable budget-rung certificate before work, allocation, cache
  insertion, or telemetry for a later rung. An observer can return `Stop`; the
  returned report is then the exact deterministic prefix and no later side
  effect occurs. `run_anytime` is the collector wrapper whose observer always
  continues. Total work never exceeds the final budget, no certificate appears
  before its cumulative cost is affordable, and retaining the best checked
  certificate makes tightening MONOTONE. The first affordable rung is the
  IMMEDIATE wide certified interval, and every step carries its guarded
  Proposal-3 color plus verifier family/flux identity and a PRICED "what would
  tighten this" hint
  (`tighten_hint`: gap extrapolation naming the next menu move and the
  hot region where refinement concentrated; cold telemetry degrades to
  the generic priced form). REFUSAL semantics: an undischargeable
  query returns the achieved certified interval, the price of the gap,
  and the explicit no-point-estimate clause — never a silent
  best-effort number. `run_anytime -> Result<AnytimeReport, PlanError>` rejects
  empty, non-finite, non-positive, or non-increasing budget ladders before work.
  A valid rung too small for the initial solve contributes no trajectory point;
  if it is the final rung the report explicitly says that no certified interval
  or color exists. `tighten_hint` is likewise fallible and cannot emit NaN/∞
  gap prices. Hints use the cost table at emission time and cannot consume
  telemetry from future work. Budget-ladder length is resource-capped and its
  trajectory allocation is fallible. Replays reproduce trajectories,
  certificate identities, and observer-selected prefixes bit-for-bit (G5).

## Invariants

- Machine entity roles are non-confusable at both the Rust type and canonical
  identity-schema levels: identical key bytes under `BodyId` and
  `SurfacePatchId` produce different strong identities. Closed
  `MachineElementId` erasure retains the role tag.
- Lineage relation targets must preserve their source's nominal entity role.
  Sources and targets are duplicate-free. Split/fracture, merge, and wear have
  explicit cardinality laws; remesh permits declared one-to-many topology but
  cannot silently transport a live attachment across it.
- Lineage and invalidation identities bind the event kind, complete canonical
  relation set, dependent class/key/source, and every admitted target. Caller
  relation/dependent order is normalized once before identity publication. A
  refusal identity retains every considered dependent as well as the exact
  invalidated subset, so changing an otherwise unambiguous attachment cannot
  alias the refused attempt.
- Machine graph identity binds the graph-schema version and every canonical
  clock, subsystem/model, terminal semantic type and shape, port energy role,
  relation mode/policy/state, material-card reference, and oriented interface
  endpoint. Every caller collection is sorted by its durable ID; duplicates
  refuse rather than last-write-wins. A `Dimensional(Dims)` terminal is an
  explicit no-kind claim and never aliases a `Semantic(SemanticType)` terminal
  with the same six exponents.
- Internal input terminals have exactly one directed producer; explicitly
  external inputs have none. Each declared state slot has exactly one stateful
  writer owned by the target subsystem. Algebraic relations with no policy
  form the structural feedthrough graph; a cycle in that graph refuses, while
  a stateful edge or named solve-policy boundary cuts it without claiming
  numerical convergence.
- A role-oriented interface resolves two distinct ports, requires exactly one
  producer and one consumer independently for each effort/flow pair plus
  complementary energy direction, and names both relations in their actual
  directed causality. Its aligned/opposed declaration
  governs relation-orientation compatibility; it does not synthesize an
  implicit relation or silently close an input.

1. Isomorphism: `parse(print(x))` has the same shape as `x`, per syntax
   and across syntaxes (property-tested on generated programs and the
   Appendix C fixtures).
2. Both parsers are total: any input yields a value or a structured error
   with an in-bounds span; recursion is depth-capped (no stack overflow).
3. No silent reinterpretation: numeric-leading tokens either fully parse
   as int/float/quantity/count or refuse; non-finite literals refuse.
4. Count authority is exact end to end: decimal text cannot round into a
   different byte/core claim, checked unit scaling precedes admission, and
   `SessionCapability` carries integer memory/core grants without an `f64`
   projection.
5. Lowering is explicit, inspectable, and idempotent; the trace names
   every injected default. Admission binds the raw and fully lowered canonical
   versioned identities before authority checks. Invalid raw or expanded ASTs
   refuse before those checks and cannot mint a missing identity.
6. A derived/Machine-IR crosswalk candidate cannot detach its nominal mapping
   artifacts from the sealed derived geometry: geometry identity, subject,
   immutable model version, frame, and unit system must all match the supplied
   admitted object exactly, and every retained field is identity-semantic.

- Admission determinism: same study + context → byte-identical
  `diagnosis()`; findings sorted (check, span).
- Admission latency is milliseconds-class on Appendix C studies (six
  checks timed individually; conformance logs and bounds the total).
- Zero false admits on the violation zoo; missing verifiers (no Router,
  no RegimeReport) degrade to WARN verification-gap findings, never to
  silent admits of violations they could not check.

## Error model

`MachineIdError` names the applicable role, segment, byte offset, and
bounded-key rule that refused an entity/dependent/graph key. `MachineReferenceError`
refuses noncanonical external namespaces and zero semantic digests.
`LineageRefusal` distinguishes
shape/resource limits, cross-role or duplicate endpoints, missing attachment
sources, duplicate dependents, canonical-identity failure, and semantic
ambiguity. Both expose stable rule codes, and `LineageAdmissionDecision` retains
accepted and refused outcomes for structured tracing without printing from the
core library. The ambiguity variant owns the deterministic
invalidation receipt; callers never need to reconstruct the invalidation set
from prose.
`MachineGraphRefusal` owns a nonempty, stable-sorted finding list whose closed
rule vocabulary and typed offending/related subjects are directly suitable for
structured logs. Canonical identity errors are retained only on the identity
rule; no admitted receipt escapes any refusal.

Syntax/study/lowering APIs return `IrError` (span, stable
`IrErrorKind::code()`, detail, hint). Feature-gated planner/anytime APIs return
`PlanError`, and valid but under-budget queries return structured
`PlanOutcome` refusals. `IrNonCanonical` identifies a semantically parseable
versioned artifact whose bytes are not its canonical persisted identity.
`admit` and lowered-receipt binding use the fallible `try_current` boundary;
caller-forged invalid atoms and trees that fit only the bare, not versioned,
depth envelope become deterministic lowering refusals. Neither boundary panics
on malformed caller data.
`DerivedMachineModelCrosswalkCandidateErrorV1` is the feature-gated typed
refusal surface for schema/IR-version drift, zero selectors or artifacts,
raw-to-sealed geometry mismatch, redundant derived-selector mismatch,
cancellation, and canonical identity failure. No partial crosswalk token escapes.

## Determinism class

Machine entity, lineage-record, invalidation, and admitted-graph identities are
bit-stable for the same schema version and semantic inputs. Caller ordering of relations,
targets, and dependents is not semantic. Event kind, target identity, dependent
kind/key/source, entity role, or canonical key changes the appropriate ID.
Machine graph collection order and each subsystem's owned-element order are
non-semantic. Model/card/policy namespace, schema version, semantic digest,
clock declaration, terminal kind/form/dimensions/shape/frame/orientation,
causality, relation endpoint/mode, state slot, port energy role, material
target, and oriented interface endpoint are semantic and move the graph ID.

Parsing, printing, and lowering are pure functions of their input text.
Planner replay is deterministic for the same family, query, ladders, cache
contents, and learned cost table: fixed operator ordering, exact cache-key
framing, coordinate-ordered prolongation, and deterministic tie-breaking.
For the same crosswalk candidate, same sealed derived object, and same current
IR version, crosswalk admission emits the same typed receipt. Every selector,
mapping artifact, aggregate declaration, and no-authority artifact occupies one
fixed canonical field; there is no map-order or platform-dependent input.

## Cancellation behavior

Parsing is bounded by source size and the depth cap. The feature-gated planner
does perform numerical work. Its operational anytime API is synchronously
stoppable between certified budget rungs: `Stop` prevents later planner work,
allocation, cache insertion, and telemetry. Sub-operator cancellation still
lands with fs-exec integration; no claim is made that a running solve or
verification can be interrupted inside its admitted coefficient-by-cell work
envelope.
Derived-crosswalk admission polls at entry, during canonical identity
construction, and immediately before publication. Its work is a fixed 17-field
envelope; cancellation publishes no partial token.
Machine entity, lineage, and graph admission are synchronous bounded metadata
operations (128-byte keys; lineage at most 4,096 relations/dependents and 8,192
target endpoints; graph at most 1,024 clocks/subsystems, 4,096 terminals,
2,048 ports/interfaces, 8,192 relations, 4,096 material bindings, and 16,384
owned elements). They use explicit collection and canonical byte/item
envelopes and do not claim cancellable long-running work.

## Unsafe boundary

None. Safe Rust only.

## Feature flags

- `ladder-planner` [F] (default OFF) — the greedy ladder-walk planner
  (`dep:fs-verify`); disabled until its Gauntlet tier + kill metric are
  green. Gates the `planner`, `plancal`, and `anytime` targets. All
  other current IR behavior is `[S]` default-path.
- `derived-crosswalk` [M] (default OFF) — nominal L6 Machine-IR selectors
  bound to exact admitted L2 derived geometry. It enables
  `fs-geom/derived-geometry` and the fs-exec cancellation dependency.
  `fs-blake3` is now a default dependency for the [S] Machine-IR identity and
  admitted-graph kernel. This feature does not bind its older selectors to an
  admitted `MachineGraphIdV1`.

## Conformance tests

`tests/conformance.rs`: Appendix C spout + frame studies as verbatim
fixtures (names, seeds, locks, lets, and typed-noun counts asserted);
isomorphism property over 200 generated programs plus the fixtures
(s-expr, JSON, and cross-syntax cycles); 8000-parse garbage battery with
in-bounds spans and non-empty hints plus 100k-deep nesting rejections;
span-accuracy cases (bad seed, bad quantity); verb lowering explicitness,
trace content, idempotence, forged-AST validation, and structured refusal;
version-pin round-trip, exact envelope-slot errors, strict canonical persisted
bytes, RFC 8259/control/surrogate cases, and tagged-object uniqueness through
both syntaxes.

`tests/query.rs` (suite `fs-ir/query`, addendum Proposal 8): the wedge QoI
menu is expressible with correct metadata; `value_dims` follows the
functional; well-posed queries admit; the FIVE ill-posed classes each reject
on a distinct check with a teaching fix (zero budget, past deadline, 100%
confidence, field-absent-from-design, self-contradictory dimensions), plus
off-dimension exceedance thresholds and integral-tolerance-needs-volume-dims;
multiple faults are reported together; admission is deterministic (identical
  verdict on replay); and every QoI/target combination round-trips through both
  versioned syntaxes with exact nested float bits and identical canonical
  identities, with a teaching error on a non-query form.

`tests/admission.rs` (suite `fs-ir/admission`): ad-001 Appendix C admits
cleanly + ms-class latency + determinism; ad-001b/c raw/explicit lowering
equivalence and atomic malformed-shorthand refusal; ad-001d forged-AST,
raw-envelope-depth, and lowered-envelope-depth binding refusals without panic;
ad-002 five-study violation zoo
(all rejected on the right dimension, fixes attached); ad-002b malformed,
negative, non-finite, empty, and duplicate resource grants/pillars fail closed,
and self-contained operator grants constrain the study; ad-003 dimensional
spans pinpoint the offending operand, products stay legal; ad-004
BudgetInfeasible with ranked cost-derived fixes + fix-quality harness
(applying fixes admits); ad-005 Router-backed feasibility; ad-006 regime
gating with alternatives + policy grading; ad-007 2000 mutants + all
truncation prefixes never panic (a fuzz-found scanner panic became a
structured refusal).

`tests/planner.rs` + planner unit tests (`ladder-planner`, G0/G3/G5): existing
accuracy, kill-ratio, cache, refusal, calibration, canonical family/cache
identity equivalence and semantic-mutation separation, and replay checks plus
empty/zero/non-monotone rungs; non-finite theta/tolerance/budget; malformed
family/mesh/candidates; poisoned cost samples; unaffordable initial solves;
independent replay of a falsely certified cache answer; non-dyadic
prolongation; adaptive-to-uniform coordinate interpolation; pessimistic learned
costs that cannot veto affordable exact work; bounded family/rung/cell and
combined coefficient-by-cell resource drivers; pre-allocation budget refusal;
verifier authority retained on every finite audit bound; and aborted
transitions that do not enter observed telemetry.

`tests/anytime.rs` (`ladder-planner`, G0/G5): monotone verified trajectories,
priced refusal/hints, cache termination, empty/zero/non-monotone budget/rung
ladders, malformed scalar/hint inputs, and explicit no-interval/no-color output
when the final budget cannot fund one solve. A counting-cache regression proves
that an entire budget ladder executes the planner once. Operational observer
regressions prove callback order, actual rung/spend receipts, contemporaneous
hints without future telemetry, verifier-family/flux identity retention, and
that `Stop` prevents later work telemetry and cache insertion while retaining
telemetry for a completed speculative transition.

`derived_crosswalk` unit tests (`derived-crosswalk`, G0/G3/G4/G5): schema and
receipt replay, exact accessors, independent identity movement for every
Machine-IR selector, derived selector, mapping artifact, aggregate declaration,
and no-authority field, including direct identity movement for the otherwise
version-gated IR-language field; stale schema/IR version refusal; raw-to-sealed geometry
and subject/version/frame/unit mismatch refusal; all-zero refusal for every
selector/artifact role; and entry, identity-construction, and pre-publication
cancellation with no token publication.

`tests/machine.rs` (Machine-IR E0 PR-1, G0/G3): repeatable and role-separated
entity IDs with pairwise schema separation across all six durable roles;
canonical-key refusal with exact diagnostic positions and boundary admission;
caller-order-invariant remesh receipts and diagnostics; structured accepted and
refused decision records; unique dependent rebinding; split recording without
attachments; fail-closed ambiguity binding complete relations, considered
dependents, and cache/contact/winding/adjoint invalidations; split, merge,
fracture, and wear laws; duplicate and public resource-limit refusal; a
maximum-endpoint record plus maximum-relation/dependent invalidation identity
envelopes; and identity movement by event, target, dependent class, complete
event context, and considered inputs.

`tests/machine_graph.rs` (Machine-IR E0 PR-2, G0/G3): structured successful
and refused decisions; full-collection permutation invariance; graph-local ID
schema separation; pressure/stress and absolute/difference-temperature
separation despite equal dimensions; scalar-form, periodic-phase, causality,
clock, frame, orientation, port-shape/power, interface-relation, source-closure,
state-ownership/writer, algebraic
cycle, material, interface, energy-role, and public-resource refusals; stateful
and explicit-policy cycle breaking; and independent graph-identity movement by
model version, clock period, terminal semantic kind, solve policy, and external
material/interface card digests.

`tests/machine_semantics.rs` (Machine-IR E0 PR-3, G0/G3): fully populated
behavior-overlay admission; complete state/initial/boundary/body-motion closure;
quantity, shape, clock, frame, causality, event-history, guard-dependency,
event-clock, reset-write, tolerance, and dependence refusals; canonical
signed-zero scalar handling; public-resource refusal before graph work;
collection and nested-set permutation invariance; refusal permutation
invariance; role-separated event/tolerance IDs; and identity movement by base
graph, guard artifact, and adjacent finite tolerance values.

`tests/machine_assurance.rs` (Machine-IR E0 PR-4, G0/G3): exact admitted-V&V
case and receipt closure; graph-visible and experiment-only sensors; instrument
mapping; aggregate Context-qualified QoIs; hazard assumption/fault coverage;
signed accounting targets, balance dimensions, intervals, and loss ownership;
exact graph-model baselines; terminating same-subsystem fidelity escalation;
public and aggregate-V&V resource refusal; rich outer/nested permutation
invariance; behavior/graph mismatch refusal; and identity movement through an
admitted V&V receipt, escalation trigger, and fixed replay reference.

`tests/machine_lowering.rs` (Machine-IR E0 PR-5, G0/G3/G4/G5): complete
durable-source enumeration; one-way scenario/external-domain crosswalk
admission; missing/duplicate/foreign source, invalid locator,
orphan artifact, invalid scenario, and pre-cancellation refusals; caller-order
invariance; independent identity movement by scenario bytes, external artifact
hash, mapping law, and projection policy; exact scenario/manifest/payload
replay; portable `fs-package::SemanticWitness` JSON transport with tamper
refusal; and deterministic `fs-ledger` replay/mismatch evidence. Package,
ledger, witness, plain content, and typed Machine identities are asserted in
their separate domains rather than conflated. Hard scenario-record admission,
empty external selectors, structured decision counts/codes, and outer/nested
portable framing tamper are covered directly.

## No-claim boundaries

- No operator catalog or per-operator semantic versions — gp3.6; the
  `IR_VERSION` constant covers the language only.
- JSON `\uXXXX` escapes decode scalar values and valid high/low surrogate
  pairs; isolated or malformed surrogates refuse with a structured error.
- The verb table is v1-small (optimize-shape, simulate-pour); verbs are
  data to extend, not a framework.
- Qty literals must be written in units fs-qty accepts; information
  units are Counts, not quantities, by design.
- Admission's dimensional pass covers arithmetic/comparison heads;
  verb-signature dimension contracts (per-operator expected dims) land
  with the operator registry.
- Chart requirements are supplied by lowering/callers; admission does not
  yet derive them from raw study text.
- Machine behavior admission is structural. Opaque condition, motion, guard,
  reset, distribution, tolerance, correlation, witness, and no-claim references
  are bound exactly but neither authenticated nor executed here. Admission does
  not prove state-model compatibility, PDE/DAE well-posedness, swept geometry,
  true-flow event coverage, root completeness, reset regularity, saltation
  derivatives, absence of grazing/simultaneous/Zeno ambiguity, covariance PSD,
  cross-clock simultaneity/synchronization, or physical/statistical validity.
  Domain execution and one-way lowering into `fs-scenario`, motion, time, and
  UQ artifacts remain PR-5 work. PR-4 sensors/experiments, hazards/faults,
  ContextOfUse, accounting, and fidelity/escalation policy are a separate
  admitted overlay and do not strengthen this PR-3 receipt by implication.
- Machine-assurance admission is structural. Exact V&V receipt binding does
  not authenticate a laboratory, calibrate an instrument, or turn synthetic
  or code-to-code evidence into physical validation. Hazard admission does not
  establish FMEA/FTA completeness, fault probability, containment,
  reachability, safety approval, standards conformance, or regulatory
  certification. Accounting declarations do not prove conservation,
  passivity, loss correctness, or numerical balance closure. Fidelity
  declarations do not validate applicability, rank model accuracy, prove
  crosswalk commutation, guarantee evidence monotonicity, choose an optimal
  model, or authorize runtime promotion. Opaque references are bound exactly
  but not executed or authenticated here. Runtime routing and scientific
  execution remain outside Machine assurance; PR-5 can bind exact projection
  artifacts but does not retroactively strengthen this receipt.
- Machine-domain projection is explicit crosswalk admission, not automatic
  semantic inference. PR-2--4 model, material, value, history, distribution,
  motion, guard, reset, correlation, evidence, and policy references do not
  contain enough concrete data to synthesize domain values. The caller owns
  materialization and supplies a versioned mapping law; PR-5 binds those exact
  inputs and validates the resulting scenario but neither authenticates the
  external law/artifact nor proves semantic equivalence, physical fidelity,
  inverse reconstruction, crosswalk commutation, execution correctness, or
  evidence-color promotion.
- The current FrankenScript codec covers admitted Machine graph and behavior
  syntax, including the behavior-to-graph identity binding. The assurance
  overlay still has no literal AST grammar, and assurance requires
  caller-supplied admitted V&V authority that must never be minted from
  embedded hashes. PR-5 therefore proves `FrankenScript -> admitted Machine
  graph -> admitted Machine behavior` plus the separately constructed
  `admitted Machine stack -> explicit projection -> canonical scenario/domain
  package -> exact transport/replay`; it does not yet claim one uninterrupted
  literal FrankenScript-to-Machine-stack leg. `fs-package` structural verification proves
  exact witness transport, not family semantics, and ledger replay compares
  recorded operation/output identities rather than re-executing the lowerer.
  `fs-package::SemanticWitness` independently caps witness payloads at 256 KiB,
  so only projections whose portable payload fits that budget are packageable
  through that witness type; PR-5 makes no universal packageability claim.
  The standalone decoder treats artifact and crosswalk row bodies as opaque
  framed bytes: it proves exact transport and scenario/hash binding, not row
  grammar, canonical ordering, mapping-law authority, or crosswalk semantics.
  Full package/ledger receipt migration remains the separate dependent Bead.
- Router certification is currently a validated declaration on
  `ConverterSpec`, not an authenticated checker/ledger receipt. Admission
  refuses explicitly estimated routes, but full opaque admitted-converter
  authority remains part of the scientific-evidence migration; callers must not
  interpret the declaration Boolean as independent proof.
- `SessionCapability` is admission's view of a token; issuance,
  revocation, and idempotency keys are fs-session's bead (gp3.7).
  A self-contained `(capability ...)` clause supports static planning and
  source-level admission only; it does not mint runtime authority. Plan §11.3's
  session token remains mandatory before execution.
- IR v2 changes Count identity from binary-float-backed count atoms to exact
  integer/decimal atoms. V1 canonical artifacts must be reparsed and re-emitted
  under v2 before their new identity is recorded; no silent v1 hash migration
  is claimed.
- IR v3 changes quantity semantics from five SI base exponents to six
  `[m, kg, s, K, A, mol]`. V1/v2 envelopes refuse at the persisted boundary;
  callers must explicitly reparse and re-emit legacy source under v3 so the
  new semantic identity is visible rather than silently reusing an old hash.
- Bare `sexpr::parse`/`json::parse` intentionally do not infer an artifact
  version. Persisted or replayed programs must use `VersionedProgram`; callers
  that ledger a bare AST have no version-binding claim.
- The `derived-crosswalk` token does not prove that any nominal Machine-IR
  selector identifies a real, admitted, or canonical model. It does not decode
  or inspect an admitted `MachineGraphIdV1`; execute the
  subject/version/frame/unit mapping artifacts; prove semantic or physical
  preservation; provide an inverse or composition law; construct a
  `VersionedProgram`; convert into a
  derived morphism/equivalence; or transport/strengthen evidence. PR-1 supplies
  strong entity identities and PR-2 supplies a structural admitted graph, but
  the feature's older selector bytes do not automatically become those IDs or
  that graph identity.
  A successor crosswalk must explicitly bind the strong IDs before those
  questions can be admitted. The mandatory no-authority artifact records this
  boundary; it is not a proof by itself.
- Machine-IR PR-2 is a structural graph admission boundary, not an executable
  multiphysics model. Matching quantity dimensions, frames, orientations,
  clocks, effort/flow power dimensions, and interface roles does not prove
  conservation, passivity, constitutive validity, geometric compatibility, or
  scheduler behavior. Algebraic-cycle detection is not a DAE-index proof, and
  naming a solve policy does not prove convergence. Opaque model/material/
  interface/policy references are not authenticated or inspected. PR-3 owns
  IC/BC/motion/event/reset and tolerance/correlation declarations; PR-4 owns
  sensors, hazards, ContextOfUse, accounting, and fidelity policy; PR-5 admits
  caller-materialized scenario/domain crosswalks and stable-ID/hash round
  trips without inventing executable meanings. Lineage persistence and
  authenticated crosswalk authority are also absent. Admission-decision
  summaries are not canonical digests of early-refused drafts and are not by
  themselves replayable ledger records.
- The query language is v0: a FIXED QoI menu (max/integral/exceedance), not
  a general program surface. `Query::admit` type-checks well-posedness and
  dimensions ONLY — it does NOT plan, cost, or execute a query (the greedy
  fidelity-ladder planner and the anytime/refusal result semantics are
  separate addendum beads). Field dimensions come from a caller-supplied
  `FieldRegistry` (the design's typed fields, Proposal 13); this module does
  not itself derive fields from geometry. `budget_usd` is a priced dollar
  budget distinct from the wall/memory/core grants of the `(budget …)` study
  clause. The returned answer's COLOR (verified/validated/estimated) is
  attached by the query result, not here.

## No-claim boundaries (planner)

- v0 discharges the verifier's 1-D elliptic kernel class; the 2-D
  cutfem DWR (fs-dwr) and real physics kernels plug into the same walk
  as the ladder registry grows rungs.
- Cost units are solved cells (the flywheel's telemetry currency);
  wall-clock costs arrive with the perf-CI lane.
- The v0 synchronous numerical envelope additionally refuses when polynomial
  coefficients times mesh cells exceeds `MAX_POLYNOMIAL_CELL_WORK`; this is a
  deterministic resource guard, not a wall-time certificate.
- Cache storage/transport authentication remains the content-addressed store's
  responsibility. This planner treats cache data as untrusted and re-verifies
  its numerical claim, but does not authenticate who wrote the entry.
- The v0 family boundary checks finite polynomial structure and homogeneous
  endpoints; it does not prove that arbitrary caller-supplied polynomial
  semantics represent the intended physical model.
- Confidence targets (`Target::Confidence`) are the e-process beads'
  contract; v0 discharges tolerance targets.
- The kill measurement (>=2x vs mid-rung+uniform; measured 4.31x on the
  steep-feature fixture) is per-fixture evidence, not a universal claim
  — the wedge query set re-measures it as kernels land.

## No-claim boundaries (anytime)

- The hint's price is an O(h) extrapolation from the achieved bound —
  an estimate for teaching, not a certified cost bound; Proposal C's
  full value-of-information ranking replaces it when C lands (the soft
  dependency the bead names).
- Operational interruption is rung-granular. A callback can stop before the
  next operator and receive a clean deterministic prefix, but sub-operator
  cancellation lands with the fs-exec tile integration.
