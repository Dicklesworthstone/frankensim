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
  caller order. Downstream cache/contact/winding/adjoint/manufacturing-state
  attachments are
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
- `machine::causalization` (I02.1 equation-variable hypergraph and receipt
  schema, [S]) is an additive post-Machine admission boundary. A
  `CausalGraphDraft` binds one exact `AdmittedMachineGraph`, an optional exact
  `AdmittedMachineBehavior` when state coordinates are present, an SI-base
  unit convention, complete or explicitly partial extraction scope, and an
  extractor/budget/capability/seed/determinism context. Equations, variables,
  activation conditions, and incidences are bounded and canonicalized rather
  than assigned collection-position IDs. `EquationId`, `VariableId`, and
  `IncidenceId` retain complete canonical identity receipts behind shared
  ownership so repeated endpoint references do not discard collision-
  adjudication metadata. Equality, ordering, and hashing use the semantic root
  plus the independent canonical-preimage/schema/metric adjudication tuple;
  nominal-digest duplicates are still terminal ambiguity. Adjudication-aware
  equality propagates through admitted graph, receipt, theorem-claim, decision,
  and migration wrappers, so evidence-only encoder
  limits cannot reappear as user-visible semantic inequality at an outer API.
  Public raw `IdentityReceipt` getters remain evidence-object views whose own
  equality includes encoder limits; callers that need map/set/cache keys must
  use causal adjudication rather than raw receipt equality. The typed public
  key that makes this distinction unambiguous is the required I02.3 predecessor
  `frankensim-leapfrog-2026-program-i94v.1.2.10`.
  Pre-canonicalization resource telemetry is itself a cancellation-polled pass:
  top-level vector lengths remain exact, nested aggregates carry an explicit
  `complete` bit, and cancellation or the first early cap refusal returns a
  partial, visibly
  incomplete count record without beginning graph admission.
  Top-level caps refuse immediately after that first cancellation checkpoint;
  nested counting short-circuits at the first aggregate/per-row cap instead of
  traversing the remainder of a hostile payload.
- Standalone incidence identity limits realize the full public activation
  sub-envelope rather than imposing a hidden smaller byte cap. The graph-
  artifact field limit is independently derived to contain every max-key
  audited-bridge artifact row at the public incidence-count cap. Because an
  `IncidenceSpec` is intentionally inspectable and mutable before admission,
  admission rederives its normalized identity with cancellation-polled row
  materialization and a cancellation-aware encoder. The normalized graph
  preserves recursive schema typing with Machine/incidence child fields and
  separately composes fixed-size full-receipt adjudication tuples for the
  Machine graph, every equation/variable endpoint, and every incidence. The
  provenance artifact does the same for structure, behavior, lineage, and
  incidence rows. Complete receipts remain attached to their public handles;
  parent identities never duplicate a potentially huge child canonical
  preimage, but they also do not collapse composition to a child digest alone.
- The canonical byte/count caps are not yet a peak-memory certificate. The
  current encoder callers retain `Vec<Vec<u8>>` rows and deterministic
  map/set/sort workspaces, and those temporary allocations are not yet
  fallibly admitted against the active operation-memory lease. A legal
  maximum-incidence graph can therefore exceed available process memory before
  returning a typed refusal. I02.1 remains blocked on the allocation-free
  encoder leaf `frankensim-sj31i.52.1.1` and lease-admitted causal workspace
  leaf `frankensim-leapfrog-2026-program-i94v.1.2.9`; neither this contract nor
  an intermediate receipt claims bounded peak memory until those leaves are
  green. The target public envelope is unchanged.
- Canonical graph and receipt ordering uses a deterministic stable index-map
  sort with cancellation checkpoints through index initialization, cumulative
  bottom-up merge work, inverse-permutation construction, fixed-point scans,
  and payload movement. Long fallible lexicographic row comparisons poll
  internally as well. Identity row materialization is likewise polled at
  bounded aggregate strides. Maximum-size node support/activation rows,
  condition rows, derived-lineage rows, incidence rows, receipt domains, and
  conditional-outcome rows precompute their payload capacity so a late legal
  reallocation cannot evade a poll boundary. Within those canonical ordering
  and identity stages, cancellation therefore returns a typed refusal before
  structure, artifact, outcome, or evidence identity publication rather than
  being deferred across one maximum-size legal row or collection.
- The four top-level causal identity families are runtime-typed and have
  independent candidate domains, versions, field layouts, and quotient-law
  tests, but their central identity-governance closure is **not yet certified**.
  A declaration DAG cannot honestly encode the self-recursive and mutually
  recursive schema closure formed by conditional child outcomes and complete
  receipts. The current byte-schema fingerprint also does not prove the actual
  `FIELDS` table plus encoder closure, and the child-authority, nested-exclusion,
  family-local guard, mutation-matrix, generated-registry, and golden migration
  obligations must move through the minimal SCC bundle set and its condensation
  DAG without co-bundling or version-locking acyclic structure/artifact
  families with recursive normalized or full-evidence families. Those
  obligations are the
  P0 I02 blocker `frankensim-leapfrog-2026-program-i94v.1.2.12`; no partial
  owner declaration, central required-ID row, or golden is an authority claim,
  and I02.1 cannot close until that SCC bundle set and its condensation
  governance are green.
- Causal graph identity has two deliberate axes. `CausalStructureIdV1` is the
  normalized, producer-independent equation-variable structure bound to the
  exact Machine graph. `CausalGraphArtifactIdV1` is the provenance-bearing
  artifact identity: it additionally binds the behavior overlay, extraction
  context, generated/derived/audited source lineage, condition-source audits,
  PortSchema crosswalk audits, and clock-bridge audits. Presentation labels
  affect neither identity. Admission requires each caller-declared audit target
  to equal its source coordinate (`source == audited_source`, or the
  projection/bridge equivalent), so the published artifact tuple cannot
  silently name two different coordinates. This equality check does not
  authenticate the opaque audit artifact, inspect its contents, or prove that
  the referenced checker actually audited that source.
- The current `v1` labels name a **candidate runtime grammar**, not a durable
  schema freeze. Early centralized tests minted only ephemeral conformance
  values without recursive governance, migration closure, lease-admitted peak-
  memory proof, ledger persistence, or release authority; those values are not
  durable migration predecessors. Stabilization requires the resource,
  migration, and recursive-governance blockers
  `frankensim-leapfrog-2026-program-i94v.1.2.9`,
  `frankensim-leapfrog-2026-program-i94v.1.2.11`, and
  `frankensim-leapfrog-2026-program-i94v.1.2.12` plus a fresh central
  proof/golden boundary. That ratification must inventory every
  candidate-era semantic addition, including `EmptyProjection`, and either
  freeze it coherently as the first governed version or assign a new version
  with explicit migration. After ratification, every canonical-row semantic
  change requires a version and migration ratchet.
- Every structural incidence names an exact equation, base variable,
  derivative order, solve participation, coefficient dimensions, term signal,
  optional transform/operator, clock relation, and activation domain. Unit
  closure is checked as
  `variable dims - derivative_order * time dims + coefficient dims = residual dims`;
  the complete residual quantity/shape/clock/frame contract must also agree.
  Cross-clock incidences require an endpoint-exact audited bridge. Equation
  role and solve participation remain distinct: matching equations may read
  known inputs, while known-closure/condition-only equations cannot silently
  contribute unknown matching vertices. A condition-only guard occurrence may
  read a variable solved elsewhere without turning that guard edge into a
  matching edge. A retained state with an active unknown derivative is matched
  through that derivative vertex; isolated/base algebraic and port unknowns
  retain their order-zero vertex. A `ModeDependent` base declaration must be
  covered in every active finite-domain cell by at least one concrete
  `Unknown` or `KnownInput` incidence on a matching equation. Auxiliary
  `ConditionOnly` reads may coexist but cannot discharge that solve-status
  totality obligation.
- Activation is canonical finite-domain DNF: each cube is a conjunction with
  at most one branch per condition and a domain is a sorted disjunction of
  bounded duplicate-free cubes. Condition tables bind the complete finite
  branch set and dependency set. Because v1 conditions are global axes with no
  conditional reachability domain of their own, every dependency must be
  `Always` available. A guard-backed condition binds an exact
  `Guard`/`ConditionOnly` equation plus its
  simultaneous/root-solve obligation, requires that guard and its incidences
  to be `Always` available, and requires the guard incidence dependency set to
  equal the declared condition dependencies. An audited external predicate may
  read only always-available known/condition-only dependencies and binds its
  exact source and audit.
  Incidence activation must logically imply both endpoint domains. Admission
  decides this with a deterministic finite-domain symbolic counterexample
  search, not a merely syntactic cube-cover test; the same engine proves total
  `ModeDependent` participation and disjoint switching. It polls cancellation
  and charges domain rows, branches, partial assignments, cubes, selections,
  and branch-table searches against one explicit proof-work budget. Its
  in-place depth-first assignment and backtrack stack retain memory linear in
  the declared condition count instead of cloning an exponential frontier.
  DNF identity currently preserves redundant but non-duplicate cubes, so
  general Boolean-minimal/BDD-canonical identity remains a future theorem and
  canonicalization ratchet rather than an implied claim.
- `CausalizationReceiptDraft` records three separately encoded axes: determination
  (empty-projection/well/under/over/mixed/unknown), generic block-structural rank
  (full-relative-to-min-side/deficient/not-applicable/unknown), and
  conditionality (unconditional/conditional/unknown). Determination and rank
  are not allowed to contradict the same bipartition: Full pairs with
  Well/Under/Over, Deficient pairs with Mixed, and NotApplicable pairs with
  Under/Over, and concrete `EmptyProjection` pairs only with `NotApplicable`;
  an Unknown axis honestly declines its side of this implication.
  The same law is enforced both when constructing a uniform-mode theorem and
  defensively during receipt admission. Matching pairs bind the
  exact admitted `IncidenceId`, equation, and derivative-variable endpoint;
  unmatched sets must be the exact complements. A domain with neither active
  matching equations nor unknown vertices uses concrete
  `EmptyProjection`/`NotApplicable`, allowing an honest off/disengaged mode
  without minting vacuous Well/Full authority. One-sided empty bipartitions use
  Under/Over with `NotApplicable`. Min-side saturation
  proves maximum directly. Any non-saturating maximum claim requires an opaque-
  field, constructor-validated `MaximumMatchingBinding` tied to the exact
  graph, inhabited domain, canonical matching set, complete set-identity
  receipt, certificate, and checker-referenced evidence.
- Receipt identity also has two axes. `CausalOutcomeIdV1` is a normalized
  producer-independent result identity over structure/domain/axes/matching/
  complements/conditional semantics. `CausalizationReceiptIdV1` is the full
  evidence artifact and additionally binds the graph artifact, analyzer Five
  Explicits, typed theorem commitments, unknown-axis reasons/checkpoints,
  evidence state, and normalized outcome child. Structure, artifact, and
  normalized-outcome composition retains both typed child fields and sibling
  full-receipt adjudication tuples. This permits conformant analyzers to compare
  outcomes without erasing their distinct provenance or silently dropping
  collision evidence at a parent boundary.
- A mode-cell receipt requires exactly one branch for every graph condition and
  projects active rows before constructing matching vertices. Hybrid summaries
  carry no union-graph matching. `ConditionalCausalOutcome` can be constructed
  only from an admitted mode-cell receipt, so its graph, assignment, axes,
  unknown-axis reasons/checkpoints, and typed receipt identity cannot be
  substituted independently. Explicit
  conditional coverage binds the exact graph and the complete Cartesian set
  of every declared condition/branch cell, not merely a nonempty unique subset;
  children must have concrete determination/rank axes, and products beyond the
  public child/selection envelope refuse before theorem publication. A future
  reachability theorem may soundly shrink that product, but v1 does not infer
  unreachable cells. A uniform theorem instead binds exact concrete summary
  axes. Both paths require checker-referenced evidence. A
  `Conditionality::Unknown` summary may retain a bound partial set containing
  incomplete or cancelled children, but it carries no coverage commitment and
  cannot enter the complete-coverage constructor. Heterogeneous concrete child
  axes force the corresponding top axis to `Unknown` with
  `NonUniformAcrossModes`; missing or Unknown child evidence is not itself a
  contradictory concrete value. Invalid coverage diagnostics canonicalize
  children and apply fixed resource/foreign/non-concrete precedence before
  reporting canonical indices, so caller permutation cannot select the error.
  Cancelled or budget-exhausted unknown axes require a deterministic resume
  checkpoint; other reasons forbid one. Receipt resource telemetry uses the
  same cancellation-aware completeness
  convention as graph telemetry and applies outer caps before walking or
  cloning nested child assignments. Coverage and maximum-matching binding
  constructors likewise cap child sets, match sets, and mode assignments before
  cancellation-polled cloning or canonical sorting. Their private theorem-set
  commitments retain and revalidate complete identity receipts, not only
  semantic roots. Public child projection likewise copies a bounded mode
  assignment and its bounded progress state under an explicit `Cx`.
- Invalid-draft diagnostics are themselves resource-bounded. Graph and receipt
  admission retain at most `MAX_CAUSAL_GRAPH_FINDINGS` and
  `MAX_CAUSAL_RECEIPT_FINDINGS` detailed rows respectively; crossing either
  public budget returns one deterministic `ResourceLimit` sentinel and no
  identity instead of allocating attacker-controlled millions of findings.
- `CausalSchemaMigrationDraft` retains explicit schema lineage without moving
  the native target identity. Historical receipts and migration drafts have
  private fields; historical family/version/identity/preimage/schema/frame/
  field/collection metadata are constructor-gated with nonzero-digest and
  canonical-metadata completeness checks, and the target is a private typed
  native receipt enum. Admission therefore derives rather than trusts target
  family, schema version, identity, preimage, and receipt metrics, and only
  admits strictly older same-family predecessors. Admission and identity
  publication accept an explicit `Cx`, checkpoint before and after minting, and
  publish nothing after observed cancellation. It does not authenticate the
  truth of caller-supplied legacy metrics. V1 currently covers the top-level
  structure, graph-artifact, and causalization-receipt families; the required
  child/outcome/matching migration closure is tracked by
  `frankensim-leapfrog-2026-program-i94v.1.2.11`, and I02.1 cannot close before
  that typed family closure lands.
- Graph and receipt admission require an explicit `fs_exec::Cx`, poll at
  bounded validation/identity boundaries, publish no identity after observed
  cancellation, expose structured decision counts/codes, and keep all library
  output silent. Resource caps cover outer collections plus aggregate supports,
  condition dependencies/branches, DNF selections/cubes, derived parents,
  labels, matching vertices, conditional children, and conditional selections.
  Structure identity, graph-artifact identity, unmatched-equation, and
  unmatched-variable failures retain distinct diagnostic subjects.
  Invalid/non-incidence or endpoint-duplicating matching pairs never enter the
  matched witness sets, so they cannot suppress complement or downstream axis
  diagnostics.
  Duplicate nominal identities are a terminal ambiguity and refuse before an
  ID-keyed map can select a caller-order-dependent payload.
- `machine::codec` and `machine::assurance_codec` are the additive
  FrankenScript-to-Machine graph/behavior/assurance-v1 boundaries. One
  current-version `VersionedProgram` body decodes through a closed positional
  grammar into an authority-free draft; only the existing graph, behavior, or
  assurance admission may publish semantic identity. Behavior syntax
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
  references. Assurance syntax binds exact graph and behavior identities and
  carries complete V&V receipt commitments, but those rows are never evidence
  authority: readmission requires caller-supplied `AdmittedVvCase` values,
  verifies them through the existing assurance boundary, and only then compares
  the derived bindings exactly. Its derived envelope is 5,832,832 nodes and
  256 MiB aggregate atom text with a 256-byte per-atom limit, with every PR-4
  top-level and 65,536 aggregate nested cap preflighted. Each admitted-value
  writer runs the same preflights. All resource passes precede recursive AST
  validation and decoded-vector allocation.
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
- `machine::interop` is the additive E7 workflow/interchange seed. A
  `MachineWorkflowDraftV1` must bind one exact admitted Machine graph and the
  ten plan-7.9 engineering stages exactly once in canonical order. Every stage
  retains a bounded, nonzero, versioned content-hash coordinate, and
  `MachineWorkflowIdV1` binds the current `IR_VERSION`, graph identity, stage
  tags, and complete ordered artifact coordinates.
- `ForeignExecutionDraftV1::admit_against` records a bounded FMI 3.0.2 or SSP
  2.0 boundary against one exact admitted workflow and Machine graph. It binds
  caller-supplied model-description, adapter, isolation-receipt, and opaque
  output-artifact coordinates; closes every output over a declared terminal or
  state slot; rejects duplicate names or targets; and canonicalizes output
  order. Foreign results are sealed at `ColorRank::Estimated` with infinite
  dispersion and a mandatory no-native-authority policy marker. Requests for
  `Validated` or `Verified` rank refuse instead of being demoted or promoted,
  and the admitted type exposes neither `Color` nor `AdmittedColor`.
- `machine::manufacturing` is the additive E7 as-built process-lineage seed.
  `MachineManufacturingDraftV1::admit_against` binds one exact admitted Machine
  graph, one exact external process-tolerance correlation-model coordinate,
  and at most 4,096 durable-body process steps. Each step binds a stable key,
  declared casting/forging/machining/additive/heat-treatment/coating/assembly
  family, predecessor, process specification, input material state, and the
  resulting microstructure, residual-stress, and property-state artifact
  coordinates. Admission requires every body to exist in the graph and every
  body's submitted steps to form exactly one finite linear predecessor chain;
  duplicate/dangling/cross-body IDs, forks, cycles, disconnected histories,
  and zero artifact hashes refuse before identity publication. Caller step
  order is non-semantic; canonical body-chain order is retained.
- An admitted manufacturing state exposes one `ManufacturingState` lineage
  dependent per process step. The existing `LineageRecord` law therefore
  transports the complete process history across a unique one-to-one body
  morphism and returns a complete `LineageInvalidation` on an ambiguous split
  or remesh instead of choosing a descendant.
- `machine::manufacturing::tolerance_axis` is the additive one-way L6 bridge to
  `fs-toleralloc`. It binds one exact admitted behavior and manufacturing-state
  identity to the complete non-forgeable `CorrelatedStackReceipt`: external
  model namespace/version/digest, dimension, every lower-factor binary64 bit,
  row-norm defect, every ordered term/name/sensitivity/color/standard-deviation
  bit, and all five published first-order moment bits. Machine dependence order
  is authoritative; every positional term label must exactly equal its durable
  `ToleranceId` key. V1 admits only tolerance-only correlated groups whose
  retained specifications are random, scalar, body-targeted, and attached to a
  body with manufacturing process history. It rejects mixed random conditions,
  vector/tensor axes, subsystem/surface/feature targets, and missing history
  rather than inventing a projection, component basis, or owner-to-body map.
- Manufacturing correlation references are content coordinates while the
  `fs-toleralloc` model carries a caller-supplied semantic digest. The crosswalk
  does not equate those 32-byte domains. It retains a separate nonzero,
  versioned `correlation-coordinate-link` artifact naming the caller's explicit
  mapping policy, and makes that policy identity-semantic.
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
- A version-one admitted engineering workflow has exactly ten steps in the
  declared plan order and is inseparable from its exact Machine graph. Stage
  tags and every artifact namespace, schema version, and content hash are
  identity-semantic; no partial workflow publishes a receipt.
- A foreign execution receipt is inseparable from its graph, workflow,
  standard version, model-description, adapter, isolation coordinate,
  no-authority policy, color-algebra version, and complete output bindings.
  Output names and associated targets are unique, targets must exist in the
  graph, and caller output order is non-semantic. Its only representable
  evidence rank is `Estimated`; no foreign output inherits native certificate
  authority.
- A manufacturing-state receipt is inseparable from its graph, manufacturing
  schema and FrankenScript IR versions, exact correlation-model namespace/
  version/content hash, canonical per-body process order, process families, and
  all five artifact coordinates carried by every step. A step may name only a
  graph-owned `BodyId`; its predecessor must name exactly one earlier step on
  that same body. Manufacturing dependents obey the same unique-rebind/
  ambiguous-invalidate law as other durable Machine attachments.
- A manufacturing-tolerance crosswalk is inseparable from its graph, behavior,
  manufacturing state, explicit correlation-coordinate link, exact factor and
  stack receipt, and ordered `(full model position, ToleranceId, BodyId)` rows.
  Every bound body has retained manufacturing history. Generic lineage records
  move or invalidate its `ManufacturingToleranceAxis` attachments; they do not
  mutate the crosswalk or mint a successor. Every graph/body/behavior/
  manufacturing transition requires explicit crosswalk readmission. Partial
  ambiguity invalidates only attachments on ambiguous source bodies, but the
  old aggregate receipt remains bound to its original endpoints and cannot be
  reused as a successor.

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
`InteropReferenceErrorV1`, `MachineWorkflowRefusalV1`, and
`ForeignExecutionRefusalV1` preserve stable rule codes and deterministic repair
hints for malformed coordinates, graph/workflow rebinding, incomplete or
reordered stages, empty/oversized/ambiguous outputs, unknown targets, evidence
laundering, and canonical-identity failure. Public collection limits are
checked before per-output validation, and no admitted workflow or foreign
receipt escapes a refusal.

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
Workflow replay is bit-stable for the same schema/IR version, graph, stages,
and artifact coordinates; stage order is semantic. Foreign-output caller order
is normalized by canonical name/target order, while every standard, artifact,
target, policy, and output coordinate is semantic and moves the receipt.

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
Workflow and foreign-output admission are likewise synchronous bounded
metadata operations: ten workflow stages and at most 4,096 opaque output
bindings. They execute no adapter or FMU and make no cancellation-latency claim
beyond these fixed public envelopes.

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

`tests/machine_causalization.rs` (I02.1, G0/G3/G4/G5): minimal graph/receipt
admission with complete node, structure, artifact, normalized-outcome, and
evidence-receipt metadata; explicit candidate four-family domains/versions/
top-level field layouts and representative semantic/provenance/analysis/
evidence quotient mutations over both typed IDs and canonical preimages; exact
state-contract/behavior-overlay binding,
missing/foreign behavior refusal, and normalized-structure versus behavior-
provenance separation; outer-collection permutation and presentation-label
invariance; normalized structure versus extractor provenance; exact incidence/
derivative matching, atomic duplicate-endpoint rejection, and complement
closure; concrete empty-projection axes without vacuous rank authority plus
directed equation-only and variable-only bipartition semantics;
graph/domain/witness/checker-bound maximum-matching commitments; exact finite-
domain `ModeDependent` coverage; admitted mode-cell projection; concrete empty
off-mode aggregation; partial Unknown child resume-state retention without a
coverage claim; mixed concrete/incomplete summaries without invented
nonuniformity plus agreeing/disagreeing concrete-only aggregation over a 2x3
partial cell set; typed hybrid children and exact coverage roots; deterministic
foreign/non-concrete diagnostic precedence under invalid child permutation;
asymmetric 2x3 Cartesian completeness under
shuffled caller order; missing/duplicate/foreign/non-concrete child refusal;
explicit Cartesian-envelope refusal; heterogeneous child-axis honesty;
unrelated-checker substitution refusal; duplicate-identity ambiguity; orphan/
matching guard refusal; foreign ownership, derivative-order, unit, and clock refusals;
single-row resource bombs before canonical sorting; oversized outer graph
collections before nested telemetry; exact finite-domain DNF consensus plus a
missing-branch counterexample; exclusion of auxiliary condition-only reads
from mode-dependent solve totality; refusal of conditionally unavailable guard/
predicate dependencies; exact structure/graph-artifact/causalization migration
metadata plus zero-digest, incomplete-history, family, and version refusals;
pre-cancelled migration publication refusal;
compatible and contradictory uniform-theorem axis tables; oversized
conditional-child and maximum-binding mode-domain refusal before nested
scans/clones/sorts; uninhabitable/wrong-condition/wrong-branch theorem-domain
refusal; nested condition-table/cube/selection permutation invariance; maximum
activation sub-envelope construction plus exact one-over aggregate, cube-count,
and per-cube-selection refusal; and pre-cancelled public incidence, child,
binding, count, graph, and receipt construction before identity publication.
Private G0/G4 unit laws additionally compare the cancellable stable sort with
the standard stable ordering across empty, singleton, poll-stride, reverse,
and duplicate-key cases, then inject deterministic errors at entry, index
initialization, cumulative merge, inverse-map, fixed-point, partial payload-
movement, and completion phases while checking payload conservation. A
long-equal-prefix fallible comparator test separately interrupts inside a row
comparison before payload publication. Graph/receipt refusal assertions
print the typed deterministic finding set so a batch-verification log retains
exact rule and subject context rather than only a Boolean test result;
constructor-boundary assertions retain their exact typed error variants.
Private G0 laws also prove receipt-adjudication Eq/Ord/Hash consistency and
admission revalidation when encoder limits differ only as evidence, while the
inclusive graph/receipt finding-budget boundaries collapse one-over inputs to
a single `ResourceLimit` row.

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
admitted V&V receipt, escalation trigger, and fixed replay reference. Codec
coverage includes canonical s-expression/JSON/version-envelope replay, exact
graph/behavior rebinding, populated nested-vocabulary round trips, receipt-commitment
tamper/duplicate/omission refusal, and proof that serialized hashes cannot mint
the live admitted V&V authority required for readmission.

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

`tests/workflow_interop.rs` (Machine-IR E7 workflow/interchange seed,
G0/G3/G5): exact ten-stage ordering and graph closure; repeatable workflow and
foreign-execution identities; identity movement by stage artifact, standard,
and output artifact; output-order normalization; bounded reference and output
sets; duplicate/unknown target refusal; and explicit refusal of foreign
`Validated` or `Verified` authority.

`tests/machine_manufacturing.rs` (Machine-IR E7 as-built process-lineage seed,
G0/G3/G5): graph-bound casting-to-machining-to-heat-treatment lineage; caller-
order-invariant receipts and canonical predecessor order; independent identity
movement by graph, correlation namespace/version/hash, step ID/body/order/kind,
and namespace/version/hash for every process/material-state coordinate; unknown-body,
duplicate, dangling, cross-body, multi-root, forked, cyclic, disconnected, and
zero-artifact refusal; refusal stability under caller reordering; complete
manufacturing-history `(key, source, target)` rebinding across one-to-one wear
versus exact fail-closed invalidation across an ambiguous body split; and the
4,096-step/128-byte-key identity boundary plus limit-plus-one refusal.

`tests/machine_manufacturing_tolerance.rs` (Machine-IR/fs-toleralloc as-built
axis bridge, G0/G3/G5): exact graph/behavior/manufacturing/coordinate-link
closure; complete factor, term, color, and first-order-moment identity movement;
deterministic replay; deliberate Machine-scale/stack-sigma inequality; exact
body-axis attachment rebind versus ambiguous invalidation; multi-body partial
ambiguity with only affected attachments invalidated; structured refusal of
graph, dependence, model, dimension, mixed-condition, non-scalar, non-body,
missing-history, and positional-name gaps; and the owner's exact 128-axis
boundary without widening it.

## No-claim boundaries

- I02.1 defines and locally validates the current candidate equation-variable
  graph and receipt runtime grammars. Peak-memory/resource closure, complete
  identity-family migration, and recursive schema governance remain explicitly
  open under `frankensim-leapfrog-2026-program-i94v.1.2.9`,
  `frankensim-leapfrog-2026-program-i94v.1.2.11`, and
  `frankensim-leapfrog-2026-program-i94v.1.2.12`; this boundary does not
  extract equations from `fs-opdsl`/`fs-couple`, compute
  a matching, Dulmage-Mendelsohn or BLT decomposition, perform Pantelides/
  dummy-derivative index reduction, choose tears, solve a system, or execute an
  end-to-end mechanism. Those source-owned adapters and algorithms are I02.2+
  consumers of this boundary. In particular, a caller-authored admitted draft
  is not evidence that extraction was complete merely because it names a
  coverage reference.
- Structural rank here is generic block incidence/cardinality, not scalar DOF
  rank, numerical matrix rank, nonsingularity, conditioning, existence,
  uniqueness, stability, convergence, DAE index, physical cause-and-effect, or
  executable scheduling authority. A matching certificate or hybrid coverage
  binding is tied to exact graph/domain/witness/checker coordinates, but this
  crate does not authenticate the checker binary, execute it, or prove its
  theorem. These are current authority boundaries, not assertions that stronger
  admitted maximum-matching, DM/index, EM-topology, or physical-causality
  theorems are impossible; the typed bindings are intentionally shaped for
  those ambitious ratchets.
- The budgeted activation engine exactly decides implication/coverage only for
  this schema's finite declared branch domains and refuses when its explicit
  symbolic-work budget is exhausted. It is not an SMT solver for continuous
  guard predicates and does not authenticate an opaque predicate or root-solve
  obligation. Redundant but non-duplicate DNF cube sets can still have distinct
  structural identities; Boolean-minimal/BDD-canonical identity is not claimed.
- Post-admission causalization cannot rediscover or authorize inter-subsystem
  algebraic loops that `MachineGraphDraft` already refused for lacking an
  explicit solve policy. It analyzes internal/generated equations and loops
  already governed by the admitted Machine graph. A future provisional pre-
  Machine structural stage must be explicit rather than bypassing that gate.
- State derivative matching treats the retained base state as integration
  memory when an active unknown derivative occurrence exists; initial-
  condition completeness remains the exact Machine behavior overlay's job.
  `StateUpdate` does not yet encode a general discrete pre/post/next/delay/stage
  temporal-occurrence algebra, and derivative order alone must not be read as
  one. General temporal occurrence is a successor schema, not an implicit v1
  meaning.
- Partial extraction still binds the complete admitted Machine graph plus an
  exact partial-boundary reference. It does not publish a freestanding Machine
  projection/seam certificate or guarantee incremental identity locality when
  unrelated portions of the parent Machine graph change.
- Causal graph composition retains the complete Machine and behavior identity
  receipts it receives. The upstream `AdmittedMachineBehavior` v1 receipt,
  however, still binds/exposes its base Machine graph by typed ID rather than a
  full base-graph receipt. Collision-adjudicable closure of that pre-existing
  behavior-to-graph edge therefore remains a successor-schema ratchet in
  `machine::semantics`; ordinary non-colliding graph mismatches already refuse.
  Likewise, the upstream Machine-graph identity composes several durable entity
  IDs by digest, so retaining its complete graph receipt here cannot
  retroactively adjudicate a collision already collapsed inside Machine v1.
  Full recursive receipt composition for those upstream identities is a
  separate hardening ratchet, not an I02.1 theorem.
- Canonical schema migration proves exact lineage binding and prevents public
  field forgery; it does not prove semantic equivalence between predecessor and
  target schemas or authorize automatic re-interpretation of legacy bytes.
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
- The workflow/interchange seed is structural receipt binding only. It does
  not parse or export XML/archives, inspect an embedded FMI/SSP edition, claim
  FMI 3.0.2 or SSP 2.0 profile conformance, execute an FMU, enforce an OS
  sandbox, authenticate adapter/isolation coordinates, or recompute a supplied
  content hash from source bytes. It does not establish unit closure, model or
  output validity, physical fidelity, semantic equivalence, third-party
  acceptance, deterministic wire export, or native V&V, safety, regulatory,
  `Validated`, or `Verified` status. Stage artifact coordinates bind the
  caller's assertions but do not perform stage-specific semantic validation;
  this seed also does not yet bind an admitted PR-5 domain lowering or prove
  the full FrankenScript-to-package round trip. Foreign outputs are opaque
  artifact coordinates rather than decoded numeric values. Parser/export,
  isolated execution, and full round-trip batteries remain follow-on work.
- The manufacturing process seed and tolerance crosswalk are structural lineage
  binding only. They do not parse
  GD&T, datum, fit, texture, process, microstructure, residual-stress, property,
  or measurement artifact bytes; authenticate their producer; prove that a
  process ran; validate a correlation population or semantic/content coordinate
  link; prove that a Machine random `scale` is a stack standard deviation;
  establish axis/QoI unit or vector-component closure; derive or validate a
  supplied sensitivity; or promote its evidence color. The crosswalk retains
  the exact already-propagated first-order stack receipt but does not itself
  perform numerical propagation. Neither lane can
  establish physical causality, assembly feasibility, joint/preload/weld/
  adhesive behavior, surface evolution, or material-property correctness.
  Version one is body-level linear process history, not the
  complete datum/GD&T/assembly schema, a nonlinear/hierarchical/mode-switching
  reliability or tail/quantile model, a population-calibration proof, or a
  gear/backlash consumer.
- The current FrankenScript codec covers admitted Machine graph, behavior, and
  assurance syntax, including both assurance base-identity bindings. Assurance
  receipt rows are only exact transport commitments: callers must retain and
  supply real `AdmittedVvCase` authority, and embedded hashes cannot mint or
  substitute for it. PR-5 can therefore consume a fully readmitted literal
  Machine stack before explicit projection, but literal syntax still does not
  materialize the caller-owned scenario or external domain artifacts.
  `fs-package` structural verification proves
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
