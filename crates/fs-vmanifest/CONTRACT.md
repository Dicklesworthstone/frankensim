# CONTRACT: fs-vmanifest

The typed VerificationManifest schema for leapfrog G1 claim/evidence freezes,
with authored draft constructors for I01 through I15 and the CP, EM, PD, and
RL portfolio aggregates.
Sibling G1-freeze beads reuse this schema and add their own instance modules;
CP/EM/RL/PD portfolio aggregates freeze the cross-instance seams.

## Purpose and layer

Layer UTIL (schema + seed data; depends only on `fs-blake3`). A G1 freeze
preregisters claims, fixtures, obligations, and waivers BEFORE any
implementation result is inspected, so tolerances, claim wording,
capability scope, and failure policy cannot drift toward favorable
outcomes. Preregistration is not proof; nothing here mints evidence
colors or promotion authority.

## Public types and semantics

Canonical manifest/obligation encoding is schema v2. V2 introduces unordered,
duplicate-free obligation-set identity and therefore uses new manifest and
obligation hash domains; unchanged claim, fixture, waiver, and authored-spec
component encodings retain their v1 subdomains.

- `ClaimSpec` — one preregistered claim: ambition lattice element
  (`Ambition::Solid/Frontier/Moonshot`), polarity (affirmative or
  refutation/falsifier lane), statement, explicit hypotheses, QoI + unit,
  `ToleranceSemantics` acceptance arithmetic, targeted `GauntletTier`,
  `OracleRoute` (identity, declared independence, declared TCB overlap),
  activation/kill criteria, fallback, and the Unknown/no-claim boundary.
- `ManifestDraft` — mutable assembly authority. Its numeric `version` is
  the sole machine-interpreted manifest-instance revision.
  `FiveExplicits::versions` is opaque, identity-bearing provenance for
  schemas, toolchains, dependencies, and data contracts; freeze hashes but
  deliberately does not semantically parse that prose. Authored instance
  constructors have an additional conformance lint against mirroring the
  numeric revision through a small, explicitly non-exhaustive set of known
  legacy semicolon-field spellings. A future structured-pin schema may remove
  this remaining human-authoring ambiguity without inventing a prose parser.
- `FixturePin` — a corpus element: authored generator spec (digest
  computed from the exact bytes) or external artifact by 64-hex digest,
  each in a `Partition` (development or held-out; the split is frozen).
- `ObligationRow` — one execution leaf's complete mapping: covered claim
  ids, unit-case classes, G0 generators/laws/shrinkers/seeds, deck ids,
  G3 relations, G4 cancellation/fault schedule, G5 determinism matrix,
  `scripts/e2e/leapfrog/*.sh` entry point, smoke/core/max tier, DSR lane,
  fs-obs event kinds, and the exact replay command.
  `claims_covered`, `unit_cases`, `decks`, `g3_relations`, and `obs_events`
  are nonempty duplicate-free sets for identity purposes, so
  presentation-only reordering cannot revoke an otherwise identical leaf.
  The five `canonical_*` accessors expose the prospective lexical view on a
  draft. They sort but do not repair duplicates in an unfrozen raw row,
  because duplicate input is a freeze refusal.
- `FrozenObligationRow` — the immutable accepted projection returned by a
  frozen manifest. Its five owned set fields are lexically sorted during
  freeze after validation; every field is private and exposed only through
  read-only getters, so external code cannot forge or mutate the accepted
  projection. `digest()` returns the exact authored-row component identity.
  Canonical execution/serialization is therefore enforced by the frozen API
  rather than left as a caller convention. This projection changes no
  schema-v2 digest bytes.
- `Waiver` — a NAMED skip: subject, narrow reason, owner, retirement
  predicate, expiry/review point, and explicit promotion effect. There
  are no unnamed skips: a claim must be covered by an obligation or a
  waiver, and an obligation deck must resolve to a pin or a waiver.
- `ManifestDraft::freeze` — fail-closed gates in documented order: all
  collection/per-row-list and cumulative UTF-8 text-byte caps (BEFORE any
  semantic scan), version, top-level blanks, required nonempty claim
  authority, per-component blanks, duplicate ids, oracle independence
  (production-oracle reuse is a refusal, not a style issue), tolerance
  validity, fixture well-formedness, orphan claim references, orphan decks,
  uncovered claims, and unused/misspelled waiver subjects. Each category is a
  whole-manifest pass, so a later-row failure in an earlier category cannot be
  preempted by an earlier-row failure in a later category. Claim
  ids and obligation-leaf ids also share one evidence-id namespace, so
  string-valued amendment invalidation can never collapse two authorities;
  an untyped waiver subject must resolve in exactly one of the claim,
  referenced-deck, or leaf namespaces. Refusals are typed
  `FreezeRefusal` values.
- `FrozenManifest` — SEALED: no public constructor, no mutating API;
  producers are `freeze` and `amend` only, so holding one proves the
  gates ran on exactly this content. `digest()` is the canonical
  identity the freeze is bound to, and `obligations()` exposes only canonical
  `FrozenObligationRow` projections.
- `FrozenManifest::amend` — the only change path: the successor must
  preserve the initiative, carry representable `version + 1`, and pass
  every gate. A predecessor claim id may not become an immediately succeeding
  obligation-leaf id, or vice versa; that adjacent-version authority-kind
  alias is a typed
  `EvidenceKindChanged` refusal. `AmendmentRecord` names exactly the predecessor
  claim and obligation-leaf authority that must be re-earned or explicitly
  rebound: direct claim edits
  reach their producer leaves; fixture/deck, obligation, and waiver
  changes follow reverse dependencies; title, Five Explicits, or
  amendment-rule changes invalidate all predecessor evidence. Filling a
  formerly waived deck slot also propagates. A leaf receipt binds both its
  execution payload and claim mapping: a mapping-only edit invalidates the
  predecessor leaf authority even though its byte-identical execution payload
  may be reused through authenticated amendment lineage. An unrelated sibling
  claim sharing that payload is not invalidated by a claim content edit or
  mapping-only removal/rename, while changing execution semantics or a deck
  invalidates every claim fed by that evidence.

`freeze` and `amend` validate manifest structure and reverse-dependency
invalidation only. A same-ID `External` digest plus waiver removal may
therefore freeze structurally, but never proves governed discharge. Only the
typed fs-vvreg receipt/envelope transaction, verified amendment lineage, and
atomic authority-head advancement grant that separately governed authority.

- `claim_digest` / `fixture_digest` / `obligation_digest` /
  `waiver_digest` — per-component canonical identities (length-framed,
  variant-tagged, exact IEEE-754 float bits, domain-separated BLAKE3).
  `ToleranceSemantics` equality follows the same float-bit rule (including
  distinct valid `-0.0`/`+0.0` bounds); valid external digest hex is compared
  by decoded bytes, so `FixtureSource` and `FixturePin` equality normalize
  hex case exactly as their digest does.
- `i01_draft()` — the I01 (multi-field equation/compiler) instance: 9
  claims (5 baseline [S], 2 [F], 2 [M], including a refutation-polarity
  falsifier lane for the completeness moonshot), 6 fixture pins (2
  held-out), 6 obligation rows covering every claim, 1 waiver for the
  not-yet-licensed external benchmark deck slot. Exposed as a draft so
  consumers freeze it themselves — no panic path hides in a static
  initializer.
- `i02_draft()` — the I02 (machine causalization and structural-index
  compiler) instance on the same schema: 9 claims (5 [S]: incidence,
  deterministic matching/DM/SCC/BLT, scoped index reduction, consistent
  initialization, block plans + repair witnesses; 2 [F]: causal-witness
  minimality/presentation-invariance, hybrid-mode structural
  completeness; 2 [M]: certified hidden-constraint discovery and a
  refutation-polarity globally-optimal-tearing falsifier lane), 7
  fixture pins (2 held-out), 6 obligation rows, 1 waiver.
- `i03_draft()` — the I03 electrostatic/EQS gate: 16 claims split 8
  baseline [S], 2 [F], and 6 [M]; 22 detailed fixture/theorem-card pins
  with 7 stage-local heldouts; 8 obligation rows split 4 Core/4 Max; and
  1 licensed-industrial-pack waiver. Exact cochain algebra is separate
  from numerical convergence. The baseline also freezes the closed-
  system versus grounded capacitance distinction, locally conservative
  current, class-specific dielectric passivity, quantitative EQS-to-
  Maxwell escalation, total-current/work ownership, and constrained
  held-variable force formulas. For a fixed port with dielectric-outward
  normal and current positive into the field, total current is primary and
  satisfies `I=-integral(J_f+partial_t D).n=dot(Q_free)-dot(Q_transfer)`;
  it equals `dot(Q_free)` only for a blocking carrier boundary. Exact
  reciprocal capacitance structure is
  attached to self-adjoint assembly/terminal-adjointness and oriented incidence
  receipts before interval PSD checks. EQS monitors carry explicit norms,
  positive gross-input denominators, and normalized QoI bounds. Fixed-charge
  force uses an affine charge space and gauge quotient, while radial coenergy
  is allowed only on a certified in-domain segment. Independent maximal lanes
  cover governed candidate-first, sealed-simultaneous raw-block commit/reveal
  IID-lot discharge calibration with exact dyadic-atom sampling and
  compatible-set width normalization; support-typed space-charge/aging;
  global-homeomorphism versus a.e.-injectivity electrostriction routes crossed
  with reduced-Schur versus mixed inf-sup stability; and machine-bound
  stationary condensation or complete filtered pronilpotent cyclic-
  L-infinity/BV force naturality, certified signed refinement defects,
  regularized topology-event jumps, and a refutation-polarity counterexample
  search. The latter targets a cardinality-proved exhaustive microgrammar,
  symbolic full-domain theorems, and non-exhaustive adversarial supergrammar,
  with canonicalization of the full decorated object. Manifest version 1
  deliberately grants neither theorem nor exhaustive-search authority from
  prose: pre-candidate successors must freeze complete machine proposition and
  grammar ASTs, definitions, translations, predicates, and rank/unrank proofs.
  Bare chain or cohomology equivalence is explicitly not claimed to preserve
  equilibria or force.
- `i04_draft()` — the I04 (conservation-defect microscope) instance:
  8 claims (4/2/2 with a counterfactual-completeness falsifier), 8
  fixture pins (2 held-out), 5 obligation rows, 1 waiver.
- `i05_draft()` — the I05 HIL/WCET/fixed-point deployment-twin compiler
  gate: 12 claims split 7 baseline `[S]`, 2 `[F]`, and 3 `[M]`; 15 fixture
  and theorem/exhaustiveness-card pins with 6 stage-local held-outs; 7
  execution leaves; and 1 production-target characterization waiver. The
  baseline freezes bounded `DeploymentTwinIR` admission, exact-rational
  fixed-point range/error containment, reproducible safe `no_std` images,
  non-confusable static/compositional/measured/unavailable timing evidence,
  interval-valued multi-clock HIL causality, typed fault/safe-state behavior,
  and directed bounded-horizon target refinement. Frontier lanes attempt
  complete interference ownership and compositional deadline/age bounds while
  reconciling authenticated measurements without relabeling samples as WCET.
  Maximal lanes separately target machine-checked floating-to-binary trace
  refinement and exact worst-case cycles for a complete finite target model,
  while a refutation lane attacks false certificates. Version 1 prose cards
  mint neither theorem nor exhaustive-search authority: successors must first
  freeze executable proposition/definition/target-semantics or finite-state
  grammar/transition/cost/coverage artifacts and independent checkers.
  Generated controller code is safe Rust; unavoidable boot/vector/HAL/
  interrupt target support is a separate content-addressed, capability-
  allowlisted, audited runtime capsule and remains an explicit TCB rather than
  being hidden inside code generation.
- `i06_draft()` — the I06 material-lot passport and substitution-impact
  gate: 12 claims split 7 baseline `[S]`, 2 `[F]`, and 3 `[M]`; 17 fixture,
  theorem, and finite-decision cards with 7 stage-local held-outs; 7 execution
  leaves; and 1 governed-industrial-lot-pack waiver. The baseline freezes
  non-confusable material/specification/supplier/heat/lot/specimen/process/
  custody identity; source-preserving, method/condition/domain-qualified
  property observations; hierarchical lot/spatial/aging posteriors without
  pseudoreplication; three-valued ContextOfUse substitution; exact selective
  invalidation relative to a declared dependency graph; anytime-valid supplier
  drift; and replayable human/agent-parity decision receipts. Property semantics
  explicitly distinguish conductivity kinds, permeability/permeance/
  diffusivity, magnetic quantities, viscosity kinds, mass/volume latent heat,
  hardness scales, tensor frames, censoring, wetting hysteresis, process/history
  and validity domains. Frontier lanes target calibrated coupled-property
  posteriors and identifiable causal transport. Moonshot lanes target a
  machine-checked transitive impact-completeness theorem and exact global robust
  decision optimality over a frozen finite grammar, while a hidden-mutant lane
  attacks forged provenance and false certificates. Version-1 prose grants no
  theorem, graph-completeness, causal, exhaustive, supplier, or decision
  authority.
- `i07_draft()` — the I07 AP242 semantic-manufacturing round-trip gate:
  15 claims split 8 baseline `[S]`, 4 `[F]`, and 3 `[M]`; 27 fixture,
  profile, theorem, and migration cards with 11 stage-local held-outs; 10
  independently scored execution leaves; and 1 narrowly scoped governed
  real-corpus/profile/tool-matrix waiver. The baseline freezes bounded Part-21
  parsing and canonical writing, occurrence-stable product/assembly/
  configuration lineage, exact B-rep geometry plus certified tessellated
  approximations, PMI/GD&T/datum/surface-texture semantics,
  materials/process/lot passports, harness/
  EWIS connectivity, canonical deterministic export, and typed semantic-loss
  receipts.
  Profile authority is the transitive tuple
  `(ProfileDefinitionSchemaRoot, ProfileKeyDigest, ProfileEntryRoot,
  ProfileRegistryRoot, RequiredAssignmentUniverseRoot,
  ProfileAssignmentRoot)`, not a filename, friendly name, inferred
  `FILE_SCHEMA`, or favorable runtime choice. The definition schema is itself
  a mandatory content-addressed canonical artifact. Canonical definitions
  carry the schema root and every schema-table field exactly once in strictly
  increasing ordinal order; the field count must equal the bound table
  cardinality. Registry, required-universe, and assignment roots use the
  frozen, domain-separated, length-framed hash algebra. The registry freezes
  the exact nineteen-field `ProfileKey` and five-field
  `AssignmentLogicalKey` tables, and the decoded assignment domain must equal,
  not merely contain, the frozen required-key universe. Many logical keys may
  intentionally share one registered profile, but missing or unavailable
  schema bytes, a mismatched schema root, duplicate profile or logical keys,
  extra assignments, omitted assignments, independently sortable role/root
  lists, or a profile absent from the exact registry are `IntegrityFailed`.
  Bootstrap parsing is independently bounded before `limit_policy_digest` can
  be trusted: UTF-8 payloads are capped at 65536 bytes, nested collections at
  65536 items and nested frame payloads at 16777216 bytes, and nesting at depth
  64. The definition schema is capped at 4194304 bytes, the canonical definition
  at 16777216 bytes and 4096 fields, and the registry at 4096 entries and
  67108864 projected bytes. The
  required universe admits at most 1048576 keys and 268435456 bytes; the exact
  assignment projection is capped at 536870912 bytes. Every length, count,
  cumulative size, and arithmetic operation is checked against remaining input
  before allocation or recursion; registry/assignment hashing is streaming.
  The nested-payload cap applies to the bytes inside a frame; the frame's
  eight-byte length prefix is included through checked addition in every
  enclosing projection cap. Each top-level profile-key, definition-schema,
  canonical-definition, registry, required-universe, and assignment projection
  is depth 0; entering a framed nested record or collection increments depth
  once, depth 64 is admitted, and depth 65 is refused before recursion or
  allocation. Concrete
  successors must provide two independent decoders plus exact-payload-cap,
  payload-cap-plus-one, exact-depth-64/depth-65, overflow, truncation,
  duplicate/order, and trailing-byte mutations before profile authority is
  used.

  Frontier lanes target safe opaque-extension capsules, kinematic-pair
  semantics, partial bidirectional edition migrations, and governed real
  industrial interoperability matrices. Moonshot lanes target sheaf/stack
  descent, proof-carrying bidirectional equivalence, and maximal falsification
  while preserving opaque unsupported content without mislabeling it
  understood. Public authored corpora prove only scoped deterministic
  mechanics; licensed conformance, vendor-wide interoperability, regulatory
  authority, and industrial-population validity require their separately
  governed evidence. Governed discharge uses an anti-cycle transaction-intent
  projection with a 25-byte ASCII header plus trailing NUL, an exact
  twenty-field tag/type table, role-addressed protected bindings,
  predecessor head/stage fencing, CandidateFrozen/profile commitments,
  explicit future-output unions, a closed two-role receipt-schema set, a closed
  singleton governed-output-schema set at tag `0x0014`, and required independent
  encoder/decoder KAT and mutation gates. The governed-output authority has one
  non-extensible matrix row:
  `GovernedRealAp242CorpusDischarge=1`,
  `RealizationCommitted=3`, scope
  `i07.governed-real-ap242-corpus-discharge`, target and exact retired subject
  `i07-governed-real-ap242-corpus`, cardinality exactly one, required protected
  binding, and role
  `GovernedRealAp242CorpusDischarge/i07-governed-real-ap242-corpus`. The role is
  derived from the fixed prefix and byte-identical retired subject; callers
  cannot supply aliases. The exact canonical matrix bytes and their digest are
  themselves inputs to the singleton schema-set root.

  The one authorized output schema is the closed 34-required-field
  `I07GovernedCorpusDischargeEnvelopeSchemaBytesV1`. Its exact schema grammar
  freezes ordinal, name, type, requiredness, and constraint AST for protocol,
  stage, scope, slot, waiver, role, intent, authorization, output-schema,
  profile, standard/AP/Part-21, conformance, custody, corpus, expectation,
  semantic, QoI/band, clause/interoperability, public/blind authority,
  CandidateFrozen component, joint-commitment, oracle/TCB, review, revocation,
  semantic-loss/component-result, redaction, and validation-axis fields. The
  17 precommit-authorization binding names are a closed, exact, case-sensitive
  namespace; missing, extra, duplicate, aliased, or differently typed bindings
  fail closed.
  Unknown, optional, omitted, duplicate, reordered, or trailing fields are not
  authorized. The realized envelope itself has a fixed header, exact field
  count, increasing ordinal/length/payload records, fixed-width scalar
  encodings, a domain-separated digest, and no extension field; the realized
  FutureArtifact digest and successor-installed same-ID `External` root must
  equal that digest. `GovernanceCommitted` freezes the exact matrix/schema bytes,
  digests, and set root before candidate execution; `CandidateFrozen` copies
  that root; the intent binds it at `0x0014`; and precommit authorization proves
  the sole FutureArtifact's exact singleton membership. Every future envelope
  record binds its slot, related waiver subject, derived role, and exact
  output-schema digest before its union can encode `Pending`. Realization
  rehashes the independently available schema bytes, decodes the envelope under
  the exact 34-field schema, and checks its committed-intent and
  CandidateFrozen equalities. A permissive schema, matrix extension, role alias
  or swap, missing/extra output, cross-protocol/stage/scope schema, forged
  membership proof, unavailable bytes, or
  `GovernanceCommitted`/`CandidateFrozen`/intent root mismatch is
  `IntegrityFailed`.

  The canonical matrix, schema-set record, membership proof, and realized
  envelope have exact derived lengths 221, 171, 227, and 1458 bytes. Bootstrap
  admission nevertheless preflight-caps untrusted matrix bytes at 4096, schema
  bytes at 4194304, membership proof at 8192, and each UTF-8 payload at 65536
  before exact-byte comparison; nesting is capped at depth 64, field count at
  exactly 34, and schema-set count at exactly one. Counts, lengths, projection
  totals, and depth are checked before allocation or hashing. Concrete
  authority use requires two independent matrix/schema/set/proof/envelope
  encoders and decoders, published exact-byte
  and exact-root KATs, and exact-cap/cap-plus-one, depth-64/depth-65,
  permissive-schema, role, membership, cardinality, stage/scope, root, overflow,
  truncation, order, and trailing-byte mutations. These are authored contracts
  and static policy locks in this crate, not claims that production codecs or
  computed KATs already exist.

  Precommit authorization binds only the predecessor/stage receipt, old
  authority head, both schema-set roots, and intent. The typed receipt, same-ID
  envelope, waiver retirement, verified amendment lineage, and authority-head
  advance then commit atomically through the acyclic intent-to-authorization-to-
  envelope-to-successor-to-amendment/commit-receipt graph. The new authority-head
  digest is derived from the old digest/generation, transaction-intent digest,
  authorization digest, realized envelope digest, final-successor digest, and
  amendment-record digest; its generation is exactly the checked old generation
  plus one. A CAS mismatch, generation overflow, receipt-role/schema mismatch,
  or caller-proposed alternative head changes nothing.

  Execution, claim, and evidence completeness are closed orthogonal axes.
  `Passed`, `Refuted`, and every favorable validation axis require a successful
  execution, complete evidence, valid integrity, admitted applicability and
  support, and their exact frozen predicate or independently reproduced
  counterexample. A complete supported predicate failure is `Failed`;
  `Unsupported` requires complete evidence of exact non-applicability.
  Partial or absent evidence forces `Unknown` and has zero promotion authority,
  while failed, cancelled, timed-out, budget-exhausted, infrastructure-failed,
  or integrity-failed execution cannot carry `Passed`, `Refuted`, or a
  favorable validation axis. Infrastructure and integrity failures never
  become scientific failure or refutation.

  Every leaf exposes request, admission/refusal, start, observation, drain,
  finalization, adjudication, and atomic-publication events. Successful
  publication is one complete atomic AP242/loss-receipt pair. Non-success
  publishes no AP242 artifact and records an explicit absent-artifact
  disposition in at most one atomic terminal/FailureBundle transaction; no
  path may publish one member of a success pair or an unrecorded absence. The
  first accepted cancellation atomically seals child admission, the
  observer-set root/count, the separately reconciled descendant-frontier
  root/count, registration epoch, and exact terminal reservation. Observations
  and descendant exits bind the primary request event, seal epoch, and
  membership proofs; repeated or late requests cannot reopen drain. Once a
  unit observes cancellation, only already reserved bounded cleanup,
  descendant joins, exit receipts, drain, and finalization may continue, and
  none may mint new scientific authority.

  The closed drain-trigger vocabulary distinguishes cancellation observation,
  a verified empty-observer seal, observation timeout, infrastructure failure,
  and non-cancellation drain. It orders effective calibrated time first, then
  the frozen causal ranks `InfrastructureFailure=0`,
  `CancellationObserved=1`, `EmptyObserverSeal=2`,
  `ObservationTimeoutDrain=3`, and `NonCancellationDrain=4`, then causal
  logical sequence and stable identity. Non-cancellation work completion,
  domain failure, budget, campaign timeout, or infrastructure termination
  drains and reconciles its admitted descendant frontier without inventing a
  cancellation seal or observation. The observation deadline is inclusive: an
  observation whose interval upper endpoint equals the deadline remains on
  time, while the first nanosecond after it is the timeout onset. Conservative
  trigger-to-drained and drained-to-finalized latencies use checked outward
  interval subtraction and
  inclusive 2-second Core or 8-second Max caps for each segment. At one
  terminal boundary the exact operational cause order is
  `InfrastructureFailed`, `TimedOut`, `Failed`, `Cancelled`,
  `BudgetExhausted`, then `Succeeded`; malformed canonical input or invalid
  integrity selects `IntegrityFailed` independently. A verified empty observer
  set has zero observation latency at seal but still drains descendants. One
  supervisor clock or a receipt-bound conservative calibration makes all SLO
  comparisons meaningful; overflow, underflow, an inward or invalid enclosure,
  or an incomparable clock domain fails closed.

  Explicit JSONL/flush/run/governor caps, deterministic precreation
  aggregation, a globally checked non-borrowable priority-terminal reserve,
  and fail-closed sink behavior prevent observability from erasing terminal
  evidence. Reservations follow exactly `Provisional -> Released` for refused
  admission, or `Provisional -> Reserved -> Consumed|Released` for admitted
  capacity; finalization emits a settlement root, and crash replay must
  reconstruct the same settlement with neither double consumption nor leakage.
  A separately serviced bounded priority writer keeps request/observation/
  drain/finalization evidence
  independent of ordinary traffic. Admission applies a checked global
  earliest-deadline-first demand-bound test over every concurrent run and
  priority segment, including nonpreemptive blocking and remaining service
  capacity; per-run byte/row capacity alone is insufficient.
- `i08_draft()` — the I08 (evidence-budget co-design planner) instance:
  8 claims (5/2/1 with a robust multi-horizon optimality falsifier), 9
  fixture pins (2 held-out), 6 obligation rows, 1 waiver.
- `i10_draft()` — the I10 constitutive identifiability and optimal coupon
  design gate: 11 claims split 6 baseline `[S]`, 2 `[F]`, and 3 `[M]`; 14
  fixture, counterexample, theorem, and finite-grammar cards with 6
  stage-local held-outs (3 Core, 3 Max); 6 execution leaves; and 1
  instrumented-lab-campaign-pack waiver. The baseline freezes non-confusable
  law/experiment/nuisance/discrepancy schema identity; adjoint sensitivities
  published only inside independent outward-rounded finite-difference
  enclosures; three-valued gauge-witnessed structural identifiability (a
  full-rank numeric FIM alone never mints Identifiable; kinematic/isotropic
  hardening under monotone proportional loading, Prony relabeling, and
  modulus-thickness products are seeded confoundings); profile-likelihood/
  sloppiness practical identifiability with preregistered coverage and
  eigen-direction ownership; conservative three-valued manufacturable
  coupon/environment/sensor admission where feasibility precedes any
  information criterion; and blind held-out design gain scored only on
  declared identifiable combinations. Frontier lanes target discrepancy-aware
  robust design (regret vs a scenario-exhaustive oracle) and anytime-valid
  adaptive law discrimination. Moonshot lanes target a machine-checked
  identifiable-combination quotient completeness theorem and certified
  zero-gap global design optimality over a frozen finite grammar, while a
  hidden-mutant lane attacks false identifiability and design certificates.
  Version-1 prose grants no theorem, completeness, discrimination, or
  design-optimality authority.
- `i12_draft()` — the I12 hybrid mode-automaton compiler instance: 10
  claims (6/2/2, including the Solid grazing false-certificate
  refutation tripwire), 7 fixture pins (3 held-out), 6 obligation rows,
  and 1 waiver.
- `i13_draft()` — the I13 cohomology-preserving electromagnetic topology-
  optimization gate: 17 claims split 6 baseline `[S]`, 5 `[F]`, and 6 `[M]`;
  28 fixture, theorem, route, machine, and adversarial cards with 12
  stage-local held-outs; 14 execution leaves; and 1 governed external
  industrial-validation waiver. The baseline freezes admissible material
  fields and energy-consistent interpolation, relative winding/cohomology
  ownership, exact-discrete complex adjoints, topology-event lineage,
  manufacturability and insulation constraints, robust QoIs, and realizable
  coil-route semantics. Manufacturing, route realization, robust-QoI, and
  multiphysics-closure evidence are independently owned rather than hidden in
  a monolithic success bit. Frontier and moonshot lanes target integer
  cohomology synthesis, differential-character invariants, certified
  topology-change continuation, topology-to-route theorems, and globally
  robust manufacturable machine designs. Public deterministic decks never
  become industrial validation; theorem prose, finite survival, and a numeric
  optimum never mint proof, global optimality, safety, or production authority.
  Governed authority advances only through `GovernanceCommitted`,
  `CandidateFrozen`, `RealizationCommitted`, `RevealedForAdjudication`, and
  `Closed`. Canonical governance framing uses checked fixed-width lengths,
  deliberately nested framing for UTF-8 payloads, exact case-sensitive bytes,
  and no implicit normalization. One governed-group projection binds four
  nonempty, sorted, duplicate-free sets—claims, leaves, fixtures, and
  strata—plus the closed protocol kind under a domain-separated root. Each set
  is capped at 4096 identifiers, each identifier at 65536 UTF-8 bytes, and the
  full projection at 4194304 bytes. Split reveal is forbidden, and two
  independent encoders must reproduce a published byte/root KAT and its
  adversarial framing, order, cardinality, cap, overflow, and trailing-byte
  mutations before group authority is used.

  Transition and protected-access chains have distinct group/protocol-bound
  genesis receipts and authority heads at generation zero. Their first durable
  append consumes the exact genesis, advances generation `0 -> 1`, and every
  later append consumes the immediately preceding receipt/head and increments
  exactly once; reuse, cross-group/protocol/epoch substitution, forks, ABA, or
  overflow are `IntegrityFailed`. The closed governance namespace contains
  exactly `governance.transition` and `governance.protected_access`; exact
  replay returns the original response receipt and never appends another
  canonical outcome.

  Every transition serializes exactly the twenty-five-field inventory
  `source_universe`, `acquisition_access_custody`, `authorized_principals`,
  `candidate_input_permissions`, `selection_exclusion`, `receipt_schema`,
  `decision_rule`, `candidate_freeze`, `schema_slot_successors`,
  `realized_root_commitments`, `hidden_opening_commitments`,
  `protected_bindings`, `retired_waivers`, `discharge_envelopes`,
  `pre_reveal_access_chain_root`, `amendment_record`,
  `realization_authority_head`, `one_shot_capability`, `reveal_receipt`,
  `reveal_access_chain_head`, `adjudication_receipt`,
  `final_access_chain_root`, `action_outcome_effect_counts`,
  `terminal_receipt`, and `failure_bundle_or_clean`. A normative two-protocol
  by five-stage matrix assigns every field exactly one of
  `Digest|Pending|NotApplicable`: a digest stays byte-identical, permanent
  non-applicability never changes, only `Pending -> Digest` is legal, and
  `Closed` has no pending field. The holdout and external-discharge protocols
  use different permanent-not-applicable sets; the matrix may not be inferred
  by carrying a previous row forward.

  Protected access uses a durable write-ahead `Prepared` row, atomic
  access-head/generation compare-and-swap, and exactly one correlated terminal
  row. Only `Prepared` carries `terminal_effect_option=Absent`; `Completed`
  requires `Performed`, `Refused` requires `NotPerformed`, and `Aborted`
  requires exactly one of `NotPerformed|Performed|MayHaveOccurred` according to
  authenticated effect evidence. A crash-dangling `Prepared` row remains
  potential access and blocks authority until conservative durable
  reconciliation; replay never
  appends a second row, and `Closed` requires zero unresolved operations.

  Industrial validation is deliberately two-phase. While the waiver is live,
  only the distinct `/governance` attempt may commit
  `GovernanceCommitted -> CandidateFrozen -> RealizationCommitted`; it has
  `NoPromotionAuthority`, may not reveal protected data or run validation
  science, and must drain and finalize. Only after its transition and access
  frontiers reconcile, its exact receipt and heads verify, and the waiver is
  absent from the committed successor may a separately admitted `/science`
  attempt enter `RevealedForAdjudication`, execute science/adjudication, enter
  `Closed`, and seek promotion. A `RealizationCommitted` transaction durably
  committed before cancellation remains committed but never auto-starts
  science.

  External discharge additionally uses an exact nineteen-field anti-cycle
  transaction intent binding protocol, predecessor stage/head,
  `CandidateFrozen` commitment, indivisible role-addressed protected records,
  mutation fence, future artifact relations, the closed two-role receipt-
  schema-set root, and `GovernanceProtocolSchemaSetRoot` at tag `0x0013`.
  Output-schema authority is not caller-selected. Exact
  `I13_GOVERNANCE_PROTOCOL_MATRIX_V1` bytes contain two and only two ordered
  rows: holdout-slot realization at `RealizationCommitted` under
  `i13.holdout-slot-realization`, with role derived from target slot and empty
  waiver; and external discharge at the same stage under
  `i13.external-discharge`, with role derived from the exact retired-waiver
  subject. Both admit 1..=4096 outputs and require a protected binding. Matrix
  header, row framing, enums, order, cardinalities, scope, role prefix, waiver
  rule, and trailing-byte closure are normative.

  Each role-specific `OutputArtifactSchemaDigest` derives from protocol, stage,
  scope, role, and complete canonical schema bytes. The flat authenticated
  schema-set projection binds the matrix digest plus sorted unique member
  records; a membership proof is exact zero-based ordinal plus framed member
  and has authority only while the complete available projection rehashes to
  the frozen root. An isolated record is not a proof. `GovernanceCommitted`
  freezes matrix and set before candidate execution, `CandidateFrozen` copies
  the root, intent tag `0x0013` byte-equals it, and each protected/future
  artifact schema equals its role-specific member digest—never the set root.
  The matrix or one schema is at most 4194304 bytes, but matrix plus all fetched
  schemas and the set projection are cumulatively capped at 67108864 bytes and
  1048576 decoded constraint-AST nodes, with depth 64 admitted and 65 refused
  before allocation. Permissive schema substitution, matrix extension, role
  alias/swap, cross-protocol/stage/scope reuse, add/drop output, forged ordinal,
  root mismatch, unavailable bytes, or cumulative cap breach is
  `IntegrityFailed`.

  Initial schema choice is not a trusted-governance loophole. The exact 50-byte
  `I13_OUTPUT_SCHEMA_SEMANTIC_AUTHORITY_V1` projection contains two and only two
  protocol-to-schema-kind rows and hashes to the semantic-authority root. Every
  exact 791-byte `I13_OUTPUT_ARTIFACT_SCHEMA_V1` wrapper embeds that root, the
  selected schema kind, the required protocol-payload schema digest, and a
  closed sixteen-field ordinal/name/type/requiredness/constraint table. The
  constraint atoms require byte equality with matrix, intent, authorization,
  future-artifact, protected-binding, payload, and realization-receipt
  authorities; a merely nonzero value is insufficient. A correctly hashed
  schema that starts permissive, optional, renamed, retyped, extended,
  role-aliased, schema-kind-swapped, or payload-schema-substituted is therefore
  a schema mismatch, not an authorized policy choice.

  Holdout outputs must reuse the exact realized-artifact payload schema frozen
  by the predecessor slot/generator contract before `GovernanceCommitted`;
  governance cannot invent a replacement. External discharge outputs use the
  exact 1422-byte, 27-field `I13_EXTERNAL_DISCHARGE_BODY_SCHEMA_V1` descriptor.
  It embeds `NestedSemanticAuthorityRoot`, not caller-selected nested bytes. The
  independently content-addressed 2615-byte
  `I13_NESTED_SEMANTIC_AUTHORITY_V1` projection binds every policy grammar,
  two distinct decoder fingerprints with independent owners/toolchains, and
  positive/negative KAT roots for the twenty nested semantic groups. The
  predecessor fixes both payload-schema authorities; `GovernanceCommitted`
  rehashes and freezes them; only then may `CandidateFrozen` copy the identical
  schema-set root. Independent encoders, decoders, and semantic conformance
  evaluators reject every field, constraint, nested-authority, decoder-set,
  payload-length, ordering, truncation, and trailing-byte mutation before any
  governance authority exists.

  Instance identity is closed as well. `I13_GOVERNED_OUTPUT_ENVELOPE_V1` uses
  exact ordinal/length/value framing for all sixteen fields, and its
  domain-separated digest is the sole `realized_artifact_digest` used by the
  future slot, successor, and `RealizedOutputRoot`. The external-discharge body
  analogously uses exact 27-field instance framing and a streamed payload digest
  that consumes exactly the frozen byte count. The output-schema proof-set root
  covers exactly one full-set-rehash membership proof per future artifact in
  transaction order. The schema-set projection is duplicate-free and all-and-
  only the protected/future derived-role image, so dormant, wildcard, extension,
  and payload-authority-free members have no authority.

  Receipt roles remain exactly `ProtocolAuthorization` and `ProtocolCommit`.
  Their exact schema projections are 697 and 805 bytes; their schema-set
  projection and membership proof are exactly 139 and 34 bytes. Authorization
  and commit receipts are independently decoded and rehashed under their exact
  domains before use. Holdout and external realization receipts have reachable
  role-derived maxima of 131556 and 196978 bytes; authorization receipts have
  reachable maxima of 320 and 314 bytes for the two exact scopes. Larger
  syntactic framing ceilings do not make impossible role/scope combinations
  valid. The exact realized-output set is a role-addressed root over every
  precommitted Pending slot and no other output. The new authority head derives
  from the old head/generation, intent, authorization-receipt digest,
  realized-output root, final successor, and verified amendment, with generation
  checked as old-plus-one. Alias/swap/absence, arbitrary nonzero schemas,
  partial output, stale CAS, overflow, or caller-proposed alternative head has
  zero authority and zero effect. These encodings and KAT suites are authored
  authority contracts and static locks in this slice; production governance
  codecs, independent decoder implementations, and executable computed KATs do
  not yet exist.

  Fresh authority uses the corrected V2 governed-subartifact and semantic-
  completeness closure. Every nonempty semantic inventory row has closed
  requiredness and semantic-role fields, a domain-separated inventory digest,
  and a predicate whose canonical constant-folded form may not be `True`.
  Terminal record keys, key lists, consumed authority receipts, and breach-
  frontier references are distinct semantic types rather than interchangeable
  raw digests. Parent-extractor ASTs have canonical child order, uniqueness,
  depth/node/byte caps, and tag-specific role checks. Each grammar/inventory pair
  binds an exact two-member independently owned evaluator set and nonempty V2
  positive/negative KAT sets; every required atom has both coverage polarities.

  The semantic matrix has the closed scopes ExternalSemanticGroup,
  TerminalEventKind, and HoldoutTargetSlot and contains all and only the frozen
  predecessor keys. Scope2 is derived from structured data only: the 294
  `obligations[*].obs_events` occurrences plus the exact 21-row
  `I13_PROTOCOL_FIXED_TERMINAL_EVENTS_V2` registry give 315 occurrences and 64
  unique event keys. The registry is an exact 1054-byte artifact; its reserved
  `/65535/<ordinal>` source paths cannot collide with ordinary manifest fields.
  Policy prose is never scanned for event authority. A canonical source catalog
  binds the registry root, and its generated member-set root precedes every
  receipt, preventing a hash cycle. Repeated lifecycle-event occurrences are
  not competing duplicate keys: exact field-path/source-record occurrence sets
  aggregate them under one unique `(scope,event)` catalog member, retain every
  preimage, and require the normalized meet of every occurrence authority.
  Ordinary non-scalar atoms carry their structural value kind and no invented
  scalar discriminator; the four specialized causal roles retain their exact
  role-payload discriminants and required Scalar/Sequence/Record shape. Two distinct governed
  generators produce byte-identical member and terminal-event catalog bytes;
  two exact 336-byte generation receipts form a 230-byte receipt-set authority
  that the final matrix itself binds and fs-vvreg independently checks against
  the predecessor. A single generator run twice, caller filtering,
  opaque evaluator/KAT roots, omitted atoms, open scope values, or normalized
  tautologies have no fresh authority. The earlier V1 and incomplete V2
  subartifact phrases remain decode/replay drafts.

  Held-out payload authority is closed end-to-end. Exact V2 generator-contract
  and 278-byte payload-schema projections bind target, alias, index range,
  executable/toolchain, grammar, numeric order, semantic inventory, acceptance,
  label sealing, evaluator/KAT sets, and byte caps. The exact V2 realization
  commitment additionally binds target, alias, predecessor, secret key, contract,
  payload schema, order, inventory, and sealing policy. Contract/order/predicate/
  matrix/receipt fields must byte-equal one another. Realized cases use a fixed
  domain-separated numeric-order Merkle tree with suffix empty-leaf padding,
  height-tagged nodes, a 90-byte root preimage, and at-most-522-byte inclusion
  proofs; no unspecified tree, duplicate-last padding, or schema-less V1
  commitment may authorize a fresh reveal. Infrastructure failure removes
  heldout-stratum and promotion authority without collapsing lifecycle
  completeness.

  The following V3 capacity-profile, physical-service/store, replay-metadata,
  shared control-ledger, bootstrap, expiry, staging, and downstream-binding
  encodings are retained only to decode and audit predecessor artifacts. They
  do not authorize a newly executed I13 attempt; the terminal V4 closure below
  governs every fresh capacity effect. Nonconflicting V3 append-slot,
  lane-content, scientific-retention, and terminal-persistence encodings retain
  their stated authority. The replay-only V3 predecessor used an acyclic
  capacity-profile proposal -> global-reservation intent -> global-reservation
  receipt -> derived lane-plan identities and physical service/store/control
  reservations -> infrastructure-capacity set -> two service-issued lane
  reservations -> admission core -> activation -> context chain. Its V3
  capacity-profile extension is exactly
  435 bytes and is a pure pre-intent proposal: it contains no committed receipt,
  infrastructure-set root, or final lane-plan digest. The global intent and
  receipt are exactly 389/534 bytes, and each lane receipt is exactly 225 bytes.
  Only the committed global receipt authorizes deriving the two plan identities
  and obtaining the physical reservation receipts; the 2542-byte infrastructure
  set then binds those receipts, so no reservation object hashes a successor
  that is needed to construct itself. The final 546-byte admission core binds
  the global receipt, lane authority set, seven-kind 304-byte operation schema
  matrix, four-route recovery registry, exact nine-member infrastructure-
  capacity set, and 4096 total-append ceiling; activation remains 214 bytes.
  Attempt, epoch, governor, reservation,
  services, retention destinations, publication services, calibration, leases,
  and ordered intervals are cross-checked at every edge. The 325-byte cyclic
  draft and intermediate 410/450/482-byte cores have no fresh authority.

  Capacity covers durable bytes, one-in-flight scratch, priority segments,
  breach persistence, phase history, and all control authority before any
  capability is issued. Each lane reserves 268435456 durable bytes, 134217728
  scratch bytes, two 269409356-byte phase-specific breach-copy slots, 10485760
  service-control bytes, and `4194304*p_i` priority-segment bytes. The exact
  lane reservation is therefore `951957656 + 4194304*p_i`. The control reserve
  is exactly `5324902 + 2801937 + 1572964 + 659555 + 126402 = 10485760`
  bytes for transition receipts, phase heads, phase leases, component plans,
  and lane-control closure respectively. Those artifacts retain the complete
  burned-slot authority and all supporting preimages. Their checked maxima sum
  exactly;
  none is charged again to the durable artifact reserve and no unbounded or
  unnamed control index remains. The final lane content consumes at
  most 261447712 bytes and leaves exactly 6987744 non-borrowable safety-slack
  bytes. That slack grants no index, metadata, or artifact authority and is
  released only at settlement; every durable byte must fit a named artifact cap.
  Equality at any named cap passes, and the first excess byte, row, proof,
  tombstone, or index entry refuses before allocation.

  Two lane services contribute `1903915312 + 8388608*p_i` bytes. The two
  publication journals and independent publication verifier contribute 8388608
  bytes together; cleanup/control authority and ordinary evidence contribute
  67108864 bytes each. Thus the final non-replica global-governor equation is
  `run_reserve_i = 2046521648 + 8388608*p_i`. Only `p_i=2..=12` is admitted:
  `p_i=2` reserves 2063298864 and leaves 84184784 bytes;
  `p_i=12` reserves 2147184944 and leaves 298704 bytes; and
  `p_i=13` requires 2155573552 bytes and is refused. Every addition and
  multiplication is checked before allocation. The ordinary-evidence
  67108864-byte hold is an exact V3 artifact with a 98-byte fixed projection and
  framed admission-frozen typed entries. `OrdinaryEvidenceStore` is the ninth
  infrastructure member: its row binds a committed generic physical-capacity
  receipt, and its final artifact root must satisfy role-12 conformance. The
  hold is released only after retained evidence and H3 control closure, never
  merely because the global arithmetic would balance. Each of
  the two pairwise
  owner/implementation/medium/failure-domain-disjoint retention destinations
  separately pre-reserves `run_reserve_i` bytes under its own physical
  cross-attempt store-capacity head. Those replica reservations are real and
  authenticated, but are not dishonestly double-counted inside the distinct
  2-GiB attempt governor. The minimum still exceeds half the governor, so this
  flagship profile admits at most one active I13 attempt; smaller profiles need
  their own proved equation.

  In that replay-only V3 protocol, every non-lane, nondestination physical
  service used the same bounded reservation protocol. Its exact 247-byte
  infrastructure-capacity intent binds
  attempt, epoch, service kind/identity, reserved bytes, capacity policy, common
  clock calibration, validity interval, idempotency, and signature set. The
  exact 121-byte active member binds intent, attempt, epoch, service kind/
  identity, reserved bytes, and policy. The replicated exact 228-byte capacity
  head binds service identity, generation, capacity, active reserved/count/root,
  permanent closed-intent count/root, and quarantined bytes/count/root; the sum
  of active members' reserved-byte fields plus quarantined bytes never exceeds
  capacity. Its proof bundle is at most 16693 bytes: an 8248-byte closed-intent
  proof, at-most-8361-byte active proof, and 8304-byte quarantine proof under the
  operation-specific one/two-proof grammar. Each exact 440-byte update receipt
  binds eleven digest fields—intent, service, member-or-zero, expected and
  committed heads, proof bundle, evidence-or-zero, policy, capability,
  idempotency, and signature—plus committed CAS and durable interval.
  `ReserveCommitted|ReserveRefused|ReserveExpiredUnknown` permanently tombstone
  the intent and respectively add the active member, add no capacity, or
  quarantine the maximum reservation. `AuthorizedRelease` requires either
  preactivation no-effect or post-H3 settlement closure; `QuarantineRelease`
  requires journal reconciliation proving zero late physical effect. Exact-key
  replay returns the stored disposition, and the replicated head fences every
  underlying allocator.

  All physical intents are derived only after the committed global reservation
  receipt. The nine infrastructure members are `PrimaryLane`, `EmergencyLane`,
  both publication journals, `PublicationVerifier`, `CleanupAuthorityStore`,
  `OrdinaryEvidenceStore`, and retention destinations A/B. Lane rows bind lane
  receipts; journal/verifier/cleanup/ordinary rows bind generic committed
  receipts; destination rows bind store `ReserveCommitted` receipts. The
  verifier's reservation includes logical-slot authority. Preactivation abort
  releases every known active hold under pre-CAS no-effect authority; an unknown
  outcome stays quarantined and unavailable. Activated branch release waits for
  the relevant retained artifacts and H3 control closure; cleanup capacity is
  last. Every refusal or expired-unknown disposition remains retained even when
  admission cannot proceed.

  Each lane authenticates capacity with one 218-byte V3 service head containing
  service identity and generation; capacity, active-reserved, active-count and
  active-root fields; a permanent closed-plan count/root; and quarantined-
  reserved, quarantined-count and quarantine-map-root fields. Every successor
  rehashes the complete roots and satisfies checked
  `active_reserved + quarantined_reserved <= service_capacity`. The 121-byte
  lane-plan identity projection binds the global intent, lane/service,
  durability, scratch, two breach slots, 10485760 control bytes, and priority
  segments. A 160-byte active member binds that plan, lane receipt, breach-slot
  reservation root, every component count, and independently recomputed total;
  its service-context 256-level sparse-map proof is exactly 8392 bytes. A
  72-byte quarantine member instead binds the complete active-member digest,
  total reserved bytes, and preactivation-catastrophe evidence in its own
  256-level sparse map; its exact proof is 8304 bytes. Exact 447-byte CAS
  receipts cover
  `IssuanceCommitted`, `IssuanceRefused`, `IssuanceExpiredUnknown`,
  `ActiveRelease`, and authenticated `QuarantineRelease`. One atomic issuance
  CAS always inserts the permanent plan tombstone; Committed also inserts the
  active member, Refused consumes no capacity, and ExpiredUnknown quarantines
  the byte-identical full active-member total. Issuance, quarantine, and release
  all rehash that same member; only authenticated recovery may release a
  quarantine. The exact successor service head, never an unattached side root,
  enters preactivation no-effect and settlement authority. A late response
  cannot resurrect a closed plan, quarantined bytes cannot be allocated, and
  compaction is legal only through a signed archival checkpoint that forever
  preserves membership. The 165-byte service head and asserted 249-byte
  side-root quarantine are replay-only; the distinct 125-byte global-governor
  head remains authoritative.

  The exact 217-byte lane breach-slot reservation binds both phase-ordered slot
  receipts, genesis persistence heads, and maximum artifact bytes before lane
  issuance. Active release additionally requires an exact 356-byte authorization
  binding settlement arm, branch evidence, quiescence, retention-conformance/
  destination quorum or canonical zero for preactivation, breach-slot
  settlement, expected service head, active member, capability, idempotency,
  signature, and interval. Thus no service update receipt can self-authorize
  early release.

  Settlement is a guarded three-arm graph rather than one accidental
  conjunction. Preactivation revokes only capabilities actually issued and
  selects exactly one per-lane `Committed|Refused|ExpiredUnknown` cleanup arm;
  it carries zero scientific-retention digest but retains the no-effect,
  tombstone, quarantine/refusal/release, attempt-set, abort, and settlement
  evidence in the control ledger. Canonical final and catastrophe each use
  their own activation-aware revocation set, quiescence receipt, retention
  requirement set, shard partition, two-store transfer, service/store/global
  releases, successor active set and governor head, and exact 582-byte
  settlement receipt. Catastrophe additionally requires the exact all-four-
  route ledger and 601-byte catastrophe evidence; any Pending or Succeeded
  route forbids catastrophe. Every arm is copied into the separately
  pre-reserved control-plane settlement ledger. The exact order is H0 -> append
  branch evidence E -> H1 -> append the signed 280-byte marker -> H2 -> commit
  the 582-byte global settlement receipt binding H2 and the marker -> append
  that receipt as the post-settlement payload -> H3. The global run reserve is
  not released before H2; physical cleanup capacity is released only after H3,
  whose durable chain retains the receipt without a self-hash cycle. Early
  release, leaked service/global/store
  capacity, stale or forked head, omitted copy, surviving route, incomplete
  shard partition, artifact loss, branch flattening, or a second settlement is
  fail-closed.

  All fourteen leaves expose a once-only lifecycle sufficient to audit
  request-to-observation, drain, finalization, and atomic publication. Terminal
  selection freezes one observer catalog and closes three independently
  reconciled frontiers: descendants, governance transitions, and protected
  access. Every registered mutation reconciles before `execution.drained`;
  dangling `Prepared` operations, unresolved transitions, missing membership,
  head forks, or post-trigger effects are `IntegrityFailed`. A racing
  registration either commits membership and reservation before the closure
  compare-and-swap or has no membership/effect. These catalog/frontier
  memberships are distinct from trigger candidates.

  Every eligible terminal cause contributes a 76-byte canonical candidate with
  calibrated lower time, causal rank, logical sequence, stable source-event
  digest, trigger kind, and terminal subcause. Cancellation ingress first
  persists a nonterminal `cancellation.proposed` event; the proposal neither
  seals nor cancels, and its sequence, interval, and digest provide the
  cancellation candidate without depending on the later request. Its canonical
  projection has the exact `I13_CANCELLATION_PROPOSAL_V1` 28-byte header plus
  NUL, framed scope, epoch/reason/idempotency/sequence/outward-interval/
  registration/capability fields, a 65709-byte maximum, and a domain-separated
  digest. Proposal-key replay returns the immutable original proposal receipt;
  a separate non-authoritative status envelope may progress from Pending to its
  bound winning or losing selection result.

  Candidate bytes are lexicographically sorted only to canonicalize roots. The
  arbiter instead minimizes the decoded typed key: unsigned numeric calibrated
  lower time, unsigned numeric rank, unsigned numeric logical sequence, then raw
  digest bytes. U64LE byte order has no selection authority. Frozen ranks are
  infrastructure failure 0, cancellation 1, timeout 2, budget exhaustion 3,
  and successful work 4; those ranks bijectively determine kind/subcause, so a
  nonidentical candidate with the same key is an integrity error. The 1..=65536
  candidate projection is exactly `16 + 84*count`, capped at 5505040 bytes; its
  cap-plus-one twin is exactly 5505041 bytes. The distinct nonwinner projection
  is `40 + 84*(count-1)`, capped at 5504980 bytes; its independent cap-plus-one
  twin is exactly 5504981 bytes. Each trailing-byte twin refuses against its own
  cap. Two encoders, exact-byte/root KATs, both cap twins, U64LE 255-versus-256
  counterexamples, malformed/set-mutation tests, and an independent typed
  reference arbiter are prerequisites for authority.

  Candidate, witness, and terminal-event registration frontiers use one exact
  fixed-depth authenticated-set algebra. A record key is a domain-separated
  BLAKE3 digest; a 256-level sparse-Merkle path uses digest bits MSB-first, with
  depth-tagged empty and internal-node hashes. The root wraps domain/context,
  count, and tree root. An update proof is old/new count, framed record, and
  exactly 256 sibling digests. Its reachable maxima are 8292 bytes for a
  76-byte candidate, 8283 for a 67-byte witness, and 8568 for the maximum
  352-byte terminal-event record. A valid insert proves old leaf empty, new leaf
  equal to the exact record leaf, and count incremented once; deletion,
  replacement, or occupied-leaf insertion is forbidden. Verification always
  performs 256 bounded hashes, and insertion never sorts or rehashes the set. A
  full-record collision at one record key is `IntegrityFailed`, never a bucket
  alias.
  Closure authenticates every hash-key-ordered leaf, then performs one bounded,
  cancellable external radix sort of full record bytes to stream any required
  flat lexicographic certification root. Hash-key traversal is not mislabeled
  as lexical order. Proof, collision-injection, spill, cancellation, cap, and
  reversed-insertion KATs are mandatory.

  Before any source capability exists, admission freezes exact
  `I13_UNAVOIDABLE_SOURCE_UNIVERSE_V1` bytes: 1..=256 sorted unique 43-byte rows
  bind pre/post-selection phase, source kind, authority identity, and exactly
  one canonical fact. Its at-most-11045-byte projection and root are
  independently rehashed; phase/kind/identity swaps, duplicate authority under
  a new id, or a maximum other than one fail closed. Checked finite-universe and
  internal-source maxima must fit every candidate, witness, terminal-event, and
  mirrored failure-lane reservation before execution.

  Eligible causes enter the attempt/epoch candidate sparse-Merkle frontier in
  the same atomic transaction as their source event. Its context binds both that
  universe root and the exact five-row
  `I13_CANDIDATE_SOURCE_CLASS_MATRIX_V1`: successful-work and
  budget-exhaustion sources are `InternalSchedulable`; independently
  authenticated cancellation ingress, supervisor timeout, and receipt-qualified
  infrastructure failure are `UnavoidableExternal`. No caller, worker, plugin,
  or recovery branch may relabel the class. Its genesis is the empty tree at
  count/generation zero; generation is count, not logical time. Exact-key replay
  inserts nothing and equality at 65536 is valid. An internal cap-plus-one is
  `UnsupportedOperation` before proof allocation/effect and proves admission
  unsound. No internal source may be armed without its ordinary slot; if an
  already authenticated internal fact nevertheless materializes at the cap,
  the supervisor-authorized admission-invariant breach path settles it rather
  than silently dropping or relabeling it.

  A full ordinary frontier opens a separately reserved authenticated
  `BreachFactFrontierV1`; it never collapses a race to one arbitrarily observed
  source. Each exact fact is 153 or 229 bytes and binds phase, source kind,
  closed reason, authority, retained source projection, logical sequence,
  interval, receipt, and the optional exact candidate. Its sparse-Merkle update
  proof is at most 8445 bytes. The frontier progresses
  `Open -> BreachCollecting -> Closed`, contains 1..=257 facts, fences and drains
  the finite source-operation frontier, and reconciles every already registered
  eligible source exactly once before closure. External facts require exact
  unavoidable-universe membership; internal facts require the admission
  receipt, missing-slot contradiction, and supervisor invariant capability.

  Pre-selection breach authority is the closed frontier, not its first opening
  record. The deterministic typed-minimum fact yields a domain-separated
  aggregate digest, while every nonselected fact remains authenticated evidence.
  Its 535-byte prepared selection payload binds that aggregate, closed fact
  root/count, full candidate frontier, expected terminal head, all four mutation
  frontiers, reservation, selection boundary, and the exact maximum sequence of
  `SelectionBatchV1`. Fresh authority appends the 32-byte breach-durability-
  bundle digest, making the exact receipt 567 bytes. The outer mirrored transaction
  is all-or-none, but it has two explicit topological sub-batches: selection and
  drain-start finish first; only then does `PostSelectionWitnessBatchV1` insert
  the receipt-qualified infrastructure witness at a sequence strictly greater
  than the selection-batch maximum. The older singular
  `CandidateCapacityBreachV1` and 501-byte receipt describe only the opening
  fact and legacy draft; they grant no V2 selection authority. A persisted V2
  breach drains and finalizes exactly one non-promoting FailureBundle rather
  than disappearing, wedging, or becoming catastrophe. Cancellation after a
  normal selection remains `AlreadyDraining`; later timeout or infrastructure
  facts become typed witnesses.

  The exact breach artifact makes that retention executable: its fixed header
  binds attempt/epoch/phase, closed fact root/count, and unavoidable-source
  universe; each record-key-ordered entry frames the complete fact and the exact
  independently decoded source projection whose digest/count it carries. The
  257-fact maximum is 269409356 bytes per phase per lane. Each lane issues an
  exact 360-byte persistence receipt; an independent 342-byte mirror certificate
  verifies byte-identical artifacts, and a 185-byte phase/frontier-typed bundle
  binds root, count, artifact, universe, and certificate. Selection or
  finalization cannot use the frontier before this bundle is durable. Add/drop/
  reorder, projection substitution, phase/universe mismatch, asymmetric
  persistence, and cap-plus-one fail closed.

  Normal selection atomically snapshots the current published durability commit,
  candidate root, terminal/index joint head, generations, counts, and authenticated
  traversals. Its causal boundary is
  checked one-plus the maximum snapshotted event logical sequence, not one-plus
  event count. The flat `TriggerCandidateSetRoot` is streamed once from the
  authenticated candidate traversal. The `Arbitrated` receipt binds that first
  boundary and `terminal_selection_transaction_max_logical_sequence`, the
  maximum of every event in the deterministic atomic batch. Selection closes
  the exact candidate root/count and prepares the next joint head; it becomes
  visible only through the mirrored durability-commit CAS. Any
  candidate/noncandidate insert or `Open -> Breached` race forces recomputation.
  Reversing concurrent arrival changes neither root nor winner.

  `I13_TERMINAL_SELECTION_RECEIPT_V2` is the closed, domain-separated union
  `Arbitrated|CapacityBreach`. Its exact arms contain the 592-byte typed-minimum
  arbitrated receipt or the 567-byte aggregate-plus-durability breach receipt,
  making the complete union 635 or 610 bytes. The legacy 501/535-byte breach and
  `TerminalSelectionReceiptV1` are decode/replay-only and cannot arm a fresh V2
  witness, anchor, core, envelope, event, or downstream receipt.
  One fenced `Arbitrated` transaction validates the typed minimum
  and compare-and-swaps observer-catalog, descendant-frontier, transition-
  frontier, and access-frontier heads, roots/counts, and exact reservations. A
  successful transaction commits selection, closes work/child/mutation and
  candidate admission, and inserts `execution.drain.started`. If and only if
  cancellation wins, it also inserts the sole primary
  `cancellation.requested`, seals the scope, freezes all four memberships,
  mints the primary identity/sequence, and inserts `execution.cancelled`. Crash
  or replay exposes all or none. Both alternatives bind the deterministic
  transaction-maximum logical sequence used by every later witness ancestry and
  sequence check. A losing proposal remains in the rejected-set
  root and creates no cancellation effect; non-cancellation winners insert no
  cancellation event.

  Selection opens a separate attempt/epoch/selection-bound witness sparse-
  Merkle frontier whose context includes the selection receipt, its selection-
  batch maximum logical sequence, the unavoidable-source universe, and the
  frozen post-selection source-class matrix. A receipt-qualified witness and
  its source commit atomically, causally descend from the exact selection
  record, have sequence strictly greater than the selection-batch maximum, and
  precede cutoff; integer separation without authenticated ancestry has no
  authority. Equality at 65536 is valid. An internal cap-plus-one is
  `UnsupportedOperation` before effect and proves admission unsound.

  Unavoidable post-selection SLO, infrastructure, or admitted-invariant facts
  also use a `BreachFactFrontierV1`. Closure retains all 1..=257 reconciled
  facts and deterministically selects only a representative ordering key; it
  never discards concurrent nonselected facts. The exact V2 finalization union
  is `Absent=0` or `Present=1 || breach_durability_bundle_digest ||
  U16LE(fact_count)`, hence one or 35 bytes. After closure and mirrored artifact
  certification, one separately reserved 262-byte
  `I13_POST_SELECTION_BREACH_ANCHOR_V1` event causally references the selection
  record and the authenticated breach-fact frontier. Its logical sequence is
  one-plus the maximum of the selection parent and every external source
  dependency. The anchor is included in the terminal snapshot and necessarily
  precedes cutoff; it forces infrastructure bit 4, `IntegrityFailed`, and one
  non-promoting FailureBundle. The older singleton digest union and out-of-band
  no-sequence source model have no V2 authority.

  Ordinary witness records use the exact 25-byte
  `I13_TERMINAL_CONDITION_V1` header plus NUL, closed kind, logical sequence,
  and nonzero receipt digest—67 bytes total. Closure verifies the all-and-only
  root/count and streams the flat witness root once. Final cutoff sequence is
  snapshot maximum plus barrier-level count plus one, never witness count or
  event generation. The V2 breach anchor supplies real causal identity instead
  of assigning a fabricated sequence to an out-of-band digest.

  Terminal event semantics are predecessor-governed, not caller-selected. Every
  canonical event schema is a V2 projection-grammar member selected by the
  semantic matrix. Its exact semantic-conformance object binds inventory,
  parent profile, two-evaluator set, and positive/negative KAT roots; the V2
  terminal authority binds the independently reproduced all-and-only catalog.
  The event digest hashes the activated terminal context, semantic-authority and
  causal-schema roots, parent set, exact kind, schema and conformance identities,
  canonical event bytes, and idempotency key. It never hashes the superseded
  cyclic admission digest.

  Parent extractors admit only typed explicit event references, consumed
  authority receipts, authenticated breach-frontier references, `None` where
  the governed schema proves no dependency, or canonical nonempty `AllOf`.
  Consuming a receipt resolves its 133-byte binding and inserts the producer's
  terminal record key into the actual causal parent set; sequence comparison is
  not a substitute for authenticated ancestry. Preselection breach context is
  exactly candidate-frontier/matrix plus a zero selection digest; postselection
  context is exactly witness- or terminal-frontier/matrix plus the nonzero V2
  selection. The global terminal-frontier context is a fixed 112-byte activated
  projection and contains neither mutable receipt-index state nor a selection
  digest, avoiding a selection self-cycle.

  Receipt lookup is a fixed-context 256-level sparse-Merkle map. Single
  insertions have at-most-8349-byte proofs; deterministic batches contain at
  most 1024 sequential proofs, 8557690 bytes, and 262144 hashes. A corrected
  127-byte joint head binds terminal and receipt-index roots, generations, and
  counts, with exact empty genesis and generation=count invariants. Replacement,
  second producer, duplicate batch key, split root/index CAS, all-`None`
  laundering, generic-digest role substitution, or omitted external dependency
  is `IntegrityFailed`. V1 causal/schema/event domains are replay-only.

  The lane-agnostic terminal authority is an order-independent sparse-Merkle
  set. A terminal record contains causal logical sequence, parent-set count/root,
  framed event kind, and event digest, at most 352 bytes. Parent identities are
  0..=64 sorted unique record keys in an at-most-2056-byte projection. Their
  canonical multiproof is reconstructed over the union trie of all target paths:
  shared nodes hash once, canonical empty siblings are omitted, and exactly one
  entry represents each required nonempty target-free sibling subtree. The
  proof has at most 16384 entries, 1081352 bytes, and 16384 node hashes. Missing,
  redundant, default-empty, duplicate, target-containing, or inconsistent
  entries fail closed, as do self/future/omitted/extra parents and cycles. Event
  generation remains set count, not logical time; incomparable events may share
  a sequence and full record bytes break canonical ties. Primary and emergency
  stores converge on this logical root while excluding lane-operational roots
  from canonical identity.

  V2 attempts retain all 177407 ordinary nonbarrier slots and reserve event
  177408 for the breach anchor. A maximum 177408-event snapshot produces 2816
  barriers in two levels and one final record, so the exact absolute cap is
  180225. Admission statically proves every frozen extractor, source universe,
  parent proof, barrier, ordinary route, breach route, and mirrored reservation
  fits before issuing a capability; equality passes and cap-plus-one refuses.

  Finalization causally dominates an arbitrarily wide valid snapshot through a
  deterministic 64-ary barrier DAG. Sorted snapshot keys are chunked by 64;
  each chunk produces one exact `I13_TERMINAL_CAUSAL_BARRIER_V1` record binding
  attempt/epoch, snapshot root/generation/count, derived level, chunk, parent
  root/count, and the complete sorted parent-key list. Its size is exactly
  `151 + 32*parent_count`, at most 2199 bytes. The event-kind extractor selects
  that typed list and independently rehashes it; root/count alone cannot
  self-assert missing parents. Its sole idempotency key is domain-derived from
  those exact bytes; no caller or lane key exists. Generated
  keys are rechunked until at most 64 final parents remain. The barrier-plan
  digest binds snapshot root/generation/count, every level/chunk record, level/
  event counts, and the exact final-parent set. Those final parent keys and root
  are embedded byte-for-byte in the applicable finalization envelope so the
  finalized event's extractor can reproduce its all-and-only parents. Zero
  through 64 snapshot events need no barrier; 65 needs one level/two events. A
  legacy V1 ordinary-only
  snapshot of 177407 events reached 180224 after 2816 barriers and finalization.
  For a V2 attempt, the reserved breach anchor may extend the snapshot to
  177408; the same two levels still contain 2816 barriers, so finalization reaches
  the exact 180225 V2 cap. Cutoff sequence is snapshot maximum plus barrier
  levels plus one. Witness/nonbarrier admission closes first, then all barriers
  and finalization commit in one all-or-none deterministic topological batch.
  Caller grouping/order, missing/extra barriers, or a cap-plus-one plan fails
  closed.

  Finalization remains deliberately acyclic. After `execution.drained`, checked
  arithmetic derives drain latency, inclusive Finalize and Total deadlines, and
  `AlreadyExceededByDrain` when total subtraction underflows. Causal subtraction
  underflow and deadline-add overflow remain integrity errors. Tickets use the
  exact `I13_TERMINAL_DEADLINE_TICKET_V1` header and are 202 bytes for an
  absolute deadline or 210 bytes for `AlreadyExceededByDrain`. Both arm
  atomically after drain, match attempt/epoch/profile/capability, and are never
  caller-selected. Exceeded requires an authenticated pre-cutoff interval lower
  strictly beyond the inclusive deadline; a straddling interval remains pending.

  `I13_FINALIZATION_CORE_V2` uses the exact V1 post-header field order but a V2
  header/domain, V2 barrier-plan digest, and 1-or-35-byte V2 breach union. The
  core freezes attempt/epoch, authenticated lane-epoch
  head/value, idempotency, selection/drain,
  expected terminal and witness roots/generations/counts, tickets, two-entry
  resolution template, exact barrier-plan digest/level/event counts, proposed
  cutoff, witness root/count, the framed 1-or-35-byte post-selection breach union,
  condition/disposition, evidence digest, exact
  record size, both service reservations, and capability epoch. Pending-OnTime
  is the one-byte union arm with no
  payload; only the eventual resolution is fixed at 102 bytes. The core excludes
  the grant, OnTime receipt, resolution root, finalized digest, and downstream
  receipts.

  Before serialization, the selected sink issues one authenticated 436-byte
  `CommitBeforeGrantV1`. It binds core, tickets, expected roots/counts,
  lane/lane-epoch (exactly Primary/0 or Emergency/1),
  reservation/service identity, ready upper, explicit prepare-and-submit,
  queue, durable-service, and clock bounds, checked certified upper, expiration,
  one-shot permit, and idempotency. Prepare-and-submit covers validation,
  resolution/envelope construction, canonical serialization, and enqueue; it
  cannot be hidden in queue or service. Pending OnTime is legal only when the
  pre-existing certified upper meets the inclusive deadline. Otherwise no grant
  is issued and a definitely expired ticket first receives its durable Exceeded
  receipt/witness.

  `OnTimeResolutionReceiptV1` has an exact 142-byte encoding and domain-
  separated digest over ticket, grant, cutoff sequence, decision interval, and
  idempotency. Exactly two 102-byte deadline resolutions exist, Finalize then
  Total; their set root has an exact 292-byte projection. Exceeded is after
  arming and before cutoff; OnTime uses the cutoff/finalization sequence.
  Unknown/swapped/missing/duplicate/cross-ticket records, cycles, bad intervals,
  bad sequence, root/count mismatch, cap-plus-one, or split commit fail closed.

  `I13_FINALIZATION_ENVELOPE_V2` likewise has a distinct V2 header/domain and
  exactly binds the core/grant identities,
  lane epoch, selection/drain, expected roots/generations/counts, immutable
  barrier plan/counts, the exact 0..=64 final-parent keys/root, cutoff,
  all-and-only ordinary witness and resolution roots/counts/records, framed
  `Absent|Present` capacity-breach union, condition/disposition/integrity,
  evidence, idempotency, and capability epoch. Its fixed overhead excluding
  parent keys and union payload is 763 bytes; the union payload is exactly one
  or 35 bytes, and the checked evidence bound is
  `1047813 - union_bytes - 32*final_parent_count`. Thus the complete envelope is
  764..=1048576 bytes without weakening the one-mebibyte record cap; a
  zero-parent Present V2 envelope admits at most 1047778 evidence bytes.
  Its digest is domain-separated. It contains no successor root, durability
  observation, or downstream receipt, eliminating a timing/content cycle.
  `I13_EMERGENCY_TERMINAL_SLOT_V2` retains the 368-byte overhead and 1048944-byte
  maximum but embeds only the V2 envelope/digest under its own V2 domain. Stable
  grant, OnTime, finalization-commit, and journal receipt layouts may carry V2
  digests only after independently decoding the V2 objects. V1 core/envelope/
  slot bytes are historical replay artifacts, never fresh authority.

  The conditional sink consumes the permit and either makes the exact envelope,
  terminal-set insertion, witness closure, cutoff, and `execution.finalized`
  durable by the certified upper, or has zero canonical finalization effect.
  Successor terminal generation is checked predecessor plus barrier-event count
  plus one; successor root is recomputed by inserting every barrier record and
  the finalized record into the sparse-Merkle set, not by an arrival-ordered
  predecessor hash. A separate exact 388-byte
  `FinalizationCommitReceiptV1` later attests lane/lane-epoch,
  core/grant/envelope, predecessor/
  successor roots, interval, reservation, permit, idempotency, and authority.
  Missing recovery blocks adjudication; late/forged evidence never rewrites the
  already durable axes.

  Downstream receipts live in a separate pre-reserved authenticated ledger. A
  precommit marker binds grant/permit/envelope/expected root. Crash recovery
  commits one conservative interval whose upper is the recovery-observation
  upper; it never backdates or narrows. Exact replay returns byte-identical
  receipt bytes. Finalized latency ends at durable envelope plus terminal-root
  insertion, not at the receipt's own persistence; receipt return has a separate
  bounded SLO.

  The selected trigger supplies one base bit; the all-and-only ordinary witness
  root adds receipt-derived resource, SLO, and infrastructure bits, while a
  Present capacity-breach union necessarily adds infrastructure bit 4 and
  `IntegrityFailed`. Completed/Cancelled remain base-only. Tags are exactly
  `Completed=0`, `Cancelled=1`, `BudgetExhausted=2`, `TimedOut=3`, and
  `InfrastructureFailed=4`; the derived tag is the maximum set-bit index in the
  five-bit condition byte. The pure decoder classifies all 31 syntactically
  bounded nonzero inputs, the zero mask, and all 256 possible tag bytes.
  Lifecycle acceptance is stricter: exactly one earned base bit plus only earned
  post-condition bits. Final disposition is therefore independent of
  registration/serialization order. No mutable setter exists; later adjudication
  maps admissible SLO overrun to claim `Unknown`, and evidence integrity remains
  an orthogonal axis.

  Each primary/emergency lane retains five cumulative precommit artifacts plus
  a post-publication metadata chain. Structural chunks have contiguous
  zero-based ordinals, exact predecessor/successor roots, 1..=1024 globally
  ordered records, and at most 63827439 bytes at 180225 events. A separately
  capped 67108864-byte terminal-payload artifact retains every record key,
  idempotency digest, and canonical event bytes required to rehash and decode
  the typed event. A key-ordered receipt-index artifact retains every 133-byte
  binding and is at most 31179051 bytes. These three total at most 162115354
  bytes. Structural chunks have no hidden outer frame; lane identity belongs to
  their surrounding artifact and bundle authority. Each event is capped at
  1048576 canonical bytes while the complete payload remains capped at
  67108864, and every declared count/total equals its checked row sum. The
  receipt-index artifact retains the exact 112-byte context preimage rather
  than only its digest. An opaque logical root without source payload and
  enumerable index bytes has no recovery authority.

  The V3 operation-payload matrix has exactly seven closed kinds:
  `OrdinaryEventBatch`, `CandidateSelection`, `CapacityBreachSelection`,
  `BarrierBatch`, `Finalization`, `Failover`, and `Recovery`. Its exact matrix is
  304 bytes. Each typed operation authority is exactly `84 + payload_byte_count`
  and at most 4194304 bytes, binding kind, its matrix-selected schema digest,
  exact count, and framed canonical payload. The cumulative operation-payload
  artifact has a 98-byte fixed projection followed by generation-contiguous
  rows that bind generation,
  digest, and the exact typed operation authority; it is capped at 4096 rows
  and 67108864 bytes. `OrdinaryEventBatch` cannot alias `BarrierBatch`, and the
  byte-identical selected tagged arm must feed event preparation, the cumulative
  artifact, and commit for the six top-level arms `OrdinaryEventBatch`,
  `CandidateSelection`, `CapacityBreachSelection`, `Finalization`, `Failover`,
  and `Recovery`; `BarrierBatch` is the exact matrix-validated nested
  subprojection of the single `Finalization` operation when derived barriers are
  nonempty; it does not consume another durability generation, append-artifact
  row, or Closing slot. The enclosing Finalization arm binds its complete bytes
  and digest, and a caller cannot submit BarrierBatch as a standalone append.
  The top-level tagged sum consequently has six arms while the schema matrix has
  seven rows. This preserves the nested seven-kind grammar without making the
  4096-append phase budget internally inconsistent. The cumulative append-delta
  artifact has a 102-byte fixed projection plus one eight-byte frame per exact
  ordered delta and is at most 11849894 bytes. The post-publication append-
  metadata artifact has a
  96-byte fixed projection plus exact generation-ordered 4974-byte rows, each
  comprising generation, published-authority digest, and the framed exact
  4926-byte request/receipt/ack/prefix/commit/head/journal/certificate bundle;
  it is at most 20373600 bytes. A metadata successor is built only after the
  publication certificate, must exist before the next append or settlement,
  and chains from the complete predecessor metadata artifact. It is never
  hashed backward into the request, commit, or certificate it contains.

  The exact 257-byte `LaneContentBundleProjectionV3` binds terminal context,
  attempt, epoch, lane, generation, and five root/byte-count pairs in fixed
  structural, terminal-payload, receipt-index, operation-payload, append-delta
  order. Its domain-separated `LaneContentBundleRootV3` is committed before
  publication, which produces the independently chained metadata
  root afterward. Thus final per-lane source-retaining durable use is at most
  `162115354 + 67108864 + 11849894 + 20373600 = 261447712` bytes, leaving
  exactly 6987744 bytes of the 268435456-byte durable reservation as
  non-borrowable, non-authoritative safety slack. No index, metadata, or hidden
  artifact may consume it. Historical snapshots, wrappers, and
  proofs are deterministic reconstructible views of those retained source
  artifacts plus append boundaries, not uncharged duplicate blobs; inability
  to reconstruct any historical preimage is `IntegrityFailed`.

  The exact V3 published durability head is 218 bytes and carries phase
  `Open|Closing|Finalized` plus `remaining_closing_slots`. Open admits at most
  4090 nonclosing `OrdinaryEventBatch` appends. Exactly one
  `CandidateSelection|CapacityBreachSelection` append atomically changes Open
  to Closing. Closing reserves six potential appends: that selection, at most
  four exact replay/failover route outcomes, and exactly one Finalization.
  Inapplicable route slots are retired only by the all-four-route ledger when
  Finalization changes Closing to Finalized and sets remaining slots to zero.
  Selection's published successor has remaining count five. Each distinct route
  ordinal may consume at most one physical Failover/Recovery append and
  decrement it once; Exhausted/Inapplicable with zero canonical append burns
  the ordinal in the route ledger without changing the published counter. A
  direct Finalization consumes the dedicated last slot and retires unused route
  reserves. If a successful route publishes the byte-identical Finalization,
  that one physical append consumes both its route ordinal and the dedicated
  final slot; a second finalization generation is forbidden.
  Finalized rejects every append permanently; capability revocation can only
  remove remaining authority and can never reopen it. At
  every Open append the governor preserves, including row/frame overhead, at
  least 25166112 operation-payload bytes, 13117392 terminal-payload bytes,
  2810954 structural bytes, 451054 append-delta bytes, 1062912 receipt-index
  bytes, and 29844 metadata bytes for the closing tail. The finalization share
  independently preserves at least 4194352 operation-payload bytes, 7464912
  terminal-payload bytes, 997689 structural bytes, 122989 append-delta bytes,
  at most 1024 final receipt-index entries, and one 4974-byte metadata entry.
  Equality passes; exhaustion by one byte, row, record, or slot refuses the
  nonclosing append before it consumes authority.

  Fresh publication is V3 end to end. The journal receipt and independent
  publication certificate use exact 342-byte and 410-byte V3 projections and
  bind V3 head/authority digests. The published-authority digest hashes the
  exact 218-byte head digest together with primary and emergency journal-
  receipt digests; the independently valid certificate remains a separate
  visibility prerequisite. Append metadata stores only those V3 preimages. A
  V2 head, journal receipt, certificate, authority digest, or metadata root is
  replay-only and cannot enter a fresh V3 append.

  Scratch lifecycle is separately authenticated. The exact 163-byte append-
  phase head carries attempt, epoch, lane, generation, state
  `Idle=0|Acquired=1|Prepared=2|Quarantined=3`, lease and idempotency identities,
  reserved/live scratch counts, and plan-or-prepared root. Idle requires every
  optional field and count to be zero. An exact 376-byte pre-CAS phase lease
  binds attempt/epoch/lane, primary and optional secondary append-slot ordinals,
  top-level operation, context, expected phase-authority head, the exact global
  logical-slot allocation receipt, expected published authority, operation and
  component-plan digests, maximum scratch, capability, idempotency, deadline,
  and signer set. Its exact 153-byte closed component plan binds the same
  identities and six component maxima in structural, payload, receipt-index,
  operation, delta, metadata order; prepublication scratch is the checked first-
  five sum and at most 134217728.

  Logical slots are allocated once globally, not independently by the two lane
  services. The exact 258-byte allocation intent binds attempt, epoch, lane-
  neutral context, operation kind, primary/secondary ordinals, expected logical-
  slot authority head, operation payload, idempotency, capability, and signer.
  The exact 164-byte verifier-owned logical-slot authority head binds its service,
  attempt, epoch, generation, slot-ledger root, burned count, and latest
  allocation receipt. Its exact 335-byte CAS receipt binds the intent, expected
  and committed authority heads, old/new generations, old/new slot-ledger roots,
  capability, idempotency, signature, CAS disposition, and durable interval.
  This replicated verifier-owned CAS commits before either lane may `Acquire`;
  both lane leases must byte-equal the same allocation receipt. Failure of one
  lane after the allocation leaves the global slot burned and forces
  deterministic paired abort/quarantine reconciliation. It never authorizes
  divergent reuse, and publication still requires both lane phase authorities
  to bind that allocation.

  The exact 633-byte slot ledger contains attempt, epoch, lane-neutral context,
  generation, burned count, next ordinary ordinal, selection state, route
  bitmap, final state, a complete 4096-slot bitmap, and predecessor root.
  Ordinary uses 1..=4090, selection 4091, routes 1..=4 use 4092..=4095, and
  direct finalization uses 4096; a route-published finalization names its route
  slot primary and 4096 secondary. No other pair is valid. The exact 348-byte
  per-lane phase-authority head binds current inner head, the shared logical-slot
  authority head, lease/plan/transition/head artifacts and counts, latest
  publication, and latest metadata root. It is the sole linearizable CAS for
  that lane's phase state; the verifier-owned logical-slot authority is the sole
  CAS for global slot allocation. Thus neither lane history nor burned-slot
  state can lag its governing authority. Historical logical-slot heads, ledgers,
  intents, and allocation receipts must reconstruct the global allocation chain,
  while the complete lease/plan/transition/head artifact prefixes plus published
  metadata reconstruct each outer lane-authority chain. Every reconstructed
  root must match the next governing receipt or lease and the final closure in
  bounded streaming passes, never quadratic prefix replay. `Acquire`
  permanently burns its classed ordinal
  even if the phase later aborts or quarantines; exact-key replay adds nothing.
  `Acquire` changes Idle to
  Acquired and binds the exact plan and maximum component bytes. `SealPrepared` binds the
  rehashed components and exact live bytes. `CommitRelease` requires the
  published metadata successor and returns Prepared to Idle. `AbortReclaim`
  proves zero canonical effect before returning Acquired or Prepared to Idle;
  only `Quarantine` and authenticated `RecoveryRelease` traverse Quarantined.
  Exact 317-byte transition receipts bind the transition tag, expected and
  committed heads, old/new generations, expected publication authority,
  plan-or-operation digest, requested-or-live scratch count, capability,
  idempotency, signature set, CAS disposition, and durable interval.

  `AppendPhaseReceiptArtifactBytesV3` has a 102-byte fixed projection and at most
  16384 framed 317-byte receipts, each consuming 325 bytes; it is exactly
  `102 + 325*receipt_count` and at most 5324902 bytes. A successful append uses
  exactly Acquire, SealPrepared, and CommitRelease; any burned slot uses at most
  four non-replay transitions and replay adds no row, so retries cannot overflow
  the 16384-receipt cap. `AppendPhaseHeadArtifactBytesV3` binds attempt, epoch,
  lane, context, head count, total head bytes, the exact 163-byte head width,
  and `completeness=Closed`; it retains the exact genesis and every committed
  successor in transition order. Its count is exactly `receipt_count + 1`, at
  most 16385; with a 102-byte fixed projection and 171 bytes per framed head,
  its maximum is exactly
  `102 + 171*16385 = 2801937`. Every receipt's expected/committed digest must
  rehash adjacent retained heads. `AppendPhaseLeaseArtifactBytesV3` retains at
  most 4096 framed 376-byte leases and is at most 1572964 bytes;
  `AppendComponentPlanArtifactBytesV3` retains the corresponding framed 153-byte
  plans and is at most 659555 bytes. All four histories chain complete
  predecessor bytes. `LaneControlClosureArtifactBytesV3`, capped at 126402
  bytes, binds the final phase-authority head, shared logical-slot authority head
  and ledger, all four artifact roots/counts, active member, service updates,
  breach-slot settlement, lane release, and framed typed preimages. The exact
  five-artifact control partition totals 10485760 bytes. An exact 178-byte two-
  lane phase-authority set binds both final heads and their byte-identical shared
  logical-slot-authority-head digest. Each
  lane has one 134217728-byte scratch reservation and at most one in-flight
  phase; both lanes execute the same logical append plan. Crash/replay resumes
  the exact key or quarantines and authentically reclaims bounded bytes before
  any later acquire. A stale loser, abandoned Prepared bytes, post-Finalized
  request, second simultaneous lease, missing history receipt, or control-cap
  excess has zero publication authority.

  Prepared terminal/index state remains invisible. Both lanes bind the same
  predecessor published authority and exact logical successor; independent
  acknowledgements and prefix verification precede an all-or-none durability
  commit. Two independent journal receipts plus the verifier certificate then
  publish the sole consumable head. Selection, barriers, finalization, failover,
  and recovery consume only this head-plus-certificate authority, never a
  terminal-only, one-journal, one-lane, Prepared, or metadata-incomplete state.
  Counts and reservations never grow after activation, and release occurs only
  through the exact branch settlement CAS. The 216-byte head, six-kind matrix,
  4924-byte metadata bundle, 20168704-byte metadata cap, 182284058-byte final
  total, structure-only lane budget, and 256-MiB attempt governor are historical
  replay values with no fresh V3 authority.

  The emergency writer is a complete mirrored typed lifecycle lane, so primary
  failure before selection or drain remains recoverable. `LaneEpochHeadV1` has
  an attempt/reservation-bound Primary/0 genesis. The only pre-settlement
  transition is Primary/0 to Emergency/1 through exact 327-byte V2 intent bytes.
  They extend the 295-byte legacy projection with a digest of an exact 310-byte
  `I13_LANE_COMMON_PREFIX_CERTIFICATE_V1`. Two independently durable append acks
  use an exact 240-byte `I13_LANE_APPEND_ACK_V1` encoding and domain-separated
  digest. Each acknowledgement binds lane, reservation, the lane-agnostic head/
  generation/count, its rehashed artifact/count, append idempotency key, durable
  interval, and independently serviced operation receipt. Both acknowledgements
  must decode to the same contiguous ordered terminal records and logical head,
  while retaining their distinct lane, artifact, reservation, and receipt
  identities. Certificate interval is the deterministic hull of the two ack
  intervals; its idempotency key and digest are derived from the exact ack/head
  projection, never caller-selected.

  The operation and verifier evidence are themselves exact. Each 370-byte
  durable-append receipt binds expected predecessor and committed successor
  joint-head identities/generation/count, the exact lane reservation,
  idempotency, lane-bundle root/count, V2 append request, independent service/
  capability, committed CAS, and interval. The 404-byte prefix receipt binds
  both acknowledgements/bundles, common logical/index content, independently
  owned verifier implementation/toolchain, verdict, interval, signatures, and
  idempotency. The certificate's historical
  `independent_prefix_verifier_capability_digest` field byte-equals this receipt
  digest; neither is an opaque assertion.

  A V1 intent without that certificate has no failover authority. V2 uses the
  distinct `org.frankensim.i13.lane-epoch-transition.v2` digest domain and
  derives the successor lane head from that V2 digest; the legacy V1 digest is
  not an alias. The exact 230-byte transition receipt binds those V2 transition
  and successor-head values. Uncertain, missing, gapped, overlapping, or
  divergent prefixes remain `TerminalPersistencePending`; the implementation
  may not guess, skip, or replay. Emergency finalization equality-fences the
  already committed Emergency/1 lane head while compare-and-swapping journal and
  logical terminal heads; it does not invent an undefined post-Emergency lane
  successor. Exact replay returns the receipt; stale/ABA/conflicting transitions
  fail. A recovered primary never regains canonical append authority for that
  attempt; copyback is non-authoritative.

  Replay and catastrophe authority are closed over an admission-bound route
  universe, not an open-ended recovery suggestion. The exact 742-byte failure
  matrix has eleven ordered stages:
  `AuthorityAndFence`, `PrimaryAppend`, `MirrorAndPrefixVerification`,
  `DurabilityCommit`, `DualJournalPublication`,
  `IndependentPublicationVerification`, `FailoverTransition`,
  `EmergencyPrefix`, `EmergencyFinalization`,
  `DownstreamReceiptRecovery`, and `ArtifactReconstruction`. Failure kinds are
  numbered exactly `CapabilityRevoked=1`, `CapabilityExpired=2`,
  `ReservationUnavailable=3`, `ServiceTerminalRefusal=4`,
  `ServiceUnavailablePastDeadline=5`, `DurableWriteNotCommitted=6`,
  `ConflictingCasCommit=7`, `RequiredArtifactAbsent=8`,
  `RequiredArtifactInvalid=9`, `PrefixNotEquivalent=10`,
  `RequiredReceiptAbsent=11`, `RequiredReceiptInvalid=12`,
  `ReconstructionImpossible=13`, `PublicationQuorumUnavailable=14`, and
  `IndependentVerifierRejected=15`. In stage order, the exact allowed code sets
  are `{1,2,3}`, `{1,2,3,4,5,6,7,8,9}`,
  `{1,2,3,4,5,6,7,8,9,10}`, `{1,2,3,4,5,6,7,8,9,10,11,12}`,
  `{1,2,3,4,5,6,7,11,12,14}`, `{1,2,8,9,11,12,15}`,
  `{1,2,3,4,5,7,11,12}`, `{1,2,3,4,5,6,7,8,9,10}`,
  `{1,2,3,4,5,6,7,8,9,10,11,12}`, `{1,2,3,4,5,8,9,11,12}`, and
  `{1,2,3,8,9,11,12,13,15}`. Their cardinalities sum to exactly 99; a
  stage/failure pair outside that matrix has no exhaustion authority.

  The exact 656-byte route registry binds that matrix and exactly four ordinal-
  zero routes. `PrimaryExactReplay` traverses stages `[1,2,3,4,5,6]`;
  `MirroredPublicationReplay` traverses `[1,5,6]`; `EmergencyFailover`
  traverses `[1,7,8,9,4,5,6]`; and `DurablePrefixRecovery` traverses
  `[1,10,11,4,5,6]`. Both publication journals plus the independent verifier
  certificate are one indivisible publication stage. Four exact 399-byte route
  specs bind attempt, epoch, route kind/ordinal, terminal context, activation,
  starting published authority, target operation, frozen precondition and
  program/success schemas, capability set, clock/deadline, and idempotency.
  Their exact all-and-only universe is 1814 bytes and its root, the failure-
  matrix root, and the registry root are all admission-bound.

  An exact 521-byte failure receipt binds universe/spec, route, scope
  `Global|Primary|Emergency`, allowed stage/failure pair, effect
  `NoCanonicalEffect|DurablePrefixOnly`, affected artifact, request, evidence
  schema and evidence preimage, pre/post state, service, capability, clock,
  idempotency, independent verifier, signature, and durable interval. An exact
  296-byte inapplicability proof can establish only that the frozen precondition
  is permanently false; transient service failure is never Inapplicable. Each
  exact 492-byte outcome receipt is one of `Exhausted|Inapplicable|Succeeded`
  with effect `NoCanonicalEffect|DurablePrefixOnly|CanonicalPublished`.
  Succeeded binds exactly one resulting authority and no failure/inapplicability
  digest; Exhausted binds one valid failure receipt; Inapplicable binds one
  valid proof; absence is Pending. The 385-byte replay ledger contains exactly
  one framed row for each registered route. Each V3 failure inventory has scope
  `Global|Primary|Emergency` and exact size `44 + 74*failure_count`; a fetched
  receipt must byte-equal that scope and occurs in exactly one inventory.
  Catastrophe evidence is exactly 601
  bytes and additionally binds the global failure inventory, route registry,
  and admissible universe. It accepts only all-four Exhausted/Inapplicable
  outcomes, each failure in exactly one scope inventory, all capabilities
  quiesced, and no success. Any DurablePrefixOnly result forces
  DurablePrefixRecovery and retention of every surviving prefix byte. A fifth
  route, missing row, caller route, one-journal success, malformed verifier,
  duplicate-scoped failure, invalid inapplicability, or success relabeled as
  catastrophe remains Pending or `IntegrityFailed`.

  Canonical terminal evidence is also V3-complete. The exact 567-byte
  `CanonicalFinalizationEvidenceBytesV3` upgrades every identity, authority,
  commit, envelope, record, and quiescence digest in the 375-byte V2-form body
  and additionally binds operation-payload, append-delta, append-metadata, and
  two-lane phase-authority-set roots plus conditional replay-route-ledger and
  breach-evidence-bundle roots. Each nonzero conditional root must fetch exact
  bytes; each zero root requires its closed predicate to prove inapplicability.

  Retention is typed, aggregate-complete, shard-aware, capacity-backed, and
  branch-specific. An arm requirement set has a 56-byte fixed projection plus
  one framed 49-byte row per role, so its exact size is
  `56 + 57*requirement_count`. A row binds its role, the canonical source-
  artifact-set root, and exact required artifact and byte counts. The thirteen
  closed roles are `TerminalStructural`, `TerminalPayload`, `ReceiptIndex`,
  `AppendMetadataAndDelta`, `BreachEvidence`, `FailureRecoveryEvidence`,
  `FinalEnvelope`, `PublicationEvidence`, `QuiescenceEvidence`,
  `SettlementEvidence`, `OperationPayload`, `OrdinaryObservability`, and
  `CleanupAuthority`. Each retained member is the exact 43-byte
  lane/role/digest/byte-count projection. A V3 retained set has kind
  `PreEvidenceSurvivors|TransferSet|RequirementSource|SourceShard|
  ControlEvidence`, a 62-byte fixed projection plus 51 bytes per framed lexical
  member, and globally unique artifact digests. Every requirement source root
  resolves to an exact
  RequirementSource set whose members all have that row's role and whose
  count/bytes equal the row. Publication evidence
  includes every head, journal receipt, certificate, and service receipt;
  quiescence evidence includes issued-capability catalog, state heads, every
  revocation receipt, revoked set, and quiescence receipt; failure/recovery
  evidence includes all three failure inventories, every failure,
  inapplicability and outcome preimage, the route matrix/registry/universe/
  ledger, and every affected or surviving-prefix artifact.

  The admission-frozen retention policy binds the exact 1880-byte, thirteen-row
  `RetentionRoleInventoryMatrixBytesV3`. Each row fixes its producer-set root,
  lane rule (`LaneNeutralUnique|LaneSpecificAll|SurvivingLaneSpecific|
  ControlOnly`), conditional predicate, count and byte formulas, and overlap
  rule (`DigestDisjoint|OwnedAggregatePreimages`). Roles 1,2,3,4,7,11,12 are
  lane-neutral named final roots; role 5 owns the breach bundle and all of its
  logical artifacts, lane receipts, fact-durability bundles, and certificates;
  role 6 owns all route matrices, registries, universes, ledgers, inventories,
  receipts, and affected artifacts; roles 8,9,10 own publication, quiescence,
  and the singleton terminal-evidence preimage respectively. Role 13 is
  `ControlOnly`: it is all-and-only the separately retained control-ledger
  inventory and is forbidden from scientific RequirementSource or TransferSet
  bytes. Shared subpreimages have one owning aggregate role, while top-level
  source-set digests are otherwise disjoint. Manifest generation evaluates the
  frozen producer/count/byte formulas; callers cannot select a favorable
  inventory.

  Every committed artifact producer atomically emits an exact 328-byte
  retention-producer insertion receipt binding attempt, epoch, branch arm,
  role, producer schema, artifact digest/byte count, lane rule, source identity,
  expected/committed role-registry heads, old/new generations, idempotency,
  signature, committed CAS, and durable interval. Quiescence freezes every role
  registry's root, count, and byte sum. RequirementSource sets are generated
  all-and-only from those frozen receipts; scientific conformance bijects every
  required receipt to one TransferSet member and every member to one receipt.
  `ControlOnly` receipts instead generate set-kind `ControlEvidence=5`, which
  appears all-and-only in pre-settlement E and never in the scientific
  TransferSet.

  Canonical final requires roles 1,2,3,4,7,8,9,10,11,12; it additionally
  requires role 5 if either phase's breach union is Present and role 6 if any
  route outcome is Exhausted or Succeeded rather than permanently Inapplicable,
  or any Failover/Recovery lease was acquired.
  Thus every actual replay, failover, or recovery attempt is retained, while an
  all-Inapplicable no-recovery execution does not fabricate failure evidence.
  Its SettlementEvidence source set is the
  exact canonical singleton containing the finalization-evidence artifact,
  rather than type-confusing an artifact digest with a set root. Catastrophe
  pre-evidence requires roles 1,2,3,4,6,8,11,12 and role 5 exactly when at
  least one breach preimage survives. A logical Present union whose two
  physical copies were both lost is proved through role-6 failure evidence; it
  cannot require the unavailable artifact and thereby make catastrophe
  unrepresentable. The pre-evidence set excludes not-yet-existing quiescence/
  catastrophe evidence; after the exact 601-byte
  catastrophe evidence exists, the transfer adds role 9 and a role-10 singleton
  source set containing that evidence. No role may name the eventual transfer
  or settlement receipt and create a hash cycle.

  An exact 305-byte V3 retention-conformance receipt binds attempt, epoch, arm,
  requirement and TransferSet roots, policy, roster, requirement/artifact
  counts, evaluator, idempotency, signature, `AllAndOnly` verdict, and interval.
  It proves every required source member occurs exactly once in the TransferSet,
  every TransferSet member is authorized by exactly one requirement, and every
  aggregate preimage is available. A favorable total without this constructive
  all-and-only proof has no transfer or settlement authority.

  Breach identity is logical and lane-neutral even though persistence is
  independently mirrored. Each Present phase therefore contributes one unique
  logical breach-artifact digest, two lane-specific persistence receipts, one
  fact-durability bundle, and one mirror certificate; Absent contributes none.
  The exact V3 breach-evidence bundle is `47 + 177*phase_count` bytes for phase
  count 0..=2 and binds each phase, logical digest, byte count, both lane
  receipts, fact-durability-bundle digest, and mirror certificate. Every bound
  durability-bundle preimage is retained. Final `BreachEvidence` consequently
  contains exactly 0, 1, or 2 unique logical artifact digests, 0, 2, or 4 lane-
  specific persistence receipts, 0, 1, or 2 fact-durability bundles, and 0, 1,
  or 2 mirror certificates. It never attempts to insert duplicate equal
  digests for the two copies into a unique-digest retained-artifact set. Each of
  the four physical lane/phase slots is nevertheless capacity-reserved and its
  release proves either exactly one consumption by the bound logical artifact
  or its canonical generation-zero unconsumed head. Catastrophe retains the
  surviving copies and receipts named by its failure/survivor inventory.

  Each exact 138-byte lane breach-slot settlement member binds phase, slot
  reservation receipt, current persistence head, disposition
  `UnusedGenesis|Consumed`, logical artifact or zero, lane persistence receipt
  or zero, and byte count. The exact 378-byte lane settlement contains both
  phase-ordered rows and the active-capacity-member digest. Unused requires the
  frozen genesis head and zero artifact/receipt/bytes; Consumed requires the
  exact successor and nonzero bound values. The exact 467-byte lane-release
  receipt binds this settlement root and the service-capacity update receipt.
  Retention conformance binds the same root, preventing a consumed slot from
  disappearing or an unused slot from being fabricated.

  The exact 345-byte `RetentionSourceShardReceiptBytesV3` binds attempt, epoch,
  branch arm, `store_role=Source`, admission-frozen source-store identity,
  TransferSet and source-shard-set roots, expected/committed source-inventory
  heads, capability, idempotency, signature set, exact shard artifact/byte
  counts, committed CAS, and durable interval. Source identity resolves through
  the frozen infrastructure roster to owner, implementation, medium, and
  failure domain. `RetentionSourceShardReceiptSetBytesV3` binds attempt, epoch,
  arm, TransferSet root, receipt count, total artifacts and bytes, followed by
  framed receipt digests; its exact size is `109 + 40*receipt_count`. The
  ordered duplicate-free source-shard sets must partition the TransferSet
  exactly: every required artifact occurs once, checked counts/bytes sum, and
  gaps, overlaps, unknown sources, or altered preimages fail. Each of exactly
  two pairwise owner/implementation/medium/failure-domain-disjoint destination
  receipts independently rehashes and persists the complete identical
  TransferSet. The infrastructure reservation row for each destination binds
  both the exact nonzero retention-policy digest and frozen destination-roster
  digest; the same two raw32 fields are canonical zero in every nondestination
  row. With those fields each framed member is 265 bytes and the exact nine-
  member set is 2542 bytes.

  Each exact 137-byte destination-retention member extends the prior 105-byte
  destination identity/receipt projection with the digest of that store's
  `CommitRetained` capacity update. The exact two-member destination set is
  `60 + 145*member_count = 350` bytes and binds both independently retained
  replicas. The exact 516-byte V3 durable-transfer receipt binds the V3 source-
  shard receipt set, destination-receipt set, retention-conformance receipt,
  and both lane-release receipts. Its acyclic order is source-shard and
  destination persistence -> all-and-only retention conformance -> lane-release
  authorizations and two committed lane releases -> durable-transfer receipt.
  No downstream `CommitRetained` digest is hashed backward into the 377-byte
  storage-write receipt.

  Historical V3 physical destination capacity remains exactly decodable for
  replay, but grants no fresh capacity authority.
  `StoreCapacityReservationIntentBytesV3` is exactly 257 bytes and binds
  attempt, epoch, store identity, maximum reserved bytes, retention policy,
  destination roster, signed lease, idempotency, and signature set. The closed-
  intent proof carries and rehashes that exact preimage rather than accepting an
  opaque digest.
  `StoreCapacityReservationMemberBytesV3` binds reservation intent, attempt,
  epoch, tagged state, store, retention policy, destination roster, maximum held
  bytes, actual used bytes, expiry ticket, and idempotency. Its Held arm is 225
  bytes; its 289-byte Retained arm additionally binds the exact TransferSet root
  and retention-store receipt digest. The exact 211-byte V3 store head binds
  store/generation/capacity; active reserved/count/root; permanent closed-intent
  count/root; and quarantined reserved/count/root. Its active and quarantine
  roots are 256-level maps and its closed root is a permanent sparse tombstone
  set. Every head proves that Held maxima plus Retained used bytes plus
  quarantined bytes equal reserved bytes and do not exceed capacity.

  The exact proof bundle has one or two typed, framed, domain-ordered proofs and
  is at most 17150 bytes. ClosedIntent proof is 8248 bytes, ActiveCapacity proof
  is at most 8818 bytes for old/new 225/289-byte arms, and Quarantine proof is
  8304 bytes. `ReserveCommitted` carries closed+active,
  `ReserveExpiredUnknown` closed+quarantine, `ReserveRefused` closed only, and
  each ordinary later operation exactly its one changed-map proof;
  `QuarantineReconcileCommitted` atomically carries active and quarantine
  proofs because it moves the fetched physical state between both maps. Each
  exact 455-byte update receipt binds one of `ReserveCommitted|ReserveRefused|
  ReserveExpiredUnknown|CommitRetained|AbortRelease|ExpiryRelease|
  QuarantineRelease|QuarantineReconcileCommitted`, twelve fixed digest fields—
  intent, store, member-or-zero, expected/committed heads, proof bundle,
  evidence-or-zero, policy, roster, capability, idempotency, and signature—plus
  committed CAS and interval. Reserve always inserts the permanent
  tombstone and selects active Held, no capacity, or maximum quarantine;
  CommitRetained derives/replaces the unique Held predecessor, charges exact
  used bytes, and releases only `maximum-used`; ExpiryRelease is legal only
  after governed expiry/legal-hold/clock authorization; quarantine release
  requires authenticated recovery proving zero physical effect.
  `QuarantineReconcileCommitted` instead converts quarantine to the exact
  independently fetched Held/Retained physical state; a late physical commit
  can never be de-accounted while bytes remain held. `ReserveExpiredUnknown`
  evidence binds the signed lease, independent clock quorum, and last observed
  store journal. Before activation each destination reserves full
  `run_reserve_i`. The acyclic order is Reserve capacity update; exact
  377-byte V3 retention-store write receipt, which extends the 345-byte
  attempt-local form with only the Reserve update-receipt digest; CommitRetained
  update, whose Retained member binds that store receipt and TransferSet; then
  source/destination receipts and conformance, which bind both capacity and
  write receipts; lane releases; then the 516-byte transfer receipt. Store
  selection after
  quiescence, unreserved staging,
  same-failure-domain replicas, stale CAS, overcommit, policy/roster swap,
  early used-byte release, or a byte-free member is `IntegrityFailed`.

  The historical V3 control-plane settlement ledger is a separately
  capacity-headed, append-only artifact capped and pre-reserved at 67108864
  bytes. It is replay-only where the terminal V4 closure supplies a conflicting
  staging, transaction, or downstream-binding rule. Each exact
  entry is `131 + payload_byte_count` bytes and binds attempt, epoch, entry kind
  `PreSettlementEvidence=1|SettlementMarker=2|PostSettlementReceipt=3`,
  predecessor control head, payload digest/count, and exact framed payload. The
  exact 121-
  byte ledger head binds cleanup-service identity, generation, entry count,
  total bytes, and latest entry digest. The complete pre-settlement payload E
  retains attempt/issuance, service/quarantine, capability/revocation/
  quiescence, requirement/source-set, shard/store/conformance, lane/store-
  release, and branch evidence. Exactly 4096 bytes remain reserved for the
  marker, settlement receipt, H1/H2/H3 head preimages, generic cleanup-capacity
  receipts/proofs, and closure, so E is at most 67104768 bytes.

  The exact 280-byte marker binds E, expected/committed control heads, old/new
  generations, capability, idempotency, signature, committed CAS, and interval;
  it does not bind the not-yet-created settlement receipt. Relative to the exact
  518-byte V2 settlement form, the exact 582-byte V3 receipt appends the H2
  control-ledger-head digest and marker digest. The only acyclic order is H0 ->
  append E as `PreSettlementEvidence` -> H1 -> append the marker(E) as
  `SettlementMarker` -> H2 -> commit global governor settlement(H2, marker) ->
  append the exact receipt bytes as `PostSettlementReceipt` -> H3. No record
  hashes its own successor. The global run reserve cannot release before H2.
  The cleanup capacity head is not a 131-byte shadow: it is the exact 228-byte
  generic infrastructure-service head for `CleanupAuthorityStore`, and its
  physical hold cannot release before H3. H3 plus the governor service's durable
  receipt journal preserve the settlement preimage even if cleanup return later
  fails. Preactivation uses the same ledger despite zero scientific transfer.
  Cleanup artifacts never masquerade as scientific results; inability to append
  after a valid reserve is retained infrastructure catastrophe and the reserve
  remains held.

  **Final V4 finite-capacity, WAL, and closed-evidence authority.** The final
  durable-allocation, coordinator/staging-capacity, and wire/accounting
  corrections in the authored I13 policy are the terminal fresh authority for
  I13 capacity. They supersede every conflicting V3 description and every
  earlier V4 draft value, including the 16384-slot decision journal, 394/604
  request/receipt pair, unpublished aggregate outcome, cyclic orphan member,
  330-byte transaction-chain step, and unpaged staging log. The V3 append-slot,
  lane-content, scientific-retention, and terminal-persistence clauses remain
  authoritative only where they do not conflict with this corrected V4 closure.

  The frozen finite roster has exactly fifteen logical authorities in fixed
  order: one Global lifecycle authority; seven `ReplayMetadata` authorities for
  Cleanup, the four Generic services, and Destination A/B; four physical
  Generic authorities for the two publication journals, publication verifier,
  and `OrdinaryEvidenceStore`; and three physical Store authorities for Cleanup
  and Destination A/B. Lane reservations remain separately governed lane
  authorities and are not a hidden sixteenth or seventeenth capacity authority.
  Each exact 161-byte roster row binds the authority kind and identity, both
  physical-copy identities, operation-protocol root, and transaction
  coordinator. The exact fifteen-row roster matrix is 2609 bytes. Its closed
  operation matrix uses exact 202-byte rows and the size equation
  `77 + 210*row_count`; each row freezes response/WAL maxima, repair budget,
  predecessor/successor dispositions, retry rule, and terminality rule before
  effect.

  The closed disposition machine is `Unused`, `Prepared`, `IntentDurable`,
  `Committed`, `AbortedNoEffect`, `ConflictLost`, `ReconciledCommitted`,
  `ReconciledAborted`, or `PendingRepair`. An exact-key retry returns the same
  terminal row, while a conflicting payload becomes `ConflictLost`. A one-copy
  append is `PendingRepair`; its reserved repair capability must write the
  byte-identical missing copy before transaction intent. Once either immutable
  Prepare is durable, it can never be discarded or inferred away: it must reach
  exactly one terminal disposition. If neither copy durably admits Prepare
  before any map effect, the API may return a structured, explicitly nondurable
  infrastructure refusal, but it grants no exact-replay or known-no-effect
  authority.

  Every logical authority is first allocated from two independent physical
  index copies. The exact 266-byte reservation subject binds its 262144-byte
  response-log quota, 131072-byte transaction-WAL quota, operation-protocol
  root, coordinator, attempt/epoch, policy, and idempotency. Each copy owns
  4096 nonborrowable entries of `262144 + 131072 = 393216` bytes, exactly
  1610612736 response/WAL bytes, and 24576 canonical-zero-padded 32768-byte
  decision slots, exactly 805306368 bytes. It additionally owns one
  canonical-zero-padded 131072-byte allocation-exception payload slot and one
  32768-byte authenticated sidecar reservation for each of the 4096 held
  ordinals, exactly 536870912 and 134217728 bytes. The fixed
  response/WAL/decision baseline is 2415919104 bytes and the complete per-copy
  allocation substrate is exactly 3087007744 bytes.

  The corrected 228-byte decision head's historical
  `reserved_subject_count` field means `held_subject_count`: the cardinality of
  the union of replay-derived reservation subjects and immutable
  `PreEffectAbsent` keys, in `0..=4096`. The first ordinary phase 1 or first
  atomic authenticated absence transaction for a subject assigns ordinal
  `old_held_subject_count`
  and reserves that ordinal, all six decision positions, and its exception
  slot. Recovery reuses the held ordinal. Joining both subject sources requires
  byte-equal ordinals on overlap, distinct ordinals for distinct subjects, and
  exact coverage `0..held_subject_count-1` with no hole. The exact invariants
  are `used_bytes = record_count*32768`,
  `record_count = sum(popcount(prefix_closed_phase_bitmap))`,
  `record_count <= 6*held_subject_count`, and
  `generation = record_count + absence_tombstone_count`. A legal zero-phase
  absence hold means `held_subject_count <= record_count` is not an invariant.
  Equality at every cap passes; cap-plus-one refuses before effect.

  The replay-derived reservation row is exactly 420 bytes: subject, held
  ordinal, prefix-closed six-bit phase bitmap, disposition
  `Pending|LiveAllocated|OrphanDominated|TerminalNoEffect`, and six ordered
  `(decision-entry digest, idempotency digest)` pairs, with canonical zeroes for
  unset phases. The reservation root is rebuilt from the retained at-most-24576
  decision entries and rejects duplicate subjects, ordinals, idempotencies,
  holes, and row/entry disagreement. The coordinator's historical
  `idempotency-map root` field is the domain-separated wrapper over that replay
  root, the authenticated absence-map root, and the fixed-depth exception-slot
  vector root, plus their checked counts and the slot vector's total used
  bytes. The combined reservation/absence subject-to-ordinal map is an exact
  bijection to `0..held_subject_count-1` and hence to the exception slots.
  The exact 272-byte copy-coordinator head is the sole mutable `CopyAtomicHead`
  CAS authority. Decision, index, orphan, absence, and exception-slot heads are
  immutable content-addressed projections selected by that one successor, not
  independently committed authorities.

  Availability is replayable from retained preimages and an explicit writer
  authority. The capability, signature, quiescence, storage, exception, and
  capacity figures in the following intermediate-design paragraphs are retained
  only to explain the defects corrected by the terminal closure below; they are
  non-normative wherever that closure differs. Each exact 163-byte
  writer-registry row binds logical authority,
  copy, writer ordinal and identity, capability-service identity, role, and
  policy. The frozen roster/protocol cross-product generates all and only
  `R <= 480` rows; its registry is exactly `136 + 171*R` bytes. The capability
  service has an exact 193-byte head carrying generation, capability epoch,
  registry root, active-roster root/count, policy, and sealed state. Every
  accepted write rechecks a capability bound to that service head and epoch.
  After all-and-only writer quiescence, the exact 437-byte epoch-fence receipt
  commits `generation+1` and `capability_epoch+1` in one service-head CAS,
  invalidating every stale writer capability.

  Every unqualified V4 `signature`, evaluator signature, or signature-set field
  is exactly a nonzero `raw32(CanonicalSignatureSetRootV4)`, never an inline
  64-byte signature. The complete canonical signer/key/algorithm/threshold/
  policy/revocation/object-binding set bytes remain content-addressed,
  independently rehashed, and separately admitted against the
  `OrdinaryEvidenceStore` reservation. Missing or unbudgeted preimages refuse;
  their physical bytes are not hidden in the maxima below.

  A 99-byte registry-projected runtime writer member and one 370-byte
  quiescence receipt per writer, one 364-byte stable-storage state, and one
  closed capability disposition form the exact availability bundle. The
  disposition is a 97-byte `NoCapabilityIssued` arm binding request, expected
  copy head, and current capability-service head, or a maximum 446-byte
  `EpochFenced` arm framing the 437-byte receipt. The admission-refused bundle
  has exact size `150 + FRAME(364) + FRAME(97) = 627` bytes at `W=0`; an
  unavailable bundle has exact size `976 + 485*W` for `1 <= W <= 16`, with
  maximum 8736 bytes. Writer and
  quiescence roots are separately domain-separated; equal counts plus an
  all-and-only ordinal bijection prove that every registry-selected writer
  quiesced. Admission refusal requires `W=0`, a zero request capability, no
  phase-1 replay, and storage not attempted. Unavailability requires `1..=16`
  writers, a revoked or expired old epoch, and the matching closed stable-state
  arm. The signed 508-byte fence evidence is accepted only when all roots,
  counts, service-head transitions, and digests recompute from these exact
  preimages.

  A `PreEffectAbsent` hold is one single-`CopyAtomicHead` transaction, never a
  pre-CAS append or a multi-head atomicity assertion. Its candidate payload is
  `tag || proposed_ordinal || request_digest || FRAME(fence) ||
  FRAME(preimage_bundle)`, exactly `559+B` and at most 9295 bytes. It contains
  expected inputs but no resulting absence member, candidate digest, committed
  head, or committed bit. The exact commit envelope frames that candidate, the
  8401-byte absence-map update proof, and the old 228/272-byte decision/copy
  heads; it is `9130 + candidate_bytes` and at most 18425 bytes.

  Exception allocation uses a canonical depth-12 vector with 4096 leaves. Its
  exact 180-byte nonempty member carries copy, held ordinal, the closed
  PreEffect/Phase3/Phase6 bitmap, `u32` used bytes, chain root, latest payload,
  policy, and sealed state. The self-describing vector update proof includes
  logical authority, copy, ordinal, policy, old empty-or-present arm, exact new
  member, old/new roots, and twelve ordinal-directed siblings. It is exactly
  786 bytes for `Empty -> Present` and 974 bytes for
  `Present -> Present`. The envelope is staged first, so its digest can enter
  the new member and proof without a self hash.

  The exact 384-byte pure commit receipt binds payload and vector-proof digests,
  expected/committed copy heads, generations, held counts, slot-used bytes,
  slot-chain roots, idempotency, and committed/completeness bytes. The durable
  absence artifact is `FRAME(envelope) || FRAME(vector_proof) || FRAME(receipt)`,
  at most 19619 bytes, under its own content digest. Envelope, proof, receipt,
  and artifact are durably staged before the sole copy-head CAS. A winner
  atomically owns the ordinal and updates the absence, decision, slot-vector,
  count/byte, and coordinator projections; a loser gives all staged bytes zero
  authority and changes no state. Crash recovery independently reconstructs
  the exact committed successor from the old head and both proofs.

  Phase-3/6 append unavailability retains the exact 508-byte fence, 431-byte
  deterministic projection, intended 756-byte outer decision entry, and
  availability-preimage bundle. Its record is `1728+B`, at most 10464 bytes.
  It is staged before a 786- or 974-byte vector proof and the same 384-byte
  receipt; their exact composite append artifact is at most 11846 bytes. One
  framed maximum PreEffect envelope plus one maximum unavailable record for
  each of phases 3 and 6 consumes
  `18433 + 10472 + 10472 = 39377` payload bytes, within the nonborrowable
  65536-byte ordinal slot. Vector proofs, receipts, and composite artifacts are
  separately charged evidence bytes. The slot recurrence and receipt bind the
  complete payload and proof; phase duplicates, wrong order, missing preimages,
  or cap-plus-one refuse. The matrix has no zero-durable initial row, so at most
  one phase-3 append can be unavailable globally; both then-local phase-6
  appends may be unavailable and repaired.

  Each copy's exact 272-byte coordinator head commits its index head, decision
  head, orphan head, expanded coordinator-state root, and generation. Every
  decision append is one local CAS of that sole coordinator head; the corrected
  228-byte decision head is an immutable successor projection selected by it,
  not a second CAS authority. Its exact 335-byte append receipt binds
  expected/committed projections and
  generations for both, old/new record and reserved-subject counts, old/new used
  bytes, and the exact entry. Every fresh phase requires both decision and
  copy-coordinator generation to equal the checked predecessor plus one; replay
  advances neither. Prepare and aggregate phases preserve complete index/orphan
  heads. A local-terminal mutation advances the selected component generation
  by exactly one with its proved root/count/bytes; a no-effect component remains
  byte-identical. Local-terminal phases may publish only the heads and proof
  frozen by their immediately preceding Prepare. The exact 395-byte local decision has an
  explicit intended orphan head and proof-bundle digest. Its unpadded entry is
  723 bytes; an exact 428-byte aggregate-outcome entry is 756 bytes. The
  decision-chain genesis/successor and every entry, head, coordinator,
  append-receipt, request, receipt, bundle, and outcome domain are normative.

  Allocation requests and receipts are exactly 426 and 636 bytes and bind the
  response-log, transaction-WAL, and decision-log identities independently. The
  subject-map member is 98 bytes and its update proof is 8330 bytes. The
  operation-matrix bundle has checked variable length. A clean two-copy Reserved
  trace has 39 framed items and exact bundle/pair/authority-set sizes
  30573/30788/31422 bytes. An ordinary six-phase Reserved trace has 63 items and
  exact sizes 40179/40394/41028, but those are not universal maxima. The final
  causal-fencing/exception closure below recomputes the recovered base,
  exception histories, signature material, direct snapshots, and page charges
  from their complete retained preimages. Its terminal maxima are a 74-item,
  198055-byte Reserved bundle, 198270-byte pair, 198904-byte authority set, and
  198542-byte response-log entry-zero wrapper. Terminal no-effect bundles have
  separately generated all-and-only rosters, a hard 262144-byte cap, and no
  pair, authority-set, or effect eligibility. Intermediate traces use checked
  exact sums, never a clean-size substitution. Profiles and every later
  consumer bind the authority-set digest, never the pair alone. The V4 terminal
  capacity-profile extension remains 531 bytes and carries three raw bootstrap
  authority-set digests.

  Allocation outcomes are the exact 428-byte closed sum `Reserved`,
  `RefusedNoEffect`, `PartialPending`, `PartialReconciledReserved`,
  `PartialReconciledNoEffect`, or `ExpiredUnknown`. Each copy projection is a
  `LocalDecision`-or-`PreEffectAbsent` sum and at least one must be a local
  decision. A separate three-arm legal-combination sum admits only both local,
  A-local/B-absent, or A-absent/B-local; two nondurable absences can never be
  laundered into durable evidence. A second A-only/B-only/both phase-3 receipt
  sum is the outcome-durability authority consumed by bootstrap evidence; a raw
  signed outcome alone is insufficient. The exact signed aggregate outcome is
  phase 3 in every available
  preprovisioned decision journal before any newly allocated response log can be
  trusted, preserving byte-identical first-Global replay after a crash. Global-
  log copies are bounded redundancy only. Partial or unknown allocation remains
  maximally charged until repair succeeds or an authenticated orphan-ordinal
  tombstone proves the unusable slot. The 132-byte orphan member binds the
  phase-1 AllocationPrepare decision and phase-3 initial pending outcome, never
  its proof-dependent phase-4/5 decisions or phase-6 terminal outcome. The
  member/proof/candidate root feed the durably appended phase-4 reconciliation
  Prepare; only its append receipt can enable phase 5, and only the phase-5
  append receipt can enable phase 6. Orphan heads remain 174 bytes and insertion
  proofs 8500 bytes. This ordering is acyclic; ordinals are monotone and never
  reused.

  The generated outcome matrix has all 882 rows in lexicographic
  `(phase, presence, local_A, local_B, durability_mask)` order and is exactly
  21271 bytes. Row ordinal is
  `((((phase_index*3 + presence_index)*7 + local_A)*7 + local_B)*3 +
  mask_index)`. Evidence code 2 means the complete durable absence projection,
  not a bare proof. Evidence code 4 means `RecoveredPrepareChain`, not the
  superseded recovery CAS/receipt: it contains the durable projection, exact
  368-byte authorization, 532-byte current-head request, 1784-byte recovered
  outer entry, committed phase-1 heads and ordinary receipt, and tagged direct
  or exception-repaired phase-3 completion with byte-equal nested fields.

  Recovery equality is defined by canonical resolvers, not by pretending every
  transitive field is serialized directly. `ResolveDurableAbsence(D)` fetches
  and rehashes the exact envelope/vector-proof/receipt artifact, rebuilds the
  old and committed decision/coordinator/absence/slot states, and checks one
  attempt, epoch, authority, subject, copy, role, ordinal, policy, and candidate.
  Recovery capability is deterministically derived from those direct context
  fields, including attempt, epoch, and role, plus `D`, the candidate, and the
  historical committed absence-state copy head.

  `ResolveAuthorization` checks only its directly encoded 368-byte fields and
  resolves `D` and the candidate through their authenticated digests.
  `ResolveRecoveredRequest` checks its direct context, capability, current
  actual predecessor, and `D`/authorization/candidate references. A harmless
  head refresh receives a new head-derived request idempotency; it need not
  equal the authorization's idempotency. `ResolveRecoveredEntry` checks its
  exact framed request/authorization, the nested LocalDecision request digest,
  retained durable references, and the phase-1 insertion at the same ordinal.
  Intervals nest and contain the winning CAS rather than being byte-identical
  across refreshed proposals. Per-copy replay guards reject a second committed
  entry carrying the capability; stale and failed proposals consume nothing.
  That CAS changes only decision/coordinator state. Phase 2 alone may publish
  the prepared index effect.

  The exact 173-byte transition-rules artifact contains six literal seven-byte
  rows, in allowed-next order:
  `Refuse=[0,0,0,0,0,0,0]`,
  `TerminalReserved=[1,1,0,0,0,0,1]`,
  `TerminalNoEffect=[2,1,0,0,0,0,1]`,
  `ReconcileRepair=[3,0,4,1,1,1,1]`,
  `ReconcileOrphan=[4,0,5,2,2,2,1]`, and
  `AwaitAppendRepair=[5,0,7,2,2,3,1]`. The byte positions are allowed-next,
  zero-final permission, required-final aggregate, final-presence rule,
  final-durability rule, evidence class, and terminality; every enum and
  NotApplicable filler code is closed. G0 exhausts all `882*883` initial-row /
  final-row-or-zero pairs and each evidence mutation.

  Append repair preserves semantic intent and exact bytes. An otherwise-
  terminal phase-3 outcome with only A or B durably appended retains the exact
  outcome and intended entry, including declared `Reserved` or
  `RefusedNoEffect`, copy projections, local outcomes, heads,
  orphan/reconciliation evidence, idempotency, signature, and interval. Only
  the derived matrix durability is one-sided and terminality is Nonterminal.
  `AwaitAppendRepair=[5,0,7,2,2,3,1]` means same declared aggregate, same
  presence, all then-local durable, retained append-repair evidence, and a
  Terminal repaired target. Exact replay of the byte-identical phase-3 entry
  changes only derived receipt mask and terminality. `PartialPending` and
  `ExpiredUnknown` remain genuine semantic uncertainty or mixed-effect states,
  not aliases for incomplete logging. Reserved repair reaches only the Reserved
  bundle/pair/authority path; repaired `RefusedNoEffect` reaches only terminal
  no-effect. An incomplete phase-6 append creates no admitted Final row; after
  repair, the original phase-3-to-phase-6 rule is evaluated normally.

  Phase 6 derives an all-then-local mask `ALocalBAbsent=1`,
  `AAbsentBLocal=2`, or `BothLocal=3`. Entry, head, and eventual-receipt sets
  contain exactly `popcount(mask)` objects in copy order. The repaired mask is
  a subset and selects exactly one history:
  `DirectAllThenLocal=0|RepairedA=1|RepairedB=2|RepairedBoth=3`; the last is
  legal only for `BothLocal`. Every repaired bit retains its maximum
  10464-byte record, vector proof, 384-byte slot receipt, and eventual ordinary
  receipt; every direct bit carries an explicit no-exception proof from the
  same slot-vector snapshot. Both Reserved and orphan/no-effect terminal
  authorities consume one complete history, so a two-copy final cannot
  terminalize after repairing only one of two unavailable appends and a
  one-copy final never fabricates its absent copy. The full state-specific
  capacity evidence bundle is a mandatory parent of the allocation-outcome set
  and bootstrap manifest, never an optional sink.

  **Terminal V4 causal-fencing and exception-publication closure.** This
  subsection is the final authority for capability, quiescence, storage,
  availability, exception publication, direct-history snapshots, and their
  capacity accounting. It supersedes the intermediate 193-byte service-head,
  fence-after-quiescence, single-signature-set, 8736-byte availability,
  65536-byte slot, 786/974-byte proof, and 96945-byte bundle descriptions
  above. Established allocation request/decision/receipt/outcome/fence and
  recovered-entry sizes remain unchanged.

  Capability and storage heads belong to independent authenticated services;
  the 272-byte `CopyAtomicHead` remains the sole mutable allocation-state head.
  Both service vectors are keyed by the stable admitted subject-slot ordinal,
  never the proposed held ordinal. Each capability subject has a dense
  sixteen-row issuance vector whose exact 195-byte members distinguish
  `NeverIssued`, `Active`, `Revoked`, and `Expired`. The exact terminal sizes
  are 343 bytes for a subject head, 217 for a service head, 515 for a subject
  update proof, 437 for a signed fence intent, 6571 for the complete sixteen-
  row issuance transition, 72 for nonissuance, `329+120*C <= 1169` for a
  three-to-seven-clock expiry quorum, and 581 for a nonrecursive fence commit.
  The complete `NoIssued`, `Revoked`, and `Expired` dispositions are exactly
  2790, 9289, and 10466 bytes. The capability CAS commits before quiescence;
  quiescence then proves every accepted old-epoch write drained. A commit
  record binds a candidate successor and `PublishOnHead`, never its own final
  service head; the final head is derived from the full record digest.

  Signatures are causal strata rather than one cyclic assertion. An exact
  `192+43*N` object set commits unsigned intent/authorization objects, an exact
  `195+77*K` canonical set commits at most eight policy-selected signers, and
  a 67-byte binding joins their roots. NoIssued/Revoked histories contain
  exactly four layers—capability intent, quiescence, storage intent, and final
  availability—while Expired prepends one clock layer. Their layer sets are
  489/564 bytes and complete material is bounded by 5949/7043 bytes. A layer
  signs only its pre-CAS intent or aggregate unsigned object. Proofs, commit
  records, and successor heads are derived afterward and become inputs only to
  later layers. The 508-byte availability evidence zeroes only its own direct
  signature field for signing; the enclosing 601-byte bundle embeds no object
  set and is assembled after the final signature exists.

  Digest authority is closed at the schema boundary. The 69-byte canonical-
  signature member and 178-byte coordinator core have their own literal,
  disjoint V4 domains rather than inheriting their enclosing set/root domains.
  Six deliberately retained short domain names—availability preimage,
  capability/storage commit, durable absence, signature material, and direct
  slot snapshot—have an explicit schema-to-domain alias table in the authored
  policy. Every other named V4 byte schema maps by its exact schema stem, and
  no two semantically distinct schema digests may share a domain.

  Writer roster members are 158 bytes and their all-and-only artifact is
  `317+166*W <= 2973`. Quiescence receipts are 477 bytes and the ordered set is
  `230+485*Q <= 7990`. Storage uses an exact 187-byte member, 212-byte head,
  841-byte Present-to-Present proof, 495-byte signed intent, reconstructible
  422-byte stable-extent projection, 691-byte nonrecursive commit, and
  2492-byte disposition. The preimage bundle contains the exact 601-byte
  context plus framed request, roster, capability disposition, quiescence set,
  storage disposition, and final signature-object set; its independently
  checked conservative cap is 26564 bytes. Every digest resolves to these
  retained bytes. The stable-extent digest is recomputed from the already
  retained old/new members and unsigned intent, not an unpriced opaque receipt.

  The exact coordinator core is 178 bytes: reservation root/count, absence
  root/count, exception-slot inner root/nonempty count/payload bytes,
  PreEffect-receipt-vector inner root/count, and global exception-chain
  root/count. Its contextual root preimage is 283 bytes. Slot and receipt
  vectors are distinct depth-12 trees with distinct nonzero empty/leaf/node/
  root domains. A 771-byte dual vacancy witness opens both old empty leaves.
  Compact slot proofs are 640 bytes for `Empty` and 828 for `Present`; the exact
  old 178-byte core opens their aggregate counts and payload total against an
  already-retained exact copy head.

  The PreEffect candidate is at most 27123 bytes. Its exact envelope adds the
  old core and framed dual witness and is at most 37218 bytes. The generic
  phase-tagged 384-byte exception receipt binds the expected old copy head and
  candidate nonrecursive coordinator state, never the receipt-dependent final
  head. Receipt digest then supplies the PreEffect receipt leaf and the exact
  199-byte global exception-chain step; only those successors produce the final
  coordinator root and one `CopyAtomicHead` CAS. Phase 3/6 preserve the
  PreEffect receipt vector and update the slot/global commitment. The durable
  absence artifact is at most 38266 bytes. Unavailable records are at most
  28292 bytes; core-bearing Empty/Present append artifacts are at most
  29526/29714 bytes. Every direct phase-3/6 arm carries a 453/641-byte slot
  snapshot proof plus the exact old core, for an at-most-835-byte artifact.

  Deterministic mode sequences logical copy mutations by frozen phase and
  subject/held-ordinal order; worker scheduling cannot change the global chain.
  Fast mode ledgers its order and makes no cross-schedule identity claim. The
  exception stage lease freezes request, stable subject slot, proposed held
  ordinal, and expected copy/decision heads. It becomes recoverably sticky
  after an external service commit until the same request publishes or an
  authenticated abort/tombstone completes. Exact retries reuse external
  records; conflicting retries refuse; no unbounded fence history is allowed.
  A genuine nonterminal phase-3 aggregate has a distinct append-repair arm that
  preserves `PartialPending|ExpiredUnknown` and remains nonterminal before
  phases 4–6. The otherwise-terminal repair arm still terminates. Legal initial
  presence has at least one local copy, so a Reserved attempt has at most one
  recovered copy and double absence issues no recovery authority.

  Each held ordinal has a 131072-byte payload slot and a 32768-byte
  content-addressed `OrdinaryEvidenceStore` sidecar under an exact charge
  lease. Worst payload use is 93826 bytes. A deliberately branch-independent
  sidecar ledger is 32144 bytes, including proofs, receipts, old cores, direct
  snapshots and retained heads, three signature-material/charge records, and
  composite framing. Referenced sidecars cannot be reclaimed; unreferenced
  losers drain before reuse. The resulting per-copy substrate is exactly
  `2415919104 + 4096*131072 + 4096*32768 = 3087007744` bytes.

  The conservative Reserved admission cap starts from the 40179-byte ordinary
  trace, adds one recovered-entry delta (1061), two historical head frames
  (516), one framed durable-absence artifact (38274), one Empty plus two Present
  exception-history items (88978), four framed maximum signature-material
  objects (28204), and one framed direct snapshot artifact (843). Therefore
  the bundle is at most 198055 bytes and 74 items; pair, authority set,
  response wrapper, raw page row, and framed page row are respectively
  198270, 198904, 198542, 199063, and 199071 bytes. The component caps are
  intentionally conservative and need not be jointly attainable. Terminal
  no-effect remains capped at 262144 bytes and is never effect-eligible.

  Each allocated copy then exposes a finite response log and transaction WAL.
  Response heads, entry fixed overhead, and pure append receipts are exactly
  214, 272, and 257 bytes; the checked prefix must remain at or below 262144
  bytes. WAL heads, entry fixed overhead, and pure append receipts are exactly
  185, 128, and 260 bytes; that checked prefix must remain at or below 131072
  bytes. Exact retained prefixes reconstruct every historical head and pure
  receipt, so no receipt-of-receipt sink exists. V4 does not inherit a V3
  aggregate replay-slot subtotal: each generated operation row admits its exact
  response and WAL prefix, including every entry wrapper and every
  manifest-located head preimage required by that row. The generated matrix's
  equality/cap-plus-one tests, not a legacy aggregate inventory, are
  authoritative; this contract asserts no additional universal per-authority
  byte total.

  Logical visibility uses one exact 221-byte transaction-coordinator head. Its
  separately domain-separated component and idempotency maps plus exact
  346-byte transaction-chain step produce the candidate head before WAL
  construction, without any backward or self hash. Mutation-component tags are
  `LogicalCapacityHead=1`, `ClosureJournalHead=2`,
  `GlobalGovernorHead=3`, and `ControlEvidenceStagingHead=4`. Ordinary
  operations carry exactly one tag-1 or tag-3 component; staging insertion and
  seal carry one tag-4 component under the already-rostered Cleanup
  `ReplayMetadata` authority; Global closure carries exactly tag 2 then tag 3.
  Response/WAL durability is a prerequisite to, not a fake component of, the
  one coordinator-head CAS. Readers accept only a state proved from the current
  coordinator component map; every lower-layer head remains Prepared data with
  zero effect until that sole CAS succeeds.

  The V4 Prepared payload is exactly `140 + primary_bytes`. Its distinct
  PreparedOnly primary sizes are 389 (`ReplayMetadata`), 440 (Generic), 455
  (Store), 431 (`ReplayArchive`), 566 (`GlobalAdmission`), 614 (`H2Settlement`),
  536 (`GlobalClosure`), 455 (`GlobalExpiry`), 460 (staging insertion), and 354
  (staging seal). Every historical CAS field in those primary projections is
  zero, including both historical Global-closure CAS positions; no V1/V2/V3
  decoder may reinterpret zero as Committed. Both response logs receive the
  byte-identical Prepared payload and all manifest-selected precommit evidence
  before the exact commit intent. The intent is the final signed and clocked
  nondeterministic object, exactly `370 + 138*component_count` bytes: 508 bytes
  normally and 646 for two-component Global closure.

  The chain step carries old/new committed-transaction counts as well as
  generations. A fresh CAS requires both successors to equal the checked old
  value plus one; exact-key replay advances neither. The exact update proof is
  at most 25100 bytes for two components. Both WALs receive the deterministic
  `676 + proof_bytes` record before the coordinator
  CAS. That record binds expected and candidate coordinator heads but no
  post-CAS fact. After the single CAS, the exact 359-byte transaction receipt is
  pure deterministic recovery evidence. The materialized durability-set
  projection has exact size
  `3299 + primary_bytes + 138*component_count + proof_bytes`; maxima are 20582
  for `ReplayMetadata`, 20633 Generic, 20648 Store, 20624 `ReplayArchive`,
  20759 `GlobalAdmission`, 20807 H2, 20648 expiry, 20653 staging insertion,
  20547 staging seal, and 29211 for two-component Global closure. It is a
  verification projection of exact response/WAL ranges, not duplicate physical
  storage.

  Precommit and postcommit evidence manifests are respectively
  `137 + 94*object_count` and `136 + 94*object_count` bytes. The precommit form
  names only already-existing objects. The postcommit form excludes its own
  durability set; its root is stored in that set, preventing a future-object
  cycle. Every effectful consumer field is generated all-and-only into exact
  133-byte authority-binding rows; the binding matrix size is
  `75 + 141*row_count`. `AdmissionCore` and issuance consume Global-admission
  durability, `CapacityOperationEvidence` consumes `ReplayMetadata` durability,
  physical and archive consumers consume their terminal durability sets, the
  normal Global-control closure consumes H2 durability, and every retention,
  release, archive, staging, and expiry field consumes its matching terminal
  set. A PreparedOnly projection or raw digest can never substitute.

  Bootstrap is derived from the roster rather than a guessed fixed count. For
  attempted roster prefix `1 <= K <= 15`, at most two immutable allocation
  outcomes exist per authority, so `A <= 2K <= 30`; preactivation mutation
  manifests satisfy `M <= 1 + 2*(K-1) <= 29`; and four fixed rows bind the
  profile, quiescence, lane-cleanup set, and dependency-trace root. Therefore
  `item_count = A + M + 4 <= 63`. Exact locators and items are 242 and 291
  bytes. The locator's two decision ordinal/digest positions select the
  per-copy phase-3-or-phase-6 aggregate-outcome entry (or zero only for a proved
  `PreEffectAbsent` copy), never a local-decision digest. At least one is
  nonzero, and its append receipt plus committed decision/coordinator heads
  resolve the complete state-specific trace. The evidence bundle is
  `61 + 299*item_count`, at most 18898 bytes.
  Bootstrap payload and Global-closure entry are at most 19321 and 19454 bytes.
  The first Global allocation is rooted in its two pre-provisioned decision
  journals; later outcomes remain authoritative in their decision journals and
  enter Global response logs only when a finite protocol row charged the
  redundant append. Bootstrap bypasses normal E/H2/H3 and feeds the two-component Global
  closure directly. Only canonical finalization and terminal-persistence
  catastrophe use the normal settlement ledger.

  Closure retention is executable rather than digest-only. A mutation-evidence
  archive is `135 + total_framed_object_bytes` and at most 65536 bytes. The
  exact Global closure-evidence artifact frames the at-most-19454-byte closure
  entry, old and new 179-byte closure-journal heads, at-most-29211-byte closure
  durability set, and at-most-65536-byte archive. Its total is
  `456 + closure_entry_bytes + closure_durability_bytes + archive_bytes`, at
  most 114657 bytes, and is self-contained for either Bootstrap or Normal.
  Normal control settlement uses the exact 1298-byte V4 postsettlement bundle,
  1429-byte entry, and 2489-byte H1/marker/H2/postreceipt/H3/footer tail, leaving
  63047 bytes. H3 cannot precede H2 durability, and Finalization/Catastrophe
  discriminators must byte-equal across every settlement projection.

  Expiry moves—not erases—closure evidence. Active `RunHeld` and
  `ClosureJournalHeld` entries have zero tickets; `Retained` carries its
  nonzero, branch-typed Bootstrap/Normal ticket. The Global lifecycle authority
  owns a nonborrowable `65536*131072 = 8589934592`-byte closed-evidence map.
  Its exact 176-byte head admits at most 65536 zero-padded 131072-byte slots.
  An exact closed member is `165 + closure_artifact_bytes`, at most 114822
  bytes; only its unpadded canonical bytes are hashed. The active-delete proof
  is at most 8457 bytes and the closed-insert proof at most 123062 bytes; their
  exact proof bundle is at most 131669 bytes. One Global-governor component
  successor commits active deletion and closed insertion in the same expiry
  CAS. Active bytes release only after the full executable closure artifact is
  in the committed closed slot; transaction durability and coordinator
  idempotency preserve exact replay.

  ControlEvidence staging is likewise finite and authority-bearing only through
  V4 durability. Exact 155-byte producer rows form the generated inventory
  `77 + 163*row_count`. Each conditional row contributes exactly one 116-byte
  `Extent` or one 100-byte `ZeroDisposition`. Extents distinguish the unique
  `PhysicalOwner` from any `ReferenceOnly` members; a separate owner map proves
  exactly one owner and every reference's membership. Zero dispositions never
  acquire positive physical layout. The exact 215-byte staging head is Open or
  Closed and binds member/owner roots, counts, physical-unique bytes, policy,
  and inventory. Member and owner proofs are bounded by 8348 and 8280 bytes;
  their bundle is at most 16716 bytes.

  Each exact 460-byte staging-insertion primary is PreparedOnly and gains
  authority solely from its at-most-20653-byte durability set. The exact
  prepared-payload, staging-proof-bundle, staging-artifact, and page-row headers
  and domains are normative; implementations may not infer them from prose.
  The staging artifact retains rows in generation order, never digest order.
  Each at-most-21261-byte row materializes the member, PreparedOnly projection,
  and durability set so replay survives response/WAL expiry. The exact 354-byte
  seal is also PreparedOnly; its at-most-20547-byte durability set is the sole
  Open-to-Closed authority. The pure 225-byte cutoff binds the closed head,
  artifact, seal durability set, inventory, and retained set.

  Staging uses finite deterministic log pages rather than overflowing one
  Cleanup response/WAL prefix. Each exact page row is
  `159 + allocation_authority_set_bytes`, at most 199063 bytes, and embeds that
  page's complete at-most-198904-byte allocation authority set, reserved response/
  WAL prefix maxima, and starting post-entry-zero heads—not post-transaction
  final heads. Its framed row is
  `167 + allocation_authority_set_bytes`, at most 199071 bytes. The raw page table
  is `10 + sum(framed_page_rows)` and its enclosing frame has
  `L = 18 + sum(framed_page_rows)` bytes. That exact frame immediately precedes
  generation rows in the artifact, whose size is `106 + L + 21261*N`. Let
  `Lmax` be the inventory-derived maximum of `L`.

  The generated inventory/protocol cross-product partitions the `Nmax` inserts
  into nonempty insertion pages and allocates exactly one distinct final
  zero-insertion `SealOnly` page before artifact construction. Starting with the
  first unassigned generation, an insertion
  page is the longest nonempty consecutive prefix satisfying both byte and count
  caps; the boundary is the earliest next invocation that would violate a cap.
  Each response prefix begins with its at-most-198542-byte allocation-pair entry;
  WAL accounting begins at zero and never charges that response-only entry. On
  each copy, decoded `old_reserved_subject_count + page_count <= 4096`, including
  every previously burned no-effect ordinal, and all six-entry decision
  reservations must exist. The seal page separately reserves entry zero plus the
  seal transaction; the pre-seal artifact contains only its starting heads and
  authority, while the downstream seal durability set/cutoff retains resulting
  heads. Every insertion transaction maps to exactly one nonzero page range.
  Sharing the seal with an insertion page, post-seal heads in the artifact,
  caller partitions, favorable ties, and a single invocation that cannot fit
  all caps refuse before effect.

  Control E lays out the staging artifact, seal durability projection, and
  cutoff as three special PhysicalOwner extents in addition to staged owner
  bytes. The no-page-table lower equation is
  `Smax + 21319*Nmax + 20547 + 66505 <= 67042864`; the authoritative paged
  equation is
  `Smax + 21319*Nmax + Lmax + 20547 + 66505 <= 67042864`.
  Equality passes and every byte, count, page, index, decision-slot, or prefix
  cap-plus-one refuses before producer effect. The coefficient 21319 is one
  21261-byte artifact row plus its 58-byte E-manifest row; 66505 includes fixed
  artifact/cutoff/E bytes, the three special manifest rows, and the
  nonborrowable 65536-byte post-E reserve.

  This V4 closure remains `[M]`. No transaction, closure, expiry, staging, or
  all-and-only certificate is promoted until G0/G3/G4/G5 prove exact bytes and
  domains; empty/first/max/max-plus-one finite structures; component and
  binding all-and-only generation; crash at every index/decision/response/WAL/
  coordinator boundary; one-sided repair and every Prepare terminal;
  phase-3 aggregate-outcome survival before response-log creation; the exact
  phase-3 -> orphan member/proof/candidate root -> phase-4 -> phase-5 -> phase-6
  order;
  reserved-subject partition/accounting; partial-allocation/orphan recovery;
  absence of backward hashes; bootstrap
  bounds and separation from normal H2; final/catastrophe arm isolation;
  closure replay after active-log expiry; staging ownership, zero arms,
  page cross-product/allocation authority, generation replay, and post-seal
  refusal; and every capacity equality and cap-plus-one twin.

  Each attempt also reserves one terminal segment. Its single emergency slot
  has the exact `I13_EMERGENCY_TERMINAL_SLOT_V2` header and fixed 368-byte
  overhead around the exact 764..=1048576-byte final envelope, so maximum slot
  size is 1048944. It cross-binds attempt/epoch/lane-epoch/idempotency, emergency
  lane, reservation/permit/capability, primary failure, selection/drain, expected
  terminal and witness roots/counts, core/grant/envelope, and frozen recovery/
  redaction policies. Missing failure witness, mixed lane/core/policy, or stale
  identity aborts without consuming the only slot.

  Under an Emergency grant, one indivisible transaction commits slot, emergency
  operational head, logical terminal root, witness closure, and
  `execution.finalized`, or none has canonical or journal effect. The downstream
  emergency receipt is 288 bytes. Maximum slot plus the 288-byte journal and
  388-byte finalization receipts is 1049620 bytes, leaving 3144684 bytes in the
  terminal segment. Receipt fields bind directly or transitively through the
  slot; both receipts remain downstream and recover under the separate ledger.

  Admission budgets the full primary and emergency finalization routes. Primary
  includes prepare/submit, queue, and durable service. Emergency additionally
  charges spent primary work, qualified failure detection/receipt, lane fence/
  core rebuild, emergency prepare/submit, queue, atomic terminal service, and
  clock uncertainty. Later primary copyback has a separate availability/
  retention bound and no OnTime/finalization authority. A qualified primary
  failure forces a fresh emergency core and grant. Before the all-or-none
  emergency commit, state is `TerminalPersistencePending`; afterward it is
  canonical and there is no reconstruction/recommit phase. Only failure of both
  mirrored prefixes, the emergency transaction, and every authenticated replay
  route permits `TerminalPersistenceCatastrophe`, never a favorable result.
- `i14_draft()` — the I14 multirung EMC/harness gate: 28 claims split 3
  baseline `[S]`, 12 `[F]`, and 13 `[M]`; 53 fixture, convention, theorem,
  falsifier, laboratory, and population cards with 20 stage-local held-outs;
  28 one-claim execution leaves (13 Core and 15 Max); and 7 independently
  scoped governed waivers. Native HarnessGraph and synthetic AP242-adapter
  mechanics, RLGC admission and MTL propagation, source/probe and victim-mode
  semantics, UQ mechanics and governed reliability, BEM formulation and FMM
  acceleration, laboratory EMC and bearing-population evidence all have
  separate authority. One `EmConventionCard` binds phasor, phase, RMS/peak,
  material, outgoing-wave, port, Poynting, PML, BEM, probe, and adjoint signs.
  Acceptance arithmetic keeps exact hard predicates disjoint from explicitly
  soft tolerances and requires certified outward enclosures. Maximal ratchets
  target FEEC stability, certified fidelity descent, passive-causal
  sheaf/cosheaf composition, hypercohomology-obstruction localization, cover-
  refinement naturality, a KYP/sheaf passivity bridge, robust mitigation, and
  exhaustive countermodel discovery. Governed blind/physical/standards/
  population authority requires an atomic staged transaction that installs a
  same-ID typed external discharge envelope for every retired waiver and a
  separately verified receipt; an arbitrary digest or waiver deletion is only
  a syntax change and cannot promote. Blind and population evidence use the
  strict `GovernanceCommitted` -> `CandidateFrozen` -> independent custodian
  realization -> atomic `RealizationCommitted` ->
  `RevealedForAdjudication` -> `Closed` lifecycle; commitment-slot replacement,
  same-ID envelope installation, receipt verification, waiver retirement, and
  authority-head advancement are one transaction. Physical validation uses a
  stronger actor-scoped protocol. Before candidate-visible calibration or
  model-input access, an independent custodian/fs-vvreg capability creates
  hiding calibration/model-input commitments, a hiding validation source-
  universe commitment, disjoint-membership commitment, exact selection
  algorithm, and pre-candidate seed/VRF commitment. Candidate builders,
  fitters, checker/threshold owners, and their transitive capabilities receive
  only opaque commitment identities before `CandidateFrozen`: no validation
  bytes, membership witnesses, labels, aggregates, derived statistics,
  selection outputs, or opening material. Named builders may use only the
  committed calibration/model-input strata under a complete derivation/access
  ledger; the frozen candidate is then followed by a contamination receipt
  naming every audited candidate-side principal/capability, independent
  validation realization, membership/disjointness/non-adaptive-selection
  proofs, atomic joined-root discharge, controlled reveal, and closure.
  Licensed standards use the explicit `GovernanceCommitted` ->
  `AuthorizedConstruction` -> `CandidateFrozen` ->
  `StandardsAuthorityCommitted` -> `StandardsAdjudicated` -> `Closed`
  lifecycle with scope-before-access construction and a pre-adjudication
  same-ID envelope/waiver transaction. Every lifecycle is cancellation-safe:
  protected access, reveal, publication, adjudication, and authority-head
  advancement are forbidden after cancellation unless that exact atomic stage
  had already committed. Receipts bind predecessor plus transaction intent,
  while the amendment record binds predecessor and final successor, avoiding
  a content-hash cycle. The 23-field I14 intent projection represents every
  not-yet-realized envelope, governed output, final successor, and amendment
  record through an explicit `Pending` digest union; it never fabricates a
  digest sentinel. It may bind the already-existing predecessor-stage receipt,
  but excludes every newly created discharge/commit receipt, realized output,
  final-successor, amendment-record, and new-authority-head digest.

  The public I14 terminal-status, selection, lifecycle, canonical-result, and
  telemetry slice described below is executable. The governed 23-field
  transaction intent, schema-set/output/head/commit-receipt algebra, authority
  CAS, and atomic-ledger path remain authored contracts with
  `NoPromotionAuthority` until their separately stated decoder, independent-
  encoder, receipt-authentication, and atomic-commit proof gates are supplied.
  `I14TerminalStatusV1` and `i14_evaluate_terminal_status_v1` retain the raw
  producer tuple, expose every fail-closed normalization action, normalize all
  3,600 eight-axis tuples, and pin raw tags, normalized tags, normalization
  bits, and exit projection under a domain-separated table digest. The scoped
  `i14_select_terminal_boundary_v1` cause selector caps requests, scope depth,
  and observer tiles; validates globally unique causal sequences, timestamps,
  deadlines, strictly prior boundary observations, acyclic scope ancestry,
  and admitted observer identities; canonicalizes request and observer-catalog
  presentation order before malformed-input diagnosis; and implements the
  frozen cause precedence: `InfrastructureFailed > TimedOut > Cancelled >
  BudgetExhausted > Completed`. An earlier pending cancellation defers only a
  normal completion/budget candidate, not a boundary with no terminal cause.
  This V1 path proves only local cause arbitration at one caller-supplied
  boundary. It is retained for old-ledger readability and has no
  first-terminal, cancellation-SLO, lifecycle, or promotion authority.
  `I14CanonicalTerminalResultInputV1`,
  `i14_canonical_terminal_result_v1`, and
  `i14_canonical_terminal_result_digest_v1` refuse a nonterminal boundary or a
  normalized-receipt/causal-disposition mismatch. The resulting canonical
  identity binds the strict pre-boundary logical trace cut, cause candidates,
  selected cause, request identities/sequences, boundary ordinals, raw and
  normalized terminal axes, normalization actions, exit projection, and an
  explicit content-addressed semantic payload. Raw monotonic times, deadlines,
  watchdog arrivals, and clock calibration are excluded. Within the bounded
  recorded trace, valid or malformed post-boundary records cannot rewrite the
  frozen cut.
  `i14_telemetry_envelope_digest_v1` separately validates the complete recorded
  trace and binds all raw timing/deadline fields, at most
  `I14_MAX_WATCHDOG_OBSERVATIONS_V1=4096` canonically ordered watchdog records,
  and the clock-calibration artifact to the canonical result digest.
  `I14_CANONICAL_TERMINAL_RESULT_V1_KAT_HEX` and
  `I14_TELEMETRY_ENVELOPE_V1_KAT_HEX` pin both version-1 byte layouts. G5 bit
  stability and promotion do not rely on those legacy layouts.

  The schema-authoritative V2 path closes the omitted-prefix, self-selected-
  deadline, and caller-selected resource-cause gaps. Schema authority is not
  promotion authority: `fs-vmanifest` validates and hashes a projection, while
  the consuming HELM/ledger gate must authenticate the bound issuer,
  capability, trust policy, revocation state, exact card, and semantic-work
  verification receipt before promotion. `i14_admit_cancellation_card_v2` is
  the sole constructor for `I14CancellationCardV2`: it admits Core, ordinary
  Max, or the explicit `MaxTheoremFalsifier` Max subtype and binds the semantic
  work unit, wall and hard logical-memory ceilings, exact count-resource kind,
  total ceiling, logical tile quantum, resource authority, deterministic
  partition, execution-environment fingerprint, response bounds, and optional
  exact external-child catalog. It preserves the frozen 90-minute/18-hour/24-
  hour and 32/96-GiB envelopes; enforces the per-kind Core ceilings
  4096/16384/1024/256 and tile quanta 64/256/16/4 (four times each for Max);
  and refuses zero, over-tier, unit-inconsistent, or external-policy-
  inconsistent contracts. Core/Max/theorem-Max cards fix observer and in-
  flight-child caps at 32/128, request-to-observation at 250/1000 ms, watchdog
  quantum at 25/100 ms, and drain/finalize SLOs at 2/8 seconds. A multi-kind
  campaign uses one card per independently governed execution leaf.
  `i14_select_first_terminal_boundary_v2` requires a nonempty trace beginning
  at genesis ordinal zero, at most 4096 contiguous boundaries, strictly
  increasing logical sequences, nondecreasing calibrated times, an immutable
  scope/observer catalog, exact card-derived request deadlines, bounded child
  and watchdog state, monotone resource consumption/work-frontier state, and
  no record after the first selected terminal. Cumulative consumption may not
  exceed the hard ceiling. `BudgetExhausted` is derived only when frozen work
  remains and a receipt-bound nonzero next-work quantum no wider than the card
  tile quantum is rejected because it would cross the ceiling; completion
  exactly at the ceiling remains `Completed`. The reference implementation
  refuses more than `I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2=1048576` boundary/
  request pairs before its repeated V1 arbitration, preventing an admitted-
  maxima sort/alloc denial of service while an incremental equivalent selector
  remains proof-pending. It
  returns either an opaque first-terminal selection or an opaque nonterminal
  frontier certificate; callers cannot manufacture either proof type or turn a
  frontier into a terminal result. Temporal precedence across boundaries is
  therefore separate from same-boundary cause precedence.

  Every successful trace recomputes its clock-free logical root from canonical
  boundary, resource, work-frontier, decision, request, and observation
  semantics; no caller-supplied hash can substitute for those records. The
  request-inclusive terminal-prefix digest also binds the independent
  verification-receipt identity. `I14CanonicalTerminalResultInputV2` cannot
  bypass the trace selector. It also
  requires `I14TerminalLifecycleTraceV2`, whose typed fields make execution
  start, drain start, drained, and finalized events mandatory; bind active and
  drained child counts; split the unconditional child lifecycle into a
  canonical semantic root, telemetry-only raw root, and independent canonical
  verification-receipt identity; and carry the typed
  `I14DrainTriggerV2::{CancellationObserved, ObservationTimeoutDrain,
  InfrastructureFailure, NonCancellationDrain}`. The trigger is derived
  deterministically from the earliest effective in-scope, pre-drain on-time
  observation, the first nanosecond after a missed inclusive observation
  deadline, or a receipt-bound infrastructure-failure onset whose receipt is
  authenticated by HELM/ledger. Equal effective
  times are broken by the frozen causal tie rank
  `Infrastructure=0, Observation=1, Timeout=2` (which is deliberately distinct
  from wire tags), then causal logical sequence and stable identity. A non-cancellation drain
  starts its SLO clock at drain start. A receipt-bound infrastructure onset may
  coalesce with drain start when both logical sequence and calibrated timestamp
  are exactly equal; it remains the infrastructure trigger rather than being
  relabeled non-cancellation. Canonical bytes bind only the derived
  variant plus request ID or infrastructure-onset logical sequence and, for an
  infrastructure witness, its closed source tag and independent verification-
  receipt digest; calibrated trigger/onset times remain telemetry. Local
  evidence-derived sources must agree with the corresponding failed evidence;
  generic supervisor/authentication/drain/publication sources are structurally
  bound here and authenticated only by HELM/ledger.

  `I14SpawnFrontierEvidenceV2` is mandatory on every terminal path. Its request
  identity names the earliest actual pre-drain observation, so it may differ
  from a timeout drain trigger; it is `None` when drain start closes a frontier
  with no prior observation. Request submission and a calibrated deadline do
  not synthesize a logical spawn cut. The typed last-spawn event participates
  in global logical/time ordering, its scheduler semantic root must equal the
  child semantic root, and its count reports every spawn at or after frontier
  closure. Watchdog, tile-poll, and external-heartbeat evidence likewise split
  canonical semantic roots from telemetry-only raw roots and independent
  verification-receipt identities. Typed summaries cover before-item-zero and
  before/after-tile polling, per-catalog heartbeats, external terminate/drain
  acknowledgements, and atomic no-partial publication.

  This layer checks structural/internal consistency, caps, trigger arithmetic,
  failure reflection, and—when `watchdog_samples_complete=true`—the complete
  multi-kind raw watchdog stream against its poll count, endpoints, maximum
  gap, and versioned raw-root derivation. A bounded diagnostic subset has no
  completeness authority. This layer structurally binds semantic roots and
  verification-receipt identities but neither authenticates the corresponding
  receipts nor grants external semantic authority. The HELM/ledger promotion
  gate must authenticate those receipts and verify issuer authority, membership,
  sequence completeness, child/external
  acknowledgements, publication atomicity, and every raw-event derivation.
  A receipt-bound first infrastructure-failure onset latch must select the
  first logical boundary on or after onset. Timeout onset is never caller-
  selected: it is the earliest first nanosecond outside the inclusive trigger-
  to-drained or drained-to-finalized cap. A boundary at the inclusive deadline
  itself remains on time; the first boundary at or after the derived
  deadline-plus-one-nanosecond onset must latch timeout. The canonical timeout
  logical field is that exact latch-
  boundary sequence, not a synthetic caller-selected onset event. Real watchdog,
  tile-poll, external-child, descendant-drain, or spawn-after-frontier failure
  remains admissible evidence only when the terminal cause honestly reports
  `InfrastructureFailed`; real trigger-to-drained or drained-to-finalized SLO
  failure remains admissible only when `TimedOut` or higher-priority
  `InfrastructureFailed` reports it.

  The V2 canonical digest binds the recomputed clock-free logical root and an
  independent verification-receipt identity whose corresponding receipt the
  HELM/ledger gate authenticates. It also binds the request-inclusive prefix,
  cancellation card, every genesis-to-terminal
  decision, selected V1 semantic projection, lifecycle logical sequences and
  typed trigger, child semantic/verification identities, mandatory spawn-
  frontier audit, clock-free tile/watchdog/external semantic/verification
  identities, the receipt-bound infrastructure-onset latch, the locally derived
  timeout latch, and every derived failure bit. Child/watchdog/tile/heartbeat
  raw roots, raw poll/heartbeat counts, and
  calibrated times remain telemetry-only.
  `I14CanonicalLifecycleProjectionV2` and the canonical result expose immutable
  getters for the validated lifecycle semantics rather than hiding them. The
  V2 telemetry digest separately binds every raw campaign, boundary, request,
  observation, lifecycle, watchdog, external-count, and calibration field; a
  completeness bit distinguishes a fully reconciled multi-kind watchdog stream
  from a bounded diagnostic subset. `I14LateEventTailV2` retains a post-terminal
  logical tail bound to an independent verification-receipt identity in
  telemetry without permitting it to rewrite final disposition; receipt
  authentication remains a HELM/ledger responsibility. Thus timing shifts that
  preserve all semantic
  relations preserve canonical identity but change telemetry; crossing a
  deadline or causal boundary changes both. Four checked-in fixture digest
  KATs pin their byte layouts through the current local implementation:
  `I14_CANCELLATION_CARD_V2_KAT_HEX`,
  `I14_TERMINAL_PREFIX_V2_KAT_HEX`,
  `I14_CANONICAL_TERMINAL_RESULT_V2_KAT_HEX`, and
  `I14_TELEMETRY_ENVELOPE_V2_KAT_HEX`. Three additional encoding/digest KATs
  pin otherwise easy-to-drift unions: `I14_DRAIN_TRIGGER_ENCODING_V2_KAT_HEX`
  exhausts all four trigger tags/payloads,
  `I14_INFRASTRUCTURE_FAILURE_ONSET_ENCODING_V2_KAT_HEX` pins the present
  presence/sequence/source/receipt form, and
  `I14_WATCHDOG_RAW_TRACE_V2_KAT_HEX` pins the known-answer digest of the exact
  629-byte complete 35-record multi-kind raw-watchdog fixture. These constants
  and same-implementation checks are drift detectors, not evidence that an
  independent encoder exists or agrees. Before promotion authority, two
  independently implemented encoders must first match the exhaustive terminal-
  table bytes and then reproduce every exact encoding/digest KAT; this version
  does not claim that gate has passed. G5
  bit stability applies only to this V2 canonical result under an identical
  logical event/cause trace, identical bound verification-receipt identities,
  and identical topology/ISA/toolchain fingerprint.
  `I14ArtifactCategoryV1` and `i14_retention_rule_v1` exhaustively assign all
  sanitized evidence/failure and raw licensed, secret, specimen, governed-
  holdout, derived-sensitive, and diagnostic categories to typed durability,
  sanitization, encryption/capability, access-ledger, and retention/erasure
  requirements.
- `i15_draft()` — the I15 (executable standards compiler) instance:
  9 claims (5/2/2 with a transitive-impact completeness falsifier and a
  human-authority preservation moonshot), 8 fixture pins (2 held-out), 6
  obligation rows, 1 waiver. LICENSING BOUNDARY: all I15 fixtures are
  synthetic standard-shaped packs; no licensed standard text is embedded;
  real editions enter only through the waived external slot pinned via
  fs-vvreg.
- `rl_draft()` — the RL reality-loop PORTFOLIO aggregate (the first
  CP/EM/RL/PD aggregate): 12 claims split 7 baseline `[S]`, 2 `[F]`, and 3
  `[M]`; 12 fixture, corpus, theorem, and protocol cards with 5 stage-local
  held-outs (3 Core, 3 Max ranges with the mutant block spanning
  135168..=143359); 5 execution leaves; and 1 governed-physical-loop-pack
  waiver. RL freezes the SEAMS of the lot-to-experiment-to-deployment loop
  composed from the I05/I06/I10/I11 instance authorities, never their
  internals: one identity spine (stage boundaries bind exact content
  addresses; names/paths/timestamps carry no authority); end-to-end
  blind-partition custody with side channels in scope; exactly-once
  uncertainty ownership whose twin-level decomposition recomposes exactly;
  a frozen calibration/measurement graph with cone-exact quarantine;
  bitwise whole-loop replay from the EvidenceRetentionReceipt chain;
  declared-graph selective reproof (a consumed instance-manifest successor
  version invalidates every receipt bound to the predecessor digest); and
  weakest-wins deployment gating with a closed typed-refusal taxonomy.
  Frontier lanes target the anytime-valid adaptive OED loop and
  target-exact HIL/timing composition (a measured sample maximum never
  promotes to a WCET bound anywhere in the chain). Moonshot lanes target a
  machine-checked end-to-end composition theorem and physical-campaign
  receipt parity, while a hidden-mutant lane attacks forged cross-stage
  certificates. Version-1 prose grants no composition, parity, or
  deployment authority.

## Invariants

- FAIL-CLOSED FREEZE: a manifest missing any load-bearing field or required
  list item, carrying no claims, with
  duplicate ids, a non-independent oracle, an invalid tolerance, a
  malformed fixture, a duplicate claim mapping, an orphan reference or
  waiver, an uncovered claim, or a claim/leaf evidence-id collision
  cannot freeze; the refusal names the gate.
- SEALED AUTHORITY: `FrozenManifest` is immutable by construction;
  "alter flags after freeze" has no code path. Change is amendment,
  amendment is a same-initiative successor version, exhaustion is a
  typed refusal rather than integer wrap, and the record names the exact
  reverse-dependency invalidation set — nothing else is invalidated.
- VERSION AUTHORITY: the numeric manifest revision advances exactly once
  per amendment and is the only revision field interpreted by the schema.
  Authored constructors are checked only against the documented,
  non-exhaustive known legacy semicolon-field spellings; arbitrary public
  drafts receive no semantic natural-language parsing claim.
  A version-only successor can carry an empty invalidation set; an
  identity-bearing campaign-authority title or campaign-policy change cannot
  masquerade as version-only.
- CANONICAL IDENTITY: components sort into one total order with content
  tie-breaks; assembly/input order and presentation order of declared set
  fields can never move `digest()`; valid external hex normalizes to raw bytes
  so case cannot fork identity; exact floating-point bits govern both digest
  and component equality; `FrozenManifest` equality is digest equality; every
  semantic field of every component is mutation-sensitive. Frozen manifests
  expose accepted obligation rows only through the canonical owned-set
  projection, never raw authored order.
- LATTICE SEPARATION: solid/frontier/moonshot claims are distinct
  elements; a weaker receipt closes its own element and is never
  relabeled as the stronger theorem (the I01 maximal lanes activate only
  after the baseline closes, per their activation fields).
- BOUNDED GATES: collection, per-row list, and cumulative UTF-8 text-byte caps
  are checked before any semantic content scan (`MAX_CLAIMS`/`MAX_FIXTURES`/
  `MAX_OBLIGATIONS`/`MAX_WAIVERS`/`MAX_ROW_ITEMS`/
  `MAX_MANIFEST_TEXT_BYTES`). The byte cap covers text payload only, not an
  overclaimed bound on allocator or canonical-framing overhead.

## Error model

Total functions; no panics in library paths. Freeze failures are typed
`FreezeRefusal` values; amendment failures are `AmendmentRefusal` values
(initiative change, exhausted/wrong version, adjacent-version evidence-kind
alias, or the successor's own refusal). Generic structural seed defects are
freeze refusals; semantic defects specific to an initiative, such as an
incorrect Philox alias/range, are focused conformance or campaign-integrity
failures, never silently normalized and never crash a library path.

## Determinism class

Fully deterministic: seed data and governed-sampling protocol text are
`const`/static; digests are
domain-separated BLAKE3 over length-framed canonical bytes with exact
float bit patterns. Byte-stable across runs and thread counts on the
same ISA (the G5 test); cross-ISA stability of the digest is expected
but not yet claimed — see no-claim boundaries.

## Cancellation behavior

None; freezing is synchronous and pure. Cost is `O(n log n)` sorting plus at
most two content-digest evaluations per claim/fixture/waiver and three per
obligation (manifest grouping, content-tie ordering, then sealed-projection
digest caching), with bounded cross-reference scans that are
worst-case quadratic in capped component/list counts. Amendment propagation
uses a deterministic set of predecessor authority ids, so fan-out cannot
amplify duplicate owned strings before finalization; its reverse-dependency
scans remain within the same caps. Chunked in-memory manifest assembly
identity is covered by clone-boundary equivalence tests: identity depends only
on frozen content, never on append history. This is a G4 precursor, not
durable checkpoint encoding, process restart, corruption recovery, or runtime
request-drain-finalize proof.

## Unsafe boundary

None. Workspace `deny(unsafe_code)` applies.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`: the I01 draft freezes with the declared lattice
split (5/2/2), held-out partitions, full coverage, and a refutation-
polarity falsifier lane; the freeze-gate battery in documented order
(all four collection caps, all six per-row/list caps, and cumulative text cap
before semantic scans, version, blanks, required
nonempty claims, duplicates, category-global oracle and reference precedence,
tolerance cases, malformed fixtures, orphan refs/decks, uncovered claims,
waiver-covered claims accepted); G5 input-order invariance of the
digest; chunked in-memory assembly identity as a G4 precursor; the G3 mutation battery
(loosened bands, weakened hypotheses, swapped held-out fixtures move the
identity; production-oracle reuse fails closed; post-freeze alteration
has no code path); compile-time-maintained, every-field component mutation
sensitivity for claims/oracles, fixtures, obligations, waivers, and every
tolerance payload/tag, including exact-bit signed-zero distinction; direct
external-hex case-normalized equality; canonical draft accessors plus sealed
canonical frozen obligation projection/digest; amendment semantics (initiative/version/
overflow/evidence-kind gates, defective successors, claim/fixture/filled-waiver-slot/
waiver discharge/obligation/waiver/title/global-policy propagation,
version-only and new-authority preservation, and mapping-set order
invariance); the I02 draft freezes with its own 5/2/2 lattice, 2 held-out
partitions, full coverage,
refutation-polarity tearing falsifier, an identity distinct from I01's,
and I02 input-order invariance; and sibling instance smoke checks cover
I04, I08, I12, and I15.

`tests/i05.rs`: exact I05 7/2/3 lattice, single maximal refutation lane,
15-fixture/6-held-out corpus and once-only seven-leaf claim map; all nine unit
case classes, campaign-policy participation, structured lifecycle events,
entry/replay/DSR bindings, sealed canonical obligation projections, and
stage-specific tiers; fixed-point signedness/scale/rounding/overflow semantics;
the four non-confusable timing-evidence kinds; whole-interval clock alignment;
fault-fidelity/safe-state and regulatory no-claim boundaries; directed target
refinement; theorem/proposition and finite-target/exhaustiveness successor
ratchets; G3 assumption/oracle/tolerance/holdout/partition/policy mutations;
G4 chunked assembly equivalence; G5 input-order invariance; and targeted versus
global amendment invalidation.

`tests/i06.rs`: exact I06 7/2/3 lattice, single maximal refutation lane,
17-fixture/7-held-out corpus, one governed industrial waiver, and once-only
seven-leaf claim map; all nine unit classes, campaign-policy participation,
structured lifecycle/events/replay/DSR bindings and sealed canonical
obligations; non-confusable property kinds, domains, units, frames, methods,
censoring and calibration; authenticity-versus-physical-truth, uncertainty-
ownership, pseudoreplication, three-valued substitution, declared-graph and
NoAlarm no-claim boundaries; disjoint single-consumer Core/Max ranges; causal/
impact theorem and finite-decision global-optimality successor ratchets; G3
hypothesis/oracle/tolerance/holdout/partition/policy mutations; G4 chunked
assembly equivalence; G5 order invariance; and targeted/global amendment
invalidation.

`tests/i07.rs`: exact I07 8/4/3 lattice, 27-fixture/11-held-out corpus,
one governed-real-corpus waiver, and once-only ten-leaf claim map; policy locks
requiring a complete profile-definition schema/table plus exact locks for the
key/entry/registry/required-universe/assignment-root formulas, versioned
framing, bootstrap caps, schema cardinality, domain equality, required successor
decoder/KAT batteries, and no runtime profile choice; native occurrence,
geometry, PMI,
material, harness, kinematic, migration, semantic-loss, and industrial-matrix
no-claim boundaries; exact split/merge provenance incidence; public versus
governed evidence axes; the closed execution/claim/completeness ontology and
cross-axis validity rules; certified outward geometry-score arithmetic;
single-consumer holdouts, bounded corpus/tool cancellation, and logical-event
versus telemetry determinism; policy locks for five-trigger/cause-arbitration
and inclusive-cap arithmetic; success-pair versus explicit absent-artifact
publication; exact NUL-terminated twenty-field governed transaction grammar,
one-row/single-output AP242 discharge schema matrix, closed 34-field
role-derived envelope schema and singleton schema-set root, governance/
candidate/intent root equality, authorization membership proof, cap/depth, and
permissive-schema/role/stage/scope/add/drop mutations; required receipt-schema
and derived-new-head successor tests; role-addressed anti-swap bindings;
primary-request observer/child seals with
empty-set and calibrated-clock rules; bounded lossless observability with
reservation settlement and global EDF priority-service admission;
sheaf/stack/effective-quotient theorem ratchets; G3 semantic/profile/oracle/
partition mutations; G4 chunked assembly and cancellation contracts; G5 order
invariance; and targeted/global amendment invalidation.

`tests/i10.rs`: exact I10 6/2/3 lattice, single maximal refutation lane,
14-fixture/6-held-out corpus, one instrumented-lab waiver, and once-only
six-leaf claim map; all nine unit classes, campaign-policy participation,
structured lifecycle/events/replay/DSR bindings and sealed canonical
obligations; gauge-confounding semantics (numeric-FIM no-mint, seeded
hardening/Prony/product confoundings, invariance-witness replay,
sloppiness-versus-structural wording separation); enclosure-checked adjoint,
feasibility-precedes-optimality, identifiable-combination-only gain, and
Undecided-never-equivalence no-claim boundaries; disjoint single-consumer
Core/Max ranges; quotient-theorem and finite-grammar design-optimality
successor ratchets; G3 hypothesis/oracle/tolerance/holdout/partition/policy
mutations; G4 chunked assembly equivalence; G5 order invariance; and
targeted/global amendment invalidation.

`tests/i13.rs`: exact I13 6/5/6 lattice, single maximal refutation lane,
28-fixture/12-held-out corpus, one governed industrial-machine waiver, and
once-only fourteen-leaf claim map; exact BLAKE3/Philox key-counter-lane
algebra and disjoint semantic streams; independently owned manufacturing,
route-realization, robust-QoI, and multiphysics-closure decks/heldouts;
relative winding/cohomology, differential-character, topology-event,
energy/tangent, adjoint, uncertainty, manufacturing, insulation, route, and
machine no-claim boundaries; theorem/proposition and finite-grammar
falsification ratchets; G3 topology/units/representation/partition mutations;
G4 cancellation/checkpoint/drain contracts; G5 order/thread/shard invariance
requirements; closed two-kind governance namespace; policy/table locks for
canonical governed-group framing/root/KAT/caps and transition/access genesis;
the exact twenty-five-field two-protocol-by-five-stage value matrix;
Prepared/terminal access correlation with the explicit terminal-effect option;
exhaustive nineteen-field transaction tag/type table; distinct governance/
science industrial attempts; sealed cancellation replay and four-frontier
reconciliation; deterministic trigger-arbitration and exact drain/finalize/
total arithmetic policies; full poll/enqueue/queue/service/calibration
logging-latency composition; corrected enumerable single-governor and dual-
service reservation/settlement state machines; closed semantic inventory/
evaluator/KAT/catalog authority with repeated-source occurrence aggregation;
exact holdout
commitment and realized Merkle algebra; mirrored breach, structural-event,
canonical-payload, and receipt-index durability; typed operation/batch authority;
bounded append metadata and phase history; dual-journal append-then-publish
authority; branch-specific retention and settlement evidence; exact route/
failure and requirement-role universes; V4 finite index/decision/response/WAL/
single-CAS capacity authority; machine-checked fresh-V2/V3/V4 prerequisite-DAG
acyclicity plus tagged-
sum, required-edge, and required-path locks; bounded priority reservations; and
local versus global amendment invalidation.

The I13 proof program is cumulative, not satisfied by prose-substring tests.
Static locks must independently parse every still-authoritative V2/V3
append/terminal clause and every terminal V4 capacity clause: header/NUL/domain,
enum, field order, exact byte count, checked size equation, cap, residual, root
projection, route/stage/failure set, requirement-role set, state transition,
tagged-sum discriminator, and authority-DAG edge. The nonconflicting append
locks still recompute the exact 21/315/64 event-registry counts, 99 route/failure
pairs, 4096 append limit, 261447712 durable total, 6987744 nonborrowable slack,
the 5324902/2801937/1572964/659555/126402 component maxima and exact 10485760
sum, `84+payload` operation authority, 257-byte lane bundle, 376-byte phase
lease, 153-byte plan, 633-byte slot ledger, 258-byte logical-slot intent,
phase-authority history, branch-specific retention, and every route and
terminal bound.

The V4 locks additionally recompute, rather than merely substring-match:
the exact fifteen-row/2609-byte authority roster; `77+210n` operation matrix;
4096-by-393216 response/WAL and 24576-by-32768 decision baselines and their
2415919104-byte per-copy sum; the 4096-by-131072 exception-payload bank,
4096-by-32768 authenticated evidence sidecar bank, and complete
3087007744-byte per-copy substrate; 180/228/272-byte index, decision, and copy-
coordinator heads; reserved-subject partition and six-phase budget;
328-plus-payload decision entries; 395-byte local decisions and 723/756-byte
local/aggregate entries; 335-byte coordinator-binding append receipts;
426/636-byte allocation request/receipt; 98-byte map members and 8330-byte
proofs; clean 30573/30788/31422-byte and maximally reconciled
40179/40394/41028-byte bundle/pair/authority values; 40666-byte maximum pair
response entry; the 163-byte writer row and `136+171R` registry; 195-byte
dense issuance member; 343/217-byte capability subject/service heads;
515-byte subject-vector proof; 437-byte signed fence intent; fixed 6571-byte
all-sixteen-row issuance proof; 72-byte no-issuance proof; `329+120C`/1169-byte
expiry quorum; 581-byte nonrecursive fence commit; and exact
2790/9289/10466-byte NoIssued/Revoked/Expired dispositions. They lock the
`317+166W`/2973-byte all-and-only writer roster; `192+43N`/1568-byte signature
object set; `195+77K`/811-byte canonical signature set; 67-byte layer binding;
489/564-byte four/five-layer set; and exact 5949/7043-byte causal signature
material. They also lock the 477-byte writer-quiescence receipt,
`230+485Q`/7990-byte set, 187/212/841/495/422/691/2492-byte storage member,
head, proof, intent, reconstructible stable-extent projection, nonrecursive
commit, and disposition; 601-byte availability bundle; and independently
replayable 26564-byte availability-preimage cap.

The terminal exception-publication locks recompute the 178-byte coordinator
core and 283-byte contextual root preimage; distinct depth-12 slot and
PreEffect-receipt trees; 771-byte dual vacancy witness; compact 640/828-byte
inner-root proofs opened by the exact old core; generic phase-tagged 384-byte
conditional receipt without a committed-head backedge; and 199-byte global
exception-chain step. Candidate, envelope, durable-absence, unavailable-record,
Empty/Present exception-artifact, direct-snapshot artifact, and worst slot use
are respectively at most 27123, 37218, 38266, 28292, 29526/29714, 835, and
93826 bytes. The evidence store uses exact 387-byte charge receipts and a
32144-byte branch-independent sidecar overbound. The conservative admitted
Reserved maxima are exactly 198055/198270/198904/198542 bytes for
bundle/pair/authority/response wrapper with at most 74 bundle items, while the
raw/framed staging-page row maxima are 199063/199071. These are checked caps
from independent per-object maxima, not a claim that every component maximum
is jointly attainable; and the 531-byte profile;
428-byte outcomes and 132/174/8500-byte orphan member/head/proof; 214/272/257
response and 185/128/260 WAL head/entry-overhead/receipt values; 221-byte
transaction coordinator, 346-byte chain step with old/new transaction counts,
`140+primary`, `370+138*components`, 25100-byte update-proof cap,
`676+proof`, 359-byte receipt, and
`3299+primary+138*components+proof` durability formula with every
stated per-operation maximum. They also lock `137+94n`/`136+94n` evidence
manifests; `75+141n` bindings; the derived `A<=30`, `M<=29`, item-count<=63,
`61+299n` bootstrap and 18898/19321/19454 maxima; 65536-byte mutation archive;
114657-byte closure artifact; 1298/1429/2489-byte normal control tail and
63047 slack; 176-byte/65536-slot/131072-byte closed-evidence map, 114822-byte
member and 131669-byte expiry-proof bundle; `77+163n` staging inventory;
116/100-byte member arms, 215-byte staging head, 16716-byte proof bundle,
`106+L+21261n` artifact, 460/354-byte PreparedOnly operations, 20653/20547-byte
durability maxima, 225-byte cutoff, at-most-199071-byte framed page rows, exact page
table/cross-product/index/decision-prefix gates, and
`Smax+21319*Nmax+Lmax+20547+66505<=67042864`. The generated protocol matrix must
derive every exact response/WAL prefix and equality/cap-plus-one twin; no V3
aggregate subtotal may be substituted. Two independently authored encoder/
decoder implementations must reproduce the same positive KAT bytes and roots
and reject truncation, trailing bytes, unknown tags, noncanonical order,
duplicate/drop/add/swap, overflow, cross-attempt/epoch/lane/domain substitution,
old-version relabeling, and unavailable hash preimages.

Focused unit tests must retain every append-phase, route, breach, terminal, and
scientific-retention transition test, while replacing V3-as-fresh capacity tests
with the V4 finite protocol. They must cover empty/first/equality/max/max-plus-
one index, map, decision, response, WAL, orphan, active, closed-evidence, and
staging states; all legal and illegal dispositions; exact-key retry and
conflict; stale generation, ABA, arithmetic overflow, two-absence refusal,
one-sided repair, per-subject six-phase/reserved-slot accounting, phase-3
aggregate recovery before response-log creation, acyclic phase-4/5 orphan
reconciliation, phase-6 terminal outcome, and every immutable Prepare terminal;
all component/binding/manifest omissions, additions, duplicates, permutations,
and PreparedOnly substitutions; and crash/restart before and after each local
decision append, response append, intent append, WAL append, coordinator CAS,
and pure-receipt reconstruction. Bootstrap tests must enumerate K=1..15 and
both known-no-effect and pending partial outcomes, prove the 63-item bound, and
prove bootstrap cannot enter normal H2. Closure/expiry tests must replay both
branches from the materialized artifact after active-log expiry and prove that
active deletion cannot precede executable closed-member insertion.

Capacity-index unit/property tests must use independently authored encoders to
lock every capability issuance state, subject/service head, dense transition,
nonissuance/expiry proof, fence intent/commit/disposition, roster, causal
signature object/set/layer/material, quiescence receipt/set, storage
member/head/proof/intent/commit/disposition, availability bundle, charge
receipt, coordinator core, dual vacancy witness, compact slot proof, generic
phase receipt, exception-chain step, direct snapshot, and composite artifact.
They independently recompute every exact size, frame, equality, conservative
maximum, and maximum-plus-one. Mutation campaigns cover header/NUL/domain,
canonical-zero rules, signer/object/layer order, writer identity/ordinal/role,
issuance state and epoch, clock failure-domain independence, signature-set
resolution/revocation, same-layer backedges, service-head recursion, stable
subject versus proposed held ordinal, both Merkle sibling paths, old-core
aggregate opening, receipt-vector/global-chain recurrence, storage prefix,
payload/proof/receipt joins, sidecar charge/reclaim, and staged-before-CAS
preconditions.

Resolver properties distinguish direct from transitive fields, permit harmless
current-head request refreshes, require nested interval containment, and reject
capability reuse across copy, role, subject, ordinal, attempt, or entry. The
`882*883` matrix test covers byte-identical terminal phase-3 repair for Reserved
and RefusedNoEffect plus the distinct nonterminal repair arm for genuine
PartialPending/ExpiredUnknown; enumerates every legal
`(required_local_mask,repaired_mask)` pair including `RepairedBoth`; and rejects
wrong-copy, repaired-bit-without-artifact, direct-bit-without-snapshot,
direct-bit-with-exception, double recovery, and partially repaired two-copy
final histories. A graph lock requires explicit phase-3 and phase-6 direct
snapshot parents as well as acyclicity and tagged-sum reachability.

Control tests must generate every producer row all-and-only; prove one physical
owner for every extent, every reference's owner membership, and zero-
disposition nonphysicality; generate the finite staging-page cross-product,
materialize every page allocation authority set, prove all response/WAL/index/
decision caps, and replay every staging generation from materialized durability
sets; reject insertion after the seal durability set; reproduce the three
special E extents and exact paged admission equation; walk the V4 normal H2
durability and H3 order; and reject every omitted preimage, arm splice, early
release, layout gap/overlap, duplicate charge, equality-plus-one, or backward
hash.

The no-mock I13 end-to-end scripts must eventually exercise real persisted
implementations of the governor, both lane stores, both publication journals,
independent verifier/logical-slot allocator, all fifteen V4 logical authorities,
both physical copies of every index/decision/response/WAL allocation, the
single-CAS coordinators, generic cleanup and ordinary-evidence stores, two
physical retention stores, closed-evidence map, staging service, and control
settlement ledger. They
must cover normal selection/finalization, preactivation Committed/Refused/
ExpiredUnknown cleanup, preselection and postselection breach independently and
together, primary exact replay, mirrored publication replay, emergency
failover, durable-prefix recovery, exhaustive catastrophe, cancellation at
every phase, and process crash/restart before and after every prepare, CAS,
journal, verifier, metadata, transfer, release, settlement, and closure point.
They must additionally crash around every local allocation outcome, terminal
and nonterminal phase-3 repair, both-copy phase-6 repair, capability-intent
signature, capability-service CAS, writer drain, quiescence-set signature,
storage-intent signature, storage-service CAS, final availability signature,
sidecar charge, dual-vacancy stage, conditional receipt, receipt-vector/global-
chain derivation, durable-artifact staging, exception-vector update, sticky-
lease recovery/abort, orphan reconciliation,
Prepare/Intent/WAL/coordinator transition,
bootstrap closure, staging generation/seal, normal H2/H3 transition, closed-map
expiry move, and store late-commit reconciliation.
Each run emits deterministic structured logs containing attempt/epoch, logical
generation, state transition, predecessor/successor roots, reservation and live
bytes, capability/idempotency identity, route/stage/failure, shard/store,
checked interval, disposition, capability/storage generations, signature layer,
subject/held ordinals, slot/receipt/global roots, sidecar charge, and artifact/
receipt digests; secrets and protected payloads remain content-addressed rather
than printed. On failure the script preserves a replay bundle, exact failing
seed/step, deterministic mutation-sequence position, and all independently
rehashable preimages needed to diagnose the first divergence.

Gauntlet evidence is explicit: G0 covers encode/decode inverses, set/map/root
laws, exact sums, tagged unions, state-machine safety, and checked arithmetic;
G1 covers manufactured append/retention histories with independently recomputed
roots and conservation of bytes, counts, generations, and capabilities; G2
replays canonical normal, breach, failover, recovery, and catastrophe vectors;
G3 mutates every identity, edge, role, route, preimage, order, cap, and branch;
G4 runs cancellation storms, injected crashes, unavailable services, torn
publication attempts, quarantine recovery, retention loss, leak/deadlock
checks, and exact restart equivalence; G5 permutes worker count, thread, shard,
arrival, replay, and serialization order on the same ISA and requires identical
canonical authority, artifact roots, receipts, disposition, and logs after
declared nondeterministic telemetry is removed. Authority remains disabled until
the intended tier is green with replayable evidence.

`tests/i14.rs`: exact I14 3/12/13 lattice, one maximal refutation lane,
53-fixture/20-held-out corpus, seven independently scoped governed waivers,
and 28 one-claim leaves split 13 Core/15 Max; exact seed derivation and
EmConventionCard sign/phase/RMS/outgoing-wave closure; hard-versus-soft
componentwise acceptance arithmetic with certified outward enclosures;
separate native/AP242, RLGC/MTL, source/victim, UQ/population, BEM/FMM,
laboratory/bearing, standards, safety, and blind-mitigation authority;
atomic governed stage/discharge semantics and exact waiver blast radii;
passive-composition, hypercohomology-obstruction, cover-naturality, KYP/sheaf,
fidelity, and countermodel theorem ratchets including AuthorityContradiction;
typed exhaustive terminal-status normalization with retained raw tuples,
explicit normalization actions, well-formed witnesses for every exit, and
pinned-table identity; count-defined logical tiles separated from wall-time
watchdog telemetry; bounded typed terminal-cause selection with canonical
malformed refusals, scope/observer/causal validation, every terminal and
nonterminal decision path, pending-request deferral, and cause precedence;
schema-authoritative V2 Core/Max/theorem-Max card admission, typed count/memory/
tile/resource derivation, exact deadline/cap enforcement, genesis/contiguous
request-inclusive first-terminal prefix proofs, bounded arbitration work,
temporal-before-cause precedence, opaque nonterminal frontiers, mandatory
drain/finalize, watchdog, tile-poll, child, and external-heartbeat coverage,
failure-preserving mandatory typed spawn-frontier evidence, typed receipt-bound
infrastructure onset, derived timeout onset with exact/cap-plus-one boundary
tests, exhaustive trigger/source wire KATs, immutable lifecycle views,
verification-receipt-bound telemetry-only late-event tails whose authentication
is deferred to HELM/ledger, explicit verifier-authority no-claim boundaries, and
canonical-versus-telemetry timing metamorphisms;
typed exhaustive artifact-retention rules; actor-scoped physical-validation,
licensed-standards, blind-mitigation, and governed-population authority
lifecycles with exact-once event-vocabulary membership, pinned causal-order
contracts, and cancellation-safe atomic stages; replay/DSR bindings,
single-consumer holdouts, G3 mutation relations, canonical-result versus
noncanonical-telemetry G5 boundaries, and targeted/global amendment
invalidation.

`tests/rl.rs`: exact RL 7/2/3 portfolio lattice, single maximal refutation
lane, 12-fixture/5-held-out corpus, one governed-physical-loop waiver, and
once-only five-leaf claim map; all nine unit classes, campaign-policy
participation, lifecycle/events/replay/DSR bindings and sealed canonical
obligations; seam semantics (content-address identity spine, side-channel
custody scope, exactly-once ownership with exact recomposition, closed
typed-refusal gate taxonomy, weakest-wins wording, target-exact timing with
the no-WCET-promotion rule, instance-digest amendment propagation);
disjoint single-consumer Core/Max ranges with loop seeds isolated from
instance-manifest domains; composition-theorem and physical-protocol
successor ratchets; G3 hypothesis/oracle/tolerance/holdout/partition/policy
mutations; G4 chunked assembly equivalence; G5 order invariance; and
targeted/global amendment invalidation.

`tests/i03.rs`: exact I03 8/2/6 claim lattice and once-only 4-Core/
4-Max leaf mapping; exact nine unit-case classes; content-bound campaign
policy, FailureBundle, and independent-adjudication obligations; exact
22-fixture/7-heldout stage-local partitions with exact deterministic fixture
aliases, inclusive development/Core/Max Philox ranges, collision-free derived
keys, and one public governed statistical-holdout protocol whose realized
nonce/lots/features/predictions/labels remain non-public until adjudication and
whose exact IID authority is conditional on candidate-first semantics, an
authenticated peer-withheld 1024-block join-DAG transcript, an audited honest
uniform component, exact dyadic sampling, candidate-input isolation, and
order/external-id invariance; separately pinned theorem-axiom and formal-
projection no-authority gates; target-only M0 grammar formalization gate;
orthogonal execution/predicate/claim/completeness/integrity/support/observable/
promotion axes, including `DeclaredDivergent` only as an observable disposition;
exact per-leaf deck/event/entry/replay/DSR maps, and one-consumer holdouts;
explicit common cancellation/checkpoint lifecycle events, two declared ISA
families with bitwise comparison confined to identical ISA fingerprints, and
I03-specific chunked clone-boundary assembly identity; per-leaf G4 fault,
no-partial-publication, drain/resume, and durable-retention authority; per-leaf
G5 exact-output authority; G5 order stability; G3 hypothesis/oracle/tolerance/
holdout/policy mutations; and targeted versus global amendment invalidation.

## No-claim boundaries

- A frozen manifest asserts NOTHING about implementation correctness:
  preregistration is not proof, and no evidence color, receipt, or
  promotion authority is minted here.
- The named `scripts/e2e/leapfrog/*.sh` entry points and DSR lanes are
  preregistered locations; their existence and behavior are verified at
  campaign time, not freeze time.
- Authored fixture text is the immutable generator/theorem contract, not
  proof that an executable generator or dataset already exists. Campaign
  receipts must additionally bind implementation/toolchain identity and
  the exact generated artifact bytes before evidence can promote.
- Named independent BEM, Maxwell, interval, and Lean oracle routes are
  preregistered identities and separation contracts. Freeze does not prove
  that those implementations exist, are correct, or satisfy the claimed
  independence; campaign receipts and certifier-adversary lanes must do so.
- I14's checked-in terminal/card/prefix/telemetry/watchdog KAT constants are
  known answers exercised through the current local implementation, not
  independent reproductions. Its governed transaction authority is also
  authored contract rather than implemented authority: this crate does not yet
  provide two independent encoders for the exact 23-field `P`, the closed
  three-role schema set, role-addressed receipt/output roots, derived successor
  head, or 361-byte commit receipt; production bounded schema/receipt decoders;
  authenticated schema/receipt registries; or the atomic idempotent fenced CAS
  ledger. Those exact KAT, cap/cap-plus-one, schema/role/output-swap,
  stale/alternative-head, generation-overflow, ABA, replay, crash-before/
  after-commit, and partial-state gates plus HELM/ledger authentication must
  pass before any I14 discharge or promotion authority. Prose locks and local
  known-answer equality alone grant none.
- The manifest digest golden constant and I07/I13 stream-key/Philox
  known-answer vectors are deliberately NOT frozen in this crate yet: per
  `docs/GOLDEN_POLICY.md` those pins require committed-tree, two-mode
  reproduction scheduled with the batch-verify lane. The I07/I13 policies
  require concrete successors to freeze independently reproduced development
  and heldout-format KATs before evidence generation; this draft crate itself
  freezes no computed I07/I13 Philox output words. Current tests lock seed
  syntax, ranges, derivation domain/endianness prose, KAT admission
  requirements, and collision structure, not a computed Philox known-answer
  value. This version likewise freezes the exact I07/I13 transaction-intent
  grammar and mutation requirements but no computed intent KAT: the governed
  successor/discharge lane must publish independently reproduced vectors before
  authority use. Cross-ISA digest stability is likewise expected but unproven
  until the two-host campaign runs.
- I07/I13 profile, governed-group, governed output/protocol-schema authority,
  transaction, authority-head, lifecycle, trigger arbitration, source-class and
  candidate/post-selection capacity-breach settlement, sparse-Merkle candidate/
  witness/event frontiers, causal-parent schema/proofs and barrier parent lists,
  semantic inventory/evaluator/KAT/catalog closure, receipt-index batches,
  lane-epoch failover, deadline ticket/resolution, mirrored breach/structural/
  payload/index/operation/delta artifacts, authenticated append-phase history,
  reconstructible append metadata, replicated append-then-publish durability,
  closed replay-route/failure authority, shard-complete scientific retention,
  persistence-recovery heads, typed quiescence/finalization/catastrophe
  evidence, holdout commitment/Merkle, and latency clauses are authored schema
  contracts in this freeze. I13 retains final V2/V3 append, route, terminal, and
  scientific-retention authority only where nonconflicting. Its finite capacity
  index/decision journals, response logs, WALs, transaction coordinator,
  bootstrap, closure evidence, control staging, expiry, and downstream binding
  are the terminal V4 schema authority and remain `[M]`.

  The current Rust surface locks vocabulary, field tables, dependency edges,
  exact sizes, formulas, caps, and mandatory proof suites; it is not runtime
  proof that those protocols exist or work. In particular this crate does not
  yet supply the two independent bounded V4 encoders/decoders and computed KATs;
  generated fifteen-authority roster, operation, producer, component, and
  all-and-only binding matrices; provisioned index/decision/response/WAL
  extents; local copy coordinators and decision reconciliation; global
  transaction coordinator maps and sole atomic CAS; prefix reconstruction,
  repair, idempotency, and crash recovery; materialized pre/postcommit manifests
  and durability projections; authenticated orphan, active, and closed-evidence
  maps; bounded bootstrap derivation; self-contained closure archive; atomic
  expiry move; staging member/owner maps, generation artifact, seal, cutoff and
  physical layout; or HELM/ledger authentication. Nor does it instantiate the
  still-required production catalog/evaluator/generator, canonical append
  stores, append-phase lease/history service, breach mirrors, lane services,
  publication journals/verifier, shard/retention stores, route/failure/failover/
  recovery machinery, or checked runtime arithmetic.

  Historical V1 artifacts and superseded V2/V3 capacity machinery are replay
  specifications only except where a final nonconflicting clause explicitly
  wraps and rehashes their exact envelope. A legacy capacity receipt, a
  PreparedOnly primary, a raw digest, a static substring, or favorable
  arithmetic cannot silently acquire V4 committed authority. Before authority
  use, concrete successors must supply every specified independent
  implementation; exact-byte/domain KAT; exact-cap/cap-plus-one and overflow
  test; adversarial semantic/type/role/component/binding mutation; empty/first/
  max finite-state test; append genesis/phase/reconstruction proof; two-lane,
  two-copy, two-response-log, and two-WAL crash matrix; repair and orphan
  reconciliation proof; every Prepare terminal; bootstrap and normal-settlement
  arm isolation; closure replay after active-log expiry; atomic closed-map
  expiry; staging ownership/generation/seal proof; all 99 route pairs; every
  settlement arm; shard/replica failure; early-release/retention-loss test; G5
  permutation; and authenticated runtime receipt with structured deterministic
  logs. Passing static locks alone grants no execution, durability,
  determinism, discharge, closure, expiry, theorem, or promotion authority.
- Waivers record discipline; they do not verify that owners discharge
  them by expiry — that policing belongs to the governance lane.
- Amendment invalidation names predecessor claim/obligation authority;
  the downstream ledger/governance layer must actually revoke or refuse
  stale receipts carrying those identities. Unchanged component evidence
  may be rebound across a version-only or targeted successor only through
  authenticated amendment lineage and byte-identical component digests;
  whole-manifest digest inequality alone neither revokes everything nor
  authorizes blind reuse.
- An invalidated obligation leaf is a mapping-bound authority whose receipt
  must be reissued or explicitly rebound. `AmendmentRecord` is exact at
  claim/leaf authority granularity; it is not by itself a numerical-kernel
  rerun schedule, and authenticated lineage may reuse byte-identical execution
  payloads where the surviving claim authority permits it.
- Every authored held-out fixture SPEC is frozen here. Deterministic holdouts
  also freeze public replay ranges and therefore carry no statistical
  untouched/IID authority. The discharge max holdout deliberately has no
  public seed. Candidate/model/toolchain and every campaign semantic freeze
  first; only then do three named custodians enter a peer-withheld,
  all-three-fixed commitment phase over exact lot-major 1024x256-bit vectors.
  Strict RFC8032-authenticated transcripts, an audited at-least-one-honest IID
  uniform vector independent of the already-fixed candidate and adversarial
  masks, exact coordinate-wise XOR, 257-bit half-open cumulative intervals for
  the pinned dyadic atom law, a later 256-bit order-challenge beacon that never
  enters sampling, and complete
  candidate-input isolation are all authority-bearing. Exactness does not rely
  on computational hiding or a short-seed PRG. Its IID receipt, external-id
  relabeling/order invariance, access control, one-shot adjudication, and later
  raw-draw reveal must be established by campaign integrity receipts. This
  crate preregisters that protocol but cannot enforce the external governance.
- I03 synthetic discharge holdouts can support calibration evidence only;
  physical/industrial authority remains blocked until the waived,
  independently governed experimental pack is admitted. I03's exact
  topology/force theorems apply only under the complete frozen
  variational or regularized-event premises; generic cohomology
  preservation, bare chain equivalence, refinement, or topology change
  receives no free invariance claim.
- I03 theorem target cards in manifest version 1 are not machine propositions.
  They mint no theorem color until a pre-proof successor freezes canonical
  proposition AST and definition bytes/digests, a total runtime-premise map,
  and a deterministic AST-to-Lean translation with structural round-trip
  checks. The exact axiom policy admits only `propext`, `Quot.sound`, and
  `Classical.choice`; its digest and the complete transitive axiom closure are
  receipt-bound, while `sorryAx`, custom postulates, and unsafe/native-oracle
  shortcuts fail integrity.
- I03's version-1 `M0` prose similarly freezes the ambitious
  16x16x16x16/N=65536 target but has no exhaustive authority. A pre-search
  successor must freeze the full decorated-record grammar, canonical
  encodings/domains, validity/stratum/tag/parameter semantics, explicitly
  encoded event primitives, total enumeration and exclusion order,
  rank/unrank/sharding algorithms, source digests, independent decoder and
  bijection proofs, cost preflight, and Merkle completeness root.
- I05 timing evidence never promotes a `MeasuredSampleMaximum` into a WCET
  upper bound. Static and compositional bounds carry authority only for their
  exact binary, target, clock/power state, task context, microarchitectural and
  interference assumptions. HIL clock containment, injected safe-state
  behavior, and bounded trace agreement are likewise not regulatory safety,
  physical-plant validation, unbounded stability, or authority on another
  target/rig profile.
- I05 version-1 theorem and finite-target cards are ambitious targets, not
  machine propositions or executable exhaustive grammars. A pre-proof
  successor must freeze canonical proposition/definition/decoded-target ASTs,
  total runtime premises, translation and axiom closure. A pre-search successor
  must freeze complete instruction/microstate/input grammar, initial/valid/
  transition/cost semantics, unsupported-state closure, enumeration or sound
  symbolic coverage, quotient obligations, independent decoding, preflight and
  completeness root. Even a model-exact maximum becomes silicon WCET only
  after independent model-fidelity qualification on the exact hardware.
- I06 cryptographic and custody validity authenticates records and transitions,
  not physical composition, label-to-specimen correspondence, supplier honesty,
  method execution, calibration, specification conformance, qualification, or
  legal authority. Synthetic posterior coverage and drift-certifier tests do
  not establish a real supplier population; `NoAlarm` never proves no drift.
  ContextOfUse `Compatible` is neither universal interchangeability nor an
  automatic procurement/manufacturing/safety approval.
- I06 exact impact closure is initially only relative to the frozen declared
  graph. Version-1 causal/impact theorem prose and finite decision grammar mint
  no completeness, causal, exhaustive, or optimality authority. Pre-proof
  successors must freeze machine graph/SCM/authority/proposition/definition
  ASTs, complete adapters/open-world boundaries, runtime premises, translation,
  axiom closure, nonvacuity and kernel replay. Pre-search successors must freeze
  the complete finite action/scenario/test grammar, feasibility/interval utility,
  ambiguity/risk semantics, validity/exclusion, enumeration or verified bounds,
  rank/unrank/sharding, independent decoding, preflight and completeness root.
  Global optimality requires coincident checker-verified lower and feasible
  upper bounds; a positive gap remains honest non-optimal/Unknown authority.
- I10 structural verdicts hold for exact model structure and noise-free
  observables and never certify practical recoverability under noise or
  discrepancy; synthetic profile coverage and held-out gain are not real-
  laboratory authority, `Feasible` is not fabrication approval, `Undecided`
  never proves rival-law equivalence, and the quotient theorem and zero-gap
  design optimality are version-1 prose targets: pre-proof successors must
  freeze law/observable/group ASTs, premises, translation, axiom closure,
  nonvacuity and kernel replay; pre-candidate successors must freeze the
  finite design grammar, bounds, rank/unrank/sharding, independent decoding,
  preflight and completeness root. All I10 fixtures are synthetic
  coupon-shaped data; real instrumented-lab campaigns enter only through the
  waived governed pack.
- RL portfolio authority is seam bookkeeping, never stage truth: identity
  continuity, custody, ownership conservation, replay, reproof, and gating
  do not accredit laboratories, suppliers, devices, targets, or models
  outside the declared ContextOfUse, and a green deployment gate is not
  physical validity or safety authority. The end-to-end composition theorem
  and physical receipt parity are version-1 prose targets behind pre-proof
  and waiver-discharge successors; a composed guarantee is conditional on
  every stage premise and validates none of them. All RL fixtures are
  synthetic loop worlds; real physical campaigns enter only through the
  waived governed pack under the rl-physical-protocol-card.
- The adjacent-version `EvidenceKindChanged` guard prevents an immediate
  claim/leaf kind swap. `FrozenManifest` carries no lineage-wide tombstone set;
  a ledger spanning nonadjacent versions must key authority by typed kind plus
  manifest/version digest (or add an authenticated tombstone schema) rather
  than treating a raw string id as globally unique forever.
- The clone-boundary assembly tests do not encode or decode durable
  checkpoints, restart a process, detect corruption, or exercise runtime
  cancellation. Those remain executable G4 campaign obligations.

## VerificationManifest v1 identity core (bead i94v.7.1.1)

`v1` is the stable-identity and source-authority layer: `ClaimId` names a
conceptual lineage while `ClaimRevisionId` content-addresses one exact
statement (kind, quantifiers, units/conventions, hypotheses, domain,
code/contract surface, no-claim boundary, supersession pointer) — distinct
revisions cannot collide, identical content is idempotent, supersession
appends and never mutates. `CaseId`/`JourneyId` are the stable
case/journey identities. `SourceAuthority` is a total lattice
(GeneratedArtifact < TestSource < Contract < BeadObligation <
FrozenSnapshot); conflicts resolve upward by re-pinning, and
equal-authority conflicts refuse with ranked fixes. Typed
`ClaimRelationReceipt`s (implication, refinement, restriction,
counterexample, certified equivalence) carry direction, checker/TCB,
quantifier variance, and policy version; promotion never transfers along
counterexample or quantifier-strengthening edges; directed cycles refuse
unless certified-equivalent, in which case the SCC canonicalizes to its
smallest member without erasing members. `NormalizedGraph` digests are
input-order invariant and the human/JSON/ledger renderings are tested
semantic projections of that one digest. Migration is additive or
breaking-with-mandatory-lossy-report; the 22-row `MANIFEST_RECORD_FIELDS`
registry declares units, cardinality, authority, default visibility, and
migration semantics per field, in data.

### No-claim boundaries (v1 identity core)

- The manifest is metadata and obligation authority, not proof; nothing
  here adjudicates a scientific claim.
- The frozen inventory compiler is V.1.2 scope; the lint battery V.1.3;
  ledger persistence fs-obs/fs-ledger scope.
- The field registry SPECIFIES the record; the typed full-record wire
  codec lands with the V.1.2 compiler against this registry.
- Relation soundness is structural (orientation, variance, cycles,
  contradiction); checker/TCB strings are recorded identities, not
  re-verified proofs.
