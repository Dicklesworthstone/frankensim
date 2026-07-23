# CONTRACT: fs-govern

The addendum's governance as machine-readable data: the design principles
(P1–P8), the governance rules, the nineteen proposals (with kill metrics +
owning beads), the original Part V runtime risks (R1–R10), the distinct
expansion-program risks (PR-001–PR-012) — each with a CI-gateable completeness
surface — plus the EXECUTABLE one-bet-per-lane admission state machine
(`lanes` module, bead rjoq.6) and the generated requirement-to-evidence
registry (`traceability` module, bead
`frankensim-ext-traceability-ledger-ziso`).

## Purpose and layer

Layer UTIL. Pure data + audit, with `fs-blake3` as its only dependency for
canonical content identities. Encodes the doctrine, proposals, original ten
runtime risks, and twelve expansion-program risks, audits that nothing
survives unmeasured (design principle P8 / Governance Rule 2), and enforces
the one-active-unproven-mechanism-per-independently-falsifiable-proof-lane
rule as an atomic, replayable admission ledger. It also owns the pure typed
source registry and fail-closed generator for the extension charter's B1–B14
and RQ-* requirement-to-evidence catalog, plus the Phase 0B-C descriptive
support/threat graph and deterministic allocation-candidate planner.

## Phase 0B-A evidence-contract algebra (`evidence_contract` module)

`AUTHORITY_ALGEBRA_VERSION = 2` is the pure common contract consumed by later
schema, graph, checker, ledger, and runtime layers. It owns canonical semantic
identity and authority composition; it performs no filesystem/network I/O,
signature verification, durable persistence, or graph traversal. Constructors,
canonical sets, JSON rendering, and diagnostics use bounded heap allocation.
Those adapters must preserve this algebra rather than reconstructing it from
booleans or prose. Every public evidence, state, edge, adjudication, tombstone,
and checker-decision constructor creates descriptive candidate data only. The
opaque grant, authenticated decision, current-head, and runtime-admission types
have no public minting path in Phase 0B-A; Phase 0B-B must authenticate durable
receipts before it can use their private boundary.

### Exact objects and identity bindings

The code-owned `AUTHORITY_CATALOG_ROWS` table is the source of
`authority_catalog_json()`. The contract drift suite requires every row below
and its domain to remain represented here. The catalog serialization schema,
the authority algebra used to interpret a row, and the row's own identity
schema are independent version axes. Catalog serialization v3 is an intentional
pre-freeze source-format correction: the row-level identity-schema field was
added while the outer marker accidentally remained v2, leaving that label
ambiguous across incompatible row shapes. The
`validate_authority_catalog_schema_version` gate accepts exactly v3; it refuses
retired catalog v2 and unknown/future versions instead of silently
reinterpreting them. In particular, the validated proof-lane identity remains
v1 while it participates in the v2 authority algebra; neither the table nor its
JSON/Markdown renderers may mislabel that identity as v2.

| Object kind | Authority algebra version | Identity schema version | Identity domain | Identity sources | Binding | No claim |
| --- | --- | --- | --- | --- | --- | --- |
| `claim-statement` | `2` | `2` | `frankensim.fs-govern.claim-statement.v2` | canonical conjunction clauses | clause-order invariant; clause mutation moves identity | does not prove statement truth |
| `quantified-domain` | `2` | `2` | `frankensim.fs-govern.quantified-domain.v2` | ordered quantifier blocks, explicit block commutativity, domain predicates | block order is semantic; only declared-commutative intra-block variables sort | does not prove satisfiability or nonvacuity |
| `assumption-set` | `2` | `2` | `frankensim.fs-govern.assumption-set.v2` | canonical assumption conjunction | assumption-order invariant; semantic mutation moves identity | does not discharge assumptions |
| `semantic-claim` | `2` | `2` | `frankensim.fs-govern.semantic-claim.v2` | statement, domain, assumptions, exact units, no-claim | semantic root excludes execution budget/seed/version/capability context | semantic identity is not an exact execution instance |
| `claim-lane-binding` | `2` | `2` | `frankensim.fs-govern.claim-lane-binding.v2` | statement/domain/assumption roots, validated lane, binder, artifact | claim instances reject a binding minted for another structured claim | binding data does not authenticate a dishonest binder |
| `claim-instance` | `2` | `2` | `frankensim.fs-govern.claim-instance.v2` | semantic claim, claim-lane binding, Five Explicits | semantic and exact-instance roots are distinct | content identity is not admission |
| `proof-lane` | `2` | `1` | `frankensim.fs-govern.proof-lane.v1` | validated LaneCharter | reuses non-forgeable lanes::ProofLaneId | lane identity is not proof |
| `evidence-ref` | `2` | `2` | `frankensim.fs-govern.evidence-ref.v2` | kind, exact claim, artifact, checker, schema | conclusion polarity is encoded before typed wrapper construction | reference is not authenticated authority |
| `satisfiable-evidence` | `2` | `2` | `frankensim.fs-govern.satisfiable-evidence.v2` | satisfiable evidence reference | type and identity both bind positive satisfiability polarity | candidate evidence is not authenticated satisfiability authority |
| `unsatisfiable-evidence` | `2` | `2` | `frankensim.fs-govern.unsatisfiable-evidence.v2` | unsatisfiable evidence reference | type and identity both bind negative satisfiability polarity | candidate evidence is not authenticated unsatisfiability authority |
| `nonvacuous-evidence` | `2` | `2` | `frankensim.fs-govern.nonvacuous-evidence.v2` | nonvacuous evidence reference, strength kind, context, fibre | positive polarity and exact strength are type- and identity-bound | one strength class cannot widen without an inference rule |
| `vacuous-evidence` | `2` | `2` | `frankensim.fs-govern.vacuous-evidence.v2` | vacuous evidence reference, strength kind, context, fibre | negative polarity and exact strength are type- and identity-bound | one strength class cannot widen without an inference rule |
| `reproduction-failed-evidence` | `2` | `2` | `frankensim.fs-govern.reproduction-failed-evidence.v2` | reproduction-failed evidence reference | type and identity both bind failed-reproduction polarity | failure evidence is not authenticated adjudication |
| `reproduced-evidence` | `2` | `2` | `frankensim.fs-govern.reproduced-evidence.v2` | reproduced evidence reference | type and identity both bind successful-reproduction polarity | reproduction is distinct from proof and admission |
| `evidence-state` | `2` | `2` | `frankensim.fs-govern.evidence-state.v2` | exact evidence reference, predecessor, lifecycle/cancellation fields | exclusive transitions replace the token; terminal states cannot revive | lifecycle completion does not establish statement truth |
| `authority-state` | `2` | `2` | `frankensim.fs-govern.authority-state.v2` | exact claim and all orthogonal authority axes | split scientific/exact refinement, clear runtime substitution, deny-biased meet | descriptive classifications are not authenticated authority |
| `exact-instance-decision` | `2` | `2` | `frankensim.fs-govern.exact-instance-decision.v2` | claim, policy, checker, verdict, artifact, schema | runtime validation requires admitted polarity and exact field agreement | decision candidate is not an authenticated admission receipt |
| `invalidation-binding` | `2` | `2` | `frankensim.fs-govern.invalidation-binding.v2` | target claim/state/head/generation and verified tombstone | private constructor binds the exact live predecessor atomically | binding is not cryptographic receipt authentication |
| `inference-rule` | `2` | `2` | `frankensim.fs-govern.inference-rule.v2` | name, version, definition artifact | default rule set is empty | registered rule is not assumed sound |
| `support-edge` | `2` | `2` | `frankensim.fs-govern.support-edge.v2` | source state, target claim/lane, rule, evidence | exact endpoint and rule identities | candidate edge is neither authenticated nor proof of graph acyclicity |
| `attack-edge` | `2` | `2` | `frankensim.fs-govern.attack-edge.v2` | candidate, target claim/lane, evidence | candidate target/domain must match | attack is not adjudication |
| `counterexample-candidate` | `2` | `2` | `frankensim.fs-govern.counterexample.v2` | target claim/domain and counterexample evidence | exact candidate-to-domain identity | candidate is not a refutation |
| `counterexample-adjudication` | `2` | `2` | `frankensim.fs-govern.counterexample-adjudication.v2` | candidate, target, verdict, adjudication evidence | only genuine verdict can derive a tombstone candidate | candidate adjudication cannot advance an authoritative head |
| `revocation-tombstone` | `2` | `2` | `frankensim.fs-govern.revocation-tombstone.v2` | target state, genuine adjudication, reason, evidence | permanent exact-state invalidation | descriptive tombstone is not an authenticated revocation receipt |
| `verified-revocation-tombstone` | `2` | `2` | `frankensim.fs-govern.verified-revocation-tombstone.v2` | authenticated tombstone candidate, target head and generation | opaque wrapper has no public minting constructor | wrapper alone makes no cryptographic-authority claim |
| `capability-policy` | `2` | `2` | `frankensim.fs-govern.capability-policy.v2` | axis/strength requirements, capabilities, accepted assumptions/no-claims | every guard changes policy identity | policy data is neither capability possession nor checker authority |
| `checker-decision` | `2` | `2` | `frankensim.fs-govern.checker-decision.v2` | claim, authority, policy, checker, verdict, artifact, cancellation | public candidate is exact data; opaque decision has no public Phase 0B-A mint | candidate verdict is neither authentication nor statement truth |
| `authority-head` | `2` | `2` | `frankensim.fs-govern.authority-head.v2` | claim, exact state, invalidation, generation, predecessor head | atomic advancement preserves permanent invalidation and replaces the head token | durable single-head authentication is Phase 0B-B scope |
| `runtime-admission` | `2` | `2` | `frankensim.fs-govern.runtime-admission.v2` | claim, authority, policy, checker decision, current head identity/generation | positive typestate plus exact product-policy and live-head validation | does not widen claim scope or survive revocation |
| `authority-migration` | `2` | `2` | `frankensim.fs-govern.authority-migration.v2` | legacy schema/record/rank/booleans and explicit demotions | v0 ambiguity demotes to v2 Unknown axes; persisted v1 is refused | migration never restores legacy positive authority |

All canonical encodings are typed, tagged and length-prefixed before
domain-separated BLAKE3. Text sets canonicalize whitespace, order and duplicate
spellings. `QuantifiedDomain` stores an ordered sequence of `QuantifierBlock`s.
Block order is semantic (`forall/exists` is not `exists/forall`), adjacent
same-kind blocks are not merged, and only a block explicitly constructed as
commutative sorts its variables. Ordered blocks retain variable order; names
must be globally unique across all blocks. Predicates remain a canonical
commutative conjunction. Unit equivalence is exact and
structural: rational scales reduce by GCD and like factor exponents combine;
there is no heuristic alias conversion. Reordering equivalent clauses,
assumptions, declared-commutative block variables, predicates, factors, versions or capabilities
preserves identity. Changing a quantifier, assumption, unit exponent, seed,
budget, version, capability, lane or no-claim boundary moves the appropriate
semantic/exact root.

### Orthogonal authority product

`AuthorityState` keeps these descriptive candidate axes distinct and private:

- truth: `Unknown | ConditionalProof | Proved | Refuted`;
- statement satisfiability: `Unknown | Satisfiable | Unsatisfiable`, carrying
  polarity-distinct `SatisfiableEvidence` or `UnsatisfiableEvidence`;
- nonvacuity: `Unknown | Nonvacuous | Vacuous`, carrying polarity-distinct
  `NonvacuousEvidence` or `VacuousEvidence`, each bound to an exact
  point/open-family/positive-measure/scale-family/custom strength, context and
  fibre, with no cross-polarity or satisfiability conversion;
- exact-instance admission: `NotEvaluated | Refused | Admitted`, bound to the
  typed `ExactInstanceDecisionCandidate` identity over exact claim, policy,
  checker, verdict, artifact and schema;
- proof-kernel checking: `NotChecked | KernelChecked`;
- empirical scale qualification: `NotQualified | ScaleQualified`;
- reproduction: `NotAttempted | Failed | Reproduced`, carrying distinct
  `ReproductionFailedEvidence` or `ReproducedEvidence`;
- invalidation: a private-representation `Clear` or exact
  `InvalidationBinding` state that public fields cannot construct.

The truth partial order is `Unknown <= ConditionalProof <= Proved`, with
`Refuted` on a separate branch sharing the `Unknown` bottom. Other positive and
negative branches likewise share their explicit unknown/not-evaluated bottom
and require exact evidence identity when satisfying a stronger requirement.
`scientific_evidence_refines` compares only clear descriptive scientific axes;
`exact_decision_refines` compares the exact-decision axis;
`is_safe_runtime_substitute_for` requires both relations and two clear states.
`deny_biased_meet` never drops invalidation, but is explicitly not a
substitutability greatest lower bound. Lattice-law claims apply only to the
clear descriptive evidence product. Different exact claims refuse composition;
distinct invalidation histories return a structured conflict. Invalidated
candidate data cannot pass public runtime-candidate assessment or form a
support edge.

Invalid combinations refuse during construction. In particular, an explicit
assumption-bearing proof is `ConditionalProof`, not unqualified `Proved`;
`Unsatisfiable + Nonvacuous`, `Refuted + Admitted`, `Unsatisfiable + Admitted`,
`Vacuous + Admitted`, and `Invalidated + Admitted` cannot exist as an
`AuthorityState`. Every evidence axis rechecks evidence kind and exact claim
identity.

`AuthorityGrant<S>` is a sealed typestate view, but positive and refuted grant
constructors are module-private. `CheckerDecision` likewise wraps a candidate
only after private receipt authentication. `AuthorityHead` is non-`Clone`,
non-`Copy`, invalidation-, generation- and predecessor-bound. Public
`AuthorityState::new` constructs only clear states. A descriptive tombstone
cannot invalidate a state: only an opaque `VerifiedRevocationTombstone` at the
private receipt boundary can atomically bind and advance the exact current
head. Wrong target, stale/replayed tombstones, invalidated-head transitions and
generation overflow fail before mutation. This algebra makes no cryptographic
authority claim. `RuntimeAdmission` is
non-`Clone`/non-`Copy`, binds the exact
head identity and generation, and must validate against the authoritative
current head on every private consumption. Phase 0B-A exposes these object
shapes and internal tests, not a public issuer or persistent single-head store.
An authenticated checker decision is state/policy-bound rather than
head-generation-bound, so an idempotent same-state head refresh may reuse that
decision while minting a distinct generation-bound admission.

`assess_runtime_candidate()` is the public, explicitly non-authoritative policy
surface. It validates exact claim/state/policy binding, an accepting candidate
verdict, an admitted exact-instance decision with exact policy/checker/artifact
agreement, current descriptive admission/invalidation, truth, accepted
conditional assumptions, satisfiability, exact nonvacuity strength,
kernel/scale/reproduction guards, declared Five-Explicits capabilities and
accepted no-claim boundaries. Its `RuntimeAssessment` result means only
`eligible candidate` and has no conversion to a grant, authenticated decision,
head, receipt or admission. Capability declarations are not possession;
Phase 0B-B's authenticated checker receipt must verify actual possession.

### Evidence lifecycle, adjudication, migration and extension rules

`EvidenceLifecycle` binds an exact evidence reference to `EvidenceState` under
`frankensim.fs-govern.evidence-state.v2`. The state advances only
`Proposed -> Checked -> Adjudicated`, or to a terminal `Failed`/`Cancelled`
state. Cancellation requires nonzero identities for request, drain and
finalize; terminal states cannot transition or revive. Lifecycle tokens are
non-`Clone`/non-`Copy`; validated transitions atomically replace the exclusive
token, and successor identity binds the predecessor. Recreating a descriptive
Proposed root remains possible, so durable single-writer/CAS enforcement is
explicitly Phase 0B-B.

A `CounterexampleCandidate` and `AttackEdge` bind the exact target and domain
but grant no refutation authority. Only a descriptive `GenuineCounterexample`
`CounterexampleAdjudication` can derive a `RevocationTombstone` candidate.
Out-of-domain, artifact-defect and indeterminate verdicts cannot derive one.
The candidate itself cannot create an invalidated successor. Private receipt
authentication produces the opaque verified wrapper consumed by the atomic
current-head transition; no public path turns raw hashes into revocation.

V0 historical records mixed rank/admission/reproduction booleans without exact
evidence bindings. `migrate_legacy_v0` accepts only schema 0, binds the retained
source record and every legacy field, and demotes all positive axes to
`Unknown`/`NotEvaluated`/`NotAttempted`. Persisted v1 is explicitly refused
because its flattened mixed quantifiers and polarity-ambiguous evidence cannot
be reinterpreted soundly; source evidence must be regenerated under v2.
Unknown or future schema versions refuse. Future migrations must be explicit
versioned functions and may never widen historical authority.

New inference machinery enters only through a positive-version
`InferenceRule` bound to its exact definition artifact. `DEFAULT_INFERENCE_RULES`
is empty: merely landing a theorem or rule does not add it to default
composition semantics.

`authority_log_json` is a non-authoritative diagnostic over state/policy/checker
candidate context. It validates that supplied context binds the exact state and
policy, names object/source/state/policy/checker identities, algebra/policy
versions, every axis, no-claim boundaries and a stable ranked remedy code. It
refuses rather than returning a record over `MAX_AUTHORITY_LOG_BYTES`; that cap
bounds the returned record, not transient `String` allocation during rendering.

### Evidence and no-claim boundary

G0/G3 coverage owns canonical identity equivalence/sensitivity, partial-order
and deny-biased clear-product laws, invalid product states, exact binding,
schema refusal, migration demotion, typed satisfiability/nonvacuity separation,
terminal cancellation, adjudication/revocation, policy guards and bounded log
shape. Compile-fail documentation prevents descriptive state/assessment values
from widening into opaque grants/admissions. The central batch must still
execute those suites and catalog drift checks; this module alone does not
implement Phase 0B-B wire decoding, durable admission receipts/current-head
storage, package/checker/ledger adapters, signatures, or production runtime
consumption. The separate `evidence_graph` module below supplies descriptive
Phase 0B-C graph/planning candidates without widening this authority boundary.

## Phase 0B-C support/threat planning (`evidence_graph` module)

`evidence_graph` is a pure, bounded planning layer over the exact descriptive
identities established by `evidence_contract`. It builds one immutable,
content-addressed graph snapshot before scoring any work. The snapshot keeps
claim instances, authority-state candidates, assumptions, evidence references,
checker/falsifier declarations, consumers, support edges, and attack edges
distinct. Each evidence node retains the exact evidence kind, artifact,
checker, schema version, claim, and content identity; standalone checker and
falsifier nodes are descriptive declarations, not authenticated executions.
Registration binds every authority state to its exact claim and every
counterexample to its exact target. Every counterexample and support/attack edge
must also name its exact in-snapshot evidence node bound to the target claim.
Unknown or mismatched endpoints, duplicate identities, self-support, and support
cycles refuse before a snapshot can exist. Attack edges may express mutual
constraints because an attack is not support and is not an adjudicated
refutation.

`EVIDENCE_GRAPH_VERSION = 2` is an identity-schema boundary. Every graph-node
preimage begins with tagged graph version field `0`, tagged authority-algebra
version field `250`, then the node-kind tag and its exact payload. The literal
node-kind tags are `Claim=1`, `AuthorityState=2`, `AssumptionSet=3`,
`Evidence=4`, `Checker=5`, `Falsifier=6`, `Consumer=7`, and
`Counterexample=8`; each preimage binds the exact payload fields documented by
its public variant. Snapshot preimages begin with the same two version fields,
then canonically sorted node identities under repeated field tag `1`, support
edge identities under tag `2`, and attack edge identities under tag `3`. The
production domains are
`frankensim.fs-govern.evidence-graph-node.v2` and
`frankensim.fs-govern.evidence-graph-snapshot.v2`. The separate test-fixture
namespace is `frankensim.fs-govern.test-evidence-graph.v2`; it creates fixture
inputs only and is never a production identity domain. Changing any node or
snapshot tag/payload contract requires a new graph identity schema.

Support reachability is evaluated only over the validated support DAG.
Consequence is the checked sum of explicitly declared downstream authority and
consumer weights reachable from a claim. Authority history cannot multiply a
claim: at most the maximum declared authority consequence is counted per exact
claim, while each consumer content identity is counted once. Consequence is not
inferred from graph degree or issue priority. An attack on an upstream claim
propagates once per attack edge to
every support-dependent claim. Doubt remains a product of named, bounded
components: calibrated uncertainty, uncovered attack surface, independence
shortfall, and unresolved assumptions. Because this layer has no authenticated
adjudication state, any graph-visible unadjudicated attack conservatively raises
the affected score's doubt to one. No favorable independence, correlation, or
attack-coverage value is inferred automatically: caller declarations are
content-bound and selected candidates sharing a declared independence or
correlation class are excluded, but their scientific adequacy remains
unauthenticated. Exact duplicate nodes/edges refuse. All arithmetic is checked
integer/fixed-point arithmetic, so overflow is a structured refusal rather than
a ranking change.

Allocation candidates bind the exact snapshot, claim, proof lane,
independence class, work kind, cost, priors/correlations, and a descriptive
anytime-accounting artifact with a nonzero observation count. Under optional
stopping, a policy may require that artifact; this crate does not manufacture or
validate an e-process. The deterministic planner ranks consequence-times-doubt
while enforcing an explicit no-action reserve, holdout/diversity/independence/
exploration floors, one selected bet per proof lane and declared independence
and correlation class, a selection cap, and a bounded work-unit budget.
Canonical identity ordering is the final tie-break. A bounded exact feasibility
search prevents greedy lane/correlation choices from falsely declaring a
feasible floor set impossible; exhausting its explicit search-state cap refuses
without claiming infeasibility. A successful decision candidate exposes the
graph/policy roots, utility model, sensitivity artifact, complete ranked rows
with the four doubt inputs, priors, correlations, and full descriptive anytime
state, selected reservations, used budget, no-action reserve, unused allocatable
budget, total unallocated budget, and whether no positive/floor-required action
was selected. The no-action model in this version is that explicit reserve plus
an empty-selection result; it is not a separately utility-ranked candidate.

The planner is transactional and synchronous: it validates and computes into
local candidate data, polls the supplied cancellation boundary at bounded
passes, and publishes nothing on refusal or cancellation. It does not mutate
`lanes::PortfolioLedger`; its scalar cost is only the work-unit projection of a
future reservation. A downstream authenticated adapter must separately admit
the selected work under the ledger's four-axis `ResourceEnvelope`, mechanism
cap, and comparison rules, then persist the graph/decision receipts.

Authority boundary: every public graph, anytime-accounting, score, allocation,
and decision object in Phase 0B-C is descriptive candidate data. BLAKE3 roots
provide deterministic content identity, not signatures, authorized review,
statistical validity, scientific truth, calibrated probabilities, durable
budget reservation, or runtime admission. Phase 0B-B and later checker/ledger
adapters must authenticate the exact referenced artifacts before any favorable
authority or capability can be consumed. Cohomological/tropical summaries and
automatic e-process execution remain future versioned inference-rule work;
they are not implied by this solid graph core.

## Proof-lane admission (`lanes` module, bead rjoq.6)

- `LaneCharter::new(statement, admissible_domain, assumptions,
  target_authority, baseline, falsifier_family, independence_class)`
  canonicalizes (whitespace collapse; assumptions sorted + deduped; every
  field non-empty and bounded before canonicalization allocates) and `lane_id()` mints the validated
  `ProofLaneId` (tagged, length-prefixed fields under BLAKE3 domain
  `frankensim.fs-govern.proof-lane.v1`). There is NO public id constructor
  from a raw hash — an id always corresponds to a validated charter
  (anti-spoofing), and cosmetic re-spellings collapse to one lane (the
  split-gaming canonicalization). `mechanism_id(name, version)` mints a
  `MechanismId` that retains both its originating `ProofLaneId` and its
  content identity; admission, comparisons, and supersession refuse a
  mechanism presented under another lane.
- `PortfolioLedger::new(PortfolioPolicy { global, max_active_mechanisms })`
  is the admission state machine (`LANE_POLICY_VERSION = 2`). Semantics:
  multiple active unproven mechanisms across independently falsifiable
  lanes; a second active mechanism in the SAME lane refuses atomically;
  lanes DECLARING the same independence class share one bet (the ledger
  retains the complete active set, so finalizing one comparison candidate
  cannot hide a surviving candidate); a comparison on another lane cannot
  evade that backstop; the four-axis global
  envelope (work / memory / reviewer slots / falsification capacity) and
  the mechanism cap bind across all lanes.
- `HeadToHeadCharter::new(&lane_charter, candidates, shared, preregistration)` is
  the ONLY carve-out: a preregistered comparison (declared BEFORE any
  admission in the lane, 2..=8 distinct candidates, non-zero
  preregistration artifact) admits exactly its declared candidates under
  its own bounded shared envelope.
- `FinalizationReceipt::new(mechanism, kind, superseded_by, ledger_artifact)`
  (kinds: refuted / tombstoned / withdrawn / superseded-with-successor) is
  the ONLY release path: `finalize` releases the slot and its reservation
  EXACTLY ONCE against an identity-consistent receipt; Unknown or stalled
  work never releases (there is deliberately no timeout path); terminal is
  permanent (no re-admission).
- Every request carries an `IdempotencyKey`: a byte-identical retry
  replays the recorded decision without double-charging or re-recording; a
  different request under a used key records one conflict row naming the
  original sequence; an exact retry of that conflicting request replays the
  same refusal without another row. Every method validates completely BEFORE
  governed state mutates; refusals may append their explicit audit decision but
  never partially charge lane/resource state. `PortfolioLedger` is deliberately
  non-`Clone`; exclusive `&mut self` access to that authority value is the
  in-process concurrency contract.
- The retained record is fail-closed and hard-bounded by
  `MAX_RETAINED_DECISIONS` and `MAX_RETAINED_DECISION_BYTES`; the primary and
  conflict-idempotency maps share the decision cap. Records are never silently
  evicted because eviction would permit an old retry to execute again. Instead,
  `RetentionCapacityExceeded` refuses before cloning variable-size request data.
  One row/key and a conservative byte allowance remain reserved for every
  active mechanism's future finalization.
- Every `AdmissionDecision` retains the complete policy and canonical replay
  preimage: statement/domain/assumptions/authority/baseline/falsifier/class,
  comparison candidates/artifact/shared envelope or admission reservation,
  terminal kind/successor/artifact/receipt identity/released envelope, plus
  lane/mechanism ids, idempotency key, request digest, verdict, and ranked
  remedy. `decisions_json(limit)` adds explicit skipped and retained-cap
  metadata; the stored request, rather than an opaque digest alone, is
  sufficient for deterministic replay.

## Moonshot portfolio (`moonshots` module, bead f85xj.16.3)

- `MoonshotDeclaration` is the typed policy record for one `[M]` lane. Its
  constructor requires a Bead id, independently falsifiable lane, named owner,
  `NamedFalsifier` observation and decision rule, positive effort cap, real
  calendar deadline, next review no later than that deadline, required review
  evidence, a boring fallback derived from the falsifier's decision rule, and
  a non-empty critical-path-disjointness boundary.
  Text inputs are bounded by `MAX_MOONSHOT_FIELD_BYTES` before the underlying
  `LaneCharter` and `MechanismId` are minted.
- `MoonshotPortfolio::new` projects those declarations into the existing
  four-axis `PortfolioLedger`. The v1 status-quo ceiling is the immutable six
  active lanes captured by `MOONSHOT_V1_INITIAL_CAP`; the supplied cap must
  equal the complete declaration count and may shrink but never exceed six.
  The sum of admitted effort caps becomes the portfolio work envelope, while
  each lane charges one reviewer and one falsification-capacity slot.
- New work has no independent admission method. `ReplacementAdmission` names
  exactly one active `DisplacementRecord` with the exact `completed`,
  `falsified`, or `shelved-with-state` disposition, retained state/evidence
  locator, and reason. `assess_replacement` consumes the portfolio, finalizes
  the predecessor through `PortfolioLedger`, admits the successor, verifies
  the active count stayed fixed, retains the exact domain disposition, and
  permanently refuses revival of a terminal Bead. Consuming `self` prevents a
  failed successor admission from exposing a partially released authority
  value to the caller. `liquidate` uses the same consumed-value transaction
  without a successor and decrements the wrapper cap permanently; later work
  can still enter only by replacing another currently active lane, so the
  released slot cannot be refilled.
- `MOONSHOT_DISPOSITION_IDENTITY_DOMAIN` binds the exact domain disposition,
  retained-state locator, reason, and displaced Bead into the finalization
  evidence root. This preserves the completed-versus-shelved distinction even
  though both conservatively map to the generic ledger's `Withdrawn` terminal
  kind; falsification maps to `Refuted`.
- The module is pure synchronous policy with no ambient clock, filesystem,
  tracker, or dependency-graph access. `xtask check-moonshots` is the explicit
  world adapter: it binds the immutable six-Bead baseline and disposition
  history to live Beads labels/status/assignment and checks direct/transitive
  reachability into `vertical-capability-graph.json`.
- No-claim boundary: constructing or admitting a declaration proves only that
  the governance fields are complete and resource-bounded. It does not prove
  the disjointness prose, authenticate a Bead owner or artifact, evaluate the
  falsifier, validate the science, or promote `[M]` work. Those remain live
  graph, evidence, checker, and review responsibilities.
- `tests/moonshots.rs` provides G0 construction/date/cap refusals, one-for-one
  disposition retention, terminal anti-revival, effort-envelope refusal, and
  cap-shrinking liquidation followed by replacement at the smaller cap.

## Crate registry (`crates` module)

- `addendum_crates() -> &[AddendumCrate]` — the seven net-new crates the
  addendum introduced, each `{ name, purpose, owning_proposal, layer, no_claim
  }`. `crate_audit() -> CrateAudit` confirms each declares a purpose, an owner,
  and a no-claim boundary (the AGENTS.md contract discipline made
  governance-legible); `crates_json()` emits the deterministic record. Actual
  `CONTRACT.md` file presence is enforced separately by `xtask check-contracts`.

## Requirement-to-evidence traceability (`traceability` module)

- `requirements()` is the closed canonical B1–B14 and sixteen-RQ source
  registry derived from the extension charter seed. Every row carries the
  exact generated-ledger columns: stable requirement id, capability/property,
  concrete blocker, owner/artifact, prerequisite phase, milestone, forcing
  flagship, benchmark/data evidence route, proof-obligation links, honest
  claim boundary, and declarative status. Status is governance declaration,
  not scientific proof and not a mirror of whether an owning Bead is closed;
  tracker `closed` is deliberately not an accepted scientific lifecycle value.
- `proof_obligations()` is the complete ordered PO-1 through PO-25 index. Each
  entry links its executable summary to the full owner Bead ids. A row may
  cite only entries in this closed index.
- `audit_traceability(rows, obligations)` accumulates structured
  `TraceabilityDiagnostic { requirement_id, field, reason }` refusals. An empty
  scope, any missing or oversized field, missing owner, missing milestone gate,
  missing benchmark/data evidence route, missing claim boundary, missing or
  dangling PO link, duplicate requirement, noncanonical PO alias, incomplete
  PO index, or duplicate PO owner is non-green. Diagnostics are sorted by
  requirement, field, and reason, so source enumeration order cannot alter the
  refusal artifact.
- `generate_traceability_ledger(rows, obligations)` validates first and emits
  no partial output on failure. Successful rows use the closed charter order,
  PO links numerically, PO definitions numerically, and owner ids lexically;
  `traceability_ledger_json()` applies it to the canonical code registry under
  schema `frankensim-requirement-traceability-v1`. Every artifact explicitly
  emits `authority: "declaration-only"` and `source_snapshot: null`; this API
  refuses `verified` and `validated` because only a future source-bound
  admission receipt may mint scientific promotion status.
- `TraceabilitySourceSnapshot::new(sources)` is the pure adapter boundary for
  live Beads/contracts/registry extraction. It requires at least one immutable
  artifact from each of those three classes, rejects blank/oversized or
  duplicate locators and the all-zero missing-content sentinel, canonicalizes
  source order, and derives a domain-separated BLAKE3 snapshot identity from
  source kind, locator, and exact content identity. It reads no filesystem and
  trusts no tracker status implicitly; the caller owns bounded parsing of the
  named bytes.
- `load_traceability_source_snapshot(root, specs, limits)` is the concrete,
  synchronous filesystem boundary. It accepts only strict root-relative
  locators, canonicalizes every regular file beneath the configured root,
  refuses final symlinks, bounds both individual and aggregate reads, rejects
  metadata/read length races, and hashes the exact admitted bytes. Its receipt
  records the adapter version, effective limits, canonical locator/class,
  byte count, exact BLAKE3 identity, Beads record count and canonical id index
  where applicable, and the resulting source-snapshot identity. Contract and
  registry sources must be nonempty NUL-free UTF-8. Beads sources additionally
  receive a bounded lexical JSON-object-envelope audit and a unique canonical
  top-level string `id` audit before snapshot admission; the same id appearing
  in two bound Beads sources refuses as an ambiguous index.
- `audit_traceability_owner_join(obligations, loaded)` checks that every Bead id
  named by the complete PO index is lexically present in those exact loaded
  Beads bytes. Its sorted diagnostics name the missing PO and owner id.
  `generate_traceability_ledger_from_loaded_sources(rows, obligations, loaded)`
  validates declarations first, performs that presence join, and returns no
  partial ledger on either failure. A green `TraceabilityOwnerJoinReceipt`
  binds the join version, owner-reference/distinct-owner/source/index counts,
  and exact source-snapshot identity. It does not inspect or grant meaning to
  tracker status, assignment, dependencies, closure, or evidence fields.
- `generate_traceability_ledger_from_snapshot(rows, obligations, snapshot)`
  retains `authority: "declaration-only"` while adding three explicit roots:
  the admitted source-snapshot identity, the canonical unbound declaration
  identity, and a binding identity over that pair. A source-byte mutation moves
  the snapshot/binding roots; a row/PO mutation moves the
  declaration/binding roots. None of these roots proves adapter correctness,
  scientific adequacy, authorized review, or that Bead closure discharges a
  proof obligation.
- Inputs are hard-bounded before rendering: at most 256 requirements, exactly
  the closed 25-PO definition space, 25 PO links per requirement, 16 owners per
  PO, 16 KiB per scalar field, 512 source artifacts, and 4 KiB per source
  locator. The core `traceability` module remains pure. `traceability_fs` owns
  bounded standard-filesystem reads only; it does not open the Beads database,
  perform network I/O, or persist immutable artifacts. `traceability_join`
  owns only canonical owner-id presence. Higher tooling owns lifecycle,
  dependency, contract-clause, evidence, and persistence semantics and must use
  these validators instead of hand-maintaining a dashboard.

## Vertical ratification (`ratification` module, bead f85xj.1.4)

- `VERTICAL_RATIFICATION_V1` is the program's one ratified vertical decision
  record: chosen vertical (`thermal-design-assurance`), runner-up and retained
  minority report, the recorded scoring totals and inventory pin they rest on,
  the kill criterion bound to the measured cycle-time baseline, three named
  mechanically-evaluable `Falsifier`s, the next review date, honored
  principles, and the downstream beads gated on the record id.
- The record is fail-closed. `ratified_vertical()` re-validates before
  handing the record out: required fields non-empty, falsifiers present and
  complete, kill parameters equal to `fs_wedge::CHT_BASELINE`'s
  (`RatificationError::KillCriterionUnbound` otherwise), the baseline's
  provenance non-placeholder and its record complete
  (`BaselineNotMeasured`), and the scoring table recomputed via
  `fs_wedge::default_recommendation()` with any disagreement refused as
  `ScoringDrift` (recommendation, runner-up, totals, and the inventory
  revision on every comparison candidate). `ratification_json()` renders
  every record in `decision_records()` fail-closed and deterministically.
- The conclusion is mirrored in `fs-wedge` (`RATIFIED_VERTICAL`,
  `RATIFICATION_RECORD_ID`) and cross-checked by this crate's battery so the
  two crates cannot fork the story. `fs-govern` therefore carries a
  UTIL-to-UTIL path dependency on the dependency-free `fs-wedge`.
- No-claim boundary: ratification is a COMMERCIAL decision record over the
  measured comparison. It confers no technical maturity, does not upgrade any
  capability level, and does not validate the thermal market thesis — the
  customer-pain factor remains a declared assumption until interviews or
  workflow traces exist. A failed quarterly falsifier forces re-selection of
  the wedge, not a change to the platform architecture.

## Measured wedge decision audit (`wedge_audit` module, bead f85xj.1.5)

- `WedgeDecisionAuditRequest::new` requires a positive finite measured
  FrankenSim cycle time and a non-empty bounded evidence locator. The locator
  is retained exactly; no default measurement or anonymous scalar can enter
  the audit.
- `build_wedge_decision_audit` fails closed across the complete current
  decision surface: `fs_wedge::audit`, the measured comparison and exhaustive
  one-factor sensitivity render, the sourced cycle-time baseline and
  conservative kill verdict, and `ratified_vertical` / `ratification_json`.
  Success emits one JSON artifact under schema
  `frankensim-wedge-decision-audit-v1`, with a domain-separated BLAKE3 identity
  over the exact payload and a 1 MiB publication cap.
- Ordered `WedgeAuditLog` JSON-lines rows retain the validated request, every
  fs-wedge self-check, the complete measured-input/baseline manifest, scoring
  report and sensitivities, ratification record, kill derivation, and final
  artifact. Warning rows carry a stable refusal code and an actionable fix.
  The `wedge-audit` binary writes the artifact to stdout and logs to stderr;
  argument refusals exit 2 and decision-boundary refusals exit 3 without
  publishing a partial artifact. Its explicit
  `--seed-fault missing-cycle-time-evidence` drill proves that removing the
  measurement provenance fails the lane.
- The artifact is a reproducible commercial-decision audit, not an authority
  upgrade. The caller-supplied evidence locator is not opened, authenticated,
  content-hashed, or proven to name a production timing run. A successful
  render does not prove market demand, technical maturity, scientific
  validity, or that a production run met the kill criterion.

## Doctrine and proposals (`doctrine`, `proposals` modules)

- `principles() -> &[Principle]` — the eight design principles P1–P8 (id, name,
  statement); `rules() -> &[GovernanceRule]` — the four governance rules
  (number, name, statement).
- `proposals() -> &[Proposal]` — the nineteen proposals in composite (Mean)
  order, each `{ id, name, phase, mean, kill_metric, owning_bead, receipt }`.
  `governance_audit() -> GovernanceAudit` enforces that every proposal
  DECLARES a kill metric AND an owning bead (Governance Rule 2), counting how
  many are instrumented; `proposals_json()` emits the deterministic
  machine-readable record.

## Expansion-program register (`program_risks` module)

- `ProgramRiskId` is the namespace `PR-001` through `PR-012`; it is deliberately
  disjoint from the original `RiskId::R1` through `R10` runtime register.
- Every `ProgramRisk` row carries a named workstream and owning Bead, likelihood
  and impact, a leading indicator, a numeric comparator/threshold/unit/typed
  domain/minimum-sample trigger, mitigation, contingency (kill/refuse/escalate),
  residual likelihood and impact, and one E0–E7 review gate.
- `program_risks()` returns the twelve rows in canonical id order;
  `program_risk(id)` performs total typed lookup; and
  `program_risk_register_json()` emits the deterministic ledger artifact.
- `assess_program_risks(observations)` always emits twelve ordered rows. Missing,
  duplicate, non-finite, unit-mismatched, out-of-domain, and under-sampled
  evidence is non-green. Only one finite, correctly unit-tagged, in-domain,
  sufficiently sampled observation below its declared trip condition is
  `Clear`. Caller units are retained only as a 64-byte UTF-8-safe preview plus
  their exact original byte length.

## Public types and semantics

- `GraphNode` has typed constructors for exact claims, authority-state
  candidates, assumption sets, evidence, checkers, falsifiers, consumers, and
  counterexamples. `GraphSnapshot::new` canonicalizes those nodes plus the
  Phase 0B-A support/attack edges and publishes one root only after endpoint,
  exact counterexample/edge-evidence, duplicate, and support-DAG validation.
- `FixedRate` and `DoubtProfile` retain the exact calibrated-uncertainty,
  attack-coverage, independent-support, and assumption-resolution inputs.
  `DoubtProfile::combined` conservatively rounds the union of named doubt
  sources outward on a one-million-part fixed-point lattice.
- `AnytimeAccountingCandidate` binds a named method/version, observation
  count (nonzero), state root, and evidence artifact without claiming the method
  or artifact is valid. `AllocationCandidate::new` additionally requires the
  exact graph, claim, and validated `LaneCharter`; it refuses a charter that
  does not mint the claim's proof lane and derives the independence class from
  that charter rather than a caller-supplied bare class hash.
- `plan_allocations` and `plan_allocations_with_cancel` return an immutable
  `AllocationDecisionCandidate` containing the complete deterministic ranking,
  selected reservation candidates, exact used/reserved/unused/unallocated
  budget decomposition, explicit no-action result, inspectable
  doubt/prior/correlation/full-anytime state, policy/model/sensitivity roots,
  and content identity.
- `RiskId` (`R1`..`R10`) with `RiskId::ALL` and `code()`.
- `InstrumentationReceipt::new(subject, dashboard, verifier,
  evidence_artifact, verified_day)` validates mandatory provenance and returns
  an opaque receipt. Its private fields prevent accidental identity drift;
  accessors expose the dashboard, verifier, evidence-artifact content hash,
  verification day, and receipt identity. `receipt_identity()` is the replay
  oracle.
- `Risk { id, name, description, mitigation, early_warning, threshold, owner,
  receipt }` — `early_warning` is the metric that makes the risk visible before
  it is fatal; `owner` is the bead that owns the mitigation; `receipt` is the
  optional evidence-bearing instrumentation assertion (`None` is the honest
  baseline).
- `register() -> &'static [Risk]` — the canonical R1–R10 in order;
  `risk(id) -> &'static Risk` for lookup.
- `audit(today_day) -> RiskAudit` / `audit_slice(&[Risk], today_day) ->
  RiskAudit` — checks every
  risk has a non-empty early-warning metric AND an owner, counts how many are
  instrumented, and separately lists schema gaps and operational receipt gaps.
  `declared_schema_ok()` and `operationally_managed()` deliberately expose the
  two different verdicts. Empty audit scopes and duplicate risk ids fail closed;
  success requires the declared/verified counts to equal the nonzero total.
- `to_json(today_day) -> String` — a deterministic machine-readable JSON array
  with JSON-escaped strings for dashboards / CI gates. Every row carries the
  unambiguous instrumentation status and either `receipt:null` or the complete
  receipt provenance (dashboard, day, verifier, evidence artifact, identity).

## Trust boundary: declaration vs live operation (bead xpck.9)

The audits report TWO verdicts, never one: `declared_schema_ok()`
(every entry names a metric and an owner — pure schema) and
`operationally_managed()` (every metric VERIFIED live). The former
single `ok()` collapsed these and rendered the zero-instrumented
registry as green — the false-green this bead removed. Instrumentation
is EVIDENCE, not a flag: an entry counts as verified only through an
`InstrumentationReceipt` that binds the subject id, dashboard locator,
verifier identity, supporting evidence-artifact content hash, and verification
day. The canonical encoding uses tagged, `u64`-length-prefixed fields under
BLAKE3 derive-key domain
`frankensim.fs-govern.instrumentation-receipt.v1`; all receipt fields are
private. Subject replay, a future verification date, missing provenance, stale
evidence, or an inconsistent identity fails closed (`BadReceipt`/`Stale`).
Audits and JSON take `today_day` (days since 2026-01-01) explicitly, so verdicts
are deterministic and replayable.

The BLAKE3 root is an **unkeyed content identity, not an authentication tag or
signature**. It provides collision-resistant identity and accidental-tamper
detection; it does not prove that the dashboard was live, that `verifier` was
authorized, or that the evidence artifact is scientifically adequate. The
canonical registry remains code-reviewed governance data, and issuer trust /
artifact checking are deployment policy. Calling a public hash an
"authentication fingerprint" would overstate this crate's security contract.

## Invariants

- The register declares all ten risks while honestly remaining operationally
  red until receipts are installed.
- `register()` and `RiskId::ALL` share the same order.
- `to_json()` and `audit()` are deterministic.
- Receipt identities change when any semantic field changes, are bound to one
  governed subject, and cannot be mutated through the safe public API.
- The program register has exactly twelve unique rows in `ProgramRiskId::ALL`
  order, and every trigger contains a finite number, explicit unit, typed
  numeric domain, and positive sampling floor.
- Program assessment is input-order independent and cannot become all-clear by
  omission, duplication, malformed units, invalid numeric domains, or NaN.
- Graph snapshot identity is input-order invariant and moves when any node or
  support/attack edge identity moves. Edge evidence must be present and bound to
  the exact target. Consequence counts at most the maximum authority weight per
  reachable exact claim and each consumer identity once; support decomposition
  or authority history cannot multiply it.
- Attacks propagate through downstream support reachability exactly once per
  attack edge and cannot reduce the affected consequence-times-doubt score.
- Allocation rankings are input-order invariant. Selected work contains at
  most one candidate per exact proof lane, validated independence class, and
  declared correlation class; selected cost never enters the no-action reserve
  or exceeds the bounded work-unit budget. The full four-axis portfolio cap is
  a separate `PortfolioLedger` admission and is not asserted by this planner.
- Refusal or cancellation returns no partial allocation decision. Fixed-point
  doubt and consequence/score/budget arithmetic either remain representable or
  fail with a structured overflow.
- Every mechanism is admitted only under the lane that minted it; every
  surviving member of a comparison independently keeps its declared
  independence class occupied.
- Retained decisions, primary/conflict idempotency bindings, and their
  variable-size canonical payloads remain within fixed caps. Capacity reserved
  for active terminal transitions is unavailable to unrelated traffic.
- The traceability registry contains exactly the thirty charter requirement
  ids and the complete canonical PO-1..PO-25 index. Generation is deterministic
  and cannot emit a partial ledger when any owner, gate, evidence route, proof
  link, boundary, or status is absent.

## Error model

`InstrumentationReceipt::new` returns `ReceiptError` for an empty subject,
dashboard, verifier, or an all-zero missing-evidence artifact sentinel. Audits
report missing, stale, future-dated, subject-mismatched, and otherwise
inconsistent receipts as data and fail closed; they do not panic or silently
promote coverage.

The expansion-program evaluator is total over its bounded slice input: evidence
defects are represented by `AssessmentStatus`, not panics. The caller remains
responsible for bounding the number of supplied aggregate observations; output
and retained unit text remain bounded by the fixed twelve-row register and unit
preview limit.

Proof-lane APIs return `LaneError` for empty/oversized preflighted inputs,
lane/mechanism mismatches, occupied lanes/classes, envelope/cap failures,
invalid terminal receipts, idempotency conflicts, and exhausted retained-log
capacity. Capacity exhaustion never evicts replay authority and never mutates
governed admission state.

Evidence-graph constructors and planners return `GraphError` for missing or
oversized inputs, invalid rates/policies, duplicate semantic identities,
unknown or mismatched endpoints, self/circular support, arithmetic overflow,
snapshot/lane mismatch, duplicated work, absent required anytime-accounting
data, infeasible floors, and cancellation. Errors never fall back to an empty
graph, zero consequence, or partially published decision. A bounded feasibility
search that reaches its state cap returns an explicit search-limit error rather
than a false infeasibility verdict.

Traceability generation returns a `TraceabilityAudit` rather than output when
the source is empty, oversized, duplicated, incomplete, or contains a dangling
reference. Diagnostics name the exact requirement and field; no missing-field
case is collapsed into an unstructured string or silently defaulted.

## Determinism class

Fully deterministic — pure functions over `const` data, no RNG or I/O.

## Cancellation behavior

Most modules are synchronous pure functions with no cancellation surface.
`plan_allocations_with_cancel` is the bounded exception: it polls the supplied
pure cancellation predicate during input validation, every active floor pass,
each feasibility-search state, ranking, and both before and after receipt
construction. Cancellation returns
`GraphError::Cancelled` and no decision candidate. This is a planning-boundary
protocol, not an asupersync execution scope or durable request-drain-finalize
receipt.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/evidence_graph.rs` (Phase 0B-C G0/G2/G3/G4): canonical snapshot and
decision permutation invariance; exact identity mutation; duplicate/unknown/
mismatched endpoint and circular-support refusals; consequence reachability,
diamond deduplication and support decomposition; attack accounting; checked
overflow; fixed-rate doubt laws; score/tie determinism; budget conservation,
floors, no-action, one-lane/independence/correlation backstops and anytime
accounting; seeded consequence-times-doubt priority without starvation; prior,
utility and correlation sensitivity; and cancellation with no partial result.

`tests/register.rs` (Part V, 10 cases): all ten risks present + ordered;
every risk has a metric/owner/mitigation; owners are real bead ids; lookup;
the canonical audit is complete with an honest zero-instrumented baseline;
`audit_slice` detects missing metric AND owner on an incomplete entry, rejects
empty scopes and duplicate ids (the audit is not vacuous), and fails closed on
subject-replay/future/stale receipts;
missing provenance is rejected; changing subject, dashboard, verifier,
evidence artifact, or day changes the identity; JSON includes explicit receipt
provenance; determinism.

`tests/governance.rs` (8 cases): eight principles P1–P8; four rules numbered
1–4 (Rule 2 = kill-criteria enforcement); all nineteen proposals present with
unique ids and in descending composite order; every proposal declares a kill
metric + bead-id owner; the governance audit is complete with a zero
instrumented baseline; owner mapping spot-checks; `proposals_json` is
well-formed + deterministic.

`tests/program_risks.rs` (G0): exactly twelve ordered unique rows; an exact
twelve-row owner/trigger/review-gate mapping lock; complete quantitative
schemas and review gates; exact comparator boundaries; fail-closed
missing/duplicate/non-finite/unit/domain/sample handling; integer/fraction
domain checks; input-order-independent assessment JSON; and deterministic
numeric-threshold register serialization.

`tests/traceability.rs` (G0): exact thirty-requirement and PO-1..PO-25 closed
inventories; canonical completeness and link coverage; one named diagnostic
for every required row field; deliberately orphaned owner detail; duplicate
requirement and dangling/repeated PO refusals; missing/noncanonical/duplicate-
owner PO index failures; success and complete failure-audit permutation
invariance; empty-scope and exact/plus-one row, PO-link, owner, scalar, id,
summary, and reference caps; declaration-only authority enforcement; and exact
generated JSON column/index coverage.

`tests/traceability_fs.rs` (G0/G3): exact byte/hash receipts, canonical source
order, source-mutation identity sensitivity, declaration-only ledger binding,
bounded Beads record/line/nesting envelope and unique-id refusals, individual
and aggregate byte caps, strict relative paths, regular-file metadata, complete
source-class coverage, limit validation, duplicate-locator/cross-source-id
refusal, complete PO-owner presence, exact missing-PO/owner diagnostics,
declaration-before-join failure order, and proof that tracker status changes
move the bound source identity without changing lexical owner admission.

`tests/lanes_e2e.rs` (bead rjoq.6 slice 2): the cross-crate no-mock
composition — fs-govern admission persisted into a FrankenSQLite-backed
fs-ledger (events + content-addressed preregistration/refutation/
decision-log artifacts), fs-package claims re-checked solver-free by
fs-checker incl. a mismatched-root refusal, fresh-ledger byte-for-byte
replay, and an idempotent full-retry pass.

`tests/wedge_audit.rs` (bead f85xj.1.5, G0/G3/e2e): byte-identical artifact
and log replay, complete measured-input/scoring/sensitivity/baseline/kill/
ratification composition, typed missing-evidence refusal, deterministic binary
stdout/stderr, a seeded missing-evidence warning with an actionable repair, and
stable duplicate-argument refusal with no partial artifact.

`tests/lanes.rs` (bead rjoq.6, cases lane-001..lane-009): G0 identity and
state-machine laws (canonical collapse, per-field mutation sensitivity,
lane-bound mechanism/comparison/successor authority, same-lane refusal with
unchanged state, terminal exactly-once release, receipt validation incl. zero
evidence and self-supersession); G3
split/merge adversaries (independence-class backstop, comparison-evasion
refusal, undeclared-candidate refusal, surviving-candidate class retention);
G4 crash/retry idempotency (replay without double-charge, one-row stable
key-conflict refusal, refusal replay); G0 global mechanism-cap and independent
work/memory/reviewer/falsification-capacity boundaries for both global and
comparison envelopes; G5 whole-ledger and acceptance-complete JSON replay with
explicit truncation/cap metadata; bounded-retention refusal with a terminal
slot held for active work. Each same-lane, global-cap, terminal-release, and
identity guard has a test that fails if the guard is removed.

## No-claim boundaries

- Evidence-graph nodes, edges, anytime-accounting records, doubt inputs,
  correlation classes, scores, and allocation decisions are descriptive
  candidates. Their public constructors and BLAKE3 identities do not prove
  endpoint evidence is authentic, uncertainty is calibrated, correlations are
  honestly declared, optional-stopping validity holds, a selected reservation
  is durable, or a claim is true. The planner does not execute `fs-eproc`,
  mutate `PortfolioLedger`, or widen into the Phase 0B-A private authority
  types. Those are explicit downstream checker/ledger responsibilities.
- This crate encodes the risk register as governance DATA; it does not itself
  measure an early-warning metric, fetch an evidence artifact, authenticate a
  verifier, or prove dashboard liveness. A dashboard/CI supplies that evidence
  and deployment policy establishes issuer authority. The audit enforces that
  each risk declares a metric + owner and fails closed when receipt evidence is
  absent or malformed; it cannot establish the truth of an issuer's assertion.
- The original R1–R10 register and the PR-001–PR-012 expansion-program register
  are separate authorities. Neither namespace silently substitutes for the
  other, and neither measures its own leading indicators.
- Program-risk thresholds are governance trip points over caller-supplied
  aggregates. They do not establish scientific validity, automatically execute
  a contingency, or prove that the named Bead owner has reviewed the result.
- Bead-id owners in governance declarations remain string references. The
  filesystem adapter may read exact caller-selected Beads JSONL bytes, but it
  validates only a conservative lexical object envelope and unique canonical
  top-level ids. The owner join additionally confirms that every declared PO
  owner id is present. Neither path opens the Beads database nor interprets
  status, dependencies, closure, assignment, ownership truth, or scientific
  evidence semantics.
- The canonical traceability rows are governance declarations. They do not
  attest that an owner exists in the current tracker snapshot, that a benchmark
  ran, that a milestone closed, or that a scientific claim is verified. A
  filesystem receipt binds exact Beads/contracts/registry bytes and the owner
  join proves lexical owner-id presence only. Higher tooling must join admitted
  lifecycle, dependency, contract-clause, and V&V evidence semantics before
  persisting an immutable generated artifact; a closed Bead must never be
  translated directly into scientific proof status.
- Filesystem loading uses bounded synchronous `std` I/O and has no `Cx`; it
  makes no cancellation-latency, hostile-directory race resistance, sandbox,
  signature, authorization, or durable-persistence claim. Deployments that
  admit adversarially mutable roots need a stronger capability-scoped adapter.
- The lane ledger is the admission STATE MACHINE, not durable storage: the
  `FinalizationReceipt`'s `ledger_artifact` and the comparison's
  preregistration artifact are content references whose durable
  finalization, issuer authority, and scientific adequacy are established
  by fs-ledger/fs-package/fs-checker integration and deployment policy.
  The cross-crate no-mock E2E lives in `tests/lanes_e2e.rs` (dev-deps
  only): two independent theorem lanes plus one preregistered in-lane
  comparison drive fs-govern admission with every decision row persisted
  as a real fs-ledger event, the preregistration and refutation
  artifacts content-addressed in the ledger (the finalization receipt
  seals a hash that actually exists there), the outcome packaged as
  fs-package claims and re-checked solver-free by fs-checker (with a
  mismatched-root refusal probe), and the whole request sequence
  replayed byte-for-byte on a fresh portfolio ledger plus an idempotent
  full-retry pass. The G4 persistence-fault drill (same file) crashes
  the process mid-sequence — in-memory portfolio dropped, design-ledger
  handle closed — reopens the same path, proves pre-fault decision
  artifacts survive byte-for-byte, re-persists the recovered prefix as
  a DEDUPE (no double write), and converges the full re-execution to
  the never-crashed control bytes with idempotent retries absorbed.
  In-crate retention exhaustion (the governed store's own fault
  surface) is pinned by lane-009. NOT claimed: kernel-level I/O fault
  injection, torn/partial-write simulation, and media corruption —
  those live with fs-ledger's own crash-recovery battery and
  deployment policy.
- Independence classes are DECLARED. Canonicalization defeats cosmetic
  splits and the class backstop defeats partition gaming among honestly
  labeled lanes, but the crate cannot algorithmically prove that two
  falsifier families are genuinely independent — adversarial mislabeling
  is a governance-review matter, bounded in damage by the global caps.
- Non-`Clone` `&mut self` atomicity covers one in-process authority value and
  its retries; a
  multi-process admission service needs an external serialization or the
  ledger-backed successor.
