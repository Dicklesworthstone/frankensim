# fs-euler-disc-e2e contract

Status: [S] structural scientific-contract infrastructure plus a bounded,
deterministic numerical/software campaign slice. The slice composes real
squat-disc line/arc geometry and mass properties, a profile unilateral-sticking
micro-trajectory, a one-way flexible-base load response, encoded reduced-decay
ablations, and a one-way exterior passivity probe. It is not a closed coupled
Euler-disc model, an experiment or video validation, or support for claims
about Steve Mould's reported observations.

## Purpose and layer

fs-euler-disc-e2e is an L6 HELM leaf for the Euler-disc emergent-prediction
flagship. Its first schema freezes what questions may eventually be asked
before geometry, solver, calibration, or experiment choices can bias those
questions.

The committed `euler_disc_campaign` binary is a deterministic JSONL runner for
the bounded slice. It is deliberately composition-only: each rung consumes a
declared snapshot or one-way input from its owning production operator. It does
not solve mutual contact/base/exterior feedback or identify a physical decay
law.

The crate composes 15 artifact roles across four ownership surfaces. Their
rows live in a separately versioned, domain-separated owner-matrix routing
registry. The registry is a leaf identity schema: each source-schema string is
an opaque routing address whose bytes are frozen by the registry, not an
identity-dependency edge to the addressed implementation. The matrix is
deliberately exhaustive: adding a public protocol artifact without adding its
exact owner/schema/ceiling row is a schema error.

| Role | Exact owner | Source schema | Local authority ceiling |
| --- | --- | --- | --- |
| Context of Use, QoIs, criteria, applicability | fs-evidence | org.frankensim.fs-evidence.vv-artifact.v3 | structural context declaration |
| generic V&V evidence artifact | fs-evidence | org.frankensim.fs-evidence.vv-artifact.v3 | structural evidence reference only |
| generic V&V schema-admission receipt | fs-evidence | org.frankensim.fs-evidence.vv-schema-admission-receipt.v2 | structural schema admission only |
| claim IDs and directed evidence-use vocabulary | fs-ir | fs-ir:experiment-campaign-schema-v1 | campaign vocabulary only |
| canonical no-claim set | fs-govern | fs-govern:authority-algebra-v2 | canonical prohibition only |
| hypothesis-source declaration | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.hypothesis-source-declaration.v1 | hypothesis only |
| Euler claim DAG | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.claim-graph.v1 | structural claim policy only |
| complete scientific contract | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.scientific-contract.v1 | candidate eligibility only |
| evidence-reference packet | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.claim-evidence-packet.v1 | candidate eligibility only |
| direct prerequisite receipt | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.prerequisite-assessment-receipt.v1 | structural dependency only |
| claim-policy assessment | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.claim-policy-assessment.v1 | candidate eligibility only |
| contract-check receipt | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.contract-check-receipt.v1 | structural check only |
| claim-policy assessment log | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.claim-policy-assessment-log.v1 | diagnostic retention only |
| owner-matrix routing registry | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.owner-matrix.v1 | structural routing registry only |
| Context-bound aggregate-QoI derivation receipt | fs-euler-disc-e2e | org.frankensim.fs-euler-disc-e2e.aggregate-qoi-derivation-receipt.v1 | structural evidence reference only |

Lower layers must never depend on this crate. This crate does not redefine
generic V&V artifacts, instantiate a fictional ValidationPlan, construct a
full experiment campaign, or mint a governance grant. Those operations need
real, content-bound experiment, split, diagnostic, custody, and review
artifacts from later work.

## Bounded campaign runner

Run the committed executable through its retaining wrapper with new output
paths:

```bash
scripts/e2e/euler_disc_campaign.sh \
  --output target/euler-disc-campaign/campaign.jsonl \
  --stderr-log target/euler-disc-campaign/campaign.stderr.log
```

The wrapper requires committed, clean `Cargo.toml`, `Cargo.lock`, and `crates/`
inputs, runs the binary through strict remote RCH, refuses output or stderr
paths that already exist, and rejects partial or malformed JSONL. Output and
log paths are exclusively created and identity-bound for the run. A narrowly
handled `RCH-E309` may recover exact LF-framed producer records only when one
ordered, worker-consistent receipt chain proves that the submitted remote
campaign command exited zero before retrieval of its separately built artifact
failed. The retained receipt still marks artifact retrieval incomplete and its
local transcript authority as non-cryptographic. Other nonzero exits refuse;
ordinary RCH success with empty direct stdout also refuses rather than inferring
records from diagnostics. The wrapper emits, in order: sharp and
1-mm-filleted squat-disc line/arc geometry-plus-mass records;
a conservative steady oracle; a profile unilateral-sticking micro-trajectory;
a one-way reduced flexible-base record; contour-only, boundary-layer-only, and
combined reduced-decay records; a one-way reduced exterior-wrench passivity
record; and a campaign manifest. The manifest is a deterministic digest of the
preceding records, not a validation certificate.

This contract deliberately records no output path, digest, or numerical value
until an actual retained run supplies them. The runner's receipts demonstrate
only the encoded numerical/software composition and its input/record checks.

## Public types and semantics

- OwnerMatrix is an independently versioned canonical v1 transport with its own
  domain-separated identity and strict encode/decode fixed point. Every role,
  owner, opaque source-schema routing address, and authority ceiling is bound.
  It declares no downstream schema dependencies; consumers resolve role meaning
  through the registry rather than constructing recursive identity edges.
- EulerScientificContract binds exact generic ContextOfUse canonical bytes and
  hash, an Euler-only population/risk/source extension, the nine-claim DAG, the
  seven mandatory core no-claims, and both the exact owner-matrix canonical
  bytes and its independently derived identity.
- EulerClaimKind is closed over numerical trajectory verification, calibrated
  reproduction, blind trajectory prediction, event/crossover prediction,
  qualitative direction, ranking, nonlinear optimum interval, energy-channel
  attribution, and mechanism attribution.
- EulerAcceptanceFamily distinguishes the decision semantics of those nine
  claims. The generic Context contains only aggregate decision QoIs. Detailed
  observables and their authority classes belong to a separately versioned,
  Context-bound measurement/score artifact owned by the downstream
  QoI/event/scoring contract. That contract (bead t6314.1.3) owns the
  measurement equations, covariance treatment, windows and filters,
  event classification, censoring, tie rules, multiplicity control, and exact
  formulas that derive aggregate Context QoIs. Its content-bound derivation
  receipt must bind the exact Context identity, complete detailed-observable
  registry, design-set identity, claim, and aggregate-QoI scoring scope. Its
  admission checker must cross-check the receipt's design-set identity against
  the consuming packet. This v1 leaf reserves the exact route and binds a
  nonzero receipt identity, but cannot construct, interpret, admit, or execute
  that artifact; its absence or non-admission prevents physical promotion.
- Numeric v1 QoIs whose names already denote a discrepancy, absolute error,
  interval width, or residual use `ClosedRange` as a predicate on that derived
  QoI itself. Non-negative magnitudes have a zero lower bound; the signed
  work-energy residual has symmetric bounds. `AbsoluteErrorAtMost` and
  `RelativeErrorAtMost` would instead mean a second discrepancy *in* the named
  QoI and are therefore deliberately not used here. This structural contract
  does not define or evaluate the downstream observation and scoring recipes
  that derive those QoIs.
- Event/crossover and nonlinear-optimum claims each require two aggregate QoIs,
  not merely a time error or interval width. `event-class-disposition` records
  whether the preregistered terminal/crossover class matched, while
  `event-time-error` bounds timing conditional on the downstream event rule.
  `optimum-containment-disposition` records whether the exact downstream score
  places the preregistered optimum inside the retained interval, while
  `optimum-interval-width` bounds its width. The categorical companion prevents
  a small numeric value from silently standing in for the required event class
  or optimum-containment decision.
- EvidenceRequirement is a categorical policy, not a score. Each claim kind
  has an exact minimum set and exact QoI set under
  EULER_CLAIM_POLICY_SCHEMA_VERSION.
- HypothesisSource is hypothesis-only by type. Its hash identifies the source
  declaration, not the source's authenticity, custody, truth, or validation
  status. An evidence packet cannot reuse that hash to satisfy an evidence
  role.
- EvidenceRecord binds one nonzero referenced artifact hash to one contract,
  claim, QoI set, evidence role, generic V&V container kind, referenced
  schema-admission receipt hash, caller-declared access class, authority
  declaration, and independence declaration. `DeclaredEvidenceAccessClass` is
  policy metadata, not the generic `fs_evidence::vv::EvidencePartition`; it
  does not prove observation membership or blind release. Numerical roles
  carry a declaration-only Color;
  categorical process roles carry a receipt hash. Neither is an admitted
  authority object. Every logical hash across the artifact, schema-admission,
  and categorical role-receipt slots is globally unique, and a declared
  hypothesis-source hash is forbidden in every slot. Reusing or laundering a
  hash refuses because no v1 typed composite-evidence receipt exists.
- A declaration-only Validated color may contain at most
  `MAX_VALIDITY_DOMAIN_AXES` (32) regime axes and at most
  `MAX_VALIDITY_DOMAIN_CANONICAL_BYTES` (8 KiB) for the canonical axis-count
  and bounds payload. Admission computes that exact shared Color-v2 wire size
  with checked arithmetic before color serialization. These are local resource
  ceilings, not evidence that 32 axes or 8 KiB is a scientifically adequate
  model of the physical validity domain; callers must refuse or use a future
  versioned external artifact rather than truncate a richer regime.
- Validated physical coverage is conjunctive: every numeric Context axis must
  be present in the ValidityDomain and covered by the exact point, and every
  bound axis carried by the ValidityDomain (including axes richer than the
  Context) must be present, finite, and in range at the exact point. Bounds are
  inclusive. Unchecked extra regime axes cannot silently widen validation.
- ClaimEvidencePacket carries explicit units, seed status, resource budgets,
  a campaign-anchor applicability point, a nonzero content identity for the
  complete comparison/configuration/search design set, a distinct nonzero
  content identity for the downstream aggregate-QoI derivation receipt,
  target-fitting state, no-claim acceptance, expected protocol result, and a
  caller-reported positive/negative/inconclusive scientific disposition. The
  single applicability point is an anchor for Context-domain checking; it does
  not enumerate or identify the relational design set. `ProtocolBudget` names
  its dimensionless computational control field
  `normalized_accuracy_limit`; that field is not a QoI, claim threshold,
  observed discrepancy, or criterion result. The packet contains identities
  and declarations, not raw target observations or an admitted
  criterion-evaluation artifact.
- ClaimPolicyAssessment returns one of ReferenceCompleteCandidate,
  DemotedCandidate, RetainedTerminal, or Refused. “Reference complete” means
  only that this local state machine found no missing/weak declared reference
  and received all exact direct-DAG prerequisite receipt identities. The
  referenced generic artifact, whole-case schema-admission receipt,
  ObservationSelection, blind release, producer, and independence declaration
  remain unresolved and unre-admitted. This is not evaluated criterion
  satisfaction, physical validity, mechanism truth, maturity, runtime
  admission, or release authority.
- ContractCheckReceipt separately re-decodes the generic Context, Euler graph,
  and whole composite transport, recomputes their hashes, and compares them to
  three literal-frozen v1 digest anchors. This detects builder drift unless the
  anchors are deliberately changed too; no independent review provenance is
  claimed.
- StructurallyAdmittedEulerContract is constructible only from a passing exact-
  contract check. It deliberately exposes no conversion to fs-govern authority
  types.
- ClaimPolicyAssessmentLog is a deterministic, bounded, policy-redacted JSON
  Lines record for one local assessment. It deliberately does not claim the
  campaign-wide evidence-event/progress/diagnostic/retained-log namespace owned
  by the later evidence-log contract. It binds
  protocol, case, packet, contract, claim, exact design-set identity,
  aggregate-QoI derivation-receipt identity, units, seed, budgets,
  expected/observed disposition, caller-reported scientific outcome, exact
  campaign-anchor applicability state, first evaluator divergence,
  authority/no-claim state,
  role-and-slot-labelled logical artifact identities, and a focused checker-
  smoke reproduction command. Every evidence row retains its artifact and
  schema-admission-receipt identities and the bounded source id/schema/kind of
  every referenced evidence row; structural/process rows additionally
  retain their categorical role-receipt identity, while numerical rows omit
  that slot because it is inapplicable. Prerequisite identities are labelled
  with the opaque prerequisite-receipt routing address obtained from the
  owner-matrix registry, avoiding a recursive log/prerequisite schema edge.
  The design set and derivation receipt are each retained once in the relative-
  artifact registry; the latter is labelled with its owner-matrix route.
  “Redacted” here means that the schema has no raw
  observation payload, source-prose, host-path, environment-dump, or artifact-
  byte field; stable case/source IDs and structured protocol metadata are
  intentionally retained and are not anonymous. This is not a secret scrubber:
  callers must not encode secrets or personal data inside machine identifiers.
  V1 does not persist or resolve the referenced artifacts themselves, so that
  command verifies the log/receipt contract and cannot replay the named physical
  or numerical case.
- `ports` is an Euler-local composition and accounting surface. A
  `PortDeclaration` binds one typed channel (gravity, normal, tangential,
  rolling/contour/spin, impact, base, exterior gas, or gas film) to an explicit
  effort/flow kind, active/inactive/unavailable state, law and source
  identities, canonical surface pair, patch interval, clock interval, and
  generalized-coordinate basis/frame/sign binding. `EulerPortRegistry` admits
  only deterministic, non-ambiguous ownership. `EulerEnergyLedger` records
  caller-supplied finite energy changes exactly once by contribution identity,
  with identity-bound checkpoints and rollback.

## Invariants

1. V1 contains exactly all nine claim kinds. Each kind carries the exact stable
   claim ID, acceptance family, QoI set, and minimum evidence set defined by
   claim-policy v1.
2. Every claim QoI exists in the exact generic Context, and every Context QoI
   is covered by at least one claim. Every QoI unit and every numeric
   applicability-axis unit occurs in the Context Five-Explicits header.
3. Physical applicability has both numeric and categorical axes. The
   Euler-declared observation frame must be admitted by the exact
   observation-frame categorical axis. V1 refuses, rather than silently
   extrapolating, outside the domain.
4. Claim dependencies are acyclic, have existing endpoints, are globally
   unique by endpoint pair, and cannot relabel one edge as both calibration and
   validation use. Evidence-gap IDs are globally unique.
5. The seven CORE_NO_CLAIMS statements are mandatory. Additional statements
   may narrow a constructed contract up to the exact public and decoder limit
   `MAX_EULER_NO_CLAIMS` (72). A 73rd row refuses before identity publication.
   A narrowed contract is not the exact frozen instance and fails
   check_frozen_contract.
6. Generic owner rows use exact crate, opaque source-schema routing address,
   and closed authority ceiling values. Caller prose cannot widen them. The
   complete registry has its own schema version, canonical transport, identity,
   byte cap, fixed-point decoder, append-only role tags, and no predecessor
   migration. Adding the aggregate-QoI route preserves the pre-existing owner-
   matrix registry tag 14 and appends the new role as tag 15. The composite
   contract embeds both its bytes and matching identity.
7. Canonical collection order is explicit. Graph and contract decoders enforce
   magic, version, bounds, known tags, canonical ordering, UTF-8/text rules,
   complete consumption, and encode/decode fixed points.
8. Context, graph, owner matrix, composite contract, hypothesis declarations,
   packets, prerequisite receipts, assessments, checker receipts, assessment
   logs, and downstream aggregate-QoI derivation receipts use distinct BLAKE3
   domains or separately routed schemas. Their identities are never
   interchangeable.
9. Physical numerical declarations require a structurally valid Validated
   color whose complete numeric regime and every numeric Context axis cover the
   exact case point; numerical
   verification declarations require Verified. Categorical process facts use
   a nonzero receipt reference instead of a Color. These checks validate local
   payload shape only: they do not re-admit the referenced generic V&V case,
   authenticate a producer, or turn a plain Color into AdmittedColor. Before
   shared Color serialization, every Validated regime is preflighted against
   the exact 32-axis and 8-KiB v1 limits with checked length arithmetic.
10. Calibration evidence cannot satisfy blind evidence. Target fitting
    terminally refuses every claim whose kind forbids it. Numerical
    verification cannot satisfy physical validation. Energy closure alone
    cannot satisfy mechanism discrimination.
11. Negative and inconclusive scientific outcomes are retained as terminal
    non-promotions when evidence is complete. They are never converted to
    missing data, deleted, or reported as positive candidates.
12. One logical hash cannot occupy two evidence slots of any kind in v1, and a
    hypothesis-source declaration hash cannot occupy an artifact, schema, or
    role-receipt slot. This prevents relabeling code proof as physical evidence,
    laundering a transcript through a receipt-shaped field, or treating energy
    closure as mechanism discrimination.
13. No constructor accepts an all-zero content identity. A packet's design-set
    and aggregate-QoI derivation-receipt identities must also be distinct from
    one another and from every evidence artifact/schema/role-receipt slot.
14. A dependent claim consumes every exact direct claim-DAG edge. A receipt
    binds contract, prerequisite, dependent, EvidenceUse, source packet,
    source assessment, exact design-set identity, and campaign-anchor
    applicability point; indirect ancestry, relabelled use, duplicates,
    missing edges, a different design set, or a different anchor point refuse.
    Assessment accepts at most `MAX_PREREQUISITE_RECEIPTS` (nine in v1), so the
    direct-receipt slice and retained reason set are bounded before sorting.
15. The assessment and log retain only a reported scientific disposition and
    structurally bind the unreadmitted aggregate-QoI derivation receipt.
    Acceptance-criterion execution is explicitly deferred to that downstream
    preregistered-analysis artifact and cannot be inferred from either value.
16. `negative-result-erasure` maps to the explicit
    `retain-as-terminal-non-promotion` decision alternative, never to candidate
    review. This response is part of the frozen risk registry and contract
    identity.
17. Active ports that overlap in canonical surface pair, patch region, clock
    interval, and generalized-coordinate identity refuse unless both carry the
    same exact additive decomposition receipt covering both named ports. Surface
    order and coordinate-frame/sign disagreement cannot create distinct
    ownership domains.
18. The ledger accepts a contribution only once, only from an active declared
    port of the matching channel, and only at a timestamp in that port's
    half-open interval. Entries are canonicalized by clock, tick, and stable
    contribution identity; invalid entries and invalid checkpoints leave the
    visible ledger unchanged. Kinetic, potential, recoverable, and numerical
    terms are signed; dissipated, heat, and unresolved magnitudes are finite
    and non-negative.

## Error model

All public fallible operations return ContractError with a stable machine code
and bounded actionable detail. Structurally required collections refuse empty,
duplicate, or oversized inputs. An evidence packet may intentionally carry zero
evidence rows so assessment can retain an exact `missing-evidence` refusal; an
empty packet can never yield a non-refused disposition. Construction also
refuses malformed text; invalid versions; dangling references; cycles; owner or
routing-address substitution; unit/frame/QoI mismatch; missing core no-claims;
more than 72 no-claim rows; zero hashes; design-set/derivation cross-role
aliasing; cross-role evidence aliasing; and
malformed color payloads. A maximum-plus-one validity axis refuses with
`EulerProtocolValidityDomainCardinality`; a maximum-plus-one canonical regime
payload or checked-length overflow refuses with
`EulerProtocolValidityDomainTooLarge`, before color serialization.

Assessment distinguishes hard refusals from weakest-evidence demotions without
computing a numeric authority score. Reasons are canonicalized and the first
divergence is retained in the log. A transported passed boolean is never
trusted: the contract-check receipt decoder enforces its result/issue invariant,
while `verify_subject` re-runs the structural checker, re-decodes and recomputes
the subject, compares it to literal-frozen digests, and demands exact equality
with that fresh receipt. A passing checker receipt establishes only exact
structural equality.

The scientific contract and owner-matrix registry are independently versioned.
Each v1 has no predecessor: its migration_policy(1) accepts the current schema,
while version 0 and every other value refuse. Protocol v1 has the same explicit
no-predecessor rule through `protocol_migration_policy`; the contract-check
receipt has a strict fixed-point decoder and the assessment log has a bounded
complete canonical reader. That log reader admits exactly the v1 field order,
JSON types, closed values, canonical escaping, resource bounds, and local
cross-field bindings that are decidable from the line alone. In particular, it
requires `contract_identity` to equal `FROZEN_CONTRACT_IDENTITY_HEX`
(`e95ae98859836b49370bc0a75749f7c6687cd1552a73ae8177fcfafbcb3d5e60`),
retains the packet's separately declared `packet_contract_identity`, and binds
an inequality between those two identities bidirectionally to the exact
`contract-identity-mismatch` reason and a refused disposition,
parses the closed evidence-source and relative-artifact row grammars, requires
each retained source role's artifact and schema-admission slots, and enforces
the role-receipt slot shape against the observed authority class carried by the
closed weak-authority reason grammar (or the role's required class when no such
reason exists). Each retained evidence-slot hash is bound bidirectionally to
the corresponding frozen hypothesis-source collision state: a collision
requires its exact slot-specific reason and refusal, while a noncollision
forbids that reason. Non-refused rows retain both the exact required evidence roles and the
exact per-claim prerequisite receipt-row counts implied by incoming frozen-DAG
edges. An absent expected prerequisite claim row requires its exact
`missing-prerequisite-receipt` edge reason and a refused disposition;
receipt-specific reasons that name a claim require that claim's retained row,
while a malformed-receipt reason requires at least one retained prerequisite
row and admits only the prerequisite receipt verifier's publicly producible
identity-mismatch code. Each log also requires exactly one top-level and relative-artifact binding
for the packet's design set and downstream aggregate-QoI derivation receipt,
refuses aliasing them with each other or any evidence slot, and admits the
closed `prerequisite-design-set-mismatch` reason only with a retained receipt
row and refusal/missing-edge semantics. Reasons contradicted by line-local
frozen claim or accepted source kind/schema state are rejected, as are all
reasons outside the closed v1 evaluator grammar. It does not treat a
version/domain envelope as sufficient. It does not re-evaluate or authenticate
the absent packet, prerequisites, referenced artifacts, producer, or scientific
claim.
Packet, prerequisite, and assessment
canonical byte functions are exact identity preimages, not promised public
rehydration transports in this leaf. A future version must define an explicit
semantic migration and receipt or continue refusing predecessor bytes. It must
never reinterpret v1 approximately.

The local port surface refuses malformed surface/patch/time domains, duplicate
port or contribution identities, active ownership overlap, unproved additive
overlap, unavailable/inactive or channel-mismatched energy records,
out-of-window or out-of-order receipts, non-finite/negative constrained energy
terms, cumulative non-finiteness, and checkpoints from another ledger, registry,
or retained prefix. It has no implicit channel default or fallback.

## Determinism class

Contract and owner-matrix construction, canonical encoding, decoding, hashing,
graph checking, evidence assessment, reason ordering, and assessment logs are
deterministic for the same exact input bytes and dependency versions. Set-like
inputs are sorted by explicit schema tags or stable machine identities. Direct
prerequisite receipts are ordered by their complete semantic tuples before any
refusal is observed, so caller slice permutations cannot change the retained
first divergence, log identity, or assessment identity. Floating criteria,
applicability bounds, and the dimensionless `normalized_accuracy_limit`
computational budget retain their canonical IEEE semantics from the owning
schemas; the local budget constructor normalizes signed zero. The packet encoder
also canonicalizes either sign of an applicability-point zero to positive zero;
that exact nonsemantic bit exclusion
is declared and mutation-tested in the packet identity schema.

The crate makes no cross-ISA claim about later floating-point simulation.
Deterministic metadata proof is not deterministic solver proof and neither is
physical validation.

Euler-local port declarations are sorted by stable declaration identity and
their remaining domain tuple; additive receipts sort contributor identities;
energy receipts are retained only in canonical clock/tick/identity order.
The resulting registry and accepted ledger order are deterministic for equal
inputs and dependency versions, but do not establish a deterministic contact or
multiphysics solver.

## Cancellation behavior

V1 contract operations are bounded metadata transformations and contain no
admitted long-running work, task spawning, blocking I/O, solver iteration, or
partial publication. Therefore no asynchronous cancellation surface is
exposed.
Caller-controlled cardinality limits are checked before sorting. Before the
packet writer allocates its canonical output buffer, it computes the complete
packet length with checked arithmetic, including every nested Color-v2 payload,
and refuses a result above `MAX_EVIDENCE_PACKET_BYTES`. Caller-owned validity
axes are likewise counted and their exact Color-v2 canonical length is checked
without cloning or serializing the regime.

The current campaign binary is synchronous and step-bounded, but exposes no
`Cx` or cancellation API. It therefore does not claim the project's general
hot-kernel cancellation invariant or conformance to it. A future cancellable
campaign must run under explicit asupersync scopes, poll at bounded tile/run
boundaries, drain cancellation, and publish no partial evidence or authority.

Port and ledger operations are bounded in-memory metadata operations and expose
no `Cx`. Record and rollback validate their complete next state before mutation,
so refusal does not publish a partial receipt or partial checkpoint rollback.

## Unsafe boundary

The crate forbids unsafe code. It uses only safe Rust and the canonical/hash
facilities of lower FrankenSim crates. There is no FFI, foreign numerical
kernel, memory mapping, platform-specific dispatch, or unchecked deserializer.

## Feature flags

There are no feature flags. The v1 scientific contract is a solid structural
spine, and the committed runner is a bounded numerical/software slice rather
than default-on experimental physics. Its profile contact, flexible-base, and
exterior records remain one-way bounded operators. Closed fluid/solid/contact
coupling, rarefied-gas or fluid-film claims, calibration, and inverse-model
capabilities remain out of scope and require their owning contracts and
evidence.

## Conformance tests

The G0/G3 contract suite must cover:

- a hand-maintained review oracle for all nine claims, exact QoIs, evidence
  minima, dependencies, seven no-claims, and all 15 artifact roles across four
  ownership surfaces;
- direct closed-range tests for each numeric derived QoI and exact categorical-
  criterion tests for all five event-class, qualitative-direction,
  configuration-ranking, optimum-containment, and rival-mechanism companions;
- whole Context, graph, owner-matrix, contract, and checker-receipt canonical
  identity and decode fixed points, plus assessment-log byte/identity determinism;
- every owner-matrix role's role tag, owner, opaque source-schema routing
  address, and authority ceiling, including independent version/domain/framing/
  ordering mutation and contract-embedded registry identity mismatch;
- every claim's synthetic caller-reported-positive structural fixture, removal
  of every required role, weakest-authority demotion, wrong declared access class,
  target-fitting guard, unre-admitted-reference candidate boundary, and
  negative/inconclusive terminal outcomes;
- evidence reuse across roles, stale contract/claim/QoI bindings,
  hypothesis-source substitution, calibration/blind leakage, and software-to-
  physical or energy-to-mechanism authority laundering;
- semantic identity mutation of decision, QoI/unit/criterion, applicability,
  user/population/environment/frame, risk, hypothesis source, claim policy,
  dependency, no-claim, and owner matrix;
- one-field identity-preimage batteries for every declared local schema field,
  including alternate domain/version/framing/order/endianness and upstream
  color/artifact-kind codec rules. Inadmissible alternate encodings prove hash
  binding only; constructor and hostile-decoder tests separately prove
  admission and refusal behavior;
- packet/assessment identity movement for exact design-set and aggregate-QoI
  derivation-receipt identities; zero/cross-role alias refusal; and direct-
  prerequisite receipt identity movement plus dependent refusal when the source
  design set differs despite an identical anchor point;
- canonical permutation invariance for valid and invalid direct-prerequisite
  receipts, stable first-divergence/log/assessment identities, and duplicate
  refusal;
- malformed magic/version/tag/UTF-8/length/truncation/trailing bytes, empty,
  maximum, and maximum-plus-one bounds, including the claim-graph transport
  cap, evidence-packet canonical byte cap (including a nested Color payload),
  public direct-prerequisite receipt cap, 72-row no-claim cap, assessment-log
  byte cap, and both validity-domain axis and exact canonical-byte caps;
- strict assessment-log refusal for missing, extra, duplicate, reordered, or
  mistyped fields; noncanonical integers and string escapes; unknown closed
  values; invalid identity widths; unsorted arrays; cross-field mismatches;
  trailing objects; and missing or embedded line terminators;
- strict assessment-log refusal for frozen-contract substitution; non-refused
  evidence-role or incoming frozen-DAG prerequisite-row mismatch; an absent
  expected prerequisite row without its exact missing-receipt reason and
  refusal; a receipt-specific reason without a line-decidable retained row; a
  hypothesis-source/evidence-slot collision not bound bidirectionally to its
  exact reason and refusal; and a claim/source mismatch reason contradicted by
  the retained line itself;
- strict assessment-log refusal for missing, zero, aliased, or route-substituted
  design-set/aggregate-QoI derivation bindings and for a design-set prerequisite
  mismatch not represented by its exact retained receipt and refusal grammar,
  including paired admission of the reachable receipt-identity mismatch and
  rejection of producer-unreachable malformed-receipt codes;
- missing/out-of-range extra validity axes and inclusive lower/upper coverage,
  exact three-slot evidence-log retention, owner-routed prerequisite labels,
  and the terminal negative-result response;
- compile/manifest direction through workspace layer/dependency gates and an
  API boundary in which a generic lookalike cannot enter the concrete Context
  constructor.
- port-registry permutation invariance; all eight typed channels; canonical
  action/reaction surface pairs; duplicate identities; partial temporal/patch
  overlap; contradictory coordinate frame/sign bindings; and additive receipt
  exact-domain/contributor refusal;
- exactly-once ledger receipt, inactive/unavailable/channel/timestamp refusal,
  canonical receipt-order refusal without partial mutation, checkpoint identity
  and prefix rollback, finite-energy validation, and explicit no-closure
  dispositions when a channel is unavailable or merely absent.

Focused tests are necessary software evidence only. Repository DSR, separately
custodied blind physical evidence, model comparison, and maturity promotion are
independent gates.

## Focused and closure-candidate runner

`scripts/ci/euler_disc_contract_e2e.sh` is the retained orchestration surface
for this leaf. Its default invocation is deliberately refused until the caller
declares where Cargo will execute. Examples:

```bash
FSIM_EULER_DISC_E2E_EXECUTOR=local \
FSIM_EULER_DISC_E2E_ALLOW_LOCAL=1 \
scripts/ci/euler_disc_contract_e2e.sh --profile focused

rch exec -- env \
CARGO_TARGET_DIR="${RCH_TARGET_BASE:-${TMPDIR:-/tmp}}/rch_target_frankensim_test" \
FSIM_EULER_DISC_E2E_EXECUTOR=rch \
scripts/ci/euler_disc_contract_e2e.sh --profile focused

scripts/ci/euler_disc_contract_e2e.sh --verify-bundle \
target/euler-disc-contract-e2e/<run-directory>
```

The focused profile may inspect a dirty candidate and therefore records source-
manifest authority as `NO_DATA`. The closure profile requires the complete root
index, tracked worktree, and non-ignored untracked set to be clean at both
bookends, requires the recorded HEAD to remain fixed, checks the constellation,
checks the whole HEAD-derived source manifest, and verifies exact membership of
the declared Euler paths. Even then it reports a clean-HEAD-bookended candidate,
not “exact HEAD”: bookends cannot disprove a transient concurrent edit. DSR is
still the repository-level authority and is not run or impersonated here.
Both profiles run the exact `xtask check-consolidation` structural gate, and
`consolidation-review.json` is an explicit snapshot/source-membership scope
path. A content change is therefore visible even when its Git porcelain state
remains `M` at both bookends.
Closure preflight exit 1 is the sole policy-refusal transition to `NO_DATA`; its
trace is exactly proof boundary, preflight attempt, and unevaluated source
manifest. Launch, supervision, timeout, signal, and every other preflight error
remain `FAIL`, stop the closure branch, and cannot mint manifest `NO_DATA`.
Each Git repository, index, object, and worktree-hash command is checked at its
real exit boundary: a non-1 infrastructure failure is retained as `FAIL`, not
collapsed into the exit-1 policy refusal.

The run derives one absolute monotonic aggregate deadline before setup; normal
and self-test supervisors share that deadline, and a lane with no remaining
launch budget refuses before spawning. Prepublication verification, the atomic
publication operation, and post-publication verification are independently
supervised against that same deadline. A success finalizer that finishes late,
times out, or fails must retract any renamed public candidate before a terminal
success is admitted; bounded process drain, retraction, and terminal failure
sealing may finish later. The success-deadline linearization point is the final
successful monotonic observation after post-publication semantic verification
and binding equality. Only signal-latch builtins and the atomic Bash success
commit follow it. This does not prove that `SEALED=1` itself executes before the
monotonic deadline instant under arbitrary scheduler preemption; no stronger
real-time commit claim is made. The retained hidden finalization logs make those
transitions diagnosable. Human-facing finalizer-log notices and failure replay
are queued until after a terminal bundle is sealed. A blocked or broken sink
during those deferred replays therefore cannot hold a renamed candidate in a
preterminal state or change the authoritative exit. Before any success rename,
the publisher opens the actual log-root directory with the same read-only,
directory-only, no-follow flags used by retraction and exercises the exact
same-filesystem atomic-exchange primitive relative to that descriptor. A
symlink final component, insufficient directory access, or unsupported kernel,
libc, or filesystem behavior therefore fails before publication. An explicit
retraction failure is a harness-integrity failure, not a success claim. Every
external lane also has its own monotonic deadline.
The supervisor re-samples signal state and monotonic time after process
inspection, after bounded output polling, and immediately before admitting
completion. Its no-Cargo self-tests deterministically delay that final
classification across a deadline and latch TERM at that boundary, proving that
an already-exited child cannot be admitted late or ahead of a consumed signal.

Each successful exchange preflight retains one hidden
`.retraction-exchange-preflight-*` sibling containing the two exchanged empty
directories, and records that exact path in the bounded publication-finalizer
log and deferred post-seal diagnostics. This retained capability evidence is
not part of the proof bundle or its byte budget. The harness never auto-deletes
it, so repeated runs accumulate these small trees; callers that require bounded
long-term storage must provide an isolated log root and an explicitly governed
artifact-lifecycle policy.

The supervisor creates a new session/process group and keeps the leader
unreaped while inspecting every live non-zombie member of the owned session,
drains combined output beyond the bounded retention cap, and uses TERM then
KILL with bounded grace. A descendant remains in scope if it changes its
process group with `setpgid`; cleanup signals every observed session group and
PID and independently checks that the session is empty. Success requires leader
exit, pipe EOF, and an empty owned session. A leader that leaves descendants is
a failure even if cleanup succeeds. Deliberately daemonized processes that
create or enter another session, power loss, PID-reuse races outside an observed
scan, and kernel failure are outside this containment claim.

Timeout and retained-byte settings are canonical positive decimal integers
validated outside Bash arithmetic. Lane and aggregate timeouts are each at most
604800 seconds. A lane's complete retained log is at most 64 MiB, all lane logs
together are at most 512 MiB, and normal scheduling reserves 256 KiB of that
aggregate for terminal diagnostics. Before every command launch, the producer
derives its remaining aggregate allowance and reduces the child-output
allowance without weakening the complete-log ceiling. Each provenance snapshot
is streamed through a 16 MiB cap; maximum-plus-one drains and fails without
publishing authoritative snapshot bytes. A snapshot enters the proof directory
and receives a seal hash if and only if its lane passes. Failed partial
candidates are retained outside the sealed bundle for diagnosis.

Each run is first written beneath a hidden candidate directory. It retains per-
lane `command.v1` metadata, canonical JSON argv, exact authority, exact log
locator, exact exit disposition, log byte count and SHA-256, and canonical
`supervisor-result.v1` metadata cross-checked against the verdict,
before/after snapshots, the exact nonterminal `verdicts-prefix.jsonl`, JSONL
verdicts, and one terminal proof seal. The complete candidate is strictly
verified before publication. The verifier emits a binding over the no-follow
root directory identity, exact relative inventory, file sizes, and verified
file digests. The publisher reopens that same identity, recomputes the binding
immediately before and after the same-parent rename, and the strict semantic
verifier runs again at the published path. Publication uses an atomic exclusive
no-replace operation (`renamex_np(RENAME_EXCL)` on Darwin or
`renameat2(RENAME_NOREPLACE)` on Linux) and fails closed on unsupported
platforms or an occupied destination. A detected gap or post-rename
mutation refuses publication or retracts the published path. Retraction opens
the expected bundle and an empty placeholder through no-follow directory
descriptors, checks both identities, and uses `renameatx_np(RENAME_SWAP)` on
Darwin or `renameat2(RENAME_EXCHANGE)` on Linux so the expected bundle lands
directly at a hidden retained sibling path. After a successful exchange and
identity check, the exact empty placeholder remains at the former public name
as a non-proof tombstone; it is never interpreted as a terminal bundle.
Identity mismatch is preserved rather than moved, while an exchange or post-
exchange integrity failure fails closed. This protocol assumes
ordinary custody of the log-root namespace. It does not claim to defeat a
continuously hostile actor that can replace names during retraction, nor does it
claim crash-durable rollback. An invalid staged success is retained separately
and a fresh minimal incomplete candidate is sealed, so contaminated success
bytes are never repurposed as terminal proof.
`summary.json` is byte-identical to the final `verdicts.jsonl` seal line.
`--verify-bundle` reruns no Cargo work: bounded, no-follow reads reject
oversized inputs, duplicate JSON keys, floating/nonfinite constants,
noncanonical bytes, unknown fields, unsafe paths, multiply linked files,
unreferenced files/directories, and record or directory-inventory overflow. It
also checks that the prefix file is byte-identical to the verdict stream before
the seal, the prefix hash, unique lanes/logs, exact ordered producer traces and
terminal transition/readiness matrices, exact no-claims, every log size/hash,
snapshot hashes, derived counts, and the byte-identical terminal seal. An
interrupt or incomplete outcome is admitted only when its terminal lane is last,
is `FAIL`, and has the exact authority, locator, and log bytes implied by the
terminal exit.

Verification is rooted at one no-follow directory descriptor. It inventories
entry identities and metadata before semantic reads, inventories again before
the final pass, reopens and rehashes every seal-derived file relative to that
descriptor, then performs a final inventory/root-identity check. This refuses
mutation observed during the invocation, including path replacement between
semantic validation and final rehash. “Verified” is descriptor-stable for that
completed invocation, not an immutability, lock, custody, crash-durability, or
after-return mutation claim. The publication rename is atomic on its current
filesystem; a reported parent-directory fsync failure explicitly leaves power-
loss durability unestablished. HUP, INT, and TERM are not masked across the
publication window. When retraction infrastructure remains available, a signal
observed by the bounded publisher either prevents candidate admission before
rename or moves an already-renamed candidate back under a hidden retained path.
The final `SEALED=1, SEALING=0` transition is one Bash arithmetic builtin, which
is the signal linearization point: a signal latched before that commit is
rechecked and revokes the candidate, while a signal processed after commit may
terminate the wrapper but does not rewrite an already sealed terminal success.
Retracted bytes are resealed as the exact `INTERRUPTED` outcome at a fresh public
name; the old name remains the empty tombstone. If placeholder allocation,
exchange, or an identity check fails after rename, the candidate may remain at
its former public-looking path. The failed transition leaves the wrapper's
`SEALED` state zero for that candidate and does not admit it; the directory's
already-written proof-seal bytes do not establish wrapper admission. The wrapper
returns harness-integrity code 125, and recovery may separately seal a fresh
`INCOMPLETE` bundle. No interrupted/success admission is reported for the
failed candidate; its path is retained failure evidence rather than admitted
proof. If a wrapper signal was concurrently latched, its number and transition
context remain in the retraction diagnostic even though integrity exit 125
takes precedence over an `INTERRUPTED` terminal status. A dedicated integrity-
recovery guard is established before the sealing critical section is released;
the signal trap may retain the first late signal but cannot re-enter interrupt
sealing, and the EXIT trap snapshots and clears that control latch only after
masking further wrapper signals. When the wrapper consumes a signal, it masks
that signal and clears the forwarding latch before launching any finalizer, so
a newly spawned supervisor cannot receive the same signal a second time. The
self-test-only active-supervisor readiness marker is likewise a one-shot
handshake: the harness captures the environment binding into shell-local state,
removes it from helper environments, and passes a pending marker only as
supervisor argv; the created marker file remains retained evidence. Before
finalization, the consumed active marker is replaced by a fresh shell-local
finalizer-marker path. The first finalizer creates that marker only after its
signal handlers are installed, then remains live until the wrapper has exercised
the forwarding latch and atomically creates the paired release marker with the
exact consumed signal number. The supervisor validates that no-follow payload,
then exclusively creates a PID-and-signal-bound acknowledgement; the wrapper
requires that acknowledgement before accepting the handshake. Readiness and
release share a 30-second supervisor bound, acknowledgement has a 10-second
bound, and bounded process-tree drain plus that handshake remain below the
60-second termination and 120-second enclosing signal-self-test budgets. Marker
payloads are completely written and synced into prepared regular files before
an exclusive same-directory hard link makes the final marker name observable;
the prepared evidence link is retained. The test arm and pending marker are
then consumed so later
verifier and publication helpers neither demand another one-shot marker nor
inherit the test hook. Armed-handler regressions cover TERM
in an ordinary lane plus dedicated HUP, INT, and TERM, require the direct
interrupted seal, and explicitly reject both staged-bundle fallback and an
`internal-incomplete` lane on these successful signal paths.

`--self-test` itself does not resolve Cargo. Each nested normal-harness fixture
uses a 20-second lane bound and 90-second aggregate run bound inside an exact
120-second finite outer containment budget. Those inner limits keep a
launched-but-never-started fixture interpreter from occupying the outer lane
for an hour; they do not change the deliberately injected half-second
success-deadline semantics and are not runtime or performance evidence. The complete recursive
candidate-corruption matrix has an exact 600-second bound: the four nested
fixture allowances plus one fixture-sized headroom interval. These are
containment budgets, not formal worst-case runtime bounds. Its 79-lane v1 matrix exercises
ordinary and nonzero exits, reserved-code child exits, child signals, spawn
failure, timeout/drain, bounded output loss,
runtime-output assertions that exclude escaped command/argv fixture text,
the exact safe log-cap boundary; exact-maximum, maximum-plus-one, and huge-text
numeric settings; exact and overflowing snapshot byte boundaries; snapshot
hash-read failure propagation; incremental aggregate-log exhaustion; exact
closure refusal versus real and injected
infrastructure failure; unexpected-supervisor cleanup with an independent
process-table check; a same-session `setpgid` escape attempt with an external
post-run session scan; separate assume-unchanged, skip-worktree, and
fsmonitor-valid Git concealment; simultaneous output truncation and timeout;
zero-test sentinel refusal; a synthetic structural gate that rejects an invalid
Euler consolidation disposition; a
`consolidation-review.json` mutation whose Git status remains `M` at both
bookends; deterministic aggregate-deadline expiry in the verified
prepublication window with success-publication refusal and later incomplete
sealing; independently supervised expiry during post-publication verification,
public-candidate retraction, and later incomplete sealing; exact focused,
closure `READY_FOR_DSR`, and `SELF_TEST_PASS` fixtures;
consistently rehashed authority, detail, log-locator, required-lane, and snapshot
mutations; hostile/malformed/oversized proof bundles; a maximum-plus-one
proof-bundle inventory that refuses before collecting or sorting entries;
wrong normal-lane authority
and argv metadata; supervisor-result, containment, writer-flag, retained-log-cap,
and aggregate-deadline contradictions; exact proof-boundary and source-manifest
body mutations; continued execution after a failed control lane; premature
snapshot authority; concurrent verifier and verifier-to-publication mutation;
prefix and complete-inventory refusals; hostile terminal-status matrices; a
final-component symlink log-root preflight refusal; deterministic post-exchange
public-tombstone substitution with zero success admission; prepublication
candidate corruption, destination-collision races, and signal injection both
before rename and on the pre-commit side of the final atomic `SEALED`
transition; fresh incomplete fallback; arbitrary-lane interrupted commit
boundaries; incomplete exits; post-rename fsync failure; and real HUP/INT/TERM
delivery both during publication and active lanes. Auxiliary fixture
repositories and failed snapshot candidates are retained as siblings,
never smuggled into the sealed proof directory, and all evidence is retained
instead of deleted. The shell `command.v1`, `verdict.v1`, and `proof-seal.v2`
formats are intentionally compatibility-free diagnostic harness evidence: they
must be consumed by the strict co-versioned verifier and are neither ledger
authority nor separately governed public identity schemas. The retained Rust
`reproduction_command` is only the non-vacuous
one-test checker smoke used to verify the log/receipt schema; it does not replay
the recorded case, resolve its logical artifacts, or establish physical
evidence.

## No-claim boundaries

The binding v1 statements are:

1. Transcript and publication sources generate hypotheses only; they are not
   validation evidence.
2. Fitting or selecting against protected target outcomes is calibrated
   reproduction, not emergent prediction.
3. Agreement in an exponent, event time, or stop time does not identify an
   energy-loss mechanism.
4. Geometric similarity does not establish dynamic similarity across scale,
   material, support, or environment.
5. Deterministic software verification does not establish physical validation.
6. A successful blind case is local to its exact declared Context of Use and
   applicability domain.
7. Negative and inconclusive results are retained terminal outcomes and are
   never erased or promoted.

Additionally, this crate does not claim that the video transcript is complete
or authoritative, that literature hypotheses are correct, that the chosen
domain is scientifically adequate, that the frozen tolerances are empirically
justified, that a Validated color authenticates laboratory custody, that the
32-axis or 8-KiB local validity-domain resource ceilings establish scientific
adequacy or license silent truncation, that the generic V&V container by itself
proves sample independence, multiplicity
control, energy closure, or rival-mechanism exclusion, that a logical artifact
identity has been persisted and replay-resolved, that caller-reported
scientific disposition is an executed acceptance criterion, that a nonzero
aggregate-QoI derivation-receipt hash means the downstream measurement/score
artifact exists or has been admitted, that a single anchor point identifies a
complete relational design set, that `normalized_accuracy_limit` is a
scientific acceptance threshold, that the
structural checker is an organizationally independent reviewer, or that any
counterintuitive Euler-disc effect will emerge from simulation. The all-claim
positive fixtures are synthetic state-machine tests only. Those questions
require later preregistration, exact measurement/custody artifacts, admitted
criterion-evaluation receipts, solver and model verification, blind
experiments, and candid comparison against negative or inconclusive outcomes.

The bounded campaign additionally makes no claim that its one-way snapshot
composition is closed multiphysics; that its encoded reduced decay exponents or
channel laws were identified from data; that any result ranks sharp versus
filleted edges, glass versus steel, rings versus discs, or other configurations;
or that it predicts a spin time. Its geometry, contact, base, decay, and
exterior records are numerical/software rungs, not experiment-, video-, or
Mould-backed evidence.

The Euler-local port registry does not implement any gravity, contact,
partial-slip, rolling/contour/spin, impact, base, exterior-gas, or gas-film
physics. Its structural additive receipt is not a signed decomposition proof;
it neither authenticates a source nor establishes action/reaction balance. Its
ledger neither derives energy from forces or impulses nor closes an energy
window: unavailable and undeclared channels always retain an explicit
no-closure boundary. It does not satisfy the blocked generic `PortSchema`,
manifest, constraint/impact-law, RATTLE/generalized-alpha/nonholonomic, or DSR
lanes.
