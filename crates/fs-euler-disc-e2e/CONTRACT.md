# fs-euler-disc-e2e contract

Status: [S] structural scientific-contract infrastructure plus a bounded,
deterministic numerical/software campaign slice. The slice composes real
squat-disc line/arc geometry and mass properties, a profile-native reduced
rigid-body/contact/base/gas trajectory, transactional finite-patch normal,
partial-slip, rolling, exterior-air, thin-gap gas-film, and reduced-base ports,
plus encoded reduced-decay ablations. The higher-fidelity ports now compose
atomically for bounded smooth-contact trajectory prefixes; impact and
separation still hand off through typed refusals. This is not an experiment or
video validation and does not support claims about Steve Mould's reported
observations.

## Purpose and layer

fs-euler-disc-e2e is an L6 HELM leaf for the Euler-disc emergent-prediction
flagship. Its first schema freezes what questions may eventually be asked
before geometry, solver, calibration, or experiment choices can bias those
questions.

The committed `euler_disc_campaign` binary is a deterministic JSONL runner for
the bounded slice. Its closed reduced lane evolves profile support, rigid-body
motion, unilateral contact, a one-mode base, rolling loss, and reduced gas drag
together. Opening and reimpact crossings are sought on four deterministic
subintervals, while terminal-inclination endpoint crossings are localized by
re-evolving the actual reduced step. Body, contact, base, rolling, and gas
channels share an explicit-midpoint stage, eliminating the former staggered
start-force/base-response lag. A separate Estimate-only
production composition exercises the higher-fidelity transactional adapters in
smooth contact, but the campaign does not yet use that composition as its
event-loop backend. Neither lane identifies a physical decay law.

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
records from diagnostics. The wrapper emits the original bounded component
records and manifest, followed by a v3 closed-lane block. That block contains
12 profile-native controlled cases (solid, fillet, chamfer, fixed-density and
sharp equal-mass rings, a matched 1-mm-outer-fillet equal-mass ring, symmetric
taper, scale, density, base, gas, and rolling ablations), a five-rung numerical
convergence record, a five-rung censor-aware matched-fillet equal-mass
ring/solid ordering record, a typed physical-calibration NO-DATA
record, and its own manifest. The closed cases use a declared 2/4/8/16/32 s
continuation policy. Any case still live at 32 s remains right-censored; that
bound is never reported as a spin duration. The finest predeclared 20/10/5
microsecond triplet owns the convergence decision while 80 and 40 microseconds
remain visible as pre-asymptotic sentinels. Both manifests are deterministic
digests, not validation certificates.

This contract deliberately records no output path, digest, or numerical value
until an actual retained run supplies them. The runner's receipts demonstrate
only the encoded numerical/software composition and its input/record checks.

The reduced event locator scans four fixed subintervals before bisection. An
empty event list does not certify absence of an open/reimpact excursion wholly
between adjacent scan nodes. The production smooth-contact helper
is restartable and atomic across its accepted adapters, but it is intentionally
not an impact/event owner. Thin-gap gas-film work has exact-once transactional
ownership and typed continuum/applicability refusals; it is not rarefied-gas or
resolved body-base FSI evidence.

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

`scientific-contract` is enabled by default and gates the `fs-ir`-backed
contract and protocol surfaces. `cinematic-render` is an opt-in [F] feature
that enables the `fs-render/tracer` dependency, animated scene bridge, and
focused render E2E without promoting the frontier tracer into the default
scientific-contract path. `render-checkpoint-ledger` implies
`cinematic-render` and adds the optional L6 `fs-ledger` persistence adapter;
it remains absent from the default graph. The committed runner remains a bounded
numerical/software slice rather than default-on experimental physics. Its
profile contact, flexible-base, and exterior records remain one-way bounded
operators. Closed fluid/solid/contact coupling, rarefied-gas or fluid-film
claims, calibration, and inverse-model capabilities remain out of scope and
require their owning contracts and evidence.

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

## Parameterized material specimen and contact ingress

`DiscProfileSpec` remains the geometry authority for solid, annular, filleted,
chamfered, and tapered axisymmetric specimens. `resolve_with_material_state`
derives density, mass, centroid, and inertia from the same complete
`IsotropicSolidStatePoint` later consumed by contact adapters;
its identity binds the exact chart, mass properties, material card, physical
query point, selected property values, and usage receipts. Equal density does
not make copper, steel, gold, wood, ruby, or another named state interchangeable.

Structural vibration instead consumes `ResolvedElasticDiscProfile`, which
binds the same geometry/mass authority to either the minimal isotropic elastic
state or a full oriented orthotropic tensor. The material-axis rotation and its
evidence identity are load-bearing inputs. They flow into the actual 3-D
stiffness matrix and therefore into modal frequencies, shapes, BEM radiation,
and sound; neither chemistry names nor an isotropic fallback select them.
The cinematic specimen accepts both symmetry families as numerical data and
binds the resolved mechanical and visible-optical bundles from one immutable
material card and thermodynamic query point. Thus anisotropy changes the real
modal/acoustic solve, while its chemistry label remains identity metadata only.

The normal-contact ingress accepts a `BoundNormalContactModel`: both ordered
bulk states, the complete `InterfaceSystemCard`, its resolved adhesion datum,
and one executable constitutive-model card must share one exact physical state
point. The model card—not the Euler fixture—selects rate-independent Hertz
versus the point-contact Hunt--Crossley rung and owns damping, characteristic
rate, applicability ratios, and finite temperature validity. Local profile
geometry independently selects a sphere/plane or true two-curvature elliptic
Hertz coefficient, so dissipation cannot replace an elliptic rim by a sphere.
The Euler adapter
supplies only geometry-owned half-space/layer extents and propagated
uncertainty. It cannot retype density, modulus, Poisson ratio, yield stress,
temperature, adhesion, damping, rate limits, surface order, or law identity.
This is a generic card-to-contact composition; it contains no Euler outcome
targets or material-name presets.

`ResolvedPhaseDiscProfile` binds the reference profile and mass properties to
the generic `fs-material` specific-enthalpy phase state. A fully solid state can
enter fixed-topology mechanics only when material-card identity, temperature,
and density exactly match the independently resolved elastic state. The first
nonzero liquid fraction refuses that rung as `EvolvingPhaseRequired`; stale
rigid geometry, modes, sound, or optics cannot silently survive the transition.

`ResolvedDiscProfile::thermal_geometry` derives complete boundary area and the
whole-body Biot length `V/A` from that same exact line/arc chart and resolved
volume. Its geometry identity is material-independent, while a later thermal
body identity must additionally bind mass, phase curve, conductivity,
emissivity, and boundary transfer. It is a valid ingress to a whole-boundary
isothermal rung only; partially insulated or spatially varying boundaries must
escalate to a partitioned thermal mesh rather than substituting a hand-entered
effective area. The cinematic critique manifest records these measures, but
their presence alone does not claim that thermal evolution is coupled into the
current source-bound trajectory.

`ResolvedPhaseDiscProfile::bind_lumped_enthalpy_march` is the reduced
thermal-to-geometry handoff. It accepts the generic `fs-conduction` implicit
enthalpy march only when that march's body uses the specimen's bit-exact mass,
complete profile-derived area, `V/A` length, phase curve, and initial state.
The preferred body transport queries temperature-dependent conductivity and
hemispherical emissivity from that phase curve's same immutable material card;
its query receipts and temperature grid participate in the body identity.
Every accepted thermal boundary is then converted through
`mass_conserving_state`: reference mass is held invariant, the volume required
by the current equilibrium density is computed, and one of three explicit
geometry regimes is returned: unchanged reference geometry, solid
thermomechanical update required, or evolving free surface required. This
prevents a convenient surrogate cylinder from driving another specimen. The
volume is a conservation constraint, not by itself an isotropic scale or liquid
shape; downstream deformation/remeshing must satisfy it before updated
mechanics, acoustics, or optics can be admitted.

For a still-solid, homogeneous, unconstrained, isothermal body,
`UniformIsotropicFreeExpansionLaw` provides one explicit bounded constitutive
rung from required volume to geometry. It applies the cube-root volume scale to
every profile length (including bores, fillets, chamfers, and taper),
reintegrates mass and inertia at the evolved equilibrium density, and checks
mass/volume closure independently. Its authority identity and strain-validity
ceiling are content-bound; material names never select the law. The result is
geometry authority only: temperature-dependent elasticity, contact, damping,
optics, and structural modes must still be re-resolved before use. Its contact
and structural-acoustic binders require the re-queried solid/elastic state to
match the phase state's exact card, temperature, and density before reusing the
evolved chart. Any liquid
fraction, anisotropic or constrained expansion, spatial temperature field, or
excessive strain refuses this rung and requires a higher-fidelity solver.

The admitted reduced rung evolves whole-body enthalpy under uniform ambient
convection, radiation, and internal power only while its Biot gate supports an
isothermal body. It does not resolve spatial thermal gradients, deform or
remesh a mushy/liquid specimen, recompute structural modes after topology
change, or derive phase-dependent optical/acoustic properties. Those effects
require the generic spatial thermal, thermomechanical/free-surface,
structural-acoustic, and optical couplings; the cinematic fixture must not
simulate them with scripted geometry or material-name switches.

## Animation-grade Euler render trajectory v3

`render_trajectory` defines the accepted public state boundary shared by later
image and sound pipelines. It exists because `CoupledSample` alone cannot place
or orient a disc: the complete `RigidBodyState` otherwise survives only in the
restart checkpoint. The v3 trajectory retains, at every accepted time, the
center-of-mass pose, canonical body-to-world unit quaternion, world linear
momentum, principal-body angular momentum bound to exact mass properties,
symmetry-axis diagnostic, contact branch and geometry, localized contact
transitions, the exact start and contact-activity flag of each accepted
interval, complete one-mode base displacement/velocity, channel
wrenches/work, total energy/defect, redundant Euler QoIs, and an explicit final
terminal, censor, or numerical-refusal disposition. A localized inclination
event is required exactly for the corresponding terminal disposition.

Top-level metadata binds the resolved specimen profile and chart, mass
properties, initial state, base and physics models, full configuration,
the nominal reduced-base frame, explicit per-channel availability, restart
fingerprint, fixed timestep, producer version, applicability, and mandatory
no-claims. V3 admits only right-handed Cartesian `+z`-up world
coordinates and SI/radian units. Those declarations are repeated on each raw
sample and must match. Quaternion inputs must already be finite unit
quaternions; admission canonicalizes the `q/-q` double cover. Times are finite,
non-negative, and strictly increasing; each interval start is exact, must equal
the preceding endpoint after the first sample, and may equal its endpoint only
for an interval-data-free initial point. A positive interval endpoint may not
exceed `interval_start + metadata.timestep_s` advanced by 32 nonnegative
binary64 ULPs. This admits producer-shaped addition at large absolute clocks
without a cancellation-prone endpoint subtraction; it is an admission
tolerance only and does not rewrite either retained time. Contact geometry is
present exactly on
the closed branch, contact normals and redundant symmetry axes are unit vectors,
localized transition times are ordered inside the retained interval and
alternate branches, and redundant QoIs must agree with the authoritative state
and bound mass properties. A channel declared unavailable must contain an exact
zero payload, while an available zero payload remains distinguishable. The
interval contact-activity flag cannot be inferred from force magnitude because
a localized reimpact root may be active at zero penetration and zero force. On
a positive-duration interval, an active-contact interval that ends open must
retain an opening transition; an inactive interval may end closed only at an
exact endpoint reimpact. Across retained adjacent samples, an interval with no
transition must preserve the preceding endpoint branch, while its first
transition must be reimpact from open or opening from closed. The first
positive-duration sample remains explicit preroll: because its segment-start
branch is not retained, v3 makes no claim that its first transition is bound to
that unavailable state.

The reduced coupled runner now publishes those fields only after an accepted
macro interval or an accepted localized terminal/contact boundary. `CoupledRun`
retains the configuration's original rigid/base state, exact mass properties,
and macro timestep even for resumed segments, so trajectory admission can bind
metadata without reconstructing the initializer. Every published center-of-
mass velocity is checked against the retained momentum and those exact mass
properties. A prohibited reimpact is retained at its localized root with the
reimpact transition, closed post-root branch, updated event count, checkpoint,
energy ledger, and numerical-refusal disposition; no positive-duration closed
mechanics is evolved past that boundary. Restarting that terminal checkpoint
publishes no duplicate sample. Profile-backed runs poll their execution scope
before setup and immediately before each checkpoint/sample commit; cancellation
at either boundary returns a typed refusal and publishes no partial `CoupledRun`.

The source-bound reduced-decay render bridge v2 retains its grounded no-slip
kinematics and declares a separately available normal-load scalar while the
aggregate contact wrench/work channel remains unavailable. Its initial
zero-duration point remains contact-inactive with an exact zero interval
payload. For every later interval it reconstructs the duration-mean normal
reaction required by endpoint vertical impulse balance under the reduced
model's gravity-plus-support vertical-force closure,
`N_bar = m * (g + (v_z1 - v_z0) / dt)`, and refuses a non-finite or negative
reaction. This scalar has no associated authoritative torque or work. It is a
kinematically implied interval quantity, not a resolved contact patch,
subinterval force history, angular-impulse solution, tangential traction,
measured force, deformable-contact model, or acoustic-radiation model.

The transactional production-coupling bridge publishes accepted smooth-contact
prefixes and also exposes one shared-checkpoint open/contact compliant driver.
Its profile-native initializer resolves the actual support point and both
principal curvatures, evaluates the caller-selected normal law, derives the
finite normal patch used to initialize the caller-selected tangential law,
starts rolling history, and seals the explicitly selected gas channel into one
checkpoint. It does not infer any constitutive coefficient or contact radius
from a specimen name. Product drivers therefore share the same initialization
boundary instead of recreating synthetic normal-patch state by hand.
Product drivers may also initialize rolling motion from explicitly declared
precession and spin rates through the same resolved profile. That generic
initializer derives mass, principal inertia, ground support, linear momentum,
and angular momentum from the chart and enforces zero material-point velocity
at the initial contact. The optional small-angle helper is only one analytical
way to choose those rates; neither the generic initializer nor the production
coupling model imposes that approximation.
For long-form product trajectories, the moving-contact modal base port consumes
the certified rectangular plate basis directly. It projects each actual
base-frame contact point through a C1 rectangular Hermite reconstruction of
the plate solve's retained DKT displacement and slope degrees of freedom,
advances the same
mass-normalized modes used by structural acoustics, publishes local surface
displacement and velocity back to contact kinematics, and retains work,
stored-energy, viscous-loss, and closure accounting. Checkpoints keep a fixed
size step-lineage root instead of cloning every prior step identity. The
production coupling owner accepts either this resolved moving-contact backend
or the older bounded one-mode estimate without changing disc/contact laws.
This smooths the numerical moving-load field but does not yet integrate the
finite normal-pressure distribution over the plate modes: the structural load
is still the admitted point resultant at the resolved contact location. The
finite-gap solver's pressure field and pressure moments remain authoritative
for contact admission, while distributed structural projection is an explicit
outstanding fidelity rung rather than a claimed property of this path.
Every contact sample comes from the exact `fs-mbd`
`StepReceipt` and retains its accepted rigid state, profile-native contact
feature and point, normal, selected structural-base endpoint, and disc mechanical-energy
diagnostic. Because the force law is evaluated at the interval start while a
render sample describes the endpoint, the bridge independently re-queries the
resolved profile support at the accepted endpoint pose and base displacement;
an endpoint beyond the plane is refused as an unlocalized opening instead of
being mislabeled closed. It checks checkpoint integrity plus exact state, time,
base, model, and specimen lineage before admission. Its normal scalar is tagged
`AppliedSubstepZeroOrderHold`: this is the normal-law evaluation that the
transactional mechanics step applies as a constant world force over the
accepted interval. It is therefore the exact duration mean of the discretized
forcing and may drive the normal-force excitation seam, but its approximation
to the continuously varying physical force and its usable bandwidth remain
timestep-dependent convergence questions. The bridge leaves all aggregate
channel-work availability false until the production composition publishes one
shared cross-channel work ledger. During an open interval, normal, tangential,
and rolling state is retained without fabricating a zero-force contact port;
the disc advances under gravity plus the selected gas law and the support
advances with exactly zero applied contact load. A later approaching or
impact-candidate state may re-enter the same finite-patch law only through the
explicit time-resolved compliant regime and its declared rate, pressure,
strain, temperature, and geometry applicability envelope. Branch changes carry
fixed-grid time brackets, not exact event-time claims; timestep refinement must
shrink those brackets. No restitution coefficient or synthetic impact impulse
is inferred. Its energy defect follows the actual
`fs-mbd` convention of world-fixed force and body-fixed torque, but remains a
disc-only residual rather than total disc/base/contact/gas closure. Reaching
the requested step budget is horizon censoring; a source refusal remains
an explicit backend-specific numerical refusal. The older smooth-prefix render
bridge never invents separation, impact, reimpact, or terminal continuation
after the source stops; the event-aware trajectory requires its own
render/control-stream admission.

An optional production substep may additionally carry two owned material-frame
surface-height traces and their actual path coordinates/speeds. The accepted
Hertz footprint filters those traces through `fs-tribo::surface_excitation`,
and the normal law's consistent tangent converts the admitted small height into
one action/reaction force perturbation. That perturbation is included in both
the `fs-mbd` wrench and the moving-base load before the step commits, and its
receipt remains available to structural acoustics. It is therefore not a
post-hoc audio oscillator. This rung is explicitly first-order: it preserves
the nominal action line and patch used by tangential/rolling laws and refuses
outside the caller's bounded height/approach fraction or if the perturbation
would open contact. Large topography requires nonlinear contact re-resolution;
no material name, audible frequency, or renderer setting selects this channel.

The cinematic configuration may supply those periodic tracks either as
explicit Fourier coefficients or as an `fs-tribo` band-limited self-affine
profile. The latter retains explicit RMS height, one-dimensional Hurst
exponent, inclusive spatial-cycle cutoffs, phase seed, source identity, and
authority ceiling. The seed changes spatial phase only; actual contact-path
kinematics convert wavelength to time, after which the same Hertz footprint,
normal tangent, structural dynamics, and acoustic radiation operate unchanged.
The built-in critique specimen defaults to an explicitly ideal-smooth empty
spectrum rather than inventing unmeasured roughness. A caller may supply an
`Estimated` or measured spectrum through the same parameterized surface path,
but it is not profilometry, a material-name finish preset, an areal roughness
field, or permission to amplify roughness until a desired soundtrack appears.
A realistically rougher specimen that violates the height/approach limit must
escalate to nonlinear rough-contact resolution rather than reuse this tangent
perturbation.

Localization brackets describe evaluated uncertainty bounds, not additional
accepted states. A terminal event's retained time must exactly equal its final
sample time, and a `ReimpactLimitExceeded` refusal's final reimpact must likewise
equal the refusal-sample time. Consequently, the final bracket may straddle a retained root
estimate when the sample itself is exactly that terminal/refusal root. Such an
overhang is admitted only for the final relevant event, only within one declared
macro timestep, and never relaxes strict retained-sample time ordering.

This schema retains accepted public state rather than every integrator stage.
It performs no interpolation, hidden-state reconstruction, calibration, or
physical-authority promotion. Its `f64` times are exact producer values whose
bit patterns are preserved by `render_trajectory_codec`. The focused G0/G3
tests cover valid construction, quaternion double-cover canonicalization,
quaternion/time/frame/unit/contact/base/terminal/QoI refusals, localized
brackets, component identities, and rigid world translation plus `+z` rotation
invariance of intrinsic QoIs. Direct
runner E2E coverage additionally checks full-state/checkpoint equality,
configuration binding, uninterrupted-versus-resumed sample equality, a resolved
1 mm filleted profile, localized terminal admission, prohibited-reimpact
publication, and refusal-checkpoint restart behavior.

`render_trajectory_codec` v3 is the canonical durable transport for that
admitted boundary. It binds a caller-supplied nonzero campaign/operation
content identity, the complete trajectory metadata and raw accepted sample
inputs, and ordered producer-declared timeline seams. It does not serialize
derived visualization points, audio intervals, resampled poses, or artistic
controls: those are deterministically regenerated from the decoded trajectory
under the control-stream and timeline-resampler versions pinned in the wire
header. This single-source rule prevents derived data from silently disagreeing
with the authoritative accepted state.

The wire format uses raw little-endian binary64 bits, a fixed 1,024-sample
chunk policy, length-prefixed sample records, a domain-separated fingerprint
for each chunk, an embedded domain-separated fingerprint for the complete
prefix, and an out-of-band domain-separated identity for the complete artifact
including its fingerprint trailer. The fixed header binds the trajectory,
control-stream, timeline-resampler, interpolation, floating-point, and chunking
versions; sample/transition/chunk counts; exact first and last sample times;
terminal disposition; frame, units, channel availability; and source campaign
identity. Unknown versions, tags, nonzero reserved fields, noncanonical chunk
order or sizing, trailing bytes, truncation, digest mismatch, and a caller's
expected-root mismatch refuse with typed errors.

Callers provide explicit byte, sample, aggregate-transition, and aggregate-text
budgets beneath hard schema ceilings. Seekable decoding first verifies the
complete envelope and every bounded chunk without retaining the sample corpus,
then allocates and decodes admitted sample inputs, re-runs
canonical `RenderTrajectory` admission without renormalizing already-canonical
quaternion bits, validates every declared seam against exact source sample
times, and streams a canonical re-encoding against the original bytes.
Only that exact byte fixed point becomes an `EulerRenderTrajectoryArtifact`.
The decoder leaves no partially accepted artifact on cancellation or refusal;
streaming encode/decode still uses bounded per-chunk buffers and polls the
execution scope at chunk, seam, semantic-admission, allocation, and final
publication boundaries. Convenience in-memory encoding is only for callers
whose explicit artifact budget makes a monolithic buffer suitable.

The transport preserves `SimulationEvidence`; integrity and content identity
are not experimental calibration, cryptographic authentication, physical
validation, cross-ISA floating-point equivalence, or proof that the reduced
contact/aerodynamic/base models match a real Euler disc. A producer restart is
not automatically declared a discontinuity: the composer supplies a seam only
when continuity is not warranted, allowing semantically identical uninterrupted
and checkpoint-resumed trajectories to share canonical bytes when their raw
accepted states and declared composition are identical.

`control_stream` v3 derives synchronized raw rendering and sound controls from
that admitted boundary. `VisualizationControlPoint` is point sampled at an
exact accepted endpoint and carries the exact pose, center-of-mass velocity,
body/world angular velocity, symmetry axis, reduced-base pose/velocity,
post-interval branch, gap, QoIs, and disposition. When and only when the
endpoint branch is closed, it also carries the exact retained contact point and
normal in world, disc-body, and displaced-base coordinates together with the
disc, base, and relative material-point velocities. The reduced base translates
along its local `+z`; v3 admits only base frames whose local `+z` coincides with
world `+z`, while allowing a declared origin and yaw so the controls remain
equivariant under admissible horizontal rigid transforms.

`AudioControlInterval` owns one exact positive-duration accepted interval. For
each explicitly available channel it exposes the duration-weighted mean force
and torque, the exact retained signed work, signed mean work rate `work / dt`,
and force-time/torque-time measures. Positive work means energy transferred into
the simulated disc/body under the producer channel convention. The aggregate
contact channel is not tangential-only work; the reduced-base channel is base
damping work rather than total contact work into the base; and the exterior-gas
channel is body work rather than relative gas dissipation. The separately
retained `interval_normal_force_n` has an explicit
`RenderNormalForceSampling` tag in `RenderChannelAvailability`: `Unavailable`,
`FirstAcceptedSubintervalMidpoint`, `IntervalMean`, or
`AppliedSubstepZeroOrderHold`. The coupled runner uses the midpoint tag;
reduced-decay bridge v2 uses the documented interval-mean tag; and the
transactional production-prefix bridge uses the applied zero-order-hold tag.
Full contact requires a non-unavailable tag. A normal-only midpoint is
diagnostic-only and **cannot** be promoted into a sound force measure; an
applied zero-order hold is a discrete interval mean and remains explicitly
distinguishable from an analytically or quadrature-derived `IntervalMean`.
Duration-mean normal-load authority comes from the full contact-channel mean
projected onto base `+z` whenever that channel is available, or from the scalar
only when its tag is `IntervalMean` or `AppliedSubstepZeroOrderHold` and full
contact is unavailable. Sound mapping for `ContactNormalForce` fails closed
(unavailable) in every other case.
Available numerical zero is never treated as missing.

Each audio interval also declares its exact visualization endpoint coverage.
If the first retained sample owns a positive-duration interval, only its closing
visual state exists; that interval is retained as endpoint-only preroll and is
emitted alone by coarsening. It is excluded from the exposed common positive-
duration audio/visual horizon. Later intervals are fully bracketed by adjacent
retained visualization points. Consumers must not reconstruct preroll motion
from `metadata.initial_state`, because that value describes the original run
configuration and is not necessarily the start state of a resumed segment.
Exact preroll animation requires a future source schema that retains the
segment-start rigid, base, branch, and contact state.

Opening and reimpact records retain their exact source interval, class, time,
and localization bracket. Their excitation measure is explicitly `TimingOnly`:
the trajectory contains neither an event-specific impulse nor a resolved force
history, so raw controls do not invent an amplitude or assign interval work to
one side of a root. Raw controls likewise perform no artistic normalization,
gain, saturation, clamping, spectral shaping, or conversion of mechanical QoIs
into an acoustic frequency.

The only v1 reduction is `WholeIntervalBoxcarV1`. It integrates complete source
intervals before decimation, sums signed work, duration-weights mean wrench and
normal-load controls, recomputes rates over the resulting duration, and treats
every event-bearing source interval as a one-interval barrier. Thus it cannot
blend across a declared contact transition or fabricate subinterval timing.
This is a deterministic preview/control filter, not a fixed-rate band-limited
audio resampler; fractional event placement and 48 kHz synthesis remain owned
by the later audio pipeline.

Both raw derivation and coarsening are transactional, poll their execution scope
at bounded source-interval work, and return no partially published stream.
Finite source values that overflow a derived rate, measure, aggregate, or work
reconciliation refuse with a typed error. Raw and coarsened controls expose
signed-work integral checks. They borrow and pointer-bind the exact admitted
trajectory and inherit its `SimulationEvidence` ceiling; durable canonical
bytes and content identity remain owned by the trajectory codec. Focused G0/G3
coverage exercises resumed clocks, exact frame/velocity derivation, zero and
unavailable channels, zero-force reimpact, opening/event barriers, signed-work
and force-time conservation, alternating-signal boxcar cancellation, admissible
rigid transforms, extreme finite values, derived overflow, deterministic replay,
pre-cancellation, and deterministic injected cancellation inside continuity,
aggregation, post-aggregation, and work-reconciliation loops. These controls
are neither calibrated sound nor physical validation.

`audio_excitation` v2 maps the admitted control stream into a content-identified,
source-clock generalized-force artifact. Canonically ordered controls may select
base-local normal contact force or signed aggregate-contact, rolling-resistance,
reduced-base-damping, and exterior-gas body-work rates. A declared scale converts
each selector's explicit SI unit into a component generalized force. For every
selected interval the mapper retains contact, rolling, base, and gas stems,
their duration means, their force-time measures, availability distinct from
numerical zero, and per-stem arithmetic reconciliation. The work-rate path is a
dimensionally valid transfer-model proxy; it is not a claim that signed work
alone determines acoustic phase, an energy-conserving mechanics force, radiated
energy, calibrated pressure, or absolute SPL.

Normal-force excitation may come either from the full contact-channel mean or
from the reduced-decay bridge v2 normal-load-only authority. In the latter case
the force-time measure is exactly the declared interval mean times duration;
the unavailable aggregate contact wrench/work channel remains unavailable, so
the mapper does not synthesize contact work or torque from that scalar.

Contact and rolling stems are explicitly localized; base and gas stems are
explicitly distributed. Azimuthal mode-shape factors derive from exact disc-body
or displaced-base contact coordinates and apply only to localized stems. Smooth
raw intervals may interpolate exact endpoint factors. A single opening holds the
exact contact-side start factor and a single reimpact holds the exact contact-side
end factor; multiple events or missing required event-side geometry refuse rather
than substituting unity coupling. Because boxcar reduction before a varying
mode-shape projection does not conserve the spatially weighted measure, v2
refuses `WholeIntervalBoxcarV1` together with nonuniform contact-coordinate
participation. Static participation remains compatible with measure-first
boxcar reduction and its event barriers.

Opening and reimpact records remain `TimingOnly` and carry an exact zero physical
impulse. An optional deterministic rolling-noise envelope and randomized reimpact
impulse are separately typed as artistic, domain-separated by source/event/seed,
and bounded before publication. They do not promote the soundtrack's physical
authority. The mapper does not produce `ModalDriveFrame` values or 48 kHz samples:
its interval measures, event brackets, spatial partitions, and explicit
`RequiresBandLimitedResampling` status are the input to the later multirate stage,
which alone owns anti-alias reconstruction, filter latency, fractional-delay
events, and audio/video clock alignment.

Mapper identity transitively binds the durable trajectory, control-stream schema,
modal model, mapping/scaling, reduction, spatial rules, artistic seed/config,
deterministic math and RNG semantics, selected grid, and all budgets. Raw source
work is capped before optional coarsening; each transactional chunk separately
caps intervals, events, and interval-by-mode envelopes, polls cancellation in
preflight and inner bounded loops, and publishes intervals and successor
checkpoint only together. Checkpoints are bound to mapper/source, exact next
interval, preceding end time, and compensated cumulative per-stem measures.
Numerical-refusal trajectories cannot mint excitation; horizon-censored sources
retain their accepted interval history without inventing terminal events.
Focused G0/G3/G4/G5 tests cover every admitted selector and SI transfer unit,
signed stem/measure reconciliation, available-zero versus unavailable sources,
alternating-signal pre-decimation cancellation, event barriers and zero physical
impulse, deterministic artistic replay, exact/event-side spatial factors,
resource refusals, split/resume equivalence, wrong-checkpoint refusal, and
injected cancellation before atomic publication.

`audio_resampling` v1 is the deterministic, offline bridge from those complete
source-clock intervals to the exact 48 kHz modal-drive clock. It admits an exact
24/1 Hz video clock and 48,000/1 Hz audio clock with identical rational start and
exclusive-end instants. Every video boundary must therefore be an integral audio
boundary: the frozen master has exactly 2,000 audio frames per video frame, an
explicit final endpoint marker, and zero integer-clock endpoint drift. Binary64
source-grid endpoints may differ from the exact audio clock by at most
`1e-6` audio frame; this tolerance is only an interchange boundary and does not
replace the source timeline's bit-exact interval continuity.

Source intervals are finite-volume cells, not point samples. The resampler
integrates each retained force-time measure `[N s]` over its overlap with an
audio cell and divides by that cell's exact duration to obtain a held force
`[N]`. For localized contact and rolling drive it first integrates the product
of component participation and the interval's linearly varying, signed modal
location factor; only that already-participated per-mode force is filtered.
Filtering force and location independently is not equivalent and is forbidden
when location varies inside the filter support. Distributed base-damping and
exterior-gas drive remains in component coordinates and uses each mode's static
declared participation later. Output localized component fields are exactly
zero; row-major `(audio_frame, canonical_mode)` arrays carry localized force
`[N]` and left-boundary impulse `[N s]` into modal synthesis's
`PreparticipatedLocalizedDrive` path.

Continuous mechanics-derived controls use a normalized odd-length, centered,
exactly symmetric Blackman-Harris-4 windowed-sinc low-pass whose cutoff is the
arithmetic midpoint of the declared passband and stopband edges. Its odd tap
count is `2h + 1` for declared half-length `h >= 4`. The declared source
bandwidth must fit both the source cadence's conservative Nyquist ceiling and
the passband; ordered
pass/stop edges must also lie below the 24 kHz output Nyquist. Admission measures
the exact versioned coefficients on a declared grid of 8,192 through 32,768
intervals, with at least eight intervals per filter half-length, and refuses a
requested passband contract weaker than 0.1 dB ripple, a stopband contract weaker
than 80 dB attenuation, or coefficients that miss the requested limits. This is
a sampled response audit, not a proof between grid frequencies and not evidence
that the mechanics source contains all acoustic bandwidth. Optional artistic
rolling texture uses a separate content-identified band filter and carries no
physical-filter response authority.

The FIR uses half-sample even reflection at the complete admitted horizon, never
at chunk seams. Centered evaluation explicitly compensates its `half_length`
frame group delay, publishes zero alignment offset, and requires the same number
of future frames as lookahead. It is therefore an offline reconstruction policy,
not a causal or live-processing latency claim; reflection is a deterministic
boundary convention, not modeled physical pre-roll or post-roll. Each chunk
recomputes its global-horizon halo from the immutable complete source payload, so
chunk size and restart boundaries do not reset filter state.

The same module exposes a geometry- and material-independent generalized-force
variant for physical structural solvers. Each anonymous coordinate carries its
own interval force-time measure, is conservatively rasterized by exact temporal
overlap, and is filtered by the identical admitted physical low-pass before
fixed-rate integration. Its identity binds every source measure, clock, filter,
bandwidth declaration, and output value. This prevents a mechanics update clock
from becoming an acoustic tone while avoiding any Euler-object, material-name,
or modal-family preset in the reconstruction layer.

The cinematic physical fixture does not treat a smooth or empty surface spectrum
as a bandwidth certificate for the coupled normal-force and base response. An
interval-measure force stream without producer-supplied band-limit authority is
admitted only against its full conservative source-grid Nyquist. Consequently,
the current 48 kHz force grid refuses the fixture's 18 kHz reconstruction
passband; no physical listening master is claimed until that authority exists.

`AudioResamplingCrop` is the only admitted way to publish a nonzero-offset
subrange while preserving that global boundary condition. Its half-open source
range must begin and end on exact video/audio alignment markers, must fit the
full source horizon, and must have exactly the same duration as its rebased
output clocks. Its identity binds the complete full-horizon resampler identity,
both source offsets, and both output clocks. A cropped `SoundSynthesisConfig`
must name that derived crop identity while retaining the full source
excitation, modal model, filter, and algorithm versions; an independently
restarted short-horizon resampler is rejected. This is an exact provenance and
continuation guarantee for the modeled signal, not a claim that history before
the retained source preroll is physically reconstructed.

Opening and reimpact events retain their source bracket and requested binary64
sample coordinate. A tolerance of `1e-9` frame may snap only near-integral
roundoff. The `LinearTwoBoundaryV1` rule divides an admitted artistic impulse
between neighboring frame boundaries with weights `1-f` and `f`, retaining the
resulting centroid error in its receipt. Nonzero impulses that would require the
exclusive endpoint refuse; zero-impulse timing receipts may name that endpoint.
Physical event impulse remains exactly zero, an opening can never acquire an
artistic reimpact impulse, and artistic event impulses bypass the continuous
control FIR. Event ownership is unique to one half-open output chunk even when
its two signal contributions straddle a chunk boundary.

Resampler identity binds the excitation and modal identities, exact source
payload, deterministic-math and algorithm versions, physical and optional
artistic filter identities, clocks, bandwidth declaration, boundary and event
policies, total horizon, and all resource ceilings. Admission requires that the
excitation mapper was itself constructed against that exact modal identity;
structurally similar modes from another model cannot be substituted. Checkpoints
bind that model and the next absolute audio frame. High-level sound admission additionally
requires exact excitation, modal, resampler, filter, mode-list, version, and
clock agreement. Source intervals, events, synchronization markers, total/chunk
frames, filter taps (hard-capped at 4,097), row-major output and raw-halo mode
values, and the combined raster plus physical and optional-artistic
convolution-work estimate are checked before result storage; allocations owned
by this boundary have typed refusal. Cancellation is polled in admission and at
bounded inner raster, mode, filter, event, identity, and synchronization-marker
work. A chunk, its event and synchronization receipts, and its successor
checkpoint publish atomically, so refusal or cancellation leaves the predecessor
valid and exposes no partial drive. The chunk's `synthesize_modal` handoff
verifies the bound modal identity and exact next-sample checkpoint, then selects
`PreparticipatedLocalizedDrive` itself, preventing a caller from reordering
chunks or silently dropping the already-projected contact/rolling channels.
Absolute frame indexing, immutable source data, versioned deterministic math,
compensated summation, and global boundary handling are the basis of replay and
split/resume equivalence.

Before this boundary supports a release claim, focused G0/G3/G4/G5 evidence must
cover constant-measure conservation; localized force-factor product ordering;
distributed-drive immunity to spatial modulation; pass/stop admission; global
reflection and compensated delay; integral, fractional, opening, reimpact, and
terminal event rules; exact A/V markers; budget and allocation refusals;
pre-cancellation and bounded in-work cancellation; deterministic replay; and
bit-identical one-shot versus split/resume output. Even with that evidence the
artifact is a model-derived modal drive, not calibrated pressure, radiated
acoustic energy, absolute SPL, room response, structural/acoustic coupling, or
physical validation.

`spatial_audio` v1 is the deterministic offline stereo transform between dry
modal/source frames and `audio_artifact`. It accepts an ordered, bounded set of
identified mono point sources. A source may borrow generic mono samples or one
disc/glass/base field directly from `ModalStemFrame`, so the cinematic path does
not copy three full-rate stem arrays. Every source declares a finite
nonnegative linear gain, static or one-per-frame world position `[m]`, and
`PhysicallyParameterized` or `Artistic` parameter provenance. Listener input is
likewise static or one pose per emission frame, with finite position `[m]` and
unit, mutually orthogonal forward/right axes. All sources, per-frame tracks,
and the listener have one exact nonzero common frame horizon; mismatches refuse.

Propagation uses the declared sample rate and positive finite speed of sound.
For distance `d`, the frozen attenuation law is
`d_min / max(d,d_min)` with positive finite `d_min`; it cannot boost above one.
At exact source/listener coincidence, v1 freezes zero delay, centred pan,
front-axis microphone gain, and clamped unit attenuation rather than inventing
a direction or dividing by zero. `IntegerCeiling` places an arrival at
`n + ceil(d fs/c)` and therefore never advances it; `LinearFloorCeil` splits it
between the adjacent sampled arrival frames. The latter is a causal sampled
linear interpolator, not a band-limited fractional-delay filter. Source and
listener poses are evaluated at emission time. V1 does not reconstruct a
retarded listener pose or Doppler shift.

The stereo law is fixed equal-power panning:
`L=sqrt((1-p)/2)`, `R=sqrt((1+p)/2)`, where `p` is the source direction dotted
with the listener's right axis and clamped to `[-1,1]`. Microphone directivity
is either omnidirectional or a first-order cardioid with an explicit rear-axis
amplitude floor in `[0,1]`; it is not an HRTF or head model. Contributions are
accumulated in caller order with deterministic compensated sums. The explicit
source gain is applied exactly once before attenuation, directivity, and pan;
there is no inferred stem balance, normalization, compressor, limiter, hidden
gain, or clipping. Any nonfinite input/derived value or final sample above the
declared absolute ceiling refuses the whole transaction.

`PreserveTail` publishes the complete sampled propagation tail and optional
room-response tail. `ClampToInputFrames` instead publishes exactly the common
input horizon, deterministically discarding deposits and convolution output at
or beyond that exclusive boundary. This choice is bound into configuration
identity. Diagnostics retain both the published count and the complete natural
final count plus discarded-tail count, so a clamped prefix cannot be relabelled
as an untruncated result. Each call is a fresh offline transaction: cuts/reset
are expressed by starting a new call, and v1 never carries a delay line or room
tail implicitly across calls.

An optional stereo room impulse response owns equal-length finite left/right
taps at an exact sample rate and derives its identity from rate, authority, and
every tap bit. Its rate and tap count must agree with the renderer and its
budget. Convolution is deterministic and channel-wise; it is a pure response,
so tap zero must explicitly contain any desired direct path. No gain
normalization or hidden dry signal is added. A room response marked `Artistic`
makes the output `Artistic`; otherwise authority is the conservative
combination of configuration, every source, and the optional response.

Configuration identity binds the algorithm version, exact sample rate, speed
of sound, minimum distance, delay and output-horizon policies, microphone law,
parameter authority, and all resource/amplitude ceilings. Input identity binds
that configuration, ordered upstream identities, selected modal fields, source
gains, every sample bit, every resolved source/listener pose, and optional room
identity. Output identity additionally binds exact stereo sample bits, final
authority, and room identity. The distinct dry-bypass transaction validates
rate, identity, budgets, finiteness, and amplitude, then copies already-stereo
frames bit-for-bit (including signed zero) without invoking any spatial law.

Sources, total input frames, output frames, room taps, owned sample bytes, and a
checked deterministic work estimate have explicit caller ceilings plus hard
source/rate/tap limits. Owned buffers use fallible reservation. Admission,
source traversal, accumulation finalization, convolution, validation, and
identity traversal poll the supplied `Cx` at no more than 256 frames/taps
between checkpoints. Output, identities, and diagnostics publish only after
all work and a final checkpoint succeed; refusal or cancellation exposes no
partial result.

Focused G0/G3/G4/G5 evidence covers static and moving left/right geometry,
centred balance, near/far attenuation, listener motion, co-location clamping,
speed-of-sound delay scaling, rear microphone response, integer/fractional
impulses, room tails and rate/identity checks, preserve-versus-clamp horizon
identity and diagnostics, unity/non-unity source gains, modal stem selection,
exact dry bypass, nonfinite/clipping refusal, pre- and in-work cancellation,
and bit-stable replay. Even with that evidence, output is a deterministic
presentation transform, not BEM, HRTF, occlusion/diffraction, measured room
acoustics, calibrated pressure, radiated power, absolute SPL, or perceptual
validation. `PhysicallyParameterized` records input provenance only and does
not promote those no-claims.

`audio_artifact` v1 turns complete modal or separately spatialized samples into
the deterministic stereo artifact consumed by the cinematic finalizer. Its
production format surface is deliberately narrow: RIFF/WAVE at exactly 48 kHz
and two channels, using either packed signed little-endian PCM24 or IEEE-754
binary32. Float32 WAV is the authoritative `euler-disc-v1` audio master; PCM24
is a deterministic quantized derivative and cannot silently replace it.
The writer emits one canonical order: `RIFF/WAVE`, `fmt `, `fact` for float32,
an optional bounded `LIST/INFO/ICMT`, then `data`. V1 metadata admits only
printable ASCII plus line feed; it does not guess a RIFF text encoding or emit
an undeclared code page. Chunk and RIFF sizes use
checked arithmetic, odd payloads receive one zero pad byte excluded from their
declared size, and RF64, compression, extensible layouts, arbitrary channel
orders, and unknown chunks are explicit no-claims. The strict reader accepts
only this emitted subset, exact internal rates/widths/counts and exact EOF;
well-formed features outside it return `Unsupported` rather than being guessed.

Dry input is the canonical disc, glass-plate, base-assembly stem order. Each
mono stem has an explicit decibel gain and pan `p` in `[-1,1]`; equal-power
coefficients are `sqrt((1-p)/2)` left and `sqrt((1+p)/2)` right, evaluated by
the versioned deterministic math core. Three terms are accumulated in fixed
order with compensation, followed by one explicit master gain. A distinct
channel-layout receipt identifies already-spatialized stereo and requires a
nonzero spatialization identity; that path is copied without repanning. Dry and
spatialized inputs cannot be blended implicitly, and the dry mix has no camera,
listener, room-response, or renderer dependency.

There is no limiter, compressor, normalizer, dither, DC blocker, or hidden gain.
All input must be finite. The configured headroom ceiling is
`10^(-headroom_db/20)` and both stored-sample peak and the declared intersample
estimate must fit it or the whole transaction refuses. PCM24 uses round-to-even
at scale `2^23` followed by signed-24-bit clamping, so `-1` is exact and `+1`
decodes as `1-2^-23`; float32 preserves finite converted bits including signed
zero. Final meters are recomputed over decoded format semantics, not the
pre-quantized `f64` mix.

Meters report exact stored-sample peak; population stereo RMS and per-channel
DC estimates; a four-times Lanczos-8 windowed-sinc intersample estimate under
half-sample-even boundary extension;
and a BS.1770-derived digital programme-loudness diagnostic. Loudness uses the
frozen 48 kHz two-biquad K weighting, complete 400 ms blocks at 100 ms hops,
stereo weights of one, a strict `-70 LUFS` absolute gate, and a strict relative
gate ten loudness units below the absolute-gated energy mean. A programme under
400 ms, silence, or an empty final gate returns unavailable loudness rather
than NaN or infinity. The Lanczos peak is not the BS.1770 Annex-2 filter, a
continuous-time supremum, or a standards certificate; the loudness result is
not an EBU normalization target, perceived-loudness proof, calibrated SPL, or
radiated acoustic power.

The typed manifest binds the complete `SoundSynthesisReceipt` (and therefore
trajectory, excitation, model, timeline, authority, assumptions, resampler and
filter configuration), an exactly matching caller-declared source-synthesis
receipt, dry/spatialized channel receipt, source-sample and
optional mix identities, exact 24/48 kHz start/end ticks, 2,000 audio frames per
video frame, stored sample count/encoding, metadata identity, complete WAV byte
count/hash, configured headroom, and final decoded-sample meters. Binary typed
fields define manifest identity; JSON is only a deterministic view. WAV bytes
do not embed their own final hash or manifest identity, avoiding a hash cycle.
The matching receipt catches accidental cross-configuration relabeling, but
raw sample slices are not an independently replayed proof that upstream modal
synthesis or spatialization actually produced them. Configuration authority is
declared unchanged and is never minted by mixing, metering, encoding, or a byte
hash; this boundary alone therefore cannot establish calibrated acoustics.

Frames, bytes, metadata and deterministic work have caller-controlled ceilings
under the RIFF `u32` hard limit. Owned buffers use fallible reservation.
Cancellation is polled at fixed sample boundaries and 64-KiB byte-hash
boundaries in source/WAV hashing, mixing, peak interpolation, K weighting,
block accumulation and gating, encode and decode. Aggregate work is preflighted
before high-level frame traversal and allocation. Public output
is an in-memory transaction published only after the complete WAV, strict
round-trip, meters and manifest agree; cancellation exposes no partial artifact.
Focused G0/G3/G4/G5 evidence must cover known RIFF bytes, zero/one/budget edges,
PCM extrema and near-positive-full-scale rounding, exact float bits/nonfinite
refusal, metadata padding, every proper prefix, deterministic mutation/junk,
dry pan/gain and spatialized separation, no-normalization/headroom refusal,
meter and loudness boundaries, exact A/V count, byte replay, verification and
cancellation. External decoder inspection is compatibility evidence only, not
authority or physical validation.

`modal_synthesis` v2 is the deterministic 48 kHz damped-resonator runtime for
the reconstructed audio controls. Every canonical mode solves
`m q'' + 2 zeta m omega q' + m omega^2 q = F`: a declared boundary impulse is
applied first as an exact velocity jump, then a constant generalized force is
integrated over the frame by an exact zero-order-held transition. Stable
underdamped, critical, overdamped, and small-step coefficient branches use the
versioned deterministic `fs-math` core. Modes are sorted by unique nonzero ID;
distinct degenerate modes remain distinct.

Each drive frame separates localized and distributed component force `[N]` and
left-boundary impulse `[N s]`. `Declared` participation applies each mode's
static signed component vector to both classes. `PerFrameModeFactors` additionally
multiplies only localized drive by one finite normalized factor in `[-1, 1]` per
`(frame, canonical_mode)`; distributed drive is invariant to that factor. The
resampling-safe `PreparticipatedLocalizedDrive` instead accepts localized
per-mode force and impulse arrays directly, while still applying static modal
participation to distributed fields. Its localized component fields must be
exactly zero or synthesis refuses, preventing accidental double routing. This
per-frame factor array, when selected, and each direct array must have exactly
`frames * canonical_modes` finite row-major values. The direct path is required
when anti-alias reconstruction spans a changing contact location; the modal
runtime does not reconstruct or filter either operand.
Each mode's radiation is assigned to one dry component stem. Off-component
participation is explicit source routing, not a coupled structural solve.

The model identity binds the algorithm and deterministic-math versions, sample
rate, complete canonical modes, routing class, output convention, and every
resource/state limit. A checkpoint binds that identity, the next sample frame,
and every canonical `(q, q')` state. Chunk synthesis preflights the complete
drive and spatial-factor envelope, reserves all result storage, polls the
execution scope at bounded frame intervals, and publishes samples and successor
state only together. Displacement, velocity, per-mode energy, summed energy,
and dry-output ceilings refuse rather than clamp; boundary impulses are checked
before decay can hide an excessive kick. Outputs include end-boundary mono
samples, component stems, per-frame internal modal energy, final per-mode
kinetic/elastic energies, and explicit population RMS and absolute peak.

Focused G0/G3/G4/G5 coverage compares impulse, step, and harmonic response to
independent analytic or exact-discrete oracles; exercises zero, critical, and
high damping, Nyquist guard refusal, unforced energy decay, mode permutation,
duplicates, degeneracy, zero and cross-component participation, signed spatial
factors, direct preparticipated localized routing, distributed-drive immunity to
spatial factors, bounded state, cancellation before and during work, exact
replay, and bit-identical split/resume. Built-in
tungsten/stainless-disc plus glass/base presets are
`RepresentativeUncalibrated`: their frequencies, damping, masses, participation,
and radiation gains are plausible film parameters, not measured eigenmodes,
radiated acoustic energy, absolute SPL, structural/acoustic calibration, or
physical validation.

`structural_acoustics` is the physical replacement path for those
representative presets. A resolved elastic specimen is tetrahedralized by the
shared rounded-cylinder mesh primitive; `fs-solid` assembles physical 3-D
`(K,M)` operators from its global-frame Mandel stiffness tensor; `fs-modal`
returns certified mass-normalized modes; the
actual body-frame contact point and force are projected through their P1 mode
shapes; and `fs-bem` maps the same boundary-normal mode shapes into SI pressure
at an explicit gas-state-dependent observer. Low acoustic size uses plain CBIE
to avoid the documented low-`ka` Burton--Miller resistance artifact; the
higher-frequency arm uses Burton--Miller for fictitious-frequency protection.
Any negative outgoing power refuses.

`PhysicalModalAudioModel` evaluates material-model loss at those structural
frequencies, projects each interval mean to a generalized force-time measure,
conservatively reconstructs those measures on the audio clock, advances every
mass-normalized oscillator by an exact constant-force transition over each
reconstructed audio cell, and emits unmastered pressure in pascals. The current
spatial approximation is declared explicitly: an interval mean force uses the
closing retained contact point/body orientation, or the opening endpoint when
the closing endpoint is open. Observer-independent BEM far fields are projected
into caller-bounded spherical-harmonic tables. A world-fixed microphone then
evaluates spherical spreading, modal propagation phase, and body-frame
directivity from the actual resampled rigid pose at every audio boundary;
simultaneous microphones share one oscillator state so stereo phase cannot
drift. No material name selects a frequency, decay constant, radiation gain,
pan curve, or digital level.

The broadband structural source stem reuses the certified enrichment modes of
`StructuralResidualFlexibilityEstimateBasis` as its fixed scalar generalized-force basis; it never
adds the static-residual response. The producer binds exact basis, operator, material, damping, gas,
frequency-partition, and solver identities. On disjoint training and withheld grids it uses the
`exp(-i omega t)` relation `v_n/a = +i phi_n/omega`, performs one multi-right-hand-side BEM solve per
frequency, fits body-frame real-tesseral spherical-harmonic transfers only from training data, and
retains direct BEM far fields for withheld validation.

The rigid-disc companion uses the same closed specimen mesh and gas/BEM/SH
pipeline for three rigid translations. Accepted-step momentum differences are
converted to linear acceleration at mechanics cadence, projected into the
material frame, and passed through the staged anti-alias decimator before
reaching 48 kHz. Rotational coordinates are excluded because their
low-frequency boundary-work estimates are not passivity-admissible on the
production mesh. This is a linear low-Mach far-field estimate: it does not
claim moving-boundary/FW-H, rotational or convective radiation, near field,
room response, two-way fluid loading, calibrated SPL, or accuracy beyond the
declared BEM/directivity/fit evidence. Its fitted rigid-body bank trains through
10 kHz and is held out through 10.5 kHz because the production surface mesh does
not meet the six-panels-per-wavelength admission floor at 12 kHz; the independently resolved support-plate
path retains the complete 18 kHz reconstruction passband.

At runtime, existing P1 point-force projection and conservative audio-cell reconstruction drive the
exact modal transition with zero pressure transfer. The closing state supplies
`a = Q - 2 zeta omega qdot - omega^2 q`; a persistent passive bank emits frame-major,
channel-ordered `Pa*m` coefficients. These are source coefficients, not listener pressure or a
`PhysicalPressureSignal`. This `EstimateOnly` path makes no claim of calibrated absolute SPL,
listener propagation, room response, nonlinear structure, or two-way acoustic loading.

The generic `InteriorOnly` retarded observer converts that body-frame source
stem to simultaneous world-fixed pressure signals. Frame zero is explicitly at
the closing boundary `start + 1/fs`. At each output boundary, deterministic
bisection solves `tau + |x_observer-X(tau)|/c = t`; all real-tesseral source
coefficients use one complete 16-tap Lanczos-8 stencil at `tau`, and their
direction is rotated by the emission pose before physical `1/r` spreading.
All observers share the intersection of their complete arrival-time horizons
and are returned transactionally in caller order. Admission enforces exact
source/basis/gas/clock identities, the conservative broadband far-field gate,
and a caller-bounded surface Mach number. This remains `EstimateOnly` under
`RETARDED_FAR_FIELD_OBSERVER_NO_CLAIM`; in particular it is not moving-boundary
FW-H, exact Doppler amplitude, near-field, room/head, or calibrated SPL.

The modal initial condition is explicit. `Zero` is admissible only when the
retained horizon begins before excitation; a cropped horizon that begins under
a held contact force uses `StaticEquilibriumAtFirstHeldForce` so that truncating
unavailable prehistory does not invent a start-up impulse. Simultaneously
radiating bodies are combined by `superpose_pressure_signals` in SI pascals
before any presentation mastering. The superposition validates a common
observer, clock, sample count, and contact sampling convention, canonicalizes
component identities and summation order, and emits aggregate structural,
radiation, and damping identities. Exactly one explicit pressure-to-digital
gain may then be applied to the composite field. The cinematic listening
derivative uses one caller-declared fixed `FS / Pa` calibration across the
whole artifact; it is never derived from that clip's peak or loudness. If the
fixed calibration exceeds the declared true-peak ceiling, publication refuses
instead of normalizing, limiting, or compressing the pressure signal.

The cinematic product path additionally simulates a caller-selected, integral
number of 24 Hz frame intervals before the synchronized picture/sound cut.
This is coupled physical history: rigid disc, contact laws, gas, rolling loss,
and the structural base all advance through the ordinary production step, and
the exact same prefix warms acoustic reconstruction before it is cropped. The
default is 48 frames (2 s), with an admitted range of 1 through 240 frames.
Its purpose is to move the media cut past the artificial free vibration caused
by initializing a moving-load problem from instantaneous static equilibrium;
it is not an audio-only fade, a fabricated sound source, or timestep-convergence
evidence. The manifest records the selected preroll horizon.

The production cinematic default mechanics clock is 1.536 MHz. Stability is
horizon- and state-dependent: in the retained `0.02 rad` qualifier, 384 kHz
produced contact-branch chatter while 768 kHz and 1.536 MHz stayed on the
smooth-contact branch, but a later `0.01 rad` six-second qualifier produced
4,528 branch transitions at 768 kHz and zero at 1.536 MHz, with more than a
200-fold difference in peak radiated pressure. The 768 kHz member is therefore
rejected for terminal audio. A caller may explicitly request the admitted
3.072 MHz and 6.144 MHz diagnostic rungs; exact 12.288 MHz is additionally
admitted only as their bounded next-refinement diagnostic. No current default,
pair, or single artifact establishes asymptotic timestep convergence.

The fixed supported-plate path does not reuse a natural-frequency Helmholtz
transfer for arbitrary forced motion. It advances the physical modal state
under the contact reaction using the same C1 rectangular reconstruction of
retained DKT nodal displacement and slope degrees of freedom as the live modal
support, evaluates each mode's instantaneous normal
acceleration, and applies a causal Rayleigh-I surface integral at every
triangle's finite sound-travel delay. Linear fractional-sample delay taps are
fixed by panel geometry, observer position, gas sound speed, and the declared
sample rate. Thus a sub-resonant moving contact radiates at its simulated
forcing harmonics instead of being relabeled as a plate eigenfrequency. This
remains linear small-displacement radiation into an infinite rigid baffle; it
does not include the housing cavity, edge diffraction, room response, or
two-way radiation impedance.

That C1 rectangular reconstruction is a smooth moving-load evaluation field,
not a claim that bicubic Hermite polynomials are the exact triangular DKT
element shape functions. Distributed projection of the admitted finite contact
pressure through the native element basis remains a higher-fidelity rung; the
current structural port applies the accepted point resultant at the resolved
contact location.

Plate support is geometry, not an apparatus-name preset. A request may constrain
the full perimeter or declare three centered-frame pin locations. Three-point
pins are deterministically resolved to distinct, non-collinear structured-mesh
nodes under an explicit snap tolerance; the requested locations, resolved node
indices, and maximum snap error are identity-bound and reported. A pin fixes
transverse displacement while leaving DKT rotations free. This is an ideal
kinematic boundary condition, not an elastic foot, housing, or table-impedance
model; those require explicit mass, stiffness, damping, and contact state.

No-claim boundary: both the legacy body-fixed transfer and the pose-dependent
directivity realization retain only each structural mode's natural-frequency
Helmholtz response. The latter applies that narrow-band propagation phase at
the current pose; it is not a broadband moving-boundary retarded-time solve,
Doppler model, or room response. Those require the passive rational frequency-
response and propagation-history rungs tracked separately. The current small-
strain solid basis also cannot cross yield, large deformation, or phase change;
a hot lead request must refuse when its evidence-bearing solid state leaves its
validity domain. The current enthalpy/phase ingress names the required
mass-conserving free-surface escalation, but no liquid deformation, remeshing,
flow, phase-dependent acoustics, or phase-dependent optics is yet executed.

`timeline_resampling` v1 reconstructs render/audio query times from this
admitted state without mutating the source artifact. Center-of-mass translation
and base displacement use cubic Hermite interpolation with the accepted
endpoint velocities. World linear momentum, body angular momentum, and base
velocity remain finite continuous reconstructions. Orientation uses normalized
shortest-arc quaternion SLERP, including explicit small-angle and `q/-q`
handling. Exact source times reproduce the accepted continuous state rather
than passing it through the interpolator. The method/version is exposed as
`CubicHermiteSlerpV1`; the later codec owns composition identity. Every declared
continuation/producer discontinuity must coincide with an accepted source
sample; otherwise admission refuses because neither one-sided reconstruction
has a trustworthy state boundary.

Contact branch, terminal/refusal disposition, and event class are discrete.
They are never averaged. Queries exactly at a localized event require an
explicit left- or right-limit policy, and retain the original event time and
localization bracket. Shutter intervals can either subdivide at every retained
contact, terminal, continuation-seam, or producer-declared boundary, or refuse
cross-event blur. Query sequences must be finite, strictly increasing, within
the accepted sample horizon, and below the trajectory resource ceiling; there
is no silent extrapolation. Non-representable reconstruction refuses instead
of emitting NaN/inf state.

This is a visualization and control-signal reconstruction assumption, not a
solver refinement, contact impulse reconstruction, or physical validation
claim. Declared continuation seams supply partition metadata only; they do not
assert continuity. G0/G3 coverage includes exact endpoints, analytic constant
translation/base motion, tiny and near-pi shortest-arc rotations, quaternion
double-cover invariance, event-side branch/terminal semantics, shutter
subdivision/refusal, declared seams, strict query refusal, frame-rate nesting,
time/rigid-translation equivariance, determinism, and extreme finite-time
refusal.

`render_motion_bridge` v1 is the L6 visualization adapter from that resampler
to L5 render-time contracts. It resolves each frame shutter inside the exact
accepted trajectory horizon, applies explicit `Subdivide` or `Refuse` policy
before ray evaluation, and treats a genuinely zero-width shutter as one static
time. A positive requested duration that collapses at the absolute-time
binary64 resolution refuses. Prepared shutters are bound to the exact source
trajectory, declared-discontinuity model, and event policy and retain source
authority unchanged. A `TimedRay` must retain the identical admitted shutter;
matching one coincident absolute time from a different exposure is insufficient.
Multi-segment `Subdivide` shutters cannot pass through the global-coordinate or
global-ray sampling APIs. The caller must select an explicit segment-local
shutter; segment endpoints deterministically use the left limit at an interior
closing event and the right limit at the following segment opening, and each
segment exposes its positive duration weight.
Each query is reconstructed by `TimelineResampler`, then maps the mechanics
quaternion from `(w,x,y,z)` into the renderer's `(x,y,z,w)` convention and maps
center-of-mass world metres directly into the proper-rigid translation. Tests
pin endpoints, quaternion double-cover equivalence, cross-shutter and
cross-event-model refusal, enforced segment selection and one-sided event evaluation,
collapsed-exposure refusal, horizon refusal, zero-width equivalence, replay,
and the unchanged `SimulationEvidence` authority ceiling.

The adapter does not promote interpolated visualization poses into accepted
solver states, assert continuity across a returned event partition, or add
mechanical bandwidth. The prepared global-coordinate and global-ray APIs refuse
unsplit sampling of a multi-segment exposure. The raw absolute-time resampling
API remains available for expert non-exposure queries and does not apply
shutter admission or event policy. Downstream weighted image accumulation over
explicitly selected segments remains a separate composition responsibility.

## Animated Euler scene bridge v1

The opt-in `cinematic-render` feature composes one decoded
`EulerRenderTrajectoryArtifact`, its exact `ResolvedDiscProfile`, a physical
camera timeline, and declared visual configuration into a renderable scene.
Admission recomputes strong domain-separated identities for the complete
ordered line/arc meridian, density-bound profile, and resolved mass properties;
all three must exactly match trajectory metadata. The bridge also checks SI
metres/radians, the center-of-mass mass reference, object-ID uniqueness, and
camera coverage over the accepted trajectory horizon. Every retained signed
gap is recomputed from the bound exact chart at the retained pose and displaced
base plane; closed samples must additionally reproduce the exact support point
and profile feature and align their normal with the base. Camera depth admission
uses conservative subject bounds over the full animated horizon, including
cubic-Hermite Bezier translation controls and arbitrary rotation, rather than
endpoint positions alone. The bridge also applies an explicit temporal
angular-sampling ceiling. The angular guard takes the larger
of endpoint angular-speed times sample duration and the retained quaternion's
shortest-arc separation. It is an alias diagnostic, not a proof that the source
contains every intervening revolution or mechanical frequency.

`AxisymmetricChart` intentionally carries no certified tracing claim, and the
production tracer correctly refuses to interpret it as certified ray geometry.
The bridge therefore performs one bounded, deterministic visualization-only
surface-of-revolution conversion from the chart's exact retained line and
circular-arc segments. Circular arcs use a fixed declared number of equal-angle
chords; azimuth uses a fixed declared ring resolution. The output mesh is
centered at the resolved center of mass before the trajectory transform is
applied, so construction origin and COM cannot be silently conflated. Its
receipt binds the exact source chart, resolution, canonical vertices and
triangles, topology counts, local bounds, chord-sagitta diagnostics, and a BVH
layout fingerprint. These are approximation and replay diagnostics, not a
watertightness, Hausdorff, shading-normal, or direct-chart intersection
certificate. Mechanics, contact, support, mass, and inertia continue to use the
original resolved chart rather than the render mesh.

The emitted beauty scene uses stable semantic indices and object IDs for the
animated disc and base plate, static housing, an optional finite studio support
surface, and one rectangular emitter. The optional support is a closed
base-local box whose top is declared relative to the housing bottom; its exact
dimensions, gap, material, geometry identity, and object identity are bound to
the configuration and scene identities. It is visual studio geometry rather
than mechanics support, and is deliberately excluded from the animated-subject
camera-coverage bound. A
configured contact marker has a separate identity/receipt and is absent from
all beauty scene and frame APIs; only the explicitly named static diagnostic
scene/render APIs append it after the beauty primitives. The
plate follows the reduced base's local `+z` displacement and velocity composed
with its nominal rigid frame; the housing remains attached to the nominal base
frame. Plate and housing dimensions are declared visual inputs and remain
representative unless separately bound to measured apparatus geometry. The
reference plate is a closed, outward-wound box rendered as homogeneous spectral
dielectric glass. Its Cauchy dispersion, Beer-Lambert absorption, polished-GGX
boundary, and explicit `RepresentativeCrownV1` provenance are bound into the
scene identity. The preset is a look-development starting point, not measured
stock or experimental calibration data. Tracer paths preserve raw boundary
sidedness, apply absorption over physical in-medium segment length, and enforce
strict path-local LIFO entry/exit semantics. Encountered reversed winding,
non-nested overlap, or escape from a declared closed medium refuses rather than
silently choosing an index of refraction. This local runtime check is not a
global watertightness, orientation, or non-overlap certificate for the mesh.
`EulerMaterialStyle::Conductor` carries the L5 renderer's admitted complex-IOR
optics and explicit isotropic roughness. The cinematic cross-domain binding
requires bulk mechanics and conductor/dielectric optics to resolve from the
same immutable material card and thermodynamic state; surface finish retains a
separate identity because it is not implied by bulk chemistry. Scene look
presets remain uncalibrated look-development inputs, not predictions of a
product, alloy, oxide/passive film, machining, or finish.
Layered engraving, anisotropic brushing, environment lighting, measured
specimen fitting, and final studio look development remain downstream
capabilities. The renderer has exact next-event connections through admitted
parallel slabs and isolated convex planar two-interface smooth dielectrics,
including the reference plate's finite bevel/side faces. Those are generic
unbiased manifold proposals with reverse MIS, not straight undeviating shadow
rays. The tracer still makes no general bidirectional, vertex-merging, curved
specular-manifold, nested-media, or arbitrary caustic claim; paths outside the
admitted connector classes may converge slowly.

Frame requests resolve their shutter and contact/event policy before tracing.
A zero-width request renders one exact time. A positive exposure crossing a
retained event either refuses or yields ordered event-delimited shutters with
explicit duration weights. Rigid-geometry cinematic frames may instead admit
contact opening/reimpact boundaries because those change a discrete solver
branch while the rendered pose and base state remain continuous; the events
remain retained in the exposure provenance, and terminal or producer-declared
discontinuities still refuse. Event-delimited films remain separate: v1 will
not fabricate sample-count or time-mode provenance by silently combining them.
Producer-declared continuation/discontinuity seams currently refuse scene
construction because this instance representation cannot encode trustworthy
one-sided transform keyframes. Camera gaps and poses wholly outside the trajectory
horizon are irrelevant; source times and in-horizon camera keyframes must meet
the declared near/far depth admission. Those depths are validation requirements
over scene bounds, not tracer clipping planes.

The serial segment/frame methods remain the reference path.
`render_segment_with_execution` and `render_frame_with_execution` expose the
same admitted cinematic shutter through an explicit tile shape, worker count,
operation-memory ceiling, and logical run identity. They return the renderer's
tile, executor, and memory-admission report and preserve `RenderExecutionError`
as a structured `EulerSceneError` variant rather than flattening it into a
generic scene refusal. For animation batches,
`render_segment_with_parked_scope` and `render_frame_with_parked_scope` accept
the callback-scoped crew from `RenderWorkerPool::with_parked_crew_local`, so
successive jobs reuse one structurally joined worker crew. Each parked job
retains its own tile shape, memory ceiling, and logical run identity; its worker
count, scheduling weights, and execution mode must match the parked crew.
`begin_segment_render` and `begin_frame_render` bind the exact Euler scene,
camera exposure, settings, layout, execution mode, compute budget, and run
identity into an opaque single-film `PendingRender`; cancellation returns a
resumable suspension with row-prefix progress but never exposes partially
accumulated pixels.
The corresponding `render_*_adaptive_with_execution`,
`render_*_adaptive_with_parked_scope`, and `begin_*_adaptive_render` methods
bind the same cinematic scene and shutter to fs-render's versioned raw-sample
adaptive policy. Their returned `AdaptiveFilm` retains raw sums, Welford
moments, sample counts, and terminal decisions; their opaque pending jobs retain
only complete tile-row prefixes across in-process cancellation. The dispersion
proxy remains a scene/profile-specific render heuristic, not a statement about
Euler mechanics, physical fidelity, perceptual error, or a universal 4K preset.
Event-delimited multi-film composition remains explicit in every execution
mode.

The separate opt-in `render-checkpoint-ledger` feature connects those opaque
pending jobs to `fs-ledger` without changing their rendering semantics. A v1
binding names the exact source trajectory artifact, the complete admitted
scene-builder configuration, the resulting beauty scene, and a canonical
event-delimited frame identity. That frame identity covers cut ownership,
segment index, resolved shutter endpoints/convention/distribution, and duration
weight. The render-job identity additionally covers fixed-spp versus adaptive
kind, every `Settings` field, every `RenderExecutionConfig` field including
raw scheduling weights, the admitted execution mode and exact `Budget`, the
runtime ISA and sorted detected feature set, the complete cinematic shutter
mode and shot identity, and every adaptive-policy field plus its estimator
semantics version. `begin_uniform_checkpoint_job` and
`begin_adaptive_checkpoint_job` atomically admit the raw renderer job and
freeze the canonical Euler frame identity from the same prepared segment.
Persistence accepts only the resulting sealed, non-cloneable Euler job; it
never accepts a raw pending render or a caller-created binding. The wrapper
privately carries both the evolving pending state and its latest durable head
through bounded advances, suspensions, parked-crew attempts, and restores.
Store derives generation zero when no head exists and otherwise derives the
next generation and predecessor digest from that private head. Safe code
therefore cannot replace the state while retaining successor authority.
Producer build and claim identities are explicit and nonzero and may change
between generations, while source, configuration, scene, frame, and render-job
identities remain fixed. A separate expectation type can reconstruct a root or
typed successor binding only inside restore; it cannot publish fresh state.
Successor restore reopens the typed predecessor artifact as well as the proposed
successor, strictly revalidates the predecessor receipt, and asks fs-render to
prove that every predecessor-committed accumulator/AOV bit is retained, every
tile-row prefix and attempt count is nondecreasing, and all newly uncommitted
state remains canonical zero. Merely naming a valid predecessor digest in a
canonical lower-level checkpoint therefore cannot mint state-ancestry
authority. This is a safe-adapter construction guarantee, not cryptographic
producer authentication, scheduler ownership, or concurrent single-writer
claim arbitration.

Checkpoint codec v1 bytes stream directly through
`Ledger::artifact_writer/write/finish`. Dropping or failing that writer rolls
back its transaction, so a prior immutable content-addressed checkpoint
remains readable. Before `finish`, the adapter counts every streamed byte and
requires exact agreement with fs-render's sealed receipt. After successful
commit, constructing the typed stored receipt and replacing the wrapper's head
are infallible; the adapter cannot report an ordinary validation error after
publication. Restore checks ledger metadata for the exact
`euler-render-checkpoint-v1` artifact kind, then uses `get_artifact_bounded`
with the caller's explicit byte ceiling and invokes fs-render's strict v1
decoder against the expected binding. For a root the ceiling covers that one
artifact; for a successor it covers the aggregate predecessor-plus-successor
bytes that must coexist during extension verification. A successful restore returns a sealed
Euler job whose private head is the typed, verified stored checkpoint, so a
post-restart advance and store derives its successor without reconstructing or
injecting lineage scalars. Encode and decode receive the caller's `Cx` and poll at
bounded codec chunks/tiles; cancellation refuses before publication or restored
state escapes. The raw BLAKE3 ledger artifact key covers the complete sealed
bytes, while fs-render's separately domain-separated receipt digest covers the
canonical body named by the seal; both identities are retained and are not
falsely equated. Their byte counts must agree. Reopening a ledger does not
weaken those checks. This adapter does not invent a universal checkpoint-size
constant or assert that any fixed memory limit admits 4K uniform or adaptive
state; callers must budget each job and fs-render remains the concrete
memory/codec admission authority.

The opt-in `render-sharding-ledger` feature adds a bounded L6 coordinator for
uniform fixed-SPP cinematic work. `EulerUniformRenderPlan` canonically sorts
stable frame ordinals, retains every event-delimited prepared segment as a
separate film, and partitions each segment's complete tile/sample cell space
exactly once. Its externally pinnable identity covers the sequence, source
trajectory, scene configuration and built scene, complete render settings,
tile/sample partition, finishing-neighbor radius, every exact checkpoint-derived
frame-segment identity, and all resource caps. Plan admission checks frame,
shard, canonical-plan-byte, per-shard path, per-shard result-byte, and per-segment
aggregate-result limits before returning the plan. Its strict canonical codec
recomputes ordering, coverage, identities, derived counts, and byte length under
an external plan pin. Construction derives the exact segment, shard, and plan
byte counts, then proves every per-shard path/result and per-segment aggregate
cap before allocating the canonical frame-order map, so a certainly over-cap
plan cannot force that retained allocation; decoded bytes alone do not authorize
a render.

Each worker must re-present the original scene-bound prepared frame, from which
the coordinator reconstructs the exact shutter, cut side, frame identity, and
generic renderer shard spec. Workers may exchange only the renderer's bounded
canonical result bytes. A single coordinator strictly decodes those bytes and
stores them under the dedicated immutable ledger kind; workers do not share or
write the coordinator ledger. Segment merge first validates the complete
reference set without artifact I/O: exact duplicate references are idempotent,
an unexpected logical shard refuses, and differing artifact references for the
same logical shard return the lowest conflicting logical identity independent
of arrival order. After that structural pass, it preflights every unique
selected artifact kind and length together against the aggregate cap before
reading any payload. Payload materialization uses the ledger's controlled
64-KiB tiles and polls the caller's `Cx` between tiles. Missing, foreign,
corrupt, aliased, or conflicting work publishes no film. Raw segments remain
independently renderable; `finishing_neighbors` exposes only explicit
frame-position dependency metadata for a later temporal finishing pass and does
not blend films. This feature is a portable file/byte boundary and local
multi-worker plan, not a network scheduler, lease/claim protocol, distributed
database, adaptive-sampling sharder, or cluster-fault-tolerance claim.

Focused G0/G3/E2E coverage builds the scene from a real 1 mm circular-filleted
disc, checks deterministic scene and mesh identities, COM/base transforms,
fillet chord retention, stable primitive binding, authority preservation,
restyling identity separation, and irrelevant camera-history isolation. It
also exercises asset/unit/angular/discontinuity/cancellation refusals and
explicit pre/post-event segment renders. A tiny real spectral-path-traced frame
must contain finite illuminated spatial variation, reproduce byte-exactly on
the same cinematic path, agree bit-for-bit between the serial and explicit
tile-parallel cinematic paths, and render two different frame jobs bit-for-bit
against their serial oracles while reusing one parked crew and preserving each
declared run identity. The adaptive cinematic surface must likewise produce
bit-identical raw sums, moments, counts, and decisions through one-shot,
parked-crew, and opaque-pending frame entry points while preserving the
declared adaptive run identity. Separately, the uniform cinematic frame must
agree tightly with the materialized static seam and round-trip through the
in-house floating-point EXR codec; adaptive-film EXR serialization is not
claimed here. Changing only
the plate phase-index magnitude from representative crown glass to a
same-absorption, same-dispersion-class lower-index control must produce a
finite, quantitative multi-pixel radiance difference. Keeping both materials
dispersive prevents packet-collapse policy from masquerading as optical
activity in this comparison.
This evidence proves the software composition path, not the physical fidelity
or 4K quality of a finished Euler-disc film.

The ledger-checkpoint G3/E2E additionally advances every uniform and adaptive
tile to a nonzero but strictly partial row-atomic safe point, persists that
state, and requires an exact final-buffer match after closing and reopening the
ledger and resuming. It also requires refusal after any binding or job identity
mutation, the same-shutter/different-duration-weight frame collision, and an
artifact-kind mismatch; cancellation refusal without artifact publication or
restored state; and preservation of an earlier checkpoint when a later store
refuses. A fresh sealed job can publish only generation zero. A persisted root
must remain byte-identical and independently readable after ledger reopen,
typed restore, state advance, automatic successor publication, and a second
reopen/restore. A separately serialized canonical generation that names the
root but regresses its committed state must refuse, as must an aggregate
successor-recovery budget one byte below the exact artifact pair. These are
replay, state-ancestry, and transactional storage
claims only; they do not establish concurrent scheduler-claim arbitration,
render convergence, 4K capacity, or scientific authority beyond the source
trajectory.

## Deterministic cinematic job conductor v1

The opt-in `cinematic-orchestration` feature owns the smallest film-specific
dependency graph spanning the pinned trajectory, independent uniform render
shards, per-segment raw merge, complete temporal-neighbor finishing frontier,
image-sequence seal, audio controls/excitation/resampling/modal synthesis/WAV,
independent bundle verification, and optional quarantined mux derivative. It is
not a general workflow engine. Node identities cover the exact configuration
partition, implementation/checker identity, stage-local typed inputs, logical
render work, and sorted dependency identities. The public constructor requires
the original `EulerCinematicScene` and scene-bound `EulerRenderFrameInput`
inventory used to build the render plan. It rechecks every prepared-segment
identity and admits every shutter through the scene's `AnimatedCamera`; the
resulting `CameraExposure::shot_id` is the sole camera-shot oracle. Temporal
neighbors are restricted to that same continuous shot, so a hard camera cut
cannot leak denoising history across the boundary. Bundle verification also
binds an independent finalization-expectation identity; changing it invalidates
only bundle verification and muxing. Resource ceilings do not change content
node identity, but the complete plan identity retains every active-node
work/output ceiling and graph limit for admission and replay.

Plan construction checks node, dependency, output, event, and snapshot bounds
before retaining dependency-sized graph storage. The canonical resume snapshot
is a bounded, self-hashed, fixed-width node-to-artifact index. A snapshot is
only a lookup hint: reuse additionally requires exact node, artifact-family,
expected-output, content, byte-hash, and byte-length bindings plus the stage
owner's independent decoder/checker. The backend discovery seam covers the
crash window after atomic artifact publication but before snapshot persistence.
A checked publication that differs from, or lacks, its exact snapshot hint is
marked recovered: the node remains complete, but every descendant is rebuilt
because v1 receipts do not yet carry dependency-publication content lineage. A
missing, corrupt, recovered, or rebuilt node therefore cannot leave a stale
downstream artifact in the current sequence or bundle.

The conductor polls the caller's `Cx` before and between discovery, staging,
checking, and publication. A checked candidate does not become complete until
atomic publication returns the identical bounded descriptor. Cancellation
after publication retains that exact record for resume; cancellation before
publication exposes no completion. Backend methods are synchronous structured
boundaries: any renderer/audio crews they start must be drained before return
or unwind. Cancellation or deterministic deadline expiry observed while a
backend returns an error takes precedence over classifying the node as failed,
so interruption remains resumable rather than becoming a false backend fault.
Typed failures and contained synchronous panics otherwise block descendants
while independent completed work remains reusable. Progress reports monotone
logical nodes, scene-derived shots, frames, event-delimited finished segments,
render shards, unique frame-segment tiles (complete only after every sample
shard), tile-sample cells, exact retained path counts, and audio stages, plus an
estimated remaining-work count. None is a wall-clock ETA or throughput promise.

G0/G3 tests pin canonical topology and parallel frontiers, snapshot codec
self-verification, checker-before-publish ordering, exact no-work replay,
camera/material/image, sound/audio, trajectory, mux-only, and
finalization-expectation invalidation; missing, corrupt, and changed-publication
rebuilding; oversized-snapshot rejection with live discovery; failure/panic
containment; retry; pre-stage/check/publish cancellation;
cancellation/deadline precedence; exact terminal resume; admission refusal;
sample-split tile progress; and domain progress. A real three-frame
`EulerUniformRenderPlan` fixture exercises the public scene/prepared-frame
authority checks, temporal-neighbor translation, public cancellation/resume,
exact replay, and one-short dependency preflight. The graph names real stage
contracts but does not claim that the still-separate
cinematic-AOV shard renderer, temporal denoiser, mux adapter, or final 4K
qualification exists. Those implementations must land and pass their own
owner-specific gates before this conductor can produce a film.

## Cinematic finalization gate v1

The opt-in `cinematic-finalization` feature provides the read-only terminal
integrity gate for the current uniform fixed-SPP sharded cinematic path. A
`CinematicFinalizationPlan` reconstructs the expected frame, segment, role,
format, dimensions, channels, sampling, source-lineage, configuration, build,
profile, clock, cut, and audio-event contract from the admitted composition,
quality profile, brief, `EulerUniformRenderPlan`, exact `EulerCinematicScene`,
original prepared exposures, and independently retained audio expectations.
The scene, not a caller-provided palette or scalar shot assertion, supplies the
sorted object/material palette and each segment's admitted camera exposure and
shot identity. Producer completion flags and the produced sequence inventory
are inputs to check, never the oracle that defines success.

The constructor accepts only a complete zero-based master range. A partial
quality-profile range is refused because pairing it with the brief's complete
WAV, clock endpoints, cuts, and event table would silently produce a false A/V
alignment claim; a future partial-range path needs its own typed range clock.
The configuration's render/audio budget references must bind the exact quality
profile identity and identity-schema version, its timeline must bind the exact
brief identity and version, and the render plan, sound configuration, scene,
trajectory, and prepared exposures must agree at their available typed
identity boundaries.

`verify_cinematic_bundle` consumes explicit persisted bytes and external pins;
it opens no path and mutates neither the bundle nor its artifacts. It strictly
decodes the finalized sequence, checks canonical order and complete inventory,
then verifies each frame's byte count and content hash before inspecting EXR or
PNG structure. Raw EXRs must match the exact dimensions, AOV layout and
custom-attribute set, finite payload constraints, and validity-bit domain
reconstructed by the plan; an unexpected vendor attribute is not silently
accepted. Final-diagnostic raw masters additionally prove the exact uniform
sample count from the per-pixel `samples` plane, one-based object/material IDs
within the scene-derived palette lengths, the reserved v2 validity bit clear,
and exact per-pixel validity relations for IDs, primary coverage, and sample
contribution. Profiles without those planes prove only their exact renderer
SPP declaration and may not promote themselves to final delivery. Derived artifacts must retain their expected
raw-source lineage. Independently pinned authority records must bind each
artifact to its expected role, source, transformation, configuration, and
disclosures. The audio side checks canonical manifest bytes and identity,
source-signal/channel-layout/mix identities, strict float32 WAV structure and
sample count, exact video/audio clocks, a synchronization marker at every cut
boundary (not a distinct cut-flag receipt), and
the ordered canonical resampled-event receipts. Sequence, audio-manifest, WAV,
alignment, event, and authority identities remain distinct; one digest cannot
stand in for another proof obligation.

The report has five closed dispositions: `Pass`, `Incomplete` for absent or
unfinished required evidence, `Refused` for explicit resource or cancellation
boundaries, `Corrupt` for malformed, noncanonical, truncated, or hash-
inconsistent bytes, and `Incompatible` for well-formed artifacts that disagree
with the independent plan. Verification order, first-divergence coordinate,
ranked repair advice, verified counts, retained no-claims, and the report's
domain-separated identity are deterministic. Only `Pass` for an admitted
`Final4k` target sets `final_delivery_eligible`; that bit is artifact-delivery
eligibility, not scientific, release, or aesthetic authority.

Aggregate and per-codec byte/cardinality ceilings are checked before
image-sized materialization, arithmetic is checked, authority envelopes have
a fixed 64-KiB parser ceiling in addition to the caller ceiling, palette strings are shared
between frame expectations, and image, hash, manifest, receipt, inventory, and
event traversals poll the supplied `Cx` at declared intervals; canonical sorts
are cardinality-capped and poll immediately before and after sorting.
Verifier-owned bundle indexes reserve fallibly and borrow artifact bytes rather
than duplicating them. Cancellation or an insufficient caller budget returns `Refused` without
changing input bytes or publishing a replacement artifact. Focused G0 tests
cover the canonical
contracts and boundary cases; G3 hostile fixtures mutate pins, lineage,
inventory, headers, payloads, clocks, cuts, and events; G4 covers bounded
resource refusal and cancellation; and G5 covers stable report bytes and
identity for equal inputs. Tiny real EXR/PNG/WAV/sequence E2E fixtures exercise
the cross-codec composition, but are not 4K throughput or visual-quality
evidence.

V1 reconstructs only the existing `EulerUniformRenderPlan`. That plan carries
neither the quality profile's adaptive stopping result nor a receipt proving
the requested number of temporal shutter samples. Consequently its production
constructor always returns `NonFinal`, including for a nominal 3840x2160
`Final4k` profile: a passing report is useful integrity evidence but
`final_delivery_eligible` remains false. V1 does not serialize an adaptive
film, admit an adaptive sharding plan, or certify adaptive stopping, estimator
convergence, temporal accumulation completion, or a per-pixel sampling
optimum. Palette material hashes are renderer material identities, not measured
material truth. Exact source-signal and event-receipt integrity also does not
prove that any event made an audible causal contribution to the encoded
waveform. Nor does the gate independently resolve the composition's abstract
camera/geometry asset references into the supplied scene; it binds both exact
inputs and verifies all identity relations currently exposed by their typed
APIs. More generally, byte integrity and internal compatibility do not
establish physical fidelity, acoustic validity, perceptual quality, material
calibration, experimental agreement, mechanism truth, or aesthetic approval.

## Offline cinematic denoising controls

`euler_cinematic_denoise` always treats the first requested frame as a temporal
history cut. A range beginning after frame zero is admitted only with the
explicit `--initial-cut` acknowledgement; continuity is then checked strictly
between every subsequent frame. The CLI exposes only the spatial à-trous pass
count and scene-linear RGB sigma. Reprojection, history weight/length, and
coverage, depth, normal, object, and material compatibility remain frozen
correspondence-safety controls rather than noise-strength knobs. Every output
is a biased display derivative. In particular, primary-surface guides on glass
do not certify the motion of reflected or refracted secondary transport, so a
smoother preview is not evidence of raw-estimator convergence or transport
correctness.

The decoder admits uniform, adaptive-stopping, and
`independent-pilot-fixed-v1` FinalDiagnostic sample planes. For an independent
pilot it validates the complete canonical pilot-plan header, requires every
retained count to remain within the declared minimum/maximum allocation, and
requires the exact count-plane sum to equal the frame's declared retained path
total. Temporal sequence identity binds the frame-invariant pilot seed,
sampler, thresholds, safety factor, and allocation bounds while deliberately
excluding that per-frame total. The discarded pilot observations are not
available in the EXR, so neither sample inspection nor denoising upgrades the
allocation into an image-error or convergence certificate.

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

The bounded campaign additionally makes no claim that its reduced point-contact
closure is resolved multiphysics or that its encoded decay exponents and
channel laws were identified from data. Its convergence receipt may support a
configuration ordering *inside that declared reduced numerical model* in
either direction; it cannot promote that result into a physical ranking of
sharp versus filleted edges, glass versus steel, rings versus discs, or any
other real configuration. It does not predict a physical spin time. Its
geometry, contact, base, decay, and exterior records are numerical/software
rungs, not experiment-, video-, or Mould-backed evidence.

The Euler-local port registry does not implement any gravity, contact,
partial-slip, rolling/contour/spin, impact, base, exterior-gas, or gas-film
physics. Its structural additive receipt is not a signed decomposition proof;
it neither authenticates a source nor establishes action/reaction balance. Its
ledger neither derives energy from forces or impulses nor closes an energy
window: unavailable and undeclared channels always retain an explicit
no-closure boundary. It does not satisfy the blocked generic `PortSchema`,
manifest, constraint/impact-law, RATTLE/generalized-alpha/nonholonomic, or DSR
lanes.
