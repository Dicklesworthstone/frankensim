# CONTRACT: fs-package

Machine-checkable evidence packages (plan addendum, Proposal 12): a
content-addressed bundle of color-typed claims a standalone checker can
re-verify without solvers.

## Purpose and layer

Layer L6. Depends on `fs-evidence` (UTIL — `Color`, `ColorRank`,
`ValidityDomain`), `fs-crosswalk` (the static standards vocabulary used by
coverage reports), and dependency-free `fs-blake3` (the shared content-hash
owner). The structural core performs no I/O and has no solver dependency.
Injected verifier implementations are caller-owned synchronous capabilities.

## Public types and semantics

- `Claim` — a sealed claim plus its epistemic color and typed `ClaimOrigin`.
  Construction is through `from_certificate`, `from_portable_certificate`,
  `anchored`, `estimated`, `derived`, or `waived`; callers cannot assemble an
  origin-free claim from a public `Color`. Claim ids are canonical,
  non-placeholder, and unique within a
  package; statements are non-blank and cannot be reserved placeholder text.
  Raw declarations are exposed only by explicitly named `*_unverified`
  accessors. Scientific callers consume `VerifiedPackage::admitted_claims`,
  where waiver-dependent descendants have no `scientific_color()`.
- Package-side color admission (bead 6pf9, stage S2) — `scientific_color()`
  still returns a raw, copyable `Color`; downstream positive-evidence
  consumers should instead convert through the admission bridge:
  `VerifiedPackage::claim_admission_receipt(claim_id)` mints an
  `fs_evidence::AdmissionReceipt` whose node identity is the claim's
  domain-separated declaration hash (`fs-package:v8:claim` surface,
  `CLAIM_DECLARATION_IDENTITY_VERSION`), refusing unknown, waiver-tainted
  (transitively), and non-positive claims, and refusing everything when the
  retained package/report binding no longer re-derives.
  `PackageColorAdmissionVerifier` is the injected
  `fs_evidence::AdmissionVerifier` that authenticates a (candidate, receipt)
  pair by re-deriving it from the retained pair: exact claim address, exact
  canonical color bytes, exact schema/algebra versions, this authority's own
  policy fingerprint, and a re-validated binding. The policy fingerprint
  (`package_color_admission_policy_fingerprint`) is deliberately distinct
  from the ledger authority's, so package-admitted and ledger-admitted
  colors remain separately auditable and receipts cannot cross authorities.
  NO-CLAIM: admission authenticates that the checker's policy admitted the
  claim inside this retained pair — it is capability injection, not
  cryptography, and it does not upgrade the origin/dataset trust recorded by
  the verification receipt.
- `SemanticWitness { family, schema_version, canonical_payload }` — a sealed,
  portable envelope whose family-owned canonical bytes can be interpreted by a
  standalone semantic-checker plugin. `content_hash()` binds field lengths,
  family/version, and payload through schema-v8 derive-key domains without
  copying the payload. Package structure authenticates the envelope address;
  only a plugin can admit its mathematical meaning.
- `Provenance { code_version, constellation_lock }`.
- `VerificationCapabilities` — explicit source-certificate, anchoring-dataset,
  falsifier-artifact, derivation-artifact, waiver, and detached-signature
  verification capabilities. `deny_all()` is the default. A
  `SourceCertificateVerifier` receives a typed `SourceCertificateRequest`
  containing the exact claim, package provenance and root, claim subject hash,
  index, producer, parsed artifact hash, and optional semantic witness. An
  `AnchoredSourceVerifier` receives the exact claim provenance, index,
  statement, validity regime, dataset id, and parsed dataset hash; its current
  request does not carry the package root or source-certificate subject hash. A
  `FalsifierVerifier` and `DerivationVerifier` receive the exact package root,
  claim subject, ordered records/parents, and parsed proof-artifact hash. A
  `WaiverVerification` pairs a `WaiverVerifier` with an explicit Unix-day clock
  context. Every callback atomically returns acceptance plus a stable
  `PolicyFingerprint` through a sealed `VerificationDecision`; rejecting and
  conflicting fingerprints survive in structured refusals.
- `EvidencePackage { format_version, claims, provenance, signature }` —
  builder: `new(prov).with_claim(..).signed(..)`.
  - `try_merkle_root() -> Result<ContentHash, PackageError>` — a bounded 32-byte
    BLAKE3 Merkle root over the
    package identity: format version, claim count, provenance, and ordered
    claims. Header, claim, and internal-node hashes use standard BLAKE3
    derive-key domains. Detached signatures are excluded. Any reproducibility
    provenance or claim change changes the root.
  - `verify()` uses `VerificationCapabilities::deny_all()`. Only empty packages
    and ungated Estimated-source claims without falsifier records can pass.
    Source certificates, anchored sources, derivations, waivers, and attached
    falsifier artifacts all require their exact capability. Detached signature
    bytes remain `Unverified` without a signature capability.
  - `verify_with(&capabilities)` performs structural verification and then
    authenticates every capability-gated origin before returning a positive
    `PackageReport`.
  - `verify_structural_integrity()` is the callback-free boundary: it checks
    schema, transport limits, root-bound claim structure, origin/witness
    consistency, and returns the recomputed root without granting scientific
    authority.
  - `color_breakdown()` and `color_breakdown_with(..)` return a budget pie only
    through their corresponding verification path.
  - `waiver_message(index) -> Result<Vec<u8>, PackageError>` constructs bounded,
    package-owned, domain-separated authorization bytes; `with_waiver_mac`
    installs the final authenticator without changing those bytes.
  - `to_json() -> Result<String, PackageError>` is bounded and deterministic.
    `from_json()` reads the deterministic FrankenSim JSON profile, checks
    structure, declared budgets, and the root while retaining unauthenticated
    origins; `from_json_with()` returns a `VerifiedPackage` after authenticating
    every gated artifact.
- `ColorBreakdown { verified, validated, estimated, waived }` counts admitted
  scientific claims in the first three buckets. A direct waiver and its entire
  derived descendant cone appear only in `waived`.
- `Claim::declared_is_release_scientific_evidence_unverified()` distinguishes
  finite informative `Verified` intervals and `Validated` claims from vacuous
  infinite enclosures. It is a raw shape predicate only; release authority still
  requires the package receipt and checker gate.
- `PackageReport` and `VerificationReceipt` have sealed fields and read-only
  accessors. The receipt binds the package root, every invoked policy
  fingerprint, waiver day, signature status, and ordered `ClaimAdmission`
  decisions into its own domain-separated `receipt_hash`. Direct waiver ids are
  interned once; descendants store immediate waiver-dependent parent indices,
  preventing quadratic string amplification while preserving the complete DAG.
- `AuthenticatedSignature` has private fields and read-only accessors, so safe
  downstream code cannot substitute a different signature or purpose into an
  authenticated status. Positive status has authority only inside a sealed
  receipt.
- `PackagePresenceReport` and `PackageCoverageReport` retain the package receipt
  and carry their own decision hashes. They do not implement consuming
  `IntoIterator` or slice `Deref`. Bare row booleans/statuses are explicitly
  non-authoritative when detached; authority requires validating the sealed
  report hash and receipt.
- `PackageError` — structured refusals for incomplete provenance, invalid or
  duplicate claim ids, blank/placeholder claim statements, malformed color
  payloads, unsupported formats, receipt mismatches/parents, malformed
  falsifier/anchor records, refuted claims, transport limits, rejected policies,
  and policy-identity drift. It implements `Display` and `Error`.

## Invariants

- COMPLETENESS: reproducibility provenance fields, origin identities, and claim
  ids are canonical machine identities: no padding and no reserved placeholder
  tokens. Claim ids and waiver ids are unique within their namespaces, and
  claim statements are meaningful rather than blank
  or one of the reserved placeholders (`TODO`, `TBD`, `placeholder`, `N/A`/`NA`,
  `none`, `not run`, `pending`, `unknown`, `-`, or `?`, case-insensitive). A
  `Validated` claim must have a non-empty regime
  (`regime.bounds()` non-empty) whose axis names are non-blank and whose bounds
  are finite and ordered, plus a non-blank anchoring `dataset`. A `Verified`
  claim must carry an ordered, non-NaN `[lo <= hi]` interval. Ordered infinite
  endpoints are a sound but vacuous enclosure, not a decision-grade bound. An `Estimated` claim needs a
  non-blank estimator identity and a non-negative, non-NaN dispersion.
  Positive infinity is preserved as the lower-layer algebra's explicit
  no-quantitative-spread-claim sentinel; it is distinct from finite subtotal
  overflow, which verification rejects. An honest all-estimated package
  remains valid.
- FALSIFIER EVIDENCE: a record has a non-blank, non-placeholder stable
  falsifier identity, at least one executed attempt, and a non-blank,
  non-placeholder outcome detail. A refuted record still rejects its claim and
  package. This structural rule does not assert that the recorded work ran;
  source authentication remains a no-claim below.
- DATASET ANCHORS: every attached anchor has a non-blank stable dataset id and
  an exactly 64-character, lowercase hexadecimal content hash. Crosswalk
  anchoring coverage for an `AnchoredSource` requires an attached record whose
  dataset id AND content hash exactly equal the origin tuple, and whose dataset
  id equals the `Validated` color. An unrelated canonical hash does not count.
  Positive admission additionally requires an injected anchor verifier to
  accept the complete typed subject; a matching hash is only an address. A
  derived `Validated` claim must carry at least one exact dataset-id match, and
  every such matching anchor is authenticated independently. Derivation
  authority cannot substitute a new dataset hash. The invoked anchor-policy
  fingerprint and root-bound anchor list participate in the receipt hash and
  release-admission context.
- SOURCE CERTIFICATES: a canonical certificate hash is only an artifact
  address. Positive verification requires an injected verifier to establish
  the exact typed claim request. Merely naming a producer and 64-hex hash never
  produces a report or coverage. The request's source-specific subject hash
  omits every external artifact address, including that certificate hash, so a
  content-addressed certificate can embed the subject without a fixed-point
  cycle. Portable subjects retain the semantic family/schema identity but omit
  the witness payload address; the typed request and package root still bind
  the full witness and complete declaration.
- PORTABLE SEMANTIC WITNESSES: only a `Verified` source-certificate claim may
  carry one. Its source `certificate_hash` must exactly equal the witness
  content hash. Family identity is canonical and at most 128 UTF-8 bytes;
  schema version is positive; payload is nonempty and at most 256 KiB. A
  package carries at most 4,096 witnesses and 8 MiB of aggregate decoded
  payload, with checked arithmetic. JSON requires `semantic_witness` on every
  claim (`null` or the exact closed object `family`, `schema_version`,
  `payload_hex`); payload hex is lowercase, even-length, and decoded limits are
  checked incrementally before payload allocation.
- WAIVERS: authorization requires an injected verifier and explicit date. The
  MAC message binds a domain tag, package provenance, ordered authorization
  context, target index, complete target claim, waiver id, and expiry. It
  excludes detached signatures and MAC fields so the message is stable before
  installing the authenticator. Expired, replayed, duplicated, or rejected
  grants fail closed. Authentication records an administrative exception; it
  never promotes the waived color. Every derived child with any
  waiver-dependent parent stores the immediate tainted-parent edge and remains
  in the `waived` bucket; direct ids live once in the receipt waiver registry.
  Such claims contribute nothing to scientific rank,
  certificate coverage, regime/dataset coverage, or numeric magnitude totals;
  `MagnitudeBudget::waived_unquantified` counts them explicitly.
- SIGNATURES: raw detached bytes are `Unverified`. An installed verifier must
  authenticate the canonical `signature_subject_hash` for the typed purpose or
  verification refuses. Generic package-root attestation is integrity evidence
  only. Signature coverage requires a `ReleaseApproval` purpose bound to an
  explicit checker protocol, expected root, and `release_admission_context`
  covering every non-signature policy fingerprint, waiver day, admission, and
  compact waiver edge, plus a separate caller-supplied `semantic_context` that
  identifies the approved semantic checker/plugin set. Producers obtain the
  admission context from an unsigned
  verification receipt, sign the public canonical subject hash, attach the
  bytes, and run the final gate. Coverage records authenticated signature intent;
  it does not establish that `fs-checker` admitted the package. No signer
  identity, role, or authorization is claimed.
  A release purpose naming any root other than the recomputed package root is
  refused before the callback, even if a permissive verifier would accept it;
  an approval from another policy, waiver day, or semantic context has a
  different subject hash.
- POLICY RECEIPTS: invoked verifier fingerprints and waiver day are decision inputs,
  not ambient process state. They, the root, signature result, origin class,
  admission class, compact waiver registry, and immediate taint edges are bound
  into the receipt hash.
- CONTENT-ADDRESSING: `try_merkle_root` is deterministic and tamper-evident across
  format version, claim count, provenance, and ordered claims. Derive-key modes
  separate package leaves/nodes from plain artifact hashes and each other; a
  detached signature does not change the root.
- Coverage is capability-aware. The plain coverage functions use deny-all
  verification and suppress every concept when a gated origin is
  unauthenticated; `_with` variants accept explicit capabilities. Color rank,
  certificate, regime, falsifier, and dataset coverage use only scientifically
  admitted claims. `Certificate` requires admitted `Verified` evidence with
  finite endpoints (an all-estimated or solely vacuous-Verified package has
  none). `ClaimOrigin` excludes Estimated declarations,
  waived claims, and waiver-dependent descendants. `Signature` requires
  release-purpose authentication but does not imply checker release admission.
  `WaiverAuthorization` remains an explicit administrative concept.

## Error model

Structured `PackageError` values (refusals that teach), never panics. Public
root computation, JSON serialization, and waiver-message construction are
fallible at the transport boundary. Panics from injected verifier callbacks are
caught and converted to structured refusals. The JSON
parser maps structural refusals into `ParseError`; `from_json_with` additionally
maps capability refusals. Plain `from_json` is intentionally an integrity and
structure boundary, not an authentication verdict. The untrusted JSON
boundary is bounded before schema mapping: 64 MiB input, depth 64, one million
values, 100,000 members per container, 1 MiB decoded strings, and 128-byte
number tokens. In-memory verification enforces the corresponding transport
envelope before a package can pass, so a verified object remains serializable
and checkable under those bounds. Accounting aborts at the first exceeded
byte/node budget, including inside one large claim. Semantic decoded-byte and
witness-count budgets are enforced separately and incrementally. Limit violations are
structured refusals.

## Determinism class

The fallible Merkle root, JSON profile, strict parsing, structural verification, and deny-all
verification are deterministic pure functions of the package (bit-exact on
float certificate payloads via `to_bits`). `verify_with` additionally depends
on caller-supplied capability decisions; reproducible deployments must provide
deterministic verifiers over pinned artifact stores, stable policy
fingerprints, and an explicit date.

## Cancellation behavior

No internal asynchronous work. Injected verifiers are synchronous callbacks;
their latency and cancellation behavior are outside this crate's control and
must be bounded by the caller's implementation.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/package.rs` (Proposal 12): complete mixed-color package with injected
source verification;
all-estimated boundary (valid + round-trips); validated-missing-regime and
validated-missing-dataset completeness failures; verified bad-interval
failure; blank/placeholder statement and falsifier refusal; Merkle determinism
and claim/provenance tamper detection;
unsupported-format rejection; optional detached signature; deterministic JSON
carrying the root; in-memory/serialized semantic parity; exact full-width
falsifier attempt-count round trips with overflow refusal; deny-all source and
waiver behavior; exact typed source-subject checks; package-context waiver
replay resistance; duplicate waiver ids; capability-aware coverage; exact
falsifier/derivation artifact subjects; policy-fingerprint drift; signed-zero
receipt identity; compact waiver-taint DAGs; sealed verified-package binding;
portable-witness construction, request binding, strict JSON, tamper detection,
shape and aggregate caps; semantic-context signature substitution; stale-v7
refusal; and origin transport/identity limits. The compile-fail battery proves that an
authenticated signature payload cannot be constructed downstream.

## Schema v8: portable semantic witnesses and semantic release context

Format 8 rotates every active `fs-package:v8:*` derive-key domain. Claims add a
mandatory JSON `semantic_witness` field and bind the witness content address
into the canonical body, Merkle root, origin certificate hash, authorization
context, source-verifier request, and transport accounting. Witness payloads
remain opaque to this crate; semantic plugins own interpretation and admission.
Release-approval signatures now bind `semantic_context` independently from the
package-computed origin-policy `release_admission_context`. Format 7 and earlier
are refused by the one-version contract.

### Schema-v8 identity owners

The crate owns the following closed identities. Every public domain constant
names the exact primary `fs-package:v8:*` context consumed by its encoder;
secondary family, payload, node, and authorization-context constants likewise
name complete contexts. Declaring these existing contexts does not change any
schema-v8 digest or authorization bytes.

| Identity | Canonical encoder and exact context | Semantic source | Intentional exclusions | Schema dependencies |
| --- | --- | --- | --- | --- |
| `semantic-witness` | `SemanticWitness::content_hash`; `fs-package:v8:semantic-witness`, `...:semantic-witness-family`, `...:semantic-witness-payload` | family UTF-8 and length, witness schema, payload bytes and length | none | none |
| `claim-declaration` | `Claim::declared_content_hash_unverified`; `fs-package:v8:claim` | complete claim declaration, including every artifact address and witness address | none | `semantic-witness` |
| `claim-verification-subject` | `Claim::declared_verification_subject_hash_unverified`; `fs-package:v8:claim-verification-subject` | claim subject and witness address | receipt/falsifier artifact addresses and waiver MAC output | `claim-declaration` |
| `source-certificate-subject` | `Claim::declared_source_certificate_subject_hash_unverified`; `fs-package:v8:source-certificate-subject` | claim/color metadata, address-free receipt/falsifier/attached-anchor metadata, portable family/schema, claim origin with source-certificate address and waiver MAC omitted | source-certificate, receipt, falsifier, attached-anchor, and witness payload addresses; witness payload; waiver MAC output | none; shared canonicalization and BLAKE3 closure are fingerprinted directly |
| `package-root` | `EvidencePackage::try_merkle_root`; `fs-package:v8:header`, dependent claim leaves, and `fs-package:v8:node` | format, claim count/order, declaration hashes, provenance, odd-node carry rule | detached package signature | `claim-declaration` |
| `waiver-authorization-subject` | `EvidencePackage::waiver_message`; exact `fs-package:v8:waiver-authorization\|...` bytes plus `fs-package:v8:authorization-context` digest | package header/provenance, ordered address-free claims, target index/body, waiver id, expiry | detached signature and every waiver MAC output | `claim-declaration` |
| `signature-subject` | `signature_subject_hash`; `fs-package:v8:signature-subject` | package root, purpose tag, and every release-purpose axis | detached signature bytes | `package-root`, `release-admission-context` |
| `verification-receipt` | `VerificationReceipt::recomputed_hash`; inner tag and outer `fs-package:v8:verification-receipt` context | root, all six policy slots, waiver day, signature status/purpose, ordered admissions and waiver registry | stored receipt hash is derived | `package-root` |
| `release-admission-context` | `VerificationReceipt::release_admission_context`; `fs-package:v8:release-admission-context` | root, five non-signature policy slots, waiver day, ordered admissions and waiver registry | signature status, signature policy slot, and stored receipt hash | `package-root`; provisional receipt encoding is fingerprinted directly |
| `presence-decision` | sealed presence-report digest; `fs-package:v8:presence-decision` | row count/order, concept, presence bit, rationale, optional receipt hash | stored decision hash is derived | `verification-receipt` |
| `coverage-decision` | sealed coverage-report digest; `fs-package:v8:coverage-decision` | standard, crosswalk/package versions, row count/order, concept, status, rationale, optional receipt hash | stored decision hash is derived | `verification-receipt` |

Retained digest transport is exactly 32 bytes paired out of band with identity
version 8. A stale version or any short/extended digest fails closed. Waiver
authorization is an exact canonical byte subject rather than a fixed-width
digest: retention requires version 8, byte equality with a freshly recomputed
subject, and the exact schema-v8 prefix. Package JSON is retained only through
the strict closed parser; accepted canonical JSON is a textual fixed point under
decode then encode. None of these integrity admissions grants artifact trust,
scientific truth, waiver authority, or signer authority.

## Schema v7: algebra-versioned derivations (historical)

Format 7 rotates current `fs-package:v7:*` domains and adds the exact
`fs-evidence` color-algebra version to every composition receipt. Both strict
transport parsing and in-memory verification refuse a stale algebra before
re-deriving a color. This prevents a receipt created under older identity,
rounding, or composition rules from being interpreted as current evidence.
Formats 6 and earlier are refused by the one-version contract.

## Schema v6: admission receipts and non-launderable waivers (historical)

Format 6 uses `fs-package:v6:*` domains and a closed JSON magnitude shape with
`waived_unquantified`. Source certificates, anchoring datasets, falsifier
artifacts, derivation artifacts, waivers, and signatures are explicit
fingerprinted capabilities over exact typed subjects. Artifact hashes are
addresses, never self-authentication. Verification emits a policy-bound receipt
with a topological admission decision for every claim. Direct waiver identities
are interned once and descendants retain immediate tainted-parent edges. Direct
waivers and every one-parent, multi-hop, or multi-parent descendant remain
waiver-dependent and cannot enter scientific color/magnitude/coverage summaries.
Format 6 and earlier are no longer accepted by current readers.

## Schema v5: sealed origins and explicit capabilities (historical)

Every claim carries one content-addressed `ClaimOrigin`. `AnchoredSource` and
`EstimatedSource` must agree exactly with the color; `Derived` must be paired
with a backward-only composition receipt that re-derives the exact color.
`SourceCertificate` and `AuthenticatedWaiver` are capability-gated. Plain
verification, breakdown, and coverage are deny-all, so certificate-shaped bytes
cannot create a positive result. Waiver messages are package-owned and
domain-separated, and source verifiers receive typed requests rather than raw
strings. Origins and all their strings are included in the in-memory transport
envelope. Schema v6 supersedes these domains and capability semantics.

## Schema v4: mode-separated BLAKE3 roots (beads 7uq9, t7x3)

The content address is a 32-byte `ContentHash`. A header leaf binds the format,
ordered claim count, and reproducibility provenance; each claim is a separate
canonical leaf; internal nodes combine exact child bytes. Header, claim, and
node hashing use distinct standard BLAKE3 derive-key contexts, which also
separate package identities from plain ledger artifact hashes. Odd tree nodes
carry upward unchanged, while the header's claim count prevents tree-shape or
duplicate-tail ambiguity. Readers refuse v3's 16-hex FNV root by version and
root width; checker protocol v2 is the matching standalone ABI.

## Schema v3: receipts, falsifiers, anchors (bead xfxq)

Claims optionally carry a COMPOSITION RECEIPT (parent claim indices +
the ledger op + the exact color-algebra version): `verify` re-runs
`fs_evidence::compose` over the parents in order and requires exact IEEE-754
bit identity (including signed zero). Readers and in-memory verification refuse
receipts from any other algebra version before interpreting their derived
identity namespace. A claim whose
color cannot be re-derived is `ReceiptMismatch`, so laundering a
Verified from Estimated parents is caught SEMANTICALLY, not just by
the content address; parents must point strictly backwards
(`BadReceiptParent` otherwise — a DAG by construction). FALSIFIER
RECORDS (name, attempts, refuted, detail) travel with the claim; any
`refuted: true` fails verification outright (`RefutedClaim`). ANCHOR
RECORDS give validated claims content-hashed dataset identities. All
three field families bind into the content address and round-trip
through the strict parser (booleans added to the closed grammar);
crosswalk coverage now reads validated falsifier logs and matching anchors from
the actual fields. Coverage fails closed for an invalid package, and raw
detached-signature presence never counts as authenticated sign-off. The
checker stays solver-free: `compose` lives in fs-evidence, already in its
dependency cone. Migration: v2 readers refuse v3 by version (the one-version
contract); in-tree constructors gained fields with builders, no call-site
changes.

## Schema v2: round trip and fail-closed parsing (bead qmao.6.1)

Fallible `to_json` emits the COMPLETE color payloads (floats as bit-exact hex
bits), provenance, signature, magnitude budget, and the content root;
`from_json` is a STRICT parser — unknown fields, missing fields, wrong
types, duplicate keys, bad hex, NaN certificates, inverted
intervals, negative dispersions, and unknown color kinds each refuse
with a structured `ParseError`. The parser re-derives the magnitude
budget from the parsed claims and RECOMPUTES the content root from the
parsed fields: a package whose embedded root does not recompute is
tampered in transit and never loads. This is an integrity check, not the
schema-v8 external-origin admission (inherited from v6) named above.
Decode-encode is both
semantically and textually stable (tested). Ordered infinite Verified endpoints
are accepted as rigorous but vacuous enclosures; coverage and magnitude
accounting treat them as unbounded, never as a finite decision-grade result.
The magnitude budget
attributes ERROR MAGNITUDES (verified interval widths, estimated
dispersions) and counts validated claims as unquantified regional
trust — never numerified. JSON number tokens are retained as decimal text until
they are converted into their target integer type; full-width `u64` falsifier
attempt counts therefore round-trip exactly, and overflow refuses instead of
rounding through `f64`. The reader implements the deterministic FrankenSim JSON
profile emitted by `to_json`; it does not claim to accept every optional escape
spelling available in the general JSON grammar.

## No-claim boundaries

- The Merkle hash is the shared in-house BLAKE3 implementation, with standard
  derive-key mode separation for typed package roots. This establishes content
  integrity, not authorship or scientific truth. A cryptographic SIGNATURE is
  detached and OPTIONAL. Signature decisions are made only by explicitly
  injected policy; the crate ships no default signing primitive or signer
  identity model.
- The crate does not fetch artifacts or choose trust roots. Callers provide
  source, anchor, falsifier, derivation, signature, and waiver verifiers; those capabilities may
  retrieve and re-check addressed artifacts, or refuse. `fs-package` supplies
  typed exact requests, policy-bound receipts, and deny-all defaults.
- The certificate payloads live in `fs-evidence::Color`; this crate bundles
  and content-addresses them, it does not produce them.
- A validated dataset hash proves content identity, not experimental quality or
  custodial authenticity. Those stronger properties require an external
  evidence policy. Likewise, successful waiver authentication proves only that
  the configured policy accepted the exact context through the stated date; it
  does not identify an authorized person or convert the waived claim into an independently
  reproduced result.
