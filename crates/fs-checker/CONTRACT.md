# CONTRACT: fs-checker

The standalone evidence-package checker (plan addendum, Proposal 12): an
independently distributable verifier — "don't trust us; here is the checker."

## Purpose and layer

Layer L6. Its sole direct dependency is `fs-package`; that package's production
cone contains `fs-evidence`, dependency-free `fs-blake3`, and the static
`fs-crosswalk` vocabulary. A HARD
distribution constraint (Proposal 12): NO solver stack, geometry kernel, or
license gate anywhere in the graph. By construction the checker cannot run a
solve. It carries `CHECKER_PROTOCOL_VERSION = 6` for the admission-receipt,
closed semantic-registry, and semantic-context ABI
(distributed independently). `CHECKER_SUPPORTED_PACKAGE_FORMAT = 8` is an
explicit protocol literal with a compile-time assertion against
`fs_package::FORMAT_VERSION`, so a package schema bump cannot silently retain
an incompatible checker ABI.

## Semantic identity contracts

Four owner-local typed-binary identities seal every checker authority surface:

- `fs-checker:decision-report` is identity version 8 under
  `fs-package:v8:checker-decision`. Its identity version is a literal tied by
  compile-time assertion to `CHECKER_SUPPORTED_PACKAGE_FORMAT`; the independent
  `CHECKER_PROTOCOL_VERSION = 6` is also fingerprinted and is not substituted
  for the package-format identity version. The decision binds every
  `CheckReport` and `Finding` field, nested semantic-context and verification
  receipt hashes, signature purpose/payload, collection length/order, and the
  derived decision hash rule.
- `fs-checker:semantic-plugin` remains identity version 1 under
  `fs-checker:semantic-plugin:v1`. Its bytes bind the exact family/schema,
  payload cap, implementation revision, family-specific arithmetic limits, and
  registry revision.
- `fs-checker:semantic-registry` remains identity version 1 under
  `fs-checker:semantic-registry:v1`. Its bytes bind the registry revision,
  implementation version, every global resource limit, and the exact compiled
  plugin count, order, descriptors, and fingerprints.
- `fs-checker:semantic-report` remains identity version 1 under
  `fs-checker:semantic-report:v1`. Its bytes bind the package root, registry,
  package status, aggregate charges, ordered claim/failure counts and rows, all
  nested receipt/failure fields, and every closed status/failure variant;
  `context_hash` is derived and never trusted as an independent input.

The three established v1 logical domains and canonical byte order are retained
unchanged. Retained decision, plugin, registry, and report digests are admitted
only at their exact declared identity version and only as exactly 32 bytes.
Stale, future, truncated, and extended transports refuse; there is no implicit
cross-version migration or best-effort reinterpretation. A migration must
re-verify the source artifact under the new implementation and emit a newly
versioned identity.

## Public types and semantics

- `check(&EvidencePackage) -> CheckReport` — re-verify a package with all
  external origin capabilities denied.
- `check_against_root(&EvidencePackage, expected_root) -> CheckReport` — also
  confirm the content address matches (tamper / substitution detection), with
  external origins still denied. A bounded-root refusal or mismatch terminates
  before any injected verifier can observe the package.
- `check_with_capabilities(package, expected_root, signature_verifier,
  capabilities)` — the in-memory entry point for explicitly authenticated
  source certificates, anchoring datasets, falsifier artifacts, derivation
  artifacts, waivers, and signatures. The separate signature argument selects
  the exact checker-purpose context and overrides a signature capability in the
  set.
- `check_json(...)` and `check_json_with_capabilities(...)` — strict schema-v8
  transport counterparts. Plain `check_json` denies external origins; the
  capability-aware form authenticates them after structural parsing.
- `check_release_preflight(&EvidencePackage, expected_root, verifier)` — a
  structurally non-admitting blocker inventory. It uses
  `CheckPolicy::ReleasePreflight`, always returns `Fail` with a
  `release-preflight-only` finding, and cannot become release authority if a
  future ungated origin is added. A bounded, structurally valid package still
  receives declaration-level falsifier, anchor, signature, and scientific-rank
  findings after its expected deny-all capability refusal. Malformed and
  oversized packages are not rescanned or amplified.
- `check_for_release_with_capabilities(...)` — release admission with explicit
  source-certificate, anchoring-dataset, falsifier, derivation, and waiver
  capabilities in addition to the mandatory signature verifier. It requires at
  least one scientifically admitted finite `Verified` interval or authenticated
  `Validated` claim. Half-infinite and whole-line `Verified` enclosures remain
  transportable no-claim artifacts but cannot satisfy this minimum; degenerate
  point-at-infinity intervals enclose no real value and no longer parse or
  verify at all (bead frankensim-extreal-program-f85xj.2.6). After proving the expected
  content root, the gate performs bounded structural inspection and refuses
  conclusive declaration-shape blockers before any external capability runs.
- `check_json_release_preflight(...)` — the non-admitting transport preflight.
- `check_json_for_release_with_capabilities(...)` — the strict-parser release
  entry point with explicit origin capabilities.
- `CheckReport` is sealed and exposes read-only verdict, bounded recomputed root,
  `IntegrityStatus`, `SemanticReport`, `OriginStatus`, breakdown, signature,
  receipt, findings, policy, expected root, and `decision_hash`. The three
  stages cannot be confused: semantic refusal leaves integrity `Verified` and
  origin `NotRun`; an origin refusal preserves the independently computed
  semantic transcript. The hash binds checker protocol, every stage status,
  semantic context, policy, expected root, package root/receipt, signature
  purpose, summary, verdict, and findings.
  `receipt` is `Some` only after successful package verification. `passed()` is
  policy-local; release consumers use `release_admitted()`, which additionally
  requires the ReleaseAdmission policy, receipt, admissible semantic status,
  authenticated origin, and valid decision/transcript hashes.
  `release_independently_verified()` additionally requires
  `SemanticStatus::Verified`; witnessless legacy packages may still satisfy
  `release_admitted()` with the explicit `NotProvided` status.
- `verify_portable_semantics(&EvidencePackage) -> SemanticReport` performs the
  callback-free semantic stage for producers and auditors. The sealed report
  exposes ordered per-claim decisions and failures, registry fingerprint,
  package root, aggregate resource charges, and `context_hash`. The exact
  context is required when producing a release-approval signature.
- `semantic_plugin_registry()` is a closed, compiled-in registry. Protocol v6
  implements exact family/version dispatch for
  `frankensim/exact-interval@1` and
  `frankensim/bounded-linf-residual@1`; unknown families and versions refuse.
- `Verdict { Pass, Fail }`;
  `SignatureStatus { Unsigned, Refused, Unverified, Authenticated(payload) }`;
  the authenticated payload has private fields and read-only accessors.
  `Finding { kind, detail }`;
  `SemanticStatus { NotProvided, Verified, Refused, NotRun }`;
  `SemanticFailureKind { StructuralIntegrity, UnknownFamily,
  UnsupportedVersion, MalformedPayload, ResourceLimit, ClaimMismatch,
  VerifierPanic }`.
- Any package-verification or origin-capability refusal carries a zeroed
  breakdown, so unauthenticated evidence cannot retain a normal-looking
  positive pie alongside the failure finding.
- Re-exports `EvidencePackage`, `ContentHash`, `ColorBreakdown`,
  `MagnitudeBudget`, `PackageError`, `ParseError`,
  `VerificationCapabilities`, `VerificationReceipt`, admission types, and all
  six verifier interfaces.

## What it re-verifies

The authority order is fixed and fail-closed:

1. The content address through bounded `try_merkle_root`, optionally checked
   against an expected value before any source, anchor, falsifier, derivation,
   waiver, or signature capability dispatch. A transport refusal uses a zero
   refusal sentinel in the sealed report and never hashes or clones rejected
   oversized bytes. A mismatched bounded package reports its actual and expected
   roots but carries no receipt or admitted breakdown.
2. Callback-free format, transport, content-binding, per-claim completeness,
   sealed origin/color consistency, receipt structure, and portable-witness
   envelope validation through `EvidencePackage::verify_structural_integrity`.
3. Every attached portable witness through the exact compiled family/version.
   The checker recomputes the mathematical result from canonical primitive
   inputs and requires bitwise equality with the claim's ordered finite
   `Verified` interval. Any false, unsupported, malformed, over-budget, or
   panicking plugin refuses the whole semantic stage and suppresses every
   external callback. Witnessless claims receive an explicit `NotProvided`
   decision; they are never described as independently verified.
4. Source-certificate, anchoring-dataset, falsifier-artifact,
   derivation-artifact, waiver, and signature decisions only through exact
   typed `VerificationCapabilities`. Plain integrity entry points use
   `deny_all()`. Package verification rechecks structural semantics, but no
   callback is reachable until stages 1–3 succeed.
5. Signature validity only through an injected `SignatureVerifier` over a typed
   purpose. Integrity uses `PackageRootAttestation`; release uses
   `ReleaseApproval { checker_protocol, expected_root, admission_context,
   semantic_context }`. The contexts separately bind every non-signature policy
   fingerprint, waiver day, admission, compact waiver edge, package root,
   plugin implementation/limits, and complete semantic transcript. Policy,
   clock, witness, plugin-set, or semantic-context replay changes the canonical
   signature subject and refuses. Purpose substitution refuses. No signer
   identity or role is inferred.
6. A policy-bound verification receipt: package root, policy fingerprints,
   waiver day, signature status, and ordered origin/admission/waiver decisions.
7. For explicit release admission only: at least one scientific finite
   `Verified` interval or authenticated `Validated` claim, purpose-bound
   approval, authenticated per-certificate falsifiers, and exact authenticated
   per-Validated anchors. Empty, unsigned, all-waiver-dependent, unpaired
   certificate-class, and unanchored `Validated` declarations are conclusive
   structural blockers: they return a zeroed breakdown, no receipt, and at most
   raw `Unverified` signature bytes without invoking a verifier. Structurally
   complete candidates may be rehashed by the package verifier, but every such
   pass remains behind the same bounded transport envelope and precedes or
   accompanies the exact authority decision it protects.

### Built-in canonical witness schemas

- `exact-interval@1` is a little-endian, iterative straight-line program with
  at most 4,096 nodes. Leaves are exact signed integers in the IEEE-754 exact
  `i53` range or finite ordered binary64 intervals. Closed tags implement add,
  subtract, multiply, divide, negate, and hull over strictly backward node
  references, followed by one in-range result index and no trailing bytes.
  Exact integer arithmetic remains exact while representable; other arithmetic
  expands endpoints by one adjacent binary64 value. Division by an interval
  containing zero refuses.
- `bounded-linf-residual@1` is a little-endian dense residual witness: norm tag
  zero, nonzero row/column counts, row-major finite binary64 `A`, finite `x`,
  and finite `b`, with no trailing bytes. The checker independently encloses
  `b - A*x` using local outward interval arithmetic and returns
  `[+0, max_i |r_i|]`. Each dimension is at most 128 and the matrix is at most
  16,384 entries.
- Per-witness payload is at most 256 KiB; a package has at most 4,096 witnesses,
  8 MiB aggregate witness bytes, and 1,000,000 charged primitive operations.
  Checker caps reuse the public `fs-package` envelope caps with compile-time
  equality assertions. All family-specific and aggregate limits plus the
  implementation version are bound into plugin/registry fingerprints.

## Invariants

- No solver / license in the build graph (enforced by the dependency set).
- A package that fails `verify_with` (incomplete claim, unsupported format, or
  refused capability) yields `Verdict::Fail` with a matching finding and a
  zeroed breakdown; a content-address mismatch fails.
- `check`, `check_against_root`, `check_json`, `check_with`, and both release
  preflight entry points use `VerificationCapabilities::deny_all()`.
  Certificate-shaped bytes and waiver fields never authorize themselves. Only
  the explicitly capability-bearing release functions can grant admission.
- Source-certificate verification receives the complete typed request: package
  provenance and recomputed root, claim index/id/subject hash, statement,
  interval, producer, artifact hash, and optional portable witness. Waiver
  verification receives the package-owned authorization message and an explicit
  date context. Anchored-source verification binds the exact validity regime,
  dataset identity, and parsed dataset hash, but — unlike the
  source-certificate request above — it does NOT receive the package root or
  the claim subject hash.
- Integrity, semantic, origin, and signature authority are separately visible
  and hash-bound. Semantic failure never invokes an external verifier; origin
  or signature refusal never erases a completed semantic transcript.
- Every semantic context binds the exact package root, registry fingerprint,
  ordered claim/content/witness/plugin identities, per-claim status and failure,
  operation charges, aggregate counters, and package status. Stored context
  hashes are revalidated before a checker decision is considered valid.
- Direct and transitively derived waiver-dependent claims appear only in the
  fourth pie bucket and never in scientific rank or magnitude summaries.
- Every successful package verification retains its `VerificationReceipt`;
  parse and capability refusals retain none. Rejected callback fingerprints and
  both identities in a fingerprint-drift event remain in the finding hashed by
  the checker decision.
- Signature verification is independent from scientific-origin verification.
  It is optional outside release admission and mandatory at release admission.
- An empty package verifies vacuously and renders a "no claims" pie.
- Release preflight never passes, even if every current blocker is absent. An
  empty package and all-estimated, all-waived, or solely vacuous-Verified
  packages never pass actual release admission; ordinary integrity, preflight,
  and admission are distinct hash-bound policies.
- Verified and Validated claims never pass release admission without
  authenticated, content-addressed falsifier artifacts. Validated claims
  additionally require an exact matching canonical dataset anchor authenticated
  against the anchor request. That request is NOT the complete package subject:
  it omits the package root and the claim subject hash, so an anchor decision
  is portable across packages sharing its fields. Package binding comes from
  the root and policy fingerprint in the receipt (bead
  frankensim-extreal-program-f85xj.2.13).
- A derived `Validated` declaration with no `Validated` parent fails before
  capability dispatch and renders the stable
  `validation-authority-promotion-refused` finding. Thus an accepting
  derivation callback, a small numerical residual, a domain intersection, or
  any number of `Verified` parents cannot mint model-form or experimental
  validation authority through the standalone checker.
- Oversized in-memory builders and expected-root mismatches are refused before
  every external verifier callback and before per-claim release diagnostics.
  Rejected oversized signature bytes are not retained; a bounded mismatched
  package may retain its detached bytes only as explicitly `Unverified`.
- Actual release admission dispatches no source, anchor, falsifier, derivation,
  waiver, or signature callback when a structurally inspectable declaration has
  a conclusive release-shape blocker. Malformed declarations still traverse the
  structural verifier for a precise refusal, which itself precedes callbacks.
- Release preflight distinguishes capability refusal from structural refusal:
  only `EvidencePackage::is_structurally_inspectable_unverified()` inputs are
  scanned for independent declaration-level blockers when no receipt exists.
- `render_pie`, reports, and decision hashes are deterministic; pie arithmetic
  widens counts to `u128` before multiplication.
- Canonical semantic parsers are iterative, reject trailing/truncated input,
  forward/self references, nonfinite residual scalars, invalid dimensions, and
  arithmetic singularities, and cannot allocate beyond validated bounds.

## Error model

The checker does not error — it REPORTS: failures become bounded `Finding`s in
a sealed `CheckReport` with `Verdict::Fail`. External verifier panics become
structured package refusals; a built-in semantic verifier panic becomes a
callback-free `VerifierPanic` refusal. No rejected transport is reserialized or
scanned again for release diagnostics.

## Determinism class

The registry order, iterative payload decoders, arithmetic, semantic transcript,
deny-all report, rendered pie, and checker decision hash are deterministic pure
functions of the package and explicit gate context on the same binary64
implementation. Capability-aware reports additionally depend on atomic verifier
decisions and explicit waiver date; reproducible deployments must use pinned
deterministic policies.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/checker.rs`, `tests/plugins.rs`, `tests/epi_e2e.rs`, plus crate unit tests
(Proposal 12): clean pass with no findings;
incomplete-validated-claim failure; content-address (Merkle) tamper detection;
including provenance tamper; malformed falsifier refusal with fail-closed pie;
signature-presence and verifier-capability reporting; deterministic budget-pie
rendering; empty-package pie; protocol version; determinism; and release-gate
admission/refusal batteries for empty, unpaired, unanchored, unsigned, and
wrong-root packages through both in-memory and strict JSON entry points. The
battery also checks positive and negative source-certificate and waiver
authentication through in-memory, JSON, release, and JSON-release paths;
ordinary positive paths remain unsigned, capability refusals zero the
breakdown, source verifiers bind the exact typed claim, and waiver verifiers
bind the complete package-owned authorization message. The battery also locks
all-estimated/all-waived release refusal, purpose-bound release signatures,
scientific-policy and waiver-clock replay refusal, structurally non-admitting
preflight policy, checker decision-hash mutation coverage, oversized-builder
diagnostic bounds, and zero-callback expected-root refusal across every injected
capability for both in-memory and strict-JSON integrity/release paths. A
six-case release-shape battery proves callback-free refusal for empty, unsigned,
all-waived, unpaired, and unanchored inputs, then proves a complete mixed-origin
candidate dispatches every capability with identical in-memory and JSON reports.
The semantic battery covers both positive built-in families through in-memory
and schema-v8 JSON paths, including hand-derived expected-bit goldens for
outward-rounded interval addition, multiplication, and division plus a nonzero
rectangular dense residual; witnessless status separation; false exact and
dense residual mathematics; unknown family/version; payload/root tamper with a
newly computed root; trailing, truncated, nonfinite, divide-by-zero,
forward/self reference, node/dimension/element overflow; exact aggregate
operation-budget overflow; callback suppression; report/context determinism;
and release signature replay under a different semantic context. The EPI
round-trip carries and independently recomputes a portable interval witness.

## Independent re-verification (bead qmao.6.1)

`check_json` is the deny-all third-party entry point: strict parse (root
recomputation, structural semantics, and budget re-derivation happen in the
parser), then semantic re-verification, optionally against an expected root and
a `SignatureVerifier` capability. `check_json_with_capabilities` adds explicit
source-certificate, anchoring-dataset, falsifier, derivation, waiver, and
signature authentication.
Signature validity is asserted only
when a supplied capability accepts the signature over the canonical typed
subject hash; the in-tree `NoSignatureVerifier` accepts nothing (the no-crypto
no-claim — presence is recorded as `Unverified`, and supplying a
capability that rejects raises a `signature-invalid` finding). The
magnitude budget must reconcile with its parts. The normal dependency graph is
`fs-package -> {fs-blake3, fs-crosswalk, fs-evidence -> fs-obs}`: it contains no
solver and the checker cannot run a solve by construction.

## No-claim boundaries

- Identity validation proves exact versioned byte binding only. It does not
  authenticate a signer, make a semantic refusal scientifically true, or
  upgrade package evidence. A stale identity is refused rather than migrated;
  refusal itself supplies no claim about the old producer beyond
  incompatibility with the current schema.
- This crate ships no cryptographic primitive or signer registry. It records
  only that the injected policy accepted an exact typed signature subject; it
  does not establish signer identity, organizational role, or authorship.
- Composition receipts are re-run, but the checker itself does not produce or
  fetch source certificates or anchoring datasets. Injected verifiers may
  retrieve and independently validate addressed artifacts; the checker only
  supplies exact typed subjects and fails closed without those capabilities.
- Schema v8 carries portable witnesses and binds release approval to their
  independent semantic transcript. It inherits schema v7's binding of the
  color-algebra version into derived receipts; the checker refuses any stale
  algebra before re-running composition. Schema v6
  sealed every claim behind a typed origin and emitted a policy-bound
  admission receipt. Content addressing proves
  package integrity, not scientific truth. Successful source verification
  means only that the caller's configured verifier accepted the exact artifact
  subject; successful waiver verification means only that the configured
  policy accepted that exact package context through the stated date.
  Waiver-dependent claims remain visible but never become scientific evidence.
- Release admission adds authenticated falsifier, anchor, and purpose-bound
  signature obligations but does not re-run the source solver or independently
  establish experimental quality. Preflight is an inventory only and is never
  release authority.
- `exact-interval@1` proves only the result of its supplied bounded expression
  under the documented binary64 enclosure rules; it does not prove that the
  expression models the claimed physical system. `bounded-linf-residual@1`
  proves only a dense finite binary64 L-infinity residual enclosure for the
  supplied `A`, `x`, and `b`; it does not establish conditioning, uniqueness,
  discretization error, convergence, or model fidelity. Origin authentication
  remains a separate mandatory authority boundary for source certificates.

## Public certificate-plugin facade (bead frankensim-checker-semantic-plugins-9e8n)

`fs_checker::plugins` is a discoverability namespace over the one authoritative
package-bound registry described above. It re-exports the exact descriptors,
limits, identity versions/domains, report types, registry fingerprints, and
`verify_portable_semantics` implementation used by ordinary checking, strict
JSON checking, and release admission. It owns no second registry, raw
`(family, version, bytes)` checker, producer encoder, arithmetic path, verdict
type, or payload-only hash.

The sole positive independent-semantic authority path is:

1. family-owned canonical bytes are enclosed in `SemanticWitness`;
2. `Claim::from_portable_certificate` binds the witness's domain-separated
   BLAKE3 identity into the source-certificate origin and claim declaration;
3. `EvidencePackage` binds the complete claim into its content root;
4. the closed registry recomputes every attached witness after structural and
   expected-root checks, sealing package root, claim/witness/plugin identities,
   resource charges, and refusals into `SemanticReport`;
5. origin policies independently authenticate the exact typed source subject;
6. release approval signs both the scientific admission context and semantic
   context.

Unknown families, unsupported schema versions, malformed or over-budget
payloads, arithmetic mismatches, and plugin panics refuse at step 4 before any
external capability runs. A public caller cannot construct a positive semantic
receipt or bypass package binding by checking detached payload bytes.

### Facade and identity no-claim boundaries

- Re-exporting a registry API does not create a second authority or relax the
  closed family/version dispatch. New semantics require a new versioned family
  and corresponding fingerprint/identity review; an existing parser is never
  widened silently.
- `SemanticWitness::content_hash`, plugin fingerprints, the registry
  fingerprint, and semantic context are domain-separated BLAKE3 integrity
  identities. They are not signatures, provenance, solver authentication, or
  scientific truth.
- The retained-identity admission helpers check only the declared identity
  version and exact 32-byte transport width. An exact-width foreign digest is
  still merely an opaque digest and gains no authority until compared with the
  recomputed compiled/package identity in its consuming protocol.
- A semantic refusal refutes the attached certificate under the compiled
  witness semantics; it need not refute the underlying mathematical statement.
  Conversely, semantic verification proves only the documented bounded
  recomputation, not model fidelity, conditioning, uniqueness, discretization
  error, convergence, material realism, or experimental validity.
- The facade adds no dependency. The checker remains solver-free, geometry-free,
  license-free, synchronous, safe Rust, and suitable for the same standalone
  and WASM distribution cone.

`tests/plugins.rs` independently pins hand-authored v1 canonical bytes and their
v8 witness, v1 plugin, and v1 registry identities; exercises both families
through `SemanticWitness -> Claim -> EvidencePackage -> strict JSON`; proves
ordinary and release-gate use of the same sealed semantic report; and covers
expected-root, family, schema, payload, registry-identity transport, and release
semantic-context substitution/refusal.
