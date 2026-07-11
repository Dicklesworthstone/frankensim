# CONTRACT: fs-checker

The standalone evidence-package checker (plan addendum, Proposal 12): an
independently distributable verifier — "don't trust us; here is the checker."

## Purpose and layer

Layer L6. Its sole direct dependency is `fs-package`; that package's production
cone contains `fs-evidence`, dependency-free `fs-blake3`, and the static
`fs-crosswalk` vocabulary. A HARD
distribution constraint (Proposal 12): NO solver stack, geometry kernel, or
license gate anywhere in the graph. By construction the checker cannot run a
solve. It carries `CHECKER_PROTOCOL_VERSION = 4` for the schema-v6 admission
receipt ABI (distributed independently). `CHECKER_SUPPORTED_PACKAGE_FORMAT = 6` is an
explicit protocol literal with a compile-time assertion against
`fs_package::FORMAT_VERSION`, so a package schema bump cannot silently retain
an incompatible checker ABI.

## Public types and semantics

- `check(&EvidencePackage) -> CheckReport` — re-verify a package with all
  external origin capabilities denied.
- `check_against_root(&EvidencePackage, expected_root) -> CheckReport` — also
  confirm the content address matches (tamper / substitution detection), with
  external origins still denied.
- `check_with_capabilities(package, expected_root, signature_verifier,
  capabilities)` — the in-memory entry point for explicitly authenticated
  source certificates, anchoring datasets, falsifier artifacts, derivation
  artifacts, waivers, and signatures. The separate signature argument selects
  the exact checker-purpose context and overrides a signature capability in the
  set.
- `check_json(...)` and `check_json_with_capabilities(...)` — strict schema-v6
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
  least one scientifically admitted Verified or Validated claim.
- `check_json_release_preflight(...)` — the non-admitting transport preflight.
- `check_json_for_release_with_capabilities(...)` — the strict-parser release
  entry point with explicit origin capabilities.
- `CheckReport` is sealed and exposes read-only verdict, bounded recomputed root,
  breakdown, signature, receipt, findings, policy, expected root, and
  `decision_hash`. The hash binds checker protocol, policy, expected root,
  package root/receipt, signature purpose, summary, verdict, and findings.
  `receipt` is `Some` only after successful package verification. `passed()` is
  policy-local; release consumers use `release_admitted()`, which additionally
  requires the ReleaseAdmission policy, receipt, and valid decision hash.
- `Verdict { Pass, Fail }`;
  `SignatureStatus { Unsigned, Refused, Unverified, Authenticated(payload) }`;
  the authenticated payload has private fields and read-only accessors.
  `Finding { kind, detail }`.
- Any package-verification or origin-capability refusal carries a zeroed
  breakdown, so unauthenticated evidence cannot retain a normal-looking
  positive pie alongside the failure finding.
- Re-exports `EvidencePackage`, `ContentHash`, `ColorBreakdown`,
  `MagnitudeBudget`, `PackageError`, `ParseError`,
  `VerificationCapabilities`, `VerificationReceipt`, admission types, and all
  six verifier interfaces.

## What it re-verifies

1. Format support, per-claim completeness, sealed origin/color consistency, and
   receipt re-derivation (delegated to `EvidencePackage::verify_with` — no
   solver).
2. Source-certificate, anchoring-dataset, falsifier-artifact,
   derivation-artifact, waiver, and signature decisions only through exact
   typed `VerificationCapabilities`. Plain integrity entry points use
   `deny_all()`.
3. The content address through bounded `try_merkle_root`, optionally checked
   against an expected value. A transport refusal uses a zero refusal sentinel
   in the sealed report and never hashes or clones rejected oversized bytes.
4. Signature validity only through an injected `SignatureVerifier` over a typed
   purpose. Integrity uses `PackageRootAttestation`; release uses
   `ReleaseApproval { checker_protocol, expected_root, admission_context }`.
   The context binds every non-signature policy fingerprint, waiver day,
   admission, and compact waiver edge; policy or clock replay changes the
   canonical signature subject and refuses. Purpose substitution refuses. No
   signer identity or role is inferred.
5. A policy-bound verification receipt: package root, policy fingerprints,
   waiver day, signature status, and ordered origin/admission/waiver decisions.
6. For explicit release admission only: non-vacuity, at least one scientific
   Verified/Validated claim, purpose-bound approval, authenticated
   per-certificate falsifiers, and exact authenticated per-Validated anchors.

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
  provenance, claim index and id, statement, interval, producer, and artifact
  hash. Waiver verification receives the package-owned authorization message
  and an explicit date context. Anchored-source verification additionally binds
  the exact validity regime, dataset identity, and parsed dataset hash.
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
  empty package and all-estimated or all-waived packages never pass actual
  release admission; ordinary integrity, preflight, and admission are distinct
  hash-bound policies.
- Verified and Validated claims never pass release admission without
  authenticated, content-addressed falsifier artifacts. Validated claims
  additionally require an exact matching canonical dataset anchor authenticated
  against the complete typed subject.
- Oversized in-memory builders are refused before root/signature canonicalization
  and before per-claim release diagnostics. Rejected raw signature bytes are not
  retained in a refusal report.
- Release preflight distinguishes capability refusal from structural refusal:
  only `EvidencePackage::is_structurally_inspectable_unverified()` inputs are
  scanned for independent declaration-level blockers when no receipt exists.
- `render_pie`, reports, and decision hashes are deterministic; pie arithmetic
  widens counts to `u128` before multiplication.

## Error model

The checker does not error — it REPORTS: failures become bounded `Finding`s in
a sealed `CheckReport` with `Verdict::Fail`. External verifier panics become
structured package refusals. No rejected transport is reserialized or scanned
again for release diagnostics.

## Determinism class

The deny-all report, rendered pie, and checker decision hash are deterministic
pure functions of the package and explicit gate context. Capability-aware
reports additionally depend on atomic verifier decisions and explicit waiver
date; reproducible deployments must use pinned deterministic policies.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/checker.rs` plus crate unit tests (23 cases, Proposal 12): clean pass with no findings;
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
preflight policy, checker decision-hash mutation coverage, and oversized-builder
diagnostic bounds.

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

- This crate ships no cryptographic primitive or signer registry. It records
  only that the injected policy accepted an exact typed signature subject; it
  does not establish signer identity, organizational role, or authorship.
- Composition receipts are re-run, but the checker itself does not produce or
  fetch source certificates or anchoring datasets. Injected verifiers may
  retrieve and independently validate addressed artifacts; the checker only
  supplies exact typed subjects and fails closed without those capabilities.
- Schema v6 seals every claim behind a typed origin and emits a policy-bound
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
