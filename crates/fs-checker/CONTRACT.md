# CONTRACT: fs-checker

The standalone evidence-package checker (plan addendum, Proposal 12): an
independently distributable verifier — "don't trust us; here is the checker."

## Purpose and layer

Layer L6. Its sole direct dependency is `fs-package`; that package's production
cone contains `fs-evidence`, dependency-free `fs-blake3`, and the static
`fs-crosswalk` vocabulary. A HARD
distribution constraint (Proposal 12): NO solver stack, geometry kernel, or
license gate anywhere in the graph. By construction the checker cannot run a
solve. It carries `CHECKER_PROTOCOL_VERSION = 2` for the schema-v4 32-byte-root
ABI (distributed independently). `CHECKER_SUPPORTED_PACKAGE_FORMAT = 4` is an
explicit protocol literal with a compile-time assertion against
`fs_package::FORMAT_VERSION`, so a package schema bump cannot silently retain
an incompatible checker ABI.

## Public types and semantics

- `check(&EvidencePackage) -> CheckReport` — re-verify a package.
- `check_against_root(&EvidencePackage, expected_root) -> CheckReport` — also
  confirm the content address matches (tamper / substitution detection).
- `check_for_release(&EvidencePackage, expected_root, verifier) -> CheckReport`
  — the stronger no-falsifier-no-ship admission gate: requires a non-empty
  package, authenticated detached signature, falsifiers on every Verified or
  Validated claim, and a matching content-hash dataset anchor on every
  Validated claim.
- `check_json_for_release(text, expected_root, verifier) -> CheckReport` — the
  strict-parser release entry point; malformed transports fail before
  admission.
- `CheckReport { verdict, merkle_root, breakdown, signature, findings }` —
  `passed()` and `render_pie()` (a deterministic text budget pie).
- `Verdict { Pass, Fail }`; `SignatureStatus { Unsigned, Unverified, Valid }`;
  `Finding { kind, detail }`.
- Invalid packages carry a zeroed breakdown, so a refused claim cannot retain a
  normal-looking positive evidence pie alongside the failure finding.
- Re-exports `EvidencePackage`, `ContentHash`, `ColorBreakdown`,
  `MagnitudeBudget`, `PackageError`, and `ParseError`.

## What it re-verifies

1. Format support + per-claim completeness (delegated to
   `EvidencePackage::verify` — no solver).
2. The content address: the Merkle root, recomputed independently and
   (optionally) checked against an expected value.
3. Signature validity only through an injected `SignatureVerifier` over the
   recomputed root. Presence without a verifier remains `Unverified`.
4. For explicit release admission only: non-vacuity, authenticated signature,
   per-certificate falsifier pairing, and per-Validated-claim dataset anchors.

## Invariants

- No solver / license in the build graph (enforced by the dependency set).
- A package that fails `verify` (incomplete claim, unsupported format) yields
  `Verdict::Fail` with a matching finding; a content-address mismatch fails.
- An empty package verifies vacuously and renders a "no claims" pie.
- An empty package never passes `check_for_release`; ordinary integrity
  checking and release admission are deliberately distinct verdicts.
- Verified and Validated claims never pass release admission without attached,
  structurally valid falsifier records. Validated claims additionally require
  an exact matching canonical dataset anchor.
- `render_pie` and the report are deterministic.

## Error model

The checker does not error — it REPORTS: failures become `Finding`s in a
`CheckReport` with `Verdict::Fail`. No panics.

## Determinism class

Fully deterministic: the report and rendered pie are pure functions of the
package.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/checker.rs` (Proposal 12): clean pass with no findings;
incomplete-validated-claim failure; content-address (Merkle) tamper detection;
including provenance tamper; malformed falsifier refusal with fail-closed pie;
signature-presence and verifier-capability reporting; deterministic budget-pie
rendering; empty-package pie; protocol version; determinism; and release-gate
admission/refusal batteries for empty, unpaired, unanchored, unsigned, and
wrong-root packages through both in-memory and strict JSON entry points.

## Independent re-verification (bead qmao.6.1)

`check_json` is the third-party entry point: strict parse (root
recomputation and budget re-derivation happen in the parser), then
semantic re-verification, optionally against an expected root and a
`SignatureVerifier` capability. Signature VALIDITY is asserted only
when a supplied capability accepts the signature over the RECOMPUTED
root; the in-tree `NoSignatureVerifier` accepts nothing (the no-crypto
no-claim — presence is recorded as `Unverified`, and supplying a
capability that rejects raises a `signature-invalid` finding). The
magnitude budget must reconcile with its parts. The normal dependency graph is
`fs-package -> {fs-blake3, fs-crosswalk, fs-evidence -> fs-obs}`: it contains no
solver and the checker cannot run a solve by construction.

## No-claim boundaries

- This crate ships no cryptographic primitive. It can assert signature validity
  only when a caller injects a `SignatureVerifier`; the default authenticates
  nothing.
- Composition receipts are re-run, but source-certificate production is not.
  The certificates are CARRIED in the package; the checker validates their
  structure and derivation receipts without running a solver.
- Schema v4 does not yet encode a non-forgeable source origin for a raw
  `Verified` claim. Content addressing proves package integrity, not scientific
  truth; schema-v5 ClaimOrigin work is tracked separately. Release admission
  makes falsifier/anchor/signature obligations structural, but cannot prove
  that a publicly constructible record corresponds to work that actually ran.
