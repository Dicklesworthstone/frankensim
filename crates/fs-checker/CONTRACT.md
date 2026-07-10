# CONTRACT: fs-checker

The standalone evidence-package checker (plan addendum, Proposal 12): an
independently distributable verifier — "don't trust us; here is the checker."

## Purpose and layer

Layer L6. Its ENTIRE dependency graph is `fs-package` → `fs-evidence`. A HARD
distribution constraint (Proposal 12): NO solver stack, NO geometry kernel, NO
license gate anywhere in the graph — by construction the checker cannot run a
solve. It carries its own `CHECKER_PROTOCOL_VERSION` (distributed
independently).

## Public types and semantics

- `check(&EvidencePackage) -> CheckReport` — re-verify a package.
- `check_against_root(&EvidencePackage, expected_root) -> CheckReport` — also
  confirm the content address matches (tamper / substitution detection).
- `CheckReport { verdict, merkle_root, breakdown, signature, findings }` —
  `passed()` and `render_pie()` (a deterministic text budget pie).
- `Verdict { Pass, Fail }`; `SignatureStatus { Unsigned, Present(sig) }`;
  `Finding { kind, detail }`.
- Re-exports `EvidencePackage`, `ColorBreakdown`, `PackageError`.

## What it re-verifies

1. Format support + per-claim completeness (delegated to
   `EvidencePackage::verify` — no solver).
2. The content address: the Merkle root, recomputed independently and
   (optionally) checked against an expected value.
3. Signature PRESENCE (recorded; cryptographic assertion awaits a Franken
   signature primitive — a present signature is not silently treated as valid).

## Invariants

- No solver / license in the build graph (enforced by the dependency set).
- A package that fails `verify` (incomplete claim, unsupported format) yields
  `Verdict::Fail` with a matching finding; a content-address mismatch fails.
- An empty package verifies vacuously and renders a "no claims" pie.
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

`tests/checker.rs` (Proposal 12, 9 cases): clean pass with no findings;
incomplete-validated-claim failure; content-address (Merkle) tamper detection;
including provenance tamper; signature-presence reporting; deterministic
budget-pie rendering; empty-package pie; protocol version; determinism.

## Independent re-verification (bead qmao.6.1)

`check_json` is the third-party entry point: strict parse (root
recomputation and budget re-derivation happen in the parser), then
semantic re-verification, optionally against an expected root and a
`SignatureVerifier` capability. Signature VALIDITY is asserted only
when a supplied capability accepts the signature over the RECOMPUTED
root; the in-tree `NoSignatureVerifier` accepts nothing (the no-crypto
no-claim — presence is recorded as `Unverified`, and supplying a
capability that rejects raises a `signature-invalid` finding). The
magnitude budget must reconcile with its parts. The dependency graph
remains fs-package → fs-evidence only: the checker cannot run a solve
by construction.

## No-claim boundaries

- Cryptographic SIGNATURE verification is not performed (no Franken-compliant
  primitive yet); the bundle is trusted by CONTENT ADDRESS and a present
  signature is recorded, not asserted valid. Wiring a primitive is later work.
- v1 re-verifies flat claims: completeness + content address. Re-running
  CONTRACT COMPOSITION (when packages carry composition DAGs, cf. fs-contract)
  and per-interval arithmetic re-derivation are follow-ons.
- The certificates are CARRIED in the package; the checker re-checks their
  structural validity, it does not re-derive them (that is the whole point —
  no solver).
