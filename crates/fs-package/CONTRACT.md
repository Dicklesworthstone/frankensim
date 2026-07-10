# CONTRACT: fs-package

Machine-checkable evidence packages (plan addendum, Proposal 12): a
content-addressed bundle of color-typed claims a standalone checker can
re-verify without solvers.

## Purpose and layer

Layer L6. Depends only on `fs-evidence` (UTIL — `Color`, `ColorRank`,
`ValidityDomain`). Pure, deterministic; no I/O.

## Public types and semantics

- `Claim { id, statement, color }` — a claim plus its epistemic color (which
  carries the certificate payload).
- `Provenance { code_version, constellation_lock }`.
- `EvidencePackage { format_version, claims, provenance, signature }` —
  builder: `new(prov).with_claim(..).signed(..)`.
  - `merkle_root() -> u64` — an FNV-1a Merkle root over the package identity:
    format version, provenance, and ordered claims. Detached signatures are
    excluded. Any reproducibility provenance or claim change changes it.
  - `verify() -> Result<PackageReport, PackageError>` — re-verify WITHOUT a
    solver: the format must be `FORMAT_VERSION`, and every claim must be
    complete for its color.
  - `color_breakdown() -> ColorBreakdown` — the by-color budget pie.
  - `to_json()` — deterministic self-describing JSON (carries the root hex).
- `PackageReport { merkle_root, breakdown, claims }`.
- `PackageError` — `IncompleteValidatedClaim { claim, missing }` /
  `IncompleteVerifiedClaim { claim }` / `UnsupportedFormat { found }`.

## Invariants

- COMPLETENESS: a `Validated` claim must have a non-empty regime (`regime.
  bounds()` non-empty) AND a non-blank anchoring `dataset`; a `Verified` claim
  must carry a finite `[lo <= hi]` interval. An `Estimated` claim needs no
  certificate — an all-estimated package is valid and round-trips.
- CONTENT-ADDRESSING: `merkle_root` is deterministic and tamper-evident across
  format version, provenance, and claims; a detached signature does not change it.
- `verify` runs no solver — pure structural re-verification (the checker's
  core).

## Error model

Structured `PackageError` values (refusals that teach), never panics.

## Determinism class

Fully deterministic: the Merkle root and JSON are pure functions of the
package (bit-exact on float certificate payloads via `to_bits`).

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/package.rs` (Proposal 12, 10 cases): complete mixed-color package;
all-estimated boundary (valid + round-trips); validated-missing-regime and
validated-missing-dataset completeness failures; verified bad-interval
failure; Merkle determinism + claim/provenance tamper detection;
unsupported-format rejection; optional detached signature; deterministic JSON
carrying the root.

## Schema v3: receipts, falsifiers, anchors (bead xfxq)

Claims optionally carry a COMPOSITION RECEIPT (parent claim indices +
the ledger op): `verify` re-runs `fs_evidence::compose` over the
parents in order and requires EXACT color equality — a claim whose
color cannot be re-derived is `ReceiptMismatch`, so laundering a
Verified from Estimated parents is caught SEMANTICALLY, not just by
the content address; parents must point strictly backwards
(`BadReceiptParent` otherwise — a DAG by construction). FALSIFIER
RECORDS (name, attempts, refuted, detail) travel with the claim; any
`refuted: true` fails verification outright (`RefutedClaim`). ANCHOR
RECORDS give validated claims content-hashed dataset identities. All
three field families bind into the content address and round-trip
through the strict parser (booleans added to the closed grammar);
crosswalk coverage now reads falsifier logs and anchors from the
actual fields (present only when records exist). The checker stays
solver-free: `compose` lives in fs-evidence, already in its dependency
cone. Migration: v2 readers refuse v3 by version (the one-version
contract); in-tree constructors gained fields with builders, no
call-site changes.

## Schema v2: round trip and fail-closed parsing (bead qmao.6.1)

`to_json` emits the COMPLETE color payloads (floats as bit-exact hex
bits), provenance, signature, magnitude budget, and the content root;
`from_json` is a STRICT parser — unknown fields, missing fields, wrong
types, duplicate keys, bad hex, non-finite certificates, inverted
intervals, negative dispersions, and unknown color kinds each refuse
with a structured `ParseError`. The parser re-derives the magnitude
budget from the parsed claims and RECOMPUTES the content root from the
parsed fields: a package whose embedded root does not recompute is
tampered or forged and never loads — which is what makes a
structurally-shaped forged Verified claim fail. Decode-encode is both
semantically and textually stable (tested). The magnitude budget
attributes ERROR MAGNITUDES (verified interval widths, estimated
dispersions) and counts validated claims as unquantified regional
trust — never numerified.

## No-claim boundaries

- The Merkle hash is an in-house FNV-1a (Franken-compliant, pure Rust); a
  production build swaps in fs-ledger's BLAKE3-class hash. A cryptographic
  SIGNATURE is detached and OPTIONAL — the bundle is verifiable by content
  address regardless; wiring a Franken signature primitive is later work.
- `verify` checks STRUCTURAL completeness + the content address; it does not
  re-run solvers to re-derive the certificates (that is the point — the
  certificates are carried). The standalone distributable checker (a separate
  bead) wraps this crate.
- The certificate payloads live in `fs-evidence::Color`; this crate bundles
  and content-addresses them, it does not produce them.
