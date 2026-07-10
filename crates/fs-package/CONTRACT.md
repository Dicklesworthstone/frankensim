# CONTRACT: fs-package

Machine-checkable evidence packages (plan addendum, Proposal 12): a
content-addressed bundle of color-typed claims a standalone checker can
re-verify without solvers.

## Purpose and layer

Layer L6. Depends on `fs-evidence` (UTIL — `Color`, `ColorRank`,
`ValidityDomain`), `fs-crosswalk` (the static standards vocabulary used by
coverage reports), and dependency-free `fs-blake3` (the shared content-hash
owner). Pure, deterministic; no I/O and no solver dependency.

## Public types and semantics

- `Claim { id, statement, color }` — a claim plus its epistemic color (which
  carries the certificate payload). Claim ids are non-blank and unique within
  a package.
- `Provenance { code_version, constellation_lock }`.
- `EvidencePackage { format_version, claims, provenance, signature }` —
  builder: `new(prov).with_claim(..).signed(..)`.
  - `merkle_root() -> ContentHash` — a 32-byte BLAKE3 Merkle root over the
    package identity: format version, claim count, provenance, and ordered
    claims. Header, claim, and internal-node hashes use standard BLAKE3
    derive-key domains. Detached signatures are excluded. Any reproducibility
    provenance or claim change changes the root.
  - `verify() -> Result<PackageReport, PackageError>` — re-verify WITHOUT a
    solver: the format must be `FORMAT_VERSION`, and every claim must be
    complete for its color.
  - `color_breakdown() -> ColorBreakdown` — the by-color budget pie.
  - `to_json()` — deterministic self-describing JSON (carries the root hex).
- `PackageReport { merkle_root, breakdown, claims }`.
- `PackageError` — structured refusals for incomplete provenance, invalid or
  duplicate claim ids, malformed color payloads, unsupported formats, receipt
  mismatches/parents, malformed falsifier/anchor records, and refuted claims.

## Invariants

- COMPLETENESS: reproducibility provenance fields and claim ids are non-blank;
  claim ids are unique. A `Validated` claim must have a non-empty regime
  (`regime.bounds()` non-empty) whose axis names are non-blank and whose bounds
  are finite and ordered, plus a non-blank anchoring `dataset`. A `Verified`
  claim must carry a finite `[lo <= hi]` interval. An `Estimated` claim needs a
  non-blank estimator identity and a non-negative, non-NaN dispersion.
  Positive infinity is preserved as the lower-layer algebra's explicit
  no-quantitative-spread-claim sentinel; it is distinct from finite subtotal
  overflow, which verification rejects. An honest all-estimated package
  remains valid.
- FALSIFIER EVIDENCE: a record has a non-blank stable falsifier identity, at
  least one executed attempt, and a non-blank outcome detail. A refuted record
  still rejects its claim and package.
- DATASET ANCHORS: every attached anchor has a non-blank stable dataset id and
  an exactly 64-character, lowercase hexadecimal content hash. Crosswalk
  anchoring coverage requires a valid anchor whose dataset id exactly matches
  the `Validated` claim's named dataset; unrelated anchors do not count.
- CONTENT-ADDRESSING: `merkle_root` is deterministic and tamper-evident across
  format version, claim count, provenance, and ordered claims. Derive-key modes
  separate package leaves/nodes from plain artifact hashes and each other; a
  detached signature does not change the root.
- `verify` runs no solver — pure structural re-verification (the checker's
  core).

## Error model

Structured `PackageError` values (refusals that teach), never panics. The JSON
parser maps the same package-level semantic refusals into `ParseError`, so a
package cannot pass one entry point and fail the other. The untrusted JSON
boundary is bounded before schema mapping: 64 MiB input, depth 64, one million
values, 100,000 members per container, 1 MiB decoded strings, and 128-byte
number tokens. In-memory verification enforces the corresponding transport
envelope before a package can pass, so a verified object remains serializable
and checkable under those bounds. Limit violations are structured refusals.

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

`tests/package.rs` (Proposal 12): complete mixed-color package;
all-estimated boundary (valid + round-trips); validated-missing-regime and
validated-missing-dataset completeness failures; verified bad-interval
failure; Merkle determinism + claim/provenance tamper detection;
unsupported-format rejection; optional detached signature; deterministic JSON
carrying the root; in-memory/serialized semantic parity; and exact full-width
falsifier attempt-count round trips with overflow refusal.

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
crosswalk coverage now reads validated falsifier logs and matching anchors from
the actual fields. Coverage fails closed for an invalid package, and raw
detached-signature presence never counts as authenticated sign-off. The
checker stays solver-free: `compose` lives in fs-evidence, already in its
dependency cone. Migration: v2 readers refuse v3 by version (the one-version
contract); in-tree constructors gained fields with builders, no call-site
changes.

## Schema v2: round trip and fail-closed parsing (bead qmao.6.1)

`to_json` emits the COMPLETE color payloads (floats as bit-exact hex
bits), provenance, signature, magnitude budget, and the content root;
`from_json` is a STRICT parser — unknown fields, missing fields, wrong
types, duplicate keys, bad hex, non-finite certificates, inverted
intervals, negative dispersions, and unknown color kinds each refuse
with a structured `ParseError`. The parser re-derives the magnitude
budget from the parsed claims and RECOMPUTES the content root from the
parsed fields: a package whose embedded root does not recompute is
tampered in transit and never loads. This is an integrity check, not the
schema-v5 source-origin proof named below. Decode-encode is both
semantically and textually stable (tested). The magnitude budget
attributes ERROR MAGNITUDES (verified interval widths, estimated
dispersions) and counts validated claims as unquantified regional
trust — never numerified. JSON number tokens are retained as decimal text until
they are converted into their target integer type; full-width `u64` falsifier
attempt counts therefore round-trip exactly, and overflow refuses instead of
rounding through `f64`.

## No-claim boundaries

- The Merkle hash is the shared in-house BLAKE3 implementation, with standard
  derive-key mode separation for typed package roots. This establishes content
  integrity, not authorship or scientific truth. A cryptographic SIGNATURE is
  detached and OPTIONAL; wiring a default Franken signature primitive is later
  work.
- `verify` checks STRUCTURAL completeness + the content address; it does not
  re-run solvers to re-derive the certificates (that is the point — the
  certificates are carried). The standalone distributable checker (a separate
  bead) wraps this crate.
- The certificate payloads live in `fs-evidence::Color`; this crate bundles
  and content-addresses them, it does not produce them.
- Schema v4's public `Claim` and `Color` fields do not prove how a source claim
  was obtained. A party can construct a fresh package around a structurally
  valid raw `Verified` interval and recompute its root. The root establishes
  integrity, not epistemic origin; sealed ClaimOrigin transport is tracked for
  schema v5.
