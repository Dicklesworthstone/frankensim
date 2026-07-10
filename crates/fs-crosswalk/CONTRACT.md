# CONTRACT: fs-crosswalk

Regulatory vocabulary crosswalk (plan addendum, Proposal 12): maps
evidence-package fields onto the regulator's existing standards language.

## Purpose and layer

Layer UTIL (pure data + audit; no dependencies). It owns risk R9
(standards-body latency): every field is mapped or explicitly flagged, so the
package doubles as internal-QA / B2B diligence collateral regardless of
standards-body pace. `CROSSWALK_VERSION = 2` identifies this vocabulary and
`SUPPORTED_PACKAGE_FORMAT = 4` makes package compatibility explicit without a
dependency cycle.

## Public types and semantics

- `PackageConcept` (10) — the evidence-package fields (three colors,
  certificate, falsifier log, regime tag, anchoring dataset, provenance, Merkle
  root, signature); `Standard` (4) — ASME V&V 10 / 20 / 40 and FAA/EASA CbA.
  Both expose `ALL` and `label()`; `Standard::full_name()`.
- `CrosswalkEntry { concept, standard, counterpart }`; `Counterpart` is
  `Mapped { clause, note }` or `NoCounterpart { reason }`.
- `crosswalk()` — all 40 rows. `for_concept` / `for_standard` slices;
  `lookup(concept, standard)` — the single row.
- `audit() -> CrosswalkAudit` — every `(concept, standard)` pair must have
  exactly one row (mapped or an explicit no-counterpart); `ok()` iff no gaps.
- `to_json()` — deterministic machine-readable record carrying both the
  vocabulary version and supported package format.

## Invariants

- COMPLETENESS: all 10 × 4 = 40 pairs are covered; `mapped + no_counterpart ==
  40`, no gaps.
- HONESTY: a concept with no named counterpart in a standard is flagged
  `NoCounterpart` (with a reason), never force-mapped — the crosswalk contains
  real no-counterpart rows (content-addressed integrity in the ASME standards;
  adversarial falsification logs in V&V 10/20).

## Error model

Total functions; no panics.

## Determinism class

Fully deterministic: a pure `const` table; `to_json` reproduces byte-for-byte.

## Cancellation behavior

None.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/crosswalk.rs` (Proposal 12, 8 cases): explicit package/vocabulary
compatibility versions; full coverage + no silent gaps; honesty (explicit
no-counterpart rows exist, verified-color does map);
per-concept (×4) and per-standard (×10) slices; a representative
validated→validation-metric lookup; unique labels; deterministic JSON.

## Package-grounded coverage (bead qmao.6.1)

The static concept↔standard table is a MAPPING, never coverage.
`package_presence` judges each concept against fields actually present
in a parsed package, and `package_coverage` reports Covered only for
the intersection of (mapped, evidenced); a mapped concept with absent
evidence is `MappedButAbsent`, and concepts unrepresentable in schema
v2 (falsifier logs) can never read Covered (tested across all four
standards).

## No-claim boundaries

- The mappings are a FIRST-PARTY engineering reading of the standards to frame
  the artifact — NOT an official ASME/FAA/EASA determination. The kill
  criterion (R9) is market: if no auditor engages the machine-checkable format
  even as supplementary evidence in the first regulated-vertical cycle, the
  package stays internal-QA / B2B collateral and crosswalk investment pauses.
- This crate is DATA + audit; it does not parse a package or drive the checker
  (fs-package / fs-checker do). A tool that renders a package AS a given
  standard's evidence dossier is a downstream consumer.
- Clause references are indicative concept names, not verbatim clause numbers
  (which move between standard editions); the crosswalk is versioned so it can
  track a specific edition.
