# CONTRACT: fs-crosswalk

Regulatory vocabulary crosswalk (plan addendum, Proposal 12): maps
evidence-package fields onto the regulator's existing standards language.

## Purpose and layer

Layer UTIL (pure data + audit; no dependencies). It owns risk R9
(standards-body latency): every field is mapped or explicitly flagged, so the
package doubles as internal-QA / B2B diligence collateral regardless of
standards-body pace. `CROSSWALK_VERSION = 4` identifies this vocabulary and
`SUPPORTED_PACKAGE_FORMAT = 6` makes package compatibility explicit without a
dependency cycle.

## Public types and semantics

- `PackageConcept` (12) — the evidence-package fields (three colors,
  certificate, falsifier log, regime tag, anchoring dataset, provenance, Merkle
  root, signature, re-verified claim origin, and authenticated waiver
  authorization); `Standard` (4) — ASME V&V 10 / 20 / 40 and FAA/EASA CbA.
  Both expose `ALL` and `label()`; `Standard::full_name()`.
- `CrosswalkEntry { concept, standard, counterpart }`; `Counterpart` is
  `Mapped { clause, note }` or `NoCounterpart { reason }`.
- `crosswalk()` — all 48 rows. `for_concept` / `for_standard` slices;
  `lookup(concept, standard)` — the single row.
- `audit() -> CrosswalkAudit` — every `(concept, standard)` pair must have
  exactly one row (mapped or an explicit no-counterpart); `ok()` iff no gaps.
- `to_json()` — deterministic machine-readable record carrying both the
  vocabulary version and supported package format.

## Invariants

- COMPLETENESS: all 12 × 4 = 48 pairs are covered; `mapped + no_counterpart ==
  48`, no gaps.
- HONESTY: a concept with no named counterpart in a standard is flagged
  `NoCounterpart` (with a reason), never force-mapped — the crosswalk contains
  real no-counterpart rows (content-addressed integrity in the ASME standards;
  adversarial falsification logs in V&V 10/20; typed claim origins in V&V
  10/20; authenticated waivers in all ASME vocabularies).

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
per-concept (×4) and per-standard (×12) slices; representative validation,
claim-origin, and waiver decisions; unique labels; deterministic JSON.

## Vocabulary v4: schema-v6 authenticated concepts

Package format 6 makes scientific claim origin and waiver authorization
independently observable proof states. `ClaimOrigin` maps to
credibility-evidence traceability in ASME V&V 40 and analysis-data traceability
in FAA/EASA CbA; Estimated declarations, administrative waivers, and
waiver-dependent descendants do not establish this concept. V&V 10/20 have
explicit no-counterpart rows. `WaiverAuthorization` maps only to FAA/EASA's
approved deviations/limitations vocabulary. Signature rows refer specifically
to a policy-authenticated signature carrying release-purpose intent bound to
checker protocol and package root. This mapping does not establish checker
admission and does not claim signer identity, role, authorization, or authorship.

## Package-grounded coverage (bead qmao.6.1)

The static concept↔standard table is a MAPPING, never coverage.
`package_presence` judges each concept after deny-all package verification, and
`package_presence_with` accepts explicit origin/artifact capabilities.
`PackagePresenceReport` and `PackageCoverageReport` are sealed, retain the exact
package receipt, and bind rows plus mapping context into decision hashes.
Coverage is `Covered` only for the intersection of mapped and authenticated
evidence; a mapped concept with absent evidence is `MappedButAbsent`.
Claim-origin presence requires scientific admission and a successfully invoked
source/anchor/derivation policy. Waiver authorization requires an authenticated,
unexpired waiver. Signature coverage requires release-purpose authentication,
but remains a signature-intent observation rather than checker admission.
Provenance strings remain absent coverage until a provenance-artifact verifier
exists. Raw declarations never count.

## No-claim boundaries

- The mappings are a FIRST-PARTY engineering reading of the standards to frame
  the artifact — NOT an official ASME/FAA/EASA determination. The kill
  criterion (R9) is market: if no auditor engages the machine-checkable format
  even as supplementary evidence in the first regulated-vertical cycle, the
  package stays internal-QA / B2B collateral and crosswalk investment pauses.
- This crate is DATA + audit; it does not parse a package or drive the checker
  (fs-package / fs-checker do). A tool that renders a package AS a given
  standard's evidence dossier is a downstream consumer.
- A mapped waiver authorization means only that the configured policy accepted
  the bounded exception through the recorded day. It does not identify a
  responsible person and is not a certificate, validation result, or proof of
  the waived claim.
- Clause references are indicative concept names, not verbatim clause numbers
  (which move between standard editions); the crosswalk is versioned so it can
  track a specific edition.
