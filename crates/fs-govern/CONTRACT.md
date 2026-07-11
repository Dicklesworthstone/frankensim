# CONTRACT: fs-govern

The addendum's governance as machine-readable data: the design principles
(P1–P8), the governance rules, the nineteen proposals (with kill metrics +
owning beads), and the risk register (Part V, R1–R10) — each with a
CI-gateable completeness audit.

## Purpose and layer

Layer UTIL. Pure data + audit, with `fs-blake3` as its only dependency for
canonical content identities. Encodes the doctrine, the proposals, and the ten
named risks, and audits that nothing survives unmeasured (design principle P8 /
Governance Rule 2).

## Crate registry (`crates` module)

- `addendum_crates() -> &[AddendumCrate]` — the seven net-new crates the
  addendum introduced, each `{ name, purpose, owning_proposal, layer, no_claim
  }`. `crate_audit() -> CrateAudit` confirms each declares a purpose, an owner,
  and a no-claim boundary (the AGENTS.md contract discipline made
  governance-legible); `crates_json()` emits the deterministic record. Actual
  `CONTRACT.md` file presence is enforced separately by `xtask check-contracts`.

## Doctrine and proposals (`doctrine`, `proposals` modules)

- `principles() -> &[Principle]` — the eight design principles P1–P8 (id, name,
  statement); `rules() -> &[GovernanceRule]` — the four governance rules
  (number, name, statement).
- `proposals() -> &[Proposal]` — the nineteen proposals in composite (Mean)
  order, each `{ id, name, phase, mean, kill_metric, owning_bead, receipt }`.
  `governance_audit() -> GovernanceAudit` enforces that every proposal
  DECLARES a kill metric AND an owning bead (Governance Rule 2), counting how
  many are instrumented; `proposals_json()` emits the deterministic
  machine-readable record.

## Public types and semantics

- `RiskId` (`R1`..`R10`) with `RiskId::ALL` and `code()`.
- `InstrumentationReceipt::new(subject, dashboard, verifier,
  evidence_artifact, verified_day)` validates mandatory provenance and returns
  an opaque receipt. Its private fields prevent accidental identity drift;
  accessors expose the dashboard, verifier, evidence-artifact content hash,
  verification day, and receipt identity. `receipt_identity()` is the replay
  oracle.
- `Risk { id, name, description, mitigation, early_warning, threshold, owner,
  receipt }` — `early_warning` is the metric that makes the risk visible before
  it is fatal; `owner` is the bead that owns the mitigation; `receipt` is the
  optional evidence-bearing instrumentation assertion (`None` is the honest
  baseline).
- `register() -> &'static [Risk]` — the canonical R1–R10 in order;
  `risk(id) -> &'static Risk` for lookup.
- `audit(today_day) -> RiskAudit` / `audit_slice(&[Risk], today_day) ->
  RiskAudit` — checks every
  risk has a non-empty early-warning metric AND an owner, counts how many are
  instrumented, and separately lists schema gaps and operational receipt gaps.
  `declared_schema_ok()` and `operationally_managed()` deliberately expose the
  two different verdicts.
- `to_json(today_day) -> String` — a deterministic machine-readable JSON array
  with JSON-escaped strings for dashboards / CI gates. Every row carries the
  unambiguous instrumentation status and either `receipt:null` or the complete
  receipt provenance (dashboard, day, verifier, evidence artifact, identity).

## Trust boundary: declaration vs live operation (bead xpck.9)

The audits report TWO verdicts, never one: `declared_schema_ok()`
(every entry names a metric and an owner — pure schema) and
`operationally_managed()` (every metric VERIFIED live). The former
single `ok()` collapsed these and rendered the zero-instrumented
registry as green — the false-green this bead removed. Instrumentation
is EVIDENCE, not a flag: an entry counts as verified only through an
`InstrumentationReceipt` that binds the subject id, dashboard locator,
verifier identity, supporting evidence-artifact content hash, and verification
day. The canonical encoding uses tagged, `u64`-length-prefixed fields under
BLAKE3 derive-key domain
`frankensim.fs-govern.instrumentation-receipt.v1`; all receipt fields are
private. Subject replay, a future verification date, missing provenance, stale
evidence, or an inconsistent identity fails closed (`BadReceipt`/`Stale`).
Audits and JSON take `today_day` (days since 2026-01-01) explicitly, so verdicts
are deterministic and replayable.

The BLAKE3 root is an **unkeyed content identity, not an authentication tag or
signature**. It provides collision-resistant identity and accidental-tamper
detection; it does not prove that the dashboard was live, that `verifier` was
authorized, or that the evidence artifact is scientifically adequate. The
canonical registry remains code-reviewed governance data, and issuer trust /
artifact checking are deployment policy. Calling a public hash an
"authentication fingerprint" would overstate this crate's security contract.

## Invariants

- The register declares all ten risks while honestly remaining operationally
  red until receipts are installed.
- `register()` and `RiskId::ALL` share the same order.
- `to_json()` and `audit()` are deterministic.
- Receipt identities change when any semantic field changes, are bound to one
  governed subject, and cannot be mutated through the safe public API.

## Error model

`InstrumentationReceipt::new` returns `ReceiptError` for an empty subject,
dashboard, or verifier. Audits report missing, stale, future-dated,
subject-mismatched, and otherwise inconsistent receipts as data and fail
closed; they do not panic or silently promote coverage.

## Determinism class

Fully deterministic — pure functions over `const` data, no RNG or I/O.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/register.rs` (Part V, 10 cases): all ten risks present + ordered;
every risk has a metric/owner/mitigation; owners are real bead ids; lookup;
the canonical audit is complete with an honest zero-instrumented baseline;
`audit_slice` detects missing metric AND owner on an incomplete entry (the
audit is not vacuous); subject-replay/future/stale receipts fail closed;
missing provenance is rejected; changing subject, dashboard, verifier,
evidence artifact, or day changes the identity; JSON includes explicit receipt
provenance; determinism.

`tests/governance.rs` (8 cases): eight principles P1–P8; four rules numbered
1–4 (Rule 2 = kill-criteria enforcement); all nineteen proposals present with
unique ids and in descending composite order; every proposal declares a kill
metric + bead-id owner; the governance audit is complete with a zero
instrumented baseline; owner mapping spot-checks; `proposals_json` is
well-formed + deterministic.

## No-claim boundaries

- This crate encodes the risk register as governance DATA; it does not itself
  measure an early-warning metric, fetch an evidence artifact, authenticate a
  verifier, or prove dashboard liveness. A dashboard/CI supplies that evidence
  and deployment policy establishes issuer authority. The audit enforces that
  each risk declares a metric + owner and fails closed when receipt evidence is
  absent or malformed; it cannot establish the truth of an issuer's assertion.
- The addendum's design principles (P1–P8), the four governance rules, and the
  per-proposal kill criteria are the sibling doctrine bead's scope; this crate
  is the risk register (R1–R10) only.
- Bead-id owners are string references; this crate does not read the beads
  database (that coupling is deliberately avoided).
