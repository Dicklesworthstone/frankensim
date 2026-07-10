# CONTRACT: fs-govern

The addendum's governance as machine-readable data: the design principles
(P1–P8), the governance rules, the nineteen proposals (with kill metrics +
owning beads), and the risk register (Part V, R1–R10) — each with a
CI-gateable completeness audit.

## Purpose and layer

Layer UTIL. Pure data + audit — no dependencies. Encodes the doctrine, the
proposals, and the ten named risks, and audits that nothing survives
unmeasured (design principle P8 / Governance Rule 2).

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
  order, each `{ id, name, phase, mean, kill_metric, owning_bead, instrumented
  }`. `governance_audit() -> GovernanceAudit` enforces that every proposal
  DECLARES a kill metric AND an owning bead (Governance Rule 2), counting how
  many are instrumented; `proposals_json()` emits the deterministic
  machine-readable record.

## Public types and semantics

- `RiskId` (`R1`..`R10`) with `RiskId::ALL` and `code()`.
- `Risk { id, name, description, mitigation, early_warning, threshold, owner,
  instrumented }` — `early_warning` is the metric that makes the risk visible
  before it is fatal; `owner` is the bead that owns the mitigation;
  `instrumented` is whether that metric is live on a dashboard (default
  `false` — the honest baseline).
- `register() -> &'static [Risk]` — the canonical R1–R10 in order;
  `risk(id) -> &'static Risk` for lookup.
- `audit() -> RiskAudit` / `audit_slice(&[Risk]) -> RiskAudit` — checks every
  risk has a non-empty early-warning metric AND an owner, counts how many are
  instrumented, and lists `(RiskId, reason)` gaps. `RiskAudit::ok()` is true
  iff there are no gaps.
- `to_json() -> String` — a deterministic machine-readable JSON array (one
  object per risk: id, name, early_warning, threshold, owner, instrumented,
  mitigation) with JSON-escaped strings, for dashboards / CI gates.

## Trust boundary: declaration vs live operation (bead xpck.9)

The audits report TWO verdicts, never one: `declared_schema_ok()`
(every entry names a metric and an owner — pure schema) and
`operationally_managed()` (every metric VERIFIED live). The former
single `ok()` collapsed these and rendered the zero-instrumented
registry as green — the false-green this bead removed. Instrumentation
is EVIDENCE, not a flag: an entry counts as verified only through an
`InstrumentationReceipt` whose FNV fingerprint authenticates over
(subject id, dashboard, verified day), that is not dated in the
future, and that is younger than `MAX_RECEIPT_AGE_DAYS` (45) — stale
or fingerprint-inconsistent receipts DEMOTE coverage
(`Stale`/`BadReceipt` in the exact-gap lists) and can never be
outweighed by flipping a boolean. Audits and JSON take `today_day`
(days since 2026-01-01) explicitly, so verdicts are deterministic and
replayable. Migration: `Risk::instrumented`/`Proposal::instrumented`
(bool) became `receipt: Option<InstrumentationReceipt>`; JSON emits
`"instrumentation":"verified|stale|bad-receipt|uninstrumented"` in
place of the boolean; no other crate consumed the old fields
(verified at migration).

## Invariants

- The register is complete: `audit().ok()` is true and `audit().complete == 10`.
- `register()` and `RiskId::ALL` share the same order.
- `to_json()` and `audit()` are deterministic.

## Error model

None (pure data + total functions); the audit reports gaps as data, it does
not error.

## Determinism class

Fully deterministic — pure functions over `const` data, no RNG or I/O.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/register.rs` (Part V, 8 cases): all ten risks present + ordered;
every risk has a metric/owner/mitigation; owners are real bead ids; lookup;
the canonical audit is complete with an honest zero-instrumented baseline;
`audit_slice` detects missing metric AND owner on an incomplete entry (the
audit is not vacuous); JSON is well-formed + complete; determinism.

`tests/governance.rs` (8 cases): eight principles P1–P8; four rules numbered
1–4 (Rule 2 = kill-criteria enforcement); all nineteen proposals present with
unique ids and in descending composite order; every proposal declares a kill
metric + bead-id owner; the governance audit is complete with a zero
instrumented baseline; owner mapping spot-checks; `proposals_json` is
well-formed + deterministic.

## No-claim boundaries

- This crate encodes the risk register as governance DATA; it does not itself
  measure the early-warning metrics — a dashboard/CI wires them and flips
  `instrumented`. The audit enforces that each risk DECLARES a metric + owner,
  not that the metric is currently green.
- The addendum's design principles (P1–P8), the four governance rules, and the
  per-proposal kill criteria are the sibling doctrine bead's scope; this crate
  is the risk register (R1–R10) only.
- Bead-id owners are string references; this crate does not read the beads
  database (that coupling is deliberately avoided).
