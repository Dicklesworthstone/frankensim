# CONTRACT: fs-govern

The addendum's governance as machine-readable data: the design principles
(P1–P8), the governance rules, the nineteen proposals (with kill metrics +
owning beads), the original Part V runtime risks (R1–R10), the distinct
expansion-program risks (PR-001–PR-012) — each with a CI-gateable completeness
surface — plus the EXECUTABLE one-bet-per-lane admission state machine
(`lanes` module, bead rjoq.6) and the generated requirement-to-evidence
registry (`traceability` module, bead
`frankensim-ext-traceability-ledger-ziso`).

## Purpose and layer

Layer UTIL. Pure data + audit, with `fs-blake3` as its only dependency for
canonical content identities. Encodes the doctrine, proposals, original ten
runtime risks, and twelve expansion-program risks, audits that nothing
survives unmeasured (design principle P8 / Governance Rule 2), and enforces
the one-active-unproven-mechanism-per-independently-falsifiable-proof-lane
rule as an atomic, replayable admission ledger. It also owns the pure typed
source registry and fail-closed generator for the extension charter's B1–B14
and RQ-* requirement-to-evidence catalog.

## Proof-lane admission (`lanes` module, bead rjoq.6)

- `LaneCharter::new(statement, admissible_domain, assumptions,
  target_authority, baseline, falsifier_family, independence_class)`
  canonicalizes (whitespace collapse; assumptions sorted + deduped; every
  field non-empty and bounded before canonicalization allocates) and `lane_id()` mints the validated
  `ProofLaneId` (tagged, length-prefixed fields under BLAKE3 domain
  `frankensim.fs-govern.proof-lane.v1`). There is NO public id constructor
  from a raw hash — an id always corresponds to a validated charter
  (anti-spoofing), and cosmetic re-spellings collapse to one lane (the
  split-gaming canonicalization). `mechanism_id(name, version)` mints a
  `MechanismId` that retains both its originating `ProofLaneId` and its
  content identity; admission, comparisons, and supersession refuse a
  mechanism presented under another lane.
- `PortfolioLedger::new(PortfolioPolicy { global, max_active_mechanisms })`
  is the admission state machine (`LANE_POLICY_VERSION = 2`). Semantics:
  multiple active unproven mechanisms across independently falsifiable
  lanes; a second active mechanism in the SAME lane refuses atomically;
  lanes DECLARING the same independence class share one bet (the ledger
  retains the complete active set, so finalizing one comparison candidate
  cannot hide a surviving candidate); a comparison on another lane cannot
  evade that backstop; the four-axis global
  envelope (work / memory / reviewer slots / falsification capacity) and
  the mechanism cap bind across all lanes.
- `HeadToHeadCharter::new(&lane_charter, candidates, shared, preregistration)` is
  the ONLY carve-out: a preregistered comparison (declared BEFORE any
  admission in the lane, 2..=8 distinct candidates, non-zero
  preregistration artifact) admits exactly its declared candidates under
  its own bounded shared envelope.
- `FinalizationReceipt::new(mechanism, kind, superseded_by, ledger_artifact)`
  (kinds: refuted / tombstoned / withdrawn / superseded-with-successor) is
  the ONLY release path: `finalize` releases the slot and its reservation
  EXACTLY ONCE against an identity-consistent receipt; Unknown or stalled
  work never releases (there is deliberately no timeout path); terminal is
  permanent (no re-admission).
- Every request carries an `IdempotencyKey`: a byte-identical retry
  replays the recorded decision without double-charging or re-recording; a
  different request under a used key records one conflict row naming the
  original sequence; an exact retry of that conflicting request replays the
  same refusal without another row. Every method validates completely BEFORE
  governed state mutates; refusals may append their explicit audit decision but
  never partially charge lane/resource state. `PortfolioLedger` is deliberately
  non-`Clone`; exclusive `&mut self` access to that authority value is the
  in-process concurrency contract.
- The retained record is fail-closed and hard-bounded by
  `MAX_RETAINED_DECISIONS` and `MAX_RETAINED_DECISION_BYTES`; the primary and
  conflict-idempotency maps share the decision cap. Records are never silently
  evicted because eviction would permit an old retry to execute again. Instead,
  `RetentionCapacityExceeded` refuses before cloning variable-size request data.
  One row/key and a conservative byte allowance remain reserved for every
  active mechanism's future finalization.
- Every `AdmissionDecision` retains the complete policy and canonical replay
  preimage: statement/domain/assumptions/authority/baseline/falsifier/class,
  comparison candidates/artifact/shared envelope or admission reservation,
  terminal kind/successor/artifact/receipt identity/released envelope, plus
  lane/mechanism ids, idempotency key, request digest, verdict, and ranked
  remedy. `decisions_json(limit)` adds explicit skipped and retained-cap
  metadata; the stored request, rather than an opaque digest alone, is
  sufficient for deterministic replay.

## Crate registry (`crates` module)

- `addendum_crates() -> &[AddendumCrate]` — the seven net-new crates the
  addendum introduced, each `{ name, purpose, owning_proposal, layer, no_claim
  }`. `crate_audit() -> CrateAudit` confirms each declares a purpose, an owner,
  and a no-claim boundary (the AGENTS.md contract discipline made
  governance-legible); `crates_json()` emits the deterministic record. Actual
  `CONTRACT.md` file presence is enforced separately by `xtask check-contracts`.

## Requirement-to-evidence traceability (`traceability` module)

- `requirements()` is the closed canonical B1–B14 and sixteen-RQ source
  registry derived from the extension charter seed. Every row carries the
  exact generated-ledger columns: stable requirement id, capability/property,
  concrete blocker, owner/artifact, prerequisite phase, milestone, forcing
  flagship, benchmark/data evidence route, proof-obligation links, honest
  claim boundary, and declarative status. Status is governance declaration,
  not scientific proof and not a mirror of whether an owning Bead is closed;
  tracker `closed` is deliberately not an accepted scientific lifecycle value.
- `proof_obligations()` is the complete ordered PO-1 through PO-25 index. Each
  entry links its executable summary to the full owner Bead ids. A row may
  cite only entries in this closed index.
- `audit_traceability(rows, obligations)` accumulates structured
  `TraceabilityDiagnostic { requirement_id, field, reason }` refusals. An empty
  scope, any missing or oversized field, missing owner, missing milestone gate,
  missing benchmark/data evidence route, missing claim boundary, missing or
  dangling PO link, duplicate requirement, noncanonical PO alias, incomplete
  PO index, or duplicate PO owner is non-green. Diagnostics are sorted by
  requirement, field, and reason, so source enumeration order cannot alter the
  refusal artifact.
- `generate_traceability_ledger(rows, obligations)` validates first and emits
  no partial output on failure. Successful rows use the closed charter order,
  PO links numerically, PO definitions numerically, and owner ids lexically;
  `traceability_ledger_json()` applies it to the canonical code registry under
  schema `frankensim-requirement-traceability-v1`. Every artifact explicitly
  emits `authority: "declaration-only"` and `source_snapshot: null`; this API
  refuses `verified` and `validated` because only a future source-bound
  admission receipt may mint scientific promotion status.
- `TraceabilitySourceSnapshot::new(sources)` is the pure adapter boundary for
  live Beads/contracts/registry extraction. It requires at least one immutable
  artifact from each of those three classes, rejects blank/oversized or
  duplicate locators and the all-zero missing-content sentinel, canonicalizes
  source order, and derives a domain-separated BLAKE3 snapshot identity from
  source kind, locator, and exact content identity. It reads no filesystem and
  trusts no tracker status implicitly; the caller owns bounded parsing of the
  named bytes.
- `load_traceability_source_snapshot(root, specs, limits)` is the concrete,
  synchronous filesystem boundary. It accepts only strict root-relative
  locators, canonicalizes every regular file beneath the configured root,
  refuses final symlinks, bounds both individual and aggregate reads, rejects
  metadata/read length races, and hashes the exact admitted bytes. Its receipt
  records the adapter version, effective limits, canonical locator/class,
  byte count, exact BLAKE3 identity, Beads record count where applicable, and
  the resulting source-snapshot identity. Contract and registry sources must
  be nonempty NUL-free UTF-8. Beads sources additionally receive a bounded
  lexical JSON-object-envelope audit and a unique canonical top-level string
  `id` audit before snapshot admission.
- `generate_traceability_ledger_from_snapshot(rows, obligations, snapshot)`
  retains `authority: "declaration-only"` while adding three explicit roots:
  the admitted source-snapshot identity, the canonical unbound declaration
  identity, and a binding identity over that pair. A source-byte mutation moves
  the snapshot/binding roots; a row/PO mutation moves the
  declaration/binding roots. None of these roots proves adapter correctness,
  scientific adequacy, authorized review, or that Bead closure discharges a
  proof obligation.
- Inputs are hard-bounded before rendering: at most 256 requirements, exactly
  the closed 25-PO definition space, 25 PO links per requirement, 16 owners per
  PO, 16 KiB per scalar field, 512 source artifacts, and 4 KiB per source
  locator. The core `traceability` module remains pure. `traceability_fs` owns
  bounded standard-filesystem reads only; it does not open the Beads database,
  perform network I/O, or persist immutable artifacts. Higher tooling owns the
  semantic join and persistence and must use these validators instead of
  hand-maintaining a dashboard.

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

## Expansion-program register (`program_risks` module)

- `ProgramRiskId` is the namespace `PR-001` through `PR-012`; it is deliberately
  disjoint from the original `RiskId::R1` through `R10` runtime register.
- Every `ProgramRisk` row carries a named workstream and owning Bead, likelihood
  and impact, a leading indicator, a numeric comparator/threshold/unit/typed
  domain/minimum-sample trigger, mitigation, contingency (kill/refuse/escalate),
  residual likelihood and impact, and one E0–E7 review gate.
- `program_risks()` returns the twelve rows in canonical id order;
  `program_risk(id)` performs total typed lookup; and
  `program_risk_register_json()` emits the deterministic ledger artifact.
- `assess_program_risks(observations)` always emits twelve ordered rows. Missing,
  duplicate, non-finite, unit-mismatched, out-of-domain, and under-sampled
  evidence is non-green. Only one finite, correctly unit-tagged, in-domain,
  sufficiently sampled observation below its declared trip condition is
  `Clear`. Caller units are retained only as a 64-byte UTF-8-safe preview plus
  their exact original byte length.

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
  two different verdicts. Empty audit scopes and duplicate risk ids fail closed;
  success requires the declared/verified counts to equal the nonzero total.
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
- The program register has exactly twelve unique rows in `ProgramRiskId::ALL`
  order, and every trigger contains a finite number, explicit unit, typed
  numeric domain, and positive sampling floor.
- Program assessment is input-order independent and cannot become all-clear by
  omission, duplication, malformed units, invalid numeric domains, or NaN.
- Every mechanism is admitted only under the lane that minted it; every
  surviving member of a comparison independently keeps its declared
  independence class occupied.
- Retained decisions, primary/conflict idempotency bindings, and their
  variable-size canonical payloads remain within fixed caps. Capacity reserved
  for active terminal transitions is unavailable to unrelated traffic.
- The traceability registry contains exactly the thirty charter requirement
  ids and the complete canonical PO-1..PO-25 index. Generation is deterministic
  and cannot emit a partial ledger when any owner, gate, evidence route, proof
  link, boundary, or status is absent.

## Error model

`InstrumentationReceipt::new` returns `ReceiptError` for an empty subject,
dashboard, verifier, or an all-zero missing-evidence artifact sentinel. Audits
report missing, stale, future-dated, subject-mismatched, and otherwise
inconsistent receipts as data and fail closed; they do not panic or silently
promote coverage.

The expansion-program evaluator is total over its bounded slice input: evidence
defects are represented by `AssessmentStatus`, not panics. The caller remains
responsible for bounding the number of supplied aggregate observations; output
and retained unit text remain bounded by the fixed twelve-row register and unit
preview limit.

Proof-lane APIs return `LaneError` for empty/oversized preflighted inputs,
lane/mechanism mismatches, occupied lanes/classes, envelope/cap failures,
invalid terminal receipts, idempotency conflicts, and exhausted retained-log
capacity. Capacity exhaustion never evicts replay authority and never mutates
governed admission state.

Traceability generation returns a `TraceabilityAudit` rather than output when
the source is empty, oversized, duplicated, incomplete, or contains a dangling
reference. Diagnostics name the exact requirement and field; no missing-field
case is collapsed into an unstructured string or silently defaulted.

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
`audit_slice` detects missing metric AND owner on an incomplete entry, rejects
empty scopes and duplicate ids (the audit is not vacuous), and fails closed on
subject-replay/future/stale receipts;
missing provenance is rejected; changing subject, dashboard, verifier,
evidence artifact, or day changes the identity; JSON includes explicit receipt
provenance; determinism.

`tests/governance.rs` (8 cases): eight principles P1–P8; four rules numbered
1–4 (Rule 2 = kill-criteria enforcement); all nineteen proposals present with
unique ids and in descending composite order; every proposal declares a kill
metric + bead-id owner; the governance audit is complete with a zero
instrumented baseline; owner mapping spot-checks; `proposals_json` is
well-formed + deterministic.

`tests/program_risks.rs` (G0): exactly twelve ordered unique rows; an exact
twelve-row owner/trigger/review-gate mapping lock; complete quantitative
schemas and review gates; exact comparator boundaries; fail-closed
missing/duplicate/non-finite/unit/domain/sample handling; integer/fraction
domain checks; input-order-independent assessment JSON; and deterministic
numeric-threshold register serialization.

`tests/traceability.rs` (G0): exact thirty-requirement and PO-1..PO-25 closed
inventories; canonical completeness and link coverage; one named diagnostic
for every required row field; deliberately orphaned owner detail; duplicate
requirement and dangling/repeated PO refusals; missing/noncanonical/duplicate-
owner PO index failures; success and complete failure-audit permutation
invariance; empty-scope and exact/plus-one row, PO-link, owner, scalar, id,
summary, and reference caps; declaration-only authority enforcement; and exact
generated JSON column/index coverage.

`tests/traceability_fs.rs` (G0/G3): exact byte/hash receipts, canonical source
order, source-mutation identity sensitivity, declaration-only ledger binding,
bounded Beads record/line/nesting envelope and unique-id refusals, individual
and aggregate byte caps, strict relative paths, regular-file metadata, complete
source-class coverage, limit validation, and duplicate-locator refusal.

`tests/lanes_e2e.rs` (bead rjoq.6 slice 2): the cross-crate no-mock
composition — fs-govern admission persisted into a FrankenSQLite-backed
fs-ledger (events + content-addressed preregistration/refutation/
decision-log artifacts), fs-package claims re-checked solver-free by
fs-checker incl. a mismatched-root refusal, fresh-ledger byte-for-byte
replay, and an idempotent full-retry pass.

`tests/lanes.rs` (bead rjoq.6, cases lane-001..lane-009): G0 identity and
state-machine laws (canonical collapse, per-field mutation sensitivity,
lane-bound mechanism/comparison/successor authority, same-lane refusal with
unchanged state, terminal exactly-once release, receipt validation incl. zero
evidence and self-supersession); G3
split/merge adversaries (independence-class backstop, comparison-evasion
refusal, undeclared-candidate refusal, surviving-candidate class retention);
G4 crash/retry idempotency (replay without double-charge, one-row stable
key-conflict refusal, refusal replay); G0 global mechanism-cap and independent
work/memory/reviewer/falsification-capacity boundaries for both global and
comparison envelopes; G5 whole-ledger and acceptance-complete JSON replay with
explicit truncation/cap metadata; bounded-retention refusal with a terminal
slot held for active work. Each same-lane, global-cap, terminal-release, and
identity guard has a test that fails if the guard is removed.

## No-claim boundaries

- This crate encodes the risk register as governance DATA; it does not itself
  measure an early-warning metric, fetch an evidence artifact, authenticate a
  verifier, or prove dashboard liveness. A dashboard/CI supplies that evidence
  and deployment policy establishes issuer authority. The audit enforces that
  each risk declares a metric + owner and fails closed when receipt evidence is
  absent or malformed; it cannot establish the truth of an issuer's assertion.
- The original R1–R10 register and the PR-001–PR-012 expansion-program register
  are separate authorities. Neither namespace silently substitutes for the
  other, and neither measures its own leading indicators.
- Program-risk thresholds are governance trip points over caller-supplied
  aggregates. They do not establish scientific validity, automatically execute
  a contingency, or prove that the named Bead owner has reviewed the result.
- Bead-id owners in governance declarations remain string references. The
  filesystem adapter may read exact caller-selected Beads JSONL bytes, but it
  validates only a conservative lexical object envelope and unique canonical
  top-level ids. It does not open the Beads database or interpret status,
  dependencies, closure, ownership, or scientific evidence semantics.
- The canonical traceability rows are governance declarations. They do not
  attest that an owner exists in the current tracker snapshot, that a benchmark
  ran, that a milestone closed, or that a scientific claim is verified. A
  filesystem receipt binds exact Beads/contracts/registry bytes but does not
  perform that semantic join. Higher tooling must join the admitted sources and
  V&V evidence before persisting an immutable generated artifact; a closed Bead
  must never be translated directly into scientific proof status.
- Filesystem loading uses bounded synchronous `std` I/O and has no `Cx`; it
  makes no cancellation-latency, hostile-directory race resistance, sandbox,
  signature, authorization, or durable-persistence claim. Deployments that
  admit adversarially mutable roots need a stronger capability-scoped adapter.
- The lane ledger is the admission STATE MACHINE, not durable storage: the
  `FinalizationReceipt`'s `ledger_artifact` and the comparison's
  preregistration artifact are content references whose durable
  finalization, issuer authority, and scientific adequacy are established
  by fs-ledger/fs-package/fs-checker integration and deployment policy.
  The cross-crate no-mock E2E lives in `tests/lanes_e2e.rs` (dev-deps
  only): two independent theorem lanes plus one preregistered in-lane
  comparison drive fs-govern admission with every decision row persisted
  as a real fs-ledger event, the preregistration and refutation
  artifacts content-addressed in the ledger (the finalization receipt
  seals a hash that actually exists there), the outcome packaged as
  fs-package claims and re-checked solver-free by fs-checker (with a
  mismatched-root refusal probe), and the whole request sequence
  replayed byte-for-byte on a fresh portfolio ledger plus an idempotent
  full-retry pass. The G4 persistence-fault drill (same file) crashes
  the process mid-sequence — in-memory portfolio dropped, design-ledger
  handle closed — reopens the same path, proves pre-fault decision
  artifacts survive byte-for-byte, re-persists the recovered prefix as
  a DEDUPE (no double write), and converges the full re-execution to
  the never-crashed control bytes with idempotent retries absorbed.
  In-crate retention exhaustion (the governed store's own fault
  surface) is pinned by lane-009. NOT claimed: kernel-level I/O fault
  injection, torn/partial-write simulation, and media corruption —
  those live with fs-ledger's own crash-recovery battery and
  deployment policy.
- Independence classes are DECLARED. Canonicalization defeats cosmetic
  splits and the class backstop defeats partition gaming among honestly
  labeled lanes, but the crate cannot algorithmically prove that two
  falsifier families are genuinely independent — adversarial mislabeling
  is a governance-review matter, bounded in damage by the global caps.
- Non-`Clone` `&mut self` atomicity covers one in-process authority value and
  its retries; a
  multi-process admission service needs an external serialization or the
  ledger-backed successor.
