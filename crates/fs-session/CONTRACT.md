# fs-session — CONTRACT

Sessions, capability tokens, and the resource GOVERNOR (plan §11.3):
budgets are ENFORCED, not advisory — plus the agent-proofing trio:
idempotency keys, `estimate()` dry runs, and errors as guidance.

Ambition tags: tokens/governor/idempotency/estimate [F per the bead
label; shipped surface tested to [S] discipline].

## Purpose and layer

Layer **L6** (HELM). Runtime deps: `std`, fs-ir (admission bridge +
study parsing), fs-exec (CancelGate/SolverState/TilePool), fs-la (production
GEMM), fs-ledger (persistence/content hash), fs-blake3 (domain-separated
receipt identity), fs-plan (cost models), fs-obs, fs-qty.
Consumers: the P2 marquee demo, the HELM e2e suite (gp3.11).

## Public types and semantics

- `CapabilityToken { session, ops globs, core_s, mem_bytes, wall_s,
  cores, ledger_scope }` — the explicit grant every IR program executes
  under; fallible `to_admission()` validates bounded canonical operator grants
  before cloning them into fs-ir static admission. The
  governor independently validates the declaration at registration and
  continuously meters core-seconds, peak memory, and wall time. Operator names
  and concurrent cores are static-admission fields today, not dynamically
  leased execution resources. Memory
  bytes and concurrent cores remain exact `u64` values through the bridge;
  source count scaling, token comparison, and governor enforcement never
  project those authorities through `f64`. `ledger_scope` is exact authority,
  admitted only as 1..=128 ASCII graphic bytes: whitespace, controls, Unicode
  normalization aliases, and oversized namespaces fail before registration.
  Invalid-scope diagnostics retain only a UTF-8-safe 128-byte preview plus the
  exact input length, so refusing an oversized authority string is itself
  memory-bounded. Operator authority is likewise structural: at most 256
  unique canonical exact names or namespace wildcards ending in `.*`, each
  1..=128 bytes and at most 8 KiB in aggregate, are admitted under fs-ir's
  shared grammar before any governor state changes. Validated strings and the
  operator vector are rebuilt before retention, so caller-chosen spare capacity
  cannot bypass byte accounting.
- `Governor` — `Send + Sync`; hot paths are mutex-guarded in-memory
  state. `open_session` returns a private-field `ScopeFlushPermit` bound to the
  exact governor and immutable ledger scope; callers cannot manufacture scope
  authority from strings. It rejects invalid ledger scope, non-finite, or
  negative floating grants before registration and rejects an already-open `SessionId`
  without replacing its immutable token or inheriting its meters, gate, pause
  state, or idempotency receipts. A pre-requested gate is refused before
  registration. Token and optional gate registration are one atomic critical
  section. `charge(session, Charge)` rejects non-finite,
  negative, or overflowing deltas before mutating meters, then meters
  core-seconds / peak memory / wall and returns `Enforcement`: `Ok` →
  `Throttled` (at the grant) → `Paused` (past 1.2× the grant, with a
  teaching resume hint). Memory admission compares integer bytes exactly,
  including above f64's 53-bit integer precision, and evaluates the hard
  threshold as `used * 5 > granted * 6`; diagnostic f64 fields do not drive the
  verdict. The governor NEVER silently kills.
- `submit_once(session, idem_key, work)` — exactly-once execution:
  the first caller in that session runs and is charged. A concurrent or
  same-thread reentrant caller observing `Pending` returns `InFlight`
  immediately without executing or waiting; a repeat after terminal
  publication receives `Duplicate` with the SAME receipt and NO charge.
  Caller panics and invalid returned charges are contained as a terminal
  `Failed { receipt, evidence }`; no charge is committed, and retry requires an
  explicit new key. A private-field `SubmissionReceipt` is a domain-separated BLAKE3
  v2 identity of the owning session, immutable ledger scope, exact key, and
  terminal charge plus enforcement decision or full failure evidence.
  Caller-controlled evidence retains at most 16 KiB of UTF-8-safe preview,
  while its exact byte length and BLAKE3 digest bind the complete original.
  Admission reserves three retained key copies plus worst-case terminal
  evidence before caller work begins, then releases unused terminal capacity
  when the result is published. A panic or invalid charge therefore cannot
  strand an unpublishable `Pending` generation at a byte boundary.
  `Executed` and every later `Duplicate` expose the same `Enforcement`, so a
  throttled/paused charge is never hidden behind a generic success.
  `idempotency_key(agent_key, program)` accepts at most 1 MiB per input,
  separately domain-hashes them, then binds exact lengths and the two digests
  into a fixed-memory v3 key; blank or oversized supplied execution keys are
  refused length-first.
- `apply_memory_pressure(session, level)` — the DECLARED
  degradation ladder (`LADDER`: spill coldest arenas → coarsen
  adaptively → pause-serialize-resume) fires steps `1..=level` in order;
  the pause step resolves the session-owned `CancelGate` bound by
  `open_session_gated` and requests it so the solver checkpoints at its next
  tile boundary (P7). Pause gates are generational: level 3 mints a
  private-field `PauseRequestId` bound to this governor, session, old gate
  generation, and request ordinal. `acknowledge_pause(request_id, claim)`
  consumes only that exact generation; while it remains the session's latest
  completed generation, identical response replay returns the same event and
  `Arc<CancelGate>`, while conflicting evidence fails closed. If an external
  owner requests that gate before activation, `activate_resume` refuses it and
  identical acknowledgement replay replaces the never-activated `Arc` at the
  same generation. The old acknowledgement becomes stale by exact pointer
  identity; repeated recovery replay returns the one current replacement.
  Completing a later generation invalidates older replay authority rather than
  retaining an unbounded gate history.
  `PauseAcknowledgement` has private fields and activation compares its full
  event, request, generation, and exact gate identity with retained state, so a
  caller-altered acknowledgement cannot activate work.
  The fresh gate remains `ReadyToResume` until
  `activate_resume(&acknowledgement)` records that resumed workers adopted it;
  all pressure levels refuse while a request is pending or a gate awaits
  activation. The old gate stays permanently requested so drained workers
  cannot re-enter. Request and completion events bind the old generation.
  Every event carries attribution and a deterministic ordinal. Admitting level
  3 reserves the mandatory future
  completion row, and all other event admission counts outstanding
  reservations so the completion cannot be starved at the cap. The request
  also reserves worst-case retained completion evidence and one global event
  ordinal before it requests the gate. Other event admission must preserve all
  outstanding completion reservations. Pause completion accepts at most 1 MiB
  of checkpoint-claim input and retains it under the bounded
  preview/full-length/full-digest evidence model, releasing unused reserved
  capacity atomically with completion.
  `events_page` is the only event reader and returns at most 1,024 rows under
  the permit's exact scope.
- `estimate(study, cost_models, cores)` — the dry run: p10/p50/p90 wall
  from fs-plan quantile models over `:dof`/`:size` features, declared
  memory ask, energy (p50 × cores × 45 W/core), and an HONEST
  `unmodeled_ops` coverage list (never silent gaps). The result is fallible:
  cores and derived wall/energy values must stay finite and non-negative, and a
  memory Count must scale to an exact positive whole-byte `u64`; zero, negative,
  fractional-byte, wrong-unit, and overflowed asks are structured refusals.
  Operation discovery includes namespaced verbs and any undotted verb present
  in the supplied model registry; registry-backed work cannot silently vanish
  from a dry-run estimate because of its spelling.
  Memory is read only from the recognized study's direct `(budget (mem ...))`
  clause: duplicate or malformed memory clauses fail closed, and a body call
  named `mem` has no authority to grant a budget.
  `unmodeled_ops` means no model exists; a present model that refuses its input,
  emits invalid numbers, or reverses its quantiles is an error, not silently
  relabeled as a coverage gap. An explicitly declared `:dof`, `:size`, or
  `:modes` feature must have one numeric value, and duplicate size features are
  refused instead of silently falling back to the unit-size default.
- `CalibrationReport` — estimate-vs-actual rows, ratio quantiles, and a
  content-addressed ledger artifact (`estimate-calibration`): the cost
  models' own report card. Zero-prediction rows (bead gp3.21) are
  EXCLUDED from ratio quantiles (no invented ratios) but never hidden:
  `zero_prediction_summary()` and the JSON's `zero_predictions` object
  carry their count — split into true-zero (fully modeled) vs unmodeled
  (coverage gap) — plus the raw actual-time quantiles; rows serialize as
  `[predicted, actual, fully_modeled]` triples. `health(&policy)` judges
  the zero-prediction fraction against the governance-configurable
  `CalibrationPolicy` threshold (default 0.25; non-finite or out-of-[0,1]
  thresholds refuse) and returns Healthy or Degraded, never a silent
  pass. Row ingestion rejects negative/non-finite values and
  non-finite ratios before mutation, so its canonical JSON cannot be poisoned
  by `NaN`/infinity spellings.
- `Guidance { code, diagnosis, fixes }` — errors as teaching:
  `from_finding` lifts fs-ir admission findings (the canonical §11.3
  `BudgetInfeasible` fixture) with their cost-model-ranked fixes intact.
- `flush_scope_to_ledger(&ScopeFlushPermit, &Ledger)` — changed consumption
  snapshots plus new degradation and terminal idempotency receipts for sessions
  whose immutable token grants that exact scope, persisted exactly once as an
  atomic bounded `session.*` event chunk. Foreign permits fail closed. Every payload
  carries the exact JSON-escaped `ledger_scope`; consumption/idempotency schemas
  are v3/v4 and degradation is v4. `FlushReport` names appended rows, encoded
  bytes, and whether more state remains dirty; each call admits at most 1,024
  rows and 4 MiB of conservatively encoded payload. An unchanged repeated call
  is a no-op; failed persistence leaves every selected generation-aware cursor
  dirty for retry.
  The call refuses an already-open ledger transaction because it cannot know
  whether the owner will commit or roll back. Preparation reserves one scope
  under the governor mutex, database I/O runs after releasing that mutex, and
  cursor commit validates the reservation plus generation/revision snapshot.
  A concurrent same-scope flush returns `ScopeFlushInFlight`; unrelated hot
  paths remain live. Each scope owns indexed dirty meter/receipt sets, an
  independent event cursor, sink, revision, and flush generation. A rotating
  three-lane start order prevents sustained meter traffic from starving
  receipts or events. Its first successful non-empty flush binds that scope to
  the ledger's opaque persisted `LedgerInstanceId`, revalidated against the
  live schema before every authority-bearing flush;
  aliases and moved handles remain the same sink, while a replacement file at
  the same path or independent memory ledger is refused before writes.
- Deterministic hard caps bound retained governor state and public
  materialization: 4,096 sessions/governor, 1,024 sessions/scope, 4,096
  idempotency keys/session, 65,536 degradation events/scope, 16 KiB evidence
  previews, 64 MiB of retained caller-controlled payload/scope, 256 MiB of
  retained caller-controlled payload/governor, 1,024 event-page and flush rows,
  4 MiB flush bytes, and checked signed-ledger ordinals. Counts bound fixed
  structure overhead; byte budgets conservatively count duplicated key strings,
  token text, event attribution/evidence, and terminal reservations. Limit+1
  refuses before partial mutation.
- `gemm_f64_session_with_pool(tuner, cache_policy, pool, gate, m, n, k, α,
  a, b, β, c)`
  — the production GEMM autotune loop (bead yqug): measure → cache →
  model → cancellation-correct dispatch through one caller-owned, reusable
  `TilePool`. `gemm_f64_session(..., threads, ...)` is the compatibility
  wrapper that constructs an unpinned host pool. The `*_budgeted` forms accept
  an explicit `GemmMemoryEnvelope`; legacy wrappers pass the explicit unbounded
  sentinel. The scoped key binds fs-la's bit semantics version, power-of-two
  shape class, requested/normalized thread budget, memory limit, exact capped
  probe dims (M/K ≤ 512, N ≤ 2048), resolved SIMD tier,
  the executing pool's canonical topology/mode/weights/arena/pin-groups
  identity, implementation version, and generated compiler/profile/codegen
  build fingerprint, plus `GEMM_TUNER_SCHEMA_VERSION`, which must bump whenever
  the producer lattice, probe/sample policy, ranking, or plan mapping changes;
  the ledger key also binds the machine fingerprint.
  `gemm_tune_build_evidence()` exposes that exact build fingerprint together
  with `GemmGraphEvidenceClass`, the fingerprint-bound class identity, and the
  optional full canonical dependency receipt + domain digest. This lets root
  orchestration require and retain the receipt artifact before publication.
  `GEMM_DEPGRAPH_RECEIPT_DOMAIN` is re-exported so the root can recompute the
  digest from retained bytes without copying a private domain string.
  `OperatorObservedReceipt` means strict receipt structure was validated, not
  that fs-session independently authenticated the operator or reconstructed
  Cargo's invocation-exact unit graph. `DevelopmentEquivalenceSalt` is never
  verified graph evidence.
  Shape/overflow and
  pre-requested-cancel checks happen before any tune mutation; one-thread,
  small-M, and no-product calls bypass tuning. Pinned plans skip measurement;
  else an exact tuner/ledger row is used; else an up-to-4×2 lattice is
  deduplicated by the `(mc,nc)` values fs-la will actually execute and sampled
  three times. Probe A, B, candidate C, and the exact-bit reference are
  fallibly reserved and jointly charged to the session envelope; bounded sweeps
  pass only the remaining child ceiling to fs-la. Every output word from every
  repeat is compared by `f64::to_bits`
  (signed zero and NaN payloads included); drift fails closed. The declared
  model is argmin of minimum wall time with lattice-order tie breaking, not a
  confidence claim. `GemmTuneCache` makes durable access explicit:
  `Disabled`, `ReadOnly`, or `ReadWrite`. Read-only callers may adopt an
  existing validated row but cannot publish during speculative work. A newly
  measured row is sealed as `ValidatedGemmTuneRow` inside `GemmDispatch`; its
  private fields cannot be forged or altered. `receipt_json()` is its canonical
  kernel/shape/machine/params/measured/memory-limit/probe-buffer-plan preimage;
  a public globally unique
  derive-key domain hashes those exact bytes, and
  `publish_to_ledger` participates in an already-open wider transaction.
  Cache adoption returns the same sealed identity on its first dispatch so
  downstream evidence can bind adopted and freshly measured rows uniformly.
  `publish_if_absent_or_identical` lets evidence populate an independent ledger
  but refuses to overwrite a different row already stored under that key.
  `replace_cache_row` is the distinct mutable-cache operation: only a sealed,
  remeasured row can replace stale or malformed dispatch state, and exact
  read-back is required. Replacement is never authority inherited by a delayed
  or cloned benchmark receipt. Read-write mode uses that repair path and
  durably writes the same sealed row before committing it locally. A decision is recorded
  only after fs-la drains and successfully commits the staged output.
  `GemmDispatch.run` carries the final compute counts and every real per-panel
  fs-exec `RunReport`; `execution_receipt()` projects kernel, mode, deterministic
  declared panel ordinal, completed/total counts, and the deterministic memory
  plan into `GemmExecutionReceipt`, explicitly excluding steal, latency,
  worker-distribution, peak-use, and refusal measurements from replay identity.
  `GemmDispatch.kernel` is the exact replay key; replay pins the recorded key
  and params rather than reconstructing a weaker base key.

## Invariants

1. **Enforcement is structured**: every over-grant outcome is `Throttled`
   or `Paused` with resource, used, granted, and a resume hint — no kill
   path exists in the API.
2. **Exactly-once within the owning session and scope**: for any `(session,
   idempotency-key)` pair, `work` runs at most once; all callers in that scope
   observe either non-blocking `InFlight` or the same terminal content-derived
   receipt, and consumption is charged exactly once (16-thread and reentrant
   race-tested). The receipt binds the immutable ledger scope. The same caller
   key in another session is independent and produces a different receipt, so
   one tenant cannot suppress another tenant's work.
3. **The ladder order is the contract**: spill before coarsen before pause,
   always; pause requests one gate generation, reserves its completion event,
   exact opaque acknowledgement rotates to a fresh gate, and explicit
   activation precedes the next pressure transition. Identical acknowledgement
   replay is idempotent while that completion is latest; a requested
   never-activated resume gate is recoverable by same-generation replay without
   another ledger event, and the replaced acknowledgement cannot activate.
   Activation replay is idempotent while its gate generation remains current,
   including after the next pause requests that gate. `SolverState` snapshots round-trip
   losslessly across repeated pause-resume cycles.
4. **Estimates state their coverage**: unmodeled ops are listed, their
   wall is excluded, nothing is silently assumed.
5. **Meters are exact under storm**: concurrent charges accumulate
   without loss (32-way storm test asserts exact totals).
6. **Every owned idempotency execution terminates**: success or caller panic
   transitions `Pending` exactly once and carries one shared terminal receipt;
   failed work never charges and same-key retry never executes implicitly.
   Successful receipts bind bit-exact charge fields and the resulting
   enforcement verdict; duplicates replay that verdict without recharging.
7. **Invalid resources fail closed**: NaN, infinities, negative values, and
   accumulated floating-point overflow are rejected before any token or meter
   mutation. Landing exactly on a grant returns `Throttled`.
8. **GEMM tuning cannot create phantom success**: malformed shapes and
   no-op/serial routes cannot create rows or decisions; cancellation or bit
   drift during measurement cannot create a row; cache failure, foreign
   execution identity, and params/body disagreement cannot install a row.
   Read-only mode performs no ledger writes and exports only an unforgeable
   validated row. In read-write mode, ledger persistence precedes local row
   commit; successful compute precedes decision commit. Cancellation during
   final dispatch may retain its already validated measured row, but records no
   successful decision and leaves caller `C` bitwise intact.
9. **Session identity is immutable**: opening an existing `SessionId` is a
   structured `SessionAlreadyOpen` refusal. The original capability, meter,
   cancellation gate, pause state, and terminal idempotency generations remain
   unchanged.
10. **Scoped flush is append-once and retryable**: one terminal submission
    generation, degradation event, or distinct meter snapshot is appended at
    most once to the sink bound by its token's exact ledger scope. Interleaved
    scopes have independent sink, generation, and degradation cursors. One
    scope's atomic success advances only its cursors; wrong-sink,
    foreign-permit, in-flight refusal, or persistence failure advances none,
    while successful unchanged repeats append zero rows. Each bounded chunk
    commits only its prepared cursors and reports whether another chunk or
    concurrently-created state remains dirty.
11. **All retained collections are bounded**: session, scope-session,
    operator-authority, idempotency-key, degradation-event plus reserved
    completion row/ordinal, checkpoint-input, evidence-preview,
    per-scope/per-governor retained payload, pagination, flush row/byte, and
    ordinal limits are checked
    before their corresponding state transition. Limit refusal never runs
    caller work or partially advances a cursor, gate, meter, event stream, or
    authority binding.

## Error model

`SessionError`: `UnknownSession`, `SessionAlreadyOpen`, `InvalidLedgerScope`,
`InvalidOperatorGrant`, `DuplicateOperatorGrant`,
`UnknownLedgerScope`, `LedgerScopeSinkMismatch`, `ScopePermitMismatch`,
`ScopeFlushInFlight`, `LimitExceeded`, `PauseAlreadyPending`,
`PreRequestedGate`, `PauseRequestMismatch`, `PauseAcknowledgementConflict`,
`ResumeNotActivated`, `ResumeAcknowledgementMismatch`,
`ResumeGateAlreadyRequested`, `InvalidResource`, `Submission`, `Persistence`.
`GemmTuneError`: cancellation with
the drained numerical report when dispatch began, structured TilePool failure
with tile provenance and its full report, typed tuner refusal, ledger refusal,
exact-bit drift with candidate and repeat, `MemoryRefused` with the outer
session peak plus any fs-la report, or `MemoryPlanOverflow` before unsafe
allocation. Cancellation observed before or between numerical dispatches has
no fs-la report but still retains the outer session peak; cancellation returned
by fs-la and every executor/memory refusal preserve the full drained report.
All such paths leave caller-visible `C` unchanged. Refusals that teach travel
as `Guidance` values with ranked fixes.
A caller-work panic is data, not an unwind across the governor API:
`SubmitOutcome::Failed` records its receipt and bounded retained evidence.

## Determinism class

Governor state transitions are deterministic given the call order;
event ordinals are logical (no wall clocks in ledgered payloads).
Concurrency outcomes (who wins a race) are scheduling-dependent by
nature — the INVARIANTS above are what is guaranteed. GEMM numerical bits are
independent of the selected MC/NC plan; the wall-clock winner is inherently
environment-dependent and therefore travels as scoped evidence plus an exact
replayable decision.

## Cancellation behavior

The governor is itself a cancellation SOURCE (pause step → CancelGate).
Its own operations are short, bounded critical sections. GEMM sweep and final
dispatch use the same caller-owned pool and poll the same gate inside bounded
packing/microtile work. fs-la stages `C`, stops claiming M-band tiles, drains
all workers and Cx arenas, and commits only after the final poll; cancellation
returns compute plus TilePool progress and leaves caller `C` bit-for-bit
unchanged. The pool's worker lifetime is not yet an asupersync child scope; the
precise no-claim and follow-up live in fs-exec's L0 contract.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-session/conformance`):
ss-001 token→admission bridge end-to-end; ss-002 Ok→Throttled→Paused
with exact meters and structured unknown-session errors; ss-003
16-thread idempotency race (one execution, one charge, non-blocking in-flight
observations or shared terminal receipt, plus same-thread reentrancy and
independent keys); ss-004 estimate p10/p50/p90 + energy + declared mem +
honest coverage, calibration ratio quantiles, ledgered artifact
round-trip; ss-004b invalid estimator/calibration numeric domains fail closed
without poisoning JSON; ss-005 ladder order + gate request + toy-SolverState
snapshot equality + attributed ordinal-ordered events; ss-006 the
canonical BudgetInfeasible finding surfacing as ranked `Guidance`;
ss-007 32-way adversarial-grant storm with exact meters and structured
outcomes only; ss-008 seeded caller panic with eight concurrent duplicates
returning `InFlight`, bounded full-digest evidence, one terminal failure
receipt, and zero charge; ss-009 NaN/infinite/negative grant and charge refusal
with no-mutation checks;
ss-010 the exact-grant throttle boundary and atomic accumulated-overflow
refusal; ss-002b duplicate session registration preserving the original token,
meters, gate, and terminal idempotency state; ss-003d atomic incremental ledger
flush, unchanged-call no-op behavior, and dirty-cursor retry after transaction
refusal; ss-003e canonical scope validation plus foreign-permit and page-bound
refusals; ss-003f two-scope interleaving, independent sink binding, exact scoped
payload escaping, cross-scope wrong-sink refusal, and per-scope cursor retry.
ss-003g proves moved handles retain opaque sink identity while independent
ledgers differ; ss-003h drains a limit+1 batch through bounded atomic chunks;
ss-011 covers pending/ready pressure refusal, cross-governor request authority,
bounded checkpoint evidence, identical response replay, same-generation
replacement of an externally requested never-activated gate, stale-
acknowledgement refusal, conflicting evidence refusal, explicit idempotent
activation, and repeated generations; ss-012
covers exact session/scope/key boundaries and atomic limit+1 refusal. In-module
tests cover exact flush row/byte limits, event/ordinal/retained-byte caps,
concurrent acknowledgement replay, altered-acknowledgement refusal, multiple
simultaneous reservations,
pagination, and same-scope in-flight reservation refusal.

`tests/gemm_tune.rs` (bead yqug drills): cold start sweeps once and
matches serial bits; ledger warm start seeds a fresh session without
re-measuring; stale (foreign-fingerprint) and invalid cache rows are
refused and re-measured; a requested gate cancels the sweep with no row
and untouched output; pre-requested warm and pinned paths leave output and
decision logs untouched; non-canonical/off-lattice/mis-keyed pins are
structured refusals; replay uses the actual recorded scoped decision and
reproduces the live bits; execution identity separates thread budgets and
exact probes even inside one shape bucket; serial/no-op/small and invalid-shape
paths cannot mutate tuning state; every lattice plan matches serial bits.
An n=640 producer test executes both the two-panel NC=512 route and the
single-panel wider route, proving the NC axis reaches fs-la rather than only
changing evidence labels. Caller-pool drills prove the same pool executes
measurement and final dispatch, the receipt contains real TilePool traversal,
legacy std-thread placement keys are refused, and pinned/unpinned placement
policies cannot share tune rows.
In-module injected Gauntlet tests force exact signed-zero/NaN-payload drift in
each repeat, candidate collapse, between-repeat cancellation, cache-write
failure/retry, params/body disagreement, wrong-probe adoption, and cancelled
dispatch without a success decision. The oracle lane
(`--ignored`, release) asserts the live choice is within the declared
25% tolerance of the exhaustive best-of-3 oracle at the real problem
size and reports its machine — measured 1.000/1.000/1.062 on
macos-aarch64 under ambient load; the second-ISA (x86) counterpart is
armed and runs when an x86 host picks it up.

## No-claim boundaries

- **The governor meters what it is TOLD** (`Charge` deltas from the
  executor); OS-level resource sampling and per-thread accounting are
  the executor/observability beads' territory.
- **Degradation steps are orchestration events**: actual arena spilling and
  adaptive coarsening are fs-alloc/solver behaviors triggered by these events,
  not implemented here. Their v4 phase is therefore `declared`, never
  `applied`. Pause requests, generational gates, exact request
  authority, replay, and activation are wired, but the checkpoint claim passed
  to `acknowledge_pause` is still operator-asserted. fs-session does not yet
  verify a ledger-minted solver-state artifact or an executor drain receipt, so
  `Complete` means "the owning orchestrator acknowledged this request," not an
  independently attested request-drain-finalize proof.
- **Capability issuance is not authenticated yet**: `CapabilityToken` and the
  fs-ir projection are caller-constructible declarations. Registration makes a
  declaration immutable inside one governor, but does not prove an external
  issuer, entitlement, expiry/revocation policy, dynamic operator binding, or a
  shared concurrent-core lease. This authority boundary is tracked by
  `frankensim-authenticate-session-capability-issuance-aeq7`.
- **Mutator retry and causal replay are incomplete**: `submit_once` and pause
  acknowledgement are keyed, but raw `charge`, pressure, and session-open calls
  are not. A lost successful open response can strand the caller without its
  private flush permit; exact token/gate replay and conflict detection remain
  required.
  Submission rows currently carry admission order while meter enforcement
  commits in completion order; they bind the returned outcome but do not yet
  contain the commit ordinal and pre/post meter snapshots needed to re-earn a
  concurrent throttle decision from row order. `submit_once` also accepts raw
  string keys without retaining an independent request/program digest, so
  callers must use the canonical helper and the governor cannot yet diagnose a
  same-key/different-request conflict. This is tracked by
  `frankensim-retry-idempotent-session-mutations-ujhp`.
- **Energy is a declared-constant model** (45 W/core), not measured
  power telemetry; the calibration channel is where reality lands.
- **Idempotency persistence is flush-based**: in-process registry +
  session-and-scope-bound ledgered success/failure receipts; cross-process
  replay reconstruction belongs to the HELM e2e/crash-recovery bead (gp3.11).
- **Each exact ledger scope flushes to one owning sink**: per-scope in-memory
  cursors prevent duplicate appends, while a later different sink for that
  scope is refused by opaque `LedgerInstanceId` before writes rather than
  receiving a partial history. Paths are not sink identity.
  Different scopes can bind independent sinks. Cross-ledger replication of one
  scope remains an event-log concern above this API.
- **Two-lane executor integration** (interactive vs batch lanes with
  core quotas) is deferred to gp3.11's study-scale batteries; the
  enforcement/idempotency/estimate surfaces here are what it composes.
- fs-session exposes fs-la's exact dependency receipt and trust class but does
  not upgrade it: correspondence between an operator-supplied receipt and the
  invoking Cargo build remains operator-trusted. Root publication from a
  development-salt build is not claimed as receipt-backed evidence.
- A mutex self-deadlock in the calibration renderer was found by the
  conformance run and fixed (single lock scope). That renderer remains
  non-reentrant; governor idempotency separately guarantees that same-thread
  duplicate submission returns `InFlight` rather than deadlocking.
- GEMM minimum-wall-time ranking is a deterministic selection rule over the
  recorded samples, not statistical confidence. The x86 oracle lane remains
  armed rather than claimed as measured until it runs on the reference host.
- The session envelope covers the four dominant numeric tune buffers, every
  fs-la-owned dispatch reservation, AND the session tune-metadata plan (bead
  wf9.15.1): candidate/ranking/observation collections (BTreeSet dedup was
  replaced by a bounded linear scan — tree-node overhead is not honestly
  accountable), one reused sample buffer plus per-observation exact copies,
  canonical plan labels, and the sealed-row strings, all with documented
  caps ENFORCED at observation/seal time. `run_sweep` charges the plan
  constant after the probe buffers clear the envelope alone (refusal
  `what: "tune-metadata-plan"`, before any allocation, never losing a
  validated fs-la report) and the child envelope excludes probe + plan
  bytes. The plan is a pure constant of the sweep lattice and schema caps —
  `gemm_tune_metadata_plan_bytes()` / `GEMM_TUNE_METADATA_PLAN_SCHEMA` — and
  every sealed row binds it into `receipt_json` as `metadata_plan`, so a
  freshly measured row and the same row adopted later derive the identical
  `receipt_identity`. Introducing the field rotates pre-plan row identities
  once (stored ledger cache rows re-tune on first touch — the same rotation
  class as a build-fingerprint change). Generic TilePool metadata remains
  the separate `frankensim-epic-substrate-wf9.16` boundary.
