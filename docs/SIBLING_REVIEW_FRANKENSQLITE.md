# Independent Sibling Review — FrankenSQLite durability boundary

Bead: `frankensim-extreal-program-f85xj.13.5`. Review id: **SREV-2026-07-B**.
Reviewer: GoldSnow, 2026-07-25. Drills:
`crates/fs-ledger/tests/sibling_review_durability.rs`.

Companion to [`SIBLING_REVIEW_ASUPERSYNC.md`](SIBLING_REVIEW_ASUPERSYNC.md),
which reviews the cancellation half and states the shared method and
independence grading. Both apply here unchanged: contract-first, minimum
credible independence (no authorship history on the target), not a genuinely
external audit.

## Claims under review

From `crates/fs-ledger/CONTRACT.md` — FrankenSim's usage assumptions about the
engine underneath — scoped by the `f85xj.13.1` P1 list to WAL durability,
transaction boundaries, migration refusal, and blob/checkpoint behaviour:

| Claim | Why it is load-bearing |
|---|---|
| `FutureSchema`: a newer file is **refused, never clobbered** | a data-loss claim: a ledger owned by a newer client must survive an old client touching it |
| Seals are atomic, idempotent for the same op, and immutable | every replay and provenance claim rests on seals being real |
| A refused seal leaves no seal | a dangling seal row would assert authorship that was rejected |
| Unknown subjects remain `NotFound` | a missing subject must not read as an accepted one |
| Transaction boundaries hold | work inside a rolled-back transaction must vanish |

### Sharp edge deliberately targeted

The **prepare-time subquery folding** incident (`904cdef1`, beads cgt07 /
lnbzs): FrankenSQLite never binds parameters referenced inside WHERE-clause
subqueries of a FROM-less `SELECT`, so `seal_artifact_output`'s guarded
`INSERT ... SELECT ?1, ?2 WHERE EXISTS(... ?1 ...)` **silently inserted zero
rows for every valid sole producer**. The seal returned `Ok` and sealed nothing.
Regressions of fixed bugs are prime review targets, and this one is unusually
instructive — see D2 below.

## Drills and results

`cargo test -p fs-ledger --test sibling_review_durability` → **5 passed, 0
failed**.

| Drill | Claim attacked | Result |
|---|---|---|
| D1 | `FutureSchema` refused, durable state byte-identical | pass |
| D2 | a valid seal is **observably persisted** (seeded regression) | pass |
| D2b | negative control: sealing a non-producer refuses, leaves no seal | pass |
| D3 | unknown artifact refuses, no dangling row | pass |
| D4 | rolled-back transaction leaves no trace | pass |

### Why D2 is shaped the way it is

A test that only checked "an invalid seal refuses" would have **passed during
the bug**: zero rows inserted means there is nothing to violate. The failure
mode was *silent success*, so the only discriminating assertion is that the
**valid** case leaves an observable seal — the drill seals, then reads the seal
back and requires `Some(op)`.

This generalises. When a historical defect was a silent no-op, the regression
drill must assert a **positive observable**, not merely that the error path
still errors. Testing only the negative path is how a silent-success bug hides
from its own regression suite.

## Findings

**No defects found at the reviewed boundary.**

One finding about *this review's own method*, recorded because it nearly
produced a false defect report:

> **The first version of D1 was wrong and reported a defect that does not
> exist.** It stamped a future schema version by patching byte offset 60 of the
> main database file — the SQLite header's `user_version` field. The drill then
> observed the ledger opening the file happily and was one step from being
> written up as "FutureSchema refusal is not enforced".
>
> The real explanation: this ledger runs in **WAL mode**. A 519 KB
> `ledger.db-wal` sidecar held a newer copy of page 1 — the page containing the
> header — so the engine read `user_version` from the WAL and never saw the
> patched byte. The stamp must go through the engine (`PRAGMA user_version`) so
> it lands wherever the engine will actually look.
>
> A false certificate is worse than an ordinary wrong answer, and a false
> *defect report* against a sibling is the review-shaped version of that. The
> corrected helper carries this explanation in its doc comment so the next
> reviewer does not rediscover it the same way.

The same discovery hardened D1's other half: the "never clobbered" assertion now
snapshots the **whole database** — main file plus every sidecar — because
comparing only the main file would let a refusal that journals into the WAL pass
as untouched.

## Residual concerns — what this review does NOT cover

1. **Only the fs-ledger adapter surface was drilled.** Nothing here exercises
   FrankenSQLite's pager, WAL replay, B-tree, or locking internals.
2. **No crash or interruption drills.** The charter called for torn-write
   recovery, mid-transaction `kill -9`, and lock-contention storms. None were
   run: they need process-level orchestration, not in-process tests. These are
   the drills most likely to find a real durability defect, and they remain
   undone.
3. **No concurrency drills.** `Busy` / write-conflict retry semantics are
   claimed by the contract and untested here.
4. **Checkpoint behaviour is untested.** D1 observes that a WAL exists; it does
   not test checkpointing, WAL truncation, or recovery from a partial
   checkpoint.
5. **The engine is not certified by green drills.** They certify that the
   boundary FrankenSim consumes behaves as documented at the observable surface.

## Status

This completes the second of the two priority reviews the bead requires. Both
halves found no defect at their reviewed boundaries, and both recorded what they
tried and what they could not reach — which, per the bead, is the deliverable
when a review finds nothing.
