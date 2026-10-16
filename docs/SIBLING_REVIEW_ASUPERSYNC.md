# Independent Sibling Review — asupersync cancellation boundary

Bead: `frankensim-extreal-program-f85xj.13.5`. Review id: **SREV-2026-07-A**.
Reviewer: GoldSnow, 2026-07-24. Drills:
`crates/fs-exec/tests/sibling_review_cancellation.rs`.

## Why this review exists

Common authorship across the constellation means an assumption shared between
FrankenSim and a sibling is invisible to every internal test on both sides. Both
suites can be green and both can be wrong in the same direction. For components
*under* the certificates — asupersync's cancellation protocol, frankensqlite's
durability — independent scrutiny is the only cure.

## Independence grading — honest

The bead defines the minimum credible standard as review by an agent with **no
authorship history on the target**, working **against the target's own claimed
contracts**, with adversarial drills written **from the contract alone before
reading the implementation**.

This review meets that standard, with one caveat stated plainly:

| Criterion | Status |
|---|---|
| No authorship history on asupersync | **Met** — I have never written asupersync code |
| Claims taken from the target's own docs | **Met** — `asupersync_v4_formal_semantics.md` §9.2 and the rule index |
| Drills written before reading the implementation | **Met** — claims enumerated first; only the *public API signatures* of `CancelGate`/`Cx` were read, to make the drills runnable |
| Genuinely external reviewer | **Not met** — I am another agent in the same fleet, under the same operator |

The caveat matters. This is the *minimum credible* bar, not the ideal one. A
reviewer sharing the fleet's habits of mind can still share its blind spots.

## Claims under review

From asupersync's own formal semantics (§9.2 proof obligations and the
cancellation rule index at §1):

| Claim | Source |
|---|---|
| `request → drain → finalize` is **monotone and idempotent** | obligation 5; rule index #1–#4 |
| `inv.cancel.idempotence` | rule index #5 |
| `inv.cancel.propagates_down` | rule index #6 |
| `rule.cancel.checkpoint_masked` | rule index #10 |
| No obligation leaks: `TaskCompleted t → Held(t) = ∅` | obligation 3 |
| Loser drain: race completion implies all losers completed | obligation 4 |
| Cleanup budget bounds drain/finalize work | §545 |
| Bounded cancel fairness: non-cancel task dispatched within `L+1` | obligation 7 |

From FrankenSim's usage assumptions (`crates/fs-exec/CONTRACT.md`):

| Assumption | Note |
|---|---|
| The latency lane inherits the region state machine **unmodified** | the load-bearing inheritance |
| Poll order is fixed: deadline observation → one poll spend → cancellation observation | |
| Deadline expiry is retained as the **first** terminal failure | |
| A successful receipt **cannot** cross an observed absolute deadline | strongest falsifiable claim |
| Clock-free manual gates produce an empty latency sample set and **make no latency claim** | |

### Known sharp edges deliberately targeted

- **Caller-owned cancellation gates.** A historically fixed defect: a race,
  pause, or memory-pressure response could manufacture *private* cancellation
  state the owner could not observe. Regressions of fixed bugs are prime review
  targets, so this is drilled directly (D5) with a negative control (D5b).
- **First-timestamp retention.** `CancelGate::request` claims the *first*
  request's timestamp is what latency histograms measure from. A re-stamping
  regression would corrupt every cancel-latency measurement in the workspace
  while still producing plausible numbers.

## Drills and results

All seven pass. `cargo test -p fs-exec --test sibling_review_cancellation` →
**7 passed, 0 failed**.

| Drill | Claim attacked | Result |
|---|---|---|
| D1 | idempotence; first-timestamp retention | pass — 64 repeat requests do not move the retained stamp |
| D2 | `rule.cancel.checkpoint_masked` | pass |
| D3 | monotonicity | pass — no revert across 10,000 polls |
| D4 | clock-free gates make no latency claim | pass — no fabricated or accruing measurement |
| D5 | caller-owned gate (seeded historical defect) | pass |
| D5b | **negative control** for D5 | pass — a detached gate is observably different |
| D6 | `inv.cancel.propagates_down` across contexts | pass |

D5b exists because a suite in which everything passes proves nothing on its own.
It seeds the regression's observable — a context bound to a gate that is not the
caller's — and confirms the boundary behaves measurably differently, so D5's
assertion is demonstrably the discriminating one. If D5b ever starts reporting
cancellation, the two situations have become indistinguishable and D5 has
silently stopped testing anything.

## Findings

**No defects found at the reviewed boundary.** Per the bead, a review that found
nothing must still record what it *tried*; the drill table above is that record.

One observation, not a defect: `CancelGate::new_clock_free` documents an
internal sentinel value of `1` for a request marker while `requested_at_ns`
remains the accessor. The distinction between "sentinel" and "measurement" is
carried by documentation rather than by the type. It is correct today and D4
pins the behaviour, but a typed `LatencySample::None` would make the no-claim
boundary unforgeable rather than conventional. Not filed as a sibling defect —
it is a hardening suggestion at FrankenSim's own adapter.

## Residual concerns — what this review does NOT cover

Stating these is the point of the exercise; the drills above are narrow.

1. **Only the adapter surface was drilled.** Every claim was tested through
   `fs-exec`'s public `CancelGate`/`Cx`. Nothing here exercises asupersync's
   scheduler, its region state machine internally, or its own task graph.
2. **Four claims are UNTESTED here.** No obligation leaks (obligation 3), loser
   drain (obligation 4), cleanup-budget bounding, and bounded cancel fairness
   (obligation 7) all require driving a real runtime with concurrent tasks.
   They are the claims most likely to hide a defect, and they remain unreviewed.
3. **The strongest FrankenSim assumption is untested.** "A successful receipt
   cannot cross an observed absolute deadline" needs a deadline expiring
   precisely at publication — a timing-race drill, not a state drill.
4. **No crash/interruption drills.** No `kill -9` mid-transaction, no
   panic-during-finalize, no nested-scope leak hunt.
5. **Passing drills do not certify the sibling.** They certify that the boundary
   FrankenSim consumes behaves as documented at the observable surface.

## frankensqlite half — BLOCKED, not skipped

The bead orders frankensqlite first. It could not be reviewed this session:
**frankensqlite does not build.** `cargo build -p fs-ledger` fails with E0053
trait-signature mismatches on `get_page` / `write_page` / `write_page_data` /
`restore_staged_page_data`, an in-flight async migration in `fsqlite-pager`
against the moving asupersync sibling. Contract-first durability drills (torn
write recovery, mid-transaction kill, lock-contention storms) all require a
working build.

The charter for that half is nonetheless fixed in advance, which is the useful
part of writing it now — the claims are chosen before the implementation can
influence them:

- WAL durability and crash recovery; transaction boundary semantics; migration
  refusal behaviour; blob/checkpoint behaviour (the `f85xj.13.1` P1 list).
- Sharp edge to target: the **prepare-time subquery folding / `IS NULL`**
  incident class (`params=None`, fixed upstream `2bd64d114`). A regression of a
  fixed bug is a prime review target, and this one silently returned wrong rows
  rather than failing.

Tracked as the open remainder of this bead.
