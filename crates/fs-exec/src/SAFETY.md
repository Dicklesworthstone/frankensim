# SAFETY: fs-exec/crew.rs parked-crew job-borrow erasure (bead tkr7)

## What the unsafe does
One `mem::transmute` in `Crew::dispatch` erasing the lifetime of a job
borrow: `&'run (dyn Fn(usize, Option<&CpuCx<Caps>>) + Sync)` →
`&'static (dyn ...)`, so the borrow can sit in the crew's `Mutex`-guarded
job slot while parked workers (spawned by the OWNER into its own scope)
dereference it.

## Why it exists
`TilePool` respawns scoped workers on every run; N-D FFT axis passes run
the pool 4–6 times per transform, and the spawn/join cost collapses
small-size attainment (measured 0.011 at 256×256 on a 5995WX, bead 27d3).
Parked workers need a job channel; a job that borrows the run's stack
frame (kernel, deques, slots) cannot cross a channel without erasing the
lifetime — the same argument rayon's scoped registry makes.

## Invariants
1. The erased borrow is dereferenced ONLY between its publication in
   `dispatch` and the completion report of the last worker executing it.
2. `dispatch` BLOCKS until `remaining == 0`, and the last worker clears
   the job slot (`st.job = None`) before signalling `done_cv` — so the
   slot never holds the borrow after `dispatch` returns, and the referent
   (the caller's stack frame) outlives every dereference.
3. A worker executes a job at most once: it captures the job only when
   `epoch` advances past its `seen` counter, and `epoch` only advances
   under the state lock in `dispatch`.
4. At most one job is in flight: `dispatch` asserts `job.is_none() &&
   remaining == 0` before publishing (overlap is a structured panic, not
   UB — the slot is typed, never dangling).

## Aliasing assumptions
The job is a `&dyn Fn` — shared, immutable, `Sync`. Workers alias it
freely by design; no `&mut` exists anywhere in the capsule. Interior
state the job touches (deques, slots, atomics) carries its own
synchronization in `pool.rs`, outside this capsule.

## Alignment assumptions
None beyond what the borrow already guarantees; the capsule performs no
raw-pointer arithmetic and no allocation.

## Lifetime assumptions
`'run` (the true lifetime of the job borrow, the duration of one
`dispatch` call) is erased to `'static` and dynamically re-established by
the completion latch: publish → workers use → last worker revokes →
dispatch returns. The referent lives in the frame of `TilePool`'s
`run_inner`, which is blocked inside `dispatch` for the entire window.

## Panic behavior
Job calls run under `catch_unwind` in `park_loop`; a panic is recorded
(first worker + payload) and the worker still decrements `remaining`, so
the latch NEVER hangs on a panicked job and the borrow is still revoked
by the last completion. Kernel panics never reach this capsule — the
pool contains them per tile — so a captured panic here is a pool
invariant failure, surfaced by `dispatch`'s return value. The
bookkeeping between lock acquisitions is panic-free (integer ops,
`Option` writes), so the state mutex cannot poison in practice; `expect`
converts pathological poisoning into a structured panic, not UB.

## Cancellation behavior
Workers park on a condvar, NOT on asupersync checkpoints; cancellation
of the owning task drains RUNNING jobs at tile boundaries (each worker's
own `CpuCx` is passed into every job call — the lx0e bridge, outside
this capsule). Parked workers are released by `Crew::shutdown`, which
the `CrewShutdown` guard runs on BOTH normal return and unwind of the
owner's scope body, so a cancelled/unwinding owner can always join its
scope. No capsule state exists at a cancellation point that could tear:
the job slot is either published-and-latched or `None`.

## Concurrency behavior
All shared state (`job`, `epoch`, `remaining`, `shutdown`,
`invariant_panic`) lives in ONE `Mutex<CrewState>`; both condvars
(`work_cv`, `done_cv`) wait on that mutex with predicate re-check loops,
so lost and spurious wakeups are handled by construction. The mutex's
release/acquire edges order the job publication before any worker
dereference, and every dereference before the slot clear that `dispatch`
observes. `CrewJob<'static, Caps>` is `Send`/shareable because the
underlying `dyn Fn` is `Sync` (bound in the type alias).

## Miri coverage
The crew tests (`every_worker_runs_each_dispatched_job_and_borrows_locals`,
`job_panic_is_captured_and_the_latch_still_completes`) exercise the
erasure + latch under Miri when the fs-exec Miri lane runs; the capsule
uses no intrinsics or FFI, so Miri sees every access. OS-thread condvar
timing under Miri is slow but sound (bounded rounds: 3 dispatches).

## Model-checking status and no-claim
The single-mutex predicate loops above are an invariant argument, not an
exhaustive enumeration of schedules. Neither this parked-crew protocol nor
the surrounding `TilePool` protocol, arena-lease handoff, or fault-drain path
currently has exhaustive schedule/model-checking coverage.

Current source-level evidence is bounded runtime coverage: the crew unit tests
exercise publication, completion, panic, and shutdown cases; the pool storms
exercise cancellation drain, panic containment, and reuse; and
`tests/fault_storm.rs` declares a fixed 16-seed battery whose receipts identify
the selected seed, tile, and touch. Those cases can falsify defects on the
schedules they reach, but they do not enumerate every mutex/condvar
interleaving.

asupersync's own model-checking results apply only to the protocols represented
in its models. They do not transitively cover FrankenSim's `Crew`, `TilePool`,
arena leases, or injected-fault drains. A future G4 model must represent job
publication and revocation, both condition variables, owner unwind, lease
release, and drain finalization before this capsule can claim exhaustive
model-checking evidence. Until then, the claim is limited to the bounded
seeded/runtime batteries above.

## Fuzz/property coverage
The safe facade is exercised end-to-end by the pool's G5 determinism
audits (parked lane bitwise-equal to the spawned lane across worker
counts and reruns) — any erasure misuse that corrupted a run would break
bitwise equality or the slot-completeness check (`RunError::Incomplete`).

## Proof obligations discharged by callers
None. The facade (`dispatch`) upholds the latch discipline internally;
callers cannot retain the borrow (it is consumed by the call), cannot
overlap dispatches (asserted), and cannot leak workers past their scope
(the owner's scope joins them; `CrewShutdown` guarantees release).
