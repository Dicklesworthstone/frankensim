//! Parked worker crew (bead tkr7): N OS threads park on a condvar INSIDE
//! an owning scope (a `std::thread::scope` or an asupersync
//! `Cx::scoped_cpu` — bead lx0e) and execute per-run jobs handed over by
//! reference, so repeated `TilePool` runs reuse warm workers instead of
//! respawning per run (the measured N-D FFT attainment collapse: every
//! axis pass respawned the pool, 0.011 attainment at 256×256 on a
//! 5995WX).
//!
//! CAPSULE: the ONE unsafe operation is the lifetime erasure of the job
//! borrow in [`Crew::dispatch`] — the rayon-latch argument: dispatch
//! publishes the borrow, wakes the crew, and BLOCKS until every worker
//! has reported completion and the slot is cleared, so the erased borrow
//! provably outlives every dereference. See SAFETY.md beside this file.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use asupersync::cx::CpuCx;
use std::sync::{Condvar, Mutex};

/// The job one run hands to the crew: called exactly once per crew worker
/// with the worker ordinal and that worker's own park-time task context
/// (the bridge that lets task cancellation/budget drain the run — lx0e).
pub(crate) type CrewJob<'a, Caps> = &'a (dyn Fn(usize, Option<&CpuCx<Caps>>) + Sync);

struct CrewState<Caps: 'static> {
    /// The published job. `'static` is ERASED, not true: the borrow is
    /// only valid until the dispatch that published it returns (SAFETY.md).
    job: Option<CrewJob<'static, Caps>>,
    /// Job generation; a worker executes a job at most once.
    epoch: u64,
    /// Workers still executing the current job.
    remaining: usize,
    /// Set once by [`Crew::shutdown`]; parked workers exit.
    shutdown: bool,
    /// First panic that escaped a job call (worker ordinal + payload).
    /// Kernel panics never reach here — the pool contains them per tile —
    /// so an entry is a pool invariant failure.
    invariant_panic: Option<(usize, String)>,
}

/// A crew of parked workers. The crew itself spawns NOTHING: the owner
/// spawns workers into its own scope and each runs [`Crew::park_loop`],
/// so worker lifetimes belong to that scope (P7) and the crew cannot
/// outlive or leak past it.
pub(crate) struct Crew<Caps: 'static> {
    workers: usize,
    state: Mutex<CrewState<Caps>>,
    /// Workers wait here for a job or shutdown.
    work_cv: Condvar,
    /// Dispatch waits here for the crew to finish the published job.
    done_cv: Condvar,
}

impl<Caps: 'static> Crew<Caps> {
    pub(crate) fn new(workers: usize) -> Self {
        Crew {
            workers,
            state: Mutex::new(CrewState {
                job: None,
                epoch: 0,
                remaining: 0,
                shutdown: false,
                invariant_panic: None,
            }),
            work_cv: Condvar::new(),
            done_cv: Condvar::new(),
        }
    }

    pub(crate) fn workers(&self) -> usize {
        self.workers
    }

    /// The worker half: park until a job (execute it, report completion)
    /// or shutdown (return, letting the owning scope join). Runs on a
    /// thread the OWNER spawned into its scope; `task_cx` is that
    /// worker's own scoped-CPU context when the owner is a task.
    pub(crate) fn park_loop(&self, w: usize, task_cx: Option<&CpuCx<Caps>>) {
        let mut seen = 0u64;
        loop {
            let job = {
                let mut st = self.state.lock().expect("crew state");
                loop {
                    if st.shutdown {
                        return;
                    }
                    if st.epoch != seen
                        && let Some(job) = st.job
                    {
                        seen = st.epoch;
                        break job;
                    }
                    st = self.work_cv.wait(st).expect("crew state");
                }
            };
            // The job runs OUTSIDE the state lock. Kernel panics are
            // already contained per tile inside the pool's worker loop;
            // anything escaping here is a pool invariant failure,
            // captured so the crew's completion latch can never hang.
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                job(w, task_cx);
            }));
            let mut st = self.state.lock().expect("crew state");
            if let Err(payload) = outcome {
                let message = payload
                    .downcast_ref::<&str>()
                    .map(ToString::to_string)
                    .or_else(|| payload.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "non-string panic payload".to_string());
                if st.invariant_panic.is_none() {
                    st.invariant_panic = Some((w, message));
                }
            }
            st.remaining -= 1;
            if st.remaining == 0 {
                // Revoke the erased borrow BEFORE releasing the caller:
                // dispatch returns only after this store (SAFETY.md).
                st.job = None;
                self.done_cv.notify_all();
            }
        }
    }

    /// The caller half: publish one job, wake the crew, and BLOCK until
    /// every worker has finished it. Returns the first invariant panic
    /// (worker ordinal + message) if one escaped a job call.
    ///
    /// # Panics
    /// On overlapping dispatches or dispatch after shutdown (programmer
    /// errors — the crew executes exactly one job at a time).
    pub(crate) fn dispatch(&self, job: CrewJob<'_, Caps>) -> Option<(usize, String)> {
        // SAFETY: the erased borrow is published under the state lock,
        // dereferenced by workers only between this publish and their
        // completion report, and revoked (slot cleared) by the LAST
        // completion — which this function blocks on before returning.
        // The borrow therefore outlives every dereference; no worker can
        // observe it after dispatch returns (a fresh job requires a new
        // epoch published by a later dispatch). See SAFETY.md.
        let erased: CrewJob<'static, Caps> =
            unsafe { std::mem::transmute::<CrewJob<'_, Caps>, CrewJob<'static, Caps>>(job) };
        {
            let mut st = self.state.lock().expect("crew state");
            assert!(
                st.job.is_none() && st.remaining == 0,
                "crew dispatch overlap (programmer error)"
            );
            assert!(
                !st.shutdown,
                "crew dispatch after shutdown (programmer error)"
            );
            st.job = Some(erased);
            st.epoch = st.epoch.wrapping_add(1);
            st.remaining = self.workers;
        }
        self.work_cv.notify_all();
        let mut st = self.state.lock().expect("crew state");
        while st.remaining > 0 {
            st = self.done_cv.wait(st).expect("crew state");
        }
        st.invariant_panic.take()
    }

    /// Wake every parked worker and tell it to exit. Idempotent. Called
    /// by [`CrewShutdown`] on drop — including when the owner's scope
    /// body unwinds — so the owning scope's join can never hang on
    /// parked workers.
    pub(crate) fn shutdown(&self) {
        let mut st = self.state.lock().expect("crew state");
        st.shutdown = true;
        drop(st);
        self.work_cv.notify_all();
    }
}

/// Shuts the crew down on drop. The owner declares this INSIDE its scope
/// body (after spawning the workers, before running user code) so both
/// normal return and unwind release the parked workers for join.
pub(crate) struct CrewShutdown<'a, Caps: 'static>(pub(crate) &'a Crew<Caps>);

impl<Caps: 'static> Drop for CrewShutdown<'_, Caps> {
    fn drop(&mut self) {
        self.0.shutdown();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::cx::cap;
    use std::sync::atomic::{AtomicU64, Ordering};

    fn crewed<R>(workers: usize, f: impl FnOnce(&Crew<cap::All>) -> R) -> R {
        let crew: Crew<cap::All> = Crew::new(workers);
        std::thread::scope(|s| {
            let _guard = CrewShutdown(&crew);
            for w in 0..workers {
                let crew = &crew;
                s.spawn(move || crew.park_loop(w, None));
            }
            f(&crew)
        })
    }

    #[test]
    fn every_worker_runs_each_dispatched_job_and_borrows_locals() {
        crewed(4, |crew| {
            let hits: Vec<AtomicU64> = (0..4).map(|_| AtomicU64::new(0)).collect();
            for round in 1..=3u64 {
                crew.dispatch(&|w, _cpu| {
                    hits[w].fetch_add(round, Ordering::Relaxed);
                });
                for h in &hits {
                    assert_eq!(h.load(Ordering::Relaxed), round * (round + 1) / 2);
                }
            }
        });
    }

    #[test]
    fn job_panic_is_captured_and_the_latch_still_completes() {
        crewed(3, |crew| {
            let refusal = crew.dispatch(&|w, _cpu| {
                assert!(w != 1, "crew invariant probe");
            });
            let (worker, message) = refusal.expect("panic captured");
            assert_eq!(worker, 1);
            assert!(message.contains("crew invariant probe"), "{message}");
            // The crew survives an invariant panic: next dispatch is clean.
            assert!(crew.dispatch(&|_w, _cpu| {}).is_none());
        });
    }

    #[test]
    fn unwinding_through_the_owner_scope_releases_parked_workers() {
        let outcome = std::panic::catch_unwind(|| {
            crewed(2, |crew| {
                crew.dispatch(&|_w, _cpu| {});
                panic!("owner body unwinds mid-session");
            })
        });
        // The proof is termination: CrewShutdown released the parked
        // workers during unwind, so the scope joined instead of hanging.
        assert!(outcome.is_err(), "the owner panic propagates");
    }

    #[test]
    fn shutdown_is_idempotent_and_workers_exit_cleanly() {
        crewed(2, |crew| {
            crew.dispatch(&|_w, _cpu| {});
            crew.shutdown();
            crew.shutdown();
        });
    }
}
