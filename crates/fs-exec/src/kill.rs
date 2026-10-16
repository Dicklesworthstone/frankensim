//! Statistical preemption: the kill-handle registry (plan §5.2 behavior 3,
//! Bet 8's machinery). The e-process layer (fs-eproc) holds handles to
//! candidate evaluation scope-trees; the moment an elimination e-value
//! crosses threshold it kills the candidate's ENTIRE tree mid-flight —
//! tile-pool runs, races, and solver drives sharing the candidate's gate
//! all drain at their next poll point, and arena scoping reclaims their
//! memory. Cancel-correctness guarantees no torn state; the reclaim-latency
//! histogram is measured per kill and ledgered (never assumed).

use crate::cx::CancelGate;
use core::fmt;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

/// A candidate's logical identity (the e-racing layer's key: generation ×
/// member, hashed however the optimizer likes).
pub type CandidateId = u64;

/// Registry of live candidate kill-handles. `Sync`; cheap to share.
/// Everything a candidate evaluates — kernels, races, solver drives — runs
/// under the gate registered here, so one `kill` reaches the whole tree.
#[derive(Debug, Default)]
pub struct KillRegistry {
    entries: Mutex<BTreeMap<CandidateId, Arc<CancelGate>>>,
}

/// A kill aimed at a candidate the registry does not know (wf9.8.1):
/// the structured missing-handle error — a claimed elimination that
/// exercised NO live evaluation tree must be visible, not an ignorable
/// `false`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnregisteredKill {
    /// The unknown candidate.
    pub id: CandidateId,
}

impl fmt::Display for UnregisteredKill {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "kill({}) found no registered handle: the candidate's evaluation tree was \
             never wired to this registry (or was already released) — register ids \
             before racing so eliminations reach live work",
            self.id
        )
    }
}

impl core::error::Error for UnregisteredKill {}

impl KillRegistry {
    /// Empty registry.
    #[must_use]
    pub fn new() -> Self {
        KillRegistry::default()
    }

    /// Panic-total lock: a poisoned registry degrades to its data
    /// (poisoning can only come from a `kill_where` predicate panic;
    /// gate state is atomic and never torn).
    fn lock(&self) -> MutexGuard<'_, BTreeMap<CandidateId, Arc<CancelGate>>> {
        self.entries.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// [`KillRegistry::kill`] for candidates that MUST be wired: the
    /// missing-handle case is a structured error, not a droppable
    /// boolean.
    ///
    /// # Errors
    /// [`UnregisteredKill`] when no handle is registered for `id`.
    pub fn kill_registered(&self, id: CandidateId) -> Result<(), UnregisteredKill> {
        self.registered_gate(id)?.request();
        Ok(())
    }

    /// Fetch a candidate gate without creating one. Admission paths use this
    /// to distinguish a caller-wired evaluation tree from a dummy gate that a
    /// tournament could otherwise create for itself.
    ///
    /// # Errors
    /// [`UnregisteredKill`] when no handle is registered for `id`.
    pub fn registered_gate(&self, id: CandidateId) -> Result<Arc<CancelGate>, UnregisteredKill> {
        self.lock().get(&id).cloned().ok_or(UnregisteredKill { id })
    }

    /// Register (or fetch) the kill-handle gate for a candidate. The
    /// candidate's ENTIRE evaluation tree must run under this gate for the
    /// kill to reach all of it.
    #[must_use]
    pub fn register(&self, id: CandidateId) -> Arc<CancelGate> {
        Arc::clone(
            self.lock()
                .entry(id)
                .or_insert_with(|| Arc::new(CancelGate::new())),
        )
    }

    /// Kill a candidate: request cancellation on its gate. Idempotent.
    /// Returns false for unknown ids (a structured non-event, not an
    /// error — the candidate may already be finished and released).
    pub fn kill(&self, id: CandidateId) -> bool {
        match self.lock().get(&id) {
            Some(gate) => {
                gate.request();
                true
            }
            None => false,
        }
    }

    /// Kill every candidate the predicate condemns (e-BH style batch
    /// elimination). Returns the killed ids, ascending (deterministic).
    ///
    /// The batch is a snapshot of the registry at call entry. A concurrently
    /// released gate can still be signalled through its in-flight `Arc`, while
    /// a candidate registered after the snapshot waits for the next batch.
    pub fn kill_where(&self, mut condemn: impl FnMut(CandidateId) -> bool) -> Vec<CandidateId> {
        // SAFETY: never invoke caller-controlled code while holding the registry mutex. A
        // predicate is allowed to query or mutate this registry without self-deadlocking.
        let entries: Vec<_> = self
            .lock()
            .iter()
            .map(|(&id, gate)| (id, Arc::clone(gate)))
            .collect();
        let mut killed = Vec::new();
        for (id, gate) in entries {
            if condemn(id) {
                gate.request();
                killed.push(id);
            }
        }
        killed
    }

    /// Release a finished candidate's handle (the gate lives on in any
    /// `Arc` still held by in-flight work; the registry just forgets it).
    /// Returns false for unknown ids.
    pub fn release(&self, id: CandidateId) -> bool {
        self.lock().remove(&id).is_some()
    }

    /// Live (registered, unreleased) candidate count.
    #[must_use]
    pub fn live(&self) -> usize {
        self.lock().len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;
    use std::time::Duration;

    #[test]
    fn register_kill_release_lifecycle_is_idempotent_and_structured() {
        let reg = KillRegistry::new();
        assert!(!reg.kill(7), "unknown id is a non-event");
        let gate = reg.register(7);
        let same = reg.register(7);
        assert!(Arc::ptr_eq(&gate, &same), "one gate per candidate");
        assert_eq!(reg.live(), 1);
        assert!(!gate.is_requested());
        assert!(reg.kill(7));
        assert!(reg.kill(7), "kill is idempotent");
        assert!(gate.is_requested());
        assert!(reg.release(7));
        assert!(!reg.release(7));
        assert!(!reg.kill(7), "released candidates are unknown");
        // The Arc'd gate outlives the registry entry (in-flight holders).
        assert!(gate.is_requested());
    }

    #[test]
    fn unregistered_kill_is_a_structured_error() {
        let reg = KillRegistry::new();
        let err = reg.kill_registered(42).expect_err("nothing wired");
        assert_eq!(err, UnregisteredKill { id: 42 });
        assert!(err.to_string().contains("no registered handle"));
        let _ = reg.register(42);
        assert!(reg.registered_gate(42).is_ok());
        reg.kill_registered(42).expect("wired candidates kill");
    }

    #[test]
    fn batch_elimination_kills_deterministically_by_ascending_id() {
        let reg = KillRegistry::new();
        for id in [5u64, 1, 9, 3] {
            let _ = reg.register(id);
        }
        let killed = reg.kill_where(|id| id % 2 == 1 && id > 2);
        assert_eq!(killed, vec![3, 5, 9], "ascending, deterministic");
        assert!(!reg.register(1).is_requested(), "survivor untouched");
        assert!(reg.register(3).is_requested());
    }

    #[test]
    fn batch_predicate_can_reenter_registry_without_deadlocking() {
        for _ in 0..10 {
            let reg = Arc::new(KillRegistry::new());
            let _ = reg.register(1);
            let worker_reg = Arc::clone(&reg);
            let (tx, rx) = mpsc::channel();
            std::thread::spawn(move || {
                let killed = worker_reg.kill_where(|id| {
                    let _ = worker_reg.register(2);
                    id == 1
                });
                let _ = tx.send(killed);
            });
            assert_eq!(
                rx.recv_timeout(Duration::from_secs(1))
                    .expect("reentrant predicate must complete"),
                vec![1]
            );
            assert_eq!(reg.live(), 2);
        }
    }
}
