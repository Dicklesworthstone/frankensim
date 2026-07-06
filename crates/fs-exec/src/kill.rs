//! Statistical preemption: the kill-handle registry (plan §5.2 behavior 3,
//! Bet 8's machinery). The e-process layer (fs-eproc) holds handles to
//! candidate evaluation scope-trees; the moment an elimination e-value
//! crosses threshold it kills the candidate's ENTIRE tree mid-flight —
//! tile-pool runs, races, and solver drives sharing the candidate's gate
//! all drain at their next poll point, and arena scoping reclaims their
//! memory. Cancel-correctness guarantees no torn state; the reclaim-latency
//! histogram is measured per kill and ledgered (never assumed).

use crate::cx::CancelGate;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

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

impl KillRegistry {
    /// Empty registry.
    #[must_use]
    pub fn new() -> Self {
        KillRegistry::default()
    }

    /// Register (or fetch) the kill-handle gate for a candidate. The
    /// candidate's ENTIRE evaluation tree must run under this gate for the
    /// kill to reach all of it.
    #[must_use]
    pub fn register(&self, id: CandidateId) -> Arc<CancelGate> {
        Arc::clone(
            self.entries
                .lock()
                .expect("kill registry")
                .entry(id)
                .or_insert_with(|| Arc::new(CancelGate::new())),
        )
    }

    /// Kill a candidate: request cancellation on its gate. Idempotent.
    /// Returns false for unknown ids (a structured non-event, not an
    /// error — the candidate may already be finished and released).
    pub fn kill(&self, id: CandidateId) -> bool {
        match self.entries.lock().expect("kill registry").get(&id) {
            Some(gate) => {
                gate.request();
                true
            }
            None => false,
        }
    }

    /// Kill every candidate the predicate condemns (e-BH style batch
    /// elimination). Returns the killed ids, ascending (deterministic).
    pub fn kill_where(&self, mut condemn: impl FnMut(CandidateId) -> bool) -> Vec<CandidateId> {
        let entries = self.entries.lock().expect("kill registry");
        let mut killed = Vec::new();
        for (&id, gate) in entries.iter() {
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
        self.entries
            .lock()
            .expect("kill registry")
            .remove(&id)
            .is_some()
    }

    /// Live (registered, unreleased) candidate count.
    #[must_use]
    pub fn live(&self) -> usize {
        self.entries.lock().expect("kill registry").len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
