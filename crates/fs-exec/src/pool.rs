//! The throughput lane: a work-stealing fork-join tile pool (plan §5.2).
//!
//! Semantics first, lock-freedom later: workers own per-worker deques
//! (`CachePadded<Mutex<VecDeque>>`) seeded with contiguous, weight-
//! proportional tile ranges; an empty worker steals HALF a victim's deque,
//! visiting same-CCD victims before cross-CCD ones (plan §5.1 consequence
//! 3). The protocol — weighted quanta, CCD-local-first stealing, fixed-slot
//! reductions, drain-on-cancel, panic containment — is the contract; the
//! Chase–Lev lock-free deque is a later optimization gated on roofline
//! evidence (CONTRACT no-claims).
//!
//! Determinism (P2): every tile's output lands in its OWN slot and slots
//! fold in ascending tile order, so results are bit-identical across worker
//! counts and steal schedules by construction. RNG stream keys derive from
//! logical identity only.

use crate::cx::{CancelGate, Cx, ExecMode, StreamKey};
use crate::kernel::{Reduce, TileKernel};
use core::fmt;
use core::ops::ControlFlow;
use fs_alloc::CachePadded;
use fs_substrate::affinity::CcdTopology;
use std::collections::VecDeque;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

/// Pool configuration. Normalized (not rejected) by [`TilePool::new`]:
/// `workers` is clamped to at least 1 and `quantum_weights` is resized to
/// `workers` (missing entries take weight 1, zero weights are raised to 1).
#[derive(Debug, Clone)]
pub struct PoolConfig {
    /// Worker count (defaults to available parallelism at the call site's
    /// discretion; the pool itself never probes).
    pub workers: usize,
    /// CCD/cluster shape used to derive the steal order (fixtures or
    /// `CcdTopology::from_probe`).
    pub topo: CcdTopology,
    /// Per-worker initial-share weights — the P/E asymmetry hook: E-core
    /// workers get proportionally smaller tile quanta instead of being
    /// ignored or stalling joins. Weights come from the autotuner
    /// eventually; explicit until then.
    pub quantum_weights: Vec<u32>,
    /// Study seed (the Five Explicits' seed pillar) for stream keys.
    pub seed: u64,
    /// Execution mode, stamped on reports and events.
    pub mode: ExecMode,
    /// Arena configuration for per-tile scope arenas.
    pub arena: fs_alloc::ArenaConfig,
}

impl PoolConfig {
    /// A sane default: `workers` workers, weight 1 each, deterministic mode.
    #[must_use]
    pub fn new(workers: usize, topo: CcdTopology, seed: u64) -> Self {
        PoolConfig {
            workers,
            topo,
            quantum_weights: Vec::new(),
            seed,
            mode: ExecMode::Deterministic,
            arena: fs_alloc::ArenaConfig::default(),
        }
    }
}

/// Structured run failure (Decalogue P10). Cancellation and panics are
/// OUTCOMES, never process aborts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RunError {
    /// The run's cancel gate was requested; workers drained cleanly.
    Cancelled {
        /// Kernel name.
        kernel: &'static str,
        /// Tiles that completed before the drain.
        completed: u64,
        /// Total tiles planned.
        total: u64,
    },
    /// A tile panicked; siblings were cancelled and drained; the pool
    /// remains usable.
    TilePanicked {
        /// Kernel name.
        kernel: &'static str,
        /// The offending tile (full provenance for the ledger).
        tile: u64,
        /// The panic payload's message, when it carried one.
        message: String,
        /// Tiles that completed despite the failure.
        completed: u64,
    },
    /// Defensive: a slot was missing at fold time (executor bug, reported
    /// structurally rather than panicking across the boundary).
    Incomplete {
        /// Kernel name.
        kernel: &'static str,
        /// First missing tile slot.
        tile: u64,
    },
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RunError::Cancelled {
                kernel,
                completed,
                total,
            } => write!(
                f,
                "kernel `{kernel}` cancelled after {completed}/{total} tiles; partial work was \
                 reclaimed with the scope arenas (request -> drain -> finalize)"
            ),
            RunError::TilePanicked {
                kernel,
                tile,
                message,
                completed,
            } => write!(
                f,
                "kernel `{kernel}` tile {tile} panicked: {message} ({completed} sibling tiles \
                 completed; siblings were cancelled, the pool remains usable)"
            ),
            RunError::Incomplete { kernel, tile } => write!(
                f,
                "kernel `{kernel}` finished without output for tile {tile}: executor invariant \
                 violation — please report this"
            ),
        }
    }
}

impl core::error::Error for RunError {}

/// Measured facts about one run: steal statistics and the cancel-latency
/// samples (ns between the first cancel request and each worker OBSERVING
/// it at a tile boundary). Measurements only — results never depend on them.
#[derive(Debug, Clone, Default)]
pub struct RunReport {
    /// Kernel name.
    pub kernel: &'static str,
    /// Execution mode of the run.
    pub mode: &'static str,
    /// Tiles completed.
    pub completed: u64,
    /// Tiles planned.
    pub total: u64,
    /// Successful steal operations.
    pub steals: u64,
    /// Steals whose victim sat on another CCD (should stay the minority
    /// under the CCD-local-first order).
    pub cross_ccd_steals: u64,
    /// Per-worker cancel-observation latencies in ns (empty when the run
    /// was not cancelled).
    pub cancel_latencies_ns: Vec<u64>,
}

impl RunReport {
    /// The p99-ish latency sample (max of the sorted lower 99%; exact max
    /// for fewer than 100 samples). `None` when the run wasn't cancelled.
    #[must_use]
    pub fn cancel_latency_p99_ns(&self) -> Option<u64> {
        if self.cancel_latencies_ns.is_empty() {
            return None;
        }
        let mut v = self.cancel_latencies_ns.clone();
        v.sort_unstable();
        let idx = ((v.len() as f64) * 0.99).ceil() as usize;
        Some(v[idx.saturating_sub(1).min(v.len() - 1)])
    }

    /// Canonical JSON (deterministic field order; latency samples included
    /// verbatim — they are measurements, envelope-class like `wall_ns`).
    #[must_use]
    pub fn to_json(&self) -> String {
        use std::fmt::Write as _;
        let mut s = String::with_capacity(160);
        let _ = write!(
            s,
            "{{\"kernel\":\"{}\",\"mode\":\"{}\",\"completed\":{},\"total\":{},\"steals\":{},\
             \"cross_ccd_steals\":{},\"cancel_latencies_ns\":[",
            self.kernel, self.mode, self.completed, self.total, self.steals, self.cross_ccd_steals
        );
        for (i, l) in self.cancel_latencies_ns.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{l}");
        }
        s.push_str("]}");
        s
    }
}

/// Compute worker `w`'s CCD index under `topo` for `workers` total workers:
/// contiguous blocks, so workers `[k*W/C, (k+1)*W/C)` share CCD `k`.
fn ccd_of_worker(w: usize, workers: usize, topo: CcdTopology) -> usize {
    let ccds = (topo.ccds as usize).max(1);
    (w * ccds) / workers.max(1)
}

/// The steal victim order for worker `w`: same-CCD workers first (ring
/// order after `w`), then the rest (ring order). Pure and deterministic —
/// this function IS what workers use, so verifying it on fixture
/// topologies verifies the runtime behavior.
#[must_use]
pub fn victim_order(w: usize, workers: usize, topo: &CcdTopology) -> Vec<usize> {
    let my_ccd = ccd_of_worker(w, workers, *topo);
    let ring = (1..workers).map(|d| (w + d) % workers);
    let mut same: Vec<usize> = Vec::new();
    let mut other: Vec<usize> = Vec::new();
    for v in ring {
        if ccd_of_worker(v, workers, *topo) == my_ccd {
            same.push(v);
        } else {
            other.push(v);
        }
    }
    same.extend(other);
    same
}

/// Split `0..tiles` into contiguous per-worker ranges proportional to
/// `weights` (largest-remainder rounding; deterministic).
#[must_use]
pub fn weighted_ranges(tiles: u64, weights: &[u32]) -> Vec<core::ops::Range<u64>> {
    let total_w: u64 = weights.iter().map(|&w| u64::from(w.max(1))).sum();
    let mut ranges = Vec::with_capacity(weights.len());
    let mut start = 0u64;
    let mut acc = 0u64;
    for (i, &w) in weights.iter().enumerate() {
        acc += u64::from(w.max(1));
        let end = if i + 1 == weights.len() {
            tiles
        } else {
            (tiles * acc) / total_w
        };
        ranges.push(start..end);
        start = end;
    }
    ranges
}

/// The throughput-lane pool. Workers are scoped per run (spawned at `run`,
/// joined before it returns) so kernel borrows need no `'static` — the
/// persistent-parked-worker optimization is deferred with the lock-free
/// deques (CONTRACT no-claims).
pub struct TilePool {
    config: PoolConfig,
    arenas: fs_alloc::ArenaPool,
    iteration: AtomicU64,
}

impl TilePool {
    /// Build a pool (normalizes the config — see [`PoolConfig`]).
    #[must_use]
    pub fn new(config: PoolConfig) -> Self {
        let mut config = config;
        config.workers = config.workers.max(1);
        config.quantum_weights.resize(config.workers, 1);
        for w in &mut config.quantum_weights {
            *w = (*w).max(1);
        }
        let arenas = fs_alloc::ArenaPool::new(config.arena.clone());
        TilePool {
            config,
            arenas,
            iteration: AtomicU64::new(0),
        }
    }

    /// The arena pool backing per-tile scopes (leak oracle for G4 tests).
    #[must_use]
    pub fn arena_pool(&self) -> &fs_alloc::ArenaPool {
        &self.arenas
    }

    /// Run a kernel to completion with an internal gate (no external
    /// cancellation source).
    ///
    /// # Errors
    /// [`RunError`] on cancellation (kernel-initiated), tile panic, or
    /// executor invariant violation.
    pub fn run<K: TileKernel>(&self, kernel: &K) -> Result<K::Out, RunError> {
        self.run_with_gate(kernel, &CancelGate::new()).0
    }

    /// Run a kernel under an external cancel gate; returns the outcome and
    /// the measured [`RunReport`].
    // One coherent protocol (seed deques -> worker loops -> fold + report);
    // splitting it would scatter the drain/containment invariants the
    // storm suite audits as a unit.
    #[allow(clippy::too_many_lines)]
    pub fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let plan = kernel.tiles();
        let kernel_id = plan.kernel_id();
        let n = plan.tiles;
        let iteration = self.iteration.fetch_add(1, Ordering::Relaxed);
        let workers = self.config.workers.min(n.max(1) as usize).max(1);

        // Fixed-slot reduction storage: one slot per tile, written once.
        let slots: Vec<Mutex<Option<K::Out>>> = (0..n).map(|_| Mutex::new(None)).collect();
        // Per-worker deques seeded with weight-proportional contiguous runs.
        let ranges = weighted_ranges(n, &self.config.quantum_weights[..workers]);
        let deques: Vec<CachePadded<Mutex<VecDeque<u64>>>> = ranges
            .iter()
            .map(|r| CachePadded::new(Mutex::new(r.clone().collect())))
            .collect();
        let victims: Vec<Vec<usize>> = (0..workers)
            .map(|w| victim_order(w, workers, &self.config.topo))
            .collect();

        let steals = AtomicU64::new(0);
        let cross_steals = AtomicU64::new(0);
        let panic_box: Mutex<Option<(u64, String)>> = Mutex::new(None);
        let observed: Vec<CachePadded<AtomicU64>> = (0..workers)
            .map(|_| CachePadded::new(AtomicU64::new(0)))
            .collect();

        std::thread::scope(|s| {
            for w in 0..workers {
                let deques = &deques;
                let slots = &slots;
                let victims = &victims[w];
                let steals = &steals;
                let cross_steals = &cross_steals;
                let panic_box = &panic_box;
                let observed = &observed[w];
                let arenas = &self.arenas;
                let config = &self.config;
                s.spawn(move || {
                    loop {
                        // Tile boundary: the drain point (P7). Record the
                        // observation timestamp once for the histogram.
                        if gate.is_requested() {
                            let _ = observed.get().compare_exchange(
                                0,
                                gate.now_ns().max(1),
                                Ordering::AcqRel,
                                Ordering::Acquire,
                            );
                            break;
                        }
                        // Own deque first (front: preserve locality runs).
                        let mut tile = deques[w].get().lock().expect("deque").pop_front();
                        if tile.is_none() {
                            // Steal HALF from the first non-empty victim,
                            // same-CCD victims first.
                            for &v in victims {
                                let mut vd = deques[v].get().lock().expect("deque");
                                let take = vd.len().div_ceil(2);
                                if take == 0 {
                                    continue;
                                }
                                let split_at = vd.len() - take;
                                let stolen: VecDeque<u64> = vd.split_off(split_at);
                                drop(vd);
                                steals.fetch_add(1, Ordering::Relaxed);
                                if ccd_of_worker(v, workers, config.topo)
                                    != ccd_of_worker(w, workers, config.topo)
                                {
                                    cross_steals.fetch_add(1, Ordering::Relaxed);
                                }
                                let mut mine = deques[w].get().lock().expect("deque");
                                *mine = stolen;
                                tile = mine.pop_front();
                                break;
                            }
                        }
                        let Some(tile) = tile else {
                            break; // every deque empty: run complete
                        };
                        let key = StreamKey {
                            seed: config.seed,
                            kernel_id,
                            tile,
                            iteration,
                        };
                        let outcome = arenas.scope(|arena| {
                            let cx = Cx::new(
                                gate,
                                arena,
                                key,
                                asupersync::types::Budget::INFINITE,
                                config.mode,
                            );
                            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                kernel.run(tile, &cx)
                            }))
                        });
                        match outcome {
                            Ok(ControlFlow::Continue(out)) => {
                                *slots[tile as usize].lock().expect("slot") = Some(out);
                            }
                            Ok(ControlFlow::Break(_cancelled)) => {
                                // Kernel observed the gate (or self-cancelled):
                                // make it global and drain.
                                gate.request();
                            }
                            Err(payload) => {
                                let message = payload
                                    .downcast_ref::<&str>()
                                    .map(ToString::to_string)
                                    .or_else(|| payload.downcast_ref::<String>().cloned())
                                    .unwrap_or_else(|| "non-string panic payload".to_string());
                                let mut pb = panic_box.lock().expect("panic box");
                                if pb.is_none() {
                                    *pb = Some((tile, message));
                                }
                                drop(pb);
                                gate.request();
                            }
                        }
                    }
                });
            }
        });

        let completed = slots
            .iter()
            .filter(|s| s.lock().expect("slot").is_some())
            .count() as u64;
        let requested_at = gate.requested_at_ns();
        let report = RunReport {
            kernel: plan.kernel,
            mode: self.config.mode.name(),
            completed,
            total: n,
            steals: steals.load(Ordering::Relaxed),
            cross_ccd_steals: cross_steals.load(Ordering::Relaxed),
            cancel_latencies_ns: requested_at.map_or_else(Vec::new, |req| {
                observed
                    .iter()
                    .filter_map(|o| match o.get().load(Ordering::Acquire) {
                        0 => None,
                        t => Some(t.saturating_sub(req)),
                    })
                    .collect()
            }),
        };

        if let Some((tile, message)) = panic_box.into_inner().expect("panic box") {
            return (
                Err(RunError::TilePanicked {
                    kernel: plan.kernel,
                    tile,
                    message,
                    completed,
                }),
                report,
            );
        }
        if gate.is_requested() {
            return (
                Err(RunError::Cancelled {
                    kernel: plan.kernel,
                    completed,
                    total: n,
                }),
                report,
            );
        }
        // Fixed-shape fold: ascending tile order, always.
        let mut acc = K::Out::identity();
        for (i, slot) in slots.into_iter().enumerate() {
            match slot.into_inner().expect("slot") {
                Some(out) => acc = acc.merge(out),
                None => {
                    return (
                        Err(RunError::Incomplete {
                            kernel: plan.kernel,
                            tile: i as u64,
                        }),
                        report,
                    );
                }
            }
        }
        (Ok(acc), report)
    }
}

impl fmt::Debug for TilePool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TilePool")
            .field("workers", &self.config.workers)
            .field("mode", &self.config.mode.name())
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::TilePlan;

    struct SumKernel {
        tiles: u64,
    }

    impl TileKernel for SumKernel {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/sum", self.tiles)
        }

        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(crate::Cancelled);
            }
            let buf = cx
                .arena()
                .alloc_slice_fill(fs_alloc::Site::named("test/sum"), 64, tile)
                .expect("arena alloc");
            ControlFlow::Continue(buf.iter().sum::<u64>() / 64 + 1)
        }
    }

    fn pool(workers: usize) -> TilePool {
        TilePool::new(PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, 0x5EED))
    }

    #[test]
    fn completeness_across_worker_and_tile_counts() {
        for workers in [1, 2, 4, 8] {
            for tiles in [0u64, 1, 7, 64, 513] {
                let p = pool(workers);
                let got = p.run(&SumKernel { tiles }).expect("run");
                let want: u64 = (0..tiles).map(|t| t + 1).sum();
                assert_eq!(got, want, "workers={workers} tiles={tiles}");
                assert!(p.arena_pool().stats().quiescent(), "arena leak");
            }
        }
    }

    #[test]
    fn weighted_ranges_are_contiguous_and_proportional() {
        let r = weighted_ranges(100, &[2, 1, 1]);
        assert_eq!(r, vec![0..50, 50..75, 75..100]);
        let r = weighted_ranges(7, &[1, 1]);
        assert_eq!(r, vec![0..3, 3..7]);
        let r = weighted_ranges(0, &[1, 1]);
        assert_eq!(r, vec![0..0, 0..0]);
    }

    #[test]
    fn victim_order_prefers_the_local_ccd() {
        // 8 workers on the Apple fixture (2 CCDs): workers 0..4 on ccd 0.
        let order = victim_order(1, 8, &CcdTopology::APPLE_M_CLASS);
        assert_eq!(order.len(), 7);
        assert_eq!(&order[..3], &[2, 3, 0], "same-CCD ring first");
        assert_eq!(&order[3..], &[4, 5, 6, 7], "cross-CCD after");
    }

    #[test]
    fn kernel_initiated_cancellation_is_a_structured_outcome() {
        struct SelfCancel;
        impl TileKernel for SelfCancel {
            type Out = u64;

            fn tiles(&self) -> TilePlan {
                TilePlan::new("test/self-cancel", 64)
            }

            fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
                if tile == 5 {
                    ControlFlow::Break(crate::Cancelled)
                } else {
                    ControlFlow::Continue(1)
                }
            }
        }
        let p = pool(4);
        let (res, report) = p.run_with_gate(&SelfCancel, &CancelGate::new());
        match res {
            Err(RunError::Cancelled { total: 64, .. }) => {}
            other => panic!("expected Cancelled, got {other:?}"),
        }
        assert_eq!(report.total, 64);
        assert!(
            p.arena_pool().stats().quiescent(),
            "cancelled work must reclaim"
        );
    }

    #[test]
    fn panics_are_contained_with_tile_provenance_and_pool_survives() {
        struct Bomb;
        impl TileKernel for Bomb {
            type Out = u64;

            fn tiles(&self) -> TilePlan {
                TilePlan::new("test/bomb", 32)
            }

            fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
                assert!(tile != 9, "tile 9 exploded");
                ControlFlow::Continue(1)
            }
        }
        let p = pool(4);
        let err = p.run(&Bomb).expect_err("must fail");
        match &err {
            RunError::TilePanicked {
                tile: 9, message, ..
            } => {
                assert!(message.contains("exploded"), "{message}");
            }
            other => panic!("expected TilePanicked{{tile:9}}, got {other:?}"),
        }
        assert!(err.to_string().contains("pool remains usable"));
        // The pool is not poisoned: a healthy kernel still runs.
        let ok = p.run(&SumKernel { tiles: 16 }).expect("pool survives");
        assert_eq!(ok, (0..16).map(|t| t + 1).sum::<u64>());
        assert!(p.arena_pool().stats().quiescent());
    }
}
