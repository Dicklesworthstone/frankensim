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

use crate::cx::{CancelGate, Cx, ExecMode, RunId, StreamKey};
use crate::kernel::TileKernel;
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
    /// OPT-IN OS pinning (fz2.2): worker `w` is pinned to
    /// `pin_groups[ccd_of_worker(w) % len]` — pass the measured L3
    /// groups so each shard's workers stay inside their cache island
    /// (measured on a 5995WX: unpinned threads migrate across CCDs and
    /// lose 8.35x on cache-resident sweeps). Empty = no pinning
    /// (default). ADVISORY and timing-only (P2): pin failures are
    /// ignored by design — results are bit-identical either way, and
    /// the ccd_ab harness verifies the mechanism separately.
    pub pin_groups: Vec<Vec<u32>>,
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
            pin_groups: Vec::new(),
        }
    }

    /// Construct an unpinned deterministic configuration from the host's
    /// topology probe. The probe is a scheduling hint, not a hardware claim;
    /// callers that already hold measured topology should use [`Self::new`].
    #[must_use]
    pub fn for_host(workers: usize, seed: u64) -> Self {
        let probe = fs_substrate::CapabilityProbe::topology_only();
        Self::new(workers, CcdTopology::from_probe(&probe), seed)
    }

    /// Enable CCD pinning from the MEASURED L3 topology where the
    /// platform exposes it (Linux sysfs); a no-op elsewhere — callers
    /// can inspect `pin_groups.is_empty()` to ledger which they got.
    #[must_use]
    pub fn with_measured_pinning(mut self) -> Self {
        let groups = fs_substrate::affinity::measured_l3_groups();
        if let Some(topo) = CcdTopology::from_l3_groups(&groups) {
            self.topo = topo;
            self.pin_groups = groups;
        }
        self
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
    /// The operating system refused to create a scoped worker. Already-started
    /// workers were cancelled and drained before this outcome was returned.
    WorkerSpawn {
        /// Kernel name.
        kernel: &'static str,
        /// Lowest worker index whose creation failed.
        worker: usize,
        /// Operating-system diagnostic.
        message: String,
    },
    /// A user-defined deterministic reduction merge panicked after every tile
    /// had completed. The unwind was contained at the pool boundary.
    ReductionPanicked {
        /// Kernel name.
        kernel: &'static str,
        /// Panic payload's message, when it carried one.
        message: String,
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
            RunError::WorkerSpawn {
                kernel,
                worker,
                message,
            } => write!(
                f,
                "kernel `{kernel}` worker {worker} could not be created: {message}; started workers were cancelled and drained"
            ),
            RunError::ReductionPanicked { kernel, message } => write!(
                f,
                "kernel `{kernel}` deterministic reduction panicked: {message}; the unwind was contained and the pool remains usable"
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

fn push_json_string(out: &mut String, value: &str) {
    use core::fmt::Write as _;

    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Measured facts about one run: steal statistics and the cancel-latency
/// samples (ns between the first cancel request and each worker OBSERVING
/// it at a tile boundary). Measurements only — results never depend on them.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RunReport {
    /// Kernel name.
    pub kernel: &'static str,
    /// Execution mode of the run.
    pub mode: &'static str,
    /// Caller-declared logical run identity used as every tile stream's
    /// iteration component.
    pub declared_run: RunId,
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
    /// Tiles completed per worker (fz2.2): the measured per-class
    /// throughput signal — on heterogeneous cores, slow-class workers
    /// complete measurably fewer tiles under work-stealing.
    pub tiles_by_worker: Vec<u64>,
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
        s.push_str("{\"kernel\":");
        push_json_string(&mut s, self.kernel);
        s.push_str(",\"mode\":");
        push_json_string(&mut s, self.mode);
        let _ = write!(
            s,
            ",\"declared_run\":{},\"completed\":{},\"total\":{},\"steals\":{},\
             \"cross_ccd_steals\":{},\"cancel_latencies_ns\":[",
            self.declared_run.0, self.completed, self.total, self.steals, self.cross_ccd_steals
        );
        for (i, l) in self.cancel_latencies_ns.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{l}");
        }
        s.push_str("],\"tiles_by_worker\":[");
        for (i, completed) in self.tiles_by_worker.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{completed}");
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
/// cumulative weights. Each interior boundary is
/// `floor(tiles * prefix_weight / total_weight)`; the implementation evaluates
/// that ratio exactly without a fixed-width intermediate product.
#[must_use]
pub fn weighted_ranges(tiles: u64, weights: &[u32]) -> Vec<core::ops::Range<u64>> {
    let total_w: u128 = weights.iter().map(|&w| u128::from(w.max(1))).sum();
    let tiles = u128::from(tiles);
    let mut ranges = Vec::with_capacity(weights.len());
    let mut start = 0u64;
    let mut acc = 0u128;
    for (i, &w) in weights.iter().enumerate() {
        acc += u128::from(w.max(1));
        let end = if i + 1 == weights.len() {
            u64::try_from(tiles).expect("u64 tile count widened losslessly")
        } else {
            mul_ratio_floor(
                u64::try_from(tiles).expect("u64 tile count widened losslessly"),
                acc,
                total_w,
            )
        };
        ranges.push(start..end);
        start = end;
    }
    ranges
}

fn mul_ratio_floor(value: u64, numerator: u128, denominator: u128) -> u64 {
    debug_assert!(denominator > 0 && numerator <= denominator);
    // A realizable &[u32] has total weight below 2^96 on 64-bit targets.
    // Maintaining the division remainder keeps every step below 3*denominator
    // instead of forming the potentially 160-bit `value * numerator` product.
    let mut quotient = 0u128;
    let mut remainder = 0u128;
    for bit in (0..u64::BITS).rev() {
        quotient *= 2;
        remainder *= 2;
        if (value >> bit) & 1 == 1 {
            remainder += numerator;
        }
        quotient += remainder / denominator;
        remainder %= denominator;
    }
    u64::try_from(quotient).expect("a ratio no greater than one cannot exceed its u64 multiplicand")
}

/// The throughput-lane pool. Workers are scoped per run (spawned at `run`,
/// joined before it returns) so kernel borrows need no `'static` — the
/// persistent-parked-worker optimization is deferred with the lock-free
/// deques (CONTRACT no-claims).
pub struct TilePool {
    config: PoolConfig,
    arenas: fs_alloc::ArenaPool,
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
        TilePool { config, arenas }
    }

    /// Construct a deterministic, unpinned pool from the host topology probe.
    #[must_use]
    pub fn for_host(workers: usize, seed: u64) -> Self {
        Self::new(PoolConfig::for_host(workers, seed))
    }

    /// Normalized maximum worker budget used by this pool.
    #[must_use]
    pub const fn workers(&self) -> usize {
        self.config.workers
    }

    /// Canonical placement/configuration identity for tune rows and replay
    /// keys. The readable prefix records topology, mode, and pinning intent;
    /// the derive-key BLAKE3 suffix binds normalized weights, arena policy,
    /// the pool's recorded hugepage decision, and exact pin groups without an
    /// unbounded key.
    ///
    /// Pinning is advisory at execution time, but requesting it changes the
    /// timing population and therefore must select a distinct tune key even
    /// on a host where the OS rejects the affinity request.
    #[must_use]
    pub fn placement_identity(&self) -> String {
        let digest = placement_digest(&self.config, self.arenas.hugepage_decision());
        let pinning_intent = if self.config.pin_groups.is_empty() {
            "pin-unrequested"
        } else {
            "ccd-pin-requested"
        };
        format!(
            "fs-exec-tilepool-v2-{pinning_intent}-ccd{}x{}-mode-{}-cfg-{digest}",
            self.config.topo.ccds,
            self.config.topo.cores_per_ccd,
            self.config.mode.name(),
        )
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

    /// Run a kernel under an explicit, caller-ledgered [`RunId`] (bead
    /// wf9.7.1): re-running the SAME kernel with a DIFFERENT logical
    /// run (a new generation, trial, or restart) diverges its streams
    /// by declared identity. `run`/`run_with_gate` are the fixed
    /// `RunId(0)` convenience — bit-identical no matter how much
    /// unrelated or concurrent work the pool has executed.
    pub fn run_declared<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
    ) -> (Result<K::Out, RunError>, RunReport) {
        self.run_inner(kernel, gate, run)
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
        self.run_inner(kernel, gate, RunId::default())
    }

    // One coherent protocol (seed deques -> worker loops -> fold + report);
    // splitting it would scatter the drain/containment invariants the
    // storm suite audits as a unit.
    #[allow(clippy::too_many_lines)]
    fn run_inner<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let plan = kernel.tiles();
        let kernel_id = plan.kernel_id();
        let n = plan.tiles;
        // Stream identity is DECLARED, never scheduled (wf9.7.1): the
        // former pool-global counter made keys depend on unrelated
        // prior runs and on concurrent invocation order.
        let iteration = run.0;
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
        let done_by: Vec<CachePadded<AtomicU64>> = (0..workers)
            .map(|_| CachePadded::new(AtomicU64::new(0)))
            .collect();

        let mut spawn_failure = None;
        std::thread::scope(|s| {
            for w in 0..workers {
                let deques = &deques;
                let slots = &slots;
                let victims = &victims[w];
                let steals = &steals;
                let cross_steals = &cross_steals;
                let panic_box = &panic_box;
                let observed = &observed[w];
                let done_by = &done_by[w];
                let arenas = &self.arenas;
                let config = &self.config;
                let spawned = std::thread::Builder::new().spawn_scoped(s, move || {
                    if !config.pin_groups.is_empty() {
                        let g = ccd_of_worker(w, workers, config.topo) % config.pin_groups.len();
                        // Advisory (see PoolConfig::pin_groups docs).
                        let _ =
                            fs_substrate::os_affinity::pin_current_thread(&config.pin_groups[g]);
                    }
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
                                done_by.get().fetch_add(1, Ordering::Relaxed);
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
                                if pb
                                    .as_ref()
                                    .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
                                {
                                    *pb = Some((tile, message));
                                }
                                drop(pb);
                                gate.request();
                            }
                        }
                    }
                });
                if let Err(error) = spawned {
                    spawn_failure = Some((w, error.to_string()));
                    gate.request();
                    break;
                }
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
            declared_run: run,
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
            tiles_by_worker: done_by
                .iter()
                .map(|c| c.get().load(Ordering::Relaxed))
                .collect(),
        };

        if let Some((worker, message)) = spawn_failure {
            return (
                Err(RunError::WorkerSpawn {
                    kernel: plan.kernel,
                    worker,
                    message,
                }),
                report,
            );
        }
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
        // Fixed-shape fold: the pairwise tree over ascending tile order
        // (shape a pure function of the tile count — plan §5.4).
        let mut outs: Vec<K::Out> = Vec::with_capacity(slots.len());
        for (i, slot) in slots.into_iter().enumerate() {
            match slot.into_inner().expect("slot") {
                Some(out) => outs.push(out),
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
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            crate::reduce::pairwise_fold(outs)
        })) {
            Ok(out) => (Ok(out), report),
            Err(payload) => {
                let message = payload
                    .downcast_ref::<&str>()
                    .map(ToString::to_string)
                    .or_else(|| payload.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "non-string panic payload".to_string());
                (
                    Err(RunError::ReductionPanicked {
                        kernel: plan.kernel,
                        message,
                    }),
                    report,
                )
            }
        }
    }
}

fn placement_digest(config: &PoolConfig, hugepage: &fs_alloc::HugepageDecision) -> String {
    const DOMAIN: &str = "org.frankensim.fs-exec.tilepool-placement.v2";
    let mut payload = Vec::new();
    append_placement_usize(&mut payload, config.workers);
    payload.extend_from_slice(&config.topo.ccds.to_le_bytes());
    payload.extend_from_slice(&config.topo.cores_per_ccd.to_le_bytes());
    payload.push(match config.mode {
        ExecMode::Deterministic => 0,
        ExecMode::Fast => 1,
    });
    append_placement_usize(&mut payload, config.quantum_weights.len());
    for weight in &config.quantum_weights {
        payload.extend_from_slice(&weight.to_le_bytes());
    }
    append_placement_usize(&mut payload, config.arena.chunk_bytes);
    append_placement_usize(&mut payload, config.arena.max_chunk_bytes);
    match config.arena.limit_bytes {
        Some(limit) => {
            payload.push(1);
            append_placement_usize(&mut payload, limit);
        }
        None => payload.push(0),
    }
    append_placement_usize(&mut payload, config.arena.free_list_max_bytes);
    payload.push(match config.arena.hugepage {
        fs_alloc::HugepagePolicy::Auto => 0,
        fs_alloc::HugepagePolicy::Never => 1,
    });
    append_placement_bytes(&mut payload, hugepage.to_json().as_bytes());
    append_placement_usize(&mut payload, config.pin_groups.len());
    for group in &config.pin_groups {
        append_placement_usize(&mut payload, group.len());
        for cpu in group {
            payload.extend_from_slice(&cpu.to_le_bytes());
        }
    }
    fs_blake3::hash_domain(DOMAIN, &payload).to_hex()
}

fn append_placement_usize(payload: &mut Vec<u8>, value: usize) {
    payload.extend_from_slice(
        &u64::try_from(value)
            .expect("TilePool placement dimension exceeds u64")
            .to_le_bytes(),
    );
}

fn append_placement_bytes(payload: &mut Vec<u8>, bytes: &[u8]) {
    append_placement_usize(payload, bytes.len());
    payload.extend_from_slice(bytes);
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
    use crate::kernel::{Reduce, TilePlan};

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

    struct MultiPanicKernel {
        tiles: u64,
        barrier: std::sync::Barrier,
    }

    impl TileKernel for MultiPanicKernel {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/multi-panic", self.tiles)
        }

        fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            self.barrier.wait();
            panic!("simultaneous panic from tile {tile}");
        }
    }

    #[derive(Debug, PartialEq, Eq)]
    struct MergeBomb(u64);

    impl Reduce for MergeBomb {
        fn identity() -> Self {
            Self(0)
        }

        fn merge(self, _other: Self) -> Self {
            panic!("reduction merge exploded")
        }
    }

    struct ReductionPanicKernel;

    impl TileKernel for ReductionPanicKernel {
        type Out = MergeBomb;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/reduction-panic", 2)
        }

        fn run(&self, _tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, MergeBomb> {
            ControlFlow::Continue(MergeBomb(1))
        }
    }

    fn pool(workers: usize) -> TilePool {
        TilePool::new(PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, 0x5EED))
    }

    #[test]
    fn run_report_json_escapes_identity_and_retains_worker_counts() {
        let report = RunReport {
            kernel: "test/\"kernel\\line\n",
            mode: "deterministic",
            declared_run: RunId(7),
            completed: 3,
            total: 4,
            steals: 2,
            cross_ccd_steals: 1,
            cancel_latencies_ns: vec![11, 13],
            tiles_by_worker: vec![2, 1],
        };

        assert_eq!(
            report.to_json(),
            "{\"kernel\":\"test/\\\"kernel\\\\line\\n\",\"mode\":\"deterministic\",\"declared_run\":7,\"completed\":3,\"total\":4,\"steals\":2,\"cross_ccd_steals\":1,\"cancel_latencies_ns\":[11,13],\"tiles_by_worker\":[2,1]}"
        );
    }

    #[test]
    fn simultaneous_panics_report_the_lowest_logical_tile() {
        for workers in [2usize, 4] {
            for _ in 0..16 {
                let kernel = MultiPanicKernel {
                    tiles: workers as u64,
                    barrier: std::sync::Barrier::new(workers),
                };
                let error = pool(workers)
                    .run(&kernel)
                    .expect_err("every in-flight tile panics");
                match error {
                    RunError::TilePanicked { tile, message, .. } => {
                        assert_eq!(tile, 0, "panic provenance must not depend on arrival order");
                        assert_eq!(message, "simultaneous panic from tile 0");
                    }
                    other => panic!("expected TilePanicked, got {other:?}"),
                }
            }
        }
    }

    #[test]
    fn reduction_panics_are_structured_and_the_pool_survives() {
        let pool = pool(2);
        let error = pool
            .run(&ReductionPanicKernel)
            .expect_err("the merge deliberately panics");
        assert_eq!(
            error,
            RunError::ReductionPanicked {
                kernel: "test/reduction-panic",
                message: "reduction merge exploded".to_string(),
            }
        );
        assert_eq!(
            pool.run(&SumKernel { tiles: 17 })
                .expect("reuse after panic"),
            (1_u64..=17).sum::<u64>(),
            "a contained reduction panic must not poison the pool"
        );
        assert!(pool.arena_pool().stats().quiescent());
    }

    #[test]
    fn placement_identity_tracks_the_requested_pinning_intent() {
        let unpinned = pool(0);
        assert_eq!(unpinned.workers(), 1, "worker budgets are normalized");
        let unpinned_identity = unpinned.placement_identity();
        assert!(
            unpinned_identity.starts_with("fs-exec-tilepool-v2-pin-unrequested-ccd"),
            "{unpinned_identity}"
        );

        let mut config = PoolConfig::new(3, CcdTopology::APPLE_M_CLASS, 0x5EED);
        config.pin_groups = vec![vec![9999]];
        let pinned = TilePool::new(config);
        assert_eq!(pinned.workers(), 3);
        let pinned_identity = pinned.placement_identity();
        assert!(
            pinned_identity.starts_with("fs-exec-tilepool-v2-ccd-pin-requested-ccd"),
            "{pinned_identity}"
        );
        assert_ne!(pinned_identity, unpinned_identity);

        let mut weighted = PoolConfig::new(1, CcdTopology::APPLE_M_CLASS, 0x5EED);
        weighted.quantum_weights = vec![2];
        let weighted_identity = TilePool::new(weighted).placement_identity();
        assert_ne!(weighted_identity, unpinned_identity);
        assert!(weighted_identity.len() <= 256);
        assert!(
            weighted_identity
                .bytes()
                .all(|byte| { byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'.' | b'_') })
        );
    }

    #[test]
    fn placement_identity_binds_the_recorded_hugepage_outcome() {
        let pool = pool(1);
        let decision = |outcome| fs_alloc::HugepageDecision {
            policy: fs_alloc::HugepagePolicy::Auto,
            outcome,
            detail: "deterministic fixture detail".to_string(),
        };
        let aligned = decision(fs_alloc::HugepageOutcome::AlignedForThp);
        let unsupported = decision(fs_alloc::HugepageOutcome::UnsupportedPlatform);

        let aligned_digest = placement_digest(&pool.config, &aligned);
        assert_eq!(
            aligned_digest,
            placement_digest(&pool.config, &aligned),
            "the same recorded decision must produce the same placement digest"
        );
        assert_ne!(
            aligned_digest,
            placement_digest(&pool.config, &unsupported),
            "different realized hugepage outcomes must not share tune rows"
        );
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
    fn pinning_is_bit_invariant_and_advisory() {
        // P2: pinning changes timing, never bits — pinned (measured
        // topology where available, garbage groups otherwise) must
        // produce exactly the unpinned result; on targets without
        // pinning support the advisory path is a no-op that still
        // completes the run.
        let tiles = 513u64;
        let want = pool(4).run(&SumKernel { tiles }).expect("unpinned run");
        let measured = TilePool::new(
            PoolConfig::new(4, CcdTopology::APPLE_M_CLASS, 0x5EED).with_measured_pinning(),
        );
        assert_eq!(
            measured.run(&SumKernel { tiles }).expect("pinned run"),
            want
        );
        // Deliberately hostile pin groups (cpu ids that may not exist):
        // advisory pinning must never fail the run or change the bits.
        let mut hostile = PoolConfig::new(4, CcdTopology::APPLE_M_CLASS, 0x5EED);
        hostile.pin_groups = vec![vec![9999], vec![0]];
        assert_eq!(
            TilePool::new(hostile)
                .run(&SumKernel { tiles })
                .expect("hostile-pin run"),
            want
        );
    }

    #[test]
    fn weighted_ranges_are_contiguous_and_proportional() {
        let r = weighted_ranges(100, &[2, 1, 1]);
        assert_eq!(r, vec![0..50, 50..75, 75..100]);
        let r = weighted_ranges(7, &[1, 1]);
        assert_eq!(r, vec![0..3, 3..7]);
        let r = weighted_ranges(0, &[1, 1]);
        assert_eq!(r, vec![0..0, 0..0]);

        let maximal = weighted_ranges(u64::MAX, &[1, 1, u32::MAX, 7]);
        assert_eq!(maximal.first().map(|range| range.start), Some(0));
        assert_eq!(maximal.last().map(|range| range.end), Some(u64::MAX));
        assert!(
            maximal.windows(2).all(|pair| pair[0].end == pair[1].start),
            "maximum-domain partition must have neither gaps nor overlap: {maximal:?}"
        );
        assert!(
            maximal.iter().all(|range| range.start <= range.end),
            "maximum-domain boundaries must be monotonic: {maximal:?}"
        );
        assert_eq!(mul_ratio_floor(u64::MAX, 1, 2), u64::MAX / 2);
        for value in [0, 1, 7, 1024, u64::from(u32::MAX)] {
            for (numerator, denominator) in [(0, 1), (1, 3), (2, 3), (17, 19), (1, 1)] {
                assert_eq!(
                    mul_ratio_floor(value, numerator, denominator),
                    u64::try_from(u128::from(value) * numerator / denominator)
                        .expect("small oracle fits")
                );
            }
        }
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
