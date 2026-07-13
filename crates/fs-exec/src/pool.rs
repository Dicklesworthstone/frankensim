//! The throughput lane: a work-stealing fork-join tile pool (plan §5.2).
//!
//! Semantics first, lock-freedom later: each worker owns one contiguous
//! tile run (`CachePadded<Mutex<TileRun>>`, bead wf9.16.2 — ownership is
//! two u64s, so stealing allocates NOTHING after launch) seeded with
//! contiguous, weight-proportional ranges; an empty worker steals the BACK
//! HALF of a victim's run, visiting same-CCD victims before cross-CCD ones
//! (plan §5.1 consequence 3). The protocol — weighted quanta, CCD-local-first stealing, fixed-slot
//! reductions, drain-on-cancel, panic containment — is the contract; the
//! Chase–Lev lock-free deque is a later optimization gated on roofline
//! evidence (CONTRACT no-claims).
//!
//! Determinism (P2): every tile's output lands in its OWN slot and slots
//! fold in ascending tile order, so results are bit-identical across worker
//! counts and steal schedules by construction. RNG stream keys derive from
//! logical identity only.

use crate::cx::{Budget, CancelGate, Cx, ExecMode, RefusalSink, RunId, StreamKey, TileFailure};
use crate::kernel::TileKernel;
use asupersync::cx::{CpuCx, ScopedCpuError};
use core::fmt;
use core::ops::ControlFlow;
use fs_alloc::CachePadded;
use fs_substrate::affinity::CcdTopology;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

/// Semantic version of the tile-pool placement/tuning identity.
///
/// Version 2 is the already-shipped key and payload. Making the version
/// explicit does not rotate or re-key it.
pub const TILEPOOL_PLACEMENT_IDENTITY_VERSION: u32 = 2;

/// BLAKE3 derive-key domain for the exact v2 placement payload.
pub const TILEPOOL_PLACEMENT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-exec.tilepool-placement.v2";

const TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM: &str = "fs-exec-tilepool-v";

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const TILEPOOL_PLACEMENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-exec:tilepool-placement",
    "version_const=TILEPOOL_PLACEMENT_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-exec.tilepool-placement.v2",
    "domain_const=TILEPOOL_PLACEMENT_IDENTITY_DOMAIN",
    "encoder=TilePool::placement_identity",
    "encoder_helpers=placement_identity_with_schema,placement_digest_with_domain,encode_tilepool_placement,PlacementCounts::from_inputs,append_placement_usize,append_placement_bytes",
    "schema_constants=TILEPOOL_PLACEMENT_IDENTITY_VERSION,TILEPOOL_PLACEMENT_IDENTITY_DOMAIN,TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM,crates/fs-blake3/src/lib.rs#IV,crates/fs-blake3/src/lib.rs#MSG_PERMUTATION,crates/fs-blake3/src/lib.rs#BLOCK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_START,crates/fs-blake3/src/lib.rs#CHUNK_END,crates/fs-blake3/src/lib.rs#PARENT,crates/fs-blake3/src/lib.rs#ROOT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_CONTEXT,crates/fs-blake3/src/lib.rs#DERIVE_KEY_MATERIAL,crates/fs-blake3/src/lib.rs#MAX_DEPTH",
    "schema_functions=crates/fs-exec/src/cx.rs#ExecMode::name,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-blake3/src/lib.rs#ContentHash::to_hex,crates/fs-blake3/src/lib.rs#g,crates/fs-blake3/src/lib.rs#round,crates/fs-blake3/src/lib.rs#permute,crates/fs-blake3/src/lib.rs#compress,crates/fs-blake3/src/lib.rs#words_from_block,crates/fs-blake3/src/lib.rs#first_8_words,crates/fs-blake3/src/lib.rs#Output::chaining_value,crates/fs-blake3/src/lib.rs#Output::root_hash,crates/fs-blake3/src/lib.rs#parent_output,crates/fs-blake3/src/lib.rs#ChunkState::new,crates/fs-blake3/src/lib.rs#ChunkState::len,crates/fs-blake3/src/lib.rs#ChunkState::start_flag,crates/fs-blake3/src/lib.rs#ChunkState::update,crates/fs-blake3/src/lib.rs#ChunkState::output,crates/fs-blake3/src/lib.rs#Blake3::new_internal,crates/fs-blake3/src/lib.rs#Blake3::push_stack,crates/fs-blake3/src/lib.rs#Blake3::pop_stack,crates/fs-blake3/src/lib.rs#Blake3::add_chunk_chaining_value,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize",
    "schema_dependencies=fs-alloc:hugepage-decision",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=PoolConfig,TilePoolPlacementTopologyFields,TilePoolPlacementArenaFields,TilePoolPlacementHugepageFields,PlacementCounts",
    "source_fields=PoolConfig.workers:semantic,PoolConfig.topo:derived:expanded-into-exact-topology-fields,PoolConfig.quantum_weights:semantic,PoolConfig.seed:nonsemantic:logical-stream-identity-not-placement,PoolConfig.mode:semantic,PoolConfig.arena:derived:expanded-into-exact-arena-fields,PoolConfig.pin_groups:semantic,TilePoolPlacementTopologyFields.ccds:semantic,TilePoolPlacementTopologyFields.cores_per_ccd:semantic,TilePoolPlacementArenaFields.chunk_bytes:semantic,TilePoolPlacementArenaFields.max_chunk_bytes:semantic,TilePoolPlacementArenaFields.limit_bytes:semantic,TilePoolPlacementArenaFields.free_list_max_bytes:semantic,TilePoolPlacementArenaFields.hugepage:semantic,TilePoolPlacementHugepageFields.policy:semantic,TilePoolPlacementHugepageFields.outcome:semantic,TilePoolPlacementHugepageFields.detail:semantic,PlacementCounts.workers:derived:exact-count-of-normalized-workers,PlacementCounts.quantum_weights:derived:exact-count-of-normalized-quantum-weights,PlacementCounts.hugepage_json_bytes:derived:exact-byte-count-of-canonical-hugepage-json,PlacementCounts.pin_groups:derived:exact-count-of-requested-pin-groups,PlacementCounts.pin_cpus:derived:ordered-exact-counts-of-cpus-per-requested-pin-group",
    "source_bindings=PoolConfig.workers>workers,PoolConfig.quantum_weights>quantum-weight-count+quantum-weights-in-order,PoolConfig.mode>mode-tag,PoolConfig.pin_groups>pinning-intent+pin-group-count+pin-cpu-counts+pin-cpu-ids-in-order,TilePoolPlacementTopologyFields.ccds>topology-ccds,TilePoolPlacementTopologyFields.cores_per_ccd>topology-cores-per-ccd,TilePoolPlacementArenaFields.chunk_bytes>arena-chunk-bytes,TilePoolPlacementArenaFields.max_chunk_bytes>arena-max-chunk-bytes,TilePoolPlacementArenaFields.limit_bytes>arena-limit-presence+arena-limit-bytes,TilePoolPlacementArenaFields.free_list_max_bytes>arena-free-list-max-bytes,TilePoolPlacementArenaFields.hugepage>arena-hugepage-policy-tag,TilePoolPlacementHugepageFields.policy>hugepage-decision-policy,TilePoolPlacementHugepageFields.outcome>hugepage-decision-outcome,TilePoolPlacementHugepageFields.detail>hugepage-json-byte-count+hugepage-decision-detail-json",
    "external_semantic_fields=digest-domain,identity-prefix-stem,identity-version",
    "semantic_fields=digest-domain,identity-prefix-stem,identity-version,workers,topology-ccds,topology-cores-per-ccd,mode-tag,quantum-weight-count,quantum-weights-in-order,arena-chunk-bytes,arena-max-chunk-bytes,arena-limit-presence,arena-limit-bytes,arena-free-list-max-bytes,arena-hugepage-policy-tag,hugepage-decision-policy,hugepage-decision-outcome,hugepage-json-byte-count,hugepage-decision-detail-json,pinning-intent,pin-group-count,pin-cpu-counts,pin-cpu-ids-in-order",
    "excluded_fields=pin-success:observed-timing-fact-not-requested-placement",
    "consumers=TilePool::placement_identity,TilePool::admit_retained_placement_identity,fs-exec::tuner,replay-and-tune-rows",
    "mutations=digest-domain:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,identity-prefix-stem:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,identity-version:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,workers:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,topology-ccds:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,topology-cores-per-ccd:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,mode-tag:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,quantum-weight-count:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,quantum-weights-in-order:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,arena-chunk-bytes:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,arena-max-chunk-bytes:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,arena-limit-presence:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,arena-limit-bytes:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,arena-free-list-max-bytes:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,arena-hugepage-policy-tag:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,hugepage-decision-policy:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,hugepage-decision-outcome:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,hugepage-json-byte-count:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,hugepage-decision-detail-json:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,pinning-intent:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,pin-group-count:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,pin-cpu-counts:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently,pin-cpu-ids-in-order:crates/fs-exec/src/pool.rs#tilepool_placement_identity_fields_move_independently",
    "nonsemantic_mutations=PoolConfig.seed:crates/fs-exec/src/pool.rs#tilepool_placement_seed_is_nonsemantic,pin-success:crates/fs-exec/src/pool.rs#pinning_is_bit_invariant_and_advisory",
    "field_guard=classify_tilepool_placement_identity_fields",
    "transport_guard=TilePool::admit_retained_placement_identity",
    "version_guard=crates/fs-exec/src/pool.rs#tilepool_placement_identity_versions_fail_closed",
    "coupling_surface=fs-exec:tilepool-placement",
];

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

#[derive(Debug, Clone, Copy)]
struct TilePoolPlacementTopologyFields {
    ccds: u32,
    cores_per_ccd: u32,
}

#[derive(Debug, Clone, Copy)]
struct TilePoolPlacementArenaFields {
    chunk_bytes: usize,
    max_chunk_bytes: usize,
    limit_bytes: Option<usize>,
    free_list_max_bytes: usize,
    hugepage: fs_alloc::HugepagePolicy,
}

#[derive(Debug, Clone, Copy)]
struct TilePoolPlacementHugepageFields<'a> {
    policy: fs_alloc::HugepagePolicy,
    outcome: fs_alloc::HugepageOutcome,
    detail: &'a str,
}

#[allow(dead_code)]
fn classify_tilepool_placement_identity_fields(
    config: &PoolConfig,
    topology_fields: TilePoolPlacementTopologyFields,
    arena_fields: &TilePoolPlacementArenaFields,
    hugepage_fields: &TilePoolPlacementHugepageFields<'_>,
    counts: &PlacementCounts,
    hugepage_decision: &fs_alloc::HugepageDecision,
) {
    let PoolConfig {
        workers,
        topo,
        quantum_weights,
        seed,
        mode,
        arena,
        pin_groups,
    } = config;
    let CcdTopology {
        ccds,
        cores_per_ccd,
    } = topo;
    let fs_alloc::ArenaConfig {
        chunk_bytes,
        max_chunk_bytes,
        limit_bytes,
        free_list_max_bytes,
        hugepage,
    } = arena;
    let TilePoolPlacementTopologyFields {
        ccds: identity_ccds,
        cores_per_ccd: identity_cores_per_ccd,
    } = topology_fields;
    let TilePoolPlacementArenaFields {
        chunk_bytes: identity_chunk_bytes,
        max_chunk_bytes: identity_max_chunk_bytes,
        limit_bytes: identity_limit_bytes,
        free_list_max_bytes: identity_free_list_max_bytes,
        hugepage: identity_hugepage,
    } = arena_fields;
    let TilePoolPlacementHugepageFields {
        policy,
        outcome,
        detail,
    } = hugepage_fields;
    let PlacementCounts {
        workers: counted_workers,
        quantum_weights: counted_quantum_weights,
        hugepage_json_bytes,
        pin_groups: counted_pin_groups,
        pin_cpus,
    } = counts;
    let fs_alloc::HugepageDecision {
        policy: recorded_policy,
        outcome: recorded_outcome,
        detail: recorded_detail,
    } = hugepage_decision;
    let _ = (
        workers,
        ccds,
        cores_per_ccd,
        quantum_weights,
        seed,
        mode,
        chunk_bytes,
        max_chunk_bytes,
        limit_bytes,
        free_list_max_bytes,
        hugepage,
        pin_groups,
        identity_ccds,
        identity_cores_per_ccd,
        identity_chunk_bytes,
        identity_max_chunk_bytes,
        identity_limit_bytes,
        identity_free_list_max_bytes,
        identity_hugepage,
        policy,
        outcome,
        detail,
        counted_workers,
        counted_quantum_weights,
        hugepage_json_bytes,
        counted_pin_groups,
        pin_cpus,
        recorded_policy,
        recorded_outcome,
        recorded_detail,
    );
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
    /// A tile returned a typed refusal; siblings were cancelled and drained.
    TileFailed {
        /// Kernel name.
        kernel: &'static str,
        /// Lowest logical tile that reported a refusal before drain completed.
        tile: u64,
        /// Typed refusal suitable for upstream policy and ledger handling.
        failure: TileFailure,
        /// Tiles that completed despite the refusal.
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
    /// The operation memory lease refused the pool's root metadata BEFORE
    /// worker launch (bead wf9.16); nothing ran and no root metadata was
    /// allocated.
    /// Mid-run per-tile refusals surface as [`RunError::TileFailed`] with an
    /// allocation failure instead.
    MemoryRefused {
        /// Kernel name.
        kernel: &'static str,
        /// Component that was refused.
        what: &'static str,
        /// Bytes the component requested.
        requested_bytes: u64,
        /// Lease bytes already in use at refusal time.
        used_bytes: u64,
        /// The lease limit in force.
        limit_bytes: u64,
    },
    /// A root-metadata dimension or byte total cannot be represented on this
    /// target. Refused before lease mutation, allocation, or worker launch.
    MemoryPlanOverflow {
        /// Kernel name.
        kernel: &'static str,
        /// First root component whose checked sizing overflowed.
        what: &'static str,
    },
    /// The global allocator refused fallible root-metadata reservation before
    /// worker launch. The operation-lease charge is rolled back on return.
    MemoryAllocationRefused {
        /// Kernel name.
        kernel: &'static str,
        /// Root component whose backing allocation was refused.
        what: &'static str,
        /// Logical bytes requested for that component.
        requested_bytes: u64,
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
            RunError::TileFailed {
                kernel,
                tile,
                failure,
                completed,
            } => write!(
                f,
                "kernel `{kernel}` tile {tile} refused: {failure} ({completed} sibling tiles \
                 completed; siblings were cancelled and drained, the pool remains usable)"
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
            RunError::MemoryRefused {
                kernel,
                what,
                requested_bytes,
                used_bytes,
                limit_bytes,
            } => write!(
                f,
                "kernel `{kernel}` refused before launch: `{what}` needs {requested_bytes} B \
                 with {used_bytes} B of the {limit_bytes} B operation memory lease already in \
                 use; nothing ran and no root metadata was allocated"
            ),
            RunError::MemoryPlanOverflow { kernel, what } => write!(
                f,
                "kernel `{kernel}` refused before launch: checked sizing for root component \
                 `{what}` exceeds this target's representable memory domain; reduce the tile \
                 or worker count"
            ),
            RunError::MemoryAllocationRefused {
                kernel,
                what,
                requested_bytes,
            } => write!(
                f,
                "kernel `{kernel}` refused before launch: the global allocator could not reserve \
                 {requested_bytes} B for root component `{what}`; the lease charge was rolled back"
            ),
        }
    }
}

impl core::error::Error for RunError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::TileFailed { failure, .. } => Some(failure),
            _ => None,
        }
    }
}

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

fn prelaunch_report(kernel: &'static str, mode: &'static str, run: RunId, total: u64) -> RunReport {
    RunReport {
        kernel,
        mode,
        declared_run: run,
        completed: 0,
        total,
        steals: 0,
        cross_ccd_steals: 0,
        cancel_latencies_ns: Vec::new(),
        tiles_by_worker: Vec::new(),
    }
}

/// Compute worker `w`'s CCD index under `topo` for `workers` total workers:
/// contiguous blocks, so workers `[k*W/C, (k+1)*W/C)` share CCD `k`.
fn ccd_of_worker(w: usize, workers: usize, topo: CcdTopology) -> usize {
    let ccds = (topo.ccds as usize).max(1);
    (w * ccds) / workers.max(1)
}

/// Checked conservative logical bytes for one run's tracked root metadata
/// (bead wf9.16): slots, deque headers and initial tile-id entries, range-plan entries,
/// victim-table headers/final entries/construction temporaries, per-worker
/// cache-padded atomics, retained pairwise-fold buffers, and report vectors.
/// Thread stacks, allocator bookkeeping, and heap owned by arbitrary kernel
/// outputs are explicit no-claims (CONTRACT): this is an enforceable tracked
/// envelope, not a full-process byte census.
fn root_metadata_bytes<K: TileKernel>(n: u64, workers: usize) -> Result<u64, &'static str> {
    let workers = u64::try_from(workers).map_err(|_| "worker-count")?;
    let slot = size_of::<Mutex<Option<K::Out>>>() as u64;
    let deque_header = size_of::<CachePadded<Mutex<TileRun>>>() as u64;
    let range = size_of::<core::ops::Range<u64>>() as u64;
    let victim_header = size_of::<Vec<usize>>() as u64;
    let report_value = size_of::<u64>() as u64;
    let atomic = size_of::<CachePadded<AtomicU64>>() as u64;
    let out = size_of::<K::Out>() as u64;
    let victim_entries = root_mul(workers, workers.saturating_sub(1), "victim-table-entries")?;
    // victim_order builds one final vector while its `other` partition is
    // still live. The checked constructor reserves workers-1 entries for that
    // temporary, so the peak is final tables plus one extra partition.
    let victim_temporary_entries = workers.saturating_sub(1);
    // pairwise_fold recursively split_offs right halves while parent buffers
    // remain allocated: n + floor tree rights = at most 2n-1 elements.
    let fold_elements = if n == 0 {
        0
    } else {
        root_mul(n, 2, "fold-buffer-elements")?
            .checked_sub(1)
            .ok_or("fold-buffer-elements")?
    };

    let components = [
        ("slot-table", root_mul(n, slot, "slot-table")?),
        // No deque-entries component (bead wf9.16.2): worker ownership is a
        // TileRun of two u64s inside the header, never per-tile storage.
        (
            "deque-headers",
            root_mul(workers, deque_header, "deque-headers")?,
        ),
        ("range-plans", root_mul(workers, range, "range-plans")?),
        (
            "victim-table-headers",
            root_mul(workers, victim_header, "victim-table-headers")?,
        ),
        (
            "victim-table-entries",
            root_mul(
                victim_entries,
                size_of::<usize>() as u64,
                "victim-table-entries",
            )?,
        ),
        (
            "victim-order-temporary",
            root_mul(
                victim_temporary_entries,
                size_of::<usize>() as u64,
                "victim-order-temporary",
            )?,
        ),
        (
            "worker-counters",
            root_mul(
                root_mul(workers, 2, "worker-counters")?,
                atomic,
                "worker-counters",
            )?,
        ),
        (
            "fold-buffers",
            root_mul(fold_elements, out, "fold-buffers")?,
        ),
        (
            "report-vectors",
            root_mul(
                root_mul(workers, 2, "report-vectors")?,
                report_value,
                "report-vectors",
            )?,
        ),
    ];
    components
        .into_iter()
        .try_fold(0_u64, |total, (what, bytes)| {
            total.checked_add(bytes).ok_or(what)
        })
}

fn root_mul(a: u64, b: u64, what: &'static str) -> Result<u64, &'static str> {
    a.checked_mul(b).ok_or(what)
}

fn allocation_bytes<T>(capacity: usize) -> u64 {
    u64::try_from(capacity)
        .ok()
        .and_then(|count| count.checked_mul(size_of::<T>() as u64))
        .unwrap_or(u64::MAX)
}

fn try_reserve_root_vec<T>(
    values: &mut Vec<T>,
    capacity: usize,
    kernel: &'static str,
    what: &'static str,
) -> Result<(), RunError> {
    values
        .try_reserve_exact(capacity)
        .map_err(|_| RunError::MemoryAllocationRefused {
            kernel,
            what,
            requested_bytes: allocation_bytes::<T>(capacity),
        })
}

struct RunRoot<T> {
    slots: Vec<Mutex<Option<T>>>,
    _ranges: Vec<core::ops::Range<u64>>,
    deques: Vec<CachePadded<Mutex<TileRun>>>,
    victims: Vec<Vec<usize>>,
    observed: Vec<CachePadded<AtomicU64>>,
    done_by: Vec<CachePadded<AtomicU64>>,
    cancel_latencies_ns: Vec<u64>,
    tiles_by_worker: Vec<u64>,
    outs: Vec<T>,
}

fn allocate_run_root<K: TileKernel>(
    n: u64,
    n_usize: usize,
    workers: usize,
    weights: &[u32],
    topo: CcdTopology,
    kernel: &'static str,
) -> Result<RunRoot<K::Out>, RunError> {
    let mut slots = Vec::new();
    try_reserve_root_vec(&mut slots, n_usize, kernel, "slot-table")?;
    for _ in 0..n_usize {
        slots.push(Mutex::new(None));
    }

    let active_weights = weights.get(..workers).ok_or(RunError::MemoryPlanOverflow {
        kernel,
        what: "worker-weights",
    })?;
    let ranges = try_weighted_ranges(n, active_weights, kernel)?;
    let mut deques = Vec::new();
    try_reserve_root_vec(&mut deques, workers, kernel, "deque-headers")?;
    for range in &ranges {
        // One contiguous run per worker (bead wf9.16.2): ownership is two
        // u64s, so there is no per-tile entry storage to reserve and the
        // steal protocol allocates nothing after launch.
        deques.push(CachePadded::new(Mutex::new(TileRun::from_range(range))));
    }

    let mut victims = Vec::new();
    try_reserve_root_vec(&mut victims, workers, kernel, "victim-table-headers")?;
    for worker in 0..workers {
        victims.push(try_victim_order(worker, workers, topo, kernel)?);
    }

    let mut observed = Vec::new();
    let mut done_by = Vec::new();
    try_reserve_root_vec(&mut observed, workers, kernel, "worker-counters")?;
    try_reserve_root_vec(&mut done_by, workers, kernel, "worker-counters")?;
    for _ in 0..workers {
        observed.push(CachePadded::new(AtomicU64::new(0)));
        done_by.push(CachePadded::new(AtomicU64::new(0)));
    }

    let mut cancel_latencies_ns = Vec::new();
    let mut tiles_by_worker = Vec::new();
    try_reserve_root_vec(&mut cancel_latencies_ns, workers, kernel, "report-vectors")?;
    try_reserve_root_vec(&mut tiles_by_worker, workers, kernel, "report-vectors")?;

    let mut outs = Vec::new();
    try_reserve_root_vec(&mut outs, n_usize, kernel, "fold-buffers")?;

    Ok(RunRoot {
        slots,
        _ranges: ranges,
        deques,
        victims,
        observed,
        done_by,
        cancel_latencies_ns,
        tiles_by_worker,
        outs,
    })
}

/// One worker's owned work: a contiguous ascending run of logical tile ids
/// (bead wf9.16.2). The pool's stealing protocol maintains a structural
/// invariant that makes this exact: deques are seeded with contiguous
/// weighted ranges, workers only ever pop the FRONT, and a (necessarily
/// empty) thief wholesale-adopts the victim's BACK half — which is itself
/// contiguous. Ownership transfer is therefore pure `Copy` arithmetic on
/// two `u64`s: ZERO allocation after launch, and the peak storage is
/// exactly one cache-padded slot per worker, admitted pre-launch.
#[derive(Debug, Clone, Copy)]
struct TileRun {
    /// Next tile to execute (front).
    next: u64,
    /// One past the last owned tile.
    end: u64,
}

impl TileRun {
    fn from_range(range: &core::ops::Range<u64>) -> Self {
        TileRun {
            next: range.start,
            end: range.end,
        }
    }

    fn len(self) -> u64 {
        self.end.saturating_sub(self.next)
    }

    fn pop_front(&mut self) -> Option<u64> {
        if self.next < self.end {
            let tile = self.next;
            self.next += 1;
            Some(tile)
        } else {
            None
        }
    }

    /// Split off the BACK `ceil(len/2)` tiles — the exact `take`
    /// arithmetic of the previous `VecDeque::split_off` protocol, so the
    /// tile→worker transfer is preserved verbatim, not just semantically.
    fn steal_back_half(&mut self) -> Option<TileRun> {
        let take = self.len().div_ceil(2);
        if take == 0 {
            return None;
        }
        let stolen = TileRun {
            next: self.end - take,
            end: self.end,
        };
        self.end -= take;
        Some(stolen)
    }
}

/// The steal victim order for worker `w`: same-CCD workers first (ring
/// order after `w`), then the rest (ring order). Pure and deterministic —
/// this function IS what workers use, so verifying it on fixture
/// topologies verifies the runtime behavior.
#[must_use]
pub fn victim_order(w: usize, workers: usize, topo: &CcdTopology) -> Vec<usize> {
    let capacity = workers.saturating_sub(1);
    let mut same = Vec::with_capacity(capacity);
    let mut other = Vec::with_capacity(capacity);
    partition_victims(w, workers, *topo, &mut same, &mut other);
    same.extend(other);
    same
}

fn try_victim_order(
    w: usize,
    workers: usize,
    topo: CcdTopology,
    kernel: &'static str,
) -> Result<Vec<usize>, RunError> {
    let capacity = workers.saturating_sub(1);
    let mut same = Vec::new();
    let mut other = Vec::new();
    try_reserve_root_vec(&mut same, capacity, kernel, "victim-table-entries")?;
    try_reserve_root_vec(&mut other, capacity, kernel, "victim-order-temporary")?;
    partition_victims(w, workers, topo, &mut same, &mut other);
    same.extend(other);
    Ok(same)
}

fn partition_victims(
    w: usize,
    workers: usize,
    topo: CcdTopology,
    same: &mut Vec<usize>,
    other: &mut Vec<usize>,
) {
    let my_ccd = ccd_of_worker(w, workers, topo);
    for d in 1..workers {
        let v = (w + d) % workers;
        if ccd_of_worker(v, workers, topo) == my_ccd {
            same.push(v);
        } else {
            other.push(v);
        }
    }
}

/// Split `0..tiles` into contiguous per-worker ranges proportional to
/// cumulative weights. Each interior boundary is
/// `floor(tiles * prefix_weight / total_weight)`; the implementation evaluates
/// that ratio exactly without a fixed-width intermediate product.
#[must_use]
pub fn weighted_ranges(tiles: u64, weights: &[u32]) -> Vec<core::ops::Range<u64>> {
    let mut ranges = Vec::with_capacity(weights.len());
    fill_weighted_ranges(tiles, weights, &mut ranges);
    ranges
}

fn try_weighted_ranges(
    tiles: u64,
    weights: &[u32],
    kernel: &'static str,
) -> Result<Vec<core::ops::Range<u64>>, RunError> {
    let mut ranges = Vec::new();
    try_reserve_root_vec(&mut ranges, weights.len(), kernel, "range-plans")?;
    fill_weighted_ranges(tiles, weights, &mut ranges);
    Ok(ranges)
}

fn fill_weighted_ranges(tiles: u64, weights: &[u32], ranges: &mut Vec<core::ops::Range<u64>>) {
    let total_w: u128 = weights.iter().map(|&w| u128::from(w.max(1))).sum();
    let tiles = u128::from(tiles);
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

/// Worker-lifetime strategy for one run: spawn into a fresh std scope,
/// spawn as scoped-CPU children of the calling task (bead lx0e), or
/// dispatch to an already-parked crew (bead tkr7). All three drive
/// [`worker_loop`], so results are bitwise-identical across strategies
/// by construction (P2).
enum Launch<'a, Caps: 'static> {
    OwnScope,
    TaskScope(&'a asupersync::Cx<Caps>),
    Crew(&'a crate::crew::Crew<Caps>),
}

// Manual impls: every variant is a reference (or unit), so Launch is Copy
// regardless of whether Caps itself is — the derive would demand
// `Caps: Copy` spuriously.
impl<Caps: 'static> Clone for Launch<'_, Caps> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<Caps: 'static> Copy for Launch<'_, Caps> {}

/// Caps stand-in for launches that carry no task context.
type NoTask = asupersync::cx::cap::All;

/// Everything one worker's loop touches, bundled so the launch
/// harnesses — `std::thread::scope`, asupersync's `Cx::scoped_cpu`
/// (bead lx0e), and the parked crew (bead tkr7) — drive the IDENTICAL
/// protocol: seed deques, steal-half, drain-on-cancel, per-tile panic
/// containment. One loop, three worker-lifetime strategies.
struct WorkerCtx<'a, K: TileKernel> {
    kernel: &'a K,
    kernel_id: u64,
    iteration: u64,
    workers: usize,
    budget: Budget,
    gate: &'a CancelGate,
    lease: &'a fs_alloc::OperationMemoryLease,
    arenas: &'a fs_alloc::ArenaPool,
    config: &'a PoolConfig,
    deques: &'a [CachePadded<Mutex<TileRun>>],
    slots: &'a [Mutex<Option<K::Out>>],
    victims: &'a [Vec<usize>],
    observed: &'a [CachePadded<AtomicU64>],
    done_by: &'a [CachePadded<AtomicU64>],
    steals: &'a AtomicU64,
    cross_steals: &'a AtomicU64,
    panic_box: &'a Mutex<Option<(u64, String)>>,
    refusal_sink: &'a RefusalSink,
}

/// The worker protocol, shared verbatim by both launch harnesses. When
/// `task_cx` is present (the asupersync lane), the CALLING task's
/// cancellation and budget bound the run: every tile boundary checkpoints
/// the task context, and a failed checkpoint converts into a gate request
/// so the pool's normal drain protocol — including its cancel-latency
/// histogram — applies unchanged (P7: one drain semantics, two signals).
fn worker_loop<Caps, K: TileKernel>(
    ctx: &WorkerCtx<'_, K>,
    w: usize,
    task_cx: Option<&CpuCx<Caps>>,
) {
    if !ctx.config.pin_groups.is_empty() {
        let g = ccd_of_worker(w, ctx.workers, ctx.config.topo) % ctx.config.pin_groups.len();
        // Advisory (see PoolConfig::pin_groups docs).
        let _ = fs_substrate::os_affinity::pin_current_thread(&ctx.config.pin_groups[g]);
    }
    loop {
        // Tile boundary: the drain point (P7). Bridge the calling task's
        // cancellation/budget first (charged once per boundary), then
        // record the observation timestamp once for the histogram.
        if let Some(task_cx) = task_cx
            && !ctx.gate.is_requested()
            && task_cx.checkpoint().is_err()
        {
            ctx.gate.request();
        }
        if ctx.gate.is_requested() {
            let _ = ctx.observed[w].get().compare_exchange(
                0,
                ctx.gate.now_ns().max(1),
                Ordering::AcqRel,
                Ordering::Acquire,
            );
            break;
        }
        // Own deque first (front: preserve locality runs).
        let mut tile = ctx.deques[w].get().lock().expect("deque").pop_front();
        if tile.is_none() {
            // Steal HALF from the first non-empty victim,
            // same-CCD victims first.
            for &v in &ctx.victims[w] {
                let mut vd = ctx.deques[v].get().lock().expect("deque");
                let Some(stolen) = vd.steal_back_half() else {
                    continue;
                };
                drop(vd);
                ctx.steals.fetch_add(1, Ordering::Relaxed);
                if ccd_of_worker(v, ctx.workers, ctx.config.topo)
                    != ccd_of_worker(w, ctx.workers, ctx.config.topo)
                {
                    ctx.cross_steals.fetch_add(1, Ordering::Relaxed);
                }
                let mut mine = ctx.deques[w].get().lock().expect("deque");
                *mine = stolen;
                tile = mine.pop_front();
                break;
            }
        }
        let Some(tile) = tile else {
            break; // every deque empty: run complete
        };
        let key = StreamKey {
            seed: ctx.config.seed,
            kernel_id: ctx.kernel_id,
            tile,
            iteration: ctx.iteration,
        };
        // Every tile arena charges the shared operation
        // lease while its chunks are held (bead wf9.16).
        let outcome = ctx.arenas.scope_leased(ctx.lease, |arena| {
            let cx = Cx::new_with_refusal_sink(
                ctx.gate,
                arena,
                key,
                ctx.budget,
                ctx.config.mode,
                ctx.refusal_sink,
                ctx.lease,
            );
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| ctx.kernel.run(tile, &cx)))
        });
        match outcome {
            Ok(ControlFlow::Continue(out)) => {
                *ctx.slots[tile as usize].lock().expect("slot") = Some(out);
                ctx.done_by[w].get().fetch_add(1, Ordering::Relaxed);
            }
            Ok(ControlFlow::Break(_cancelled)) => {
                // Kernel observed the gate (or self-cancelled):
                // make it global and drain.
                ctx.gate.request();
            }
            Err(payload) => {
                let message = payload
                    .downcast_ref::<&str>()
                    .map(ToString::to_string)
                    .or_else(|| payload.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "non-string panic payload".to_string());
                let mut pb = ctx.panic_box.lock().expect("panic box");
                if pb
                    .as_ref()
                    .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
                {
                    *pb = Some((tile, message));
                }
                drop(pb);
                ctx.gate.request();
            }
        }
    }
}

/// The throughput-lane pool. Workers are scoped per run (spawned at `run`,
/// joined before it returns) so kernel borrows need no `'static`; callers
/// with many small runs park a crew once instead
/// ([`TilePool::with_parked_crew`], bead tkr7) — the lock-free deques
/// remain deferred (CONTRACT no-claims). Three launch harnesses share one
/// worker protocol: the std lane (`run`/`run_declared*`), the asupersync
/// lane ([`TilePool::run_scoped`], bead lx0e) where workers are scoped
/// CPU children of the calling task via `Cx::scoped_cpu`, and the parked
/// lane ([`ParkedTilePool`]) where runs dispatch to workers already
/// parked inside their owner's scope.
pub struct TilePool {
    config: PoolConfig,
    arenas: fs_alloc::ArenaPool,
}

impl TilePool {
    /// Current producer version for placement and tune-row identities.
    pub const PLACEMENT_IDENTITY_VERSION: u32 = TILEPOOL_PLACEMENT_IDENTITY_VERSION;

    /// Current BLAKE3 derive-key domain for placement identities.
    pub const PLACEMENT_IDENTITY_DOMAIN: &str = TILEPOOL_PLACEMENT_IDENTITY_DOMAIN;

    /// Normalized worker count — preflight sizing for callers that
    /// budget per-worker scratch (bead wf9.15).
    #[must_use]
    pub const fn workers(&self) -> usize {
        self.config.workers
    }

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
        let pinning_intent = if self.config.pin_groups.is_empty() {
            "pin-unrequested"
        } else {
            "ccd-pin-requested"
        };
        let hugepage = self.arenas.hugepage_decision();
        let counts = PlacementCounts::from_inputs(&self.config, hugepage);
        placement_identity_with_schema(
            &self.config,
            hugepage,
            TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM,
            TILEPOOL_PLACEMENT_IDENTITY_VERSION,
            TILEPOOL_PLACEMENT_IDENTITY_DOMAIN,
            pinning_intent,
            &counts,
        )
    }

    /// Fail-closed admission for a retained placement/tuning identity.
    ///
    /// Only the current explicit producer version and the exact identity
    /// recomputed from this normalized pool are admitted. A stale/future
    /// version or any byte mismatch is refused; callers must migrate old
    /// tune rows deliberately rather than treating them as current.
    ///
    /// # Errors
    /// Returns a stable refusal message when the producer version is not v2
    /// or the retained identity differs from the current normalized pool.
    pub fn admit_retained_placement_identity(
        &self,
        declared_version: u32,
        retained_identity: &str,
    ) -> Result<(), &'static str> {
        if declared_version != TILEPOOL_PLACEMENT_IDENTITY_VERSION {
            return Err("tile-pool placement identity version is unsupported");
        }
        if retained_identity != self.placement_identity() {
            return Err("tile-pool placement identity does not match normalized configuration");
        }
        Ok(())
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
        self.run_inner(
            kernel,
            gate,
            run,
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::<NoTask>::OwnScope,
        )
    }

    /// Run a kernel under explicit logical identity and asupersync budget.
    /// Every tile receives the exact same budget slice in its [`Cx`]; kernels
    /// remain responsible for consuming or interpreting its quota dimensions.
    pub fn run_declared_budgeted<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
    ) -> (Result<K::Out, RunError>, RunReport) {
        self.run_inner(
            kernel,
            gate,
            run,
            budget,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::<NoTask>::OwnScope,
        )
    }

    /// [`TilePool::run_declared_budgeted`] under a shared operation memory
    /// lease (bead wf9.16): root metadata is reserved fallibly BEFORE worker
    /// launch, and every tile arena's chunks charge the lease while held.
    /// The caller keeps the lease and reads `lease.receipt()` for the
    /// canonical accounting of that admission trace. Thread stacks, allocator
    /// bookkeeping, and arbitrary heap owned directly by kernels or their
    /// outputs are explicitly not claimed.
    /// The output bound is the sealed admission contract (bead wf9.16.1):
    /// `K::Out` must be [`crate::LeaseAdmittedOut`], so a heap-bearing
    /// custom output whose payload is invisible to `size_of` FAILS TO
    /// COMPILE here. List-shaped outputs use [`crate::Concat`] over
    /// [`fs_alloc::LeasedVec`]; legacy unleased entries stay unconstrained.
    pub fn run_declared_leased_budgeted<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> (Result<K::Out, RunError>, RunReport)
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_inner(kernel, gate, run, budget, lease, Launch::<NoTask>::OwnScope)
    }

    /// Run a kernel under a LIVE asupersync task context (bead lx0e): the
    /// workers are scoped CPU children of the calling task via
    /// `Cx::scoped_cpu`, so task cancellation and budget exhaustion drain
    /// the pool exactly like a gate request, and the scope tree stays
    /// honest (P7) — the pool cannot outlive or leak past the calling
    /// task, which remains blocked here until every worker joins.
    ///
    /// `budget` remains the per-tile slice stamped into each tile's
    /// [`Cx`] (fs-exec vocabulary), independent of the CALLING task's
    /// asupersync budget, which bounds the run itself: each worker
    /// checkpoints `task_cx` at every tile boundary (charging poll quota
    /// per boundary), and cancellation or exhaustion converts to a drain.
    ///
    /// # Errors
    /// Everything [`TilePool::run_declared_leased_budgeted`] can return,
    /// plus [`RunError::Cancelled`] when the calling task is cancelled or
    /// budget-exhausted at entry (nothing runs), mid-run (drain), or at
    /// exit (completed results are refused fail-closed: a cancelled task
    /// must not admit work finished under it).
    ///
    /// # Panics
    /// Pool-invariant panics (a worker dying OUTSIDE per-tile
    /// containment) propagate, exactly like the std-scope lane. OS-level
    /// worker-spawn failure also panics in this lane (upstream
    /// `scoped_cpu` spawns through `std::thread::Scope::spawn`), unlike
    /// the std lane's structured [`RunError::WorkerSpawn`].
    pub fn run_scoped<Caps, K: TileKernel>(
        &self,
        task_cx: &asupersync::Cx<Caps>,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> (Result<K::Out, RunError>, RunReport)
    where
        Caps: Send + Sync + 'static,
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_inner(kernel, gate, run, budget, lease, Launch::TaskScope(task_cx))
    }

    /// Park a crew of exactly [`TilePool::workers`] workers as scoped CPU
    /// children of the CALLING task (bead tkr7) and run `f` with a
    /// [`ParkedTilePool`] whose runs dispatch to those parked workers
    /// instead of spawning — the per-run spawn/join cost that collapses
    /// small-kernel attainment (measured on N-D FFT axis passes, bead
    /// 27d3) drops to a condvar wake/sleep.
    ///
    /// The scope tree stays honest (P7): the crew lives inside this
    /// task's `Cx::scoped_cpu` scope, the calling task blocks here until
    /// every worker joins, and a shutdown guard releases parked workers
    /// on BOTH normal return and unwind of `f`, so the join can never
    /// hang. Task cancellation and budget exhaustion drain RUNNING
    /// kernels at tile boundaries through each worker's own scoped-CPU
    /// context, exactly like [`TilePool::run_scoped`].
    ///
    /// # Errors
    /// [`CrewScopeError::Cancelled`] when the calling task is cancelled
    /// or budget-exhausted at the crew scope's entry (nothing runs, `f`
    /// is never called) or exit (fail closed: a cancelled task must not
    /// admit results computed under it).
    ///
    /// # Panics
    /// Pool-invariant panics (a parked worker dying outside job
    /// containment) propagate, with spawned-lane parity.
    pub fn with_parked_crew<Caps, R, F>(
        &self,
        task_cx: &asupersync::Cx<Caps>,
        f: F,
    ) -> Result<R, CrewScopeError>
    where
        Caps: Send + Sync + 'static,
        F: FnOnce(&ParkedTilePool<'_, Caps>) -> R,
    {
        let crew = crate::crew::Crew::new(self.config.workers);
        match task_cx.scoped_cpu(self.config.workers, |scope| {
            let _shutdown = crate::crew::CrewShutdown(&crew);
            for w in 0..crew.workers() {
                let crew = &crew;
                scope
                    .spawn(move |cpu| crew.park_loop(w, Some(cpu)))
                    .expect("crew spawns exactly its own cap of workers");
            }
            f(&ParkedTilePool {
                pool: self,
                crew: &crew,
            })
        }) {
            Ok(out) => Ok(out),
            Err(ScopedCpuError::Cancelled(_)) => Err(CrewScopeError::Cancelled),
            // Parked workers contain job panics inside the crew; a panic
            // escaping park bookkeeping is a pool invariant failure.
            Err(ScopedCpuError::ChildPanicked { child, message }) => std::panic::panic_any(
                format!("parked crew worker {child} panicked outside job containment: {message}"),
            ),
            Err(ScopedCpuError::WorkerCapExceeded { cap }) => std::panic::panic_any(format!(
                "parked crew launch exceeded its own worker cap {cap}"
            )),
        }
    }

    /// [`TilePool::with_parked_crew`] for callers with NO ambient
    /// asupersync task (perf lanes, tests, batch tools): the crew parks
    /// inside a plain `std::thread::scope`, which is itself scope-sound
    /// (function-blocks-until-join). Cancellation still flows through
    /// each run's [`CancelGate`]; there is simply no task to bridge.
    pub fn with_parked_crew_local<R, F>(&self, f: F) -> R
    where
        F: FnOnce(&ParkedTilePool<'_, NoTask>) -> R,
    {
        let crew: crate::crew::Crew<NoTask> = crate::crew::Crew::new(self.config.workers);
        std::thread::scope(|s| {
            let _shutdown = crate::crew::CrewShutdown(&crew);
            for w in 0..crew.workers() {
                let crew = &crew;
                s.spawn(move || crew.park_loop(w, None));
            }
            f(&ParkedTilePool {
                pool: self,
                crew: &crew,
            })
        })
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
        self.run_inner(
            kernel,
            gate,
            RunId::default(),
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::<NoTask>::OwnScope,
        )
    }

    // One coherent protocol (seed deques -> worker loops -> fold + report);
    // splitting it would scatter the drain/containment invariants the
    // storm suite audits as a unit.
    #[allow(clippy::too_many_lines)]
    fn run_inner<Caps, K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
        launch: Launch<'_, Caps>,
    ) -> (Result<K::Out, RunError>, RunReport)
    where
        Caps: Send + Sync + 'static,
    {
        let plan = kernel.tiles();
        let kernel_id = plan.kernel_id();
        let n = plan.tiles;
        // Stream identity is DECLARED, never scheduled (wf9.7.1): the
        // former pool-global counter made keys depend on unrelated
        // prior runs and on concurrent invocation order.
        let iteration = run.0;
        let Ok(n_usize) = usize::try_from(n) else {
            return (
                Err(RunError::MemoryPlanOverflow {
                    kernel: plan.kernel,
                    what: "tile-count",
                }),
                prelaunch_report(plan.kernel, self.config.mode.name(), run, n),
            );
        };
        let workers = self.config.workers.min(n_usize.max(1)).max(1);

        // Root metadata is reserved fallibly BEFORE any of it is allocated
        // and BEFORE worker launch (bead wf9.16). The charge covers slots,
        // deque headers/initial entries, range plans, victim tables plus their
        // construction temporary, per-worker atomics, the retained
        // pairwise-fold buffers, and report vectors. Thread stacks, allocator
        // bookkeeping, and arbitrary kernel/output-owned heap are explicit
        // no-claims. The guard holds until the run returns (including unwinds).
        let root_bytes = match root_metadata_bytes::<K>(n, workers) {
            Ok(bytes) => bytes,
            Err(what) => {
                return (
                    Err(RunError::MemoryPlanOverflow {
                        kernel: plan.kernel,
                        what,
                    }),
                    prelaunch_report(plan.kernel, self.config.mode.name(), run, n),
                );
            }
        };
        let _root_charge = match lease.reserve("tilepool-root-metadata", root_bytes) {
            Ok(charge) => charge,
            Err(refusal) => {
                return (
                    Err(RunError::MemoryRefused {
                        kernel: plan.kernel,
                        what: refusal.what,
                        requested_bytes: refusal.requested_bytes,
                        used_bytes: refusal.used_bytes,
                        limit_bytes: refusal.limit_bytes,
                    }),
                    prelaunch_report(plan.kernel, self.config.mode.name(), run, n),
                );
            }
        };

        let RunRoot {
            slots,
            _ranges,
            deques,
            victims,
            observed,
            done_by,
            mut cancel_latencies_ns,
            mut tiles_by_worker,
            mut outs,
        } = match allocate_run_root::<K>(
            n,
            n_usize,
            workers,
            &self.config.quantum_weights,
            self.config.topo,
            plan.kernel,
        ) {
            Ok(root) => root,
            Err(error) => {
                return (
                    Err(error),
                    prelaunch_report(plan.kernel, self.config.mode.name(), run, n),
                );
            }
        };

        let steals = AtomicU64::new(0);
        let cross_steals = AtomicU64::new(0);
        let panic_box: Mutex<Option<(u64, String)>> = Mutex::new(None);
        let refusal_sink = RefusalSink::default();

        let ctx = WorkerCtx {
            kernel,
            kernel_id,
            iteration,
            workers,
            budget,
            gate,
            lease,
            arenas: &self.arenas,
            config: &self.config,
            deques: &deques,
            slots: &slots,
            victims: &victims,
            observed: &observed,
            done_by: &done_by,
            steals: &steals,
            cross_steals: &cross_steals,
            panic_box: &panic_box,
            refusal_sink: &refusal_sink,
        };
        let mut spawn_failure = None;
        let mut scope_refusal = None;
        match launch {
            Launch::OwnScope => {
                std::thread::scope(|s| {
                    for w in 0..workers {
                        let ctx = &ctx;
                        let spawned = std::thread::Builder::new()
                            .spawn_scoped(s, move || worker_loop::<Caps, K>(ctx, w, None));
                        if let Err(error) = spawned {
                            spawn_failure = Some((w, error.to_string()));
                            gate.request();
                            break;
                        }
                    }
                });
            }
            Launch::Crew(crew) => {
                // The parked lane (bead tkr7): no spawns at all — the job
                // is dispatched to workers already parked inside their
                // owner's scope, and dispatch blocks until every one of
                // them reports done, so run-local borrows in `ctx` outlive
                // every use (the crew capsule's latch argument). Task
                // cancellation/budget bridging rides each worker's own
                // park-time CpuCx, exactly like the scoped lane. Crew
                // workers beyond this run's normalized count no-op: the
                // run-local tables are sized to `ctx.workers`.
                let job = |w: usize, cpu: Option<&CpuCx<Caps>>| {
                    if w < ctx.workers {
                        worker_loop(&ctx, w, cpu);
                    }
                };
                if let Some((worker, message)) = crew.dispatch(&job) {
                    // Parity with the spawned lanes: a panic escaping a
                    // worker (kernel panics are contained per tile) is a
                    // pool invariant failure — propagate.
                    std::panic::panic_any(format!(
                        "tile-pool worker {worker} panicked outside tile containment: {message}"
                    ));
                }
            }
            Launch::TaskScope(cx) => {
                // The asupersync lane (bead lx0e): workers are scoped CPU
                // children of the CALLING task via `Cx::scoped_cpu`, which
                // blocks here until every child joins — the scope tree is
                // honest by construction (the region cannot close under the
                // workers because its task is inside this call).
                match cx.scoped_cpu(workers, |scope| {
                    for w in 0..workers {
                        let ctx = &ctx;
                        if let Err(error) = scope.spawn(move |cpu| worker_loop(ctx, w, Some(cpu))) {
                            // Structurally unreachable (exactly `workers`
                            // spawns under a cap of `workers`); kept as the
                            // same failure class as an OS spawn refusal.
                            spawn_failure = Some((w, error.to_string()));
                            gate.request();
                            break;
                        }
                    }
                }) {
                    Ok(()) => {}
                    Err(refusal) => scope_refusal = Some(refusal),
                }
            }
        }
        if let Some(refusal) = scope_refusal {
            match refusal {
                // Entry refusal (nothing ran) or exit checkpoint (the
                // calling task was cancelled or exhausted its budget, even
                // if every tile completed first): fail closed as a
                // cancellation — the drain is already complete because the
                // scope joins all workers before returning.
                ScopedCpuError::Cancelled(_) => gate.request(),
                // Parity with the std lane: kernel panics are contained per
                // tile INSIDE the worker loop, so a panic escaping a worker
                // is a pool invariant failure — propagate, exactly as
                // `std::thread::scope` would have.
                ScopedCpuError::ChildPanicked { child, message } => std::panic::panic_any(format!(
                    "tile-pool worker {child} panicked outside tile containment: {message}"
                )),
                // Unreachable by construction; refuse loudly rather than
                // misreport.
                ScopedCpuError::WorkerCapExceeded { cap } => std::panic::panic_any(format!(
                    "tile-pool scoped launch exceeded its own worker cap {cap}"
                )),
            }
        }

        let completed = slots
            .iter()
            .filter(|s| s.lock().expect("slot").is_some())
            .count() as u64;
        let requested_at = gate.requested_at_ns();
        if let Some(requested_at) = requested_at {
            for observed_at in &observed {
                match observed_at.get().load(Ordering::Acquire) {
                    0 => {}
                    observed_at => {
                        cancel_latencies_ns.push(observed_at.saturating_sub(requested_at));
                    }
                }
            }
        }
        for completed_by_worker in &done_by {
            tiles_by_worker.push(completed_by_worker.get().load(Ordering::Relaxed));
        }
        let report = RunReport {
            kernel: plan.kernel,
            mode: self.config.mode.name(),
            declared_run: run,
            completed,
            total: n,
            steals: steals.load(Ordering::Relaxed),
            cross_ccd_steals: cross_steals.load(Ordering::Relaxed),
            cancel_latencies_ns,
            tiles_by_worker,
        };

        // Stable failure-class precedence preserves legacy panic containment
        // while keeping typed refusals distinct from ordinary cancellation.
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
        if let Some((tile, failure)) = refusal_sink.take() {
            return (
                Err(RunError::TileFailed {
                    kernel: plan.kernel,
                    tile,
                    failure,
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

/// Structured refusal from [`TilePool::with_parked_crew`]: the calling
/// task was cancelled or budget-exhausted at the crew scope's entry or
/// exit checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrewScopeError {
    /// Entry refusal (`f` never ran) or exit fail-closed (a cancelled
    /// task must not admit results computed under it).
    Cancelled,
}

impl fmt::Display for CrewScopeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CrewScopeError::Cancelled => write!(
                f,
                "parked-crew scope refused: the calling task was cancelled or exhausted its \
                 budget at the scope boundary"
            ),
        }
    }
}

impl core::error::Error for CrewScopeError {}

/// A [`TilePool`] view whose runs dispatch to an already-parked worker
/// crew (bead tkr7) instead of spawning per run. Created by
/// [`TilePool::with_parked_crew`] / [`TilePool::with_parked_crew_local`];
/// same run surface and the SAME worker protocol, so results are
/// bitwise-identical to the spawned lanes by construction (P2) — only
/// the worker-lifetime strategy differs.
pub struct ParkedTilePool<'a, Caps: 'static> {
    pool: &'a TilePool,
    crew: &'a crate::crew::Crew<Caps>,
}

impl<Caps: Send + Sync + 'static> ParkedTilePool<'_, Caps> {
    /// Normalized worker count (the crew's size — the same value the
    /// spawned lanes normalize to).
    #[must_use]
    pub fn workers(&self) -> usize {
        self.pool.workers()
    }

    /// The arena pool backing per-tile scopes (leak oracle for G4 tests).
    #[must_use]
    pub fn arena_pool(&self) -> &fs_alloc::ArenaPool {
        self.pool.arena_pool()
    }

    /// [`TilePool::run`] on the parked crew.
    ///
    /// # Errors
    /// As [`TilePool::run`].
    pub fn run<K: TileKernel>(&self, kernel: &K) -> Result<K::Out, RunError> {
        self.run_with_gate(kernel, &CancelGate::new()).0
    }

    /// [`TilePool::run_with_gate`] on the parked crew.
    pub fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        self.pool.run_inner(
            kernel,
            gate,
            RunId::default(),
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::Crew(self.crew),
        )
    }

    /// [`TilePool::run_declared`] on the parked crew.
    pub fn run_declared<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
    ) -> (Result<K::Out, RunError>, RunReport) {
        self.pool.run_inner(
            kernel,
            gate,
            run,
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::Crew(self.crew),
        )
    }

    /// [`TilePool::run_declared_leased_budgeted`] on the parked crew.
    pub fn run_declared_leased_budgeted<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> (Result<K::Out, RunError>, RunReport) {
        self.pool
            .run_inner(kernel, gate, run, budget, lease, Launch::Crew(self.crew))
    }
}

impl<Caps: Send + Sync + 'static> crate::kernel::KernelRunner for ParkedTilePool<'_, Caps> {
    fn workers(&self) -> usize {
        ParkedTilePool::workers(self)
    }

    fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        ParkedTilePool::run_with_gate(self, kernel, gate)
    }
}

impl<Caps: 'static> fmt::Debug for ParkedTilePool<'_, Caps> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ParkedTilePool")
            .field("workers", &self.crew.workers())
            .finish_non_exhaustive()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PlacementCounts {
    workers: usize,
    quantum_weights: usize,
    hugepage_json_bytes: usize,
    pin_groups: usize,
    pin_cpus: Vec<usize>,
}

impl PlacementCounts {
    fn from_inputs(config: &PoolConfig, hugepage: &fs_alloc::HugepageDecision) -> Self {
        Self {
            workers: config.workers,
            quantum_weights: config.quantum_weights.len(),
            hugepage_json_bytes: hugepage.to_json().len(),
            pin_groups: config.pin_groups.len(),
            pin_cpus: config.pin_groups.iter().map(Vec::len).collect(),
        }
    }
}

fn placement_identity_with_schema(
    config: &PoolConfig,
    hugepage: &fs_alloc::HugepageDecision,
    prefix_stem: &str,
    version: u32,
    domain: &str,
    pinning_intent: &str,
    counts: &PlacementCounts,
) -> String {
    let digest = placement_digest_with_domain(config, hugepage, domain, counts);
    format!(
        "{prefix_stem}{version}-{pinning_intent}-ccd{}x{}-mode-{}-cfg-{digest}",
        config.topo.ccds,
        config.topo.cores_per_ccd,
        config.mode.name(),
    )
}

#[cfg(test)]
fn placement_digest(config: &PoolConfig, hugepage: &fs_alloc::HugepageDecision) -> String {
    let counts = PlacementCounts::from_inputs(config, hugepage);
    placement_digest_with_domain(
        config,
        hugepage,
        TILEPOOL_PLACEMENT_IDENTITY_DOMAIN,
        &counts,
    )
}

fn placement_digest_with_domain(
    config: &PoolConfig,
    hugepage: &fs_alloc::HugepageDecision,
    domain: &str,
    counts: &PlacementCounts,
) -> String {
    let payload = encode_tilepool_placement(config, hugepage, counts);
    fs_blake3::hash_domain(domain, &payload).to_hex()
}

fn encode_tilepool_placement(
    config: &PoolConfig,
    hugepage: &fs_alloc::HugepageDecision,
    counts: &PlacementCounts,
) -> Vec<u8> {
    let topology = TilePoolPlacementTopologyFields {
        ccds: config.topo.ccds,
        cores_per_ccd: config.topo.cores_per_ccd,
    };
    let arena = TilePoolPlacementArenaFields {
        chunk_bytes: config.arena.chunk_bytes,
        max_chunk_bytes: config.arena.max_chunk_bytes,
        limit_bytes: config.arena.limit_bytes,
        free_list_max_bytes: config.arena.free_list_max_bytes,
        hugepage: config.arena.hugepage,
    };
    let mut payload = Vec::new();
    append_placement_usize(&mut payload, counts.workers);
    payload.extend_from_slice(&topology.ccds.to_le_bytes());
    payload.extend_from_slice(&topology.cores_per_ccd.to_le_bytes());
    payload.push(match config.mode {
        ExecMode::Deterministic => 0,
        ExecMode::Fast => 1,
    });
    append_placement_usize(&mut payload, counts.quantum_weights);
    for weight in &config.quantum_weights {
        payload.extend_from_slice(&weight.to_le_bytes());
    }
    append_placement_usize(&mut payload, arena.chunk_bytes);
    append_placement_usize(&mut payload, arena.max_chunk_bytes);
    match arena.limit_bytes {
        Some(limit) => {
            payload.push(1);
            append_placement_usize(&mut payload, limit);
        }
        None => payload.push(0),
    }
    append_placement_usize(&mut payload, arena.free_list_max_bytes);
    payload.push(match arena.hugepage {
        fs_alloc::HugepagePolicy::Auto => 0,
        fs_alloc::HugepagePolicy::Never => 1,
    });
    let hugepage_json = hugepage.to_json();
    append_placement_bytes(
        &mut payload,
        hugepage_json.as_bytes(),
        counts.hugepage_json_bytes,
    );
    append_placement_usize(&mut payload, counts.pin_groups);
    for (index, group) in config.pin_groups.iter().enumerate() {
        append_placement_usize(&mut payload, counts.pin_cpus[index]);
        for cpu in group {
            payload.extend_from_slice(&cpu.to_le_bytes());
        }
    }
    payload
}

fn append_placement_usize(payload: &mut Vec<u8>, value: usize) {
    payload.extend_from_slice(
        &u64::try_from(value)
            .expect("TilePool placement dimension exceeds u64")
            .to_le_bytes(),
    );
}

fn append_placement_bytes(payload: &mut Vec<u8>, bytes: &[u8], declared_len: usize) {
    append_placement_usize(payload, declared_len);
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

    struct BudgetProbe {
        tiles: u64,
    }

    impl TileKernel for BudgetProbe {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/budget-probe", self.tiles)
        }

        fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            ControlFlow::Continue(cx.budget().remaining_cost().unwrap_or(u64::MAX))
        }
    }

    struct SimultaneousAllocationRefusal {
        tiles: u64,
        barrier: std::sync::Barrier,
    }

    impl TileKernel for SimultaneousAllocationRefusal {
        type Out = ();

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/allocation-refusal", self.tiles)
        }

        fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, ()> {
            self.barrier.wait();
            match cx
                .arena()
                .alloc_slice_fill(fs_alloc::Site::named("test/refusal"), 1, 0_u8)
            {
                Ok(_) => ControlFlow::Continue(()),
                Err(error) => ControlFlow::Break(cx.refuse(TileFailure::Allocation(error))),
            }
        }
    }

    struct MixedPanicAndRefusal {
        barrier: std::sync::Barrier,
    }

    impl TileKernel for MixedPanicAndRefusal {
        type Out = ();

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/mixed-panic-refusal", 2)
        }

        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, ()> {
            self.barrier.wait();
            assert!(tile != 1, "mixed failure panic");
            match cx
                .arena()
                .alloc_slice_fill(fs_alloc::Site::named("test/mixed-refusal"), 1, 0_u8)
            {
                Ok(_) => ControlFlow::Continue(()),
                Err(error) => ControlFlow::Break(cx.refuse(TileFailure::Allocation(error))),
            }
        }
    }

    struct NoAllocation;

    impl TileKernel for NoAllocation {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/no-allocation", 1)
        }

        fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            if cx.checkpoint().is_err() {
                ControlFlow::Break(crate::Cancelled)
            } else {
                ControlFlow::Continue(1)
            }
        }
    }

    struct UnrepresentablePlan;

    impl TileKernel for UnrepresentablePlan {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/unrepresentable-root", u64::MAX)
        }

        fn run(&self, _tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            panic!("an unrepresentable plan must be refused before launch")
        }
    }

    fn pool(workers: usize) -> TilePool {
        TilePool::new(PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, 0x5EED))
    }

    /// Run `f` inside a REAL asupersync root task (the canonical shape from
    /// `latency.rs`): the scoped lane demands a live task context, and using
    /// the runtime — not a synthetic Cx — is the P7-honest harness.
    fn in_task<R, F>(f: F) -> R
    where
        F: FnOnce(asupersync::Cx) -> R + Send + 'static,
        R: Send + 'static,
    {
        let lane = crate::LatencyLane::new(1).expect("lane");
        let root = lane.runtime().handle().spawn(async move {
            let cx = asupersync::Cx::current().expect("task Cx");
            f(cx)
        });
        lane.block_on(root)
    }

    fn run_scoped_simple<K: TileKernel>(
        p: &TilePool,
        cx: &asupersync::Cx,
        kernel: &K,
    ) -> (Result<K::Out, RunError>, RunReport)
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        p.run_scoped(
            cx,
            kernel,
            &CancelGate::new(),
            RunId(0),
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
        )
    }

    /// A kernel that cancels the CALLING asupersync task from inside a tile
    /// — the G4 mid-run cancellation storm's trigger.
    struct CancelTaskAt {
        tiles: u64,
        at: u64,
        task: asupersync::Cx,
    }

    impl TileKernel for CancelTaskAt {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/cancel-task-at", self.tiles)
        }

        fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            if tile == self.at {
                self.task.set_cancel_requested(true);
            }
            ControlFlow::Continue(1)
        }
    }

    struct PanicAt {
        tiles: u64,
        at: u64,
    }

    impl TileKernel for PanicAt {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/panic-at", self.tiles)
        }

        fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            assert!(tile != self.at, "scoped containment probe");
            ControlFlow::Continue(1)
        }
    }

    /// G5 (lx0e): the asupersync lane is deterministic across reruns and
    /// bitwise-identical to the std lane — one worker protocol, two scopes —
    /// with every tile accounted to a worker and arenas quiescent after.
    #[test]
    fn scoped_lane_is_deterministic_and_matches_std_lane_bitwise() {
        let expected = pool(4).run(&SumKernel { tiles: 257 }).expect("std lane");
        let (first, second, report) = in_task(|cx| {
            let p = pool(4);
            let (first, report) = run_scoped_simple(&p, &cx, &SumKernel { tiles: 257 });
            let (second, _) = run_scoped_simple(&p, &cx, &SumKernel { tiles: 257 });
            assert!(
                p.arena_pool().stats().quiescent(),
                "arenas quiescent after scoped runs"
            );
            (
                first.expect("scoped lane"),
                second.expect("scoped rerun"),
                report,
            )
        });
        assert_eq!(first, expected, "scoped lane bitwise-matches the std lane");
        assert_eq!(second, expected, "scoped lane deterministic across reruns");
        assert_eq!(report.completed, 257);
        assert_eq!(
            report.tiles_by_worker.iter().sum::<u64>(),
            257,
            "every tile accounted to a worker"
        );
    }

    /// G4 (lx0e): a pre-cancelled task refuses at ENTRY — no worker spawns,
    /// nothing runs — and the pool stays usable once the task is live again.
    #[test]
    fn scoped_lane_refuses_pre_cancelled_task_at_entry_and_pool_survives() {
        in_task(|cx| {
            let p = pool(2);
            cx.set_cancel_requested(true);
            let (out, report) = run_scoped_simple(&p, &cx, &SumKernel { tiles: 64 });
            match out {
                Err(RunError::Cancelled {
                    completed: 0,
                    total: 64,
                    ..
                }) => {}
                other => panic!("entry refusal must be Cancelled with zero work: {other:?}"),
            }
            assert_eq!(report.completed, 0, "no tile ran under a cancelled task");
            assert!(p.arena_pool().stats().quiescent(), "nothing to leak");
            cx.set_cancel_requested(false);
            let (out, _) = run_scoped_simple(&p, &cx, &SumKernel { tiles: 64 });
            let rerun = out.expect("pool usable after an entry refusal");
            assert_eq!(
                rerun,
                pool(2).run(&SumKernel { tiles: 64 }).expect("std lane"),
                "post-refusal rerun bitwise-matches the std lane"
            );
        });
    }

    /// G4 (lx0e): cancelling the calling TASK mid-run converts into the
    /// pool's own drain protocol — workers stop at tile boundaries, the run
    /// fails closed as Cancelled, arenas quiesce, and the pool survives.
    #[test]
    fn scoped_lane_drains_on_mid_run_task_cancel_and_fails_closed() {
        in_task(|cx| {
            let p = pool(4);
            let kernel = CancelTaskAt {
                tiles: 16_384,
                at: 0,
                task: cx.clone(),
            };
            let (out, report) = run_scoped_simple(&p, &cx, &kernel);
            match out {
                Err(RunError::Cancelled {
                    completed, total, ..
                }) => {
                    assert_eq!(total, 16_384);
                    assert!(
                        completed < total,
                        "a tile-0 task cancel must drain before completion \
                         (completed {completed})"
                    );
                }
                other => panic!("task cancel must fail closed as Cancelled: {other:?}"),
            }
            assert!(
                report.completed < 16_384,
                "report agrees the drain preempted completion"
            );
            assert!(
                p.arena_pool().stats().quiescent(),
                "drained workers leaked no arena chunks"
            );
            cx.set_cancel_requested(false);
            let (out, _) = run_scoped_simple(&p, &cx, &SumKernel { tiles: 64 });
            out.expect("pool usable after a mid-run task cancel");
        });
    }

    /// G4 (lx0e): per-tile panic containment holds unchanged in the scoped
    /// lane — the panic is localized to its tile, siblings drain, the scope
    /// joins, and the pool survives.
    #[test]
    fn scoped_lane_contains_tile_panics_and_survives() {
        in_task(|cx| {
            let p = pool(4);
            let (out, _report) = run_scoped_simple(&p, &cx, &PanicAt { tiles: 512, at: 7 });
            match out {
                Err(RunError::TilePanicked { tile, message, .. }) => {
                    assert_eq!(tile, 7, "the panic is localized to its tile");
                    assert!(
                        message.contains("scoped containment probe"),
                        "payload survives: {message}"
                    );
                }
                other => panic!("a tile panic must surface as TilePanicked: {other:?}"),
            }
            assert!(
                p.arena_pool().stats().quiescent(),
                "containment leaked no arena chunks"
            );
            let (out, _) = run_scoped_simple(&p, &cx, &SumKernel { tiles: 64 });
            out.expect("pool usable after tile panic containment");
        });
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
    fn declared_budget_reaches_every_tile_without_changing_legacy_wrappers() {
        for workers in [1, 4] {
            let pool = pool(workers);
            let gate = CancelGate::new();
            let budget = Budget::new().with_cost_quota(65_536);
            let probe = BudgetProbe {
                tiles: workers as u64,
            };
            let (result, report) = pool.run_declared_budgeted(&probe, &gate, RunId(17), budget);
            assert_eq!(result.expect("budgeted probe"), 65_536 * workers as u64);
            assert_eq!(report.declared_run, RunId(17));
            assert_eq!(
                pool.run(&probe).expect("legacy probe"),
                u64::MAX.wrapping_mul(workers as u64)
            );
        }
    }

    #[test]
    fn root_metadata_plan_counts_fold_and_victim_construction_peaks() {
        let n = 9_u64;
        let workers = 4_u64;
        let slot = size_of::<Mutex<Option<u64>>>() as u64;
        let deque_header = size_of::<CachePadded<Mutex<TileRun>>>() as u64;
        let range = size_of::<core::ops::Range<u64>>() as u64;
        let victim_header = size_of::<Vec<usize>>() as u64;
        let atomic = size_of::<CachePadded<AtomicU64>>() as u64;
        // No per-tile deque-entries term (bead wf9.16.2).
        let expected = n * slot
            + workers * deque_header
            + workers * range
            + workers * victim_header
            + (workers * (workers - 1) + (workers - 1)) * size_of::<usize>() as u64
            + workers * 2 * atomic
            + (2 * n - 1) * size_of::<u64>() as u64
            + workers * 2 * size_of::<u64>() as u64;
        assert_eq!(
            root_metadata_bytes::<SumKernel>(n, workers as usize),
            Ok(expected)
        );
    }

    #[test]
    fn unrepresentable_root_plan_is_refused_before_lease_or_launch() {
        let pool = pool(2);
        let lease = fs_alloc::OperationMemoryLease::unbounded();
        let (result, report) = pool.run_declared_leased_budgeted(
            &UnrepresentablePlan,
            &CancelGate::new(),
            RunId(29),
            Budget::INFINITE,
            &lease,
        );
        assert!(
            matches!(
                result,
                Err(RunError::MemoryPlanOverflow {
                    what: "fold-buffer-elements",
                    ..
                })
            ),
            "got {result:?}"
        );
        assert_eq!(report.completed, 0);
        let receipt = lease.receipt();
        assert_eq!(receipt.requested_bytes, 0);
        assert_eq!(receipt.used_bytes, 0);
        assert_eq!(receipt.refusals, 0);
    }

    #[test]
    fn simultaneous_typed_refusals_report_lowest_tile_and_drain() {
        for workers in [2, 4] {
            let mut config = PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, 0xFA11);
            config.arena.limit_bytes = Some(0);
            let pool = TilePool::new(config);
            let gate = CancelGate::new();
            let kernel = SimultaneousAllocationRefusal {
                tiles: workers as u64,
                barrier: std::sync::Barrier::new(workers),
            };
            let (result, report) = pool.run_declared_budgeted(
                &kernel,
                &gate,
                RunId(23),
                Budget::new().with_cost_quota(1 << 20),
            );
            match result {
                Err(RunError::TileFailed {
                    tile: 0,
                    failure:
                        TileFailure::Allocation(fs_alloc::AllocError::Exhausted {
                            limit_bytes: 0, ..
                        }),
                    completed: 0,
                    ..
                }) => {}
                other => panic!("expected deterministic allocation refusal, got {other:?}"),
            }
            assert!(gate.is_requested());
            assert_eq!(report.completed, 0);
            assert_eq!(report.total, workers as u64);
            assert!(pool.arena_pool().stats().quiescent());
            assert_eq!(pool.run(&NoAllocation).expect("pool remains reusable"), 1);
        }
    }

    #[test]
    fn panic_precedence_over_typed_refusal_is_explicit_and_drained() {
        let mut config = PoolConfig::new(2, CcdTopology::APPLE_M_CLASS, 0xFA12);
        config.arena.limit_bytes = Some(0);
        let pool = TilePool::new(config);
        let gate = CancelGate::new();
        let kernel = MixedPanicAndRefusal {
            barrier: std::sync::Barrier::new(2),
        };
        let (result, report) = pool.run_declared_budgeted(
            &kernel,
            &gate,
            RunId(24),
            Budget::new().with_cost_quota(1 << 20),
        );
        match result {
            Err(RunError::TilePanicked {
                tile: 1,
                message,
                completed: 0,
                ..
            }) => assert!(message.contains("mixed failure panic"), "{message}"),
            other => panic!("panic class must precede typed refusal, got {other:?}"),
        }
        assert!(gate.is_requested());
        assert_eq!(report.completed, 0);
        assert!(pool.arena_pool().stats().quiescent());
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

    fn placement_identity_fixture() -> (PoolConfig, fs_alloc::HugepageDecision) {
        let mut config = PoolConfig::new(
            3,
            CcdTopology {
                ccds: 3,
                cores_per_ccd: 5,
            },
            0xA110_CAFE,
        );
        config.quantum_weights = vec![2, 3, 5];
        config.mode = ExecMode::Fast;
        config.arena = fs_alloc::ArenaConfig {
            chunk_bytes: 2 << 20,
            max_chunk_bytes: 32 << 20,
            limit_bytes: Some(96 << 20),
            free_list_max_bytes: 48 << 20,
            hugepage: fs_alloc::HugepagePolicy::Auto,
        };
        config.pin_groups = vec![vec![3, 1], vec![7, 9, 11]];
        let hugepage = fs_alloc::HugepageDecision {
            policy: fs_alloc::HugepagePolicy::Auto,
            outcome: fs_alloc::HugepageOutcome::ThpNotEnabled,
            detail: "fixture detail alpha".to_string(),
        };
        (config, hugepage)
    }

    fn fixture_placement_identity(
        config: &PoolConfig,
        hugepage: &fs_alloc::HugepageDecision,
        counts: &PlacementCounts,
    ) -> String {
        placement_identity_with_schema(
            config,
            hugepage,
            TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM,
            TILEPOOL_PLACEMENT_IDENTITY_VERSION,
            TILEPOOL_PLACEMENT_IDENTITY_DOMAIN,
            "ccd-pin-requested",
            counts,
        )
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn tilepool_placement_identity_fields_move_independently() {
        let (config, hugepage) = placement_identity_fixture();
        let counts = PlacementCounts::from_inputs(&config, &hugepage);
        let canonical = fixture_placement_identity(&config, &hugepage, &counts);
        let assert_moves = |field: &str, changed: String| {
            assert_ne!(
                changed, canonical,
                "semantic placement field {field} did not move the identity"
            );
        };

        assert_moves(
            "digest-domain",
            placement_identity_with_schema(
                &config,
                &hugepage,
                TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM,
                TILEPOOL_PLACEMENT_IDENTITY_VERSION,
                "org.frankensim.fs-exec.tilepool-placement.w2",
                "ccd-pin-requested",
                &counts,
            ),
        );
        assert_moves(
            "identity-prefix-stem",
            placement_identity_with_schema(
                &config,
                &hugepage,
                "xs-exec-tilepool-v",
                TILEPOOL_PLACEMENT_IDENTITY_VERSION,
                TILEPOOL_PLACEMENT_IDENTITY_DOMAIN,
                "ccd-pin-requested",
                &counts,
            ),
        );
        assert_moves(
            "identity-version",
            placement_identity_with_schema(
                &config,
                &hugepage,
                TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM,
                TILEPOOL_PLACEMENT_IDENTITY_VERSION + 1,
                TILEPOOL_PLACEMENT_IDENTITY_DOMAIN,
                "ccd-pin-requested",
                &counts,
            ),
        );

        let mut changed_counts = counts.clone();
        changed_counts.workers += 1;
        assert_moves(
            "workers",
            fixture_placement_identity(&config, &hugepage, &changed_counts),
        );
        let mut changed = config.clone();
        changed.topo.ccds += 1;
        assert_moves(
            "topology-ccds",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
        let mut changed = config.clone();
        changed.topo.cores_per_ccd += 1;
        assert_moves(
            "topology-cores-per-ccd",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
        let mut changed = config.clone();
        changed.mode = ExecMode::Deterministic;
        assert_moves(
            "mode-tag",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );

        let mut changed_counts = counts.clone();
        changed_counts.quantum_weights += 1;
        assert_moves(
            "quantum-weight-count",
            fixture_placement_identity(&config, &hugepage, &changed_counts),
        );
        let mut changed = config.clone();
        changed.quantum_weights.swap(0, 1);
        assert_moves(
            "quantum-weights-in-order",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );

        let mut changed = config.clone();
        changed.arena.chunk_bytes += 4096;
        assert_moves(
            "arena-chunk-bytes",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
        let mut changed = config.clone();
        changed.arena.max_chunk_bytes += 4096;
        assert_moves(
            "arena-max-chunk-bytes",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
        let mut changed = config.clone();
        changed.arena.limit_bytes = None;
        assert_moves(
            "arena-limit-presence",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
        let mut changed = config.clone();
        changed.arena.limit_bytes = changed.arena.limit_bytes.map(|limit| limit + 4096);
        assert_moves(
            "arena-limit-bytes",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
        let mut changed = config.clone();
        changed.arena.free_list_max_bytes += 4096;
        assert_moves(
            "arena-free-list-max-bytes",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
        let mut changed = config.clone();
        changed.arena.hugepage = fs_alloc::HugepagePolicy::Never;
        assert_moves(
            "arena-hugepage-policy-tag",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );

        let mut changed_hugepage = hugepage.clone();
        changed_hugepage.policy = fs_alloc::HugepagePolicy::Never;
        assert_moves(
            "hugepage-decision-policy",
            fixture_placement_identity(&config, &changed_hugepage, &counts),
        );
        let mut changed_hugepage = hugepage.clone();
        changed_hugepage.outcome = fs_alloc::HugepageOutcome::AlignedForThp;
        assert_moves(
            "hugepage-decision-outcome",
            fixture_placement_identity(&config, &changed_hugepage, &counts),
        );
        let mut changed_counts = counts.clone();
        changed_counts.hugepage_json_bytes += 1;
        assert_moves(
            "hugepage-json-byte-count",
            fixture_placement_identity(&config, &hugepage, &changed_counts),
        );
        let mut changed_hugepage = hugepage.clone();
        changed_hugepage.detail = "fixture detail omega".to_string();
        assert_eq!(changed_hugepage.detail.len(), hugepage.detail.len());
        assert_moves(
            "hugepage-decision-detail-json",
            fixture_placement_identity(&config, &changed_hugepage, &counts),
        );

        assert_moves(
            "pinning-intent",
            placement_identity_with_schema(
                &config,
                &hugepage,
                TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM,
                TILEPOOL_PLACEMENT_IDENTITY_VERSION,
                TILEPOOL_PLACEMENT_IDENTITY_DOMAIN,
                "pin-unrequested",
                &counts,
            ),
        );
        let mut changed_counts = counts.clone();
        changed_counts.pin_groups += 1;
        assert_moves(
            "pin-group-count",
            fixture_placement_identity(&config, &hugepage, &changed_counts),
        );
        let mut changed_counts = counts.clone();
        changed_counts.pin_cpus[0] += 1;
        assert_moves(
            "pin-cpu-counts",
            fixture_placement_identity(&config, &hugepage, &changed_counts),
        );
        let mut changed = config.clone();
        changed.pin_groups[0].swap(0, 1);
        assert_moves(
            "pin-cpu-ids-in-order",
            fixture_placement_identity(
                &changed,
                &hugepage,
                &PlacementCounts::from_inputs(&changed, &hugepage),
            ),
        );
    }

    #[test]
    fn tilepool_placement_seed_is_nonsemantic() {
        let (mut first, _) = placement_identity_fixture();
        first.arena.hugepage = fs_alloc::HugepagePolicy::Never;
        let mut second = first.clone();
        second.seed ^= u64::MAX;
        assert_eq!(
            TilePool::new(first).placement_identity(),
            TilePool::new(second).placement_identity(),
            "the scheduling-stream seed must not partition placement tune rows"
        );
    }

    #[test]
    fn tilepool_placement_identity_versions_fail_closed() {
        let (mut config, _) = placement_identity_fixture();
        config.arena.hugepage = fs_alloc::HugepagePolicy::Never;
        let pool = TilePool::new(config);
        let identity = pool.placement_identity();
        let admitted =
            pool.admit_retained_placement_identity(TILEPOOL_PLACEMENT_IDENTITY_VERSION, &identity);
        assert_eq!(admitted, Ok(()));
        for version in [
            TILEPOOL_PLACEMENT_IDENTITY_VERSION - 1,
            TILEPOOL_PLACEMENT_IDENTITY_VERSION + 1,
        ] {
            assert!(
                pool.admit_retained_placement_identity(version, &identity)
                    .is_err(),
                "retained producer version {version} must fail closed"
            );
        }
        let mut tampered = identity;
        tampered.push('x');
        let refused =
            pool.admit_retained_placement_identity(TILEPOOL_PLACEMENT_IDENTITY_VERSION, &tampered);
        assert!(refused.is_err());
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

    /// G5 (tkr7): the parked lane is deterministic across reruns on ONE
    /// crew and bitwise-identical to the spawned std lane — one worker
    /// protocol, three lifetime strategies — including runs with fewer
    /// tiles than crew workers (excess workers no-op).
    #[test]
    fn parked_local_lane_matches_std_lane_bitwise_across_reruns() {
        let p = pool(4);
        let expected = p.run(&SumKernel { tiles: 257 }).expect("std lane");
        let expected_small = p.run(&SumKernel { tiles: 2 }).expect("std lane small");
        p.with_parked_crew_local(|parked| {
            let first = parked.run(&SumKernel { tiles: 257 }).expect("parked");
            let second = parked.run(&SumKernel { tiles: 257 }).expect("parked rerun");
            assert_eq!(first, expected, "parked lane bitwise-matches the std lane");
            assert_eq!(second, expected, "parked lane deterministic across reruns");
            let small = parked
                .run(&SumKernel { tiles: 2 })
                .expect("fewer tiles than crew workers");
            assert_eq!(small, expected_small, "excess crew workers no-op cleanly");
            let (_, report) = parked.run_with_gate(&SumKernel { tiles: 257 }, &CancelGate::new());
            assert_eq!(report.completed, 257);
            assert_eq!(
                report.tiles_by_worker.iter().sum::<u64>(),
                257,
                "every tile accounted to a worker"
            );
        });
        assert!(p.arena_pool().stats().quiescent(), "arenas quiescent");
    }

    /// G4 (tkr7): per-tile panic containment holds unchanged on the
    /// parked lane, and the SAME crew keeps serving runs afterwards.
    #[test]
    fn parked_lane_contains_tile_panics_and_the_crew_survives() {
        let p = pool(4);
        p.with_parked_crew_local(|parked| {
            let err = parked
                .run(&PanicAt { tiles: 512, at: 7 })
                .expect_err("tile panic surfaces");
            match &err {
                RunError::TilePanicked {
                    tile: 7, message, ..
                } => {
                    assert!(message.contains("scoped containment probe"), "{message}");
                }
                other => panic!("expected TilePanicked{{tile:7}}, got {other:?}"),
            }
            let ok = parked
                .run(&SumKernel { tiles: 64 })
                .expect("crew survives a contained tile panic");
            assert_eq!(ok, (0..64).map(|t| t + 1).sum::<u64>());
        });
        assert!(p.arena_pool().stats().quiescent());
    }

    /// G4 (tkr7): a gate request mid-run drains a parked run exactly like
    /// a spawned run, and the crew serves the next run.
    #[test]
    fn parked_lane_drains_on_gate_request_and_reuses_the_crew() {
        struct SelfCancel {
            tiles: u64,
        }
        impl TileKernel for SelfCancel {
            type Out = u64;

            fn tiles(&self) -> TilePlan {
                TilePlan::new("test/self-cancel", self.tiles)
            }

            fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
                if tile == 0 {
                    return ControlFlow::Break(crate::Cancelled);
                }
                ControlFlow::Continue(1)
            }
        }
        let p = pool(4);
        p.with_parked_crew_local(|parked| {
            let (out, _report) =
                parked.run_with_gate(&SelfCancel { tiles: 16_384 }, &CancelGate::new());
            match out {
                Err(RunError::Cancelled {
                    completed, total, ..
                }) => {
                    assert_eq!(total, 16_384);
                    assert!(completed < total, "drain preempted completion");
                }
                other => panic!("expected Cancelled, got {other:?}"),
            }
            parked
                .run(&SumKernel { tiles: 64 })
                .expect("crew serves runs after a drained cancellation");
        });
        assert!(p.arena_pool().stats().quiescent());
    }

    /// G4+G5 (tkr7): the parked crew under a REAL task scope — bitwise
    /// equality with the spawned lanes, mid-run task cancellation drains
    /// through each parked worker's own scoped-CPU context, and a
    /// pre-cancelled task refuses the whole crew scope at entry.
    #[test]
    fn parked_task_crew_bridges_cancellation_and_matches_other_lanes() {
        let expected = pool(4).run(&SumKernel { tiles: 257 }).expect("std lane");
        in_task(move |cx| {
            let p = pool(4);
            let out = p
                .with_parked_crew(&cx, |parked| {
                    let first = parked
                        .run(&SumKernel { tiles: 257 })
                        .expect("parked task lane");
                    assert_eq!(first, expected, "parked task lane bitwise-matches");

                    // Mid-run task cancel: drains at tile boundaries via the
                    // park-time CpuCx bridge, then the task is revived and
                    // the SAME crew serves the next run.
                    let kernel = CancelTaskAt {
                        tiles: 16_384,
                        at: 0,
                        task: cx.clone(),
                    };
                    let (out, _) = parked.run_with_gate(&kernel, &CancelGate::new());
                    match out {
                        Err(RunError::Cancelled {
                            completed, total, ..
                        }) => {
                            assert_eq!(total, 16_384);
                            assert!(completed < total, "task cancel drained the parked run");
                        }
                        other => panic!("expected Cancelled, got {other:?}"),
                    }
                    cx.set_cancel_requested(false);
                    parked
                        .run(&SumKernel { tiles: 64 })
                        .expect("crew serves runs after task revival")
                })
                .expect("crew scope completes");
            assert_eq!(out, (0..64).map(|t| t + 1).sum::<u64>());
            assert!(p.arena_pool().stats().quiescent(), "arenas quiescent");

            // Entry refusal: a pre-cancelled task parks nothing and runs
            // nothing — f is never called.
            cx.set_cancel_requested(true);
            let refused = p.with_parked_crew(&cx, |_parked| {
                panic!("f must not run under a pre-cancelled task")
            });
            assert_eq!(refused, Err(CrewScopeError::Cancelled));
            cx.set_cancel_requested(false);
        });
    }
}
