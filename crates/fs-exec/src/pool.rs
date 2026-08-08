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
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

/// Semantic version of the tile-pool placement/tuning identity.
///
/// Version 2 is the already-shipped key and payload. Making the version
/// explicit does not rotate or re-key it.
pub const TILEPOOL_PLACEMENT_IDENTITY_VERSION: u32 = 2;

/// BLAKE3 derive-key domain for the exact v2 placement payload.
pub const TILEPOOL_PLACEMENT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-exec.tilepool-placement.v2";

const TILEPOOL_PLACEMENT_IDENTITY_PREFIX_STEM: &str = "fs-exec-tilepool-v";

/// Schema version of executor-minted TilePool completion evidence.
///
/// Version 2 binds the exact normalized placement identity, pool seed,
/// call-local replay root, request phase, full parked-crew callback set, and
/// lossless `u128` cumulative lease counters.
pub const TILEPOOL_COMPLETION_WITNESS_VERSION: u32 = 2;

const TILEPOOL_COMPLETION_WITNESS_DOMAIN: &str =
    "org.frankensim.fs-exec.tilepool-completion-witness.v2";
const TILEPOOL_COMPLETION_PLAN_DOMAIN: &str = "org.frankensim.fs-exec.tilepool-completion-plan.v2";
const TILEPOOL_COMPLETION_CALL_REPLAY_DOMAIN: &str =
    "org.frankensim.fs-exec.tilepool-call-replay.v1";

/// Claims deliberately excluded from every TilePool completion witness.
///
/// The witness covers only the run-local executor lifecycle that mints it.
/// In particular, it cannot discover work that was never admitted to this
/// pool run and it says nothing about publication outside the executor.
const TILEPOOL_COMPLETION_WITNESS_BASE_NO_CLAIMS: &[&str] = &[
    "external-or-unregistered-threads",
    "cancel-gate-instance-identity",
    "wall-clock-or-cancellation-latency",
    "scheduler-trace-or-steal-replay",
    "application-state-publication",
    "scientific-output-correctness",
    "durable-invocation-admission",
];

/// No-claim set for standalone witnessed calls that do not consume an affine
/// invocation permit.
pub const TILEPOOL_COMPLETION_WITNESS_NO_CLAIMS: &[&str] = &[
    "external-or-unregistered-threads",
    "cancel-gate-instance-identity",
    "wall-clock-or-cancellation-latency",
    "scheduler-trace-or-steal-replay",
    "application-state-publication",
    "scientific-output-correctness",
    "durable-invocation-admission",
    "cross-call-uniqueness-without-affine-invocation-permit",
];

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
    /// A parked crew already has an admitted dispatch. Concurrent or
    /// recursive dispatch is refused before root allocation or worker wakeup,
    /// leaving the active run and crew reusable.
    ParkedCrewBusy {
        /// Kernel whose attempted dispatch was refused.
        kernel: &'static str,
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
            RunError::ParkedCrewBusy { kernel } => write!(
                f,
                "kernel `{kernel}` was not dispatched: this parked crew already has an active run"
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

/// Terminal class bound into an executor-minted completion witness.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TilePoolCompletionDisposition {
    /// Every admitted tile completed and reduction returned normally.
    Completed,
    /// The run's cancellation gate was requested and the worker set drained.
    Cancelled,
    /// At least one tile returned a typed refusal.
    TileFailed,
    /// At least one tile panic was contained.
    TilePanicked,
    /// A scoped worker could not be created; already-launched workers joined.
    WorkerSpawnRefused,
    /// An overlapping or recursive parked-crew dispatch was refused.
    ParkedCrewBusy,
    /// The operation lease refused root metadata before worker launch.
    MemoryRefused,
    /// Checked root planning refused an unrepresentable run.
    MemoryPlanOverflow,
    /// Fallible root backing allocation was refused.
    MemoryAllocationRefused,
    /// Every tile completed, but the deterministic reduction panicked.
    ReductionPanicked,
    /// The executor found a missing output slot after worker join.
    Incomplete,
}

impl TilePoolCompletionDisposition {
    /// Stable canonical spelling used by completion logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Completed => "completed",
            Self::Cancelled => "cancelled",
            Self::TileFailed => "tile-failed",
            Self::TilePanicked => "tile-panicked",
            Self::WorkerSpawnRefused => "worker-spawn-refused",
            Self::ParkedCrewBusy => "parked-crew-busy",
            Self::MemoryRefused => "memory-refused",
            Self::MemoryPlanOverflow => "memory-plan-overflow",
            Self::MemoryAllocationRefused => "memory-allocation-refused",
            Self::ReductionPanicked => "reduction-panicked",
            Self::Incomplete => "incomplete",
        }
    }

    #[must_use]
    const fn tag(self) -> u64 {
        match self {
            Self::Completed => 0,
            Self::Cancelled => 1,
            Self::TileFailed => 2,
            Self::TilePanicked => 3,
            Self::WorkerSpawnRefused => 4,
            Self::ParkedCrewBusy => 5,
            Self::MemoryRefused => 6,
            Self::MemoryPlanOverflow => 7,
            Self::MemoryAllocationRefused => 8,
            Self::ReductionPanicked => 9,
            Self::Incomplete => 10,
        }
    }
}

/// Stable reason a retained completion witness failed verification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TilePoolCompletionWitnessError {
    /// The producer schema is not understood by this verifier.
    UnsupportedVersion {
        /// Encountered schema version.
        found: u32,
    },
    /// The retained plan root does not bind the declared run and tile plan.
    PlanRootMismatch,
    /// The retained witness root does not bind the retained fields.
    RootMismatch,
    /// One lifecycle/count/quiescence invariant was violated.
    Invariant {
        /// Stable invariant label.
        name: &'static str,
    },
    /// The result/report/witness bundle disagrees internally.
    BundleInvariant {
        /// Stable bundle-invariant label.
        name: &'static str,
    },
}

impl fmt::Display for TilePoolCompletionWitnessError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedVersion { found } => {
                write!(
                    formatter,
                    "unsupported TilePool completion witness version {found}"
                )
            }
            Self::PlanRootMismatch => {
                formatter.write_str("TilePool completion witness plan root mismatch")
            }
            Self::RootMismatch => formatter.write_str("TilePool completion witness root mismatch"),
            Self::Invariant { name } => {
                write!(formatter, "TilePool completion witness violated `{name}`")
            }
            Self::BundleInvariant { name } => {
                write!(formatter, "TilePool witnessed run violated `{name}`")
            }
        }
    }
}

impl core::error::Error for TilePoolCompletionWitnessError {}

/// First request phase retained by a completion witness.
///
/// This is a logical ordering relative to executor checkpoints, not a
/// wall-clock timestamp.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TilePoolRequestPhase {
    /// No request existed when the terminal bundle was sealed.
    NotRequested,
    /// The gate was already requested when the run entered the executor.
    BeforeEntry,
    /// The request appeared after entry but before terminal outcome selection.
    BeforeTerminalDecision,
    /// The request appeared only after terminal outcome selection but before
    /// the immutable witness was sealed.
    AfterTerminalDecision,
}

impl TilePoolRequestPhase {
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::NotRequested => "not-requested",
            Self::BeforeEntry => "before-entry",
            Self::BeforeTerminalDecision => "before-terminal-decision",
            Self::AfterTerminalDecision => "after-terminal-decision",
        }
    }

    #[must_use]
    const fn tag(self) -> u64 {
        match self {
            Self::NotRequested => 0,
            Self::BeforeEntry => 1,
            Self::BeforeTerminalDecision => 2,
            Self::AfterTerminalDecision => 3,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CompletionArenaSnapshot {
    live: u64,
    reserved_bytes: u64,
    free_bytes: u64,
    quiescent: bool,
}

impl CompletionArenaSnapshot {
    fn capture(pool: &fs_alloc::ArenaPool) -> Self {
        let stats = pool.stats();
        Self {
            live: u64::try_from(stats.arenas_live).unwrap_or(u64::MAX),
            reserved_bytes: u64::try_from(stats.reserved_bytes).unwrap_or(u64::MAX),
            free_bytes: u64::try_from(stats.free_bytes).unwrap_or(u64::MAX),
            quiescent: stats.quiescent(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CompletionLeaseSnapshot {
    limit_bytes: Option<u64>,
    requested_bytes: u128,
    peak_bytes: u64,
    used_bytes: u64,
    refusals: u128,
    release_invariant_violations: u128,
}

impl CompletionLeaseSnapshot {
    fn capture(lease: &fs_alloc::OperationMemoryLease) -> Self {
        let receipt = lease.receipt();
        Self {
            limit_bytes: receipt.limit_bytes,
            requested_bytes: receipt.requested_bytes,
            peak_bytes: receipt.peak_bytes,
            used_bytes: receipt.used_bytes,
            refusals: receipt.refusals,
            release_invariant_violations: receipt.release_invariant_violations,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CompletionScopeIdentity {
    kind: &'static str,
    parent_region_id: Option<u64>,
    parent_task_id: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CompletionBefore {
    scope: CompletionScopeIdentity,
    pool_placement_identity: String,
    arena: CompletionArenaSnapshot,
    lease: CompletionLeaseSnapshot,
    cancellation_requested_at_entry: bool,
    planned_crew_callbacks: u64,
}

impl CompletionScopeIdentity {
    const STD_SCOPED: Self = Self {
        kind: "std-thread-scope",
        parent_region_id: None,
        parent_task_id: None,
    };

    const STD_PARKED: Self = Self {
        kind: "std-thread-parked-crew",
        parent_region_id: None,
        parent_task_id: None,
    };

    fn task_scoped<Caps>(cx: &asupersync::Cx<Caps>) -> Self {
        Self {
            kind: "asupersync-task-scope",
            parent_region_id: Some(cx.region_id().as_u64()),
            parent_task_id: Some(cx.task_id().as_u64()),
        }
    }

    fn task_parked<Caps>(cx: &asupersync::Cx<Caps>) -> Self {
        Self {
            kind: "asupersync-task-parked-crew",
            parent_region_id: Some(cx.region_id().as_u64()),
            parent_task_id: Some(cx.task_id().as_u64()),
        }
    }
}

/// Immutable, versioned evidence minted by the executor only after the real
/// TilePool launch/join path has closed.
///
/// Private fields and the absence of a public constructor prevent callers
/// from manufacturing a `drained=true` token. The witness binds the declared
/// run and plan, actual worker entry/exit counts, exact claimed tile outcomes,
/// run-local arena-scope closure, root-charge release, before/after allocator
/// observations, the selected terminal error, and an explicit no-claim set.
///
/// `failed_tiles()` counts panic outcomes plus the one typed refusal retained
/// with exact provenance by the existing lowest-tile `RunError` contract.
/// Additional simultaneous `ControlFlow::Break` values are retained exactly
/// in `break_tiles()` but are not relabelled as independently proven faults.
/// Consequently `cancelled_tiles()` is the conservative remainder without
/// retained failure provenance, not a claim that every such tile personally
/// observed the cancellation gate.
///
/// Global arena/lease before-after values are observations, not a claim that
/// no concurrent caller exists. Run-local authority comes from zero live
/// worker guards and tile scopes plus explicit release of the executor's root
/// charge. The witness never attests external threads, wall-clock latency,
/// scheduler replay, application publication, or scientific correctness.
///
/// The type is intentionally non-forgeable and immutable outside this module:
///
/// ```compile_fail
/// use fs_exec::TilePoolCompletionWitness;
///
/// fn forge_or_mutate(witness: &mut TilePoolCompletionWitness) {
///     witness.version = 99;
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[must_use]
pub struct TilePoolCompletionWitness {
    version: u32,
    producer_version: &'static str,
    pool_placement_identity_version: u32,
    pool_placement_identity: String,
    pool_seed: u64,
    kernel: &'static str,
    kernel_id: u64,
    declared_run: RunId,
    mode: &'static str,
    scope: CompletionScopeIdentity,
    planned_tiles: u64,
    plan_root: [u8; 32],
    call_replay_root: [u8; 32],
    affine_invocation_permit_root: Option<[u8; 32]>,
    admission_completed: bool,
    admitted_tiles: u64,
    unadmitted_tiles: u64,
    claimed_tiles: u64,
    completed_tiles: u64,
    break_tiles: u64,
    panicked_tiles: u64,
    planned_workers: u64,
    launched_workers: u64,
    joined_workers: u64,
    worker_admission_closed: bool,
    live_worker_guards_at_seal: u64,
    planned_crew_callbacks: u64,
    entered_crew_callbacks: u64,
    exited_crew_callbacks: u64,
    tile_scopes_opened: u64,
    live_tile_scopes_at_seal: u64,
    cancellation_requested_at_entry: bool,
    cancellation_requested_at_terminal: bool,
    cancellation_requested: bool,
    request_phase: TilePoolRequestPhase,
    cancellation_observed_workers: u64,
    root_metadata_bytes: u64,
    root_charge_admitted: bool,
    root_charge_released: bool,
    arena_before: CompletionArenaSnapshot,
    arena_after: CompletionArenaSnapshot,
    lease_before: CompletionLeaseSnapshot,
    lease_after: CompletionLeaseSnapshot,
    disposition: TilePoolCompletionDisposition,
    first_failure_kind: Option<&'static str>,
    first_failure_tile: Option<u64>,
    terminal_error: Option<RunError>,
    root: [u8; 32],
}

/// Affine authority to execute one permit-bound witnessed TilePool call.
///
/// The root invocation layer constructs this token only after atomically
/// reserving its one execution attempt. The token is deliberately neither
/// [`Clone`] nor [`Copy`], and every permit-bound run method consumes it.
/// Its per-run permit root is bound into both the call-replay root and
/// terminal witness. The root invocation layer derives that permit root from
/// the invocation occurrence plus the declared run ordinal; it is not the
/// bare invocation root.
///
/// Standalone `*_witnessed` compatibility methods do not need a permit and
/// explicitly retain the cross-call-uniqueness no-claim.
///
/// ```compile_fail
/// use fs_exec::TilePoolInvocationPermit;
///
/// fn consume(_: TilePoolInvocationPermit) {}
///
/// fn cannot_reuse(permit: TilePoolInvocationPermit) {
///     consume(permit);
///     consume(permit);
/// }
/// ```
#[derive(Debug, PartialEq, Eq)]
#[must_use]
pub struct TilePoolInvocationPermit {
    permit_root: [u8; 32],
}

impl TilePoolInvocationPermit {
    /// Mint one permit from a per-run root derived by invocation authority.
    ///
    /// This constructor stays crate-private so ordinary callers cannot mint
    /// one-shot authority. `InvocationBudget` is responsible for deriving the
    /// root from its invocation occurrence plus run ordinal and refusing a
    /// second mint for that ordinal.
    #[must_use]
    pub(crate) const fn from_permit_root(permit_root: [u8; 32]) -> Self {
        Self { permit_root }
    }

    /// Per-run permit root bound into the witnessed call.
    #[must_use]
    pub const fn permit_root(&self) -> [u8; 32] {
        self.permit_root
    }

    const fn into_root(self) -> [u8; 32] {
        self.permit_root
    }
}

/// One coherent executor result, measured report, and completion witness.
///
/// Private fields prevent callers from pairing a valid witness with a report
/// or outcome from another run. Use [`Self::into_parts`] to consume the
/// executor-minted bundle.
#[derive(Debug)]
#[must_use]
pub struct WitnessedRun<Out> {
    outcome: Result<Out, RunError>,
    report: RunReport,
    witness: TilePoolCompletionWitness,
}

impl<Out> WitnessedRun<Out> {
    /// Borrow the terminal kernel outcome.
    #[must_use]
    pub const fn outcome(&self) -> &Result<Out, RunError> {
        &self.outcome
    }

    /// Borrow the measured, non-semantic scheduling report.
    #[must_use]
    pub const fn report(&self) -> &RunReport {
        &self.report
    }

    /// Borrow the immutable semantic completion witness.
    #[must_use]
    pub const fn witness(&self) -> &TilePoolCompletionWitness {
        &self.witness
    }

    /// Consume the bundle without weakening its prior verification.
    #[must_use]
    pub fn into_parts(self) -> (Result<Out, RunError>, RunReport, TilePoolCompletionWitness) {
        (self.outcome, self.report, self.witness)
    }

    /// Recheck the witness plus exact report/outcome coupling.
    pub fn verify_bundle(&self) -> Result<(), TilePoolCompletionWitnessError> {
        self.witness.verify()?;
        if self.report.kernel != self.witness.kernel {
            return Err(completion_bundle_invariant("report-kernel"));
        }
        if self.report.mode != self.witness.mode {
            return Err(completion_bundle_invariant("report-mode"));
        }
        if self.report.declared_run != self.witness.declared_run {
            return Err(completion_bundle_invariant("report-declared-run"));
        }
        if self.report.completed != self.witness.completed_tiles {
            return Err(completion_bundle_invariant("report-completed"));
        }
        if self.report.total != self.witness.planned_tiles {
            return Err(completion_bundle_invariant("report-total"));
        }
        let reported_completed = self
            .report
            .tiles_by_worker
            .iter()
            .try_fold(0_u64, |total, count| total.checked_add(*count))
            .ok_or_else(|| completion_bundle_invariant("report-worker-count-overflow"))?;
        if self.witness.admission_completed {
            if u64::try_from(self.report.tiles_by_worker.len()).ok()
                != Some(self.witness.planned_workers)
            {
                return Err(completion_bundle_invariant("report-worker-cardinality"));
            }
            if reported_completed != self.witness.completed_tiles {
                return Err(completion_bundle_invariant("report-worker-conservation"));
            }
        } else if !self.report.tiles_by_worker.is_empty() || reported_completed != 0 {
            return Err(completion_bundle_invariant("prelaunch-report-workers"));
        }
        match (&self.outcome, self.witness.terminal_error()) {
            (Ok(_), None) => {}
            (Err(outcome), Some(retained)) if outcome == retained => {}
            _ => return Err(completion_bundle_invariant("terminal-outcome")),
        }
        Ok(())
    }
}

/// Runner surface for consumers that require executor-minted completion
/// evidence rather than the legacy outcome/report pair.
pub trait CompletionKernelRunner {
    /// Normalized worker count used for preflight sizing.
    fn workers(&self) -> usize;

    /// Execute under an explicit gate and return one verified bundle.
    fn run_with_gate_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>;

    /// Consume one affine invocation permit and execute the exact call once.
    fn run_with_gate_witnessed_once<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>;
}

impl TilePoolCompletionWitness {
    /// Current witness schema version.
    pub const VERSION: u32 = TILEPOOL_COMPLETION_WITNESS_VERSION;

    /// Fixed explicit no-claim set.
    ///
    /// This is the compatibility set used by standalone calls. Permit-bound
    /// witnesses expose the narrower run-specific set through
    /// [`Self::no_claims`].
    pub const NO_CLAIMS: &'static [&'static str] = TILEPOOL_COMPLETION_WITNESS_NO_CLAIMS;

    /// Explicit no-claim set applicable to this exact witnessed call.
    #[must_use]
    pub fn no_claims(&self) -> &'static [&'static str] {
        if self.affine_invocation_permit_root.is_some() {
            TILEPOOL_COMPLETION_WITNESS_BASE_NO_CLAIMS
        } else {
            TILEPOOL_COMPLETION_WITNESS_NO_CLAIMS
        }
    }

    /// Retained producer schema version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// fs-exec crate version that minted this witness.
    #[must_use]
    pub const fn producer_version(&self) -> &'static str {
        self.producer_version
    }

    /// Placement-identity schema used by the exact normalized pool.
    #[must_use]
    pub const fn pool_placement_identity_version(&self) -> u32 {
        self.pool_placement_identity_version
    }

    /// Exact normalized placement identity of the pool that executed this
    /// call.
    #[must_use]
    pub fn pool_placement_identity(&self) -> &str {
        &self.pool_placement_identity
    }

    /// Declared stream seed from the executing pool configuration.
    #[must_use]
    pub const fn pool_seed(&self) -> u64 {
        self.pool_seed
    }

    /// Kernel name from the actual `TilePlan`.
    #[must_use]
    pub const fn kernel(&self) -> &'static str {
        self.kernel
    }

    /// Stable kernel identity from the actual `TilePlan`.
    #[must_use]
    pub const fn kernel_id(&self) -> u64 {
        self.kernel_id
    }

    /// Caller-declared logical run identity.
    #[must_use]
    pub const fn declared_run(&self) -> RunId {
        self.declared_run
    }

    /// Execution-mode name bound into the witness.
    #[must_use]
    pub const fn mode(&self) -> &'static str {
        self.mode
    }

    /// Worker-lifetime/scope strategy used by this exact launch.
    #[must_use]
    pub const fn scope_kind(&self) -> &'static str {
        self.scope.kind
    }

    /// Owning asupersync region when this run was launched inside one.
    #[must_use]
    pub const fn parent_region_id(&self) -> Option<u64> {
        self.scope.parent_region_id
    }

    /// Owning asupersync task when this run was launched inside one.
    #[must_use]
    pub const fn parent_task_id(&self) -> Option<u64> {
        self.scope.parent_task_id
    }

    /// Tiles declared by the actual plan.
    #[must_use]
    pub const fn planned_tiles(&self) -> u64 {
        self.planned_tiles
    }

    /// Whether checked root allocation/admission completed before launch.
    #[must_use]
    pub const fn admission_completed(&self) -> bool {
        self.admission_completed
    }

    /// Plan tiles admitted to the worker protocol.
    #[must_use]
    pub const fn admitted_tiles(&self) -> u64 {
        self.admitted_tiles
    }

    /// Plan tiles refused before worker-protocol admission.
    #[must_use]
    pub const fn unadmitted_tiles(&self) -> u64 {
        self.unadmitted_tiles
    }

    /// Tiles actually removed from an owned/stolen run.
    #[must_use]
    pub const fn claimed_tiles(&self) -> u64 {
        self.claimed_tiles
    }

    /// Tiles that returned `ControlFlow::Continue` and populated a slot.
    #[must_use]
    pub const fn completed_tiles(&self) -> u64 {
        self.completed_tiles
    }

    /// Tiles that returned `ControlFlow::Break`.
    #[must_use]
    pub const fn break_tiles(&self) -> u64 {
        self.break_tiles
    }

    /// Tiles whose panic was contained by the worker loop.
    #[must_use]
    pub const fn panicked_tiles(&self) -> u64 {
        self.panicked_tiles
    }

    /// Typed-refusal tiles retained with exact provenance (zero or one).
    #[must_use]
    pub fn retained_refusal_tiles(&self) -> u64 {
        if matches!(
            self.terminal_error.as_ref(),
            Some(RunError::TileFailed { .. })
        ) {
            1
        } else {
            0
        }
    }

    /// Failures with exact executor-retained provenance.
    #[must_use]
    pub fn failed_tiles(&self) -> u64 {
        self.panicked_tiles + self.retained_refusal_tiles()
    }

    /// Conservative non-success remainder without retained failure provenance.
    #[must_use]
    pub fn cancelled_tiles(&self) -> u64 {
        let unclaimed = self.admitted_tiles.saturating_sub(self.claimed_tiles);
        let nonfailure_breaks = self
            .break_tiles
            .saturating_sub(self.retained_refusal_tiles());
        unclaimed.saturating_add(nonfailure_breaks)
    }

    /// Workers selected after tile-count normalization.
    #[must_use]
    pub const fn planned_workers(&self) -> u64 {
        self.planned_workers
    }

    /// Workers that actually entered the shared worker loop.
    #[must_use]
    pub const fn launched_workers(&self) -> u64 {
        self.launched_workers
    }

    /// Entered workers that exited before the launch harness returned.
    #[must_use]
    pub const fn joined_workers(&self) -> u64 {
        self.joined_workers
    }

    /// Whether run-local worker admission was permanently closed before seal.
    #[must_use]
    pub const fn worker_admission_closed(&self) -> bool {
        self.worker_admission_closed
    }

    /// Worker guards still live when the evidence was sealed.
    #[must_use]
    pub const fn live_worker_guards_at_seal(&self) -> u64 {
        self.live_worker_guards_at_seal
    }

    /// Number of parked-crew callbacks required for this dispatch.
    ///
    /// This is the complete parked crew, which may be larger than
    /// [`Self::planned_workers`] for a short tile plan.
    #[must_use]
    pub const fn planned_crew_callbacks(&self) -> u64 {
        self.planned_crew_callbacks
    }

    /// Parked-crew callbacks that entered this exact dispatch.
    #[must_use]
    pub const fn entered_crew_callbacks(&self) -> u64 {
        self.entered_crew_callbacks
    }

    /// Entered parked-crew callbacks that exited before sealing.
    #[must_use]
    pub const fn exited_crew_callbacks(&self) -> u64 {
        self.exited_crew_callbacks
    }

    /// Tile-arena scopes opened by claimed work.
    #[must_use]
    pub const fn tile_scopes_opened(&self) -> u64 {
        self.tile_scopes_opened
    }

    /// Run-local tile-arena scopes still live when sealed.
    #[must_use]
    pub const fn live_tile_scopes_at_seal(&self) -> u64 {
        self.live_tile_scopes_at_seal
    }

    /// Whether the run's actual gate was requested by seal time.
    #[must_use]
    pub const fn cancellation_requested(&self) -> bool {
        self.cancellation_requested
    }

    /// Whether the gate was already requested at executor entry.
    #[must_use]
    pub const fn cancellation_requested_at_entry(&self) -> bool {
        self.cancellation_requested_at_entry
    }

    /// Whether the gate was requested at the terminal-decision checkpoint.
    #[must_use]
    pub const fn cancellation_requested_at_terminal(&self) -> bool {
        self.cancellation_requested_at_terminal
    }

    /// First logical request phase observed by the executor.
    #[must_use]
    pub const fn request_phase(&self) -> TilePoolRequestPhase {
        self.request_phase
    }

    /// Workers that actually reached a tile boundary and observed the
    /// requested gate before exiting.
    ///
    /// This is distinct from [`Self::cancellation_requested`]: a refusal can
    /// close launch before any worker exists, and a late external request can
    /// race after workers have already exhausted their deques.
    #[must_use]
    pub const fn cancellation_observed_workers(&self) -> u64 {
        self.cancellation_observed_workers
    }

    /// Checked logical root-metadata charge for this plan.
    #[must_use]
    pub const fn root_metadata_bytes(&self) -> u64 {
        self.root_metadata_bytes
    }

    /// Whether the executor admitted its root lease charge.
    #[must_use]
    pub const fn root_charge_admitted(&self) -> bool {
        self.root_charge_admitted
    }

    /// Whether an admitted executor root charge was explicitly released.
    #[must_use]
    pub const fn root_charge_released(&self) -> bool {
        self.root_charge_released
    }

    /// Run-local transient quiescence (workers/scopes/root charge only).
    #[must_use]
    pub const fn executor_transients_quiescent(&self) -> bool {
        self.worker_admission_closed
            && self.live_worker_guards_at_seal == 0
            && self.live_tile_scopes_at_seal == 0
            && (!self.root_charge_admitted || self.root_charge_released)
    }

    /// Global ArenaPool live-arena observation before the run.
    #[must_use]
    pub const fn arena_live_before(&self) -> u64 {
        self.arena_before.live
    }

    /// Global ArenaPool live-arena observation after run-local drain.
    #[must_use]
    pub const fn arena_live_after(&self) -> u64 {
        self.arena_after.live
    }

    /// Whether the global ArenaPool happened to be quiescent before the run.
    #[must_use]
    pub const fn arena_pool_quiescent_before(&self) -> bool {
        self.arena_before.quiescent
    }

    /// Whether the global ArenaPool happened to be quiescent after this run.
    #[must_use]
    pub const fn arena_pool_quiescent_after(&self) -> bool {
        self.arena_after.quiescent
    }

    /// Shared operation-lease bytes live before this run.
    #[must_use]
    pub const fn lease_used_before(&self) -> u64 {
        self.lease_before.used_bytes
    }

    /// Shared operation-lease bytes live after executor-transient release.
    #[must_use]
    pub const fn lease_used_after(&self) -> u64 {
        self.lease_after.used_bytes
    }

    /// Shared operation-lease refusal count before this run.
    #[must_use]
    pub const fn lease_refusals_before(&self) -> u128 {
        self.lease_before.refusals
    }

    /// Shared operation-lease refusal count after this run.
    #[must_use]
    pub const fn lease_refusals_after(&self) -> u128 {
        self.lease_after.refusals
    }

    /// Cumulative granted lease bytes before this call.
    #[must_use]
    pub const fn lease_requested_before(&self) -> u128 {
        self.lease_before.requested_bytes
    }

    /// Cumulative granted lease bytes after this call.
    #[must_use]
    pub const fn lease_requested_after(&self) -> u128 {
        self.lease_after.requested_bytes
    }

    /// Cumulative lease release-invariant violations before this call.
    #[must_use]
    pub const fn lease_release_invariant_violations_before(&self) -> u128 {
        self.lease_before.release_invariant_violations
    }

    /// Cumulative lease release-invariant violations after this call.
    #[must_use]
    pub const fn lease_release_invariant_violations_after(&self) -> u128 {
        self.lease_after.release_invariant_violations
    }

    /// Terminal disposition bound to the selected `RunError`, if any.
    #[must_use]
    pub const fn disposition(&self) -> TilePoolCompletionDisposition {
        self.disposition
    }

    /// Stable disposition spelling.
    #[must_use]
    pub const fn disposition_name(&self) -> &'static str {
        self.disposition.name()
    }

    /// Selected terminal tile-failure class, when the terminal error is a
    /// typed refusal or contained tile panic.
    #[must_use]
    pub const fn first_failure_kind(&self) -> Option<&'static str> {
        self.first_failure_kind
    }

    /// Lowest logical tile retained under the selected terminal failure
    /// class, independently cross-checked against `terminal_error`.
    #[must_use]
    pub const fn first_failure_tile(&self) -> Option<u64> {
        self.first_failure_tile
    }

    /// Exact selected terminal error under existing `RunError` precedence.
    #[must_use]
    pub const fn terminal_error(&self) -> Option<&RunError> {
        self.terminal_error.as_ref()
    }

    /// Exact plan-root bytes.
    #[must_use]
    pub const fn plan_root_bytes(&self) -> [u8; 32] {
        self.plan_root
    }

    /// Lowercase plan-root hex.
    #[must_use]
    pub fn plan_root_hex(&self) -> String {
        completion_hex(&self.plan_root)
    }

    /// Exact replay root for the declared call identity.
    ///
    /// This root is intentionally reproducible. It is not a uniqueness token
    /// for standalone calls: identical declared calls have the same root.
    /// Permit-consuming calls also bind their one-shot invocation root.
    #[must_use]
    pub const fn call_replay_root_bytes(&self) -> [u8; 32] {
        self.call_replay_root
    }

    /// Lowercase call-replay-root hex.
    #[must_use]
    pub fn call_replay_root_hex(&self) -> String {
        completion_hex(&self.call_replay_root)
    }

    /// Affine invocation-permit root for a permit-consuming call.
    ///
    /// Standalone compatibility entry points retain `None` and therefore
    /// make no cross-call uniqueness claim.
    #[must_use]
    pub const fn affine_invocation_permit_root(&self) -> Option<[u8; 32]> {
        self.affine_invocation_permit_root
    }

    /// Whether this call was bound to an affine invocation permit.
    #[must_use]
    pub const fn has_affine_invocation_permit(&self) -> bool {
        self.affine_invocation_permit_root.is_some()
    }

    /// Exact witness-root bytes.
    #[must_use]
    pub const fn root_bytes(&self) -> [u8; 32] {
        self.root
    }

    /// Lowercase witness-root hex.
    #[must_use]
    pub fn root_hex(&self) -> String {
        completion_hex(&self.root)
    }

    /// Recompute both roots and verify count, join, admission, quiescence,
    /// terminal-error, and disposition semantics.
    ///
    /// # Errors
    /// Returns the first stable violated invariant.
    pub fn verify(&self) -> Result<(), TilePoolCompletionWitnessError> {
        verify_completion_witness(self)
    }

    /// Convenience integrity verdict.
    #[must_use]
    pub fn verifies_integrity(&self) -> bool {
        self.verify().is_ok()
    }

    /// Canonical semantic JSON with deterministic field order.
    ///
    /// Timing samples, steal counts, stdout/stderr, and caller publication
    /// state are deliberately absent.
    #[must_use]
    pub fn to_canonical_json(&self) -> String {
        use core::fmt::Write as _;

        let mut out = String::with_capacity(2304);
        out.push_str("{\"schema\":\"fs-exec-tilepool-completion-witness-v2\"");
        let _ = write!(out, ",\"version\":{}", self.version);
        out.push_str(",\"producer_version\":");
        push_json_string(&mut out, self.producer_version);
        let _ = write!(
            out,
            ",\"pool_placement_identity_version\":{}",
            self.pool_placement_identity_version
        );
        out.push_str(",\"pool_placement_identity\":");
        push_json_string(&mut out, &self.pool_placement_identity);
        let _ = write!(out, ",\"pool_seed\":{}", self.pool_seed);
        out.push_str(",\"kernel\":");
        push_json_string(&mut out, self.kernel);
        let _ = write!(
            out,
            ",\"kernel_id\":{},\"declared_run\":{},\"mode\":",
            self.kernel_id, self.declared_run.0
        );
        push_json_string(&mut out, self.mode);
        out.push_str(",\"scope\":{\"kind\":");
        push_json_string(&mut out, self.scope.kind);
        out.push_str(",\"parent_region_id\":");
        completion_push_optional_u64_json(&mut out, self.scope.parent_region_id);
        out.push_str(",\"parent_task_id\":");
        completion_push_optional_u64_json(&mut out, self.scope.parent_task_id);
        out.push('}');
        let _ = write!(
            out,
            ",\"planned_tiles\":{},\"plan_root\":\"{}\",\"call_replay_root\":\"{}\",\
             \"admission_completed\":{},\"admitted_tiles\":{},\"unadmitted_tiles\":{},\
             \"claimed_tiles\":{},\"completed_tiles\":{},\"break_tiles\":{},\
             \"panicked_tiles\":{},\"retained_refusal_tiles\":{},\"failed_tiles\":{},\
             \"cancelled_tiles\":{},\"planned_workers\":{},\"launched_workers\":{},\
             \"joined_workers\":{},\"worker_admission_closed\":{},\
             \"live_worker_guards_at_seal\":{},\"planned_crew_callbacks\":{},\
             \"entered_crew_callbacks\":{},\"exited_crew_callbacks\":{},\
             \"tile_scopes_opened\":{},\"live_tile_scopes_at_seal\":{},\
             \"cancellation_requested_at_entry\":{},\
             \"cancellation_requested_at_terminal\":{},\"cancellation_requested\":{},\
             \"request_phase\":\"{}\",\
             \"cancellation_observed_workers\":{},\
             \"root_metadata_bytes\":{},\"root_charge_admitted\":{},\
             \"root_charge_released\":{},\"executor_transients_quiescent\":{}",
            self.planned_tiles,
            self.plan_root_hex(),
            self.call_replay_root_hex(),
            self.admission_completed,
            self.admitted_tiles,
            self.unadmitted_tiles,
            self.claimed_tiles,
            self.completed_tiles,
            self.break_tiles,
            self.panicked_tiles,
            self.retained_refusal_tiles(),
            self.failed_tiles(),
            self.cancelled_tiles(),
            self.planned_workers,
            self.launched_workers,
            self.joined_workers,
            self.worker_admission_closed,
            self.live_worker_guards_at_seal,
            self.planned_crew_callbacks,
            self.entered_crew_callbacks,
            self.exited_crew_callbacks,
            self.tile_scopes_opened,
            self.live_tile_scopes_at_seal,
            self.cancellation_requested_at_entry,
            self.cancellation_requested_at_terminal,
            self.cancellation_requested,
            self.request_phase.name(),
            self.cancellation_observed_workers,
            self.root_metadata_bytes,
            self.root_charge_admitted,
            self.root_charge_released,
            self.executor_transients_quiescent(),
        );
        out.push_str(",\"affine_invocation_permit_root\":");
        match self.affine_invocation_permit_root {
            Some(root) => push_json_string(&mut out, &completion_hex(&root)),
            None => out.push_str("null"),
        }
        let _ = write!(
            out,
            ",\"arena_before\":{{\"live\":{},\"reserved_bytes\":{},\"free_bytes\":{},\
             \"quiescent\":{}}},\"arena_after\":{{\"live\":{},\"reserved_bytes\":{},\
             \"free_bytes\":{},\"quiescent\":{}}}",
            self.arena_before.live,
            self.arena_before.reserved_bytes,
            self.arena_before.free_bytes,
            self.arena_before.quiescent,
            self.arena_after.live,
            self.arena_after.reserved_bytes,
            self.arena_after.free_bytes,
            self.arena_after.quiescent,
        );
        completion_push_lease_json(&mut out, "lease_before", self.lease_before);
        completion_push_lease_json(&mut out, "lease_after", self.lease_after);
        out.push_str(",\"disposition\":");
        push_json_string(&mut out, self.disposition.name());
        out.push_str(",\"first_failure\":{\"kind\":");
        match self.first_failure_kind {
            Some(kind) => push_json_string(&mut out, kind),
            None => out.push_str("null"),
        }
        out.push_str(",\"tile\":");
        completion_push_optional_u64_json(&mut out, self.first_failure_tile);
        out.push('}');
        out.push_str(",\"terminal_error\":");
        match &self.terminal_error {
            Some(error) => {
                out.push_str("{\"kind\":");
                push_json_string(&mut out, completion_error_kind(error));
                out.push_str(",\"detail\":");
                push_json_string(&mut out, &error.to_string());
                out.push('}');
            }
            None => out.push_str("null"),
        }
        out.push_str(",\"no_claims\":[");
        for (index, no_claim) in self.no_claims().iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            push_json_string(&mut out, no_claim);
        }
        let _ = write!(out, "],\"witness_root\":\"{}\"}}", self.root_hex());
        out
    }

    /// Deterministic JSONL envelope for focused/e2e evidence logging.
    ///
    /// `case`, `sequence`, and `source_build_identity` are logging-envelope
    /// metadata and therefore do not alter the immutable semantic witness.
    /// stdout/stderr remains outside semantic identity.
    #[must_use]
    pub fn to_jsonl(&self, case: &str, sequence: u64, source_build_identity: &str) -> String {
        self.to_jsonl_with_reuse(case, sequence, source_build_identity, None)
    }

    /// Deterministic JSONL envelope with an optional post-run pool-reuse
    /// verdict.
    ///
    /// Reuse is necessarily observed after this immutable witness was minted,
    /// so it is explicit envelope metadata rather than executor-attested
    /// semantic identity. This method never changes [`Self::root_bytes`].
    #[must_use]
    pub fn to_jsonl_with_reuse(
        &self,
        case: &str,
        sequence: u64,
        source_build_identity: &str,
        reuse_verdict: Option<bool>,
    ) -> String {
        use core::fmt::Write as _;

        let mut out = String::with_capacity(2560);
        out.push_str("{\"schema\":\"fs-exec-tilepool-completion-e2e-v2\",\"case\":");
        push_json_string(&mut out, case);
        let _ = write!(out, ",\"sequence\":{sequence},\"source_build_identity\":");
        push_json_string(&mut out, source_build_identity);
        out.push_str(",\"reuse_verdict\":");
        match reuse_verdict {
            Some(true) => out.push_str("true"),
            Some(false) => out.push_str("false"),
            None => out.push_str("null"),
        }
        out.push_str(",\"witness\":");
        out.push_str(&self.to_canonical_json());
        out.push('}');
        out
    }
}

fn completion_invariant(name: &'static str) -> TilePoolCompletionWitnessError {
    TilePoolCompletionWitnessError::Invariant { name }
}

fn completion_bundle_invariant(name: &'static str) -> TilePoolCompletionWitnessError {
    TilePoolCompletionWitnessError::BundleInvariant { name }
}

fn completion_disposition(error: Option<&RunError>) -> TilePoolCompletionDisposition {
    match error {
        None => TilePoolCompletionDisposition::Completed,
        Some(RunError::Cancelled { .. }) => TilePoolCompletionDisposition::Cancelled,
        Some(RunError::TileFailed { .. }) => TilePoolCompletionDisposition::TileFailed,
        Some(RunError::TilePanicked { .. }) => TilePoolCompletionDisposition::TilePanicked,
        Some(RunError::WorkerSpawn { .. }) => TilePoolCompletionDisposition::WorkerSpawnRefused,
        Some(RunError::ParkedCrewBusy { .. }) => TilePoolCompletionDisposition::ParkedCrewBusy,
        Some(RunError::MemoryRefused { .. }) => TilePoolCompletionDisposition::MemoryRefused,
        Some(RunError::MemoryPlanOverflow { .. }) => {
            TilePoolCompletionDisposition::MemoryPlanOverflow
        }
        Some(RunError::MemoryAllocationRefused { .. }) => {
            TilePoolCompletionDisposition::MemoryAllocationRefused
        }
        Some(RunError::ReductionPanicked { .. }) => {
            TilePoolCompletionDisposition::ReductionPanicked
        }
        Some(RunError::Incomplete { .. }) => TilePoolCompletionDisposition::Incomplete,
    }
}

fn completion_first_failure(error: Option<&RunError>) -> (Option<&'static str>, Option<u64>) {
    match error {
        Some(RunError::TileFailed { tile, .. }) => (Some("tile-failed"), Some(*tile)),
        Some(RunError::TilePanicked { tile, .. }) => (Some("tile-panicked"), Some(*tile)),
        _ => (None, None),
    }
}

fn verify_completion_witness(
    witness: &TilePoolCompletionWitness,
) -> Result<(), TilePoolCompletionWitnessError> {
    if witness.version != TILEPOOL_COMPLETION_WITNESS_VERSION {
        return Err(TilePoolCompletionWitnessError::UnsupportedVersion {
            found: witness.version,
        });
    }
    if witness.pool_placement_identity_version != TILEPOOL_PLACEMENT_IDENTITY_VERSION {
        return Err(completion_invariant("pool-placement-identity-version"));
    }
    if !completion_placement_identity_is_well_formed(&witness.pool_placement_identity) {
        return Err(completion_invariant("pool-placement-identity-shape"));
    }
    if completion_plan_root(
        witness.pool_placement_identity_version,
        &witness.pool_placement_identity,
        witness.pool_seed,
        witness.kernel,
        witness.kernel_id,
        witness.declared_run,
        witness.mode,
        witness.planned_tiles,
    ) != witness.plan_root
    {
        return Err(TilePoolCompletionWitnessError::PlanRootMismatch);
    }
    if completion_call_replay_root(
        witness.plan_root,
        witness.scope,
        witness.affine_invocation_permit_root,
    ) != witness.call_replay_root
    {
        return Err(completion_invariant("call-replay-root"));
    }
    if completion_witness_root(witness) != witness.root {
        return Err(TilePoolCompletionWitnessError::RootMismatch);
    }
    if witness.cancellation_requested_at_entry && !witness.cancellation_requested_at_terminal {
        return Err(completion_invariant("request-entry-monotonic"));
    }
    if witness.cancellation_requested_at_terminal && !witness.cancellation_requested {
        return Err(completion_invariant("request-terminal-monotonic"));
    }
    let request_phase = completion_request_phase(
        witness.cancellation_requested_at_entry,
        witness.cancellation_requested_at_terminal,
        witness.cancellation_requested,
    )
    .ok_or_else(|| completion_invariant("request-phase-observations"))?;
    if witness.request_phase != request_phase {
        return Err(completion_invariant("derived-request-phase"));
    }
    completion_verify_arena_snapshot(witness.arena_before, "arena-before-internal")?;
    completion_verify_arena_snapshot(witness.arena_after, "arena-after-internal")?;
    completion_verify_lease_snapshot(witness.lease_before, "lease-before-internal")?;
    completion_verify_lease_snapshot(witness.lease_after, "lease-after-internal")?;
    if witness.lease_before.limit_bytes != witness.lease_after.limit_bytes {
        return Err(completion_invariant("lease-limit-stable"));
    }
    if witness.lease_after.requested_bytes < witness.lease_before.requested_bytes {
        return Err(completion_invariant("lease-requested-monotonic"));
    }
    if witness.lease_after.peak_bytes < witness.lease_before.peak_bytes {
        return Err(completion_invariant("lease-peak-monotonic"));
    }
    if witness.lease_after.refusals < witness.lease_before.refusals {
        return Err(completion_invariant("lease-refusals-monotonic"));
    }
    if witness.lease_after.release_invariant_violations
        != witness.lease_before.release_invariant_violations
    {
        return Err(completion_invariant(
            "no-run-observed-lease-release-violation",
        ));
    }
    let lease_requested_delta = witness
        .lease_after
        .requested_bytes
        .checked_sub(witness.lease_before.requested_bytes)
        .ok_or_else(|| completion_invariant("lease-requested-delta"))?;
    let lease_refusal_delta = witness
        .lease_after
        .refusals
        .checked_sub(witness.lease_before.refusals)
        .ok_or_else(|| completion_invariant("lease-refusal-delta"))?;
    if witness.root_charge_admitted
        && lease_requested_delta < u128::from(witness.root_metadata_bytes)
    {
        return Err(completion_invariant("root-charge-request-observed"));
    }
    if witness.root_charge_admitted != witness.root_charge_released {
        return Err(completion_invariant("root-charge-release"));
    }
    let parked_scope = match (
        witness.scope.kind,
        witness.scope.parent_region_id,
        witness.scope.parent_task_id,
    ) {
        ("std-thread-scope", None, None) => false,
        ("asupersync-task-scope", Some(_), Some(_)) => false,
        ("std-thread-parked-crew", None, None) => true,
        ("asupersync-task-parked-crew", Some(_), Some(_)) => true,
        _ => return Err(completion_invariant("scope-identity")),
    };
    if witness.admitted_tiles.checked_add(witness.unadmitted_tiles) != Some(witness.planned_tiles) {
        return Err(completion_invariant("plan-admission-conservation"));
    }
    if witness.admission_completed {
        if witness.admitted_tiles != witness.planned_tiles || witness.unadmitted_tiles != 0 {
            return Err(completion_invariant("admitted-plan"));
        }
        if !witness.root_charge_admitted {
            return Err(completion_invariant("root-charge-admitted"));
        }
    } else if witness.admitted_tiles != 0
        || witness.claimed_tiles != 0
        || witness.completed_tiles != 0
        || witness.break_tiles != 0
        || witness.panicked_tiles != 0
        || witness.launched_workers != 0
        || witness.joined_workers != 0
        || witness.entered_crew_callbacks != 0
        || witness.exited_crew_callbacks != 0
        || witness.tile_scopes_opened != 0
        || witness.cancellation_observed_workers != 0
    {
        return Err(completion_invariant("prelaunch-zero-work"));
    }
    if parked_scope {
        if witness.planned_crew_callbacks == 0
            || witness.planned_workers > witness.planned_crew_callbacks
        {
            return Err(completion_invariant("parked-crew-plan"));
        }
        if witness.admission_completed
            && (witness.entered_crew_callbacks != witness.planned_crew_callbacks
                || witness.exited_crew_callbacks != witness.entered_crew_callbacks)
        {
            return Err(completion_invariant("parked-crew-callback-drain"));
        }
    } else if witness.planned_crew_callbacks != 0
        || witness.entered_crew_callbacks != 0
        || witness.exited_crew_callbacks != 0
    {
        return Err(completion_invariant("nonparked-crew-callbacks"));
    }
    if witness.claimed_tiles > witness.admitted_tiles {
        return Err(completion_invariant("claimed-within-admitted"));
    }
    if witness
        .completed_tiles
        .checked_add(witness.break_tiles)
        .and_then(|count| count.checked_add(witness.panicked_tiles))
        != Some(witness.claimed_tiles)
    {
        return Err(completion_invariant("claimed-terminal-conservation"));
    }
    if witness.tile_scopes_opened != witness.claimed_tiles {
        return Err(completion_invariant("one-scope-per-claimed-tile"));
    }
    if witness.retained_refusal_tiles() > witness.break_tiles {
        return Err(completion_invariant("retained-refusal-is-break"));
    }
    if witness.launched_workers > witness.planned_workers {
        return Err(completion_invariant("launched-within-plan"));
    }
    if witness.claimed_tiles != 0 && witness.launched_workers == 0 {
        return Err(completion_invariant("claimed-work-needs-worker"));
    }
    if witness.admission_completed {
        match witness.terminal_error.as_ref() {
            Some(RunError::WorkerSpawn { .. }) => {}
            Some(RunError::Cancelled { .. }) => {
                if witness.launched_workers != 0
                    && witness.launched_workers != witness.planned_workers
                {
                    return Err(completion_invariant("cancelled-launch-set"));
                }
            }
            _ if witness.launched_workers != witness.planned_workers => {
                return Err(completion_invariant("all-planned-workers-entered"));
            }
            _ => {}
        }
    }
    if witness.joined_workers != witness.launched_workers {
        return Err(completion_invariant("all-launched-workers-joined"));
    }
    if !witness.worker_admission_closed {
        return Err(completion_invariant("worker-admission-closed"));
    }
    if witness.live_worker_guards_at_seal != 0 {
        return Err(completion_invariant("no-live-worker-guards"));
    }
    if witness.live_tile_scopes_at_seal != 0 {
        return Err(completion_invariant("no-live-tile-scopes"));
    }
    if witness.cancellation_observed_workers > witness.launched_workers {
        return Err(completion_invariant("request-observation-within-launch"));
    }
    if witness.cancellation_observed_workers != 0 && !witness.cancellation_requested_at_terminal {
        return Err(completion_invariant("observed-request-before-terminal"));
    }
    if witness.disposition != completion_disposition(witness.terminal_error.as_ref()) {
        return Err(completion_invariant("derived-disposition"));
    }
    if (witness.first_failure_kind, witness.first_failure_tile)
        != completion_first_failure(witness.terminal_error.as_ref())
    {
        return Err(completion_invariant("retained-first-failure"));
    }
    if witness
        .first_failure_tile
        .is_some_and(|tile| tile >= witness.planned_tiles)
    {
        return Err(completion_invariant("first-failure-within-plan"));
    }
    match &witness.terminal_error {
        None => {
            if !witness.admission_completed
                || witness.completed_tiles != witness.planned_tiles
                || witness.claimed_tiles != witness.planned_tiles
                || witness.break_tiles != 0
                || witness.panicked_tiles != 0
                || witness.cancellation_requested_at_terminal
            {
                return Err(completion_invariant("completed-run"));
            }
        }
        Some(RunError::Cancelled {
            kernel,
            completed,
            total,
        }) => {
            if *kernel != witness.kernel
                || *completed != witness.completed_tiles
                || *total != witness.planned_tiles
                || !witness.cancellation_requested_at_terminal
            {
                return Err(completion_invariant("cancelled-error"));
            }
        }
        Some(RunError::TilePanicked {
            kernel, completed, ..
        })
        | Some(RunError::TileFailed {
            kernel, completed, ..
        }) => {
            if *kernel != witness.kernel
                || *completed != witness.completed_tiles
                || !witness.cancellation_requested_at_terminal
                || witness.failed_tiles() == 0
            {
                return Err(completion_invariant("tile-failure-error"));
            }
        }
        Some(RunError::WorkerSpawn { kernel, worker, .. }) => {
            if *kernel != witness.kernel
                || !witness.admission_completed
                || !witness.cancellation_requested_at_terminal
                || u64::try_from(*worker).ok() != Some(witness.launched_workers)
            {
                return Err(completion_invariant("worker-spawn-error"));
            }
        }
        Some(RunError::ParkedCrewBusy { kernel }) => {
            if *kernel != witness.kernel
                || witness.admission_completed
                || !parked_scope
                || witness.root_metadata_bytes != 0
                || witness.root_charge_admitted
            {
                return Err(completion_invariant("parked-crew-busy-error"));
            }
        }
        Some(RunError::MemoryRefused { kernel, .. }) => {
            if *kernel != witness.kernel
                || witness.admission_completed
                || witness.root_charge_admitted
                || lease_refusal_delta == 0
            {
                return Err(completion_invariant("memory-refused-error"));
            }
        }
        Some(RunError::MemoryPlanOverflow { kernel, .. }) => {
            if *kernel != witness.kernel || witness.admission_completed {
                return Err(completion_invariant("memory-plan-overflow-error"));
            }
        }
        Some(RunError::MemoryAllocationRefused { kernel, .. }) => {
            if *kernel != witness.kernel
                || witness.admission_completed
                || !witness.root_charge_admitted
            {
                return Err(completion_invariant("memory-allocation-refused-error"));
            }
        }
        Some(RunError::ReductionPanicked { kernel, .. }) => {
            if *kernel != witness.kernel
                || !witness.admission_completed
                || witness.completed_tiles != witness.planned_tiles
                || witness.cancellation_requested_at_terminal
            {
                return Err(completion_invariant("reduction-error"));
            }
        }
        Some(RunError::Incomplete { kernel, tile }) => {
            if *kernel != witness.kernel
                || !witness.admission_completed
                || *tile >= witness.planned_tiles
                || witness.cancellation_requested_at_terminal
            {
                return Err(completion_invariant("incomplete-error"));
            }
        }
    }
    Ok(())
}

fn completion_placement_identity_is_well_formed(identity: &str) -> bool {
    // Producer (placement_identity_with_schema):
    //   fs-exec-tilepool-v{2}-{pinning_intent}-ccd{ccds}x{cores}-mode-{mode}-cfg-{digest}
    // The invariant must accept exactly that shape: the v2 prefix, a
    // non-empty descriptor section, the "-cfg-" separator, and the exact
    // 64-char lowercase-hex BLAKE3 digest. An earlier version of this
    // check demanded prefix+digest only, so EVERY legitimately minted
    // witness failed "pool-placement-identity-shape" (kh5tf).
    const PREFIX: &str = "fs-exec-tilepool-v2-";
    let Some(rest) = identity.strip_prefix(PREFIX) else {
        return false;
    };
    let Some((descriptor, digest)) = rest.rsplit_once("-cfg-") else {
        return false;
    };
    !descriptor.is_empty()
        && digest.len() == 64
        && digest
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn completion_request_phase(
    requested_at_entry: bool,
    requested_at_terminal: bool,
    requested_at_seal: bool,
) -> Option<TilePoolRequestPhase> {
    match (requested_at_entry, requested_at_terminal, requested_at_seal) {
        (false, false, false) => Some(TilePoolRequestPhase::NotRequested),
        (true, true, true) => Some(TilePoolRequestPhase::BeforeEntry),
        (false, true, true) => Some(TilePoolRequestPhase::BeforeTerminalDecision),
        (false, false, true) => Some(TilePoolRequestPhase::AfterTerminalDecision),
        _ => None,
    }
}

fn completion_verify_arena_snapshot(
    snapshot: CompletionArenaSnapshot,
    invariant: &'static str,
) -> Result<(), TilePoolCompletionWitnessError> {
    if snapshot.free_bytes > snapshot.reserved_bytes
        || snapshot.quiescent
            != (snapshot.live == 0 && snapshot.reserved_bytes == snapshot.free_bytes)
    {
        return Err(completion_invariant(invariant));
    }
    Ok(())
}

fn completion_verify_lease_snapshot(
    snapshot: CompletionLeaseSnapshot,
    invariant: &'static str,
) -> Result<(), TilePoolCompletionWitnessError> {
    if snapshot.used_bytes > snapshot.peak_bytes
        || snapshot
            .limit_bytes
            .is_some_and(|limit| snapshot.used_bytes > limit || snapshot.peak_bytes > limit)
    {
        return Err(completion_invariant(invariant));
    }
    Ok(())
}

fn completion_hex(bytes: &[u8; 32]) -> String {
    use core::fmt::Write as _;

    let mut out = String::with_capacity(64);
    for byte in bytes {
        let _ = write!(out, "{byte:02x}");
    }
    out
}

fn completion_push_optional_u64_json(out: &mut String, value: Option<u64>) {
    use core::fmt::Write as _;

    match value {
        Some(value) => {
            let _ = write!(out, "{value}");
        }
        None => out.push_str("null"),
    }
}

fn completion_push_lease_json(out: &mut String, name: &str, snapshot: CompletionLeaseSnapshot) {
    use core::fmt::Write as _;

    out.push_str(",\"");
    out.push_str(name);
    out.push_str("\":{\"limit_bytes\":");
    match snapshot.limit_bytes {
        Some(limit) => {
            let _ = write!(out, "{limit}");
        }
        None => out.push_str("null"),
    }
    let _ = write!(
        out,
        ",\"requested_bytes\":{},\"peak_bytes\":{},\"used_bytes\":{},\
         \"refusals\":{},\"release_invariant_violations\":{}}}",
        snapshot.requested_bytes,
        snapshot.peak_bytes,
        snapshot.used_bytes,
        snapshot.refusals,
        snapshot.release_invariant_violations,
    );
}

fn completion_error_kind(error: &RunError) -> &'static str {
    match error {
        RunError::Cancelled { .. } => "cancelled",
        RunError::TilePanicked { .. } => "tile-panicked",
        RunError::TileFailed { .. } => "tile-failed",
        RunError::WorkerSpawn { .. } => "worker-spawn-refused",
        RunError::ParkedCrewBusy { .. } => "parked-crew-busy",
        RunError::MemoryRefused { .. } => "memory-refused",
        RunError::MemoryPlanOverflow { .. } => "memory-plan-overflow",
        RunError::MemoryAllocationRefused { .. } => "memory-allocation-refused",
        RunError::ReductionPanicked { .. } => "reduction-panicked",
        RunError::Incomplete { .. } => "incomplete",
    }
}

fn completion_hash_field(hasher: &mut fs_blake3::DomainHasher, name: &'static str, value: &[u8]) {
    let name_len = u64::try_from(name.len()).unwrap_or(u64::MAX);
    let value_len = u64::try_from(value.len()).unwrap_or(u64::MAX);
    hasher.update(&name_len.to_le_bytes());
    hasher.update(name.as_bytes());
    hasher.update(&value_len.to_le_bytes());
    hasher.update(value);
}

fn completion_hash_bool(hasher: &mut fs_blake3::DomainHasher, name: &'static str, value: bool) {
    completion_hash_field(hasher, name, &[u8::from(value)]);
}

fn completion_hash_u64(hasher: &mut fs_blake3::DomainHasher, name: &'static str, value: u64) {
    completion_hash_field(hasher, name, &value.to_le_bytes());
}

fn completion_hash_u128(hasher: &mut fs_blake3::DomainHasher, name: &'static str, value: u128) {
    completion_hash_field(hasher, name, &value.to_le_bytes());
}

fn completion_hash_usize(hasher: &mut fs_blake3::DomainHasher, name: &'static str, value: usize) {
    completion_hash_u64(hasher, name, u64::try_from(value).unwrap_or(u64::MAX));
}

fn completion_hash_optional_u64(
    hasher: &mut fs_blake3::DomainHasher,
    presence_name: &'static str,
    value_name: &'static str,
    value: Option<u64>,
) {
    completion_hash_bool(hasher, presence_name, value.is_some());
    if let Some(value) = value {
        completion_hash_u64(hasher, value_name, value);
    }
}

fn completion_hash_optional_root(
    hasher: &mut fs_blake3::DomainHasher,
    presence_name: &'static str,
    value_name: &'static str,
    value: Option<[u8; 32]>,
) {
    completion_hash_bool(hasher, presence_name, value.is_some());
    if let Some(value) = value {
        completion_hash_field(hasher, value_name, &value);
    }
}

fn completion_run_error_tag(error: &RunError) -> u64 {
    match error {
        RunError::Cancelled { .. } => 0,
        RunError::TilePanicked { .. } => 1,
        RunError::TileFailed { .. } => 2,
        RunError::WorkerSpawn { .. } => 3,
        RunError::ParkedCrewBusy { .. } => 4,
        RunError::MemoryRefused { .. } => 5,
        RunError::MemoryPlanOverflow { .. } => 6,
        RunError::MemoryAllocationRefused { .. } => 7,
        RunError::ReductionPanicked { .. } => 8,
        RunError::Incomplete { .. } => 9,
    }
}

fn completion_hash_run_error(hasher: &mut fs_blake3::DomainHasher, error: &RunError) {
    completion_hash_u64(
        hasher,
        "terminal-error.tag",
        completion_run_error_tag(error),
    );
    match error {
        RunError::Cancelled {
            kernel,
            completed,
            total,
        } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_u64(hasher, "terminal-error.completed", *completed);
            completion_hash_u64(hasher, "terminal-error.total", *total);
        }
        RunError::TilePanicked {
            kernel,
            tile,
            message,
            completed,
        } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_u64(hasher, "terminal-error.tile", *tile);
            completion_hash_field(hasher, "terminal-error.message", message.as_bytes());
            completion_hash_u64(hasher, "terminal-error.completed", *completed);
        }
        RunError::TileFailed {
            kernel,
            tile,
            failure,
            completed,
        } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_u64(hasher, "terminal-error.tile", *tile);
            completion_hash_tile_failure(hasher, failure);
            completion_hash_u64(hasher, "terminal-error.completed", *completed);
        }
        RunError::WorkerSpawn {
            kernel,
            worker,
            message,
        } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_u64(
                hasher,
                "terminal-error.worker",
                u64::try_from(*worker).unwrap_or(u64::MAX),
            );
            completion_hash_field(hasher, "terminal-error.message", message.as_bytes());
        }
        RunError::ParkedCrewBusy { kernel } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
        }
        RunError::MemoryRefused {
            kernel,
            what,
            requested_bytes,
            used_bytes,
            limit_bytes,
        } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_field(hasher, "terminal-error.what", what.as_bytes());
            completion_hash_u64(hasher, "terminal-error.requested-bytes", *requested_bytes);
            completion_hash_u64(hasher, "terminal-error.used-bytes", *used_bytes);
            completion_hash_u64(hasher, "terminal-error.limit-bytes", *limit_bytes);
        }
        RunError::MemoryPlanOverflow { kernel, what } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_field(hasher, "terminal-error.what", what.as_bytes());
        }
        RunError::MemoryAllocationRefused {
            kernel,
            what,
            requested_bytes,
        } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_field(hasher, "terminal-error.what", what.as_bytes());
            completion_hash_u64(hasher, "terminal-error.requested-bytes", *requested_bytes);
        }
        RunError::ReductionPanicked { kernel, message } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_field(hasher, "terminal-error.message", message.as_bytes());
        }
        RunError::Incomplete { kernel, tile } => {
            completion_hash_field(hasher, "terminal-error.kernel", kernel.as_bytes());
            completion_hash_u64(hasher, "terminal-error.tile", *tile);
        }
    }
}

fn completion_hash_tile_failure(hasher: &mut fs_blake3::DomainHasher, failure: &TileFailure) {
    match failure {
        TileFailure::Allocation(error) => {
            completion_hash_u64(hasher, "terminal-error.failure.tag", 0);
            completion_hash_alloc_error(hasher, error);
        }
        TileFailure::InjectedFault {
            plan_version,
            plan_seed,
            tiles,
            touches_per_tile,
            touch,
        } => {
            completion_hash_u64(hasher, "terminal-error.failure.tag", 1);
            completion_hash_u64(
                hasher,
                "terminal-error.failure.plan-version",
                u64::from(*plan_version),
            );
            completion_hash_u64(hasher, "terminal-error.failure.plan-seed", *plan_seed);
            completion_hash_u64(hasher, "terminal-error.failure.tiles", *tiles);
            completion_hash_u64(
                hasher,
                "terminal-error.failure.touches-per-tile",
                u64::from(*touches_per_tile),
            );
            completion_hash_u64(hasher, "terminal-error.failure.touch", u64::from(*touch));
        }
    }
}

#[allow(clippy::too_many_lines)]
fn completion_hash_alloc_error(hasher: &mut fs_blake3::DomainHasher, error: &fs_alloc::AllocError) {
    match error {
        fs_alloc::AllocError::Exhausted {
            site,
            requested_bytes,
            reserved_bytes,
            limit_bytes,
        } => {
            completion_hash_u64(hasher, "terminal-error.failure.alloc.tag", 0);
            completion_hash_field(hasher, "terminal-error.failure.alloc.site", site.as_bytes());
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.requested-bytes",
                *requested_bytes,
            );
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.reserved-bytes",
                *reserved_bytes,
            );
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.limit-bytes",
                *limit_bytes,
            );
        }
        fs_alloc::AllocError::OutOfMemory {
            site,
            requested_bytes,
        } => {
            completion_hash_u64(hasher, "terminal-error.failure.alloc.tag", 1);
            completion_hash_field(hasher, "terminal-error.failure.alloc.site", site.as_bytes());
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.requested-bytes",
                *requested_bytes,
            );
        }
        fs_alloc::AllocError::LeaseExhausted {
            site,
            requested_bytes,
            used_bytes,
            limit_bytes,
        } => {
            completion_hash_u64(hasher, "terminal-error.failure.alloc.tag", 2);
            completion_hash_field(hasher, "terminal-error.failure.alloc.site", site.as_bytes());
            completion_hash_u64(
                hasher,
                "terminal-error.failure.alloc.requested-bytes",
                *requested_bytes,
            );
            completion_hash_u64(
                hasher,
                "terminal-error.failure.alloc.used-bytes",
                *used_bytes,
            );
            completion_hash_u64(
                hasher,
                "terminal-error.failure.alloc.limit-bytes",
                *limit_bytes,
            );
        }
        fs_alloc::AllocError::LayoutOverflow {
            site,
            len,
            elem_bytes,
        } => {
            completion_hash_u64(hasher, "terminal-error.failure.alloc.tag", 3);
            completion_hash_field(hasher, "terminal-error.failure.alloc.site", site.as_bytes());
            completion_hash_usize(hasher, "terminal-error.failure.alloc.len", *len);
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.elem-bytes",
                *elem_bytes,
            );
        }
        fs_alloc::AllocError::ReservationOverflow {
            site,
            base_bytes,
            additional_bytes,
        } => {
            completion_hash_u64(hasher, "terminal-error.failure.alloc.tag", 4);
            completion_hash_field(hasher, "terminal-error.failure.alloc.site", site.as_bytes());
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.base-bytes",
                *base_bytes,
            );
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.additional-bytes",
                *additional_bytes,
            );
        }
        fs_alloc::AllocError::ReclaimedChunkCorrupted {
            site,
            poison_version,
            poison_seed,
            chunk_bytes,
            offset,
            expected,
            actual,
        } => {
            completion_hash_u64(hasher, "terminal-error.failure.alloc.tag", 5);
            completion_hash_field(hasher, "terminal-error.failure.alloc.site", site.as_bytes());
            completion_hash_u64(
                hasher,
                "terminal-error.failure.alloc.poison-version",
                u64::from(*poison_version),
            );
            completion_hash_u64(
                hasher,
                "terminal-error.failure.alloc.poison-seed",
                *poison_seed,
            );
            completion_hash_usize(
                hasher,
                "terminal-error.failure.alloc.chunk-bytes",
                *chunk_bytes,
            );
            completion_hash_usize(hasher, "terminal-error.failure.alloc.offset", *offset);
            completion_hash_u64(
                hasher,
                "terminal-error.failure.alloc.expected",
                u64::from(*expected),
            );
            completion_hash_u64(
                hasher,
                "terminal-error.failure.alloc.actual",
                u64::from(*actual),
            );
        }
    }
}

fn completion_plan_root(
    pool_placement_identity_version: u32,
    pool_placement_identity: &str,
    pool_seed: u64,
    kernel: &'static str,
    kernel_id: u64,
    run: RunId,
    mode: &'static str,
    planned_tiles: u64,
) -> [u8; 32] {
    let mut hasher = fs_blake3::DomainHasher::new(TILEPOOL_COMPLETION_PLAN_DOMAIN);
    completion_hash_u64(
        &mut hasher,
        "pool-placement-identity-version",
        u64::from(pool_placement_identity_version),
    );
    completion_hash_field(
        &mut hasher,
        "pool-placement-identity",
        pool_placement_identity.as_bytes(),
    );
    completion_hash_u64(&mut hasher, "pool-seed", pool_seed);
    completion_hash_field(&mut hasher, "kernel", kernel.as_bytes());
    completion_hash_u64(&mut hasher, "kernel-id", kernel_id);
    completion_hash_u64(&mut hasher, "declared-run", run.0);
    completion_hash_field(&mut hasher, "mode", mode.as_bytes());
    completion_hash_u64(&mut hasher, "planned-tiles", planned_tiles);
    *hasher.finalize().as_bytes()
}

fn completion_call_replay_root(
    plan_root: [u8; 32],
    scope: CompletionScopeIdentity,
    affine_invocation_permit_root: Option<[u8; 32]>,
) -> [u8; 32] {
    let mut hasher = fs_blake3::DomainHasher::new(TILEPOOL_COMPLETION_CALL_REPLAY_DOMAIN);
    completion_hash_field(&mut hasher, "plan-root", &plan_root);
    completion_hash_field(&mut hasher, "scope.kind", scope.kind.as_bytes());
    completion_hash_optional_u64(
        &mut hasher,
        "scope.parent-region-id.present",
        "scope.parent-region-id",
        scope.parent_region_id,
    );
    completion_hash_optional_u64(
        &mut hasher,
        "scope.parent-task-id.present",
        "scope.parent-task-id",
        scope.parent_task_id,
    );
    completion_hash_optional_root(
        &mut hasher,
        "affine-invocation-permit-root.present",
        "affine-invocation-permit-root",
        affine_invocation_permit_root,
    );
    *hasher.finalize().as_bytes()
}

fn completion_witness_root(witness: &TilePoolCompletionWitness) -> [u8; 32] {
    let mut hasher = fs_blake3::DomainHasher::new(TILEPOOL_COMPLETION_WITNESS_DOMAIN);
    completion_hash_u64(&mut hasher, "version", u64::from(witness.version));
    completion_hash_field(
        &mut hasher,
        "producer-version",
        witness.producer_version.as_bytes(),
    );
    completion_hash_u64(
        &mut hasher,
        "pool-placement-identity-version",
        u64::from(witness.pool_placement_identity_version),
    );
    completion_hash_field(
        &mut hasher,
        "pool-placement-identity",
        witness.pool_placement_identity.as_bytes(),
    );
    completion_hash_u64(&mut hasher, "pool-seed", witness.pool_seed);
    completion_hash_field(&mut hasher, "kernel", witness.kernel.as_bytes());
    completion_hash_u64(&mut hasher, "kernel-id", witness.kernel_id);
    completion_hash_u64(&mut hasher, "declared-run", witness.declared_run.0);
    completion_hash_field(&mut hasher, "mode", witness.mode.as_bytes());
    completion_hash_field(&mut hasher, "scope.kind", witness.scope.kind.as_bytes());
    completion_hash_optional_u64(
        &mut hasher,
        "scope.parent-region-id.present",
        "scope.parent-region-id",
        witness.scope.parent_region_id,
    );
    completion_hash_optional_u64(
        &mut hasher,
        "scope.parent-task-id.present",
        "scope.parent-task-id",
        witness.scope.parent_task_id,
    );
    completion_hash_u64(&mut hasher, "planned-tiles", witness.planned_tiles);
    completion_hash_field(&mut hasher, "plan-root", &witness.plan_root);
    completion_hash_field(&mut hasher, "call-replay-root", &witness.call_replay_root);
    completion_hash_optional_root(
        &mut hasher,
        "affine-invocation-permit-root.present",
        "affine-invocation-permit-root",
        witness.affine_invocation_permit_root,
    );
    completion_hash_bool(
        &mut hasher,
        "admission-completed",
        witness.admission_completed,
    );
    completion_hash_u64(&mut hasher, "admitted-tiles", witness.admitted_tiles);
    completion_hash_u64(&mut hasher, "unadmitted-tiles", witness.unadmitted_tiles);
    completion_hash_u64(&mut hasher, "claimed-tiles", witness.claimed_tiles);
    completion_hash_u64(&mut hasher, "completed-tiles", witness.completed_tiles);
    completion_hash_u64(&mut hasher, "break-tiles", witness.break_tiles);
    completion_hash_u64(&mut hasher, "panicked-tiles", witness.panicked_tiles);
    completion_hash_u64(&mut hasher, "planned-workers", witness.planned_workers);
    completion_hash_u64(&mut hasher, "launched-workers", witness.launched_workers);
    completion_hash_u64(&mut hasher, "joined-workers", witness.joined_workers);
    completion_hash_bool(
        &mut hasher,
        "worker-admission-closed",
        witness.worker_admission_closed,
    );
    completion_hash_u64(
        &mut hasher,
        "live-worker-guards-at-seal",
        witness.live_worker_guards_at_seal,
    );
    completion_hash_u64(
        &mut hasher,
        "planned-crew-callbacks",
        witness.planned_crew_callbacks,
    );
    completion_hash_u64(
        &mut hasher,
        "entered-crew-callbacks",
        witness.entered_crew_callbacks,
    );
    completion_hash_u64(
        &mut hasher,
        "exited-crew-callbacks",
        witness.exited_crew_callbacks,
    );
    completion_hash_u64(
        &mut hasher,
        "tile-scopes-opened",
        witness.tile_scopes_opened,
    );
    completion_hash_u64(
        &mut hasher,
        "live-tile-scopes-at-seal",
        witness.live_tile_scopes_at_seal,
    );
    completion_hash_bool(
        &mut hasher,
        "cancellation-requested-at-entry",
        witness.cancellation_requested_at_entry,
    );
    completion_hash_bool(
        &mut hasher,
        "cancellation-requested-at-terminal",
        witness.cancellation_requested_at_terminal,
    );
    completion_hash_bool(
        &mut hasher,
        "cancellation-requested",
        witness.cancellation_requested,
    );
    completion_hash_u64(&mut hasher, "request-phase", witness.request_phase.tag());
    completion_hash_u64(
        &mut hasher,
        "cancellation-observed-workers",
        witness.cancellation_observed_workers,
    );
    completion_hash_u64(
        &mut hasher,
        "root-metadata-bytes",
        witness.root_metadata_bytes,
    );
    completion_hash_bool(
        &mut hasher,
        "root-charge-admitted",
        witness.root_charge_admitted,
    );
    completion_hash_bool(
        &mut hasher,
        "root-charge-released",
        witness.root_charge_released,
    );
    completion_hash_u64(&mut hasher, "arena-before.live", witness.arena_before.live);
    completion_hash_u64(
        &mut hasher,
        "arena-before.reserved-bytes",
        witness.arena_before.reserved_bytes,
    );
    completion_hash_u64(
        &mut hasher,
        "arena-before.free-bytes",
        witness.arena_before.free_bytes,
    );
    completion_hash_bool(
        &mut hasher,
        "arena-before.quiescent",
        witness.arena_before.quiescent,
    );
    completion_hash_u64(&mut hasher, "arena-after.live", witness.arena_after.live);
    completion_hash_u64(
        &mut hasher,
        "arena-after.reserved-bytes",
        witness.arena_after.reserved_bytes,
    );
    completion_hash_u64(
        &mut hasher,
        "arena-after.free-bytes",
        witness.arena_after.free_bytes,
    );
    completion_hash_bool(
        &mut hasher,
        "arena-after.quiescent",
        witness.arena_after.quiescent,
    );
    completion_hash_lease_before(&mut hasher, witness.lease_before);
    completion_hash_lease_after(&mut hasher, witness.lease_after);
    completion_hash_u64(&mut hasher, "disposition", witness.disposition.tag());
    completion_hash_bool(
        &mut hasher,
        "first-failure-kind-present",
        witness.first_failure_kind.is_some(),
    );
    if let Some(kind) = witness.first_failure_kind {
        completion_hash_field(&mut hasher, "first-failure-kind", kind.as_bytes());
    }
    completion_hash_optional_u64(
        &mut hasher,
        "first-failure-tile-present",
        "first-failure-tile",
        witness.first_failure_tile,
    );
    completion_hash_bool(
        &mut hasher,
        "terminal-error-present",
        witness.terminal_error.is_some(),
    );
    if let Some(error) = &witness.terminal_error {
        completion_hash_run_error(&mut hasher, error);
    }
    completion_hash_u64(
        &mut hasher,
        "no-claim-count",
        u64::try_from(witness.no_claims().len()).unwrap_or(u64::MAX),
    );
    for no_claim in witness.no_claims() {
        completion_hash_field(&mut hasher, "no-claim", no_claim.as_bytes());
    }
    *hasher.finalize().as_bytes()
}

fn completion_hash_lease_before(
    hasher: &mut fs_blake3::DomainHasher,
    lease: CompletionLeaseSnapshot,
) {
    completion_hash_optional_u64(
        hasher,
        "lease-before.limit-bytes.present",
        "lease-before.limit-bytes",
        lease.limit_bytes,
    );
    completion_hash_u128(
        hasher,
        "lease-before.requested-bytes",
        lease.requested_bytes,
    );
    completion_hash_u64(hasher, "lease-before.peak-bytes", lease.peak_bytes);
    completion_hash_u64(hasher, "lease-before.used-bytes", lease.used_bytes);
    completion_hash_u128(hasher, "lease-before.refusals", lease.refusals);
    completion_hash_u128(
        hasher,
        "lease-before.release-invariant-violations",
        lease.release_invariant_violations,
    );
}

fn completion_hash_lease_after(
    hasher: &mut fs_blake3::DomainHasher,
    lease: CompletionLeaseSnapshot,
) {
    completion_hash_optional_u64(
        hasher,
        "lease-after.limit-bytes.present",
        "lease-after.limit-bytes",
        lease.limit_bytes,
    );
    completion_hash_u128(hasher, "lease-after.requested-bytes", lease.requested_bytes);
    completion_hash_u64(hasher, "lease-after.peak-bytes", lease.peak_bytes);
    completion_hash_u64(hasher, "lease-after.used-bytes", lease.used_bytes);
    completion_hash_u128(hasher, "lease-after.refusals", lease.refusals);
    completion_hash_u128(
        hasher,
        "lease-after.release-invariant-violations",
        lease.release_invariant_violations,
    );
}

struct CompletionWitnessFields {
    pool_placement_identity: String,
    pool_seed: u64,
    affine_invocation_permit_root: Option<[u8; 32]>,
    kernel: &'static str,
    kernel_id: u64,
    declared_run: RunId,
    mode: &'static str,
    scope: CompletionScopeIdentity,
    planned_tiles: u64,
    admission_completed: bool,
    admitted_tiles: u64,
    claimed_tiles: u64,
    completed_tiles: u64,
    break_tiles: u64,
    panicked_tiles: u64,
    planned_workers: u64,
    launched_workers: u64,
    joined_workers: u64,
    planned_crew_callbacks: u64,
    entered_crew_callbacks: u64,
    exited_crew_callbacks: u64,
    tile_scopes_opened: u64,
    live_tile_scopes_at_seal: u64,
    cancellation_requested_at_entry: bool,
    cancellation_requested_at_terminal: bool,
    cancellation_requested: bool,
    cancellation_observed_workers: u64,
    root_metadata_bytes: u64,
    root_charge_admitted: bool,
    root_charge_released: bool,
    arena_before: CompletionArenaSnapshot,
    arena_after: CompletionArenaSnapshot,
    lease_before: CompletionLeaseSnapshot,
    lease_after: CompletionLeaseSnapshot,
    terminal_error: Option<RunError>,
}

fn mint_completion_witness(
    fields: CompletionWitnessFields,
) -> Result<TilePoolCompletionWitness, TilePoolCompletionWitnessError> {
    let plan_root = completion_plan_root(
        TILEPOOL_PLACEMENT_IDENTITY_VERSION,
        &fields.pool_placement_identity,
        fields.pool_seed,
        fields.kernel,
        fields.kernel_id,
        fields.declared_run,
        fields.mode,
        fields.planned_tiles,
    );
    let affine_invocation_permit_root = fields.affine_invocation_permit_root;
    let call_replay_root =
        completion_call_replay_root(plan_root, fields.scope, affine_invocation_permit_root);
    let request_phase = completion_request_phase(
        fields.cancellation_requested_at_entry,
        fields.cancellation_requested_at_terminal,
        fields.cancellation_requested,
    )
    .ok_or_else(|| completion_invariant("mint-request-phase-observations"))?;
    let admission_completed = fields.admission_completed;
    let admitted_tiles = fields.admitted_tiles;
    let launched_workers = fields.launched_workers;
    let joined_workers = fields.joined_workers;
    let (first_failure_kind, first_failure_tile) =
        completion_first_failure(fields.terminal_error.as_ref());
    let mut witness = TilePoolCompletionWitness {
        version: TILEPOOL_COMPLETION_WITNESS_VERSION,
        producer_version: env!("CARGO_PKG_VERSION"),
        pool_placement_identity_version: TILEPOOL_PLACEMENT_IDENTITY_VERSION,
        pool_placement_identity: fields.pool_placement_identity,
        pool_seed: fields.pool_seed,
        kernel: fields.kernel,
        kernel_id: fields.kernel_id,
        declared_run: fields.declared_run,
        mode: fields.mode,
        scope: fields.scope,
        planned_tiles: fields.planned_tiles,
        plan_root,
        call_replay_root,
        affine_invocation_permit_root,
        admission_completed,
        admitted_tiles,
        unadmitted_tiles: fields.planned_tiles.saturating_sub(admitted_tiles),
        claimed_tiles: fields.claimed_tiles,
        completed_tiles: fields.completed_tiles,
        break_tiles: fields.break_tiles,
        panicked_tiles: fields.panicked_tiles,
        planned_workers: fields.planned_workers,
        launched_workers,
        joined_workers,
        worker_admission_closed: true,
        live_worker_guards_at_seal: launched_workers.saturating_sub(joined_workers),
        planned_crew_callbacks: fields.planned_crew_callbacks,
        entered_crew_callbacks: fields.entered_crew_callbacks,
        exited_crew_callbacks: fields.exited_crew_callbacks,
        tile_scopes_opened: fields.tile_scopes_opened,
        live_tile_scopes_at_seal: fields.live_tile_scopes_at_seal,
        cancellation_requested_at_entry: fields.cancellation_requested_at_entry,
        cancellation_requested_at_terminal: fields.cancellation_requested_at_terminal,
        cancellation_requested: fields.cancellation_requested,
        request_phase,
        cancellation_observed_workers: fields.cancellation_observed_workers,
        root_metadata_bytes: fields.root_metadata_bytes,
        root_charge_admitted: fields.root_charge_admitted,
        root_charge_released: fields.root_charge_released,
        arena_before: fields.arena_before,
        arena_after: fields.arena_after,
        lease_before: fields.lease_before,
        lease_after: fields.lease_after,
        disposition: completion_disposition(fields.terminal_error.as_ref()),
        first_failure_kind,
        first_failure_tile,
        terminal_error: fields.terminal_error,
        root: [0; 32],
    };
    witness.root = completion_witness_root(&witness);
    witness.verify()?;
    Ok(witness)
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

#[cfg(test)]
thread_local! {
    /// One-shot, thread-local backing-allocation refusal used only to prove
    /// prelaunch sealing and lease-charge rollback. Root allocation happens
    /// on the calling thread, so parallel tests cannot consume each other's
    /// injected refusal.
    static TEST_ROOT_ALLOCATION_REFUSAL: std::cell::Cell<Option<&'static str>> =
        const { std::cell::Cell::new(None) };

    /// One-shot scoped-thread spawn refusal. The hook is consumed on the
    /// caller's launch thread before attempting the selected worker, while
    /// every earlier real worker still traverses the normal join/drain path.
    static TEST_WORKER_SPAWN_REFUSAL: std::cell::Cell<Option<usize>> =
        const { std::cell::Cell::new(None) };
}

#[cfg(test)]
fn refuse_next_root_allocation(what: &'static str) {
    TEST_ROOT_ALLOCATION_REFUSAL.with(|refusal| {
        assert!(
            refusal.replace(Some(what)).is_none(),
            "root allocation refusal already armed on this test thread"
        );
    });
}

#[cfg(test)]
fn refuse_next_worker_spawn(worker: usize) {
    TEST_WORKER_SPAWN_REFUSAL.with(|refusal| {
        assert!(
            refusal.replace(Some(worker)).is_none(),
            "worker spawn refusal already armed on this test thread"
        );
    });
}

#[cfg(test)]
fn take_worker_spawn_refusal(worker: usize) -> bool {
    TEST_WORKER_SPAWN_REFUSAL.with(|refusal| {
        if refusal.get() == Some(worker) {
            refusal.set(None);
            true
        } else {
            false
        }
    })
}

fn try_reserve_root_vec<T>(
    values: &mut Vec<T>,
    capacity: usize,
    kernel: &'static str,
    what: &'static str,
) -> Result<(), RunError> {
    #[cfg(test)]
    if TEST_ROOT_ALLOCATION_REFUSAL.with(|refusal| {
        if refusal.get() == Some(what) {
            refusal.set(None);
            true
        } else {
            false
        }
    }) {
        return Err(RunError::MemoryAllocationRefused {
            kernel,
            what,
            requested_bytes: allocation_bytes::<T>(capacity),
        });
    }

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
    Crew {
        crew: &'a crate::crew::Crew<Caps>,
        scope: CompletionScopeIdentity,
        dispatch_admission: &'a AtomicBool,
    },
}

impl<'a, Caps: 'static> Launch<'a, Caps> {
    fn completion_scope(self) -> CompletionScopeIdentity {
        match self {
            Self::OwnScope => CompletionScopeIdentity::STD_SCOPED,
            Self::TaskScope(cx) => CompletionScopeIdentity::task_scoped(cx),
            Self::Crew { scope, .. } => scope,
        }
    }

    fn planned_crew_callbacks(self) -> u64 {
        match self {
            Self::OwnScope | Self::TaskScope(_) => 0,
            Self::Crew { crew, .. } => u64::try_from(crew.workers()).unwrap_or(u64::MAX),
        }
    }

    fn dispatch_admission(self) -> Option<&'a AtomicBool> {
        match self {
            Self::OwnScope | Self::TaskScope(_) => None,
            Self::Crew {
                dispatch_admission, ..
            } => Some(dispatch_admission),
        }
    }
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

struct ParkedDispatchAdmission<'a> {
    flag: &'a AtomicBool,
}

impl<'a> ParkedDispatchAdmission<'a> {
    fn try_acquire(flag: &'a AtomicBool) -> Option<Self> {
        flag.compare_exchange(false, true, Ordering::AcqRel, Ordering::Acquire)
            .ok()
            .map(|_| Self { flag })
    }
}

impl Drop for ParkedDispatchAdmission<'_> {
    fn drop(&mut self) {
        self.flag.store(false, Ordering::Release);
    }
}

/// Capability set used by locally parked crews that carry no ambient
/// asupersync task. Exported so higher-level scoped runners can name the
/// callback-only [`ParkedTilePool`] type without depending on asupersync's
/// concrete capability module.
pub type LocalTaskCaps = asupersync::cx::cap::All;

type NoTask = LocalTaskCaps;

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
    claimed: &'a AtomicU64,
    breaks: &'a AtomicU64,
    panics: &'a AtomicU64,
    workers_entered: &'a AtomicU64,
    workers_exited: &'a AtomicU64,
    tile_scopes_opened: &'a AtomicU64,
    tile_scopes_live: &'a AtomicU64,
    cancellation_observed_workers: &'a AtomicU64,
    steals: &'a AtomicU64,
    cross_steals: &'a AtomicU64,
    panic_box: &'a Mutex<Option<(u64, String)>>,
    refusal_sink: &'a RefusalSink,
}

struct WorkerCompletionGuard<'a> {
    exited: Option<&'a AtomicU64>,
}

struct CrewCallbackCompletionGuard<'a> {
    exited: Option<&'a AtomicU64>,
}

impl Drop for CrewCallbackCompletionGuard<'_> {
    fn drop(&mut self) {
        if let Some(exited) = self.exited {
            exited.fetch_add(1, Ordering::Release);
        }
    }
}

impl Drop for WorkerCompletionGuard<'_> {
    fn drop(&mut self) {
        if let Some(exited) = self.exited {
            exited.fetch_add(1, Ordering::Release);
        }
    }
}

struct TileScopeCompletionGuard<'a> {
    live: Option<&'a AtomicU64>,
}

impl Drop for TileScopeCompletionGuard<'_> {
    fn drop(&mut self) {
        if let Some(live) = self.live {
            live.fetch_sub(1, Ordering::Release);
        }
    }
}

/// The worker protocol, shared verbatim by both launch harnesses. When
/// `task_cx` is present (the asupersync lane), the CALLING task's
/// cancellation and budget bound the run: every tile boundary checkpoints
/// the task context, and a failed checkpoint converts into a gate request
/// so the pool's normal drain protocol — including its cancel-latency
/// histogram — applies unchanged (P7: one drain semantics, two signals).
fn worker_loop<const COMPLETION: bool, Caps, K: TileKernel>(
    ctx: &WorkerCtx<'_, K>,
    w: usize,
    task_cx: Option<&CpuCx<Caps>>,
) {
    if COMPLETION {
        ctx.workers_entered.fetch_add(1, Ordering::Release);
    }
    let _worker_completion = WorkerCompletionGuard {
        exited: COMPLETION.then_some(ctx.workers_exited),
    };
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
            if COMPLETION {
                ctx.cancellation_observed_workers
                    .fetch_add(1, Ordering::Relaxed);
            }
            if let Some(observed_at_ns) = ctx.gate.latency_now_ns() {
                let _ = ctx.observed[w].get().compare_exchange(
                    0,
                    observed_at_ns.max(1),
                    Ordering::AcqRel,
                    Ordering::Acquire,
                );
            }
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
        if COMPLETION {
            ctx.claimed.fetch_add(1, Ordering::Relaxed);
        }
        let key = StreamKey {
            seed: ctx.config.seed,
            kernel_id: ctx.kernel_id,
            tile,
            iteration: ctx.iteration,
        };
        // Every tile arena charges the shared operation
        // lease while its chunks are held (bead wf9.16).
        if COMPLETION {
            ctx.tile_scopes_opened.fetch_add(1, Ordering::Relaxed);
            ctx.tile_scopes_live.fetch_add(1, Ordering::Release);
        }
        let tile_scope_completion = TileScopeCompletionGuard {
            live: COMPLETION.then_some(ctx.tile_scopes_live),
        };
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
        drop(tile_scope_completion);
        match outcome {
            Ok(ControlFlow::Continue(out)) => {
                *ctx.slots[tile as usize].lock().expect("slot") = Some(out);
                ctx.done_by[w].get().fetch_add(1, Ordering::Relaxed);
            }
            Ok(ControlFlow::Break(_cancelled)) => {
                if COMPLETION {
                    ctx.breaks.fetch_add(1, Ordering::Relaxed);
                }
                // Kernel observed the gate (or self-cancelled):
                // make it global and drain.
                ctx.gate.request();
            }
            Err(payload) => {
                if COMPLETION {
                    ctx.panics.fetch_add(1, Ordering::Relaxed);
                }
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

    /// [`Self::run`] plus executor-minted completion evidence from the same
    /// actual launch/join path.
    pub fn run_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.run_with_gate_witnessed(kernel, &CancelGate::new())
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

    /// [`Self::run_declared`] plus executor-minted completion evidence.
    pub fn run_declared_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.run_inner_witnessed(
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

    /// [`Self::run_declared_budgeted`] plus executor-minted completion
    /// evidence.
    pub fn run_declared_budgeted_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.run_inner_witnessed(
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

    /// Run a kernel under the exact cancellation authority and budget of an
    /// ambient executor [`Cx`], while retaining an explicit logical run
    /// identity and operation memory lease.
    ///
    /// This is the nested-throughput bridge for callers that already execute
    /// under an `fs-exec` context: it deliberately reuses `outer`'s gate
    /// instead of creating a second cancellation authority. The worker lane,
    /// memory admission, structured failures, full drain, and [`RunReport`]
    /// semantics are exactly those of
    /// [`Self::run_declared_leased_budgeted`].
    pub fn run_declared_leased_with_cx<K: TileKernel>(
        &self,
        outer: &Cx<'_>,
        kernel: &K,
        run: RunId,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> (Result<K::Out, RunError>, RunReport)
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_declared_leased_budgeted(kernel, outer.cancel_gate(), run, outer.budget(), lease)
    }

    /// [`Self::run_declared_leased_budgeted`] plus executor-minted completion
    /// evidence. The witness records the shared lease before and after
    /// executor-transient release; caller-owned retained output may
    /// legitimately keep the shared lease nonzero.
    pub fn run_declared_leased_budgeted_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_inner_witnessed(kernel, gate, run, budget, lease, Launch::<NoTask>::OwnScope)
    }

    /// Permit-consuming form of
    /// [`Self::run_declared_leased_budgeted_witnessed`].
    pub fn run_declared_leased_budgeted_witnessed_once<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_inner_witnessed_with_permit(
            kernel,
            gate,
            run,
            budget,
            lease,
            Some(permit.into_root()),
            Launch::<NoTask>::OwnScope,
        )
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

    /// [`Self::run_scoped`] plus executor-minted completion evidence.
    pub fn run_scoped_witnessed<Caps, K: TileKernel>(
        &self,
        task_cx: &asupersync::Cx<Caps>,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        Caps: Send + Sync + 'static,
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_inner_witnessed(kernel, gate, run, budget, lease, Launch::TaskScope(task_cx))
    }

    /// Permit-consuming form of [`Self::run_scoped_witnessed`].
    pub fn run_scoped_witnessed_once<Caps, K: TileKernel>(
        &self,
        task_cx: &asupersync::Cx<Caps>,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        Caps: Send + Sync + 'static,
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_inner_witnessed_with_permit(
            kernel,
            gate,
            run,
            budget,
            lease,
            Some(permit.into_root()),
            Launch::TaskScope(task_cx),
        )
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
        let dispatch_admission = AtomicBool::new(false);
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
                scope: CompletionScopeIdentity::task_parked(task_cx),
                dispatch_admission: &dispatch_admission,
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
        let dispatch_admission = AtomicBool::new(false);
        std::thread::scope(|s| {
            let _shutdown = crate::crew::CrewShutdown(&crew);
            for w in 0..crew.workers() {
                let crew = &crew;
                s.spawn(move || crew.park_loop(w, None));
            }
            f(&ParkedTilePool {
                pool: self,
                crew: &crew,
                scope: CompletionScopeIdentity::STD_PARKED,
                dispatch_admission: &dispatch_admission,
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

    /// [`Self::run_with_gate`] plus executor-minted completion evidence.
    pub fn run_with_gate_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.run_inner_witnessed(
            kernel,
            gate,
            RunId::default(),
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::<NoTask>::OwnScope,
        )
    }

    /// Consume one affine invocation permit and run under an external gate.
    ///
    /// Unlike the standalone compatibility method, the returned witness
    /// binds the permit root and therefore can participate in a one-shot
    /// invocation proof.
    pub fn run_with_gate_witnessed_once<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.run_inner_witnessed_with_permit(
            kernel,
            gate,
            RunId::default(),
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Some(permit.into_root()),
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
        let (result, report, witness) =
            self.run_inner_core::<false, Caps, K>(kernel, gate, run, budget, lease, None, launch);
        debug_assert!(witness.is_none());
        (result, report)
    }

    // This is deliberately one lifecycle: plan/admit -> launch -> join ->
    // consume/drop staging -> release root charge -> seal. Moving witness
    // construction into a caller wrapper would lose the authority boundary.
    fn run_inner_witnessed<Caps, K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
        launch: Launch<'_, Caps>,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        Caps: Send + Sync + 'static,
    {
        self.run_inner_witnessed_with_permit(kernel, gate, run, budget, lease, None, launch)
    }

    fn run_inner_witnessed_with_permit<Caps, K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
        affine_invocation_permit_root: Option<[u8; 32]>,
        launch: Launch<'_, Caps>,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        Caps: Send + Sync + 'static,
    {
        let (result, report, witness) = self.run_inner_core::<true, Caps, K>(
            kernel,
            gate,
            run,
            budget,
            lease,
            affine_invocation_permit_root,
            launch,
        );
        let witness =
            witness.ok_or_else(|| completion_bundle_invariant("completion-witness-missing"))??;
        let bundle = WitnessedRun {
            outcome: result,
            report,
            witness,
        };
        bundle.verify_bundle()?;
        Ok(bundle)
    }

    // `COMPLETION=false` is the legacy zero-evidence lane. Const propagation
    // removes witness counters, allocator snapshots, and hashing so the
    // additive evidence API does not tax existing hot-kernel callers.
    #[allow(clippy::too_many_lines)]
    fn run_inner_core<const COMPLETION: bool, Caps, K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
        affine_invocation_permit_root: Option<[u8; 32]>,
        launch: Launch<'_, Caps>,
    ) -> (
        Result<K::Out, RunError>,
        RunReport,
        Option<Result<TilePoolCompletionWitness, TilePoolCompletionWitnessError>>,
    )
    where
        Caps: Send + Sync + 'static,
    {
        let plan = kernel.tiles();
        let kernel_id = plan.kernel_id();
        let n = plan.tiles;
        let mode = self.config.mode.name();
        let completion_before = COMPLETION.then(|| CompletionBefore {
            scope: launch.completion_scope(),
            pool_placement_identity: self.placement_identity(),
            arena: CompletionArenaSnapshot::capture(&self.arenas),
            lease: CompletionLeaseSnapshot::capture(lease),
            cancellation_requested_at_entry: gate.is_requested(),
            planned_crew_callbacks: launch.planned_crew_callbacks(),
        });
        let prelaunch = |error: RunError,
                         planned_workers: u64,
                         root_metadata_bytes: u64,
                         root_charge_admitted: bool,
                         root_charge_released: bool| {
            let report = prelaunch_report(plan.kernel, mode, run, n);
            let cancellation_requested_at_terminal = gate.is_requested();
            let cancellation_requested = gate.is_requested();
            let witness = completion_before.as_ref().map(|before| {
                mint_completion_witness(CompletionWitnessFields {
                    pool_placement_identity: before.pool_placement_identity.clone(),
                    pool_seed: self.config.seed,
                    affine_invocation_permit_root,
                    kernel: plan.kernel,
                    kernel_id,
                    declared_run: run,
                    mode,
                    scope: before.scope,
                    planned_tiles: n,
                    admission_completed: false,
                    admitted_tiles: 0,
                    claimed_tiles: 0,
                    completed_tiles: 0,
                    break_tiles: 0,
                    panicked_tiles: 0,
                    planned_workers,
                    launched_workers: 0,
                    joined_workers: 0,
                    planned_crew_callbacks: before.planned_crew_callbacks,
                    entered_crew_callbacks: 0,
                    exited_crew_callbacks: 0,
                    tile_scopes_opened: 0,
                    live_tile_scopes_at_seal: 0,
                    cancellation_requested_at_entry: before.cancellation_requested_at_entry,
                    cancellation_requested_at_terminal,
                    cancellation_requested,
                    cancellation_observed_workers: 0,
                    root_metadata_bytes,
                    root_charge_admitted,
                    root_charge_released,
                    arena_before: before.arena,
                    arena_after: CompletionArenaSnapshot::capture(&self.arenas),
                    lease_before: before.lease,
                    lease_after: CompletionLeaseSnapshot::capture(lease),
                    terminal_error: Some(error.clone()),
                })
            });
            (Err(error), report, witness)
        };

        // Stream identity is DECLARED, never scheduled (wf9.7.1): the
        // former pool-global counter made keys depend on unrelated
        // prior runs and on concurrent invocation order.
        let iteration = run.0;
        let Ok(n_usize) = usize::try_from(n) else {
            return prelaunch(
                RunError::MemoryPlanOverflow {
                    kernel: plan.kernel,
                    what: "tile-count",
                },
                0,
                0,
                false,
                false,
            );
        };
        let workers = self.config.workers.min(n_usize.max(1)).max(1);
        let planned_workers = u64::try_from(workers).unwrap_or(u64::MAX);
        let _parked_dispatch_admission = match launch.dispatch_admission() {
            Some(admission) => match ParkedDispatchAdmission::try_acquire(admission) {
                Some(guard) => Some(guard),
                None => {
                    return prelaunch(
                        RunError::ParkedCrewBusy {
                            kernel: plan.kernel,
                        },
                        planned_workers,
                        0,
                        false,
                        false,
                    );
                }
            },
            None => None,
        };

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
                return prelaunch(
                    RunError::MemoryPlanOverflow {
                        kernel: plan.kernel,
                        what,
                    },
                    planned_workers,
                    0,
                    false,
                    false,
                );
            }
        };
        let root_charge = match lease.reserve("tilepool-root-metadata", root_bytes) {
            Ok(charge) => charge,
            Err(refusal) => {
                return prelaunch(
                    RunError::MemoryRefused {
                        kernel: plan.kernel,
                        what: refusal.what,
                        requested_bytes: refusal.requested_bytes,
                        used_bytes: refusal.used_bytes,
                        limit_bytes: refusal.limit_bytes,
                    },
                    planned_workers,
                    root_bytes,
                    false,
                    false,
                );
            }
        };

        let root = match allocate_run_root::<K>(
            n,
            n_usize,
            workers,
            &self.config.quantum_weights,
            self.config.topo,
            plan.kernel,
        ) {
            Ok(root) => root,
            Err(error) => {
                drop(root_charge);
                return prelaunch(error, planned_workers, root_bytes, true, true);
            }
        };

        let (
            result,
            report,
            terminal_error,
            claimed_tiles,
            break_tiles,
            panicked_tiles,
            launched_workers,
            joined_workers,
            entered_crew_callbacks,
            exited_crew_callbacks,
            tile_scopes_opened,
            live_tile_scopes_at_seal,
            cancellation_requested_at_terminal,
            cancellation_observed_workers,
        ) = {
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
            } = root;

            let claimed = AtomicU64::new(0);
            let breaks = AtomicU64::new(0);
            let panics = AtomicU64::new(0);
            let workers_entered = AtomicU64::new(0);
            let workers_exited = AtomicU64::new(0);
            let crew_callbacks_entered = AtomicU64::new(0);
            let crew_callbacks_exited = AtomicU64::new(0);
            let tile_scopes_opened_counter = AtomicU64::new(0);
            let tile_scopes_live = AtomicU64::new(0);
            let cancellation_observed_workers = AtomicU64::new(0);
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
                claimed: &claimed,
                breaks: &breaks,
                panics: &panics,
                workers_entered: &workers_entered,
                workers_exited: &workers_exited,
                tile_scopes_opened: &tile_scopes_opened_counter,
                tile_scopes_live: &tile_scopes_live,
                cancellation_observed_workers: &cancellation_observed_workers,
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
                            #[cfg(test)]
                            if take_worker_spawn_refusal(w) {
                                spawn_failure =
                                    Some((w, "test-injected worker spawn refusal".to_string()));
                                gate.request();
                                break;
                            }
                            let spawned = std::thread::Builder::new().spawn_scoped(s, move || {
                                worker_loop::<COMPLETION, Caps, K>(ctx, w, None)
                            });
                            if let Err(error) = spawned {
                                spawn_failure = Some((w, error.to_string()));
                                gate.request();
                                break;
                            }
                        }
                    });
                }
                Launch::Crew { crew, .. } => {
                    // Every parked worker enters the callback, even when the
                    // active tile plan uses a smaller worker subset. Dispatch
                    // returns only after this complete callback set exits.
                    let job = |w: usize, cpu: Option<&CpuCx<Caps>>| {
                        if COMPLETION {
                            crew_callbacks_entered.fetch_add(1, Ordering::Release);
                        }
                        let _callback_completion = CrewCallbackCompletionGuard {
                            exited: COMPLETION.then_some(&crew_callbacks_exited),
                        };
                        if w < ctx.workers {
                            worker_loop::<COMPLETION, Caps, K>(&ctx, w, cpu);
                        }
                    };
                    if let Some((worker, message)) = crew.dispatch(&job) {
                        std::panic::panic_any(format!(
                            "tile-pool worker {worker} panicked outside tile containment: {message}"
                        ));
                    }
                }
                Launch::TaskScope(cx) => {
                    match cx.scoped_cpu(workers, |scope| {
                        for w in 0..workers {
                            let ctx = &ctx;
                            if let Err(error) = scope.spawn(move |cpu| {
                                worker_loop::<COMPLETION, Caps, K>(ctx, w, Some(cpu))
                            }) {
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
                    ScopedCpuError::Cancelled(_) => gate.request(),
                    ScopedCpuError::ChildPanicked { child, message } => {
                        std::panic::panic_any(format!(
                            "tile-pool worker {child} panicked outside tile containment: {message}"
                        ));
                    }
                    ScopedCpuError::WorkerCapExceeded { cap } => {
                        std::panic::panic_any(format!(
                            "tile-pool scoped launch exceeded its own worker cap {cap}"
                        ));
                    }
                }
            }

            // Every launch harness above is a blocking join boundary. These
            // counters come from the actual worker-loop entry/exit guards,
            // rather than from a caller-authored tracker.
            let launched_workers = if COMPLETION {
                workers_entered.load(Ordering::Acquire)
            } else {
                0
            };
            let joined_workers = if COMPLETION {
                workers_exited.load(Ordering::Acquire)
            } else {
                0
            };
            let entered_crew_callbacks = if COMPLETION {
                crew_callbacks_entered.load(Ordering::Acquire)
            } else {
                0
            };
            let exited_crew_callbacks = if COMPLETION {
                crew_callbacks_exited.load(Ordering::Acquire)
            } else {
                0
            };
            let claimed_tiles = if COMPLETION {
                claimed.load(Ordering::Acquire)
            } else {
                0
            };
            let break_tiles = if COMPLETION {
                breaks.load(Ordering::Acquire)
            } else {
                0
            };
            let panicked_tiles = if COMPLETION {
                panics.load(Ordering::Acquire)
            } else {
                0
            };
            let tile_scopes_opened = if COMPLETION {
                tile_scopes_opened_counter.load(Ordering::Acquire)
            } else {
                0
            };
            let live_tile_scopes_at_seal = if COMPLETION {
                tile_scopes_live.load(Ordering::Acquire)
            } else {
                0
            };
            let cancellation_observed_workers = if COMPLETION {
                cancellation_observed_workers.load(Ordering::Acquire)
            } else {
                0
            };

            let completed = slots
                .iter()
                .filter(|slot| slot.lock().expect("slot").is_some())
                .count() as u64;
            if let Some(requested_at) = gate.requested_at_ns() {
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
                mode,
                declared_run: run,
                completed,
                total: n,
                steals: steals.load(Ordering::Relaxed),
                cross_ccd_steals: cross_steals.load(Ordering::Relaxed),
                cancel_latencies_ns,
                tiles_by_worker,
            };

            let panic_failure = panic_box.into_inner().expect("panic box");
            let typed_failure = refusal_sink.take();
            let cancellation_requested_at_terminal = gate.is_requested();
            // Stable failure-class precedence preserves the existing public
            // RunError semantics; the exact selected value is cloned into the
            // immutable witness and rechecked by its verifier.
            let result = if let Some((worker, message)) = spawn_failure {
                Err(RunError::WorkerSpawn {
                    kernel: plan.kernel,
                    worker,
                    message,
                })
            } else if let Some((tile, message)) = panic_failure {
                Err(RunError::TilePanicked {
                    kernel: plan.kernel,
                    tile,
                    message,
                    completed,
                })
            } else if let Some((tile, failure)) = typed_failure {
                Err(RunError::TileFailed {
                    kernel: plan.kernel,
                    tile,
                    failure,
                    completed,
                })
            } else if cancellation_requested_at_terminal {
                Err(RunError::Cancelled {
                    kernel: plan.kernel,
                    completed,
                    total: n,
                })
            } else {
                let mut missing = None;
                for (index, slot) in slots.into_iter().enumerate() {
                    match slot.into_inner().expect("slot") {
                        Some(out) => outs.push(out),
                        None => {
                            missing = Some(index as u64);
                            break;
                        }
                    }
                }
                match missing {
                    Some(tile) => Err(RunError::Incomplete {
                        kernel: plan.kernel,
                        tile,
                    }),
                    None => match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        crate::reduce::pairwise_fold(outs)
                    })) {
                        Ok(out) => Ok(out),
                        Err(payload) => {
                            let message = payload
                                .downcast_ref::<&str>()
                                .map(ToString::to_string)
                                .or_else(|| payload.downcast_ref::<String>().cloned())
                                .unwrap_or_else(|| "non-string panic payload".to_string());
                            Err(RunError::ReductionPanicked {
                                kernel: plan.kernel,
                                message,
                            })
                        }
                    },
                }
            };
            let terminal_error = result.as_ref().err().cloned();
            (
                result,
                report,
                terminal_error,
                claimed_tiles,
                break_tiles,
                panicked_tiles,
                launched_workers,
                joined_workers,
                entered_crew_callbacks,
                exited_crew_callbacks,
                tile_scopes_opened,
                live_tile_scopes_at_seal,
                cancellation_requested_at_terminal,
                cancellation_observed_workers,
            )
        };

        // The output/result and report have left staging. All other root
        // vectors dropped at the block boundary above; now release the
        // executor's own root lease charge before observing/sealing evidence.
        drop(root_charge);
        let cancellation_requested = gate.is_requested();
        let witness = completion_before.map(|before| {
            mint_completion_witness(CompletionWitnessFields {
                pool_placement_identity: before.pool_placement_identity,
                pool_seed: self.config.seed,
                affine_invocation_permit_root,
                kernel: plan.kernel,
                kernel_id,
                declared_run: run,
                mode,
                scope: before.scope,
                planned_tiles: n,
                admission_completed: true,
                admitted_tiles: n,
                claimed_tiles,
                completed_tiles: report.completed,
                break_tiles,
                panicked_tiles,
                planned_workers,
                launched_workers,
                joined_workers,
                planned_crew_callbacks: before.planned_crew_callbacks,
                entered_crew_callbacks,
                exited_crew_callbacks,
                tile_scopes_opened,
                live_tile_scopes_at_seal,
                cancellation_requested_at_entry: before.cancellation_requested_at_entry,
                cancellation_requested_at_terminal,
                cancellation_requested,
                cancellation_observed_workers,
                root_metadata_bytes: root_bytes,
                root_charge_admitted: true,
                root_charge_released: true,
                arena_before: before.arena,
                arena_after: CompletionArenaSnapshot::capture(&self.arenas),
                lease_before: before.lease,
                lease_after: CompletionLeaseSnapshot::capture(lease),
                terminal_error,
            })
        });
        (result, report, witness)
    }
}

impl CompletionKernelRunner for TilePool {
    fn workers(&self) -> usize {
        TilePool::workers(self)
    }

    fn run_with_gate_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        TilePool::run_with_gate_witnessed(self, kernel, gate)
    }

    fn run_with_gate_witnessed_once<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        TilePool::run_with_gate_witnessed_once(self, kernel, gate, permit)
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
    scope: CompletionScopeIdentity,
    dispatch_admission: &'a AtomicBool,
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

    /// [`TilePool::run_witnessed`] on the parked crew.
    pub fn run_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.run_with_gate_witnessed(kernel, &CancelGate::new())
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
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
        )
    }

    /// [`TilePool::run_with_gate_witnessed`] on the parked crew.
    pub fn run_with_gate_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.pool.run_inner_witnessed(
            kernel,
            gate,
            RunId::default(),
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
        )
    }

    /// Permit-consuming [`Self::run_with_gate_witnessed`] on the parked crew.
    pub fn run_with_gate_witnessed_once<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.pool.run_inner_witnessed_with_permit(
            kernel,
            gate,
            RunId::default(),
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Some(permit.into_root()),
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
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
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
        )
    }

    /// [`TilePool::run_declared_witnessed`] on the parked crew.
    pub fn run_declared_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        self.pool.run_inner_witnessed(
            kernel,
            gate,
            run,
            Budget::INFINITE,
            &fs_alloc::OperationMemoryLease::unbounded(),
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
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
        self.pool.run_inner(
            kernel,
            gate,
            run,
            budget,
            lease,
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
        )
    }

    /// Run on the parked crew under the exact cancellation authority and
    /// budget of an ambient executor [`Cx`]. This is the parked counterpart
    /// of [`TilePool::run_declared_leased_with_cx`].
    pub fn run_declared_leased_with_cx<K: TileKernel>(
        &self,
        outer: &Cx<'_>,
        kernel: &K,
        run: RunId,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> (Result<K::Out, RunError>, RunReport)
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        self.run_declared_leased_budgeted(kernel, outer.cancel_gate(), run, outer.budget(), lease)
    }

    /// [`TilePool::run_declared_leased_budgeted_witnessed`] on the parked
    /// crew.
    pub fn run_declared_leased_budgeted_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        self.pool.run_inner_witnessed(
            kernel,
            gate,
            run,
            budget,
            lease,
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
        )
    }

    /// Permit-consuming
    /// [`Self::run_declared_leased_budgeted_witnessed`] on the parked crew.
    pub fn run_declared_leased_budgeted_witnessed_once<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        run: RunId,
        budget: Budget,
        lease: &fs_alloc::OperationMemoryLease,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError>
    where
        K::Out: crate::LeaseAdmittedOut,
    {
        self.pool.run_inner_witnessed_with_permit(
            kernel,
            gate,
            run,
            budget,
            lease,
            Some(permit.into_root()),
            Launch::Crew {
                crew: self.crew,
                scope: self.scope,
                dispatch_admission: self.dispatch_admission,
            },
        )
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

impl<Caps: Send + Sync + 'static> CompletionKernelRunner for ParkedTilePool<'_, Caps> {
    fn workers(&self) -> usize {
        ParkedTilePool::workers(self)
    }

    fn run_with_gate_witnessed<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        ParkedTilePool::run_with_gate_witnessed(self, kernel, gate)
    }

    fn run_with_gate_witnessed_once<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
        permit: TilePoolInvocationPermit,
    ) -> Result<WitnessedRun<K::Out>, TilePoolCompletionWitnessError> {
        ParkedTilePool::run_with_gate_witnessed_once(self, kernel, gate, permit)
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

    macro_rules! witnessed_parts {
        ($run:expr) => {{
            let bundle = $run.expect("executor must seal a valid witnessed-run bundle");
            bundle
                .verify_bundle()
                .expect("fresh executor bundle must verify");
            bundle.into_parts()
        }};
    }

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

    struct BlockingParkedKernel {
        entered: std::sync::Arc<std::sync::Barrier>,
        release: std::sync::Arc<std::sync::Barrier>,
    }

    impl TileKernel for BlockingParkedKernel {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/blocking-parked-dispatch", 1)
        }

        fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            self.entered.wait();
            self.release.wait();
            if cx.checkpoint().is_err() {
                ControlFlow::Break(crate::Cancelled)
            } else {
                ControlFlow::Continue(1)
            }
        }
    }

    struct RecursiveParkedKernel<'pool, 'crew> {
        parked: &'pool ParkedTilePool<'crew, NoTask>,
        nested: Mutex<Option<WitnessedRun<u64>>>,
    }

    impl TileKernel for RecursiveParkedKernel<'_, '_> {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/recursive-parked-dispatch", 1)
        }

        fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(crate::Cancelled);
            }
            let nested = self
                .parked
                .run_witnessed(&NoAllocation)
                .expect("recursive dispatch refusal must itself seal a valid witness");
            self.nested
                .lock()
                .expect("recursive result")
                .replace(nested);
            ControlFlow::Continue(1)
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

    fn with_outer_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        const OUTER_STREAM: StreamKey = StreamKey {
            seed: 0x4F55_5445_525F_4358,
            kernel_id: 0x4252_4944_4745,
            tile: 3,
            iteration: 5,
        };

        let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let result = arenas.scope(|arena| {
            let outer = Cx::new(gate, arena, OUTER_STREAM, budget, ExecMode::Deterministic);
            f(&outer)
        });
        assert!(
            arenas.stats().quiescent(),
            "ambient Cx arena must be quiescent after the nested run"
        );
        result
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

    const CANCEL_HANDSHAKE_WORKERS: usize = 4;

    /// A kernel that cancels the CALLING asupersync task from inside a tile
    /// — the G4 mid-run cancellation storm's trigger. The two barrier
    /// generations make the trigger schedule-independent: every worker has
    /// entered exactly one tile before `at` requests cancellation, and no
    /// first-wave tile can return until that request is visible. The workers'
    /// next REAL tile boundary must therefore bridge the task request into the
    /// pool gate and drain.
    struct CancelTaskAt {
        tiles: u64,
        at: u64,
        task: asupersync::Cx,
        rendezvous: std::sync::Barrier,
    }

    impl CancelTaskAt {
        fn synchronized(tiles: u64, at: u64, task: asupersync::Cx) -> Self {
            Self {
                tiles,
                at,
                task,
                rendezvous: std::sync::Barrier::new(CANCEL_HANDSHAKE_WORKERS),
            }
        }
    }

    impl TileKernel for CancelTaskAt {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("test/cancel-task-at", self.tiles)
        }

        fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
            self.rendezvous.wait();
            if tile == self.at {
                self.task.set_cancel_requested(true);
            }
            self.rendezvous.wait();
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
            let p = pool(CANCEL_HANDSHAKE_WORKERS);
            let kernel = CancelTaskAt::synchronized(16_384, 0, cx.clone());
            let (out, report) = run_scoped_simple(&p, &cx, &kernel);
            match out {
                Err(RunError::Cancelled {
                    completed, total, ..
                }) => {
                    assert_eq!(total, 16_384);
                    assert_eq!(
                        completed, CANCEL_HANDSHAKE_WORKERS as u64,
                        "only the synchronized first wave may complete before the next \
                         boundary drains (completed {completed})"
                    );
                }
                other => panic!("task cancel must fail closed as Cancelled: {other:?}"),
            }
            assert_eq!(
                report.completed, CANCEL_HANDSHAKE_WORKERS as u64,
                "report agrees that the first post-request boundaries drained"
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
    fn declared_leased_with_cx_inherits_budget_and_run_identity() {
        let pool = pool(4);
        let gate = CancelGate::new();
        let budget = Budget::new().with_cost_quota(65_536);
        let lease = fs_alloc::OperationMemoryLease::unbounded();

        let (result, report) = with_outer_cx(&gate, budget, |outer| {
            pool.run_declared_leased_with_cx(
                outer,
                &BudgetProbe { tiles: 4 },
                RunId(0x4358),
                &lease,
            )
        });

        assert_eq!(
            result.expect("ambient-Cx run"),
            4 * 65_536,
            "every tile must receive the outer context's exact budget"
        );
        assert_eq!(report.declared_run, RunId(0x4358));
        assert_eq!((report.completed, report.total), (4, 4));
        assert!(!gate.is_requested());
        assert!(pool.arena_pool().stats().quiescent());
        let receipt = lease.receipt();
        assert!(receipt.requested_bytes > 0, "root metadata was admitted");
        assert_eq!(receipt.used_bytes, 0, "all transient charges were released");
    }

    #[test]
    fn declared_leased_with_cx_honors_ambient_pre_cancel() {
        struct MustNotRun;

        impl TileKernel for MustNotRun {
            type Out = ();

            fn tiles(&self) -> TilePlan {
                TilePlan::new("test/ambient-pre-cancel", 8)
            }

            fn run(&self, _tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, ()> {
                panic!("a pre-cancelled ambient context must not execute a tile")
            }
        }

        let pool = pool(2);
        let gate = CancelGate::new_clock_free();
        gate.request();
        let lease = fs_alloc::OperationMemoryLease::unbounded();

        let (result, report) = with_outer_cx(&gate, Budget::INFINITE, |outer| {
            pool.run_declared_leased_with_cx(outer, &MustNotRun, RunId(0xCA11), &lease)
        });

        assert!(matches!(
            result,
            Err(RunError::Cancelled {
                completed: 0,
                total: 8,
                ..
            })
        ));
        assert_eq!(report.declared_run, RunId(0xCA11));
        assert_eq!((report.completed, report.total), (0, 8));
        assert!(report.cancel_latencies_ns.is_empty());
        assert!(pool.arena_pool().stats().quiescent());
        assert_eq!(lease.receipt().used_bytes, 0);
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
        let (result, report, witness) =
            witnessed_parts!(pool.run_declared_leased_budgeted_witnessed(
                &UnrepresentablePlan,
                &CancelGate::new(),
                RunId(29),
                Budget::INFINITE,
                &lease,
            ));
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
        assert_eq!(witness.verify(), Ok(()));
        assert_eq!(
            witness.disposition(),
            TilePoolCompletionDisposition::MemoryPlanOverflow
        );
        assert!(!witness.admission_completed());
        assert_eq!(witness.admitted_tiles(), 0);
        assert_eq!(witness.claimed_tiles(), 0);
        assert_eq!(witness.launched_workers(), 0);
        assert_eq!(witness.joined_workers(), 0);
        assert!(!witness.root_charge_admitted());
        assert!(!witness.root_charge_released());
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
            let (result, report, witness) = witnessed_parts!(pool.run_declared_budgeted_witnessed(
                &kernel,
                &gate,
                RunId(23),
                Budget::new().with_cost_quota(1 << 20),
            ));
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
            assert_eq!(witness.verify(), Ok(()));
            assert_eq!(
                witness.disposition(),
                TilePoolCompletionDisposition::TileFailed
            );
            assert_eq!(witness.first_failure_kind(), Some("tile-failed"));
            assert_eq!(witness.first_failure_tile(), Some(0));
            assert_eq!(witness.claimed_tiles(), workers as u64);
            assert_eq!(witness.break_tiles(), workers as u64);
            assert_eq!(witness.failed_tiles(), 1);
            assert_eq!(witness.cancelled_tiles(), workers.saturating_sub(1) as u64);
            assert_eq!(witness.launched_workers(), workers as u64);
            assert_eq!(witness.joined_workers(), workers as u64);
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
                let p = pool(workers);
                let (result, _, witness) = witnessed_parts!(p.run_witnessed(&kernel));
                let error = result.expect_err("every in-flight tile panics");
                match error {
                    RunError::TilePanicked { tile, message, .. } => {
                        assert_eq!(tile, 0, "panic provenance must not depend on arrival order");
                        assert_eq!(message, "simultaneous panic from tile 0");
                    }
                    other => panic!("expected TilePanicked, got {other:?}"),
                }
                assert_eq!(witness.verify(), Ok(()));
                assert_eq!(
                    witness.disposition(),
                    TilePoolCompletionDisposition::TilePanicked
                );
                assert_eq!(witness.first_failure_kind(), Some("tile-panicked"));
                assert_eq!(witness.first_failure_tile(), Some(0));
                assert_eq!(witness.claimed_tiles(), workers as u64);
                assert_eq!(witness.panicked_tiles(), workers as u64);
                assert_eq!(witness.failed_tiles(), workers as u64);
                assert_eq!(witness.cancelled_tiles(), 0);
                assert_eq!(witness.launched_workers(), workers as u64);
                assert_eq!(witness.joined_workers(), workers as u64);
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
    fn clock_free_gate_never_mints_cancel_latency_samples() {
        let p = pool(2);
        let gate = CancelGate::new_clock_free();
        gate.request();

        let (result, report) = p.run_with_gate(&SumKernel { tiles: 64 }, &gate);
        assert!(matches!(result, Err(RunError::Cancelled { .. })));
        assert!(report.cancel_latencies_ns.is_empty());
        assert_eq!(report.cancel_latency_p99_ns(), None);
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
        struct RequestGateAfterFirstWave<'a> {
            tiles: u64,
            gate: &'a CancelGate,
            rendezvous: std::sync::Barrier,
        }
        impl TileKernel for RequestGateAfterFirstWave<'_> {
            type Out = u64;

            fn tiles(&self) -> TilePlan {
                TilePlan::new("test/request-gate-after-first-wave", self.tiles)
            }

            fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
                self.rendezvous.wait();
                if tile == 0 {
                    self.gate.request();
                }
                self.rendezvous.wait();
                ControlFlow::Continue(1)
            }
        }
        let p = pool(CANCEL_HANDSHAKE_WORKERS);
        p.with_parked_crew_local(|parked| {
            let gate = CancelGate::new();
            let kernel = RequestGateAfterFirstWave {
                tiles: 16_384,
                gate: &gate,
                rendezvous: std::sync::Barrier::new(CANCEL_HANDSHAKE_WORKERS),
            };
            let (out, report) = parked.run_with_gate(&kernel, &gate);
            match out {
                Err(RunError::Cancelled {
                    completed, total, ..
                }) => {
                    assert_eq!(total, 16_384);
                    assert_eq!(
                        completed, CANCEL_HANDSHAKE_WORKERS as u64,
                        "the gate drains the parked run after its first wave"
                    );
                }
                other => panic!("expected Cancelled, got {other:?}"),
            }
            assert_eq!(report.completed, CANCEL_HANDSHAKE_WORKERS as u64);
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
                    let kernel = CancelTaskAt::synchronized(16_384, 0, cx.clone());
                    let (out, _) = parked.run_with_gate(&kernel, &CancelGate::new());
                    match out {
                        Err(RunError::Cancelled {
                            completed, total, ..
                        }) => {
                            assert_eq!(total, 16_384);
                            assert_eq!(
                                completed, CANCEL_HANDSHAKE_WORKERS as u64,
                                "task cancel drains the parked run after its first wave"
                            );
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

    fn assert_completion_invariant(witness: &TilePoolCompletionWitness, expected: &'static str) {
        assert_eq!(
            witness.verify(),
            Err(TilePoolCompletionWitnessError::Invariant { name: expected }),
            "mutated witness should fail `{expected}`: {}",
            witness.to_canonical_json()
        );
    }

    #[test]
    fn completion_witness_binds_the_real_asupersync_parent_scope() {
        in_task(|cx| {
            let expected_region = cx.region_id().as_u64();
            let expected_task = cx.task_id().as_u64();
            let p = pool(2);
            let (result, _, witness) = witnessed_parts!(p.run_scoped_witnessed(
                &cx,
                &SumKernel { tiles: 8 },
                &CancelGate::new_clock_free(),
                RunId(301),
                Budget::INFINITE,
                &fs_alloc::OperationMemoryLease::unbounded(),
            ));
            assert_eq!(result, Ok((1..=8).sum::<u64>()));
            assert_eq!(witness.verify(), Ok(()));
            assert_eq!(witness.scope_kind(), "asupersync-task-scope");
            assert_eq!(witness.parent_region_id(), Some(expected_region));
            assert_eq!(witness.parent_task_id(), Some(expected_task));
            println!(
                "{}",
                witness.to_jsonl("asupersync-parent-scope", 0, "inline-real-task")
            );
        });
    }

    #[test]
    fn affine_invocation_permit_is_consumed_and_binds_call_identity() {
        const CROSS_CALL_NO_CLAIM: &str = "cross-call-uniqueness-without-affine-invocation-permit";

        let p = pool(2);
        let gate = CancelGate::new_clock_free();
        let ordinal_zero_permit_root = [0xA5; 32];
        let permit = TilePoolInvocationPermit::from_permit_root(ordinal_zero_permit_root);
        assert_eq!(permit.permit_root(), ordinal_zero_permit_root);

        let bound = p
            .run_with_gate_witnessed_once(&NoAllocation, &gate, permit)
            .expect("permit-bound run");
        bound.verify_bundle().expect("permit-bound bundle");
        assert_eq!(bound.outcome(), &Ok(1));
        assert!(bound.witness().has_affine_invocation_permit());
        assert_eq!(
            bound.witness().affine_invocation_permit_root(),
            Some(ordinal_zero_permit_root)
        );
        assert!(
            !bound.witness().no_claims().contains(&CROSS_CALL_NO_CLAIM),
            "a consumed affine permit discharges the standalone uniqueness no-claim"
        );

        let standalone = p
            .run_with_gate_witnessed(&NoAllocation, &CancelGate::new_clock_free())
            .expect("standalone witnessed run");
        assert!(!standalone.witness().has_affine_invocation_permit());
        assert!(
            standalone
                .witness()
                .no_claims()
                .contains(&CROSS_CALL_NO_CLAIM)
        );
        assert_eq!(
            bound.witness().plan_root_bytes(),
            standalone.witness().plan_root_bytes(),
            "the declared plan is unchanged by invocation authority"
        );
        assert_ne!(
            bound.witness().call_replay_root_bytes(),
            standalone.witness().call_replay_root_bytes(),
            "the permit root must distinguish the authority-bearing call"
        );

        let ordinal_one_permit_root = [0x5A; 32];
        let second = p
            .run_with_gate_witnessed_once(
                &NoAllocation,
                &CancelGate::new_clock_free(),
                TilePoolInvocationPermit::from_permit_root(ordinal_one_permit_root),
            )
            .expect("second independent permit-bound run");
        assert_eq!(
            second.witness().affine_invocation_permit_root(),
            Some(ordinal_one_permit_root)
        );
        assert_eq!(
            bound.witness().plan_root_bytes(),
            second.witness().plan_root_bytes(),
            "run ordinal authority changes the permit root, not the declared executor plan"
        );
        assert_ne!(
            bound.witness().call_replay_root_bytes(),
            second.witness().call_replay_root_bytes()
        );
        assert_ne!(bound.witness().root_bytes(), second.witness().root_bytes());

        let mut permit_corruption = bound.witness().clone();
        permit_corruption.affine_invocation_permit_root = Some(ordinal_one_permit_root);
        permit_corruption.root = completion_witness_root(&permit_corruption);
        assert_completion_invariant(&permit_corruption, "call-replay-root");
    }

    #[test]
    fn witnessed_run_bundle_rejects_cross_run_report_and_outcome_tampering() {
        let p = pool(2);

        let mut report_kernel = p
            .run_witnessed(&NoAllocation)
            .expect("fresh witnessed bundle");
        report_kernel.report.kernel = "test/not-the-executed-kernel";
        assert_eq!(
            report_kernel.verify_bundle(),
            Err(TilePoolCompletionWitnessError::BundleInvariant {
                name: "report-kernel"
            })
        );

        let mut report_completed = p
            .run_witnessed(&NoAllocation)
            .expect("fresh witnessed bundle");
        report_completed.report.completed = 0;
        assert_eq!(
            report_completed.verify_bundle(),
            Err(TilePoolCompletionWitnessError::BundleInvariant {
                name: "report-completed"
            })
        );

        let mut worker_conservation = p
            .run_witnessed(&NoAllocation)
            .expect("fresh witnessed bundle");
        worker_conservation.report.tiles_by_worker[0] =
            worker_conservation.report.tiles_by_worker[0].saturating_add(1);
        assert_eq!(
            worker_conservation.verify_bundle(),
            Err(TilePoolCompletionWitnessError::BundleInvariant {
                name: "report-worker-conservation"
            })
        );

        let mut terminal_outcome = p
            .run_witnessed(&NoAllocation)
            .expect("fresh witnessed bundle");
        terminal_outcome.outcome = Err(RunError::Incomplete {
            kernel: "test/no-allocation",
            tile: 0,
        });
        assert_eq!(
            terminal_outcome.verify_bundle(),
            Err(TilePoolCompletionWitnessError::BundleInvariant {
                name: "terminal-outcome"
            })
        );
    }

    /// G5: both retained roots are checked before lifecycle semantics, then a
    /// self-consistent rehash still cannot turn impossible join/quiescence
    /// states into executor-minted evidence.
    #[test]
    fn completion_witness_verifier_rejects_mutated_lifecycle_states() {
        let p = pool(2);
        let (result, _, base) = witnessed_parts!(p.run_witnessed(&NoAllocation));
        assert_eq!(result, Ok(1));
        assert_eq!(base.verify(), Ok(()));

        let mut plan_corruption = base.clone();
        plan_corruption.plan_root[0] ^= 0x80;
        assert_eq!(
            plan_corruption.verify(),
            Err(TilePoolCompletionWitnessError::PlanRootMismatch)
        );

        let mut root_corruption = base.clone();
        root_corruption.joined_workers = root_corruption.joined_workers.saturating_sub(1);
        assert_eq!(
            root_corruption.verify(),
            Err(TilePoolCompletionWitnessError::RootMismatch)
        );

        let mutations: &[(fn(&mut TilePoolCompletionWitness), &'static str)] = &[
            (
                |witness| {
                    witness.joined_workers = witness.joined_workers.saturating_sub(1);
                },
                "all-launched-workers-joined",
            ),
            (
                |witness| {
                    witness.launched_workers = witness.launched_workers.saturating_sub(1);
                    witness.joined_workers = witness.joined_workers.saturating_sub(1);
                },
                "claimed-work-needs-worker",
            ),
            (
                |witness| witness.worker_admission_closed = false,
                "worker-admission-closed",
            ),
            (
                |witness| witness.live_worker_guards_at_seal = 1,
                "no-live-worker-guards",
            ),
            (
                |witness| witness.live_tile_scopes_at_seal = 1,
                "no-live-tile-scopes",
            ),
            (
                |witness| witness.root_charge_released = false,
                "root-charge-release",
            ),
            (
                |witness| {
                    witness.completed_tiles = witness.completed_tiles.saturating_sub(1);
                },
                "claimed-terminal-conservation",
            ),
            (
                |witness| witness.disposition = TilePoolCompletionDisposition::Cancelled,
                "derived-disposition",
            ),
            (
                |witness| witness.scope.parent_region_id = Some(1),
                "scope-identity",
            ),
            (
                |witness| {
                    witness.cancellation_observed_workers =
                        witness.launched_workers.saturating_add(1);
                },
                "request-observation-within-launch",
            ),
            (
                |witness| witness.cancellation_requested_at_entry = true,
                "request-entry-monotonic",
            ),
            (
                |witness| witness.cancellation_requested_at_terminal = true,
                "request-terminal-monotonic",
            ),
            (
                |witness| witness.request_phase = TilePoolRequestPhase::BeforeEntry,
                "derived-request-phase",
            ),
            (
                |witness| {
                    witness.arena_before.quiescent = !witness.arena_before.quiescent;
                },
                "arena-before-internal",
            ),
            (
                |witness| {
                    witness.arena_after.free_bytes =
                        witness.arena_after.reserved_bytes.saturating_add(1);
                },
                "arena-after-internal",
            ),
            (
                |witness| {
                    witness.lease_after.limit_bytes = Some(
                        witness
                            .lease_after
                            .peak_bytes
                            .max(witness.lease_after.used_bytes),
                    );
                },
                "lease-limit-stable",
            ),
            (
                |witness| {
                    witness.lease_before.requested_bytes =
                        witness.lease_after.requested_bytes.saturating_add(1);
                },
                "lease-requested-monotonic",
            ),
            (
                |witness| {
                    witness.lease_before.peak_bytes =
                        witness.lease_after.peak_bytes.saturating_add(1);
                },
                "lease-peak-monotonic",
            ),
            (
                |witness| {
                    witness.lease_before.refusals = witness.lease_after.refusals.saturating_add(1);
                },
                "lease-refusals-monotonic",
            ),
            (
                |witness| {
                    witness.lease_after.release_invariant_violations = witness
                        .lease_before
                        .release_invariant_violations
                        .saturating_add(1);
                },
                "no-run-observed-lease-release-violation",
            ),
            (
                |witness| {
                    witness.lease_after.requested_bytes = witness.lease_before.requested_bytes;
                },
                "root-charge-request-observed",
            ),
        ];
        for (mutate, expected) in mutations {
            let mut witness = base.clone();
            mutate(&mut witness);
            witness.root = completion_witness_root(&witness);
            assert_completion_invariant(&witness, expected);
        }

        let cancelled_gate = CancelGate::new_clock_free();
        cancelled_gate.request();
        let (_, _, mut cancelled) =
            witnessed_parts!(p.run_with_gate_witnessed(&NoAllocation, &cancelled_gate));
        assert!(cancelled.cancellation_observed_workers() > 0);
        cancelled.cancellation_requested = false;
        cancelled.root = completion_witness_root(&cancelled);
        assert_completion_invariant(&cancelled, "request-terminal-monotonic");

        let mut config = PoolConfig::new(1, CcdTopology::APPLE_M_CLASS, 0xFA13);
        config.arena.limit_bytes = Some(0);
        let fault_pool = TilePool::new(config);
        let (_, _, mut first_failure) =
            witnessed_parts!(fault_pool.run_declared_budgeted_witnessed(
                &SimultaneousAllocationRefusal {
                    tiles: 1,
                    barrier: std::sync::Barrier::new(1),
                },
                &CancelGate::new_clock_free(),
                RunId(302),
                Budget::INFINITE,
            ));
        match first_failure.terminal_error.as_mut() {
            Some(RunError::TileFailed { tile, .. }) => *tile = 7,
            other => panic!("expected retained typed refusal, got {other:?}"),
        }
        first_failure.root = completion_witness_root(&first_failure);
        assert_completion_invariant(&first_failure, "retained-first-failure");
    }

    #[test]
    fn completion_request_phase_model_and_wide_lease_counters_are_lossless() {
        let p = pool(1);
        let (_, _, base) = witnessed_parts!(p.run_witnessed(&NoAllocation));
        assert_eq!(base.request_phase(), TilePoolRequestPhase::NotRequested);

        let mut after_terminal = base.clone();
        after_terminal.cancellation_requested = true;
        after_terminal.request_phase = TilePoolRequestPhase::AfterTerminalDecision;
        after_terminal.root = completion_witness_root(&after_terminal);
        assert_eq!(
            after_terminal.verify(),
            Ok(()),
            "a request first observed after terminal selection is distinct from run cancellation"
        );

        let mut wide = base;
        wide.lease_before.requested_bytes = u128::from(u64::MAX) + 41;
        wide.lease_after.requested_bytes = wide
            .lease_before
            .requested_bytes
            .saturating_add(u128::from(wide.root_metadata_bytes));
        wide.root = completion_witness_root(&wide);
        assert_eq!(wide.verify(), Ok(()));
        let json = wide.to_canonical_json();
        assert!(
            json.contains(&wide.lease_before.requested_bytes.to_string()),
            "canonical evidence must retain the full u128 counter"
        );
        assert!(
            json.contains(&wide.lease_after.requested_bytes.to_string()),
            "canonical evidence must retain the full u128 delta"
        );
    }

    /// Bounded lifecycle model: across every join/admission/live-guard state
    /// reachable in this two-worker fixture, the verifier accepts exactly the
    /// closed state where every entered worker has exited and no guard remains.
    #[test]
    fn completion_witness_bounded_join_model_accepts_only_closed_drained_states() {
        let p = pool(2);
        let (_, _, base) = witnessed_parts!(p.run_witnessed(&SumKernel { tiles: 2 }));
        assert_eq!(base.launched_workers, 2);

        let mut cases = 0_u64;
        for joined in 0..=base.launched_workers {
            for live in 0..=base.launched_workers {
                for admission_closed in [false, true] {
                    cases += 1;
                    let mut witness = base.clone();
                    witness.joined_workers = joined;
                    witness.live_worker_guards_at_seal = live;
                    witness.worker_admission_closed = admission_closed;
                    witness.root = completion_witness_root(&witness);

                    let should_verify =
                        joined == base.launched_workers && live == 0 && admission_closed;
                    assert_eq!(
                        witness.verify().is_ok(),
                        should_verify,
                        "bounded model case joined={joined}, live={live}, \
                         admission_closed={admission_closed}: {}",
                        witness.to_canonical_json()
                    );
                }
            }
        }
        assert_eq!(cases, 18);
    }

    /// G4: a real fallible root-backing refusal seals before launch, records
    /// explicit charge rollback, and leaves the same pool reusable.
    #[test]
    fn completion_witness_root_allocation_refusal_seals_and_pool_reuses() {
        let p = pool(2);
        refuse_next_root_allocation("slot-table");
        let (result, report, witness) = witnessed_parts!(p.run_witnessed(&SumKernel { tiles: 8 }));

        assert!(matches!(
            result,
            Err(RunError::MemoryAllocationRefused {
                what: "slot-table",
                ..
            })
        ));
        assert_eq!(report.completed, 0);
        assert_eq!(report.total, 8);
        assert_eq!(witness.verify(), Ok(()));
        assert_eq!(
            witness.disposition(),
            TilePoolCompletionDisposition::MemoryAllocationRefused
        );
        assert!(!witness.admission_completed());
        assert_eq!(witness.admitted_tiles(), 0);
        assert_eq!(witness.unadmitted_tiles(), 8);
        assert_eq!(witness.claimed_tiles(), 0);
        assert_eq!(witness.launched_workers(), 0);
        assert_eq!(witness.joined_workers(), 0);
        assert!(witness.root_charge_admitted());
        assert!(witness.root_charge_released());
        assert_eq!(witness.lease_used_after(), witness.lease_used_before());
        assert!(witness.executor_transients_quiescent());
        println!(
            "{}",
            witness.to_jsonl("root-backing-allocation-refusal", 0, "inline-private-hook")
        );

        assert_eq!(
            p.run(&SumKernel { tiles: 8 })
                .expect("pool remains reusable"),
            (1..=8).sum::<u64>()
        );
        assert!(p.arena_pool().stats().quiescent());
    }

    /// G4: a refusal after one real worker was admitted closes admission,
    /// requests drain, joins that worker, and seals the structured launch
    /// failure without fabricating an entry for the refused worker.
    #[test]
    fn completion_witness_worker_spawn_refusal_joins_started_workers() {
        let p = pool(2);
        refuse_next_worker_spawn(1);
        let (result, report, witness) =
            witnessed_parts!(p.run_witnessed(&SumKernel { tiles: 16_384 }));

        assert!(matches!(
            result,
            Err(RunError::WorkerSpawn {
                worker: 1,
                ref message,
                ..
            }) if message == "test-injected worker spawn refusal"
        ));
        assert_eq!(witness.verify(), Ok(()));
        assert_eq!(
            witness.disposition(),
            TilePoolCompletionDisposition::WorkerSpawnRefused
        );
        assert!(witness.admission_completed());
        assert_eq!(witness.planned_workers(), 2);
        assert_eq!(witness.launched_workers(), 1);
        assert_eq!(witness.joined_workers(), 1);
        assert!(witness.worker_admission_closed());
        assert!(witness.cancellation_requested());
        assert_eq!(
            witness.claimed_tiles(),
            witness.completed_tiles() + witness.break_tiles() + witness.panicked_tiles()
        );
        assert_eq!(report.completed, witness.completed_tiles());
        assert!(witness.executor_transients_quiescent());
        assert_eq!(
            p.run(&SumKernel { tiles: 8 })
                .expect("pool remains reusable after spawn refusal"),
            (1..=8).sum::<u64>()
        );
        assert!(p.arena_pool().stats().quiescent());
        println!(
            "{}",
            witness.to_jsonl_with_reuse(
                "worker-spawn-refusal",
                0,
                "inline-private-hook",
                Some(true),
            )
        );
    }

    #[test]
    fn parked_crew_overlap_is_refused_before_dispatch_and_then_reuses() {
        let p = pool(4);
        p.with_parked_crew_local(|parked| {
            let entered = std::sync::Arc::new(std::sync::Barrier::new(2));
            let release = std::sync::Arc::new(std::sync::Barrier::new(2));

            std::thread::scope(|scope| {
                let active_kernel = BlockingParkedKernel {
                    entered: std::sync::Arc::clone(&entered),
                    release: std::sync::Arc::clone(&release),
                };
                let active = scope.spawn(move || {
                    parked
                        .run_witnessed(&active_kernel)
                        .expect("active parked run must seal")
                });

                entered.wait();
                // Capture, do not assert, while this thread is the barrier
                // partner: a panic here would otherwise leave the active
                // run's worker blocked at `release` forever, turning a
                // clean failure into a suite-hanging join deadlock (kh5tf).
                let busy = parked.run_witnessed(&NoAllocation);
                release.wait();
                let active = active.join().expect("active parked caller");
                let busy = busy.expect("overlap refusal must seal a valid bundle");
                busy.verify_bundle().expect("busy refusal bundle");
                assert!(matches!(
                    busy.outcome(),
                    Err(RunError::ParkedCrewBusy {
                        kernel: "test/no-allocation"
                    })
                ));
                assert_eq!(
                    busy.witness().disposition(),
                    TilePoolCompletionDisposition::ParkedCrewBusy
                );
                assert!(!busy.witness().admission_completed());
                assert_eq!(busy.witness().root_metadata_bytes(), 0);
                assert!(!busy.witness().root_charge_admitted());
                assert_eq!(busy.witness().planned_crew_callbacks(), 4);
                assert_eq!(busy.witness().entered_crew_callbacks(), 0);
                assert_eq!(busy.witness().exited_crew_callbacks(), 0);

                active.verify_bundle().expect("active parked bundle");
                assert_eq!(active.outcome(), &Ok(1));
                assert_eq!(active.witness().planned_workers(), 1);
                assert_eq!(active.witness().launched_workers(), 1);
                assert_eq!(active.witness().joined_workers(), 1);
                assert_eq!(active.witness().planned_crew_callbacks(), 4);
                assert_eq!(active.witness().entered_crew_callbacks(), 4);
                assert_eq!(active.witness().exited_crew_callbacks(), 4);

                let mut missing_callback = active.witness().clone();
                missing_callback.entered_crew_callbacks = 3;
                missing_callback.root = completion_witness_root(&missing_callback);
                assert_completion_invariant(&missing_callback, "parked-crew-callback-drain");
            });

            let reused = parked
                .run_witnessed(&NoAllocation)
                .expect("crew must remain reusable after overlap refusal");
            assert_eq!(reused.outcome(), &Ok(1));
            assert_eq!(reused.witness().entered_crew_callbacks(), 4);
            assert_eq!(reused.witness().exited_crew_callbacks(), 4);
        });
        assert!(p.arena_pool().stats().quiescent());
    }

    #[test]
    fn parked_crew_recursive_dispatch_is_refused_without_poisoning_outer_run() {
        let p = pool(3);
        p.with_parked_crew_local(|parked| {
            let kernel = RecursiveParkedKernel {
                parked,
                nested: Mutex::new(None),
            };
            let outer = parked
                .run_witnessed(&kernel)
                .expect("outer recursive-dispatch fixture must seal");
            outer.verify_bundle().expect("outer bundle");
            assert_eq!(outer.outcome(), &Ok(1));
            assert_eq!(outer.witness().planned_crew_callbacks(), 3);
            assert_eq!(outer.witness().entered_crew_callbacks(), 3);
            assert_eq!(outer.witness().exited_crew_callbacks(), 3);

            let nested = kernel
                .nested
                .lock()
                .expect("recursive result")
                .take()
                .expect("kernel must retain its nested refusal");
            nested.verify_bundle().expect("nested busy bundle");
            assert!(matches!(
                nested.outcome(),
                Err(RunError::ParkedCrewBusy {
                    kernel: "test/no-allocation"
                })
            ));
            assert_eq!(nested.witness().planned_crew_callbacks(), 3);
            assert_eq!(nested.witness().entered_crew_callbacks(), 0);
            assert_eq!(nested.witness().exited_crew_callbacks(), 0);

            assert_eq!(
                parked
                    .run(&NoAllocation)
                    .expect("crew must remain reusable after recursive refusal"),
                1
            );
        });
        assert!(p.arena_pool().stats().quiescent());
    }
}
