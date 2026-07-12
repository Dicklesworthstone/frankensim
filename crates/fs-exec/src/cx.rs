//! The `Cx` contract: what every tile sees (plan §5.2, Appendix B).
//!
//! A `Cx` carries the cancellation gate (poll at tile boundaries), the
//! tile-scoped bump arena, the counter-based RNG stream key derived from
//! LOGICAL identity `(seed, kernel_id, tile, iteration)` — never from the
//! worker that happens to run the tile (Decalogue P2) — the budget slice,
//! and the execution mode. The ledger handle joins when fs-ledger lands;
//! until then accounting flows through fs-obs events.

pub use asupersync::types::Budget;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

/// Execution mode, part of every run's provenance (plan §5.4). Both modes
/// currently use the same fixed-shape slot reduction; `Fast` reserves the
/// right to relax reduction shape once the roofline harness can price it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExecMode {
    /// Bit-stable across runs and worker counts on the same ISA.
    #[default]
    Deterministic,
    /// May relax determinism for throughput; the mode is recorded in every
    /// event so a result always knows how it was made.
    Fast,
}

impl ExecMode {
    /// Stable lowercase name for events and ledger rows.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            ExecMode::Deterministic => "deterministic",
            ExecMode::Fast => "fast",
        }
    }
}

/// An explicit, caller-ledgered logical run identity (bead wf9.7.1):
/// the generation/iteration axis of [`StreamKey`]. Execution order is
/// NOT identity — a reused or concurrently shared pool must never
/// perturb streams, so the run id comes from the caller's ledger
/// (trial index, generation number, restart ordinal), never from a
/// pool-global counter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct RunId(pub u64);

impl RunId {
    /// Derive a deterministic child run for one phase/ordinal of this logical
    /// operation. The domain must name the producer and schema; derivation is
    /// independent of scheduling and pool history.
    #[must_use]
    pub fn derive(self, domain: &str, ordinal: u64) -> Self {
        let mut preimage = [0_u8; 16];
        preimage[..8].copy_from_slice(&self.0.to_le_bytes());
        preimage[8..].copy_from_slice(&ordinal.to_le_bytes());
        let hash = fs_blake3::hash_domain(domain, &preimage);
        let mut id = [0_u8; 8];
        id.copy_from_slice(&hash.as_bytes()[..8]);
        Self(u64::from_le_bytes(id))
    }
}

/// Logical RNG stream identity (plan §5.2): results must be independent of
/// which worker ran which tile, so streams are keyed by WHAT the work is,
/// not WHERE it ran. fs-rand's Philox generator consumes [`StreamKey::key128`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StreamKey {
    /// Study/run seed (the Five Explicits' seed pillar).
    pub seed: u64,
    /// Kernel identity (stable hash of the kernel name).
    pub kernel_id: u64,
    /// Logical tile identity within the kernel's plan.
    pub tile: u64,
    /// Iteration/generation counter for kernels re-run over the same tiles.
    pub iteration: u64,
}

impl StreamKey {
    /// Pack the logical identity into the 128-bit Philox key domain.
    /// Mixing is a fixed xor/multiply fold — deterministic across platforms,
    /// with each field influencing both halves.
    #[must_use]
    pub const fn key128(self) -> u128 {
        // SplitMix64-style finalizer per field, then lane packing.
        const fn mix(mut z: u64) -> u64 {
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        }
        let hi = mix(self.seed ^ mix(self.kernel_id));
        let lo = mix(self.tile ^ mix(self.iteration).rotate_left(17));
        ((hi as u128) << 64) | lo as u128
    }
}

/// Shared cancellation gate for one kernel run: request → drain → finalize.
/// Setting the gate is a REQUEST; workers drain (finish the current tile,
/// then stop claiming) and the run finalizes with structured outcomes —
/// never a silent drop (Decalogue P7).
///
/// The gate owns a monotonic origin so request/observation timestamps share
/// one domain for latency histograms. Timestamps feed REPORTS only — they
/// never influence results (determinism discipline, like fs-substrate's
/// measured bandwidth).
#[derive(Debug)]
pub struct CancelGate {
    requested: AtomicBool,
    origin: std::time::Instant,
    /// Nanoseconds after `origin` of the FIRST request. 0 = unset.
    requested_at_ns: AtomicU64,
}

impl Default for CancelGate {
    fn default() -> Self {
        CancelGate {
            requested: AtomicBool::new(false),
            origin: std::time::Instant::now(),
            requested_at_ns: AtomicU64::new(0),
        }
    }
}

impl CancelGate {
    /// Fresh, un-requested gate.
    #[must_use]
    pub fn new() -> Self {
        CancelGate::default()
    }

    /// Nanoseconds since this gate's origin (shared timestamp domain for
    /// latency accounting).
    #[must_use]
    pub fn now_ns(&self) -> u64 {
        self.origin.elapsed().as_nanos() as u64
    }

    /// Request cancellation (idempotent; the first request's timestamp is
    /// the one latency histograms measure from).
    pub fn request(&self) {
        let now = self.now_ns().max(1);
        let _ = self
            .requested_at_ns
            .compare_exchange(0, now, Ordering::AcqRel, Ordering::Acquire);
        self.requested.store(true, Ordering::Release);
    }

    /// Cheap synchronous poll (the tile-boundary check).
    #[must_use]
    pub fn is_requested(&self) -> bool {
        self.requested.load(Ordering::Acquire)
    }

    /// Timestamp of the first request (ns after origin), if any.
    #[must_use]
    pub fn requested_at_ns(&self) -> Option<u64> {
        match self.requested_at_ns.load(Ordering::Acquire) {
            0 => None,
            t => Some(t),
        }
    }
}

/// Marker returned by kernels that observed cancellation at a tile boundary
/// (plan Appendix B: `run(...) -> ControlFlow<Cancelled, Out>`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cancelled;

/// A typed tile-local refusal that requests sibling drain without pretending
/// the condition was either ordinary cancellation or a panic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TileFailure {
    /// The tile's scoped arena refused an allocation.
    Allocation(fs_alloc::AllocError),
}

impl core::fmt::Display for TileFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Allocation(error) => write!(f, "tile arena allocation refused: {error}"),
        }
    }
}

impl core::error::Error for TileFailure {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Allocation(error) => Some(error),
        }
    }
}

#[derive(Default)]
pub(crate) struct RefusalSink {
    lowest: Mutex<Option<(u64, TileFailure)>>,
}

impl RefusalSink {
    fn record(&self, tile: u64, failure: TileFailure) {
        let mut lowest = self
            .lowest
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if lowest
            .as_ref()
            .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
        {
            *lowest = Some((tile, failure));
        }
    }

    pub(crate) fn take(&self) -> Option<(u64, TileFailure)> {
        self.lowest
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .take()
    }
}

/// The per-tile execution context handed to [`crate::TileKernel::run`].
/// Lifetime-scoped: the arena (and everything allocated from it) cannot
/// outlive the tile's scope — fs-alloc's lifetime discipline does the
/// enforcement.
pub struct Cx<'s> {
    gate: &'s CancelGate,
    arena: &'s fs_alloc::Arena,
    key: StreamKey,
    budget: Budget,
    mode: ExecMode,
    refusals: Option<&'s RefusalSink>,
    lease: Option<&'s fs_alloc::OperationMemoryLease>,
}

impl<'s> Cx<'s> {
    /// Assemble a context (the executor's job; kernels only consume it).
    #[must_use]
    pub fn new(
        gate: &'s CancelGate,
        arena: &'s fs_alloc::Arena,
        key: StreamKey,
        budget: Budget,
        mode: ExecMode,
    ) -> Self {
        Cx {
            gate,
            arena,
            key,
            budget,
            mode,
            refusals: None,
            lease: None,
        }
    }

    pub(crate) fn new_with_refusal_sink(
        gate: &'s CancelGate,
        arena: &'s fs_alloc::Arena,
        key: StreamKey,
        budget: Budget,
        mode: ExecMode,
        refusals: &'s RefusalSink,
        lease: &'s fs_alloc::OperationMemoryLease,
    ) -> Self {
        Cx {
            gate,
            arena,
            key,
            budget,
            mode,
            refusals: Some(refusals),
            lease: Some(lease),
        }
    }

    /// The run's operation memory lease, when this context came from a
    /// pool run (bead wf9.16.1). Kernels admit output payload storage
    /// through it (`fs_alloc::LeasedVec::with_capacity(cx.lease()?, ..)`);
    /// manually constructed contexts have none.
    #[must_use]
    pub fn lease(&self) -> Option<&'s fs_alloc::OperationMemoryLease> {
        self.lease
    }

    /// Poll the cancellation gate: the MANDATORY call at tile boundaries
    /// (and at bounded strides inside long tiles). Returns `Err(Cancelled)`
    /// when a request is pending; the kernel converts that into
    /// `ControlFlow::Break(Cancelled)` and returns promptly.
    ///
    /// # Errors
    /// [`Cancelled`] when cancellation has been requested.
    pub fn checkpoint(&self) -> Result<(), Cancelled> {
        if self.gate.is_requested() {
            Err(Cancelled)
        } else {
            Ok(())
        }
    }

    /// Non-consuming form of the poll for `while !cx.is_cancel_requested()`
    /// loop shapes.
    #[must_use]
    pub fn is_cancel_requested(&self) -> bool {
        self.gate.is_requested()
    }

    /// Record a typed tile refusal, request sibling drain, and return the
    /// marker expected by [`crate::TileKernel::run`]'s `ControlFlow::Break`.
    ///
    /// Pool-created contexts surface the lowest observed logical tile as
    /// [`crate::RunError::TileFailed`]. Manually constructed contexts have no
    /// pool refusal sink, but still request their gate so the failure cannot be
    /// ignored accidentally.
    #[must_use]
    pub fn refuse(&self, failure: TileFailure) -> Cancelled {
        if let Some(refusals) = self.refusals {
            refusals.record(self.key.tile, failure);
        }
        self.gate.request();
        Cancelled
    }

    /// The tile-scoped bump arena (O(chunks) reclaim at scope end; escapes
    /// are compile errors — see fs-alloc's contract).
    #[must_use]
    pub fn arena(&self) -> &'s fs_alloc::Arena {
        self.arena
    }

    /// This tile's logical RNG stream identity.
    #[must_use]
    pub fn stream_key(&self) -> StreamKey {
        self.key
    }

    /// The budget slice this tile runs under (asupersync's vocabulary;
    /// enforcement beyond cancellation lands with the session governor —
    /// see CONTRACT no-claims).
    #[must_use]
    pub fn budget(&self) -> Budget {
        self.budget
    }

    /// The execution mode (provenance, P2).
    #[must_use]
    pub fn mode(&self) -> ExecMode {
        self.mode
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn child_run_ids_bind_parent_domain_and_ordinal() {
        let parent = RunId(17);
        let first = parent.derive("org.frankensim.test.phase-a.v1", 0);
        assert_eq!(first, parent.derive("org.frankensim.test.phase-a.v1", 0));
        assert_ne!(first, parent.derive("org.frankensim.test.phase-a.v1", 1));
        assert_ne!(first, parent.derive("org.frankensim.test.phase-b.v1", 0));
        assert_ne!(first, RunId(18).derive("org.frankensim.test.phase-a.v1", 0));
    }

    #[test]
    fn stream_keys_depend_on_every_logical_field_and_nothing_else() {
        let base = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let k0 = base.key128();
        assert_eq!(k0, base.key128(), "pure function of logical identity");
        for (i, variant) in [
            StreamKey { seed: 9, ..base },
            StreamKey {
                kernel_id: 9,
                ..base
            },
            StreamKey { tile: 9, ..base },
            StreamKey {
                iteration: 9,
                ..base
            },
        ]
        .iter()
        .enumerate()
        {
            assert_ne!(k0, variant.key128(), "field {i} must matter");
        }
    }

    #[test]
    fn gate_records_first_request_and_polls_cheaply() {
        let gate = CancelGate::new();
        assert!(!gate.is_requested());
        assert_eq!(gate.requested_at_ns(), None);
        gate.request();
        let first = gate.requested_at_ns().expect("stamped");
        gate.request(); // later request must not overwrite the first stamp
        assert!(gate.is_requested());
        assert_eq!(gate.requested_at_ns(), Some(first));
        assert!(gate.now_ns() >= first);
    }

    #[test]
    fn cx_checkpoint_observes_the_gate() {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(cx.checkpoint().is_ok());
            gate.request();
            assert_eq!(cx.checkpoint(), Err(Cancelled));
            assert!(cx.is_cancel_requested());
            assert_eq!(cx.mode().name(), "deterministic");
        });
    }
}
