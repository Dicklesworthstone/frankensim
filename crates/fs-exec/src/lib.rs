//! fs-exec — Two-lane executor: latency lane on asupersync, work-stealing
//! tile lane, `Cx` contract. Layer: L0.
//!
//! The beating heart of L0 (plan §5.2), in two lanes:
//!
//! - **Latency lane** ([`LatencyLane`]): asupersync's native async
//!   scheduling for orchestration — ledger I/O, IR interpretation, progress
//!   streaming. Thin by design; asupersync's request → drain → finalize
//!   cancellation protocol applies unmodified.
//! - **Throughput lane** ([`TilePool`]): a work-stealing fork-join pool
//!   whose units of work are TILES. Per-worker deques seeded with
//!   weight-proportional contiguous tile runs (the P/E-asymmetry hook),
//!   CCD-local-first steal order derived from fs-substrate topology, and
//!   FIXED-SLOT reductions folded in tile order — results are bit-identical
//!   across worker counts and steal schedules by construction (Decalogue
//!   P2).
//!
//! The contract every kernel programs against (plan Appendix B):
//! [`TileKernel`] (`tiles() -> TilePlan`, `run(tile, &Cx) ->
//! ControlFlow<Cancelled, Out>` with `Out: Reduce`) and [`Cx`] — the
//! cancellation gate polled at tile boundaries, the tile-scoped fs-alloc
//! arena, the [`StreamKey`] RNG identity keyed by `(seed, kernel_id, tile,
//! iteration)` (never by worker), the budget slice, and the [`ExecMode`].
//!
//! Failure containment: a panicking tile is caught, siblings drain via the
//! gate, and the run returns a structured [`RunError`] with full tile
//! provenance — never a process abort mid-campaign (Decalogue P10/P7).
//!
//! On top of the lanes, the three executor behaviors most HPC runtimes
//! lack (plan §5.2): speculative races with loser-cancellation and
//! deterministic victory ([`Racer`]), resumable/forkable solvers with
//! bit-exact pause-serialize-resume ([`solver`]), and the statistical-
//! preemption kill-handle registry ([`KillRegistry`], Bet 8's machinery).
//!
//! See CONTRACT.md for invariants, determinism class, cancellation
//! behavior, and no-claim boundaries.

mod cx;
mod kernel;
mod kill;
mod latency;
mod pool;
mod race;
pub mod reduce;
pub mod solver;

pub use cx::{CancelGate, Cancelled, Cx, ExecMode, StreamKey};
pub use kernel::{Reduce, TileKernel, TilePlan};
pub use kill::{CandidateId, KillRegistry};
pub use latency::{LaneError, LatencyLane};
pub use pool::{PoolConfig, RunError, RunReport, TilePool, victim_order, weighted_ranges};
pub use race::{BranchOutcome, BranchReport, NoWinner, RaceBranch, RaceRun, Racer, RacerConfig};

/// Crate version, re-exported for provenance stamping (the Five Explicits'
/// "versions" pillar reaches down to individual crates).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
