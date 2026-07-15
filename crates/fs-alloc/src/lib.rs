//! fs-alloc — Scope arenas with O(1) cancel reclaim, 128-byte alignment,
//! hugepage-eligible chunks, and sharded object pools. Layer: L0.
//!
//! The memory discipline of plan §5.3, as types:
//!
//! - [`ArenaPool`] + [`Arena`]: bump arenas scoped 1:1 to units of work.
//!   Completion or cancellation drops the arena and reclaims EVERYTHING at
//!   a cost independent of allocation count (Decalogue P7: cancellation
//!   reclaims arenas without walking them). Cross-scope escapes are compile
//!   errors — lifetimes do the enforcement (see the `compile_fail` battery
//!   on [`ArenaPool::scope`] and [`Arena::alloc`]).
//! - [`ALLOC_ALIGN`] = 128 bytes UNCONDITIONALLY (superset of Apple's
//!   128-byte and x86-64's 64-byte cache lines): every arena allocation is
//!   128-byte aligned, and [`CachePadded`] pads shared slots so they never
//!   false-share.
//! - [`HugepagePolicy`]: chunks of 2 MiB and up become THP-*eligible* on
//!   Linux where the kernel allows it; the decision — including every
//!   fallback — is recorded in [`HugepageDecision`] (never claimed
//!   silently; see CONTRACT.md no-claims).
//! - [`ShardedPool`]: recycling pools for recurring same-shape allocations
//!   with per-CCD sharding hooks and first-touch construction.
//! - [`Site`] tags + [`SiteReport`]: allocation-site accounting whose JSON
//!   is deterministic and diffable between runs, feeding the Ledger through
//!   `fs-obs` events.
//! - [`ReclaimPoison`]: opt-in seeded G4 mode that poisons retained chunks,
//!   verifies them before reuse, and quarantines stale-write corruption with
//!   a structured receipt.
//!
//! Error model: every fallible path returns a structured, teaching
//! [`AllocError`] (Decalogue P10). Out-of-memory is a `Result`, not an
//! abort.
//!
//! Unsafe boundary: two registered capsules behind this crate's safe facade:
//! `src/raw/mod.rs` (bump-pointer core and chunk alloc/dealloc) and
//! `src/raw/poison/mod.rs` (diagnostic byte fill/scan). See each adjacent
//! SAFETY.md.
//!
//! See CONTRACT.md for invariants, determinism class, cancellation
//! behavior, and no-claim boundaries.

mod arena;
mod hugepage;
mod lease;
mod lease_vec;
mod pool;
mod raw;

pub use arena::{
    AllocError, Arena, ArenaConfig, ArenaPool, ArenaStats, PoolStats, RECLAIM_POISON_VERSION,
    ReclaimPoison, ReclaimPoisonMutation, Site, SiteReport, SiteStats,
};
pub use hugepage::{HUGEPAGE_BYTES, HugepageDecision, HugepageOutcome, HugepagePolicy};
pub use lease::{LeaseCharge, LeaseReceipt, LeaseRefusal, OperationMemoryLease};
pub use lease_vec::LeasedVec;
pub use pool::{PoolItem, ShardStats, ShardedPool, ShardedPoolStats};

/// Crate version, re-exported for provenance stamping (the Five Explicits'
/// "versions" pillar reaches down to individual crates).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The unconditional allocation alignment policy (plan §5.3): 128 bytes, the
/// superset of both reference targets' cache lines (Apple aarch64 128 B,
/// x86-64 64 B). Every arena allocation is aligned to this; getting it wrong
/// on M-series silently halves effective bandwidth on contended structures.
pub const ALLOC_ALIGN: usize = 128;

/// Pads and aligns `T` to [`ALLOC_ALIGN`] so adjacent instances (e.g. one
/// per worker/shard) never share a cache line on either reference target.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[repr(align(128))]
pub struct CachePadded<T>(T);

impl<T> CachePadded<T> {
    /// Wrap a value.
    #[must_use]
    pub const fn new(value: T) -> Self {
        CachePadded(value)
    }

    /// Shared access to the padded value.
    #[must_use]
    pub const fn get(&self) -> &T {
        &self.0
    }

    /// Exclusive access to the padded value.
    pub const fn get_mut(&mut self) -> &mut T {
        &mut self.0
    }

    /// Unwrap.
    #[must_use]
    pub fn into_inner(self) -> T {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn alignment_policy_is_128() {
        assert_eq!(ALLOC_ALIGN, 128);
        assert!(ALLOC_ALIGN.is_power_of_two());
    }

    #[test]
    fn cache_padded_pads_and_aligns() {
        assert_eq!(align_of::<CachePadded<u8>>(), 128);
        assert_eq!(size_of::<CachePadded<u8>>(), 128);
        assert_eq!(size_of::<CachePadded<[u8; 200]>>(), 256);
        let mut c = CachePadded::new(41u32);
        *c.get_mut() += 1;
        assert_eq!(*c.get(), 42);
        assert_eq!(c.into_inner(), 42);
    }
}
