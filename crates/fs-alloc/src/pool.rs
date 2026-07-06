//! Sharded object pools for recurring same-shape allocations (tiles, small
//! matrices, packing buffers — plan §5.3).
//!
//! Sharding is the per-CCD hook: callers pass a shard index (fs-exec will
//! pass the worker's CCD/core-class id), each shard has its own lock and
//! free list padded to the 128-byte policy so shards never false-share, and
//! an object is constructed by the ACQUIRING worker (`make` runs on the
//! caller — the first-touch placement hook from plan §5.1 consequence 3).

use core::fmt;
use core::ops::{Deref, DerefMut};
use std::fmt::Write as _;
use std::sync::Mutex;

use crate::CachePadded;

struct Shard<T> {
    free: Vec<T>,
    created: u64,
    recycled: u64,
    detached: u64,
    live: u64,
}

impl<T> Shard<T> {
    const fn new() -> Self {
        Shard {
            free: Vec::new(),
            created: 0,
            recycled: 0,
            detached: 0,
            live: 0,
        }
    }
}

/// A fixed-shard object pool. `Sync`: shards are independently locked.
pub struct ShardedPool<T> {
    shards: Box<[CachePadded<Mutex<Shard<T>>>]>,
}

impl<T> ShardedPool<T> {
    /// Build a pool with `shards` shards (clamped to at least 1). Size it to
    /// the machine's CCD/cluster count (fs-substrate's probe) for the
    /// intended locality effect.
    #[must_use]
    pub fn new(shards: usize) -> Self {
        let shards = shards.max(1);
        ShardedPool {
            shards: (0..shards)
                .map(|_| CachePadded::new(Mutex::new(Shard::new())))
                .collect(),
        }
    }

    /// Number of shards.
    #[must_use]
    pub fn shard_count(&self) -> usize {
        self.shards.len()
    }

    /// Take an object from `shard` (index taken modulo the shard count),
    /// constructing one with `make` if the shard's free list is empty.
    /// `make` runs on the calling thread — that is the first-touch hook.
    /// The item returns to ITS shard when the guard drops.
    pub fn acquire_with(&self, shard: usize, make: impl FnOnce() -> T) -> PoolItem<'_, T> {
        let home = &self.shards[shard % self.shards.len()];
        let value = {
            let mut s = home.get().lock().expect("fs-alloc shard poisoned");
            s.live += 1;
            match s.free.pop() {
                Some(v) => {
                    s.recycled += 1;
                    Some(v)
                }
                None => {
                    s.created += 1;
                    None
                }
            }
        };
        // Construction happens OUTSIDE the lock, on the acquiring thread.
        let value = value.unwrap_or_else(make);
        PoolItem {
            value: Some(value),
            home: home.get(),
        }
    }

    /// Deterministic per-shard accounting (index order).
    #[must_use]
    pub fn stats(&self) -> ShardedPoolStats {
        ShardedPoolStats {
            shards: self
                .shards
                .iter()
                .map(|padded| {
                    let s = padded.get().lock().expect("fs-alloc shard poisoned");
                    ShardStats {
                        created: s.created,
                        recycled: s.recycled,
                        detached: s.detached,
                        live: s.live,
                        free: s.free.len() as u64,
                    }
                })
                .collect(),
        }
    }

    /// True when every acquired item has been returned or detached — the
    /// pool-side leak oracle for G4 storms.
    #[must_use]
    pub fn quiescent(&self) -> bool {
        self.stats().shards.iter().all(|s| s.live == 0)
    }
}

impl<T> fmt::Debug for ShardedPool<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ShardedPool")
            .field("stats", &self.stats().to_json())
            .finish_non_exhaustive()
    }
}

/// RAII guard for a pooled object; derefs to `T` and returns the object to
/// its home shard on drop.
pub struct PoolItem<'p, T> {
    value: Option<T>,
    home: &'p Mutex<Shard<T>>,
}

impl<T> PoolItem<'_, T> {
    /// Permanently remove the object from the pool ("graduation" of a
    /// long-lived artifact — plan §5.3). The pool records the detachment;
    /// the object no longer counts as live.
    #[must_use]
    pub fn detach(mut self) -> T {
        let value = self.value.take().expect("value present until drop/detach");
        let mut s = self.home.lock().expect("fs-alloc shard poisoned");
        s.live -= 1;
        s.detached += 1;
        value
    }
}

impl<T> Deref for PoolItem<'_, T> {
    type Target = T;

    fn deref(&self) -> &T {
        self.value.as_ref().expect("value present until drop/detach")
    }
}

impl<T> DerefMut for PoolItem<'_, T> {
    fn deref_mut(&mut self) -> &mut T {
        self.value.as_mut().expect("value present until drop/detach")
    }
}

impl<T> Drop for PoolItem<'_, T> {
    fn drop(&mut self) {
        if let Some(value) = self.value.take() {
            let mut s = self.home.lock().expect("fs-alloc shard poisoned");
            s.live -= 1;
            s.free.push(value);
        }
    }
}

/// Accounting for one shard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShardStats {
    /// Objects constructed (free-list misses).
    pub created: u64,
    /// Acquisitions served from the free list.
    pub recycled: u64,
    /// Objects permanently removed via [`PoolItem::detach`].
    pub detached: u64,
    /// Objects currently checked out.
    pub live: u64,
    /// Objects parked in the free list.
    pub free: u64,
}

/// Deterministic pool accounting (shards in index order).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShardedPoolStats {
    /// Per-shard stats, index order.
    pub shards: Vec<ShardStats>,
}

impl ShardedPoolStats {
    /// Canonical JSON object (deterministic field and shard order).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::from("{\"shards\":[");
        for (i, s) in self.shards.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"created\":{},\"recycled\":{},\"detached\":{},\"live\":{},\"free\":{}}}",
                s.created, s.recycled, s.detached, s.live, s.free
            );
        }
        out.push_str("]}");
        out
    }

    /// Package as an `fs-obs` event payload for the ledger pipeline.
    #[must_use]
    pub fn to_event_kind(&self) -> fs_obs::EventKind {
        fs_obs::EventKind::Custom {
            name: "fs-alloc-sharded-pool-stats".to_string(),
            json: self.to_json(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acquire_reuse_and_detach_account_exactly() {
        let pool: ShardedPool<Vec<u8>> = ShardedPool::new(2);
        {
            let mut a = pool.acquire_with(0, || vec![0u8; 64]);
            a[0] = 7;
            assert_eq!(a.len(), 64);
        } // returned to shard 0
        {
            let b = pool.acquire_with(0, || vec![0u8; 64]);
            assert_eq!(b[0], 7, "same buffer must come back (LIFO reuse)");
            let owned = b.detach();
            assert_eq!(owned.len(), 64);
        } // detached: not returned
        let stats = pool.stats();
        assert_eq!(stats.shards[0].created, 1);
        assert_eq!(stats.shards[0].recycled, 1);
        assert_eq!(stats.shards[0].detached, 1);
        assert_eq!(stats.shards[0].free, 0);
        assert_eq!(stats.shards[1].created, 0);
        assert!(pool.quiescent());
        assert_eq!(
            stats.to_json(),
            "{\"shards\":[{\"created\":1,\"recycled\":1,\"detached\":1,\"live\":0,\"free\":0},\
             {\"created\":0,\"recycled\":0,\"detached\":0,\"live\":0,\"free\":0}]}"
        );
    }

    #[test]
    fn shard_index_wraps_and_zero_shards_clamp() {
        let pool: ShardedPool<u32> = ShardedPool::new(0);
        assert_eq!(pool.shard_count(), 1);
        let x = pool.acquire_with(17, || 5);
        assert_eq!(*x, 5);
    }

    #[test]
    fn shards_do_not_false_share_by_construction() {
        // Structural check for the padding policy: each shard slot is
        // 128-byte sized/aligned via CachePadded.
        assert_eq!(align_of::<CachePadded<Mutex<Shard<u64>>>>() % 128, 0);
        assert_eq!(size_of::<CachePadded<Mutex<Shard<u64>>>>() % 128, 0);
    }

    #[test]
    fn parallel_hammer_stays_quiescent() {
        let pool: ShardedPool<Vec<u64>> = ShardedPool::new(4);
        std::thread::scope(|s| {
            let pool = &pool;
            for t in 0..8 {
                s.spawn(move || {
                    for i in 0..500 {
                        let mut buf = pool.acquire_with(t * 31 + i, || vec![0u64; 32]);
                        buf[i % 32] = i as u64;
                    }
                });
            }
        });
        assert!(pool.quiescent(), "{}", pool.stats().to_json());
        let total_created: u64 = pool.stats().shards.iter().map(|s| s.created).sum();
        assert!(total_created <= 8 * 4, "reuse must bound construction");
    }
}
