//! Scope arenas: the safe facade over the bump-pointer capsule.
//!
//! The model (plan §5.3): an [`ArenaPool`] owns configuration, accounting,
//! and a chunk free list; each unit of scoped work gets its own [`Arena`]
//! (bump allocator). Completing or cancelling the scope drops the arena,
//! which reclaims every allocation at once — cost proportional to the
//! CHUNK count (O(log bytes) from geometric growth), independent of the
//! allocation count — and recycles the chunks for the next scope. Escapes
//! are compile errors: allocations borrow the arena, and the closure-based
//! scope APIs are higher-ranked so nothing tied to the scope's lifetime can
//! leave it.
//!
//! Binding arenas 1:1 to *asupersync* scopes is fs-exec's contract (its
//! `Cx` will carry `&'scope Arena`); fs-alloc supplies the lifetime
//! discipline and the O(1)-per-chunk reclaim primitive that make that
//! binding leak-free under cancellation (G4).

use core::cell::{Cell, RefCell};
use core::fmt;
use core::mem::{align_of, size_of};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use crate::ALLOC_ALIGN;
use crate::hugepage::{HugepageDecision, HugepagePolicy};
use crate::raw::{Chunk, RawArena};

/// An allocation-site tag: a static name under which bytes are accounted,
/// so memory usage feeds the Ledger and regressions are diffable between
/// runs (plan §5.3). Cheap to copy; compare by name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Site(&'static str);

impl Site {
    /// Tag for a static site name (convention: `"crate/what"`, e.g.
    /// `"fs-la/gemm-packing"`).
    #[must_use]
    pub const fn named(name: &'static str) -> Self {
        Site(name)
    }

    /// The site name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.0
    }
}

/// Per-site accounting: cumulative padded payload bytes and allocation count.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SiteStats {
    /// Cumulative 128-byte-padded payload bytes allocated under this site.
    /// Saturates at `u64::MAX` rather than wrapping.
    pub bytes: u64,
    /// Cumulative allocation count under this site.
    /// Saturates at `u64::MAX` rather than wrapping.
    pub allocations: u64,
}

/// Pool configuration. Values are normalized (not rejected) by
/// [`ArenaPool::new`]: `chunk_bytes` is clamped to at least 4 KiB and
/// `max_chunk_bytes` to at least `chunk_bytes`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArenaConfig {
    /// First-chunk size for each arena; chunks grow geometrically from here.
    pub chunk_bytes: usize,
    /// Ceiling for the geometric chunk growth.
    pub max_chunk_bytes: usize,
    /// Budget on total OS-reserved bytes (in-use + free-listed) across the
    /// pool. `None` = unlimited. Exceeding it yields a structured
    /// [`AllocError::Exhausted`], never an abort (Decalogue P4/P10).
    pub limit_bytes: Option<usize>,
    /// Cap on bytes parked in the chunk free list; chunks beyond it are
    /// returned to the OS on arena drop.
    pub free_list_max_bytes: usize,
    /// Hugepage intent; the probe outcome is recorded per pool.
    pub hugepage: HugepagePolicy,
}

impl Default for ArenaConfig {
    fn default() -> Self {
        ArenaConfig {
            chunk_bytes: 1 << 20,
            max_chunk_bytes: 32 << 20,
            limit_bytes: None,
            free_list_max_bytes: 64 << 20,
            hugepage: HugepagePolicy::default(),
        }
    }
}

/// Structured allocation failure (Decalogue P10: errors are guidance).
/// Never a panic, never an abort.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AllocError {
    /// The pool's `limit_bytes` budget cannot cover a needed chunk.
    Exhausted {
        /// Site requesting the allocation.
        site: &'static str,
        /// Chunk bytes that were requested.
        requested_bytes: usize,
        /// OS-reserved bytes at refusal time (in-use + free list).
        reserved_bytes: usize,
        /// The configured budget.
        limit_bytes: usize,
    },
    /// The global allocator refused the chunk request.
    OutOfMemory {
        /// Site requesting the allocation.
        site: &'static str,
        /// Chunk bytes that were requested.
        requested_bytes: usize,
    },
    /// The operation memory lease refused the chunk (bead wf9.16). The
    /// pool budget may still have room — this is the per-operation gate.
    LeaseExhausted {
        /// Site requesting the allocation.
        site: &'static str,
        /// Chunk bytes the lease was asked for.
        requested_bytes: u64,
        /// Lease bytes in use at refusal time.
        used_bytes: u64,
        /// The lease limit in force.
        limit_bytes: u64,
    },
    /// `len * size_of::<T>()` overflows `usize`.
    LayoutOverflow {
        /// Site requesting the allocation.
        site: &'static str,
        /// Requested element count.
        len: usize,
        /// Element size in bytes.
        elem_bytes: usize,
    },
    /// Arithmetic needed to size or account for a chunk exceeds `usize`.
    ReservationOverflow {
        /// Site requesting the reservation.
        site: &'static str,
        /// Existing bytes participating in the failed sum.
        base_bytes: usize,
        /// Additional bytes participating in the failed sum.
        additional_bytes: usize,
    },
}

impl fmt::Display for AllocError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AllocError::Exhausted {
                site,
                requested_bytes,
                reserved_bytes,
                limit_bytes,
            } => write!(
                f,
                "arena budget exhausted at site `{site}`: requested {requested_bytes} B for a \
                 new chunk with {reserved_bytes} B of the {limit_bytes} B limit already \
                 reserved. Fixes (ranked): (1) raise ArenaConfig::limit_bytes; (2) end or \
                 shrink concurrent scopes holding arenas; (3) lower ArenaConfig::chunk_bytes \
                 so growth is finer-grained"
            ),
            AllocError::OutOfMemory {
                site,
                requested_bytes,
            } => write!(
                f,
                "global allocator refused {requested_bytes} B at site `{site}`: the process \
                 is near its memory ceiling. Fixes (ranked): (1) set ArenaConfig::limit_bytes \
                 so pressure degrades to a budget error earlier; (2) split the workload into \
                 smaller scopes"
            ),
            AllocError::LeaseExhausted {
                site,
                requested_bytes,
                used_bytes,
                limit_bytes,
            } => write!(
                f,
                "operation memory lease exhausted at site `{site}`: requested {requested_bytes} B \
                 for a chunk with {used_bytes} B of the {limit_bytes} B operation lease in use. \
                 Fixes (ranked): (1) raise the operation lease; (2) reduce concurrent tile \
                 working sets; (3) lower ArenaConfig::chunk_bytes so growth is finer-grained"
            ),
            AllocError::LayoutOverflow {
                site,
                len,
                elem_bytes,
            } => write!(
                f,
                "slice layout overflows at site `{site}`: {len} elements x {elem_bytes} B \
                 exceeds the address space; validate the length upstream"
            ),
            AllocError::ReservationOverflow {
                site,
                base_bytes,
                additional_bytes,
            } => write!(
                f,
                "arena reservation arithmetic overflows at site `{site}`: {base_bytes} base B + \
                 {additional_bytes} additional B exceeds the address space; reduce the allocation \
                 or arena chunk size"
            ),
        }
    }
}

impl core::error::Error for AllocError {}

/// Chunk free list; byte total lives under the same lock as the chunks so
/// they can never disagree.
struct FreeList {
    chunks: Vec<Chunk>,
    bytes: usize,
}

/// State shared by a pool and all its arenas.
struct PoolShared {
    config: ArenaConfig,
    hugepage: HugepageDecision,
    free: Mutex<FreeList>,
    /// Total OS-reserved bytes (in-use + free-listed).
    reserved_bytes: AtomicUsize,
    arenas_live: AtomicUsize,
    chunks_created: AtomicU64,
    chunks_recycled: AtomicU64,
    sites: Mutex<BTreeMap<&'static str, SiteStats>>,
}

/// A chunk admitted by both ledgers but not yet installed in an arena.
/// Until `install_into` completes, dropping this value releases the operation
/// charge first and returns the chunk to the pool. This closes the unwind gap
/// between acquisition and arena ownership transfer.
struct ChunkAcquisition<'a> {
    shared: &'a PoolShared,
    chunk: Option<Chunk>,
    charge: Option<crate::LeaseCharge>,
}

impl ChunkAcquisition<'_> {
    fn len(&self) -> usize {
        self.chunk.as_ref().expect("pending chunk is present").len()
    }

    fn install_into(mut self, arena: &Arena, next_chunk_bytes: usize) {
        let chunk_bytes = self.len() as u64;
        let next_leased_bytes = self.charge.as_ref().map(|_| {
            arena
                .leased_bytes
                .get()
                .checked_add(chunk_bytes)
                .expect("live arena chunks cannot exceed the address space")
        });
        let chunk = self.chunk.take().expect("pending chunk is present");
        arena.raw.install_chunk(chunk);
        arena.next_chunk_bytes.set(next_chunk_bytes);
        if let Some(next) = next_leased_bytes {
            arena.leased_bytes.set(next);
        }
        if let Some(charge) = self.charge.take() {
            charge.commit_to_manual_release();
        }
    }
}

impl Drop for ChunkAcquisition<'_> {
    fn drop(&mut self) {
        // A rejected/unwinding acquisition must stop belonging to the
        // operation before another operation can recycle the chunk.
        drop(self.charge.take());
        if let Some(chunk) = self.chunk.take() {
            self.shared.release_chunk(chunk);
        }
    }
}

impl PoolShared {
    fn claim_new_chunk_bytes(&self, want: usize, site: Site) -> Result<(), AllocError> {
        loop {
            let reserved = self.reserved_bytes.load(Ordering::Acquire);
            let next = reserved
                .checked_add(want)
                .ok_or(AllocError::ReservationOverflow {
                    site: site.name(),
                    base_bytes: reserved,
                    additional_bytes: want,
                })?;
            if self.config.limit_bytes.is_some_and(|limit| next > limit) {
                // Cached chunks count against the hard pool limit. Return them
                // to the OS once before deciding that this request cannot fit.
                let mut free = self
                    .free
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                let drained = core::mem::take(&mut free.chunks);
                let free_bytes = core::mem::replace(&mut free.bytes, 0);
                drop(free);
                // Keep the counter conservative until deallocation completes:
                // another claimant must never observe capacity that is still
                // physically owned by these chunks.
                drop(drained);
                self.reserved_bytes.fetch_sub(free_bytes, Ordering::AcqRel);

                let reserved = self.reserved_bytes.load(Ordering::Acquire);
                let Some(next) = reserved.checked_add(want) else {
                    return Err(AllocError::ReservationOverflow {
                        site: site.name(),
                        base_bytes: reserved,
                        additional_bytes: want,
                    });
                };
                if let Some(limit) = self.config.limit_bytes
                    && next > limit
                {
                    return Err(AllocError::Exhausted {
                        site: site.name(),
                        requested_bytes: want,
                        reserved_bytes: reserved,
                        limit_bytes: limit,
                    });
                }
                continue;
            }
            if self
                .reserved_bytes
                .compare_exchange_weak(reserved, next, Ordering::AcqRel, Ordering::Acquire)
                .is_ok()
            {
                return Ok(());
            }
        }
    }

    /// Get a chunk of at least `min_bytes` (free list first, then the OS),
    /// enforcing the pool budget.
    fn acquire_chunk(
        &self,
        min_bytes: usize,
        want: usize,
        site: Site,
        lease: Option<&crate::OperationMemoryLease>,
    ) -> Result<ChunkAcquisition<'_>, AllocError> {
        let want = want.max(min_bytes);
        let mut retried_after_pool_pressure = false;
        loop {
            {
                let mut free = self
                    .free
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                if let Some((i, chunk_bytes)) = free
                    .chunks
                    .iter()
                    .enumerate()
                    .filter_map(|(i, chunk)| (chunk.len() >= min_bytes).then_some((i, chunk.len())))
                    .min_by_key(|(_, bytes)| *bytes)
                {
                    // A cached chunk can be larger than this arena's normal
                    // fresh chunk. Do not let that historical cache shape
                    // cause a false operation refusal when the smaller fresh
                    // request fits.
                    let prefer_fresh = lease.is_some_and(|lease| {
                        chunk_bytes > want
                            && !lease.can_reserve_now(chunk_bytes as u64)
                            && lease.can_reserve_now(want as u64)
                    });
                    if !prefer_fresh {
                        let charge = lease
                            .map(|lease| lease.reserve("arena-chunk", chunk_bytes as u64))
                            .transpose()
                            .map_err(|refusal| AllocError::LeaseExhausted {
                                site: site.name(),
                                requested_bytes: refusal.requested_bytes,
                                used_bytes: refusal.used_bytes,
                                limit_bytes: refusal.limit_bytes,
                            })?;
                        let chunk = free.chunks.swap_remove(i);
                        free.bytes -= chunk.len();
                        self.chunks_recycled.fetch_add(1, Ordering::Relaxed);
                        return Ok(ChunkAcquisition {
                            shared: self,
                            chunk: Some(chunk),
                            charge,
                        });
                    }
                }
            }
            let charge = lease
                .map(|lease| lease.reserve("arena-chunk", want as u64))
                .transpose()
                .map_err(|refusal| AllocError::LeaseExhausted {
                    site: site.name(),
                    requested_bytes: refusal.requested_bytes,
                    used_bytes: refusal.used_bytes,
                    limit_bytes: refusal.limit_bytes,
                })?;
            match self.claim_new_chunk_bytes(want, site) {
                Ok(()) => {}
                Err(error @ AllocError::Exhausted { .. }) if !retried_after_pool_pressure => {
                    // A concurrent arena can publish a suitable chunk after
                    // the pressure drain's final capacity check. Drop this
                    // tentative operation charge and recheck the free list
                    // exactly once before returning the refusal.
                    drop(charge);
                    let free = self
                        .free
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner);
                    let cached_bytes = free
                        .chunks
                        .iter()
                        .filter_map(|chunk| (chunk.len() >= min_bytes).then_some(chunk.len()))
                        .min();
                    drop(free);
                    let cached_is_admissible = cached_bytes.is_some_and(|chunk_bytes| {
                        lease.is_none_or(|lease| {
                            let prefer_fresh = chunk_bytes > want
                                && !lease.can_reserve_now(chunk_bytes as u64)
                                && lease.can_reserve_now(want as u64);
                            !prefer_fresh && lease.can_reserve_now(chunk_bytes as u64)
                        })
                    });
                    if cached_is_admissible {
                        retried_after_pool_pressure = true;
                        continue;
                    }
                    return Err(error);
                }
                Err(error) => return Err(error),
            }
            let align = self.hugepage.chunk_align(want);
            let Some(chunk) = Chunk::allocate(want, align) else {
                self.reserved_bytes.fetch_sub(want, Ordering::AcqRel);
                return Err(AllocError::OutOfMemory {
                    site: site.name(),
                    requested_bytes: want,
                });
            };
            debug_assert_eq!(chunk.len(), want);
            self.chunks_created.fetch_add(1, Ordering::Relaxed);
            return Ok(ChunkAcquisition {
                shared: self,
                chunk: Some(chunk),
                charge,
            });
        }
    }

    fn release_chunk(&self, chunk: Chunk) {
        let mut free = self
            .free
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(retained) = free
            .bytes
            .checked_add(chunk.len())
            .filter(|&bytes| bytes <= self.config.free_list_max_bytes)
        {
            free.bytes = retained;
            free.chunks.push(chunk);
        } else {
            let chunk_bytes = chunk.len();
            drop(chunk);
            self.reserved_bytes.fetch_sub(chunk_bytes, Ordering::AcqRel);
        }
    }

    /// Return chunks for reuse, respecting the free-list cap.
    fn release_chunks(&self, chunks: Vec<Chunk>) {
        let mut free = self
            .free
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        for chunk in chunks {
            if let Some(retained) = free
                .bytes
                .checked_add(chunk.len())
                .filter(|&bytes| bytes <= self.config.free_list_max_bytes)
            {
                free.bytes = retained;
                free.chunks.push(chunk);
            } else {
                let chunk_bytes = chunk.len();
                // Deallocate before advertising the capacity to concurrent
                // claimants. Temporary over-accounting is safe; temporary
                // under-accounting would violate the hard pool limit.
                drop(chunk);
                self.reserved_bytes.fetch_sub(chunk_bytes, Ordering::AcqRel);
            }
        }
    }
}

/// Factory and accounting hub for scope arenas. Cheap to clone (`Arc`);
/// `Send + Sync` so worker threads can draw per-tile arenas from one pool.
#[derive(Clone)]
pub struct ArenaPool {
    shared: Arc<PoolShared>,
}

impl ArenaPool {
    /// Build a pool, probing and RECORDING the hugepage decision once.
    #[must_use]
    pub fn new(config: ArenaConfig) -> Self {
        let mut config = config;
        config.chunk_bytes = config.chunk_bytes.max(4096);
        config.max_chunk_bytes = config.max_chunk_bytes.max(config.chunk_bytes);
        let hugepage = HugepageDecision::probe(config.hugepage, config.chunk_bytes);
        ArenaPool {
            shared: Arc::new(PoolShared {
                hugepage,
                free: Mutex::new(FreeList {
                    chunks: Vec::new(),
                    bytes: 0,
                }),
                reserved_bytes: AtomicUsize::new(0),
                arenas_live: AtomicUsize::new(0),
                chunks_created: AtomicU64::new(0),
                chunks_recycled: AtomicU64::new(0),
                sites: Mutex::new(BTreeMap::new()),
                config,
            }),
        }
    }

    /// Create a fresh arena for one unit of scoped work. Prefer
    /// [`ArenaPool::scope`] unless the arena's lifetime is managed by an
    /// executor (fs-exec's use case).
    #[must_use]
    pub fn arena(&self) -> Arena {
        self.build_arena(None)
    }

    /// Create a fresh arena whose chunks charge `lease` while held (bead
    /// wf9.16). The pool's own `limit_bytes` remains an independent gate:
    /// both must admit a chunk, and a refusal names whichever refused.
    #[must_use]
    pub fn arena_leased(&self, lease: &crate::OperationMemoryLease) -> Arena {
        self.build_arena(Some(lease.clone()))
    }

    fn build_arena(&self, lease: Option<crate::OperationMemoryLease>) -> Arena {
        self.shared.arenas_live.fetch_add(1, Ordering::AcqRel);
        Arena {
            shared: Arc::clone(&self.shared),
            raw: RawArena::new(),
            next_chunk_bytes: Cell::new(self.shared.config.chunk_bytes),
            allocated_bytes: Cell::new(0),
            allocation_count: Cell::new(0),
            sites: RefCell::new(Vec::new()),
            lease,
            leased_bytes: Cell::new(0),
        }
    }

    /// Bytes a fresh arena must reserve for its first slice allocation.
    ///
    /// This is an allocation-free preflight using the pool's normalized first
    /// chunk size and the exact alignment slack used by [`Arena::grow`]. Empty
    /// and zero-sized slices require no chunk. A free-list hit may hand the
    /// arena a larger already-reserved chunk, but cannot increase the pool's
    /// OS reservation.
    ///
    /// # Errors
    /// [`AllocError::LayoutOverflow`] or [`AllocError::ReservationOverflow`]
    /// when the request cannot be represented.
    pub fn reservation_bytes_for_slice<T>(
        &self,
        site: Site,
        len: usize,
    ) -> Result<usize, AllocError> {
        let bytes = slice_bytes::<T>(site, len)?;
        if bytes == 0 {
            return Ok(0);
        }
        Ok(self
            .shared
            .config
            .chunk_bytes
            .max(reservation_min_bytes(bytes, align_of::<T>(), site)?))
    }

    /// Run `f` with a scope arena; the arena and every allocation in it are
    /// reclaimed when `f` returns (or unwinds). The higher-ranked closure
    /// bound makes escaping an allocation a COMPILE ERROR:
    ///
    /// ```compile_fail
    /// let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    /// let escaped = pool.scope(|a| {
    ///     a.alloc(fs_alloc::Site::named("t/escape"), 7u64).unwrap()
    /// }); // ERROR: allocation cannot outlive the scope
    /// ```
    ///
    /// Smuggling through an outer binding fails the same way:
    ///
    /// ```compile_fail
    /// let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    /// let mut smuggled: Option<&mut u64> = None;
    /// pool.scope(|a| {
    ///     smuggled = Some(a.alloc(fs_alloc::Site::named("t/smuggle"), 7u64).unwrap());
    /// }); // ERROR: borrow escapes the higher-ranked scope lifetime
    /// ```
    pub fn scope<R>(&self, f: impl for<'a> FnOnce(&'a Arena) -> R) -> R {
        let arena = self.arena();
        f(&arena)
    }

    /// [`ArenaPool::scope`] with the scope's chunks charged to `lease`
    /// while held (bead wf9.16). Same escape discipline.
    pub fn scope_leased<R>(
        &self,
        lease: &crate::OperationMemoryLease,
        f: impl for<'a> FnOnce(&'a Arena) -> R,
    ) -> R {
        let arena = self.arena_leased(lease);
        f(&arena)
    }

    /// Snapshot of pool-level accounting.
    #[must_use]
    pub fn stats(&self) -> PoolStats {
        let free_bytes = self
            .shared
            .free
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .bytes;
        PoolStats {
            arenas_live: self.shared.arenas_live.load(Ordering::Acquire),
            reserved_bytes: self.shared.reserved_bytes.load(Ordering::Acquire),
            free_bytes,
            chunks_created: self.shared.chunks_created.load(Ordering::Relaxed),
            chunks_recycled: self.shared.chunks_recycled.load(Ordering::Relaxed),
            hugepage: self.shared.hugepage.clone(),
        }
    }

    /// The hugepage decision recorded at pool construction.
    #[must_use]
    pub fn hugepage_decision(&self) -> &HugepageDecision {
        &self.shared.hugepage
    }

    /// Deterministic per-site accounting report, accumulated from every
    /// arena this pool has dropped so far (sorted by site name — diffable
    /// between runs).
    #[must_use]
    pub fn site_report(&self) -> SiteReport {
        let sites = self
            .shared
            .sites
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        SiteReport {
            sites: sites.iter().map(|(k, v)| (*k, *v)).collect(),
        }
    }
}

impl fmt::Debug for ArenaPool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ArenaPool")
            .field("stats", &self.stats().to_json())
            .finish_non_exhaustive()
    }
}

/// A bump arena for one unit of scoped work. `Send` (movable to a worker)
/// but deliberately `!Sync` — tiles get their own arenas instead of sharing
/// one:
///
/// ```compile_fail
/// fn assert_sync<T: Sync>() {}
/// assert_sync::<fs_alloc::Arena>(); // ERROR: Arena is intentionally !Sync
/// ```
pub struct Arena {
    shared: Arc<PoolShared>,
    raw: RawArena,
    next_chunk_bytes: Cell<usize>,
    allocated_bytes: Cell<u64>,
    allocation_count: Cell<u64>,
    sites: RefCell<Vec<(&'static str, SiteStats)>>,
    /// Operation lease this arena's chunks charge while held (bead wf9.16).
    /// `None` = legacy unleased arena. The charge covers the hold interval
    /// only: chunks returned to the pool free list on drop stop being the
    /// operation's live set, so a recycled chunk is never double-charged.
    lease: Option<crate::OperationMemoryLease>,
    leased_bytes: Cell<u64>,
}

impl Arena {
    /// Allocate one value. The reference borrows the arena, so it cannot
    /// outlive it:
    ///
    /// ```compile_fail
    /// let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    /// let r;
    /// {
    ///     let a = pool.arena();
    ///     r = a.alloc(fs_alloc::Site::named("t/outlive"), 7u64).unwrap();
    /// } // ERROR: `a` dropped while `r` still borrows it
    /// let _ = *r;
    /// ```
    ///
    /// Types with destructors are rejected at compile time (bump arenas
    /// never run `Drop`):
    ///
    /// ```compile_fail
    /// let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    /// pool.scope(|a| {
    ///     let _ = a.alloc(fs_alloc::Site::named("t/drop"), String::from("no"));
    /// }); // ERROR: post-monomorphization assert — String needs Drop
    /// ```
    ///
    /// # Errors
    /// [`AllocError::Exhausted`] / [`AllocError::OutOfMemory`] when the pool
    /// cannot supply a chunk.
    // &self -> &mut is the arena shape: returns are disjoint by the bump
    // discipline (capsule SAFETY.md), so the exclusive borrows never alias.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc<T>(&self, site: Site, value: T) -> Result<&mut T, AllocError> {
        match self.raw.try_place(value) {
            Ok(r) => {
                self.note(site, size_of::<T>(), align_of::<T>());
                Ok(r)
            }
            Err(value) => {
                self.grow(size_of::<T>(), align_of::<T>(), site)?;
                // Defensive fallback: grow() sizes the window for the request.
                let Ok(r) = self.raw.try_place(value) else {
                    return Err(self.exhausted_defensively(site, size_of::<T>()));
                };
                self.note(site, size_of::<T>(), align_of::<T>());
                Ok(r)
            }
        }
    }

    /// Allocate a slice of `len` copies of `fill` (128-byte aligned).
    ///
    /// # Errors
    /// [`AllocError::LayoutOverflow`] on `len * size_of::<T>()` overflow,
    /// otherwise as [`Arena::alloc`].
    // &self -> &mut: see Arena::alloc.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc_slice_fill<T: Copy>(
        &self,
        site: Site,
        len: usize,
        fill: T,
    ) -> Result<&mut [T], AllocError> {
        let bytes = slice_bytes::<T>(site, len)?;
        if let Ok(s) = self.raw.try_place_slice_fill(len, fill) {
            self.note(site, bytes, align_of::<T>());
            return Ok(s);
        }
        self.grow(bytes, align_of::<T>(), site)?;
        // Defensive fallback: grow() sizes the window for the request.
        let Ok(s) = self.raw.try_place_slice_fill(len, fill) else {
            return Err(self.exhausted_defensively(site, bytes));
        };
        self.note(site, bytes, align_of::<T>());
        Ok(s)
    }

    /// Allocate a slice built element-by-element from `f(i)` (128-byte
    /// aligned). If `f` panics midway, the arena stays usable and the
    /// reserved bytes are reclaimed with the scope (no leak, no torn state).
    ///
    /// # Errors
    /// As [`Arena::alloc_slice_fill`].
    // &self -> &mut: see Arena::alloc.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc_slice_with<T>(
        &self,
        site: Site,
        len: usize,
        mut f: impl FnMut(usize) -> T,
    ) -> Result<&mut [T], AllocError> {
        let bytes = slice_bytes::<T>(site, len)?;
        if let Ok(s) = self.raw.try_place_slice_with(len, &mut f) {
            self.note(site, bytes, align_of::<T>());
            return Ok(s);
        }
        self.grow(bytes, align_of::<T>(), site)?;
        // Defensive fallback: grow() sizes the window for the request.
        let Ok(s) = self.raw.try_place_slice_with(len, &mut f) else {
            return Err(self.exhausted_defensively(site, bytes));
        };
        self.note(site, bytes, align_of::<T>());
        Ok(s)
    }

    /// Run `f` with a CHILD scope arena drawn from the same pool; the child
    /// and all its allocations are reclaimed when `f` returns. (Chunks
    /// recycle through the pool, so tight scope loops do not thrash the OS
    /// allocator.)
    pub fn scope<R>(&self, f: impl for<'c> FnOnce(&'c Arena) -> R) -> R {
        self.shared.arenas_live.fetch_add(1, Ordering::AcqRel);
        // A child scope inherits the parent's operation lease (bead
        // wf9.16): sub-scopes of a leased tile stay inside the same
        // operation live set.
        let child = Arena {
            shared: Arc::clone(&self.shared),
            raw: RawArena::new(),
            next_chunk_bytes: Cell::new(self.shared.config.chunk_bytes),
            allocated_bytes: Cell::new(0),
            allocation_count: Cell::new(0),
            sites: RefCell::new(Vec::new()),
            lease: self.lease.clone(),
            leased_bytes: Cell::new(0),
        };
        f(&child)
    }

    /// Snapshot of this arena's accounting.
    #[must_use]
    pub fn stats(&self) -> ArenaStats {
        ArenaStats {
            allocated_bytes: self.allocated_bytes.get(),
            allocation_count: self.allocation_count.get(),
            chunk_count: self.raw.chunk_count(),
        }
    }

    /// Verification hook: base address of the current chunk (used by the
    /// conformance suite to check the recorded hugepage alignment; addresses
    /// are never part of deterministic reports).
    #[must_use]
    pub fn current_chunk_base(&self) -> Option<usize> {
        self.raw.last_chunk_base()
    }

    /// Install a fresh chunk sized for a `bytes`/`align` request.
    fn grow(&self, bytes: usize, align: usize, site: Site) -> Result<(), AllocError> {
        let min_bytes = reservation_min_bytes(bytes, align, site)?;
        let acquisition = self.shared.acquire_chunk(
            min_bytes,
            self.next_chunk_bytes.get(),
            site,
            self.lease.as_ref(),
        )?;
        // Admission precedes recycled-chunk removal and fresh allocation, so
        // the operation lease covers the complete hold interval. Pool/OS
        // failures roll that charge back inside acquire_chunk.
        let next_chunk_bytes = acquisition
            .len()
            .checked_mul(2)
            .unwrap_or(self.shared.config.max_chunk_bytes)
            .min(self.shared.config.max_chunk_bytes);
        acquisition.install_into(self, next_chunk_bytes);
        Ok(())
    }

    /// Record an allocation under `site`. Accounted bytes are the payload
    /// padded to `max(align, 128)` — the window consumption for the common
    /// `align <= 128` case (higher alignments may consume extra window
    /// padding, visible only as reserved-vs-allocated slack).
    fn note(&self, site: Site, bytes: usize, align: usize) {
        let padded = padded_bytes(bytes, align);
        self.allocated_bytes
            .set(self.allocated_bytes.get().saturating_add(padded));
        self.allocation_count
            .set(self.allocation_count.get().saturating_add(1));
        let mut sites = self.sites.borrow_mut();
        if let Some((_, stats)) = sites.iter_mut().find(|(name, _)| *name == site.name()) {
            stats.bytes = stats.bytes.saturating_add(padded);
            stats.allocations = stats.allocations.saturating_add(1);
        } else {
            sites.push((
                site.name(),
                SiteStats {
                    bytes: padded,
                    allocations: 1,
                },
            ));
        }
    }

    fn exhausted_defensively(&self, site: Site, bytes: usize) -> AllocError {
        debug_assert!(false, "grow() must size the window for the request");
        AllocError::Exhausted {
            site: site.name(),
            requested_bytes: bytes,
            reserved_bytes: self.shared.reserved_bytes.load(Ordering::Acquire),
            limit_bytes: self.shared.config.limit_bytes.unwrap_or(usize::MAX),
        }
    }
}

impl fmt::Debug for Arena {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Arena")
            .field("stats", &self.stats().to_json())
            .finish_non_exhaustive()
    }
}

impl Drop for Arena {
    fn drop(&mut self) {
        let chunks = self.raw.take_chunks();
        // Chunks returned to the pool free list leave this operation's live
        // set. Drop the operation charge BEFORE publishing those chunks to
        // the shared free list, so a concurrent recycler can never observe a
        // free chunk that is still charged to its previous operation.
        if let Some(lease) = &self.lease {
            let held = self.leased_bytes.get();
            if held > 0 {
                let _released = lease.release_raw(held);
            }
        }
        self.shared.release_chunks(chunks);
        let local = core::mem::take(self.sites.get_mut());
        if !local.is_empty() {
            let mut sites = self
                .shared
                .sites
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            for (name, stats) in local {
                let entry = sites.entry(name).or_default();
                entry.bytes = entry.bytes.saturating_add(stats.bytes);
                entry.allocations = entry.allocations.saturating_add(stats.allocations);
            }
        }
        self.shared.arenas_live.fetch_sub(1, Ordering::AcqRel);
    }
}

/// Slice payload bytes, or a structured overflow error.
fn slice_bytes<T>(site: Site, len: usize) -> Result<usize, AllocError> {
    size_of::<T>()
        .checked_mul(len)
        .ok_or(AllocError::LayoutOverflow {
            site: site.name(),
            len,
            elem_bytes: size_of::<T>(),
        })
}

fn reservation_min_bytes(bytes: usize, align: usize, site: Site) -> Result<usize, AllocError> {
    if bytes == 0 {
        return Ok(0);
    }
    // Slack covers worst-case in-window padding when T's alignment exceeds
    // the chunk's unconditional 128-byte base alignment.
    let slack = if align > ALLOC_ALIGN { align } else { 0 };
    bytes
        .checked_add(slack)
        .ok_or(AllocError::ReservationOverflow {
            site: site.name(),
            base_bytes: bytes,
            additional_bytes: slack,
        })
}

/// Payload bytes padded to the accounting granularity.
fn padded_bytes(bytes: usize, align: usize) -> u64 {
    if bytes == 0 {
        return 0;
    }
    let pad = align.max(ALLOC_ALIGN) as u64;
    (bytes as u64).div_ceil(pad) * pad
}

/// Per-arena accounting snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArenaStats {
    /// Cumulative padded payload bytes allocated by this arena.
    /// Saturates at `u64::MAX` rather than wrapping.
    pub allocated_bytes: u64,
    /// Cumulative allocation count.
    /// Saturates at `u64::MAX` rather than wrapping.
    pub allocation_count: u64,
    /// Chunks currently owned by this arena.
    pub chunk_count: usize,
}

impl ArenaStats {
    /// Canonical JSON object (deterministic field order).
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"allocated_bytes\":{},\"allocation_count\":{},\"chunk_count\":{}}}",
            self.allocated_bytes, self.allocation_count, self.chunk_count
        )
    }
}

/// Pool-level accounting snapshot. `quiescent()` is the G4 leak oracle.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PoolStats {
    /// Arenas currently alive.
    pub arenas_live: usize,
    /// Total OS-reserved bytes (in-use + free-listed).
    pub reserved_bytes: usize,
    /// Bytes parked in the chunk free list.
    pub free_bytes: usize,
    /// Chunks ever obtained from the OS.
    pub chunks_created: u64,
    /// Chunk acquisitions served from the free list.
    pub chunks_recycled: u64,
    /// The recorded hugepage decision.
    pub hugepage: HugepageDecision,
}

impl PoolStats {
    /// True when no arena is alive and every reserved byte is parked in the
    /// free list — i.e. nothing leaked. The G4 storm asserts this after
    /// 10^6 random cancellations.
    #[must_use]
    pub fn quiescent(&self) -> bool {
        self.arenas_live == 0 && self.reserved_bytes == self.free_bytes
    }

    /// Canonical JSON object (deterministic field order; no addresses, no
    /// clocks — safe for G5 comparisons).
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"arenas_live\":{},\"reserved_bytes\":{},\"free_bytes\":{},\
             \"chunks_created\":{},\"chunks_recycled\":{},\"hugepage\":{}}}",
            self.arenas_live,
            self.reserved_bytes,
            self.free_bytes,
            self.chunks_created,
            self.chunks_recycled,
            self.hugepage.to_json()
        )
    }

    /// Package as an `fs-obs` event payload for the ledger pipeline.
    #[must_use]
    pub fn to_event_kind(&self) -> fs_obs::EventKind {
        fs_obs::EventKind::Custom {
            name: "fs-alloc-pool-stats".to_string(),
            json: self.to_json(),
        }
    }
}

/// Deterministic per-site accounting report (sorted by site name).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteReport {
    /// `(site, stats)` pairs, sorted by site name.
    pub sites: Vec<(&'static str, SiteStats)>,
}

impl SiteReport {
    /// Canonical JSON object (deterministic order — DIFFABLE between runs,
    /// which is the point of site tracking).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::from("{\"sites\":[");
        for (i, (name, stats)) in self.sites.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let escaped_name = crate::lease::json_escape(name);
            let _ = write!(
                s,
                "{{\"site\":\"{escaped_name}\",\"bytes\":{},\"allocations\":{}}}",
                stats.bytes, stats.allocations
            );
        }
        s.push_str("]}");
        s
    }

    /// Package as an `fs-obs` event payload for the ledger pipeline.
    #[must_use]
    pub fn to_event_kind(&self) -> fs_obs::EventKind {
        fs_obs::EventKind::Custom {
            name: "fs-alloc-site-report".to_string(),
            json: self.to_json(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_pool(limit: Option<usize>) -> ArenaPool {
        ArenaPool::new(ArenaConfig {
            chunk_bytes: 4096,
            max_chunk_bytes: 16 * 4096,
            limit_bytes: limit,
            free_list_max_bytes: 1 << 20,
            hugepage: HugepagePolicy::Never,
        })
    }

    #[test]
    fn first_slice_reservation_uses_normalized_chunk_and_checked_layout() {
        #[derive(Clone, Copy)]
        #[repr(align(256))]
        struct Wide(u8);

        let default_pool = ArenaPool::new(ArenaConfig::default());
        assert_eq!(
            default_pool
                .reservation_bytes_for_slice::<f64>(Site::named("t/plan-default"), 8 * 256)
                .expect("A micro-panel plan"),
            1 << 20,
            "a fresh default arena reserves its 1 MiB first chunk"
        );

        let pool = small_pool(None);
        assert_eq!(
            pool.reservation_bytes_for_slice::<u8>(Site::named("t/plan-small"), 1024)
                .expect("small plan"),
            4096
        );
        assert_eq!(
            pool.reservation_bytes_for_slice::<u8>(Site::named("t/plan-large"), 8192)
                .expect("large plan"),
            8192
        );
        assert_eq!(
            pool.reservation_bytes_for_slice::<u64>(Site::named("t/plan-empty"), 0)
                .expect("empty plan"),
            0
        );

        let error = pool
            .reservation_bytes_for_slice::<Wide>(
                Site::named("t/plan-overflow"),
                usize::MAX / core::mem::size_of::<Wide>(),
            )
            .expect_err("alignment slack must overflow structurally");
        assert!(
            matches!(error, AllocError::ReservationOverflow { .. }),
            "{error:?}"
        );
        let _ = Wide(0).0;
    }

    #[test]
    fn concurrent_chunk_claims_cannot_cross_the_pool_limit() {
        const THREADS: usize = 8;
        let pool = std::sync::Arc::new(small_pool(Some(4096)));
        let start = std::sync::Arc::new(std::sync::Barrier::new(THREADS));
        let hold = std::sync::Arc::new(std::sync::Barrier::new(THREADS));
        let successes = std::sync::atomic::AtomicUsize::new(0);

        std::thread::scope(|scope| {
            for _ in 0..THREADS {
                let pool = std::sync::Arc::clone(&pool);
                let start = std::sync::Arc::clone(&start);
                let hold = std::sync::Arc::clone(&hold);
                let successes = &successes;
                scope.spawn(move || {
                    start.wait();
                    pool.scope(|arena| {
                        if arena.alloc(Site::named("t/concurrent-limit"), 1_u8).is_ok() {
                            successes.fetch_add(1, Ordering::AcqRel);
                        }
                        // Keep the successful arena live until every competing
                        // reservation has reached the hard-limit check.
                        hold.wait();
                    });
                });
            }
        });

        assert_eq!(successes.load(Ordering::Acquire), 1);
        let stats = pool.stats();
        assert!(stats.reserved_bytes <= 4096, "{}", stats.to_json());
        assert!(stats.quiescent(), "{}", stats.to_json());
    }

    #[test]
    fn happy_path_alloc_and_reclaim() {
        let pool = small_pool(None);
        let sum = pool.scope(|a| {
            let x = a.alloc(Site::named("t/x"), 40u64).expect("alloc");
            let ys = a
                .alloc_slice_fill(Site::named("t/ys"), 100, 0.5f64)
                .expect("slice");
            let zs = a
                .alloc_slice_with(Site::named("t/zs"), 4, |i| i as u64)
                .expect("slice_with");
            assert_eq!(a.stats().allocation_count, 3);
            *x + ys.iter().sum::<f64>() as u64 + zs.iter().sum::<u64>()
        });
        assert_eq!(sum, 40 + 50 + 6);
        assert!(pool.stats().quiescent(), "{}", pool.stats().to_json());
    }

    #[test]
    fn growth_crosses_chunks_and_recycles() {
        let pool = small_pool(None);
        for _ in 0..3 {
            pool.scope(|a| {
                for i in 0..64 {
                    let s = a
                        .alloc_slice_fill(Site::named("t/grow"), 256, i as u8)
                        .expect("fits");
                    assert_eq!(s[0], i as u8);
                }
                assert!(a.stats().chunk_count >= 2, "growth must add chunks");
            });
        }
        let stats = pool.stats();
        assert!(stats.quiescent());
        assert!(
            stats.chunks_recycled > 0,
            "later scopes must reuse chunks: {}",
            stats.to_json()
        );
    }

    #[test]
    fn budget_exhaustion_is_structured_not_fatal() {
        let pool = small_pool(Some(8192));
        let err = pool.scope(|a| {
            a.alloc_slice_fill(Site::named("t/budget"), 1 << 20, 0u8)
                .expect_err("must exceed the 8 KiB limit")
        });
        match &err {
            AllocError::Exhausted {
                site, limit_bytes, ..
            } => {
                assert_eq!(*site, "t/budget");
                assert_eq!(*limit_bytes, 8192);
            }
            other => panic!("wrong error variant: {other:?}"),
        }
        let msg = err.to_string();
        assert!(
            msg.contains("limit_bytes") || msg.contains("limit"),
            "{msg}"
        );
        // The pool remains fully usable after refusal.
        pool.scope(|a| {
            a.alloc(Site::named("t/after"), 1u8)
                .expect("small alloc still fine");
        });
        assert!(pool.stats().quiescent());
    }

    #[test]
    fn leased_recycled_chunk_is_admitted_before_removal() {
        let pool = small_pool(None);
        pool.scope(|arena| {
            arena
                .alloc(Site::named("t/seed-recycled"), 1_u8)
                .expect("seed free list");
        });
        let seeded = pool.stats();
        assert_eq!(seeded.free_bytes, 4096);

        let refusing = crate::OperationMemoryLease::bounded(0);
        let error = pool.scope_leased(&refusing, |arena| {
            arena
                .alloc(Site::named("t/refuse-recycled"), 1_u8)
                .expect_err("operation lease refuses before chunk removal")
        });
        assert!(matches!(error, AllocError::LeaseExhausted { .. }));
        assert_eq!(refusing.receipt().used_bytes, 0);
        let after_refusal = pool.stats();
        assert_eq!(after_refusal.free_bytes, seeded.free_bytes);
        assert_eq!(after_refusal.chunks_recycled, seeded.chunks_recycled);

        let admitted = crate::OperationMemoryLease::unbounded();
        pool.scope_leased(&admitted, |arena| {
            arena
                .alloc(Site::named("t/retry-recycled"), 1_u8)
                .expect("unchanged free-list chunk remains reusable");
        });
        let after_retry = pool.stats();
        assert_eq!(after_retry.chunks_created, seeded.chunks_created);
        assert_eq!(after_retry.chunks_recycled, seeded.chunks_recycled + 1);
        assert_eq!(admitted.receipt().used_bytes, 0);
    }

    #[test]
    fn oversized_cached_chunk_does_not_force_a_false_lease_refusal() {
        let pool = small_pool(None);
        pool.scope(|arena| {
            arena
                .alloc_slice_fill(Site::named("t/seed-oversized"), 8192, 0_u8)
                .expect("seed an oversized cached chunk");
        });
        let seeded = pool.stats();
        assert_eq!(seeded.free_bytes, 8192);

        let lease = crate::OperationMemoryLease::bounded(4096);
        pool.scope_leased(&lease, |arena| {
            arena
                .alloc(Site::named("t/fresh-below-oversized-cache"), 1_u8)
                .expect("the 4 KiB fresh chunk fits the lease");
        });
        let receipt = lease.receipt();
        assert_eq!(receipt.requested_bytes, 4096);
        assert_eq!(receipt.used_bytes, 0);
        assert_eq!(receipt.refusals, 0);
        let after_fresh = pool.stats();
        assert_eq!(after_fresh.chunks_created, seeded.chunks_created + 1);
        assert_eq!(after_fresh.chunks_recycled, seeded.chunks_recycled);
        assert_eq!(after_fresh.free_bytes, 8192 + 4096);

        let recycled = crate::OperationMemoryLease::bounded(4096);
        pool.scope_leased(&recycled, |arena| {
            arena
                .alloc(Site::named("t/recycle-smallest"), 1_u8)
                .expect("the smallest sufficient cached chunk is selected");
        });
        assert_eq!(recycled.receipt().requested_bytes, 4096);
        assert_eq!(
            pool.stats().chunks_recycled,
            after_fresh.chunks_recycled + 1
        );
    }

    #[test]
    fn pending_chunk_unwind_rolls_back_both_ledgers() {
        let pool = small_pool(None);
        let lease = crate::OperationMemoryLease::bounded(4096);
        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let pending = pool
                .shared
                .acquire_chunk(1, 4096, Site::named("t/pending-unwind"), Some(&lease))
                .expect("both gates admit the pending chunk");
            assert_eq!(pending.len(), 4096);
            assert_eq!(lease.receipt().used_bytes, 4096);
            panic!("unwind before arena installation");
        }));
        assert!(panicked.is_err());
        let receipt = lease.receipt();
        assert_eq!(receipt.used_bytes, 0);
        assert_eq!(receipt.release_invariant_violations, 0);
        let stats = pool.stats();
        assert_eq!(stats.free_bytes, 4096);
        assert!(stats.quiescent(), "{}", stats.to_json());
    }

    #[test]
    fn poisoned_accounting_locks_do_not_break_acquire_recycle_or_cleanup() {
        let pool = small_pool(None);
        let poisoned = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _free = pool.shared.free.lock().expect("initial lock");
            panic!("poison the free-list mutex without changing its state");
        }));
        assert!(poisoned.is_err());
        assert!(pool.shared.free.is_poisoned());
        let poisoned_sites = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _sites = pool.shared.sites.lock().expect("initial site lock");
            panic!("poison the site-table mutex without changing its state");
        }));
        assert!(poisoned_sites.is_err());
        assert!(pool.shared.sites.is_poisoned());

        let first = crate::OperationMemoryLease::bounded(4096);
        pool.scope_leased(&first, |arena| {
            arena
                .alloc(Site::named("t/poison-fresh"), 1_u8)
                .expect("fresh acquisition recovers the poisoned lock");
        });
        assert_eq!(first.receipt().used_bytes, 0);
        assert!(pool.stats().quiescent());

        let recycled_before = pool.stats().chunks_recycled;
        let second = crate::OperationMemoryLease::bounded(4096);
        pool.scope_leased(&second, |arena| {
            arena
                .alloc(Site::named("t/poison-recycled"), 1_u8)
                .expect("recycled acquisition recovers the poisoned lock");
        });
        assert_eq!(pool.stats().chunks_recycled, recycled_before + 1);
        assert_eq!(second.receipt().used_bytes, 0);
        assert!(pool.stats().quiescent());
        let report = pool.site_report();
        assert_eq!(report.sites.len(), 2);
        assert_eq!(report.sites[0].0, "t/poison-fresh");
        assert_eq!(report.sites[1].0, "t/poison-recycled");
    }

    #[test]
    fn pool_refusal_rolls_back_prior_operation_admission() {
        let pool = small_pool(Some(0));
        let lease = crate::OperationMemoryLease::unbounded();
        let error = pool.scope_leased(&lease, |arena| {
            arena
                .alloc(Site::named("t/pool-refusal-rollback"), 1_u8)
                .expect_err("zero-byte pool gate refuses")
        });
        assert!(matches!(error, AllocError::Exhausted { .. }));
        let receipt = lease.receipt();
        assert_eq!(receipt.requested_bytes, 4096);
        assert_eq!(receipt.used_bytes, 0, "failed pool gate must roll back");
        assert_eq!(receipt.refusals, 0, "the operation gate admitted");
        assert!(pool.stats().quiescent());
    }

    #[test]
    fn layout_overflow_is_reported() {
        let pool = small_pool(None);
        pool.scope(|a| {
            let err = a
                .alloc_slice_fill(Site::named("t/overflow"), usize::MAX / 4, 0u64)
                .expect_err("layout must overflow");
            assert!(matches!(err, AllocError::LayoutOverflow { .. }), "{err:?}");
        });
    }

    #[test]
    fn zst_and_empty_allocations_cost_nothing() {
        let pool = small_pool(None);
        pool.scope(|a| {
            let unit = a.alloc(Site::named("t/zst"), ()).expect("zst");
            assert_eq!(*unit, ());
            let empty = a
                .alloc_slice_fill(Site::named("t/empty"), 0, 0u8)
                .expect("empty");
            assert!(empty.is_empty());
            assert_eq!(a.stats().allocated_bytes, 0);
            assert_eq!(a.stats().allocation_count, 2);
        });
        assert!(pool.stats().quiescent());
    }

    #[test]
    fn over_128_alignment_is_honored() {
        #[repr(align(512))]
        #[derive(Clone, Copy)]
        struct Big([u8; 512]);
        let pool = small_pool(None);
        pool.scope(|a| {
            let b = a
                .alloc(Site::named("t/align512"), Big([7; 512]))
                .expect("alloc");
            assert_eq!(core::ptr::from_mut(b) as usize % 512, 0);
            assert_eq!(b.0[511], 7);
        });
        assert!(pool.stats().quiescent());
    }

    #[test]
    fn panic_mid_fill_leaves_arena_usable_and_leak_free() {
        let pool = small_pool(None);
        pool.scope(|a| {
            let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let _ = a.alloc_slice_with(Site::named("t/panic"), 16, |i| {
                    assert!(i < 8, "boom at element 8");
                    i as u64
                });
            }))
            .is_err();
            assert!(panicked, "the fill closure must panic");
            // Arena still works after the unwind.
            let x = a.alloc(Site::named("t/after-panic"), 3u32).expect("alloc");
            assert_eq!(*x, 3);
        });
        assert!(pool.stats().quiescent(), "{}", pool.stats().to_json());
    }

    #[test]
    fn nested_scopes_account_independently_and_merge_sites() {
        let pool = small_pool(None);
        pool.scope(|outer| {
            outer
                .alloc_slice_fill(Site::named("t/outer"), 64, 1u8)
                .expect("outer");
            outer.scope(|inner| {
                inner
                    .alloc_slice_fill(Site::named("t/inner"), 64, 2u8)
                    .expect("inner");
                assert_eq!(inner.stats().allocation_count, 1);
            });
            assert_eq!(outer.stats().allocation_count, 1);
        });
        let report = pool.site_report();
        assert_eq!(report.sites.len(), 2);
        assert_eq!(report.sites[0].0, "t/inner", "sorted by name");
        assert_eq!(report.sites[1].0, "t/outer");
        assert!(pool.stats().quiescent());
    }

    #[test]
    fn site_report_json_is_deterministic_and_diffable() {
        let run = || {
            let pool = small_pool(None);
            pool.scope(|a| {
                a.alloc_slice_fill(Site::named("t/b"), 10, 0u8).expect("b");
                a.alloc(Site::named("t/a"), 1u64).expect("a");
                a.alloc_slice_fill(Site::named("t/b"), 10, 0u8).expect("b2");
            });
            pool.site_report().to_json()
        };
        let (r1, r2) = (run(), run());
        assert_eq!(r1, r2);
        assert_eq!(
            r1,
            "{\"sites\":[{\"site\":\"t/a\",\"bytes\":128,\"allocations\":1},\
             {\"site\":\"t/b\",\"bytes\":256,\"allocations\":2}]}"
        );
    }

    #[test]
    fn site_counters_saturate_instead_of_wrapping() {
        let pool = small_pool(None);
        let site = Site::named("t/saturating-site");
        let arena = pool.arena();
        arena.allocated_bytes.set(u64::MAX - 64);
        arena.allocation_count.set(u64::MAX);
        arena.sites.borrow_mut().push((
            site.name(),
            SiteStats {
                bytes: u64::MAX - 64,
                allocations: u64::MAX,
            },
        ));
        pool.shared.sites.lock().expect("site table").insert(
            site.name(),
            SiteStats {
                bytes: u64::MAX - 64,
                allocations: u64::MAX,
            },
        );

        arena.note(site, 1, 1);
        assert_eq!(arena.allocated_bytes.get(), u64::MAX);
        assert_eq!(arena.allocation_count.get(), u64::MAX);
        assert_eq!(arena.sites.borrow()[0].1.bytes, u64::MAX);
        assert_eq!(arena.sites.borrow()[0].1.allocations, u64::MAX);
        drop(arena);

        let report = pool.site_report();
        assert_eq!(report.sites[0].1.bytes, u64::MAX);
        assert_eq!(report.sites[0].1.allocations, u64::MAX);
    }

    #[test]
    fn site_report_json_escapes_hostile_static_names() {
        let pool = small_pool(None);
        let hostile = Site::named("t/quote\"slash\\line\ncontrol\u{001f}");
        pool.scope(|arena| {
            arena.alloc(hostile, 1_u8).expect("hostile-name allocation");
        });
        let json = pool.site_report().to_json();
        assert!(!json.contains('\n'));
        assert!(json.contains("t/quote\\\"slash\\\\line\\ncontrol\\u001f"));

        let mut emitter = fs_obs::Emitter::new("fs-alloc-test", "hostile-site");
        let line = emitter
            .emit(
                fs_obs::Severity::Info,
                pool.site_report().to_event_kind(),
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).unwrap_or_else(|error| panic!("{line}: {error}"));
    }

    #[test]
    fn stats_events_validate_against_the_obs_schema() {
        let pool = small_pool(None);
        pool.scope(|a| {
            a.alloc(Site::named("t/evt"), 1u8).expect("alloc");
        });
        let mut em = fs_obs::Emitter::new("fs-alloc-test", "arena-tests");
        for kind in [
            pool.stats().to_event_kind(),
            pool.site_report().to_event_kind(),
        ] {
            let line = em.emit(fs_obs::Severity::Info, kind, None).to_jsonl();
            fs_obs::validate_line(&line).unwrap_or_else(|e| panic!("{line}: {e}"));
        }
    }
}
