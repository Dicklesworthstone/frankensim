//! The bump-pointer core: fs-alloc's primary unsafe capsule.
//!
//! Chunk allocation/deallocation against the global allocator and placement
//! of values into bump-allocated memory live here. The separate reclaim-poison
//! byte-access capsule is in `poison`; the safe facade (`crate::arena`)
//! contains zero `unsafe` tokens. xtask `check-unsafe` enforces registration,
//! size, and SAFETY.md coverage for both capsules.
#![allow(unsafe_code)]
// registered capsule — see SAFETY.md beside this file
// `&self -> &mut` is the load-bearing arena shape (as in bumpalo): every
// call returns a DISTINCT region by the monotone-bump invariant, so the
// returned exclusive borrows never alias (SAFETY.md, Aliasing assumptions).
#![allow(clippy::mut_from_ref)]

use core::cell::{Cell, RefCell};
use core::mem::{align_of, needs_drop, size_of};
use core::ptr::NonNull;
use std::alloc::{Layout, alloc, dealloc};

mod poison;

/// One block of memory obtained from the global allocator. Exclusively owned:
/// the block is freed exactly once, in [`Drop`].
pub(crate) struct Chunk {
    base: NonNull<u8>,
    layout: Layout,
}

impl Chunk {
    /// Allocate `bytes` with the given base alignment (power of two).
    /// Returns `None` when the layout is invalid (zero/overflowing size,
    /// bad alignment) or the OS refuses the allocation. Never aborts.
    pub(crate) fn allocate(bytes: usize, base_align: usize) -> Option<Chunk> {
        if bytes == 0 {
            return None; // empty chunks are never useful; refuse early
        }
        let layout = Layout::from_size_align(bytes, base_align).ok()?;
        // SAFETY: `layout` has non-zero size (checked above).
        let p = unsafe { alloc(layout) };
        NonNull::new(p).map(|base| Chunk { base, layout })
    }

    /// Usable size in bytes.
    pub(crate) fn len(&self) -> usize {
        self.layout.size()
    }

    /// Base address (for alignment assertions in tests; never dereferenced
    /// outside this capsule).
    pub(crate) fn base_addr(&self) -> usize {
        self.base.as_ptr() as usize
    }
}

impl Drop for Chunk {
    fn drop(&mut self) {
        // SAFETY: `base` came from `alloc` with exactly `self.layout`, is
        // freed exactly once (Chunk is not Clone/Copy), and no references
        // into the block outlive the owning arena (see SAFETY.md, Lifetime
        // assumptions).
        unsafe { dealloc(self.base.as_ptr(), self.layout) }
    }
}

// SAFETY: a Chunk exclusively owns its block; moving it to another thread
// transfers that ownership. `Sync` is deliberately NOT implemented — the
// bump offset lives outside the chunk, in RawArena, which is !Sync.
unsafe impl Send for Chunk {}

/// "The current window is full" — the facade reacts by installing a larger
/// chunk and retrying. Carries no payload; sizing context lives in the
/// facade where the budget is known.
#[derive(Debug)]
pub(crate) struct Full;

/// Multi-chunk bump arena core. Single-threaded by construction (`Cell`,
/// `RefCell`), so the whole type is `!Sync`; it is `Send` because `Chunk`
/// ownership transfers cleanly.
///
/// Invariant (window): `cur <= end`, and when `cur < end` both lie within
/// the LAST chunk pushed onto `chunks`. Placement never touches any earlier
/// chunk; earlier chunks only stay alive so outstanding references into
/// them remain valid.
pub(crate) struct RawArena {
    chunks: RefCell<Vec<Chunk>>,
    cur: Cell<usize>,
    end: Cell<usize>,
}

impl RawArena {
    /// New arena with an empty window; every placement fails with [`Full`]
    /// until [`Self::install_chunk`] provides memory.
    pub(crate) fn new() -> Self {
        RawArena {
            chunks: RefCell::new(Vec::new()),
            cur: Cell::new(0),
            end: Cell::new(0),
        }
    }

    /// Make `chunk` the current bump window. Remaining space in the previous
    /// window is abandoned (bounded waste; the facade grows geometrically).
    ///
    /// Pushing moves the `Chunk` STRUCT (pointer + layout), never the block
    /// it owns, so references previously handed out remain valid.
    pub(crate) fn install_chunk(&self, chunk: Chunk) {
        let base = chunk.base_addr();
        self.cur.set(base);
        self.end.set(base + chunk.len());
        self.chunks.borrow_mut().push(chunk);
    }

    /// Number of chunks currently owned.
    pub(crate) fn chunk_count(&self) -> usize {
        self.chunks.borrow().len()
    }

    /// Base address of the most recently installed chunk (test/verification
    /// hook for the hugepage-alignment decision). `None` before any chunk.
    pub(crate) fn last_chunk_base(&self) -> Option<usize> {
        self.chunks.borrow().last().map(Chunk::base_addr)
    }

    /// Strip all chunks for recycling. `&mut self` is the load-bearing part:
    /// the borrow checker guarantees no reference handed out by any
    /// `try_place*` call (all tied to `&self`) is still alive, so the blocks
    /// may be reused or freed. The window is reset to empty.
    pub(crate) fn take_chunks(&mut self) -> Vec<Chunk> {
        self.cur.set(0);
        self.end.set(0);
        core::mem::take(self.chunks.get_mut())
    }

    /// Bump-allocate `size` bytes at `align` (raised to the crate-wide
    /// 128-byte policy floor) from the current window.
    fn bump(&self, size: usize, align: usize) -> Result<NonNull<u8>, Full> {
        debug_assert!(align.is_power_of_two());
        let align = align.max(crate::ALLOC_ALIGN);
        let cur = self.cur.get();
        let aligned = cur.checked_add(align - 1).ok_or(Full)? & !(align - 1);
        let next = aligned.checked_add(size).ok_or(Full)?;
        if next > self.end.get() {
            return Err(Full);
        }
        self.cur.set(next);
        // SAFETY: `aligned..next` lies inside the current chunk (window
        // invariant), `aligned` is non-zero (chunk bases are non-null and
        // in-window addresses are >= base), and the monotone bump pointer
        // guarantees the range was never handed out before.
        Ok(unsafe { NonNull::new_unchecked(aligned as *mut u8) })
    }

    /// Place one value. On success the exclusive reference lives as long as
    /// the `&self` borrow; on a full window the value is handed back so the
    /// caller can grow and retry. Types needing `Drop` are rejected at
    /// compile time: bump arenas never run destructors.
    pub(crate) fn try_place<T>(&self, value: T) -> Result<&mut T, T> {
        const {
            assert!(
                !needs_drop::<T>(),
                "arena-placed types must not need Drop (bump arenas never run destructors)"
            );
        }
        if size_of::<T>() == 0 {
            // ZSTs occupy no arena bytes; a well-aligned dangling pointer is
            // the canonical valid representation.
            core::mem::forget(value); // no-op: T is a ZST without Drop
            // SAFETY: NonNull::dangling is aligned and valid for ZST access.
            return Ok(unsafe { NonNull::<T>::dangling().as_mut() });
        }
        match self.bump(size_of::<T>(), align_of::<T>()) {
            Ok(p) => {
                let p = p.cast::<T>();
                // SAFETY: `p` is fresh (never aliased), aligned for T, and
                // sized for T by the bump contract above.
                unsafe {
                    p.as_ptr().write(value);
                    Ok(&mut *p.as_ptr())
                }
            }
            Err(Full) => Err(value),
        }
    }

    /// Place a slice of `len` copies of `fill`.
    pub(crate) fn try_place_slice_fill<T: Copy>(
        &self,
        len: usize,
        fill: T,
    ) -> Result<&mut [T], Full> {
        // Copy implies !Drop; no const assert needed.
        self.place_slice_raw(len, |dst: *mut T, i| {
            // SAFETY (closure contract): `dst` is the reserved, aligned,
            // in-bounds element slot `i`, written exactly once.
            unsafe { dst.add(i).write(fill) }
        })
    }

    /// Place a slice built element-by-element from `f(i)`. If `f` panics
    /// midway the reserved bytes stay bumped-but-unreferenced (no leak of
    /// Drop types — T cannot need Drop — and no UB; see SAFETY.md).
    pub(crate) fn try_place_slice_with<T>(
        &self,
        len: usize,
        f: &mut dyn FnMut(usize) -> T,
    ) -> Result<&mut [T], Full> {
        const {
            assert!(
                !needs_drop::<T>(),
                "arena-placed types must not need Drop (bump arenas never run destructors)"
            );
        }
        self.place_slice_raw(len, |dst: *mut T, i| {
            // SAFETY (closure contract): as in try_place_slice_fill.
            unsafe { dst.add(i).write(f(i)) }
        })
    }

    /// Shared slice-placement plumbing. `write_at(dst, i)` must write element
    /// `i` exactly once; it runs for i = 0..len in order.
    fn place_slice_raw<T>(
        &self,
        len: usize,
        mut write_at: impl FnMut(*mut T, usize),
    ) -> Result<&mut [T], Full> {
        if len == 0 || size_of::<T>() == 0 {
            // Zero-size cases: valid dangling slice, no arena bytes.
            // SAFETY: dangling is aligned; len elements of a ZST (or zero
            // elements of any T) are valid at any aligned non-null address.
            return Ok(unsafe {
                core::slice::from_raw_parts_mut(NonNull::<T>::dangling().as_ptr(), len)
            });
        }
        let bytes = size_of::<T>().checked_mul(len).ok_or(Full)?;
        let p = self.bump(bytes, align_of::<T>())?.cast::<T>();
        for i in 0..len {
            write_at(p.as_ptr(), i);
        }
        // SAFETY: all `len` elements were just initialized in the freshly
        // bumped, aligned, exclusively-held range.
        Ok(unsafe { core::slice::from_raw_parts_mut(p.as_ptr(), len) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refuses_bad_layouts() {
        assert!(Chunk::allocate(0, 128).is_none());
        assert!(Chunk::allocate(64, 3).is_none()); // not a power of two
        assert!(Chunk::allocate(usize::MAX - 7, 128).is_none());
    }

    #[test]
    fn window_math_and_alignment() {
        let raw = RawArena::new();
        assert!(raw.try_place(1u8).is_err(), "empty window must be Full");
        raw.install_chunk(Chunk::allocate(4096, crate::ALLOC_ALIGN).expect("small chunk"));
        let a = raw.try_place(7u8).expect("fits");
        let b = raw.try_place(9u64).expect("fits");
        assert_eq!(*a, 7);
        assert_eq!(*b, 9);
        assert_eq!(core::ptr::from_mut(a) as usize % crate::ALLOC_ALIGN, 0);
        assert_eq!((core::ptr::from_mut(b) as usize) % crate::ALLOC_ALIGN, 0);
        let s = raw.try_place_slice_fill(8, 1.5f64).expect("fits");
        assert_eq!(s, &[1.5; 8]);
        assert_eq!(s.as_ptr() as usize % crate::ALLOC_ALIGN, 0);
    }

    #[test]
    fn full_window_hands_value_back() {
        let raw = RawArena::new();
        raw.install_chunk(Chunk::allocate(256, crate::ALLOC_ALIGN).expect("chunk"));
        let big = [0u8; 512];
        let Err(returned) = raw.try_place(big) else {
            panic!("512B cannot fit a 256B window");
        };
        assert_eq!(returned.len(), 512);
    }

    #[test]
    fn zero_size_cases_consume_nothing() {
        let raw = RawArena::new();
        // No chunk installed: ZST and empty-slice placements still succeed.
        let unit = raw.try_place(()).expect("ZST needs no memory");
        assert_eq!(*unit, ());
        let empty = raw.try_place_slice_fill(0, 0u8).expect("empty slice");
        assert!(empty.is_empty());
        let zsts = raw.try_place_slice_with(3, &mut |_| ()).expect("ZST slice");
        assert_eq!(zsts.len(), 3);
    }
}
