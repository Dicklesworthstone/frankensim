# SAFETY: fs-alloc/src/raw/mod.rs

> The bump-pointer core named by Decalogue P1 as a sanctioned unsafe zone
> ("arena allocators"). Registered in unsafe-capsules.json; enforced by
> `cargo run -p xtask -- check-unsafe` (registration, <300 lines, this file).

## Invariants
1. A `Chunk` exclusively owns one live global-allocator block described
   exactly by its stored `Layout`; it is freed exactly once, in `Drop`.
2. The bump window `[cur, end)` always lies inside the most recently
   installed chunk; `cur` is monotone non-decreasing between
   `install_chunk`/`take_chunks` calls, so no byte range is ever handed out
   twice.
3. Every placement pointer is aligned to `max(align_of::<T>(), 128)` and
   in-bounds by the window arithmetic (checked_add, no wrap).
4. Placed types never need `Drop` — enforced at compile time by
   `const { assert!(!needs_drop::<T>()) }` in `try_place` /
   `try_place_slice_with` and by the `Copy` bound in
   `try_place_slice_fill`. The arena never runs destructors, so this makes
   "reclaim without drops" sound rather than leaky.

## Aliasing assumptions
Each successful `bump` returns a byte range disjoint from every previously
returned range (monotone bump, invariant 2), so the `&mut` references handed
out never alias each other. The arena's own bookkeeping (`Cell` offsets,
`RefCell<Vec<Chunk>>`) is only accessed inside this capsule's methods and is
never referenced from placed memory.

## Alignment assumptions
Chunk base alignment is whatever the facade requests (>= 128, or 2 MiB for
THP-eligible chunks); `Layout::from_size_align` validates it. In-window
placement alignment is computed by `bump` (invariant 3). This capsule
ESTABLISHES the crate's 128-byte policy floor; nothing upstream is relied on.

## Lifetime assumptions
The one deliberate erasure: references returned by `try_place*` borrow the
`RawArena` (`&self`), while the bytes they point into are owned by `Chunk`s
in the arena's `Vec`. This is sound because chunks are only dropped or
recycled by `take_chunks(&mut self)` or by dropping the arena — both require
exclusive access, which the borrow checker grants only after every
handed-out `&mut` (tied to `&self`) is dead. Pushing to the `Vec` moves
`Chunk` structs (pointer + layout), never the blocks they own, so growth
does not invalidate outstanding references.

## Panic behavior
No unsafe block brackets a user-code call except slice element construction
in `try_place_slice_with`: if `f(i)` panics mid-fill, the bumped range stays
reserved-but-unreferenced (offset already advanced), no reference has
escaped, and `T` cannot need `Drop` — so unwinding leaks at most arena bytes
until scope reclaim, never touching freed or uninitialized memory through a
live reference. `Chunk::allocate` returns `None` on failure rather than
calling `handle_alloc_error` (no aborts mid-campaign — CONVENTIONS.md).

## Cancellation behavior
No poll points: every method is a bounded, lock-free, allocation-free
operation (chunk allocation aside, which is a single global-allocator call).
Cancellation interacts with arenas only through scope teardown in the safe
facade, which drops or recycles whole arenas; invariant 1 makes that
reclamation correct regardless of when cancellation lands.

## Concurrency behavior
`RawArena` is `!Sync` by construction (`Cell`/`RefCell`), so all placement
is single-threaded; there is nothing to order. `unsafe impl Send for Chunk`
is justified because a chunk is an exclusively owned heap block with no
thread affinity — sending it transfers ownership, exactly like `Box<[u8]>`.
`RawArena` derives `Send` from its fields for the per-worker-arena pattern.

## Miri coverage
The whole capsule is plain Rust (no intrinsics, no FFI, no syscalls), so
Miri interprets every path. `cargo miri test -p fs-alloc` exercises the
in-module tests plus the facade unit tests and doctests that route through
this capsule (placement, growth, ZSTs, full-window handback, reclaim).
Known limitation: none specific to this capsule.

## Model-checking coverage
N/A (single-threaded capsule; the concurrent structures in fs-alloc —
chunk free list, sharded pools — are plain `Mutex`/atomic safe code outside
this capsule, covered by the G4 storm and multi-thread conformance cases).

## Fuzz/property coverage
The G0 shadow-model battery in tests/conformance.rs drives randomized
(seeded LCG) interleavings of placements, scope creation, and teardown
through the safe facade, asserting accounting equivalence and pairwise
non-overlap of returned ranges. The G4 storm runs 10^6 randomized
cancel/complete cycles and asserts quiescence (zero leaked bytes/arenas).

## Proof obligations discharged by callers
None. The facade-visible API of this capsule is safe: misuse that could
cause UB (overlapping placements, premature chunk reuse) is prevented by
the capsule's own bump discipline and by `&mut`-gated reclamation, not by
caller promises. The only caller-visible contracts are compile-time
(`!needs_drop`) or ordinary panics/`Result`s.
