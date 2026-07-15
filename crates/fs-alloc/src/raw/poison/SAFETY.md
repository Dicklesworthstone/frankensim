# SAFETY: fs-alloc/src/raw/poison/mod.rs

> Raw byte access for the opt-in reclaimed-chunk poison detector. Registered
> separately from the bump-pointer capsule in `unsafe-capsules.json` and kept
> behind `ArenaPool`'s safe, mutex-serialized diagnostic facade.

## Invariants
1. Every `Chunk` exclusively owns one live allocation for exactly `len()`
   bytes.
2. Poisoning occurs only after the chunk has left an arena, or while an
   exclusively owned pending acquisition is being returned before installation.
   Fault injection and verification occur only while the chunk is exclusively
   owned by the pool free list.
3. Verification reads exactly `len()` initialized poison bytes. Injection
   writes exactly one seed-selected in-bounds byte.

## Aliasing assumptions
The safe facade holds the free-list mutex while verifying or injecting. Arena
reclamation requires exclusive access to the arena, so every typed reference
previously returned from the chunk is dead before poisoning begins. Concurrent
free-list reuse is excluded by the mutex.

## Alignment assumptions
Byte reads/writes require alignment one. No stronger alignment is assumed;
the allocation's original alignment remains unchanged for later reuse.

## Lifetime assumptions
Raw slices exist only for the duration of `reclaimed_poison_mismatch` and never
escape the capsule. The owning `Chunk` stays live and cannot be moved or
dropped while borrowed.

## Panic behavior
The raw operations do not call user code. Injection asserts only in debug
builds against the impossible empty-allocated-chunk state. A detected mismatch
is returned as data; the safe facade quarantines the chunk and returns a
structured `AllocError` rather than panicking.

## Cancellation behavior
Poisoning is enabled only by the explicit diagnostic constructor and runs
during whole-arena reclaim or rollback of an exclusively owned pending
acquisition. There are no partial published states: the chunk is poisoned
before insertion under the free-list mutex.

## Concurrency behavior
All calls are serialized by the pool free-list mutex. This capsule introduces
no atomics and implements neither `Send` nor `Sync`; `Chunk` ownership transfer
is justified by the primary raw capsule.

## Miri coverage
The capsule is plain pointer/slice access over global-allocator blocks, so Miri
can interpret every path. The intended lane is
`cargo miri test -p fs-alloc --test conformance alloc_010_seeded_reclaim_poison_detects_and_quarantines_corruption`;
it routes poison, injection, mismatch detection, quarantine, and healthy reuse
through the safe facade. That integration lane is batch-proof-pending and is
not executed by the existing `cargo miri test -p fs-alloc --lib` command, so no
current green Miri claim is made for this capsule.

## Model-checking coverage
N/A: the unsafe operations are mutex-serialized. The surrounding plain
mutex/atomic pool protocol is exercised by seeded G4 runtime tests, not an
exhaustive schedule model checker.

## Fuzz/property coverage
`alloc-010` uses a retained seed to select the poison pattern and corrupted
offset, and asserts the exact structured mismatch receipt plus quiescence.

## Proof obligations discharged by callers
None for public callers. These methods are crate-private and invoked only while
the safe facade exclusively owns a chunk during reclaim/rollback, or owns a
free-listed chunk under the free-list mutex.
