# SAFETY: fs-substrate prefetch capsule

## What is unsafe here

Two one-line `unsafe` blocks:
1. `slice.as_ptr().add(idx)` with `idx < slice.len()` checked
   immediately above — an in-bounds pointer into a live slice.
2. The prefetch instruction itself: `_mm_prefetch::<_MM_HINT_T0>`
   (x86-64 intrinsic marked unsafe by signature) or inline
   `prfm pldl1keep` (aarch64 `asm!`).

## Why it is sound

Prefetch instructions have NO architectural effect on memory state:
they never write, never fault (invalid or unmapped addresses are
ignored by the CPU by specification), never touch flags
(`preserves_flags` declared on the asm), and cannot change any
program-visible value — the only possible effect is cache-line
residency, i.e. timing. The pointer passed is additionally always
in-bounds of a live borrow, which is stronger than the instruction
requires.

## Blast radius

None on correctness by construction (determinism P2: results can
never depend on the hint). Worst case is a useless prefetch —
wasted bandwidth, measured and tuned by the sweep harness rather
than assumed.
