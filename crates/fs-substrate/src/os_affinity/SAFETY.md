# SAFETY: fs-substrate os_affinity capsule

## What is unsafe here

One `unsafe fn syscall6` (inline `asm!` issuing the Linux x86-64
`syscall` instruction) plus three `unsafe { }` call sites, all in the
`imp` module compiled only for `target_os = "linux", target_arch =
"x86_64"`. Every other target is 100% safe stub code returning
`Unsupported`.

## Why it is sound

- The `asm!` block declares the exact Linux x86-64 syscall ABI:
  arguments in rdi/rsi/rdx/r10/r8/r9, number/result in rax, rcx and
  r11 declared clobbered, `nostack`.
- All pointers passed to the kernel name live, correctly-sized Rust
  stack/heap allocations owned by the calling frame, and every call is
  synchronous — the kernel does not retain the pointers past return.
  - `sched_setaffinity`: `&[u64; 16]` (1024-bit cpu_set_t) read by the
    kernel; size passed explicitly.
  - `getcpu`: two `&mut u32` out-slots; tcache NULL (modern form).
  - `move_pages` (QUERY form, nodes = NULL, flags = 0): a `Vec<u64>`
    of page addresses (read) and a `Vec<i32>` status array (written)
    both sized to `count`; nothing is migrated.
- Failures are returned as negative errnos and surfaced as structured
  `OsAffinityError::Syscall` — no error state is swallowed.
- Affinity is ADVISORY: results of computations never depend on
  placement (determinism P2); a lost race with an external affinity
  change degrades timing only.

## Blast radius

Wrong register constraints would corrupt at most the calling thread's
state at the call site (caught by the capsule tests, which pin to the
reported CPU and re-query). The capsule never transmutes, never frees,
never aliases: kernel-owned copies only.
