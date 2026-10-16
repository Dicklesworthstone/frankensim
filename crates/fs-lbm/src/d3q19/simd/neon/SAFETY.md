# SAFETY: fs-lbm D3Q19 axial-BGK NEON capsule

## Invariants

`kernel_neon` uses AArch64 NEON load, store, and arithmetic intrinsics over one
frozen 64-cell D3Q19 tile. The adjacent safe selected thunk is its only caller.
NEON is architectural on AArch64; no optional CPU feature is assumed. Offsets
are exactly `0,2,...,62`, so each two-lane access is in bounds.
The safe thunk asserts tile divisibility by two before entering unsafe code,
making future tile-extent drift fail closed.

## Aliasing assumptions

Inputs are immutable Rust borrows and all output fields are exclusive mutable
borrows. Constructing the `[&mut [f64; 64]; 19]` output proves the fields do not
alias for the duration of the call.

## Alignment assumptions

None. Loads and stores are explicitly unaligned; the tile's stronger 128-byte
alignment is a performance fact, never a soundness precondition.

## Lifetime assumptions

No lifetime is erased or reconstructed. All pointers are derived from fixed-size
array references whose Rust lifetimes cover the complete intrinsic call.

## Panic behavior

Tau, force, incoming populations, density, and velocity are checked before the
capsule runs. Outputs are checked in canonical lane/direction order after the
body and return the existing typed collision error on non-finite output. The
destination is private post-collision scratch until the check succeeds. The
safe façade panics before unsafe code if the compile-time tile extent stops
being divisible by two; the intrinsic body contains no deliberate panic path.

## Cancellation behavior

One fixed 64-cell tile is allocation-free and bounded. No cancellation point
exists inside it; the caller polls only between complete tiles.

## Concurrency behavior

The capsule has no shared or static mutable state. Rust borrows give one caller
exclusive output ownership, so independent tiles may execute concurrently.

## Miri coverage

Miri excludes this intrinsic module and the one-shot dispatcher selects the safe
scalar twin. The scalar/reference differential battery supplies compensating
arithmetic coverage.

## Model-checking coverage

N/A: the capsule is a bounded, single-call arithmetic kernel with no shared
state, synchronization, atomics, or scheduler interaction.

## Fuzz/property coverage

The safe scalar twin uses the same per-cell operation tree. A seeded battery
compares every output bit and retains the first differing lane/direction plus
all 19 input population bits. The frozen end-to-end Duct golden is an additional
gate.

## Proof obligations discharged by callers

None. The safe façade validates parameters and incoming macroscopic state;
fixed-size Rust borrows discharge bounds, aliasing, and lifetime obligations.
Architecture admission is compile-time gated to AArch64.
