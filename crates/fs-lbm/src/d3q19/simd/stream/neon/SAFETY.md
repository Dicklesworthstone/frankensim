# SAFETY: fs-lbm D3Q19 pull-stream NEON capsule

## Invariants

`kernel_neon` moves one four-cell x row in two AArch64 NEON registers. Loads,
extracts, duplicates, and stores perform no floating-point arithmetic, so every
population bit is retained. The safe façade validates positive 4-aligned
dimensions and equal tile counts for all 19 input and output fields. The
private selected thunk additionally asserts the frozen row extent before
entering the capsule. NEON is architectural on AArch64. All derived tile
indices remain inside the validated grid, and each pointer access covers one
2-element half-row inside a 64-element tile.

## Aliasing assumptions

The input field is an immutable Rust borrow and the output field is a distinct
exclusive mutable borrow. Safe Rust cannot construct the call with overlapping
storage, and the capsule does not manufacture aliases between output rows.

## Alignment assumptions

No alignment beyond a valid `f64` pointer is required by AArch64 loads and
stores. The tile's 128-byte alignment is a performance property, not a
soundness precondition.

## Lifetime assumptions

No lifetime is erased or reconstructed. Every raw pointer is derived from an
indexed tile borrow and consumed before that borrow can end; no pointer escapes
the capsule.

## Panic behavior

Dimension and field-shape errors panic before unsafe code. Given those checked
preconditions, the source map yields only in-bounds tiles and rows and periodic
z always yields a source. A logic-regression panic during the synchronous move
can leave the private destination field partially written, but unwinding is
memory-safe and publishes no aliased reference.

## Cancellation behavior

The synchronous Duct step has no cancellation point inside the stream sweep.
The capsule holds no resource requiring cleanup; a future scoped caller must
poll between complete steps or split the sweep at deterministic tile rows.

## Concurrency behavior

The capsule has no shared or static mutable state. Rust grants exclusive access
to the complete output field, so this implementation is single-caller; a future
parallel schedule must partition destination rows into disjoint mutable slices.

## Miri coverage

Miri excludes this intrinsic module and the one-shot dispatcher selects the
safe scalar stencil. The scalar/active differential battery is the compensating
memory-access and source-map check.

## Model-checking coverage

N/A: this is a bounded synchronous data-move capsule with no atomics,
synchronization, scheduler interaction, or shared mutable state.

## Fuzz/property coverage

A seeded G0 battery compares every output population bit against the retained
scalar stencil on single-tile and asymmetric 2x3x4-tile domains. Independent
route anchors cover inter-tile x/y/z moves, periodic z, corner bounce, and wall
precedence over a simultaneous periodic source.

## Proof obligations discharged by callers

None. The safe façade validates dimensions and field lengths, Rust borrows
discharge aliasing and lifetime obligations, and architecture gating admits
this capsule only on AArch64.
