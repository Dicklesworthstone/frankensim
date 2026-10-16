# SAFETY: fs-lbm D3Q19 pull-stream AVX2 capsule

## Invariants

`kernel_256` moves one four-cell x row per AVX2 register. Loads, permutations,
blends, and stores perform no floating-point arithmetic, so every population
bit is retained. The safe façade validates positive 4-aligned dimensions and
equal tile counts for all 19 input and output fields. The private selected
thunk additionally asserts the frozen four-cell row extent before entering the
capsule. Its one-shot selector publishes the thunk only after AVX2 detection.
All derived tile indices remain within the validated tile grid, and every
pointer access covers exactly one 4-element row inside a 64-element tile.

## Aliasing assumptions

The input field is an immutable Rust borrow and the output field is a distinct
exclusive mutable borrow. Safe Rust cannot construct the call with overlapping
storage, and the capsule does not manufacture aliases between output rows.

## Alignment assumptions

None beyond valid `f64` pointers. The capsule uses unaligned AVX2 loads and
stores; the tile's 128-byte alignment is a performance property, not a
soundness precondition.

## Lifetime assumptions

No lifetime is erased or reconstructed. Every raw pointer is derived from an
indexed tile borrow and is consumed before that borrow can end; no pointer
escapes the capsule.

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
discharge aliasing and lifetime obligations, and the private selector performs
CPU-feature admission before the capsule can run.
