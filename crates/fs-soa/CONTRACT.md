# fs-soa CONTRACT

## Purpose and layer

Layer: **L0 SUBSTRATE** (deps: fs-soa-derive only). The structure-of-
arrays runtime `#[derive(Soa)]` targets (plan §5.3): per-field growable
buffers aligned to 128 bytes, AoS gather/scatter, chunked SIMD-friendly
access with explicit masked tails, zero-copy strided view descriptors
for the FrankenNumpy membrane (§12), and chunk-quantum grouping as the
tile-identity hook. "SoA everywhere hot" is a memory-discipline pillar:
batched small dense LA (6ys.4) and LBM lattices need SIMD lanes running
ACROSS elements.

## Public types and semantics

- `Soa` — re-export of the derive macro (`use fs_soa::Soa;` is the
  whole user surface).
- `FieldBuf<T: Copy>` — one field's growable aligned buffer: `new`,
  `with_capacity` (hint, honored at first push), `len`, `is_empty`,
  `capacity`, `as_slice`, `as_mut_slice`, `clear` (keeps allocation),
  `reserve`, `push`, `view(name) -> RawView`. Amortized-doubling
  growth; clear/reuse never shrinks.
- `RawView` — leaf view descriptor: dotted `name`, `addr` (0 when
  unallocated), `len`, `elem_bytes`, `stride_bytes` (== `elem_bytes`,
  dense), `achieved_align`, `dtype` (type name, auditability only).
  `descr()` renders an address-free JSON line (deterministic logs).
- `SoaAble` (`type Soa`) and `SoaContainer<T>` (`c_new`,
  `c_with_capacity`, `c_len`, `c_push`, `c_get`, `c_set`, `c_clear`,
  `c_reserve`, `c_views`, `c_layout`) — the composition traits behind
  `#[soa(nested)]`; generated containers implement both.
- `SOA_ALIGN = 128` (matches `fs_alloc::ALLOC_ALIGN` doctrine),
  `DEFAULT_CHUNK_QUANTUM = 512` (E8 tile volume, 8³),
  `chunks_with_tail`/`chunks_with_tail_mut` (chunks_exact + explicit
  remainder), `chunk_count`, `view_name`, `leaf_layout`.

## Invariants

- NO unsafe anywhere: alignment comes from over-allocating the backing
  `Vec` by 128 slack elements and starting the payload at the first
  whole-element offset reaching 128 bytes (`align_offset`); slack slots
  hold copies of pushed values and are never exposed.
- Alignment guarantee: 128 bytes for every element size s where a
  solution of (base + k·s) ≡ 0 (mod 128) exists — always true for
  primitives (s ≤ 16 given ≥16-byte allocator bases) and for s = 24
  ([f64; 3]); for exotic sizes the buffer degrades to the best
  reachable power of two and `RawView.achieved_align` REPORTS it —
  consumers never guess.
- Within a generated container all field buffers have equal length;
  gather(get)/scatter(set/push) preserve bit patterns exactly.
- `chunks_with_tail(s, w)` covers every element exactly once:
  `len/w` full chunks plus a `len % w` tail.

## Error model

Panics on out-of-bounds indices (slice indexing semantics) and on
`width == 0`/`quantum == 0` for the chunk helpers. Allocation failure
aborts as std does. No Result-based paths; containers are plain data.

## Determinism class

Fully deterministic value semantics: layout offsets depend on allocator
base addresses, but no VALUE ever depends on them; view `descr()` and
`layout_descr()` exclude addresses so logs and goldens are stable
across runs and ISAs. No floating-point arithmetic in this crate.

## Cancellation behavior

All operations are synchronous, bounded, and allocation-light (one
buffer copy on growth). Nothing long-running; batch chunking to tile
quanta with Cx poll points is the CONSUMER's job (recorded for 6ys.4).

## Unsafe boundary

None. `unsafe_code = "deny"` via workspace lints; no capsules. This is
load-bearing: the crate proves SoA needs no unsafe.

## Feature flags

None.

## Conformance tests

`tests/soa_battery.rs` (12 cases, RFC 8259 JSON-line logging): 128-byte alignment
asserted for every leaf view on a 4-field fixture ([f64;3]/f64/u32);
AoS↔SoA round-trip bitwise over 500 random elements + iterator
equivalence vs a `Vec` reference; scatter/column-mutation/clear-reuse;
capacity-hint no-regrowth (address-stable across 256 pushes); nested
containers (dotted view paths `inner.a`, drill-down accessors);
generic struct with bounds + const param (`GenericSoa<f32, 4>`);
Qty-typed columns stay dimensionally typed; chunked access with masked
tail (128×8 + 6, mutable pass scales the tail too, quantum groups);
exact address-free layout/view descriptions; RFC 8259 escaping of hostile
descriptor strings and of the battery envelope's dynamic case/verdict/detail
fields (including a descriptor nested as JSON-string detail); 20k-op random
property battery mirrored against `Vec` (fs-rand keyed streams).
`tests/compile_fail.rs`: 8-case diagnostics battery via the in-house
offline harness (scratch cargo project, path deps, no trybuild) —
tuple/unit/enum/lifetime/generic-default/zero-field/reserved-name/
unknown-attribute all rejected with our structured messages.

## No-claim boundaries

- Descriptors only for FrankenNumpy: `RawView` is the membrane shape,
  but LIVE franken_numpy wiring ships when the membrane crate lands
  (no fnp dependency exists anywhere in the workspace yet).
- No arena-backed variant (fs-alloc `Arena` returns borrowed slices;
  growable owned buffers don't fit a bump arena) — revisit if a
  fixed-capacity arena SoA is needed.
- No Morton SORTING of contents; chunk-quantum grouping provides tile
  identities, spatial ordering belongs to consumers (fs-substrate
  morton is the tool).
- No autotuned interleave width (that is 6ys.4's autotuner) and no
  SIMD kernels here — this crate provides the layout, fs-simd/6ys.4
  provide the compute.
- Alignment above 128 bytes (hugepage-granularity placement) is
  fs-alloc territory.
