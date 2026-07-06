# CONTRACT: fs-sparse

## Purpose and layer
Sparse matrix formats (CSR, BSR, SELL-C-σ), deterministic COO assembly,
SpMV/SpMM kernels, and pattern algebra (transpose, symmetrize, Gustavson
SpGEMM — the building block of AMG's Galerkin triple product). Layer:
**L1 BEDROCK**. Zero dependencies (pure std arithmetic). Plan §6.2.

v1 is the correctness core on scalar kernels. The roofline lane (≥85% of
measured STREAM, per-CCD sharding, prefetch autotuning, fs-tilelang SIMD
bodies, FrankenNetworkx graph interop) is the recorded follow-up bead,
gated on fs-tilelang and the autotuner — see No-claim boundaries.

## Public types and semantics
- `Coo` — triplet staging; duplicates ACCUMULATE (FEM element-assembly
  contract). `assemble()` produces canonical CSR.
- `Csr` — the canonical format. INVARIANT: within each row, columns are
  strictly ascending, no duplicates. Every constructor establishes this
  (`from_parts` validates and panics otherwise); every algorithm may rely on
  it. `row(r)` doubles as the graph neighbor view (CSR IS the adjacency
  structure). `spmv`, `spmm` (dense row-major right-hand sides), `to_dense`
  (oracle use), `identity`, `get`.
- `Bsr` — block CSR with fixed r×c blocks for FEM vector unknowns.
  `from_csr` requires divisible dimensions (padding is a modeling decision,
  never invented by a conversion). `to_csr` drops exact-zero fill.
- `Sell` — SELL-C-σ, stable-sorted by descending row length within σ-row
  windows, lane-fastest layout. Stores TRUE per-row lengths; pad slots exist
  physically but are never read. `to_csr` is bitwise lossless.
- `ops::{transpose, symmetrize, spgemm}` — pattern algebra on canonical CSR.

## Invariants
1. **Cross-format bitwise SpMV equality**: CSR, BSR, and SELL SpMV produce
   BIT-IDENTICAL outputs. Mechanism: every kernel accumulates each row in
   ascending-global-column order with fused `mul_add` from a +0.0 start; BSR
   fill zeros are provably inert (a fused `0·x + acc` is exactly `acc`, and
   `acc` cannot become −0.0 from a +0.0 start under round-to-nearest); SELL
   never reads pads. Tested on FEM, random, and skewed fixtures.
2. **Assembly canonicalization**: `Coo::assemble` output is a pure function
   of the triplet multiset ordered by (row, col, insertion sequence) —
   stream/tile interleavings that preserve each (row, col)'s own contribution
   order produce bitwise-identical matrices (tested with shuffled streams).
   Contribution order within one (row, col) pair is LOGICAL identity (e.g.
   element id); callers parallelizing assembly must preserve it.
3. `spmm` output equals column-by-column `spmv` bitwise (tested).
4. `transpose` is a bitwise involution: `(Aᵀ)ᵀ = A` exactly; values are
   moved, never recomputed.
5. `symmetrize` output is bitwise symmetric (`Sᵀ = S` exactly; IEEE
   `midpoint` commutes) and fixes symmetric inputs.
6. `spgemm(A, I) = A` and `spgemm(I, A) = A` bitwise; contributions to each
   C[i][j] accumulate in ascending-k order (deterministic).

## Error model
Structural violations panic with structured messages: out-of-range COO
indices, non-canonical `from_parts` input, dimension mismatches in
spmv/spmm/spgemm/symmetrize, indivisible BSR block shapes. These are
programmer errors; silently proceeding would void determinism claims. No
other fallible paths; no allocation-failure handling beyond std's.

## Determinism class
**Bit-deterministic cross-ISA by construction**: kernels are fixed-order
+, ×, mul_add; no libm, no data-dependent reassociation, no threading in v1.
Evidence: the conformance battery (three-matrix zoo × three formats ×
transpose/symmetrize/SpGEMM) folds all output bits into FNV-64 golden hash
`0xbcf5_52b6_c5bf_aed6`, recorded on aarch64-apple (M4 Pro) and required to
match on x86-64 (Threadripper). Golden-evidence policy applies.

## Cancellation behavior
v1 kernels are single-tile and uninterruptible; the executor-tiled parallel
lanes (follow-up bead) will poll at row-range boundaries per Decalogue P7.

## Unsafe boundary
None. `unsafe_code` denied; no capsules.

## Feature flags
None.

## Conformance tests
`tests/conformance.rs`: cross-format bitwise battery + golden hash. In-crate
suites: assembly canonicalization + stream-order invariance, SpMV vs dense
oracle, linearity, adversarial patterns (empty rows, dense row, single
column, empty matrix), BSR/SELL round-trips, SELL padding economics,
transpose involution, symmetrize bitwise symmetry, SpGEMM vs dense oracle +
Laplacian-square pattern sanity, structured rejections. Any reimplementation
must pass the conformance battery bit-for-bit.

## No-claim boundaries
- **No performance claims yet**: scalar reference kernels; the ≥85% STREAM
  target, CCD sharding, prefetch, and SIMD belong to the perf follow-up.
- No parallel assembly implementation (the CONTRACT for its accumulation
  order is stated in Invariant 2; the tiled implementation is follow-up).
- BSR `to_csr` is only structurally lossless for matrices without stored
  exact-zero values (fill is dropped by value test); the dense expansion is
  always bitwise faithful.
- SpGEMM uses a dense SPA per row (O(ncols) scratch); hash-SPA for very wide
  matrices is unclaimed.
- No AMG components yet (SpGEMM is the substrate; AMG has its own bead).
- No FrankenNumpy/FrankenNetworkx interop views yet (follow-up).
- Indices are `usize` (compact u32 indices are a recorded perf-bead item).
