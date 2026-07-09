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
- `ops::{transpose, symmetrize, spgemm, spgemm_sparse_spa}` — pattern
  algebra on canonical CSR. `spgemm` is the dense-SPA reference path;
  `spgemm_sparse_spa` uses a per-row BTree accumulator for very wide products
  without O(ncols) scratch, while preserving the reference accumulation order
  and ascending-column output.
- `interop::{csr_to_graph_snapshot, graph_snapshot_to_csr, WEIGHT_KEY,
  InteropError}` (feature `fnx-interop`, bead gtql) — FrankenNetworkx bridge.
  `Csr → GraphSnapshot`: square adjacency only (else `NotSquare`); node `i` ⇄
  key `"i"`, stored `(r, c, v)` ⇄ directed `EdgeSnapshot` `r→c` with weight `v`
  under `WEIGHT_KEY`, emitted in canonical row-major/ascending-column order.
  `GraphSnapshot → Csr`: nodes map to indices BY ORDER (keys may be arbitrary
  strings), weights read from `WEIGHT_KEY` (`Float`/`Int` coerced, else `1.0`),
  parallel `(r, c)` edges SUM (multigraph → simple). Necessarily COPIES (fnx's
  `Graph` is string-keyed and owned — no zero-copy path; the honest deviation
  from "wrap, do not copy"). Round-trip `Csr → snapshot → Csr` is the identity,
  bitwise on values.
- `precond::{Precond, IdentityPrecond, Chebyshev, Ilu0/ilu0/IluBreakdown,
  SaAmg, pcg/PcgReport, lambda_max_estimate}` — the solver-stack toolkit.
  `Precond::apply(r, z)` is the operator interface. Chebyshev: degree-k
  three-term recurrence on the band [λmax/α, λmax], λmax by deterministic
  power iteration (fixed libm-free start, 1.1 safety — enclosure tested
  vs analytic spectra). ILU(0): zero-fill on the CSR pattern, sequential
  v1, typed `IluBreakdown{row}` (shift-retry is caller policy). SaAmg:
  symmetric strength graph, GREEDY INDEX-ORDER aggregation with
  lowest-index tie-breaks (P2 on setup), Jacobi-smoothed prolongator
  (ω = 4/3·λmax(D⁻¹A)), Galerkin triple product via in-crate SpGEMM,
  V-cycle with Chebyshev smoothing and ILU-PCG coarsest solve;
  `operator_complexity()` and `level_sizes` are the memory-honesty
  evidence. `pcg` is a REFERENCE driver (solver stack supersedes it).

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
7. `spgemm_sparse_spa(A, B) = spgemm(A, B)` bitwise on tested random and
   very-wide products; it changes accumulator storage only, not numerical
   order or output canonicalization.
8. Preconditioner setup and solves are rerun-deterministic BITWISE
   (hierarchy shapes, iteration counts, and solutions — tested), and the
   spectral-bound estimate ENCLOSES the true λmax on tested fixtures
   (over-estimation safe by construction, safety factor 1.1).
9. AMG on 2D Poisson: near-grid-independent PCG iterations (32² vs 64²
   within a tested band), operator complexity < 2, and the anisotropic
   ε = 1e-3 fixture converges (tested).

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
`0xbcf5_52b6_c5bf_aed6`; the preconditioner battery (Chebyshev apply +
ILU-PCG + AMG-PCG solutions + hierarchy shapes) hashes to
`0x752f_215a_26e3_2fea`. Both recorded on aarch64-apple (M4 Pro) and
verified identical on x86-64 (Threadripper). Golden-evidence policy
applies. NO platform libm feeds any solver state (workspace contract
rule).

## Cancellation behavior
v1 kernels are single-tile and uninterruptible; the executor-tiled parallel
lanes (follow-up bead) will poll at row-range boundaries per Decalogue P7.

## Unsafe boundary
None. `unsafe_code` denied; no capsules.

## Feature flags
- `fnx-interop` (default OFF) — pulls the optional `fnx-classes` + `fnx-runtime`
  path deps and enables the `interop` module (Csr ⇄ FrankenNetworkx
  `GraphSnapshot`). Off by default so the L1 crate stays dependency-lean; the
  numeric core never pulls a constellation crate.

## Conformance tests
`tests/conformance.rs`: cross-format bitwise battery + golden hash. In-crate
suites: assembly canonicalization + stream-order invariance, SpMV vs dense
oracle, linearity, adversarial patterns (empty rows, dense row, single
column, empty matrix), BSR/SELL round-trips, SELL padding economics,
transpose involution, symmetrize bitwise symmetry, SpGEMM vs dense oracle +
Laplacian-square pattern sanity, sparse-SPA SpGEMM vs dense-SPA reference on
random and 2e6-column-wide products, structured rejections. Any
reimplementation must pass the conformance battery bit-for-bit.

## No-claim boundaries
- **No performance claims yet**: scalar reference kernels; the ≥85% STREAM
  target, CCD sharding, prefetch, and SIMD belong to the perf follow-up.
- No parallel assembly implementation (the CONTRACT for its accumulation
  order is stated in Invariant 2; the tiled implementation is follow-up).
- BSR `to_csr` is only structurally lossless for matrices without stored
  exact-zero values (fill is dropped by value test); the dense expansion is
  always bitwise faithful.
- `spgemm_sparse_spa` is a deterministic scalar BTree-accumulator path, not a
  tuned hash-SPA throughput path; wide-matrix memory shape is improved, but no
  speedup claim is made.
- ILU(0) is sequential (level scheduling recorded); IC(0)-specific
  symmetric storage is unclaimed (ILU covers SPD use). Supernodal
  Cholesky deferred per its own scope cap. AMG coarsest solve is
  ILU-PCG (dense direct coarse solve joins solver-stack integration).
  No 1e8-DOF scaling claims yet (release-mode scaling lane).
- **Interop scope (bead gtql)**: the `fnx-interop` feature ships the
  FrankenNetworkx `GraphSnapshot` bridge (round-trip tested). FrankenNumpy
  array views (item c) are a scoped follow-up — the owned f64 array lives in
  `fnp-ufunc::UFuncArray`, which transitively pulls `rayon` + `fnp-linalg`, a
  heavy dependency for an L1 numeric core; the borrow-vs-convert decision and a
  separate `fnp-interop` feature are deferred until a consumer needs it.
- No FrankenNumpy/FrankenNetworkx interop views yet (follow-up).
- Indices are `usize` (compact u32 indices are a recorded perf-bead item).
