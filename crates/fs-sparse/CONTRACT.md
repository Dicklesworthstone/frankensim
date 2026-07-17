# CONTRACT: fs-sparse

## Purpose and layer
Sparse matrix formats (CSR, BSR, SELL-C-σ), deterministic COO assembly,
SpMV/SpMM kernels, and pattern algebra (transpose, symmetrize, Gustavson
SpGEMM — the building block of AMG's Galerkin triple product). Layer:
**L1 BEDROCK**. The default numeric core uses only `std`; the off-by-default
`fnx-interop` and `fnp-interop` features add the documented FrankenNetworkx and
FrankenNumpy path dependencies. Plan §6.2.

v1 includes the scalar correctness core plus compact-index CSR, deterministic
sharded CSR/SELL SpMV, tiled parallel COO assembly, blocked SpMM, sparse-SPA
SpGEMM, and a runtime-dispatched x86 AVX2+FMA code-generation capsule. The
ignored release-only roofline harness reports attainment against measured
STREAM and enforces its 85% gate only when `FS_SPARSE_ROOFLINE_GATE=1`;
results are machine-specific, not a universal throughput guarantee.

## Public types and semantics
- `Coo` — triplet staging; duplicates ACCUMULATE (FEM element-assembly
  contract). `assemble()` produces canonical CSR.
- `Csr` — the canonical format. INVARIANT: within each row, columns are
  strictly ascending, no duplicates. Every constructor establishes this
  (`from_parts` validates and panics otherwise); every algorithm may rely on
  it. `try_from_parts_with_checkpoint` validates the same invariant without a
  panic and invokes a caller callback once per row and stored column, allowing
  cancellation-aware callers without adding an executor dependency. `row(r)`
  doubles as the graph neighbor view (CSR IS the adjacency structure). `spmv`,
  `spmm` (dense row-major right-hand sides), `to_dense` (oracle use),
  `identity`, `get`.
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
allocation-failure handling beyond std's. The checkpointed CSR constructor
instead returns `Ok(None)` for malformed canonical parts and passes through a
caller checkpoint error as `Err`.

## Determinism class
**Bit-deterministic cross-ISA by construction**: kernels are fixed-order
+, ×, mul_add with no libm or data-dependent reassociation; parallel assembly
and sharded kernels preserve their serial twins' logical accumulation order.
Evidence: the conformance battery (three-matrix zoo × three formats ×
transpose/symmetrize/SpGEMM) folds all output bits into FNV-64 golden hash
`0xbcf5_52b6_c5bf_aed6`; the preconditioner battery (Chebyshev apply +
ILU-PCG + AMG-PCG solutions + hierarchy shapes) hashes to
`0x752f_215a_26e3_2fea`. Both recorded on aarch64-apple (M4 Pro) and
verified identical on x86-64 (Threadripper). Golden-evidence policy
applies. NO platform libm feeds any solver state (workspace contract
rule).

## Cancellation behavior
The serial and `std`-threaded numeric kernels do not accept an asupersync `Cx`
and are uninterruptible once called. Scoped parallel lanes join/drain their
workers before returning. After top-level shape checks,
`Csr::try_from_parts_with_checkpoint` invokes a caller checkpoint once per row
and once per visited stored column. No bounded-latency cancellation claim is
made for SpMV, SpMM, SpGEMM, parallel assembly, or preconditioners.

## Unsafe boundary
One registered capsule exists at `src/fma/mod.rs`, with its calling invariants
recorded in `src/fma/SAFETY.md` and `unsafe-capsules.json`. On x86-64, safe
dispatchers verify AVX2+FMA support before calling private `#[target_feature]`
functions; the unsafe scope is the target-feature calling contract, while the
kernel bodies remain safe slice arithmetic. Workspace-level `unsafe_code =
"deny"` remains in force outside that module-local capsule.

## Feature flags
- `fnx-interop` (default OFF) — pulls the optional `fnx-classes` + `fnx-runtime`
  path deps and enables the `interop` module (Csr ⇄ FrankenNetworkx
  `GraphSnapshot`). Off by default so the L1 crate stays dependency-lean; the
  numeric core never pulls a constellation crate.
- `fnp-interop` (default OFF) — pulls the optional `fnp-ufunc` + `fnp-dtype`
  path deps and enables the `interop_fnp` module (Csr ⇄ FrankenNumpy
  `UFuncArray`). Same dependency-lean rationale.

## Conformance tests
`tests/conformance.rs`: cross-format bitwise battery + golden hash. In-crate
suites: assembly canonicalization + stream-order invariance, SpMV vs dense
oracle, linearity, adversarial patterns (empty rows, dense row, single
column, empty matrix), BSR/SELL round-trips, SELL padding economics,
checkpointed CSR publication and malformed canonical-part refusal,
transpose involution, symmetrize bitwise symmetry, SpGEMM vs dense oracle +
Laplacian-square pattern sanity, sparse-SPA SpGEMM vs dense-SPA reference on
random and 2e6-column-wide products, structured rejections. Any
reimplementation must pass the conformance battery bit-for-bit.
Bead 4nh8 adds 600 shrink-armed generated 8×8 cases (seed `0x5A_5001`): up
to 64 integer COO triplets, including duplicates and stored zeros, are applied
to an integer vector through CSR, BSR 4×4, and SELL-C-σ `(8, 32)` and compared
bitwise. The fixed cross-ISA golden `0xbcf5_52b6_c5bf_aed6` is unchanged.

Bead 6ys.18.7 adds the cheap structured PR subset through `fs-casebook`: an
exact out-of-order COO assembly KAT with duplicate chains, exact CSR/BSR/SELL
SpMV agreement against an independent literal result, checkpointed canonical
publication plus typed and malformed-input refusals, pinned input digests, and
a disclosed one-bit seeded corruption proving structured failure reporting and
merge-gate refusal. It preserves the aggregate golden and the broader
property, WSBF, preconditioner, release-performance, and cross-ISA lanes as
separate evidence; this subset awards no fresh oracle, performance, nightly,
or dual-ISA claim.

`tests/frankenscipy_oracle_casebook.rs` supplies the missing test-only
FrankenScipy sparse baseline. Identical canonical CSR storage is admitted by
`fs-sparse::Csr` and `fsci_sparse::CsrMatrix`; a 4×5 dyadic KAT pins both
outputs bit-for-bit, including an empty-row positive zero. Twelve disclosed
LCG fixtures (root `0x5A25_5A25_D15C_0001`, 84 output rows) compare the fused
FrankenSim path with FrankenScipy's separately rounded matvec under the
declared absolute oracle-agreement bound `32·f64::EPSILON`. Canonical
little-endian frames retain every row pointer, column, value, vector,
implementation/API version, tolerance, and expected KAT bits. The input
frame is sufficient to regenerate both outputs; the passing structured record
retains separate FNV digests of all 84 computed output-bit rows, while any
failing row logs both implementations' exact bits and the full frame. The
input digests are `d7c1b8ae30777f71` (KAT), `735ea2d324e967bc` (seeded), and
`3d1185b7e82a131f` (shape-refusal policy). Seed `0x5A25_0000` flips one
reference bit, producing one replay-identical red record at digest
`323f3591ddcb0d6e` and an `assert_green` refusal. The structural case records
fs-sparse's documented actionable programmer-error panic and FrankenScipy's
typed `IncompatibleShape` refusal, with the fs-sparse output unchanged.
`fsci-sparse` remains a development dependency only; this tranche is central
package-proof pending.

## FrankenNumpy interop (bead gtql item c, feature `fnp-interop`)

Scout verdict recorded: fnp-ndarray holds only layout metadata; the
owned array type is `fnp_ufunc::UFuncArray` (owned `Vec<f64>` +
shape), so these are CONVERSIONS by necessity — no zero-copy view
exists into another crate's owned storage; the borrowed views remain
`Csr::row`/`Csr::to_dense` on this side. `csr_to_dense_array`
densifies (O(nrows·ncols), overflow-refused) with unstored entries as
+0.0 — explicit stored zeros are documented lossiness (round trip is
the identity exactly for CSRs without them). `dense_array_to_csr`
accepts rank-2 only, reads the f64 value plane, drops ±0.0, and
REFUSES non-finite entries with their position (fail closed). Off by
default; the L1 core pulls no constellation crate unless opted in.

## No-claim boundaries
- **Performance scope**: `tests/roofline.rs` reports machine-specific STREAM
  attainment and asserts the ≥85% all-core gate only under
  `FS_SPARSE_ROOFLINE_GATE=1`; no every-host throughput guarantee is made.
  Software-prefetch/autotuner integration and fs-tilelang-generated bodies are
  not present. The shipped x86 acceleration is runtime AVX2+FMA codegen for the
  existing fixed-order bodies.
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
- **Interop scope**: `fnx-interop` copies between square CSR adjacency and the
  owned `GraphSnapshot`; `fnp-interop` converts between CSR and the owned dense
  `UFuncArray`, requiring O(nrows·ncols) memory when densifying and losing the
  distinction between explicit and unstored ±0.0. `Csr::row` remains the
  borrowed sparse neighbor/value view; no zero-copy FrankenNetworkx or
  FrankenNumpy view type is claimed.
- Canonical `Csr`, `Bsr`, and `Sell` use `usize` indices. `CsrCompact` narrows
  column indices to `u32` and refuses dimensions outside that space; compact
  BSR/SELL indices are not claimed.
- The FrankenScipy Casebook proves agreement for its fixed canonical KAT and
  disclosed seeded finite fixture family. It does not establish Python SciPy
  execution, arbitrary-matrix equivalence, performance parity, fresh
  dual-ISA execution, or the umbrella's full G5 sweep. The `32·EPSILON`
  margin is an absolute oracle-agreement bound for the fixture's five terms per row and
  `|a_ij|, |x_j| ≤ 1/2`; it is not a dimension-free tolerance for callers.
