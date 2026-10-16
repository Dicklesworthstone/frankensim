# fs-fmm CONTRACT

## Purpose and layer

Layer: L2. Kernel-independent black-box FMM (plan §8.3 [F], bead
tfz.20): the Fong–Darve Chebyshev scheme — every translation operator
is polynomial interpolation, so any kernel smooth away from the
diagonal works unchanged, and accuracy is controlled by one number
(the interpolation order p).

## Public types and semantics

- `Kernel` (trait): `eval(x, y)` — smooth off the diagonal; the
  diagonal convention is the caller's (Laplace3d returns 0 there).
- `Laplace3d`: `1/(4π|x−y|)`.
- `Fmm::new(kernel, points, order, leaf_cap)`: fallible UNIFORM-depth octree
  (depth from N/leaf_cap, empty cells omitted, ancestors registered;
  `leaf_cap` targets uniform-cloud occupancy rather than imposing a strict
  clustered-leaf capacity) —
  on a uniform tree, "adjacent leaves run P2P, first-separated
  ancestors run M2L" partitions every source–target pair EXACTLY ONCE
  (no gap, no double count). `potentials(charges)` runs P2M/M2M
  (anterpolation to p³ Chebyshev grids), M2L (kernel evaluated between
  well-separated same-level grids whose parents are adjacent), L2L/L2P
  and the direct near field. `direct(charges)` is the O(N²) oracle;
  `stats()` the octree ledger row.

## Invariants

1. Accuracy is controlled by p: the order sweep against the direct
   oracle falls monotonically, reaching < 1e-5 relative L2 at p = 7 on
   uniform clouds (fmm-001; curves ledgered).
2. G3 translation invariance: a rigidly shifted cloud reproduces
   potentials to < 1e-9 relative (fmm-002).
3. Scaling: time-vs-N fitted exponent < 1.6 over a 4096→32768
   doubling ladder (fmm-003; O(N log N)-class trend — the 10⁷-point
   wall-clock target is the perf lanes' scope, stated, not silently
   skipped).
4. Determinism: BTree-keyed tree, fixed traversal orders.

## Error model

`FmmError` rejects empty/non-finite/oversized clouds, invalid order or leaf
capacity, charge mismatches/non-finite charges, non-finite kernel output, and
conservative coefficient-slot, cell-pair scan, M2L translation, direct-oracle,
and distribution-aware near-field pair work envelopes before the corresponding
large allocation or pass.
No silent accuracy degradation: the order is explicit at construction and the
battery curves are the evidence. These are admission guarantees, not a claim
that process-level allocator exhaustion is recoverable: `BTreeMap` and several
internal pass buffers still use Rust's infallible allocation APIs after the
request has passed its bounded envelope.

## Determinism class

Bit-deterministic across runs on a platform. Cross-ISA goldens not
yet recorded.

## Cancellation behavior

The current passes are synchronous and do **not** accept a `Cx`, poll
cancellation, or serialize mid-pass resume state. The explicit work envelope
bounds admitted calls, but is not a cancellation-latency claim. Cross-crate
Cx/resume integration is tracked separately under `frankensim-ccmn`.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None ([F] realized as a standalone crate per the crate-granular
gating rule).

## Conformance tests

`tests/battery.rs`: fmm-001 order sweep vs oracle; fmm-002
translation invariance; fmm-003 scaling trend; invalid geometry, tuning,
charge, and work-envelope refusal.

Completed battery cases emit schema-validated `fs-obs` JSONL under suite
`fs-fmm/battery`. Aggregate identities are `fmm-001` through `fmm-003`, with
distinct `fmm-00N/measurement` companions that retain the order/error/tree,
translation-error, and timing/exponent rows. Their literal input roots are
`0x1001_2026_0708_0001`, `0x1001_2026_0708_0002`, and
`0x1001_2026_0708_0003`; the last root is intentionally replayed independently
for every size in the scaling ladder. Measurements and aggregate verdicts are
emitted before the original terminal assertion, so a failed gate remains an
ordinary Rust test failure as well as a canonical error event.

## No-claim boundaries

- Adaptive trees (U/V/W/X interaction lists) — the uniform-depth
  partition ships; adaptivity is the recorded successor.
- SIMD-batched near field and precomputed M2L tables (perf lanes).
- Gradient outputs (consumers run per-component kernels today;
  a fused gradient pass is follow-up).
- Periodic/boundary-image variants; high-kD directional oscillatory kernels.
- Recoverable allocator exhaustion for internal tree/map/pass buffers. Work is
  bounded before those allocations, but allocator failure remains process-level.
- The `fmm-003/measurement` timing rows are diagnostic observations, not
  portable performance authority; performance claims still require the
  machine-fingerprinted perf lanes and their acceptance bands.
- Retired from v1 (owner decision 2026-09-01, bead frankensim-rc-root-q61wp.25):
  the adaptive kernel-independent (bbFMM) tree that would drive an FMM-accelerated
  vortex particle method is NOT scheduled; the promise is preserved in deferred
  epic frankensim-rc-root-q61wp.43 and un-parks only behind a consumer journey.

## Complex Helmholtz extension (IRP M2.5-W5a)

`ComplexKernel::eval -> fs_math::c64::C64` and `FmmC64` provide complex charges
and potentials. The existing real `Kernel`, `Laplace3d` and `Fmm` arithmetic is
unchanged. The complex engine reuses the existing uniform-depth tree builder
at order 2, which allocates no expansion coefficients, then applies its own
order-dependent admission and complex passes. It uses the same exactly-once
near/far pair partition. `direct` remains a separate, bounded quadratic oracle;
`potentials` never falls back to it.

`FmmC64::near_field_blocks()` borrows each nonempty leaf's original point
indices in deterministic tree order. The blocks partition all input indices
exactly once. Callers can split these lists to build bounded local
preconditioners; the method adds no kernel evaluation or change to the pair
partition. The `near_blocks_partition_indices_in_fixed_order` unit test checks
this partition directly; BEM exercises it through block-Jacobi solves.

The outgoing convention is `exp(-i omega t)`, `G=exp(ikr)/(4 pi r)`, matching
`fs-bem::helmholtz`. Shipped kernels are `HelmholtzSingle`, `HelmholtzDouble`
(`dG/dn_y`), `HelmholtzAdjointDouble` (`dG/dn_x`), `HelmholtzGradient` (one
Cartesian `grad_x G` component), and `HelmholtzHypersingular` (the mixed normal
derivative **away from the diagonal**). Self terms, jumps, finite parts and
singular/near-panel quadrature remain BEM's responsibility. These point kernels
alone are not a complete Burton–Miller panel operator.

Double/adjoint/mixed constructors accept one unit normal per point, squared
norm tolerance `1e-12`. Source normals weight three P2M charge fields; target
normals contract three L2P fields. Normals are not interpolated or evaluated at
synthetic grid points. No automatic normalization changes caller geometry.
All new kernel and interpolation transcendentals use `fs_math::det`; inherited
tree-depth selection still uses the real builder's platform `ln`. The existing
real interpolation and Laplace `powi` remain untouched. The claim is repeated
execution on the same platform/profile, without a new cross-ISA claim.

The requested interpolation-order ladder is `4,6,8,10,12` against the complex
direct oracle for root-box `kD=1,5,10` (D includes the builder's `1e-9` padding).
`helmholtz_order_for` conservatively recommends p=12 only within `kD<=10` and
`1e-6<=tol<1`. Higher kD or tighter requested tolerances return the named
`FmmError::HelmholtzRegime`, route `M2.5-W5b`. Actual Helmholtz constructors
also enforce the kD ceiling. Generic custom kernels must report their
Helmholtz wavenumber through `ComplexKernel::wavenumber` when applicable.
Explicit lower p values are diagnostic approximations. The order selector is
an empirically tested domain recommendation, **not** a universal relative-error
certificate for arbitrary charge cancellation or untested point distributions.
Consumers needing such a certificate must independently bound their operator
error; this extension does not mint that authority. The new numeric ladders
remain pending central native verification until the bead's run receipt exists.

`MAX_COMPLEX_COEFFICIENT_SLOTS=MAX_COEFFICIENT_SLOTS/2` keeps the original byte
budget. Admission counts source and target channels, every actual tree cell,
and two extra cells for transient coefficient buffers. M2L work counts actual
interaction pairs times `p^6` times source/target channel counts, with the
existing `MAX_TRANSLATION_WORK` ceiling. Tree/cell-pair/point/near-field caps
remain in force. This avoids the real path's overly conservative M2L estimate
rejecting p=12 even on a small sparse tree; it does not raise any cap.
Allocation failure and cancellation boundaries remain as described above.

`tests/battery.rs` includes `battery/helmholtz.rs` in the existing test target.
It emits fs-obs JSONL rows for the single/double order ladders, adjoint/gradient/
mixed checks, `1e-12` translation, high-kD refusal (including an unguarded
test-only p12 failure), coefficient cap/cap+1 and bitwise replay. Unit tests
independently finite-difference the single-layer kernel to check all derivative
signs, check analytic outgoing phase, and check pair partition and input errors.
`laplace_bitwise_unchanged` freezes the original real operation-tree source at
`ee8e613806eaa9ac3979a8b139af5e1528b699d5`; it is a source guard, not a fabricated
historical runtime golden. Existing fmm-001..003 stay intact. Native comparison
against a baseline run is still required to claim historic measurement identity;
elapsed timing rows are diagnostic and cannot be bitwise reproducibility claims.
