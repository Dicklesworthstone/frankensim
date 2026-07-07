# CONTRACT: fs-rand

## Purpose and layer
Counter-based Philox streams keyed by LOGICAL work identity + deterministic
distributions (plan §6.7; P2's seed pillar). Layer: L1.

## Public types and semantics
- `philox::philox4x32_10(ctr, key)` — the Random123 block function, KAT-pinned.
- `StreamKey { seed: u64, kernel: u32, tile: u32 }` — the Cx-carried logical
  identity; field widths are contract (2⁶⁴ draws per stream, 2³² kernels/tiles).
- `Stream` — sequential view with RANDOM ACCESS (`Stream::at(key, index)`),
  checkpoint/resume by index, `Copy` (forks diverge by IDENTITY, not state).
- Draws: `next_u64`, `next_f64` (53-bit, [0,1)), `next_below` (Lemire,
  deterministic rejection consumption), `next_normal` (Box–Muller on
  fs-math strict fns — cross-ISA deterministic SAMPLES), `next_exponential`
  (inversion), `fill_f64`.

### Extended distributions (bead 6ys.19, module `dist`)
- `Stream::{next_gamma, next_beta, next_dirichlet, next_truncated_normal,
  next_truncated_exponential, next_vmf3}` + `dist::AliasTable`.
- CONSUMPTION CONTRACTS: rejection samplers (gamma via Marsaglia–Tsang,
  truncated normal via Robert) advance the index on every proposal —
  consumed count is a pure function of stream content (replay-tested,
  including mid-stream interleaving). Fixed-consumption samplers are
  documented and TESTED as such: truncated exponential 1 draw, vMF 2
  draws (Ulrich inversion — no rejection), alias sampling 1 draw.
- AliasTable construction is DETERMINISTIC (index-order worklists,
  P2 on setup): same weights, same table, bitwise.
- All arithmetic routes through fs-math strict kernels (incl. the wf9.14
  pow for the α < 1 gamma boost) — sampled VALUES are cross-ISA
  bit-deterministic, golden-hashed (`0x4224_6e28_56de_673c`, verified on
  both reference ISAs).
- vMF mean-resultant lengths match the analytic coth(κ) − 1/κ at
  κ ∈ {1, 10, 100} (tested); truncated-normal mean matches the analytic
  hazard ratio via erfc (tested).

## Invariants
- A draw is a pure function of (seed, kernel, tile, index) — never of
  thread/worker/order (shuffle-invariance is a test).
- Random access ≡ sequential access (tested bitwise).
- Rejection sampling advances the index deterministically (replay-safe;
  consumed-count is content-determined — tested).
- Integer core is trivially cross-ISA; float distributions inherit fs-math's
  proven cross-ISA determinism.

## Error model
`next_below(0)` panics (programmer error). Everything else total.

## Determinism class
Deterministic CROSS-ISA (integer core + fs-math-strict distributions).

## Cancellation behavior
Pure computation, O(1) per draw; no poll points needed.

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
Random123 KATs (3 vectors), avalanche battery, random-access≡sequential,
16-tile×3-order shuffle invariance, adjacent-identity decorrelation,
chi-square/moment gates (uniform/normal/exponential), Lemire bias +
rejection-replay, checkpoint-resume equality.

## QMC (qmc module)
- `Sobol::new(dim)` / `Sobol::scrambled(dim, seed)` — base-2 Sobol,
  embedded Joe-Kuo head (dims 1..=10, preconditions asserted at load),
  Gray-code + RANDOM ACCESS `point(n)`; TRUE Owen nested-uniform scrambling
  via Philox-derived lazy random tree (zero storage, seed-replayable,
  net-preserving — all tested). Verified: exact per-dim stratification
  (m 1..=8) and 2D elementary intervals; scrambled-Sobol RMSE 3.17e-6 vs MC
  1.53e-3 at n=4096/dim=5 on Genz product-peak (~480x).
- `Lattice::cbc(n, dim)` — rank-1 CBC in the gamma=1 Korobov space (B2
  kernel), `korobov_error_sq` diagnostic (verified decay 4.46e-4@257 ->
  5.26e-5@1031, beats naive vectors), `baker` periodization.

## No-claim boundaries
- Sobol dims > 10 (full Joe-Kuo table import = recorded follow-up).
- Owen scrambling performance (correct lazy-tree v1 is 32 Philox calls per
  point-dim; hash-based fast path = recorded follow-up).
- Gamma/beta/Dirichlet/categorical-alias/von-Mises–Fisher/truncated
  distributions: follow-up bead (consumer-driven: UQ/BO/rendering).
- Ziggurat normal (perf) — Box–Muller chosen v1 FOR strict-mode determinism.
- SIMD bulk generation lanes; PractRand/TestU01-class nightly battery.
