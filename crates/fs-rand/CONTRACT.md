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

## No-claim boundaries
- Sobol/Owen/lattice QMC: the separate fs-rand-qmc bead.
- Gamma/beta/Dirichlet/categorical-alias/von-Mises–Fisher/truncated
  distributions: follow-up bead (consumer-driven: UQ/BO/rendering).
- Ziggurat normal (perf) — Box–Muller chosen v1 FOR strict-mode determinism.
- SIMD bulk generation lanes; PractRand/TestU01-class nightly battery.
