# CONTRACT: fs-rand

## Purpose and layer
Counter-based Philox streams keyed by LOGICAL work identity + deterministic
distributions (plan §6.7; P2's seed pillar). Layer: L1.

## Public types and semantics
- `philox::philox4x32_10(ctr, key)` — the Random123 block function, KAT-pinned.
- `StreamKey { seed: u64, kernel: u32, tile: u32 }` — the Cx-carried logical
  identity; field widths are contract (2⁶⁴ draws per stream, 2³² kernels/tiles).
- `Stream` — sequential view with RANDOM ACCESS (`Stream::at(key, index)`),
  `Copy` (forks diverge by IDENTITY, not state). Persisted replay uses
  `StreamCheckpoint { checkpoint_version, stream_semantics_version, key,
  index }`; `Stream::resume` accepts only both exact current versions.
- Canonical retained checkpoints are exactly 83 bytes from
  `StreamCheckpoint::to_canonical_le_bytes`: the 8-byte `FSRCKPT\0` magic,
  the exact 43-byte `org.frankensim.fs-rand.stream-checkpoint.v1` domain, and
  little-endian checkpoint version, stream-semantics version, seed, kernel,
  tile, and next index. `StreamCheckpoint::from_canonical_le_bytes` and
  `Stream::resume_retained` refuse truncation, trailing bytes, foreign magic or
  domain, and past or future versions before replay.
- Draws: `next_u64`, `next_f64` (53-bit, [0,1)), `next_below` (Lemire,
  deterministic rejection consumption), `next_normal` (Box–Muller on
  fs-math strict fns — cross-ISA deterministic SAMPLES), `next_normal_ziggurat`
  (bead 1za9: the ZIGGURAT perf path — deterministic `det`-generated table +
  deterministic rejection consumption; FAST-MODE-ONLY until a cross-ISA bitwise
  proof admits it to strict mode), `next_exponential` (inversion). BULK
  generation (bead 1za9): `fill_f64` / `fill_u64` fill via 8-lane batched
  `philox4x32_10_batch` (structure-of-arrays, auto-vectorizable) and are
  BITWISE-IDENTICAL to the sequential draws (index advances by `len`).

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
- Beta and Dirichlet use direct gamma normalization when its total is finite
  and positive. If valid extreme shapes underflow every gamma to zero (or make
  the direct total overflow), they replay the same checkpointed samples as
  scaled log weights. The replay advances no caller-visible stream state, uses
  only strict arithmetic, and guarantees finite normalized outputs without
  changing ordinary-shape bits.
- vMF mean-resultant lengths match the analytic coth(κ) − 1/κ at
  κ ∈ {1, 10, 100} (tested); truncated-normal mean matches the analytic
  hazard ratio via erfc (tested).

## Invariants
- A draw is a pure function of (seed, kernel, tile, index) — never of
  thread/worker/order (shuffle-invariance is a test).
- Stream counters advance modulo 2⁶⁴, matching the Philox counter width, so
  debug and release builds agree even at the end of the counter space. Callers
  needing more than 2⁶⁴ draws must change `StreamKey` identity rather than
  relying on a longer single stream.
- Random access ≡ sequential access (tested bitwise).
- Rejection sampling advances the index deterministically (replay-safe;
  consumed-count is content-determined — tested).
- A retained checkpoint cannot bypass version admission: stale/future
  checkpoint transports and stale/future stream semantics return structured
  `StreamReplayError` before any draw.
- Canonical transport binds every replay field at one documented little-endian
  offset. Independent mutations prove that checkpoint version, stream version,
  seed, kernel, tile, and index move only their own fixed-width byte range;
  magic/domain mutations are refused as cross-type/cross-domain input.
- Integer core is trivially cross-ISA; float distributions inherit fs-math's
  proven cross-ISA determinism.

## Error model
Invalid distribution parameters panic as programmer errors (`next_below(0)`,
non-positive/non-finite gamma-family shapes, empty or length-mismatched
Dirichlet outputs, and invalid truncation/vMF parameters). Within those
documented domains, operations are total; in particular, every valid beta is
finite in `[0,1]` and every valid Dirichlet result is a finite non-negative
simplex point even at `f64::MIN_POSITIVE` shapes.
Untrusted retained replay state returns `StreamReplayError` for a non-canonical
length (including trailing data), foreign magic/domain, or an unknown past or
future checkpoint/stream-semantics version; it is never interpreted under the
current Philox mapping by guesswork.

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
rejection-replay, versioned checkpoint-resume equality, stale-version refusal,
the exact 83-byte canonical checkpoint KAT, truncation/trailing refusal,
independent transport-field mutations, the nonzero Random123 `Stream::at`
mapping KAT, and independent low/high seed/index plus kernel/tile draw-identity
mutations. Bead 1za9: `tests/ziggurat.rs`
(ziggurat moments + bit-determinism + two-sample KS vs Box–Muller) and
`tests/stream_battery.rs` (DEV-ONLY: uniform χ², lag-1 serial correlation,
monobit balance, the fixed 8-stream × 100,000 inter-stream correlation
matrix, and bead 4nh8's 128-case shrink-armed logical-identity battery).
The generated battery checks distinct blocks at counters 0, 1, and 2⁶⁴−1
plus a 4,096-sample `|correlation| < 0.10` smoke band; seed
`0xF5_AA_0001` is the replay root.

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
- The canonical checkpoint frame is deterministic type/domain separation, not
  a cryptographic authenticator. Retained storage must provide its own content
  digest or authenticated envelope when corruption or hostile modification is
  in scope; a changed but well-formed seed/kernel/tile/index is a different
  valid replay identity.
- Sobol dims > 10 (full Joe-Kuo table import = recorded follow-up).
- Owen scrambling performance (correct lazy-tree v1 is 32 Philox calls per
  point-dim; hash-based fast path = recorded follow-up).
- Gamma/beta/Dirichlet/categorical-alias/von-Mises–Fisher/truncated
  distributions: follow-up bead (consumer-driven: UQ/BO/rendering).
- Ziggurat normal STRICT-MODE admission (bead 1za9): `next_normal_ziggurat`
  ships FAST-MODE-ONLY — deterministic table + rejection consumption and
  accuracy-gated (moments + two-sample KS vs Box–Muller), but not admitted to
  strict mode until a cross-ISA bitwise-equal run on both reference machines
  lands (the `trj` pipeline); Box–Muller stays the strict default.
- SIMD bulk generation (bead 1za9): the batched SoA primitive
  `philox4x32_10_batch` + `fill_f64`/`fill_u64` ship and are gated
  bitwise-equivalent to the scalar stream (`tests/bulk.rs`); a HAND-WRITTEN
  NEON/AVX u32 kernel and its throughput gate remain a follow-up — `fs-simd`
  exposes only float ops and the workspace is stable-Rust (no `portable_simd`),
  so the vectorized capsule is a separate hardware-gated effort.
- A FULL PractRand/TestU01 port: `tests/stream_battery.rs` ships a
  representative dev-only subset (χ²/serial/monobit/inter-stream); the complete
  suite + CI wiring remain a follow-up. Finite sampled correlation bands are
  defect detectors, not proofs of statistical independence.

## Exec key bridge (bead wf9.7.1, v1)

`StreamKey::from_exec_parts(seed, kernel_id, tile, iteration)` is the
CHECKED bridge from fs-exec's four-u64 logical key: seed is lossless,
kernel and tile must fit their u32 slots, and iteration must be 0 (the
draw index is a within-stream counter, not an identity axis — callers
with generation-diverging streams ledger the generation upstream).
REFUSAL, never truncation: a lossy mapping would alias distinct
logical streams (collision-tested at every boundary). Field widths and
rules are versioned (`EXEC_KEY_BRIDGE_VERSION`); bump only with
recorded justification, since ledgered replay depends on them.

## Perf-lane evidence (bead 1za9, measured)

- Release, macos-aarch64 (Mac16,11): ziggurat 84.2M normals/s vs
  Box–Muller 29.4M/s — 2.86× — the perf path justified by measurement.
- Bulk SoA fill measured ~0.9× the scalar loop: bitwise-equivalent but
  NOT yet faster — the speedup claim awaits the hand-written NEON/AVX
  Philox capsule (the recorded resource-gated no-claim); the perf-lane
  gate trips only on pathological regression.
- The cross-ISA GOLDEN HASH (995960fe709f00bc over 100k draws) is the
  ready-to-run strict-mode admission instrument: it reproduces across
  independent aarch64 environments; the x86-64 run completes the trj
  proof the moment an x86 reference machine is available (the rch
  fleet is currently ARM-only, verified by census).
