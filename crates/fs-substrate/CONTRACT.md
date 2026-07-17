# CONTRACT: fs-substrate

## Purpose and layer
Hardware capability probes, stable machine fingerprints, topology facts,
one-shot SIMD dispatch resolution (plan §5.1, patch Rev Q: measured facts,
never static claims), and the Morton/tile-major layout layer — tile
identities, tiled field containers, halo gathers, and shard-by-CCD affinity
maps (plan §5.3). Layer: L0.

## Public types and semantics
Probe surface:
- `CACHE_LINE: usize` — compile-time cache-line const (128 Apple aarch64 /
  64 elsewhere), runtime-VERIFIED against the OS report by the test suite.
- `CapabilityProbe` — `run()` (topology + ~100–300 ms bandwidth measurement)
  / `topology_only()` (fast); fields: ISA, brand, sorted feature list,
  cache line, page size, memory, logical CPUs, Apple perf/eff core counts,
  per-cluster L2, Linux NUMA node count, OPTIONAL gpu fact (separate class,
  unused by the CPU build), `Measured { single_thread_gbs, all_core_gbs }`.
- `fingerprint() -> u64` — FNV-1a over STABLE facts only; measured bandwidth
  excluded by design (jitter must not mint new machine identities).
- `to_json()` — ledger `capability_probes` row shape; caller supplies
  probe_time (this crate is clock-free outside the bandwidth measurement).
- `dispatch_tier() -> SimdTier` (Scalar/Neon/Avx2/Avx512) — resolved EXACTLY
  once (OnceLock; `dispatch_resolution_count()` test hook proves it).
- `bandwidth::measure` — STREAM-triad sweep (24 B/elem accounting, 64 MiB
  arrays, safe 16-element unroll plus exact scalar tail), best-of-5
  single-thread + best-of-3 all-core aggregate. The fixed-width body exposes
  four AVX2 vectors or eight Neon vectors of independent work to the optimizer
  without ISA-specific code or a different accounting convention.

Tile-layout surface:
- `morton::{morton3_encode, morton3_decode}` — 3D Morton codec, 21 bits per
  axis; backends (magic-bits reference; BMI2 PDEP/PEXT capsule on x86-64)
  resolved once into a function table, `morton_backend()` reports which.
- `tile::TileEdge` (4/8/16, default 8), `tile::TileCoord`,
  `tile::TileGrid` — ceil-division tile cover of a cell domain with
  structured `TileError` refusals (zero dims, Morton-domain overflow);
  deterministic iteration orders (`iter_linear`, `iter_zorder`,
  `iter_boundary_first`), cell↔tile maps, `zorder_ranks` storage slots.
- `tile::TileId { grid, code }` — THE stable logical tile identity
  (grid-geometry hash + Morton code); `stream_key() -> u128` keys Philox
  streams, reduction slots, and ledger events (Decalogue P2 foundation).
- `field::TiledField<T>` — dense tile-major storage in z-order slot order,
  x-fastest rows within tiles; per-tile `TileMeta` (first-touch owner tag,
  dirty flag, occupancy hint); `Boundary::{Clamp, Constant, Periodic}` halo
  policies; `gather_halo` (reference) and `gather_halo_fast` (row-copy body,
  per-cell ghost shell) — bit-identical by conformance law;
  `shard_views_mut` splits the field into per-shard disjoint `ShardViewMut`
  windows for first-touch initialization; `TiledField::<f32>::eval_f64` is
  the f32-storage/f64-evaluate SDF convention.
- `affinity::CcdTopology` — fixtures (TR 7995WX 12×8, EPYC 16×8, Apple M
  2-cluster) + `from_probe` heuristic (a recorded HINT, not a claim);
  `affinity::AffinityMap` — balanced (±1) contiguous z-order ranges per
  shard, `shard_of_slot`/`slots_of`/`to_json` (deterministic table).

## Invariants
- Fingerprint stable across probes on the same machine; changes iff stable
  topology changes; never depends on measured numbers.
- Dispatch resolution (SIMD tier, Morton backend) happens at most once per
  process; no hot-path re-detection.
- Probe never panics; missing OS facts degrade to None/defaults, loudly
  visible in the JSON row.
- Morton codec is bijective over `[0, 2^21)³` and backend-equivalent
  bit-for-bit (G0 law, 200k-case seeded battery + exhaustive small cube).
- Tile iteration orders are permutations of the tile set; z-order and
  boundary-first orders, tile ids, and affinity tables are deterministic
  (G5).
- `gather_halo_fast` ≡ `gather_halo` on every tile (including partial edge
  tiles) under every boundary policy.
- Affinity ranges cover every tile exactly once, balanced within one tile;
  per-shard core ranges never straddle a CCD boundary on the fixtures.
- Tile bases of `TiledField` storage are 128-byte aligned whenever
  `edge³ · size_of::<T>()` is a multiple of 128 (all f32/f64 fields).

## Error model
Absent OS facts → `None` fields (no fabrication). Subprocess/file read
failures degrade gracefully. Tile-grid construction returns structured
teaching `TileError`s. In-domain cell access is a programmer-error contract
(asserts before arithmetic, like fs-simd's length contract); boundary-aware
reads go through `get_bc`/halo APIs instead.

## Determinism class
Topology facts: deterministic per machine/boot. Bandwidth: measurement
(jittery by nature, quarantined from identity). Tile layer: fully
deterministic — orders, identities, affinity tables, and halo contents are
bit-stable across runs and thread counts (parallel first-touch writes land
in disjoint slots; conformance sub-006 checks parallel ≡ serial).

## Cancellation behavior
Probe: startup-time, latency-lane work; longest operation ~300 ms; no Cx.
Tile layer: bounded pure computation per call (a halo gather touches one
tile + ghost shell); tile-granular cancellation polling is the consuming
kernel's contract (fs-exec), for which this crate supplies the TileId
boundaries.

## Unsafe boundary
One registered capsule: `src/morton/bmi2/mod.rs` — BMI2 PDEP/PEXT bit
interleave (register-only intrinsics, no memory access), behind runtime
`is_x86_feature_detected!` dispatch with the magic-bits twin as reference
and Miri fallback; SAFETY.md beside it. Everything else: feature detection
via std::arch macros; OS facts via subprocess/file I/O (no FFI, per P1).

## Feature flags
None. SME/AMX-class matrix units: deliberately NOT claimed; future probes
gate the fs-simd `frontier-sme2` tier.

## Conformance tests
tests/conformance.rs, cases sub-001..sub-006. Each case's aggregate verdict
uses the canonical fs-obs `ConformanceCase` schema; randomized cases carry
their input seed. Morton bijection + backend equivalence, tile/world map
bijection + iteration permutations, halo fast-vs-reference equivalence
across boundary policies and partial tiles, affinity fixtures (balance,
coverage, CCD-respecting core ranges, deterministic tables), the tiled-vs-
linear stencil smoke (bitwise agreement + timings emitted separately as
`benchmark_result` events — documentation, not claims), and parallel
first-touch ≡ serial with owner tags recorded. The first-touch case has
deterministic input and output but makes no exact thread-interleaving replay
claim. Assertions that abort before the aggregate verdict remain ordinary
Rust test diagnostics. In-module suites cover the probe battery, Morton
known answers, tile geometry edge cases, halo boundary semantics, and
128-byte tile-base alignment.

## No-claim boundaries
- Per-core-CLASS bandwidth (P vs E pinning) — needs QoS/affinity outside safe
  std; per-class COUNTS are reported, aggregate bandwidth measured; real
  pinning lands with fs-exec or an audited capsule.
- x86-64/Linux paths (BMI2 capsule, THP-adjacent behavior) are exercised on
  the Threadripper CI runner where available; AVX-512-class hosts remain
  unverified until that runner lands.
- GPU/Neural facts are descriptive strings only; no capability claims.
- `CcdTopology::from_probe` is a HEURISTIC hint (Apple cluster counts and
  x86 16-logical-CPUs-per-CCD guesses); real CCD/L3-domain discovery needs
  the x86 machine's cache-topology sysfs and lands with the autotuner. The
  fixtures, not the heuristic, are the tested contract.
- NO stencil performance claim: sub-005 timings are measurements written to
  events for the record; roofline verdicts (targets, bands, regressions)
  belong to the perf-harness bead.
- Sparse occupancy is a HOOK (`TileMeta.occupied`); FrankenVDB-style masked
  tiles are a MORPH-layer bead.
- Thread pinning, NUMA binding, and actual first-touch scheduling are
  fs-exec's contract; this crate only produces maps, views, and tags.
