# CONTRACT: fs-archive

Quality-diversity archives: MAP-Elites / CVT illumination + novelty.

## Purpose and layer

Layer L4 (ASCENT). No dependencies — pure Rust (`BTreeMap` for deterministic
niche storage).

## Public types and semantics

- `Elite { solution, descriptor, fitness }` — the best solution in a niche;
  descriptors and fitness values are finite, and fitness is non-negative so
  QD-score remains monotone when new niches are filled.
- `MapElites::new(lo, hi, bins)` — a gridded behavior space; `cell_of`
  (discretize, clamped), `add` (insert iff empty niche or strict improvement,
  returns whether it became an elite), `capacity`, `num_elites`, `coverage`,
  `qd_score`, `elite_at`, `best`, `elites`.
- `CvtArchive::new(centroids)` — the same illumination over a Voronoi
  tessellation; `nearest_centroid`, `add`, `capacity`, `num_elites`, `coverage`,
  `qd_score`, `best`.
- `novelty(descriptor, others, k)` — mean distance to the `k` nearest
  neighbours (`+∞` for an empty set).

## Invariants

- MAP-Elites keeps exactly one elite per niche: a worse candidate in a filled
  niche is rejected, a strictly better one replaces it.
- Descriptor dimensions match the archive dimension exactly; mismatches are
  rejected instead of being silently truncated in distance or cell math.
- ILLUMINATION preserves DIVERSITY: a low-fitness solution in an otherwise-empty
  niche is retained (unlike single-objective optimization).
- `coverage` and `qd_score` are monotone non-decreasing under `add`.
- `cell_of` / `nearest_centroid` assign descriptors correctly (grid clamps
  out-of-range descriptors).
- `novelty` grows with distance from the archive.

## Error model

Constructors and descriptor-entry methods panic on malformed configuration or
inputs: dimension mismatch, zero bins, empty centroids, non-increasing or
non-finite bounds, non-finite descriptors, and negative or non-finite fitness.

## Determinism class

Fully deterministic: `BTreeMap` gives a fixed niche order; identical add
sequences yield identical archives.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/archive.rs` (8 cases): descriptors map to the right cells; MAP-Elites
keeps the best per niche; illumination preserves diversity; coverage + QD-score
are monotone; the CVT archive assigns to the nearest centroid; novelty rewards
distance; malformed dimensions and fitness are rejected; determinism.

`tests/conformance.rs` registers the load-bearing QD laws with the shared
`fs-casebook` harness in a fixed order. Its exact records exercise MAP-Elites
strict niche replacement plus coverage/QD monotonicity, CVT nearest-centroid
ties plus strict replacement, and novelty known answers. Canonical little-endian
input frames bind every operation, descriptor, fitness, expected acceptance,
centroid, query, and expected result to literal FNV-1a digests; each record
separately declares its exact tolerance. Failure details retain the reversible
inputs and exact result bits. A disclosed seeded acceptance-oracle corruption
proves both the typed report and `assert_green` merge gate turn red.

## No-claim boundaries

- v0 is the ARCHIVE data structure + its metrics + novelty scoring. The full QD
  ALGORITHM loop (variation operators, batch selection/emitters, CMA-ME-class
  emitters) and the descriptor computations from layer-native definitions
  (behavior descriptors sourced from geometry/physics/optimization layers) are
  the fuller deliverable, staged.
- The CVT centroids are supplied by the caller; generating them by k-means /
  Lloyd relaxation over a sampled descriptor space is a follow-on.
- Novelty is brute-force k-NN; a spatial index for large archives is a
  performance follow-on.
- Casebook records are deterministic diagnostics, not authenticated ledger
  receipts or a benchmark-performance claim. A local green suite does not
  establish dual-ISA study-scale G5; the central conformance lane must compare
  the same input digests and exact verdict records on both reference ISAs.
