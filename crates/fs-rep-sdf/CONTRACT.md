# CONTRACT: fs-rep-sdf

## Purpose and layer
Signed-distance-field charts (plan §7.2): dense Morton-tiled grids with C¹
reconstruction, the FrankenVDB sparse hierarchical tile tree, adaptive
octree SDFs, and narrow-band level sets — all presenting through fs-geom's
`Chart` so agreement checking and conversion receipts apply unchanged.
Layer: L2. Depends on fs-geom, fs-substrate, fs-exec, fs-evidence,
fs-alloc, fs-obs.

## Public types and semantics
- `TiledSdf` — built from ANY chart with a certified Lipschitz bound
  (`build(source, target_h, cx)`, teaching refusals: `NoLipschitzBound`,
  `ResolutionTooFine` with the 512-samples/axis cap). f32 storage on
  `fs_substrate::TiledField`, f64 evaluation via triquadratic B-spline
  (C¹: continuous analytic gradients — shape optimization differentiates
  through samples). Declared error enclosure = `2·L·h_max + f32-quant`,
  conservative by construction. Chart Lipschitz claim =
  `L + 2·bound/h_min` (sample-slope bound, derivation in code).
  `measure_eikonal` (seeded |∇φ|−1 statistics: EVIDENCE, labeled),
  `mean_curvature_estimate` (Estimate-grade stencil), `raycast` —
  sphere tracing with steps `(sd − bound)/lipschitz`, cancellable.
- `VdbGrid<T: Copy>` — FrankenVDB: BTreeMap root → 32³ internal nodes
  (sparse BTreeMap children) → 8³ bitmasked leaves with full value
  arrays. `set/get/is_active/deactivate/active_count`, DETERMINISTIC
  `iter_active` (root order → internal order → leaf linear order; no
  HashMap anywhere near results), face-6 `dilate` (new voxels copy the
  activating value)/`erode`, `memory_stats` (layout-derived footprint +
  bytes-per-active, ledgered).
- `NarrowBand` — band of `half_width_cells` around a chart's zero set on
  the VDB: cancellable `from_chart`, trilinear `interpolate` (None
  outside band), semi-Lagrangian `advect` (first-order, band
  dilate+trim per step), Godunov-upwind `reinitialize` sweeps with the
  smoothed sign function, `stats` (active count, mean eikonal residual
  over interior voxels, max |φ|).
- `AdaptiveSdf` — octree with per-cell trilinear corner fits; refinement
  splits cells whose residual at 7 probe points exceeds tolerance (down
  to max depth); `stats` ledger the cells/depth/worst residual; error
  model is ESTIMATE-grade (probed, not enclosed) and says so.

## Invariants
1. Fixture reproduction (rsdf-001): |TiledSdf − source| ≤ the declared
   enclosure over 12k seeded points on sphere/box/torus.
2. C¹ seamlessness (rsdf-002): values and gradients vary continuously
   across tile boundaries (the B-spline never sees tile seams — storage
   layout is invisible to reconstruction); gradients match central FD.
3. Eikonal honesty (rsdf-003): |∇φ|−1 statistics are measured, ledgered,
   and NEVER promoted to a certificate.
4. VDB oracle equivalence (rsdf-004): set/get/iterate/miss agree EXACTLY
   with a BTreeMap oracle on clustered+stray actives; iteration order is
   total and deterministic; dilate grows and erode shrinks the active
   set; footprint stats ledgered.
5. Narrow-band evolution (rsdf-005): translation drift of the zero
   crossing stays within the stated fixture bound; reinitialization does
   not worsen the mean eikonal residual; band stats ledgered.
6. Sphere tracing (rsdf-006): steps respect the chart's own bound and
   Lipschitz claim (never tunnel by construction); hits match analytic
   intersections within bound-scaled tolerance; misses miss; polls
   cancellation.

## Error model
Structured teaching errors (`SdfBuildError`), `fs_exec::Cancelled` for
interrupted builds/checks, Estimate-vs-Enclosure honesty carried in every
`ChartSample`. No panics across the boundary.

## Determinism class
Deterministic: seeded probes, BTreeMap orders everywhere, no clocks or
addresses in results. Float behavior inherits scalar-arithmetic classes.

## Cancellation behavior
`TiledSdf::build` polls through the source chart's eval; `NarrowBand::
from_chart` polls per row; `AdaptiveSdf::build` polls per cell; `raycast`
polls per step. Bounded work between polls (P7).

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
None. `[S]` solid-tier throughout.

## Conformance tests
tests/conformance.rs, cases rsdf-001..rsdf-006 (JSON-line verdicts;
seeded cases carry seeds) covering invariants 1–6 with fs-obs-validated
evidence events (eikonal stats, VDB footprint, band stats, adaptive
residuals).

## No-claim boundaries
- Eikonal deviation and adaptive residuals are MEASURED statistics
  (Estimate-grade), not certificates; fs-ivl interval-verified sampling
  promotes them later.
- The hash-grid variant is deferred to the point-cloud chart bead
  (wqd.6) where its consumer lives.
- Narrow-band advection is first-order semi-Lagrangian with fixture-
  grade drift bounds; WENO advection, fast-iterative redistancing, and
  velocity extension are the topo-levelset bead's (9.5), which owns the
  production evolution claims.
- NO performance claims: sample-throughput rooflines belong to the perf
  harness; VDB "O(1)-ish" access is verified structurally (two bounded
  map hops), not benchmarked here.
- The 10⁹-voxel footprint target is scaled in CI (~12.5k actives with
  the per-active overhead ledgered); the absolute-scale run is a
  nightly/perf-lane job.
- Curvature is a labeled estimate; certified stencils arrive with
  fs-ivl integration.
