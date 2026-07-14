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
  `ResolutionTooFine` with the 512-samples/axis cap). Unbounded charts use
  `build_clipped(source, target_h, clip, cx)`, which samples the geometric
  intersection `source ∩ clip`, not merely a replacement extent. f32 storage on
  `fs_substrate::TiledField`, f64 evaluation via triquadratic B-spline
  (C¹: continuous analytic gradients — shape optimization differentiates
  through samples). `bound()` is retained as an explicitly NOMINAL
  compatibility accessor: like `nominal_field_bound()`, it bounds
  reconstruction relative to the sampled chart field by the source Lipschitz
  constant times an outward Euclidean radius to the farthest actual
  spline-support node, plus f32 storage and interpolation-roundoff terms. The
  storage term is the outward-rounded maximum of
  `max_abs_sample·f32::EPSILON` and half the minimum positive f32 subnormal, so
  subnormal spacing and underflow-to-zero cannot escape the band. Dense axes
  store every actual representable sample coordinate and require them to be
  strictly increasing. Evaluation locates those actual cells and uses a
  monotone C¹ physical-to-cardinal-index warp; a translated domain with too few
  floats for the requested nodes receives `DenseLatticeUnrepresentable`
  before source evaluation. The spline radius uses the maximum actual gap on
  every axis, so nonuniform one-ulp placement cannot escape the theorem.
  Nonnegative products, sums, quotients, norms, source-certificate radii,
  interpolation error, and abstract-bound composition are rounded outward.
  The nominal bound is not by itself authority relative to the abstract
  region's signed distance.
  `abstract_distance_kind()` and
  `abstract_distance_bound()` aggregate every sampled `ChartSample::error`:
  Exact/Enclosure remain rigorous only when the source also opts into the
  global `TraceStepClaim::ExactDistance` theorem; interpolation demotes Exact
  to Enclosure. `ExactDistance` imposes a global Lipschitz floor of 1, so a
  contradictory local `Some(L < 1)` hint cannot shrink the rigorous
  reconstruction band or derived chart bound. Sampled LOCAL Lipschitz maxima
  without that theorem make the result at best Estimate and are not re-exposed
  as a certified Lipschitz value. Estimate remains weak, NoClaim absorbs,
  malformed finite certificates fail closed to NoClaim, and the maximum finite
  source certificate radius is added to the nominal reconstruction bound. This
  authority discipline also applies outside support: an `ExactDistance`
  source propagates from the clamped point with its theorem constant `L=1`,
  independent of a looser sampled `Some(L)` hint. The outside norm and interval
  arithmetic are outward, the published nominal value is included in the
  certificate hull, and non-finite queries or overflowed arithmetic fail
  closed to NoClaim. Chart Lipschitz claim =
  `L + 2·nominal_field_bound/h_min` (sample-slope bound, derivation in code).
  Both the nominal bound and this derived Lipschitz value are preflighted as
  finite before publication. Axis construction uses convex span fractions
  (with exact endpoints), not `min + index·h`, then retains the resulting
  strictly ordered nodes as representation state. `downgrade_abstract_distance_authority(kind)`
  is a monotone weakening-only hook: strengthening requests are ignored and
  NoClaim clears only the abstract bound, retaining the finite nominal field.
  `measure_eikonal` (seeded |∇φ|−1 statistics: EVIDENCE, labeled;
  polls every 256 probes and returns structured cancellation),
  `resample_box` (finite dirty boxes are validated and intersected with the
  field before index arithmetic; incremental re-sampling occurs at EXACTLY the
  original positions and is staged transactionally before field mutation —
  bit-identical to a full rebuild over the same source; source authority,
  maximum certificate radius, Lipschitz, quantization slack, and the field
  commit together, while cancellation/refusal leaves all prior state
  unchanged; the converter beads' G5 law),
  `mean_curvature_estimate` (Estimate-grade stencil), `raycast` —
  cancellable sphere tracing with steps
  `(sd − abstract_distance_bound)/lipschitz` only for rigorous
  abstract-distance fields. A ray is intersected with the stored finite AABB
  before the reconstruction is sampled, so an outside origin cannot silently
  reuse a clamped boundary value; non-finite rays and overflowed slab/step
  arithmetic fail closed. Its compatibility `Option` API returns no hit for
  Estimate/NoClaim fields rather than laundering the nominal bound.
- `VdbGrid<T: Copy>` — FrankenVDB: BTreeMap root → 32³ internal nodes
  (sparse BTreeMap children) → 8³ bitmasked leaves with full value
  arrays. `set/get/is_active/deactivate/active_count`, DETERMINISTIC
  `iter_active` (root order → internal order → leaf linear order; no
  HashMap anywhere near results), face-6 `dilate` (new voxels copy the
  activating value)/`erode`, `memory_stats` (layout-derived footprint +
  bytes-per-active, ledgered).
- `NarrowBand` — band of `half_width_cells` around a chart's zero set on
  the VDB: cancellable `from_chart` / `from_chart_clipped`, with a
  16,777,216-point deterministic pre-sparsification scan cap and checked signed
  VDB coordinate range. Its exactly uniform lattice uses a floor-derived count
  and max-anchored origin, so every node stays within the admitted inflated
  support and the final node hits the finite maximum without endpoint
  overshoot, including near `f64::MAX`. Every source result must be finite as
  both `f64` and stored `f32` before it enters the sparse grid. Trilinear
  `interpolate` returns None outside the band; semi-Lagrangian `advect` is
  first-order with band dilate+trim per step; Godunov-upwind `reinitialize`
  uses the smoothed sign function; `stats` records active count, mean eikonal
  residual over interior voxels, and max |φ|.
- `AdaptiveSdf` — octree with per-cell trilinear corner fits; refinement
  splits cells whose residual at 7 probe points exceeds tolerance (down
  to max depth); `build_clipped` samples `source ∩ clip`; the worst-case full
  octree is checked against a 1,000,000-node cap before evaluation or
  allocation. Before any branch is created, all three `f64` midpoints must be
  finite and strictly interior; adjacent-float axes receive a structured
  `AdaptiveSubdivisionUnrepresentable` refusal instead of degenerate children.
  interpolation is overflow-stable for opposite-sign finite corners, and any
  non-finite fit or probe residual is a structured reconstruction refusal;
  `stats` ledger the cells/depth/worst residual. Its probed fit is at best
  Estimate-grade, composed with the maximum valid source-certificate radius;
  source NoClaim or malformed certificates remain NoClaim inside and outside
  support. The nominal and abstract bands are separately exposed.

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
7. Sampling admission (rsdf-007): default dense/adaptive/band builders reject
   unresolved extended supports before chart evaluation; paired clipped APIs
   sample the actual geometric intersection; invalid spacing and checked
   count/work-limit refusals also precede evaluation; translated source+clip
   pairs give translated fields (G3). Narrow-band construction rejects
   non-finite or non-`f32`-representable samples and its max-anchored lattice
   remains finite and within an extreme admitted support near `f64::MAX`;
   adaptive refinement refuses unrepresentable adjacent-float subdivisions.
8. Authority preservation (rsdf-008): dense and adaptive reconstruction add
   the maximum finite source-certificate radius to their nominal fit band,
   retain Estimate as weak authority, refuse to promote sampled local
   Lipschitz maxima into a global enclosure, and make explicit or malformed
   NoClaim absorbing inside and outside support; an ExactDistance source cannot
   understate global L below 1; a tiny translated ExactDistance fixture proves
   f32 underflow error remains inside the advertised enclosure; translated
   few-ulp domains either retain strict actual nodes with a radius derived from
   their nonuniform gaps or refuse a collapsed lattice; outside propagation
   uses theorem `L=1`, contains its nominal value, and fails closed on
   non-finite arithmetic; raycast refuses weak authority.

## Error model
Structured teaching errors (`SdfBuildError`) wrap `SamplingDomainError` and
name invalid spacing, count/coordinate overflow, deterministic work caps, and
cancellation, plus non-finite/unrepresentable samples and non-finite
constructed bounds. Unrepresentable dense axes and adaptive midpoints are
structured refusals, not coincident nodes or degenerate child cells.
Estimate-vs-Enclosure-vs-NoClaim honesty is carried in every `ChartSample`;
invalid derived Lipschitz claims are structured refusals. No panics across the
boundary.

## Determinism class
Deterministic: seeded probes, BTreeMap orders everywhere, no clocks or
addresses in results. Float behavior inherits scalar-arithmetic classes.

## Cancellation behavior
`TiledSdf::build` polls directly per sampled row and through its eikonal
probes; `TiledSdf::measure_eikonal` polls every 256 probes and before
publication; `TiledSdf::resample_box` polls per staged row and once before its
mutation-only commit; `NarrowBand::from_chart` polls at most every 256 source
evaluations and once more before publishing the completed grid;
`AdaptiveSdf::build` polls per cell; `raycast` polls per step. Bounded work
between polls (P7).

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
None. `[S]` solid-tier throughout.

## Conformance tests
tests/conformance.rs, cases rsdf-001..rsdf-008 (JSON-line verdicts;
seeded cases carry seeds) covering invariants 1–8 with fs-obs-validated
evidence events (eikonal stats, VDB footprint, band stats, adaptive
residuals).

## No-claim boundaries
- Eikonal deviation and adaptive residuals are MEASURED statistics
  (Estimate-grade), not certificates; fs-ivl interval-verified sampling
  promotes them later.
- A finite nominal-field reconstruction bound remains useful for clipped and
  other NoClaim sources, but confers no abstract-region distance or safe-ray
  authority. The separate accessors and ChartSample kind enforce that line.
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
