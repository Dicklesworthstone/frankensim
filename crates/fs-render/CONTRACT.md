# CONTRACT: fs-render

Unbiased spectral path-tracing core: the verifiable Monte-Carlo foundations.

## Purpose and layer

Layer L5 (LUMEN). The Monte-Carlo core depends on deterministic `fs-math`;
the default chart backends consume lower-layer `fs-evidence`, `fs-exec`,
`fs-geom`, and `fs-rep-nurbs`. Optional differentiable, volume, and tracer
surfaces add their declared lower-layer dependencies (`fs-ad` for the
differentiable lift). Pure Rust throughout.

## Public types and semantics

- `radical_inverse(base, i)` / `halton(dim, i)` — deterministic low-discrepancy
  coordinates (an image is as replayable as a solve).
- `cosine_sample_hemisphere(u1, u2) -> (dir, pdf)` — cosine-weighted hemisphere
  sample (`pdf = cosθ/π`).
- `Lambertian { albedo }` — `brdf` (`ρ/π`); `furnace_radiance(incident,
  samples)` — the FURNACE Monte-Carlo estimate (exactly `albedo·incident`).
- `balance_heuristic` / `power_heuristic` — MIS weights; `mis_weight_sum(pf,
  pg)` — the weight-sum audit (nominally `1`).
- `mis_integrate_unit(f, n)` — an unbiased MIS estimate of `∫₀¹ f` combining
  uniform + linear-importance strategies.
- `hero_wavelengths(hero, count, min, max)` / `spectral_integral(spectrum, min,
  max, samples)` — hero-wavelength spectral integration.

- `charts` module (plan §10.2, beads qfx.2 + 8ll9; [S], default-on through
  `chart-backends`): render charts that opt into a typed trace theorem, WITHOUT
  conversion; other chart types remain explicit no-claim previews until they
  supply the theorem their error model needs.
  `sphere_trace` steps within the endpoint closest to zero of a rigorous trace
  evaluation enclosure: for `ExactDistance` this defaults to the sample's
  abstract-distance certificate; for `LipschitzImplicit` the chart must supply a
  separate `trace_value_enclosure` of its real implicit field, whose zero-nearest
  magnitude is divided outward by the chart's CERTIFIED `L`. The sign cannot
  flip inside either radius, so the marcher provably never tunnels (audited:
  `TraceAudit.worst_step_ratio`);
  over-relaxation uses the standard certified fallback (retreat when
  spheres fail to overlap). `ray_intersect_nurbs` is grid-seeded 3×3
  Newton on `S(u,v) − o − t·d` with the `[S_u, S_v, −d]` Jacobian.
  Certification requires the chart-level typed `TraceStepClaim`; a sample
  carrying `Some(L)` cannot upgrade the default `NoClaim`. Exact-distance hits
  use world-space distance tolerance. A Lipschitz-implicit chart authorizes a
  geometric `Hit` with either a rigorous singleton-zero field enclosure or a
  short sign bracket. At a strict-sign residual limit, the marcher may inspect
  one cancellation-aware witness no more than `2*eps` ahead in actual
  outward-rounded Euclidean distance. Rigorously opposite endpoint signs plus
  the claim's finite-segment continuity prove a boundary in that segment; only
  an evaluated midpoint within `eps` of both endpoints becomes the hit. The
  witness is evidence only and is never adopted as march state. A nonzero
  normalized residual `|f|/L` by itself certifies step safety, not proximity,
  and stops as `ResidualLimit` with no `Hit`. Pending over-relaxed
  endpoints are validated before either hit or miss acceptance. An
  outward-rounded working limit is classified at the caller's actual endpoint:
  hits require a validated endpoint sample, and misses require an
  epsilon-clear safe-ball bridge across any coordinate-rounding gap. Normalized
  working parameters size steps, but every chart and overlap evaluation uses
  the caller ray's own `Ray::at` arithmetic so a certificate is never returned
  for a numerically different point. This is chart-backend bit-semantics v5.
  `TraceAudit`
  states whether every marched sample supplied a positive finite certified
  bound and compatible rigorous trace-value certificate, counts retreats to the
  safe endpoint, and
  distinguishes geometric hit, residual limit, clean miss, step-limit,
  invalid-input, and invalid-sample termination. Returned chart gradients are normalized before becoming hit
  normals. `trace_scene` and the spectral/differentiable renderers accept chart
  terminal results only when the full trace stayed certified; an uncertified
  miss is not evidence of empty geometry. The uncertified `L = 1` fallback is
  a direct-call preview surface, never a production geometry decision.
  Mixed-scene tracing returns `Result`, propagates cancellation and chart
  refusal, and enforces `t_max` uniformly across charts, NURBS, and meshes.
  `TriMesh` is Möller–Trumbore over a deterministic median-split BVH with
  outward-rounded slab pruning;
  `bvh_fingerprint` is a stable diagnostic receipt over its sorted layout.
  `trace_scene` mixes all three backend kinds by closest hit.

- `volumes` module (bead qfx.3, feature `volumes`): [`VolumeGrid`]
  BORROWS its density buffer (zero-copy: live simulation fields render
  in place), [`MajorantGrid`] per-block maxima, Woodcock delta
  tracking (`woodcock_transmittance`, unbiased for ANY bound ≥ max σ;
  the tile stage thins field lookups), the collision emission
  estimator with Planck spectral weights, HG/Rayleigh phase sampling
  (Rayleigh via exact Cardano inversion), Beer–Lambert fast path, and
  a deterministic per-pixel-stream orthographic transmittance
  renderer.

## Invariants

- FURNACE: `furnace_radiance` returns exactly `albedo·incident` (energy
  conservation; cosine importance sampling gives zero variance).
- MIS WEIGHT-SUM: the two balance weights at a sample sum to `1` (no energy lost
  or gained at strategy boundaries).
- MIS integration is unbiased (converges to `∫f`).
- Hero-wavelength integration is exact on a constant spectrum and accurate on a
  ramp; `cosine_sample_hemisphere` returns unit vectors in the upper hemisphere.
- Everything is deterministic (low-discrepancy sequences, no RNG here).

- Volumes (vol-001..006): homogeneous slabs match exp(−σL) within
  3σ_stat; heterogeneous means are invariant under a 3× LOOSE
  majorant (48.8k vs 229.3k null collisions ledgered — looseness
  costs work, never bias) and match a deterministic fine-quadrature
  reference; HG E[cosθ] = g (a sign error in the inversion was CAUGHT
  by this gate: −0.5995 measured before the fix) and Rayleigh
  E[cos²θ] = 2/5; spectral emission matches B_λ(T)(1 − e^(−σL)) to
  0.5% at three hero wavelengths; the live LBM dam-break binding
  renders bitwise-replayably through a borrowed buffer with the free
  surface visible (0.917 vs 0.167 transmittance); per-pixel streams
  make any pixel recomputable standalone to bitwise equality.

## Error model

`TraceTermination` reports invalid input/sample, cancellation, iteration-limit,
residual-only stop, miss, or geometrically authorized hit without conflation.
Differentiable rendering returns
`RenderError` for cancellation, invalid parameters/configuration/targets,
backend refusal, uncertified traces, and singular implicit/boundary
derivatives. The tracer returns `TracerError`, preserving cancellation,
invalid dimensions/film buffers/progressive ranges, backend refusal,
uncertified traces, and missing normals. `halton`
panics only on `dim >= 8` (out of the prime table).

## Determinism class

Fully deterministic: the sampling is low-discrepancy, keyed by sample index.

## Cancellation behavior

`sphere_trace` polls its `Cx` before and after each chart evaluation and before
terminal success. Cancellable NURBS seeding/Newton and BVH traversal poll before
and after each bounded seed/iteration/node/triangle. Differentiable scanline rendering
polls at entry, per search iteration, row, pixel, and loss-reduction element and
propagates `RenderError::Cancelled`. The spectral tracer polls per row, sample,
bounce, and primitive, and copies progressive staging buffers in checked
chunks; it propagates `TracerError::Cancelled`. A failed or reversed range
leaves both film sums and `spp_done` unchanged so retry cannot double-count.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

`chart-backends` is a DEFAULT feature. Bead 8ll9 requires its thin-feature
falsifier, deterministic-BVH, workspace/default-matrix, nested-Wasm, and
four-quadrant tracer-golden gates before closeout. The wider SIMD BVH and
ray-rate claims remain evidence-gated successors; default-on does not promote
those claims.

`volumes` [F] gates the volumetric media stack (fs-rand dependency).

`differentiable` (bead qfx.5) gates the edge-aware differentiable renderer
(fs-ad + fs-evidence + fs-math dependencies) and explicitly co-enables the
default chart backend surface. Its primal silhouette and hit decisions use the
same `Chart`/certified-sphere-trace backend; dual lanes lift those decisions by
the implicit hit equation.

`tracer` (bead 872c) gates the spectral path tracer v1
(chart-backends + fs-rand + fs-img): hero-wavelength (4-packet)
NEE+MIS path integration, Lambertian + GGX with spectral reflectance
(the `spectral` module's bounded sigmoid lift; round-trip RGB error
pinned under 1e-3), one rect area light per scene, CIE-XYZ film →
Bradford-adapted linear sRGB → byte-exact EXR. Streams are
counter-based and keyed (pixel, sample, dimension) — Philox for path
decisions, optional Owen-scrambled Sobol' for pixel dimensions
(measured at 64 spp on the Cornell fixture: variance ratio 0.676 vs
iid, ledgered on bead 872c) — so images are bitwise invariant to any
pixel/tile scheduling and progressive checkpoints resume bitwise.
Radiance-path transcendentals go through `fs_math::det`; the Cornell
golden (`fs-render:cornell` in golden-couplings.json) reproduced
identically in all four ISA/profile quadrants at freeze. v1 no-claims:
single NEE light, no volumetric coupling, no environment light, no
Russian roulette, GGX samples the NDF (VNDF is a recorded follow-up),
emitters do not reflect.

## Conformance tests

`tests/render.rs` (7 cases): radical inverse known values; cosine samples are
unit vectors with the right pdf; the furnace test conserves energy exactly; MIS
weights sum to one (+ heuristic ordering); MIS integration is unbiased;
hero-wavelength integration exact on a constant / accurate on a ramp;
determinism.

`tests/diff_battery.rs` (bead qfx.5, feature `differentiable`): edge-aware
gradient vs central FD, the frozen-crossing negative control, shrinking
quadrature bias, inverse rendering, a combined appearance/physics objective,
bitwise primal/gradient replay through the shared backend, a smooth-min seam
derivative regression, and fail-closed cancellation. Numerical receipts are
emitted by the current-tree run; this contract does not carry stale measurements
across backend-semantic changes.

`tests/charts.rs` (beads qfx.2 + 8ll9, default feature): four distinct
thin-shell/scaling falsifiers that all defeat the naive unit-bound marcher while
the certified `d/L` path closes a short first-root bracket; pending-overlap
regressions at both a far shell boundary and `t_max`; witness cancellation,
caller-limit clipping, non-adoption, loose-L, indeterminate-enclosure,
same-sign ultra-thin-shell, and tangent fail-closed controls; 120 additional
grazing-biased rays against a micro-step oracle; explicit no-certificate
behavior when a bound is withheld; analytic NURBS hits; one BVH fingerprint and
bit-identical hit receipt across 1/2/4/8 concurrent builders; mixed-backend and
translated-scene consistency; and honestly labeled throughput telemetry. The
tracer's Cornell EXR golden composes both F-rep sphere tracing and the mesh BVH;
its prior 872c freeze was four-quadrant, and 8ll9 requires current-tree replay.

## No-claim boundaries

- v0 includes the scalar-BVH spectral path tracer. Wide-BVH SIMD traversal,
  watertight ray-triangle tests, a LIGHT-BVH, media coupling, ray-stream
  sorting, and progressive tile streaming to HELM remain staged.
- The spectral pipeline here integrates a spectrum; the radiometrically correct
  spectra→XYZ→display transforms and layered measured-spectrum materials are
  staged.
- `mis_integrate_unit` is a 1-D demonstrator of the balance heuristic; the
  production MIS lives in the path integrator across BSDF/light strategies.

## No-claim boundaries (differentiable)

- Smoke tier is DETERMINISTIC QUADRATURE on SDF scenes (scanline with
  analytic horizontal antialiasing; primal crossings and hits through the
  default chart backend, interior derivatives lifted from the certified hit by
  the implicit equation, and boundary terms through explicit crossing
  velocities with Danskin's envelope at the z-argmin). The
  Monte-Carlo/reparameterized estimators for path-traced integration,
  FrankenTorch-bridged learned BSDFs, heterogeneous differentiable
  `charts::Backend` scenes, fs-xform θ→Region chart perturbations, and fs-opt-ir
  term registration are the RECORDED SUCCESSORS (the loss term's (value,
  gradient) shape is already compatible).
- Vertical antialiasing is sub-row averaging (piecewise constant in
  y): FD steps that push a silhouette tangency across a sub-row line
  see an O(subrow²) kink — fixtures sit away from tangency rows; the
  bias battery measures the induced error honestly.
- The smoke fixture uses deterministic ternary closest-approach search. A
  general proof for arbitrary separated, multi-modal parameter sets is not
  claimed; a certified global 1-D minimum/uniqueness diagnostic is required
  before extending the exact-gradient claim beyond the conformance domain.
- `render_grad(…, edge_terms = false)` exists ONLY as the battery's
  negative control; it is documented WRONG for real gradients.

## No-claim boundaries (charts)

- The tunneling guarantee holds only when a chart opts into
  `TraceStepClaim::{ExactDistance,LipschitzImplicit}` and every sample supplies
  a positive finite Lipschitz bound plus a compatible `Exact` or `Enclosure`
  trace-value certificate containing the reported field value. An `Estimate`
  in `ChartSample.error` remains valid abstract-distance honesty, but it cannot
  substitute for that trace certificate. Charts using the default `NoClaim` may retain an `L = 1` preview
  fallback, but `TraceAudit::certified` is false and production render paths
  reject every terminal result from that trace, including a miss. Malformed
  claims stop as `InvalidSample`.
- A nonzero `LipschitzImplicit` normalized residual is not a Euclidean
  distance-to-boundary certificate. Without a rigorous singleton zero or the
  short opposite-sign witness above, it returns `ResidualLimit` and no `Hit`.
  Same-sign and indeterminate witnesses — including generic tangencies and
  even-contact intervals — remain explicit no-claim outcomes until a chart
  supplies a proximity or first-ray-root certificate. Exact-distance charts
  retain world-space tolerance semantics.
- The mesh BVH is the interim scalar backend; the 8-wide SIMD BVH and
  ray streams are qfx.1's ledgered follow-up scope.
- Ray-rate NUMBERS are measured and ledgered per build/machine; the
  Mray/s TARGETS (80/120) are release-build perf-CI gates (fz2.4), not
  claims this module makes.
- Trimmed-NURBS awareness rides fs-rep-nurbs trim classification; the
  intersection here treats the full patch (no-claim on trimmed holes).

## No-claim boundaries (volumes)

- FrankenVDB tile-maxima majorants: no fvdb crate exists in-workspace;
  [`MajorantGrid`] builds per-block maxima from dense grids, and the
  per-tile-rate DDA traversal (rather than lookup thinning under a
  global bound) is the recorded successor alongside the FVDB wiring.
- Progressive live tiles with ledger artifact pinning (frame-consistent
  snapshots of evolving fields) — staged with the vessel flagship's
  render lane; the smoke tier renders a paused simulation's buffer.
- Refractive free-surface rendering (fill-fraction interface
  reconstruction) and MIS integration of phase functions into the full
  tracer — successors; the phase samplers and their moment gates ship
  now.
- The zero-copy claim at smoke tier is BORROW SEMANTICS (the API takes
  `&[f64]`; the battery binds a live `FreeSurface` mass buffer); the
  FrankenNumpy membrane view protocol is the fuller deliverable.
