# CONTRACT: fs-render

Unbiased spectral path-tracing core: the verifiable Monte-Carlo foundations.

## Purpose and layer

Layer L5 (LUMEN). No dependencies — pure Rust.

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

- `charts` module (plan §10.2, bead qfx.2; [F], behind
  `chart-backends`): render whatever chart exists, WITHOUT conversion.
  `sphere_trace` steps `|f(p)|/L` with the chart's CERTIFIED Lipschitz
  bound — the sign cannot flip within that radius, so the marcher
  provably never tunnels (audited: `TraceAudit.worst_step_ratio`);
  over-relaxation uses the standard certified fallback (retreat when
  spheres fail to overlap). `ray_intersect_nurbs` is grid-seeded 3×3
  Newton on `S(u,v) − o − t·d` with the `[S_u, S_v, −d]` Jacobian.
  `TriMesh` is Möller–Trumbore over a deterministic median-split BVH.
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

Total functions; `halton` panics only on `dim >= 8` (out of the prime table).

## Determinism class

Fully deterministic: the sampling is low-discrepancy, keyed by sample index.

## Cancellation behavior

None here; the production tracer polls `Cx` at tile boundaries (a render is a
budgeted, cancellable study).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None (the `frontier-polarization` Mueller-calculus path is staged).

`volumes` gates the volumetric media stack (fs-rand dependency) per
the same Ambition-Tag rule as `chart-backends`.

`differentiable` (bead qfx.5) gates the edge-aware differentiable
renderer (fs-ad + fs-math dependencies), same Ambition-Tag rule.

## Conformance tests

`tests/render.rs` (7 cases): radical inverse known values; cosine samples are
unit vectors with the right pdf; the furnace test conserves energy exactly; MIS
weights sum to one (+ heuristic ordering); MIS integration is unbiased;
hero-wavelength integration exact on a constant / accurate on a ramp;
determinism.

`tests/diff_battery.rs` (bead qfx.5, feature `differentiable`, 6
cases): dr-001 edge-aware gradient of the full-image L2 loss vs
central FD of the RENDER — worst rel 5.9e-6 over all 9 parameters
(centers, radii, blend); dr-002 the NEGATIVE CONTROL — freezing the
silhouette crossings (what naive pointwise autodiff computes) is off
by 1.4e-1 relative where the edge-aware gradient sits at 1.6e-8 (the
boundary term is ~9e6× the honest error budget: measured, not
asserted); dr-003 deterministic-quadrature bias shrinks monotonically
with supersampling (2.8e-1 → 1.6e-2); dr-004 inverse rendering
recovers sphere position + radius to 7e-7 from a target image;
dr-005 combined appearance+physics objective (image L2 + volume
budget) optimized end-to-end through the shared gradient path
(1.1e-2 → 8.9e-12, budget met); dr-006 bitwise replay of render AND
gradient.

## No-claim boundaries

- v0 is the verifiable Monte-Carlo core (sampling, furnace, MIS, spectral
  integration). The full unidirectional PATH TRACER — wide-BVH SIMD traversal,
  watertight ray-triangle tests, next-event estimation with a LIGHT-BVH,
  Beer–Lambert media, ray-stream sorting, progressive tile streaming to HELM,
  per-tile Philox keyed by (seed, frame, tile), and `Cx` cancellation — is the
  fuller deliverable, staged.
- The spectral pipeline here integrates a spectrum; the radiometrically correct
  spectra→XYZ→display transforms and layered measured-spectrum materials are
  staged.
- `mis_integrate_unit` is a 1-D demonstrator of the balance heuristic; the
  production MIS lives in the path integrator across BSDF/light strategies.

## No-claim boundaries (differentiable)

- Smoke tier is DETERMINISTIC QUADRATURE on SDF scenes (scanline with
  analytic horizontal antialiasing; interior terms through converged
  dual sphere tracing, boundary terms through explicit crossing
  velocities with Danskin's envelope at the z-argmin). The
  Monte-Carlo/reparameterized estimators for path-traced integration,
  FrankenTorch-bridged learned BSDFs, unification with
  `charts::Backend`, fs-xform θ→Region chart perturbations, and
  fs-opt-ir term registration are the RECORDED SUCCESSORS (the loss
  term's (value, gradient) shape is already compatible).
- Vertical antialiasing is sub-row averaging (piecewise constant in
  y): FD steps that push a silhouette tangency across a sub-row line
  see an O(subrow²) kink — fixtures sit away from tangency rows; the
  bias battery measures the induced error honestly.
- `render_grad(…, edge_terms = false)` exists ONLY as the battery's
  negative control; it is documented WRONG for real gradients.

## No-claim boundaries (charts)

- The tunneling guarantee holds for charts whose `lipschitz` claim is
  certified (Frep/exact SDFs); charts reporting no bound default to
  L = 1, which is only safe for true distance fields.
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
