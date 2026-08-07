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
- `dielectric` module (feature `tracer`): validated Cauchy phase-index laws in
  vacuum nanometres, homogeneous Beer-Lambert absorption in inverse metres,
  provenance-bearing representative glass presets, exact unpolarized Fresnel,
  Snell refraction and total internal reflection, plus smooth-delta and
  single-scattering isotropic-GGX reflection/transmission. Radiance transport
  uses the documented `(eta_i/eta_t)^2` transmission convention.
- `conductor` module (feature `tracer`): fixed-capacity visible-band tables of
  absolute complex refractive index `eta(lambda) + i k(lambda)` in vacuum
  nanometres, strict inclusive 380--780 nm support, deterministic linear
  interpolation, exact unpolarized complex-IOR Fresnel against the path's active
  incident medium, and validated isotropic-GGX roughness. Built-in tungsten and
  stainless-steel tables are explicitly representative look-development data;
  table/source/material identities prove replay binding, not specimen fidelity.

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
  Newton on `S(u,v) − o − t·d` with the `[S_u, S_v, −d]` Jacobian. Its
  allocation-bearing legacy seed grid is checked, capped, and fallibly
  reserved before evaluation; defensive sealed-source invariant failure,
  invalid settings, resource refusal, cancellation, and bounded-search
  nonconvergence remain distinct typed outcomes. The current Newton path has no exclusion theorem,
  so exhausted starts fail closed as `IterationLimit`; `Ok(None)` is reserved
  for a future certified miss. Jacobian columns and residuals are normalized
  before the dimensionless angular-volume condition test, and residual/normal
  norms use scale-safe finite arithmetic, so knot-domain or ray scaling cannot
  manufacture convergence, a miss, or a `Some(NaN)` normal.
  Certification requires the chart-level typed `TraceStepClaim`; a sample
  carrying `Some(L)` cannot upgrade the default `NoClaim`. Exact-distance hits
  use world-space distance tolerance. A Lipschitz-implicit chart authorizes a
  geometric `Hit` with either a rigorous singleton-zero field enclosure or a
  short sign bracket. At each strict-sign residual sample, the marcher may
  inspect one cancellation-aware witness no more than `2*eps` ahead in actual
  outward-rounded Euclidean distance. Rigorously opposite endpoint signs plus
  the claim's finite-segment continuity prove a boundary in that segment. When
  that bounded witness is same-sign or otherwise supplies no root evidence, it
  remains non-adopted: the marcher may advance only to the independently
  certified `|f|/L` safe endpoint and retry under the caller-admitted step
  budget. The public default is 4096 steps; explicit requests are hard-capped
  at 16384. Residual refinement never launches an over-relaxed endpoint.
  The marcher uses bounded evaluated refinement when binary64 `Ray::at` rounding
  puts the first midpoint microscopically outside an endpoint-distance guard.
  It preserves the original root-free prefix endpoint and may pull back only a
  rigorously opposite-sign or singleton-zero far endpoint; same-sign,
  indeterminate, and no-progress refinement fails closed rather than discarding
  possible earlier even crossings. Only an evaluated representative within
  `eps` of both retained endpoints becomes the hit. The witness is evidence
  only and is never adopted as march state. A nonzero normalized residual
  `|f|/L` by itself certifies only the next root-free safe step, not proximity.
  Reaching the caller boundary while it remains residual, exhausting bounded
  residual refinement, or losing finite forward progress without a sign
  bracket stops as `ResidualLimit` with no `Hit`. Pending over-relaxed
  endpoints are validated before either hit or miss acceptance. An
  outward-rounded working limit is classified at the caller's actual endpoint:
  hits require a validated endpoint sample, and misses require an
  epsilon-clear safe-ball bridge across any coordinate-rounding gap. Normalized
  working parameters size steps, but every chart and overlap evaluation uses
  the caller ray's own `Ray::at` arithmetic so a certificate is never returned
  for a numerically different point. Mesh hits additionally retain original
  triangle indices and barycentric coordinates; exact-distance feature ties
  select the lowest original triangle index. This is chart-backend
  bit-semantics v10.
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

- `instances` module (bead `frankensim-h7xu5.3.1`, default-on through
  `chart-backends`): `RigidTransform` admits only finite proper-rigid
  body-to-world placements (a near-unit quaternion and translation in metres),
  canonicalizes the quaternion double cover, and cannot represent scale,
  shear, or reflection. `SharedGeometry` keeps one immutable chart or mesh
  allocation behind `Arc`; `GeometryInstance` binds it to a nonzero stable
  object ID, a caller-supplied geometry content identity, and a placement.
  Intersection transforms the world ray into local coordinates, delegates to
  the existing certified chart or deterministic mesh backend, preserves the ray
  parameter and backend audit, and rotates geometric/shading normals, tangents,
  and surface derivatives back to world space. `InstanceScene` rejects object-ID
  collisions, orders instances by object ID, and content-identifies the ordered
  placements with a length-framed streaming hash. Material and emission remain
  properties of the tracer `Primitive`, separate from shared geometry and pose.

- `tracer` tile execution (bead `frankensim-h7xu5.5.1`):
  `RenderExecutionConfig` makes logical tile shape, worker ceiling, operation
  memory, scheduling weights, and `RunId` explicit. The
  `render_*_with_execution` convenience functions preserve the serial APIs as
  bitwise oracles while executing fixed row-major tiles through `fs-exec`.
  `RenderWorkerPool::with_parked_crew_local` is the animation/batch surface:
  one scoped worker crew serves every `ParkedRenderScope` job and joins
  structurally when the callback exits. Each success reports requested versus
  admitted workers, layout, setup/traversal/compute/merge/publication timings,
  executor drain diagnostics, and the operation-memory receipt.

- `tracer` deterministic adaptive sampling (bead
  `frankensim-h7xu5.5.2`): `AdaptiveSamplingConfig` declares a minimum,
  maximum through `Settings::spp`, fixed decision-batch spacing, per-channel
  absolute and relative dispersion thresholds, and a dark-channel scale.
  `AdaptiveFilm` privately owns the raw sequential XYZ sum, Welford mean and
  second central moment, exact sample count, and terminal decision for every
  pixel. It binds sampler, stream seed, policy, maximum SPP, camera/time mode,
  and `ADAPTIVE_SAMPLING_SEMANTICS_VERSION`. This is diagnostic estimator
  provenance, not a complete replay identity: the current film does not bind
  full render settings or scene content, which remain obligations of the
  checkpoint/artifact layer. `beauty_mean_xyz` always divides the raw sum,
  preserving the uniform renderer as the oracle;
  `estimator_mean_xyz` exposes the separately retained Welford mean used by
  stopping. Decisions occur only at the declared fixed checkpoints, require
  all three XYZ channels, and record threshold success in preference to the
  hard ceiling when both occur at the final checkpoint. Static, legacy-motion,
  cinematic-camera, scoped-worker, parked-crew, and opaque resumable entry
  points share the same estimator.

- `aov` module (bead `frankensim-h7xu5.6.1`): opt-in cinematic diagnostic and
  denoising accumulation is deliberately parallel to legacy `Film`; ordinary
  RGB rendering, RGB EXR output, checkpoint bytes, and shard bytes neither
  allocate nor carry these planes. `CinematicAovConfig` freezes one profile,
  L6-supplied finite ordered current/previous/next frame times and nonzero
  trajectory, scene, and composition assertions, plus explicit pixel, retained
  payload, EXR-plane, EXR-metadata, EXR-encoder-scratch, encoded-EXR-output,
  and palette limits. Every limit is part of configuration identity and uniform
  checkpoint state. Those supplied identities are retained provenance
  assertions, not scene or trajectory authority minted from borrowed renderer
  data. `BeautyOnly` exports the three linear-sRGB beauty channels;
  `DailyCore` additionally exports linear-sRGB reflectance albedo, the
  normalized face-forwarded world-space surface frame actually consumed by
  beauty, axial-metre depth, `primary.coverage`, unbiased raw CIE-Y
  sample variance, and previous-minus-current raster-pixel motion (14 `FLOAT`
  channels total); `FinalDiagnostic` additionally exports world-unit geometric
  normal, direct/indirect/emission linear-sRGB beauty splits, object/material
  palette indices, per-pixel sample count, and `diagnostic.validity` (30
  `FLOAT` channels total). Surface-filtered values average only qualifying
  primary hits; no primary hit produces zero surface values and zero
  `primary.coverage`. `primary.coverage` is the primary-hit fraction, so it is
  the DailyCore validity indicator; FinalDiagnostic additionally names valid
  primary, albedo, previous-motion, object-ID, material-ID, and
  contribution-split observations in its bit mask. The authored-shading-normal
  validity bit is reserved and remains clear in AOV semantics v2 because the
  beauty transport does not yet consume backend-authored shading normals.
  Palette
  index zero always means background or unavailable; nonzero object IDs and
  material content hashes are stored losslessly in sorted EXR metadata palettes.
  Final categorical IDs select the accepted primary closest to the pixel centre,
  with deterministic ties by absolute sample, primitive index, object palette,
  then material palette. A direct mesh/chart hit can therefore retain material
  identity while correctly leaving instance object identity and motion
  unavailable; it is never assigned a fabricated instance ID.

  `render_cinematic_with_aovs` and transactional
  `render_cinematic_range_with_aovs` obtain beauty, primary, normals, material,
  contribution split, sampled shutter time, depth, and motion from the same
  accepted path traversal; they do not retrace to populate diagnostics. A range
  stages the complete beauty/AOV state and swaps it only after success. Its
  process-local continuity guard rejects a progressive append when borrowed
  scene, lighting, geometry, or camera inputs changed; chart allocation
  addresses participate only in that in-process guard, which is not exported
  and is not an authority-bearing content identity. `render_cinematic_adaptive_with_aovs`
  retains the adaptive film's exact per-pixel terminal prefix and matching AOV
  observations. AOV samples never participate in the adaptive stopping decision
  and cannot change its terminal count. `FinalDiagnostic` exports the actual
  per-pixel count in `samples`; narrower frozen profiles retain counts
  internally and set `frankensim.render.spp=per-pixel-unexported-by-profile`
  instead of naming a nonexistent channel.

  `to_exr` emits deterministic `FLOAT` planes and raw-estimate metadata naming
  the AOV schema/profile/configuration hash, channel and invalid-value
  semantics, material domain, frame/provenance assertions, render sampler and
  strategy, uniform or adaptive sample mode and policy, exact shutter
  convention/distribution/strata, renderer versions, and exact palettes. The
  `frankensim.aov.authority=raw-estimate` marker is
  intentional: these images do not certify material truth, target-frame
  visibility, physical trajectory truth, mechanical state, or an image-error
  bound. Retained-payload admission covers owned beauty/AOV accumulator payloads
  (and the uniform progressive full-frame staging copy); export admission covers
  four independently refused logical categories: uncompressed `f32` plane
  storage, channel-descriptor plus metadata payloads, EXR canonical-order
  reference scratch, and the exact encoded EXR artifact length. Scene-sized
  palette metadata is bounded before construction; the encoder sizes a 4K
  layout without allocating its planes, reserves its output once, and returns
  no partial artifact on refusal. The public defaults admit the 30-channel 4K
  profile (995,328,000 plane bytes and an encoded artifact below 1 GiB), while
  preserving distinct ceilings so callers can choose a smaller profile or
  envelope. These logical payload limits do not claim allocator bookkeeping,
  allocator usable-size rounding, process RSS, or throughput; the current
  returned-`Vec` path necessarily holds staged planes and encoded output
  concurrently during finalization.

  `CinematicAovFilm::denoise_guides` exposes the retained `DailyCore` guide
  subset as owned row-major `f32` planes: previous-frame motion, axial depth,
  shading normal, primary coverage, and raw-luminance variance. It uses the
  exact same averaging, normal normalization, finite narrowing, invalid
  background-zero convention, and signed-zero canonicalization as EXR export;
  the focused battery compares every returned plane bit-for-bit with its
  decoded EXR counterpart. `BeautyOnly` refuses because it retained no primary
  observations. Stable IDs are intentionally omitted rather than converting
  per-frame palette indices into a stronger identity claim.

  Progressive uniform `CinematicAovFilm` state has a distinct `FSRAOVC1`
  checkpoint codec; legacy `FSRCP001` bytes remain unchanged. The AOV codec
  streams at most 64 KiB payload chunks plus one 32-byte seal under an explicit
  encoded-byte ceiling. It stores every raw beauty/AOV accumulator field,
  committed uniform sample prefix, exact profile/configuration/provenance,
  all four export ceilings, settings, shutter/shot/cut-side binding, sorted
  palettes, and the process-local
  continuity fingerprint field-by-field in little-endian form. Restore checks
  the domain-separated seal, exact profile-dependent length, semantic versions,
  canonical finite values, count relationships, palette ordering/ranges, and
  render binding before publishing a film. The caller must also supply the
  independently expected complete AOV configuration and exact uniform render
  binding, including its committed sample prefix; an internally coherent
  checkpoint for another job or an older prefix is refused.
  Restoring a committed prefix and
  continuing through `render_cinematic_range_with_aovs` is bitwise equivalent
  to straight-through accumulation. The seal detects accidental corruption or
  unresealed substitution; it provides neither authentication nor adversarial
  integrity, and the supplied L6 identities
  retain their assertion-only authority. Because a direct chart currently has
  no content identity, its continuity guard contains a process-local allocation
  address and a cross-process chart resume correctly refuses rather than
  pretending identity.

  Fresh uniform AOV films also have an explicit deterministic tile path:
  `render_cinematic_with_aovs_execution` and
  `ParkedRenderScope::render_cinematic_with_aovs`. Each logical tile owns
  disjoint pixels and accumulates every pixel in ascending absolute-sample
  order before copying the complete tile into one private full-frame staging
  film. The path reuses `RenderExecutionConfig`, the ambient `Cx`, the existing
  operation lease, and the ordinary/parked `RenderPoolRunner`; worker count,
  stealing, and tile completion order cannot change AOV arithmetic. The fresh
  path allocates one staged AOV film rather than constructing an empty film and
  cloning it as the progressive compatibility API must. Cancellation, AOV
  refusal, memory refusal, or a contained worker panic drains the crew and
  publishes no film.

  There is currently no AOV-specific tile-pending state, adaptive checkpoint
  codec/resume state, shard specification, shard result, or shard merge API.
  The existing tracer checkpoint and uniform-shard surfaces still apply only to
  their legacy RGB films and must not be represented as adaptive persistence or
  distributed merge support for cinematic AOV films.

- `tracer` crash-recovery codec (bead `frankensim-h7xu5.5.3`):
  `PendingRender` and `PendingAdaptiveRender` stream a schema-versioned,
  domain-hashed checkpoint in at most 64 KiB body chunks and one fixed seal.
  The canonical body binds uniform/adaptive kind, every bit-affecting tracer
  semantics version, complete settings/time/execution/adaptive policy, the
  caller's exact execution `Budget`, the runtime ISA and sorted detected
  feature set, exact tile-row prefixes and attempt count, raw binary64
  accumulators/AOVs, and an L6-supplied
  source/configuration/scene/frame/job/build/claim/generation chain. The
  renderer derives and validates the render-job identity from those owned
  inputs; an L6 binding cannot substitute a caller-invented job digest.
  Generation zero has no predecessor; every later generation names a nonzero
  prior renderer-content digest. Emission accepts a caller byte ceiling and
  fallible sink. It snapshots at most one committed tile row while holding the
  pending-state lock and never calls the external sink under that lock. Restore
  consumes a freshly admitted opaque job, verifies the complete seal and
  binding before decoding, refuses malformed numeric or row-prefix state, and
  additionally verifies adaptive raw-sum/mean/count consistency plus canonical
  nonnegative second moments. Successor restore consumes an already restored
  predecessor job and refuses unless tile prefixes and attempts are
  nondecreasing and every predecessor-committed uniform accumulator or adaptive
  sum/moment/count/decision is retained bit-for-bit. It returns no partially
  restored job on error.
  The codec never publishes a film or durable artifact by itself; only a
  completed film may enter final image/manifest publication.

- `tracer` uniform render sharding (bead `frankensim-h7xu5.5.4`):
  `UniformRenderShardSpec` admits one nonempty rectangle of the canonical
  `(row-major tile, absolute sample)` space under explicit path and encoded-byte
  caps. Its identity binds the external plan and frame identities, stable frame
  ordinal, complete fixed-SPP settings and time mode, tile layout and ranges,
  renderer semantic versions, and runtime ISA/features. Workers trace only that
  rectangle and return finite immutable XYZ partial sums; random streams remain
  keyed by absolute pixel/sample identity. The strict canonical result codec
  requires external plan and shard pins and rejects truncation, trailing bytes,
  corruption, non-finite payloads, or a foreign header before returning a
  result. `merge_uniform_shards` derives the exact complete-set input envelope
  from expected specs before allocating merge indexes, validates exact
  nonoverlapping full coverage, and ignores submitted exact duplicates without
  changing that envelope. Expected and submitted order are nonsemantic:
  reference semantics, diagnostics, conflicts, and accumulation use canonical
  ordering or fixed error precedence. Missing, foreign, corrupt, or conflicting
  work and aggregate-input/output-film cap violations publish no partial film.
  Full-SPP tile-only plans retain the legacy serial accumulation order and bits.
  A sample-split result is bit-stable across worker counts and arrival orders
  for the same frozen plan, but is not claimed bit-identical to the monolithic
  accumulator because binary64 addition is non-associative. This surface does
  not implement adaptive sample repartition, a process launcher, remote
  transport, leases, or cluster fault tolerance.

- `volumes` module (bead qfx.3, feature `volumes`): [`VolumeGrid`]
  BORROWS its density buffer (zero-copy: live simulation fields render
  in place), [`MajorantGrid`] per-block maxima, Woodcock delta
  tracking (`woodcock_transmittance`, unbiased for ANY bound ≥ max σ;
  the tile stage thins field lookups), the collision emission
  estimator with Planck spectral weights, HG/Rayleigh phase sampling
  (Rayleigh via exact Cardano inversion), Beer–Lambert fast path, and
  deterministic per-pixel-stream orthographic renderers. The
  transfer-function path uses a validated piecewise-linear
  `TransferFunction` from scalar value to nonnegative extinction and
  linear-RGB source radiance. `render_transfer_emission` derives the
  mapped global majorant from the borrowed field and admits explicit
  field-scan, pixel, primary-sample, and per-sample null-collision
  budgets before rendering. Its bit contract is
  `TRANSFER_RENDER_SEMANTICS_VERSION = 1`.

## Invariants

- FURNACE: `furnace_radiance` returns exactly `albedo·incident` (energy
  conservation; cosine importance sampling gives zero variance).
- G3 RADIANCE RESCALING: coherently rescaling incident radiance rescales the
  furnace result by the same factor under an explicit numeric tolerance.
- MIS WEIGHT-SUM: the two balance weights at a sample sum to `1` (no energy lost
  or gained at strategy boundaries).
- MIS integration is unbiased (converges to `∫f`).
- Hero-wavelength integration is exact on a constant spectrum and accurate on a
  ramp; `cosine_sample_hemisphere` returns unit vectors in the upper hemisphere.
- Smooth dielectric Fresnel branch probabilities partition unity; a lossless
  entry/exit slab cancels its reciprocal radiance eta factors. Rough dielectric
  evaluation and sampling share the Walter solid-angle Jacobian and remain
  finite and nonnegative over the admitted grazing/roughness grid. Homogeneous
  attenuation composes as `T(L1 + L2) = T(L1) T(L2)` to floating-point
  tolerance.
- Sampling is deterministic: analytic and low-discrepancy helpers are
  stateless, while stochastic render paths use explicitly keyed counter-based
  streams rather than ambient mutable RNG state.
- Proper-rigid placement preserves lengths and therefore preserves the local
  backend's ray parameter. Exact-distance instance ties select the lowest
  object ID, independent of caller insertion order. A geometry identity stays
  unchanged when only its placement changes; the frame identity binds object,
  geometry, and the canonical placement.

- Volumes (vol-001..009): homogeneous slabs match exp(−σL) within
  3σ_stat; heterogeneous means are invariant under a 3× LOOSE
  majorant (48.8k vs 229.3k null collisions ledgered — looseness
  costs work, never bias) and match a deterministic fine-quadrature
  reference; HG E[cosθ] = g (a sign error in the inversion was CAUGHT
  by this gate: −0.5995 measured before the fix) and Rayleigh
  E[cos²θ] = 2/5; spectral emission matches B_λ(T)(1 − e^(−σL)) to
  0.5% at three hero wavelengths; the live LBM dam-break binding
  renders bitwise-replayably through a borrowed buffer with the free
  surface visible (0.917 vs 0.167 transmittance); per-pixel streams
  make any pixel recomputable standalone to bitwise equality. Transfer
  construction rejects unordered, non-finite, or negative optical
  knots; interpolation and endpoint clamping are deterministic; a
  mapped homogeneous slab matches
  `source_rgb * (1 - exp(-extinction * length))` channel by channel;
  replay is bitwise; insufficient work budgets, non-finite fields,
  insufficient majorants, and exhausted tracking limits return typed
  refusal instead of partial or biased images.

## Error model

`TraceTermination` reports invalid input/sample, cancellation, iteration-limit,
residual-only stop, miss, or geometrically authorized hit without conflation.
Differentiable rendering returns
`RenderError` for cancellation, invalid parameters/configuration/targets,
backend refusal, uncertified traces, and singular implicit/boundary
derivatives. `RenderCfg.max_trace_steps` makes its per-ray work envelope
explicit; zero or values above the hard 16384-step ceiling are invalid. The
tracer returns `TracerError`, preserving cancellation,
invalid dimensions/film buffers/progressive ranges, backend refusal,
uncertified traces, missing normals, invalid/colliding rigid instances,
dielectric evaluation refusal, LIFO medium-stack mismatch/overflow, and a miss
while still inside a declared closed medium.
Instance construction rejects zero identities, invalid transforms, and any
scene count that cannot fit the canonical identity encoding;
intersection rejects malformed rays, non-positive/non-finite limits,
non-finite hit data, and missing geometric normals rather than manufacturing
shading data. `halton`
panics only on `dim >= 8` (out of the prime table).
Transfer construction and direct-volume admission return `DvrError`.
The high-level renderer reserves its complete private image buffer before
sampling and drops it on any tracking refusal; no error path returns a partial
image.
The tile renderer additionally reserves retained/staging film payloads, one
shared three-dimensional Sobol direction table for Owen-Sobol jobs, its worst-case concurrent
tile-pixel scratch envelope, and `fs-exec` root metadata before dispatch. Tile
scratch uses fallible raw buffers only inside that already-held aggregate
charge, so scheduling overlap cannot change lease admission. The dielectric
medium stack is fixed-capacity inline state. Thread stacks, allocator
usable-size overhead, and heap owned by lower chart implementations remain
outside the operation-memory claim.
Fresh cinematic AOV tile jobs charge their complete private beauty/AOV film,
the exact profile-dependent maximum concurrent tile accumulators, sampler
state, and executor metadata before dispatch; successful publication releases
every transient charge while the returned film retains its owned storage.
Adaptive jobs reserve three XYZ binary64 planes (raw sum, Welford mean, and
second moment), one `u32` count plane, one terminal-decision plane, sampler
state when applicable, and tile or row scratch before dispatch. Published AOV
allocations leave the operation lease while remaining owned by the returned
film. This is allocation accounting, not a 4K resident-set or throughput claim.

## Determinism class

SAME-ISA bit-deterministic: the sampling is low-discrepancy, keyed by sample
index. The complete `RenderCfg`, including `max_trace_steps`, is replay-critical
input because a different bounded work envelope may change a typed refusal into
a certified hit and therefore change image bytes. The transfer renderer uses a
distinct per-pixel Philox domain, so
pixel/tile execution order does not alter its samples. The claim is scoped
to one ISA/libm build (determinism-tier policy, bead frankensim-lyms):
volume tracking, phase functions, and disk sampling call platform libm
(`exp`/`ln`/`cos`/`sin`/`powf`/`cbrt`), which is not correctly rounded and
may differ by last-ULP across ISAs or libm versions. Cross-ISA bitwise
replay of renders is deliberately NOT claimed; promote by routing through
`fs_math::det` and registering in `check-libm` if a lane ever needs it.
Instance ordering and exact-hit ties are deterministic by stable object ID.
Canonical quaternion sign and signed-zero normalization ensure equivalent
proper-rigid inputs produce the same transform and frame identities.
Tile geometry is independent of worker count. Each tile owns disjoint pixels,
and each pixel evaluates samples in ascending absolute sample order with the
same `(pixel, sample, dimension)` stream and floating-point expression as the
serial oracle. No cross-tile floating-point reduction exists, so worker counts,
steals, scheduling weights, and parked-versus-spawned worker lifetime do not
change film bits on the same ISA.
Adaptive pixels likewise consume the unchanged absolute sample prefix
`0..terminal_count`; active-mask shape, batch checkpoint spacing, tile shape,
worker count, scheduling weights, and resume attempts never reseed or compact a
continuing pixel's stream. The raw sum preserves the uniform tracer's
sequential per-pixel addition order. Welford state is updated transactionally
across all XYZ channels, so a non-finite sample or invalid intermediate moment
publishes none of that sample.

## Cancellation behavior

`sphere_trace` polls its `Cx` before and after each chart evaluation and before
terminal success. Cancellable NURBS seeding/Newton and BVH traversal poll before
and after each bounded seed/iteration/node/triangle. Differentiable scanline rendering
polls at entry, per search iteration, row, pixel, and loss-reduction element and
propagates `RenderError::Cancelled`. Spectral scene preflight polls while
validating each light and primitive and while materializing canonical light
order; both progressive and shard entry points map that refusal to
`TracerError::Cancelled` without publishing a partial result. The spectral
tracer also polls per row, sample, bounce, and primitive, and copies progressive
staging buffers in checked chunks. A failed or reversed range
leaves both film sums and `spp_done` unchanged so retry cannot double-count.
Tile execution derives the exact gate and budget from the caller's `Cx`, polls
inside path traversal as well as at tile boundaries, contains panics, and fully
drains the worker lane before returning. A public film is swapped only after a
successful complete executor report and a final caller-gate checkpoint.
Parked crews use the same run protocol and join on callback exit or unwind.
`PendingAdaptiveRender` commits one complete tile row at a time into opaque
private raw-sum/moment/count/decision buffers. Cancellation or a contained
worker panic discards the uncommitted row, retains completed row prefixes, and
recomputes the row from the same absolute sample IDs under a fresh `Cx`.
Policy, estimator version, sampler/seed, scene/camera borrow, shutter mode,
layout, `RunId`, execution mode, and budget cannot be substituted on resume.
Instance traversal polls between objects and delegates to the backend's existing
bounded cancellation points; cancellation is propagated without a partial hit.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

`chart-backends` is a DEFAULT feature. Bead 8ll9 requires its thin-feature
falsifier, deterministic-BVH, workspace/default-matrix, nested-Wasm, and
four-quadrant tracer-golden gates before closeout. The wider SIMD BVH and
ray-rate claims remain evidence-gated successors; default-on does not promote
those claims.

`volumes` [F] gates the volumetric media stack (fs-rand dependency).
The v1 transfer renderer claims emission/absorption only. It does not claim
scattering, preintegrated transfer functions, adaptive ray integration,
display encoding, tone mapping, or ledger/EXR provenance embedding; those are
separate composition and evidence lanes.

`differentiable` (bead qfx.5) gates the edge-aware differentiable renderer
(fs-ad + fs-evidence + fs-math dependencies) and explicitly co-enables the
default chart backend surface. Its primal silhouette and hit decisions use the
same `Chart`/certified-sphere-trace backend; dual lanes lift those decisions by
the implicit hit equation.

`tracer` (bead 872c) gates the spectral path tracer
(chart-backends + fs-rand + fs-img): hero-wavelength (4-packet)
NEE+MIS path integration, Lambertian + reflective GGX with spectral
reflectance, exact-Fresnel spectral conductors, and provenance-bearing
smooth/rough spectral dielectric glass
(the `spectral` module's bounded sigmoid lift; round-trip RGB error
pinned under 1e-3), deterministic multiple rectangular emitters plus an
optional canonical lat-long environment, CIE-XYZ film →
Bradford-adapted linear sRGB → byte-exact EXR. Streams are
counter-based and keyed (pixel, sample, dimension) — Philox for path
decisions, optional Owen-scrambled Sobol' for pixel dimensions
(measured at 64 spp on the Cornell fixture: variance ratio 0.676 vs
iid, ledgered on bead 872c) — so images are bitwise invariant to any
pixel/tile scheduling and progressive checkpoints resume bitwise.
Radiance-path transcendentals go through `fs_math::det`; the Cornell
golden (`fs-render:cornell` in golden-couplings.json) reproduced
identically in all four ISA/profile quadrants at freeze. The lighting-v1
extension admits each rectangle only when its named primitive resolves to the
same two-triangle world-space quad (directly or through a static mesh instance),
then orders semantic light identities by cancellation-controlled keyed
insertion independently of construction order. It
selects among incident solid-angle/luminance rectangle weights and the
environment's integrated luminance. Environment texels use exact spherical-cell
solid angles and a two-level row/column CDF, so selection is logarithmic in map
dimensions; rotation is around world +Y, with increasing columns running from
+X toward +Z. A semantic pixel hash drives deterministic sample ordering while
the separate source/provenance hashes retain container lineage without changing
the sampling stream. Both NEE and BSDF-hit paths include the same selection
probability in their solid-angle PDFs. At a lighting-v1 path's final permitted
bounce, NEE has weight one because no competing BSDF continuation is evaluated.
The exact legacy one-rectangle/no-environment lighting branch retains its
original light-selection draws and estimator. Tracer bit-semantics v2
deliberately changes reflective GGX sample/PDF bits from NDF sampling to
isotropic visible-normal sampling while preserving the same BSDF integral in
expectation.
Current no-claims: no volumetric coupling, no Russian roulette, reflective GGX
and conductor sampling is isotropic-only, rough dielectric sampling remains an
NDF-plus-Fresnel-branch path, conductor and dielectric GGX are single-scattering
without multiple-scattering compensation, and emitters do not reflect. On the
first sampled transmission through a dispersive dielectric, the shared
four-wavelength geometric packet collapses once to its uniformly selected hero
lane with the matching factor-four estimator weight; companion lanes are zeroed
instead of being biased along the hero wavelength's refracted direction.
Current-vertex NEE and sampled reflection retain all four lanes because their
geometric directions are wavelength-independent; this preserves the same
spectral transport while avoiding a gratuitous fourfold reduction in spectral
samples on camera-visible glass reflections.
Absorption uses unshifted physical segment length. Smooth events have zero
solid-angle query density and receive delta-correct MIS treatment. A strict
path-local stack mutates only after sampled transmission and supports nested,
closed, consistently outward-oriented media.

## Conformance tests

`tests/instances.rs` (bead `frankensim-h7xu5.3.1`, feature
`chart-backends`, with one `tracer`-gated E2E case): transform admission and
quaternion canonicalization; inverse/composition round trips; identity and
placed-mesh equivalence; world-space differential frames; certified chart
authority preservation; shared-allocation reuse; collision rejection and
object-ID tie order; missing-normal, no-claim, and cancellation refusal; stable
geometry/frame identities; and a production spectral render whose visible
result changes when the immutable mesh instance moves out of view. Exact chart
tangency is pinned as a bounded fail-closed backend outcome, not presented as a
certified hit. A second tracer E2E regression proves exact-tie output is
independent of primitive insertion order and that duplicate object IDs refuse
before sampling.

`tests/dielectric_battery.rs` (feature `tracer`): equal-IOR null slabs across all
direct-light strategies; a translation-invariant 0.1 mm slab at `x = 10^9`;
exact thin/thick Beer-Lambert scaling; nested active-medium accounting;
reversed-winding and non-LIFO transactional refusal; finite rough transmitted
NEE with target-medium attenuation; cancellation; and bitwise progressive
replay. Inline analytic tests add independent Fresnel, Snell/critical-angle,
signed-vector refraction, Walter rough-transmission, eta-factor, signed-zero,
grazing, pole-frame, adjacent-IOR, reflected-packet retention, and one-time
post-transmission dispersive packet-collapse fixtures. The Cornell tracer
battery is deliberately re-frozen when tracer bit
semantics change and remains the opaque-path non-regression gate.

`conductor` and `tracer` inline tests plus the Euler scene-bridge E2E fixture
(feature `tracer`, and `cinematic-render` downstream) cover complex-Fresnel
known answers and an independent real-form grid, normal/grazing/common-index
relations, complete-table and wavelength refusal, exact-knot interpolation,
packet permutation, source/material fingerprint sensitivity, roughness bounds,
GGX visible-normal sampling/PDF agreement at grazing incidence and both pole
frames, directional-PDF mass versus sampler acceptance, reciprocity,
active-incident-medium dependence, a numerical no-energy-gain furnace, an
independent equal-expectation furnace check, equal-draw grazing variance against
the former NDF sampler, exact replay, and visibly distinct representative
tungsten/stainless renders under one neutral light. The Cornell golden remains
the legacy-material bit gate.

`lighting` inline tests and `tests/lighting_battery.rs` (feature `tracer`) cover
rectangle admission and identity order, zero/duplicate emitters, exact
solid-angle selection, constant and concentrated environment maps, spherical
PDF normalization, rotation/seam/pole conventions, supported linear-float EXR
ingestion and source/provenance binding, deterministic replay, mixed-emitter
MIS accounting, no-light refusal, cancellation, and finite high-dynamic-range
transport. Focused assertion messages include whichever seed, sample count,
light identities, PDFs, and radiance deltas are applicable to reproducing the
case.

`tests/tile_parallel_battery.rs` (feature `tracer`) covers odd and clipped tile
partitioning, one-tile worker admission, pre-dispatch film and scratch-envelope
memory refusal, exact serial equality across 1/2/4/8 workers and skewed
schedules, progressive partition changes, pre/mid-run cancellation, panic
containment, unchanged-film retry, and multiple bit-identical jobs through one
parked crew. Its owned `PendingRender` cases cancel after a strict committed-row
prefix, resume under a fresh authority without double-counting, remain below
two film payloads, refuse scratch before dispatch with zero progress, and retain
a precompleted zero-sample private job when publication starts under a
cancelled authority. A contained one-time worker panic also retains a strict
row prefix and resumes bit-exactly. Mode- or budget-changing retries refuse
before dispatch; cancellation/resume uses four workers and multiple tiles; and
a pool placement seed is not observable by scene charts. IID and zero-sample
jobs prove that they neither allocate nor charge unused Sobol direction state. Failures name
run, tile policy, pixel, channel, and binary64 bits; success logs setup,
traversal, compute, merge, idle, and memory measurements. The Euler cinematic
bridge adds animated-camera, geometry, dielectric, and direct-light
serial-versus-tiled equality through both one-shot and parked-crew paths.

`tests/aov_battery.rs` additionally covers complete cinematic AOV-film and EXR
byte equality against the serial oracle across 1/2/8 workers, both IID and
Owen-Sobol sampling, clipped tile layouts, a skewed schedule, and multiple AOV
profiles; two different jobs on one parked crew; pre-dispatch memory and AOV
configuration refusal; mid-trace cancellation; contained worker panic; and
exact denoise-guide/EXR plane equivalence including background signed zero.

`tests/adaptive_battery.rs` plus `tracer` inline tests (feature `tracer`) cover
policy admission and signed-zero canonicalization; fixed and truncated
checkpoint schedules; exact threshold equality and final-checkpoint precedence;
Welford variance/dispersion arithmetic, dark channels, one-noisy-channel
behavior, HDR offsets, power-of-two scaling, NaN/later-channel overflow/count
overflow rollback; constant-scene minimum stopping; exact path counts and
per-tile summaries; private AOV shape and identity getters; IID and Owen-Sobol
full-ceiling equality to uniform raw sums; exact decisions, moments, sums, and
counts across worker counts, tile shapes, skewed schedules, and parked crews;
and cancellation after committed rows followed by bit-exact parked resume. A
sparse GGX fixture records heterogeneous path allocation and compares adaptive
and rounded-up uniform cost against a disjoint-seed high-SPP reference. That
fixture is a regression for this scene/profile, not a universal quality claim.

`tests/checkpoint_battery.rs` plus `tracer::checkpoint` inline tests (feature
`tracer`) cover strict nonzero partial uniform and adaptive safe points;
bit-exact serialize/restore/finalize equivalence; preservation of raw sums,
moments, counts, and adaptive decisions; every-prefix truncation; body and seal
corruption; stale job, binding, execution-budget, and runtime-environment
refusal; well-sealed but semantically malformed tile/AOV state; one-byte-short
read/write ceilings; fallible and re-entrant sinks; cancellation before
emission and during the final seal; and no-retrace completion of already
finished uniform and adaptive jobs with the correct kernel identity. These
tests establish deterministic codec and resume semantics, not filesystem
durability or concurrent scheduler-claim arbitration.

`tests/render.rs` (7 cases): radical inverse known values; cosine samples are
unit vectors with the right pdf; the furnace test conserves energy exactly; MIS
weights sum to one (+ heuristic ordering); MIS integration is unbiased;
hero-wavelength integration exact on a constant / accurate on a ramp;
determinism.

`tests/metamorphic.rs`: relation
`lambertian-radiance-scale-equivariance` applies the shared shrinkable G3
unit-rescaling harness to production
`fs-render::Lambertian::furnace_radiance` (seed `0x2ACE_0006`, 384 cases,
absolute-or-relative tolerance `max_abs = 2e-12`,
`max_relative = 2e-12`). Each case draws albedo in `[0.05, 0.95)`, incident
radiance in `[0.125, 32)`, and a jointly shrinkable exponent in `[-3, 3]`; the
transform and expected output both use the exact power-of-two scale `2^e`.
This code-first declaration is proof-pending until the owning-crate batch suite
passes. The existing translated-scene frame-invariance pin in
`tests/charts.rs` remains independent.

`tests/diff_battery.rs` (bead qfx.5, feature `differentiable`): edge-aware
gradient vs central FD, the frozen-crossing negative control, a retained-fixture
quadrature-bias refinement regression, inverse rendering, a combined appearance/physics objective,
bitwise primal/gradient replay through the shared backend, a smooth-min seam
derivative regression, and fail-closed cancellation. Numerical receipts are
emitted by the current-tree run; this contract does not carry stale measurements
across backend-semantic changes. Every change to chart hit certification or
termination control flow must prove both the complete `tests/charts.rs` battery
and all nine `tests/diff_battery.rs` cases with the `differentiable` feature;
a charts-only proof set is narrower than the shared backend's blast radius.

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

- `GeometryInstance::try_new` validates that the supplied geometry identity is
  nonzero but does not derive or independently verify that identity against the
  chart or mesh bytes. Callers remain responsible for supplying the immutable
  geometry artifact's authoritative content identity. Instance transforms are
  visualization placements only: they do not alter mass properties, contact
  geometry, mechanics state, material parameters, or emission. Only proper
  rigid transforms are supported. Time-varying proper-rigid placements and
  motion blur are implemented under the contracts below; deformation, scaling,
  and time-varying emission or animated NEE-light sampling remain successor
  work. Emissive geometry may move with its instance, but its emission itself is
  time-invariant.
- v1 includes the scalar-BVH spectral path tracer. Wide-BVH SIMD traversal,
  watertight ray-triangle tests, a LIGHT-BVH for large emitter populations,
  heterogeneous volume/surface coupling, ray-stream sorting, and progressive
  tile streaming to HELM remain staged.
- The free `render_*_with_execution` functions and
  `render_cinematic_with_aovs_execution` are one-job convenience APIs and
  construct a scoped worker lane per call. Sustained frame sequences must use
  `RenderWorkerPool::with_parked_crew_local`. Timing fields are diagnostic
  wall-clock observations, not a universal scaling or 4K-attainment claim;
  bead `frankensim-h7xu5.5.5` owns representative 1080p/4K qualification.
- Transactional `&mut Film` compatibility calls necessarily retain the
  caller's committed film and one complete private staging film. Fresh
  `PendingRender` instead owns the one eventual film and atomic per-tile row
  prefixes; a suspension exposes counts and reports but never partial pixels.
  It binds the execution mode and compute budget at admission and refuses a
  retry under another mode or budget. Tracer-visible child contexts rebind the pool's placement seed to the
  public `Settings::seed`, so changing a parked pool's scheduling seed cannot
  change chart or camera semantics.
  `PendingAdaptiveRender` owns the corresponding raw sum, Welford, count, and
  decision AOVs with the same row-prefix rule. Both opaque states can emit and
  restore the versioned canonical checkpoint codec described below. `fs-render`
  itself still owns no filesystem, transaction, replacement, or scheduler-claim
  policy: crash recovery exists only when L6 streams a successfully sealed
  checkpoint through a transactional artifact store and re-admits the exact
  bound job. The codec also cannot certify identity against interior mutation
  inside a borrowed `Chart`; L6 must bind the immutable scene/source identities.
- Adaptive dispersion is an IID standard-error estimate only for IID streams;
  it is a within-stream heuristic for a single Owen-Sobol scramble. It is not a
  confidence interval under adaptive stopping, a formal image-error
  certificate, or evidence that a threshold transfers between scenes,
  exposures, samplers, materials, or output transforms. Statistically
  meaningful randomized-QMC error estimation requires independent scrambles.
  Denoised and postprocessed pixels cannot enter the accumulator or stopping
  API. Uniform rendering remains the final-quality fallback.
- Operation-memory receipts cover the named film, progress, Sobol, tile
  scratch, executor metadata, and arena charges. They do not cover OS thread
  stacks, allocator usable-size overhead, the small lighting-admission
  candidate/index structures created during preflight, or arbitrary heap owned
  by chart implementations. They are not process-RSS certificates. For an
  in-memory resumable job, lease requests, refusals, and peak usage are
  cumulative from job admission through the reported attempt; executor and
  timing fields are scoped to the most recent attempt.
- Environment radiance uses a deterministic bounded linear-sRGB spectral lift
  and piecewise-constant lat-long texels. It does not claim measured spectra,
  arbitrary HDR decoder compatibility, texture filtering, sun/sky delta
  models, portal sampling, or calibrated real-studio illumination. The current
  exact-CDF light selector is intended for the small studio rigs in scope; it
  is not the planned many-thousand-emitter light BVH. The standalone
  `AdmittedLighting::try_new` compatibility/analysis constructor has no `Cx`;
  production tracer entry points use its cancellation-aware counterpart.
- Dielectric support is homogeneous, non-polarizing geometric optics. It does
  not claim polarization, coherence, fluorescence, birefringence, measured
  preset fidelity, camera-start-inside-medium support, arbitrary overlapping
  media, or a topology certificate. GGX transmission is single-scattering and
  does not claim multiple-scattering furnace closure at appreciable roughness.
  Shadow rays stop at the first intervening surface rather than travelling
  undeviated through glass. Difficult focused caustics therefore remain a
  slow-convergence case for this unidirectional integrator; no
  bidirectional/manifold transport, denoising, or unbiased firefly-clamping
  claim is made. Fixed-metric, representability-aware ray offsets are robust
  engineering, not a certified positional-error enclosure.
- Conductor support is homogeneous, opaque, non-polarizing geometric optics
  with isotropic single-scattering GGX. The representative presets do not claim
  material identification, alloy/lot or product calibration, a measured finish,
  oxide/passivation or thin-film layers, anisotropic brushing/machining grooves,
  temperature dependence, polarization, interfacet multiple scattering, or a
  conversion from measured `Ra`/`Rq` to GGX alpha. Caller-asserted measured
  source metadata is retained but not independently authenticated. A content
  hash establishes exact replay inputs, not physical truth.
- The tracer implements spectrum→CIE XYZ→Bradford-adapted linear sRGB and
  floating-point EXR output. Display-referred tone/color management and layered
  measured-spectrum materials remain staged.
- `mis_integrate_unit` is a 1-D demonstrator of the balance heuristic; the
  production MIS lives in the path integrator across BSDF/light strategies.
- The G3 furnace adopter covers only bounded positive scalar albedo/radiance,
  exact power-of-two unit rescaling, and the `16`-sample furnace call. It does
  not establish a runtime physical-unit system, arbitrary-scale conditioning,
  Monte Carlo error bounds, spectral/path-integrator rescaling, or the separate
  translated-scene frame-invariance claim. Its declaration is not evidence of
  a passing batch run.

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
  short opposite-sign witness above, it can authorize only another certified
  root-free safe step; it never mints a `Hit`. If bounded continuation cannot
  reach such evidence, the production trace returns `ResidualLimit`.
  An independent oracle may therefore locate a root beyond the bounded
  continuation budget while the production trace honestly remains
  `ResidualLimit`; that is unresolved completeness, not a certified miss.
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

## Motion-time foundation

`motion` defines a timed-ray envelope for every existing spatial ray backend.
Frame shutters resolve centered, front-loaded, or back-loaded exposure
intervals inside explicit shot bounds. Normalized shutter coordinates are
finite and in `[0,1]`; `TimedRay` retains the complete admitted shutter,
normalized coordinate, and absolute seconds so a downstream adapter can reject
mixing rays from coincident but different exposures. `UniformCounterV1` and
`StratifiedCounterV1` use a dedicated stable counter domain, independently mix
the explicit render-stream seed, logical pixel, and absolute sample identities,
and never consume the tracer's wavelength, pixel, light, lens, or BSDF streams.
Every consecutive full-cycle window of a stratified stream visits each temporal stratum once,
while the keyed permutation spreads incomplete groups across the exposure
instead of pinning low sample IDs near shutter open. Zero-width shutters map
every sample bit-exactly to the static frame time. Positive requested durations
that collapse to one binary64 endpoint at the declared absolute-time scale
refuse instead of silently becoming static; explicit normalized endpoints
return the stored open/close bits exactly.

`animated_instances` binds shared immutable geometry to strictly ordered
proper-rigid keyframes. Translation uses cubic Hermite interpolation with
declared endpoint velocities; orientation uses shortest-arc SLERP with a
normalized-linear fallback for nearly coincident quaternions.
Evaluation is by a `TimedRay`'s absolute time, never extrapolates, retains object
and geometry identities, and binds the returned frame identity to the evaluated
pose. It polls cancellation before and after trajectory materialization and
backend intersection. Equal camera/object translations preserve the relative
hit, while holding the camera fixed changes it as expected.

`motion_bounds` validates finite local boxes and returns conservative finite
world AABBs. Arbitrary proper rotation is enclosed without angular sampling by
an outward-rounded body-origin sphere, so endpoint quaternions cannot hide
additional full spins. Linear segments retain their declared endpoint
translation envelope and constant identity rotation uses the tighter analytic
translated-box envelope. The runtime-trajectory path converts every overlapping
cubic-Hermite translation segment to its equivalent Bezier controls and bounds
their convex hull, including velocity-driven interior overshoot; it keeps
shutter-clipped boundary segments whole and adds an outward numerical margin.
This is a safe primitive for the later animated TLAS, not a TLAS itself.

`temporal_accumulation` evaluates absolute logical sample IDs in fixed order and
transactionally appends three-channel linear RGB or XYZ sums. A checkpoint
retains shutter, shutter-stream identity, pixel identity, color interpretation,
sum, accepted count, and next sample ID. One-shot and contiguous progressive
partitions execute the same floating additions and are bit-identical.
Cancellation, non-finite evaluator output, range mismatch, and arithmetic
overflow leave accepted state unchanged.

With feature `tracer`, `Shape::AnimatedInstance` and
`render_motion`/`render_motion_range` make the time contract part of actual
spectral rendering. One shutter time keyed by render seed, pixel, and absolute
logical sample is drawn per camera path and retained for every secondary and
shadow ray. A progressive `Film` checkpoints the exact shutter plus shutter
stream seed and refuses mixed exposure histories. Animated geometry on the
legacy static entry point refuses with `MissingRayTime`; a shutter outside any
animated trajectory refuses before film publication. Legacy static rendering
does not draw the motion dimension, and rendering static geometry through the
motion entry point is bit-identical in XYZ output. A zero-width shutter at an animated
keyframe is bit-identical to the equivalent static instance. G0/G3/G5 coverage
includes malformed/full/zero exposures, deterministic strata and progressive
partitions, trajectory admission and interpolation, moving-camera/object
relative intersection, conservative sampled bounds, transactional
cancellation, and high-rate spectral renders matching analytic
constant-velocity and constant-spin occupancy envelopes.

Rectangular NEE lights are static metadata. A scene that names an animated
instance as any such light refuses with `AnimatedLightUnsupported`; it is never
rendered with a stale light-sampling transform. Environment maps are immutable
linear-sRGB radiance artifacts in a declared Y-up lat-long layout. Their
semantic identity binds dimensions, canonical f32 pixel bits, color/layout
interpretation, and rotation; separate source and provenance identities retain
the supported-EXR import lineage.

## Cinematic physical camera

`camera` adds an opt-in, validated ideal thin-lens camera without changing the
legacy `tracer::Camera` or its frozen image stream. `CameraProjection` accepts
either focal length plus **vertical** active sensor height in metres, an
explicit vertical FOV in radians, or the legacy vertical half tangent. Sensor
format is never guessed: in particular, there is no hidden 36x24 mm default.
The render dimensions still declare film aspect ratio.

`PhysicalCamera` admits a finite eye, a positive focus distance, and an
orthonormal right/up/forward frame. Look-at construction uses scaled
normalization so representable very large or small vectors do not overflow or
underflow before validation. An explicit up reference whose scale-free sine
against the view axis is below the contract threshold refuses with ranked
fixes; the camera never silently chooses a different roll. Camera local axes
are the proper rotation `+X=right`, `+Y=up`, `-Z=forward`.

Focus distance is the positive axial distance from the lens plane. For raster
vector `v = forward + x*right + y*up`, all thin-lens rays converge on
`eye + focus_distance*v`. A declarative world-point focus track projects the
identity-resolved point onto the evaluated optical axis; a target on or behind
the lens refuses. It is not heuristic autofocus.

`Aperture` is constructible only through validating constructors. A circular
aperture uses the area-preserving concentric-square map. A regular bladed
aperture declares its circumradius and canonical rotation, precomputes 3--64
vertices, and samples equal-area centre triangles uniformly; its closing
triangle reuses vertex zero bit-exactly. Radius zero canonicalizes to a true
pinhole. `Aperture::try_from_f_number` uses
`radius = focal_length / (2*f_number)`.

`CameraShot` linearly interpolates eye position, uses deterministic
shortest-arc quaternion SLERP for orientation, and linearly interpolates either
axial focus distance or an explicitly moving world focus point. It returns
stored endpoint bits exactly and holds the first/last keyframe only inside the
declared shot bounds. Projection, aperture, exposure metadata, and focus policy
may not change inside one continuous shot. `AnimatedCamera` never extrapolates
or blends across hard cuts. A positive-width shutter must belong wholly to one
shot; a zero-width exposure exactly on a cut declares `CutSide::Before` or
`CutSide::After`. The admitted `CameraExposure` binds subsequent evaluation to
that shot so a rare exact boundary sample cannot jump across the cut.

With feature `tracer`, `render_cinematic`/`render_cinematic_range` evaluate the
camera and animated geometry at the same one absolute `PathTime`. Lens U/V use
the separately versioned Philox camera-lens domain and never consume the
pixel/wavelength, shutter, light, or BSDF streams. The pinhole branch performs
the legacy camera arithmetic directly rather than reconstructing a focus point,
so an equivalent aperture-zero camera is bit-identical to the legacy render.
Camera validation, shutter admission, evaluation, and ray generation poll
`Cx`; failures leave a progressive film transactionally unchanged.

Determinism is same-ISA/toolchain bit replay under the existing tracer policy.
The camera model makes no claim for lens distortion, chromatic aberration,
diffraction, rolling shutter, vignetting, autofocus, sensor irradiance, or
photometrically calibrated exposure. `ExposureMetadata` is retained but is not
applied by the current radiance-averaging integrator. `FilmTimeMode` binds the
render-path kind, shutter, stream seed, and cinematic shot ID, so progressive
appends cannot cross a cut or switch between legacy-motion and cinematic rays.
Like the pre-existing progressive tracer, it does not yet content-bind the
scene or full camera content; composition-level artifact identity owns that
broader provenance closure.

The motion path reconstructs only the supplied interpolation model. It does not
invent mechanical bandwidth, unwrap rotations absent from keyframes, smooth a
declared contact/terminal discontinuity, validate the source dynamics, or claim
that motion blur recovers unresolved terminal physics. Event-aware subdivision
or refusal is owned by the Euler trajectory adapter; final animated TLAS and
cinematic-camera composition remain their dependent beads.

## Stable hit correspondence and motion vectors

`motion_vectors` (bead `frankensim-h7xu5.3.5`) retains the original triangle
index, barycentric coordinates, accepted local point, and distinct local
geometric/shading normals from the same instance hit used for shading. Stable
hit identity binds the nonzero object ID, caller-supplied immutable-geometry
identity, deterministic material identity, and original mesh triangle. The
low-level constructor accepts an explicit nonzero material digest; the tracer's
aligned cinematic-sample path derives that digest from the complete `Material`
value (excluding separately modeled emission). Pose or frame identity is
deliberately separate because it changes while the same material feature moves.
The geometry digest remains a caller assertion; this module does not re-hash
mutable public mesh storage.

Motion evaluates one retained local mesh point under explicit previous,
current, and next proper-rigid object transforms and physical cameras. The
current frame must match the accepted hit's frame identity, reference times are
nondecreasing, and every frame must name the same object and geometry. Camera
projection uses the deterministic optical centre, independent of thin-lens
samples. NDC is `+x` right/`+y` up; row-major pixel coordinates are `+x`
right/`+y` down. Each displacement is explicitly `target - current`, both in
NDC per declared frame step and pixels per declared frame step. It is neither
metres per second nor an observed mechanical velocity.

The raster extent is explicit and uses the same half-open pixel-edge convention
as the tracer, including 3840x2160 without a shape-specific path. Off-screen
in-front projections retain their finite displacement. Points on or behind the
lens plane have no vector. Different continuous-shot IDs produce `CameraCut`
rather than a cross-cut numeric vector. Generic chart hits return
`ChartHasNoStableParameter`; a local implicit hit point is not silently promoted
to a stable material coordinate. The Euler cinematic bridge's current
axisymmetric surface is a deterministic immutable preview mesh, so it uses the
mesh correspondence path; this does not promote that chordal visualization
mesh into mechanics authority.

Projection alone makes no visibility claim. `validate_reprojection` compares a
frontmost target observation using full categorical identity, axial depth, and
mesh barycentrics under caller-declared finite tolerances. It distinguishes
visible, nearer occluder, absent/farther disocclusion, hard cut, off-screen,
behind-camera, identity mismatch, depth mismatch, topology ambiguity, and a
missing local witness. Target observations must come from an aligned AOV or an
equivalent target-frame trace; the motion module does not fabricate them.

`trace_cinematic_pixel_sample` is the opt-in lossless producer seam: it returns
the unaccumulated beauty XYZ, exact sampled shutter time, and primary record
from one traversal. Ordinary beauty rendering does not compute or hash this
record when it will be discarded. The dependent film/AOV layer can therefore
accumulate aligned records without retracing and without imposing the feature's
cost on legacy renders.

G0/G3/G5 coverage pins analytic translation, rigid rotation, camera pan,
common camera/object motion, quaternion double-cover equivalence, aperture
independence, hard cuts, off-screen and behind-camera endpoints, mesh-edge tie
ownership, chart refusal, reprojection classifications, deterministic replay,
and exact 4K raster indexing. Aligned film buffers, checkpoint payloads, and EXR
channels are intentionally owned by dependent bead `frankensim-h7xu5.6.1` and
are not claimed by this kernel alone.
