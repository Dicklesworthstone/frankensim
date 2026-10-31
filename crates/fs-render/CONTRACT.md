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
- `spectral::SpectralReflectance` and `Material::SpectralDiffuse` retain a
  declared passive visible-band reflectance curve. The native diffuse BSDF
  evaluates `rho(lambda)/pi` at every wavelength, including exact black/white,
  without assigning a temperature or an emission source. See the calibrated
  diffuse-material contract below for table admission and accuracy limits.
- `balance_heuristic` / `power_heuristic` — MIS weights; `mis_weight_sum(pf,
  pg)` — the weight-sum audit (nominally `1`).
- `mis_integrate_unit(f, n)` — an unbiased MIS estimate of `∫₀¹ f` combining
  uniform + linear-importance strategies.
- `hero_wavelengths(hero, count, min, max)` /
  `hero_wavelength_packet::<N>(hero, min, max)` /
  `spectral_integral(spectrum, min, max, samples)` — heap-backed and fixed-size
  hero-wavelength spectral integration surfaces.
- `dielectric` module (feature `tracer`): validated Cauchy and three-oscillator
  Sellmeier laws share `DielectricIor`, in vacuum nanometres. Sellmeier strengths
  and squared resonance wavelengths are nonnegative; active resonances inside
  the inclusive 380--780 nm band refuse. Monotonicity between poles makes the
  endpoint admission meaningful across the whole transparent band, with
  checked binary64 evaluation rather than an interval certificate. Both laws
  expose the narrowband group index `n - lambda dn/dlambda`; every existing
  smooth/rough, nested-medium, direct-connector and bidirectional transport path
  continues to use the phase index. The new law is allocation-free and requires
  three oscillator evaluations per wavelength, with validation at construction.
  Existing Cauchy material/scene identities remain unchanged; Sellmeier identities
  bind a distinct law tag and all six canonical coefficients.
  Coefficients must describe absolute indices versus vacuum wavelength;
  air-relative catalog data need an explicit conversion/refit before physical
  accuracy is claimed. There is no built-in measured-glass preset or claim for
  anomalous/absorptive resonance, transient propagation or pulse broadening.
  `TabulatedIor` extends the same phase-law owner with nine ordered absolute
  index samples spanning exactly 380--780 nm, with gaps at most 100 nm and
  bounded linear interpolation. Every knot is retained in material identity.
  Its group-index diagnostic uses the declared interval slope and refuses
  unequal-slope interior knots, where that interpolant has no derivative.
  Spectral polarized Fresnel channels and phase-only waveplates accept the
  same Cauchy, Sellmeier and tabulated laws, converting once per spectral loop.
  Homogeneous Beer-Lambert absorption remains in inverse metres, with
  provenance-bearing representative glass presets, exact unpolarized Fresnel,
  Snell refraction and total internal reflection, plus smooth-delta and
  single-scattering isotropic-GGX reflection/transmission. Radiance transport
  uses the documented `(eta_i/eta_t)^2` transmission convention.
  Native spectral bulk absorption retains 2–16 declared vacuum-wavelength
  knots covering exactly 380–780 nm, with nonnegative intensity coefficients
  in inverse metres and deterministic linear interpolation. Imaginary-index
  inputs convert knotwise as `alpha = 4*pi*k/lambda_m`, retaining the stated
  linear-alpha approximation and weakly absorbing dielectric interface model.
  No RGB fit, opacity floor, or extrapolation is introduced. Thin-path absorbed
  fractions use `expm1`, and material identities retain every admitted knot.
- **Dielectric launch media:** camera pupils, cinematic/moving origins, area-only
  photometric probes, and both BDPT endpoints initialize the physical host
  before their first segment. Mesh edge closure/orientation is prepared once;
  local support bounds prune queries, and compensated oriented solid-angle
  winding determines mesh membership. Charts supply a finite certified strict
  sign. Static legacy eyes cache their admitted boundary IDs; aperture and
  moving origins use the same immutable prepared topology in their own frame.
  Multiple containing bodies require a conservative nesting witness: the inner
  transformed support lies in all outer mesh inward half-spaces, or in a
  certified negative chart step ball. Unproven concave nesting, overlaps,
  reversed/open/unresolved winding and origins on a boundary refuse rather than
  choosing the host from a first back-facing hit. These checks are not a global
  self-intersection/topology certificate. Pure emitter materials do not invent
  enclosing volumes. The original LIFO stack then supplies absorption, nested
  restoration, incident index, total internal reflection, participating-field
  selection and endpoint connection attenuation. Exterior origins retain their
  existing random streams and radiance arithmetic. Calibrated point-fixture
  probes still require their separate refractive direct connector; area-only
  embedded probes use host absorption in both NEE and continuation paths.
- `gradient_index` module (feature `tracer`): allocation-free adaptive
  Dormand--Prince 5(4) propagation of the geometric-optics ray equation
  `d(n t)/ds = grad(n)`. An immutable field supplies local phase index,
  analytic spatial gradient, true power absorption and optional group index.
  Position, direction, parallel-transported transverse frame, phase/group
  optical lengths and absorption depth evolve along physical arc length.
  Explicit local tolerances, metric step bounds and a hard attempted-step cap
  govern execution; cancellation is polled before every attempted step,
  including rejections. Local estimates are not global error certificates.
  Native homogeneous spectral glass has a one-sample analytic fast path;
  a bounded affine field retains the same dispersion/absorption owner.
  Missing or nonpositive tabulated group derivatives suppress the delay claim
  while retaining supported phase propagation. `PolarizedBeam` consumes the
  propagator directly: all Stokes power components receive integrated
  Beer--Lambert attenuation in the transported frame, and subsequent native
  smooth/rough interfaces use the endpoint's actual phase index. Its
  `BeamMedium` distinguishes homogeneous laws from position/wavelength-bound
  field samples, refusing straight advancement of a gradient medium and
  mismatched local state. Calibrated fixtures can emit into either medium.
  `gradient_index::event` connects the same adaptive state to a bounded
  clearance/proximity oracle. A true free-ball radius bounds traveled arc
  length; every internal Runge--Kutta stage must remain inside it before the
  field is sampled. Upper surface proximity and a transverse approach are
  separate contact requirements. Numerical contacts retain their resolution
  and integrated endpoint; they are not certified global crossing enclosures.
  An admitted local planar patch permits positive-length departure, including
  its optical lengths and attenuation, without ignoring a whole object or
  skipping its opposite wall. The selected field must support the explicit
  numerical contact collar. Unsupported patches, grazing ambiguity and
  exhausted step/query budgets refuse. Spatially uniform fields reuse one
  optical sample and take analytic steps limited by geometry. The caller
  retains the same accepted field and geometry snapshots.
  No scalar camera/BDPT curved visibility, ray-bundle Jacobians,
  scattering, diffraction, birefringence or full-wave accuracy is claimed.
  G1 checks include analytic affine-index trajectories, Snell invariants,
  phase/group separation, absorption and tolerance refinement; G3 checks
  cover reverse paths, transported frames and composed segments; G4 checks
  include cancellation during rejected adaptive trials.
- `conductor` module (feature `tracer`): fixed-capacity visible-band tables of
  absolute complex refractive index `eta(lambda) + i k(lambda)` in vacuum
  nanometres, strict inclusive 380--780 nm support, deterministic linear
  interpolation, exact unpolarized complex-IOR Fresnel against the path's active
  incident medium, and validated isotropic-GGX roughness. Built-in tungsten and
  stainless-steel tables are explicitly representative look-development data;
  table/source/material identities prove replay binding, not specimen fidelity.
  `ConductorOptics::with_coating(CoherentDielectricCoating)` additionally admits
  one explicitly source-bound transparent isotropic film: the shared visible
  Cauchy, Sellmeier or tabulated phase-index law and a uniform thickness in metres, from
  zero through 100 micrometres. At each wavelength the same phase index sets
  both the interface admittance and the coherent optical thickness; group
  index is never substituted. Existing Cauchy coating identities retain their
  original bytes; Sellmeier binds a separate law tag and all six coefficients,
  while tabulated laws bind each wavelength/index pair. Separate complex s/p amplitudes retain substrate reflection
  phase and sum every internal reflection before unpolarized power averaging.
  A fixed-size characteristic-matrix evaluation handles the film's critical
  angle by its analytic sinc limit and scales evanescent transfer entries to
  preserve frustrated total internal reflection without exponential overflow.
  Both isotropic and brushed conductor BSDFs, including their existing direct,
  camera, photometric and bidirectional consumers, use this same Fresnel owner.
  The film follows the local microfacet normal; geometric proposals remain
  unchanged. Bare presets and bare optical identities remain unchanged. A
  coating is never inferred from a metal name or material-state binding.
  `ConductorOptics::with_absorbing_coating(CoherentAbsorbingCoating)` admits an
  explicitly sourced passive film with nine absolute `eta + i k` knots over
  the same visible support. `CoatingIorSample` accepts exactly zero extinction;
  an opaque substrate continues to require positive extinction. Both index
  components interpolate linearly without extrapolation. Complex extinction
  participates in the interface reflection phase and propagation before power
  averaging; it is not a post-interference tint. Complex propagation entries
  share an `exp(-Im(delta))` scale, so highly absorbing and evanescent films
  never require growing exponentials. An analytic critical-angle series retains
  the transparent and weak-loss limits. Zero thickness, a film matching the
  substrate, or a transparent film matching the incident index retain bare
  Fresnel arithmetic. Each attachment replaces the existing film; it does not
  silently create a multilayer. Storage is one fixed spectral table, evaluation
  is bounded and allocation-free, and all native conductor BSDF consumers use
  the same owner. Opaque substrate absorption remains implicit in `1-R`; this
  interface does not resolve heat deposition between film and substrate.
  The caller declares sufficient coherence and locally parallel, laterally
  uniform boundaries. This is local single-film interference, with no arbitrary
  multilayer stack, thickness distribution, polarization
  transport between bounces, diffraction, or scene-wide coherent field claim.
  G0/G1/G3 tests cover independent complex Airy answers, quarter-wave
  cancellation, half-wave/zero-thickness/bare limits, incident-medium effects,
  critical/evanescent continuity, passivity and both native conductor BSDFs.
  Sellmeier regressions additionally cover independent Airy spectra, constant-law
  bit equivalence, sensitivity to all six coefficients and frozen Cauchy identity.
  Absorbing-film regressions cover independent complex Airy answers, weak loss
  near a destructive null, critical/evanescent transparent limits, opaque-film
  screening of the substrate, retained spectral knots, passivity across optical
  depth, and both native conductor BSDFs. A supplied visible-band table does not
  establish a broadband causal dispersion law or measured specimen fidelity.
  The derivation follows [Byrnes, *Multilayer optical calculations*,
  sections 2--5 and appendix D](https://arxiv.org/abs/1603.02720).

- `material_state` (opt-in `material-state`, including `tracer` and L3
  `fs-material`/`fs-conduction`): `MaterialOpticalBinding` converts resolved visible conductor
  or dielectric data directly to the existing tracer `Material`. Mechanics and
  optics must share one exact card and query point; mechanical requirements
  remain consumer-owned. Surface identity and roughness are independent inputs.
  Cauchy coefficients convert from SI to micrometre powers, while wavelength
  evaluation uses nanometres and absorption distance uses metres. No fallback
  preset or chemistry-name dispatch is used. Admission is synchronous, bounded
  to nine optical knots, safe Rust, and does not mutate accepted physical or
  scene state. Underlying optical model refusals remain structured errors.
  Binding identities establish replay consistency, not measurement authority;
  conductor source status is `MaterialStateClaims`, preserving query identity
  without labeling authored or unmeasured inputs as measured.
  `try_tabulated_dielectric` and the generic family selector additionally bind
  eighteen material-state-selected visible phase-index/absorption samples to
  the existing dielectric transport. No Cauchy fit or RGB conversion is used
  on that path; wavelengths are vacuum nm and intensity absorption is 1/m.
  State-dependent dispersion and absorption therefore reach Snell refraction
  and Beer-Lambert attenuation directly, with source identities preserved.
  Dielectrics retain native nine-knot intensity absorption with typed
  power-nepers/metre semantics;
  the separately selected legacy RGB schema retains its lifting approximation.
  The material resolver refuses ambiguous or incomplete absorption schemas.
  `MaterialOpticalBinding::at_material_state` refreshes the exact original
  optical source to a new accepted state's complete condition through that
  shared pinned-claim owner. It retains the admitted optical family, smooth or
  GGX finish, surface identity and spectral-versus-RGB absorption schema, and
  returns both the native material binding and refreshed source for subsequent
  endpoints. Stale or substituted source states, foreign cards, changed axis
  meanings and unsupported original claims refuse before either input changes.
  The existing conductor Fresnel and dielectric Beer-Lambert owners consume the
  refreshed values; no temperature-color curve or new optical law is introduced.
  Focused G0/G1/G3 material-state tests cover those refusals, continued spectral
  transmittance, independent complex-index Fresnel and actual reflected pixels.
  This optical binding adds no thermal
  emission, birefringence, phase evolution, or measured-spectrum validation.
  Focused G0/G1 regressions in the
  existing Euler scene-adapter unit seam exercise the generic optical families,
  minimal string versus disc mechanical requirements, Fresnel and absorption
  formulas, SI conversion, material replacement, and mismatch refusals.

- `thermal_surface` (feature `tracer`): `GrayDiffuseThermalSurface` declares an
  opaque gray/diffuse local-equilibrium approximation explicitly. Its source,
  temperature and emissivity identity drives SI Planck spectral radiance in
  W/(m² sr nm), with reciprocal diffuse reflectance `1-epsilon`. This material
  emits on camera-path hits and continues reflected-light transport; it is not
  an RGB temperature-color curve or a second debit to thermal storage.
  `MaterialThermalBinding::try_gray_diffuse` retains the original
  `fs-conduction::radiation::SurfaceEmissivity` sample and requires its exact
  card, complete query point and coordinate quantity meanings to match an
  accepted `ResolvedMaterialStatePoint`. The physical temperature comes from
  the radiation owner, including explicitly typed coordinates not named `T`;
  pressure and phase conditions cannot be dropped when binding appearance.
  Changed thermal endpoints require rebinding. The declaration does not infer
  spectral/directional emissivity from a total measurement or visible n+k,
  and replaces rather than adds to the independent conductor/dielectric law.
  Like the existing opaque material convention, each side emits into its
  outward hemisphere; `epsilon*sigma*T^4` is power per emitting side's area.
  Emission into an active dielectric medium refuses: this law describes a
  vacuum/exterior-air side, with no inferred refractive-index correction.
  The observer integrates the visible band under the existing normalized CIE
  convention; conduction retains the full bolometric cooling budget. Actual
  static thermal mesh faces, including rigid mesh instances and accepted body
  boundaries, participate in NEE with their original Planck spectra. Source
  selection caches visible CIE luminance per thermal state and weights each
  triangle by its subtended solid angle; the conditional proposal is exact
  uniform area. Forward and reverse mixture PDFs include all overlapping
  same-primitive triangle proposals. Thermal hits use the corresponding MIS
  weight and continue through the reciprocal reflection law. The last admitted
  bounce gives NEE full weight because no competing continuation is evaluated.
  Charts and paths through smooth dielectric
  connectors retain camera-path emission without a fictitious NEE competitor.
  Separate `Primitive::emission` is refused for this material. BDPT refuses it
  before film publication until a physical thermal
  light-subpath proposal exists. Focused regressions cover the Planck integral
  against Stefan--Boltzmann, emitted pixel intensity/color, retained reflection,
  original cooling flux, material swaps and mismatched-state refusals. The
  thermal lighting battery also compares diffuse illumination against independent
  area quadrature under NEE, BSDF, balance MIS and power MIS, and verifies
  occlusion, emissivity scaling and construction-order/progressive replay. These
  are numerical consistency checks, not experimental spectral validation.

- `SpectralDiffuseThermalSurface` extends the same opaque diffuse transport
  with an explicit visible-band Kirchhoff law. `SpectralEmissivity` admits
  2–16 strictly increasing linear knots spanning exactly 380–780 nm with
  emissivity in `[0,1]`; evaluation never extrapolates. Each wavelength emits
  `epsilon(lambda) * B_lambda(T)` and reflects `1-epsilon(lambda)`. Source
  identity binds the actual numerical curve, temperature and declared source.
  A visible table does not supply bolometric exitance, infrared continuation,
  a cooling debit or authenticated specimen measurements. The original gray
  model and its all-wavelength source arithmetic remain available unchanged.
  Camera paths, opaque BSDFs, stationary and moving thermal mesh proposals,
  physical-medium phase lighting and accepted scene exposures use the same
  spectral law. AOV albedo caches a separate visible reflectance projection.
  Prepared triangles and light samples borrow a compact source view; they do
  not copy a 16-knot table per triangle or path. Visible proposal quadrature
  splits at the emissivity knots and retains positive narrow-band support
  even when the old fixed 5 nm nodes miss it. Physical transport still evaluates
  the original curve at each path wavelength. BDPT and the static photometric
  scene adapter retain their explicit thermal-source limitations. Numerical
  tests cover independent SI Planck/Kirchhoff equilibrium, gray limits,
  narrow-band proposals, source PDFs, direct-estimator energy agreement,
  moving accepted source states and deterministic execution replay.

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
  bit-semantics v11. Mesh BVH traversal accepts an inclusive scene-level
  current-best distance, transforms it outward into the normalized-ray
  parameterization, and prunes farther nodes without changing the selected hit
  or exact-distance tie rule. `Hit.steps` remains an honest count of work
  actually performed and can therefore decrease under this pruning.
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
  `TriMesh` uses Woop/Benthin/Wald shear-space intersections over a
  deterministic median-split BVH with outward-rounded slab pruning;
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

- `tracer` finite-light bidirectional transport (opt-in, [F]):
  `render_cinematic_bidirectional` constructs paired camera and rectangle-light
  subpaths, evaluates every admissible `(s,t)` connection through the declared
  depth, and combines reciprocal area-measure densities with power-heuristic
  path-space MIS. Radiance and importance walks use their respective
  refractive transport factors. Two-sided rectangle emission has identical
  support and radiance in emitter-hit, direct-light, and launch strategies.
  Each logical sample evaluates the established randomized four-wavelength
  stratification and averages its lanes, retaining unbiased spectral
  integration while reducing single-hero chromatic variance. Spectral lane is
  part of the canonical identity of every cross-pixel splat.
  The first smooth camera-side dielectric interface is evaluated as the exact
  sum of its Fresnel reflection and transmission terms; later interfaces stay
  sampled, bounding live suffixes while removing the dominant plate-entry
  roulette. Camera-independent `t=1` light splats are evaluated once, not once
  per enumerated camera suffix.
  For an all-specular camera prefix followed by a connectible surface, finite
  rectangle lighting through one admitted parallel smooth slab reuses the
  exact two-interface Snell/Jacobian connector. That proposal and the matching
  camera-walk emitter hit carry complementary power-heuristic weights. General
  manifold-path MIS with a non-delta camera prefix is not claimed and does not
  silently inherit this restricted proposal.
  `render_cinematic_bidirectional_with_execution` and the parked-crew method
  partition work into fixed eight-sample full-raster blocks because `t=1`
  light paths may splat into arbitrary pixels; blocks publish only after full
  executor drain and reduce in ascending sample order, so worker count and
  arrival order do not affect bits. Owen-Sobol full-path mode assigns disjoint
  nested-scramble domains to camera walks, light walks, and each connector.
  Current support is a pinhole cinematic camera, finite rectangle emitters,
  opaque materials, and smooth/rough dielectrics. Environment illumination
  remains a separate camera-path/NEE estimator rather than a fictitious finite
  launch surface. Finite apertures, adaptive BDPT sampling, aligned BDPT AOVs,
  vertex merging, and participating-media bidirectional transport are explicit
  no-claim/refusal boundaries.

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

- `tracer` independent-pilot fixed production allocation: an
  `IndependentPilotSamplePlan` converts a completed pilot film's per-channel
  sample variances into a fixed row-major production count map using an
  explicit absolute/relative error scale, dark floor, safety factor, and
  minimum/maximum interval. The pilot and production counter-stream seeds must
  differ, pilot radiance is discarded, and every production pixel consumes its
  complete preselected prefix before any retained production sample is
  inspected. Conditional on the plan, this is therefore the ordinary
  fixed-count path estimator rather than data-dependent stopping. The
  allocation remains a quality/cost heuristic: a single Owen scramble's
  within-stream variance is not a confidence interval. Version 1 supports a
  fresh cinematic beauty film, exactly co-sampled variable-count cinematic
  AOVs, count-aware denoising guides, and linear-float EXR carrying the count
  map and pilot-plan provenance. It does not yet claim checkpoint/resume or
  shard merge for a fixed-plan production pass.

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

  `CinematicAovFilm::denoise_guides` exposes the retained guide subset as owned
  row-major planes: previous-frame motion, axial depth, shading normal, primary
  coverage, and raw-luminance variance. `FinalDiagnostic` additionally exposes
  the exact nearest-primary object/material palette indices as optional `u64`
  categorical labels; `DailyCore` leaves those labels unavailable. The palette
  and scene binding make the indices sequence-stable only while that exact
  binding is preserved. The bridge uses the same averaging, normal
  normalization, finite narrowing, invalid background-zero convention, and
  signed-zero canonicalization as EXR export; the focused battery compares
  every returned plane bit-for-bit with its decoded EXR counterpart.
  `BeautyOnly` refuses because it retained no primary observations.

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

  Fresh uniform and adaptive AOV films have explicit deterministic tile paths:
  `render_cinematic_with_aovs_execution`,
  `render_cinematic_adaptive_with_aovs_execution`, and their corresponding
  `ParkedRenderScope` methods. Each logical tile owns disjoint pixels and
  accumulates every pixel in ascending absolute-sample order before copying
  the complete tile into one private full-frame staging film. The adaptive path
  co-stages its raw XYZ Welford state and AOV observations, and retains only
  the exact `0..terminal_count` path prefix selected for each pixel; AOV state
  cannot influence the stopping decision. The paths reuse
  `RenderExecutionConfig`, the ambient `Cx`, the existing operation lease, and
  the ordinary/parked `RenderPoolRunner`; worker count, stealing, and tile
  completion order cannot change estimator, decision, or AOV arithmetic. The
  fresh paths allocate private staging rather than constructing an empty film
  and cloning it as the progressive compatibility API must. Cancellation, AOV
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

- `volumes::transport` (features `volumes,tracer`) connects borrowed simulation
  density and temperature snapshots to physical scene transport. `AbsorbingGrid`
  requires kg/m³ density, metre geometry and an explicit source-domain mass
  absorption spectrum in m²/kg; raw lattice density is not silently interpreted
  as a physical optical coefficient. Each crossed constant cell integrates
  Beer--Lambert attenuation and, when a kelvin field is supplied, local-equilibrium
  emission `sigma_a * B_lambda(T)` using the same SI Planck owner as thermal
  surfaces. Front-to-back emission has correct intervening absorption, stable
  optically thin `expm1` evaluation and a repeated-temperature packet cache.
  Admission validates and hashes borrowed fields once. Ray traversal allocates
  nothing, visits intersected cells through DDA, checks cancellation and enforces
  the caller's cell-scan/segment-step limits. Shadow transport skips Planck
  evaluation and performs one final exponential per wavelength.
  `render_absorbing` and `render_absorbing_range` use the existing surface
  integrator, including its spectral BSDFs, lights and counter-based sampling;
  camera, reflected, RGB/thermal direct-light and environment segments all
  consume the same medium. This is a static-camera, absorption/emission,
  refractive-index-one ambient model. Closed dielectric bodies displace the
  grid and retain their material's own interior absorption; vacuum grids preserve
  the ordinary tracer's exact results. With nonzero media, final-vertex NEE keeps
  its complete weight because no competing BSDF continuation is evaluated,
  including a single rectangular light. `AbsorbingFilm` binds the medium snapshot
  and all sample settings, caps progressive samples and commits whole ranges
  transactionally. Callers retain the same scene between appends. This surface
  adds no dynamic-field interpolation, thermal-state debit or measured-property
  authority. G1/G3/G4/G5 tests compare
  slabs, shadow paths and hot-cell emission to independent integrals, exercise
  boundary/refinement/unit changes, and verify replay and refusal rollback.
- **Source-bounded UV/IR transport (`tracer,volumes`):**
  `MassAbsorptionSpectrum::source_table` and
  `MassScatteringSpectrum::source_table` retain 2–16 positive vacuum-nm knots,
  their exact inclusive support and a nonzero caller-supplied source identity.
  That identity binds source/condition data, not measurement authority. The
  retained interpolation is linear; unresolved narrow features and spectral
  quadrature remain caller responsibilities. Existing `table`/`constant`
  constructors and camera sampling keep their 380–780 nm domain. Every query,
  including a zero-power beam lane, refuses outside its actual source support.
  Scattering requires matching coefficient domains and complete phase support;
  `PhaseFunction::source_spectral_table` retains normalized source-domain
  angular laws without extrapolation. Beam deposition preserves source bounds
  in its report, and `OpticalHeating` admits UV/IR absorbed watts into the same
  native enthalpy cells with those bounds intact. The shared SI Planck law and
  emission-direction proposal cover the admitted medium domain, including cold
  far-IR sources; no second spectral probability is introduced. No material
  dataset, unprovided infrared continuation, or bolometric heat-loss closure is
  inferred. Focused G0/G1 tests cover source/domain refusals, UV/IR Beer–Lambert
  heat, separate scattering, cold far-IR emission, and native thermal energy.
- **Physical participating cells (`tracer,volumes`):**
  `volumes::transport::scattering::ScatteringGrid` extends the same accepted
  borrowed density/temperature snapshot with independent mass-scattering
  coefficients and a reciprocal isotropic, Henyey--Greenstein, unpolarized
  Rayleigh or tabulated phase law. Spectral extinction is `rho * (kappa_a + kappa_s)`.
  Exact cell crossings invert exponential free flight without fixed marching
  steps, density copies, global majorants or null-collision loops. A uniformly
  selected packet lane proposes distance; complete spectral-mixture densities
  weight every collision and escape, including zero-extinction lanes. LTE
  emission uses absorption as its source coefficient and total extinction for
  attenuation. Its complete segment integral is added before stochastic path
  weighting, so a sampled collision cannot erase an unsampled source term.
  `render_scattering` and its progressive range entry integrate surface and
  multiple-scattering vertices in the same path loop. Rectangle, environment
  and thermal-mesh NEE evaluate the phase law and exact shadow transmittance;
  light-hit and direct-sample MIS include the same packet escape probability.
  Final vertices retain full direct weight. `ScatteringFilm` binds field,
  coefficients, phase and settings, and publishes each sample range atomically.
  Zero-scattering snapshots retain deterministic absorption; vacuum retains
  the original surface path. This is one finite, index-one grid with global
  spectral coefficients/phase, not overlapping-medium composition, polarization,
  or thermal-energy mutation. With nonzero ambient grids,
  every recursive Fresnel/spectral branch retains that same field after glass
  exit. Specialized surface-only refractive NEE and its reverse competitor are
  disabled together; complete interface random walks supply those paths.
  Ordinary rough-interface NEE uses the correct ambient or glass side, and
  rough Fresnel events have an independent random dimension from free flight.
  Thermal sources embedded in glass remain explicitly unsupported. Scattering
  semantics v4 adds ambient-side rough-interface source sampling while
  preserving cold/vacuum streams.
  Focused regressions cover independent single-scatter source integrals,
  thermal spectra, surface shadow extinction, spectral furnace conservation,
  phase normalization/moments, exact progressive replay and bounded refusals.
  Both physical-media renderers also expose fresh and progressive
  `*_with_execution` entries. They share the original tile pool, render-seeded
  pixel streams, per-pixel accumulation order, operation memory leases and
  execution reports. Workers borrow the same accepted fields. A fresh render
  allocates one film payload plus bounded tile scratch; appends stage privately
  and retain the original field/settings binding. Serial, tile-layout,
  worker-count and split-range comparisons assert exact binary64 pixels.
  Memory, traversal, cancellation and nonfinite-radiance refusals publish no
  partial film. This adds parallel static-field rendering, not adaptive media,
  dynamic-field shutter interpolation or a separate executor.
- **Participating dielectric interiors (`tracer,volumes`):**
  `SceneMedia` binds borrowed absorbing/scattering grids to sorted dielectric
  primitive IDs, together with an optional ambient grid. The top of the existing
  boundary stack selects one physical owner; nested entry displaces the outer
  field and exit restores it. A bound field replaces the host glass's complete
  homogeneous absorption, even when the field is vacuum. The host's spectral
  phase index and smooth/rough Fresnel law remain in force. Grids use world
  metres and vacuum wavelengths; their finite support is clipped at each next
  surface, with vacuum outside the grid. No density/temperature buffer is copied.
  Exact DDA attenuation, free flight, phase scattering, source integration,
  rectangle/environment lighting and both MIS directions use the selected
  region. Every recursive spectral/Fresnel branch retains the binding.
  LTE source radiance inside a host is `n(lambda)^2 * B_lambda(T)`, consistently
  paired with the existing radiance-mode interface factors. Rough boundaries
  sample emitting incident and transmitted regions on disjoint hemispheres;
  each continuation competes only with its outgoing source proposal. Opaque
  receivers evaluate conductor Fresnel response against the actual host index.
  Specialized surface-only refractive NEE is disabled with its reverse
  competitor for these jobs; complete boundary random walks retain the paths.
  Admission bounds binding/primitive/source-triangle metadata, hashes numeric
  scene and field identities, and prepares lighting once. Region lookup is
  allocation-free binary search. `SceneMediaFilm` binds the scene, all supplied
  fields and settings, and fresh/progressive serial/tile entries reuse the
  existing bounded executor and transactional publication. Tests cover replaced
  absorption, explicit vacuum, heterogeneous and nested paths, refractive LTE,
  phase transport, direct-source laws, execution replay and refusal rollback.
  This is static geometry with held world-space fields and spatially constant
  host index. Camera launches use the dielectric admission above. It does not
  certify boundary closure, admit arbitrary overlaps, implement gradient-index rays, or debit
  thermal energy. Surface thermal-emission laws inside glass remain refused,
  including phase-to-surface connections; the implemented interior thermal
  source is the participating field. The `n^2` source convention follows
  [Siegel and Spuckler (1992)](https://ntrs.nasa.gov/citations/20150021009).
- **Accepted medium history (`tracer,volumes`):**
  `MediumTimeline` borrows contiguous physical-state intervals and their
  already admitted absorbing/scattering grids without copying fields or
  allocating during admission or lookup. Content identity binds intervals,
  model kinds and numerical states. Binary search selects the later state at
  an internal discontinuity and the final state at the final endpoint.
  Complete shutter coverage is required; gaps and overlaps are named refusals.
  The exposure model holds each accepted state over its declared interval,
  including explicit grid/model changes. It does not interpolate temperatures,
  average optical depth before rendering, or extrapolate missing simulation
  state. One selected state remains fixed for the complete optical path.
  `MediumExposure` admits that history together with a shutter and optional
  thin-lens/keyframed camera. `render_medium_exposure` selects the snapshot at
  the same absolute shutter time used by the camera and animated geometry;
  camera, shadow, scattering and dielectric continuation rays retain it.
  This integrates radiance over held simulation states rather than rendering
  averaged density or temperature. Vacuum snapshots retain ordinary motion and
  cinematic rendering, and zero-scattering states retain deterministic
  absorption. `MediumExposureFilm` binds the exact history, shutter/sample
  policy, settings and numeric cinematic-camera content, including changes
  made under the same shot identifier. Scene geometry/materials/lights and the
  legacy scene camera remain the caller's fixed owners between appends.
  Fresh and progressive execution entries reuse the existing bounded tile
  pool and transactional publication. Focused tests cover duration-weighted
  extinction/emission, per-sample mixed-state shadow transport, motion/camera
  coherence, exact split-range/tile replay and refusal rollback. This remains
  quasi-static exposure integration: no finite-light-travel-time effects,
  interpolation across missing simulation states, or thermal-energy debit.
- **Accepted scene history (`tracer`):**
  `SceneExposure` borrows contiguous complete `SceneFrame` states and prepares
  each selected state's stationary and moving light proposals once. Explicit
  aggregate limits bound frame, primitive and source-triangle metadata before
  preparation. Each animated trajectory is admitted over its own frame's
  intersection with the shutter. Binary search selects the later accepted
  state at an exact discontinuity and includes the final history endpoint.
  Geometry, topology, thermal temperature/emissivity, optical materials,
  environment and light proposals remain coherent for the entire sampled path,
  including shadow rays and recursive spectral/Fresnel continuations.
  Optional cinematic cameras and physical medium histories use that same
  absolute shutter time, even when their state boundaries differ. Integration
  holds accepted states; it never interpolates missing topology or averages
  temperature before evaluating emission. Scenes and physical fields remain
  borrowed, and lookup allocates no memory or light hierarchy per path.
  `SceneExposureFilm` binds numeric scene/camera/medium history, shutter policy,
  settings and sample prefix. Chart allocations retain the existing explicitly
  process-local continuity guard, not a portable authority-bearing identity.
  Fresh, progressive and bounded tile entries reuse the same path kernel and
  atomic publication. Focused tests exercise independent Planck/reflection
  integration, per-sample frozen-state lighting oracles, topology changes,
  medium/moving-source clock coherence, exact execution replay and rollback.
  With `volumes`, `with_scene_media` attaches one already-admitted `SceneMedia`
  owner to every supplied frame. Primitive-index bindings must borrow that
  frame's exact scene; mismatched geometry and competing ambient-timeline
  ownership are refused. The same shutter lookup selects geometry, ambient and
  interior density/temperature, coefficients, phase functions and refractive
  owners. Accepted frames may change topology or primitive numbering without
  carrying the prior frame's field bindings into the next state. The exposure
  borrows each owner's prepared lighting and releases its duplicate hierarchy
  metadata. Media continuity is hashed once at attachment; lookup and progressive
  append do not rehash cell fields or rebuild proposals. Tests compare every
  sampled state with independent frozen-scene paths, exercise exact internal
  and final endpoints, and verify cinematic-camera, tile and progressive replay
  plus changed-field/cancellation rollback. These per-frame owners hold static
  world-space fields; continuously moving interior grids remain unsupported.
  This is quasi-static uniform exposure rendering; adaptive/AOV history and
  finite light travel time are outside this entry point's contract.
- **Physical volume-source lighting (`tracer,volumes`):**
  Absorbing and scattering snapshots cache an emission-direction distribution
  over their original temperature/density buffers. A source-band Planck proxy
  guides cell selection, with an explicit uniform component retaining cold
  spectral-tail support. The bounded CDF stores one f64 per cell plus one;
  fields with no representable source in that band allocate no proposal. Sampling
  selects a source point only to propose a direction. Exact cell crossings
  marginalize radius into the complete solid-angle PDF, including cells behind
  occluders, while source radiance is integrated only up to the first surface.
  This removes the point estimator's inverse-square near-field singularity
  without a distance clamp. Opaque and phase vertices sample their active
  region's source alongside existing lights; rough dielectric interfaces sample
  each emitting side's explicit region. Rough reflection and transmission use
  the actual wavelength-dependent Fresnel, visible-normal density and radiance
  eta factors. Dispersive suffixes
  retain each lane's own competing density; an interior-facing connection never
  integrates the displaced ambient grid. Complementary balance/power MIS weights
  the source integral before the stochastic free-flight outcome; surface-light survival
  PDFs retain their separate meaning. Terminal vertices receive full direct
  weight. Camera, smooth/null-index dielectric and BSDF-only paths retain full
  source-integral support without a continuous competitor.
  `ScatteringGrid::transport` shares its exact LTE cell integral with free
  flights: absorption creates radiation, total extinction attenuates it.
  A disjoint counter stream preserves existing flight/phase/light identities.
  Physical exposures select the matching source distribution with each accepted
  snapshot, and serial/tiled/progressive entries all use the same estimator.
  Independent cube densities, hot-box irradiance, rough-interface Fresnel/eta
  limits and single-scattering integrals cover normalization, visibility, source
  energy and replay. This does not add refractive multi-interface source
  connectors or participating fields inside a glass body.

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
- Terminal rough-dielectric direct lighting uses unit MIS weight when no BSDF
  continuation is evaluated. Reflection and transmission agree bitwise with
  NEE-only for the same samples, including the one-rectangle lighting path.
  Opaque terminal receivers in dielectric scenes follow the same rule after
  glass exit; the frozen opaque-only estimator remains separately versioned.
- The first two smooth-dielectric boundaries encountered along every path
  suffix are evaluated by exact forced Fresnel splitting: reflection and
  transmission suffixes are both traced with their physical spectral weights,
  rather than selecting one event and dividing by its probability. Every child
  consumes one unit of the private path-local budget before it can continue.
  This is an unbiased sum of the same transport integral, removes roulette at
  the primary interface and its immediate internal reflection/exit decision,
  and still prevents an unbounded specular tree. Third and later smooth
  boundaries retain the stochastic Fresnel estimator.
- The smooth parallel-slab direct-light specialization is an unbiased raw
  estimator, not a straight-through visibility shortcut. For a finite sampled
  rectangle point it solves the two-interface Snell connection, differentiates
  the slab's direction-dependent lateral displacement to obtain
  `dA_light/dOmega_source`, and converts the light's original area density into
  the connected source solid-angle density. For an environment, equal ambient
  media and parallel interfaces make the external angular map the identity.
  Both paths apply the two exact Fresnel BTDF radiance factors and physical
  Beer-Lambert length. MIS compares the NEE density with the opaque BSDF density
  times the event probabilities that the competing BSDF strategy actually
  sampled: a deterministically split first entry contributes probability one,
  while a roulette-sampled entry and every later sampled exit contribute their
  Fresnel probabilities. Physical radiance throughput retains every Fresnel
  factor independently of that sampling density, and the reciprocal weights
  partition unity. Reverse MIS replays the exact forward straight shadow-ray
  spawn from the retained opaque-source geometric normal. Its competing NEE
  density is zero unless that ray first meets the same admitted local slab face
  pair; an unsupported first hit on that same slab remains a typed refusal.
  Failure to classify an ordinary sampled dielectric path as this optional slab
  proposal does not forbid BSDF transport or invent an NEE competitor: it keeps
  normal delta-emitter MIS semantics and therefore unit weight.
- The finite smooth-dielectric manifold specialization extends that estimator
  to an isolated convex, outward-oriented instanced triangle solid with two
  nonparallel planar transmission faces. For each sampled finite emitter point
  it solves the four-variable stationary Fermat path with a damped Newton method
  and analytic Hessian. For an environment direction, exact inverse Snell
  transport fixes both segment directions and their support-plane intersections;
  the induced source density uses the reciprocal two-interface etendue Jacobian,
  with the ambient-to-glass and glass-to-ambient eta-squared factors cancelling.
  Both forms validate the finite face hits and final target against actual scene
  traversal and apply exact spectral Snell, Fresnel, Beer-Lambert, and
  radiance-mode eta transport. The finite light-area to source-solid-angle
  density is the implicit derivative of the same Fermat system, not a
  finite-difference approximation. Forward NEE deterministically
  prefers the analytic slab specialization when its wavelength-specific finite
  connection exists, then falls back to the manifold specialization if the
  refracted path leaves that slab face through a bevel. Reverse MIS replays that
  per-wavelength selection and the complete manifold density, so the BSDF and
  NEE weights partition unity. Multiple distinct admitted
  stationary transmission paths are summed; coplanar triangulation duplicates
  are direction-deduplicated. Non-convergence, a non-convex/inward-wound mesh,
  rough or nested media, and charts simply have no optional manifold proposal
  and retain ordinary BSDF support.
- For finite rectangle endpoints, that same isolated-convex-planar admission
  additionally covers exactly one internal specular reflection (`T-R-T`). The
  six-variable stationary Fermat solve, its analytic Hessian, and implicit
  target-area Jacobian define the forward density; actual entry, reflection,
  exit, and light hits must match the named finite support planes. Coplanar mesh
  triangles share one deterministic lowest-index interface identity in forward
  enumeration, the sampled random walk, and reverse replay. Throughput applies
  exact spectral `T*R*T`, radiance-mode eta transport, and Beer-Lambert length.
  Before Newton iteration, a necessary angular-support test rejects only strict
  impossibilities: air-to-glass Snell transmission bounds the incident internal
  direction to `sin(theta) = 1/eta`, specular reflection maps that cone exactly,
  and glass-to-air transmission requires overlap with the equal exit cone. The
  implementation compares `cos(2*theta) = 1 - 2/eta^2` algebraically and retains
  both equality boundaries under a roundoff guard. This cull changes neither
  finite-path support nor density; finite-face traversal remains authoritative.
  Reverse MIS replays the same reflection support plane and stationary solution,
  so its event-probability density partitions with the forward proposal. This
  is deliberately not arbitrary `T-R^k-T`, an environment reflection connector,
  or a general manifold-next-event solver.
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
while still inside a declared closed medium. A transparent blocker outside the
admitted analytic-slab and finite-convex-planar-manifold NEE classes returns
`UnsupportedSlabNee`; it is never silently reclassified as opaque visibility.
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
shared three-dimensional Sobol direction table for pixel-Owen jobs plus a
six-dimensional reusable path-block table for full-path-Owen jobs, its worst-case concurrent
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
iid, ledgered on bead 872c), or pixel Sobol' plus independently
nested-uniform-scrambled six-dimensional light/BSDF/dielectric-event blocks
per bounce. The latter is a randomized blocked-QMC estimator, not one
monolithic high-dimensional Sobol' net; its quality ordering is measured per
scene at equal fixed SPP and is not assumed. Streams remain bitwise invariant to any
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
  probability in their solid-angle PDFs and use one shared `1e-9` rectangle
  grazing-support boundary, including the formerly mismatched `1e-12..1e-9`
  band. This support repair is lighting tracer bit-semantics v2. Lighting
  bit-semantics v3 replaces the within-rectangle uniform-area draw with the
  area-preserving spherical-rectangle map of Urena et al. (EGSR 2013), whose
  conditional density is exactly `1 / solid_angle`. Projections below `1e-3`
  steradians retain the exact uniform-area proposal because the spherical
  inverse map is ill-conditioned there; forward sampling and reverse MIS use
  that same branch and area-to-solid-angle Jacobian. At a
  lighting-v1 path's final permitted bounce, NEE has weight one because no
  competing BSDF continuation is evaluated.
Lighting bit-semantics v4 adds physical thermal mesh proposals and their matched
hit weights. Scenes without these proposals retain their previous numerical
samples. Thermal spectra never pass through an RGB approximation.
Lighting bit-semantics v5 adds a deterministic thermal light hierarchy above
  eight positive mesh faces. Cancellable Morton ordering and balanced spatial
  partitions produce bounded-depth sampling; node luminance/area/distance
  heuristics choose proposals only, never approximate emitted spectra. Matched
  forward and reverse densities use the same branch probabilities. Bounds-pruned
  endpoint queries sum every overlapping same-primitive proposal; shared edges,
  coincident proposals and null coplanar subtrees retain their explicit support.
  The small-source sampler and all nonthermal image streams remain unchanged.
  Sampling is logarithmic in thermal face count; reverse queries are spatially
  pruned, with worst-case linear work for fully overlapping source bounds.
Lighting bit-semantics v6 adds direct sampling of rigidly moving thermal mesh
  sources in the existing motion, cinematic, adaptive and physical-medium
  exposure paths. Each body's local-space hierarchy is prepared once; a
  stable-object-ordered area/luminance CDF selects a body with logarithmic
  lookup and a uniform support component. Only the selected rigid placement is
  evaluated at the path's absolute shutter time, reusing the geometry pose
  cache when available. Source points, visibility and reverse hit PDFs share
  that time. Proper-rigid motion preserves area and solid angle; neither mesh
  positions nor light hierarchies are rebuilt per ray. Moving sources have a
  separate light technique and counter domain, so static-only light/BSDF/lens
  streams are unchanged. Opaque and rough-interface receivers, and physical
  medium phase vertices, use matching direct/hit MIS with complete final-vertex
  weight. Static and moving source sets are disjoint. Spectral radiance remains
  the accepted graybody Planck law; transforms and observation do not update
  temperature or debit thermal energy. Focused tests include independent
  area/time integrals, mixed-source energy, moving occlusion, physical extinction,
  matched source densities and exact tiled/progressive replay.
Camera-path transport exposes both balance-heuristic and squared-power-
heuristic MIS. Both use the same reciprocal proposal densities and partition
unity; the power variant is an unbiased variance policy, not radiance clipping
or a change to the underlying light-transport integral. Its quality ordering
must be measured per scene at equal SPP.
The exact legacy one-rectangle/no-environment lighting branch retains its
original light-selection draws and estimator. Tracer bit-semantics v2
deliberately changes reflective GGX sample/PDF bits from NDF sampling to
isotropic visible-normal sampling while preserving the same BSDF integral in
expectation.
Current scalar scene-integrator no-claims: no gradient-index propagation or arbitrary overlapping media,
no Russian roulette, reflective GGX
and conductor sampling is isotropic-only, conductor and rough-dielectric GGX
are single-scattering without multiple-scattering compensation, and standalone
RGB emitters do not reflect. Thermal emitters retain their reciprocal diffuse
reflection. Dielectric tracer bit-semantics v14 retains v13's isotropic
view-conditioned visible-normal rough-dielectric proposal, exact environment
endpoint, finite-rectangle one-internal-reflection manifold connector,
coplanar-interface canonicalization, and support-preserving `T-R-T` cone cull.
It raises the private path-local smooth-dielectric forced-split budget from one
to two: the first and immediate second smooth boundaries on every suffix use
the exact Fresnel reflection/transmission sum, every child consumes one budget
unit, and the third and later boundaries remain sampled. Version 13
canonicalizes every coplanar triangulated interface to its lowest triangle ID
for both forward enumeration and reverse BSDF replay; this removes duplicate
Fermat solves and prevents a quad-diagonal hit from losing its competing NEE
density. It also rejects strictly disjoint `T-R-T` Snell/specular direction
cones before Newton iteration, using the same support-preserving test in forward
enumeration and reverse replay. Physical Fresnel transport remains distinct from
the actual competing BSDF event density: either deterministically split
interface contributes
physical `T` to throughput but probability one to MIS, whereas a stochastically
selected entry contributes `T` to both. At the
first two smooth dielectric boundaries, reflection is retained as a spectral
packet and transmission fans out into single-wavelength continuations only when
its Snell directions require it; their deterministic Fresnel-weighted sum
retains support when one wavelength undergoes TIR and another does not. After
the forced budget is exhausted, a later first wavelength-dependent dielectric
boundary can still fan an unsplit packet into four correlated stochastic
continuations. A child wavelength lane cannot fan out spectrally again, fixed
lane/branch-order accumulation is deterministic, and the combined
forced-Fresnel plus spectral-dispersion path-work bound is at most sixteen times
the declared maximum depth. Wavelength-independent transmission retains the
packet.
Current-vertex NEE uses each lane's native competing BSDF density.
Absorption uses unshifted physical segment length. Smooth events have zero
solid-angle query density and receive delta-correct MIS treatment. A strict
path-local stack mutates only after sampled transmission and supports nested,
  closed, consistently outward-oriented media. Opaque-source direct lighting has
  one additional exact specialization when its first shadow hit supplies two
  smooth, homogeneous, locally parallel, piecewise-planar instanced-mesh face
  witnesses satisfying the mesh thin-axis support admission: finite rectangles use the refracted endpoint
connection and area/solid-angle Jacobian, while environment samples use the
unchanged external direction. Geometry is checked again on the connected ray;
rough, nonparallel, non-slab, overlapping, and nested direct-light cases refuse
rather than being treated as opaque or undeviating.

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
grazing, pole-frame, adjacent-IOR, packet retention at nondispersive boundaries,
one-time dispersive fan-out, wavelength-straddled TIR support, exact primary and
immediate-secondary forced-split energy, and two-unit child-budget/traversal-bound
fixtures. The Cornell tracer
battery is deliberately re-frozen when tracer bit
semantics change and remains the opaque-path non-regression gate.
Inline slab-NEE tests additionally pin the equal-IOR reduction, the closed-form
normal-incidence paraxial Jacobian, an independent finite-difference ray
differential, reciprocal MIS weights including both Fresnel event
  probabilities, two-interface plus post-slab visibility, and typed refusal of
  rough, wedge/nonparallel, nested, smooth chart/sphere, and chamfer/side-face
  cases while admitting the central thin-axis faces of a beveled plate. Reverse-MIS
  regressions prove that a finite-edge straight miss and an opaque first blocker
  produce zero competing NEE density and therefore unit BSDF weight. Lighting
  tests pin identical forward/reverse rectangle support within and above the old
  grazing mismatch band.

`conductor` and `tracer` inline tests plus the Euler scene-bridge E2E fixture
(feature `tracer`, and `cinematic-render` downstream) cover complex-Fresnel
known answers and an independent real-form grid, normal/grazing/common-index
relations, complete-table and wavelength refusal, exact-knot interpolation,
packet permutation, source/material fingerprint sensitivity, roughness bounds,
GGX visible-normal sampling/PDF agreement at grazing incidence and both pole
frames, directional-PDF mass versus sampler acceptance, reciprocity,
active-incident-medium dependence, a numerical no-energy-gain furnace, an
independent equal-expectation furnace check, equal-draw grazing variance against
the former NDF sampler, rough-dielectric visible-normal sample/PDF replay,
exact replay, and visibly distinct representative tungsten/stainless renders
under one neutral light. The Cornell golden remains the legacy-material bit
gate.

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
jobs prove that they neither allocate nor charge unused Sobol direction state;
pixel-only and full-path modes charge exactly three and nine direction tables.
Failures name
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
  motion blur are implemented under the contracts below. Rigidly moving thermal
  meshes support direct-light proposals; deforming emission, scaling, moving
  chart-source proposals and time-varying surface emission remain successor
  work. Emissive geometry may move with its instance, but its material emission
  itself is time-invariant.
- v1 includes the scalar-BVH spectral path tracer. Wide-BVH SIMD traversal,
  watertight ray-triangle tests, a LIGHT-BVH for large emitter populations,
  moving participating-grid interiors, ray-stream sorting, and
  progressive tile streaming to HELM remain staged.
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
  Full-path Owen sampling is unbiased only under its fixed-SPP randomized-QMC
  interpretation; data-dependent stopping does not inherit that statement.
  Independent-pilot allocation preserves the fixed-count interpretation only
  because the pilot is discarded and a distinct production seed is required;
  it does not certify that the selected count meets the requested image error.
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
  preset fidelity, arbitrary overlapping
  media, or a topology certificate. GGX transmission is single-scattering and
  does not claim multiple-scattering furnace closure at appreciable roughness.
  Outside the explicitly checked smooth homogeneous local parallel-face-pair
  specialization and the isolated convex planar two-transmission-interface
  manifold (optionally with one finite-rectangle internal reflection) on an
  instanced mesh, shadow rays stop at the first intervening surface; they never
  travel undeviated through glass. Curved, rough, layered, overlapping, and
  nested manifold connections are refused, not approximated.
  The manifold connector covers transmission through one finite convex planar
  solid with zero or one internal reflection; it is not general bidirectional
  path tracing, vertex merging, or a
  solver for arbitrary specular chains. Difficult focused caustics outside that
  class can still converge slowly. No denoising or unbiased firefly-clamping
  claim is made.
  The support-face/thin-axis extent check, two instanced triangle witnesses,
  and path-local probes do not constitute
  a global closed-solid, planarity, thickness-uniformity, or topology
  certificate; support is deliberately limited to the locally observed planar
  face pair.
  Fixed-metric, representability-aware ray offsets are robust engineering, not
  a certified positional-error enclosure.
- Conductor support is an opaque complex-index substrate with optional explicit
  local coherent transparent or absorbing-film interference and single-scattering GGX.
  The representative presets do not infer an oxide/passivation film or claim
  material identification, alloy/lot or product calibration, a measured finish,
  temperature dependence, polarization transport, interfacet multiple scattering,
  or a conversion from measured `Ra`/`Rq` to GGX alpha. Caller-asserted measured
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

## Incident optical power and thermal absorption

`AbsorbingGrid::deposit_beam_power` and `ScatteringGrid::deposit_beam_power`
retain the original exact cell crossings and return directed transmission,
true absorption and scattering separately, in quadrature-weighted watts per
spectral lane. The optional `optical-heating` feature maps only the absorbed
cell powers through `volumes::heating::OpticalHeating` into the existing
`fs-thermal::enthalpy::DepositionObserver` source hook. Stable, explicit
cell-index bindings permit aggregation into a resolved thermal control volume;
missing, duplicate, inactive or nonconservative assignments refuse.

The native `advance_observed` owns caloric storage (including latent heat),
conduction, sources, event alignment and transactional publication. Absorbed
power is counted once in its external-heat audit. The optical observer retains
incident, transmitted, absorbed and scattered joules over those same accepted
substep durations; cancellation or domain refusal discards both histories.
Prepared heat rates require no optical traversal during thermal substeps.
The caller must start from active receiver cells and recompute deposition at
accepted endpoints when beam powers, optics or geometry change. This fixed
declared-band incident-beam observer does not infer net radiative cooling,
unsupported infrared properties, multiple-scattered reabsorption or
material-temperature feedback. Scattered power remains radiation and is never
converted to heat.

### Temperature-dependent optical heating

`volumes::heating::feedback::TemperatureOpticalHeating` retains source-bound
mass-absorption tables versus kelvin and explicit original-grid/thermal-cell
bindings. It caches the existing DDA's physical cell lengths/densities and
samples each original spectral law at the beam wavelengths once. Each native
thermal substep then obtains actual temperatures from the registered caloric
laws and recomputes serial Beer attenuation, including changing upstream
shielding. Every non-vacuum path cell needs a binding even when initially
transparent. Unsupported temperature/wavelength states refuse without
extrapolation; spectral rows and their source identities remain explicit.

The shared thermal owner continues to integrate all conduction, latent storage,
authored sources and scheduled events. Optical energies integrate the actual
changing watts with compensated sums over those exact accepted durations.
`advance_refined` compares complete native coarse/fine/finer trajectories and
returns only the finest original state after decreasing discrepancies satisfy
declared enthalpy, temperature and optical-energy tolerances. Two comparisons
at the documented binary64 roundoff floor may also pass those same tolerances
without a measured-order claim, so exact constant-source solutions do not
refuse merely because their final rounding is nonmonotone. Relative scales
use changes, and optical comparisons retain compensated energy differences;
large formation-energy datums or historic fluence cannot hide discrepancies.
All trials share a whole-call thermal work ceiling. Domain, cancellation,
conservation, budget or convergence refusal publishes no state or partial heat.

This is first-order explicit incident heating with measured time discrepancies,
not a certified global error bound or an implicit nonlinear solve. Geometry,
density, beam powers and spectral quadrature remain fixed; data validity,
spatial/spectral accuracy, scattering, thermal emission and thermochemical
changes require their respective physical owners. Tests compare an independent
closed-form nonlinear Beer-heating trajectory and its refinement order, dynamic
two-cell shielding, native latent/event/split behavior, infrared source support,
initially transparent cells and transactional refusal.

### Multiple-scattering absorbed power

`ScatteringGrid::deposit_multiple_scattering_power` follows the same native
exact-cell extinction free flights and wavelength-resolved phase laws through
the full open, index-one grid. Independent monochromatic source histories use
implicit capture: each collision deposits only its true absorption fraction
and continues its scattering fraction. Absorption-only lanes reuse the exact
directed-beam path. Source powers are already quadrature-weighted watts; LTE
emission is independent and is neither credited nor debited by this operation.

Every history retains ballistic escape, diffuse escape, cell-resolved absorbed
power and unresolved radiation. Collision caps and below-representable packet
splits retain the surviving packet in the unresolved account, including a
separate count of precision-limited histories. No roulette or hidden power
cutoff is used. Segment-work and distinct absorbing-cell limits are explicit;
hard work/storage exhaustion, cancellation or invalid arithmetic returns no
partial report. Sparse cell storage is bounded by the declared output count;
source/lane/history counter streams and ordered reductions are deterministic
independently of hash-table iteration and worker identity.

Compensated survival updates retain cumulative absorption at near-unit albedo;
absorption-dominant splits preserve their small tails through direct coefficient
ratios. Scaled Welford moments avoid squaring raw SI powers. Repeat visits to a cell
are summed within each history before calculating its variance, and total
absorption uses complete history sums rather than incorrectly summing cell
variances. Retained standard errors are fixed-sample estimates, not confidence
intervals, optional-stopping guarantees or certified bounds. Unobserved rare
events can have zero empirical error. Mean unresolved power estimates the
omitted radiative tail; its own sampling error is retained.

Under `optical-heating`, `volumes::heating::scattering::ScatteringHeating`
maps only the absorbed estimates to the original enthalpy observer protocol.
Admission requires explicit allowances for unresolved power, total sampling
error, and receiving-cell sampling error. Many-to-one optical/thermal mappings
conservatively sum the contributing optical-cell standard errors, retaining
possible spatial covariance. Native accepted time integrates radiation and
heat together; fixed-field sampling errors scale linearly with duration across
substeps. Authored heat sources remain additive and absorbed power enters the
native energy audit once. Mapping, source admission and advance refusals
publish neither heat nor radiation accounts. This model does not include
reflecting/refracting grid boundaries, surfaces, fluorescence, moving geometry,
temperature feedback or net thermal radiation exchange.

G1 tests compare exact Beer transport and heterogeneous spectral
first-collision expectations. A G2 reference independently samples analog
absorption in a homogeneous cube, with separate ray/box geometry and HG
sampling, for backward, isotropic and forward phase laws. G0/G3/G4/G5 tests
cover repeated-visit covariance, high-albedo truncation, precision-limited
tails, source replay and bounded refusal. Native thermal integration tests
cover actual enthalpy/temperature rise, many-to-one uncertainty admission,
pure scattering without heat, accepted-time error propagation and rollback.

## No-claim boundaries (volumes)

- FrankenVDB tile-maxima majorants: no fvdb crate exists in-workspace;
  [`MajorantGrid`] builds per-block maxima from dense grids, and the
  per-tile-rate DDA traversal (rather than lookup thinning under a
  global bound) is the recorded successor alongside the FVDB wiring.
- Progressive live tiles with ledger artifact pinning (frame-consistent
  snapshots of evolving fields) — staged with the vessel flagship's
  render lane; the smoke tier renders a paused simulation's buffer.
- Refractive free-surface reconstruction from fill fractions and specialized
  refractive volume-source connectors remain successors. The shared path tracer
  implements borrowed absorption/emission/scattering in ambient and explicitly
  bound dielectric regions, with phase/source MIS; interface random walks
  supply paths that cross boundaries.
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
keyframe is bit-identical to the equivalent static instance when no moving thermal
NEE technique is active. Moving thermal sources have equal physical expectation
under the same pose, but their dedicated proposal stream is different from the
static mixture. G0/G3/G5 coverage
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
The geometric camera model makes no claim for lens distortion, chromatic aberration,
rolling shutter, vignetting, autofocus, sensor irradiance, or
photometrically calibrated exposure. `ExposureMetadata` is retained but is not
applied by the current radiance-averaging integrator. `FilmTimeMode` binds the
render-path kind, shutter, stream seed, and cinematic shot ID, so progressive
appends cannot cross a cut or switch between legacy-motion and cinematic rays.
Like the pre-existing progressive tracer, it does not yet content-bind the
scene or full camera content; composition-level artifact identity owns that
broader provenance closure.

### Focused pupil diffraction

`PhysicalCamera::with_in_focus_diffraction(FocusedDiffraction)` opts the
cinematic spectral tracer into scalar, paraxial, incoherent pupil imaging.
It requires a nonzero circular or regular 3–64-blade aperture,
`lambda/(2*inradius) <= 0.01` (the usual radius for a circle), and an explicit
maximum pupil-rim defocus phase between zero and 0.1 radians. A physical
focal-length declaration requires focus beyond that focal length, admitting
only a positive conjugate image distance. First-hit surfaces must have an
opaque diffuse BSDF (including diffuse thermal emission); specular, conductor
and dielectric first hits refuse even when they also emit, because their
virtual-image depth cannot be inferred from the surface position. Paths after
an admitted diffuse first hit retain ordinary radiance transport.
At each primary visible surface the tracer checks
`pi*a^2/lambda_min * abs(1/z - 1/focus) <= maximum_phase`, using axial depth.
An emitting environment uses `1/z = 0`; an empty black background needs no
depth claim. Refusals preserve the film transactionally. Camera animation
retains the policy and render-continuity hashes distinguish it from geometric
imaging. Existing geometric render bits and random streams are unchanged.

For a circular pupil, the in-focus image-plane PSF is
`[pi*D/(lambda*q)]^2/(4*pi) * [2*J1(v)/v]^2`,
`v = pi*D*r/(lambda*q)`, where `q` is the actual image distance. Its integral
is one, not its peak. Since this camera's raster is angular (`r=q*theta`),
the angular Airy scale is `lambda/(pi*D)` at finite conjugates as well.
The existing `fs-math` Bessel implementation supplies the rings. A shared
two-variate Cauchy proposal has the same inverse-cubic radial-density tail,
with no rejection loop, image-space convolution buffer, or ring cutoff.
Each wavelength receives its own PSF/proposal weight before spectral XYZ
integration; no wavelength is treated as a representative achromatic blur.

Regular bladed pupils use the analytic Fourier transform of their actual
polygon edges. With unit-circumradius pupil area `A`, transform `F(q)` and
`q=2*pi*circumradius*theta/lambda`, the normalized density is
`abs(F(q))²/(4*pi²*A)`. Edge integrals use sinc limits, compensated sums and
a second-moment DC limit to avoid removable singularities. This retains
orientation-dependent diffraction spikes and full tails without a pupil
raster, FFT aliasing, an image convolution buffer, or per-ray allocation.
A mixture of blade-aligned product-Cauchy proposals resolves narrow
inverse-square ridge tails as well as generic inverse-fourth-power tails.
Quarter-turn-equivalent proposal frames are shared. The full mixture PDF is
used for every spectral lane; work is bounded by the admitted 64 blades,
and cancellation is polled between wavelengths. Exact random endpoints use
nearest open-interval floating-point representatives. The prior circular
sampler's arithmetic and random domain are retained.

`in_focus_diffraction_density` exposes this same analytic PSF per unit
angular-raster area. Its central density is actual pupil area divided by
wavelength squared, and its integral is one. The existing circumradius
defocus-phase admission remains conservative over the whole polygon pupil.
The edge-integral derivation follows Gallatin, *J. Math. Phys.* 53, 013509
(2012), equations (5) and (21),
<https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=909668>.

The focused mode traces the optical-centre image and does not additionally
sample geometric pupil blur. Its phase budget is a model admission condition,
not a certified image-error bound. It does not solve defocused wavefronts,
coherence between scene points, pupil-dependent occlusion, aberrations,
vector/polarization diffraction, high numerical apertures, or arbitrary
edge diffraction in scene transport. Primary volumes and calibrated
photometric fixture paths refuse this mode; bidirectional rendering retains
its existing pinhole-only admission. Wavelength-unaware camera-ray callers
receive `DiffractionRequiresSpectralPath`, rather than silently losing the
selected optical model.

G2/G3 regressions in `tests/camera.rs` check the first Airy zero, wavelength
and pupil scaling, separate packet ring positions, unit total power (0.2%
quadrature band), the 0.837784869 central-disk energy (1e-7 absolute), phase
admission and animation. `tests/tracer_battery.rs` checks actual image
contributions outside geometrical support, flat-environment conservation
(4% finite Monte Carlo band), exact progressive replay, and atomic defocus
refusal. Pure polygon-kernel tests compare the independent separable square
aperture sinc² solution (including zero and edge-parallel limits), verify
rotation and SI scaling, and integrate all three packet wavelengths with
fixed Hammersley cubature for 3, 4, 6, 9 and 64 blades (0.4% power band).
Camera tests bind pupil area, wavelength scaling and sampled orientation to
that kernel. These are numerical model checks, not measured lens validation.

### Finite-defocus circular-pupil imaging

`PhysicalCamera::with_defocused_diffraction(DefocusedDiffraction)` selects
scalar propagation of a circular pupil at a declared object plane. The
policy takes a positive axial object distance in metres (`Some(z0)`) or
infinity (`None`), plus a residual depth-phase budget in `[0,0.1]` radians.
It replaces the focused policy rather than combining the two. The camera's
focus still sets the sensor conjugate. At each wavelength, the modeled
pupil-rim phase is `phi=pi*a²/lambda*(1/z0-1/focus)` with `1/z0=0` for
infinity; `abs(phi)<=16` is required. Actual first-hit depths must separately
meet `pi*a²/lambda_min*abs(1/z_actual-1/z0)<=maximum_depth_phase_error`.
Thus finite defocus is explicitly propagated, while the residual budget
only admits a small mismatch from the declared source plane.

The real spectral renderer preflights a deliberately bounded scene model.
Forward geometry must be static triangle meshes (including fixed rigid
mesh instances) in one frontoparallel opaque diffuse object plane. Vertex
coplanarity uses coordinate-roundoff tolerance, independently of the optical
residual-phase budget. Finite-plane scenes have black misses and no competing
environment. Other geometry may lie wholly behind the pupil plane, where it
can illuminate the object without clipping a primary pupil-to-object path.
Infinity mode permits an environment with all geometry behind the pupil.
The camera eye and basis stay fixed throughout the owning shot; focus racks
remain supported. Unknown charts, animated geometry, prepared scene histories,
primary media, calibrated photometric fixtures and specular primary surfaces
refuse. The mesh scan occurs once at render admission, and actual sampled
hits retain the depth check. Camera-only consumers must establish equivalent
clear-pupil and pupil-independent-radiance conditions themselves: a first-hit
depth check cannot detect an obstacle that clips only part of a pupil.

For dimensionless image radius `q=2*pi*a*abs(theta)/lambda`, the propagated
amplitude is `A(q,phi)=2*integral_0^1 r*J0(q*r)*exp(i*phi*r²) dr`; the
two-dimensional normalized density is `abs(A)²/(4*pi)`. Its angular density
integrates to one. On axis, `A=exp(i*phi/2)*sinc(phi/2)`, retaining central
nulls and changing ring structure through focus. At zero phase it recovers
the existing Airy function. Opposite defocus signs have the same intensity
for this uniform circular pupil. Each wavelength has its own phase, angular
scale and importance weight; neither an achromatic blur nor a convolution
of focused Airy and geometric DOF is substituted. Rays originate at the
optical centre and the geometric lens draw is not applied again.

The stack-only numerical owner evaluates a 40-term Bauer/Zernike expansion,
whose convergence depends on the admitted phase rather than sampled radius.
The omitted amplitude is at most `1.849e-23` for `abs(phi)<=16`, before
floating-point error. Spherical coefficients use bounded backward recurrence
with a nonvanishing normalization anchor. Cylindrical Bessel values use
rescaled Miller recurrence and Neumann normalization through `q=128`, then
stable upward recurrence from Hankel order-zero and shared `fs-math` order-one
starts. Small arguments use explicit ratio series; angle addition preserves
Hankel phase at large arguments. No pupil raster, FFT, runtime dependency,
per-ray allocation, rejection loop or radial cutoff is introduced. A Cauchy
proposal of scale `max(2,2*abs(phi_reference))` covers the broadened core and
inverse-cubic tails; packet wavelengths share a ray and receive independent
PSF/proposal weights. Cancellation is polled between spectral lanes.
`defocused_diffraction_density` exposes the same physical angular density.
The expansion follows Braat et al., *Progress in Optics* 51 (2008), equations
2.45 and 2.61, <https://nijboerzernike.nl/_PDF/BraatProgOpt2008_20140413.pdf>,
with the coefficient bound from NIST DLMF 10.54.1 and 10.60.7.

Pure kernel tests check the axial sinc law, exact focused limit, sign
symmetry, independent adaptive-integral landmarks across all numerical
branches and into the far tail, and three-wavelength unit power with 65,536
deterministic radial draws (0.4% quadrature band). The last representable
open-unit draw retains finite importance weights. Camera tests cover physical
depth and SI scaling, chromatic axial nulls, lens-sample independence,
shortest-wavelength phase refusal, infinity and focus-rack retention. The
phase cap bounds numerical work; it does not certify the paraxial, high-NA,
wide-angle, lens-aberration or measured-lens error. General multi-depth
occlusion, defocused polygon pupils, virtual images, coherent scene transport
and vector diffraction remain outside this mode.

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

For an aligned cinematic AOV sample, previous-frame transforms and cameras are
evaluated at the same shutter phase one presentation cadence before that
sample's absolute mechanics time. The previous presentation timestamp defines
the cadence; it is not substituted for shutter-sample time. This is
load-bearing for front- and back-loaded exposures because the temporal denoiser
consumes the immediately preceding shutter-integrated image, not an image at
the encoded presentation timestamp. The aligned AOV exports no next-frame
motion and does not evaluate mechanics beyond a terminal trajectory cut.

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

## Accepted material-body geometry (MR41, fixed-topology rung)

The optional `material-geometry` ingress maps `DeformedBody` (including the
original `ThermalBody::body`) to the original `TriMesh` and deterministic BVH.
`TetLoadGeometry3` owns finite/quality/incidence admission and exterior-face
extraction. Original shared vertex IDs and accepted metre coordinates are
preserved; internal faces are omitted and winding follows the solid owner's
outward normals. Ray hits, silhouettes, shading normals and visibility use the
actual accepted shape. No visual displacement, alternate solver or thermal
state is introduced. The caller replaces a scene shape only after success.

The original solid node/tet/face-entry bounds apply before extraction; checked
renderer index bounds and finite cross-product norms greater than 1e-12 square
metres admit the current triangle backend. Cancellation polls extraction and
each copied vertex/face, then before and after the original bounded synchronous
BVH build. Refusal/cancellation returns no partial scene mutation.

The source is a conforming, non-overlapping P1 volume mesh; face-incidence
admission does not establish global injectivity or T-junction freedom. Existing
instances place the body's metre frame. Appearance is explicit through the
existing material bindings; a combined thermoelastic fixture binds the same
accepted temperature/card to gray-diffuse thermal emission. Regional optical
assignment is explicit through the regional ingress below. Remeshing,
melt/free-surface evolution, and live trajectory sampling remain outside this
fixed-topology rung.

`material_geometry::regional::deformed_body_regions` partitions the globally
extracted exterior by the original region-to-element map retained by
`DeformedBody`. Canonical material-region order does not imply solver element
order. Each patch retains its original node IDs and four bitwise accepted cell
coordinates; only exterior triangles enter the original mesh/BVH. Bonded
internal faces stay absent, and storage is linear in the admitted cell count.

With `material-state`, `thermal_body_primitives` matches explicit surface bindings
by region key and requires their exact accepted constitutive-state identities.
Each exterior cell can therefore carry its own conductor optics or gray/diffuse
thermal emission at its actual accepted temperature. Interior cells do not gain
optical-data requirements. Stale, foreign, missing and duplicate assignments
refuse before returning any candidate scene replacement. This observer adds no
thermal debit or state update. Dielectric bindings refuse: an exterior-only
regional partition does not supply the internal optical interfaces and shared
medium-boundary tracking required for transparent multi-material bodies.

`continue_thermal_body_primitives` connects later accepted local solid states to
the same native regional geometry/appearance owner. It reevaluates each exterior
region's original optical or total-emissivity claims once at that cell's complete
accepted query and quantity schema, preserving its optical family, claim pins,
independent surface identity and finish. Competing later claims cannot replace
the originals. Geometry extraction and regional BVH construction occur once;
the shared primitive builder validates the complete candidate. Successor sources
and primitives return together for atomic caller publication and the next
endpoint. Unsupported source domains, stale/foreign pairs and cancellation expose
no partial replacement. This remains a fixed-topology opaque endpoint observer,
with no extra heat debit or inferred surface history. The regional continuation
regression drives real thermoelastic deformation, conductor reflection and
gray/diffuse thermal pixels, including later overlapping claims and late refusal.

## Calibrated point-fixture photometry

The opt-in `photometry` feature composes `tracer` and `volumes` for a headless
SI reference calculation. `PhotometricSource` takes the existing proper-rigid
fixture transform, emitted lumens, linear dimmer, angular distribution and
visible spectrum. Local +Z is the optical axis; callers can use a right-handed
Z-up world without coordinate conversion. Isotropic, uniform-cone, cosine-power,
explicit beam/field/cutoff and tabulated theta/phi laws integrate to unit flux.
Tables interpolate in cos(theta)/azimuth with exact cell solid-angle weights;
this is the declared represented distribution, not an interpolation-error bound.
Passive spectral filters reduce flux instead of renormalizing away attenuation.

`PhotometricScene` borrows the original tracer primitives, native mesh BVHs and
rigid mesh instances. Chart geometry refuses: a proximity-certified hit does
not determine whether a surface is before a finite light-connection endpoint,
and must not silently turn a chart receiver into its own shadow.
`IlluminanceProbe`, `CalculationGrid` and `IlluminanceJob` return the
direct, unscattered component in lux using inverse-square, receiving cosine,
visibility and Beer–Lambert extinction. Point-source near-field singularities,
out-of-range connections, transparent materials, animated geometry and legacy
uncalibrated emission refuse. Existing opaque material-state bindings remain
usable. Receiver endpoint exclusion is an explicit metre-valued tolerance and
is not substituted for the physical transport distance.

`HomogeneousMedium` reuses the volume transport/phase kernels with nonnegative
inverse-metre absorption/scattering and isotropic, HG, Rayleigh or tabulated
scattering. Its public `PhaseFunction` is the native volume owner's type.
The single-scattering integrator analytically integrates camera-path extinction
within midpoint panels and refines under an explicit panel budget. Its returned
quadrature discrepancy is an estimate, never an enclosure. The scene image
entry point partitions a camera segment at fixture cone/cosine cutoffs and
tabulated polar/azimuth boundaries, including the azimuth seam. Since photometry
semantics version 6, it also splits at every interior tabulated phase-law cosine
knot for each active point fixture. Those phase crossings use photon propagation
toward the camera, independent of the fixture's optical frame, and retain the
admitted camera direction's actual norm. Isotropic and full-sphere fixtures
still contribute phase crossings. Collinear sources outside the admitted
segment have constant phase and need no phase cuts; sources on the segment
remain subject to the existing singularity refusal.

Beam intersections use the original rigid fixture frames and represented angular
laws. Opposite
cone nappes and azimuth half-planes are excluded by their geometric signs,
retaining ambiguous zeros; distinct roots are never epsilon-merged. Source-distance
and angular boundaries share one exactly deduplicated partition. Each
interval then receives its own midpoint refinements, so a narrow represented
beam or phase lobe cannot evade all sampling merely by lying between the
original uniform grid's nodes.

Since photometry semantics version 8, scene images also partition near-source
distance scales: closest approach to the finite camera segment and offsets of
`d_min`, `2*d_min`, `4*d_min`, and so on, with the actual ray-direction norm
retained. Collinear sources outside the segment use their distance to its nearest
endpoint. A segment shorter than one natural source-distance scale needs no added
distance cuts. These intervals capture inverse-square variation even when a
fixture's angular law adds no boundaries; they do not validate the near-field
point-source approximation or fixture calibration.

Since photometry semantics version 9, smooth cosine-power beams also contribute
their natural intensity contours: `exp(-1)`, `exp(-2)`, `exp(-4)`, and so on
relative to the on-axis intensity. The source caches the first contour cosine;
successive squaring reaches the actual cutoff or floating-point zero without
discarding a chosen tail fraction. These contours reuse the native cone-crossing
geometry, exact deduplication and shared panel cap. A positive exponent whose
first contour rounds to cosine one is refused as `SourceError::NumericalRange`.
An exponent of zero introduces no smooth-profile contours. This sampling
resolution rule is not an integration enclosure or a measured-fixture claim.

The intervals share `ScatterSettings::max_panels`, with at least one panel per
interval and equal shares of the remaining capacity. Exceeding that capacity
returns `ScatteringPartitionBudgetExceeded`; no required distance or angular
interval is discarded.
Total callback evaluations remain less than four times the total cap. A missing
refinement comparison in any interval leaves the aggregate discrepancy unknown
and `tolerance_met` false. Otherwise the discrepancy is the sum of transmitted
interval discrepancies, checked against the requested global tolerance.
Photometry semantics version 3 records the angular partitioning.

Since photometry semantics version 7, optically thick intervals use geometric
optical-depth panel boundaries. For extinction `sigma`, length `L` and `N`
panels, the boundaries are `expm1(i*ln(1+sigma*L)/N)/sigma`, with endpoints
exactly zero and `L`. Metric midpoints retain analytic Beer panel weights.
This resolves foreground source variation that uniform-distance midpoint
grids could miss completely, while retaining the full ray interval rather
than truncating it at a chosen number of optical depths. Intervals with
`sigma*L <= 1` retain their existing uniform arithmetic.

Refinement cannot report `tolerance_met` until the first panel spans at most
one optical depth. A cap too small to satisfy this prerequisite returns the
computed estimate and any actual refinement discrepancy, with the tolerance
flag false. The partitioned wrapper preserves this prerequisite for every
interval. Panel and callback caps, cancellation, deterministic ordering and
allocation limits remain unchanged. This resolution check is not a numerical
enclosure and does not prove arbitrary incident-field variation resolved.

This partition is not a certified intersection or integration enclosure.
Unresolved shadow boundaries, untabulated analytic phase peaks, smooth
within-cell variation and floating-point geometry can still defeat a
refinement estimate. The public
callback-only integrator has no fixture-profile metadata and retains its stated
aliasing limitation. Omitted multiple scattering has no finite bound. Gray
coefficients are an explicit wavelength-independent approximation.

Source spectra use the published CIE 1931 2-degree observer (CIE dataset
`10.25039/CIE.DS.xvudnb9b`) linearly interpolated at 80 five-nanometre midpoints
over 380--780 nm, with the 683 lm/W photopic conversion. Since photometry semantics
version 2 and path semantics version 9, this observer is used consistently for
source normalization, filtering, direct response, all path XYZ and probe area
NEE. The legacy display renderer retains its analytic fit. The original CIE
CSV, metadata, published checksum and separately licensed CC BY-SA 4.0 adapted
table are retained under `data/cie`; `photometry::observer::OBSERVER_ID` names
the data and interpolation convention. This removes the observer-fit mismatch
at represented wavelengths, not spectral-bin/truncation or metrological error.

Input calibration, geometry/material/medium uncertainty and general omitted
transport are not certified. The restricted diffuse depth bounds below
concerns only that model's omitted reflection orders. Direct and path estimates
carry `NumericalKind::NoClaim`,
not a fictitious zero-width certificate; unknown total error is `None`.
Sampling standard error and haze quadrature discrepancy are separate estimates,
neither a confidence interval nor a total error bound. A zero estimate does not
exclude unobserved paths, depth bias or unresolved beams. There is no claim of
CIE 171 conformance or measured-fixture agreement. G0/G1/G3/G4/G5 tests exercise flux, spectral filtering,
inverse square/cosine laws, real mesh visibility, homogeneous reference integrals,
frame invariance, cancellation and deterministic progressive partitions.

Jobs borrow immutable scene/receiver state, reserve their bounded result storage
fallibly and publish only whole completed receivers. Every source, intersection
and scattering sample polls `Cx`. A fresh uncancelled context can resume the
retained prefix; partitioning changes no sums. Determinism is same-ISA/toolchain
and fixed caller source order. No GPU, window, file parser or random stream is
required. Cross-ISA bit equality and large-rig throughput are separate claims.

`LuminanceJob` uses the existing pinhole `PhysicalCamera` and the native
diffuse (including tabulated spectral reflectance), GGX and spectral conductor
BRDFs. Its row-major pixels report
cd/m² before exposure or tone mapping, with direct surface and single-scattered
haze components separately available. Primary opaque geometry terminates both
components at the same hit; fixture-to-surface/haze and surface/haze-to-camera
extinction are each applied once. Camera segments that enter an active point
source's exclusion radius in scattering media refuse before quadrature, even
if no midpoint would sample the singularity. Finite apertures refuse.

This image model omits interreflection, scattering onto surfaces, multiple
scattering and finite emitter appearance. Pixel centres are point samples, not
pixel-area integrals. Neither a haze refinement estimate nor a completed image
certifies those omitted terms. G1 image tests check rho*E/pi, spectral reflection,
both extinction legs and HG propagation direction; G3/G4/G5 tests check shadows,
primary truncation, singularity refusal and progressive cancellation/replay.
The `photometric_stage` example accepts ordinary Rust fixture/scene data and
prints an illuminance summary and luminance diagnostics, or a lux CSV with
`--csv`; its fixture profiles and haze coefficients are declared examples.
Both CSV modes name unknown total bounds. Stderr records the observer and
semantics versions, units, numerical claim, sampling interpretation/seed/depth,
and separate unknown calibration, discretization and physical-model errors.

### Independent physical reference tests

The following are independently specified mathematical regression references,
with fixed test tolerances; passing them does not certify arbitrary scenes.
Existing inverse-square/incidence, disk refinement, spectral extinction and
progressive replay tests remain in the photometry batteries.

- `observer` unit tests independently parse the original CIE CSV and check all
  80 midpoint XYZ triples. `photometric_path::spectral_reference` compares
  selective filters, colored reflection and finite-source measurements with
  the published data. Two equal-power blue/red bins expose an observer error
  that an unfiltered lumen-normalization round trip would cancel.
- `photometric_path::area` evaluates inclined, off-axis receivers against a
  spherical-polygon boundary integral of projected solid angle. A closed,
  uniformly emitting diffuse cube tests the finite-depth geometric reflection
  series, camera/probe units, terminal MIS weighting and absorbed-power balance.
  Its analytic vacuum tail is a property of that benchmark only.
- `photometric_path::point_energy` instead lights passive cube walls with two
  off-center point fixtures totaling 1000 lm. Integrating lux over all six
  walls checks source injection and repeated diffuse reflection against
  `Phi * sum(rho^k, k=0..depth)`. The absorbed flux plus the known finite-depth
  residual must recover the emitted flux. Black walls, two spatial quadrature
  orders and three reflection depths distinguish direct, sampled and omitted
  contributions. This is an area-integrated balance; local wall illuminance
  is nonuniform. Its gray closed-cavity tail and fixed regression tolerances
  are not general error bounds.
- `photometric_image::scattering_reference` compares three non-collinear
  fixtures with an independent angle-space airlight integral for isotropic,
  forward/backward HG and Rayleigh phases. It checks passive extinction and SI
  scale covariance. Independent narrow-cone, annular polar-profile and azimuth
  seam integrals check angular partitioning and rigid-frame equivalence. A raw
  callback counterexample retains the distinction between zero refinement
  discrepancy and a genuine enclosure.
- A conservative-scattering diffuse-cube furnace tests all-order radiance and
  energy balance at two path depths. Its acceptance tolerance does not provide
  a general bound on scattering-order truncation.
- `photometric_path::equilibrium` couples colored diffuse walls, calibrated
  wall emission and a spectrally absorbing, emitting, scattering gas at 3000 K.
  Independent SI Planck constants and original CIE CSV rows define the
  all-order equilibrium radiance and the probe reference `E = pi * L`.
  Uniform and layered densities, a nonuniform tabulated phase, camera XYZ and
  two differently oriented probes exercise absorption/emission balance and
  conservative redistribution together. The benchmark separates a derived
  model-tail allowance from its fixed IID sampling tolerance; it does not
  assert continuous-spectrum accuracy or measured-material agreement.
- `photometric_path::depth_remainder` checks the restricted public model bound
  against an emitting cavity's independent geometric series, including the
  native depth convention and lux conversion. It checks global colored
  reflectance maxima, progress-independent metadata, exact black tails and
  refusal of unsupported models while retaining numerical `NoClaim`.
- `photometric_path::media_remainder` checks the joint surface/volume bound
  against an independent colored-albedo maximum at an absorption-only knot.
  It covers nonemitting absorption, varied density layers, retained phase
  normalization, both headless jobs and refusal of conservative or LTE media.
- `photometric_path::multiple_emitter_reference` combines finite fixtures with
  unequal area, flux, dimming and disjoint spectra over a measured colored
  diffuse receiver. Independent projected polygon solid angles and original
  CIE rows supply incident lux and reflected XYZ. Different direct-light
  strategies and reversed emitter order test the multi-source mixture.
- `photometric_path::uncertainty_reference` measures one unresolved finite
  panel and two disjoint panels with different luminances. Independently
  flattened photopic spectra make each entire path a member of a known
  discrete population. Exact pixel-area probabilities then give both mean
  luminance and population variance, independently of the renderer's sample
  moments. Fixed sample prefixes check absolute cd/m², estimated IID standard
  error and its inverse-square-root sample-count scaling. Single-sample and
  Owen–Sobol results retain no IID error claim. These fixed regression
  tolerances are neither confidence certificates nor total-error bounds.

- `photometric_image::phase_reference` independently integrates narrow forward
  and backward tabulated phase lobes in angle space. For a declared 100 cd
  fixture, 0.05 m⁻¹ absorption and 0.1 m⁻¹ scattering, the two references are
  approximately 2.04655970 and 2.04655733 cd/m². Both lobes lie entirely between
  the original two/four-panel midpoint nodes, which otherwise yield zero and
  apparent convergence. The regressions also check rigid placement, SI scaling,
  shared multi-fixture partition budgets, and unknown discrepancy when no
  actual refinement comparison fits the budget.
- `photometric_image::optical_depth_reference` checks a nearby point fixture
  behind the camera in an optically thick medium. An independent dimensionless
  source-leg integral gives approximately 6495.34889475 cd/m² where the former
  default uniform grids returned almost zero with an apparent tolerance pass.
  The reference also checks SI scale covariance and insufficient-budget
  refusal. Pure panel tests cover numerical range and optical resolution;
  generic-callback regressions retain constant-field and distant-support cases.
- `photometric_image::near_source_reference` checks a 100 cd isotropic fixture
  0.01 m from a 10 m camera ray in a 0.01 m⁻¹ scattering medium. An independent
  angle-space reference gives approximately 23.71377356 cd/m²; the former
  32/64-panel estimates agreed within the requested discrepancy near 11.75 cd/m².
  Distance partitions recover that contribution under the shared cap. Tests
  also check SI/rigid transforms, exact interval budgets and typed refusal.
  Results retain `NumericalKind::NoClaim` and unknown total error: resolving these
  source scales supplies neither a numerical enclosure nor a calibration claim.
- `photometric_image::smooth_beam_reference` checks an approximately 0.954-degree
  full-width-at-half-maximum cosine-power beam against an independent angle-space
  airlight integral. The reference is approximately 0.5734469215 cd/m² while the
  former distance-only grids miss almost the entire beam. The regression also
  checks rigid/SI transformations, tilted optical axes, shared work caps and
  unchanged unknown physical-error fields. The native partition regression uses
  the actual production roots with a separately evaluated radiometric integrand.
- `photometric_path::penumbra_reference` defines independent clipped-emitter
  polygon references across five receiver positions and four direct-light
  strategies. Its comparisons cover probe lux and reflected camera `rho*E/pi`,
  with fixed sampling-regression tolerances and unchanged `NoClaim` reporting.
- `photometric_path::brushed` defines an absolute tilted-normal point-source
  radiance oracle for the existing law `f_R = f_s * (ns.wi)/(ng.wi)`. It checks
  weighted local reciprocity after removing that radiance projection, preserving
  the main shading-normal law and unknown total-error status. These entries
  describe regression coverage, not completed execution results.

`photometric_path::cie171_analytical` independently reconstructs the diffuse
point-source subcase of CIE 171:2006 section 5.2 from the factual setup in the
[DDCI/Lighting Analysts validation report, page 17](https://www.agi32.com/Downloads/TechnicalDocs/Report%20on%20AGI%2032%20validation%20of%20CIE%20171_Compiled_070620.pdf#page=17).
A downward `1000*cos(theta)` cd source sits 3 m above ten horizontal probes.
Expected lux is evaluated from the original expression
`9000 / (9 + x*x + y*y)^2`, combining source and receiver cosines with inverse
square, rather than copying rounded reference values. Both the direct scene
and progressive path-probe APIs are checked, including their unknown total
error and numerical `NoClaim`.

Only this analytical subcase is reconstructed; the other original scenes above
are independent references. Neither the CIE publication nor the validation
report's text, diagrams or reference tables is redistributed. The
[engine author's publication note](https://discourse.radiance-online.org/t/radiance-validation-challenge-cie-171-2006/1981)
distinguishes independently implemented test models from the publication.
This selected case does not establish full CIE 171 conformance or physical
fixture accuracy. The open observer dataset is a separate CIE reference.
Real-fixture accuracy still requires licensed measurements, material/medium
characterization and calibration uncertainties.

### Passive spectral diffuse materials

`Material::SpectralDiffuse { reflectance }` admits
`spectral::SpectralReflectance::constant(rho)` or `::table(knots)` as a native
opaque Lambertian surface. Tables contain 2--16 strictly increasing wavelength
knots in nanometres, with endpoints exactly 380 and 780 nm and finite fractions
in [0,1]. Evaluation is linear in wavelength; malformed data and out-of-band
evaluation refuse without sorting, clipping or extrapolation. Exact zero and
unit reflectance are admitted. Values are stored directly, preserving tiny
reflectance that a rounded `1 - emissivity` representation could erase.

The fixed-capacity table validation/interpolation owner is shared with existing
spectral thermal emissivity, while the public reflectance and emissivity types
retain distinct meanings. Existing thermal identity bytes and arithmetic are
preserved. The passive material has no temperature, LTE source or implicit
emission. Its material identity includes the original reflectance knots;
identities establish declaration/replay binding, not measurement authority.

Direct photometric images, spectral path images, illuminance-probe paths and
the native bidirectional diffuse branches evaluate the original `rho(lambda)`.
The existing cosine-weighted sampler and matched PDF are reused; each diffuse
bounce preserves the wavelength's own reflected fraction. A cached RGB AOV
projection is display metadata and never replaces the physical path spectrum.
Photometry version 5 and path version 11 add this material admission without
changing existing source spectra or the arithmetic of existing material families.

Independent regressions use published CIE weights for two-band illumination,
exact black/white limits and a closed emitting cube with wavelength-dependent
reflectance. The cube's finite-depth reflection series averages powers of the
original spectral reflectance, exposing premature gray averaging. The table
does not certify measurement calibration, spectral resolution/interpolation,
angular Lambertian accuracy or omitted transport. Numerical `NoClaim`, unknown
total error and the separate interpretation of sampling uncertainty remain.

### Borrowed spectral physical media

`PhotometricScene::try_with_medium` accepts `PhotometricMedium`, borrowing the
same `AbsorbingGrid` or `ScatteringGrid` used by the spectral volume integrator.
Original density buffers, mass-coefficient spectra, physical cell bounds,
traversal budgets and cancellation remain authoritative. Direct lux integrates
the 80 wavelength transmittances before the photopic conversion; a colored
medium is not replaced by one attenuation factor on an already integrated lux
value. Gray and vacuum segments retain the homogeneous calculation's bits.

`direct_surface_luminance_along_ray` explicitly returns only fixture light
reflected once by the first opaque surface and transmitted without scattering
to the camera. Both leg transmittances enter its spectral BRDF integral.
The deterministic `LuminanceJob` accepts passive absorption grids but refuses
physical-grid scattering and LTE emission; those require the path estimator
below rather than a zero contribution presented as a complete image.

### Tabulated angular scattering

`volumes::transport::scattering::TabulatedPhase::try_new` accepts 2--4096
strictly increasing `(cos_theta, relative_density)` knots covering exactly
[-1,1]. Values are finite and nonnegative, with nonzero total mass. The density
is linear in propagation-direction cosine, azimuthally symmetric and normalized
with the physical solid-angle measure `2*pi*d(cos_theta)`. It is not a histogram
of angular counts or interpolation in degrees. `PhaseFunction::tabulated`
constructs the shared law used by both gray quadrature and spectral paths.
Scattering strength still comes from gray `sigma_s` or native `rho*kappa_s`;
the table supplies angular redistribution in sr^-1, not an extinction coefficient.

Compensated trapezoid integration normalizes the represented angular law;
analytic quadratic CDF inversion samples it with the same per-steradian PDF.
Zero-density intervals retain their zero response and receive no samples.
Unrepresentable positive mass, unresolved active intervals and invalid input
refuse rather than being assigned a density floor. Floating-point endpoint
repairs stay in the selected positive interval and evaluate the actual returned
direction's PDF. The immutable table is shared through `Arc`; native phase
sampling and point/finite-source/thermal NEE borrow it without per-vertex clones
or allocation. Native medium identities include the normalized PDF and represented
CDF. Existing analytic phase arithmetic, random addresses and identity bytes
remain unchanged. Native scattering semantics is version 6 after the spectral
angular-law extension below.
`PhaseFunction`, `HomogeneousMedium` and `PhotometricMedium` are now `Clone`
rather than `Copy`; clone an owned medium explicitly when constructing two scenes.

Photometry version 4 and path version 10 add this table support. The headless
example accepts `--haze-phase CSV` in gray or spectral modes; the ordinary Rust
interface needs no file parser. Tables are wavelength-independent declared
angular data, not evidence of a measured aerosol or a particle-size model.
Measurement, interpolation, angular-resolution and polarization errors remain
unbounded; normalizing a table cannot strengthen `NoClaim` or turn sampling
standard error into a total error bound. Focused tests use independent linear-law
moments/CDFs, HG refinement, gray and two-layer spectral single-scattering
integrals, finite-panel solid-angle integration, and progressive image/lux replay.

### Spectral angular transport and physical particle scattering

`PhaseFunction::spectral_table` admits 2--16 visible wavelength knots covering
380--780 nm. Each interval is the convex mixture of two admitted normalized
angular laws. Packet continuation samples the uniform mixture of all lane laws
and multiplies each lane by its own phase response divided by that complete
mixture PDF. Rectangle/environment, finite photometric and thermal-source NEE
use the same proposal for MIS, while physical responses remain wavelength
resolved. No phase table is cloned or allocated at a path vertex. Achromatic
laws retain their previous direction/PDF arithmetic. The gray deterministic
photometry owner explicitly refuses spectral phase tables.

`PhaseFunction::rayleigh_depolarized` implements the normalized scalar molecular
phase for the stated unpolarized 90-degree depolarization convention; zero
depolarization retains the existing Rayleigh law. `MassScatteringSpectrum::rayleigh`
evaluates the exact declared inverse-fourth-power wavelength law, rather than
a linear interpolation between its visible endpoints. Extinction admission and
the volume contraction enclosure include these coefficient and angular laws.

`material_state::gas_scattering::MaterialGasScatteringBinding::bind_refractive_field`
connects the same accepted molecular optical sources and gas composition to
`gradient_index` propagation. The caller supplies an immutable smooth physical
density reconstruction and its actual spatial gradient, with an explicit
snapshot/reconstruction identity. The material owner evaluates the
Lorentz--Lorenz mixture index and exact density derivative; the renderer uses
`grad(n) = (dn/d rho) grad(rho)` without introducing a gas preset or fitting a
visual effect. Its dilute visible-band domain, constant-polarizability
approximation and source identities remain attached. Phase and group indices
coincide under that declared nondispersive model. No cell-average smoothing,
composition gradient, measured dispersion, or new thermodynamic law is inferred.
The refractive field has zero true absorption: Rayleigh scattering remains a
separate owner and is not converted to heating. Combining curved propagation
with scene scattering/visibility is not claimed by this adapter.

`MaterialGasScatteringBinding::bind_cell_grid` borrows accepted gas states and
their exact density buffer in the native `VolumeGrid` cell order. Temperature
and pressure may vary by cell while composition and molecular optical sources
remain fixed. Each original state identity survives even when two states have
equal density. The same exact inverse-fourth spectrum and depolarized phase
feed native traversal without per-cell optical tables or a replacement gas law.
This elastic scattering adapter creates no absorption or thermal emission.
Thermochemical gas-field samples can supply the same `RayleighGasState` view.

`volumes::transport::scattering::mie::MieSphere` computes passive homogeneous,
nonmagnetic sphere cross sections, particle-mass coefficients, asymmetry and
the retained unpolarized Mie angular law. Its admitted domain is
`1e-4 <= 2*pi*n_host*r/lambda_vacuum <= 128` and `|m*x| <= 512`, with a real
nonabsorbing host. Downward logarithmic-derivative recurrence is checked from
two starts; small real arguments use a convergent Riccati--Bessel series.
Exact index matching retains zero scattering. Recurrence, multipole and angular
work have explicit budgets and preparation polls cancellation.

`MieSolution::phase_table` adaptively prepares the native matched angular sampler
and reports observed normalization, moment and interpolation discrepancies.
`MieSphere::visible_spectrum` prepares native absorption, scattering and spectral
phase owners and compares the retained wavelength interpolation with fresh
midpoint solves. These diagnostics are sampled numerical evidence, not continuous
error bounds; insufficient angular resolution refuses. Grid density means
dispersed particle mass per mixture volume. A nonunit host requires independently
matching dielectric-interior ray ownership. The model does not claim polarized
transport, particle interactions, nonspherical particles, measured material
dispersion, or an inferred aerosol distribution. A constant-index visible spectrum
does not acquire material dispersion from wavelength sampling alone.

`MiePopulation` admits up to 256 explicitly number-weighted sphere bins in one
real host. Cross sections add by number fraction; mass coefficients divide by
the number-mean particle mass; phase and asymmetry add by scattered power.
Different radii, material densities and complex indices therefore retain their
distinct physical weights. Zero-scattering particles still contribute mass.
Compensated reductions, positive-contribution underflow refusal and per-sphere,
total-multipole and bin budgets bound preparation. The monodisperse case retains
the original sphere coefficients and angular arithmetic. Shared phase/spectral
preparation reduces populations to the existing native tables, so transport
vertices do not loop over radius bins or solve Mie recurrences.

`ParticleSizeDistribution` accepts explicit radius/count bins or a lognormal
number distribution on a caller-declared finite radius interval. Lognormal
weights integrate each log-radius bin's probability; number and volume tail
fractions remain separate from angular and spectral diagnostics. They do not
bound omitted optical scattering. Radius quadrature needs independent refinement;
the constructor never clips particles to the Mie domain. Population regressions
cover monodisperse reduction, mixed density/composition, scattered-power phase
weights, index-matched mass dilution and truncated Rayleigh-moment refinement.

`material_state::particle_scattering::MaterialParticleScattering` binds an
accepted material's explicit scalar bulk-density property, its matching optical
state and current radius/number distribution. Surface and particle consumers
share one conversion of conductor n/k, tabulated dielectric, Cauchy and Sellmeier
sources. Dielectric bulk intensity absorption becomes `k = alpha*lambda/(4*pi)`
with vacuum wavelength in metres; RGB-derived absorption refuses. Every spectral
knot and independent midpoint uses the source's wavelength-dependent index.
Temperature-dependent source values enter when the material is resolved and
bound; a changed state or particle geometry requires rebinding.

`MaterialParticleSpectrum::bind_ambient_grid` supplies native absorption,
scattering and optional LTE emission from the prepared spectrum and borrowed
dispersed-particle density. It admits only an index-one ambient host, enforces a
caller-declared particle volume-fraction ceiling, and requires LTE temperatures
to equal the accepted material temperature. These checks do not certify dilution,
independent-particle accuracy, unmeasured material data or continuous spectral
error. No material, radius or Mie solve runs at a path vertex. Regressions cover
all four optical source families, temperature-dependent source changes and
two-cell Beer--Lambert transport through the original native grid.

Tests include published complex Mie coefficients and efficiencies, Rayleigh
size/wavelength limits, passive energy balance, index matching, angular
normalization and first moments through size parameter 128, spectral refinement,
per-lane phase moments, disjoint angular support, and achromatic replay. Native
path regressions additionally exercise wavelength-resolved single-scattering
integrals and a conservative multiple-scattering environment furnace.

### Calibrated spectral fixture paths

`PathLuminanceJob` connects the same `PhotometricSource` states to the original
spectral surface and volume path integrator. Fixture connections are evaluated
at each opaque surface and volume interaction. These point sources have delta
position distributions: the continuous BSDF/phase hit technique cannot reach
them, so their separately summed next-event terms have unit MIS weight. They
do not alter the continuous-light sampling mixture or its random draws.

The admitted scene uses native meshes/rigid mesh instances and spectral
materials, with the finite-source glass extension described below. Relative
area/environment/material emission refuses at this SI
boundary. Borrowed physical grids supply spectral attenuation, scattering and
their already calibrated LTE radiance. Existing constructors retain the pinhole
`Scene::camera`; `PathLuminanceJob::new_with_camera` instead borrows an explicit
native `PhysicalCamera` and ignores the unused legacy camera. Calibrated hero
wavelengths use the same 80-bin midpoints and tabulated observer as source
normalization, and XYZ accumulation uses 683 rather than
the legacy display-white normalization. Y therefore reports cd/m². The old
entry points pass no fixture context and retain their previous arithmetic.

The explicit physical camera supports focal length plus active vertical sensor
height, axial focus distance, and circular or regular bladed pupils through the
original thin-lens ray generator. Image dimensions declare aspect ratio. Lens
draws use the native independent Philox camera-lens domain at the same absolute
pixel/sample address; they do not consume wavelength, fixture, phase or BSDF
draws. An equivalent legacy-compatible pinhole preserves the original image
bits. The complete sampled ray starts at its actual pupil point, including its
medium attenuation, geometry visibility and subsequent spectral transport.
Every pupil origin uses the shared admitted dielectric membership and nesting
rules above, including a camera placed inside glass. An unresolved boundary
point or unsupported containing topology refuses before publication.

These values are pixel- and pupil-averaged luminance in cd/m², without an
aperture-area or f-number brightness multiplier. Focal length, sensor size,
focus and pupil shape affect the measured rays; ISO and exposure compensation
remain metadata and do not scale absolute luminance. Sensor irradiance,
photographic exposure, vignetting, diffraction, camera animation and a shutter
integral are not supplied by this static calibrated entry. Whole-path IID
standard error includes pupil sampling variability, while the existing unknown
total-error and Owen–Sobol no-IID-claim policies remain in force. The borrowed
camera stays fixed across progressive passes, cancellation and worker counts.
Fixture-path semantics version 7 adds this explicit camera admission while
retaining the existing constructor streams and source-connection policies.

The default point-fixture policy remains the complete source sum. The optional
receiver-dependent connection policy below samples that sum without changing
the underlying source, material or medium models.

The image advances in row-major sample passes and reuses the existing Welford
accumulator. Each complete path and moment update commits atomically. IID
estimates expose a sample-based standard error after two paths; Owen–Sobol
camera sampling exposes no IID standard-error claim. The blocked full-path
Sobol mode refuses here: its repeated base dimensions have not qualified for
these spectral physical references. A finite vertex ceiling omits higher-order
transport, so neither sampling error nor apparent convergence bounds total
image error. Fixture calibration, observer/discretization, geometry and that
finite-depth remainder remain unknown.

Every physical path segment and source connection respects the declared metre
range. A geometry miss uses the finite physical-grid exit, preserving its
foreground transport without an arbitrary truncated tail. Scattering segments
refuse point-source exclusion-radius crossings before free-flight sampling.
The `photometric_stage --paths` example shares ordinary fixture/scene data with
the direct calculation and renders reflections plus multiple scattering through
a borrowed heterogeneous, spectrally absorbing field.
Appending `--camera-lens 50 24 2.8 12` supplies a 50 mm lens, 24 mm active sensor
height, f/2.8 circular pupil and 12 m axial focus distance to that same path job.

### Surface illuminance with reflected and volume transport

`PathIlluminanceJob` borrows the same scene, fixtures and native physical grids
and accepts `IlluminanceProbe` arrays or `CalculationGrid` without a probe-array
copy. `Scene::camera` is unused. `IlluminanceSampling` fixes the sample count,
maximum vertex depth, sampler and seed independently of image dimensions.
The probe measures incoming light without applying a receiver material.

The original direct solver supplies the unscattered fixture term once per
probe. Cosine-weighted hemisphere paths supply `pi * L` in lux using the shared
spectral integrator and its existing surface/volume interactions. The measured
hemisphere consumes no vertex depth. `continuous_lux` includes reflected and
scattered fixture light and physical LTE volume emission; delta point-source
directions occur only in `direct_lux`. Their sum is `illuminance_lux`, without
double counting or converting an outgoing image back to an assumed irradiance.

The initial geometric ray is offset outward, while physical-medium transport
retains the exact receiver endpoint, including adjacent thin cells. Both poles
of the probe hemisphere use the existing stable orthonormal basis. Calibrated
mesh and rigid-instance paths discover the actual first hit before enforcing
the caller's metre range; the legacy instance range is unchanged.

Probe-order sample passes publish complete moment updates. Direct lux and the
first successful path publish together; a failed path exposes no incomplete
probe. The original Welford accumulator estimates IID sampling standard error
in lux after two paths. The exact direct term contributes no sampling variance.
Owen–Sobol results carry no IID standard-error claim. Total error remains
unknown. Finite-depth model truncation is unbounded except in the restricted
diffuse cases below, and physical-model uncertainty remains unknown.
`photometric_stage --total-csv` writes positions, direct/continuous/total lux,
estimated IID standard error, sample count, explicit unknown total bounds and
a named sampling-uncertainty state for a stage grid using heterogeneous
spectral haze and reflected light.

### Interpreting path sampling uncertainty

Both path jobs expose `calculation_context()`, returning the actual admitted
sampler, seed and maximum surface/volume vertex depth, together with the
published observer identifier and photometry/path semantics versions.
`PathCalculationContext` records declared settings for result interpretation;
it is neither a complete scene/replay identity nor an accuracy certificate.
Geometry, fixture, material and medium inputs remain the caller's retained
scene data. A named observer does not establish fixture calibration.

`sampling_uncertainty_cd_m2(pixel)` and `sampling_uncertainty_lux(probe)` return
`SamplingUncertainty`: `InsufficientIidSamples`, `EstimatedIidStandardError`,
or `CorrelatedSamples`, with the accepted whole-path count. The standard-error
value has the method's named SI unit. Uncomputed and out-of-range outputs
return outer `None`. The typed estimate retains the existing Welford value
bit-for-bit; legacy optional standard-error fields remain available. This
distinguishes a one-path IID result from an Owen–Sobol result without inventing
an IID estimate for correlated samples. Zero estimated sampling error leaves
numerical `NoClaim` and any finite-depth model bound unchanged; calibration,
discretization and total error remain unknown.

The headless example derives path metadata from these job accessors. Its total
lux CSV appends `sampling_uncertainty` to identify the reason for a missing
standard-error value; path-image summaries report an unavailable estimate
when no IID values exist instead of substituting zero.

### Restricted diffuse-vacuum depth remainder

`PathLuminanceJob::diffuse_vacuum_remainder()` and the corresponding probe-job
accessor return an optional `DiffuseVacuumRemainder`. Admission requires
explicit vacuum, no active point fixtures, only `Material::SpectralDiffuse`
surfaces and uniform Lambertian finite emitters without projector bindings.
Every primitive participates, including hidden or unlit geometry. Let `q` be
the largest reflectance knot over all materials; the piecewise-linear laws
then satisfy `0 <= rho(lambda) <= q < 1`. A white surface with `q = 1`, any
other material, a physical medium (even an empty grid) or any projector binding
returns `None` for this metadata without rejecting an otherwise valid job.

For this ideal model, the positive diffuse reflection operator has maximum
spectral-radiance norm at most `q`. A conservative source envelope sums the
admitted emitter radiance laws. Its photopic ceiling `S` is evaluated from
each retained emitter coefficient, source-spectrum bin, published observer Y
weight, 5 nm bin width and 683 lm/W factor. It does not assume that a stored
spectrum's floating-point photopic normalization is exactly one. The omitted
reflection orders at native depth `D` therefore satisfy
`B_L = S * q^D / (1 - q)` in cd/m² and `B_E = pi * B_L` in lux. BSDF-only
paths retain orders zero through `D - 1`; terminal next-event estimation can
include an additional order, making this common bound conservative. A probe's
initial hemisphere direction does not consume a transport vertex.

The bound calculation uses outward interval arithmetic for the spectral sum,
geometric tail and mathematical-pi conversion. Exponentiation takes
logarithmic work even at the largest native depth. Exact black or zero-source
tails can be zero; a positive underflow retains a nonzero outward upper bound.
A bound that cannot be represented finitely returns `None`. Getters retain
the depth, maximum reflectance, source ceiling and both SI bounds. This
immutable scene/depth metadata is available before sampling and does not
change as a job advances.

This bounds only the ideal model's omitted reflection orders using its
admitted radiance coefficients and represented spectral laws. It does not
enclose sample deviation, floating-point transport/intersection arithmetic,
errors in geometry or source flux/area normalization, spectral discretization,
fixture calibration or disagreement with physical materials and fixtures.
It is not an interval around a sampled result. In particular, a zero model
tail does not change `NumericalKind::NoClaim`, estimated IID sampling error or
the `None` total-error fields. Point-source, projector and non-diffuse
truncation bounds remain unavailable.

### Nonemitting participating-media depth remainder

Both path jobs also expose `diffuse_transport_remainder()`, returning an
optional `DiffuseTransportRemainder` for the same finite emitters and diffuse
surfaces in vacuum or a nonemitting native absorption/scattering grid. This
broader accessor does not change `diffuse_vacuum_remainder()`: every physical
grid still returns `None` from the older accessor, even when its density is
zero. The source and surface envelope is shared by both calculations.

For a nonemitting free flight with transmittance `T`, the continuation operator
has norm at most `T*r + (1-T)*w`, where `r` bounds surface reflectance and `w`
bounds spectral single-scattering albedo times the phase integral. Therefore
`q = max(r,w)` bounds the combined continuation; adding `r+w` would incorrectly
treat surface and volume events as independent alternatives on one flight.
The source ceiling and the outward `S*q^D/(1-q)` and `pi*S*q^D/(1-q)` formulas
remain unchanged. Absorption alone cannot increase this bound. Density may be
heterogeneous: the native model multiplies absorption and scattering by the
same density, so density cancels in their ratio.

The albedo maximum is enclosed at the union of absorption and scattering
spectral knots. Between consecutive combined knots, the ratio is fractional
linear and monotone wherever extinction is positive; shared zero endpoints
are covered by their neighboring limits. Coefficient interpolation and the
ratio use outward intervals. The phase integral accounts for retained native
analytic coefficients or the trapezoidal integral of retained tabulated
densities; normalization rounding is not assumed absent. Its gain envelope
is at least one, so downward normalization roundoff cannot create artificial
absorption. The final `q` must remain strictly below one without clipping.

Metadata retains the maximum surface reflectance, scattering gain, combined
contraction, depth, source ceiling and SI bounds. Zero actual scattering has
zero volume gain. Thermal emission, conservative scattering without a strict
global contraction, unsupported materials/sources or nonrepresentable bounds
return `None`; these cases require further source or escape arguments. The
bound still concerns only omitted ideal-model orders. Sampling deviation,
floating-point transport and normalization error, geometry, spectral
discretization and calibration remain outside it, with numerical `NoClaim`
and unknown total-error fields unchanged.

### Native CPU execution of calibrated paths

Both calibrated path jobs expose `advance_with_execution`, using the existing
`RenderExecutionConfig`, native `RenderWorkerPool`, parked crews and caller `Cx`.
One crew serves a complete advance call. Work tiles contain contiguous absolute
sample addresses, with `tile_width * tile_height` paths per tile. Worker count,
scheduling weights and tile size do not enter the spectral sample identity.

Each wave privately stages at most 65,536 paths, also bounded by twice the worker
count times the tile area and by the caller's requested work. Its moment updates
replay in ascending absolute address order using the original accumulator.
Only a fully successful wave publishes; cancellation or numerical refusal keeps
the previously accepted prefix. Subsequent serial or parallel calls resume from
that prefix. Direct probe lux remains cached with the first accepted sample.
Completed sample sets retain the serial XYZ, lux and estimated-error bits.

The existing operation-memory lease accounts for retained moment capacity,
bounded wave samples, private candidate moments and native executor metadata.
Advancing a job does not copy its geometry, fixtures or fields. Prepared emitter
proposal tables, allocator bookkeeping and native thread stacks remain outside
this payload accounting. Memory/executor
refusals retain their structured errors; a physical transport refusal is not
relabelled as ordinary cancellation.

The example accepts `--paths --workers N` and `--total-csv --workers N`; the CSV
schema and physical values are shared with serial execution. This execution
bridge changes scheduling and publication granularity, not the transport model
or its finite-depth and uncertainty limitations.

### Finite calibrated surface emitters

`PhotometricAreaEmitter` binds total lumens, a linear dimmer and the original
80-bin source spectrum to a native mesh primitive. `PhotometricEmitters` combines
these bindings with existing point fixtures. Both path jobs accept this set;
static rigid instances share the original geometry and material abstractions.
One-sided luminance is `flux/(pi*area)` and two-sided luminance is
`flux/(2*pi*area)`, using the actual world-space triangle area and winding.
Polygonal apertures use their represented mesh area, without claiming an exact
curved aperture. Emission adds to the existing material's reflected radiance.

`PhotometricAreaEmitter::filtered` applies a passive 80-bin spectral transmission
using the same normalization owner as point-fixture filtering. The photopic
transmission scales both clear and emitted lumens; only the surviving spectrum
is normalized. Gray filters retain the original spectrum, opaque filters
produce zero flux, and repeated filters multiply their spectral transmission.
The filter is uniform over the source and independent of direction. Attaching
or replacing achromatic projector optics retains the filtered clear calibration;
it does not restore absorbed power. This is not a spatially varying color gobo
or an angle-dependent dichroic model.

Finite-source proposals select emitters by flux and triangles by area. Source
hits and surface/phase/probe proposals use complementary MIS with the original
BSDF/phase density and physical-medium survival probability. A terminal path
vertex gives its direct proposal full weight. Visibility tests the actual
endpoint; another face of the emitting primitive can occlude it. The complete
source spectrum participates in the original absorption and scattering paths.

Probe `continuous_lux` includes directly observed finite sources. Its area NEE
and hemisphere-path contributions are combined before moment accumulation, so
the sampling-error estimate describes the complete sample. The cached
`direct_lux` term remains point-fixture-only. Point-only paths retain their
original sample streams. `BsdfOnly` admits area-only scenes but refuses active
point fixtures, whose delta directions cannot be sampled by a continuous BSDF.

The stage example accepts `--paths --panels` or `--total-csv --panels`, optionally
with `--workers N`. Default bindings remain uniform Lambertian patches. The
optional projector model below supplies a position/direction-dependent source
law. Calibration, finite depth, spectral and geometric uncertainty retain
their explicit unknown total-error bounds.

### Native glass and brushed surfaces in calibrated paths

Both path jobs accept the existing smooth/rough `Material::Dielectric` for
finite-source scenes. Closed outward-oriented native meshes and static rigid
instances retain the spectral renderer's boundary stack, wavelength-dependent
Snell/Fresnel transport, internal Beer–Lambert absorption and radiance eta
factors. Ambient physical grids are displaced inside homogeneous glass. An
opaque finite emitter may sit inside a separate glass solid; the emitting
primitive itself cannot also be a dielectric boundary. Camera and calculation
points must start in ambient space. General nested-medium initialization at
interior sensor positions is not supplied by these jobs.

Rough-interface source proposals evaluate native reflection/transmission
response and continuation density separately at each packet wavelength. The
shadow direction selects the actual incident or transmitted medium. Smooth
interfaces use the native refractive walk, including its spectral/delta splits;
there is no calibrated area-source manifold proposal, so those source hits
retain full weight. Optional legacy rectangle/environment connector discovery
is skipped for calibrated jobs, which cannot contain those sources.

Active point-fixture-plus-glass scenes refuse explicitly. Their delta-position
emitters cannot be found by continuous refractive walks; a general specular
point-connection estimator remains absent. The public direct-only
`PhotometricScene` retains its opaque-material admission rather than presenting
straight-through or missing refractive transport as total illumination.

`Material::BrushedConductor` uses its existing body-local tangent frame for
calibrated point/area connections and BSDF sampling, including rigid placement.
The original optical table, surface detail and medium-index conventions remain
authoritative. Surface-detail identity v2 rotates the bitangent together with
the perturbed shading normal, keeping local/world reflection frames orthonormal.
Brushed reflection sampling and its direct-light MIS density both exclude
directions below either the geometric or shading horizon. Focused regressions
cover admitted nonzero tilts, unit reflected rays, and matching sample densities.
Brushed BSDF semantics version 2 applies the radiance projection
`f_effective = f_shading * |ns.wi| / |ng.wi|`, following PBRT 3e
§16.1.3, “Non-symmetry Due to Shading Normals.” Native and calibrated
consumers retain their geometric incoming cosine. This removes artificial
grazing-view gain while preserving the reciprocal local GGX law; the effective
radiance BSDF itself is generally non-symmetric. The existing bidirectional
owner still refuses brushed conductors rather than treating this radiance law
as its importance adjoint. Furnace regressions cover the bounded reflected
radiance and an independent angular-quadrature reference. Decorative normal
mapping remains a shading approximation, not resolved microgeometry.
The new Rust regressions remain natively unexecuted in the implementation environment.
Reference: [PBRT, The Path-Space Measurement Equation](https://www.pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation).
This does not promote representative optical data or decorative normal/mark
parameters into measured material claims. Path semantics version 4
adds these material paths; previously admitted opaque sample streams remain
unchanged. Sampling estimates still exclude finite-depth and physical-model
error, including geometric-optics and surface-model approximations.

### Finite-aperture thin-lens projectors

`PhotometricAreaEmitter::with_projector` binds `ThinLensProjector` to the same
native mesh source. Its optical frame is relative to the original primitive;
static rigid placement moves the aperture and optics together. The pupil must
lie in optical `z=0` with winding toward `+Z`. Plane/normal checks allow only
coordinate-roundoff tolerance, which is not a geometry-error certificate.
The original material still reflects incoming light, and source rays use the
same spectral paths, shadows, absorption, scattering and glass as other finite
emitters. No independent projection renderer is introduced.

This is an explicitly ideal exit-pupil light field: clear output carries
uniform ray power over pupil area `A` and rectangular gobo area `G`. With gobo
distance `s`, positive thin-lens focal length `f`, pupil point `a`, and outgoing
slope `u = omega.xy / omega.z`, the source evaluates
`g = a*(1-s/f) - s*u`. Thus the real image distance satisfies
`1/f = 1/s + 1/q`, with image magnification `-q/s`. Different pupil points
produce different rays outside that focus plane; defocus is determined by
the actual aperture geometry, not an angular texture blur.

`GoboTransmission` is a centred rectangular array of piecewise-constant gray
transmittances in `[0,1]`, in metres, row-major from negative to positive Y.
Outside the rectangle transmission is zero. If `M(g)` is its value, emitted
luminance is `Phi_clear*s²*M(g)/(A*G*cos(theta)^4)`. The projected-plane
Jacobian `dOmega = cos(theta)^3*dg/s²` makes
`L*cos(theta)*dA*dOmega = Phi_clear*M(g)*dA*dg/(A*G)`.
Consequently `flux_lm()` reports `clear_flux_lm()*mean(M)` after the dimmer;
blocked power is removed, never redistributed into the transmitting texels.
These are total outgoing source lumens before later scene-medium/glass losses.
The original spectrum scales this luminance without changing its color.

Both source hits and direct source proposals evaluate this same light field.
Emitter selection uses transmitted flux; the proposal still samples the actual
mesh aperture by area, so blocked gobo rays contribute zero with the original
endpoint density. Its complementary MIS weights therefore stay valid. Small
transmitting features or severe defocus can require more samples; no adaptive
pattern-resolution or error-per-second guarantee is claimed. Whole path
samples, cancellation, deterministic parallel accumulation and estimated IID
standard error retain the existing job semantics. Path semantics version 5
adds this source model without changing existing Lambertian/point streams.

The model does not solve lamp illumination or a multi-element lens train and
does not include aberrations, diffraction, wavelength-dependent focal length,
physical prism refraction. Its thin-lens/paraxial law and uniform pupil/gobo ray-power
assumption are declared approximations, with unknown physical-model error.
Finite input and radiance representability are checked at admission; numerical
agreement with the ideal model does not certify a manufactured fixture.
`photometric_stage --paths --projector-focus 5 --workers 4` exercises this
source together with native material reflection and heterogeneous spectral
haze. The same option works with `--total-csv` for surface lux.

### Calibrated multi-beam projector optics

`ParaxialBeamSplitter` supplies at most 32 `ProjectorBeam` branches to
`ThinLensProjector::with_beam_splitter`. Each branch declares a dimensionless
optical X/Y slope offset `d_k` and a passive flux fraction `t_k`. Fractions are
finite and nonnegative, and their compensated sum must not exceed one. They
are not renormalized. An absorbing splitter with zero total transmission is
valid. Splitter replacement replaces the optical state rather than multiplying
it into a previously applied splitter.

For outgoing slope `u = omega.xy/omega.z`, branch `k` reads the same gobo at
`g_k = a*(1-s/f) - s*(u-d_k)`. The shared source evaluator uses
`L = Phi_clear*s²*sum(t_k*M(g_k))/(A*G*cos(theta_out)^4)`. Translating slopes
does not change the pupil/gobo ray-power measure, so total emitted lumens are
`Phi_clear*mean(M)*sum(t_k)`. Overlapping branches add their actual light;
absorbed power is not redistributed among the surviving branches. At the
existing image distance `q`, each projected copy shifts by `q*d_k` while
retaining magnification `-q/s` and the same geometric pupil defocus.

Source hits, surface/probe connections and volume connections evaluate this
one light field. The existing native mesh-area proposal and complementary MIS
remain in force; all bounded branches are evaluated without another random
draw or a separate branch-selection estimator. Admission includes the shifted
angular support in its finite source-radiance bound. An unsplit projector keeps
its original arithmetic. Reflections, absorption and multiple scattering still
use the original material and physical-grid owners.

This is an ideal, achromatic paraxial exit-pupil model. It does not resolve a
faceted prism's spatial pupil partition, wavelength-dependent Snell refraction,
Fresnel/internal reflections, diffraction or aberrations. A supplied branch
pattern and throughput are declared source data, not inferred manufacturer
measurements. Sampling uncertainty still excludes those model/calibration
errors and the finite-depth remainder. Path semantics version 8 adds the
splitter and finite-source filtering while preserving existing unsplit streams.
The headless example accepts `--projector-focus 5 --projector-split 3 4` for
three beams on a 4-degree ring with 90% total splitter transmission, with either
`--paths` or `--total-csv` and the existing worker/progressive execution.

### Passive paraxial projector diffusion

`ThinLensProjector::with_diffuser` accepts `ParaxialBoxDiffuser`, declaring
positive optical X/Y half angles below pi/2, a scattered fraction `alpha` and
an absolute transmission `tau`, both in [0,1]. Scattered rays add independent
uniform slope offsets in `[-tan(half_angle),tan(half_angle)]`; the remaining
fraction keeps the original outgoing slope. This is an ideal achromatic
exit-pupil ray-power model, not a reciprocal surface BSDF or measured frost.

With the original backprojected coordinate `g`, the source mask becomes
`tau*((1-alpha)*M(g)+alpha*average(M(g+s*delta)))`. The rectangular average
integrates the represented piecewise-constant gobo, including its opaque
exterior. The existing outgoing cosine Jacobian remains in force. Thus total
lumens are `clear_flux*mean(M)*splitter_transmission*tau`; widening the beam or
changing its scattered fraction does not restore absorbed power. At focus
distance `q`, diffuser half widths on the receiving plane are `q*tan(half_angle)`.
Pupil defocus and splitter translations compose with this redistribution.
Replacing the diffuser replaces its optical state rather than compounding loss.

A private positive-sum rectangular range tree integrates the gobo without
subtracting large nearly equal prefix sums. Its optional immutable cache is
shared by clones, capped at `4*MAX_GOBO_TEXELS` doubles including power-of-two
padding, and allocated fallibly only for active diffusion. Each evaluation
uses bounded logarithmic range queries with no heap allocation or extra random
draw. The clear and identity-diffuser paths retain the existing arithmetic.
Source admission includes broadened angular support in its radiance bound.
Unrepresentable positive contributions refuse; numerical rounding is not certified.

Source hits and surface, probe, and volume connections evaluate this same law
through the existing finite-emitter estimator and native materials/media. Path
semantics version 12 adds diffusion; cancellation, progressive ordering and
IID sampling-error interpretation are retained. The headless example accepts
`--projector-focus 5 --projector-frost 1.5 3 0.8 0.9` in either path mode,
including filters, splitters, physical cameras and heterogeneous haze. Physical
diffuser calibration, wavelength-dependent scattering, diffraction and the
paraxial/model discrepancy remain unknown; total error retains `NoClaim`.

### Receiver-dependent point-fixture sampling

`PointFixtureSampling::Importance { connections_per_vertex }` is an opt-in
policy on both calibrated path jobs, selected with `with_point_sampling` or
`set_point_sampling` before any path has committed. The retained policy is
queryable; invalid counts or changes after progress refuse. The mutating setter
preserves the job and its accepted prefix on refusal. `All` remains the default
and retains the original source sum, arithmetic and sample streams.

At a surface vertex, proposal weight is the actual fixture intensity toward
the receiver divided by squared physical distance, times receiving cosine.
At a volume vertex the phase response replaces that cosine. Thus off-axis
beam structure, fixture orientation and receiver position affect selection.
All active source distances, exclusion radii and angular intensities are
examined, including sources whose proposal response is zero. Whole-segment
point-source exclusion in scattering media remains unchanged.

Between 1 and 32 independent weighted reservoirs select connections with
replacement. For proposal marginal `p_j`, a whole path includes the average
of `C_j/p_j` over its connection draws. `C_j` uses the original visibility,
spectral source, physical Beer attenuation and BSDF/phase response. Point
sources still have unit MIS weight because continuous paths cannot hit their
delta positions. Duplicate selected identities share one expensive connection
evaluation with the corresponding multiplicity. Finite-area/projector source
proposals and their complementary MIS are unchanged.

Reservoir storage is fixed-size and path-local; no per-vertex heap allocation
or shared mutable sampler is introduced. Independent Philox domains depend on
logical path, vertex and source/draw identities, without consuming camera,
wavelength, BSDF, phase or area-source draws. Same-ISA serial/parallel and
progressive replay retain identical values for the same policy and seed.
Positive importance that cannot retain represented proposal support refuses
explicitly; callers can use `All` for such extreme mixtures. Finite arithmetic
and random-word resolution are not rigorous numerical error enclosures.

The direct unscattered point-fixture lux cached by each illuminance probe
always uses the complete original sum. Importance sampling affects only its
reflected/scattered paths, and point contributions to camera paths. All
connection draws belong to one path sample, so the existing IID standard error
includes their selection variance. Owen–Sobol retains its no-IID-error-claim
policy, and total model/transport error remains unknown.

The policy bounds costly shadow, medium and spectral shading evaluations by
the connection count. It still performs `O(N*K)` inexpensive reservoir work for
`N` sources and at most `K=32` draws; it is not a light-tree or sublinear source
traversal claim. The proposal does not predict visibility, colored attenuation
or sharp material lobes, so variance and error per second remain scene-dependent.
Path semantics version 6 adds this policy. The headless stage example accepts
`--point-connections K` with either path mode; omitted means `All`.

### Prepared geometry for calibrated transport

Each immutable photometric scene prepares a spatial index over its admitted
native meshes and rigid mesh instances. Calibrated path jobs share this index
with their direct-lux scene. Direct/image queries, point-fixture connections,
finite-emitter visibility, path continuation and thermal-volume source queries
use the same original geometry and material witnesses.

The index is a broad phase: surviving primitive indices are visited in their
original order through the existing native intersection and tie rules. This
retains closest-hit arithmetic, coincident mesh/instance behavior, material
frames and per-hit backend diagnostics. The original mesh BVHs remain borrowed;
no mesh copy, per-query index allocation or random draw is introduced.

Admission and traversal poll cancellation. Metadata allocation is fallible.
Tiny scenes and scenes beyond the bounded candidate-bitmap capacity use the
linear traversal; uncertain bounds remain candidates, and rays outside the
index's numerical admission use the native traversal. Bounds reuse the rigid
motion enclosure and outward-rounded mesh slab test with additional broad-phase
roundoff padding. This is not a numerical error certificate.

The physical estimator, source-selection policy and reported uncertainty are
unchanged. Bitwise comparisons against the native linear traversal cover direct
lux, image/probe paths, rigid instances, endpoint/tie cases and progressive
parallel replay. The paired dense-stage benchmark reports index construction
separately from repeated connection evaluation, with identical lux outputs.
