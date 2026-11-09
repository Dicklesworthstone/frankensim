## Spectral diffuse cooling on the native enthalpy timeline

`thermal_surface::diffuse_exchange::DiffuseExchangeResponse` integrates the
original piecewise-linear spectral emissivity against the shared stable SI
Planck difference. Positive four-point quadrature in log wavelength splits at
every source knot; the diffuse hemisphere factor is analytic. Prepared grids
are shared across state-dependent emissivity refreshes. Evaluation allocates
nothing and preserves heating/cooling reversal and the equilibrium tangent.

`heating::diffuse_exchange::NativeDiffuseExchange` resolves the SAME native
temperature/enthalpy-dependent law for applied heat and its accepted renderer
surface. The existing thermal observer owner advances storage, conduction,
convection and events and applies each accepted-start radiation rate once.
Explicit budgets bound source work, spectral work and retained quadrature.
The frozen-emissivity secant supplies the native step restriction; state-law
and temporal accuracy still require independent time refinement.

Unknown wavelengths are never filled with guessed emissivity. Integrated Wien
and Rayleigh--Jeans majorants, capped by complete blackbody power, bound omitted
exchange using the hotter body/bath. Area-safe logarithmic evaluation uses the
same SI constants and outward floating-point padding. This is an analytic model
bound, not a formal interval or material-uncertainty certificate. Each source
may require an absolute omitted-watt allowance; partial-band heat remains
explicitly labeled when no allowance is requested. The bound adds no heat and
does not include wavelength-quadrature or temporal error.

Three tests of the exact production tail helper passed against independent
Planck integration, physical scaling and numerical extremes, linked to the
actual pinned fs-math implementation. Three response and four native integration
regressions cover spectral refinement, thermal balance, latent spectral changes,
coverage refusal and ambient events. These integration tests remain unexecuted
because the checkout lacks required constellation siblings; formatting passed.

## Nonuniform incandescent emission and physical exposure

`ThermalVisibleEmission::observe` combines original diffuse thermal surfaces
only after integrating each source spectrum. Explicit per-face areas yield
absolute visible watts, lumens and flux-scaled XYZ; hot spots and differing
spectral emissivities survive. No representative temperature or white
normalization enters the physical result.

`ThermalVisibleExposure::integrate` clips explicitly held accepted emission
intervals against the existing rectangular shutter. It returns visible joules,
lumen-seconds, integrated XYZ and mean flux. This is gross emitted energy;
camera collection still belongs to the original geometry/transport owner.
The held reconstruction is an explicit time approximation and adds no heat.

Both paths use bounded work and cancellation, compensated sums and no heap
allocation. Unsupported spectral coverage, nonpositive areas, time gaps,
overlaps, missing shutter coverage and numerical overflow/underflow refuse.
Three focused regression groups define hotspot/area scaling, clipped pulse
exposure and refusal behavior. Pinned formatting and source/API checks passed;
native integration execution remains pending the missing sibling dependencies.

# CONTRACT: fs-render

Unbiased spectral path-tracing core: the verifiable Monte-Carlo foundations.

## Purpose and layer

Layer L5 (LUMEN). The Monte-Carlo core depends on deterministic `fs-math`;
the default chart backends consume lower-layer `fs-evidence`, `fs-exec`,
`fs-geom`, and `fs-rep-nurbs`. Optional differentiable, volume, and tracer
surfaces add their declared lower-layer dependencies (`fs-ad` for the
differentiable lift). Pure Rust throughout.

## Public types and semantics

All native GGX reflection uses height-correlated Smith masking-shadowing,
`G2 = 1/(1 + Lambda_o + Lambda_i)`, with directional `G1 = 1/(1 + Lambda)`.
The same facet height governs visibility from both directions; the two
marginal masking probabilities are no longer multiplied as independent events.
This is the statistical single-scattering law of
[PBRT 4e section 9.6.2](https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory).
Legacy spectral GGX, exact conductor and body-locked anisotropic conductor
evaluation share it with both hemispheres of rough dielectric transport.
Fresnel, coherent films, the original visible-normal proposal and its PDFs,
rejected-facet probability mass and radiance eta factors retain their owners.

Polarized metal samples and finite-direction connections use the same
`G2/G1(wo)` weight, including the existing adjoint for tilted shading normals.
That ratio is evaluated directly from admitted G1 probabilities; a weak
Fresnel channel or narrow-lobe density is never divided out. Shared G2
evaluation orders its bounded factors to retain representable grazing power.
Work is constant with no allocation or additional random samples. Source
semantics change to tracer 9, dielectric 20, conductor 2 and brushed 3.
This still omits inter-facet multiple scattering: the missing single-scattered
energy is not a physical absorptance or thermal debit. Smooth conductor
emission/absorption and the separate EON rough-diffuse thermal law retain their
original constitutive owners. No rough-metal thermal closure is inferred.
Focused G1/G3 controls compare independent slope-form Smith/Walter equations,
grazing reflection, sampled reflection/transmission weights, reciprocity and
angular power integrals; native scalar/polarized controls exercise shared use.

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
- `rough_diffuse::RoughDiffuseSurface` and `Material::RoughDiffuse` provide
  isotropic, depolarizing energy-preserving Oren--Nayar (EON) reflection.
  The normalized roughness in `[0,1]` is a model parameter, not a measured
  height or slope angle; spectral table values are diffuse-facet
  single-scattering albedos. The reciprocal multiple-return closure uses the
  exact analytic projected albedo, with a factored smooth limit and stable
  grazing deficit. White facets conserve integrated reflected power, colored
  facets absorb on additional returns, and zero roughness is exactly the
  existing Lambertian law. Scalar camera/BDPT, calibrated photometric and
  native polarized camera/connection paths evaluate the same angular law at
  the original wavelength. Native sampling mixes the cosine hemisphere with
  a view-aligned uniform-elevation, cosine-azimuth lobe. The unchanged BRDF's
  positive grazing term determines the mixture mass; the cosine component
  retains full support. Two original random coordinates select both the branch
  and direction without allocation or rejection. Scalar, thermal and polarized
  paths use the same proposal. Reciprocal connections and BDPT evaluate the
  actual view-conditioned reverse density, which need not equal the forward
  density. Zero roughness retains the original Lambertian samples and PDFs.
  Material identities bind sampling semantics version 2 for replay/checkpoints.
  A single path weight above one is not an integrated energy gain.
  G0/G1/G3 checks cover admission, reciprocity, exact limits, independent
  projected-angle white/colored furnace quadrature and native consumers.
  Actual-source focused tests integrate the proposal mass and sampling moments
  independently and retain white/colored furnace means. For roughness 1 and
  view cosines 0.1 and 0.01, second-moment quadrature gives approximately 8.8x
  and 18x lower variance per path than cosine-only reflection, for facet
  albedos 0.4 and 1.0; refinement changes the new moments by less than 4.6e-7.
  This measures constant-illumination estimator variance, not wall-clock speed
  or a universal improvement under arbitrary directional illumination.
  This is the analytic approximation of Portsmouth, Kutz and Hill,
  <https://jcgt.org/published/0014/01/06/>; it does not infer specimen
  calibration, subsurface transport, coherent diffraction, polarization
  retention or a thermal emission law.
- `balance_heuristic` / `power_heuristic` — MIS weights; `mis_weight_sum(pf,
  pg)` — the weight-sum audit (nominally `1`). Finite nonnegative active
  densities retain complementary weights even when count products, sums or
  squares overflow, or squares underflow. Ordinary-range arithmetic is
  unchanged. Zero counts disable their unused density, two zero masses return
  zero, and invalid active densities return NaN. G0/G3 tests in `tests/render.rs`
  check count-weighted analytic ratios, subnormal and extreme density scales,
  and a unit-power integral under a common change of measure. These scalar
  weights do not replace the caller's geometric PDF/Jacobian validation.
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
  `SellmeierIor::try_new_in_band` retains an explicit positive finite vacuum-
  wavelength interval for original oscillator coefficients and analytic group
  index. The closed interval excludes every active pole; monotone endpoint
  checks retain ordinary phase indices `n >= 1` with squared-index arithmetic
  headroom `n² <= f64::MAX/64`. Nonrepresentable
  squared-micrometre bounds and queries outside the source interval refuse.
  No table fit or extrapolation is introduced through native fields, interfaces,
  coherent films or thermal-source consumers. Effective glass support still
  intersects original absorption. The legacy constructor retains its visible
  arithmetic and identities; nondefault bands enter material, camera, film and
  thermal optical identities. This does not broaden the CIE observer or confer
  measurement authority on supplied coefficients.
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
  single-scattering GGX reflection/transmission with isotropic or explicitly
  authored anisotropic slopes. Radiance transport
  uses the documented `(eta_i/eta_t)^2` transmission convention.
  Bare reflection rationalizes the original index contrast and obtains
  transmission from normal Poynting flux, so a rounded unit reflectance does
  not extinguish a positive transmitted channel. Scalar and polarized paths
  share these coefficients and the factored Snell normal wave number. Only
  cancellation at a critical ray uses `fs-math` double-double arithmetic;
  ordinary interfaces remain allocation-free binary64 calculations. Focused
  G1/G3 tests retain near-matched-index refraction, critical-angle side,
  grazing transmission, polarized power closure and finite sampled weights.
  Native spectral bulk absorption retains 2–16 declared vacuum-wavelength
  knots covering exactly 380–780 nm, with nonnegative intensity coefficients
  in inverse metres and deterministic linear interpolation. Imaginary-index
  inputs retain their original dimensionless `k` knots, interpolate `k`, then
  evaluate `alpha = 4*pi*k/lambda_m` at the requested vacuum wavelength. They
  do not replace the source law with linear interpolation of converted alpha.
  The weakly absorbing dielectric interface approximation is unchanged.
  No RGB fit, opacity floor, or extrapolation is introduced. Thin-path absorbed
  fractions use `expm1`, and material identities retain every admitted knot.
  Explicit `SpectralAbsorption::try_new_in_band` and Beer-Lambert band
  constructors retain 2–16 original positive vacuum-wavelength knots on their
  actual inclusive support, including infrared data. The legacy constructors
  still require 380–780 nm and retain identical canonical parameters. The
  imaginary-index band constructor retains the same source parameterization;
  canonical parameters and material identities distinguish it from a supplied
  linear-alpha spectrum. Clear and constant
  extinction are explicit achromatic declarations; lifted RGB retains visible
  support and a source table is never extrapolated. `DielectricGlass` exposes
  the checked intersection of its phase-index and absorption bands, refusing
  disjoint sources while retaining a shared endpoint for monochromatic use.
  No extra allocation or per-query spectral resampling is introduced. G1/G3
  regressions cover original infrared knots, inverse-wavelength absorption for
  constant `k`, affine-`k` knot refinement, path composition, intensity versus
  amplitude attenuation, common-band/refusal boundaries, weak absorption and
  exact visible-constructor parameter agreement. These declarations do not
  establish measurement accuracy, a causal complex-index fit or validity of
  weak-absorption geometric optics for strongly absorbing media.
  Bare `DielectricSurface::interface` accepts positive finite vacuum
  wavelengths, including UV and infrared, using the supplied actual adjacent
  phase indices. Smooth/GGX sampling, evaluation and reverse densities retain
  the same native Fresnel and radiance arithmetic at every admitted wavelength.
  Compact Cauchy coatings remain visible-only, including zero-thickness films;
  borrowed original stacks continue to enforce their own source support.
- **Anisotropic rough glass:** `DielectricSurface::try_anisotropic_rough`
  admits two GGX slope scales in `[1e-4, 1]` and a finite nonzero authored
  tangent axis. `AnisotropicDielectricRoughness` retains both scales and the
  normalized axis; its projection at a contact defines alpha-x, and the
  perpendicular tangent defines alpha-y. A parallel axis explicitly refuses.
  Native scalar scenes interpret the axis in geometry-local coordinates and
  reuse the existing rigid instance/animated placement at the physical path
  time. No tangent is inferred from the camera, material name or measured Ra.
  Equal scales canonicalize to the original isotropic law and material identity.
  Distinct scales and their authored axis enter the new material identity tag.

  The original Walter reflection/transmission owner retains Fresnel, coherent
  coating order, TIR, null-interface measures, radiance eta-squared factors
  and refractive Jacobians. The anisotropic distribution and directional Smith
  masking use the same geometric frame as the original Heitz visible-normal
  warp, generalized to its two-axis stretch. Sampling, forward/reverse PDFs,
  scalar camera packets and dispersive continuations, both BDPT transport
  directions, and thermal/photometric/volume direct-source receivers consume
  the same contact-bound law. Rejected facets retain their original probability
  mass without resampling. Preparation/evaluation has bounded constant work
  and no allocation; no performance benchmark or exact microgeometry claim
  is made. Isotropic paths preserve their original arithmetic and identities.

  Eight G0/G1/G3 native optical tests compare independent anisotropic Walter
  equations, reflection/transmission reciprocity, radiance versus power,
  independent hemispherical PDF/power quadrature, rotated/swapped axes,
  absorbing-film survival, extreme slopes, TIR and exact null measures.
  Together with 130 existing optical/spectral controls they pass in a focused
  original-source build using the original vector operations and optical
  kernels. Three additional native scene/BDPT/packet/refusal regressions are
  authored and pass a real-dependency Cargo metadata check with tracer,
  volumes and photometry enabled; scene-test execution remains unverified
  in the implementation environment. Polarized sampling and reciprocal connections consume the same body-bound
  slope law and the original facet Fresnel Mueller response. This remains
  single-scattering GGX: no inter-facet multiple scattering, finite coherence
  averaging, diffuse substrate, inferred roughness spectrum or thermal law.
- **Coherent coatings on transmitting dielectrics:** `DielectricSurface` owns
  an optional source-bound `DielectricCoating`, separate from homogeneous glass
  and bulk Beer-Lambert absorption. It admits the existing one-to-four-law
  Cauchy unit cell and at most 64 expanded physical layers, with passive
  wavelength-independent extinction for each compact law. The borrowed
  `DielectricCoatingStack` additionally evaluates arbitrary ordered full Cauchy,
  Sellmeier or tabulated transparent laws; `try_mixed` borrows those together
  with original nine-knot absorbing-film laws without fitting, resampling or
  heap-backed scene materials. Authored order runs outside-to-inside and reverses on exits and
  reverse transmission proposals. Indices always use vacuum-wavelength phase
  dispersion. Storage and work stay fixed and bounded, with no per-hit allocation.
  Both dielectric and conductor films share the same passive layer propagation
  and analytic critical-angle limit. The transmitting branch composes incident
  and reflected wave coordinates with factored original index contrasts, so
  weak reflection does not subtract nearly equal total fields. Opaque reflection
  retains its field/flux specialization. Logarithmic normalization scales keep
  evanescent barriers and strongly absorbing films on decaying exponentials;
  transmitted complex amplitudes retain their phase.
  Separate s/p transmitted amplitudes are normalized by outgoing/incident normal
  Poynting flux before Mueller construction. Lossless layers check and normalize
  bounded roundoff in `R+T=1`; absorbing layers retain their physical `R+T<1`.
  Only an above-one excess within the stated roundoff allowance is normalized;
  larger gain refuses. `FresnelDielectric` exposes explicit transmittance,
  preserving weak Bragg leakage when reflectance rounds to one. Its signed
  `absorption_residual()` is a balance diagnostic, not a positive layer-resolved
  heat source or a certificate for sub-ULP losses. Absorbing coatings can remove
  power under external total internal reflection while transmitted power remains
  zero. Existing material identities are retained; bare-interface arithmetic
  follows the shared stable Fresnel law described above.
  `DielectricCoatingStack::layer_power_balance` additionally supplies positive
  total and per-film s/p absorption from that same native transmitting solve.
  It retains factored incident port waves at each boundary and contracts the
  original analytical even/odd dissipative modes; no heat is computed by
  subtracting reflected/transmitted fluxes. Weak film losses survive a rounded
  unit `R+T`, including external TIR with zero transmitted power. The original
  Mueller/Snell result remains unchanged. Reverse incidence reverses optical
  traversal while returned per-layer heat stays in authored outside-to-inside
  order. Lossless/zero-thickness films contribute exactly zero; the real-index
  exit half-space receives transmitted radiation, never film heat. Storage and
  work remain bounded by 64 layers and evaluation allocates nothing. Tests
  compare independent oblique forward/reverse Poynting balances, matched-index
  Beer--Lambert weak loss, opaque limits, exact lossless behavior and native
  source/grazing refusals. This adds no bulk-heat double counting, finite
  coherence averaging, rough-film scattering or inferred optical data.
  Smooth/GGX sampling, hero packet weights and forced splits, forward/reverse
  BSDF PDFs, bidirectional transport, slab/manifold connectors and direct
  photometric/thermal/volume-light receivers consume the same coating response.
  Snell directions, geometric Jacobians, masking and radiance eta factors retain
  their original owners. Polarized beam traversal uses the same reflected and
  transmitted Mueller responses and outward-normal layer orientation.
  A smooth equal-bulk-index coating retains its etalon reflection/transmission;
  empty, zero-thickness and index-matched coatings retain null limits. A
  nontrivial rough coating between exactly equal bulk indices explicitly refuses
  the unsupported mixed continuous-reflection/delta-transmission measure.
  Exact grazing of a nonzero coating also refuses zero incoming normal power.
  Film boundaries are isotropic, locally parallel and laterally uniform; the
  caller declares sufficient coherence across the entire stack. This does not
  introduce scene-path coherence, layer-resolved heating, rough inter-facet
  multiple scattering or thickness statistics. G0/G1/G3 regressions cover
  independent complex Airy R/T phases, quarter-wave antireflection, etalons,
  internal critical/evanescent and terminal-TIR limits, reciprocal flux and
  polarized power closure, deep Bragg transmission, native sampling/PDF replay,
  packet/connector consumption and polarized beam transport.
  Original arbitrary spectral stacks additionally expose `interface` through
  that same scalar smooth/GGX sampling and polarized response owner.
  `DielectricStackBinding` and the straight/curved scenes'
  `with_dielectric_coating_stacks` attach borrowed original layers to exact
  existing uncoated dielectric object IDs. Bounded preparation stores only
  per-boundary slice descriptors; duplicate, unknown, nondielectric and
  conflicting compact/single-film attachments refuse. The actual incident and
  target local bulk laws supply indices at every contact and wavelength.
  Native ordinary/molecular camera continuation, smooth differentiated source
  bundles and shared rough point/area/multiple-image receivers consult that
  same table. Reverse transmission proposals reverse both media and layer
  order; their conditional probabilities remain distinct from physical film
  power. The original Snell directions, bundle area, receiver eta factors and
  bulk Beer absorption are unchanged. Empty bindings preserve existing paths;
  empty physical stacks still retain their declared wavelength support.
  This is a collapsed locally parallel film interface, with no geometric film
  travel or group delay, scene-path coherence, inferred optical fit, or change
  to scalar Copy scene-material storage. Independent spectral Airy controls
  exercise asymmetric absorbing stacks in both native scene directions,
  physical bundle powers/reverse densities, polarized GGX evaluation and actual
  calibrated camera/source connections through the original refractive area.
  At a stochastic nonzero coating, an unsplit spectral packet takes the existing
  one-time per-wavelength continuation path before conditioning on an event.
  This retains companion-wavelength reflection/transmission at a hero-lane
  interference null even when all Snell directions agree. Each child is already
  monochromatic and cannot split again; the existing sixteen-suffix work bound
  is unchanged. Forced smooth splits already carry each lane's exact R/T weight.
  Coated Monte Carlo events place the smaller R/T probability at the low end of
  the variate, retaining weak transmitted support even when R rounds to one.
  Smooth weights, GGX PDFs, polarized throughput and importance transport use
  the same normalized branch probabilities. For absorbing layers, these are
  conditional survival proposals `R/(R+T)` and `T/(R+T)`: sampled throughput
  retains physical power divided by its proposal. Deterministic slab/manifold
  TT and TRT connectors carry physical R/T separately from those MIS densities,
  and completed random-walk paths replay the same conditional probabilities.
  A representably black coating has a zero-weight reflection proposal. Uncoated
  event ordering is unchanged except when bare R rounds to one with positive T:
  that weak transmitted channel also takes the low end of the variate. The
  scalar camera model remains unpolarized at
  each interaction; repeated scalar R/T products do not carry Stokes history.
  Native polarized dielectric sampling, reciprocal connections, ray bundles
  and conductor beams apply their coherent s/p response through the diagonal
  Jones identity `R_s R_p = |r_s conj(r_p)|²`. The weaker channel is recovered
  from the retained U/V block instead of subtracting rounded Mueller I/Q
  terms. Thus a near-Brewster reflected p wave survives even when its power
  is smaller than one ulp of the s wave. General or depolarizing Mueller
  application retains its original semantics; this does not infer coherence
  between independent paths or resolve powers below representable range.
  The exterior longitudinal wave is taken from native Snell/TIR classification,
  preventing a separate squared-index rounding from inventing a propagating
  transmitted mode at a native critical ray; internal film roots remain complex.
  Absorbing-film regressions compare s/p Mueller entries against independent
  interface-scattering/Airy recursion over angle, wavelength and extinction,
  check original mixed spectral knots and reversed layer/media reciprocity,
  and exercise equal-index loss, absorbing TIR, smooth/GGX survival powers,
  weak-loss reflection, opaque/zero-thickness limits and connector MIS density.
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
- **Thin dielectric mesh boundaries:** camera and bidirectional continuations,
  rough-interface connections to rectangle/environment, thermal surface,
  calibrated area and volume-emission sources, and slab/manifold replay
  reconstruct mesh origins from the original triangle barycentrics and project
  componentwise position roundoff onto the geometric normal, with outward
  representable-coordinate rounding. Static and animated rigid mesh instances
  include placement arithmetic. A fixed one-micrometre displacement no longer
  skips a nearby exit and extends its material absorption into the exterior.
  Unoffset barycentric vertices determine absorption lengths and path-space
  densities without long-ray cancellation; instanced world/local witnesses
  remain consistent and traversal retains its original ray parameter.
  Numerical origins never shorten the declared optical path. Implicit charts
  use the native departures described below. This is geometric-optics
  transport, not coherent interference across a geometrically thin body or a
  general certified minimum separation. Analytic regressions cover submicron
  absorption, distant camera origins, oblique Snell/Fresnel power, metre
  scaling and rigid placement.

  The same reconstructed mesh witnesses and bounded normal offsets also
  govern opaque reflection: diffuse, GGX, smooth and brushed conductor paths,
  their direct-light visibility, bidirectional vertices, and finite thermal
  or volume-source connections. Static, rigid and animated mesh placement
  retain their original nearest-hit ray parameters. A numerical launch must
  neither skip a nearby reflector/emitter nor jump over an occluder. Tracer
  bit semantics v5 records bounded origins; v6 also aims finite visibility
  rays at the physical emitter point and validates reconstructed endpoints.
  Physical directions and distances still own BSDFs, PDFs and attenuation.
  Focused analytical controls cover both sides of submicrometre mirror/diffuse
  gaps, instancing, translated scenes, metre-scale light gaps at large world
  coordinates, and environment-light occlusion. Calibrated photometric probes
  retain their authored geometric resolution limit.
- **Native chart departures:** straight camera, bidirectional and shadow paths
  retain physical chart hit points and prepare side-specific ray origins from
  the original signed field, native normal-neighborhood certificates and
  prepared native surface clearance. Direct, rigid and shutter-evaluated
  animated charts share this path. Rigid candidates also clear the exact
  inverse-transform arithmetic used by native intersection. Normal-cone bounds
  cover the complete skipped chord; endpoint clearance alone does not establish
  self-separation. Nearby unresolved separation refuses only the hemisphere
  actually used, so an unused side does not suppress terminal emission.
  Charts without a supplied local normal-neighborhood certificate retain
  legacy self-departure resolution, with certified exclusion of other
  primitives; that compatibility path does not certify uniqueness of their
  own surface sheet. Shared meshes validate once at admission, native
  proximity work remains bounded, and cached box pruning is separately bounded
  by the known primitive count. Physical directions, PDFs and optical distances
  remain attached to the original hit. Camera and light continuations, thermal,
  volume-emission and calibrated-area shadows, and slab/manifold source replay
  retain the same prepared origin. Tracer, dielectric and bidirectional semantic
  versions are respectively 8, 19 and 12. Focused regressions specify both sides
  of submicrometre gaps, direct/rigid/animated charts, large rotated placement,
  blocker visibility, unused-side refusals and unsupported-chart compatibility.
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
  Vacuum wavelengths must be positive and finite; the immutable field owns
  their spectral support independently of the visible color observer. Native
  affine and quadratic preparation checks the common phase/absorption band,
  its endpoints and every contained phase knot throughout the admitted spatial
  bounds. Original infrared tables therefore reach scalar, differential-bundle,
  boundary-event, LTE-emission and polarized curved-beam propagation without
  extrapolation or resampling. The same local wavelength drives phase delay,
  group delay and intensity extinction. Uniform propagation retains its
  one-sample fast path and the adaptive integration adds no new allocation.
  Synthetic infrared fixtures check independent homogeneous/affine path laws,
  quadratic focusing, Planck emission through the same absorbing path, Stokes
  attenuation and common-band refusal. This does not broaden visible-only
  analytic materials, color observers or compact coatings.
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
  An admitted local planar patch or regular signed-distance neighborhood
  permits positive-length departure, including its optical lengths and
  attenuation, without ignoring a whole object or skipping its opposite wall.
  Smooth departures query the actual signed field and use a bounded cone of
  surface normals throughout the neighborhood. The same cone strengthens
  arrival transversality so a near-tangent approach cannot become a false
  crossing. Preparation, arrival and signed-stage queries retain cumulative
  query/work accounting; cancellation and budget errors cannot become retries.
  The selected field must support the explicit numerical contact collar.
  Unsupported neighborhoods, grazing ambiguity and exhausted step/query
  budgets refuse. Spatially uniform fields reuse one
  optical sample and take analytic steps limited by geometry. The caller
  retains the same accepted field and geometry snapshots.
  No scalar camera/BDPT curved visibility, scattering, diffraction,
  birefringence or full-wave accuracy is claimed.
  G1 checks include analytic affine-index trajectories, Snell invariants,
  phase/group separation, absorption and tolerance refinement; G3 checks
  cover reverse paths, transported frames and composed segments; G4 checks
  include cancellation during rejected adaptive trials.

- `gradient_index::bundle` differentiates the same physical-arc-length ray
  equation with respect to two orthonormal launch angles. Immutable
  `SecondOrderIndexField` samples supply the analytic symmetric spatial Hessian
  together with the optical sample; no finite-difference Hessian is inferred.
  Native homogeneous and affine fields provide exact zero Hessians. A bounded
  `QuadraticIndexField` retains native spectral dispersion, group index and
  absorption, admitting its entire spatial box with interval arithmetic,
  including interior extrema. Its complete interval shares the native
  real-index and squared-index headroom checks, so supported high-index
  infrared optics retain quadratic focusing and defocusing without a stale
  visible-glass cap. Analytic axial ray-tube checks retain phase/group path
  length, absorption, footprint area and irradiance for these indices.
  The allocation-free augmented Dormand--Prince
  state controls both angular-derivative columns alongside the central ray,
  differentiates direction normalization and retains explicit step budgets
  and cancellation. Uniform fields use a one-sample analytic fast path.
  Endpoint-plane projection accounts for neighboring rays' different arrival
  lengths and returns signed/absolute area per launch solid angle. Point-source
  irradiance divides spectral intensity by this area and applies central
  Beer--Lambert attenuation once; it needs no additional cosine or index-square
  factor. Explicit grazing, footprint-area and conditioning thresholds refuse
  unresolved caustic results. This is a local differential of one supplied ray,
  not a connecting-ray search, visibility test, finite-pixel caustic integral
  or MIS density. Analytic homogeneous obliquity, quadratic focusing/defocusing and
  independently perturbed affine/mixed-Hessian trajectories check the result.
  `BundleContinuation` preserves the original emission ray, wavelength and
  launch-angle basis across smooth legs and explicit interface maps. The same
  solver accepts arbitrary validated endpoint columns; its uniform fast path
  advances `X` by `distance*V`. Distance, phase/group optical lengths, absorption
  and checked work counts accumulate. The next field must match the retained
  endpoint sample, and mapped endpoints cannot relocate the ray. Integration
  budgets remain per leg with cumulative counts available to the path owner.
  `bundle::interface` differentiates surface arrival, the geometric normal and
  vector Snell/reflection laws. Both incident and transmitted phase-index
  spatial gradients enter the ratio derivative, including equal local indices
  with different gradients. Grazing incidence, critical transmission and
  forbidden total-internal-reflection transmission refuse. The caller supplies
  the admitted normal Jacobian and actual target sample at the existing hit;
  this fixed-branch map does not find boundaries or supply Fresnel power.
  Independent affine-index rays and tangential optical-momentum identities
  check the gradient terms. Normal Jacobians and local integration tolerances
  do not turn a numerical contact collar into a global derivative certificate.

- `photometry::source::bundle` consumes angular bundles in
  the native calibrated point-source path. `PlanarIrradianceReceiver` binds a
  finite point, a one-sided normal and an explicit numerical endpoint-position
  tolerance. `PhotometricSource::spectral_irradiance_from_bundle` verifies the
  retained source origin, wavelength, launch direction and receiver endpoint,
  then uses the source's actual spectral/angular intensity and the bundle's
  receiving-area Jacobian. This bare-bundle entry point rejects paths with
  unaccounted discrete interface events. The refractive entry point and meter's
  `record_refractive` consume `PolarizedRayBundle`, retaining the true native
  interface transfer through an explicit specular path. `BundleIrradianceMeter` borrows one immutable source
  and receiver, admits exactly one supplied connecting bundle at each of the
  80 native spectral midpoints, and uses compensated sums for W/m² and native
  CIE-weighted lux. Wrong source/receiver/bin bindings, duplicate wavelengths
  and incomplete spectra refuse without changing an admitted reading. The
  receiver tolerance is not a finite aperture or a global numerical error
  certificate. This evaluates one supplied branch at each wavelength; it does
  not search for connections, establish visibility, sum all optical branches
  or average a caustic over a pixel.
  Independent full-spectrum checks recover homogeneous inverse-square
  obliquity/attenuation and dispersive quadratic-field focusing in both radiant
  and photopic units. A full-spectrum finite-thickness ball-lens check also
  verifies both curved interfaces, native dispersion, absorption and true
  Fresnel transmission through the actual calibrated meter.

- `polarization::dielectric::transport::bundle::PolarizedRayBundle` composes
  those geometric continuations with the existing native Stokes-frame and
  Fresnel Mueller owner. It retains unit-intensity polarization and logarithmic
  interface power separately from the integrated bulk optical depth. Subsequent
  oblique interfaces therefore see the actual partial/circular polarization;
  deterministic branch contributions contain neither proposal-probability
  division nor a radiance eta-square multiplier. Native transparent coherent
  coatings reuse their existing response and collapsed-interface model,
  including nonzero coated reflection between equal-index bulk media.
  Physical transfer, Beer absorption and receiving-area conversion are applied
  once, with log-domain recovery for representable attenuated irradiance.
  Completely dark polarized branches retain zero power without invalid Stokes
  normalization. A native geometry contact must match the exact endpoint.
  The caller supplies each segment and interface; no connecting-ray search,
  visibility, rough-interface differential, branch completeness or finite-pixel
  caustic claim is added. Independent thick-lens ABCD, off-axis analytic Snell,
  slab, total-reflection, critical-angle, quarter-wave-film and Brewster-null
  references exercise the complete path.

- `polarization::dielectric::scene::curved` connects the field propagator,
  native instance boundaries, spectral interfaces and calibrated receiving
  meters into one forward ray-power path. `CurvedPolarizedScene` borrows
  homogeneous glass or explicitly identified immutable gradient-index fields.
  Each interface uses the actual local phase indices and existing smooth/GGX
  Stokes sampler. Reflection retains its medium; transmission selects the
  explicitly bound opposite medium. Smooth opaque metal bindings borrow the
  same conductor/coating tables as the native polarized scene and retain the
  local gradient medium's index, reflected phase and attenuated power.
  The returned path accumulates physical
  distance, phase/group optical lengths, absorption depth and numerical contact
  resolution. A private beam continuation consumes each solved segment once.
  Surface-event, total attempted-step and geometry-work budgets are cumulative;
  cancellations and geometry/field refusals do not publish partial paths.
  With `photometry`, `SpectralPowerMeter::record_curved` shares the original 80-bin observer and
  all-launch normalization. Receiver and analyzer readings retain contact
  resolution; truncated surviving paths remain explicitly unresolved.
  Geometry/frame identities are checked through separate numerical contacts,
  without constructing certified straight-ray hits. Fresh launches on
  unresolved surfaces refuse, and limited results are not resume tokens.
  This core walk is numerical forward power transport with a required
  material-field contact collar. The separate `camera` consumer below supplies
  scalar radiance; general curved BDPT/MIS remains outside this walk.
  Analytic slab-to-finite-receiver and weighted reflection/transmission
  regressions check exit position/direction, optical lengths, attenuation,
  received watts and Stokes components through the actual native scene.
  Native unmeshed sphere regressions additionally cover off-axis dispersive
  ball-lens entry/exit, internal total reflection with Stokes phase, an affine
  gradient-index trajectory refracting out of a sphere, and exact-tangent
  refusal before interface sampling. Their receiving geometry and optical
  lengths use independent analytic references.

- `curved::bundle::BundleSurfaceSource` admits emission from a named original
  surface on an explicit geometric side. Every connection trial and molecular
  survival replay queries the unchanged native boundary oracle at the actual
  launch point, validates the source object and adjacent medium, and enters
  its existing departure neighborhood. The source is not removed from later
  visibility and the launch point is not offset. Emission already belongs to
  the selected side, so departure adds no source Fresnel/BSDF factor or surface
  event. Native source contact and source-side branch identity remain in the
  result; source admission consumes the original geometry/query budgets.
  Surface-source connection and image APIs share the original differentiated
  shooting, deduplication and survival owners. Duplicate seeds retain one
  source replay, including exact ratio-tracking zeros that still resolve the
  receiver. Invalid source roles, side/medium mismatches, unresolved contacts,
  unsupported departure neighborhoods and work failures refuse explicitly.

- `curved::bundle::connection::BundleConnectionTarget::try_diffuse_receiver`
  selects an original native planar depolarizing diffuse surface as an incident
  power receiver. The same differentiated scene walk and bounded local Newton
  shooting retain actual finite geometry, source and medium identities, optical
  branch plans, unsnapped endpoint resolution and the receiving-area Jacobian.
  The selected contact stops before reflectance or BSDF sampling. A dedicated
  result marker and branch-signature tag distinguish it from an opaque blocker;
  the declared arriving medium and complete prescribed branch plan are checked.
  Other diffuse or rough encounters remain unsupported. Missing, blocked,
  caustic and unresolved connections refuse; this local solver does not discover
  every optical image. Existing detector targets retain their original result.
  Incident point-source irradiance includes native Fresnel power, Beer loss and
  the Jacobian once, independently of the receiver's reflectance. Four native
  regression groups cover detector parity, an off-axis absorbing ball lens,
  branch plans, visibility, finite support, medium identity, contact resolution,
  planarity, work bounds and cancellation.

- `curved::bundle::connection::BundleConnectionTarget::try_volume_receiver`
  connects a source to a point inside a named original optical medium. Its
  numerical target section is composed with the original boundary oracle;
  nearer geometry, native departure neighborhoods, cancellation and work limits
  remain authoritative. The endpoint creates no scene object, material, BSDF or
  physical surface event. A distinct result marker and branch tag retain its
  original medium, and Newton/replay must resolve the actual unsnapped point.
  The section normal controls endpoint solving only: molecular radiometry uses
  the area perpendicular to the actual arriving ray, with no surface cosine.
  Source transmittance replays this same volume target, including complete
  endpoint validation after a legitimate zero ratio sample. Detector and
  diffuse targets retain their existing semantics. This remains a bounded
  declared-branch solver, not a certificate that every optical image was found.

- With `volumes`, `curved::participating::transmittance` connects a selected
  native source branch with its molecular scattering survival. It first solves
  the original deterministic bundle connection, then replays the solved launch
  once through the same differentiated geometry walk. The replay must retain
  the original ordered branch and resolve the same actual receiver. It supplies
  only scattering survival; the original solved bundle remains authoritative
  for Fresnel power, true absorption, focusing and the receiving-area Jacobian.
  Uniform scattering uses exact exponential survival over actual physical arc
  lengths without random draws. Bounded heterogeneous density and spectral
  molecular fields use Poisson ratio tracking on that native curved arc, with
  logarithmic products of `1 - beta_s/majorant`. Scattering is never folded into
  the absorption field or evaluated on an invented straight shadow segment.
  Original source/snapshot identities and wavelength-specific rate/phase
  admission remain binding. A zero ratio factor is a legitimate stochastic zero;
  floating underflow of a positive weight does not establish exact darkness.
  Candidate/contact ordering inside unresolved geometry collars refuses.
  Independent source samples, draw/candidate bounds and native cumulative
  integration/query limits remain explicit; connection trace budgets include
  the final replay. Reported work includes both solving and replay. The estimator
  does not remove native numerical trajectory error, establish branch/image
  completeness or convert budget failures into dark samples. Six native tests
  cover homogeneous survival, prescribed heterogeneous factors on curved paths,
  stochastic zeros, multi-medium refraction, invalid laws and work/cancellation.

- `curved::bundle::connection::images::BundleImageWorkspace` prepares bounded
  image and seed storage for repeated native discovery. The allocating and
  workspace entry points share the same shooting, branch-history, angular
  resolution, parity and ambiguity rules. Failed discovery invalidates the
  workspace result; capacity failures never trigger hidden growth while tracing.
  With `volumes`, `transmittance::images::CurvedTransmittanceImageWorkspace`
  discovers and deduplicates the complete authored seed list before source
  tracking. Each retained image then replays its solved launch once through the
  same original molecular scene, without another Newton solve. The original
  bundle remains authoritative for geometry, Fresnel power, absorption and
  focusing; only independent scattering survival is added. Duplicate seeds
  consume geometric search work but no duplicate power or tracking draws.
  Search and replay share one total trace limit; source draw/candidate limits
  are cumulative across all retained images. A failed replay publishes no
  partial survival catalog. Finite seeds, recorded local misses and numerical
  caustic refusals do not establish image completeness or exact darkness.

- `BundleImageCatalog::validate_proposal_support` admits a deterministic
  finite light proposal paired with the original continuous camera estimator.
  Completed source-side exclusions, proven escape, opaque/detector obstruction,
  another receiver and the wrong receiving hemisphere may leave no represented
  light image. The complementary camera path must use this same catalog and
  retain full emitter weight for an unrepresented image. Empty proposal support
  does not certify zero physical illumination. Missing branch decisions,
  numerical contact/refinement failures, singular/grazing images, exhausted
  work and cancellation remain errors. Strict standalone irradiance keeps its
  all-seeds-resolved requirement because it has no complementary camera term.
  `CurvedBundleSceneError::SourceSideExcluded` is produced only after the
  original source geometry, emitting role and selected-side medium are admitted;
  it retains the actual completed launch/query work for catalog accounting.
  Three authored native regressions cover source-side/medium distinctions,
  geometric obstruction/escape, work accounting and numerical refusal. Rust
  execution remains unavailable without a toolchain/build runner.

- `polarization::dielectric::scene::curved::camera::CurvedCameraScene` consumes
  native curved visibility, full smooth/GGX reciprocal Fresnel/frame transport and
  Beer integrals as scalar unpolarized camera radiance. It retains the native
  branch importance weights and applies one actual local phase-index factor
  `(n_camera/n_source)^2` separately for each encountered source. No source-distance
  or bundle-area factor is added. Validated finite spectral radiance binds
  original native object IDs with explicit emitting sides. Emission is added
  from the arriving beam before that surface reflects, transmits or terminates,
  so a calibrated source can also reflect its surroundings. Each source keeps
  its own Beer prefix and incident local index; the current surface's Fresnel,
  cosine and sampling PDF do not multiply its own emission. Original geometric
  normals determine emitting sides, including brushed surfaces. The private
  callback shares the original walk and uses stack-only accumulation; public
  forward-power APIs retain their original behavior. Source spectra remain
  independently calibrated, with no inferred emissivity or thermal coupling.
  `RadianceSpectrum::try_new_in_band` retains 2--1024 ordered original source
  knots over their inclusive positive vacuum-wavelength band; interpolation
  refuses extrapolation. Native monochromatic camera paths can therefore use
  infrared surface emission with original dispersive media and absorption.
  Legacy visible constructors retain their exact numerical behavior. RGB
  environment reconstruction and CIE raster output remain visible-only;
  extending an optical material does not create missing emitted source data.
  Environment maps require caller radiometric
  calibration and are sampled only after proven homogeneous escape. Any
  positive absorption on the remaining infinite tail extinguishes external
  environment radiance before source evaluation. Native rough dielectric/conductor and anisotropic brushing
  retain their original directional proposal, Fresnel, masking and
  single-scattering approximation. Authored brushed-normal tilt now uses the
  scalar camera law `F_mueller*G2/G1(wo)`, selected before forming the forward-power
  adjoint's correction `(ns.wo/ng.wo)*(ng.wi/ns.wi)`. Both laws keep the same
  body-local frame, mark, spectral/host optics, PDF and horizon support. Public
  forward-beam and molecular-power APIs retain their original importance law;
  surface, raster and molecular camera consumers select the radiance law in
  the shared surface driver. The camera's retained reciprocal beam therefore
  carries an estimator weight, not a forward-power/heat-deposition measurement.
  Decorative normals remain a shading approximation. Rejected
  microfacets have explicit exact-zero continuation without resampling; any
  earlier surface emission remains in the radiance sample. Geometric
  cosine and refractive Jacobians are not applied twice. Detector interactions
  refuse. Numerical limits and indistinguishable
  native extinction/underflow are unresolved
  errors. This covers reciprocal isotropic stationary media and passive
  intensity observations, including the explicit source/analyzer Stokes laws
  described below. Terminal transport weight describes only the final continuation,
  while returned radiance sums all encountered sources. Unresolved continuation,
  overflow and cancellation still return errors without publishing a partial sum.
  `camera::raster::CurvedCameraRaster` prepares the original physical
  pinhole/thin-lens camera and traces each of four stratified wavelengths along
  its own curved path. Raster, lens, wavelength and boundary streams retain
  absolute pixel/sample identities across scheduling and batches. The returned
  sample retains SI spectral radiances, native paths and renderer-normalized
  CIE XYZ, using fixed stack storage and once-prepared color normalization.
  XYZ represents pixel/aperture-averaged scene radiance, not sensor irradiance
  or exposure. Diffraction refuses; every wavelength retains its native path
  limits and cancellation. Twelve focused regressions cover curved visibility,
  phase-index radiance, Beer, branch averaging, emitting sides, calibrated
  directional backgrounds, reciprocal noncoplanar interfaces, independent
  dispersive paths, physical camera sampling and refusal boundaries. Four
  additional rough-surface regressions compare oblique bidirectional glass
  samples to native scalar radiance weights, bare/coated metal to Fresnel and
  Smith masking, body-locked anisotropic brushing under rigid motion, and
  exact-zero rejection without an extra random draw.
  Focused extensions compare the tilted camera weight directly with the native
  scalar BSDF and the forward importance adjoint, retain coated-metal/absorbing-
  film/GRIN references, and exercise spectral raster and molecular camera paths.
  These added native integration regressions remain unexecuted in the
  implementation environment; no multiple-scattering or geometry-detail claim
  is added by admitting the existing bounded normal perturbation.
  Six additive-emission regressions cover analytic diffuse/coated-metal and
  refractive/GRIN sums, exact-black termination, sidedness/null proposals,
  spectral raster and all molecular camera routes including a real collision,
  preserved public power results, limits, overflow and late cancellation.
  They are authored native integration regressions; dependency compilation
  constraints still prevent execution in the implementation environment.

- `transport::diffuse::DiffuseBeamSample` and
  `Interaction::DepolarizingDiffuse { medium, reflectance }` add explicit
  two-sided ideal depolarizing Lambertian reflection to the original straight
  and curved scene walkers. Both `try_new_with_materials` constructors borrow
  one immutable `SpectralReflectance` table; bindings retain indices rather
  than copying spectra into every surface record. Geometric-normal cosine
  sampling uses deterministic math and `M00=rho/pi`, yielding exactly `rho*I`
  with zero Q/U/V. Position, wavelength and the exact local surrounding medium
  remain unchanged. Original departure witnesses, absorption, phase/group
  optical lengths and cumulative work/surface/distance limits continue across
  the bounce. Exact black terminates without a random draw; positive-power
  underflow refuses. The existing scalar camera, raster and molecular-volume
  wrappers consume this same native walk and its source-local phase-index
  conversion. Calibrated camera emitters can also bind diffuse surfaces, adding
  emission before the diffuse interaction, including exact black. Diffuse
  reflection does not infer thermal emission, subsurface transport, material
  polarization, direct-light sampling or MIS. Differentiated
  ray bundles explicitly refuse continuation through the nondelta interaction. Nine regressions
  cover hemisphere power and moments, pole-safe frames, spectral diffuse
  forward detectors in both walkers, gray/white/black camera radiance, analytic
  curved finite-emitter visibility with Beer and group delay, table/medium
  admission, underflow, cumulative limits and bundle refusal.

- With `photometry`, `camera::point_lights::CurvedPointLightCameraScene`
  borrows the exact native camera and at most 64 explicitly bound original
  point-fixture/continuous-receiver branches. Diffuse, rough dielectric and
  rough/brushed conductor receivers preserve their original emission and
  continuation. The original reciprocal polarized surface response contracts
  the arriving source Stokes state with the complete camera analyser; diffuse
  response retains `reflectance * E/pi`. The camera prefix and
  `(n_camera/n_contact)^2` are applied once. Native differentiated source
  connection and calibrated angular spectrum supply incident irradiance E;
  source solid-angle/receiver-area Jacobian, Fresnel transfer and Beer absorption
  occur once. Optional `with_source_polarizations` borrows at most 64 original
  `PointLightSourcePolarization` bindings. Each binding requires pointer identity
  between its scalar fixture and `PolarizedSource::intensity_source()`; a clone
  or differently calibrated source refuses, as do duplicate source pointers.
  The retained source-optics chain supplies normalized Stokes at the exact
  piecewise-constant 80-bin wavelength selected by the original intensity
  spectrum, including 380 and 780 nm. No fraction interpolation or restored
  polarizer loss is introduced. Dark bins and the original source's zero
  out-of-band emission keep zero intensity and use an unpolarized unit
  bookkeeping state; native geometry/visibility checks still apply. The original branch/image
  `source_transverse_axis` declares the polarization frame and must describe
  the same physical source across its declarations. Unlisted sources preserve
  the original unpolarized path. Bindings may precede later molecular or image
  attachments; all three launch paths share the same source lookup. This is
  incoherent Stokes emission, not mutual phase coherence between fixtures or
  separate optical images. Rough glass admits reflection and transmission from either actual
  adjacent medium; its native radiance BSDF already owns local transmission eta.
  Irradiance projects onto the actual source-facing geometric receiving plane.
  Smooth delta receivers remain outside this continuous proposal. Source and
  both adjacent receiver media must support second-order bundles.
  One prescribed branch and seed per original source/receiver pair selects a
  local contribution, without claiming all-image discovery or area-light MIS.
  Duplicate bindings, unsupported media, visibility/numerical failures and
  exhausted cumulative connection/trace limits refuse. Exact black reflectance
  and globally dark sources skip solves; positive underflow is not black.
  The stack-only accumulator consumes no camera random draws, restarts no paths
  and preserves original camera limits. Source-leg traces retain their native
  limits and expose cumulative field/query counts. Participating-media source
  legs remain outside this consumer. `CurvedCameraRaster::sample_pixel_point_lights`
  independently connects every native wavelength lane, preserves physical-camera
  and Philox streams, and feeds combined spectral radiance to original XYZ
  integration. It retains each lane's fixture contribution and work counts;
  any lane's failed source solve refuses the complete pixel with its original
  error. Native regressions cover calibrated cosine/Beer radiometry, quadratic
  GRIN focusing, affine-index endpoint factors, multiple diffuse bounces and
  emission, source/receiver admission, numerical darkness, cumulative budgets,
  cancellation, empty-fixture pixel parity and spectral pixel radiometry.
  Continuous-receiver controls cover bidirectional rough-glass reflection and
  transmission with radiance reciprocity, complex-index and anisotropic metal
  response, non-coplanar camera polarization, and a calibrated narrow-beam source
  through a glass slab with an independent refracted bundle-area reference.
  Source-polarization regressions compare rough Brewster reflection and
  conditional Rayleigh lighting against independent analytic formulas,
  including discovered-image duplicate seeds and polarizer throughput once.

- With `photometry` and `volumes`,
  `camera::point_lights::participating::CurvedRayleighPointLightCameraScene`
  composes the exact original molecular camera with its admitted point fixtures.
  At matching continuous surface arrivals, source scattering survival multiplies
  the existing calibrated polarized source contribution once. Camera-prefix
  survival remains exclusively in the native molecular free-flight law; there
  is no second camera extinction factor, phase weighting or collision rate in
  the surface connection. Rough-receiver controls compare the same original
  molecular camera walk with independent source-only scattering attenuation. Original surface emission, continuation and all camera random
  streams remain unchanged. Source tracking has independent counter streams
  and cumulative draw/candidate limits across all connections and bounces.
  Legitimate zero ratio samples retain camera emission and continuation without
  retrying a source path. Source failures preserve their native errors.
  `CurvedCameraRaster::sample_pixel_rayleigh_point_lights` carries this behavior
  through each original wavelength lane and native XYZ integration, retaining
  separate camera and source counts. Failed lanes refuse the complete pixel.
  Optional `MolecularPointLightBranch` bindings also supply direct point-fixture
  contributions at accepted molecular collisions in their original receiver
  medium. Null candidates contribute nothing. The observer reads the incoming
  camera path and local admitted phase before continuation, so its contribution
  includes neither the continuation-phase PDF nor a second collision-rate or
  camera-survival factor. Each native source connection retains its original
  Mueller/Fresnel transfer, absorption, transverse bundle area and independently
  tracked scattering survival. The normalized Rayleigh contraction uses all
  Stokes components, shared transverse frames and reciprocal direction reversal;
  multiplying source and camera intensities alone is insufficient. The local
  `(n_camera/n_vertex)^2` factor is applied once. Reciprocal scalar camera
  response is supported; arbitrary nonreciprocal optics and polarization-resolved
  sensor channels are not inferred. Surface and volume connections share the
  same cumulative work budgets and source ordinals without changing the camera
  random streams. Existing constructors keep empty volume bindings. Native
  single-scattering pixel tests compare independent spectral radiometry and XYZ
  while requiring identical original camera paths. General area-light MIS and
  complete optical-image discovery remain outside this consumer.

- With `photometry` and `volumes`,
  `camera::point_lights::participating::images::CurvedRayleighImageLightCameraScene`
  admits finite image groups for original continuous surfaces or molecular
  receiver media. `PointLightImageReceiver::SurfaceObject` admits diffuse,
  rough dielectric and rough/brushed conductor laws; the original `DiffuseObject`
  spelling remains equivalent, including duplicate source/receiver identity.
  At each matching original camera vertex it discovers distinct native
  optical images, tracks their source survival separately and adds their
  incoherent radiometric contributions. Each image keeps its own calibrated
  angular intensity, Mueller/Fresnel transfer, absorption and bundle focusing.
  Each image uses its original reciprocal polarized surface or molecular
  contraction. Rough-glass images may arrive from either adjacent medium;
  source-facing irradiance projection and the native radiance transmission law
  retain their original index factors exactly once. Camera continuation and surface/flight
  random streams remain unchanged. Source ordinals and work budgets are shared
  with any separately bound single connections; duplicate original
  source/receiver bindings refuse. Empty catalogs, recorded local misses,
  ambiguity and exhausted budgets refuse the sample instead of silently
  dropping unresolved light. A caller-owned workspace is reused across
  vertices and wavelengths without a mutable scene cache.
  `CurvedCameraRaster::sample_pixel_rayleigh_image_lights` carries these
  contributions through all four native wavelength lanes and original XYZ
  integration, preserving complete-pixel failure and independent source
  counters. Finite authored seeds do not prove all optical images were found;
  coherent interference, singular caustic regularization and area-light MIS
  remain outside this consumer. G1 controls sum independent TT/TRRT slab
  throughput, physical absorption/scattering lengths and unfolded bundle areas
  at rough reflecting/transmitting glass and rough/brushed metal receivers.
  G3 checks duplicate-seed suppression, compatibility-name identity, exact trace
  budgets and dispersive image sums in all native spectral pixel lanes.

- Native polarized bundles retain `reverse_branch_log_probability`
  separately from physical Mueller power and Beer absorption. This is the
  exact scalar proposal for the selected smooth chain traversed backward;
  transmitted absorbing multilayers use the opposite local indices and
  reversed authored layer order. Smooth conductors contribute probability
  one. Camera surface observers can retain the original sampled direction,
  event and discrete/solid-angle proposal without resampling the native walk.
- `camera::area_lights::mis` evaluates complementary power-heuristic weights
  for a regular source-to-diffuse image in one receiver solid-angle measure.
  The source area density uses the reciprocal native footprint and actual
  endpoint phase indices; the competitor is the native cosine PDF times the
  reverse discrete branch proposals. Per-image enumeration adds no image-count
  divisor. Original-face hit lookup retains the represented area density and
  unsnapped camera point, refusing edge or owner ambiguity. These local weights
  do not establish membership of a separate camera path in a finite seed
  catalog. G1/G3 cases cover oblique receiving area, common refractive host,
  scalar versus polarized lens power and an asymmetric absorbing coating.

- With `volumes`, the native molecular observer also retains successful
  surface and molecular departures, including the original incoming beam,
  sampled outgoing beam and scalar phase or surface proposal. Null tracking
  candidates produce no physical vertex and do not reset a connection suffix.
  Existing arrival-only consumers share the same driver and original random
  draws, optical arithmetic and path limits.
  `CurvedAreaEmitter::volume_mis_weights` uses the actual ray-transverse bundle
  footprint, original endpoint phase indices and reverse smooth-interface
  probabilities. An arbitrary numerical target section introduces no physical
  receiver cosine. Participating consumers use complementary deterministic
  angular/branch weights; these omit no-collision survival and do not claim
  full marginalized path densities or variance-optimal volumetric weights.
  Physical source survival remains independent of the nonlinear heuristic;
  a noisy ratio-tracking estimate is never substituted into its denominator.
  G1/G3 regressions cover oblique target-section invariance, common-index and
  absorption independence, native refractive proposals, and unchanged density
  and local-field paths through null and accepted collisions. These native
  regressions are authored but unexecuted in the implementation environment.

- `camera::area_lights::CurvedAreaEmitter` prepares area sampling on original
  native geometry and its original calibrated or thermal emission law,
  directly from the surface scene or an admitted camera emitter. Mesh sources may include all
  original faces or an explicit strictly increasing face selection through
  `try_new_faces` / `from_camera_faces`; disconnected, noncoplanar and closed
  meshes retain original face winding and per-side radiance. `AreaEmitterLimits`
  bounds the whole mesh, preparation work and each hit-membership query.
  Admission reuses the original BVH and exact coplanarity/orientation
  predicates to reject positive-area coplanar overlap, including overlap with
  excluded faces. Shared edges and vertices remain zero-area contacts; no
  mesh-union or multiplicity correction is inferred. The triangle-area CDF
  retains original geometry, frame and the exact emission law. Samples
  retain the exact prepared support, original face index and represented
  conditional area PDF; proper-rigid placement preserves area. Native camera
  hit lookup reuses the original BVH, checks actual three-dimensional plane
  membership and strict triangle interior, and keeps the unsnapped camera
  point. Edges, multiple memberships, numerical range and exhausted bounds
  refuse. Sampling and hit lookup introduce no per-call allocation, chart
  conversion or surrogate point fixture. Image seeds remain original world
  directions. Native surface-source discovery supplies each
  distinct regular image. Its sampled differential intensity is
  `Le * emitting_cosine / area_pdf`; the original polarized bundle supplies
  Fresnel power, absorption, focusing and receiver obliquity exactly once.
  Volume targets use the actual arriving ray's transverse area. Optional
  molecular survival applies independently to each image on the same scene.
  Aggregate work includes all authored seed searches and retained replays.
  This standalone irradiance consumer still requires every authored seed to
  resolve a retained image; an excluded seed or empty catalog refuses. Finite
  seeds do not establish image completeness. It does not add a second
  emission term to camera continuation. Singular caustic integration and
  coherent interference remain outside it.

- Thermal area proposals retain the complete `DiffuseThermalEmitter` snapshot.
  The thermal constructors and original `from_camera` constructors reuse native
  geometry sampling and image discovery. NEE and complementary emitter hits
  evaluate `epsilon(lambda)*B_lambda(T)` at the actual wavelength, without a
  tabulated surrogate, extra pi/side factor or second heat step. Exact Kirchhoff
  material agreement, source-band support and local index-one emission remain
  required; camera matching includes temperature, emissivity and source identity.
  `refresh_thermal_emitter` validates and atomically replaces only that snapshot,
  preserving prepared mesh CDF/chart metadata and work limits. Its mutable
  borrow excludes stale samples, and camera bindings must match the refreshed
  state. `source_emitter()` returns an optional calibrated owner;
  `source_thermal_emitter()` exposes the original thermal state. Native surface
  and participating MIS share this source path. Focused native regression
  fixtures cover Planck radiometry, complementary hits, thermal QE response,
  geometry reuse and refusal; full renderer execution remains unavailable in
  this implementation environment.

- Native chart area emitters use the original owner's explicit
  `Chart::{surface_area, sample_surface_area, surface_area_pdf_at}` capability.
  `SphereChart` provides closed-form uniform area sampling and matching hit
  density without tessellation. `CurvedAreaEmitter::sample_area_with_context`
  samples either backend; legacy `sample_area` retains mesh behavior and
  returns `ContextRequired` for a chart. `AreaEmitterSample::triangle_index`
  now returns `Option<usize>`: an original mesh face or `None` for a chart.
  Samples retain the exact prepared support identity, native outward normal
  and area PDF under proper-rigid placement. Chart preparation and hit lookup
  each require one bounded native operation; zero work limits refuse. Native
  face-subset constructors reject charts, and overlapping whole-chart bindings
  remain duplicate support. Camera-hit density queries retain the unsnapped
  integrated point and validate the original normal and object/frame identity.
  Surface and molecular cameras invoke the same contextual sampler and
  reciprocal MIS path. Each emitting side uses its original adjacent medium
  and calibrated per-side radiance, with no source-interface Fresnel event.
  Six authored renderer regressions cover native sphere measure/placement,
  source ownership and work bounds, actual refractive surface/molecular camera
  paths, matching hit PDFs and both sides of an emitting dielectric sphere.
  Independent numerical Snell/radiometry controls passed; Rust execution
  remains unavailable without Rust/DSR/RCH.

- `BundleConnectionTarget::try_surface_receiver` stops at the original
  continuous scattering material before its receiving BSDF or local branch
  choice. Rough glass, rough/brushed conductors and diffuse surfaces retain
  native geometry, including curved charts. The actual geometric approach side
  selects the glass's original incident medium; equal IOR never substitutes for
  that identity. The target normal defines numerical shooting coordinates only.
  `BundleConnection::receiver_projection`, image parity/resolution and source
  irradiance use the actual endpoint tangent facing that image's arrival.
  Full unsnapped three-dimensional residuals, finite support, preceding smooth
  refraction, opacity and complete branch plans remain authoritative. No
  derivative through a rough scattering event is inferred.
  `CurvedAreaEmitter::surface_mis_weights` converts the original area PDF into
  receiver solid angle using native etendue and the actual light-side phase
  index. Its competing PDF is the native continuous BSDF proposal, including
  local reflection/transmission selection, times the original reverse smooth
  branch probability. Physical Mueller power and absorption never substitute
  for proposal probability. Five authored regressions cover curved receiving
  derivatives/area on both sides, material/media admission, lens Fresnel/Beer
  transport and complementary directional weights. Rust execution remains
  unavailable in the implementation environment.

- `transport::connection::SurfaceConnectionResponse` evaluates the original
  continuous dielectric, rough-metal or body-locked brushed-metal response at
  a fixed source direction. It retains the scalar native proposal separately
  from the full reciprocal Stokes contraction, including Q/U/V and coherent
  coating phase. A dielectric connection reuses one original coherent solve
  for its BSDF and proposal; no direction or branch is resampled. The other
  original dielectric medium determines Fresnel contrast even on reflection.
  The radiance response already includes the local transmission eta factor.
  Consumers multiply camera-prefix intensity, actual receiving irradiance and
  `(n_camera/n_camera_arrival)^2` exactly once. Light-side etendue belongs to
  the proposal conversion, not another physical index factor. The brushed
  response uses the native camera radiance law and geometric receiving cosine,
  retaining instance-local detail, marks and the original pole/zero fallback.
  A positive coherency overlap preserves crossed-polarization extinction with
  a positive scalar PDF. Smooth/null delta transport remains outside continuous
  densities; coated equal-index rough mixtures retain native refusal. These
  evaluations allocate nothing and consume no random draws. Nine authored
  regressions compare native samples, scalar laws, both glass sides, coherent
  coatings, polarization/frame changes and brushed radiance/adjoint behavior.
  Independent numerical coherency and index-composition checks passed; Rust
  regressions remain unexecuted without a toolchain/build runner.

- `camera::area_lights::camera::CurvedAreaLightCameraScene` composes
  original mesh-face area samples with the native surface camera at declared
  continuous receivers: diffuse, GGX glass, rough metal and body-locked brushed
  metal. Each regular refractive image receives its own power-heuristic weight
  before physical irradiance, the original polarized radiance BSDF and the
  camera-prefix index ratio are composed. The camera retains its original
  continuation draws and native path. At an emitter hit, only the retained
  last-continuous-surface suffix is reversed: original ordered object/geometry/frame/
  feature, media and selected branches must agree, as must resolved directions
  and actual scalar continuation proposals. Rough anchors retain the original
  incoming polarized beam, sampled outgoing direction and post-scatter medium.
  Their sharp GGX proposal is checked at that exact recorded direction, while
  native image resolution separately validates reversed geometry. Each matching
  image evaluates its own directional PDF. The same authored light-seed
  catalog is reevaluated at that source point without adding the known camera
  direction as a light seed. Native angular-resolution, separation and parity
  rules decide overlap; an unrepresented image keeps full emitter-hit weight.
  Resolved geometric exclusions use `validate_proposal_support` in both
  next-event and emitter-hit evaluation; an empty light catalog contributes no
  next-event term while the native continuation keeps its full emission.
  Disjoint face selections may share an emitter/receiver/side binding; each
  retains its own conditional area PDF and authored seeds and is sampled once,
  without a selection-count divisor. Overlapping selections refuse. Emitter
  hits select the containing original face group, independently of binding
  order; uncovered faces retain full native emission. Ambiguous or numerical
  failures still refuse the complete sample. Direct
  emitter visibility and unsupported nondelta suffixes keep original emission.
  One caller-owned bounded workspace retains source images and up to 64 delta
  events; at most 64 source/receiver/side bindings are admitted. Cumulative
  area-sample and trace limits are charged before new work, including reverse
  reconstruction and overlap discovery. Sampling and tracing add no per-ray
  allocation. This consumer has no molecular/free-flight proposal binding.
- `camera::raster::area_lights` integrates that surface estimator into the
  original physical-camera pixel and XYZ quadrature. All four wavelengths
  independently trace their camera and source paths; raster, lens, wavelength
  and continuation streams retain their identities. Area variates use a
  separate absolute pixel/sample/event/binding/lane Philox domain. Pixel
  results retain per-lane radiance and actual source work, and a failed lane
  refuses the whole pixel. G0/G1 regressions cover analytic paired estimators,
  unbound/direct sources, an unseeded ball-lens image, source ownership,
  cumulative limits, aperture/sample identity and independent spectral NEE.
  Four additional authored regressions cover sharp rough/anisotropic metal
  paired estimators, noncoplanar polarized prefixes, heterogeneous rough-glass
  source survival and independent dispersive transmission in all four physical
  camera wavelength lanes. Native camera paths and random draws remain
  identical. The Rust regressions are unexecuted without a toolchain/build runner.

- With `volumes`,
  `camera::area_lights::participating::CurvedRayleighAreaLightCameraScene`
  adds native mesh-face next-event sampling at original continuous surfaces and
  real molecular collisions. It borrows the exact admitted camera, face selections
  and molecular scene. Each retained refractive image reuses its differentiated
  native geometry and independently tracked scattering survival. Molecular
  radiometry contracts both subpaths' full normalized Stokes states with the
  local reciprocal phase; source absorption, Fresnel transfer, focusing and
  the camera-to-vertex phase-index ratio occur once. No second camera survival,
  collision rate or sampled continuation-phase divisor is introduced.
  Continuous surfaces reuse the same original receiving BSDF, geometric
  tangent, adjacent-media validation and camera radiance law as the clear
  camera. Surface and molecular paths share original-source validation and
  reverse native history/image-membership implementation. Real collisions
  replace the angular anchor, including in unbound media; null candidates
  preserve it. Full emitter weight remains for unrepresented optical images
  and uncovered source faces. Disjoint selections of one emitter may bind the
  same molecular medium with their own conditional area PDFs and image seeds;
  the same geometric proposal-support rules apply to clear and molecular paths.
  Source tracking and area draws use independent streams and cumulative work
  limits. Caller-owned bounded image/history storage is reused without per-ray
  allocation. The deterministic angular partition deliberately omits survival
  from the heuristic, as documented above; noisy source survival remains a
  separate radiometric estimator.
  `CurvedCameraRaster::sample_pixel_rayleigh_area_lights` integrates all four
  wavelength lanes into the original physical-camera and XYZ sample, retaining
  separate camera and source work. Any failed lane refuses the whole pixel.
  G0/G1 regressions cover complementary Rayleigh/depolarized emitter hits,
  two polarized collisions, null-history retention, independently varied
  source ratios, cumulative area/trace limits, unchanged native pixels and
  independently evaluated spectral source loss. Rust execution is unavailable
  in the implementation environment; these regressions remain unexecuted.
  Additional authored regressions cover nonplanar/closed mesh sources, exact
  support ownership, conditional face densities, disjoint bindings on surface
  and molecular paths, excluded source sides and partial occlusion. Independent
  numerical references and 3,000 exact planar-overlap comparisons passed; Rust
  compilation and regressions remain unexecuted without Rust/DSR/RCH. Charts
  lacking explicit area sampling or native departure authority, singular caustic
  integration and global optical-image
  completeness remain outside scope.

- `camera::IsothermalVolumeEmitter` explicitly binds a frozen temperature state
  to an original curved medium through `with_isothermal_volume_emitters`.
  Positive kelvin temperature, nonzero source identity, valid medium ownership
  and unique bindings are required. The same native true absorption supplies
  Kirchhoff LTE emission; no opacity, thermal emissivity or temperature is
  inferred from unrelated data. The shared Planck owner supplies vacuum
  spectral radiance. Each isothermal leg contributes
  `n_camera² * I_departure * B_lambda(T) * (1-exp(-Delta_tau_abs))`,
  with a stable thin-optical-depth evaluation. The local source phase index
  cancels in reduced radiance; group index and extra Fresnel/geometry factors
  do not enter this source. Spatially varying phase index and true absorption
  follow the original curved integrator without extra field/geometry queries.
  Ordinary, molecular, point/area-light and spectral raster cameras use the
  same accumulator. Arrivals close a leg before scattering, and departures
  retain the native Mueller/PDF throughput. Null candidates do not reset it.
  Scattering survival is already in the flight proposal; surface-source MIS
  weights do not attenuate thermal emission. Planck evaluation is reused
  across successive collisions in the same medium, with stack-only ray state.
  A proven homogeneous escape with positive absorption supplies the finite
  remaining equilibrium source while extinguishing the external environment.
  Unresolved and extinguished paths retain their existing refusals. Empty
  source bindings preserve existing camera results. Eight authored native
  regressions cover finite/transparent/thin slabs, spatial absorption on a
  curved GRIN ray, phase/group normalization, Fresnel transmission, infinite
  thermal equilibrium, source admission, spectral pixels, two polarized
  collisions and null tracking. Spatially varying temperature and thermal
  direct-source sampling remain outside this isothermal capability.

- `camera::SpatialVolumeEmitter` binds a borrowed, immutable `TemperatureField`
  through `with_spatial_volume_emitters`. Its nonzero source identity and medium
  ownership are explicit; overlapping spatial/isothermal laws refuse. The
  original optical Dormand–Prince stages also integrate
  `alpha_abs(s) * B_lambda(T(r(s))) * exp(-tau_abs(s))` in reduced radiance.
  Temperature queries follow the same stage clearance/departure guards. The
  embedded source estimate joins the native optical acceptance decision using
  a positive per-step absolute reduced-radiance tolerance and the original
  relative tolerance. Rejected increments are never published, the accepted
  endpoint derivative is reused, and original optical, geometry, cancellation
  and work budgets remain binding. Uniform optical samples stay prepared.
  Before segment Beer attenuation, the native current Mueller intensity weights
  its source integral; the common camera accumulator applies `n_camera²` once.
  Physical and null collisions, boundaries, source MIS and every existing
  ordinary/molecular/point/area/raster consumer retain those same semantics.
  With `material-state`, `AxialGasTemperatureField` embeds an original
  `fs_material::gas_field::AxialGasField` along an explicit world axis and queries
  its source-resolved `state_at`; the original model and phase domains remain
  authoritative. Cellwise reconstructions and axial extrapolation refuse.
  Linear primitive station spans are continuous and piecewise smooth; callers
  must resolve/refine those spans and all other source structure with the step
  cap. Embedded estimates are not certified global error bounds or guarantees
  across missed features. Unknown absorbing infinite thermal tails refuse;
  globally transparent homogeneous paths skip temperature queries and preserve
  the original analytic path. Empty bindings preserve existing bits, and
  isothermal sources retain their closed-form fast path. Nine authored native
  regressions cover adaptive source refinement, prepared-optics/FSAL reuse,
  transparent/disabled equivalence, straight and GRIN temperature gradients,
  phase/group and interface normalization, isothermal reduction, source-domain,
  cancellation/work/infinite-tail refusals, polarized collisions and null/zero
  scattering. Thermal direct-source sampling and energy feedback remain outside
  these passive camera source bindings.

- With `volumes`, `camera::participating::CurvedRayleighCameraScene`
  composes the prepared scalar camera with the exact original
  `CurvedRayleighScene` surface owner, checked by pointer identity. Native
  arc-length flights, heterogeneous majorant/null tracking, complete
  Rayleigh/depolarized-Rayleigh Mueller/PDF response and absorption are reused
  unchanged. The same source-local radiance conversion applies once per
  encounter; reflective emitters use the shared arriving-surface callback. No
  duplicate scattering survival, collision rate, albedo or phase weight is
  introduced. Samples retain real collision and free-flight counters. Zero
  scattering preserves surface-camera results bit for bit. Positive-majorant
  tails do not become vacuum because one density query is zero. All native
  work/distance/volume/cancellation refusals propagate. Optional isothermal or
  spatial LTE emission is accumulated on these same original absorption paths.
  General particle polarization inference and polarized camera adjoints remain
  outside this consumer.
  `CurvedCameraRaster::sample_pixel_rayleigh` uses the same physical camera,
  spectral stratification and color integration, with separate per-wavelength
  surface/flight streams and a fourth independent delta-acceptance variate.
  Every lane retains its volume counters and path budgets. Nine regressions
  cover analytic uncollided GRIN transport, source visibility after scattering,
  polarization-dependent multiple collisions, null tracking, owner identity,
  zero-scattering pixel equivalence, conditional Beer weights and refusals.

- With `volumes`, `CurvedRayleighScene::try_new_with_fields` additionally binds
  local spectral molecular rates and depolarization to the exact original GRIN
  snapshot identity, optionally alongside unchanged density-scaled laws.
  Immutable wavelength majorants drive the same native curved-arc delta tracking;
  local rates and phases are validated, accepted events reuse native Mueller/PDF
  response, and null events retain beam weight. One stack cache reuses the last
  field medium's sampled majorant without allocating per-path tables. Positive
  field majorants require an independent fourth tracking variate; exact zero
  consumes no flight RNG. Existing camera/raster owners inherit the field law
  with original local source radiometry, geometry, absorption, cancellation and
  cumulative budgets. Local vacuum never proves escape under a positive global
  majorant. Eight regressions cover analytic arcs, changing local phases through
  actual camera radiance, nulls/cache reuse, legacy bits, source identity, invalid
  inputs and native limits. No particle-polarization inference, volume emission,
  extra extinction weight, direct-light estimator or MIS is added.

- `polarization::photometry::SpectralOptics` (feature `photometry`) connects
  native metal and coherent multilayer reflection to calibrated source optics.
  `conductor_reflection` retains the original conductor's attached coatings;
  `conductor_stack_reflection` consumes an arbitrary borrowed ordered spectral
  stack with the native replacement semantics. `coated_dielectric_channel`
  selects original reflected or power-normalized transmitted Mueller response,
  including absorption and explicitly reversed layer order with the actual
  adjacent media. Each constructor evaluates the original phase-index and
  complex-field owners at the existing photometric quadrature wavelengths.
  Invalid angles or unsupported source wavelengths refuse construction.
  The resulting fixed-size optic retains no layer borrow and adds no per-ray
  solve. The existing `PolarizedSource::through` chain preserves Q/U/V and
  absolute surviving power relative to the original source normalization, so
  mirror phase and coating interference affect downstream analyzers, lumens,
  candela and spectral source intensity. Caller-owned s/p frames and outgoing
  placement remain explicit; no radiance eta factor, extra absorption or
  stochastic event probability is applied. Native scalar rendering consumes
  the final intensity source; it does not inherit the chain's Stokes history.
  Source chains retain two positive, independent coherency factors internally,
  applying precomputed native Jones responses without interference between
  populations. This preserves finite-extinction leakage and weak Brewster
  populations through later analyzers even when intermediate binary64 Stokes
  I/Q round to equality. Public Stokes and calibrated power remain the original
  interfaces. Exact authored cardinal axes avoid a trigonometric leakage floor;
  arbitrary near-null field cancellation retains ordinary binary64 limits.
- `conductor` module (feature `tracer`): fixed-capacity visible-band tables of
  absolute complex refractive index `eta(lambda) + i k(lambda)` in vacuum
  nanometres, strict inclusive 380--780 nm support, deterministic linear
  interpolation, exact unpolarized complex-IOR Fresnel against the path's active
  incident medium, and validated isotropic-GGX roughness. Built-in tungsten and
  stainless-steel tables are explicitly representative look-development data;
  table/source/material identities prove replay binding, not specimen fidelity.
  `ConductorSurface::SMOOTH` selects exact delta reflection without a numerical
  roughness floor. Camera packets, bidirectional radiance/importance paths and
  the polarized surface adapter retain the same spectral/coated Fresnel power
  and the actual incident host index. The mirror event has unit discrete
  probability, zero continuous BSDF/PDF and no finite-direction NEE competitor;
  emitter/environment hits keep their full delta-path contribution. Rough
  surfaces and their content identities retain their original semantics.
  Material-state conductor bindings admit explicit zero alpha and preserve it
  across optical-state refresh. This geometric-optics limit supplies no
  diffraction, finite-roughness approximation or mirror-caustic connector.
  Scalar complex square roots compute the larger component first and recover
  the smaller from `2uv=y`, retaining weak absorption on the total-reflection
  side. The shared scalar/polarized bare amplitudes factor index contrast and
  rationalize weak reflected waves. Bare, single-film, multilayer and per-film
  heating evaluations retain the incident normal wave in the longitudinal
  wavevector even when a grazing `cos²` is smaller than machine epsilon.
  Optical semantics version 3 distinguishes this evaluation from cancellation
  that could turn an almost matched transparent interface into a mirror;
  authored tables, coatings and source metadata retain their original values.
  Independent high-precision Fresnel/Airy references cover weak absorption on
  both sides of critical incidence, matched grazing media, reflection phase,
  one-ulp film contrast and coherent reflection/absorption energy closure.
  `ConductorOptics::with_coating(CoherentDielectricCoating)` additionally admits
  one explicitly source-bound transparent isotropic film: the shared visible
  Cauchy, Sellmeier or tabulated phase-index law and a uniform thickness in metres, from
  zero through 100 micrometres. At each wavelength the same phase index sets
  both the interface admittance and the coherent optical thickness; group
  index is never substituted. Cauchy coatings retain their existing identity
  encoding; Sellmeier binds a separate law tag and all six coefficients,
  while tabulated laws bind each wavelength/index pair. Separate complex s/p amplitudes retain substrate reflection
  phase and sum every internal reflection before unpolarized power averaging.
  A fixed-size characteristic-matrix evaluation handles the film's critical
  angle by its analytic sinc limit and scales evanescent transfer entries to
  preserve frustrated total internal reflection without exponential overflow.
  Both isotropic and brushed conductor BSDFs, including their existing direct,
  camera, photometric and bidirectional consumers, use this same Fresnel owner.
  The film follows the local microfacet normal; geometric proposals remain
  unchanged. Bare presets retain their authored optical coefficients. A
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
  `ConductorOptics::with_periodic_coating(CoherentPeriodicCoating)` additionally
  integrates explicitly ordered multilayers into those same scalar BSDF and
  Mueller consumers. One to four `CauchyCoatingLayer` laws form a unit cell,
  repeated coherently up to 64 total physical layers. Each layer declares its
  shared Cauchy phase dispersion, constant visible-band extinction coefficient,
  and normal thickness; the whole coating carries a source/status binding.
  These are authored optical laws, not a silent constant-index or constant-loss
  fit to arbitrary measured data. One repeat admits any ordering of four such
  layers; multiple repeats support Bragg stacks. Cell propagation is prepared
  once per wavelength/angle and reused across repeats. Storage fits inside the
  existing absorbing-film payload, preserving Copy material size and avoiding
  per-hit allocation. Identity binds physical ordering, repetition count, every
  coefficient, thickness and source; adding multilayer support preserves
  existing bare and single-film identity encodings.
  `CoherentCoatingStack` separately borrows up to 64 arbitrary ordered transparent
  or absorbing layers, retaining each full Sellmeier/Cauchy/tabulated law or
  nine-knot eta/k spectrum and source. It evaluates scalar reflectance and the
  signed Mueller response through the same kernel. This general borrowed API
  replaces any coating on its supplied substrate for that evaluation; it does
  not attach arbitrary heap-backed spectral stacks to Copy scene materials.
  Backward boundary-field propagation normalizes each field/flux pair after
  every layer and shares the passive scaled propagation and analytic critical
  limit with single films. No exponentially growing fields, intensity-only
  mixing, effective-medium replacement, or amplitude-based layer pruning is
  used. Zero-thickness layers disappear exactly. Bare, incident-matched prefix
  and substrate-matched suffix limits retain established scalar/Mueller bits.
  The caller declares sufficient coherence over the full stack and locally
  parallel, laterally uniform boundaries. This is local multilayer interference;
  it does not infer thickness distributions, spatial diffraction or scene-wide
  coherent fields, or make the scalar tracer transport polarization between
  bounces. Absorption remains the total opaque-system `1-R`, with no claim to
  layer-resolved heat deposition or broadband thermal emissivity.
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
  Multilayer regressions additionally cover independent nested complex Airy
  and signed Mueller answers, closed-form quarter-wave Bragg admittance in both
  orders, full eta/k spectral interpolation, expanded/repeated-cell equivalence,
  critical and evanescent weak-loss limits, zero/split-layer reductions, opaque
  front-layer screening, bounded storage/work, and native isotropic/brushed BSDFs.
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

  `try_spectral_conductor` and `try_spectral_tabulated_dielectric` additionally
  consume the explicit finite-band source snapshots from `fs-material`. Each
  original wavelength converts once from SI metres to vacuum nm; nine original
  `(eta,k)` or `(n,alpha)` knots reach the native conductor or tabulated dielectric
  owners without visible-band projection. Conductor knots and identity come from
  the existing `SpectralComplexIndexStatePoint` projection exposed by the SI
  source wrapper; no second n/k law is introduced. Dielectric absorption retains
  original linear-alpha interpolation. The first and last source knots are
  the actual closed support. A positive finite caller-declared maximum knot gap
  is checked by those owners; it declares interpolation admission, not accuracy.
  Smooth and GGX finishes remain independent, and the native opaque-conductor
  and weakly absorbing real-index dielectric domains still apply.
  `at_spectral_material_state` retains the original claims, schema identity,
  maximum knot gap and surface finish while continuing to a new accepted state.
  Source-selected wavelength changes must pass that same interpolation limit;
  SI conversion overflow, unsupported states and invalid grids refuse atomically.
  Schema identity also disambiguates two families adopted from one complete
  source bundle. Regression sources compare infrared complex Fresnel, phase and
  group indices, Snell refraction and intensity Beer-Lambert attenuation against
  independent formulae, including temperature refresh and a narrow absorption
  feature. These material bindings do not infer thermal emission or unsupported
  spectral tails.

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
  with an explicit source-band Kirchhoff law. `SpectralEmissivity` retains
  2–16 strictly increasing positive vacuum-wavelength knots with emissivity in
  `[0,1]`. Its existing `table` and `constant` constructors keep exact 380–780 nm
  support; `table_in_band` and `constant_in_band` retain the actual finite band,
  including infrared. `SpectralReflectance` shares this bounded numerical owner
  and the same explicit-band constructors. Evaluation never extrapolates.
  Visible-observer light proposals require full visible source coverage and
  restrict their quadrature to that observer band even for wider source tables.
  Each wavelength emits
  `epsilon(lambda) * B_lambda(T)` and reflects `1-epsilon(lambda)`. Source
  identity binds the actual numerical curve, temperature and declared source.
  A finite-band table does not supply bolometric exitance, out-of-band continuation,
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

- With `photometry`, `DiffuseThermalSurface::visible_observables` measures the
  original gray or borrowed spectral source in absolute units: 380–780 nm
  radiant exitance in W/m², CIE XYZ with Y in cd/m², Lambertian luminous
  exitance in lm/m², and optional xy chromaticity. Source identity and the
  published observer ID are retained. Four-point Gauss–Legendre integration
  splits at every emissivity knot and every existing 5 nm observer boundary;
  narrow accepted emission peaks therefore survive even when all original
  midpoint samples miss them. The operation allocates no heap storage and
  performs at most 376 pointwise source evaluations. Its CIE response remains
  piecewise constant inside the existing observer bins; tracer packet semantics
  do not change. Missing visible-band support and nonfinite outputs refuse.
  Zero represented XYZ has no invented chromaticity. Exitance is gross emission
  per area of one emitting face, not net exchange or an additional thermal debit.
  No display-white normalization, exposure, out-of-band cooling law, calibrated
  specimen or certified quadrature error is inferred. G1/G3 numerical tests
  compare independent SI references, emissivity scaling, gray/spectral limits
  and sub-bin emission. The `incandescent_emission` example reports watts,
  lumens and source color for an explicitly authored gray model.

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
  seed evaluation grid is checked and capped; a stable streaming top-k retains
  only the admitted Newton starts (six by default), with fallible reservation
  before evaluation. `NurbsRayBudget` and `ray_intersect_nurbs_with_budget`
  admit caller work units, start/iteration counts and retained-seed bytes.
  The work charge conservatively prices the source structure, full grid,
  retained Newton attempts, hull scan and seed ranking; the seed-byte allowance
  does not include the native evaluator's separately bounded scratch or source
  storage. Defensive sealed-source invariant failure,
  invalid settings, resource refusal, cancellation, and bounded-search
  nonconvergence remain distinct typed outcomes. Outward interval projection
  of the stored homogeneous controls encloses the entire positive-weight
  rational surface. Interval slab intersection against the original positive
  ray proves global misses, returned as `Ok(None)`, without a finite cutoff or
  a lost small direction component. A ray intersecting this control hull still
  needs Newton; exhausted starts fail closed as `IterationLimit`. The hull test
  does not certify Newton hit existence, nearest-hit completeness or trimmed
  geometry. Jacobian columns and residuals are normalized
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

  `charts::proximity` adds bounded point-to-surface queries over that same
  native mesh BVH. Interval separating-plane bounds provide conservative
  clearance, while convex surface witnesses independently bound proximity;
  a supporting-plane bound retains precision close to broad planar faces.
  Prepared `MeshPlanarPatches` borrows immutable geometry, validates retained
  BVH bounds once and shares them across repeated queries. Exact-coordinate
  duplicate vertices and consistently oriented coplanar neighbors form one
  patch, so internal triangulation diagonals do not block refraction.
  Patch edges and every other face, including a thin object's opposite wall,
  bound its departure neighborhood. Direct unprepared queries validate their
  authoring arrays before using cached bounds. Exact-distance charts can
  provide upper proximity; Lipschitz implicit fields provide only clearance.
  Optional `Chart::trace_normal_neighborhood` supplies an analytic regularity
  witness over an entire closed ball: a unit-gradient cone for an exact signed
  distance, excluding creases, medial axes and opposing sheets. `SphereChart`
  provides this witness with interval rounding and a curvature bound; its
  finite support is outward rounded. A point gradient alone never supplies
  this capability.
  The separate optional `Chart::trace_normal_differential` supplies the local
  analytic spatial Jacobian of its geometric unit normal. `SphereChart` uses
  `(I-nnᵀ)/|p-center|` at the queried point, including numerical contact collars;
  missing derivatives, creases and medial axes are unsupported. This is a
  numerical first derivative, not a global derivative-error certificate.
  The queries retain deterministic feature ties, cumulative work limits and
  cancellation. They return separate surface witnesses, never invented ray hits.

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

  `instances::proximity::BoundaryGeometry` prepares a static top-level BVH,
  shared borrowed mesh patches and interval affine transforms for native
  curved-path queries. Inverse-transform uncertainty and singular-value
  enclosures convert local clearance/proximity to world metres; stored
  binary64 quaternions are not assumed to be exact isometries. Queries keep
  original object/geometry/frame identities and one cumulative work budget.
  Resolved contacts exclude other objects and noncoplanar patches; coincident
  interfaces, unresolved finite edges and features within the contact collar
  refuse. A planar arrival's incidence-scaled unresolved tail must fit inside
  its exclusion ball. This prevents finite-detector near misses from becoming
  surface events. The prepared owner admits finite support. Smooth chart
  neighborhoods use interval inverse-transpose normal-cone transformation and
  world-metric signed-distance enclosures, exclude other objects, and refine
  within a fixed 32-query ceiling. Existing planar-witness fields retain their
  meaning; smooth regularity is a separate optional capability, currently
  native for `SphereChart`. Unsupported exact-distance charts retain numerical
  arrival proximity and refuse departure. These are separate witnesses, not
  extensions of the certified straight-ray `InstanceHit` authority.
  `normal_differential_with_cx` resolves a retained object/frame and exact query
  point through the same accelerated proximity owner, contact resolution,
  ambiguity checks and cumulative work budget. Planar mesh-patch interiors
  supply zero curvature. Chart derivatives transform with inverse coordinates,
  the inverse-transpose normal map and its normalization derivative; stored
  quaternion coefficients are not assumed to be exact isometries. The result
  retains its exact endpoint and geometry/frame/feature identities. Stale
  frames, patch edges, competing contacts and unsupported charts refuse; the
  path oracle exposes the same query using its admitted contact and budget.

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
  `TRANSFER_RENDER_SEMANTICS_VERSION = 3`.

  `render_absorption_transmittance` is the exact orthographic specialization
  for nonnegative, piecewise-constant intensity absorption in inverse world
  lengths. It sums the original column's optical depths with compensated
  arithmetic, using each z cell's nominal thickness, then evaluates the shared
  deterministic Beer--Lambert exponential once per pixel. There is no
  collision noise or optical-depth-dependent work. Explicit pixel, field-scan
  and pixel/cell-visit budgets admit the complete operation before allocation.
  All cells are validated, including cells outside sampled image columns.
  Pixel centres use exact integer cell ownership, independent of grid origin;
  lateral subpixel averaging, scattering and refraction are not supplied.
  Focused G0/G1/G3 tests cover heterogeneous layered columns, thin/vacuum/opaque
  limits, length-unit rescaling, layer refinement, translation, cell-face
  ownership and invalid-field/budget refusal. The original Monte Carlo APIs
  retain their sampling semantics and remain available.

- `volumes::transport` (features `volumes,tracer`) connects borrowed simulation
  density and temperature snapshots to physical scene transport. `AbsorbingGrid`
  requires kg/m³ density, metre geometry and an explicit source-domain mass
  absorption spectrum in m²/kg; raw lattice density is not silently interpreted
  as a physical optical coefficient. Each crossed constant cell integrates
  Beer--Lambert attenuation and, when a kelvin field is supplied, local-equilibrium
  emission `sigma_a * B_lambda(T)` using the same SI Planck owner as thermal
  surfaces. Front-to-back emission has correct intervening absorption, stable
  optically thin `expm1` evaluation and a repeated-temperature packet cache.
  Public and observed transport share compensated optical-depth accumulation.
  Thermal radiance is scaled before applying subnormal survival weights;
  subnormal cell absorption retains the original density, mass coefficient
  and path-length factors, preserving representable emission. Public APIs
  retain transparent/opaque limits, while observed APIs additionally require
  representable finite depth and positive emitted contributions.
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

## Material thermo-optic refraction

With `material-state`, `SchottThermoOpticIor` also carries the original
`fs-material::thermo_optic::SchottThermoOpticLaw` through homogeneous native
glass. Cauchy, Sellmeier, tabulated and explicitly air-corrected Sellmeier
reference laws retain their exact coefficients, source coordinates and
reference atmosphere. The shared material owner evaluates the frozen
temperature correction and its analytic spectral derivative; no heated-index
fit or finite difference enters transport. `schott_glass_at_temperature`
retains reference absorption and admits its common source band.

The same phase law reaches Snell/Fresnel surfaces, nested media, thin and
sequential lenses, coherent films and native bulk paths; group propagation
uses the separately evaluated group index. Temperature refreshes start from
the immutable reference, so a thermal cycle restores the same optical state.
The original thermal source, reference dispersion and frozen T enter native
material, scene, lens and film identities. Phase-only queries remain usable at
tabulated kinks while undefined group derivatives refuse. The source bands and
temperature bounds are not extrapolated; applying a second thermal correction
to an already heated law refuses. Endpoint admission and checked queries are
not a whole-band interval certificate. This composition supplies no thermal
expansion, evolving absorption, stress optics or new measurement authority.
G1/G3 tests cover independently differentiated TIE-19 response, critical-angle
motion, lens power, Airy film reflection, bulk phase/group/Beer propagation,
air-reference separation and temperature-cycle identity.

`thermo_optic::binding::MaterialThermoOpticBinding` carries this same law
through accepted material observations. Its original dielectric binding must
match the supplied reference state at exactly the law's reference temperature.
Every successor retains the same card, complete coordinate meanings and all
non-temperature coordinates, and reads the explicitly named typed static
absolute-temperature axis. The resulting ordinary `MaterialOpticalBinding`
carries the actual accepted-state identity into native tracer materials and
slab emission. Source identity includes the original binding and thermal law;
surface finish and reference absorption remain fixed explicit assumptions.
Heating/cooling always reevaluates the immutable original reference. An `Arc`
shares that source across endpoints; wavelength queries allocate nothing.
This does not admit dielectric volumes from open regional exterior meshes.

With `material-state`, `silica_refraction::corning7980_glass` binds NASA CHARMS
Corning 7980 fused-silica coefficients directly to the native Sellmeier glass
owner. Lenses, surface refraction, polarization and group propagation retain
the original absolute-index fit, its 30--300 K and 400--2600 nm domain, and
explicitly independent Beer absorption. The source class identifies only the
index fit; it does not authenticate a specimen or its absorption spectrum.

`MaterialSilicaRefractiveField` composes the same material law with the original
thermal snapshot and its analytic gradient/Hessian. Native rays and bundles
receive both `n_T Hess(T)` and `n_TT grad(T) grad(T)^T`; a zero central gradient
does not discard thermal focusing. Wavelength preparation is allocation-free
and caches absorption while retaining spatial temperature queries at every
stage. Bounded fields keep their original domain refusals. These derivatives
come from the published global index fit, not the source's separately processed
thermal derivative estimates; their endpoint discrepancy remains explicit.

`P1ConductionMesh::propagate_silica` and `propagate_silica_bundle` carry that
same source law through the original accepted transient conduction mesh. The
launch binding must borrow the exact mesh and nodal temperature slice. Shared
nodes preserve continuous phase index; element gradient jumps retain the native
crossing-time and bundle corrections without adding artificial Fresnel events.
The original mesh walker owns face adjacency, transition collars, full-path
budgets, cancellation, phase/group distance and absorption. Propagation adds no
heap workspace. Exterior or duplicated contact-node faces remain native
boundaries, and the CHARMS temperature/wavelength support is not extended.
Regressions cover an actual backward-Euler cooling endpoint, spectral ray and
bundle agreement with an independent affine field, gradient-jump sensitivities,
and contact/snapshot/source/work refusals.

With `photometry` and `material-state`, the silica constructors on
`NativeMaterialFieldBinding` bind these original fields to the existing native
scene, molecular-ambient scene, camera and differentiated-bundle observers.
Admission requires the complete released Corning glass at an explicit reference
temperature and exactly the field's independent absorption; an equal index or
provenance tag alone is insufficient. Actual local temperatures replace the
reference medium during observation. First/second-order and wavelength-prepared
bindings retain the same complete source/reference/snapshot identity. Geometry,
medium ownership, source domains and native transport remain authoritative;
no spatial reconstruction or contact smoothing is inferred. Three regressions
cover analytic thermal focusing at visible/IR wavelengths, temperature-driven
camera visibility/radiance, and source/reference/cancellation refusals. They are
authored and formatted; execution remains blocked by the upstream build resources.

With `material-state`, `MaterialThermoOpticField` binds the original reference
glass and `fs-material::thermo_optic::SchottThermoOpticLaw` to one immutable
thermal snapshot. The default constructor intersects reference phase/absorption
and thermo-optic support. Preparing a wavelength evaluates reference dispersion
and the spectral material coefficients once; optical stages still query the
actual spatial temperature and its analytic derivatives.

The original gradient-index and ray-bundle owners consume
`grad(n) = n_T grad(T)` and
`Hess(n) = n_T Hess(T) + n_TT grad(T) grad(T)^T`, including mixed terms.
Phase and available group optical distance retain the original material law.
The default absorption remains the explicitly declared reference law. The same
bound source supplies the native temperature trait. No finite-difference field,
moving state, expansion or stress response is invented. Each material/contact
discontinuity requires a separate region and the native boundary event.

`MaterialThermoOpticField::try_with_absorption_source` additionally borrows an
original `SpectralDielectricAbsorptionTemperaturePath`. Local refraction and
intensity absorption read the same temperature; native ray and bundle optical
depth integrates that alpha along the actual curved path. The source preserves
its original selected claims, auxiliary requirements, temperature breakpoints
and moving spectral knots. Preparation checks continuous source support once;
optical stages perform bounded searches with no allocation or material query.
The joined temperature domain and common spectral band remain explicit,
including for transparent spectra. Replaced reference absorption does not
restrict the independently retained phase and current-absorption bands.

The exact absorption source survives P1 mesh transitions, bundle gradient-jump
collars and curved contact reflection/transmission collars. Retained contact
arrivals require the same source; native scene identities bind its prepared
identity, while the default reference-loss identity is unchanged. The caller
associates phase reference, thermal law and absorption with the same specimen
and held environment. This retains the original weak-loss geometric-optics
model and frozen thermal snapshot; no complex-interface law, thermal evolution
or constitutive data is inferred. Regressions compare curved-ray/bundle optical
depth and mesh gradient jumps with independent analytic integrals, and check
source bands, scene identity, cancellation and the retained reference policy.

`thermo_optic::conduction::PlanarConductionTemperatureField` borrows one original
`GeneratedLayeredWallSolution` profile, with an explicit world origin and unit
axis. The conduction owner supplies temperature, Fourier gradient and
curvature; the bridge only rotates them. Prepared profile queries are O(1),
allocation-free and bounded by the original layer and source temperature span.
Contact jumps and exact stored face values survive. Even a constant-temperature
bounded layer cannot activate an unbounded uniform-field shortcut.

Controls run original conduction solves through analytic infrared dispersion,
thermal refraction and native bundle focusing, with independent path/footprint,
phase/group, Beer attenuation, mixed-Hessian, contact and domain references.

## Absorption-aware thermal source sampling

Temperature-field volume emission proposals include each cell's density,
maximum mass absorption over its original spectral knots, and in-band Planck
source proxy. Fully transparent cells receive zero probability before Planck
evaluation and do not dilute the uniform support mixture. The absorption
maximum is cached once per palette model; the existing 1/16 uniform mixture
over active absorbing cells retains support for cold and weak sources.

Sampling and density evaluation use the same normalized CDF and complete
radial-to-angular marginal. Physical LTE emission, Beer attenuation, occlusion,
and multiple-importance weights retain their original owners. A spatially
uniform non-indexed absorption scale cancels, preserving that path's previous
sampling arithmetic. Proposal semantics version 3 records the change.
The maximum-product weight is a bounded-cost proposal proxy, not an exact
spectral integral or a variance-optimal, visibility-aware distribution.

Controls cover transparent-cell exclusion, original interior spectral knots,
probability/PDF normalization, infrared and extreme-ratio support, unchanged
LTE transport, and the uniform-opacity compatibility path.

## Finite-thickness glass opacity and transmission measurements

`DielectricGlass::incoherent_slab_response` combines the original spectral
index and absorption laws with normal thickness, refracted path length and
every internal reflection between smooth parallel faces with identical exterior
media. `incoherent_slab_response_between` additionally admits distinct lossless
incident and exit indices, using the original incident ray to retain the
conserved tangential wave number at the rear interface. Both sum s and p power
independently before averaging. Front total reflection admits no internal heat;
rear total reflection retains absorption during every forward and returning
traversal. Reflection, transmission and true bulk absorption remain separate;
direct Poynting transmissions and `expm1` retain weak escape and heating. Equal
exterior indices preserve the symmetric calculation. Exterior indices are
caller-supplied at the actual wavelength; original bulk support remains enforced.

`BeerLambertAbsorption` can import spectral internal transmission at a known
path length, or normal-incidence external slab transmission under that same
uncoated incoherent model. The latter accounts for both surface reflections
and all internal-bounce orders when recovering the native intensity absorption
law. External reflection must not be relabeled as bulk heat. Original
wavelength support remains bounded and alpha interpolates linearly; zero
transmission is a measurement limit rather than an infinite inferred
coefficient. Impossible above-clear values refuse. The calculation adds no
measurement or uncertainty authority.

Controls compare independent asymmetric interface bounce sums, power closure,
polarization, reciprocal transmission, material-bound infrared opacity after a
supported temperature refresh, thickness scaling, weak/opaque limits and inverse
measurements propagated through the original Beer law at new thicknesses and angles. This retains the
native real-index, weak-absorption interface approximation: no coherent
thin-film interference, strongly absorbing complex-index faces, rough faces,
scattering, coatings, lateral exit location or traversal delay. See
[SCHOTT TIE-35, section 1](https://media.schott.com/api/public/content/110b56518128431b9f3134627d7f4c4e?v=012d3e7e)
for the distinction between internal and external transmittance.

### Slab emission and native thermal exchange

`DielectricGlass::prepare_incoherent_slab_absorption_between` resolves the same
spectral and interface laws once for allocation-free depth-interval queries.
`PreparedIncoherentSlabAbsorption::absorption_between` integrates positive
forward and backward Beer losses with every internal reflection, separately
for s/p incident power. Normal-depth partitions create no optical interfaces.
Unequal exteriors and rear total reflection retain their original power laws;
front total reflection, clear materials and empty intervals absorb zero.
`expm1` retains weak cell heating and positive path decomposition preserves a
finite entrance skin even when the full slab is opaque. Opposite-face queries
retain the authored interval width without subtracting two reversed depths,
so a very thin physical cell is not rounded away inside a thick slab. G1/G3 tests compare
independent bounce integrals, additive partitions, metre scaling and opaque
skin localization. Optical coefficients remain homogeneous and frozen; no
thermal storage, within-slab material variation or coherent profile is inferred.

`SlabDepthPowerResponse` in `thermal_surface::heating::slab_depth` compiles
complete authored normal-depth partitions into the original spectral R/T/A
account. Reverse illumination changes traversal while retaining authored cell
order. `SurfaceOpticalHeating::try_new_slab_depths` assigns those positive
fractions to existing native thermal cells and reuses scheduled intercepted
watts, receiver material/temperature checks, source accounting and rollback.
Depth intervals may share a receiver; incident, reflected and transmitted
energy are still counted once per physical source. Receiver mapping and cell
volumes are caller-authored; generic receivers do not certify slab geometry.
Explicit per-response and whole-source binding ceilings bound storage;
preparation polls cancellation.
Focused G1/G3/G4 tests cover independent local heat, weak-loss enthalpy,
reversed localization, finite ramp schedules, shared receivers and atomic
refusals. Optics remain homogeneous and frozen on declared receiver support;
only prompt thermalized absorption is deposited, with no inferred cooling.

`CoatingThermalEmission::try_depth_resolved_incoherent_slab` observes one
temperature per authored slab depth interval. Its outgoing s/p radiance sums
each interval's reciprocal absorption times its own Planck source, retaining
the original homogeneous Fresnel/Beer owner without an averaged temperature or
invented internal optical interfaces. The native constructor snapshots each
bound thermal receiver's actual temperature and material. Depths are ordered
from the observed vacuum face; at most 64 intervals are admitted, with no
allocation during spectral/angular queries. Temperature refreshes preserve
optical cache identity while changing emission identity. Independent local
Beer/Planck controls, isothermal reduction, reversed temperature profiles,
snapshot immutability and domain/cancellation refusals cover this path.
This remains finite-band outgoing appearance at frozen homogeneous optics;
internal radiative exchange, net cooling and thermal geometry are not inferred.

`CoatingThermalEmission::try_incoherent_slab` retains one uniform-temperature
bulk-glass receiver and the original finite-band phase/absorption laws. Native
slab s/p absorptance supplies reciprocal emissivity; transmitted power remains
radiation, and transparent slabs retain a thermal channel with zero emission.
The selected shared exchange port remains vacuum. The separate
`DielectricGlass::incoherent_slab_thermal_response_between` observation supports
either lossless exterior with that medium's index-squared radiance convention
([CIE S 017:2020, 17-24-005 note 2](https://cie.co.at/eilvterm/17-24-005)).
Existing one-port, layer-channel and two-vacuum-port exchange owners reuse the
same reflection, transmission and absorption.

`MaterialOpticalBinding::incoherent_slab_emission` requires the exact accepted
state, an existing typed absolute-temperature axis, and smooth uncoated glass.
`IsothermalCoatingStack::IncoherentSlab` and
`NativeCoatingEmission::try_incoherent_slab` connect the same opacity to native
caloric storage, step restrictions, cached quadrature and signed energy
accounting. These fixed-optics entry points retain their original optical laws
over the admitted temperature interval; material-dependent optical changes use
`MaterialSlabLaw`, described below. No bolometric continuation,
coherent/strongly absorbing interface model, rough face or geometry-derived
thickness is inferred. Focused
regressions exercise directional emission, transparent/thin/black limits,
material-state refresh, two exposed faces, equilibrium and native energy closure.

## High-index transparent dielectric transport

Explicit finite-band Sellmeier and tabulated phase laws admit ordinary real
indices above three, including high-index infrared windows and lenses. The
native scalar/polarized Fresnel, smooth/GGX transmission, GRIN sample and
transmitting-coating exterior checks share `n >= 1` and
`n² <= f64::MAX/64`; this upper condition reserves arithmetic headroom and
supplies no material or transparency claim. Original wavelength support,
absorption laws, and each coherent-film owner's additional limits remain.
The visible Cauchy, sample, table and legacy Sellmeier constructors retain
their original `[1,3]` domain and existing accepted arithmetic/identities.

Exact normal incidence retains the transmitted normal square directly,
without cancellation against a much larger incident index. Evanescent Fresnel
phase ratios scale only when their intermediate squares are unrepresentable,
so total internal reflection remains a finite polarization transformation.
No metal, anomalous-index or complex-Snell model is inferred from this range
extension. Controls cover independent high-index Snell/Fresnel and Brewster
values, reciprocal rough transmission, quarter-wave matching in both
directions, native analytic GRIN paths/group distance/Beer attenuation, source
bounds and normal/evanescent extreme-contrast arithmetic.

## Refraction from nonlinear conduction

`material_state::thermo_optic::nonlinear_conduction::PlanarNonlinearConductionTemperatureField`
embeds one original `NonlinearGeneratedLayeredWallSolution` layer using an
explicit origin, unit axis and execution context. It consumes the conduction
owner's temperature and Fourier derivatives, including the original `dk/dT`
curvature term, and rotates the gradient/Hessian into world coordinates. The
existing thermo-optic, GRIN and bundle owners then retain phase bending,
chromatic focusing, available group distance and the declared Beer law.

First-order optical queries remain usable at continuous conductivity knots.
Second-order queries require the physical curvature capability of that source;
they refuse a nonzero-flux interior slope ambiguity, admit zero-flux extrema,
and use inward derivatives at layer faces. Every query retains the original
layer/source bounds and cancellation. There is no extra thermal solve, spatial
fit, contact smoothing or inferred transient evolution. Controls start with an
actual nonlinear `k(T)=T/300` solve and check native infrared focusing at three
wavelengths against analytic temperature and ray-variation references.

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
terminal success. Cancellable NURBS tracing carries the same `Cx` through native
surface admission, basis/contraction evaluation and partial derivatives, checks
each control during outward hull construction, and polls before and after each
bounded seed/iteration. Final cancellation wins over either a hit or a proved
miss. BVH traversal polls before and after each bounded node/triangle.
Differentiable scanline rendering
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
expectation. Tracer bit-semantics v4 and bidirectional semantics v10 retain
ordinary emitting geometry's authored surface scattering. Source radiance is
added at a hit with its existing emitter MIS weight; reflection/transmission
then uses the original diffuse, GGX, conductor/coating, brushing or dielectric
owner. An emitting surface can be an intermediate vertex on either BDPT
subpath; its emission is a separate source endpoint contribution for the
corresponding camera prefix. Emitting dielectric solids also participate in
initial-medium discovery and Beer absorption. An exact black
`SpectralDiffuse` material supplies an opaque emission-only surface. These
authored relative-radiance sources do not infer an equilibrium thermal law;
thermal and calibrated photometric source ownership remains unchanged.
The five emitter-continuation regressions are authored but not yet executed:
the native feature build was killed while compiling the pinned runtime under
the available memory limit. The Cornell reference remains pinned to tracer
v3 and must be reproduced in debug and release before a justified re-freeze;
neither its digest nor its evidence pin has been replaced without that proof.

Tracer bit-semantics v7 uses the same all-sphere orthonormal frame for diffuse
and reflective GGX sampling. The diffuse cosine warp retains the actual unit
surface normal near the negative Z pole; it neither substitutes fixed tangents
nor repairs a distorted distribution by normalizing the result. Returned rays
stay unit length to floating-point roundoff, and their sampled/evaluated
solid-angle densities agree. This preserves the orientation-independent
Lambertian furnace weight in ordinary and spectral diffuse materials, their
thermal reflection counterparts, and native camera and BDPT consumers.
The frame has constant work and no allocation or rejection loop. Focused
regressions exercise both sides of the former pole cutoff, grazing directions,
an independent rotated quadratic-illumination integral, and forward/reverse
BDPT densities. The existing Cornell digest and its evidence pin are unchanged.

Current scalar scene-integrator no-claims: no gradient-index propagation or arbitrary overlapping media,
no Russian roulette, reflective GGX
and conductor sampling is isotropic-only, conductor and rough-dielectric GGX
are single-scattering without multiple-scattering compensation. Thermal
emitters retain their reciprocal diffuse reflection. Dielectric tracer
bit-semantics v14 retains v13's isotropic
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

Directed beam power uses compensated cumulative optical depth, including
absorption and scattering, before applying attenuation to the segment's
original incident power. Refining a layer into cells whose individual survival
factors round to one therefore retains the resolved total loss. Per-cell heat
and scattered watts still use their separate positive `expm1` partitions;
small transmission fractions are scaled without discarding representable watts.

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

### Prompt fluorescence and absorbed heat

`volumes::heating::fluorescence::{FluorescenceSource, FluorescentOpticalHeating}`
adds prompt Stokes fluorescence to original exact beam deposition and the same
native enthalpy observer. The supplied source retains its identity, excitation
and temperature domains, photon quantum yield, and 1–16 discrete emission
nodes with normalized photon weights. Escaping power uses each emitted photon's
energy relative to the absorbed pump photon. Positive nonradiative and Stokes
losses alone heat mapped cells, preserving small heat without subtracting
nearly equal radiant powers. The full absorbed pump account remains separate;
native accepted substeps integrate heat and escaping radiation transactionally
on the same clock, without repeating optical traversal.

The absorption-only grid must describe the complete fixed, open, index-one,
non-scattering optical environment. Every occupied cell must have exactly zero
absorption at every emitted wavelength, including cells outside the pump ray.
An absorption view stripped from a reflecting, refracting or scattering scene
does not establish this model. Self-absorption, unsupported wavelengths,
anti-Stokes emission, mismatched medium identity, scattered pump power,
excessive preparation counts, cancellation and temperature-domain escape
refuse. This consumer does not infer fluorescence lifetime or excited-state
storage, saturation, bleaching, cascades or calibrated source data. Independent
photon/energy closure, native temperature response, thin/opaque and yield limits,
one-ULP Stokes heat, domain refusals and transactional rollback have focused
regressions in `tests/optical_heating/fluorescence.rs`.

`volumes::transport::luminescence::LuminescenceGrid` observes original prompt
or delayed cell emission as nonthermal SI radiance. It borrows the original
absorption grid and prepares a bounded sparse map of 1–16 discrete-node powers,
with source, wavelength and numerical snapshot identities. Each cell declares
uniform isotropic emission `j = P/(4*pi*V)`; the native DDA integrates emission
and foreground Beer attenuation analytically, including the transparent limit.
Output is line/bin-integrated W/(m² sr), without an invented spectral width or
temperature. Rays allocate nothing; preparation and original cell traversal
retain explicit limits and cancellation. Thin/opaque limits, grid refinement,
direction-dependent shielding, extreme but representable radiance, and native
fluorescence power/heat consistency have G1/G3 tests in
`tests/optical_heating/luminescence.rs`. This is a stationary non-scattering,
index-one observation component; it does not feed absorbed secondary light
back to heat or solve fluorescence cascades, subcell pump structure, dipole
anisotropy, surfaces or light-travel delay. Thermal emission remains separate.

### Accepted molecular opacity in luminescence observations

`LuminescenceGrid::with_additional_absorption` binds a sparse original-cell
map of absorption coefficients in inverse metres at the declared emission
nodes. This admits the accepted `sigma * ground_count / volume` opacity of
finite emitters alongside the host's passive mass absorption. Emission and
background attenuation use the same `rho*kappa + alpha_active` coefficient;
absorbing cells need not emit. Host and active contributions must not overlap.
The cell integral retains the transparent and opaque limits, including finite
radiance when total extinction overflows its direct floating-point product.
Each ray retains its original traversal budget and allocates no opacity map.

The binding validates original cells, nonnegative finite coefficients, a source
identity and a preparation cap before publishing an immutable snapshot.
Source powers and opacity must describe the same accepted physical state;
separate interval means do not establish a time-integrated exposure when
opacity changes. The observer evolves no populations or heat. The multiple-
scattering adapter includes the same additional opacity in its flight density,
spectral continuation weights and complete source integral; see the scattered
nonthermal observation contract below.

Three focused G0/G1/G4 regression groups cover exact layered attenuation,
unchanged inactive spectral lanes, representable opaque-cell radiance,
nonemitting foreground absorbers, invalid bindings and cancellation. They are
source-formatted and reviewed, but remain unexecuted: the genuine renderer
build was killed in `asupersync` before renderer typechecking.

### Camera images of nonthermal luminescence

`volumes::transport::luminescence::render_luminescence_camera` samples original
physical pinhole or ideal thin-lens rays and the same native cell integral.
A complete image retains original wavelengths and source/field identities,
with per-node mean integrated radiance, transmission, thin opacity and effective
negative-log survival. Constant-footprint backgrounds compose spectrally,
including finite bright backgrounds behind underflowed unit transmission.
Pixel and pupil samples average transported light independently at each node;
they do not average density or optical depth before evaluating transmission.
Shared image sampling limits and an actual cumulative DDA-cell allowance remain
binding, including transparent emitting cells, with no partial output.

`approximate_emission_xyz` uses the existing analytic visible observer fits and
relative-renderer normalization without an inferred bin-width factor. Coarse
supplied bins carry no spectral accuracy claim; UV/IR data remain available in
the primary spectral image and refuse the visible-only projection. Dispersive,
diffractive and sequential lenses retain native refusal. Thermal emission,
scattering, surfaces and secondary fluorescence remain outside this observer.
Eight focused regressions are in
`volumes/transport/luminescence/camera/tests.rs`; native compilation and execution
of this new adapter remain unverified because the execution environment became
unavailable. The exact adapter and parent change received independent source
review against the existing physical-camera, cell-traversal and spectral APIs.

### Spectral absorption and thermal camera images

`volumes::transport::camera::render_absorbing_camera` connects the original
absorption grid and its accepted thermal or independently thermalized
constituent sources to physical pinhole and ideal thin-lens images. Each ray
uses the same native cell integral, preserving original wavelength support,
foreground shielding, Kirchhoff emission and original field identity. No
mean constituent temperature, emissivity, spectrum or RGB absorption is
invented. Per-node foreground and composed background radiance remain
W/(m² sr nm), with separate transmission, thin opacity and effective depth.
Background composition requires a constant background within each averaged
pixel/pupil footprint. Optional approximate XYZ requires explicitly supplied
finite nonnegative quadrature weights in nm; the caller owns spectral support,
resolution and integration accuracy. UV/IR spectral output remains available.

The raster/pupil sampler, stable spectral average, background composition and
visible projection are shared with luminescence. Its public API, integrated-bin
radiance units and original random-stream domain/draw order remain intact.
Actual cell counting is independent of the native integrator's strict observed-
column policy: camera observations retain the public thin/opaque limits while
every existing strict numerical refusal remains unchanged. Both cameras retain
native lens refusals, cancellation, per-ray crossing limits and cumulative
whole-image work bounds, with no partial image on refusal. They remain
stationary, non-scattering, index-one observations without surfaces or exposure.
Six new regressions cover independent SI Planck/Kirchhoff equilibrium,
constituent temperatures, foreground order, opaque-limit policy, explicit nm
quadrature and cumulative limits. The existing eight luminescence fixtures are
retained apart from import updates. All changed sources parsed with the pinned
formatter and received independent interface/physics review; native tracer
compilation and execution remain unverified because its dependency closure is
unavailable in the storage-constrained execution environment.

### Finite-lifetime luminescence on the native heat clock

`volumes::heating::luminescence::{LuminescenceSource, LuminescentOpticalHeating}`
extends the same admitted pump/transparent-escape model with a separate,
initially unexcited energy store for each absorbing optical cell. An authored
emitting-state wavelength fixes the stored energy per excitation; it must lie
between the entire pump band and every emitted photon energy. Pump-to-state
relaxation heats the mapped native cell immediately. The L3
`fs_thermal::luminescence::ExcitationDecay` then integrates constant pumping and
decay analytically over the final accepted interval. Radiative branching retains
the original photon fractions and photon-energy ratios; nonradiative decay and
the remaining radiative relaxation energy heat that same cell once.

Excitation is not thermal enthalpy. Cumulative absorbed energy equals escaped
light plus native heat plus the remaining excitation, within explicit joule
roundoff allowances. Every candidate retains signed global and cumulative
absolute local defects. Native explicit steps, implicit interval retries, event
alignment and observer composition keep acceptance with the original thermal
owner. Retrying an interval replaces tentative work from the unchanged accepted
store. Cancellation and source-temperature refusal publish no partial history.

The existing `SpectralPowerSchedule` supplies bounded rectangular pulses and
dark gaps; its transitions split native steps exactly. Ramps refuse because a
fluence-only mean would misrepresent excitation timing. Scheduled lanes remain
within the prepared linear Beer-law envelope; unrepresentable small scaling
refuses rather than losing transmitted or absorbed energy. Pump-off intervals
continue to emit and heat from the retained excitation. Optional authored
Arrhenius quenching adds a nonradiative rate and changes lifetime and yield
together. Rates freeze at each native start temperature, so thermal feedback
still requires independent timestep refinement.

Accepted endpoint cell-emission watts and last-interval mean watts are separate
observations. Both retain original optical cell indices and integrate discrete
emission nodes, allowing `transport::luminescence::LuminescenceGrid` to observe
the actual excitation history without an invented temperature or bin width.
This remains one linear unsaturated emitting state with fixed pump optics,
transparent escape and supplied material data. It does not infer a spin/trap
phosphorescence mechanism, saturation, bleaching, self-absorption/cascades,
temperature-dependent absorption or measured rate constants. Focused G0/G1/G3/G4
regressions cover analytic pumping and afterglow, photon/heat/storage closure,
prompt limits, quenching, temporal partition/restart, implicit retry replacement,
shared thermal receivers, source-domain refusal and transactional cancellation.

### Singlet/triplet phosphorescence on the same thermal clock

`LuminescentOpticalHeating::with_singlet_triplet` explicitly binds a
`phosphorescence::SingletTripletSource` to the existing singlet pump owner.
The authored lifetime form retains independent intrinsic lifetimes, photon
yields, quenching laws and normalized photon spectra. Their common discrete wavelength
grid can contain zero-weight nodes in either spectrum; only populated spectral
nodes must lie below that state's excitation energy. Neither yield is applied
twice: radiative, dark and intersystem-crossing transitions compete in the
existing L3 `fs_thermal::luminescence::transfer::SingletTripletDecay` model.

Absorbed pump photons populate the singlet. Forward crossing transfers actual
populations to the separate triplet store and releases the authored energy gap
as heat. Explicit Arrhenius reverse crossing removes that gap from the native
thermal bath and permits delayed singlet emission. Both stores continue to
evolve in dark intervals. Each state's radiative release is converted through
its own photon fractions and photon energies, so the observed spectrum can
change during the afterglow without an invented temperature or RGB fade.

The source exposes separate accepted `[singlet, triplet]` joules and reports
downhill crossing heat and uphill bath absorption separately. Native external
heat includes the former and subtracts the latter. Absorbed pump energy equals
total stored excitation plus escaped/unresolved light plus net native heat,
within the supplied absolute local and global joule allowances. Passive
emission-band absorption uses the existing directed transport operator. The
equal-energy scalar resonant recycler remains exclusive from singlet/triplet
coupling. `with_singlet_triplet_recycling` instead admits explicit spectral
photon return through `PhosphorescenceReceiver` probabilities. Every positive
receiver branch must pay the destination-state energy, and one absorbed photon
creates at most one excitation. Explicitly mapped unpumped receivers acquire
real stores; omitted receivers stay passive. Source Stokes heat, incoming
stored excitation and receiver excess/passive heat remain separate. Repeated
returns evolve within the same interval without a bounce cutoff. The current
native 64-state bound permits 32 active two-state cells; passive absorbers add
no kinetic state. Immutable transport and prepared plans are shared; actual
temperature-refreshed rates and interval duration determine plan reuse. The
source identity binds the original transport snapshot, thermal mapping and
explicit receiver probabilities. Absolute transport and receiver-partition
defects are charged before they can cancel between distinct excitation cells.
Focused cases in `tests/optical_heating/phosphorescence_recycling.rs` cover
self-return, unpumped active receivers, passive receiver heat, closed photon
cycles, energetically forbidden activation, retry and cancellation.

Fixed-rate intervals use the existing bounded positive matrix-exponential
owner, with a prepared plan shared across clones and reused for unchanged
rates and duration. Rates freeze at each native start temperature; changing
thermal feedback still requires timestep refinement. Pulse boundaries, retry
replacement, cancellation and final acceptance remain owned by the original
thermal observer. `emission_source_identity()` binds both source identities,
spectra, state energies, intrinsic decay/quenching and crossing laws.

G1/G3/G4 regressions in `tests/optical_heating/phosphorescence.rs` compare
population/energy transfer with independent Bateman solutions, observe dark
afterglow through the existing radiance grid, exercise reverse-crossing bath
exchange and thermal quenching, and check partition/restart, retry and refusal.
The unextended crossing mechanism is dilute and first order. Optional
`with_triplet_annihilation` binds a `PhosphorescenceAnnihilation` law to each
explicitly selected active optical cell. Its pair-event coefficient, active
population volume and singlet yield remain authored inputs; the renderer
never substitutes thermal mass or volume. Each event consumes two triplets
and creates at most one singlet, with per-event energy admission. The
separate `annihilation_heat_j` account also enters the actual native thermal
source and global energy closure.

Both independent cells and full spectral photon-recycling networks use the
native controlled nonlinear owner. Original linear rates still refresh at
each cell's actual temperature, and absent or zero-rate pair laws retain the
original linear path. Independent cells incur no coupled-network state cap;
nonlinear accuracy/work controls apply per independent cell or once to the
complete coupled network. Bind photon recycling before annihilation when
unpumped active receivers are needed. Law getters and source identity retain
the explicit coefficients, volume, yield, cell mapping and numerical controls.
The density and productive-afterglow regressions include independently
derived steady-state and external ODE references; these are manufactured
physics fixtures, not measured material data.

The `phosphorescent_afterglow` example accepts `--triplet-pairs` and reports
pair heat together with singlet and triplet radiance. None of these paths
claims finite ground-state depletion, trap-mediated persistent phosphors,
coherent spin dynamics, oxygen transport, measured material constants or an
automatically enforced detailed-balance law.

With `material-state`, `SingletTripletSource::try_material_temperature`
binds a complete `fs_material::phosphorescence::PhosphorescenceStatePoint`
and its immutable original card directly to this same owner. The six native
rates remain authoritative, including exact zero rates and stable stores.
Spectrum descriptors supply conditional photon fractions and validity domains;
their legacy scalar quantum yields are not applied to the native rate model.
No finite lifetime or branching yield is fabricated. The returned singlet
descriptor must be paired with its original coupling.

The native thermal owner controls one explicitly typed static absolute-
temperature coordinate. It refreshes the original selected material claims at
the actual cell temperature while retaining every other supplied condition,
including quencher concentration. The same original state energies must remain
fixed: selected energy models depending on that temperature coordinate refuse
at binding, and every refreshed snapshot rechecks the original energies.
Changing an occupied energy level would require separate parameter-work
accounting. Source-domain or identity mismatches refuse.

Temperature snapshots are cached within an interval with a cell-count-bounded
cache; numerical plans still reuse the original positive-network propagation.
The cache carries no physical history and cannot change accepted populations,
heat or source authority. The material card, initial complete context and owned
temperature coordinate contribute to the emission source identity. Manufactured
integration regressions check stable triplet storage, actual-temperature rate
refresh, fixed quencher conditions affecting both afterglow and heat, and
unaccounted energy-work or mismatched-singlet refusal. The runnable
`phosphorescent_afterglow` example demonstrates a synthetic finite pump pulse,
nanosecond singlet and second-scale triplet through native heat and radiance.

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

### Photon momentum on native volume histories

`ScatteringGrid::deposit_multiple_scattering_momentum` and the corresponding
`CellScatteringGrid` methods return the original power report together with
per-cell radiation force and ballistic, diffuse and unresolved photon momentum
rates. The same native free flights and wavelength/cell phase samples transfer
`[P_abs d_in + P_continue (d_in - d_out)] / c`; elastic scattering exerts force
without creating heat. Absorption-only lanes retain exact spatial Beer transport.
The force query samples the final capped collision's outgoing direction and
retains that packet as unresolved. It introduces no subsequent force or heat
from that unresolved continuation.

Signed cell/history reductions and the mean momentum residual use compensated
sums. Total-force statistics retain correlations between cells; component
standard errors are fixed-sample estimates, not confidence bounds. One bounded
cell lookup handles three components. Heating-only calls allocate no force
tallies, and independent history streams preserve their original power results.
Cancellation and force-storage refusal publish neither account.

The physical domain is stationary, open, index-one scalar transport with
quadrature-weighted incident watts. Refractive-host stress, moving-medium Doppler
work, thermal-emission recoil, spin angular momentum and mechanical time
integration are outside this API. Seven focused regressions cover exact
absorption, elastic isotropy, forward/backward scattering, individual-history
momentum closure, local spectral phases, diffuse escape, repeatability and
transactional refusal.

### Polarized multiple-scattering heat and force

With `tracer`, `ScatteringGrid::deposit_polarized_multiple_scattering_power`
and `deposit_polarized_multiple_scattering_momentum` accept a common physical
reference frame and quadrature-weighted full-Stokes watts at each wavelength.
Every history retains its actual polarization through the native Rayleigh,
depolarized-Rayleigh, retained spherical-Mueller, spectral and mixture phase
owners. Angular proposals use the existing scalar density and reject against
the full Mueller intensity with envelope `1 + hypot(Q,U)/I`; accepted packets
retain their original power and carry the normalized outgoing Stokes state.
Circular and unpolarized inputs need no rejection draw. The original native
free flights, implicit capture, original-cell heat identities, force transfer,
compensated statistics and directed unresolved accounts remain the owners.
The common power report feeds `ScatteringHeating` without another heat debit.

`deposit_polarized_continuous_absorption_power` and
`deposit_polarized_continuous_absorption_momentum` use continuous absorption
weighting through that same source admission and matched Stokes phase owner.
The admitted total scattering coefficient is polarization independent, so only
scattering determines free flights. Isotropic absorption deposits positive Beer
loss along every original segment while leaving normalized Stokes state intact;
it is not charged again at a scattering event. The common history reducer
carries the evolving polarization into every angular continuation. Its momentum
variant also scores segment absorption and actual polarized elastic recoil,
including the outgoing direction at a capped event. Scalar and polarized
continuous queries share one DDA, conservative accounts and refusal paths.
Tests compare the two absorption estimators, weak/no-collision heat, physical
orthogonal-dipole receiver patterns, spectral replay and proposal/cancel refusal.

`PolarizedDepositionConfig::max_phase_proposals` bounds each angular draw.
Cancellation, invalid/nonphysical states, unsupported scattering laws and
proposal exhaustion return no partial heat or force report. Absorption-only
lanes retain exact Beer deposition, including scalar-only phase declarations
that are never used. Scope is stationary, index-one transport with isotropic
extinction and randomly oriented molecular/spherical scatterers; aligned or
dichroic media, refractive-host stress, optical spin torque and moving-medium
work remain outside these APIs. Focused regressions specify analytic dipole
angular moments/PDFs, dark-proposal rejection, exact absorption, spectral
multicollision heat/momentum closure, orthogonal-polarization spatial heating,
transactional refusal and transfer to original thermal enthalpy cells.

### Orbital torque on native volume histories

`deposit_multiple_scattering_wrench` on both `ScatteringGrid` and
`CellScatteringGrid` adds orbital angular momentum to the same native power and
linear-momentum reports. Its explicit finite world-coordinate reference point
is in meters. Actual collision locations transfer
`(r_collision - r_reference) × F_collision`; every original cell receives the
sum of its physical event moments, not its center crossed with net force.
Absorption-only lanes retain exact Beer power and the source line of action:
free-flight displacement is parallel to momentum, so its cross product
vanishes and no absorption-position approximation is needed.

Ballistic escape, diffuse escape and directed unresolved packets each retain
orbital momentum flux. The incident-minus-outgoing-minus-torque account closes
about the same reference. Capped radiation never contributes untraced torque.
FMA-compensated cross products and signed history/cell sums retain small
moments under cancellation; total standard errors combine each history before
averaging. Original force/power queries allocate no orbital cell storage, and
the wrench query uses the existing work and spatial budgets without additional
phase draws. Invalid references, unrepresentable moments, cancellation or
budget refusal publish no partial account.

This resolves extrinsic orbital torque of scalar packets and off-axis beams,
not polarization/spin torque, phase-resolved optical vortices, moving-medium
work or mechanical time integration. Six additional regressions cover exact
off-axis absorption, actual collision locations and individual-history closure,
reference translation across escape/unresolved channels, local spectral cells,
small cross-product cancellation residuals and transactional refusal.

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

`ScatteringGrid::deposit_continuous_absorption_power` and the corresponding
`CellScatteringGrid` methods provide continuous absorption weighting for scalar
heat. Scattering alone determines the sampled free flight. Each crossed native
cell contributes positive analytic Beer absorption along the portion actually
traversed, before the surviving radiation samples the cell's original phase.
No absorption/albedo debit is applied again at the scattering event. A single
DDA pass per flight retains the residual scattering optical depth across cell
boundaries; it needs no temporary cell list or copied medium. Compensated
cumulative absorption depth preserves weak losses over many cells and flights.
This estimator follows Hayakawa et al., JOSA A 31 (2014), 301-311,
DOI:10.1364/JOSAA.31.000301.

The report retains `AbsorptionEstimator::{Collision, Continuous}` independently
of the physical medium identity. Exact absorption-only lanes, source streams,
history-level covariance, work limits and atomic refusal share the original
reducer. Surviving power at the last admitted scattering event remains
unresolved; truncation is not credited as heat. The same native
`ScatteringHeating` consumer admits either report. Collision weighting remains
available, and the existing force/polarized queries retain that estimator.
The continuous momentum/wrench variants integrate positive absorption impulse
along every traversed native cell segment and elastic recoil at each original
scattering position, within the same histories. Distributed absorption retains
the actual ray line of action for orbital torque; no cell-center approximation
or extra force traversal is used. The final capped event samples its outgoing
phase when momentum is observed, closing recoil with directed unresolved
radiation. Power reports match the continuous heat query, and force/torque
retain original per-history covariance, reference geometry and finite-rate
refusals. Native heat admission consumes only the absorbed-power account.
Continuous weighting can reduce heat noise in thin absorbers; no universal variance or
runtime improvement is claimed. Focused regressions cover the exact
absorption-only limit, layered first-scatter integrals, an independent analog
cube reference, native-cell subdivision/replay, prepared cell sources,
budget/cancellation refusal and accepted native enthalpy input.

### Cell-dependent absorption, scattering and phase

`scattering::CellScatteringGrid` borrows the original density/geometry and fixed
spectral packets of mass absorption and mass scattering for every physical
cell. A shared phase law or indexed phase palette retains the original angular
owners without copying large tables into cells. Explicit wavelength support
and source identities bind the prepared field; no constitutive data or spectral
extrapolation is inferred. Density stays in kg/m³ and coefficients in m²/kg.

Uniform and cell-dependent media share one monomorphized DDA/free-flight kernel,
one stable directed-power integrator and one implicit-capture history reducer.
Each collision uses its actual cell's albedo and phase. Shadow extinction,
local LTE emission and directed heat use those same coefficients. Absorption-
only lanes retain exact Beer transport and zero sampling error even when other
lanes scatter. Solver checkpoint adapters preserve the original refusal; all
existing traversal, history, segment and output-cell budgets remain effective.

Tests compare bitwise uniform-media reductions, independent layered spectral
collision laws, infrared support, local LTE absorption/extinction balance and
an independent analog two-layer calculation with distinct anisotropic phases.
The prepared packet is fixed for one call; temperature/state evolution belongs
to the hosting constitutive and thermal owners.

### Temperature-dependent multiple-scattering heating

Under `optical-heating`,
`volumes::heating::scattering_feedback::TemperatureScatteringHeating` couples
local temperature-dependent absorption **and** scattering to the original
`DepositionSimulation` clock. It owns one immutable original geometry/density
snapshot, shared compiled source spectra, a phase palette, and complete
original-grid/material/thermal mappings. Every positive-density cell requires
a binding, including initially transparent cells and cells outside the incident
beam. `TemperatureScatteringTable` retains actual mass-scattering spectra;
it shares the existing bounded temperature interpolation with
`TemperatureAbsorptionTable` without relabelling absorption as scattering.
Spectral laws are sampled once at the fixed source wavelengths, while the
native caloric laws supply actual kelvin values at every thermal substep.

The immutable trajectory can select `AbsorptionEstimator::Continuous` through
`try_new_with_estimator` or `from_field_with_estimator`. It samples scattering
free flights and scores conservative Beer absorption on every traversed native
cell, so a weakly absorbing path contributes heat even without a collision.
The original constructors retain collision capture. Both choices share the same
compiled physical source tables, temperature preparation, thermal admission,
cancellation and aggregate work bounds; estimator identity stays explicit on
each original deposition report. Replica clones, changed history counts and
whole-trajectory refinement retain the chosen estimator. Continuous absorption
does not universally reduce variance, alter the admitted material laws or
remove the separate time/sampling/unresolved-energy accuracy gates.

`from_polarized_field` binds an `IncidentPolarizedScatteringBeam` to the same
compiled field and native thermal observer. Original frame/Stokes admission and
molecular, depolarized-molecular, spherical or supported mixture Mueller laws
supply the physical angular response. `CellScatteringGrid::deposit_polarized_power`
and its hosting-checkpoint variant share those histories with either native
absorption estimator. Isotropic absorption changes watts while actual scattering
carries normalized Stokes state; no extra angular weight, albedo debit,
dichroism or response inferred from a scalar phase is introduced. The source
and phase-proposal bound remain immutable across clones and trajectory
refinement. Only phase entries with actual scattering in prepared physical
cells require polarized support. If heating activates an unsupported phase,
admission refuses before source histories and the native advance publishes no
partial thermal or optical state. Existing heat/error admission, unresolved
radiation, work budgets and exact native cancellation remain authoritative.

The original cell transport kernel follows all scattering continuations through
that prepared local field. True absorption supplies the original thermal cells,
including off-beam receivers; ballistic escape, diffuse escape and unresolved
radiation stay separate. The existing `ScatteringHeating` owner admits each
sampled report with explicit conservation, global and receiving-cell empirical
sampling-error allowances, and an unresolved-power allowance. Native conduction,
authored sources, latent enthalpy and events retain their owners. Accepted
physical durations integrate the actual changing radiation/heat with compensated
sums, and the native external-heat audit receives that heat once. Absorption-only
fields take the exact deterministic Beer path. Source-domain, native cancellation,
conservation and whole-advance source-evaluation/segment-budget refusals publish
neither thermal state nor partial optical energy.

Within a trajectory each accepted physical start time selects a fresh deterministic
counter stream from the original source seed/identity and replica seed. The key
excludes changing temperature/coefficient state, history count and worker identity.
Split advances with matching physical steps replay the same photons; paired
studies retain matching physical-time streams and photon-history prefixes.
This avoids recycling one set of packet draws throughout a thermally changing
trajectory. The most recent report and maximum conditional source standard
errors remain diagnostics in watts. They are **not** propagated final thermal
uncertainty, and changing-field energies do not claim the fixed-field linear
standard-error scaling of `ScatteringHeating`.

`advance_ensemble_refined` evaluates independent complete native trajectories,
with paired coarse/fine/finer time comparisons and separate final ensemble
sampling and H/2H photon-sensitivity gates. Both consecutive time comparisons
must actually increase native thermal substep counts; unchanged CFL/event-clipped
time grids cannot establish refinement. Paired temperature/enthalpy requests
below native coordinate resolution refuse without enlarging user tolerances.
Every thermal-cell enthalpy and actual caloric temperature, and every
radiation/heat energy channel, is assessed.
Energy comparisons use interval increments and compensated differences so old
fluence or arbitrary enthalpy datums cannot hide current discrepancies. The
result retains actual finest native trajectories and separate statistical
summaries; it does not manufacture a thermal state by averaging temperature or
enthalpy. All discarded trials share explicit total thermal-step, source-query,
traced-segment and retained-state ceilings. These empirical fixed-replica
assessments are sensitivity/standard-error diagnostics, not confidence sequences,
certified nonlinear bias bounds, rare-path guarantees or validation of supplied
material data. Selecting a resolution after looking at them confers no
optional-stopping guarantee.

Geometry, density, phase laws, incident power and spectral quadrature remain
prescribed. The incident-heating model does not infer missing material spectra,
reflecting/refracting interfaces, moving geometry, temperature-dependent phase
functions, fluorescence, chemistry or net thermal-radiation exchange. Native
integration tests cover off-beam heat, thermal activation of scattering, changing
upstream absorption and downstream heat, exact absorption-only reduction,
latent/event/split clocks, stream pairing, source-domain refusal and aggregate
work/cancellation rollback. Ensemble tests separately exercise deterministic
reduction, nonlinear stochastic trajectories and statistical/work refusal.

## No-claim boundaries (volumes)

- FrankenVDB tile-maxima majorants: no fvdb crate exists in-workspace;
  [`MajorantGrid`] builds per-block maxima from dense grids, and the
  per-tile-rate DDA traversal (rather than lookup thinning under a
  global bound) is the recorded successor alongside the FVDB wiring.
- Progressive live tiles with ledger artifact pinning (frame-consistent
  snapshots of evolving fields) — staged with the vessel flagship's
  render lane; the smoke tier renders a paused simulation's buffer.
- Specialized refractive volume-source connectors remain successors. Detached
  native liquid snapshots have the opt-in reconstruction described below. The shared path tracer
  implements borrowed absorption/emission/scattering in ambient and explicitly
  bound dielectric regions, with phase/source MIS; interface random walks
  supply paths that cross boundaries.
- The zero-copy claim at smoke tier is BORROW SEMANTICS (the API takes
  `&[f64]`; the battery binds a live `FreeSurface` mass buffer); the
  FrankenNumpy membrane view protocol is the fuller deliverable.

## Native free-surface optical geometry

The opt-in `free-surface` feature observes immutable `FreeSurface3` snapshots
through their actual unclamped interface mass/density ratios. Cell-center
samples of `0.5 - fill` pass through the original `fs-viz` marching-tetrahedra
owner, then the original renderer `TriMesh` BVH. Outward liquid/gas normals and
physical metre coordinates let the ordinary dielectric owner apply spectral
Snell/Fresnel transport, nested regions and Beer attenuation. The material-flow
entry point takes its cell width from the original SI scaling; neither entry
point infers refractive indices or absorption from numerical LBM density.

The caller supplies dense-sample and triangle budgets. Sampling, extraction and
facet admission poll cancellation; the existing BVH construction remains one
synchronous section bounded by the triangle limit. Empty interior sparse tiles
are gas. Nonzero unlocated carry, raw fills outside `[0,1]`, solids, wet domain
boundaries, unresolved crossings, degenerate facets and nonclosed or inconsistently
oriented edges refuse without changing accepted fluid state or returning a
partial mesh. Detached liquid requires an empty domain collar.

The default retains the half-fill geometry. The explicit
`volume_matched_free_surface_mesh` and material-flow counterpart instead match
the compensated original sum(fill) times SI cell volume. One global isovalue is
solved using the native mesher's allocation-free clipped-tetrahedron volume,
with an explicit pass ceiling and absolute/relative volume tolerance. Sampling
occurs once, in cell coordinates; only the accepted level creates a mesh.
After SI placement, translated compensated signed-tetrahedron volume checks the
actual closed mesh against the same tolerance before its BVH is constructed.
The returned report retains target and enclosed volumes, signed residual,
selected fill threshold and actual volume-query count. Unattainable constant
fill plateaus, exhausted passes and lost world-coordinate precision refuse.

The `component_volume_matched_free_surface_mesh` entry and its material-flow
counterpart assign a separate accepted-volume target to each disconnected
positive-fill support.
Connectivity uses the native tetrahedral edge graph on raw negative fill,
before half-minus-fill rounding; equal-zero nodes never join supports. IDs are
ordered by first original cell and describe only the current finite-resolution
PL support. They do not infer physical drops or track identities over time.
Each component uses one masked crop with a full gas-node halo and original
cell coordinates, so overlapping crop boxes cannot duplicate another support.
The same private allocation-free solver selects each level, each final crop is
meshed once, and one combined original TriMesh/BVH is built.

The original total absolute/relative tolerance is shared in proportion to
accepted component volume, with finite-positive allowances admitted before any
volume queries. Underflow refuses. Individual physical meshes and the actual
combined remapped closed mesh must satisfy their component and original total
allowances. Query and triangle ceilings are aggregate; explicit component-count
and summed-crop-node limits also bound support allocation and field work. Raw
positive fills remain labeled even when the half-minus-fill geometric model
cannot resolve them; unresolved supports refuse the entire result.

These are piecewise-linear visualization models: no PLIC, per-cell volume
conservation, inferred physical/temporal drop topology, subcell shape recovery,
smooth-normal reconstruction, global vertex-link manifold certification,
moving-exposure sampling or wet-wall closure is claimed. The numerical volume
tolerance is not a physical shape-error bound. Controls cover analytic
single-cell and detached unequal-drop volumes, masked overlapping crops,
single-support reduction, SI placement, immutable density-independent replay,
native Snell/Beer paths, and aggregate work/plateau/precision refusals.

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

`Aperture::try_annular(a,b)` declares a concentric circular obstruction with
outer radius `a` and inner radius `b`. Its inner/outer ratio must be representable
and at most 0.95; zero inner radius canonicalizes to the existing circle or
pinhole. The 0.95 ceiling is a numerical and sampling admission, not a physical
claim that narrower rings cannot exist. Geometric lens samples are uniform in
`r²=b²+u*(a²-b²)` with an independent uniform azimuth. A scaled `hypot` form
avoids squaring physical radii and preserves the inner endpoint. No rejection
loop or sampled origin inside the obstruction is substituted. The conventional
`LensSample::CENTER` is a unit-square midpoint, not a ray through this blocked
centre. `inner_radius_m`, `obscuration_ratio` and `is_annular` expose the shape.

`transmitting_area_m2` returns the actual collection area for each pupil:
zero for an ideal pinhole, `pi*a²` for a circle, `pi*(a-b)*(a+b)` for an annulus,
and the represented vertex area for a polygon. Nonzero areas that cannot be
represented as finite positive `f64` refuse. This is separate from the camera's
unit-power radiance average: no hidden aperture-area exposure multiplier is
introduced. Annular geometric rays also use the existing photometric path's
pupil-averaged radiance convention.

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

`Aperture::try_rectangular(width_m, height_m, rotation_rad)` admits a finite
rectangular or slit pupil with two positive full dimensions, finite positive
area and aspect ratio at most 1024. The circumscribed radius supplies the
existing defocus phase; the shorter dimension supplies scalar-wavelength
admission. Geometric rays uniformly sample the same rotated rectangle.
Focused propagation evaluates the signed separable sinc amplitude exactly
within the scalar Fraunhofer model, normalized by the actual transmitting
area. A product-Cauchy proposal in the two aperture-phase coordinates retains
the full anisotropic tails with its inverse-transpose Jacobian. No image-grid
FFT, slit-width blur fit or tail cutoff is introduced. Finite defocus reuses
the existing bounded polygon Fresnel boundary integral, with the rectangle's
actual moments at small phase. The same pupil enters source-depth layers,
dispersive lenses, coherent/Jones amplitudes and their finite-pixel/exposure
owners. Existing phase, numerical-aperture and visibility limits still apply;
zero-width or subwavelength slits and vector/high-NA diffraction are outside
this model. `tests/camera/rectangular.rs` compares signed lobes, physical
normalization and rotated Fresnel fields against independent separable
integrals; the pupil unit test integrates spectral power over complete tails.

`Aperture::try_rectangular_pair(width_m, height_m, center_separation_m,
rotation_rad)` admits two identical nonoverlapping rectangular openings on
one coherent pupil reference. The complete pupil span divided by the smaller
opening dimension may not exceed 1024. The circumscribed radius includes the
separation; transmitting area includes only the two openings. Geometric
samples select each opening with equal probability. Focused propagation is
the signed single-rectangle sinc product times the Young cosine factor,
normalized by total open area. Finite defocus restores both the translated
quadratic phase and shifted spatial frequency for each opening, then sums
complex amplitudes before squaring. Each component keeps its own radius for
numerical conditioning. Full-support mixture proposals retain the same
source-independent sampling contract across focused, defocused, layered and
dispersive paths. Coherent/Jones point, pixel and exposure owners share that
transfer. The touching-opening limit equals a single rectangle of twice the
width; overlap, zero-width openings and unresolved scalar wavelengths refuse.
Render continuity includes both component dimensions and separation. The
pair tests check physical samples, signed Young fringes, absolute Fresnel
phase, normalization and the touching limit; all existing scalar/paraxial,
phase and visibility restrictions continue to apply.

Finite rectangular gratings extend the same pupil owner through
`Aperture::try_rectangular_array(width_m, height_m, center_separation_m,
openings, rotation_rad)`. Counts 1..=256 are admitted; one opening
canonicalizes to the original rectangle and two share the original pair API.
The centre pitch must be at least the width, and total span divided by the
smaller opening dimension remains bounded by 1024. Geometry samples the actual
openings uniformly. Transmitting area is `N*width*height`; gaps contribute no
area or light.

At focus the signed field is the single-opening sinc envelope times
`sin(N*t)/sin(t)`, normalized by total open area, with `t=pi*pitch*theta/lambda`.
A division-free binary composition evaluates the finite sum in O(log N),
including the removable principal-order limits and destructive-interference
zeros. Finite defocus sums each translated original Fresnel field with its
quadratic and linear phases retained separately. Components keep their own
circumradius; DD phase reduction and accumulation preserve cancellation.
The bounded sum and mixture need O(N) work, with N <= 256. The full-support
Cauchy proposals retain all orders and tails; their peak importance weight
can grow with N, so adding openings may require more samples.

The array enters the original focused, defocused, source-depth, dispersive,
coherent/Jones and finite-detector paths through their shared pupil transfer.
Progressive rendering binds count, pitch and opening dimensions. This is
scalar thin-aperture diffraction within the existing phase, wavelength,
paraxial and visibility restrictions; it does not introduce a material,
blazed, subwavelength or high-NA grating model. Numerical tests compare
independently enumerated phasors, physical grating order widths, direct
opening integrals, complete spectral power and the touching-array limit.
The diffraction-array theorem follows
[TU Delft, Scalar Diffraction Optics, equations 403–413](https://interactivetextbooks.tudelft.nl/interactive-optics/content/Chap6_Diffraction/DiffractiveOptics_2022_01Clean.html).

`PhysicalCamera::with_in_focus_diffraction(FocusedDiffraction)` opts the
cinematic spectral tracer into scalar, paraxial, incoherent pupil imaging.
It requires a nonzero circular, annular, rectangular or regular 3–64-blade aperture,
`lambda/(2*inradius) <= 0.01` for circles/polygons or
`lambda/(a-b) <= 0.01` for an annulus, and an explicit
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

A concentric annulus retains coherent interference between its two curved
edges. With `e=b/a`, its focused amplitude is
`B(q)=2*(J1(q)-e*J1(e*q))/q` and its unit-integral dimensionless density is
`abs(B)²/(4*pi*(1-e²))`. Its focused peak is the actual transmitting area
divided by wavelength squared. Removing the centre changes the rings and the
central peak; it does not merely darken an otherwise unchanged Airy function.
Small image radii use a factored series. The finite-defocus composition below
uses the same pupil and normalization. Its complete-support proposal is an
equal mixture of isotropic Cauchy densities with scales 2 and
`max(2,2*abs(phi_reference),1/(1-e))`. This keeps the narrow central lobe while
covering the wider ring groups of an obscured pupil. Every draw evaluates the
full mixture PDF for each wavelength; both components retain inverse-cubic
area-density tails without an angular cutoff or rejection loop.

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

### Finite-defocus pupil imaging

`PhysicalCamera::with_defocused_diffraction(DefocusedDiffraction)` selects
scalar propagation of a circular, annular or regular 3–64-blade pupil at a declared object plane. The
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
The camera basis and axial eye coordinate stay fixed throughout the owning
shot; lateral translation parallel to the object plane and focus racks
remain supported. Unknown charts, animated geometry, prepared scene histories,
primary media, calibrated photometric fixtures and specular primary surfaces
refuse. The mesh scan occurs once at render admission, and actual sampled
hits retain the depth check. Camera-only consumers must establish equivalent
clear-pupil and pupil-independent-radiance conditions themselves: a first-hit
depth check cannot detect an obstacle that clips only part of a pupil.
For an annular aperture, the declared concentric obstruction is part of the
propagated pupil; these conditions forbid additional scene-dependent clipping.

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

The circular numerical owner evaluates a stack-only 40-term Bauer/Zernike expansion,
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

Annular finite defocus splits that same complex pupil integral at the inner
radius. With `e=b/a` and `delta=(1-e)*(1+e)`, its propagated amplitude is
`B(q,phi)=A(q,phi)-e²*A(e*q,e²*phi)`, and its normalized density is
`abs(B)²/(4*pi*delta)`. Both disk amplitudes restore their different common
`exp(i*phi/2)` factors before subtraction. The signed phase is retained;
negative phase conjugates the complete amplitude. Subtracting disk intensities
or their phase-stripped complex series would give a different optical model.
The on-axis limit is evaluated directly as
`delta*sinc(phi*delta/2)²/(4*pi)`, including its nulls. The existing outer-rim
phase cap of 16 applies to both constituent disk evaluations; a small annular
phase spread does not relax that cap. The clear-area fraction is at least
0.0975 at the admitted maximum ratio, bounding the cancellation associated
with subtracting nearly coincident disks. This uses two bounded circular
expansions, with no additional quadrature, pupil raster or radius-dependent
work. A zero-width ring is refused: it has no transmitting area and its
formal `J0²` limit is not a normalizable two-dimensional PSF.

Bladed defocus instead propagates the actual polygon by exact Fresnel boundary
identities. Its normalized density is `abs(F(q,phi))²/(4*pi²*area)`, where
`F=integral_pupil exp(i*phi*abs(r)²-i*q.r) dA`; the circumradius remains the
phase/coordinate scale and the actual polygon area remains the normalization.
Near the stationary point, the analytically regularized `exprel` edge formula
retains its finite limit even on an edge or vertex. In the exterior region,
`F=i*integral_boundary [(q-2*phi*r).n/abs(q-2*phi*r)²]*exp(iS) ds` avoids
large completed-square phases. At very small image radius and defocus, exact
polygon second/fourth radial moments give an intensity expansion with a
fourth-order remainder below `1.07e-19`. Zero phase calls the existing focused
polygon owner exactly.

Bounded tangential frequencies use positive 16- or 32-point Gauss rules,
choosing the cheaper composite rule from an explicit derivative-remainder
bound. Positive recurrences bound the quadratic exponential and the separated
rational edge factor. Subdivision is capped at 64 segments per edge; unresolved
work refuses rather than publishing a fallback blur. For `abs(q)>=8*abs(phi)+32`
(with actual-vertex roundoff clearance) and tangential phase at least 96,
a fixed 32-term endpoint integration-by-parts expansion replaces oscillatory
quadrature. Its derivative remainder is bounded separately and must meet the
same `2e-13` normalized edge budget. These analytic truncation bounds exclude
floating-point evaluation error. Small-phase subtraction, compensated sums and
separate bounded defocus/large linear phases protect the removable limits and
far tails. This adds no runtime dependency, pupil raster, Fourier grid,
rejection loop, per-ray allocation or tail cutoff.

The existing blade-aligned product-Cauchy mixture widens with defocus to cover
the geometric pupil footprint and its directional diffraction ridges. Each
packet lane evaluates the actual signed phase and its own PSF/proposal weight.
For odd blade counts, changing defocus sign mirrors the intensity pattern:
`I(q,-phi)=I(-q,phi)`; using `abs(phi)` would be physically wrong. Circular, annular
and even-blade pupils also have central symmetry. The boundary identities
follow Barnett, *arXiv:2010.05978*, Appendix A, equations 19 and 22,
<https://arxiv.org/pdf/2010.05978>, and Huang, Christian and McDonald,
*JOSA A* 23, 2768 (2006), equations 14–16,
<https://www.salfordphysics.com/gsmcdonald/32_JOSAA_v23_p2768_2006.pdf>.

Pure kernel tests check the axial sinc law, exact focused limit, sign
symmetry, independent adaptive-integral landmarks across all numerical
branches and into the far tail, and three-wavelength unit power with 65,536
deterministic radial draws (0.4% quadrature band). The last representable
open-unit draw retains finite importance weights. Camera tests cover physical
depth and SI scaling, chromatic axial nulls, lens-sample independence,
shortest-wavelength phase refusal, infinity and focus-rack retention. The
polygon kernel adds independent separable-square Fresnel and direct triangle
area-integral landmarks through stationary edge/vertex points, branch
transitions and far tails; signed odd-blade patterns; rotation; focused-bit
preservation; extreme full-support draws; and three-wavelength unit power with
16,384 Hammersley samples for triangular and heptagonal pupils (1% integration
band). Annular tests compare direct adaptive integrals over the opening
against coherent disk composition, including recurrence transitions, both
phase signs, 0.95 obscuration and image radius 1000. They check axial nulls,
exact clear-disk limits, finite extreme draws and unit power at three
wavelengths using 65,536 radial mixture draws (0.6% integration band).
Camera tests bind these kernels to physical aperture area, wavelength
scaling, geometric annulus support and area moments, radial-width admission,
rotated rays and separate packet weights. The
phase cap bounds numerical work; it does not certify the paraxial, high-NA,
wide-angle, lens-aberration or measured-lens error. General multi-depth
occlusion, virtual images, coherent scene transport
and vector diffraction remain outside this mode.

### Complex pupil transfer for authored coherent fields

`PhysicalCamera::in_focus_diffraction_amplitude` and
`defocused_diffraction_amplitude` expose the same circular, annular and bladed
pupils before the intensity square. With unit-circumradius pupil P, area A,
`s=lambda_medium/(2*pi*a)` and `q=tangent_offset/s`, the transfer is
`integral_P exp(i*phi*|r|²-i*q.r) dA / (2*pi*sqrt(A)*s)`.
Its squared magnitude matches the original angular density, and its focused
on-axis value is positive `sqrt(physical pupil area)/lambda_medium`.

Circular and annular transfers retain their signed piston phase. Polygon
transfers restore the stationary formula's `exp(-i*|q|²/(4*phi))` and the
exterior formula's common `i`. Near DC, a cubic complex exponential expansion
retains odd moments as well as quadratic phase; bounded 4x4 Duffy/Gauss area
rules integrate that polynomial exactly. Its Taylor remainder per unit area
is below `5e-18` in the focused branch and `7e-21` in the tiny-defocus branch,
excluding floating-point error. The original density arithmetic is retained.
A prepared internal owner reuses wavelength admission and pupil coefficients
across observations, borrows original polygon vertices and allocates nothing.

These queries share the original medium, material-band, clear-pupil, paraxial
and 16-radian phase limits. Their reference is the pupil plane; external source
propagation phases and source coherence must be authored separately. They do
not change the radiance renderer, support layered source planes, or claim
polarized vector diffraction. G1/G3 tests compare absolute complex transfers
with independent square Fresnel and triangle area integrals, both sides of the
stationary/exterior transition, endpoint asymptotics, annular axial phase,
small odd moments, density parity, conjugation, physical medium wavelength and
original material-band refusal.

### Coherent and partially coherent pupil illumination

`PhysicalCamera::prepare_coherent_diffraction` composes the original normalized
complex pupil transfer with authored `CoherentPupilWave` fields. Every wave
shares one entrance-pupil phase reference, scalar polarization, vacuum
wavelength and declared source-plane curvature. A pupil tilt
`a_j exp(+i k theta_j.u)` contributes `a_j h(theta-theta_j)` at the observation;
there is no separately shifted sensor chirp or inferred object-to-pupil phase.
Coefficients are sqrt(W), normalized so their squared magnitudes are isolated
transmitted powers through the original pupil. Results are power per tangent
coordinate area, not detector irradiance per square metre.

The deterministic mean and each explicit covariance factor are reduced as
complex fields before squaring. Their powers then add as
`|H mean|² + sum_r |H factor_r|²`, with no cross-factor coherence and no final
power renormalization. Factor weights belong in amplitudes as square roots.
This positive factorization permits partial coherence without constructing a
dense source-pair matrix. DD complex products and sums preserve weak residuals
under destructive interference within the supported floating-point range.

Preparation borrows immutable source banks and caches the original pupil and
wavelength admission. Each point query has O(total authored waves) work,
constant scratch and no allocation. Explicit wave/factor limits apply before
preparation; cancellation is checked for every wave and factor and before
returning a complete sample. Source and observation tangent magnitudes must
not exceed 0.1. This is a scalar paraxial model boundary, not an error bound.
The original clear-pupil, material-domain, physical-medium and finite-defocus
limits remain. Layered source planes, scene-derived coherence, vector/polarized
imaging and finite-shutter beat integration are outside this query; the
existing radiance renderer retains its original semantics.

`PreparedCoherentDiffraction::integrate_pixel` collects mean, fluctuation and
total power in watts over an admitted `CoherentDiffractionPixel` rectangle.
Bounds use the original tangent coordinates. For a perpendicular physical
image plane at distance d, lateral bounds map as x/d and the irradiance and
area Jacobians cancel; no solid-angle cosine is inserted. It integrates each
source bank's complete field square, then adds independent covariance powers,
so a dark centre never replaces bright fringes elsewhere in the pixel.

The original medium wavelength and pupil radius set the diffraction scale s.
The pupil's intensity Fourier support lies within [-2,2] on each normalized
axis, including defocus and all tilted-source cross terms. Positive 4-by-4
Gauss quadrature uses subcells with `2*width/s` no greater than the caller's
`maximum_phase_span_rad` in (0,1]. This bandwidth-guided resolution is not an
error certificate; reducing the span assesses quadrature error independently
of the scalar/paraxial model. Counts for all samples and original pupil-wave
evaluations are checked before the first field query and returned with the
result. The fixed traversal uses DD accumulation and constant scratch with no
allocation, polls at every sample and wave, and leaves no partial publication
on cancellation. Empty illumination returns zero without field evaluation.

`tests/coherent_diffraction.rs` checks the physical four/two/zero interference
limits, signed Airy lobes, weak cancellation residuals, common-phase and
covariance-basis invariance, independent complex pupil integration, bounded
admission, cancellation/reuse and numerical overflow.
Finite-pixel regressions additionally compare an independent separable square-
pupil integral, retain exact coherent nulls and independent-factor powers,
check pixel partition additivity and weak residuals, and refuse work before
evaluation. The runnable
`examples/coherent_pupil.rs` emits a CSV comparison of coherent, partially
coherent and incoherent illumination with the same pupil and isolated powers.

### Two-arm detection of mutual spatial coherence

`PreparedCoherentDiffraction::observe_interferometer` combines the original
monochromatic scalar camera field at two sensor coordinates through passive
arm power fractions and a balanced lossless combiner. Each mean/covariance
bank propagates separately to both coordinates; its two fields interfere
before independent bank powers add. Coincident samples cannot duplicate
source power: the two input fractions must sum to at most one.

The receiver restores the relative physical sensor Fresnel chirp and accepts
an explicit optical path difference (including its phase index) plus an
authored interface/combiner phase. The original physical sensor-conjugate,
medium-wavelength, pupil and paraxial admissions remain. Both output powers,
arm powers, complex mutual coherence and fringe visibility are returned as
point densities per tangent-coordinate area. Equal-area unit-magnification
arm relays are assumed; resolved interferometer optics, finite detector
integration, polarization and spectral exposure are outside this receiver.

Complex DD fields retain destructive output differences before squaring.
The lossless output powers sum to the admitted arm powers; unequal arm
powers and partial spatial coherence reduce fringe visibility. Original
pupil-wave work is preflighted, returned and allocation-free; cancellation
leaves the preparation reusable. Focused regressions cover path-length
fringes, output closure, weak dark ports, partial coherence, complex defocus
with sensor chirp, and passivity/sensor/work/cancellation admission.

### Finite relay pixels of the two-arm receiver

`PreparedCoherentDiffraction::integrate_interferometer_pixel` collects both
ports over one rectangular relay-offset footprint, translated about each
receiver coordinate. Both relays use the same tangent-area measure and unit
magnification; full translated corners retain the original paraxial admission.
The original positive spatial rule integrates six DD moments: both arm powers,
both port powers and complex mutual coherence. Fields interfere before these
moments, and independent covariance banks add afterward. Collected values are
watts, not point densities; the integrated complex coherence determines the
visibility of an achromatic phase scan across the complete detector.

Resolution includes the original pupil intensity bandwidth plus the relative
sensor-chirp gradient in each axis. This prevents a large arm separation or
sensor distance from hiding unresolved spatial fringes. Active-arm pupil work
and all quadrature samples are admitted before field evaluation. Repeated
queries allocate nothing, preserve weak dark-port signals and leave preparation
reusable after refusal/cancellation. Refinement assesses quadrature convergence,
not an optical error certificate. Independent square-pupil/Simpson references,
spatial contrast loss, partition/factor additivity, work refusal and finite
large-power cases supplement the original point-interferometer regressions.

### Polarized interferometer arms and full output coherency

`PreparedPolarizedCoherentDiffraction::observe_interferometer` carries both
original transverse Jones fields through independent passive arm optics and
the shared two-arm receiver. `CoherentJonesArm` caches at most eight ordered,
uniform native `CoherentJonesPupilLayer` elements, preserving their principal
transmissions, retardance and axes. Spatial gradients refuse; these uniform
arm elements act after the pupil-to-sensor propagation. Their passive losses
remain in the returned arm powers and are never renormalized away.

Both combiner outputs retain full coherency and native Stokes, together with
complex trace mutual coherence, total-power fringe visibility and pupil work.
Each source bank interferes before the outer product; independent covariance
banks add afterward. Orthogonal arm fields can change output polarization
without modulating total intensity. An analyzer can reveal existing coherence
but cannot create it between independent source banks. DD field differences
precede output powers, retaining aligned weak populations and dark ports;
arbitrary rotated cancellations retain the existing finite-precision limits.
The shared scalar sensor/phase/passivity admissions remain authoritative.
Queries allocate nothing, and ordered arm work is bounded by the layer cap and
original source-bank limits. Focused regressions cover polarization/analyzer
laws, circular phase, independent factors, scalar parity, weak populations,
finite bright outputs and refusal/cancellation reuse.

### Stationary finite-linewidth interferometry

`InterferometerTemporalCoherence` supplies an explicit Lorentzian or Gaussian
power-spectrum FWHM in Hz and a separately authored complete relative group
delay in seconds. Scalar and Jones `observe_interferometer_with_linewidth`
reuse the original source banks, carrier phases, passive splits and arm optics.
The normalized field correlations are `exp(-pi*FWHM*abs(delay))` and
`exp(-(pi*FWHM*delay)^2/(4*ln(2)))`. Both ports use positive mixtures of their
coherent and independent-arm outer products; cross-arm coherency changes but
arm energy does not. `expm1` retains linewidth-induced leakage at a dark port
even when the rounded correlation is one. The scalar finite-pixel counterpart
evaluates an optional authored linear delay gradient inside the original positive
quadrature and admits resolution/work for the temporal envelope too.

This is a stationary separable narrowband ensemble mean: pupil transfer and
Jones optics are frozen at the original wavelength, and the authored delay
contains all relative propagation/optics. Phase OPD is not substituted for
group delay. Admission caps FWHM at 1% of carrier frequency; this is not an
accuracy certificate. Higher-order dispersion, broad chromatic response,
source turn-on and finite-shutter fluctuations remain outside this model.
The existing explicitly resolved line-spectrum exposure path is unchanged.
Focused G1/G3/G4 tests compare independent spectral Fourier integrals,
finite-pixel square-pupil integrals, weak dark leakage, polarization/spatial
independence, power closure and refusal/cancellation reuse.

### Spectral interferometry on the original finite shutter

`PreparedCoherentDiffractionExposure::observe_interferometer` recombines both
arms at every original frequency before the analytic temporal Gram contraction.
Equal-frequency entries coalesce, and each covariance column retains its shared
identity across frequencies and arms. Results include shutter-averaged arm/port
power densities, complex mutual coherence, visibility and exposure-energy
densities. Visibility refers to an achromatic combiner-phase scan. The authored
OPD is one nondispersive optical length; additional arm dispersion is not inferred.

Propagation retains authored Hz and the original phase index. Frequency-length
products reduce modulo exact c before division and conversion to radians;
relative optical frequencies form before propagation. This preserves close-tone
beats across large source epochs and optical paths without inverting a rounded
pupil wavelength. Spatial pupil admission still uses its original wavelength.
No optical carrier or shutter raster is introduced. Fixed observation scratch,
active-arm pupil limits and explicit frequency-contraction limits bound work:
B*N*(4*N+1) for two arms, or the original B*N*(N+1)/2 for one active arm.
Closed/empty observations perform no pupil or pair work. Failure and cancellation
leave preparation reusable. Six focused groups cover scalar equivalence,
independent analytic beat washout, covariance identity, equal-line cancellation,
weak dark ports, long-path/source-clock phase and exact work/admission behavior.

### Spectral interferograms over finite relay pixels

`PreparedCoherentDiffractionExposure::integrate_interferometer_pixel` collects
both scalar ports, arm powers and complex mutual coherence over the original
finite shutter and a common rectangular relay offset. `CoherentInterferometerPixelExposure`
returns average watts, shutter-integrated joules and exact spatial/pupil/pair
work. Source banks recombine as DD fields before analytic temporal contraction;
the native positive spatial rule integrates those unrounded moments. Visibility
is derived from the integrated complex coherence, preserving spatial and temporal
fringe washout.

The shared work admission covers every translated corner, the shortest pupil
scale and all active cross-frequency/cross-arm sensor-phase gradients. A range
of endpoint gradients and extreme optical frequencies bounds those pairings in
O(N) preprocessing. Common gradient subtraction and center-plus-offset DD
propagation retain close-tone beats even for offsets smaller than one binary64
ULP of an arm center. All sample, pupil and contraction limits precede field
work; closed/empty receivers do none. Queries allocate nothing and remain reusable
after cancellation. Equal-area unit-magnification relays, nondispersive arm OPD,
paraxial propagation and cycle-averaged detection remain the model boundaries;
phase-span refinement is resolution control, not an error certificate.

Five focused G1/G3/G4 groups cover independent square-pupil/shutter integrals,
monochromatic parity and weak ports, covariance/partition additivity, large
path/epoch invariance, sub-ULP offset beats and exact work/cancellation admission.

### Polarized spectral interference at both output ports

`PreparedCoherentJonesDiffractionExposure::observe_interferometer` retains both
Jones components through passive arm optics and recombines each original
frequency before the native full-coherency shutter contraction. Both output
coherencies and native Stokes, arm powers, complex trace mutual coherence,
visibility and energy densities remain available. Original covariance-column
identities span frequencies, polarization and arms; equal-frequency source
entries share one pupil response and combine before products. Independent banks
add afterward, so an analyzer reveals existing coherence without creating it.

Native patterned chromatic pupils, authored-Hz propagation, strict Jones-arm
products and physical coherency admission are reused. Uniform arm sequences are
explicitly frequency-independent. DD fields preserve weak aligned populations
and nearly dark ports; no arbitrary-null relative certificate is implied.
Observation allocates nothing and borrows the cached arm sequences internally.
For N frequencies and B source banks, exact contraction work is
B*2*N*(6*N+1) with two active arms, B*N*(3*N+1) with one, and zero when closed
or empty. Explicit receiver limits preflight work independently of preparation's
single-detector limit. Cancellation/refusal preserves the preparation.

Seven focused groups cover monochromatic parity, analytic polarization beats
and analyzers, source independence, equal-line and weak-channel cancellation,
patterned-pupil work, exact limits/cancellation and finite bright-port range.

### Joint polarized spectral interferograms on finite pixels

`PreparedCoherentJonesDiffractionExposure::integrate_interferometer_pixel`
combines polarization, spectral exposure and physical detector area in one
observation. `PolarizedCoherentInterferometerPixelExposure` returns both output
coherencies/Stokes in watts and joules, collected arm powers/energies, complex
mutual coherence, visibility and exact sample/pupil/contraction counts. A mixed
spatial/temporal observation is not assigned a single Jones field.

The receiver integrates twelve unrounded DD moments with the native positive
area rule. Source-bank fields recombine before temporal coherency contraction;
complex mutual coherence integrates before visibility is formed. It reuses the
scalar spectral receiver's complete footprint and cross-arm/frequency phase
bound, including DD relay offsets. All work is preflighted, scratch is fixed,
and each sample borrows the uniform arm sequences without allocating or copying
their full cache. Wholly underflowed nonzero diagonal populations are refused.
The original equal-area relays, paraxial cycle-averaged detector, frequency-specific
pupil responses and frequency-independent uniform arm matrices define scope.

Four focused regression groups cover native patterned-pupil pixel parity,
independent separable square-pupil and analytic-shutter polarization beats,
Mueller/analyzer and pixel-partition laws, weak dark channels, energy closure,
exact budgets, closed receivers and deterministic cancellation/reuse. The
combined original-source optical harness passed all 76 selected tests; native
full-workspace and tracer-feature integration were not exercised in that run.

### Prepared frequency-specific interferometer arm optics

`CoherentJonesInterferometerResponse` supplies a pair of passive Jones-arm
sequences for one exact original source frequency. The exposure owner's
`prepare_interferometer` binds a receiver and a complete response table into
`PreparedCoherentJonesInterferometer`, whose `observe` and `integrate_pixel`
methods use the same full-coherency, shutter and spatial moment evaluators as
the uniform-arm APIs. Each frequency retains its own transmission, retardance,
principal axes and common phase before interference. Passive losses remain
physical; no missing response, interpolation or material dispersion is inferred.

Preparation accepts arbitrary response order but requires exactly one response
per distinct authored frequency, including shared equal-frequency source entries.
Missing, duplicate and foreign identities refuse. Cache payload bytes are admitted
before allocation as N*size_of::<[CoherentJonesArm;2]>(), excluding Vec/allocator
metadata. Existing 64-frequency and eight-layer-per-arm limits bound the cache.
Bounded matching occurs once; observation uses direct group indexing without
allocation or frequency search. The immutable source borrow prevents changing
its pupil/source owner while the receiver is live. Failure/cancellation publishes
no receiver and preserves the source.

Five focused G1/G3/G4 groups cover bitwise uniform-arm parity including patterned
pupils, independent chromatic Jones/shutter references, retained common spectral
phase versus independent sources, response/source order and equal-frequency
identity, exact cache bounds and cancellation/reuse. The combined original-source
optical harness passed all 81 selected tests. Native full-workspace and tracer
feature integration remain unverified by that harness.

### Finite-exposure coherent camera pixels

`PhysicalCamera::prepare_coherent_diffraction_exposure` accepts an authored
positive-frequency `CoherentDiffractionTone` spectrum, an explicit sensor axial
distance and a uniform `CoherentDiffractionShutter`. Frequencies use a shared
pupil-centre clock and `exp(-i*2*pi*f*t)`. Each line independently reuses the
original pupil at vacuum wavelength c/f, including its signed defocus, physical
medium and material-band admission. Equal frequencies coalesce as DD fields
before any power or cross moment. Covariance column r names the same source
coefficient at every frequency; absent contributions use empty wave slices.
Distinct columns remain uncorrelated. Matching wavelengths or nearby spectral
lines never establish source coherence automatically.

The measured quantity is the squared analytic scalar optical field, preserving
the original sqrt(W) amplitude and monochromatic power normalization. It is
cycle-averaged photodetection, not subcycle real-field measurement; optical
sum-frequency terms and detector/material nonlinear response are excluded.
The existing `fs-fft` difference-frequency kernel integrates every cross term
analytically over the actual shutter. Its real-peak normalization is explicitly
removed for analytic-field power. Computation does not step through optical
carrier cycles or beat periods, even for long exposures.

At the physical sensor coordinate x=d*theta, the field includes the common
Fresnel phase `exp(i*k_f*d*(1+|theta|²/2))`. Longitudinal propagation delay and
quadratic detector chirp cancel for one frequency but cannot be omitted across
frequencies. They use the original medium phase index and retained DD frequency
differences. Source-to-pupil and additional material/lens piston phases remain
authored. A physical focal-length camera enforces its original sensor conjugate;
an FOV-only camera explicitly declares its otherwise unknown sensor distance.
This is the same scalar paraxial model, with source/observation tangent norm at
most 0.1, and does not enable coherent scene transport or vector diffraction.

`PreparedCoherentDiffractionExposure::observe` returns shutter-average power
density. `integrate_pixel` returns average collected power and integrated joules
using the original positive spatial Gauss rule. Resolution includes both the
shortest physical wavelength's pupil bandwidth and the frequency-dependent
sensor-chirp gradient across that pixel. It is a bandwidth-guided resolution,
not a quadrature or optical error certificate. Original pupil evaluations,
spatial sample counts and total frequency-pair contractions are preflighted.
The spectrum has at most 64 authored tones; all source/factor counts are also
caller-bounded. Cached temporal kernels and pupil transfers make observations
allocation-free with bounded scratch and cancellation checks. Negative or
unrepresentable contracted powers are refused, never silently clipped; DD
moments do not certify relative accuracy at arbitrary interference nulls.

Focused tests compare analytic beat visibility, shutter-clock and sensor-delay
phase, independent wavelength/defocus transfers, a separable square-pupil pixel
integral with strong detector chirp, a positive sub-ulp destructive-beat energy,
shared versus independent spectral covariance columns, monochromatic power and
energy normalization, equal-frequency weak residuals, spatial refinement,
physical-conjugate refusal, work preflight and cancellation/reuse.

### Polarization-resolved coherent camera detection

`PhysicalCamera::prepare_coherent_jones_diffraction_exposure` propagates authored
`CoherentJonesPupilWave` spectra through the original wavelength-specific pupil,
defocus, physical medium and finite-shutter optics. Each vector wave has two
complex sqrt(W) amplitudes in the same orthonormal paraxial pupil x/y basis;
their squared magnitudes sum to its isolated clear-pupil power, before any
optional Jones-pupil layers described below. One tone is
the monochromatic case. Frequencies retain the common pupil-centre clock,
negative-time analytic-field convention and original chromatic sensor phase.

Jones covariance column r denotes the same random coefficient across both
polarizations and all frequencies. Different columns remain uncorrelated.
Independent horizontal and vertical fluctuations therefore need separate
columns; a deterministic Jones vector is not an unpolarized source. Equal
frequencies coalesce as DD vector fields before detector projection and power.
The camera never creates phase coherence from matching polarization, similar
frequencies, Stokes radiance or scene geometry.

`CoherentJonesDetector::POLARIZATION_INSENSITIVE` measures |Ex|² + |Ey|².
Orthogonal components cannot cancel. `CoherentJonesDetector::analyzer(mode)`
normalizes a finite nonzero Jones direction and measures |a† E|². The ideal
analyzer is frequency-independent and selects one mode; its rejected
orthogonal output is not another input to this detector. Linear, circular and
elliptical modes are specified by their complex components, with no inferred
viewing-direction handedness. Mode magnitude is deliberately not attenuation
or gain. DD projection happens before the analytic temporal contraction,
preserving weak analyzer residuals and permitting interference between
coherent orthogonal inputs after projection. An analyzer never establishes
coherence between independent source columns.

Point observations return shutter-average tangent-coordinate power density;
finite pixels return average watts and integrated joules. They reuse the
original positive spatial quadrature, physical pupil bandwidth, cross-frequency
sensor-chirp bound and analytic shutter kernels. Each unmasked vector source
evaluates its pupil once per spatial sample for both components; patterned
layers multiply this count by their retained Fourier modes. Source budgets count
vector waves. Preparation preflights both polarization components so a cached
query supports either detector; each pixel charges its actual one or two
analyzed components. Observation uses bounded scratch, no allocations and
bounded cancellation checks. Errors publish no partial measurement, and the
same source/pupil preparation is reusable with another analyzer or after
cancellation. Original numerical refusals and resolution limits remain.

The default is two-component paraxial diffraction through an isotropic scalar
pupil; the optional passive patterned layers below extend its Jones response.
Longitudinal/high-NA vector diffraction, nonparallel ray-frame transport and
coherent scene rendering remain outside this query. An
arbitrary finite-pixel polarization state is not inferred from its integrated
power. Scalar/paraxial and quadrature errors are not certified by this API.

The coherent-camera tests verify orthogonal and matched interference, Malus's
law, complex circular-mode conjugation, shared versus independent spectral
coherence, wavelength-specific defocused beats, independent square-pupil pixel
integrals, complementary analyzer power, exposure energy, DD weak residuals,
basis rotation, exact logical work bounds and cancellation/reuse.

### Passive patterned Jones pupils for spectral camera interference

`PreparedCoherentJonesDiffractionExposure::set_pupil_responses` installs one
explicit `CoherentJonesPupilResponse` per distinct source frequency. Repeated
tones at that frequency share the same pupil. Each response is an ordered
stack of passive thin `CoherentJonesPupilLayer` elements on the same pupil
plane: `J(u)=R(theta) diag(t0,t1) R(theta)^T`, with
`theta=theta0+g.u` and `tj=sqrt(Tj) exp(i phi_j)`, `0<=Tj<=1`.
Uniform layers represent linear retarders, polarizers and diattenuators;
a rotating axis represents a geometric-phase polarization grating. Authored
frequency-specific loss and both principal phases retain the common phase
needed for finite-exposure beats. No spectral interpolation or missing
material response is inferred. Input amplitudes remain normalized to the
clear pupil, so passive loss is never renormalized away.

Each layer has at most three exact spatial Fourier modes. Coplanar composition
preserves matrix order, convolves those modes and merges identical gradients
in DD; only exact-zero matrices are removed. Each mode applies its Jones
matrix before coherent field addition and queries the ORIGINAL chromatic
pupil at `theta_observation-theta_source-q/k_medium`. Circular, annular and
polygon geometry, signed defocus, medium wavelength and common sensor phase
remain owned by the original transfer. Deterministic means and shared spectral
covariance columns retain their existing coherent/incoherent semantics. The
shared field owner supplies detector power, analyzer projection and full
spectral coherency/Stokes through the same pupil; sensor phase stays at the
actual observation coordinate, not at a shifted Fourier-mode coordinate.
The original analytic shutter and positive finite-pixel integration remain
unchanged. Multiplying the physical pupil by a Jones layer does not enlarge
its spatial support, so the original intensity-bandwidth resolution remains
applicable; it is still not a quadrature certificate.

`CoherentJonesPupilLimits` bounds total authored layers, completed frequency
modes plus the active intermediate expansion, and original pupil evaluations
per point. Each frequency has an additional hard cap of 64 intermediate
modes; expansion refuses instead of truncating coefficients. Pixel work
includes every compiled mode for every wave/factor and is preflighted before
evaluation. Shifted incident directions retain the original paraxial admission.
Wholly underflowed products of nonzero complex fields/transfers or layer
coefficients refuse rather than becoming false extinction. Preparation is
transactional under data, work, range or cancellation failure, while
point/pixel queries allocate nothing and retain bounded cancellation checks.
Empty layer lists restore the clear pupil. This model excludes propagation
between layers, multiple reflections, arbitrary pupil textures and high-NA
longitudinal fields; it does not extend the separate monochromatic Stokes owner.

Five focused regression groups check independent physical square-aperture
Jones integration with a 32/64/128 refinement ladder and both defocus signs,
noncommuting layers, weak aligned transmission, chromatic finite-exposure
phase, polarization-grating steering, finite-pixel energy, and transactional
bounds/cancellation. The focused actual-source coherent camera suite passes
36 tests, including the existing scalar, spectral, monochromatic coherency
and full spectral coherency tests. Patterned coherency traces and analyzer
contractions agree with detector power, and both pixel APIs account for modes. This bounded native execution is not a full workspace/tracer build or
experimental validation of authored material responses.

### Source-depth-resolved diffraction across separated planes

`PhysicalCamera::with_layered_diffraction(LayeredDiffraction)` admits one to
sixteen explicitly declared positive axial source depths, strictly ordered
nearest to farthest, with the same `[0,0.1]` residual pupil-rim phase budget.
It replaces the focused or single-plane policy. Circular, annular and regular
bladed pupils reuse their existing propagated complex-pupil PSFs, physical
area normalization, scalar feature-width gates and 16-radian outer phase cap.
Depth varies between visible sources, rather than forcing the complete image
to share one declared object plane.

`sample_layered_diffraction_ray` returns a `LayeredCameraRay` that borrows the
exact camera and retains its wavelength packet and complete proposal. It
exposes a source ray, with no provisional spectral weights. The caller traces
that ray and consumes `resolve_source` with its actual first-hit point, or
`None` for a black miss. The latter returns zero weights. The older one-stage
spectral API refuses this policy with `DiffractionRequiresSourceDepth`.
The spectral tracer caches that one primary intersection and supplies it to
the normal first-vertex material, radiance and AOV paths; it does not perform
a duplicate visibility traversal or bypass volume/fixture refusals.

The source direction is sampled independently of its unknown depth. For
reference coordinates `q` and wavelength ratio `r=lambda_ref/lambda`, each
lane receives `r²*p(r*q,phi(z_actual,lambda))/proposal(q)`. Its signed phase
uses the sampled source's actual axial depth, after admitting that depth
against the declared set. It never chooses a kernel from the unshifted
destination pixel's hit. Each source PSF has unit integrated power, so the
image is the incoherent sum of its independently blurred source contributions.
`layered_diffraction_density` exposes the same physical angular-raster density
for an admitted source depth and wavelength. There is no image postprocessor,
intermediate layer image, PSF crop or wavelength-independent blur.

The proposal is a full-support equal mixture of narrow and broad components.
The broad component covers the largest declared phase magnitude at the current
focus and reference wavelength. That maximum occurs at one of the ordered
depth endpoints because phase is affine in inverse depth; choosing it takes
two phase evaluations regardless of the number of layers. Circular and annular
pupils retain the radial Cauchy mixture, including annular ring spreading;
polygons retain the blade-aligned product-Cauchy mixture at narrow and broad
scales. Every draw evaluates the full mixture PDF and every source depth uses
that same PDF. The numerical kernels retain complete tails and bounded work.

The renderer preflights a conservative full-pupil geometric visibility proof.
Every forward primitive must be a static triangle mesh or fixed rigid mesh
instance, entirely in one declared frontoparallel plane to coordinate
roundoff, with an opaque diffuse material. All its vertices must belong to
the same layer; a mesh connecting layers is refused. Other geometry must be
wholly behind the pupil plane. Finite layers require a black background;
environments, unknown charts, animated geometry, prepared scene histories,
primary media and calibrated fixture paths refuse. The residual optical phase
budget does not replace the independent geometric coplanarity requirement.

For every occupied nearer/farther layer pair, `zf<zb`, a farther source with
reference-eye angular coordinate `theta_b` crosses the nearer plane at
`theta_b + a*(1/zf-1/zb)`, where `a` includes the complete pupil and the
lateral camera displacement. Each layer's projected vertex rectangle is
bounded with outward-rounded coordinate/depth intervals. The farther rectangle
is expanded over the entire signed pupil/trajectory rectangle and must be
disjoint from the nearer rectangle along at least one axis. Thus a centre ray
that misses the foreground does not hide clipping of an outer pupil segment.
The enclosing outer radius is conservative for annular and bladed openings.
Overlapping or numerically unresolved footprints refuse; the method can reject
safe arrangements whose bounding rectangles overlap.

This proof applies to the whole owning camera shot, including every linear
lateral eye segment between keyframes. Fixed orientation and axial eye
coordinate are required, up to the existing coordinate-roundoff gate; the
admitted axial-roundoff excursion is also included in source-footprint bounds.
Lateral pans and focus racks remain supported when the entire swept pupil
passes clearance. The renderer rebuilds admission for each render call and
does not reuse it across changed scene geometry or prepared histories. Its
work is one bounded scene scan and at most 120 layer-pair checks, with no
per-ray layer-visibility scan. Camera-only consumers must establish equivalent
pupil-independent visibility and source radiance themselves.

The native numerical checks retain the existing three pupil kernels and add
unit-power integration of distinct source phases under one shared polygon
proposal (32,768 Hammersley samples; 1% numerical integration band), plus
full-pupil footprint cases that distinguish clear centre rays from clipped
outer rays and include both signs of lateral camera motion. Camera regressions
check actual-source chromatic nulls, common proposals, black misses, phase
admission and focus interpolation. Renderer regressions check incoherent
addition of separated emitters, progressive replay during a lateral focus
rack, and transactional refusal of partial pupil clipping or axial motion.

This extends the clear-pupil scalar imaging model to separated source depths.
It does not simulate arbitrary 3D occlusion, clipped-pupil wavefronts, hidden
background reconstruction, inter-object edge-wave propagation, coherent
interference between scene points, polarization or high numerical aperture.
The distinction between depth-dependent PSFs and occlusion modeling is also
discussed by Ikoma et al., *ICCP* (2021), section 3.2,
<https://www.computationalimaging.org/wp-content/uploads/2021/04/DeepDfD_ICCP2021.pdf>;
the explicit clearance proof above restricts this mode to a class that does
not require their approximate layer-occlusion composition.

### Native material dispersion in thin camera lenses

`PhysicalCamera::project_geometric_spectral_point` analytically inverts one
original geometric material-dispersive ray at a supplied wavelength and pupil
sample. Fixed-size double-double reciprocal-basis contractions preserve the
authored world frame. The zero-power-difference limit uses the original finite
world-space focus construction and actual rounded pupil origin. Every result
is checked against one newly generated native ray under an explicit positive
metre residual tolerance; no iterative solve or per-query allocation is added.

The inverse retains original material/exterior source bands, optical-opening
and full-pupil paraxial gates. It refuses behind-camera targets, unsupported
wave/sequential modes and unresolved residuals. Off-screen coordinates remain
valid. This is a conditional ideal-phase-lens image, not a pupil average,
visibility certificate or finite-thickness optical model.

`compute_geometric_dispersive_motion_vectors` uses that inverse on the original
stable material-point correspondence, fixing the actual wavelength and pupil
draw across frames. It shares original identity/time admission, normals,
camera cuts, target-minus-current raster signs and behind-camera categories.
Reusing the current frame as an unused endpoint reuses the accepted projection.
Three native fixture groups cover rotated frames, circular/annular/polygonal
pupils, reference and dispersed wavelengths, source support, cancellation,
correspondence and motion signs. The changed Rust was parsed and independently
source/API reviewed; these integration fixtures were not executed here.


With `tracer`, `ThinLensElement` and `DispersiveThinLens` explicitly enable
longitudinal chromatic aberration in explicitly selected geometric,
finite-defocus and layered camera paths. A prescription contains one to eight coincident ideal thin
elements, each retaining its exact native `DielectricIor` (Cauchy, Sellmeier or
tabulated phase index) and signed front/back surface curvatures in inverse
metres. The default exterior index is one and wavelengths are vacuum nm within
the intersection of the original material laws’ source bands, including
explicit infrared support; group index is never substituted for phase index. No IOR fit, RGB conversion, external optics library or new pupil raster
is introduced. The immutable bounded prescription is shared across evaluated
camera states; per-ray material evaluation and pupil propagation allocate
nothing.

For reference wavelength `lambda_ref`, element power adds at the common pupil:

`P(lambda) = sum_j (n_j(lambda)-1)*(C_front,j-C_back,j)`.

An explicit `try_new_in_medium` prescription retains the surrounding native
phase-index law. Its power uses `n_j(lambda)/n_exterior(lambda)-1`, and its
stable chromatic difference subtracts double-double cross products of the
native material/exterior indices. Wavelength support is their intersection;
no surrounding material is guessed from the camera or scene.

The declared reference focus `s` sets the positive sensor conjugate
`t = 1/(P(lambda_ref)-1/s)`. At actual source depth `z`, the propagated pupil
phase is

`phi(lambda,z) = pi*a^2/lambda * [1/z - 1/s + P(lambda_ref)-P(lambda)]`.

Infinity has `1/z=0`. The chromatic power difference is evaluated as the
compensated sum of each element's **index difference**, before multiplication
and summation, avoiding subtraction of two large aggregate powers. Positive
and negative elements can therefore represent an achromatic doublet while
retaining its secondary spectrum. The native circular, annular and polygonal
kernels receive this same signed total phase; annular inner/outer amplitude
interference and odd-blade phase orientation remain intact. Single-plane and
source-resolved density queries and actual spectral ray weights use the same
law, including the inverse-square wavelength density Jacobian.

`PhysicalCamera::with_dispersive_lens` requires a nonzero pupil and explicit
focal-length/sensor-height projection. Its focal length must agree with the
reference lens power within 64 relative machine epsilons. The finite-defocus
or layered diffraction policy must already be selected; a geometric or
focused-Airy camera refuses the lens, and an attached lens cannot silently
switch to focused Airy imaging. A source at reference focus uses the
finite-defocus policy with that declared plane. The original `CameraProjection`
metadata is retained, while `raster_vertical_half_tan()` returns `H/(2*t)`.
Both cinematic ray construction and optical-centre reprojection use this
active scale. Focus racks therefore change finite-conjugate magnification for
the fixed physical sensor height. All wavelengths share that sensor and chief
ray mapping; this models longitudinal dispersion, not an artificial lateral
color displacement. Lens prescriptions cannot change within one shot and
their full native material coefficients/tables enter AOV camera continuity.
Without an explicitly attached lens, previous raster arithmetic, sampling
paths and continuity bytes are unchanged.

`PhysicalCamera::with_geometric_dispersive_lens` separately selects the
geometric limit of the same prescription and original pupil sampler. The
fixed sensor conjugate gives each wavelength its own tangent
`u + p*(P(reference)-1/s-P(lambda))`, evaluated from the original compensated
material index difference. `generate_geometric_spectral_rays_from_tangent_offsets`
returns those separate rays with unit aperture-average camera weight. Native
curved raster XYZ, QE and participating paths consume each lane’s own ray,
retaining absorption, refraction, visibility and photon-energy weighting.
Their actual launch phase index must equal the prescription’s native exterior
law (index one by default); a mismatch refuses the complete pixel. Source bands/gaps,
`lambda_vac/(n_exterior*opening_width) <= 0.01`, and the full pupil’s outgoing paraxial slope
remain admission conditions. This mode omits diffraction; shared-ray entry
points refuse it, and curved raster still refuses wave diffraction without
the original source-depth and clear-pupil proof. Achromatic cameras retain
their existing ray arithmetic and Philox streams.

The main cinematic tracer also consumes these distinct geometric material rays.
One native pupil draw launches four original spectral path states, with unit
ideal-phase camera weights and the original wavelength measure, path RNG,
nested-medium state, scattering, absorption and source owners. Real sequential
lens launches share only this small fixed-array adapter and retain their
original clipping/throughput. Every accepted launch checks the exact native
exterior phase index against the actual starting transport medium. A calibrated
photometric packet keeps its original wavelengths rather than reconstructing
a display observer packet.

All five aligned AOV paths retain the actual lane-zero wavelength and pupil
witness for geometric as well as sequential optics. Geometric motion uses the
native analytic inverse with a 1e-8 m forward-ray residual tolerance. Unsupported
local inverse projections clear motion validity; missing internal witnesses
remain alignment errors and cancellation remains cancellation. AOV depth stays
in the authored camera frame. Existing admitted camera output semantics and
stream identities are unchanged.

Two focused native integration groups specify four distinct chromatic emitter
hits and their CIE integral, actual hero/pupil retention, analytical translated
camera motion, public beauty/AOV parity, exterior-medium mismatch, full-field
paraxial refusal and continued shared-ray API refusal. The full modified tracer
and all changed test/adapter sources passed Rust parsing and independent
source/API review; these renderer integration groups were not executed here.

Paraxial admission requires `a*abs(C_surface) <= 0.05`, object- and image-side
pupil tangent slopes at most 0.1, and nominal raster chief-ray tangent radius
at most 0.1. The tracer checks the entire raster corner, including its aspect
ratio, independently of pixel jitter. The nominal raster limit is checked **before** drawing the PSF;
the scalar kernel's complete mathematical tails are never clipped or
renormalized. Existing scalar opening-width admission and the absolute
16-radian total pupil-phase limit still apply. Every packet wavelength and
both source-depth endpoints determine a common proposal envelope, using
`max(abs(phi(lambda,z))*lambda/lambda_proposal)` in the shared reference
coordinates. A narrow/broad full-support mixture retains both a sharp
reference wavelength and wider defocused colors. Actual source depths are
also admitted at every packet wavelength; the shortest wavelength alone
still bounds only the separate residual geometric depth error.

The existing clear-pupil scene and motion admission is unchanged: this lens
does not bypass partial-pupil occlusion, primary-material, medium, prepared
history, or finite/infinite source restrictions. This is a unit-transmission
quadratic phase model, not a finite-thickness lens ray trace. It makes no
claim for separated elements/principal-plane shifts, surface spherical
aberration, coatings/reflection, bulk absorption, sensor irradiance or
photometric exposure. It preserves normalized radiance averaging and adds
no brightness multiplier. The lensmaker and quadratic phase conventions
follow Kirchauer, TU Wien, section 4.1.2, equations (4.15)–(4.18),
<https://www.iue.tuwien.ac.at/phd/kirchauer/node51.html>; finite-conjugate
imaging and magnification follow the thin-lens construction in MIT's
*Foundations of Computer Vision*, chapter 6,
<https://visionbook.mit.edu/lenses.html>.

Focused G1/G3 regressions cover native Cauchy/Sellmeier/table units, phase
versus group index, signed achromat cancellation and secondary spectrum,
reference focus and sensor magnification, all-pupil actual-source phase,
packet proposal width, unit power and complete tails, focus-rack state, and
strict unsupported-input admission. Renderer regressions exercise native
glass dispersion for all pupil families, separate source-depth contributions
during a lateral focus rack, and exact progressive replay. A temporary native
numerical slice can execute the exact IOR/lens numerical bodies and existing
pupil kernels independently of the runtime; that check is not a substitute
for compiling or running the full camera/tracer integration tests.

### Sequential refractive camera

With `tracer`, `LensSurface` and `SequentialLens` admit up to 32 explicitly
positioned spherical caps, conic/even aspheres or planes, each with a clear
radius and the original scene-side `DielectricGlass`. The original constructor retains vacuum ports.
`try_new_with_sensor_medium` admits an authored sensor-side medium and the final
surface's authored scene exterior, including immersed optical systems.
Disjoint axial cap bounds preserve surface order; an unchanged-medium plane
acts as an internal aperture stop. Native phase-index and absorption laws keep
their original spectral bands and query refusals. Construction never invents
an optical prescription or measured glass data.

`PhysicalCamera::with_sequential_lens` places the sensor at the camera eye and
derives reference focal/focus diagnostics from the thick-system ABCD matrix,
using the actual sensor and scene indices for the physical focal distances.
Actual rays use the authored surface intersections and vector Snell refraction,
retaining chromatic and spherical aberration, finite thickness, skew rays, clipping, TIR,
and absorption over the actual distance inside each glass. Cumulative Stokes
transport retains interface-generated polarization when computing unpolarized
direct-image transmission. It does not multiply independent Fresnel averages.

`SequentialLens::paraxial_focus_at` evaluates the same reduced-angle matrix at
any wavelength in the complete original lens band, returning sensor-side
effective focal length and real object focus measured from the sensor.
`reference_focus` delegates at the retained reference wavelength. Phase-index
knots need no group derivative; chromatic queries do not refocus the camera.
The result describes direct paraxial geometry; coating phase, finite-aperture
aberrations and diffraction remain outside that diagnostic. Afocal/diverging
and virtual-object results refuse. Focused regressions compare independent
thick-lens formulas and exact near-axis Snell rays in vacuum and dispersive
immersion media, plus spectral support and phase-only knot behavior.

`SequentialLens::paraxial_focus_sensitivity_at` propagates the analytic first
vacuum-wavelength derivative through the same reduced-angle transfer, in one
allocation-free traversal. Signed focal and object-conjugate derivatives are in
metres per nanometre and include both original immersion exteriors. The original
phase/group laws supply the derivatives without wavelength perturbation or
refitting. Undefined source-table derivatives refuse while ordinary phase-only
focus remains available. Coating phase, finite-aperture aberrations, ghosts and
diffraction retain their separate models.

`SequentialLens::refocused_for_object` explicitly translates the complete
authored prescription to image a prescribed finite axial plane with a fixed
sensor. The caller supplies inclusive mechanical travel and an absolute
paraxial object-focus tolerance. The original thick-system matrix yields both
conjugate positions; double-double coefficients and stable quadratic roots
retain cancellation-sensitive terms. The nearest admissible position wins,
with the lower signed shift breaking an exact tie. Positive power, finite
arithmetic, original spectral support, ordered cap bounds wholly between
sensor and object, and the actual focus recomputed from rounded translated
coordinates are required. A currently virtual conjugate may be repositioned
into a real one; an unreachable target explicitly refuses.

The returned `SequentialLensFocus` owns the complete new prescription and
reports requested/achieved distance and common axial shift. Its reference
wavelength is the explicit focusing wavelength, so the existing
`PhysicalCamera::with_sequential_lens` directly consumes that same focus.
Surface shapes, original media, stops, compact/spectral coatings and ghost
sampling remain intact; the input prescription is unchanged. This operation
moves the whole lens, not selected internal groups, and does not choose scene
targets or optimize finite-aperture spot size. Regression cases cover vacuum
and dispersive immersion, the native exact-Snell camera path, both conjugates,
travel/refusal boundaries, retained optical metadata and physical length
scaling. Chromatic-gradient cases compare closed thick-lens laws and near-axis
exact rays, including nonvisible source support and derivative-domain refusals.

`LensSurface::with_coating` binds the existing coherent dielectric coating
owner at each physical hit, with layers explicitly ordered sensor to scene.
Its original Fresnel/Mueller response retains spectral interference, s/p phase
and film absorption. Adjacent bulk indices continue to determine Snell geometry;
the locally parallel films are collapsed onto the surface. Compact Cauchy
coatings retain their 380--780 nm source support, including zero-thickness
stacks. Active films refuse timed generation because coherent residence delay
is not modeled; zero-thickness films preserve the original bulk timing.
Each effective coating enters camera continuity, including its source identity.
Focused regressions cover quarter-wave anti-reflection transmission, unchanged
geometric refraction, bulk Beer loss, asymmetric absorbing stacks, polarization
correlation, source-band admission and the explicit timing boundary.

`LensCoating` owns immutable shared original transparent and absorbing films;
`SequentialLens::with_spectral_coating` attaches up to 64 ordered layers to a
face without changing `LensSurface`'s compact value contract. Sellmeier and
tabulated phase laws and tabulated complex-index films retain their exact
source data, authority and support, including nonvisible bands. The effective
coating replaces that face's compact coating. The prescription reintersects
all effective film and bulk domains; replacing the last visible-only film can
restore the original bulk law's broader support. Rays allocate no coating
data and use the native borrowed mixed-stack evaluator for coherent loss and
signed s/p transmission. Only effective coating identities enter camera
continuity. Active effective films refuse timed transport before clipping;
zero-thickness replacements retain their source domain and bulk transit.
Focused regressions cover nonvisible mixed-film transport against the original
owner, near-infrared quarter-wave cancellation, replacement domains, and
effective-only timing and continuity.

`LensGhostPath` and `SequentialLens::trace_ghost_ray` trace a deterministic
sensor-to-scene path with exactly two interface reflections: first at the
larger surface index, then at the smaller one. Other encounters transmit.
The shared traversal uses the original Snell/Fresnel/Mueller owners, reverses
both adjacent media and authored film order on return legs, tests every
revisited clear aperture, and retains the actual repeated Beer, geometric,
phase-index and qualified group-index path lengths. Spheres, conics and even
aspheres use the same intersection owners in both axial directions; polynomial
uniqueness remains admitted by its global slope bound. A pair visits
`N + 2*(first-second)` surfaces, at most 94 for 32 faces, with bounded work and
no per-ray heap allocation. Exact zero-reflection support contributes zero
sample mass; a weak channel retained by native Jones coherence is not erased.

`with_two_reflection_ghosts(direct_probability)` explicitly samples the
incoherent sum of the direct image and all `N*(N-1)/2` first-order ghost pairs.
The supplied probability controls proposal variance; original optics fixes
ghost strength. The first uniform lens coordinate is partitioned into path
intervals and remapped within the chosen interval, retaining the original
conditional pupil measure. Its inverse path probability enters the camera's
geometric/measure weight exactly once, including clipped zero-contribution
paths. The proposal enters camera continuity; `without_ghosts` restores the
original direct-only proposal and identity. Higher-reflection paths and
coherent interference between bulk paths are excluded. Timed paths retain
actual repeated bulk transit; active films and encountered TIR reflections
refuse timing because their coherent interface delay is not modeled.
Regressions cover analytic slab Fresnel/three-pass Beer and group delay,
signed Stokes transport through reversed spectral films, bidirectional conic
and polynomial hits, proposal measure, and the explicit direct-plus-pair sum.

`LensSurface::try_new_aspheric` retains the original vertex curvature, a
dimensionless conic constant and seven SI even-polynomial coefficients A4
through A16. Excluding A2 keeps the paraxial reference curvature unchanged.
The conventional conic sag plus polynomial departure supplies both the exact
normal and the rear projected-area/solid-angle Jacobian. Normalized polynomial
coefficients and conservative full-disk sag/slope bounds are prepared once;
surface order is admitted from those bounds, not sampled radial extrema.
Pure conics use analytic quadratic/linear intersections. Polynomial profiles
use a safeguarded solve bounded to 96 iterations after the complete slope bound
proves one forward crossing. Unqualified grazing or multiple-crossing cases
and exhausted/numerically unresolved solves refuse rather than becoming optical
misses. Zero conic and zero coefficients retain the original sphere/plane
arithmetic. Spectral Snell, Stokes, Fresnel, absorption and optional group transit
all reuse the same original transport through the resulting physical hits.

`generate_sequential_spectral_rays_from_tangent_offsets` returns separate
optional rays and weights. Uniform sampling of the projected rear cap uses its
exact projected-area/solid-angle Jacobian; the planar limit is cos^4(theta).
Clipped wavelengths retain zero sample mass. The native cinematic path tracer
continues every surviving wavelength independently and sums the original packet
measure, including native surface/volume transport and contribution AOVs. AOV
primary geometry follows the hero lane, including an absent clipped hero.
Camera continuity includes every surface, all authored aspheric parameters and
its full native glass identity, plus the complete sensor-medium law when
explicitly authored. Cinematic paths
validate the outgoing wavelength's exterior against their original initial
medium stack before transport.
Depth is axial distance from the sensor, and surface normal/albedo AOVs remain
available. Aligned motion output consumes the native conditional inverse below, using
the accepted beauty hero's actual wavelength, pupil/ghost draw and sensor
tangent. A clipped or absent hero never borrows another wavelength's geometry.

`PhysicalCamera::project_sequential_point` supplies a local inverse for an
explicit original wavelength and `LensSample`. It replays the original steady
spectral lens generator at trial sensor coordinates; the fixed original sample
retains both the direct/ghost selection and remapped conditional pupil. Snell
geometry, spectral/coating support, full optical throughput and clipping remain
with their existing owners. Group timing is not required. The target must lie
beyond the complete front cap, and every accepted residual uses a nonnegative
forward distance along an actual native ray to the target's axial plane.

`SequentialProjectionOptions` specifies the local sensor seed, physical miss
tolerance, finite-difference step, maximum sensor step and total native trace
budget. Fixed-size double-double plane residuals and a scaled 2x2 solve drive
bounded monotone backtracking. Clipped probes may shorten the step or supply
a one-sided local direction; source errors propagate. The difference stencil
guides the solve and makes no certified-gradient claim. Clipped seeds, singular
or unresolved solves, exhausted work and cancellation return explicit
`SequentialProjectionError` values. Every ray evaluation, including failed
probes, counts against the original caller's cap.

The returned `SequentialProjection` retains sensor tangents, NDC, axial depth,
achieved physical residual and actual trace count. Off-screen coordinates are
valid. This is one locally converged conditional image; it does not establish
visibility, uniqueness, all ghost images, or a pupil-averaged projection of a
defocused point. Three G1/G3 integration fixtures specify native spectral/pupil/
ghost and world-pose round trips, an independent closest-point residual,
chromatic geometry, steady absorbing coatings, and work/clipping/cancellation
refusals. The changed Rust passed parser and source/API review; the full
renderer integration fixtures were not executed in this editing environment.

`motion_vectors::compute_sequential_motion_vectors` consumes that conditional
projection for the original stable local material-point correspondence.
It holds the actual wavelength and normalized pupil/ghost draw fixed across
previous, current and next camera/object frames. The current converged tangent
seeds both target frames; each target uses its own native lens prescription.
The original time ordering, object/geometry/pose identities, normals,
camera-cut handling, target-minus-current displacement and beauty raster
arithmetic are shared with the optical-centre motion owner.

The explicit trace budget applies to each inverse. Reusing the same current
frame as an unused endpoint reuses its accepted projection without another
lens trace. Behind-camera endpoints remain categorical; clipped or unresolved
local inverses return typed errors. The result is neither a pupil average nor
a visibility certificate. Original target observations and reprojection
validation remain necessary for sample reuse. No query allocation is added.

Two native correspondence fixtures cover chromatic full-ray reconstruction,
pixel signs, cuts, behind-camera points, identity refusal and per-frame work
limits. They were parsed and independently source/API reviewed; they were not
executed in this editing environment.

All five serial and worker cinematic AOV routes retain that private optical
witness from the actual sequential beauty launch. Photometric packets use
their original lane-zero wavelength. No wavelength or pupil draw is
reconstructed, consumed again, or promoted from a different lane. The
previous reference preserves the accepted shutter phase and the current frame
serves as the cached unused next endpoint.

The complete cinematic AOV film now exports previous-motion values for admitted
native lens correspondences. Clipped seeds, unresolved local inverses and
exhausted inverse budgets clear motion validity while preserving other valid
AOVs. Missing internal witnesses or malformed options refuse; cancellation
retains the transactional tracer cancellation category. Camera cuts,
behind-camera points and unavailable correspondence remain categorical.
Ordinary camera projection and beauty integration retain their original
arithmetic. AOV semantics version 4 prevents mixing these changed motion
outputs with prior progressive/configuration/checkpoint identities.

Four original-owner fixtures cover exact beauty hero/pupil/tangent retention,
shutter-aligned nonzero motion, public film validity/EXR guides, clipping,
missing witnesses, local refusal, cuts, and ordinary-camera/beauty bit parity.
Both complete changed source files and the new fixtures passed Rust parsing,
with independent API/source review. The full renderer tests were not executed.

The film convention is sensor irradiance divided by
`pi*rear_clear_radius^2/rear_vertex_z^2`; multiplying by that factor recovers
irradiance from a radiometric scene. Exposure, ISO and QE remain explicit.
The direct-only proposal remains the default. This finite-thickness model
omits diffraction and autofocus. Its unweighted ray and optical-centre
projection entry points refuse; a distorted lens has no single projective
centre. Pose animation retains the authored prescription and fixed focus.
The rear sampling aperture is not an entrance-pupil measurement, so no f-number
is inferred.

`SequentialCameraBeam` retains one immutable wavelength, world ray/frame,
complete reciprocal detector Stokes state, outgoing scene phase index, and a
separate geometric weight. Generated detector power includes the endpoint
`(n_sensor/n_scene)^2` radiance factor exactly once; it may exceed one without
violating passive power transmission. Raw `trace_ray` transmission excludes
this endpoint factor and remains power-normalized.
Identity and analyzed spectral generators share the original lens traversal.
An authored passive analyzer is placed at the physical sensor before the
lens. Reciprocal U reversal and the camera basis's separate U/V handedness
conversion are both retained; no scalar intensity projection precedes scene
interfaces. Surviving beams whose positive transmission underflows refuse,
while aperture and TIR clipping remain explicit absent samples.

`CurvedCameraScene::sample_sequential_ray` and its original-variate counterpart
continue those complete beams through native curved refraction, scene Fresnel
responses, Beer paths, polarized sources and thermal emission. The original
scene must resolve the authored exterior index at the actual front-surface
launch point and wavelength.
With `volumes`, `CurvedRayleighCameraScene` retains the same lens state through
native molecular Mueller scattering and uniform or heterogeneous free-flight
tracking. These scene samples include lens/analyzer transmission but exclude
the rear-cap geometric factor; callers apply that factor exactly once during
pixel integration. Retained scene paths start at the lens exit.

`CurvedAreaLightCameraScene::sample_sequential_ray_with_samples` continues
an original steady sequential beam through finite-area NEE and complementary
emitter-hit MIS. The complete Stokes state, original coating/Fresnel response,
bulk Beer attenuation and endpoint normalization stay in that beam. The
steady `sample_pixel_sequential_area_lights` and
`sample_pixel_sensor_sequential_area_lights` APIs, with constant and spectral
analyzer variants, reuse the original CIE/QE lens integrators and native area
work accumulator. Analyzers are admitted over the original integration band
and installed at the physical sensor before lens traversal. Group derivatives
are not requested: active steady coatings and phase-only laws retain their
original optical domain. The separate pupil geometry factor is applied once;
clipped lanes retain the original four-lane proposal measure and zero source
work. Source errors retain their original lane and cause. Focused G1 controls
cover quarter-wave-coated axial area radiometry; G0 controls retain the native
steady CIE/QE path, detector and clipping semantics.

The explicit `generate_sequential_spectral_timed_beams*` entry points retain
sensor-to-front-surface group transit in DD seconds. Every finite lens segment
uses its original positive group index and actual geometric length, divided by
the exact SI light speed. The final exterior has no lens segment and does not
acquire a lens group-index requirement. Steady generation remains independent
of group derivatives, including at an authored phase-table cusp. Timed
generation admits required group laws before optical clipping and refuses
undefined or unrepresentable transit. Phase optical length is never substituted
for group length; interface residence times and coherent pulse reshaping remain
outside this geometric-envelope model.

`TransientCurvedCameraScene::sample_sequential_ray_at_time` and its
original-variate counterpart continue timed beams through the original native
source-history walk. Each source encounter adds its own scene group delay to
the lens delay exactly once, then retards the original DD sensor observation
clock. Source callbacks retain that sensor clock; the retained scene path and
its group delay still begin at the physical front lens surface. Steady beams
without qualified lens timing refuse. A raw prepared beam already owns its
sensor response, so an independently configured scene analyzer refuses instead
of being applied after the lens.

`CurvedCameraRaster::sample_pixel_sequential_transient_at_time` and
`sample_pixel_sensor_sequential_transient_at_time` integrate the same histories
into visible XYZ or original-band QE and photoelectron response. Their analyzer
is installed at the sensor during timed beam generation. Original spatial,
pupil, spectral and surface variate addresses are retained; clipped wavelengths
keep zero contribution and their complete proposal mass without history calls.
Even a fully clipped packet admits the sensor clock and cancellation. Original
source failures abort the complete observation. These entry points cover
instantaneous finite-source histories.

`sample_sequential_ray_during_exposure` and its original-variate counterpart
retain a rectangular `ExposureWindow` at the physical sensor. Each emitting
hit requests its original history's interval average after subtracting that
hit's lens-plus-scene group delay from both shutter endpoints. There is one
native path per wavelength, with no temporal quadrature or extra random draws.
`sample_pixel_sequential_during_exposure` and
`sample_pixel_sensor_sequential_during_exposure` reuse the original visible
and QE integrations, including optional clipped lanes and sensor-side analyzer
placement. QE energy and expected photoelectrons use the explicit DD shutter
duration and the common detector rate integrator. Default generic exposure
payloads preserve the original mandatory-ray API. Different per-pixel windows
can express authored row timing; moving geometry, time-varying optical laws,
and nonlinear sensor response are not inferred. Unsupported source interval
averages refuse rather than substituting sampled instantaneous values.
Timed sequential beams also feed the original finite-area source NEE/MIS
owner through `TransientCurvedAreaLightCameraScene::sample_sequential_ray_with_samples_at_time`
and `sample_sequential_ray_with_samples_during_exposure`. Each camera-hit or
connected-image history receives its native scene group delay plus the complete
qualified lens DD prefix exactly once, at the unchanged sensor observation
clock. Full Stokes, lens Beer/Fresnel factors and endpoint normalization stay
in the reciprocal beam; rear-cap geometric measure remains separate. Raw-beam
calls refuse a separately configured analyzer, which sequential raster callers
install at the sensor before the lens. The shutter method uses the original
source interval average on one frozen camera and connection walk.

`CurvedCameraRaster::sample_pixel_sequential_area_lights_at_time` and
`sample_pixel_sensor_sequential_area_lights_at_time`, together with their
`_during_exposure` counterparts, reuse the same native area sampler/work
accumulator and optional real-lens CIE/QE lane integrators. The original
analyzer is admitted over the integration band and installed at the sensor
before the lens. Chromatically clipped lanes contribute zero while retaining
the original four-lane proposal denominator. Fully clipped packets still
admit the clock and cancellation without source work. Source failures retain
the original lane and cause and refuse the complete pixel. Finite shutters
observe each source interval on one frozen optical path; QE energy and
photoelectrons apply the explicit shutter duration once after spectral
integration. Active-film timing restrictions are unchanged.

`CurvedCameraRaster::sample_pixel_sequential` and
`sample_pixel_sensor_sequential` integrate the same full optical seeds into
visible XYZ or original-band QE and photon responses, including infrared.
Their analyzed variants place the original passive analyzer at the sensor.
The shared raster, pupil and spectral proposal owner retains the original
counter streams. Each new result stores optional rays and native paths:
physically clipped lanes have neither, contribute zero and retain their full
one-quarter wavelength measure. No ray is fabricated, retried or normalized
by the number of surviving colors. Only the separate geometric weight is
applied after scene transport. Existing mandatory-ray pixel APIs preserve
their original behavior and refuse sequential lenses explicitly.

With `volumes`, the visible and QE sequential Rayleigh entry points reuse
those same proposal and integration helpers. Their shared tracking callback
continues complete lens/analyzer Stokes states through original molecular
Mueller events and heterogeneous acceptance draws. The generic
`CurvedSequentialRayleighPixelSample` retains per-wavelength volume and flight
counts; clipped lanes retain zero work. Zero scattering preserves the complete
surface-only sample. Native volume/work limits abort the whole observation.

With `volumes`, `TransientCurvedRayleighCameraScene::sample_sequential_ray_at_time`
and its uniform/heterogeneous original-variate counterparts continue timed
full-Stokes lens seeds through the existing molecular history walk. Every
surface source retains the physical sensor clock and its own lens-plus-scene
group delay along the actual curved, multiply scattered path. Rayleigh events
remain instantaneous; null candidates add neither history queries nor extra
collision weights and do not reset the accumulated delay or polarization.
`sample_pixel_sequential_rayleigh_transient_at_time` and
`sample_pixel_sensor_sequential_rayleigh_transient_at_time` reuse the same
optional-lane visible/QE integrations and independent surface/free-flight
counter domains. Their `TransientSequentialRayleighPixelSample` retains the
original per-wavelength collision and flight counts; clipped lanes have zero
work and retain their proposal mass. Analyzer placement remains at the sensor,
and source/cancellation failures abort the complete observation. Frozen media
and finite surface-source histories remain the original scope; fluorescent
lifetimes and evolving volume sources are not supplied.

The same molecular camera's `sample_sequential_ray_during_exposure` and its
uniform/heterogeneous original-variate counterparts integrate source histories
on one frozen refracted/scattered path. Every emitting hit receives its original
interval-average query with both sensor shutter endpoints retarded by that
hit's own lens-plus-scene group delay. Source integration adds no temporal
quadrature, transport draws or path replays. The visible
`sample_pixel_sequential_rayleigh_during_exposure` and original-band
`sample_pixel_sensor_sequential_rayleigh_during_exposure` retain the complete
sensor-side analyzer, optional clipped lanes, proposal mass and per-wavelength
collision/flight counts. Their defaulted generic exposure payloads preserve
existing mandatory-ray callers; QE energy and expected photoelectrons share
the original DD-duration integrator. Source, work-limit and cancellation errors
abort the observation. Frozen media, instantaneous molecular events and finite
surface histories retain the scope above. G1/G3/G4 regressions cover two-turn
chromatic source intervals at a large DD epoch, clipped analyzed packets, and
shared infrared QE paths with integrated energy/electron response.

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

### Multi-emitter fixture assembly

`photometry::fixture::PhotometricFixture` groups original `PhotometricSource`
and `PhotometricAreaEmitter` values under one proper-rigid fixture-to-world
transform and a shared `FixtureOpticalState`. Point emitters retain their own
local transform, angular distribution, spectrum and calibrated/dimmed/filtered
lumens. Each `FixtureSurfaceEmitter` binds its source's primitive index to a
local slot in the fixture's `geometry` list. Its original sidedness, spectrum,
clear-output calibration and optional thin-lens/gobo/splitter/diffuser law stay
with that native surface. Projector frames remain body-local; native instance
placement composes body and fixture transforms exactly once.

Shared state is an additional linear dimmer and passive 80-bin spectral
transmission applied through the existing source owners. For finite emitters,
both clear and emitted flux are attenuated, preserving passive projector loss;
replacing optics cannot restore absorbed power. Each point or surface entry
contributes independent calibrated power. There is no fixture-wide flux
normalization, RGB conversion, redistribution between emitters or extra emissive
proxy for a point's already-counted lumens. Positive attenuation that underflows
the represented flux or spectral filter computation refuses as `NumericalRange`
instead of silently becoming physical opacity.

`PhotometricRig::try_new` retains fixture input order and each fixture's emitter
order in separate point and finite-source lists. Dark emitters and empty
fixtures retain identities and ranges. Fixture IDs are nonzero and unique in
the rig; emitter IDs are nonzero and unique across both emitter types within
their fixture. Different fixtures may reuse a local emitter ID.
`FixtureEmitterIdentity` explicitly carries both. Fixture metadata reports
point/surface/geometry ranges and the compensated sum of post-control lumens,
adding point entries followed by surface entries. Admission limits bound the
combined emitter count, fixture count and native geometry count. Assembly checks
cancellation per fixture/emitter/body part and publishes no partial rig. Source
data is never modified; control changes assemble a fresh snapshot from the
original calibration. Point-only fixture literals set `surfaces: Vec::new()`.

Fixture bodies retain original native `GeometryInstance` handles, geometry
identities and materials. `append_geometry_to` mounts these instances once into
the caller's native primitive list without copying meshes or changing existing
primitive indices. Native static/animated object-ID collisions refuse. The
append is staged, so cancellation/refusal leaves primitive contents unchanged.

`replacement.replace_geometry_in(&previous, primitives, cx)` updates an already
mounted rig between discrete lighting snapshots. Both rigs must contain the
same native object-ID set; fixture/body ordering may differ. Every previous
body, including non-emitting housing, must still be mounted exactly once with
its original static instance, shared geometry handle, geometry identity,
placement and material, and without a legacy emission override. Missing,
duplicated, animated or stale targets and added/removed object IDs refuse
before any scene edit. Refusals identify the previous fixture and native object
where applicable.

Replacement instances are staged with fallible allocation and cancellation
checkpoints, followed by a final checkpoint before an infallible assignment
phase. Cancellation observed during preparation exposes no partial update;
cancellation arriving after the final checkpoint does not interrupt publication.
The scene keeps its length and primitive ordering, unrelated primitives and
supplemental area-emitter indices. Returned scene indices follow the replacement
rig's geometry order. Native geometry handles remain shared rather than copying
mesh buffers. Original and replacement rigs remain immutable, and the caller
admits new lighting jobs against the replacement rig and updated scene through
the usual geometry/transport owners.

This operation does not interpolate poses or controls over an exposure, update
an already prepared job, or claim incremental acceleration-structure refitting.
It permits a new material or shared geometry under an unchanged native object
ID; subsequent native job admission still validates the resulting optical and
geometric inputs. Regression cases compare moved multi-emitter snapshots with
fresh mounting, preserve independent controls and supplemental bindings, check
atomic refusal/cancellation and exercise a representative 300-fixture update.
These are functional cases, not throughput or physical error qualifications.

Rig path constructors resolve each finite source's rig-local geometry slot to
the actual scene slot by native object ID, independent of insertion offset or
later primitive ordering. Missing bodies, duplicate object IDs and mismatched
shared geometry, placement or material refuse; every finite source binds one
distinct fixture geometry slot. Area admission still validates the original
mesh area, aperture, material and source law. Bound sources remain borrowed;
the temporary object-index map creates no second renderer or mesh collection.
Fixture surface bindings precede supplemental scene-level area bindings, and
the existing native owner rejects duplicate emission on one primitive.

`PathIlluminanceJob` and `PathLuminanceJob` include rig-owned finite emission
through the existing source-hit, next-event, MIS, visibility, spectral and
physical-medium paths. Both retain separate point and surface emitter identity
maps. The surface map identifies the rig-owned prefix before supplemental area
bindings; it does not claim per-emitter decomposition of sampled indirect light.
Finite emission contributes to probe `continuous_lux`, while `direct_lux`
remains the deterministic point-source component. Public direct-only
`PhotometricScene` rig admission refuses active finite emitters explicitly.
Use the native path jobs for mixed or finite-source rig results.

For point rigs, direct illuminance breakdowns retain emitter/fixture
contributions in lux through the original connection/visibility calculation.
The total preserves source-order compensated summation; separately rounded
fixture totals need not reproduce its last bit. Direct and single-scattered
images retain the original `LuminanceJob` path, and physical-camera spectral
images use the original path integrator.

The fixture layer adds no GDTF/MVR parser, shared lens solver, inter-emitter
coupling or new transport/error model. Points remain far-field delta-position
sources. Finite apertures obey the represented Lambertian or ideal-projector
law; they do not establish real curved-lens glow, diffraction, aberrations or
manufacturer calibration. Existing glass/point-source refusals, spectral
discretization and unknown total-error qualifications remain in force.
Regression cases cover transformed finite/point geometry, independent and
shared controls, native camera/probe parity, visible luminance, stale binding
refusal, haze, cancellation/replay and a 300-fixture point-emitter rig.
These functional cases do not establish throughput or general physical bounds.

#### Shared physical lens controls

`FixtureOpticalState::with_lens_scale` applies positive finite multipliers to
each emitter's original gobo separation and focal length. The existing native
thin-lens model owns the coupled field-angle/focus response, including real
images, exact collimation and virtual images. These are physical distances,
not calibrated DMX zoom/focus channels. Nonidentity controls refuse point and
ordinary Lambertian emitters, including dark emitters; unit scales preserve
their existing behavior and empty fixtures remain valid.

`ThinLensProjector::with_lens_distances` retains each source's optical frame,
gobo, splitter, diffuser law and passive transmission. Focal-only changes share
the prepared diffuser unchanged; changing gobo distance updates its two metric
widths while sharing the existing positive integral tree. No texel traversal,
tree rebuild or new per-ray operation is introduced by lens refresh.
`PhotometricAreaEmitter::with_lens_distances` preserves clear and emitted
lumens bit-for-bit, including prior spectral attenuation and projector losses.
Fixture assembly reports unrepresentable distances and unsupported controls
with the original fixture/emitter identity. The original snapshots stay intact.

Native regressions check analytical focus/radiance, prepared-asset sharing,
diffuser support, identity controls and refusals. Fixture regressions compare
shared controls against independently configured native camera/probe paths.
These checks retain the existing ideal-projector and unknown calibration/model
error scope; no measured performance or manufactured-lens accuracy is claimed.

#### Shared projector component controls

`FixtureOpticalState` can replace a gobo and independently retain, remove or
replace each projector's native beam splitter and diffuser through
`FixtureOpticOverride`. Retained components stay individually authored.
Replacement gobos use their declared physical dimensions in each emitter's
optical frame; they are not implicitly fitted to the former mask. Frames,
geometry, spectra and clear-output calibration remain emitter-specific.

The final gobo, splitter, diffuser and scaled lens distances are applied
atomically by `ThinLensProjector::with_optical_state`. No obsolete intermediate
diffuser is prepared. Existing `PhotometricAreaEmitter::with_projector` derives
emitted lumens from retained clear-output flux and the final passive optical
transmission, then ordinary shared dimming and spectral filtering apply.
Opening an optic can restore its blocked power, but cannot restore individual
or shared dimmer/filter losses. Default component policies retain the original
assembly arithmetic. Unsupported point and Lambertian emitters refuse with
fixture/emitter identity, including dark sources; cancellation publishes no
partial rig.

Each immutable `GoboTransmission` shares a lazily prepared positive integral
across projector clones, diffuser laws and lens distances. Successful
preparation is cached; allocation failures are not. A different gobo owns a
different cache, and semantic equality ignores preparation state. Source
evaluation has no cache access or lock. Cached storage remains alive with the
gobo; concurrent first preparations can temporarily build duplicate trees.
This is an allocation-sharing property, not a measured throughput claim.

The optical laws retain their existing scope: achromatic piecewise-constant
gobos, ideal paraxial beam replication and passive box-slope diffusion.
Gobo-plane irises and framing shutters use the sampled-mask extension below.
Wavelength-dependent prism refraction, measured frost BSDFs and calibrated
hardware control channels remain outside these models. Regressions
exercise transformed mixed-calibration projectors through native camera and
probe paths, component power accounting, cache sharing, refusal/cancellation,
and a 300-fixture rig with 24 emitters per fixture. Calibration, discretization
and total physical error remain unqualified.

#### Polygonal iris and framing shutters

`GoboMask` combines a centered `PolygonalIris` with up to four angled
`FramingShutter` half-planes in the existing optical gobo plane. Iris inradius
and shutter offsets are metres; rotations are radians. The iris has an
explicit 3-to-64 blade count and is a regular polygon, not a circular-aperture
claim. A shutter passes `normal dot position <= offset`; its offset may be
negative. Zero iris inradius is an explicitly closed aperture.

`GoboTransmission::masked` subdivides each original texel by positive integer
X/Y factors, retaining the original physical rectangle. It clips each output
cell against the admitted half-planes and multiplies the original transmission
by the open-area fraction. It uses fixed-size convex-polygon scratch and area
integration, not point classification at the cell centre. The represented
integrated flux follows those cell averages without renormalizing blocked
power. Spatial mask edges remain piecewise constant; neither clipping
roundoff nor the spatial approximation has a certified error bound.

Output dimensions obey `MAX_GOBO_TEXELS`; invalid subdivisions and
unrepresentable scales or positive transmission products refuse. Allocation
is fallible and preparation checks cancellation within bounded cell work,
publishing only a complete validated gobo. A mask with no boundaries and
unit subdivision retains its original map.

`FixtureOpticalState::with_gobo_mask` applies the mask to each emitter's
authored or shared replacement gobo before the existing atomic projector
update. Lens distances, splitter/frost state, spectrum and clear calibration
retain their native owners. Assembly reuses a result for shared gobo/mask
storage and equal subdivisions. These transient lookup identities never
determine source order, numerical sums or published fixture/emitter IDs.
Different source maps, masks or resolutions cannot reuse stale preparation.
Closed projectors retain their identities; unsupported point and Lambertian
sources refuse even when dark.

Native regressions use analytical polygon/shutter areas, passivity and cache
isolation. These masks represent planar opaque blade geometry sampled into
the native gobo law; three-dimensional blade placement, diffraction,
manufactured-iris calibration and total physical error remain unqualified.

End-to-end fixture regressions exercise transformed camera/probe paths,
closed apertures, cancellation/replay, and 300 fixtures with 24 finite emitters
each. The large-rig case uses an independently specified one-eighth opening
and 759.375 lm total after optical and emitter controls. These are regression
references, not measured throughput or total physical-error bounds.

#### Bounded direct-area references for fixture rigs

`UnoccludedAreaReference::new_with_rig` and
`OccludedAreaReference::new_with_rig` admit the mounted rig's finite emitters,
followed by supplemental scene-level bindings. The occluded owner's
`new_with_rig_and_medium` also prepares the existing homogeneous or native
spectral-grid extinction model. These references bound the selected finite
component; point-source illuminance remains a separate calculation.

A shared borrowed-binding preparation resolves the original native scene
indices once. Mesh admission, model hashing, source-face exclusion and spectral
weight lookup all use those resolved indices. Calibrated source objects and
geometry are not cloned or replaced. Prepared `AreaReferenceBinding` records
retain each native primitive and optional fixture/emitter identity, including
dark sources; supplemental/flat sources have no fixture identity. The original
reference and progressive-job `source_bindings()` accessors expose this metadata
without duplicating it in every receiver result.

All interval area integration, visibility classification, adaptive patch
refinement and spectral extinction remain with their existing owners. Flat
inputs retain the same numeric traversal and model identity. Existing
preparation/output/work budgets and complete-receiver cancellation publication
remain in effect. Missing/stale mounts and duplicate emission bindings refuse.

Projectors still refuse, including dark projectors: the analytical source law
is uniform Lambertian emission. The occluded reference includes direct
absorption and out-scattering, with unresolved shadows retained as intervals;
it omits in-scattering, reflection, refraction and thermal emission. Exhausted
work budgets retain their existing termination and enclosure rather than a
tolerance pass. Source-construction, calibration, discretization and total
physical error remain unqualified. Regressions compare transformed/reordered
rigs against analytical and flat references, spectral-grid source mapping,
shadow/haze calculations, explicit refusals and cancellation/replay.

#### Finite-emitter and fixture illuminance attribution

Both area-reference owners provide `illuminance_breakdown` at a receiver.
`AreaIlluminanceBreakdown<T>` retains the original native aggregate result and
adds ordered `EmitterIlluminance` and `FixtureIlluminance` records. Emitter
indices resolve through the reference's `source_bindings()`; fixture and local
emitter identities remain distinct. Dark emitters stay present. Every original
rig fixture retains a record, including empty and point-only fixtures whose
finite-source contribution is zero. Supplemental sources contribute to the
total and emitter list with no invented fixture identity.

Attribution runs alongside one existing numerical solve. The unoccluded owner
collects its original analytical face intervals. The occluded owner retains a
separate raw interval and intersected enclosure for each source, updating them
with the original parent-replacement operation whenever a patch is split.
Resolved leaves remain accounted for even when they leave the adaptive
frontier. Global accumulation, patch priority, work/depth budget, termination,
diagnostic counters and numerical-model identity retain their original values.
Each source update uses its original root's stored source slot; it performs no
additional visibility, angular-integral or extinction query.

The new query checks `1 + emitters + fixtures` against the preparation budget's
`max_outputs` before allocating attribution or adaptive-frontier storage.
Reservations are fallible, cancellation publishes no partial breakdown, and
the original aggregate/progressive queries retain their existing behavior.
Fixture labels and surface ranges are copied from the already-admitted rig;
they do not alter the numeric model hash or clone source/mesh data.

Component values are interval midpoints with conditional numerical enclosures.
Their separately rounded midpoints and fixture reductions need not reproduce
the aggregate's last bit or endpoints. They share the global adaptive solve
and do not claim an independent tolerance pass per emitter or fixture. These
results cover only the selected uniform-Lambertian finite direct component;
point light, in-scattering, indirect transport and projector optics retain
their existing separate calculations or refusals. Calibration, spectral/field
discretization and total physical error remain unknown. Functional regressions
exercise transformed rigs, dark/empty/supplemental identities, adaptive shadows,
spectral extinction, budget refusal, cancellation and representative large rigs
without throughput or general physical-error claims.

### Direct and native transport

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

The same binding implements `SecondOrderIndexField` when its reconstruction
also implements `GasDensitySecondOrderField`. One allocation-free density-jet
query supplies the value, gradient and symmetric Hessian in SI units. The
material owner supplies `n'` and `n''`; the renderer applies both chain-rule
terms, `Hess(n) = n' Hess(rho) + n'' grad(rho) grad(rho)^T`. An affine density
therefore retains its nonlinear index curvature, and an on-axis zero gradient
does not erase focusing by a curved density profile. This feeds the existing
angular bundles, polarized differentiable scene media and detector projections
without a second optical law or inferred numerical Hessian. G1 checks compare
spatial gradient differences, exact on-axis focusing/defocusing, independent
neighboring rays and physical point-source irradiance. Existing smooth-field,
fixed-composition, molecular-scattering and unresolved-caustic limits remain.

`MaterialGasCompositionRefractiveField` binds a prepared molecular model to
an immutable `GasCompositionGradientField` with exact model/order identity,
component count, explicit density basis and nonzero snapshot/reconstruction
identity. Every query supplies all partial-density values and gradients;
unwritten entries refuse. The adapter computes
`grad(L) = sum_i K_i grad(rho_i)` and uses the material phase derivative, so
composition changes bend rays even at constant total mass density.
`GasCompositionSecondOrderField` supplies actual component Hessians, retaining
both terms `Hess(n) = n_L sum_i K_i Hess(rho_i) + n_LL grad(L) grad(L)^T`.
Queries use bounded stack storage without heap allocation and reuse one
material wavelength preparation for phase/group, gradient and Hessian. Native
curved rays and angular bundles bend and focus with phase index, while group
optical path uses the separately differentiated full-mixture group index.
Constant-polarizability fields retain phase/group equality and existing
arithmetic. Source-supported Cauchy fits permit chromatic refraction and
focusing at fixed total density without fitting renderer coefficients or
averaging constituent group indices. True absorption remains zero.
Source mismatches and provider/numerical failures refuse. No EOS
reconstruction, cell smoothing or chemical-species inference is supplied. Focused regressions cover spatial derivatives in both
density bases, constant-total-density bending and Snell invariance, analytic
on-axis focusing/defocusing, neighboring rays, detector irradiance, vacuum and
invalid/source/missing-row refusals. Two additional spectral regressions cover
blue/red bending, spatial Hessians, analytic lens focusing and homogeneous
and curved group-delay accumulation under the same explicit optical sources.

`MaterialGasCompositionRefractiveField::at_wavelength` retains the material's
stack-owned spectral preparation across native curved-ray stages, angular-bundle
stages, boundary segments and connecting-ray trials. Its borrowed
`PreparedMaterialGasCompositionRefractiveField` implements both existing field
interfaces without a cache, lock or new optical law. Spatial values and jets
still query the original provider; source/snapshot identity remains unchanged.
A different wavelength refuses explicitly, so a prepared blue field cannot
silently supply red refraction or delay. Focused regressions compare complete
ray/bundle results exactly in both density bases, with and without sourced
dispersion, and exercise the prepared field in the native curved molecular scene.

`MaterialGasCompositionRefractiveField::bind_scattering_field` borrows that exact
refractive model, provider, density basis and snapshot for native curved molecular
transport. A stack-owned componentwise upper-density vector must satisfy the
original dilute material law at the source endpoints; exact-wavelength evaluation
of the nonnegative ordered molecular sum supplies its tracking majorant. Local
candidates query the original partial-density reconstruction and refuse unwritten
or nonfinite rows, negative densities, changed source/units and exceeded component
bounds. The material owner supplies both the dispersive rate and the effective
depolarization, including the complete scattering-weighted Mueller response.
Native attachment uses the original GRIN source identity; a distinct adapter
content identity additionally binds its bounds. Bounds are caller-declared global
envelopes, not certification of unsampled field values, including a zero envelope.
Preparation and queries allocate no heap memory and retain the original component
limit. No wavelength fit, gray phase, independent gas law, absorption or thermal
emission is introduced. Seven focused regressions cover spectral/spatial variation,
all Mueller entries, density bases, source/bound identities, envelope/row/domain
refusals, vacuum and a real native curved collision with the original snapshot.

`AxialGasCompositionField` binds one original `AxialGasField` linear-primitive
span to the same molecular refraction, dispersion, angular-bundle and scattering
owners. Canonical species and exact molar masses must match the explicitly
assigned optical catalogue, including species that enter from zero at a span
endpoint. The original `state_at` supplies temperature, active mole fractions,
source/phase admission and EOS molar volume. At common static pressure the
adapter forms `N_i = x_i*N_A/V_m`, with analytic axial derivatives
`N_i' = N*(x_i' - x_i*T'/T)` and `N_i'' = -2*(T'/T)*N_i'`; its explicit world
axis embeds those derivatives as a gradient and symmetric Hessian. No independent
density reconstruction, calorically perfect replacement, numerical Hessian or
absorbing-gas law is inferred. The same borrowed adapter implements
`TemperatureField`, while the original elastic molecular optical law retains
exactly zero true absorption and therefore supplies no LTE emission.
`bind_refraction()` reuses the existing prepared wavelength owner and its exact
phase/group distinction. Original endpoint extrema give outward-rounded
component number-density envelopes for `bind_scattering_field`; original
spectral/dilute admission and runtime bound checks remain binding, without an
interval-certificate claim. Derivative coefficients and envelopes are prepared
once, but every query still calls the original thermochemical source evaluator,
including its allocations and intermediate-state validity checks.
Only `LinearPrimitives` with `FrozenGas` are admitted. A selected span is an
explicit closed domain with one-sided endpoint derivatives; queries beyond it
refuse even when another original station span exists. This prevents an implicit
smooth Hessian across primitive-gradient jumps. Paths and numerical interface
collars must stay within source support; multiple spans need explicit domain
handling. Four focused native regressions cover original-state/temperature/EOS
jet agreement, zero-species entrance derivatives, coherent wavelength/scattering
and bound reuse, bent native bundles against neighboring rays, and knot/domain,
phase/reconstruction/source admission. The `thermochemical_refraction` example
loads actual normalized species/NASA-9/transport triples and explicit optical
rows, feeds the original primitive-column reader and this adapter, and reports
native spectral paths, group delay, local molecular scattering and angular
footprints without replacing the source laws or treating scattering as heat.

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

Absorbed power is accumulated from each multipole's Riccati--Bessel Wronskian,
`-Im(D/m)/|Delta_a|²` and `-Im(mD)/|Delta_b|²`, without subtracting extinction and
scattering or clamping their cancellation error into heat. The original downward
recurrence retains its tail to evaluate the small-sphere magnetic loss without
cancelling its leading terms. No second recurrence, fit or transport model is
introduced. Negative dissipative flux refuses; an exactly real particle index
retains zero absorption. G1/G3 regressions compare direct 150/200-digit Bessel
references, resolve loss indices through `1e-30`, retain strong-absorption energy
balance and carry the weak coefficient through native spectral beam heating.
Finite multipole truncation and the existing host/model domain are unchanged.

`MieSolution::phase_table` adaptively prepares the native matched angular sampler
and reports observed normalization, moment and interpolation discrepancies.
`MieSphere::visible_spectrum` prepares native absorption, scattering and spectral
phase owners and compares the retained wavelength interpolation with fresh
midpoint solves. These diagnostics are sampled numerical evidence, not continuous
error bounds; insufficient angular resolution refuses. Grid density means
dispersed particle mass per mixture volume. A nonunit host requires independently
matching dielectric-interior ray ownership. The scalar preparation APIs retain intensity-only phases; full-Stokes source
preparation is specified below. The model does not claim particle interactions,
nonspherical particles, inferred material dispersion, or an inferred aerosol
distribution. A constant-index visible spectrum
does not acquire material dispersion from wavelength sampling alone.

`MieSphere::source_spectrum` and `MiePopulation::source_spectrum` prepare the
same positive absorption, scattering and wavelength-resolved angular owners on
2--16 explicitly sourced UV/visible/infrared knots. All three laws retain exactly
the supplied wavelength support and nonzero source identity; neither absorption
nor phase extrapolates into an unsupported band. The existing visible entry
points preserve their numerical responses and unsourced table identities.

`MieOpticalSpectrum::from_population_source` queries the original optical source
at every knot and fresh midpoint before solving the native population. Complex
indices may disperse; physical radii, densities, normalized number fractions and
the real host must remain fixed. Source errors remain distinguishable from Mie
numerical/work failures. The source identity is an association, not measurement
authority. Preparation retains the original recurrence, multipole, angular and
population limits, with no source query or Mie solve during transport. Sampled
midpoint discrepancies still do not certify unobserved spectral resonances.

Focused regressions verify infrared Rayleigh absorption in m²/kg, population
mass weighting, dispersive source queries, visible numerical parity and strict
support/source/geometry/budget/cancellation refusal. An infrared particle beam
also reaches the original native enthalpy timeline: only true absorbed watts
become heat, while elastically scattered watts remain radiative power. Fixed
indices must be physically supported over their declared band; this extension
does not invent missing infrared material data or relax the native Mie domain.

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

`MaterialParticleScattering::try_spectral` consumes the original
`SpectralComplexIndexStatePoint` with explicit wavelength/full-property-key
bindings. Arbitrary supported UV, visible and infrared bands enter without
relabeling a visible source or fitting a replacement material. Optical and
bulk-density bundles must share the exact card, complete physical condition and
quantity conventions. The binding identity includes the original wavelength/key
map as well as the accepted state, density, host and radius-number distribution.
An exactly zero extinction index remains zero particle absorption.

`source_spectrum` uses the SAME native Mie preparation owner for 2--16 caller
knots inside that original support. It evaluates material n/k at every knot and
fresh midpoint, then retains the exact output support in absorption, scattering
and angular tables. Material-state refresh keeps original claim pins; rebinding
after temperature changes updates optical response and particle bulk density
together. Existing cell palettes, constituent mixtures, LTE emission and
constituent heat deposition consume these spectra without new transport loops.
Source-knot interpolation and prepared Mie interpolation remain separate
approximations; midpoint diagnostics do not certify unsampled resonances.
Synthetic G0/G1/G3 regressions cover infrared Rayleigh absorption, temperature
refresh, source-map identity, two-temperature Kirchhoff transport, and exclusion
of purely scattering constituents from heat, plus source/state/work refusals.

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

### Cell-resolved optical composition

`CellPalette<T>` retains prepared optical models and a shared immutable `u32`
index per native physical cell, with explicit model/cell budgets and bounded
cancellation checks. Absorption and scattering can share the same map or use
different maps on the same original grid. `AbsorbingGrid::try_new_indexed` and
`ScatteringGrid::try_new_indexed` preserve borrowed density/temperature buffers
and require every declared model to cover the same exact wavelength interval.
Neither table order nor an unused model supplies an implicit physical material.

Exact DDA attenuation, spectral-channel free flight and ordered LTE source
integration use each traversed cell's coefficients. A stack-local cache reuses
wavelength evaluations across runs of the same material pair; there is no
per-ray palette allocation or constitutive solve. Collision likelihoods use the
complete accumulated spectral optical depth and the local extinction/scattering
ratio. The original collision cell supplies the phase law to continuation and
all direct-light/MIS techniques, avoiding position reclassification at a face.
Global `phase()`/`phase_ref()` queries return `None` for indexed scattering;
`phase_at_cell` is the authoritative angular query. Tracer semantics version 7
adds cell selection while retaining uniform-model arithmetic and random streams.

Beam power deposition and multiple-scattering heating use the same local
absorption/scattering split; elastic scattering never becomes immediate heat.
Absorption-only spectral lanes are recognized across every occupied cell, not
from a representative entry. Photometric contraction envelopes maximize over
actual occupied material pairs and their own angular integrals. Vacuum and LTE
support likewise depend on mapped, positive-density cells. Models remain
piecewise constant within each cell; no interpolation across material identities
or uncertainty certificate is inferred. Analytical regressions cover layered
spectral attenuation/source integration, collision-weighted angular response,
native direct lighting, local heating, uniform replay and admission refusals.

`MaterialGasScatteringCellGrid::try_new` admits accepted gas cells with different
compositions, molecular references and depolarization. It deduplicates exact
material model identities in deterministic first-cell order, prepares one
inverse-fourth spectrum and phase per model, and retains the original gas states
and exact density buffer. Pressure/temperature changes within one model reuse
its optical entry. The existing fixed-composition binding remains available;
both paths describe elastic, nonemitting molecular scattering.

`MaterialParticleCellGrid` binds a palette of original prepared particle spectra,
allowing material, accepted temperature and size distribution to vary by cell.
Absorption and scattering reuse the same immutable material-index map. Dilution
is checked with each mapped particle material's own bulk density; an optional
LTE field must equal that material's accepted temperature endpoint in every
cell. All declared hosts remain index one. These adapters preserve source and
state identities and feed the existing native transport; they do not add dense
particle correlations, a refracting carrier interface, or inferred gas chemistry.

### Coexisting optical constituents

`MassAbsorptionSpectrum::weighted_sum`, `MassScatteringSpectrum::weighted_sum`
and `ScatteringOptics::weighted_sum` prepare coexisting independent populations.
Weights are explicit kg of constituent per kg of the density field's reference
mass; numerical preparation never silently normalizes them. With total mixture
density, the coefficient is `rho_total * sum(w_i * kappa_i(lambda))`. Sources
must share exact wavelength endpoints. Active linear spectra retain their exact
union of knots, up to the native 16-knot capacity; a larger union refuses without
smoothing or discarding spectral features. Molecular inverse-fourth terms remain
analytic and combine with the linear contribution in a fixed-size, `Copy`
spectrum. Coefficient preparation uses fixed-order compensated sums and binds
the original weights, source identities and analytical references.

`SpectralPhaseMixture` retains each constituent spectrum and original angular
law. At every requested wavelength its phase is
`sum(w_i*kappa_s_i*p_i) / sum(w_i*kappa_s_i)`. Averaging refractive indices,
asymmetry parameters, or already-normalized phase endpoints is not this law.
The matched sampler conditionally remaps the existing random variates through
packet-lane, component and constituent-phase selection. Evaluation and MIS use
the complete packet-mixture PDF. No heap allocation, material solve, angular
table cloning or additional random draw occurs at a scattering vertex. A lane
with exactly zero scattering uses an isotropic proposal; its physical collision
weight remains zero. `OpticalMixtureBudget` limits supplied components, with a
hard cap of 64 and no nested phase mixtures. Gray-only photometry refuses these
spectral phases; native spectral transport, tracing and deposition consume them.

The photometric albedo maximum remains at the combined spectral knots: for
affine absorption `a` and `s = affine + K/lambda^4`, `K >= 0`, the derivative
numerator obeys `(s' a - s a')' = s'' a >= 0`. Interior stationary points cannot
be maxima. Outward coefficient intervals include the molecular addend; angular
integral bounds use the maximum constituent bound under convex mixing.

`MaterialOpticalMixture` borrows the original accepted gas/particle owners and
the original `fs_qty::Composition` in mass-fraction form. A recipe can contain
one molecular gas composition and multiple prepared particle materials or size
distributions. `MaterialOpticalMixtureGrid` binds a recipe palette to borrowed
TOTAL mass density and original accepted gas cells. Particle volume fractions
add as `phi = sum(rho_total*w_i/rho_particle_i)`. A present gas occupies the
remaining volume: `rho_total*w_gas = (1-phi)*rho_gas`, where `rho_gas` is the
original intrinsic EOS density. Admission checks that relation within the
documented floating-point allowance and never changes the accepted state.
Thus source composition, gas pressure/temperature, particle densities and bulk
mass/volume accounting must agree. The declared particle-volume ceiling applies
to the sum of all populations, not to each population independently.

Absorption and scattering palettes share the original immutable cell map.
An optional LTE field must match every constituent's accepted temperature;
different-temperature sources are not replaced by an averaged temperature.
Omitting the field explicitly selects nonemitting transport. Only true absorption
deposits heat. All particle hosts remain index one; volume accounting does not
establish dense-particle scattering accuracy or add coherent effective indices,
refracting carrier interfaces, dependent scattering or particle correlations.
Independent optical-depth, angular-moment, LTE, beam-power and volume/mass
regressions exercise the existing native owners; they are not experimental
material qualification or a performance benchmark.

### Constituent thermal radiation

`ThermalEmissionModel` binds original component mass weights, absorption
spectra and accepted kelvin temperatures. Each absorber is internally in thermal
equilibrium; different constituents need not share a temperature. Its emission
coefficient is `j_lambda = rho * sum(w_i*kappa_ai(lambda)*B_lambda(T_i))`.
`AbsorbingGrid::try_new_thermal_sources` derives absorption from those same
sources and shares their immutable cell-index map. The deterministic DDA and
prepared/free-flight scattering kernels use the absorption-weighted source
`sum(w_i*kappa_ai*B_i)/sum(w_i*kappa_ai)`, multiplied by the existing absorption
fraction and exact cell attenuation. No average temperature or extra emissivity
is introduced. Equal-temperature absorbers reduce to the original Planck
evaluation; a transparent wavelength emits exactly zero. Missing temperature
is permitted only for zero-weight or identically nonabsorbing components.

Source evaluation retains original spectra and fixed-order compensated sums,
with no per-ray allocation. Consecutive cells sharing a source model reuse the
last source radiance per spectral lane. The existing emission direction/PDF
owner caches importance weights from density, maximum combined absorption and
the hottest active absorber, in log space with a uniform support floor. The
hottest temperature affects only the proposal, never physical emission.
Unused models and zero-density cells cannot create a proposal. Shadows and
incident-power deposition never evaluate thermal sources.

A `CellScatteringGrid` over constituent sources must retain exactly their
wavelength support and absorption coefficients at its packet wavelengths;
incompatible replacements refuse with `ThermalSource`. Its scattering remains
independently supplied. Existing uniform-temperature and nonemitting
constructors retain their source semantics and identities. These finite-band
sources describe observed radiation, not a thermal-state debit, net radiation
exchange, fluorescence, or an inferred spectrum beyond admitted support.

### Material temperatures and constituent heating

`MaterialOpticalMixtureGrid::try_new_with_emission` selects
`MaterialMixtureEmission::None`, `CommonLte` or `ConstituentTemperatures`.
The original constructor forwards to the first two modes. Constituent mode
binds each prepared particle population's original accepted temperature and
absorption to the native thermal-source model. The current molecular gas law
absorbs no radiation, so it contributes no emission and retains its independent
accepted thermodynamic state. Gas EOS phase-volume closure, aggregate particle
dilution and source identity checks remain unchanged. An absorbing particle
without an accepted temperature refuses. No per-cell temperature copy is needed.

`deposit_beam_power_by_component` performs the existing native beam transport
once, then partitions each cell's true absorption using its original spectral
shares `w_i*kappa_ai(lambda)/sum_j(w_j*kappa_aj(lambda))`. The immutable
`ConstituentBeamDeposition` retains the complete original beam report,
mixture identity and sparse cell/component/source rows. Nonabsorbers and
zero-power cells receive exactly zero heat. Explicit row and component-work
budgets bound partitioning; cancellation or a nonrepresentable positive share
refuses the complete report. Compensated reductions check closure and retain
the native-minus-partition roundoff residual. Scattered power remains radiation
in the unchanged beam balance.

Under `optical-heating`, `ConstituentBeamDeposition::prepare_heating` binds
each retained `(optical cell, component)` exactly once to an active
`fs_thermal::CellId`. Constituents may heat different cells, and multiple
contributions may feed one cell. Admission checks complete power accounting
against the caller's tolerance. The returned owner is the existing
`OpticalHeating`, with its native thermal clock, accepted-substep energy
accounts, caloric laws and transactional observer composition. The caller owns
physical control-volume/mass assignments and constituent thermal contacts.
Optics and source power remain a fixed accepted snapshot; changed states need
rebinding. This bridge adds incident absorption, not an emission debit or a
net radiative-equilibrium solve.

### Constituent absorption after multiple scattering

`MaterialOpticalMixtureGrid::deposit_multiple_scattering_power_by_component`
runs the original native collision estimator once and partitions its retained
cell heat by the same fixed spectral absorption shares used for beam deposition.
`deposit_continuous_absorption_power_by_component` instead scores Beer heat
through native scattering free flights and uses that SAME constituent partition.
A traversed weakly absorbing cell can therefore contribute without a sampled
collision. `ConstituentScatteringDeposition` retains the complete unchanged
native report, including estimator identity, ballistic and diffuse escape,
unresolved radiation, correlated whole-history totals, sample counts, work
counts, seed and source configuration. The partition introduces no retracing,
additional random draws or thermal debit. Both reports can move directly through
the existing `into_heating` adapter with the original receiving-error and
unresolved-power admission gates.

Rows retain native sparse cell order and original component/source identities.
At a fixed cell and wavelength, every history's constituent heat is the same
constant fraction of the native cell heat; the reported mean and estimated
standard error therefore scale by that fraction. A sampled absorbing lane with
no observed deposition reports unknown constituent error (`None`), while
deterministic or structurally zero lanes retain `Some(0)`. An observed sample
can also have zero empirical variance; this remains a fixed-sample estimate,
not a confidence bound. Absent sampled cells remain unobserved, not proven
unheated. No aggregate constituent standard error across cells is invented:
the native cell report does not retain the covariance needed for that sum.

The existing constituent row/work budgets bound only retained output, without
expanding sparse histories to a full-grid inventory. Compensated sums check
partition closure against native cell and total mean powers and retain the
roundoff residual. Original collision caps and precision-limited histories
remain visible as unresolved radiation; neither scattering nor an unresolved
tail is credited as heat. This report does not create a stochastic constituent
thermal observer or assert net radiative equilibrium.

### Native thermal observer for scattered constituent heat

Under `optical-heating`, `ConstituentScatteringDeposition::into_heating`
consumes the already sampled constituent report and maps each retained
`(optical cell, component)` row exactly once to an active native thermal cell.
It moves the unchanged native `ScatteringPowerDeposition` into the existing
`ScatteringHeating` observer, so conduction, caloric storage, latent events,
authored sources, accepted-time energy accounts and rollback retain one owner.
No photon history is rerun and no escaped or unresolved radiation becomes heat.

Per-recipient heat is the compensated sum of constituent means. Its admission
error is the conservative sum of every contributing row/lane standard error;
this deliberately retains correlation between constituent shares of the same
native cell. The original whole-history total absorption error and unresolved
radiation thresholds are still checked independently. A sampled absorbing lane
with no observed native deposition has no numerical constituent error estimate
and therefore refuses thermal admission instead of silently becoming a zero-
uncertainty source. The constituent partition residual must fit the same
roundoff allowance used for the native radiation balance.

The report is consumed at admission so its potentially large sparse stochastic
state is moved rather than cloned. Different constituents may heat distinct
thermal cells, and many rows may feed one cell. Fixed-snapshot semantics remain:
changed optical/material state requires a fresh sampled report and observer.
This bridge is incident-radiation heating, not an emission debit, confidence
interval, rare-event guarantee or net radiative-equilibrium solve.

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
height, axial focus distance, and circular, annular or regular bladed pupils through the
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

Photometry semantics v10 and path semantics v14 also refuse active point
fixtures with `ConductorSurface::SMOOTH`. A mirror's delta reflection has no
finite BSDF density for the straight point-source connection, and a sampled
mirror continuation cannot hit an ideal point emitter. Omitting that reflected
light must not produce an accepted calibrated luminance or illuminance result.
Rough metal retains its ordinary point-source connections; dark point fixtures
do not impose this restriction. Smooth metal with finite emitters retains
native reflected paths and the original wavelength-dependent Fresnel loss.

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
radiance BSDF itself is generally non-symmetric. Bidirectional vertices now
retain the same body-local frame at the sampled shutter time. Camera walks and
connections evaluate the existing radiance law; light walks and camera splats
use its adjoint `f_importance(wo, wi) = f_radiance(wi, wo)`. Both directions
retain the native GGX visible-normal proposal, including geometric/shading
horizon rejection, mark roughness, and isotropic pole fallback. No rejected
sample is resampled, and the importance correction is not clamped to one.
Focused regressions compare sample/PDF/evaluation consistency, coated-metal
forward power in a dielectric host, rigidly placed frames on both sides, and
finite-emitter camera connections/light splats. Other material paths retain
their existing arithmetic. Furnace regressions cover the bounded reflected
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

#### Streaming inverse-CDF point selection

`PointFixtureSampling::ImportanceCdf { connections_per_vertex }` is an additional
opt-in policy on the same camera/probe jobs. It uses the same receiver-dependent
beam, inverse-square and cosine/phase weights, with two complete source sweeps.
The first admits every candidate and its scaled weight sum. It then generates
1-to-32 independent 53-bit integer draws in a separate Philox domain, sorts a
fixed-size stack array, and sweeps those draws through recomputed prefix weights.
Original draw slots determine final accumulation order. No source or fixture
identity is merged, and no per-vertex heap allocation is introduced.

Each CDF endpoint becomes an integer threshold `ceil(CDF * 2^53)`. The proposal
probability is the actual integer interval width divided by `2^53`, including
rounding of the CDF. Every positive candidate must retain a positive interval;
zero or non-increasing intervals refuse, even after all requested draws have
selected a source. The final endpoint must cover the entire random grid.
These probabilities affect proposal selection and inverse-probability weights;
physical intensity, source spectra and calibrated lumens retain their original
owners. Floating-point evaluation and the physical model remain unqualified.

Selected sources use the same native visibility, Beer attenuation, spectral
BSDF/phase response and multiplicity coalescing as the reservoir policy. All
active source distances are checked before angular/receiver filtering. The
direct point-fixture lux cached by a probe remains the full source sum. New
selection draws remain part of one whole path sample, so existing IID standard
errors and correlated-sample qualifications retain their meanings. Cancellation
and serial/parallel/progressive replay use the existing path commit boundary.

Selection uses `O(N + K log K)` work, `O(K)` random draws and `O(K)` bounded stack
storage, while evaluating source intensity and receiver response twice. It is
still a full source traversal and can cost more than reservoirs for small `K`
or expensive angular/phase laws. Native duplicate coalescing retains bounded
`O(K^2)` comparisons. No throughput or variance-per-time improvement is claimed.
Default `All` and existing `Importance` retain their arithmetic and
random addresses. Path semantics version 13 identifies this extension; CSV
reports use `importance-cdf`, and the headless stage example accepts
`--point-cdf-connections K` in either path mode, mutually exclusive with
`--point-connections K`.

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

### Reciprocal volumetric radiative exchange

`volumes::transport::scattering::exchange` estimates a temperature-independent
`RadiativeExchangeOperator<N>` from the original absorption/scattering geometry,
coefficient packet, supported wavelength window and scalar reciprocal phase
laws. Uniform and cell-resolved optical owners share the same implicit-capture
history primitive as incident beam deposition. Emission histories sample actual
extinction collisions even in a purely absorbing medium; the renderer's
no-scattering deterministic attenuation shortcut is not an emission sampler.
Existing incident-beam draw order and physical deposition remain unchanged.

Every absorbing cell/lane launches independent uniform-volume, isotropic
histories. With `A_i = 4*pi*rho_i*kappa_a_i*V_i`, a directed captured fraction
gives `A_i*p_ij`. The retained unordered coefficient is
`C_ij = (A_i*p_ij + A_j*p_ji)/2`, in m² sr. Its nonnegativity and shared value
in both directions make a Planck-difference heat transfer reciprocal and
passive for every retained sample. Self return cancels from net heat. The
escape coefficient `A_i*p_i,escape` couples to an external isotropic black bath.
Rows are not renormalized after symmetrization: finite-sample reciprocity and
directed emissive-capacity closure are different properties.

Reports retain all positive-density material-cell indices separately from sparse
emitter rows, including pure scatterers and currently transparent cells. Their
metadata has an independent storage limit so thermal consumers can check every
fixed optical snapshot's temperature domain, even outside the emitting subset.
Reports also retain cell volume, source identity, wavelengths,
per-emitter histories, empirical coefficient standard errors, directed
reciprocity residuals, self return, escape and unresolved tails. Admission
reconstructs directed absorption rows and checks their energy balance before
thermal consumption. Repeated receiver visits within one history are reduced
before statistical moments. Sparse pair storage, emitter count and total
segments across every emitter/lane have explicit limits and cancellation.
An unresolved collision/precision tail remains unresolved, never invented heat.
Empirical errors do not certify rare-path coverage or provide confidence bounds.

The operator is a volume-launched net-exchange formulation on fixed, index-one
geometry with elastic scattering. It does not implement boundary-launch
acceleration, reflecting/refracting interfaces, polarization, fluorescence or
temperature-varying optical coefficients. Spectral integration requires explicit
positive wavelength weights over supported source data; visible support never
implies infrared or bolometric support.

The original `thermal_surface` SI Planck law also supplies stable signed
temperature differences, analytic tangents, and positive secants. A combined
evaluation applies exchange coefficient and quadrature weight in logarithmic
space before conversion to watts and W/K. Nearly equal temperatures, overflowing
gross radiances, and underflowed intermediate products therefore do not erase
representable net heat. Equal temperatures produce exactly zero exchange.
Invalid domains or unrepresentable final physical rates are explicit refusals.

### Native enthalpy coupling for volume emission and reabsorption

Under `optical-heating`, `volumes::heating::exchange::VolumeRadiativeExchange<N>`
applies a fixed reciprocal optical operator to the original deposition thermal
state. The immutable operator and preaggregated receiver pairs are shared in
`Arc`; temperatures, prepared rates and energy accounts remain candidate-local.
Every positive-density material cell requires an explicit active thermal-cell
binding and optical coefficient-valid temperature interval. This includes pure
scatterers and transparent cells whose temperatures can change through other
native sources. Multiple optical cells may share one native receiver; their
domains intersect and transfers within that receiver cancel. An empty native
thermal state is refused because its time shortcut supplies no observer steps.

`RadiativeBath::Vacuum` is an explicit zero incident radiance. A blackbody bath
has a positive supplied temperature. Material temperatures always come from the
native caloric law, including latent storage. One signed Planck difference per
pair creates an original `HeatTransfer`; only the signed bath source enters the
native external-heat audit. Uniform material/bath temperature has exactly zero
sampled net heat. Spectral bath energy and internal-transfer magnitudes are
integrated only over actually accepted native steps using compensated accounts.

The same source evaluation supplies nonnegative Planck secant conductances.
Their receiving-row sums enter the original additional-conductance timestep
limiter, retaining native mass, caloric slope and harmonic composition with
conduction/boundary decay. A nonzero watt result with an unrepresentable zero
secant refuses before heat is applied. This is an explicit stability control,
not an accuracy claim for large nonlinear cooling steps. Source-domain,
conservation, allocation, work-budget or checkpoint refusal rolls back both
thermal state and radiation accounts through the existing observer transaction.

The source gate admits explicit per-cell/total coefficient-propagation errors
and unresolved heat allowances. Unresolved net-heat L1 scale uses the maximum
Planck difference from each emitter to a possible absorbing receiver or bath;
it vanishes at equal-temperature equilibrium. Its mean and empirical error are
sample diagnostics, not certified expectation bounds. Watt error propagation
treats the realized temperature vector as fixed: the same sampled operator
drives prior evolution, so its correlated errors do not become independent
per-step errors or final temperature confidence intervals. The snapshot remains
fixed within its declared temperature domains; dynamic optical resampling and
refractive/polarized exchange are separate constitutive extensions.

`VolumeRadiativeExchange::advance_refined` measures temporal discrepancy for
that fixed sampled operator. Coarse/fine/finer trajectories all use the original
native conduction, enthalpy and event owners. Two successive comparisons must
meet the supplied h/T/radiative-energy tolerances, with actual native substep
counts increasing on both refinements. A requested maximum dt that changes no
CFL-limited steps cannot demonstrate convergence. Resolved discrepancies must
decrease; roundoff only relaxes that ordering check, never the requested error
allowances. A native-coordinate resolution floor refuses sub-ULP accuracy
claims, and arbitrary caloric formation energy does not enlarge relative scales.

Shared total substep and source-evaluation limits include discarded trials.
Radiation accounts are compared as compensated interval increments before
rounding away small new heat against historic fluence. The returned state is
the finest actual native trajectory, with its conservation audit, measured
discrepancies, observed order when resolved, and all-trial work totals. Its
original continuation budget is restored. This is conditional time-convergence
evidence; it neither resamples the operator nor certifies spatial, spectral,
constitutive or Monte Carlo trajectory errors. Independent-operator ensembles
would be required to assess sampling variability of the thermal trajectory.
### Shared temperature-dependent optical material fields

`volumes::heating::temperature_field::TemperatureOpticalField<N>` is the shared
constitutive owner for local absorption and scattering versus native thermal
temperature. It compiles the original source spectra at the supplied wavelength
packet once and reuses `feedback::CompiledTable` for interpolation in kelvin.
The neutral field contains geometry, density, source identities, phases and
thermal/material bindings; it requires no incident beam, thermal bath or
invented optical data. Declared spectral and temperature support never expands.

Preparations read each distinct native thermal temperature once, then evaluate
each distinct receiver/table pair once. Candidate-local coefficient packets are
updated only when their actual coefficient bits change. A temperature change
on a flat constitutive table therefore preserves the optical packet. Every
positive-density cell retains its domain binding, including transparent cells
and pure scatterers. A packet from another field is refused before transport.
The original cell transport owner remains responsible for coefficient, phase
and constituent-source/Kirchhoff admission.

`TemperatureScatteringHeating::from_field` accepts an `Arc` of this field so
incident heating and volume emission can share the same compiled source laws.
Existing `try_new` callers retain their original source identity serialization,
photon streams, physical rates and energy accounts. `TemperatureScatteringTable`
and `ThermalScatteringCell` remain available at their original public paths.
The field owns neither a second temperature evolution nor a second caloric law;
immutable inputs are shared and mutable candidate packets roll back with the
native observer transaction.
### Temperature feedback in reciprocal volume radiation

`volumes::heating::exchange::feedback::TemperatureRadiativeExchange<N>` accepts
the shared `Arc<TemperatureOpticalField>` and evaluates local absorption and
scattering at each actual native source-start temperature. Emission capacities,
self absorption, escape and reciprocal cell exchange therefore respond to the
evolving material state. Fixed geometry, density and reciprocal phase laws stay
prescribed; coefficient source domains remain explicit. Construction validates
the field and sampling controls without launching photon histories.

Each replica retains common photon stream keys across temperature/time updates;
the history count is absent from those keys, so H/2H studies share prefixes.
Changed coefficient bits or a changed replica/history policy rebuild the
original reciprocal optical operator. An unchanged packet reuses its cached
transport while the original fixed-exchange owner still evaluates new Planck
rates. Every source query consumes a source evaluation, including cache hits;
actual optical rebuilds and photon segments have separate whole-advance counts.
An exact cached continuation can use a zero remaining photon allowance.

The fixed exchange owner remains the only Planck/rate/CFL/energy implementation.
Replacement snapshots inherit its compensated accepted energy history, while
their current powers come from the current native state. Domains are checked
at starts, events and endpoints, including transparent/pure-scattering cells.
The original thermal transaction owns conservation, latent storage, cancellation
and rollback. The coefficient-frozen secant limiter is not a bound on dκ/dT or
nonlinear feedback accuracy; time and photon refinement must assess that error.

`examples/temperature_radiative_balance.rs` composes incident scattering heat,
temperature-dependent emission/reabsorption, an external black bath and native
conduction through the existing `CoupledObservers`. Both radiative observers
share one compiled material field and one accepted native timeline. It exports
cumulative physical energy and temperature CSV with an explicit balance check.
Its stated synthetic source tables and three 8–16 micrometre bins demonstrate
integration, not measured material behavior or bolometric cooling. Run with
`cargo run -p fs-render --example temperature_radiative_balance --no-default-features --features optical-heating`.

### Whole-trajectory assessment of temperature-dependent radiation

`TemperatureRadiativeExchange::advance_ensemble_refined` evolves independent
optical replicas through complete native enthalpy trajectories from the same
accepted input state. Each replica shares photon streams across temperature
updates and paired time/history resolutions. This preserves correlations from
earlier radiative heat while permitting paired comparisons; existing input-state
uncertainty is conditional, not resampled. Returned trajectories are actual
native states. Mean caloric temperature is measured from each trajectory's
`T(h)`, never inferred by applying a nonlinear caloric law to mean enthalpy.

Acceptance requires two successive real native time refinements, a separate
final-state sampling gate, and a separate H/2H photon-sensitivity gate. Actual
substep counts must increase on both time refinements. Volume emission remains
sampled in a purely absorbing medium, so H/2H is never skipped just because
scattering is absent. The caller supplies independent h, T and radiative-energy
tolerances for the three gates. Formation enthalpy does not enlarge relative
allowances; a native-coordinate floor refuses sub-ULP accuracy claims instead
of silently relaxing them.

Centered, scaled scalar statistics are shared with the existing incident-heating
ensemble owner. Compensated interval energy differences retain small new heat
against large accepted historical accounts. Final sampling uses empirical
standard errors across complete trajectories; paired widths use the mean
difference plus the requested multiple of its empirical standard error. These
are finite-sample diagnostics, not confidence coverage, rare-path guarantees,
optional-stopping guarantees or certified bounds on nonlinear finite-H bias.

All accepted and discarded trials share substep, source-evaluation and traced-
segment ceilings. Retained trajectory, observable-cell and conservative operator
pair/cell-slot limits bound admission, including vacuum coefficient packets and
temporary operator replacement/construction overlap. Slot limits count retained
representations, not allocator bytes. The finest accepted trajectories retain
their refined photon policy and regain the original per-advance continuation
budgets. Native physics, source-domain, conservation and cancellation refusals
remain transactional; only unresolved assessment gates trigger further bounded
refinement. Geometry, spectral quadrature and constitutive-model errors remain
outside these trajectory assessments.

### Coherent transparent coatings on dielectric ray boundaries

The polarized forward-power path accepts an original
`CoherentDielectricCoating` on a smooth dielectric boundary. Its admitted
visible Cauchy, Sellmeier or sampled phase-index law, physical thickness and
source assertion remain intact. The adapter borrows that film through the shared
lossless-stack owner used by transmitting multilayers and conductor propagation;
explicit scalar transmittance, signed Mueller flux and native critical-angle
classification therefore come from the same coherent solution. There is no
separate single-film interference solver in this adapter, and
no synthetic permittivity/permeability or material-name preset is introduced.

`fresnel_coated_dielectric_polarized` retains complex s/p phase in both Mueller
matrices and normal Poynting transmission. Equal surrounding indices do not
remove a contrasting film. Zero thickness and index-matched films recover the
bare response; critical film waves, frustrated total internal reflection and
an opaque evanescent-film limit remain finite. Nontrivial exact grazing refuses
because its incident normal power vanishes. Boundary indices retain the native
real-index domain, with original bulk absorption applied separately by segment
transport.

`sample_smooth_coated_dielectric_polarized` reuses native direction admission,
reflection and Snell refraction. It chooses from the two positive optical powers,
placing the smaller branch first so weak reflection or transmission is not
erased by a rounded complement. `PolarizedBeam::sample_coated_interface` rotates
the actual ray frames, transfers forward power and changes the medium only on
transmission. It removes the camera radiance Jacobian exactly once.

Straight and curved scenes expose `try_new_with_optics` and an indexed
`Interaction::CoatedDielectric`; existing constructors retain their behavior.
Both owners borrow one original coating table and use the same boundary sampler.
Curved contacts evaluate the arriving and target local phase indices while
retaining the original field snapshot, geometry witness and detector route.
The existing spectral receiver/analyzers consume these path outcomes directly.

This is a coherent smooth, isotropic, locally parallel film model with both film
faces collapsed onto the geometric boundary. Film thickness changes optical
interference; no nanometre displacement, group delay, rough coated interface,
scene-wide coherent phase is inferred. This borrowed original-film adapter does
not attach arbitrary full-spectrum films to Copy camera-path materials; those
materials separately admit the compact source-bound multilayer cells described
above. Regression
cases cover independent Airy phase/work, antireflection, reversed media, critical
and evanescent films, original spectral laws, weak positive branch support, both
native scene consumers and calibrated detector power/handedness.


### Distributed illumination on temperature-dependent material fields

`volumes::heating::illumination::TemperatureIlluminationHeating<N>` admits a
bounded set of source-identified incident rays with explicit quadrature-weighted
spectral watts. It prepares one native-temperature coefficient packet and one
borrowed physical field for the whole illumination. Sources are reduced in
canonical ID order; reordering inputs or adding an inactive source does not
change existing photon streams. Ray weights describe caller-supplied spatial/
angular quadrature; no finite-emitter or beam-profile accuracy is inferred.

Every source uses the original multiple-scattering deposition and fixed heating
owners. Absorption heats native cells; ballistic/diffuse escape and unresolved
energy remain separate. The existing compensated account owner integrates time
once per accepted native interval and sums the physical source energies. Sampling,
receiving-cell error, truncation and power-balance allowances apply to the whole
illumination as well as its original per-source admission. Source standard errors
remain conditional diagnostics, not final temperature uncertainty.

Independent source/replica streams are common across temperature, time and
history-count comparisons, with common H/2H prefixes. The entire admitted
transport set is shared immutably and reused when coefficient bits are unchanged;
Planck or other coupled physical rates still evolve through their own owners.
This sampled response field may have nonlinear finite-history bias. Empty
illumination and unchanged cached continuation require no photon segments.

Source count, combined retained deposition cells, source evaluations and all-source
photon segments have explicit limits. Segment allowances cover every source and
rebuild in one native advance; no per-ray reset replenishes them. Every material
binding, including transparent/scattering-only cells, retains source-domain
checks at actual starts, events and endpoints. Original native thermal cloning
owns cancellation, conservation and atomic rollback of the entire candidate.


### Joint visible illumination and infrared thermal exchange

`volumes::heating::joint::JointOpticalHeating<NI, NR>` couples distributed
illumination and temperature-dependent volume radiation through the original
native `CoupledObservers`. Its two wavelength packets, supported source tables,
and optical grids remain independently specified. For example, visible incident
watts may heat the same native enthalpy cells whose supported infrared data
control emission, reabsorption and bath exchange. Equal wavelength counts,
identical field pointers or invented visible-to-infrared extrapolation are not
admission requirements. Consistent specimen/source bindings remain explicit
model inputs; equal material names or densities do not establish that consistency.

The adapter forwards every native source, transfer, timestep, interval, event
and acceptance hook. Both owners see the same start state and accepted duration;
conduction and latent storage retain the existing thermal implementation. Only
absorbed illumination and signed bath heat enter external optical heat. Internal
radiative transfers and conduction remain conservative native exchanges. Clock,
material-domain, allocation, work or cancellation refusal discards the whole
native candidate with both optical histories.

`examples/distributed_optical_heating.rs` evolves four conducting material cells
under six explicit weighted visible rays totaling 1300 W and a separate three-bin
8–16 micrometre radiation packet. The two fields use independent supported,
synthetic source tables and one native thermal state. Its CSV records per-cell
temperatures, absorbed/bath/internal radiation accounts and a direct energy
balance. Prescribed ray weights, spectral bins and time steps demonstrate
integration; measured-material accuracy, spatial quadrature accuracy and
bolometric cooling are outside that example's claims. Run with
`cargo run -p fs-render --example distributed_optical_heating --no-default-features --features optical-heating`.


### Joint thermal trajectories with separate incident and emission sampling

`JointOpticalHeating::advance_ensemble_refined` evolves independent complete
joint trajectories using the original native solver and shared scalar-statistics
owner. Each replica retains common source streams across temperatures and paired
resolutions. Two consecutive time comparisons must meet physical tolerances with
actual native substep counts increasing. Final sampling variability and the two
photon sensitivities have independent caller-supplied gates.

Photon comparisons lower illumination histories while holding radiation fine,
then lower radiation histories while holding illumination fine. Opposing heating
and cooling discrepancies therefore cannot cancel a simultaneous-only photon
check. Emitted volume radiation is always sampled, including pure absorption;
only exact non-scattering incident transport skips its photon comparison. The
spatial source quadrature remains the caller's responsibility.

Results preserve actual native trajectories, actual mean caloric T(h), separate
incident and radiation energy accounts, and covariance in net optical heat.
Compensated signed component accounts are combined before rounding interval or
paired differences, retaining representable small net heat against large
opposing historical fluences. Time, sampling and photon widths are empirical
finite-sample diagnostics; they do not certify confidence coverage, rare-path
resolution, optional-stopping behavior or nonlinear finite-history bias.

Separate total source/segment ceilings for each optical owner cover every trial;
a shared substep limit covers the native solver. The same remaining work is never
granted twice. Conservative retained pair/cell/report slot limits include both
geometry packets, vacuum cells, per-source headers and construction/replacement
overlap. Coarser/low-history ensembles are released before another is allocated.
Accepted trajectories regain both original continuation budgets, including the
cached radiation owner's effective source limit. Existing input-state uncertainty
is conditional rather than resampled, and native refusal remains transactional.

### Absorption and scattering through latent-enthalpy plateaus

`TemperatureOpticalField::try_new_enthalpy` compiles explicitly supplied mass
absorption and scattering spectra on the original native specific-enthalpy
coordinate. The existing temperature constructor, interpolation, transport,
source accounting and thermal solver remain their owners. This permits optical
coefficients to evolve during isothermal melting or freezing, where temperature
alone cannot determine the material's phase inventory.

Each table carries its optical source identity, caloric-law identity, finite
ordered enthalpy knots and explicit positive kelvin support. Zero and negative
enthalpy are legal reference choices. Every binding verifies the actual native
material ID, enthalpy storage basis and nonzero caloric identity; an unavailable
or different law identity refuses. Neither enthalpy nor temperature extrapolates.
Accepted endpoint/event checks include transparent and nonemitting material.

The shared field deduplicates native-state reads and coefficient interpolation.
Only changed coefficient bits invalidate transport: changing enthalpy on a
constant optical branch reuses cached transport, while changing coefficients
on a constant-temperature plateau rebuilds it. Runtime enthalpy does not enter
photon-stream identity. Incident illumination, reciprocal infrared exchange and
joint trajectory assessment consume the same field without another heat or
latent-energy source. Both absorption and scattering use the selected coordinate.

G0/G1 integration checks compare a sharp melting plateau with an independent
analytic Beer–Lambert feedback solution, exercise energy conservation and time
refinement, and check the equal-temperature radiative-bath limit while opacity
changes. G3/G4 checks cover enthalpy-datum changes, bound-law/domain refusals,
cache continuation and transactional cancellation. These are numerical/model
checks, not experimental validation. Linear interpolation of supplied effective
coefficients does not establish a universal solid/liquid optical mixing law.
Fixed density, geometry, phase functions and spectral quadrature remain explicit
inputs; this extension adds no deformation, phase kinetics or hysteresis.

`phase_change_optical_heating` exercises this path with the original
`fs-material::phase::EquilibriumEnthalpyPhaseCurve`: a synthetic 1 g specimen,
a sharp 350 K fusion plateau, independently authored visible/infrared spectra,
continuous 40 W illumination and an explicit convection event that refreezes
the specimen. The original native timeline and optical energy histories remain
continuous. Output includes actual liquid mass fraction, enthalpy, temperature,
absorption and separate optical/convection energy accounts. The example checks
conservation and traversal of both the melting and freezing plateaus. Geometry
and density are fixed, data are synthetic, and infrared quadrature covers only
6–18 micrometres. Run it with
`cargo run -p fs-render --example phase_change_optical_heating --no-default-features --features optical-heating,material-state`.

### Scheduled spectral absorption on the native physical clock

`illumination::schedule::SpectralPowerSchedule` admits bounded, ordered,
nonoverlapping linear power segments with independent wavelength-lane watts.
Segments can meet with a jump or leave dark gaps. Power is zero outside the
explicit support. Immutable identities include all times and endpoint powers;
clones share their segment storage. Interval means integrate local overlaps,
avoiding cancellation between large historical fluences and short late pulses.
They use the represented clock span and the original native integration duration
for energy, including endpoint-rounding cases.

`TemperatureIlluminationHeating::from_scheduled_field` retains the existing
temperature/enthalpy field, scattering transport, heat receivers and joint
infrared coupling. Source transitions cut the native step after all stability
limits. Interval source heat enters once, after that final cut. Native events at
the same instant precede preparation of the following interval; no synthetic
event consumes the authored event budget.

The same illumination owner now selects `AbsorptionEstimator` explicitly through
`from_field_with_estimator` and `from_scheduled_field_with_estimator`; original
constructors retain collision capture and the original source counter streams.
`from_polarized_field` admits full-Stokes source sets, while
`from_scheduled_polarized_field` accepts `ScheduledPolarizedIlluminationSource`
with physical unit-intensity Stokes states and the original spectral watt
profiles. Each lane's polarization remains fixed while pulse power scales its
admitted native response. Actual scattering uses the original Mueller paths;
isotropic absorption retains the same conservative heat owner. Canonical model
identity includes estimator, source IDs, physical frames and Stokes data,
including authored polarization in permanently dark lanes. Source ordering
does not alter histories or reductions. Power transitions reuse the response;
changed temperature-dependent coefficients rebuild it. Original phase/proposal
admission, envelope error allowances, native interval timing, joint thermal
composition and transactional cancellation remain authoritative. This adds no
coherent optical interference, propagation delay or inferred scalar-phase
polarization response.

Each wavelength/source is admitted at its complete-schedule power envelope.
Absolute uncertainty and unresolved-power allowances therefore apply at a
conservative simultaneous peak. The same immutable transport reports are scaled
by interval-mean power without reseeding, retracing or copying cell reports when
only illumination power changes. Whole dark intervals draw no incident photons,
retain reusable transport and still validate the native material domain. Changed
material coefficients rebuild transport when illumination resumes. Separate
incident, absorbed, ballistic/diffuse escape and unresolved energy accounts
retain their spectral lanes; reference sampling errors retain the admitted
envelope scale. Thermal input uses the actual combined
receiver map. Pending energy publishes only with the accepted native endpoint.

The existing work limits count every interval query, including darkness and
cache reuse. An additional total segment limit bounds source-profile storage
and at most twice as many transition entries. Whole-trajectory refinement and
joint visible/infrared coupling retain the same schedule and stream identities.
Geometry, density, direction, spectrum locations and source quadrature remain
prescribed. Exact profile fluence does not remove temporal error from freezing
nonlinear material coefficients during a step; existing refinement remains
necessary. These waveforms do not add laser propagation delay, pulse coherence,
ultrafast carrier physics, temperature-dependent source hardware or experimental
material validation.

G0/G1/G3/G4 checks cover rectangular/triangular spectral fluences despite a larger
requested timestep, independent Beer absorption, a latent-plateau fluence
solution and time refinement, shared conduction/infrared pulse boundaries,
continued dark cooling, cached scattering through darkness, split continuation,
bounded admission and transactional cancellation/domain refusal.

`examples/pulsed_optical_heating.rs` drives the original equilibrium enthalpy
phase chart with a 320 J blue rectangular pulse followed by a 200 J red triangle.
The independently authored blue/red and infrared coefficients evolve with the
actual phase inventory of a synthetic 1 g specimen. At 8.5 s, illumination ends
and an authored convection quench begins on the same native clock. The specimen
fully melts and then refreezes in darkness while infrared exchange continues.
CSV output retains temperature, liquid fraction, native specific enthalpy,
per-lane incident energy, absorption, escape, bath/convection heat and the energy
residual. Assertions check pulse fluence, full phase traversal, dark transport
work and conservation. Fixed geometry/density and finite 6–18 micrometre
infrared quadrature remain the same explicit model limits as the continuous
phase-change example. Run with
`cargo run -p fs-render --example pulsed_optical_heating --no-default-features --features optical-heating,material-state`.

### Opaque surface absorption from electromagnetic power flux

`ConductorOptics::power_balance` and `CoherentCoatingStack::power_balance`
return separate s/p reflected and absorbed fractions. They reuse the original
scaled electromagnetic boundary-field solve for the source-bound complex-index
substrate and its coherent transparent, absorbing or periodic coatings. The
explicit ordered-stack API replaces any coating attached to the substrate,
matching the existing stack reflection API.

Absorptance comes directly from normal Poynting flux,
`4*Y_i*Re(F*conj(G))/|Y_i*F+G|^2`, with compensated products for the normal flux.
It is not obtained by subtracting reflected power from one. Consequently a
representable small absorption survives a reflectance rounded to exactly one.
Numerically unresolved/nonpassive flux refuses; exact grazing gives zero heat.
Existing reflection/transmission callers retain their original arithmetic and
response storage through static output selection and perform no added absorption
calculations.

This is aggregate coating plus opaque-half-space absorption, without an inferred
individual-film, penetration-depth or thermal-cell partition. Original visible
wavelength support, positive substrate eta/k, incident-medium constraints,
source qualification and local planar/coherence assumptions remain in force.
Representative constants do not become measured material data, and small
roundoff closure does not certify optical uncertainty. Independent normal-flux
and Airy references cover weak absorption, single films, periodic stacks, s/p
response, grazing and exact zero-thickness/matched-layer identities.

### Native thermal intake for incident surface absorption

Under `optical-heating`, `thermal_surface::heating::SurfacePowerResponse<N>`
compiles reflected, transmitted and absorbed fractions from the original smooth
conductor, coherent coating, dielectric or diffuse emissivity owners. Incident
s/p power fractions resolve polarization at the supplied incidence geometry.
Opaque conductors use direct boundary flux. Clear dielectric transmission is a
separate outgoing power channel, never heat; subsequent bulk absorption belongs
to volume transport. Transmitting absorbing films use their original
power-normalized Mueller R/T solution without the camera's radiance eta-squared
factor. Positive coherent film-work integrals retain weak absorption independently
of rounded R/T; lossless films receive exactly zero heat. The thermal observer
retains even sub-roundoff absorbed energy in native enthalpy storage. Its power
closure allowance is not a certified material/numerical error bound.

Explicit gray diffuse absorption equals the original constant emissivity at
every supplied wavelength. Spectral diffuse absorption uses the original
finite-band emissivity curve and refuses extrapolation. A visible spectrum does
not establish an infrared or bolometric law. Optical source identities and
ordered wavelength packets remain part of response identity.

`SurfaceOpticalHeating<N>` admits bounded scheduled illumination sources, each
bound to an active native cell, its actual material ID and a declared temperature
domain. A source supplies physical watts already intercepted and weighted by its
external spatial/angular/spectral quadrature. This adapter applies no additional
area, cosine, inverse-square or exposure factor. Source IDs canonicalize sums;
multiple sources may heat one cell, whose temperature support is their
intersection. All sources use the same ordered wavelength packet. Immutable
responses, schedules and receiver maps are shared across candidate clones.

The original `SpectralPowerSchedule` cuts the native clock at exact source
transitions and supplies interval-mean watts after the final stability cut.
Absorbed power enters each native receiver exactly once. The original enthalpy
owner retains phase storage, conduction, boundary cooling and event ordering.
Dark intervals and perfectly reflecting/transparent sources still validate
material/temperature support. No optical ray retracing or coating solve occurs
during these fixed-response thermal steps.

Accepted histories separately retain per-lane incident, reflected, transmitted,
absorbed and signed residual joules, plus actual combined native heat input.
Compensated reductions reuse the existing optical heating sum owner. Arithmetic
power residuals and receiver-map reduction must satisfy the caller's absolute
watt allowance; that allowance does not remove or reassign energy. Positive
physical power/energy products that underflow to zero refuse, as do nonfinite
or overflowing sums. A weak direct-flux absorption may coexist with reflectance
rounded to one; its representability discrepancy remains visible in the signed
residual instead of being removed from heat.

Sources, retained profile segments and whole-advance interval evaluations are
explicitly bounded. Pending energy and receiver maps publish only with the
accepted native endpoint; cancellation, budgets and domain refusals preserve the
original thermal and optical state. This fixed constitutive response does not
infer temperature/phase-dependent optical constants, emitting-surface feedback,
scene visibility, multiple surface bounces, coherent pulse dynamics, per-layer
heat or absorption depth. Coating and opaque-substrate losses heat one assigned
cell. Existing source qualification and optical support remain authoritative.

### Native scene mirror paths coupled to thermal storage

`SurfaceIlluminationSource<1>::from_mirror_path` compiles actual
`PolarizedScene::trace_with_mirror_balances` hits into existing surface heating
sources. `SurfaceOpticalHeating::try_new_from_mirror_path` binds them directly
to native enthalpy cells. The retained positive launch power normalizes the
prescribed launch schedule; each hit receives its actual incident/launch ratio,
including earlier reflection and Beer attenuation. Original positive absorbed
power supplies its heat fraction, even when reflectance rounds to one.

`MirrorPathIllumination` pairs each history with its physically weighted launch
schedule. `SurfaceIlluminationSource<N>::from_mirror_paths` and
`SurfaceOpticalHeating<N>::try_new_from_mirror_paths` combine many histories and
wavelengths on one native thermal clock. Each exact launch wavelength selects
its lane in the common ordered packet; other lanes carry zero power and neutral
response padding. No untraced optical spectrum, wavelength interpolation or
coherence between histories is inferred. Source IDs follow path/hit order and
source/segment budgets cover the whole batch. Splitting a path requires splitting
its scheduled watts; the compiler does not divide by the number of rays.

Geometry object IDs have explicit native material/cell/temperature bindings.
Repeated visits and mirrors sharing a cell sum without a second Fresnel factor.
Optical identities retain the original conductor/coating source, actual geometry
and frame, launch beam and counter-stream path. Source IDs are consecutive in hit
order; several paths may use disjoint ranges in one existing thermal observer.
Local incident/reflected energy totals count mirror returns; they are not the
external source's one-pass energy budget. Only actual absorption enters enthalpy.
Capped surviving light, bulk absorption and other surface interactions acquire
no extra thermal credit. Sampled histories retain their supplied estimator
weights; ensemble normalization belongs in the launch schedules.

Each history remains monochromatic transport through frozen geometry/optics,
with aggregate coating/substrate heat assigned to each mirror's receiver.
Original thermal schedule cuts, caloric/latent storage, conduction, temperature domains and
rollback remain authoritative. Zero or unrepresentable power normalization,
missing receiver bindings and exceeded budgets refuse. The absorbed fraction
reconstructed from reference watts must agree with the retained original s/p
boundary-flux fraction; unresolved reference loss requires rescaling the launch
beam and retracing. Native regressions cover coated-mirror heating after
attenuation, exact pulse fluence, repeated reflection into a latent plateau,
source identity and weak positive enthalpy increments. Spectral batch checks
compare two-color coating heat and native temperature against independent Airy
and Beer formulas and preserve heat when weighted histories are split.

`SurfaceOpticalHeating::try_new_from_stateful_mirror_paths` additionally binds
each object to a `mirror_state::StatefulMirrorReceiver`. Its original conductor
optics are fixed or read that exact native heat receiver through an existing
`SurfaceConductorLaw`, including enthalpy-coordinate laws on latent plateaus.
The native path retains each physical propagation distance and surface-event
index. Replay admits only consecutive smooth mirrors, optionally ending at a
detector or opaque object, or at an explicit surface/distance cap. Intervening
dielectric, rough or diffuse events and extinguished/collision-limited histories
refuse. Geometry and homogeneous host media remain fixed.

At the accepted thermal interval start, changed n/k invalidate the complete
affected itinerary. Original Beer propagation and complex polarized reflection
recompute all downstream incident watts, reference-frame rotations, U/V phase
and positive absorption. Only local reflectance scaling would be insufficient.
Immutable receiver/law bindings are deduplicated, unchanged histories reuse
packets, and darkness checks domains without replay. Source transitions, native
caloric storage, conduction, endpoint admission and rollback remain unchanged.
The explicit start-state approximation requires time refinement; no implicit
optical solve, moving geometry, coherent interference between path visits or
extra heat from surviving capped light is inferred. Native regressions match
an independent two-cell Fresnel ODE and its first-order refinement, retain
enthalpy-driven changes at constant temperature, reproduce freshly traced
noncoplanar polarized histories and reject unsupported native endpoints.

Native integration checks cover independent conductor/film R/T/A references,
exact rectangular/triangular pulse fluence into a latent plateau, clear-glass and
lossless-stack zero heat, polarization, shared receivers, conduction and dark
cooling, Kirchhoff gray/spectral absorption, weak loss, canonical continuation,
bounded admission and transactional cancellation/domain/underflow refusals.

`examples/surface_optical_heating.rs` compares three synthetic 1 g receivers
under the same intercepted 192 J blue box and 160 J red triangular pulse: a
smooth opaque metal, the same substrate with a coherent 100 nm transparent
coating, and a clear dielectric interface. The original equilibrium caloric
chart supplies a 350 K latent plateau. Bare and coated receivers both melt;
after source switch-off at 8.5 s, native convection refreezes both on the same
clock. The clear interface reflects/transmits its illumination and remains at
300 K. CSV output includes each temperature, metal liquid fractions, separate
R/T/A energies, native convection and stored-energy residual. Assertions check
independent pulse fluence, phase traversal, dark heat history and conservation.
Fixed synthetic optical constants, prescribed incidence/intercepted power and
aggregate coating/substrate heat are explicit assumptions; radiation and bulk
absorption behind the clear interface are outside this example. Run with
`cargo run -p fs-render --example surface_optical_heating --no-default-features --features optical-heating,material-state`.

### Surface absorption follows native temperature and latent enthalpy

`thermal_surface::heating::feedback::SurfaceConductorLaw` admits explicitly
authored complex-index rows versus either native temperature [K] or specific
enthalpy [J/kg]. Every row retains its original nine-knot spectrum and source.
Rows must share one wavelength grid and exactly the same attached transparent,
absorbing or periodic coating. Interpolation acts on substrate n/k before the
original Fresnel/coherent-film solve; reflected and absorbed fractions are never
interpolated. Exact endpoints/constants and an endpoint-near FMA preserve small
positive indices without introducing another optical kernel.

An enthalpy law names the exact original caloric identity and energy datum,
requires native enthalpy storage and checks both h and T support. Negative
enthalpies remain valid when the named datum permits them. Two states on one
constant-temperature latent plateau may therefore have different optical
constants and absorbed powers. No phase fraction, solid/liquid mixing rule or
complete material-card query point is inferred. The separately declared law
source and all row identities remain bound; representative endpoint spectra
cannot become caller-asserted measured optics through this constructor.

`SurfaceOpticalHeating::try_new_with_stateful` admits fixed and evolving sources
into the original surface energy/source path. `SurfaceIncidence` retains the
prescribed wavelength packet, incident medium, angle and s/p power fractions.
The existing cell/material and source temperature domains remain active. The
additional knot budget counts state rows across distinct immutable law identities,
including shared laws only once; each row has the fixed nine spectral knots.

Receiver/material/law bindings are deduplicated. Immutable laws and prepared
optical packets share storage across candidate clones. Each source interval
resolves n/k once per distinct binding; unchanged optical values keep their
Fresnel packets. Changed packets rebuild only when nonzero illumination needs
them. Dark intervals still check native material, caloric identity and T/h
support, including at accepted endpoints and events. `response_evaluations`
reports actual packet rebuilds during an advance; admission builds are excluded.

The optical response is held at the accepted interval's START state. Native
thermal storage, events, conduction and cooling remain authoritative, and all
existing pending-history, energy reduction and rollback rules apply. Exact
scheduled fluence does not remove this explicit constitutive time error. The
feature has no fully endpoint-implicit optical claim and does not supply
temperature-dependent film thickness, coating optics, roughness, geometry or a
broadband thermal-emission law.

G1/G3/G4 checks cover independent n/k-before-Fresnel response, convergence to a
separated one-cell nonlinear ODE, actual latent-enthalpy feedback, retained
coatings, constant-table bit parity and zero redundant Fresnel builds, shared-law
budgets, canonical split continuation, and transactional dark-domain/caloric/
cancellation refusals. The original caloric owner remains the sole phase and
temperature model in every coupled test.

`examples/phase_surface_optical_heating.rs` compares an evolving conductor law
with the same initial optics held fixed under identical 320 J blue and 200 J red
pulses. Both synthetic 1 g specimens use the original equilibrium enthalpy chart
and 350 K latent plateau. The authored n/k rows follow actual native enthalpy,
so extinction changes while temperature stays at 350 K. Independent optical
histories share one native clock through `CoupledObservers`; an authored
convection change then refreezes both in darkness. The fixed control checks its
analytic normal-incidence absorption. Output retains n/k, temperatures, liquid
fractions, absorbed/reflected energies, convection, storage and Fresnel rebuilds.
Assertions require phase traversal, pulse fluence, conservation and zero unused
Fresnel rebuilds during dark cooling. The example's authored interpolation and
fixed geometry/coating assumptions are explicit; these are not measured metal
or phase-mixture data. Run with
`cargo run -p fs-render --example phase_surface_optical_heating --no-default-features --features optical-heating,material-state`.


### Native-state diffuse absorptivity

`thermal_surface::heating::diffuse_feedback::SurfaceDiffuseLaw` binds the
original finite-band spectral emissivity/absorptivity law to either native
temperature or the original specific-enthalpy coordinate. Enthalpy binding
requires the exact caloric identity and enthalpy storage basis, so two states on
one isothermal latent plateau may retain different absorptivity. Every state row
shares one authored wavelength grid; interpolation occurs in the declared state
coordinate before the existing wavelength interpolation. State refresh retains
the source band and uses fixed inline knot storage without a per-refresh heap
allocation. Infrared heating and reciprocal emission use the same native
temperature or enthalpy state and emissivity. No RGB conversion,
solid/liquid mixture rule, infrared continuation or measured-data authority is
invented.

`SurfaceOpticalHeating::try_new_with_stateful_diffuse` reuses the existing
pulse schedule, receiver, compensated energy and transactional thermal owners.
Distinct receiver/law pairs resolve once per accepted interval and invalidate
compiled spectral responses only when emissivity values change. Constant rows
retain fixed-response parity; dark intervals continue to validate source domains
without rebuilding optical packets. This coupling is start-state explicit and
therefore retains the existing temporal-refinement requirement.

### Coherent coating and substrate absorption

`ConductorOptics::layer_power_balance` and
`CoherentCoatingStack::layer_power_balance` retain the original boundary-flux
total and separately return each authored film's absorbed s/p incident-power
fraction plus the opaque substrate's fraction. Zero-thickness and lossless
films retain their positions with zero heat; compact periodic coatings expand
in physical incident-to-substrate order. The native 64-layer bound remains.

The query reuses the original scaled coherent field solution and its source,
wavelength, incidence and passive-index admission. Positive symmetric and
antisymmetric slab-loss integrals retain standing-wave interference and weak
absorption without subtracting nearly equal interface powers. Small-depth
series preserve the thin-film loss; logarithmic field-scale recovery avoids
growing propagation exponentials in opaque stacks. Output and intermediate
history use bounded stack storage. Existing reflection and aggregate-power
queries do not retain the new boundary history.

Per-film and substrate absorption must independently close the original total
within the declared floating-point allowance; unresolved or nonpassive
arithmetic refuses. These are whole-film powers, not a penetration-depth mesh,
thermal-emission model, rough-film/scattering law, broadband continuation, or
experimental validation of representative optical data. G1/G3 regressions
exercise independent matched-index Beer-Lambert layers, weak/zero loss,
oblique-polarization power closure, opaque shielding and periodic order.

### Coherent coating heat on the native thermal timeline

`CoatedSurfacePowerResponse` prepares polarization-weighted spectral film and
substrate fractions through the original coherent absorption owner. The optical
solve occurs once during preparation. `SurfaceOpticalHeating::try_new_coated`
binds every authored film, including zero-loss films, and the opaque substrate
to explicit active native thermal cells and their material identities. Multiple
receivers may share a cell; their powers add while incident, reflected and total
absorbed energy are accounted once per physical source.

The existing source schedule, compensated reductions, native enthalpy storage,
external heat sources, thermal stepping and rollback remain authoritative.
Source permutation retains deterministic identity and reductions. Every
receiver declares a positive Kelvin interval supporting its frozen optical
coefficients; repeated cells use the intersection of all declared intervals.
Initial and accepted-step validation includes film cells even when their
instantaneous absorption is zero. Failed domains, cancellation or invalid
bindings publish neither partial thermal state nor partial optical energy.

Each receiver represents whole-film power. Callers justify this thermal lumping
and the optical/material association; source identity alone does not certify
them. Absorbed energy is assumed to thermalize within the resolved interval.
Within-film penetration, temperature-varying coating optics, fluorescence,
carrier extraction and delayed nonthermal stores require their own owners.
G1/G3 tests compare distinct-cell heat against an independent matched-index
Fresnel/Beer-Lambert solution, retain existing external heat, and cover shared
cells, source permutation, receiver refusal and late transactional failure.

### Finite-film optical constants follow native thermal state

`SurfaceCoatingLaw` binds absolute film n/k to explicit temperature or specific
enthalpy rows. The native `CoherentAbsorbingCoating` still owns wavelength
interpolation, passive propagation and absorption. State interpolation precedes
the coherent solve; reflected or absorbed intensities are never mixed. All rows
retain the same nine-knot wavelength grid and exact thickness, so no expansion,
phase-mixture rule or broadband continuation is inferred. Exact transparent
endpoints and representable subnormal extinction remain intact.

Native resolution checks active cell/material identity and the declared Kelvin
domain. Enthalpy laws additionally bind the original caloric identity, energy
datum and enthalpy storage basis, allowing distinct optical states at the same
latent-plateau temperature. Neither state coordinate is extrapolated. Immutable
rows are shared; evaluation uses a bounded binary lookup and reconstructs the
original film without allocation. Source identities retain every endpoint and
cannot promote representative data to caller-asserted measured optics.
Focused G0/G1/G3 tests exercise original coherent absorption after interpolation,
weak/transparent/subnormal limits, exact latent-state binding and domain refusal.

### Independent film and substrate states drive coherent absorption

`SurfaceOpticalHeating::try_new_coated_with_stateful` binds each film and the
opaque substrate to its own original native heat receiver and constitutive law.
Fixed absorbing films, original dispersive clear films and state-dependent
films may coexist. Each temperature or enthalpy law reads that receiver's actual
material/caloric state. The complete coherent stack is then solved through the
original optical owner; local changes may redistribute heat in every layer.
Original source schedules, total R/T/A accounting, material/domain admission,
native storage/conduction and transactional acceptance remain authoritative.

Receiver/material/law pairs resolve once per interval. Exactly matching ordered
stack, receiver and incidence bindings share immutable response packets across
illumination sources. Only changed optical values invalidate dependent packets;
constant rows retain static behavior, and dark intervals defer optical solves
while still checking all state domains at accepted endpoints. Shared packet
identity never combines different physical thermal receivers. The program
identity includes the laws and geometry rather than its initial temperature.
The existing source budgets bound work, and the additional row budget counts
each distinct film/substrate law once.

The coefficients use each accepted interval's start state. This is an explicit
coupling discretization requiring time refinement, with fixed film geometry and
the original whole-film thermal lumping and prompt-thermalization assumptions.
No phase-mixture rule, thermal expansion or endpoint-implicit optics is inferred.
Four integration regressions cover bit-identical constant-row/static behavior,
independent film/substrate temperatures and coherent first-step work, three
temporal refinements against an independently integrated RK4 thermal ODE,
equal-temperature latent states, shared packet counts, permutation/continuation,
dark-domain refusal and late cancellation followed by identical retry.

### Transmitting coating heat retains authored film receivers

`TransmittingSurfacePowerResponse::from_dielectric_stack` compiles the original
dielectric coherent owner's reflected/transmitted power and positive finite-film
work once per wavelength. `SurfaceOpticalHeating::try_new_transmitting` binds
exactly those films to native thermal cells. Transmitted radiation remains an
outgoing optical channel; no substrate absorption, radiance eta-squared factor,
extra area/cosine factor or implicit exit-medium heat destination is introduced.

Film fractions and receivers retain authored outside-to-inside order under both
forward and reverse incidence. The caller supplies the appropriate incident/exit
indices and incidence angle. Original dispersive lossless films remain lossless;
weak positive film work survives below the rounding scale of R+T, and evanescent
fields under total internal reflection can heat genuinely absorbing finite films.
Native optical spectral, coherence, passivity and grazing restrictions remain
binding. Fixed responses retain the existing whole-film thermal lumping,
temperature support and prompt-thermalization assumptions.

The common source admission and thermal observer own budgets, schedules,
compensated energy reductions, material/domain checks and rollback. Empty film
stacks have no thermal receivers and can account for a bare interface on an empty
native thermal state. Native interval observers still honor waveform transitions,
work budgets and cancellation in that case; the ordinary no-observer empty-state
fast-forward remains unchanged. Existing opaque source identities and numerical
operation order are retained. G1/G3 regressions cover asymmetric reverse routing,
TIR, weak stored enthalpy, lossless dispersion, empty-state pulse/ramp schedules,
receiver/budget refusal and late cancellation without partial publication.

### Stateful transmitting films compose with opaque coating sources

`StatefulTransmittingSurfaceIlluminationSource` binds every finite film to its
own native thermal receiver and fixed, dispersive-lossless or temperature/
enthalpy optical model. `CoatingSurfaceSource` and
`SurfaceOpticalHeating::try_new_coating_stacks` compose fixed and stateful opaque
and transmitting sources in one observer. Existing opaque constructors delegate
to the same preparation and retain their physical arithmetic and identities.
Transmitting exits have fixed real spectral indices and no thermal destination;
all film/receiver associations remain in authored order under reverse incidence.

The common receiver/law bindings and coherent packet cache resolve independent
film states, including distinct enthalpies at the same latent temperature.
Direction, surrounding media, incidence and physical receivers distinguish
packets. Identical bindings share their next solve; a changed film invalidates
all dependent opaque or transmitting stacks. Native material/caloric identities,
row/source/work budgets, dark-domain validation, accepted-endpoint checks and
transactional rollback remain binding. Reflected and transmitted energy track
the same updated coherent solution as the positive film heat.

Coupling retains fixed geometry, whole-film thermal lumping, prompt
thermalization and accepted-start coefficients; it requires time refinement and
does not infer expansion, phase mixing or temperature-dependent surrounding
media. Four G1/G3 regressions cover all four source kinds with exact static
parity, reverse illumination of independently evolving films against a separate
RK4 integration of native work and outgoing R/T energy over three refinements,
latent forward/reverse cache sharing, canonical permutation/split advances,
domain/caloric refusal and late cancellation followed by an identical retry.

### Reciprocal coating emission retains each physical layer temperature

`CoatingThermalEmission` owns the original ordered absorbing/lossless film laws
and one uniform thermodynamic temperature per film. Opaque stacks additionally
retain the substrate's complex-index law and uniform half-space temperature.
The selected outgoing port is vacuum, with radiance in W/(m² sr nm) per vacuum
wavelength. A transmitting stack retains the other port's original real
dispersive index; that port has no assigned temperature or emission source.
Reverse observation traverses the optical stack in reverse while preserving
every authored film/temperature association. The observed physical port must
be vacuum for either selected side.

For a mode at wavelength 380--780 nm and outward cosine in `(0,1]`, the native
coherent owner resolves positive per-film work under reciprocal illumination
from vacuum. Each layer contributes `A_s B_lambda(T)/2` and
`A_p B_lambda(T)/2`; the original Planck owner already includes both independent
polarizations. This is the layer integral of Greffet, Bouchon, Brucoli and
Marquier's local Kirchhoff law, Phys. Rev. X 8, 021008 (2018), equations (5)--(6),
<https://doi.org/10.1103/PhysRevX.8.021008>. Isotropic reciprocal local optical
response and local thermal equilibrium within each uniform layer are required.
The coherent fields determine extraction, while independent thermal sources add
radiances. Positive work preserves weak emission below the rounded R+T deficit.

The sample returns original film contributions, optional substrate contribution,
and physical Stokes radiance with `I=L_s+L_p`, `Q=L_s-L_p`, and U=V=0 in the
outgoing interface frame. Reciprocal R/T/A fractions remain separate from emitted
radiance; transmitted light is never reassigned as emissivity. Immutable storage
copies at most 64 full laws once and shares them on clone; each query uses fixed
reference storage and allocates nothing. Temperature/source/domain admission
precedes publication, with no extrapolation, averaged temperature, infrared or
bolometric emissivity, interior-medium index correction, roughness average,
near-field exchange or inferred thermal cooling. This is a physical spectral
observation API; existing diffuse scene-material/light-sampling behavior is
unchanged.

Four G1/G3 regressions cover a nonisothermal matched-index stack against
independent Beer--Lambert depth losses and SI Planck radiance, polarized
isothermal Kirchhoff closure with the original dispersive optical owners,
asymmetric reverse temperature routing, positive emission at k=1e-20 when R+T
rounds to one, and empty/lossless/zero-thickness and support-refusal limits.

### Native coating emission resolves the original accepted thermal state

`NativeCoatingEmission` prepares an immutable spectral observation from an
accepted `DepositionSimulation`. Every `StatefulCoatingLayer` retains its own
active cell, material identity and declared temperature support. Opaque emission
also binds the original substrate receiver and fixed or state-dependent complex
index law. Transmitting emission retains only its physical finite-film receivers;
a bare lossless interface requires no cells and emits zero.

Preparation obtains temperature and exact specific enthalpy from the native
thermal owner and resolves the same fixed, temperature or caloric-enthalpy optical
laws used by absorption heating. A latent plateau may keep Planck temperature
unchanged while its optical response evolves with enthalpy. Original caloric
identity, storage basis and law domains remain binding. Repeated cell/law pairs
share one resolution; distinct optical-law rows consume the declared preparation
budget once. The existing 64-film ceiling applies before allocation. Cancellation
or any binding/domain failure publishes no snapshot and changes no thermal state.

The snapshot owns resolved optics and temperatures, records the accepted native
time, and hashes ordered physical bindings, source identity, native state and
resolved emitter identity. Later angle/wavelength queries reuse that immutable
emitter without reading a changed simulation or resolving laws again. Observation
needs no illumination schedule and neither advances nor debits thermal storage;
internal exchange and radiative cooling require a separate thermal transport
owner. The reciprocal emitter's vacuum, wavelength, polarization and local
equilibrium assumptions remain unchanged.

Three G1/G3 integration groups compare independently evolving film/substrate
states with direct coherent absorption and independent SI Planck radiance, check
snapshot immutability and identities, preserve distinct equal-temperature latent
states under reverse observation and shared-law budgets, and exercise material,
caloric, domain and cancellation refusals without altering native heat.

<!-- coating-exchange:start -->
### Finite-band isothermal coating exchange reuses reciprocal absorption

`CoatingExchangeResponse::from_emitter` prepares one physical vacuum-facing
port of an existing `CoatingThermalEmission`. All finite films and the opaque
substrate must have exactly the same temperature; a nonisothermal stack needs
internal radiative exchange before its outward emission can define full cooling.
An empty transmitting interface has no emitter and remains exactly dark.
Transmitting responses account for only the selected port; another exposed port
is a separately authored physical exchange.

`CoatingExchangeQuadrature` declares composite-midpoint wavelength and cosine
counts over exactly 380--780 nm and outward cosine `(0,1)`. Preparation sums the
original positive film/substrate absorptances and integrates
`2*pi*mu*(A_s+A_p)/2 dmu` at each wavelength. The azimuthal factor, projected area
and polarization average enter exactly once; a black interface recovers pi
steradians. Original dispersive laws, source identity, film order and observed
side remain binding. Explicit spectral-storage and coherent-solve budgets are
checked before allocation, and checkpoints precede each bounded optical solve.
Both quadrature counts require refinement for coherent spectral/angular features;
they are approximation controls, not a certified integration error.

The immutable response then returns signed outward watts and a positive W/K
secant from the original stable scaled Planck difference owner. The body has one
common temperature; the environment is an isotropic blackbody hemisphere or an
explicit empty incoming vacuum bath. Equal temperatures give exactly zero heat
and the analytic tangent. Reversing two temperatures reverses heat while retaining
the secant for the same frozen optics. Changing temperatures performs no optical
solve or allocation. Separate or balanced positive scale factors preserve
representable exchange when an intermediate area/absorptance product would
underflow or overflow; unrepresentable heat or stability slope refuses.

The response's optical identity excludes temperature while retaining the original
optical laws and physical port declaration. A state-dependent law must be resolved
at the accepted state before comparing that identity and reusing a response.
The secant is a frozen-optics exchange conductance, without an inferred derivative
of the constitutive optical law. Existing emission identity bytes are preserved.
No infrared continuation, bolometric emissivity, roughness average, internal
nonisothermal exchange or exposure/display factor is introduced.

Five G1/G3 groups verify the independent finite-band blackbody/hemisphere limit,
integrated reciprocal emission through asymmetric transmitting ports, fixed-optics
reuse, positive weak-loss exchange at adjacent temperatures, representable heat
and secant with underflowed scale products, three spectral/angular refinements,
and nonisothermal/domain/budget/cancellation refusals without partial publication.
<!-- coating-exchange:end -->

<!-- native-coating-exchange:start -->
### Native coating exchange applies signed finite-band heat exactly once

`thermal_surface::heating::coating_exchange::NativeCoatingExchange` binds each
physical exposed patch to one original isothermal thermal receiver. All its
films and opaque substrate share that receiver's actual temperature and specific
energy; the authored specimen mass and caloric law must represent their combined
storage. Separate film temperatures require an internal radiative-exchange owner.
An opaque stack or one selected vacuum port of a transmitting stack exchanges
only the original 380--780 nm band with vacuum, a fixed blackbody bath, or the
native boundary's evolving ambient temperature.

Original temperature and caloric-enthalpy optical laws resolve before coherent
quadrature. A changed optical identity rebuilds the response; unchanged optics
reuse it while Planck rates follow native temperature. Equal-temperature latent
states retain their optical enthalpy dependence. Repeated law tables share one
declared row budget; source and coherent-mode work have whole-advance limits.
Physical patches have canonical unique IDs and must fit current exposed area
after reservoir/contact coverage. Native material, caloric identity, storage
basis, geometry and mass remain bound. The caller supplies disjoint patches and
must not attach a second radiation owner to the same area and wavelength band;
ambient convection may coexist on that area.

Signed negative outward watts enter the original source hook once. The matching
positive Planck secants use the native additional-conductance limiter, composing
with conduction rather than taking independent stability minima. The observer
charges accepted-start rates before refreshing endpoint optics, and reports
per-patch outward energy separately from native external heat already booked.
Boundary events, exposed-area changes, domain failures and late cancellation
retain native ordering and atomic rollback. Frozen-step monotonicity does not
establish nonlinear time accuracy; refine thermal time and optical quadrature.

Six G1/G3 integration groups cover stiff cooling/warming/vacuum and equilibrium,
fixed-optics reuse, state-dependent optics and ambient events, latent feedback,
coverage/domain/budget/cancellation refusal, independent direct coherent/SI
Planck RK4 refinement with native convection/reservoir heat, and illumination
balanced by cooling through the original `CoupledObservers` on one timeline.
<!-- native-coating-exchange:end -->

<!-- optical-band:start -->
### Explicit optical support extends the original kernels beyond the visible band

`ConductorOptics`, `CoherentAbsorbingCoating` and `TabulatedIor` admit original
nine-knot data through opt-in `try_new_in_band` constructors. The first and last
positive finite vacuum-wavelength knots define the exact closed support; each
table retains a positive finite caller-declared maximum interpolation gap. Knots
must increase strictly and every interval must satisfy that ceiling. This keeps
the original bounded linear interpolation and source/passivity/index limits;
the gap declaration is not an interpolation-error certificate or causal spectral
fit. Coating extinction may be exactly zero; an opaque substrate retains its
strictly positive extinction. No constructor extends or extrapolates an existing
visible data table.

Actual support and retained gap limits are exposed through
`wavelength_bounds_nm()` and `max_knot_span_nm()`. Existing visible constructors
keep their exact 380/780 nm endpoints and 100 nm gap policy. Original visible
identity bytes and arithmetic are preserved; source identities/statuses remain
unchanged. Opt-in band boundaries are carried by the original samples, and
nondefault gap policies enter optical, coating and tabulated scene/lens identity.
The structs remain bounded inline Copy values. Original Cauchy and Sellmeier
laws retain visible validity; nonvisible clear films and real exterior media use
explicitly supported real-index tables rather than unchecked analytic extension.

The original complex-interface, coherent phase, polarized reflection/transmission
and positive film-work kernels consume the actual vacuum wavelength in nm.
No wavelength remapping or replacement Fresnel solver is introduced. Every
consumed film and substrate must admit the query, including attached coatings,
zero-thickness films, matched-index and grazing shortcuts. An explicit conductor
stack continues to replace any coating attached to its supplied substrate.
`DielectricCoatingStack::try_new_in_band` and `try_mixed_in_band` retain a requested
band covered by all finite films; empty/all-zero stacks still enforce that band
before their bare-interface shortcuts. A singleton requested band admits one
spectral packet lane. Scalar exterior indices remain caller inputs at that
wavelength; their original exterior laws must also admit it. Integrated consumers
must intersect all consumed material supports and require positive band width.

Four G1/G3 groups cover retained IR knots and declared gaps, visible compatibility,
independent 8--14 micrometre Fresnel and Beer--Lambert losses, an IR quarter-wave
impedance oracle, weak transmitting absorption below the R+T rounding scale,
TIR and reverse film partitions, and support refusal through attached, empty,
zero-depth, matched and grazing paths. These tests verify implementation and
model boundaries; their synthetic optical values establish no measured-specimen
accuracy or authority outside the supplied source data.
<!-- optical-band:end -->

<!-- source-band-thermal-exchange:start -->

## Source-bounded coating emission and native infrared heat

`CoatingThermalEmission` retains the positive-width intersection of every
original film law and the opaque bulk substrate or transmitting opposite-medium
law. `wavelength_bounds_nm()` reports that exact vacuum-nanometre interval.
Transparent and zero-thickness films still constrain it. An empty transmitting
interface retains its explicitly supplied real-index law's support and emits
zero; an empty opaque stack retains its original substrate's support. Disjoint
source bands are refused before preparation. An explicit coating stack continues
to replace any coating attached to its opaque substrate.

The same original positive reciprocal s/p absorption supplies emitted radiance
and the `CoatingExchangeResponse` angular integral. The latter now integrates
exactly that retained common band, with wavelength widths in nm, physical area,
projected solid angle and the native stable Planck differences/secants. Its
reported band and optical identity bind cache reuse. The explicit midpoint
counts remain convergence controls, not a quadrature error certificate. Only
supported radiation contributes; finite-band power is not bolometric cooling.

`SurfaceCoatingLaw` and `SurfaceConductorLaw` retain each row's source band and
caller-declared maximum knot span while interpolating absolute n/k in native
T or exact-caloric specific enthalpy. Rows must have identical wavelength grids
and gap policies. No wavelength/state extrapolation or upgrade of source data
status occurs. Packet-based transmitting heating admits only its supplied
ordered wavelength nodes and requires every original film to cover their range;
explicit real surrounding indices remain supplied at those exact nodes.

`NativeCoatingEmission` and `NativeCoatingExchange` consume these same resolved
laws, so supported infrared absorption, nonisothermal emitted appearance and
isothermal external cooling share the original optical owners. Existing native
storage, caloric identity, accepted-step work, cancellation, area and budget
contracts remain binding. Fixed optics reuse their angular/spectral response;
state-dependent optics invalidate it only when the resolved optical identity
changes. Independently emitting internal layers do not constitute a complete
model of internal or near-field heat exchange.

Focused physical witnesses cover independent Beer--Lambert depth work with SI
Planck sources at infrared wavelengths, common-domain and dark-layer refusals,
finite-band Planck integral convergence, original gap-policy retention, 300 K
native cooling against independent SI Planck/RK4 integration, and infrared
latent-enthalpy feedback at 310 K. The pre-existing visible native paths remain
in the same focused regression suite.

<!-- source-band-thermal-exchange:end -->

<!-- coating-gap-exchange:start -->
### Mutual coating exchange preserves directional spectral absorption

`CoatingGapResponse::from_emitters` prepares the facing vacuum ports of two
opaque, internally isothermal coherent stacks. The bodies may have different
temperatures. The geometry is the parallel, coextensive-plate approximation
with unit view factor: each emitted/reflected mode reaches the opposite face
at the same wavelength, cosine and s/p polarization. Gap-crossing phases are
averaged, while the films inside each body remain fully coherent. This is an
explicit far-field/phase-averaged approximation, with no inferred applicability
from a gap distance, evanescent transfer, gap interference or finite-edge escape.

The original owners supply positive reciprocal absorptances A1 and A2 for
each matched mode. Summing repeated reflected powers gives
`A1*A2/(A1+A2-A1*A2)`. Its scaled evaluation avoids both an underflowed numerator
and subtraction from nearly unit reflectances. Each polarization retains this
nonlinear exchange factor before the angular integration; averaging the two
surface absorptances first is physically different. The classic specular plate
limit is Howell and Siegel, NASA SP-164, volume II (1969), section 4.3.1,
equation (4-1), with spectral/directional treatment in sections 5.3 and 5.6:
<https://ntrs.nasa.gov/citations/19690020939>.

The two emitters' common positive-width source band constrains all integration.
Both original opaque and isothermal admissions remain binding; neither surface
is silently replaced by a black bath. The shared `CoatingExchangeResponse`
preparation owns quadrature, physical-area/solid-angle factors, conditioned
spectral scales and stable Planck differences/secants. Preparation charges two
original optical solves per node and checks cancellation between them; immutable
temperature queries perform no optical solves or allocation. Ordered original
optical identities distinguish the two physical sides and exclude temperature
for fixed-optics reuse. Positive watts flow from the first side to the second;
a native conservative owner must apply that same transfer to both receivers.

Four focused G1/G4 groups cover independent repeated-ray gray limits, positive
weak exchange with representable watts and W/K despite underflowed intermediate
products, the infrared blackbody integral and spectral refinement, coherent
spectral/polarized fields against independent bounce histories, temperature
reversal/equilibrium/cache identities, and opaque/isothermal/domain/work-budget/
late-cancellation refusals. No bolometric continuation of finite source data or
internal exchange between independently heated films is claimed.
<!-- coating-gap-exchange:end -->

<!-- native-coating-gap:start -->

## Conservative mutual coating exchange on native caloric stores

`thermal_surface::heating::coating_exchange::gap::NativeCoatingGapExchange`
connects independent opaque coating specimens through declared physical vacuum
gaps. Each `CoatingGapLink` has two distinct native receivers, one common area,
original source identities and explicit spectral/angular quadrature. Each side
is internally isothermal; different sides retain independent native temperature
and exact-caloric enthalpy-dependent optical laws. The original
`CoatingGapResponse` owns their supported-band, polarized multiple-reflection
exchange in the phase-averaged far-field, parallel/full-view-factor model.

The original `DepositionObserver::heat_transfers` channel applies the SAME
signed watts to both caloric stores. Radiation is internal transfer, never
external source heat and never a conductive contact. Per-link accepted energy
and per-cell net radiation are diagnostic accounts of that same heat, not
additional storage. Native masses, geometry, caloric/storage identities,
conduction, other boundaries and events remain authoritative.

All gap areas assigned to each receiver must fit its current native exposed
area. Distinct links declare disjoint physical patches, including when sharing
a receiver. Separate radiation owners must not double-own these patches/bands.
The owner rechecks coverage and every material/temperature domain at events and
accepted states. Geometry, occlusion, finite-edge escape and nonunit view factors
are not inferred from cell coordinates.

The two source-state keys retain both T and h bits. Unchanged optical identities
reuse the immutable joint angular/spectral response; unchanged source states
also reuse its Planck rates. Unique law rows share one preparation budget across
all sides. Link/rate/optical budgets bound work, with two original coherent
solves per quadrature node. The native `DepositionBudget.max_links` independently
bounds submitted internal transfers. Fixed-coefficient secants accumulate at
both endpoints and compose with the native conduction/source timestep bound;
nonlinear state-law and time accuracy require timestep refinement.

Accepted-start pair watts are charged before endpoint optical refresh. Native
candidate cloning retains the original simulation/observer on domain, coverage,
budget or late-cancellation refusal. Focused integration witnesses cover direct
SI Planck power, zero external heat and unequal caloric masses, equilibrium and
fixed-optics cache reuse, independent two-body RK4 refinement with native
conduction, shared latent-state optical laws, canonical multi-gap ordering,
events/continuation, physical area limits and both work-budget channels.

<!-- native-coating-gap:end -->

<!-- coating-external-bath:start -->
### Exterior baths exchange heat with each original coating layer

`CoatingExternalBathResponse` prepares the contribution of one selected exterior
vacuum port to every authored film and the optional opaque substrate. Each
channel retains its own uniform local temperature. Its signed outward watts are
the channel's escaped emission minus absorption of that port's isotropic black
bath, integrated only over the original common wavelength band. An explicitly
empty incoming vacuum bath has zero incident radiance; it does not authorize a
zero material temperature. The positive reciprocal work of the matched incoming
mode supplies both terms. This is the layer integral of Greffet et al.'s local
Kirchhoff relation, Phys. Rev. X 8, 021008 (2018), equations (5)--(6),
<https://doi.org/10.1103/PhysRevX.8.021008>.

This operator supplies EXTERIOR boundary terms. It does not compute radiation
emitted by one layer and absorbed by another. Such internal transfers cancel in
the specimen total but may change individual layer temperatures; accurate total
external energy alone is not complete nonisothermal radiative evolution. The
native thermal caller remains responsible for any needed internal redistribution.
An internally isothermal specimen or a single absorbing thermal subsystem has no
additional internal radiative redistribution in its aggregate thermal account.

One original coherent solve at each wavelength/cosine supplies all film channels
simultaneously. Independent s/p work enters the same projected hemispherical
integral as the existing isothermal owner. Angular data are reduced immediately;
only `(film_count + opaque_substrate_count) * wavelength_bins` spectral slots are
retained, including explicit lossless/zero-thickness channels. A distinct budget
bounds that total storage and the coherent solve count before allocation. Queries
use fixed bounded channel storage, original scaled Planck differences/secants and
compensated sums, with no allocation or optical solve. Temperatures do not change
the prepared identity; changed resolved optical laws invalidate it.

Transmitting stacks retain authored layer order under reverse observation. Each
physical vacuum-facing side has its own external port and may have an independent
bath. A material-side port is not turned into vacuum by reversing traversal, and
transmitted radiation is never a thermal destination. The existing finite-band,
local-equilibrium, planar, coherent, reciprocal and passive assumptions remain.
No diffuse view factor, infrared extrapolation or internal-transfer coefficient
is inferred.

Five focused groups verify independent oblique Fresnel/Beer depth work with SI
Planck sources and quadrature refinement, isothermal reduction and asymmetric
reverse routing, a complete two-bath single-film reference using independent
Fresnel amplitude sums, positive weak-loss/subnormal conductances, and empty/dark
channels, exact slot/work budgets, temperature shape and late cancellation.
<!-- coating-external-bath:end -->

<!-- native-coating-external-bath:start -->
### Native per-layer exterior radiative boundary

`NativeCoatingExternalBathExchange` composes the original layer-resolved
exterior response with `DepositionSimulation::advance_observed`. Each
`CoatingExternalBathSource` supplies one selected vacuum-facing port, its bath,
area and quadrature, plus an `ExternalCoatingStack` of the original
`StatefulCoatingLayer` receivers and optional opaque-substrate receiver.
Every channel retains its original material, caloric identity, storage basis,
mass and geometry. Films read their own accepted temperature and exact native
specific enthalpy, including optical changes on a latent-temperature plateau.
Reverse transmitting traversal preserves authored film/receiver association.

The source's `surface_cell` owns the physical patch area and its optional native
ambient bath. Area is charged once per port to that cell's CURRENT exposed
surface. Buried film/substrate receivers need no exposed area. Repeated receiver
cells sum their channel watts and positive frozen-optics secants. Several ports
share the original-law row budget; retained channel-by-wavelength slots and
cumulative optical/rate work have separate explicit bounds. IDs are canonical,
and duplicated source IDs refuse. Scalar area checks do not prove geometric
non-overlap or detect another radiation owner on the same patch/band.

One coherent solve supplies all channels of a port. Immutable quadrature
responses are shared by clones and reused while resolved optical identity is
unchanged. Temperature/enthalpy state keys invalidate native optical snapshots;
unchanged states and baths reuse rates. Native ambient events refresh the bath
even when optics are unchanged. Secants combine with native conduction/source
step bounds. Accepted-start watts and accepted duration own the interval before
endpoint refresh. Domain, area, budget or cancellation refusal leaves the
accepted native simulation and observer unchanged through transactional advance.

Signed channel heat enters native `external_heat_j` once. Port energy and
`energy_by_cell_j` report that already-applied heat; they are diagnostic accounts,
not additional thermal stores. Compose independent supported transfers through
the existing native observer machinery.

This is an EXTERIOR boundary contribution. Internal radiative redistribution
between separately heated absorbing layers remains a separate operator; total
external closure alone does not certify complete individual-layer trajectories.
Both ports of a vacuum/vacuum specimen must be authored when exposed. A complete
isothermal specimen, or a single absorbing thermal subsystem surrounded by
lossless films, then has its complete far-field boundary within the admitted
source band. Material-side mode density, near-field transport, bolometric tails
and internal layer-to-layer exchange are not inferred.

The four native `coating_exchange/external.rs` integration groups check distinct
hot/cold receivers behind occupied buried faces, repeated-cell isothermal
reduction, two-sided single-absorber temperature feedback against independent
SI-Planck/native-field integration and RK4, latent enthalpy optics, ambient
events and split continuation, global slot/law/work limits, and atomic late
cancellation/coverage refusal. The low-level exterior tests supply independent
thin-film amplitude and layer-depth absorption references.
<!-- native-coating-external-bath:end -->

<!-- coating-transmitting-pair:start -->
### Coherent transmitting specimens joined by a phase-averaged radiation gap

`CoatingTransmittingPairResponse` closes the four-terminal finite-band radiation
network of two isothermal transmitting specimens and two independent exterior
baths. The selected side of each original `CoatingThermalEmission` faces the
gap; reversing that same authored stack supplies its exterior response. Both
half spaces must be identically vacuum throughout their retained original law,
not merely at selected quadrature wavelengths. Asymmetric films retain distinct
gap-side and exterior absorptances. Empty transparent stacks exchange no heat
while passing exterior bath radiation.

At each matched vacuum wavelength, cosine and s/p polarization, solve the two
incoherent gap-intensity equations using the original coherent reflection,
transmission and positive absorption. Their elimination yields one reciprocal
body/body channel, four body/bath channels and one bath/bath channel. The gap
denominator is `(A1+T1)+(A2+T2)-(A1+T1)*(A2+T2)`, evaluated after dividing by the
larger positive escape fraction. No rounded `1-R1*R2` residual, product of
independently averaged spectra/polarizations or invented opaque approximation
replaces the original matched modes. Perfectly trapped modes carry no gap
exchange. Each edge multiplies its own stable original Planck difference and
secant; global thermal equilibrium gives exact zero heat.

Four coherent solves per spectral/angular node serve all six channels. Budgets
admit `6*N_lambda` retained channel nodes and `4*N_lambda*N_mu` optical solves
before allocation or solve; original per-stack 64-film limits remain. Subsequent
temperature queries allocate nothing and preserve representable weak heat and
conductance through the original logarithmic Planck scaling. Positive body
temperatures and supplied bath temperatures are required even for dark ports;
an absent bath denotes zero incoming radiance. Ordered optical identities,
common source band, physical facing area and quadrature enter response identity.

This is unit-view-factor, parallel/coextensive, specular far-field exchange;
the caller must justify phase averaging between stacks. Interference inside
each specimen remains coherent. Neither gap interference, evanescent transfer,
finite edges, rough scattering, unsupported spectral tails nor nonisothermal
internal film exchange is implied. The intensity construction follows Byrnes,
*Multilayer optical calculations*, arXiv:1603.02720v5 section 6; emission retains
the original local Kirchhoff owner and its source references.

Six focused low-level groups check independent unreduced intensity balances,
original asymmetric coherent infrared fields and four-terminal conservation,
exact equilibrium and mirror reciprocity, transparent bath throughput against
independent SI-Planck quadrature, opaque recovery of the original pairgap and
exterior owners, weak positive work below a rounded reflectance deficit,
full-band vacuum/source admission, budgets, cancellation and dark temperatures.
<!-- coating-transmitting-pair:end -->

<!-- native-coating-transmitting-pair:start -->
### Complete native radiation for two transmitting bodies and two baths

`heating::coating_exchange::transmitting_pair::NativeCoatingTransmittingExchange`
connects two independently stored isothermal coating specimens across an
explicit phase-averaged far-field vacuum gap. The original coherent optical
stacks determine both their gap-facing and opposite exterior responses. Their
common supported source band, matched wavelength/cosine/polarization, unit view
factor and parallel/coextensive geometry remain binding assumptions. Each
`CoatingTransmittingPair` declares two gap-facing transmitting stacks, one common
physical face area and two independent exterior baths. Both half spaces of each
stack must be vacuum throughout its retained band. Each body may have its own
temperature/enthalpy optical laws and native caloric/storage identity.

The six-channel `CoatingTransmittingPairResponse` owns all repeated optical
reflections. One equal-and-opposite native `HeatTransfer` carries the body pair
exchange. Four signed body/bath channels enter native external heat exactly
once. A bath can reach either body through the other transmitting specimen;
`NativeAmbient` follows the original boundary outside its corresponding body.
Direct bath-to-bath throughput is an exterior diagnostic only, never a thermal
source or a subtraction from either body's positive absorption. Transparent
bodies therefore retain zero heat even with nonzero bath throughput.

Both physical faces belong to this owner: each body is charged `2 * area_m2`
against its CURRENT native exposed area. Native contact/reservoir events can
therefore invalidate coverage, including at an advance endpoint. Distinct pairs
declare disjoint patches; scalar area cannot prove non-overlap or detect another
owner on the same physical ports and wavelength band. No conductive contact,
finite-edge view factor, gap interference, near-field transfer or bolometric
tail is inferred. Internal radiation between nonisothermal films remains outside
this internally isothermal specimen model.

Preparation bounds pair count, shared original-law rows, retained six-channel
spectral storage, source resolutions and four coherent solves per angular/
spectral node. Each body caches both native T and h. Changed original optical
identities rebuild the common response; fixed optics reuse it as temperatures
change, and unchanged states/baths reuse their Planck rates. The native source
hook prepares both external sources and internal transfers once for the same
accepted start state. Positive frozen-optics secants from pair and bath channels
combine with the original conductive/source timestep limit. Nonlinear optical
state laws and time accuracy still require timestep refinement.

Accepted-start watts and actual accepted duration own the interval before
endpoint refresh. Per-pair and per-cell energy maps report that same applied
heat, with exterior throughput separate and internal cancellation explicit.
Clones share immutable laws/responses and keep candidate state/energy private;
domain, work, coverage or late-cancellation refusal leaves both accepted owners
unchanged. The native tests include an independent unreduced gap-intensity and
SI-Planck reference with two asymmetric evolving film stacks, unequal caloric
masses, ordinary conduction and three timestep refinements against RK4. Further
focused checks cover transparent throughput, two-face coverage and contact
events, original ambient events, split continuation, fixed-optics reuse, shared
enthalpy laws on a latent plateau, canonical pairs and atomic refusal.
<!-- native-coating-transmitting-pair:end -->

<!-- coating-radiative-stack:start -->
### Mutual screening and exchange across multiple coated bodies

`CoatingRadiativeStackResponse` prepares one through sixteen ordered,
vacuum/vacuum, internally isothermal reciprocal coating bodies. Every selected
emitter port faces LEFT. Each body retains its original coherent film optics;
inter-body gaps carry phase-averaged far-field powers. Two original coherent
face solves per body/wavelength/cosine serve all `N*(N-1)/2` body pairs, `2*N`
body/bath channels and one exterior-throughput channel. Angular and polarization
modes are matched before integration, not combined from independent averages.

A positive normalized prefix-source elimination and backward sweep carry escape
independently of rounded unit reflection, retain representable weak coupling,
verify reciprocity and handle source-free trapped gaps without inventing fields.
Optical preparation costs O(N^2) per matched mode, with bounded scratch storage.
Only common source-band midpoint spectral/angular integrals are retained in
shared immutable storage. Temperature queries allocate nothing and use the
original stable Planck differences and secants. Canonical i<j transfers apply
equal/opposite watts; direct bath throughput heats no body. Equal temperatures
give exact zero exchange, including for near-trapped modes.

Preparation enforces hard and caller body bounds, retained-channel and optical
solve budgets, full-band vacuum/isothermal/common-band admission and bounded
cancellation. No partial response is published. Six G1/G3/G4 groups check
independent dense directional-intensity balances for 1, 2, 3, 5 and 16 bodies,
original coherent three-body spectra and conservation, reversal, single/two-body
reductions, opaque shielding, transparent limits, sub-roundoff positive loss,
stable weak Planck work and admission/cancellation. Model scope follows Byrnes,
arXiv:1603.02720v5 section 6, and the existing local Kirchhoff emission owner.
The caller must justify phase averaging; no finite-edge, rough, near-field,
coherent-gap, internal-nonisothermal or bolometric extrapolation is implied.
<!-- coating-radiative-stack:end -->

<!-- native-coating-radiative-stack:start -->
### Native thermal evolution of ordered radiative shields and coated bodies

`heating::coating_exchange::radiative_stack::NativeCoatingRadiativeStack`
connects one through sixteen independently stored isothermal bodies in each
`CoatingRadiativeStack`. Bodies remain in authored physical left-to-right order;
every original emitter's selected port is its LEFT face. Only assembly IDs are
canonicalized. The common finite source band, vacuum on both sides of each body,
parallel/coextensive geometry and phase-averaged propagating gaps are explicit.
Coherent film interference remains inside the original optical owners.

One complete `CoatingRadiativeStackResponse` includes screening and repeated
reflections through all bodies. All body/body channels, including nonadjacent
bodies reached through transmitting shields, enter conservative native
`HeatTransfer`s. Body/bath channels alone enter external thermal heat. Direct
bath throughput stays separate and cannot heat transparent bodies. Each body
keeps its original mass, geometry, material and caloric/storage identity, with
independent temperature/enthalpy optical laws. Native ambient baths follow the
first and last bodies' exterior boundary temperatures, including events.

Both faces charge `2 * area_m2` against each body's current exposed area;
coverage is rechecked after native contact and boundary events. Bodies within an
assembly need distinct receivers. Different assemblies may share receivers on
explicitly disjoint patches, subject to their summed area. This owner neither
proves geometric non-overlap nor detects another owner on the same ports/band.
It adds no contact, finite-edge view factor, near-field exchange, coherent gap
phase, unsupported spectral tails or internal nonisothermal-film radiation.

Admission bounds assembly count, total body count, shared original optical-law
rows, retained channel nodes, source evaluations and original face solves before
copying or solving the model. Each assembly retains
`N_lambda * (N_body + 1) * (N_body + 2) / 2` spectral slots and needs
`2 * N_body * N_lambda * N_mu` coherent solves when its optics change. Body T/h
caches preserve original state; immutable optical responses are shared across
clones and reused when optical identities are unchanged. Source and transfer
hooks reuse one prepared native start state. Positive frozen-optics conductances
compose with native conduction limits; nonlinear state laws require timestep
refinement rather than an unconditional stability claim.

Accepted-start watts times actual accepted duration enter per-assembly pair,
exterior-throughput and per-cell energy maps before endpoint refresh. Internal
heat cancels in the whole-system account. Native candidate ownership preserves
accepted states on domain, coverage, budget or cancellation refusal. Focused
G1/G3 tests compare a three-body asymmetric infrared assembly to an independent
unreduced gap-intensity solve and thermal RK4 refinement with unequal masses and
ordinary conduction. Further tests exercise transparent intermediates, opaque
screening, bath events, continuation, fixed-optics reuse, enthalpy-driven optics
on a latent plateau, coverage and atomic refusal.
<!-- native-coating-radiative-stack:end -->

<!-- material-temperature-conductor:start -->
### Original material claims drive native infrared absorption and cooling

With `material-state` and `optical-heating`,
`SurfaceConductorLaw::try_material_temperature` binds an immutable material card
and its complete `SpectralConductorStatePoint` to one explicitly selected static
absolute-temperature coordinate. The original native caloric owner supplies T;
all other coordinates and quantity conventions remain fixed. Each evaluation
refreshes the original selected claims, including source wavelength properties,
all n/k knots, density and auxiliary properties. It does not choose replacement
claims, fit a new temperature table, infer phase mixing or widen wavelength
support. The original material source remains available through `material_source`.

The existing material optical binding admits the refreshed knots into the same
native conductor used by surface illumination, coated substrate absorption,
`NativeCoatingEmission` and native coating radiative exchange. Absorption uses
the original normal electromagnetic power flux. Reciprocal thermal emission
uses that absorption and the same cell temperature through the original Planck
owner. Only absorbed incident power enters heat; signed finite-band radiative
exchange is debited once from native storage. Explicit coatings continue to use
their original film laws and independent receivers.

This is an explicitly opaque, isothermal, semi-infinite substrate model. Its
native admission requires n and k strictly positive at every source knot;
exactly transparent knots refuse, including when transparency first occurs at a
heated endpoint. Positive k alone does not justify treating a finite specimen
as opaque; thickness-dependent transmission and deposition remain with the
existing finite-film and volume owners. No roughness law, penetration depth,
bolometric tail or thermal material/card association is inferred.

All retained property queries count against the source-property limit and the
existing shared optical-law preparation budget. Laws and immutable cards share
storage across native candidates. Native endpoint validation rechecks complete
source support and native optical admissibility even in darkness. Every retained
bare/coated incident wavelength is checked against the freshly admitted source
band without another source or Fresnel solve, so moving support, late auxiliary
expiry or loss of opacity refuses the whole candidate. Caller T
bounds are additional restrictions, not certification of an intervening path.
The existing accepted-start feedback remains explicit and needs time refinement.
Five focused regression groups cover source pin retention and nonaffine T laws,
analytic Fresnel heating with matching Planck emission, native IR bath heating,
cooling and equilibrium, strict transparent limits, budgets and dark rollback.
<!-- material-temperature-conductor:end -->

<!-- material-temperature-coating:start -->
### Original material claims drive transparent and absorbing finite films

With `material-state` and `optical-heating`,
`SurfaceCoatingLaw::try_material_temperature` binds the complete original
`SpectralConductorStatePoint` n/k schema to an explicitly typed native absolute
temperature and fixed physical film thickness. This source schema admits k=0;
the existing finite-film owner preserves exact transparency and transmission,
including when a film becomes absorbing after external heating. Each refresh
requeries the original selected wavelength, n/k and auxiliary claims together.
Neither a new temperature table nor a replacement claim or alpha-to-k fit is
inferred. Other source coordinates and quantity meanings remain fixed.

The original `CoherentAbsorbingCoating` supplies interference, positive layer
absorption and reflected/transmitted power in opaque or transmitting stacks.
Native heat receives only the film's actual absorbed watts. The same law feeds
`NativeCoatingEmission` and the existing coating exchange families, retaining
directional/spectral reciprocity and original per-film temperatures. Transparent
films emit and absorb exactly zero in this model; transmitted power never
becomes a thermal credit. Complete material-property counts use the existing
shared-law budgets, and immutable card/source data share storage across clones.

Endpoint validation refreshes source support even during dark intervals and
checks every bound incident wavelength against the refreshed film before any
thermal candidate is published. Auxiliary source expiry, moving bands, invalid
native constants or cancellation preserve the original accepted state. Caller
temperature bounds restrict evaluation without certifying an intervening path.
Five focused tests cover original claim pins and nonlinear temperature laws,
transparent-to-absorbing heating against independent two-interface slab powers,
matching SI Planck emission, infrared bath heating/cooling/equilibrium, complete
property budgets and dark opaque/transmitting rollback.

This retains the native finite, coherent, locally planar, isotropic film model
and its thickness/index limits. A film is isothermal in its assigned native
cell; thermalization is prompt and accepted-start optical feedback is explicit.
No thermal expansion, arbitrary depth profile, phase-mixture rule, bulk thermal
card association, unsupported spectral tail or new caloric integrator is added.
<!-- material-temperature-coating:end -->

### Original dielectric claims drive native slab opacity and thermal exchange

With `material-state` and `optical-heating`, `MaterialSlabLaw` retains an
original tabulated dielectric or analytic Sellmeier source and its independent
spectral absorption. Each accepted native temperature refreshes all selected
claims, including auxiliary properties, wavelength coordinates and fit limits.
The existing material optical binding admits the resulting glass; other state
coordinates, claim selection and quantity conventions remain fixed. Declared
temperature bounds further restrict queries without certifying source coverage.

`MaterialSlabLaw::try_schott_with_reference_absorption` selects the original
shared SCHOTT law instead of refreshing phase-index claims at each temperature.
Its spectral source stays at T0, including explicitly frozen absorption. A
separate accepted mechanical reference supplies the exact card, complete query
and selected physical requirements; only those requirements are re-resolved
at the native receiver's actual temperature. Optical claims may therefore be
valid only at T0 while the independent thermal law supports a wider interval.
No supplied mechanical requirement is dropped, and its original domain remains
binding. The original reference phase is corrected once through
`MaterialThermoOpticBinding`, preserving thermal cycles and independent source
identities. Optical plus mechanical requirement counts share the existing
distinct-law budget. This mode reaches the same native emission, radiative
exchange and material-slab illumination below.

`MaterialSlabLaw::try_schott_with_absorption_source` instead retains an
independent `SpectralDielectricAbsorptionStatePoint` from the same reference
condition. Only its original wavelength, intensity-absorption and auxiliary
claims are refreshed at each accepted temperature; reference phase claims
remain at T0. `MaterialThermoOpticBinding::at_material_state_with_absorption`
composes that current spectrum with the original thermal phase law into one
native material binding. Complete card/query/quantity agreement is required.
The nine original loss knots retain their interpolation and finite band,
including for an all-zero spectrum; phase and absorption support intersect.
Native Beer attenuation, slab illumination, emission and radiative exchange
therefore use the same sourced alpha(T) and n(T), without a second loss fit or
thermal debit. All three source bundles count toward the distinct-law property
budget, and both original loss policy and current source identity distinguish
cache entries. Wavelength queries reuse inline prepared absorption storage.

`StatefulMaterialSlabIlluminationSource` enters `SurfaceOpticalHeating` through
`CoatingSurfaceSource::MaterialIncoherentSlab`. One isothermal native cell
supplies the accepted temperature and receives the complete positive slab
absorption. The original source law updates phase and absorption together; SCHOTT mode
retains reference phase claims and uses either explicitly frozen loss or
independently refreshed current absorption. Both real-index interfaces, refracted Beer path lengths and the
complete incoherent internal-reflection series use that same glass snapshot.
Independent s/p powers are summed before polarization mixing. Exterior media,
positive normal thickness and interception geometry stay fixed; reflected and
transmitted watts never become deposited heat. The same owner admits mixed
coherent-coating and incoherent-slab sources without adding an optical interface.

Original material laws, receiver/law resolution, identical response packets,
source budgets, pulse transitions, energy accounts and transactional native
advancement reuse the existing coating observer. Wavelength and temperature
support remain binding at accepted endpoints and during darkness, including
independent absorption support. This is accepted-start optical feedback and
requires temporal refinement. It does not infer a homogeneous index from
different depth-cell temperatures, add spatial thermal-lens transport, move
geometry or duplicate separately composed radiative cooling.

`NativeCoatingEmission::try_material_incoherent_slab` and
`IsothermalCoatingStack::MaterialIncoherentSlab` feed the existing finite-band
emission and native exchange owners. Material resolution occurs outside optical
quadrature. Complete source-property counts share distinct-law budgets; accepted
source identity governs optical-cache reuse. Endpoint refusal publishes neither
candidate heat nor candidate emission/cache state.

The original weak-absorption, real-interface, incoherent slab model remains
binding. Its geometry and lossless opposite medium are explicitly fixed; one
native receiver owns its isothermal storage. Only true absorption exchanges heat,
and transparent slabs emit zero. No temperature/emissivity fit, phase mixture,
thermal expansion, unsupported spectral continuation or new caloric solver is
introduced. Accepted-start thermal feedback remains explicit and requires time
refinement.


## Source-bounded quantum-efficiency detector channels

`spectral::detector::QuantumEfficiencySensor<C>` retains 1–16 external
quantum-efficiency channels on 2–1024 original positive vacuum-wavelength
knots. Efficiencies are collected electrons per incident photon without gain,
between zero and one. Piecewise-linear source curves, dark channels and dark
spectral gaps remain explicit; there is no extrapolation or CIE/RGB extension.
An analytic sum-of-QE proposal costs O(C*K) to prepare and O(log K + C) per
wavelength with no query allocation or rejection sampling. Every draw supplies
the physical QE/pdf radiance weight and QE*lambda_vac/(h*c*pdf) electron weight,
using the exact SI defining constants. Zero-density endpoints use the analytic
one-sided ratio. Unrepresentable nonzero support or response refuses instead of
silently becoming a dark channel.

`CurvedCameraRaster::sample_pixel_sensor` integrates those original channels
through four wavelength-specific native camera paths. A stratified CDF packet
shares four paths across every channel; channel count does not multiply scene
transport. Curved geometry, wavelength-dependent index and absorption, source
emission, surface events, path budgets and absolute sample identities remain
owned by the original camera. The common raster/lens helper preserves the
existing visible XYZ path and random-counter domains.
`sample_pixel_sensor_rayleigh` also retains each wavelength's original
scattering spectrum, phase law, free-flight draws and volume-event budget.
Zero scattering preserves the complete surface detector sample bit for bit.
Failure or cancellation publishes no partial channels or renormalized lanes.

`sample_pixel_sensor_area_lights` and `sample_pixel_sensor_rayleigh_area_lights`
observe finite-area emission through the same native source connections,
complementary emitter-hit MIS and independent source-survival tracking as the
color observer. Their original detector proposal supplies each wavelength;
all channels share four paths and receive their own physical QE and photon
weights only after complete path radiance is assembled. Generic area-pixel
results retain the detector sample and original per-lane work accounting.
Unsupported source bands, exhausted work or cancellation refuse the complete
observation. No CIE conversion or additional light/survival weight is applied.

Outputs are QE-weighted radiance in W/(m² sr) and collected-electron rate per
m² sr, before pixel etendue, exposure and electronics. The existing aperture
average is retained. These observations do not infer sensor irradiance,
calibrated counts, gain, shot noise, saturation, or display colors. Authored
response curves are not empirical detector calibration.

Focused native detector tests exercise analytic energy/photon moments, sloped
and disjoint responses, weak-response scaling and range refusals. Camera
regressions cover source-bounded infrared emission, independent dispersive GRIN
paths and absorption, deterministic replay, dark channels, molecular transport
parity, volume limits and cancellation. The detector numerical module was
executed against real fs_math; full native camera/material integration tests
remain unexecuted without the complete constellation dependency checkouts.

## Analytic spectral material bindings

`MaterialOpticalBinding::try_spectral_sellmeier` carries the original finite-band
material fit directly into `SellmeierIor::try_new_in_band`, converting source
m² to native µm² and vacuum metres to nm exactly once. Pole and transparent-index
admission covers the entire original fitted interval before intersection with
the separately retained Beer intensity-absorption support. No sampled phase-index
approximation replaces the analytic law or its group derivative. Native
Snell/Fresnel, phase/group path and Beer consumers share that same accepted law.

Bindings require matching complete physical state and retain explicit surface
state. Refresh re-evaluates pinned source claims, readmits the original family
and reconstructs the binding atomically. Focused controls include high-index
infrared sources, full-fit pole refusal outside absorption support, original
claim identity and source-state evolution.


### Native diffuse thermal surfaces on spectral camera paths

`CurvedCameraScene::with_diffuse_thermal_emitters` binds immutable
`DiffuseThermalEmitter` values to original native geometry. Each value retains
one `SpectralDiffuseThermalSurface`; its `reflectance_spectrum()` supplies the
same finite-band opaque Kirchhoff law to the scene's depolarizing diffuse
interaction. Binding checks the complete numerical knot/value identity, or
admits an opaque terminating surface only for emissivity identically one.
Unknown geometry, duplicate thermal bindings and calibrated emitters on the
same object refuse. Both physical faces follow the original thermal owner's
per-side convention. Actual emitting contacts require phase index exactly one,
so a vacuum Planck spectrum is never inferred for an immersed surface.

The existing emission accumulator evaluates the original Planck/emissivity
law at each actual path wavelength before reflection. It preserves source-band
refusals, reciprocal radiance conversion and native scattering weights. QE
pixels therefore observe thermal infrared energy and wavelength-correct
photoelectron-rate integrals without tabulating or fitting Planck radiance.
Preparation stores one bounded indexed source table; per-ray queries allocate
nothing. Observation does not advance, cool or mutate the original thermal
state. These sources currently use emitter-hit transport; no thermal area-light
proposal or additional cooling owner is inferred. Detector surface transport
remains available without `volumes`; only the scattering-specific detector
APIs and corresponding tests require that feature.

## Broadband molecular refraction through native GRIN transport

The original material gas-composition field consumes explicitly supported
broadband molecular Cauchy sources through its existing prepared-wavelength path.
Phase index, true group index, composition gradient and Hessian keep the same
Lorentz--Lorenz mixture; no phase-index table, extra gas law or per-stage spectral
fit is introduced. Material wavelength refusal remains authoritative and a
prepared field admits only its retained wavelength. Optical temperature/density
and composition continue to come from the original frozen reconstruction.

Focused native controls cover 300, 1550 and 2500 nm phase/group paths and analytic
ray-bundle footprints, with independent material scattering and source-domain
controls. These are explicit synthetic fits, not measured-gas accuracy claims.
Molecular scattering stays distinct from true absorption.

## Polarized anisotropic dielectric reflection and transmission

`PolarizedBeam::sample_interface` and `evaluate_dielectric_connection`, including their original borrowed-stack variants, consume the same scalar anisotropic GGX distribution, directional masking, visible-normal proposal and refractive measure. Raw beam queries receive the authored tangent in the beam's world/query frame. Straight and curved native scenes bind the original geometry-local axis through the owning instance placement before scattering; the shared continuous camera receiver applies that identical binding before source connections and recorded-PDF replay. Scene snapshots remain immutable during a trace.

The finish changes facet statistics. Every facet retains the original isotropic bulk/coating Fresnel Mueller law, its actual s/p plane and coherent relative phase; no depolarization or birefringent material response is invented. Film order, survival, adjacent phase indices and forward-power versus radiance conversion retain their existing owners. Coated original stacks retain anisotropy instead of reducing to alpha-x. Isotropic surfaces bypass axis binding and preserve their former numerical path; invalid contact axes refuse.

Five focused tests pass against the original numerical owners, checking independently reconstructed Walter densities and full Stokes power, normal-facet analytic power/probability, Brewster extinction with positive proposal density, reciprocal source/frame changes, original-film equivalence and reverse incidence, and isotropic/null limits. Three native scene/receiver regressions cover body rotations, both boundary sides and original stack paths, plus consistent camera receiver/MIS densities. The actual renderer and unit-test code pass a real-dependency Cargo metadata check with tracer, volumes and photometry enabled; the scene tests remain unexecuted. A wider original-source optical/beam build passes 189 of 190 tests: the existing deep-film sampler-cache test fails identically on the pre-edit source. This remains allocation-free local single-scattering transport; no inter-facet multiple scattering, coherent interference between scene rays, measured finish calibration, or bulk anisotropic dielectric tensor is claimed.

## Geometric thin-lens refraction in an explicit exterior medium

`DispersiveThinLens::try_new_in_medium` retains one native absolute phase-index law for the homogeneous medium on both sides of coincident paraxial thin elements. Optical power uses the lens/exterior index contrast; double-double cross-product differences preserve weak chromatic changes. Every element and the exterior keep their original source band and table-gap refusals. Native raster transport must return the same exterior phase index at each sampled vacuum wavelength. Opening-scale admission uses the wavelength in the medium, while detector photon conversion retains vacuum wavelength. Preparation is bounded to eight elements and spectral ray evaluation allocates nothing.

This is an ideal transparent, unit-transmission geometric lens; existing scene transport still owns exterior absorption, scattering and radiance conversion. No thick-lens geometry, surface aberrations, coating/Fresnel lens losses, or unprovided source provenance is inferred. An explicitly declared exterior supports the finite-defocus and actual-source-depth layered wave policies in the same homogeneous medium; the in-focus policy remains unsupported for a dispersive lens. Focused actual-source controls cover the paraxial limit of exact two-surface Snell refraction, dispersive pupils, common-index scaling, cancellation-sensitive color, source/domain and medium mismatch refusals, and native geometry-to-Beer-to-QE/photoelectron evaluation.

## Transient P1 conduction snapshots in native thermo-optic transport

`material_state::thermo_optic::transient::P1ConductionTemperatureSnapshot` borrows the original `ConductionMesh` and an accepted temperature-based nodal state. Snapshot admission checks all vertices once without copying or allocating; each prepared tetrahedron reuses its original barycentric gradients in constant work and storage. Temperature and its gradient retain the original one-sided P1 reconstruction, including duplicated contact-node discontinuities. The thermal Hessian is zero inside an element, while `MaterialThermoOpticField` still supplies the nonlinear optical term `n_TT grad(T) grad(T)^T` from the original material law and caller-retained snapshot identity. Prepared sampling allocates nothing.

`ConductionElementBoundary` reuses the native interval triangle/edge proximity owner for the four original outward-oriented faces. Curved transport can reach the first face without skipping interpolation or contact discontinuities. A conservative finite-face radius guards shallow incident tails; it does not authorize departure. Edges, competing faces and roundoff-ambiguous contacts refuse. The adapter supplies no global locator, automatic interelement crossing, reflected/transmitted departure, geometry deformation, temporal interpolation, or substitution of nodal temperatures for a pointwise enthalpy reconstruction. Other scene surfaces remain the responsibility of the original scene geometry owner.

Five focused actual-source controls verify both tetrahedral orientations, original contact jumps, admission/cancellation and finite-face/departure refusals, original backward-Euler cooling through analytic curved-ray phase/group paths and Beer absorption, and nonlinear optical Hessians through native bundle derivatives checked against independent primal rays. These controls use the original finite-element, conduction, material, interval and optical numerical owners.

## Continuous thermo-optic rays across the original P1 thermal mesh

`thermo_optic::transient::mesh::P1ConductionMesh` prepares adjacency from the original vertex-ID faces and propagates primal rays under one original reference glass, thermo-optic law and accepted nodal-temperature snapshot. Shared nodal values give continuous temperature/index while each side retains its original gradient. A private two-cell departure is rebound to the native arrival witness; interval proximity excludes every other face of both tetrahedra. The existing event integrator accounts for the complete transition collar without moving the endpoint or dropping arc length, phase/group path or Beer optical depth. Wavelength preparation is shared, propagation allocates nothing, and compensated segment accumulation and whole-path step/query/face/cancellation budgets remain explicit.

The caller supplies the containing launch element and the intended conforming mesh. Exterior and duplicated contact-node faces stop at native events; no material-index jump, other scene surface, global locator or global embedding proof is inferred. The mesh also supports native launch-angle bundles through the continuous-index face-coordinate continuation described below, including the crossing-time variation caused by gradient jumps. Six new actual-source controls plus five original transient controls pass, including a twelve-tetrahedron column in both orientations, independent piecewise-affine ray solutions, partial-collar endpoints, original two-cell backward-Euler cooling, narrow-neighbor/contact refusals and global budgets.

## Original water state in shared dielectric and gradient-index transport

With `material-state`, `DielectricIor::Water(IapwsWaterIor)` retains the original IAPWS R9-97 phase/group law at the supplied temperature and mass density, without a Cauchy fit or wavelength table. Native Fresnel/Snell, total internal reflection, dispersion and homogeneous group delay use the existing optical owners. Material, camera, coating, thermal-coating and Euler scene identities retain the released law and exact frozen state; old variants retain their identity bytes.

`material_state::water_refraction::MaterialWaterRefractiveField` borrows one immutable ordinary-water T/rho reconstruction and an explicit independent spectral absorption law. Preparation caches wavelength terms and intersects source domains. Analytic phase gradients and Hessians include both source gradients/Hessians, both pure second partials and both mixed outer products; the separate group index supplies narrowband delay. Prepared ray and bundle queries allocate nothing and reject foreign wavelengths. A bounded provider cannot acquire global uniformity from a locally zero gradient.

Eight focused actual-source controls pass for native Snell/Fresnel/TIR, phase/group/Beer paths, state/material identities, source/domain refusals, full three-dimensional finite-difference derivatives, conserved tangential optical momentum and density-lens bundle footprints against analytic solutions. The caller retains homogeneous ordinary-water purity, phase validity, state correspondence, spatial smoothness and boundary ownership. No EOS, phase-mixture interpolation, absorption coefficient, or arbitrary mesh crossing is inferred.

## Immersed finite-defocus and source-depth wave optics

The finite-defocus and layered camera owners now use the explicit exterior's physical wavelength `lambda_vac/n_exterior` for angular scale, quadratic pupil phase, spectral sampling Jacobians and the strongest depth-error gate. Relative lens power remains the sum of `(n_lens/n_exterior - 1)` times each original curvature difference. Native material, source and sensor queries retain their original vacuum wavelength; a synthetic proposal coordinate never becomes a source query. The original bounded circular, annular and polygonal pupil kernels are reused without per-sample allocation.

Cinematic camera continuity includes the complete exterior phase law and frozen state while retaining legacy index-one identity bytes. Launch admission compares each original wavelength with the phase index of the actual original medium stack. Thirty focused actual-source controls pass: twelve lens/wave controls (including independent complex pupil quadrature, integrated radiance, actual source depth, nonmonotone exterior dispersion, source-coordinate retention and the native IAPWS water state), two exact camera-identity controls and sixteen existing camera regressions.

The model retains one homogeneous medium on both sides of up to eight coincident paraxial ideal thin elements, scalar incoherent propagation, clear pupil and pupil-independent radiance/visibility. It does not add thick-lens or Fresnel lens losses, high-NA/vector optics, coherent scene-point interference or clipped-pupil reconstruction. The existing general cinematic scene admission still refuses unsupported enclosing dielectric geometry, and CurvedCameraRaster still refuses wave imaging. The positive capability is the native sampling, density and actual-source-resolution API under these existing optical conditions.

## Water thermal lenses from original isobaric IF97 states

With `material-state`, `material_state::water_refraction::isobaric::IsobaricWaterField` borrows an existing `TemperatureGradientField` and, when available, its original `TemperatureSecondOrderField` capability. The caller explicitly chooses one fixed SI pressure, IF97 Region 1 or 2, pure ordinary-water applicability and an accepted snapshot identity. Every query evaluates the original EOS at that field's actual temperature. Analytic density gradients and Hessians retain rho_T * grad(T) and rho_T * Hess(T) + rho_TT * grad(T) outer grad(T), including all mixed spatial entries.

`bind_refraction` supplies those original states to the existing prepared IAPWS water optical owner and an independent absorption law. Phase/group transport and ray/bundle integration stay with their existing owners. Source/pressure/branch/snapshot identity, bounded-field behavior and exact EOS/optical source-domain refusals remain explicit. Preparation and queries allocate nothing; only the second-order path requests the additional Gibbs third partial.

Four focused actual-source controls cover the complete three-dimensional EOS-plus-optics derivative chain against independent scalar differences, native bundle focusing against an analytic on-axis solution, a solved conduction profile driving native curved water transport with conserved tangential optical momentum, and source/phase/spectral-domain refusals. No pressure field, phase-mixture interpolation, moving-fluid/Doppler optics, free-surface geometry or mesh-face crossing is inferred. A saturation or regional crossing is refused by the original selected branch; the adapter does not widen IF97 or IAPWS optical uncertainty.

## Prepared immersed wave images through the native scene and sensor

`polarization::dielectric::scene::curved::camera::raster::CurvedWaveCameraRaster` binds the existing finite-defocus/layered diffraction camera to its original `CurvedCameraScene` and an explicit clear homogeneous initial medium. The original medium's full phase law must match the camera exterior; equality at one reference wavelength is insufficient. Implicit-vacuum cameras require a nondispersive index-one medium. Original rigid source meshes must be opaque, lie on the declared finite axial planes, and pass the existing interval `ProjectedPupilFootprint` full-pupil separation proof.

Each wavelength traces and retains its own original path, emission and first-contact witness. Geometric axial depth comes from that witness, never from n*distance or group optical length. `LayeredCameraRay::resolve_spectral_sources` reuses the original proposal and pupil kernel for separate spectral first hits and black misses. Native XYZ and multichannel QE/photoelectron accumulators apply each optical weight once; exposed spectral samples remain the original unweighted transport radiances, and photon energy retains vacuum wavelength. The observer adds no per-sample workspace allocation or alternate ray integrator.

Eighteen focused actual-source checks pass: eight new physical/admission controls and ten existing raster/sensor/geometric-dispersion regressions. They include positive IAPWS Water imaging at 800–1000 nm, independent phase/group lengths and photon energy, the vacuum limit, integrated uniform-source normalization, per-lane layered depths, pupil-edge occlusion, and source/medium/budget/cancellation refusals. The generic geometric raster still rejects diffraction. Arbitrary source geometry, refracting or reflecting primary surfaces, GRIN propagation, absorbing-pupil apodization, environment/volume emission, unresolved clipping and coherent scene interference remain outside this prepared scalar paraxial incoherent model.

## Native thermal mesh ray bundles across continuous-index faces

`material_state::thermo_optic::transient::mesh::P1ConductionMesh::propagate_bundle` now carries the native launch-angle family through the same accepted P1 thermal snapshot and original Schott thermo-optic material. The original one-sided temperature gradient and zero in-cell temperature Hessian feed the complete nonlinear material chain rule; the refractive-index Hessian is not assumed zero. Wavelength preparation and the original angular basis are retained across cells, along with phase/group optical length, Beer optical depth and receiving-plane footprint/caustic admission.

Inside each original planar face collar, a temporary fixed-plane-coordinate representation carries both variational columns and their crossing-time derivatives through the actual central trajectory. This recovers the continuous-index saltation at the gradient jump without applying it prematurely at the approximate arrival point. Endpoints on either side of the finite collar are converted to derivatives at the same total geometric arc length. The native Dormand–Prince stages and embedded error estimates include the crossing-time columns and their effect on physical position/direction errors; ordinary smooth propagation retains its original fast path. No finite-difference bundles, neighbor-ray pretraces, smoothed gradients, endpoint offsets or per-ray workspace allocations are introduced.

The prepared mesh and material must retain the original mesh/nodal-buffer snapshot identity and one continuous material phase law. Each stage requires the declared transverse approach bound. A roundoff-unresolved actual endpoint on a genuine temperature-gradient-jump plane refuses the derivative. Native exterior and duplicated contact-node faces still terminate. Compensated distance accumulation, whole-path step/face/query budgets and cancellation remain explicit. Local adaptive ODE/variation estimates are not a global accuracy certificate; this adds no discontinuous-index material interface, curved-face continuation, moving thermoelastic mesh or global geometry/point-location proof.

Forty-two focused actual-source tests pass: seven new bundle controls, eleven retained thermal/transient controls and twenty-four original native bundle regressions. They cover many-cell footprints in both orientations, independent exact piecewise-affine ray references, equal-arc endpoints on both sides inside the collar, nonlinear material curvature and distinct group path, an accepted backward-Euler thermal state, contact handling, original source identity and global budgets. In the genuine-gradient-jump refinement control, the maximum position/direction column error against independent neighboring analytic rays falls from 2.1625e-5 to 5.0302e-9 to 9.0960e-11 over three step/tolerance settings.

## Water refraction from joint pressure and temperature fields

`material_state::water_refraction::pressure_temperature::PressureTemperatureWaterField` borrows one immutable joint reconstruction of absolute thermodynamic pressure and temperature. Each query evaluates the selected original IF97 Region 1/2 state, retaining compressibility, expansion and their mixed response. The spatial density Hessian contains both field Hessians, both pure gradient outer products and both mixed pressure/temperature outer products. Ordinary ray queries reuse the stored first EOS responses; only bundle queries request the full analytic density jet. Preparation and sampling allocate nothing.

`bind_refraction` uses the existing IAPWS phase/group law, absorption law and native ray/bundle transport. Focused controls cover liquid/vapor scalar finite differences through the complete EOS-plus-optics chain, exact constant-pressure compatibility, pressure-driven bending with conserved tangential optical momentum, and analytic pressure-lens focusing with distinct phase/group path lengths. Snapshot/branch identity, bounded-field behavior, nonfinite range and native EOS/optical refusals remain explicit. A frozen compressive snapshot must supply its actual pressure and temperature; the adapter supplies no acoustic thermal closure, phase-mixture tangent, moving-fluid/Doppler effect, regional continuation or mesh transition.

## Original water refraction through continuous thermal meshes

With `material-state`, `P1ConductionMesh::propagate_water` and `propagate_water_bundle` carry an original isobaric IF97/IAPWS water binding through the accepted P1 conduction snapshot. The source mesh and nodal buffer must be identical to the prepared mesh's snapshot. Pressure, selected IF97 branch, separate inner-state/outer-optical snapshot identities and independent absorption are retained. One original phase law and shared nodal temperatures give continuous refractive index across internal faces; these faces change curvature without introducing Fresnel events.

Both glass and water now call the same private prepared-field traversal owners. Water's complete nonlinear density/optical Hessian feeds the existing variational integrator even though the in-element temperature Hessian is zero. Native gradient-jump crossing-time corrections, transition collars, phase/group path lengths, Beer loss, receiving-plane projection, cancellation and whole-path work limits remain in force. Wavelength preparation is reused; propagation adds no workspace allocation. Contacts and exterior faces still stop, and EOS/optical domain crossings refuse.

Four focused controls cover original smooth-field equivalence through a twelve-tetrahedron column in both orientations and both liquid/vapor branches, gradient-jump bundles against neighboring primal rays, accepted backward-Euler cooling, and contact/phase/snapshot/work refusals. This route uses a frozen pure-water thermal state at fixed pressure; it adds no fluid closure, mixture interpolation, varying-pressure mesh reconstruction, moving-fluid optics, global locator or free-surface interface law.

## Polarized refraction at original thermal-contact faces

`material_state::thermo_optic::transient::mesh::contact::P1ContactInterface`
explicitly pairs two coincident original triangular faces in the same accepted
P1 conduction snapshot. Each side retains its own reference glass and Schott
thermo-optic law. The four original nodal values must prove each element
isothermal; distinct temperatures and material indices across the contact are
retained. Faces may share vertex IDs or use duplicated thermal-contact nodes.
Exact geometric pairing and opposite interiors are checked without smoothing
the jump, inventing an air gap or replacing a gradient by an element average.

`split` traces an initially unpolarized interior ray to the actual first face
and returns both nonzero Fresnel branches. `split_polarized` accepts a physical
Stokes state with unit input intensity and the declared transverse reference.
The existing dielectric owner supplies Snell directions, power-normalized
Mueller matrices and the relative reflection phase under total internal
reflection. Each contact rotates the actual incident frame into its local s/p
basis. The complete Stokes state survives subsequent contacts; unpolarized
scalar Fresnel averages remain diagnostic coefficients, not cumulative power
weights. Returned Stokes components and `power_fraction()` are relative to the
original unit launch and may be scaled by the caller's incident spectral power.

Each branch propagates within its original receiving element to the next
native face or the explicit outgoing-distance limit. Phase optical length,
group optical length, geometric length and Beer optical depth accumulate from
the original wavelength-dependent laws; each interface and bulk attenuation
factor is applied once. Native watertight triangle hits retain the actual
endpoint, without an unintegrated position offset or a GRIN arrival collar.
`split_arrival` composes a caller-selected retained branch with the next
explicitly paired contact, checking its private original-owner/face witness.
The API supplies no implicit branching tree or unbounded reflection loop.

`with_coating_stack` binds an original borrowed `DielectricCoatingStack` in
side-0-to-side-1 order. Reverse incidence reverses the physical film traversal
and uses the actual reversed thermal indices and Snell angle. Transparent and
absorbing films retain their original wavelength domains, dispersion,
extinction, thickness and source identities, including zero-thickness films.
One native coherent solve supplies the signed Fresnel matrices and positive
dissipative film work. `coating_absorbed_power_fraction` weights that work by
the actual incident s/p powers and retains losses below the rounded R+T
residual. It refers to this contact relative to the original launch, excluding
bulk attenuation and earlier contact losses. Native diagonal-Jones application
also retains weak reflected polarization channels near Brewster incidence.

The coating is a local planar scattering boundary. Its phase affects R/T and
polarization, but the ray is not displaced across its geometric thickness and
no film group delay or heat deposition into the conduction state is inferred.
The unchanged bulk phase/group integrals therefore describe the two original
elements, not pulse transit or resonant residence inside the films.

Preparation stores eight native one-triangle charts. A split uses bounded
work, at most twelve triangle-intersection queries, no query allocation and
the original cancellation checkpoints. Optional coating work is linear in the
native limit of 64 films; bare interfaces retain their direct Fresnel path.
Edges, ambiguous faces, foreign
snapshots/materials, unsupported wavelengths, invalid launches and nonuniform
elements refuse. The model is smooth-interface polarized geometrical power
transport: it adds no coherent interference between paths, roughness,
discontinuous-index GRIN continuation, angular-bundle derivatives or global
mesh locator. Thermal contact resistance is not assigned optical absorption;
bulk absorption remains with each original glass.

### Curved refraction at nonuniform thermal contacts

`contact::curved::P1CurvedContactInterface` accepts the same original paired
faces without requiring uniform nodal temperatures. Incoming and both nonzero
outgoing branches use the native clearance-controlled GRIN event integrator.
Each optical law retains its original reference glass, thermo-optic coefficients,
spectral support and accepted P1 snapshot. Temperature and gradient are evaluated
at every optical stage; phase, group delay, Beer loss and parallel-transported
polarization follow the actual curved arc. The exact and curved contact owners
share one Fresnel/coherent-coating response, including positive film absorption.

A private numerical collar continues the selected element's original affine P1
polynomial through a thin contact slab inside a ball excluding the face edges
and every other face of both cells. The public temperature element keeps its
strict domain. No field average, position projection or unintegrated offset is
introduced. Temperature/source-domain refusals remain active inside the collar.
This follows the native numerical-contact model: the result reports
`contact_resolution_m`, not an exact interface location or a certified global
optical-path error. Tightening the event tolerance controls this localization.

Both attempted ODE steps (including rejection) and boundary queries share one
request budget across arrival, collar preparation and both outgoing branches.
Queries allocate nothing. A private original-binding/native-witness branch can
continue through `split_arrival`; its existing history is retained without
charging it again. All four faces return to the boundary search after departure,
so a gradient can bend a ray back to the same interface. The model remains
isotropic geometrical power transport on one immutable accepted thermal state;
it adds no birefringence, moving-interface physics or film transit delay.

Six focused controls cover independent Snell/Fresnel and distinct phase/group
paths, Beer loss, temperature-jump refraction and reverse incidence, TIR and
equal-index limits, composition through three original elements with nonparallel
incidence planes, analytic TIR polarization phase, frame covariance and
source/domain/cancellation refusals.

## Original P1 pressure and temperature as one water-optical source

`thermo_optic::transient::mesh::water::pressure_temperature` supplies an immutable `P1PressureSnapshot` and `P1PressureTemperatureSnapshot` on the original native mesh. Pressure is an independently supplied absolute thermodynamic field: `NodalPressureValues` borrows either dimensioned pressure quantities or explicitly declared Pa values without allocation, conversion or a pressure-as-temperature surrogate. The mesh supplies geometry and original P1 gradients only. The joint snapshot requires the exact same mesh owner and combines the accepted pressure and temperature source identities under its reconstruction identity.

Each joint element uses one original barycentric admission and the same coordinates for both fields, retaining exact vertex values and the separate original buffers. `bind_water` connects the combined identity and source to the existing IF97/IAPWS owner. Both input Hessians are zero inside a P1 element, while the complete optical Hessian retains pressure, temperature and mixed nonlinear responses. Four focused controls cover typed/raw-Pa equivalence, both orientations, liquid/vapor scalar derivative references with nonparallel source gradients, native in-element rays/bundles and source/domain refusals. This source adapter declares original nodal P1 fields; it does not turn cell averages, DG traces or gauge-pressure perturbations into an absolute pressure field or add a fluid-pressure solve.

## Water rays and bundles through joint pressure and temperature meshes

`P1ConductionMesh::propagate_pressure_temperature_water` and `propagate_pressure_temperature_water_bundle` carry the original joint P1 water source through the native face traversal. Every stage makes one original-element/barycentric admission for both pressure and temperature, including the departure collar. A two-entry pressure-polynomial cache follows that admitted element without choosing a side or interpolating between traces. The canonical joint state identity includes both accepted source identities; a mismatched inner water binding refuses, while the separate outer optical identity and independent absorption law are preserved.

The existing IF97/IAPWS jets retain compressibility, thermal expansion and the full mixed optical Hessian. Continuous original nodal pressure and temperature preserve phase index across an internal face, so the shared crossing-time/saltation transport applies without a Fresnel event. A pressure-gradient jump now participates in the native unresolved-face bundle-endpoint refusal even when the temperature gradient is unchanged. Glass and fixed-pressure water retain their existing predicate. Original phase/group paths, Beer loss, face/step/query budgets and cancellation remain under the shared transport owners; the joint route adds no per-ray workspace allocation.

Five focused controls cover liquid/vapor many-cell propagation with nonparallel source gradients in both tetrahedron orientations, pressure-only gradient-jump bundles against neighboring primal rays, exact same-side departure-collar field pairing, the native pressure-only unresolved-face guard, and source/contact/phase/work refusals. Exterior and duplicated-contact faces still terminate; a true pressure jump does not become a continuous face transition. This frozen-snapshot route supplies neither a fluid closure nor a pressure solve, moving-fluid optics, cross-phase continuation, global point locator or optical law for discontinuous contacts.

## Smooth dielectric coats over diffuse substrates

`SmoothDiffuseCoating` declares constant dimensionless indices
`n_exterior >= 1` and `n_coating >= 1` for a smooth, geometrically collapsed
dielectric above the original wavelength-resolved, depolarizing Lambertian
substrate. Its actual surrounding phase index must agree. Both reflection
branches remain in that medium; this material never changes the bulk-medium
stack or invents a thermal-emission law.

The original native Fresnel owner supplies lossless entry and exit transmission.
The clear constructor caches internal escape
`E = (n_exterior/n_coating)^2 * 2*integral(mu*T(mu), mu=0..1)`
using bounded Gauss-Legendre refinement on the lower-index hemisphere. For
`n_coating < n_exterior`, the directly integrated internal hemisphere supplies
`E`; external directions beyond the critical angle have total internal
reflection and exactly zero diffuse response. Repeated
internal returns use the cancellation-resistant
`D = (1-rho) + rho*E`, giving continuous reciprocal reflection
`rho*(n_exterior/n_coating)^2*T_in*T_out/(pi*D)`, alongside native mirror
Fresnel. This is the physical nonlinear smooth-plastic model described in
[Mitsuba's material reference](https://mitsuba.readthedocs.io/en/latest/src/generated/plugins_bsdfs.html#smooth-plastic-material-plastic).
The cached quadrature is a convergence check, not a certified enclosure or
measured-material claim. Ray queries perform no angular integration or allocation.

`with_absorption(thickness_m, absorption_per_m)` adds an explicit homogeneous
gray intensity-absorption law in the locally parallel slab. The physical
internal angle supplies `A(mu)=exp(-alpha*thickness/mu)`, including the initial
descent, final escape and both legs of each internally reflected return.
Preparation caches `E=<A*T>` and coating loss `C=<(1-A)*(1+A*R)>` in the
internal Lambertian flux measure. Returned flux is `B=<A^2*R>`; `B+E+C=1`.
The shared denominator is `D=(1-rho)+rho*(E+C)`, and the diffuse BRDF is
`rho*(n_exterior/n_coating)^2*T_in*T_out*A_in*A_out/(pi*D)`. Separate trapped
and escaping angular domains retain TIR. Bounded adaptive Gauss panels use
`expm1` for weak absorption and a factored normal transmission for nearly
opaque escape; unresolved convergence or unrepresentable optical depth refuses.
The exact zero-depth path retains the original clear calculation.

`directional_power` reports reflected power and distinct coating/substrate
absorptances. Initial descent absorption is counted before the first base hit;
repeated layer loss and substrate absorption share the original return sum.
No arithmetic balance residual is assigned to either material as a heat source.
This optical partition does not itself evolve temperature or authorize emission.
Native scalar, BDPT and polarized sampling/connections retain both Beer paths;
their proposal densities remain the original explicit mixed measures. The
material identity includes declared thickness and gray absorption coefficient,
and the substrate still uses its original reflectance at each path wavelength.

Native polarized scene/camera transport retains full entry s/p response,
substrate depolarization and polarization on escape. Scalar
`Material::CoatedDiffuse`, spectral path continuation and BDPT share the same
law, with a discrete mirror branch and unconditional continuous density
`diffuse_mass*cos(theta)/(pi*F)` on the accessible exterior cone, where
`F = min(1, (n_coating/n_exterior)^2)`. The cosine sampler scales its squared
radius by `F`; it does not waste samples outside the transmitted cone. Its
conditional weight cancels the cone's reciprocal etendue analytically, before
forming large coefficients. For the original high-index coats, `F=1` and the
original full-hemisphere law remains. Forward/reverse densities, emitting-surface
addition, material identity, actual-medium admission and existing receiver
lighting remain explicit. The smaller sampling branch lies next to zero, so
a mirror mass rounded to one cannot erase weak diffuse support. The polarized
sampler normalizes small-probability transmission before applying it to Stokes
power, retaining finite grazing and subnormal conditional power.

After spectral splitting, only active coated material responses are queried;
the original complete packet still determines every medium/source query,
free-flight proposal and MIS escape density. Calibrated jobs refuse active
point fixtures with a nonzero coating mirror because specular point-source
connections are absent. Clear matched-index coats reduce to the original
diffuse law; matched absorbing coats retain their actual Beer paths.
Intermediate analytic bundle mirror plans still refuse this mixed
material; ordinary sampled continuation and continuous bundle receivers apply.

Independent numerical controls cover flux integrals, white-furnace conservation,
internal returns, reciprocity, matched/black limits, full Stokes frames,
conditional grazing/subnormal power and mixed forward/reverse PDFs. Lower-index
coats additionally test independent internal flux integrals, external TIR for
black and white substrates, escape-cone densities, scalar/polarized agreement
and polarized white-furnace conservation. Absorbing coats add independent
internal-angle flux/return integrals, separate coating/substrate energy,
matched-index Beer limits, physical depth scaling, reciprocity, weak and opaque
limits, wavelength-specific substrate evaluation and native adapter checks.
Native scene, scalar Material and BDPT integration regressions are also present.
The model supplies no coherent interference, complex-index interface response,
roughness, volume scattering, spectral coating absorption, optical flight time
or lateral displacement. Its geometric surface remains collapsed even when a
physical thickness supplies Beer attenuation.

## Selected point-fixture transmission through ordinary scene glass

`PhotometricSource::transmitted_slab_irradiance` measures one selected image
through the two smooth parallel broad faces of a named static mesh instance in
the ordinary `tracer::Scene`. The original slab owner solves each of the 80
native wavelengths and supplies receiving area per source solid angle; the
actual source intensity, native s/p Mueller transmission history, original
surface coatings and Beer loss then supply spectral irradiance and CIE lux.
The first face's induced polarization survives the second face. Geometry and
its acceleration are prepared once; spectral queries allocate no workspace.

Body-only slab discovery cannot make unrelated chord or normal-probe blockers
into optical occlusion. All three solved physical segments retain original
scene visibility and endpoint exclusion. Source and receiver must be in vacuum
on opposite sides; the full refracted distance is budgeted. Opaque blockers
produce zero bins; unsupported media, roughness, nonparallel/edge exits,
back-facing receivers and cancellation refuse the complete reading.
The returned spectrum is the selected two-transmission image, with no claim to
total illuminance, reflected/internal-return images, bulk interference,
scattering, calibration bounds or a certified numerical enclosure. Existing
path jobs retain their stricter point-fixture/glass admission. Focused controls
cover analytic axial spreading, independent oblique Snell rays and angular
finite differences, retained s/p power, dispersion, native visibility and
physical-domain/work refusals.

## Incoherent transmitted slab history for calibrated point fixtures

`PhotometricSource::transmitted_slab_history_irradiance` extends that selected
image to the broad-face `T(RR)^mT` family for all nonnegative internal round-trip
orders, subject to an explicit omitted-image tolerance. Each order unfolds
`(2m+1)h` of glass and translates the virtual receiver by `2mh`, preserving the
actual exterior distances. The existing Snell solver and analytic area
Jacobian supply that order's wavelength-dependent source direction and
spreading. Native geometry then traces every folded physical reflection and
transmission. Original s/p power, reversed surface-film order and Beer loss
along the summed internal segments remain intact. No radiance eta factor,
roulette probability or first-image geometric multiplier enters source power.

The actual signed receiver projection applies after the axis-plane solve.
Back-facing/grazing orders, zero source intensity and opaque blockers zero
only that order; later images can still illuminate the receiver. Existing
selected-image back-face/grazing refusals remain unchanged. Additional media,
nonparallel/edge encounters and unresolved geometry still refuse a result.

The source's global angular peak bounds unevaluated orders, including hidden
or unlit earlier images. Bare glass uses monotone s-polarized Fresnel
reflectance to bound both polarization channels by a geometric series.
Passive coated glass uses a decreasing inverse-square spreading series and
integral bound; this can require looser tolerances or more work. Minimum bulk
absorption tightens both. Bound products are evaluated in log space so optical
underflow cannot erase a remainder that a bright source makes representable.
Both slab APIs share source/Jacobian scaling with a log s/p power fallback;
representable irradiance survives when normalized Beer/return throughput
underflows, while ordinary finite arithmetic remains the fast path.
Returned per-bin bounds and their native radiant/CIE integrals describe
truncation of this model, not roundoff, geometric, calibration or measurement
uncertainty. Bulk coherent interference, side-face paths and other scene
illumination remain outside this observable.

Each bin must satisfy its absolute-plus-relative tolerance before the complete
80-bin spectrum is published. Image count, physical distance and total reserved
interface work have independent hard caps; order m reserves `2m+2` interactions
before tracing, including early terminated orders. Cancellation is checked at
each bin, order and physical return. A work cap without convergence is an
explicit error, never a silently truncated complete reading. Geometry is
prepared once and individual histories use fixed local storage.

## Incoherent reflected slab histories for calibrated point fixtures

`PhotometricSource::reflected_slab_history_irradiance` measures the complementary
same-side broad-face family: the front mirror `R`, then `T R^(2k-1) T` for
positive k. The caller supplies a point on the original source-facing broad
face. Body-only discovery verifies that witness, the parallel opposite face,
and the original closed interior; it does not guess a face from the same-side
endpoint chord or treat the witness as the physical reflection point. Source
and receiver remain in vacuum in front of that face. Colocated endpoints are
valid because direct source-to-receiver light is outside this observable.

The front mirror uses its exact virtual-distance area map. Each internal image
mirrors the receiver across the front face and unfolds `2kh` of glass before
using the original dispersive inverse Snell solve and analytic Jacobian. The
actual reflected arrival direction supplies the signed receiver projection.
Native geometry then folds and checks every source, internal and outgoing
segment against the complete original scene. The same shared s/p transport
retains front-versus-inside coatings, `2k-1` physical internal reflections,
and Beer absorption along their actual lengths. Dark, back-facing or blocked
early images cannot remove subsequent illumination. A discovery-only blocker
has no optical effect. Internal normalized-power underflow retains the existing
logarithmic source-scaling fallback.

The front contribution has its own peak-source bound, with factored logarithmic
area scaling. The remaining family reuses the passive spreading bound with
`D_k=H+2kh/n`, `H=a+b`, and the tighter bare-glass bound with Fresnel exponent
`2k-1`; it never applies a negative exponent to the front mirror. Per-bin
absolute-plus-relative tolerances and complete-spectrum cancellation remain
shared with transmission. Reflected order k reserves `2k+1` interface events,
so N evaluated orders consume `N²` events per wavelength. Physical distance,
image and total interface caps refuse incomplete readings. Geometry/BVH
preparation is shared once and individual images allocate no workspace.

`SlabHistorySpectrum` exposes the selected `SlabHistoryFamily`; the transmitted
and reflected result names alias this common reading and its spectral/radiant/
photopic omitted-image bounds. These are analytic model-truncation bounds, not
roundoff, mesh, calibration or measurement enclosures. Side-face paths,
rough/curved/moving surfaces, bulk coherent interference, other illumination,
and general point-fixture specular path-job admission remain outside this API.
Focused controls include axial series and retroreception, independently folded
oblique Snell rays and finite-difference areas, angular source support,
visibility, tilted/grazing reception, a separate complex-amplitude absorbing-
film reference, rescaled extreme attenuation and domain/work refusals.

## Original tracer scenes in calibrated curved transport

`photometry::scene::NativeOpticalScene` admits an ordinary static `tracer::Scene`
to the existing curved polarized transport and its calibrated point connections,
image discovery and finite-detector measurements. Original mesh/chart objects
move into shared immutable geometry without tessellation or a replacement lens;
placed instances retain their object/geometry/frame identities and proper-rigid
transforms. `render_scene()` and `into_scene()` retain the original materials,
emission, lights, environment and camera with the same shared geometry. Raw
meshes/charts require caller-supplied object IDs and immutable geometry hashes.
Animated shapes require an explicit static snapshot.

Every primitive has an explicit material or detector role. Material exteriors
and source initial media use `Ambient` or `InteriorOf(original_primitive_index)`.
Every material dielectric retains a distinct medium slot, including equal-index
bodies with different absorption. Invalid references and self/cyclic containment
refuse; parent validation is linear in the object count. The caller remains
responsible for geometrically correct nesting and outward boundary normals.
Original native medium-transition checks still refuse inconsistent paths.

Dielectric surfaces retain the original dispersion, bulk absorption, roughness,
anisotropy and compact coherent coating. Smooth/rough/brushed conductors retain
their optical tables/coatings and body-local detail. Native spectral diffuse and
coated diffuse laws are unchanged. Lifted Lambertian reflectance is evaluated
from its original coefficients at each physical wavelength throughout the native
380–780 nm band, without fitting or a new table allocation. Legacy Schlick GGX
currently refuses rather than becoming a surrogate absorber. Original gray and
spectral diffuse thermal materials retain their reciprocal reflection and exact
physical source, as specified in the thermal camera extension below.
Detectors explicitly intercept the chosen original geometry. All other original
geometry remains present, including physical occluders.

Object allocations are bounded and checked; `with_curved_scene` builds the native
boundary accelerator once per observation batch under its existing work budget.
Per-path and per-solve limits, cancellation, actual spectral Snell/Jacobian
transport, Stokes transfer, and native visibility remain unchanged. Selected
connections retain selected-branch semantics; image seed discovery does not
prove completeness. Curved source observations measure the explicitly supplied
source contribution, not all emitters retained in the ordinary rendering view.
Regressions connect the original dispersive ball chart to independent axial
ray-transfer/Fresnel/Beer/CIE values, preserve embedded nested absorption even
at equal phase indices, retain black blockers and shared instance identity, and
refuse malformed topology, unsupported laws and exhausted build budgets.

`Interaction::DepolarizingLiftedDiffuse` also carries that exact law through
ordinary polarized paths, curved/participating walks, selected diffuse bundle
endpoints and reciprocal camera receivers. It shares the native depolarizing
cosine sampler and the same calibrated point/area-light response and proposal
PDFs as tabulated diffuse. Finite coefficients are checked at admission;
out-of-band wavelengths or unrepresentable reflectance refuse at evaluation.
Exactly black response consumes no direction sample. Thermal Kirchhoff bindings
still require their existing exact supported thermal curve. Regressions compare
native straight/curved power and depolarization to the original polynomial, and
connected fixture/camera radiance to independent cosine/Beer/inverse-square
values at multiple wavelengths while retaining the same continuation path.

## Original lifted emission in native reciprocal camera transport

`RadianceSpectrum::try_lifted` retains the original `LiftedSpectrum` coefficients
and an explicit local radiance scale in W/(m² sr nm). Its support is exactly
380–780 nm; no knots are fitted or allocated. `lifted_source()` reports that
original pair and `knots()` is empty for the analytic representation. Existing
tabulated spectra retain their original knots, support and interpolation
arithmetic. Nonfinite coefficients or invalid calibration refuse at admission;
out-of-band queries and unrepresentable evaluated radiance refuse at observation.
Zero scale is exact black within its declared source band.

The existing `SurfaceEmitter` owns this representation, so camera-hit emission,
native mesh/chart area proposals, complementary light-sampling MIS and molecular
camera consumers query the same source at the actual physical wavelength. The
original emitting side, local-medium calibration, emission-before-scattering,
reciprocal index conversion and Beer attenuation remain unchanged. A lifted
RGB spectrum remains authored source data; calibration does not recover a
measured spectrum or impose a thermal emission law. Regressions retain exact
interior-wavelength evaluation and compare native camera radiance and original
area-source irradiance against independent area/distance²/Beer values, including
exact-black source behavior and explicit source-domain refusals.

## Original scene cameras and emitters through native refraction

`NativeOpticalScene::with_curved_camera` connects an admitted original tracer
scene to the existing reciprocal camera and its raster/area-light consumers.
The caller supplies independent surface and environment conversions from
authored radiance units to W/(m² sr nm) using `NativeRadianceCalibration`.
Geometry already uses metres. Each original emissive primitive registers once,
with its exact lifted spectrum, authored strength and original two-sided law.
Material reflectance does not multiply local emission: native camera transport
adds emission on arrival before sampling the original scattering material.
The original environment is borrowed with its existing RGB spectral
reconstruction and declared visible support; calibration does not recover a
measured spectrum. Zero calibration is exact black for these authored sources
in the declared band; original thermal materials already have SI source units.

The original camera eye, orthonormal forward/up basis and half-tangent projection
are retained by the existing exact legacy-compatible pinhole branch. Callers
select the launch medium with the optical owner's original-primitive references;
camera containment is not inferred. Actual refracted paths, reciprocal
`n_camera²/n_source²`, Beer attenuation, native proposal probabilities and
unresolved-continuation errors remain owned by the same curved transport.
`CurvedCameraRaster` traces each of its four wavelengths independently, with the
existing deterministic sample identities and normalized XYZ integration.
Detector-role surfaces remain unsupported by the reciprocal camera owner.

Legacy rectangles remain NEE metadata naming original physical sources. They
undergo existing rectangle/duplicate admission and the same primitive emission
and geometry predicate as ordinary render admission. They add no second source
or geometry. `CurvedAreaEmitter::from_camera` can prepare each shared original
surface and exact spectrum for native area-light MIS under explicit face,
geometry and image budgets. Legacy proposals do not infer refracted image seeds,
an automatic NEE strategy or image completeness.

Camera/source preparation checkpoints cancellation and bounds source storage by
the admitted object count; native acceleration uses the caller's build-work
budget. Invalid camera bases, source coefficients/strengths, unrepresentable
calibration products, mismatched light ownership and exhausted budgets refuse
before an observation callback. Observation errors remain distinct from
admission errors. Regressions cover original emission through a glass sphere,
immersed-camera index scaling and actual Beer lengths, two-sided emission before
diffuse/environment continuation, raster spectral/XYZ replay, exact source and
geometry identity for area proposals, and absence of duplicate rectangle power.

## Original thermal surfaces in refractive camera observations

`NativeOpticalScene` now retains original `GrayDiffuseThermal` and
`SpectralDiffuseThermal` materials instead of refusing them. Native diffuse
reflection is exactly the reciprocal `1-epsilon` law. Original gray sources
retain their accepted temperature, declared gray model and complete source
identity; the constant native reflectance is prepared over their original
surrounding medium's phase-law band. Spectral sources retain every original
emissivity knot and their finite source band. Neither preparation fits Planck
radiance, changes emissivity or extends a material/medium domain.

`DiffuseThermalEmitter::try_gray_in_band` also makes this original gray source
available directly to native camera and area-source users. `surface()` returns
the existing `DiffuseThermalSurface` view, preserving the original gray or
spectral value and identity. Exact source/material agreement remains required;
only a blackbody may use native opaque termination. Source storage is bounded
by the original object count and prepared once; path and spectral-lane queries
reuse the original source law without allocations or new geometry.

`with_curved_camera` attaches those accepted physical sources to the same native
camera used by raster, sensor, area-light and molecular observations. Thermal
radiance is already W/(m² sr nm), so neither authored surface calibration nor
environment calibration scales it. The original two-sided law, emission before
scattering, refracted path, Beer attenuation and camera/source index ratio are
unchanged. Native area sampling and complementary emitter-hit MIS retain the
same thermal snapshot and original geometry. A thermal material with any
`Primitive::emission` metadata refuses at scene admission, including a zero
strength, rather than creating two emission owners.

The original vacuum/index-one emission restriction remains local to the actual
source contact. A camera inside glass can observe a thermal source outside that
glass; a source immersed in an unsupported host still refuses. This observes a
frozen accepted thermal state and does not advance heat, infer bolometric power
from finite spectral data, invent temperature for other materials, or supply
coherent thermal fields. Four focused G0/G1/G3 regression groups cover analytic
Planck radiance through glass with actual absorption lengths and immersed-camera
scaling, gray/spectral equilibrium reflection and complementary area MIS,
infrared photon-weighted detector response with replay and temperature changes,
and source/domain/material/work refusals. Native Rust execution remains
unavailable in the implementation environment; scoped formatting and source
review are separate from those authored regressions.

## Thermal gas fields driving native refraction

`material_state::gas_scattering::isobaric::IsobaricGasField` borrows an original
single-component molecular gas model and an immutable temperature reconstruction.
It evaluates the original gas EOS at each actual temperature and fixed absolute
pressure. Exact density gradients and Hessians include both the thermal Hessian
and the gradient outer product; mixed spatial derivatives are retained. The
existing molecular Lorentz--Lorenz owner then supplies wavelength-dependent
phase/group indices and their nonlinear density derivatives to native rays and
angular bundles. Preparation binds the original source coefficients, explicit
reference state and accepted snapshot. Admitted stage queries allocate nothing.

The explicit reference temperature prepares optical coefficients; it never
replaces a spatial temperature. The original source T/p domain, water label and
thermal-provider cancellation/refusals remain authoritative. The adapter admits
one component or declared pseudocomponent, fixed composition and fixed pressure;
it does not supply flow, pressure/composition gradients, gas-phase existence,
condensation, thermal-contact traversal or true absorption. Focused G1/G3 tests
connect an original generated conduction solution to rotated density jets,
independent Lorentz--Lorenz derivatives, analytic divergent-lens footprints,
curved-ray Snell invariants, phase/group paths and source-domain refusals.

## Original scene media driven by accepted thermal material fields

With `photometry` and `material-state`, `NativeThermoOpticBinding` connects an
original `Ambient` or `InteriorOf` medium to its accepted
`MaterialThermoOpticField`. `with_thermo_optic_scene` and
`with_thermo_optic_camera` use that actual spatial field throughout native
propagation and local interface scattering. The original reference glass must
match in all dispersion, absorption and provenance values. Distinct physical
medium slots remain distinct even when their reference indices agree; each
slot accepts at most one field, and unlisted slots retain their original law.

The borrowed field supplies temperature-derived phase and group indices,
analytic spatial gradient, reference absorption and, when explicitly admitted,
the analytic Hessian. First-order bindings support rays and cameras but refuse
differentiated bundle connections. The differentiated binding carries the
original chain rule through curved boundary arrival and Snell/Fresnel maps,
so native point/area connections and detector footprints see actual thermal
focusing. It does not infer a Hessian from a locally zero gradient. The native
field identity binds the complete reference optics, thermo-optic constitutive
law and frozen temperature/reconstruction snapshot; derivative capability does
not alter that physical identity.

Positions are world metres, temperatures kelvin, absorption inverse metres and
wavelength vacuum nanometres. The accepted snapshot stays fixed throughout an
observation batch. Field source domains, optical step resolution and native
boundary contact collars remain authoritative; an unsupported spatial point,
temperature or wavelength refuses the observation without nominal-glass
fallback. Internal discontinuities still require their original separate
region/event treatment. Thermal deformation, missing constitutive response,
heat evolution and thermal source emission are not inferred from this binding.

Both entry points reuse the original immutable geometry, topology and optical
materials. Camera preparation is shared with `with_curved_camera`, including
the original physical camera, calibrated authored emission, environment and
exact SI thermal sources. The ordinary `render_scene()` retains the original
reference-glass rendering view. Medium preparation is linear in the original
medium/binding counts, allocates one bounded medium table, and borrows fields
without copies; native geometry is prepared once per callback. Existing query
budgets and cancellation remain in force.

Focused regressions compare spectral thermal bending inside the original ball
and refraction onto an off-axis original emitter against a conserved-momentum
analytic path, including phase/group lengths, Beer depth and immersed-camera
radiance. A quadratic-temperature field compares differentiated focusing after
the actual spherical interface against the analytic paraxial area map. Further
controls retain the homogeneous reference limit, physical snapshot identity,
explicit Hessian capability and source/domain/work refusals. These native Rust
regressions are authored; this environment permits scoped parser/format checks
and source/math review but does not provide native Cargo execution.

## Original scenes with source-backed molecular gas ambient

With `photometry` and `material-state`, `NativeGasOpticalScene` admits original
scene geometry and materials with an actual molecular ambient. Its
`NativeGasAmbient` can be constructed only from existing
`MaterialGasRefractiveField`, `MaterialGasCompositionRefractiveField` or
`MaterialAirRefractiveField` owners, including an `IsobaricGasField` refraction
binding or explicit Ciddor pressure/temperature reconstruction. Material identity and source
band come from those owners; no nominal glass, vacuum placeholder, fitted index
or caller-asserted material identity is introduced. The distinct scene type
does not expose ordinary `render_scene`/`into_scene` views that would lose the
gas physics. Existing glass-scene APIs retain their behavior.

Every gas scene/camera observation requires an explicit accepted ambient field.
Its material identity and source band must match scene admission, even for
numerically equal optics or local vacuum. Different accepted frozen snapshots
of that same material may serve different batches. The native field identity
remains the exact original source/snapshot identity. Interior homogeneous glass
and optional reference-matched thermo-optic fields retain their original slots,
surfaces and absorption; another ambient binding is duplicate ownership and
refuses. Gas ambient has no fabricated reference glass.

Fixed-composition density fields retain the original visible, nondispersive
Lorentz--Lorenz approximation. Composition fields use their explicit source
bands, including supported infrared dispersion, with the full mixture phase
and group derivatives. First-order bindings refuse angular-bundle queries;
differentiable bindings retain the actual density/component Hessians and both
nonlinear chain-rule terms. A typed prepared-composition binding verifies its
exact source snapshot and reuses the existing stack-owned molecular wavelength
coefficients across stages, boundary segments and connection trials. Its one
selected wavelength remains mandatory; source-band metadata is not a license
to reuse prepared coefficients at another wavelength. No optical heap workspace,
cache, finite-difference Hessian or independent gas EOS is added.

Camera/source preparation is shared with the original native camera. Authored
radiance calibration, exact thermal sources, reflection, source spectral support
and emitting-host restrictions remain intact. These refraction bindings
explicitly select zero true absorption. For molecular-composition bindings,
the same composition scattering field attaches
directly to the existing `CurvedRayleighScene` and `CurvedRayleighCameraScene`
using its unchanged refractive identity, so elastic collision and survival
sampling are not double-counted as Beer absorption. Refraction alone does not
claim the complete scattered radiance. Fields use world metres and their
declared SI density bases/derivatives. Snapshot, dilute/EOS, spatial and spectral
domains, contact collars, cancellation and native work limits remain binding.

Four focused regression groups cover an isobaric thermal gradient bending the
original camera onto an off-axis original source against an independent exact
conserved-momentum solution; visible/infrared composition focusing combined with
thermal glass through the actual original ball interfaces and detector;
prepared/unprepared agreement and source-coupled elastic molecular camera
collisions; and material, snapshot, capability, domain, duplicate, work and
cancellation refusals. Native Rust execution is not available in the current
implementation environment; parser/format checks and independent source/math
review are recorded separately from the authored regressions.

## Atmospheric air through original state fields and native cameras

`material_state::air_refraction::MaterialAirRefractiveField` borrows the existing
joint pressure/temperature reconstruction and the original Ciddor composition.
Its immutable identity binds the source snapshot and the complete material law.
Analytic optical gradients and Hessians retain both original state Hessians and
the full TT, PP and mixed TP chain rule. Constant temperature does not suppress
pressure-driven bending. The prepared form caches wavelength coefficients outside
ray/bundle stages and refuses any different wavelength.

`try_nonabsorbing` explicitly selects a zero-loss approximation; refractive-index
data do not establish atmospheric transparency or absence of scattering.
`NativeGasAmbient` accepts prepared or unprepared first- and second-order Ciddor
fields on the existing original-scene camera path. Phase, group path, focusing,
source identity and refractive radiance scaling use the existing transport
owners. Variable-composition spatial gradients, phase equilibrium, droplets,
moving-fluid effects and attenuation/scattering require separate source owners.
Focused exact-source controls cover joint-state derivatives, analytic
pressure-lens focusing and source-domain refusals. The original off-axis scene
regression requires the complete native camera feature build.

## Absolute atmospheric and air-relative glass indices

With `material-state`, `CiddorAirIor` retains the original atmospheric composition,
absolute pressure and temperature as a native phase/group law. The independent
`DielectricGlass` absorption source remains authoritative.

`AirRelativeSellmeierIor` retains original relative-index oscillator coefficients
and source support, together with an explicit index-reference atmosphere and an
independently declared vacuum or air-wavelength coordinate. For reference-air
index A, coordinate-air index B and relative source R, absolute phase is
A R(lambda/B). Group index includes both air dispersion and the coordinate
chain rule; equal atmospheres reduce to the product of their group indices.
Changing the actual scene atmosphere leaves the measurement reference intact.

Source-coordinate endpoints are mapped inward into vacuum support once at
construction, intersecting the original Ciddor band. Ray queries do no root
solving, fitting, interpolation or allocation. Native Fresnel/lens, bulk
phase/group transport, coating and material identities consume the same absolute
law; existing identity tags are unchanged. No automatic catalog convention,
source authentication, glass temperature evolution or source extrapolation is
inferred. Focused exact-source controls exercise all coordinate conventions,
independent group differences, source support, Fresnel/lens power and identities.

## Joint pressure and temperature driving molecular gas refraction

`material_state::gas_scattering::pressure_temperature::PressureTemperatureGasField`
borrows the original single-component gas model and the existing joint absolute
pressure/temperature reconstruction. Each stage uses one original field query
and allocation-free native gas-state evaluation. Exact density gradients and
Hessians include both pressure and temperature curvature and the mixed p/T
response; the existing molecular Lorentz–Lorenz owner supplies optical jets and
native ray/bundle propagation. An original P1 joint element is accepted directly,
without a new reconstruction, copied source buffers or a renderer-side EOS.
Original gas validity, water labeling, provider cancellation, source support and
dilute-optics refusal remain in force. This frozen single-component or declared
pseudocomponent model does not establish phase equilibrium, real-gas corrections,
moving-medium optics, new dispersion data or scattered radiance. It does not
smooth gradient jumps or cross native mesh faces by itself.
G1/G3 tests cover independent native P1 scalar-source derivatives, pressure-lens
ray/bundle focusing, correlated p/T constant-density cancellation and original
source/domain/cancellation refusals. Native execution is reported separately.

## Native mesh traversal of molecular gas refraction

Under `material-state,volumes`, `thermo_optic::transient::mesh::gas` extends the
original `P1ConductionMesh` with isobaric and joint pressure/temperature gas ray
and angular-bundle propagation. Each stage retains the original scalar gas EOS,
molecular optical model and full density/index Hessian. The existing mesh owns
face admission, the complete departure collar and crossing-time corrections;
continuous shared nodal state introduces no artificial Fresnel jump. Pressure-only
gradient jumps participate in the native one-sided bundle endpoint checks.
The joint path reuses the existing two-element pressure cache and one barycentric
query for both fields. It accepts only the original mesh and temperature buffer,
and requires the original joint p/T snapshot identity. Contacts with duplicated
nodes and exterior faces terminate at the native boundary; source, numerical,
face-budget and cancellation refusals remain explicit. These frozen gas fields
add no pressure solve, flow coupling, phase-equilibrium claim or scattered radiance.
Four G1/G3 regressions compare an unsplit analytic thermal field, pressure-jump
bundle variations, correlated p/T constant-index transport and source/contact/work
refusals. Native execution is reported separately from source and format checks.

## Original water scene media driven by accepted thermodynamic state

With `photometry` and `material-state`, `NativeMaterialFieldBinding` connects an
original IAPWS water medium to its existing `MaterialWaterRefractiveField`.
`NativeThermoOpticBinding` remains a direct alias with its existing thermal-glass
constructors. The same original scene, gas-ambient scene and camera entry points
now accept either material owner; medium topology, original geometry, authored
radiance calibration and exact thermal sources retain their existing handling.

Water admission requires an actual original `DielectricIor::Water` reference
with the released IAPWS law and exactly the field's independently owned
absorption. The complete reference, including its state and provenance, must
match the selected original medium. A fitted equal-index glass, changed
absorption or second binding to the same slot refuses. Native identity includes
that complete original reference, optical law and accepted source snapshot.
Prepared fields expose their actual source owner for this admission; equal
snapshot labels cannot conceal different absorption. Prepared and unprepared
bindings have the same physical identity, while prepared queries retain their
one exact vacuum wavelength and reuse the original wavelength coefficients.

The original field supplies the actual local temperature and density, or the
original IF97 isobaric/joint pressure-temperature reconstruction and its selected
single-phase branch. Reference state values are never substituted or added as
an index correction. The existing optical and EOS chain rules supply phase and
group indices, spatial gradients and, for explicitly differentiated bindings,
complete state Hessians. Pressure and temperature curvature both contribute to
native bundle focusing and boundary Snell/Fresnel maps. First-order bindings
support rays and cameras and explicitly refuse differentiated bundle queries;
a locally zero gradient does not imply a zero Hessian.

Positions remain world metres, pressure absolute pascals, temperature kelvin,
density kilograms per cubic metre, absorption inverse metres and vacuum
wavelength nanometres. The source stays frozen throughout the observation.
Original optical, EOS, spatial, spectral and contact-collar domains remain
authoritative. No absorption, phase crossing, fluid motion, thermal deformation
or source emission is inferred from water refraction. Unsupported local states
refuse instead of falling back to the reference. Existing bounded medium
preparation and native ray, event, bundle and camera solvers are reused without
an additional propagation solver or field copies.

Four focused regression groups exercise joint pressure/temperature focusing
through the original spherical water interface against its analytic paraxial
area map; an isobaric temperature gradient bending the original camera onto an
off-axis original source against independent conserved-momentum quadrature,
including phase/group paths, Beer depth and immersed radiance; water interiors
and an actual molecular ambient in the same camera; and original reference,
absorption, prepared wavelength, capability, domain and work refusals. Native
Rust execution is unavailable in this environment; these authored regressions
are distinguished from scoped parser/format checks and source/math review.

## Original gas composition on native tetrahedral meshes

Under `material-state,volumes`,
`thermo_optic::transient::mesh::gas::composition` borrows original nodal
partial densities in the molecular model's exact order and declared mass or
molecule-number basis. One native barycentric query supplies every component
and its unsmoothed P1 gradient. Original molecular polarizabilities, dispersion,
phase/group indices and nonlinear optical Hessians drive native rays and
angular bundles, including composition-driven bending at constant total density.

Preparation validates the borrowed source once. Each propagation prepares one
wavelength and retains a private two-element polynomial cache; stage queries
allocate nothing. Original mesh events, complete collars, crossing-time
corrections, spectral domains, work budgets and cancellation remain authoritative.
Continuous internal faces create no Fresnel event; duplicated contact topology
and exterior faces stop. Zero species and vacuum stay explicit. No fraction
normalization, new EOS, cell-average reconstruction, moving-medium model,
absorption or scattered-radiance integration is inferred.

Five G1/G3 regression groups cover both density bases, original P1 derivatives
and optical Hessians, chromatic bending and group paths, concentration-jump
bundle variations, constant-index/vacuum limits, and source/topology/work
refusals. Native execution is reported separately from source/format checks.

## Finite-source histories on refractive camera paths

`polarization::dielectric::scene::curved::camera::transient` observes
prescribed finite `SurfaceEmitter` radiance histories on the original camera
paths. Each encountered source is evaluated at the observation time minus
that encounter's accumulated group optical length divided by the exact vacuum
light speed. Earlier source contributions never inherit the terminal path's
delay. Phase index still owns bending, Fresnel response and radiance conversion.

The observation clock is double-double seconds at the launched ray/pupil.
Histories retain that clock until subtracting their own physical epoch. Missing
group information, unresolved delays, invalid histories and cancellation refuse
the complete sample. Histories multiply the original radiance spectrum, consume
no random variates, and leave the original path, Beer attenuation and reciprocal
polarization response intact. Exact unit histories preserve static arithmetic.

This is frozen-geometry, frozen-medium, incoherent spectral envelope transport
with a narrow-band group-delay approximation. It supplies no carrier phase,
within-bin pulse reshaping, moving-medium frequency shift, fluorescence,
unmodeled lens transit, interface/coating residence time, or temporal environment,
volume history. The finite-area and point-fixture adapters below support timed
clear-medium and molecular connections, including native refracted source
images. Unbound sources remain explicitly steady.

Seven G1/G3/G4 tests cover native static/RNG parity, distinct per-source delays,
chromatic pulse arrivals, an independent curved-path group integral, large-epoch
invariance, source/domain refusals and late cancellation. Numerical execution
remains separate from source and formatting checks.

## Transient spectral pixels and quantum-efficiency channels

`camera::raster::transient` connects the finite-source history observer to the
original raster, pupil, CIE/XYZ and quantum-efficiency sensor integrations.
`sample_pixel_transient_at_time` and `sample_pixel_sensor_transient_at_time`
trace each original wavelength independently and evaluate its own refracted
emission times. All detector channels share the same four spectral paths;
adding QE channels adds no propagation or source-history calls.

Original random-counter domains, wavelength proposals, spectral support,
photon-energy conversion and aperture-average conventions are retained.
A failed source or wavelength lane refuses the complete pixel, preserving its
original error without retrying or renormalizing surviving channels. Results
are instantaneous radiance and photoelectron rates at the declared DD launch
clock, before exposure, etendue, lens transit and electronics. The underlying
finite-source and frozen-medium assumptions remain unchanged.

Three G1/G3/G4 tests cover static pixel/sensor parity and replay, channel-count
independent transport work, chromatic pulse energy/photon references, and
failed-lane/cancellation behavior. These are authored numerical regressions;
execution is reported separately.

## Passive analyzer measurements through native refracting cameras

`CameraAnalyzer` supplies a declared passive detector intensity row and a world
reference axis to the existing reciprocal curved camera. Authored rows are
admitted only when `0 < a0 <= 1` and `norm([aQ,aU,aV]) <= min(a0,1-a0)`;
their response is consequently between zero and the incident intensity for
every physical Stokes vector. Ideal linear and circular constructors retain
their physical one-half unpolarized transmission. All-zero rows refuse under
the existing extinction/underflow semantics. A polarized analyzer whose axis
cannot define the queried transverse frame also refuses.

For the original opposite-travel frame convention, the reciprocal seed is
`[a0,aQ,-aU,aV]`. The full Stokes vector traverses the original Fresnel, TIR,
coating and gradient-index owners. The analyzer is therefore applied inside
transport, with the original branch proposal division and endpoint `n^2`
conversion exactly once. It is not an additional scalar postprocessing loss.
An identity row takes the original scalar launch path exactly. Nonpolarizing
rows do not require a distinguished transverse axis.

`CurvedCameraScene::sample_ray_analyzed` and
`sample_ray_analyzed_with_samples` observe the original surface emissions,
including optional `SourcePolarization`, and unpolarized environment and
supported thermal emissions. The geometric spectral raster
adds `sample_pixel_analyzed` and `sample_pixel_sensor_analyzed`; the original
four wavelength-dependent paths, counter streams, CIE or multichannel QE
responses, lens/pupil weights and all-or-error pixel behavior are retained.
No scene-table clone, extra traversal, per-query allocation or second spectral
integration is introduced. Native phase/group paths, Beer absorption, medium
identity, cancellation and work limits remain authoritative.

The original `CameraAnalyzer` is an ideal wavelength-independent detector row;
the spectral extension below retains a source-owned chromatic row. Neither
provides a finite-thickness or angle-dependent polarizer, bulk
birefringence or analyzer overloads for explicit emitter-connection consumers.
The molecular-scattering extension is described below.

Focused regressions specify passive-row admission and independent Malus
values, a native Brewster/rotated-TIR path resolving both U and V, a curved
gradient-index path that reaches an emitter missed by its straight launch,
and the original XYZ and multichannel QE raster. The numerical expectations
include native proposal probabilities, Beer attenuation, endpoint indices,
identity replay and failure semantics. The implementation passed pinned Rust
formatter/parser checks; these native regressions were added but were not
executed in the shared editing environment.

## Passive analyzers through native molecular-scattering cameras

With `volumes`, `CurvedRayleighCameraScene` accepts the same original passive
`CameraAnalyzer` through `sample_ray_analyzed`,
`sample_ray_analyzed_with_samples` and
`sample_ray_analyzed_with_tracking_samples`. The reciprocal seed enters the
original Rayleigh/depolarized-Rayleigh Mueller walk once. Frame rotations,
curved index propagation, discrete refractive interfaces and molecular
collisions retain the complete Stokes state. The original scalar scattering
proposal divides the Mueller response once; scattering survival remains in
the original free-flight law. No second analyzer transmission or scalar
scattering extinction is applied after transport.

Uniform and heterogeneous density use their original flight and independent
acceptance variates. Rejected tracking candidates preserve polarization.
Observed surface and supported thermal emission, original phase/group paths,
medium identities, local refractive indices, Beer attenuation and endpoint
`n^2` conversion remain in their existing owners. Identity analysis follows
the exact original launch. The admitted molecular law samples direction from
the original incoming direction and world frame, independently of the
analyzer's polarization reference axis. Original counter streams, work,
volume-event and distance limits, cancellation and unresolved-tail failures
are retained.

The spectral raster adds `sample_pixel_rayleigh_analyzed` and
`sample_pixel_sensor_rayleigh_analyzed`. Their original four wavelength paths
feed the original XYZ or multichannel QE/photoconversion integration; each
lane retains its volume-event and free-flight counts. Additional detector
channels do not add propagation, and a failed lane publishes no partial pixel.
The implementation shares the original launched-beam traversal bodies and
result construction without cloning scene tables or adding path traversals.

The physical scope remains the original admitted reciprocal molecular
Mueller laws and the original source models, including authored surface
polarization. The detector is the declared ideal local analyzer, without
resolved polarizer/lens propagation, bulk birefringence, or analyzer overloads for explicit emitter
connections or the separate molecular transient-history consumers.

Three focused regression groups specify right-angle Rayleigh/depolarized
responses and rotated analyzer frames; a closed-form gradient-index path with
an actual collision, discrete refraction, Beer/index/Fresnel factors and a
rejected tracking candidate; and original XYZ/QE pixels with identity replay
and complementary analyzer channels reconstructing the scalar signal. The
fixed raster fixture's native Philox first-flight variates were independently
checked to include an interior collision. Formatter/parser checks passed;
these native Rust regressions were added but were not executed in the shared
editing environment.


## Retarded finite emitters on refracted molecular-scattering paths

With `volumes`, `camera::transient::participating::TransientCurvedRayleighCameraScene`
composes the original admitted surface histories with the native curved Rayleigh
and depolarized-Rayleigh walk. Both borrow the same original scene. Each actual
emitting encounter queries its source at the observation clock minus that
encounter's accumulated group path divided by c, including all preceding curved
legs and real molecular turns. Rejected tracking candidates preserve that clock.

The original free-flight, delta-acceptance and surface streams, full Mueller
weights, Beer absorption, endpoint index conversion, thermal emission, work
limits and cancellation remain authoritative. History queries consume no random
variates and introduce no additional collision, extinction or albedo factor.
Seeded, uniform-density custom and heterogeneous-tracking custom paths retain
their original sampling contracts. Source failure or an unresolved native path
refuses the complete sample; finite-source histories remain separate from steady
volume and environment emission. Geometry and media remain frozen, molecular
events instantaneous, and spectral envelopes incoherent under the original
group-delay approximation.

Four G1/G3/G4 regression groups specify two molecular turns with an independent
polarization/absorption/chromatic-delay reference, the analytic curved group
integral through a rejected tracking candidate, static and zero-scattering
replay, scene identity and history/domain/work/cancellation refusals. Pinned
Rust formatting/parser checks passed. The native Rust regressions were authored
but not executed in the shared editing environment.


## Spectral pixels and QE rates with retarded molecular paths

With `volumes`, the native raster adds `sample_pixel_rayleigh_transient_at_time`
and `sample_pixel_sensor_rayleigh_transient_at_time`. Each wavelength reaches
finite source histories along its actual refracted and scattered path; all
sensor channels share those same four original trajectories. Original CIE and
QE proposals, photon-energy conversion, pupil/diffraction weights, surface
Philox domain `0x6372_7340 + lane`, and tracking domain `0x6372_7350 + lane`
remain unchanged. Additional QE channels add neither traces nor history calls.

Samples retain their normalized DD observation clock and per-lane native
volume-event and free-flight counts. Histories, unresolved transport or cancelled
lanes refuse the complete pixel while preserving the original source error;
remaining lanes are not renormalized. Outputs are instantaneous aperture-average
radiance and photoelectron rates under the original frozen-medium, finite-surface
history scope, before exposure, etendue and sensor electronics.

Three regression groups use a closed emitting geometry with actual seeded
molecular collisions: static CIE/QE parity and replay, channel-independent work,
per-wavelength scattered group delay against independent energy/photon sums,
and late-history/work/cancellation refusal. Pinned Rust formatting/parser checks
passed; these native tests were authored but not executed in the shared editing
environment.


## Finite shutters over original retarded surface-source histories

`camera::transient::exposure` adds an explicit positive rectangular
`ExposureWindow` on the original DD launch/pupil clock. Each encountered emitter
receives both shutter endpoints shifted by its own native accumulated group delay.
For frozen geometry, media and surface laws, linearity permits one original trace
per spectral sample. `SurfaceEmissionHistory::average` supplies the source-owned
interval mean; its default returns `ExposureUnsupported` rather than guessing an
integral from instantaneous samples. Existing instantaneous observations retain
the same shared source identity, group-delay admission and clock checks.

Seeded and caller-sampled camera observations preserve original path and random
work, yielding mean spectral radiance and its time integral in J/(m² sr nm).
Unbound surface, volume/LTE and environment emission remains steady. Every
source must complete a nonnegative finite interval mean; late source failure or
cancellation refuses the entire observation. This owner does not infer motion,
a shutter response curve, lens transit, fluorescence or sensor electronics.

`PiecewiseLinearHistory` borrows the original strictly ordered nonnegative
source knots and retains their DD epoch. It integrates every crossed segment
of the declared interpolant, including clipped endpoints and narrow interior
pulses, without resampling, query allocation, extrapolation or subtracting large
cumulative integrals. Point lookup is logarithmic; interval work scales with
crossed knots and polls cancellation. Normalized DD contributions and a scaled
operation order preserve representable pulse area under extreme duration and
intensity ranges. Exact interpolant integration does not assert measured-source
bandwidth or physical calibration.

Four G1/G3/G4 regression groups specify narrow/clipped pulse areas at ordinary
and large epochs, extreme representable pulse means, original seeded/caller
path replay, different delayed windows at two emitting encounters, invalid
shutters/knots/domains, unsupported interval integration and late failures.
Pinned Rust formatting/parser and source review passed. Native Rust execution
was not performed in the shared editing environment.


## Original CIE and QE detector observations over finite shutters

`CurvedCameraRaster::sample_pixel_during_exposure` observes mean source radiance
through the original four wavelength paths and CIE integration.
`sample_pixel_sensor_during_exposure` shares those paths and source means across
all QE channels, then applies the explicit DD duration to the native spectral
energy and electron rates. Results expose mean detector samples, QE-weighted
energy in J/(m² sr), and expected photoelectrons per m² sr. Original pupil,
diffraction, wavelength proposal, QE and photon-energy weights remain unchanged.

This is an aperture-average observation before pixel etendue, gain, saturation,
read noise or lens transit. No additional temporal samples, random draws, source
interpolation or path integrations are introduced. Each wavelength uses every
emitter's own delayed shutter interval. A failed source integral or wavelength
lane refuses the complete observation without retrying or renormalizing the
remaining detector channels.

Three regression groups specify exact constant-source CIE/QE parity, four shared
paths independent of channel count, explicit integrated units, independent
chromatic triangular-pulse area and Beer/QE/photon references at ordinary and
large DD epochs, and transactional failure of a later spectral lane. Pinned Rust
formatting/parser checks passed; the authored native tests were not executed in
the shared editing environment.


## Finite exposure through the original refracted molecular walk

The transient molecular camera and raster now compose the same explicit shutter
and source-owned interval means with their original participating-medium paths.
Seeded, uniform custom-flight and heterogeneous tracking exposure entry points
trace each spectral path once. Every actual surface source uses the accumulated
group delay after all preceding refractive legs and molecular turns. Native
volume/LTE emission, Beer attenuation, Mueller weights, collision proposals,
null-candidate acceptance, path limits and cancellation remain unchanged.

`sample_pixel_rayleigh_during_exposure` retains the original CIE mean observation;
`sample_pixel_sensor_rayleigh_during_exposure` shares the original four paths
across all QE channels and uses the common exposure conversion to J/(m² sr)
and photoelectrons per m² sr. Both retain per-lane volume/free-flight counters,
original random domains and complete-or-error output semantics. There are no
new temporal samples, repeated paths or unsupported point-sampling fallbacks.
The original frozen-medium, instantaneous molecular-event and finite-surface
history boundaries continue to apply.

Two focused regression groups specify a narrow triangular pulse after two
molecular turns at a large DD epoch, one source mean and unchanged native path/
flight work; and actual seeded scattering with steady CIE parity, independent
chromatic integrated energy/photon references, one-versus-three-channel shared
work and unsupported-history refusal. Pinned Rust formatting/parser checks
passed for both module trees; native test execution was not performed in the
shared editing environment.

## Polarization-resolved transient camera and shutter observations

`TransientCurvedCameraScene::with_analyzer` configures the existing
`CameraAnalyzer` at the original launch plane. The optional setting feeds the
same complete-Stokes reciprocal seed into every clear and Rayleigh
instantaneous/exposure launch, including seeded, explicit-flight and tracking
observations. The existing transient CIE XYZ and QE raster adapters inherit
that prepared configuration; no additional tracing or detector wrapper is
introduced. Default construction retains the original analyzer-free behavior.

Every original source keeps its own DD retarded group clock and finite-shutter
interval mean. Native refractive/Fresnel and molecular Mueller responses,
proposal probabilities, Beer loss, volume emission, endpoint radiance
conversion and random domains remain with the existing transport owners.
Configuration moves the original history bindings and adds no samples,
path integrations, random draws or source-history allocations. An invalid
projected analyzer axis refuses the complete observation.

The constant analyzer remains an ideal passive wavelength-independent launch-plane response;
`with_spectral_analyzer` instead borrows the chromatic detector law described below,
resolved at each original wavelength before launch. Both operate
over the original emitters, frozen geometry/media and narrowband
envelope/group-delay assumptions. A prescribed intensity history scales the
source's complete spectral Stokes state together. It introduces no lens/analyzer
transit, time-varying polarization fractions or coherent optical carrier evolution.

Three focused regression groups specify native oblique-Fresnel analyzer
parity and original history identity; complementary complete-Stokes molecular
observations of a picosecond pulse inside a finite shutter at a large DD
epoch; and chromatic XYZ/QE observation, shared four-path work, zero-response
channels and unsupported-source-integral refusal. Pinned Rust formatting/
parser checks and bounded source/API/physics review passed. The native Rust
regressions were authored but not executed in the shared editing environment.

## Source-owned spectral polarization detectors

`SpectralCameraAnalyzer` retains 2..=1024 ordered positive vacuum-wavelength
knots and a common world reference axis. Every knot is an absolute passive
Mueller intensity row: `0 < a0 <= 1` and `norm([aQ,aU,aV]) <= min(a0,1-a0)`.
All-zero rows retain the existing camera's extinction/underflow refusal.
The source declaration is immutable; no measured accuracy or missing spectral
support is inferred. Queries return exact authored knots and linearly
interpolate absolute rows, including their throughput, without interpolating
angles or normalizing Stokes fractions. Generated roundoff is directed inward
to the convex passive set; authored inputs are never repaired. Constant rows
retain the original values bit for bit. Lookup is bounded logarithmic work
and allocates nothing.

Geometric `sample_pixel_spectrally_analyzed` and
`sample_pixel_sensor_spectrally_analyzed` resolve each row before the original
scene walk and integrate the measured radiances through the existing XYZ/QE
owners. Steady molecular counterparts retain the original free-flight counters,
scattering Mueller response and shared per-channel paths. Their sequential
counterparts install each wavelength's detector row
at the sensor BEFORE the real lens. The original lens geometry, basis changes,
Fresnel/Beer response, clipping and proposal mass remain authoritative. The
shared lens generator accepts an explicitly resolved row array; no wavelength
borrows the hero's response and no extra analyzer factor is applied afterwards.

`TransientCurvedCameraScene::with_spectral_analyzer` borrows the same source law
for instantaneous and finite-exposure clear, molecular and finite-area paths.
Native sequential raster observers place the row before lens propagation and
retain the lens's original group delay. Public raw-lens continuations refuse a
separately configured constant OR spectral analyzer, preventing accidental
application after the lens. Visible and QE raster owners require the complete
integration band before tracing; unsupported rows never become dark samples.
The existing constant APIs and analyzer-free paths are retained.

This is a local, time-independent detector response, not resolved bulk
birefringence, an angle-dependent polarizer, additional group delay or coherent
spectral phase. Focused G0/G1 regressions cover passive interpolation and support,
independent circular-source spectral powers through XYZ/QE/exposure, orthogonal
channel superposition through a dispersive real lens and curved scene, and
constant-response parity. Native execution results are reported separately.

## Retarded finite-area source connections and finite shutters

`camera::transient::area_lights::TransientCurvedAreaLightCameraScene` couples original finite-surface histories to the native clear-medium area-light camera. Each discovered refracted connection image uses the original camera-prefix plus connection group path, summed in DD before conversion to delay. It samples the original source identity at the actual sampled emitter endpoint and its own retarded time. Complementary camera emitter hits use their original accumulated group delay and MIS weight. Distinct optical images never share an averaged clock.

The adapter borrows exactly the same prepared camera as the area-source owner and retains the native analyzer launch, continuation samples, image discovery, Jacobians, refractive transport, source PDFs, MIS, cancellation, and bounded work. One modulation callback changes only source radiance. Logarithmic range recovery retains representable contributions across attenuation and source scaling. Unbound sources, supported thermal emission, and environment emission remain steady; thermal terms receive neither surface-history scaling nor surface MIS weights. Missing group delay on a contributing bound source, invalid source scale, unresolved clock, or a source-history failure refuses the complete observation.

`sample_ray_with_samples_at_time` returns the complete native radiance and per-image work at a DD observation clock. `sample_ray_with_samples_during_exposure` requests each source's mean over its individually retarded shutter interval on that same frozen transport, then returns mean radiance and integrated radiance [J/(m² sr nm)]. No temporal path sampling, extra image tracing, midpoint source substitution, or additional geometric delay is introduced. Geometry/optical coefficients stay frozen; moving media, pulse reshaping, interface residence time, point-source connections, and molecular-scattering connections are outside this adapter.

Four focused regression groups cover unity/native-sample and work parity, dispersive group clocks at a large epoch, triangular shutter integration on both MIS strategies, distinct refracted-image clocks, and complete refusal for missing connection group delay or source failures. They were authored and formatter/parser checked; native numerical execution was unavailable in the editing environment.

### Retarded finite-area sources through complete Smith glass

`TransientCurvedAreaLightCameraScene::with_smith_dielectric` prepares a
borrowed `TransientSmithDielectricAreaLightCameraScene` over the exact original
histories and full opaque/two-sided dielectric Smith finite-area estimator.
Instantaneous and finite-shutter observations retain every original source
point, polarization channel, direct-source endpoint partition, facet walk,
refractive/Beer response, medium conversion, six RNG domains, cancellation
boundary and source-work budget.

A camera emitter hit uses its original accumulated group path. A finite-area
connection uses the original camera-prefix and connection group paths, summed
in DD by the existing temporal owner before division by c. Its source-history
query is shared by its original physical polarization channels and admitted
before a stochastic fixed-direction response can vanish. The history changes
only its original source radiance; logarithmic range recovery retains the
native weighting, and thermal emission keeps its existing owner.

Qualified sequential beams contribute their original DD lens delay exactly
once to each source clock. Raw sequential observations refuse a separately
configured analyzer, which must precede lens propagation. Finite shutters
request each source's interval mean on the same frozen transport and return
mean radiance and its physical time integral. No temporal sampling, extra
transport, query allocation or inferred group law is introduced by the adapter.

The native direct finite-area proposal and equal-index mixed-measure policy
remain unchanged. The Smith response is local: statistical height is not a
physical flight length, and no interface residence time, pulse reshaping,
moving optics, molecular connection or new refractive-image proposal is
claimed. Two focused native regression groups cover both source sides,
wavelength and lens group clocks, steady path/work parity, pulse energy and
complete source or lens admission failures. These renderer fixtures were
parsed and source/API reviewed; they were not executed in this environment.

## Retarded point-fixture connections

`camera::transient::point_lights::TransientCurvedPointLightCameraScene` binds
`PointEmissionHistory` to the exact borrowed `PhotometricSource` and an explicit
nonzero content identity. Instantaneous and finite-shutter queries observe the
original accepted point-light connection, with its calibrated source direction,
camera and connection group paths, and complete sequential-lens DD delay prefix.
Surface-emission histories retain their original object identities. Unbound
fixtures are steady and require no new group-delay evidence.

The observer scales the native point-light intensity before logarithmic range
admission; a unit history preserves the existing arithmetic. Point fixtures
remain additive delta-endpoint next-event estimates. Clipped, zero-survival and
otherwise noncontributing connections acquire no fabricated source event.
Missing group evidence, source-history failures and cancellation refuse the
whole requested observation. Steady raw-beam and timed sequential-lens methods
retain the native polarization analyzer and optical transport.

These histories assume frozen geometry and optics. Environment, LTE and volume
emission remain steady. Molecular point connections and discovered source
images use the additional adapters below. Three regression fixtures cover separate
surface/fixture clocks, native work parity, timed analyzed lens pulses and
shutters, and admission/refusal cases. They were source reviewed and parsed;
the complete renderer test target was not executed in the editing environment.

## Retarded molecular point fixtures and refracted images

With `volumes`, `TransientCurvedRayleighPointLightCameraScene` and
`TransientCurvedRayleighImageLightCameraScene` apply the same exact-source
histories at native surface and accepted molecular vertices. Prescribed
connections and every deduplicated refracted source image retain their actual
camera-prefix and connection group paths, source polarization, photometric
normalization, and independent source-survival realization. Duplicate search
seeds do not duplicate source events or source streams.

Instantaneous and finite-shutter methods accept ideal rays or complete timed
sequential-lens beams. The DD lens prefix enters the clock once. Original
zero-survival, zero-phase and noncontributing guards precede history evaluation;
null molecular candidates create no invented emission event. A steady
raw-beam launch remains available without requiring group timing. Surface,
flight and source-counter streams, original geometry searches, trace budgets
and cumulative work accounting remain with their native workers.

History/source identity and camera ownership are admitted against every actual
attached surface, volume and image fixture. Unbound fixtures remain steady and
bypass new group eligibility; missing group evidence and original history
failures preserve typed refusal of the entire observation. Frozen geometry and
optics, steady environment, LTE and volume emission remain explicit limits.

Three additional regression fixtures cover heterogeneous null/real candidates
and source stochastic zero, distinct polarized images and shutter means with
exact native trace budgets, and chromatic analyzed-lens timing at a real
molecular vertex. They were parsed and independently source/API/math reviewed;
renderer execution was unavailable in this editing environment.

## Steady coated-lens point-fixture pixels

`camera::raster::point_lights::sequential` adds native steady CIE/XYZ and
multi-channel QE observations for clear point fixtures, molecular point
connections and discovered refracted source images. Each family exposes
ordinary, constant-analyzer and source-supported spectral-analyzer variants.
The eighteen entry points share the original point collectors and complete
steady beam methods; molecular consumers retain all three original random
domains and per-lane source/camera work.

Steady observation requests no group clock, so nontrivial absorbing lens
coatings remain usable through their original full-Stokes response. A detector
row is installed once at the physical sensor before the lens. Existing source
calibration, scene-medium admission, source polarization, Rayleigh response,
source survival and discovered-image accounting remain with the native owners.
The CIE/QE integrators apply original geometry weights, spectral probabilities
and photon conversion once. Adding detector channels creates no additional
transport path. Clipped lanes keep their proposal mass and zero scene work;
missing analyzer/source support and original transport failures refuse.

The final native-owner regression exercises all six clear detector modes and
molecular/image spectral-analyzer CIE/QE observations, positive independent
fixture radiometry, exact native lane/work replay, active-coating timing
refusal, closed-pupil accounting, spectral support and cancellation. Source
and fixtures were parsed and independently source/API reviewed. Numerical
renderer execution remains unperformed in this editing environment.

## Retarded point-fixture CIE and quantum-efficiency pixels

`camera::raster::transient::point_lights` connects the admitted point histories
to the original physical-camera integrators. Twenty-four instantaneous and
finite-shutter entry points cover CIE/XYZ and multi-channel QE observations,
ideal rays and complete timed sequential-lens beams, and clear, molecular and
discovered-image source owners. QE channels reuse the same four spectral
paths, original wavelength support and photon-energy conversion.

Shared point collectors now accept the original launch, pixel payload and
typed transport error without changing their existing public default types.
Surface, camera-flight and source-survival counters retain their three
original Philox domains and connection ordinals. Per-lane fixture radiance,
connection/search work, molecular events, flight draws and source candidates
survive payload mapping into shutter-energy/electron observations.

Timed lens observations install the configured detector row before the lens
and use the complete qualified DD lens prefix. Clipped lanes retain their
proposal mass, zero contribution and zero scene work; they are never replaced.
Clock admission precedes optional lens lanes. An original source/history
failure retains its lane and refuses the complete pixel. No extra source
survival, delta-light MIS weight or fabricated group delay is introduced.
Mean and integrated exposure products use the existing native finish helpers.

Four additional fixtures cover independent chromatic CIE/QE radiometry,
steady bit parity, full-lens clock/analyzer/lane replay, clipped mass,
molecular collisions, deduplicated images, all three random domains, and
finite-shutter energy/electron work. All changed Rust was parsed and
independently source/API reviewed; these renderer fixtures were not executed.

## Retarded molecular finite-area source connections

With `volumes`, `camera::transient::area_lights::participating::TransientCurvedRayleighAreaLightCameraScene`
borrows the original finite-source histories and molecular area-light owner.
Instantaneous and shutter observations retain native surface and molecular
NEE, complementary emitter hits, refracted image discovery, Rayleigh Mueller
transport and independent source-survival tracking. The existing source-clock
calculation observes each image's actual camera-prefix and connection group
paths; null collisions and optical images do not acquire surrogate delays.
All five original random streams, work limits and per-sample accounts remain
with the molecular estimator. Contributing bound sources without original
group data refuse; exact survival zeros need no history query.

Steady raw sequential beams now enter that native molecular area owner without
a group-timing requirement. Timed sequential observations add their qualified
lens DD prefix once and retain the original sensor clock, full Stokes state,
lens absorption and optical normalization. Raw-beam temporal calls require any
analyzer to be installed at the sensor before the lens. Finite shutters use
each source's interval mean on one frozen camera/collision/source realization.
Spatial volume/LTE and unbound emission remain steady. Frozen media and
instantaneous molecular events do not imply fluorescence, moving-medium
frequency shifts, interface residence time or coherent pulse reshaping.

Focused G0/G1 source regressions cover null candidates and independent source
survival, distinct image clocks through an analyzed dispersive lens, complete
pulse shutters, untimed/analyzer admission and original group/history errors.



## Retarded area-light raster and QE exposure observations

`camera::raster::transient::area_lights` adds `sample_pixel_area_lights_at_time`, `sample_pixel_sensor_area_lights_at_time`, `sample_pixel_area_lights_during_exposure`, and `sample_pixel_sensor_area_lights_during_exposure`. Each of the four native wavelength lanes retains its own camera path, refracted source images, original group clock, history evaluation, pupil weight and MIS. The same lanes feed the existing CIE/XYZ or multichannel QE integration; detector channel count adds no transport or source samples.

Steady and transient consumers share the original area-raster lane worker and typed Philox counter addresses. Surface and area-source random domains and native per-lane image, trace, field and query work are unchanged. The defaulted source-error type preserves the existing steady API while temporal failures retain their original lane and typed source cause. Any failed lane refuses the complete observation; surviving lanes are not renormalized or exposed as a partial pixel.

Finite exposures use each source's retarded interval mean on frozen transport, then the existing detector/exposure owners supply mean CIE pixels, mean QE rates, integrated spectral energy and photoelectron counts per m² sr. Pixel etendue, gain, saturation and noise remain separate. The area-sample wrapper preserves original per-lane radiance and work across detector conversion. No extra spectral paths, temporal quadrature traces or delay reconstruction are introduced.

Two focused regression groups cover complete static CIE/QE parity, one-versus-three-channel shared work, large-epoch chromatic camera-plus-connection group clocks, native QE energy/photon conversion, late source failure with deterministic replay, and cancellation. Formatter/parser checks passed; native numerical execution was unavailable in the editing environment.


### Steady coated-lens molecular area-light pixels

The original molecular raster also exposes
`sample_pixel_sequential_rayleigh_area_lights` and
`sample_pixel_sensor_sequential_rayleigh_area_lights`, each with
`_analyzed` and `_spectrally_analyzed` variants. These use the original
steady lens integrators with timing disabled, so coated or absorbing lens
optics need no group-delay model to produce CIE or multichannel QE images.

A detector row is applied at the physical sensor before the full polarized
lens. Native lens transmission, Beer absorption, molecular scattering,
finite-area source images and complementary MIS retain their original
correlation. All six methods reuse the original five sampling streams and
molecular raster worker. Wavelengths clipped by the pupil contribute zero
with the full proposal mass; detector channels share paths and work.
Spectral analyzer support is admitted before tracing, including fully clipped
packets. No source clock or additional coating phase derivative is inferred.

One native fixture specifies positive finite-area radiometry with independent
molecular source survival through absorbing coatings, all six detector modes,
CIE/photon integration and unchanged work. It also checks that the same
active coating refuses an unsupported timed-lens query while steady rendering
succeeds, plus complete clipping, missing analyzer support and cancellation.
The four changed Rust files passed parser and source/API review; this renderer
integration fixture was not executed in the editing environment.

### Timed real-lens molecular area-light pixels

`camera::raster::transient::area_lights::participating` exposes
`sample_pixel_sequential_rayleigh_area_lights_at_time`,
`sample_pixel_sensor_sequential_rayleigh_area_lights_at_time`, and their
`_during_exposure` counterparts. They connect the original molecular
finite-area scene to the original optional real-lens CIE and QE integrators.
The analyzer is applied at the sensor before refracting through the lens.
Each surviving wavelength retains its complete Stokes beam, qualified lens
group delay, molecular camera path, independently tracked source survival and
the emitted source channel's own retarded history.

The shared molecular raster worker retains the original five Philox counter
domains and per-lane work. Clipped wavelengths contribute zero with the original
four-lane proposal mass; channels share the same paths. The original integrators
apply lens geometry, spectral density and QE/photon weights once. A typed
source-clock failure identifies its lane and refuses the whole observation.
Clock admission and cancellation are checked even when every pupil lane clips.

Finite shutters integrate each history on its frozen retarded interval without
adding temporal traces. The existing exposure owners convert the completed mean
to spectral energy and electrons per m² sr once, retaining all molecular and
source-work counts. These remain frozen-media, instantaneous-scattering source
histories; spatial volume/LTE emission remains steady.

Two native regression fixtures independently specify source-survival radiometry,
double-double source clocks, replay of all five sampling streams, CIE/photon
integration, partial and complete clipping, and source/clock/cancellation
refusals. The changed Rust files passed parser checks and source/API review;
the renderer integration fixtures were not executed in this editing environment.

### Exact scalar-transfer opacity and physical density calibration

The volumes feature exposes render_transfer_transport for orthographic
piecewise-constant scalar fields. Each cell contributes its analytically
integrated source radiance through the actual foreground optical depth.
TransferTransportPixel retains emitted linear RGB, transmittance, optical
depth and independently evaluated expm1 opacity. Its background compositor
uses L_out = L_emitted + T L_background; source and background radiance remain
representable behind subnormal unit-radiance transmission where possible.
There is no sampling seed, null-collision budget or opacity cutoff in this
path. Work is proportional to admitted pixel/cell intersections, independent
of extinction magnitude. The existing named Woodcock APIs retain their
sampling semantics.

With `tracer,volumes`, `volumes::opacity::TransferRayGrid` admits one immutable
borrowed scalar field and transfer function, then integrates arbitrary clipped
and perspective rays with the shared native cell DDA. Non-unit ray parameters
are converted to actual world length; emission accumulates in camera order.
Admission scans each field value once, polling cancellation every 1024 values;
each ray allocates nothing and has an explicit crossing limit. Cancellation and
budget refusal return no partial result, and optical thickness never changes
the work limit. Orthographic and arbitrary-ray paths share the analytic cell
solution, retaining representable source-scaled light when separate tiny opacity
or foreground transmission underflows. Transfer rendering semantics version 3
records the changed tiny-source arithmetic. This remains cellwise-constant
prescribed RGB source plus direct extinction, without inferred LTE, in-scattering,
refraction, continuous-field interpolation or pixel-area filtering. Regressions
exercise camera-generated perspective rays, analytic ordered layers, clipped and
rescaled rays, existing orthographic parity, extreme radiance and refusals.

With `tracer,volumes`, `volumes::opacity::render_transfer_camera` renders one
borrowed transfer snapshot through `PhysicalCamera` perspective or ideal
finite-aperture rays. One sample uses pixel/lens centres; larger counts use
independent uniform pixel/pupil variates from per-pixel Philox streams. Returned
pixels retain mean emitted linear RGB, mean transmission and stable mean opacity.
Effective depth represents negative log mean transmission, using thin-opacity
and thick log-sum evaluation. Misses remain transparent, and normalized
pupil/radiance averaging applies no exposure or aperture gain. Pixel/sample
admission and actual whole-image cell-crossing limits fail transactionally with
cancellation. Existing spectral/weighted-camera refusals remain binding.
Per-pixel background composition assumes a background constant over the sampled
footprint; prescribed scalar RGB source/extinction does not infer LTE,
in-scattering or surface refraction. Physical regressions compare perspective
slabs and defocused/refocused opaque-strip coverage against analytic references.

render_density_transmittance borrows the original density buffer and uses
explicit kg/m³ per field unit, metres per grid unit and m²/kg mass absorption.
It produces monochromatic direct-beam transmission without allocating an
intermediate coefficient field. No wavelength response, density conversion,
surface coverage, scattering, refraction or material calibration is inferred.

Both orthographic renderers admit output pixels, input cells and full cell visits before
integration, validate finite geometry/fields, and use the existing integer
pixel-centre ownership. Transfer opacity uses expm1 for thin cells, optical
depth and emission use compensated sums, and extreme positive products have
a logarithmic range-recovery path. Analytic regressions cover ordered slabs,
radiance balance, thin/opaque/vacuum limits, spatial subdivision, translations,
SI unit changes, bright sources behind thick foregrounds, and admission
refusals. These are representation-level arithmetic checks; this piecewise-constant
path makes no claim of continuous-field preintegration, lateral antialiasing,
experimental optical calibration, scattering transport or measured performance.

`render_density_transport_box`, `render_density_transmittance_box` and
`render_transfer_transport_box` integrate the complete orthographic pixel
footprint over the original cellwise columns. Exact integer overlap areas work
at nondivisible resolutions and arbitrary world translations. Transport is
evaluated before filtering: a half-clear, half-opaque footprint transmits half
its background, and layered source emission retains its original correlation
with foreground attenuation. This is lateral filtering of the represented
columns, not reconstruction of an unknown interface inside a cell.

The box and physical-camera paths share one compensated light accumulator.
It averages premultiplied radiance, transmission and thin opacity separately,
retaining finite effective log attenuation for bright backgrounds behind thick
fields. Box rendering needs only output storage and constant scratch; it does
not allocate a converted field, high-resolution image or supersampling buffer.
Its exact column-cell work is `(nx+w-gcd(nx,w))*(ny+h-gcd(ny,h))*nz`, admitted
before field scans and allocation against the existing cell-visit budget.
These orthographic entry points remain synchronous without a cancellation hook.
Eight focused G0/G1/G3 regressions cover analytic coverage, nondivisible 2-D
footprints, layered emission, SI/refinement invariance, thin and opaque limits,
radiance balance, exact work admission and integer coverage. They establish
the represented box average up to floating-point error, not a measured speed
or optical-data accuracy claim.

### Wavelength-resolved physical opacity

`render_density_spectral_box` and `render_density_spectral_transmittance_box`
extend the same borrowed density model to 1–16 explicitly supplied increasing
vacuum wavelengths and matching mass absorption coefficients. The wavelength
samples retain their original order; there is no spectral interpolation,
implicit RGB conversion, bin-width inference or material-name calibration.
A single original-column density reduction supplies every wavelength, and
exact integer footprints average each wavelength's transported intensity.
The reduction retains a finite scale when column mass would overflow, so a
small coefficient can still produce finite optical depth. SI units, physical
thickness, thin opacity and effective log survival retain the same meanings as
the monochromatic owner. The full image preserves wavelengths and permits
stable attenuation of a constant-footprint spectral background.

The direct-transmission entry point writes pixel-major, wavelength-fastest
scalars without a full-statistics intermediate image, RGB emission work or
opacity/log-survival evaluation. Both paths share geometry and density reads;
neither constructs a converted field or one image per wavelength. In addition
to the existing exact cell-visit budget, `max_spectral_evaluations` admits
`(nx+w-gcd(nx,w))*(ny+h-gcd(ny,h))*wavelengths` before scanning the field or
allocating output. This also bounds retained sample count. These synchronous
absorption views do not infer scattering, emission or unresolved subcell shape.
Eight spectral regressions cover SI slabs, wavelength-dependent partial
coverage, refinement, nondivisible footprints, monochromatic parity, finite
depth behind an overflowing mass sum, thin/opaque backgrounds, flat layout and
work/input refusal. The targeted native volumes-enabled source build passed
117 tests, including these eight; this is not a full-workspace or tracer gate.

### Accepted material-state absorption

With both `material-state` and `volumes` enabled,
`material_state::absorption::MaterialAbsorptionBinding` prepares physical mass
absorption from the original accepted material state, optical source, and an
explicit scalar mass-density property. The accepted and optical states must
agree on material, source card, complete query point, and axis quantity
conventions. The binding borrows the original accepted bundles and density
receipt; it does not construct a replacement state or infer a density basis.

Native visible dielectric intensity absorption is divided directly by intrinsic
material density. Visible conductor and generic spectral complex-index sources
evaluate their original extinction index at each requested wavelength before
converting it with `kappa_a = 4*pi*k / (lambda_m*rho_material)`. In particular,
sampling converted coefficients at source knots and then interpolating them
would generally change the source law and is not performed. The original
closed wavelength support and spectral wavelength/property-key interpretation
remain part of the binding. RGB-derived attenuation is refused.

`MaterialAbsorptionSamples` prepares one through sixteen finite, positive,
strictly ordered wavelength nodes into a fixed coefficient array. It retains
the requested wavelength slice and accepted/source/binding identities, and
`as_density_spectrum` borrows those samples directly into the lean density
renderer. The coefficient slice can also supply indexed material palettes and
prepared mass-mixture inputs. Preparation introduces no per-ray constituent
evaluation, converted field, or interpolated replacement material law.
Intrinsic material density normalizes the mass coefficient; the rendered
field's density-unit and length-unit scales remain explicit caller inputs.

Changing material conditions requires refreshing the original accepted and
optical sources and rebinding. Stale state pairs, absent or incorrectly typed
density, out-of-support wavelengths, malformed sample nodes, and positive
coefficients that cannot be represented are refused. Ordinary complex-index
conversion takes a short arithmetic path; a logarithmic fallback handles
recoverable intermediate range loss, while representable positive subnormal
results remain valid. Particle scattering uses the same extracted scalar
density admission predicate without changing its existing admission or receipt
semantics.

Five analytic tests cover original-law sampling through actual density opacity,
temperature/density refresh and stale-state refusal, native dielectric
absorption and RGB refusal, closed support for lossless sources, and extreme
coefficient conversion. These material-state tests received pinned Rust parser
validation and two source/physics reviews; they were not compiled or executed
because the material-state/tracer dependency closure is unavailable in this
workspace. The separately executed lean volume tests do not cover this feature
gate. This bridge adds no measured optical dataset and makes no scattering,
interface-reflection, or emission inference from absorption alone.

### Finite-band dielectric absorption sources

`MaterialAbsorptionBinding` accepts the existing finite-band tabulated
dielectric and Sellmeier state-point families directly, with a unified
`try_spectral_optical` constructor for the existing explicit optical enum.
The material owner's allocation-free `absorption_per_m_at` query evaluates
original tabulated alpha; analytic dielectrics reuse their original alpha
evaluator. Both divide intensity absorption directly by intrinsic material
density. Conductor dispatch retains the original n/k binding and identity.

Absorption-only Sellmeier observations use the independently declared absorption
support without evaluating or truncating to the refractive-index fit. Refraction
consumers still require that separate index domain. SI wavelength support must
have a finite, positive, ordered nanometre representation. Closed-band admission
precedes endpoint roundtrip correction, so no out-of-band input is extrapolated.
Source identity, density receipts, complete state matching and refresh behavior
remain those of the original accepted claims. This adds no new measured data.

### Accepted constituent opacity mixtures

`MaterialAbsorptionComponent` pairs an original accepted absorption binding
with explicit kg-of-constituent per kg-of-reference-mass loading.
`MaterialAbsorptionMixture` borrows bounded component declarations, retains
original source identities and support, and never normalizes weights.
`MaterialAbsorptionMixtureField::try_new_mixtures` samples every original law
at the exact immutable wavelength nodes. Its mass absorption is
`kappa_mix = sum(w_i*kappa_i)`; accepted LTE source radiance is
`sum(w_i*kappa_i*B_lambda(T_i))/kappa_mix`. Constituent temperatures remain
distinct. The common-temperature case uses the original single Planck value.

The same prepared `MaterialAbsorptionField`, cameras and transport methods
serve individual bindings and mixtures through a default model parameter.
Existing single-binding API, numerical preparation and identity ordering are
preserved. Final storage is O(recipes times nodes); source evaluation,
constituent loops and Planck evaluation occur during preparation only.
Preparation charges every component/node pair, including unused models and
zero weights, before allocation, and retains the existing 64-component cap.
Original compensated positive sums and range-aware weighted products are
shared unchanged. A positive coefficient that cannot be represented refuses.

Zero-weight components still participate in source admission, support and
identity, but need no temperature and create no absorption or emission.
An all-zero recipe is transparent. This is additive homogeneous independent
absorption on an explicit reference-mass basis; it does not imply effective
refractive indices, interfaces, scattering, chemical interactions or thermal
equilibration. New accepted states or wavelengths require new preparation.

Four G0/G3 tests cover exact original spectra, independent 300/400 K Planck
sources in native camera radiance, beam power closure, unnormalized loadings,
common-temperature and single-component parity, zero weights, complete work
budgets, source support, coefficient range and cancellation.

### Material beam deposition and thermal input

`MaterialAbsorptionField::deposit_beam_power` accepts incident watts already
integrated over beam area and the caller's spectral widths at the field's
fixed nodes. Individual materials and mixtures reuse the existing native
DDA and conservative `deposit_prepared` kernel; no new traversal or partition
owner is introduced. Original density cells, accepted identity, source-band
intersection, finite endpoints, traversal/output limits and cancellation
remain authoritative. The shared native kernel retains its existing
compensated cumulative-depth implementation.

The result is the existing `BeamPowerDeposition`, directly admitted by
`OpticalHeating::try_new` with an explicit optical-to-thermal cell map and that
owner's original source-support requirements. It reports absorbed incident
watts even with accepted LTE emission enabled. Cooling, source/state refresh,
thermalization assumptions and material feedback remain explicit; the
deposition call does not silently advance temperature.

Three native G0/G3 tests exercise original source coefficients, per-cell and
reverse finite-beam accounting, incident-only heating, thin/opaque representable
power, vacuum, budgets and cancellation. One `optical-heating` integration test
feeds original thermal cells and the actual native enthalpy clock while
retaining boundary heat and energy accounts. These four tests and the four
mixture tests passed pinned Rust parsing and independent source review.
They were **not compiled or executed** because the material/tracer/heating
dependency closure is unavailable here; no runtime claim follows for that
feature configuration.

### Accepted material opacity in native transport

`material_state::absorption::native::MaterialAbsorptionField` prepares one
immutable packet of one through sixteen wavelength nodes per original material
binding. It retains the original cell palette and borrows the existing admitted
SI geometry/density carrier, following the native prepared-field convention.
The carrier's previous absorption and temperature sources are explicitly not
used. Coefficients evaluate the original material law at each exact node;
there is no off-node query or interpolation of converted mass coefficients.

Preparation storage is proportional to models times wavelengths, not cells
times wavelengths. Tracing performs neither material-source evaluation nor
Planck evaluation. Original material identities, cell mapping, exact nodes and
explicit emission mode enter the prepared identity. The source-support
intersection records scope without widening any supplied source.

`MaterialAbsorptionEmission::None` needs no accepted temperature.
`AcceptedLte` uses each absorbing material's own unique accepted Kelvin
temperature and the native Planck law, linked to its same absorption
coefficient. It introduces no averaged temperature or independent emissivity.
The extracted temperature predicate is shared with particle optics without
changing existing particle admission. Intrinsic density normalizes kappa;
the borrowed grid still declares selected-material mass per bulk cubic metre.

`transmittance`, `transport`, `render_camera` and
`render_camera_with_background` reuse one original native cell integral and
the original camera sampler. DDA crossings, compensated optical depth,
emission arithmetic, opaque limits, finite endpoints, cell accounting and
cancellation remain shared. Per-ray illumination is composed before pixel/pupil
averaging. Only mapped positive-density cells determine vacuum or emission;
an unused absorbing model cannot create haze.

Preparation admits model/node budgets and map length before output allocation,
checks every original source and polls cancellation. Source, temperature,
arithmetic, traversal or camera refusal publishes no partial field/image.
This is the existing non-scattering, index-one continuum model, without
surface reflection, refraction, unresolved spatial coverage or thermal-state
evolution inferred from absorption.

Two material-owner tests, four finite-band renderer tests and five native-field
tests cover source-law interpolation, IR alpha transport, independent Sellmeier
support, refreshed state/density/wavelength coordinates, layered 300/400 K
emission, reversed layer order, finite endpoints, temperature-free absorption,
unused models, budgets and cancellation. Pinned Rust parsing and source/API/
physics review passed, including independent review of the native extraction.
These material-state/tracer tests were not compiled or executed because their
complete dependency closure is unavailable. The separate 141-test lean volume
run does not establish execution coverage of these feature configurations.

### Indexed material opacity

`IndexedDensityAbsorptionSpectrum` borrows the original x-fastest density
field, one `u32` material/state index per cell, and a material-major,
wavelength-fastest table of mass absorption in m²/kg on common 1–16 positive
increasing vacuum wavelengths. The indexed full-image and direct-transmission
entry points integrate each cell's density, selected coefficient and physical
length before averaging transported intensity over exact pixel footprints.
Distinct layers and adjacent materials retain their original assignments;
there is no converted extinction grid or inferred mixing law.

Indexed spectral work counts every intersected original cell times wavelengths;
uniform spectral work remains overlapping columns times wavelengths and keeps
its single density-column reduction. A separate `max_materials` bounds every
palette row, including unused rows. The complete palette and every map index,
including zero-density cells, are validated before output allocation. Outputs
retain the existing pixel-major spectral layout, SI scaling, thin opacity and
effective log-survival semantics. Six regressions cover independent layered
Beer–Lambert references, colored subpixel coverage/refinement, uniform and SI
parity, range limits, palette/map refusal and exact spectral work admission.

### Prepared mass-mixture opacity

`prepare_density_absorption_mixtures` prepares one point-sampled absorption row
per declared composition for the existing indexed renderer. Constituent input
is component-major/wavelength-fastest in m²/kg of constituent; supplied weights
are mixture-major/component-fastest in kg of constituent per kg of the density
field's declared reference mass. Output is mixture-major/wavelength-fastest in
m²/kg of that reference mass. Total-density mass fractions and carrier-relative
loadings retain their explicit meaning; weights are never normalized.

The independent-constituent model sums mass-weighted absorption at each exact
supplied wavelength. Homogeneous composition within one cell remains distinct
from the renderer's spatial average of transported intensities across different
rays. No index, phase, temperature, interpolation, material identity or optical
measurement authority is inferred. Mixtures are prepared once; rendering keeps
the existing indexed hot loop and never evaluates constituents per ray.

Component count, mixture count and the full mixture/component/wavelength
product are admitted before complete input scans or output allocation, including
zero weights and unused components. All inputs are finite and nonnegative.
Positive coefficient products that underflow to zero and overflowing prepared
coefficients refuse explicitly, because later density/path scaling could make
the lost absorption observable. Only the final palette is allocated.

The browser `prepare_vessel_spectral_absorption_mixtures` function and wasm
export forward this same owner with limits of 64 component spectra, 4,096
recipes and 1,048,576 total evaluations; invalid requests return an empty vector.
The targeted volumes-enabled native source build passed all 134 tests, including
four mixture regressions and the previously typechecked indexed/observer cases.
This is separate from material-state, tracer and WASM feature validation.

### Spectral illumination and visible observation

`SpectralOpacityImage::approximate_transmitted_xyz` and its per-image-background
variant project the original wavelength-resolved transmission under an explicit
incident radiance spectrum. `AbsorbingCameraImage` and `LuminescenceCameraImage`
provide corresponding `approximate_composited_xyz` methods that include native
foreground emission. Uniform or exactly pixel-matched backgrounds are supported;
each background spectrum must be constant within its own filtered pixel/pupil
footprint. Averaged opacity cannot reconstruct unresolved background correlation.

The lean and thermal observers consume radiance per nanometre plus explicit
nonnegative quadrature widths in nanometres. Luminescence consumes integrated
bin radiance and uses unit weights; it never invents bin widths. The shared
allocation-free lane accumulator prepares the existing analytic CIE fits once
per image and projects emission and transmitted illumination separately into
compensated XYZ sums. It avoids an overflowing intermediate spectral sum,
retains original factors when an observer-weight product underflows, and uses
retained optical depth to recover finite bright-background transmission. Only
the resulting XYZ image is allocated. Physical-camera cancellation returns no
partial image, and sampling/native transport behavior is preserved.

These are approximate relative visible observations using the existing
380–780 nm observer and `1/integral(y_bar)` normalization, not calibrated
luminance, broadband radiative power, exposure or display RGB. Unsupported
wavelengths, malformed weights/background layout, negative/nonfinite radiance
and unrepresentable XYZ are refused. Primary spectra remain available outside
the visible band. Seven lean/common regressions are typechecked by the targeted
native volumes-enabled metadata build; five additional physical-camera
regressions have parser and independent source review but have not run.

### Spectral beam power on original density cells

The `volumes` feature exposes `deposit_density_spectral_beam` and its indexed
material counterpart. They integrate explicit incident spectral irradiance in
W/(m² nm), caller-supplied nonnegative quadrature widths in nm, and each original
column's physical XY area. Sources may be uniform or original-column resolved;
propagation may enter either Z face. This is collimated, normal-incidence,
independent absorption: camera resolution, RGB observation and exposure do not
define incident energy.

The result retains exact nodes and weights, per-node incident/absorbed/transmitted
watts, original-cell absorbed watts, and original-column transmitted watts.
Cell and column maps sum wavelengths and retain original X-fastest ordering.
The original density, material index, SI density scale and SI length scale
participate in every optical depth; coefficients and illumination are never
spatially averaged before the energy calculation. Compensated cumulative depth
and `expm1` retain thin-cell absorption and deep-column transmission; logarithmic
fallbacks preserve representable watts when intermediate optical depth, area or
transmission underflows. Nonrepresentable power refuses.

Grid cells, complete cells-times-wavelength work, indexed material rows, and
cell-plus-column output values are admitted before allocation. The report owns
one final buffer: three global power totals, the cell map, then the column map.
`into_power_buffer` transfers that allocation to downstream consumers. Fixed-size
per-node totals require no additional heap allocation.

G0/G3 checks cover independent layered Beer integrals, both directions, original
material/source correlation, energy accounting, spatial refinement, unit changes,
10,000 sub-ULP optical depths, extreme representable powers, zero quadrature and
input/work/storage refusals. The real pinned-compiler `fs-render` test target with
`volumes` passed **148 tests**, including all seven new beam checks. This establishes
the tested finite spectral quadrature, not broadband integration accuracy or a
measured performance target. Deposited watts are absorbed incident power; this
owner does not add scattering, subtract LTE cooling, infer a thermalization
fraction, or advance temperature.

### Direct illuminated density observations

`SpectralIllumination` supplies either one explicit spectral radiance or one
spectrum per original XY density column, with X fastest and wavelength innermost.
The column layout is independent of output resolution. Light comes from behind
the minimum-Z face; each supplied spectrum is constant over its original column
footprint. Inputs are finite nonnegative W/(m² sr nm). No source, exposure,
emission, scattering or subcolumn geometry is inferred.

`render_density_spectral_radiance_box` and its indexed-material counterpart
integrate each original column, apply that column's incident spectrum, then
average transported intensity over the exact integer pixel footprint. The
result is pixel-major, wavelength-fastest radiance. This retains illumination
and opacity correlation that separately averaged backgrounds and transmission
would lose. Uniform material spectra retain one shared density reduction per
column; indexed spectra retain the original ordered per-cell coefficient sums.

The corresponding `render_density_spectral_xyz_box` functions project directly
to one XYZ triple per pixel, using explicit nonnegative nanometre quadrature
weights, the existing analytic 380--780 nm CIE fits and relative
`1 / integral(y_bar)` normalization. Spatial weights enter the shared observer
before range-losing intermediate rounding. A representable pixel result can
therefore survive both an overflowing unweighted column and underflowed
intermediate transmitted radiance. Raw spectral output retains UV/IR support;
XYZ is an approximate visible observation, without display or luminance claims.

A single statically dispatched column traversal owns admission, exact footprint
geometry and absorption for the full-statistics, transmission-only, radiance
and XYZ paths. Only the selected final image is allocated. Radiance uses
compensated online means; XYZ uses compensated weighted sums. All original
geometry/spectral budgets apply before background scans and output allocation;
invalid illumination, palettes, observer weights or arithmetic return no partial
image. The original scalar/opacity arithmetic and traversal order are retained.

Seven analytic tests cover correlated backgrounds, layered indexed media,
nondivisible footprints, uniform-observer parity, extreme optical ranges,
spatial refinement and refusals. The targeted native `volumes` configuration
was compiled against the real fs-math/fs-rand libraries and executed:
141 tests passed, zero failed. This lean configuration does not execute the
material-state, native tracer camera or browser consumer feature tests.

### Native and browser pour beam-power ports

Uniform and indexed `PourSpectralOpticalView` owners expose `deposit_beam`;
`PourOutcome` provides corresponding spectral beam entry points over its
accepted mass field. They return the shared core's typed
`SpectralBeamDeposition`, preserving original cell order and per-node power
accounts. The original SI lattice spacing defines physical area.

The `fs-wasm` functions and browser exports `deposit_vessel_spectral_beam`
and `deposit_vessel_indexed_spectral_beam` consume that same core. Callers
supply W/(m² nm) irradiance, explicit nm quadrature widths, a uniform or
original-column layout, and the entering Z face. The returned buffer is
`[incident_total_W, absorbed_total_W, transmitted_total_W, cell_absorbed_W...,
column_transmitted_W...]`; it transfers the core allocation without another
map copy. Browser admission caps cells and complete spectral work at 1,048,576,
cell-plus-column map values at 2,097,152, and indexed material rows at 4,096.
Camera resolution and visible-light observation are absent from these APIs.

Three native and two browser G0/G3 regressions cover independent
22.4 = 13.6 + 8.8 W accounting, original-cell and spectral totals, indexed
sources, lateral refinement, density-unit changes, UV/IR nodes, flat output
and source/weight/work/storage refusals. These consumer tests passed pinned
Rust parsing and source review but were not compiled or executed in this
environment; the shared `volumes` core separately passed 148 tests.

The frozen two-dimensional pour has no accepted caloric/enthalpy state.
These ports expose physically dimensioned absorbed incident power for an
explicit thermal coupling; they do not invent temperature evolution or infer
an optical-to-thermal map.

### Illuminated native and browser pour frames

`PourSpectralOpticalView` and `PourIndexedSpectralOpticalView` now expose
`render_radiance` and `render_xyz`, with matching `PourOutcome` entry points.
They borrow the captured original mass buffer and use its explicit SI density,
spacing and extruded slab depth. Supplied illumination is uniform or indexed
by original XY column, independently of output image resolution. The shared
column kernel applies source radiance and absorption before footprint filtering.

Browser functions `render_vessel_spectral_radiance`,
`render_vessel_spectral_xyz`, `render_vessel_indexed_spectral_radiance` and
`render_vessel_indexed_spectral_xyz` expose the same physical calculation.
The explicit `background_per_column` flag selects one wavelength spectrum or
`nx*ny*wavelengths` source values; shape never silently changes meaning.
Radiance remains W/(m² sr nm). XYZ requires explicit nanometre integration
weights and retains the shared approximate visible observer and normalization.
No illuminant, exposure, RGB transform, scattering or hidden liquid spectrum is
introduced. XYZ is flattened without copying its final allocation.

The browser admits at most 1,048,576 pixels, original cells, cell visits and
spectral evaluations, and at most 4,096 indexed material rows. Any invalid
input or exceeded limit returns the existing empty-vector refusal, with no
partial image. Native callers receive the original typed optical errors.

Six consumer regressions cover source/opacity covariance, uniform and indexed
paths, output-resolution independence, explicit XYZ widths/layout, native/browser
parity and invalid inputs/work limits. These received pinned Rust parser and
source/API review only; the full fs-vessel/fs-wasm dependency configurations
were not compiled or executed here. The 141 passing lean renderer tests cover
the shared transport kernel, not these consumer feature configurations.

### Correlated spectral camera illumination

`render_absorbing_camera_with_background` and
`render_luminescence_camera_with_background` evaluate caller-owned incident
spectra once at every actual admitted camera ray, including its sampled pupil
origin and direction. The callback receives the original wavelengths and
execution context, and returns finite nonnegative spectral radiance plus the
nonnegative endpoint distance in metres. Positive infinity denotes distant
illumination; zero is the native zero-path case. Native counted transport clips
absorption and emission at finite endpoints before composing that ray's
foreground and transmitted background.

The resulting typed radiance images average E + T*B over the original sample
stream, preserving illumination/opacity covariance that independently filtered
background and transmission images lose. Thermal/absorption radiance remains
W/(m² sr nm); luminescence remains integrated W/(m² sr). Their visible observation
methods reuse the shared XYZ owner with explicit nm quadrature or unit bin
weights respectively, without a second attenuation factor. Field identities
describe the accepted medium; callback-owned illumination is independent.

One private statically dispatched sampler preserves original kernel IDs, draw
order, camera refusals, cancellation and whole-image cell work. Full-statistics
images retain their original accumulation. Radiance-only outputs use compensated
online means and skip unused averaged opacity, transmission and logarithmic
survival calculations. Native per-ray log depth still preserves bright opaque
backgrounds. There is no additional grid traversal or intermediate spectral
image. Callback work is caller-owned and bounded in invocation count by the
admitted camera samples. Every composed ray must be representable; any callback,
transport or cancellation refusal returns no partial image.

Five regressions cover a defocused half-pupil covariance counterexample,
finite/zero endpoint clipping, callback and budget refusals, constant-background
parity, integrated-bin units and direct visible projection. They have pinned
parser and independent source/physics review; the tracer dependency closure
has not been compiled or executed in this environment.

### Continuous scalar-field transfer preintegration

`volumes::preintegrated` adds bounded source/extinction integration for linear
scalar segments and orthographic rendering of borrowed node-centered grids.
The renderer reconstructs the scalar field before applying the existing
transfer function, splits segments at every crossed transfer knot, and
integrates optical depth analytically for that reconstruction. Node spacing
defines an extent of `(dims - 1) * spacing`; no extra cell thickness is invented.
Narrow transfer peaks remain visible even when both sampled endpoints are
transparent.

Constant-source and constant-extinction source terms have analytic paths.
Remaining source transport uses deterministic adaptive 8/16-point Gauss
quadrature over cumulative absorption probability, resolving optically thick
ramps without stepping through tiny physical mean free paths. Pixel, node,
segment, transfer-interval and quadrature limits apply to the whole call.
Failed admission or exhausted work publishes no partial image.

G1/G3 controls cover analytic Beer-Lambert slabs, source gradients, thin layers,
narrow peaks, thick-ramp asymptotics, an independent physical-distance Simpson
integral, endpoint reversal and node-grid refinement. Returned source-error
values are quadrature estimates and exclude roundoff, reconstruction and
transfer-data errors. This original orthographic entry point does not claim arbitrary-ray
reconstruction, lateral antialiasing, scattering, refraction, measured optical
authority, or a per-call cancellation hook. The separate continuous-ray and
camera entry points below add their own bounded reconstruction and cancellation. The source convention remains `j = sigma * S`.
See [Engel, Kraus and Ertl (2001)](https://doi.org/10.1145/383507.383515)
for preintegration and interpolation-before-classification.

With `tracer,volumes`, `PreintegratedRayGrid` transports the original
node-centered trilinear field along clipped physical rays. Each native cell
crossing becomes a cubic Bezier scalar profile; its convex hull and all
intersected transfer intervals bound classification error, so endpoint or
midpoint agreement cannot conceal an interior opacity feature. Bounded adaptive
subdivision admits linear chords for the existing knot-resolving preintegrator.

Returned depth and source reconstruction bounds exclude floating-point and
source-data uncertainty and remain separate from Gauss source-error estimates.
The emitted-source bound includes changed foreground attenuation via maximum
source-profile deviation plus maximum source radiance times integrated
extinction-error bound. This is controlled cubic reconstruction, not an exact
integral or a certified floating-point enclosure.

`opacity::render_preintegrated_camera` reuses the original pupil/raster sampler
and radiance/transmission averages. Maximum sample depth error bounds effective
pixel depth; source bounds and estimates average. These do not bound pixel or
pupil sampling error. Cell, subdivision, transfer-interval and quadrature
budgets are shared across the complete image; cancellation or refusal exposes
no partial output. The model supplies continuous prescribed RGB source and
extinction, without inferred LTE, in-scattering or refraction.

Focused physical controls include analytic cubic depth, narrow interior opacity
peaks, independently integrated varying source radiance, grid/unit refinement,
clipping, finite-aperture camera parity, and complete-image work limits.


## Vacuum Doppler observations of thermal emission

`spectral::doppler::VacuumDoppler` prepares an allocation-free local Lorentz
transform from a physical photon's world propagation direction and separate
emitter and receiver SI velocities. Both velocities must be finite and strictly
subluminal. The direction points from emission toward reception. The transform
retains aberrated photon directions in both endpoint frames, the observed/source
frequency ratio D, reciprocal wavelength mapping and the solid-angle Jacobian.
It includes transverse time dilation. Spectral radiance per vacuum nanometre
transforms as `I_receiver(lambda) = D^5 I_source(D lambda)`; the corresponding
solid-angle ratio is D^-2 and must not be applied to radiance a second time.

With `tracer`, `thermal_surface::doppler::ObservedThermalSurface` borrows the
original spectral thermal surface or retains the original gray source value.
Its explicit `try_new_directional` also borrows the original rough EON source,
using its actual source-rest-frame outgoing cosine and facet absorptance at the
mapped source wavelength. The direction-free constructor still refuses rough
sources. Zero shift delegates to the unchanged directional source law, preserving
its radiance bits. No hemispherical average or extra projected cosine is used.
Its emissivity is evaluated at the mapped SOURCE wavelength. Actual temperature,
source identity and finite-band law remain unchanged. The equivalent Planck
evaluation `epsilon(D lambda) B_lambda(lambda,D T)` avoids prematurely
underflowing the rest-frame Wien tail; D T is only directional apparent color
temperature. Closed finite source support is mapped into the receiver frame;
out-of-band queries are refused rather than extrapolated.

`with_sensor` prepares `ObservedThermalDetector` only after checking the complete
original `QuantumEfficiencySensor` band against shifted source support. Its
`sample` uses the detector's existing receiver-frame wavelength proposal and
channel weights. It returns receiver and source wavelengths, physical spectral
radiance, QE-weighted channel radiance and collected-electron rates per m² sr.
Photon energy is h c / lambda_RECEIVER. No extra PDF Jacobian is introduced for
a proposal already defined in receiver wavelengths. Averaging draws estimates
the observed channel integrals. Per-query storage is fixed and no spectrum,
response grid or source table is resampled or allocated.

These are local directional vacuum observations before visibility, projected
area, etendue and exposure. The caller supplies the physical propagation ray and
endpoint velocities at emission/reception; the adapter does not solve retarded
moving geometry, transform normals, infer dielectric/moving-medium shifts,
propagate polarization, debit thermal storage or claim gravitational shifts.
Analytical regressions specify longitudinal/transverse motion, aberration and
inverse boosts, D^5 spectral/D^4 integrated scaling, original source identity,
finite-band admission, detector-frame photon conversion, nonuniform receiver
proposals, stationary parity and amplified Wien-tail range recovery.


## Directed transport of luminescent volume emission

`scattering::luminescence::LuminescenceTransferOperator` supplies dimensionless
directed transfer fractions for isotropic, uniform emission within explicitly
selected original cells. It reuses native exact extinction and implicit-capture
histories, supports uniform, palette and cell-local optics, retains self-absorption
and every original-cell receiver, and separates escaped from unresolved radiant
energy. Unlike LTE exchange, it neither weights emitters by absorptivity nor
symmetrizes directed probabilities. Fractions apply to bin-integrated watts or
joules without wavelength quadrature or a second solid-angle factor.

Immutable optical/absorption identities and every positive-density material cell
bind later thermal use. Independent deterministic emitter/lane/history streams
and global segment, material, emitter and directed-link limits bound work.
Public-report admission checks sorted mappings, nonnegative finite estimates,
receiver/total absorption and unit-energy closure against an explicit fraction
allowance. Fixed-sample standard errors are diagnostics, not confidence bounds;
tails remain radiation. Source kinetics, reexcitation, heat conversion and
wavelength redistribution belong to the caller. Transport is quasistatic,
elastic, index-one and open-boundary.

Seven focused regression groups specify transparent/opaque limits, independent
finite-box and neighbor-absorption integrals, heterogeneous directed transfer,
scattering tails, deterministic replay and conservation/budget/cancellation
refusals. Pinned formatter/parser checks passed; native test execution was
unavailable in the shared editing environment.


## Luminescent-light absorption on the native thermal timeline

`LuminescentOpticalHeating::try_new_with_emission_transport` couples the directed
emission operator to the original finite-lifetime excitation stores and native
enthalpy cells. An `AbsorbingLuminescenceEnvironment` explicitly declares passive
emission-band absorbers. Every generated radiant joule is assigned once to its
original receiving cell, escape, or a retained unresolved tail; self-return heats
the emitter without deducting excitation twice. Total emitted energy remains a
diagnostic, separate from external energy sinks. This is parasitic absorption,
not an inferred fluorescence recycling channel; those mechanisms are physically
distinct (Abebe et al., https://arxiv.org/abs/1804.02230, Fig. 1).

The coupling uses analytically integrated emission over each accepted interval,
not endpoint power multiplied by duration. It preserves exact pulse boundaries,
dark-afterglow heating, native implicit retries and transactional rollback.
Immutable prepared transfer data is shared across observer clones, with original
receiver IDs compiled once. All positive-density optical cells carry explicit
native thermal bindings and temperature-validity domains, including unpumped
receivers and scatterers. The complete absorption/scattering/phase identity is
checked; scattering at pump wavelengths refuses the absorption-only pump model,
while elastic scattering of emitted light remains supported.

Cumulative accounts separately expose emitted/escaped/unresolved radiation,
passive absorption heat, prompt relaxation, nonradiative loss, radiative Stokes
heat and excitation storage. Caller allowances separately admit probability-row
roundoff, fixed-sample statistical diagnostics, cumulative unresolved joules and
local/global energy defects. No tail is relabeled heat or escape. Quasistatic
open index-one transport and fixed unsaturated pump absorption remain explicit.

Nine focused regression groups specify two-cell and self-return energy closure,
unpumped receiver heating during dark afterglow, split-clock agreement, actual
opaque-cell transport, tail admission, endpoint temperature rollback, original
report/binding integrity, pump-scattering refusal, work bounds and cancellation.
Pinned formatter/parser checks and independent source review passed. These new
native integration tests have not been executed in the shared editing environment.

## Retarded vacuum geometry and material proper-time histories

`spectral::doppler::retarded::UniformVacuumWorldline` admits one uniformly
translating material point with paired reference coordinate/proper clocks and
an explicit closed proper-time source history. `VacuumReception` describes an
independent receiver event and its instantaneous subluminal velocity in the
same world frame. The receiver need not move uniformly away from that event.

`observe` solves the unique past-light-cone quadratic analytically. Its scaled
discriminant and rationalized recession branch avoid squaring large ranges or
subtracting almost equal roots. Coordinate and material clocks use the
original `fs_math::dd::Dd` owner. The emitted position, path length and null
residual are recomputed from the published DD emission event; a delay lost at
an unresolved clock epoch is refused. Residual tolerance is an explicit
absolute distance in metres. It is a numerical admission criterion, not an
interval certificate for uncertain geometry. The same tolerance separately
bounds the loss from publishing the material proper-clock offset, converted to
equivalent coordinate seconds times c. An unresolved material epoch therefore
cannot silently select the wrong source history. A source-history miss, zero-flight
coincidence, invalid frame or unrepresentable arithmetic remains an error.

`RetardedVacuumRay` retains the physical propagation direction, emission and
reception events, source proper clock, flight time, range and the original
`VacuumDoppler`. Its coordinate-clock rate is dt_source/dt_receiver; its proper
clock rate is d(tau_source)/d(tau_receiver) and equals the observed/source
frequency ratio. SI velocity admission, c²-v² and the gamma-factored Lorentz
boost share DD arithmetic so near-light-speed source clocks and photon frames
use the same physical velocity. Zero-velocity and identical-frame paths avoid
unneeded frame transforms. Evaluation has fixed work and no allocation.

With `tracer`, `thermal_surface::doppler::retarded::RetardedThermalSurface`
calls the original source-history callback once at the admitted emission
proper time. Invalid geometry and expired history do not access that callback.
`ThermalEmissionState` carries the original diffuse thermal source and a normal
in its material rest frame. The aberrated source-frame photon direction checks
the emitting hemisphere and selects the original rough source's directional
emissivity. The returned cosine does not multiply radiance again. `RetardedThermalDetector` retains that ray and delegates directly
to the existing observed surface and `QuantumEfficiencySensor`, preserving
source identity, finite wavelength support, receiver-frame proposals and
receiver photon-energy conversion. No shifted source table is allocated.

G1/G3 regressions cover stationary and approaching/receding events, the
oblique retarded ray of a presently transverse source, independent clock
derivatives along a moving receiver, DD epoch limits, velocities one SI float
below c, source-history admission, proper-time thermal state selection,
aberrated emitting sides and the original QE rate path. This route supplies
neither scene visibility/surface intersection nor apparent projected areas,
finite-source geometry, smooth accelerated-worldline integration, moving dielectric transport,
gravitational shifts, polarization or a second debit to thermal storage.

`PiecewiseVacuumWorldline` adds a finite continuous history from authored
world-coordinate time/position knots. Each inertial interval derives its SI
velocity in DD, checks strict subluminality, and accumulates its own
`dt * sqrt(1-v²/c²)` proper-clock increment. The source's initial proper epoch
is kept separate until observation, and lost prefix/epoch increments remain
subject to the original clock-resolution tolerance. Preparation has an explicit
segment cap and fallible allocation; observation allocates nothing and locates
the causal segment with O(log N) monotone knot-residual queries before reusing
the original analytic null-ray and Lorentz owners. There is no time stepping,
constant-gamma substitution over the full history, or extrapolation.

At an exact interior knot, the outgoing interval supplies the right-hand
velocity and Doppler derivative; the final endpoint uses the incoming interval.
These velocity jumps are idealized impulses, not resolved acceleration. Smooth
motion requires caller-controlled trajectory refinement. `RetardedThermalSurface::observe_piecewise`
resolves the complete event before its single source-history query, then uses
the unchanged thermal spectrum, emitting hemisphere and receiver QE path.
G1/G3 controls cover changing gamma, reversal, knot semantics, finite support,
uniform subdivision, DD clock limits and proper-history-to-electron-rate
integration. The existing geometry, visibility and material no-claims remain.

Both worldline owners also provide `reception_at_proper_time`: a detector can
sample on its own clock without independently authoring coordinate time,
position and velocity. Uniform motion inverts the admitted `d(tau)=dt/gamma`;
piecewise motion locates the accumulated proper-time interval by binary search,
then inverts only that interval's rate. Earlier time-dilation increments are
retained. Receiver queries use the same outgoing/interior and incoming/final
knot convention, explicit closed history and allocation-free evaluation as
source queries. Repeated knots at an unresolved proper epoch cannot choose
between distinct velocities and are refused.

Receiver geometry now retains DD positions through the original light-cone
solver, preserving small detector displacements about large world origins.
Clock inversion is checked by reconstructing the proper interval from the
published coordinate event. Prefix and epoch losses remain subject to the
absolute distance tolerance after conversion by the receiver's inverse gamma.
An observation rechecks this receiver residual at its own tolerance, so an
event prepared with a loose allowance cannot bypass a stricter observation.
`VacuumReception` and `RetardedVacuumRay` retain the receiver proper clock and
its separate residual. Manually authored reception events retain `None` for
proper time; no detector epoch is invented.

`RetardedThermalSurface::observe_at_receiver_proper_time` and
`observe_piecewise_at_receiver_proper_time` feed these events through the
original retarded thermal and QE owners. Receiver admission precedes the
single source-history callback. The resulting collected-electron rate is per
receiver proper second; no second Doppler or coordinate-clock factor applies
when integrating exposure on that clock. Mixed uniform/piecewise endpoint
pairs can compose `reception_at_proper_time` with either existing observe
method. G1/G3 controls compare independent longitudinal/transverse laws and
finite differences of emitted proper time against the Lorentz frequency
ratio, and cover common motion, changing detector gamma, knot/history
semantics, uniform subdivision, DD receiver geometry, clock refusal and
proper-time thermal-to-electron-rate integration. These instantaneous APIs supply
detector events, not apparent-area geometry or smooth acceleration.

Both worldlines also provide `receive_emission` for the future reception of a
fixed DD emission event. The receiver uses the same scaled, rationalized null
quadratic as the past-source owner; piecewise motion first locates the causal
receiver segment by a monotone binary search. Public DD position pairs are
normalized before geometry. Published null geometry and receiver proper-clock
loss are re-admitted against the original distance tolerance, finite history
and outgoing/interior, incoming/final knot convention. No receiver motion is
extrapolated into the returned event.

`spectral::doppler::retarded::exposure` integrates nonnegative per-ray rates over
authored detector-proper-time panels. Uniform and piecewise endpoint owners
share four-point Gauss--Legendre temporal quadrature, DD accumulation, fixed
channel caps and an explicit ray-query budget. Both interval endpoints are
admitted before the first observable callback; the subluminal causal clock is
monotone, so a missing source-history tail cannot hide beyond the interior
quadrature nodes. Numerical node collapse, nonfinite/negative rates, original
geometry errors and cancellation refuse without publishing a partial integral.
Observable callbacks receive the same cancellation checkpoint and must read,
rather than advance, the source history. Quadrature is an estimate: authored
panel edges and refinement must resolve source discontinuities and velocity
changes. There is no certified error bound or claim about unsampled geometry.

`integrate_piecewise_with_knots` adds detector turns and received source turns
to those authored panels. Source knots are filtered against the two admitted
endpoint emission clocks before their future arrivals are queried, so unrelated
knots outside the exposure neither require receiver history nor consume ray
queries. Three ordered edge streams merge with exact duplicate removal and one
fallible allocation capped by the total ray budget. That budget includes two
planning endpoints, each relevant source-knot arrival, the two integration
endpoints and four quadrature rays per resulting panel. Planning, allocation,
history and complete panel-budget admission precede all observable callbacks;
cancellation remains available during admission and merging. The original
authored-only functions retain their behavior. Material-history discontinuities
and smooth-panel refinement still require caller control.

`thermal_surface::doppler::retarded::exposure::RetardedThermalExposure` connects
that owner to the original thermal history and `QuantumEfficiencySensor`.
At each temporal node it resolves the original emitting side and source state,
then checks the complete receiver sensor band against shifted source support.
It merges original QE knots and Doppler-shifted emissivity or rough-facet albedo
knots with a bounded streaming cursor, splitting at narrow spectral features
without a resampled table or heap allocation. Rough-source absorptance uses the
same retarded material-frame angle at every wavelength; the source law remains
nonlinear in facet albedo rather than interpolating transformed emissivities. Explicit wavelength subdivisions refine Gauss4
integration inside these smooth intervals without changing physical QE data.
Every spectral point consumes the caller's total Planck/QE-evaluation budget
and checks cancellation. Fixed-size DD channel sums retain expected collected
electrons and QE-weighted energy in joules; the latter is a response integral,
not an absorbed-heat debit. Photon energy uses receiver wavelength exactly once.
Its `piecewise_with_knots` entry uses the same motion-aligned plan before any
thermal-history query, sharing the original spectrum, QE and collection rules.

The supplied finite nonnegative etendue is constant in the receiver rest frame
and already includes projected area, solid angle and throughput. No extra
Doppler, clock, incidence-cosine or inverse-square multiplier is applied.
Changing optical collection geometry requires caller-controlled exposure
panels/etendue or a separate geometry owner. No visibility, stochastic charge,
shutter profile, electrical gain, saturation or material calibration is inferred.
Time and spectral quadratures remain estimates; sampled source-band admission
does not certify arbitrary unsampled source-state changes. G1/G3/G4 controls
cover analytic proper-clock integrals, receiver turns, refinement, DD epochs,
history/budget admission and cancellation, together with absolute SI thermal
energy/photon integrals, narrow shifted source/QE features, wavelength refinement
and dark-etendue/original-band refusal in the native thermal adapter tests.


## Coupled resonant photon recycling on the native excitation timeline

`LuminescentOpticalHeating<N,1>::with_resonant_recycling` explicitly changes
emission-band absorption from passive heat to fresh equal-energy excitation.
The sole emission wavelength must equal the stored-state wavelength and lie
inside the original source's excitation domain. Every positive-density material
cell receives a state and a directed emission row, including cells never reached
by the pump. No excitation band is broadened and no forbidden spectrum is clipped.

The original `fs-thermal` positive reaction/source propagator evolves all stores
simultaneously over the accepted thermal interval. It retains repeated self-return
and cross-cell reemission without delayed next-step injection, a generation cap,
a steady-state inverse or an artificial leakage floor. Quenching uses each native
start temperature. Prepared topology and interval maps are shared across observer
clones; identical decay parameters and duration reuse the cached map even when
pumping or stored excitation changes. The existing limits bound the complete
network to 64 states and 4,096 directed edges, with explicit numerical work limits.

Local accounts include incoming recycled energy, while global accounts exclude
internal emission/reabsorption from external energy sinks. Resolved recycling
never also enters the passive heat route. `recycled_excitation_j` exposes internal
circulation; the emitted/recycled totals can exceed external pump input. Dark
heat, escape, unresolved radiation, primary Stokes heat and current excitation
remain separate. Candidate local/global energy defects and cumulative unresolved
energy retain the original explicit allowances. Pulse boundaries, implicit retry
replacement, native temperature domains and complete rollback remain binding.

Six integration regression groups specify a one-cell effective-lifetime solution,
two-cell pulse/dark eigenmodes, asymmetric source/receiver orientation, a closed
unit-return system, supported-source/row/budget admission and transactional replay.
The closed fixture retains 32 J of excitation while recording 256 J of internal
emission/reabsorption. These are analytical regression expectations, not executed
results. Pinned formatter/parser checks and independent source review passed;
native execution was blocked by dependency-build storage exhaustion. General
wavelength-changing recycling, saturation and finite photon flight delays remain
outside this explicitly resonant, linear, quasistatic model.

## Polarized finite surface sources through native cameras and area sampling

`camera::SourcePolarization` retains constant or wavelength-dependent normalized
Stokes emission, with I exactly one and an explicitly authored geometry-local
positive-Q axis. The original `SurfaceEmitter` retains the absolute radiance
spectrum and sidedness. Source tables contain 2–1024 ordered positive vacuum
wavelengths, interpolate fractional Q/U/V by convex linear weights, and must
cover the complete radiance band. Queries never extrapolate or fit polarization
angles. Invalid states, axes and spectral support refuse.

Camera preparation binds source axes through the immutable native instance
placement once. At each actual emitter hit, the physical outward frame and the
original positive-coherency reciprocal contraction retain all Stokes components
and the A=diag(1,1,-1,1) convention. Exact source/analyzer extinction contributes
zero without terminating the original material continuation. The existing
prefix, source history, MIS weight, Beer loss and local phase-index radiance
conversion each apply once; logarithmic range recovery includes the source
polarization response. Sources without a polarization declaration preserve the
original arithmetic.

Original finite-area proposals launch that same source state and body axis into
the existing differentiated bundle/image owner. An area source's declared axis
also owns unpolarized spectral knots, so preparation and shooting use the same
frame throughout its spectrum. Rough and molecular receivers
therefore receive the original polarized illumination. Reverse emitter-hit
membership remains an unpolarized geometric probe, preserving proposal support
at a polarized dark port. No new propagator, random draws, optical solves or
per-query allocation is introduced; constant source evaluation is constant work
and table queries use bounded binary search.

Ordinary, analyzed, molecular, transient/history, shutter and spectral camera
observations share the same source contribution owner. Histories scale the full
state's intensity; they do not rotate its polarization in time. Focused native
regressions cover physical-cone admission, linear/circular/partial analyzer
responses, exact extinction, body rotation, spectral changes, independently
expanded total-internal-reflection phase, delayed source exposures and molecular
scattering. These are authored source laws and numerical model checks, not
measured emitter calibration, coherent interference between scene paths,
bulk birefringence or a resolved physical polarizer. The original opt-in
tracer/volume feature boundaries remain in force.

## Directional polarized coating heat reaches native camera observations

`camera::CoatingThermalEmitter` binds one original geometric side of an existing
smooth conductor to `CoatingThermalEmission`. Positive and negative faces may
retain independent thermal snapshots; duplicate same-side sources and competing
calibrated/diffuse thermal emitters refuse before publication. The source must
match the native conductor's complete bulk knots, source/status, provenance,
interpolation policy and exact zero/one absorbing or lossless film. Agreement
at a few sampled wavelengths is insufficient. The thermal owner's explicit
films retain their existing precedence over a stored substrate attachment.

Each actual hit queries the original positive layer-work/Planck owner at its
vacuum wavelength and outward geometric cosine, preserving every film's own
temperature and physical `I=L_s+L_p`, `Q=L_s-L_p`, `U=V=0`. The emission sample
also retains the original positive s/p radiances, before their sum/difference
rounding. Camera contraction combines the two pure-channel responses with
positive radiance fractions; a weak channel is not reconstructed by subtracting
nearly equal Stokes components. Exact extinction stays zero; an unrepresentable
positive response refuses rather than becoming false black. The physical s
axis is perpendicular to the outgoing ray/normal plane; the transported camera
frame resolves its exact normal-incidence degeneracy. The original reciprocal
contraction retains polarization through earlier Fresnel/scattering events
and camera analyzers. Existing source/radiance weights and logarithmic
range recovery apply once, with no averaged temperature, scalar emissivity,
tabulated lamp fit, extra cosine factor or thermal debit.

Binding allocates only a surface-indexed pair of optional face slots. Queries
reuse shared immutable optical data and the existing allocation-free coherent
source owner. Ordinary, analyzed, molecular, raster and transient paths share
the same source accumulator. A `NativeCoatingEmission` snapshot can supply its
original emitter directly; rebinding changes the observed frozen state without
advancing or mutating heat. No source implies no new thermal table or solve.

The observed mode retains the original local-index-one exterior approximation
and finite optical support. Rough/periodic surfaces, unsupported arbitrary
multilayer/transmitting native scene representations and nonvacuum emitting
hosts refuse. Directional thermal sources enter both camera hits and the native
finite-area proposals described below. Time-varying thermal states and bulk
birefringence remain separate work.

## Directional and spectral polarized environment observations

`EnvironmentRadiance::with_polarization` borrows one explicit source Stokes law;
`with_polarization_map` instead borrows one independently authored spectral law
per original row-major texel. Full visible RGB-reconstruction support is checked
before publication. Fractional Q/U/V retain the source owner's physical-cone
admission and convex spectral interpolation. The native map still determines
scalar intensity, spectral reconstruction, angular support and calibration;
color and brightness never imply polarization or measured spectral accuracy.

`EnvironmentEvaluation::texel_index` exposes the texel already selected by the
original direction/rotation lookup. Polarization uses that same result without
a second trigonometric lookup, resampling or interpolation across differently
oriented Stokes frames. Source axes live in the map's Cartesian coordinates;
the exact map rotation is prepared once and applied before transverse
projection. Source travel is opposite the escaped camera ray, and the original
reciprocal contraction retains A=diag(1,1,-1,1). Circular-only states have no
linear orientation and use the existing transverse camera axis at degeneracies;
a singular declared axis with a nonzero linear component refuses.

The shared camera endpoint observes these states through native dielectric,
conductor and molecular paths, camera analyzers, spectral/raster observations,
and the existing steady environment term of transient observations. Physical
camera-prefix intensity and the local phase-index conversion remain separate
and apply once. Exact source/analyzer extinction is zero emitted contribution,
while unresolved traversal and positive-response underflow remain refusals.
The original proven-escape rule applies; infinite true absorption is exactly
zero before source lookup. Unpolarized declarations preserve scalar arithmetic.

Map admission polls cancellation at bounded texel intervals and borrows storage;
queries allocate nothing, select one texel and perform at most one bounded
spectral lookup. Four focused numerical regressions execute the actual source
evaluation, original texel lookup and reciprocal polarization owner, alongside
the prior optical/thermal checks. They cover independent full-Stokes detector
rows, spectral/directional maps, map-axis rotation, circular-frame degeneracy,
identity, support and range behavior. Native scene regressions additionally
specify total-internal-reflection phase, actual escaped-path admission and
infinite absorption; those scene regressions have not yet been executed.

This is an explicitly authored, piecewise-constant angular source law. It adds
no atmospheric solution, polarized EXR decoding, time-varying polarization or
coherence between distinct paths. Its native direct-source proposal and
complementary camera weighting are specified below.

## Intensity-dependent two-photon beam absorption and native heat

`volumes::transport::deposition::nonlinear` binds explicit uniform or native-cell
bulk beta coefficients to the original absorbing grid. Beta is in m/W and
already accounts for the energy of both photons. Each source declares its
constant-beta wavelength band, irradiance ceiling and nonzero identity; prepared
models retain the intersection with the original linear spectrum. Numerical
identity binds every model and native cell index without copying density fields.
Zero-density cells remain vacuum. Bulk beta is applied at the authored material
state, with no inferred density, temperature or polarization scaling.

For a monochromatic beam of constant area A, each physical DDA crossing solves
`dI/ds = -alpha I - beta I²`, `I=P/A`, `alpha=rho*kappa_a` analytically.
The exact solution is `P_out=P*exp(-alpha*L)/(1+beta*(P/A)*L_eff)`, with
`L_eff=-expm1(-alpha*L)/alpha` and its continuous zero-alpha limit L. See
[Sheik-Bahae et al. (1990), equations (23)–(25)](https://doi.org/10.1109/3.53394).
Positive power partitions avoid subtracting nearly equal beam powers; a rare
logarithmic fallback retains finite watts when optical fractions underflow or
intermediate products overflow. Runs of beta-zero cells share the compensated
Beer-law depth accumulator; a nonlinear cell starts the next linear run at its
transmitted power, preserving the exact linear limit under cell refinement.
Only crossed cells are visited after preparation; output allocation and traversal
are explicitly bounded, cancellation is polled through the original DDA, and a
late failure publishes no partial report.

`BeamPowerDeposition<1>` carries the absorbed watts into the existing
`OpticalHeating` and enthalpy timeline. This requires the explicitly declared
passive, promptly thermalizing channel: the model does not infer retained
carrier/excitation energy, fluorescence, bleaching, saturation, nonlinear
refraction, diffraction, reflection or cross-frequency absorption. Input power
is already inside the medium and uniform across the given area. Pulses require
their instantaneous power; replacing a pulse by mean power or adding independently
deposited overlapping beams is generally invalid for this nonlinear law.

G0/G1/G3/G4 regressions cover an independent Runge–Kutta attenuation solution,
linear and pure-quadratic limits, weak and opaque cells, extreme representable
power, subdivision, area scaling, non-superposition, heterogeneous source maps,
source support, budgets and cancellation. A native enthalpy regression checks
local temperatures and incident/transmitted/thermal energy closure. These are
numerical checks of the declared constitutive law, not measured beta calibration
or evidence that other high-intensity processes can be neglected.

## Directional polarized coating heat through finite-area light proposals

`CurvedAreaEmitter` now admits the original `CoatingThermalEmitter` snapshot
on its selected positive or negative source face. Camera-derived preparation
selects an explicit face when both sides emit; complementary-hit matching
retains exact source identity, side, geometry and body placement. Frozen source
refresh retains the prepared area support but rechecks the original native
smooth-conductor optical law. The existing calibrated and diffuse thermal
source paths retain their radiance and sampling arithmetic.

Shooting, image discovery and reverse MIS membership remain geometric and
unpolarized. Each converged image evaluates the original coherent layer-work
and Planck owner once at the actual emitted direction and vacuum wavelength.
The retained positive s/p radiances stay separate through native transport and
receiver response, preserving weak power that a rounded I/Q difference loses.
The physical outgoing-ray/normal plane determines the source frame; a caller's
shooting chart does not rotate the thermal polarization. No initial-angle
evaluation or fixed lamp spectrum substitutes for the directional source law.

Unequal nonzero s/p channels each use one fixed-ray native replay of the
original image's branch policy. The original connection owner must accept its
first trace; its existing replay validator checks the same branch signature
and receiver residual. Exact-zero native interface power continues with the
existing negative-infinity log-power marker and contributes zero; generic
extinction/errors are not treated as successful dark paths. Representable
near-Brewster residues remain positive. Equal channels reuse the original
unpolarized connection, zero channels skip work, and all replay traces, field
evaluations and geometry queries count against the original work limits.

Standalone irradiance, surface NEE, molecular NEE and complementary camera
hits share this source law. MIS densities and image counts belong to the
geometric image, not to its two polarization replays. A molecular proposal
reuses its original scattering-survival realization for both deterministic
channels, without additional stochastic free-flight draws. This extension
retains the thermal source's original finite spectral support, local-index-one
emission domain and exactly matched native smooth-conductor/film restrictions.
Focused native regressions cover a nonisothermal absorbing film, converged-angle
irradiance, replay budgets, weak near-Brewster filtering, rotated emission
frames, face/snapshot refusal, complementary hit/NEE MIS and clear molecular
parity. They are numerical model regressions, not measured thermal calibration.

Validation in this continuation is limited to pinned-nightly parsing of the
changed Rust leaves, formatting of the new leaves, and an executed actual-source
Fresnel/Jones probe retaining the finite near-Brewster pure-p channel. The three
native scene regression groups are authored but unrun under the current build
and storage limits; the numerical probe does not establish scene replay or NEE
integration.

### Temperature-dependent two-photon heating and shielding

`TemperatureOpticalHeating<1>::try_new_two_photon` binds independently sourced
bulk beta(T) tables to the original non-vacuum optical cells. Each table carries
a constant-beta wavelength band, an incident-irradiance ceiling and ordered
kelvin knots. Linear interpolation does not extrapolate; bulk beta remains in
m/W without density scaling. Linear and nonlinear source assignments remain
independent, while both use the original receiving thermal cells.

The existing cached DDA, linear temperature tables, positive analytic
two-photon attenuation kernel, native thermal observer and complete-trajectory
refinement remain the owners. Each native step refreshes alpha(T) and beta(T),
so upstream heating changes downstream nonlinear shielding. Applied-source and
endpoint checks enforce the local irradiance domain, including final endpoints
with no subsequent source query. The endpoint transport pass allocates no source
map and consumes the original source-evaluation allowance. Linear and nonlinear
tables and knots share the declared preparation budgets.

The monochromatic constructor requires a constant uniform beam area and power
already inside the medium. Prompt thermalization is explicit. It does not infer
cross-frequency absorption, carrier storage, saturation, nonlinear refraction
or changing density/geometry. Coefficients remain frozen over one native
explicit step; exact spatial attenuation does not remove temporal error.

Focused regressions cover an independent nonlinear heating solution and time
refinement, changing shielding, weak/dark/opaque limits, area scaling, native
latent/event/restart behavior, and transactional source-domain/work refusal.
These are numerical constitutive checks, not measured-material qualification
or a certified global time-error bound.

## Polarized environment next-event sampling on native clear paths

`camera::environment::sampling::CurvedEnvironmentCameraScene` borrows the
original camera and its calibrated environment. Each continuous native surface
receiver takes one original environment-map direction proposal, including its
existing luminance-times-solid-angle distribution and conditional PDF. The
sampler retains the exact selected row-major texel; a rounded pole or seam is
not re-inverted into another texel's intensity or polarization. Original
spectral Stokes fractions, map placement and physical incoming-light frames
remain intact. Polarized power does not replace or modify the scalar proposal.

Diffuse, coated-diffuse, rough dielectric and rough/brushed conductor receiving
responses reuse the original native reciprocal BSDF evaluator. The actual
incoming direction selects the original adjacent source-side medium. Passive
camera analyzers act at launch, before any interaction. Source radiance, the
physical receiving cosine, the native solid-angle BSDF, camera prefix, original
local-index conversion and proposal PDF are each applied once. Positive
products recover representable values through logarithms when an intermediate
overflows or underflows; exact optical zeros stay zero, while an unresolved
positive range failure refuses the complete sample.

A direct proposal represents only a homogeneous, nonabsorbing source segment
with no intervening surface. The original first-boundary event owner propagates
the actual receiving contact's departure witness without a ray-origin offset
or removal of that object. Any subsequent surface blocks this direct strategy.
At a finite distance limit, the original outward-bounded geometry owner must
also prove that the remaining infinite halfline escapes. A budget ending inside
possible geometry is unresolved visibility, not an assumed clear or dark ray.
Source samples and cumulative geometry work have explicit caps; optical field
evaluations and shadow-query counts remain observable. No per-ray allocation
or new geometric representation is introduced.

The original camera continuation and its Philox domain are unchanged; source
directions use a disjoint counter domain or explicit independent samples.
Complementary power-heuristic weights compare the same two solid-angle
proposals. Only a last continuous departure followed by a supported clear
escape receives the environment endpoint weight. An intervening delta event
clears that competitor, so unrepresented refractive or reflected chains retain
their camera contribution. Encountered surface and thermal emission retain
their original accumulator and are never multiplied by the environment weight.
A homogeneous absorbing infinite tail has the original exact-zero environment
response; curved-index source segments retain camera transport without an
invented angular Jacobian.

This surface-only view has no participating-medium law. With `volumes`, the
participating environment consumer below composes the same source sampler and
visibility owner with the original finite-grid scattering support. Arbitrary
active fields and unresolved infinite tails remain outside that proposal.
Focused native regressions specify original-texel spectral/frame consistency,
analytic diffuse MIS with independent surface emission, rough-metal linear and
circular polarization, native blockers, unresolved distance limits and delta
chains. The actual source sampler, source-frame and MIS/range arithmetic passed
a bounded numerical execution; the full native scene regressions were authored
but not executed in the constrained shared build environment.

## Source-bound full-Stokes Mie scattering and particle transport

`MieSolution::amplitudes` retains the original complex S1/S2 series in the
documented exp(-i omega t), right-handed s/p/travel convention. Its differential
response stores both nonnegative diagonal Jones powers and the complex cross
product before conversion to I/Q. `SphericalPhaseValue` and
`SphericalPhaseTable` admit positive semidefinite coherency; convex angular
and spectral interpolation and independent populations mix coherency, never
complex particle amplitudes. Population weights use the same original
scattering cross sections and number distribution as the scalar owner.

`polarized_phase_table` refines all four independent real response components
against the same direct multipoles, using bounded angular knots, cancellation,
normalization and first-moment checks. Its scalar CDF and M00 share the exact
same admitted table. The reported discrepancy is sampled numerical evidence;
it does not certify unprobed angular or spectral resonances. Existing scalar
preparation remains intensity-only and preserves its numerical arithmetic.

`MieOpticalSpectrum::from_population_source_polarized` prepares the original
dispersive population at each declared knot and diagnostic midpoint.
`MaterialParticleScattering::source_spectrum_polarized` uses the SAME native
material query and sphere solver, retaining accepted state, original n/k laws,
temperature, particle density, host and radius/number distribution. Absorption
and total scattering remain the original scalar coefficients. Every declared
polarized knot requires positive scattering; source support is never
extrapolated or replaced by an invented isotropic Mueller matrix.

The original `PhaseFunction` dispatcher retains this full response, including
spectral tables and independent species mixed with their original
mass-loading-times-scattering weights. Scalar evaluation, CDF sampling,
contraction and enclosure consumers retain the original M00 law. Full-Stokes
collision sampling uses the actual wavelength and M/PDF response in the
physical scattering frame; Rayleigh/depolarized ordinary arithmetic is
preserved. Positive s/p channels and a bounded binary-exponent fallback retain
weak responses that rounded M00-M01 or an intermediate quotient could erase.
Exact optical extinction stays zero; unresolved positive range failures refuse.

The borrowed `RayleighGrid` adapter and density-bound curved scattering media
admit these retained particle laws without changing native free-flight,
absorption, density/majorant or geometry ownership. Camera point and area-light
connections use the common reciprocal full-Stokes contraction, including
complex phase and circular polarization. Area-light continuation anchors
borrow cheap shared phase tables and evaluate the original wavelength-specific
PDF. Sharp particle lobes replay the original sampled direction for the PDF
consistency guard; independent geometric checks and each image's actual-arrival
MIS density remain separate. No Mie recurrence, radius-bin solve or material
query runs at a transport vertex.

Scalar-only HG/isotropic/intensity tables cannot determine polarization and
refuse full-Stokes transport. The original independent nonmagnetic homogeneous
sphere domain and real nonabsorbing host restrictions remain. Nonspherical or
aligned nonspherical particles, dependent scattering, coherent interference
between distinct particles and nonreciprocal magneto-optics are not supplied.
Active-scattering environment NEE still requires an original infinite-tail or
finite-support law; this extension does not infer one.

Executed focused checks cover direct Jones/Mueller agreement, the Rayleigh
limit, incoherent population weights, all-component angular refinement,
original dispersive spectral preparation, every input polarization's sampled
power, frame covariance, reciprocal connections, mixed-species wavelength
PDFs, exact zeros and representable weak channels. The narrow-peak camera
anchor arithmetic was executed separately with actual native phase owners.
Native material-state and complete camera/grid scene integration are not
claimed by those bounded executions; the infrared material regression is
authored for the full native build.

## Polarized environment connections through native finite density grids

With `volumes`,
`camera::environment::sampling::participating::CurvedParticipatingEnvironmentCameraScene`
borrows the already identity-checked `CurvedRayleighCameraScene`. Continuous
surface vertices and accepted physical volume collisions use the original
environment direction/texel sampler. Rayleigh, depolarized Rayleigh, retained
spherical Mie tables, spectral tables, and admitted physical mixtures reuse
`PhaseFunction::polarized_connection_factor` and its wavelength-specific
`eval_at` PDF. The physical source frame and full Stokes state enter the
original reciprocal phase response. No scalar particle law is promoted to
Mueller polarization, and no additional cosine is applied at a volume vertex.

The original `try_new_volume_grid` binding retains its exact density grid and
finite support. In a homogeneous host, the native cell traversal supplies the
forward intersection interval and integrates the piecewise-constant mass
column. An unobstructed source connection uses
`log(T_s) = -kappa_s(lambda) * sum_cells(rho_cell * length_cell)`.
Its finite extent is a property of the original grid constructor, not an
inference from one sampled zero or from a caller's distance limit. Erasing the
grid behind an arbitrary density-field binding does not retain this admission.
Absent or globally zero scattering needs no cell visits; a known positive
uniform infinite tail has exact zero unscattered source radiance. GRIN hosts
and other active fields remain unsupported by this straight source proposal.

Camera majorant flights use the same retained grid interval: they skip only
the proven empty prefix and stop at the native far face when the sampled
majorant flight leaves support. Within support, the original exponential,
independent acceptance draw, local density, polarized phase sampler and
surface-versus-collision arbitration remain authoritative. This avoids an
unbounded exterior null-collision tail without changing an active arbitrary
field or treating a rejected candidate as a physical scattering event.
Scattering survival is already represented by these camera flights; it is
never multiplied into the camera endpoint a second time.

The source leg retains the original native first-boundary visibility query.
Surface receivers carry their real departure witness; volume receivers start
at their actual collision point without a fictitious surface or offset. Any
intervening interface blocks this direct technique. Only homogeneous geometric source
segments with a proven clear infinite geometric tail participate. The host
remains nonabsorbing; explicitly bound native finite absorption grids contribute
their original true-absorption column in index-one hosts as specified below.
Refractive images, absorbing infinite hosts and curved angular maps retain their
existing camera paths.

Both techniques use a deterministic complementary angular partition from the
same environment and BSDF/phase solid-angle PDFs. Source survival multiplies
the direct radiance only; neither a stochastic survival estimate nor a new
collision rate enters the MIS weight. This does not claim the full marginal
path-density power heuristic or variance optimality. Physical volume vertices
replace the continuous anchor, null tracking candidates preserve it, and delta
surface events clear it. Independent environment streams use native surface
ordinals or accepted volume ordinals; original surface and flight Philox
domains are unchanged. Encountered surface and thermal emission remain in
their original accumulator and are never reweighted by the environment PDF.

Source-grid visits, environment samples and geometric work have cumulative
explicit caps. Finite log survival is combined with radiometric factors before
a potentially underflowing exponential can erase a representable bright
source. Exact optical zeros and the declared infinite-tail zero stay zero;
unresolved positive range, support, work or cancellation failures refuse the
complete sample. There is no per-ray heap allocation or copied grid/phase/CDF.

Five focused consumer regression groups specify exact zero-law surface
parity, heterogeneous source columns and null-collision MIS retention, native
Rayleigh and spectral-Mie volume connections, refusal to infer arbitrary-field
support, and positive radiometric recovery from a very small survival factor.
The native integration regressions and new Rust parsing/format checks were
not executed in this continuation because the execution environment became
unavailable. These source-reviewed tests are not a claim of end-to-end build,
render convergence, or experimental model validation.

Four additional native-cell regression groups specify directional face ownership, layered and oblique density columns, finite support, zero-work misses, crossing budgets, cancellation and numerical-range refusal. A fifth support-distance regression specifies preserved vacuum prefixes and far-face truncation. These five core regression groups are also unexecuted after the workspace interruption; finite-column overflow refuses rather than being relabeled as an infinite uniform tail.

### Pulsed temperature-dependent optical absorption

`TemperatureOpticalHeating::with_power_schedule` binds the existing immutable
`SpectralPowerSchedule` before the observer's first advance. Original beam powers
bound each spectral lane. Native intervals stop at affine-segment boundaries,
jumps and dark-gap boundaries; the original accepted duration owns every energy
credit. A pending interval publishes powers, cell heat and cumulative accounts
only with the matching accepted thermal step.

Frozen linear transport uses the original kernel at the exact interval mean.
Monochromatic two-photon transport integrates the entire affine waveform
through the cached heterogeneous path: each prefix is a passive rational map,
and positive analytic first/product moments yield local heat without
substituting mean intensity into the nonlinear law. Peak transport checks
every interior ramp irradiance; actual left/right endpoint powers check changing
shielding and jumps, including final endpoints. Repeated evaluations retain
O(path cells × spectral lanes) work and source-evaluation/cancellation budgets.

Focused numerical regressions cover independent time quadrature, nearly equal
shape parameters, weak-loss and opaque limits, heterogeneous cascades and
interval subdivision. Public regressions cover spectral fluence, native
events/restart, nonlinear ramp energy, local shielding, temperature refinement
and transactional endpoint refusal. Exact waveform integration does not remove
temperature feedback error: the same native explicit state and
complete-trajectory refinement remain authoritative. Geometry, density and beam
area stay prescribed; carrier storage, material saturation, nonlinear refraction
and ultrafast relaxation remain unmodeled.

### Common-frame Jones pupil diffraction

`PhysicalCamera::prepare_polarized_coherent_diffraction` accepts the existing `CoherentJonesPupilWave`, `CoherentJonesPupilIllumination` and
`JonesPupilCoherenceFactor` inputs owned by Jones spectral exposure. These
same authored two-component waves propagate through the focused or
declared-plane defocused circular, annular or polygon pupil transfer. Every wave uses one
complex scalar transfer for both components in the same declared transverse
pupil frame and with the same pupil-centre phase reference. The axes are
camera right/up for common incoming propagation along -camera.forward,
forming the same right-handed convention used by native Stokes transport.
The deterministic mean and each authored covariance factor interfere internally before their
coherencies are formed. Independent factors add coherencies without cross
terms. Equal wavelength alone never establishes mutual coherence; source
weights remain square-root amplitudes and no isolated-power renormalization
is applied.

Point queries return the deterministic Jones field and separate mean,
fluctuation and total coherencies. Each coherency retains its two diagonal
powers and complex cross moment independently of its native Stokes result:
`I=Cxx+Cyy, Q=Cxx-Cyy, U=2 Re(Cxy), V=2 Im(Cxy)`, with
`Cxy=Ex conj(Ey)`. DD field and moment accumulation precedes rounding, so a
weak diagonal remains directly available when rounded Stokes I and Q agree.
Exact coherent cancellation remains dark. A product of nonzero complex
wave/transfer values that underflows completely, and a represented nonzero
field whose axis power underflows, are refused rather than reported as exact
extinction.

Finite physical-pixel queries integrate the four coherency moments with the
existing positive area quadrature, after interference at each sample. They
retain spatial depolarization instead of assigning one Jones field to an
intensity-integrating detector. The existing paraxial, material/wavelength,
pupil, source/factor and exact pixel-work budgets remain authoritative.
Each authored wave costs one pupil evaluation per quadrature sample shared
by both Jones components. Immutable preparations allocate nothing during
queries and remain reusable after cancellation or refusal. The generic
moment integrator preserves the original scalar/exposure arithmetic order.

This is monochromatic paraxial propagation in one common pupil frame, using
the original scalar transfer. It does not implement vector high-NA
diffraction, a polarization-dependent pupil, inferred scene coherence, or
cross-frequency Stokes exposure. The existing Jones spectral-exposure query
continues to own cross-frequency detector-power integration with these shared
inputs. Quadrature refinement controls
resolution; it is not a certified model-error bound. No relative-error
guarantee is made at a dark fringe.

Six native regression functions cover scalar parity across the original pupils,
coherent versus independent Jones illumination, polarization and covariance
basis changes, finite-pixel depolarization, cancellation and weak-axis range,
and admission/work/cancellation bounds. They were source-reviewed but could
not be executed because the execution environment was unavailable.

### Full spectral coherency over the original finite shutter

`PreparedCoherentJonesDiffractionExposure::observe_coherency` returns the
shutter-averaged full coherency and native Stokes state in the original camera
right/up frame, for incoming propagation along -forward. It uses the existing
Jones source, covariance-column identities, equal-frequency field grouping,
wavelength-specific pupil caches, chromatic sensor phases and analytic
`ProperHarmonicProduct` kernels. The matrix is the average outer product of
the analytic field; it is not the outer product of a time-averaged field.
Deterministic multicolor illumination can therefore depolarize during a finite
shutter. Independent source columns add coherencies without cross terms.

The four DD moments are Cxx, Cyy, Re(Cxy) and Im(Cxy), with Cxy=Ex conj(Ey).
Their direct temporal contraction preserves native I=Cxx+Cyy, Q=Cxx-Cyy,
U=2 Re(Cxy), V=2 Im(Cxy). No analyzer-intensity differencing is used. Axis
moments and the complex cross term remain available independently of rounded
Stokes I/Q. The result retains separate deterministic-mean, covariance and
total contributions. Diagonal range, native Stokes admission and direct
weak-channel correlation checks refuse invalid contractions without clipping.

`integrate_pixel_coherency` integrates these same moments with the original
positive spatial rule, shortest-wavelength pupil bandwidth and mixed-frequency
sensor-chirp bound. It returns spatially integrated average coherency in watts
and shutter-integrated coherency in joules; all signed Stokes components share
those units. This retains temporal and spatial depolarization. Each unmasked
vector wave evaluates its physical pupil once per spatial sample; optional
patterned Jones pupils charge every retained Fourier mode in that same count.

For n distinct frequencies and each mean/covariance bank, full coherency costs
four diagonal real moments plus six per unequal-frequency pair: n*(3*n+1).
The original preparation's pair limit preflights the full query, while existing
pixel limits bound total contractions, samples and pupil evaluations. A budget
that admits detector power may refuse the more expensive full matrix. Original
detector-power arithmetic and grouping are preserved. Shared source products
now refuse a wholly underflowed product of nonzero complex wave/transfer values
instead of reporting false exact extinction. Queries use bounded scratch,
allocate nothing, and remain reusable after refusal or cancellation.

The original uniform-shutter, cycle-averaged analytic-field detector and scalar
paraxial pupil assumptions remain. There is no optical-carrier time raster,
subcycle real-field detector, inferred coherence or vector high-NA diffraction.
Passive patterned Jones pupils use the shared field owner described above.
Spatial resolution is not an error
certificate, and DD moments do not guarantee relative accuracy at an arbitrary
interference null. Point results are tangent-coordinate power densities; this
does not change physical scene radiance or imply a sensor QE/electronics model.

Six native regression groups cover analytic polarization beats, analyzer agreement,
shutter depolarization, covariance identity, frame rotation, finite-pixel chromatic
integration and energy, weak-axis retention, independent time-quadrature
convergence, and work/cancellation/range admission. All six groups passed in the native `coherent_diffraction`
target with `photometry`; all 31 tests in that target also passed, including
the prior monochromatic coherency and original detector/pixel regressions.

### Direct-beam absorption through native latent enthalpy states

`TemperatureOpticalHeating::try_new_enthalpy` consumes the existing
`EnthalpyAbsorptionTable` and explicit `EnthalpyOpticalCell` bindings. The
original enthalpy optical compiler and native binding checks remain the owners
of spectral interpolation, caloric identity, material identity, storage basis
and h support. An independent kelvin domain also remains binding. Zero and
negative h retain the original source datum; no phase fraction or optical
mixing rule is inferred.

The same cached DDA, positive attenuation kernels, native accepted-time
accounting, power schedules and complete-trajectory refinement serve both
temperature and enthalpy coordinates. Native h enters the optical snapshot
identity and source queries, including dark intervals and final/event endpoint
checks. The solver uses precompiled spectral lanes with bounded per-cell state-table
lookups, without new geometry traversal or spectral source searches during an advance.

`try_new_enthalpy_two_photon` composes authored linear kappa(h) with the
existing monochromatic bulk beta(T). Both coefficients read the same native
start state, and the existing pulse owner integrates each affine waveform.
Beta does not acquire an inferred phase dependence. Combined source table/knot
limits, local irradiance support, cancellation and transactional publication
remain binding.

Focused regressions cover native latent storage, caloric datum invariance,
original-source refusals, transparent-cell endpoints and temperature-path
parity. Public integration regressions additionally cover analytic latent
heating, downstream shielding, nonlinear pulses, events/restart and native
time refinement. Geometry, density and beam area remain prescribed; hysteresis,
phase mixing, carrier storage and material saturation are not inferred.

## Native finite absorption and constituent LTE on polarized participating paths

With `volumes`, `CurvedRayleighScene::with_absorbing_grids` borrows one optional
original `AbsorbingGrid` per native scene medium. The existing grid retains its
finite support, density cells, absorption spectrum/palette, constituent source
law, temperature state and content identity. No density, spectrum or source
table is copied, fitted or inferred. A nonempty binding table must match the
original medium table. Bound hosts must be homogeneous and, at each actual
transported wavelength, index-one and clear; unsupported source wavelengths
and GRIN carriers refuse. Scalar absorption preserves polarization and does
not imply dichroism, birefringence or aligned-particle extinction.

The original curved walk supplies each segment's actual start ray and accepted
physical arc length. In a bound homogeneous host, the original finite-cell
traversal integrates true absorption and optional constituent LTE on that
straight segment. Every Stokes component receives the same scalar attenuation;
the native endpoint, frame, medium, scattering response and proposal are
retained. An original geometric escape proof also integrates the remaining
finite grid column, including a grid beyond the finite geometry integration
distance. The clear unbounded distance after its support is not added to the
reported finite geometric/optical lengths. Positive finite optical-depth
overflow and unrepresentable nonzero attenuated components refuse; they never
become an exact extinction marker.

Reciprocal camera paths accumulate the grid's original unpolarized LTE source
against the actual arriving full Mueller intensity before applying that
segment's absorption. The existing reduced-volume source accumulator applies
the camera-index conversion once. Original scattering flights already supply
their own survival probability, including null candidates and native finite
support exit; no scattering Beer factor or extra albedo is applied to camera
segments or to their local thermal source. Forward source-power paths evaluate
no Planck source. A medium bound to a finite absorption grid cannot also bind
an isothermal or spatial camera temperature emitter, which would otherwise
reinterpret the same absorption under a second source law.

Environment source connections add the exact finite absorption depth to the
original deterministic scattering log survival. The shared angular MIS
partition is unchanged, and camera endpoint weights receive no second source
extinction factor. Explicit environment NEE contains environment radiance
only; local thermal emission remains in the original camera accumulator.
Native point/area source connections reuse their same-branch bundle replay,
whose observer now retains the exact segment-start ray. Scalar finite
absorption remains separate from geometry, interface coherency and receiving
Jacobian; `log_source_survival` combines it with the original scattering
estimate for all participating point/area consumers. The existing
`log_scattering_survival` accessor retains its original meaning, and a separate
absorption accessor retains finite deep columns even when their isolated
exponential underflows. Only the original ratio-tracking exact zero or the
declared positive uniform infinite scattering tail supplies negative infinity.

The finite-grid adapter exposes actual absorption-cell visits and a cumulative
per-path/per-source-replay cap, bounded again by each grid's original segment
cap. Environment queries count both scattering and absorption cell traversals
against their existing cumulative source-cell allowance. Original cancellation,
geometry, integration, free-flight and collision budgets remain binding. No
per-ray heap workspace or additional random stream is introduced.

Seven native regression groups cover full-Stokes finite segments and escape
tails, temperature ordering and analyzer throughput across a null flight,
the original real Rayleigh Mueller collision, carrier/domain/work/source
conflicts, deterministic source replay with separate absorption/scattering,
retained finite depth 800, and polarized environment NEE with complementary
camera MIS. All seven passed with native `tracer,volumes` execution, along
with the existing equilibrium control selected by the same filter (eight total).
The native cell-transport group passed all 14 cases, including the four new
observed-column controls; all six constituent thermal-source cases and all 15
source-connection cases also passed. The source-connection group includes the
two new replay cases. Four additional public-API end-to-end cases passed against
the actual compiled renderer with real dependencies and an independent SI Planck
oracle; no production implementation was copied or stubbed.

The shared observed grid integral retains compensated finite optical depth and
actual native cell work without changing the original public transport arithmetic.
It reuses the original palette and thermal source cache. A logarithmic fallback
combines a bright local source with its foreground survival when that intermediate
weight underflows; a nonzero unrepresentable contribution refuses explicitly.
Four focused observed-column groups cover ordered thermal radiance, opaque finite
log columns, source-range recovery, native/cumulative cell budgets and cancellation.
Original constituent-source coverage also checks the observed transfer against the
same retained source palette. Native execution uses the required `tracer,volumes`
features; unrelated photometry-only unit-test compilation failures are separate.


## Rough diffuse thermal reflection and emission

`RoughDiffuseThermalSurface` binds the original finite-band spectral **facet
single-scattering albedo**, `RoughDiffuseSurface`, positive source temperature,
and nonzero physical-source identity. The facet spectrum is an explicit input;
a measured hemispherical reflectance or emissivity is not fitted or reinterpreted
as facet data. Both reciprocal reflection and directional thermal emission use
the same EON roughness and multiple-return closure. Source identity includes the
original knots, roughness, temperature, source hash and constitutive-law domain.

Kirchhoff emissivity is the original EON directional absorptance, evaluated at
the actual outward cosine. Spectral radiance is that emissivity times the shared
vacuum SI Planck source in W/(m² sr nm). The stable absorptance factors `1-rho`
and retains representable weak absorption near a white reflector. Exact black
facets emit a blackbody, exact white facets emit zero, and zero roughness retains
the original spectral Lambertian radiance arithmetic. Hemispherical mean
emissivity and one-face spectral exitance use the same analytic EON angular
integral. They do not supply an unmeasured out-of-band tail or advance cooling.
The admitted source remains isothermal, opaque, depolarizing and two-sided;
each face emits into its own hemisphere in an index-one medium.

`Material::RoughDiffuseThermal` is the sole thermal emission owner on its
primitive. Scalar camera hits, thermal next-event proposals and volume-phase
source connections evaluate directional emission at the source's actual normal
and path direction. Reflected continuation and direct receiver evaluation use
the same original EON BRDF and explicit facet spectrum as nonthermal rough
surfaces. Prepared visible light selection uses hemispherical mean radiance only
as a proposal weight; it never substitutes that mean for transported radiance.
The wavelength support remains the declared source support, and visible proposal
preparation refuses spectra that do not cover its complete visible band.

`NativeOpticalScene` retains the original thermal owner and binds its exact
facet spectrum and roughness to a native `RoughDiffuse` interaction. Native
emitter admission rejects a substituted facet curve, substituted roughness or
Lambertian reinterpretation. Native camera hits and converged area-source images
share the same geometric directional-radiance method. Rough area emission is
evaluated on the original converged source ray; its unpolarized source reuses the
original transported bundle without an extra image search or polarization replay.
Area geometry, projected Jacobians, proposal densities and complementary MIS
remain unchanged. Smooth, black and nonemitting sources retain the earlier
uniform-source fast path. Existing refractive-host refusals remain binding.

Direction-free borrowed-source radiance and reflectance queries refuse rough
sources explicitly. Direction-free Doppler and Lambertian visible-observable
entry points also refuse them rather than silently using an angular average.
The explicit directional Doppler view, retarded thermal observations and finite
receiver-proper-time exposures retain the original rough-source angle and
shifted facet spectrum. Focused G1/G3 regressions cover rest-frame aberration,
D^5 radiance, receiver-photon conversion, exact zero-shift source parity, narrow
shifted facet features, source-band admission and spectral-work refusal. Scalar thermal BDPT light subpaths and calibrated-fixture
thermal emission remain outside their existing supported contracts. This does
not claim measured rough-material authority, thermal polarization, resolved
surface microgeometry, nonisothermal rough cavities or a bolometric heat debit.

Validation probes independently integrate the original BRDF over the projected
hemisphere to verify that reflected blackbody illumination plus emitted radiance
equals the same Planck source in each tested direction. A separate angular
emission quadrature verifies spectral exitance. Controls also cover exact
smooth/black/white limits, positive near-white absorption, original-source
identity, finite-band and angle refusals, exact native facet/roughness admission,
and both geometric sides.

## Native sourced room-temperature water absorption

With `material-state`, `BeerLambertAbsorption::hale_querry_water_25c()` consumes the complete, shared 169-knot Hale–Querry (1973) imaginary-index source from `fs-material`. It retains 200–200000 vacuum-nm absorption support, source temperature 298.15 K, linear-k interpolation before SI intensity conversion, existing zero/thin/deep Beer semantics, and source-sensitive material identity. No per-material table enlargement or allocation is introduced. This historical, fixed-temperature compilation supplies no pressure, temperature, scattering or mixture dependence; strongly absorbing infrared interfaces need their separate complex-index owner.

`material_state::water_refraction::hale_querry_water_25c()` binds that absorption to original IF97 Region 1 density at explicitly chosen 101325 Pa and 298.15 K and original IAPWS phase/group refraction. The combined medium admits only their common 200–1100 nm band. Native tracer and Euler-disc material identities retain the full source identity and preserve existing tags; the Euler-disc exhaustive consumer also admits original imaginary-index tables and CHARMS glass provenance.

Four focused native regressions cover source Beer semantics, reference-state/source-band binding, obliquely refracted physical-path absorption with phase/group length, and source identity. They are authored and formatted; their execution remains unverified because the existing original-renderer build was killed in upstream `asupersync` before `fs-render`. The source material's three independent-control tests and strict Clippy passed separately.

## Native mechanically induced Corning 7980 polarization

With `material-state`, `PolarizedBeam::advance_corning7980_stress_optic` consumes the original uniform Cauchy-stress tensor for a physical segment, derives its transverse differential index in the current right-handed beam frame, and applies its relative polarization phase. It reuses the original Beer advance, medium, metric path and Stokes validation. One half-angle sine/cosine evaluation retains weak-retardance mixing without matrix allocation. Successive stressed pieces compose in traversal order, preserving noncommuting principal axes.

The method requires the native homogeneous CHARMS 7980 provenance and one measured stress-source line. Their present intersection is 633 nm; retained 193/248 nm stress measurements do not extend the scalar phase law. Positive stress-optic coefficient and tension-positive stress select the higher-index slow axis; horizontal tension sends +45-degree input toward positive V under `exp(-i*omega*t)` and `V=2 Im(E_h conj(E_v))`. The source calibrates differential stress retardation in laboratory uniaxial tests, not the scalar CHARMS temperature response or a general multi-axial/temperature envelope. Common stress-induced phase/group shifts, ray splitting, walk-off and new interface physics are not supplied.

Four native regressions are authored and formatted for independent signed Jones/Stokes controls, crossed-polarizer signals with original Beer loss, rotated/reverse frames, segment composition/order and source/material refusals. They remain unexecuted because the established full-renderer dependency build is blocked. The material source's three original assertion bodies executed separately; this is not a full native-transport gate.


## Composition-resolved atmospheric refraction

With `material-state`, `MaterialAirCompositionRefractiveField` binds the original pressure/temperature reconstruction and an explicit water-mole-fraction/dry-air-CO2 reconstruction at the same world point. Its nonzero snapshot identity binds both providers and their common physical time. The Ciddor model identity remains independent of the varying composition. Prepared fields cache the original wavelength poles once and retain the full source domain; they do not freeze composition at a reference point.

Phase gradients and spatial Hessians apply the complete four-variable material chain rule, including all six mixed state derivatives and both composition Hessians. Native first-order and differentiable ambient constructors, including their prepared forms, reuse the existing curved camera, group optical path and ray-bundle transport. The uniform fast path requires both original providers to guarantee spatial uniformity. Invalid state, wavelength, derivatives and zero snapshot identity refuse explicitly.

Seven focused source-level G1/G3 checks passed using the original material, geometry, scalar GRIN and bundle implementations. They cover mixed spatial derivatives, humidity-only and CO2-only bending/focusing at constant temperature and pressure, the four original ambient bindings, and the emitting-target fixture's path and group delay. The complete native camera/scene regression is authored but has not executed end-to-end in this environment. This remains an explicit nonabsorbing-air approximation; composition refraction does not establish aerosol scattering, cloud microphysics, condensation, atmospheric absorption or time evolution.


## Shared absolute catalog-glass material law

`AirRelativeSellmeierIor` delegates source admission, phase/group evaluation and mapped wavelength support to `fs-material::catalog_refraction::AirRelativeSellmeier`. Its existing constructors, source conventions and renderer identity encoding are retained. `material_law()` exposes the same immutable original law to native electromagnetic consumers; no renderer dependency is introduced below LUMEN. Five existing focused renderer checks passed with this delegation, covering independent reference/coordinate atmospheres, inward support mapping, original Fresnel/thin-lens behavior and identity sensitivity. The catalog reference atmosphere remains independent of the actual scene medium and glass thermal state.

## Finite emitter inventories on the native thermal clock

`heating::luminescence::saturation::{PrescribedPhotonPump,
SaturableEmitterCell, SaturableLuminescentHeating}` connects the original
finite-population kinetic owner to the original enthalpy timeline. Each
binding supplies mean counts [ground, singlet, triplet], an optical cell,
a thermal cell, two molecular cross sections and two prescribed local photon
fluxes. Hazards are cross section times photon flux, in inverse seconds per
available ground-state emitter. The shared positive reaction owner integrates
actual absorbed pump events, ground-state recovery, triplet shelving and all
radiative, nonradiative and intersystem-crossing channels. It does not clamp a
dilute solution or deposit a separate Beer-law pump budget.

The original scalar or material-bound constitutive rates are resolved at the
actual thermal temperature. Rates are fixed over each accepted thermal interval;
step refinement assesses this feedback approximation. Intrinsic lifetimes do
not force lifetime-sized thermal timesteps. Immutable prepared propagators are
reused per emitter and across consecutive matching cells when the model,
duration, pump hazards and both photon energies agree. A pass-local,
constant-memory temperature cache also avoids repeated constitutive queries
for consecutive cells at the same physical temperature. Each cell retains its
own population, thermal lookup and cancellation checkpoint. A bounded pump
window splits native steps at both pulse edges. Initial excitation enters
the conservation account explicitly. Failed or retried thermal steps publish
neither tentative populations nor their accumulated heat/light channels.

Actual pump absorption plus initial excitation equals final excitation, net
native heat, escaped light and explicitly unresolved light. Reverse crossing
draws its state-energy gap from the bath. Spectral relaxation and passive
absorption use the existing original owners. Absorbing emission requires an
explicit passive transport environment; passive absorption creates heat rather
than new emitter populations. Accepted endpoint and last-interval mean spectra
remain bin-integrated watts and feed the existing luminescence observers.

The default prescribed pump field must be justified as optically thin or
externally maintained. The optional depleting beam below supplies its own
incident/transmitted photon balance. Stimulated emission, excited-state
absorption and irreversible photobleaching remain outside this adapter.
Triplet-pair laws and active photon recycling require the explicit bindings below.
Cross sections and rates are authored
inputs, not a material calibration certificate. The example documents synthetic
coefficients, a 1 us pulse and maximum ground-state optical depth 1e-3. Focused
integration tests cover analytic saturation, triplet shelving, exact pulse
edges and afterglow, precharged stores, retry/cancellation/refusal and original
temperature-dependent material rates. Native execution of these renderer
fixtures remains separate from the passing lower kinetic-owner tests.

### Depleting spectral pumps and material cross sections

`saturation::beam::BeamEmitterCrossSections` and
`SaturableLuminescentHeating::with_depleting_beam` bind one directed,
quadrature-weighted two-wavelength beam to the existing finite inventories.
Each wavelength can excite both singlet and triplet states. The original DDA
retains every positive-length cell in ray order, including separately authored
inventories at zero host density. Molecular opacity uses cross section times
ray length divided by the original cell volume, then the complete ground
inventory. Host mass density enters only separately declared passive absorption.
Stable original arithmetic avoids intermediate underflow/overflow in both
geometric factors. No beam area, emitter count or material dataset is inferred.

The shared thermal column owner solves pumping and ground residence together,
with bounded temporal refinement of the complete ordered column. It feeds
actual transmitted photons downstream and keeps wavelength-order beam energies
separate from destination-order excitation energies. Accepted incident,
transmitted, molecular and passive energies have an independent cumulative
audit. Passive pump absorption appears exactly once as external input and once
as native heat; host coefficients must exclude the same active emitters. Every
passive absorbing path cell requires its original thermal receiver, which may
be bound with zero emitter counts. Emitters outside every beam retain their
prescribed local excitation and intrinsic decay.
Passive host coefficients retain the original optical snapshot; the molecular
and kinetic material refresh below does not infer a host temperature law.

With `material-state`, cross sections come from the original
`GroundStatePumpStatePoint` claims at each actual pump wavelength and thermal
cell temperature. Other source coordinates and selected claims remain pinned.
Immutable prepared columns are shared across snapshots while all cell
temperatures agree exactly; changed endpoint temperatures re-query original
support before acceptance. Populations and photon histories are never cached
as constitutive data. Pump windows, cancellation and failed thermal trials keep
the same transactional semantics as the original observer.

This is a cell-uniform, quasistatic, one-pass pump model. Spatial refinement
resolves longitudinal/transverse saturation; native thermal interval refinement
resolves temperature feedback separately from the inner population/opacity
solve. Pump scattering, reflection, stimulated emission and excited-state
absorption are outside its scope. Simultaneous beams use the shared extension
below. Existing
passive emitted-light transport and camera observation remain authoritative.
The synthetic `saturable_afterglow --depleting` example supplies an initial
optical depth of 10 and a 10 microsecond pulse, then reports both afterglow
spectra, native heat, remaining ground opacity and the full pump ledger.

Focused integration fixtures cover dense downstream shielding and reverse ray
order, competing excitation destinations, passive host heating, pulse-off
afterglow, exact pulse boundaries, cancellation/retry and missed/empty beams.
Material fixtures compare both wavelengths with independently authored cold
and warm cross sections, then check endpoint domain refusal and clean retry.
These renderer tests are authored and formatted; their native execution is not
established by the independently passing lower-owner and material-source tests.

### Simultaneous and counterpropagating molecular pumps

`saturation::beam::DepletingPumpBeam` and
`SaturableLuminescentHeating::with_depleting_beams` bind multiple directed
spectral pumps to the same unique finite emitter inventories. Each declared
beam retains its original wavelength nodes, molecular source, native ordered
cell path and optional half-open pulse window. Different paths can overlap,
cross or oppose one another. A cell's complete ground/singlet/triplet state
advances once under all incident channels together; independent per-beam
population updates are not summed. Distinct rays may share a source identity
without deduplicating their physical powers.

The lower thermal owner handles the simultaneous optical/kinetic solve. The
renderer compiles the original geometry and refreshes the original molecular
claims at each beam wavelength and actual shared cell temperature. Its cache
contains only immutable constitutive data, keyed by every unique cell's exact
temperature. Pulse activity, populations and photon histories are excluded
from the cache. The original one-beam entry point delegates to the same binding
and retains its one-column evolution path when its cells have no prescribed
local excitation; co-driven cells use the shared finite-inventory solve.

Each beam's pulse intersects the existing global pump window. Every strictly
future edge splits the native thermal interval, so staggered pulses, overlap
and dark recovery share one accepted clock. Populations, all beam accounts,
heat/light channels and changed-temperature constitutive data commit together;
a failed material query, coupling solve, resource cap or cancellation publishes
none of the trial history.

`depleting_pump_report()` returns per-beam incident, transmitted, molecular and
passive energy in original declaration order with explicit wavelengths and
source identities. Shared kinetic work and temporal discrepancy are reported
once for the complete system. `depleting_pump_energy()` retains its original
single-beam meaning and returns no ambiguous summed spectral account for a
multi-beam system. Passive host heat from all paths is accumulated once into
each original thermal receiver.

These beams are mutually incoherent, with additive intensities. Optical phase,
interference, coherent population dynamics and pump scattering/reflection are
outside this rate model. Active photon return uses the joint extension below. Within a
cell, every beam addresses the same homogeneous inventory; separated beam
footprints or a resolved saturation front require spatial refinement. Material
and native thermal interval limits remain distinct from inner optical coupling
and temporal-refinement controls.

Focused integration tests compare split and combined illumination, mirrored
counterpropagating colors and declaration order, exact per-beam/global pulse
intersections, passive host heat and transactional retry. The existing
`saturable_afterglow --counterpropagating` example uses two original cells,
400/410 nm opposed pumps, staggered pulses, separate beam accounts and observed
afterglow. These new renderer tests and the coupled example are unexecuted:
the native build received SIGKILL and the executor became unavailable during
this extension. Earlier lower-owner passes do not establish this new gate.

### Finite-inventory triplet-pair quenching and delayed fluorescence

`SaturableLuminescentHeating::with_triplet_annihilation` binds the original
`PhosphorescenceAnnihilation` laws to the same finite ground/singlet/triplet
inventories used by prescribed photon fields and one or several depleting
beams. Each pair event consumes two triplets; the authored productive yield
creates at most one singlet, with all remaining emitters returning to ground.
Pair recovery participates inside each controlled kinetic/optical trial,
including during dark intervals and in cells outside every incident beam.
It is not a pair kick applied after a completed optical history.

The pair-event coefficient is in m³/s and the active volume is explicitly
supplied by the original law. Neither beam area, thermal mass nor material
name supplies a replacement volume or coefficient. A productive event must
supply its singlet energy from two triplets even when its probability is small.
`LuminescenceEnergy::annihilation_heat_j` retains physical pair heat separately
from ordinary nonradiative decay and reverse-crossing thermal withdrawal.
Temporary positive-reaction coordinates never become fictitious bath work.
Integrated singlet emission includes delayed fluorescence after the pump stops
and feeds the same accepted endpoint/mean observers and passive light routing.

`with_triplet_annihilation_sources` also accepts
`saturation::annihilation::EmitterTripletAnnihilation::from_material`. This
preserves the original `TripletAnnihilationStatePoint` claims and fixed source
conditions while resolving the pair coefficient and productive yield at each
actual native temperature. The source must own an absolute-temperature axis;
the active population volume remains caller-owned. Unchanged cell temperatures
reuse accepted laws, and only a successful endpoint replaces cached laws,
populations and energy accounts. Unsupported endpoints, cancellation and
failed nonlinear/optical trials publish no partial history.

Zero coefficients retain the original linear results. Existing nonlinear and
optical budgets bound work; componentwise temporal refinement is an accuracy
indicator, not a certified continuum error bound. Exhaustion refuses instead
of changing a tolerance or clipping populations. These laws do not establish
spatial encounter/diffusion physics, stimulated emission or measured
coefficients. The existing example accepts
`--counterpropagating --annihilation`, retaining its synthetic coefficients,
original cell volumes and separate beam and pair-heat accounts.

Four authored renderer regression groups cover independent radiative/quadratic
triplet decay, prescribed-field equilibrium, combined/overlapping beam pumping
and delayed-light observation, and unilluminated cells with zero-law parity,
retry/cancellation and resource refusal. Two sourced groups cover native
temperature refresh and endpoint refusal with clean retry. These renderer
tests are source-formatted; native renderer execution remains unverified after
the unchanged `asupersync` dependency was killed before renderer typechecking.

### Finite-ground emitted-photon recycling and accepted active opacity

`SaturableLuminescentHeating::with_photon_recycling` connects original emitted
photons, changing ground-state absorption, external depleted beams and optional
triplet-pair encounters in one native history. The coupled constructor
`try_new_with_photon_recycling` also admits passive host absorption directly,
without first constructing an unused fixed-probability transfer operator.
Every positive-density host cell needs its original thermal receiver; a zero
emitter inventory can supply that binding. Prescribed local fields, declared
depleted beams and initial stores may supply excitation together, with each
original emitter inventory advanced once.

`PhotonRecyclingGeometry` explicitly bounds sampled rays and the complete
path/segment allocation, seed, logical source and passive-host temperature
domain. The original uniform-cell isotropic sampler and original DDA compile
complete paths to the open boundary. A traversal refusal discards the binding;
no unfinished ray becomes escape or heat. Spatial/angular samples are reused
across physical stages and wavelengths. Their finite quadrature has its own
ray-count/seed sensitivity; it is not an unbiased nonlinear estimator or an
empirical confidence bound. This ballistic binding does not infer scattering,
reflecting/refracting interfaces, coherent interference or photon flight time.
The explicit elastic-scattering binding below retains its own original optical
snapshot and wavelength-dependent paths.

`EmitterAbsorptionCrossSections<B>` is the same molecular source for external
and internally emitted photons. The existing `BeamEmitterCrossSections` name
is its two-lane alias; `RecyclingEmitterCrossSections<M>` follows the actual
emission nodes. Original material claims are queried at those vacuum
wavelengths and each native cell temperature. Source excitation support and
destination-state energy are both checked. The passive host must exclude the
same active molecules. Molecular opacity is `N_ground * sigma / V` on the
original cell; neither host mass density nor an invented beam area replaces it.

`SaturableRecyclingLimits::maximum_secant_history` bounds the inverse-Broyden
history (default eight; zero selects safeguarded physical source iteration).
Coupled volumes may exceed sixty-four radiating states while retaining the
original per-cell kinetic limits and aggregate work budgets. Accelerator
storage and application scale with source count times retained history;
geometry and spectral quadrature retain their separate explicit bounds.

The existing shared-beam kinetic owner integrates all unique inventories once.
Emitted spectral/directional rates are solved from that same history's
radiative events, and absorption uses their actual receiver ground
availability. External and recycled absorption, passive heat and excitation
excess are attributed from original channel events. Only external photons
enter the pump ledger. Source Stokes heat, receiver relaxation, ordinary dark
decay, pair heat and bath withdrawal retain separate physical meanings.
`recycled_excitation_j` counts actual internal excitation-store input; gross
emission/recycling can exceed external energy without becoming extra supply.

Prepared caches contain immutable paths and constitutive data keyed by source
identity and every exact native temperature. Changing beams or pulse bindings
invalidates preparation without caching photon traffic or populations. Native
acceptance updates populations, energy, per-external-beam accounts and source
temperatures together. Candidate retry, cancellation, missing material support
and exhausted work caps publish none of that history.

`recycling_emission_grid` binds accepted endpoint powers to accepted molecular
absorption and the same original host for the existing native camera integral.
Both foreground attenuation and emission within each cell use this opacity.
Endpoint opacity with interval-mean powers is not claimed to integrate an
exposure through changing opacity; resolve observation times explicitly. The
extended example accepts `--counterpropagating --recycling --annihilation`,
uses explicit synthetic spectra/cross sections and reports internal excitation,
secondary heat and escaped energy separately from each external beam.

Three authored renderer groups cover transparent/zero-active analytic decay,
passive host heat, transfer into an unpumped original cell, accepted-opacity
observations, shared beams with pair recovery, and cancellation/retry/budgets.
Two material groups cover original temperature-dependent emission absorption
and endpoint domain refusal with clean retry. These renderer fixtures and the
example are source-formatted and source-reviewed; native execution remains
unverified because the unchanged `asupersync` dependency was killed before
renderer typechecking. Lower-owner native tests do not establish this gate.

### Elastic scattering coupled to finite-ground photon return

`SaturableLuminescentHeating::{try_new_with_scattered_photon_recycling,
with_scattered_photon_recycling}` binds the original `ScatteringGrid`, including
its cell-specific elastic spectra and reciprocal phase laws, to the same finite
emitter/thermal history. Each emission wavelength samples complete independent
elastic walks from the existing uniform-cell isotropic source. The original
native free-flight and phase samplers are reused; prepared optical palettes are
amortized over all histories. Zero scattering delegates to the original
ballistic binding and retains its arithmetic, samples and source identity.

The geometry proposal contains scattering only. Passive and evolving molecular
absorption act as survival factors along the retained walk, so no ground-state
absorption probability is frozen into the geometry. This decomposition uses
one monochromatic proposal per path: native collision/phase weights are one,
with no cross-wavelength mixture correction. Every clipped original-cell
traversal remains ordered, including zero-host-density cells and repeated
returns. A repeated cell uses its original volume and single shared population;
nonadjacent visits are never collapsed across intervening passive absorption.
Source singlet/triplet spectra lower separately onto those walks, with their
authored photon fractions and distinct source/destination energy accounts.

The explicit per-history scattering-event limit and native DDA-visit limits
refuse unfinished walks. Aggregate path/segment bounds include every wavelength
and state copy plus external beams at native admission. Reaching a cap does not
make a surviving path escape or deposit its residual as heat. Cancellation or
any refused binding publishes no kinetic state. The complete scattering/phase
snapshot, collision bound, sample seed and geometry enter source identity.
Its declared temperature domain covers the entire immutable passive optical
snapshot; this adapter does not invent a temperature law for scattering.

Fixed sampled geometry is reused through the nonlinear kinetics. Angular,
spectral and spatial sampling sensitivity remains separate from kinetic and
thermal error control; the nonlinear result is not an unbiased estimator or
a certified Monte Carlo interval. Elastic scattering does not alter photon
energy. Interfaces, coherent propagation, directional dipole emission and
finite flight delays remain outside this open index-one volume model.

The existing `saturable_afterglow` example accepts `--scattered-recycling` with
the prescribed, depleted or counterpropagating pump modes and optional local
background/pair kinetics. Its explicitly synthetic spectrum has zero elastic
scattering at both pump colors and nonzero scattering at the two emission
wavelengths. The same snapshot drives photon return and the accepted-opacity
endpoint observation; direct pump transport therefore remains consistent with
that declared spectrum. The endpoint uses 128 bounded native photon histories,
while recycling uses 32 complete walks per cell and wavelength. Those sample
counts demonstrate the API and do not certify sampling convergence.

Authored focused fixtures check wavelength-dependent collision distances,
clipped returns through vacuum, independence of the walk from passive
absorption, exact zero-scattering parity, secondary excitation/native heat,
deterministic replay and aggregate/event-budget refusal. The original recycling
fixture now explicitly authors passive spectral support through its 800 nm
emission node. These Rust fixtures are source-reviewed and unexecuted: Rust,
DSR and RCH are unavailable in the current workspace. Independent ordered-path
fluence/energy equations and source review do not establish the renderer gate.

### Maintained local excitation with beams and photon return

Existing `PrescribedPhotonPump` bindings may now remain active when depleted
beams or finite-ground recycling are bound. The renderer passes their original
hazards and photon energies to the shared native population owner before
advancing the cell. Local fields and every transported pump compete for the
same available ground molecules, including optional triplet-pair reactions.
A prescribed field must represent distinct optically thin or externally
maintained illumination; photons already declared in a beam must not also
appear in that field. No beam area or physical illumination is inferred.

The existing global pump window gates local fields and external beams. Each
beam's own narrower window affects only that beam. Every pulse edge splits
the same native thermal timeline; cells outside all beam paths continue to
receive their declared local excitation. Dark continuation closes both
external supplies while retained excitation, recycling and pair reactions
continue on the original clock. Empty local input preserves the previous
numerical routes, including the original single-column fast path.

`absorbed_pump_energy_j()` includes consumed local and external beam photons,
ordered by destination state. Each beam's incident/transmitted/absorbed report
continues to describe that beam alone. The recycling adapter adds the native
`LocalPumpStep` absorption and relaxation directly into external totals;
it never recovers a weak local contribution by subtracting large recycled
traffic. Combined cell kinetics supply heat once. Populations, accepted
material temperatures, source powers and all cumulative accounts retain the
same transactional acceptance and retry behavior.

The synthetic example accepts `--recycling` alone with its existing thin
prescribed pump, or
`--counterpropagating --background-pump --recycling --annihilation` for
maintained local excitation alongside opposed beams. The additional background
is explicitly authored as 400 nm, 5e24 photons/(m² s), giving a 5e4/s hazard per
available ground emitter under the global pulse. Output labels total external
absorption separately from per-beam columns. Coefficients remain synthetic.

Three new renderer fixture groups cover analytic excitation outside beam
paths with global/per-beam pulse intersections, analytic local-only transparent
recycling, and local triplet pumping with a depleted singlet beam and pair
recovery both with and without active return, including retry and dark
continuation. These fixtures and example are formatted and source-reviewed
against the exact new native APIs; full renderer execution remains unverified
because its genuine `asupersync` dependency build is killed before renderer
typechecking. The lower owner separately passed all 28 focused pump tests,
including three new independent physical oracles; those passes do not replace
this full renderer gate.

## Multiply scattered nonthermal luminescence observations

`volumes::transport::luminescence::{ScatteringLuminescenceGrid,
render_scattered_luminescence_camera}` borrows accepted sparse cell-emission
powers and the same complete absorption snapshot extended by its original
elastic scattering and reciprocal phase laws. The native spectral-mixture free
flight and actual collision-cell phase supply continuation; a source hook in
that same DDA integrates the complete ray's `j=P/(4*pi*V)` under total extinction
before applying continuation weights. Original discrete nodes remain
bin-integrated W/(m² sr), with no wavelength width, absorptivity or LTE source
applied to the authored emission. Logarithmic path weights and the existing
logarithmic source fallback preserve representable weighted source radiance.

Accepted sparse molecular absorption is included in the same flight proposal
and complete source integral. In each original cell, the total linear
extinction is `rho*(kappa_absorption+kappa_scattering)+alpha_active`;
the scattering numerator remains `rho*kappa_scattering`. Additional
absorption stays in inverse metres, including in zero-host-density cells.
Spectral continuation divides by the full mixture collision density, so adding
absorption cannot silently retain a host-only proposal. A sampled absorption
event with no surviving continuation still completes the independent source
ray, including emitters behind that event. Direct and scattered observers
share logarithmic total-extinction arithmetic for representable opaque limits.
The sparse snapshot is borrowed without a per-ray opacity-map allocation;
original phase laws, cell identities, cancellation and work bounds remain.

The pinhole/thin-lens image shares original raster/pupil sampling and returns
spectral radiance. Zero scattering retains the original exact integral and draw
order. Logical path IDs and seeds determine independent streams. Per-history
collision and actual DDA-visit limits, plus the whole-image visit limit, are
binding; an unresolved surviving continuation at a cap refuses the entire
sample/image. Refused histories must not be discarded or replaced with zero.
Sampling uncertainty is not certified by a completed image.

This is an open, stationary/quasistatic, index-one observation with uniform
isotropic cell emission; surfaces, flight delay, wavelength-changing scattering
and kinetic/thermal updates remain outside this observer. Seven regression
groups in `luminescence/scattered/tests.rs` cover zero-scattering and camera
parity, independent two-lane one-scatter integrals, actual cell-selected phase,
lossless boundary flux, source arithmetic range, deterministic replay and
admission/cancellation/work limits. Three additional integration groups cover asymmetric spectral absorption,
a zero-density absorbing foreground with a complete source ray, and exact
zero-scattering parity without random draws; the identity fixture also checks
accepted active opacity and rejects a different host snapshot.

The exact production `scattering/extinction.rs` numerical leaf passed all three
focused native tests in the existing actual-source harness: spectral mixture
weights/event distances, independent absorption at zero host density, and
representable distance/weight/source scale when the direct extinction product
overflows. No test was ignored. This does not establish execution of the full
adapter groups: the genuine renderer dependency build remains blocked by
`asupersync` compilation being killed before renderer typechecking. Source
formatting and interface review are reported separately from that native gate.

## Time-resolved luminescence camera exposure

`volumes::transport::luminescence::render_luminescence_exposure` observes an
explicit ordered set of accepted field snapshots through one frozen native
`PhysicalCamera`. `LuminescenceExposureNode` binds each declared physical time
to a direct `LuminescenceGrid` or `ScatteringLuminescenceGrid` and a positive
dimensionless quadrature weight. Weights must sum to one within 64 machine
epsilons; finite increasing shutter bounds and strictly increasing in-window
times are required. The original source identity and exact wavelength order
must agree, while complete field identities may change with accepted powers,
opacity and material state. The caller supplies the accepted snapshots and
their times; this observer does not create or certify a kinetic history.

Each pixel/pupil ray transports every node's actual source and absorption
before averaging. Transporting independently averaged powers and opacity is
not substituted for this operation. The original raster kernel, pupil draws,
path seeds and native direct/scattered owners remain in use. Common random
numbers correlate observations across time without altering their individual
path measures; constant fields retain the frozen camera's samples and result.
One output image and a bounded temporal weight vector are allocated, rather
than a full intermediate image for every time node.

`LuminescenceCameraExposure` retains the camera, sampling settings, shutter
bounds, original source/wavelengths and an ordered field-history identity. Its
mean is bin-integrated radiance in W/(m² sr). The optional time integral is in
J/(m² sr), and the existing approximate visible observer projects the mean.
These are radiance quantities, without inferred sensor irradiance, aperture
gain or display calibration. The compensated positive weighted mean retains
both old and new fractions independently, preserving a short bright pulse
even when the dominant fraction rounds to one. Mean radiance remains usable
when its time integral overflows or positive energy underflows to zero; the
integral query explicitly refuses those unrepresentable results.

`LuminescenceExposureBudget` admits temporal node storage and the full
pixels-times-primary-samples-times-nodes query count. Every actual temporal DDA
visit consumes the same whole-image allowance; checked totals preserve the
native per-history bounds and prevent visit-counter overflow. Cancellation or
any refused observation returns no partial exposure and changes no field.
The positive quadrature is an estimate with caller-owned temporal resolution,
not an adaptive integration or temporal error certificate. Backgrounds, camera
motion, photon flight delay and sensor response remain separate models.

Four authored fixture groups cover midpoint convergence to an independent
exponential light integral with changing emission and foreground opacity;
constant direct/scattered camera parity with actual multiplied work; bright
short pulses in either weight order and integrated-radiance representability;
and source/spectrum, time, quadrature, budget and cancellation refusals. Scalar
calculations independently verify second-order midpoint convergence and the
pulse/correction arithmetic. These new Rust fixtures are source-reviewed but
unexecuted because Rust, DSR and RCH are unavailable in the current workspace.

## Bounded opaque GGX reflection sampling

Opaque legacy GGX, spectral conductors and body-locked brushed conductors use
the bounded spherical-cap VNDF of Eto and Tokuyoshi, equations 5, 6 and 8:
<https://gpuopen.com/download/Bounded_VNDF_Sampling_for_Smith-GGX_Reflections.pdf>.
The cap contains every contributing reflection; anisotropic slopes use the
conservative smaller-slope bound. Sampling and both forward/reverse PDFs share
the same view-conditioned cap. Polarized samples and reciprocal connections
retain its normalization analytically as `cap_mass * G2/G1(wo)`, so the original
BRDF, optical phase, Fresnel law and shading-normal adjoint remain unchanged.
Residual rejected facets contribute zero without resampling. The warp uses the
same two uniforms, bounded constant work and no allocation. Transmitting glass
retains its original full VNDF. Replay/checkpoint semantics advance to tracer
10, conductor 3 and brushed 4; dielectric semantics remain 20.

Five G1/G3 regression groups in `tracer/ggx_sampling_tests.rs` pass using the
original production helper sources, original `Vec3` operations and native
`fs-math`. Independent angular quadrature checks proposal mass, directional
moments, anisotropic furnace power and reciprocal normalization. In the
unit-Fresnel, unit-roughness, normal-view single-scattering furnace, the mean
remains the analytic `1 - ln(2)`. Equal-draw variance falls from approximately
0.1332526 to 0.0195470 (6.82 times lower), and rejection falls from half of draws
to zero. This is an estimator-variance result, not a wall-clock speed claim or
inter-facet multiple-scattering closure. Existing native scalar/polarized scene
regressions also exercise the changed owners; their execution remains pending
the unavailable constellation dependency checkouts in this environment.

## Positive narrow-band temporal interference

Scalar exposure, full Jones coherency, analyzer projections and scalar/Jones
interferometers share `fs-fft::harmonic_window::PreparedUniformHarmonicBand`
when the original uniform shutter admits its bounded narrow-band factor.
The same temporal modes determine diagonal powers, complex cross-arm moments
and polarized coherency. Original source covariance columns remain independent;
equal-frequency fields still combine before projection. This retains small
positive multi-tone null signals without subtracting nearly equal powers or
clamping a negative result. Wide-band and one/two-tone observations retain the
existing analytic pair kernels. Spatial quadrature and authored propagation
phases are unchanged. The shared owner documents its absolute truncation bound;
this is not a relative floating-point accuracy certificate at arbitrary nulls.

Preparation preflights the selected representation before temporal allocation.
For N distinct frequencies and M retained modes, scalar power costs S=M*(N+1),
complex mutual coherence C=M*(2*N+1), and full Jones coherency H=M*(2*N+4) per
source bank. Two active scalar arms cost 4*S+C; two active Jones arms cost
4*S+2*H+2*C. A single active arm costs S or H respectively. Counts include all
independent source banks and every spatial quadrature sample. Observations use
fixed scratch, no allocation, and bounded cancellation; no carrier-time raster
is introduced. Regression fixtures retain positive fourth-order three-tone
nulls through the original chromatic pupil, source-column identity, full Stokes,
passive port closure, finite-pixel energy and exact work-budget boundaries.

## Anisotropic Smith reflection in native brushed scenes

The opt-in Smith forward and reciprocal camera paths now apply the complete
anisotropic height walk to `BrushedConductor`, including its original body-locked
radial/circular axes and roughness mark. The finish is evaluated once per
geometric encounter; the existing Smith kernel and original complex
Fresnel/Mueller owner handle every hidden facet. Zero detail and the undefined
radial pole reuse the isotropic Smith path exactly. Microscopic heights do not
move the geometric endpoint or add path length, delay, absorption or a coating.

For coincident shading/geometric normals the anisotropic walk is passive.
Authored normal tilt remains the original empirical shading approximation:
geometric clipping applies only after complete escape, and forward power
receives the single endpoint adjoint
`(ns.wo/ng.wo)*(ng.wi/ns.wi)`. Reciprocal camera continuation uses the native
radiance weight. Internal facets receive no extra G factor, horizon clipping
or adjoint. Clipped escape contributes zero without retries; cancellation,
invalid data and facet exhaustion refuse. This is not a resolved-groove heat
model or an outgoing marginal PDF suitable for MIS.

Six focused Rust tests (three existing and three new) pass with the actual Smith
kernel, native mirror/frame and conductor/transparent-film/absorbing-film
Mueller math, original finish/rigid-transform math and `fs-math`. They cover
facet replay, spectra/media, body rotation, mark and exact fallback, endpoint
adjoint, final clipping, cancellation and budgets. The focused executable
uses minimal execution/hash and homogeneous-medium wiring; it is not a full
workspace build. A native camera regression additionally covers analyzer,
Beer attenuation, source-before-scatter ordering and anisotropic scene dispatch;
execution of that full scene test remains pending dependency availability.

## Point-fixture illumination with complete Smith reflection

`conductor::microsurface::SmithGgx` supplies fixed-direction collision
connections and an analytic height-integrated first reflection. The polarized
`evaluate_smith_conductor_connection` and
`evaluate_instanced_smith_brushed_conductor_connection` APIs retain that exact
first-order term and sample the higher orders with the original Smith walk.
Every connection uses its own smooth complex Fresnel/Mueller response on the
retained incoming optical prefix. Its escape probability is applied once.
Native coherent films, local phase index, source polarization and frame
reciprocity remain authoritative. Brushed axes/marks retain the original
instance chart. An authored tilted shading normal keeps the existing empirical
radiance law and its geometric/shading cosine conversion once at the endpoints.

`point_lights::SmithPointLightCameraScene` connects the original calibrated
fixtures on the complete Smith camera continuation. Native bundle shooting,
physical irradiance, optical source polarization, Beer/GRIN propagation,
receiver geometry, and camera n² conversion are reused. Analyzer and resolved
sequential-lens seeds enter before every facet. Source and thermal emission
remain on the same native camera path, including its post-microsurface prefix.
The added Smith evaluations use fixed local state with no heap allocation;
work is linear in the actual bounded facet count. Existing bundle owners retain
their workspace policy. Finish frames and optical tables are not rebuilt per facet.

Point fixtures cannot be hit by continuous path continuation, so their selected
connections add without a competing endpoint strategy. The directional
response is a stochastic BRDF contraction, NOT an outgoing marginal PDF.
Neither a sampled response nor the hidden-path joint density is used in MIS.
Surface decisions, continuing facets and each source-connection evaluation
have independent logical counter domains. The first-order connection is
integrated rather than counted again as a random collision contribution.

Facet limits and original cumulative source/shooting limits bound the whole
sample. Any unresolved path, source connection, numerical representation
failure or cancellation refuses it; averaging only successful histories would
bias the result. This implements the original declared point-source branches,
not all optical images, area/environment NEE, BDPT or molecular scattering.
Statistical micro-heights add no metric travel, delay or coherent path phase.

Twenty-three focused Rust tests pass with the actual geometric walk,
conductor/film Fresnel-Mueller mathematics, native beam/frame and reciprocal
contraction bodies, original finish/placement mathematics, fs-math and Philox.
Eight new groups cover independent analytic first order, hidden-height
integration, unit-furnace and absorption agreement with sampled escapes,
anisotropic reciprocity, actual Mueller prefixes, frame/scale/fallback behavior
and transactional limits. Minimal execution/hash/homogeneous-medium wiring
supplies the isolated dependency closure. Two full native point-camera
regressions are authored for physical irradiance, Beer attenuation, independent
connection/continuation draws and refusal; those scene tests remain unexecuted
with the missing constellation checkouts.

Two additional compiled integration probes execute the actual new point-camera
and receiver modules and original fixture accumulation/photometry seams with
signature fixtures for the driver, bundle and evaluator. They verify isotropic
and brushed dispatch, radiometric composition, independent sample indices,
work counters, contact refusal and cancellation. These are integration-shape
checks, not execution of the full native scene/geometry solver.


### Original tabulated thermal refraction in native rays, meshes and scenes

- [S] The material-state adapter now exposes
  `MaterialSpectralRefractiveField` and a borrowed wavelength preparation.
  The shared `SpectralDielectricTemperaturePath` owns the original selected
  phase, wavelength and intensity-absorption claims and their exact temperature
  interpolation. A local sample uses all three at the same accepted temperature;
  no additive reference-index correction, source refit or resampling is applied.
- [S] The native reference glass must equal the complete original tabulated
  spectrum, including phase/loss/provenance and the interpolation-gap policy.
  The same gap ceiling is checked over the full source temperature interval,
  including interior wavelength breakpoints. Wavelength support is the exact
  common original band, and snapshot identities must be nonzero.
- [S] Prepared ray stages reuse the original smooth thermal cell without
  allocation, material queries or finite differences. Analytic moving-grid
  derivatives compose as grad(n)=n_T grad(T) and
  Hess(n)=n_T Hess(T)+n_TT grad(T) grad(T)^T, including a nonzero Hessian when
  the central temperature gradient vanishes. Phase optical length, available
  group delay and Beer loss use the same original state.
- [S] `P1ConductionMesh::propagate_spectral` and
  `propagate_spectral_bundle` reuse the native element cursor, prepared
  walkers, geometric budgets, transition collars and gradient-jump mapping.
  They require the exact accepted mesh and nodal slice; equal-valued foreign
  snapshots are refused. No additional propagation heap workspace or solver
  is introduced.
- [S] Four native material-field scene constructors retain gradient or full
  Hessian capability in prepared/unprepared form. Their identity binds the
  complete reference, original material path and accepted snapshot, and the
  existing full reference-medium comparison remains authoritative.
- [S] Added native regressions cover exact affine curved rays and bundles with
  phase/group/current loss, moving-grid nonlinear Hessians, original source
  cusps and gap/domain refusals, and four mesh orientation/gradient-jump cases
  with exact piecewise-affine ray references. Independent scalar reference
  arithmetic and bounded source/API review were performed. Rust compilation,
  native regressions and rustfmt were unavailable in this session.
- [F] Smooth preparation refuses genuine active source cusps or interior moving
  spectral-knot crossings; callers must split the declared thermal interval.
  The accepted geometry and state stay frozen during optical transport. This
  adapter infers no source extrapolation, scattering, stress optics, expansion,
  moving-medium correction or causal transient constitutive storage.

## Physical triangular shutters for coherent optical observations

`CoherentDiffractionShutter::try_with_profile` and
`try_with_profile_from_clock_dd` admit Uniform and Triangular physical shutters
with unit peak transmission. `profile()` exposes the authored shape and
`effective_duration_s()` gives integrated transmission: T for Uniform, T/2
for Triangular. A zero unrepresentable transmission integral is refused.

Reported powers and coherencies are normalized transmission-weighted means.
Every scalar/Jones pixel and point/pixel interferometer energy uses the same
physical transmission integral. This includes analyzer outputs, signed Stokes
energies and both passive interferometer ports. A triangular detector's
sensitivity is applied once to the observable, rather than squared as a field
taper or multiplied by the full elapsed duration.

Shared profile-aware positive temporal factors now serve scalar, Jones,
analyzer, mutual-coherence and interferometer queries. Original source columns,
chromatic pupil/sensor phase, precise source epoch, spatial quadrature and
explicit work limits remain authoritative. Unsupported bands retain the original
analytic pair kernels; no time raster or new optical transport model is added.

All 101 functions in the current coherent_diffraction target pass: 96 existing
tests and five new groups for physical throughput across every energy output,
independent polarization-beat integration, fourth-order scalar/Jones nulls
through 1 ps, retained DD epochs, and admission/cancellation. The execution
uses the actual complete non-tracer camera/diffraction/coherent source graph,
original geometry vectors and actual math/harmonic code, with narrow Cx/Ray/
ShutterInterval boundary carriers. Full tracer/material-dispersion, moving-camera
and workspace/DSR integration are outside this executed validation.

### Complete opaque Smith finite-area camera connections

`area_lights::camera::SmithAreaLightCameraScene` borrows an admitted area-light
camera and opts rough/brushed conductors into complete Smith continuation and
fixed-direction reflection evaluation. Each selected receiver/source pair must
share one homogeneous medium. It samples each original area/face selection and
emitting side with its native area PDF, then launches the unique straight
source-to-contact direction. `BundleConnectionTarget::try_direct_surface_receiver`
retains the original receiving contact and physical bundle projection while
every intervening native surface blocks this direct proposal, including smooth
mirrors, lossless interfaces and rough scattering surfaces.

The estimator assigns the represented DIRECT source endpoints entirely to area
sampling. Only the immediately following continuation hit of the same original
source object, selected face set, physical side and medium is suppressed.
An intervening surface retires this partition, so reflected/refracted source
chains retain their complete continuation contribution. Primary surface
emission, unselected sources/faces/sides and LTE on all camera segments remain
with the original camera emission owner. Existing source admission preserves
the exact calibrated spectral/polarization owner or directional thermal/coating
snapshot and rejects overlapping duplicate source selections.

The fixed-direction full Smith factor is stochastic; it and the hidden-walk
joint density are never used as outgoing PDFs or nonlinear MIS inputs.
Independent counter domains separate continuing facets, area samples, bindings
and source polarization channels. Native source s/p replay, complex per-facet
Fresnel/Mueller transfer, Beer attenuation, local refractive-index scaling and
work receipts remain unchanged. Analyzed and resolved real-lens beams retain
their original reciprocal Stokes seed. Bounded geometry, source replay or
facet exhaustion refuses the entire sample, as does cancellation.

This direct source proposal does not cover GRIN legs, refractive/mirror images
or dielectric microsurface multiple scattering; such geometric source chains
remain on the original continuation path. Tilted brushed normals retain the
existing empirical shading law, not resolved groove geometry or a new thermal
emissivity model.

### Full Smith reflection under directional environment illumination

`environment::sampling::SmithEnvironmentCameraScene` borrows the admitted
native environment camera. Rough and body-locked brushed conductors use the
same complete polarized Smith law for source connections and continuation;
ordinary materials retain their original optical law and exact-PDF MIS.
Analyzed and resolved real-lens launches preserve the reciprocal Stokes seed.

The complete microsurface outgoing PDF is unavailable. The new geometry-only
queries report the analytic density q1 of escape after exactly one reflection
at unit facet albedo. Its integral is the single-reflection probability, not
one; it is never substituted for full-walk throughput or claimed as a marginal
PDF. Given the original environment proposal p, the light contribution uses
p²/(p²+q1²) and the escaped continuation uses the complementary direction-only
weight. Both are deterministic functions of the arriving macrostate and source
direction. Linearity therefore retains the complete BRDF integral while keeping
narrow-lobe continuation useful; no sampled BRDF or hidden joint path density
enters a nonlinear weight. This is not a claim of globally optimal variance.

The original environment texel owns direction, angular density, calibrated
spectrum and source polarization. Shared native visibility, contact ownership,
radiometric/index scaling and work bounds remain in the original source owner.
Only direct homogeneous nonabsorbing legs compete. Unsupported Smith legs
consume no environment/connection samples. Every later surface retires the
previous partition; smooth delta chains keep their complete endpoint term.
Surface and LTE emission remain unweighted by this environment partition.
The opaque source wrappers do not silently opt rough glass into a different
BSDF; full dielectric microsurface transport has its separate explicit owner.

The full polarized connection evaluator also retains isolated grazing escape
terms through a logarithmic geometry fallback and positive log-tail summation.
Normal contributions keep their compensated linear fast path. Subnormal
positive terms accumulate before linear rounding; only a complete positive
response that is still unrepresentable refuses. Intermediate underflow is not
declared physical extinction. The density endpoint query performs no optical
solve, random walk, allocation or material query.

Ordinary-surface, environment, continuing-facet and evaluating-facet streams
have independent counter domains and zero-based surface/facet ordinals.
Original source/query limits and per-estimate facet limits bound all work;
cancellation or exhaustion refuses the complete sample. Averaging only
successful histories after a refusal would bias the result.

Verification: 27 focused native numerical tests passed, using original
fs_math/Philox, Smith geometry, complex conductor/coating Fresnel-Mueller
responses, beam frames and anisotropy. Four new groups cover the exact
one-reflection integral, brushed-frame deterministic density, 32,768 independent
complete polarized path/source estimates, and grazing/subnormal accumulation.
Three compiled integration probes additionally exercised the complete modified
environment sampling module, new view, receiver adapter and original source/
weight helpers with explicit scene/visibility/emission/signature fixtures:
ordinary/Smith dispatch, endpoint partition resets, primary/LTE preservation,
blocked and unsupported legs, cancellation, and facet/geometry refusals.
Four native full-scene regression groups are authored but unrun. Full workspace
Cargo and resolved-geometry runtime validation remain unavailable.

### Full two-sided rough dielectric Smith transport

`PolarizedBeam::sample_smith_dielectric` and
`sample_smith_stack_dielectric` add the Heitz et al. (2016) dielectric random
walk to the original native Smith geometry. Every microscopic collision uses
the original spectral smooth Fresnel/coating Mueller response and ray-frame
owner. Refraction changes the active medium and reverses the same statistical
height: for `d = -ln(C1(h))`, the opposite side is
`-ln(1-exp(-d))`, evaluated without losing small complementary probabilities.
Both sides retain the authored anisotropic slope distribution. Transmission
therefore has the shared-height beta-function masking law from the original
model; it is not approximated by the reflection masking expression used in
the older single-scattering dielectric path. TIR phase, multiple reflection
and transmission returns, reversed compact or arbitrary passive film order,
and local incident/exit phase indices remain in the native optical owner.

The caller supplies independent `[flight, facet_u, facet_v, Fresnel_branch]`
variates. The finite facet budget and cancellation refuse incomplete samples.
Only final escape determines the macroscopic event and scene medium identity.
Each microscopic response carries ray power; no radiance eta multiplier,
albedo compensation, extra Smith G factor, fictitious metric distance, delay,
or bulk Beer attenuation is applied inside the statistical surface. Coherent
film absorption remains local to each facet. The joint log density is in
per-collision depth and facet-normal measure with branch/escape atoms. It is
never an outgoing marginal PDF. A bare equal-index boundary is an exact null
event with no draws and preserved frame. At a coated equal-index boundary,
paths with no facet reflection retain the straight-through delta atom of the
mixed distribution instead of claiming a continuous rough transmission PDF.

`SmithDielectricCameraScene` is an explicit path-only camera owner for full
opaque and dielectric Smith transport. The curved scene binds physical body
axes once per hit, retains source/GRIN/Beer propagation and optical distance,
and changes topology only on net macroscopic transmission. The existing
camera's reciprocal Stokes launch and endpoint n-squared conversion remain
the sole radiance conversion. Seeded spectral, passive-analyzer, resolved
lens, and four-wavelength/pupil raster paths use independent logical facet
streams; variable microscopic work cannot shift other surface or pixel
draws. Existing opaque point, environment and area-light connection wrappers
retain their original dielectric law until they have a compatible full
dielectric connection estimator. No hidden-path PDF is inserted into MIS.

Regression coverage includes the independently integrated transmission beta
law, two-sided polarized white furnaces, cosine-integrated n-squared
transmission reciprocity, genuine refractive returns, analytic TIR ellipticity,
reversed lossy films, arbitrary-stack/compact-film agreement with rotated
anisotropy, equal-index mixed measures, deterministic replay and explicit
limits/cancellation. Native scene regressions cover final-medium/Beer/radiance
and geometric-distance accounting plus smooth camera/raster parity.

Verification in this session: seven production optical groups passed with the
original dielectric/coating/anisotropy, interface Mueller, smooth beam scatter,
Smith geometry, fs_math and Philox sources; external execution/container/hash
wiring used explicit fixtures. An eighth focused probe executed the exact new
curved dispatch block with fixture context containers and native glass optics.
Actual camera/raster modules typechecked against signature fixtures. The two
full native scene regressions are authored but unrun; the absent sibling
workspace graph prevents full Cargo, resolved geometry/GRIN, raster, clippy,
material-state or remote-runner validation. These are numerical-model checks,
not experimental validation of arbitrary measured rough glass.

## Original optical finishes across accepted material states

`MaterialOpticalBinding::with_surface_finish` binds an explicitly authored native
finish while retaining the exact accepted bulk optical law and original source
receipts. Visible and spectral accepted-state refreshes retain the full
dielectric tangent, anisotropic GGX law and compact coherent coating, or the
conductor's isotropic finish and transparent, absorbing or periodic coating.
Fresh bulk coefficients continue to come from the original pinned material
claims at the complete accepted state.

Independently authored coating laws and thicknesses stay fixed until explicitly
replaced; bulk temperature changes do not infer surface evolution. Original
per-query wavelength support, tangent-frame admission, scalar reflection and
polarized transport remain authoritative, including zero-thickness coating
support. The existing material content identity covers the complete finish.
Bulk-law, optical-family, source and zero-surface-identity substitutions refuse
without changing either the original binding or accepted physical state. The
Euler scene bridge preserves this refusal through its existing material-state
error.

Focused native fixtures cover visible and spectral endpoint continuation,
coated anisotropic dielectric scalar/polarized reflection, all three original
conductor coating families, retained source and interpolation policy, and
independent narrower coating support. They do not establish measured finish
calibration, bulk birefringence, inter-facet multiple scattering, or inferred
thermal evolution of the coatings.

### Full rough-glass fixed-direction polarized connections

`PolarizedBeam::evaluate_smith_dielectric_connection` and its arbitrary-stack
variant evaluate every reflection and transmission order of the original
two-sided Smith walk toward one specified source-arrival frame. The first
height is integrated analytically: reflection uses the correlated Smith sum;
transmission uses Beta(1+Lambda_i,1+Lambda_t). A bounded positive-argument
log-Beta recurrence and combined Stirling expression retain grazing arguments
without subtracting large log-Gammas. Later sampled collisions use the native
full visible-normal law, the reflection/refraction solid-angle Jacobian and
escape from the appropriate side of the same height. Every forced and sampled
facet retains native coherent spectral Mueller transfer, reversed film order,
TIR phase and authored anisotropic axes.

Microscopic continuation carries ray power with its original scalar branch
proposal. The completed connection includes the net
(n_camera_arrival/n_source_arrival)^2 radiance conversion exactly once; the
camera's prefix index conversion remains separate. The stochastic response has
no outgoing PDF. Neither it nor the hidden-path joint density may enter a
nonlinear MIS weight. The explicit capability query admits only rough
unequal-index interfaces; smooth/bare index-matched surfaces have zero
continuous response, and rough index-matched coatings/stacks refuse direct
fixed-direction evaluation with UnsupportedMixedMeasure. Their existing full
mixed-measure path sampler remains available.

The shared positive accumulator retains logarithmic/subnormal terms until
complete summation. A wholly unrepresentable positive result still refuses
NumericalRange; one deterministic inside-to-outside grazing replay is covered.
First-order analytic integration removes the ordinary grazing first-height
underflow. Tests never discard refused histories: all-order integral comparisons
use the same declared source support abs(cosine)>=0.02 for both strategies.
No universal finite-range, globally optimal variance, measured-surface or
arbitrary grazing-success claim is made. Work is bounded and cancellable; facet
exhaustion refuses the complete estimate, not a truncated dark contribution.

Verification: 41 focused tests passed against original optical/math/Philox
code with explicit execution/vector/container/hash boundary fixtures: seven new
connection groups, seven existing full-glass groups, and the previous 27 opaque
Smith/geometry/Philox groups. New checks include independent Beta quadrature,
262144 paths per incident side with separately compared reflected/transmitted
polarized source integrals, the matching deterministic native-proposal angular
partition, scalar n-squared reciprocity, lossy layered source integrals and
polarized linearity, compact/stack equality under rotated anisotropy, near-index
contrast, high-precision log-Beta references, cancellation, replay, budgets and
the explicit grazing range boundary. Full workspace Cargo/DSR and native scene
geometry execution remain outside this focused probe.

### Calibrated point fixtures with complete rough-glass response

The original full-glass raster now exposes six additional detector routes:
ordinary and resolved-lens QE energy/photon measurements, their spectral-analyzer
variants, and ordinary/resolved-lens spectral-analyzer CIE measurements.
All channels share the original four wavelength paths. The original five
surface/facet/source RNG domains, source receipts, work accounting and
receiving-medium index conversions are unchanged. A detector's full wavelength
band is admitted before tracing, including a fully clipped lens packet.

A spectral analyzer is installed at the sensor before the coated lens and
rough-glass reciprocal polarization response. Steady complete-lens beams do
not request an unavailable group clock. Detector and pupil weights are applied
once; a clipped lane keeps its original quarter proposal mass and zero source
work. Two native integration groups cover anisotropic glass from both media,
source polarization, replay/work parity, nonconstant QE/photon arithmetic,
coated-lens analyzer ordering, clipping and band refusal. Changed Rust was
parsed and independently source/API reviewed; these full-scene fixtures were
not executed in this editing environment.


`SmithDielectricPointLightCameraScene` explicitly selects complete opaque and
two-sided dielectric Smith transport for the original point-fixture camera.
Its shared receiving-contact adapter validates both original adjacent medium
identities, binds anisotropic body axes once, and retains the actual compact or
arbitrary spectral coating stack. Reflection and transmission source legs use
their own native receiving plane, source polarization, bundle irradiance and
Beer/GRIN propagation. The full-glass connection owns its net local radiance
index ratio; the original camera-prefix ratio is applied exactly once.

The original source solver, cumulative work, surface/LTE emission and finishing
owner are shared with the opaque point-light view. Seeded and caller-owned
ordinary, analyzed and resolved-lens entry points use separate ordinary-surface,
opaque/glass continuation and opaque/glass connection streams. This new API's
surface and facet identities are zero-based, independently of microscopic work.
A point fixture has no competing continuous endpoint proposal; no stochastic
response or hidden-path density enters MIS. Selected equal-index coated rough
receivers return the core's explicit mixed-measure refusal. Smooth and bare
null boundaries have no continuous receiving response.

Unresolved connections, cancellation, facet exhaustion or numerical range
failure refuse the complete camera sample. This covers the original declared
source branches, not all caustic images or participating-medium source legs.
Four focused integration checks passed using the actual old/new point wrappers
and native fixture accumulator with explicit scene, photometry and optical
response fixtures; two are new glass dispatch/stream/refusal groups and two
retain the existing opaque behavior. Two full native scene regression groups
for both source media, analyzer/Beer/index composition and cancellation are
authored but unrun. The actual new shared receiver also compiled and ran in the
environment integration probe. Full Cargo/DSR/RCH remain unavailable.

### Full Smith rough-glass finite-area illumination

The full-glass finite-area raster now adds calibrated multi-channel QE energy
and photoelectron-rate measurements, with constant or spectral analyzers,
through both ordinary and complete sequential cameras. Spectral-analyzer CIE
views complete the same six additional detector routes as point fixtures.
All channels reuse each original wavelength's full source/facet path, native
logical streams and caller-owned area workspace. Original source-side selection,
endpoint ownership, receiving index ratios and per-lane work arrays survive.

The complete detector/analyzer band is admitted before tracing, including an
entirely clipped lens packet. A spectral analyzer precedes the coated lens once;
steady optics request no group clock. Pupil and QE/photon weights are applied
once, with original quarter packet mass and zero source work for clipped lanes.
One native fixture group covers both source sides, nonconstant QE and detector
units, native replay/receipt parity, channel-count independence, original
receiver emission, coated-lens order, work exhaustion, clipping and refusal.
It was parsed and independently source/API reviewed, not executed here.


`SmithDielectricAreaLightCameraScene` borrows the admitted native area camera
and reuses its original direct-source solver, face/side selection, spectral
radiance, polarization-channel replay, Beer attenuation and work accounting.
Bound receivers may be rough dielectrics or opaque rough/brushed conductors;
each source's declared medium must be one of the receiver's original adjacent
homogeneous media. The source's actual approach side remains authoritative.
Other geometric surfaces block this straight proposal before scattering.

Both connection and continuation use the complete two-sided Smith glass law.
The glass connection owns its net local radiance index ratio and the original
receiver accumulator applies the camera-prefix ratio exactly once. Finite
sources use a direct light-only partition: only the immediately following
emitter hit represented by the same receiver, source face selection, side and
FINAL departure medium is suppressed. Arriving-medium identity is insufficient
after net transmission. Intervening events retire that partition; primary
surface emission and LTE emission remain in their original owner. No sampled
BSDF or hidden joint path density enters MIS.

At equal-index rough glass, including coated mixed continuous/delta laws, this
view retains complete full-glass continuation and endpoint weight one. It draws
neither an area source nor a connection sample there. Smooth interfaces retain
their original delta path. This does not claim continuous NEE for mixed
equal-index laws, GRIN direct legs, source images through other boundaries,
or globally optimal variance.

Separate zero-based logical streams own ordinary surfaces, opaque and glass
continuations, area points, and opaque/glass source estimates; source binding,
polarization channel and facet identity remain explicit. Seeded, passive
analyzer, resolved-lens and corresponding caller-sample entry points share the
same estimator. Cancellation and exhausted optical/source/query work refuse
the entire sample.

Verification: seven compiled integration probes passed (four new glass groups
plus three existing opaque groups), executing the actual modified shared area
estimator and new wrapper with explicit scene, geometry, emission, response and
source fixtures. They cover both-side dispatch, final-medium endpoint ownership,
face/side support, source channels, equal-index no-draw fallback, intervening
delta reset, admission, cancellation/refusal, replay and analyzed entry points.
Three native full-scene regression groups are authored but unrun. Full workspace
Cargo, resolved geometry runtime and remote DSR/RCH checks were unavailable.

### Full Smith rough-glass environment connections

The full-glass environment raster adds the same six ordinary/resolved-lens
QE/photon and spectral-analyzer CIE/QE routes. Each wavelength retains its
complete original environment receipt, texel polarization and deterministic
complementary endpoint partition. Detector channels share source/facet work;
no replacement directional proposal, stochastic-response MIS density or
wavelength resampling is introduced. Typed source failures retain their lane.

Spectral analyzers are admitted over the complete detector band and applied
before the native coated lens. Steady optics need no group clock; original
geometry and detector weights occur once. Clipped lanes retain quarter proposal
mass and no source receipt. One native fixture group covers rotated polarized
environment texels, independent QE/photon arithmetic, original seeded replay and
receipts, receiver emission, channel-count independence, full-lens weighting,
work exhaustion, clipping and band/cancellation refusal. It was parsed and
independently source/API reviewed, not executed here.


`SmithDielectricEnvironmentCameraScene` borrows the original environment camera
and extends the existing full opaque Smith transport with two-sided rough
dielectric continuation and full fixed-direction source responses. The original
source spectrum/polarization, native contact visibility, adjacent medium IDs,
surface/LTE emission, real-lens/analyzer launch and single camera-prefix index
conversion are retained. Separate Philox domains and caller callbacks index
ordinary surfaces, environment directions, opaque continuation/connections and
dielectric continuation/connections by zero-based surface/facet ordinals.

The original polarization-independent single-facet dielectric proposal is used
only as a deterministic angular guide q, with environment map density p:
source weight p²/(p²+q²), escaped-camera weight q²/(p²+q²). This is a matching
partition of the full optical integrand, not a claim that q is the complete
Smith outgoing marginal PDF. The stochastic evaluated response enters only
linearly. No hidden path density or sampled Fresnel/BSDF response enters the
weights. A later contact retires the anchor; delta chains, primary environment
and unrepresented angular laws retain their original endpoint. Emission already
encountered along the camera path is never scaled by this endpoint weight.

The connection supports rough distinct-index dielectrics. Equal-index rough
glass, including coated mixed straight-through/continuous laws, retains the
complete path-only Smith walk with endpoint weight one; no legacy single-facet
NEE is substituted. Both adjacent media are checked before an environment draw.
After a direction is chosen, its own medium is checked before connection-facet
sampling. Infinity connections retain the original homogeneous, nonabsorbing
medium requirement; GRIN/absorbing infinity legs are not approximated as clear.
A facet, visibility-work or cancellation refusal discards the complete sample.

Focused verification: seven scoped integration groups passed (three existing
opaque/ordinary regressions and four glass groups), compiling the actual whole
sampling owner, old/new wrapper modules, source/weight arithmetic and shared
native receiver modules. Scene propagation, visibility and Smith optical
responses were explicit fixtures; these checks prove dispatch, source-side
admission, deterministic complementary weights, independent streams, emission
ownership, resets, replay and complete refusal. Three full native scene groups
were authored for both-side glass/eta composition, equal-index bare/coated
path-only measures and per-side infinity admission, but were not run because
the sibling workspace dependencies and DSR/RCH were unavailable. They do not
constitute full native geometry, optics, GRIN or workspace validation.

### Native pixels with full rough-glass direct lighting

The original `CurvedCameraRaster` now exposes ordinary and resolved-sequential
pixel entry points for full Smith dielectric point, finite-area and environment
lighting. All six adapters use the original four hero wavelengths, raster/pupil
sampling, wavelength-specific rays, reciprocal lens Stokes seed and native XYZ
integration. Passive analyzers enter before transport (before the lens for a
resolved lens); lens geometry factors remain in the raster owner exactly once.
The original source/visibility/optical budgets are retained independently in each
spectral lane. No new light integrator or camera radiometric convention is added.

Point and area adapters reuse both original native lane-accounting owners and
their generic pixel receipts, retaining per-lane source radiance and work for
ordinary and sequential lenses. Environment pixels retain each complete native
`EnvironmentCameraSample`, including its endpoint partition and source/query
work. Clipped sequential wavelengths keep
zero work and their original proposal mass; environment receipts are `None`
there. Original transport errors retain their concrete type and failed lane;
any failed live lane refuses the whole pixel.

Seeded light wrappers retain their own optical/source counter domains. Four
lanes occupy the low Philox key's top two bits; original source binding/channel
addresses occupy the high key, and original fixed low-key domains remain
distinct below the lane bits. This prevents a lane address from cancelling a
source binding or thermal channel. Pixel/sample occupy the complete 64-bit
path address. Variable facet/source work cannot shift pupil, wavelength or
another lane's draws. Replay is independent of pixel scheduling.

Verification: four compiled integration probes passed using all six actual new
adapter methods, the original point and area lane-accounting owners, and the
actual environment receipt code. Explicit spectral/lens/source signature fixtures exercise
four-lane receipts and keys, exact error/lane propagation, deterministic replay,
analyzer placement, clipped-lane work and reciprocal-beam handoff. They do not
establish resolved geometry, physical lens or full native image convergence;
full workspace/DSR execution remains unavailable.

A native regression protects the previously observed point-binding and
area-channel lane alias. Two additional focused checks execute the actual
three seeded wrapper callback bodies and original Smith/Philox helpers:
they reproduce the old alias and verify distinct effective counter/key
addresses and blocks across four lanes, 64 source bindings, both area channels
and three seeds (5,520 addressed draws). All three checks and the four raster
adapter probes pass. This verifies address separation, not statistical or
physical convergence of an image.

### Matched-index rough-film continuous source connections

The full Smith dielectric fixed-direction methods now evaluate the continuous
part of rough coated interfaces with equal bulk phase index. Each non-atomic
path is partitioned at its unique last facet reflection: native continuation
samples its prefix, a forced reflection connects toward the source, and raw
native Mueller transmission weights integrate later direction-preserving
crossings. Film order, polarization frame, physical medium and complementary
statistical height are retained. The zero-transmission suffix uses exact escape
and the analytic first-height reflection term. Factored logarithmic optical
power protects lossy tails before final positive summation.

The straight-through atom remains exclusively in the original path sampler.
The whole-law continuous capability query stays false for matched media,
preserving area/environment endpoint accounting; point-fixture NEE uses the
new continuous response. Neither response nor hidden-path density is a
marginal PDF. Callback ordinals are consecutive across prefix and tail work.
One shared facet cap bounds all sampled collisions, with at most twice that
cap plus one draw, fixed storage, cancellation before each draw, and
whole-estimate refusal on exhaustion.

Five new G1/G3/G4 groups cover an analytic first reflection, independent
131072-path-per-side polarized continuous/atom integrals, transmission
reciprocity, compact/stack and frame covariance, polarization linearity,
replay, shared budgets and cancellation. All 46 optical/math/Philox probe
groups pass using native optical owners with explicit execution/vector/
container/hash fixtures. This does not establish full native scene or
resolved-lens convergence. Refs: frankensim-h7xu5.4 and frankensim-h7xu5.4.5.

### Thermal emission through complete Smith conductor reflection

`thermal_surface::SmithCoatingThermalEmission` binds the original opaque
conductor optical law to retained film/substrate LTE temperatures. At each
native Smith collision, positive layer-resolved absorption supplies separate
s/p Planck radiance in the reciprocal facet mode. The arriving Mueller
importance contracts that source before original reflection. Summing positive
sources retains weak loss even when reflected power rounds to one; no
effective albedo or cancelling one-minus-reflectance approximation is fitted.
The same walk returns its unmodified reflected beam and microscopic receipt.

Isotropic and physical-tangent anisotropic GGX share the existing geometry.
The smooth limit consumes no random draws. Exact matching currently admits
bare substrates and one original lossless or absorbing film; unmatched,
transmitting and periodic representations refuse. The unit input camera
keeps physical prefix and index conversion with the caller. Local phase
index one, original spectral support and per-layer temperatures are required.
Queries use fixed storage, no heap allocation, and bounded cancellable work;
no metric microsurface delay, temperature advance, infrared extrapolation or
bolometric cooling is supplied.

Seven new G1/G3/G4 numerical groups pass using actual coating-emission,
positive absorption, Planck, Fresnel/Mueller, Smith, fs_math and Philox owners
with explicit external execution/container/hash fixtures. They include 4608
per-realization isothermal source-plus-reflection checks, unchanged native
continuation, weak loss, independent layer temperatures, polarization,
anisotropic frame covariance and whole-sample refusals. Six original opaque
Smith groups also pass. Native `cargo check -p fs-render --features
tracer,volumes,photometry --lib --locked -j 1` passes with one compiler thread
and debug info disabled; full native scene-test execution remains separate.
Refs: frankensim-h7xu5.4 and frankensim-material-reality-ui2if.6.4.

### Complete Smith surfaces inside physical participating transport

`participating::microsurface::SmithDielectricRayleighScene` and
`camera::participating::SmithDielectricRayleighCameraScene` compose original
full Smith metal/glass with the native molecular/spectral-particle tracker.
Physical free flights, rejected null candidates, local Mueller phase laws,
finite-cell absorption/LTE, GRIN propagation and actual volume-event limits
retain their original owners. The completed microscopic escape selects the
physical medium for the next volume flight. Extended departures have no
marginal direction PDF and refresh the thermal Mueller prefix.

Separate surface, physical-flight, opaque-facet and dielectric-facet streams
retain zero-based logical addresses; variable work in one domain cannot
advance another. Seeded zero-scattering transport preserves the existing
Smith path stream. Forward source power and reciprocal ordinary/analyzed/
resolved-lens cameras are exposed without inventing statistical-height
distance, absorption or delay.

The native raster adds four full-Smith participating entry points: visible
and multichannel QE observations, each with ordinary and sequential lenses.
Original CIE, QE/photon proposals, pupil geometry and lens Stokes seeds are
reused. Volume events and tracking draws are retained per wavelength; clipped
lanes keep zero work and their original probability mass. Every live-path
failure refuses the complete pixel.

Native library checking passes with `tracer,volumes,photometry`, one compiler
thread and no debug info. Four focused native integration groups are authored
for zero-scattering parity, real Rayleigh-to-multiple-facet glass paths, both
physical media, independent Beer/eta/LTE composition, source-power parity,
detector measure, resolved lenses, and whole-sample limits/cancellation.
Their native test binary is still building at this commit; no executed scene
test pass is claimed here. This is path-only source observation, without
participating NEE or new caustic image coverage. Ref: frankensim-h7xu5.4.
