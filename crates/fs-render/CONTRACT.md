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
  inputs convert knotwise as `alpha = 4*pi*k/lambda_m`, retaining the stated
  linear-alpha approximation and weakly absorbing dielectric interface model.
  No RGB fit, opacity floor, or extrapolation is introduced. Thin-path absorbed
  fractions use `expm1`, and material identities retain every admitted knot.
  Explicit `SpectralAbsorption::try_new_in_band` and Beer-Lambert band
  constructors retain 2–16 original positive vacuum-wavelength knots on their
  actual inclusive support, including infrared data. The legacy constructors
  still require 380–780 nm and retain identical canonical parameters. The
  imaginary-index band constructor converts intensity absorption with the same
  vacuum-wavelength law before linear-alpha interpolation. Clear and constant
  extinction are explicit achromatic declarations; lifted RGB retains visible
  support and a source table is never extrapolated. `DielectricGlass` exposes
  the checked intersection of its phase-index and absorption bands, refusing
  disjoint sources while retaining a shared endpoint for monochromatic use.
  No extra allocation or per-query spectral resampling is introduced. G1/G3
  regressions cover original infrared knots, path composition, intensity versus
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
  in the implementation environment. Polarized scene sampling and direct connections explicitly
  refuse this new finish until they carry its body-frame binding. This remains
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
  retain their separate proximity-based offset. This is geometric-optics
  transport, not coherent interference across a geometrically thin body or a
  general certified minimum separation. Analytic regressions cover submicron
  absorption, distant camera origins, oblique Snell/Fresnel power, metre
  scaling and rigid placement.
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
  including interior extrema. The allocation-free augmented Dormand--Prince
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
  scalar camera law `F_mueller*G1(wi)`, selected before forming the forward-power
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
  errors. This covers reciprocal isotropic stationary media and intensity-only
  endpoints, not general polarized adjoints. Terminal transport weight describes only the final continuation,
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
  occur once. Rough glass admits reflection and transmission from either actual
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
  native geometry and its exact `SurfaceEmitter`, directly from the surface
  scene or an admitted camera emitter. Mesh sources may include all
  original faces or an explicit strictly increasing face selection through
  `try_new_faces` / `from_camera_faces`; disconnected, noncoplanar and closed
  meshes retain original face winding and per-side radiance. `AreaEmitterLimits`
  bounds the whole mesh, preparation work and each hit-membership query.
  Admission reuses the original BVH and exact coplanarity/orientation
  predicates to reject positive-area coplanar overlap, including overlap with
  excluded faces. Shared edges and vertices remain zero-area contacts; no
  mesh-union or multiplicity correction is inferred. The triangle-area CDF
  retains original geometry, frame and calibrated spectral radiance. Samples
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

## Material thermo-optic refraction

With `material-state`, `MaterialThermoOpticField` binds the original reference
glass and `fs-material::thermo_optic::SchottThermoOpticLaw` to one immutable
thermal snapshot. Its usable wavelengths are the intersection of reference
phase/absorption and thermo-optic support. Preparing a wavelength evaluates
reference dispersion and the spectral material coefficients once; optical
stages still query the actual spatial temperature and its analytic derivatives.

The original gradient-index and ray-bundle owners consume
`grad(n) = n_T grad(T)` and
`Hess(n) = n_T Hess(T) + n_TT grad(T) grad(T)^T`, including mixed terms.
Phase and available group optical distance retain the original material law;
absorption remains the explicitly declared reference law. The same bound
source supplies the native temperature trait. No finite-difference field,
temperature-dependent absorption, moving state, expansion or stress response
is invented. Each material/contact discontinuity requires a separate region
and the native boundary event, not a smooth integration across the jump.

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
v3 and must be reproduced in debug and release before a justified v4 re-freeze;
neither its digest nor its evidence pin has been replaced without that proof.
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

`PhysicalCamera::with_in_focus_diffraction(FocusedDiffraction)` opts the
cinematic spectral tracer into scalar, paraxial, incoherent pupil imaging.
It requires a nonzero circular, annular or regular 3–64-blade aperture,
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

With `tracer`, `ThinLensElement` and `DispersiveThinLens` explicitly enable
longitudinal chromatic aberration in explicitly selected geometric,
finite-defocus and layered camera paths. A prescription contains one to eight coincident ideal thin
elements, each retaining its exact native `DielectricIor` (Cauchy, Sellmeier or
tabulated phase index) and signed front/back surface curvatures in inverse
metres. The exterior index is exactly one and wavelengths are vacuum nm within
the intersection of the original material laws’ source bands, including
explicit infrared support; group index is never substituted for phase index. No IOR fit, RGB conversion, external optics library or new pupil raster
is introduced. The immutable bounded prescription is shared across evaluated
camera states; per-ray material evaluation and pupil propagation allocate
nothing.

For reference wavelength `lambda_ref`, element power adds at the common pupil:

`P(lambda) = sum_j (n_j(lambda)-1)*(C_front,j-C_back,j)`.

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
Their actual launch phase index must equal the prescription’s index-one
exterior; a mismatch refuses the complete pixel. Source bands/gaps,
`lambda/opening_width <= 0.01`, and the full pupil’s outgoing paraxial slope
remain admission conditions. This mode omits diffraction; shared-ray entry
points refuse it, and curved raster still refuses wave diffraction without
the original source-depth and clear-pupil proof. Achromatic cameras retain
their existing ray arithmetic and Philox streams.

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
matching dielectric-interior ray ownership. The model does not claim polarized
transport, particle interactions, nonspherical particles, measured material
dispersion, or an inferred aerosol distribution. A constant-index visible spectrum
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
runs the original native estimator once and partitions its retained cell heat
by the same fixed spectral absorption shares used for beam deposition.
`ConstituentScatteringDeposition` retains the complete unchanged native report,
including ballistic and diffuse escape, unresolved radiation, correlated
whole-history totals, sample counts, work counts, seed and source configuration.
The partition introduces no retracing, additional random draws or thermal debit.

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
