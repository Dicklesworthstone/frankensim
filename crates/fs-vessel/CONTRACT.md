# fs-vessel CONTRACT

Flagship 3 (plan §15.3, bead mye.4): the LAMINAR-POUR VESSEL, smoke
tier. Design a carafe lip so the pouring film stays laminar across a
FAMILY of fluids — stability objective, free-surface validation,
CVaR robustification, e-raced screening, and a render that is the
physics.

## Physical SI free-surface film studies

`stability::laminar::LaminarFilmStudy` is the physical alternative to the
historical U-shaped channel proxy. Each section declares its actual wetted
width and effective acceleration tangent to and toward the supporting wall.
An operating point supplies total volume flow and the actual constant liquid
density, kinematic viscosity and surface tension. These are SI quantities;
no arbitrary factor converts them into the two-wall transition benchmark.

The original `fs_lbm::film::NusseltFilm` owner solves the no-slip/shear-free
base profile and the isothermal Benney long-wave approximation. Continuity
holds section by section: `q=Q/(width*scale)`, `Re=q/nu`. The caller supplies
physical wavenumbers and explicit `kh` and `Re*kh` budgets. Every section and
mode must be admitted; failure returns its exact section/mode, never a clipped
Reynolds number, skipped mode, or invented stable score. The report identifies
the worst sampled growth in 1/s and the observed nondimensional ranges.
Evaluation is allocation-free and linear in section count times mode count.

`LaminarFilmTransitStudy` additionally transports one inlet frequency through a
finite, slowly varying attached-film path. Each section carries an actual length.
Continuity fixes the Nusselt kinematic-wave speed at `c=dq/dh=3*q/h`; the local
wavenumber is `omega/c`, and signed log amplification is the sum of `s*L/c`.
The relative film-height/flow disturbance retains no extra geometric factor,
because `width*c*h=3*Q`. Liquid residence time is reported separately and is
three times the wave transit time. The report retains the largest accumulated
amplification along the path, so downstream decay does not hide prior growth.

This is the leading weak-growth temporal-to-spatial conversion of the same
Benney dispersion. Every section must pass both native long-wave limits and an
explicit `|s/omega|` budget. It does not solve the full spatial eigenproblem,
reflection at abrupt geometry changes, nonlinear saturation or transition to
turbulence. Geometry must vary slowly enough for the local parallel assumption;
piecewise sections represent quadrature of that path, not resolved sharp joints.
Evaluation needs constant storage and work linear in section count. The G1/G3
transit tests use independent dimensional coefficients, fixed-frequency
transport, path ordering, length scaling and homogeneous section partitioning.

`LaminarFilmTransitBandStudy` maximizes amplification over a continuous inlet
frequency interval and the entire finite path. The original film owner exposes
its dimensional dispersion `s(k)=a*k^2-b*k^4`. Each path prefix therefore has
gain `A*omega^2-B*omega^4`, with `A=sum(a*L/c^3)` and `B=sum(b*L/c^5)>=0`.
Endpoints and the positive stationary point exhaust its possible maxima.
Section ends exhaust the spatial maxima because gain is linear in distance
within each locally constant section. Evaluation uses normalized frequency,
compensated prefix sums, constant storage and O(section count) work.

The report distinguishes the worst outlet frequency from the frequency with
the largest earlier amplification. An inlet peak has value zero, section count
zero and `None` for its frequency: all inlet frequencies tie there. A zero
frequency endpoint is the neutral mass-conservation limit; a point band
recovers the existing single-frequency model. Every local wavelength in the
band must satisfy native long-wave admission. The weak-growth budget also
checks the interior extrema of `|s/omega|`, a cubic in frequency, so neutral
band endpoints cannot hide an invalid interior conversion. The independent
regressions in `tests/laminar_transit_band.rs` cover distinct interior optima,
upstream decay, both full-band admission conditions and section subdivision.

`robust::robustify_laminar` uses the same physical evaluator over an authored
width-candidate set and fluid band, then invokes canonical empirical CVaR.
It returns both candidates' band CVaR and worst sampled growth. Equal nominal
and robust choices, ties and boundary optima are legitimate outcomes; the
algorithm does not force an interior optimum or a strict robustness benefit.
A domain failure at any candidate/band point refuses the complete study.

`race::screen_laminar_lips` sends native growth rates to the existing public
race owner. Its fixed demonstration jitter has units 1/s here and remains
synthetic; it is not measured fluid uncertainty. Candidate refusals produce
no race verdict. Historical smoke APIs, browser fixtures and frozen smoke
goldens retain their original explicitly documented proxy.

Scope: constant-property Newtonian, locally planar, slowly varying,
fully developed, laterally wide films with no slip and negligible gas shear.
Section widths alone do not establish lateral-wall, curvature or entrance
validity. The long-wave budgets are declared approximation restrictions,
not error certificates. Positive growth predicts a surface-wave mechanism,
not turbulence; negative sampled growth is not proof over all wavelengths.
This route does not replace or calibrate the transient LBM pour, contact-line
physics, or measured-fluid laws.

`tests/laminar.rs` covers the independent dimensional dispersion relation,
conserved flux/viscosity scaling, finite-path disturbance transport, mode-domain
refusals, physical CVaR candidate selection and the original race integration.

## Purpose and layer

Layer **L6 (HELM)**. Composes battle-tested lower crates end-to-end:
fs-cheb (profile + Orr–Sommerfeld), fs-lbm (free-surface pour,
rheology), fs-race/fs-exec (e-raced screening), fs-robust (canonical
empirical CVaR), and fs-render `volumes` (exact slab transmittance
deliverable). No new numerics live here; the crate's claim is the
COMPOSITION with certificates at every joint.

## Public types and semantics

- `stability::VesselProfile { radius: Cheb1, lip_width: f64 }` — a
  vessel of revolution: Chebyshev radius r(z) on z ∈ [0, 1] (base →
  lip) plus a scalar lip-channel width. `carafe(lip_width)` builds the
  smooth default profile.
- `VesselProfile::film_reynolds(rate, viscosity, stations)` — the
  quasi-steady film Reynolds proxy along the lip run z ∈ [0.7, 1.0].
  **U-shaped in lip width by design**: a thickness branch
  (film/viscosity, wide lips) plus a velocity branch (2.5/film,
  viscosity-blind, narrow lips), scaled 1000× so the design range
  straddles the plane-Poiseuille transition. A SMOKE PROXY, documented
  as such — not a calibrated film model.
- `stability::growth_objective(profile, rate, visc, stations, modes)`
  — min-max certified modal growth: max over stations and the first
  `modes` Orr–Sommerfeld modes (α = 1.02056, n = 32 collocation) of
  the real growth rate. Negative = every mode decays (laminar proxy).
- `pour::PourRig` — the weir-tank fixture (48×28 lattice, lip column,
  rest-to-rest rigid tilt schedule, and explicit nonnegative lattice surface
  tension). Defaults use g0 = 6e-4, sigma = 0.002 and tilt 0.7 rad over
  900 steps. Set sigma explicitly to zero for a no-capillarity comparison.
  SI conversion is sigma_lu = sigma_SI dt² / (rho_SI dx³); the fixture
  default does not imply a sourced water or other fluid calibration.
- `pour::run_pour(rig, contact, law) -> Result<PourOutcome, PourError>` — one physical rigid-frame pour;
  outcome carries `mass_drift` (worst relative ledger drift),
  `fragments` (2-D connected-component count), `poured_mass` (mass right of
  the weir), `dribble_cells` (wet outer-lip cells), and `mass_field`
  (the same bytes the render binds).
- `pour::advance_pour_liquid` owns the checked constitutive advance shared by
  native and browser pours. Every step uses the original `update_tau` and
  refuses before collision/streaming/VOF if any liquid or interface cell
  would need a viscosity below its numerical floor or above its ceiling.
  The native counts do not include gas or walls. The error reports the
  zero-based step and both counts;
  the caller must rescale the lattice instead of accepting a different fluid.
  This adds no population pass or per-step allocation. A refused trial retains
  its prepared tau buffer but has not advanced populations, topology or mass;
  both consumers discard it and publish no successful pour or captured frames.
- `pour::run_pour_in_frame(rig, contact, law, pivot)` selects an explicit
  lattice pivot; the default is the centre of the halfway-bounce-back floor.
  A quintic clockwise tilt starts and ends with zero angular velocity and
  acceleration. Its midpoint angle/rates drive gravity, centrifugal, Euler
  and Coriolis terms through `fs_lbm::core2::frame::RigidFrame2` on the
  original populations. Force replacement preserves physical velocity and
  stress before collision. The native viscosity, capillary, wetting and mass
  operators remain active. Support motion is prescribed; no free vessel,
  gas inertia, 3-D geometry or calibrated jet-breakup claim is added.
- `pour::run_pour_rotating_gravity` explicitly retains the historical linear
  gravity-only smoke model. It omits frame inertia and is not a physically
  rotating vessel. The existing cross-ISA smoke golden selects this legacy
  model explicitly; default and MID pours select the rigid-frame model. It
  returns the same `Result` and also refuses clipped rheology.
- `pour::render_pour(outcome, nx, ny, res)` — fs-render exact Beer--Lambert
  slab transmittance bound ZERO-COPY to `outcome.mass_field`, with explicit
  image/field/work admission. Lattice mass is displayed as absorption per
  lattice length through an `ny`-thick slab. This is an absorption
  visualization, not sourced liquid optics or refracting 3-D geometry.
- `robust::robustify(beta) -> RobustReport` — nominal (band-center)
  vs CVaR_β lip optimization over the fluid band (rates ×
  viscosities), both evaluated on the OFF-NOMINAL corners.
  `robust::empirical_cvar` and `robust::cvar` directly re-export the canonical
  `fs-robust` report and scalar surfaces, respectively. The report retains
  deterministic VaR/minimizer and fractional-boundary metadata;
  empty/non-finite losses and beta outside `(0,1)` are structured refusals.
  `robustify` generates its own fixed, finite fluid band and
  treats a canonical refusal there as an internal programmer-contract
  defect.
- `race::screen_lips(lips, seed) -> Result<LipScreenReport, ScreenError>`
  — the vessel's PUBLIC e-raced candidate screen: score each lip with
  `growth_objective` at the nominal fluid (`SCREEN_RATE = 1.0`,
  `SCREEN_VISCOSITY = 1.0`, `SCREEN_STATIONS = 3`, `SCREEN_MODES = 3`),
  then race the table through `fs-race` under the vessel's DECLARED
  convention. `race::race_base_losses(base, seed)` is the same race over
  an EXPLICIT loss table, so an auditor outside this crate can drive the
  vessel's convention over a table it did not generate.
- **The vessel's declared racing convention**, owned here and nowhere
  else: losses handed to the race are `SCREEN_SCALE x (base + jitter)`
  with `SCREEN_SCALE = 200.0`; the jitter is a hashed counter of
  `(candidate, round, seed)` with total width
  `SCREEN_JITTER_WIDTH = 1e-4` in unscaled objective units; and the
  declared paired-loss support is DATA-DERIVED,
  `race::declared_span(base) = SCREEN_SCALE x (fixture spread +
  SCREEN_JITTER_WIDTH)`. The equality of the span's slack term and the
  jitter width is the SOUNDNESS of the declaration — a paired difference
  of jitters cannot exceed `SCREEN_JITTER_WIDTH` — and drifting the two
  apart is the failure mode the fs-flagship-e2e fe2e-006 audit exists to
  catch. This convention was inlined in `tests/battery.rs` until bead
  `frankensim-extreal-program-f85xj.2.31`; a convention that lives only
  in a test cannot be audited by anyone else, which is why it moved.
- `race::LipScreenReport { winner, eliminated, evaluations_used,
  fixed_n_equivalent, losses, declared_span }` — the screen's evidence
  row, carrying the declared span it actually raced under.
- `race::ScreenError` — structured refusals: `TooFewCandidates`,
  `NonFiniteLoss { candidate, value_bits }`, `InvalidSpan { span_bits }`,
  and `Race(fs_race::RaceError)`. A loss table that cannot support a
  declared span produces no verdict; it never produces a forged one.

## Invariants

1. The mass ledger is STRICT on every shipped pour: relative drift
   < 1e-10 at every step (gated by vsl-002/003/004).
2. The render reads the simulation's own mass buffer (zero-copy
   borrow) and replays BITWISE (vsl-006).
3. The growth objective is deterministic: same profile and fluid →
   the same eigensolve → the same objective, byte-for-byte.
4. The contact-line bracket (Neutral vs Wetting) is a REPORTED
   first-class output; it is never folded into a single pretended
   number.
5. `film_reynolds` floors at Re = 50 (a defensive guard; the shipped
   scale keeps the design range at Re ≈ 3200–6000 nominal, so the
   floor never binds).
6. The e-raced screen's declared span is derived from the loss table it
   races, and its slack term equals the jitter width, so every paired
   difference the race observes lies inside the declared support. If it
   cannot (non-finite loss, non-finite derived span), the screen returns
   a `ScreenError` and no verdict.

## Error model

Direct empirical-CVaR calls preserve `fs-robust`'s typed `RobustError`
refusals. The e-raced screen returns `race::ScreenError` (never a
fabricated verdict) when the loss table cannot support a declared span
or when the shared race core refuses. Fixture-scale orchestration errors
panic (`expect`) — including
an impossible rejection of the internally generated fluid-band losses,
the eigensolver failing to converge, `Cheb1::build` failing on the smooth
carafe, or render dimensions disagreeing with the outcome buffer. The
flagship has no recovery story for those internally generated smoke-tier
invariants and converts their failure to a loud stop.

All three public pour runners now return `PourError` for a zero-step schedule,
unrepresentable local viscosity, or the physical frame owner's refusal. A
`RheologyClipped` result does not mean the fluid itself is invalid: its requested
viscosity is outside the native numerical range at the selected scale. A
shear-dependent law is checked at every evolving state, not only at rest.
Malformed constitutive parameters and the existing native geometry/capillarity
programmer contracts still panic; no broad input-recovery API is claimed.

## Determinism class

Fully deterministic. The tilt schedule is closed-form, fs-lbm
free-surface stepping is deterministic, the race jitter is a hashed
counter (no RNG state), and the exact absorption render uses the original
deterministic Beer--Lambert owner without random sampling. Bitwise replay is
gated (vsl-006); image values intentionally change from the previous 24-sample
collision estimates to their analytic values.

## Cancellation behavior

None. All entry points run to completion at fixture scale (seconds).
The e-raced screen consumes an `fs_exec::KillRegistry` and records
kills, but fs-vessel itself never blocks.

## Unsafe boundary

`#![forbid(unsafe_code)]` — none.

## Feature flags

None. The fs-render dependency enables `volumes` in addition to fs-render's
default certified chart backends; nothing in fs-vessel itself is gated.

## Conformance tests

`tests/battery.rs`, verdict-JSON rows (all gates strict unless marked
REPORTED):

- **vsl-001** stability objective: moderate lip (0.6) grows slower
  than wide lip (2.4) through the thickness flank, AND the Orszag
  bracket σ(5000) < 0 < σ(6000) re-gated through this crate's call
  path.
- **vsl-002** the capillary-enabled pour crosses the lip while retaining
  a strict mass ledger under the rigid-frame tilt schedule.
- **vsl-003** contact bracket: Neutral vs Wetting must change the solved
  mass field while both mass ledgers remain strict. A companion native
  pour test requires a nonzero response to turning surface tension off.
  These are numerical sensitivity checks, not calibrated contact-line
  or quantitative cylindrical-jet breakup validation.
- **vsl-004** Carreau band under the rigid-frame schedule: three fluids
  pour [2.04, 5.94, 6.53], spread 4.49; all ledgers remain strict.
- **vsl-005** robustification: robust lip 0.92 stays stable off-band
  (−0.00044) while the nominal lip 1.11 goes UNSTABLE off-band
  (+0.00023) — the flagship claim, gated. The e-raced screen, driven
  through the public `race::screen_lips` wrapper, recovers the
  deterministic argmin under jitter and eliminates all 4 dominated lips
  in 136 evals vs fixed-N 2000 under declared span 1.49318. Invalid
  direct CVaR inputs return exact structured refusals rather than
  panicking or reporting fake risk.
- **vsl-005b** (`vsl_005_race_wrapper_owns_the_declared_convention`) the
  extracted wrapper IS the convention the battery used to inline: the
  declared span matches `200 x (spread + 1e-4)` bit-for-bit, the base
  loss table matches `screening_losses` bit-for-bit,
  `race_base_losses` over that table reproduces `screen_lips` exactly,
  and degenerate tables (one candidate, a NaN loss, an overflowing
  spread) return structured refusals instead of a forged span.
- **vsl-006** deliverable: bitwise render replay, transmittance range
  1.000.

Measured design rejections (kept for the record): a clamped
single-branch Reynolds proxy pinned every design at Re = 50 (constant
objective); a monotone unclamped branch pinned both optimizers at the
lip-grid floor (no robustification story); g0 = 1.2e-4 never overtopped
the weir in 700 steps (poured mass 0.00); unscaled growth losses
(~1e-3 gaps) starved the PairwiseRace betting process (0 eliminations
in 400 rounds) — losses are now scaled by 200 for power, while the exact
fixture-wide paired difference plus jitter is declared separately as the
checked `LossSpan`. A support breach produces no race verdict.

## No-claim boundaries

- **Level-set lip topology change** (fs-topols wiring): the scalar
  lip-channel width is the smoke knob; reshaping the lip as a level
  set is the RECORDED SUCCESSOR, not claimed.
- **Cumulant collision**: pours run BGK-family collision via fs-lbm;
  no cumulant-collision claim.
- **Recorded-fluid calibration**: Carreau parameters are plausible
  bands, NOT calibrated against measured liquids.
- **Appendix C study-program runner**: no ledger artifacts, no
  experiment orchestration — the battery rows are the evidence here.
- **Fine-mesh validation**: 48×28 fixture lattice only; no
  grid-convergence claim for the pour.
- **Film model fidelity**: `film_reynolds` is a documented U-shaped
  proxy; no claim that its Reynolds numbers match a resolved thin-film
  solution.

## Rigid-frame momentum checks

`PourRig::prepare_step` is the shared native/browser entry for one scheduled
force update. It admits the zero-based step before using the same quintic
midpoint kinematics and original D2Q9 rigid-frame owner. Preparing a step
changes no physical momentum; the original collision supplies the held impulse.
An empty or exhausted schedule refuses before changing the liquid. Consumers
retain their own geometry, capture policy and subsequent native flow step.

The original native flow tests retain mass, wetting/capillary response and
Carreau sensitivity with the physical schedule. New tests verify its zero
endpoint spin/acceleration and that moving the pivot changes the native flow
even at exactly zero gravity. That zero-gravity response comes from frame
inertia, not a rotated body-force proxy. A native comparison also verifies the
explicit historical smoke path is bitwise identical to its prior implementation.
The core's independent momentum/stress and pure-Coriolis tests cover the force
operator; this consumer does not claim global energy exactness for LBM/VOF.

The G0/G3 pour regressions additionally reject both too-small and too-large
Newtonian viscosities on every public motion path, reject a Carreau law only
when the evolving flow reaches the numerical floor, and prove that refusal
precedes changes to populations, interface mass and conversion state. The
existing Newtonian/contact/Carreau pour tests require successful admission.
The browser uses the same checked advance and drops a refused trial's frames.

## Continuous-band physical optimization

`LaminarFilmBandStudy` asks the original Nusselt/Benney owner for the maximum
over every wavenumber in an authored SI interval at each actual section.
Its evaluation is allocation-free and O(sections), independent of wavelength
sampling density. It reports the worst section/wavenumber and full-band
nondimensional extrema. An unsupported band tail refuses the candidate even
when its maximizing mode would be admissible.

`LaminarGrowthObjective` lets the existing `robustify_laminar` and
`screen_laminar_lips` consume either sampled or continuous-band physical studies.
Risk still uses canonical empirical CVaR; racing retains its declared synthetic
jitter. Robust reports now retain every evaluated candidate's nominal growth,
band CVaR and band maximum, so product consumers can plot actual objectives
without repeating the candidate/fluid sweep. Existing sampled studies keep
their previous semantics. Continuous-band maxima do not establish global,
finite-length, nonlinear or turbulent stability.


## Monochromatic physical pour opacity

PourOpticalView renders an accepted mass field or PourOutcome through the
shared fs-render density-opacity integrator. The view explicitly supplies
vacuum wavelength, mass absorption in m²/kg, kg/m³ per lattice mass unit,
in-plane lattice spacing and an independent viewing depth in metres. Field
storage is borrowed; no converted extinction image or new simulation is
needed. Every pixel now integrates its complete rectangular footprint by
averaging original-column transmission, so a thin stream cannot disappear
because its pixel centre lands in gas. Doubling depth squares each column's
direct transmittance before filtering; a mixed pixel's mean transmission is
not squared. Refining the representation or changing density units without
changing the physical field preserves the filtered result. The fs-wasm render_vessel_transmittance export exposes
the same computation for captured browser frames with fixed work limits.

An interface cell's liquid mass is treated as homogeneous along the viewing
ray. This is a declared cell-average slab model: it does not reconstruct
free-surface geometry or subcell coverage, wall occlusion, refraction, scattering,
or a liquid's measured spectrum. The coefficient must be supplied at the
declared wavelength and material state. The older render_pour visualization
retains its explicitly authored lattice-opacity convention.

The original cell-work budget now covers every overlapping source column.
Native/browser regressions specify subpixel opaque-stream coverage, analytic
mixed-column transmission, lateral refinement and one-short work refusal.
The shared renderer's exact footprint kernel passed native tests; these
downstream integration tests were added after the execution environment became
unavailable and remain unexecuted.

### Spectral pour observation

`PourSpectralOpticalView` borrows 1–16 explicitly supplied increasing vacuum
wavelengths and their mass absorption coefficients at the observed material
state. It uses the same accepted mass field, SI density conversion, lattice
spacing and independently supplied slab depth as `PourOpticalView`. Its full
image retains wavelength labels, transmission, thin opacity and effective log
attenuation. `render_transmittance` and the `PourOutcome` convenience method
write the direct pixel-major, wavelength-fastest buffer without constructing
the full-statistics image. Every wavelength is attenuated before exact spatial
filtering, and original geometry/density traversal is shared across wavelengths.
No spectral interpolation, bin widths, illumination or detector response is
inferred from the supplied samples.

The existing fs-wasm flagship bindings now expose
`render_vessel_spectral_transmittance` with the same original mass bytes and SI
parameters. The browser retains its 1,048,576 pixel/grid/cell-work ceilings and
adds a 1,048,576 column/wavelength-evaluation ceiling across the complete
spectrum. Invalid inputs or exhausted budgets return an empty vector before
publishing partial output. This remains bulk absorption without refraction,
wall occlusion, scattering, emission or an inferred subcell free surface.
Four native/browser regressions specify spectral filament coverage, one-band
parity, original-output forwarding, sample layout and total-work refusal.
The shared renderer passed its 117-test targeted native suite; these four
consumer tests received parser/API/physics review but remain unexecuted because
the actual fs-vessel/fs-wasm dependency closure is unavailable here.

### Indexed material spectra in native and browser pours

`PourIndexedSpectralOpticalView` adds one borrowed `u32` material/state index per
original mass cell and a material-major, wavelength-fastest mass-absorption
palette on common supplied wavelengths. The native full-image and direct
transmission methods, plus `PourOutcome::render_indexed_spectral_transmittance`,
reuse the same SI geometry validation, original mass bytes and independently
declared slab depth. Each cell uses its own coefficient before exact footprint
filtering. The full spectral image also supports explicit-illumination visible
observation through its original renderer API. No material identity, mixture
law, spectral data or unresolved interface is inferred from density.

`render_vessel_indexed_spectral_transmittance` and its wasm export expose this
path directly to browser consumers with the same pixel-major, wavelength-fastest
layout. Browser admission allows at most 4,096 material rows (65,536 coefficients
at 16 wavelengths), 1,048,576 pixels/grid cells/cell visits, and 1,048,576 total
cell-by-wavelength evaluations. The complete palette and map are validated,
including unused rows and vacuum cells, before output allocation. Invalid input
or exhausted budgets returns an empty vector. Four native/browser regression
fixtures cover distinct material colors and footprint coverage, original-outcome
forwarding, palette/index refusal, uniform parity and total work limits. They
have pinned parser and independent API/physics review; native and WASM runtime
execution remains unavailable in this environment.

## Cylindrical free-jet capillary screen

`jet::FreeJetStudy<J>` composes an original `fs_lbm::jet::JetDispersion`
owner with an authored cylindrical radius, volume flow, flight length,
continuous physical wavenumber band and disturbance amplitude/threshold.
Continuity fixes U=Q/(pi R²). Native temporal growth is observed over L/U;
the result retains the fastest native mode, natural-log amplification,
log amplitude/radius, and the linear time and advected distance to the
declared threshold. Neutral or oscillatory bands below threshold return
no threshold time. Large scores remain logarithmic, avoiding exponential
overflow. Evaluation is allocation-free and independent of wavelength
sampling density.

The inviscid `CapillaryJet` includes azimuthal curvature and its stable
capillary-wave branch. The Newtonian `ViscousSlenderJet` includes extensional
viscous damping in the linearized one-dimensional cylindrical model and
requires an explicit maximum kR in (0,1]. It refuses a band outside that
budget. Its zero-viscosity limit is the slender approximation, not the exact
inviscid cylindrical dispersion. Both models maximize the continuous admitted
band with at most three mode evaluations; no sampling or runtime dispatch is
introduced.

`MaterialCapillaryJet` and `MaterialViscousSlenderJet` implement the same
trait and can be carried directly by a study. They borrow the exact resolved
liquid source, preserving its density, viscosity, tension, query coordinates,
property selections and receipts. The inviscid binding requires a declared
Ohnesorge allowance; the viscous binding consumes the source viscosity
directly. Neither silently freezes an evolving non-Newtonian rheology source.

This screen is a translating uniform cylinder with an authored normal-mode
amplitude; it does not infer a jet radius or noise spectrum from the 2-D pour,
replace its mass field, solve inlet-forced spatial instability, or resolve
gravity stretching, gas inertia/viscosity, surfactants or nonlinear pinch-off.
A threshold crossing is a linear extrapolation, not a simulated breakup or
measured droplet size. The caller must justify the selected model and
small-disturbance approximation. A long-wave budget declares admitted
wavelengths; it does not certify their approximation error.

The module's G0/G1/G3 tests cover the independently evaluated cylindrical
dispersion reaching the final score, conserved volume flux, flow/length
sensitivity, stable and neutral bands, zero transit, logarithmic large scores
and structured refusal. The existing pour connectivity result keeps its
original 2-D meaning.

## Source-backed SI free-surface pours

`pour::material::MaterialPour` binds an original resolved Newtonian liquid
to the same rotating-frame weir solver used by `PourRig`. The shared native
free-surface planner converts the source density, viscosity and surface tension
at the requested cell size and timestep. The original source and its property
selections remain borrowed by the accepted outcome; no material is reselected
or numerical viscosity substituted.

The authored physical scenario supplies gravity, pivot, tilt duration, grid
geometry, contact-model bracket, allocation limit and velocity/density limits.
Preparation enforces the native Mach, capillary and relaxation admissions.
The observer checks actual vessel-frame wet-cell velocities, including the Guo
half-force correction, and density excursions at initialization and every
completed step. Cancellation or a failed envelope returns no partial outcome.
The legacy lattice fixture uses a no-op observer on that same integration path.

The accepted report exposes physical duration, vessel-frame speed and mass per metre of
out-of-plane depth. The mass conversion is rho * dx²: this two-dimensional
simulation does not infer a vessel depth or report kilograms. The original
mass field remains the rendering field; physical mass values are available
through a lazy iterator without a second field allocation.

This binding holds one Newtonian source at its resolved conditions and refuses
sources requiring evolving rheology. It does not solve thermal or composition
evolution, certify laminar transition or substep extrema, calibrate the contact
model, or turn two-dimensional connectivity into three-dimensional breakup.
The focused regressions compare the exact native trajectory and render field,
source-temperature response, physical unit conversion, cancellation and
declared speed/density refusal.

## Nonlinear cylindrical free flight

`jet::nonlinear::NonlinearFreeJet<J>` composes the original periodic nonlinear
slender-jet owner or its material-bound implementation. Initial radii are at
cell centers and laboratory axial velocities at faces. The study derives
transport speed from initial axial momentum divided by liquid mass; it accepts
no contradictory second flow rate and applies no hidden reference-frame shift.

The authored flight distance and that conserved mean speed define a nominal
residence target. Each call advances the native PDE for an explicitly bounded
time chunk, capped by remaining residence time. Reports distinguish ongoing
flight, completion, and the original neck/slope/resolution stop. A domain stop
retains the last admitted native state and terminates that flight; a subsequent
call observes it without inventing further evolution. Numerical/cancellation/
work-budget errors retain the pre-call state. Native volume, momentum,
kinetic/surface energy and separate physical/numerical losses reach the result.

The original owner and admitted initial profiles remain available. A material
binding therefore preserves the same immutable source through all accepted
chunks. This temporal-to-distance mapping uses bulk speed; it does not track
each parcel's residence, impose nozzle/end conditions, resolve gravity-induced
stretching or continue through pinch-off. A neck limit is a model stop.

Four focused native G0/G1/G3 tests exercise nonlinear shape/velocity changes,
volume and momentum conservation, velocity-dependent residence and necking,
early neck-limit publication, zero flight, budget rollback and successful retry.

## Source-driven viscosity during a physical pour

`MaterialRheologyPour` adds generalized-Newtonian source viscosity to the
same admitted SI weir pour. An explicit `IsothermalRheologyLaw` retains the
original material card, selected claims, fixed temperature and context.
Each incoming wet-cell Guo-corrected stress closes physical shear against
that source through the shared native bounded root. The resulting per-cell
tau drives the existing BGK collision and conservative VOF integration.

Geometry, capillary admission, prescribed rotating-frame motion, endpoint
speed/density limits and SI mass-per-depth reporting are shared with
`MaterialPour`. The simulation mass buffer remains the rendering field.
A read-only preparation can be retried; cancellation or source/native refusal
returns no partial pour. Frozen `MaterialPour` retains its existing rejection
of evolving sources.

The accepted result exposes `pour()` and `last_collision()`. Complete
material receipts and physical rates describe the last collision input,
before streaming and conversion, not a second source sweep of the final field.
Rates are objective physical shear; the speed envelope is in the vessel frame.
Pressure-dependent rheology, temperature evolution, viscoelastic memory and
yield history remain outside this explicit fixed-temperature shear owner.

Focused synthetic regressions specify exact constant-source parity under both
contact models and actual changes to the native pour field and speed when the
original viscosity curve varies with physical shear. They do not establish
measured-fluid calibration, a laminar transition threshold or 3-D breakup.


## Implicit-viscosity nonlinear capillary flight

The existing generic `NonlinearFreeJet` consumes the integration policy of its
original native or material-bound `SlenderJetEvolution`. Selecting implicit
viscosity before flight allows stiff Newtonian jets to advance without the
explicit viscous timestep restriction, while retaining the same conservative
geometry, liquid source, initial profile, mean transport and residence time.
Flight reports retain the native method and local embedded error indicators.
No flight-specific PDE, retuned amplification score or replacement source is used.

The high-viscosity flight regression covers two resumed chunks, physical
residence/distance, volume and axial momentum conservation, retained initial
profiles and a work budget that refuses the equivalent explicit run.
Capillary/advection stability and maximum-step accuracy requirements still
apply. A slender-domain stop is not a prediction of resolved pinch-off or
post-breakup droplets.


## Finite-wavelength free-surface film stability

`LaminarFilmSpectralStudy` composes the native Nusselt base with
`fs_cheb::orr_sommerfeld::film::PreparedFilmSpectrum`. Actual section widths,
volume flow, material properties and effective acceleration determine h, mean
speed, Re, kh, gravity and capillarity. Rates return in 1/s and rad/s on the
same h/mean-speed clock. The caller's h/width budget explicitly admits the
laterally wide approximation.

One prepared pair of spectral grids serves the entire path/wavenumber query.
The report identifies the worst sampled mode, physical operating range, native
residual and dimensional order difference. Finite kh does not invoke the
Benney long-wave approximation. Native numerical or physical refusal publishes
no partial report, and the controlled entry point observes cancellation.

The study implements the existing `LaminarGrowthObjective` for physical risk
and racing. Direct evaluation retains full native refusal details; the common
objective error preserves the failed section and wavenumber. Resolution
differences are numerical indicators, not certified error bounds. A negative
sampled result does not establish global stability, nonlinear laminarity or
stability of unexamined wavelengths. This local temporal study does not extend
the separate weak-growth finite-path transport approximation.

SI regressions cover finite-wavelength decay, physical clock/dimensional
scaling, the shared growth-objective seam and geometry/wavelength/cancellation
refusal. No material coefficients or physical kernel are replaced in the
vessel adapter.

## Original thermal Carreau film stability

`CarreauFilmSpectralStudy` uses `CarreauPourConditions` to retain the original
thermal Carreau-Yasuda law, actual temperature, density, surface tension and
declared source-rate ceiling. Each actual path section resolves its nonlinear
thickness at prescribed flux through `fs_lbm::film::carreau::CarreauFilm`.
Both native spectral grids then sample that same accepted velocity, shear,
curvature, apparent viscosity and stress tangent. No zero-shear parabolic base
or single effective viscosity is substituted.

The native reference Reynolds number is `rho*q_actual/mu0(T)`; both viscosity
profiles use that same original `mu0(T)`. Source profiles and derivatives are
prepared once per section and reused for every requested SI wavelength. Native
free-surface eigenvalues return through the actual `Ubar/h` clock. Reports
retain physical thickness ranges, actual/requested flux residual, native profile
quadrature indicators, and the original spectral residual/resolution evidence.

Nonconstant sources currently require transition exponent `a>=2`, admitting
regular viscosity coefficient derivatives at the zero-shear surface. Their exact
surface stress-tangent gradient is zero. Exact constant-law limits retain the
original Newtonian behavior at the actual temperature for every admitted source
exponent. Unsupported regularity, source support, lateral-width applicability,
work exhaustion or cancellation returns no partial study.

This is isothermal, planar, locally parallel, two-dimensional temporal stability
of sampled finite wavelengths. It does not infer spatial amplification, global
stability, thermal or viscoelastic response, or turbulent transition. Nonlinear
operating conditions remain explicit rather than being passed through the
constant-property Newtonian objective interface. Focused tests check original
temperature, nonlinear thickness/spectrum against independent stress inversion
and uneliminated-surface QZ, physical similarity including constitutive time,
source limits and cancellation.
