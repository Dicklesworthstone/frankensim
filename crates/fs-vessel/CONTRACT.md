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
conserved flux/viscosity scaling, mode-domain refusals, physical CVaR candidate
selection and the original race integration. Native Rust execution remains
unavailable in this session; numerical equation checks are distinct evidence.

## Purpose and layer

Layer **L6 (HELM)**. Composes battle-tested lower crates end-to-end:
fs-cheb (profile + Orr–Sommerfeld), fs-lbm (free-surface pour,
rheology), fs-race/fs-exec (e-raced screening), fs-robust (canonical
empirical CVaR), and fs-render `volumes` (Woodcock transmittance
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
- `pour::run_pour(rig, contact, law) -> PourOutcome` — one physical rigid-frame pour;
  outcome carries `mass_drift` (worst relative ledger drift),
  `fragments` (2-D connected-component count), `poured_mass` (mass right of
  the weir), `dribble_cells` (wet outer-lip cells), and `mass_field`
  (the same bytes the render binds).
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
  model explicitly; default and MID pours select the rigid-frame model.
- `pour::render_pour(outcome, nx, ny, res)` — fs-render Woodcock
  transmittance bound ZERO-COPY to `outcome.mass_field`.
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

## Determinism class

Fully deterministic. The tilt schedule is closed-form, fs-lbm
free-surface stepping is deterministic, the race jitter is a hashed
counter (no RNG state), the render uses fixed Philox streams keyed by
a constant seed. Bitwise replay is gated (vsl-006).

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
