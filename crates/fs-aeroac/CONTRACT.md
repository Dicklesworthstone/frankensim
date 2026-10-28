# fs-aeroac — Contract

## Purpose and layer

Layer L3. Low-Mach aeroacoustic source models (bead
frankensim-fsim-aeroacoustic-sources-9ok02, IN PROGRESS — this crate
is the first slice): the acoustic-analogy half of the honest hybrid
over fs-lbm's incompressible base flows. Everything here is a
SHAPE/SCALING authority; ABSOLUTE SPL IS NEVER CLAIMED (the crate's
scope law, pinned as data and by test).

## Public types and semantics

- `bessel::{j0, j1, y0, y1, hankel0_outgoing, hankel1_outgoing}` —
  deterministic Bessel/Hankel functions built EXCLUSIVELY from
  fs_math primitives (det transcendentals, Dd double-double series
  accumulation below the x = 18 crossover, recurrence-GENERATED
  asymptotic sums above it — deliberately no transcribed coefficient
  tables). `e^{-i omega t}` convention: the outgoing cylindrical wave
  is `H^(1) = J + iY`. Y at non-positive arguments is NaN (real
  domain).
- `curle2d::dipole_pressure` — 2D frequency-domain Curle compact
  dipole over the outgoing Hankel Green's function
  `G = (i/4) H0^(1)(kr)`: `p = (ik/4) H1^(1)(kr) (rhat . F)` per unit
  span. The 2D-vs-3D Green's-function trap from the bead's polish
  round is structurally excluded (frequency domain, Hankel kernel);
  every output embeds [`SCOPE_STATEMENT`].
- `bickley::bickley_rayleigh_mode` — Rayleigh-equation shooting
  solver (RK4 + complex secant, half-line with symmetry boundary
  conditions and the exact `phi' = -alpha phi` far-field decay) for
  the Bickley jet `U = sech^2(y)`: the inviscid instability oracle
  the fs-lbm jet validation chain runs against. Non-convergence
  REFUSES with the residual disclosed; no partial eigenvalue.
- `bickley::rayleigh_residual_closed_form` — the self-verification
  surface for the analytic pins.
- `jetlab::run_jet_labium` — the edge-tone base-flow fixture: an
  fs-lbm slot jet (smoothed top-hat) impinging on a two-cell splitter
  plate (DISCLOSED simplification of a wedge — the sharp edge, not
  the wedge angle, drives the oscillation class), in a PERIODIC-x
  domain closed by a per-row-profile FRINGE layer (the spectral-DNS
  fringe method: the sponge re-conditions outflow to the authored
  inflow profile, and doubles as the measured acoustic absorber).
  Records the plate force series via momentum exchange; returns
  per-run diagnostics (max Mach, plate-vs-fringe plane fluxes,
  Reynolds) and embeds the scope statement.
  `jetlab::dipole_spectrum_line` radiates a caller-FFT'd force line
  through the 2D Curle dipole.
- `jetlab::transverse_force_peak` — the shared Hann-periodogram peak
  instrument (mean removal, guard-band skip, prominence vs median
  admitted power, unwindowed RMS floor value); reports the peak BIN
  so the +-1/(2 bin) quantization travels with every Strouhal claim.
- Lattice destabilization is a TYPED refusal: a blow-up inside
  fs-lbm's stepping (its density assert panics on the first bad
  cell, so no pre-scan can precede it) is contained at the
  settle/ramp/record boundary via `catch_unwind` and returned as a
  named `NonFinite` refusal; the rig state is discarded on that
  path. Executed trigger: plain-BGK destabilizes at cell Reynolds
  `Re/delta` in (36, 48) on this rig family (measured at the fixed
  rig speed u = 0.08 across delta 6 and 12; the Mach/u dependence of
  the boundary was not probed and is not claimed).
- `jetlab::run_adiabatic_ramp` — the hysteresis-following protocol
  for the MULTI-STABLE rig: after one base lock-in, Reynolds is
  swept quasi-statically (per-step linear relaxation-time ramp at
  FIXED jet speed, so Mach, fringe profile, nozzle slit, and seed
  are all unchanged) up to `reynolds_end` and back, with a spectral
  rung measured after every transition. Returns `2 rungs - 1`
  direction-tagged measurements; adiabaticity (ramp steps >> the
  oscillation period) is the caller's declared obligation. Refuses
  degenerate spans and a top-of-ramp tau at or below the 0.5005
  stability floor.
- `SCOPE_STATEMENT` — the no-absolute-SPL / 2D-to-3D-span-correction
  law as data (the marketing-mutation guard asserts it).
- `regime` — bead l011o: pinned CentralMoment 2D tonal-lock rows
  (`PINNED_2D_CENTRAL_MOMENT_TONAL`), the named inverse-cascade mechanism,
  `admit_broadband_spectrum` (refuses flatness at or below
  `TONAL_FLATNESS_CEILING`), `measure_spectral_flatness` /
  `classify_spectrum` (geometric/arithmetic power-mean; no table is
  minted), and `evaluate_slot_jet_3d_operator` (D3Q19 CentralMoment
  smoke). `SlotJet3dOperatorSmoke::incorporate_measured_spectrum`
  raises `broadband_demonstrated` only from a real spectrum that
  clears the floor. The smoke alone does not mint a 3-D broadband
  table.
- `slot_jet_3d` — bead frankensim-music-v8-root-3ez8g.10.1: the 3-D
  slot-jet rig on the D3Q19 central-moment operator, the executed
  follow-up lane to the recorded 2-D tonal lock. Periodic domain +
  right-side per-y-row fringe (`Fringe3`, the Sponge2 port; its 3-D
  reflection coefficient is MEASURED by the test battery with a
  bounce-back-wall control), binary-slit nozzle wall (the executed
  2-D lesson kept), two-cell splitter plate across the full span,
  momentum-exchange plate force, `run_slot_jet_3d` /
  `classify_rung` / `SlotJet3dRung::to_jsonl` receipts with the FFT
  bin disclosure and the `FORCE_RMS_AMPLITUDE_FLOOR` amplitude law.
  Chunked checkpoints fail closed on unknown format versions, files beyond the
  configured size envelope, inconsistent settle/record/force counters, any
  truncated or trailing state bytes, non-finite force/population values, and
  non-positive fluid density before resumed stepping can mutate the rig. Their
  identity covers the full run configuration, including settle and record
  horizons; terminal checkpoints are retained for idempotent replay.
  The Re actuator is the central-moment second-order rate at fixed
  `u_jet`; independent rungs, not adiabatic continuation. Both
  outcomes are honest wins: a broadband rung above the tonal
  ceiling or a quantified tonal refusal. Driver:
  `src/bin/slot_jet_3d_sweep.rs` +
  `scripts/e2e/music/slot_jet_3d_sweep.sh`. No experimental,
  video-backed, or absolute-level claim survives this lane; the
  crate scope statement applies to every spectral quantity.
- `jetcard::{JetCard, JetCardClaim, CardAuthority, mint_tonal_interim_card,
  mint_broadband_card, mint_refusal_boundary_card}` — minted jet
  cards (bead 3ez8g.10.2), the typed boundary object between this
  lab lane and the audio-rate performance lane under the downstream
  "no card -> refuse" menu law. Authority law enforced at minting
  and revalidated on parse: `X-Struct` requires a residual issued
  against lab data, `X-Est` otherwise; a broadband claim requires
  16-band noise content; every card embeds the scope statement and
  refuses queries outside its demonstrated Reynolds range. The
  serialization is a pinned-field-order v1 JSON document with a
  strict fail-closed parser (foreign schema/claim kinds, tampered
  scope, trailing bytes all refuse by name); `content_hash` is
  FNV-1a over those canonical bytes — a deterministic correlation
  identity, NOT a cryptographic or authenticity anchor (the
  workspace `ProvenanceHash` caveat applies verbatim). The tonal
  interim card is minted from the recorded stage-I staging point
  (St 0.03662 vs Brown 0.03554, +3.0% inside the ±6% bin; ramp
  hysteresis and 0.0366/0.0458 multi-stability recorded) with the
  exact staging-rig FNV fingerprint as provenance.
- `forced_jet::{ForcedJetResolution, ForcedJet, JetStepResponse,
  bickley_step_response, FRONT_MARGIN}` — beads pexfp.135/.138/.321,
  ruling 7376: the causal response of a planar Bickley jet to a uniform
  transverse acoustic velocity at a baffled flue exit, integrated in
  time from the linearised vorticity equation (Thom wall at the exit:
  the jet leaves undisplaced and parallel; the exact discrete harmonic
  exterior beyond the lateral band; a downstream wall kept
  `FRONT_MARGIN` beyond the disturbance front). Units: `b`, `U0`,
  `b / U0`, `Re = U0 b / nu`. `bickley_step_response` records the
  centreline displacement at one station for a unit velocity step from
  rest over the declared window and holds it afterwards; `transfer(w)`
  is that held response's exact Fourier transform and `hat_weights`
  the audio-rate convolution weights in seconds for a mouth velocity
  linear between samples (`g_0` multiplies the current sample, the
  weights end one sample after the window and sum to the held gain).
  `leading_hat_weights` returns their count and only the first `count`
  of them, bit for bit, for a slow jet whose caller reads only its
  leading lags.
  `interval_hat_weights` integrates this SAME interpolant's derivative
  against both physical-time endpoint hats on a supplied convective interval.
  It returns newer/older weights in seconds, leaves the held tail inactive,
  and refuses invalid ages, advances or durations. At uniform spacing,
  adjacent interval weights recompose the existing hat weights. The consumer
  owns the accepted clock history; this method adds no variable-Re flow law.
  Refuses non-power-of-two or short boxes (`length < window + 8`), an
  explicit viscous number `4 dt / (Re dx^2)` above 2, stations outside
  the box, and non-finite results.

## Invariants

### Shared material integration

`material::MaterialJetLabium` binds geometry and physical speed/dx/dt inputs
to the shared `GasState` and `fs-lbm::material::MaterialScaling` owner. Gas
coherence is admitted before compilation. The existing jet-labium solver
receives material-derived viscosity; its actual surface force is converted
to N/m with `rho * dx^3 / dt^2`. Physical frequency and the same gas's sound
speed supply `k = 2 pi f / c` to the existing Curle observer. No independent
rho, viscosity, Reynolds number or acoustic speed can override the binding.
Pure, humid, mixed and thermochemical gas evaluations use the same path.
`MaterialJetLabiumRun::observe_spectral_line` transforms its own accepted
force history with `fs-fft`, mean removal and the existing symmetric Hann
window. Coherent-gain normalization preserves resolved sinusoid amplitude;
conjugation makes the phase convention `exp(-i omega t)` relative to the first
recorded sample. Physical bin frequency and same-medium pressure are returned
together. An independent direct DFT checks force amplitude and phase.

The physical Mach envelope is strictly below 0.3 and is checked separately
from numerical lattice Mach. Grid allocation is budgeted before construction.
An observed speed-envelope violation refuses output; diagnostics sample the
existing two planes every 16 steps and do not prove the full-field bound.
Observer frequencies at or above the force record's Nyquist limit refuse.
The medium is homogeneous and frozen during a run. Compressible feedback,
reaction, relaxation/absorption, global flow validity and quantitative
experimental accuracy are not established by this adapter; the crate's
two-dimensional shape/scaling scope remains unchanged.

`tests/material.rs` verifies the actual lattice-to-SI force/frequency mapping,
same-state radiation, independent physical/numerical Mach admission, and
data-only pressure/temperature/mixture changes through executed flow and
observer outputs (G0/G3). These short transient tests establish integration,
not a converged tonal spectrum or absolute sound level.

`material::MaterialSlotJet3d` applies the same material/scaling owner to the
existing D3Q19 central-moment slot-jet rig. Its second-order relaxation rate is
derived from the gas's actual kinematic viscosity; only higher-order numerical
relaxation remains an independent control. Three-dimensional momentum exchange
converts with `rho * dx^4 / dt^2` into N across one explicitly retained periodic
span. Both material adapters share one Hann/coherent-gain/phase FFT extraction.
The physical line frequency and acoustic wavenumber use the same timestep and
gas as the flow. A three-dimensional periodic-cell force is not passed to the
two-dimensional Curle Green function or relabeled as microphone pressure.

The 3-D initial-speed admission includes the rig's 1.5 maximum spanwise seed
factor. Final-plane speed diagnostics use all three velocity components
(`momentum / density`), and plane mass flux is momentum itself. These are
final-plane samples, not a global speed bound or maximum over the recording.
`tests/material_slot_jet_3d.rs` executes pressure, temperature and composition
substitutions through the original solver, independently checks SI force and
DFT amplitude/phase, and checks necessary material/Mach/allocation refusals.
This frozen homogeneous low-Mach adapter supplies no compressible feedback,
chemical evolution, acoustic radiation, turbulence-convergence or SPL claim.

1. EXACT Wronskian: `J1 Y0 - J0 Y1 = 2/(pi x)` holds within 5e-13
   over a log grid 1e-3..3e3 including the series/asymptotic
   crossover from both sides. Honest limits (review-measured against
   a 120-digit reference): the identity is blind to COMMON phase
   errors, and pointwise accuracy beyond the certified grid degrades
   as ~x*eps from f64 phase representation (~2e-11 at x = 1e6) —
   stated in the module doc; intended kr ranges sit inside the
   certified band.
2. Independent-path derivatives: central-difference `J0'` matches
   `-J1` and `Y0'` matches `-Y1` at FD accuracy (1e-8) across both
   regimes.
3. Cross-implementation oracle: fsci-special's cephes-heritage
   j0/j1/y0/y1 agree within 5e-12 of envelope (fsci's plain-f64
   Y-series carry ~1e-12 cancellation near their own crossover —
   measured; the Wronskian is the precision arbiter).
4. Small-argument limits exact (J0(0) = 1, J1 ~ x/2, bounded Y0 log
   constant); far-field `|H0| -> sqrt(2/(pi x))` within the O(x^-2)
   envelope correction.
5. Curle dipole physics: far-field decay exponent is -1/2
   (CYLINDRICAL spreading; the 3D-Green's-swap mutation reads -1 and
   fails the gate), the perpendicular directivity null is EXACT, the
   45-degree amplitude matches cos(theta) within 1e-3, and the scope
   statement rides on every output.
6. Bickley analytic pins are SELF-VERIFIED per run: `phi = sech^2` at
   `(alpha, c) = (2, 2/3)` and `phi = sech tanh` at `(1, 2/3)` drive
   the closed-form Rayleigh residual to machine zero (< 1e-14) at
   every probe point, and the residual is proven LIVE (perturbing
   alpha or c by 0.01-0.05 moves it above 1e-3). No literature value
   is transcribed on trust.
7. The shooting solver reproduces the pins from the unstable side
   (alpha = 1.95 sinuous: c within 0.02 of 2/3 with small positive
   Im c; alpha = 0.95 varicose: within 0.03), shows positive
   mid-band growth exceeding the near-neutral growth, and the
   eigenvalue is grid-converged (2048 vs 4096 RK4 steps within
   1e-8).
8. VERIFICATION CHAIN link 1 (dev-dep composition with fs-lbm): the
   D2Q9 Bickley jet's measured linear growth rate CONVERGES toward
   the inviscid Rayleigh oracle under simultaneous grid + Reynolds
   refinement — measured -21.4% (b = 8 lu, Re 240) closing to -12.0%
   (b = 12 lu, Re 432), both biased LOW as finite-Re physics demands,
   with clean-exponential window consistency, low-Mach and linearity
   diagnostics green per run, and gates authored just outside the
   measured values (a fixed loose envelope would hide a broken
   convergence direction). Executed lesson in the test doc: the free
   jet's own viscous diffusion (~nu/b^2 per step) decays the base
   flow and drags the measured rate down over long runs — fit early,
   at high Re.
9. Jet-labium fixture: THE VACUOUS-OSCILLATION TRAP (executed and
   now gated): the mirror-symmetric rig, unseeded, preserves
   symmetry to roundoff and its force spectrum shows high-prominence
   structure in ~1e-15 amplified machine noise — prominence-style
   ratios cannot detect it, so every amplitude-bearing test now
   REQUIRES the deterministic sinuous seed (`seed_amplitude`) and a
   force-RMS floor (1e-6 lattice). Seeded, the default configuration
   (Re 240, edge at 3 slot heights) saturates into a REAL limit
   cycle (Fy rms >> floor, essentially pure tone), max Mach 0.07,
   flux imbalance under 1% (the outlet-reflection pathology read
   6%). Geometry/regime refusals typed; bitwise deterministic.
10. EDGE-TONE STAGING (ignored heavy test, executed on record;
    seed-provenance trio in the test doc): with a nozzle wall at the
    jet root and the 0.005 sinuous seed, the SATURATED rig at
    Re 144, h/delta = 10 locks stage-I St = 0.03662 vs Brown's
    (1937) 0.03554 — +3.0%, within the record's +-6% bin
    quantization and INSIDE the published spread (Vaik/Paal exp
    0.03723, CFD 0.03497; two fetched sources). The unseeded run
    selects the SAME frequency in amplified roundoff
    (frequency-selection only); a 0.02 seed lands a NEIGHBORING
    locked state (St 0.0458) — edge-tone multi-stability, recorded.
    STRUCTURAL FINDING: without the nozzle the rig locks to the free
    jet's own mode (St 0.46, 12x the ladder) — the jet-root
    receptivity edge closes Brown's loop. SLIT-SHOULDER FIX
    (executed): with a nozzle the fringe target is the BINARY slit
    profile matching the wall opening (a smooth target slammed into
    the shoulder every wrap and fed a slit-lip mode at St_delta 2-3
    that blocked lattice refinement); the fine-lattice regression
    test (delta = 7.5, geometric similarity, recorded St 0.0458 =
    1.29x Brown) gates the LADDER BAND [0.7, 1.4] x Brown.
    MULTI-STABILITY (executed): attractor selection varies with seed
    and resolution (recorded states 0.0366 and 0.0458), so strict
    cross-resolution convergence of the selected attractor is not
    claimable from fixed-parameter runs; the adiabatic ramp
    (invariant 13) is the follow-up instrument.
13. ADIABATIC RAMP PROTOCOL (ignored heavy test, executed on record
    at TWO ramp rates — 15 and 3 stage-I periods per transition —
    with qualitative agreement): branch-following over
    Re 144 -> 264 -> 144 at h/delta = 10 stays inside the stage-I
    NEIGHBORHOOD [0.6, 1.7] x Brown on every rung of both legs
    (executed extremes 0.717/1.546; the slit-lip >5x, free-jet >12x,
    and stage-II ~2.3x signatures never appear), wanders between the
    neighboring 0.0366/0.0458 locked families (mode competition),
    and returns to Re 144 on a DIFFERENT attractor than it started
    from (bin 9 vs bin 7 on the committed protocol) — REAL,
    bitwise-reproducible hysteresis. NO stage-II lock occurs by
    Re 264 despite the literature onset ~220-250 for a wedge: a
    recorded NO-CLAIM boundary of this two-cell-plate periodic-
    fringe 2D rig, not a validated staging result.
11. NOISE TABLES (the product deliverable): `noisetable` fits
    Strouhal-band power-DENSITY shapes (record-length-independent —
    a band-SUM convention broke the synth round trip by the log-band
    width ratio, executed) over a velocity sweep with per-entry
    regime gates and real-amplitude floors; the JSON export embeds
    the scope statement (marketing-mutation guard) and the MEASURED
    power exponent as data (the saturated tonal limit cycle is NOT
    monotone in u near mode switches — multi-stability, reported not
    prescribed); the demo synth consumes the table and its output's
    band densities match within 6 dB (measured worst 2.5 dB).
12. Refusals typed by name (non-positive wavenumber, coincident
   observer/source, non-finite inputs, bad solver parameters,
   non-convergence with residual); Hankel functions outside their
   domain are FULLY NaN (both components — a half-valid complex once
   leaked, review-caught); the Rayleigh solver's result is
   documented GUESS-DEPENDENT (any eigenvalue of the symmetry class
   zeroes the mismatch); determinism bitwise.
14. FORCED JET (`tests/forced_jet.rs`, executed on trj lane q): at
   station 11 b and Re 600, the step kernel's transfer equals
   independent harmonic forcing from rest within 0.1 % at w 0.5. At
   w 1-2, holding the response after window 90 cuts a slowly decaying
   tail. That costs 2.5-3.9 % (1.2 % at w 3; the test holds w 2 to
   4 %). At w 2, window 150 halves it, while the harmonic side moves
   1e-3, and spacing and time step leave it at 3.2-3.8 %;
   the kernel converges at second order on the ladder dx 0.2/0.1/0.05
   (window 90 b/U0; change 0.1 -> 0.05 is 0.07-0.11 % at w 0.03-0.16,
   0.7-1.9 % at w 0.5-1, held gain 14.8625 -> 14.8663 b/U0) and is bit
   identical at 1 and 5 workers; doubling the box (102.4 -> 204.8 b)
   or narrowing the band (8 -> 6 b) moves it at most 3.3e-3 / 1.9e-3
   at w >= 0.06. Inviscid, the downstream wave carries the published
   sinuous spatial mode at w 0.8 (Sabatini & Bailly 2015,
   1.449709 - 0.134109i): wavenumber within 0.1 %, growth rate
   converging from below with the critical-layer resolution (4.4 % low
   at dx 0.1 with window 125; a dx 0.05 prototype is 1.2 % low). The
   hat weights for a recorder-like jet (22.5 m/s, 0.4 mm, 44.1 kHz)
   are 73 causal taps whose current-sample weight is the early
   `s ~ t` response (0.5 % from half a period) and whose sum is the
   held gain to 1e-9. SnowySummit's causality admission: the transfer
   sampled on [0, 10] b/U0 and inverted puts less energy and L1 before
   time zero than any of the probe's pure-delay controls (5, 10, 15,
   20 b) on the same grid: energy 1.1e-5 against 4.2e-4 and L1 2.1e-3
   against 0.15 untapered, 3.0e-9 against 2.6e-6 and 6.8e-5 against
   0.015 under the 0.25 Gaussian taper (prototype kernels hold the same
   level over windows 30-150 b/U0). Windows stay at
   or below ~150 b/U0: beyond ~200 the downstream packet's growth
   brings FFT round-off into the global Poisson solve and every
   station blows up at once; from 90 to 200 the transfer moves at most
   3 % and 4 degrees at w >= 0.04.

## Error model

Typed `AeroacError` (non-finite, invalid-parameter, not-converged
with disclosed residual). No silent degradation, no fabricated
eigenvalues.

## Determinism class

Deterministic: fs_math::det transcendentals + Dd arithmetic only, no
platform libm, no RNG. `forced_jet` splits its rows across
`std::thread::scope` workers; every row of both transforms is computed
whole by one worker, so the result is bit-identical for every worker
count (tested at 1 and 5).

## Cancellation behavior

Synchronous, milliseconds-class. No `Cx` integration (workspace
`frankensim-ccmn`). Exception: `bickley_step_response` is
seconds-to-minutes class (window 90 b/U0 at dx 0.1: 100-175 s on a
loaded 8-worker trj; dx 0.2: 52 s) and is not interruptible; callers
poll their gate before and after it.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/bickley_lbm.rs` (1): the LBM-vs-Rayleigh growth-rate
convergence fixture (invariant 8).

`tests/edgetone_staging.rs` (4, ignored heavy): stage-I Strouhal vs
Brown + Vaik/Paal published values; fine-lattice slit-fix regression
(invariant 10); adiabatic ramp hysteresis protocol (invariant 13);
turbulent-regime probe (stability-boundary bracket + tonal-window +
typed-destabilization pins).

`tests/noisetable.rs` (2): sweep + export + synth round trip;
typed refusals.

`tests/jetlab.rs` (5): edge-tone oscillation + diagnostics +
radiation with scope; typed refusals; bitwise determinism; ramp
structure + rung grid + bitwise determinism; ramp + peak-instrument
typed refusals.

`tests/jetcard.rs` (13): tonal interim card mints the recorded
stage-I point with X-Struct residual + narrow claim + validity
refusals; exact round trip + stable content hash; broadband/refusal
minting laws (no broadband rung -> refuse; out-of-regime rung ->
refuse; boundary needs >= 2 admitted all-tonal rungs); the
no-residual X-Est cap and the forged-authority refusal; parser
refusals on foreign schema, laundered scope, truncation, and
trailing bytes; momentum-thickness tanh-edge law + refusals; pinned
schema constant; the sweep-receipt reader (`SlotJet3dRung::from_jsonl`,
`parse_sweep_receipts`): bitwise round trip incl. a non-finite flux
imbalance, pinned-order/schema/trailing-byte/non-finite refusals,
header-first and foreign-schema file refusals, and every archived
`tests/receipts/slot-jet-3d-re-sweep*.jsonl` parsing; the
temporal-Nyquist law (`peak_at_nyquist_edge`, `NYQUIST_EDGE_BINS`):
a rung whose raw peak sits within 8 bins of the record's Nyquist edge
is the bounce-back parity artifact (executed at rate 1.60: bin 4095 of
4096) and no minter admits it; the parity filter
(`parity_filtered_force_series`, `classify_rung_parity_filtered`):
non-overlapping pair averaging annihilates the period-2 line exactly,
keeps a tone's bin index and Strouhal (record duration unchanged; the
pipeline Strouhal is halved), refuses odd/short records.

`tests/forced_jet.rs` (6, full suite minutes-class): the published spatial mode
downstream; step kernel against harmonic forcing; the dx 0.2/0.1/0.05
ladder and worker-count bits; box and band invariance; hat weights and
the causality admission against a pure delay (invariant 14); local interval
weights against uniform hats at three clocks, tail support and invalid clocks.

`tests/aeroac.rs` (11): Wronskian identity; derivative cross-checks;
fsci-special oracle; small-argument limits; Hankel far-field
amplitude; Curle spreading + directivity + scope guard; Curle typed
refusals; Bickley pin self-verification + falsifiers; shooting solver
vs pins + grid convergence; Bickley typed refusals; bitwise
determinism.

## No-claim boundaries (the bead's remaining scope — OPEN)

- Edge-tone staging BEYOND stage I: only the stage-I point
  (Re 144, h/delta = 10) is validated. Hysteresis and stage-I-band
  confinement under adiabatic ramping ARE now measured (invariant
  13), but no stage-II/III LOCK exists on this rig by Re 264 —
  stage-ladder claims above stage I remain refused.
- Grid-convergence of source spectra across refinement levels: not
  implemented (single-resolution v1).
- Noise tables catalog THIS rig's tonal limit cycle in 2D; they are
  not turbulent flute-noise spectra. QUANTIFIED (executed probe,
  delta = 12 lu): the plain-BGK stable window ends at cell Reynolds
  `Re/delta` in (36, 48), and the spectrum is still an essentially
  pure tone (flatness ~1e-15) at the last stable rung. The fs-lbm
  `CollisionModel2::CentralMoment` operator (bead 3zkcr) removes the
  stability wall — the same rig runs through Re 2304 (Re_cell 192),
  now selectable via `JetLabiumConfig::collision` — but every
  executed rung REMAINS an essentially pure tone (flatness ~1e-18;
  free-jet-scale St ~ 0.467 at Re >= 1152). Bead l011o records that
  landscape as [`crate::regime::PINNED_2D_CENTRAL_MOMENT_TONAL`] and
  [`crate::regime::two_d_broadband_refusal`]: cataloging those rows as
  broadband flute noise is a typed refusal. The D3Q19 central-moment
  operator is live ([`crate::regime::evaluate_slot_jet_3d_operator`])
  but no 3D Re-sweep has been executed; a demonstrated 3D broadband
  regime remains a no-claim. The 2D inverse-cascade mechanism is the
  recorded physics boundary.
- Jet cards: the tonal interim card claims EDGE-TONE CLASS ONLY —
  it is not flute broadband, carries no band noise content, no
  absolute SPL, and no oscillation-amplitude authority; consuming
  it as a noise source is a type error by construction. The
  broadband and refusal-boundary minters are law + schema today:
  no card of either kind has been minted from executed 3-D sweep
  receipts (the .10.1 sweep is the pending producer). The card
  content hash is FNV-1a: deterministic correlation identity only,
  never authenticity.
- Quadrupole (volume) sources: dipoles dominate at low Mach over
  rigid surfaces; the Lighthill volume term is out of scope v1.
- Only orders 0 and 1 of the Bessel functions are provided (all the
  2D monopole/dipole machinery needs); higher orders are not
  implemented.
- The Rayleigh oracle is INVISCID (no critical-layer viscous
  correction); exactly-neutral real-c modes are ill-posed for
  shooting and are validated by approach from the unstable side.
- Forced jet: the base flow is parallel (its own viscous spreading
  over the window is neglected; the disturbance is viscous); 2-D (no
  span ends); linear (saturation belongs to the consumer); the exit is
  a baffle wall (no channel profile development, no chamfer); the
  labium does not react back on the jet, which runs on freely past
  the station; the response is held after the window, so content
  below w ~ 2 pi / window is the held limit, not computed, and the
  slowly decaying tail the hold cuts costs 2.5-3.9 % at w 1-2 at
  window 90 (invariant 14); the
  spatial growth rate is low by the critical-layer resolution (4.4 %
  at dx 0.1, invariant 14).
