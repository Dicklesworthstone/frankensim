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
- `SCOPE_STATEMENT` — the no-absolute-SPL / 2D-to-3D-span-correction
  law as data (the marketing-mutation guard asserts it).

## Invariants

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
    HONEST OPEN SCOPE: the rig is MULTI-STABLE — attractor selection
    varies with seed and resolution (recorded states 0.0366 and
    0.0458) — so strict cross-resolution convergence of the selected
    attractor and a stage-II lock both need a hysteresis-following
    (adiabatic ramp) protocol, recorded on the bead.
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

## Error model

Typed `AeroacError` (non-finite, invalid-parameter, not-converged
with disclosed residual). No silent degradation, no fabricated
eigenvalues.

## Determinism class

Deterministic: fs_math::det transcendentals + Dd arithmetic only, no
platform libm, no RNG.

## Cancellation behavior

Synchronous, milliseconds-class. No `Cx` integration (workspace
`frankensim-ccmn`).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/bickley_lbm.rs` (1): the LBM-vs-Rayleigh growth-rate
convergence fixture (invariant 8).

`tests/edgetone_staging.rs` (2, ignored heavy): stage-I Strouhal vs
Brown + Vaik/Paal published values; fine-lattice slit-fix regression
(invariant 10).

`tests/noisetable.rs` (2): sweep + export + synth round trip;
typed refusals.

`tests/jetlab.rs` (3): edge-tone oscillation + diagnostics +
radiation with scope; typed refusals; bitwise determinism.

`tests/aeroac.rs` (11): Wronskian identity; derivative cross-checks;
fsci-special oracle; small-argument limits; Hankel far-field
amplitude; Curle spreading + directivity + scope guard; Curle typed
refusals; Bickley pin self-verification + falsifiers; shooting solver
vs pins + grid convergence; Bickley typed refusals; bitwise
determinism.

## No-claim boundaries (the bead's remaining scope — OPEN)

- Edge-tone staging BEYOND stage I: only the stage-I point
  (Re 144, h/delta = 10) is validated; the stage II/III ladder and
  hysteresis are not exercised.
- Grid-convergence of source spectra across refinement levels: not
  implemented (single-resolution v1).
- Noise tables catalog THIS rig's tonal limit cycle at low Re in
  2D; they are not turbulent flute-noise spectra (that regime needs
  higher Re than the rig currently reaches).
- Quadrupole (volume) sources: dipoles dominate at low Mach over
  rigid surfaces; the Lighthill volume term is out of scope v1.
- Only orders 0 and 1 of the Bessel functions are provided (all the
  2D monopole/dipole machinery needs); higher orders are not
  implemented.
- The Rayleigh oracle is INVISCID (no critical-layer viscous
  correction); exactly-neutral real-c modes are ill-posed for
  shooting and are validated by approach from the unstable side.
