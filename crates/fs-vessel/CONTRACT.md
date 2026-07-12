# fs-vessel CONTRACT

Flagship 3 (plan §15.3, bead mye.4): the LAMINAR-POUR VESSEL, smoke
tier. Design a carafe lip so the pouring film stays laminar across a
FAMILY of fluids — stability objective, free-surface validation,
CVaR robustification, e-raced screening, and a render that is the
physics.

## Purpose and layer

Layer **L6 (HELM)**. Composes battle-tested lower crates end-to-end:
fs-cheb (profile + Orr–Sommerfeld), fs-lbm (free-surface pour,
rheology), fs-race/fs-exec (e-raced screening), fs-render `volumes`
(Woodcock transmittance deliverable). No new numerics live here; the
crate's claim is the COMPOSITION with certificates at every joint.

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
  rotating-gravity tilt schedule). Defaults pour decisively
  (g0 = 6e-4, tilt 0.7 rad over 900 steps).
- `pour::run_pour(rig, contact, law) -> PourOutcome` — one pour;
  outcome carries `mass_drift` (worst relative ledger drift),
  `fragments` (Plateau–Rayleigh score), `poured_mass` (mass right of
  the weir), `dribble_cells` (wet outer-lip cells), and `mass_field`
  (the same bytes the render binds).
- `pour::render_pour(outcome, nx, ny, res)` — fs-render Woodcock
  transmittance bound ZERO-COPY to `outcome.mass_field`.
- `robust::robustify(beta) -> RobustReport` — nominal (band-center)
  vs CVaR_β lip optimization over the fluid band (rates ×
  viscosities), both evaluated on the OFF-NOMINAL corners.
  `robust::empirical_cvar` is the exact Rockafellar–Uryasev empirical
  tail mean; empty/non-finite losses and beta outside `(0,1)` panic
  as programmer-contract defects instead of returning fake risk.

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

## Error model

Fixture-scale programmer errors panic (`expect`) — the eigensolver
failing to converge, `Cheb1::build` failing on the smooth carafe, or
render dimensions disagreeing with the outcome buffer. No `Result`
surface: the flagship consumes typed results from lower crates and
converts failure to a loud stop; there is no recovery story at the
smoke tier.

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
- **vsl-002** the pour POURS with a strict ledger: drift 7.96e-14,
  poured mass 46.11, under the rotating-gravity tilt schedule.
- **vsl-003** contact bracket: Neutral vs Wetting band REPORTED
  (poured band 0.00, dribble band 0 on the default rig — the models
  do not differentiate this fixture; the honest row says so); gate is
  strictness of both ledgers only.
- **vsl-004** Carreau band: three fluids pour [7.04, 16.71, 18.13] —
  spread 11.09 (the family responds); all ledgers strict.
- **vsl-005** robustification: robust lip 0.92 stays stable off-band
  (−0.00044) while the nominal lip 1.11 goes UNSTABLE off-band
  (+0.00023) — the flagship claim, gated. The e-raced screen recovers
  the deterministic argmin under jitter and eliminates all 4 dominated
  lips in 76 evals vs fixed-N 2000.
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
