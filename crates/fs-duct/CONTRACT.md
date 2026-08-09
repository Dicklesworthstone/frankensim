# fs-duct — Contract

## Purpose and layer

Layer L3 (FLUX). 1D viscothermal duct/horn acoustics
(bead frankensim-fsim-duct-acoustics-zdmm1): transfer matrices over
cylinder/cone segment chains with Zwikker–Kosten wide-tube losses,
low-`ka` radiation terminations, and input impedance. Generic by
doctrine: instrument bores, mufflers, HVAC runs and waveguides are the
same object, and every medium property derives from
`fs_material::gas::GasState` — no hardcoded air.

## Public types and semantics

- `Segment::{Cylinder, Cone}` — geometry in SI; a cone is a linear
  taper propagated with exact 1D spherical waves, its viscothermal
  correction evaluated at the mean radius (documented standard
  treatment; refine by subdivision).
- `LossModel::{Lossless, WideTube}` — the exact closed-form arm and
  the first-order ZK wide-tube arm.
- `Termination::{Closed, IdealOpen, UnflangedOpen, FlangedOpen}` —
  ideal limits plus the classic low-`ka` radiation fits (end
  corrections 0.6133a / 0.8216a) with a named `ka` ceiling.
- `segment_wave` — complex `k` and characteristic impedance from the
  gas state: `k = k0 [1 + (1+i)/(sqrt2 rv)(1 + (gamma-1)/sqrt Pr)]`,
  `Zc` with the `(1 - (gamma-1)/sqrt Pr)` factor; shear number
  `rv = r sqrt(rho omega / mu)` reported, refused below 10.
- `input_impedance` / `impedance_sweep` / `impedance_peaks` —
  `Z_in = p/U` at the inlet with per-solve diagnostics (minimum shear
  number, mouth `ka`).
- `DuctError` — stable `FS-DUCT-*` refusals: bad parameter, too-narrow
  (wide-tube floor), radiation-`ka` ceiling, empty duct, singular.

Time convention `e^{-i omega t}` (matching `fs_bem::helmholtz`):
`Im k > 0` decays, mass-like reactance is negative imaginary, closed
pipe `Z_in = +i Zc cot(kL)`.

## Invariants

1. Segment 2-ports are built numerically from exact analytic basis
   solutions (plane/spherical waves with the transmission-line
   relation `U = S p'/(i k z_specific)`) — never transcribed matrices.
   Oracles: lossless closed/open cylinder matches `+i Zc cot(kL)` /
   `-i Zc tan(kL)` to 1e-9; the cone reproduces the cylinder in the
   zero-taper limit; the (p, U) determinant is exactly 1 (constant
   Wronskian, `r` proportional to apex distance); half-cone
   composition equals the whole cone to 1e-9.
2. A near-complete closed-open cone resonates at ALL harmonics
   `n c / (2 L_apex)` (the classic conical-bore result, within 3% on
   the committed fixture).
3. The lossless unflanged quarter-wave ladder lands on
   `(2n-1) c / (4 (L + 0.6133 a))` within 0.5%; the viscothermal
   ladder sits below it by the independently computed dispersion
   deficit (measured ratio 1.006/0.998/1.000 across three peaks).
4. Resonance Q of a closed viscothermal cylinder matches the
   closed-form Kirchhoff wall-loss `Q = k/(2 alpha)` within 3%
   (measured 0.6%) — the loss model validated through the full TMM
   machinery against an independent expression.
5. The thermal loss share for air is the `(gamma-1)/sqrt(Pr)` factor
   (measured 1.475x viscous): a dropped thermal term fails loudly.
6. `Re Z_in >= 0` across the band for radiating viscothermal ducts
   (sign pin); lossy peaks flatten below lossless ones (dispersion
   sign pin).
7. Ambient parameterization: lossless resonances scale with
   `sqrt(T_hot/T_cold)` against the independently computed constant.
8. Repeat evaluations are bitwise identical.

## Error model

Typed `DuctError` refusals with stable codes; validity is refused by
name (narrow tubes, radiation `ka`), never silently degraded.
Diagnostics (`min_shear_number`, `mouth_ka`) ride on every response.

## Determinism class

Bit-deterministic across runs on a platform: sequential fixed-order
assembly, `fs_math::det` transcendentals, dense 2x2 complex solves via
`fs_la::eigen_complex::lu_complex`.

## Cancellation behavior

Pure short-running functions; sweeps are caller-chunkable per
frequency. No `Cx` integration (tracked with the workspace-wide
`frankensim-ccmn` effort).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

Inline `tests` module: lossless closed forms; quarter-wave end
correction + dispersion deficit; four cone oracles; Kirchhoff-Q;
passivity/flattening/thermal-share pins; hot-duct sqrt(T) law; named
refusals; bitwise determinism.

## No-claim boundaries

- Wide-tube first order ONLY: narrow tubes (`rv < 10`) refuse; the
  full Bessel/Kelvin ZK solution is the recorded follow-up (needs
  fs-math Bessel functions, shared with the piston closed form).
- Straight rigid smooth isothermal walls: no wall compliance, no
  roughness, no porous liners, no mean flow, no nonlinear (high
  amplitude) propagation.
- Cone losses use the mean radius (documented approximation); strongly
  flaring horns may need the multimodal expansion — the recorded
  trigger is a bell mismatch beyond authored tolerance in a future
  validation against measured instrument impedance.
- Radiation fits are low-`ka` (ceiling 1.0, refused by name) for
  unflanged/flanged circular mouths; BEM-computed loads are the
  recorded successor (`fs_bem::helmholtz`).
- No tone holes / side branches / junction elements yet — slice 3 of
  the bead, along with the provenance-pinned measured-impedance
  validation of a full fingering lattice.
- The gas state inherits `fs_material::gas` boundaries (ideal gas,
  phase validity unchecked, calorically perfect gamma).
