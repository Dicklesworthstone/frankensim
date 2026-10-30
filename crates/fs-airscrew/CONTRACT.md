# fs-airscrew — CONTRACT

Bead: frankensim-wf-root-guzez.5.12.1 (E4.5-i, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN §5.3 (ROUND 6 steady state). Evidence:
prop-geometry-v1 (E1.6) under the frozen registry.

## Purpose and layer

L3 (sibling of fs-wing; NEITHER depends on the other — fs-flyer owns the
coupling). E4.5-i ships the BEMT rotor kernel in the w-formulation (valid
through J = 0), Prandtl tip/root loss, the warm-started bounded station
solve with per-station convergence receipts, the declared engine torque
curve, and the rotor spin-up step. E4.5-ii adds fs-flyer's
CoupledPropAirframeStep (Aitken candidate A).

## Public types and semantics

- `BladeStation { r_over_r, chord_m, beta_rad }`, `Rotor { radius,
  n_blades, camber, stations }` — provenance lives with the caller
  (prop-geometry-v1 rules).
- `bemt_solve(rotor, rho, V_axial, omega) -> BemtSolution { thrust,
  torque, ct, cq, j, station receipts }` — CT/CQ/J in the
  prop-geometry-v1 conventions.
- `engine_torque_at_prop_nm(omega)` — 12 hp @ 1025 engine rpm through the
  23:8 chain, flat-torque-then-power-limited (DISCLOSED approximation).
- `rotor_spinup_step(I, omega, Q_eng, Q_prop, dt)`.

## Invariants

- The w-formulation momentum balance (V+w)·w = dT/(4πrρF) stays valid at
  V = 0 (static bench states are first-class, not a special case).
- Station marching warm-starts from the inboard neighbor; nonconvergence
  is a TYPED refusal, never a clamped answer.
- The normal quadratic root is rationalized as `2t/(V + sqrt(V² + 4t))`,
  with signed annular loading `t = dT/(4πrρF)`, to retain weak induction
  continuously through zero. Negative inner-annulus loading is compatible
  with positive total rotor thrust; no accepted thrust is clipped to zero.
  Accepted thrust, torque, angle, Prandtl factor and induction are evaluated
  at one identical station state. Both the original velocity tolerance and
  the dimensional annular momentum equation must close at that state.
- Algebraic normal throughflow requires `V >= 0` and `w >= -V/2`, giving
  disk speed `V+w >= V/2` and nonnegative far-wake speed `V+2w`. Empirical
  admission is stricter: accepted braking requires `a=-w/V <= 0.4` for
  `V>0`, exposed as `MAX_BRAKING_INDUCTION`; static states retain `w>=0`.
  The ordinary BEM/high-induction transition at about `a=0.4` is described
  in Moriarty and Hansen, [AeroDyn Theory Manual, NREL/TP-500-36881](https://www.nlr.gov/docs/fy05osti/36881.pdf),
  printed pp. 6–7. A real quadratic root alone does not qualify stronger
  braking. This entry does not import an empirical high-induction correction.
- The ordinary relaxed Picard path retains its warm-start branch. At
  evaluation 96, three nearby, monotonically contracting samples may justify
  a small local secant probe when their measured contraction predicts a
  budget overrun. A sign bracket, compatible local slopes and monotone
  bisection then finish within the SAME 120 force evaluations. The reported
  iteration count includes probes and bisections. No global root search,
  tolerance increase or arbitrary selection among stalled-section roots is
  introduced; ambiguous or unresolved continuation refuses.
- A non-real trial quadratic target is not proof that the nonlinear station
  equation has no physical root: section loading changes with induction.
  Such a trial only proposes a bounded half-retreat toward `w=-V/2`, sharing
  the original evaluation budget. It cannot pass convergence or supply the
  real-target contraction evidence required by the local finish. Unresolved
  continuation refuses numerically; a converged root beyond the empirical
  braking limit receives a distinct wake-domain refusal. No reversed-flow,
  turbulent-wake, vortex-ring or universal turbine model is inferred.
- Prandtl F ∈ (0, 1], falling toward the tip.

## Error model

`rotor-invalid` (stations cap AND cap+1, ordering, ranges),
`operating-point-invalid` (nonfinite/unrepresentable state, V < 0, ω ≤ 0),
`station-did-not-converge` (names the station and unresolved numerical branch),
`station-wake-domain` (high induction or non-normal throughflow),
`station-numerical-invalid` (unrepresentable station arithmetic),
`rotor-dynamics-invalid`.

## Determinism class

Deterministic (det:: transcendentals, fixed iteration order); golden
pinned under `org.frankensim.fs-airscrew.v03-golden.v1`.

## Cancellation behavior

Synchronous pure functions.

## Unsafe boundary

Workspace `deny(unsafe_code)`; none.

## Feature flags

None.

## Conformance tests

`tests/bemt_battery.rs` — V-03 core: the E1.6 HOLDOUT static anchor
(1903 reconstruction at 350 rpm vs the 285 N bench, factor-2 trend band);
CT strictly falling over J 0.1–0.9; η rising to a 0.5–0.95 peak in the
Wright J-class then falling; Prandtl limits + tip bite; spin-up from rest
to the 250–500 rpm equilibrium class against the declared engine curve;
caps at cap AND cap+1; golden at the rail state (J ≈ 0.71).

`tests/station_momentum.rs` independently reconstructs the blade-element
loads and annular momentum from public accepted station states. It covers
the actual flight query that formerly exhausted Picard at positive induction,
retention of the lower branch where three positive roots exist, bounded force
evaluations, and static/advancing states with positive, negative and mixed
annular loading. A constructed `a=0.35` root is recovered despite a non-real
first target; a real `a=0.41` root refuses the high-induction domain. Nonfinite
and unresolved non-real static states refuse explicitly. This is numerical
closure of the declared model, not new validation of the reconstructed rotor
or a universal branch selector.

## No-claim boundaries

- The 1903 radial table is an ESTIMATED reconstruction (prop-geometry-v1
  1903-absence rule); quantitative CT/CQ validation awaits real tables.
- Swirl enters via the section angle only (small-swirl DISCLOSED);
  nonuniform/skew dynamic inflow, blade unsteadiness, and inflow harmonics remain E4.5-ii+.
- The engine curve is a declared flat-torque approximation; the real
  torque-vs-rpm curve is load-bearing and remains an E1 follow-on.
- No airframe coupling here (fs-flyer owns CoupledPropAirframeStep).


## Uniform axial dynamic inflow

`dynamic_inflow::UniformAxialInflow` adds physical axial induction memory
using the Carpenter--Fridovich apparent-air-mass model:
`m_a dv_i/dt = T_blade - 2 rho A v_i (V_axial + v_i)`, with the ideal
normal-disk apparent mass `m_a = (8/3) rho R^3`. The model follows the
normal-flow approximation described in Chen and Hindson, NASA TM-88327,
equation (1), https://ntrs.nasa.gov/citations/19870015049. It is an approximate
uniform-disk closure, not a resolved wake or a fitted arbitrary time constant.

Blade thrust and torque are evaluated at the stored induced velocity using
the same native section force law and radial widths as steady BEMT. No steady
annular momentum solve is substituted for this evaluation. Uniform momentum
transport replaces the steady annular Prandtl closure, so identical steady
radial distributions between the two models are not claimed. Existing BEMT
and station goldens remain unchanged.

A held-input implicit-midpoint step solves the scalar wake quadratic through
its rationalized normal root. It reports induced-coordinate thrust impulse,
transported momentum, apparent-fluid kinetic storage, induced work, transport
work, and numerical residuals. Signed loads and moderate braking remain
explicit; leaving normal throughflow refuses rather than clipping. The model
requires constant positive density and matching rotor radius. Apparent wake
mass is never added to rigid aircraft mass or inertia.

Six native cases establish analytic hover thrust-step and unpowered-decay
refinement, equilibrium/weak-load behavior, signed normal-branch admission,
physical blade-force dependence on stored induction, and impulse/storage/work
balances. The method is second order for held smooth admitted inputs. These
axial reduced-state balances do not establish whole-fluid/body passivity,
skew/crossflow forces, azimuthal harmonics, swirl, wake geometry transport,
blade unsteadiness, or historical calibration.


## Oblique rigid-blade force and hub moment

`oblique::{DiskFrame, AzimuthQuadrature}` supplies a separately selected,
revolution-averaged blade-element law for uniform incident crossflow and
carrier tilt. Each section uses the existing radial widths and native
attached/separated section closure. Its material velocity relative to incident
air is `V_hub + omega_tilt cross r + Omega_absolute r t`. The axial carrier
rate is already included in the absolute shaft speed and is not added again.
The force is `dT n - (dQ/r) t`; its actual arm gives the transverse hub
moment. Positive thrust, physical spin direction and material radial reference
are distinct inputs. Frame admission checks orthonormality without silently
changing geometry. The construction follows the inclined-propeller section
kinematics discussed with figures 93–95 of NASA CR-182121
(<https://ntrs.nasa.gov/citations/19880016976>).

`ObliqueLoads` returns full force, transverse hub moment EXCLUDING axial shaft
torque, required drive torque, profile-drag dissipation, actual section
kinematic power, induced work and a floating-point power residual. The admitted
section law satisfies `F dot V_hub + M_perp dot omega_tilt - Q Omega_absolute
+ T v_i + D_profile = 0` in real arithmetic. Required torque retains its sign;
the aerodynamic axial moment is opposite physical spin. These are local
section-power identities, not a complete fluid/wake passivity result.

Quadrature has 8–512 points in multiples of four, with deterministic antipodal
pairs and trigonometry prepared once. Evaluation allocates no storage and costs
O(radial stations times azimuth points). Nonpositive local tangential speed,
reversed section throughflow, and `|alpha| > 1.3` refuse instead of extrapolating
the native angle clamp. Normal-branch induction admission remains explicit.
The maximum ignored spanwise speed is reported: a two-dimensional strip law
does not model spanwise force. There is no blade phase, flapping, skew-induced
wake distribution, compressibility, swirl, cyclic pitch or section memory.

Five focused native tests in `tests/oblique.rs` cover the original axial limit,
an independent attached-section midpoint quadrature and azimuth refinement,
physical handedness parity, full frame rotation and shaft/carrier power, and
bounded work/reversed-flow/angle/overflow refusals. This establishes numerical
conformance of the declared Estimated section model, not propeller calibration.

## Normal-flow first harmonic wake memory and distributed blade induction

`dynamic_inflow::first_harmonic::{FirstHarmonicInflow, HarmonicInduction}`
adds a physical mean plus two first spatial modes:
`v_i(x,y) = u0 + hx x/R + hy y/R`, in metres per second. The radial
reference and `axis cross reference` define the cosine/sine directions;
they are independent of blade handedness and do not follow blade phase.
Dimensional states avoid artificial wake-memory rescaling when shaft RPM changes.

The selected transport closure is the **normal-flow, nonrigid-wake** limit
of Chen, NASA TM-102219, equations (13b), (14), (17), (18) and (21):
<https://ntrs.nasa.gov/api/citations/19900006622/downloads/19900006622.pdf>.
The apparent modal masses are `m0=(8/3) rho R^3` and
`m1=(16/45) rho R^3`. Mean transport remains the native
`2 rho A u0 (V+u0)`; first-mode transport is
`rho A (V+2u0) h/2`. This explicitly selects the nonrigid-wake harmonic
gain; the alternative rigid-wake assumption has different gains/time scales.
It does not implement the source's general skewed-flow matrix.

Actual transverse hub moments supply the power-conjugate loads
`Q=[T, -M dot e2/R, M dot e1/R]`, all in newtons. This follows directly
from `M=sum(r cross n dT)` and `Q dot u=sum(dT v_i)`. The mean and two
moment channels retain separate generalized impulse residuals. They are
reduced modal balances, not a full world-fluid linear/angular momentum account.
No apparent air mass is added to the solid rotor or aircraft inertia.

The moving disk reference contributes `m1 omega_n J h`, where
`J[hx,hy]=[-hy,hx]` and `omega_n` is the **carrier** axial rate, not rotor
shaft speed. This is tangent-frame coordinate transport with zero power in
real arithmetic. Its impulse and measured work are exposed separately; it
is not a physical source. This tangent-plane transport does not represent
far-wake geometry, tilt-induced wake distortion or a resolved air domain.

`step_held` solves the original mean quadratic followed by the coupled
two-mode implicit midpoint equations in closed form, with constant storage
and work. It exposes actual modal kinetic storage, driving/transport/frame
impulses, induced work, transport work, frame work and numerical residuals.
Transport work is nonnegative on the admitted normal branch. Midpoint is
second order for smooth held data and is not L-stable; a coarse transient
still requires refinement. Initial, midpoint and final whole-disk normal
throughflow must be nonnegative; overflow or inadmissible harmonic amplitude
refuses rather than clipping or reverting to uniform flow.

`oblique::AzimuthQuadrature::evaluate_first_harmonic` evaluates the SAME
native section law at each actual radial/azimuth induction value. Induced
work is summed as `sum(dT v_i)`, including both moment ports, instead of
using thrust times mean induction. The zero-harmonic path preserves the
original oblique arithmetic. Prepared 8–512 point quadrature remains
allocation-free during evaluation, with O(radial stations times azimuth)
work. `DiskFrame::rebase_induction` changes the reference of the same
oriented disk; `wake_force_n` provides the conjugate load normalization.
All original positive tangential-flow, normal-flow and angle-domain guards
remain. Higher harmonics, general skew transport, swirl, flapping, cyclic
pitch, wake contraction, ground effects, radial flow force and source
calibration remain outside this model.

Six native G0/G1/G3 cases in `tests/first_harmonic.rs` establish the analytic
exponential modal response, second-order refinement, zero-work rotating-frame
transport, original zero-mode arithmetic, independent attached-section force
and work integration, full-scene rotation, disk-reference re-expression and
domain/late-endpoint refusals. Harmonic response errors at 16/32/64 steps are
`2.4097506965e-4`, `6.0192955857e-5`, `1.5045067779e-5`. These focused native
results do not establish a full Cargo/DSR, cross-platform or physical
calibration gate.
