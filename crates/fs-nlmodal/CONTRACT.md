# fs-nlmodal — Contract

## Purpose and layer

Layer L3. Geometric nonlinearity for thin structures in modal
coordinates (bead frankensim-fsim-vonkarman-u5nme): von Karman plates
and Kirchhoff-Carrier strings as port-Hamiltonian systems with a
quartic SUM-OF-SQUARES potential
`H = 1/2 sum (p^2 + w^2 q^2) + 1/4 sum_j c_j (q^T E_j q)^2` — the
Airy-stress-eliminated modal form. Pitch glide, mode coupling, and
the crash cascade with structurally exact energy accounting.

## Public types and semantics

- `SosModalStorage` (implements fs-phs `Storage`): frequencies +
  `StressChannel`s (`c_j >= 0`, symmetric `E_j`). Gradient is the
  EXACT analytic gradient of the coded `H`; state layout matches
  fs-phs `modal_bank` (`[q, p]` interleaved).
  `validate()` owns frequency/channel admission for `assemble` and
  consumers that retain this storage inside multiport/full-viscosity systems.
  Its Gonzalez midpoint energy remainder is evaluated directly as
  `sum_j c_j/4 (d^T E_j d) (m^T E_j d)`, where `d = q1-q0` and
  `m = (q1+q0)/2`; every quadratic term contributes exactly zero.
  This retains small resolved restoring forces beside large stored modal
  energy without subtracting total endpoint energies. It uses one triangle
  traversal per symmetric stress channel and no per-evaluation allocation,
  through the existing fs-phs residual and analytic-tangent interface.
  Admission-sized coupling asymmetry retains its small energy-gradient
  remainder as well.
- `assemble(storage, zetas, strike_weights)` — the pHS: symplectic
  pair blocks, per-mode viscous damping `R = diag(0, 2 zeta w)`
  (zetas come from the caller; the visco-damping facility's per-mode
  output slots in here — no second damping representation invented),
  one strike port. REFUSES asymmetric couplings, negative
  coefficients/damping, non-positive/non-finite frequencies, and
  non-finite couplings by name (NaN-proof negated comparisons);
  duplicate modes in either constructor list are refused.
- `assemble_storage_with_damping` retains the full symmetric positive-
  semidefinite modal viscosity C [1/s] in the momentum block of R, using
  the same storage, signed force ports and fs-phs stepping as the diagonal
  `assemble_storage` path. All supplied coordinates must use the same
  mass-normalized basis. The full matrix owns all selected viscous loss;
  no separate damping ratios are added. Off-diagonal terms, exact zero rows
  and coupled lossless directions remain intact, without jitter, eigenvalue
  clipping or modal decoupling. Finite entries and port projections are
  admitted here; symmetry and PSD admission use the existing fs-phs owner
  and its numerical tolerances. Hereditary material loss still requires
  explicit constitutive states and is not represented by this fixed C.
- `von_karman_ss_plate` — simply-supported rectangle with ANALYTIC
  sine modes for displacement (mass-normalized) and Airy stress
  (unit-normalized; both are biharmonic eigenfunctions on the SS
  rectangle, in-plane movable edges). Channel coefficients
  `E h / (2 xi_j^4)`. Optional isotropic pretension `T` [N/m]
  raises the linear spectrum `ω² = (D k⁴ + T k²)/(ρ h)` without
  changing the Airy channels. The Hessian bracket is separated exactly into
  products of one-dimensional sine/sine/sine and sine/cosine/cosine
  integrals: O(Q) work per mode triple and order instead of a Q-by-Q grid.
  Reflection parity makes forbidden channels exactly zero. Signed
  off-diagonal interactions and mass/stress normalizations are retained.
  Coupling integrals use Gauss-Legendre quadrature
  whose order SCALES with the highest half-wave **index sum**
  `m+n` (not `max(m,n)` — a 0.20×0.15×2 mm panel with three odd–odd
  modes sat at 1.3e-8). Certified by a second
  independent order judged against max(channel scale, 1e-12 x global
  scale) — entrywise relative comparison falsely refuses
  analytically-zero entries, and an ALL-zero channel's own scale is
  pure roundoff (both executed); the residual is returned, not
  hidden. The two orders share the quadrature engine, so the
  certificate is a CONVERGENCE check, not an independent-route one. Stress-mode count is a SEPARATE truncation from the
  displacement count (both explicit inputs).
  G1 compares the separated kernel with independent pointwise two-dimensional
  area integration and checks the assembled high-mode tensor against analytic
  product-to-sum integrals, including both signs of off-diagonal coupling,
  symmetry, and forbidden-channel zeros. These verify the same retained
  continuum model; they do not close modal or Airy-basis truncation error.
- `von_karman_sampled_plate` / `SampledPlateMode` — same Airy
  channels from a regular sampled `w` grid. Second-order FD is the
  production stencil (FE `w` is C⁰). Second vs fourth Hessians
  must agree to `0.9` of the (dimensional) channel scale or
  construction refuses; a sine-grid reprint of the analytic Airy
  channel is the accuracy pin. Clamped and orthotropic bending
  enter as FE-sampled displacement; the membrane channel stays the
  analytic sine Airy basis (in-plane movable). Callers with
  orthotropic `E1 ≠ E2` pass `(E1+E2)/2` as `params.young`.
- `kirchhoff_carrier_string` — one diagonal channel
  `E[k,k] = (k pi/L)^2 * 2/(mu L)`, coefficient `E A L / 8`: exactly
  the averaged-tension Kirchhoff-Carrier form in mass-normalized
  coordinates (hand-derivation pinned in the battery). A guitar
  string is one filling; the law is a prestressed cable.
- `kirchhoff_carrier_moving_end` — the same channel on free-fixed
  cosines `κ = (k−½)π/L`. `φ(0) ≠ 0`, so the attachment is a
  Dirac port. Fixed-fixed sines stay the one-way bridge.
- `prestressed_beam_omega` — linear angular frequency of sine mode
  `k` on a uniform pinned-pinned Euler–Bernoulli beam,
  `ω² = κ² (T + EI κ²)/μ`, `κ = kπ/L`. `EI = 0` is the flexible
  taut-string limit; zero or vanishing tension retains the pure-bending
  `k²` spectrum without division by the preload. Subcritical compression
  lowers the spectrum, critical modal stiffness gives zero, and negative
  modal stiffness returns NaN instead of fabricating a stable frequency.
  Invalid dimensions, mode index, or bending stiffness also return NaN.
  Cables, stays, and strings share this law; it is not a free-free or
  clamped beam formula. G1 checks the bending and Euler-load limits;
  G3 checks mass/stiffness and geometric scaling, including representable
  frequencies whose squares lie outside the f64 range.
  The shared `prestressed_beam_omega_at_wavenumber` accepts each boundary
  image's spatial eigenvalue; the moving-end Galerkin image uses the same law
  at half-integer wavenumbers. For nonnegative tension, a scaled norm of
  tension and bending frequencies avoids squaring large frequencies or
  dividing by tension. Additional G0/G1 tests cover phase/group-speed ratios,
  limiting spectra and invalid data. The slender-beam law makes no claim for
  transverse shear, rotary inertia, postbuckling or nonlinear stability.
- `longitudinal::LongitudinalString` — the same moderate-slope axial
  material strain `u_x + (y_x² + z_x²)/2` with explicit finite longitudinal
  inertia. The caller supplies one or two cosine slope series, including
  the moving transverse chord, and strictly ordered longitudinal sine
  indices. Mass-normalized axial frequencies are
  `omega_j = (j pi/L) sqrt(EA/mu)`. Exact cosine products give the geometric
  mean `s0` and coefficients `sj`; storage is the positive sum
  `EA L s0²/2 + EA L/4 sum_retained(dj uj + sj)²`. Unretained zero-mean
  strain channels are statically equilibrated, without their inertia.
  `energy` includes retained axial quadratic stiffness; it excludes all
  kinetic energy and the pre-existing transverse linear potential.
  `discrete_gradient` gives exact finite-step work, without dividing by
  state differences. `coupling_discrete_gradient` removes only the axial
  linear gradient for a caller whose existing propagator already owns it.
  `strain_bounds` and `strain_at` expose actual reduced material strain
  for caller-owned material/slack admission. No method selects damping,
  advances state, or allocates during energy/force evaluation. Preparation
  admits at most 512 transverse modes per direction, two directions,
  and 32 explicit axial modes with indices at most 1024; it never drops
  a mode to satisfy a frequency cutoff. The caller owns the cutoff,
  resolved bandwidth, moving-end pullback, clock, and finite-strain domain.
- `longitudinal::relaxation::StrainRelaxation` retains each original positive
  Maxwell arm's mean/cosine relaxed strains. Its trial calls the canonical
  `fs_material::visco::PronyRamp` for a piecewise-linear strain history,
  preparing once per source arm and reusing coefficients across its
  mean/cosine coordinates. It returns exact endpoint memory, interval-mean
  stress and independently integrated nonnegative dashpot heat.
  `stage_with_arm_dissipation` performs that same update while retaining
  each original arm's physical heat for finite-body routing. Mean/cosine
  orthogonality weights enter both storage and dissipation.
  No midpoint relaxation or dt/tau <= 2 limit remains in this update.
  Work is mean stress times the endpoint strain increment; consumers pull it
  back through their quadratic strain map at midpoint coordinates. This
  preserves storage balance without treating endpoint-stress trapezoids as
  physical work. Exact history integration does not make nonlinear mechanics
  exact: the hosting coupled solve, source bandwidth, strain domain and other
  physical participants retain their time-resolution requirements. Accepted
  memory is never mutated by a trial. Zero/infinite rate use the
  elastic/fully-relaxed limits.
- `duffing_backbone` / `single_mode_beta` — the analytic
  perturbation pins.
- The `assemble*` systems retain fs-phs's original-energy discrete-gradient
  stepping, striking, damping, and work ledger.
- `sav::{SavStepper, SavConfig, SavState, SavRecord}` is a separate opt-in
  exact-linear, non-iterative stepper for the complete retained SOS potential.
  It introduces ONE `r=sqrt(V+C0)` for the total nonlinear potential; all
  stress channels, signed couplings and physical linear coefficients remain.
  The exact scalar `PreparedExactZoh` owner supplies held-force trajectories,
  same-interval mean velocities and independently integrated positive heat.
  A half-step predictor uses the current extended state, then one scalar
  Sherman-Morrison solve closes the held nonlinear force. Diagonal KC
  channels cost O(N); dense channels retain O(N^2) multiplication. Preparation
  and hot work grow with the retained modes without a dense 2N-state solve;
  warm stepping allocates nothing. Damping input is `C_i` [1/s], not zeta.
  The declared mode cap is at most 4096, and every retained omega*dt must
  remain within a positive caller band strictly below pi. Exact nonnegative
  displacement gains make the rank-one denominator at least one; invalid
  data or arithmetic refuses, without dropping modes or clipping a gain.
  This is a second-order exponential SAV method (Jiang/Wang/Cai, JCP 419,
  109690, doi:10.1016/j.jcp.2020.109690; scalar string formulation in
  Russo/Bilbao/Ducceschi, IFAC 2024). Its invariant is the MODIFIED energy
  `H_linear+r^2-C0`, not the original `H_linear+V`. Each accepted record
  reports both, their discrepancy, physical held-force work and independent
  viscous heat. The balance must close to 1e-12 of the physical/evolved energy
  scale; the arbitrary offset cannot inflate that allowance. The shifted
  auxiliary representation `r-sqrt(C0)` avoids losing weak nonlinear energy
  beside the offset. `SavState` contains the entire auxiliary checkpoint;
  reinitializing r from q during restart would change the method.

## Invariants

1. Quadrature certificate: two independent Gauss orders agree to
   1e-8 of the channel scale, or construction refuses
   (`QuadratureMismatch`).
2. Coupling symmetry is an admission requirement
   (`AsymmetricCoupling`), and the tensor is non-vacuous for the SS
   plate.
3. Duffing backbone: the measured amplitude-dependent frequency of a
   single-mode system tracks `w0 (1 + 3 beta a^2 / (8 w0^2))` within
   10% of the SHIFT at small amplitude.
4. Kirchhoff-Carrier exactness: `single_mode_beta` equals the
   hand-derived `E A (k pi/L)^4 / (2 mu^2 L)` to 1e-12, and the
   backbone reproduces the classical glide
   `dw/w = (3/32)(EA/T0)(k pi A/L)^2` in physical amplitude.
5. Energy conservation: undamped struck-level free run holds `H` to
   1e-9 relative over 20k steps.
6. ARCHITECTURAL FINDING (executed, recorded): under Gonzalez
   discrete-gradient stepping the correction term forces
   `dg.(x1-x0) = H(x1)-H(x0)` for ANY gradient function, so
   force-vs-energy divergence (the unsymmetrized-tensor bug class)
   CANNOT break energy conservation — unlike force-side integrators.
   The guards here are admission (2) and the TRAJECTORY: a divergent
   mutant conserves its coded `H` to ~1e-7 while departing the true
   dynamics by > 5% — both pinned.
7. Amplitude-scaling covariance: the nonlinear frequency shift
   scales as amplitude^2 (measured ratio 4 +/- 10%).
8. Mode coupling cross-check: the two-mode tension-modulated string
   transfers energy into a seeded second mode (weak, second-order
   parametric — correct string physics) and the discrete-gradient
   trajectory agrees with an independent in-test RK4 integration of
   the same vector field to 2%.
9. Cascade onset: a fundamental-only initial condition leaks energy
   into other modes only through the nonlinearity; the leaked
   fraction grows strongly with `w/h` (0.1 -> 0.6 -> 2.0), ledger
   closed throughout.

## Error model

Typed `NlModalError` (parameter, asymmetric coupling, quadrature
mismatch, pHS passthrough). No silent degradation.
The opt-in stepper adds `SavError` for incompatible state/inputs, the exact
mode and angular-step band, physical linear-owner refusal, modified-energy
balance failure, and cancellation. Refused steps publish no physical or
auxiliary state; retry and checkpoint replay retain the same trajectory.

## Determinism class

Deterministic: Newton-on-recurrence Gauss nodes, fixed orders, fs-phs
stepping; no RNG or time dependence.

## Cancellation behavior

Construction is milliseconds at analytic-mode scale; stepping is the
long pole (implicit Newton per step). No `Cx` integration (workspace
`frankensim-ccmn`).
`SavStepper::step_controlled` polls the supplied cancellation callback between
modal operations and immediately before publication. Its state and auxiliary
coordinate are transactional; scratch may change on refusal.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/sav.rs` independently checks exact 1/5/10 kHz linear phase at 48 kHz,
bitwise equality with the existing damped/forced exact owner when V is zero,
unchanged nonlinear static compliance, and a forced/damped modified-energy
ledger. A signed two-channel model transfers energy into an initially silent
mode and converges to an independent continuous RK4 oracle; both state error
and physical-versus-auxiliary energy defect refine at second order. An AGM
elliptic Duffing oracle checks a 5 kHz mode at nonlinear ratios .01, .1 and .5,
with 1/2/4/8/16 substeps and finer asymptotic-order rungs. At 16 substeps the
frequency error is below 0.1 cent for these fixtures; this is not an automatic
production clock or a wider nonlinear-band guarantee. The suite also retains
all 128 KC modes with one rank-one solve, tests band/capacity/offset refusals,
mid-step cancellation, unchanged rejected state, and bitwise checkpoint replay.

`tests/nlmodal.rs::full_viscous_operator_preserves_local_damper_nullspace_and_forced_response`
checks an independent common/differential analytic oscillator solution,
rank-one damping's lossless direction and physical port work to 3e-12.
`full_damping_diagonal_limit_is_the_same_system` pins identical state and
work for the old diagonal image. `full_damping_nonlinear_trajectory_converges_and_closes_work`
checks the coupled quartic trajectory against independent RK4 on three time
rungs (second-order refinement; finest state error <3e-5) with total
work/energy closure <2e-10. Active, asymmetric and nonfinite operators
refuse; zero loss remains exactly zero. These are G0/G1/G3 numerical model
checks, not experimental material-damping validation.

The nonlinear modal battery also checks a resolved quartic force beside an
unchanged high-energy mode against an independently factored scalar divided
difference, and coupled-channel finite-step work against two-node Gauss
integration of the ordinary cubic gradient, including time reversal and
the zero-increment limit.

`src/longitudinal/relaxation.rs` tests compare original-arm heat to independent
spatial integrals, retain identical ordinary updates, verify work/storage
balance and the second-order relaxation-time ladder, and check refusal/retry
without mutating accepted memory.

`src/longitudinal_tests.rs`: independent sampled-span strain integration
including statically eliminated channels; ordinary force versus energy
finite differences; exact finite-step work and removal of only the axial
linear force; the existing Kirchhoff--Carrier condensed limit; analytic
single-sine doubling of the axial drive index; moving-chord reciprocal
pullback and whole-span strain bounds; finite and declared-cap admission.
These are G0/G1/G3 model tests, not measured piano calibration.

`src/longitudinal/relaxation.rs` checks independent spatial Maxwell stress
and Simpson-integrated heat; exact fixed-strain decay under 1/4/16/64
partitions including stiff material times; linear-strain work, heat and memory
partition invariance; and transactional repeated/invalid trials. The immediate
string and fiber-ribbon consumers check their actual nonlinear geometric
pullbacks against independent spatial and temporal integrations. These G0/G1/G3
tests establish the coded constitutive law and work semantics, not experimental
loss calibration or unrestricted nonlinear time-step accuracy.

`tests/nlmodal.rs` (10): construction certificates (quadrature
residual, mode counts, assembly); non-vacuous coupling + asymmetry
refusal; Duffing backbone vs perturbation formula; Kirchhoff-Carrier
beta + glide hand formulas; energy conservation; the
unsymmetrized-tensor architectural falsifier (conservation survives,
trajectory diverges); amplitude-scaling covariance; parametric
coupling + RK4 cross-check; plate cascade casebook (JSON lines with
modal energies, leaked fractions, tensor size); typed refusals.

## No-claim boundaries

- The opt-in SAV API does not migrate existing render consumers or close the
  full M1.1-W3/M5.2-W1 work items. Conserving auxiliary energy alone certifies
  neither physical nonlinear energy nor accurate harmonics. The caller must
  refine temporal resolution and inspect state, pitch and auxiliary-energy
  defect; nonlinear harmonics can exceed the retained linear band. Offset
  sensitivity, modal/stress truncation and physical model error remain separate
  convergence axes. Full matrix viscosity, contact, constitutive memory and
  reciprocal termination coupling remain with their existing pHS owners.
- Moderate-rotation von Karman regime only: no damage, plasticity,
  wrinkling, or full shell curvature (quadratic coupling from
  curvature is the shells follow-up).
- Simply-supported in-plane-movable boundary in v1; FE-mode input is
  the follow-up (trigger: consumers with non-analytic geometry —
  which is also when coupling-tensor CACHING starts to matter; at
  analytic-mode scale the tensor build is milliseconds and a cache
  manifest would be ceremony, deliberately not built).
- The Airy potential form `V = (Eh/8) sum S_j^2/xi_j^4` is taken
  from the standard modal literature (Ducceschi-Touze class); its
  validation here is through the Duffing/KC analytic pins and the
  structural energy identities, not a from-scratch continuum
  re-derivation.
- Literature-value comparisons for internal-resonance energy
  exchange are deferred (measured tables need verified
  transcription); the executed cross-check is an independent
  INTEGRATOR over the same coded vector field — it validates the
  stepper, not the tensor physics (review-corrected wording).
- Viscous ratios or full damping matrices are supplied by the caller.
  Frequency-dependent loss remains the shared material facility's;
  the longitudinal history adapter consumes its actual source Maxwell arms.
- The finite longitudinal string image fixes longitudinal endpoint
  displacements while allowing caller-supplied transverse chord motion.
  Omitted axial coordinates are statically condensed, not resolved waves;
  this omits their inertial correction and does not certify broadband
  convergence. No tension-dependent winding slip, axial bridge compliance,
  full finite-strain geometry, or identified material loss is implied.
