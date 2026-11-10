# fs-nlmodal — Contract

## Purpose and layer

Layer L3. Geometric nonlinearity for thin structures in modal
coordinates (bead frankensim-fsim-vonkarman-u5nme): von Karman plates
and Kirchhoff-Carrier strings as port-Hamiltonian systems with a
quartic SUM-OF-SQUARES potential
`H = 1/2 sum (p^2 + w^2 q^2) + 1/4 sum_j c_j (q^T E_j q)^2` — the
Airy-stress-eliminated modal form. Pitch glide, mode coupling, and
the crash cascade with structurally exact energy accounting. The linear
case is also a composition-kernel Component (`component`), with a joint
loaded basis for modes coupled through shared air (`loaded`), for which the
crate depends on fs-modal (the exact damped-mode step and the dense pencil
solve), fs-ports and fs-exec.

## Public types and semantics

- `SosModalStorage` (implements fs-phs `Storage`): frequencies +
  `StressChannel`s (`c_j >= 0`, symmetric `E_j`). Gradient is the
  EXACT analytic gradient of the coded `H`; state layout matches
  fs-phs `modal_bank` (`[q, p]` interleaved).
  `validate()` owns frequency/channel admission for `assemble` and
  consumers that retain this storage inside multiport/full-viscosity systems.
  Exactly symmetric channels own the complete average-vector-field discrete
  gradient through the existing fs-phs constitutive interface. With
  `d=q1-q0` and midpoint `m`, each quartic channel integrates its original
  cubic force exactly as
  `c [(m^T E m + d^T E d/12) E m + (m^T E d/6) E d]`.
  The kinetic efforts stay equal to the mean momenta; a channel's null
  directions and unrelated modes receive no artificial elastic correction.
  The source also supplies the analytic endpoint action
  `integral_0^1 t Hessian(H)(a+t(b-a)) v dt`, distinct from the continuous
  tangent. Existing direct-sum composition and prepared analytic Newton use
  this action without a separate model, numerical quadrature, finite
  differences, dense Hessian construction, or evaluation allocation.
  Both force and endpoint action cost O(stress count × mode count²) with
  fixed traversal and constant temporary memory. Public model data and
  endpoint extents are checked; invalid/nonfinite arithmetic refuses.
  Admission-sized asymmetry retains the existing Gonzalez path, because
  its coded effort is not exactly the Hamiltonian gradient. No source
  coefficients are silently symmetrized to manufacture constitutive support.
  Independent energy accounting retains the complete midpoint remainder.
  Let `o=(q0-m)+(q1-m)`, `D=d^T E d`, `O=o^T E o`, `C=d^T E o`,
  `M=m^T E d`, `N=m^T E o`, and `S=m^T E m`. Each symmetric quartic
  channel contributes `c [(M+C/2)(N+(D+O)/4)+S C/2]`; quadratic modal
  terms contribute their represented-midpoint offsets as well. At an exact
  midpoint this reduces to `c D M/4`, with zero quadratic remainder.
  This retains small resolved work beside large stored modal energy without
  subtracting total endpoint energies. It uses one triangle traversal per
  stress channel and no per-evaluation allocation. Admission-sized coupling
  asymmetry retains its small energy-gradient remainder as well.
  `hessian_action(state, direction, out)` supplies the exact continuous
  directional derivative of this same effort law for harmonic-balance and
  variational solvers. Each symmetric channel contributes
  `c [(q^T E q) E v + 2 (q^T E v) E q]`, plus the original linear modal
  stiffness and unit momentum derivative. The admitted roundoff asymmetry is
  retained through `v^T E q + q^T E v`; no symmetrization changes the coded
  force law. Evaluation uses fixed-order matrix-vector contractions, constant
  temporary memory, and no allocation, dense Hessian or finite differences.
  It validates public model data and input extents, refuses nonfinite
  arithmetic, and costs O(stress count × mode count²). G1 tests compare the
  full tangent with independent effort derivatives for native string/plate
  models and scalar/cross-mode analytic stiffness; G3 checks reciprocity and
  preservation of admitted asymmetry. The action is not a discrete-gradient
  step Jacobian, a harmonic truncation bound, or an orbit-stability claim.
  The same action is exposed through `Storage::continuous_tangent`, including
  after assembly into `Box<dyn Storage>` and supported pHS composition. Ordinary
  harmonic-balance and variational adapters can therefore differentiate the
  actual assembled law without maintaining a second force model or callback.
- `PreloadedModalStorage::new(source, equilibrium_displacements)` translates
  the SAME conservative SOS law under the required held generalized force
  `f0 = grad V(a)`. Incremental storage is
  `p^2/2 + V(a+q)-V(a)-f0.q`. Each source channel retains its quadratic,
  cubic and quartic increments through the factored form
  `c (s*u/2 + t^2 + t*u + u^2/4)`, with `s=a^T E a`, `t=(E a)^T q`,
  and `u=q^T E q`. Small motion is not lost by subtracting nearly equal
  source energies or forces. The complete local stiffness
  `diag(omega^2) + sum c [s E + 2 (E a)(E a)^T]` includes its signed
  off-diagonal entries; existing pHS PSD admission refuses negative tangent
  curvature without altering it. Semidefinite curvature does not establish
  local stability: cubic terms can destabilize a null direction. Exact stress symmetry is required rather
  than silently repairing the source. `preload_forces()` is already
  accounted in the translated potential; incremental applied drive excludes
  that held force. `Storage::quadratic_part` supplies this full local stiffness
  to existing exact-linear nonlinear stepping, and the continuous tangent
  retains every nonlinear interaction. The same source owns a complete shifted
  average-vector-field discrete gradient and its analytic endpoint action.
  With cached `r=E a`, `s=a^T E a`, incremental midpoint `m`, difference `d`,
  `t=r^T m`, `u=m^T E m`, `D=d^T E d` and `Z=r^T d+m^T E d`, the original
  channel's mean restoring force is
  `c [(2t+u+D/12) r + (s+2t+u+D/12) E m + (Z/6) E d]`.
  It includes the preload-induced quadratic/cubic interactions and cancels
  the held dead load algebraically. It never reconstructs `a+q`, subtracts
  large absolute forces, or transfers elastic work into momentum efforts.
  The endpoint derivative differentiates this same complete polynomial;
  coincidence gives half the full continuous Hessian. Both capabilities
  propagate through the existing boxed/composed pHS interface and prepared
  analytic Newton. The immutable source admission remains the owner of exact
  stress symmetry and local-curvature refusal. Evaluations reuse diagonal
  channel metadata, cost O(modes) per diagonal channel or O(modes²) per dense
  channel, and allocate no temporary state or Hessian. The stable midpoint remainder retains
  the represented midpoint offset `o=(left-m)+(right-m)`, including quadratic
  kinetic and elastic terms. For each channel, with `z=a+m`, it evaluates
  `c [(z^T E d + d^T E o/2)(z^T E o + (d^T E d + o^T E o)/4)
  + (z^T E z)(d^T E o)/2]` through compensated factored contractions of
  source-reference coefficients, increments and offsets; it never forms `z`
  or subtracts endpoint energies. At an exact midpoint this reduces to
  `c/4 (d^T E d) [(E a)^T d + m^T E d]`. Warm calls allocate nothing;
  diagonal stress channels cost O(modes), dense channels O(modes^2).
  This is a caller-declared static displacement under a dead load, not a
  pressure projection, equilibrium solve, curved-shell theory, local/global
  stability certificate, or permission to ignore generated-harmonic bandwidth.
  `tests/preload.rs` checks the exact biased Duffing law and its DC/even/odd
  force harmonics, full coupled translation and tangent, tiny resolved motion
  beside a large static reference and kinetic background, independent
  Gauss-integrated force components/work, endpoint polynomial derivatives,
  native/prepared mechanics against an independent scalar solve with separate
  spring work, external work and physical viscous heat, the
  centered limit, represented-midpoint offsets, cubic instability in an
  admitted Hessian nullspace, and negative-curvature/asymmetric refusals.
- `PreloadedModalStorage::invariant_axis_spectrum` derives the conservative
  orbit about that held preload from the original cached tangent and stress
  contractions. Every active channel must preserve the chosen coordinate,
  including preload-induced cross-mode forces; no coupling threshold or
  cross-channel cancellation manufactures an invariant axis. The scalar
  potential must be globally convex and confining. The query retains its
  unequal turning points, physical frequency, mean incremental displacement,
  every cosine order through H, signed physical velocity sine amplitudes and
  acceleration cosine amplitudes. The static reference is not added twice.
  A turning-point substitution cancels the energy integral's endpoint square
  roots before composite Simpson integration of time and Fourier coefficients.
  Two consecutive grid refinements must resolve both relative frequency and
  the normalized, h-squared-weighted coefficient changes through order 2H.
  Caller budgets bound work and memory; exhausted budgets refuse. Symmetric
  and zero-amplitude limits reuse the exact Jacobi kinematic owner. Binary
  product scaling preserves representable physical derivatives without an
  intermediate squared-frequency overflow. The measured H+1..2H guard-band
  norms are finite omitted content, NOT an infinite-tail bound; refinement
  estimates are not interval enclosures. This is an opt-in source/solver seed,
  not a change to coupled dynamics, damping, rendering or stability admission.
  `tests/preloaded_spectrum.rs` compares the waveform, period, DC/even/odd
  coefficients and velocity with independent RK4 integration of the original
  `Storage::gradient`; it checks original energy/force, preload-reflection
  parity, symmetric limits, extreme physical-unit scaling, finite guard-band
  refinement and exact coupling/work-budget refusals.
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
  changing the Airy channels. The Hessian bracket is separated into exact
  integer trigonometric products. For odd `s+a+b`, the triple-sine integral
  on [0,1] is `-4 s a b / (pi (s+a+b)(s+a-b)(s-a+b)(s-a-b))`; for even
  sums it is zero. The sine/cosine/cosine integral follows by multiplying
  this result by `(a²+b²-s²)/(2ab)`. These are product-to-sum identities,
  not fitted coefficients. Signed integer differences preserve near-resonant
  overlaps; the combined normal/shear polynomial uses exact i128 arithmetic
  when representable and f64 otherwise. Reflection parity makes forbidden
  channels exactly zero. Signed off-diagonal interactions and mass/stress
  normalizations are retained. Assembly costs O(stress count × displacement
  count²), with no quadrature order increasing with the highest mode index.
  Each index must be exactly representable in f64 (at most 2^53), and
  nonfinite material, geometry, or assembled storage is refused.
  `quadrature_residual` is zero for this analytic constructor; that denotes
  absence of quadrature truncation, not a bound on all numerical errors.
  Stress-mode count remains a SEPARATE truncation from the displacement
  count (both explicit inputs).
  G1 compares the exact kernel with independent pointwise two-dimensional
  area integration at two Gauss orders and checks the assembled tensor against
  direct product-to-sum integration. Regressions retain far-separated small
  overlaps, near-resonance, signed weak shear/normal cancellation, exact parity
  zeros, permutation symmetry, and the known odd-integer scaling of a sparse
  high-overtone model. These verify the same retained continuum model; they
  do not close modal or Airy-basis truncation, time-step convergence, or the
  material's physical high-frequency validity band.
- `von_karman_sampled_plate` / `SampledPlateMode` — same Airy
  channels from a regular sampled `w` grid. Second-order FD is the
  production stencil (FE `w` is C⁰). All interior nodes participate,
  including the first ring next to the boundary; the Airy trace is zero
  on the boundary itself. Displacement Hessians are prepared once and
  reused across pairs and Airy channels. Airy indices must be strictly
  below the grid Nyquist limit to prevent a zero or aliased stress field;
  this admission does not certify the resolution of displacement products.
  Second vs fourth Hessians
  must agree to `0.9` of the (dimensional) channel scale or
  construction refuses; signed sine-mode interactions must converge
  quadratically to the analytic Airy overlaps under grid refinement.
  Nonfinite inputs and assembled storage are refused. Clamped bending
  enters as FE-sampled displacement; this constructor retains isotropic
  membrane elasticity and the analytic sine Airy basis (in-plane movable).
- `von_karman_sampled_plate_with_membrane` / `VkSampledPlateParams` —
  the same sampled kinematics with an admitted uniform `SheetMatrix` [N/m]
  and structural mass per area. The complete plane-stress tensor supplies
  the complementary Airy matrix `C_jk = integral b_j^T A^-1 b_k`, with
  `b_j = [Psi_j,yy, Psi_j,xx, -Psi_j,xy]`. Exact sine/cosine integrals
  retain the off-diagonal Airy coupling from rotated normal/shear compliance.
  Construction factors `C = L L^T` and transforms each geometric coupling
  by `L^-1`, so `Vnl = 1/8 e^T C^-1 e` remains positive sum-of-squares
  storage with no time-step matrix solve or added evaluation allocation.
  The isotropic limit is exactly `C_jj = xi_j^4/(E h)`. Independent
  area integration in material axes checks rotated compliance; hand
  orthotropic formulas check shear/Poisson dependence, and covariance
  tests check stiffness and mass normalization. Singular, nonfinite or
  numerically unresolved compliance refuses. This is a uniform membrane
  model; it does not add spatially varying material, laminate extension/
  bending coupling, nonlinear constitutive elasticity, or in-plane inertia.
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
- `longitudinal::LongitudinalModalStorage` supplies the complete conservative
  fixed-position modal storage for that SAME finite-axial-inertia strain law.
  Its interleaved state retains all transverse modes plane-major, followed by
  the authored axial indices. One or two physical bending rigidities select
  transverse planes; their linear frequencies come from T, EI, length and
  solid mass, while axial frequencies come from EA and the same mass. The
  mass-normalized sine basis gives a reciprocal slope/force pullback. Nonzero
  EI means pinned bending supports, not clamped endpoints. Exact cosine
  products retain spatial sum/difference coupling without runtime quadrature;
  no requested axial mode is removed to satisfy a time-step cutoff.
  Original kinetic energy, tension/bending storage, geometric mean strain and
  every retained nonuniform strain channel contribute to one Hamiltonian.
  The source's one-state strain contractions are reused once per force call.
  The analytic continuous Hessian differentiates both strain and its geometric
  pullback, including reciprocal axial/transverse and cross-polarization
  terms. The original exact two-point strain gradient and its analytic finite
  endpoint Jacobian also survive the generic Storage interface; a midpoint
  Hessian never replaces that endpoint derivative. Warm callbacks allocate
  nothing and use bounded stack scratch, with O(planes * (transverse modes + 1) *
  (axial_modes + 1)) work. Generic pHS assembly, harmonic balance, variational
  and implicit mechanics can consume this complete source directly. Damping,
  material memory, moving supports and generated-harmonic bandwidth remain
  caller-owned; unretained axial channels are still statically equilibrated.
  `tests/longitudinal_modal.rs` checks independent resolved spatial energy and
  force quadrature, original KC condensation, complete tangent/endpoint
  derivatives and work, a circular two-polarization orbit with axial DC, and
  a sum-frequency axial resonance with three-rung temporal refinement.
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
- `duffing_backbone_exact` evaluates the full conservative Duffing period by
  the energy-integral/AGM identity, using the original physical turning-point
  amplitude. Hardening, pure-quartic motion and softening below the separatrix
  share one bounded, allocation-free owner. Binary scaling and fs-math
  quad-double products retain extreme units and near-separatrix cancellation;
  the final AGM uses deterministic square roots. Invalid/nonperiodic inputs
  and unrepresentable frequencies refuse. The old perturbation function keeps
  its small-amplitude meaning. `SosModalStorage::invariant_axis_backbone`
  derives the frequency from the actual source coefficients, requiring each
  active stress channel to leave the selected axis invariant. It never drops
  coupling to manufacture a scalar resonance. These are conservative-orbit
  frequencies, useful as periodic-solver seeds, not forced peaks or transverse
  stability results; turning amplitude differs from the pinned Fourier
  amplitude in harmonic balance. G1/G3 tests independently integrate the energy
  period and check the weak/pure-quartic/separatrix and unit-scaling limits.
- `duffing_backbone_directional_derivative` supplies the analytic local
  frequency change for any direction in the original linear frequency,
  cubic stiffness and turning amplitude. It differentiates the full-period
  AGM, including a nonlinear-stiffness direction at a linear baseline and the
  divergent softening slope near the separatrix. Independent binary scaling
  of source and derivative products retains extreme units; quad-double
  contractions retain cancellation. No finite-difference step or
  small-amplitude approximation is used. Nonperiodic sources and
  unrepresentable derivatives refuse. The storage-bound
  `invariant_axis_backbone_amplitude_derivative` reuses the same original
  source reduction and channel-by-channel invariance admission. G1/G3 tests
  compare an independently differentiated energy-period integral and verify
  linear, pure-quartic, separatrix, homogeneity and unit-scaling limits.
- `SosModalStorage::gradient_parameter_derivative` differentiates the complete
  original effort at fixed mass-normalized modal state for a supplied direction
  in the original angular frequencies and every stress channel's coefficient
  and coupling matrix. It retains signed/off-diagonal coupling and the cubic
  derivative of a channel whose baseline coefficient is zero. The exact
  channel chain rule uses two matrix-vector passes and constant scratch; no
  perturbed model or finite-difference step is formed. This is the direct
  material term, distinct from the state Hessian action; periodic response and
  physical output derivatives need both. Source/direction admission and finite
  arithmetic remain mandatory, and a changing modal basis or mass normalization
  needs an additional coordinate chain rule. Tests compare independently
  perturbed original gradients, verify the cubic small-amplitude and exact
  linear limits, and refuse incomplete, asymmetric or nonfinite directions.
- `duffing_spectrum_exact` / `DuffingSpectrum` retain an explicit count of odd
  cosine harmonics of that same conservative orbit. The Jacobi `cn` (hardening)
  and `cd` (softening) nome series retain physical amplitudes and the alternating
  softening signs. DC and even orders are zero by the original equation's
  symmetry. The fundamental is not silently equated to the turning amplitude.
  A separately converged infinite-series normalization is independent of the
  requested retained order: refinement preserves every earlier coefficient.
  `displacement_tail_bound` numerically evaluates an analytic absolute-sum
  bound on the omitted displacement at all phases, excluding floating-point
  evaluation error; it is not an interval certificate or a coupled-model error
  bound. The work/storage cap is 4096 odd harmonics. Binary scaling retains
  extreme physical units; split nome multiplication preserves representable
  harmonics even when the dimensionless nome underflows. Unrepresentable
  coefficients, moduli or bounds refuse. `SosModalStorage::invariant_axis_spectrum`
  uses all original stress channels and the same invariant-axis admission as
  the frequency query; coupled axes refuse. These are explicit conservative
  references and nonlinear solver seeds, not changes to any render path.
  G1/G3 tests compare full waveforms and harmonics against independent RK4
  trajectories using an energy-quadrature period, check weak-nonlinearity and
  quartic/softening limits, truncation/refinement and physical-unit scaling.
- `duffing_kinematic_spectrum_exact` / `DuffingKinematicSpectrum` extend the
  same prepared Jacobi series with signed physical velocity sine coefficients
  and acceleration cosine coefficients. Scalar frequency/nome preparation is
  shared once, and displacement remains the existing spectrum. Binary-scaled
  source amplitudes retain derivative harmonics even when their displacement
  rounds below the subnormal range; no intermediate squared frequency is
  required. Closed weighted geometric sums report separate uniform omitted
  velocity and acceleration bounds in their original physical units. These
  estimate harmonic truncation, excluding coefficient roundoff; they are not
  interval enclosures. A displacement tolerance alone does not establish
  derivative resolution. This opt-in query refuses unrepresentable derivatives
  or bounds without narrowing the displacement query. The storage adapter
  `invariant_axis_kinematic_spectrum` preserves all source channels and the
  same coupling refusal. Original-ODE, source-force, truncation and extreme
  physical-scale tests cover the complete interleaved kinematic seed.
- The `assemble*` systems retain fs-phs's original-energy discrete-gradient
  stepping, striking, damping, and work ledger.
- `sav::{SavStepper, SavConfig, SavState, SavRecord}` is a separate opt-in
  exact-linear, non-iterative stepper for the complete retained SOS potential.
  It introduces ONE `r=sqrt(V+C0)` for the total nonlinear potential; all
  stress channels, signed couplings and physical linear coefficients remain.
  The existing `PreparedExactZoh` owner supplies held-force trajectories,
  same-interval mean velocities and independently integrated positive heat.
  A half-step predictor uses the current extended state, then one scalar
  Sherman-Morrison solve closes the held nonlinear force. Diagonal KC
  channels cost O(N); dense channels retain O(N^2) multiplication. Diagonal
  viscosity uses one scalar exact cache per mode and O(N) linear work, with
  a declared mode cap of at most 4096. `new_with_damping_matrix` instead
  admits the complete symmetric PSD `C` [1/s] in the original mass-normalized
  coordinates, retaining cross-mode losses and their lossless directions.
  Its complete linear state and independent force ports use the same exact
  dense owner; bounded cold preparation is followed by allocation-free
  O(N^2) linear warm work. The existing augmented-exponential cap admits
  at most `MAX_SAV_DAMPING_MATRIX_MODES = 42` non-diagonally damped modes;
  excess work refuses before preparation. An exactly diagonal matrix reuses
  the scalar path bitwise and keeps its larger capacity. Both paths reuse all
  buffers, retain one nonlinear scalar solve, and leave physical stiffness
  unchanged. Damping input is viscosity C, not zeta. Every retained omega*dt
  must remain within a positive caller band strictly below pi. The complete
  held-force displacement map B enters `1 + g^T B g / 4`; passivity makes
  this denominator at least one. Invalid data or arithmetic refuses, without
  dropping modes, replacing C by its diagonal, or clipping a gain.
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

### `component`: the linear modal structure on the composition kernel

`component::ModalStructure` is the linear case of the modal structures above (no quartic term)
as an `fs_ports::component::Component` of contract v0. It wraps an admitted
`fs_modal::time::ModalAcousticTimeModel` and takes that owner's exact held-force step,
`step_duration_into`, for every proposal: nothing is integrated here. `from_modes` admits the
owner's model from `(angular frequency, viscous damping ratio)` pairs with no observer
transfer; `new` takes an admitted model and refuses one that carries an observer transfer,
because an observer is a probe of its own.

- Ports, in order: one power port per `ModalAttachment`, receiving effort, then `heat`, a heat
  port receiving temperature. Row `r` of an attachment has one weight `s_rk` per mode: the
  generalized force on mode `k` is `sum_r s_rk e_r` and the port returns
  `sum_k s_rk (q1_k - q0_k) / dt`, the mean flow of the step. Point forces on a structure are
  rows of mode shapes at the points, a mouth pressure on an air column is a row of `1/sqrt(M_k)`
  (`M_k` the mode's acoustic mass per unit volume flow [kg/m^4]), and identity weights are the
  modal port. The heat port returns
  `-D / dt`, where `D` is the owner's viscous heat of the step.
- State: `Vec<ModalAcousticState>`, one entry per mode, owned by the system.
  `stored_energy_j` is `ModalAcousticTimeModel::energy_j`, the owner's per-mode law and
  compensated sum.
- Receipt: stored change `E1 - E0` by that function; dissipation `D`; numerical loss
  `W - D - (E1 - E0)`, the owner's work-heat residual, which the owner itself refuses above its
  allowance; entropy `D / T`. The energy tolerance covers the kernel's port products against
  the owner's work (sums reordered, endpoint differences rounding at the size of the endpoints),
  the heat port's division and multiplication by `dt`, and the subtractions, at
  `4 (rows + modes + 16)` units in the last place of their sizes. The owner's heat allowance
  is not in it: the same `D` enters the audit through the heat port and the numerical loss, so
  its error cancels.
- Class `Smooth`: the mechanical rows are affine in the efforts, the heat row is quadratic. No
  tangent tier.
- `ModalStructureRefusal`: the owner's error unchanged, an observer transfer, weights of the
  wrong count or not finite, a refused identifier or port list, a state of the wrong length, a
  temperature that is not positive and finite, cancellation.

### `loaded`: the joint loaded modal basis

`loaded::LoadedModalBasis::new(w, M_a)` takes `n` modes, each in its own
mass-normalized coordinates (unit modal mass, so `K = diag(w_k^2)` and `M = I`),
and the symmetric added mass `M_a` of the air they share, in the same
coordinates. It solves them as one pencil `(K, M + M_a)` with fs-modal's dense
solve (`eigh_gen_dense`: Cholesky reduction, then Jacobi) and refuses
(`LoadedMass`) unless `M + M_a` is positive definite. It returns:
- the loaded angular frequencies, ascending;
- fs-modal's certified residual of each `w'^2`;
- the basis `Phi`, with `Phi^T (M + M_a) Phi = I` and
  `Phi^T K Phi = diag(w'^2)`; the unloaded coordinates are `x = Phi y`.

Methods:
- `congruent(A)` is `Phi^T A Phi` for a symmetric `A`, symmetric bit for bit.
  The loss `C` is one such matrix: each mode's own loss on the diagonal, the
  radiation resistance between modes off it.
- `loaded_attachment` carries an attachment written in the unloaded
  coordinates to the loaded ones (`s' = s Phi`) and keeps each row's flow;
  `unloaded(y)` is `Phi y`.
- `structure(rate, budget, modal, kind, attachments)` builds the loaded modes
  as one undamped `ModalStructure`. Its ports are the modal port first (width
  `n`, identity weights), then the carried attachments, then `heat`.
  fs-system's `elements::Resistance` with `Phi^T C Phi` joins the modal port by
  a common-flow junction, so the definiteness of the whole loss is checked on
  the sum.
- `held_loss_factors(dt)` is `sinc^2(w'_j dt / 2)` for each loaded mode. The
  structure steps exactly under a force held over the step, and the resistor
  answers the step's mean velocity. So for one underdamped mode with loss `c`
  the step's determinant is `(1 - e) / (1 + e)`, with
  `e = (c dt / 2) sinc^2(w dt / 2)`: the energy envelope decays by that factor
  per step, while within one step the energy ratio depends on the phase,
  between its square and 1. The envelope decays at the factor times `c` to
  first order in `e`. The factor is stated with a render, not corrected; it is
  NaN for a step that is not positive and finite.

Refusals, all on `ModalStructureRefusal`:
- no mode;
- a frequency that is not positive and finite (unloaded, or a loaded square);
- more modes than fs-la's dense symmetric eigensolver admits within its work
  cap, refused before the matrices are read;
- a matrix that is not `n x n`, finite and symmetric (a signed zero equals its
  mirror);
- a loaded mass that is not positive definite;
- the refusals of the carried attachments (each with its index in the list
  given, a weight that overflows when carried included) and of `from_modes`.

No claim: `M_a` and `C` are constants, so a frequency-dependent loading enters
at the frequency it was sampled at. The basis is dense, `O(n^3)`.

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

The full-viscosity SAV regressions independently resolve a local dashpot into
lossless common motion and damped relative motion. Closed-form forced and
unforced trajectories cover underdamping, critical damping, and overdamping
with damping*dt = 2600; state error stays below 5e-12 and accumulated physical
heat agrees within 2e-11 J. The lossless common motion generates less than
1e-23 J of numerical heat. With signed nonlinear stress channels and
noncommuting stiffness/viscosity, state, independently integrated heat and
auxiliary defect converge at second order to continuous RK4; finest weighted
state error is below 5e-6, heat error below 6e-6 J, and maximum auxiliary
defect below 2e-5 J. Modified energy/work/heat closes within 2e-11 J over the
run. Exact diagonal compatibility is bitwise; full matrix static compliance,
complete checkpoint replay, cancellation after candidate propagation,
nonfinite/overflowing drive and matrix/capacity refusals are also pinned.

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
difference, and every component of the coupled-channel finite-step effort
against two-node Gauss integration of the ordinary cubic gradient, including
time reversal and the zero-increment limit. Rank-one channels preserve their
null directions, uncoupled modal forces and physical momentum efforts beside
large kinetic storage. Five-point polynomial derivatives check the complete
endpoint tangent; coincidence recovers half the continuous Hessian. Actual
prepared analytic steps retain the kinematic velocity relation, independent
source work and viscous heat without a fallback tangent. Adjacent-float
midpoint tests pin both elastic and kinetic remainders. These G0/G1/G3 checks
concern the retained model and stepping equations, not unrestricted timestep
accuracy, convergence of the mode/stress truncation, or experimental fidelity.

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

`tests/modal_structure.rs` refuses a model that carries an observer transfer and admits the
same model without it (an integration test: the crate's lib test target does not compile at
origin, where `fd_probe` imports `fd_coupling`, removed in ca7b47dbf8). The composition-kernel slice of `ModalStructure` is
`crates/fs-system/tests/modal_structure_slice.rs` (`--features modal-slices`); see the
fs-system contract.

`tests/loaded_basis.rs` checks `LoadedModalBasis` against cases solved by hand:
- one mode under added mass `mu` sounds at `w0 / sqrt(1 + mu)` to 4 eps, for four
  `(f0, mu)`;
- two equal modes split by a mutual mass `m` sound at `w0^2 / (1 + mu + m)` and
  `w0^2 / (1 + mu - m)`, with in-phase and out-of-phase vectors;
- two unequal modes sit at the roots of their pencil's quadratic (its
  discriminant formed free of cancellation), within fs-modal's certificate
  widened by the roundoff of the residual itself;
- a three-mode basis is normalized in `M + M_a` and diagonalizes `K`, within
  bounds derived from fs-modal's own mass normalization and, off the diagonal,
  from its certified residuals over the gap between modes (the identity
  `(lambda_i - lambda_j) phi_i^T M phi_j = phi_i^T r_j - phi_j^T r_i`), each
  widened by the rounding of the products that form it;
- `held_loss_factors` is `sinc^2(w' dt / 2)` bit for bit and NaN for a step
  that is not positive and finite;
- a carried attachment keeps its rows' flows, within the rounding of two nested
  inner products on each side;
- the structure has its ports and runs the loaded modes undamped; one step
  from rest shows the modal port moving each loaded mode alone, and the strike
  row pushing and reading each through its carried weight;
- the refusals, among them a signed zero admitted as its own mirror, a basis
  beyond fs-la's work cap, and a weight that overflows when carried.

Its composition-kernel slice is `crates/fs-system/tests/loaded_modal_slice.rs`
(`--features modal-slices`).

## No-claim boundaries

- `ModalStructure` offers no port tangent, clones the owner's model and allocates its
  workspace in every proposal, and claims nothing the owner's step does not: exact for a held
  force, modes admitted below the Nyquist guard of the model's nominal clock, the owner's
  budgets. It adds no observer, radiation or contact law.
- Its port flows are endpoint differences `q1 - q0`; the owner books work on the increment it
  added and does not publish it, so the energy tolerance carries `|f| (|q0| + |q1|)` per mode.
  Under a statically held load or a free coordinate drifting under a constant force that
  tolerance is wide against the step's work, and a leak below it shows only in the ledger's
  accepted-defect account.
- The opt-in SAV API does not migrate existing render consumers or close the
  full M1.1-W3/M5.2-W1 work items. Conserving auxiliary energy alone certifies
  neither physical nonlinear energy nor accurate harmonics. The caller must
  refine temporal resolution and inspect state, pitch and auxiliary-energy
  defect; nonlinear harmonics can exceed the retained linear band. Offset
  sensitivity, modal/stress truncation and physical model error remain separate
  convergence axes. Contact, constitutive memory and reciprocal termination
  coupling remain with their existing pHS owners.
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
