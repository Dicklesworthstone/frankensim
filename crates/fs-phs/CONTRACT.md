# fs-phs — Contract

## Purpose and layer

Layer L2. Port-Hamiltonian systems (bead
frankensim-fsim-port-hamiltonian-8cx9i): passivity as a property of
the FORMULATION — `dx/dt = (J - R) grad H + G u`, `y = G^T grad H`
with `J` skew and `R` PSD — so any power-conjugate interconnection is
passive by construction, discrete gradients give energy-exact time
stepping, and Galerkin reduction preserves the structure.

## Public types and semantics

### Continuous harmonic response of the complete linear system

`harmonic::HarmonicModel::{from_ode,from_descriptor}` prepares the original
quadratic pHS once. `factor(omega)` factors its continuous AC operator once per
frequency, and the returned `HarmonicFactor::response` reuses that LU for any
number of source phasors. Peak amplitudes follow `exp(-i omega t)`; cycle means
therefore carry the appropriate one-half power and one-quarter energy factors.
No integration clock, frequency warping, peak fitting or artificial damping is
part of this path.

For a descriptor, the solved equation is
`[-i omega diag(I_diff,0) - (J-R) diag(Q_diff,I_alg)] x = G u`.
All algebraic multipliers, mutual storage, internal relaxation states and
original port directions are retained. Algebraic variables enter as efforts,
including source currents and transformer reactions, but store zero energy.
The singular descriptor mass is never inverted, and individual resonant
components are not eliminated through a singular impedance. A genuinely
singular complete operator still refuses. Nonlinear or undeclared storage
refuses; this is not an implicit tangent approximation.

The response exposes state, effort and port-output phasors, mean stored energy,
per-port supplied power, and independently evaluated `Re(e* R e)/2` loss.
Loss minus supply is reported without clipping or inferring heat from balance.
Independent row/column equilibration handles mixed units; the original equations
must pass a componentwise backward-residual ceiling, including small algebraic
constraints. This residual is not a forward-error, linewidth or pole certificate.

Dimension (hard cap 1024), internal numerical payload and dense work are admitted
before matrix allocation through `HarmonicBudget`; no state is truncated to fit.
The payload bound excludes source-owned storage, allocator overhead and process
RSS. Cancellation is polled during assembly, scaling, the existing complex LU,
response and diagnostics; source and prepared factors are immutable on refusal.
Repeat execution on the same host follows the existing deterministic LU owner.
No cross-ISA proof, transient-decay guarantee, nonlinear sounding threshold or
physical validity beyond the supplied coefficient laws is inferred.

G1/G3 regressions in `tests/harmonic_response.rs` compare mechanical amplitude,
phase and power with the damped-oscillator solution; Maxwell dispersion/loss
with its analytic complex stiffness; a high-Q RLC with its Kirchhoff impedance;
mutual-winding resonance splitting in both polarities; pure algebraic
transformer power; and coordinate scales `1e-120..1e120`. True lossless poles,
capacity limits, unsupported storage and cancelled reusable solves are tested.

### Finite resonance poles of the complete continuous model

`HarmonicModel::resonance_pole` refines one caller-seeded finite free mode of
the original `s diag(I_diff,0) - (J-R) diag(Q_diff,I_alg)` pencil through the
single residual-gated, cancellable bordered-Newton owner in fs-la. All mutual
storage, Maxwell memory states, algebraic efforts and original port directions
remain present. The singular descriptor matrix is never inverted. Fixed row
conditioning uses the declared spectral scale and retains the analytic
derivative and original state coordinates; an accepted result must also meet
the model's existing componentwise residual tolerance in the ORIGINAL split
pencil.

`HarmonicPole` returns the Laplace pole `s`, `omega=i*s` under the existing
`exp(-i omega t)` convention, original state and effort vectors, and original
`G^T effort` observations. Negative `Re(s)` describes decay; growth is never
clamped. The vector has unit numerical 2-norm, not unit physical energy.
Reported `H = x* Q x / 2` and `D = Re(e* R e)` share that amplitude, and
`2 Re(s) H + D` is evaluated independently from the physical storage and
resistance rather than set to zero. Algebraic efforts store no energy.

Separate `PoleOptions` admit dimension, dense numerical payload and total
iteration work before the lazy pencil callback allocates. Cancellation reaches
assembly, bordered LU, iteration and physical observations. Refusal leaves the
prepared model and existing reusable AC factors immutable. The budget excludes
the already prepared source model, allocator overhead and process RSS.

This is local isolated-pole refinement, not a complete spectrum search, contour
mode count, nonlinear limit cycle or global convergence certificate. Infinite,
repeated, defective or unresolved poles can refuse. A finite internal mode may
be unobservable at the selected ports; it is not automatically a transfer
function pole. Physical validity remains that of the supplied linear laws.

The four native regressions in `tests/harmonic_poles.rs` compare mechanical
frequency, damping, original coordinates and port output with the closed form;
a singular-descriptor RLC pole with its Kirchhoff constraints; a hereditary
mode with independently eliminated Maxwell dispersion, memory and loss; and
cancelled/budget-refused solves with deterministic retry and an unchanged
reusable harmonic response. The eight existing harmonic-response regressions
also pass with this owner.

### Stable constitutive energy changes

Quadratic storage evaluates its midpoint remainder analytically, retaining
represented-midpoint rounding and admitted matrix asymmetry without subtracting
large endpoint energies. Direct sums retain this capability only when both
children provide a complete remainder. Native, descriptor and prepared stepping
use the same constitutive endpoint difference in their energy reports; the
prepared path reuses its scratch buffers. This calculation remains independent
of source work, physical dissipation and the solved update equation. It does not
force a residual to zero or turn numerical error into heat. The regression
source checks a quartic system shifted by an arbitrary 1e30 J datum, composition,
large quadratic states, represented midpoint rounding and finite extreme Q.

A resistor whose terminals are the same node has zero incidence and contributes
neither dissipation nor a ground connection. Its conductance is not stamped, and
it cannot hide a floating node. The circuit regression compares the actual
source work and Joule heat with the independent V²/R result.

### Exact-Bessel circular wall memory

`ZkDiffusive` converts even truncations of the exact cylindrical Bessel
continued fractions into positive viscous and thermal Foster arms. The existing
Jacobi eigensolver owns the spectral decomposition. `series_drag_per_s` is
Poiseuille resistance divided by inviscid inertance; each `ZkDiffusiveTerm`
supplies `gain_per_s * s/(s + rate_per_s)`, normalized by the corresponding
inviscid inertance or compliance. No acoustic baseline storage is duplicated.
The low-frequency total inertance is 4/3 of the inviscid value and total
compliance is gamma times the adiabatic value at every pole count.

Pole counts 1..=64 are an explicit diffusion resolution axis. Construction
refuses nonpositive rates/weights or unresolved spectral decomposition, and
checks both complex excess wall operators and their real losses against the
existing exact Bessel owner at 129 logarithmic points in the supplied band.
The user supplies the error ceiling; no automatic clipping or fallback occurs.
The reported maximum errors are sampled evidence, not interval bounds. This
owner supplies immutable coefficients, not a time integrator. Focused regression
source in `zk/diffusive/tests.rs` covers a 2/4/8/16 pole ladder over viscous shear
1..30, static storage, physical gas scaling, and refusal boundaries.

### Opt-in linear evolution (IRP M5.2-W1)

`linear_exact::{split_quadratic, step_split}` requires explicit
`Storage::quadratic_part` data; gradient sampling cannot prove a storage linear.
The complete Jacobi energy basis is retained, including zero-energy coordinates.
`LinearSplitBudget.max_modes` counts scalar basis coordinates (two per ordinary
second-order oscillator), and the augmented exponential is capped at dimension
256. No modal truncation is performed. Nonfinite operators, stale splits,
unsupported storage, capacity exhaustion and unrepresentable work refuse.

`linear_exact::exponential_with_budget` exposes the same dense Taylor/scaling
algorithm to explicit cold preparation callers. `DenseExponentialBudget`
declares matrix-dimension, simultaneous dense-byte and multiply-term ceilings;
the opt-in dimension limit is 4096. Checked preflight reserves five square
matrices including the supplied input, plus 32 Taylor products and 48 squarings
before workspace allocation. Overflow, insufficient budgets, nonfinite input,
unconverged series and unrepresentable output refuse. The ordinary
`exponential`, split and ZOH paths retain their original 256-dimensional cap
and arithmetic. This entry permits complete-operator preparation; it supplies
no truncation, hot-path timing or physical-validity guarantee. Its regression
source checks larger explicit admission, unchanged default refusal and exact
byte/work boundaries. Native execution is pending.

`ExactZoh` advances the full coupled `(J-R)Q` operator with held input. Port
outputs are step averages. Physical heat independently integrates
`grad H^T R grad H` along that same complete trajectory. It retains a positive
square-root factor instead of expanding a Gramian whose terms can cancel near
a force-supported equilibrium. Initial loss effort and the original force
imbalance supply numerically stable coordinates; no equilibrium solve, state
reduction or change of physical generator occurs.

Positive eight-point Gauss quadrature on a norm-bounded short interval and QR
semigroup doubling preserve the heat factor. The existing fs-la Householder QR
owns factorization, and the existing bounded exponential owns every transition.
Absolute resistance scaling retains even very small physical damping; a negative
computed resistance eigenvalue refuses instead of being clipped. Exactly zero
resistance gives exactly zero heat. No inverse-growing passive propagator is
formed. Separate resistance channels still observe the same full original
trajectory. The original state/port admission and 48-doubling limit remain.

Source work and constitutive storage change are computed independently, and
the original 1e-12 relative energy-balance gate remains. Heat is neither an
endpoint discrete-gradient approximation nor an inferred energy remainder.
Generic signed quadratic observations retain their existing Gramian integration.
This is floating-point evidence, not an interval certificate or cooperative
cancellation guarantee; callers poll at their existing block boundary. Purely
quadratic storage is mandatory. `PreparedExactZoh` caches the same heat factor,
and its successful `step_into` reuses all numerical buffers.

Four physical regressions in `tests/exact_heat.rs` compare heat with independent
physical-equation RK4 and analytic limits: exact static silence and 1e-10 motion
with resistance down to 1e-200; high-Q and million-per-second relaxation;
relative-velocity damping nullspaces; and complete off-diagonal damping plus
Maxwell memory, clock subdivision and atomic refusal. Existing prepared-runtime
allocation, exact-flow, stiff-energy and multirate regressions pass in native
numerical testing. This does not claim measured wall-time acceleration or full
workspace qualification.

`FrequencyMatched` replaces Q by a per-mode matched quadratic energy while
retaining every term of the original nonlinear Hamiltonian. For lossless
energy-coordinate eigenfrequency omega its scale is
`tan(omega*dt/2)/(omega*dt/2)`. The ordinary Gonzalez discrete-gradient identity
then applies to the MODIFIED energy. Both the caller's angular-step bound and
`omega*dt < pi` are mandatory; no clamping at the tangent pole is allowed.
Static compliance changes, and damped decay rates are not exact in this mode.
Band admission does not waive finite-representation or Newton admission.
Newton refuses at its strict residual tolerance rather than accepting a best
iterate. The original `step` and prepared step remain unchanged.

G1/G0 tests in `linear_exact::tests` cover analytic held-force damped evolution
over 10,000 steps (1e-12 absolute at unit scale), lossless phase (0.01 cents),
exact damping radius (1e-12 relative), independent energy balance (1e-12
relative), a live quartic remainder, cap/cap+1, the tangent-pole refusal,
static compliance, free inertia, stale splits and unsupported storage.
`linear_exact::gramian::tests` additionally checks million-relaxation-time
decay, independently integrated stiff source work, nilpotent transport and
noncommuting damping/rotation. `tests/stiff_energy.rs` exercises the complete
public propagator across quadratic storage scales, separate physical heat
channels and driven thermal/electrical-style steady states. These retain the
independent 1e-12 energy-balance gate rather than widening its tolerance.

`admit_midpoint_band(rate, dt)` is the one-cent admission for every midpoint
(Gonzalez/Cayley) linear block. The caller supplies a bound on ALL its linear
pole magnitudes, damping and coupled stiffness included. The block is admitted
only if rate*dt <= `MIDPOINT_ONE_CENT_ANGULAR_STEP` = 0.04 rad, which keeps the
Cayley phase error below one cent even near critical damping. Outside the band
it refuses with `ResolvedBandError::Outside`, naming the largest admitted
frequency; invalid inputs refuse with `Invalid`. Coefficients are never retuned
or clamped: the remedy is a finer clock or an exact propagator. `step_render`
uses `ExactZoh` for exactly quadratic storage and otherwise steps midpoint
only inside this band. Tested by
`linear_exact::tests::midpoint_resolved_band_cap_cap_plus_one_and_damped_poles`
(cap, cap+1e-12, invalid inputs). Each render path's integrator and band is
inventoried and measured in fs-couple `tests/integrator_warping_audit.rs`.

- `Storage` trait (`hamiltonian`, analytic `gradient`);
  `QuadraticStorage` (`H = x^T Q x / 2`, `Q` admitted symmetric PSD),
  `SeparableStorage` (per-state scalar laws), `SumStorage` (direct
  sum).
- `PortHamiltonian::new` — admission: `J` skew, `R` symmetric PSD
  (Jacobi eigenvalue check), dimensions; refusal, never repair.
  `from_raw_parts` bypasses admission for trusted callers and mutation
  batteries — the supply-rate audit is what catches violated trust.
- `with_relaxation_branches` appends Maxwell internal strains to any admitted
  storage/port system. Original states and ports retain their order; each arm
  appends `z = sqrt(k) epsilon_viscous`, adds `(sqrt(k) b.x-z)^2/2` to H,
  and adds `1/tau` to R. Thus `z_dot=(sqrt(k) b.x-z)/tau`; arm dissipation
  is `(sqrt(k) b.x-z)^2/tau`, and the existing base H contains equilibrium
  stiffness. Coefficients must be finite, stiffness and relaxation rate positive,
  and the strain projection nonzero. The caller supplies loading history:
  zero z denotes unrelaxed internal strain; fully relaxed z is `sqrt(k) b.x`.
  G1 tests check two-arm exponential stress, analytical dissipated energy,
  second-order convergence, and the hereditary ramp integral through the actual
  stress port. The same step and interconnection machinery applies. No material
  identification, frequency-domain validity or thermal evolution is inferred.
  The composed storage keeps the original nonlinear storage's stable midpoint
  remainder and adds every Maxwell arm's exact quadratic remainder directly
  from coordinate increments, including represented-midpoint rounding. This
  preserves small restoring forces and material work beside a large stored
  energy or arbitrary energy datum; no extra state or hot-path allocation is
  introduced. An unavailable base remainder retains the original fallback.
  `fs-nlmodal/tests/relaxation.rs` checks identical nonlinear motion and work
  after adding a 1e24 J energy datum, nested material histories, and independent
  polarized stored-energy, force-work and Maxwell/Kelvin heat identities.

- `discrete_gradient` — the PINNED Gonzalez midpoint formula (order
  2): satisfies `dg . (b - a) = H(b) - H(a)` identically.
  `Storage::midpoint_energy_remainder` optionally supplies a cancellation-stable
  evaluation of its original correction numerator. It must account for EVERY
  stored-energy contribution and the supplied actual midpoint; `None` retains
  endpoint-energy subtraction. Both the canonical residual and analytic Newton
  use this same value. This changes numerical evaluation, never the Hamiltonian,
  loss law, solver request, rounding-floor policy or iteration cap. fs-dcontact's
  power-law implementation requires a current borrowed quadratic inner matrix
  (or a complete stable inner callback) and retains the total contact storage's
  nonlinear status.
  Native HC regression and independent double-double remainder tests are the
  qualification gates; this hook is not authority to replace an unknown energy
  contribution by zero or waive algebraic admission.
- `step` — implicit discrete-gradient step (approximate-Newton with
  central-difference Jacobian and state-relative target 1e-13). On
  stagnation the best iterate must pass the same quantitative algebraic
  admission as the final returned state; there is no 1e-6 fallback.
  Returns `StepRecord` with the
  ledger: `delta_h`, `dissipated = dt dg^T R dg`, `supplied = dt u^T
  y`, and `solver_residual` (the ACHIEVED residual, disclosed).
  `balance_residual()` restates the update equation for an admitted
  skew-J system (solver diagnostic, NOT a structure audit);
  `supply_defect() = delta_h - supplied` is the independent passivity
  audit, valid to the DISCLOSED band `~n * |dg|_inf *
  solver_residual` — audit thresholds must scale with it (the crate
  states the blind band rather than hiding it).
- `interconnect` — canonical skew pairing `u_a = -y_b`, `u_b = +y_a`
  over chosen port pairs; unpaired ports stay external (a's first);
  the composite is RE-ADMITTED (no false passivity). Associative on
  structure matrices for disjoint pairings.
- `reduce_galerkin` — `V^T J V` / `V^T R V` / `V^T G` with
  `H_r(x_r) = H(V x_r)`; skewness and PSD-ness survive Galerkin
  automatically and are re-admitted anyway. v1 requires quadratic
  storage (probed and refused otherwise); the reduced energy deficit
  at t = 0 equals the energy outside the basis exactly.
- Zoo: `mass_spring_damper`, `helmholtz_resonator` (pressure-driven
  lumped acoustic 1-DOF with unflanged neck correction),
  `helmholtz_resonator_flow` (dual port: injected volume velocity
  in, pressure out), `helmholtz_resonator_radiating` (damper =
  compact-mouth `Re Z_rad(ω₀)`), `lc_ladder`,
  `lc_ladder_terminated` (resistive far-end load — compact
  radiation on a discrete waveguide), `acoustic_cylinder`
  (inviscid `ρ,c` cylinder as that ladder; two inlets share the
  mouth so a blow and a transformer body can coexist; open end
  uses compact `Re Z_rad` and the same fit's mass `L = −X/ω`
  — unflanged `Δℓ = 0.6133 a`, flanged `0.8216 a` via
  `acoustic_chain_mouth` — at the quarter-wave pin; a flanged
  mouth above the fit ceiling is the Rayleigh piston;
  closed stays the lossless last inertance),
  `acoustic_waveguide` / `AcousticTap` (`open_fraction` σ:
  compact open neck inertance, or a long neck
  (`kℓ > 0.2` or `ℓ > 4b`, `ℓ` = chimney + Dalmont inner)
  as the same 2-cell LC line plus flanged mouth the TMM
  chimney is; closed cavity `C`, or the mix
  `Y = σ Y_open + (1−σ) Y_closed`; a `ViscothermalPin`
  puts series `R` on a lumped neck and thermal `G` on
  remaining `C`, or the bore's per-cell wall law on a
  line; `foster_branches > 0` adds the bore's Foster
  series on each lumped neck, series plus thermal Foster
  on each chimney cell of a line, and the thermal Foster
  shunt on each remaining pad `C`),
  `acoustic_chain` / `AcousticSection` (concatenated LC runs
  with an area jump at each interface — a muffler, a
  constriction, a cone sliced into cylinders),
  `ChainOutlet` / `acoustic_chain_outlet_wall` (the far end:
  `PressureRelease` is `p = 0` behind the last flux, the
  quarter-wave stop family that `mouth: None` has always
  built; `Rigid` cuts the last flux, so `U = 0` and a
  flow-driven inlet sees the half-wave family; `Radiating`
  is the compact mouth load. G1
  `rigid_outlet_converges_to_the_closed_closed_resonance`,
  frankensim-gs17m),
  `WallPin` / `acoustic_chain_mouth_wall` (locally reacting
  wall: per-cell shunt `L = σ/A_w`, `C = A_w/K`,
  `R = r/A_w` on `A_w = 2π a dx` — `×√(1+s²)` on a
  taper; TMM `input_impedance_wall` adds the same
  `Y' = 2π a slant / Z'` to the gas shunt; `None` is
  rigid; not a plate and not mean flow),
  `ViscothermalPin` (exact Bessel ZK pin, with series R/inertance and
  thermal G/compliance evaluated at the same physical frequency;
  zero `μ` is the lossless mutation; `foster_branches`
  collocates Foster networks to Bessel Zwikker–Kosten
  `F(r_v)` at every shear number, and on an open mouth
  to `Re Z_rad(ω)` of the compact / piston load — the
  pin `X` still supplies the radiation mass; tabulated
  `R(ω)` remains the FIR path),
  `zwikker_kosten_f` / `zk_one_minus_f` (also the fs-duct Bessel TMM),
  `foster_sqrt_omega_terms` (positive by construction: the log-trapezoid
  rule of the exact Stieltjes form of the wall law `√ω(1 - i)`, tails lumped
  into the end poles; no collocation and no fallback split),
  `foster_match_re`, `slice_linear_taper`,
  `spherical_cone` (a widening cone as a 1-D wave on `ψ = x p`.
  Physical ports are `p = ψ/x` and `U = x U_ψ + (α/ρ)∫ψ`. The
  near-field term enters only where ψ meets physical flow, so the
  inlet port carries the one shunt `L = ρ x₁/α`; an interior tap
  joins at `1/x` with no near-field term. The open mouth is the same
  compact `(R, X)` on the last ψ-flux with `Z_ψ = x_out² Z_phys`.
  Narrowing tapers are refused. A lone widening taper in
  `acoustic_chain` is this object. Any other taper, alone or in a
  multi-section chain, is frustum LC cells with the exact
  `ρ dx/(π r_a r_b)` inertance and frustum volume.
  G1 `tapered_chains_converge_to_the_spherical_wave_resonance` checks
  the lossless p = 0 resonance against `k L + atan(k x1) = π` and
  the cylinder+cone root. The former per-cell shunt diverged,
  frankensim-d1j6s), `modal_bank`
  (mass-normalized modes -> canonical pHS: the bridge from the
  eig/plate beads), `duffing_oscillator` (non-quadratic exercise).
- `compact_radiation_impedance` — low-`ka` Levine–Schwinger /
  flanged piston `(R, X)` under `e^{-iωt}` (mass-like `X < 0`),
  refused above `ka = 1`. Same numbers as `fs_duct::Termination`.
- `baffled_piston_impedance` — the sole uniform circular-piston owner,
  exact Bessel–Struve closed form in an infinite rigid baffle, acoustic
  `p/U = (rho c/S) [1 - 2J1(2ka)/(2ka) - i 2H1(2ka)/(2ka)]`.
  Area converts to mean p/v; area squared converts to force/velocity.
  No rings parameter or numerical surface discretization. Nonpositive,
  nonfinite or unrepresentable derived scales refuse with NotPsd.
  G2 `piston_closed_form_certified` targets 1e-13 relative at independent
  analytic landmarks, resistance's 1e-12 low-ka band and the high-ka limit.
  At ka=1e-3 the cubic reactance correction is 2.26e-10: the leading-term
  1e-12 absolute band in .101 is inconsistent with the exact solution;
  the test retains the cubic term and explicitly detects this difference.
  Native qualification is batch-pending. No unbaffled-pipe claim.
- `side_hole_inner_length` / `side_hole_neck_length` — isolated
  disc `8/(3π) b`, or Dalmont–Nederveen–Joly `d = b/a` on a
  finite bore, plus wall-flanged `0.8216 b`. Same `t_eff` as
  the fs-duct open-hole shunt. A compact open `AcousticTap`
  uses this length and flanged `Re Z_rad`. A long neck uses
  chimney plus inner on the cells and `0.8216 b` on the last
  flux.
- `side_hole_series_length` — Nederveen `t_s = −0.37 b²/a`. The
  TMM open hole is the T-junction `series(Z_s/2)·shunt·series(Z_s/2)`;
  the ODE tap subtracts that inertance from the station cell.
- `side_hole_mutual_length` — extra series on hole `i` from
  neighbor `j`: `t_m = −0.37 (b_i b_j / a) e^{−s/a}`. Added
  into the same T-junction `t_s` (TMM and ODE).
- `series_impedance_ports` — ODE series of two 1-port impedance
  systems (same `u`, `y = y_a + y_b`).
- `common_effort_capacitor` — two-port `C` with shared pressure
  `p = q/C` and `q̇ = U₁ + U₂`. That is the ODE image of a
  Kirchhoff effort junction.
- `kirchhoff_parallel_step` — index-1 Newton split of an external
  flow across two 1-port impedance systems so their pressures
  match. This is Kirchhoff current law as a time-step.
- `common_effort_dirac` / `common_flow_dirac` /
  `join_port` (1-junction on a named port pair; leftover
  ports stay external — a bow, a blow, a side load),
  `common_effort_star` / `DescriptorPortHamiltonian` /
  `step_descriptor` — the true composite Dirac structure of a
  Kirchhoff star. Two members are a 0-junction
  (`U_a = λ`, `U_b = U_ext − λ`, `y_a = y_b`). `N` members
  use `N−1` multipliers. On admittance ports the same `J` is the
  mechanical 1-junction (common velocity, forces split). Gonzalez
  on the differential block; algebraic rows are the common-output
  residual. The capacitor and the Newton split are regularizations.
  `DescriptorPortHamiltonian::new` admits explicit finite J/R/G matrices and
  a differential prefix, enforcing zero algebraic storage by passing only that
  prefix to the supplied energy law. Read-only `structure()` and `storage()`
  let coupled material owners reuse original circuit stamps and declared Q
  rather than build another circuit solver. Differential quadratic metadata is
  lifted with an exactly zero algebraic block. The descriptor Gonzalez norm,
  cancellation guard and corresponding rounding estimate exclude algebraic
  coordinates; actual multiplier efforts remain their original midpoints.
  Nonfinite state, storage, effort, flow or ledger refuses as `DescriptorDomain`;
  a NaN residual cannot count as convergence. Centered Jacobian probes and Newton
  updates retain their original size when finite, with at most 24 halvings to
  regain the declared storage support. No state clamping, material extrapolation
  or convergence-budget relaxation occurs.
  `tests/descriptor_storage.rs` exercises a nonlinear RC energy defect, changed
  algebraic coordinate scales, native RLC reconstruction, source-law admission,
  a 100-micrometre charged plate with a 100-volt algebraic source, explicit
  nonfinite refusals and a supported nonlinear root requiring update backtracking.
- `transformer(a, b, port_a, port_b, n)` — power-conserving
  `u_a = n y_b`, `u_b = −n y_a`. A plate area, a hydraulic ram,
  and a lever are this object.
- `moving_end_waveguide` — free-fixed taut eigenfunctions with a
  1-port at the free end (`φ(0) ≠ 0`). Fixed-fixed sines cannot
  Dirac-join a body. `modal_bank_ports` is the same bank with
  one drive column per port.
- Three-pHS string–plate–cavity is
  `common_flow_dirac(moving_end_waveguide, transformer(plate, cavity))`.
- Memoryless dissipative ports: `bernoulli_volume_flow` is the
  two-sided jet `U = w h sgn(Δp) √(2|Δp|/ρ)` (dissipative on both
  branches); `quasistatic_aperture_opening` is the zero-mass
  reduction of a linearly restoring slit; `regularized_coulomb`
  is `F = −μ N tanh(v/v_reg)` (dissipative, stick-slip by
  construction). A beating reed, a vocal fold, a relief valve, a
  leaflet, a bow, and a brake are fillings of those ports plus a
  1-DOF `mass_spring_damper` or a `modal_bank`. Music is not a
  special case.

## Invariants

The cylindrical ZK owner uses ζ=rv√(-i) in +iω telegraph bookkeeping:
Im F<0 and Re(iωL/(1-F))>0. Double-double defining series evaluate both
F and its cancellation-free deficit below |ζ|=12; the equivalent I1/I0
Lentz fraction is used above it. The admitted range is [1e-6,1e6]; a
16392-iteration cap and the delta stop bound work and refuse exhaustion.
The old 1e-12 stop arithmetic is preserved above the switch to avoid
gratuitous WideTube-chimney golden changes. G0/G2 `zk_f_certified_all_regimes`
targets 1e-12 agreement against a dd backward fraction on 400 log points,
1e-13 switch continuity and 1e-6 two-term wide-limit agreement. The
48th retained series term is checked below 1e-17 of J0 at the switch.
`zk_one_minus_f_small_rv_no_cancellation` targets 1e-14; domain/cap tests
also certify the extended endpoints. `zk_pin_keeps_positive_resistance_and_thermal_compliance`
checks both series/shunt quadratures and refusals. These new tests await
central native qualification. PhsError has no domain variant, so named
NotPsd refusals are retained for this item; the enum label is historical.

1. Admission rejects non-skew `J`, asymmetric or non-PSD `R`, non-PSD
   `Q`, and bad dimensions by name.
2. Lossless conservation: undriven `R = 0` systems hold `H` to
   1e-10 relative over thousands of steps (state-relative Newton
   tolerance; an absolute tolerance leaks — executed lesson).
3. Damped ledger exactness: per-step balance residual is solver-zero
   and the summed dissipation ledger equals the total `H` drop.
4. Gonzalez order 2 on non-quadratic `H` (Richardson-verified) with
   `H` STILL conserved exactly despite trajectory truncation error.
5. Interconnection reproduces the hand-assembled monolithic system's
   trajectories to 1e-10 and is associative on structure.
6. Supply audit fires on smuggled mutations: symmetrized `J` and
   sign-flipped `R` both violate `supply_defect <= 0` observably; a
   one-sided (broken) coupling map is refused at admission.
7. Reduction preserves structure (re-admits) and the reduced modal
   bank's impedance certifies passive under the INDEPENDENT fs-vfit
   Hamiltonian test.

## Error model

DEF-05 algebraic admission is shared by `step`, `step_descriptor`, prepared
FD/analytic/dense/condensed stepping and the Kirchhoff split. `SolveCertificate`
records the site, residual, scale, iterations, stagnant iterations,
`requested_budget`, `rounding_floor`, `effective_budget` and `floor_raised`.
The strict requested target remains 1e-13 times the existing step scale;
Kirchhoff uses its existing 1e-8 target, including its final assessment.
Residuals at the cap accept; next-up residuals refuse as `NewtonRefused`.
An explicit `SolveBudget` below the evaluated floor refuses as
`BudgetBelowRoundingFloor`, retaining both numbers. Prepared declarations use
`set_solve_budget`; `set_iteration_limit` retains its update-cap semantics.
The default APIs explicitly record an effective ceiling
`max(requested_budget, rounding_floor)` and the corresponding `floor_raised`.
This resolves the bead's contradictory below-floor/default-noise clauses under
MagentaBadger's 2026-10-16 decision; it is never an implicit relaxed tolerance.
Dropping the floor from a raised certificate makes it internally inconsistent
and refuses. Returned Gonzalez records carry `Some(certificate)`; the separate
linear propagators currently carry `None`, with no Gonzalez-admission claim.

The arithmetic model uses `gamma_k = k*EPSILON/(1-k*EPSILON)`, with
`k = 4*n+8` for each dense residual row. It bounds the sum of the evaluated
absolute operands: `|x1|+|x0|+|dt|*(sum((|J|+|R|)*|dg|)+|Gu|+|D|)` for a
differential row; algebraic rows omit the state terms and use factor one.
When the Gonzalez quotient is active under its existing cancellation guard,
`gamma_(4*n+8)*(|H0|+|H1|+sum(|dg_i*dx_i|))` is additionally propagated through
`|dx_i|/sum(dx^2)` and the absolute flow operator. These operands and operation
counts compute the pointwise floor; neither an observed residual nor a fitted
noise constant defines it. At a vanishing increment this pointwise model alone
is insufficient: the best iterate can be `x1=x0`, with no quotient term, while
every nearby Newton trial amplifies energy-subtraction roundoff.

For storage exposing a borrowed exact quadratic matrix through
`Storage::exact_quadratic_matrix`, and no nonlinear dissipation callback, a local
attainable-residual envelope also applies. Write `E` for the energy roundoff
above, `B_i = factor_i*sum_j(|J_ij|+|R_ij|)`, and
`L_i = 1_(i differential) + B_i*||Q||_inf/2`. Descriptor multiplier efforts
replace the corresponding `Q` rows with identity rows, so their effort bound
uses `max(||Q||_inf,1)/2`. The ideal quadratic residual is affine with row
Lipschitz bound `L_i`. Moving a distance at most `rho` from its root changes that
row by at most `L_i*rho`, while a Gonzalez correction evaluated at increment
length `rho` contributes at most `E*B_i/rho`. Minimizing their sum gives
`rho_i^2=E*B_i/L_i` and `2*sqrt(L_i*E*B_i)`. Only when the actual increment
satisfies `||x1-x0||_2^2 <= rho_i^2` is this local envelope added to that row's
pointwise floor. Thus the `x1=x0` endpoint has a finite attainable allowance;
ordinary steps and nonquadratic laws retain the pointwise model. The borrowed
matrix hook preserves prepared stepping's no-allocation guarantee. A nonlinear
dissipation callback disables this extension because its slope is not supplied
by `Q`. Explicit requests below either computed floor still refuse.

Kirchhoff uses `gamma_3*(|y_a|+|y_b|)` for its subtraction.
This conditional floating-point envelope assumes conforming, read-only storage
and dissipation callbacks. It does not bound callback evaluation error,
constitutive error, finite-difference truncation, conditioning of a trajectory,
or discretization error, and it mints no physical validation authority.

G0/G4 `tests/newton_refusal.rs` covers all four sites, a consistent staircase
Hamiltonian, cap/next-up, explicit below-floor, exact update caps, the real reed
lamella's full stationary tail and visible default floor, a dropped-floor mutant,
and unchanged caller buffers after refusal. The native regression and music-baseline lanes remain
required before any closure or golden change.

Typed `PhsError` (symmetry class, PSD, dimensions, Newton stall, port
pairing); no silent degradation.

## Determinism class

Deterministic: fixed iteration caps, deterministic Jacobi/LU kernels
from fs-la, no RNG or time dependence.

## Cancellation behavior

Pure synchronous computation; per-step cost is one small dense Newton
solve. No `Cx` integration (workspace `frankensim-ccmn` effort).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

- Circuit DAE battery (`src/circuit.rs`, cases kd-001..kd-005):
  series RLC vs the analytic impedance at three frequencies straddling
  resonance with the supply-defect audit below 1e-9 throughout; the
  pure LC loop (the index-relevant case) holding energy to ~1e-16
  over 1e5 discrete-gradient steps; the ideal-transformer two-loop
  fixture matching the reflected-load current analytically with the
  voltage law pinned; refusals by name; and the audit non-vacuity
  falsifier (a sign-flipped supply fails by >1e3x).

- Device battery (`src/device.rs`, cases dv-000..dv-005): the
  descriptor READ-LAW regression (algebraic efforts are midpoints of
  the stored coordinate — the consistent-IC raw read is 2x and a
  direct potential read is a bug by construction); single-triode
  small-signal gain matching the plate-curve-derived
  `mu*RL/(RL+rp)` from the card's own finite differences (60.9 vs
  60.9); the EMERGENT harmonic ladder (H2 rising 0.26%->3.2% across
  drives 0.05->0.6 with H2 >> H3 throughout — asymmetry from the
  card + bias, nothing shaped); the diode clipper vs an in-test
  bisection oracle at residuals <=1e-19; refusals by name (validity
  region, bias admission); bitwise determinism across re-runs.

## Conformance tests

`tests/phs.rs` (9): admission rejections; Gonzalez identity +
lossless conservation; damped ledger exactness; order-2 Richardson on
Duffing; interconnection vs monolithic; structural associativity;
driven power balance + supply audit; mutation battery (symmetrized J,
sign-flipped R via `from_raw_parts`, one-sided coupling refused at
admission); reduction (structure re-admission, t = 0 energy
exactness, fs-vfit passivity certification of the reduced impedance,
realized H-error inside an authored 2% envelope with the measurement
recorded).

Inline `valve_ports` (2): two-sided Bernoulli jet is odd, dissipative,
and scales as `√Δp`; the quasistatic aperture closes at the named
pressure.

`tests/reed_casebook.rs` (1): single-reed exciter as pHS components —
reed lamella (msd pHS) + modal bore bank + TWO-SIDED Bernoulli valve
(backflow reverses with dp, so `dp * U >= 0` is a property of the
law, not a guard); per-step component supply audits; an integral
mouth-work bound from state-recomputed energies; the per-step closed
accounting is labeled as the consistency restatement it is (implied
by component passivity + valve dissipativity); quasi-static
equilibrium pinned against the hand-derived analytic solution
(`q* = Sr pm / k`, zero bore DC pressure, `U* = w h* sqrt(2 pm /
rho)`) — the casebook's independent evidential content. JSON-lines
summary. The damped-ledger battery test carries its own independent
oracle: trapezoidal quadrature of `c v^2` from recorded states
matches the dissipation ledger within 2%.

## No-claim boundaries

- THE DAE DEFERRAL IS CONSUMED (music bead 3ez8g.9.1; the
  electric-guitar track was the named first consumer): the `circuit`
  module assembles R/L/C/source/ideal-transformer graphs into
  `DescriptorPortHamiltonian` systems where KCL/KVL are algebraic
  Dirac constraints. Each constraint enters the composite `J` once
  and ANTISYMMETRY manufactures its dual (the inductor's KVL row
  transposes into its KCL current; the capacitor branch row into both
  charge dynamics and node currents; the transformer's single voltage
  row into the power-exact current pair) — passivity by construction,
  never bookkeeping. Admission is loud (floating node, shorted or
  parallel ideal sources, non-physical elements refuse BY NAME);
  consistent ICs start with the dt=0 descriptor step. If redundant ideal
  constraints leave that Newton matrix singular, the existing fs-la SVD owner
  solves only algebraic efforts with the original charge and flux held exactly.
  The represented midpoint must satisfy the original unscaled rows under the
  same native solve budget; incompatible clamps still refuse. The minimum-norm
  multiplier choice does not claim uniquely determined initial branch currents.
  Source-free circuits expose an exactly empty n-by-zero G, so their native
  stamps reconstruct through the same descriptor constructor.
  `tests/circuit_initial_constraints.rs` checks compatible and conflicting
  voltage clamps, a subsequent physical RLC step, and source-free reconstruction.
  `assemble_circuit_with_mutual_inductance` accepts signed reciprocal links
  between original inductor branch indices. It admits the complete positive-
  definite L after scaling by sqrt(L_ii), prepares its inverse once, and stores
  H_m = phi^T L^-1 phi / 2 in the existing flux coordinates. The represented
  squared Cholesky pivots must exceed 16*n*EPSILON in that dimensionless matrix;
  this refuses unresolved leakage instead of inventing it through roundoff.
  The represented
  dimensionless inverse must satisfy max|C C^-1 - I| <= 1e-10; nonfinite or
  unresolved inverses refuse. Each pair is unique and names distinct inductors;
  pairwise coupling limits do not replace whole-matrix energy admission.
  `CircuitDae::inductor_currents` evaluates the actual storage gradient,
  including mutual terms. No source, capacitor, Kirchhoff stamp, integrator or
  energy ledger is duplicated. Empty links preserve the uncoupled path exactly.
  G0/G1 tests in `tests/circuit_mutual_inductance.rs` check interleaved branch
  mapping, independent coupled-RL KVL/heat, winding polarity equivariance,
  unequal winding scales, and the finite magnetizing/leakage transient of a
  loaded transformer against its closed form with second-order refinement.
  Negative-energy multiwinding data and perfect coupling refuse; exact ideal
  transformers retain the algebraic TransformerLink representation. Authored
  or measured L/M values carry no inferred core saturation, hysteresis,
  eddy-current loss, moving geometry or temperature law.
  The ledger stays exact
  through constraints, and `step_audited` returns the SUPPLY-DEFECT
  residual |dH + dissipated - supplied| — the independent check
  (balance alone is tautological under sign mutations; the audit's
  non-vacuity is itself tested by a sign-flip falsifier). Remaining
  DAE no-claims: no event switching in v1 (diodes/triodes arrive as
  smooth device laws in the device bead; switching topologies are a
  later, named extension); index analysis is structural (admission
  catches the classic pathologies; exotic graphs that stall Newton
  refuse loudly).
- DEVICE CARDS + ISLANDS (music bead 3ez8g.9.3): the `device` module
  carries Koren-class triode and Shockley-class diode-pair CARDS
  (typed parameters, named provenance, validity regions in the matdb
  shape) and couples them to `circuit` DAEs as nonlinear ISLANDS —
  the split-circuit image: the linear network steps through the
  descriptor machinery while the device iterates against it
  (Gauss-Seidel sweeps with a load-line-damped island update; fixed
  caps; analytic Jacobians only — the laws are smooth). Bias
  admission is a DC Newton refusal, not silent garbage. Distortion is
  CONSTITUTIVE: the harmonic ladder emerges from curvature at the
  operating point. No-claims: the split image is the BUDGET solver —
  the full-DAE monolith is the authority image and the bake-off
  belongs to the electric-gates bead; parameter provenance is
  authored-typical (Koren 12AX7 fit constants, silicon-class diode),
  not a measured-tube ingest; no matdb migration yet (schema is
  migration-shaped by construction).
- ACOUSTIC COMPONENTS ARE QUIESCENT-MEDIUM (decision record
  `docs/DECISION_MUSIC_MEAN_FLOW.md`, music bead 3ez8g.10.5): the lumped
  Helmholtz resonator, the `AcousticTap`s, and the `acoustic_chain`
  inviscid ladders model a medium at rest. Validity boundary: peak
  local acoustic Mach < 0.01 AND zero through-flow; forte-regime
  impedance shifts, vortex-shedding end corrections near tone holes,
  and jet-drive interaction regions are OUTSIDE every claim these
  surfaces mint. Escalation (a named convected uniform-flow card, then
  shear-flow/CAA territory) follows the ladder in that record and fires
  only on registered-corpus measured deviation — never as a silent
  rescale of an existing quiescent result.
- Non-quadratic storage through `reduce_galerkin` is deferred with the
  same trigger; the refusal is typed.
- A-priori trajectory/H-error bounds for reduction are a no-claim; the
  certified statement is the t = 0 energy deficit plus the measured
  realized error under an authored envelope.
- Oscillation-regime reed validation (threshold pressures, limit
  cycles) is not a pHS claim. The casebook validates the
  quasi-static regime where an exact analytic answer exists. A
  time-domain aperture ⊕ traveling-wave composition lives in
  fs-couple as a consumer of these ports, not as an instrument
  algebra.
- The jet and aperture maps are nameless fluid/solid ports. They
  do not mint a clarinet, voice, or valve product.
- fs-time strategy wiring is a follow-up owned at L3 (fs-time depends
  down, not up).

### Frozen storage during algebraic initialization

`CircuitDae::consistent_initial_state` admits a converged zero-duration Newton
result only when every differential coordinate retains its original bits. Any
roundoff drift retries the existing rank-aware algebraic solve at the original
storage, with the original unscaled constraint rows and solve budget. No charge,
flux, or other differential coordinate is projected or newly created to make
an ideal source fit. The source-onset regression covers an interleaved RLC
network with zero stored coordinates and nonzero instantaneous source/capacitor
currents; incompatible ideal clamps still refuse.

### Prepared complete flow and physical multirate coupling

`linear_exact::PreparedExactZoh` wraps the existing `PreparedQuadraticFlow` for
a fixed original quadratic owner and clock. It caches the complete exponential,
original resistance heat integral and all conjugate port integrals; no state,
mode, damping term or port is removed. Warm calls perform matrix-vector products
and quadratic contractions, with allocating result vectors. Heat is independently
integrated from the original loss law, never recovered as an energy remainder.

`multirate::MultiratePair` advances both complete physical states on independently
declared clocks of 1 through 16 native steps per common window. Its held interface
inputs satisfy `u_left = -average(y_right)` and `u_right = average(y_left)` using
the original integrators' actual window-average conjugate outputs. Consequently
interface works cancel over that SAME physical duration. Authored external inputs,
each owner's physical heat, and interface exchange are reported separately. The
whole-window constitutive energy increment is also checked independently against
the compensated sum of native increments; a supplied nonadditive remainder cannot
hide merely because each individual update closes its own equation.

Quadratic pairs prepare the complete affine interface response and its LU once.
Nonlinear pairs retain the original Gonzalez `StepWorkspace`, with bounded
interface Newton updates and at most 32 paired ports. Caller-declared local pole
rate and optional contact-resolution bounds are enforced before preparation.
All subsystem states are returned together only on success; refusal/cancellation
does not publish a state or advance a clock. Numerical scratch may change.

Relative energy gates use a conservative coordinate-dependent numerical scale
from `sum(abs(x_i*gradient_i(H)))`, actual changes, input work and heat. This scale
is invariant under a constant Hamiltonian datum; it is not itself thermodynamic
stored energy. Explicit absolute allowances are never silently enlarged.

This is a held-interface time discretization with second-order convergence on
smooth converged owners, not the exact interconnected flow. Internal subcycling
cannot eliminate common-window phase error: both clocks AND the coupling window
require refinement. Local exact flow does not certify a macro frequency band,
and no observation-rate anti-aliasing, nonlinear tangent-bound inference,
million-step drift certification or full multirate-bead closure is claimed. The
nine focused integration regressions exercise 10,000 lossless windows across four
clock ratios, physical heat/source separation, three-rung full-flow convergence,
nonlinear subcycle convergence, multiple ports and numerical scaling, explicit
clock limits, cancellation retry, stiff analytic heat/work and datum-invariant
rejection of an inconsistent storage law.

## Owned, reusable exact quadratic evolution

`PreparedExactZoh::new_owned` moves the original immutable physical system and
its fixed-clock cache together; `system()` exposes only read access. The
existing borrowed constructor remains. `step_into` reuses complete augmented
coordinates, endpoint, port-integral and energy scratch. No matrix exponential
or owner allocation occurs on a successful warm call with allocation-free
storage callbacks; both caller outputs remain unchanged on every refusal.
The allocating API shares these contractions and the original independent
heat, port-work and storage-change admission. No state or off-diagonal
coupling is dropped. `LinearSplit::is_exact_quadratic` exposes the original
storage classification for consumer dispatch.

`prepared_exact_owned` checks 1,000 warm calls with zero allocations, complete
coupled damping and Maxwell states, bitwise borrowed/owned work and motion,
direct exact propagation, energy-datum invariance, capacity and atomic refusal.
All nine existing multirate tests pass against the changed owner. These are
native numerical and allocation checks, not measured wall-time acceleration or
full workspace qualification.

### Resonance spectrum of the original descriptor model

`HarmonicModel::resonance_spectrum` discovers finite free modes inside an
explicit `ResonanceRegion` ellipse in the Laplace s-plane [1/s], using
`BeynOptions`. It shares the exact same lazy `sE-(J-R)W` pencil, fixed row
conditioning, bordered-Newton refinement and physical observation owner as
`resonance_pole`. This complete linear pencil is entire even when E is
singular. All algebraic constraints, mutual storage and material memory remain;
finite modes with zero port output are retained. E is never inverted.

Every returned `HarmonicPole` passes the unchanged original split-pencil
componentwise residual gate and reports original coordinates/efforts/ports,
physical energy and dissipation, and independent decay-power balance. The
spectrum reports numerical argument count, resolved rank, unrefined values,
nested-quadrature indicators and contour margin. Physical result storage and
observation work are reserved in addition to the shared numerical budget.
Cancellation or refusal leaves the source model and reusable AC factors intact.

Four native integrations in `tests/harmonic_spectrum.rs` recover an unobservable
mechanical resonance, all three poles of a Maxwell hereditary oscillator,
only the finite poles of an RLC descriptor with original Kirchhoff constraints,
and exact replay after cancellation with unchanged reusable AC output. The
four isolated-pole and eight harmonic-response integrations also pass against
the complete production fs-phs source. These finite-region results provide
numerical convergence and original-equation evidence; they do not certify an
entire spectrum, nonlinear sounding threshold, observable transfer pole, or
resolution of repeated/defective modes.

## Explicit physical load ports beside independent drive ports

`PortHamiltonian::with_parallel_load_columns` accepts the complete finite host-state-by-load-port matrix B. The load receives B^T grad(H_host), and the host receives the reciprocal reaction -B y_load through the existing `interconnect` owner. Every original state, constitutive law and external drive port remains in order. This permits a surface-volume-flow or distributed-velocity port to differ from a point-force drive without changing internal work. Signed, repeated and nonorthogonal columns are retained. The existing indexed `with_parallel_load` delegates to the same owner.

Callers supply physical units and explicit work/state capacities; this is no new integrator, reduction or impedance fit. Three focused native G0/G1 tests check signed reciprocal work, independent force work and total storage/loss, undriven-coordinate reaction, malformed projections, and exact structure equivalence of the indexed API. Full-workspace qualification is not asserted.

## Material histories with zero mechanical projection

`with_relaxation_branches` admits an exactly zero finite projection while still
requiring strictly positive finite branch stiffness and relaxation time. Such a
branch is an independent internal energy coordinate: its initial history stores
energy and decays into its own physical heat through the original relaxation
law. It adds no mechanical stiffness or work. Complete tensor material sources
can contain these directions; a zero observed force does not authorize dropping
their storage or histories. Dimensions, finite coefficients, relaxation-rate
representability and existing time/energy admission remain unchanged.

The authored exact-flow regression checks exponential decay of a nonzero
uncoupled history, unchanged mechanical coordinates, independently integrated
heat and the original energy balance. Native execution remains unverified for
this increment.

## Immutable constitutive sharing for transactional preparation

`PortHamiltonian::share_storage` copies the original J/R/G matrices during cold
preparation and shares its immutable constitutive callbacks. Evolving physical
history remains in the caller's explicit state vector. Hamiltonian, gradient,
declared quadratic data and stable nonlinear energy remainder are delegated
unchanged; this neither approximates the law nor infers a quadratic model.
`storage` exposes the original read-only law. The optional `Storage::shared_clone`
hook lets the internal sharing wrapper duplicate its handle without nesting
another callback wrapper on every cache preparation. Other storage owners keep
the default and require no implementation change. Read-only trial semantics and
accepted-state synchronization remain the caller's responsibility.

The focused `prepared::port_load` regression makes 1,024 cold copies of an
original oscillator with Maxwell memory and checks stable wrapper identity,
complete J/R/G, energy, effort, quadratic data, nonlinear-remainder delegation
and actual step/heat/work equivalence. Existing physical load-port regressions
remain unchanged. This is a composition and numerical-cache reuse seam, not a
new state integrator, threading guarantee or measured performance certificate.
