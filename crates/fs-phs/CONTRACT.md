# fs-phs — Contract

## Purpose and layer

Layer L2. Port-Hamiltonian systems (bead
frankensim-fsim-port-hamiltonian-8cx9i): passivity as a property of
the FORMULATION — `dx/dt = (J - R) grad H + G u`, `y = G^T grad H`
with `J` skew and `R` PSD — so any power-conjugate interconnection is
passive by construction, discrete gradients give energy-exact time
stepping, and Galerkin reduction preserves the structure.

## Public types and semantics

- `Storage` trait (`hamiltonian`, analytic `gradient`);
  `QuadraticStorage` (`H = x^T Q x / 2`, `Q` admitted symmetric PSD),
  `SeparableStorage` (per-state scalar laws), `SumStorage` (direct
  sum).
- `PortHamiltonian::new` — admission: `J` skew, `R` symmetric PSD
  (Jacobi eigenvalue check), dimensions; refusal, never repair.
  `from_raw_parts` bypasses admission for trusted callers and mutation
  batteries — the supply-rate audit is what catches violated trust.
- `discrete_gradient` — the PINNED Gonzalez midpoint formula (order
  2): satisfies `dg . (b - a) = H(b) - H(a)` identically.
- `step` — implicit discrete-gradient step (approximate-Newton with
  central-difference Jacobian, state-relative tolerance; on
  stagnation the best iterate is accepted iff its residual sits six
  orders below both the iterate scale and the initial residual, else
  a typed `NewtonStalled` refusal). Returns `StepRecord` with the
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
  uses compact `Re Z_rad` at the quarter-wave pin),
  `acoustic_waveguide` / `AcousticTap` (open side-branch neck
  inertance shunted to atmosphere — a tone hole),
  `acoustic_chain` / `AcousticSection` (concatenated LC runs
  with an area jump at each interface — a muffler, a
  constriction, a cone sliced into cylinders),
  `ViscothermalPin` (all-regime pin: wide-tube ZK for
  `r_v ≥ 10`, Poiseuille + isothermal-tending shunt below;
  zero `μ` is the lossless mutation; `foster_branches`
  collocates Foster networks to Bessel Zwikker–Kosten
  `F(r_v)` at every shear number),
  `zwikker_kosten_f` (also the fs-duct `LossModel::Bessel` TMM),
  `foster_sqrt_omega_terms`,
  `foster_match_re`, `slice_linear_taper`,
  `spherical_cone` (1-D wave on `ψ = x p`, physical ports
  `p = ψ/x` and `U = x U_ψ + (α/ρ)∫ψ` via a per-cell
  shunt `L = ρ|x|/α`; a single taper in `acoustic_chain`
  is this object; mixed cylinder+taper runs stitch it to
  the LC ladder with transformer `x` at each interface), `modal_bank`
  (mass-normalized modes -> canonical pHS: the bridge from the
  eig/plate beads), `duffing_oscillator` (non-quadratic exercise).
- `compact_radiation_impedance` — low-`ka` Levine–Schwinger /
  flanged piston `(R, X)` under `e^{-iωt}` (mass-like `X < 0`),
  refused above `ka = 1`. Same numbers as `fs_duct::Termination`.
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

- ODE-form only: constrained/DAE Dirac structures (rigid
  interconnection, Kirchhoff laws) are DEFERRED with the stated
  trigger (first consumer needing constraints).
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
