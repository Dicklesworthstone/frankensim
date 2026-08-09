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
  central-difference Jacobian, state-relative tolerance, stagnation
  acceptance at the best iterate). Returns `StepRecord` with the exact
  ledger: `delta_h`, `dissipated = dt dg^T R dg`, `supplied = dt u^T
  y`. `balance_residual()` restates the update equation (solver
  diagnostic); `supply_defect() = delta_h - supplied` is the
  INDEPENDENT passivity audit (`<= 0` for a true pHS).
- `interconnect` — canonical skew pairing `u_a = -y_b`, `u_b = +y_a`
  over chosen port pairs; unpaired ports stay external (a's first);
  the composite is RE-ADMITTED (no false passivity). Associative on
  structure matrices for disjoint pairings.
- `reduce_galerkin` — `V^T J V` / `V^T R V` / `V^T G` with
  `H_r(x_r) = H(V x_r)`; skewness and PSD-ness survive Galerkin
  automatically and are re-admitted anyway. v1 requires quadratic
  storage (probed and refused otherwise); the reduced energy deficit
  at t = 0 equals the energy outside the basis exactly.
- Zoo: `mass_spring_damper`, `lc_ladder`, `modal_bank`
  (mass-normalized modes -> canonical pHS: the bridge from the
  eig/plate beads), `duffing_oscillator` (non-quadratic exercise).

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

`tests/reed_casebook.rs` (1): single-reed exciter as pHS components —
reed lamella (msd pHS) + modal bore bank + Bernoulli valve as a
memoryless dissipative port (`dp * U >= 0` asserted every step);
per-step component supply audits + closed whole-instrument energy
accounting bounded by mouth supply; quasi-static equilibrium pinned
against the hand-derived analytic solution (`q* = Sr pm / k`, zero
bore DC pressure, `U* = w h* sqrt(2 pm / rho)`). JSON-lines summary.

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
  cycles) belongs to the distributed-contact/exciter follow-up
  (q6nmy); the casebook validates the quasi-static regime where an
  exact analytic answer exists.
- The bead's "reproduces the existing reed solver" premise was STALE:
  no reed solver exists in-tree (verified 2026-08-09); the analytic
  quasi-static pin replaces it.
- fs-time strategy wiring is a follow-up owned at L3 (fs-time depends
  down, not up).
