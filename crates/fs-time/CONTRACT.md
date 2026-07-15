# fs-time CONTRACT

## Purpose and layer

Layer: **L3 FLUX** (deps: fs-ad, fs-ga, fs-la, fs-math — all L0–L2).
Structure-preserving time integration (plan §8.5): integrators that
preserve what the physics preserves — symplectic (Störmer–Verlet with
its discrete-Lagrangian equivalence tested), Lie-group SO(3) via
exponential-map updates and SE(3) on fs-ga PGA motors (exp-map lane
plus a free-body discrete Euler–Poincaré variational lane),
generalized-α with controllable dissipation, IMEX and exponential
integrators for stiffness, and embedded-pair adaptivity with a PI
controller. The two universal obligations (P7 + §8.7) — resumable
state machines and discrete adjoints of the stepper — are discharged
where claimed below.

## Public types and semantics

- `symplectic::verlet_step(q, p, h, force, scratch)` — one
  kick–drift–kick step for separable H = ½‖p‖² + V(q), unit mass
  (scale q/p externally otherwise). `force` writes F(q) = −∇V.
- `symplectic::verlet_adjoint(q0, p0, h, steps, force, force_jvp,
  (bar_q_N, bar_p_N)) -> (bar_q0, bar_p0)` — discrete adjoint OF THE
  STEPPER (transposed tangent map, reverse order), propagated with
  `fs_ad::revolve` checkpointing (O(log N) memory).
  **Precondition: symmetric ∂F/∂q** (conservative forces): the reverse
  pass evaluates Jᵀv through the user's `force_jvp` as Jᵀv = Jv.
- `lie::{quat_mul, quat_exp, quat_exp_step, quat_rotate}` — (w, x, y, z)
  unit quaternions; `quat_exp` uses a series branch below θ = 1e−6 (no
  sinc cancellation). `lie::rigid_body_step(q, ω, I, h)` — torque-free
  Euler equations via midpoint RK2 + CG2 exp-map attitude update at the
  midpoint ω; returns (q′, ω′).
- `galpha::GeneralizedAlpha::new(m, c, k, n, h, rho_inf)` +
  `galpha_step(ga, q, v, a, f_next)` — Chung–Hulbert: αm = (2ρ−1)/(ρ+1),
  αf = ρ/(ρ+1), γ = ½ − αm + αf, β = ¼(1 − αm + αf)²; one prefactored
  LU of the effective matrix per (M, C, K, h); `f_next` is the load at
  t + (1−αf)h. Newmark correctors update (q, v, a) in place.
- `stiff::Imex2::new(l, n, h)` + `imex2_step` — ARS(2,2,2),
  γ = 1 − 1/√2, diagonally implicit on L (one LU reused for both
  stages), explicit N, stiffly accurate, R(∞) = 0. Second order in
  BOTH parts: the explicit weights are (δ, 1−δ) with δ = 1 − 1/(2γ);
  trapezoidal (½, ½) is only first order in N (caught during
  construction, locked by the convergence test).
- `stiff::ExpEuler::new(a, n, h)` + `.step(u, nonlin)` — exponential
  Euler for u′ = Au + N(u), **symmetric A** via the fs-la Jacobi
  eigenbasis; φ₁(x) = expm1(x)/x (cancellation-free). Exact for N ≡ 0.
- `adaptive::{AdaptiveState, PiController, rk45_adaptive}` —
  Dormand–Prince 5(4), max-norm error against atol + rtol·|u|, PI
  update h ← h·safety·err^{−0.14}·err_prev^{0.08} clamped to
  [0.2, 5.0]; rejections shrink with the classical −1/5 exponent.
- `se3` module (bead 3ol0): `Twist` (body ω then v);
  `se3_exp_step(M, twist, h)` — `M ← M ∘ exp(−(h/2)·B(ω, v))` through
  fs-ga's screw exponential, returning the CANONICAL double-cover
  representative (`canonicalize_motor`: first nonzero even component
  in `EVEN_BLADES` order made positive; `M` and `−M` canonicalize
  bit-identically; all-zero refuses). `se3_exp_step_renorm` adds
  drift-controlled versor renormalization gated by `RenormPolicy` and
  returns a `RenormReceipt` (defect before, whether renormalized,
  reported drift) — drift is ledger fodder, never silently absorbed.
  `se3_rigid_body_step` — free-body Euler equations + body transport
  of the spatially constant velocity, midpoint RK2 in the algebra with
  one exp update at the midpoint twist. `dep_free_step` — discrete
  Euler–Poincaré free-body step in body-momentum form (fixed-point
  solve controlled by `DepSolveParams`, per-step `DepStepReceipt`);
  spatial angular momentum is conserved EXACTLY by construction.
  `run_dep_free` returns a `BalanceReceipt` whose `Se3ClaimClass` is
  decided by `claim_for(declaration, all_solves_converged)`: the
  conservative variational theorem class only for declared smooth,
  conservative, regular-constraint, fixed-step fixtures with every
  solve converged — dissipation, adaptivity, or divergence demote to
  `MeasuredOnly` with measured drift receipts. `dep_momentum_adjoint`
  pulls a terminal cotangent through the transposed implicit-function
  tangent of each step's ACTUAL residual. `RattleProjection` is the
  constraint hook fs-mbd's constrained lanes plug into
  (`Unconstrained` is the trivial impl).

- `slabs` module (addendum Proposal 4, bead bk0o.7; [F], behind
  `time-slabs`): time slabs as CELLS. `SlabEntry` carries the TEMPORAL
  COCYCLE — the split step's defect against the monolithic residual
  over the slab — and `SlabLedger::attribute` is the budget pie
  pointed at time ("your error is in the coupling handoff at
  t ∈ [2.0, 2.5]"). `march_adaptive` doubles subcycles where the
  cocycle exceeds tolerance (cap 64). `activation_report` encodes the
  Proposal-4 SEQUENCING gate: splitting error under 20% of budget →
  INSTRUMENT-ONLY (measure it, don't control it). Slab constructors
  and marchers fail fast on zero substeps/slabs, invalid time spans,
  non-finite states/couplings, and non-positive total error budgets;
  ledger reporting rejects hand-built malformed entries before
  attribution or JSON emission.

## Invariants

- Verlet is symplectic ⇒ **bounded** (non-secular) energy error;
  its KDK positions satisfy the discrete Euler–Lagrange recurrence
  q_{k+1} = 2q_k − q_{k−1} + h²F(q_k) (Marsden–West).
- Exp-map attitude updates keep ‖q‖ = 1 by construction — no
  renormalization anywhere in the crate.
- Generalized-α: high-frequency per-step contraction → ρ∞; order 2
  across the whole ρ∞ range.
- `AdaptiveState` is COMPLETE: checkpoint = `clone()`; split runs are
  **bitwise-equal** to straight runs, controller memory (`err_prev`)
  and counters included (P7). A step shortened to land exactly on
  t_end does not feed the controller (the clamped h would poison the
  h carried into a later resumed segment).

## Error model

Construction panics with a structured message when a required
factorization is singular (`GeneralizedAlpha::new`, `Imex2::new`):
these are modeling errors (h incompatible with the operators), not
runtime conditions. `verlet_adjoint` inherits `fs_ad::revolve`'s
budget assertion. Steppers themselves are panic-free on finite input;
NaN/Inf propagate as NaN/Inf (garbage-in, garbage-out, never UB).

## Determinism class

Bit-deterministic across ISAs BY CONSTRUCTION: no platform libm in any
solver path (workspace law) — transcendentals go through
`fs_math::det` (`sin`, `cos`, `sqrt`, `exp`, `expm1`, `pow`), including
the PI controller's `pow`, because the adaptive step SEQUENCE is part
of the contract. Test-side oracles may use std (disjoint-path rule).
Golden FNV-64 over Verlet, rigid-body, generalized-α, IMEX, ExpEuler
and RK45 trajectories (controller state and counters included):
`0xeae8_ccec_5e2e_cf41`, recorded on Apple M4 Pro (aarch64), verified
identical on Threadripper (x86_64).

## Cancellation behavior

All entry points are synchronous, allocation-light, and run to
completion; long integrations are resumable via `AdaptiveState`
(interrupt between calls, `clone()` to checkpoint, continue bitwise).
No async, no internal threading, no I/O.

## Unsafe boundary

None. `unsafe_code = "deny"` via workspace lints; no capsules.

## Feature flags

- `time-slabs` — [F] time slabs as temporal cells and splitting-error
  ledger/controller instrumentation. No optional dependencies.

## Conformance tests

`tests/time_battery.rs` (16 cases, JSON logging): 10⁶-step harmonic
energy boundedness vs RK4's secular decay at the same h + e = 0.6
Kepler orbit (~16 revolutions); discrete Euler–Lagrange residual
≤ 1e−12 on a nonlinear potential; quaternion norm drift ≈ 1e−12 over
10⁵ steps with no renormalization; gyroscope battery (ω₃ constant to
1e−9, analytic precession phase Ω = (I₃−I₁)/I₁·ω₃ to 1e−3, energy and
spatial angular momentum to O(h²)); generalized-α spectral radius → ρ∞
at ωh = 10³ to 2% with ρ∞ = 0 annihilation and ρ∞ = 1 energy
preservation; IMEX hλ = −100 monotone contraction + measured order 2
on the logistic equation; ExpEuler exact (1e−13) on a linear system
against a disjoint-path oracle + measured ETD1 order 1; RK45 accuracy
tracking rtol, rejection recovery from an absurd h₀, and bitwise
split-run resumability at 4 cut points; Verlet adjoint gradcheck vs
central FD ≤ 1e−7 relative; the cross-ISA golden hash.
`tests/galpha_probe.rs`: order-2 sweep across ρ∞ ∈ {0, .3, .5, .8, 1}
at two horizons — kept from the probe that diagnosed the period-point
metric blindness (at t = 2π, cos′ = 0: a q-only error measures phase
error quadratically and fakes order ≈ 4; the honest metric is
max(q, v) error).
`tests/slabs.rs` under `time-slabs`: temporal cocycle consistency,
budget-pie localization, adaptive-vs-uniform cost, G3 repartition
envelope, activation gate, and fail-fast validation for invalid slab
counts/substeps/tolerances/budgets, malformed public ledger entries,
and non-finite couplings.
`tests/se3.rs` (bead 3ol0, printed measurements on every gate):
SO(3)-lane agreement for pure rotations (se3-001); constant-twist
one-parameter composition exactness (se3-002); double-cover
canonicalization determinism, `M`/`−M` bit-equality, and bitwise
replay across a scalar-zero crossing (se3-003); free-body DEP spatial
momentum to roundoff plus bounded energy over 10⁴ steps with the
theorem claim class earned (se3-004); adjoint-vs-central-FD gate on
the DEP momentum map, four directions, rel ≤ 1e−6 (se3-005); the
honesty fixture — a damped run demotes to `MeasuredOnly` with nonzero
measured drift despite converged solves (se3-006); 10⁵-step
renormalization receipts bounding the final unit defect (se3-007);
SE(3) rigid-body agreement with the SO(3) lane and spatial
free-velocity drift at the measured-order level (se3-008).

## No-claim boundaries

- SE(3) motor states ship the exp-map lane and the FREE-BODY
  variational lane only: no forced/damped/constrained variational
  steps here (constrained lanes are fs-mbd's, through
  `RattleProjection`); the backward-error/modified-energy THEOREM is
  claimed only for the declared smooth conservative fixed-step class
  with converged solves — impacts, adaptivity, and dissipation get
  measured `BalanceReceipt`s, never the theorem. The
  `dep_momentum_adjoint` residual Jacobians are central differences of
  the actual residual (3×3), not analytic tangents; the analytic
  derivation and revolve checkpointing (O(N) memory today) are
  follow-up work. The variational lane's group update solves the
  midpoint transported-momentum fixed point; no claim of equivalence
  to any OTHER published DEP variant. No higher-order compositions
  (Yoshida/Suzuki), no splitting beyond kick–drift–kick.
- No Krylov φ-actions for large nonsymmetric A (needs Arnoldi);
  `ExpEuler` is dense-symmetric only. No Rosenbrock/SDIRK families,
  no BDF/multistep.
- No dense output / continuous extension for RK45; no stiffness
  detection; no event location.
- Adjoints ship for Verlet only (the template); generalized-α/IMEX/
  RK45 adjoints are the fs-ad integration lane (o3ui).
- `Imex2`/`GeneralizedAlpha` take dense row-major operators; sparse
  variants belong to the fs-sparse integration lane.

## No-claim boundaries (slabs)

- Visibility and control of SPLITTING error only — NO claim of
  coupling STABILITY: added-mass FSI instabilities and stiff
  time-scale pathologies are per-coupling analysis problems (the
  proposal's own honest scope, verbatim).
- The monolithic reference is a fine RK4 over the slab — a numerical
  reference, not an exact solution; defects below its own error floor
  are not resolvable.
- The fixture family is the linear two-field testbed; PDE couplings
  ride fs-couple's port-Hamiltonian interconnection when both land.
- Parallel-in-time coarse propagators (the BDDC-pattern extension) are
  explicitly deferred.
