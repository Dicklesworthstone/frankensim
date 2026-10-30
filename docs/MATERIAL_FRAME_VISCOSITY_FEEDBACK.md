# Prepared source-viscosity feedback for material frames

`MaterialFluidFrame::prepare` retains the original full Hamiltonian, loaded
modal inertia, source sections, and every fluid-memory momentum. It uses the
existing native `StepWorkspace` with analytic Hessian/resistance actions and
exact independent-memory condensation. No modes or fluid coordinates are
removed. The factored resistance callback evaluates each station velocity
once, then projects its accumulated reaction once. Native algebra retains
its cancellation, residual, iteration and safe full-LU fallback gates.

Preparation allocates. Native stepping and its physical callbacks allocate
nothing after preparation; returned reports and material evaluations can
allocate. No elapsed-time, real-time, Rust benchmark or full-workspace
qualification is claimed. Complete source-frame state/dense limits still
apply, and the native solver additionally owns its bounded workspace.

`LiquidViscosityFrame::new` binds the exact original `LiquidTemperatureLaw`
for every loaded member. `step_at_temperatures` evaluates those pinned claims
at prescribed temperatures and changes the same physical resistances.
A validity-wide scalar constant-density source is mandatory: matching two
density samples is not enough. At fixed density/radius and dimensionless
spectrum, c_j/lambda_j is viscosity-independent. The original kinetic storage
and all momenta therefore remain unchanged when viscosity varies.
Dynamic and kinematic source-viscosity conventions use the existing material
owner; unsupported source domains and unowned rate/frequency laws refuse.

```rust,ignore
let mut prepared = model.prepare(&gate)?;
let next = prepared.step_candidate(&state, &nodal_forces, dt, energy_budget_j, &gate)?;
state = next.state; // only after the enclosing transaction admits the proposal

let mut viscous = LiquidViscosityFrame::new(&model, original_liquid_laws, &gate)?;
let next = viscous.step_at_temperatures(
    &state, &nodal_forces, &prescribed_temperatures_k, dt, energy_budget_j, &gate,
)?;
```

Prescribed temperature is NOT a solved finite thermal history. In particular,
assigning every unconfined Stokes loss to an arbitrary finite liquid volume
would not resolve spatial fluid heating. This path makes no such claim.
Variable density, thermal expansion, phase change, finite boundaries and
spatially varying fluid temperature require their corresponding physical owners.

Five native Rust regressions cover dense/condensed parity, source-only
viscosity changes, cold-hot operator comparison, density-source refusal and
retry after source/cancellation failures. They are authored but UNEXECUTED:
Rust, Cargo, RCH and DSR are unavailable in this editing environment.
Independent NumPy/SciPy matrix equations passed ten one-/four-mode cases over
viscosity ratios 0.01 through 100. Maximum normalized resistance-action
mismatch was 6.66e-17, dense/Schur state discrepancy 2.78e-17, and work residual
6.35e-17 J. These reference calculations do not execute or compile Rust.
