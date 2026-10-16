# Total fan-speed sensitivity and bounded target sizing

The experimental `cooling-network` command now connects a steady temperature
objective to fan speed through the declared fan, quadratic hydraulic network,
transported heat capacity, Reynolds-dependent convection and shared solid.
It reuses the existing implicit coupled adjoint, including matching contacts
and the material K'(T) Jacobian. It does not finite-difference the solver or
hold the air temperatures artificially fixed.

```bash
frankensim --json cooling-network examples/cooling-network/size-fan-adjoint-hotspot.json
```

The example has a one-watt localized load, two materials, a cooling branch and
bypass which merge before a downstream exchange. Its inputs are illustrative,
not measured hardware. Remove `fan_speed_design` to evaluate just its declared
speed. `objective.gradient=true` requests the existing conditional controls and
the new total fan-speed response together.

## What the result means

`fan_speed_sensitivity.status="available"` carries:

- `dobjective_dlog_speed_ratio_k`: total dJ/dln(s), in kelvin.
- `capacity_contribution_k`: air-capacity response at fixed coefficients.
- `convection_contribution_k`: sum of the total dJ/dln(h) controls times the
  actual convection cards' dln(Nu)/dln(Re). The two contributions sum to total.

For dJ/ds, divide the logarithmic derivative by the actual returned fan speed
ratio, not by an unrelated base-request speed. A small fractional speed change
uses dJ approximately equal to the logarithmic derivative times that fraction.
This is a local numerical derivative, not a bound on finite changes or safety.
For discrete maxima, it retains the existing objective's active-vertex/tie
semantics; it does not assert smoothness across a hotspot switch.

A gradient-disabled, transient, or prescribed-pressure result has a null fan
sensitivity. The rectangular developing table currently reports `unavailable`
instead of inventing a smooth slope or quietly treating its h as constant.
The original conditional inlet/h/contact gradients remain available. Smooth
supported cards are the circular and rectangular fully developed CWT limits,
Hausen, Dittus-Boelter, and Gnielinski, within the original validity domains.

## Why a complete fan derivative needs no extra solid solve

The admitted single fan bank obeys affinity scaling p(q,s)=s^2 p(q/s,1).
Every passive branch obeys delta-p=R q|q| at fixed R. Hence all signed branch
flows scale with s, pressures scale with s^2, and mixing fractions do not
change. Every branch and external capacity rate and each Reynolds number
therefore scale with s. This remains true for reversed edge orientations and
fixed series/parallel bank configurations.

At fixed walls, scaling both capacity and hA by the same factor preserves
transported temperatures and scales watt outputs by that factor. The air-only
pullback, seeded with the converged interface adjoint, therefore yields the
capacity derivative by subtracting its log(hA) controls from its explicit
heat-functional term. Subtracting the total solid-plus-air h controls would be
wrong. The extension performs one additional air-only reverse sweep, not an
additional solid solve. Its general library heat functionals retain the
explicit watt term even when a numerical heat residual is merely small.

This does not implement derivatives of arbitrary independent branch losses,
flow-dependent resistance, multiple independently controlled fans, geometry,
fan heat, fluid properties, or transient trajectories.

## Target search uses derivatives only as proposals

When an available gradient was requested, the existing fan target search tries
a log-speed Newton proposal strictly inside the central 80% of its evaluated
passing/failing bracket. Otherwise it uses the original midpoint. Thus either
surviving bracket contracts, and a zero, positive, unavailable, nonfinite or
out-of-bracket slope never sends a speculative design outside that bracket.

Every proposed speed reruns the real hydraulics, correlations and solid/air
solve. Model refusals are errors, not hot/cold observations. Both the original
speed-width and passing-temperature tolerances must hold. Results retain the
actual passing field, failed lower endpoint, full history and `newton_trials`;
per-trial history now includes the evaluated logarithmic speed derivative.
There are no hidden adjoints when `objective.gradient=false`, no increased
budgets, and no promise that requesting derivatives reduces total compute.
A passing bracket is not a global minimum-speed or hardware certificate.

## Verification scope

```bash
cargo test -p fs-airflow --test mixed_network_flow_scale
cargo test -p fs-cli --bin frankensim network_command::fan_gradient
cargo test -p fs-cli --bin frankensim network_command::fan_speed
cargo test -p fs-cli --test cooling_fan_sensitivity
```

Added tests compare the real graph/FEM/transport producers against an
independent slab/NTU result, preserve explicit heat-functional scaling, and
refuse interrupted/unconverged adjoints. Actual-binary tests perturb speed in
the complete cooling model, retain contact/nonlinear-material feedback, reverse
an edge, and re-evaluate the selected design with gradients disabled.

Independent NumPy P1 FEM checks covered eighteen steady cases, with worst
central-difference discrepancy 4.97e-8 K per log-speed change and worst
chain-rule versus direct implicit discrepancy 1.25e-14 K. The baseline
correlated hotspot reference is approximately 302.457593084 K, with total
log-speed derivative -0.600305330 K (capacity -0.116572917 K and convection
-0.483732413 K). These are reference calculations, not Rust executions.
On this example at a 302.4 K target, the independent search mirror used 21
complete evaluations with three Newton proposals versus 19 for bisection;
the new derivative capability is not presented as a measured speedup.

The authoring environment did not provide a Rust toolchain. Compilation and
Rust test execution remain unverified. The workflow remains a nominal,
fixed-discretization JSON product, not physical validation or a native .fsim
ledger-backed study. Executable-bound UQ checkpoints from a different binary
continue to refuse rather than silently mixing old and new sample producers.
