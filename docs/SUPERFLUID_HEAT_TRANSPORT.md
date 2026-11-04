# Superfluid heat transport and vortex friction

The `superfluid` feature connects the finite-temperature Landau two-fluid
material, the original two-sound channel, native rigid pistons, and the existing
finite thermal timeline. A heater can now launch counterflow and second sound
while its own finite caloric storage changes the heat it supplies.

## Run the finite heater example

```bash
cargo run -p fs-couple --features superfluid --example superfluid_heater
```

This uses the existing sourced helium-II point at **1.8 K and saturated vapor
pressure**. It prints heater temperature, channel temperature at two locations,
heat flow, and independent reference-heat balances. The heater capacitance,
interface conductance, geometry, source power and frozen vortex population are
explicit demonstration inputs, not fitted measurements of an apparatus.

The example propagates both sound modes, retains thermal expansion, and applies
isotropic vortex friction on the original channel. Exported friction heat and
wave availability are printed separately from heat transported by counterflow.

## Finite caloric contact

`fs_couple::superfluid::counterflow::advance_counterflow_heat` connects a selected
channel end to an active `fs_thermal::enthalpy::CellId`. The opposite end is
closed. The contact has an authored **total conductance in W/K**, including area.
Both endpoints carry zero net mass current; heat flow is positive from the
channel into the finite thermal cell.

For each interval it:

1. Prepares the original channel's reciprocal midpoint temperature response.
2. Advances the original thermal sources, links and scheduled events once.
3. Solves the interface law `Q = G * (T_channel_mid - T_cell_mid)` against the
   original cell's caloric chart, with bounded scalar probes. The channel's
   wave impedance and the finite cell's heat capacity both feed back into Q.
4. Applies exactly `Q * dt` through native `apply_heat_j`, preserving original
   material domains and any latent plateau.
5. Publishes both histories only after physical limits, transfer checks and
   cancellation admission succeed.

The native thermal network followed by the contact is a **first-order split**;
refine the timestep for coupling accuracy. The contact itself uses midpoint
temperatures. Probes do not repeatedly solve the wave equation or clone the
thermal network: they query its original caloric chart at candidate energies.
An out-of-chart trial is shortened while bracketing; it never extrapolates.

The source and sink must remain in the caller's declared small-signal range
about the channel's frozen reference temperature. No nonlinear helium equation
of state, phase crossing, evolving material coefficients, or uncalibrated
Kapitza conductance is inferred.

## Three distinct thermal quantities

For reference temperature `T0` and channel temperature perturbation `theta`:

| Quantity | Outward rate | Meaning |
| --- | --- | --- |
| Reference-state counterflow heat | `Q` [W] | Heat exchanged with the caloric owner |
| Reference entropy flow | `Q / T0` [W/K] | Linearized entropy transport |
| Wave-availability power | `theta * Q / T0` [W] | Boundary power in the quadratic wave model |

These rates are not interchangeable. The report separately retains native
caloric energy, reference heat and entropy residuals, the wave port residual,
and the contact's linearized availability loss. That availability loss is not
extra heat to deposit. The coupled reference-heat account is not a nonlinear
total-internal-energy theorem for helium or an arbitrary enthalpy-based solid.

Optional mutual-friction losses remain in the channel's original per-cell heat
field and report. This boundary contact does not redirect all distributed
frictional heating into the heater. Existing closed-channel and piston thermal
compositions remain available for explicitly mapped friction-heat destinations.

## Thermal boundary response for other couplings

`TwoFluidChannel::counterflow_boundary_velocities([Q_left, Q_right])` converts
authored total watts into opposite normal and superfluid velocities.
`counterflow_ends_response` returns the full reciprocal 2x2 midpoint-temperature
response, retaining propagation between both ends, both sound modes, body forces,
thermal expansion and the original friction law. It reuses the cached native
tridiagonal factors and endpoint columns in O(cells) work and storage.

The response is an algebraic proposal. It never publishes a field, clock or heat
history, and does not replace the final step's physical amplitude admission.

## Vortex-dependent mutual friction

`fs_material::superfluid::VortexMutualFriction` requires explicit circulation
magnitude, vortex-line density, dissipative Hall–Vinen B and reactive B-prime.
`TwoFluidState::relax_hall_vinen` analytically advances homogeneous flow through
a frozen aligned straight-vortex population. It leaves parallel counterflow
unchanged, damps and rotates transverse counterflow, preserves barycentric
momentum, and exports heat from dissipative transverse motion only. Pure
reactive rotation produces no heat.

The scalar channel supports `MutualFriction::IsotropicVortexTangle`. Its force
average is `gamma = (rho_n*rho_s/rho) * B*kappa*L/3`; it uses the existing exact
linear-drag step and heat destinations. A one-dimensional channel does not
resolve oriented transverse vortex motion. Neither API evolves or nucleates
vortices, models tension/reconnections, or predicts a critical velocity.

The coefficient conventions and angular average follow
[Jou, Mongiovi and Sciacca, equations 2.5–2.9](https://arxiv.org/abs/0806.1889).
The helium reference and B/alpha convention follow
[Donnelly and Barenghi (1998), section 11](https://doi.org/10.1063/1.556028).
The latter source documents frequency and velocity dependence; the API therefore
does not silently assign a universal B value to a material name.

## Focused regression targets

```bash
cargo test -p fs-material --test superfluid
cargo test -p fs-duct --features superfluid --test superfluid_channel --test superfluid_ports
cargo test -p fs-couple --features superfluid --test superfluid_counterflow --test superfluid_heat
```

The contact tests include independent finite-capacity/two-cell equations,
heating and cooling at either end, nonlinear caloric feedback, original sources,
latent storage, reference heat and linear availability balances, zero-contact
identity, and atomic refusal/cancellation/retry. Vortex regressions cover
direction, frame covariance, opposite impulses, decay, reactive rotation,
weak/strong damping, the angular force average and channel equivalence.
