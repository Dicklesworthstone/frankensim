# Endpoint heat exchange for constant-pressure chemistry

The public APIs live under
`cavity_phs::piston::reactive::implicit::coupled::isobaric::thermal`.
They extend the original `GasReaction` and `ClosedGasParcel` owners; they do not
introduce another equation of state, chemistry integrator or phase model.

## Shared endpoint, signed contacts

`GasReaction::isobaric_thermal_ports_candidate` closes the existing
constant-pressure implicit mechanism against signed thermal contacts. Every
heat proposal starts at the original accepted gas and solves

```
H_gas1 = H_gas0 + Q_gas_external + Q
p_gas1 = p_gas0
n1 = n0 + S xi
xi_r = dt V1 rate_r(T1, p0, n1)
Q = sum_i dt G_i(T_port_i - T_gas1)
```

The SAME endpoint composition, temperature, pressure and volume determine all
reaction rates and contact laws. This is nonlinear endpoint feedback, not a
reaction-then-cooling split. Formation energy is already in gas enthalpy and is
never supplied again as heat. Callbacks are repeatable material/interface laws,
not permission to mutate an external reservoir during a numerical trial.

Ports are canonicalized by original ID, with compensated heat and entropy
reductions. Every `ThermalPortSample` retains its signed heat and concrete
source result, including opposing transfers at zero net gas heat. Missing
sources refuse even on an isolated contact. A finite receiving network must
apply the reciprocal transfers and publish all histories together. The returned
gas already contains Q; depositing it again double-counts energy. Per-contact
law values and their aggregate source residual remain explicit.

Independent acceptance checks the heat equation and both

```
Delta H_gas - Q_gas_external - Q_evaluated
Delta U_gas + p0 (V1-V0) - Q_gas_external - Q_evaluated
```

The contact-entropy field sums the actual endpoint laws
`Q_i (1/T_gas1 - 1/T_port_i)` before heat cancellation. It is not an entropy
account for unspecified reservoirs or a global passivity theorem.

## Reciprocal finite body and latent heat

`GasReaction::isobaric_thermal_body_candidate` additionally evolves a
`UniformEnthalpyStepInput` on its ORIGINAL equilibrium phase chart:

```
h_body1 = h_body0 + (Q_body_external - Q) / m_body
T_body1, phase_body1 = original_chart(h_body1)
Q = dt G(gas1, body1) (T_body1 - T_gas1)
```

The body is not an infinite temperature reservoir. Its temperature and phase
fraction feed back into the gas chemistry and the interface source at the same
endpoint. A latent plateau can change fraction while its temperature remains
fixed; no apparent heat capacity or extra latent-energy credit is invented.

`IsobaricBodyCandidate` returns the complete gas, original `UniformEnthalpyStep`
and joint report together. The latter independently checks

```
Delta H_gas + Delta H_body - Q_gas_external - Q_body_external
Delta U_gas + p0 Delta V + Delta H_body - Q_gas_external - Q_body_external
```

The body boundary heat is exactly `-Q`; it is already included in the returned
body state. Both original histories remain unchanged on source refusal,
unsupported phase state, unresolved heat or chemistry, exhausted work, or
cancellation, including late rejection of the joint entropy account.

The body entropy is integrated by the original phase chart as `integral dh/T`,
including `Delta h/T_transition` on a latent plateau. `IsobaricBodySolve` supplies
a separate absolute entropy allowance. The report distinguishes

```
Delta S_body - (Q_body_external-Q)/T_body1
Q (1/T_gas1 - 1/T_body1)
Delta S_gas + Delta S_body - Q_gas_external/T_gas1 - Q_body_external/T_body1
```

The applied-Q contact entropy is distinct from the independently evaluated
contact-law entropy. Negative values outside the declared allowance refuse.
These endpoint residuals include backward-Euler numerical entropy; they are not
isolated chemical production or measured heat-reservoir entropy. A physical
external reservoir needs its own entropy account.

`GasReaction::isobaric_thermal_layers_candidate` reuses
`MaterialConductionLayers`. Gas conductivity comes from the actual reacting
endpoint; body conductivity and phase storage must use the SAME material card.
Declared thermal/fraction axes and original source selection survive. The
existing series resistance is `G = A/(d_gas/k_gas + d_body/k_body)`; the patch
cannot exceed the body's exposed area. This remains a reduced conduction model,
not a resolved through-wall profile or an inferred convection correlation.

## Work, ownership and scope

`IsobaricHeatSolve` retains the original heat bracket, absolute residual
tolerances and scalar iteration cap. It adds a strict cumulative
`maximum_chemical_evaluations` across accepted and discarded heat trials. Only
the small numerical policy is cloned once. Each original Newton solve receives
the lesser of its original evaluation limit and the remaining total budget,
so the final heat trial cannot overspend. Thermochemical source models are not
cloned merely to rebind that numerical limit. Existing per-inversion and
per-mechanism work bounds remain in force.

This is first-order, closed homogeneous ideal-gas chemistry under a quasistatic
constant-pressure load. The finite body is uniform, stationary/fixed-geometry,
with negligible pressure-volume work. Its melting does not implicitly change
geometry, inertia, flow mobility or mechanical supports. No flame, open-flow
reactor, pressure-controller dynamics, experimental material qualification,
phase-stability certificate or arbitrary-step convergence claim is added.
Nonlinear feedback may have multiple roots. Every queried bracket proposal
must remain in the original gas, body and source domains; temporal refinement
and adequate material data remain required.

## Verification

```
cargo test -p fs-couple --test reactive_isobaric_heat
```

Eleven authored Rust regression groups exercise the production owners:
constant-Cp pressure-work/heat, canonical opposing contacts, exothermic endpoint
rate feedback, exact cumulative work limits, missing sources/cancellation,
finite sensible-wall energy/entropy, latent fraction changes, exothermic finite
wall feedback, analytical first-order exchange refinement, original layer
source identity/conductivity, and finite-body domain/work/late rollback.

**Native Rust compilation, Cargo tests and rustfmt were not executed in this
authoring environment: the Rust toolchain is unavailable.** The following
executed independent Python equations are supplementary reference checks, not
execution of the Rust implementation:

| Independent check | Observed result |
| --- | --- |
| 300 nonlinear reactor/contact endpoints, nested versus simultaneous equations | Maximum energy/source residual 3.80e-14 J; reaction residual 3.19e-17 mol; endpoint discrepancy at most 2.28e-13 K or mol |
| 300 nonlinear reacting finite-wall endpoints | Maximum total energy residual 1.09e-12 J; heat-law residual 2.46e-14 J; all endpoint heat-entropy residuals nonnegative |
| Latent-wall reference | 1.9774011299435 J supplied to gas; wall remains at 600 K while liquid fraction decreases from 0.5 to 0.4980225988701 |
| Exothermic finite-wall trajectory versus DOP853, 32/64/128 steps | Error-reduction ratios 1.99675 and 1.99838 |
| Independent two-capacity exchange, 8/16/32 steps | Error-reduction ratios 1.98036 and 1.99006 |

No native throughput benchmark, full-workspace pass or cross-ISA qualification
is claimed. Energy convention: Cantera 3.2 constant-pressure reactor equation,
https://www.cantera.org/3.2/reference/reactors/constant-pressure-reactor.html .
