# Two material gases, one material-bound mechanism, one finite thermal store

`fs_couple::cavity_phs::piston::double_acting` connects two finite gas volumes
to the same translating body. `piston::thermal::double_acting` additionally
connects both gases to one finite enthalpy/phase store. These are shared physical
states, not separately advanced pistons or two withdrawals from an unchanged
reservoir.

## One owner for each physical law

The existing `ClosedGasParcel::midpoint_heat_volume` supplies each gas's complete
source-resolved caloric, pressure, transport and acoustic midpoint/endpoint.
Composition and supported caloric data remain inputs; there is no nominal air,
constant-gamma substitution or second gas model in the coupling.

`MaterialDoubleActingGasPiston` retains `fs_mbd_material::ResolvedBody`. Its
canonical spatial inertia supplies the only moving mass. Regional geometry and
density records are not discarded after extracting a number. An existing
`MaterialGasPiston` can become double-acting through `with_opposing_gas` without
resetting its primary gas state or accepted momentum.

The original `EquilibriumEnthalpyPhaseCurve` supplies body temperature and phase
from enthalpy. Both faces use the same enthalpy midpoint. A latent plateau can
therefore accept energy while changing phase rather than inventing a temperature
rise. `MaterialConductionLayers` supplies both interface conductances from the
actual gas conductivity, source-bound body conductivity and caller geometry.
Missing properties, unsupported states and mismatched body cards remain errors,
including on an explicitly insulated face.

Changing a supported density card before resolving the mechanical body changes
acceleration and both gas trajectories. Changing a supported body conductivity
card changes heat transfer, gas pressure, common motion and the updated acoustic
states. These are different initialized material scenarios, not unaccounted
changes of live mass or chemical energy.

## Reciprocal equations

Positive displacement expands the primary gas and compresses the opposing gas.
The positive projected face areas may differ:

```text
dVa = Aa * dx                    dVb = -Ab * dx
dx  = dt * (v0 + v1) / 2
M * (v1-v0) = dt * (Aa*pa_mid - Ab*pb_mid + F_external)

dUa = Ha + Qa - pa_mid*dVa       dUb = Hb + Qb - pb_mid*dVb
m_body * (h1-h0) = H_body - Qa - Qb
Qi = dt * Gi * (T_body_mid - T_gas_i_mid), i in {a,b}
```

`Ha`, `Hb` and `H_body` are prescribed integrated heat sources in joules. `Qa`
and `Qb` are internal exchanges, not additional external supplies. The combined
balance independently checks

```text
dUa + dUb + dK + m_body*(h1-h0) - Ha - Hb - H_body - F_external*dx
```

There is one kinetic storage and one thermal-body debit. Opposing pressure work
is internal. Both original discrete pressure derivatives enter the mechanical
Newton tangent. The existing bounded thermal feedback solver enforces both
interface residuals, with an additional `maximum_joint_trials` cap across all
nested proposals. Energy closure is not proof of entropy production, source
calibration, uniqueness or unconditional timestep stability.

## Source-driven calling pattern

The inputs below are supplied by the existing material, geometry and control
owners. No per-material physics callback is necessary for series conduction:

```rust
use fs_couple::cavity_phs::piston::material::MaterialDoubleActingGasPiston;
use fs_couple::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;

let actuator = MaterialDoubleActingGasPiston::new(
    primary_gas, opposing_gas, resolved_body,
    [primary_area_m2, opposing_area_m2], initial_velocity_m_s,
)?;
let candidate = actuator.thermal_layers_candidate(
    gas_heat_and_force_input,
    finite_body_enthalpy_input,
    paired_gas_and_momentum_controls,
    DoubleHeatExchangeSolve {
        faces: [primary_heat_controls, opposing_heat_controls],
        maximum_joint_trials: total_trial_budget,
    },
    [primary_material_layers, opposing_material_layers],
    &mut cancelled,
)?;
let (next_actuator, thermal_step, balances) = candidate.into_parts();
// Publish next_actuator and thermal_step.state together in the owning simulation.
// Do not debit the returned boundary heat again.
```

For other interface laws, `thermal_candidate` accepts one repeatable resolver
that sees both gas midpoints and the shared body phase state. Both returned
source evaluations survive in `thermal_step.report`. Callbacks must not debit
live state. Source, geometry, balance, work-budget and cancellation failures
return no partially advanced system. The material-bound wrapper retains the
original specimen and checks cancellation again before returning it.

## Boundaries and verification

The gases are uniform, closed and frozen-composition; the moving specimen is
rigid, and areas are constant. The finite thermal store has fixed geometry and
negligible pressure-volume/deformation work. It may be a separate wall or store;
its phase change does not silently deform or re-resolve the moving specimen.
There are no valves, combustion, mass exchange, contact stops or crankshafts in
this component. It does not yet connect the shared mechanism into the existing
independent-wall electrical/thermal network or a live audio renderer. Acoustic
consumers must use the exposed updated gas states.

Twenty-four added Rust regressions cover shared motion, material inertia, both
heat laws, phase storage, source-layer substitution, original-owner parity,
repeated energy balance, invalid sources and final cancellation. They are in
the existing piston target:

```sh
cargo test -p fs-couple --lib cavity_phs::piston:: -- --test-threads=1
```

Native Rust validation is not claimed. The authoring environment has no Rust
toolchain, and the observed GitHub run stopped this target at dependency
resolution (`fsci-special = 0.4.0` required, pinned sibling `0.2.0` available).
Independent synthetic reference-equation checks passed, including 72 paired
thermal cases compared with a separate two-variable root solve and timestep
refinement. Those checks neither execute the Rust code nor calibrate real
material behavior.
