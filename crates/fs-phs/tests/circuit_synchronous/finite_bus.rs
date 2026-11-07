//! Actual converter, finite buses and nonlinear machine in one native solve.
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::synchronous::{
    SynchronousMotor, SynchronousMotorParameters, phase_voltage_to_alpha_beta,
};
use fs_phs::{InterconnectionBudget, SolveBudget, StepWorkspace, interconnect_mapped};

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * scale,
        "{actual:.17e} != {expected:.17e} (scale {scale:e})"
    );
}

#[test]
fn held_three_phase_bridge_motors_and_regenerates_on_finite_native_buses() {
    for direction in [-1.0, 1.0] {
        for generating in [false, true] {
            let bus = BuckBoostParameters {
                inductance_h: 0.1,
                bus_capacitance_f: [0.2, 0.3],
                winding_resistance_ohm: 0.04,
                input_on_resistance_ohm: [0.0; 2],
                output_on_resistance_ohm: [0.0; 2],
            };
            let p = SynchronousMotorParameters {
                d_inductance_h: 0.05,
                q_inductance_h: 0.08,
                winding_resistance_ohm: 0.02,
                pm_flux_linkage_wb: 0.8,
                pole_pairs: 2,
                rotor_inertia_kg_m2: 0.2,
                stator_inertia_kg_m2: 0.3,
                bearing_drag_nm_s_per_rad: 0.001,
            };
            let converter = AveragedBuckBoost::new(
                bus,
                BridgeModulation {
                    input_high_fraction: 0.75,
                    output_high_fraction: 1.0,
                },
            )
            .unwrap();
            let machine = SynchronousMotor::new(p).unwrap();
            let bus_initial = converter.state(0.0, [4.0, 3.0]).unwrap();
            let speed = if generating { direction * 10.0 } else { 0.0 };
            let machine_initial = machine
                .state(
                    [0.0; 2],
                    -core::f64::consts::FRAC_PI_4,
                    [
                        speed,
                        -p.rotor_inertia_kg_m2 / p.stator_inertia_kg_m2 * speed,
                    ],
                )
                .unwrap();
            // Actual complementary bridge states 100 and 011. Their common
            // mode has no current path in this balanced three-wire machine.
            let gates = if direction > 0.0 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 1.0]
            };
            let modulation = phase_voltage_to_alpha_beta(gates).unwrap();
            let system = interconnect_mapped(
                converter.into_system(),
                machine.into_system(),
                &[1],
                &[0, 1],
                &modulation,
                InterconnectionBudget {
                    max_states: 8,
                    max_total_ports: 6,
                    max_dense_bytes: 1 << 20,
                    max_coupling_multiply_terms: 1_000_000,
                },
            )
            .unwrap();
            assert_eq!((system.state_dim(), system.port_dim()), (8, 3));
            let mut initial = [0.0; 8];
            initial[..3].copy_from_slice(&bus_initial);
            initial[3..].copy_from_slice(&machine_initial);
            // Independent storage reconstruction, including rotor-angle
            // dependent magnetic energy and both actual finite inertias.
            let energies = |x: &[f64; 8]| {
                let (s, c) = (f64::from(p.pole_pairs) * x[5]).sin_cos();
                let d = c * x[3] + s * x[4] - p.pm_flux_linkage_wb;
                let q = -s * x[3] + c * x[4];
                [
                    0.5 * (x[0] * x[0] / bus.inductance_h
                        + x[1] * x[1] / bus.bus_capacitance_f[0]
                        + x[2] * x[2] / bus.bus_capacitance_f[1]),
                    0.5 * (d * d / p.d_inductance_h + q * q / p.q_inductance_h),
                    0.5 * (x[6] * x[6] / p.rotor_inertia_kg_m2
                        + x[7] * x[7] / p.stator_inertia_kg_m2),
                ]
            };
            let before = energies(&initial);
            near(system.hamiltonian(&initial), before.iter().sum(), 30.0);
            let mut workspace = StepWorkspace::new(&system).unwrap();
            workspace
                .set_solve_budget(SolveBudget {
                    relative: 1e-11,
                    absolute: 1e-13,
                    max_iterations: 30,
                })
                .unwrap();
            workspace.set_backtracking_limit(8).unwrap();
            let mut state = initial;
            let mut heat = 0.0;
            let dt = 0.001;
            for _ in 0..20 {
                let (mut next, mut outputs) = ([0.0; 8], [0.0; 3]);
                // No prescribed bus voltage: surviving ports are source
                // current, rotor torque and stator torque, all held at zero.
                let record = workspace
                    .step_into_controlled(
                        &system,
                        &state,
                        &[0.0; 3],
                        dt,
                        &mut next,
                        &mut outputs,
                        || false,
                    )
                    .unwrap();
                let effort = workspace.last_step_discrete_effort().unwrap();
                let source_heat = dt
                    * (bus.winding_resistance_ohm * effort[0].powi(2)
                        + p.winding_resistance_ohm * (effort[3].powi(2) + effort[4].powi(2))
                        + p.bearing_drag_nm_s_per_rad * (effort[6] - effort[7]).powi(2));
                near(record.dissipated, source_heat, 1.0);
                assert!(record.dissipated >= 0.0);
                assert_eq!(record.supplied, 0.0);
                near(record.delta_h + record.dissipated, 0.0, 30.0);
                // Composition must preserve the machine's physical component
                // discrete gradient, including exact quadratic shaft efforts.
                near(
                    outputs[1],
                    0.5 * (state[6] + next[6]) / p.rotor_inertia_kg_m2,
                    20.0,
                );
                near(
                    outputs[2],
                    0.5 * (state[7] + next[7]) / p.stator_inertia_kg_m2,
                    20.0,
                );
                near(next[5] - state[5], dt * (outputs[1] - outputs[2]), 1.0);
                near(next[6] + next[7], initial[6] + initial[7], 10.0);
                // Independent physical phase-current reconstruction verifies
                // that the held gate vector and DC load are power conjugate.
                let abc = [
                    (2.0_f64 / 3.0).sqrt() * effort[3],
                    -effort[3] / 6.0_f64.sqrt() + effort[4] / 2.0_f64.sqrt(),
                    -effort[3] / 6.0_f64.sqrt() - effort[4] / 2.0_f64.sqrt(),
                ];
                let dc_load = gates
                    .iter()
                    .zip(abc)
                    .map(|(gate, current)| gate * current)
                    .sum::<f64>();
                near(
                    dc_load,
                    modulation[0] * effort[3] + modulation[1] * effort[4],
                    10.0,
                );
                near(abc.iter().sum(), 0.0, 10.0);
                heat += record.dissipated;
                state = next;
            }
            let after = energies(&state);
            near(
                after.iter().sum::<f64>() - before.iter().sum::<f64>() + heat,
                0.0,
                30.0,
            );
            assert!(heat > 0.0);
            assert_eq!(after[2] < before[2], generating);
            assert_eq!(after[0] > before[0], generating);
            assert_eq!(state[1] > initial[1], generating);
            assert_eq!(state[2] > initial[2], generating);
            if !generating {
                assert!(direction * state[6] > 0.0);
                assert!(direction * state[7] < 0.0);
            }
        }
    }
}
