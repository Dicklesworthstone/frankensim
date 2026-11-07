//! Original RLC storage and physical electrode currents on native 3D fields.
use super::*;
use fs_em::spatial::vector_potential::induction::terminals::circuit::{
    SpatialTerminalCircuit, SpatialTerminalCircuitConfig, SpatialTerminalCircuitControls,
    SpatialTerminalCircuitError, SpatialTerminalCircuitInitial,
};
use fs_phs::circuit::{Branch, CircuitGraph};

#[path = "circuit/graph_resistors.rs"]
mod graph_resistors;
#[path = "circuit/thermal.rs"]
mod thermal;

fn config(reference_ohm: f64) -> SpatialTerminalCircuitConfig {
    SpatialTerminalCircuitConfig {
        reference_resistance_ohm: reference_ohm,
        maximum_native_states: 128,
        maximum_interface_scalars: 100_000,
    }
}
fn circuit_controls() -> SpatialTerminalCircuitControls {
    SpatialTerminalCircuitControls {
        field: controls(),
        voltage_probe_v: 0.5,
        current_probe_a: 0.25,
        maximum_field_trials: 18,
        relative_terminal_tolerance: 2e-10,
        voltage_tolerance_v: 5e-8,
        current_tolerance_a: 5e-8,
        energy_tolerance_j: 2e-7,
    }
}
fn terminal<'m>(
    problem: &'m SpatialVectorMagnetostatics,
    sources: Vec<Option<ResolvedConductorState>>,
    b: [f64; 3],
) -> SpatialTerminalInduction<'m> {
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let initial = problem
        .solve(
            &zero,
            &trace(problem, b, false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let native =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    SpatialTerminalInduction::new(native, &electrodes(problem), preparation(), &mut || false)
        .unwrap()
}
fn initial(
    fluxes_wb: Vec<f64>,
    charges_c: Vec<f64>,
    external_inputs: Vec<f64>,
) -> SpatialTerminalCircuitInitial {
    SpatialTerminalCircuitInitial {
        fluxes_wb,
        charges_c,
        external_inputs,
    }
}

#[test]
fn native_series_resistor_drives_exact_field_transient_and_reference_cancels() {
    // Independent one-edge FE reduction: adot = (-K a + f U)/(M+K dt/2),
    // I = sigma U - f adot. Kirchhoff closes U + R I = E. The field
    // self-induction remains a dynamic state; it is not a DC resistor model.
    let (problem, sources) = cube(1, None);
    let field = terminal(&problem, sources, [0.0; 3]);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::VoltageSource { port: 0 }),
            (2, 1, Branch::Resistor { ohms: 0.4 }),
        ],
        transformers: vec![],
    };
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let (sigma, mass, stiffness, f, dt, source): (f64, f64, f64, f64, f64, f64) =
        (5.0, 1.0, 2.0, 5.0 / 6.0, 0.05, 1.2);
    let mut outcomes = vec![];
    for reference in [0.001, 1.0, 1000.0] {
        let mut circuit = SpatialTerminalCircuit::new(
            field.clone(),
            graph.clone(),
            vec![1, 0],
            vec![],
            initial(vec![], vec![], vec![0.0]),
            config(reference),
            &mut || false,
        )
        .unwrap();
        let mut a = 0.0;
        for step in 0..24 {
            let d = mass + stiffness * dt / 2.0;
            let y = sigma - f * f / d;
            let i0 = f * stiffness * a / d;
            let u = (source - 0.4 * i0) / (1.0 + 0.4 * y);
            let current = y * u + i0;
            a += dt * (-stiffness * a + f * u) / d;
            let candidate = circuit
                .step_candidate(
                    &[source],
                    &zero,
                    &velocity,
                    dt,
                    circuit_controls(),
                    &mut || false,
                )
                .unwrap();
            let result = candidate.step();
            near(result.terminal_voltage_v[0], u, 2e-8);
            near(result.terminal_current_a[0], current, 2e-8);
            near(
                candidate.circuit().field().potential_wb()[diagonal(&problem)],
                a,
                2e-8,
            );
            near(result.resistor_heat_j, dt * 0.4 * current * current, 2e-8);
            near(result.external_work_j, dt * source * current, 2e-8);
            near(result.field.terminal_work_j, dt * u * current, 2e-8);
            near(result.energy_residual_j, 0.0, 2e-8);
            assert_eq!(result.reused_response, step != 0);
            assert_eq!(result.field_trials, if step == 0 { 3 } else { 2 });
            assert!(result.maximum_kirchhoff_current_defect_a < 5e-8);
            circuit = candidate.into_parts().0;
        }
        let mut warm_only = circuit_controls();
        warm_only.maximum_field_trials = 2;
        assert!(
            circuit
                .step_candidate(&[source], &zero, &velocity, dt, warm_only, &mut || false)
                .is_ok()
        );
        assert!(matches!(
            circuit.step_candidate(
                &[source],
                &zero,
                &velocity,
                dt / 2.0,
                warm_only,
                &mut || false
            ),
            Err(SpatialTerminalCircuitError::FieldTrialBudget)
        ));
        let mut different_field_accuracy = warm_only;
        different_field_accuracy.field.absolute_residual_a *= 0.5;
        assert!(matches!(
            circuit.step_candidate(
                &[source],
                &zero,
                &velocity,
                dt,
                different_field_accuracy,
                &mut || false
            ),
            Err(SpatialTerminalCircuitError::FieldTrialBudget)
        ));
        outcomes.push(circuit.field().potential_wb()[diagonal(&problem)]);
    }
    for &a in &outcomes[1..] {
        near(a, outcomes[0], 3e-8);
    }
}

#[test]
fn original_capacitor_discharge_and_field_storage_follow_independent_midpoint_equations() {
    let (problem, sources) = cube(1, None);
    let field = terminal(&problem, sources, [0.0; 3]);
    let capacitance = 0.4;
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(
            1,
            0,
            Branch::Capacitor {
                farads: capacitance,
            },
        )],
        transformers: vec![],
    };
    let mut circuit = SpatialTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        initial(vec![], vec![capacitance], vec![]),
        config(0.7),
        &mut || false,
    )
    .unwrap();
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let (dt, sigma, stiffness, f): (f64, f64, f64, f64) = (0.02, 5.0, 2.0, 5.0 / 6.0);
    let (mut a, mut charge, mut heat) = (0.0, capacitance, 0.0);
    for _ in 0..32 {
        let d = 1.0 + stiffness * dt / 2.0;
        let y = sigma - f * f / d;
        let i0 = f * stiffness * a / d;
        let current = (y * charge / capacitance + i0) / (1.0 + y * dt / (2.0 * capacitance));
        let voltage = charge / capacitance - dt * current / (2.0 * capacitance);
        charge -= dt * current;
        a += dt * (-stiffness * a + f * voltage) / d;
        let candidate = circuit
            .step_candidate(&[], &zero, &velocity, dt, circuit_controls(), &mut || false)
            .unwrap();
        let result = candidate.step();
        near(result.terminal_voltage_v[0], voltage, 2e-8);
        near(result.terminal_current_a[0], current, 2e-8);
        near(
            result.circuit_stored_after_j,
            charge * charge / (2.0 * capacitance),
            2e-8,
        );
        near(result.field.stored_after_j, stiffness * a * a / 2.0, 2e-8);
        near(result.external_work_j, 0.0, 1e-15);
        near(result.resistor_heat_j, 0.0, 1e-15);
        heat += result.field.heat_j;
        near(
            result.circuit_stored_after_j + result.field.stored_after_j + heat,
            capacitance / 2.0,
            2e-8,
        );
        circuit = candidate.into_parts().0;
    }
}

#[test]
fn external_inductor_retains_its_original_flux_and_does_not_replace_field_induction() {
    let (problem, sources) = cube(1, None);
    let field = terminal(&problem, sources, [0.0; 3]);
    let inductance = 0.3;
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (2, 0, Branch::VoltageSource { port: 0 }),
            (
                2,
                1,
                Branch::Inductor {
                    henries: inductance,
                },
            ),
        ],
        transformers: vec![],
    };
    let mut circuit = SpatialTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        initial(vec![0.0], vec![], vec![0.0]),
        config(1.1),
        &mut || false,
    )
    .unwrap();
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    let (dt, source, f): (f64, f64, f64) = (0.03, 1.4, 5.0 / 6.0);
    let (mut a, mut flux, mut heat, mut supplied) = (0.0, 0.0, 0.0, 0.0);
    for _ in 0..24 {
        let d = 1.0 + dt;
        let y = 5.0 - f * f / d;
        let i0 = 2.0 * f * a / d;
        let current = (y * (source + 2.0 * flux / dt) + i0) / (1.0 + 2.0 * inductance * y / dt);
        let voltage = source - 2.0 * (inductance * current - flux) / dt;
        flux = 2.0 * inductance * current - flux;
        a += dt * (-2.0 * a + f * voltage) / d;
        let candidate = circuit
            .step_candidate(
                &[source],
                &zero,
                &velocity,
                dt,
                circuit_controls(),
                &mut || false,
            )
            .unwrap();
        let result = candidate.step();
        near(result.terminal_current_a[0], current, 2e-8);
        near(result.branches[1].current_a, current, 2e-8);
        near(result.terminal_voltage_v[0], voltage, 2e-8);
        near(
            result.circuit_stored_after_j,
            flux * flux / (2.0 * inductance),
            2e-8,
        );
        near(result.field.stored_after_j, a * a, 2e-8);
        near(result.external_work_j, dt * source * current, 2e-8);
        heat += result.field.heat_j;
        supplied += result.external_work_j;
        near(
            result.circuit_stored_after_j + result.field.stored_after_j + heat,
            supplied,
            3e-8,
        );
        circuit = candidate.into_parts().0;
    }
}

#[test]
fn motion_generates_real_load_heat_and_consistent_mechanical_work() {
    let (problem, sources) = cube(2, None);
    let field = terminal(&problem, sources, [0.0, 0.0, 1.0]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::Resistor { ohms: 0.8 })],
        transformers: vec![],
    };
    let circuit = SpatialTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        initial(vec![], vec![], vec![]),
        config(0.5),
        &mut || false,
    )
    .unwrap();
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0, 0.7, 0.0]; problem.mesh().vertex_count];
    let candidate = circuit
        .step_candidate(&[], &zero, &velocity, 0.03, circuit_controls(), &mut || {
            false
        })
        .unwrap();
    let result = candidate.step();
    assert!(result.resistor_heat_j > 0.0 && result.field.heat_j > 0.0);
    assert!(result.field.mechanical_work_j < 0.0 && result.field.terminal_work_j < 0.0);
    assert!(result.terminal_voltage_v[0] < 0.0 && result.terminal_current_a[0] > 0.0);
    near(result.field.terminal_work_j, -result.resistor_heat_j, 3e-8);
    near(result.external_work_j, 0.0, 1e-15);
    near(result.energy_residual_j, 0.0, 3e-8);
    let cells: Vec<_> = (0..problem.mesh().tets.len()).collect();
    let loads = candidate
        .lorentz_loads(
            &cells,
            SpatialLorentzLoadLimits {
                max_cells: cells.len(),
                max_vertices: problem.mesh().vertex_count,
            },
            &mut || false,
        )
        .unwrap();
    near(
        0.03 * loads.power_w(&velocity, &mut || false).unwrap(),
        result.field.mechanical_work_j,
        1e-10,
    );
}

#[test]
fn native_ideal_current_and_voltage_sources_both_close_without_fictitious_loads() {
    let (problem, sources) = cube(1, None);
    let field = terminal(&problem, sources, [0.0; 3]);
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    for branch in [
        Branch::CurrentSource { port: 0 },
        Branch::VoltageSource { port: 0 },
    ] {
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, branch)],
            transformers: vec![],
        };
        let circuit = SpatialTerminalCircuit::new(
            field.clone(),
            graph,
            vec![1, 0],
            vec![],
            initial(vec![], vec![], vec![0.0]),
            config(0.37),
            &mut || false,
        )
        .unwrap();
        let candidate = circuit
            .step_candidate(
                &[0.9],
                &zero,
                &velocity,
                0.05,
                circuit_controls(),
                &mut || false,
            )
            .unwrap();
        let result = candidate.step();
        match branch {
            Branch::CurrentSource { .. } => near(result.terminal_current_a[0], 0.9, 2e-8),
            Branch::VoltageSource { .. } => near(result.terminal_voltage_v[0], 0.9, 2e-8),
            _ => unreachable!(),
        }
        near(result.resistor_heat_j, 0.0, 1e-15);
        near(result.external_work_j, result.field.terminal_work_j, 2e-8);
        near(result.energy_residual_j, 0.0, 2e-8);
        assert!(result.cancelled_reference_heat_j > 0.0);
    }
}

#[test]
fn electrically_shorted_electrodes_use_zero_voltage_without_a_native_self_loop_source() {
    let (problem, sources) = cube(2, None);
    let field = terminal(&problem, sources, [0.0, 0.0, 1.0]);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let velocity = vec![[0.0, 0.7, 0.0]; problem.mesh().vertex_count];
    let expected = field
        .step_voltage_candidate(&[0.0, 0.0], &zero, &velocity, 0.03, controls(), &mut || {
            false
        })
        .unwrap();
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::Resistor { ohms: 0.8 })],
        transformers: vec![],
    };
    // Both electrodes join node 1 by ideal wires. The separate external
    // resistor sees no net current, while real short-circuit current flows
    // between the two electrodes and heats the moving body.
    let circuit = SpatialTerminalCircuit::new(
        field,
        graph,
        vec![1, 1],
        vec![],
        initial(vec![], vec![], vec![]),
        config(0.7),
        &mut || false,
    )
    .unwrap();
    let mut one_trial = circuit_controls();
    one_trial.maximum_field_trials = 1;
    let candidate = circuit
        .step_candidate(&[], &zero, &velocity, 0.03, one_trial, &mut || false)
        .unwrap();
    let result = candidate.step();
    assert_eq!(result.field_trials, 1);
    assert!(!result.reused_response);
    assert_eq!(
        candidate.circuit().field().potential_wb(),
        expected.state().potential_wb()
    );
    for (actual, expected) in result
        .terminal_current_a
        .iter()
        .zip(&expected.step().terminals)
    {
        near(*actual, expected.incoming_current_a, 1e-12);
    }
    assert!(result.terminal_current_a[0] > 0.0 && result.field.heat_j > 0.0);
    near(result.field.terminal_work_j, 0.0, 1e-15);
    near(result.resistor_heat_j, 0.0, 1e-15);
    near(result.cancelled_reference_heat_j, 0.0, 1e-15);
    near(result.energy_residual_j, 0.0, 3e-8);
}

#[test]
fn storage_topology_budgets_bad_drives_and_cancellation_preserve_both_histories() {
    let (problem, sources) = cube(1, None);
    let field = terminal(&problem, sources, [0.0; 3]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::VoltageSource { port: 0 })],
        transformers: vec![],
    };
    for cfg in [
        config(0.0),
        SpatialTerminalCircuitConfig {
            maximum_native_states: 0,
            ..config(1.0)
        },
        SpatialTerminalCircuitConfig {
            maximum_interface_scalars: 0,
            ..config(1.0)
        },
    ] {
        assert!(
            SpatialTerminalCircuit::new(
                field.clone(),
                graph.clone(),
                vec![1, 0],
                vec![],
                initial(vec![], vec![], vec![0.0]),
                cfg,
                &mut || false
            )
            .is_err()
        );
    }
    let incompatible = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 0, Branch::Capacitor { farads: 1.0 }),
        ],
        transformers: vec![],
    };
    assert!(
        SpatialTerminalCircuit::new(
            field.clone(),
            incompatible,
            vec![1, 0],
            vec![],
            initial(vec![], vec![0.0], vec![1.0]),
            config(1.0),
            &mut || false
        )
        .is_err()
    );
    let mixed_units = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    };
    assert!(
        SpatialTerminalCircuit::new(
            field.clone(),
            mixed_units,
            vec![1, 0],
            vec![],
            initial(vec![], vec![], vec![0.0]),
            config(1.0),
            &mut || false
        )
        .is_err()
    );
    let circuit = SpatialTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        initial(vec![], vec![], vec![0.0]),
        config(1.0),
        &mut || false,
    )
    .unwrap();
    let native_before = circuit.native_state().to_vec();
    let field_before = circuit.field().potential_wb().to_vec();
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    for inputs in [vec![], vec![f64::NAN], vec![1.0, 0.0]] {
        assert!(
            circuit
                .step_candidate(
                    &inputs,
                    &zero,
                    &velocity,
                    0.05,
                    circuit_controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    let mut limited = circuit_controls();
    limited.maximum_field_trials = 2;
    assert!(matches!(
        circuit.step_candidate(&[1.0], &zero, &velocity, 0.05, limited, &mut || false),
        Err(SpatialTerminalCircuitError::FieldTrialBudget)
    ));
    let mut polls = 0;
    assert!(matches!(
        circuit.step_candidate(
            &[1.0],
            &zero,
            &velocity,
            0.05,
            circuit_controls(),
            &mut || {
                polls += 1;
                polls > 25
            }
        ),
        Err(SpatialTerminalCircuitError::Cancelled)
    ));
    assert_eq!(circuit.native_state(), native_before);
    assert_eq!(circuit.field().potential_wb(), field_before);
    assert_eq!(circuit.time_s(), 0.0);
}

#[test]
fn field_response_probes_do_not_bypass_original_recoil_operating_envelopes() {
    let (problem, sources) = cube(1, Some(0.4));
    let field = terminal(&problem, sources, [0.0; 3]);
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::VoltageSource { port: 0 })],
        transformers: vec![],
    };
    let circuit = SpatialTerminalCircuit::new(
        field,
        graph,
        vec![1, 0],
        vec![],
        initial(vec![], vec![], vec![0.0]),
        config(1.0),
        &mut || false,
    )
    .unwrap();
    let mut probe = circuit_controls();
    probe.voltage_probe_v = 1000.0;
    let zero = vec![[0.0; 3]; 6];
    let velocity = vec![[0.0; 3]; 8];
    assert!(matches!(
        circuit.step_candidate(&[0.1], &zero, &velocity, 0.1, probe, &mut || false),
        Err(SpatialTerminalCircuitError::Field(_))
    ));
    assert_eq!(circuit.time_s(), 0.0);
}
