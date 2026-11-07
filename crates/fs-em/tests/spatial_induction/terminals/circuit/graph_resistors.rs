//! Source-backed resistor updates preserve the original field and RLC storage.
use super::super::thermal::{source, temperature_problem, temperature_sources};
use super::*;
use fs_em::winding::nonlinear::resistors::GraphResistorMaterial;
use fs_material::conductor::ElectricalPropertyBasis;

fn resistor(
    card: &MaterialCard,
    temperature: f64,
    branch: usize,
    length: f64,
) -> GraphResistorMaterial {
    let basis = ElectricalPropertyBasis::Resistivity;
    let bundle = source::bundle_at(card, &source::absolute_point(temperature), basis);
    GraphResistorMaterial {
        branch,
        conductor: ResolvedConductorState::try_from_material_with_basis(bundle, basis)
            .unwrap()
            .into_uniform_conductor(length, 1.0)
            .unwrap(),
    }
}

#[test]
fn resistor_temperature_and_field_material_refresh_preserve_rlc_history_and_separate_heat() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let field_card = source::absolute_affine_magnetic_card(0.2, 0.002, 2.0, 0.003, basis);
    let resistor_card = source::absolute_affine_magnetic_card(0.4, 0.001, 1.0, 0.0, basis);
    let (problem, sources) = temperature_problem(&field_card, basis);
    let (inductance, capacitance, dt): (f64, f64, f64) = (0.3, 0.2, 0.03);
    let (mut a, mut flux, mut charge, mut k0) = (0.0, 0.05, 0.14, 2.0);
    let mut circuit = SpatialTerminalCircuit::new(
        terminal(&problem, sources, [0.0; 3]),
        CircuitGraph {
            node_count: 4,
            branches: vec![
                (3, 0, Branch::VoltageSource { port: 0 }),
                (3, 2, Branch::Resistor { ohms: 9.0 }),
                (
                    2,
                    1,
                    Branch::Inductor {
                        henries: inductance,
                    },
                ),
                (
                    1,
                    0,
                    Branch::Capacitor {
                        farads: capacitance,
                    },
                ),
            ],
            transformers: vec![],
        },
        vec![1, 0],
        vec![],
        initial(vec![flux], vec![charge], vec![0.0]),
        config(0.7),
        &mut || false,
    )
    .unwrap();
    assert!(circuit.graph_resistor_materials().is_empty());
    for (field_temperature, resistor_temperature, source_voltage) in [
        (300.0, 300.0, 1.2),
        (600.0, 800.0, -0.4),
        (450.0, 350.0, 0.8),
        (800.0, 600.0, 1.1),
    ] {
        let graph_sources = vec![resistor(&resistor_card, resistor_temperature, 1, 1.0)];
        let resistance = 0.4 * (1.0 + 0.001 * (resistor_temperature - 300.0));
        near(
            graph_sources[0].conductor.resistance_ohm(),
            resistance,
            1e-14,
        );
        let retained_native = circuit.native_state().to_vec();
        let retained_field = circuit.field().potential_wb().to_vec();
        let rebound = circuit
            .with_graph_resistors(graph_sources.clone(), &mut || false)
            .unwrap();
        assert_eq!(rebound.native_state(), retained_native);
        assert_eq!(rebound.field().potential_wb(), retained_field);
        assert_eq!(rebound.time_s(), circuit.time_s());
        assert_eq!(rebound.terminal_nodes(), circuit.terminal_nodes());
        assert_eq!(rebound.graph_resistor_materials(), graph_sources);
        let (materials, conductors) = temperature_sources(&field_card, field_temperature, basis);
        let sigma = 1.0 / conductors[0].as_ref().unwrap().resistivity_ohm_m();
        let (f, k1) = (
            sigma / 6.0,
            4.0 / (2.0 + 0.003 * (field_temperature - 300.0)),
        );
        let d = sigma / 5.0 + 0.5 * dt * k1;
        let y = sigma - f * f / d;
        let i0 = f * (k0 + k1) * a / (2.0 * d);
        // Independent midpoint KVL/KCL with the endpoint material R(T).
        let z = resistance + 2.0 * inductance / dt;
        let u = (source_voltage + 2.0 * flux / dt - z * (i0 - 2.0 * charge / dt))
            / (1.0 + z * (y + 2.0 * capacitance / dt));
        let incoming = y * u + i0;
        let branch_current = incoming + 2.0 * (capacitance * u - charge) / dt;
        let (old_a, old_flux, old_charge) = (a, flux, charge);
        a += dt * (f * u - 0.5 * (k0 + k1) * a) / d;
        flux = 2.0 * inductance * branch_current - flux;
        charge = 2.0 * capacitance * u - charge;
        let candidate = rebound
            .step_with_materials_candidate(
                &[source_voltage],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                materials,
                conductors,
                assembly(),
                circuit_controls(),
                &mut || false,
            )
            .unwrap();
        let step = candidate.step();
        near(
            candidate.circuit().field().potential_wb()[diagonal(&problem)],
            a,
            3e-8,
        );
        near(step.terminal_voltage_v[0], u, 3e-8);
        near(step.terminal_current_a[0], incoming, 3e-8);
        near(
            step.circuit_stored_before_j,
            old_flux * old_flux / (2.0 * inductance)
                + old_charge * old_charge / (2.0 * capacitance),
            3e-8,
        );
        near(
            step.circuit_stored_after_j,
            flux * flux / (2.0 * inductance) + charge * charge / (2.0 * capacitance),
            3e-8,
        );
        near(step.parameter_work_j, 0.5 * (k1 - k0) * old_a * a, 3e-8);
        near(
            step.resistor_heat_j,
            dt * resistance * branch_current * branch_current,
            3e-8,
        );
        near(
            step.external_work_j,
            dt * source_voltage * branch_current,
            3e-8,
        );
        near(step.heat_j, step.field.heat_j + step.resistor_heat_j, 3e-8);
        near(step.energy_residual_j, 0.0, 3e-8);
        assert_eq!(step.graph_resistors.len(), 1);
        assert_eq!(step.graph_resistors[0].branch, 1);
        near(step.graph_resistors[0].heat_j, step.resistor_heat_j, 3e-8);
        assert_eq!(step.graph_materials, graph_sources);
        assert_eq!(
            candidate.circuit().graph_resistor_materials(),
            graph_sources
        );
        assert_eq!(circuit.native_state(), retained_native);
        assert_eq!(circuit.field().potential_wb(), retained_field);
        circuit = candidate.into_parts().0;
        k0 = k1;
    }
}

#[test]
fn identical_resistor_bindings_reuse_response_and_changed_resistance_requires_rebuild() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.4, 0.001, 1.0, 0.0, basis);
    let (problem, sources) = cube(1, None);
    let cold = vec![resistor(&card, 300.0, 1, 1.0)];
    let state = SpatialTerminalCircuit::new(
        terminal(&problem, sources, [0.0; 3]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (2, 0, Branch::VoltageSource { port: 0 }),
                (2, 1, Branch::Resistor { ohms: 0.4 }),
            ],
            transformers: vec![],
        },
        vec![1, 0],
        vec![],
        initial(vec![], vec![], vec![0.0]),
        config(1.0),
        &mut || false,
    )
    .unwrap()
    .with_graph_resistors(cold.clone(), &mut || false)
    .unwrap()
    .step_candidate(
        &[1.2],
        &[[0.0; 3]; 6],
        &[[0.0; 3]; 8],
        0.05,
        circuit_controls(),
        &mut || false,
    )
    .unwrap()
    .into_parts()
    .0;
    let mut bounded = circuit_controls();
    bounded.maximum_field_trials = 2;
    let unchanged = state
        .with_graph_resistors(cold.clone(), &mut || false)
        .unwrap();
    let cached = unchanged
        .step_candidate(
            &[1.2],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.05,
            bounded,
            &mut || false,
        )
        .unwrap();
    assert!(cached.step().reused_response);
    assert_eq!(cached.step().graph_materials, cold);
    let hot = vec![resistor(&card, 800.0, 1, 1.0)];
    let changed = state
        .with_graph_resistors(hot.clone(), &mut || false)
        .unwrap();
    assert!(matches!(
        changed.step_candidate(
            &[1.2],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.05,
            bounded,
            &mut || false,
        ),
        Err(SpatialTerminalCircuitError::FieldTrialBudget)
    ));
    let actual = changed
        .step_candidate(
            &[1.2],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.05,
            circuit_controls(),
            &mut || false,
        )
        .unwrap();
    let old_a = state.field().potential_wb()[diagonal(&problem)];
    let d = 1.05;
    let y = 5.0 - (5.0 / 6.0_f64).powi(2) / d;
    let offset = (5.0 / 6.0) * 2.0 * old_a / d;
    let r = hot[0].conductor.resistance_ohm();
    let voltage = (1.2 - r * offset) / (1.0 + r * y);
    near(actual.step().terminal_voltage_v[0], voltage, 3e-8);
    near(
        actual.step().terminal_current_a[0],
        y * voltage + offset,
        3e-8,
    );
    assert!(!actual.step().reused_response);
    assert_eq!(actual.step().field_trials, 3);
    assert_eq!(actual.circuit().graph_resistor_materials(), hot);
    assert_eq!(state.graph_resistor_materials(), cold);
    assert_eq!(state.time_s(), 0.05);
}

#[test]
fn shorted_electrodes_retain_sorted_resistor_sources_and_reject_incomplete_maps() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.4, 0.001, 1.0, 0.0, basis);
    let (problem, sources) = cube(1, None);
    let state = SpatialTerminalCircuit::new(
        terminal(&problem, sources, [0.0; 3]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (2, 0, Branch::VoltageSource { port: 0 }),
                (2, 1, Branch::Resistor { ohms: 9.0 }),
                (1, 0, Branch::Resistor { ohms: 9.0 }),
            ],
            transformers: vec![],
        },
        vec![0, 0],
        vec![],
        initial(vec![], vec![], vec![0.0]),
        config(1.0),
        &mut || false,
    )
    .unwrap();
    let sorted = vec![
        resistor(&card, 300.0, 1, 1.0),
        resistor(&card, 300.0, 2, 2.0),
    ];
    let bound = state
        .with_graph_resistors(vec![sorted[1].clone(), sorted[0].clone()], &mut || false)
        .unwrap();
    assert_eq!(bound.graph_resistor_materials(), sorted);
    let before = bound.native_state().to_vec();
    let mut controls = circuit_controls();
    controls.maximum_field_trials = 1;
    let candidate = bound
        .step_candidate(
            &[1.2],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.1,
            controls,
            &mut || false,
        )
        .unwrap();
    assert_eq!(candidate.step().field_trials, 1);
    assert_eq!(candidate.step().graph_materials, sorted);
    assert_eq!(candidate.circuit().graph_resistor_materials(), sorted);
    near(candidate.step().graph_resistors[0].heat_j, 0.04, 3e-9);
    near(candidate.step().graph_resistors[1].heat_j, 0.08, 3e-9);
    near(candidate.step().resistor_heat_j, 0.12, 3e-9);
    near(candidate.step().parameter_work_j, 0.0, 1e-15);
    assert!(
        bound
            .with_graph_resistors(vec![sorted[0].clone()], &mut || false)
            .is_err()
    );
    assert!(
        bound
            .with_graph_resistors(vec![sorted[0].clone(), sorted[0].clone()], &mut || false)
            .is_err()
    );
    assert!(
        bound
            .with_graph_resistors(
                vec![resistor(&card, 300.0, 0, 1.0), sorted[1].clone()],
                &mut || false
            )
            .is_err()
    );
    assert!(
        bound
            .with_graph_resistors(sorted.clone(), &mut || true)
            .is_err()
    );
    assert_eq!(bound.native_state(), before);
    assert_eq!(bound.graph_resistor_materials(), sorted);
    assert_eq!(bound.time_s(), 0.0);
}
