//! The native RC/RLC states and both magnetic endpoints survive source refresh.
use super::super::thermal::{source, temperature_problem, temperature_sources};
use super::*;
use fs_material::conductor::ElectricalPropertyBasis;

#[test]
fn capacitor_discharge_refresh_retains_charge_and_invalidates_only_the_old_response() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.002, 2.0, 0.003, basis);
    let (problem, sources) = temperature_problem(&card, basis);
    let (capacitance, dt): (f64, f64) = (0.4, 0.04);
    let mut circuit = SpatialTerminalCircuit::new(
        terminal(&problem, sources, [0.0; 3]),
        CircuitGraph {
            node_count: 2,
            branches: vec![(
                1,
                0,
                Branch::Capacitor {
                    farads: capacitance,
                },
            )],
            transformers: vec![],
        },
        vec![1, 0],
        vec![],
        initial(vec![], vec![capacitance], vec![]),
        config(0.7),
        &mut || false,
    )
    .unwrap();
    let (mut a, mut charge, mut k0) = (0.0, capacitance, 2.0);
    for (index, (temperature, refresh)) in [
        (300.0, false),
        (600.0, true),
        (600.0, false),
        (350.0, true),
        (350.0, false),
    ]
    .into_iter()
    .enumerate()
    {
        let (materials, conductors) = temperature_sources(&card, temperature, basis);
        let sigma = 1.0 / conductors[0].as_ref().unwrap().resistivity_ohm_m();
        let (f, k1) = (sigma / 6.0, 4.0 / (2.0 + 0.003 * (temperature - 300.0)));
        // The old material enters the history offset; only the endpoint
        // conductivity and stiffness enter the finite-stage response.
        let d = sigma / 5.0 + 0.5 * dt * k1;
        let y = sigma - f * f / d;
        let i0 = f * (k0 + k1) * a / (2.0 * d);
        let incoming = (y * charge / capacitance + i0) / (1.0 + y * dt / (2.0 * capacitance));
        let voltage = charge / capacitance - dt * incoming / (2.0 * capacitance);
        let old_a = a;
        let old_charge = charge;
        a += dt * (f * voltage - 0.5 * (k0 + k1) * a) / d;
        charge -= dt * incoming;
        let before = circuit.native_state().to_vec();
        if refresh {
            let mut bounded = circuit_controls();
            bounded.maximum_field_trials = 2;
            assert!(matches!(
                circuit.step_with_materials_candidate(
                    &[],
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    dt,
                    materials.clone(),
                    conductors.clone(),
                    assembly(),
                    bounded,
                    &mut || false,
                ),
                Err(SpatialTerminalCircuitError::FieldTrialBudget)
            ));
        }
        let candidate = if refresh {
            circuit.step_with_materials_candidate(
                &[],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                materials,
                conductors,
                assembly(),
                circuit_controls(),
                &mut || false,
            )
        } else {
            circuit.step_candidate(
                &[],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                dt,
                circuit_controls(),
                &mut || false,
            )
        }
        .unwrap();
        let step = candidate.step();
        assert_eq!(step.reused_response, index != 0 && !refresh);
        assert_eq!(step.field_trials, if step.reused_response { 2 } else { 3 });
        near(
            candidate.circuit().field().potential_wb()[diagonal(&problem)],
            a,
            3e-9,
        );
        near(step.terminal_voltage_v[0], voltage, 3e-9);
        near(step.terminal_current_a[0], incoming, 3e-9);
        near(
            step.circuit_stored_before_j,
            old_charge * old_charge / (2.0 * capacitance),
            3e-9,
        );
        near(
            step.circuit_stored_after_j,
            charge * charge / (2.0 * capacitance),
            3e-9,
        );
        near(step.field.stored_before_j, 0.5 * k0 * old_a * old_a, 3e-9);
        near(
            step.field.parameter_work_j,
            0.5 * (k1 - k0) * old_a * a,
            3e-9,
        );
        near(step.field.stored_after_j, 0.5 * k1 * a * a, 3e-9);
        near(step.energy_residual_j, 0.0, 3e-9);
        assert_eq!(circuit.native_state(), before);
        near(
            circuit.field().potential_wb()[diagonal(&problem)],
            old_a,
            3e-9,
        );
        circuit = candidate.into_parts().0;
        k0 = k1;
    }
}

#[test]
fn thermal_rlc_closure_retains_original_flux_charge_and_external_source_work() {
    let basis = ElectricalPropertyBasis::Conductivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.002, 2.0, 0.003, basis);
    let (problem, sources) = temperature_problem(&card, basis);
    let (resistance, inductance, capacitance, dt): (f64, f64, f64, f64) = (0.4, 0.3, 0.2, 0.03);
    let graph = CircuitGraph {
        node_count: 4,
        branches: vec![
            (3, 0, Branch::VoltageSource { port: 0 }),
            (3, 2, Branch::Resistor { ohms: resistance }),
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
    };
    let field = terminal(&problem, sources, [0.0; 3]);
    let mut endpoints = Vec::new();
    for reference in [0.01, 1.0, 100.0] {
        let (mut a, mut flux, mut charge, mut k0) = (0.0, 0.05, 0.14, 2.0);
        let mut circuit = SpatialTerminalCircuit::new(
            field.clone(),
            graph.clone(),
            vec![1, 0],
            vec![],
            initial(vec![flux], vec![charge], vec![0.0]),
            config(reference),
            &mut || false,
        )
        .unwrap();
        for (temperature, source_voltage) in
            [(300.0, 1.2), (650.0, -0.4), (450.0, 0.8), (800.0, 1.1)]
        {
            let (materials, conductors) = temperature_sources(&card, temperature, basis);
            let sigma = 1.0 / conductors[0].as_ref().unwrap().resistivity_ohm_m();
            let (f, k1) = (sigma / 6.0, 4.0 / (2.0 + 0.003 * (temperature - 300.0)));
            let d = sigma / 5.0 + 0.5 * dt * k1;
            let y = sigma - f * f / d;
            let i0 = f * (k0 + k1) * a / (2.0 * d);
            // Native midpoint KVL: (R+2L/dt)*I_L + U = E+2*flux0/dt.
            // Native KCL: I_L = (Y+2C/dt)*U + I0 - 2*q0/dt.
            let z = resistance + 2.0 * inductance / dt;
            let voltage = (source_voltage + 2.0 * flux / dt - z * (i0 - 2.0 * charge / dt))
                / (1.0 + z * (y + 2.0 * capacitance / dt));
            let incoming = y * voltage + i0;
            let source_current = incoming + 2.0 * (capacitance * voltage - charge) / dt;
            let (old_a, old_flux, old_charge) = (a, flux, charge);
            a += dt * (f * voltage - 0.5 * (k0 + k1) * a) / d;
            flux = 2.0 * inductance * source_current - flux;
            charge = 2.0 * capacitance * voltage - charge;
            let candidate = circuit
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
            assert!(!step.reused_response);
            assert_eq!(step.field_trials, 3);
            near(step.terminal_voltage_v[0], voltage, 3e-8);
            near(step.terminal_current_a[0], incoming, 3e-8);
            near(step.branches[1].current_a, source_current, 3e-8);
            near(step.branches[2].current_a, source_current, 3e-8);
            near(step.branches[3].current_a, (charge - old_charge) / dt, 3e-8);
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
            near(step.field.stored_before_j, 0.5 * k0 * old_a * old_a, 3e-8);
            near(step.field.stored_after_j, 0.5 * k1 * a * a, 3e-8);
            near(
                step.field.parameter_work_j,
                0.5 * (k1 - k0) * old_a * a,
                3e-8,
            );
            near(
                step.resistor_heat_j,
                dt * resistance * source_current * source_current,
                3e-8,
            );
            near(
                step.external_work_j,
                dt * source_voltage * source_current,
                3e-8,
            );
            near(step.energy_residual_j, 0.0, 3e-8);
            circuit = candidate.into_parts().0;
            k0 = k1;
        }
        endpoints.push(circuit.field().potential_wb()[diagonal(&problem)]);
    }
    for &a in &endpoints[1..] {
        near(a, endpoints[0], 3e-8);
    }
}

#[test]
fn shorted_circuit_refresh_uses_one_field_solve_and_refusals_preserve_both_histories() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.002, 2.0, 0.003, basis);
    let (problem, sources) = temperature_problem(&card, basis);
    let field = terminal(&problem, sources, [0.0; 3])
        .step_voltage_candidate(
            &[1.0, 0.0],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.1,
            controls(),
            &mut || false,
        )
        .unwrap()
        .into_parts()
        .0;
    let circuit = SpatialTerminalCircuit::new(
        field,
        CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, Branch::Resistor { ohms: 0.7 })],
            transformers: vec![],
        },
        vec![0, 0],
        vec![],
        initial(vec![], vec![], vec![]),
        config(1.0),
        &mut || false,
    )
    .unwrap();
    let old_native = circuit.native_state().to_vec();
    let old_field = circuit.field().potential_wb().to_vec();
    let mut bounded = circuit_controls();
    bounded.maximum_field_trials = 1;
    let (materials, conductors) = temperature_sources(&card, 600.0, basis);
    let candidate = circuit
        .step_with_materials_candidate(
            &[],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.1,
            materials.clone(),
            conductors.clone(),
            assembly(),
            bounded,
            &mut || false,
        )
        .unwrap();
    assert_eq!(candidate.step().field_trials, 1);
    assert_eq!(candidate.step().field.terminal_work_j, 0.0);
    assert!(candidate.step().field.parameter_work_j < 0.0);
    near(candidate.step().energy_residual_j, 0.0, 3e-9);
    let mut missing = conductors.clone();
    missing[0] = None;
    assert!(
        circuit
            .step_with_materials_candidate(
                &[],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.1,
                materials.clone(),
                missing,
                assembly(),
                bounded,
                &mut || false,
            )
            .is_err()
    );
    assert!(
        circuit
            .step_with_materials_candidate(
                &[],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.1,
                materials,
                conductors,
                assembly(),
                bounded,
                &mut || true,
            )
            .is_err()
    );
    assert_eq!(circuit.native_state(), old_native);
    assert_eq!(circuit.field().potential_wb(), old_field);
    assert_eq!(circuit.time_s(), 0.1);
}
