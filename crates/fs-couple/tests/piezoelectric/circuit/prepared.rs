//! Independent electrical/mechanical ODE checks on authored material data.
use super::{loaded, near, rod, solve, source};
use fs_couple::material_circuit::MaterialCircuitError;
use fs_couple::material_circuit::piezoelectric::PiezoelectricCircuit;
use fs_couple::material_circuit::piezoelectric::prepared::PreparedPiezoelectricCircuit;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};

fn budgets() -> (LinearSplitBudget, DenseExponentialBudget) {
    (
        LinearSplitBudget {
            max_modes: 32,
            max_omega_dt: 3.0,
        },
        DenseExponentialBudget {
            max_dimension: 64,
            max_dense_bytes: 64 * 1024 * 1024,
            max_multiply_terms: 1_000_000_000,
        },
    )
}
fn prepare(model: PiezoelectricCircuit, ports: &[f64], h: f64) -> PreparedPiezoelectricCircuit {
    let (split, dense) = budgets();
    PreparedPiezoelectricCircuit::new(model, ports, h, split, dense, || false).unwrap()
}

fn rlc(d: f64) -> PiezoelectricCircuit {
    // Real source -> material resistor -> inductor -> piezo electrode,
    // with a second capacitor before the inductor and a harvesting resistor.
    PiezoelectricCircuit::new(
        CircuitGraph {
            node_count: 4,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: f64::NAN }),
                (2, 3, Branch::Inductor { henries: 0.5 }),
                (3, 0, Branch::Capacitor { farads: f64::NAN }),
                (3, 0, Branch::Resistor { ohms: f64::NAN }),
                (2, 0, Branch::Capacitor { farads: 0.3 }),
            ],
            transformers: vec![],
        },
        vec![source(4, 1.7), source(1, 0.8)],
        vec![0.03],
        vec![0.012, 0.02],
        3,
        rod(d, 3.0),
        0.01,
        0.003,
        &[0.05],
        || false,
    )
    .unwrap()
}

// Independent physical variables: flux, piezo charge, other charge, extension,
// velocity, supplied electrical work, each resistor heat, mechanical work.
// No native matrix, descriptor index, storage callback or propagator is used.
fn derivatives(y: &[f64; 9], d: f64, voltage: f64, force: f64) -> [f64; 9] {
    let theta = 2.5 * d;
    let capacitance = 0.25 * (2.0 - 10.0 * d * d);
    let piezo_voltage = (y[1] - theta * y[3]) / capacitance;
    let other_voltage = y[2] / 0.3;
    let inductor_current = y[0] / 0.5;
    let source_current = (voltage - other_voltage) / 0.8;
    [
        other_voltage - piezo_voltage,
        inductor_current - piezo_voltage / 1.7,
        source_current - inductor_current,
        y[4],
        force - 2.5 * y[3] + theta * piezo_voltage, // original consistent rod mass = 1
        voltage * source_current,
        (voltage - other_voltage) * source_current,
        piezo_voltage * piezo_voltage / 1.7,
        force * y[4],
    ]
}
fn reference(y: &mut [f64; 9], d: f64, voltage: f64, force: f64, h: f64) {
    y[5..].fill(0.0);
    let dt = h / 1024.0;
    for _ in 0..1024 {
        let a = derivatives(y, d, voltage, force);
        let b = derivatives(
            &core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]),
            d,
            voltage,
            force,
        );
        let c = derivatives(
            &core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]),
            d,
            voltage,
            force,
        );
        let e = derivatives(
            &core::array::from_fn(|i| y[i] + dt * c[i]),
            d,
            voltage,
            force,
        );
        for i in 0..9 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + e[i]) / 6.0;
        }
    }
}

#[test]
fn g1_prepared_actual_rlc_retains_all_states_and_independent_heat_and_source_work() {
    for d in [-0.2, 0.0, 0.2] {
        let model = rlc(d);
        let original_rod = model.rod().clone();
        let sources = model.circuit().resistors().to_vec();
        let mut runtime = prepare(model, &[0.05], 0.04);
        let mut expected = [0.03, 0.012, 0.02, 0.01, 0.003, 0.0, 0.0, 0.0, 0.0];
        assert!(runtime.last_step().is_none());
        for index in 0..8 {
            let voltage = if index % 2 == 0 { 0.05 } else { -0.03 };
            let force = 0.015 + 0.001 * index as f64;
            reference(&mut expected, d, voltage, force, 0.04);
            let report = runtime.step(&[voltage], force, solve(), || false).unwrap();
            let model = runtime.model();
            near(model.circuit().fluxes_wb()[0], expected[0], 2e-11);
            near(model.circuit().charges_c()[0], expected[1], 2e-11);
            near(model.circuit().charges_c()[1], expected[2], 2e-11);
            near(model.extension_m(), expected[3], 2e-11);
            near(model.velocity_m_s(), expected[4], 2e-11);
            near(report.electrical_supply_j, expected[5], 2e-12);
            near(report.mechanical_supply_j, expected[8], 2e-12);
            let heat = runtime.resistor_heat().unwrap();
            assert_eq!((heat[0].branch, heat[1].branch), (1, 4));
            near(heat[0].heat_j, expected[6], 2e-12);
            near(heat[1].heat_j, expected[7], 2e-12);
            near(report.heat_j, expected[6] + expected[7], 2e-12);
            near(
                runtime.mean_electrical_outputs().unwrap()[0],
                expected[5] / (voltage * 0.04),
                2e-11,
            );
            near(
                heat[0].mean_current_a,
                expected[5] / (voltage * 0.04),
                2e-11,
            );
            assert!(
                heat.iter()
                    .all(|h| h.heat_j >= 0.0 && h.rms_current_a + 1e-12 >= h.mean_current_a.abs())
            );
            near(report.energy_residual_j, 0.0, 2e-12);
            near(
                report.other_electrical_energy_j[1],
                expected[0] * expected[0] / 1.0 + expected[2] * expected[2] / 0.6,
                2e-12,
            );
            assert_eq!(model.rod(), &original_rod);
            assert_eq!(model.circuit().resistors(), sources.as_slice());
            // Original native order: five differential coordinates, then
            // node potentials 1, 2 and 3, followed by constraint currents.
            let full = runtime.descriptor_state();
            assert_eq!(full.len(), 11);
            near(full[5], voltage, 2e-12);
            near(full[6], expected[2] / 0.3, 2e-11);
            near(
                full[7],
                (expected[1] - 2.5 * d * expected[3]) / (0.25 * (2.0 - 10.0 * d * d)),
                2e-11,
            );
        }
        let last = runtime.model().clone();
        assert_eq!(runtime.into_model(), last);
    }
}

#[test]
fn g1_prepared_real_source_resistor_retains_pure_input_feedthrough_heat() {
    let model = PiezoelectricCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (2, 0, Branch::VoltageSource { port: 0 }),
                (2, 0, Branch::Resistor { ohms: f64::NAN }),
            ],
            transformers: vec![],
        },
        vec![source(2, 2.0)],
        vec![],
        vec![0.0],
        0,
        rod(0.2, 3.0),
        0.0,
        0.0,
        &[3.0],
        || false,
    )
    .unwrap();
    let original = model.clone();
    let mut runtime = prepare(model, &[3.0], 0.25);
    for voltage in [3.0, -3.0] {
        let report = runtime.step(&[voltage], 0.0, solve(), || false).unwrap();
        let heat = &runtime.resistor_heat().unwrap()[0];
        near(heat.heat_j, 1.125, 2e-12);
        near(heat.mean_voltage_v, voltage, 2e-12);
        near(heat.mean_current_a, voltage / 2.0, 2e-12);
        near(heat.rms_current_a, 1.5, 2e-12);
        near(heat.mean_power_w, 4.5, 2e-12);
        near(report.electrical_supply_j, 1.125, 2e-12);
        near(report.heat_j, 1.125, 2e-12);
        near(report.stored_before_j, 0.0, 0.0);
        near(report.stored_after_j, 0.0, 0.0);
        near(report.energy_residual_j, 0.0, 2e-12);
        assert_eq!(runtime.model(), &original);
    }
}

#[test]
fn g1_prepared_current_source_uses_integrated_terminal_voltage_and_rc_loss() {
    let model = PiezoelectricCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::Resistor { ohms: f64::NAN }),
                (1, 0, Branch::CurrentSource { port: 0 }),
            ],
            transformers: vec![],
        },
        vec![source(1, 2.0)],
        vec![],
        vec![0.01],
        0,
        rod(0.0, 3.0),
        0.0,
        0.0,
        &[0.03],
        || false,
    )
    .unwrap();
    let h = 0.2;
    let mut runtime = prepare(model, &[0.03], h);
    let mut charge = 0.01;
    for current in [0.03, -0.02, 0.0] {
        // C = 0.5 F, R = 2 ohm, hence tau = 1 s.
        let equilibrium_v = 2.0 * current;
        let transient_v = charge / 0.5 - equilibrium_v;
        let decay = fs_math::det::exp(-h);
        let integral_v = equilibrium_v * h + transient_v * (1.0 - decay);
        let expected_heat = 0.5
            * (equilibrium_v * equilibrium_v * h
                + 2.0 * equilibrium_v * transient_v * (1.0 - decay)
                + 0.5 * transient_v * transient_v * (1.0 - decay * decay));
        charge = 0.5 * (equilibrium_v + transient_v * decay);
        let report = runtime.step(&[current], 0.0, solve(), || false).unwrap();
        near(runtime.model().charge_c(), charge, 2e-12);
        near(
            runtime.mean_electrical_outputs().unwrap()[0],
            integral_v / h,
            2e-12,
        );
        near(report.electrical_supply_j, current * integral_v, 2e-12);
        near(report.heat_j, expected_heat, 2e-12);
        near(report.energy_residual_j, 0.0, 2e-12);
    }
}

#[test]
fn g3_prepared_late_cancellation_and_material_refusal_publish_nothing_and_retry_exactly() {
    let model = loaded(0.2, 3.0, 0.8, 0.01, 0.003, 0.004);
    let mut reference = prepare(model.clone(), &[], 0.02);
    let mut polls = 0;
    let expected = reference
        .step(&[], 0.015, solve(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 4);
    let mut runtime = prepare(model.clone(), &[], 0.02);
    let original_full = runtime.descriptor_state().to_vec();
    let mut seen = 0;
    assert!(matches!(
        runtime.step(&[], 0.015, solve(), || {
            seen += 1;
            seen == polls
        }),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert_eq!(seen, polls);
    assert_eq!(runtime.model(), &model);
    assert_eq!(runtime.descriptor_state(), original_full.as_slice());
    assert!(runtime.last_step().is_none());
    assert!(runtime.resistor_heat().is_none());
    assert!(runtime.mean_electrical_outputs().is_none());
    let retry = runtime.step(&[], 0.015, solve(), || false).unwrap();
    assert_eq!(retry, expected);
    assert_eq!(runtime.model(), reference.model());
    let accepted_model = runtime.model().clone();
    let accepted_full = runtime.descriptor_state().to_vec();
    let accepted_heat = runtime.resistor_heat().unwrap().to_vec();
    assert!(runtime.step(&[], 1e8, solve(), || false).is_err());
    assert!(runtime.step(&[0.0], 0.0, solve(), || false).is_err());
    assert_eq!(runtime.model(), &accepted_model);
    assert_eq!(runtime.descriptor_state(), accepted_full.as_slice());
    assert_eq!(runtime.resistor_heat().unwrap(), accepted_heat.as_slice());
    assert_eq!(runtime.last_step(), Some(&retry));
}

#[test]
fn g3_prepared_refuses_higher_index_source_constraints_and_exhausted_cold_budget() {
    let model = PiezoelectricCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        vec![],
        vec![],
        vec![0.005],
        0,
        rod(0.2, 3.0),
        0.01,
        0.0,
        &[0.0],
        || false,
    )
    .unwrap();
    let (split, dense) = budgets();
    assert!(matches!(PreparedPiezoelectricCircuit::new(
        model,&[0.0],0.02,split,dense,|| false,
    ),Err(MaterialCircuitError::Invalid(message)) if message.contains("index-one")));
    let model = loaded(0.2, 3.0, 0.8, 0.01, 0.0, 0.004);
    assert!(
        PreparedPiezoelectricCircuit::new(
            model.clone(),
            &[],
            0.02,
            split,
            DenseExponentialBudget {
                max_dimension: 1,
                ..dense
            },
            || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedPiezoelectricCircuit::new(model, &[], 0.02, split, dense, || true,),
        Err(MaterialCircuitError::Cancelled)
    ));
}
