//! Exact constrained-oscillator references with original electrode work.
use super::*;
use fs_couple::piezoelectric::prepared::PreparedPiezoelectricVoltage;

fn prepared_voltage(body: PiezoelectricBody, dt: f64) -> PreparedPiezoelectricVoltage {
    PreparedPiezoelectricVoltage::new(
        body,
        dt,
        LinearSplitBudget {
            max_modes: 2,
            max_omega_dt: 1.0,
        },
        DenseExponentialBudget {
            max_dimension: 4,
            max_dense_bytes: 8_000_000,
            max_multiply_terms: 1_000_000_000,
        },
        || false,
    )
    .unwrap()
}

fn voltage_reference(body: &PiezoelectricBody, voltage: f64, force: f64, t: f64) -> [f64; 2] {
    let rod = body.rod();
    let omega = (rod.short_circuit_stiffness_n_m() / rod.modal_mass_kg()).sqrt();
    let equilibrium = (force + rod.coupling_c_m() * voltage) / rod.short_circuit_stiffness_n_m();
    let offset = body.extension_m() - equilibrium;
    let (sn, cs) = (omega * t).sin_cos();
    [
        equilibrium + offset * cs + body.velocity_m_s() * sn / omega,
        -omega * offset * sn + body.velocity_m_s() * cs,
    ]
}

#[test]
fn g1_exact_short_circuit_retains_the_original_mechanical_frequency_and_state() {
    for (d, density) in [(0.2, 3.0), (-0.2, 3.0), (0.0, 3.0), (0.2, 6.0)] {
        let source = rod(d, density);
        let charge = source.coupling_c_m() * 0.001;
        let original = PiezoelectricBody::new(source, 0.001, 0.0003, charge).unwrap();
        let omega =
            (original.rod().short_circuit_stiffness_n_m() / original.rod().modal_mass_kg()).sqrt();
        let initial_energy = original.stored_energy_j().unwrap();
        let dt = 0.75 / omega;
        let mut runtime = prepared_voltage(original.clone(), dt);
        assert_eq!(runtime.body(), &original);
        assert_eq!(runtime.duration_s(), dt);
        near(runtime.voltage_v(), 0.0, 1e-15);
        for i in 1..=128 {
            let expected = voltage_reference(&original, runtime.voltage_v(), 0.0, i as f64 * dt);
            let step = runtime.step(0.0, 1e-10, controls(), || false).unwrap();
            near(runtime.body().extension_m(), expected[0], 2e-12);
            near(runtime.body().velocity_m_s(), expected[1], 2e-12);
            near(
                runtime.body().voltage_v().unwrap(),
                runtime.voltage_v(),
                2e-14,
            );
            near(
                runtime.body().stored_energy_j().unwrap(),
                initial_energy,
                2e-13,
            );
            near(step.electrical_supply_j, 0.0, 2e-15);
            assert_eq!(step.mechanical_supply_j, 0.0);
            near(step.charge_residual_c, 0.0, 2e-14);
        }
        assert_eq!(runtime.body().rod(), original.rod());
        let accepted = runtime.body().clone();
        assert_eq!(runtime.into_body(), accepted);
    }
}

#[test]
fn g1_biased_clamp_retains_forced_equilibrium_charge_anchor_and_both_source_works() {
    for d in [0.2, -0.2, 0.0] {
        let source = rod(d, 3.0);
        let theta = source.coupling_c_m();
        let c = source.axially_clamped_capacitance_f();
        let initial =
            PiezoelectricBody::new(source, 0.001, -0.0003, c * 0.002 + theta * 0.001).unwrap();
        let force = 0.0002;
        let dt = 0.3;
        let initial_energy = initial.stored_energy_j().unwrap();
        let mut runtime = prepared_voltage(initial.clone(), dt);
        assert_eq!(runtime.body(), &initial);
        let voltage = initial.voltage_v().unwrap();
        assert_eq!(runtime.voltage_v(), voltage);
        let mut previous_x = initial.extension_m();
        for i in 1..=128 {
            let expected = voltage_reference(&initial, voltage, force, i as f64 * dt);
            let step = runtime.step(force, 1e-10, controls(), || false).unwrap();
            near(runtime.body().extension_m(), expected[0], 2e-12);
            near(runtime.body().velocity_m_s(), expected[1], 2e-12);
            assert_eq!(
                runtime.body().charge_c(),
                initial.charge_c() + theta * (runtime.body().extension_m() - initial.extension_m())
            );
            near(
                runtime.body().charge_c(),
                initial.charge_c() + theta * (expected[0] - initial.extension_m()),
                2e-12,
            );
            near(
                step.mechanical_supply_j,
                force * (expected[0] - previous_x),
                2e-15,
            );
            near(
                step.electrical_supply_j,
                voltage * theta * (expected[0] - previous_x),
                2e-15,
            );
            near(runtime.body().voltage_v().unwrap(), voltage, 2e-14);
            assert_eq!(step.mean_voltage_v, voltage);
            near(
                runtime.body().stored_energy_j().unwrap() - initial_energy,
                (force + theta * voltage) * (expected[0] - initial.extension_m()),
                2e-14,
            );
            previous_x = expected[0];
        }
        let equilibrium = (force + theta * voltage) / initial.rod().short_circuit_stiffness_n_m();
        let resting = PiezoelectricBody::new(
            initial.rod().clone(),
            equilibrium,
            0.0,
            c * voltage + theta * equilibrium,
        )
        .unwrap();
        let mut stationary = prepared_voltage(resting.clone(), dt);
        stationary.step(force, 1e-10, controls(), || false).unwrap();
        near(stationary.body().extension_m(), equilibrium, 2e-14);
        near(stationary.body().velocity_m_s(), 0.0, 2e-14);
        near(stationary.body().charge_c(), resting.charge_c(), 2e-14);
    }
}

#[test]
fn g4_fixed_voltage_refuses_a_hidden_interior_strain_excursion() {
    let initial = PiezoelectricBody::new(rod(0.2, 3.0), 0.0, 0.0, 0.0).unwrap();
    let omega =
        (initial.rod().short_circuit_stiffness_n_m() / initial.rod().modal_mass_kg()).sqrt();
    let dt = 2.0 * core::f64::consts::PI / omega;
    near(voltage_reference(&initial, 0.0, 1.0, dt)[0], 0.0, 1e-12);
    let mut runtime = prepared_voltage(initial.clone(), dt);
    assert!(runtime.step(1.0, 1e-10, controls(), || false).is_err());
    assert_eq!(runtime.body(), &initial);
}

#[test]
fn g4_fixed_voltage_cancellation_and_invalid_controls_keep_full_state_retryable() {
    let initial = PiezoelectricBody::new(rod(0.2, 3.0), 0.001, -0.0003, 0.0013).unwrap();
    let mut runtime = prepared_voltage(initial.clone(), 0.1);
    let mut calls = 0;
    assert!(matches!(
        runtime.step(0.0002, 1e-10, controls(), || {
            calls += 1;
            calls == 2
        }),
        Err(PiezoelectricCouplingError::Cancelled)
    ));
    assert_eq!(runtime.body(), &initial);
    assert!(runtime.step(f64::NAN, 1e-10, controls(), || false).is_err());
    assert!(runtime.step(0.0002, 0.0, controls(), || false).is_err());
    assert_eq!(runtime.body(), &initial);
    let step = runtime.step(0.0002, 1e-10, controls(), || false).unwrap();
    let mut clean = prepared_voltage(initial, 0.1);
    let clean_step = clean.step(0.0002, 1e-10, controls(), || false).unwrap();
    assert_eq!(runtime.body(), clean.body());
    assert_eq!(step.electrical_supply_j, clean_step.electrical_supply_j);
    assert_eq!(step.energy_residual_j, clean_step.energy_residual_j);
}
