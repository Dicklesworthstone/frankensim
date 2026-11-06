//! G0/G1/G4 source-preserving exact motion, charge and source-work checks.
#[path = "prepared/voltage.rs"]
mod voltage;

use super::{near, rod};
use fs_couple::piezoelectric::prepared::{
    ExactPiezoelectricControls, PreparedPiezoelectricCurrent,
};
use fs_couple::piezoelectric::{PiezoelectricBody, PiezoelectricCouplingError};
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};

fn controls() -> ExactPiezoelectricControls {
    ExactPiezoelectricControls {
        energy_tolerance_j: 1e-12,
        charge_tolerance_c: 1e-12,
    }
}

fn prepared(body: PiezoelectricBody, dt: f64) -> PreparedPiezoelectricCurrent {
    PreparedPiezoelectricCurrent::new(
        body,
        dt,
        LinearSplitBudget {
            max_modes: 3,
            max_omega_dt: 1.0,
        },
        DenseExponentialBudget {
            max_dimension: 32,
            max_dense_bytes: 8_000_000,
            max_multiply_terms: 1_000_000_000,
        },
        || false,
    )
    .unwrap()
}

// Independent forced-oscillator solution with Q(t)=Q0+I*t. The mean voltage
// integrates the continuous motion, without endpoint trapezoids.
fn reference(body: &PiezoelectricBody, force: f64, current: f64, t: f64) -> [f64; 4] {
    let rod = body.rod();
    let m = rod.modal_mass_kg();
    let k = rod.open_circuit_stiffness_n_m();
    let c = rod.axially_clamped_capacitance_f();
    let theta = rod.coupling_c_m();
    let omega = (k / m).sqrt();
    let angle = omega * t;
    let (sn, cs) = angle.sin_cos();
    let equilibrium = (force + theta * body.charge_c() / c) / k;
    let ramp = theta * current / c / k;
    let offset = body.extension_m() - equilibrium;
    let x = equilibrium + offset * cs
        + body.velocity_m_s() * sn / omega
        + ramp * (t - sn / omega);
    let velocity = -omega * offset * sn
        + body.velocity_m_s() * cs
        + ramp * (1.0 - cs);
    let charge = body.charge_c() + current * t;
    let x_integral = equilibrium * t + offset * sn / omega
        + body.velocity_m_s() * (1.0 - cs) / (omega * omega)
        + ramp * (0.5 * t * t + (cs - 1.0) / (omega * omega));
    let q_integral = body.charge_c() * t + 0.5 * current * t * t;
    [x, velocity, charge, (q_integral - theta * x_integral) / c / t]
}

#[test]
fn g1_exact_current_flow_matches_independent_motion_and_voltage_integral() {
    for (d, density) in [(0.2, 3.0), (-0.2, 3.0), (0.0, 3.0), (0.2, 6.0)] {
        let body = PiezoelectricBody::new(rod(d, density), 0.001, -0.0003, 0.0001).unwrap();
        let original = body.clone();
        let (force, current, dt) = (0.0002, 0.00003, 0.3);
        let expected = reference(&body, force, current, dt);
        let before = body.stored_energy_j().unwrap();
        let mut runtime = prepared(body, dt);
        let step = runtime.step(force, current, controls(), || false).unwrap();
        near(runtime.body().extension_m(), expected[0], 2e-13);
        near(runtime.body().velocity_m_s(), expected[1], 2e-13);
        near(runtime.body().charge_c(), expected[2], 2e-14);
        near(step.mean_voltage_v, expected[3], 2e-13);
        near(
            step.mechanical_supply_j,
            force * (expected[0] - original.extension_m()),
            2e-15,
        );
        near(step.electrical_supply_j, current * dt * expected[3], 2e-15);
        near(
            runtime.body().stored_energy_j().unwrap() - before,
            step.mechanical_supply_j + step.electrical_supply_j,
            2e-15,
        );
        near(step.charge_residual_c, 0.0, 2e-14);
        assert_eq!(runtime.body().rod().material(), original.rod().material());
    }
}

#[test]
fn g1_repeated_prepared_open_circuit_preserves_frequency_and_initial_charge() {
    let initial = PiezoelectricBody::new(rod(0.2, 3.0), 0.001, 0.0003, 0.0001).unwrap();
    let omega = (initial.rod().open_circuit_stiffness_n_m()
        / initial.rod().modal_mass_kg()).sqrt();
    // This interval deliberately resolves the original continuous linear flow
    // rather than relying on a small midpoint phase-error approximation.
    let dt = 0.75 / omega;
    let initial_energy = initial.stored_energy_j().unwrap();
    let mut runtime = prepared(initial.clone(), dt);
    for i in 1..=256 {
        let step = runtime.step(0.0, 0.0, controls(), || false).unwrap();
        let expected = reference(&initial, 0.0, 0.0, i as f64 * dt);
        near(runtime.body().extension_m(), expected[0], 3e-11);
        near(runtime.body().velocity_m_s(), expected[1], 3e-11);
        near(runtime.body().charge_c(), initial.charge_c(), 2e-13);
        near(runtime.body().stored_energy_j().unwrap(), initial_energy, 2e-13);
        assert_eq!(step.mechanical_supply_j, 0.0);
        assert_eq!(step.electrical_supply_j, 0.0);
    }
}

#[test]
fn g4_exact_flow_refuses_hidden_interior_excursion_and_preserves_state() {
    let initial = PiezoelectricBody::new(rod(0.2, 3.0), 0.0, 0.0, 0.0).unwrap();
    let omega = (initial.rod().open_circuit_stiffness_n_m()
        / initial.rod().modal_mass_kg()).sqrt();
    let dt = 2.0 * std::f64::consts::PI / omega;
    // The ideal solution returns to zero at one period, but its intermediate
    // strain exceeds the material envelope. Endpoint-only admission misses it.
    let endpoint = reference(&initial, 1.0, 0.0, dt);
    near(endpoint[0], 0.0, 1e-12);
    let mut runtime = prepared(initial.clone(), dt);
    assert!(runtime.step(1.0, 0.0, controls(), || false).is_err());
    assert_eq!(runtime.body(), &initial);
}

#[test]
fn g4_prepared_current_cancellation_is_atomic_and_retry_is_deterministic() {
    let initial = PiezoelectricBody::new(rod(0.2, 3.0), 0.001, 0.0, 0.0001).unwrap();
    let mut runtime = prepared(initial.clone(), 0.1);
    let mut calls = 0;
    let result = runtime.step(0.0002, 0.00003, controls(), || {
        calls += 1;
        calls == 2
    });
    assert!(matches!(result, Err(PiezoelectricCouplingError::Cancelled)));
    assert_eq!(runtime.body(), &initial);
    runtime.step(0.0002, 0.00003, controls(), || false).unwrap();
    let mut clean = prepared(initial, 0.1);
    clean.step(0.0002, 0.00003, controls(), || false).unwrap();
    assert_eq!(runtime.body(), clean.body());
    let accepted = runtime.body().clone();
    assert!(runtime.step(f64::NAN, 0.0, controls(), || false).is_err());
    assert_eq!(runtime.body(), &accepted);
}
