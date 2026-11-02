//! G1/G3: real-area force and torque share one static set and one body response.

use core::f64::consts::PI;
use fs_tribo::adhesive::spinning::contact::*;

const RADIUS: f64 = 0.2;

fn law(static_capacity: f64, kinetic_capacity: f64) -> AdhesiveWrenchContact {
    AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
        contact_radius_m: RADIUS,
        static_shear_strength_pa: static_capacity / (PI * RADIUS * RADIUS),
        kinetic_shear_strength_pa: kinetic_capacity / (PI * RADIUS * RADIUS),
    })
    .unwrap()
}

fn controls() -> AdhesiveWrenchConfig {
    AdhesiveWrenchConfig {
        maximum_static_iterations: 96,
        maximum_iterations: 512,
        maximum_evaluations: 4096,
        force_tolerance_n: 1e-8,
        moment_tolerance_nm: 1e-9,
        displacement_tolerance_m: 1e-12,
        energy_tolerance_j: 1e-9,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 100.0,
        maximum_spin_moment_nm: 100.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tolerance {tolerance:e})"
    );
}

fn product(c: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    c.map(|row| row.into_iter().zip(x).map(|(a, b)| a * b).sum())
}

fn input(free: [f64; 3], c: [[f64; 3]; 3]) -> AdhesiveWrenchInput {
    AdhesiveWrenchInput {
        free_displacement_m: [free[0], free[1]],
        free_spin_angle_rad: free[2] / RADIUS,
        scaled_compliance_m_per_n: c,
    }
}

#[test]
fn pure_spin_uses_actual_disk_first_moment_and_retains_body_recoil() {
    let law = law(15.0, 10.0);
    let c = [[0.01, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.03]];
    let result = law
        .advance(&input([0.0, 0.0, 1.0], c), controls(), || false)
        .unwrap();
    assert_eq!(result.mode, AdhesiveWrenchMode::Sliding);
    let expected_moment = -(2.0 / 3.0) * 10.0 * RADIUS;
    close(result.spin_moment_nm, expected_moment, 1e-14);
    close(
        result.spin_angle_rad,
        (1.0 + c[2][2] * expected_moment / RADIUS) / RADIUS,
        1e-13,
    );
    assert_eq!(result.force_n, [0.0; 2]);
    close(
        result.dissipated_energy_j,
        -expected_moment * result.spin_angle_rad,
        1e-13,
    );
    close(result.work_residual_j, 0.0, 1e-13);
}

#[test]
fn actual_coupled_sliding_obeys_all_cross_axis_responses_and_native_work() {
    let law = law(6.0, 4.0);
    let c = [
        [0.020, 0.003, -0.004],
        [0.003, 0.015, 0.002],
        [-0.004, 0.002, 0.012],
    ];
    let actual = [0.11, -0.07, 0.13];
    let native = law
        .kinetic_patch()
        .evaluate([actual[0], actual[1]], actual[2] / RADIUS)
        .unwrap();
    let expected = [
        native.force_n[0],
        native.force_n[1],
        native.spin_moment_nm / RADIUS,
    ];
    let recoil = product(c, expected);
    let free = core::array::from_fn(|i| actual[i] - recoil[i]);
    let result = law
        .advance_phase(
            &input(free, c),
            AdhesiveWrenchPhase::Kinetic,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result.mode, AdhesiveWrenchMode::Sliding);
    for i in 0..2 {
        close(result.force_n[i], native.force_n[i], 1e-8);
        close(result.slip_displacement_m[i], actual[i], 3e-10);
    }
    close(result.spin_moment_nm, native.spin_moment_nm, 1e-9);
    close(result.spin_angle_rad, actual[2] / RADIUS, 3e-9);
    close(result.contact_work_j, native.dissipated_power_w, 1e-9);
    assert!(result.force_n[0].hypot(result.force_n[1]) < 4.0);
    assert!(result.spin_moment_nm.abs() < (2.0 / 3.0) * 4.0 * RADIUS);
}

#[test]
fn static_force_and_moment_share_the_full_junction_capacity() {
    let law = law(10.0, 8.0);
    let c = [
        [0.02, 0.003, -0.004],
        [0.003, 0.015, 0.002],
        [-0.004, 0.002, 0.012],
    ];
    let requested = [2.0, -1.0, 1.5];
    let result = law
        .advance(&input(product(c, requested), c), controls(), || false)
        .unwrap();
    assert_eq!(result.mode, AdhesiveWrenchMode::Sticking);
    close(result.force_n[0], -requested[0], 1e-13);
    close(result.force_n[1], -requested[1], 1e-13);
    close(result.spin_moment_nm, -RADIUS * requested[2], 1e-13);
    assert_eq!(result.slip_displacement_m, [0.0; 2]);
    assert_eq!(result.spin_angle_rad, 0.0);
    assert_eq!(result.dissipated_energy_j, 0.0);

    // Each marginal limit admits these components separately, but their joint
    // demand lies outside the actual uniform-disk static surface.
    let impossible = [9.5, 0.0, 0.95 * (2.0 / 3.0) * 10.0];
    let mixed = law
        .advance(&input(product(c, impossible), c), controls(), || false)
        .unwrap();
    assert_eq!(mixed.mode, AdhesiveWrenchMode::Sliding);
    assert!(mixed.dissipated_energy_j > 0.0);
}

#[test]
fn moving_phase_cannot_gain_the_larger_static_twisting_capacity() {
    let law = law(20.0, 5.0);
    let c = [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]];
    let trial = input([0.0, 0.0, 0.08], c);
    let at_rest = law.advance(&trial, controls(), || false).unwrap();
    let moving = law
        .advance_phase(&trial, AdhesiveWrenchPhase::Kinetic, controls(), || false)
        .unwrap();
    assert_eq!(at_rest.mode, AdhesiveWrenchMode::Sticking);
    assert_eq!(moving.mode, AdhesiveWrenchMode::Sliding);
    close(moving.spin_moment_nm, -(2.0 / 3.0) * 5.0 * RADIUS, 1e-14);
    assert!(moving.spin_angle_rad > 0.0);
}

#[test]
fn covariance_zero_strength_and_atomic_refusals_preserve_the_model() {
    let law = law(6.0, 4.0);
    let c = [[0.02, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.03]];
    let a = law
        .advance(&input([0.3, -0.1, 0.2], c), controls(), || false)
        .unwrap();
    let b = law
        .advance(&input([0.1, 0.3, 0.2], c), controls(), || false)
        .unwrap();
    close(b.force_n[0], -a.force_n[1], 1e-8);
    close(b.force_n[1], a.force_n[0], 1e-8);
    close(b.spin_moment_nm, a.spin_moment_nm, 1e-9);
    close(b.dissipated_energy_j, a.dissipated_energy_j, 1e-9);
    let trial = input([0.3, -0.1, 0.2], c);
    let before = law;
    assert_eq!(
        law.advance(&trial, controls(), || true),
        Err(AdhesiveWrenchError::Cancelled)
    );
    let mut limited = controls();
    limited.maximum_spin_moment_nm = 0.0;
    assert!(matches!(
        law.advance(&trial, limited, || false),
        Err(AdhesiveWrenchError::EffortBudget { .. })
    ));
    assert_eq!(law, before);
    assert_eq!(law.advance(&trial, controls(), || false).unwrap(), a);
    let free = crate::law(0.0, 0.0)
        .advance(&trial, controls(), || false)
        .unwrap();
    assert_eq!(free.mode, AdhesiveWrenchMode::Frictionless);
    assert_eq!(free.force_n, [0.0; 2]);
    assert_eq!(free.spin_moment_nm, 0.0);
    assert_eq!(free.slip_displacement_m, trial.free_displacement_m);
    close(free.spin_angle_rad, trial.free_spin_angle_rad, 1e-14);
}
