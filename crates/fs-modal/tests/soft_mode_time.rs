//! G1/G3 near-rigid modal limits: finite motion, stored energy and port work.

use fs_math::{c64::C64, det};
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel, ModalAcousticWorkspace, advance_exact_zoh,
};

fn mode(omega: f64, zeta: f64) -> ModalAcousticMode {
    ModalAcousticMode {
        angular_frequency_rad_s: omega,
        damping_ratio: zeta,
        pressure_per_modal_velocity: C64::ZERO,
    }
}

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn model(omega: f64, zeta: f64, maximum_energy: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new_with_free_coordinates(
        8,
        vec![mode(omega, zeta)],
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 2.0e200,
            maximum_abs_velocity_m_sqrt_kg_per_s: 1.0e16,
            maximum_total_energy_j: maximum_energy,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= 16.0 * f64::EPSILON * scale,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_soft_undamped_modes_recover_free_motion_and_signed_work() {
    let initial = state(0.3, -0.4);
    let force = 1.7;
    let h = 0.125;
    // The first omitted restoring-force contribution is below 1e-280 here.
    // Thus the independent free-coordinate polynomial is the rounded limit.
    let delta = h * (initial.velocity_m_sqrt_kg_per_s + 0.5 * h * force);
    let expected = state(
        initial.displacement_m_sqrt_kg + delta,
        initial.velocity_m_sqrt_kg_per_s + h * force,
    );
    for omega in [0.0, 1.0e-140, 1.0e-155, 1.0e-160, 1.0e-200, 1.0e-300] {
        let mut whole = model(omega, 0.0, 1.0);
        whole.restore_states(&[initial]).unwrap();
        let mut split = whole.clone();
        let mut workspace = ModalAcousticWorkspace::new(&whole);
        let direct = advance_exact_zoh(mode(omega, 0.0), initial, force, h);
        let cached = whole.advance_exact_zoh_step(0, initial, force);
        assert_eq!(
            direct, cached,
            "cached and direct owner disagree at {omega:e}"
        );
        let frame = whole.step_into(&[force], &mut workspace).unwrap();
        close(
            whole.states()[0].displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            1.0,
        );
        close(
            whole.states()[0].velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            1.0,
        );
        close(frame.input_work_j, force * delta, 1.0);
        close(
            frame.total_modal_energy_j,
            0.5 * expected.velocity_m_sqrt_kg_per_s.powi(2),
            1.0,
        );
        assert_eq!(frame.viscous_dissipation_j, 0.0);
        let first = split.step_duration(&[force], 0.5 * h).unwrap();
        let second = split.step_duration(&[force], 0.5 * h).unwrap();
        close(
            split.states()[0].displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            1.0,
        );
        close(
            split.states()[0].velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            1.0,
        );
        close(
            first.input_work_j + second.input_work_j,
            frame.input_work_j,
            1.0,
        );
    }
}

#[test]
fn g1_overflowing_static_equilibrium_keeps_a_finite_attack_and_work() {
    // omega^2 is normal, but F/omega^2 overflows. Neither the finite-time
    // displacement nor the physical work requires that static equilibrium.
    let (omega, force, h) = (1.0e-140, 1.0e30, 1.0e-15);
    let mut actual = model(omega, 0.0, 1.0e32);
    let expected_velocity = h * force;
    let expected_displacement = 0.5 * h * expected_velocity;
    let frame = actual.step_duration(&[force], h).unwrap();
    close(
        actual.states()[0].displacement_m_sqrt_kg,
        expected_displacement,
        1.0,
    );
    close(
        actual.states()[0].velocity_m_sqrt_kg_per_s,
        expected_velocity,
        expected_velocity,
    );
    close(
        frame.input_work_j,
        force * expected_displacement,
        frame.total_modal_energy_j,
    );
    close(
        frame.total_modal_energy_j,
        0.5 * expected_velocity * expected_velocity,
        frame.total_modal_energy_j,
    );
    assert_eq!(frame.viscous_dissipation_j, 0.0);
}

#[test]
fn g1_tiny_stiffness_retains_finite_stored_energy_and_restoring_force() {
    let omega = 1.0e-200;
    let initial = state(1.0e200, 0.0);
    let h = 0.125;
    // omega*q=1: H=1/2 J and q''(0)=-omega, despite omega^2 underflowing.
    // Damping corrections over this interval are below representable accuracy.
    for zeta in [0.0, 0.5, 1.0, 2.0] {
        let mut actual = model(omega, zeta, 1.0);
        actual.restore_states(&[initial]).unwrap();
        let frame = actual.step(&[0.0]).unwrap();
        assert_eq!(
            actual.states()[0].displacement_m_sqrt_kg,
            initial.displacement_m_sqrt_kg
        );
        close(
            actual.states()[0].velocity_m_sqrt_kg_per_s,
            -h * omega,
            h * omega,
        );
        close(frame.total_modal_energy_j, 0.5, 1.0);
        assert_eq!(frame.input_work_j, 0.0);
        assert_eq!(frame.viscous_dissipation_j, 0.0); // Physical heat is below f64 range.
    }
}

#[test]
fn g0_soft_mode_energy_budget_cannot_be_bypassed_by_stiffness_underflow() {
    let mut actual = model(1.0e-200, 0.0, 0.25);
    let before = actual.states().to_vec();
    let error = actual.restore_states(&[state(1.0e200, 0.0)]).unwrap_err();
    assert!(matches!(
        error,
        ModalAcousticTimeError::BudgetExceeded {
            what: "total modal energy",
            ..
        }
    ));
    assert_eq!(actual.states(), before);
}

#[test]
fn g3_ordinary_undamped_transitions_keep_the_existing_operation_tree() {
    let initial = state(0.0002, -0.3);
    let force = 1.7;
    for omega in [0.1, 7.0, core::f64::consts::TAU * 440.0] {
        for h in [1.0 / 48_000.0, 2.0e-4, 0.02] {
            let angle = omega * h;
            let sine = det::sin(angle);
            let half_sine = det::sin(0.5 * angle);
            let one_minus_cosine = 2.0 * half_sine * half_sine;
            let q = initial.displacement_m_sqrt_kg;
            let v = initial.velocity_m_sqrt_kg_per_s;
            let expected = state(
                q + (v / omega) * sine + (force / (omega * omega) - q) * one_minus_cosine,
                v + ((force / omega - omega * q) * sine - v * one_minus_cosine),
            );
            // Preserve the original increment's parentheses as well as its terms.
            let expected_q =
                q + ((v / omega) * sine + (force / (omega * omega) - q) * one_minus_cosine);
            let actual = advance_exact_zoh(mode(omega, 0.0), initial, force, h);
            assert_eq!(
                actual.displacement_m_sqrt_kg.to_bits(),
                expected_q.to_bits()
            );
            assert_eq!(
                actual.velocity_m_sqrt_kg_per_s.to_bits(),
                expected.velocity_m_sqrt_kg_per_s.to_bits()
            );
        }
    }
}
