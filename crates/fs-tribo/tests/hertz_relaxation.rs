//! Finite Hertz arrest from independent pressure moments and native work.

use core::f64::consts::PI;
use fs_tribo::sliding_patch::contact::relaxation::{HertzRelaxationError, HertzRelaxationInput};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};

fn config() -> HertzWrenchConfig {
    HertzWrenchConfig {
        maximum_static_iterations: 1,
        maximum_iterations: 1,
        maximum_evaluations: 1,
        force_tolerance_n: 1e-12,
        moment_tolerance_nm: 1e-12,
        displacement_tolerance_m: 1e-12,
        energy_tolerance_j: 1e-12,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 100.0,
        maximum_spin_moment_nm: 100.0,
    }
}
fn law(radius: f64, static_mu: f64, kinetic_mu: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 10.0,
        contact_radius_m: radius,
        static_mu,
        kinetic_mu,
    })
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-2),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}
fn coupled_mobility() -> [[f64; 3]; 3] {
    // At |u|=a*|omega|, the exact Hertz wrench is proportional to [4,0,1].
    // Every matrix row participates and W*[4,0,1] = [8.4,0,8.4].
    [[2.0, 0.04, 0.4], [0.04, 3.0, -0.16], [0.4, -0.16, 6.8]]
}

#[test]
fn translation_and_spin_stop_at_kinetic_time_and_dwell_without_extra_heat() {
    let radius = 0.07;
    let law = law(radius, 0.9, 0.4);
    let equal_bounds = self::law(radius, 0.4, 0.4);
    for spin in [false, true] {
        let rate: f64 = if spin { -3.4 * radius } else { 0.3 };
        let response = if spin { 11.0 } else { 2.0 };
        // Integrating r against Hertz pressure gives 3*pi*a/16, not the
        // uniform-disk first moment 2*a/3 and not an RMS patch radius.
        let resistance = if spin { 4.0 * 3.0 * PI / 16.0 } else { 4.0 };
        let stop = rate.abs() / (response * resistance);
        for fraction in [0.25_f64, 1.0, 4.0] {
            let active = fraction.min(1.0);
            let input = HertzRelaxationInput {
                velocity_m_s: if spin { [0.0; 2] } else { [rate, 0.0] },
                spin_rate_rad_s: if spin { rate / radius } else { 0.0 },
                scaled_mobility_per_kg: [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 11.0]],
                dt_s: fraction * stop,
            };
            let result = law.relax_velocity(&input, config(), || false).unwrap();
            // A stronger static coefficient cannot alter a moving phase.
            assert_eq!(
                result,
                equal_bounds
                    .relax_velocity(&input, config(), || false)
                    .unwrap()
            );
            let heat = rate * rate / response * (active - 0.5 * active * active);
            close(result.dissipated_energy_j, heat);
            close(result.contact_work_j, heat);
            close(result.work_residual_j, 0.0);
            close(result.moving_duration_s, active * stop);
            close(result.stationary_duration_s, (fraction - active) * stop);
            if spin {
                close(result.spin_rate_after_rad_s, rate / radius * (1.0 - active));
                close(
                    result.spin_angular_impulse_nms,
                    -rate.signum() * resistance * radius * active * stop,
                );
                close(
                    result.spin_increment_rad,
                    rate / radius * stop * (active - 0.5 * active * active),
                );
                assert_eq!(result.force_impulse_ns, [0.0; 2]);
            } else {
                close(result.velocity_after_m_s[0], rate * (1.0 - active));
                close(result.force_impulse_ns[0], -resistance * active * stop);
                close(
                    result.slip_increment_m[0],
                    rate * stop * (active - 0.5 * active * active),
                );
                assert_eq!(result.spin_angular_impulse_nms, 0.0);
            }
            if fraction > 1.0 {
                close(result.kinetic_phase_end_s.unwrap(), stop);
                assert_eq!(result.velocity_after_m_s, [0.0; 2]);
                assert_eq!(result.spin_rate_after_rad_s, 0.0);
            }
        }
    }
}

#[test]
fn combined_hertz_rim_motion_retains_the_full_off_diagonal_response() {
    let radius = 0.03;
    let law = law(radius, 0.8, 0.4);
    // Independent closed Hertz integrals at epsilon=1.
    let force = 4.0 * 9.0 * PI / 32.0;
    let moment_over_radius = 4.0 * 9.0 * PI / 128.0;
    let rate = 0.7;
    let acceleration = 2.0 * force + 0.4 * moment_over_radius;
    let stop = rate / acceleration;
    let result = law
        .relax_velocity(
            &HertzRelaxationInput {
                velocity_m_s: [rate, 0.0],
                spin_rate_rad_s: rate / radius,
                scaled_mobility_per_kg: coupled_mobility(),
                dt_s: 3.0 * stop,
            },
            config(),
            || false,
        )
        .unwrap();
    close(result.initial_wrench.force_n[0], -force);
    close(
        result.initial_wrench.spin_moment_nm,
        -radius * moment_over_radius,
    );
    close(result.force_impulse_ns[0], -force * stop);
    close(
        result.spin_angular_impulse_nms,
        -radius * moment_over_radius * stop,
    );
    close(result.slip_increment_m[0], 0.5 * rate * stop);
    close(result.spin_increment_rad, 0.5 * rate * stop / radius);
    close(result.kinetic_phase_end_s.unwrap(), stop);
    close(
        result.dissipated_energy_j,
        0.5 * (force + moment_over_radius) * rate * stop,
    );
    close(result.contact_work_j, result.dissipated_energy_j);
    assert!(result.ray_residual_relative < 2e-15);
    assert!(
        result
            .scaled_velocity_residual_m_s
            .iter()
            .all(|v| v.abs() < 1e-15)
    );
}

#[test]
fn splitting_before_arrest_preserves_impulse_travel_and_heat() {
    let radius = 0.1;
    let law = law(radius, 0.8, 0.4);
    let mut input = HertzRelaxationInput {
        velocity_m_s: [0.2, 0.0],
        spin_rate_rad_s: 2.0,
        scaled_mobility_per_kg: coupled_mobility(),
        dt_s: 10.0,
    };
    let whole = law.relax_velocity(&input, config(), || false).unwrap();
    input.dt_s = 0.4 * whole.kinetic_phase_end_s.unwrap();
    let first = law.relax_velocity(&input, config(), || false).unwrap();
    assert!(first.kinetic_phase_end_s.is_none());
    input.velocity_m_s = first.velocity_after_m_s;
    input.spin_rate_rad_s = first.spin_rate_after_rad_s;
    input.dt_s = 10.0 - input.dt_s;
    let second = law.relax_velocity(&input, config(), || false).unwrap();
    close(
        first.force_impulse_ns[0] + second.force_impulse_ns[0],
        whole.force_impulse_ns[0],
    );
    close(
        first.spin_angular_impulse_nms + second.spin_angular_impulse_nms,
        whole.spin_angular_impulse_nms,
    );
    close(
        first.slip_increment_m[0] + second.slip_increment_m[0],
        whole.slip_increment_m[0],
    );
    close(
        first.spin_increment_rad + second.spin_increment_rad,
        whole.spin_increment_rad,
    );
    close(
        first.dissipated_energy_j + second.dissipated_energy_j,
        whole.dissipated_energy_j,
    );
    input.velocity_m_s = second.velocity_after_m_s;
    input.spin_rate_rad_s = second.spin_rate_after_rad_s;
    let rest = law.relax_velocity(&input, config(), || false).unwrap();
    assert_eq!(rest.force_impulse_ns, [0.0; 2]);
    assert_eq!(rest.spin_angular_impulse_nms, 0.0);
    assert_eq!(rest.slip_increment_m, [0.0; 2]);
    assert_eq!(rest.dissipated_energy_j, 0.0);
    assert_eq!(rest.kinetic_phase_end_s, Some(0.0));
    assert_eq!(rest.stationary_duration_s, input.dt_s);
}

#[test]
fn zero_kinetic_capacity_preserves_motion_and_rest_generates_no_heat() {
    let input = HertzRelaxationInput {
        velocity_m_s: [0.2, -0.3],
        spin_rate_rad_s: 2.0,
        scaled_mobility_per_kg: coupled_mobility(),
        dt_s: 4.0,
    };
    let result = law(0.1, 0.8, 0.0)
        .relax_velocity(&input, config(), || false)
        .unwrap();
    assert_eq!(result.velocity_after_m_s, input.velocity_m_s);
    assert_eq!(result.spin_rate_after_rad_s, input.spin_rate_rad_s);
    assert_eq!(result.force_impulse_ns, [0.0; 2]);
    assert_eq!(result.spin_angular_impulse_nms, 0.0);
    assert_eq!(result.dissipated_energy_j, 0.0);
    assert_eq!(result.kinetic_phase_end_s, None);
    assert_eq!(result.slip_increment_m, [0.8, -1.2]);
    assert_eq!(result.spin_increment_rad, 8.0);
    let rest = law(0.1, 0.8, 0.4)
        .relax_velocity(
            &HertzRelaxationInput {
                velocity_m_s: [0.0; 2],
                spin_rate_rad_s: 0.0,
                ..input
            },
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(rest.kinetic_phase_end_s, Some(0.0));
    assert_eq!(rest.stationary_duration_s, 4.0);
    assert_eq!(rest.dissipated_energy_j, 0.0);
}

#[test]
fn curved_motion_and_cancellation_do_not_publish_an_approximate_event() {
    let law = law(0.1, 0.8, 0.4);
    let input = HertzRelaxationInput {
        velocity_m_s: [0.2, 0.0],
        spin_rate_rad_s: 2.0,
        scaled_mobility_per_kg: coupled_mobility(),
        dt_s: 3.0,
    };
    for dt_s in [3.0, 1e-30] {
        let mut curved = HertzRelaxationInput { dt_s, ..input };
        curved.scaled_mobility_per_kg[2][2] = 7.0;
        assert!(matches!(
            law.relax_velocity(&curved, config(), || false),
            Err(HertzRelaxationError::CurvedTrajectory { .. })
        ));
    }
    let expected = law.relax_velocity(&input, config(), || false).unwrap();
    for cancel_at in 1..=3 {
        let mut calls = 0;
        assert!(matches!(
            law.relax_velocity(&input, config(), || {
                calls += 1;
                calls == cancel_at
            }),
            Err(HertzRelaxationError::Cancelled)
        ));
        assert_eq!(
            law.relax_velocity(&input, config(), || false).unwrap(),
            expected
        );
    }
    let mut limited = config();
    limited.maximum_spin_moment_nm = 0.0;
    assert!(matches!(
        law.relax_velocity(&input, limited, || false),
        Err(HertzRelaxationError::EffortBudget)
    ));
}
