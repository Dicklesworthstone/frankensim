//! Exact finite arrest of the original uniform disk, with independent work.

use core::f64::consts::PI;
use fs_tribo::adhesive::spinning::contact::{
    AdhesiveWrenchConfig, AdhesiveWrenchContact, AdhesiveWrenchParameters,
};
use fs_tribo::adhesive::spinning::relaxation::{AdhesiveRelaxationError, AdhesiveRelaxationInput};

fn config() -> AdhesiveWrenchConfig {
    AdhesiveWrenchConfig {
        maximum_static_iterations: 64,
        maximum_iterations: 128,
        maximum_evaluations: 256,
        force_tolerance_n: 1e-12,
        moment_tolerance_nm: 1e-12,
        displacement_tolerance_m: 1e-12,
        energy_tolerance_j: 1e-12,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 100.0,
        maximum_spin_moment_nm: 100.0,
    }
}
fn law(radius: f64, strength: f64) -> AdhesiveWrenchContact {
    AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
        contact_radius_m: radius,
        static_shear_strength_pa: 7.0 * strength,
        kinetic_shear_strength_pa: strength,
    })
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        // Independent stopping-time arithmetic can place a nominal boundary
        // a few ulps to either side; the law must not clip to the test's time.
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-2),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}
fn matrix() -> [[f64; 3]; 3] {
    // All three coordinate rows participate. The true uniform-disk rim wrench
    // is proportional to [3,0,1]; this maps it to [6.4,0,6.4].
    [[2.0, 0.04, 0.4], [0.04, 3.0, -0.12], [0.4, -0.12, 5.2]]
}

#[test]
fn pure_translation_and_spin_cross_the_actual_stop_without_extra_work() {
    let radius = 0.07;
    let strength = 5.0;
    let law = law(radius, strength);
    let capacity = strength * PI * radius * radius;
    for spin in [false, true] {
        let rate: f64 = if spin { -3.4 * radius } else { 0.3 };
        let response = if spin { 11.0 } else { 2.0 };
        let resistance = if spin { 2.0 * capacity / 3.0 } else { capacity };
        let stop = rate.abs() / (response * resistance);
        for fraction in [0.25_f64, 1.0, 4.0] {
            let active = fraction.min(1.0);
            let result = law
                .relax_velocity(
                    &AdhesiveRelaxationInput {
                        velocity_m_s: if spin { [0.0; 2] } else { [rate, 0.0] },
                        spin_rate_rad_s: if spin { rate / radius } else { 0.0 },
                        scaled_mobility_per_kg: [
                            [2.0, 0.0, 0.0],
                            [0.0, 3.0, 0.0],
                            [0.0, 0.0, 11.0],
                        ],
                        dt_s: fraction * stop,
                    },
                    config(),
                    || false,
                )
                .unwrap();
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
                close(result.force_impulse_ns[0], -capacity * active * stop);
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
fn combined_rim_motion_uses_the_original_surface_and_every_mobility_entry() {
    let radius = 0.03;
    let strength = 17.0;
    let law = law(radius, strength);
    let capacity = strength * PI * radius * radius;
    let f = capacity * 8.0 / (3.0 * PI);
    let m_over_a = capacity * 8.0 / (9.0 * PI);
    let rate = 0.7;
    let acceleration = 2.0 * f + 0.4 * m_over_a;
    let stop = rate / acceleration;
    let result = law
        .relax_velocity(
            &AdhesiveRelaxationInput {
                velocity_m_s: [rate, 0.0],
                spin_rate_rad_s: rate / radius,
                scaled_mobility_per_kg: matrix(),
                dt_s: 3.0 * stop,
            },
            config(),
            || false,
        )
        .unwrap();
    close(result.initial_wrench.force_n[0], -f);
    close(result.initial_wrench.spin_moment_nm, -radius * m_over_a);
    close(result.force_impulse_ns[0], -f * stop);
    close(result.spin_angular_impulse_nms, -radius * m_over_a * stop);
    close(result.slip_increment_m[0], 0.5 * rate * stop);
    close(result.spin_increment_rad, 0.5 * rate * stop / radius);
    close(result.kinetic_phase_end_s.unwrap(), stop);
    close(
        result.dissipated_energy_j,
        0.5 * (f + m_over_a) * rate * stop,
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
fn time_partition_and_stationary_remainder_preserve_impulse_travel_and_loss() {
    let radius = 0.1;
    let law = law(radius, 4.0);
    let mut input = AdhesiveRelaxationInput {
        velocity_m_s: [0.2, 0.0],
        spin_rate_rad_s: 2.0,
        scaled_mobility_per_kg: matrix(),
        dt_s: 100.0,
    };
    let whole = law.relax_velocity(&input, config(), || false).unwrap();
    let stop = whole.kinetic_phase_end_s.unwrap();
    input.dt_s = 0.4 * stop;
    let first = law.relax_velocity(&input, config(), || false).unwrap();
    assert!(first.kinetic_phase_end_s.is_none());
    input.velocity_m_s = first.velocity_after_m_s;
    input.spin_rate_rad_s = first.spin_rate_after_rad_s;
    input.dt_s = 100.0 - input.dt_s;
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
fn tangent_rotation_and_full_motion_reversal_preserve_the_event_and_heat() {
    let radius = 0.04;
    let law = law(radius, 6.0);
    let input = AdhesiveRelaxationInput {
        velocity_m_s: [0.3, 0.0],
        spin_rate_rad_s: 0.3 / radius,
        scaled_mobility_per_kg: matrix(),
        dt_s: 30.0,
    };
    let original = law.relax_velocity(&input, config(), || false).unwrap();
    let index = [1, 0, 2];
    let sign = [-1.0, 1.0, 1.0];
    let rotated = AdhesiveRelaxationInput {
        velocity_m_s: [-input.velocity_m_s[1], input.velocity_m_s[0]],
        scaled_mobility_per_kg: core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                sign[i] * sign[j] * input.scaled_mobility_per_kg[index[i]][index[j]]
            })
        }),
        ..input
    };
    let rotated = law.relax_velocity(&rotated, config(), || false).unwrap();
    close(rotated.force_impulse_ns[0], -original.force_impulse_ns[1]);
    close(rotated.force_impulse_ns[1], original.force_impulse_ns[0]);
    close(
        rotated.spin_angular_impulse_nms,
        original.spin_angular_impulse_nms,
    );
    close(rotated.dissipated_energy_j, original.dissipated_energy_j);
    assert_eq!(rotated.kinetic_phase_end_s, original.kinetic_phase_end_s);
    let reversed = law
        .relax_velocity(
            &AdhesiveRelaxationInput {
                velocity_m_s: input.velocity_m_s.map(|v| -v),
                spin_rate_rad_s: -input.spin_rate_rad_s,
                ..input
            },
            config(),
            || false,
        )
        .unwrap();
    close(reversed.force_impulse_ns[0], -original.force_impulse_ns[0]);
    close(
        reversed.spin_angular_impulse_nms,
        -original.spin_angular_impulse_nms,
    );
    close(reversed.spin_increment_rad, -original.spin_increment_rad);
    close(reversed.dissipated_energy_j, original.dissipated_energy_j);
    assert_eq!(reversed.kinetic_phase_end_s, original.kinetic_phase_end_s);
}

#[test]
fn true_curvature_is_refused_even_for_a_tiny_interval_and_roundoff_keeps_work_gates() {
    let law = law(0.1, 4.0);
    let mut input = AdhesiveRelaxationInput {
        velocity_m_s: [1.0, 0.0],
        spin_rate_rad_s: 10.0,
        scaled_mobility_per_kg: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dt_s: 1e-100,
    };
    assert!(matches!(
        law.relax_velocity(&input, config(), || false),
        Err(AdhesiveRelaxationError::CurvedTrajectory { .. })
    ));
    input.spin_rate_rad_s = 0.0;
    input.scaled_mobility_per_kg[0][1] = 1e-14;
    input.scaled_mobility_per_kg[1][0] = 1e-14;
    input.dt_s = 100.0;
    let strict = AdhesiveWrenchConfig {
        displacement_tolerance_m: 0.0,
        ..config()
    };
    assert!(matches!(
        law.relax_velocity(&input, strict, || false),
        Err(AdhesiveRelaxationError::Unresolved { .. })
    ));
}

#[test]
fn zero_kinetic_strength_and_tiny_motion_do_not_select_static_capture() {
    let frictionless = AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
        contact_radius_m: 0.1,
        static_shear_strength_pa: 1000.0,
        kinetic_shear_strength_pa: 0.0,
    })
    .unwrap();
    let input = AdhesiveRelaxationInput {
        velocity_m_s: [1e-15, -2e-15],
        spin_rate_rad_s: 3e-15,
        scaled_mobility_per_kg: matrix(),
        dt_s: 4.0,
    };
    let free = frictionless
        .relax_velocity(&input, config(), || false)
        .unwrap();
    assert_eq!(free.velocity_after_m_s, input.velocity_m_s);
    assert_eq!(free.spin_rate_after_rad_s, input.spin_rate_rad_s);
    assert_eq!(free.kinetic_phase_end_s, None);
    assert_eq!(free.stationary_duration_s, 0.0);
    assert_eq!(
        free.slip_increment_m,
        input.velocity_m_s.map(|v| v * input.dt_s)
    );
    assert_eq!(free.dissipated_energy_j, 0.0);
    let moving = AdhesiveRelaxationInput {
        velocity_m_s: [1e-15, 0.0],
        spin_rate_rad_s: 0.0,
        scaled_mobility_per_kg: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ..input
    };
    let stopped = law(0.1, 1.0)
        .relax_velocity(&moving, config(), || false)
        .unwrap();
    assert!(stopped.kinetic_phase_end_s.unwrap() > 0.0);
    assert!(stopped.force_impulse_ns[0] < 0.0);
    assert!(stopped.dissipated_energy_j > 0.0);
}

#[test]
fn cancellation_effort_and_original_matrix_refusals_leave_retry_identical() {
    let law = law(0.1, 4.0);
    let input = AdhesiveRelaxationInput {
        velocity_m_s: [0.2, 0.0],
        spin_rate_rad_s: 2.0,
        scaled_mobility_per_kg: matrix(),
        dt_s: 10.0,
    };
    let expected = law.relax_velocity(&input, config(), || false).unwrap();
    for at in 0..3 {
        let mut calls = 0;
        assert_eq!(
            law.relax_velocity(&input, config(), || {
                let refuse = calls == at;
                calls += 1;
                refuse
            }),
            Err(AdhesiveRelaxationError::Cancelled)
        );
        assert_eq!(
            law.relax_velocity(&input, config(), || false).unwrap(),
            expected
        );
    }
    assert_eq!(
        law.relax_velocity(
            &input,
            AdhesiveWrenchConfig {
                maximum_force_n: 0.0,
                ..config()
            },
            || false
        ),
        Err(AdhesiveRelaxationError::EffortBudget)
    );
    let mut invalid = input;
    invalid.scaled_mobility_per_kg[0][1] += 0.01;
    assert!(matches!(
        law.relax_velocity(&invalid, config(), || false),
        Err(AdhesiveRelaxationError::Contact(_))
    ));
    invalid.scaled_mobility_per_kg = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(matches!(
        law.relax_velocity(&invalid, config(), || false),
        Err(AdhesiveRelaxationError::Contact(_))
    ));
}
