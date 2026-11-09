//! Physical gait kinematics, native wake coupling and bounded screen execution.

use fs_bem::BemError;
use fs_bem::wake2d::{PanelWakeAdaptiveOptions, PanelWakeTransportOptions, WakeSim};
use fs_ornith::OrnithCandidate;
use fs_ornith::screen::{GaitScreenControls, GaitTransport, SectionGait, screen_section_gait};

fn candidate() -> OrnithCandidate {
    OrnithCandidate::from_genes(&[0.5, 0.4, 0.5, 0.5, 0.5])
}

fn gait() -> SectionGait {
    SectionGait {
        heave_amplitude_chord: 0.025,
        pitch_amplitude_rad: 0.035,
        reduced_frequency: 1.4,
        heave_phase_rad: 0.2,
        pitch_phase_rad: -0.4,
        pivot_chord: 0.25,
    }
}

fn controls(substeps: usize) -> GaitScreenControls {
    GaitScreenControls {
        panels: 24,
        intervals: 10,
        interval_nondimensional: 0.025,
        core_radius_chord: 0.05,
        transport: GaitTransport::Fixed(PanelWakeTransportOptions {
            substeps,
            max_coupled_solves: 2 * substeps + 1,
            max_displacement_fraction: None,
        }),
        maximum_coupled_solves: 20000,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

fn position_error(a: &WakeSim, b: &WakeSim) -> f64 {
    let a = a.wake_in_world().unwrap();
    let b = b.wake_in_world().unwrap();
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(a, b)| (a.pos[0] - b.pos[0]).hypot(a.pos[1] - b.pos[1]))
        .fold(0.0, f64::max)
}

#[test]
fn physical_gait_has_the_declared_pivot_phase_reduced_time_and_material_velocity() {
    let gait = gait();
    for time in [0.0, 0.13, 0.7, 2.4] {
        let frame = gait.frame_at(0.09, time).unwrap();
        let phase = 2.0 * gait.reduced_frequency * time;
        let heave = gait.heave_amplitude_chord * (phase + gait.heave_phase_rad).sin();
        let heave_rate = gait.heave_amplitude_chord
            * 2.0
            * gait.reduced_frequency
            * (phase + gait.heave_phase_rad).cos();
        let (s, c) = frame.angle_rad.sin_cos();
        let pivot = [gait.pivot_chord * c, gait.pivot_chord * s];
        close(frame.origin_world[0] + pivot[0], gait.pivot_chord, 2e-15);
        close(frame.origin_world[1] + pivot[1], heave, 2e-15);
        close(
            frame.origin_velocity_world[0] - frame.angular_velocity * pivot[1],
            0.0,
            2e-15,
        );
        close(
            frame.origin_velocity_world[1] + frame.angular_velocity * pivot[0],
            heave_rate,
            2e-15,
        );
        close(
            frame.angle_rad,
            -0.09 + gait.pitch_amplitude_rad * (phase + gait.pitch_phase_rad).sin(),
            2e-15,
        );
        // Independent trajectory differences test the actual LE derivative,
        // including rotation about a nonzero material pivot.
        let h = 1e-5;
        let before = gait.frame_at(0.09, time - h).unwrap();
        let after = gait.frame_at(0.09, time + h).unwrap();
        for axis in 0..2 {
            close(
                (after.origin_world[axis] - before.origin_world[axis]) / (2.0 * h),
                frame.origin_velocity_world[axis],
                2e-11,
            );
        }
        close(
            (after.angle_rad - before.angle_rad) / (2.0 * h),
            frame.angular_velocity,
            2e-11,
        );
    }
    let held = SectionGait {
        reduced_frequency: 0.0,
        ..gait
    };
    assert_eq!(
        held.frame_at(0.09, 0.0).unwrap(),
        held.frame_at(0.09, 100.0).unwrap()
    );
    assert_eq!(
        held.frame_at(0.09, 1.0).unwrap().origin_velocity_world,
        [0.0; 2]
    );
}

#[test]
fn zero_motion_matches_original_native_incidence_and_independent_transport_metric() {
    let candidate = candidate();
    let config = controls(4);
    let held = SectionGait {
        heave_amplitude_chord: 0.0,
        pitch_amplitude_rad: 0.0,
        ..gait()
    };
    let report = screen_section_gait(&candidate, held, config, || false).unwrap();
    let mut fixed = WakeSim::new_panel_coupled(
        &candidate.section(config.panels),
        candidate.alpha,
        config.interval_nondimensional,
        config.core_radius_chord,
    )
    .unwrap();
    let GaitTransport::Fixed(options) = config.transport else {
        unreachable!()
    };
    let mut old = Vec::<([f64; 2], f64)>::new();
    let (s, c) = candidate.alpha.sin_cos();
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for _ in 0..config.intervals {
        fixed
            .step_panel_coupled_with_transport(candidate.alpha, options, || false)
            .unwrap();
        for ((position, gamma), new) in old.iter().zip(fixed.wake()) {
            let dx_world = c * (new.pos[0] - position[0]) + s * (new.pos[1] - position[1]);
            numerator += gamma.abs() * dx_world.abs();
            denominator += gamma.abs() * config.interval_nondimensional;
        }
        old = fixed.wake().iter().map(|v| (v.pos, v.gamma)).collect();
    }
    assert_eq!(report.wake.wake().len(), fixed.wake().len());
    for (a, b) in report.wake.wake().iter().zip(fixed.wake()) {
        close(a.pos[0], b.pos[0], 2e-12);
        close(a.pos[1], b.pos[1], 2e-12);
        close(a.gamma, b.gamma, 2e-12);
    }
    close(
        report.wake.panel_state().unwrap().bound_circulation,
        fixed.panel_state().unwrap().bound_circulation,
        2e-12,
    );
    close(
        report.circulation_weighted_streamwise_speed,
        numerator / denominator,
        3e-12,
    );
    close(report.circulation_time_weight, denominator, 3e-12);
    assert_eq!(
        report.coupled_solves,
        config.intervals * (2 * options.substeps + 1)
    );
    assert_eq!(
        report.transport_substeps,
        config.intervals * options.substeps
    );
    assert!(report.maximum_kelvin_residual < 2e-12);
}

#[test]
fn physical_gait_changes_accepted_wake_without_reinterpreting_legacy_gain_genes() {
    let candidate = candidate();
    let config = controls(4);
    let physical = screen_section_gait(&candidate, gait(), config, || false).unwrap();
    let mut different_legacy = candidate;
    different_legacy.flap_amp = 0.0;
    different_legacy.flap_freq = 0.2;
    different_legacy.inlet_x = 0.75;
    let repeated = screen_section_gait(&different_legacy, gait(), config, || false).unwrap();
    assert_eq!(
        physical.wake.trace_json(1).unwrap(),
        repeated.wake.trace_json(1).unwrap()
    );
    assert_eq!(position_error(&physical.wake, &repeated.wake), 0.0);
    assert_eq!(
        physical.circulation_weighted_streamwise_speed,
        repeated.circulation_weighted_streamwise_speed
    );
    let held = SectionGait {
        heave_amplitude_chord: 0.0,
        pitch_amplitude_rad: 0.0,
        ..gait()
    };
    let stationary = screen_section_gait(&candidate, held, config, || false).unwrap();
    assert!(position_error(&physical.wake, &stationary.wake) > 1e-4);
    assert!(
        (physical.wake.panel_state().unwrap().bound_circulation
            - stationary.wake.panel_state().unwrap().bound_circulation)
            .abs()
            > 1e-4
    );
    assert!(physical.maximum_kelvin_residual < 2e-12);
    assert!(physical.circulation_time_weight > 0.0);
    assert!(physical.circulation_weighted_streamwise_speed > 0.0);
    assert_eq!(physical.wake.history().len(), config.intervals);
    assert_eq!(
        physical.wake.panel_frame(),
        Some(
            gait()
                .frame_at(
                    candidate.alpha,
                    config.intervals as f64 * config.interval_nondimensional
                )
                .unwrap()
        )
    );
    // A packet first observed after the only interval has no sampled lifetime.
    let first = screen_section_gait(
        &candidate,
        gait(),
        GaitScreenControls {
            intervals: 1,
            ..config
        },
        || false,
    )
    .unwrap();
    assert_eq!(first.circulation_time_weight, 0.0);
    assert_eq!(first.circulation_weighted_streamwise_speed, 0.0);
}

#[test]
fn native_refinement_and_adaptive_accuracy_reach_the_physical_consumer() {
    let candidate = candidate();
    let reference = screen_section_gait(&candidate, gait(), controls(64), || false).unwrap();
    let mut errors = Vec::new();
    for substeps in [2, 4, 8] {
        let report = screen_section_gait(&candidate, gait(), controls(substeps), || false).unwrap();
        errors.push(position_error(&report.wake, &reference.wake));
    }
    assert!(errors[0] > 1e-9, "resolved temporal error: {errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.4 * pair[0], "midpoint convergence: {errors:?}");
    }
    let config = GaitScreenControls {
        transport: GaitTransport::Adaptive(PanelWakeAdaptiveOptions {
            initial_substeps: 1,
            max_substeps: 128,
            max_coupled_solves: 2048,
            absolute_position_tolerance: 1e-8,
            relative_position_tolerance: 0.0,
            max_displacement_fraction: None,
        }),
        ..controls(4)
    };
    let result = screen_section_gait(&candidate, gait(), config, || false).unwrap();
    assert!(position_error(&result.wake, &reference.wake) < errors[2]);
    assert!(result.coupled_solves > 8 * config.intervals);
    assert!(result.maximum_kelvin_residual < 2e-12);
    assert_eq!(result.wake.history().len(), config.intervals);
    assert_eq!(result.wake.wake().len(), reference.wake.wake().len());
    // The complete rejected-trial cost is admitted across the WHOLE screen.
    assert!(matches!(
        screen_section_gait(
            &candidate,
            gait(),
            GaitScreenControls {
                maximum_coupled_solves: result.coupled_solves - 1,
                ..config
            },
            || false
        ),
        Err(BemError::WorkEnvelopeExceeded { .. })
    ));
}

#[test]
fn invalid_gaits_work_exhaustion_and_mid_screen_cancellation_return_no_report() {
    let candidate = candidate();
    for invalid in [
        SectionGait {
            reduced_frequency: -1.0,
            ..gait()
        },
        SectionGait {
            pitch_amplitude_rad: f64::NAN,
            ..gait()
        },
        SectionGait {
            heave_phase_rad: f64::INFINITY,
            ..gait()
        },
        SectionGait {
            pivot_chord: 1.1,
            ..gait()
        },
    ] {
        assert!(screen_section_gait(&candidate, invalid, controls(4), || false).is_err());
    }
    let mut invalid_candidate = candidate;
    invalid_candidate.thickness = f64::NAN;
    assert!(screen_section_gait(&invalid_candidate, gait(), controls(4), || false).is_err());
    let config = controls(4);
    assert!(matches!(
        screen_section_gait(
            &candidate,
            gait(),
            GaitScreenControls {
                maximum_coupled_solves: 9 * config.intervals - 1,
                ..config
            },
            || false
        ),
        Err(BemError::WorkEnvelopeExceeded { .. })
    ));
    let mut visits = 0;
    screen_section_gait(&candidate, gait(), config, || {
        visits += 1;
        false
    })
    .unwrap();
    let mut checked = 0;
    assert!(matches!(
        screen_section_gait(&candidate, gait(), config, || {
            checked += 1;
            checked == visits / 2
        }),
        Err(BemError::Cancelled { .. })
    ));
    assert_eq!(checked, visits / 2);
}
