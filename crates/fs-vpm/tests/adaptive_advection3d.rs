//! Independent material-wake trajectories, temporal refinement and rollback.
use fs_vpm::farfield3d::WakeCoreEvolutionMode;
use fs_vpm::filament3d::{
    FilamentWake, Refusal, SelfInducedShed, WakeAdaptiveReport, WakeBoundary, WakeRateCertificate,
    WakeTimeSpec, WakeTreeSpec,
};
use fs_vpm::images3d::ImagePlane;

const SPAN: f64 = 1.4;
const CORE: f64 = 0.2;

fn ring(gamma: f64) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 2.0,
            n_stations: 1,
            max_rows: 128,
        },
        vec![[0.0, -SPAN / 2.0, 0.0], [0.0, SPAN / 2.0, 0.0]],
    )
    .unwrap();
    wake.shed(&[gamma], [0.0; 3]).unwrap();
    for p in &mut wake.rows[0].nodes {
        p[0] = -1.0;
    }
    wake.rows[0].core2_m2 = CORE * CORE;
    wake
}

fn settings(substeps: usize) -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [0.0; 3],
        emitted_core_radius_m: CORE,
        core_mode: WakeCoreEvolutionMode::Frozen,
        substeps,
        max_displacement_fraction: 1.0,
        max_work: 10_000_000,
    }
}

fn time(tolerance: f64) -> WakeTimeSpec {
    WakeTimeSpec {
        absolute_position_tolerance_m: tolerance,
        relative_position_tolerance: 0.0,
        min_step_s: 1e-10,
        max_trials: 8192,
    }
}

fn advance(wake: &mut FilamentWake, tolerance: f64) -> WakeAdaptiveReport {
    wake.shed_self_induced_adaptive(
        &[7.5],
        settings(1),
        time(tolerance),
        WakeBoundary::default(),
        None,
        |_, _| Ok([0.0; 3]),
        || false,
    )
    .unwrap()
}

// Integrate the two opposite straight edges at a free rectangle corner.
// Its two incident edges induce zero; the remaining velocity rotates the
// free span about the fixed lifting line at this constant angular rate.
fn angular_rate(gamma: f64) -> f64 {
    gamma * SPAN / (4.0 * std::f64::consts::PI * (1.0 + SPAN * SPAN + CORE * CORE).sqrt())
        * (1.0 / (1.0 + CORE * CORE) + 1.0 / (SPAN * SPAN + CORE * CORE))
}

fn expected(gamma: f64, y: f64) -> [f64; 3] {
    let angle = angular_rate(gamma) / 2.0;
    [-angle.cos(), y, angle.sin()]
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}

fn assert_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    assert!(distance(a, b) <= tolerance, "{a:?} != {b:?}");
}

#[test]
fn tolerance_refines_actual_roll_up_and_preserves_kelvin_balance_and_emission() {
    let original = ring(4.0);
    let mut errors = Vec::new();
    let mut counts = Vec::new();
    for tolerance in [1e-5, 1e-8, 1e-11] {
        let mut wake = original.clone();
        let report = advance(&mut wake, tolerance);
        errors.push(distance(wake.rows[0].nodes[0], expected(4.0, -SPAN / 2.0)));
        counts.push(report.accepted_substeps);
        assert!(report.rejected_trials > 0, "{report:?}");
        assert!(report.max_accepted_position_error_m <= tolerance);
        assert!(report.max_accepted_error_ratio <= 1.0);
        assert_eq!(wake.line_nodes, original.line_nodes);
        assert_eq!(wake.rows.len(), 2);
        assert_eq!(wake.ticks, original.ticks + 1);
        assert_eq!(wake.rows[0].gamma, [4.0]);
        assert_eq!(wake.rows[0].age_ticks, 1.0);
        assert_eq!(wake.rows[0].core2_m2, CORE * CORE);
        assert_eq!(wake.rows[1].gamma, [7.5]);
        assert_eq!(wake.rows[1].age_ticks, 0.0);
        assert_eq!(wake.rows[1].nodes, wake.line_nodes);
        for row in 0..=wake.rows.len() {
            for edge in 0..=wake.cert.n_stations {
                assert_eq!(wake.node_circulation_residual(row, edge), Some(0.0));
            }
        }
    }
    assert!(
        errors[1] < errors[0] / 20.0 && errors[2] < errors[1] / 20.0,
        "{errors:?}"
    );
    assert!(errors[2] < 5e-10, "{errors:?}");
    assert!(counts[0] < counts[1] && counts[1] < counts[2], "{counts:?}");
}

#[test]
fn accepted_third_order_endpoints_converge_at_the_expected_order() {
    let mut errors = Vec::new();
    for substeps in [2, 4, 8] {
        let mut wake = ring(4.0);
        let report = wake
            .shed_self_induced_adaptive(
                &[4.0],
                settings(substeps),
                time(1.0),
                WakeBoundary::default(),
                None,
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap();
        assert_eq!(report.accepted_substeps, substeps as u64);
        assert_eq!(report.rejected_trials, 0);
        errors.push(distance(wake.rows[0].nodes[0], expected(4.0, -SPAN / 2.0)));
    }
    for pair in errors.windows(2) {
        assert!((7.0..9.0).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn adaptive_evolution_resolves_the_fixed_step_endpoint_collapse_under_strain() {
    let mut original = ring(1e-30);
    original.cert.shed_hz = 1.0;
    original.line_nodes = vec![[1.0, -1.0, 2.0], [1.0, 1.0, 2.0]];
    original.rows[0].nodes = vec![[0.0, -1.0, 2.0], [0.0, 1.0, 2.0]];
    original.rows[0].core2_m2 = 4.0;
    let mut config = settings(1);
    config.emitted_core_radius_m = 2.0;
    // This divergence-free incident flow has y(t)=y(0)*exp(-t^2).
    let strain = |p: [f64; 3], t: f64| Ok([2.0 * t * p[0], -2.0 * t * p[1], 0.0]);
    let mut fixed = original.clone();
    assert_eq!(
        fixed
            .shed_self_induced_with_environment(
                &[1e-30],
                config,
                WakeBoundary::default(),
                strain,
                || false,
            )
            .unwrap_err()
            .code,
        "wake-advection-core-invalid"
    );
    assert_eq!(fixed, original);
    let mut adaptive = original;
    let report = adaptive
        .shed_self_induced_adaptive(
            &[1e-30],
            config,
            time(1e-10),
            WakeBoundary::default(),
            None,
            strain,
            || false,
        )
        .unwrap();
    assert!(report.rejected_trials > 0 && report.accepted_substeps > 1);
    for (p, sign) in adaptive.rows[0].nodes.iter().zip([-1.0, 1.0]) {
        assert_close(*p, [0.0, sign * (-1.0_f64).exp(), 2.0], 3e-8);
    }
}

#[test]
fn coarse_image_plane_stage_crossings_refine_boundary_compatible_compression() {
    let mut wake = ring(1e-30);
    wake.cert.shed_hz = 1.0;
    wake.line_nodes = vec![[1.0, -1.0, 0.001], [1.0, 1.0, 0.001]];
    wake.rows[0].nodes = vec![[0.0, -1.0, 0.001], [0.0, 1.0, 0.001]];
    wake.rows[0].core2_m2 = 4.0;
    let mut config = settings(1);
    config.emitted_core_radius_m = 2.0;
    // div u=0 and u_z=0 at the plane. Exact z(t)=0.001*exp(-4t)
    // remains positive. The first h=1 and h=1/2 midpoint stages cross or
    // reach the plane despite full displacement being far below core/edge
    // limits; these are numerical trial failures, not invalid initial data.
    let report = wake
        .shed_self_induced_adaptive(
            &[1e-30],
            config,
            time(1e-12),
            WakeBoundary {
                line_endpoint_m: None,
                image_plane: Some(ImagePlane { z0_m: 0.0 }),
            },
            None,
            |p, _| Ok([4.0 * p[0], 0.0, -4.0 * p[2]]),
            || false,
        )
        .unwrap();
    assert!(report.rejected_trials >= 2, "{report:?}");
    for (p, sign) in wake.rows[0].nodes.iter().zip([-1.0, 1.0]) {
        assert!(p[2] > 0.0);
        assert_close(*p, [0.0, sign, 0.001 * (-4.0_f64).exp()], 2e-10);
    }
}

#[test]
fn rapid_roll_up_and_prescribed_line_motion_refine_the_resolution_guard() {
    let original = ring(16.0);
    let velocity = [-1.0, 0.0, 0.25];
    let endpoint: Vec<_> = original
        .line_nodes
        .iter()
        .map(|p| std::array::from_fn(|c| p[c] + 0.5 * velocity[c]))
        .collect();
    let boundary = WakeBoundary {
        line_endpoint_m: Some(&endpoint),
        image_plane: None,
    };
    let mut config = settings(1);
    config.ambient_velocity_mps = velocity;
    let mut fixed = original.clone();
    assert_eq!(
        fixed
            .shed_self_induced_with_environment(
                &[16.0],
                config,
                boundary,
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap_err()
            .code,
        "wake-advection-displacement"
    );
    let mut wake = original;
    let report = wake
        .shed_self_induced_adaptive(
            &[16.0],
            config,
            time(1e-10),
            boundary,
            None,
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    assert!(report.rejected_trials >= 2, "{report:?}");
    assert!(report.max_accepted_step_s < 0.125);
    for p in &wake.rows[0].nodes {
        let rotated = expected(16.0, p[1]);
        assert_close(
            *p,
            std::array::from_fn(|c| rotated[c] + 0.5 * velocity[c]),
            5e-8,
        );
    }
    assert_eq!(wake.line_nodes, endpoint);
    assert_eq!(wake.rows[1].nodes, endpoint);
}

#[test]
fn relative_time_accuracy_does_not_depend_on_the_world_origin() {
    let original = ring(4.0);
    let shift = [1048576.0, -1048576.0, 1048576.0];
    let mut shifted = original.clone();
    for p in shifted
        .line_nodes
        .iter_mut()
        .chain(shifted.rows.iter_mut().flat_map(|r| &mut r.nodes))
    {
        for c in 0..3 {
            p[c] += shift[c];
        }
    }
    let control = WakeTimeSpec {
        relative_position_tolerance: 1e-4,
        ..time(1e-12)
    };
    let run = |wake: &mut FilamentWake| {
        wake.shed_self_induced_adaptive(
            &[4.0],
            settings(1),
            control,
            WakeBoundary::default(),
            None,
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap()
    };
    let mut base = original;
    let a = run(&mut base);
    let b = run(&mut shifted);
    assert!(a.accepted_substeps > 2);
    assert_eq!(a.accepted_substeps, b.accepted_substeps);
    assert_eq!(a.rejected_trials, b.rejected_trials);
    for (a, b) in base.rows[0].nodes.iter().zip(&shifted.rows[0].nodes) {
        assert_close(*a, std::array::from_fn(|c| b[c] - shift[c]), 2e-8);
    }
}

#[test]
fn uniform_transport_reuses_endpoint_fields_and_emits_only_once() {
    let mut wake = ring(0.0);
    wake.rows[0].age_ticks = 2.25;
    let original = wake.clone();
    let mut config = settings(10);
    config.ambient_velocity_mps = [0.125, -0.5, 0.25];
    let mut samples = 0;
    let report = wake
        .shed_self_induced_adaptive(
            &[3.0],
            config,
            time(1e-20),
            WakeBoundary::default(),
            None,
            |_, _| {
                samples += 1;
                Ok([0.5, 0.0, 0.125])
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.accepted_substeps, 10);
    assert_eq!(report.rejected_trials, 0);
    assert_eq!(report.max_accepted_position_error_m, 0.0);
    assert_eq!(samples, 2 * (1 + 3 * 10));
    assert_eq!(report.advection.direct_interactions, 0);
    for (old, new) in original.rows[0].nodes.iter().zip(&wake.rows[0].nodes) {
        assert_close(
            *new,
            [old[0] + 0.3125, old[1] - 0.25, old[2] + 0.1875],
            1e-15,
        );
    }
    assert_eq!(wake.rows.len(), 2);
    assert_eq!(wake.ticks, original.ticks + 1);
    assert_eq!(wake.rows[0].age_ticks, 3.25);
    assert_eq!(wake.rows[1].age_ticks, 0.0);
    assert_eq!(wake.rows[1].gamma, [3.0]);
}

#[test]
fn rejected_trial_work_is_cumulative_and_late_failures_publish_nothing() {
    let original = ring(4.0);
    let control = time(1e-10);
    let mut completed = original.clone();
    let mut full_checks = 0;
    let report = completed
        .shed_self_induced_adaptive(
            &[4.0],
            settings(1),
            control,
            WakeBoundary::default(),
            None,
            |_, _| Ok([0.0; 3]),
            || {
                full_checks += 1;
                false
            },
        )
        .unwrap();
    assert!(report.rejected_trials > 0);
    let mut config = settings(1);
    config.max_work = report.advection.work_used;
    let mut exact_cap = original.clone();
    let repeated = exact_cap
        .shed_self_induced_adaptive(
            &[4.0],
            config,
            control,
            WakeBoundary::default(),
            Some(WakeTreeSpec {
                absolute_velocity_error_mps: 0.0,
                leaf_capacity: 8,
            }),
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    assert_eq!(exact_cap, completed);
    assert_eq!(repeated, report);
    config.max_work -= 1;
    let mut exhausted = original.clone();
    assert_eq!(
        exhausted
            .shed_self_induced_adaptive(
                &[4.0],
                config,
                control,
                WakeBoundary::default(),
                None,
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap_err()
            .code,
        "wake-advection-work-exhausted"
    );
    assert_eq!(exhausted, original);
    let mut cancelled = original.clone();
    let mut checks = 0;
    assert_eq!(
        cancelled
            .shed_self_induced_adaptive(
                &[4.0],
                settings(1),
                control,
                WakeBoundary::default(),
                None,
                |_, _| Ok([0.0; 3]),
                || {
                    checks += 1;
                    checks == full_checks
                },
            )
            .unwrap_err()
            .code,
        "wake-advection-cancelled"
    );
    assert_eq!(cancelled, original);
    let mut capped = original.clone();
    let cap = WakeTimeSpec {
        max_trials: report.accepted_substeps + report.rejected_trials - 1,
        ..control
    };
    assert_eq!(
        capped
            .shed_self_induced_adaptive(
                &[4.0],
                settings(1),
                cap,
                WakeBoundary::default(),
                None,
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap_err()
            .code,
        "wake-advection-trials-exhausted"
    );
    assert_eq!(capped, original);
}

#[test]
fn invalid_time_controls_and_minimum_steps_refuse_atomically() {
    let original = ring(4.0);
    let base = time(1e-12);
    for control in [
        WakeTimeSpec {
            absolute_position_tolerance_m: 0.0,
            ..base
        },
        WakeTimeSpec {
            absolute_position_tolerance_m: f64::NAN,
            ..base
        },
        WakeTimeSpec {
            relative_position_tolerance: -1.0,
            ..base
        },
        WakeTimeSpec {
            min_step_s: 0.0,
            ..base
        },
        WakeTimeSpec {
            min_step_s: 1.0,
            ..base
        },
        WakeTimeSpec {
            max_trials: 0,
            ..base
        },
    ] {
        let mut wake = original.clone();
        assert_eq!(
            wake.shed_self_induced_adaptive(
                &[4.0],
                settings(1),
                control,
                WakeBoundary::default(),
                None,
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap_err()
            .code,
            "wake-advection-time-invalid"
        );
        assert_eq!(wake, original);
    }
    let mut wake = original.clone();
    assert_eq!(
        wake.shed_self_induced_adaptive(
            &[4.0],
            settings(1),
            WakeTimeSpec {
                min_step_s: 0.5,
                ..base
            },
            WakeBoundary::default(),
            None,
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap_err()
        .code,
        "wake-advection-time-step-exhausted"
    );
    assert_eq!(wake, original);
}

#[test]
fn original_incident_refusals_never_turn_into_internal_retries() {
    let original = ring(0.0);
    let refusal = Refusal {
        code: "wake-advection-displacement",
        message: "original incident source unavailable".into(),
        ranked_repairs: vec!["restore original incident source".into()],
    };
    let mut wake = original.clone();
    let mut config = settings(4);
    config.ambient_velocity_mps = [0.125, 0.0, 0.0];
    let error = wake
        .shed_self_induced_adaptive(
            &[0.0],
            config,
            time(1e-12),
            WakeBoundary::default(),
            None,
            |_, elapsed| {
                if elapsed >= 0.2 {
                    Err(refusal.clone())
                } else {
                    Ok([0.0; 3])
                }
            },
            || false,
        )
        .unwrap_err();
    assert_eq!(error, refusal);
    assert_eq!(wake, original);
}
