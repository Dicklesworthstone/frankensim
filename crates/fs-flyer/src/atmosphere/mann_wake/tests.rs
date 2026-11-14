use super::*;
use fs_atmo::mann::MannParams;
use fs_atmo::mann_field::{HomogeneousMannField, MannQuadrature, MannSpectralModel};
use fs_vpm::farfield3d::WakeCoreEvolutionMode;
use fs_vpm::filament3d::WakeRateCertificate;
use fs_vpm::images3d::ImagePlane;

#[path = "tests/affine.rs"]
mod affine;

#[path = "tests/rigid.rs"]
mod rigid;

const MODES: usize = 8;
const BYTES: usize = 4 * 88 * MODES;
const DT: f64 = 0.05;

fn source() -> FixedStepMannField {
    let model = MannSpectralModel::build(
        927,
        MannParams {
            alpha_eps23: 0.04,
            length_m: 1.0,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: MODES,
            min_kl: 0.03,
            max_kl: 3.0,
        },
        [-0.4, 0.05, 0.0],
    )
    .unwrap();
    let snapshot = HomogeneousMannField::at_anchor(model, 17, vec![0.02; MODES]).unwrap();
    FixedStepMannField::from_snapshot(snapshot, 91, DT, 88 * MODES).unwrap()
}

fn atmosphere_config() -> MannWakeConfig {
    MannWakeConfig {
        bridge_seed: 218,
        maximum_bridge_depth: 128,
        wall: None,
        mean: None,
    }
}

fn budget() -> MannWakeBudget {
    MannWakeBudget {
        maximum_normal_draws: 10_000_000,
        maximum_field_samples: 1_000_000,
    }
}

fn settings() -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [-0.2, 0.0, 0.0],
        emitted_core_radius_m: 0.2,
        core_mode: WakeCoreEvolutionMode::Frozen,
        substeps: 2,
        max_displacement_fraction: 1.0,
        max_work: 100_000_000,
    }
}

fn wake() -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 1.0 / DT,
            n_stations: 2,
            max_rows: 64,
        },
        vec![[0.0, -0.5, 1.0], [0.0, 0.0, 1.0], [0.0, 0.5, 1.0]],
    )
    .unwrap();
    for i in 0..4 {
        wake.shed(&[0.12 + i as f64 * 0.01, 0.18], [-0.2, 0.0, 0.0])
            .unwrap();
        wake.rows.last_mut().unwrap().core2_m2 = 0.04;
    }
    wake
}

fn adaptive() -> MannWakeIntegration {
    MannWakeIntegration::Adaptive {
        time: WakeTimeSpec {
            absolute_position_tolerance_m: 2e-5,
            relative_position_tolerance: 0.0,
            min_step_s: 1e-9,
            max_trials: 4096,
        },
        tree: None,
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}

fn assert_kelvin(wake: &FilamentWake) {
    for row in 0..=wake.rows.len() {
        for edge in 0..=wake.cert.n_stations {
            assert!(wake.node_circulation_residual(row, edge).unwrap().abs() < 2e-16);
        }
    }
}

// Independent wiring of the ORIGINAL native bridge and wake callback, with
// no stage cache: each spatial probe resamples the same conditional history.
fn bridge_reference(
    bridge: &MannOuBridge,
    frame: &mut MannOuBridgeFrame,
    point: [f64; 3],
    elapsed: f64,
) -> Result<[f64; 3], Refusal> {
    let fraction = elapsed / DT;
    let physical = if fraction == 1.0 {
        bridge.end_time_s()
    } else {
        bridge.start_time_s() + fraction * (bridge.end_time_s() - bridge.start_time_s())
    };
    bridge
        .sample_frame(physical, frame, 100_000, || false)
        .map_err(atmosphere)?;
    frame.sample(point).map(|flow| flow.u).map_err(atmosphere)
}

#[test]
fn native_midpoint_uses_one_spatial_history_and_changes_real_wake_trajectories() {
    let source = source();
    let original = wake();
    let mut actual = original.clone();
    let mut wind =
        MannWakeWind::prepare(source.clone(), &actual, atmosphere_config(), BYTES).unwrap();
    let bridge = MannOuBridge::prepare(&source, 218, 128, 88 * MODES, 4 * MODES, || false).unwrap();
    let mut frame = bridge.frame(88 * MODES, || false).unwrap();
    let mut expected = original.clone();
    expected
        .shed_self_induced_with_environment(
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            |point, elapsed| bridge_reference(&bridge, &mut frame, point, elapsed),
            || false,
        )
        .unwrap();
    let report = wind
        .shed_near(
            &mut actual,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            MannWakeIntegration::Midpoint { tree: None },
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(report.stage_frames, 3); // initial frame plus three changed stages
    assert!(report.field_samples > 3 * report.stage_frames);
    assert!(report.normal_draws > 8 * MODES);
    assert_eq!(report.start_time_s, source.time_s());
    let mut endpoint = source.clone();
    endpoint.advance_to(1, 4 * MODES, || false).unwrap();
    assert_eq!(wind.checkpoint().field, endpoint.checkpoint());
    let probe = [-0.3, 0.2, 1.1];
    assert_eq!(wind.sample(probe).unwrap(), endpoint.sample(probe).unwrap());
    assert_kelvin(&actual);
    for (row, before) in actual.rows.iter().zip(&original.rows) {
        assert_eq!(row.gamma, before.gamma);
        assert_eq!(row.core2_m2, before.core2_m2);
    }
    let mut frozen = original;
    frozen
        .shed_self_induced_with_environment(
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            |point, _| source.sample(point).map(|flow| flow.u).map_err(atmosphere),
            || false,
        )
        .unwrap();
    assert!(distance(frozen.rows[0].nodes[0], actual.rows[0].nodes[0]) > 1e-5);
}

#[test]
fn rejected_adaptive_stages_revisit_the_original_conditioned_path() {
    let source = source();
    let mut actual = wake();
    let mut expected = actual.clone();
    let mut wind =
        MannWakeWind::prepare(source.clone(), &actual, atmosphere_config(), BYTES).unwrap();
    let bridge = MannOuBridge::prepare(&source, 218, 128, 88 * MODES, 4 * MODES, || false).unwrap();
    let mut frame = bridge.frame(88 * MODES, || false).unwrap();
    let MannWakeIntegration::Adaptive { time, tree } = adaptive() else {
        unreachable!()
    };
    let mut query_order = Vec::new();
    let original_report = expected
        .shed_self_induced_adaptive(
            &[0.17, 0.2],
            settings(),
            time,
            WakeBoundary::default(),
            tree,
            |point, elapsed| {
                query_order.push(elapsed);
                bridge_reference(&bridge, &mut frame, point, elapsed)
            },
            || false,
        )
        .unwrap();
    assert!(original_report.rejected_trials > 0);
    assert!(query_order.windows(2).any(|pair| pair[1] < pair[0]));
    let report = wind
        .shed_near(
            &mut actual,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            adaptive(),
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(
        report.advection,
        MannWakeAdvection::Adaptive(original_report)
    );
    assert_eq!(report.field_samples, query_order.len());
    assert!(report.stage_frames < report.field_samples);
    assert!(report.maximum_bridge_depth_used > 2);
    assert_kelvin(&actual);
}

#[test]
fn near_and_far_sources_share_the_same_bridge_and_replay() {
    let mut actual = wake();
    let mut far =
        FarField::aggregate(&mut actual, 2, 1, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let source = source();
    let mut wind =
        MannWakeWind::prepare(source.clone(), &actual, atmosphere_config(), BYTES).unwrap();
    let mut expected = actual.clone();
    let mut expected_far = far.clone();
    let bridge = MannOuBridge::prepare(&source, 218, 128, 88 * MODES, 4 * MODES, || false).unwrap();
    let mut frame = bridge.frame(88 * MODES, || false).unwrap();
    let MannWakeIntegration::Adaptive { time, tree } = adaptive() else {
        unreachable!()
    };
    let original_report = expected_far
        .shed_self_induced_adaptive(
            &mut expected,
            &[0.17, 0.2],
            settings(),
            time,
            WakeBoundary::default(),
            tree,
            |point, elapsed| bridge_reference(&bridge, &mut frame, point, elapsed),
            || false,
        )
        .unwrap();
    let report = wind
        .shed_hybrid(
            &mut actual,
            &mut far,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            adaptive(),
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(far, expected_far);
    assert_eq!(
        report.advection,
        MannWakeAdvection::Adaptive(original_report)
    );

    let mut resumed = MannWakeWind::resume(wind.checkpoint(), BYTES).unwrap();
    let mut replay_near = actual.clone();
    let mut replay_far = far.clone();
    let integration = MannWakeIntegration::Midpoint {
        tree: Some(WakeTreeSpec {
            absolute_velocity_error_mps: 0.0,
            leaf_capacity: 2,
        }),
    };
    // At 2 substeps this shed moves a wake line 2.246e-2 m in one substep,
    // beyond fs-vpm's core/edge resolution limit of 1.041e-2 m, which refuses
    // it and asks for more substeps. The replay runs at a resolved count,
    // the same for the original and the resumed wind.
    let resolved = SelfInducedShed {
        substeps: 8,
        ..settings()
    };
    let a = wind
        .shed_hybrid(
            &mut actual,
            &mut far,
            &[0.15, 0.2],
            resolved,
            WakeBoundary::default(),
            integration,
            budget(),
            || false,
        )
        .unwrap();
    let b = resumed
        .shed_hybrid(
            &mut replay_near,
            &mut replay_far,
            &[0.15, 0.2],
            resolved,
            WakeBoundary::default(),
            integration,
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(actual, replay_near);
    assert_eq!(far, replay_far);
    assert_eq!(wind.checkpoint(), resumed.checkpoint());
}

#[test]
fn cancellation_and_both_work_limits_publish_neither_clock_nor_wake() {
    let original = wake();
    let source = source();
    let integration = MannWakeIntegration::Midpoint { tree: None };
    let mut reference = original.clone();
    let mut wind =
        MannWakeWind::prepare(source.clone(), &reference, atmosphere_config(), BYTES).unwrap();
    let mut polls = 0;
    wind.shed_near(
        &mut reference,
        &[0.17, 0.2],
        settings(),
        WakeBoundary::default(),
        integration,
        budget(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 100);
    for cancel_at in [1, 50, polls] {
        let mut near = original.clone();
        let mut wind =
            MannWakeWind::prepare(source.clone(), &near, atmosphere_config(), BYTES).unwrap();
        let before = wind.checkpoint();
        let mut seen = 0;
        assert!(
            wind.shed_near(
                &mut near,
                &[0.17, 0.2],
                settings(),
                WakeBoundary::default(),
                integration,
                budget(),
                || {
                    seen += 1;
                    seen == cancel_at
                },
            )
            .is_err()
        );
        assert_eq!(near, original);
        assert_eq!(wind.checkpoint(), before);
        wind.shed_near(
            &mut near,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            integration,
            budget(),
            || false,
        )
        .unwrap();
        assert_eq!(near, reference);
    }
    for limited in [
        MannWakeBudget {
            maximum_normal_draws: 8 * MODES,
            ..budget()
        },
        MannWakeBudget {
            maximum_field_samples: 1,
            ..budget()
        },
    ] {
        let mut near = original.clone();
        let mut far =
            FarField::aggregate(&mut near, 2, 1, WakeCoreEvolutionMode::Frozen, DT).unwrap();
        let mut wind =
            MannWakeWind::prepare(source.clone(), &near, atmosphere_config(), BYTES).unwrap();
        let (old_near, old_far, old_wind) = (near.clone(), far.clone(), wind.checkpoint());
        assert!(
            wind.shed_hybrid(
                &mut near,
                &mut far,
                &[0.17, 0.2],
                settings(),
                WakeBoundary::default(),
                integration,
                limited,
                || false,
            )
            .is_err()
        );
        assert_eq!(near, old_near);
        assert_eq!(far, old_far);
        assert_eq!(wind.checkpoint(), old_wind);
    }
}

#[test]
fn wall_and_log_mean_share_the_native_image_plane_and_analytic_gradient() {
    let source = source();
    let wall = MannPlaneWall::new(0.2).unwrap();
    let mean = FlatSiteLogLaw {
        scenario_effective_z0_m: 0.001,
        displacement_height_m: 0.0,
        reference_height_m: 2.0,
        reference_speed_mps: 0.5,
    };
    let config = MannWakeConfig {
        wall: Some(wall),
        mean: Some(mean),
        ..atmosphere_config()
    };
    let mut near = wake();
    let mut wind = MannWakeWind::prepare(source.clone(), &near, config, BYTES).unwrap();
    let point = [-0.3, 0.2, 0.9];
    let mut expected = source.sample_wall_blocked(point, wall).unwrap();
    expected.u[0] += mean.speed(point[2] - wall.height_m());
    expected.grad[0][2] += mean.dspeed_dh(point[2] - wall.height_m());
    assert_eq!(wind.sample(point).unwrap(), expected);
    assert_eq!(wind.sample([0.2, 0.3, 0.2]).unwrap().u[2], 0.0);
    let before = (near.clone(), wind.checkpoint());
    let mut native = settings();
    native.ambient_velocity_mps = [0.0; 3];
    let integration = MannWakeIntegration::Midpoint { tree: None };
    assert_eq!(
        wind.shed_near(
            &mut near,
            &[0.17, 0.2],
            native,
            WakeBoundary::default(),
            integration,
            budget(),
            || false,
        )
        .unwrap_err()
        .code,
        "mann-wake-boundary"
    );
    let boundary = WakeBoundary {
        line_endpoint_m: None,
        image_plane: Some(ImagePlane { z0_m: 0.2 }),
    };
    assert_eq!(
        wind.shed_near(
            &mut near,
            &[0.17, 0.2],
            settings(),
            boundary,
            integration,
            budget(),
            || false,
        )
        .unwrap_err()
        .code,
        "mann-wake-boundary"
    );
    assert_eq!(near, before.0);
    assert_eq!(wind.checkpoint(), before.1);
    wind.shed_near(
        &mut near,
        &[0.17, 0.2],
        native,
        boundary,
        integration,
        budget(),
        || false,
    )
    .unwrap();
    assert!(
        near.rows
            .iter()
            .flat_map(|row| &row.nodes)
            .all(|point| point[2] > 0.2)
    );
    assert_kelvin(&near);
}

#[test]
fn mismatched_clocks_and_storage_refuse_without_publication() {
    let mut near = wake();
    let source = source();
    assert_eq!(
        MannWakeWind::prepare(source.clone(), &near, atmosphere_config(), BYTES - 1)
            .unwrap_err()
            .code,
        "mann-wake-budget"
    );
    let mut wind =
        MannWakeWind::prepare(source.clone(), &near, atmosphere_config(), BYTES).unwrap();
    let before = wind.checkpoint();
    near.cert.shed_hz *= 2.0;
    assert_eq!(
        MannWakeWind::prepare(source, &near, atmosphere_config(), BYTES)
            .unwrap_err()
            .code,
        "mann-wake-clock"
    );
    let old_near = near.clone();
    assert_eq!(
        wind.shed_near(
            &mut near,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            MannWakeIntegration::Midpoint { tree: None },
            budget(),
            || false,
        )
        .unwrap_err()
        .code,
        "mann-wake-clock"
    );
    assert_eq!(wind.checkpoint(), before);
    assert_eq!(near, old_near);
    near.cert.shed_hz *= 0.5;
    near.ticks += 1;
    assert_eq!(
        wind.shed_near(
            &mut near,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            MannWakeIntegration::Midpoint { tree: None },
            budget(),
            || false,
        )
        .unwrap_err()
        .code,
        "mann-wake-clock"
    );
    assert_eq!(wind.checkpoint(), before);
}
