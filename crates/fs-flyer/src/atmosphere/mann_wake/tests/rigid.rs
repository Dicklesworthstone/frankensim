//! Real rotating native sources driven by the original conditioned wind.
use super::*;
use fs_vpm::filament3d::RigidLineMotion;

fn boundary() -> RigidWakeBoundary {
    RigidWakeBoundary {
        line_motion: RigidLineMotion {
            pivot_start_m: [0.0, 0.0, 1.0],
            pivot_velocity_mps: [0.025, -0.015, 0.01],
            angular_velocity_world_rad_s: [0.4, -0.3, 4.0],
        },
        image_plane: None,
    }
}

fn resolved() -> SelfInducedShed {
    SelfInducedShed { substeps: 32, ..settings() }
}

fn midpoint(tree: Option<WakeTreeSpec>) -> MannWakeIntegration {
    MannWakeIntegration::Midpoint { tree }
}

#[test]
fn zero_rigid_motion_preserves_the_original_adapter_transaction() {
    let original = wake();
    let mut legacy = original.clone();
    let mut rigid = original;
    let mut old_wind = MannWakeWind::prepare(source(), &legacy, atmosphere_config(), BYTES).unwrap();
    let mut new_wind = MannWakeWind::prepare(source(), &rigid, atmosphere_config(), BYTES).unwrap();
    let old_report = old_wind.shed_near(
        &mut legacy, &[0.17, 0.2], settings(), WakeBoundary::default(),
        midpoint(None), budget(), || false,
    ).unwrap();
    let new_report = new_wind.shed_near_rigid(
        &mut rigid, &[0.17, 0.2], settings(), RigidWakeBoundary::default(),
        midpoint(None), MannWakeSampling::Exact, budget(), || false,
    ).unwrap();
    assert_eq!(rigid, legacy);
    assert_eq!(new_report, old_report);
    assert_eq!(new_wind.checkpoint(), old_wind.checkpoint());
}

#[test]
fn rotating_near_wake_matches_native_stage_wind_and_changes_material_trajectories() {
    let source = source();
    let original = wake();
    let motion = boundary();
    for tree in [
        None,
        Some(WakeTreeSpec { absolute_velocity_error_mps: 0.0, leaf_capacity: 2 }),
    ] {
        let mut actual = original.clone();
        let mut expected = original.clone();
        let mut wind = MannWakeWind::prepare(
            source.clone(), &actual, atmosphere_config(), BYTES,
        ).unwrap();
        let bridge = MannOuBridge::prepare(
            &source, 218, 128, 88 * MODES, 4 * MODES, || false,
        ).unwrap();
        let mut frame = bridge.frame(88 * MODES, || false).unwrap();
        let native_report = expected.shed_self_induced_with_rigid_line(
            &[0.17, 0.2], resolved(), motion, tree,
            |point, elapsed| bridge_reference(&bridge, &mut frame, point, elapsed),
            || false,
        ).unwrap();
        let report = wind.shed_near_rigid(
            &mut actual, &[0.17, 0.2], resolved(), motion, midpoint(tree),
            MannWakeSampling::Exact, budget(), || false,
        ).unwrap();
        assert_eq!(actual, expected);
        let expected_report = match tree {
            Some(_) => MannWakeAdvection::TreeMidpoint(native_report),
            None => MannWakeAdvection::DirectMidpoint,
        };
        assert_eq!(report.advection, expected_report);
        assert_eq!(report.stage_frames, 2 * resolved().substeps - 1);
        assert_eq!(
            report.field_samples,
            2 * resolved().substeps * original.rows.len() * (original.cert.n_stations + 1),
        );
        assert_eq!(report.affine_frames, 0);
        let mut endpoint = source.clone();
        endpoint.advance_to(1, 4 * MODES, || false).unwrap();
        assert_eq!(wind.checkpoint().field, endpoint.checkpoint());
        for (&start, &end) in original.line_nodes.iter().zip(&actual.line_nodes) {
            assert_eq!(end, motion.line_motion.point_and_velocity(start, DT).unwrap().0);
        }
        assert_eq!(actual.rows.last().unwrap().nodes, actual.line_nodes);
        assert_kelvin(&actual);

        // Hold the original line in the SAME evolving wind. The changed
        // free-vertex trajectory therefore comes from the rigid source motion.
        let mut fixed = original.clone();
        let mut fixed_wind = MannWakeWind::prepare(
            source.clone(), &fixed, atmosphere_config(), BYTES,
        ).unwrap();
        fixed_wind.shed_near(
            &mut fixed, &[0.17, 0.2], resolved(), WakeBoundary::default(),
            midpoint(tree), budget(), || false,
        ).unwrap();
        assert!(distance(actual.rows[0].nodes[2], fixed.rows[0].nodes[2]) > 1e-7);
        assert_eq!(wind.checkpoint(), fixed_wind.checkpoint());
    }
}

#[test]
fn rotating_hybrid_adaptive_trials_share_the_native_history_and_replay() {
    let mut near = wake();
    let mut far = FarField::aggregate(
        &mut near, 2, 1, WakeCoreEvolutionMode::Frozen, DT,
    ).unwrap();
    let initial_near = near.clone();
    let initial_far = far.clone();
    let source = source();
    let mut wind = MannWakeWind::prepare(
        source.clone(), &near, atmosphere_config(), BYTES,
    ).unwrap();
    let checkpoint = wind.checkpoint();
    let mut expected = near.clone();
    let mut expected_far = far.clone();
    let bridge = MannOuBridge::prepare(
        &source, 218, 128, 88 * MODES, 4 * MODES, || false,
    ).unwrap();
    let mut frame = bridge.frame(88 * MODES, || false).unwrap();
    let MannWakeIntegration::Adaptive { time, tree } = adaptive() else { unreachable!() };
    let mut samples = 0;
    let native_report = expected_far.shed_self_induced_adaptive_with_rigid_line(
        &mut expected, &[0.17, 0.2], settings(), time, boundary(), tree,
        |point, elapsed| {
            samples += 1;
            bridge_reference(&bridge, &mut frame, point, elapsed)
        },
        || false,
    ).unwrap();
    let report = wind.shed_hybrid_rigid(
        &mut near, &mut far, &[0.17, 0.2], settings(), boundary(), adaptive(),
        MannWakeSampling::Exact, budget(), || false,
    ).unwrap();
    assert_eq!(near, expected);
    assert_eq!(far, expected_far);
    assert_eq!(report.advection, MannWakeAdvection::Adaptive(native_report));
    assert_eq!(report.field_samples, samples);
    assert!(native_report.rejected_trials > 0);
    assert!(report.stage_frames < report.field_samples);
    assert_ne!(far, initial_far);

    let mut replay_near = initial_near;
    let mut replay_far = initial_far;
    let mut replay_wind = MannWakeWind::resume(checkpoint, BYTES).unwrap();
    let repeated = replay_wind.shed_hybrid_rigid(
        &mut replay_near, &mut replay_far, &[0.17, 0.2], settings(), boundary(), adaptive(),
        MannWakeSampling::Exact, budget(), || false,
    ).unwrap();
    assert_eq!(repeated, report);
    assert_eq!(replay_near, near);
    assert_eq!(replay_far, far);
    assert_eq!(replay_wind.checkpoint(), wind.checkpoint());
}

#[test]
fn rotating_adaptive_affine_sampling_reuses_each_native_stage_and_keeps_its_error_budget() {
    let source = source();
    let mut actual = wake();
    let mut expected = actual.clone();
    let mut wind = MannWakeWind::prepare(
        source.clone(), &actual, atmosphere_config(), BYTES,
    ).unwrap();
    let anchor = [0.0, 0.0, 1.0];
    let accuracy = MannAffineAccuracy {
        maximum_velocity_error_mps: 10.0,
        maximum_gradient_error_per_s: 10.0,
    };
    let sampling = MannWakeSampling::Affine { anchor_m: anchor, accuracy };
    let bridge = MannOuBridge::prepare(
        &source, 218, 128, 88 * MODES, 4 * MODES, || false,
    ).unwrap();
    let mut frame = bridge.frame(88 * MODES, || false).unwrap();
    let MannWakeIntegration::Adaptive { time, tree } = adaptive() else { unreachable!() };
    let mut maximum = MannAffineError::default();
    let mut samples = 0;
    let native_report = expected.shed_self_induced_adaptive_with_rigid_line(
        &[0.17, 0.2], settings(), time, boundary(), tree,
        |point, elapsed| {
            samples += 1;
            bridge_reference(&bridge, &mut frame, point, elapsed)?;
            // Deliberately uncached: rebuild the original affine observation
            // on every actual native query, including revisited trial times.
            let observation = frame.affine_observation(anchor, None, || false).map_err(atmosphere)?;
            let sample = observation.sample(point, accuracy).map_err(atmosphere)?;
            maximum.velocity_mps = maximum.velocity_mps.max(sample.error.velocity_mps);
            maximum.gradient_per_s = maximum.gradient_per_s.max(sample.error.gradient_per_s);
            Ok(sample.flow.u)
        },
        || false,
    ).unwrap();
    let report = wind.shed_near_rigid(
        &mut actual, &[0.17, 0.2], settings(), boundary(), adaptive(),
        sampling, budget(), || false,
    ).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(report.advection, MannWakeAdvection::Adaptive(native_report));
    assert_eq!(report.field_samples, samples);
    assert_eq!(report.maximum_affine_error, maximum);
    assert_eq!(report.affine_frames, report.stage_frames + 1);
    assert!(report.field_samples > 3 * report.affine_frames);
    assert!(native_report.rejected_trials > 0);
    assert_kelvin(&actual);
}

#[test]
fn rigid_motion_work_spatial_boundary_and_final_cancellation_refusals_are_joint() {
    let mut near = wake();
    let mut far = FarField::aggregate(
        &mut near, 2, 1, WakeCoreEvolutionMode::Frozen, DT,
    ).unwrap();
    let initial_near = near.clone();
    let initial_far = far.clone();
    let mut wind = MannWakeWind::prepare(source(), &near, atmosphere_config(), BYTES).unwrap();
    let checkpoint = wind.checkpoint();
    let limited = MannWakeBudget { maximum_field_samples: 1, ..budget() };
    assert_eq!(wind.shed_hybrid_rigid(
        &mut near, &mut far, &[0.17, 0.2], resolved(), boundary(), midpoint(None),
        MannWakeSampling::Exact, limited, || false,
    ).unwrap_err().code, "mann-wake-work");
    assert_eq!(near, initial_near);
    assert_eq!(far, initial_far);
    assert_eq!(wind.checkpoint(), checkpoint);

    let exact_only = MannWakeSampling::Affine {
        anchor_m: near.rows[0].nodes[0],
        accuracy: MannAffineAccuracy {
            maximum_velocity_error_mps: 0.0,
            maximum_gradient_error_per_s: 0.0,
        },
    };
    assert_eq!(wind.shed_hybrid_rigid(
        &mut near, &mut far, &[0.17, 0.2], resolved(), boundary(), midpoint(None),
        exact_only, budget(), || false,
    ).unwrap_err().code, "mann-affine-accuracy");
    assert_eq!(near, initial_near);
    assert_eq!(far, initial_far);
    assert_eq!(wind.checkpoint(), checkpoint);

    let mismatch = RigidWakeBoundary {
        image_plane: Some(ImagePlane { z0_m: 0.0 }),
        ..boundary()
    };
    assert_eq!(wind.shed_hybrid_rigid(
        &mut near, &mut far, &[0.17, 0.2], resolved(), mismatch, midpoint(None),
        MannWakeSampling::Exact, budget(), || false,
    ).unwrap_err().code, "mann-wake-boundary");
    assert_eq!(near, initial_near);
    assert_eq!(wind.checkpoint(), checkpoint);

    let mut polls = 0;
    let accepted = wind.shed_hybrid_rigid(
        &mut near, &mut far, &[0.17, 0.2], resolved(), boundary(), midpoint(None),
        MannWakeSampling::Exact, budget(), || { polls += 1; false },
    ).unwrap();
    assert!(accepted.field_samples > 1);
    let mut refused_near = initial_near.clone();
    let mut refused_far = initial_far.clone();
    let mut cancelled = MannWakeWind::resume(checkpoint.clone(), BYTES).unwrap();
    let mut seen = 0;
    assert!(cancelled.shed_hybrid_rigid(
        &mut refused_near, &mut refused_far, &[0.17, 0.2], resolved(), boundary(), midpoint(None),
        MannWakeSampling::Exact, budget(), || { seen += 1; seen == polls },
    ).is_err());
    assert_eq!(seen, polls);
    assert_eq!(refused_near, initial_near);
    assert_eq!(refused_far, initial_far);
    assert_eq!(cancelled.checkpoint(), checkpoint);
}
