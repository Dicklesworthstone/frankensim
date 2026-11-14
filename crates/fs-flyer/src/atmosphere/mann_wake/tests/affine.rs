//! Actual native wake transactions under a budgeted spatial approximation.
use super::*;

fn sampling(anchor_m: [f64; 3]) -> MannWakeSampling {
    MannWakeSampling::Affine {
        anchor_m,
        accuracy: MannAffineAccuracy {
            maximum_velocity_error_mps: 10.0,
            maximum_gradient_error_per_s: 10.0,
        },
    }
}

// Independently rebuild on EVERY actual native query, with no affine cache.
// The original bridge addresses the same conditional history even after a
// rejected trial returns to an earlier stage. This also checks the spatial
// error estimate at the actual moving native vertices, not only fixed probes.
fn reference(
    bridge: &MannOuBridge,
    frame: &mut MannOuBridgeFrame,
    point: [f64; 3],
    elapsed: f64,
    config: MannWakeConfig,
    sampling: MannWakeSampling,
    maximum: &mut MannAffineError,
) -> Result<[f64; 3], Refusal> {
    bridge_reference(bridge, frame, point, elapsed)?;
    let MannWakeSampling::Affine { anchor_m, accuracy } = sampling else {
        unreachable!()
    };
    let observation = frame
        .affine_observation(anchor_m, config.wall, || false)
        .map_err(atmosphere)?;
    let sample = observation.sample(point, accuracy).map_err(atmosphere)?;
    let exact = match config.wall {
        Some(wall) => frame.sample_wall_blocked(point, wall),
        None => frame.sample(point),
    }
    .map_err(atmosphere)?;
    assert!(distance(exact.u, sample.flow.u) <= sample.error.velocity_mps + 2e-13);
    maximum.velocity_mps = maximum.velocity_mps.max(sample.error.velocity_mps);
    maximum.gradient_per_s = maximum.gradient_per_s.max(sample.error.gradient_per_s);
    add_mean(sample.flow, point, config).map(|flow| flow.u)
}

#[test]
fn native_adaptive_affine_cache_matches_uncached_bridge_including_rejected_stages() {
    let source = source();
    let config = atmosphere_config();
    let sampling = sampling([0.0, 0.0, 1.0]);
    let mut actual = wake();
    let mut expected = actual.clone();
    let mut wind = MannWakeWind::prepare(source.clone(), &actual, config, BYTES).unwrap();
    let bridge = MannOuBridge::prepare(
        &source,
        config.bridge_seed,
        128,
        88 * MODES,
        4 * MODES,
        || false,
    )
    .unwrap();
    let mut frame = bridge.frame(88 * MODES, || false).unwrap();
    let MannWakeIntegration::Adaptive { time, tree } = adaptive() else {
        unreachable!()
    };
    let mut query_order = Vec::new();
    let mut maximum = MannAffineError::default();
    let direct_report = expected
        .shed_self_induced_adaptive(
            &[0.17, 0.2],
            settings(),
            time,
            WakeBoundary::default(),
            tree,
            |point, elapsed| {
                query_order.push(elapsed);
                reference(
                    &bridge,
                    &mut frame,
                    point,
                    elapsed,
                    config,
                    sampling,
                    &mut maximum,
                )
            },
            || false,
        )
        .unwrap();
    let report = wind
        .shed_near_with_sampling(
            &mut actual,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            adaptive(),
            sampling,
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(report.advection, MannWakeAdvection::Adaptive(direct_report));
    assert_eq!(report.maximum_affine_error, maximum);
    assert!(direct_report.rejected_trials > 0);
    assert!(query_order.windows(2).any(|times| times[1] < times[0]));
    assert_eq!(report.affine_frames, report.stage_frames + 1);
    assert!(report.field_samples > 3 * report.affine_frames);
    assert_eq!(report.field_samples, query_order.len());
    assert_kelvin(&actual);
    let mut endpoint = source;
    endpoint.advance_to(1, 4 * MODES, || false).unwrap();
    assert_eq!(wind.checkpoint().field, endpoint.checkpoint());
}

#[test]
fn wall_hybrid_affine_sampling_retains_exact_mean_and_checkpoint_replay() {
    let mut near = wake();
    let mut far = FarField::aggregate(&mut near, 2, 1, WakeCoreEvolutionMode::Frozen, DT).unwrap();
    let config = MannWakeConfig {
        wall: Some(MannPlaneWall::new(0.2).unwrap()),
        mean: Some(FlatSiteLogLaw {
            scenario_effective_z0_m: 0.001,
            displacement_height_m: 0.0,
            reference_height_m: 2.0,
            reference_speed_mps: 0.5,
        }),
        ..atmosphere_config()
    };
    let sampling = sampling([0.0, 0.0, 0.2]);
    let boundary = WakeBoundary {
        line_endpoint_m: None,
        image_plane: Some(ImagePlane { z0_m: 0.2 }),
    };
    let native = SelfInducedShed {
        ambient_velocity_mps: [0.0; 3],
        ..settings()
    };
    let integration = MannWakeIntegration::Midpoint { tree: None };
    let source = source();
    let mut wind = MannWakeWind::prepare(source.clone(), &near, config, BYTES).unwrap();
    let mut expected = near.clone();
    let mut expected_far = far.clone();
    let bridge = MannOuBridge::prepare(
        &source,
        config.bridge_seed,
        128,
        88 * MODES,
        4 * MODES,
        || false,
    )
    .unwrap();
    let mut frame = bridge.frame(88 * MODES, || false).unwrap();
    let mut maximum = MannAffineError::default();
    expected_far
        .shed_self_induced_with_environment(
            &mut expected,
            &[0.17, 0.2],
            native,
            boundary,
            |point, elapsed| {
                reference(
                    &bridge,
                    &mut frame,
                    point,
                    elapsed,
                    config,
                    sampling,
                    &mut maximum,
                )
            },
            || false,
        )
        .unwrap();
    let report = wind
        .shed_hybrid_with_sampling(
            &mut near,
            &mut far,
            &[0.17, 0.2],
            native,
            boundary,
            integration,
            sampling,
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(near, expected);
    assert_eq!(far, expected_far);
    assert_eq!(report.maximum_affine_error, maximum);
    assert!(report.affine_frames > 1);
    assert!(report.field_samples > report.affine_frames);

    // At 2 substeps this shed moves a wake line 2.755e-2 m in one substep,
    // beyond fs-vpm's core/edge resolution limit of 2.720e-2 m, which refuses
    // it and asks for more substeps. The replay runs at a resolved count,
    // the same for the original and the resumed wind.
    let resolved = SelfInducedShed {
        substeps: 8,
        ..native
    };
    let mut replay = MannWakeWind::resume(wind.checkpoint(), BYTES).unwrap();
    let mut replay_near = near.clone();
    let mut replay_far = far.clone();
    let a = wind
        .shed_hybrid_with_sampling(
            &mut near,
            &mut far,
            &[0.15, 0.19],
            resolved,
            boundary,
            integration,
            sampling,
            budget(),
            || false,
        )
        .unwrap();
    let b = replay
        .shed_hybrid_with_sampling(
            &mut replay_near,
            &mut replay_far,
            &[0.15, 0.19],
            resolved,
            boundary,
            integration,
            sampling,
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(near, replay_near);
    assert_eq!(far, replay_far);
    assert_eq!(wind.checkpoint(), replay.checkpoint());
}

#[test]
fn late_affine_refusal_and_final_cancellation_preserve_the_native_transaction() {
    let original = wake();
    // Discover the native first material query through the original callback;
    // its deliberate refusal must leave the original wake unchanged.
    let mut probe = original.clone();
    let mut anchor = None;
    assert_eq!(
        probe
            .shed_self_induced_with_environment(
                &[0.17, 0.2],
                settings(),
                WakeBoundary::default(),
                |point, _| {
                    anchor = Some(point);
                    Err(error(
                        "capture-first-query",
                        "test observes the native first material query",
                    ))
                },
                || false,
            )
            .unwrap_err()
            .code,
        "capture-first-query"
    );
    assert_eq!(probe, original);
    let anchor = anchor.unwrap();
    let mut near = original.clone();
    let mut wind = MannWakeWind::prepare(source(), &near, atmosphere_config(), BYTES).unwrap();
    let before = wind.checkpoint();
    let integration = MannWakeIntegration::Midpoint { tree: None };
    let zero = MannWakeSampling::Affine {
        anchor_m: anchor,
        accuracy: MannAffineAccuracy {
            maximum_velocity_error_mps: 0.0,
            maximum_gradient_error_per_s: 0.0,
        },
    };
    // The first query is exact at the anchor; a later, distinct vertex exceeds
    // the zero spatial budget. No accepted wind or native geometry may escape.
    assert_eq!(
        wind.shed_near_with_sampling(
            &mut near,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            integration,
            zero,
            budget(),
            || false,
        )
        .unwrap_err()
        .code,
        "mann-affine-accuracy"
    );
    assert_eq!(near, original);
    assert_eq!(wind.checkpoint(), before);

    let admitted = sampling([0.0, 0.0, 1.0]);
    let mut polls = 0;
    let accepted = wind
        .shed_near_with_sampling(
            &mut near,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            integration,
            admitted,
            budget(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut replay_wind = MannWakeWind::resume(before.clone(), BYTES).unwrap();
    let mut replay = original.clone();
    let repeated = replay_wind
        .shed_near_with_sampling(
            &mut replay,
            &[0.17, 0.2],
            settings(),
            WakeBoundary::default(),
            integration,
            admitted,
            budget(),
            || false,
        )
        .unwrap();
    assert_eq!(accepted, repeated);
    assert_eq!(near, replay);
    assert_eq!(wind.checkpoint(), replay_wind.checkpoint());
    let mut cancelled = MannWakeWind::resume(before.clone(), BYTES).unwrap();
    let mut refused = original.clone();
    let mut seen = 0;
    assert!(
        cancelled
            .shed_near_with_sampling(
                &mut refused,
                &[0.17, 0.2],
                settings(),
                WakeBoundary::default(),
                integration,
                admitted,
                budget(),
                || {
                    seen += 1;
                    seen == polls
                },
            )
            .is_err()
    );
    assert_eq!(refused, original);
    assert_eq!(cancelled.checkpoint(), before);
}
