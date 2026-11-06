//! Native histories whose shared velocity previously implied false bandwidth.
//! No sample-rate, bandwidth, root, or source-history gate is bypassed here.

use fs_aeroac::curle3d::moving::flow::UniformFlowFrame;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, PreparedReflectedSourceTrajectory, PreparedSourceTrajectory,
    ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    TimedCompactSource, UniformSourceMotion,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn path(
    first: f64,
    last: f64,
    offset: f64,
    laboratory_speed: f64,
    wind: f64,
) -> PreparedKinematicTrajectory {
    let knots = [first, last].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: SourceKinematics {
            position_m: [offset + laboratory_speed * source_time_s, 0.0, 0.0],
            velocity_m_s: [laboratory_speed, 0.0, 0.0],
            acceleration_m_s2: [0.0; 3],
        },
    });
    PreparedKinematicTrajectory::in_uniform_flow(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        &knots,
        UniformFlowFrame::new([wind, 0.0, 0.0], 0.0, [0.0; 3]).unwrap(),
        0.95,
        TrajectoryBudget {
            maximum_segments: 4,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}
fn pair<'a>(
    source: &'a PreparedKinematicTrajectory,
    receiver: &'a PreparedKinematicTrajectory,
) -> PreparedMovingReceiverSource<'a> {
    PreparedMovingReceiverSource::new(
        source,
        receiver,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.0, 0.125],
            minimum_distance_m: 0.01,
            maximum_time_residual_s: 1e-12,
        },
        MovingReceiverBudget {
            maximum_segment_pairs: 16,
            maximum_bound_nodes: 65_536,
            maximum_subdivision_depth: 24,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.0,
        sample_rate_hz: 512,
        samples: 65,
        max_block: 65,
        receiver_bandwidth_hz: 40.1,
        maximum_time_residual_s: 1e-12,
    }
}
fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 2,
        maximum_source_samples: 2048,
        maximum_output_samples: 4096,
        maximum_payload_bytes: 262_144,
        maximum_source_receiver_pairs: 1_000_000,
    }
}
fn strength(t: f64) -> CompactSourceStrength {
    let omega = core::f64::consts::TAU * 32.0;
    CompactSourceStrength {
        volume_flow_m3_s: 1e-5 * (omega * t).sin(),
        volume_flow_rate_m3_s2: 1e-5 * omega * (omega * t).cos(),
        ..CompactSourceStrength::default()
    }
}
fn history(first: i32, last: i32) -> Vec<TimedCompactSource> {
    (first..=last)
        .map(|i| {
            let source_time_s = f64::from(i) / 512.0;
            TimedCompactSource { source_time_s, strength: strength(source_time_s) }
        })
        .collect()
}

#[test]
fn high_common_speed_renders_original_bandwidth_without_inventing_doppler() {
    // All three configurations have exactly 1/32 s travel time, so every
    // retarded receiver sample lies on an original source-history knot.
    // The first has common motion in still fluid, the second fixed laboratory
    // endpoints in wind, and the third a common laboratory velocity boost.
    for (speed, wind, separation) in [
        (240.0, 0.0, 1.875),
        (0.0, 240.0, 16.875),
        (1024.0, 1264.0, 16.875),
    ] {
        let source = path(-0.0625, 0.1875, 0.0, speed, wind);
        let receiver = path(0.0, 0.125, separation, speed, wind);
        let motion = pair(&source, &receiver);
        let old_bound = (1.0 + receiver.maximum_mach()) / (1.0 - source.maximum_mach());
        assert!(40.0 * old_bound > 256.0); // Would exceed receiver Nyquist.
        assert!(40.0 * motion.maximum_source_time_rate() < config().receiver_bandwidth_hz);
        let accepted = history(-32, 96);
        let sources = [CompactSourceHistory {
            motion: &motion,
            samples: &accepted,
            source_bandwidth_hz: 40.0,
        }];
        let build = || {
            MovingCompactRenderer::with_moving_receiver(&sources, config(), budget(), &mut || false)
                .unwrap()
        };
        let mut renderer = build();
        let mut pressure = [0.0; 65];
        renderer.block(&mut pressure).unwrap();
        // Independent axial moving-monopole formula, including the convective
        // Q term. The microphone's observed rate is one, NOT 1/(1-M_source).
        let m = (speed - wind) / 300.0;
        let kappa = 1.0 - m;
        let r = 300.0 / 32.0;
        for (i, &actual) in pressure.iter().enumerate() {
            let t = renderer.receiver_time_at(i as u64).unwrap();
            let q = strength(t - 1.0 / 32.0);
            let expected = 1.2 / (4.0 * core::f64::consts::PI)
                * (q.volume_flow_rate_m3_s2 / (r * kappa * kappa)
                    + q.volume_flow_m3_s * 300.0 * (m - m * m)
                        / (r * r * kappa * kappa * kappa));
            assert!((actual - expected).abs() < 2e-10,
                "pressure={actual:e}, oracle={expected:e}");
        }
        let mut partitioned = build();
        let mut chunks = [0.0; 65];
        for chunk in chunks.chunks_mut(7) {
            partitioned.block(chunk).unwrap();
        }
        assert_eq!(pressure.map(f64::to_bits), chunks.map(f64::to_bits));

        let mut cancelled = build();
        let mut pending = [123.0; 65];
        let mut polls = 0;
        assert!(cancelled.block_controlled(&mut pending, &mut || {
            polls += 1;
            polls == 12
        }).is_err());
        assert_eq!(pending, [123.0; 65]);
        assert_eq!(cancelled.samples_rendered(), 0);
        cancelled.block(&mut pending).unwrap();
        assert_eq!(pressure.map(f64::to_bits), pending.map(f64::to_bits));
    }
}

#[test]
fn shared_motion_does_not_bypass_source_or_receiver_sampling_limits() {
    let source = path(-0.0625, 0.1875, 0.0, 0.0, 240.0);
    let receiver = path(0.0, 0.125, 16.875, 0.0, 240.0);
    let motion = pair(&source, &receiver);
    let accepted = history(-32, 96);
    let sources = [CompactSourceHistory {
        motion: &motion, samples: &accepted, source_bandwidth_hz: 40.0,
    }];
    for bandwidth in [39.9, 256.0] {
        let config = MovingPressureConfig { receiver_bandwidth_hz: bandwidth, ..config() };
        assert!(MovingCompactRenderer::with_moving_receiver(&sources, config, budget(), &mut || false).is_err());
    }
    let sparse = [accepted[0], accepted[64], accepted[128]];
    let unresolved = [CompactSourceHistory {
        motion: &motion, samples: &sparse, source_bandwidth_hz: 40.0,
    }];
    assert!(MovingCompactRenderer::with_moving_receiver(&unresolved, config(), budget(), &mut || false).is_err());
    assert!(MovingCompactRenderer::with_moving_receiver(&sources, config(), budget(), &mut || true).is_err());
}

#[test]
fn opposing_motion_still_requires_its_actual_doppler_bandwidth() {
    let source = path(-0.5, 0.5, 0.0, 120.0, 0.0);
    let receiver = path(0.0, 0.125, 60.0, -120.0, 0.0);
    let motion = pair(&source, &receiver);
    let accepted = history(-256, 256);
    let sources = [CompactSourceHistory {
        motion: &motion, samples: &accepted, source_bandwidth_hz: 40.0,
    }];
    assert!(motion.maximum_source_time_rate() >= 7.0 / 3.0);
    assert!(MovingCompactRenderer::with_moving_receiver(&sources, config(), budget(), &mut || false).is_err());
    let resolved = MovingPressureConfig { receiver_bandwidth_hz: 100.0, ..config() };
    assert!(MovingCompactRenderer::with_moving_receiver(&sources, resolved, budget(), &mut || false).is_ok());
}

#[test]
fn receding_native_history_renders_at_its_contracted_receiver_bandwidth() {
    let source = path(-0.5, 0.5, 0.0, -120.0, 0.0);
    let receiver = path(0.0, 0.125, 100.0, 0.0, 0.0);
    let motion = pair(&source, &receiver);
    let accepted = history(-256, 256);
    let sources = [CompactSourceHistory {
        motion: &motion,
        samples: &accepted,
        source_bandwidth_hz: 40.0,
    }];
    let config = MovingPressureConfig {
        receiver_bandwidth_hz: 29.0,
        ..config()
    };
    assert!(40.0 * motion.maximum_source_time_rate() < config.receiver_bandwidth_hz);
    let build = || {
        MovingCompactRenderer::with_moving_receiver(&sources, config, budget(), &mut || false)
            .unwrap()
    };
    let mut renderer = build();
    let mut pressure = [0.0; 65];
    renderer.block(&mut pressure).unwrap();
    // Independent axial pressure including convection/induction. The accepted
    // Hermite history converges to this analytic tone, it is not a pitch filter.
    let mut error_square = 0.0;
    let mut reference_square = 0.0;
    for (i, &actual) in pressure.iter().enumerate() {
        let t = i as f64 / 512.0;
        let tau = (300.0 * t - 100.0) / 420.0;
        let r = 100.0 + 120.0 * tau;
        let q = strength(tau);
        let m = -0.4_f64;
        let kappa = 1.0 - m;
        let expected = 1.2 / (4.0 * core::f64::consts::PI)
            * (q.volume_flow_rate_m3_s2 / (r * kappa * kappa)
                + q.volume_flow_m3_s * 300.0 * (m - m * m)
                    / (r * r * kappa * kappa * kappa));
        error_square += (actual - expected).powi(2);
        reference_square += expected * expected;
    }
    assert!((error_square / reference_square).sqrt() < 1e-3);
    let mut partitioned = build();
    let mut chunks = [0.0; 65];
    for chunk in chunks.chunks_mut(7) {
        partitioned.block(chunk).unwrap();
    }
    assert_eq!(pressure.map(f64::to_bits), chunks.map(f64::to_bits));
    // A contracted clock still needs enough receiver bandwidth for its tone.
    let unresolved = MovingPressureConfig {
        receiver_bandwidth_hz: 28.0,
        ..config
    };
    assert!(MovingCompactRenderer::with_moving_receiver(
        &sources, unresolved, budget(), &mut || false,
    ).is_err());
}

fn fixed_path_with_later_approach(
    observer: MovingCompactObserver,
    source_y: f64,
) -> PreparedSourceTrajectory {
    // The first span recedes at 120 m/s. The next turns smoothly and ends
    // approaching at 120 m/s, but none of its emissions reaches this window.
    let knots = [(-0.5, 60.0, -120.0), (0.0, 0.0, -120.0), (1.0, 0.0, 120.0)].map(
        |(source_time_s, x, speed)| TimedSourceKinematics {
            source_time_s,
            kinematics: SourceKinematics {
                position_m: [x, source_y, 0.0],
                velocity_m_s: [speed, 0.0, 0.0],
                acceleration_m_s2: [0.0; 3],
            },
        },
    );
    PreparedSourceTrajectory::new(
        observer,
        &knots,
        TrajectoryBudget {
            maximum_segments: 4,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

#[test]
fn fixed_curved_history_admits_only_the_requested_receding_window() {
    let observer = MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        [100.0, 0.0, 0.0],
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let motion = fixed_path_with_later_approach(observer, 0.0);
    assert!(motion.maximum_source_time_rate() > 1.0);
    let accepted = history(-256, 512);
    let sources = [CompactSourceHistory {
        motion: &motion,
        samples: &accepted,
        source_bandwidth_hz: 40.0,
    }];
    let config = MovingPressureConfig {
        receiver_bandwidth_hz: 29.0,
        ..config()
    };
    let build =
        || MovingCompactRenderer::new(observer, &sources, config, budget(), &mut || false).unwrap();
    let mut renderer = build();
    let mut actual = [0.0; 65];
    renderer.block(&mut actual).unwrap();

    // An independent analytic retarded inversion of the same straight first
    // span must produce the same physical pressure; later motion is irrelevant.
    let uniform = UniformSourceMotion::new(0.0, [0.0; 3], [-120.0, 0.0, 0.0]).unwrap();
    let reference_sources = [CompactSourceHistory {
        motion: uniform,
        samples: &accepted,
        source_bandwidth_hz: 40.0,
    }];
    let mut reference =
        MovingCompactRenderer::new(observer, &reference_sources, config, budget(), &mut || {
            false
        })
        .unwrap();
    let mut expected = [0.0; 65];
    reference.block(&mut expected).unwrap();
    for (actual, expected) in actual.iter().zip(expected) {
        assert!(
            (actual - expected).abs() < 1e-10,
            "{actual:e} vs {expected:e}"
        );
    }
    let mut partitioned = build();
    let mut chunks = [0.0; 65];
    for chunk in chunks.chunks_mut(7) {
        partitioned.block(chunk).unwrap();
    }
    assert_eq!(actual.map(f64::to_bits), chunks.map(f64::to_bits));
    let insufficient = MovingPressureConfig {
        receiver_bandwidth_hz: 28.0,
        ..config
    };
    assert!(
        MovingCompactRenderer::new(observer, &sources, insufficient, budget(), &mut || false,)
            .is_err()
    );
}

#[test]
fn fixed_reflected_history_keeps_each_receding_delay_and_wall_sign() {
    let medium = FluidAcousticState::try_new(1.2, 300.0).unwrap();
    let limits = MovingCompactLimits {
        maximum_mach: 0.8,
        minimum_distance_m: 0.01,
    };
    let observer = MovingCompactObserver::new(medium, [100.0, 10.0, 0.0], limits).unwrap();
    let trajectory = fixed_path_with_later_approach(observer, 10.0);
    let accepted = history(-256, 512);
    let uniform = UniformSourceMotion::new(0.0, [0.0, 10.0, 0.0], [-120.0, 0.0, 0.0]).unwrap();
    let reference_sources = [CompactSourceHistory {
        motion: uniform,
        samples: &accepted,
        source_bandwidth_hz: 40.0,
    }];
    let config = MovingPressureConfig {
        receiver_bandwidth_hz: 31.0,
        ..config()
    };
    for (boundary, sign) in [
        (IdealPlaneBoundary::Rigid, 1.0),
        (IdealPlaneBoundary::PressureRelease, -1.0),
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], [0.0, 1.0, 0.0], boundary).unwrap();
        let images = IdealReflectionObserver::new(observer, &[plane], 2).unwrap();
        let reflected = PreparedReflectedSourceTrajectory::new(
            &trajectory,
            &images,
            ReflectedTrajectoryBudget {
                maximum_payload_bytes: 131_072,
                maximum_bound_nodes: 4096,
                maximum_subdivision_depth: 20,
            },
            &mut || false,
        )
        .unwrap();
        let sources = [CompactSourceHistory {
            motion: &reflected,
            samples: &accepted,
            source_bandwidth_hz: 40.0,
        }];
        let mut renderer =
            MovingCompactRenderer::with_ideal_reflections(&sources, config, budget(), &mut || {
                false
            })
            .unwrap();
        let mut actual = [0.0; 65];
        renderer.block(&mut actual).unwrap();

        let virtual_observer =
            MovingCompactObserver::new(medium, [100.0, -10.0, 0.0], limits).unwrap();
        let mut direct =
            MovingCompactRenderer::new(observer, &reference_sources, config, budget(), &mut || {
                false
            })
            .unwrap();
        let mut image = MovingCompactRenderer::new(
            virtual_observer,
            &reference_sources,
            config,
            budget(),
            &mut || false,
        )
        .unwrap();
        let mut direct_pressure = [0.0; 65];
        let mut image_pressure = [0.0; 65];
        direct.block(&mut direct_pressure).unwrap();
        image.block(&mut image_pressure).unwrap();
        for i in 0..65 {
            let expected = direct_pressure[i] + sign * image_pressure[i];
            assert!(
                (actual[i] - expected).abs() < 1e-10,
                "image sign {sign}: {} vs {expected}",
                actual[i]
            );
        }
    }
}
