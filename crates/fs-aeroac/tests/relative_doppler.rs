//! Observed-clock admission must retain shared velocity, not assume opposition.

use fs_aeroac::curle3d::moving::flow::UniformFlowFrame;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::SourceKinematics;
use fs_material::fluid_acoustics::FluidAcousticState;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}
fn budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 16,
        maximum_payload_bytes: 262_144,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}
fn knots(times: [f64; 2], motion: impl Fn(f64) -> SourceKinematics) -> [TimedSourceKinematics; 2] {
    times.map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: motion(source_time_s),
    })
}
fn straight(t: f64, position: [f64; 3], velocity: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m: core::array::from_fn(|i| position[i] + velocity[i] * t),
        velocity_m_s: velocity,
        acceleration_m_s2: [0.0; 3],
    }
}
fn path(samples: &[TimedSourceKinematics], maximum_mach: f64) -> PreparedKinematicTrajectory {
    PreparedKinematicTrajectory::new(medium(), samples, maximum_mach, budget(), &mut || false)
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
            maximum_segment_pairs: 256,
            maximum_bound_nodes: 65_536,
            maximum_subdivision_depth: 24,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn check_window(owner: &PreparedMovingReceiverSource<'_>) {
    let upper = owner.maximum_source_time_rate();
    for i in 0..=64 {
        let sample = owner.prepare_receiver(f64::from(i) / 512.0, &mut || false).unwrap();
        assert!(sample.observed_source_time_rate <= upper,
            "rate={} exceeds bound={upper}", sample.observed_source_time_rate);
    }
}

#[test]
fn common_translation_does_not_demand_opposing_motion_bandwidth() {
    for speed in [-240.0, 0.0, 240.0] {
        let velocity = [speed, 0.0, 0.0];
        let source = path(&knots([-0.5, 0.5], |t| straight(t, [0.0; 3], velocity)), 0.95);
        let receiver = path(&knots([0.0, 0.125], |t| straight(t, [12.0, 0.0, 0.0], velocity)), 0.95);
        let owner = pair(&source, &receiver);
        // The old independent bound is about nine at |v|/c=0.8. The actual
        // two-clock map is a translation, irrespective of approach/recession.
        if speed != 0.0 {
            let old = (1.0 + receiver.maximum_mach()) / (1.0 - source.maximum_mach());
            assert!(old > 8.9);
        } else {
            assert_eq!(owner.maximum_source_time_rate(), 1.0);
        }
        assert!(owner.maximum_source_time_rate() < 1.0 + 1e-9);
        check_window(&owner);
        let sample = owner.prepare_receiver(0.0625, &mut || false).unwrap();
        close(sample.observed_source_time_rate, 1.0);
        close(sample.pressure.source_time_s(), 0.0625 - 12.0 / (300.0 - speed));
    }
}

#[test]
fn common_wind_and_laboratory_boost_preserve_the_tight_bound() {
    let mut expected: Option<(f64, f64)> = None;
    for boost in [[0.0; 3], [1024.0, -512.0, 256.0]] {
        let wind = core::array::from_fn(|i| [240.0, -32.0, 16.0][i] + boost[i]);
        let frame = UniformFlowFrame::new(wind, 0.0, [0.0; 3]).unwrap();
        let source_knots = knots([-0.5, 0.5], |t| straight(t, [0.0; 3], boost));
        let receiver_knots = knots([0.0, 0.125], |t| straight(t, [12.0, 4.0, 2.0], boost));
        let source = PreparedKinematicTrajectory::in_uniform_flow(
            medium(), &source_knots, frame, 0.95, budget(), &mut || false,
        ).unwrap();
        let receiver = PreparedKinematicTrajectory::in_uniform_flow(
            medium(), &receiver_knots, frame, 0.95, budget(), &mut || false,
        ).unwrap();
        let owner = pair(&source, &receiver);
        assert!(owner.maximum_source_time_rate() < 1.0 + 1e-9);
        check_window(&owner);
        let sample = owner.prepare_receiver(0.0625, &mut || false).unwrap();
        close(sample.observed_source_time_rate, 1.0);
        if let Some((time, rate)) = expected {
            close(sample.pressure.source_time_s(), time);
            close(owner.maximum_source_time_rate(), rate);
        } else {
            expected = Some((sample.pressure.source_time_s(), owner.maximum_source_time_rate()));
        }
    }
}

#[test]
fn real_opposing_velocity_doppler_is_not_relaxed() {
    let source = path(&knots([-0.5, 0.5], |t| straight(t, [0.0; 3], [120.0, 0.0, 0.0])), 0.8);
    let receiver = path(&knots([0.0, 0.125], |t| straight(t, [60.0, 0.0, 0.0], [-120.0, 0.0, 0.0])), 0.8);
    let owner = pair(&source, &receiver);
    let expected = 420.0 / 180.0;
    assert!(owner.maximum_source_time_rate() >= expected);
    close(owner.maximum_source_time_rate(), expected);
    check_window(&owner);
    for t in [0.0, 0.0625, 0.125] {
        let sample = owner.prepare_receiver(t, &mut || false).unwrap();
        close(sample.pressure.source_time_s(), (420.0 * t - 60.0) / 180.0);
        close(sample.observed_source_time_rate, expected);
    }
}

#[test]
fn equal_simultaneous_accelerating_velocities_are_not_retarded_comotion() {
    let motion = |t: f64, offset: f64| SourceKinematics {
        position_m: [offset + 64.0 * t - 64.0 * t * t, 0.0, 0.0],
        velocity_m_s: [64.0 - 128.0 * t, 0.0, 0.0],
        acceleration_m_s2: [-128.0, 0.0, 0.0],
    };
    let source = path(&knots([-0.25, 0.25], |t| motion(t, 0.0)), 0.6);
    let receiver = path(&knots([0.0, 0.125], |t| motion(t, 12.0)), 0.6);
    let owner = pair(&source, &receiver);
    check_window(&owner);
    for t in [0.0, 0.0625, 0.125] {
        let sample = owner.prepare_receiver(t, &mut || false).unwrap();
        // Independent bisection of the analytic trajectory's sound cone.
        let (mut lo, mut hi) = (-0.25, t);
        for _ in 0..80 {
            let tau = 0.5 * lo + 0.5 * hi;
            let distance = motion(t, 12.0).position_m[0] - motion(tau, 0.0).position_m[0];
            if tau + distance / 300.0 < t { lo = tau; } else { hi = tau; }
        }
        let tau = 0.5 * lo + 0.5 * hi;
        let rate = (300.0 - motion(t, 12.0).velocity_m_s[0])
            / (300.0 - motion(tau, 0.0).velocity_m_s[0]);
        close(sample.pressure.source_time_s(), tau);
        close(sample.observed_source_time_rate, rate);
        assert!(rate > 1.01);
        assert!(owner.maximum_source_time_rate() >= rate);
    }
}

#[test]
fn interior_velocity_extrema_and_admitted_subdivisions_enter_the_bound() {
    // Endpoint velocities agree with the microphone, but the source reaches
    // 124 m/s inside its quintic. Its initial Bernstein velocity control is
    // 224 m/s, above the explicit 135 m/s cap: native subdivision is required.
    // A knot-only bound would incorrectly say one. Keeping the parent hull
    // would lose the accepted subdivision's useful tighter bound.
    let source_samples = [
        TimedSourceKinematics { source_time_s: -0.5, kinematics: straight(0.0, [-32.0, 0.0, 0.0], [64.0, 0.0, 0.0]) },
        TimedSourceKinematics { source_time_s: 0.5, kinematics: straight(0.0, [64.0, 0.0, 0.0], [64.0, 0.0, 0.0]) },
    ];
    let source = path(&source_samples, 0.45);
    let receiver = path(&knots([0.0, 0.125], |t| straight(t, [50.0, 0.0, 0.0], [64.0, 0.0, 0.0])), 0.45);
    let owner = pair(&source, &receiver);
    assert!(owner.maximum_source_time_rate() > 1.25);
    assert!(owner.maximum_source_time_rate() < 1.5);
    check_window(&owner);
    assert!(owner.prepare_receiver(0.0625, &mut || false).unwrap().observed_source_time_rate > 1.25);
    // Observation and refusal do not change the prepared bound or native path.
    let retained = owner.maximum_source_time_rate();
    assert!(owner.prepare_receiver(0.0625, &mut || true).is_err());
    assert!(owner.prepare_receiver(0.126, &mut || false).is_err());
    assert_eq!(owner.maximum_source_time_rate(), retained);
}

#[test]
fn directional_bound_preserves_receding_redshift_on_every_axis() {
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut position = [0.0; 3];
            let mut velocity = [0.0; 3];
            position[axis] = sign * 100.0;
            velocity[axis] = -sign * 120.0;
            let source = path(&knots([-0.5, 0.5], |t| straight(t, [0.0; 3], velocity)), 0.8);
            let receiver = path(&knots([0.0, 0.125], |t| straight(t, position, [0.0; 3])), 0.8);
            let owner = pair(&source, &receiver);
            let expected = 300.0 / 420.0;
            assert!(owner.maximum_source_time_rate() >= expected);
            close(owner.maximum_source_time_rate(), expected);
            assert!(owner.maximum_source_time_rate() < 1.0);
            check_window(&owner);
            for t in [0.0, 0.0625, 0.125] {
                let sample = owner.prepare_receiver(t, &mut || false).unwrap();
                close(sample.observed_source_time_rate, expected);
                close(sample.pressure.source_time_s(), (300.0 * t - 100.0) / 420.0);
            }
        }
    }
}

#[test]
fn moving_microphone_recession_keeps_its_own_clock_numerator() {
    let source = path(&knots([-0.5, 0.5], |t| straight(t, [0.0; 3], [0.0; 3])), 0.8);
    let receiver = path(
        &knots([0.0, 0.125], |t| straight(t, [100.0, 0.0, 0.0], [120.0, 0.0, 0.0])),
        0.8,
    );
    let owner = pair(&source, &receiver);
    close(owner.maximum_source_time_rate(), 0.6);
    check_window(&owner);
    let sample = owner.prepare_receiver(0.0625, &mut || false).unwrap();
    close(sample.observed_source_time_rate, 0.6);
    // Source radiation is still Eulerian; receiver Doppler is not an
    // amplitude multiplier or a replacement for the source's kappa.
    close(sample.pressure.source_time_rate(), 1.0);
}

#[test]
fn unused_future_source_spans_do_not_impose_their_approach_bandwidth() {
    let source_samples = [
        TimedSourceKinematics { source_time_s: -0.5, kinematics: straight(-0.5, [0.0; 3], [-120.0, 0.0, 0.0]) },
        TimedSourceKinematics { source_time_s: 0.0, kinematics: straight(0.0, [0.0; 3], [-120.0, 0.0, 0.0]) },
        TimedSourceKinematics { source_time_s: 0.5, kinematics: straight(0.0, [0.0; 3], [120.0, 0.0, 0.0]) },
    ];
    let source = path(&source_samples, 0.95);
    let receiver = path(
        &knots([0.0, 0.125], |t| straight(t, [100.0, 0.0, 0.0], [0.0; 3])),
        0.8,
    );
    let owner = pair(&source, &receiver);
    close(owner.maximum_source_time_rate(), 300.0 / 420.0);
    check_window(&owner);
}

#[test]
fn dense_histories_admit_only_causally_possible_span_pairs() {
    let trajectory_budget = TrajectoryBudget {
        maximum_segments: 512,
        maximum_payload_bytes: 2_000_000,
        ..budget()
    };
    let source_knots: Vec<_> = (0..=384)
        .map(|i| TimedSourceKinematics {
            source_time_s: -0.5 + f64::from(i) / 256.0,
            kinematics: straight(0.0, [0.0; 3], [0.0; 3]),
        })
        .collect();
    let receiver_knots: Vec<_> = (0..=256)
        .map(|i| TimedSourceKinematics {
            source_time_s: f64::from(i) / 256.0,
            kinematics: straight(0.0, [9.375, 0.0, 0.0], [0.0; 3]),
        })
        .collect();
    let source = PreparedKinematicTrajectory::new(
        medium(), &source_knots, 0.8, trajectory_budget, &mut || false,
    ).unwrap();
    let receiver = PreparedKinematicTrajectory::new(
        medium(), &receiver_knots, 0.8, trajectory_budget, &mut || false,
    ).unwrap();
    let mut polls = 0;
    let owner = PreparedMovingReceiverSource::new(
        &source,
        &receiver,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.0, 1.0],
            minimum_distance_m: 0.01,
            maximum_time_residual_s: 1e-12,
        },
        MovingReceiverBudget {
            // 384 * 256 = 98,304 combinations would fail. At most three
            // original source spans touch each receiver span's sound cone.
            maximum_segment_pairs: 1024,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 24,
            maximum_root_iterations: 96,
        },
        &mut || { polls += 1; false },
    ).unwrap();
    assert!(polls < 12_000, "noncausal history scan: {polls} polls");
    assert_eq!(owner.maximum_source_time_rate(), 1.0);
    for i in 0..=256 {
        let t = f64::from(i) / 256.0;
        let sample = owner.prepare_receiver(t, &mut || false).unwrap();
        close(sample.pressure.source_time_s(), t - 1.0 / 32.0);
        close(sample.observed_source_time_rate, 1.0);
    }
}

#[test]
fn causal_brackets_evaluate_clipped_receiver_windows_inside_original_spans() {
    let trajectory_budget = TrajectoryBudget {
        maximum_segments: 512,
        maximum_payload_bytes: 2_000_000,
        ..budget()
    };
    let source_knots: Vec<_> = (0..=256)
        .map(|i| {
            let t = -0.5 + f64::from(i) / 256.0;
            TimedSourceKinematics {
                source_time_s: t,
                kinematics: straight(t, [0.0; 3], [-120.0, 0.0, 0.0]),
            }
        })
        .collect();
    let source = PreparedKinematicTrajectory::new(
        medium(), &source_knots, 0.8, trajectory_budget, &mut || false,
    ).unwrap();
    let receiver = path(
        &knots([-1.0, 2.0], |t| straight(t, [100.0, 0.0, 0.0], [60.0, 0.0, 0.0])),
        0.8,
    );
    let owner = PreparedMovingReceiverSource::new(
        &source,
        &receiver,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.03125, 0.0625],
            minimum_distance_m: 0.01,
            maximum_time_residual_s: 1e-12,
        },
        MovingReceiverBudget {
            maximum_segment_pairs: 16,
            maximum_bound_nodes: 256,
            maximum_subdivision_depth: 24,
            maximum_root_iterations: 96,
        },
        &mut || false,
    ).unwrap();
    close(owner.maximum_source_time_rate(), 240.0 / 420.0);
    for i in 0..=32 {
        let t = 0.03125 + f64::from(i) / 1024.0;
        let sample = owner.prepare_receiver(t, &mut || false).unwrap();
        close(sample.pressure.source_time_s(), (240.0 * t - 100.0) / 420.0);
        assert!(sample.observed_source_time_rate <= owner.maximum_source_time_rate());
    }
}
