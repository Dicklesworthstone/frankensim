use super::*;
use fs_aeroac::curle3d::moving::SourceKinematics;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, MovingReceiverBudget, MovingReceiverConfig,
    PreparedKinematicTrajectory, PreparedMovingReceiverSource, PreparedReflectedSourceTrajectory,
    PreparedSourceTrajectory, ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};

fn trajectory_budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 8,
        maximum_payload_bytes: 131_072,
        maximum_bound_nodes: 4096,
        maximum_subdivision_depth: 16,
        maximum_root_iterations: 64,
    }
}
fn knots(motion: impl Fn(f64) -> SourceKinematics) -> [TimedSourceKinematics; 3] {
    [-0.25, 0.1, 0.5].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: motion(source_time_s),
    })
}
fn curved(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [0.2 + 10.0 * t + 2.0 * t * t, 0.5 + 3.0 * t - t * t, 0.2],
        velocity_m_s: [10.0 + 4.0 * t, 3.0 - 2.0 * t, 0.0],
        acceleration_m_s2: [4.0, -2.0, 0.0],
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
// Original analytic curve and independent bisection, without the prepared
// quintics, retarded solver, gradient coefficients or strength interpolation.
fn potential(x: [f64; 3], time: f64) -> f64 {
    let mut lo: f64 = -0.25;
    let mut hi: f64 = 0.5;
    for _ in 0..100 {
        let tau = lo.midpoint(hi);
        let source = curved(tau);
        let r: [f64; 3] = core::array::from_fn(|i| x[i] - source.position_m[i]);
        if tau + dot(r, r).sqrt() / 300.0 < time {
            lo = tau;
        } else {
            hi = tau;
        }
    }
    let tau = lo.midpoint(hi);
    let source = curved(tau);
    let r: [f64; 3] = core::array::from_fn(|i| x[i] - source.position_m[i]);
    let radius = dot(r, r).sqrt();
    let kappa = 1.0 - dot(r.map(|r| r / radius), source.velocity_m_s) / 300.0;
    -strength(tau, 1.0).volume_flow_m3_s / (4.0 * PI * radius * kappa)
}
fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
}

#[test]
fn curved_near_fields_match_independent_retarded_potential_space_and_time_derivatives() {
    let position = [12.0, 2.0, 1.0];
    let observer = observer(position);
    let trajectory =
        PreparedSourceTrajectory::new(observer, &knots(curved), trajectory_budget(), &mut || false)
            .unwrap();
    let samples = history(1.0);
    let sources = [CompactSourceHistory {
        motion: &trajectory,
        samples: &samples,
        source_bandwidth_hz: 1.0,
    }];
    let mut renderer = MovingMonopoleRenderer::new(
        observer,
        &sources,
        FirstOrder::default(),
        config(),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut fields = [CompactFieldSample::default(); 81];
    renderer.field_block(&mut fields).unwrap();
    for field in fields.iter().step_by(10) {
        let time = field.receiver_time_s;
        let pressure = -1.2 * derivative(|dt| potential(position, time + dt), 1e-5);
        close(field.pressure_pa, pressure, 2e-8);
        for axis in 0..3 {
            let velocity = derivative(
                |dx| {
                    let mut shifted = position;
                    shifted[axis] += dx;
                    potential(shifted, time)
                },
                0.001,
            );
            close(field.particle_velocity_m_s[axis], velocity, 2e-8);
            close(
                field.instantaneous_intensity_w_m2[axis],
                pressure * velocity,
                4e-8,
            );
        }
    }
}

#[test]
fn moving_microphone_changes_emission_clock_without_rescaling_particle_velocity() {
    let source_speed = 60.0;
    let source_path = PreparedKinematicTrajectory::new(
        medium(),
        &knots(|t| SourceKinematics {
            position_m: [source_speed * t, 0.0, 0.0],
            velocity_m_s: [source_speed, 0.0, 0.0],
            acceleration_m_s2: [0.0; 3],
        }),
        0.3,
        trajectory_budget(),
        &mut || false,
    )
    .unwrap();
    let samples = history(1.0);
    for receiver_speed in [-20.0, 0.0, 40.0] {
        let receiver = PreparedKinematicTrajectory::new(
            medium(),
            &knots(|t| SourceKinematics {
                position_m: [12.0 + receiver_speed * t, 0.0, 0.0],
                velocity_m_s: [receiver_speed, 0.0, 0.0],
                acceleration_m_s2: [0.0; 3],
            }),
            0.3,
            trajectory_budget(),
            &mut || false,
        )
        .unwrap();
        let pair = PreparedMovingReceiverSource::new(
            &source_path,
            &receiver,
            MovingReceiverConfig {
                receiver_time_interval_s: [0.05, 0.13],
                minimum_distance_m: 0.1,
                maximum_time_residual_s: 1e-12,
            },
            MovingReceiverBudget {
                maximum_segment_pairs: 64,
                maximum_bound_nodes: 4096,
                maximum_subdivision_depth: 16,
                maximum_root_iterations: 64,
            },
            &mut || false,
        )
        .unwrap();
        let sources = [CompactSourceHistory {
            motion: &pair,
            samples: &samples,
            source_bandwidth_hz: 1.0,
        }];
        let mut renderer = MovingMonopoleRenderer::new(
            &receiver,
            &sources,
            FirstOrder::default(),
            config(),
            budget(),
            &mut || false,
        )
        .unwrap();
        let mut fields = [CompactFieldSample::default(); 81];
        renderer.field_block(&mut fields).unwrap();
        for field in fields {
            let (pressure, velocity) =
                collinear(field.receiver_time_s, source_speed, receiver_speed);
            close(field.pressure_pa, pressure, 2e-9);
            close(field.particle_velocity_m_s[0], velocity, 2e-9);
            close(
                field.instantaneous_intensity_w_m2[0],
                pressure * velocity,
                4e-9,
            );
        }
        assert_eq!(renderer.maximum_retarded_iterations(), 81 * 64);
    }
}

#[test]
fn moving_reflected_histories_preserve_rigid_and_pressure_release_vector_boundary_conditions() {
    let physical = observer([0.0, 3.0, 1.0]);
    let trajectory = PreparedSourceTrajectory::new(
        physical,
        &knots(|t| SourceKinematics {
            position_m: [2.0, 1.0 + 2.0 * t, 0.5],
            velocity_m_s: [0.0, 2.0, 0.0],
            acceleration_m_s2: [0.0; 3],
        }),
        trajectory_budget(),
        &mut || false,
    )
    .unwrap();
    let samples = history(1.0);
    let normal_pickup = FirstOrder::new(0.0, [1.0, 0.0, 0.0]).unwrap();
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let walls = IdealReflectionObserver::new(
            physical,
            &[IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], boundary).unwrap()],
            2,
        )
        .unwrap();
        let reflected = PreparedReflectedSourceTrajectory::new(
            &trajectory,
            &walls,
            ReflectedTrajectoryBudget {
                maximum_bound_nodes: 4096,
                maximum_subdivision_depth: 16,
            },
            &mut || false,
        )
        .unwrap();
        let sources = [CompactSourceHistory {
            motion: &reflected,
            samples: &samples,
            source_bandwidth_hz: 1.0,
        }];
        let mut renderer = MovingMonopoleRenderer::new(
            &walls,
            &sources,
            normal_pickup,
            config(),
            budget(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(renderer.source_count(), 1);
        assert_eq!(renderer.source_path_count(), 2);
        assert_eq!(renderer.maximum_retarded_iterations(), 2 * 81 * 64);
        let mut fields = [CompactFieldSample::default(); 81];
        renderer.field_block(&mut fields).unwrap();
        for field in fields {
            match boundary {
                IdealPlaneBoundary::Rigid => {
                    assert!(field.pressure_pa.abs() > 1e-8);
                    assert_eq!(field.particle_velocity_m_s[0], 0.0);
                    assert_eq!(field.microphone_pressure_pa, 0.0);
                    assert_eq!(field.instantaneous_intensity_w_m2[0], 0.0);
                }
                IdealPlaneBoundary::PressureRelease => {
                    assert_eq!(field.pressure_pa, 0.0);
                    assert_eq!(&field.particle_velocity_m_s[1..], &[0.0; 2]);
                    assert!(field.particle_velocity_m_s[0].abs() > 1e-8);
                    close(
                        field.microphone_pressure_pa,
                        -360.0 * field.particle_velocity_m_s[0],
                        2e-12,
                    );
                    assert_eq!(field.instantaneous_intensity_w_m2, [0.0; 3]);
                }
            }
        }
    }
}
