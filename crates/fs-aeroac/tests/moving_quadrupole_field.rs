//! Full moving-stress fields against stationary fields, force derivatives and momentum.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedMovingReceiverSource,
    PreparedReflectedMovingReceiverSource, PreparedReflectedSourceTrajectory,
    PreparedSourceTrajectory, ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    UniformSourceMotion,
};
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, RetardedQuadrupoleField, TransientQuadrupoleSource,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;

const TENSOR: [f64; 6] = [0.7, -0.2, 0.4, 0.3, -0.1, 0.23];
const QUARTIC: [f64; 6] = [-0.11, 0.13, 0.09, -0.05, 0.07, 0.02];
const C: f64 = 4.0;
const RHO: f64 = 1.2;
const START: f64 = -1.0;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(RHO, C).unwrap()
}

fn observer(x: [f64; 3], medium: FluidAcousticState) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium,
        x,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 1e-5,
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, scale: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale.abs().max(1e-25),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}, tolerance={relative:e}"
    );
}

fn norm(v: [f64; 3]) -> f64 {
    v.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn compare_fields(
    actual: RetardedQuadrupoleField,
    expected: RetardedQuadrupoleField,
    tolerance: f64,
) {
    close(
        actual.pressure_pa,
        expected.pressure_pa,
        expected.pressure_pa,
        tolerance,
    );
    let scale = norm(expected.particle_velocity_m_s);
    for i in 0..3 {
        close(
            actual.particle_velocity_m_s[i],
            expected.particle_velocity_m_s[i],
            scale,
            tolerance,
        );
        close(
            actual.instantaneous_intensity_w_m2[i],
            actual.pressure_pa * actual.particle_velocity_m_s[i],
            actual.pressure_pa.abs() * scale,
            5e-15,
        );
    }
}

fn tensor(q: [f64; 6]) -> [[f64; 3]; 3] {
    [[q[0], q[3], q[4]], [q[3], q[1], q[5]], [q[4], q[5], q[2]]]
}

fn strength(time: f64) -> TransientQuadrupoleSource {
    let u = time - START;
    if u <= 0.0 {
        return TransientQuadrupoleSource::default();
    }
    TransientQuadrupoleSource {
        stress_integral_n_m_s: core::array::from_fn(|i| {
            TENSOR[i] * u.powi(4) / 4.0 + QUARTIC[i] * u.powi(5) / 5.0
        }),
        stress_moment_n_m: core::array::from_fn(|i| TENSOR[i] * u.powi(3) + QUARTIC[i] * u.powi(4)),
        stress_rate_n_m_s: core::array::from_fn(|i| {
            3.0 * TENSOR[i] * u.powi(2) + 4.0 * QUARTIC[i] * u.powi(3)
        }),
        stress_acceleration_n_m_s2: core::array::from_fn(|i| {
            6.0 * TENSOR[i] * u + 12.0 * QUARTIC[i] * u.powi(2)
        }),
    }
}

#[test]
fn stationary_limit_retains_all_four_velocity_orders_without_changing_pressure() {
    for c in [317.0, 1490.0] {
        let medium = FluidAcousticState::try_new(RHO, c).unwrap();
        for range in [0.003, 0.2, 12.0] {
            for n in [[1.0, 0.0, 0.0], [0.6, 0.0, 0.8], [-0.48, 0.64, 0.6]] {
                let x = n.map(|v| range * v);
                let obs = observer(x, medium);
                let motion = UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap();
                let full = obs
                    .prepare_uniform_quadrupole_field(0.137, &motion, 1e-12)
                    .unwrap();
                let original = obs
                    .prepare_uniform_quadrupole(0.137, &motion, 1e-12)
                    .unwrap();
                let fixed =
                    PreparedTransientQuadrupoleObserver::new(medium, x, &[[0.0; 3]], 1).unwrap();
                assert_eq!(full.source_time_s(), original.source_time_s());
                for order in 0..4 {
                    let mut q = TransientQuadrupoleSource::default();
                    match order {
                        0 => q.stress_integral_n_m_s = TENSOR,
                        1 => q.stress_moment_n_m = TENSOR,
                        2 => q.stress_rate_n_m_s = TENSOR,
                        _ => q.stress_acceleration_n_m_s2 = TENSOR,
                    }
                    // Geometry is constant, so its kernel applies to the
                    // actual tensor integral directly in this stationary limit.
                    let history = obs
                        .quadrupole_history_integrand([0.0; 3], q.stress_integral_n_m_s)
                        .unwrap();
                    let actual = full.evaluate(q, history).unwrap();
                    assert_eq!(actual.pressure_pa, original.pressure_pa(q).unwrap());
                    compare_fields(actual, fixed.evaluate_at_retarded(&[q]).unwrap(), 2e-12);
                }
            }
        }
    }
}

#[test]
fn spatial_history_kernel_has_exact_isotropic_null_and_independent_force_derivative() {
    let x = [1.1, -0.8, 1.7];
    let source = [0.2, -0.1, 0.3];
    let obs = observer(x, medium());
    assert_eq!(
        obs.quadrupole_history_integrand(source, [1e15, 1e15, 1e15, 0.0, 0.0, 0.0])
            .unwrap(),
        [0.0; 3]
    );
    let expected = obs.quadrupole_history_integrand(source, TENSOR).unwrap();
    let q = tensor(TENSOR);
    let errors = [0.04, 0.02, 0.01].map(|h| {
        let mut derivative = [0.0; 3];
        for j in 0..3 {
            let column = core::array::from_fn(|i| q[i][j]);
            for (shift, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let mut at = x;
                at[j] += shift * h;
                let force = observer(at, medium())
                    .force_history_integrand(source, column)
                    .unwrap();
                for i in 0..3 {
                    derivative[i] -= weight * force[i] / (12.0 * h);
                }
            }
        }
        norm(core::array::from_fn(|i| derivative[i] - expected[i]))
    });
    assert!(
        errors[1] < 0.1 * errors[0] && errors[2] < 0.1 * errors[1],
        "{errors:?}"
    );
    assert!(
        errors[2] < 2e-7 * norm(expected),
        "{errors:?}, expected={expected:?}"
    );
    let doubled = core::array::from_fn(|i| source[i] + 2.0 * (x[i] - source[i]));
    let distant = observer(doubled, medium())
        .quadrupole_history_integrand(source, TENSOR)
        .unwrap();
    for i in 0..3 {
        close(16.0 * distant[i], expected[i], norm(expected), 3e-14);
    }
}

#[derive(Clone, Copy)]
struct Motion {
    acceleration: [f64; 3],
    jerk: [f64; 3],
}

impl Motion {
    fn at(self, t: f64) -> SourceKinematics {
        let y = [-0.2, 0.1, 0.3];
        let v = [0.65, -0.35, 0.25];
        SourceKinematics {
            position_m: core::array::from_fn(|i| {
                y[i] + v[i] * t
                    + self.acceleration[i] * t * t / 2.0
                    + self.jerk[i] * t * t * t / 6.0
            }),
            velocity_m_s: core::array::from_fn(|i| {
                v[i] + self.acceleration[i] * t + self.jerk[i] * t * t / 2.0
            }),
            acceleration_m_s2: core::array::from_fn(|i| self.acceleration[i] + self.jerk[i] * t),
        }
    }
}

fn motions() -> [Motion; 2] {
    [
        Motion {
            acceleration: [0.0; 3],
            jerk: [0.0; 3],
        },
        Motion {
            acceleration: [0.6, 0.2, -0.3],
            jerk: [-0.4, 0.6, 0.2],
        },
    ]
}

// Independent light-cone bisection; no native inverse or retarded derivative.
fn emission(x: [f64; 3], time: f64, motion: Motion) -> f64 {
    let (mut lo, mut hi) = (START, time);
    for _ in 0..100 {
        let tau = 0.5 * (lo + hi);
        let source = motion.at(tau);
        let range = norm(core::array::from_fn(|i| x[i] - source.position_m[i]));
        if tau + range / C > time {
            hi = tau;
        } else {
            lo = tau;
        }
    }
    0.5 * (lo + hi)
}

// Fixed high-resolution Simpson integration is only an independent test
// oracle. Production history integration belongs to the accepted-history owner.
fn integral(end: f64, mut value: impl FnMut(f64) -> [f64; 3]) -> [f64; 3] {
    const STEPS: usize = 1024;
    let h = (end - START) / STEPS as f64;
    let mut sums = [Dd::ZERO; 3];
    for i in 0..=STEPS {
        let weight = if i == 0 || i == STEPS {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let sample = value(START + i as f64 * h);
        for axis in 0..3 {
            sums[axis] = sums[axis] + Dd::from_f64(weight) * Dd::from_f64(sample[axis]);
        }
    }
    sums.map(|sum| (sum * Dd::from_f64(h / 3.0)).to_f64())
}

fn moving_field(x: [f64; 3], time: f64, motion: Motion) -> RetardedQuadrupoleField {
    let tau = emission(x, time, motion);
    let obs = observer(x, medium());
    let sample = obs
        .prepare_emission_quadrupole_field(tau, motion.at(tau), motion.jerk)
        .unwrap();
    close(sample.propagation().receiver_time_s(), time, time, 2e-14);
    let history = integral(tau, |s| {
        obs.quadrupole_history_integrand(motion.at(s).position_m, strength(s).stress_moment_n_m)
            .unwrap()
    });
    sample.evaluate(strength(tau), history).unwrap()
}

fn force_column_field(x: [f64; 3], time: f64, motion: Motion, column: usize) -> [f64; 3] {
    let tau = emission(x, time, motion);
    let obs = observer(x, medium());
    let sample = obs.prepare_force_emission(tau, motion.at(tau)).unwrap();
    let q = strength(tau);
    let force = tensor(q.stress_moment_n_m);
    let rate = tensor(q.stress_rate_n_m_s);
    let history = integral(tau, |s| {
        let force = tensor(strength(s).stress_moment_n_m);
        obs.force_history_integrand(
            motion.at(s).position_m,
            core::array::from_fn(|i| force[i][column]),
        )
        .unwrap()
    });
    sample
        .evaluate(
            CompactSourceStrength {
                force_on_fluid_n: core::array::from_fn(|i| force[i][column]),
                force_rate_on_fluid_n_s: core::array::from_fn(|i| rate[i][column]),
                ..CompactSourceStrength::default()
            },
            history,
        )
        .unwrap()
        .particle_velocity_m_s
}

#[test]
fn translating_and_accelerating_velocity_converges_to_independent_force_field_divergence() {
    let x = [2.0, -1.3, 1.8];
    for motion in motions() {
        let expected = moving_field(x, 0.7, motion).particle_velocity_m_s;
        let errors = [0.04, 0.02, 0.01].map(|h| {
            let mut velocity = [0.0; 3];
            for column in 0..3 {
                for (shift, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                    let mut at = x;
                    at[column] += shift * h;
                    let force = force_column_field(at, 0.7, motion, column);
                    for i in 0..3 {
                        velocity[i] -= weight * force[i] / (12.0 * h);
                    }
                }
            }
            norm(core::array::from_fn(|i| velocity[i] - expected[i]))
        });
        assert!(
            errors[1] < 0.12 * errors[0] && errors[2] < 0.12 * errors[1],
            "{errors:?}"
        );
        assert!(
            errors[2] < 2e-7 * norm(expected),
            "{errors:?}, velocity={expected:?}"
        );
    }
}

fn moving_pressure(x: [f64; 3], time: f64, motion: Motion) -> f64 {
    let tau = emission(x, time, motion);
    observer(x, medium())
        .prepare_emission_quadrupole(tau, motion.at(tau), motion.jerk)
        .unwrap()
        .pressure_pa(strength(tau))
        .unwrap()
}

#[test]
fn full_fields_close_all_three_source_free_linear_momentum_equations() {
    let x = [2.0, -1.3, 1.8];
    for motion in motions() {
        let errors = [0.02, 0.01, 0.005].map(|h| {
            let mut acceleration = [0.0; 3];
            for (shift, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                let field = moving_field(x, 0.7 + shift * h, motion);
                for i in 0..3 {
                    acceleration[i] += RHO * weight * field.particle_velocity_m_s[i] / (12.0 * h);
                }
            }
            let gradient: [f64; 3] = core::array::from_fn(|axis| {
                let mut sum = 0.0;
                for (shift, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                    let mut at = x;
                    at[axis] += shift * h;
                    sum += weight * moving_pressure(at, 0.7, motion) / (12.0 * h);
                }
                sum
            });
            norm(core::array::from_fn(|i| acceleration[i] + gradient[i])) / norm(gradient)
        });
        assert!(
            errors[1] < 0.15 * errors[0] && errors[2] < 0.15 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 2e-8, "{errors:?}");
    }
}

fn trajectory_budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 8,
        maximum_payload_bytes: 131_072,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}

fn reflection_budget() -> ReflectedTrajectoryBudget {
    ReflectedTrajectoryBudget {
        maximum_payload_bytes: 1_000_000,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
    }
}

fn receiver_budget() -> MovingReceiverBudget {
    MovingReceiverBudget {
        maximum_segment_pairs: 1024,
        maximum_bound_nodes: 65_536,
        maximum_subdivision_depth: 24,
        maximum_root_iterations: 96,
    }
}

fn knots(motion: impl Fn(f64) -> SourceKinematics) -> Vec<TimedSourceKinematics> {
    [-2.0, 0.0, 2.0]
        .into_iter()
        .map(|source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect()
}

fn fixed(position_m: [f64; 3]) -> SourceKinematics {
    SourceKinematics {
        position_m,
        ..SourceKinematics::default()
    }
}

fn walls() -> [IdealAcousticPlane; 2] {
    [
        IdealAcousticPlane::new([0.0; 3], [0.6, 0.8, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
        IdealAcousticPlane::new(
            [0.0; 3],
            [-0.8, 0.6, 0.0],
            IdealPlaneBoundary::PressureRelease,
        )
        .unwrap(),
    ]
}

fn compare_sum(pressure: Dd, velocity: [Dd; 3], expected: RetardedQuadrupoleField) {
    close(
        pressure.to_f64(),
        expected.pressure_pa,
        expected.pressure_pa,
        5e-10,
    );
    for i in 0..3 {
        close(
            velocity[i].to_f64(),
            expected.particle_velocity_m_s[i],
            norm(expected.particle_velocity_m_s),
            5e-10,
        );
    }
}

#[test]
fn native_source_trajectory_and_oblique_images_preserve_full_retarded_fields() {
    let x = [1.0, 3.0, 2.0];
    let source = [0.2, 0.8, 0.5];
    let obs = observer(x, medium());
    let trajectory = PreparedSourceTrajectory::new(
        obs,
        &knots(|_| fixed(source)),
        trajectory_budget(),
        &mut || false,
    )
    .unwrap();
    let planes = walls();
    let images = IdealReflectionObserver::new(obs, &planes, 4).unwrap();
    let reflected = PreparedReflectedSourceTrajectory::new(
        &trajectory,
        &images,
        reflection_budget(),
        &mut || false,
    )
    .unwrap();
    let stationary = PreparedTransientQuadrupoleObserver::new(medium(), x, &[source], 1)
        .unwrap()
        .with_orthogonal_planes(&planes, 4)
        .unwrap();
    for time in [0.4, 0.7, 1.0] {
        let direct = trajectory
            .prepare_quadrupole_field_receiver(time, 1e-13, &mut || false)
            .unwrap();
        let original = trajectory
            .prepare_quadrupole_receiver(time, 1e-13, &mut || false)
            .unwrap();
        assert_eq!(
            direct
                .pressure()
                .pressure_pa(strength(direct.source_time_s()))
                .unwrap(),
            original
                .pressure_pa(strength(original.source_time_s()))
                .unwrap()
        );
        let mut pressure = Dd::ZERO;
        let mut velocity = [Dd::ZERO; 3];
        for path in 0..images.path_count() {
            let full = reflected
                .prepare_quadrupole_field_path(path, time, 1e-13, &mut || false)
                .unwrap();
            let (virtual_receiver, _, _) = images.stationary_field_path(path, source).unwrap();
            let q = strength(full.source_time_s());
            let history = virtual_receiver
                .quadrupole_history_integrand(source, q.stress_integral_n_m_s)
                .unwrap();
            let (p, u) = full.evaluate_dd(q, history.map(Dd::from_f64)).unwrap();
            pressure = pressure + p;
            for i in 0..3 {
                velocity[i] = velocity[i] + u[i];
            }
        }
        compare_sum(
            pressure,
            velocity,
            stationary
                .evaluate_at_time(time, |_, t| Ok(strength(t)))
                .unwrap(),
        );
    }
    assert!(matches!(
        trajectory.prepare_quadrupole_field_receiver(0.7, 1e-13, &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        reflected
            .prepare_quadrupole_field_path(4, 0.7, 1e-13, &mut || false)
            .is_err()
    );
}

#[test]
fn moving_microphone_and_its_images_use_the_current_spatial_history_kernel() {
    let source = [0.2, 0.8, 0.5];
    let receiver = |t: f64| SourceKinematics {
        position_m: [1.0 + 0.1 * t, 3.0 - 0.03 * t, 2.0 + 0.04 * t],
        velocity_m_s: [0.1, -0.03, 0.04],
        acceleration_m_s2: [0.0; 3],
    };
    let source_path = PreparedKinematicTrajectory::new(
        medium(),
        &knots(|_| fixed(source)),
        0.8,
        trajectory_budget(),
        &mut || false,
    )
    .unwrap();
    let receiver_path = PreparedKinematicTrajectory::new(
        medium(),
        &knots(receiver),
        0.8,
        trajectory_budget(),
        &mut || false,
    )
    .unwrap();
    let config = MovingReceiverConfig {
        receiver_time_interval_s: [0.4, 1.0],
        minimum_distance_m: 1e-5,
        maximum_time_residual_s: 1e-13,
    };
    let pair = PreparedMovingReceiverSource::new(
        &source_path,
        &receiver_path,
        config,
        receiver_budget(),
        &mut || false,
    )
    .unwrap();
    let planes = walls();
    let images =
        MovingReflectionObserver::new(&receiver_path, &planes, 4, reflection_budget(), &mut || {
            false
        })
        .unwrap();
    let reflected = PreparedReflectedMovingReceiverSource::new(
        &source_path,
        &images,
        config,
        reflection_budget(),
        receiver_budget(),
        &mut || false,
    )
    .unwrap();
    for time in [0.4, 0.7, 1.0] {
        let at = receiver(time).position_m;
        let direct = pair
            .prepare_quadrupole_field_receiver(time, &mut || false)
            .unwrap();
        let q = strength(direct.source_time_s());
        let history = observer(at, medium())
            .quadrupole_history_integrand(source, q.stress_integral_n_m_s)
            .unwrap();
        let fixed = PreparedTransientQuadrupoleObserver::new(medium(), at, &[source], 1).unwrap();
        compare_fields(
            direct.evaluate(q, history).unwrap(),
            fixed
                .evaluate_at_time(time, |_, t| Ok(strength(t)))
                .unwrap(),
            5e-10,
        );
        let original = pair
            .prepare_quadrupole_receiver(time, &mut || false)
            .unwrap();
        assert_eq!(
            direct
                .pressure()
                .observed_source_time_rate(receiver(time).velocity_m_s)
                .unwrap(),
            original
                .observed_source_time_rate(receiver(time).velocity_m_s)
                .unwrap()
        );
        let mut pressure = Dd::ZERO;
        let mut velocity = [Dd::ZERO; 3];
        for path in 0..4 {
            let field = reflected
                .prepare_quadrupole_field_path(path, time, &mut || false)
                .unwrap();
            let (virtual_receiver, _, _) = reflected.force_history_path_at(path, time).unwrap();
            let q = strength(field.source_time_s());
            let history = virtual_receiver
                .quadrupole_history_integrand(source, q.stress_integral_n_m_s)
                .unwrap();
            let (p, u) = field.evaluate_dd(q, history.map(Dd::from_f64)).unwrap();
            pressure = pressure + p;
            for i in 0..3 {
                velocity[i] = velocity[i] + u[i];
            }
        }
        compare_sum(
            pressure,
            velocity,
            fixed
                .with_orthogonal_planes(&planes, 4)
                .unwrap()
                .evaluate_at_time(time, |_, t| Ok(strength(t)))
                .unwrap(),
        );
    }
    assert!(matches!(
        pair.prepare_quadrupole_field_receiver(0.7, &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(matches!(
        reflected.prepare_quadrupole_field_path(1, 0.7, &mut || true),
        Err(AeroacError::Cancelled { .. })
    ));
}

#[test]
fn explicit_history_changes_transverse_velocity_and_density_only_scales_velocity() {
    let x = [2.0, -1.3, 1.8];
    let motion = motions()[1];
    let tau = emission(x, 0.7, motion);
    let source = strength(tau);
    let prepare = |rho| {
        observer(x, FluidAcousticState::try_new(rho, C).unwrap())
            .prepare_emission_quadrupole_field(tau, motion.at(tau), motion.jerk)
            .unwrap()
    };
    let sample = prepare(RHO);
    let history = [0.3, -0.2, 0.5];
    let a = sample.evaluate(source, history).unwrap();
    let b = prepare(2.0 * RHO).evaluate(source, history).unwrap();
    assert_eq!(a.pressure_pa, b.pressure_pa);
    for i in 0..3 {
        assert_eq!(a.particle_velocity_m_s[i], 2.0 * b.particle_velocity_m_s[i]);
    }
    let empty = sample.evaluate(source, [0.0; 3]).unwrap();
    assert_eq!(a.pressure_pa, empty.pressure_pa);
    for i in 0..3 {
        close(
            a.particle_velocity_m_s[i] - empty.particle_velocity_m_s[i],
            history[i] / (4.0 * std::f64::consts::PI * RHO),
            norm(a.particle_velocity_m_s),
            4e-15,
        );
    }
    let null = sample
        .evaluate(TransientQuadrupoleSource::default(), history)
        .unwrap();
    assert_eq!(null.pressure_pa, 0.0);
    assert!(norm(null.particle_velocity_m_s) > 0.0);
    assert_eq!(null.instantaneous_intensity_w_m2, [0.0; 3]);
}

#[test]
fn nonfinite_input_refuses_and_scalar_integral_cannot_replace_worldline_history() {
    let obs = observer([2.0, -1.3, 1.8], medium());
    let motion = motions()[1];
    let sample = obs
        .prepare_emission_quadrupole_field(0.1, motion.at(0.1), motion.jerk)
        .unwrap();
    let mut source = strength(0.1);
    let expected = sample.evaluate(source, [0.0; 3]).unwrap();
    assert!(sample.evaluate(source, [f64::NAN, 0.0, 0.0]).is_err());
    assert!(
        sample
            .evaluate_dd(source, [Dd::from_f64(f64::INFINITY); 3])
            .is_err()
    );
    source.stress_rate_n_m_s[4] = f64::NAN;
    assert!(sample.evaluate(source, [0.0; 3]).is_err());
    source = strength(0.1);
    source.stress_integral_n_m_s = [f64::NAN; 6];
    assert_eq!(sample.evaluate(source, [0.0; 3]).unwrap(), expected);
    assert!(
        obs.quadrupole_history_integrand(obs.observer_position_m(), TENSOR)
            .is_err()
    );
    assert!(
        obs.quadrupole_history_integrand([f64::NAN, 0.0, 0.0], TENSOR)
            .is_err()
    );
    assert!(
        obs.quadrupole_history_integrand([0.0; 3], [f64::INFINITY; 6])
            .is_err()
    );
    assert!(
        obs.prepare_emission_quadrupole_field(0.1, motion.at(0.1), [f64::NAN; 3])
            .is_err()
    );
}
