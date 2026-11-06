//! G1/G3 accelerating stress potentials, original trajectory jets and clocks.

use fs_aeroac::curle3d::moving::trajectory::{
    PreparedSourceTrajectory, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    MovingCompactLimits, MovingCompactObserver, SourceKinematics, UniformSourceMotion,
};
use fs_aeroac::curle3d::quadrupole::transient::TransientQuadrupoleSource;
use fs_material::fluid_acoustics::FluidAcousticState;

const TENSOR: [f64; 6] = [0.7, -0.2, 0.4, 0.3, -0.1, 0.23];
const ENTRIES: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
const C: f64 = 343.0;

fn observer(x: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, C).unwrap(),
        x,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}

fn path(
    observer: MovingCompactObserver,
    times: &[f64],
    motion: impl Fn(f64) -> SourceKinematics,
) -> PreparedSourceTrajectory {
    let samples: Vec<_> = times
        .iter()
        .map(|&source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect();
    PreparedSourceTrajectory::new(
        observer,
        &samples,
        TrajectoryBudget {
            maximum_segments: 16,
            maximum_payload_bytes: 100_000,
            maximum_bound_nodes: 8192,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}

fn strength(t: f64, omega: f64) -> TransientQuadrupoleSource {
    let (s, c) = (omega * t).sin_cos();
    TransientQuadrupoleSource {
        stress_moment_n_m: TENSOR.map(|v| v * c),
        stress_rate_n_m_s: TENSOR.map(|v| -omega * v * s),
        stress_acceleration_n_m_s2: TENSOR.map(|v| -omega * omega * v * c),
        ..Default::default()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[derive(Clone, Copy)]
struct Cubic {
    v: [f64; 3],
}
impl Cubic {
    const JERK: [f64; 3] = [8000.0, 3000.0, -4000.0];
    fn at(self, t: f64) -> SourceKinematics {
        let p = [-0.2, 0.1, 0.3];
        let a = [700.0, -500.0, 400.0];
        SourceKinematics {
            position_m: core::array::from_fn(|i| {
                p[i] + self.v[i] * t + a[i] * t * t / 2.0 + Self::JERK[i] * t * t * t / 6.0
            }),
            velocity_m_s: core::array::from_fn(|i| {
                self.v[i] + a[i] * t + Self::JERK[i] * t * t / 2.0
            }),
            acceleration_m_s2: core::array::from_fn(|i| a[i] + Self::JERK[i] * t),
        }
    }
}

#[test]
fn stationary_and_uniform_trajectories_recover_the_original_prepared_pressure() {
    let obs = observer([4.0, -2.0, 3.0]);
    for velocity in [[0.0; 3], [80.0, 0.0, 0.0], [-50.0, 20.0, 0.0]] {
        let motion = UniformSourceMotion::new(0.0, [-0.2, 0.1, 0.3], velocity).unwrap();
        let trajectory = path(obs, &[-0.04, 0.01, 0.06], |t| {
            motion.kinematics_at(t).unwrap()
        });
        for time in [0.01, 0.03, 0.05] {
            let original = obs
                .prepare_uniform_quadrupole(time, &motion, 1e-13)
                .unwrap();
            let curved = trajectory
                .prepare_quadrupole_receiver(time, 1e-13, &mut || false)
                .unwrap();
            near(curved.source_time_s(), original.source_time_s(), 2e-13);
            assert_eq!(curved.propagation().receiver_time_s(), time);
            for omega in [0.0, 250.0] {
                let expected = original
                    .pressure_pa(strength(original.source_time_s(), omega))
                    .unwrap();
                near(
                    curved
                        .pressure_pa(strength(curved.source_time_s(), omega))
                        .unwrap(),
                    expected,
                    2e-9 * expected.abs().max(1e-5),
                );
            }
        }
    }
}

// Independently bisect the physical sound cone and then evaluate only the
// scalar retarded potential. No production inverse or derivative is used.
fn potential(x: [f64; 3], t: f64, motion: Cubic, omega: f64) -> f64 {
    let mut lo = -0.04;
    let mut hi = 0.06;
    for _ in 0..100 {
        let tau = 0.5 * (lo + hi);
        let state = motion.at(tau);
        let range = (0..3)
            .map(|i| (x[i] - state.position_m[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        if tau + range / C > t {
            hi = tau;
        } else {
            lo = tau;
        }
    }
    let tau = 0.5 * (lo + hi);
    let state = motion.at(tau);
    let r: [f64; 3] = core::array::from_fn(|i| x[i] - state.position_m[i]);
    let range = r.iter().map(|v| v * v).sum::<f64>().sqrt();
    let radial = (0..3).map(|i| r[i] * state.velocity_m_s[i]).sum::<f64>() / range;
    (omega * tau).cos() / (4.0 * std::f64::consts::PI * range * (1.0 - radial / C))
}

fn hessian(x: [f64; 3], t: f64, motion: Cubic, omega: f64, h: f64) -> f64 {
    let mut pressure = 0.0;
    for (entry, (i, j)) in ENTRIES.into_iter().enumerate() {
        let derivative = if i == j {
            let mut plus = x;
            plus[i] += h;
            let mut minus = x;
            minus[i] -= h;
            (potential(plus, t, motion, omega) - 2.0 * potential(x, t, motion, omega)
                + potential(minus, t, motion, omega))
                / (h * h)
        } else {
            let mut sum = 0.0;
            for (a, b) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
                let mut shifted = x;
                shifted[i] += a * h;
                shifted[j] += b * h;
                sum += a * b * potential(shifted, t, motion, omega);
            }
            2.0 * sum / (4.0 * h * h)
        };
        pressure += TENSOR[entry] * derivative;
    }
    pressure
}

#[test]
fn accelerating_pressure_converges_to_an_independent_retarded_potential_hessian() {
    let x = [4.0, -2.0, 3.0];
    for velocity in [[80.0, 0.0, 0.0], [-80.0, 20.0, 0.0], [0.0, 100.0, -30.0]] {
        let motion = Cubic { v: velocity };
        let trajectory = path(observer(x), &[-0.04, 0.01, 0.06], |t| motion.at(t));
        let prepared = trajectory
            .prepare_quadrupole_receiver(0.03, 1e-13, &mut || false)
            .unwrap();
        for omega in [0.0, 250.0] {
            let expected = prepared
                .pressure_pa(strength(prepared.source_time_s(), omega))
                .unwrap();
            let errors =
                [0.02, 0.01, 0.005].map(|h| (hessian(x, 0.03, motion, omega, h) - expected).abs());
            assert!(
                errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
                "v={velocity:?}, omega={omega}, errors={errors:?}"
            );
            assert!(
                errors[2] < 1e-5 * expected.abs().max(1e-4),
                "pressure={expected:e}, errors={errors:?}"
            );
        }
    }
}

#[test]
fn trajectory_jerk_is_an_analytic_polynomial_derivative_with_explicit_c2_knot_sides() {
    let polynomial = |t: f64| SourceKinematics {
        position_m: [
            0.2 + 10.0 * t + 20.0 * t * t - 15.0 * t.powi(3) + 7.0 * t.powi(4) + 9.0 * t.powi(5),
            0.0,
            0.0,
        ],
        velocity_m_s: [
            10.0 + 40.0 * t - 45.0 * t * t + 28.0 * t.powi(3) + 45.0 * t.powi(4),
            0.0,
            0.0,
        ],
        acceleration_m_s2: [40.0 - 90.0 * t + 84.0 * t * t + 180.0 * t.powi(3), 0.0, 0.0],
    };
    let trajectory = path(
        observer([8.0, 5.0, 3.0]),
        &[-0.2, -0.03, 0.11, 0.4],
        polynomial,
    );
    for t in [-0.2, -0.12, -0.03, 0.0, 0.11, 0.28, 0.4] {
        let (state, jerk) = trajectory.kinematics_with_jerk_at(t).unwrap();
        assert_eq!(state, trajectory.kinematics_at(t).unwrap());
        near(jerk[0], -90.0 + 168.0 * t + 540.0 * t * t, 1e-8);
        assert_eq!([jerk[1], jerk[2]], [0.0; 2]);
    }
    let piecewise = path(observer([4.0, 2.0, 1.0]), &[-0.1, 0.0, 0.1], |t| {
        let t = t.max(0.0);
        SourceKinematics {
            position_m: [1000.0 * t * t * t, 0.0, 0.0],
            velocity_m_s: [3000.0 * t * t, 0.0, 0.0],
            acceleration_m_s2: [6000.0 * t, 0.0, 0.0],
        }
    });
    assert_eq!(piecewise.kinematics_with_jerk_at(0.0).unwrap().1, [0.0; 3]);
    near(
        piecewise.kinematics_with_jerk_at(0.01).unwrap().1[0],
        6000.0,
        1e-8,
    );
    assert!(piecewise.kinematics_with_jerk_at(-0.11).is_err());
    assert!(piecewise.kinematics_with_jerk_at(f64::NAN).is_err());
}

#[test]
fn source_jerk_changes_constant_stress_radiation_and_rotates_with_the_tensor() {
    let obs = observer([4.0, -2.0, 3.0]);
    let source = Cubic {
        v: [80.0, 0.0, 0.0],
    }
    .at(0.01);
    let base = obs
        .prepare_emission_quadrupole(0.01, source, [0.0; 3])
        .unwrap();
    let with_jerk = obs
        .prepare_emission_quadrupole(0.01, source, Cubic::JERK)
        .unwrap();
    let radius = base.propagation().distance_m();
    let n: [f64; 3] =
        core::array::from_fn(|i| (obs.observer_position_m()[i] - source.position_m[i]) / radius);
    let radial_jerk = (0..3).map(|i| n[i] * Cubic::JERK[i]).sum::<f64>();
    let projected_stress = ENTRIES
        .into_iter()
        .enumerate()
        .map(|(e, (i, j))| TENSOR[e] * n[i] * n[j] * if i == j { 1.0 } else { 2.0 })
        .sum::<f64>();
    let kappa = 1.0 - base.propagation().radial_mach();
    let expected = radial_jerk * projected_stress
        / (4.0 * std::f64::consts::PI * radius * C.powi(3) * kappa.powi(4));
    let base_pressure = base.pressure_pa(strength(0.0, 0.0)).unwrap();
    near(
        with_jerk.pressure_pa(strength(0.0, 0.0)).unwrap() - base_pressure,
        expected,
        // Subtracting the much larger acceleration/near-field pressure sets
        // the floating-point scale of this isolated small jerk contribution.
        64.0 * f64::EPSILON * base_pressure.abs().max(expected.abs()),
    );

    let rotate = |x: [f64; 3]| [x[2], x[0], x[1]];
    let rotated = observer(rotate(obs.observer_position_m()))
        .prepare_emission_quadrupole(
            0.01,
            SourceKinematics {
                position_m: rotate(source.position_m),
                velocity_m_s: rotate(source.velocity_m_s),
                acceleration_m_s2: rotate(source.acceleration_m_s2),
            },
            rotate(Cubic::JERK),
        )
        .unwrap();
    let q = strength(0.01, 250.0);
    let rotate_tensor = |t: [f64; 6]| [t[2], t[0], t[1], t[4], t[5], t[3]];
    let rq = TransientQuadrupoleSource {
        stress_moment_n_m: rotate_tensor(q.stress_moment_n_m),
        stress_rate_n_m_s: rotate_tensor(q.stress_rate_n_m_s),
        stress_acceleration_n_m_s2: rotate_tensor(q.stress_acceleration_n_m_s2),
        ..Default::default()
    };
    near(
        rotated.pressure_pa(rq).unwrap(),
        with_jerk.pressure_pa(q).unwrap(),
        1e-15,
    );
}

#[test]
fn invalid_jerk_out_of_path_receiver_times_and_cancellation_refuse() {
    let obs = observer([4.0, -2.0, 3.0]);
    let motion = Cubic {
        v: [80.0, 0.0, 0.0],
    };
    assert!(
        obs.prepare_emission_quadrupole(0.01, motion.at(0.01), [f64::NAN, 0.0, 0.0])
            .is_err()
    );
    assert!(
        obs.prepare_emission_quadrupole(0.01, motion.at(0.01), [0.0, f64::INFINITY, 0.0])
            .is_err()
    );
    let trajectory = path(obs, &[-0.04, 0.06], |t| motion.at(t));
    assert!(
        trajectory
            .prepare_quadrupole_receiver(1.0, 1e-12, &mut || false)
            .is_err()
    );
    assert!(
        trajectory
            .prepare_quadrupole_receiver(0.03, 1e-12, &mut || true)
            .is_err()
    );
}
