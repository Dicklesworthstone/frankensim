//! Full moving-force fields retain the retarded endpoint and its past worldline.

use fs_aeroac::curle3d::CompactForceObserver;
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 20.0).unwrap()
}

fn observer(position: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium(),
        position,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 1e-4,
        },
    )
    .unwrap()
}

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn force(time: f64) -> ([f64; 3], [f64; 3]) {
    let quadratic = [2.0, -1.0, 0.7];
    let cubic = [1.0, 0.3, -0.4];
    (
        core::array::from_fn(|i| time * time * (quadratic[i] + time * cubic[i])),
        core::array::from_fn(|i| time * (2.0 * quadratic[i] + 3.0 * time * cubic[i])),
    )
}

fn motion(time: f64) -> SourceKinematics {
    let start = [0.2, -0.5, 0.3];
    let velocity = [8.0, -5.0, 3.0];
    let acceleration = [-4.0, 6.0, -2.0];
    SourceKinematics {
        position_m: core::array::from_fn(|i| {
            start[i] + time * velocity[i] + 0.5 * time * time * acceleration[i]
        }),
        velocity_m_s: core::array::from_fn(|i| velocity[i] + time * acceleration[i]),
        acceleration_m_s2: acceleration,
    }
}

fn integrate(end: f64, integrand: impl Fn(f64) -> [f64; 3]) -> [f64; 3] {
    // Composite Simpson quadrature is an independent test oracle, never the
    // production history interpolator or its adaptive integration algorithm.
    const PANELS: usize = 1024;
    let mut sum = [0.0; 3];
    for index in 0..=PANELS {
        let weight = if index == 0 || index == PANELS {
            1.0
        } else if index % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let value = integrand(end * index as f64 / PANELS as f64);
        for axis in 0..3 {
            sum[axis] += weight * value[axis];
        }
    }
    sum.map(|v| v * end / (3.0 * PANELS as f64))
}

#[test]
fn stationary_limit_recovers_original_force_owner_and_coherent_monopole() {
    let position = [4.1, 1.9, -0.8];
    let source = [0.3, -0.7, 1.1];
    let original = CompactForceObserver::new(medium(), position, source).unwrap();
    let moving = observer(position);
    let kinematics = SourceKinematics {
        position_m: source,
        velocity_m_s: [0.0; 3],
        acceleration_m_s2: [0.0; 3],
    };
    let impulse = [0.7, -0.4, 0.6];
    let strength = CompactSourceStrength {
        volume_flow_m3_s: 0.13,
        volume_flow_rate_m3_s2: -0.17,
        force_on_fluid_n: [2.0, -1.0, 3.0],
        force_rate_on_fluid_n_s: [-4.0, 0.5, 6.0],
    };
    // At rest the geometric Hessian is constant, so its history integral
    // equals the integrand applied to the ORIGINAL force impulse.
    let history = moving.force_history_integrand(source, impulse).unwrap();
    let actual = moving
        .prepare_force_emission(0.4, kinematics)
        .unwrap()
        .evaluate(strength, history)
        .unwrap();
    let dipole = original
        .field_at_retarded(
            impulse,
            strength.force_on_fluid_n,
            strength.force_rate_on_fluid_n_s,
        )
        .unwrap();
    let prepared = moving.prepare_emission(0.4, kinematics).unwrap();
    let monopole = prepared
        .evaluate_monopole_field(strength.volume_flow_m3_s, strength.volume_flow_rate_m3_s2)
        .unwrap();
    let pressure = dipole.pressure_pa + monopole.sample.pressure_pa;
    close(actual.sample.pressure_pa, pressure, 2e-13);
    close(
        actual.sample.pressure_pa,
        prepared.evaluate(strength).unwrap().pressure_pa,
        2e-13,
    );
    for axis in 0..3 {
        let velocity = dipole.particle_velocity_m_s[axis] + monopole.particle_velocity_m_s[axis];
        close(actual.particle_velocity_m_s[axis], velocity, 3e-13);
        close(
            actual.instantaneous_intensity_w_m2[axis],
            pressure * velocity,
            4e-13,
        );
    }
}

// Independent vector potential B(x,t)=integral_0^tau F(s)/|x-y(s)| ds.
// Both the root and the integral are recomputed at each displaced receiver;
// freezing tau in the finite differences would omit the retarded boundary.
fn potential(position: [f64; 3], receiver_time: f64) -> [f64; 3] {
    let distance = |time| {
        let source = motion(time).position_m;
        norm(core::array::from_fn(|i| position[i] - source[i]))
    };
    let mut lo = 0.0;
    let mut hi = 1.5;
    assert!(lo + distance(lo) / 20.0 < receiver_time);
    assert!(hi + distance(hi) / 20.0 > receiver_time);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        if mid + distance(mid) / 20.0 < receiver_time {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    integrate(0.5 * (lo + hi), |time| {
        let radius = distance(time);
        force(time).0.map(|value| value / radius)
    })
}

fn potential_hessian_velocity(position: [f64; 3], time: f64, step: f64) -> [f64; 3] {
    let center = potential(position, time);
    core::array::from_fn(|i| {
        let mut value = 0.0;
        for j in 0..3 {
            if i == j {
                let mut plus = position;
                let mut minus = position;
                plus[i] += step;
                minus[i] -= step;
                value += (potential(plus, time)[j] - 2.0 * center[j] + potential(minus, time)[j])
                    / (step * step);
            } else {
                let shifted = |a: f64, b: f64| {
                    let mut point = position;
                    point[i] += a * step;
                    point[j] += b * step;
                    potential(point, time)[j]
                };
                value += (shifted(1.0, 1.0) - shifted(1.0, -1.0) - shifted(-1.0, 1.0)
                    + shifted(-1.0, -1.0))
                    / (4.0 * step * step);
            }
        }
        value / (4.0 * core::f64::consts::PI * medium().density_kg_m3())
    })
}

#[test]
fn accelerating_transverse_fields_match_independent_retarded_potential_hessian() {
    let source_time = 0.65;
    let (force, rate) = force(source_time);
    for position in [[10.0, 5.0, -3.0], [16.0, -6.0, 4.0], [-8.0, 4.0, -2.0]] {
        let receiver = observer(position);
        let history = integrate(source_time, |time| {
            receiver
                .force_history_integrand(motion(time).position_m, self::force(time).0)
                .unwrap()
        });
        let actual = receiver
            .prepare_force_emission(source_time, motion(source_time))
            .unwrap()
            .evaluate(
                CompactSourceStrength {
                    force_on_fluid_n: force,
                    force_rate_on_fluid_n_s: rate,
                    ..CompactSourceStrength::default()
                },
                history,
            )
            .unwrap();
        let mut errors = [0.0_f64; 3];
        for (index, step) in [0.04, 0.02, 0.01].into_iter().enumerate() {
            let reference =
                potential_hessian_velocity(position, actual.sample.receiver_time_s, step);
            errors[index] = norm(core::array::from_fn(|axis| {
                actual.particle_velocity_m_s[axis] - reference[axis]
            }));
        }
        let scale = norm(actual.particle_velocity_m_s);
        assert!(
            errors[2] < 3e-6 * scale,
            "errors={errors:?}, scale={scale:e}"
        );
        assert!(errors[1] < 0.4 * errors[0], "errors={errors:?}");
        assert!(errors[2] < 0.5 * errors[1], "errors={errors:?}");
    }
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    let k = core::f64::consts::FRAC_1_SQRT_2;
    [(v[0] - v[1]) * k, (v[0] + v[1]) * k, v[2]]
}

#[test]
fn full_fields_and_worldline_integrand_are_rotation_covariant() {
    let position = [10.0, 5.0, -3.0];
    let kinematics = motion(0.65);
    let (force, rate) = force(0.65);
    let strength = CompactSourceStrength {
        volume_flow_m3_s: 0.2,
        volume_flow_rate_m3_s2: -0.3,
        force_on_fluid_n: force,
        force_rate_on_fluid_n_s: rate,
    };
    let receiver = observer(position);
    let rotated = observer(rotate(position));
    let history = integrate(0.65, |time| {
        receiver
            .force_history_integrand(motion(time).position_m, self::force(time).0)
            .unwrap()
    });
    let expected_integrand = rotate(
        receiver
            .force_history_integrand(kinematics.position_m, force)
            .unwrap(),
    );
    let actual_integrand = rotated
        .force_history_integrand(rotate(kinematics.position_m), rotate(force))
        .unwrap();
    let actual = receiver
        .prepare_force_emission(0.65, kinematics)
        .unwrap()
        .evaluate(strength, history)
        .unwrap();
    let transformed = rotated
        .prepare_force_emission(
            0.65,
            SourceKinematics {
                position_m: rotate(kinematics.position_m),
                velocity_m_s: rotate(kinematics.velocity_m_s),
                acceleration_m_s2: rotate(kinematics.acceleration_m_s2),
            },
        )
        .unwrap()
        .evaluate(
            CompactSourceStrength {
                force_on_fluid_n: rotate(force),
                force_rate_on_fluid_n_s: rotate(rate),
                ..strength
            },
            rotate(history),
        )
        .unwrap();
    close(
        transformed.sample.pressure_pa,
        actual.sample.pressure_pa,
        3e-13,
    );
    let velocity = rotate(actual.particle_velocity_m_s);
    let intensity = rotate(actual.instantaneous_intensity_w_m2);
    for axis in 0..3 {
        close(actual_integrand[axis], expected_integrand[axis], 3e-13);
        close(
            transformed.particle_velocity_m_s[axis],
            velocity[axis],
            5e-13,
        );
        close(
            transformed.instantaneous_intensity_w_m2[axis],
            intensity[axis],
            6e-13,
        );
    }
}

#[test]
fn nonfinite_history_and_source_collisions_refuse_without_partial_fields() {
    let receiver = observer([10.0, 5.0, -3.0]);
    let kinematics = motion(0.65);
    let prepared = receiver.prepare_force_emission(0.65, kinematics).unwrap();
    assert!(
        prepared
            .evaluate(CompactSourceStrength::default(), [f64::NAN, 0.0, 0.0])
            .is_err()
    );
    assert!(
        prepared
            .evaluate(
                CompactSourceStrength {
                    force_on_fluid_n: [f64::INFINITY, 0.0, 0.0],
                    ..CompactSourceStrength::default()
                },
                [0.0; 3]
            )
            .is_err()
    );
    assert!(
        receiver
            .force_history_integrand(receiver.observer_position_m(), [1.0; 3])
            .is_err()
    );
    assert!(
        receiver
            .force_history_integrand([f64::NAN, 0.0, 0.0], [1.0; 3])
            .is_err()
    );
    assert!(
        receiver
            .force_history_integrand(kinematics.position_m, [f64::NAN; 3])
            .is_err()
    );
    assert!(
        receiver
            .prepare_force_emission(f64::NAN, kinematics)
            .is_err()
    );
    assert!(
        receiver
            .prepare_force_emission(
                0.65,
                SourceKinematics {
                    velocity_m_s: [20.0, 0.0, 0.0],
                    ..kinematics
                }
            )
            .is_err()
    );
    let zero = prepared
        .evaluate(CompactSourceStrength::default(), [0.0; 3])
        .unwrap();
    assert_eq!(zero.sample.pressure_pa, 0.0);
    assert_eq!(zero.particle_velocity_m_s, [0.0; 3]);
    assert_eq!(zero.instantaneous_intensity_w_m2, [0.0; 3]);
}
