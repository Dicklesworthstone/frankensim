//! G1/G3 moving-source checks against stationary owners and retarded potentials.

use fs_aeroac::curle3d::CompactForceObserver;
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    UniformSourceMotion,
};
use fs_aeroac::monopole::{CompactMonopoleObserver, MonopoleLimits};
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

fn medium(rho: f64, c: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, c).unwrap()
}

fn observer(rho: f64, c: f64, x: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium(rho, c),
        x,
        MovingCompactLimits {
            maximum_mach: 0.9999999,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn stationary_limit_preserves_original_monopole_and_force_observers() {
    let x = [3.0, 4.0, 0.0];
    let fluid = medium(1.2, 300.0);
    let strength = CompactSourceStrength {
        volume_flow_m3_s: 0.008,
        volume_flow_rate_m3_s2: 0.12,
        force_on_fluid_n: [2.0, -3.0, 9.0],
        force_rate_on_fluid_n_s: [200.0, 150.0, -600.0],
    };
    let moving = observer(1.2, 300.0, x)
        .prepare_emission(0.2, SourceKinematics::default())
        .unwrap()
        .evaluate(strength)
        .unwrap();
    let force = CompactForceObserver::new(fluid, x, [0.0; 3])
        .unwrap()
        .pressure_at_retarded(strength.force_on_fluid_n, strength.force_rate_on_fluid_n_s)
        .unwrap();
    let monopole = CompactMonopoleObserver::new(
        fluid,
        5.0,
        100.0,
        MonopoleLimits {
            maximum_ka: 0.1,
            minimum_distance_radii: 10.0,
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap()
    .sample(0.2, 0.01, strength.volume_flow_rate_m3_s2)
    .unwrap();
    close(moving.monopole_pressure_pa, monopole.pressure_pa, 2e-14);
    close(moving.dipole_pressure_pa, force, 2e-14);
    close(moving.pressure_pa, monopole.pressure_pa + force, 2e-14);
    assert_eq!(moving.source_time_rate, 1.0);
    assert_eq!(moving.receiver_time_s, monopole.arrival_time_s);
}

#[test]
fn uniform_approach_and_recession_have_correct_clock_and_amplitude() {
    let c = 300.0;
    let rho = 1.2;
    let receiver_time = 0.02;
    let microphone = observer(rho, c, [12.0, 0.0, 0.0]);
    for speed in [-180.0, 60.0, 240.0] {
        let motion = UniformSourceMotion::new(0.0, [0.0; 3], [speed, 0.0, 0.0]).unwrap();
        let prepared = microphone
            .prepare_uniform(receiver_time, &motion, 2e-14)
            .unwrap();
        // Before crossing, t = tau + (12-v*tau)/c has this linear inverse.
        let kappa = 1.0 - speed / c;
        let source_time = (receiver_time - 12.0 / c) / kappa;
        let distance = 12.0 - speed * source_time;
        let sample = prepared
            .evaluate(CompactSourceStrength {
                volume_flow_rate_m3_s2: 2.0,
                ..CompactSourceStrength::default()
            })
            .unwrap();
        close(sample.source_time_s, source_time, 2e-14);
        close(sample.distance_m, distance, 2e-14);
        close(sample.source_time_rate, 1.0 / kappa, 2e-14);
        close(
            sample.pressure_pa,
            rho * 2.0 / (4.0 * PI * distance * kappa * kappa),
            3e-14,
        );
        assert_eq!(sample.receiver_time_s, receiver_time);
        assert!(sample.retarded_time_residual_s <= 2e-14);
        let next = microphone
            .prepare_uniform(receiver_time + 1e-5, &motion, 2e-14)
            .unwrap();
        close(
            (next.source_time_s() - sample.source_time_s) / 1e-5,
            1.0 / kappa,
            2e-11,
        );
    }
}

#[test]
fn nearly_sonic_recession_keeps_the_small_causal_root() {
    let speed = -0.999999 * 300.0;
    let motion = UniformSourceMotion::new(0.0, [0.0; 3], [speed, 0.0, 0.0]).unwrap();
    let sample = observer(1.2, 300.0, [2.0, 0.0, 0.0])
        .prepare_uniform(0.0, &motion, 2e-14)
        .unwrap();
    let delay = 2.0 / (300.0 - speed);
    close(sample.source_time_s(), -delay, 2e-14);
    close(sample.propagation_delay_s(), delay, 2e-14);
    assert!(sample.source_time_s() < sample.receiver_time_s());
}

#[test]
fn transverse_motion_has_convective_pressure_at_a_stationary_dipole_null() {
    let rho = 1.2;
    let c = 300.0;
    let mach = 0.4;
    let radius = 4.0;
    let sample = observer(rho, c, [0.0, radius, 0.0])
        .prepare_emission(
            0.0,
            SourceKinematics {
                velocity_m_s: [mach * c, 0.0, 0.0],
                ..SourceKinematics::default()
            },
        )
        .unwrap()
        .evaluate(CompactSourceStrength {
            volume_flow_m3_s: 0.002,
            force_on_fluid_n: [3.0, 0.0, 0.0],
            ..CompactSourceStrength::default()
        })
        .unwrap();
    close(
        sample.dipole_pressure_pa,
        -3.0 * mach / (4.0 * PI * radius * radius),
        2e-14,
    );
    close(
        sample.monopole_pressure_pa,
        -rho * 0.002 * c * mach * mach / (4.0 * PI * radius * radius),
        2e-14,
    );
    assert_eq!(sample.radial_mach, 0.0);
    assert_eq!(sample.source_time_rate, 1.0);
}

// Independent manufactured history: the reference below differentiates the
// retarded potentials, never the production coefficients or analytic root.
fn trajectory(t: f64) -> SourceKinematics {
    let position = [0.2, -0.3, 0.1];
    let velocity = [20.0, -15.0, 5.0];
    let acceleration = [400.0, 250.0, -300.0];
    SourceKinematics {
        position_m: core::array::from_fn(|i| {
            position[i] + velocity[i] * t + 0.5 * acceleration[i] * t * t
        }),
        velocity_m_s: core::array::from_fn(|i| velocity[i] + acceleration[i] * t),
        acceleration_m_s2: acceleration,
    }
}

fn history(t: f64) -> CompactSourceStrength {
    let (s, c) = (70.0 * t).sin_cos();
    let (s2, c2) = (31.0 * t).sin_cos();
    CompactSourceStrength {
        volume_flow_m3_s: 0.003 * s + 0.001 * c2,
        volume_flow_rate_m3_s2: 0.21 * c - 0.031 * s2,
        force_on_fluid_n: [2.0 * s, -3.0 * c, 0.5 * s2],
        force_rate_on_fluid_n_s: [140.0 * c, 210.0 * s, 15.5 * c2],
    }
}

fn distance(x: [f64; 3], y: [f64; 3]) -> f64 {
    ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt()
}

fn retarded_time(x: [f64; 3], receiver_time: f64) -> f64 {
    // This whole bracket has |v|<c and therefore a unique monotone arrival map.
    let mut lo = -0.2;
    let mut hi = receiver_time;
    assert!(lo + distance(x, trajectory(lo).position_m) / 300.0 < receiver_time);
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if mid + distance(x, trajectory(mid).position_m) / 300.0 < receiver_time {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn potentials(x: [f64; 3], receiver_time: f64) -> [f64; 4] {
    let tau = retarded_time(x, receiver_time);
    let state = trajectory(tau);
    let r = distance(x, state.position_m);
    let radial_velocity: f64 = (0..3)
        .map(|i| (x[i] - state.position_m[i]) * state.velocity_m_s[i] / r)
        .sum();
    let denominator = 4.0 * PI * r * (1.0 - radial_velocity / 300.0);
    let source = history(tau);
    [
        source.volume_flow_m3_s / denominator,
        source.force_on_fluid_n[0] / denominator,
        source.force_on_fluid_n[1] / denominator,
        source.force_on_fluid_n[2] / denominator,
    ]
}

fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
}

#[test]
fn accelerating_source_matches_eulerian_derivatives_of_retarded_potentials() {
    let x = [8.0, 5.0, -3.0];
    for t in [-0.02, 0.04, 0.09] {
        let tau = retarded_time(x, t);
        let sample = observer(1.2, 300.0, x)
            .prepare_emission(tau, trajectory(tau))
            .unwrap()
            .evaluate(history(tau))
            .unwrap();
        let monopole = 1.2 * derivative(|dt| potentials(x, t + dt)[0], 2e-5);
        let mut dipole = 0.0;
        for axis in 0..3 {
            dipole -= derivative(
                |dx| {
                    let mut shifted = x;
                    shifted[axis] += dx;
                    potentials(shifted, t)[axis + 1]
                },
                2e-3,
            );
        }
        close(sample.monopole_pressure_pa, monopole, 2e-7);
        close(sample.dipole_pressure_pa, dipole, 2e-7);
        close(sample.pressure_pa, monopole + dipole, 2e-7);
        close(sample.receiver_time_s, t, 2e-13);
    }
}

#[test]
fn moving_force_remains_density_independent_and_volume_pressure_scales_with_density() {
    let x = [3.0, 4.0, 2.0];
    let kinematics = trajectory(0.01);
    let source = history(0.01);
    let evaluate = |rho, strength| {
        observer(rho, 300.0, x)
            .prepare_emission(0.01, kinematics)
            .unwrap()
            .evaluate(strength)
            .unwrap()
    };
    let baseline = evaluate(1.0, source);
    for rho in [0.125, 8.0, 1000.0] {
        let sample = evaluate(rho, source);
        assert_eq!(sample.dipole_pressure_pa, baseline.dipole_pressure_pa);
        close(
            sample.monopole_pressure_pa,
            rho * baseline.monopole_pressure_pa,
            3e-14,
        );
        assert_eq!(sample.receiver_time_s, baseline.receiver_time_s);
    }
    let reversed = evaluate(
        1.0,
        CompactSourceStrength {
            volume_flow_m3_s: -source.volume_flow_m3_s,
            volume_flow_rate_m3_s2: -source.volume_flow_rate_m3_s2,
            force_on_fluid_n: source.force_on_fluid_n.map(|v| -v),
            force_rate_on_fluid_n_s: source.force_rate_on_fluid_n_s.map(|v| -v),
        },
    );
    assert_eq!(reversed.pressure_pa, -baseline.pressure_pa);
}

#[test]
fn singular_sonic_nonfinite_and_unrepresentable_events_refuse() {
    let microphone = observer(1.2, 300.0, [1.0, 0.0, 0.0]);
    for state in [
        SourceKinematics {
            position_m: [1.0, 0.0, 0.0],
            ..SourceKinematics::default()
        },
        SourceKinematics {
            velocity_m_s: [300.0, 0.0, 0.0],
            ..SourceKinematics::default()
        },
        SourceKinematics {
            acceleration_m_s2: [f64::NAN, 0.0, 0.0],
            ..SourceKinematics::default()
        },
    ] {
        assert!(microphone.prepare_emission(0.0, state).is_err());
    }
    let static_motion = UniformSourceMotion::new(0.0, [0.0; 3], [0.0; 3]).unwrap();
    assert!(
        microphone
            .prepare_uniform(1e16, &static_motion, 1.0)
            .is_err()
    );
    assert!(
        microphone
            .prepare_uniform(0.0, &static_motion, -1.0)
            .is_err()
    );
    let prepared = microphone
        .prepare_emission(0.0, SourceKinematics::default())
        .unwrap();
    assert!(
        prepared
            .evaluate(CompactSourceStrength {
                force_rate_on_fluid_n_s: [0.0, f64::INFINITY, 0.0],
                ..CompactSourceStrength::default()
            })
            .is_err()
    );
    assert_eq!(
        prepared
            .evaluate(CompactSourceStrength::default())
            .unwrap()
            .pressure_pa,
        0.0
    );
}
