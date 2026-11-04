//! G1/G3: uniform-flow pressure against the independent convected Green function.

use fs_aeroac::curle3d::moving::flow::UniformFlowCompactObserver;
use fs_aeroac::curle3d::moving::{CompactSourceStrength, MovingCompactLimits};
use fs_aeroac::curle3d::{CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn fluid(rho: f64) -> FluidAcousticState {
    FluidAcousticState::try_new(rho, 300.0).unwrap()
}

fn observer(flow: [f64; 3], position: [f64; 3]) -> UniformFlowCompactObserver {
    UniformFlowCompactObserver::new(
        fluid(1.2),
        flow,
        position,
        [0.0; 3],
        MovingCompactLimits {
            maximum_mach: 0.9999999,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-10),
        "actual={actual:e}, expected={expected:e}, relative tolerance={tol:e}"
    );
}

fn complex_close(actual: C64, expected: C64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-12),
        "actual={actual:?}, expected={expected:?}, relative tolerance={tol:e}"
    );
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn source() -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: C64::new(0.002, -0.001),
        force_on_fluid_n: [C64::new(1.0, -0.4), C64::new(0.3, 0.7), C64::new(-0.8, 0.2)],
    }
}

// Independently differentiate G=exp(i*w*T)/(4*pi*R*), without using moving
// source kinematics or their coefficients. This is the physical source PDE
// rho(-i*w+U.grad)G Q - F.grad(G), retaining the convective Q term.
fn green_pressure(
    rho: f64,
    c: f64,
    flow: [f64; 3],
    r: [f64; 3],
    frequency: f64,
    source: CompactSourcePhasor,
) -> C64 {
    let omega = TAU * frequency;
    let mach = flow.map(|v| v / c);
    let beta2 = 1.0 - dot(mach, mach);
    let mr = dot(mach, r);
    let rstar = (beta2 * dot(r, r) + mr * mr).sqrt();
    let travel = (rstar - mr) / (c * beta2);
    let (sin, cos) = (omega * travel).sin_cos();
    let g = C64::new(cos, sin).scale(1.0 / (4.0 * PI * rstar));
    let mut pressure = g * source.volume_flow_m3_s * C64::new(0.0, -rho * omega);
    for i in 0..3 {
        let dr = (beta2 * r[i] + mr * mach[i]) / rstar;
        let dt = (dr - mach[i]) / (c * beta2);
        let gradient = g * C64::new(-dr / rstar, omega * dt);
        pressure = pressure
            + gradient
                * (source.volume_flow_m3_s.scale(rho * flow[i]) - source.force_on_fluid_n[i]);
    }
    pressure
}

fn at_time(source: CompactSourcePhasor, frequency: f64, time: f64) -> CompactSourceStrength {
    let omega = TAU * frequency;
    let (sin, cos) = (-omega * time).sin_cos();
    let phase = C64::new(cos, sin);
    let rate = C64::new(0.0, -omega);
    CompactSourceStrength {
        volume_flow_m3_s: (source.volume_flow_m3_s * phase).re,
        volume_flow_rate_m3_s2: (source.volume_flow_m3_s * rate * phase).re,
        force_on_fluid_n: source.force_on_fluid_n.map(|f| (f * phase).re),
        force_rate_on_fluid_n_s: source.force_on_fluid_n.map(|f| (f * rate * phase).re),
    }
}

#[test]
fn zero_flow_recovers_original_stationary_harmonic_pressure() {
    let position = [3.0, 4.0, -2.0];
    let wind = observer([0.0; 3], position);
    for frequency in [0.1, 125.0, 2300.0] {
        let original =
            PreparedCompactSourceObserver::new(fluid(1.2), frequency, position, &[[0.0; 3]], 1)
                .unwrap()
                .evaluate(&[source()])
                .unwrap();
        complex_close(
            wind.pressure_phasor(frequency, source()).unwrap(),
            original.pressure_pa,
            3e-13,
        );
    }
    close(
        wind.propagation_delay_s(),
        dot(position, position).sqrt() / 300.0,
        3e-15,
    );
    assert_eq!(wind.medium(), fluid(1.2));
    assert_eq!(wind.source_position_m(), [0.0; 3]);
    assert_eq!(wind.observer_position_m(), position);
}

#[test]
fn axial_and_crosswind_delays_preserve_the_lab_source_clock() {
    let flow = [180.0, 0.0, 0.0];
    for (position, delay) in [
        ([12.0, 0.0, 0.0], 12.0 / 480.0),
        ([-12.0, 0.0, 0.0], 12.0 / 120.0),
        ([0.0, 12.0, 0.0], 12.0 / (300.0 * 0.8)),
    ] {
        let wind = observer(flow, position);
        close(wind.propagation_delay_s(), delay, 3e-15);
        let first = wind.prepare_receiver(0.2, 1e-15).unwrap();
        let next = wind.prepare_receiver(0.201, 1e-15).unwrap();
        close(next.source_time_s() - first.source_time_s(), 0.001, 2e-13);
        close(first.source_time_s(), 0.2 - delay, 3e-15);
        assert_eq!(first.receiver_time_s(), 0.2);
        let sample = first.evaluate(CompactSourceStrength::default()).unwrap();
        assert_eq!(sample.pressure_pa, 0.0);
        assert!(sample.retarded_time_residual_s <= 1e-15);
    }
}

#[test]
fn complete_pressure_matches_the_independent_convected_green_function() {
    for flow in [[0.0, 0.0, 0.0], [170.0, -30.0, 50.0], [-70.0, 120.0, 10.0]] {
        for position in [[0.3, 0.4, -0.2], [3.0, -2.0, 4.0], [-30.0, 20.0, 40.0]] {
            for frequency in [0.2, 70.0, 1800.0] {
                let wind = observer(flow, position);
                // Separate volume and Cartesian force terms ensure that a
                // fortuitous cancellation cannot hide a wrong source law.
                for input in [
                    source(),
                    CompactSourcePhasor::monopole(source().volume_flow_m3_s),
                    CompactSourcePhasor::force_dipole([C64::ONE, C64::ZERO, C64::ZERO]),
                    CompactSourcePhasor::force_dipole([C64::ZERO, C64::ONE, C64::ZERO]),
                    CompactSourcePhasor::force_dipole([C64::ZERO, C64::ZERO, C64::ONE]),
                ] {
                    complex_close(
                        wind.pressure_phasor(frequency, input).unwrap(),
                        green_pressure(1.2, 300.0, flow, position, frequency, input),
                        8e-12,
                    );
                }
            }
        }
    }
}

#[test]
fn transient_and_harmonic_views_apply_the_same_physical_delay_once() {
    let wind = observer([160.0, -30.0, 50.0], [4.0, 3.0, -2.0]);
    for frequency in [0.3, 37.0, 1400.0] {
        let phasor = green_pressure(
            1.2,
            300.0,
            wind.velocity_m_s(),
            [4.0, 3.0, -2.0],
            frequency,
            source(),
        );
        for time in [0.0, 0.0004, 0.003, 0.017] {
            let sample = wind.prepare_receiver(time, 1e-15).unwrap();
            let result = sample
                .evaluate(at_time(source(), frequency, sample.source_time_s()))
                .unwrap();
            let (sin, cos) = (-TAU * frequency * time).sin_cos();
            let expected = (phasor * C64::new(cos, sin)).re;
            close(result.pressure_pa, expected, 1e-11);
            close(
                result.monopole_pressure_pa + result.dipole_pressure_pa,
                expected,
                1e-11,
            );
            assert_eq!(result.receiver_time_s, time);
        }
    }
}

#[test]
fn steady_volume_injection_retains_convective_near_field_pressure() {
    let q = 0.003;
    // D[Q delta_y] contains U.grad(delta_y) even for constant Q. Dropping
    // this term would falsely make a steady injected volume acoustically zero.
    for x in [-5.0_f64, 5.0] {
        let wind = observer([120.0, 0.0, 0.0], [x, 0.0, 0.0]);
        let result = wind
            .prepare_receiver(1.0, 1e-15)
            .unwrap()
            .evaluate(CompactSourceStrength {
                volume_flow_m3_s: q,
                ..CompactSourceStrength::default()
            })
            .unwrap();
        let expected = -1.2 * 120.0 * q * x.signum() / (4.0 * PI * x * x);
        close(result.pressure_pa, expected, 2e-14);
        assert_eq!(result.dipole_pressure_pa, 0.0);
    }
}

#[test]
fn flow_geometry_rotation_and_material_density_preserve_source_units() {
    let flow = [90.0, -40.0, 30.0];
    let position = [4.0, 2.0, -3.0];
    let rotate = |a: [f64; 3]| [-a[2], a[0], -a[1]];
    let input = source();
    let rotated = CompactSourcePhasor {
        volume_flow_m3_s: input.volume_flow_m3_s,
        force_on_fluid_n: [
            -input.force_on_fluid_n[2],
            input.force_on_fluid_n[0],
            -input.force_on_fluid_n[1],
        ],
    };
    let reference = observer(flow, position)
        .pressure_phasor(83.0, input)
        .unwrap();
    complex_close(
        observer(rotate(flow), rotate(position))
            .pressure_phasor(83.0, rotated)
            .unwrap(),
        reference,
        2e-13,
    );
    let offset = [20.0, -30.0, 80.0];
    let translated = UniformFlowCompactObserver::new(
        fluid(1.2),
        flow,
        core::array::from_fn(|i| position[i] + offset[i]),
        offset,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    complex_close(
        translated.pressure_phasor(83.0, input).unwrap(),
        reference,
        2e-13,
    );
    let heavy = UniformFlowCompactObserver::new(
        fluid(12.0),
        flow,
        position,
        [0.0; 3],
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let ordinary = observer(flow, position);
    let force = CompactSourcePhasor::force_dipole(input.force_on_fluid_n);
    complex_close(
        heavy.pressure_phasor(83.0, force).unwrap(),
        ordinary.pressure_phasor(83.0, force).unwrap(),
        2e-14,
    );
    let volume = CompactSourcePhasor::monopole(input.volume_flow_m3_s);
    complex_close(
        heavy.pressure_phasor(83.0, volume).unwrap(),
        ordinary.pressure_phasor(83.0, volume).unwrap().scale(10.0),
        2e-14,
    );
}

#[test]
fn near_sonic_downstream_delay_remains_finite_and_clearance_is_physical() {
    let speed = 0.999999 * 300.0;
    let wind = observer([speed, 0.0, 0.0], [2.0, 0.0, 0.0]);
    close(wind.propagation_delay_s(), 2.0 / (300.0 + speed), 3e-15);
    let limits = MovingCompactLimits {
        maximum_mach: 0.6,
        minimum_distance_m: 1.0,
    };
    // Acoustic rest-frame distance is 1/1.5 m; actual laboratory clearance
    // is exactly 1 m and must be accepted at its declared physical boundary.
    assert!(
        UniformFlowCompactObserver::new(
            fluid(1.2),
            [150.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            limits
        )
        .is_ok()
    );
    assert!(
        UniformFlowCompactObserver::new(
            fluid(1.2),
            [150.0, 0.0, 0.0],
            [0.999, 0.0, 0.0],
            [0.0; 3],
            limits
        )
        .is_err()
    );
}

#[test]
fn invalid_flow_history_and_unrepresentable_observations_refuse() {
    let limits = MovingCompactLimits {
        maximum_mach: 0.5,
        minimum_distance_m: 0.1,
    };
    for flow in [[150.001, 0.0, 0.0], [300.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
        assert!(
            UniformFlowCompactObserver::new(fluid(1.2), flow, [2.0, 0.0, 0.0], [0.0; 3], limits)
                .is_err()
        );
    }
    let wind = observer([100.0, 0.0, 0.0], [2.0, 3.0, 1.0]);
    assert!(wind.prepare_receiver(f64::NAN, 1e-12).is_err());
    assert!(wind.prepare_receiver(1e30, 1e-12).is_err());
    assert!(wind.prepare_receiver(0.2, -1.0).is_err());
    assert!(wind.prepare_receiver(0.2, 0.0).is_err());
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(wind.pressure_phasor(frequency, source()).is_err());
    }
    assert!(
        wind.pressure_phasor(
            100.0,
            CompactSourcePhasor::monopole(C64::new(f64::NAN, 0.0))
        )
        .is_err()
    );
    let sample = wind.prepare_receiver(0.2, 1e-15).unwrap();
    assert!(
        sample
            .evaluate(CompactSourceStrength {
                volume_flow_rate_m3_s2: f64::INFINITY,
                ..CompactSourceStrength::default()
            })
            .is_err()
    );
    // An unsuccessful observation never poisons the prepared physical owner.
    assert!(wind.pressure_phasor(100.0, source()).is_ok());
}
