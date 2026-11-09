//! G1/G3: total shaft-speed tangents and actual Omega^2-source resonance.

use fs_ga::{Mat3, So3, Vec3};
use fs_math::c64::C64;
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::resonance::{ResonanceControls, ResonanceSearchError, RotorObservation};
use fs_mbd::rotordynamics::rigid::{RadialBearing, RigidRotor, RotorUnbalance};
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn complex(actual: C64, expected: C64, tolerance: f64) {
    near(actual.re, expected.re, tolerance);
    near(actual.im, expected.im, tolerance);
}

fn times_i(value: C64) -> C64 {
    C64::new(-value.im, value.re)
}

fn rotor(bearing_damping: f64) -> RigidRotor {
    let inertia = SpatialInertia::new(
        4.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap();
    let bearings = [-0.5, 0.5].map(|station_m| RadialBearing {
        station_m,
        stiffness_n_per_m: [[18.0, 0.0], [0.0, 18.0]],
        damping_ns_per_m: [[bearing_damping, 0.0], [0.0, bearing_damping]],
    });
    RigidRotor::new(inertia, So3::identity(), &bearings, 1e-12, || false).unwrap()
}

fn source(station_m: f64) -> [RotorUnbalance; 1] {
    [RotorUnbalance {
        station_m,
        mass_eccentricity_kg_m: 0.02,
        phase_rad: 0.35,
    }]
}

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 2e-11,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: true,
    }
}

#[test]
fn g1_synchronous_tangent_keeps_source_order_and_both_gyroscopic_terms() {
    let (mass, damping, stiffness, gyro) = (1.25, 0.3, 6.0, 0.4);
    let model = RotorBearingModel::new(
        2,
        &[mass, 0.0, 0.0, mass],
        &[damping, 0.0, 0.0, damping],
        &[stiffness, 0.0, 0.0, stiffness],
        &[0.0, -gyro, gyro, 0.0],
        || false,
    )
    .unwrap();
    let unbalance = C64::new(0.025, 0.02);
    for order in [0.5, 1.0, 2.0] {
        for spin in [-2.0, -0.4, 0.0, 0.7, 2.3] {
            let omega = order * spin;
            let f = unbalance.scale(spin * spin);
            let df = unbalance.scale(2.0 * spin);
            let result = model
                .synchronous_response_with_derivative(
                    spin,
                    order,
                    &[f, -times_i(f)],
                    &[df, -times_i(df)],
                    1e-12,
                    || false,
                )
                .unwrap();
            // G*[1,-i] = i*g*[1,-i]. With omega=order*Omega,
            // the circular scalar denominator is k-(m*order^2+g*order)*Omega^2+i*c*order*Omega.
            let inertia = mass * order * order + gyro * order;
            let denominator = C64::new(stiffness - inertia * spin * spin, damping * omega);
            let slope = C64::new(-2.0 * inertia * spin, damping * order);
            let q = f / denominator;
            let dq = (df - slope * q) / denominator;
            let v = times_i(q.scale(omega));
            let dv = times_i(dq.scale(omega) + q.scale(order));
            for (actual, expected) in [
                (&result.response.displacement, [q, -times_i(q)]),
                (&result.displacement_derivative, [dq, -times_i(dq)]),
                (&result.response.velocity, [v, -times_i(v)]),
                (&result.velocity_derivative, [dv, -times_i(dv)]),
            ] {
                for (&value, reference) in actual.iter().zip(expected) {
                    complex(value, reference, 2e-12);
                }
            }
            assert!(result.derivative_relative_residual <= 1e-12);
        }
    }
}

#[test]
fn g1_rigid_unbalance_retains_phase_station_moments_reverse_and_zero_speed() {
    let body = rotor(2.0);
    let unbalances = source(0.25);
    let phase = C64::new(0.35_f64.cos(), 0.35_f64.sin()).scale(0.02);
    for spin in [-3.3, -1.1, 0.0, 0.4, 2.8, 3.3] {
        let result = body
            .unbalance_response_with_derivative(spin, &unbalances, 1e-12, || false)
            .unwrap();
        // Translation: m=4, k=36, c=4. Tilt: Jd=2, Jp=3, k=9, c=1.
        // The rotating station moment selects k+(Jp-Jd)*Omega^2+i*c*Omega.
        let translation = C64::new(36.0 - 4.0 * spin * spin, 4.0 * spin);
        let translation_slope = C64::new(-8.0 * spin, 4.0);
        let tilt = C64::new(9.0 + spin * spin, spin);
        let tilt_slope = C64::new(2.0 * spin, 1.0);
        let force = phase.scale(spin * spin);
        let force_slope = phase.scale(2.0 * spin);
        let x = force / translation;
        let dx = (force_slope - translation_slope * x) / translation;
        let angle_y = force.scale(0.25) / tilt;
        let dangle_y = (force_slope.scale(0.25) - tilt_slope * angle_y) / tilt;
        let q = [x, -times_i(x), times_i(angle_y), angle_y];
        let dq = [dx, -times_i(dx), times_i(dangle_y), dangle_y];
        for i in 0..4 {
            complex(result.response.displacement[i], q[i], 2e-13);
            complex(result.displacement_derivative[i], dq[i], 2e-13);
            complex(
                result.velocity_derivative[i],
                times_i(q[i] + dq[i].scale(spin)),
                2e-13,
            );
        }
        if spin == 0.0 {
            assert!(result.response.displacement.iter().all(|&q| q == C64::ZERO));
            assert!(
                result
                    .displacement_derivative
                    .iter()
                    .all(|&q| q == C64::ZERO)
            );
        }
    }
}

#[test]
fn g1_unbalance_peak_and_half_power_follow_the_actual_speed_dependent_force() {
    let body = rotor(2.0);
    let unbalances = source(0.0);
    let observation = RotorObservation::Displacement(&[1.0, 0.0, 0.0, 0.0]);
    let result = body
        .unbalance_resonance(
            (1.0, 8.0),
            &unbalances,
            observation,
            1e-12,
            controls(),
            || false,
        )
        .unwrap();
    let (mass, damping, stiffness) = (4.0_f64, 4.0_f64, 36.0_f64);
    // |q|^2 = U^2*x^2 / (m^2*x^2+(c^2-2mk)*x+k^2), x=Omega^2.
    let x_peak = 2.0 * stiffness * stiffness / (2.0 * mass * stiffness - damping * damping);
    let peak = x_peak.sqrt();
    near(result.peak.angular_frequency_rad_s, peak, 2e-11);
    let fixed_force_peak = (stiffness / mass - damping * damping / (2.0 * mass * mass)).sqrt();
    assert!(peak > (stiffness / mass).sqrt());
    assert!(fixed_force_peak < (stiffness / mass).sqrt());
    assert!(result.peak.angular_frequency_rad_s - fixed_force_peak > 0.1);
    let d_peak = (stiffness - mass * x_peak).powi(2) + damping * damping * x_peak;
    // The two actual half-height crossings solve this quadratic in x.
    let a = mass * mass - 2.0 * d_peak / (x_peak * x_peak);
    let b = damping * damping - 2.0 * mass * stiffness;
    let discriminant = (b * b - 4.0 * a * stiffness * stiffness).sqrt();
    let low = ((-b - discriminant) / (2.0 * a)).sqrt();
    let high = ((-b + discriminant) / (2.0 * a)).sqrt();
    let band = result.half_power.unwrap();
    near(band.lower.angular_frequency_rad_s, low, 2e-11);
    near(band.upper.angular_frequency_rad_s, high, 2e-11);
    near(band.bandwidth_rad_s, high - low, 4e-11);
    near(band.response_quality_factor, peak / (high - low), 2e-10);
    let native = body
        .unbalance_response_with_derivative(
            result.peak.angular_frequency_rad_s,
            &unbalances,
            1e-12,
            || false,
        )
        .unwrap();
    complex(result.observation, native.response.displacement[0], 0.0);
    complex(
        result.observation_derivative,
        native.displacement_derivative[0],
        0.0,
    );
    let repeated = body
        .unbalance_resonance(
            (1.0, 8.0),
            &unbalances,
            observation,
            1e-12,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result, repeated);
    let scaled = [RotorUnbalance {
        mass_eccentricity_kg_m: 0.14,
        ..unbalances[0]
    }];
    let scaled_result = body
        .unbalance_resonance((1.0, 8.0), &scaled, observation, 1e-12, controls(), || {
            false
        })
        .unwrap();
    near(scaled_result.peak.angular_frequency_rad_s, peak, 2e-11);
    complex(
        scaled_result.observation,
        result.observation.scale(7.0),
        2e-11,
    );
}

#[test]
fn g0_g4_source_singularity_input_budget_and_cancellation_refuse() {
    let body = rotor(2.0);
    let unbalances = source(0.0);
    let observation = RotorObservation::Displacement(&[1.0, 0.0, 0.0, 0.0]);
    let mut limited = controls();
    limited.maximum_evaluations = 3;
    assert!(matches!(
        body.unbalance_resonance((1.0, 8.0), &unbalances, observation, 1e-12, limited, || {
            false
        },),
        Err(ResonanceSearchError::EvaluationLimit { maximum: 3 })
    ));
    assert!(matches!(
        body.unbalance_resonance((1.0, 8.0), &[], observation, 1e-12, controls(), || false,),
        Err(ResonanceSearchError::NotBracketed(_))
    ));
    assert!(matches!(
        rotor(0.0).unbalance_resonance(
            (2.0, 4.0),
            &unbalances,
            observation,
            1e-12,
            controls(),
            || false,
        ),
        Err(ResonanceSearchError::Source(RotorError::Linear(_)))
    ));
    for order in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            body.model().synchronous_response_with_derivative(
                2.0,
                order,
                &[C64::ZERO; 4],
                &[C64::ZERO; 4],
                1e-12,
                || false,
            ),
            Err(RotorError::Invalid(_))
        ));
    }
    assert!(matches!(
        body.model().synchronous_response_with_derivative(
            2.0,
            1.0,
            &[C64::ZERO; 4],
            &[],
            1e-12,
            || false,
        ),
        Err(RotorError::Invalid(_))
    ));
    for stop_after in [0, 75] {
        let mut polls = 0;
        let result = body.unbalance_resonance(
            (1.0, 8.0),
            &unbalances,
            observation,
            1e-12,
            controls(),
            || {
                polls += 1;
                polls > stop_after
            },
        );
        assert!(matches!(
            result,
            Err(ResonanceSearchError::Cancelled
                | ResonanceSearchError::Source(RotorError::Cancelled))
        ));
    }
}
