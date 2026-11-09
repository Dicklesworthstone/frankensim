//! Independent shaft-speed oracles with original hereditary bearing memory.

use fs_ga::{Mat3, So3, Vec3};
use fs_math::c64::C64;
use fs_mbd::{
    articulated::SpatialInertia,
    rotordynamics::{
        RotorBearingModel, RotorError,
        relaxation::{MaxwellDampingBranch, ResonanceControls, ResonanceSearchError},
        resonance::RotorObservation,
        rigid::{RadialBearing, RigidRotor, RotorUnbalance},
    },
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    near(actual.re, expected.re, tolerance);
    near(actual.im, expected.im, tolerance);
}
fn times_i(value: C64) -> C64 {
    C64::new(-value.im, value.re)
}
fn branch(strain: &[f64], k: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: strain.to_vec(),
        stiffness: k,
        relaxation_time_s: tau,
    }
}
fn transfer(w: f64, k: f64, tau: f64) -> (C64, C64) {
    let denominator = C64::new(1.0, w * tau);
    (
        C64::new(0.0, k * w * tau) / denominator,
        C64::new(0.0, k * tau) / (denominator * denominator),
    )
}
fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 2e-11,
        relative_frequency_tolerance: 1e-12,
        maximum_evaluations: 512,
        measure_half_power_bandwidth: true,
    }
}
fn rotor(k: f64, c: f64) -> RigidRotor {
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
        stiffness_n_per_m: [[k, 0.0], [0.0, k]],
        damping_ns_per_m: [[c, 0.0], [0.0, c]],
    });
    RigidRotor::new(inertia, So3::identity(), &bearings, 1e-12, || false).unwrap()
}
fn source(station: f64) -> [RotorUnbalance; 1] {
    [RotorUnbalance {
        station_m: station,
        mass_eccentricity_kg_m: 0.02,
        phase_rad: 0.35,
    }]
}

#[test]
fn signed_order_and_spin_include_gyro_source_and_each_relaxation_chain_rule() {
    let (m, c, k, g) = (1.25, 0.3, 6.0, 0.4);
    let model = RotorBearingModel::new(
        2,
        &[m, 0.0, 0.0, m],
        &[c, 0.0, 0.0, c],
        &[k, 0.0, 0.0, k],
        &[0.0, -g, g, 0.0],
        || false,
    )
    .unwrap();
    let branches = [branch(&[1.0, 0.0], 2.0, 0.6), branch(&[0.0, 1.0], 2.0, 0.6)];
    let source = C64::new(0.025, 0.02);
    for order in [-2.0, -0.5, 0.5, 1.0, 2.0] {
        for spin in [-2.0, -0.4, 0.0, 0.7, 2.3] {
            let point = model
                .at_speed_with_maxwell_memory(spin, &branches, || false)
                .unwrap();
            let w = order * spin;
            let force = source.scale(spin * spin);
            let force_slope = source.scale(2.0 * spin);
            let result = point
                .synchronous_response_with_derivative(
                    order,
                    &[force, -times_i(force)],
                    &[force_slope, -times_i(force_slope)],
                    3e-12,
                    || false,
                )
                .unwrap();
            // G*[1,-i]=i*g*[1,-i], so both moving omega and Omega contribute.
            let inertia = m * order * order + g * order;
            let (h, dh) = transfer(w, 2.0, 0.6);
            let d = C64::new(k - inertia * spin * spin, c * w) + h;
            let dd = C64::new(-2.0 * inertia * spin, c * order) + dh.scale(order);
            let q = force / d;
            let dq = (force_slope - dd * q) / d;
            for (i, phase) in [C64::ONE, C64::new(0.0, -1.0)].into_iter().enumerate() {
                close(result.response.motion.displacement[i], phase * q, 3e-12);
                close(result.displacement_derivative[i], phase * dq, 3e-12);
                close(
                    result.response.motion.velocity[i],
                    phase * times_i(q.scale(w)),
                    3e-12,
                );
                close(
                    result.velocity_derivative[i],
                    phase * times_i(dq.scale(w) + q.scale(order)),
                    3e-12,
                );
                close(result.response.branch_stress[i], phase * h * q, 3e-12);
                close(
                    result.branch_stress_derivative[i],
                    phase * (h * dq + dh.scale(order) * q),
                    3e-12,
                );
                near(
                    result.response.branch_mean_heat_w[i],
                    (h * q).norm_sq() / (2.0 * 2.0 * 0.6),
                    3e-12,
                );
            }
            assert_eq!(result.response.motion.spin_rad_s, spin);
            assert_eq!(result.response.motion.frequency_rad_s, w);
            assert!(result.derivative_relative_residual < 3e-12);
            assert!(result.relaxation_derivative_relative_residual < 3e-12);
        }
    }
}

#[test]
fn rigid_unbalance_retains_original_station_phase_moments_heat_and_zero_speed() {
    let body = rotor(18.0, 2.0);
    let unbalances = source(0.25);
    let branches = [
        branch(&[1.0, 0.0, 0.0, 0.0], 6.0, 0.4),
        branch(&[0.0, 1.0, 0.0, 0.0], 6.0, 0.4),
        branch(&[0.0, 0.0, 1.0, 0.0], 2.0, 0.7),
        branch(&[0.0, 0.0, 0.0, 1.0], 2.0, 0.7),
    ];
    let u = C64::new(0.35_f64.cos(), 0.35_f64.sin()).scale(0.02);
    for spin in [-3.3, -1.1, 0.0, 0.4, 2.8, 3.3] {
        let response = body
            .unbalance_maxwell_response_with_derivative(spin, &unbalances, &branches, 3e-12, || {
                false
            })
            .unwrap();
        let (ht, dht) = transfer(spin, 6.0, 0.4);
        let (hr, dhr) = transfer(spin, 2.0, 0.7);
        let d = C64::new(36.0 - 4.0 * spin * spin, 4.0 * spin) + ht;
        let dd = C64::new(-8.0 * spin, 4.0) + dht;
        // The physical bending moment selects Jp-Jd=+1, not a removed gyro term.
        let r = C64::new(9.0 + spin * spin, spin) + hr;
        let dr = C64::new(2.0 * spin, 1.0) + dhr;
        let f = u.scale(spin * spin);
        let df = u.scale(2.0 * spin);
        let x = f / d;
        let dx = (df - dd * x) / d;
        let angle_y = f.scale(0.25) / r;
        let dangle_y = (df.scale(0.25) - dr * angle_y) / r;
        let q = [x, -times_i(x), times_i(angle_y), angle_y];
        let dq = [dx, -times_i(dx), times_i(dangle_y), dangle_y];
        for i in 0..4 {
            let (h, dh, k, tau) = if i < 2 {
                (ht, dht, 6.0, 0.4)
            } else {
                (hr, dhr, 2.0, 0.7)
            };
            close(response.response.motion.displacement[i], q[i], 3e-12);
            close(response.displacement_derivative[i], dq[i], 3e-12);
            close(
                response.velocity_derivative[i],
                times_i(q[i] + dq[i].scale(spin)),
                3e-12,
            );
            close(response.response.branch_stress[i], h * q[i], 3e-12);
            close(
                response.branch_stress_derivative[i],
                h * dq[i] + dh * q[i],
                3e-12,
            );
            let heat = (h * q[i]).norm_sq() / (2.0 * k * tau);
            near(response.response.branch_mean_work_w[i], heat, 3e-12);
            near(response.response.branch_mean_heat_w[i], heat, 3e-12);
        }
        near(response.response.mean_power_residual_w, 0.0, 3e-12);
        if spin == 0.0 {
            assert!(
                response
                    .response
                    .motion
                    .displacement
                    .iter()
                    .chain(&response.displacement_derivative)
                    .chain(&response.response.branch_stress)
                    .chain(&response.branch_stress_derivative)
                    .all(|&z| z == C64::ZERO)
            );
        }
    }
}

// The oracle is a fixed, independently derived real polynomial; it never calls
// a rotor response, a response derivative or the production resonance search.
fn polynomial_root(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    let sign = f(lo).is_sign_positive();
    assert_ne!(sign, f(hi).is_sign_positive());
    for _ in 0..90 {
        let mid = 0.5 * lo + 0.5 * hi;
        if f(mid).is_sign_positive() == sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * lo + 0.5 * hi
}

#[test]
fn hereditary_unbalance_peak_and_width_follow_actual_speed_squared_source() {
    // Translation divided by m=4 has relaxed k=1, Maxwell k=1, tau=1.
    let body = rotor(2.0, 0.0);
    let branches = [
        branch(&[1.0, 0.0, 0.0, 0.0], 4.0, 1.0),
        branch(&[0.0, 1.0, 0.0, 0.0], 4.0, 1.0),
    ];
    let unbalances = source(0.0);
    let observation = RotorObservation::Displacement(&[1.0, 0.0, 0.0, 0.0]);
    let result = body
        .unbalance_maxwell_resonance(
            (0.8, 2.6),
            &unbalances,
            &branches,
            observation,
            3e-12,
            controls(),
            || false,
        )
        .unwrap();
    // y=Omega^2: |q|^2/(U/m)^2 = y^2*(1+y)/(y^3-3y^2+2y+1).
    // Its positive maximum solves 4y^3-4y^2-5y-2=0.
    let y = polynomial_root(1.0, 2.0, |y| 4.0 * y * y * y - 4.0 * y * y - 5.0 * y - 2.0);
    let denominator = |y: f64| y * y * y - 3.0 * y * y + 2.0 * y + 1.0;
    let gain = y * y * (1.0 + y) / denominator(y);
    let half_polynomial = |x: f64| gain * denominator(x) - 2.0 * x * x * (1.0 + x);
    let lower = polynomial_root(0.64, y, half_polynomial).sqrt();
    let upper = polynomial_root(y, 6.76, half_polynomial).sqrt();
    near(result.peak.angular_frequency_rad_s, y.sqrt(), 3e-10);
    near(
        result.observation.norm_sq(),
        (0.02_f64 / 4.0).powi(2) * gain,
        3e-12,
    );
    let band = result.half_power.unwrap();
    near(band.lower.angular_frequency_rad_s, lower, 3e-10);
    near(band.upper.angular_frequency_rad_s, upper, 3e-10);
    near(band.bandwidth_rad_s, upper - lower, 3e-10);
    near(
        band.response_quality_factor,
        y.sqrt() / (upper - lower),
        3e-10,
    );
    let fixed_force_peak = polynomial_root(1.0, 2.0, |y| 2.0 * y * y * y - 6.0 * y + 1.0).sqrt();
    assert!(result.peak.angular_frequency_rad_s - fixed_force_peak > 0.03);
    let physical = body
        .unbalance_maxwell_response_with_derivative(
            result.peak.angular_frequency_rad_s,
            &unbalances,
            &branches,
            3e-12,
            || false,
        )
        .unwrap();
    assert_eq!(result.observation, physical.response.motion.displacement[0]);
    assert_eq!(
        result.observation_derivative,
        physical.displacement_derivative[0]
    );
    assert!(
        physical
            .response
            .branch_mean_heat_w
            .iter()
            .all(|&p| p > 0.0)
    );
    near(physical.response.mean_power_residual_w, 0.0, 3e-12);
    let larger = [RotorUnbalance {
        mass_eccentricity_kg_m: 0.1,
        ..unbalances[0]
    }];
    let scaled = body
        .unbalance_maxwell_resonance(
            (0.8, 2.6),
            &larger,
            &branches,
            observation,
            3e-12,
            controls(),
            || false,
        )
        .unwrap();
    near(
        scaled.peak.angular_frequency_rad_s,
        result.peak.angular_frequency_rad_s,
        3e-10,
    );
    close(scaled.observation, result.observation.scale(5.0), 3e-10);
}

#[test]
fn synchronous_peak_width_and_derivative_are_in_shaft_speed_units_for_signed_orders() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let branches = [branch(&[1.0], 1.0, 1.0)];
    let template = model
        .at_speed_with_maxwell_memory(9.0, &branches, || false)
        .unwrap();
    let frequency_peak = (3.0_f64.sqrt() - 1.0).sqrt();
    let minimum = 2.0 * 3.0_f64.sqrt() - 3.0;
    let level = 2.0 * minimum;
    let delta = ((1.0 + level).powi(2) - 4.0 * (1.0 - level)).sqrt();
    let low = ((1.0 + level - delta) / 2.0).sqrt();
    let high = ((1.0 + level + delta) / 2.0).sqrt();
    for order in [-2.0_f64, -0.5, 0.5, 2.0] {
        let result = template
            .synchronous_resonance(
                (0.05 / order.abs(), 2.0 / order.abs()),
                order,
                &[C64::ONE],
                RotorObservation::Velocity(&[1.0]),
                3e-12,
                controls(),
                || false,
            )
            .unwrap();
        near(
            result.peak.angular_frequency_rad_s,
            frequency_peak / order.abs(),
            3e-10,
        );
        let band = result.half_power.unwrap();
        near(band.lower.angular_frequency_rad_s, low / order.abs(), 3e-10);
        near(
            band.upper.angular_frequency_rad_s,
            high / order.abs(),
            3e-10,
        );
        near(
            band.response_quality_factor,
            frequency_peak / (high - low),
            3e-10,
        );
        let spin = result.peak.angular_frequency_rad_s;
        let w = order * spin;
        let (h, dh) = transfer(w, 1.0, 1.0);
        let d = C64::from_re(-w * w) + h;
        let q = C64::ONE / d;
        let dq = -(C64::from_re(-2.0 * w) + dh).scale(order) * q / d;
        close(result.observation, times_i(q.scale(w)), 3e-10);
        close(
            result.observation_derivative,
            times_i(dq.scale(w) + q.scale(order)),
            3e-10,
        );
    }
}

#[test]
fn bad_orders_sources_singular_limits_and_query_budgets_refuse_atomically() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.2], &[4.0], &[0.0], || false).unwrap();
    let branches = [branch(&[1.0], 1.0, 1.0)];
    let point = model
        .at_speed_with_maxwell_memory(0.5, &branches, || false)
        .unwrap();
    for order in [0.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(
            point
                .synchronous_response_with_derivative(
                    order,
                    &[C64::ONE],
                    &[C64::ZERO],
                    3e-12,
                    || false,
                )
                .is_err()
        );
    }
    assert!(
        point
            .synchronous_response_with_derivative(1.0, &[C64::ONE], &[], 3e-12, || false)
            .is_err()
    );
    assert!(
        point
            .synchronous_resonance(
                (0.5, 3.0),
                1.0,
                &[C64::ONE],
                RotorObservation::Velocity(&[]),
                3e-12,
                controls(),
                || false,
            )
            .is_err()
    );
    let body = rotor(2.0, 0.0);
    let memories = [
        branch(&[1.0, 0.0, 0.0, 0.0], 4.0, 1.0),
        branch(&[0.0, 1.0, 0.0, 0.0], 4.0, 1.0),
    ];
    let source = source(0.0);
    let receiver = RotorObservation::Displacement(&[1.0, 0.0, 0.0, 0.0]);
    let mut limited = controls();
    limited.maximum_evaluations = 3;
    assert!(matches!(
        body.unbalance_maxwell_resonance(
            (0.8, 2.6),
            &source,
            &memories,
            receiver,
            3e-12,
            limited,
            || false
        ),
        Err(ResonanceSearchError::EvaluationLimit { maximum: 3 })
    ));
    assert!(
        body.unbalance_maxwell_resonance(
            (0.8, 2.6),
            &[],
            &memories,
            receiver,
            3e-12,
            controls(),
            || false,
        )
        .is_err()
    );
    for stop_after in [0, 75] {
        let mut polls = 0;
        assert!(matches!(
            body.unbalance_maxwell_resonance(
                (0.8, 2.6),
                &source,
                &memories,
                receiver,
                3e-12,
                controls(),
                || {
                    polls += 1;
                    polls > stop_after
                },
            ),
            Err(ResonanceSearchError::Cancelled
                | ResonanceSearchError::Source(RotorError::Cancelled))
        ));
    }
    let native = body
        .unbalance_response_with_derivative(1.3, &source, 3e-12, || false)
        .unwrap();
    let empty = body
        .unbalance_maxwell_response_with_derivative(1.3, &source, &[], 3e-12, || false)
        .unwrap();
    close(
        empty.response.motion.displacement[0],
        native.response.displacement[0],
        3e-12,
    );
    close(
        empty.displacement_derivative[0],
        native.displacement_derivative[0],
        3e-12,
    );
    let after = body
        .unbalance_maxwell_response_with_derivative(1.3, &source, &memories, 3e-12, || false)
        .unwrap();
    let replay = body
        .unbalance_maxwell_response_with_derivative(1.3, &source, &memories, 3e-12, || false)
        .unwrap();
    assert_eq!(after, replay);
}
