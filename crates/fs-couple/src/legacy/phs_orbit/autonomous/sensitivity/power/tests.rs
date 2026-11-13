//! G1/G3: original power with complete feedback and explicit memory chains.

use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn original_feedback_supply_and_resistance_keep_distinct_complete_covectors() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let supply = prepared
        .mean_supplied_power(objective_controls(), &jacobian, || false)
        .unwrap();
    let loss = prepared
        .mean_resistance_power(objective_controls(), || false)
        .unwrap();
    let expected = cycle.energy_and_cycle_loss().1 * orbit.omega / core::f64::consts::TAU;
    close(supply.value, expected, 2e-13);
    close(loss.value, expected, 2e-13);
    assert_eq!(supply.angular_frequency_partial, 0.0);
    assert_eq!(loss.angular_frequency_partial, 0.0);
    // The rotational feedback cancels from u.y exactly, leaving a nonlinear
    // radial supply plus an original-state/memory cross term. This independent
    // formula includes the total feedback chain once, even though supply and
    // resistance have the same value on the exact cycle.
    let q2 = (cycle.g * cycle.k).powi(2);
    let a = q2 * cycle.radius_squared();
    let radial = 2.0 * q2 * (cycle.gain - 2.0 * cycle.saturation * a);
    let cross = cycle.mu * cycle.h * cycle.alpha;
    let mut supplied = vec![0.0; 28];
    let mut resistance = vec![0.0; 28];
    for n in 1..=3 {
        for state in 0..4 {
            let z = orbit.coeffs[state][n];
            let dissipation = if state < 2 {
                cycle.r * cycle.k * cycle.k
            } else {
                cycle.eta * cycle.mu * cycle.mu
            };
            let ps = if state < 2 {
                z.scale(radial) + orbit.coeffs[state + 2][n].scale(cross)
            } else {
                orbit.coeffs[state - 2][n].scale(cross)
            };
            let re = 4 + 8 * (n - 1) + state;
            let im = re + 4;
            supplied[re] = 2.0 * ps.re;
            supplied[im] = 2.0 * ps.im;
            resistance[re] = 4.0 * dissipation * z.re;
            resistance[im] = 4.0 * dissipation * z.im;
        }
    }
    for (a, b) in supply.gradient.iter().zip(supplied) {
        close(*a, b, 4e-13);
    }
    for (a, b) in loss.gradient.iter().zip(resistance) {
        close(*a, b, 4e-13);
    }
    // The explicit Hessian path applies the same original complete storage law.
    let hessian = |_: &[f64], v: &[f64], out: &mut [f64]| {
        for i in 0..4 {
            out[i] = v[i] * if i < 2 { cycle.k } else { cycle.mu };
        }
        true
    };
    let analytic = prepared
        .mean_supplied_power_analytic(objective_controls(), &jacobian, &hessian, || false)
        .unwrap();
    close(analytic.value, supply.value, 2e-14);
    for (a, b) in analytic.gradient.iter().zip(&supply.gradient) {
        close(*a, *b, 2e-14);
    }
}

#[test]
fn power_and_per_cycle_energy_adjoints_follow_independent_changed_feedback_cycles() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let slope = Cycle {
        j: 0.0,
        r: 0.0,
        g: 0.0,
        k: 0.0,
        eta: 0.0,
        mu: 0.0,
        h: 0.0,
        gain: 1.0,
        saturation: 0.0,
        shear: 0.0,
        beta: 0.0,
        alpha: 0.0,
    };
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let direction = |_: &[f64], flows: &[f64], _: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[flows[0], flows[1], 0.0, 0.0]);
        true
    };
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let supply = prepared
        .mean_supplied_power(objective_controls(), &jacobian, || false)
        .unwrap();
    let loss = prepared
        .mean_resistance_power(objective_controls(), || false)
        .unwrap();
    let epsilon = 1e-5;
    let lower = cycle.shifted(&slope, -epsilon);
    let upper = cycle.shifted(&slope, epsilon);
    let cycle_derivative =
        (upper.energy_and_cycle_loss().1 - lower.energy_and_cycle_loss().1) / (2.0 * epsilon);
    let power_derivative = (upper.energy_and_cycle_loss().1 * upper.frequency()
        - lower.energy_and_cycle_loss().1 * lower.frequency())
        / (core::f64::consts::TAU * 2.0 * epsilon);
    let source_partial = (cycle.g * cycle.k).powi(2) * cycle.radius_squared();
    for per_cycle in [false, true] {
        for (objective, direct) in [(&supply, source_partial), (&loss, 0.0)] {
            let scale = if per_cycle {
                core::f64::consts::TAU / orbit.omega
            } else {
                1.0
            };
            let gradient: Vec<_> = objective.gradient.iter().map(|g| g * scale).collect();
            let frequency = if per_cycle {
                -scale * objective.value / orbit.omega
            } else {
                0.0
            };
            let mut adjoint = prepared
                .prepare_adjoint(&gradient, frequency, 1.0, || false)
                .unwrap();
            let result = adjoint
                .gradient_parameter(
                    1.0,
                    scale * direct,
                    PhsAutonomousParameterDirection {
                        structure: Default::default(),
                        storage: None,
                        feedback: &direction,
                    },
                    || false,
                )
                .unwrap();
            close(
                result.derivative,
                if per_cycle {
                    cycle_derivative
                } else {
                    power_derivative
                },
                3e-8,
            );
        }
    }
}

#[test]
fn incomplete_power_tangents_cancellation_and_retry_preserve_the_cycle() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let baseline = prepared
        .mean_supplied_power(objective_controls(), &jacobian, || false)
        .unwrap();
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(matches!(
        prepared.mean_supplied_power(objective_controls(), &incomplete, || false),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::ObjectiveRefused { .. }
        ))
    ));
    assert!(
        prepared
            .mean_resistance_power_analytic(objective_controls(), &incomplete, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        prepared.mean_supplied_power(objective_controls(), &jacobian, || {
            polls += 1;
            polls > 11
        }),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    assert_eq!(
        baseline,
        prepared
            .mean_supplied_power(objective_controls(), &jacobian, || false)
            .unwrap()
    );
    assert!(
        prepared
            .mean_resistance_power(objective_controls(), || false)
            .unwrap()
            .value
            > 0.0
    );
}
