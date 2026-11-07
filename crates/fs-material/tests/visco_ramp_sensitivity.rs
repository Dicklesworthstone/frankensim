//! G0/G1/G3 derivatives of the original transient damping and physical heat.
//! Synthetic constitutive experiments, not measured-material identification.

use fs_material::visco::{
    GeneralizedMaxwell, PronyRamp, PronyRampLinearization, PronyRampSensitivity, PronyState,
};
use fs_math::det;

fn close(actual: f64, expected: f64, relative: f64) {
    let tolerance = relative * expected.abs() + 16.0 * f64::from_bits(1);
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn exact_jacobian_matches_independent_sensitivity_integrals() {
    let (q0, delta) = (0.002, -0.003);
    for r in [0.01, 0.0999999, 0.1, 0.1000001, 0.5, 1.0, 8.0] {
        let prepared = PronyRampSensitivity::new(r);
        let response = prepared.apply(q0, delta);
        let original = PronyRamp::new(r).apply(q0, delta);
        assert_eq!(response.response.0.to_bits(), original.0.to_bits());
        assert_eq!(response.response.1.to_bits(), original.1.to_bits());
        assert_eq!(response.response.2.to_bits(), original.2.to_bits());
        assert_eq!(prepared.ramp().apply(q0, delta), original);

        // Independently differentiate the ODE solution at physical time dt*u,
        // then integrate mean strain and r*q(u)^2. In particular the heat
        // derivative includes d r / d ln(tau) = -r, not only d q / d ln(tau).
        let n = 8192;
        let mut mean_gradient = [0.0; 3];
        let mut heat_gradient = [0.0; 3];
        let mut endpoint_gradient = [0.0; 3];
        for i in 0..=n {
            let u = f64::from(i) / f64::from(n);
            let a = det::exp(-r * u);
            let b = -det::expm1(-r * u) / r;
            let q = a * q0 + b * delta;
            let dq = [a, b, r * u * a * q0 + (b - u * a) * delta];
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for j in 0..3 {
                mean_gradient[j] += weight * dq[j];
                heat_gradient[j] += weight * 2.0 * r * q * dq[j];
            }
            heat_gradient[2] -= weight * r * q * q;
            if i == n {
                endpoint_gradient = dq;
            }
        }
        for j in 0..3 {
            close(response.jacobian[0][j], endpoint_gradient[j], 4e-12);
            close(
                response.jacobian[1][j],
                mean_gradient[j] / (3.0 * f64::from(n)),
                4e-11,
            );
            close(
                response.jacobian[2][j],
                heat_gradient[j] / (3.0 * f64::from(n)),
                4e-11,
            );
        }
    }
}

#[test]
fn weak_loss_retains_relaxation_time_and_heat_derivatives() {
    let delta = 0.002;
    for q0 in [0.0, 0.001, -0.0003] {
        for r in [1e-300, 1e-18, 1e-14, 1e-10, 1e-6] {
            let response = PronyRampSensitivity::new(r).apply(q0, delta);
            let leading_heat = r * (q0 * q0 + q0 * delta + delta * delta / 3.0);
            assert!(response.response.2 > 0.0);
            assert!(response.jacobian[2][2] < 0.0);
            let tolerance = 4.0 * r + 2e-13;
            close(response.response.2, leading_heat, tolerance);
            close(response.jacobian[0][2], r * (q0 + 0.5 * delta), tolerance);
            close(
                response.jacobian[1][2],
                r * (0.5 * q0 + delta / 6.0),
                tolerance,
            );
            close(response.jacobian[2][0], r * (2.0 * q0 + delta), tolerance);
            close(
                response.jacobian[2][1],
                r * (q0 + 2.0 * delta / 3.0),
                tolerance,
            );
            close(response.jacobian[2][2], -leading_heat, tolerance);
        }
    }
}

#[test]
fn elastic_and_stiff_limits_have_physical_derivatives() {
    let (q0, delta) = (0.001, 0.002);
    let instant = PronyRampSensitivity::new(0.0).apply(q0, delta);
    assert_eq!(
        instant.jacobian,
        [[1.0, 1.0, 0.0], [1.0, 0.5, 0.0], [0.0; 3]]
    );
    assert_eq!(instant.response, (q0 + delta, q0 + 0.5 * delta, 0.0));
    for r in [1e4, 1e100, f64::MAX] {
        let response = PronyRampSensitivity::new(r).apply(q0, delta);
        close(response.jacobian[0][2], delta / r, 3e-13);
        close(
            response.jacobian[1][2],
            (q0 + delta * (1.0 - 2.0 / r)) / r,
            3e-13,
        );
        close(
            response.jacobian[2][2],
            (delta * (q0 + delta) - 3.0 * delta * delta / r) / r,
            3e-13,
        );
    }
    let relaxed = PronyRampSensitivity::new(f64::INFINITY).apply(q0, delta);
    assert_eq!(relaxed.response, (0.0, 0.0, 0.5 * q0 * q0));
    assert_eq!(relaxed.jacobian, [[0.0; 3], [0.0; 3], [q0, 0.0, 0.0]]);
}

#[test]
fn tangent_and_pullback_preserve_energy_identity_and_duality() {
    let (q0, delta) = (0.002, -0.003);
    let direction = [0.7, -0.2, 1.3];
    let seed = [-0.4, 0.8, 2.7];
    for r in [0.0, 1e-6, 0.03, 0.1, 1.0, 20.0, 1e6, f64::INFINITY] {
        let response = PronyRampSensitivity::new(r).apply(q0, delta);
        let tangent = response.push_forward(direction);
        let pullback = response.pull_back(seed);
        let dot = |a: [f64; 3], b: [f64; 3]| a.into_iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        close(dot(seed, tangent), dot(direction, pullback), 2e-14);
        for j in 0..3 {
            let work =
                response.jacobian[1][j] * delta + if j == 1 { response.response.1 } else { 0.0 };
            let storage =
                response.response.0 * response.jacobian[0][j] - if j == 0 { q0 } else { 0.0 };
            let residual = work - storage - response.jacobian[2][j];
            assert!(
                residual.abs() < 2e-17,
                "differentiated energy residual {residual:e}"
            );
        }
    }
}

type Trace = (PronyState, Vec<Vec<PronyRampLinearization>>, Vec<f64>);

fn trace(law: &GeneralizedMaxwell, steps: &[(f64, f64)]) -> Trace {
    let mut state = law.state();
    let mut tape = Vec::new();
    let mut stresses = Vec::new();
    let mut heat = 0.0;
    for &(strain, dt) in steps {
        let delta = strain - state.strain;
        let step: Vec<_> = law
            .terms
            .iter()
            .zip(&state.q)
            .map(|(&(_, tau), &q)| PronyRampSensitivity::new(dt / tau).apply(q, delta))
            .collect();
        // The public material owner advances the sole primal history.
        law.step(&mut state, strain, dt);
        for ((response, &q), &(e, _)) in step.iter().zip(&state.q).zip(&law.terms) {
            assert_eq!(response.response.0.to_bits(), q.to_bits());
            heat += e * response.response.2;
        }
        stresses.push(state.stress);
        tape.push(step);
    }
    close(heat, state.dissipated, 3e-15);
    (state, tape, stresses)
}

fn reverse(
    law: &GeneralizedMaxwell,
    tape: &[Vec<PronyRampLinearization>],
    stress_seeds: &[f64],
    heat_seed: f64,
) -> Vec<f64> {
    let mut adjoint = vec![0.0; law.terms.len()];
    let mut gradient = adjoint.clone();
    for (step, &stress_seed) in tape.iter().zip(stress_seeds).rev() {
        for (j, &(e, _)) in law.terms.iter().enumerate() {
            let pulled = step[j].pull_back([adjoint[j] + e * stress_seed, 0.0, e * heat_seed]);
            adjoint[j] = pulled[0];
            gradient[j] += pulled[2];
        }
    }
    gradient
}

#[test]
fn multibranch_loading_reversal_adjoint_keeps_all_retained_memory_and_heat() {
    let law = GeneralizedMaxwell::new(2.0, vec![(3.0, 0.05), (7.0, 0.8)]).unwrap();
    let truth = GeneralizedMaxwell::new(2.0, vec![(3.0, 0.035), (7.0, 1.1)]).unwrap();
    let steps = [
        (0.02, 0.0),
        (0.02, 0.025),
        (-0.01, 0.075),
        (0.005, 0.1),
        (0.005, 0.7),
        (0.015, 0.2),
        (0.0, 0.05),
    ];
    let (_, _, observed) = trace(&truth, &steps);
    let (state, tape, stresses) = trace(&law, &steps);
    let residuals: Vec<_> = stresses.iter().zip(&observed).map(|(s, o)| s - o).collect();
    let heat_weight = 0.2;
    let gradient = reverse(&law, &tape, &residuals, heat_weight);
    let mut forward_q = [0.0; 2];
    let mut forward_heat = [0.0; 2];
    let mut forward_work = [0.0; 2];
    let mut forward_goal = [0.0; 2];
    let mut old_strain = 0.0;
    for ((step, &(strain, _)), &residual) in tape.iter().zip(&steps).zip(&residuals) {
        for (j, &(e, _)) in law.terms.iter().enumerate() {
            let pushed = step[j].push_forward([forward_q[j], 0.0, 1.0]);
            forward_q[j] = pushed[0];
            forward_heat[j] += e * pushed[2];
            forward_work[j] += e * pushed[1] * (strain - old_strain);
            forward_goal[j] += e * (residual * pushed[0] + heat_weight * pushed[2]);
        }
        old_strain = strain;
    }
    let objective = |model: &GeneralizedMaxwell| {
        let (state, _, stresses) = trace(model, &steps);
        0.5 * stresses
            .iter()
            .zip(&observed)
            .map(|(s, o)| (s - o).powi(2))
            .sum::<f64>()
            + heat_weight * state.dissipated
    };
    for j in 0..2 {
        close(gradient[j], forward_goal[j], 3e-13);
        close(
            forward_work[j],
            law.terms[j].0 * state.q[j] * forward_q[j] + forward_heat[j],
            3e-12,
        );
        // Only the regression perturbs parameters; production uses the exact
        // native pullback and never advances a second perturbed history.
        let h = 2e-5;
        let mut plus = law.clone();
        let mut minus = law.clone();
        plus.terms[j].1 *= det::exp(h);
        minus.terms[j].1 *= det::exp(-h);
        let finite_difference = (objective(&plus) - objective(&minus)) / (2.0 * h);
        close(gradient[j], finite_difference, 3e-8);
    }
}

#[test]
fn synthetic_relaxation_time_calibration_uses_the_transient_adjoint() {
    let (equilibrium, branch, strain, authored_tau) = (2e9, 3e9, 0.002, 0.37);
    let times = [0.0, 0.02, 0.06, 0.15, 0.35, 0.7, 1.2];
    let mut previous = 0.0;
    let steps: Vec<_> = times
        .iter()
        .map(|&time| {
            let dt = time - previous;
            previous = time;
            (strain, dt)
        })
        .collect();
    // Independent synthetic relaxation observations: a step strain followed
    // by held strain, with the authored exponential Maxwell relaxation.
    let observed: Vec<_> = times
        .iter()
        .map(|&time| equilibrium * strain + branch * strain * det::exp(-time / authored_tau))
        .collect();
    let scale = branch * strain;
    let mut log_tau = det::ln(1.1);
    let mut initial_cost = 0.0;
    let mut final_cost = 0.0;
    for iteration in 0..12 {
        let tau = det::exp(log_tau);
        let law = GeneralizedMaxwell::new(equilibrium, vec![(branch, tau)]).unwrap();
        let (state, tape, stresses) = trace(&law, &steps);
        let seeds: Vec<_> = stresses
            .iter()
            .zip(&observed)
            .map(|(s, o)| ((s - o) / scale) / scale)
            .collect();
        final_cost = 0.5
            * stresses
                .iter()
                .zip(&observed)
                .map(|(s, o)| ((s - o) / scale).powi(2))
                .sum::<f64>();
        if iteration == 0 {
            initial_cost = final_cost;
        }
        let adjoint = reverse(&law, &tape, &seeds, 0.0)[0];
        let mut curvature = 0.0;
        let mut q_tangent = 0.0;
        let mut forward_gradient = 0.0;
        for (((step, &time), &stress), &observed) in
            tape.iter().zip(&times).zip(&stresses).zip(&observed)
        {
            q_tangent = step[0].push_forward([q_tangent, 0.0, 1.0])[0];
            let normalized = branch * q_tangent / scale;
            let analytic = time / tau * det::exp(-time / tau);
            if time > 0.0 {
                close(normalized, analytic, 3e-13);
            }
            curvature += normalized * normalized;
            forward_gradient += ((stress - observed) / scale) * normalized;
        }
        if final_cost < 1e-26 {
            break;
        }
        close(adjoint, forward_gradient, 3e-12);
        // The physical cumulative heat has an independent analytic derivative
        // after the entire nonuniform hold sequence, including its old memory.
        let heat_gradient = reverse(&law, &tape, &vec![0.0; times.len()], 1.0)[0];
        let r = times[times.len() - 1] / tau;
        close(
            heat_gradient,
            -branch * strain * strain * r * det::exp(-2.0 * r),
            3e-12,
        );
        close(
            state.dissipated,
            -0.5 * branch * strain * strain * det::expm1(-2.0 * r),
            3e-13,
        );
        // A bounded scalar Gauss--Newton step, using the native reverse
        // gradient, identifies the positive relaxation time in log space.
        log_tau -= (adjoint / curvature).clamp(-1.0, 1.0);
    }
    assert!(initial_cost > 0.1);
    assert!(final_cost < 1e-26, "calibration residual {final_cost:e}");
    close(det::exp(log_tau), authored_tau, 2e-12);
}
