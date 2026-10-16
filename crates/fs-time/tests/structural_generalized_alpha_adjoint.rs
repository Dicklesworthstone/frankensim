//! G3: structural pullbacks against an independent acceleration solve,
//! forward-mode duals, and fourth-order finite differences.
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    second_order_adjoint::{SecondOrderAdjointConfig, SecondOrderVjp},
};
use std::ops::{Add, Div, Mul, Sub};

struct Model([f64; 2]);

impl Model {
    fn mass(&self) -> [f64; 4] {
        let [k, d] = self.0;
        [1.4 + 0.04 * k, 0.2 + 0.03 * d, -0.05 * k, 1.7 + 0.04 * d]
    }

    fn damping(&self) -> [f64; 4] {
        let [k, d] = self.0;
        [0.3 + 0.4 * d, 0.1 + 0.02 * k, -0.08 * d, 0.25 + 0.3 * d]
    }

    fn tangent(&self, q: &[f64]) -> [f64; 4] {
        let k = self.0[0];
        [
            2.0 + k + 0.24 * k * q[0],
            0.7,
            -0.4 + 0.08 * k * q[1],
            1.0 + 0.3 * k + 0.08 * k * q[0],
        ]
    }
}

fn multiply(a: [f64; 4], x: &[f64], y: &mut [f64]) {
    y[0] = a[0] * x[0] + a[1] * x[1];
    y[1] = a[2] * x[0] + a[3] * x[1];
}

fn transpose(a: [f64; 4], x: &[f64], y: &mut [f64]) {
    multiply([a[0], a[2], a[1], a[3]], x, y);
}

impl SecondOrderProblem for Model {
    fn dimension(&self) -> usize {
        2
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        multiply(self.mass(), input, output);
    }

    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        multiply(self.damping(), input, output);
    }

    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        let k = self.0[0];
        output[0] = (2.0 + k) * q[0] + 0.7 * q[1] + 0.12 * k * q[0] * q[0];
        output[1] = -0.4 * q[0] + (1.0 + 0.3 * k) * q[1] + 0.08 * k * q[0] * q[1];
    }

    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        multiply(self.tangent(q), direction, output);
    }
}

impl SecondOrderVjp for Model {
    fn parameter_count(&self) -> usize {
        2
    }

    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        transpose(self.mass(), seed, output);
        Ok(())
    }

    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        transpose(self.damping(), seed, output);
        Ok(())
    }

    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        transpose(self.tangent(q), seed, output);
        Ok(())
    }

    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        output[0] = (0.04 * a[0] + 0.02 * v[1] + q[0] + 0.12 * q[0] * q[0]) * seed[0]
            + (-0.05 * a[0] + 0.3 * q[1] + 0.08 * q[0] * q[1]) * seed[1];
        output[1] = (0.03 * a[1] + 0.4 * v[0]) * seed[0]
            + (0.04 * a[1] - 0.08 * v[0] + 0.3 * v[1]) * seed[1];
        Ok(())
    }
}

struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _iteration: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}

const H: f64 = 0.15;

fn method(rho: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        2,
        H,
        rho,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-12,
                relative_tolerance: 2e-13,
                linear_restart: 2,
                max_linear_cycles: 4,
                forcing_maximum: 1e-3,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 12,
        },
    )
}

fn config() -> SecondOrderAdjointConfig {
    SecondOrderAdjointConfig {
        restart: 2,
        max_cycles: 4,
        tolerance: 1e-13,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

// Test-only forward arithmetic. The independent reference differentiates
// converged direct 2x2 Newton iterations with endpoint ACCELERATION as unknown;
// the production pullback differentiates the converged displacement residual.
#[derive(Clone, Copy)]
struct Dual {
    value: f64,
    derivative: [f64; 10],
}

fn c(value: f64) -> Dual {
    Dual {
        value,
        derivative: [0.0; 10],
    }
}

impl Add for Dual {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
            derivative: std::array::from_fn(|i| self.derivative[i] + rhs.derivative[i]),
        }
    }
}

impl Sub for Dual {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            value: self.value - rhs.value,
            derivative: std::array::from_fn(|i| self.derivative[i] - rhs.derivative[i]),
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)] // Forward-mode product rule.
impl Mul for Dual {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
            derivative: std::array::from_fn(|i| {
                self.derivative[i] * rhs.value + self.value * rhs.derivative[i]
            }),
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)] // Forward-mode quotient rule.
impl Div for Dual {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        Self {
            value: self.value / rhs.value,
            derivative: std::array::from_fn(|i| {
                (self.derivative[i] * rhs.value - self.value * rhs.derivative[i])
                    / (rhs.value * rhs.value)
            }),
        }
    }
}

fn variables(input: [f64; 10]) -> [Dual; 10] {
    std::array::from_fn(|i| {
        let mut result = c(input[i]);
        result.derivative[i] = 1.0;
        result
    })
}

#[allow(clippy::many_single_char_names)] // Two-component algebraic reference.
fn oracle_step(input: &[Dual; 10], rho: f64) -> [Dual; 6] {
    let [x, y, vx, vy, ax, ay, fx, fy, k, d] = *input;
    let am = (2.0 * rho - 1.0) / (rho + 1.0);
    let af = rho / (rho + 1.0);
    let gamma = 0.5 - am + af;
    let beta = 0.25 * (1.0 - am + af).powi(2);
    let mass = [
        c(1.4) + c(0.04) * k,
        c(0.2) + c(0.03) * d,
        c(-0.05) * k,
        c(1.7) + c(0.04) * d,
    ];
    let damping = [
        c(0.3) + c(0.4) * d,
        c(0.1) + c(0.02) * k,
        c(-0.08) * d,
        c(0.25) + c(0.3) * d,
    ];
    let mut acceleration = [ax, ay];
    let endpoint = |a: [Dual; 2]| {
        [
            x + c(H) * vx + c(H * H) * (c(0.5 - beta) * ax + c(beta) * a[0]),
            y + c(H) * vy + c(H * H) * (c(0.5 - beta) * ay + c(beta) * a[1]),
            vx + c(H) * (c(1.0 - gamma) * ax + c(gamma) * a[0]),
            vy + c(H) * (c(1.0 - gamma) * ay + c(gamma) * a[1]),
        ]
    };
    for iteration in 0..10 {
        let end = endpoint(acceleration);
        let q = [
            c(af) * x + c(1.0 - af) * end[0],
            c(af) * y + c(1.0 - af) * end[1],
        ];
        let v = [
            c(af) * vx + c(1.0 - af) * end[2],
            c(af) * vy + c(1.0 - af) * end[3],
        ];
        let a = [
            c(am) * ax + c(1.0 - am) * acceleration[0],
            c(am) * ay + c(1.0 - am) * acceleration[1],
        ];
        let residual = [
            mass[0] * a[0]
                + mass[1] * a[1]
                + damping[0] * v[0]
                + damping[1] * v[1]
                + (c(2.0) + k) * q[0]
                + c(0.7) * q[1]
                + c(0.12) * k * q[0] * q[0]
                - fx,
            mass[2] * a[0] + mass[3] * a[1] + damping[2] * v[0] + damping[3] * v[1] - c(0.4) * q[0]
                + (c(1.0) + c(0.3) * k) * q[1]
                + c(0.08) * k * q[0] * q[1]
                - fy,
        ];
        if iteration == 9 {
            assert!(residual.iter().all(|r| r.value.abs() < 1e-12));
        }
        let tangent = [
            c(2.0) + k + c(0.24) * k * q[0],
            c(0.7),
            c(-0.4) + c(0.08) * k * q[1],
            c(1.0) + c(0.3) * k + c(0.08) * k * q[0],
        ];
        let jacobian: [Dual; 4] = std::array::from_fn(|i| {
            c(1.0 - am) * mass[i]
                + c((1.0 - af) * gamma * H) * damping[i]
                + c((1.0 - af) * beta * H * H) * tangent[i]
        });
        let det = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
        acceleration[0] =
            acceleration[0] - (jacobian[3] * residual[0] - jacobian[1] * residual[1]) / det;
        acceleration[1] =
            acceleration[1] - (jacobian[0] * residual[1] - jacobian[2] * residual[0]) / det;
    }
    let end = endpoint(acceleration);
    [
        end[0],
        end[1],
        end[2],
        end[3],
        acceleration[0],
        acceleration[1],
    ]
}

fn finite_difference(input: [f64; 10], index: usize, objective: impl Fn([f64; 10]) -> f64) -> f64 {
    let delta = 2e-4 * (1.0 + input[index].abs());
    let values = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
        let mut perturbed = input;
        perturbed[index] += offset * delta;
        objective(perturbed)
    });
    (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * delta)
}

#[test]
fn structural_step_matches_duals_and_five_point_differences_for_all_inputs() {
    let mut seed = 0x09a7_15c6_b927_458du64;
    let mut random = || {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        ((seed >> 11) as f64) / ((1u64 << 53) as f64)
    };
    for rho in [0.0, 0.2, 0.35, 0.6, 0.85, 1.0] {
        let input = [
            0.6 + random(),
            -0.7 + random(),
            random() - 0.5,
            random() - 0.5,
            random() - 0.5,
            random() - 0.5,
            random(),
            random(),
            0.6 + random(),
            0.2 + random(),
        ];
        let terminal: [f64; 6] = std::array::from_fn(|_| random() - 0.5);
        let method = method(rho);
        let initial = SecondOrderState::new(0.2, &input[..2], &input[2..4], &input[4..6]);
        let gradient = method
            .step_vjp(
                &initial,
                &Model([input[8], input[9]]),
                &input[6..8],
                (&terminal[..2], &terminal[2..4], &terminal[4..]),
                &Identity,
                config(),
                method.adjoint_workspace_components(2, config()).unwrap(),
                &mut || false,
            )
            .unwrap();
        assert!(gradient.primal.newton.converged);
        assert!(gradient.adjoint.converged);
        assert!(gradient.adjoint.euclidean_rel_residual().unwrap() < 1e-13);
        let endpoint = oracle_step(&variables(input), rho);
        let dual = (0..6).fold(c(0.0), |sum, i| sum + c(terminal[i]) * endpoint[i]);
        for (actual, expected) in gradient
            .q
            .iter()
            .chain(&gradient.v)
            .chain(&gradient.a)
            .zip(endpoint)
        {
            close(*actual, expected.value, 3e-11);
        }
        let actual = [
            gradient.initial_q.as_slice(),
            &gradient.initial_v,
            &gradient.initial_a,
            &gradient.forcing,
            &gradient.parameters,
        ]
        .concat();
        for i in 0..10 {
            close(actual[i], dual.derivative[i], 3e-10);
            let difference = finite_difference(input, i, |x| {
                let end = oracle_step(&x.map(c), rho);
                (0..6).map(|j| terminal[j] * end[j].value).sum()
            });
            close(actual[i], difference, 2e-9);
        }
    }
}

fn trajectory(input: &[Dual; 10], rho: f64, steps: usize) -> [Dual; 6] {
    let mut state = *input;
    for _ in 0..steps {
        let next = oracle_step(&state, rho);
        state[..6].copy_from_slice(&next);
    }
    [state[0], state[1], state[2], state[3], state[4], state[5]]
}

#[test]
fn structural_trajectory_chains_q_v_a_and_descends_stiffness_and_damping_loss() {
    let input = [0.9, -0.4, 0.2, -0.1, -0.1, 0.3, 0.3, 0.8, 0.8, 0.6];
    let rho = 0.35;
    let steps = 9;
    let weights = [1.0, 1.0, 0.2, 0.2, 0.02, 0.02];
    let mut truth = input;
    truth[8] = 1.8;
    truth[9] = 0.25;
    let target = trajectory(&truth.map(c), rho, steps).map(|x| x.value);
    let loss = |x: [Dual; 10]| {
        let end = trajectory(&x, rho, steps);
        (0..6).fold(c(0.0), |sum, i| {
            let error = end[i] - c(target[i]);
            sum + c(0.5 * weights[i]) * error * error
        })
    };
    let dual = loss(variables(input));
    let method = method(rho);
    let model = Model([input[8], input[9]]);
    let mut state = SecondOrderState::new(0.2, &input[..2], &input[2..4], &input[4..6]);
    let mut checkpoints = Vec::new();
    for _ in 0..steps {
        checkpoints.push(state.clone());
        method.step(&mut state, &model, &input[6..8]).unwrap();
    }
    let mut seed_q: Vec<f64> = (0..2)
        .map(|i| weights[i] * (state.q[i] - target[i]))
        .collect();
    let mut seed_v: Vec<f64> = (0..2)
        .map(|i| weights[i + 2] * (state.v[i] - target[i + 2]))
        .collect();
    let mut seed_a: Vec<f64> = (0..2)
        .map(|i| weights[i + 4] * (state.a[i] - target[i + 4]))
        .collect();
    let mut parameter_gradient = [0.0; 2];
    let mut forcing_gradient = [0.0; 2];
    for initial in checkpoints.iter().rev() {
        let gradient = method
            .step_vjp(
                initial,
                &model,
                &input[6..8],
                (&seed_q, &seed_v, &seed_a),
                &Identity,
                config(),
                method.adjoint_workspace_components(2, config()).unwrap(),
                &mut || false,
            )
            .unwrap();
        seed_q = gradient.initial_q;
        seed_v = gradient.initial_v;
        seed_a = gradient.initial_a;
        for i in 0..2 {
            parameter_gradient[i] += gradient.parameters[i];
            forcing_gradient[i] += gradient.forcing[i];
        }
    }
    let actual = [
        seed_q.as_slice(),
        seed_v.as_slice(),
        seed_a.as_slice(),
        &forcing_gradient,
        &parameter_gradient,
    ]
    .concat();
    for i in 0..10 {
        close(actual[i], dual.derivative[i], 3e-10);
        close(
            actual[i],
            finite_difference(input, i, |x| loss(x.map(c)).value),
            2e-9,
        );
    }
    let mut candidate = input;
    for i in 0..2 {
        assert!(parameter_gradient[i].abs() > 1e-5);
        candidate[8 + i] -= parameter_gradient[i];
    }
    assert!(loss(candidate.map(c)).value < 0.999 * dual.value);
}
