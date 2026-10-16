//! G3: discrete pullbacks against independent dense Newton, forward duals,
//! and fourth-order finite differences, including a transient inverse step.
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_time::galpha::{
    FirstOrderProblem, FirstOrderState, ImplicitSolveConfig, OperatorFirstOrderGeneralizedAlpha,
    adjoint::{FirstOrderAdjointConfig, FirstOrderVjp},
};
use std::ops::{Add, Div, Mul, Sub};

struct Model([f64; 2]);

impl Model {
    fn mass(&self) -> [f64; 4] {
        let [p, q] = self.0;
        [1.4 + 0.2 * p, 0.25 + 0.1 * q, -0.1 * p, 1.8 + 0.3 * q]
    }

    fn tangent(&self, t: f64, u: &[f64]) -> [f64; 4] {
        let [p, q] = self.0;
        [
            2.0 + p + 0.3 * q * u[0],
            0.7 + 0.05 * t,
            -0.4 + 0.1 * p * u[1] + 0.03 * t,
            1.0 + q + 0.1 * p * u[0],
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

impl FirstOrderProblem for Model {
    fn dimension(&self) -> usize {
        2
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        multiply(self.mass(), input, output);
    }

    fn internal_force(&self, t: f64, u: &[f64], output: &mut [f64]) {
        let [p, q] = self.0;
        output[0] = (2.0 + p) * u[0] + (0.7 + 0.05 * t) * u[1] + 0.15 * q * u[0] * u[0];
        output[1] = (-0.4 + 0.03 * t) * u[0] + (1.0 + q) * u[1] + 0.1 * p * u[0] * u[1];
    }

    fn tangent_apply(&self, t: f64, u: &[f64], direction: &[f64], output: &mut [f64]) {
        multiply(self.tangent(t, u), direction, output);
    }
}

impl FirstOrderVjp for Model {
    fn parameter_count(&self) -> usize {
        2
    }

    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        transpose(self.mass(), seed, output);
        Ok(())
    }

    fn tangent_transpose_apply(
        &self,
        t: f64,
        u: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        transpose(self.tangent(t, u), seed, output);
        Ok(())
    }

    fn residual_parameter_vjp(
        &self,
        _t: f64,
        u: &[f64],
        rate: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        output[0] =
            (0.2 * rate[0] + u[0]) * seed[0] + (-0.1 * rate[0] + 0.1 * u[0] * u[1]) * seed[1];
        output[1] =
            (0.1 * rate[1] + 0.15 * u[0] * u[0]) * seed[0] + (0.3 * rate[1] + u[1]) * seed[1];
        Ok(())
    }
}

struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _iteration: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}

const H: f64 = 0.08;

fn method(rho: f64) -> OperatorFirstOrderGeneralizedAlpha {
    OperatorFirstOrderGeneralizedAlpha::new(
        2,
        H,
        rho,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-13,
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

fn config() -> FirstOrderAdjointConfig {
    FirstOrderAdjointConfig {
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

// Test-only forward-mode arithmetic. The oracle has its own residual and
// solves for the endpoint RATE by a direct 2x2 inverse; production uses
// matrix-free Newton--Krylov with endpoint STATE as its unknown. Differentiating
// these converged reference Newton iterations is only an independent oracle.
#[derive(Clone, Copy)]
struct Dual {
    value: f64,
    derivative: [f64; 8],
}

fn c(value: f64) -> Dual {
    Dual {
        value,
        derivative: [0.0; 8],
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

fn variables(input: [f64; 8]) -> [Dual; 8] {
    std::array::from_fn(|i| {
        let mut result = c(input[i]);
        result.derivative[i] = 1.0;
        result
    })
}

#[allow(clippy::many_single_char_names)] // Two-component algebraic reference.
fn oracle_step(input: &[Dual; 8], t: f64, rho: f64) -> [Dual; 4] {
    let [x, y, vx, vy, fx, fy, p, q] = *input;
    let am = (3.0 - rho) / (2.0 * (1.0 + rho));
    let af = 1.0 / (1.0 + rho);
    let gamma = 0.5 + am - af;
    let time = t + af * H;
    let mass = [
        c(1.4) + c(0.2) * p,
        c(0.25) + c(0.1) * q,
        c(-0.1) * p,
        c(1.8) + c(0.3) * q,
    ];
    let mut rate = [vx, vy];
    let endpoint = |rate: [Dual; 2]| {
        [
            x + c(H) * (c(1.0 - gamma) * vx + c(gamma) * rate[0]),
            y + c(H) * (c(1.0 - gamma) * vy + c(gamma) * rate[1]),
        ]
    };
    for iteration in 0..10 {
        let u = endpoint(rate);
        let a = c(1.0 - af) * x + c(af) * u[0];
        let b = c(1.0 - af) * y + c(af) * u[1];
        let da = c(1.0 - am) * vx + c(am) * rate[0];
        let db = c(1.0 - am) * vy + c(am) * rate[1];
        let residual = [
            mass[0] * da
                + mass[1] * db
                + (c(2.0) + p) * a
                + c(0.7 + 0.05 * time) * b
                + c(0.15) * q * a * a
                - fx,
            mass[2] * da
                + mass[3] * db
                + c(-0.4 + 0.03 * time) * a
                + (c(1.0) + q) * b
                + c(0.1) * p * a * b
                - fy,
        ];
        if iteration == 9 {
            assert!(residual.iter().all(|r| r.value.abs() < 1e-12));
        }
        let tangent = [
            c(2.0) + p + c(0.3) * q * a,
            c(0.7 + 0.05 * time),
            c(-0.4 + 0.03 * time) + c(0.1) * p * b,
            c(1.0) + q + c(0.1) * p * a,
        ];
        let j: [Dual; 4] =
            std::array::from_fn(|i| c(am) * mass[i] + c(af * gamma * H) * tangent[i]);
        let det = j[0] * j[3] - j[1] * j[2];
        rate[0] = rate[0] - (j[3] * residual[0] - j[1] * residual[1]) / det;
        rate[1] = rate[1] - (j[0] * residual[1] - j[2] * residual[0]) / det;
    }
    let u = endpoint(rate);
    [u[0], u[1], rate[0], rate[1]]
}

fn finite_difference(input: [f64; 8], index: usize, objective: impl Fn([f64; 8]) -> f64) -> f64 {
    let delta = 2e-4 * (1.0 + input[index].abs());
    let values = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
        let mut perturbed = input;
        perturbed[index] += offset * delta;
        objective(perturbed)
    });
    (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * delta)
}

#[test]
fn step_cotangents_match_forward_duals_and_five_point_differences() {
    let mut seed = 0x78a4_920d_f29c_390bu64;
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
            random(),
            random(),
            0.6 + random(),
            0.3 + random(),
        ];
        let terminal = [
            random() - 0.5,
            random() - 0.5,
            random() - 0.5,
            random() - 0.5,
        ];
        let t = 0.1 + random();
        let method = method(rho);
        let initial = FirstOrderState::new(t, &input[..2], &input[2..4]);
        let gradient = method
            .step_vjp(
                &initial,
                &Model([input[6], input[7]]),
                &input[4..6],
                (&terminal[..2], &terminal[2..]),
                &Identity,
                config(),
                method.adjoint_workspace_components(2, config()).unwrap(),
                &mut || false,
            )
            .unwrap();
        assert!(gradient.primal.newton.converged);
        assert!(gradient.adjoint.converged);
        assert!(gradient.adjoint.euclidean_rel_residual().unwrap() < 1e-13);
        let endpoint = oracle_step(&variables(input), t, rho);
        let dual = (0..4).fold(c(0.0), |sum, i| sum + c(terminal[i]) * endpoint[i]);
        let actual = [
            &gradient.initial[..],
            &gradient.initial_rate[..],
            &gradient.forcing[..],
            &gradient.parameters[..],
        ]
        .concat();
        for (actual, expected) in gradient.value.iter().chain(&gradient.rate).zip(endpoint) {
            close(*actual, expected.value, 2e-12);
        }
        for i in 0..8 {
            close(actual[i], dual.derivative[i], 3e-11);
            let difference = finite_difference(input, i, |x| {
                let end = oracle_step(&x.map(c), t, rho);
                (0..4).map(|j| terminal[j] * end[j].value).sum()
            });
            close(actual[i], difference, 2e-9);
        }
    }
}

fn trajectory(input: &[Dual; 8], rho: f64, steps: usize) -> [Dual; 4] {
    let mut state = *input;
    for step in 0..steps {
        let next = oracle_step(&state, 0.2 + step as f64 * H, rho);
        state[..4].copy_from_slice(&next);
    }
    [state[0], state[1], state[2], state[3]]
}

#[test]
fn trajectory_reverses_both_state_vectors_and_reduces_inverse_loss() {
    let input = [0.9, -0.4, 0.2, -0.1, 0.3, 0.8, 0.8, 1.4];
    let rho = 0.35;
    let steps = 7;
    let mut truth = input;
    truth[6] = 1.8;
    truth[7] = 0.45;
    let target = trajectory(&truth.map(c), rho, steps).map(|x| x.value);
    let loss = |x: [Dual; 8]| {
        let end = trajectory(&x, rho, steps);
        (0..4).fold(c(0.0), |sum, i| {
            let error = end[i] - c(target[i]);
            sum + c(if i < 2 { 0.5 } else { 0.025 }) * error * error
        })
    };
    let dual = loss(variables(input));
    let method = method(rho);
    let model = Model([input[6], input[7]]);
    let mut state = FirstOrderState::new(0.2, &input[..2], &input[2..4]);
    let mut checkpoints = Vec::new();
    for _ in 0..steps {
        checkpoints.push(state.clone());
        method.step(&mut state, &model, &input[4..6]).unwrap();
    }
    let mut seed_u: Vec<f64> = (0..2).map(|i| state.u[i] - target[i]).collect();
    let mut seed_rate: Vec<f64> = (0..2)
        .map(|i| 0.05 * (state.rate[i] - target[i + 2]))
        .collect();
    let mut parameter_gradient = [0.0; 2];
    let mut forcing_gradient = [0.0; 2];
    for initial in checkpoints.iter().rev() {
        let gradient = method
            .step_vjp(
                initial,
                &model,
                &input[4..6],
                (&seed_u, &seed_rate),
                &Identity,
                config(),
                method.adjoint_workspace_components(2, config()).unwrap(),
                &mut || false,
            )
            .unwrap();
        seed_u = gradient.initial;
        seed_rate = gradient.initial_rate;
        for i in 0..2 {
            parameter_gradient[i] += gradient.parameters[i];
            forcing_gradient[i] += gradient.forcing[i];
        }
    }
    let actual = [
        seed_u.as_slice(),
        seed_rate.as_slice(),
        &forcing_gradient,
        &parameter_gradient,
    ]
    .concat();
    for i in 0..8 {
        close(actual[i], dual.derivative[i], 3e-11);
        close(
            actual[i],
            finite_difference(input, i, |x| loss(x.map(c)).value),
            2e-9,
        );
    }
    let mut candidate = input;
    for i in 0..2 {
        candidate[6 + i] -= parameter_gradient[i];
    }
    assert!(loss(candidate.map(c)).value < 0.999 * dual.value);
}
