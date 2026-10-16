//! G3/G4: consistent initialization, its IFT pullback, and the full step chain.
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_time::galpha::{
    FirstOrderGeneralizedAlpha, FirstOrderProblem, FirstOrderState, GeneralizedAlpha,
    ImplicitSolveConfig, OperatorFirstOrderGeneralizedAlpha, OperatorGeneralizedAlpha,
    SecondOrderProblem, SecondOrderState,
    adjoint::{FirstOrderAdjointConfig, FirstOrderVjp},
    first_order_galpha_step, galpha_step,
    initialization::{
        FirstOrderInitialGradient, InitialSolveConfig, InitialSolveError,
        SecondOrderInitialGradient, first_order_rate, first_order_rate_vjp,
        initial_workspace_components, second_order_acceleration, second_order_acceleration_vjp,
    },
    second_order_adjoint::{SecondOrderAdjointConfig, SecondOrderVjp},
};

const TIME: f64 = 0.37;
const H: f64 = 0.15;
const RHO: f64 = 0.35;

#[derive(Clone, Copy)]
struct Model {
    p: [f64; 2],
    fault: u8,
}
impl Model {
    fn new(p: [f64; 2]) -> Self {
        Self { p, fault: 0 }
    }
    fn mass(&self) -> [f64; 4] {
        let [p, q] = self.p;
        [1.4 + 0.2 * p, 0.25 + 0.1 * q, -0.1 * p, 1.8 + 0.3 * q]
    }
    fn damping(&self) -> [f64; 4] {
        let [p, q] = self.p;
        [0.3 + 0.4 * q, 0.1 + 0.02 * p, -0.08 * q, 0.25 + 0.3 * q]
    }
    fn stiffness(&self) -> [f64; 4] {
        [2.0 + self.p[0], 0.7, -0.4, 1.0 + self.p[1]]
    }
    fn bias(&self, t: f64) -> [f64; 2] {
        [0.03 * t * self.p[0], -0.05 * t * self.p[1]]
    }
    fn mass_vjp(&self, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        if self.fault == 1 {
            return Err("mass transpose unavailable".into());
        }
        if self.fault != 2 {
            out.copy_from_slice(&transpose(self.mass(), seed));
        }
        Ok(())
    }
    fn tangent_vjp(&self, seed: &[f64], out: &mut [f64]) {
        if self.fault != 3 {
            out.copy_from_slice(&transpose(self.stiffness(), seed));
        }
    }
}

fn multiply(a: [f64; 4], x: &[f64]) -> [f64; 2] {
    [a[0] * x[0] + a[1] * x[1], a[2] * x[0] + a[3] * x[1]]
}
fn transpose(a: [f64; 4], x: &[f64]) -> [f64; 2] {
    multiply([a[0], a[2], a[1], a[3]], x)
}
// Independent direct 2x2 solve; no production Krylov or adjoint implementation.
fn dense_solve(a: [f64; 4], rhs: [f64; 2]) -> [f64; 2] {
    let det = a[0] * a[3] - a[1] * a[2];
    [
        (a[3] * rhs[0] - a[1] * rhs[1]) / det,
        (a[0] * rhs[1] - a[2] * rhs[0]) / det,
    ]
}

impl FirstOrderProblem for Model {
    fn dimension(&self) -> usize {
        2
    }
    fn mass_apply(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&multiply(self.mass(), x));
    }
    fn internal_force(&self, t: f64, u: &[f64], out: &mut [f64]) {
        let linear = multiply(self.stiffness(), u);
        let bias = self.bias(t);
        for i in 0..2 {
            out[i] = linear[i] + bias[i];
        }
    }
    fn tangent_apply(&self, _t: f64, _u: &[f64], x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&multiply(self.stiffness(), x));
    }
}
impl FirstOrderVjp for Model {
    fn parameter_count(&self) -> usize {
        2
    }
    fn mass_transpose_apply(&self, s: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.mass_vjp(s, out)
    }
    fn tangent_transpose_apply(
        &self,
        _t: f64,
        _u: &[f64],
        s: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        self.tangent_vjp(s, out);
        Ok(())
    }
    fn residual_parameter_vjp(
        &self,
        t: f64,
        u: &[f64],
        rate: &[f64],
        s: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        if self.fault == 4 {
            return Err("parameter derivative unavailable".into());
        }
        if self.fault != 5 {
            out[0] = (0.2 * rate[0] + u[0] + 0.03 * t) * s[0] - 0.1 * rate[0] * s[1];
            out[1] = 0.1 * rate[1] * s[0] + (0.3 * rate[1] + u[1] - 0.05 * t) * s[1];
        }
        Ok(())
    }
}
impl SecondOrderProblem for Model {
    fn dimension(&self) -> usize {
        2
    }
    fn mass_apply(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&multiply(self.mass(), x));
    }
    fn damping_apply(&self, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&multiply(self.damping(), x));
    }
    fn internal_force(&self, q: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&multiply(self.stiffness(), q));
    }
    fn tangent_apply(&self, _q: &[f64], x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&multiply(self.stiffness(), x));
    }
}
impl SecondOrderVjp for Model {
    fn parameter_count(&self) -> usize {
        2
    }
    fn mass_transpose_apply(&self, s: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.mass_vjp(s, out)
    }
    fn damping_transpose_apply(&self, s: &[f64], out: &mut [f64]) -> Result<(), String> {
        if self.fault != 6 {
            out.copy_from_slice(&transpose(self.damping(), s));
        }
        Ok(())
    }
    fn tangent_transpose_apply(
        &self,
        _q: &[f64],
        s: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        self.tangent_vjp(s, out);
        Ok(())
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        s: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        if self.fault == 4 {
            return Err("parameter derivative unavailable".into());
        }
        if self.fault != 5 {
            out[0] = (0.2 * a[0] + 0.02 * v[1] + q[0]) * s[0] - 0.1 * a[0] * s[1];
            out[1] = (0.1 * a[1] + 0.4 * v[0]) * s[0]
                + (0.3 * a[1] - 0.08 * v[0] + 0.3 * v[1] + q[1]) * s[1];
        }
        Ok(())
    }
}

struct Preconditioner(f64);
impl FlexiblePreconditioner for Preconditioner {
    fn apply(&self, _iteration: usize, x: &[f64], out: &mut [f64]) {
        for (y, x) in out.iter_mut().zip(x) {
            *y = self.0 * x;
        }
    }
}
fn config() -> InitialSolveConfig {
    InitialSolveConfig {
        restart: 2,
        max_cycles: 4,
        tolerance: 1e-13,
    }
}
fn cap() -> usize {
    initial_workspace_components(2, 2, config()).unwrap()
}
fn step_config() -> ImplicitSolveConfig {
    ImplicitSolveConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-13,
            linear_restart: 2,
            max_linear_cycles: 4,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 8,
    }
}
fn first_pull<C: FnMut() -> bool>(
    model: &Model,
    x: &[f64; 6],
    seed: &[f64],
    cancel: &mut C,
) -> Result<FirstOrderInitialGradient, InitialSolveError> {
    first_order_rate_vjp(
        model,
        TIME,
        &x[..2],
        &x[2..4],
        seed,
        &Preconditioner(1.0),
        &Preconditioner(1.0),
        config(),
        cap(),
        cancel,
    )
}
fn second_pull<C: FnMut() -> bool>(
    model: &Model,
    x: &[f64; 8],
    seed: &[f64],
    cancel: &mut C,
) -> Result<SecondOrderInitialGradient, InitialSolveError> {
    second_order_acceleration_vjp(
        model,
        &x[..2],
        &x[2..4],
        &x[4..6],
        seed,
        &Preconditioner(1.0),
        &Preconditioner(1.0),
        config(),
        cap(),
        cancel,
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() < tol * (1.0 + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn check_derivatives<const N: usize>(
    actual: &[f64],
    x: [f64; N],
    objective: impl Fn([f64; N]) -> f64,
) {
    assert_eq!(actual.len(), N);
    for i in 0..N {
        let h = 2e-4 * (1.0 + x[i].abs());
        let values = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
            let mut point = x;
            point[i] += h * offset;
            objective(point)
        });
        close(
            actual[i],
            (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * h),
            2e-8,
        );
    }
}
fn first_dense(x: [f64; 6]) -> [f64; 2] {
    let model = Model::new([x[4], x[5]]);
    let internal = multiply(model.stiffness(), &x[..2]);
    let bias = model.bias(TIME);
    dense_solve(
        model.mass(),
        [x[2] - internal[0] - bias[0], x[3] - internal[1] - bias[1]],
    )
}
fn second_dense(x: [f64; 8]) -> [f64; 2] {
    let model = Model::new([x[6], x[7]]);
    let internal = multiply(model.stiffness(), &x[..2]);
    let damping = multiply(model.damping(), &x[2..4]);
    dense_solve(
        model.mass(),
        [
            x[4] - internal[0] - damping[0],
            x[5] - internal[1] - damping[1],
        ],
    )
}

#[test]
fn consistent_rates_and_accelerations_match_dense_solves_and_all_partials() {
    let first = [0.9, -0.4, 0.3, 0.8, 1.2, 0.65];
    let second = [0.9, -0.4, 0.2, -0.1, 0.3, 0.8, 1.2, 0.65];
    let model = Model::new([1.2, 0.65]);
    let seed = [0.7, -0.4];
    let rate = first_order_rate(
        &model,
        TIME,
        &first[..2],
        &first[2..4],
        &Preconditioner(1.0),
        config(),
        cap(),
        &mut || false,
    )
    .unwrap();
    let acceleration = second_order_acceleration(
        &model,
        &second[..2],
        &second[2..4],
        &second[4..6],
        &Preconditioner(1.0),
        config(),
        cap(),
        &mut || false,
    )
    .unwrap();
    assert!(rate.report.converged && acceleration.report.converged);
    for i in 0..2 {
        close(rate.value[i], first_dense(first)[i], 1e-13);
        close(acceleration.value[i], second_dense(second)[i], 1e-13);
    }
    let mass_rate = multiply(model.mass(), &rate.value);
    let mass_a = multiply(model.mass(), &acceleration.value);
    let internal = multiply(model.stiffness(), &first[..2]);
    let damping = multiply(model.damping(), &second[2..4]);
    for i in 0..2 {
        close(
            mass_rate[i] + internal[i] + model.bias(TIME)[i],
            first[2 + i],
            1e-13,
        );
        close(mass_a[i] + damping[i] + internal[i], second[4 + i], 1e-13);
    }
    let first_gradient = first_pull(&model, &first, &seed, &mut || false).unwrap();
    let second_gradient = second_pull(&model, &second, &seed, &mut || false).unwrap();
    assert_eq!(first_gradient.rate, rate.value);
    assert_eq!(second_gradient.acceleration, acceleration.value);
    assert!(first_gradient.primal.converged && first_gradient.adjoint.converged);
    assert!(second_gradient.primal.converged && second_gradient.adjoint.converged);
    check_derivatives(
        &[
            first_gradient.initial,
            first_gradient.forcing,
            first_gradient.parameters,
        ]
        .concat(),
        first,
        |x| dot(&first_dense(x), &seed),
    );
    check_derivatives(
        &[
            second_gradient.initial_q,
            second_gradient.initial_v,
            second_gradient.forcing,
            second_gradient.parameters,
        ]
        .concat(),
        second,
        |x| dot(&second_dense(x), &seed),
    );
}

fn first_endpoint(x: [f64; 6]) -> [f64; 4] {
    let model = Model::new([x[4], x[5]]);
    let dense = FirstOrderGeneralizedAlpha::new(&model.mass(), &model.stiffness(), 2, H, RHO);
    let mut u = [x[0], x[1]];
    let mut rate = first_dense(x);
    let bias = model.bias(TIME + H / (1.0 + RHO));
    first_order_galpha_step(&dense, &mut u, &mut rate, &[x[2] - bias[0], x[3] - bias[1]]);
    [u[0], u[1], rate[0], rate[1]]
}
fn second_endpoint(x: [f64; 8]) -> [f64; 6] {
    let model = Model::new([x[6], x[7]]);
    let dense = GeneralizedAlpha::new(
        &model.mass(),
        &model.damping(),
        &model.stiffness(),
        2,
        H,
        RHO,
    );
    let mut q = [x[0], x[1]];
    let mut v = [x[2], x[3]];
    let mut a = second_dense(x);
    galpha_step(&dense, &mut q, &mut v, &mut a, &x[4..6]);
    [q[0], q[1], v[0], v[1], a[0], a[1]]
}

#[test]
fn first_order_initial_consistency_chain_matches_the_complete_dense_step() {
    let x = [0.9, -0.4, 0.3, 0.8, 1.2, 0.65];
    let model = Model::new([x[4], x[5]]);
    let initial = first_order_rate(
        &model,
        TIME,
        &x[..2],
        &x[2..4],
        &Preconditioner(1.0),
        config(),
        cap(),
        &mut || false,
    )
    .unwrap();
    let state = FirstOrderState::new(TIME, &x[..2], &initial.value);
    let method = OperatorFirstOrderGeneralizedAlpha::new(2, H, RHO, step_config());
    let solve = FirstOrderAdjointConfig {
        restart: 2,
        max_cycles: 4,
        tolerance: 1e-13,
    };
    let seed = [0.3, -0.7, 0.2, 0.5];
    let step = method
        .step_vjp(
            &state,
            &model,
            &x[2..4],
            (&seed[..2], &seed[2..]),
            &Preconditioner(1.0),
            solve,
            method.adjoint_workspace_components(2, solve).unwrap(),
            &mut || false,
        )
        .unwrap();
    for (actual, expected) in step.value.iter().chain(&step.rate).zip(first_endpoint(x)) {
        close(*actual, expected, 1e-12);
    }
    let chain = first_pull(&model, &x, &step.initial_rate, &mut || false).unwrap();
    let total = [
        std::array::from_fn::<_, 2, _>(|i| step.initial[i] + chain.initial[i]),
        std::array::from_fn(|i| step.forcing[i] + chain.forcing[i]),
        std::array::from_fn(|i| step.parameters[i] + chain.parameters[i]),
    ]
    .concat();
    check_derivatives(&total, x, |point| dot(&first_endpoint(point), &seed));
    assert!(
        chain.parameters.iter().any(|v| v.abs() > 1e-4),
        "omitting initialization must change the gradient"
    );
}

#[test]
fn structural_initial_consistency_chain_matches_the_complete_dense_step() {
    let x = [0.9, -0.4, 0.2, -0.1, 0.3, 0.8, 1.2, 0.65];
    let model = Model::new([x[6], x[7]]);
    let initial = second_order_acceleration(
        &model,
        &x[..2],
        &x[2..4],
        &x[4..6],
        &Preconditioner(1.0),
        config(),
        cap(),
        &mut || false,
    )
    .unwrap();
    let state = SecondOrderState::new(TIME, &x[..2], &x[2..4], &initial.value);
    let method = OperatorGeneralizedAlpha::new(2, H, RHO, step_config());
    let solve = SecondOrderAdjointConfig {
        restart: 2,
        max_cycles: 4,
        tolerance: 1e-13,
    };
    let seed = [0.3, -0.7, 0.2, 0.5, -0.15, 0.25];
    let step = method
        .step_vjp(
            &state,
            &model,
            &x[4..6],
            (&seed[..2], &seed[2..4], &seed[4..]),
            &Preconditioner(1.0),
            solve,
            method.adjoint_workspace_components(2, solve).unwrap(),
            &mut || false,
        )
        .unwrap();
    for (actual, expected) in step
        .q
        .iter()
        .chain(&step.v)
        .chain(&step.a)
        .zip(second_endpoint(x))
    {
        close(*actual, expected, 1e-12);
    }
    let chain = second_pull(&model, &x, &step.initial_a, &mut || false).unwrap();
    let total = [
        std::array::from_fn::<_, 2, _>(|i| step.initial_q[i] + chain.initial_q[i]),
        std::array::from_fn(|i| step.initial_v[i] + chain.initial_v[i]),
        std::array::from_fn(|i| step.forcing[i] + chain.forcing[i]),
        std::array::from_fn(|i| step.parameters[i] + chain.parameters[i]),
    ]
    .concat();
    check_derivatives(&total, x, |point| dot(&second_endpoint(point), &seed));
    assert!(
        chain.parameters.iter().any(|v| v.abs() > 1e-4),
        "omitting initialization must change the gradient"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Compare both initialization families against the same refusals.
fn initialization_refuses_missing_derivatives_bad_inputs_and_failed_solves() {
    let x = [0.9, -0.4, 0.3, 0.8, 1.2, 0.65];
    let y = [0.9, -0.4, 0.2, -0.1, 0.3, 0.8, 1.2, 0.65];
    let model = Model::new([1.2, 0.65]);
    let seed = [0.7, -0.4];
    for mode in 1..=5 {
        let faulty = Model {
            fault: mode,
            ..model
        };
        let first = first_pull(&faulty, &x, &seed, &mut || false);
        let second = second_pull(&faulty, &y, &seed, &mut || false);
        if mode == 1 || mode == 4 {
            assert!(matches!(first, Err(InitialSolveError::Derivative(_))));
            assert!(matches!(second, Err(InitialSolveError::Derivative(_))));
        } else {
            assert!(matches!(first, Err(InitialSolveError::NonFiniteOutput)));
            assert!(matches!(second, Err(InitialSolveError::NonFiniteOutput)));
        }
    }
    assert!(matches!(
        second_pull(&Model { fault: 6, ..model }, &y, &seed, &mut || false),
        Err(InitialSolveError::NonFiniteOutput)
    ));
    assert!(initial_workspace_components(usize::MAX, 2, config()).is_none());
    assert!(matches!(
        first_order_rate_vjp(
            &model,
            TIME,
            &x[..2],
            &x[2..4],
            &seed,
            &Preconditioner(1.0),
            &Preconditioner(1.0),
            config(),
            cap() - 1,
            &mut || false
        ),
        Err(InitialSolveError::WorkspaceLimit { .. })
    ));
    assert!(matches!(
        first_order_rate(
            &model,
            TIME,
            &[0.9],
            &x[2..4],
            &Preconditioner(1.0),
            config(),
            cap(),
            &mut || false
        ),
        Err(InitialSolveError::InvalidInput(_))
    ));
    assert!(matches!(
        first_order_rate(
            &model,
            TIME,
            &x[..2],
            &x[2..4],
            &Preconditioner(1.0),
            InitialSolveConfig {
                restart: 0,
                ..config()
            },
            cap(),
            &mut || false
        ),
        Err(InitialSolveError::InvalidInput(_))
    ));
    assert!(matches!(
        first_order_rate(
            &model,
            TIME,
            &x[..2],
            &x[2..4],
            &Preconditioner(0.0),
            config(),
            cap(),
            &mut || false
        ),
        Err(InitialSolveError::NotConverged(_))
    ));
    assert!(matches!(
        first_order_rate_vjp(
            &model,
            TIME,
            &x[..2],
            &x[2..4],
            &seed,
            &Preconditioner(1.0),
            &Preconditioner(0.0),
            config(),
            cap(),
            &mut || false
        ),
        Err(InitialSolveError::NotConverged(_))
    ));
    assert!(matches!(
        second_order_acceleration(
            &model,
            &y[..2],
            &y[2..4],
            &[f64::NAN, 0.8],
            &Preconditioner(1.0),
            config(),
            cap(),
            &mut || false
        ),
        Err(InitialSolveError::InvalidInput(_))
    ));
    let zero = first_pull(&model, &x, &[0.0; 2], &mut || false).unwrap();
    assert_eq!(zero.initial, [0.0; 2]);
    assert_eq!(zero.forcing, [0.0; 2]);
    assert_eq!(zero.parameters, [0.0; 2]);
    assert!(zero.adjoint.converged);
}

#[test]
fn every_initialization_cancellation_boundary_refuses_and_retries() {
    let x = [0.9, -0.4, 0.3, 0.8, 1.2, 0.65];
    let y = [0.9, -0.4, 0.2, -0.1, 0.3, 0.8, 1.2, 0.65];
    let model = Model::new([1.2, 0.65]);
    let seed = [0.7, -0.4];
    let mut first_polls = 0;
    let first = first_pull(&model, &x, &seed, &mut || {
        first_polls += 1;
        false
    })
    .unwrap();
    let mut second_polls = 0;
    let second = second_pull(&model, &y, &seed, &mut || {
        second_polls += 1;
        false
    })
    .unwrap();
    for stop in 1..=first_polls {
        let mut count = 0;
        assert!(
            matches!(
                first_pull(&model, &x, &seed, &mut || {
                    count += 1;
                    count == stop
                }),
                Err(InitialSolveError::Cancelled)
            ),
            "first-order cancellation boundary {stop}"
        );
    }
    for stop in 1..=second_polls {
        let mut count = 0;
        assert!(
            matches!(
                second_pull(&model, &y, &seed, &mut || {
                    count += 1;
                    count == stop
                }),
                Err(InitialSolveError::Cancelled)
            ),
            "second-order cancellation boundary {stop}"
        );
    }
    let retry_first = first_pull(&model, &x, &seed, &mut || false).unwrap();
    let retry_second = second_pull(&model, &y, &seed, &mut || false).unwrap();
    assert_eq!(retry_first.rate, first.rate);
    assert_eq!(retry_first.initial, first.initial);
    assert_eq!(retry_first.parameters, first.parameters);
    assert_eq!(retry_first.forcing, first.forcing);
    assert_eq!(retry_second.acceleration, second.acceleration);
    assert_eq!(retry_second.initial_q, second.initial_q);
    assert_eq!(retry_second.initial_v, second.initial_v);
    assert_eq!(retry_second.parameters, second.parameters);
    assert_eq!(retry_second.forcing, second.forcing);
}
