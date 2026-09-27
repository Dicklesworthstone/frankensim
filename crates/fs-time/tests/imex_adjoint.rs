//! Gradients of the actual stiff step and a multi-step terminal objective.
use fs_solver::LinearOp;
use fs_time::stiff::{
    IdentityPreconditioner, Imex2, ImexSolveConfig, ImexSolveError, ImexState, OperatorImex2,
    adjoint::{ImexAdjointError, ImexVjp},
    imex2_step,
};
use std::cell::Cell;

struct Model {
    p: [f64; 2],
}
impl Model {
    fn matrix(&self) -> [f64; 4] {
        [-self.p[0], 1.3, 0.4, -self.p[0] - 3.0]
    }
}
impl LinearOp for Model {
    fn n(&self) -> usize {
        2
    }
    fn apply(&self, u: &[f64], out: &mut [f64]) {
        let a = self.matrix();
        out[0] = a[0] * u[0] + a[1] * u[1];
        out[1] = a[2] * u[0] + a[3] * u[1];
    }
    fn apply_transpose(&self, u: &[f64], out: &mut [f64]) {
        let a = self.matrix();
        out[0] = a[0] * u[0] + a[2] * u[1];
        out[1] = a[1] * u[0] + a[3] * u[1];
    }
}
impl ImexVjp for Model {
    fn parameter_count(&self) -> usize {
        2
    }
    fn nonlinear(&self, u: &[f64], out: &mut [f64]) {
        out[0] = self.p[1] * u[0] * u[0] + 0.1 * u[1];
        out[1] = -0.2 * u[0] * u[1] + self.p[1];
    }
    fn nonlinear_vjp(
        &self,
        u: &[f64],
        seed: &[f64],
        ub: &mut [f64],
        pb: &mut [f64],
    ) -> Result<(), String> {
        ub[0] = 2.0 * self.p[1] * u[0] * seed[0] - 0.2 * u[1] * seed[1];
        ub[1] = 0.1 * seed[0] - 0.2 * u[0] * seed[1];
        pb[0] = 0.0;
        pb[1] = u[0] * u[0] * seed[0] + seed[1];
        Ok(())
    }
    fn linear_parameter_vjp(&self, u: &[f64], seed: &[f64], pb: &mut [f64]) -> Result<(), String> {
        pb[0] = -u[0] * seed[0] - u[1] * seed[1];
        pb[1] = 0.0;
        Ok(())
    }
}

fn method() -> OperatorImex2 {
    OperatorImex2::new(
        2,
        0.08,
        ImexSolveConfig {
            tolerance: 1e-13,
            restart: 2,
            max_cycles: 4,
        },
    )
}
fn dense_endpoint(p: [f64; 2], u: &[f64], steps: usize) -> Vec<f64> {
    let model = Model { p };
    let dense = Imex2::new(&model.matrix(), 2, 0.08);
    let mut result = u.to_vec();
    for _ in 0..steps {
        imex2_step(&dense, &mut result, &|x, y| model.nonlinear(x, y));
    }
    result
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-8 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn nonsymmetric_step_matches_dense_differences_in_state_and_both_parameters() {
    let model = Model { p: [4.0, 0.7] };
    let u = [0.9, -0.4];
    let terminal = [0.3, -1.1];
    let method = method();
    let gradient = method
        .step_vjp(
            &u,
            &model,
            &IdentityPreconditioner,
            &IdentityPreconditioner,
            &terminal,
            method.adjoint_workspace_components(2).unwrap(),
            &mut || false,
        )
        .unwrap();
    let mut production = ImexState::new(0.0, &u);
    method
        .step(&mut production, &model, &IdentityPreconditioner, &|x, y| {
            model.nonlinear(x, y)
        })
        .unwrap();
    assert_eq!(production.u, gradient.value);
    for report in gradient.primal.iter().chain(&gradient.adjoint) {
        assert!(report.converged && report.rel_residual < 1e-13);
    }
    let eps = 1e-5;
    for i in 0..2 {
        let (mut up, mut um) = (u, u);
        up[i] += eps;
        um[i] -= eps;
        let fd = (dot(&dense_endpoint(model.p, &up, 1), &terminal)
            - dot(&dense_endpoint(model.p, &um, 1), &terminal))
            / (2.0 * eps);
        close(gradient.initial[i], fd);
        let (mut pp, mut pm) = (model.p, model.p);
        pp[i] += eps;
        pm[i] -= eps;
        let fd = (dot(&dense_endpoint(pp, &u, 1), &terminal)
            - dot(&dense_endpoint(pm, &u, 1), &terminal))
            / (2.0 * eps);
        close(gradient.parameters[i], fd);
        assert!(fd.abs() > 1e-3, "parameter check must be informative");
    }
}

#[test]
fn reverse_trajectory_gradient_matches_independent_resolves_and_improves_design() {
    let model = Model { p: [4.0, 0.7] };
    let method = method();
    let initial = [0.9, -0.4];
    let target = [0.3, 0.1];
    let loss = |p| {
        let end = dense_endpoint(p, &initial, 12);
        0.5 * end
            .iter()
            .zip(target)
            .map(|(u, t)| (u - t).powi(2))
            .sum::<f64>()
    };
    let mut state = ImexState::new(0.0, &initial);
    let mut checkpoints = Vec::new();
    for _ in 0..12 {
        checkpoints.push(state.u.clone());
        method
            .step(&mut state, &model, &IdentityPreconditioner, &|x, y| {
                model.nonlinear(x, y)
            })
            .unwrap();
    }
    let mut bar = state
        .u
        .iter()
        .zip(target)
        .map(|(u, t)| u - t)
        .collect::<Vec<_>>();
    let mut parameters = [0.0; 2];
    for u in checkpoints.iter().rev() {
        let gradient = method
            .step_vjp(
                u,
                &model,
                &IdentityPreconditioner,
                &IdentityPreconditioner,
                &bar,
                10_000,
                &mut || false,
            )
            .unwrap();
        for (sum, value) in parameters.iter_mut().zip(&gradient.parameters) {
            *sum += value;
        }
        bar = gradient.initial;
    }
    for i in 0..2 {
        let (mut pp, mut pm) = (model.p, model.p);
        pp[i] += 1e-5;
        pm[i] -= 1e-5;
        close(parameters[i], (loss(pp) - loss(pm)) / 2e-5);
    }
    let improved = [model.p[0] - parameters[0], model.p[1] - parameters[1]];
    assert!(loss(improved) < loss(model.p));
}

#[test]
fn cancellation_during_second_stage_leaves_checkpoint_retryable() {
    let model = Model { p: [4.0, 0.7] };
    let method = method();
    let mut state = ImexState::new(0.0, &[0.9, -0.4]);
    // Preserve a real, nonempty accepted history as well as the current state.
    method
        .step(&mut state, &model, &IdentityPreconditioner, &|x, y| {
            model.nonlinear(x, y)
        })
        .unwrap();
    let before = state.clone();
    let calls = Cell::new(0);
    let error = method
        .step_controlled(
            &mut state,
            &model,
            &IdentityPreconditioner,
            &|x, y| {
                calls.set(calls.get() + 1);
                model.nonlinear(x, y);
            },
            &mut || calls.get() == 2,
        )
        .unwrap_err();
    assert!(matches!(error, ImexSolveError::Cancelled));
    assert_eq!(state, before);
    let mut straight = before;
    method
        .step(&mut straight, &model, &IdentityPreconditioner, &|x, y| {
            model.nonlinear(x, y)
        })
        .unwrap();
    method
        .step_controlled(
            &mut state,
            &model,
            &IdentityPreconditioner,
            &|x, y| model.nonlinear(x, y),
            &mut || false,
        )
        .unwrap();
    assert_eq!(state, straight);
}

#[test]
fn bounded_sweeps_refuse_without_publishing_partial_gradients() {
    let model = Model { p: [4.0, 0.7] };
    let method = method();
    let u = [0.9, -0.4];
    let seed = [0.3, -1.1];
    let required = method.adjoint_workspace_components(2).unwrap();
    assert!(matches!(
        method.step_vjp(
            &u,
            &model,
            &IdentityPreconditioner,
            &IdentityPreconditioner,
            &seed,
            required - 1,
            &mut || false
        ),
        Err(ImexAdjointError::WorkspaceLimit { .. })
    ));
    assert!(method.adjoint_workspace_components(usize::MAX).is_none());
    assert!(matches!(
        method.step_vjp(
            &u,
            &model,
            &IdentityPreconditioner,
            &IdentityPreconditioner,
            &[f64::NAN, 0.0],
            required,
            &mut || false
        ),
        Err(ImexAdjointError::InvalidInput(_))
    ));
    let polls = Cell::new(0);
    method
        .step_vjp(
            &u,
            &model,
            &IdentityPreconditioner,
            &IdentityPreconditioner,
            &seed,
            required,
            &mut || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    // Every polling boundary, including after the last derivative, is retryable.
    for stop in 1..=polls.get() {
        let mut calls = 0;
        assert!(matches!(
            method.step_vjp(
                &u,
                &model,
                &IdentityPreconditioner,
                &IdentityPreconditioner,
                &seed,
                required,
                &mut || {
                    calls += 1;
                    calls == stop
                }
            ),
            Err(ImexAdjointError::Step(ImexSolveError::Cancelled))
        ));
    }
    let strict = OperatorImex2::new(
        2,
        0.08,
        ImexSolveConfig {
            tolerance: 1e-15,
            restart: 1,
            max_cycles: 1,
        },
    );
    assert!(matches!(
        strict.step_vjp(
            &u,
            &model,
            &IdentityPreconditioner,
            &IdentityPreconditioner,
            &seed,
            10_000,
            &mut || false
        ),
        Err(ImexAdjointError::Step(ImexSolveError::NotConverged { .. }))
    ));
    struct BrokenTransposePreconditioner;
    impl fs_solver::FlexiblePreconditioner for BrokenTransposePreconditioner {
        fn apply(&self, _iteration: usize, _rhs: &[f64], out: &mut [f64]) {
            out.fill(f64::NAN);
        }
    }
    // Both forward stages succeed; refusal belongs to the reverse solve.
    assert!(matches!(
        method.step_vjp(
            &u,
            &model,
            &IdentityPreconditioner,
            &BrokenTransposePreconditioner,
            &seed,
            required,
            &mut || false
        ),
        Err(ImexAdjointError::Step(ImexSolveError::NotConverged {
            stage: fs_time::stiff::ImexStage::Two,
            ..
        }))
    ));
}

struct Broken(Model);
impl LinearOp for Broken {
    fn n(&self) -> usize {
        self.0.n()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply(x, y);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply_transpose(x, y);
    }
}
impl ImexVjp for Broken {
    fn parameter_count(&self) -> usize {
        2
    }
    fn nonlinear(&self, x: &[f64], y: &mut [f64]) {
        self.0.nonlinear(x, y);
    }
    fn nonlinear_vjp(
        &self,
        _x: &[f64],
        _seed: &[f64],
        _xb: &mut [f64],
        _pb: &mut [f64],
    ) -> Result<(), String> {
        Ok(()) // Deliberately leaves outputs unwritten.
    }
    fn linear_parameter_vjp(
        &self,
        _x: &[f64],
        _seed: &[f64],
        _pb: &mut [f64],
    ) -> Result<(), String> {
        Err("derivative unavailable".into())
    }
}
#[test]
fn missing_derivative_outputs_cannot_be_mistaken_for_zero_gradient() {
    assert!(matches!(
        method().step_vjp(
            &[0.9, -0.4],
            &Broken(Model { p: [4.0, 0.7] }),
            &IdentityPreconditioner,
            &IdentityPreconditioner,
            &[0.3, -1.1],
            10_000,
            &mut || false
        ),
        Err(ImexAdjointError::NonFiniteDerivative)
    ));
}
