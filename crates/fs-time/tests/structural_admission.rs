//! A material envelope is enforced before numerical callbacks or publication.
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    TimeSolveError,
    initialization::{
        InitialSolveConfig, InitialSolveError, second_order_acceleration,
        second_order_acceleration_vjp,
    },
    second_order_adjoint::{SecondOrderAdjointConfig, SecondOrderAdjointError, SecondOrderVjp},
};
use std::cell::Cell;

// Free unit mass with an admitted displacement, speed, and step envelope.
struct BoundedMass {
    max_step: f64,
    last_step: Cell<f64>,
    calls: Cell<usize>,
    last_force_q: Cell<f64>,
}
impl BoundedMass {
    fn new(max_step: f64) -> Self {
        Self {
            max_step,
            last_step: Cell::new(f64::NAN),
            calls: Cell::new(0),
            last_force_q: Cell::new(f64::NAN),
        }
    }
    fn touch(&self) {
        self.calls.set(self.calls.get() + 1);
    }
    fn checked_zero(&self, q: &[f64], v: &[f64], out: &mut [f64]) {
        assert!(
            self.validate_state(q, v).is_ok(),
            "force or tangent evaluated outside the admitted material domain"
        );
        self.touch();
        out.fill(0.0);
    }
}
impl SecondOrderProblem for BoundedMass {
    fn dimension(&self) -> usize {
        1
    }
    fn validate_step(&self, step_s: f64) -> Result<(), String> {
        self.last_step.set(step_s);
        if step_s > self.max_step {
            Err("step exceeds source envelope".into())
        } else {
            Ok(())
        }
    }
    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        if !q[0].is_finite() || q[0] > 0.15 || !v[0].is_finite() || v[0].abs() > 2.0 {
            Err("state exceeds source envelope".into())
        } else {
            Ok(())
        }
    }
    fn mass_apply(&self, x: &[f64], out: &mut [f64]) {
        self.touch();
        out.copy_from_slice(x);
    }
    fn damping_apply(&self, _: &[f64], out: &mut [f64]) {
        self.touch();
        out.fill(0.0);
    }
    fn internal_force(&self, q: &[f64], out: &mut [f64]) {
        self.checked_zero(q, &[0.0], out);
    }
    fn tangent_apply(&self, q: &[f64], _: &[f64], out: &mut [f64]) {
        self.checked_zero(q, &[0.0], out);
    }
    fn resistance(&self, q: &[f64], v: &[f64], out: &mut [f64]) {
        self.checked_zero(q, v, out);
        self.last_force_q.set(q[0]);
    }
    fn resistance_displacement_apply(&self, q: &[f64], v: &[f64], _: &[f64], out: &mut [f64]) {
        self.checked_zero(q, v, out);
    }
    fn resistance_velocity_apply(&self, q: &[f64], v: &[f64], _: &[f64], out: &mut [f64]) {
        self.checked_zero(q, v, out);
    }
}
impl SecondOrderVjp for BoundedMass {
    fn parameter_count(&self) -> usize {
        0
    }
    fn mass_transpose_apply(&self, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.mass_apply(seed, out);
        Ok(())
    }
    fn damping_transpose_apply(&self, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.damping_apply(seed, out);
        Ok(())
    }
    fn tangent_transpose_apply(&self, q: &[f64], _: &[f64], out: &mut [f64]) -> Result<(), String> {
        self.internal_force(q, out);
        Ok(())
    }
    fn residual_parameter_vjp(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[f64],
        _: &[f64],
        out: &mut [f64],
    ) -> Result<(), String> {
        out.fill(0.0);
        Ok(())
    }
}
struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, input: &[f64], output: &mut [f64]) {
        output.copy_from_slice(input);
    }
}
fn method(h: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(1, h, 1.0, ImplicitSolveConfig::default())
}
fn assert_adjoint_refused(
    method: OperatorGeneralizedAlpha,
    model: &BoundedMass,
    state: &SecondOrderState,
    forcing: f64,
) {
    let error = method
        .step_vjp(
            state,
            model,
            &[forcing],
            (&[1.0], &[0.0], &[0.0]),
            &Identity,
            SecondOrderAdjointConfig::default(),
            10_000,
            &mut || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        SecondOrderAdjointError::Step(TimeSolveError::ModelRefused { .. })
    ));
}
#[test]
fn actual_step_is_admitted_before_forward_or_adjoint_numerics() {
    let model = BoundedMass::new(0.15);
    let mut state = SecondOrderState::new(0.0, &[0.0], &[0.0], &[0.0]);
    let before = state.clone();
    // rho=1 evaluates at h/2; admitting that instead of h would allow this step.
    assert!(matches!(
        method(0.2).step(&mut state, &model, &[0.0]),
        Err(TimeSolveError::ModelRefused { .. })
    ));
    assert_eq!(model.last_step.get(), 0.2);
    assert_eq!(state, before);
    assert_adjoint_refused(method(0.2), &model, &state, 0.0);
    assert_eq!(model.calls.get(), 0);
}
#[test]
fn invalid_initial_displacement_or_velocity_refuses_acceleration_and_pullback() {
    let model = BoundedMass::new(0.3);
    for (q, v) in [(0.2, 0.0), (0.0, 3.0)] {
        let config = InitialSolveConfig::default();
        assert!(matches!(
            second_order_acceleration(
                &model,
                &[q],
                &[v],
                &[0.0],
                &Identity,
                config,
                10_000,
                &mut || false
            ),
            Err(InitialSolveError::ModelRefused(_))
        ));
        assert!(matches!(
            second_order_acceleration_vjp(
                &model,
                &[q],
                &[v],
                &[0.0],
                &[1.0],
                &Identity,
                &Identity,
                config,
                10_000,
                &mut || false
            ),
            Err(InitialSolveError::ModelRefused(_))
        ));
        let mut state = SecondOrderState::new(0.0, &[q], &[v], &[0.0]);
        let before = state.clone();
        assert!(matches!(
            method(0.1).step(&mut state, &model, &[0.0]),
            Err(TimeSolveError::ModelRefused { .. })
        ));
        assert_eq!(state, before);
        assert_adjoint_refused(method(0.1), &model, &state, 0.0);
    }
    assert_eq!(model.calls.get(), 0);
}
#[test]
fn inadmissible_newton_trial_never_reaches_force_or_tangent() {
    let model = BoundedMass::new(0.3);
    let mut state = SecondOrderState::new(0.0, &[0.0], &[0.0], &[0.0]);
    let before = state.clone();
    // Admissible predictor; F=20 requires q_new=.4 and refused q_stage=.2.
    assert!(matches!(
        method(0.2).step(&mut state, &model, &[20.0]),
        Err(TimeSolveError::ModelRefused { .. })
    ));
    assert!(model.calls.get() > 0);
    assert_eq!(state, before);
    assert_adjoint_refused(method(0.2), &model, &state, 20.0);
}
#[test]
fn valid_stage_cannot_publish_an_invalid_endpoint_or_partial_history() {
    let model = BoundedMass::new(0.3);
    let mut state = SecondOrderState::new(-0.05, &[-0.05], &[1.0], &[0.0]);
    method(0.05).step(&mut state, &model, &[0.0]).unwrap();
    assert_eq!(state.steps, 1);
    assert_eq!(state.history.len(), 1);
    let before = state.clone();
    // Constant velocity gives admitted q_stage=.1 but refused q_new=.2.
    assert!(matches!(
        method(0.2).step(&mut state, &model, &[0.0]),
        Err(TimeSolveError::ModelRefused { .. })
    ));
    assert!((model.last_force_q.get() - 0.1).abs() < 1e-14);
    assert_eq!(state, before);
    assert_adjoint_refused(method(0.2), &model, &state, 0.0);
    assert_eq!(state, before);
}
