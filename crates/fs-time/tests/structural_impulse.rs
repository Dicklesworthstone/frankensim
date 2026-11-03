//! G1/G3/G4: prescribed impulses retain full inertia and restart damped motion.
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    impulse::{ImpulseControls, ImpulseError, impulse_workspace_components, second_order_impulse},
    initialization::InitialSolveConfig,
};

struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, input: &[f64], output: &mut [f64]) {
        output.copy_from_slice(input);
    }
}

/// Independent full consistent mass with C = 2 M and zero stiffness.
struct Damped {
    refuse_positive_velocity: bool,
    bad_force: bool,
}
impl SecondOrderProblem for Damped {
    fn dimension(&self) -> usize {
        2
    }
    fn validate_state(&self, _: &[f64], v: &[f64]) -> Result<(), String> {
        if self.refuse_positive_velocity && v[0] > 2.0 {
            Err("velocity outside physical domain".into())
        } else {
            Ok(())
        }
    }
    fn mass_apply(&self, x: &[f64], output: &mut [f64]) {
        output[0] = 2.0 * x[0] + x[1];
        output[1] = x[0] + 3.0 * x[1];
    }
    fn damping_apply(&self, x: &[f64], output: &mut [f64]) {
        self.mass_apply(x, output);
        for value in output {
            *value *= 2.0;
        }
    }
    fn internal_force(&self, _: &[f64], output: &mut [f64]) {
        output.fill(if self.bad_force { f64::NAN } else { 0.0 });
    }
    fn tangent_apply(&self, _: &[f64], _: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
}

fn model() -> Damped {
    Damped {
        refuse_positive_velocity: false,
        bad_force: false,
    }
}
fn controls() -> ImpulseControls {
    ImpulseControls {
        mass_solve: InitialSolveConfig {
            restart: 2,
            max_cycles: 4,
            tolerance: 1e-13,
        },
        absolute_impulse_tolerance: 1e-12,
        relative_impulse_tolerance: 1e-12,
        work_tolerance: 1e-11,
    }
}
fn capacity() -> usize {
    impulse_workspace_components(2, controls().mass_solve).unwrap()
}
fn state() -> SecondOrderState {
    let mut state = SecondOrderState::new(0.75, &[0.5, -0.25], &[1.0, -2.0], &[-2.0, 4.0]);
    state.steps = 9;
    state
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn full_mass_jump_and_post_event_force_determine_the_new_acceleration() {
    let mut state = state();
    let before = state.clone();
    let report = second_order_impulse(
        &mut state,
        &model(),
        &[5.0, 5.0],
        &[1.0, -2.0],
        &Identity,
        controls(),
        capacity(),
        &mut || false,
    )
    .unwrap();
    // Independent inverse: M^-1 = [3,-1;-1,2]/5.
    for (actual, expected) in state.v.iter().zip([3.0, -1.0]) {
        close(*actual, expected);
    }
    for (actual, expected) in state.a.iter().zip([-5.0, 1.0]) {
        close(*actual, expected);
    }
    for (actual, expected) in report.velocity_change.iter().zip([2.0, 1.0]) {
        close(*actual, expected);
    }
    for value in report.momentum_change {
        close(value, 5.0);
    }
    close(report.kinetic_energy_change, 2.5);
    close(report.supplied_work, 2.5);
    close(report.work_residual, 0.0);
    assert!(report.velocity_solve.converged && report.acceleration_solve.converged);
    assert_eq!(state.q, before.q);
    assert_eq!(state.t, before.t);
    assert_eq!(state.steps, before.steps);
    assert_eq!(state.history, before.history);
}

#[test]
fn restarted_damped_trajectory_converges_to_closed_form_without_smearing_the_strike() {
    let run = |h: f64| {
        let mut state = SecondOrderState::new(0.0, &[0.0; 2], &[0.0; 2], &[0.0; 2]);
        second_order_impulse(
            &mut state,
            &model(),
            &[5.0, 5.0],
            &[0.0; 2],
            &Identity,
            controls(),
            capacity(),
            &mut || false,
        )
        .unwrap();
        let stepper = OperatorGeneralizedAlpha::new(2, h, 1.0, ImplicitSolveConfig::default());
        for _ in 0..(0.5 / h).round() as usize {
            stepper.step(&mut state, &model(), &[0.0; 2]).unwrap();
        }
        // v(t) = [2,1] exp(-2t), q(t) = [1,1/2] (1-exp(-2t)).
        let expected = (-1.0_f64).exp();
        (state.v[0] - 2.0 * expected).abs() + (state.q[1] - 0.5 * (1.0 - expected)).abs()
    };
    let coarse = run(0.05);
    let fine = run(0.025);
    assert!(
        coarse > 3.8 * fine && coarse < 4.2 * fine,
        "{coarse:e} {fine:e}"
    );
    assert!(fine < 3e-4);
}

#[test]
fn simultaneous_impulses_include_cross_work_and_braking_extracts_energy() {
    let initial = SecondOrderState::new(0.0, &[0.0; 2], &[0.0; 2], &[0.0; 2]);
    let mut combined = initial.clone();
    let combined_report = second_order_impulse(
        &mut combined,
        &model(),
        &[5.0, 5.0],
        &[0.0; 2],
        &Identity,
        controls(),
        capacity(),
        &mut || false,
    )
    .unwrap();
    let mut split = initial;
    let first = second_order_impulse(
        &mut split,
        &model(),
        &[5.0, 0.0],
        &[0.0; 2],
        &Identity,
        controls(),
        capacity(),
        &mut || false,
    )
    .unwrap();
    let second = second_order_impulse(
        &mut split,
        &model(),
        &[0.0, 5.0],
        &[0.0; 2],
        &Identity,
        controls(),
        capacity(),
        &mut || false,
    )
    .unwrap();
    close(combined_report.supplied_work, 7.5);
    close(first.supplied_work + second.supplied_work, 7.5);
    for (a, b) in split.v.iter().zip(&combined.v) {
        close(*a, *b);
    }
    let brake = second_order_impulse(
        &mut combined,
        &model(),
        &[-5.0, -5.0],
        &[0.0; 2],
        &Identity,
        controls(),
        capacity(),
        &mut || false,
    )
    .unwrap();
    close(brake.supplied_work, -7.5);
    close(brake.kinetic_energy_change, -7.5);
    for value in combined.v {
        close(value, 0.0);
    }
}

#[test]
fn zero_impulse_restarts_finite_force_without_changing_momentum() {
    let mut state = state();
    let velocity = state.v.clone();
    let report = second_order_impulse(
        &mut state,
        &model(),
        &[0.0; 2],
        &[1.0, -2.0],
        &Identity,
        controls(),
        capacity(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(state.v, velocity);
    close(state.a[0], -1.0);
    close(state.a[1], 3.0);
    assert_eq!(report.supplied_work, 0.0);
}

#[test]
fn rounded_away_impulse_refuses_even_when_the_mass_solve_converges() {
    let mut state = SecondOrderState::new(0.0, &[0.0; 2], &[1e20; 2], &[0.0; 2]);
    let before = state.clone();
    assert!(matches!(
        second_order_impulse(
            &mut state,
            &model(),
            &[5.0, 5.0],
            &[0.0; 2],
            &Identity,
            controls(),
            capacity(),
            &mut || false,
        ),
        Err(ImpulseError::Momentum { .. })
    ));
    assert_eq!(state, before);
}

#[test]
fn failed_acceleration_admission_budget_and_cancellation_never_publish_half_an_event() {
    let mut original = state();
    // Retain real accepted-step history through both successful and refused kicks.
    OperatorGeneralizedAlpha::new(2, 0.01, 0.8, ImplicitSolveConfig::default())
        .step(&mut original, &model(), &[0.0; 2])
        .unwrap();
    let mut polls = 0;
    let mut accepted = original.clone();
    second_order_impulse(
        &mut accepted,
        &model(),
        &[5.0, 5.0],
        &[0.0; 2],
        &Identity,
        controls(),
        capacity(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(accepted.history, original.history);
    for cutoff in 1..=polls {
        let mut state = original.clone();
        let mut count = 0;
        let result = second_order_impulse(
            &mut state,
            &model(),
            &[5.0, 5.0],
            &[0.0; 2],
            &Identity,
            controls(),
            capacity(),
            &mut || {
                count += 1;
                count >= cutoff
            },
        );
        assert!(
            matches!(result, Err(ImpulseError::Cancelled)),
            "poll {cutoff}: {result:?}"
        );
        assert_eq!(state, original);
    }
    for bad in [
        Damped {
            refuse_positive_velocity: true,
            bad_force: false,
        },
        Damped {
            refuse_positive_velocity: false,
            bad_force: true,
        },
    ] {
        let mut state = original.clone();
        assert!(
            second_order_impulse(
                &mut state,
                &bad,
                &[5.0, 5.0],
                &[0.0; 2],
                &Identity,
                controls(),
                capacity(),
                &mut || false,
            )
            .is_err()
        );
        assert_eq!(state, original);
    }
    let mut state = original.clone();
    assert!(
        second_order_impulse(
            &mut state,
            &model(),
            &[5.0, 5.0],
            &[0.0; 2],
            &Identity,
            controls(),
            capacity() - 1,
            &mut || false,
        )
        .is_err()
    );
    assert_eq!(state, original);
}
