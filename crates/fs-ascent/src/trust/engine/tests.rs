use super::*;

fn rosenbrock(x: &[f64]) -> (f64, Vec<f64>) {
    let residual = x[1] - x[0] * x[0];
    (
        100.0 * residual * residual + (1.0 - x[0]).powi(2),
        vec![-400.0 * x[0] * residual - 2.0 * (1.0 - x[0]), 200.0 * residual],
    )
}

fn rosenbrock_hv(x: &[f64], v: &[f64]) -> Vec<f64> {
    vec![
        (1200.0 * x[0] * x[0] - 400.0 * x[1] + 2.0) * v[0]
            - 400.0 * x[0] * v[1],
        -400.0 * x[0] * v[0] + 200.0 * v[1],
    ]
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

fn assert_same_state(a: &TrustRegionState, b: &TrustRegionState) {
    assert_eq!(bits(&a.x), bits(&b.x));
    assert_eq!(a.f.to_bits(), b.f.to_bits());
    assert_eq!(bits(&a.g), bits(&b.g));
    assert_eq!(a.delta.to_bits(), b.delta.to_bits());
    assert_eq!(bits(&a.history), bits(&b.history));
    assert_eq!(a.iters, b.iters);
    assert_eq!(a.evals, b.evals);
    assert_eq!(a.hv_evals, b.hv_evals);
    assert_eq!(a.negative_curvature_hits, b.negative_curvature_hits);
    assert_eq!(a.stalled, b.stalled);
}

#[test]
fn trust_region_every_split_preserves_complete_state() {
    let initial = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
    let rule = StopRule::GradNorm(1e-12);
    let mut straight = initial.clone();
    straight.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 12);
    for split in 0..=12 {
        let mut segmented = initial.clone();
        segmented.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, split);
        let mut resumed = segmented.clone();
        resumed.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 12 - split);
        assert_same_state(&straight, &resumed);
    }
}

#[test]
fn trust_region_pause_is_mutation_free_and_resumable() {
    let initial = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
    let rule = StopRule::GradNorm(1e-12);
    let mut paused = initial.clone();
    paused.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 2);
    let outcome = paused.run_with_pause(
        &mut rosenbrock, &mut rosenbrock_hv, &rule, 10,
        &mut || true,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Paused);
    assert_eq!(outcome.solution.iters, 2);
    let checkpoint = paused.clone();
    let repeated = paused.run_with_pause(
        &mut |_| panic!("objective evaluated during pause"),
        &mut |_, _| panic!("Hessian evaluated during pause"),
        &rule, 10, &mut || true,
    );
    assert_eq!(repeated.progress, TrustRegionProgress::Paused);
    assert_same_state(&checkpoint, &paused);
    paused.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 10);
    let mut straight = initial;
    straight.run(&mut rosenbrock, &mut rosenbrock_hv, &rule, 12);
    assert_same_state(&straight, &paused);
}

// G4: every callback-side checkpoint of a real quadratic proposal must
// retain its work charges but publish none of the incomplete iteration.
#[test]
fn trust_region_pause_at_each_callback_boundary_is_atomic() {
    let mut fg = |x: &[f64]| ((x[0] - 2.0).powi(2), vec![2.0 * (x[0] - 2.0)]);
    let mut hv = |_: &[f64], v: &[f64]| vec![2.0 * v[0]];
    let initial = TrustRegionState::new(&[0.0], &mut fg);
    let rule = StopRule::GradNorm(1e-12);
    let mut straight = initial.clone();
    straight.run(&mut fg, &mut hv, &rule, 10);
    // Entry, before/after CG Hv, before/after final-model Hv,
    // before/after the candidate objective+gradient.
    for (boundary, spent_hv, spent_fg) in [
        (1, 0, 0), (2, 0, 0), (3, 1, 0), (4, 1, 0),
        (5, 2, 0), (6, 2, 0), (7, 2, 1),
    ] {
        let mut state = initial.clone();
        let mut polls = 0;
        let outcome = state.run_with_pause(&mut fg, &mut hv, &rule, 10, &mut || {
            polls += 1;
            polls >= boundary
        });
        assert_eq!(outcome.progress, TrustRegionProgress::Paused);
        let mut expected = initial.clone();
        expected.hv_evals += spent_hv;
        expected.evals += spent_fg;
        assert_same_state(&expected, &state);
        let checkpoint = state.clone();
        state.run_with_pause(
            &mut |_| panic!("paused objective"), &mut |_, _| panic!("paused Hessian"),
            &rule, 10, &mut || true,
        );
        assert_same_state(&checkpoint, &state);
        state.run(&mut fg, &mut hv, &rule, 10);
        let mut expected = straight.clone();
        expected.hv_evals += spent_hv;
        expected.evals += spent_fg;
        assert_same_state(&expected, &state);
    }
}

#[test]
fn trust_region_pause_during_later_krylov_product_retains_only_work() {
    use std::cell::Cell;
    let mut fg = |x: &[f64]| (0.5 * x[0] * x[0] + x[1] * x[1], vec![x[0], 2.0 * x[1]]);
    let mut state = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    let before = state.clone();
    let calls = Cell::new(0);
    let outcome = state.run_with_pause(
        &mut |_| panic!("cancelled Krylov solve reached objective"),
        &mut |_, v| { calls.set(calls.get() + 1); vec![v[0], 2.0 * v[1]] },
        &StopRule::GradNorm(1e-12), 10, &mut || calls.get() == 2,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Paused);
    let mut expected = before;
    expected.hv_evals += 2;
    assert_same_state(&expected, &state);
    state.run(&mut fg, &mut |_, v| vec![v[0], 2.0 * v[1]], &StopRule::GradNorm(1e-12), 10);
    assert!(state.report().grad_norm <= 1e-12);
}

#[test]
fn trust_region_cancelled_negative_curvature_trial_is_not_committed() {
    use std::cell::Cell;
    let mut fg = |x: &[f64]| (-x[0] * x[0], vec![-2.0 * x[0]]);
    let mut state = TrustRegionState::new(&[1.0], &mut fg);
    let mut expected = state.clone();
    let calls = Cell::new(0);
    let outcome = state.run_with_pause(
        &mut |_| panic!("cancelled model reached objective"),
        &mut |_, v| { calls.set(calls.get() + 1); vec![-2.0 * v[0]] },
        &StopRule::GradNorm(0.0), 1, &mut || calls.get() == 2,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Paused);
    expected.hv_evals += 2;
    assert_same_state(&expected, &state);
    state.run(&mut fg, &mut |_, v| vec![-2.0 * v[0]], &StopRule::GradNorm(0.0), 1);
    assert_eq!(state.report().negative_curvature_hits, 1);
    assert_eq!(state.report().x, vec![2.0]);
    assert_eq!(state.radius(), 2.0);
}

#[test]
fn trust_region_nested_all_budget_is_a_hard_callback_ceiling() {
    use std::cell::Cell;
    let calls = Cell::new(0);
    let mut fg = |x: &[f64]| { calls.set(calls.get() + 1); rosenbrock(x) };
    let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut fg);
    let rule = StopRule::All(vec![
        StopRule::GradNorm(0.0),
        StopRule::Any(vec![StopRule::Budget(3), StopRule::Budget(5)]),
    ]);
    let outcome = state.run(&mut fg, &mut rosenbrock_hv, &rule, 100);
    assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Budget));
    assert_eq!(outcome.solution.evals, 3);
    assert_eq!(calls.get(), 3);
    let checkpoint = state.clone();
    state.run(&mut fg, &mut rosenbrock_hv, &rule, 100);
    assert_eq!(calls.get(), 3);
    assert_same_state(&checkpoint, &state);
}

#[test]
fn trust_region_budget_has_priority_over_convergence() {
    let mut fg = |x: &[f64]| (x[0] * x[0], vec![2.0 * x[0]]);
    let mut state = TrustRegionState::new(&[0.0], &mut fg);
    let outcome = state.run(
        &mut fg, &mut |_, _| panic!("unexpected Hessian evaluation"),
        &StopRule::Any(vec![StopRule::GradNorm(1e-8), StopRule::Budget(1)]), 100,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Budget));
    assert_eq!(outcome.solution.iters, 0);
    assert_eq!(outcome.solution.evals, 1);
}

#[test]
fn trust_region_zero_iteration_segment_spends_no_callbacks() {
    let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
    let checkpoint = state.clone();
    let outcome = state.run(
        &mut |_| panic!("unexpected objective evaluation"),
        &mut |_, _| panic!("unexpected Hessian evaluation"),
        &StopRule::GradNorm(1e-12), 0,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::IterationCap));
    assert_same_state(&checkpoint, &state);
}

#[test]
fn trust_region_wrapper_and_state_return_the_same_solution() {
    let legacy = trust_region_newton(
        &[-1.2, 1.0], &mut rosenbrock, &mut rosenbrock_hv, 1e-8, 100,
    );
    let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
    let run = state.run(&mut rosenbrock, &mut rosenbrock_hv, &StopRule::GradNorm(1e-8), 100);
    assert_eq!(run.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
    assert!(run.solution.grad_norm <= 1e-8);
    assert_eq!(bits(&legacy.x), bits(&run.solution.x));
    assert_eq!(legacy.f.to_bits(), run.solution.f.to_bits());
    assert_eq!(legacy.grad_norm.to_bits(), run.solution.grad_norm.to_bits());
    assert_eq!(legacy.iters, run.solution.iters);
    assert_eq!(legacy.evals, run.solution.evals);
    assert_eq!(legacy.hv_evals, run.solution.hv_evals);
    assert_eq!(legacy.negative_curvature_hits, run.solution.negative_curvature_hits);
}

#[test]
#[should_panic(expected = "budget is below retained evaluations")]
fn trust_region_refuses_underwritten_resume_budget() {
    let mut state = TrustRegionState::new(&[-1.2, 1.0], &mut rosenbrock);
    state.run(
        &mut |_| panic!("unexpected objective evaluation"),
        &mut |_, _| panic!("unexpected Hessian evaluation"),
        &StopRule::Budget(0), 1,
    );
}

#[test]
fn trust_region_rejects_nonfinite_objectives_without_poisoning_state() {
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut fg = |x: &[f64]| {
            let f = if x[0] == 1.0 { invalid } else { (x[0] - 2.0).powi(2) };
            (f, vec![2.0 * (x[0] - 2.0)])
        };
        let mut hv = |_: &[f64], v: &[f64]| vec![2.0 * v[0]];
        let mut state = TrustRegionState::new(&[0.0], &mut fg);
        let run = state.run(&mut fg, &mut hv, &StopRule::GradNorm(1e-10), 1);
        assert_eq!(run.solution.x, vec![0.0]);
        assert_eq!(run.solution.f, 4.0);
        assert_eq!(state.gradient(), &[-4.0]);
        assert_eq!(state.radius(), 0.25);
        assert_eq!(state.history(), &[4.0, 4.0]);
        assert_eq!(run.solution.evals, 2);
        let recovered = state.run(&mut fg, &mut hv, &StopRule::GradNorm(1e-10), 100);
        assert_eq!(recovered.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
        assert!((recovered.solution.x[0] - 2.0).abs() < 1e-10);
    }
}

#[test]
fn trust_region_finite_objective_with_nan_gradient_is_not_accepted() {
    let mut fg = |x: &[f64]| {
        let g = if x[0] == 1.0 { f64::NAN } else { 2.0 * (x[0] - 2.0) };
        ((x[0] - 2.0).powi(2), vec![g])
    };
    let mut state = TrustRegionState::new(&[0.0], &mut fg);
    let run = state.run(
        &mut fg, &mut |_, v| vec![2.0 * v[0]], &StopRule::GradNorm(1e-10), 1,
    );
    assert_eq!(run.solution.x, vec![0.0]);
    assert_eq!(run.solution.grad_norm, 4.0);
    assert_eq!(state.radius(), 0.25);
}

#[test]
fn trust_region_nonfinite_hessian_stalls_without_objective_trials() {
    let mut fg = |x: &[f64]| (x[0] * x[0], vec![2.0 * x[0]]);
    let mut state = TrustRegionState::new(&[1.0], &mut fg);
    let outcome = state.run(
        &mut |_| panic!("invalid Hessian must not reach the objective"),
        &mut |_, _| vec![f64::NAN], &StopRule::GradNorm(1e-8), 100,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Stall));
    assert_eq!(outcome.solution.x, vec![1.0]);
    assert_eq!(outcome.solution.grad_norm, 2.0);
    assert_eq!(outcome.solution.evals, 1);
    assert_eq!(outcome.solution.hv_evals, 1);
    let checkpoint = state.clone();
    state.run(
        &mut |_| panic!("stalled state evaluated the objective"),
        &mut |_, _| panic!("stalled state evaluated the Hessian"),
        &StopRule::GradNorm(1e-8), 100,
    );
    assert_same_state(&checkpoint, &state);
}

#[test]
fn trust_region_nonpositive_model_cannot_accept_an_uphill_step() {
    let mut fg = |x: &[f64]| (x[0] * x[0], vec![2.0 * x[0]]);
    let mut state = TrustRegionState::new(&[1.0], &mut fg);
    let mut calls = 0;
    let mut inconsistent_hv = |_: &[f64], v: &[f64]| {
        calls += 1;
        let scale = if calls == 1 { 2.0 } else { 10.0 };
        vec![scale * v[0]]
    };
    let outcome = state.run(
        &mut |_| panic!("nonpositive model must not reach the objective"),
        &mut inconsistent_hv, &StopRule::GradNorm(1e-8), 100,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Stall));
    assert_eq!(outcome.solution.x, vec![1.0]);
    assert_eq!(outcome.solution.evals, 1);
    assert_eq!(outcome.solution.hv_evals, 2);
}

#[test]
fn trust_region_precision_noop_stops_instead_of_reevaluating() {
    let mut state = TrustRegionState::new(&[1e100], &mut |_| (0.0, vec![-1.0]));
    let outcome = state.run(
        &mut |_| panic!("unchanged point must not be reevaluated"),
        &mut |_, _| vec![0.0], &StopRule::GradNorm(1e-8), 100,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Stopped(StopReason::Stall));
    assert_eq!(outcome.solution.evals, 1);
    assert_eq!(outcome.solution.x, vec![1e100]);
}

#[test]
fn trust_region_large_finite_gradient_takes_a_boundary_step() {
    let mut fg = |x: &[f64]| (1e308 * (1.0 - x[0]), vec![-1e308]);
    let mut state = TrustRegionState::new(&[0.0], &mut fg);
    let outcome = state.run(&mut fg, &mut |_, _| vec![0.0], &StopRule::GradNorm(0.0), 1);
    assert_eq!(
        outcome.progress,
        TrustRegionProgress::Stopped(StopReason::IterationCap)
    );
    assert!((outcome.solution.x[0] - 1.0).abs() < 1e-15);
    assert!(outcome.solution.f.abs() / 1e308 < 1e-15);
    assert_eq!(outcome.solution.evals, 2);
    assert_eq!(outcome.solution.hv_evals, 2);
    assert_eq!(state.radius(), 2.0);
}

// G1/G3: multiplying f, g and H by a positive constant preserves the
// manufactured minimizer. Previously both ends stalled before one trial.
#[test]
fn trust_region_convex_quadratics_survive_extreme_objective_units() {
    for scale in [1e-250, 1e-200, 1.0, 1e200, 1e250] {
        let mut fg = |x: &[f64]| {
            let a = x[0] - 2.0;
            let b = x[1] + 1.0;
            (
                scale * (a * a + 4.0 * b * b),
                vec![scale * 2.0 * a, scale * 8.0 * b],
            )
        };
        let mut hv = |_: &[f64], v: &[f64]| vec![scale * (2.0 * v[0]), scale * (8.0 * v[1])];
        let mut state = TrustRegionState::new(&[0.0, 0.0], &mut fg);
        let rule = StopRule::GradNorm(scale * 1e-10);
        let outcome = state.run(&mut fg, &mut hv, &rule, 30);
        assert_eq!(
            outcome.progress,
            TrustRegionProgress::Stopped(StopReason::GradNorm),
            "objective scale {scale}"
        );
        assert!((outcome.solution.x[0] - 2.0).abs() < 1e-10);
        assert!((outcome.solution.x[1] + 1.0).abs() < 1e-10);
        assert!(outcome.solution.f / scale < 1e-20);
        assert!(outcome.solution.evals > 1);
    }
}

// G0/G3: a concave direction must still reach the correct boundary when
// the original model decrease is subnormal or exceeds the old norm range.
#[test]
fn trust_region_negative_curvature_survives_subnormal_objective_units() {
    for scale in [
        f64::from_bits(1u64 << 24),
        1e-310,
        1e-200,
        1.0,
        1e200,
        1e300,
    ] {
        let mut fg = |x: &[f64]| {
            (
                scale * (-0.5 * x[0] * x[0] + 0.5 * x[1] * x[1]),
                vec![-scale * x[0], scale * x[1]],
            )
        };
        let mut state = TrustRegionState::new(&[1.0, 0.0], &mut fg);
        let outcome = state.run(
            &mut fg,
            &mut |_, v| vec![-scale * v[0], scale * v[1]],
            &StopRule::GradNorm(0.0),
            1,
        );
        assert_eq!(
            outcome.solution.x,
            vec![2.0, 0.0],
            "objective scale {scale}"
        );
        assert_eq!(outcome.solution.negative_curvature_hits, 1);
        assert_eq!(outcome.solution.evals, 2);
        assert_eq!(outcome.solution.hv_evals, 2);
        assert_eq!(state.radius(), 2.0);
    }
}

#[test]
fn trust_region_scaled_agreement_survives_actual_decrease_overflow() {
    let mut fg = |x: &[f64]| (9e307 - 1.4e308 * x[0] - 1.4e308 * x[1], vec![-1.4e308; 2]);
    let mut state = TrustRegionState::new(&[0.0, 0.0], &mut fg);
    let outcome = state.run(
        &mut fg,
        &mut |_, _| vec![0.0; 2],
        &StopRule::GradNorm(0.0),
        1,
    );
    assert!((9e307 - outcome.solution.f).is_infinite());
    assert!(outcome.solution.f.is_finite() && outcome.solution.f < 0.0);
    for value in outcome.solution.x {
        assert!((value - 1.0 / 2.0f64.sqrt()).abs() < 1e-15);
    }
    assert_eq!(state.radius(), 2.0);
}

#[test]
fn trust_region_tiny_positive_curvature_reaches_finite_boundary() {
    let mut fg = |x: &[f64]| (-x[0] + (0.5e-320 * x[0]) * x[0], vec![-1.0 + 1e-320 * x[0]]);
    let mut state = TrustRegionState::new(&[0.0], &mut fg);
    let outcome = state.run(
        &mut fg,
        &mut |_, v| vec![1e-320 * v[0]],
        &StopRule::GradNorm(0.0),
        1,
    );
    assert_eq!(outcome.solution.x, vec![1.0]);
    assert_eq!(outcome.solution.negative_curvature_hits, 0);
    assert_eq!(state.radius(), 2.0);
}

#[test]
fn trust_region_representable_tiny_model_decrease_is_not_zero() {
    let mut fg = |x: &[f64]| {
        (
            1e-10 * x[0] + (0.5e290 * x[0]) * x[0],
            vec![1e-10 + 1e290 * x[0]],
        )
    };
    let mut state = TrustRegionState::new(&[0.0], &mut fg);
    let outcome = state.run(
        &mut fg,
        &mut |_, v| vec![1e290 * v[0]],
        &StopRule::GradNorm(1e-25),
        3,
    );
    assert_eq!(
        outcome.progress,
        TrustRegionProgress::Stopped(StopReason::GradNorm)
    );
    assert!((outcome.solution.x[0] / -1e-300 - 1.0).abs() < 1e-14);
    assert!(outcome.solution.f < 0.0 && outcome.solution.f.is_finite());
}

// G5/G4: rescaling is local model arithmetic, so checkpoints still retain
// everything necessary and cancelled normalized products remain charged.
#[test]
fn trust_region_scaled_models_replay_after_pause_and_every_split() {
    use std::cell::Cell;
    let scale = 1e200;
    let mut fg = |x: &[f64]| {
        let (f, mut g) = rosenbrock(x);
        for value in &mut g {
            *value *= scale;
        }
        (scale * f, g)
    };
    let mut hv = |x: &[f64], v: &[f64]| {
        rosenbrock_hv(x, v)
            .into_iter()
            .map(|value| scale * value)
            .collect()
    };
    let initial = TrustRegionState::new(&[-1.2, 1.0], &mut fg);
    let rule = StopRule::GradNorm(scale * 1e-12);
    let mut straight = initial.clone();
    straight.run(&mut fg, &mut hv, &rule, 12);
    for split in 0..=12 {
        let mut segmented = initial.clone();
        segmented.run(&mut fg, &mut hv, &rule, split);
        segmented.run(&mut fg, &mut hv, &rule, 12 - split);
        assert_same_state(&straight, &segmented);
    }
    let calls = Cell::new(0usize);
    let mut paused = initial;
    let outcome = paused.run_with_pause(
        &mut |_| panic!("cancelled scaled model reached objective"),
        &mut |x, v| {
            calls.set(calls.get() + 1);
            hv(x, v)
        },
        &rule,
        12,
        &mut || calls.get() == 1,
    );
    assert_eq!(outcome.progress, TrustRegionProgress::Paused);
    paused.run(&mut fg, &mut hv, &rule, 12);
    straight.hv_evals += 1;
    assert_same_state(&straight, &paused);
}

#[test]
#[should_panic(expected = "Hessian dimension mismatch")]
fn trust_region_refuses_malformed_hessian_products() {
    let mut state = TrustRegionState::new(&[1.0], &mut |_| (1.0, vec![2.0]));
    state.run(
        &mut |_| panic!("bad Hessian reached objective"),
        &mut |_, _| vec![], &StopRule::GradNorm(1e-8), 1,
    );
}

#[test]
#[should_panic(expected = "gradient dimension mismatch")]
fn trust_region_refuses_malformed_trial_gradients() {
    let mut state = TrustRegionState::new(&[1.0], &mut |_| (1.0, vec![2.0]));
    state.run(
        &mut |_| (0.0, vec![]), &mut |_, v| vec![2.0 * v[0]],
        &StopRule::GradNorm(1e-8), 1,
    );
}

#[test]
#[should_panic(expected = "initial evaluation")]
fn trust_region_refuses_nan_initial_gradients() {
    let _ = TrustRegionState::new(&[1.0], &mut |_| (1.0, vec![f64::NAN]));
}
