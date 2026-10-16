//! Native IR, exact reverse HVPs, retained tapes and real Cx integration.
use super::*;
use fs_ascent::trust::PreconditionerError;
use std::cell::Cell;

fn stiff_diagonal() -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("design", Manifold::Rn { dim: 16 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let mut root = b.konst(0.0, Dims::NONE).unwrap();
    for i in 0..16 {
        let x = b.component(r, i).unwrap();
        let square = b.powi(x, 2).unwrap();
        let weight = b.konst(0.5 * (2.0f64).powi(i as i32), Dims::NONE).unwrap();
        let term = b.mul(weight, square).unwrap();
        root = b.add(root, term).unwrap();
    }
    b.objective(root, Sense::Minimize, 1.0).unwrap();
    b.finish()
}

fn coupled() -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("xy", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap(); let y = b.component(r, 1).unwrap();
    let xx = b.powi(x, 2).unwrap(); let yy = b.powi(y, 2).unwrap();
    let xy = b.mul(x, y).unwrap();
    let a = b.konst(2.0, Dims::NONE).unwrap(); let c = b.konst(1.5, Dims::NONE).unwrap();
    let axx = b.mul(a, xx).unwrap(); let cyy = b.mul(c, yy).unwrap();
    let sum = b.add(axx, xy).unwrap(); let f = b.add(sum, cyy).unwrap();
    b.objective(f, Sense::Minimize, 1.0).unwrap();
    b.finish()
}
fn diagonal(_: &[f64], r: &[f64]) -> Vec<f64> { vec![r[0] / 4.0, r[1] / 3.0] }
fn retained(a: &ReverseNewtonStudy<'_, '_>, b: &ReverseNewtonStudy<'_, '_>) {
    let (a0, b0) = (a.snapshot(), b.snapshot());
    let bits = |x: &[f64]| x.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&a0.x), bits(&b0.x)); assert_eq!(a0.f.to_bits(), b0.f.to_bits());
    assert_eq!(bits(a.gradient()), bits(b.gradient())); assert_eq!(bits(a.history()), bits(b.history()));
    assert_eq!(a.radius().to_bits(), b.radius().to_bits()); assert_eq!(a0.iters, b0.iters);
    assert_eq!(a0.negative_curvature_hits, b0.negative_curvature_hits);
}

#[test]
fn exact_ir_hessian_uses_one_preconditioned_direction_for_stiff_design() {
    let p = stiff_diagonal(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseNewtonStudy::new(&oracle, &[0.125; 16], None).unwrap();
    let mut plain = study.clone(); let rule = StopRule::GradNorm(1e-10);
    let report = study.run_preconditioned(
        &mut |_, r| r.iter().enumerate().map(|(i, r)| r / (2.0f64).powi(i as i32)).collect(),
        &rule, 1, 2, 1, None,
    ).unwrap();
    assert_eq!(report.stop, ReverseNewtonStop::Stopped(StopReason::GradNorm));
    assert_eq!(report.solution.x, vec![0.0; 16]); assert_eq!(report.solution.f, 0.0);
    assert_eq!(report.solution.hv_evals, 2); assert_eq!(report.solution.evals, 2);
    assert_eq!(report.preconditioner_evals, 1); assert_eq!(study.preconditioner_evals(), 1);
    assert!(plain.run(&rule, 1, 100, None).unwrap().solution.hv_evals > 10);
    assert_eq!(plain.preconditioner_evals(), 0);
}

#[test]
fn discarded_ir_searches_charge_both_ceilings_before_exact_budget_restart() {
    let p = coupled(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseNewtonStudy::new(&oracle, &[0.2, -0.3], None).unwrap();
    let start = study.clone(); let rule = StopRule::GradNorm(1e-12);
    for _ in 0..2 {
        let out = study.run_preconditioned(&mut diagonal, &rule, 10, 100, 1, None).unwrap();
        assert_eq!(out.stop, ReverseNewtonStop::PreconditionerBudget);
        assert_eq!(out.preconditioner_evals, 1); assert_eq!(out.solution.hv_evals, 1);
        assert_eq!(out.solution.evals, 1); retained(&study, &start);
    }
    // The staged kernel's count is still zero, but the study has already paid
    // for one discarded direction. Only ONE further HVP may be admitted here.
    let out = study.run_preconditioned(&mut diagonal, &rule, 10, 2, 100, None).unwrap();
    assert_eq!(out.stop, ReverseNewtonStop::HessianBudget);
    assert_eq!(out.solution.hv_evals, 2); assert_eq!(out.preconditioner_evals, 2);
    retained(&study, &start);
    let out = study.run_preconditioned(&mut diagonal, &rule, 10, 5, 4, None).unwrap();
    assert_eq!(out.stop, ReverseNewtonStop::Stopped(StopReason::GradNorm));
    assert!(out.solution.x.iter().all(|x| x.abs() < 1e-15));
    assert_eq!(out.solution.evals, 2); assert_eq!(out.solution.hv_evals, 5);
    assert_eq!(out.preconditioner_evals, 4); assert_eq!(out.solution.iters, 1);
}

#[test]
fn inverse_faults_after_a_real_ir_product_retain_tape_and_actual_attempts() {
    for fault in 0..3 {
        let p = coupled(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
        let mut study = ReverseNewtonStudy::new(&oracle, &[0.2, -0.3], None).unwrap();
        let start = study.clone(); let rule = StopRule::GradNorm(1e-12);
        let calls = Cell::new(0usize);
        let out = study.run_preconditioned(&mut |x, r| {
            calls.set(calls.get() + 1);
            if calls.get() == 1 { return diagonal(x, r); }
            match fault {
                0 => vec![],
                1 => vec![f64::NAN; 2],
                _ => panic!("injected inverse failure after an exact IR product"),
            }
        }, &rule, 10, 100, 100, None);
        match (fault, out) {
            (0, Err(ReverseNewtonError::Preconditioner(PreconditionerError::Dimension { expected: 2, actual: 0 })))
            | (1, Err(ReverseNewtonError::Preconditioner(PreconditionerError::NonFinite { component: 0 })))
            | (2, Err(ReverseNewtonError::Preconditioner(PreconditionerError::CallbackPanicked))) => {},
            (_, other) => panic!("wrong failure attribution: {other:?}"),
        }
        retained(&study, &start); assert_eq!(study.snapshot().hv_evals, 1);
        assert_eq!(study.snapshot().evals, 1); assert_eq!(study.preconditioner_evals(), 2);
        let out = study.run_preconditioned(&mut diagonal, &rule, 10, 4, 4, None).unwrap();
        assert_eq!(out.stop, ReverseNewtonStop::Stopped(StopReason::GradNorm));
        assert_eq!(out.solution.hv_evals, 4); assert_eq!(out.preconditioner_evals, 4);
        assert_eq!(out.solution.evals, 2);
    }
}

#[test]
fn real_context_cancellation_inside_inverse_prevents_the_next_ir_product() {
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    for at in 1..=2 {
        let p = coupled(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
        let mut study = ReverseNewtonStudy::new(&oracle, &[0.2, -0.3], None).unwrap();
        let start = study.clone(); let calls = Cell::new(0usize);
        let gate = CancelGate::new(); let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(&gate, arena,
                StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic);
            let out = study.run_preconditioned(&mut |x, r| {
                calls.set(calls.get() + 1);
                if calls.get() == at { gate.request(); }
                diagonal(x, r)
            }, &StopRule::GradNorm(1e-12), 10, 100, 100, Some(&cx));
            assert!(matches!(out, Err(ReverseNewtonError::Cancelled)), "{out:?}");
        });
        retained(&study, &start); assert_eq!(study.preconditioner_evals(), at);
        assert_eq!(study.snapshot().hv_evals, at - 1); assert_eq!(study.snapshot().evals, 1);
        let out = study.run_preconditioned(&mut diagonal, &StopRule::GradNorm(1e-12),
            10, at + 2, at + 2, None).unwrap();
        assert_eq!(out.stop, ReverseNewtonStop::Stopped(StopReason::GradNorm));
        assert_eq!(out.preconditioner_evals, at + 2); assert_eq!(out.solution.evals, 2);
    }
}

#[test]
fn preconditioned_ir_splits_replay_and_identity_keeps_legacy_trajectory() {
    let p = rosenbrock(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let start = ReverseNewtonStudy::new(&oracle, &[-1.2, 1.0], None).unwrap();
    let rule = StopRule::GradNorm(1e-12);
    let mut legacy = start.clone(); legacy.run(&rule, 12, 1000, None).unwrap();
    for identity in [false, true] {
        let mut inverse = |x: &[f64], r: &[f64]| {
            if identity { r.to_vec() }
            else { vec![r[0] / (1200.0 * x[0] * x[0] + 400.0 * x[1].abs() + 2.0), r[1] / 200.0] }
        };
        let mut straight = start.clone();
        straight.run_preconditioned(&mut inverse, &rule, 12, 1000, 1000, None).unwrap();
        if identity { same(&straight, &legacy); }
        for cut in 0..=12 {
            let mut split = start.clone();
            split.run_preconditioned(&mut inverse, &rule, cut, 1000, 1000, None).unwrap();
            split.run_preconditioned(&mut inverse, &rule, 12 - cut, 1000, 1000, None).unwrap();
            same(&straight, &split);
            assert_eq!(straight.preconditioner_evals(), split.preconditioner_evals());
        }
    }
}

#[test]
fn objective_and_iteration_admission_do_not_spend_inverse_work() {
    let p = quadratic(1, Some(1)); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseNewtonStudy::new(&oracle, &[1.0], None).unwrap();
    let mut forbidden = |_: &[f64], _: &[f64]| panic!("unadmitted inverse callback");
    let out = study.run_preconditioned(&mut forbidden, &StopRule::GradNorm(1e-12),
        10, 0, 0, None).unwrap();
    assert_eq!(out.stop, ReverseNewtonStop::Stopped(StopReason::Budget));
    assert_eq!(out.preconditioner_evals, 0); assert_eq!(out.solution.hv_evals, 0);
    let p = quadratic(1, None); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseNewtonStudy::new(&oracle, &[1.0], None).unwrap();
    let out = study.run_preconditioned(&mut forbidden, &StopRule::GradNorm(1e-12),
        0, 100, 100, None).unwrap();
    assert_eq!(out.stop, ReverseNewtonStop::Stopped(StopReason::IterationCap));
    assert_eq!(out.preconditioner_evals, 0); assert_eq!(out.solution.hv_evals, 0);
    assert_eq!(study.run_preconditioned(&mut forbidden, &StopRule::GradNorm(1e-12),
        10, 100, 0, None).unwrap().stop, ReverseNewtonStop::PreconditionerBudget);
}
