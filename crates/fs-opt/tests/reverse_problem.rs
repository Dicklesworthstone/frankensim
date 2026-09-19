//! Solver-facing objective, constraint and manifold pullbacks.
use fs_opt::reverse::ReverseLimits;
use fs_opt::{
    ConstraintKind, Manifold, Problem, ProblemBuilder, ProblemTag,
    ReverseProblem, ReverseProblemError, Sense,
};
use fs_qty::Dims;

fn limits() -> ReverseLimits { ReverseLimits { max_nodes: 1 << 18, max_scalar_slots: 1 << 22 } }

fn fixture() -> Problem {
    let mut b = ProblemBuilder::new();
    let xv = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let sv = b.var("s", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let x = b.var_ref(xv).unwrap(); let sr = b.var_ref(sv).unwrap();
    let a = b.component(x, 0).unwrap(); let c = b.component(x, 1).unwrap();
    let s = b.component(sr, 0).unwrap(); let norm = b.norm_sq(x).unwrap();
    let sum = b.add(a, c).unwrap(); let equality = b.sub(sum, s).unwrap();
    let one = b.konst(1.0, Dims::NONE).unwrap(); let inequality = b.sub(a, one).unwrap();
    b.objective(norm, Sense::Minimize, 2.0).unwrap();
    b.objective(s, Sense::Maximize, 0.5).unwrap();
    b.constraint(equality, ConstraintKind::EqZero, "balance").unwrap();
    b.constraint(inequality, ConstraintKind::LeZero, "upper-bound").unwrap();
    b.finish()
}

#[test]
fn weighted_objective_and_full_constraint_pullbacks_use_declared_semantics() {
    let p = fixture();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut input = vec![vec![2.0, 3.0], vec![4.0]];
    let tape = oracle.evaluate(&input).unwrap();
    assert_eq!(tape.objective_values(), &[13.0, 4.0]);
    assert_eq!(tape.objective_value(), 24.0);
    assert_eq!(tape.constraint_values(), &[1.0, 1.0]);
    assert_eq!(tape.objective_gradient().unwrap(), vec![vec![8.0, 12.0], vec![-0.5]]);
    assert_eq!(tape.constraint_pullback(&[5.0, 7.0]).unwrap(), vec![vec![12.0, 5.0], vec![-5.0]]);
    assert_eq!(tape.lagrangian_gradient(&[5.0, 7.0]).unwrap(), vec![vec![20.0, 17.0], vec![-5.5]]);
    assert_eq!(oracle.problem().constraints()[0].kind, ConstraintKind::EqZero);
    assert_eq!(oracle.problem().constraints()[1].kind, ConstraintKind::LeZero);
    // Retained valuation is independent of subsequently modified input storage.
    input[0][0] = 200.0;
    input[1][0] = -100.0;
    assert_eq!(tape.objective_value(), 24.0);
    assert_eq!(tape.objective_parameter_gradient().unwrap(), vec![vec![8.0, 12.0], vec![-0.5]]);
    assert_eq!(tape.lagrangian_parameter_gradient(&[5.0, 7.0]).unwrap(), vec![vec![20.0, 17.0], vec![-5.5]]);
    assert_eq!(tape.constraint_parameter_pullback(&[5.0, 7.0]).unwrap(), tape.constraint_pullback(&[5.0, 7.0]).unwrap());
}

#[test]
fn feasible_inequalities_are_not_clipped_or_removed() {
    let p = fixture();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let input = vec![vec![-2.0, 3.0], vec![4.0]];
    let tape = oracle.evaluate(&input).unwrap();
    assert_eq!(tape.constraint_values(), &[-3.0, -3.0]);
    assert_eq!(tape.constraint_pullback(&[0.0, 1.0]).unwrap(), vec![vec![1.0, 0.0], vec![0.0]]);
}

#[test]
fn constraint_seed_errors_retain_caller_constraint_indices() {
    let p = fixture();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let input = vec![vec![2.0, 3.0], vec![4.0]];
    let tape = oracle.evaluate(&input).unwrap();
    assert!(matches!(tape.constraint_pullback(&[1.0]), Err(ReverseProblemError::SeedCount { expected: 2, actual: 1 })));
    assert!(matches!(tape.lagrangian_gradient(&[0.0, f64::INFINITY]), Err(ReverseProblemError::SeedNonFinite { index: 1, .. })));
    assert_eq!(tape.objective_value(), 24.0);
    assert_eq!(tape.objective_gradient().unwrap(), vec![vec![8.0, 12.0], vec![-0.5]]);
}

#[test]
fn mixed_manifold_gradients_use_authoritative_parameter_coordinates() {
    let mut b = ProblemBuilder::new();
    let manifolds = [Manifold::Sphere { ambient: 3 }, Manifold::So3, Manifold::Stiefel { n: 3, p: 2 }];
    for (i, manifold) in manifolds.iter().enumerate() {
        let v = b.var(&format!("v{i}"), *manifold, Dims::NONE).unwrap();
        let r = b.var_ref(v).unwrap(); let component = b.component(r, if i == 2 { 2 } else { 1 }).unwrap();
        b.objective(component, Sense::Minimize, 1.0).unwrap();
    }
    let p = b.finish();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let input = vec![vec![1.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0]];
    let tape = oracle.evaluate(&input).unwrap();
    let ambient = tape.objective_gradient().unwrap();
    let parameter = tape.objective_parameter_gradient().unwrap();
    assert_eq!(ambient[1].len(), 4);
    assert_eq!(parameter[1].len(), 3);
    assert_eq!(parameter[1], vec![0.5, 0.0, 0.0]);
    assert_eq!(tape.lagrangian_parameter_gradient(&[]).unwrap(), parameter);
    for (i, manifold) in manifolds.iter().enumerate() {
        assert_eq!(parameter[i], manifold.parameter_gradient(&input[i], &ambient[i]).unwrap());
        for j in 0..parameter[i].len() {
            let h = 1e-6;
            let mut step = vec![0.0; parameter[i].len()]; step[j] = h;
            let mut positive = input.clone(); positive[i] = manifold.retract(&input[i], &step).unwrap();
            step[j] = -h;
            let mut negative = input.clone(); negative[i] = manifold.retract(&input[i], &step).unwrap();
            let difference = (oracle.evaluate(&positive).unwrap().objective_value() - oracle.evaluate(&negative).unwrap().objective_value()) / (2.0 * h);
            assert!((difference - parameter[i][j]).abs() < 1e-8, "block {i}, lane {j}: {difference} != {}", parameter[i][j]);
        }
    }
}

#[test]
fn problem_adapter_refuses_unexecuted_structural_tags_and_missing_objectives() {
    let mut b = ProblemBuilder::new();
    let root = b.konst(1.0, Dims::NONE).unwrap();
    let no_objective = b.finish();
    assert!(matches!(ReverseProblem::new(&no_objective, limits()), Err(ReverseProblemError::NoObjectives)));
    let mut b = ProblemBuilder::new();
    let f = b.konst(1.0, Dims::NONE).unwrap();
    assert_eq!(root, f);
    b.objective(f, Sense::Minimize, 1.0).unwrap();
    b.tag(ProblemTag::MultiFidelity { levels: 2 }).unwrap();
    let tagged = b.finish();
    assert!(matches!(ReverseProblem::new(&tagged, limits()), Err(ReverseProblemError::UnsupportedProblemTags)));
}

#[test]
fn scalarization_overflow_is_not_a_finite_objective_report() {
    let mut b = ProblemBuilder::new();
    let f = b.konst(2.0, Dims::NONE).unwrap();
    b.objective(f, Sense::Minimize, f64::MAX).unwrap();
    let p = b.finish();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    assert!(matches!(oracle.evaluate(&[]), Err(ReverseProblemError::NonFiniteObjective { index: 0, .. })));
}

#[test]
fn duplicate_objective_and_constraint_roots_accumulate_with_their_own_weights() {
    let mut b = ProblemBuilder::new();
    let variable = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let reference = b.var_ref(variable).unwrap();
    let root = b.norm_sq(reference).unwrap();
    b.objective(root, Sense::Minimize, 3.0).unwrap();
    b.objective(root, Sense::Maximize, 1.0).unwrap();
    b.constraint(root, ConstraintKind::LeZero, "shared-root").unwrap();
    let p = b.finish();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let tape = oracle.evaluate(&[vec![2.0]]).unwrap();
    assert_eq!(tape.objective_values(), &[4.0, 4.0]);
    assert_eq!(tape.objective_value(), 8.0);
    assert_eq!(tape.constraint_values(), &[4.0]);
    assert_eq!(tape.objective_gradient().unwrap(), vec![vec![8.0]]);
    assert_eq!(tape.constraint_pullback(&[5.0]).unwrap(), vec![vec![20.0]]);
    assert_eq!(tape.lagrangian_gradient(&[5.0]).unwrap(), vec![vec![28.0]]);
    // Pullbacks are reusable, not destructive accumulations in the primal tape.
    assert_eq!(tape.lagrangian_gradient(&[0.0]).unwrap(), vec![vec![8.0]]);
}
