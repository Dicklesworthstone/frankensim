//! Native reverse-mode product regressions against analytic derivatives and
//! the existing value interpreter. No derivative is estimated in production.
use fs_opt::reverse::{ReverseError, ReverseLimits, ReverseProgram};
use fs_opt::{Manifold, NodeId, OptError, ProblemBuilder, eval};
use fs_qty::Dims;

const LIMITS: ReverseLimits = ReverseLimits { max_nodes: 4096, max_scalar_slots: 65536 };
const UNIT: Dims = Dims([0; 6]);

fn scalar(builder: &mut ProblemBuilder) -> NodeId {
    let var = builder.var("x", Manifold::Rn { dim: 1 }, UNIT).unwrap();
    let vector = builder.var_ref(var).unwrap();
    builder.component(vector, 0).unwrap()
}
fn close(a: f64, b: f64) { assert!((a - b).abs() <= 2e-11 * (1.0 + b.abs()), "{a} != {b}"); }

#[test]
fn shared_graph_and_duplicate_roots_accumulate_every_path() {
    let mut b = ProblemBuilder::new();
    let x = scalar(&mut b);
    let square = b.mul(x, x).unwrap();
    assert_eq!(square, b.mul(x, x).unwrap());
    let problem = b.finish();
    let roots = [square, square, x];
    let program = ReverseProgram::new(&problem, &roots, LIMITS).unwrap();
    let inputs = vec![vec![3.0]];
    let evaluation = program.evaluate(&inputs, None).unwrap();
    assert_eq!(evaluation.values(), &[9.0, 9.0, 3.0]);
    assert_eq!(evaluation.pullback(&[3.0, -2.0, 4.0], None).unwrap(), vec![vec![10.0]]);
    assert_eq!(evaluation.pullback(&[0.0, 0.0, 0.0], None).unwrap(), vec![vec![0.0]]);
    assert_eq!(evaluation.pullback(&[3.0, -2.0, 4.0], None).unwrap(), vec![vec![10.0]]);
}

#[test]
fn every_smooth_scalar_rule_matches_analytic_derivatives_and_eval_bits() {
    let mut b = ProblemBuilder::new();
    let x = scalar(&mut b);
    let two = b.konst(2.0, UNIT).unwrap();
    let roots = [b.add(x, two).unwrap(), b.sub(x, two).unwrap(), b.mul(x, two).unwrap(),
        b.div(two, x).unwrap(), b.neg(x).unwrap(), b.powi(x, -3).unwrap(),
        b.sqrt(x).unwrap(), b.exp(x).unwrap(), b.ln(x).unwrap(), b.tanh(x).unwrap()];
    let problem = b.finish();
    let program = ReverseProgram::new(&problem, &roots, LIMITS).unwrap();
    let inputs = vec![vec![1.5]];
    let evaluation = program.evaluate(&inputs, None).unwrap();
    let x = 1.5_f64;
    let expected = [1.0, 1.0, 2.0, -2.0 / (x*x), -1.0, -3.0*x.powi(-4),
        0.5/x.sqrt(), x.exp(), 1.0/x, 1.0-x.tanh().powi(2)];
    for (i, root) in roots.iter().enumerate() {
        assert_eq!(evaluation.values()[i].to_bits(), eval(&problem, *root, &inputs).unwrap().scalar().unwrap().to_bits());
        let mut weights = vec![0.0; roots.len()]; weights[i] = 1.0;
        close(evaluation.pullback(&weights, None).unwrap()[0][0], expected[i]);
    }
}

#[test]
fn vector_broadcast_reductions_and_multiple_variable_blocks() {
    let mut b = ProblemBuilder::new();
    let v = b.var("v", Manifold::Rn { dim: 3 }, UNIT).unwrap();
    let vector = b.var_ref(v).unwrap();
    let s = scalar(&mut b);
    let left = b.mul(s, vector).unwrap();
    let right = b.mul(vector, s).unwrap();
    let sum = b.add(left, right).unwrap();
    let difference = b.sub(sum, vector).unwrap();
    let negative = b.neg(difference).unwrap();
    let root = b.dot(negative, vector).unwrap();
    let norm = b.norm_sq(vector).unwrap();
    let problem = b.finish();
    let program = ReverseProgram::new(&problem, &[root, norm], LIMITS).unwrap();
    let inputs = vec![vec![1.0, -2.0, 3.0], vec![2.0]];
    let evaluation = program.evaluate(&inputs, None).unwrap();
    assert_eq!(evaluation.values(), &[-42.0, 14.0]);
    assert_eq!(evaluation.pullback(&[1.0, 0.5], None).unwrap(), vec![vec![-5.0, 10.0, -15.0], vec![-28.0]]);
    for (i, root) in [root, norm].iter().enumerate() {
        assert_eq!(evaluation.values()[i].to_bits(), eval(&problem, *root, &inputs).unwrap().scalar().unwrap().to_bits());
    }
}

#[test]
fn one_lane_vector_is_not_confused_with_scalar_broadcast() {
    let mut b = ProblemBuilder::new();
    let v = b.var("v", Manifold::Rn { dim: 1 }, UNIT).unwrap();
    let vector = b.var_ref(v).unwrap();
    let s = b.component(vector, 0).unwrap();
    let scaled = b.mul(vector, s).unwrap();
    let root = b.norm_sq(scaled).unwrap();
    let problem = b.finish();
    let program = ReverseProgram::new(&problem, &[root], LIMITS).unwrap();
    assert_eq!(program.evaluate(&[vec![2.0]], None).unwrap().pullback(&[1.0], None).unwrap(), vec![vec![32.0]]);
}

#[test]
fn unused_kinks_physics_and_domain_failures_do_not_poison_a_root() {
    let mut b = ProblemBuilder::new();
    let x = scalar(&mut b);
    let unused = b.var("unused", Manifold::Rn { dim: 2 }, UNIT).unwrap();
    let _ = b.abs(x).unwrap();
    let _ = b.pde_residual("unavailable", unused, true, UNIT).unwrap();
    let minus = b.konst(-1.0, UNIT).unwrap();
    let _ = b.ln(minus).unwrap();
    let problem = b.finish();
    let program = ReverseProgram::new(&problem, &[x], LIMITS).unwrap();
    let evaluation = program.evaluate(&[vec![2.0], vec![8.0, 9.0]], None).unwrap();
    assert_eq!(evaluation.pullback(&[1.0], None).unwrap(), vec![vec![1.0], vec![0.0, 0.0]]);
    assert!(matches!(program.evaluate(&[vec![2.0]], None), Err(ReverseError::Evaluation(OptError::BindingCount { .. }))));
}

#[test]
fn reachable_kinks_and_external_nodes_are_refused_before_execution() {
    let mut b = ProblemBuilder::new();
    let var = b.var("x", Manifold::Rn { dim: 1 }, UNIT).unwrap();
    let v = b.var_ref(var).unwrap(); let x = b.component(v, 0).unwrap();
    let kink = b.abs(x).unwrap();
    let physics = b.pde_residual("physics", var, true, UNIT).unwrap();
    let uq = b.expectation(x, "samples").unwrap();
    let problem = b.finish();
    assert!(matches!(ReverseProgram::new(&problem, &[kink], LIMITS), Err(ReverseError::Nonsmooth { .. })));
    for node in [physics, uq] { assert!(matches!(ReverseProgram::new(&problem, &[node], LIMITS), Err(ReverseError::Evaluation(OptError::Unevaluable { .. })))); }
}

#[test]
fn malformed_bindings_roots_seeds_and_storage_caps_are_typed_refusals() {
    let mut b = ProblemBuilder::new();
    let var = b.var("x", Manifold::Rn { dim: 1 }, UNIT).unwrap();
    let v = b.var_ref(var).unwrap(); let x = b.component(v, 0).unwrap();
    let problem = b.finish();
    assert!(matches!(ReverseProgram::new(&problem, &[v], LIMITS), Err(ReverseError::Evaluation(OptError::NotScalar { .. }))));
    assert!(matches!(ReverseProgram::new(&problem, &[NodeId(999)], LIMITS), Err(ReverseError::Evaluation(OptError::UnknownNode { .. }))));
    assert!(ReverseProgram::new(&problem, &[x], ReverseLimits { max_nodes: 1, ..LIMITS }).is_err());
    assert!(ReverseProgram::new(&problem, &[x], ReverseLimits { max_scalar_slots: 1, ..LIMITS }).is_err());
    let program = ReverseProgram::new(&problem, &[x], LIMITS).unwrap();
    assert!(matches!(program.evaluate(&[vec![1.0, 2.0]], None), Err(ReverseError::Evaluation(OptError::BindingLen { .. }))));
    assert!(matches!(program.evaluate(&[vec![f64::NAN]], None), Err(ReverseError::Evaluation(OptError::BindingNonFinite { .. }))));
    let evaluation = program.evaluate(&[vec![2.0]], None).unwrap();
    assert!(matches!(evaluation.pullback(&[], None), Err(ReverseError::SeedCount { .. })));
    assert!(matches!(evaluation.pullback(&[f64::INFINITY], None), Err(ReverseError::NonFiniteAdjoint { .. })));
    assert_eq!(evaluation.pullback(&[1.0], None).unwrap(), vec![vec![1.0]]);
}

#[test]
fn singular_derivative_refuses_without_corrupting_cached_primals() {
    let mut b = ProblemBuilder::new(); let x = scalar(&mut b); let root = b.sqrt(x).unwrap();
    let problem = b.finish(); let program = ReverseProgram::new(&problem, &[root], LIMITS).unwrap();
    assert!(matches!(program.evaluate(&[vec![-1.0]], None), Err(ReverseError::Evaluation(OptError::EvalNonFinite { .. }))));
    let evaluation = program.evaluate(&[vec![0.0]], None).unwrap();
    assert!(matches!(evaluation.pullback(&[1.0], None), Err(ReverseError::NonFiniteAdjoint { .. })));
    assert_eq!(evaluation.values(), &[0.0]);
    assert_eq!(evaluation.pullback(&[0.0], None).unwrap(), vec![vec![0.0]]);
}

#[test]
fn integer_power_edges_and_small_derivatives_are_not_erased() {
    let mut b = ProblemBuilder::new(); let x = scalar(&mut b);
    let roots = [b.powi(x, 0).unwrap(), b.powi(x, 1).unwrap(), b.powi(x, 2).unwrap()];
    let problem = b.finish(); let program = ReverseProgram::new(&problem, &roots, LIMITS).unwrap();
    let evaluation = program.evaluate(&[vec![1e-200]], None).unwrap();
    assert_eq!(evaluation.pullback(&[1.0, 0.0, 0.0], None).unwrap()[0][0], 0.0);
    assert_eq!(evaluation.pullback(&[0.0, 1.0, 0.0], None).unwrap()[0][0], 1.0);
    assert_eq!(evaluation.pullback(&[0.0, 0.0, 1.0], None).unwrap()[0][0], 2e-200);
    let mut b = ProblemBuilder::new(); let x = scalar(&mut b); let root = b.powi(x, i32::MIN).unwrap();
    let problem = b.finish(); let program = ReverseProgram::new(&problem, &[root], LIMITS).unwrap();
    assert_eq!(program.evaluate(&[vec![1.0]], None).unwrap().pullback(&[1.0], None).unwrap()[0][0], i32::MIN as f64);
}

#[test]
fn saturated_tanh_retains_a_representable_tail_derivative() {
    let mut b = ProblemBuilder::new(); let x = scalar(&mut b); let root = b.tanh(x).unwrap();
    let problem = b.finish(); let program = ReverseProgram::new(&problem, &[root], LIMITS).unwrap();
    let gradient = program.evaluate(&[vec![20.0]], None).unwrap().pullback(&[1.0], None).unwrap()[0][0];
    assert!(gradient > 0.0);
    assert!((gradient / (4.0 * (-40.0_f64).exp()) - 1.0).abs() < 1e-12);
}

#[test]
fn manifold_parameter_gradients_match_retraction_directional_derivatives() {
    for (manifold, point) in [
        (Manifold::Sphere { ambient: 3 }, vec![0.6, 0.8, 0.0]),
        (Manifold::So3, vec![1.0, 0.0, 0.0, 0.0]),
        (Manifold::Stiefel { n: 3, p: 2 }, vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
    ] {
        let mut b = ProblemBuilder::new(); let var = b.var("point", manifold, UNIT).unwrap();
        let v = b.var_ref(var).unwrap(); let root = b.component(v, 1).unwrap();
        let problem = b.finish(); let program = ReverseProgram::new(&problem, &[root], LIMITS).unwrap();
        let evaluation = program.evaluate(&[point.clone()], None).unwrap();
        let gradient = evaluation.parameter_pullback(&[1.0], None).unwrap();
        assert_eq!(gradient[0].len(), manifold.param_dim().unwrap() as usize);
        for j in 0..gradient[0].len() {
            let mut step = vec![0.0; gradient[0].len()]; step[j] = 1e-6;
            let plus = manifold.retract(&point, &step).unwrap(); step[j] = -1e-6;
            let minus = manifold.retract(&point, &step).unwrap();
            assert!((gradient[0][j] - (plus[1] - minus[1]) / 2e-6).abs() < 2e-7);
        }
    }
}

#[test]
fn large_vector_uses_one_graph_sweep_not_one_evaluation_per_coordinate() {
    let mut b = ProblemBuilder::new(); let var = b.var("x", Manifold::Rn { dim: 4096 }, UNIT).unwrap();
    let v = b.var_ref(var).unwrap(); let root = b.norm_sq(v).unwrap();
    let problem = b.finish(); let program = ReverseProgram::new(&problem, &[root], LIMITS).unwrap();
    assert_eq!(program.node_count(), 2); assert_eq!(program.scalar_slots(), 4097);
    let input = vec![vec![0.25; 4096]];
    let evaluation = program.evaluate(&input, None).unwrap();
    assert_eq!(evaluation.values(), &[256.0]);
    assert_eq!(evaluation.pullback(&[1.0], None).unwrap(), vec![vec![0.5; 4096]]);
}
