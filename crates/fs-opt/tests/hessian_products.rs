//! Mathematical HVPs on the actual live IR and shared primal evaluator.
use fs_opt::reverse::{HessianError, ReverseError, ReverseLimits, ReverseProgram};
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, ReverseProblemError, Sense};
use fs_qty::Dims;

fn limits() -> ReverseLimits {
    ReverseLimits { max_nodes: 100_000, max_scalar_slots: 1_000_000 }
}
fn close(actual: &[Vec<f64>], expected: &[Vec<f64>], tol: f64) {
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected) {
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(b) {
            assert!((a - b).abs() <= tol * (1.0 + b.abs()), "{a} != {b}");
        }
    }
}

#[test]
fn quartic_and_negative_curvature_have_analytic_products() {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap(); let y = b.component(r, 1).unwrap();
    let x4 = b.powi(x, 4).unwrap(); let yy = b.powi(y, 2).unwrap();
    let xy = b.mul(x, y).unwrap(); let z = b.sub(x4, yy).unwrap();
    let z = b.add(z, xy).unwrap(); b.objective(z, Sense::Minimize, 1.0).unwrap();
    let p = b.finish(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let tape = oracle.evaluate(&[vec![0.0, 1.0]]).unwrap();
    close(&tape.objective_hessian_vector_product(&[vec![2.0, 3.0]], None).unwrap(),
        &[vec![3.0, -4.0]], 0.0);
    let tape = oracle.evaluate(&[vec![2.0, 1.0]]).unwrap();
    close(&tape.objective_hessian_vector_product(&[vec![2.0, 3.0]], None).unwrap(),
        &[vec![99.0, -4.0]], 0.0);
}

#[test]
fn vector_scaling_shared_reductions_and_duplicate_roots_accumulate() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 3 }, Dims::NONE).unwrap();
    let s = b.var("s", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let x = b.var_ref(x).unwrap(); let sr = b.var_ref(s).unwrap();
    let s = b.component(sr, 0).unwrap(); let sx = b.mul(s, x).unwrap();
    let neg = b.neg(sx).unwrap(); let difference = b.sub(sx, neg).unwrap();
    let dot = b.dot(difference, x).unwrap(); // 2*s*|x|^2
    let norm = b.norm_sq(sx).unwrap(); // s^2*|x|^2
    let root = b.add(dot, norm).unwrap();
    let p = b.finish(); let program = ReverseProgram::new(&p, &[root, root], limits()).unwrap();
    let tape = program.evaluate(&[vec![1.0, 2.0, 3.0], vec![2.0]], None).unwrap();
    let d = vec![vec![3.0, -1.0, 2.0], vec![0.5]];
    // For (s^2+2s)*q: Hxx=16I, Hxs=12*x, Hss=2q=28.
    let expected = vec![vec![54.0, -4.0, 50.0], vec![98.0]];
    close(&tape.hessian_vector_product(&[2.0, -1.0], &d, None).unwrap(), &expected, 0.0);
    assert_eq!(tape.hessian_vector_product(&[1.0, -1.0], &d, None).unwrap(), vec![vec![0.0; 3], vec![0.0]]);
    close(&tape.hessian_vector_product(&[1.0, 0.0], &d, None).unwrap(), &expected, 0.0);
}

#[test]
fn all_unary_rules_and_quotient_match_gradient_differences_and_symmetry() {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap(); let x = b.component(r, 0).unwrap(); let y = b.component(r, 1).unwrap();
    let roots = [b.sqrt(x).unwrap(), b.exp(x).unwrap(), b.ln(x).unwrap(), b.tanh(x).unwrap(),
        b.powi(x, -3).unwrap(), b.div(x, y).unwrap(), b.powi(y, 3).unwrap()];
    for root in roots { b.objective(root, Sense::Minimize, 1.0).unwrap(); }
    let p = b.finish(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let x = vec![vec![1.3, 0.8]]; let v = vec![vec![0.3, -0.7]]; let u = vec![vec![0.6, 0.2]];
    let tape = oracle.evaluate(&x).unwrap();
    let hv = tape.objective_hessian_vector_product(&v, None).unwrap();
    let hu = tape.objective_hessian_vector_product(&u, None).unwrap();
    let h = 1e-5;
    let plus = vec![(0..2).map(|i| x[0][i] + h * v[0][i]).collect()];
    let minus = vec![(0..2).map(|i| x[0][i] - h * v[0][i]).collect()];
    let gp = oracle.evaluate(&plus).unwrap().objective_gradient().unwrap();
    let gm = oracle.evaluate(&minus).unwrap().objective_gradient().unwrap();
    let fd = vec![(0..2).map(|i| (gp[0][i] - gm[0][i]) / (2.0 * h)).collect()];
    close(&hv, &fd, 2e-8);
    let symmetry: f64 = (0..2).map(|i| u[0][i] * hv[0][i] - v[0][i] * hu[0][i]).sum();
    assert!(symmetry.abs() < 1e-12);
    assert_eq!(tape.objective_hessian_vector_product(&v, None).unwrap(), hv);
}

#[test]
fn objective_senses_and_nonlinear_constraint_multipliers_are_preserved() {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap(); let x = b.component(r, 0).unwrap(); let y = b.component(r, 1).unwrap();
    let xx = b.powi(x, 2).unwrap(); let yy = b.powi(y, 2).unwrap(); let xy = b.mul(x, y).unwrap();
    b.objective(xx, Sense::Minimize, 2.0).unwrap(); b.objective(yy, Sense::Maximize, 0.5).unwrap();
    b.constraint(xy, ConstraintKind::EqZero, "product").unwrap();
    b.constraint(yy, ConstraintKind::LeZero, "square").unwrap();
    let p = b.finish(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let tape = oracle.evaluate(&[vec![0.0, 0.0]]).unwrap(); let d = [vec![2.0, 3.0]];
    close(&tape.objective_hessian_vector_product(&d, None).unwrap(), &[vec![8.0, -3.0]], 0.0);
    close(&tape.lagrangian_hessian_vector_product(&[5.0, 7.0], &d, None).unwrap(), &[vec![23.0, 49.0]], 0.0);
    assert!(matches!(tape.lagrangian_hessian_vector_product(&[1.0], &d, None), Err(ReverseProblemError::SeedCount { .. })));
}

#[test]
fn integer_power_boundary_exponents_do_not_overflow_exponent_arithmetic() {
    for power in [0, 1, 2, 3, -1, i32::MIN, i32::MIN + 1] {
        let mut b = ProblemBuilder::new();
        let v = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
        let r = b.var_ref(v).unwrap(); let x = b.component(r, 0).unwrap();
        let z = b.powi(x, power).unwrap(); let p = b.finish();
        let program = ReverseProgram::new(&p, &[z], limits()).unwrap();
        let tape = program.evaluate(&[vec![1.0]], None).unwrap();
        let e = f64::from(power);
        close(&tape.hessian_vector_product(&[1.0], &[vec![0.25]], None).unwrap(),
            &[vec![0.25 * e * (e - 1.0)]], 1e-15);
    }
}

#[test]
fn zero_weight_singular_root_is_not_differentiated_but_active_singularity_refuses() {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap(); let x = b.component(r, 0).unwrap();
    let square = b.powi(x, 2).unwrap(); let sqrt = b.sqrt(x).unwrap(); let p = b.finish();
    let program = ReverseProgram::new(&p, &[square, sqrt], limits()).unwrap();
    let tape = program.evaluate(&[vec![0.0]], None).unwrap();
    assert_eq!(tape.hessian_vector_product(&[1.0, 0.0], &[vec![1.0]], None).unwrap(), vec![vec![2.0]]);
    assert!(matches!(tape.hessian_vector_product(&[0.0, 1.0], &[vec![1.0]], None), Err(HessianError::NonFiniteDerivative { .. })));
}

#[test]
fn directions_and_seeds_refuse_without_mutating_the_primal() {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap(); let z = b.norm_sq(r).unwrap(); let p = b.finish();
    let program = ReverseProgram::new(&p, &[z], limits()).unwrap();
    let tape = program.evaluate(&[vec![1.0, 2.0]], None).unwrap();
    assert!(matches!(tape.hessian_vector_product(&[], &[vec![1.0, 2.0]], None), Err(HessianError::Reverse(ReverseError::SeedCount { .. }))));
    assert!(matches!(tape.hessian_vector_product(&[1.0], &[], None), Err(HessianError::DirectionCount { .. })));
    assert!(matches!(tape.hessian_vector_product(&[1.0], &[vec![1.0]], None), Err(HessianError::DirectionLength { variable: 0, .. })));
    assert!(matches!(tape.hessian_vector_product(&[1.0], &[vec![1.0, f64::NAN]], None), Err(HessianError::NonFiniteDirection { component: 1, .. })));
    assert!(tape.hessian_vector_product(&[f64::NAN], &[vec![1.0, 2.0]], None).is_err());
    assert_eq!(tape.values(), &[5.0]);
    assert_eq!(tape.hessian_vector_product(&[1.0], &[vec![1.0, 2.0]], None).unwrap(), vec![vec![2.0, 4.0]]);
}

#[test]
fn large_vector_uses_linear_scalar_storage_and_ambient_manifold_layout() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 4096 }, Dims::NONE).unwrap();
    let q = b.var("q", Manifold::So3, Dims::NONE).unwrap();
    let xr = b.var_ref(x).unwrap(); let qr = b.var_ref(q).unwrap();
    let nx = b.norm_sq(xr).unwrap(); let nq = b.norm_sq(qr).unwrap(); let z = b.add(nx, nq).unwrap();
    let p = b.finish(); let program = ReverseProgram::new(&p, &[z], limits()).unwrap();
    assert_eq!(program.scalar_slots(), 4103);
    let tape = program.evaluate(&[vec![0.25; 4096], vec![1.0, 0.0, 0.0, 0.0]], None).unwrap();
    let hv = tape.hessian_vector_product(&[1.0], &[vec![1.0; 4096], vec![0.0, 1.0, 2.0, 3.0]], None).unwrap();
    assert_eq!(hv, vec![vec![2.0; 4096], vec![0.0, 2.0, 4.0, 6.0]]);
}

#[test]
fn cancelled_hessian_product_leaves_tape_reusable() {
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let r = b.var_ref(x).unwrap(); let z = b.norm_sq(r).unwrap(); let p = b.finish();
    let program = ReverseProgram::new(&p, &[z], limits()).unwrap();
    let tape = program.evaluate(&[vec![1.0]], None).unwrap();
    let gate = CancelGate::new(); gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        assert!(matches!(tape.hessian_vector_product(&[1.0], &[vec![1.0]], Some(&cx)),
            Err(HessianError::Reverse(ReverseError::Evaluation(fs_opt::OptError::Cancelled)))));
    });
    assert_eq!(tape.hessian_vector_product(&[1.0], &[vec![1.0]], None).unwrap(), vec![vec![2.0]]);
}
