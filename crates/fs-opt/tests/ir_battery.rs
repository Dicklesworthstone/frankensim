//! Problem-IR battery (7tv.1): graph algebra laws (CSE by construction),
//! round-trip + hash stability, seeded ill-typed rejection with named
//! nodes, differentiability-class routing, exact gradients vs central FD,
//! Riemannian descent fixtures, deterministic mutation fuzz, and the
//! cross-ISA golden hash.

use fs_opt::sexpr::{from_sexpr, problem_hash, to_sexpr};
use fs_opt::{DiffClass, Manifold, OptimizerFamily, Problem, ValidationError, descend};
use fs_qty::Dims;
use std::collections::HashMap;

/// Rosenbrock over Euclidean(2) as an IR graph.
fn rosenbrock() -> Problem {
    let mut p = Problem::new();
    let x = p.variable("x", Manifold::Euclidean(2));
    let x0 = p.component(x, 0);
    let x1 = p.component(x, 1);
    let one = p.scalar(1.0);
    let hundred = p.scalar(100.0);
    let a = p.sub(one, x0).unwrap();
    let a2 = p.powi(a, 2);
    let x0sq = p.powi(x0, 2);
    let b = p.sub(x1, x0sq).unwrap();
    let b2 = p.powi(b, 2);
    let hb2 = p.mul(hundred, b2);
    let obj = p.add(a2, hb2).unwrap();
    p.set_objective(obj);
    p
}

#[test]
fn cse_identity_and_graph_algebra() {
    let mut p = Problem::new();
    let x = p.variable("x", Manifold::Euclidean(1));
    let c0 = p.component(x, 0);
    let c0_again = p.component(x, 0);
    assert_eq!(c0, c0_again, "hash-consing must dedupe identical leaves");
    let s1 = p.powi(c0, 2);
    let s2 = p.powi(c0_again, 2);
    assert_eq!(s1, s2, "identical subexpressions share one node (CSE)");
    // Substitution law: (x+x) evaluates like 2x for all inputs.
    let two = p.scalar(2.0);
    let xx = p.add(c0, c0).unwrap();
    let tx = p.mul(two, c0);
    let mut vals = HashMap::new();
    for t in [-3.0, 0.0, 0.5, 7.25] {
        vals.insert(x, vec![t]);
        p.set_objective(xx);
        let a = p.eval(&vals).unwrap();
        p.set_objective(tx);
        let b = p.eval(&vals).unwrap();
        assert_eq!(a.to_bits(), b.to_bits(), "x+x must equal 2x bitwise at {t}");
    }
}

#[test]
fn round_trip_and_hash_stability() {
    let p = rosenbrock();
    let s1 = to_sexpr(&p);
    let q = from_sexpr(&s1).expect("canonical form parses");
    let s2 = to_sexpr(&q);
    assert_eq!(s1, s2, "serialize → parse → serialize must be identical");
    assert_eq!(problem_hash(&p), problem_hash(&q), "problem hash stable");
    // Evaluation agrees bitwise after the round trip.
    let mut x = HashMap::new();
    x.insert(fs_opt::VarId(0), vec![0.3, -1.7]);
    assert_eq!(
        p.eval(&x).unwrap().to_bits(),
        q.eval(&x).unwrap().to_bits(),
        "round-tripped problem evaluates bitwise-identically"
    );
    println!(
        "{{\"suite\":\"fs-opt\",\"case\":\"roundtrip\",\"verdict\":\"pass\",\"detail\":\"hash {:#018x}\"}}",
        problem_hash(&p)
    );
}

#[test]
fn ill_typed_problems_are_rejected_with_named_nodes() {
    // Dimension mismatch: meters + seconds.
    let mut p = Problem::new();
    let m = p.constant(1.0, Dims([1, 0, 0, 0, 0]));
    let s = p.constant(2.0, Dims([0, 0, 1, 0, 0]));
    match p.add(m, s) {
        Err(ValidationError::DimensionMismatch { op, left, right, .. }) => {
            assert_eq!(op, "add");
            assert_ne!(left, right, "diagnosis must show both unit strings");
        }
        other => panic!("expected DimensionMismatch, got {other:?}"),
    }
    // Transcendental of a dimensioned quantity.
    match p.unary("sin", m) {
        Err(ValidationError::DimensionedTranscendental { op, arg, .. }) => {
            assert_eq!(op, "sin");
            assert!(!arg.is_empty());
        }
        other => panic!("expected DimensionedTranscendental, got {other:?}"),
    }
    // Non-smooth objective fed to a smooth-only optimizer: named node.
    let mut q = Problem::new();
    let x = q.variable("x", Manifold::Euclidean(1));
    let c = q.component(x, 0);
    let zero = q.scalar(0.0);
    let hinge = q.max(c, zero).unwrap();
    q.set_objective(hinge);
    assert_eq!(q.diff_class(hinge), DiffClass::NonSmooth);
    match q.validate_for(OptimizerFamily::SmoothGradient) {
        Err(ValidationError::ClassTooRough { op, found, .. }) => {
            assert_eq!(op, "max");
            assert_eq!(found, DiffClass::NonSmooth);
        }
        other => panic!("expected ClassTooRough, got {other:?}"),
    }
    // The same objective is FINE for subgradient methods.
    q.validate_for(OptimizerFamily::Subgradient).expect("subgradient accepts kinks");
    // PDE placeholder without adjoint: NonDiff, rejected even by
    // subgradient, accepted by derivative-free.
    let mut r = Problem::new();
    let pde = r.pde_residual(42, false);
    r.set_objective(pde);
    assert!(r.validate_for(OptimizerFamily::Subgradient).is_err());
    r.validate_for(OptimizerFamily::DerivativeFree).expect("DFO accepts structure nodes");
    println!(
        "{{\"suite\":\"fs-opt\",\"case\":\"validation\",\"verdict\":\"pass\",\"detail\":\"dim mismatch + transcendental + class routing all rejected with names\"}}"
    );
}

#[test]
fn tampered_serialization_is_rejected_on_parse() {
    // Hand-craft a file whose add mixes meters and seconds: the parser
    // re-runs typed constructors and must refuse it.
    let text = r#"(problem
 (vars
 )
 (nodes
  (const 1.0 (1 0 0 0 0))
  (const 2.0 (0 0 1 0 0))
  (add 0 1)
 )
 (objective 2))"#;
    match from_sexpr(text) {
        Err(fs_opt::sexpr::ParseError::Invalid(ValidationError::DimensionMismatch {
            ..
        })) => {}
        other => panic!("tampered file must be rejected by re-validation: {other:?}"),
    }
}

#[test]
fn gradient_matches_central_differences() {
    let p = rosenbrock();
    let x = fs_opt::VarId(0);
    let mut point = HashMap::new();
    point.insert(x, vec![-0.7, 1.4]);
    let g = p.gradient(&point).unwrap();
    let h = 1e-6;
    for k in 0..2 {
        let mut plus = point.clone();
        plus.get_mut(&x).unwrap()[k] += h;
        let mut minus = point.clone();
        minus.get_mut(&x).unwrap()[k] -= h;
        let fd = (p.eval(&plus).unwrap() - p.eval(&minus).unwrap()) / (2.0 * h);
        assert!(
            (g[&x][k] - fd).abs() < 1e-5 * fd.abs().max(1.0),
            "grad[{k}] {} vs FD {fd}",
            g[&x][k]
        );
    }
}

#[test]
fn riemannian_descent_consumes_manifold_metadata() {
    // Linear objective ⟨x, a⟩ on Sphere(3) with a as CONSTANTS (the toy
    // driver moves every variable, so the fixed direction must be baked
    // in as constant coefficients): minimizer is −a/|a|.
    let target = [0.6f64, -0.64, 0.48]; // |a| = 1 exactly (0.36+0.4096+0.2304)
    let mut p2 = Problem::new();
    let xv = p2.variable("x", Manifold::Sphere(3));
    let mut acc = p2.scalar(0.0);
    for (k, &t) in target.iter().enumerate() {
        let ck = p2.component(xv, u32::try_from(k).unwrap());
        let tk = p2.scalar(t);
        let term = p2.mul(tk, ck);
        acc = p2.add(acc, term).unwrap();
    }
    p2.set_objective(acc);
    let mut start = HashMap::new();
    start.insert(xv, vec![1.0, 0.0, 0.0]);
    let (xf, f) = descend(&p2, start, 0.4, 300).unwrap();
    for (got, want) in xf[&xv].iter().zip(&[-0.6, 0.64, -0.48]) {
        assert!((got - want).abs() < 1e-6, "sphere minimizer: {:?}", xf[&xv]);
    }
    assert!((f + 1.0).abs() < 1e-9, "minimum of <x,a> on the sphere is -|a| = -1: {f}");
    // Norm stays EXACTLY 1 through retraction (unit invariant).
    let nrm: f64 = xf[&xv].iter().map(|t| t * t).sum::<f64>();
    assert!((nrm - 1.0).abs() < 1e-12, "retraction must keep the sphere: {nrm}");
    println!(
        "{{\"suite\":\"fs-opt\",\"case\":\"riemann\",\"verdict\":\"pass\",\"detail\":\"sphere descent to -a with f = {f:.9}\"}}"
    );
}

#[test]
fn deterministic_mutation_fuzz_classification() {
    // LCG-driven random graphs mixing valid and seeded-invalid builds:
    // every constructor either succeeds or returns a STRUCTURED error —
    // never panics, never silently accepts a dimension violation.
    let mut seed = 0xF422_u64;
    let mut lcg = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 33) as usize
    };
    let mut valid = 0usize;
    let mut rejected = 0usize;
    for _ in 0..2000 {
        let mut p = Problem::new();
        let dims_pool = [
            Dims([0, 0, 0, 0, 0]),
            Dims([1, 0, 0, 0, 0]),
            Dims([0, 0, 1, 0, 0]),
        ];
        let a = p.constant(1.5, dims_pool[lcg() % 3]);
        let b = p.constant(-0.5, dims_pool[lcg() % 3]);
        let r = match lcg() % 4 {
            0 => p.add(a, b).map(|_| ()),
            1 => p.sub(a, b).map(|_| ()),
            2 => p.min(a, b).map(|_| ()),
            _ => p.unary("exp", a).map(|_| ()),
        };
        match r {
            Ok(()) => valid += 1,
            Err(
                ValidationError::DimensionMismatch { .. }
                | ValidationError::DimensionedTranscendental { .. },
            ) => rejected += 1,
            Err(other) => panic!("unexpected error class: {other:?}"),
        }
    }
    assert!(valid > 200 && rejected > 200, "fuzz must exercise both: {valid}/{rejected}");
    println!(
        "{{\"suite\":\"fs-opt\",\"case\":\"fuzz\",\"verdict\":\"pass\",\"detail\":\"{valid} valid / {rejected} structured rejections, 0 panics\"}}"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x0; // placeholder: set from first run

#[test]
fn ir_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let p = rosenbrock();
    feed(problem_hash(&p) as f64);
    let x = fs_opt::VarId(0);
    let mut pt = HashMap::new();
    pt.insert(x, vec![-0.7, 1.4]);
    feed(p.eval(&pt).unwrap());
    let g = p.gradient(&pt).unwrap();
    feed(g[&x][0]);
    feed(g[&x][1]);
    // Descent trajectory endpoint bits.
    let mut q = Problem::new();
    let xv = q.variable("x", Manifold::Sphere(3));
    let c0 = q.component(xv, 0);
    let c1 = q.component(xv, 1);
    let two = q.scalar(2.0);
    let t1 = q.mul(two, c1);
    let obj = q.add(c0, t1).unwrap();
    q.set_objective(obj);
    let mut start = HashMap::new();
    start.insert(xv, vec![0.0, 0.0, 1.0]);
    let (xf, f) = descend(&q, start, 0.3, 200).unwrap();
    for &v in &xf[&xv] {
        feed(v);
    }
    feed(f);
    println!(
        "{{\"suite\":\"fs-opt\",\"case\":\"ir-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "IR bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
