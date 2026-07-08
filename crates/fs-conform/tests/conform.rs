//! Battery for the restriction-map conformance SDK (addendum Proposal 7).
//! Uses concrete matrix converters to exercise each sheaf axiom: functoriality
//! (composition + identity), adjoint consistency (honest vs a lying transpose),
//! and manufactured-solution tolerance honesty (honest vs an understated error
//! model), plus tier assignment and the R6 same-severity rule.

use fs_conform::{
    Composition, ConformanceSuite, Converter, ManufacturedCase, Tier, certify, check_adjoint,
    check_functoriality, check_identity, check_tolerance_honesty,
};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn transpose(a: &[Vec<f64>], rows: usize, cols: usize) -> Vec<Vec<f64>> {
    (0..cols)
        .map(|j| (0..rows).map(|i| a[i][j]).collect())
        .collect()
}

/// A dense linear converter: `apply = A·x`, `adjoint = At·y`.
struct Mtx {
    id: String,
    rows: usize,
    cols: usize,
    a: Vec<Vec<f64>>,
    at: Vec<Vec<f64>>,
    err: f64,
}
impl Mtx {
    /// Honest: the declared adjoint is the true transpose.
    fn honest(id: &str, a: Vec<Vec<f64>>, err: f64) -> Mtx {
        let rows = a.len();
        let cols = a[0].len();
        let at = transpose(&a, rows, cols);
        Mtx {
            id: id.into(),
            rows,
            cols,
            a,
            at,
            err,
        }
    }
    /// Custom (possibly lying) declared adjoint.
    fn with_adjoint(id: &str, a: Vec<Vec<f64>>, at: Vec<Vec<f64>>, err: f64) -> Mtx {
        let rows = a.len();
        let cols = a[0].len();
        Mtx {
            id: id.into(),
            rows,
            cols,
            a,
            at,
            err,
        }
    }
}
impl Converter for Mtx {
    fn id(&self) -> &str {
        &self.id
    }
    fn source_dim(&self) -> usize {
        self.cols
    }
    fn target_dim(&self) -> usize {
        self.rows
    }
    fn apply(&self, x: &[f64]) -> Vec<f64> {
        self.a.iter().map(|row| dot(row, x)).collect()
    }
    fn adjoint(&self, y: &[f64]) -> Vec<f64> {
        self.at.iter().map(|row| dot(row, y)).collect()
    }
    fn declared_error(&self) -> f64 {
        self.err
    }
}

fn manufactured() -> Vec<ManufacturedCase> {
    // true operator T = diag(2, 3).
    vec![
        ManufacturedCase {
            input: vec![1.0, 1.0],
            exact_output: vec![2.0, 3.0],
        },
        ManufacturedCase {
            input: vec![1.0, 0.0],
            exact_output: vec![2.0, 0.0],
        },
    ]
}

#[test]
fn adjoint_consistency_catches_a_lying_transpose() {
    let a = vec![vec![2.0, 1.0], vec![0.0, 3.0]];
    let pairs = vec![
        (vec![1.0, -2.0], vec![0.5, 4.0]),
        (vec![3.0, 1.0], vec![-1.0, 2.0]),
    ];
    // honest transpose passes.
    assert!(check_adjoint(
        &Mtx::honest("h", a.clone(), 1e-9),
        &pairs,
        1e-9
    ));
    // a wrong "adjoint" (identity instead of the transpose) fails.
    let liar = Mtx::with_adjoint("liar", a, vec![vec![1.0, 0.0], vec![0.0, 1.0]], 1e-9);
    assert!(!check_adjoint(&liar, &pairs, 1e-9));
}

#[test]
fn tolerance_honesty_catches_an_understated_error_model() {
    // honest: matrix IS the true operator diag(2,3), tiny declared error.
    let honest = Mtx::honest("honest", vec![vec![2.0, 0.0], vec![0.0, 3.0]], 1e-9);
    let (ok, measured) = check_tolerance_honesty(&honest, &manufactured(), 1e-12);
    assert!(ok && measured < 1e-9);
    // dishonest: matrix diverges from the manufactured truth but declares 1e-9.
    let liar = Mtx::honest("liar", vec![vec![2.5, 0.0], vec![0.0, 3.5]], 1e-9);
    let (ok2, measured2) = check_tolerance_honesty(&liar, &manufactured(), 1e-12);
    assert!(!ok2 && measured2 > 0.5);
}

#[test]
fn functoriality_holds_for_a_real_composition() {
    // g: A->B = diag(2,1); f: B->C = [[1,1],[0,1]]; direct = f*g.
    let g = Mtx::honest("g", vec![vec![2.0, 0.0], vec![0.0, 1.0]], 1e-9);
    let f = Mtx::honest("f", vec![vec![1.0, 1.0], vec![0.0, 1.0]], 1e-9);
    let direct = Mtx::honest("h", vec![vec![2.0, 1.0], vec![0.0, 1.0]], 1e-9); // f·g
    let probes = vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![3.0, -2.0]];
    let comp = Composition {
        after: &f,
        direct: &direct,
        probes: probes.clone(),
    };
    assert!(check_functoriality(&g, &comp, 1e-9));
    // a wrong direct converter breaks functoriality.
    let wrong = Mtx::honest("wrong", vec![vec![1.0, 0.0], vec![0.0, 1.0]], 1e-9);
    let bad = Composition {
        after: &f,
        direct: &wrong,
        probes,
    };
    assert!(!check_functoriality(&g, &bad, 1e-9));
}

#[test]
fn identity_converters_are_recognised() {
    let id = Mtx::honest("id", vec![vec![1.0, 0.0], vec![0.0, 1.0]], 1e-9);
    assert!(check_identity(
        &id,
        &[vec![3.0, 7.0], vec![-1.0, 2.0]],
        1e-12
    ));
    let not_id = Mtx::honest("scale", vec![vec![2.0, 0.0], vec![0.0, 1.0]], 1e-9);
    assert!(!check_identity(&not_id, &[vec![1.0, 1.0]], 1e-12));
}

fn full_suite<'a>() -> ConformanceSuite<'a> {
    ConformanceSuite {
        adjoint_pairs: vec![(vec![1.0, -2.0], vec![0.5, 4.0])],
        manufactured: manufactured(),
        composition: None,
        tolerance: 1e-9,
    }
}

#[test]
fn a_conformant_converter_is_certified_into_a_tier() {
    let good = Mtx::honest("good", vec![vec![2.0, 0.0], vec![0.0, 3.0]], 1e-9);
    let report = certify(&good, &full_suite());
    assert!(report.certified());
    assert_eq!(report.tier, Tier::Gold); // declared 1e-9 <= 1e-6
    assert!(report.functoriality && report.adjoint_consistent && report.tolerance_honest);
    assert!(report.findings.is_empty());
}

#[test]
fn a_dishonest_converter_is_rejected_not_certified() {
    // understated error model -> tolerance honesty fails -> Rejected.
    let liar = Mtx::honest("liar", vec![vec![2.5, 0.0], vec![0.0, 3.5]], 1e-9);
    let report = certify(&liar, &full_suite());
    assert!(!report.certified());
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.tolerance_honest);
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.contains("tolerance honesty"))
    );
}

#[test]
fn a_failing_functoriality_witness_rejects_the_converter() {
    let g = Mtx::honest("g", vec![vec![2.0, 0.0], vec![0.0, 3.0]], 1e-9);
    let f = Mtx::honest("f", vec![vec![1.0, 0.0], vec![0.0, 1.0]], 1e-9); // identity
    let wrong_direct = Mtx::honest("wrong", vec![vec![9.0, 0.0], vec![0.0, 9.0]], 1e-9);
    let mut suite = full_suite();
    suite.composition = Some(Composition {
        after: &f,
        direct: &wrong_direct,
        probes: vec![vec![1.0, 1.0]],
    });
    let report = certify(&g, &suite);
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.functoriality);
}

#[test]
fn tiers_track_the_declared_error_and_r6_severity_is_uniform() {
    // tier by declared error (all otherwise conformant).
    let mk = |err: f64| Mtx::honest("c", vec![vec![2.0, 0.0], vec![0.0, 3.0]], err);
    assert_eq!(certify(&mk(1e-9), &full_suite()).tier, Tier::Gold);
    assert_eq!(certify(&mk(1e-4), &full_suite()).tier, Tier::Silver);
    assert_eq!(certify(&mk(1e-2), &full_suite()).tier, Tier::Bronze);
    // R6: a FIRST-PARTY converter runs the IDENTICAL certify() path and is held
    // to the same bar — a dishonest first-party converter is still Rejected.
    let first_party_liar = Mtx::honest("first-party", vec![vec![5.0, 0.0], vec![0.0, 5.0]], 1e-9);
    assert_eq!(
        certify(&first_party_liar, &full_suite()).tier,
        Tier::Rejected
    );
}
