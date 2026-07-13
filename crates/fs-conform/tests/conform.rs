//! Battery for the restriction-map conformance SDK (addendum Proposal 7).
//! Uses concrete matrix converters to exercise each sheaf axiom: functoriality
//! (composition + identity), adjoint consistency (honest vs a lying transpose),
//! and manufactured-solution tolerance honesty (honest vs an understated error
//! model), plus tier assignment and the R6 same-severity rule.

use fs_conform::{
    Composition, ConformanceSuite, Converter, ManufacturedCase, Tier, certify, check_adjoint,
    check_functoriality, check_identity, check_tolerance_honesty,
};
use fs_propcheck::{Shrink, check};

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
        identity: None, // `good` is a scaling converter, not an identity
        tolerance: 1e-9,
    }
}

#[test]
fn a_false_identity_claim_is_rejected_by_certify() {
    // a converter carrying an identity witness must actually be the identity.
    // A real identity passes the axiom; a scaling map that claims identity fails.
    let probes = vec![vec![3.0, 7.0], vec![-1.0, 2.0]];
    let real_id = Mtx::honest("id", vec![vec![1.0, 0.0], vec![0.0, 1.0]], 1e-9);
    let mut suite = full_suite();
    suite.adjoint_pairs = vec![(vec![1.0, -2.0], vec![1.0, -2.0])]; // idᵀ = id
    suite.manufactured = vec![ManufacturedCase {
        input: vec![2.0, 5.0],
        exact_output: vec![2.0, 5.0],
    }];
    suite.identity = Some(probes.clone());
    let ok = certify(&real_id, &suite);
    assert!(ok.functoriality && ok.certified());

    // the SAME suite applied to a non-identity converter fails the identity axiom.
    let fake = Mtx::honest("fake-id", vec![vec![2.0, 0.0], vec![0.0, 1.0]], 1e-9);
    let bad = certify(&fake, &suite);
    assert!(!bad.functoriality);
    assert_eq!(bad.tier, Tier::Rejected);
    assert!(bad.findings.iter().any(|f| f.contains("identity")));
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

// ---------------------------------------------------------------------------
// G0 property adoption (bead frankensim-4nh8): generated functor and identity
// laws with deterministic shrinking. The fixed cases above remain regression
// pins; these exact small-integer cases cover the space between them.
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct FunctorCase([i64; 10]);

impl Shrink for FunctorCase {
    fn shrink_candidates(&self) -> Vec<Self> {
        let mut out = Vec::new();
        for (index, value) in self.0.iter().enumerate() {
            for candidate in value.shrink_candidates() {
                let mut next = self.clone();
                next.0[index] = candidate;
                out.push(next);
            }
        }
        out
    }
}

fn integer_mtx(id: &str, [a00, a01, a10, a11]: [i64; 4]) -> Mtx {
    Mtx::honest(
        id,
        vec![vec![a00 as f64, a01 as f64], vec![a10 as f64, a11 as f64]],
        0.0,
    )
}

#[test]
fn g0_generated_functoriality_holds_exactly() {
    check(
        "restriction-map-functoriality",
        0xC0F0_4A48_0001,
        512,
        |stream| FunctorCase(std::array::from_fn(|_| stream.int_in(-8, 8))),
        |case| {
            let [f00, f01, f10, f11, g00, g01, g10, g11, p0, p1] = case.0;
            let f = integer_mtx("generated-f", [f00, f01, f10, f11]);
            let g = integer_mtx("generated-g", [g00, g01, g10, g11]);
            let direct = integer_mtx(
                "generated-f-after-g",
                [
                    f00 * g00 + f01 * g10,
                    f00 * g01 + f01 * g11,
                    f10 * g00 + f11 * g10,
                    f10 * g01 + f11 * g11,
                ],
            );
            let composition = Composition {
                after: &f,
                direct: &direct,
                probes: vec![vec![p0 as f64, p1 as f64]],
            };
            check_functoriality(&g, &composition, 0.0)
        },
    );
}

#[test]
fn g0_generated_identity_holds_exactly() {
    check(
        "restriction-map-identity",
        0xC0F0_4A48_0002,
        512,
        |stream| {
            (
                stream.int_in(-1_000_000, 1_000_000),
                stream.int_in(-1_000_000, 1_000_000),
            )
        },
        |&(x, y)| {
            let identity = integer_mtx("generated-identity", [1, 0, 0, 1]);
            check_identity(&identity, &[vec![x as f64, y as f64]], 0.0)
        },
    );
}
