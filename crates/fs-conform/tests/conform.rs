//! Battery for the restriction-map conformance SDK (addendum Proposal 7).
//! Uses concrete matrix converters to exercise each sheaf axiom: functoriality
//! (composition + identity), adjoint consistency (honest vs a lying transpose),
//! and manufactured-solution tolerance honesty (honest vs an understated error
//! model), plus tier assignment and the R6 same-severity rule.

use fs_conform::{
    Composition, ConformanceSuite, Converter, ManufacturedCase, Tier, certify, check_adjoint,
    check_functoriality, check_identity, check_tolerance_honesty,
};
use fs_propcheck::{
    Shrink, check,
    metamorphic::{
        RelationCase, RelationObservation, Tolerance, check_relation, conversion_path_independence,
    },
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
        identity: None, // `good` is a scaling converter, not an identity
        tolerance: 1e-9,
    }
}

struct EmptyOutput;

impl Converter for EmptyOutput {
    fn id(&self) -> &str {
        "empty-output"
    }
    fn source_dim(&self) -> usize {
        1
    }
    fn target_dim(&self) -> usize {
        1
    }
    fn apply(&self, _x: &[f64]) -> Vec<f64> {
        Vec::new()
    }
    fn adjoint(&self, _y: &[f64]) -> Vec<f64> {
        Vec::new()
    }
    fn declared_error(&self) -> f64 {
        0.0
    }
}

struct ConstantOutput {
    value: f64,
}

impl Converter for ConstantOutput {
    fn id(&self) -> &str {
        "constant-output"
    }
    fn source_dim(&self) -> usize {
        1
    }
    fn target_dim(&self) -> usize {
        1
    }
    fn apply(&self, _x: &[f64]) -> Vec<f64> {
        vec![self.value]
    }
    fn adjoint(&self, _y: &[f64]) -> Vec<f64> {
        vec![self.value]
    }
    fn declared_error(&self) -> f64 {
        0.0
    }
}

struct ZeroDim;

impl Converter for ZeroDim {
    fn id(&self) -> &str {
        "zero-dimensional"
    }
    fn source_dim(&self) -> usize {
        0
    }
    fn target_dim(&self) -> usize {
        0
    }
    fn apply(&self, _x: &[f64]) -> Vec<f64> {
        Vec::new()
    }
    fn adjoint(&self, _y: &[f64]) -> Vec<f64> {
        Vec::new()
    }
    fn declared_error(&self) -> f64 {
        0.0
    }
}

fn scalar_suite<'a>() -> ConformanceSuite<'a> {
    ConformanceSuite {
        adjoint_pairs: vec![(vec![1.0], vec![1.0])],
        manufactured: vec![ManufacturedCase {
            input: vec![1.0],
            exact_output: vec![1.0],
        }],
        composition: None,
        identity: Some(vec![vec![1.0]]),
        tolerance: 0.0,
    }
}

#[test]
fn structurally_invalid_evidence_cannot_receive_gold() {
    let broken = EmptyOutput;
    let composition = Composition {
        after: &broken,
        direct: &broken,
        probes: vec![vec![1.0]],
    };
    let mut suite = scalar_suite();
    suite.composition = Some(composition);

    let report = certify(&broken, &suite);
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.certified());
    assert!(!report.functoriality);
    assert!(!report.adjoint_consistent);
    assert!(!report.tolerance_honest);
    assert!(report.measured_error.is_infinite());

    let identity = Mtx::honest("scalar-id", vec![vec![1.0]], 0.0);
    assert!(!check_adjoint(&identity, &[(Vec::new(), Vec::new())], 0.0));
    assert!(!check_identity(&identity, &[Vec::new()], 0.0));
    let (honest, measured) = check_tolerance_honesty(
        &identity,
        &[ManufacturedCase {
            input: Vec::new(),
            exact_output: Vec::new(),
        }],
        0.0,
    );
    assert!(!honest && measured.is_infinite());
}

#[test]
fn non_finite_evidence_and_policy_cannot_certify() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let report = certify(&ConstantOutput { value }, &scalar_suite());
        assert_eq!(report.tier, Tier::Rejected);
        assert!(!report.adjoint_consistent);
        assert!(!report.tolerance_honest);
        assert!(report.measured_error.is_infinite());

        let invalid_declared = Mtx::honest("invalid-declared", vec![vec![1.0]], value);
        let report = certify(&invalid_declared, &scalar_suite());
        assert_eq!(report.tier, Tier::Rejected);
        assert!(!report.tolerance_honest);
        assert!(report.measured_error.is_infinite());

        let valid = Mtx::honest("valid", vec![vec![1.0]], 0.0);
        let mut suite = scalar_suite();
        suite.tolerance = value;
        let report = certify(&valid, &suite);
        assert_eq!(report.tier, Tier::Rejected);
        assert!(!report.adjoint_consistent);
        assert!(!report.tolerance_honest);
    }

    let invalid_declared = Mtx::honest("negative-declared", vec![vec![1.0]], -1.0);
    assert_eq!(
        certify(&invalid_declared, &scalar_suite()).tier,
        Tier::Rejected
    );
    let valid = Mtx::honest("valid", vec![vec![1.0]], 0.0);
    let mut suite = scalar_suite();
    suite.tolerance = -1.0;
    assert_eq!(certify(&valid, &suite).tier, Tier::Rejected);

    let overflow = Mtx::honest("overflow", vec![vec![f64::MAX]], 0.0);
    let (honest, measured) = check_tolerance_honesty(
        &overflow,
        &[ManufacturedCase {
            input: vec![f64::MAX],
            exact_output: vec![0.0],
        }],
        0.0,
    );
    assert!(!honest && measured.is_infinite());
}

#[test]
fn floating_point_underflow_cannot_erase_failed_evidence() {
    let tiny = f64::MIN_POSITIVE;
    let wrong_adjoint = Mtx::with_adjoint("underflow", vec![vec![tiny]], vec![vec![0.0]], 0.0);
    assert!(!check_adjoint(
        &wrong_adjoint,
        &[(vec![1.0], vec![tiny])],
        0.0
    ));

    let smallest = f64::from_bits(1);
    let rounded_product_adjoint = Mtx::with_adjoint(
        "rounded-product-adjoint",
        vec![vec![smallest]],
        vec![vec![0.0]],
        0.0,
    );
    assert!(
        !check_adjoint(
            &rounded_product_adjoint,
            &[(vec![1.0], vec![1.25])],
            smallest
        ),
        "1.25 * min-subnormal exceeds the tolerance even when its FMA residual rounds to zero"
    );

    let tiny_output = ConstantOutput { value: tiny };
    assert!(!check_identity(&tiny_output, &[vec![0.0]], 0.0));
    let manufactured = [ManufacturedCase {
        input: vec![0.0],
        exact_output: vec![0.0],
    }];
    let (honest, measured) = check_tolerance_honesty(&tiny_output, &manufactured, 0.0);
    assert!(!honest);
    assert_eq!(measured.to_bits(), tiny.to_bits());

    let rounded_difference = [ManufacturedCase {
        input: vec![1.0],
        exact_output: vec![-f64::from_bits(1)],
    }];
    let one_output = Mtx::honest("rounded-difference", vec![vec![1.0]], 1.0);
    let (honest, measured) = check_tolerance_honesty(&one_output, &rounded_difference, 0.0);
    assert!(
        !honest,
        "DD evidence must retain the exact difference 1 + min-subnormal"
    );
    assert_eq!(measured.to_bits(), 1.0_f64.next_up().to_bits());
    let admitted_one_output = Mtx::honest(
        "admitted-rounded-difference",
        vec![vec![1.0]],
        1.0_f64.next_up(),
    );
    assert!(
        check_tolerance_honesty(&admitted_one_output, &rounded_difference, 0.0).0,
        "the next f64 bound must contain the retained subnormal residual"
    );

    let report = certify(
        &tiny_output,
        &ConformanceSuite {
            adjoint_pairs: vec![(vec![0.0], vec![0.0])],
            manufactured: manufactured.into(),
            composition: None,
            identity: Some(vec![vec![0.0]]),
            tolerance: 0.0,
        },
    );
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.functoriality);
    assert!(!report.tolerance_honest);
    assert_eq!(report.measured_error.to_bits(), tiny.to_bits());

    let absorbed_adjoint = Mtx::with_adjoint(
        "absorbed-adjoint",
        vec![vec![1.0], vec![smallest], vec![-1.0]],
        vec![vec![0.0, 0.0, 0.0]],
        0.0,
    );
    assert!(!check_adjoint(
        &absorbed_adjoint,
        &[(vec![1.0], vec![1.0, 1.0, 1.0])],
        0.0
    ));
    let reverse_absorbed_adjoint = Mtx::with_adjoint(
        "reverse-absorbed-adjoint",
        vec![vec![smallest], vec![1.0], vec![-1.0]],
        vec![vec![0.0, 0.0, 0.0]],
        0.0,
    );
    assert!(!check_adjoint(
        &reverse_absorbed_adjoint,
        &[(vec![1.0], vec![1.0, 1.0, 1.0])],
        0.0
    ));
    let honest_scale_disparate = Mtx::honest(
        "honest-scale-disparate",
        vec![vec![1.0], vec![smallest], vec![-1.0]],
        0.0,
    );
    assert!(
        check_adjoint(
            &honest_scale_disparate,
            &[(vec![1.0], vec![1.0, 1.0, 1.0])],
            smallest
        ),
        "representable roundoff belongs to the declared tolerance, not a structural refusal"
    );
    let ordinary_inexact = Mtx::honest("ordinary-inexact", vec![vec![1.1]], 0.0);
    assert!(
        check_adjoint(&ordinary_inexact, &[(vec![1.1], vec![1.1])], 0.0),
        "an exactly represented TwoProd residual inside the DD rung must not be refused"
    );
    let a = f64::from_bits(963_u64 << 52); // 2^-60
    let b = f64::from_bits((911_u64 << 52) | (1_u64 << 50)); // 5*2^-114
    let partial_low_rounding = Mtx::with_adjoint(
        "partial-low-rounding",
        vec![vec![1.0], vec![a], vec![b], vec![-1.0], vec![-a]],
        vec![vec![0.0; 5]],
        0.0,
    );
    assert!(
        !check_adjoint(
            &partial_low_rounding,
            &[(vec![1.0], vec![1.0; 5])],
            f64::from_bits(911_u64 << 52)
        ),
        "exact residual 5*2^-114 exceeds 2^-112 even if a DD low-component add rounds it down"
    );

    let small = f64::from_bits(523_u64 << 52); // 2^-500; its square is finite.
    let absorbed_distance = Mtx::honest("absorbed-distance", vec![vec![1.0], vec![small]], 1.0);
    let manufactured = vec![ManufacturedCase {
        input: vec![1.0],
        exact_output: vec![0.0, 0.0],
    }];
    let (honest, measured) = check_tolerance_honesty(&absorbed_distance, &manufactured, 0.0);
    assert!(!honest);
    assert_eq!(measured.to_bits(), 1.0_f64.next_up().to_bits());
    let reverse_absorbed_distance = Mtx::honest(
        "reverse-absorbed-distance",
        vec![vec![small], vec![1.0]],
        1.0,
    );
    let (honest, measured) =
        check_tolerance_honesty(&reverse_absorbed_distance, &manufactured, 0.0);
    assert!(!honest);
    assert_eq!(measured.to_bits(), 1.0_f64.next_up().to_bits());
    let absorbed_tolerance = small * small * 0.5;
    let (honest, measured) =
        check_tolerance_honesty(&absorbed_distance, &manufactured, absorbed_tolerance);
    assert!(
        honest,
        "the exact squared-domain decision must retain a bound that ordinary f64 addition absorbs"
    );
    assert_eq!(measured.to_bits(), 1.0_f64.next_up().to_bits());

    let square_residual_underflow = f64::from_bits((486_u64 << 52) | (1_u64 << 51));
    let rounded_square = Mtx::honest(
        "rounded-square",
        vec![vec![1.0], vec![square_residual_underflow]],
        1.0,
    );
    let (honest, measured) =
        check_tolerance_honesty(&rounded_square, &manufactured, f64::from_bits(1));
    assert!(
        !honest && measured.is_infinite(),
        "(3*2^-538)^2 is 2.25 min-subnormals, not its inward-rounded 2-min-subnormal f64 product"
    );

    let admitted_rounding = Mtx::honest(
        "admitted-rounding",
        vec![vec![1.0], vec![small]],
        1.0_f64.next_up(),
    );
    let (honest, measured) = check_tolerance_honesty(&admitted_rounding, &manufactured, 0.0);
    assert!(
        honest,
        "the next representable declared bound contains the DD norm"
    );
    assert_eq!(measured.to_bits(), 1.0_f64.next_up().to_bits());
    let report = certify(
        &absorbed_distance,
        &ConformanceSuite {
            adjoint_pairs: vec![(vec![1.0], vec![0.0, 0.0])],
            manufactured,
            composition: None,
            identity: None,
            tolerance: 0.0,
        },
    );
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.tolerance_honest);
    assert_eq!(report.measured_error.to_bits(), 1.0_f64.next_up().to_bits());

    let outside_dd_rung = Mtx::honest(
        "outside-dd-rung",
        vec![
            vec![1.0],
            vec![f64::from_bits(993_u64 << 52)], // 2^-30
            vec![f64::from_bits(923_u64 << 52)], // 2^-100
        ],
        f64::MAX,
    );
    let (honest, measured) = check_tolerance_honesty(
        &outside_dd_rung,
        &[ManufacturedCase {
            input: vec![1.0],
            exact_output: vec![0.0, 0.0, 0.0],
        }],
        0.0,
    );
    assert!(!honest && measured.is_infinite());
}

#[test]
fn dimension_topology_and_finite_overflow_fail_closed() {
    let base = Mtx::honest("base", vec![vec![1.0, 0.0], vec![0.0, 1.0]], 0.0);
    let scalar = Mtx::honest("scalar", vec![vec![1.0]], 0.0);
    let two_to_one = Mtx::honest("two-to-one", vec![vec![1.0, 0.0]], 0.0);
    let one_to_two = Mtx::honest("one-to-two", vec![vec![1.0], vec![0.0]], 0.0);
    let probe = vec![vec![1.0, 0.0]];

    for composition in [
        Composition {
            after: &scalar,
            direct: &two_to_one,
            probes: probe.clone(),
        },
        Composition {
            after: &base,
            direct: &one_to_two,
            probes: probe.clone(),
        },
        Composition {
            after: &two_to_one,
            direct: &base,
            probes: probe,
        },
    ] {
        assert!(!check_functoriality(&base, &composition, 0.0));
    }

    let far = ConstantOutput { value: f64::MAX };
    let (honest, measured) = check_tolerance_honesty(
        &far,
        &[ManufacturedCase {
            input: vec![0.0],
            exact_output: vec![-f64::MAX],
        }],
        0.0,
    );
    assert!(!honest && measured.is_infinite());

    let overflowing_bound = Mtx::honest("overflowing-bound", vec![vec![1.0]], f64::MAX);
    let mut suite = scalar_suite();
    suite.tolerance = f64::MAX;
    let report = certify(&overflowing_bound, &suite);
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.tolerance_honest);
    assert!(report.measured_error.is_infinite());
}

#[test]
fn evidence_free_suites_and_empty_optional_witnesses_cannot_certify() {
    let identity = Mtx::honest("scalar-id", vec![vec![1.0]], 0.0);
    assert!(!check_adjoint(&identity, &[], 0.0));
    assert_eq!(
        check_tolerance_honesty(&identity, &[], 0.0),
        (false, f64::INFINITY)
    );
    assert!(!check_functoriality(
        &identity,
        &Composition {
            after: &identity,
            direct: &identity,
            probes: Vec::new(),
        },
        0.0
    ));
    assert!(!check_identity(&identity, &[], 0.0));

    let report = certify(&identity, &ConformanceSuite::new(0.0));
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.adjoint_consistent);
    assert!(!report.tolerance_honest);

    let mut empty_composition = scalar_suite();
    empty_composition.composition = Some(Composition {
        after: &identity,
        direct: &identity,
        probes: Vec::new(),
    });
    let report = certify(&identity, &empty_composition);
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.functoriality);

    let mut empty_identity = scalar_suite();
    empty_identity.identity = Some(Vec::new());
    let report = certify(&identity, &empty_identity);
    assert_eq!(report.tier, Tier::Rejected);
    assert!(!report.functoriality);
}

#[test]
fn nonempty_zero_dimensional_evidence_remains_valid() {
    let zero = ZeroDim;
    let pairs = vec![(Vec::new(), Vec::new())];
    let manufactured = vec![ManufacturedCase {
        input: Vec::new(),
        exact_output: Vec::new(),
    }];
    let probes = vec![Vec::new()];
    let composition = Composition {
        after: &zero,
        direct: &zero,
        probes: probes.clone(),
    };

    assert!(check_adjoint(&zero, &pairs, 0.0));
    assert_eq!(
        check_tolerance_honesty(&zero, &manufactured, 0.0),
        (true, 0.0)
    );
    assert!(check_functoriality(&zero, &composition, 0.0));
    assert!(check_identity(&zero, &probes, 0.0));

    let report = certify(
        &zero,
        &ConformanceSuite {
            adjoint_pairs: pairs,
            manufactured,
            composition: Some(composition),
            identity: Some(probes),
            tolerance: 0.0,
        },
    );
    assert_eq!(report.tier, Tier::Gold);
    assert!(report.certified());
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
fn tiers_track_the_admitted_bound_and_r6_severity_is_uniform() {
    // Tier by the admitted bound (declared error + suite tolerance), with the
    // default suite tolerance far below the three declarations here.
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

#[test]
fn suite_tolerance_is_charged_to_the_tier_and_cannot_launder_gold() {
    let exact = Mtx::honest("exact", vec![vec![2.0, 0.0], vec![0.0, 3.0]], 0.0);

    let mut suite = full_suite();
    suite.tolerance = 1e-2;
    let bronze = certify(&exact, &suite);
    assert!(bronze.certified() && bronze.tolerance_honest);
    assert_eq!(bronze.measured_error, 0.0);
    assert_eq!(bronze.tier, Tier::Bronze);

    suite.tolerance = 1e-4;
    assert_eq!(certify(&exact, &suite).tier, Tier::Silver);

    suite.tolerance = 1e-6;
    assert_eq!(certify(&exact, &suite).tier, Tier::Gold);

    suite.tolerance = 1e-6_f64.next_up();
    assert_eq!(
        certify(&exact, &suite).tier,
        Tier::Silver,
        "one representable step beyond the Gold bound must not round back into Gold"
    );

    // Direct laundering witness: the zero map is one full unit away from the
    // manufactured truth, but a suite tolerance of one explicitly admits that
    // error. It may receive only the corresponding coarse tier, never Gold from
    // its misleading zero declaration.
    let poor = Mtx::honest("poor-but-admitted", vec![vec![0.0]], 0.0);
    let poor_report = certify(
        &poor,
        &ConformanceSuite {
            adjoint_pairs: vec![(vec![1.0], vec![1.0])],
            manufactured: vec![ManufacturedCase {
                input: vec![1.0],
                exact_output: vec![1.0],
            }],
            composition: None,
            identity: None,
            tolerance: 1.0,
        },
    );
    assert!(poor_report.certified() && poor_report.tolerance_honest);
    assert_eq!(poor_report.measured_error, 1.0);
    assert_eq!(poor_report.tier, Tier::Bronze);
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

// ---------------------------------------------------------------------------
// G3 metamorphic adoption (bead frankensim-2uce): the route transform selects
// two real Converter::apply paths for the same manufactured map. The fixed and
// generated G0 functoriality pins above remain independent regressions.
// ---------------------------------------------------------------------------

fn nonvacuous_functor_case(stream: &mut fs_propcheck::Stream) -> FunctorCase {
    let mut values = std::array::from_fn(|_| stream.int_in(-8, 8));
    if values[..4].iter().all(|value| *value == 0) {
        values[0] = 1;
    }
    if values[4..8].iter().all(|value| *value == 0) {
        values[4] = 1;
    }
    if values[8..].iter().all(|value| *value == 0) {
        values[8] = 1;
    }
    FunctorCase(values)
}

fn apply_conversion_path(case: &FunctorCase, composed: bool) -> Vec<f64> {
    let [f00, f01, f10, f11, g00, g01, g10, g11, p0, p1] = case.0;
    let f = integer_mtx("g3-f", [f00, f01, f10, f11]);
    let g = integer_mtx("g3-g", [g00, g01, g10, g11]);
    let probe = [p0 as f64, p1 as f64];
    if composed {
        f.apply(&g.apply(&probe))
    } else {
        let direct = integer_mtx(
            "g3-f-after-g",
            [
                f00 * g00 + f01 * g10,
                f00 * g01 + f01 * g11,
                f10 * g00 + f11 * g10,
                f10 * g01 + f11 * g11,
            ],
        );
        direct.apply(&probe)
    }
}

fn exact_path_observation(
    direct: &[f64],
    composed: &[f64],
    tolerance: Tolerance,
) -> RelationObservation {
    if direct.is_empty() || direct.len() != composed.len() {
        return RelationObservation::new(
            -1.0,
            "direct and composed conversion paths must return the same nonempty dimension",
        );
    }
    let margin = direct
        .iter()
        .zip(composed)
        .map(|(reference, candidate)| {
            let reference = if *reference == 0.0 { 0.0 } else { *reference };
            let candidate = if *candidate == 0.0 { 0.0 } else { *candidate };
            tolerance.evaluate_scalar(reference, candidate).margin()
        })
        .fold(0.0_f64, f64::min);
    RelationObservation::new(
        margin,
        "direct and composed Converter::apply paths agree componentwise",
    )
}

#[test]
fn g3_generated_conversion_paths_agree_exactly() {
    let relation = conversion_path_independence(
        "restriction-map-direct-vs-composed",
        Tolerance::Exact,
        |input: &(FunctorCase, i64), route: &i64| (input.0.clone(), *route),
        |direct: &Vec<f64>, composed: &Vec<f64>, _route: &i64, tolerance| {
            exact_path_observation(direct, composed, tolerance)
        },
    );
    let operator = |input: &(FunctorCase, i64)| apply_conversion_path(&input.0, input.1 != 0);

    check_relation(
        "fs-conform/converter-apply",
        0xC0F0_4A48_0003,
        512,
        |stream| RelationCase::new((nonvacuous_functor_case(stream), 0_i64), 1_i64),
        &operator,
        &relation,
    );
}
