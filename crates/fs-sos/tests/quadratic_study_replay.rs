//! G5 full-result replay for the production quadratic SOS constructor
//! (`7tv.21.33`).
//!
//! The fixture binds every bit returned by `certify_quadratic(1, -2, 3)` and
//! every derived public verdict for the exact identity
//! `x^2 - 2x + 3 = (x - 1)^2 + 2`. A retained schema-v1 root and independent
//! in-process run must both reproduce the full canonical frame. A disclosed seed
//! selects one square coefficient and one mantissa bit; the mutation is refused
//! while stale, refused against the retained reference after resealing, rejected
//! by the production certificate semantics, emitted as stable fs-obs evidence,
//! and caught by the merge gate.
//!
//! This proves only one exact dyadic quadratic and one mutation lane. It does
//! not claim arbitrary-quadratic exactness, tolerance-based global soundness,
//! general SOS search, exhaustive tamper resistance, cryptographic identity,
//! cross-ISA equality, cancellation, persistence, or performance.

use fs_obs::ident::{IDENT_SCHEMA_VERSION, IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_sos::{Poly, SosCertificate, certify_quadratic, square};
use std::panic::catch_unwind;

const SUITE: &str = "fs-sos/quadratic-study-replay";
const CASE: &str = "exact-dyadic-quadratic-full-result";
const RED_CASE: &str = "seeded-square-coefficient-corruption";
const FS_SOS_VERSION: &str = env!("CARGO_PKG_VERSION");

const A: f64 = 1.0;
const B: f64 = -2.0;
const C: f64 = 3.0;
const EXPECTED_LOWER_BOUND: f64 = 2.0;
const VERIFY_TOLERANCE: f64 = 1.0e-9;
const CERTIFIED_RADIUS: f64 = 16.0;

const MUTATION_SEED: u64 = 0x5053_4f53_0000_0118;
const MUTATION_TARGET: &str = "certificate.squares[0].coefficients[1]";
const MUTATION_SELECTOR: &str = "square=0;coefficient=(seed>>8)&1;mantissa-bit=seed&0x1f";

// Retained schema-v1 roots independently derived from the exact typed frames.
// A change requires a semantic identity review, not an automatic regeneration.
const EXPECTED_FIXTURE_ROOT: u64 = 0xf693_af68_89d0_321b;
const EXPECTED_FIXTURE_CANONICAL_BYTES: usize = 1_219;
const EXPECTED_RESULT_ROOT: u64 = 0x043e_5c9f_f9ac_b99b;
const EXPECTED_RESULT_CANONICAL_BYTES: usize = 1_295;

const _: () = assert!(A.to_bits() == 0x3ff0_0000_0000_0000);
const _: () = assert!(B.to_bits() == 0xc000_0000_0000_0000);
const _: () = assert!(C.to_bits() == 0x4008_0000_0000_0000);
const _: () = assert!(EXPECTED_LOWER_BOUND.to_bits() == 0x4000_0000_0000_0000);
const _: () = assert!((MUTATION_SEED & 0x1f) == 24);
const _: () = assert!(((MUTATION_SEED >> 8) & 1) == 1);

#[derive(Debug, Clone, PartialEq, Eq)]
struct CertificateRecord {
    inputs: [u64; 3],
    polynomial: Vec<u64>,
    certificate_present: bool,
    lower_bound: u64,
    squares: Vec<Vec<u64>>,
    residual: u64,
    coefficient_identity_pass: bool,
    global_bound: Option<u64>,
    radius_bound: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    record: CertificateRecord,
    result: ReplayIdentity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SemanticError {
    InputMismatch {
        index: usize,
        expected: u64,
        found: u64,
    },
    PolynomialLength {
        expected: usize,
        found: usize,
    },
    PolynomialCoefficient {
        index: usize,
        expected: u64,
        found: u64,
    },
    CertificateAbsent,
    LowerBound {
        expected: u64,
        found: u64,
    },
    SquareCount {
        expected: usize,
        found: usize,
    },
    SquareCoefficientCount {
        square: usize,
        expected: usize,
        found: usize,
    },
    SquareCoefficient {
        square: usize,
        coefficient: usize,
        expected: u64,
        found: u64,
    },
    Residual {
        expected: u64,
        found: u64,
    },
    CoefficientIdentityVerdict {
        expected: bool,
        found: bool,
    },
    GlobalBound {
        expected: Option<u64>,
        found: Option<u64>,
    },
    RadiusBound {
        expected: Option<u64>,
        found: Option<u64>,
    },
    PolynomialIdentity,
    AnalyticMinimum {
        expected: u64,
        found: u64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MergeRefusal {
    Admission(AdmissionError),
    Semantics(SemanticError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    square: usize,
    coefficient: usize,
    mantissa_bit: u32,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SeededCorruption {
    run: StudyRun,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    semantic_error: SemanticError,
    first_mismatch: String,
    production_verify: bool,
    production_global_bound: Option<u64>,
    production_radius_bound: Option<u64>,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn input_bits() -> [u64; 3] {
    [A.to_bits(), B.to_bits(), C.to_bits()]
}

fn polynomial_bits() -> Vec<u64> {
    vec![C.to_bits(), B.to_bits(), A.to_bits()]
}

fn polynomial() -> Poly {
    Poly::new(vec![C, B, A])
}

fn fixture_identity() -> ReplayIdentity {
    IdentityBuilder::new("fs-sos-quadratic-study-fixture-v1")
        .str("function", "fs_sos::certify_quadratic")
        .str("polynomial", "a*x^2+b*x+c")
        .str("coefficient-units", "dimensionless")
        .str("value-units", "dimensionless")
        .f64_bits("a", A)
        .f64_bits("b", B)
        .f64_bits("c", C)
        .f64_bits("coefficient-identity-tolerance", VERIFY_TOLERANCE)
        .f64_bits("certified-radius", CERTIFIED_RADIUS)
        .str("mutation-selector", MUTATION_SELECTOR)
        .u64("mutation-seed", MUTATION_SEED)
        .str("mutation-target", MUTATION_TARGET)
        .str("execution-context", "single-threaded-direct-test-no-Cx")
        .str("fs-sos-version", FS_SOS_VERSION)
        .str("fs-ivl-version", fs_ivl::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .u64("fs-obs-identity-schema", u64::from(IDENT_SCHEMA_VERSION))
        .u64("fs-obs-wire-schema", u64::from(fs_obs::SCHEMA_VERSION))
        .u64(
            "fs-obs-event-content-identity-schema",
            u64::from(fs_obs::EVENT_CONTENT_IDENTITY_VERSION),
        )
        .str(
            "no-claims",
            "arbitrary-quadratic-exactness;tolerance-as-global-theorem;general-SOS-search;multivariate-completeness;exhaustive-tamper-resistance;cryptographic-authenticity;cross-ISA;Cx;concurrency;checkpoint;persistence;external-oracle;performance",
        )
        .finish()
}

fn result_identity(fixture: &ReplayIdentity, record: &CertificateRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-sos-quadratic-study-result-v1")
        .child("fixture", fixture)
        .u64("input-count", usize_u64(record.inputs.len()));
    for (index, &bits) in record.inputs.iter().enumerate() {
        builder = builder
            .u64("input-index", usize_u64(index))
            .f64_bits("input-value", f64::from_bits(bits));
    }
    builder = builder.u64(
        "polynomial-coefficient-count",
        usize_u64(record.polynomial.len()),
    );
    for (index, &bits) in record.polynomial.iter().enumerate() {
        builder = builder
            .u64("polynomial-coefficient-index", usize_u64(index))
            .f64_bits("polynomial-coefficient", f64::from_bits(bits));
    }
    builder = builder
        .flag("certificate-present", record.certificate_present)
        .f64_bits("lower-bound", f64::from_bits(record.lower_bound))
        .u64("square-count", usize_u64(record.squares.len()));
    for (square_index, coefficients) in record.squares.iter().enumerate() {
        builder = builder
            .u64("square-index", usize_u64(square_index))
            .u64("square-coefficient-count", usize_u64(coefficients.len()));
        for (coefficient, &bits) in coefficients.iter().enumerate() {
            builder = builder
                .u64("square-coefficient-index", usize_u64(coefficient))
                .f64_bits("square-coefficient", f64::from_bits(bits));
        }
    }
    builder = builder
        .f64_bits("residual", f64::from_bits(record.residual))
        .flag(
            "coefficient-identity-pass",
            record.coefficient_identity_pass,
        )
        .flag("global-bound-present", record.global_bound.is_some());
    if let Some(bits) = record.global_bound {
        builder = builder.f64_bits("global-bound", f64::from_bits(bits));
    }
    builder = builder.flag("radius-bound-present", record.radius_bound.is_some());
    if let Some(bits) = record.radius_bound {
        builder = builder.f64_bits("radius-bound", f64::from_bits(bits));
    }
    builder.finish()
}

fn record_certificate(p: &Poly, certificate: &SosCertificate) -> CertificateRecord {
    CertificateRecord {
        inputs: input_bits(),
        polynomial: p.coeffs().iter().map(|value| value.to_bits()).collect(),
        certificate_present: true,
        lower_bound: certificate.lower_bound.to_bits(),
        squares: certificate
            .squares
            .iter()
            .map(|q| q.coeffs().iter().map(|value| value.to_bits()).collect())
            .collect(),
        residual: certificate.residual(p).to_bits(),
        coefficient_identity_pass: certificate.verify(p, VERIFY_TOLERANCE),
        global_bound: certificate.certified_bound_global(p).map(f64::to_bits),
        radius_bound: certificate
            .certified_bound_on(p, CERTIFIED_RADIUS)
            .map(f64::to_bits),
    }
}

fn run_study() -> StudyRun {
    let p = polynomial();
    let certificate = certify_quadratic(A, B, C).expect("finite positive quadratic");
    let record = record_certificate(&p, &certificate);
    let fixture = fixture_identity();
    let result = result_identity(&fixture, &record);
    StudyRun {
        fixture,
        record,
        result,
    }
}

fn polynomial_from_record(record: &CertificateRecord) -> Poly {
    Poly::new(
        record
            .polynomial
            .iter()
            .copied()
            .map(f64::from_bits)
            .collect(),
    )
}

fn certificate_from_record(record: &CertificateRecord) -> Option<SosCertificate> {
    record.certificate_present.then(|| SosCertificate {
        squares: record
            .squares
            .iter()
            .map(|coefficients| {
                Poly::new(coefficients.iter().copied().map(f64::from_bits).collect())
            })
            .collect(),
        lower_bound: f64::from_bits(record.lower_bound),
    })
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let computed = result_identity(&run.fixture, &run.record);
    if run.result == computed {
        Ok(())
    } else {
        Err(AdmissionError::PayloadIdentityMismatch {
            declared: run.result.root(),
            computed: computed.root(),
        })
    }
}

fn admit_against(run: &StudyRun, reference: &ReplayIdentity) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if &run.result == reference {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: reference.root(),
            found: run.result.root(),
        })
    }
}

#[allow(clippy::too_many_lines)] // Every public result surface is checked explicitly.
fn validate_semantics(record: &CertificateRecord) -> Result<(), SemanticError> {
    for (index, (&found, expected)) in record.inputs.iter().zip(input_bits()).enumerate() {
        if found != expected {
            return Err(SemanticError::InputMismatch {
                index,
                expected,
                found,
            });
        }
    }
    let expected_polynomial = polynomial_bits();
    if record.polynomial.len() != expected_polynomial.len() {
        return Err(SemanticError::PolynomialLength {
            expected: expected_polynomial.len(),
            found: record.polynomial.len(),
        });
    }
    for (index, (&found, &expected)) in record
        .polynomial
        .iter()
        .zip(&expected_polynomial)
        .enumerate()
    {
        if found != expected {
            return Err(SemanticError::PolynomialCoefficient {
                index,
                expected,
                found,
            });
        }
    }
    if !record.certificate_present {
        return Err(SemanticError::CertificateAbsent);
    }
    if record.lower_bound != EXPECTED_LOWER_BOUND.to_bits() {
        return Err(SemanticError::LowerBound {
            expected: EXPECTED_LOWER_BOUND.to_bits(),
            found: record.lower_bound,
        });
    }
    if record.squares.len() != 1 {
        return Err(SemanticError::SquareCount {
            expected: 1,
            found: record.squares.len(),
        });
    }
    let expected_square = [(-1.0_f64).to_bits(), 1.0_f64.to_bits()];
    if record.squares[0].len() != expected_square.len() {
        return Err(SemanticError::SquareCoefficientCount {
            square: 0,
            expected: expected_square.len(),
            found: record.squares[0].len(),
        });
    }
    for (coefficient, (&found, &expected)) in
        record.squares[0].iter().zip(&expected_square).enumerate()
    {
        if found != expected {
            return Err(SemanticError::SquareCoefficient {
                square: 0,
                coefficient,
                expected,
                found,
            });
        }
    }

    let exact_zero = 0.0_f64.to_bits();
    if record.residual != exact_zero {
        return Err(SemanticError::Residual {
            expected: exact_zero,
            found: record.residual,
        });
    }
    if !record.coefficient_identity_pass {
        return Err(SemanticError::CoefficientIdentityVerdict {
            expected: true,
            found: record.coefficient_identity_pass,
        });
    }
    let exact_bound = Some(EXPECTED_LOWER_BOUND.to_bits());
    if record.global_bound != exact_bound {
        return Err(SemanticError::GlobalBound {
            expected: exact_bound,
            found: record.global_bound,
        });
    }
    if record.radius_bound != exact_bound {
        return Err(SemanticError::RadiusBound {
            expected: exact_bound,
            found: record.radius_bound,
        });
    }

    let p = polynomial_from_record(record);
    let certificate = certificate_from_record(record).ok_or(SemanticError::CertificateAbsent)?;
    let derived_residual = certificate.residual(&p).to_bits();
    if record.residual != derived_residual {
        return Err(SemanticError::Residual {
            expected: derived_residual,
            found: record.residual,
        });
    }
    let derived_verify = certificate.verify(&p, VERIFY_TOLERANCE);
    if record.coefficient_identity_pass != derived_verify {
        return Err(SemanticError::CoefficientIdentityVerdict {
            expected: derived_verify,
            found: record.coefficient_identity_pass,
        });
    }
    let derived_global = certificate.certified_bound_global(&p).map(f64::to_bits);
    if record.global_bound != derived_global {
        return Err(SemanticError::GlobalBound {
            expected: derived_global,
            found: record.global_bound,
        });
    }
    let derived_radius = certificate
        .certified_bound_on(&p, CERTIFIED_RADIUS)
        .map(f64::to_bits);
    if record.radius_bound != derived_radius {
        return Err(SemanticError::RadiusBound {
            expected: derived_radius,
            found: record.radius_bound,
        });
    }

    let reconstructed = certificate
        .squares
        .iter()
        .fold(Poly::constant(certificate.lower_bound), |sum, q| {
            sum.add(&square(q))
        });
    if reconstructed != p {
        return Err(SemanticError::PolynomialIdentity);
    }
    let at_minimizer = p.eval(1.0).to_bits();
    if at_minimizer != EXPECTED_LOWER_BOUND.to_bits() {
        return Err(SemanticError::AnalyticMinimum {
            expected: EXPECTED_LOWER_BOUND.to_bits(),
            found: at_minimizer,
        });
    }
    Ok(())
}

fn first_record_mismatch(left: &CertificateRecord, right: &CertificateRecord) -> Option<String> {
    if left.inputs != right.inputs {
        return Some(format!(
            "inputs:{:016x?}!={:016x?}",
            left.inputs, right.inputs
        ));
    }
    if left.polynomial != right.polynomial {
        return Some(format!(
            "polynomial:{:016x?}!={:016x?}",
            left.polynomial, right.polynomial
        ));
    }
    if left.certificate_present != right.certificate_present {
        return Some(format!(
            "certificate_present:{}!={}",
            left.certificate_present, right.certificate_present
        ));
    }
    if left.lower_bound != right.lower_bound {
        return Some(format!(
            "certificate.lower_bound:0x{:016x}!=0x{:016x}",
            left.lower_bound, right.lower_bound
        ));
    }
    if left.squares.len() != right.squares.len() {
        return Some(format!(
            "certificate.squares.length:{}!={}",
            left.squares.len(),
            right.squares.len()
        ));
    }
    for (square_index, (a, b)) in left.squares.iter().zip(&right.squares).enumerate() {
        if a.len() != b.len() {
            return Some(format!(
                "certificate.squares[{square_index}].coefficients.length:{}!={}",
                a.len(),
                b.len()
            ));
        }
        for (coefficient, (&a, &b)) in a.iter().zip(b).enumerate() {
            if a != b {
                return Some(format!(
                    "certificate.squares[{square_index}].coefficients[{coefficient}]:0x{a:016x}!=0x{b:016x}"
                ));
            }
        }
    }
    if left.residual != right.residual {
        return Some(format!(
            "residual:0x{:016x}!=0x{:016x}",
            left.residual, right.residual
        ));
    }
    if left.coefficient_identity_pass != right.coefficient_identity_pass {
        return Some(format!(
            "coefficient_identity_pass:{}!={}",
            left.coefficient_identity_pass, right.coefficient_identity_pass
        ));
    }
    if left.global_bound != right.global_bound {
        return Some(format!(
            "global_bound:{:?}!={:?}",
            left.global_bound, right.global_bound
        ));
    }
    if left.radius_bound != right.radius_bound {
        return Some(format!(
            "radius_bound:{:?}!={:?}",
            left.radius_bound, right.radius_bound
        ));
    }
    None
}

fn mutation_from_seed(seed: u64) -> (usize, usize, u32) {
    let coefficient = usize::try_from((seed >> 8) & 1).expect("selector fits usize");
    let mantissa_bit = u32::try_from(seed & 0x1f).expect("selector fits u32");
    (0, coefficient, mantissa_bit)
}

fn exact_one_bit_delta(
    reference: &CertificateRecord,
    mutant: &CertificateRecord,
    mutation: Mutation,
) -> bool {
    let mut repaired = mutant.clone();
    repaired.squares[mutation.square][mutation.coefficient] = mutation.before;
    &repaired == reference
        && (mutation.before ^ mutation.after) == 1_u64 << mutation.mantissa_bit
        && (mutation.before ^ mutation.after).count_ones() == 1
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
    let (square, coefficient, mantissa_bit) = mutation_from_seed(MUTATION_SEED);
    let mut run = reference.clone();
    let before = run.record.squares[square][coefficient];
    let after = before ^ (1_u64 << mantissa_bit);
    run.record.squares[square][coefficient] = after;
    let mutation = Mutation {
        seed: MUTATION_SEED,
        square,
        coefficient,
        mantissa_bit,
        before,
        after,
    };

    let stale_error = validate_payload(&run).expect_err("stale result identity must refuse");
    let semantic_error =
        validate_semantics(&run.record).expect_err("mutated square must fail semantics");
    let mutant_certificate =
        certificate_from_record(&run.record).expect("mutation retains certificate shape");
    let p = polynomial_from_record(&run.record);
    let production_verify = mutant_certificate.verify(&p, VERIFY_TOLERANCE);
    let production_global_bound = mutant_certificate
        .certified_bound_global(&p)
        .map(f64::to_bits);
    let production_radius_bound = mutant_certificate
        .certified_bound_on(&p, CERTIFIED_RADIUS)
        .map(f64::to_bits);

    run.result = result_identity(&run.fixture, &run.record);
    let reference_error = admit_against(&run, &reference.result)
        .expect_err("resealed mutation must not match retained reference");
    let first_mismatch = first_record_mismatch(&reference.record, &run.record)
        .expect("seeded mutation changes the retained record");
    SeededCorruption {
        run,
        mutation,
        stale_error,
        reference_error,
        semantic_error,
        first_mismatch,
        production_verify,
        production_global_bound,
        production_radius_bound,
    }
}

fn merge_refusal(run: &StudyRun, reference: &ReplayIdentity) -> Result<(), MergeRefusal> {
    admit_against(run, reference).map_err(MergeRefusal::Admission)?;
    validate_semantics(&run.record).map_err(MergeRefusal::Semantics)?;
    Ok(())
}

fn assert_mergeable(run: &StudyRun, reference: &ReplayIdentity) {
    if let Err(error) = merge_refusal(run, reference) {
        panic!(
            "{RED_CASE}: seed=0x{MUTATION_SEED:016x}; target={MUTATION_TARGET}; refusal={error:?}"
        );
    }
}

fn corruption_detail(reference: &StudyRun, corruption: &SeededCorruption) -> String {
    format!(
        "fixture={}; reference={}; mutant={}; seed=0x{:016x}; selector={}; target={}; bit={}; before=0x{:016x}; after=0x{:016x}; stale_gate={:?}; reference_gate={:?}; semantic_gate={:?}; production_verify={}; production_global={:?}; production_radius={:?}; first_mismatch={}",
        reference.fixture.hex(),
        reference.result.hex(),
        corruption.run.result.hex(),
        corruption.mutation.seed,
        MUTATION_SELECTOR,
        MUTATION_TARGET,
        corruption.mutation.mantissa_bit,
        corruption.mutation.before,
        corruption.mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.semantic_error,
        corruption.production_verify,
        corruption.production_global_bound,
        corruption.production_radius_bound,
        corruption.first_mismatch,
    )
}

fn failure_event(detail: &str) -> Event {
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail: detail.to_string(),
            seed: MUTATION_SEED,
        },
        None,
    )
}

fn assert_event_admits(event: &Event) {
    fs_obs::lint_failure_record(event).expect("event retains required replay context");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("event uses the fs-obs wire schema");
    let receipt = event.content_identity_receipt();
    event
        .admit_content_identity(&receipt)
        .expect("fresh event identity admits exact content");
}

fn emit_green_receipt(run: &StudyRun) {
    let square = &run.record.squares[0];
    let global = run.record.global_bound.expect("exact global bound");
    let radius = run.record.radius_bound.expect("finite radius bound");
    let mut emitter = Emitter::new(SUITE, CASE);
    let event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "quadratic-sos-full-result-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"result_identity\":\"{}\",",
                    "\"function\":\"fs_sos::certify_quadratic\",",
                    "\"inputs\":{{\"a\":\"0x{:016x}\",\"b\":\"0x{:016x}\",\"c\":\"0x{:016x}\"}},",
                    "\"result\":{{\"present\":true,\"lower_bound\":\"0x{:016x}\",",
                    "\"squares\":[[\"0x{:016x}\",\"0x{:016x}\"]],",
                    "\"residual\":\"0x{:016x}\",\"coefficient_identity_pass\":{},",
                    "\"global_bound\":\"0x{:016x}\",\"radius_bound\":\"0x{:016x}\"}},",
                    "\"mutation\":{{\"seed\":\"0x{MUTATION_SEED:016x}\",",
                    "\"selector\":\"{MUTATION_SELECTOR}\",\"target\":\"{MUTATION_TARGET}\"}},",
                    "\"identity_schema\":{},\"wire_schema\":{},",
                    "\"event_content_identity_schema\":{},",
                    "\"versions\":{{\"fs_sos\":\"{}\",\"fs_ivl\":\"{}\",",
                    "\"fs_math\":\"{}\",\"fs_obs\":\"{}\"}},",
                    "\"execution\":\"single-threaded-direct-test-no-Cx\",",
                    "\"no_claims\":[\"arbitrary-quadratic-exactness\",",
                    "\"tolerance-as-global-theorem\",\"general-SOS-search\",",
                    "\"exhaustive-tamper-resistance\",\"cryptographic-authenticity\",",
                    "\"cross-ISA\",\"cancellation\",\"persistence\",\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.result.hex(),
                run.record.inputs[0],
                run.record.inputs[1],
                run.record.inputs[2],
                run.record.lower_bound,
                square[0],
                square[1],
                run.record.residual,
                run.record.coefficient_identity_pass,
                global,
                radius,
                IDENT_SCHEMA_VERSION,
                fs_obs::SCHEMA_VERSION,
                fs_obs::EVENT_CONTENT_IDENTITY_VERSION,
                FS_SOS_VERSION,
                fs_ivl::VERSION,
                fs_math::VERSION,
                fs_obs::VERSION,
            ),
        },
        None,
    );
    assert_event_admits(&event);
    println!("{}", event.to_jsonl());
}

fn emit_green_verdict(run: &StudyRun) -> Event {
    let detail = format!(
        "fixture={}; result={}; complete_result=bit-exact; global_bound=0x{:016x}; radius_bound=0x{:016x}",
        run.fixture.hex(),
        run.result.hex(),
        run.record.global_bound.expect("global bound"),
        run.record.radius_bound.expect("radius bound"),
    );
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    let event = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail,
            seed: 0,
        },
        None,
    );
    assert_event_admits(&event);
    println!("{}", event.to_jsonl());
    event
}

fn exercise_seeded_corruption(original: &StudyRun, replay: &StudyRun) {
    let first = seeded_corruption(original);
    let second = seeded_corruption(replay);
    assert_eq!(first, second, "seeded red state must replay exactly");
    assert_eq!(first.mutation.square, 0);
    assert_eq!(first.mutation.coefficient, 1);
    assert_eq!(first.mutation.mantissa_bit, 24);
    assert!(f64::from_bits(first.mutation.after).is_finite());
    assert!(exact_one_bit_delta(
        &original.record,
        &first.run.record,
        first.mutation
    ));
    assert_eq!(
        first.first_mismatch.split(':').next(),
        Some(MUTATION_TARGET)
    );
    assert!(matches!(
        first.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if declared == original.result.root() && computed == first.run.result.root()
    ));
    assert!(matches!(
        first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.result.root() && found == first.run.result.root()
    ));
    assert!(matches!(
        first.semantic_error,
        SemanticError::SquareCoefficient {
            square: 0,
            coefficient: 1,
            expected,
            found,
        } if expected == first.mutation.before && found == first.mutation.after
    ));
    assert!(
        !first.production_verify,
        "the production coefficient diagnostic must reject the one-bit mutant"
    );
    assert_eq!(
        first.production_global_bound, None,
        "nonzero non-constant exact residual must block a global claim"
    );
    assert!(
        first
            .production_radius_bound
            .is_some_and(|bits| f64::from_bits(bits) <= EXPECTED_LOWER_BOUND),
        "radius-scoped admission may degrade, never overstate, the true minimum"
    );

    let first_detail = corruption_detail(original, &first);
    let second_detail = corruption_detail(replay, &second);
    assert_eq!(first_detail, second_detail);
    let first_event = failure_event(&first_detail);
    let second_event = failure_event(&second_detail);
    assert_event_admits(&first_event);
    assert_event_admits(&second_event);
    assert_eq!(first_event, second_event);
    let first_event_identity = first_event.content_identity();
    let second_event_identity = second_event.content_identity();
    assert_eq!(
        first_event_identity.canonical_bytes(),
        second_event_identity.canonical_bytes(),
        "red evidence must replay byte-for-byte"
    );
    println!("{}", first_event.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first.run, &original.result))
        .expect_err("merge gate must refuse the resealed mutation");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(MUTATION_TARGET));
    assert!(message.contains("ReferenceIdentityMismatch"));
}

#[test]
fn quadratic_certificate_replays_and_seeded_false_certificate_is_refused() {
    let original = run_study();
    let replay = run_study();

    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));
    assert_eq!(validate_semantics(&original.record), Ok(()));
    assert_eq!(validate_semantics(&replay.record), Ok(()));
    assert_eq!(original.record.residual, 0.0_f64.to_bits());
    assert!(original.record.coefficient_identity_pass);
    assert_eq!(
        original.record.global_bound,
        Some(EXPECTED_LOWER_BOUND.to_bits())
    );
    assert_eq!(
        original.record.radius_bound,
        Some(EXPECTED_LOWER_BOUND.to_bits())
    );
    assert_eq!(
        first_record_mismatch(&original.record, &replay.record),
        None,
        "complete public result must replay exactly"
    );
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.result, replay.result);
    assert_eq!(original.fixture.version(), IDENT_SCHEMA_VERSION);
    assert_eq!(original.fixture.root(), EXPECTED_FIXTURE_ROOT);
    assert_eq!(
        original.fixture.canonical_bytes().len(),
        EXPECTED_FIXTURE_CANONICAL_BYTES
    );
    assert_eq!(original.result.version(), IDENT_SCHEMA_VERSION);
    assert_eq!(original.result.root(), EXPECTED_RESULT_ROOT);
    assert_eq!(
        original.result.canonical_bytes().len(),
        EXPECTED_RESULT_CANONICAL_BYTES
    );
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete result identity must replay byte-for-byte"
    );

    let q = Poly::new(vec![-1.0, 1.0]);
    assert_eq!(
        square(&q).add(&Poly::constant(EXPECTED_LOWER_BOUND)),
        polynomial(),
        "independent exact oracle: p=(x-1)^2+2"
    );
    assert_eq!(
        polynomial().eval(1.0).to_bits(),
        EXPECTED_LOWER_BOUND.to_bits()
    );
    assert_mergeable(&original, &original.result);
    emit_green_receipt(&original);
    let green = emit_green_verdict(&original);
    assert!(matches!(
        green.kind,
        EventKind::ConformanceCase { pass: true, .. }
    ));
    exercise_seeded_corruption(&original, &replay);
}
