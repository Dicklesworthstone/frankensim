//! FrankenScipy quadrature-oracle evidence for the BEDROCK spectral core.
//!
//! The exact case is a direct-coefficient Chebyshev polynomial KAT. The
//! bounded case independently integrates three smooth functions through
//! `Cheb1::build`/`Cheb1::integral` and FrankenScipy's adaptive GK15 `quad`.
//! Canonical input frames bind every declared option, domain, parameter, and
//! version; output receipts retain both implementations' bits and margins.
//!
//! This is same-binary, finite-fixture G0 agreement and replay evidence. It is
//! not a forward-error certificate, does not run Python SciPy, and makes no
//! claim for discontinuous, singular, improper, infinite, or extreme-domain
//! integrals. It does not cover ODE solvers, performance, cancellation, or
//! establish fresh cross-ISA G5 evidence.

use core::fmt::Write as _;
use std::panic::catch_unwind;

use fs_casebook::{
    CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, SuiteReport, ToleranceSpec, fnv1a64,
};
use fs_cheb::{Cheb1, VERSION as FS_CHEB_VERSION};
use fs_math::VERSION as FS_MATH_VERSION;
use fsci_integrate::{QuadOptions, quad};

const SUITE: &str = "bedrock/fs-cheb-frankenscipy-integrate-oracle-v1";
const ORACLE_VERSION: &str = "fsci-integrate/0.1.0";
const ORACLE_PIN: &str = "9e271fd734465e2b2ff755aa73ea66a7217d619b";
const ORACLE_API: &str = "fsci_integrate::{quad,QuadOptions,QuadResult}:finite-adaptive-gk15:v1";
const FRAME_ENCODING: &str =
    "field=(tag_len:u64le,tag,payload_len:u64le,payload);numbers=le;f64=bits:v1";
const CHEB_CONVENTION: &str = "first-kind-roots;stored-c0-unhalved;Cheb1::build+integral:v1";
const INTEGRAND_POLICY: &str = "exp,sin=fs_math::det;runge=fixed-order-f64;all-functions-pure:v1";

const QUAD_EPSABS: f64 = 1.0e-13;
const QUAD_EPSREL: f64 = 1.0e-13;
const QUAD_LIMIT: usize = 20;
const AGREEMENT_FLOOR: f64 = 4.0e-12;
const ORACLE_ERROR_MULTIPLIER: f64 = 8.0;
const CASEBOOK_ABSOLUTE_CEILING: f64 = 1.0e-11;

const POLY_DOMAIN: [f64; 2] = [-1.0, 1.0];
const POLY_COEFFICIENTS: [f64; 4] = [2.0, 3.0, 0.0, 5.0];
const POLY_POINTS: [f64; 4] = [-1.0, 0.0, 0.5, 1.0];
const POLY_EXPECTED_VALUES: [f64; 4] = [-7.0, 1.0, -2.5, 9.0];
const POLY_EXPECTED_INTEGRAL: f64 = 2.0;

// Filled from the canonical framing code and independently recomputed without
// executing either numerical implementation.
const POLY_FRAME_LEN: usize = 1_271;
const POLY_FRAME_FNV: u64 = 0x2689_7416_3b65_3396;
const ORACLE_FRAME_LEN: usize = 1_607;
const ORACLE_FRAME_FNV: u64 = 0x5382_bf7a_49ab_a61f;
const RED_FRAME_LEN: usize = 2_475;
const RED_FRAME_FNV: u64 = 0x050c_5b1c_a938_949b;

const CORRUPTION_SEED: u64 = 0xF5C1_0019_0000_0007;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FixtureKind {
    Exp,
    Runge,
    Sin20,
}

#[derive(Debug, Clone, Copy)]
struct FixtureSpec {
    kind: FixtureKind,
    id: &'static str,
    domain: [f64; 2],
    max_degree: usize,
    parameter_name: &'static str,
    parameter: Option<f64>,
}

const FIXTURES: [FixtureSpec; 3] = [
    FixtureSpec {
        kind: FixtureKind::Exp,
        id: "exp[-1,1]",
        domain: [-1.0, 1.0],
        max_degree: 64,
        parameter_name: "none",
        parameter: None,
    },
    FixtureSpec {
        kind: FixtureKind::Runge,
        id: "runge25[-1,1]",
        domain: [-1.0, 1.0],
        max_degree: 512,
        parameter_name: "quadratic-scale",
        parameter: Some(25.0),
    },
    FixtureSpec {
        kind: FixtureKind::Sin20,
        id: "sin20x[0,3]",
        domain: [0.0, 3.0],
        max_degree: 512,
        parameter_name: "angular-frequency",
        parameter: Some(20.0),
    },
];

impl FixtureSpec {
    fn evaluate(self, x: f64) -> f64 {
        match self.kind {
            FixtureKind::Exp => fs_math::det::exp(x),
            FixtureKind::Runge => 1.0 / (1.0 + 25.0 * x * x),
            FixtureKind::Sin20 => fs_math::det::sin(20.0 * x),
        }
    }
}

#[derive(Debug, Clone)]
struct Measurement {
    fixture: FixtureSpec,
    cheb_degree: usize,
    cheb_coefficients: Vec<u64>,
    cheb_integral: f64,
    oracle_integral: f64,
    oracle_error: f64,
    oracle_neval: usize,
    oracle_converged: bool,
}

impl Measurement {
    fn delta(&self) -> f64 {
        (self.cheb_integral - self.oracle_integral).abs()
    }

    fn margin(&self) -> f64 {
        self.derived_bound() - self.delta()
    }

    fn derived_bound(&self) -> f64 {
        AGREEMENT_FLOOR + ORACLE_ERROR_MULTIPLIER * self.oracle_error
    }
}

#[derive(Debug, Clone)]
struct Corruption {
    bit: u32,
    canonical: u64,
    corrupted: u64,
    frame: Vec<u8>,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("quadrature Casebook frame lengths fit u64"),
    );
}

fn push_field(bytes: &mut Vec<u8>, tag: &str, payload: &[u8]) {
    push_len(bytes, tag.len());
    bytes.extend_from_slice(tag.as_bytes());
    push_len(bytes, payload.len());
    bytes.extend_from_slice(payload);
}

fn push_text_field(bytes: &mut Vec<u8>, tag: &str, value: &str) {
    push_field(bytes, tag, value.as_bytes());
}

fn push_u32_field(bytes: &mut Vec<u8>, tag: &str, value: u32) {
    push_field(bytes, tag, &value.to_le_bytes());
}

fn push_u64_field(bytes: &mut Vec<u8>, tag: &str, value: u64) {
    push_field(bytes, tag, &value.to_le_bytes());
}

fn push_usize_field(bytes: &mut Vec<u8>, tag: &str, value: usize) {
    push_u64_field(
        bytes,
        tag,
        u64::try_from(value).expect("quadrature Casebook values fit u64"),
    );
}

fn push_f64_field(bytes: &mut Vec<u8>, tag: &str, value: f64) {
    push_u64_field(bytes, tag, value.to_bits());
}

fn push_f64s_field(bytes: &mut Vec<u8>, tag: &str, values: &[f64]) {
    let mut payload = Vec::with_capacity(8 + values.len() * 8);
    push_len(&mut payload, values.len());
    for value in values {
        push_u64(&mut payload, value.to_bits());
    }
    push_field(bytes, tag, &payload);
}

fn common_frame(domain: &str) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_text_field(&mut bytes, "domain", domain);
    push_text_field(&mut bytes, "encoding", FRAME_ENCODING);
    push_u32_field(
        &mut bytes,
        "casebook-record-version",
        CASEBOOK_RECORD_VERSION,
    );
    push_text_field(&mut bytes, "fs-cheb-version", FS_CHEB_VERSION);
    push_text_field(&mut bytes, "fs-math-version", FS_MATH_VERSION);
    push_text_field(&mut bytes, "oracle-version", ORACLE_VERSION);
    push_text_field(&mut bytes, "oracle-pin", ORACLE_PIN);
    push_text_field(&mut bytes, "oracle-api", ORACLE_API);
    push_f64_field(&mut bytes, "quad-epsabs", QUAD_EPSABS);
    push_f64_field(&mut bytes, "quad-epsrel", QUAD_EPSREL);
    push_usize_field(&mut bytes, "quad-limit", QUAD_LIMIT);
    push_f64_field(&mut bytes, "absolute-agreement-floor", AGREEMENT_FLOOR);
    push_f64_field(
        &mut bytes,
        "oracle-error-multiplier",
        ORACLE_ERROR_MULTIPLIER,
    );
    push_f64_field(
        &mut bytes,
        "casebook-absolute-ceiling",
        CASEBOOK_ABSOLUTE_CEILING,
    );
    push_text_field(&mut bytes, "cheb-convention", CHEB_CONVENTION);
    push_text_field(&mut bytes, "integrand-policy", INTEGRAND_POLICY);
    bytes
}

fn polynomial_inputs() -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-cheb:direct-coefficient-polynomial-kat:v1");
    push_f64s_field(&mut bytes, "domain", &POLY_DOMAIN);
    push_f64s_field(
        &mut bytes,
        "stored-c0-unhalved-coefficients",
        &POLY_COEFFICIENTS,
    );
    push_f64s_field(&mut bytes, "evaluation-points", &POLY_POINTS);
    push_f64s_field(
        &mut bytes,
        "expected-evaluation-values",
        &POLY_EXPECTED_VALUES,
    );
    push_f64_field(
        &mut bytes,
        "expected-definite-integral",
        POLY_EXPECTED_INTEGRAL,
    );
    bytes
}

fn fixture_frame(fixture: FixtureSpec) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_text_field(&mut bytes, "fixture-id", fixture.id);
    push_f64s_field(&mut bytes, "domain", &fixture.domain);
    push_usize_field(&mut bytes, "cheb-max-degree", fixture.max_degree);
    push_text_field(&mut bytes, "parameter-name", fixture.parameter_name);
    if let Some(parameter) = fixture.parameter {
        push_f64_field(&mut bytes, "parameter-value", parameter);
    }
    bytes
}

fn oracle_inputs() -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-cheb:frankenscipy-finite-quad-oracle:v1");
    push_usize_field(&mut bytes, "fixture-count", FIXTURES.len());
    for fixture in FIXTURES {
        push_field(&mut bytes, "fixture", &fixture_frame(fixture));
    }
    bytes
}

fn output_receipt(measurements: &[Measurement]) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_text_field(
        &mut bytes,
        "domain",
        "bedrock:fs-cheb:frankenscipy-finite-quad-output-receipt:v1",
    );
    push_u64_field(&mut bytes, "input-frame-fnv1a64", fnv1a64(&oracle_inputs()));
    push_usize_field(&mut bytes, "measurement-count", measurements.len());
    for measurement in measurements {
        let mut row = Vec::new();
        push_text_field(&mut row, "fixture-id", measurement.fixture.id);
        push_usize_field(&mut row, "cheb-degree", measurement.cheb_degree);
        let mut coefficient_payload =
            Vec::with_capacity(8 + measurement.cheb_coefficients.len() * 8);
        push_len(
            &mut coefficient_payload,
            measurement.cheb_coefficients.len(),
        );
        for &coefficient in &measurement.cheb_coefficients {
            push_u64(&mut coefficient_payload, coefficient);
        }
        push_field(&mut row, "cheb-coefficient-bits", &coefficient_payload);
        push_f64_field(&mut row, "cheb-integral", measurement.cheb_integral);
        push_f64_field(&mut row, "oracle-integral", measurement.oracle_integral);
        push_f64_field(&mut row, "oracle-error", measurement.oracle_error);
        push_usize_field(&mut row, "oracle-neval", measurement.oracle_neval);
        push_u32_field(
            &mut row,
            "oracle-converged",
            if measurement.oracle_converged { 1 } else { 0 },
        );
        push_f64_field(&mut row, "absolute-delta", measurement.delta());
        push_f64_field(
            &mut row,
            "derived-agreement-bound",
            measurement.derived_bound(),
        );
        push_f64_field(
            &mut row,
            "casebook-absolute-ceiling",
            CASEBOOK_ABSOLUTE_CEILING,
        );
        push_f64_field(&mut row, "remaining-margin", measurement.margin());
        push_field(&mut bytes, "measurement", &row);
    }
    bytes
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn polynomial_kat_outcome() -> CaseOutcome {
    let first = Cheb1::from_coeffs(POLY_DOMAIN[0], POLY_DOMAIN[1], POLY_COEFFICIENTS.to_vec());
    let replay = Cheb1::from_coeffs(POLY_DOMAIN[0], POLY_DOMAIN[1], POLY_COEFFICIENTS.to_vec());
    let first_values = POLY_POINTS.map(|point| first.eval(point));
    let replay_values = POLY_POINTS.map(|point| replay.eval(point));
    let first_integral = first.integral();
    let replay_integral = replay.integral();

    if bits(&first_values) != bits(&replay_values)
        || first_integral.to_bits() != replay_integral.to_bits()
    {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; first_value_bits={:016x?}; replay_value_bits={:016x?}; first_integral_bits=0x{:016x}; replay_integral_bits=0x{:016x}",
            bits(&first_values),
            bits(&replay_values),
            first_integral.to_bits(),
            replay_integral.to_bits(),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#determinism-class");
    }

    for (index, (&computed, &expected)) in
        first_values.iter().zip(&POLY_EXPECTED_VALUES).enumerate()
    {
        if computed.to_bits() != expected.to_bits() {
            return CaseOutcome::fail(format!(
                "stage=direct-coefficient-evaluation; index={index}; point_bits=0x{:016x}; computed_bits=0x{:016x}; expected_bits=0x{:016x}; coefficients={:016x?}",
                POLY_POINTS[index].to_bits(),
                computed.to_bits(),
                expected.to_bits(),
                bits(&POLY_COEFFICIENTS),
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
        }
    }

    if first_integral.to_bits() != POLY_EXPECTED_INTEGRAL.to_bits() {
        return CaseOutcome::fail(format!(
            "stage=direct-coefficient-integral; computed_bits=0x{:016x}; expected_bits=0x{:016x}; domain_bits={:016x?}; coefficients={:016x?}",
            first_integral.to_bits(),
            POLY_EXPECTED_INTEGRAL.to_bits(),
            bits(&POLY_DOMAIN),
            bits(&POLY_COEFFICIENTS),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
    }

    let mut receipt = Vec::new();
    push_f64s_field(&mut receipt, "evaluation-values", &first_values);
    push_f64_field(&mut receipt, "definite-integral", first_integral);
    CaseOutcome::pass(format!(
        "domain=[-1,1]; coefficients=[2,3,0,5]; evaluation_bits={:016x?}; integral_bits=0x{:016x}; output_receipt_len={}; output_receipt_fnv1a64=0x{:016x}",
        bits(&first_values),
        first_integral.to_bits(),
        receipt.len(),
        fnv1a64(&receipt),
    ))
    .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-cheb/CONTRACT.md#invariants")
}

fn measure_oracles() -> Result<Vec<Measurement>, String> {
    let options = QuadOptions {
        epsabs: QUAD_EPSABS,
        epsrel: QUAD_EPSREL,
        limit: QUAD_LIMIT,
    };
    let mut measurements = Vec::with_capacity(FIXTURES.len());
    for fixture in FIXTURES {
        let cheb = Cheb1::build(
            &|x| fixture.evaluate(x),
            fixture.domain[0],
            fixture.domain[1],
            fixture.max_degree,
        );
        let oracle = quad(
            |x| fixture.evaluate(x),
            fixture.domain[0],
            fixture.domain[1],
            options,
        )
        .map_err(|error| {
            format!(
                "fixture={}; oracle_refusal={error}; domain_bits={:016x?}; epsabs_bits=0x{:016x}; epsrel_bits=0x{:016x}; limit={QUAD_LIMIT}",
                fixture.id,
                bits(&fixture.domain),
                QUAD_EPSABS.to_bits(),
                QUAD_EPSREL.to_bits(),
            )
        })?;
        measurements.push(Measurement {
            fixture,
            cheb_degree: cheb.degree(),
            cheb_coefficients: cheb.coeffs().iter().map(|value| value.to_bits()).collect(),
            cheb_integral: cheb.integral(),
            oracle_integral: oracle.integral,
            oracle_error: oracle.error,
            oracle_neval: oracle.neval,
            oracle_converged: oracle.converged,
        });
    }
    Ok(measurements)
}

#[allow(clippy::too_many_lines)] // One auditable oracle row keeps checks in diagnostic order.
fn quadrature_oracle_outcome() -> CaseOutcome {
    let first = match measure_oracles() {
        Ok(measurements) => measurements,
        Err(error) => {
            return CaseOutcome::fail(format!("stage=first-measurement; {error}"))
                .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
    };
    let replay = match measure_oracles() {
        Ok(measurements) => measurements,
        Err(error) => {
            return CaseOutcome::fail(format!("stage=replay-measurement; {error}"))
                .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
    };
    let first_receipt = output_receipt(&first);
    let replay_receipt = output_receipt(&replay);
    if first_receipt != replay_receipt {
        return CaseOutcome::fail(format!(
            "stage=same-run-output-replay; first_receipt_len={}; replay_receipt_len={}; first_receipt_fnv1a64=0x{:016x}; replay_receipt_fnv1a64=0x{:016x}; first_receipt={}; replay_receipt={}",
            first_receipt.len(),
            replay_receipt.len(),
            fnv1a64(&first_receipt),
            fnv1a64(&replay_receipt),
            hex_bytes(&first_receipt),
            hex_bytes(&replay_receipt),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#determinism-class")
        .with_evidence("constellation.lock:frankenscipy-0.1.0");
    }

    for measurement in &first {
        let non_finite_coefficient = measurement
            .cheb_coefficients
            .iter()
            .position(|&value| !f64::from_bits(value).is_finite());
        if measurement.cheb_degree.checked_add(1) != Some(measurement.cheb_coefficients.len())
            || non_finite_coefficient.is_some()
        {
            return CaseOutcome::fail(format!(
                "stage=cheb-receipt-shape; fixture={}; degree={}; coefficient_count={}; first_non_finite_coefficient={non_finite_coefficient:?}",
                measurement.fixture.id,
                measurement.cheb_degree,
                measurement.cheb_coefficients.len(),
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
        }
        if !measurement.cheb_integral.is_finite()
            || !measurement.oracle_integral.is_finite()
            || !measurement.oracle_error.is_finite()
            || measurement.oracle_error < 0.0
            || measurement.oracle_neval == 0
        {
            return CaseOutcome::fail(format!(
                "stage=finite-oracle-admission; fixture={}; cheb_bits=0x{:016x}; oracle_bits=0x{:016x}; oracle_error_bits=0x{:016x}; oracle_error={}; converged={}; neval={}",
                measurement.fixture.id,
                measurement.cheb_integral.to_bits(),
                measurement.oracle_integral.to_bits(),
                measurement.oracle_error.to_bits(),
                measurement.oracle_error,
                measurement.oracle_converged,
                measurement.oracle_neval,
            ))
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
        if !measurement.oracle_converged {
            return CaseOutcome::fail(format!(
                "stage=oracle-convergence; fixture={}; converged=false; oracle_bits=0x{:016x}; oracle_error_bits=0x{:016x}; oracle_error={}; neval={}; epsabs={QUAD_EPSABS:e}; epsrel={QUAD_EPSREL:e}; limit={QUAD_LIMIT}",
                measurement.fixture.id,
                measurement.oracle_integral.to_bits(),
                measurement.oracle_error.to_bits(),
                measurement.oracle_error,
                measurement.oracle_neval,
            ))
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
        let derived_bound = measurement.derived_bound();
        if !derived_bound.is_finite() || derived_bound > CASEBOOK_ABSOLUTE_CEILING {
            return CaseOutcome::fail(format!(
                "stage=derived-bound-admission; fixture={}; oracle_error_bits=0x{:016x}; oracle_error={}; floor={}; multiplier={}; derived_bound={}; casebook_ceiling={}",
                measurement.fixture.id,
                measurement.oracle_error.to_bits(),
                measurement.oracle_error,
                AGREEMENT_FLOOR,
                ORACLE_ERROR_MULTIPLIER,
                derived_bound,
                CASEBOOK_ABSOLUTE_CEILING,
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#invariants")
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
        if !measurement.delta().is_finite() || measurement.delta() > derived_bound {
            return CaseOutcome::fail(format!(
                "stage=absolute-agreement; fixture={}; domain_bits={:016x?}; max_degree={}; cheb_degree={}; cheb_coefficient_count={}; cheb_bits=0x{:016x}; oracle_bits=0x{:016x}; oracle_error_bits=0x{:016x}; oracle_error={}; absolute_delta={}; derived_bound={}; casebook_ceiling={}; remaining_margin={}; neval={}",
                measurement.fixture.id,
                bits(&measurement.fixture.domain),
                measurement.fixture.max_degree,
                measurement.cheb_degree,
                measurement.cheb_coefficients.len(),
                measurement.cheb_integral.to_bits(),
                measurement.oracle_integral.to_bits(),
                measurement.oracle_error.to_bits(),
                measurement.oracle_error,
                measurement.delta(),
                derived_bound,
                CASEBOOK_ABSOLUTE_CEILING,
                measurement.margin(),
                measurement.oracle_neval,
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#invariants")
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
    }

    let rows = first
        .iter()
        .map(|measurement| {
            format!(
                "{}:degree={},cheb=0x{:016x},oracle=0x{:016x},error={:.3e},delta={:.3e},bound={:.3e},margin={:.3e},neval={}",
                measurement.fixture.id,
                measurement.cheb_degree,
                measurement.cheb_integral.to_bits(),
                measurement.oracle_integral.to_bits(),
                measurement.oracle_error,
                measurement.delta(),
                measurement.derived_bound(),
                measurement.margin(),
                measurement.oracle_neval,
            )
        })
        .collect::<Vec<_>>()
        .join("|");
    CaseOutcome::pass(format!(
        "fixtures={}; {rows}; output_receipt_len={}; output_receipt_fnv1a64=0x{:016x}",
        first.len(),
        first_receipt.len(),
        fnv1a64(&first_receipt),
    ))
    .with_evidence("crates/fs-cheb/CONTRACT.md#invariants")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn corruption_frame(bit: u32, canonical: u64, corrupted: u64) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-cheb:direct-polynomial-reference-corruption:v1");
    push_field(
        &mut bytes,
        "canonical-polynomial-input-frame",
        &polynomial_inputs(),
    );
    push_u64_field(&mut bytes, "corruption-seed", CORRUPTION_SEED);
    push_text_field(&mut bytes, "reference", "direct-coefficient-integral");
    push_u32_field(&mut bytes, "corrupted-bit", bit);
    push_u64_field(&mut bytes, "canonical-integral-bits", canonical);
    push_u64_field(&mut bytes, "corrupted-integral-bits", corrupted);
    bytes
}

fn reconstruct_corruption() -> Corruption {
    let bit = u32::try_from((CORRUPTION_SEED & 0xff) % 52)
        .expect("low-byte-selected mantissa bit fits u32");
    let canonical = POLY_EXPECTED_INTEGRAL.to_bits();
    let corrupted = canonical ^ (1_u64 << bit);
    Corruption {
        bit,
        canonical,
        corrupted,
        frame: corruption_frame(bit, canonical, corrupted),
    }
}

fn corruption_outcome(corruption: Corruption) -> CaseOutcome {
    let computed = Cheb1::from_coeffs(POLY_DOMAIN[0], POLY_DOMAIN[1], POLY_COEFFICIENTS.to_vec())
        .integral()
        .to_bits();
    if computed == corruption.corrupted {
        CaseOutcome::pass("disclosed polynomial-reference corruption was not detected")
    } else {
        CaseOutcome::fail(format!(
            "seed=0x{CORRUPTION_SEED:016x}; reference=direct-coefficient-integral; bit={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}; corrupted_bits=0x{:016x}; input_frame_len={}; input_frame_fnv1a64=0x{:016x}; input_frame={}",
            corruption.bit,
            corruption.canonical,
            corruption.corrupted,
            corruption.frame.len(),
            fnv1a64(&corruption.frame),
            hex_bytes(&corruption.frame),
        ))
        .with_evidence("crates/fs-cheb/tests/frankenscipy_integrate_oracle_casebook.rs#disclosed-corruption")
    }
}

fn run_red_report() -> SuiteReport {
    let corruption = reconstruct_corruption();
    assert_eq!(corruption.frame.len(), RED_FRAME_LEN);
    let digest = fnv1a64(&corruption.frame);
    assert_eq!(digest, RED_FRAME_FNV);
    Suite::new(SUITE)
        .case(
            "disclosed-polynomial-reference-bit-corruption",
            digest,
            ToleranceSpec::Exact,
            move || corruption_outcome(corruption),
        )
        .run()
}

#[test]
fn bedrock_casebook_emits_exact_and_bounded_green_records() {
    let polynomial = polynomial_inputs();
    let oracle = oracle_inputs();
    assert_eq!(polynomial.len(), POLY_FRAME_LEN);
    assert_eq!(fnv1a64(&polynomial), POLY_FRAME_FNV);
    assert_eq!(oracle.len(), ORACLE_FRAME_LEN);
    assert_eq!(fnv1a64(&oracle), ORACLE_FRAME_FNV);

    let report = Suite::new(SUITE)
        .case(
            "direct-coefficient-polynomial-kat",
            POLY_FRAME_FNV,
            ToleranceSpec::Exact,
            polynomial_kat_outcome,
        )
        .case(
            "finite-adaptive-quadrature-agreement",
            ORACLE_FRAME_FNV,
            ToleranceSpec::AbsoluteLe(CASEBOOK_ABSOLUTE_CEILING),
            quadrature_oracle_outcome,
        )
        .run();

    report.assert_green();
    assert_eq!(
        report
            .records
            .iter()
            .map(|record| record.case.as_str())
            .collect::<Vec<_>>(),
        [
            "direct-coefficient-polynomial-kat",
            "finite-adaptive-quadrature-agreement",
        ]
    );
    assert!(report.records.iter().all(|record| {
        record.version == CASEBOOK_RECORD_VERSION && record.pass && !record.evidence.is_empty()
    }));
    assert_eq!(report.records[0].tolerance, "exact");
    assert_eq!(report.records[1].tolerance, "abs<=1e-11");
}

#[test]
fn disclosed_polynomial_reference_corruption_is_stable_and_refused() {
    let first_corruption = reconstruct_corruption();
    let replay_corruption = reconstruct_corruption();
    assert_eq!(first_corruption.bit, replay_corruption.bit);
    assert_eq!(first_corruption.canonical, replay_corruption.canonical);
    assert_eq!(first_corruption.corrupted, replay_corruption.corrupted);
    assert_eq!(first_corruption.frame, replay_corruption.frame);
    assert_eq!(first_corruption.frame.len(), RED_FRAME_LEN);
    assert_eq!(fnv1a64(&first_corruption.frame), RED_FRAME_FNV);
    assert_eq!(first_corruption.bit, 7);
    assert_eq!(
        first_corruption.canonical ^ first_corruption.corrupted,
        1_u64 << first_corruption.bit,
    );

    let first = run_red_report();
    let replay = run_red_report();
    assert!(!first.all_passed());
    assert!(!replay.all_passed());
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("disclosed corruption must produce exactly one red record");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("replayed corruption must produce exactly one red record");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert_eq!(
        first_failure.case,
        "disclosed-polynomial-reference-bit-corruption"
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(
        first_failure
            .details
            .contains("reference=direct-coefficient-integral")
    );
    assert!(first_failure.details.contains("bit=7"));
    assert!(first_failure.details.contains("input_frame="));

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("assert_green must refuse the disclosed polynomial corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("disclosed-polynomial-reference-bit-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
