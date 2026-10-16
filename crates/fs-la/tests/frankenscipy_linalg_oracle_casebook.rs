//! FrankenScipy dense-solve oracle evidence for the BEDROCK factorization core.
//!
//! The exact case checks power-of-two diagonal systems through both public
//! factorization solve paths. The bounded case compares one general LU solve
//! and one positive-definite Cholesky solve with FrankenScipy's finite dense
//! solver. Canonical frames bind the fixtures, policies, versions, and oracle
//! pin; output receipts retain every solution bit and disclosed error field.
//!
//! This is finite-fixture G0 agreement and same-build replay evidence. It is
//! not Python SciPy evidence, a conditioning survey, a forward-error
//! certificate, or fresh cross-ISA G5 evidence. It does not cover singular or
//! ill-conditioned systems, QR, SVD, eigenproblems, least squares, performance,
//! cancellation, or external-checkout attestation.

use core::fmt::Write as _;
use std::panic::catch_unwind;

use fs_casebook::{
    CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, SuiteReport, ToleranceSpec, fnv1a64,
};
use fs_la::VERSION as FS_LA_VERSION;
use fs_la::factor::{FACTOR_BIT_SEMANTICS_VERSION, cholesky, lu};
use fsci_linalg::{DecompOptions, cho_factor, cho_solve, lu_factor, lu_solve};

const SUITE: &str = "bedrock/fs-la-frankenscipy-linalg-oracle-v1";
const ORACLE_VERSION: &str = "fsci-linalg/0.1.0";
const ORACLE_PIN: &str = "9e271fd734465e2b2ff755aa73ea66a7217d619b";
const ORACLE_API: &str = "fsci_linalg::{lu_factor,lu_solve,cho_factor,cho_solve,DecompOptions,SolveResult}:strict-finite:v1";
const FRAME_ENCODING: &str =
    "field=(tag_len:u64le,tag,payload_len:u64le,payload);numbers=le;f64=bits:v1";
const PRODUCTION_POLICY: &str =
    "fs_la::factor::{lu,cholesky};row-major;in-place-rhs;scaled-inf-residual:v1";

const N: usize = 3;
const SOLUTION_CEILING: f64 = 1.0e-12;
const RESIDUAL_CEILING: f64 = 1.0e-13;
const CORRUPTION_SEED: u64 = 0xF5C1_0020_0000_0007;

const EXACT_MATRIX: [[f64; N]; N] = [[4.0, 0.0, 0.0], [0.0, 16.0, 0.0], [0.0, 0.0, 64.0]];
const EXACT_RHS: [f64; N] = [8.0, -32.0, 128.0];
const EXACT_SOLUTION: [f64; N] = [2.0, -2.0, 2.0];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FactorKind {
    Lu,
    Cholesky,
}

impl FactorKind {
    const fn name(self) -> &'static str {
        match self {
            Self::Lu => "lu",
            Self::Cholesky => "cholesky",
        }
    }

    const fn oracle_factorization(self) -> &'static str {
        match self {
            Self::Lu => "lu-factor",
            Self::Cholesky => "cholesky-factor",
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Fixture {
    id: &'static str,
    kind: FactorKind,
    matrix: [[f64; N]; N],
    rhs: [f64; N],
    expected: [f64; N],
}

const FIXTURES: [Fixture; 2] = [
    Fixture {
        id: "general-binary-rational-3x3",
        kind: FactorKind::Lu,
        matrix: [[4.0, 1.0, -1.0], [2.0, 6.0, 1.0], [-1.0, 1.0, 5.0]],
        rhs: [4.0, 1.0, 12.0],
        expected: [2.0, -1.0, 3.0],
    },
    Fixture {
        id: "spd-binary-rational-3x3",
        kind: FactorKind::Cholesky,
        matrix: [[4.0, 1.0, 1.0], [1.0, 3.0, 0.0], [1.0, 0.0, 2.0]],
        rhs: [5.0, -5.0, 7.0],
        expected: [1.0, -2.0, 3.0],
    },
];

// Independently recomputed from the canonical framing code without executing
// either numerical implementation.
const EXACT_FRAME_LEN: usize = 1_084;
const EXACT_FRAME_FNV: u64 = 0x3cf8_0bdb_f20f_9f72;
const ORACLE_FRAME_LEN: usize = 1_667;
const ORACLE_FRAME_FNV: u64 = 0xd2b6_4149_bffc_2ce2;
const RED_FRAME_LEN: usize = 2_168;
const RED_FRAME_FNV: u64 = 0x7728_e77c_fa08_be27;

#[derive(Debug, Clone)]
struct Measurement {
    fixture: Fixture,
    production_solution: Vec<f64>,
    oracle_solution: Vec<f64>,
    production_residual_vector: Vec<f64>,
    oracle_residual_vector: Vec<f64>,
    production_residual: f64,
    oracle_residual: f64,
    oracle_warning: String,
    oracle_backward_error: Option<f64>,
    oracle_certificate: String,
}

impl Measurement {
    fn max_cross_delta(&self) -> f64 {
        max_abs_delta(&self.production_solution, &self.oracle_solution)
    }

    fn max_production_reference_delta(&self) -> f64 {
        max_abs_delta(&self.production_solution, &self.fixture.expected)
    }

    fn max_oracle_reference_delta(&self) -> f64 {
        max_abs_delta(&self.oracle_solution, &self.fixture.expected)
    }

    fn remaining_margin(&self) -> f64 {
        SOLUTION_CEILING
            - self
                .max_cross_delta()
                .max(self.max_production_reference_delta())
                .max(self.max_oracle_reference_delta())
    }
}

#[derive(Debug, Clone)]
struct Corruption {
    index: usize,
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
        u64::try_from(value).expect("dense oracle Casebook frame lengths fit u64"),
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
        u64::try_from(value).expect("dense oracle Casebook values fit u64"),
    );
}

fn push_bool_field(bytes: &mut Vec<u8>, tag: &str, value: bool) {
    push_u32_field(bytes, tag, if value { 1 } else { 0 });
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
    push_text_field(&mut bytes, "fs-la-version", FS_LA_VERSION);
    push_u32_field(
        &mut bytes,
        "factor-bit-semantics-version",
        FACTOR_BIT_SEMANTICS_VERSION,
    );
    push_text_field(&mut bytes, "production-policy", PRODUCTION_POLICY);
    push_text_field(&mut bytes, "oracle-version", ORACLE_VERSION);
    push_text_field(&mut bytes, "oracle-pin", ORACLE_PIN);
    push_text_field(&mut bytes, "oracle-api", ORACLE_API);
    push_text_field(&mut bytes, "oracle-mode", "strict");
    push_bool_field(&mut bytes, "oracle-check-finite", true);
    push_f64_field(&mut bytes, "solution-ceiling", SOLUTION_CEILING);
    push_f64_field(&mut bytes, "residual-ceiling", RESIDUAL_CEILING);
    bytes
}

fn flatten(matrix: &[[f64; N]; N]) -> Vec<f64> {
    matrix.iter().flatten().copied().collect()
}

fn exact_inputs() -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-la:exact-diagonal-factor-solve:v1");
    push_usize_field(&mut bytes, "dimension", N);
    push_text_field(&mut bytes, "factor-order", "lu,cholesky");
    push_f64s_field(&mut bytes, "matrix-row-major", &flatten(&EXACT_MATRIX));
    push_f64s_field(&mut bytes, "rhs", &EXACT_RHS);
    push_f64s_field(&mut bytes, "expected-solution", &EXACT_SOLUTION);
    bytes
}

fn fixture_frame(fixture: Fixture) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_text_field(&mut bytes, "fixture-id", fixture.id);
    push_text_field(&mut bytes, "production-factor", fixture.kind.name());
    push_text_field(
        &mut bytes,
        "oracle-factorization",
        fixture.kind.oracle_factorization(),
    );
    push_usize_field(&mut bytes, "dimension", N);
    push_f64s_field(&mut bytes, "matrix-row-major", &flatten(&fixture.matrix));
    push_f64s_field(&mut bytes, "rhs", &fixture.rhs);
    push_f64s_field(&mut bytes, "expected-solution", &fixture.expected);
    bytes
}

fn oracle_inputs() -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-la:frankenscipy-dense-solve-oracle:v1");
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
        "bedrock:fs-la:frankenscipy-dense-solve-output-receipt:v1",
    );
    push_u64_field(&mut bytes, "input-frame-fnv1a64", fnv1a64(&oracle_inputs()));
    push_usize_field(&mut bytes, "measurement-count", measurements.len());
    for measurement in measurements {
        let mut row = Vec::new();
        push_text_field(&mut row, "fixture-id", measurement.fixture.id);
        push_f64s_field(
            &mut row,
            "production-solution",
            &measurement.production_solution,
        );
        push_f64s_field(&mut row, "oracle-solution", &measurement.oracle_solution);
        push_f64s_field(
            &mut row,
            "production-signed-residual-vector",
            &measurement.production_residual_vector,
        );
        push_f64s_field(
            &mut row,
            "oracle-signed-residual-vector",
            &measurement.oracle_residual_vector,
        );
        push_f64_field(
            &mut row,
            "production-residual",
            measurement.production_residual,
        );
        push_f64_field(&mut row, "oracle-residual", measurement.oracle_residual);
        push_text_field(&mut row, "oracle-warning", &measurement.oracle_warning);
        push_bool_field(
            &mut row,
            "oracle-backward-error-present",
            measurement.oracle_backward_error.is_some(),
        );
        if let Some(backward_error) = measurement.oracle_backward_error {
            push_f64_field(&mut row, "oracle-backward-error", backward_error);
        }
        push_f64_field(&mut row, "max-cross-delta", measurement.max_cross_delta());
        push_f64_field(
            &mut row,
            "max-production-reference-delta",
            measurement.max_production_reference_delta(),
        );
        push_f64_field(
            &mut row,
            "max-oracle-reference-delta",
            measurement.max_oracle_reference_delta(),
        );
        push_f64_field(&mut row, "remaining-margin", measurement.remaining_margin());
        push_text_field(
            &mut row,
            "oracle-certificate",
            &measurement.oracle_certificate,
        );
        push_field(&mut bytes, "measurement", &row);
    }
    bytes
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn max_abs_delta(lhs: &[f64], rhs: &[f64]) -> f64 {
    lhs.iter()
        .zip(rhs)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max)
}

fn residual_vector(matrix: &[[f64; N]; N], solution: &[f64], rhs: &[f64]) -> Vec<f64> {
    matrix
        .iter()
        .zip(rhs)
        .map(|(row, expected)| {
            let computed = row
                .iter()
                .zip(solution)
                .fold(0.0, |sum, (coefficient, value)| {
                    coefficient.mul_add(*value, sum)
                });
            computed - expected
        })
        .collect()
}

fn scaled_inf_residual(matrix: &[[f64; N]; N], solution: &[f64], rhs: &[f64]) -> f64 {
    let residual = residual_vector(matrix, solution, rhs)
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
    let matrix_norm = matrix
        .iter()
        .map(|row| row.iter().map(|value| value.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    let solution_norm = solution.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let rhs_norm = rhs.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let scale = matrix_norm.mul_add(solution_norm, rhs_norm);
    if scale == 0.0 {
        residual
    } else {
        residual / scale
    }
}

fn production_solve(
    kind: FactorKind,
    matrix: &[[f64; N]; N],
    rhs: &[f64; N],
) -> Result<Vec<f64>, String> {
    let flat = flatten(matrix);
    let mut solution = rhs.to_vec();
    match kind {
        FactorKind::Lu => lu(&flat, N)
            .map_err(|error| format!("production-lu-refusal={error}"))?
            .solve(&mut solution),
        FactorKind::Cholesky => cholesky(&flat, N)
            .map_err(|error| format!("production-cholesky-refusal={error}"))?
            .solve(&mut solution),
    }
    Ok(solution)
}

fn measure_fixture(fixture: Fixture) -> Result<Measurement, String> {
    let production_solution = production_solve(fixture.kind, &fixture.matrix, &fixture.rhs)?;
    let matrix_rows = fixture
        .matrix
        .iter()
        .map(|row| row.to_vec())
        .collect::<Vec<_>>();
    let options = DecompOptions::default();
    let oracle = match fixture.kind {
        FactorKind::Lu => {
            let factor = lu_factor(&matrix_rows, options).map_err(|error| {
                format!("oracle-lu-factor-refusal={error}; fixture={}", fixture.id)
            })?;
            lu_solve(&factor, &fixture.rhs).map_err(|error| {
                format!("oracle-lu-solve-refusal={error}; fixture={}", fixture.id)
            })?
        }
        FactorKind::Cholesky => {
            let factor = cho_factor(&matrix_rows, options).map_err(|error| {
                format!(
                    "oracle-cholesky-factor-refusal={error}; fixture={}",
                    fixture.id
                )
            })?;
            if factor.dimension() != N {
                return Err(format!(
                    "oracle-cholesky-dimension={}; expected={N}; fixture={}",
                    factor.dimension(),
                    fixture.id,
                ));
            }
            cho_solve(&factor, &fixture.rhs).map_err(|error| {
                format!(
                    "oracle-cholesky-solve-refusal={error}; fixture={}",
                    fixture.id
                )
            })?
        }
    };
    let oracle_warning = format!("{:?}", oracle.warning);
    let oracle_backward_error = oracle.backward_error;
    let oracle_certificate = format!("{:?}", oracle.certificate);
    let oracle_solution = oracle.x;
    let production_residual_vector =
        residual_vector(&fixture.matrix, &production_solution, &fixture.rhs);
    let oracle_residual_vector = residual_vector(&fixture.matrix, &oracle_solution, &fixture.rhs);
    Ok(Measurement {
        fixture,
        production_residual: scaled_inf_residual(
            &fixture.matrix,
            &production_solution,
            &fixture.rhs,
        ),
        oracle_residual: scaled_inf_residual(&fixture.matrix, &oracle_solution, &fixture.rhs),
        production_solution,
        oracle_solution,
        production_residual_vector,
        oracle_residual_vector,
        oracle_warning,
        oracle_backward_error,
        oracle_certificate,
    })
}

fn measure_all() -> Result<Vec<Measurement>, String> {
    FIXTURES.into_iter().map(measure_fixture).collect()
}

fn exact_factor_outcome() -> CaseOutcome {
    let mut receipts = Vec::new();
    for kind in [FactorKind::Lu, FactorKind::Cholesky] {
        let first = match production_solve(kind, &EXACT_MATRIX, &EXACT_RHS) {
            Ok(solution) => solution,
            Err(error) => return CaseOutcome::fail(format!("stage=first; kind={kind:?}; {error}")),
        };
        let replay = match production_solve(kind, &EXACT_MATRIX, &EXACT_RHS) {
            Ok(solution) => solution,
            Err(error) => {
                return CaseOutcome::fail(format!("stage=replay; kind={kind:?}; {error}"));
            }
        };
        if bits(&first) != bits(&replay) || bits(&first) != bits(&EXACT_SOLUTION) {
            return CaseOutcome::fail(format!(
                "stage=exact-diagonal-solve; kind={kind:?}; first_bits={:016x?}; replay_bits={:016x?}; expected_bits={:016x?}",
                bits(&first),
                bits(&replay),
                bits(&EXACT_SOLUTION),
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics");
        }
        push_text_field(&mut receipts, "factor", kind.name());
        push_f64s_field(&mut receipts, "solution", &first);
    }
    CaseOutcome::pass(format!(
        "factors=2; solution_bits={:016x?}; output_receipt_len={}; output_receipt_fnv1a64=0x{:016x}",
        bits(&EXACT_SOLUTION),
        receipts.len(),
        fnv1a64(&receipts),
    ))
    .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-la/CONTRACT.md#invariants")
}

#[allow(clippy::too_many_lines)]
fn oracle_outcome() -> CaseOutcome {
    let input_frame = oracle_inputs();
    let input_frame_digest = fnv1a64(&input_frame);
    let input_frame_hex = hex_bytes(&input_frame);
    let first = match measure_all() {
        Ok(measurements) => measurements,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=first-measurement; {error}; input_frame_len={}; input_frame_fnv1a64=0x{input_frame_digest:016x}; input_frame={input_frame_hex}",
                input_frame.len(),
            ))
                .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
    };
    let replay = match measure_all() {
        Ok(measurements) => measurements,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=replay-measurement; {error}; input_frame_len={}; input_frame_fnv1a64=0x{input_frame_digest:016x}; input_frame={input_frame_hex}",
                input_frame.len(),
            ))
                .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
    };
    let first_receipt = output_receipt(&first);
    let replay_receipt = output_receipt(&replay);
    if first_receipt != replay_receipt {
        return CaseOutcome::fail(format!(
            "stage=same-run-output-replay; first_len={}; replay_len={}; first_fnv1a64=0x{:016x}; replay_fnv1a64=0x{:016x}; first={}; replay={}",
            first_receipt.len(),
            replay_receipt.len(),
            fnv1a64(&first_receipt),
            fnv1a64(&replay_receipt),
            hex_bytes(&first_receipt),
            hex_bytes(&replay_receipt),
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#determinism-class")
        .with_evidence("constellation.lock:frankenscipy-0.1.0");
    }

    for measurement in &first {
        let finite = measurement
            .production_solution
            .iter()
            .all(|value| value.is_finite())
            && measurement
                .oracle_solution
                .iter()
                .all(|value| value.is_finite())
            && measurement
                .production_residual_vector
                .iter()
                .all(|value| value.is_finite())
            && measurement
                .oracle_residual_vector
                .iter()
                .all(|value| value.is_finite())
            && measurement.production_residual.is_finite()
            && measurement.oracle_residual.is_finite();
        if measurement.production_solution.len() != N
            || measurement.oracle_solution.len() != N
            || measurement.production_residual_vector.len() != N
            || measurement.oracle_residual_vector.len() != N
            || !finite
            || measurement.oracle_warning != "None"
            || measurement.oracle_backward_error.is_some()
            || measurement.oracle_certificate != "None"
        {
            return CaseOutcome::fail(format!(
                "stage=finite-shape-admission; fixture={}; production_len={}; oracle_len={}; production_residual_len={}; oracle_residual_len={}; finite={finite}; warning={}; backward_error={:?}; certificate={}",
                measurement.fixture.id,
                measurement.production_solution.len(),
                measurement.oracle_solution.len(),
                measurement.production_residual_vector.len(),
                measurement.oracle_residual_vector.len(),
                measurement.oracle_warning,
                measurement.oracle_backward_error,
                measurement.oracle_certificate,
            ))
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
        if measurement.production_residual > RESIDUAL_CEILING
            || measurement.oracle_residual > RESIDUAL_CEILING
        {
            return CaseOutcome::fail(format!(
                "stage=residual-admission; fixture={}; production_residual_bits={:016x?}; oracle_residual_bits={:016x?}; production_scaled_residual={}; oracle_scaled_residual={}; ceiling={RESIDUAL_CEILING}; input_frame_len={}; input_frame_fnv1a64=0x{input_frame_digest:016x}; input_frame={input_frame_hex}",
                measurement.fixture.id,
                bits(&measurement.production_residual_vector),
                bits(&measurement.oracle_residual_vector),
                measurement.production_residual,
                measurement.oracle_residual,
                input_frame.len(),
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#invariants")
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
        let maximum_delta = measurement
            .max_cross_delta()
            .max(measurement.max_production_reference_delta())
            .max(measurement.max_oracle_reference_delta());
        if !maximum_delta.is_finite() || maximum_delta > SOLUTION_CEILING {
            return CaseOutcome::fail(format!(
                "stage=solution-agreement; fixture={}; production_bits={:016x?}; oracle_bits={:016x?}; expected_bits={:016x?}; cross_delta={}; production_reference_delta={}; oracle_reference_delta={}; ceiling={SOLUTION_CEILING}; margin={}; input_frame_len={}; input_frame_fnv1a64=0x{input_frame_digest:016x}; input_frame={input_frame_hex}",
                measurement.fixture.id,
                bits(&measurement.production_solution),
                bits(&measurement.oracle_solution),
                bits(&measurement.fixture.expected),
                measurement.max_cross_delta(),
                measurement.max_production_reference_delta(),
                measurement.max_oracle_reference_delta(),
                measurement.remaining_margin(),
                input_frame.len(),
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#invariants")
            .with_evidence("constellation.lock:frankenscipy-0.1.0");
        }
    }

    let rows = first
        .iter()
        .map(|measurement| {
            format!(
                "{}:factor={},cross_delta={:.3e},prod_scaled_residual={:.3e},oracle_scaled_residual={:.3e},margin={:.3e}",
                measurement.fixture.id,
                measurement.fixture.kind.name(),
                measurement.max_cross_delta(),
                measurement.production_residual,
                measurement.oracle_residual,
                measurement.remaining_margin(),
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
    .with_evidence("crates/fs-la/CONTRACT.md#invariants")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn corruption_frame(index: usize, bit: u32, canonical: u64, corrupted: u64) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-la:exact-solution-reference-corruption:v1");
    push_field(&mut bytes, "canonical-exact-input-frame", &exact_inputs());
    push_u64_field(&mut bytes, "corruption-seed", CORRUPTION_SEED);
    push_text_field(&mut bytes, "reference", "exact-diagonal-lu-solution");
    push_usize_field(&mut bytes, "solution-index", index);
    push_u32_field(&mut bytes, "corrupted-bit", bit);
    push_u64_field(&mut bytes, "canonical-value-bits", canonical);
    push_u64_field(&mut bytes, "corrupted-value-bits", corrupted);
    bytes
}

fn reconstruct_corruption() -> Corruption {
    let index = usize::try_from((CORRUPTION_SEED >> 8) % N as u64)
        .expect("selected exact-solution index fits usize");
    let bit = u32::try_from((CORRUPTION_SEED & 0xff) % 52).expect("selected mantissa bit fits u32");
    let canonical = EXACT_SOLUTION[index].to_bits();
    let corrupted = canonical ^ (1_u64 << bit);
    Corruption {
        index,
        bit,
        canonical,
        corrupted,
        frame: corruption_frame(index, bit, canonical, corrupted),
    }
}

fn corruption_outcome(corruption: Corruption) -> CaseOutcome {
    let computed = match production_solve(FactorKind::Lu, &EXACT_MATRIX, &EXACT_RHS) {
        Ok(solution) => solution[corruption.index].to_bits(),
        Err(error) => return CaseOutcome::fail(format!("stage=production-solve; {error}")),
    };
    if computed == corruption.corrupted {
        CaseOutcome::pass("disclosed exact-solution reference corruption was not detected")
    } else {
        CaseOutcome::fail(format!(
            "seed=0x{CORRUPTION_SEED:016x}; reference=exact-diagonal-lu-solution; index={}; bit={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}; corrupted_bits=0x{:016x}; input_frame_len={}; input_frame_fnv1a64=0x{:016x}; input_frame={}",
            corruption.index,
            corruption.bit,
            corruption.canonical,
            corruption.corrupted,
            corruption.frame.len(),
            fnv1a64(&corruption.frame),
            hex_bytes(&corruption.frame),
        ))
        .with_evidence("crates/fs-la/tests/frankenscipy_linalg_oracle_casebook.rs#disclosed-corruption")
    }
}

fn run_red_report() -> SuiteReport {
    let corruption = reconstruct_corruption();
    assert_eq!(corruption.frame.len(), RED_FRAME_LEN);
    let digest = fnv1a64(&corruption.frame);
    assert_eq!(digest, RED_FRAME_FNV);
    Suite::new(SUITE)
        .case(
            "disclosed-exact-solution-reference-bit-corruption",
            digest,
            ToleranceSpec::Exact,
            move || corruption_outcome(corruption),
        )
        .run()
}

#[test]
fn bedrock_casebook_emits_exact_and_bounded_green_records() {
    let exact = exact_inputs();
    let oracle = oracle_inputs();
    assert_eq!(exact.len(), EXACT_FRAME_LEN);
    assert_eq!(fnv1a64(&exact), EXACT_FRAME_FNV);
    assert_eq!(oracle.len(), ORACLE_FRAME_LEN);
    assert_eq!(fnv1a64(&oracle), ORACLE_FRAME_FNV);

    let report = Suite::new(SUITE)
        .case(
            "exact-diagonal-factor-solve-kat",
            EXACT_FRAME_FNV,
            ToleranceSpec::Exact,
            exact_factor_outcome,
        )
        .case(
            "finite-dense-factorization-oracle-agreement",
            ORACLE_FRAME_FNV,
            ToleranceSpec::AbsoluteLe(SOLUTION_CEILING),
            oracle_outcome,
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
            "exact-diagonal-factor-solve-kat",
            "finite-dense-factorization-oracle-agreement",
        ]
    );
    assert!(report.records.iter().all(|record| {
        record.version == CASEBOOK_RECORD_VERSION && record.pass && !record.evidence.is_empty()
    }));
    assert_eq!(report.records[0].tolerance, "exact");
    assert_eq!(report.records[1].tolerance, "abs<=1e-12");
}

#[test]
fn disclosed_exact_solution_reference_corruption_is_stable_and_refused() {
    let first_corruption = reconstruct_corruption();
    let replay_corruption = reconstruct_corruption();
    assert_eq!(first_corruption.index, replay_corruption.index);
    assert_eq!(first_corruption.bit, replay_corruption.bit);
    assert_eq!(first_corruption.canonical, replay_corruption.canonical);
    assert_eq!(first_corruption.corrupted, replay_corruption.corrupted);
    assert_eq!(first_corruption.frame, replay_corruption.frame);
    assert_eq!(first_corruption.frame.len(), RED_FRAME_LEN);
    assert_eq!(fnv1a64(&first_corruption.frame), RED_FRAME_FNV);
    assert_eq!(first_corruption.index, 2);
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
        "disclosed-exact-solution-reference-bit-corruption"
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(
        first_failure
            .details
            .contains("reference=exact-diagonal-lu-solution")
    );
    assert!(first_failure.details.contains("index=2"));
    assert!(first_failure.details.contains("bit=7"));
    assert!(first_failure.details.contains("input_frame="));

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("assert_green must refuse the disclosed solution corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("disclosed-exact-solution-reference-bit-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
