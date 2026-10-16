//! Structured portable evidence for the sparse preconditioner surface.
//!
//! The exact case exercises ILU(0) as both a direct inverse application and
//! the preconditioner for a one-step PCG solve. Bounded cases retain a
//! Chebyshev solve on a disclosed diagonal spectrum and a genuine multilevel
//! smoothed-aggregation solve on a 9x9 grid Laplacian. A structural case
//! retains typed ILU breakdowns and the documented setup panics.
//!
//! Canonical input frames bind complete CSR storage, initial buffers, solver
//! options, versions, and every numerical ceiling. Output receipts bind every
//! direct-apply and solution bit, every `PcgReport` field, signed residuals,
//! Chebyshev bands, AMG hierarchy and complexity, and refusal observations.
//!
//! This is finite-fixture G0 and same-build replay evidence. It is not an
//! external oracle, a performance or cancellation claim, a large-grid or
//! grid-independence study, proof for arbitrary SPD matrices, or fresh
//! cross-ISA G5 evidence. In particular, no AMG output bit is promoted to a
//! cross-ISA golden here.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_casebook::{
    CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, SuiteReport, ToleranceSpec, fnv1a64,
};
use fs_sparse::precond::{Chebyshev, PcgReport, Precond, SaAmg, ilu0, pcg};
use fs_sparse::{Csr, VERSION as FS_SPARSE_VERSION};

const SUITE: &str = "bedrock/fs-sparse-preconditioner-casebook-v1";
const FRAME_ENCODING: &str =
    "field=(tag_len:u64le,tag,payload_len:u64le,payload);numbers=le;f64=bits:v1";
const RESIDUAL_POLICY: &str =
    "signed-r=Ax-b;scaled-inf=maxabs(r)/(row-sum-inf*maxabs(x)+maxabs(b)):v1";

const EXACT_ROW_PTR: [usize; 4] = [0, 2, 5, 7];
const EXACT_COLUMNS: [usize; 7] = [0, 1, 0, 1, 2, 1, 2];
const EXACT_VALUES: [f64; 7] = [4.0, 2.0, 2.0, 5.0, 2.0, 2.0, 5.0];
const EXACT_RHS: [f64; 3] = [8.0, 18.0, 19.0];
const EXACT_SOLUTION: [f64; 3] = [1.0, 2.0, 3.0];
const EXACT_PCG_TOL: f64 = 0.0;
const EXACT_PCG_CAP: usize = 1;

const CHEB_DEGREE: usize = 4;
const CHEB_ALPHA: f64 = 30.0;
const CHEB_PCG_TOL: f64 = 1.0e-12;
const CHEB_PCG_CAP: usize = 12;
const CHEB_SOLUTION_CEILING: f64 = 1.0e-9;
const CHEB_RESIDUAL_CEILING: f64 = 1.0e-12;
const CHEB_BAND_HI_MIN: f64 = 4.0;
const CHEB_BAND_HI_MAX: f64 = 4.5;

const AMG_GRID: usize = 9;
const AMG_THETA: f64 = 0.08;
const AMG_SMOOTHER_DEGREE: usize = 3;
const AMG_PCG_TOL: f64 = 1.0e-10;
const AMG_PCG_CAP: usize = 100;
const AMG_SOLUTION_CEILING: f64 = 1.0e-7;
const AMG_RESIDUAL_CEILING: f64 = 1.0e-10;
const AMG_COMPLEXITY_MIN: f64 = 1.0;
const AMG_COMPLEXITY_MAX: f64 = 2.0;
const AMG_MIN_LEVELS: usize = 2;

const CORRUPTION_SEED: u64 = 0xF5A5_0021_0000_0101;
const EXACT_FRAME_LEN: usize = 1_072;
const EXACT_FRAME_FNV1A64: u64 = 0x4aed_c303_b42f_00cf;
const CHEBYSHEV_FRAME_LEN: usize = 1_191;
const CHEBYSHEV_FRAME_FNV1A64: u64 = 0x3d24_759f_a982_a6fd;
const AMG_FRAME_LEN: usize = 10_428;
const AMG_FRAME_FNV1A64: u64 = 0x66d6_7fc3_3896_1288;
const REFUSAL_FRAME_LEN: usize = 1_874;
const REFUSAL_FRAME_FNV1A64: u64 = 0x9ece_4960_5afd_e61d;
const CORRUPTION_FRAME_LEN: usize = 1_721;
const CORRUPTION_FRAME_FNV1A64: u64 = 0x9dc6_64b9_bc83_5f75;

#[derive(Debug, Clone)]
struct CsrSpec {
    nrows: usize,
    ncols: usize,
    row_ptr: Vec<usize>,
    columns: Vec<usize>,
    values: Vec<f64>,
}

impl CsrSpec {
    fn build(&self) -> Csr {
        Csr::from_parts(
            self.nrows,
            self.ncols,
            self.row_ptr.clone(),
            self.columns.clone(),
            self.values.clone(),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PcgBits {
    iters: usize,
    rel_residual: u64,
    converged: bool,
}

impl From<&PcgReport> for PcgBits {
    fn from(report: &PcgReport) -> Self {
        Self {
            iters: report.iters,
            rel_residual: report.rel_residual.to_bits(),
            converged: report.converged,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct NumericMeasurement {
    apply: Vec<u64>,
    solution: Vec<u64>,
    report: PcgBits,
    residual: Vec<u64>,
    residual_scale: u64,
    scaled_residual: u64,
    solution_delta: Vec<u64>,
    max_solution_delta: u64,
    band: Option<(u64, u64)>,
    hierarchy: Vec<usize>,
    operator_complexity: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PanicObservation {
    panicked: bool,
    message: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RefusalMeasurement {
    missing_diagonal_row: Option<usize>,
    stored_zero_row: Option<usize>,
    chebyshev_degree: PanicObservation,
    chebyshev_alpha: PanicObservation,
    nonsquare_ilu: PanicObservation,
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
        u64::try_from(value).expect("preconditioner Casebook lengths fit u64"),
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
        u64::try_from(value).expect("preconditioner Casebook values fit u64"),
    );
}

fn push_bool_field(bytes: &mut Vec<u8>, tag: &str, value: bool) {
    push_u32_field(bytes, tag, u32::from(value));
}

fn push_f64_field(bytes: &mut Vec<u8>, tag: &str, value: f64) {
    push_u64_field(bytes, tag, value.to_bits());
}

fn push_usizes_field(bytes: &mut Vec<u8>, tag: &str, values: &[usize]) {
    let mut payload = Vec::with_capacity(8 + values.len() * 8);
    push_len(&mut payload, values.len());
    for &value in values {
        push_len(&mut payload, value);
    }
    push_field(bytes, tag, &payload);
}

fn push_f64s_field(bytes: &mut Vec<u8>, tag: &str, values: &[f64]) {
    let mut payload = Vec::with_capacity(8 + values.len() * 8);
    push_len(&mut payload, values.len());
    for value in values {
        push_u64(&mut payload, value.to_bits());
    }
    push_field(bytes, tag, &payload);
}

fn push_bits_field(bytes: &mut Vec<u8>, tag: &str, values: &[u64]) {
    let mut payload = Vec::with_capacity(8 + values.len() * 8);
    push_len(&mut payload, values.len());
    for &value in values {
        push_u64(&mut payload, value);
    }
    push_field(bytes, tag, &payload);
}

fn push_csr(bytes: &mut Vec<u8>, tag: &str, spec: &CsrSpec) {
    let mut payload = Vec::new();
    push_usize_field(&mut payload, "nrows", spec.nrows);
    push_usize_field(&mut payload, "ncols", spec.ncols);
    push_usizes_field(&mut payload, "row-pointers", &spec.row_ptr);
    push_usizes_field(&mut payload, "column-indices", &spec.columns);
    push_f64s_field(&mut payload, "stored-values", &spec.values);
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
    push_text_field(&mut bytes, "fs-sparse-version", FS_SPARSE_VERSION);
    push_text_field(&mut bytes, "residual-policy", RESIDUAL_POLICY);
    bytes
}

fn exact_spec() -> CsrSpec {
    CsrSpec {
        nrows: 3,
        ncols: 3,
        row_ptr: EXACT_ROW_PTR.to_vec(),
        columns: EXACT_COLUMNS.to_vec(),
        values: EXACT_VALUES.to_vec(),
    }
}

fn chebyshev_spec() -> CsrSpec {
    CsrSpec {
        nrows: 3,
        ncols: 3,
        row_ptr: vec![0, 1, 2, 3],
        columns: vec![0, 1, 2],
        values: vec![1.0, 2.0, 4.0],
    }
}

fn chebyshev_solution() -> Vec<f64> {
    vec![1.0, -2.0, 0.5]
}

fn chebyshev_rhs() -> Vec<f64> {
    vec![1.0, -4.0, 2.0]
}

fn laplacian_2d_spec(grid: usize) -> CsrSpec {
    let n = grid * grid;
    let mut row_ptr = Vec::with_capacity(n + 1);
    let mut columns = Vec::with_capacity(5 * n);
    let mut values = Vec::with_capacity(5 * n);
    row_ptr.push(0);
    for row in 0..grid {
        for column in 0..grid {
            let node = row * grid + column;
            let mut entries = Vec::with_capacity(5);
            entries.push((node, 4.0));
            if row > 0 {
                entries.push((node - grid, -1.0));
            }
            if row + 1 < grid {
                entries.push((node + grid, -1.0));
            }
            if column > 0 {
                entries.push((node - 1, -1.0));
            }
            if column + 1 < grid {
                entries.push((node + 1, -1.0));
            }
            entries.sort_unstable_by_key(|entry| entry.0);
            for (index, value) in entries {
                columns.push(index);
                values.push(value);
            }
            row_ptr.push(columns.len());
        }
    }
    CsrSpec {
        nrows: n,
        ncols: n,
        row_ptr,
        columns,
        values,
    }
}

fn amg_solution() -> Vec<f64> {
    (0..AMG_GRID * AMG_GRID)
        .map(|index| {
            let numerator =
                i32::try_from((index * 5 + 3) % 17).expect("AMG fixture numerator fits i32") - 8;
            f64::from(numerator) * 0.125
        })
        .collect()
}

fn rhs_from_spec(spec: &CsrSpec, solution: &[f64]) -> Vec<f64> {
    (0..spec.nrows)
        .map(|row| {
            let mut value = 0.0_f64;
            for position in spec.row_ptr[row]..spec.row_ptr[row + 1] {
                value = spec.values[position].mul_add(solution[spec.columns[position]], value);
            }
            value
        })
        .collect()
}

fn exact_inputs() -> Vec<u8> {
    let spec = exact_spec();
    let mut bytes = common_frame("bedrock:fs-sparse:exact-ilu-pcg-kat:v1");
    push_csr(&mut bytes, "matrix", &spec);
    push_f64s_field(&mut bytes, "rhs", &EXACT_RHS);
    push_f64s_field(&mut bytes, "known-solution", &EXACT_SOLUTION);
    push_f64s_field(&mut bytes, "apply-initial-output", &[0.0; 3]);
    push_f64s_field(&mut bytes, "pcg-initial-iterate", &[0.0; 3]);
    push_f64_field(&mut bytes, "pcg-tolerance", EXACT_PCG_TOL);
    push_usize_field(&mut bytes, "pcg-iteration-cap", EXACT_PCG_CAP);
    push_text_field(
        &mut bytes,
        "claim",
        "ILU(0)-apply=known-solution;PCG=one-step-exact:v1",
    );
    bytes
}

fn chebyshev_inputs() -> Vec<u8> {
    let spec = chebyshev_spec();
    let rhs = chebyshev_rhs();
    let solution = chebyshev_solution();
    let mut bytes = common_frame("bedrock:fs-sparse:chebyshev-diagonal-replay:v1");
    push_csr(&mut bytes, "matrix", &spec);
    push_f64s_field(&mut bytes, "rhs", &rhs);
    push_f64s_field(&mut bytes, "known-solution", &solution);
    push_f64s_field(&mut bytes, "apply-initial-output", &[0.0; 3]);
    push_f64s_field(&mut bytes, "pcg-initial-iterate", &[0.0; 3]);
    push_usize_field(&mut bytes, "chebyshev-degree", CHEB_DEGREE);
    push_f64_field(&mut bytes, "chebyshev-alpha", CHEB_ALPHA);
    push_f64_field(&mut bytes, "band-hi-min", CHEB_BAND_HI_MIN);
    push_f64_field(&mut bytes, "band-hi-max", CHEB_BAND_HI_MAX);
    push_f64_field(&mut bytes, "pcg-tolerance", CHEB_PCG_TOL);
    push_usize_field(&mut bytes, "pcg-iteration-cap", CHEB_PCG_CAP);
    push_f64_field(
        &mut bytes,
        "solution-absolute-ceiling",
        CHEB_SOLUTION_CEILING,
    );
    push_f64_field(&mut bytes, "scaled-residual-ceiling", CHEB_RESIDUAL_CEILING);
    bytes
}

fn amg_inputs() -> Vec<u8> {
    let spec = laplacian_2d_spec(AMG_GRID);
    let solution = amg_solution();
    let rhs = rhs_from_spec(&spec, &solution);
    let mut bytes = common_frame("bedrock:fs-sparse:sa-amg-9x9-grid-replay:v1");
    push_usize_field(&mut bytes, "grid-rows", AMG_GRID);
    push_usize_field(&mut bytes, "grid-columns", AMG_GRID);
    push_csr(&mut bytes, "matrix", &spec);
    push_f64s_field(&mut bytes, "rhs", &rhs);
    push_f64s_field(&mut bytes, "known-dyadic-solution", &solution);
    push_f64s_field(&mut bytes, "apply-initial-output", &vec![0.0; spec.nrows]);
    push_f64s_field(&mut bytes, "pcg-initial-iterate", &vec![0.0; spec.nrows]);
    push_f64_field(&mut bytes, "amg-strength-theta", AMG_THETA);
    push_usize_field(
        &mut bytes,
        "amg-chebyshev-smoother-degree",
        AMG_SMOOTHER_DEGREE,
    );
    push_f64_field(&mut bytes, "pcg-tolerance", AMG_PCG_TOL);
    push_usize_field(&mut bytes, "pcg-iteration-cap", AMG_PCG_CAP);
    push_usize_field(&mut bytes, "minimum-hierarchy-levels", AMG_MIN_LEVELS);
    push_f64_field(&mut bytes, "operator-complexity-min", AMG_COMPLEXITY_MIN);
    push_f64_field(&mut bytes, "operator-complexity-max", AMG_COMPLEXITY_MAX);
    push_f64_field(
        &mut bytes,
        "solution-absolute-ceiling",
        AMG_SOLUTION_CEILING,
    );
    push_f64_field(&mut bytes, "scaled-residual-ceiling", AMG_RESIDUAL_CEILING);
    push_text_field(
        &mut bytes,
        "golden-policy",
        "same-build-byte-replay-only;do-not-promote-amg-output-bits-cross-isa:v1",
    );
    bytes
}

fn missing_diagonal_spec() -> CsrSpec {
    CsrSpec {
        nrows: 2,
        ncols: 2,
        row_ptr: vec![0, 2, 3],
        columns: vec![0, 1, 0],
        values: vec![1.0, 1.0, 1.0],
    }
}

fn stored_zero_spec() -> CsrSpec {
    CsrSpec {
        nrows: 2,
        ncols: 2,
        row_ptr: vec![0, 2, 4],
        columns: vec![0, 1, 0, 1],
        values: vec![0.0, 1.0, 1.0, 2.0],
    }
}

fn nonsquare_spec() -> CsrSpec {
    CsrSpec {
        nrows: 2,
        ncols: 3,
        row_ptr: vec![0, 1, 2],
        columns: vec![0, 1],
        values: vec![1.0, 1.0],
    }
}

fn refusal_inputs() -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-sparse:preconditioner-refusals:v1");
    push_csr(
        &mut bytes,
        "ilu-missing-diagonal-matrix",
        &missing_diagonal_spec(),
    );
    push_usize_field(&mut bytes, "expected-missing-diagonal-row", 1);
    push_csr(&mut bytes, "ilu-stored-zero-matrix", &stored_zero_spec());
    push_usize_field(&mut bytes, "expected-stored-zero-row", 0);
    push_csr(&mut bytes, "chebyshev-validation-matrix", &chebyshev_spec());
    push_usize_field(&mut bytes, "rejected-chebyshev-degree", 0);
    push_f64_field(&mut bytes, "rejected-chebyshev-alpha", 1.0);
    push_csr(&mut bytes, "nonsquare-ilu-matrix", &nonsquare_spec());
    push_text_field(
        &mut bytes,
        "required-degree-panic",
        "Chebyshev degree must be >= 1",
    );
    push_text_field(
        &mut bytes,
        "required-alpha-panic",
        "band divisor must exceed 1",
    );
    push_text_field(
        &mut bytes,
        "required-nonsquare-panic",
        "ilu0 requires a square matrix",
    );
    bytes
}

fn output_bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn finite_bits(values: &[u64]) -> bool {
    values
        .iter()
        .all(|&value| f64::from_bits(value).is_finite())
}

fn residual_evidence(spec: &CsrSpec, solution: &[f64], rhs: &[f64]) -> (Vec<u64>, u64, u64) {
    let mut residual = Vec::with_capacity(spec.nrows);
    for (row, bounds) in spec.row_ptr.windows(2).enumerate() {
        let mut value = 0.0_f64;
        for position in bounds[0]..bounds[1] {
            value = spec.values[position].mul_add(solution[spec.columns[position]], value);
        }
        residual.push(value - rhs[row]);
    }
    let residual_max = residual.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let matrix_norm = (0..spec.nrows)
        .map(|row| {
            spec.values[spec.row_ptr[row]..spec.row_ptr[row + 1]]
                .iter()
                .map(|value| value.abs())
                .sum::<f64>()
        })
        .fold(0.0, f64::max);
    let solution_norm = solution.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let rhs_norm = rhs.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let scale = matrix_norm.mul_add(solution_norm, rhs_norm);
    let scaled = if scale == 0.0 {
        residual_max
    } else {
        residual_max / scale
    };
    (output_bits(&residual), scale.to_bits(), scaled.to_bits())
}

fn solution_delta(actual: &[f64], expected: &[f64]) -> (Vec<u64>, u64) {
    let delta = actual
        .iter()
        .zip(expected)
        .map(|(actual, expected)| actual - expected)
        .collect::<Vec<_>>();
    let maximum = delta.iter().map(|value| value.abs()).fold(0.0, f64::max);
    (output_bits(&delta), maximum.to_bits())
}

/// Per-case acceptance settings for `measure`, named so call sites read as
/// declarations instead of a positional None/Vec::new() parade.
struct MeasureSettings {
    tolerance: f64,
    iteration_cap: usize,
    band: Option<(f64, f64)>,
    hierarchy: Vec<usize>,
    operator_complexity: Option<f64>,
}

fn measure<P: Precond>(
    spec: &CsrSpec,
    rhs: &[f64],
    expected: &[f64],
    preconditioner: &P,
    settings: MeasureSettings,
) -> NumericMeasurement {
    let MeasureSettings {
        tolerance,
        iteration_cap,
        band,
        hierarchy,
        operator_complexity,
    } = settings;
    let matrix = spec.build();
    let mut apply = vec![0.0; rhs.len()];
    preconditioner.apply(rhs, &mut apply);
    let mut solution = vec![0.0; rhs.len()];
    let report = pcg(
        &matrix,
        rhs,
        &mut solution,
        preconditioner,
        tolerance,
        iteration_cap,
    );
    let (residual, residual_scale, scaled_residual) = residual_evidence(spec, &solution, rhs);
    let (solution_delta, max_solution_delta) = solution_delta(&solution, expected);
    NumericMeasurement {
        apply: output_bits(&apply),
        solution: output_bits(&solution),
        report: PcgBits::from(&report),
        residual,
        residual_scale,
        scaled_residual,
        solution_delta,
        max_solution_delta,
        band: band.map(|(lo, hi)| (lo.to_bits(), hi.to_bits())),
        hierarchy,
        operator_complexity: operator_complexity.map(f64::to_bits),
    }
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| {
            payload
                .downcast_ref::<&str>()
                .map(|message| (*message).to_owned())
        })
        .unwrap_or_else(|| "non-text panic payload".to_owned())
}

fn capture<T>(stage: &str, run: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
    match catch_unwind(AssertUnwindSafe(run)) {
        Ok(result) => result,
        Err(payload) => Err(format!("stage={stage}; panic={}", panic_message(&*payload))),
    }
}

fn observe_panic(run: impl FnOnce()) -> PanicObservation {
    match catch_unwind(AssertUnwindSafe(run)) {
        Ok(()) => PanicObservation {
            panicked: false,
            message: "returned-without-panic".to_owned(),
        },
        Err(payload) => PanicObservation {
            panicked: true,
            message: panic_message(&*payload),
        },
    }
}

fn exact_measurement() -> Result<NumericMeasurement, String> {
    let spec = exact_spec();
    let matrix = spec.build();
    let preconditioner = ilu0(&matrix)
        .map_err(|error| format!("stage=exact-ilu-setup; breakdown-row={}", error.row))?;
    Ok(measure(
        &spec,
        &EXACT_RHS,
        &EXACT_SOLUTION,
        &preconditioner,
        MeasureSettings {
            tolerance: EXACT_PCG_TOL,
            iteration_cap: EXACT_PCG_CAP,
            band: None,
            hierarchy: Vec::new(),
            operator_complexity: None,
        },
    ))
}

#[expect(
    clippy::unnecessary_wraps,
    reason = "returns Result to keep the panic-capturing `capture` harness \\
              signature uniform with `exact_measurement`, which genuinely \\
              fails on ILU breakdown"
)]
fn chebyshev_measurement() -> Result<NumericMeasurement, String> {
    let spec = chebyshev_spec();
    let matrix = spec.build();
    let preconditioner = Chebyshev::new(&matrix, CHEB_DEGREE, CHEB_ALPHA);
    let band = preconditioner.band();
    Ok(measure(
        &spec,
        &chebyshev_rhs(),
        &chebyshev_solution(),
        &preconditioner,
        MeasureSettings {
            tolerance: CHEB_PCG_TOL,
            iteration_cap: CHEB_PCG_CAP,
            band: Some(band),
            hierarchy: Vec::new(),
            operator_complexity: None,
        },
    ))
}

#[expect(
    clippy::unnecessary_wraps,
    reason = "returns Result to keep the panic-capturing `capture` harness \\
              signature uniform with `exact_measurement`, which genuinely \\
              fails on ILU breakdown"
)]
fn amg_measurement() -> Result<NumericMeasurement, String> {
    let spec = laplacian_2d_spec(AMG_GRID);
    let solution = amg_solution();
    let rhs = rhs_from_spec(&spec, &solution);
    let matrix = spec.build();
    let preconditioner = SaAmg::new(&matrix, AMG_THETA, AMG_SMOOTHER_DEGREE);
    let hierarchy = preconditioner.level_sizes.clone();
    let complexity = preconditioner.operator_complexity();
    Ok(measure(
        &spec,
        &rhs,
        &solution,
        &preconditioner,
        MeasureSettings {
            tolerance: AMG_PCG_TOL,
            iteration_cap: AMG_PCG_CAP,
            band: None,
            hierarchy,
            operator_complexity: Some(complexity),
        },
    ))
}

fn refusal_measurement() -> RefusalMeasurement {
    let missing_diagonal_row = ilu0(&missing_diagonal_spec().build())
        .err()
        .map(|error| error.row);
    let stored_zero_row = ilu0(&stored_zero_spec().build())
        .err()
        .map(|error| error.row);
    let chebyshev_degree = observe_panic(|| {
        let _ = Chebyshev::new(&chebyshev_spec().build(), 0, CHEB_ALPHA);
    });
    let chebyshev_alpha = observe_panic(|| {
        let _ = Chebyshev::new(&chebyshev_spec().build(), CHEB_DEGREE, 1.0);
    });
    let nonsquare_ilu = observe_panic(|| {
        let _ = ilu0(&nonsquare_spec().build());
    });
    RefusalMeasurement {
        missing_diagonal_row,
        stored_zero_row,
        chebyshev_degree,
        chebyshev_alpha,
        nonsquare_ilu,
    }
}

fn push_pcg_report(bytes: &mut Vec<u8>, report: &PcgBits) {
    push_usize_field(bytes, "pcg-iters", report.iters);
    push_u64_field(bytes, "pcg-relative-residual-bits", report.rel_residual);
    push_bool_field(bytes, "pcg-converged", report.converged);
}

fn numeric_receipt(domain: &str, inputs: &[u8], measurement: &NumericMeasurement) -> Vec<u8> {
    let mut bytes = common_frame(domain);
    push_field(&mut bytes, "canonical-input-frame", inputs);
    push_bits_field(&mut bytes, "preconditioner-apply-bits", &measurement.apply);
    push_bits_field(&mut bytes, "pcg-solution-bits", &measurement.solution);
    push_pcg_report(&mut bytes, &measurement.report);
    push_bits_field(&mut bytes, "signed-residual-bits", &measurement.residual);
    push_u64_field(
        &mut bytes,
        "residual-scale-bits",
        measurement.residual_scale,
    );
    push_u64_field(
        &mut bytes,
        "scaled-residual-bits",
        measurement.scaled_residual,
    );
    push_bits_field(
        &mut bytes,
        "solution-delta-bits",
        &measurement.solution_delta,
    );
    push_u64_field(
        &mut bytes,
        "max-solution-delta-bits",
        measurement.max_solution_delta,
    );
    push_bool_field(
        &mut bytes,
        "chebyshev-band-present",
        measurement.band.is_some(),
    );
    if let Some((lo, hi)) = measurement.band {
        push_u64_field(&mut bytes, "chebyshev-band-lo-bits", lo);
        push_u64_field(&mut bytes, "chebyshev-band-hi-bits", hi);
    }
    push_usizes_field(&mut bytes, "amg-level-sizes", &measurement.hierarchy);
    push_bool_field(
        &mut bytes,
        "operator-complexity-present",
        measurement.operator_complexity.is_some(),
    );
    if let Some(complexity) = measurement.operator_complexity {
        push_u64_field(&mut bytes, "operator-complexity-bits", complexity);
    }
    bytes
}

fn refusal_receipt(inputs: &[u8], measurement: &RefusalMeasurement) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-sparse:preconditioner-refusal-receipt:v1");
    push_field(&mut bytes, "canonical-input-frame", inputs);
    push_bool_field(
        &mut bytes,
        "missing-diagonal-row-present",
        measurement.missing_diagonal_row.is_some(),
    );
    if let Some(row) = measurement.missing_diagonal_row {
        push_usize_field(&mut bytes, "missing-diagonal-row", row);
    }
    push_bool_field(
        &mut bytes,
        "stored-zero-row-present",
        measurement.stored_zero_row.is_some(),
    );
    if let Some(row) = measurement.stored_zero_row {
        push_usize_field(&mut bytes, "stored-zero-row", row);
    }
    for (prefix, observation) in [
        ("chebyshev-degree", &measurement.chebyshev_degree),
        ("chebyshev-alpha", &measurement.chebyshev_alpha),
        ("nonsquare-ilu", &measurement.nonsquare_ilu),
    ] {
        push_bool_field(
            &mut bytes,
            &format!("{prefix}-panicked"),
            observation.panicked,
        );
        push_text_field(
            &mut bytes,
            &format!("{prefix}-panic-message"),
            &observation.message,
        );
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

fn replay_failure(stage: &str, first: &[u8], replay: &[u8]) -> CaseOutcome {
    CaseOutcome::fail(format!(
        "stage={stage}; first_len={}; replay_len={}; first_fnv1a64=0x{:016x}; replay_fnv1a64=0x{:016x}; first={}; replay={}",
        first.len(),
        replay.len(),
        fnv1a64(first),
        fnv1a64(replay),
        hex_bytes(first),
        hex_bytes(replay),
    ))
    .with_evidence("crates/fs-sparse/CONTRACT.md#determinism-class")
}

fn receipt_pass(details: &str, receipt: &[u8]) -> CaseOutcome {
    CaseOutcome::pass(format!(
        "{details}; output_receipt_len={}; output_receipt_fnv1a64=0x{:016x}; output_receipt={}",
        receipt.len(),
        fnv1a64(receipt),
        hex_bytes(receipt),
    ))
    .with_evidence("crates/fs-sparse/CONTRACT.md#invariants")
    .with_evidence("crates/fs-sparse/tests/preconditioner_casebook.rs")
}

fn numeric_shape_and_finite(
    measurement: &NumericMeasurement,
    dimension: usize,
) -> Result<(), String> {
    if measurement.apply.len() != dimension
        || measurement.solution.len() != dimension
        || measurement.residual.len() != dimension
        || measurement.solution_delta.len() != dimension
    {
        return Err(format!(
            "stage=shape; dimension={dimension}; apply={}; solution={}; residual={}; delta={}",
            measurement.apply.len(),
            measurement.solution.len(),
            measurement.residual.len(),
            measurement.solution_delta.len(),
        ));
    }
    let scalar_bits = [
        measurement.report.rel_residual,
        measurement.residual_scale,
        measurement.scaled_residual,
        measurement.max_solution_delta,
    ];
    if !finite_bits(&measurement.apply)
        || !finite_bits(&measurement.solution)
        || !finite_bits(&measurement.residual)
        || !finite_bits(&measurement.solution_delta)
        || !finite_bits(&scalar_bits)
        || measurement
            .band
            .is_some_and(|(lo, hi)| !finite_bits(&[lo, hi]))
        || measurement
            .operator_complexity
            .is_some_and(|value| !f64::from_bits(value).is_finite())
    {
        return Err("stage=finite-admission; non-finite retained output".to_owned());
    }
    Ok(())
}

fn exact_outcome() -> CaseOutcome {
    let inputs = exact_inputs();
    let first = match capture("exact-first", exact_measurement) {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay = match capture("exact-replay", exact_measurement) {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let first_receipt =
        numeric_receipt("bedrock:fs-sparse:exact-ilu-pcg-output:v1", &inputs, &first);
    let replay_receipt = numeric_receipt(
        "bedrock:fs-sparse:exact-ilu-pcg-output:v1",
        &inputs,
        &replay,
    );
    if first_receipt != replay_receipt {
        return replay_failure("exact-same-run-replay", &first_receipt, &replay_receipt);
    }
    if let Err(error) = numeric_shape_and_finite(&first, 3) {
        return CaseOutcome::fail(error);
    }
    let expected = output_bits(&EXACT_SOLUTION);
    if first.apply != expected
        || first.solution != expected
        || first.residual != vec![0; 3]
        || first.solution_delta != vec![0; 3]
        || first.max_solution_delta != 0.0_f64.to_bits()
        || first.report.iters != 1
        || first.report.rel_residual != 0.0_f64.to_bits()
        || !first.report.converged
    {
        return CaseOutcome::fail(format!(
            "stage=exact-ilu-pcg-kat; apply={:016x?}; solution={:016x?}; expected={expected:016x?}; residual={:016x?}; delta={:016x?}; report={:?}",
            first.apply, first.solution, first.residual, first.solution_delta, first.report,
        ))
        .with_evidence("crates/fs-sparse/CONTRACT.md#invariants");
    }
    receipt_pass(
        "apply=exact; pcg=one-step-exact; same_run=byte-identical",
        &first_receipt,
    )
}

fn chebyshev_outcome() -> CaseOutcome {
    let inputs = chebyshev_inputs();
    let first = match capture("chebyshev-first", chebyshev_measurement) {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay = match capture("chebyshev-replay", chebyshev_measurement) {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let first_receipt = numeric_receipt(
        "bedrock:fs-sparse:chebyshev-diagonal-output:v1",
        &inputs,
        &first,
    );
    let replay_receipt = numeric_receipt(
        "bedrock:fs-sparse:chebyshev-diagonal-output:v1",
        &inputs,
        &replay,
    );
    if first_receipt != replay_receipt {
        return replay_failure("chebyshev-same-run-replay", &first_receipt, &replay_receipt);
    }
    if let Err(error) = numeric_shape_and_finite(&first, 3) {
        return CaseOutcome::fail(error);
    }
    let Some((lo_bits, hi_bits)) = first.band else {
        return CaseOutcome::fail("stage=chebyshev-band; missing band evidence");
    };
    let lo = f64::from_bits(lo_bits);
    let hi = f64::from_bits(hi_bits);
    let solution_delta = f64::from_bits(first.max_solution_delta);
    let residual = f64::from_bits(first.scaled_residual);
    let reported = f64::from_bits(first.report.rel_residual);
    if !(lo > 0.0
        && lo < hi
        && lo.to_bits() == (hi / CHEB_ALPHA).to_bits()
        && (CHEB_BAND_HI_MIN..=CHEB_BAND_HI_MAX).contains(&hi))
    {
        return CaseOutcome::fail(format!(
            "stage=chebyshev-band; lo={lo}; hi={hi}; expected_lo_bits=0x{:016x}; hi_range=[{CHEB_BAND_HI_MIN},{CHEB_BAND_HI_MAX}]",
            (hi / CHEB_ALPHA).to_bits(),
        ));
    }
    if !first.report.converged
        || first.report.iters > CHEB_PCG_CAP
        || reported > CHEB_PCG_TOL
        || solution_delta > CHEB_SOLUTION_CEILING
        || residual > CHEB_RESIDUAL_CEILING
    {
        return CaseOutcome::fail(format!(
            "stage=chebyshev-solve-gates; converged={}; iters={}; reported={reported}; solution_delta={solution_delta}; scaled_residual={residual}; tol={CHEB_PCG_TOL}; solution_ceiling={CHEB_SOLUTION_CEILING}; residual_ceiling={CHEB_RESIDUAL_CEILING}",
            first.report.converged, first.report.iters,
        ));
    }
    receipt_pass(
        &format!(
            "degree={CHEB_DEGREE}; alpha={CHEB_ALPHA}; band=[{lo},{hi}]; iters={}; solution_delta={solution_delta:.3e}; scaled_residual={residual:.3e}; same_run=byte-identical",
            first.report.iters,
        ),
        &first_receipt,
    )
}

fn amg_outcome() -> CaseOutcome {
    let inputs = amg_inputs();
    let first = match capture("amg-first", amg_measurement) {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay = match capture("amg-replay", amg_measurement) {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let first_receipt = numeric_receipt(
        "bedrock:fs-sparse:sa-amg-9x9-grid-output:v1",
        &inputs,
        &first,
    );
    let replay_receipt = numeric_receipt(
        "bedrock:fs-sparse:sa-amg-9x9-grid-output:v1",
        &inputs,
        &replay,
    );
    if first_receipt != replay_receipt {
        return replay_failure("amg-same-run-replay", &first_receipt, &replay_receipt);
    }
    let dimension = AMG_GRID * AMG_GRID;
    if let Err(error) = numeric_shape_and_finite(&first, dimension) {
        return CaseOutcome::fail(error);
    }
    let hierarchy_is_genuine = first.hierarchy.len() >= AMG_MIN_LEVELS
        && first.hierarchy.first() == Some(&dimension)
        && first.hierarchy.windows(2).all(|pair| pair[1] < pair[0]);
    let Some(complexity_bits) = first.operator_complexity else {
        return CaseOutcome::fail("stage=amg-complexity; missing operator complexity");
    };
    let complexity = f64::from_bits(complexity_bits);
    let solution_delta = f64::from_bits(first.max_solution_delta);
    let residual = f64::from_bits(first.scaled_residual);
    let reported = f64::from_bits(first.report.rel_residual);
    if !hierarchy_is_genuine || complexity < AMG_COMPLEXITY_MIN || complexity >= AMG_COMPLEXITY_MAX
    {
        return CaseOutcome::fail(format!(
            "stage=amg-hierarchy-complexity; hierarchy={:?}; genuine={hierarchy_is_genuine}; complexity={complexity}; range=[{AMG_COMPLEXITY_MIN},{AMG_COMPLEXITY_MAX})",
            first.hierarchy,
        ));
    }
    if !first.report.converged
        || first.report.iters > AMG_PCG_CAP
        || reported > AMG_PCG_TOL
        || solution_delta > AMG_SOLUTION_CEILING
        || residual > AMG_RESIDUAL_CEILING
    {
        return CaseOutcome::fail(format!(
            "stage=amg-solve-gates; converged={}; iters={}; reported={reported}; solution_delta={solution_delta}; scaled_residual={residual}; tol={AMG_PCG_TOL}; solution_ceiling={AMG_SOLUTION_CEILING}; residual_ceiling={AMG_RESIDUAL_CEILING}",
            first.report.converged, first.report.iters,
        ));
    }
    receipt_pass(
        &format!(
            "grid={AMG_GRID}x{AMG_GRID}; levels={:?}; complexity={complexity:.6}; iters={}; solution_delta={solution_delta:.3e}; scaled_residual={residual:.3e}; same_run=byte-identical; cross_isa_golden=no-claim",
            first.hierarchy, first.report.iters,
        ),
        &first_receipt,
    )
}

fn refusal_outcome() -> CaseOutcome {
    let inputs = refusal_inputs();
    let first = refusal_measurement();
    let replay = refusal_measurement();
    let first_receipt = refusal_receipt(&inputs, &first);
    let replay_receipt = refusal_receipt(&inputs, &replay);
    if first_receipt != replay_receipt {
        return replay_failure("refusal-same-run-replay", &first_receipt, &replay_receipt);
    }
    let admitted = first.missing_diagonal_row == Some(1)
        && first.stored_zero_row == Some(0)
        && first.chebyshev_degree.panicked
        && first
            .chebyshev_degree
            .message
            .contains("Chebyshev degree must be >= 1")
        && first.chebyshev_alpha.panicked
        && first
            .chebyshev_alpha
            .message
            .contains("band divisor must exceed 1")
        && first.nonsquare_ilu.panicked
        && first
            .nonsquare_ilu
            .message
            .contains("ilu0 requires a square matrix");
    if !admitted {
        return CaseOutcome::fail(format!(
            "stage=refusal-admission; measurement={first:?}; output_receipt={}",
            hex_bytes(&first_receipt),
        ));
    }
    receipt_pass(
        "typed_ilu_rows=missing:1,stored-zero:0; caught_panics=degree,alpha,nonsquare; same_run=byte-identical",
        &first_receipt,
    )
}

fn corruption_frame(index: usize, bit: u32, canonical: u64, corrupted: u64) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-sparse:exact-solution-reference-corruption:v1");
    push_field(&mut bytes, "canonical-exact-input-frame", &exact_inputs());
    push_u64_field(&mut bytes, "corruption-seed", CORRUPTION_SEED);
    push_text_field(&mut bytes, "reference", "exact-ilu-apply-solution");
    push_usize_field(&mut bytes, "solution-index", index);
    push_u32_field(&mut bytes, "corrupted-bit", bit);
    push_u64_field(&mut bytes, "canonical-value-bits", canonical);
    push_u64_field(&mut bytes, "corrupted-value-bits", corrupted);
    bytes
}

fn reconstruct_corruption() -> Corruption {
    let index = usize::try_from((CORRUPTION_SEED >> 8) % EXACT_SOLUTION.len() as u64)
        .expect("corruption index fits usize");
    let bit =
        u32::try_from((CORRUPTION_SEED & 0xff) % 52).expect("corruption mantissa bit fits u32");
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

fn corruption_outcome(corruption: &Corruption) -> CaseOutcome {
    let measurement = match capture("red-exact-measurement", exact_measurement) {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let computed = measurement.apply[corruption.index];
    if computed == corruption.corrupted {
        CaseOutcome::pass("disclosed exact-solution reference corruption was not detected")
    } else if computed != corruption.canonical {
        CaseOutcome::fail(format!(
            "stage=seeded-corruption-baseline-drift; seed=0x{CORRUPTION_SEED:016x}; reference=exact-ilu-apply-solution; index={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}",
            corruption.index, corruption.canonical,
        ))
    } else {
        CaseOutcome::fail(format!(
            "stage=seeded-exact-solution-reference-corruption; seed=0x{CORRUPTION_SEED:016x}; reference=exact-ilu-apply-solution; index={}; bit={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}; corrupted_bits=0x{:016x}; input_frame_len={}; input_frame_fnv1a64=0x{:016x}; input_frame={}",
            corruption.index,
            corruption.bit,
            corruption.canonical,
            corruption.corrupted,
            corruption.frame.len(),
            fnv1a64(&corruption.frame),
            hex_bytes(&corruption.frame),
        ))
        .with_evidence("crates/fs-sparse/tests/preconditioner_casebook.rs#seeded-corruption")
    }
}

fn run_red_report() -> SuiteReport {
    let corruption = reconstruct_corruption();
    assert_eq!(corruption.frame.len(), CORRUPTION_FRAME_LEN);
    assert_eq!(fnv1a64(&corruption.frame), CORRUPTION_FRAME_FNV1A64);
    Suite::new(SUITE)
        .case(
            "seeded-exact-solution-reference-bit-corruption",
            CORRUPTION_FRAME_FNV1A64,
            ToleranceSpec::Exact,
            move || corruption_outcome(&corruption),
        )
        .run()
}

#[test]
fn preconditioner_casebook_emits_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let exact = exact_inputs();
    let chebyshev = chebyshev_inputs();
    let amg = amg_inputs();
    let refusals = refusal_inputs();
    assert_eq!(exact, exact_inputs());
    assert_eq!(chebyshev, chebyshev_inputs());
    assert_eq!(amg, amg_inputs());
    assert_eq!(refusals, refusal_inputs());
    assert_eq!(exact.len(), EXACT_FRAME_LEN);
    assert_eq!(fnv1a64(&exact), EXACT_FRAME_FNV1A64);
    assert_eq!(chebyshev.len(), CHEBYSHEV_FRAME_LEN);
    assert_eq!(fnv1a64(&chebyshev), CHEBYSHEV_FRAME_FNV1A64);
    assert_eq!(amg.len(), AMG_FRAME_LEN);
    assert_eq!(fnv1a64(&amg), AMG_FRAME_FNV1A64);
    assert_eq!(refusals.len(), REFUSAL_FRAME_LEN);
    assert_eq!(fnv1a64(&refusals), REFUSAL_FRAME_FNV1A64);

    let report = Suite::new(SUITE)
        .case(
            "exact-dyadic-ilu-apply-and-one-step-pcg",
            EXACT_FRAME_FNV1A64,
            ToleranceSpec::Exact,
            exact_outcome,
        )
        .case(
            "chebyshev-diagonal-band-and-pcg-replay",
            CHEBYSHEV_FRAME_FNV1A64,
            ToleranceSpec::AbsoluteLe(CHEB_SOLUTION_CEILING),
            chebyshev_outcome,
        )
        .case(
            "genuine-multilevel-sa-amg-pcg-replay",
            AMG_FRAME_FNV1A64,
            ToleranceSpec::AbsoluteLe(AMG_SOLUTION_CEILING),
            amg_outcome,
        )
        .case(
            "preconditioner-setup-refusals",
            REFUSAL_FRAME_FNV1A64,
            ToleranceSpec::Structural,
            refusal_outcome,
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
            "exact-dyadic-ilu-apply-and-one-step-pcg",
            "chebyshev-diagonal-band-and-pcg-replay",
            "genuine-multilevel-sa-amg-pcg-replay",
            "preconditioner-setup-refusals",
        ]
    );
    assert!(report.records.iter().all(|record| {
        record.version == CASEBOOK_RECORD_VERSION
            && record.pass
            && !record.evidence.is_empty()
            && record.details.contains("output_receipt=")
    }));
    assert_eq!(report.records[0].tolerance, "exact");
    assert_eq!(report.records[1].tolerance, "abs<=1e-9");
    assert_eq!(report.records[2].tolerance, "abs<=1e-7");
    assert_eq!(report.records[3].tolerance, "structural");
}

#[test]
fn seeded_exact_solution_reference_corruption_is_stable_and_refused() {
    let first_corruption = reconstruct_corruption();
    let replay_corruption = reconstruct_corruption();
    assert_eq!(first_corruption.index, 0);
    assert_eq!(first_corruption.bit, 1);
    assert_eq!(first_corruption.index, replay_corruption.index);
    assert_eq!(first_corruption.bit, replay_corruption.bit);
    assert_eq!(first_corruption.canonical, replay_corruption.canonical);
    assert_eq!(first_corruption.corrupted, replay_corruption.corrupted);
    assert_eq!(first_corruption.frame, replay_corruption.frame);
    assert_eq!(first_corruption.frame.len(), CORRUPTION_FRAME_LEN);
    assert_eq!(fnv1a64(&first_corruption.frame), CORRUPTION_FRAME_FNV1A64);
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
        panic!("seeded corruption must produce exactly one red record");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("replayed corruption must produce exactly one red record");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert!(
        first_failure
            .details
            .contains("stage=seeded-exact-solution-reference-corruption")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(first_failure.details.contains("index=0"));
    assert!(first_failure.details.contains("bit=1"));
    assert!(first_failure.details.contains("input_frame="));

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("assert_green must refuse the disclosed reference corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("seeded-exact-solution-reference-bit-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
