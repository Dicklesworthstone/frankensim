//! FrankenScipy optimizer-oracle evidence for the shared Rosenbrock fixture.
//!
//! An exact analytic KAT checks a test-local Rosenbrock implementation against
//! pinned FrankenScipy values and proves `LbfgsState` is stationary at the exact
//! minimizer. A bounded case then compares fs-ascent L-BFGS with FrankenScipy
//! BFGS and L-BFGS-B from one disclosed global-basin start after independently
//! admitting the analytic gradient through `check_grad`.
//!
//! This is finite-fixture G0 and same-build replay evidence. FrankenScipy is a
//! pinned comparison implementation, not ground truth. Basin choice is fixed by
//! the input frame; this tranche makes no claim for arbitrary starts, global or
//! constrained optimization, stochastic methods, general tolerance
//! calibration, performance, cancellation, or fresh cross-ISA/full-G5 proof.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_ascent::{LbfgsReport, LbfgsState, StopReason, StopRule, VERSION as FS_ASCENT_VERSION};
use fs_casebook::{
    CASEBOOK_RECORD_VERSION, CaseOutcome, ReplaySpec, Suite, SuiteReport, ToleranceSpec, fnv1a64,
};
use fsci_opt::{
    ConvergenceStatus, MinimizeOptions, OptimizeMethod, OptimizeResult, check_grad, minimize,
    rosen, rosen_der,
};

const SUITE: &str = "bedrock/fs-ascent-frankenscipy-optimizer-oracle-v1";
const ORACLE_VERSION: &str = "fsci-opt/0.1.0";
const ORACLE_PIN: &str = "9e271fd734465e2b2ff755aa73ea66a7217d619b";
const PRODUCTION_API: &str = "fs_ascent::LbfgsState:new+run;strong-wolfe;infinity-gradient-stop:v1";
const ORACLE_API: &str = "fsci_opt::{rosen,rosen_der,check_grad,minimize};Bfgs+LBfgsB;strict:v1";
const FRAME_ENCODING: &str =
    "field=(tag_len:u64le,tag,payload_len:u64le,payload);numbers=le;f64=bits:v1";
const ERROR_POLICY: &str = "signed=implementation-reference;aggregate=max-absolute:v1";
const NO_CLAIM_POLICY: &str = "no-arbitrary-starts;no-basin-equivalence;no-global-or-constrained-optimization;no-stochastic-methods;no-general-tolerance-calibration;no-performance;no-cancellation;no-fresh-cross-isa:v1";

const DIMENSION: usize = 4;
const MEMORY: usize = 10;
const PRODUCTION_GRADIENT_TOLERANCE: f64 = 1.0e-9;
const PRODUCTION_MAX_ITERATIONS: usize = 4_000;
const ORACLE_TOLERANCE: f64 = 1.0e-12;
const ORACLE_MAX_ITERATIONS: usize = 5_000;
const ORACLE_GRADIENT_EPSILON: f64 = 1.0e-8;
const ORACLE_FIXTURE_ID: &str = "bedrock-rosenbrock-global-basin-v1";
const GREEN_REPLAY_COMMAND: &str = "cargo test --locked -p fs-ascent --test frankenscipy_optimizer_oracle_casebook frankenscipy_optimizer_oracle_casebook_emits_complete_green_records -- --exact --nocapture";
const RED_REPLAY_COMMAND: &str = "cargo test --locked -p fs-ascent --test frankenscipy_optimizer_oracle_casebook seeded_exact_minimizer_reference_corruption_is_stable_and_refused -- --exact --nocapture";

const GRADIENT_CHECK_CEILING: f64 = 1.0e-5;
const PRODUCTION_OBJECTIVE_CEILING: f64 = 1.0e-12;
const PRODUCTION_GRADIENT_CEILING: f64 = 1.0e-5;
const PRODUCTION_POINT_CEILING: f64 = 1.0e-4;
const ORACLE_OBJECTIVE_CEILING: f64 = 1.0e-6;
const PAIR_POINT_CEILING: f64 = 1.0e-4;
const PAIR_OBJECTIVE_CEILING: f64 = 1.0e-6;

const ZERO_POINT: [f64; DIMENSION] = [0.0; DIMENSION];
const ZERO_VALUE: f64 = 3.0;
const ZERO_GRADIENT: [f64; DIMENSION] = [-2.0, -2.0, -2.0, 0.0];
const MINIMIZER: [f64; DIMENSION] = [1.0; DIMENSION];
const MINIMUM_VALUE: f64 = 0.0;
const MINIMUM_GRADIENT: [f64; DIMENSION] = [-0.0, 0.0, 0.0, 0.0];
const GLOBAL_BASIN_START: [f64; DIMENSION] = [0.9; DIMENSION];

const CORRUPTION_SEED: u64 = 0xF5A5_0024_0000_0101;

// Filled after independent reconstruction of the framing functions without
// executing either numerical implementation.
const KAT_FRAME_LEN: usize = 1_616;
const KAT_FRAME_FNV1A64: u64 = 0x9c2f_b48a_772e_b0fa;
const ORACLE_FRAME_LEN: usize = 3_458;
const ORACLE_FRAME_FNV1A64: u64 = 0x2273_eec7_2003_e072;
const CORRUPTION_FRAME_LEN: usize = 2_902;
const CORRUPTION_FRAME_FNV1A64: u64 = 0xb809_77f8_b595_cbc2;

#[derive(Debug, Clone, PartialEq, Eq)]
struct ProductionBits {
    x: Vec<u64>,
    f: u64,
    g: Vec<u64>,
    memory: usize,
    iters: usize,
    evals: usize,
    history: Vec<u64>,
    report_reason: String,
    report_grad_norm: u64,
    report_f: u64,
    report_iters: usize,
    report_evals: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleBits {
    method: String,
    x: Vec<u64>,
    fun: Option<u64>,
    success: bool,
    status: String,
    message: String,
    nfev: usize,
    njev: usize,
    nhev: usize,
    nit: usize,
    jac: Option<Vec<u64>>,
    hess_inv: Option<Vec<Vec<u64>>>,
    maxcv: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct KatMeasurement {
    local_zero_value: u64,
    local_zero_gradient: Vec<u64>,
    oracle_zero_value: u64,
    oracle_zero_gradient: Vec<u64>,
    local_minimum_value: u64,
    local_minimum_gradient: Vec<u64>,
    oracle_minimum_value: u64,
    oracle_minimum_gradient: Vec<u64>,
    production: ProductionBits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleMeasurement {
    gradient_check: u64,
    production: ProductionBits,
    bfgs: OracleBits,
    lbfgsb: OracleBits,
}

#[derive(Debug, Clone)]
struct AgreementEvidence {
    gradient_check: f64,
    production_objective: f64,
    production_grad_norm: f64,
    production_point_error: f64,
    bfgs_objective: f64,
    bfgs_point_delta: Vec<f64>,
    bfgs_point_max: f64,
    bfgs_objective_delta: f64,
    lbfgsb_objective: f64,
    lbfgsb_point_delta: Vec<f64>,
    lbfgsb_point_max: f64,
    lbfgsb_objective_delta: f64,
}

#[derive(Debug, Clone)]
struct Corruption {
    component: usize,
    bit: u32,
    canonical: u64,
    corrupted: u64,
    frame: Vec<u8>,
}

fn local_rosen(x: &[f64]) -> f64 {
    let mut value = 0.0;
    for index in 0..x.len().saturating_sub(1) {
        let residual = x[index + 1] - x[index] * x[index];
        let offset = 1.0 - x[index];
        value += 100.0 * residual * residual + offset * offset;
    }
    value
}

fn local_rosen_gradient(x: &[f64]) -> Vec<f64> {
    if x.len() < 2 {
        return vec![0.0; x.len()];
    }
    let mut gradient = vec![0.0; x.len()];
    gradient[0] = -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]);
    for index in 1..x.len() - 1 {
        gradient[index] = 200.0 * (x[index] - x[index - 1] * x[index - 1])
            - 400.0 * x[index] * (x[index + 1] - x[index] * x[index])
            - 2.0 * (1.0 - x[index]);
    }
    let last = x.len() - 1;
    gradient[last] = 200.0 * (x[last] - x[last - 1] * x[last - 1]);
    gradient
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn bit_rows(rows: &[Vec<f64>]) -> Vec<Vec<u64>> {
    rows.iter().map(|row| bits(row)).collect()
}

fn finite_bits(values: &[u64]) -> bool {
    values
        .iter()
        .all(|&value| f64::from_bits(value).is_finite())
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().map(|value| value.abs()).fold(0.0, f64::max)
}

fn signed_delta(left: &[u64], right: &[u64]) -> Vec<f64> {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| f64::from_bits(left) - f64::from_bits(right))
        .collect()
}

fn stop_reason_name(reason: &StopReason) -> &'static str {
    match reason {
        StopReason::GradNorm => "grad-norm",
        StopReason::ObjectiveBelow => "objective-below",
        StopReason::Budget => "budget",
        StopReason::Stall => "stall",
        StopReason::Composite => "composite",
        StopReason::IterationCap => "iteration-cap",
    }
}

fn convergence_status_name(status: &ConvergenceStatus) -> &'static str {
    match status {
        ConvergenceStatus::Success => "success",
        ConvergenceStatus::MaxIterations => "max-iterations",
        ConvergenceStatus::MaxEvaluations => "max-evaluations",
        ConvergenceStatus::PrecisionLoss => "precision-loss",
        ConvergenceStatus::NanEncountered => "nan-encountered",
        ConvergenceStatus::OutOfBounds => "out-of-bounds",
        ConvergenceStatus::CallbackStop => "callback-stop",
        ConvergenceStatus::NotImplemented => "not-implemented",
        ConvergenceStatus::InvalidInput => "invalid-input",
    }
}

fn method_name(method: OptimizeMethod) -> &'static str {
    match method {
        OptimizeMethod::Bfgs => "Bfgs",
        OptimizeMethod::ConjugateGradient => "ConjugateGradient",
        OptimizeMethod::Powell => "Powell",
        OptimizeMethod::NelderMead => "NelderMead",
        OptimizeMethod::LBfgsB => "LBfgsB",
        OptimizeMethod::NewtonCg => "NewtonCg",
        OptimizeMethod::TrustExact => "TrustExact",
        OptimizeMethod::Tnc => "Tnc",
        OptimizeMethod::Slsqp => "Slsqp",
        OptimizeMethod::TrustConstr => "TrustConstr",
    }
}

fn production_bits(state: &LbfgsState, report: &LbfgsReport) -> ProductionBits {
    ProductionBits {
        x: bits(&state.x),
        f: state.f.to_bits(),
        g: bits(&state.g),
        memory: state.memory,
        iters: state.iters,
        evals: state.evals,
        history: bits(&state.history),
        report_reason: stop_reason_name(&report.reason).to_owned(),
        report_grad_norm: report.grad_norm.to_bits(),
        report_f: report.f.to_bits(),
        report_iters: report.iters,
        report_evals: report.evals,
    }
}

fn oracle_bits(method: OptimizeMethod, result: OptimizeResult) -> OracleBits {
    let OptimizeResult {
        x,
        fun,
        success,
        status,
        message,
        nfev,
        njev,
        nhev,
        nit,
        jac,
        hess_inv,
        maxcv,
    } = result;
    OracleBits {
        method: method_name(method).to_owned(),
        x: bits(&x),
        fun: fun.map(f64::to_bits),
        success,
        status: convergence_status_name(&status).to_owned(),
        message,
        nfev,
        njev,
        nhev,
        nit,
        jac: jac.map(|values| bits(&values)),
        hess_inv: hess_inv.map(|rows| bit_rows(&rows)),
        maxcv: maxcv.map(f64::to_bits),
    }
}

fn oracle_options(method: OptimizeMethod) -> MinimizeOptions {
    MinimizeOptions {
        method: Some(method),
        tol: Some(ORACLE_TOLERANCE),
        maxiter: Some(ORACLE_MAX_ITERATIONS),
        maxfev: None,
        gradient_eps: ORACLE_GRADIENT_EPSILON,
        callback: None,
        gradient: None,
        hessp: None,
        bounds: None,
        has_general_constraints: false,
        gradient_available: true,
        fixture_id: Some(ORACLE_FIXTURE_ID),
        seed: None,
        ..MinimizeOptions::default()
    }
}

fn measure_kat() -> KatMeasurement {
    let mut objective = |x: &[f64]| (local_rosen(x), local_rosen_gradient(x));
    let mut state = LbfgsState::new(&MINIMIZER, MEMORY, &mut objective);
    let report = state.run(
        &mut objective,
        &StopRule::GradNorm(PRODUCTION_GRADIENT_TOLERANCE),
        PRODUCTION_MAX_ITERATIONS,
    );
    KatMeasurement {
        local_zero_value: local_rosen(&ZERO_POINT).to_bits(),
        local_zero_gradient: bits(&local_rosen_gradient(&ZERO_POINT)),
        oracle_zero_value: rosen(&ZERO_POINT).to_bits(),
        oracle_zero_gradient: bits(&rosen_der(&ZERO_POINT)),
        local_minimum_value: local_rosen(&MINIMIZER).to_bits(),
        local_minimum_gradient: bits(&local_rosen_gradient(&MINIMIZER)),
        oracle_minimum_value: rosen(&MINIMIZER).to_bits(),
        oracle_minimum_gradient: bits(&rosen_der(&MINIMIZER)),
        production: production_bits(&state, &report),
    }
}

fn measure_oracle() -> Result<OracleMeasurement, String> {
    let gradient_check = check_grad(local_rosen, local_rosen_gradient, &GLOBAL_BASIN_START)
        .map_err(|error| format!("stage=check-grad; error={error:?}"))?;
    let mut objective = |x: &[f64]| (local_rosen(x), local_rosen_gradient(x));
    let mut state = LbfgsState::new(&GLOBAL_BASIN_START, MEMORY, &mut objective);
    let report = state.run(
        &mut objective,
        &StopRule::GradNorm(PRODUCTION_GRADIENT_TOLERANCE),
        PRODUCTION_MAX_ITERATIONS,
    );
    let bfgs = minimize(
        rosen,
        &GLOBAL_BASIN_START,
        oracle_options(OptimizeMethod::Bfgs),
    )
    .map_err(|error| format!("stage=oracle-bfgs; error={error:?}"))?;
    let lbfgsb = minimize(
        rosen,
        &GLOBAL_BASIN_START,
        oracle_options(OptimizeMethod::LBfgsB),
    )
    .map_err(|error| format!("stage=oracle-lbfgsb; error={error:?}"))?;
    Ok(OracleMeasurement {
        gradient_check: gradient_check.to_bits(),
        production: production_bits(&state, &report),
        bfgs: oracle_bits(OptimizeMethod::Bfgs, bfgs),
        lbfgsb: oracle_bits(OptimizeMethod::LBfgsB, lbfgsb),
    })
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("optimizer Casebook frame lengths fit u64"),
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
        u64::try_from(value).expect("optimizer Casebook values fit u64"),
    );
}

fn push_bool_field(bytes: &mut Vec<u8>, tag: &str, value: bool) {
    push_u32_field(bytes, tag, if value { 1 } else { 0 });
}

fn push_f64_field(bytes: &mut Vec<u8>, tag: &str, value: f64) {
    push_u64_field(bytes, tag, value.to_bits());
}

fn push_bits_field(bytes: &mut Vec<u8>, tag: &str, values: &[u64]) {
    let mut payload = Vec::with_capacity(8 + values.len() * 8);
    push_len(&mut payload, values.len());
    for &value in values {
        push_u64(&mut payload, value);
    }
    push_field(bytes, tag, &payload);
}

fn push_f64s_field(bytes: &mut Vec<u8>, tag: &str, values: &[f64]) {
    push_bits_field(bytes, tag, &bits(values));
}

fn push_bit_rows(bytes: &mut Vec<u8>, tag: &str, rows: &[Vec<u64>]) {
    let mut payload = Vec::new();
    push_len(&mut payload, rows.len());
    for row in rows {
        push_bits_field(&mut payload, "row", row);
    }
    push_field(bytes, tag, &payload);
}

fn push_optional_bits(bytes: &mut Vec<u8>, tag: &str, values: Option<&[u64]>) {
    push_bool_field(bytes, &format!("{tag}-present"), values.is_some());
    if let Some(values) = values {
        push_bits_field(bytes, tag, values);
    }
}

fn push_optional_f64(bytes: &mut Vec<u8>, tag: &str, value: Option<u64>) {
    push_bool_field(bytes, &format!("{tag}-present"), value.is_some());
    if let Some(value) = value {
        push_u64_field(bytes, tag, value);
    }
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
    push_text_field(&mut bytes, "fs-ascent-version", FS_ASCENT_VERSION);
    push_text_field(&mut bytes, "production-api", PRODUCTION_API);
    push_text_field(&mut bytes, "oracle-version", ORACLE_VERSION);
    push_text_field(&mut bytes, "oracle-pin", ORACLE_PIN);
    push_text_field(&mut bytes, "oracle-api", ORACLE_API);
    push_text_field(&mut bytes, "error-policy", ERROR_POLICY);
    push_text_field(&mut bytes, "no-claim-policy", NO_CLAIM_POLICY);
    bytes
}

fn push_production_options(bytes: &mut Vec<u8>) {
    push_usize_field(bytes, "production-memory", MEMORY);
    push_text_field(bytes, "production-stop-rule", "GradNorm");
    push_f64_field(
        bytes,
        "production-gradient-tolerance",
        PRODUCTION_GRADIENT_TOLERANCE,
    );
    push_usize_field(bytes, "production-iteration-cap", PRODUCTION_MAX_ITERATIONS);
    push_text_field(
        bytes,
        "production-line-search",
        "strong-wolfe:c1=1e-4:c2=0.9:alpha0=1:v1",
    );
}

fn push_oracle_options(bytes: &mut Vec<u8>, method: OptimizeMethod) {
    push_text_field(bytes, "oracle-method", method_name(method));
    push_bool_field(bytes, "oracle-tolerance-present", true);
    push_f64_field(bytes, "oracle-tolerance", ORACLE_TOLERANCE);
    push_bool_field(bytes, "oracle-max-iterations-present", true);
    push_usize_field(bytes, "oracle-max-iterations", ORACLE_MAX_ITERATIONS);
    push_bool_field(bytes, "oracle-max-evaluations-present", false);
    push_f64_field(bytes, "oracle-gradient-epsilon", ORACLE_GRADIENT_EPSILON);
    push_bool_field(bytes, "oracle-callback-present", false);
    push_bool_field(bytes, "oracle-gradient-callback-present", false);
    push_bool_field(bytes, "oracle-hessian-product-present", false);
    push_bool_field(bytes, "oracle-bounds-present", false);
    push_bool_field(bytes, "oracle-general-constraints-present", false);
    push_bool_field(bytes, "oracle-gradient-available", true);
    push_bool_field(bytes, "oracle-fixture-id-present", true);
    push_text_field(bytes, "oracle-fixture-id", ORACLE_FIXTURE_ID);
    push_bool_field(bytes, "oracle-seed-present", false);
    push_text_field(bytes, "oracle-runtime-mode", "Strict");
}

fn kat_inputs() -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-ascent:rosenbrock-analytic-kat-input:v1");
    push_text_field(
        &mut bytes,
        "objective",
        "sum_i(100*(x[i+1]-x[i]^2)^2+(1-x[i])^2):fixed-order:v1",
    );
    push_usize_field(&mut bytes, "dimension", DIMENSION);
    push_f64s_field(&mut bytes, "zero-point", &ZERO_POINT);
    push_f64_field(&mut bytes, "zero-expected-value", ZERO_VALUE);
    push_f64s_field(&mut bytes, "zero-expected-gradient", &ZERO_GRADIENT);
    push_f64s_field(&mut bytes, "exact-minimizer", &MINIMIZER);
    push_f64_field(&mut bytes, "minimum-expected-value", MINIMUM_VALUE);
    push_f64s_field(&mut bytes, "minimum-expected-gradient", &MINIMUM_GRADIENT);
    push_production_options(&mut bytes);
    bytes
}

fn oracle_inputs() -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-ascent:rosenbrock-optimizer-oracle-input:v1");
    push_text_field(
        &mut bytes,
        "objective",
        "sum_i(100*(x[i+1]-x[i]^2)^2+(1-x[i])^2):fixed-order:v1",
    );
    push_text_field(
        &mut bytes,
        "basin-policy",
        "single-disclosed-global-basin-start;no-basin-equivalence:v1",
    );
    push_usize_field(&mut bytes, "dimension", DIMENSION);
    push_f64s_field(&mut bytes, "start", &GLOBAL_BASIN_START);
    push_f64s_field(&mut bytes, "exact-minimizer", &MINIMIZER);
    push_production_options(&mut bytes);
    push_text_field(
        &mut bytes,
        "gradient-check",
        "fsci_opt::check_grad:forward-difference-step=1.4901161193847656e-8:v1",
    );
    push_f64_field(&mut bytes, "gradient-check-ceiling", GRADIENT_CHECK_CEILING);
    push_oracle_options(&mut bytes, OptimizeMethod::Bfgs);
    push_oracle_options(&mut bytes, OptimizeMethod::LBfgsB);
    push_f64_field(
        &mut bytes,
        "production-objective-ceiling",
        PRODUCTION_OBJECTIVE_CEILING,
    );
    push_f64_field(
        &mut bytes,
        "production-gradient-ceiling",
        PRODUCTION_GRADIENT_CEILING,
    );
    push_f64_field(
        &mut bytes,
        "production-point-ceiling",
        PRODUCTION_POINT_CEILING,
    );
    push_f64_field(
        &mut bytes,
        "oracle-objective-ceiling",
        ORACLE_OBJECTIVE_CEILING,
    );
    push_f64_field(&mut bytes, "pair-point-ceiling", PAIR_POINT_CEILING);
    push_f64_field(&mut bytes, "pair-objective-ceiling", PAIR_OBJECTIVE_CEILING);
    bytes
}

fn push_production(bytes: &mut Vec<u8>, production: &ProductionBits) {
    push_bits_field(bytes, "iterate-bits", &production.x);
    push_u64_field(bytes, "objective-bits", production.f);
    push_bits_field(bytes, "gradient-bits", &production.g);
    push_usize_field(bytes, "memory", production.memory);
    push_usize_field(bytes, "iterations", production.iters);
    push_usize_field(bytes, "evaluations", production.evals);
    push_bits_field(bytes, "objective-history-bits", &production.history);
    push_text_field(bytes, "report-reason", &production.report_reason);
    push_u64_field(
        bytes,
        "report-gradient-norm-bits",
        production.report_grad_norm,
    );
    push_u64_field(bytes, "report-objective-bits", production.report_f);
    push_usize_field(bytes, "report-iterations", production.report_iters);
    push_usize_field(bytes, "report-evaluations", production.report_evals);
}

fn push_oracle(bytes: &mut Vec<u8>, oracle: &OracleBits) {
    push_text_field(bytes, "method", &oracle.method);
    push_bits_field(bytes, "iterate-bits", &oracle.x);
    push_optional_f64(bytes, "objective-bits", oracle.fun);
    push_bool_field(bytes, "success", oracle.success);
    push_text_field(bytes, "status", &oracle.status);
    push_text_field(bytes, "message", &oracle.message);
    push_usize_field(bytes, "function-evaluations", oracle.nfev);
    push_usize_field(bytes, "gradient-evaluations", oracle.njev);
    push_usize_field(bytes, "hessian-evaluations", oracle.nhev);
    push_usize_field(bytes, "iterations", oracle.nit);
    push_optional_bits(bytes, "jacobian-bits", oracle.jac.as_deref());
    push_bool_field(bytes, "inverse-hessian-present", oracle.hess_inv.is_some());
    if let Some(rows) = &oracle.hess_inv {
        push_bit_rows(bytes, "inverse-hessian-row-bits", rows);
    }
    push_optional_f64(bytes, "max-constraint-violation-bits", oracle.maxcv);
}

fn kat_receipt(inputs: &[u8], measurement: &KatMeasurement) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-ascent:rosenbrock-analytic-kat-output:v1");
    push_field(&mut bytes, "canonical-input-frame", inputs);
    push_u64_field(
        &mut bytes,
        "local-zero-value-bits",
        measurement.local_zero_value,
    );
    push_bits_field(
        &mut bytes,
        "local-zero-gradient-bits",
        &measurement.local_zero_gradient,
    );
    push_u64_field(
        &mut bytes,
        "oracle-zero-value-bits",
        measurement.oracle_zero_value,
    );
    push_bits_field(
        &mut bytes,
        "oracle-zero-gradient-bits",
        &measurement.oracle_zero_gradient,
    );
    push_u64_field(
        &mut bytes,
        "local-minimum-value-bits",
        measurement.local_minimum_value,
    );
    push_bits_field(
        &mut bytes,
        "local-minimum-gradient-bits",
        &measurement.local_minimum_gradient,
    );
    push_u64_field(
        &mut bytes,
        "oracle-minimum-value-bits",
        measurement.oracle_minimum_value,
    );
    push_bits_field(
        &mut bytes,
        "oracle-minimum-gradient-bits",
        &measurement.oracle_minimum_gradient,
    );
    let mut production = Vec::new();
    push_production(&mut production, &measurement.production);
    push_field(&mut bytes, "production-stationary-output", &production);
    bytes
}

fn oracle_receipt(
    inputs: &[u8],
    measurement: &OracleMeasurement,
    evidence: &AgreementEvidence,
) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-ascent:rosenbrock-optimizer-oracle-output:v1");
    push_field(&mut bytes, "canonical-input-frame", inputs);
    push_u64_field(
        &mut bytes,
        "gradient-check-bits",
        measurement.gradient_check,
    );
    let mut production = Vec::new();
    push_production(&mut production, &measurement.production);
    push_field(&mut bytes, "production-output", &production);
    let mut bfgs = Vec::new();
    push_oracle(&mut bfgs, &measurement.bfgs);
    push_field(&mut bytes, "bfgs-output", &bfgs);
    let mut lbfgsb = Vec::new();
    push_oracle(&mut lbfgsb, &measurement.lbfgsb);
    push_field(&mut bytes, "lbfgsb-output", &lbfgsb);
    for (tag, value, ceiling) in [
        (
            "gradient-check",
            evidence.gradient_check,
            GRADIENT_CHECK_CEILING,
        ),
        (
            "production-objective",
            evidence.production_objective,
            PRODUCTION_OBJECTIVE_CEILING,
        ),
        (
            "production-gradient-norm",
            evidence.production_grad_norm,
            PRODUCTION_GRADIENT_CEILING,
        ),
        (
            "production-point-error",
            evidence.production_point_error,
            PRODUCTION_POINT_CEILING,
        ),
        (
            "bfgs-objective",
            evidence.bfgs_objective,
            ORACLE_OBJECTIVE_CEILING,
        ),
        (
            "bfgs-point-max",
            evidence.bfgs_point_max,
            PAIR_POINT_CEILING,
        ),
        (
            "bfgs-objective-delta",
            evidence.bfgs_objective_delta,
            PAIR_OBJECTIVE_CEILING,
        ),
        (
            "lbfgsb-objective",
            evidence.lbfgsb_objective,
            ORACLE_OBJECTIVE_CEILING,
        ),
        (
            "lbfgsb-point-max",
            evidence.lbfgsb_point_max,
            PAIR_POINT_CEILING,
        ),
        (
            "lbfgsb-objective-delta",
            evidence.lbfgsb_objective_delta,
            PAIR_OBJECTIVE_CEILING,
        ),
    ] {
        push_f64_field(&mut bytes, tag, value);
        push_f64_field(&mut bytes, &format!("{tag}-ceiling"), ceiling);
        push_f64_field(
            &mut bytes,
            &format!("{tag}-remaining-margin"),
            ceiling - value,
        );
    }
    push_f64s_field(
        &mut bytes,
        "production-minus-bfgs-point",
        &evidence.bfgs_point_delta,
    );
    push_f64s_field(
        &mut bytes,
        "production-minus-lbfgsb-point",
        &evidence.lbfgsb_point_delta,
    );
    bytes
}

fn finite_oracle(oracle: &OracleBits) -> bool {
    finite_bits(&oracle.x)
        && oracle.fun.is_none_or(|value| finite_bits(&[value]))
        && oracle.jac.as_ref().is_none_or(|values| finite_bits(values))
        && oracle
            .hess_inv
            .as_ref()
            .is_none_or(|rows| rows.iter().all(|row| finite_bits(row)))
        && oracle.maxcv.is_none_or(|value| finite_bits(&[value]))
}

fn admit_kat(measurement: &KatMeasurement) -> Result<(), String> {
    let zero_value = ZERO_VALUE.to_bits();
    let zero_gradient = bits(&ZERO_GRADIENT);
    let minimum_value = MINIMUM_VALUE.to_bits();
    let minimum_gradient = bits(&MINIMUM_GRADIENT);
    let minimizer = bits(&MINIMIZER);
    if measurement.local_zero_value != zero_value
        || measurement.oracle_zero_value != zero_value
        || measurement.local_zero_gradient != zero_gradient
        || measurement.oracle_zero_gradient != zero_gradient
        || measurement.local_minimum_value != minimum_value
        || measurement.oracle_minimum_value != minimum_value
        || measurement.local_minimum_gradient != minimum_gradient
        || measurement.oracle_minimum_gradient != minimum_gradient
    {
        return Err(format!(
            "stage=analytic-known-answer; measurement={measurement:?}; zero_value=0x{zero_value:016x}; zero_gradient={zero_gradient:016x?}; minimum_value=0x{minimum_value:016x}; minimum_gradient={minimum_gradient:016x?}"
        ));
    }
    let production = &measurement.production;
    if production.x != minimizer
        || production.f != minimum_value
        || production.g != minimum_gradient
        || production.memory != MEMORY
        || production.iters != 0
        || production.evals != 1
        || production.history != [minimum_value]
        || production.report_reason != "grad-norm"
        || production.report_grad_norm != minimum_value
        || production.report_f != minimum_value
        || production.report_iters != 0
        || production.report_evals != 1
    {
        return Err(format!(
            "stage=production-stationary-known-answer; production={production:?}; minimizer={minimizer:016x?}; minimum_gradient={minimum_gradient:016x?}"
        ));
    }
    Ok(())
}

fn oracle_objective(oracle: &OracleBits) -> Result<f64, String> {
    oracle.fun.map(f64::from_bits).ok_or_else(|| {
        format!(
            "stage=oracle-missing-objective; method={}; output={oracle:?}",
            oracle.method
        )
    })
}

fn admit_oracle_result(oracle: &OracleBits) -> Result<f64, String> {
    // At the authoritative constellation pin, BFGS can report PrecisionLoss
    // after reaching an admissible endpoint, while L-BFGS-B reports both an
    // exhausted line search and a true iteration cap as MaxIterations. Keep
    // those diagnostics in the complete receipt and let the explicit endpoint
    // bounds below decide agreement; reject every unrelated terminal state.
    let terminal_state_is_admitted = match oracle.method.as_str() {
        "Bfgs" => matches!(
            (oracle.success, oracle.status.as_str()),
            (true, "success") | (false, "precision-loss")
        ),
        "LBfgsB" => matches!(
            (oracle.success, oracle.status.as_str()),
            (true, "success") | (false, "max-iterations")
        ),
        _ => false,
    };
    let jacobian_shape_is_admitted = oracle
        .jac
        .as_ref()
        .is_some_and(|jacobian| jacobian.len() == DIMENSION);
    let curvature_shape_is_admitted = match oracle.method.as_str() {
        "Bfgs" => oracle.hess_inv.as_ref().is_some_and(|rows| {
            rows.len() == DIMENSION && rows.iter().all(|row| row.len() == DIMENSION)
        }),
        "LBfgsB" => oracle.hess_inv.is_none(),
        _ => false,
    };
    if oracle.x.len() != DIMENSION
        || !finite_oracle(oracle)
        || !terminal_state_is_admitted
        || !jacobian_shape_is_admitted
        || !curvature_shape_is_admitted
        || oracle.message.is_empty()
        || oracle.nfev == 0
        || oracle.njev == 0
        || oracle.nhev != 0
        || oracle.maxcv.is_some()
        || oracle.nit > ORACLE_MAX_ITERATIONS
    {
        return Err(format!(
            "stage=oracle-admission; method={}; output={oracle:?}",
            oracle.method
        ));
    }
    oracle_objective(oracle)
}

#[allow(clippy::too_many_lines)] // Ordered admission keeps every public field fail-closed.
fn admit_oracle(measurement: &OracleMeasurement) -> Result<AgreementEvidence, String> {
    let production = &measurement.production;
    if production.x.len() != DIMENSION
        || production.g.len() != DIMENSION
        || production.history.len() != production.iters + 1
        || !finite_bits(&production.x)
        || !finite_bits(&production.g)
        || !finite_bits(&production.history)
        || !finite_bits(&[
            production.f,
            production.report_grad_norm,
            production.report_f,
        ])
        || production.memory != MEMORY
        || production.evals == 0
        || production.iters > PRODUCTION_MAX_ITERATIONS
        || production.report_reason == "iteration-cap"
        || production.f != production.report_f
        || production.iters != production.report_iters
        || production.evals != production.report_evals
    {
        return Err(format!(
            "stage=production-admission; production={production:?}"
        ));
    }
    let observed_grad_norm = max_abs(
        &production
            .g
            .iter()
            .map(|&value| f64::from_bits(value))
            .collect::<Vec<_>>(),
    );
    if observed_grad_norm.to_bits() != production.report_grad_norm {
        return Err(format!(
            "stage=production-gradient-certificate; observed_bits=0x{:016x}; report_bits=0x{:016x}; production={production:?}",
            observed_grad_norm.to_bits(),
            production.report_grad_norm,
        ));
    }
    let gradient_check = f64::from_bits(measurement.gradient_check);
    let production_objective = f64::from_bits(production.f);
    let production_grad_norm = f64::from_bits(production.report_grad_norm);
    let production_point_error = max_abs(&signed_delta(&production.x, &bits(&MINIMIZER)));
    let bfgs_objective = admit_oracle_result(&measurement.bfgs)?;
    let lbfgsb_objective = admit_oracle_result(&measurement.lbfgsb)?;
    let bfgs_point_delta = signed_delta(&production.x, &measurement.bfgs.x);
    let bfgs_point_max = max_abs(&bfgs_point_delta);
    let bfgs_objective_delta = (bfgs_objective - production_objective).abs();
    let lbfgsb_point_delta = signed_delta(&production.x, &measurement.lbfgsb.x);
    let lbfgsb_point_max = max_abs(&lbfgsb_point_delta);
    let lbfgsb_objective_delta = (lbfgsb_objective - production_objective).abs();
    let evidence = AgreementEvidence {
        gradient_check,
        production_objective,
        production_grad_norm,
        production_point_error,
        bfgs_objective,
        bfgs_point_delta,
        bfgs_point_max,
        bfgs_objective_delta,
        lbfgsb_objective,
        lbfgsb_point_delta,
        lbfgsb_point_max,
        lbfgsb_objective_delta,
    };
    let gates = [
        ("gradient-check", gradient_check, GRADIENT_CHECK_CEILING),
        (
            "production-objective",
            production_objective,
            PRODUCTION_OBJECTIVE_CEILING,
        ),
        (
            "production-gradient",
            production_grad_norm,
            PRODUCTION_GRADIENT_CEILING,
        ),
        (
            "production-point",
            production_point_error,
            PRODUCTION_POINT_CEILING,
        ),
        ("bfgs-objective", bfgs_objective, ORACLE_OBJECTIVE_CEILING),
        ("bfgs-point", bfgs_point_max, PAIR_POINT_CEILING),
        (
            "bfgs-objective-delta",
            bfgs_objective_delta,
            PAIR_OBJECTIVE_CEILING,
        ),
        (
            "lbfgsb-objective",
            lbfgsb_objective,
            ORACLE_OBJECTIVE_CEILING,
        ),
        ("lbfgsb-point", lbfgsb_point_max, PAIR_POINT_CEILING),
        (
            "lbfgsb-objective-delta",
            lbfgsb_objective_delta,
            PAIR_OBJECTIVE_CEILING,
        ),
    ];
    if let Some((name, value, ceiling)) = gates
        .into_iter()
        .find(|(_, value, ceiling)| !value.is_finite() || value > ceiling)
    {
        return Err(format!(
            "stage=bounded-oracle-agreement; failed={name}; value={value}; ceiling={ceiling}; gates={gates:?}; measurement={measurement:?}"
        ));
    }
    Ok(evidence)
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

fn capture_kat(stage: &str) -> Result<KatMeasurement, String> {
    catch_unwind(AssertUnwindSafe(measure_kat))
        .map_err(|payload| format!("stage={stage}; panic={}", panic_message(payload.as_ref())))
}

fn capture_oracle(stage: &str) -> Result<OracleMeasurement, String> {
    match catch_unwind(AssertUnwindSafe(measure_oracle)) {
        Ok(result) => result,
        Err(payload) => Err(format!(
            "stage={stage}; panic={}",
            panic_message(payload.as_ref())
        )),
    }
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn kat_outcome() -> CaseOutcome {
    let inputs = kat_inputs();
    let first = match capture_kat("first-kat-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay = match capture_kat("replay-kat-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    if let Err(error) = admit_kat(&first) {
        return CaseOutcome::fail(error);
    }
    if let Err(error) = admit_kat(&replay) {
        return CaseOutcome::fail(format!("stage=replay-admission; {error}"));
    }
    let first_receipt = kat_receipt(&inputs, &first);
    let replay_receipt = kat_receipt(&inputs, &replay);
    if first_receipt != replay_receipt {
        return CaseOutcome::fail(format!(
            "stage=kat-output-replay; first_len={}; replay_len={}; first_fnv1a64=0x{:016x}; replay_fnv1a64=0x{:016x}; first={}; replay={}",
            first_receipt.len(),
            replay_receipt.len(),
            fnv1a64(&first_receipt),
            fnv1a64(&replay_receipt),
            hex_bytes(&first_receipt),
            hex_bytes(&replay_receipt),
        ));
    }
    CaseOutcome::pass(format!(
        "zero_value=3; zero_gradient=[-2,-2,-2,0]; exact_minimizer=[1,1,1,1]; production_stationary=true; output_receipt_len={}; output_receipt_fnv1a64=0x{:016x}; output_receipt={}",
        first_receipt.len(),
        fnv1a64(&first_receipt),
        hex_bytes(&first_receipt),
    ))
    .with_evidence("crates/fs-ascent/tests/frankenscipy_optimizer_oracle_casebook.rs#analytic-kat")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn oracle_outcome() -> CaseOutcome {
    let inputs = oracle_inputs();
    let first = match capture_oracle("first-oracle-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay = match capture_oracle("replay-oracle-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let first_evidence = match admit_oracle(&first) {
        Ok(evidence) => evidence,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay_evidence = match admit_oracle(&replay) {
        Ok(evidence) => evidence,
        Err(error) => return CaseOutcome::fail(format!("stage=replay-admission; {error}")),
    };
    let first_receipt = oracle_receipt(&inputs, &first, &first_evidence);
    let replay_receipt = oracle_receipt(&inputs, &replay, &replay_evidence);
    if first_receipt != replay_receipt {
        return CaseOutcome::fail(format!(
            "stage=oracle-output-replay; first_len={}; replay_len={}; first_fnv1a64=0x{:016x}; replay_fnv1a64=0x{:016x}; first={}; replay={}",
            first_receipt.len(),
            replay_receipt.len(),
            fnv1a64(&first_receipt),
            fnv1a64(&replay_receipt),
            hex_bytes(&first_receipt),
            hex_bytes(&replay_receipt),
        ));
    }
    CaseOutcome::pass(format!(
        "gradient_check={:.3e}; production_f={:.3e}; production_grad={:.3e}; production_point={:.3e}; bfgs_f={:.3e}; bfgs_point={:.3e}; lbfgsb_f={:.3e}; lbfgsb_point={:.3e}; output_receipt_len={}; output_receipt_fnv1a64=0x{:016x}; output_receipt={}",
        first_evidence.gradient_check,
        first_evidence.production_objective,
        first_evidence.production_grad_norm,
        first_evidence.production_point_error,
        first_evidence.bfgs_objective,
        first_evidence.bfgs_point_max,
        first_evidence.lbfgsb_objective,
        first_evidence.lbfgsb_point_max,
        first_receipt.len(),
        fnv1a64(&first_receipt),
        hex_bytes(&first_receipt),
    ))
    .with_evidence("crates/fs-ascent/tests/frankenscipy_optimizer_oracle_casebook.rs#bounded-oracle")
    .with_evidence("constellation.lock:frankenscipy-0.1.0")
}

fn corruption_frame(component: usize, bit: u32, canonical: u64, corrupted: u64) -> Vec<u8> {
    let kat = kat_inputs();
    let mut bytes = common_frame("bedrock:fs-ascent:rosenbrock-reference-corruption:v1");
    push_field(&mut bytes, "canonical-kat-input-frame", &kat);
    push_u64_field(&mut bytes, "canonical-kat-input-fnv1a64", fnv1a64(&kat));
    push_u64_field(&mut bytes, "corruption-seed", CORRUPTION_SEED);
    push_text_field(
        &mut bytes,
        "selection-policy",
        "component=(seed>>8)%dimension;bit=(seed&0xff)%52:v1",
    );
    push_text_field(
        &mut bytes,
        "reference",
        "exact-stationary-minimizer-component",
    );
    push_usize_field(&mut bytes, "component", component);
    push_u32_field(&mut bytes, "mantissa-bit", bit);
    push_u64_field(&mut bytes, "canonical-value-bits", canonical);
    push_u64_field(&mut bytes, "corrupted-value-bits", corrupted);
    bytes
}

fn reconstruct_corruption() -> Corruption {
    let component = usize::try_from((CORRUPTION_SEED >> 8) % DIMENSION as u64)
        .expect("corruption component fits usize");
    let bit = u32::try_from((CORRUPTION_SEED & 0xff) % 52).expect("corruption bit fits u32");
    let canonical = MINIMIZER[component].to_bits();
    let corrupted = canonical ^ (1_u64 << bit);
    Corruption {
        component,
        bit,
        canonical,
        corrupted,
        frame: corruption_frame(component, bit, canonical, corrupted),
    }
}

fn corruption_outcome(corruption: Corruption) -> CaseOutcome {
    let measurement = match capture_kat("red-baseline-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    if let Err(error) = admit_kat(&measurement) {
        return CaseOutcome::fail(format!("stage=red-baseline-admission; {error}"));
    }
    let computed = measurement.production.x[corruption.component];
    if computed == corruption.corrupted {
        CaseOutcome::pass("disclosed exact-minimizer corruption was not detected")
    } else if computed != corruption.canonical {
        CaseOutcome::fail(format!(
            "stage=seeded-corruption-baseline-drift; seed=0x{CORRUPTION_SEED:016x}; component={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}",
            corruption.component, corruption.canonical,
        ))
    } else {
        CaseOutcome::fail(format!(
            "stage=seeded-exact-minimizer-reference-corruption; seed=0x{CORRUPTION_SEED:016x}; component={}; bit={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}; corrupted_bits=0x{:016x}; input_frame_len={}; input_frame_fnv1a64=0x{:016x}; input_frame={}",
            corruption.component,
            corruption.bit,
            corruption.canonical,
            corruption.corrupted,
            corruption.frame.len(),
            fnv1a64(&corruption.frame),
            hex_bytes(&corruption.frame),
        ))
        .with_evidence("crates/fs-ascent/tests/frankenscipy_optimizer_oracle_casebook.rs#seeded-corruption")
    }
}

fn run_red_report() -> SuiteReport {
    let corruption = reconstruct_corruption();
    assert_eq!(corruption.frame.len(), CORRUPTION_FRAME_LEN);
    assert_eq!(fnv1a64(&corruption.frame), CORRUPTION_FRAME_FNV1A64);
    let replay = ReplaySpec::new(RED_REPLAY_COMMAND, corruption.frame.clone());
    Suite::new(SUITE)
        .case_replayable(
            "seeded-exact-minimizer-reference-bit-corruption",
            replay,
            ToleranceSpec::Exact,
            move || corruption_outcome(corruption),
        )
        .run()
}

#[test]
fn frankenscipy_optimizer_oracle_casebook_emits_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let kat = kat_inputs();
    let oracle = oracle_inputs();
    assert_eq!(kat, kat_inputs());
    assert_eq!(oracle, oracle_inputs());
    assert_eq!(kat.len(), KAT_FRAME_LEN);
    assert_eq!(fnv1a64(&kat), KAT_FRAME_FNV1A64);
    assert_eq!(oracle.len(), ORACLE_FRAME_LEN);
    assert_eq!(fnv1a64(&oracle), ORACLE_FRAME_FNV1A64);

    let report = Suite::new(SUITE)
        .case_replayable(
            "rosenbrock-analytic-and-stationary-kat",
            ReplaySpec::new(GREEN_REPLAY_COMMAND, kat.clone()),
            ToleranceSpec::Exact,
            kat_outcome,
        )
        .case_replayable(
            "rosenbrock-global-basin-optimizer-agreement",
            ReplaySpec::new(GREEN_REPLAY_COMMAND, oracle.clone()),
            ToleranceSpec::AbsoluteLe(PAIR_POINT_CEILING),
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
            "rosenbrock-analytic-and-stationary-kat",
            "rosenbrock-global-basin-optimizer-agreement",
        ]
    );
    assert!(report.records.iter().all(|record| {
        record.version == CASEBOOK_RECORD_VERSION
            && record.pass
            && !record.evidence.is_empty()
            && record.details.contains("output_receipt=")
    }));
    assert_eq!(report.records[0].tolerance, "exact");
    assert_eq!(report.records[1].tolerance, "abs<=1e-4");
    assert_eq!(report.replay_records.len(), 2);
    assert_eq!(
        report.replay_records[0]
            .verify_and_decode()
            .expect("KAT replay frame verifies"),
        kat,
    );
    assert_eq!(
        report.replay_records[1]
            .verify_and_decode()
            .expect("oracle replay frame verifies"),
        oracle,
    );
}

#[test]
fn seeded_exact_minimizer_reference_corruption_is_stable_and_refused() {
    let first_corruption = reconstruct_corruption();
    let replay_corruption = reconstruct_corruption();
    assert_eq!(first_corruption.component, 1);
    assert_eq!(first_corruption.bit, 1);
    assert_eq!(first_corruption.component, replay_corruption.component);
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
    assert_eq!(first.replay_records.len(), 1);
    assert_eq!(replay.replay_records.len(), 1);
    assert_eq!(
        first.replay_records[0]
            .verify_and_decode()
            .expect("red replay frame verifies"),
        first_corruption.frame,
    );
    assert_eq!(
        first.replay_records[0].json_line(),
        replay.replay_records[0].json_line(),
    );
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
            .contains("stage=seeded-exact-minimizer-reference-corruption")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(first_failure.details.contains("component=1"));
    assert!(first_failure.details.contains("bit=1"));
    assert!(first_failure.details.contains("input_frame="));

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("assert_green must refuse the disclosed minimizer corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("seeded-exact-minimizer-reference-bit-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
