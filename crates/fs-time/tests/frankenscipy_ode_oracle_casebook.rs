//! FrankenScipy RK45 endpoint-oracle evidence for `fs-time`.
//!
//! A stationary exact fixture proves nontrivial time advancement preserves a
//! dyadic state bit-for-bit. Decay and oscillator fixtures compare production
//! `rk45_adaptive`, pinned FrankenScipy `solve_ivp`, and independent
//! closed-form endpoints under explicit max-norm ceilings. Canonical frames
//! bind every fixture, dimensionless unit convention, and solver option;
//! receipts retain every public output field and all signed comparison
//! errors. The declared oracle pin is checked against `constellation.lock`;
//! proving that the sibling path checkout is at that pin remains the external
//! `xtask check-constellation`/DSR admission precondition.
//!
//! This is finite-fixture G0 and same-build replay evidence. FrankenScipy is a
//! pinned comparison implementation, not ground truth, and the test-side
//! standard-library analytic functions are not a cross-ISA bit oracle. No
//! claim is made for general tolerance calibration, controller-path equality,
//! dense output, events, backward integration, stiffness, cancellation,
//! checkpoint equivalence, adjoints, performance, Python SciPy, or fresh
//! cross-ISA/full-G5 execution.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_casebook::{
    CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, SuiteReport, ToleranceSpec, fnv1a64,
};
use fs_time::{AdaptiveState, PiController, VERSION as FS_TIME_VERSION, rk45_adaptive};
use fsci_integrate::{SolveIvpOptions, SolveIvpResult, SolverKind, ToleranceValue, solve_ivp};

const SUITE: &str = "bedrock/fs-time-frankenscipy-rk45-oracle-v1";
const ORACLE_VERSION: &str = "fsci-integrate/0.2.0";
const ORACLE_LOCK_VERSION: &str = "0.2.0";
const ORACLE_PIN: &str = "06108aea50fbb1541863101d371f9cd503e283d2";
const CONSTELLATION_LOCK: &str = include_str!("../../../constellation.lock");
const PRODUCTION_API: &str = "fs_time::rk45_adaptive:Dormand-Prince-5(4)+PI:max-norm:v1";
const ORACLE_API: &str = "fsci_integrate::solve_ivp:SolverKind::Rk45:t-eval-endpoints:RMS-error:v1";
const FRAME_ENCODING: &str =
    "field=(tag_len:u64le,tag,payload_len:u64le,payload);numbers=le;f64=bits:v1";
const UNIT_POLICY: &str = "time=dimensionless;state-components=dimensionless;rhs=state/time;rtol=dimensionless;atol=state:v1";
const ERROR_POLICY: &str = "signed=implementation-reference;aggregate=max-absolute-component:v1";
const ANALYTIC_POLICY: &str = "test-side-std-exp-sin-cos;finite-fixture-only;not-cross-isa-bits:v1";
const NO_CLAIM_POLICY: &str = "no-general-tolerance-calibration;no-controller-path-or-counter-equality;no-dense-output;no-events;no-backward;no-stiffness;no-cancellation;no-checkpoint-equivalence;no-adjoints;no-performance;no-python-scipy;no-fresh-cross-isa:v1";

const RTOL: f64 = 1.0e-9;
const ATOL: f64 = 1.0e-12;
const INITIAL_STEP: f64 = 0.05;
const ORACLE_MAX_STEP: f64 = f64::INFINITY;
const MAX_STEPS: usize = 100_000;
const DECAY_CEILING: f64 = 2.0e-7;
const OSCILLATOR_CEILING: f64 = 5.0e-7;
const CORRUPTION_SEED: u64 = 0xF5A5_0022_0000_0101;

// Filled from the framing code and independently reconstructed without
// executing either numerical implementation.
const STATIONARY_FRAME_LEN: usize = 2_183;
const STATIONARY_FRAME_FNV1A64: u64 = 0x3160_f7ed_184b_ddcb;
const DECAY_FRAME_LEN: usize = 2_182;
const DECAY_FRAME_FNV1A64: u64 = 0x53ef_99d7_09e0_8492;
const OSCILLATOR_FRAME_LEN: usize = 2_234;
const OSCILLATOR_FRAME_FNV1A64: u64 = 0x8d21_dcc9_ed3d_1a1d;
const CORRUPTION_FRAME_LEN: usize = 4_392;
const CORRUPTION_FRAME_FNV1A64: u64 = 0x14a2_1010_8a3d_949e;

const STATIONARY_Y0: [f64; 2] = [1.5, -0.25];
const DECAY_Y0: [f64; 1] = [2.0];
const OSCILLATOR_Y0: [f64; 2] = [0.75, -0.5];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FixtureKind {
    Stationary,
    Decay,
    Oscillator,
}

#[derive(Debug, Clone, Copy)]
struct FixtureSpec {
    kind: FixtureKind,
    id: &'static str,
    rhs_id: &'static str,
    analytic_id: &'static str,
    t0: f64,
    t_end: f64,
    y0: &'static [f64],
    ceiling: f64,
}

const STATIONARY: FixtureSpec = FixtureSpec {
    kind: FixtureKind::Stationary,
    id: "stationary-two-component-exact",
    rhs_id: "du/dt=[0,0]:v1",
    analytic_id: "u(t)=u0:exact:v1",
    t0: 0.0,
    t_end: 2.0,
    y0: &STATIONARY_Y0,
    ceiling: 0.0,
};

const DECAY: FixtureSpec = FixtureSpec {
    kind: FixtureKind::Decay,
    id: "scalar-decay-half-rate",
    rhs_id: "du/dt=-0.5*u:v1",
    analytic_id: "u(t)=u0*exp(-0.5*t):std-exp:v1",
    t0: 0.0,
    t_end: 4.0,
    y0: &DECAY_Y0,
    ceiling: DECAY_CEILING,
};

const OSCILLATOR: FixtureSpec = FixtureSpec {
    kind: FixtureKind::Oscillator,
    id: "two-component-harmonic-oscillator",
    rhs_id: "dq/dt=v;dv/dt=-q:v1",
    analytic_id: "q=q0*cos(t)+v0*sin(t);v=-q0*sin(t)+v0*cos(t):std-sin-cos:v1",
    t0: 0.0,
    t_end: core::f64::consts::TAU,
    y0: &OSCILLATOR_Y0,
    ceiling: OSCILLATOR_CEILING,
};

impl FixtureSpec {
    fn rhs_into(self, _t: f64, y: &[f64], out: &mut [f64]) {
        match self.kind {
            FixtureKind::Stationary => out.fill(0.0),
            FixtureKind::Decay => out[0] = -0.5 * y[0],
            FixtureKind::Oscillator => {
                out[0] = y[1];
                out[1] = -y[0];
            }
        }
    }

    fn rhs_vec(self, t: f64, y: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; y.len()];
        self.rhs_into(t, y, &mut out);
        out
    }

    fn analytic_endpoint(self) -> Vec<f64> {
        let dt = self.t_end - self.t0;
        match self.kind {
            FixtureKind::Stationary => self.y0.to_vec(),
            FixtureKind::Decay => vec![self.y0[0] * (-0.5 * dt).exp()],
            FixtureKind::Oscillator => {
                let (sin, cos) = dt.sin_cos();
                let q0 = self.y0[0];
                let v0 = self.y0[1];
                vec![q0.mul_add(cos, v0 * sin), (-q0).mul_add(sin, v0 * cos)]
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ProductionBits {
    t: u64,
    u: Vec<u64>,
    h: u64,
    err_prev: u64,
    accepted: usize,
    rejected: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct DenseSolutionBits {
    knots: Vec<u64>,
    values: Vec<Vec<u64>>,
    alt_segment: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleBits {
    t: Vec<u64>,
    y: Vec<Vec<u64>>,
    sol: Option<DenseSolutionBits>,
    t_events: Option<Vec<Vec<u64>>>,
    y_events: Option<Vec<Vec<Vec<u64>>>>,
    nfev: usize,
    njev: usize,
    nlu: usize,
    status: i32,
    message: String,
    success: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Measurement {
    analytic: Vec<u64>,
    production: ProductionBits,
    oracle: OracleBits,
}

#[derive(Debug, Clone)]
struct ErrorEvidence {
    production_minus_analytic: Vec<f64>,
    oracle_minus_analytic: Vec<f64>,
    production_minus_oracle: Vec<f64>,
    production_analytic_max: f64,
    oracle_analytic_max: f64,
    production_oracle_max: f64,
}

#[derive(Debug, Clone)]
struct Corruption {
    component: usize,
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
        u64::try_from(value).expect("ODE Casebook frame lengths fit u64"),
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

fn push_i32_field(bytes: &mut Vec<u8>, tag: &str, value: i32) {
    push_field(bytes, tag, &value.to_le_bytes());
}

fn push_u64_field(bytes: &mut Vec<u8>, tag: &str, value: u64) {
    push_field(bytes, tag, &value.to_le_bytes());
}

fn push_usize_field(bytes: &mut Vec<u8>, tag: &str, value: usize) {
    push_u64_field(
        bytes,
        tag,
        u64::try_from(value).expect("ODE Casebook values fit u64"),
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
    let value_bits = bits(values);
    push_bits_field(bytes, tag, &value_bits);
}

fn push_bit_rows(bytes: &mut Vec<u8>, tag: &str, rows: &[Vec<u64>]) {
    let mut payload = Vec::new();
    push_len(&mut payload, rows.len());
    for row in rows {
        push_bits_field(&mut payload, "row", row);
    }
    push_field(bytes, tag, &payload);
}

fn push_bit_tables(bytes: &mut Vec<u8>, tag: &str, tables: &[Vec<Vec<u64>>]) {
    let mut payload = Vec::new();
    push_len(&mut payload, tables.len());
    for table in tables {
        push_bit_rows(&mut payload, "table", table);
    }
    push_field(bytes, tag, &payload);
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn json_string_field<'a>(object: &'a str, field: &str) -> Option<&'a str> {
    let field = format!("\"{field}\"");
    let (_, tail) = object.split_once(field.as_str())?;
    let (_, value) = tail.split_once(':')?;
    let value = value.trim_start().strip_prefix('"')?;
    value.split_once('"').map(|(value, _)| value)
}

fn constellation_library_object(library: &str) -> Option<&'static str> {
    CONSTELLATION_LOCK.split('}').find_map(|prefix| {
        let object = prefix.rsplit_once('{').map_or(prefix, |(_, object)| object);
        (json_string_field(object, "lib") == Some(library)).then_some(object)
    })
}

fn admit_oracle_declaration() -> Result<(), String> {
    let object = constellation_library_object("frankenscipy").ok_or_else(|| {
        "stage=oracle-pin-declaration; constellation.lock has no frankenscipy object".to_owned()
    })?;
    let framed_version = ORACLE_VERSION.strip_prefix("fsci-integrate/");
    let declared_version = json_string_field(object, "version");
    let declared_pin = json_string_field(object, "git_head");
    if framed_version != Some(ORACLE_LOCK_VERSION)
        || declared_version != Some(ORACLE_LOCK_VERSION)
        || declared_pin != Some(ORACLE_PIN)
    {
        return Err(format!(
            "stage=oracle-pin-declaration; expected_version={ORACLE_LOCK_VERSION}; framed_version={framed_version:?}; declared_version={declared_version:?}; expected_pin={ORACLE_PIN}; declared_pin={declared_pin:?}"
        ));
    }
    Ok(())
}

fn bit_rows(rows: &[Vec<f64>]) -> Vec<Vec<u64>> {
    rows.iter().map(|row| bits(row)).collect()
}

fn bit_tables(tables: &[Vec<Vec<f64>>]) -> Vec<Vec<Vec<u64>>> {
    tables.iter().map(|table| bit_rows(table)).collect()
}

fn common_frame(domain: &str) -> Vec<u8> {
    let pi = PiController::default();
    let mut bytes = Vec::new();
    push_text_field(&mut bytes, "domain", domain);
    push_text_field(&mut bytes, "encoding", FRAME_ENCODING);
    push_u32_field(
        &mut bytes,
        "casebook-record-version",
        CASEBOOK_RECORD_VERSION,
    );
    push_text_field(&mut bytes, "fs-time-version", FS_TIME_VERSION);
    push_text_field(&mut bytes, "production-api", PRODUCTION_API);
    push_text_field(&mut bytes, "oracle-version", ORACLE_VERSION);
    push_text_field(&mut bytes, "oracle-pin", ORACLE_PIN);
    push_text_field(&mut bytes, "oracle-api", ORACLE_API);
    push_text_field(&mut bytes, "oracle-runtime-mode", "strict:v1");
    push_text_field(&mut bytes, "unit-policy", UNIT_POLICY);
    push_f64_field(&mut bytes, "relative-tolerance", RTOL);
    push_f64_field(&mut bytes, "absolute-tolerance", ATOL);
    push_f64_field(&mut bytes, "production-initial-step", INITIAL_STEP);
    push_usize_field(&mut bytes, "production-attempt-cap", MAX_STEPS);
    push_f64_field(&mut bytes, "pi-k-p", pi.k_p);
    push_f64_field(&mut bytes, "pi-k-i", pi.k_i);
    push_f64_field(&mut bytes, "pi-safety", pi.safety);
    push_f64_field(&mut bytes, "pi-max-growth", pi.max_growth);
    push_f64_field(&mut bytes, "pi-max-shrink", pi.max_shrink);
    push_text_field(&mut bytes, "oracle-method", "RK45");
    push_text_field(&mut bytes, "oracle-absolute-tolerance-kind", "scalar:v1");
    push_f64_field(&mut bytes, "oracle-first-step", INITIAL_STEP);
    push_f64_field(&mut bytes, "oracle-max-step", ORACLE_MAX_STEP);
    push_bool_field(&mut bytes, "oracle-dense-output", false);
    push_bool_field(&mut bytes, "oracle-events-present", false);
    push_text_field(&mut bytes, "error-policy", ERROR_POLICY);
    push_text_field(&mut bytes, "analytic-policy", ANALYTIC_POLICY);
    push_text_field(&mut bytes, "no-claim-policy", NO_CLAIM_POLICY);
    bytes
}

fn fixture_inputs(fixture: FixtureSpec) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-time:frankenscipy-rk45-endpoint-input:v1");
    push_text_field(&mut bytes, "fixture-id", fixture.id);
    push_text_field(&mut bytes, "rhs-identity", fixture.rhs_id);
    push_text_field(&mut bytes, "analytic-identity", fixture.analytic_id);
    push_usize_field(&mut bytes, "state-dimension", fixture.y0.len());
    push_f64_field(&mut bytes, "initial-time", fixture.t0);
    push_f64_field(&mut bytes, "end-time", fixture.t_end);
    push_f64s_field(&mut bytes, "initial-state", fixture.y0);
    push_f64s_field(
        &mut bytes,
        "oracle-evaluation-times",
        &[fixture.t0, fixture.t_end],
    );
    push_f64_field(&mut bytes, "max-norm-absolute-ceiling", fixture.ceiling);
    bytes
}

fn production_bits(state: &AdaptiveState) -> ProductionBits {
    ProductionBits {
        t: state.t.to_bits(),
        u: bits(&state.u),
        h: state.h.to_bits(),
        err_prev: state.err_prev.to_bits(),
        accepted: state.accepted,
        rejected: state.rejected,
    }
}

fn oracle_bits(result: SolveIvpResult) -> OracleBits {
    let SolveIvpResult {
        t,
        y,
        sol,
        t_events,
        y_events,
        nfev,
        njev,
        nlu,
        status,
        message,
        success,
    } = result;
    OracleBits {
        t: bits(&t),
        y: bit_rows(&y),
        sol: sol.map(|dense| DenseSolutionBits {
            knots: bits(&dense.knots),
            values: bit_rows(&dense.values),
            alt_segment: dense.alt_segment,
        }),
        t_events: t_events.map(|rows| bit_rows(&rows)),
        y_events: y_events.map(|tables| bit_tables(&tables)),
        nfev,
        njev,
        nlu,
        status,
        message,
        success,
    }
}

fn measure_fixture(fixture: FixtureSpec) -> Result<Measurement, String> {
    let pi = PiController::default();
    let mut state = AdaptiveState::new(fixture.t0, fixture.y0, INITIAL_STEP);
    let production_rhs = |t: f64, y: &[f64], out: &mut [f64]| fixture.rhs_into(t, y, out);
    rk45_adaptive(
        &mut state,
        &production_rhs,
        fixture.t_end,
        RTOL,
        ATOL,
        &pi,
        MAX_STEPS,
    );

    let t_eval = [fixture.t0, fixture.t_end];
    let mut oracle_rhs = |t: f64, y: &[f64]| fixture.rhs_vec(t, y);
    let oracle = solve_ivp(
        &mut oracle_rhs,
        &SolveIvpOptions {
            t_span: (fixture.t0, fixture.t_end),
            y0: fixture.y0,
            method: SolverKind::Rk45,
            t_eval: Some(&t_eval),
            dense_output: false,
            events: None,
            rtol: RTOL,
            atol: ToleranceValue::Scalar(ATOL),
            first_step: Some(INITIAL_STEP),
            max_step: ORACLE_MAX_STEP,
            ..SolveIvpOptions::default()
        },
    )
    .map_err(|error| {
        format!(
            "fixture={}; oracle-refusal={error}; pin={ORACLE_PIN}; t_span_bits=[0x{:016x},0x{:016x}]",
            fixture.id,
            fixture.t0.to_bits(),
            fixture.t_end.to_bits(),
        )
    })?;

    Ok(Measurement {
        analytic: bits(&fixture.analytic_endpoint()),
        production: production_bits(&state),
        oracle: oracle_bits(oracle),
    })
}

fn finite_bits(values: &[u64]) -> bool {
    values
        .iter()
        .all(|&value| f64::from_bits(value).is_finite())
}

fn signed_delta(left: &[u64], right: &[u64]) -> Vec<f64> {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| f64::from_bits(left) - f64::from_bits(right))
        .collect()
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().map(|value| value.abs()).fold(0.0, f64::max)
}

fn error_evidence(measurement: &Measurement) -> ErrorEvidence {
    let oracle_endpoint = measurement
        .oracle
        .y
        .last()
        .expect("oracle shape is admitted before error evidence");
    let production_minus_analytic = signed_delta(&measurement.production.u, &measurement.analytic);
    let oracle_minus_analytic = signed_delta(oracle_endpoint, &measurement.analytic);
    let production_minus_oracle = signed_delta(&measurement.production.u, oracle_endpoint);
    ErrorEvidence {
        production_analytic_max: max_abs(&production_minus_analytic),
        oracle_analytic_max: max_abs(&oracle_minus_analytic),
        production_oracle_max: max_abs(&production_minus_oracle),
        production_minus_analytic,
        oracle_minus_analytic,
        production_minus_oracle,
    }
}

#[allow(clippy::too_many_lines)] // One ordered gate keeps every public field fail-closed.
fn admit_measurement(
    fixture: FixtureSpec,
    measurement: &Measurement,
) -> Result<ErrorEvidence, String> {
    let dimension = fixture.y0.len();
    if measurement.production.u.len() != dimension {
        return Err(format!(
            "stage=production-shape; fixture={}; expected={dimension}; found={}",
            fixture.id,
            measurement.production.u.len(),
        ));
    }
    if !finite_bits(&[
        measurement.production.t,
        measurement.production.h,
        measurement.production.err_prev,
    ]) || !finite_bits(&measurement.production.u)
    {
        return Err(format!(
            "stage=production-finite; fixture={}; state={:?}",
            fixture.id, measurement.production,
        ));
    }
    if f64::from_bits(measurement.production.h) <= 0.0
        || f64::from_bits(measurement.production.err_prev) <= 0.0
    {
        return Err(format!(
            "stage=production-controller-state; fixture={}; next_step_bits=0x{:016x}; previous_error_bits=0x{:016x}",
            fixture.id, measurement.production.h, measurement.production.err_prev,
        ));
    }
    let Some(attempts) = measurement
        .production
        .accepted
        .checked_add(measurement.production.rejected)
    else {
        return Err(format!(
            "stage=production-attempt-overflow; fixture={}",
            fixture.id,
        ));
    };
    if measurement.production.accepted == 0 || attempts == 0 || attempts > MAX_STEPS {
        return Err(format!(
            "stage=production-work-admission; fixture={}; accepted={}; rejected={}; cap={MAX_STEPS}",
            fixture.id, measurement.production.accepted, measurement.production.rejected,
        ));
    }
    if measurement.production.t != fixture.t_end.to_bits() {
        return Err(format!(
            "stage=production-incomplete; fixture={}; reached_bits=0x{:016x}; expected_bits=0x{:016x}; accepted={}; rejected={}",
            fixture.id,
            measurement.production.t,
            fixture.t_end.to_bits(),
            measurement.production.accepted,
            measurement.production.rejected,
        ));
    }

    let expected_times = bits(&[fixture.t0, fixture.t_end]);
    if measurement.oracle.t != expected_times
        || measurement.oracle.y.len() != expected_times.len()
        || measurement
            .oracle
            .y
            .iter()
            .any(|row| row.len() != dimension)
    {
        return Err(format!(
            "stage=oracle-shape; fixture={}; times={:016x?}; expected_times={expected_times:016x?}; row_lengths={:?}; dimension={dimension}",
            fixture.id,
            measurement.oracle.t,
            measurement
                .oracle
                .y
                .iter()
                .map(Vec::len)
                .collect::<Vec<_>>(),
        ));
    }
    if !finite_bits(&measurement.oracle.t)
        || measurement.oracle.y.iter().any(|row| !finite_bits(row))
    {
        return Err(format!(
            "stage=oracle-finite; fixture={}; oracle={:?}",
            fixture.id, measurement.oracle,
        ));
    }
    if !measurement.oracle.success
        || measurement.oracle.status != 0
        || measurement.oracle.nfev == 0
        || measurement.oracle.njev != 0
        || measurement.oracle.nlu != 0
        || measurement.oracle.sol.is_some()
        || measurement.oracle.t_events.is_some()
        || measurement.oracle.y_events.is_some()
    {
        return Err(format!(
            "stage=oracle-admission; fixture={}; success={}; status={}; nfev={}; njev={}; nlu={}; sol_present={}; t_events_present={}; y_events_present={}; message={:?}",
            fixture.id,
            measurement.oracle.success,
            measurement.oracle.status,
            measurement.oracle.nfev,
            measurement.oracle.njev,
            measurement.oracle.nlu,
            measurement.oracle.sol.is_some(),
            measurement.oracle.t_events.is_some(),
            measurement.oracle.y_events.is_some(),
            measurement.oracle.message,
        ));
    }
    if measurement.oracle.y[0] != bits(fixture.y0) {
        return Err(format!(
            "stage=oracle-initial-state; fixture={}; observed={:016x?}; expected={:016x?}",
            fixture.id,
            measurement.oracle.y[0],
            bits(fixture.y0),
        ));
    }
    if measurement.analytic.len() != dimension || !finite_bits(&measurement.analytic) {
        return Err(format!(
            "stage=analytic-admission; fixture={}; analytic={:016x?}; dimension={dimension}",
            fixture.id, measurement.analytic,
        ));
    }

    let evidence = error_evidence(measurement);
    let maxima = [
        evidence.production_analytic_max,
        evidence.oracle_analytic_max,
        evidence.production_oracle_max,
    ];
    if maxima.iter().any(|value| !value.is_finite()) {
        return Err(format!(
            "stage=error-finite; fixture={}; maxima={maxima:?}",
            fixture.id,
        ));
    }
    if fixture.kind == FixtureKind::Stationary {
        let expected = bits(fixture.y0);
        if measurement.analytic != expected
            || measurement.production.u != expected
            || measurement.oracle.y.iter().any(|row| row != &expected)
            || maxima
                .iter()
                .any(|&value| value.to_bits() != 0.0_f64.to_bits())
        {
            return Err(format!(
                "stage=stationary-exact-kat; analytic={:016x?}; production={:016x?}; oracle={:016x?}; expected={expected:016x?}; maxima={maxima:?}",
                measurement.analytic, measurement.production.u, measurement.oracle.y,
            ));
        }
    } else if maxima.iter().any(|&value| value > fixture.ceiling) {
        return Err(format!(
            "stage=bounded-endpoint-agreement; fixture={}; production_bits={:016x?}; oracle_bits={:016x?}; analytic_bits={:016x?}; production_analytic={}; oracle_analytic={}; production_oracle={}; ceiling={}; margins=[{},{},{}]",
            fixture.id,
            measurement.production.u,
            measurement.oracle.y.last().expect("oracle shape admitted"),
            measurement.analytic,
            evidence.production_analytic_max,
            evidence.oracle_analytic_max,
            evidence.production_oracle_max,
            fixture.ceiling,
            fixture.ceiling - evidence.production_analytic_max,
            fixture.ceiling - evidence.oracle_analytic_max,
            fixture.ceiling - evidence.production_oracle_max,
        ));
    }
    Ok(evidence)
}

fn push_production(bytes: &mut Vec<u8>, production: &ProductionBits) {
    push_u64_field(bytes, "time-bits", production.t);
    push_bits_field(bytes, "state-bits", &production.u);
    push_u64_field(bytes, "next-step-bits", production.h);
    push_u64_field(bytes, "previous-error-bits", production.err_prev);
    push_usize_field(bytes, "accepted-steps", production.accepted);
    push_usize_field(bytes, "rejected-steps", production.rejected);
}

fn push_oracle(bytes: &mut Vec<u8>, oracle: &OracleBits) {
    push_bits_field(bytes, "evaluation-time-bits", &oracle.t);
    push_bit_rows(bytes, "state-rows", &oracle.y);
    push_bool_field(bytes, "dense-solution-present", oracle.sol.is_some());
    if let Some(dense) = &oracle.sol {
        let mut payload = Vec::new();
        push_bits_field(&mut payload, "knot-bits", &dense.knots);
        push_bit_rows(&mut payload, "value-rows", &dense.values);
        push_bool_field(&mut payload, "alternate-segment", dense.alt_segment);
        push_field(bytes, "dense-solution", &payload);
    }
    push_bool_field(bytes, "time-events-present", oracle.t_events.is_some());
    if let Some(events) = &oracle.t_events {
        push_bit_rows(bytes, "time-events", events);
    }
    push_bool_field(bytes, "state-events-present", oracle.y_events.is_some());
    if let Some(events) = &oracle.y_events {
        push_bit_tables(bytes, "state-events", events);
    }
    push_usize_field(bytes, "function-evaluations", oracle.nfev);
    push_usize_field(bytes, "jacobian-evaluations", oracle.njev);
    push_usize_field(bytes, "lu-factorizations", oracle.nlu);
    push_i32_field(bytes, "status", oracle.status);
    push_text_field(bytes, "message", &oracle.message);
    push_bool_field(bytes, "success", oracle.success);
}

fn output_receipt(
    fixture: FixtureSpec,
    inputs: &[u8],
    measurement: &Measurement,
    evidence: &ErrorEvidence,
) -> Vec<u8> {
    let mut bytes = common_frame("bedrock:fs-time:frankenscipy-rk45-endpoint-output:v1");
    push_field(&mut bytes, "canonical-input-frame", inputs);
    push_text_field(&mut bytes, "fixture-id", fixture.id);
    push_bits_field(&mut bytes, "analytic-endpoint-bits", &measurement.analytic);
    let mut production = Vec::new();
    push_production(&mut production, &measurement.production);
    push_field(&mut bytes, "production-output", &production);
    let mut oracle = Vec::new();
    push_oracle(&mut oracle, &measurement.oracle);
    push_field(&mut bytes, "oracle-output", &oracle);
    push_f64s_field(
        &mut bytes,
        "production-minus-analytic",
        &evidence.production_minus_analytic,
    );
    push_f64s_field(
        &mut bytes,
        "oracle-minus-analytic",
        &evidence.oracle_minus_analytic,
    );
    push_f64s_field(
        &mut bytes,
        "production-minus-oracle",
        &evidence.production_minus_oracle,
    );
    for (tag, maximum) in [
        ("production-analytic-max", evidence.production_analytic_max),
        ("oracle-analytic-max", evidence.oracle_analytic_max),
        ("production-oracle-max", evidence.production_oracle_max),
    ] {
        push_f64_field(&mut bytes, tag, maximum);
        push_f64_field(
            &mut bytes,
            &format!("{tag}-remaining-margin"),
            fixture.ceiling - maximum,
        );
    }
    push_f64_field(&mut bytes, "max-norm-absolute-ceiling", fixture.ceiling);
    bytes
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
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

fn capture_measurement(fixture: FixtureSpec, stage: &str) -> Result<Measurement, String> {
    match catch_unwind(AssertUnwindSafe(|| measure_fixture(fixture))) {
        Ok(result) => result,
        Err(payload) => Err(format!(
            "stage={stage}; fixture={}; panic={}",
            fixture.id,
            panic_message(&*payload),
        )),
    }
}

fn fixture_outcome(fixture: FixtureSpec) -> CaseOutcome {
    if let Err(error) = admit_oracle_declaration() {
        return CaseOutcome::fail(error).with_evidence("constellation.lock:frankenscipy-0.2.0");
    }
    let inputs = fixture_inputs(fixture);
    let first = match capture_measurement(fixture, "first-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay = match capture_measurement(fixture, "replay-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    let first_evidence = match admit_measurement(fixture, &first) {
        Ok(evidence) => evidence,
        Err(error) => return CaseOutcome::fail(error),
    };
    let replay_evidence = match admit_measurement(fixture, &replay) {
        Ok(evidence) => evidence,
        Err(error) => return CaseOutcome::fail(format!("stage=replay-admission; {error}")),
    };
    let first_receipt = output_receipt(fixture, &inputs, &first, &first_evidence);
    let replay_receipt = output_receipt(fixture, &inputs, &replay, &replay_evidence);
    if first_receipt != replay_receipt {
        return CaseOutcome::fail(format!(
            "stage=same-run-output-replay; fixture={}; first_len={}; replay_len={}; first_fnv1a64=0x{:016x}; replay_fnv1a64=0x{:016x}; first={}; replay={}",
            fixture.id,
            first_receipt.len(),
            replay_receipt.len(),
            fnv1a64(&first_receipt),
            fnv1a64(&replay_receipt),
            hex_bytes(&first_receipt),
            hex_bytes(&replay_receipt),
        ))
        .with_evidence("crates/fs-time/CONTRACT.md#determinism-class")
        .with_evidence("constellation.lock:frankenscipy-0.2.0");
    }

    CaseOutcome::pass(format!(
        "fixture={}; production_analytic={:.3e}; oracle_analytic={:.3e}; production_oracle={:.3e}; ceiling={:.3e}; production_steps={}/{}; oracle_nfev={}; output_receipt_len={}; output_receipt_fnv1a64=0x{:016x}; output_receipt={}",
        fixture.id,
        first_evidence.production_analytic_max,
        first_evidence.oracle_analytic_max,
        first_evidence.production_oracle_max,
        fixture.ceiling,
        first.production.accepted,
        first.production.rejected,
        first.oracle.nfev,
        first_receipt.len(),
        fnv1a64(&first_receipt),
        hex_bytes(&first_receipt),
    ))
    .with_evidence("crates/fs-time/CONTRACT.md#conformance-tests")
    .with_evidence("constellation.lock:frankenscipy-0.2.0")
}

fn corruption_frame(component: usize, bit: u32, canonical: u64, corrupted: u64) -> Vec<u8> {
    let stationary = fixture_inputs(STATIONARY);
    let mut bytes = common_frame("bedrock:fs-time:stationary-oracle-reference-corruption:v1");
    push_field(&mut bytes, "canonical-stationary-input-frame", &stationary);
    push_u64_field(
        &mut bytes,
        "canonical-stationary-input-fnv1a64",
        fnv1a64(&stationary),
    );
    push_u64_field(&mut bytes, "corruption-seed", CORRUPTION_SEED);
    push_text_field(
        &mut bytes,
        "selection-policy",
        "component=(seed>>8)%dimension;bit=(seed&0xff)%52:v1",
    );
    push_text_field(&mut bytes, "reference", "pinned-oracle-stationary-endpoint");
    push_usize_field(&mut bytes, "component", component);
    push_u32_field(&mut bytes, "mantissa-bit", bit);
    push_u64_field(&mut bytes, "canonical-value-bits", canonical);
    push_u64_field(&mut bytes, "corrupted-value-bits", corrupted);
    bytes
}

fn reconstruct_corruption() -> Corruption {
    let component = usize::try_from((CORRUPTION_SEED >> 8) % STATIONARY.y0.len() as u64)
        .expect("corruption component fits usize");
    let bit = u32::try_from((CORRUPTION_SEED & 0xff) % 52).expect("corruption bit fits u32");
    let canonical = STATIONARY.y0[component].to_bits();
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
    if let Err(error) = admit_oracle_declaration() {
        return CaseOutcome::fail(error).with_evidence("constellation.lock:frankenscipy-0.2.0");
    }
    let measurement = match capture_measurement(STATIONARY, "red-baseline-measurement") {
        Ok(measurement) => measurement,
        Err(error) => return CaseOutcome::fail(error),
    };
    if let Err(error) = admit_measurement(STATIONARY, &measurement) {
        return CaseOutcome::fail(format!("stage=red-baseline-admission; {error}"));
    }
    let computed = measurement
        .oracle
        .y
        .last()
        .expect("stationary oracle shape admitted")[corruption.component];
    if computed != corruption.canonical {
        return CaseOutcome::fail(format!(
            "stage=seeded-corruption-baseline-drift; seed=0x{CORRUPTION_SEED:016x}; component={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}",
            corruption.component, corruption.canonical,
        ));
    }

    let mut corrupted_measurement = measurement;
    corrupted_measurement
        .oracle
        .y
        .last_mut()
        .expect("stationary oracle shape admitted")[corruption.component] = corruption.corrupted;
    match admit_measurement(STATIONARY, &corrupted_measurement) {
        Ok(_) => CaseOutcome::pass(format!(
            "stage=seeded-corruption-not-detected; seed=0x{CORRUPTION_SEED:016x}; component={}; bit={}; canonical_bits=0x{:016x}; corrupted_bits=0x{:016x}",
            corruption.component,
            corruption.bit,
            corruption.canonical,
            corruption.corrupted,
        )),
        Err(refusal) if refusal.starts_with("stage=stationary-exact-kat;") => CaseOutcome::fail(format!(
            "stage=seeded-stationary-oracle-reference-corruption; seed=0x{CORRUPTION_SEED:016x}; component={}; bit={}; computed_bits=0x{computed:016x}; canonical_bits=0x{:016x}; corrupted_bits=0x{:016x}; refusal={refusal}; input_frame_len={}; input_frame_fnv1a64=0x{:016x}; input_frame={}",
            corruption.component,
            corruption.bit,
            corruption.canonical,
            corruption.corrupted,
            corruption.frame.len(),
            fnv1a64(&corruption.frame),
            hex_bytes(&corruption.frame),
        ))
        .with_evidence("crates/fs-time/tests/frankenscipy_ode_oracle_casebook.rs#seeded-corruption"),
        Err(refusal) => CaseOutcome::fail(format!(
            "stage=seeded-corruption-unexpected-refusal; seed=0x{CORRUPTION_SEED:016x}; component={}; bit={}; refusal={refusal}",
            corruption.component, corruption.bit,
        )),
    }
}

fn run_red_report() -> SuiteReport {
    let corruption = reconstruct_corruption();
    assert_eq!(corruption.frame.len(), CORRUPTION_FRAME_LEN);
    assert_eq!(fnv1a64(&corruption.frame), CORRUPTION_FRAME_FNV1A64);
    Suite::new(SUITE)
        .case(
            "seeded-stationary-oracle-reference-bit-corruption",
            CORRUPTION_FRAME_FNV1A64,
            ToleranceSpec::Exact,
            move || corruption_outcome(corruption),
        )
        .run()
}

#[test]
fn frankenscipy_rk45_oracle_casebook_emits_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let stationary = fixture_inputs(STATIONARY);
    let decay = fixture_inputs(DECAY);
    let oscillator = fixture_inputs(OSCILLATOR);
    assert_eq!(stationary, fixture_inputs(STATIONARY));
    assert_eq!(decay, fixture_inputs(DECAY));
    assert_eq!(oscillator, fixture_inputs(OSCILLATOR));
    assert_eq!(stationary.len(), STATIONARY_FRAME_LEN);
    assert_eq!(fnv1a64(&stationary), STATIONARY_FRAME_FNV1A64);
    assert_eq!(decay.len(), DECAY_FRAME_LEN);
    assert_eq!(fnv1a64(&decay), DECAY_FRAME_FNV1A64);
    assert_eq!(oscillator.len(), OSCILLATOR_FRAME_LEN);
    assert_eq!(fnv1a64(&oscillator), OSCILLATOR_FRAME_FNV1A64);

    let report = Suite::new(SUITE)
        .case(
            "stationary-two-component-exact-endpoint",
            STATIONARY_FRAME_FNV1A64,
            ToleranceSpec::Exact,
            || fixture_outcome(STATIONARY),
        )
        .case(
            "scalar-decay-pinned-rk45-endpoint",
            DECAY_FRAME_FNV1A64,
            ToleranceSpec::AbsoluteLe(DECAY_CEILING),
            || fixture_outcome(DECAY),
        )
        .case(
            "harmonic-oscillator-pinned-rk45-endpoint",
            OSCILLATOR_FRAME_FNV1A64,
            ToleranceSpec::AbsoluteLe(OSCILLATOR_CEILING),
            || fixture_outcome(OSCILLATOR),
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
            "stationary-two-component-exact-endpoint",
            "scalar-decay-pinned-rk45-endpoint",
            "harmonic-oscillator-pinned-rk45-endpoint",
        ]
    );
    assert!(report.records.iter().all(|record| {
        record.version == CASEBOOK_RECORD_VERSION
            && record.pass
            && !record.evidence.is_empty()
            && record.details.contains("output_receipt=")
    }));
    assert_eq!(report.records[0].tolerance, "exact");
    assert_eq!(report.records[1].tolerance, "abs<=2e-7");
    assert_eq!(report.records[2].tolerance, "abs<=5e-7");
}

#[test]
fn seeded_stationary_oracle_reference_corruption_is_stable_and_refused() {
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
            .contains("stage=seeded-stationary-oracle-reference-corruption")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(first_failure.details.contains("component=1"));
    assert!(first_failure.details.contains("bit=1"));
    assert!(
        first_failure
            .details
            .contains("refusal=stage=stationary-exact-kat")
    );
    assert!(first_failure.details.contains("input_frame="));

    let panic = catch_unwind(|| first.assert_green())
        .expect_err("assert_green must refuse the disclosed oracle-reference corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook refusal carries text");
    assert!(message.contains("seeded-stationary-oracle-reference-bit-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
