//! G3/G5 replay receipt for the production greedy batched q-EI path (7tv.21.35).
//!
//! A short `q = 2` normalized-Branin study retains every ordered objective
//! callback, every public `BoReport` bit, and exact batch boundaries. The
//! fixture reconstructs its Sobol initialization, checks objective values with
//! an algebraically independent reference, and recomputes best-trace bits.
//! Public-report prefix runs witness every batch boundary, and an independent
//! run must reproduce the complete frame. A disclosed acquired-point mutation
//! must fail stale-payload,
//! retained-reference, semantic, wire-evidence, and test-local merge gates.
//!
//! This target makes no optimizer-quality, convergence, analytic-qEI,
//! all-objective, all-configuration, all-seed, internal-acquisition-search,
//! cross-ISA, cancellation, persistence, authenticated-admission,
//! external-oracle, or performance claim.

#![deny(unsafe_code)]

use fs_bo::{BoConfig, BoReport, Matern, minimize};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_rand::StreamKey;

const SUITE: &str = "fs-bo/qei-study-replay-v1";
const STUDY_CASE: &str = "greedy-q-two-branin-short";
const STUDY_SEED: u64 = 17;

const DIMENSION: usize = 2;
const N_INIT: usize = 4;
const ITERATIONS: usize = 2;
const BATCH_SIZE: usize = 2;
const EXPECTED_EVALUATIONS: usize = N_INIT + ITERATIONS * BATCH_SIZE;
const EXPECTED_TRACE_POINTS: usize = ITERATIONS + 1;
const HYPER_STARTS: usize = 1;
const ACQUISITION_STARTS: usize = 1;
const ACQUISITION_EVALUATIONS: usize = 60;
const MC_SAMPLES: usize = 32;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND: f64 = 1.0;
const LOG_BOX_LOWER: f64 = -2.0;
const LOG_BOX_UPPER: f64 = 0.5;
const QEI_BANK_SEED_XOR: u64 = 0xACC5;
const ACQUISITION_SOBOL_SEED_XOR: u64 = 0x5EED;
const BRANIN_ORACLE_ABS_TOLERANCE: f64 = 1.0e-11;
const BRANIN_ORACLE_REL_TOLERANCE: f64 = 1.0e-11;

const MUTATION_SEED: u64 = 0xB0A7_5EED_0000_0035;
const MUTATION_KERNEL: u32 = 0xB035;
const MUTATION_TILE: u32 = 0;

const _: () = assert!(DIMENSION == 2);
const _: () = assert!(BATCH_SIZE > 1);
const _: () = assert!(EXPECTED_EVALUATIONS == 8);
const _: () = assert!(EXPECTED_TRACE_POINTS == 3);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    point: Vec<u64>,
    value: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BatchRange {
    iteration: usize,
    start: usize,
    end: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PrefixReceipt {
    iterations: usize,
    callbacks: Vec<EvaluationBits>,
    report_points: Vec<Vec<u64>>,
    report_values: Vec<u64>,
    best_trace: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    callbacks: Vec<EvaluationBits>,
    report_points: Vec<Vec<u64>>,
    report_values: Vec<u64>,
    best_trace: Vec<u64>,
    batches: Vec<BatchRange>,
    prefixes: Vec<PrefixReceipt>,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-bo-qei-study-output-frame-v1".to_vec();
        push_u64(&mut bytes, config_digest);
        push_len(&mut bytes, self.callbacks.len());
        for callback in &self.callbacks {
            push_u64_slice(&mut bytes, &callback.point);
            push_u64(&mut bytes, callback.value);
        }
        push_len(&mut bytes, self.report_points.len());
        for point in &self.report_points {
            push_u64_slice(&mut bytes, point);
        }
        push_u64_slice(&mut bytes, &self.report_values);
        push_u64_slice(&mut bytes, &self.best_trace);
        push_len(&mut bytes, self.batches.len());
        for batch in &self.batches {
            push_len(&mut bytes, batch.iteration);
            push_len(&mut bytes, batch.start);
            push_len(&mut bytes, batch.end);
        }
        push_len(&mut bytes, self.prefixes.len());
        for prefix in &self.prefixes {
            push_len(&mut bytes, prefix.iterations);
            push_len(&mut bytes, prefix.callbacks.len());
            for callback in &prefix.callbacks {
                push_u64_slice(&mut bytes, &callback.point);
                push_u64(&mut bytes, callback.value);
            }
            push_len(&mut bytes, prefix.report_points.len());
            for point in &prefix.report_points {
                push_u64_slice(&mut bytes, point);
            }
            push_u64_slice(&mut bytes, &prefix.report_values);
            push_u64_slice(&mut bytes, &prefix.best_trace);
        }
        bytes
    }

    #[allow(clippy::too_many_lines)] // Complete callback/report/batch accounting is the receipt.
    fn semantic_mismatch(&self) -> Option<String> {
        if self.callbacks.len() != EXPECTED_EVALUATIONS {
            return Some(format!(
                "callback-count:{}!=expected-{EXPECTED_EVALUATIONS}",
                self.callbacks.len()
            ));
        }
        if self.report_points.len() != EXPECTED_EVALUATIONS {
            return Some(format!(
                "report-point-count:{}!=expected-{EXPECTED_EVALUATIONS}",
                self.report_points.len()
            ));
        }
        if self.report_values.len() != EXPECTED_EVALUATIONS {
            return Some(format!(
                "report-value-count:{}!=expected-{EXPECTED_EVALUATIONS}",
                self.report_values.len()
            ));
        }
        if self.best_trace.len() != EXPECTED_TRACE_POINTS {
            return Some(format!(
                "best-trace-count:{}!=expected-{EXPECTED_TRACE_POINTS}",
                self.best_trace.len()
            ));
        }
        if self.batches.len() != ITERATIONS {
            return Some(format!(
                "batch-count:{}!=expected-{ITERATIONS}",
                self.batches.len()
            ));
        }
        if self.prefixes.len() != ITERATIONS + 1 {
            return Some(format!(
                "prefix-count:{}!=expected-{}",
                self.prefixes.len(),
                ITERATIONS + 1
            ));
        }

        for (iterations, prefix) in self.prefixes.iter().enumerate() {
            let expected_evaluations = N_INIT + iterations * BATCH_SIZE;
            let expected_trace_points = iterations + 1;
            if prefix.iterations != iterations {
                return Some(format!(
                    "prefix[{iterations}].iterations:{}!=expected-{iterations}",
                    prefix.iterations
                ));
            }
            if prefix.callbacks.len() != expected_evaluations
                || prefix.report_points.len() != expected_evaluations
                || prefix.report_values.len() != expected_evaluations
            {
                return Some(format!(
                    "prefix[{iterations}]-evaluation-counts:callbacks={};points={};values={};expected={expected_evaluations}",
                    prefix.callbacks.len(),
                    prefix.report_points.len(),
                    prefix.report_values.len(),
                ));
            }
            if prefix.best_trace.len() != expected_trace_points {
                return Some(format!(
                    "prefix[{iterations}]-trace-count:{}!=expected-{expected_trace_points}",
                    prefix.best_trace.len()
                ));
            }
        }

        for (iteration, batch) in self.batches.iter().enumerate() {
            let expected_start = self.prefixes[iteration].report_points.len();
            let expected_end = self.prefixes[iteration + 1].report_points.len();
            if *batch
                != (BatchRange {
                    iteration,
                    start: expected_start,
                    end: expected_end,
                })
            {
                return Some(format!(
                    "batch[{iteration}]:actual={batch:?};expected={expected_start}..{expected_end}"
                ));
            }
            if iteration > 0 && self.batches[iteration - 1].end != batch.start {
                return Some(format!("batch[{iteration}]-not-contiguous"));
            }
            if expected_end - expected_start != BATCH_SIZE {
                return Some(format!(
                    "batch[{iteration}]-observed-width:{}!=expected-{BATCH_SIZE}",
                    expected_end - expected_start
                ));
            }
        }

        for index in 0..EXPECTED_EVALUATIONS {
            let callback = &self.callbacks[index];
            let report_point = &self.report_points[index];
            if callback.point.len() != DIMENSION || report_point.len() != DIMENSION {
                return Some(format!(
                    "evaluation[{index}]-dimension:callback-{};report-{};expected-{DIMENSION}",
                    callback.point.len(),
                    report_point.len()
                ));
            }
            if callback.point != *report_point || callback.value != self.report_values[index] {
                return Some(format!("evaluation[{index}]-callback-report"));
            }
            let point = decode(report_point);
            if point
                .iter()
                .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
            {
                return Some(format!("evaluation[{index}]-point-outside-box"));
            }
            let reported = f64::from_bits(self.report_values[index]);
            let recomputed = branin_oracle(&point);
            let tolerance = BRANIN_ORACLE_ABS_TOLERANCE
                + BRANIN_ORACLE_REL_TOLERANCE * reported.abs().max(recomputed.abs());
            if !reported.is_finite()
                || !recomputed.is_finite()
                || (reported - recomputed).abs() > tolerance
            {
                return Some(format!(
                    "evaluation[{index}]-objective-oracle:recomputed={recomputed:.17e};reported={reported:.17e};tolerance={tolerance:.3e}"
                ));
            }
        }

        let sobol = fs_rand::qmc::Sobol::scrambled(DIMENSION, STUDY_SEED);
        let mut unit = vec![0.0; DIMENSION];
        for index in 0..N_INIT {
            sobol.point(
                u32::try_from(index + 1).expect("fixture initialization index fits u32"),
                &mut unit,
            );
            for (coordinate, value) in unit.iter().enumerate() {
                let expected = (UPPER_BOUND - LOWER_BOUND)
                    .mul_add(*value, LOWER_BOUND)
                    .to_bits();
                if self.report_points[index][coordinate] != expected {
                    return Some(format!(
                        "initializer[{index}].x[{coordinate}]:0x{:016x}!=0x{expected:016x}",
                        self.report_points[index][coordinate]
                    ));
                }
            }
        }

        let mut best = self.report_values[..N_INIT]
            .iter()
            .copied()
            .map(f64::from_bits)
            .fold(f64::INFINITY, f64::min);
        if best.to_bits() != self.best_trace[0] {
            return Some(format!(
                "best-trace[0]:recomputed=0x{:016x};reported=0x{:016x}",
                best.to_bits(),
                self.best_trace[0]
            ));
        }
        for batch in &self.batches {
            for &value in &self.report_values[batch.start..batch.end] {
                best = best.min(f64::from_bits(value));
            }
            if best.to_bits() != self.best_trace[batch.iteration + 1] {
                return Some(format!(
                    "best-trace[{}]:recomputed=0x{:016x};reported=0x{:016x}",
                    batch.iteration + 1,
                    best.to_bits(),
                    self.best_trace[batch.iteration + 1]
                ));
            }
        }
        for (iterations, prefix) in self.prefixes.iter().enumerate() {
            let expected_evaluations = N_INIT + iterations * BATCH_SIZE;
            let expected_trace_points = iterations + 1;
            if prefix.callbacks.as_slice() != &self.callbacks[..expected_evaluations]
                || prefix.report_points.as_slice() != &self.report_points[..expected_evaluations]
                || prefix.report_values.as_slice() != &self.report_values[..expected_evaluations]
                || prefix.best_trace.as_slice() != &self.best_trace[..expected_trace_points]
            {
                return Some(format!("prefix[{iterations}]-frame-mismatch"));
            }
        }
        None
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SealedStudy {
    config_digest: u64,
    output_digest: u64,
    record: StudyRecord,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
}

impl SealedStudy {
    fn seal(config_digest: u64, record: StudyRecord) -> Self {
        let output_digest = fnv1a64(&record.canonical_bytes(config_digest));
        Self {
            config_digest,
            output_digest,
            record,
        }
    }

    fn validate_payload(&self) -> Result<(), AdmissionError> {
        let computed = fnv1a64(&self.record.canonical_bytes(self.config_digest));
        if computed == self.output_digest {
            Ok(())
        } else {
            Err(AdmissionError::PayloadIdentityMismatch {
                declared: self.output_digest,
                computed,
            })
        }
    }

    fn admit_against(&self, reference_output_digest: u64) -> Result<(), AdmissionError> {
        self.validate_payload()?;
        if self.output_digest == reference_output_digest {
            Ok(())
        } else {
            Err(AdmissionError::ReferenceIdentityMismatch {
                expected: reference_output_digest,
                found: self.output_digest,
            })
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    evaluation: usize,
    coordinate: usize,
    mantissa_bit: u32,
    selector_draws: u64,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SeededCorruption {
    run: SealedStudy,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    first_mismatch: String,
    semantic_mismatch: String,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(bytes, usize_u64(value));
}

fn push_str(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_u64_slice(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for &value in values {
        push_u64(bytes, value);
    }
}

fn config() -> BoConfig {
    BoConfig {
        bounds: (LOWER_BOUND, UPPER_BOUND),
        family: Matern::FiveHalves,
        log_box: (LOG_BOX_LOWER, LOG_BOX_UPPER),
        hyper_starts: HYPER_STARTS,
        acq_starts: ACQUISITION_STARTS,
        acq_evals: ACQUISITION_EVALUATIONS,
        q: BATCH_SIZE,
        mc_samples: MC_SAMPLES,
        seed: STUDY_SEED,
    }
}

fn config_bytes() -> Vec<u8> {
    let config = config();
    let mut bytes = b"fs-bo-qei-study-config-v1".to_vec();
    push_str(&mut bytes, STUDY_CASE);
    push_str(&mut bytes, "dimensionless-normalized-coordinate");
    push_str(&mut bytes, "dimensionless-branin-value");
    push_str(&mut bytes, "branin-rescaled-from-standard-domain-v1");
    push_str(
        &mut bytes,
        match config.family {
            Matern::Half => "matern-half",
            Matern::ThreeHalves => "matern-three-halves",
            Matern::FiveHalves => "matern-five-halves",
        },
    );
    push_u64(&mut bytes, usize_u64(DIMENSION));
    push_u64(&mut bytes, usize_u64(N_INIT));
    push_u64(&mut bytes, usize_u64(ITERATIONS));
    push_u64(&mut bytes, config.bounds.0.to_bits());
    push_u64(&mut bytes, config.bounds.1.to_bits());
    push_u64(&mut bytes, config.log_box.0.to_bits());
    push_u64(&mut bytes, config.log_box.1.to_bits());
    push_u64(&mut bytes, usize_u64(config.hyper_starts));
    push_u64(&mut bytes, usize_u64(config.acq_starts));
    push_u64(&mut bytes, usize_u64(config.acq_evals));
    push_u64(&mut bytes, usize_u64(config.q));
    push_u64(&mut bytes, usize_u64(config.mc_samples));
    push_u64(&mut bytes, config.seed);
    push_u64(&mut bytes, QEI_BANK_SEED_XOR);
    push_u64(&mut bytes, ACQUISITION_SOBOL_SEED_XOR);
    push_u64(&mut bytes, BRANIN_ORACLE_ABS_TOLERANCE.to_bits());
    push_u64(&mut bytes, BRANIN_ORACLE_REL_TOLERANCE.to_bits());
    push_str(
        &mut bytes,
        "Sobol-init+standardize+QMC-hyperfit+greedy-joint-qEI-fixed-normal-bank+CMA-ES-v1",
    );
    push_str(
        &mut bytes,
        "public-report-prefix-witness-iterations-zero-through-two-v1",
    );
    push_str(
        &mut bytes,
        "independent-direct-arithmetic-std-cos-branin-oracle-v1",
    );
    push_str(&mut bytes, "synchronous-direct-test-no-Cx");
    push_str(&mut bytes, fs_bo::VERSION);
    push_str(&mut bytes, fs_ascent::VERSION);
    push_str(&mut bytes, fs_dfo::VERSION);
    push_str(&mut bytes, fs_la::VERSION);
    push_str(&mut bytes, fs_math::VERSION);
    push_str(&mut bytes, fs_rand::VERSION);
    push_u64(&mut bytes, u64::from(fs_rand::STREAM_SEMANTICS_VERSION));
    push_str(&mut bytes, fs_rand::STREAM_POSITION_IDENTITY_DOMAIN);
    push_u64(&mut bytes, u64::from(CASEBOOK_RECORD_VERSION));
    push_str(
        &mut bytes,
        "no-quality-convergence-analytic-qei-all-objective-all-config-all-seed-internal-search-cross-ISA-Cx-persistence-authentication-oracle-performance-claim",
    );
    bytes
}

fn branin_objective(point: &[f64]) -> f64 {
    let x1 = 15.0f64.mul_add(point[0], -5.0);
    let x2 = 15.0 * point[1];
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let inner = b.mul_add(-(x1 * x1), c.mul_add(x1, x2 - 6.0));
    inner * inner
        + 10.0 * (1.0 - 1.0 / (8.0 * core::f64::consts::PI)) * fs_math::det::cos(x1)
        + 10.0
}

fn branin_oracle(point: &[f64]) -> f64 {
    let x1 = point[0] * 15.0 - 5.0;
    let x2 = point[1] * 15.0;
    let pi = core::f64::consts::PI;
    let residual = x2 - (5.1 * x1 * x1) / (4.0 * pi * pi) + (5.0 * x1) / pi - 6.0;
    let oscillation = 10.0 * (1.0 - 1.0 / (8.0 * pi)) * x1.cos();
    residual * residual + oscillation + 10.0
}

fn evaluation_bits(point: &[f64], value: f64) -> EvaluationBits {
    EvaluationBits {
        point: point.iter().map(|value| value.to_bits()).collect(),
        value: value.to_bits(),
    }
}

fn decode(bits: &[u64]) -> Vec<f64> {
    bits.iter().copied().map(f64::from_bits).collect()
}

fn run_prefix(iterations: usize) -> PrefixReceipt {
    let expected_evaluations = N_INIT + iterations * BATCH_SIZE;
    let mut callbacks = Vec::with_capacity(expected_evaluations);
    let report = {
        let mut objective = |point: &[f64]| {
            let value = branin_objective(point);
            callbacks.push(evaluation_bits(point, value));
            value
        };
        minimize(&mut objective, DIMENSION, N_INIT, iterations, &config())
    };
    let BoReport { x, y, best_trace } = report;
    PrefixReceipt {
        iterations,
        callbacks,
        report_points: x
            .iter()
            .map(|point| point.iter().map(|value| value.to_bits()).collect())
            .collect(),
        report_values: y.iter().map(|value| value.to_bits()).collect(),
        best_trace: best_trace.iter().map(|value| value.to_bits()).collect(),
    }
}

fn run_study(config_digest: u64) -> SealedStudy {
    let full = run_prefix(ITERATIONS);
    let prefixes: Vec<_> = (0..=ITERATIONS).map(run_prefix).collect();
    let batches = prefixes
        .windows(2)
        .enumerate()
        .map(|(iteration, adjacent)| BatchRange {
            iteration,
            start: adjacent[0].report_points.len(),
            end: adjacent[1].report_points.len(),
        })
        .collect();
    let PrefixReceipt {
        iterations: _,
        callbacks,
        report_points,
        report_values,
        best_trace,
    } = full;
    let record = StudyRecord {
        callbacks,
        report_points,
        report_values,
        best_trace,
        batches,
        prefixes,
    };
    SealedStudy::seal(config_digest, record)
}

fn first_study_mismatch(left: &StudyRecord, right: &StudyRecord) -> Option<String> {
    if left.callbacks.len() != right.callbacks.len() {
        return Some(format!(
            "callbacks.length:{}!={}",
            left.callbacks.len(),
            right.callbacks.len()
        ));
    }
    for (index, (left, right)) in left.callbacks.iter().zip(&right.callbacks).enumerate() {
        if left != right {
            return Some(format!("callbacks[{index}]"));
        }
    }
    if left.report_points.len() != right.report_points.len() {
        return Some(format!(
            "report.x.length:{}!={}",
            left.report_points.len(),
            right.report_points.len()
        ));
    }
    for (index, (left, right)) in left
        .report_points
        .iter()
        .zip(&right.report_points)
        .enumerate()
    {
        if left.len() != right.len() {
            return Some(format!(
                "report.x[{index}].length:{}!={}",
                left.len(),
                right.len()
            ));
        }
        if let Some((coordinate, (left, right))) = left
            .iter()
            .zip(right)
            .enumerate()
            .find(|(_, (left, right))| left != right)
        {
            return Some(format!(
                "report.x[{index}][{coordinate}]:0x{left:016x}!=0x{right:016x}"
            ));
        }
    }
    if left.report_values.len() != right.report_values.len() {
        return Some(format!(
            "report.y.length:{}!={}",
            left.report_values.len(),
            right.report_values.len()
        ));
    }
    if let Some((index, (left, right))) = left
        .report_values
        .iter()
        .zip(&right.report_values)
        .enumerate()
        .find(|(_, (left, right))| left != right)
    {
        return Some(format!("report.y[{index}]:0x{left:016x}!=0x{right:016x}"));
    }
    if left.best_trace.len() != right.best_trace.len() {
        return Some(format!(
            "report.best_trace.length:{}!={}",
            left.best_trace.len(),
            right.best_trace.len()
        ));
    }
    if let Some((index, (left, right))) = left
        .best_trace
        .iter()
        .zip(&right.best_trace)
        .enumerate()
        .find(|(_, (left, right))| left != right)
    {
        return Some(format!(
            "report.best_trace[{index}]:0x{left:016x}!=0x{right:016x}"
        ));
    }
    if left.batches != right.batches {
        return Some(format!(
            "batches:left={:?};right={:?}",
            left.batches, right.batches
        ));
    }
    if left.prefixes != right.prefixes {
        return Some("public-report-prefix-witnesses".to_string());
    }
    None
}

fn exact_report_bit_delta(
    reference: &SealedStudy,
    mutant: &SealedStudy,
    mutation: Mutation,
) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    if reference.config_digest != mutant.config_digest || mutation.before ^ mutation.after != mask {
        return false;
    }
    let mut expected = reference.record.clone();
    let Some(point) = expected.report_points.get_mut(mutation.evaluation) else {
        return false;
    };
    let Some(bits) = point.get_mut(mutation.coordinate) else {
        return false;
    };
    if *bits != mutation.before {
        return false;
    }
    *bits = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &SealedStudy) -> SeededCorruption {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let acquired = usize::try_from(selector.next_below(usize_u64(ITERATIONS * BATCH_SIZE)))
        .expect("acquired-point index fits usize");
    let evaluation = N_INIT + acquired;
    let coordinate = usize::try_from(selector.next_below(usize_u64(DIMENSION)))
        .expect("coordinate index fits usize");
    let mantissa_bit = u32::try_from(selector.next_below(20)).expect("mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut stale = reference.clone();
    let before = stale.record.report_points[evaluation][coordinate];
    let after = before ^ (1u64 << mantissa_bit);
    stale.record.report_points[evaluation][coordinate] = after;
    let stale_error = stale
        .validate_payload()
        .expect_err("unsealed acquired-point mutation must refuse");
    let run = SealedStudy::seal(stale.config_digest, stale.record);
    let reference_error = run
        .admit_against(reference.output_digest)
        .expect_err("resealed acquired-point mutation must miss retained reference");
    let first_mismatch = first_study_mismatch(&reference.record, &run.record)
        .expect("seeded mutation changes the report frame");
    let semantic_mismatch = run
        .record
        .semantic_mismatch()
        .expect("mutated report point must disagree with callback");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            evaluation,
            coordinate,
            mantissa_bit,
            selector_draws,
            before,
            after,
        },
        stale_error,
        reference_error,
        first_mismatch,
        semantic_mismatch,
    }
}

fn case_inputs(domain: &str, config_digest: u64, digests: &[u64]) -> Vec<u8> {
    let mut bytes = domain.as_bytes().to_vec();
    push_u64(&mut bytes, config_digest);
    push_u64_slice(&mut bytes, digests);
    bytes
}

fn mutation_inputs(reference: &SealedStudy, corruption: &SeededCorruption) -> Vec<u8> {
    let mutation = corruption.mutation;
    let mut bytes = case_inputs(
        "fs-bo-qei-seeded-study-mutation-v1",
        reference.config_digest,
        &[reference.output_digest, corruption.run.output_digest],
    );
    push_u64(&mut bytes, mutation.seed);
    push_u64(&mut bytes, u64::from(mutation.kernel));
    push_u64(&mut bytes, u64::from(mutation.tile));
    push_len(&mut bytes, mutation.evaluation);
    push_len(&mut bytes, mutation.coordinate);
    push_u64(&mut bytes, u64::from(mutation.mantissa_bit));
    push_u64(&mut bytes, mutation.selector_draws);
    push_u64(&mut bytes, mutation.before);
    push_u64(&mut bytes, mutation.after);
    bytes
}

fn mutation_red_report(
    inputs_digest: u64,
    reference: &SealedStudy,
    corruption: &SeededCorruption,
) -> fs_casebook::SuiteReport {
    let mutation = corruption.mutation;
    let details = format!(
        "optimizer_seed=0x{STUDY_SEED:016x}; corruption_seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=report.x[{}][{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; reference_output=0x{:016x}; mutant_output=0x{:016x}; stale_gate={:?}; reference_gate={:?}; first_mismatch={}; semantic_mismatch={}",
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.evaluation,
        mutation.coordinate,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        reference.output_digest,
        corruption.run.output_digest,
        corruption.stale_error,
        corruption.reference_error,
        corruption.first_mismatch,
        corruption.semantic_mismatch,
    );
    let reference_error = corruption.run.admit_against(reference.output_digest);
    Suite::new(SUITE)
        .case(
            "seeded-acquired-point-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                if matches!(
                    reference_error,
                    Err(AdmissionError::ReferenceIdentityMismatch { .. })
                ) {
                    CaseOutcome::fail(details).with_evidence(
                        "crates/fs-bo/tests/qei_study_replay.rs::seeded-acquired-point-corruption",
                    )
                } else {
                    CaseOutcome::pass(format!(
                        "seeded acquired-point corruption escaped: {reference_error:?}"
                    ))
                }
            },
        )
        .run()
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| payload.downcast_ref::<&str>().map(ToString::to_string))
        .unwrap_or_else(|| "non-string panic payload".to_string())
}

#[test]
#[allow(clippy::too_many_lines)]
fn batched_qei_study_replays_and_seeded_corruption_stays_red() {
    let config_digest = fnv1a64(&config_bytes());
    let original = run_study(config_digest);
    let replay = run_study(config_digest);

    let original_mismatch = original.record.semantic_mismatch();
    let replay_mismatch = replay.record.semantic_mismatch();
    let frame_mismatch = first_study_mismatch(&original.record, &replay.record);
    let replay_pass = original_mismatch.is_none()
        && replay_mismatch.is_none()
        && frame_mismatch.is_none()
        && original.validate_payload().is_ok()
        && replay.validate_payload().is_ok()
        && original.admit_against(original.output_digest).is_ok()
        && replay.admit_against(original.output_digest).is_ok()
        && original.output_digest == replay.output_digest
        && original
            .record
            .batches
            .iter()
            .all(|batch| batch.end - batch.start == BATCH_SIZE);

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    let first_inputs = mutation_inputs(&original, &first);
    let second_inputs = mutation_inputs(&replay, &second);
    let first_inputs_digest = fnv1a64(&first_inputs);
    let second_inputs_digest = fnv1a64(&second_inputs);
    let first_red = mutation_red_report(first_inputs_digest, &original, &first);
    let second_red = mutation_red_report(second_inputs_digest, &replay, &second);
    let target = format!(
        "report.x[{}][{}]",
        first.mutation.evaluation, first.mutation.coordinate
    );
    let expected_semantic = format!("evaluation[{}]-callback-report", first.mutation.evaluation);
    let red_stable = first == second
        && first_inputs == second_inputs
        && first_red.records.len() == 1
        && second_red.records.len() == 1
        && first_red.records[0].json_line() == second_red.records[0].json_line();
    let red_is_typed_failure = !first_red.all_passed()
        && first_red.records[0].case == "seeded-acquired-point-corruption"
        && first_red.records[0].inputs_digest == format!("{first_inputs_digest:016x}")
        && first_red.records[0]
            .details
            .contains("PayloadIdentityMismatch")
        && first_red.records[0]
            .details
            .contains("ReferenceIdentityMismatch")
        && first_red.records[0].details.contains(&target)
        && first_red.records[0].details.contains(&expected_semantic);
    let merge_panic = std::panic::catch_unwind(|| first_red.assert_green())
        .expect_err("test-local merge gate must refuse seeded q-EI corruption");
    let merge_message = panic_message(merge_panic.as_ref());
    let mutation_pass = first.run.validate_payload().is_ok()
        && second.run.validate_payload().is_ok()
        && matches!(
            first.stale_error,
            AdmissionError::PayloadIdentityMismatch { declared, computed }
                if declared == original.output_digest && computed == first.run.output_digest
        )
        && matches!(
            first.reference_error,
            AdmissionError::ReferenceIdentityMismatch { expected, found }
                if expected == original.output_digest && found == first.run.output_digest
        )
        && exact_report_bit_delta(&original, &first.run, first.mutation)
        && exact_report_bit_delta(&replay, &second.run, second.mutation)
        && first.mutation.evaluation >= N_INIT
        && first.mutation.evaluation < EXPECTED_EVALUATIONS
        && first.mutation.coordinate < DIMENSION
        && (0..20).contains(&first.mutation.mantissa_bit)
        && f64::from_bits(first.mutation.after).is_finite()
        && first.first_mismatch.starts_with(&target)
        && first.semantic_mismatch == expected_semantic
        && red_stable
        && red_is_typed_failure
        && merge_message.contains("seeded-acquired-point-corruption")
        && merge_message.contains("ReferenceIdentityMismatch")
        && merge_message.contains(&target)
        && merge_message.contains(&expected_semantic);

    let accounting_inputs = case_inputs(
        "fs-bo-qei-study-accounting-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let replay_inputs = case_inputs(
        "fs-bo-qei-study-replay-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let mutation_inputs = case_inputs(
        "fs-bo-qei-study-mutation-pair-case-v1",
        config_digest,
        &[first_inputs_digest, second_inputs_digest],
    );
    let accounting_detail = format!(
        "config=0x{config_digest:016x}; output=0x{:016x}; callbacks={}; report_points={}; report_values={}; best_trace={}; batches={:?}; prefix_lengths={:?}; original_mismatch={original_mismatch:?}; replay_mismatch={replay_mismatch:?}",
        original.output_digest,
        original.record.callbacks.len(),
        original.record.report_points.len(),
        original.record.report_values.len(),
        original.record.best_trace.len(),
        original.record.batches,
        original
            .record
            .prefixes
            .iter()
            .map(|prefix| prefix.report_points.len())
            .collect::<Vec<_>>(),
    );
    let replay_detail = format!(
        "config=0x{config_digest:016x}; original=0x{:016x}; replay=0x{:016x}; frame_mismatch={frame_mismatch:?}; q={BATCH_SIZE}; iterations={ITERATIONS}",
        original.output_digest, replay.output_digest
    );
    let mutation_detail = format!(
        "config=0x{config_digest:016x}; reference=0x{:016x}; mutant=0x{:016x}; mutation_inputs=0x{first_inputs_digest:016x}; seed=0x{MUTATION_SEED:016x}; target={target}; mantissa_bit={}; selector_draws={}; first_mismatch={}; semantic_mismatch={}; red_stable={red_stable}; merge_gate_message={merge_message:?}",
        original.output_digest,
        first.run.output_digest,
        first.mutation.mantissa_bit,
        first.mutation.selector_draws,
        first.first_mismatch,
        first.semantic_mismatch,
    );

    let report = Suite::new(SUITE)
        .case(
            "complete-batched-callback-and-public-report-accounting",
            fnv1a64(&accounting_inputs),
            ToleranceSpec::Exact,
            move || {
                if original_mismatch.is_none() && replay_mismatch.is_none() {
                    CaseOutcome::pass(accounting_detail)
                } else {
                    CaseOutcome::fail(accounting_detail)
                }
                .with_evidence("crates/fs-bo/tests/qei_study_replay.rs::batch-accounting")
            },
        )
        .case(
            "same-seed-complete-qei-output-frame-replay",
            fnv1a64(&replay_inputs),
            ToleranceSpec::Exact,
            move || {
                if replay_pass {
                    CaseOutcome::pass(replay_detail)
                } else {
                    CaseOutcome::fail(replay_detail)
                }
                .with_evidence("crates/fs-bo/tests/qei_study_replay.rs::same-seed-replay")
            },
        )
        .case(
            "seeded-resealed-acquired-point-mutation-is-refused",
            fnv1a64(&mutation_inputs),
            ToleranceSpec::Structural,
            move || {
                if mutation_pass {
                    CaseOutcome::pass(mutation_detail)
                } else {
                    CaseOutcome::fail(mutation_detail)
                }
                .with_evidence("crates/fs-bo/tests/qei_study_replay.rs::seeded-red")
            },
        )
        .run();

    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    assert_eq!(
        report
            .records
            .iter()
            .map(|record| record.case.as_str())
            .collect::<Vec<_>>(),
        [
            "complete-batched-callback-and-public-report-accounting",
            "same-seed-complete-qei-output-frame-replay",
            "seeded-resealed-acquired-point-mutation-is-refused",
        ]
    );
    report.assert_green();
}
