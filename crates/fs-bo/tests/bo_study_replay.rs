//! Replay-complete production BO study receipt and seeded red self-test.
//!
//! This battery runs the same short sequential-EI Branin configuration used by
//! the crate's retained golden.  It records every objective callback and every
//! public `BoReport` field, independently reconstructs the evaluation and
//! best-so-far accounting, and compares a complete canonical output frame
//! across same-seed runs.  A disclosed mutation then changes one mantissa bit
//! in an initial reported point, reseals the canonical payload so its declared
//! identity remains valid, and proves the reference-identity admission gate
//! and Casebook merge gate both refuse it.
//!
//! The fixture establishes replay and evidence plumbing for this one finite
//! study.  It adds no optimizer-quality, all-objective, all-seed, cross-ISA,
//! cancellation, persistence, authenticated-ledger, or performance claim.

use fs_bo::{BoConfig, BoReport, Matern, minimize};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_rand::StreamKey;

const SUITE: &str = "fs-bo/ascent-study-replay-v1";
const STUDY_CASE: &str = "sequential-ei-branin-short-golden";
const STUDY_SEED: u64 = 5;
const MUTATION_SEED: u64 = 0xB0A7_5EED_0000_0024;
const MUTATION_KERNEL: u32 = 0xB024;
const MUTATION_TILE: u32 = 0;

const DIMENSION: usize = 2;
const N_INIT: usize = 6;
const ITERATIONS: usize = 2;
const BATCH_SIZE: usize = 1;
const EXPECTED_EVALUATIONS: usize = N_INIT + ITERATIONS * BATCH_SIZE;
const EXPECTED_TRACE_POINTS: usize = ITERATIONS + 1;
const HYPER_STARTS: usize = 2;
const ACQUISITION_STARTS: usize = 2;
const ACQUISITION_EVALUATIONS: usize = 150;
const MC_SAMPLES: usize = 128;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND: f64 = 1.0;
const LOG_BOX_LOWER: f64 = -2.0;
const LOG_BOX_UPPER: f64 = 0.5;

const _: () = assert!(DIMENSION == 2);
const _: () = assert!(BATCH_SIZE == 1);
const _: () = assert!(EXPECTED_EVALUATIONS == 8 && EXPECTED_TRACE_POINTS == 3);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    point: Vec<u64>,
    value: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    callbacks: Vec<EvaluationBits>,
    report_points: Vec<Vec<u64>>,
    report_values: Vec<u64>,
    best_trace: Vec<u64>,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-bo-study-output-frame-v1".to_vec();
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
        bytes
    }

    fn accounting_mismatch(&self) -> Option<String> {
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

        for index in 0..EXPECTED_EVALUATIONS {
            let callback = &self.callbacks[index];
            let point = &self.report_points[index];
            if callback.point.len() != DIMENSION || point.len() != DIMENSION {
                return Some(format!(
                    "evaluation[{index}]-dimension:callback-{};report-{};expected-{DIMENSION}",
                    callback.point.len(),
                    point.len()
                ));
            }
            if callback.point != *point || callback.value != self.report_values[index] {
                return Some(format!(
                    "evaluation[{index}]-callback-report:callback={callback:?};report-point={point:016x?};report-value=0x{:016x}",
                    self.report_values[index]
                ));
            }
            let decoded: Vec<f64> = point.iter().copied().map(f64::from_bits).collect();
            if decoded
                .iter()
                .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
            {
                return Some(format!(
                    "evaluation[{index}]-point-outside-box:{point:016x?}"
                ));
            }
            let recomputed = branin(&decoded).to_bits();
            if recomputed != self.report_values[index] {
                return Some(format!(
                    "evaluation[{index}]-objective:recomputed=0x{recomputed:016x};reported=0x{:016x}",
                    self.report_values[index]
                ));
            }
        }

        let mut best = f64::INFINITY;
        for &value in &self.report_values[..N_INIT] {
            best = best.min(f64::from_bits(value));
        }
        if best.to_bits() != self.best_trace[0] {
            return Some(format!(
                "best-trace[0]:recomputed=0x{:016x};reported=0x{:016x}",
                best.to_bits(),
                self.best_trace[0]
            ));
        }
        for iteration in 0..ITERATIONS {
            let start = N_INIT + iteration * BATCH_SIZE;
            let end = start + BATCH_SIZE;
            for &value in &self.report_values[start..end] {
                best = best.min(f64::from_bits(value));
            }
            if best.to_bits() != self.best_trace[iteration + 1] {
                return Some(format!(
                    "best-trace[{}]:recomputed=0x{:016x};reported=0x{:016x}",
                    iteration + 1,
                    best.to_bits(),
                    self.best_trace[iteration + 1]
                ));
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
    let mut bytes = b"fs-bo-study-config-v1".to_vec();
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
    push_str(
        &mut bytes,
        "Sobol-init+per-iteration-standardize+QMC-LBFGS-hyperfit+EI+CMA-ES-v1",
    );
    push_str(&mut bytes, "mc-samples-bound-but-unused-for-q-one");
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
        "no-quality-all-objective-all-config-all-seed-cross-ISA-Cx-persistence-performance-claim",
    );
    bytes
}

fn branin(point: &[f64]) -> f64 {
    let x1 = 15.0f64.mul_add(point[0], -5.0);
    let x2 = 15.0 * point[1];
    let a = 1.0;
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let r = 6.0;
    let s = 10.0;
    let t = 1.0 / (8.0 * core::f64::consts::PI);
    let inner = b.mul_add(-(x1 * x1), c.mul_add(x1, x2 - r));
    a * inner * inner + s * (1.0 - t) * fs_math::det::cos(x1) + s
}

fn evaluation_bits(point: &[f64], value: f64) -> EvaluationBits {
    EvaluationBits {
        point: point.iter().map(|value| value.to_bits()).collect(),
        value: value.to_bits(),
    }
}

fn run_study(config_digest: u64) -> SealedStudy {
    let mut callbacks = Vec::with_capacity(EXPECTED_EVALUATIONS);
    let mut objective = |point: &[f64]| {
        let value = branin(point);
        callbacks.push(evaluation_bits(point, value));
        value
    };
    let config = config();
    let report = minimize(&mut objective, DIMENSION, N_INIT, ITERATIONS, &config);
    let BoReport { x, y, best_trace } = report;
    let record = StudyRecord {
        callbacks,
        report_points: x
            .iter()
            .map(|point| point.iter().map(|value| value.to_bits()).collect())
            .collect(),
        report_values: y.iter().map(|value| value.to_bits()).collect(),
        best_trace: best_trace.iter().map(|value| value.to_bits()).collect(),
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
    for (index, (a, b)) in left.callbacks.iter().zip(&right.callbacks).enumerate() {
        if a != b {
            return Some(format!("callbacks[{index}]:left={a:?};right={b:?}"));
        }
    }
    if left.report_points.len() != right.report_points.len() {
        return Some(format!(
            "report.x.length:{}!={}",
            left.report_points.len(),
            right.report_points.len()
        ));
    }
    for (index, (a, b)) in left
        .report_points
        .iter()
        .zip(&right.report_points)
        .enumerate()
    {
        if a != b {
            if a.len() != b.len() {
                return Some(format!("report.x[{index}].length:{}!={}", a.len(), b.len()));
            }
            let coordinate = a.iter().zip(b).position(|(x, y)| x != y).unwrap_or(0);
            return Some(format!(
                "report.x[{index}][{coordinate}]:0x{:016x}!=0x{:016x}",
                a[coordinate], b[coordinate]
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
    if let Some((index, (a, b))) = left
        .report_values
        .iter()
        .zip(&right.report_values)
        .enumerate()
        .find(|(_, (a, b))| a != b)
    {
        return Some(format!("report.y[{index}]:0x{a:016x}!=0x{b:016x}"));
    }
    if left.best_trace.len() != right.best_trace.len() {
        return Some(format!(
            "report.best_trace.length:{}!={}",
            left.best_trace.len(),
            right.best_trace.len()
        ));
    }
    left.best_trace
        .iter()
        .zip(&right.best_trace)
        .enumerate()
        .find_map(|(index, (a, b))| {
            (a != b).then(|| format!("report.best_trace[{index}]:0x{a:016x}!=0x{b:016x}"))
        })
}

fn mutate_report_point(reference: &SealedStudy) -> (SealedStudy, Mutation) {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let evaluation = usize::try_from(selector.next_below(usize_u64(N_INIT)))
        .expect("initial evaluation index fits usize");
    let coordinate = usize::try_from(selector.next_below(usize_u64(DIMENSION)))
        .expect("coordinate index fits usize");
    let mantissa_bit =
        u32::try_from(selector.next_below(20)).expect("selected mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut record = reference.record.clone();
    let before = record.report_points[evaluation][coordinate];
    let after = before ^ (1u64 << mantissa_bit);
    record.report_points[evaluation][coordinate] = after;
    let mutation = Mutation {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
        evaluation,
        coordinate,
        mantissa_bit,
        selector_draws,
        before,
        after,
    };
    (SealedStudy::seal(reference.config_digest, record), mutation)
}

fn is_exact_report_bit_delta(
    reference: &SealedStudy,
    mutant: &SealedStudy,
    mutation: Mutation,
) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let Some(reference_point) = reference.record.report_points.get(mutation.evaluation) else {
        return false;
    };
    let Some(&reference_bit_pattern) = reference_point.get(mutation.coordinate) else {
        return false;
    };
    let Some(mutant_point) = mutant.record.report_points.get(mutation.evaluation) else {
        return false;
    };
    let Some(&mutant_bit_pattern) = mutant_point.get(mutation.coordinate) else {
        return false;
    };
    if reference.config_digest != mutant.config_digest
        || reference_bit_pattern != mutation.before
        || mutant_bit_pattern != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }

    let mut expected = reference.record.clone();
    expected.report_points[mutation.evaluation][mutation.coordinate] = mutation.after;
    expected == mutant.record
}

fn stale_payload_identity_is_refused(sealed: &SealedStudy) -> bool {
    let expected_computed = sealed.output_digest;
    let expected_declared = expected_computed ^ 1;
    let mut stale = sealed.clone();
    stale.output_digest = expected_declared;
    matches!(
        stale.validate_payload(),
        Err(AdmissionError::PayloadIdentityMismatch { declared, computed })
            if declared == expected_declared && computed == expected_computed
    )
}

fn case_inputs(domain: &str, config_digest: u64, digests: &[u64]) -> Vec<u8> {
    let mut bytes = domain.as_bytes().to_vec();
    push_u64(&mut bytes, config_digest);
    push_u64_slice(&mut bytes, digests);
    bytes
}

fn mutation_inputs(
    config_digest: u64,
    reference_digest: u64,
    mutant_digest: u64,
    mutation: Mutation,
) -> Vec<u8> {
    let mut bytes = case_inputs(
        "fs-bo-seeded-study-mutation-v1",
        config_digest,
        &[reference_digest, mutant_digest],
    );
    push_u64(&mut bytes, mutation.seed);
    push_u64(&mut bytes, u64::from(mutation.kernel));
    push_u64(&mut bytes, u64::from(mutation.tile));
    push_u64(&mut bytes, usize_u64(mutation.evaluation));
    push_u64(&mut bytes, usize_u64(mutation.coordinate));
    push_u64(&mut bytes, u64::from(mutation.mantissa_bit));
    push_u64(&mut bytes, mutation.selector_draws);
    push_u64(&mut bytes, mutation.before);
    push_u64(&mut bytes, mutation.after);
    bytes
}

fn mutation_red_report(
    inputs_digest: u64,
    reference_digest: u64,
    mutant: &SealedStudy,
    mutation: Mutation,
    mismatch: &str,
) -> fs_casebook::SuiteReport {
    let gate_error = mutant.admit_against(reference_digest);
    let details = format!(
        "seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=report.x[{}][{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; reference_output=0x{reference_digest:016x}; mutant_output=0x{:016x}; gate={gate_error:?}; first_mismatch={mismatch}",
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.evaluation,
        mutation.coordinate,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        mutant.output_digest,
    );
    Suite::new(SUITE)
        .case(
            "seeded-report-coordinate-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                if matches!(
                    gate_error,
                    Err(AdmissionError::ReferenceIdentityMismatch { .. })
                ) {
                    CaseOutcome::fail(details).with_evidence(
                        "crates/fs-bo/tests/bo_study_replay.rs::seeded-report-coordinate-corruption",
                    )
                } else {
                    CaseOutcome::pass(format!(
                        "seeded report corruption escaped the reference gate: {gate_error:?}"
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
fn production_bo_study_replays_and_seeded_corruption_stays_red() {
    let config_frame = config_bytes();
    let config_digest = fnv1a64(&config_frame);
    let original = run_study(config_digest);
    let replay = run_study(config_digest);

    let original_accounting = original.record.accounting_mismatch();
    let replay_accounting = replay.record.accounting_mismatch();
    let replay_mismatch = first_study_mismatch(&original.record, &replay.record);
    let replay_pass = original.validate_payload().is_ok()
        && replay.validate_payload().is_ok()
        && original.output_digest == replay.output_digest
        && replay_mismatch.is_none();

    let (mutant, mutation) = mutate_report_point(&original);
    let (mutant_replayed, mutation_replayed) = mutate_report_point(&replay);
    let mutant_mismatch = first_study_mismatch(&original.record, &mutant.record);
    let mutant_replayed_mismatch = first_study_mismatch(&replay.record, &mutant_replayed.record);
    let expected_mismatch = format!("report.x[{}][{}]", mutation.evaluation, mutation.coordinate);
    let mutation_input_frame = mutation_inputs(
        config_digest,
        original.output_digest,
        mutant.output_digest,
        mutation,
    );
    let mutation_inputs_digest = fnv1a64(&mutation_input_frame);
    let mutation_replayed_input_frame = mutation_inputs(
        config_digest,
        replay.output_digest,
        mutant_replayed.output_digest,
        mutation_replayed,
    );
    let mutation_replayed_inputs_digest = fnv1a64(&mutation_replayed_input_frame);
    let mutation_case_inputs = case_inputs(
        "fs-bo-seeded-study-mutation-pair-case-v1",
        config_digest,
        &[mutation_inputs_digest, mutation_replayed_inputs_digest],
    );
    let mutation_case_inputs_digest = fnv1a64(&mutation_case_inputs);
    let mismatch_text = mutant_mismatch.as_deref().unwrap_or("none");
    let replayed_mismatch_text = mutant_replayed_mismatch.as_deref().unwrap_or("none");
    let red_first = mutation_red_report(
        mutation_inputs_digest,
        original.output_digest,
        &mutant,
        mutation,
        mismatch_text,
    );
    let red_second = mutation_red_report(
        mutation_replayed_inputs_digest,
        replay.output_digest,
        &mutant_replayed,
        mutation_replayed,
        replayed_mismatch_text,
    );
    let red_lines_stable = red_first.records.len() == 1
        && red_second.records.len() == 1
        && mutation_input_frame == mutation_replayed_input_frame
        && red_first.records[0].json_line() == red_second.records[0].json_line();
    let red_is_typed_failure = !red_first.all_passed()
        && red_first.records[0].case == "seeded-report-coordinate-corruption"
        && red_first.records[0].inputs_digest == format!("{mutation_inputs_digest:016x}")
        && red_first.records[0]
            .details
            .contains("ReferenceIdentityMismatch")
        && red_first.records[0].details.contains(&expected_mismatch);
    let merge_gate_panic = std::panic::catch_unwind(|| red_first.assert_green());
    let merge_gate_message = merge_gate_panic
        .as_ref()
        .err()
        .map(|payload| panic_message(payload.as_ref()))
        .unwrap_or_default();
    let mutation_pass = mutant.validate_payload().is_ok()
        && mutant_replayed.validate_payload().is_ok()
        && stale_payload_identity_is_refused(&mutant)
        && mutant == mutant_replayed
        && mutation == mutation_replayed
        && is_exact_report_bit_delta(&original, &mutant, mutation)
        && is_exact_report_bit_delta(&replay, &mutant_replayed, mutation_replayed)
        && mutation.evaluation < N_INIT
        && mutation.coordinate < DIMENSION
        && (0..20).contains(&mutation.mantissa_bit)
        && mutation.before != mutation.after
        && f64::from_bits(mutation.after).is_finite()
        && (LOWER_BOUND..UPPER_BOUND).contains(&f64::from_bits(mutation.after))
        && mutant.output_digest != original.output_digest
        && matches!(
            mutant.admit_against(original.output_digest),
            Err(AdmissionError::ReferenceIdentityMismatch {
                expected,
                found
            }) if expected == original.output_digest && found == mutant.output_digest
        )
        && mutant_mismatch
            .as_deref()
            .is_some_and(|mismatch| mismatch.starts_with(&expected_mismatch))
        && red_lines_stable
        && red_is_typed_failure
        && merge_gate_message.contains("seeded-report-coordinate-corruption")
        && merge_gate_message.contains("ReferenceIdentityMismatch");

    let accounting_inputs = case_inputs(
        "fs-bo-study-accounting-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let replay_inputs = case_inputs(
        "fs-bo-study-replay-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let accounting_pass = original_accounting.is_none() && replay_accounting.is_none();
    let accounting_detail = format!(
        "config=0x{config_digest:016x}; original=0x{:016x}; replay=0x{:016x}; callbacks={}; report_points={}; report_values={}; best_trace={}; original_mismatch={original_accounting:?}; replay_mismatch={replay_accounting:?}",
        original.output_digest,
        replay.output_digest,
        original.record.callbacks.len(),
        original.record.report_points.len(),
        original.record.report_values.len(),
        original.record.best_trace.len(),
    );
    let replay_detail = format!(
        "config=0x{config_digest:016x}; original=0x{:016x}; replay=0x{:016x}; first_mismatch={replay_mismatch:?}",
        original.output_digest, replay.output_digest
    );
    let mutation_detail = format!(
        "config=0x{config_digest:016x}; reference=0x{:016x}; replay=0x{:016x}; mutant=0x{:016x}; replay_mutant=0x{:016x}; mutation_inputs=0x{mutation_inputs_digest:016x}; replay_mutation_inputs=0x{mutation_replayed_inputs_digest:016x}; seed=0x{MUTATION_SEED:016x}; target=report.x[{}][{}]; mantissa_bit={}; selector_draws={}; first_mismatch={mutant_mismatch:?}; replay_first_mismatch={mutant_replayed_mismatch:?}; red_record_stable={red_lines_stable}; merge_gate_message={merge_gate_message:?}",
        original.output_digest,
        replay.output_digest,
        mutant.output_digest,
        mutant_replayed.output_digest,
        mutation.evaluation,
        mutation.coordinate,
        mutation.mantissa_bit,
        mutation.selector_draws,
    );

    let report = Suite::new(SUITE)
        .case(
            "complete-callback-and-public-report-accounting",
            fnv1a64(&accounting_inputs),
            ToleranceSpec::Exact,
            move || {
                if accounting_pass {
                    CaseOutcome::pass(accounting_detail)
                } else {
                    CaseOutcome::fail(accounting_detail)
                }
                .with_evidence("crates/fs-bo/CONTRACT.md#conformance-tests")
            },
        )
        .case(
            "same-seed-full-output-frame-replay",
            fnv1a64(&replay_inputs),
            ToleranceSpec::Exact,
            move || {
                if replay_pass {
                    CaseOutcome::pass(replay_detail)
                } else {
                    CaseOutcome::fail(replay_detail)
                }
                .with_evidence("crates/fs-bo/CONTRACT.md#determinism-class")
            },
        )
        .case(
            "seeded-resealed-mutation-is-refused",
            mutation_case_inputs_digest,
            ToleranceSpec::Structural,
            move || {
                if mutation_pass {
                    CaseOutcome::pass(mutation_detail)
                } else {
                    CaseOutcome::fail(mutation_detail)
                }
                .with_evidence("crates/fs-bo/CONTRACT.md#no-claim-boundaries")
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
            "complete-callback-and-public-report-accounting",
            "same-seed-full-output-frame-replay",
            "seeded-resealed-mutation-is-refused",
        ]
    );
    report.assert_green();
}
