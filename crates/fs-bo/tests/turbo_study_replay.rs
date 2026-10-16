//! Replay-complete production TuRBO restart receipt and seeded red self-test.
//!
//! A constant two-dimensional plateau makes every post-initialization
//! evaluation a failure under TuRBO's strict local-improvement rule. With the
//! retained trust-region lengths and failure tolerance, each round therefore
//! performs three initial evaluations, halves twice, collapses, and restarts.
//! The ten-evaluation fixture completes exactly two such rounds. This test
//! records every objective callback and every public `TurboReport` field,
//! independently reconstructs the strict-tie best, full trace, initializer
//! points, collapse cadence, and restart count, and requires an exact second
//! run. A disclosed mutation then changes one low mantissa bit in `x_best`,
//! reseals the payload, and proves the typed reference and Casebook merge gates
//! refuse the altered receipt.
//!
//! This is finite replay/accounting evidence, not an optimizer-quality claim.
//! It does not generalize to other objectives, configurations, seeds, ISAs,
//! cancellation regimes, persistence formats, or performance envelopes.

use fs_bo::{Matern, TurboConfig, TurboReport, turbo_minimize};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_rand::StreamKey;

const SUITE: &str = "fs-bo/ascent-turbo-study-replay-v1";
const STUDY_CASE: &str = "turbo-flat-plateau-two-forced-restarts";
const STUDY_SEED: u64 = 17;
const MUTATION_SEED: u64 = 0xB0A7_5EED_0000_0031;
const MUTATION_KERNEL: u32 = 0xB031;
const MUTATION_TILE: u32 = 0;

const DIMENSION: usize = 2;
const HYPER_STARTS: usize = 1;
const N_INIT: usize = 3;
const CANDIDATES: usize = 4;
const MAX_EVALUATIONS: usize = 10;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND: f64 = 1.0;
const LOG_BOX_LOWER: f64 = -2.0;
const LOG_BOX_UPPER: f64 = 0.5;
const L_INIT: f64 = 0.5;
const L_MIN: f64 = 0.2;
const L_MAX: f64 = 1.0;
const SUCCESS_TOLERANCE: usize = 2;
const FAILURE_TOLERANCE: usize = 1;
const MAX_LOCAL: usize = 3;
const REFIT_EVERY: u64 = 2;
const PLATEAU_VALUE: f64 = 7.0;
const EXPECTED_FAILURE_STEPS_PER_ROUND: usize = 2;
const EXPECTED_EVALUATIONS_PER_ROUND: usize = N_INIT + EXPECTED_FAILURE_STEPS_PER_ROUND;
const EXPECTED_RESTARTS: usize = MAX_EVALUATIONS / EXPECTED_EVALUATIONS_PER_ROUND;

const _: () = assert!(DIMENSION == 2);
const _: () = assert!(EXPECTED_EVALUATIONS_PER_ROUND == 5);
const _: () = assert!(EXPECTED_RESTARTS == 2);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    point: Vec<u64>,
    value: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    callbacks: Vec<EvaluationBits>,
    report_x_best: Vec<u64>,
    report_f_best: u64,
    report_evals: usize,
    report_restarts: usize,
    report_trace: Vec<u64>,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-bo-turbo-study-output-frame-v1".to_vec();
        push_u64(&mut bytes, config_digest);
        push_len(&mut bytes, self.callbacks.len());
        for callback in &self.callbacks {
            push_u64_slice(&mut bytes, &callback.point);
            push_u64(&mut bytes, callback.value);
        }
        push_u64_slice(&mut bytes, &self.report_x_best);
        push_u64(&mut bytes, self.report_f_best);
        push_u64(&mut bytes, usize_u64(self.report_evals));
        push_u64(&mut bytes, usize_u64(self.report_restarts));
        push_u64_slice(&mut bytes, &self.report_trace);
        bytes
    }

    #[allow(clippy::too_many_lines)] // one complete independent accounting audit
    fn accounting_mismatch(&self) -> Option<String> {
        let config = config();
        let failure_steps = failure_steps_until_collapse(&config);
        if failure_steps != EXPECTED_FAILURE_STEPS_PER_ROUND {
            return Some(format!(
                "collapse-failure-steps:{failure_steps}!=expected-{EXPECTED_FAILURE_STEPS_PER_ROUND}"
            ));
        }
        let evaluations_per_round = config.n_init + failure_steps;
        if evaluations_per_round != EXPECTED_EVALUATIONS_PER_ROUND {
            return Some(format!(
                "evaluations-per-round:{evaluations_per_round}!=expected-{EXPECTED_EVALUATIONS_PER_ROUND}"
            ));
        }
        if !config.max_evals.is_multiple_of(evaluations_per_round) {
            return Some(format!(
                "fixture-budget-does-not-end-on-collapse:{}%{evaluations_per_round}",
                config.max_evals
            ));
        }
        let reconstructed_restarts = config.max_evals / evaluations_per_round;
        if reconstructed_restarts != EXPECTED_RESTARTS {
            return Some(format!(
                "reconstructed-restarts:{reconstructed_restarts}!=expected-{EXPECTED_RESTARTS}"
            ));
        }
        if self.callbacks.len() != config.max_evals {
            return Some(format!(
                "callback-count:{}!=configured-{}",
                self.callbacks.len(),
                config.max_evals
            ));
        }
        if self.report_trace.len() != self.callbacks.len() {
            return Some(format!(
                "report.trace-count:{}!=callbacks-{}",
                self.report_trace.len(),
                self.callbacks.len()
            ));
        }
        if self.report_x_best.len() != DIMENSION {
            return Some(format!(
                "report.x_best-dimension:{}!=expected-{DIMENSION}",
                self.report_x_best.len()
            ));
        }

        for round in 0..EXPECTED_RESTARTS {
            let round_start = round * evaluations_per_round;
            for initialization in 0..config.n_init {
                let callback_index = round_start + initialization;
                let expected = expected_initial_point(&config, round, initialization);
                if self.callbacks[callback_index].point != expected {
                    return Some(format!(
                        "callbacks[{callback_index}]-round-{round}-initializer-{initialization}:recorded={:016x?};expected={expected:016x?}",
                        self.callbacks[callback_index].point
                    ));
                }
            }
        }

        let mut reconstructed_best = f64::INFINITY;
        let mut reconstructed_x = Vec::new();
        for (index, callback) in self.callbacks.iter().enumerate() {
            if callback.point.len() != DIMENSION {
                return Some(format!(
                    "callback[{index}]-dimension:{}!=expected-{DIMENSION}",
                    callback.point.len()
                ));
            }
            let point: Vec<f64> = callback.point.iter().copied().map(f64::from_bits).collect();
            if point
                .iter()
                .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
            {
                return Some(format!(
                    "callback[{index}]-point-outside-box:{:016x?}",
                    callback.point
                ));
            }
            let recomputed = plateau(&point).to_bits();
            if recomputed != callback.value {
                return Some(format!(
                    "callback[{index}]-objective:recomputed=0x{recomputed:016x};recorded=0x{:016x}",
                    callback.value
                ));
            }
            let value = f64::from_bits(callback.value);
            if value < reconstructed_best {
                reconstructed_best = value;
                reconstructed_x.clone_from(&callback.point);
            }
            if self.report_trace[index] != reconstructed_best.to_bits() {
                return Some(format!(
                    "report.trace[{index}]:reconstructed=0x{:016x};reported=0x{:016x}",
                    reconstructed_best.to_bits(),
                    self.report_trace[index]
                ));
            }
        }

        if self.report_evals != self.callbacks.len() {
            return Some(format!(
                "report.evals:{}!=callbacks-{}",
                self.report_evals,
                self.callbacks.len()
            ));
        }
        if self.report_restarts != reconstructed_restarts {
            return Some(format!(
                "report.restarts:{}!=reconstructed-{reconstructed_restarts}",
                self.report_restarts
            ));
        }
        if self.report_f_best != reconstructed_best.to_bits() {
            return Some(format!(
                "report.f_best:0x{:016x}!=reconstructed-0x{:016x}",
                self.report_f_best,
                reconstructed_best.to_bits()
            ));
        }
        if self.report_x_best != reconstructed_x {
            let coordinate = self
                .report_x_best
                .iter()
                .zip(&reconstructed_x)
                .position(|(reported, expected)| reported != expected)
                .unwrap_or(0);
            return Some(format!(
                "report.x_best[{coordinate}]:0x{:016x}!=reconstructed-0x{:016x}",
                self.report_x_best[coordinate], reconstructed_x[coordinate]
            ));
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

fn config() -> TurboConfig {
    TurboConfig {
        bounds: (LOWER_BOUND, UPPER_BOUND),
        family: Matern::FiveHalves,
        log_box: (LOG_BOX_LOWER, LOG_BOX_UPPER),
        hyper_starts: HYPER_STARTS,
        n_init: N_INIT,
        candidates: CANDIDATES,
        max_evals: MAX_EVALUATIONS,
        l_init: L_INIT,
        l_min: L_MIN,
        l_max: L_MAX,
        succ_tol: SUCCESS_TOLERANCE,
        fail_tol: FAILURE_TOLERANCE,
        seed: STUDY_SEED,
        max_local: MAX_LOCAL,
        refit_every: REFIT_EVERY,
    }
}

fn config_bytes() -> Vec<u8> {
    let config = config();
    let mut bytes = b"fs-bo-turbo-study-config-v1".to_vec();
    push_str(&mut bytes, STUDY_CASE);
    push_str(&mut bytes, "dimensionless-normalized-coordinate");
    push_str(&mut bytes, "dimensionless-constant-objective");
    push_str(&mut bytes, "constant-plateau-seven-v1");
    push_u64(&mut bytes, usize_u64(DIMENSION));
    push_u64(&mut bytes, config.bounds.0.to_bits());
    push_u64(&mut bytes, config.bounds.1.to_bits());
    push_str(
        &mut bytes,
        match config.family {
            Matern::Half => "matern-half",
            Matern::ThreeHalves => "matern-three-halves",
            Matern::FiveHalves => "matern-five-halves",
        },
    );
    push_u64(&mut bytes, config.log_box.0.to_bits());
    push_u64(&mut bytes, config.log_box.1.to_bits());
    push_u64(&mut bytes, usize_u64(config.hyper_starts));
    push_u64(&mut bytes, usize_u64(config.n_init));
    push_u64(&mut bytes, usize_u64(config.candidates));
    push_u64(&mut bytes, usize_u64(config.max_evals));
    push_u64(&mut bytes, config.l_init.to_bits());
    push_u64(&mut bytes, config.l_min.to_bits());
    push_u64(&mut bytes, config.l_max.to_bits());
    push_u64(&mut bytes, usize_u64(config.succ_tol));
    push_u64(&mut bytes, usize_u64(config.fail_tol));
    push_u64(&mut bytes, config.seed);
    push_u64(&mut bytes, usize_u64(config.max_local));
    push_u64(&mut bytes, config.refit_every);
    push_str(
        &mut bytes,
        "strict-local-improvement+failure-halving+collapse-restart+global-best-retention-v1",
    );
    push_str(
        &mut bytes,
        "scrambled-Sobol-round-init+QMC-LBFGS-local-fit+ARD-TR+joint-Thompson-Philox-v1",
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
        "no-quality-all-objective-all-config-all-seed-cross-ISA-Cx-persistence-performance-claim",
    );
    bytes
}

fn plateau(_point: &[f64]) -> f64 {
    PLATEAU_VALUE
}

fn failure_steps_until_collapse(config: &TurboConfig) -> usize {
    let mut length = config.l_init;
    let mut failures = 0usize;
    let mut steps = 0usize;
    loop {
        steps += 1;
        failures += 1;
        if failures >= config.fail_tol {
            length *= 0.5;
            failures = 0;
        }
        if length < config.l_min {
            return steps;
        }
    }
}

fn expected_initial_point(config: &TurboConfig, round: usize, initialization: usize) -> Vec<u64> {
    let sobol = fs_rand::qmc::Sobol::scrambled(
        DIMENSION.min(fs_rand::qmc::MAX_SOBOL_DIM),
        config.seed ^ u64::try_from(round).expect("fixture round fits u64"),
    );
    let mut point = vec![0.0f64; DIMENSION];
    sobol.point(
        u32::try_from(initialization + 1).expect("few fixture initializers"),
        &mut point,
    );
    let span = config.bounds.1 - config.bounds.0;
    point
        .iter()
        .map(|unit| span.mul_add(*unit, config.bounds.0).to_bits())
        .collect()
}

fn evaluation_bits(point: &[f64], value: f64) -> EvaluationBits {
    EvaluationBits {
        point: point.iter().map(|value| value.to_bits()).collect(),
        value: value.to_bits(),
    }
}

fn run_study(config_digest: u64) -> SealedStudy {
    let mut callbacks = Vec::with_capacity(MAX_EVALUATIONS);
    let mut objective = |point: &[f64]| {
        let value = plateau(point);
        callbacks.push(evaluation_bits(point, value));
        value
    };
    let report = turbo_minimize(&mut objective, DIMENSION, &config());
    let TurboReport {
        x_best,
        f_best,
        evals,
        restarts,
        trace,
    } = report;
    let record = StudyRecord {
        callbacks,
        report_x_best: x_best.iter().map(|value| value.to_bits()).collect(),
        report_f_best: f_best.to_bits(),
        report_evals: evals,
        report_restarts: restarts,
        report_trace: trace.iter().map(|value| value.to_bits()).collect(),
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
    if left.report_x_best.len() != right.report_x_best.len() {
        return Some(format!(
            "report.x_best.length:{}!={}",
            left.report_x_best.len(),
            right.report_x_best.len()
        ));
    }
    if let Some((coordinate, (a, b))) = left
        .report_x_best
        .iter()
        .zip(&right.report_x_best)
        .enumerate()
        .find(|(_, (a, b))| a != b)
    {
        return Some(format!(
            "report.x_best[{coordinate}]:0x{a:016x}!=0x{b:016x}"
        ));
    }
    if left.report_f_best != right.report_f_best {
        return Some(format!(
            "report.f_best:0x{:016x}!=0x{:016x}",
            left.report_f_best, right.report_f_best
        ));
    }
    if left.report_evals != right.report_evals {
        return Some(format!(
            "report.evals:{}!={}",
            left.report_evals, right.report_evals
        ));
    }
    if left.report_restarts != right.report_restarts {
        return Some(format!(
            "report.restarts:{}!={}",
            left.report_restarts, right.report_restarts
        ));
    }
    if left.report_trace.len() != right.report_trace.len() {
        return Some(format!(
            "report.trace.length:{}!={}",
            left.report_trace.len(),
            right.report_trace.len()
        ));
    }
    left.report_trace
        .iter()
        .zip(&right.report_trace)
        .enumerate()
        .find_map(|(index, (a, b))| {
            (a != b).then(|| format!("report.trace[{index}]:0x{a:016x}!=0x{b:016x}"))
        })
}

fn mutate_report_x_best(reference: &SealedStudy) -> (SealedStudy, Mutation) {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let coordinate =
        usize::try_from(selector.next_below(usize_u64(DIMENSION))).expect("coordinate fits usize");
    let mantissa_bit =
        u32::try_from(selector.next_below(20)).expect("selected mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut record = reference.record.clone();
    let before = record.report_x_best[coordinate];
    let after = before ^ (1u64 << mantissa_bit);
    record.report_x_best[coordinate] = after;
    let mutation = Mutation {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
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
    let Some(&reference_bits) = reference.record.report_x_best.get(mutation.coordinate) else {
        return false;
    };
    let Some(&mutant_bits) = mutant.record.report_x_best.get(mutation.coordinate) else {
        return false;
    };
    if reference.config_digest != mutant.config_digest
        || reference_bits != mutation.before
        || mutant_bits != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }

    let mut expected = reference.record.clone();
    expected.report_x_best[mutation.coordinate] = mutation.after;
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
        "fs-bo-turbo-seeded-study-mutation-v1",
        config_digest,
        &[reference_digest, mutant_digest],
    );
    push_u64(&mut bytes, mutation.seed);
    push_u64(&mut bytes, u64::from(mutation.kernel));
    push_u64(&mut bytes, u64::from(mutation.tile));
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
        "seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=report.x_best[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; reference_output=0x{reference_digest:016x}; mutant_output=0x{:016x}; gate={gate_error:?}; first_mismatch={mismatch}",
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.coordinate,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        mutant.output_digest,
    );
    Suite::new(SUITE)
        .case(
            "seeded-turbo-best-coordinate-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                if matches!(
                    gate_error,
                    Err(AdmissionError::ReferenceIdentityMismatch { .. })
                ) {
                    CaseOutcome::fail(details).with_evidence(
                        "crates/fs-bo/tests/turbo_study_replay.rs::seeded-turbo-best-coordinate-corruption",
                    )
                } else {
                    CaseOutcome::pass(format!(
                        "seeded TuRBO corruption escaped the reference gate: {gate_error:?}"
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
fn production_turbo_study_replays_and_seeded_corruption_stays_red() {
    let config_frame = config_bytes();
    let config_digest = fnv1a64(&config_frame);
    let original = run_study(config_digest);
    let replay = run_study(config_digest);

    let original_accounting = original.record.accounting_mismatch();
    let replay_accounting = replay.record.accounting_mismatch();
    let replay_mismatch = first_study_mismatch(&original.record, &replay.record);
    let replay_pass = original.validate_payload().is_ok()
        && replay.validate_payload().is_ok()
        && original == replay
        && replay_mismatch.is_none();

    let (mutant, mutation) = mutate_report_x_best(&original);
    let (mutant_replayed, mutation_replayed) = mutate_report_x_best(&replay);
    let mutant_mismatch = first_study_mismatch(&original.record, &mutant.record);
    let mutant_replayed_mismatch = first_study_mismatch(&replay.record, &mutant_replayed.record);
    let mutant_accounting = mutant.record.accounting_mismatch();
    let mutant_replayed_accounting = mutant_replayed.record.accounting_mismatch();
    let expected_mismatch = format!("report.x_best[{}]", mutation.coordinate);
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
        "fs-bo-turbo-seeded-study-mutation-pair-case-v1",
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
        && red_first.records[0].case == "seeded-turbo-best-coordinate-corruption"
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
    let before = f64::from_bits(mutation.before);
    let after = f64::from_bits(mutation.after);
    let mutation_pass = mutant.validate_payload().is_ok()
        && mutant_replayed.validate_payload().is_ok()
        && stale_payload_identity_is_refused(&mutant)
        && mutant == mutant_replayed
        && mutation == mutation_replayed
        && is_exact_report_bit_delta(&original, &mutant, mutation)
        && is_exact_report_bit_delta(&replay, &mutant_replayed, mutation_replayed)
        && mutation.coordinate < DIMENSION
        && (0..20).contains(&mutation.mantissa_bit)
        && mutation.before != mutation.after
        && before.is_finite()
        && after.is_finite()
        && (LOWER_BOUND..=UPPER_BOUND).contains(&after)
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
        && mutant_accounting
            .as_deref()
            .is_some_and(|mismatch| mismatch.starts_with(&expected_mismatch))
        && mutant_replayed_accounting
            .as_deref()
            .is_some_and(|mismatch| mismatch.starts_with(&expected_mismatch))
        && red_lines_stable
        && red_is_typed_failure
        && merge_gate_message.contains("seeded-turbo-best-coordinate-corruption")
        && merge_gate_message.contains("ReferenceIdentityMismatch");

    let accounting_inputs = case_inputs(
        "fs-bo-turbo-study-accounting-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let replay_inputs = case_inputs(
        "fs-bo-turbo-study-replay-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let accounting_pass = original_accounting.is_none() && replay_accounting.is_none();
    let accounting_detail = format!(
        "config=0x{config_digest:016x}; original=0x{:016x}; replay=0x{:016x}; callbacks={}; x_best={:016x?}; f_best=0x{:016x}; evals={}; restarts={}; trace={}; original_mismatch={original_accounting:?}; replay_mismatch={replay_accounting:?}",
        original.output_digest,
        replay.output_digest,
        original.record.callbacks.len(),
        original.record.report_x_best,
        original.record.report_f_best,
        original.record.report_evals,
        original.record.report_restarts,
        original.record.report_trace.len(),
    );
    let replay_detail = format!(
        "config=0x{config_digest:016x}; original=0x{:016x}; replay=0x{:016x}; first_mismatch={replay_mismatch:?}",
        original.output_digest, replay.output_digest
    );
    let mutation_detail = format!(
        "config=0x{config_digest:016x}; reference=0x{:016x}; replay=0x{:016x}; mutant=0x{:016x}; replay_mutant=0x{:016x}; mutation_inputs=0x{mutation_inputs_digest:016x}; replay_mutation_inputs=0x{mutation_replayed_inputs_digest:016x}; seed=0x{MUTATION_SEED:016x}; target=report.x_best[{}]; mantissa_bit={}; selector_draws={}; first_mismatch={mutant_mismatch:?}; replay_first_mismatch={mutant_replayed_mismatch:?}; mutant_accounting={mutant_accounting:?}; replay_mutant_accounting={mutant_replayed_accounting:?}; red_record_stable={red_lines_stable}; merge_gate_message={merge_gate_message:?}",
        original.output_digest,
        replay.output_digest,
        mutant.output_digest,
        mutant_replayed.output_digest,
        mutation.coordinate,
        mutation.mantissa_bit,
        mutation.selector_draws,
    );

    let report = Suite::new(SUITE)
        .case(
            "callback-restart-and-public-report-accounting",
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
            "same-seed-full-turbo-output-frame-replay",
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
            "seeded-resealed-turbo-best-mutation-is-refused",
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
            "callback-restart-and-public-report-accounting",
            "same-seed-full-turbo-output-frame-replay",
            "seeded-resealed-turbo-best-mutation-is-refused",
        ]
    );
    report.assert_green();
}
