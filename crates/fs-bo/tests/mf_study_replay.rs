//! Replay-complete production multi-fidelity BO receipt and seeded red self-test.
//!
//! This battery runs one cost-aware allocation step after a deterministic
//! two-fidelity Branin initialization. It records every objective callback,
//! including the chosen fidelity, and every public `MfReport` field. Independent
//! accounting reconstructs objective values, fidelity counts, cumulative cost,
//! the best-high trace, and the final learned correlation from the pre-step
//! observations. A disclosed mutation changes one finite best-high trace bit,
//! reseals the canonical payload, and proves both the typed reference gate and
//! the Casebook merge gate refuse the altered receipt.
//!
//! The fixture establishes replay and evidence plumbing for this one finite
//! study. It adds no optimizer-quality, all-objective, all-configuration,
//! all-seed, cross-ISA, cancellation, persistence, authenticated-ledger, or
//! performance claim.

use fs_bo::{MfConfig, MfReport, fit_mf, mf_minimize};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_rand::StreamKey;

const SUITE: &str = "fs-bo/ascent-mf-study-replay-v1";
const STUDY_CASE: &str = "cost-aware-two-fidelity-branin-one-step";
const STUDY_SEED: u64 = 5;
const MUTATION_SEED: u64 = 0xB0A7_5EED_0000_0030;
const MUTATION_KERNEL: u32 = 0xB030;
const MUTATION_TILE: u32 = 0;

const DIMENSION: usize = 2;
const N_INIT_LOW: usize = 10;
const N_INIT_HIGH: usize = 4;
const INITIAL_EVALUATIONS: usize = N_INIT_LOW + N_INIT_HIGH;
const SEQUENTIAL_EVALUATIONS: usize = 1;
const EXPECTED_EVALUATIONS: usize = INITIAL_EVALUATIONS + SEQUENTIAL_EVALUATIONS;
const HYPER_STARTS: usize = 2;
const ACQUISITION_STARTS: usize = 2;
const ACQUISITION_EVALUATIONS: usize = 200;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND: f64 = 1.0;
const LOG_BOX_LOWER: f64 = -2.0;
const LOG_BOX_UPPER: f64 = 1.0;
const COST_LOW: f64 = 1.0;
const COST_HIGH: f64 = 10.0;
const INITIAL_COST: f64 = N_INIT_LOW as f64 * COST_LOW + N_INIT_HIGH as f64 * COST_HIGH;
const TOTAL_BUDGET: f64 = INITIAL_COST + COST_LOW;
const FIRST_ITERATION_HYPER_SEED_XOR: u64 = 0x9e37_79b9;

const _: () = assert!(DIMENSION == 2);
const _: () = assert!(INITIAL_EVALUATIONS == 14);
const _: () = assert!(EXPECTED_EVALUATIONS == 15);

#[derive(Debug, Clone, PartialEq, Eq)]
struct CallbackBits {
    point: Vec<u64>,
    fidelity: usize,
    value: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    callbacks: Vec<CallbackBits>,
    f_best_high: u64,
    evals_low: usize,
    evals_high: usize,
    cost: u64,
    trace: Vec<(u64, u64)>,
    learned_correlation: u64,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-bo-mf-study-output-frame-v1".to_vec();
        push_u64(&mut bytes, config_digest);
        push_len(&mut bytes, self.callbacks.len());
        for callback in &self.callbacks {
            push_u64_slice(&mut bytes, &callback.point);
            push_u64(&mut bytes, usize_u64(callback.fidelity));
            push_u64(&mut bytes, callback.value);
        }
        push_u64(&mut bytes, self.f_best_high);
        push_u64(&mut bytes, usize_u64(self.evals_low));
        push_u64(&mut bytes, usize_u64(self.evals_high));
        push_u64(&mut bytes, self.cost);
        push_len(&mut bytes, self.trace.len());
        for &(cost, best_high) in &self.trace {
            push_u64(&mut bytes, cost);
            push_u64(&mut bytes, best_high);
        }
        push_u64(&mut bytes, self.learned_correlation);
        bytes
    }

    #[allow(clippy::too_many_lines)] // one complete independent accounting audit
    fn accounting_mismatch(&self) -> Option<String> {
        if self.callbacks.len() != EXPECTED_EVALUATIONS {
            return Some(format!(
                "callback-count:{}!=expected-{EXPECTED_EVALUATIONS}",
                self.callbacks.len()
            ));
        }
        if self.trace.len() != EXPECTED_EVALUATIONS {
            return Some(format!(
                "report.trace-count:{}!=expected-{EXPECTED_EVALUATIONS}",
                self.trace.len()
            ));
        }

        let mut reconstructed_cost = 0.0f64;
        let mut reconstructed_best_high = f64::INFINITY;
        let mut reconstructed_low = 0usize;
        let mut reconstructed_high = 0usize;

        for (index, callback) in self.callbacks.iter().enumerate() {
            if callback.point.len() != DIMENSION {
                return Some(format!(
                    "callback[{index}]-dimension:{}!=expected-{DIMENSION}",
                    callback.point.len()
                ));
            }
            let expected_fidelity = if index < N_INIT_LOW {
                0
            } else if index < INITIAL_EVALUATIONS {
                1
            } else {
                // This retained fixture has exactly one affordable unit of
                // budget after initialization. Its strongly correlated cheap
                // model must therefore win the MFEI cost-aware allocation.
                0
            };
            if callback.fidelity != expected_fidelity {
                return Some(format!(
                    "callback[{index}]-fidelity:{}!=expected-{expected_fidelity}",
                    callback.fidelity
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
            let recomputed_value = objective(&point, callback.fidelity).to_bits();
            if recomputed_value != callback.value {
                return Some(format!(
                    "callback[{index}]-objective:recomputed=0x{recomputed_value:016x};recorded=0x{:016x}",
                    callback.value
                ));
            }

            match callback.fidelity {
                0 => {
                    reconstructed_low += 1;
                    reconstructed_cost += COST_LOW;
                }
                1 => {
                    reconstructed_high += 1;
                    reconstructed_cost += COST_HIGH;
                    reconstructed_best_high =
                        reconstructed_best_high.min(f64::from_bits(callback.value));
                }
                other => return Some(format!("callback[{index}]-invalid-fidelity:{other}")),
            }

            let (reported_cost, reported_best_high) = self.trace[index];
            if reported_cost != reconstructed_cost.to_bits() {
                return Some(format!(
                    "report.trace[{index}].cost:reconstructed=0x{:016x};reported=0x{reported_cost:016x}",
                    reconstructed_cost.to_bits()
                ));
            }
            if reported_best_high != reconstructed_best_high.to_bits() {
                return Some(format!(
                    "report.trace[{index}].best_high:reconstructed=0x{:016x};reported=0x{reported_best_high:016x}",
                    reconstructed_best_high.to_bits()
                ));
            }
        }

        if self.evals_low != reconstructed_low {
            return Some(format!(
                "report.evals_low:{}!=reconstructed-{reconstructed_low}",
                self.evals_low
            ));
        }
        if self.evals_high != reconstructed_high {
            return Some(format!(
                "report.evals_high:{}!=reconstructed-{reconstructed_high}",
                self.evals_high
            ));
        }
        if self.cost != reconstructed_cost.to_bits() {
            return Some(format!(
                "report.cost:0x{:016x}!=reconstructed-0x{:016x}",
                self.cost,
                reconstructed_cost.to_bits()
            ));
        }
        if self.f_best_high != reconstructed_best_high.to_bits() {
            return Some(format!(
                "report.f_best_high:0x{:016x}!=reconstructed-0x{:016x}",
                self.f_best_high,
                reconstructed_best_high.to_bits()
            ));
        }

        let correlation = independently_reconstructed_correlation(&self.callbacks);
        if self.learned_correlation != correlation.to_bits() {
            return Some(format!(
                "report.learned_correlation:0x{:016x}!=reconstructed-0x{:016x}",
                self.learned_correlation,
                correlation.to_bits()
            ));
        }
        let learned = f64::from_bits(self.learned_correlation);
        if !learned.is_finite() || learned.abs() > 1.0 + 1e-12 {
            return Some(format!(
                "report.learned_correlation-out-of-range:{learned:?}"
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
    trace_index: usize,
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

fn config() -> MfConfig {
    MfConfig {
        bounds: (LOWER_BOUND, UPPER_BOUND),
        log_box: (LOG_BOX_LOWER, LOG_BOX_UPPER),
        hyper_starts: HYPER_STARTS,
        n_init_low: N_INIT_LOW,
        n_init_high: N_INIT_HIGH,
        cost_low: COST_LOW,
        cost_high: COST_HIGH,
        budget: TOTAL_BUDGET,
        acq_starts: ACQUISITION_STARTS,
        acq_evals: ACQUISITION_EVALUATIONS,
        seed: STUDY_SEED,
    }
}

fn config_bytes() -> Vec<u8> {
    let config = config();
    let mut bytes = b"fs-bo-mf-study-config-v1".to_vec();
    push_str(&mut bytes, STUDY_CASE);
    push_str(&mut bytes, "dimensionless-normalized-coordinate");
    push_str(&mut bytes, "dimensionless-branin-value");
    push_str(&mut bytes, "branin-rescaled-from-standard-domain-v1");
    push_str(&mut bytes, "branin-low-smooth-linear-bias-v1");
    push_str(&mut bytes, "fidelity-0-low-fidelity-1-high");
    push_u64(&mut bytes, usize_u64(DIMENSION));
    push_u64(&mut bytes, config.bounds.0.to_bits());
    push_u64(&mut bytes, config.bounds.1.to_bits());
    push_u64(&mut bytes, config.log_box.0.to_bits());
    push_u64(&mut bytes, config.log_box.1.to_bits());
    push_u64(&mut bytes, usize_u64(config.hyper_starts));
    push_u64(&mut bytes, usize_u64(config.n_init_low));
    push_u64(&mut bytes, usize_u64(config.n_init_high));
    push_u64(&mut bytes, config.cost_low.to_bits());
    push_u64(&mut bytes, config.cost_high.to_bits());
    push_u64(&mut bytes, config.budget.to_bits());
    push_u64(&mut bytes, usize_u64(config.acq_starts));
    push_u64(&mut bytes, usize_u64(config.acq_evals));
    push_u64(&mut bytes, config.seed);
    push_str(
        &mut bytes,
        "scrambled-Sobol-init+joint-standardize+ICM-QMC-LBFGS-fit+high-EI-CMA-ES+MFEI-cost-rule-v1",
    );
    push_str(
        &mut bytes,
        "one-post-initialization-low-fidelity-allocation",
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

fn branin(point: &[f64]) -> f64 {
    let x1 = 15.0f64.mul_add(point[0], -5.0);
    let x2 = 15.0 * point[1];
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let t = 1.0 / (8.0 * core::f64::consts::PI);
    let inner = b.mul_add(-(x1 * x1), c.mul_add(x1, x2 - 6.0));
    inner * inner + 10.0 * (1.0 - t) * fs_math::det::cos(x1) + 10.0
}

fn branin_low(point: &[f64]) -> f64 {
    branin(point) + 3.0f64.mul_add(point[0], -2.0 * point[1]) + 5.0
}

fn objective(point: &[f64], fidelity: usize) -> f64 {
    match fidelity {
        0 => branin_low(point),
        1 => branin(point),
        other => panic!("fixture received invalid fidelity {other}"),
    }
}

fn callback_bits(point: &[f64], fidelity: usize, value: f64) -> CallbackBits {
    CallbackBits {
        point: point.iter().map(|value| value.to_bits()).collect(),
        fidelity,
        value: value.to_bits(),
    }
}

fn run_study(config_digest: u64) -> SealedStudy {
    let mut callbacks = Vec::with_capacity(EXPECTED_EVALUATIONS);
    let mut evaluator = |point: &[f64], fidelity: usize| {
        let value = objective(point, fidelity);
        callbacks.push(callback_bits(point, fidelity, value));
        value
    };
    let report = mf_minimize(&mut evaluator, DIMENSION, &config());
    let MfReport {
        f_best_high,
        evals_low,
        evals_high,
        cost,
        trace,
        learned_correlation,
    } = report;
    let record = StudyRecord {
        callbacks,
        f_best_high: f_best_high.to_bits(),
        evals_low,
        evals_high,
        cost: cost.to_bits(),
        trace: trace
            .iter()
            .map(|(cost, best_high)| (cost.to_bits(), best_high.to_bits()))
            .collect(),
        learned_correlation: learned_correlation.to_bits(),
    };
    SealedStudy::seal(config_digest, record)
}

fn independently_reconstructed_correlation(callbacks: &[CallbackBits]) -> f64 {
    let initial = &callbacks[..INITIAL_EVALUATIONS];
    let points: Vec<Vec<f64>> = initial
        .iter()
        .map(|callback| callback.point.iter().copied().map(f64::from_bits).collect())
        .collect();
    let fidelities: Vec<usize> = initial.iter().map(|callback| callback.fidelity).collect();
    let values: Vec<f64> = initial
        .iter()
        .map(|callback| f64::from_bits(callback.value))
        .collect();
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum::<f64>()
        / n;
    let scale = fs_math::det::sqrt(variance.max(1e-30));
    let standardized: Vec<f64> = values.iter().map(|value| (value - mean) / scale).collect();
    let config = config();
    let gp = fit_mf(
        &points,
        &fidelities,
        &standardized,
        config.log_box,
        config.hyper_starts,
        config.seed ^ FIRST_ITERATION_HYPER_SEED_XOR,
    );
    gp.kernel.correlation()
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
    if left.f_best_high != right.f_best_high {
        return Some(format!(
            "report.f_best_high:0x{:016x}!=0x{:016x}",
            left.f_best_high, right.f_best_high
        ));
    }
    if left.evals_low != right.evals_low {
        return Some(format!(
            "report.evals_low:{}!={}",
            left.evals_low, right.evals_low
        ));
    }
    if left.evals_high != right.evals_high {
        return Some(format!(
            "report.evals_high:{}!={}",
            left.evals_high, right.evals_high
        ));
    }
    if left.cost != right.cost {
        return Some(format!(
            "report.cost:0x{:016x}!=0x{:016x}",
            left.cost, right.cost
        ));
    }
    if left.trace.len() != right.trace.len() {
        return Some(format!(
            "report.trace.length:{}!={}",
            left.trace.len(),
            right.trace.len()
        ));
    }
    for (index, (a, b)) in left.trace.iter().zip(&right.trace).enumerate() {
        if a.0 != b.0 {
            return Some(format!(
                "report.trace[{index}].cost:0x{:016x}!=0x{:016x}",
                a.0, b.0
            ));
        }
        if a.1 != b.1 {
            return Some(format!(
                "report.trace[{index}].best_high:0x{:016x}!=0x{:016x}",
                a.1, b.1
            ));
        }
    }
    (left.learned_correlation != right.learned_correlation).then(|| {
        format!(
            "report.learned_correlation:0x{:016x}!=0x{:016x}",
            left.learned_correlation, right.learned_correlation
        )
    })
}

fn mutate_report_trace(reference: &SealedStudy) -> (SealedStudy, Mutation) {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let mutable_trace_count = reference.record.trace.len() - N_INIT_LOW;
    let trace_offset = usize::try_from(selector.next_below(usize_u64(mutable_trace_count)))
        .expect("trace offset fits usize");
    let trace_index = N_INIT_LOW + trace_offset;
    let mantissa_bit =
        u32::try_from(selector.next_below(20)).expect("selected mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut record = reference.record.clone();
    let before = record.trace[trace_index].1;
    let after = before ^ (1u64 << mantissa_bit);
    record.trace[trace_index].1 = after;
    let mutation = Mutation {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
        trace_index,
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
    let Some(&reference_bits) = reference
        .record
        .trace
        .get(mutation.trace_index)
        .map(|entry| &entry.1)
    else {
        return false;
    };
    let Some(&mutant_bits) = mutant
        .record
        .trace
        .get(mutation.trace_index)
        .map(|entry| &entry.1)
    else {
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
    expected.trace[mutation.trace_index].1 = mutation.after;
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
        "fs-bo-mf-seeded-study-mutation-v1",
        config_digest,
        &[reference_digest, mutant_digest],
    );
    push_u64(&mut bytes, mutation.seed);
    push_u64(&mut bytes, u64::from(mutation.kernel));
    push_u64(&mut bytes, u64::from(mutation.tile));
    push_u64(&mut bytes, usize_u64(mutation.trace_index));
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
        "seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=report.trace[{}].best_high; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; reference_output=0x{reference_digest:016x}; mutant_output=0x{:016x}; gate={gate_error:?}; first_mismatch={mismatch}",
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.trace_index,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        mutant.output_digest,
    );
    Suite::new(SUITE)
        .case(
            "seeded-report-trace-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                if matches!(
                    gate_error,
                    Err(AdmissionError::ReferenceIdentityMismatch { .. })
                ) {
                    CaseOutcome::fail(details).with_evidence(
                        "crates/fs-bo/tests/mf_study_replay.rs::seeded-report-trace-corruption",
                    )
                } else {
                    CaseOutcome::pass(format!(
                        "seeded trace corruption escaped the reference gate: {gate_error:?}"
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
fn production_mf_study_replays_and_seeded_corruption_stays_red() {
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

    let (mutant, mutation) = mutate_report_trace(&original);
    let (mutant_replayed, mutation_replayed) = mutate_report_trace(&replay);
    let mutant_mismatch = first_study_mismatch(&original.record, &mutant.record);
    let mutant_replayed_mismatch = first_study_mismatch(&replay.record, &mutant_replayed.record);
    let mutant_accounting = mutant.record.accounting_mismatch();
    let mutant_replayed_accounting = mutant_replayed.record.accounting_mismatch();
    let expected_mismatch = format!("report.trace[{}].best_high", mutation.trace_index);
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
        "fs-bo-mf-seeded-study-mutation-pair-case-v1",
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
        && red_first.records[0].case == "seeded-report-trace-corruption"
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
    let mutated_value = f64::from_bits(mutation.after);
    let mutation_pass = mutant.validate_payload().is_ok()
        && mutant_replayed.validate_payload().is_ok()
        && stale_payload_identity_is_refused(&mutant)
        && mutant == mutant_replayed
        && mutation == mutation_replayed
        && is_exact_report_bit_delta(&original, &mutant, mutation)
        && is_exact_report_bit_delta(&replay, &mutant_replayed, mutation_replayed)
        && (N_INIT_LOW..EXPECTED_EVALUATIONS).contains(&mutation.trace_index)
        && (0..20).contains(&mutation.mantissa_bit)
        && mutation.before != mutation.after
        && f64::from_bits(mutation.before).is_finite()
        && mutated_value.is_finite()
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
        && merge_gate_message.contains("seeded-report-trace-corruption")
        && merge_gate_message.contains("ReferenceIdentityMismatch");

    let accounting_inputs = case_inputs(
        "fs-bo-mf-study-accounting-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let replay_inputs = case_inputs(
        "fs-bo-mf-study-replay-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let accounting_pass = original_accounting.is_none() && replay_accounting.is_none();
    let accounting_detail = format!(
        "config=0x{config_digest:016x}; original=0x{:016x}; replay=0x{:016x}; callbacks={}; evals_low={}; evals_high={}; cost_bits=0x{:016x}; trace={}; learned_correlation_bits=0x{:016x}; original_mismatch={original_accounting:?}; replay_mismatch={replay_accounting:?}",
        original.output_digest,
        replay.output_digest,
        original.record.callbacks.len(),
        original.record.evals_low,
        original.record.evals_high,
        original.record.cost,
        original.record.trace.len(),
        original.record.learned_correlation,
    );
    let replay_detail = format!(
        "config=0x{config_digest:016x}; original=0x{:016x}; replay=0x{:016x}; first_mismatch={replay_mismatch:?}",
        original.output_digest, replay.output_digest
    );
    let mutation_detail = format!(
        "config=0x{config_digest:016x}; reference=0x{:016x}; replay=0x{:016x}; mutant=0x{:016x}; replay_mutant=0x{:016x}; mutation_inputs=0x{mutation_inputs_digest:016x}; replay_mutation_inputs=0x{mutation_replayed_inputs_digest:016x}; seed=0x{MUTATION_SEED:016x}; target=report.trace[{}].best_high; mantissa_bit={}; selector_draws={}; first_mismatch={mutant_mismatch:?}; replay_first_mismatch={mutant_replayed_mismatch:?}; mutant_accounting={mutant_accounting:?}; replay_mutant_accounting={mutant_replayed_accounting:?}; red_record_stable={red_lines_stable}; merge_gate_message={merge_gate_message:?}",
        original.output_digest,
        replay.output_digest,
        mutant.output_digest,
        mutant_replayed.output_digest,
        mutation.trace_index,
        mutation.mantissa_bit,
        mutation.selector_draws,
    );

    let report = Suite::new(SUITE)
        .case(
            "callback-fidelity-cost-and-public-report-accounting",
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
            "same-seed-full-mf-output-frame-replay",
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
            "seeded-resealed-mf-trace-mutation-is-refused",
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
            "callback-fidelity-cost-and-public-report-accounting",
            "same-seed-full-mf-output-frame-replay",
            "seeded-resealed-mf-trace-mutation-is-refused",
        ]
    );
    report.assert_green();
}
