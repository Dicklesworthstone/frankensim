//! G0/G3 replay-complete heteroscedastic-GP state receipt and red self-test.
//!
//! This target fits production `Gp::try_fit_diag` on one deterministic finite
//! study and retains all public GP fields, the complete public training view,
//! every study input, and every probe prediction. It independently regenerates
//! the stream-derived study inputs and requires a fresh same-process fit to
//! reproduce the complete frame. A disclosed mutation changes one finite alpha
//! bit, reseals the payload, and proves that payload identity alone cannot make
//! altered continuation state authoritative.
//!
//! The fixture makes no posterior-accuracy, calibration, independent-numerical-
//! equivalence, optimizer-quality, all-seed, all-configuration, cross-process,
//! cross-ISA, cancellation, persistence, authenticated-admission, or
//! performance claim.

#![deny(unsafe_code)]

use fs_bo::{Gp, Kernel, Matern};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_rand::StreamKey;

const SUITE: &str = "fs-bo/hetero-gp-study-replay-v1";
const STUDY_CASE: &str = "heteroscedastic-gp-finite-state";
const STUDY_SEED: u64 = 0x4E73_5EED_0000_0026;
const STUDY_KERNEL: u32 = 0x4E26;
const TRAINING_TILE: u32 = 0;
const PROBE_TILE: u32 = 1;
const MUTATION_SEED: u64 = 0x4E73_5EED_0000_00D4;
const MUTATION_KERNEL: u32 = 0x4ED4;
const MUTATION_TILE: u32 = 0;

const DIMENSION: usize = 1;
const TRAINING_POINTS: usize = 10;
const PROBE_POINTS: usize = 5;
const SIGNAL: f64 = 1.1;
const LENGTHSCALE: f64 = 0.3;
const LOW_NOISE_VARIANCE: f64 = 4e-4;
const HIGH_NOISE_VARIANCE: f64 = 4e-2;
const OBSERVATION_PERTURBATION: f64 = 0.03;

const _: () = assert!(DIMENSION == 1 && TRAINING_POINTS > PROBE_POINTS);

#[derive(Debug, Clone, PartialEq, Eq)]
struct KernelBits {
    family: u8,
    signal: u64,
    lengthscales: Vec<u64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PredictionBits {
    mean: u64,
    variance: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    training_points: Vec<Vec<u64>>,
    observations: Vec<u64>,
    noise_variances: Vec<u64>,
    probe_points: Vec<Vec<u64>>,
    returned_kernel: KernelBits,
    reported_noise: u64,
    training_view_points: Vec<Vec<u64>>,
    training_view_cholesky: Vec<u64>,
    training_view_alpha: Vec<u64>,
    lml: u64,
    predictions: Vec<PredictionBits>,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-bo-hetero-gp-study-output-v1".to_vec();
        push_u64(&mut bytes, config_digest);
        push_matrix(&mut bytes, &self.training_points);
        push_u64_slice(&mut bytes, &self.observations);
        push_u64_slice(&mut bytes, &self.noise_variances);
        push_matrix(&mut bytes, &self.probe_points);
        bytes.push(self.returned_kernel.family);
        push_u64(&mut bytes, self.returned_kernel.signal);
        push_u64_slice(&mut bytes, &self.returned_kernel.lengthscales);
        push_u64(&mut bytes, self.reported_noise);
        push_matrix(&mut bytes, &self.training_view_points);
        push_u64_slice(&mut bytes, &self.training_view_cholesky);
        push_u64_slice(&mut bytes, &self.training_view_alpha);
        push_u64(&mut bytes, self.lml);
        push_len(&mut bytes, self.predictions.len());
        for prediction in &self.predictions {
            push_u64(&mut bytes, prediction.mean);
            push_u64(&mut bytes, prediction.variance);
        }
        bytes
    }

    #[allow(clippy::too_many_lines)]
    fn semantic_mismatch(&self) -> Option<String> {
        let generated = generated_study();
        let expected_training_bits = matrix_bits(&generated.training_points);
        let expected_probe_bits = matrix_bits(&generated.probe_points);
        let expected_observations: Vec<u64> = generated
            .observations
            .iter()
            .map(|value| value.to_bits())
            .collect();
        let expected_noises: Vec<u64> = generated
            .noise_variances
            .iter()
            .map(|value| value.to_bits())
            .collect();
        if self.training_points != expected_training_bits {
            return Some(first_matrix_mismatch(
                "training_points",
                &expected_training_bits,
                &self.training_points,
            ));
        }
        if self.observations != expected_observations {
            return Some(first_slice_mismatch(
                "observations",
                &expected_observations,
                &self.observations,
            ));
        }
        if self.noise_variances != expected_noises {
            return Some(first_slice_mismatch(
                "noise_variances",
                &expected_noises,
                &self.noise_variances,
            ));
        }
        if self.probe_points != expected_probe_bits {
            return Some(first_matrix_mismatch(
                "probe_points",
                &expected_probe_bits,
                &self.probe_points,
            ));
        }
        let expected_kernel = kernel_bits(&kernel());
        if self.returned_kernel != expected_kernel {
            return Some(format!(
                "returned_kernel:{:?}!=expected-{expected_kernel:?}",
                self.returned_kernel
            ));
        }
        let minimum_noise = generated
            .noise_variances
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        if self.reported_noise != minimum_noise.to_bits() {
            return Some(format!(
                "reported_noise:0x{:016x}!=expected-0x{:016x}",
                self.reported_noise,
                minimum_noise.to_bits()
            ));
        }

        let gp = Gp::try_fit_diag(
            &generated.training_points,
            &generated.observations,
            kernel(),
            &generated.noise_variances,
        )
        .expect("finite fixture covariance is positive definite");
        if kernel_bits(&gp.kernel) != self.returned_kernel {
            return Some("returned_kernel:fresh-fit-mismatch".to_string());
        }
        if gp.noise.to_bits() != self.reported_noise {
            return Some("reported_noise:fresh-fit-mismatch".to_string());
        }
        if gp.lml.to_bits() != self.lml {
            return Some(format!(
                "lml:0x{:016x}!=fresh-fit-0x{:016x}",
                self.lml,
                gp.lml.to_bits()
            ));
        }
        let (view_points, cholesky, alpha) = gp.training_view();
        let expected_view_points = matrix_bits(view_points);
        let expected_cholesky: Vec<u64> = cholesky.iter().map(|value| value.to_bits()).collect();
        let expected_alpha: Vec<u64> = alpha.iter().map(|value| value.to_bits()).collect();
        if self.training_view_points != expected_view_points {
            return Some(first_matrix_mismatch(
                "training_view.points",
                &expected_view_points,
                &self.training_view_points,
            ));
        }
        if self.training_view_points != self.training_points {
            return Some(first_matrix_mismatch(
                "training_view.points-vs-input",
                &self.training_points,
                &self.training_view_points,
            ));
        }
        if self.training_view_cholesky != expected_cholesky {
            return Some(first_slice_mismatch(
                "training_view.cholesky",
                &expected_cholesky,
                &self.training_view_cholesky,
            ));
        }
        if self.training_view_alpha != expected_alpha {
            return Some(first_slice_mismatch(
                "training_view.alpha",
                &expected_alpha,
                &self.training_view_alpha,
            ));
        }
        if self.training_view_cholesky.len() != TRAINING_POINTS * TRAINING_POINTS {
            return Some(format!(
                "training_view.cholesky.length:{}!=expected-{}",
                self.training_view_cholesky.len(),
                TRAINING_POINTS * TRAINING_POINTS
            ));
        }
        if self.training_view_alpha.len() != TRAINING_POINTS {
            return Some(format!(
                "training_view.alpha.length:{}!=expected-{TRAINING_POINTS}",
                self.training_view_alpha.len()
            ));
        }
        if self.predictions.len() != PROBE_POINTS {
            return Some(format!(
                "predictions.length:{}!=expected-{PROBE_POINTS}",
                self.predictions.len()
            ));
        }
        for (probe_index, (probe, found)) in generated
            .probe_points
            .iter()
            .zip(&self.predictions)
            .enumerate()
        {
            let (mean, variance) = gp.predict(probe);
            let expected = PredictionBits {
                mean: mean.to_bits(),
                variance: variance.to_bits(),
            };
            if *found != expected {
                let field = if found.mean != expected.mean {
                    "mean"
                } else {
                    "variance"
                };
                let (found_bits, expected_bits) = if field == "mean" {
                    (found.mean, expected.mean)
                } else {
                    (found.variance, expected.variance)
                };
                return Some(format!(
                    "predictions[{probe_index}].{field}:0x{found_bits:016x}!=fresh-fit-0x{expected_bits:016x}"
                ));
            }
            if !mean.is_finite() || !variance.is_finite() || variance < 0.0 {
                return Some(format!("predictions[{probe_index}]:invalid-domain"));
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

struct GeneratedStudy {
    training_points: Vec<Vec<f64>>,
    observations: Vec<f64>,
    noise_variances: Vec<f64>,
    probe_points: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    alpha_index: usize,
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

fn push_matrix(bytes: &mut Vec<u8>, matrix: &[Vec<u64>]) {
    push_len(bytes, matrix.len());
    for row in matrix {
        push_u64_slice(bytes, row);
    }
}

fn matrix_bits(matrix: &[Vec<f64>]) -> Vec<Vec<u64>> {
    matrix
        .iter()
        .map(|row| row.iter().map(|value| value.to_bits()).collect())
        .collect()
}

fn first_slice_mismatch(label: &str, expected: &[u64], found: &[u64]) -> String {
    if expected.len() != found.len() {
        return format!("{label}.length:{}!={}", found.len(), expected.len());
    }
    if let Some((index, (&expected_bits, &found_bits))) = expected
        .iter()
        .zip(found)
        .enumerate()
        .find(|(_, (expected_bits, found_bits))| expected_bits != found_bits)
    {
        format!("{label}[{index}]:0x{found_bits:016x}!=0x{expected_bits:016x}")
    } else {
        format!("{label}:different-without-located-element")
    }
}

fn first_matrix_mismatch(label: &str, expected: &[Vec<u64>], found: &[Vec<u64>]) -> String {
    if expected.len() != found.len() {
        return format!("{label}.length:{}!={}", found.len(), expected.len());
    }
    for (row_index, (expected_row, found_row)) in expected.iter().zip(found).enumerate() {
        if expected_row != found_row {
            return first_slice_mismatch(&format!("{label}[{row_index}]"), expected_row, found_row);
        }
    }
    format!("{label}:different-without-located-row")
}

fn config_bytes() -> Vec<u8> {
    let mut bytes = b"fs-bo-hetero-gp-study-config-v1".to_vec();
    push_str(&mut bytes, STUDY_CASE);
    push_str(&mut bytes, "dimensionless-unit-interval-coordinate");
    push_str(
        &mut bytes,
        "dimensionless-polynomial-plus-seeded-observation",
    );
    push_u64(&mut bytes, STUDY_SEED);
    push_u64(&mut bytes, u64::from(STUDY_KERNEL));
    push_u64(&mut bytes, u64::from(TRAINING_TILE));
    push_u64(&mut bytes, u64::from(PROBE_TILE));
    push_len(&mut bytes, DIMENSION);
    push_len(&mut bytes, TRAINING_POINTS);
    push_len(&mut bytes, PROBE_POINTS);
    push_str(&mut bytes, "matern-five-halves");
    push_u64(&mut bytes, SIGNAL.to_bits());
    push_u64(&mut bytes, LENGTHSCALE.to_bits());
    push_u64(&mut bytes, LOW_NOISE_VARIANCE.to_bits());
    push_u64(&mut bytes, HIGH_NOISE_VARIANCE.to_bits());
    push_u64(&mut bytes, OBSERVATION_PERTURBATION.to_bits());
    push_str(
        &mut bytes,
        "try-fit-diag+public-training-view+latent-predict-v1",
    );
    push_str(&mut bytes, fs_bo::VERSION);
    push_str(&mut bytes, fs_la::VERSION);
    push_str(&mut bytes, fs_math::VERSION);
    push_str(&mut bytes, fs_rand::VERSION);
    push_u64(&mut bytes, u64::from(fs_rand::STREAM_SEMANTICS_VERSION));
    push_str(&mut bytes, fs_rand::STREAM_POSITION_IDENTITY_DOMAIN);
    push_u64(&mut bytes, u64::from(CASEBOOK_RECORD_VERSION));
    push_str(
        &mut bytes,
        "no-posterior-accuracy-calibration-independent-numerical-equivalence-optimizer-all-seed-all-config-cross-process-cross-ISA-Cx-persistence-auth-performance-claim",
    );
    bytes
}

fn kernel() -> Kernel {
    Kernel {
        family: Matern::FiveHalves,
        signal: SIGNAL,
        lengthscales: vec![LENGTHSCALE],
    }
}

fn kernel_bits(kernel: &Kernel) -> KernelBits {
    let family = match kernel.family {
        Matern::Half => 1,
        Matern::ThreeHalves => 3,
        Matern::FiveHalves => 5,
    };
    KernelBits {
        family,
        signal: kernel.signal.to_bits(),
        lengthscales: kernel
            .lengthscales
            .iter()
            .map(|value| value.to_bits())
            .collect(),
    }
}

fn target(x: f64) -> f64 {
    let centered = x - 0.6;
    centered.mul_add(centered, 0.15 * x - 0.05)
}

fn generated_study() -> GeneratedStudy {
    let mut training = StreamKey {
        seed: STUDY_SEED,
        kernel: STUDY_KERNEL,
        tile: TRAINING_TILE,
    }
    .stream();
    let mut training_points = Vec::with_capacity(TRAINING_POINTS);
    let mut observations = Vec::with_capacity(TRAINING_POINTS);
    let mut noise_variances = Vec::with_capacity(TRAINING_POINTS);
    for index in 0..TRAINING_POINTS {
        let x = training.next_f64();
        let observation = OBSERVATION_PERTURBATION.mul_add(training.next_normal(), target(x));
        let noise = if index.is_multiple_of(3) {
            HIGH_NOISE_VARIANCE
        } else {
            LOW_NOISE_VARIANCE
        };
        training_points.push(vec![x]);
        observations.push(observation);
        noise_variances.push(noise);
    }
    let mut probes = StreamKey {
        seed: STUDY_SEED,
        kernel: STUDY_KERNEL,
        tile: PROBE_TILE,
    }
    .stream();
    let probe_points = (0..PROBE_POINTS).map(|_| vec![probes.next_f64()]).collect();
    GeneratedStudy {
        training_points,
        observations,
        noise_variances,
        probe_points,
    }
}

fn run_study(config_digest: u64) -> SealedStudy {
    let generated = generated_study();
    let gp = Gp::try_fit_diag(
        &generated.training_points,
        &generated.observations,
        kernel(),
        &generated.noise_variances,
    )
    .expect("finite fixture covariance is positive definite");
    let (view_points, cholesky, alpha) = gp.training_view();
    let predictions = generated
        .probe_points
        .iter()
        .map(|probe| {
            let (mean, variance) = gp.predict(probe);
            PredictionBits {
                mean: mean.to_bits(),
                variance: variance.to_bits(),
            }
        })
        .collect();
    let record = StudyRecord {
        training_points: matrix_bits(&generated.training_points),
        observations: generated
            .observations
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        noise_variances: generated
            .noise_variances
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        probe_points: matrix_bits(&generated.probe_points),
        returned_kernel: kernel_bits(&gp.kernel),
        reported_noise: gp.noise.to_bits(),
        training_view_points: matrix_bits(view_points),
        training_view_cholesky: cholesky.iter().map(|value| value.to_bits()).collect(),
        training_view_alpha: alpha.iter().map(|value| value.to_bits()).collect(),
        lml: gp.lml.to_bits(),
        predictions,
    };
    SealedStudy::seal(config_digest, record)
}

fn mutate_alpha(reference: &SealedStudy) -> (SealedStudy, Mutation) {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let alpha_index = usize::try_from(selector.next_below(usize_u64(TRAINING_POINTS)))
        .expect("alpha index fits usize");
    let mantissa_bit =
        u32::try_from(selector.next_below(20)).expect("selected mantissa bit fits u32");
    let selector_draws = selector.index();
    let mut record = reference.record.clone();
    let before = record.training_view_alpha[alpha_index];
    let after = before ^ (1u64 << mantissa_bit);
    record.training_view_alpha[alpha_index] = after;
    let mutation = Mutation {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
        alpha_index,
        mantissa_bit,
        selector_draws,
        before,
        after,
    };
    (SealedStudy::seal(reference.config_digest, record), mutation)
}

fn is_exact_alpha_delta(reference: &SealedStudy, mutant: &SealedStudy, mutation: Mutation) -> bool {
    if reference.config_digest != mutant.config_digest
        || mutation.before ^ mutation.after != 1u64 << mutation.mantissa_bit
    {
        return false;
    }
    let mut expected = reference.record.clone();
    if expected.training_view_alpha[mutation.alpha_index] != mutation.before {
        return false;
    }
    expected.training_view_alpha[mutation.alpha_index] = mutation.after;
    expected == mutant.record
}

fn stale_payload_identity_is_refused(reference: &SealedStudy, mutant: &SealedStudy) -> bool {
    let mut stale = mutant.clone();
    stale.output_digest = reference.output_digest;
    matches!(
        stale.validate_payload(),
        Err(AdmissionError::PayloadIdentityMismatch { declared, computed })
            if declared == reference.output_digest && computed == mutant.output_digest
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
        "fs-bo-hetero-gp-seeded-alpha-mutation-v1",
        config_digest,
        &[reference_digest, mutant_digest],
    );
    push_u64(&mut bytes, mutation.seed);
    push_u64(&mut bytes, u64::from(mutation.kernel));
    push_u64(&mut bytes, u64::from(mutation.tile));
    push_len(&mut bytes, mutation.alpha_index);
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
    semantic_mismatch: &str,
) -> fs_casebook::SuiteReport {
    let gate_error = mutant.admit_against(reference_digest);
    let semantic_rejected = semantic_mismatch.contains("training_view.alpha[");
    let details = format!(
        "seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=training_view.alpha[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; reference_output=0x{reference_digest:016x}; mutant_output=0x{:016x}; gate={gate_error:?}; semantic={semantic_mismatch}",
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.alpha_index,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        mutant.output_digest,
    );
    Suite::new(SUITE)
        .case(
            "seeded-alpha-state-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                if matches!(
                    gate_error,
                    Err(AdmissionError::ReferenceIdentityMismatch { .. })
                ) && semantic_rejected
                {
                    CaseOutcome::fail(details).with_evidence(
                        "crates/fs-bo/tests/hetero_gp_study_replay.rs::seeded-alpha-state-corruption",
                    )
                } else {
                    CaseOutcome::pass(format!(
                        "seeded alpha corruption escaped a gate: identity={gate_error:?}; detail={details}"
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
fn heteroscedastic_gp_state_replays_and_seeded_alpha_corruption_stays_red() {
    let config_frame = config_bytes();
    let config_digest = fnv1a64(&config_frame);
    let original = run_study(config_digest);
    let replay = run_study(config_digest);
    let original_semantic = original.record.semantic_mismatch();
    let replay_semantic = replay.record.semantic_mismatch();
    let semantic_pass = original.validate_payload().is_ok()
        && replay.validate_payload().is_ok()
        && original_semantic.is_none()
        && replay_semantic.is_none();
    let replay_pass = original == replay;

    let (mutant, mutation) = mutate_alpha(&original);
    let (replayed_mutant, replayed_mutation) = mutate_alpha(&replay);
    let mutant_semantic = mutant.record.semantic_mismatch();
    let replayed_mutant_semantic = replayed_mutant.record.semantic_mismatch();
    let expected_path = format!("training_view.alpha[{}]", mutation.alpha_index);
    let mutation_frame = mutation_inputs(
        config_digest,
        original.output_digest,
        mutant.output_digest,
        mutation,
    );
    let replayed_mutation_frame = mutation_inputs(
        config_digest,
        replay.output_digest,
        replayed_mutant.output_digest,
        replayed_mutation,
    );
    let mutation_inputs_digest = fnv1a64(&mutation_frame);
    let replayed_mutation_inputs_digest = fnv1a64(&replayed_mutation_frame);
    let red_first = mutation_red_report(
        mutation_inputs_digest,
        original.output_digest,
        &mutant,
        mutation,
        mutant_semantic.as_deref().unwrap_or("none"),
    );
    let red_second = mutation_red_report(
        replayed_mutation_inputs_digest,
        replay.output_digest,
        &replayed_mutant,
        replayed_mutation,
        replayed_mutant_semantic.as_deref().unwrap_or("none"),
    );
    let red_lines_stable = red_first.records.len() == 1
        && red_second.records.len() == 1
        && mutation_frame == replayed_mutation_frame
        && red_first.records[0].json_line() == red_second.records[0].json_line();
    let merge_gate_panic = std::panic::catch_unwind(|| red_first.assert_green());
    let merge_gate_message = merge_gate_panic
        .as_ref()
        .err()
        .map(|payload| panic_message(payload.as_ref()))
        .unwrap_or_default();
    let mutation_pass = mutant.validate_payload().is_ok()
        && replayed_mutant.validate_payload().is_ok()
        && stale_payload_identity_is_refused(&original, &mutant)
        && mutant == replayed_mutant
        && mutation == replayed_mutation
        && is_exact_alpha_delta(&original, &mutant, mutation)
        && mutation.before != mutation.after
        && f64::from_bits(mutation.after).is_finite()
        && mutant.output_digest != original.output_digest
        && matches!(
            mutant.admit_against(original.output_digest),
            Err(AdmissionError::ReferenceIdentityMismatch { expected, found })
                if expected == original.output_digest && found == mutant.output_digest
        )
        && mutant_semantic
            .as_deref()
            .is_some_and(|mismatch| mismatch.starts_with(&expected_path))
        && red_lines_stable
        && !red_first.all_passed()
        && red_first.records[0]
            .details
            .contains("ReferenceIdentityMismatch")
        && red_first.records[0].details.contains(&expected_path)
        && merge_gate_message.contains("seeded-alpha-state-corruption")
        && merge_gate_message.contains("ReferenceIdentityMismatch");

    let semantic_inputs = case_inputs(
        "fs-bo-hetero-gp-study-semantic-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let replay_inputs = case_inputs(
        "fs-bo-hetero-gp-study-replay-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let mutation_case_inputs = case_inputs(
        "fs-bo-hetero-gp-study-mutation-case-v1",
        config_digest,
        &[mutation_inputs_digest, replayed_mutation_inputs_digest],
    );
    let semantic_detail = format!(
        "config=0x{config_digest:016x}; output=0x{:016x}; replay=0x{:016x}; training={}; probes={}; cholesky={}; alpha={}; original_mismatch={original_semantic:?}; replay_mismatch={replay_semantic:?}",
        original.output_digest,
        replay.output_digest,
        original.record.training_points.len(),
        original.record.probe_points.len(),
        original.record.training_view_cholesky.len(),
        original.record.training_view_alpha.len(),
    );
    let replay_detail = format!(
        "config=0x{config_digest:016x}; output=0x{:016x}; replay=0x{:016x}; equal={replay_pass}",
        original.output_digest, replay.output_digest,
    );
    let mutation_detail = format!(
        "config=0x{config_digest:016x}; reference=0x{:016x}; mutant=0x{:016x}; target={expected_path}; mantissa_bit={}; selector_draws={}; semantic={mutant_semantic:?}; red_record_stable={red_lines_stable}; merge_gate={merge_gate_message:?}",
        original.output_digest,
        mutant.output_digest,
        mutation.mantissa_bit,
        mutation.selector_draws,
    );

    let report = Suite::new(SUITE)
        .case(
            "complete-input-public-and-training-view-state-audit",
            fnv1a64(&semantic_inputs),
            ToleranceSpec::Exact,
            move || {
                if semantic_pass {
                    CaseOutcome::pass(semantic_detail)
                } else {
                    CaseOutcome::fail(semantic_detail)
                }
                .with_evidence("crates/fs-bo/CONTRACT.md#conformance-tests")
            },
        )
        .case(
            "same-process-full-state-frame-repetition",
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
            "seeded-resealed-alpha-state-mutation-is-refused",
            fnv1a64(&mutation_case_inputs),
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
            "complete-input-public-and-training-view-state-audit",
            "same-process-full-state-frame-repetition",
            "seeded-resealed-alpha-state-mutation-is-refused",
        ]
    );
    report.assert_green();
}
