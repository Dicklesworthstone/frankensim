//! G0/G3 replay-complete sparse-GP study receipt and seeded red self-test.
//!
//! This target runs production `SparseGp` across a short inducing-count
//! ladder. It retains every training/probe bit, selected inducing location,
//! ELBO, and public prediction, independently reconstructs farthest-point row
//! selection, and requires an independently executed same-configuration frame
//! to replay exactly. A disclosed mutation changes one finite prediction bit,
//! reseals the payload, and proves that payload identity alone cannot make a
//! semantically altered study authoritative.
//!
//! The fixture establishes same-process repetition and evidence plumbing for
//! this one finite study. It makes no ELBO-theorem, predictive-accuracy,
//! calibration, exact-GP-equivalence, monotonic-quality, optimizer, all-seed,
//! all-configuration, cross-process, cross-ISA, cancellation, persistence,
//! authenticated-admission, private-state, or performance claim.

#![deny(unsafe_code)]

use fs_bo::{Kernel, Matern, SparseGp, farthest_point_inducing};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_rand::StreamKey;

const SUITE: &str = "fs-bo/sparse-gp-study-replay-v1";
const STUDY_CASE: &str = "sparse-gp-short-inducing-ladder";
const STUDY_SEED: u64 = 0x5A12_5EED_0000_0025;
const STUDY_KERNEL: u32 = 0x5A25;
const TRAINING_TILE: u32 = 0;
const PROBE_TILE: u32 = 1;
const MUTATION_SEED: u64 = 0x5A12_5EED_0000_00D3;
const MUTATION_KERNEL: u32 = 0x5AD3;
const MUTATION_TILE: u32 = 0;

const DIMENSION: usize = 2;
const TRAINING_POINTS: usize = 12;
const PROBE_POINTS: usize = 5;
const INDUCING_LADDER: [usize; 3] = [3, 6, TRAINING_POINTS];
const SIGNAL: f64 = 1.25;
const LENGTHSCALES: [f64; DIMENSION] = [0.35, 0.45];
const NOISE: f64 = 1e-3;

const _: () = assert!(DIMENSION == 2);
const _: () = assert!(INDUCING_LADDER[INDUCING_LADDER.len() - 1] == TRAINING_POINTS);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PredictionBits {
    mean: u64,
    variance: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct KernelBits {
    family: u8,
    signal: u64,
    lengthscales: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SparseFrame {
    inducing_count: usize,
    kernel: KernelBits,
    noise: u64,
    inducing_points: Vec<Vec<u64>>,
    elbo: u64,
    predictions: Vec<PredictionBits>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    training_points: Vec<Vec<u64>>,
    training_values: Vec<u64>,
    probe_points: Vec<Vec<u64>>,
    frames: Vec<SparseFrame>,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-bo-sparse-gp-study-output-v1".to_vec();
        push_u64(&mut bytes, config_digest);
        push_matrix(&mut bytes, &self.training_points);
        push_u64_slice(&mut bytes, &self.training_values);
        push_matrix(&mut bytes, &self.probe_points);
        push_len(&mut bytes, self.frames.len());
        for frame in &self.frames {
            push_len(&mut bytes, frame.inducing_count);
            bytes.push(frame.kernel.family);
            push_u64(&mut bytes, frame.kernel.signal);
            push_u64_slice(&mut bytes, &frame.kernel.lengthscales);
            push_u64(&mut bytes, frame.noise);
            push_matrix(&mut bytes, &frame.inducing_points);
            push_u64(&mut bytes, frame.elbo);
            push_len(&mut bytes, frame.predictions.len());
            for prediction in &frame.predictions {
                push_u64(&mut bytes, prediction.mean);
                push_u64(&mut bytes, prediction.variance);
            }
        }
        bytes
    }

    #[allow(clippy::too_many_lines)]
    fn semantic_mismatch(&self) -> Option<String> {
        let expected_training = generated_points(TRAINING_POINTS, TRAINING_TILE);
        let expected_probes = generated_points(PROBE_POINTS, PROBE_TILE);
        let expected_training_bits = matrix_bits(&expected_training);
        let expected_probe_bits = matrix_bits(&expected_probes);
        if self.training_points != expected_training_bits {
            return Some(first_matrix_mismatch(
                "training_points",
                &expected_training_bits,
                &self.training_points,
            ));
        }
        if self.probe_points != expected_probe_bits {
            return Some(first_matrix_mismatch(
                "probe_points",
                &expected_probe_bits,
                &self.probe_points,
            ));
        }
        if self.training_values.len() != TRAINING_POINTS {
            return Some(format!(
                "training_values.length:{}!=expected-{TRAINING_POINTS}",
                self.training_values.len()
            ));
        }
        for (index, (point, &found)) in expected_training
            .iter()
            .zip(&self.training_values)
            .enumerate()
        {
            let expected = target(point).to_bits();
            if found != expected {
                return Some(format!(
                    "training_values[{index}]:0x{found:016x}!=0x{expected:016x}"
                ));
            }
        }
        if self.frames.len() != INDUCING_LADDER.len() {
            return Some(format!(
                "frames.length:{}!=expected-{}",
                self.frames.len(),
                INDUCING_LADDER.len()
            ));
        }

        let training_values: Vec<f64> = self
            .training_values
            .iter()
            .copied()
            .map(f64::from_bits)
            .collect();
        for (frame_index, (&expected_count, frame)) in
            INDUCING_LADDER.iter().zip(&self.frames).enumerate()
        {
            if frame.inducing_count != expected_count {
                return Some(format!(
                    "frames[{frame_index}].inducing_count:{}!=expected-{expected_count}",
                    frame.inducing_count
                ));
            }
            let expected_kernel = kernel_bits(&kernel());
            if frame.kernel != expected_kernel {
                return Some(format!(
                    "frames[{frame_index}].kernel:{:?}!=expected-{expected_kernel:?}",
                    frame.kernel
                ));
            }
            if frame.noise != NOISE.to_bits() {
                return Some(format!(
                    "frames[{frame_index}].noise:0x{:016x}!=expected-0x{:016x}",
                    frame.noise,
                    NOISE.to_bits()
                ));
            }
            let oracle_indices = independent_farthest_indices(&expected_training, expected_count);
            let oracle_points: Vec<Vec<u64>> = oracle_indices
                .iter()
                .map(|&index| expected_training_bits[index].clone())
                .collect();
            if frame.inducing_points != oracle_points {
                return Some(first_matrix_mismatch(
                    &format!("frames[{frame_index}].inducing_points"),
                    &oracle_points,
                    &frame.inducing_points,
                ));
            }
            if !f64::from_bits(frame.elbo).is_finite() {
                return Some(format!("frames[{frame_index}].elbo:not-finite"));
            }
            if frame.predictions.len() != PROBE_POINTS {
                return Some(format!(
                    "frames[{frame_index}].predictions.length:{}!=expected-{PROBE_POINTS}",
                    frame.predictions.len()
                ));
            }

            let inducing_points: Vec<Vec<f64>> = frame
                .inducing_points
                .iter()
                .map(|point| point.iter().copied().map(f64::from_bits).collect())
                .collect();
            let sparse = SparseGp::fit(
                &expected_training,
                &training_values,
                kernel(),
                NOISE,
                inducing_points,
            );
            if kernel_bits(&sparse.kernel) != frame.kernel {
                return Some(format!("frames[{frame_index}].returned_kernel"));
            }
            if sparse.noise.to_bits() != frame.noise {
                return Some(format!("frames[{frame_index}].returned_noise"));
            }
            if matrix_bits(&sparse.z) != frame.inducing_points {
                return Some(format!("frames[{frame_index}].returned_inducing_points"));
            }
            if sparse.elbo.to_bits() != frame.elbo {
                return Some(format!(
                    "frames[{frame_index}].elbo:0x{:016x}!=recomputed-0x{:016x}",
                    frame.elbo,
                    sparse.elbo.to_bits()
                ));
            }
            for (probe_index, (probe, found)) in
                expected_probes.iter().zip(&frame.predictions).enumerate()
            {
                let (mean, variance) = sparse.predict(probe);
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
                        "frames[{frame_index}].predictions[{probe_index}].{field}:0x{found_bits:016x}!=recomputed-0x{expected_bits:016x}"
                    ));
                }
                if !mean.is_finite() || !variance.is_finite() || variance < 0.0 {
                    return Some(format!(
                        "frames[{frame_index}].predictions[{probe_index}]:invalid-domain"
                    ));
                }
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
    frame: usize,
    probe: usize,
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

fn config_bytes() -> Vec<u8> {
    let mut bytes = b"fs-bo-sparse-gp-study-config-v1".to_vec();
    push_str(&mut bytes, STUDY_CASE);
    push_str(&mut bytes, "dimensionless-unit-square-coordinate");
    push_str(&mut bytes, "dimensionless-polynomial-response");
    push_u64(&mut bytes, STUDY_SEED);
    push_u64(&mut bytes, u64::from(STUDY_KERNEL));
    push_u64(&mut bytes, u64::from(TRAINING_TILE));
    push_u64(&mut bytes, u64::from(PROBE_TILE));
    push_len(&mut bytes, DIMENSION);
    push_len(&mut bytes, TRAINING_POINTS);
    push_len(&mut bytes, PROBE_POINTS);
    push_len(&mut bytes, INDUCING_LADDER.len());
    for count in INDUCING_LADDER {
        push_len(&mut bytes, count);
    }
    push_str(&mut bytes, "matern-five-halves");
    push_u64(&mut bytes, SIGNAL.to_bits());
    for lengthscale in LENGTHSCALES {
        push_u64(&mut bytes, lengthscale.to_bits());
    }
    push_u64(&mut bytes, NOISE.to_bits());
    push_str(
        &mut bytes,
        "production-dtc-sor-predict+titsias-elbo+farthest-point-v1",
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
        "no-elbo-theorem-predictive-accuracy-calibration-exact-GP-equivalence-monotonic-quality-optimizer-all-seed-all-config-cross-process-cross-ISA-Cx-persistence-auth-private-state-performance-claim",
    );
    bytes
}

fn kernel() -> Kernel {
    Kernel {
        family: Matern::FiveHalves,
        signal: SIGNAL,
        lengthscales: LENGTHSCALES.to_vec(),
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

fn generated_points(count: usize, tile: u32) -> Vec<Vec<f64>> {
    let mut stream = StreamKey {
        seed: STUDY_SEED,
        kernel: STUDY_KERNEL,
        tile,
    }
    .stream();
    (0..count)
        .map(|_| (0..DIMENSION).map(|_| stream.next_f64()).collect())
        .collect()
}

fn target(point: &[f64]) -> f64 {
    let x = point[0];
    let y = point[1];
    let bowl = (x - 0.35).mul_add(x - 0.35, (y - 0.7) * (y - 0.7));
    (0.125 * x).mul_add(y, bowl + 0.2 * x - 0.1 * y)
}

/// Separate row-index oracle for the public farthest-point rule. Production
/// retains and updates a minimum-distance vector; this version rescans every
/// selected row for every candidate on every round.
fn independent_farthest_indices(points: &[Vec<f64>], count: usize) -> Vec<usize> {
    assert!(!points.is_empty());
    assert!((1..=points.len()).contains(&count));
    let mut chosen = vec![0usize];
    while chosen.len() < count {
        let mut best_index = None;
        let mut best_distance = f64::NEG_INFINITY;
        for (candidate, point) in points.iter().enumerate() {
            if chosen.contains(&candidate) {
                continue;
            }
            let mut nearest = f64::INFINITY;
            for &selected in &chosen {
                let distance = point
                    .iter()
                    .zip(&points[selected])
                    .map(|(left, right)| {
                        let delta = left - right;
                        delta * delta
                    })
                    .sum::<f64>();
                nearest = nearest.min(distance);
            }
            if nearest > best_distance {
                best_index = Some(candidate);
                best_distance = nearest;
            }
        }
        chosen.push(best_index.expect("an unselected fixture row remains"));
    }
    chosen
}

fn run_study(config_digest: u64) -> SealedStudy {
    let training_points = generated_points(TRAINING_POINTS, TRAINING_TILE);
    let training_values: Vec<f64> = training_points.iter().map(|point| target(point)).collect();
    let probe_points = generated_points(PROBE_POINTS, PROBE_TILE);
    let frames = INDUCING_LADDER
        .iter()
        .copied()
        .map(|inducing_count| {
            let inducing_points = farthest_point_inducing(&training_points, inducing_count);
            let sparse = SparseGp::fit(
                &training_points,
                &training_values,
                kernel(),
                NOISE,
                inducing_points,
            );
            let predictions = probe_points
                .iter()
                .map(|probe| {
                    let (mean, variance) = sparse.predict(probe);
                    PredictionBits {
                        mean: mean.to_bits(),
                        variance: variance.to_bits(),
                    }
                })
                .collect();
            SparseFrame {
                inducing_count,
                kernel: kernel_bits(&sparse.kernel),
                noise: sparse.noise.to_bits(),
                inducing_points: matrix_bits(&sparse.z),
                elbo: sparse.elbo.to_bits(),
                predictions,
            }
        })
        .collect();
    let record = StudyRecord {
        training_points: matrix_bits(&training_points),
        training_values: training_values
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        probe_points: matrix_bits(&probe_points),
        frames,
    };
    SealedStudy::seal(config_digest, record)
}

fn first_matrix_mismatch(label: &str, expected: &[Vec<u64>], found: &[Vec<u64>]) -> String {
    if expected.len() != found.len() {
        return format!("{label}.length:{}!={}", found.len(), expected.len());
    }
    for (row_index, (expected_row, found_row)) in expected.iter().zip(found).enumerate() {
        if expected_row.len() != found_row.len() {
            return format!(
                "{label}[{row_index}].length:{}!={}",
                found_row.len(),
                expected_row.len()
            );
        }
        if let Some((column, (&expected_bits, &found_bits))) = expected_row
            .iter()
            .zip(found_row)
            .enumerate()
            .find(|(_, (expected_bits, found_bits))| expected_bits != found_bits)
        {
            return format!(
                "{label}[{row_index}][{column}]:0x{found_bits:016x}!=0x{expected_bits:016x}"
            );
        }
    }
    format!("{label}:different-without-located-element")
}

fn first_record_mismatch(expected: &StudyRecord, found: &StudyRecord) -> Option<String> {
    if expected.training_points != found.training_points {
        return Some(first_matrix_mismatch(
            "training_points",
            &expected.training_points,
            &found.training_points,
        ));
    }
    if expected.training_values != found.training_values {
        return Some("training_values:first-bit-mismatch".to_string());
    }
    if expected.probe_points != found.probe_points {
        return Some(first_matrix_mismatch(
            "probe_points",
            &expected.probe_points,
            &found.probe_points,
        ));
    }
    if expected.frames.len() != found.frames.len() {
        return Some(format!(
            "frames.length:{}!={}",
            found.frames.len(),
            expected.frames.len()
        ));
    }
    for (frame_index, (expected_frame, found_frame)) in
        expected.frames.iter().zip(&found.frames).enumerate()
    {
        if expected_frame.inducing_count != found_frame.inducing_count {
            return Some(format!("frames[{frame_index}].inducing_count"));
        }
        if expected_frame.kernel != found_frame.kernel {
            return Some(format!("frames[{frame_index}].kernel"));
        }
        if expected_frame.noise != found_frame.noise {
            return Some(format!("frames[{frame_index}].noise"));
        }
        if expected_frame.inducing_points != found_frame.inducing_points {
            return Some(first_matrix_mismatch(
                &format!("frames[{frame_index}].inducing_points"),
                &expected_frame.inducing_points,
                &found_frame.inducing_points,
            ));
        }
        if expected_frame.elbo != found_frame.elbo {
            return Some(format!("frames[{frame_index}].elbo"));
        }
        for (probe_index, (expected_prediction, found_prediction)) in expected_frame
            .predictions
            .iter()
            .zip(&found_frame.predictions)
            .enumerate()
        {
            if expected_prediction.mean != found_prediction.mean {
                return Some(format!(
                    "frames[{frame_index}].predictions[{probe_index}].mean"
                ));
            }
            if expected_prediction.variance != found_prediction.variance {
                return Some(format!(
                    "frames[{frame_index}].predictions[{probe_index}].variance"
                ));
            }
        }
        if expected_frame.predictions.len() != found_frame.predictions.len() {
            return Some(format!("frames[{frame_index}].predictions.length"));
        }
    }
    None
}

fn mutate_prediction(reference: &SealedStudy) -> (SealedStudy, Mutation) {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let frame = usize::try_from(selector.next_below(usize_u64(INDUCING_LADDER.len())))
        .expect("frame index fits usize");
    let probe = usize::try_from(selector.next_below(usize_u64(PROBE_POINTS)))
        .expect("probe index fits usize");
    let mantissa_bit =
        u32::try_from(selector.next_below(20)).expect("selected mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut record = reference.record.clone();
    let before = record.frames[frame].predictions[probe].mean;
    let after = before ^ (1u64 << mantissa_bit);
    record.frames[frame].predictions[probe].mean = after;
    let mutation = Mutation {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
        frame,
        probe,
        mantissa_bit,
        selector_draws,
        before,
        after,
    };
    (SealedStudy::seal(reference.config_digest, record), mutation)
}

fn is_exact_prediction_delta(
    reference: &SealedStudy,
    mutant: &SealedStudy,
    mutation: Mutation,
) -> bool {
    if reference.config_digest != mutant.config_digest
        || mutation.before ^ mutation.after != 1u64 << mutation.mantissa_bit
    {
        return false;
    }
    let mut expected = reference.record.clone();
    if expected.frames[mutation.frame].predictions[mutation.probe].mean != mutation.before {
        return false;
    }
    expected.frames[mutation.frame].predictions[mutation.probe].mean = mutation.after;
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
        "fs-bo-sparse-gp-seeded-prediction-mutation-v1",
        config_digest,
        &[reference_digest, mutant_digest],
    );
    push_u64(&mut bytes, mutation.seed);
    push_u64(&mut bytes, u64::from(mutation.kernel));
    push_u64(&mut bytes, u64::from(mutation.tile));
    push_len(&mut bytes, mutation.frame);
    push_len(&mut bytes, mutation.probe);
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
    let semantic_rejected = semantic_mismatch.contains(".mean:");
    let details = format!(
        "seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=frames[{}].predictions[{}].mean; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; reference_output=0x{reference_digest:016x}; mutant_output=0x{:016x}; gate={gate_error:?}; semantic={semantic_mismatch}",
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.frame,
        mutation.probe,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        mutant.output_digest,
    );
    Suite::new(SUITE)
        .case(
            "seeded-sparse-prediction-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                if matches!(
                    gate_error,
                    Err(AdmissionError::ReferenceIdentityMismatch { .. })
                ) && semantic_rejected
                {
                    CaseOutcome::fail(details).with_evidence(
                        "crates/fs-bo/tests/sparse_gp_study_replay.rs::seeded-sparse-prediction-corruption",
                    )
                } else {
                    CaseOutcome::pass(format!(
                        "seeded prediction corruption escaped a gate: identity={gate_error:?}; detail={details}"
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
fn production_sparse_gp_study_replays_and_seeded_prediction_corruption_stays_red() {
    let config_frame = config_bytes();
    let config_digest = fnv1a64(&config_frame);
    let original = run_study(config_digest);
    let replay = run_study(config_digest);

    let original_semantic_mismatch = original.record.semantic_mismatch();
    let replay_semantic_mismatch = replay.record.semantic_mismatch();
    let replay_mismatch = first_record_mismatch(&original.record, &replay.record);
    let semantic_pass = original.validate_payload().is_ok()
        && replay.validate_payload().is_ok()
        && original_semantic_mismatch.is_none()
        && replay_semantic_mismatch.is_none();
    let replay_pass = original.output_digest == replay.output_digest
        && original.record == replay.record
        && replay_mismatch.is_none();

    let (mutant, mutation) = mutate_prediction(&original);
    let (replayed_mutant, replayed_mutation) = mutate_prediction(&replay);
    let mutant_semantic = mutant.record.semantic_mismatch();
    let replayed_mutant_semantic = replayed_mutant.record.semantic_mismatch();
    let mutant_mismatch = first_record_mismatch(&original.record, &mutant.record);
    let expected_path = format!(
        "frames[{}].predictions[{}].mean",
        mutation.frame, mutation.probe
    );
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
        && mutation == replayed_mutation
        && mutant == replayed_mutant
        && is_exact_prediction_delta(&original, &mutant, mutation)
        && mutation.before != mutation.after
        && f64::from_bits(mutation.after).is_finite()
        && mutant.output_digest != original.output_digest
        && matches!(
            mutant.admit_against(original.output_digest),
            Err(AdmissionError::ReferenceIdentityMismatch { expected, found })
                if expected == original.output_digest && found == mutant.output_digest
        )
        && mutant_mismatch.as_deref() == Some(expected_path.as_str())
        && mutant_semantic
            .as_deref()
            .is_some_and(|mismatch| mismatch.starts_with(&expected_path))
        && red_lines_stable
        && !red_first.all_passed()
        && red_first.records[0]
            .details
            .contains("ReferenceIdentityMismatch")
        && red_first.records[0].details.contains(&expected_path)
        && merge_gate_message.contains("seeded-sparse-prediction-corruption")
        && merge_gate_message.contains("ReferenceIdentityMismatch");

    let semantic_inputs = case_inputs(
        "fs-bo-sparse-gp-study-semantic-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let replay_inputs = case_inputs(
        "fs-bo-sparse-gp-study-replay-case-v1",
        config_digest,
        &[original.output_digest, replay.output_digest],
    );
    let mutation_case_inputs = case_inputs(
        "fs-bo-sparse-gp-study-mutation-case-v1",
        config_digest,
        &[mutation_inputs_digest, replayed_mutation_inputs_digest],
    );
    let semantic_detail = format!(
        "config=0x{config_digest:016x}; output=0x{:016x}; replay=0x{:016x}; training={}; probes={}; frames={}; inducing={INDUCING_LADDER:?}; original_mismatch={original_semantic_mismatch:?}; replay_mismatch={replay_semantic_mismatch:?}",
        original.output_digest,
        replay.output_digest,
        original.record.training_points.len(),
        original.record.probe_points.len(),
        original.record.frames.len(),
    );
    let replay_detail = format!(
        "config=0x{config_digest:016x}; output=0x{:016x}; replay=0x{:016x}; first_mismatch={replay_mismatch:?}",
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
            "complete-public-frame-and-independent-selection-audit",
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
            "same-configuration-full-frame-replay",
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
            "seeded-resealed-prediction-mutation-is-refused",
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
            "complete-public-frame-and-independent-selection-audit",
            "same-configuration-full-frame-replay",
            "seeded-resealed-prediction-mutation-is-refused",
        ]
    );
    report.assert_green();
}
