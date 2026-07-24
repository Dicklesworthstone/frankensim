//! G0/G3 full-study replay for the production NSGA-III path.
//!
//! This target retains every ordered objective callback and every returned
//! first-front bit for a short three-objective study. It independently
//! reconstructs the seeded initializer, recomputes every objective, checks
//! exact evaluation accounting, and requires a fresh same-process run to
//! reproduce the complete frame. A disclosed mutation changes one returned
//! objective bit and must fail stale-payload, retained-reference, semantic,
//! stable-red-record, and merge gates.
//!
//! This target makes no convergence, reference-direction coverage,
//! hypervolume, diversity, optimizer-superiority, all-seed, all-configuration,
//! cross-process, cross-ISA, cancellation, persistence, authenticated-
//! admission, internal-selection-history, or performance claim.

#![deny(unsafe_code)]

use fs_dfo::{NsgaParams, das_dennis, nsga3};
use fs_obs::ident::ReplayIdentity;
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/nsga3-study-replay-v3";
const CASE: &str = "short-three-objective-complete-callback-and-front";
const RED_CASE: &str = "seeded-returned-objective-corruption";

const DIMENSION: usize = 4;
const OBJECTIVES: usize = 3;
const DIRECTION_DIVISIONS: usize = 3;
const EXPECTED_DIRECTIONS: usize = 10;
const POPULATION: usize = 12;
const GENERATIONS: usize = 2;
const ETA_C: f64 = 20.0;
const ETA_M: f64 = 25.0;
const MUTATION_PROBABILITY: f64 = 0.25;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND: f64 = 1.0;
const STUDY_SEED: u64 = 0xA3A3_0000_0000_0039;
const OPTIMIZER_STREAM_KERNEL: u32 = 0x05A3;
const OPTIMIZER_STREAM_TILE: u32 = 0;
const EXPECTED_EVALUATIONS: usize = POPULATION * (GENERATIONS + 1);

const MUTATION_SEED: u64 = 0xD0F0_7E1D_0000_0039;
const MUTATION_KERNEL: u32 = 0xD039;
const MUTATION_TILE: u32 = 0;

const _: () = assert!(DIMENSION > 1 && OBJECTIVES == 3);
const _: () = assert!(POPULATION.is_multiple_of(2));
const _: () = assert!(EXPECTED_EVALUATIONS == 36);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    decision: Vec<u64>,
    objectives: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    evaluations: Vec<EvaluationBits>,
    front: Vec<EvaluationBits>,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-dfo-nsga3-study-output-v2".to_vec();
        push_u64(&mut bytes, config_digest);
        push_len(&mut bytes, self.evaluations.len());
        for evaluation in &self.evaluations {
            push_u64_slice(&mut bytes, &evaluation.decision);
            push_u64_slice(&mut bytes, &evaluation.objectives);
        }
        push_len(&mut bytes, self.front.len());
        for individual in &self.front {
            push_u64_slice(&mut bytes, &individual.decision);
            push_u64_slice(&mut bytes, &individual.objectives);
        }
        bytes
    }

    #[allow(clippy::too_many_lines)]
    fn semantic_mismatch(&self) -> Option<String> {
        if self.evaluations.len() != EXPECTED_EVALUATIONS {
            return Some(format!(
                "evaluation-count:{}!=expected-{EXPECTED_EVALUATIONS}",
                self.evaluations.len()
            ));
        }
        let mut initializer = StreamKey {
            seed: STUDY_SEED,
            kernel: OPTIMIZER_STREAM_KERNEL,
            tile: OPTIMIZER_STREAM_TILE,
        }
        .stream();
        for (index, evaluation) in self.evaluations.iter().enumerate() {
            if evaluation.decision.len() != DIMENSION {
                return Some(format!(
                    "evaluations[{index}].x.length:{}!=expected-{DIMENSION}",
                    evaluation.decision.len()
                ));
            }
            if evaluation.objectives.len() != OBJECTIVES {
                return Some(format!(
                    "evaluations[{index}].f.length:{}!=expected-{OBJECTIVES}",
                    evaluation.objectives.len()
                ));
            }
            let decision = decode(&evaluation.decision);
            if decision
                .iter()
                .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
            {
                return Some(format!("evaluations[{index}].x:outside-box"));
            }
            let expected = tri_objective(&decision);
            let expected_bits: Vec<u64> = expected.iter().map(|value| value.to_bits()).collect();
            if evaluation.objectives != expected_bits {
                return Some(format!("evaluations[{index}].f"));
            }
            if index < POPULATION {
                for coordinate in 0..DIMENSION {
                    let expected = (UPPER_BOUND - LOWER_BOUND)
                        .mul_add(initializer.next_f64(), LOWER_BOUND)
                        .to_bits();
                    if evaluation.decision[coordinate] != expected {
                        return Some(format!(
                            "initial_population[{index}].x[{coordinate}]:0x{:016x}!=0x{expected:016x}",
                            evaluation.decision[coordinate]
                        ));
                    }
                }
            }
        }
        if self.front.is_empty() || self.front.len() > POPULATION {
            return Some(format!(
                "front-count:{} outside 1..={POPULATION}",
                self.front.len()
            ));
        }
        for (individual, member) in self.front.iter().enumerate() {
            if member.decision.len() != DIMENSION || member.objectives.len() != OBJECTIVES {
                return Some(format!("front[{individual}].dimensions"));
            }
            let decision = decode(&member.decision);
            if decision
                .iter()
                .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
            {
                return Some(format!("front[{individual}].x:outside-box"));
            }
            let expected = tri_objective(&decision);
            let expected_bits: Vec<u64> = expected.iter().map(|value| value.to_bits()).collect();
            if member.objectives != expected_bits {
                return Some(format!("front[{individual}].f"));
            }
            if !self
                .evaluations
                .iter()
                .any(|evaluation| evaluation == member)
            {
                return Some(format!("front[{individual}]:missing-callback"));
            }
            for (other, candidate_dominator) in self.front.iter().enumerate() {
                if individual != other
                    && independent_dominates(&candidate_dominator.objectives, &member.objectives)
                {
                    return Some(format!("front[{individual}]:dominated-by-front[{other}]"));
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
    individual: usize,
    objective: usize,
    mantissa_bit: u32,
    selector_draws: u64,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CorruptionRun {
    mutant: SealedStudy,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    first_mismatch: String,
    semantic_mismatch: String,
    red_line: String,
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

fn push_policy_identity(bytes: &mut Vec<u8>, field: &str, identity: &ReplayIdentity) {
    push_str(bytes, field);
    push_u64(bytes, u64::from(identity.version()));
    push_u64(bytes, identity.root());
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

fn parameters() -> NsgaParams {
    NsgaParams {
        pop: POPULATION,
        generations: GENERATIONS,
        eta_c: ETA_C,
        eta_m: ETA_M,
        p_mut: MUTATION_PROBABILITY,
        seed: STUDY_SEED,
    }
}

fn directions() -> Vec<Vec<f64>> {
    das_dennis(OBJECTIVES, DIRECTION_DIVISIONS)
}

fn config_bytes_with_policies(normalization: &ReplayIdentity, mating: &ReplayIdentity) -> Vec<u8> {
    let directions = directions();
    assert_eq!(directions.len(), EXPECTED_DIRECTIONS);
    let mut bytes = b"fs-dfo-nsga3-study-config-v3".to_vec();
    push_str(&mut bytes, CASE);
    push_str(&mut bytes, "fs_dfo::nsga3");
    push_policy_identity(&mut bytes, "normalization-policy", normalization);
    push_policy_identity(&mut bytes, "mating-policy", mating);
    push_str(&mut bytes, "three-objective-polynomial-tradeoff-v1");
    push_str(&mut bytes, "dimensionless");
    push_len(&mut bytes, DIMENSION);
    push_len(&mut bytes, OBJECTIVES);
    push_len(&mut bytes, DIRECTION_DIVISIONS);
    push_len(&mut bytes, directions.len());
    for direction in &directions {
        push_len(&mut bytes, direction.len());
        for value in direction {
            push_u64(&mut bytes, value.to_bits());
        }
    }
    push_len(&mut bytes, POPULATION);
    push_len(&mut bytes, GENERATIONS);
    push_u64(&mut bytes, ETA_C.to_bits());
    push_u64(&mut bytes, ETA_M.to_bits());
    push_u64(&mut bytes, MUTATION_PROBABILITY.to_bits());
    push_u64(&mut bytes, LOWER_BOUND.to_bits());
    push_u64(&mut bytes, UPPER_BOUND.to_bits());
    push_u64(&mut bytes, STUDY_SEED);
    push_u64(&mut bytes, u64::from(OPTIMIZER_STREAM_KERNEL));
    push_u64(&mut bytes, u64::from(OPTIMIZER_STREAM_TILE));
    push_u64(&mut bytes, u64::from(fs_rand::STREAM_SEMANTICS_VERSION));
    push_str(&mut bytes, fs_rand::STREAM_POSITION_IDENTITY_DOMAIN);
    push_str(&mut bytes, fs_dfo::VERSION);
    push_str(&mut bytes, fs_math::VERSION);
    push_str(&mut bytes, fs_rand::VERSION);
    push_str(
        &mut bytes,
        "no-convergence-direction-coverage-hypervolume-diversity-superiority-all-seed-all-config-cross-process-cross-ISA-Cx-persistence-auth-internal-history-performance-claim",
    );
    bytes
}

fn config_bytes() -> Vec<u8> {
    let normalization = fs_dfo::moo::NSGA3_NORMALIZATION_POLICY.replay_identity();
    let mating = fs_dfo::moo::NSGA3_MATING_POLICY.replay_identity();
    config_bytes_with_policies(&normalization, &mating)
}

fn tri_objective(decision: &[f64]) -> Vec<f64> {
    let first = decision[0];
    let second = 1.0 - decision[0];
    let radial = decision[1..]
        .iter()
        .map(|value| {
            let centered = value - 0.5;
            centered * centered
        })
        .sum::<f64>();
    vec![first, second, radial]
}

fn evaluation_bits(decision: &[f64], objectives: &[f64]) -> EvaluationBits {
    EvaluationBits {
        decision: decision.iter().map(|value| value.to_bits()).collect(),
        objectives: objectives.iter().map(|value| value.to_bits()).collect(),
    }
}

fn run_study(config_digest: u64) -> SealedStudy {
    let mut evaluations = Vec::with_capacity(EXPECTED_EVALUATIONS);
    let front = {
        let mut objective = |decision: &[f64]| {
            let objectives = tri_objective(decision);
            evaluations.push(evaluation_bits(decision, &objectives));
            objectives
        };
        nsga3(
            &mut objective,
            DIMENSION,
            (LOWER_BOUND, UPPER_BOUND),
            &directions(),
            &parameters(),
        )
    };
    let front = front
        .iter()
        .map(|individual| evaluation_bits(&individual.x, &individual.f))
        .collect();
    SealedStudy::seal(config_digest, StudyRecord { evaluations, front })
}

fn decode(bits: &[u64]) -> Vec<f64> {
    bits.iter().copied().map(f64::from_bits).collect()
}

fn independent_dominates(left: &[u64], right: &[u64]) -> bool {
    if left.len() != right.len() || left.is_empty() {
        return false;
    }
    let mut strictly_better = false;
    for (&left, &right) in left.iter().zip(right) {
        let left = f64::from_bits(left);
        let right = f64::from_bits(right);
        if left > right {
            return false;
        }
        strictly_better |= left < right;
    }
    strictly_better
}

fn first_record_mismatch(expected: &StudyRecord, found: &StudyRecord) -> Option<String> {
    if expected.evaluations.len() != found.evaluations.len() {
        return Some(format!(
            "evaluations.length:{}!={}",
            found.evaluations.len(),
            expected.evaluations.len()
        ));
    }
    for (index, (expected, found)) in expected
        .evaluations
        .iter()
        .zip(&found.evaluations)
        .enumerate()
    {
        if expected != found {
            return Some(format!("evaluations[{index}]"));
        }
    }
    if expected.front.len() != found.front.len() {
        return Some(format!(
            "front.length:{}!={}",
            found.front.len(),
            expected.front.len()
        ));
    }
    for (individual, (expected, found)) in expected.front.iter().zip(&found.front).enumerate() {
        if expected.decision != found.decision {
            return Some(format!("front[{individual}].x"));
        }
        if expected.objectives.len() != found.objectives.len() {
            return Some(format!("front[{individual}].f.length"));
        }
        if let Some((objective, (&expected, &found))) = expected
            .objectives
            .iter()
            .zip(&found.objectives)
            .enumerate()
            .find(|(_, (expected, found))| expected != found)
        {
            return Some(format!(
                "front[{individual}].f[{objective}]:0x{found:016x}!=0x{expected:016x}"
            ));
        }
    }
    None
}

fn mutate_returned_objective(reference: &SealedStudy) -> CorruptionRun {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let individual = usize::try_from(selector.next_below(usize_u64(reference.record.front.len())))
        .expect("front index fits usize");
    let objective = usize::try_from(selector.next_below(usize_u64(OBJECTIVES)))
        .expect("objective index fits usize");
    let mantissa_bit = u32::try_from(selector.next_below(20)).expect("mantissa bit fits u32");
    let selector_draws = selector.index();

    let mut stale = reference.clone();
    let before = stale.record.front[individual].objectives[objective];
    let after = before ^ (1u64 << mantissa_bit);
    stale.record.front[individual].objectives[objective] = after;
    let stale_error = stale
        .validate_payload()
        .expect_err("unsealed returned-objective mutation must refuse");
    let mutant = SealedStudy::seal(reference.config_digest, stale.record);
    let reference_error = mutant
        .admit_against(reference.output_digest)
        .expect_err("resealed mutation must miss retained reference");
    let first_mismatch = first_record_mismatch(&reference.record, &mutant.record)
        .expect("one returned objective bit changes the record");
    let semantic_mismatch = mutant
        .record
        .semantic_mismatch()
        .expect("mutated returned objective must fail recomputation");
    let mutation = Mutation {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
        individual,
        objective,
        mantissa_bit,
        selector_draws,
        before,
        after,
    };
    let red_line = format!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"{RED_CASE}\",\"pass\":false,\"config\":\"0x{:016x}\",\"reference\":\"0x{:016x}\",\"mutant\":\"0x{:016x}\",\"corruption_seed\":\"0x{:016x}\",\"kernel\":\"0x{:04x}\",\"tile\":{},\"selector_draws\":{},\"target\":\"front[{}].f[{}]\",\"mantissa_bit\":{},\"before\":\"0x{:016x}\",\"after\":\"0x{:016x}\",\"stale_gate\":\"PayloadIdentityMismatch\",\"reference_gate\":\"ReferenceIdentityMismatch\",\"first_mismatch\":\"{}\",\"semantic_mismatch\":\"{}\"}}",
        reference.config_digest,
        reference.output_digest,
        mutant.output_digest,
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.individual,
        mutation.objective,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        first_mismatch,
        semantic_mismatch,
    );
    CorruptionRun {
        mutant,
        mutation,
        stale_error,
        reference_error,
        first_mismatch,
        semantic_mismatch,
        red_line,
    }
}

fn exact_returned_bit_delta(
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
    if expected.front[mutation.individual].objectives[mutation.objective] != mutation.before {
        return false;
    }
    expected.front[mutation.individual].objectives[mutation.objective] = mutation.after;
    expected == mutant.record
}

fn assert_mergeable(candidate: &SealedStudy, reference_output_digest: u64) {
    let payload_gate = candidate.validate_payload();
    let reference_gate = candidate.admit_against(reference_output_digest);
    let semantic_gate = candidate.record.semantic_mismatch();
    assert!(
        payload_gate.is_ok() && reference_gate.is_ok() && semantic_gate.is_none(),
        "merge gate refused {RED_CASE}: payload={payload_gate:?}; reference={reference_gate:?}; semantic={semantic_gate:?}"
    );
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| payload.downcast_ref::<&str>().map(ToString::to_string))
        .unwrap_or_else(|| "non-string panic payload".to_string())
}

#[test]
fn nsga3_study_config_consumes_shared_normalization_policy_root() {
    let policy = fs_dfo::moo::NSGA3_NORMALIZATION_POLICY;
    let current = policy.replay_identity();
    let mut mutant_policy = policy;
    mutant_policy.condition_error_limit *= 2.0;
    let mutant = mutant_policy.replay_identity();
    assert_ne!(current.root(), mutant.root());

    let mating = fs_dfo::moo::NSGA3_MATING_POLICY.replay_identity();
    let current_config = config_bytes_with_policies(&current, &mating);
    let mutant_config = config_bytes_with_policies(&mutant, &mating);
    assert_ne!(current_config, mutant_config);
    assert_ne!(
        fnv1a64(&current_config),
        fnv1a64(&mutant_config),
        "the retained study configuration must consume the shared typed policy root"
    );
}

#[test]
fn nsga3_study_config_consumes_shared_mating_policy_root() {
    let policy = fs_dfo::moo::NSGA3_MATING_POLICY;
    let current = policy.replay_identity();
    let mut mutant_policy = policy;
    mutant_policy.equal_rank_policy = "lower-live-index-mutant";
    let mutant = mutant_policy.replay_identity();
    assert_ne!(current.root(), mutant.root());

    let normalization = fs_dfo::moo::NSGA3_NORMALIZATION_POLICY.replay_identity();
    let current_config = config_bytes_with_policies(&normalization, &current);
    let mutant_config = config_bytes_with_policies(&normalization, &mutant);
    assert_ne!(current_config, mutant_config);
    assert_ne!(
        fnv1a64(&current_config),
        fnv1a64(&mutant_config),
        "the retained study configuration must consume the mating-policy root"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn nsga3_full_study_replays_and_seeded_failure_is_refused() {
    let config_frame = config_bytes();
    let config_digest = fnv1a64(&config_frame);
    let original = run_study(config_digest);
    let replay = run_study(config_digest);

    assert_eq!(original.record.semantic_mismatch(), None);
    assert_eq!(replay.record.semantic_mismatch(), None);
    assert_eq!(original.validate_payload(), Ok(()));
    assert_eq!(replay.validate_payload(), Ok(()));
    assert_eq!(original.admit_against(original.output_digest), Ok(()));
    assert_eq!(replay.admit_against(original.output_digest), Ok(()));
    assert_eq!(
        first_record_mismatch(&original.record, &replay.record),
        None,
        "complete callback and returned-front frames must repeat bit-for-bit"
    );
    assert_eq!(original, replay);
    assert_mergeable(&original, original.output_digest);

    let first = mutate_returned_objective(&original);
    let second = mutate_returned_objective(&replay);
    assert_eq!(first, second, "seeded red state must repeat exactly");
    assert!(exact_returned_bit_delta(
        &original,
        &first.mutant,
        first.mutation
    ));
    assert!(f64::from_bits(first.mutation.after).is_finite());
    assert!(matches!(
        first.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if declared == original.output_digest && computed == first.mutant.output_digest
    ));
    assert!(matches!(
        first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.output_digest && found == first.mutant.output_digest
    ));
    let expected_path = format!(
        "front[{}].f[{}]",
        first.mutation.individual, first.mutation.objective
    );
    assert!(first.first_mismatch.starts_with(&expected_path));
    assert_eq!(
        first.semantic_mismatch,
        format!("front[{}].f", first.mutation.individual)
    );
    assert_eq!(first.red_line, second.red_line);
    assert!(first.red_line.contains("PayloadIdentityMismatch"));
    assert!(first.red_line.contains("ReferenceIdentityMismatch"));
    assert!(first.red_line.contains(&expected_path));
    println!("{red_line}", red_line = first.red_line);

    let merge_gate = catch_unwind(|| assert_mergeable(&first.mutant, original.output_digest));
    let message = merge_gate
        .as_ref()
        .err()
        .map(|payload| panic_message(payload.as_ref()))
        .unwrap_or_default();
    assert!(message.contains(RED_CASE));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains(&first.semantic_mismatch));

    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"{CASE}\",\"pass\":true,\"config\":\"0x{config_digest:016x}\",\"output\":\"0x{:016x}\",\"callbacks\":{},\"front\":{},\"directions\":{},\"seed\":\"0x{STUDY_SEED:016x}\",\"stream_semantics_version\":{},\"scope\":\"same-process finite fixture\"}}",
        original.output_digest,
        original.record.evaluations.len(),
        original.record.front.len(),
        EXPECTED_DIRECTIONS,
        fs_rand::STREAM_SEMANTICS_VERSION,
    );
}
