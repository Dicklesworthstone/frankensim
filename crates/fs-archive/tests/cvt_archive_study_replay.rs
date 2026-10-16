//! G0/G3 public-state full-study replay for the production CVT archive.
//!
//! This target retains every candidate, nearest-centroid decision, acceptance
//! decision, and complete publicly observable post-operation state for one
//! finite deterministic CVT illumination study. A separate vector oracle
//! regenerates the candidates and reconstructs globally scaled nearest-centroid
//! assignment, lowest-index tie handling, strict replacement, coverage,
//! QD score, and unique-best semantics. A disclosed mutation changes one
//! retained per-step best-elite fitness bit and must fail stale-payload,
//! retained-reference, independent semantic, stable-red-record, and
//! evidence-derived merge gates.
//!
//! `CvtArchive` does not expose its per-centroid incumbents, so this receipt is
//! deliberately complete only over its public state. It makes no private-state
//! equivalence, centroid-generation/Lloyd, search-optimality,
//! coverage/quality-superiority, MAP-Elites, novelty, all-seed,
//! all-configuration, extreme-scale, cross-process, cross-ISA, cancellation,
//! persistence, authenticated-admission, or performance claim.

#![deny(unsafe_code)]

use fs_archive::{CvtArchive, Elite};
use std::panic::catch_unwind;

const SUITE: &str = "fs-archive/cvt-archive-study-replay-v1";
const CASE: &str = "seeded-cvt-public-state-replay";
const RED_CASE: &str = "seeded-step-best-fitness-corruption";
const STUDY_SEED: u64 = 0xC7A5_7A2E_0000_0043;
const MUTATION_SEED: u64 = 0xC7A5_7A2E_0000_00D5;
const CANDIDATES: usize = 25;
const CONTROL_NEW: usize = 0;
const CONTROL_EQUAL: usize = 1;
const CONTROL_WORSE: usize = 2;
const CONTROL_BETTER: usize = 3;
const GENERATED_START: usize = 4;
const DESCRIPTOR_DIMENSION: usize = 2;
const SOLUTION_DIMENSION: usize = 3;
const CENTROIDS: [[f64; DESCRIPTOR_DIMENSION]; 5] = [
    [-1.0, -1.0],
    [1.0, -1.0],
    [-1.0, 1.0],
    [1.0, 1.0],
    [0.0, 0.0],
];
const CAPACITY: usize = CENTROIDS.len();
const CONTROL_DESCRIPTOR: [f64; DESCRIPTOR_DIMENSION] = [0.0, -1.25];
const GENERATOR_FRAME_DIGEST: u64 = 0x20d1_ecfc_a8bd_bc95;
const GENERATOR_GOLDENS: [(usize, u64, u64, u64); 5] = [
    (4, 0, 0xc1cb_d4bf_fe1e_3141, 0x3fe8_397a_97ff_c3c6),
    (4, 1, 0xc5b4_f67a_5e81_9007, 0x3fe8_b69e_cf4b_d032),
    (8, 0, 0xcb45_98b7_4951_9b21, 0x3fe9_68b3_16e9_2a33),
    (13, 1, 0xbc01_9f6b_1f6f_298a, 0x3fe7_8033_ed63_ede5),
    (24, 0, 0x5ebd_4680_9ed4_0278, 0x3fd7_af51_a027_b500),
];

const _: () = assert!(CANDIDATES > CAPACITY && GENERATED_START == CONTROL_BETTER + 1);

#[derive(Debug, Clone, PartialEq, Eq)]
struct CandidateBits {
    solution: Vec<u64>,
    descriptor: Vec<u64>,
    fitness: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PublicStateBits {
    capacity: usize,
    num_elites: usize,
    coverage: u64,
    qd_score: u64,
    best: Option<CandidateBits>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StepBits {
    candidate: CandidateBits,
    nearest_centroid: usize,
    accepted: bool,
    state: PublicStateBits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    steps: Vec<StepBits>,
    final_state: PublicStateBits,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-archive-cvt-study-output-v1".to_vec();
        push_u64(&mut bytes, config_digest);
        push_len(&mut bytes, self.steps.len());
        for step in &self.steps {
            push_candidate(&mut bytes, &step.candidate);
            push_len(&mut bytes, step.nearest_centroid);
            bytes.push(u8::from(step.accepted));
            push_public_state(&mut bytes, &step.state);
        }
        push_public_state(&mut bytes, &self.final_state);
        bytes
    }

    #[allow(clippy::too_many_lines)]
    fn semantic_mismatch(&self) -> Option<String> {
        let generator_digest = generator_frame_digest();
        if generator_digest != GENERATOR_FRAME_DIGEST {
            return Some(format!(
                "generator.frame:0x{generator_digest:016x}!=expected-0x{GENERATOR_FRAME_DIGEST:016x}"
            ));
        }
        for &(candidate, lane, expected_word, expected_unit_bits) in &GENERATOR_GOLDENS {
            let word = generator_word(STUDY_SEED, candidate, lane);
            if word != expected_word {
                return Some(format!(
                    "generator.word[{candidate},{lane}]:0x{word:016x}!=expected-0x{expected_word:016x}"
                ));
            }
            let unit_bits = unit_interval(word).to_bits();
            if unit_bits != expected_unit_bits {
                return Some(format!(
                    "generator.unit[{candidate},{lane}]:0x{unit_bits:016x}!=expected-0x{expected_unit_bits:016x}"
                ));
            }
        }
        if self.steps.len() != CANDIDATES {
            return Some(format!(
                "steps.length:{}!=expected-{CANDIDATES}",
                self.steps.len()
            ));
        }

        // Deliberately use a fixed index-addressed vector instead of the
        // production BTreeMap. This independently reconstructs centroid
        // ownership and the exact centroid-index order used by QD-score.
        let mut oracle: Vec<Option<CandidateBits>> = vec![None; CAPACITY];
        let mut previous_coverage = 0.0;
        let mut previous_qd_score = 0.0;
        for (index, step) in self.steps.iter().enumerate() {
            let generated = oracle_candidate(index);
            let expected_candidate = candidate_bits(&generated);
            if step.candidate != expected_candidate {
                return Some(format!("steps[{index}].candidate"));
            }

            let nearest = independent_nearest_centroid(&generated.descriptor);
            if step.nearest_centroid != nearest {
                return Some(format!(
                    "steps[{index}].nearest_centroid:{}!=expected-{nearest}",
                    step.nearest_centroid
                ));
            }
            let accepted = oracle[nearest]
                .as_ref()
                .is_none_or(|incumbent| generated.fitness > f64::from_bits(incumbent.fitness));
            if accepted {
                oracle[nearest] = Some(expected_candidate.clone());
            }
            if step.accepted != accepted {
                return Some(format!(
                    "steps[{index}].accepted:{}!=expected-{accepted}",
                    step.accepted
                ));
            }

            let expected_state = oracle_public_state(&oracle);
            if let Some(field) = public_state_mismatch(&step.state, &expected_state) {
                return Some(format!("steps[{index}].{field}"));
            }
            let expected_best_fitness = expected_state
                .best
                .as_ref()
                .map(|best| best.fitness)
                .expect("a post-add oracle state has a best elite");
            let best_multiplicity = oracle
                .iter()
                .filter_map(Option::as_ref)
                .filter(|elite| elite.fitness == expected_best_fitness)
                .count();
            if best_multiplicity != 1 {
                return Some(format!(
                    "steps[{index}].best_multiplicity:{best_multiplicity}!=expected-1"
                ));
            }
            let coverage = f64::from_bits(step.state.coverage);
            let qd_score = f64::from_bits(step.state.qd_score);
            if coverage < previous_coverage {
                return Some(format!("steps[{index}].coverage.non_monotone"));
            }
            if qd_score < previous_qd_score {
                return Some(format!("steps[{index}].qd_score.non_monotone"));
            }
            previous_coverage = coverage;
            previous_qd_score = qd_score;

            if index >= GENERATED_START {
                let intended = (index - GENERATED_START) % CAPACITY;
                if nearest != intended || !accepted {
                    return Some(format!("steps[{index}].generated-rotation-control"));
                }
            }
        }

        let first = &self.steps[CONTROL_NEW];
        let equal = &self.steps[CONTROL_EQUAL];
        let worse = &self.steps[CONTROL_WORSE];
        let better = &self.steps[CONTROL_BETTER];
        if [
            first.accepted,
            equal.accepted,
            worse.accepted,
            better.accepted,
        ] != [true, false, false, true]
        {
            return Some("controls.acceptance-vector".to_string());
        }
        let tie_distances = independent_distance_frame(&CONTROL_DESCRIPTOR);
        if tie_distances[0].to_bits() != tie_distances[1].to_bits()
            || tie_distances
                .iter()
                .skip(2)
                .any(|distance| *distance <= tie_distances[0])
        {
            return Some("controls.exact-two-way-nearest-tie".to_string());
        }
        if first.nearest_centroid != 0 || !first.accepted {
            return Some("controls.new-niche-and-lowest-index-tie".to_string());
        }
        if equal.candidate.descriptor != first.candidate.descriptor
            || equal.candidate.fitness != first.candidate.fitness
            || equal.candidate.solution == first.candidate.solution
            || equal.accepted
            || equal.state != first.state
        {
            return Some("controls.strict-equality-refusal".to_string());
        }
        if worse.candidate.descriptor != first.candidate.descriptor
            || f64::from_bits(worse.candidate.fitness) >= f64::from_bits(first.candidate.fitness)
            || worse.accepted
            || worse.state != first.state
        {
            return Some("controls.worse-refusal".to_string());
        }
        if better.candidate.descriptor != first.candidate.descriptor
            || f64::from_bits(better.candidate.fitness) <= f64::from_bits(first.candidate.fitness)
            || !better.accepted
            || better.state.num_elites != 1
            || better.state.best.as_ref() != Some(&better.candidate)
        {
            return Some("controls.strict-improvement-replacement".to_string());
        }

        if oracle.iter().any(Option::is_none) {
            return Some("fixture.not-all-centroids-filled".to_string());
        }
        let expected_final = oracle_public_state(&oracle);
        if let Some(field) = public_state_mismatch(&self.final_state, &expected_final) {
            return Some(format!("final_state.{field}"));
        }
        if self.final_state.num_elites != CAPACITY || self.final_state.coverage != 1.0_f64.to_bits()
        {
            return Some("final_state.full-public-coverage".to_string());
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

#[derive(Debug, Clone)]
struct Candidate {
    solution: Vec<f64>,
    descriptor: Vec<f64>,
    fitness: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    step: usize,
    mantissa_bit: u32,
    selector_word: u64,
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

fn push_candidate(bytes: &mut Vec<u8>, candidate: &CandidateBits) {
    push_u64_slice(bytes, &candidate.solution);
    push_u64_slice(bytes, &candidate.descriptor);
    push_u64(bytes, candidate.fitness);
}

fn push_public_state(bytes: &mut Vec<u8>, state: &PublicStateBits) {
    push_len(bytes, state.capacity);
    push_len(bytes, state.num_elites);
    push_u64(bytes, state.coverage);
    push_u64(bytes, state.qd_score);
    match &state.best {
        Some(best) => {
            bytes.push(1);
            push_candidate(bytes, best);
        }
        None => bytes.push(0),
    }
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn generator_word(seed: u64, candidate: usize, lane: u64) -> u64 {
    let ordinal = usize_u64(candidate).wrapping_mul(0xd1b5_4a32_d192_ed03);
    splitmix64(seed ^ ordinal ^ lane.wrapping_mul(0xa24b_aed4_963e_e407))
}

fn generator_frame_digest() -> u64 {
    let mut bytes = b"fs-archive-cvt-generator-frame-v1".to_vec();
    push_len(&mut bytes, CANDIDATES - GENERATED_START);
    for candidate in GENERATED_START..CANDIDATES {
        push_len(&mut bytes, candidate);
        for lane in 0..2 {
            push_u64(&mut bytes, lane);
            push_u64(&mut bytes, generator_word(STUDY_SEED, candidate, lane));
        }
    }
    fnv1a64(&bytes)
}

fn unit_interval(word: u64) -> f64 {
    (word >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
}

fn dyadic_jitter(word: u64) -> f64 {
    let numerator = i32::try_from((word >> 60) & 0x0f).expect("nibble fits i32") * 2 - 15;
    f64::from(numerator) / 128.0
}

fn control_candidate(index: usize, fitness: f64) -> Candidate {
    Candidate {
        solution: vec![index as f64, 0.0, 0.0],
        descriptor: CONTROL_DESCRIPTOR.to_vec(),
        fitness,
    }
}

fn candidate(index: usize) -> Candidate {
    match index {
        CONTROL_NEW => control_candidate(index, 4.0),
        CONTROL_EQUAL => control_candidate(index, 4.0),
        CONTROL_WORSE => control_candidate(index, 3.0),
        CONTROL_BETTER => control_candidate(index, 6.0),
        GENERATED_START..CANDIDATES => {
            let centroid = CENTROIDS[(index - GENERATED_START) % CAPACITY];
            let x_word = generator_word(STUDY_SEED, index, 0);
            let y_word = generator_word(STUDY_SEED, index, 1);
            Candidate {
                solution: vec![
                    index as f64 / 32.0,
                    unit_interval(x_word),
                    unit_interval(y_word),
                ],
                descriptor: vec![
                    centroid[0] + dyadic_jitter(x_word),
                    centroid[1] + dyadic_jitter(y_word),
                ],
                fitness: 8.0 + index as f64 / 32.0,
            }
        }
        _ => unreachable!("study index is bounded by CANDIDATES"),
    }
}

fn oracle_candidate(index: usize) -> Candidate {
    if index <= CONTROL_BETTER {
        let fitness = match index {
            CONTROL_NEW | CONTROL_EQUAL => 4.0,
            CONTROL_WORSE => 3.0,
            CONTROL_BETTER => 6.0,
            _ => unreachable!("control index is bounded"),
        };
        return Candidate {
            solution: vec![index as f64, 0.0, 0.0],
            descriptor: vec![0.0, -5.0 / 4.0],
            fitness,
        };
    }

    let centroid_index = (index - 4) % 5;
    let base = CENTROIDS[centroid_index];
    let x_word = generator_word(STUDY_SEED, index, 0);
    let y_word = generator_word(STUDY_SEED, index, 1);
    let x_nibble = ((x_word >> 60) & 15) as f64;
    let y_nibble = ((y_word >> 60) & 15) as f64;
    let denominator = 9_007_199_254_740_992.0;
    Candidate {
        solution: vec![
            index as f64 / 32.0,
            (x_word >> 11) as f64 / denominator,
            (y_word >> 11) as f64 / denominator,
        ],
        descriptor: vec![
            base[0] + (x_nibble / 64.0 - 15.0 / 128.0),
            base[1] + (y_nibble / 64.0 - 15.0 / 128.0),
        ],
        fitness: 8.0 + index as f64 / 32.0,
    }
}

fn candidate_bits(candidate: &Candidate) -> CandidateBits {
    CandidateBits {
        solution: candidate
            .solution
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        descriptor: candidate
            .descriptor
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        fitness: candidate.fitness.to_bits(),
    }
}

fn elite_bits(elite: &Elite) -> CandidateBits {
    CandidateBits {
        solution: elite.solution.iter().map(|value| value.to_bits()).collect(),
        descriptor: elite
            .descriptor
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        fitness: elite.fitness.to_bits(),
    }
}

fn independent_distance_frame(descriptor: &[f64]) -> [f64; CAPACITY] {
    let global_scale = CENTROIDS
        .iter()
        .flat_map(|centroid| centroid.iter().zip(descriptor))
        .map(|(coordinate, query)| (coordinate - query).abs())
        .fold(0.0_f64, f64::max);
    if global_scale == 0.0 {
        return [0.0; CAPACITY];
    }
    std::array::from_fn(|index| {
        let centroid = &CENTROIDS[index];
        centroid
            .iter()
            .zip(descriptor)
            .map(|(coordinate, query)| {
                let delta = (coordinate - query) / global_scale;
                delta * delta
            })
            .sum::<f64>()
    })
}

fn independent_nearest_centroid(descriptor: &[f64]) -> usize {
    let distances = independent_distance_frame(descriptor);
    let mut nearest = 0;
    let mut minimum = distances[0];
    for (index, &distance) in distances.iter().enumerate().skip(1) {
        // Strict comparison deliberately retains the lowest index on a tie.
        if distance < minimum {
            nearest = index;
            minimum = distance;
        }
    }
    nearest
}

fn oracle_public_state(oracle: &[Option<CandidateBits>]) -> PublicStateBits {
    let num_elites = oracle.iter().filter(|elite| elite.is_some()).count();
    let qd_score = oracle
        .iter()
        .filter_map(Option::as_ref)
        .map(|elite| f64::from_bits(elite.fitness))
        .sum::<f64>();
    let mut best: Option<CandidateBits> = None;
    for elite in oracle.iter().filter_map(Option::as_ref) {
        if best.as_ref().is_none_or(|incumbent| {
            f64::from_bits(elite.fitness) > f64::from_bits(incumbent.fitness)
        }) {
            best = Some(elite.clone());
        }
    }
    PublicStateBits {
        capacity: CAPACITY,
        num_elites,
        coverage: (num_elites as f64 / CAPACITY as f64).to_bits(),
        qd_score: qd_score.to_bits(),
        best,
    }
}

fn archive_public_state(archive: &CvtArchive) -> PublicStateBits {
    PublicStateBits {
        capacity: archive.capacity(),
        num_elites: archive.num_elites(),
        coverage: archive.coverage().to_bits(),
        qd_score: archive.qd_score().to_bits(),
        best: archive.best().map(elite_bits),
    }
}

fn candidate_mismatch(found: &CandidateBits, expected: &CandidateBits) -> Option<&'static str> {
    if found.solution != expected.solution {
        Some("solution")
    } else if found.descriptor != expected.descriptor {
        Some("descriptor")
    } else if found.fitness != expected.fitness {
        Some("fitness")
    } else {
        None
    }
}

fn public_state_mismatch(found: &PublicStateBits, expected: &PublicStateBits) -> Option<String> {
    if found.capacity != expected.capacity {
        return Some("capacity".to_string());
    }
    if found.num_elites != expected.num_elites {
        return Some("num_elites".to_string());
    }
    if found.coverage != expected.coverage {
        return Some("coverage".to_string());
    }
    if found.qd_score != expected.qd_score {
        return Some("qd_score".to_string());
    }
    match (&found.best, &expected.best) {
        (Some(found), Some(expected)) => {
            candidate_mismatch(found, expected).map(|field| format!("best.{field}"))
        }
        (None, None) => None,
        _ => Some("best.presence".to_string()),
    }
}

fn config_bytes() -> Vec<u8> {
    let mut bytes = b"fs-archive-cvt-study-config-v1".to_vec();
    push_str(&mut bytes, CASE);
    push_str(&mut bytes, "fs_archive::CvtArchive");
    push_str(&mut bytes, "splitmix64-indexed-dyadic-jitter-v1");
    push_str(&mut bytes, "dimensionless");
    push_u64(&mut bytes, STUDY_SEED);
    push_len(&mut bytes, CANDIDATES);
    push_len(&mut bytes, SOLUTION_DIMENSION);
    push_len(&mut bytes, DESCRIPTOR_DIMENSION);
    push_len(&mut bytes, CENTROIDS.len());
    for centroid in CENTROIDS {
        for coordinate in centroid {
            push_u64(&mut bytes, coordinate.to_bits());
        }
    }
    for coordinate in CONTROL_DESCRIPTOR {
        push_u64(&mut bytes, coordinate.to_bits());
    }
    for control in [CONTROL_NEW, CONTROL_EQUAL, CONTROL_WORSE, CONTROL_BETTER] {
        push_len(&mut bytes, control);
    }
    push_u64(&mut bytes, GENERATOR_FRAME_DIGEST);
    push_len(&mut bytes, GENERATOR_GOLDENS.len());
    for (candidate, lane, word, unit_bits) in GENERATOR_GOLDENS {
        push_len(&mut bytes, candidate);
        push_u64(&mut bytes, lane);
        push_u64(&mut bytes, word);
        push_u64(&mut bytes, unit_bits);
    }
    push_str(&mut bytes, env!("CARGO_PKG_VERSION"));
    push_str(
        &mut bytes,
        "complete-public-state-only-no-private-state-centroid-generation-optimality-quality-superiority-MAP-Elites-novelty-all-seed-all-config-extreme-scale-cross-process-cross-ISA-Cx-persistence-auth-performance-claim",
    );
    bytes
}

fn run_study(config_digest: u64) -> SealedStudy {
    let mut archive = CvtArchive::new(CENTROIDS.into_iter().map(|row| row.to_vec()).collect());
    let mut steps = Vec::with_capacity(CANDIDATES);
    for index in 0..CANDIDATES {
        let candidate = candidate(index);
        let bits = candidate_bits(&candidate);
        let nearest_centroid = archive.nearest_centroid(&candidate.descriptor);
        let accepted = archive.add(candidate.solution, candidate.descriptor, candidate.fitness);
        steps.push(StepBits {
            candidate: bits,
            nearest_centroid,
            accepted,
            state: archive_public_state(&archive),
        });
    }
    let record = StudyRecord {
        steps,
        final_state: archive_public_state(&archive),
    };
    SealedStudy::seal(config_digest, record)
}

fn first_record_mismatch(expected: &StudyRecord, found: &StudyRecord) -> Option<String> {
    if expected.steps.len() != found.steps.len() {
        return Some(format!(
            "steps.length:{}!={}",
            found.steps.len(),
            expected.steps.len()
        ));
    }
    for (index, (expected, found)) in expected.steps.iter().zip(&found.steps).enumerate() {
        if let Some(field) = candidate_mismatch(&found.candidate, &expected.candidate) {
            return Some(format!("steps[{index}].candidate.{field}"));
        }
        if expected.nearest_centroid != found.nearest_centroid {
            return Some(format!("steps[{index}].nearest_centroid"));
        }
        if expected.accepted != found.accepted {
            return Some(format!("steps[{index}].accepted"));
        }
        if let Some(field) = public_state_mismatch(&found.state, &expected.state) {
            return Some(format!("steps[{index}].{field}"));
        }
    }
    public_state_mismatch(&found.final_state, &expected.final_state)
        .map(|field| format!("final_state.{field}"))
}

fn mutate_step_best(reference: &SealedStudy) -> CorruptionRun {
    let selector_word = generator_word(MUTATION_SEED, reference.record.steps.len(), 0);
    let step = usize::try_from(selector_word % usize_u64(reference.record.steps.len()))
        .expect("step index fits usize");
    let mantissa_bit = u32::try_from((selector_word >> 32) % 20).expect("mantissa bit fits u32");
    let mut stale = reference.clone();
    let retained_best = stale.record.steps[step]
        .state
        .best
        .as_mut()
        .expect("every post-add study state has a best elite");
    let before = retained_best.fitness;
    let after = before ^ (1u64 << mantissa_bit);
    retained_best.fitness = after;
    let stale_error = stale
        .validate_payload()
        .expect_err("unsealed retained-best mutation must refuse");
    let mutant = SealedStudy::seal(reference.config_digest, stale.record);
    let reference_error = mutant
        .admit_against(reference.output_digest)
        .expect_err("resealed retained-best mutation must miss retained reference");
    let first_mismatch = first_record_mismatch(&reference.record, &mutant.record)
        .expect("one retained best-fitness bit changes the record");
    let semantic_mismatch = mutant
        .record
        .semantic_mismatch()
        .expect("mutated retained best must fail independent reconstruction");
    let mutation = Mutation {
        seed: MUTATION_SEED,
        step,
        mantissa_bit,
        selector_word,
        before,
        after,
    };
    let red_line = format!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"{RED_CASE}\",\"pass\":false,\"config\":\"0x{:016x}\",\"reference\":\"0x{:016x}\",\"mutant\":\"0x{:016x}\",\"corruption_seed\":\"0x{:016x}\",\"selector_word\":\"0x{:016x}\",\"target\":\"steps[{}].best.fitness\",\"mantissa_bit\":{},\"before\":\"0x{:016x}\",\"after\":\"0x{:016x}\",\"stale_gate\":\"PayloadIdentityMismatch\",\"reference_gate\":\"ReferenceIdentityMismatch\",\"first_mismatch\":\"{}\",\"semantic_mismatch\":\"{}\"}}",
        reference.config_digest,
        reference.output_digest,
        mutant.output_digest,
        mutation.seed,
        mutation.selector_word,
        mutation.step,
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

fn exact_best_bit_delta(reference: &SealedStudy, mutant: &SealedStudy, mutation: Mutation) -> bool {
    if reference.config_digest != mutant.config_digest
        || mutation.before ^ mutation.after != 1u64 << mutation.mantissa_bit
    {
        return false;
    }
    let mut expected = reference.record.clone();
    let Some(best) = expected.record_step_best_mut(mutation.step) else {
        return false;
    };
    if best.fitness != mutation.before {
        return false;
    }
    best.fitness = mutation.after;
    expected == mutant.record
}

impl StudyRecord {
    fn record_step_best_mut(&mut self, step: usize) -> Option<&mut CandidateBits> {
        self.steps.get_mut(step)?.state.best.as_mut()
    }
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
#[allow(clippy::too_many_lines)]
fn cvt_archive_public_state_replays_and_seeded_failure_is_refused() {
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
        "complete candidate, assignment, decision, and public-state frames must repeat bit-for-bit"
    );
    assert_eq!(original, replay);
    assert_mergeable(&original, original.output_digest);

    let first = mutate_step_best(&original);
    let second = mutate_step_best(&replay);
    assert_eq!(first, second, "seeded red state must repeat exactly");
    assert!(exact_best_bit_delta(
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
    let expected_path = format!("steps[{}].best.fitness", first.mutation.step);
    assert_eq!(first.first_mismatch, expected_path);
    assert_eq!(first.semantic_mismatch, expected_path);
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
    assert!(message.contains(&expected_path));

    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"{CASE}\",\"pass\":true,\"config\":\"0x{config_digest:016x}\",\"output\":\"0x{:016x}\",\"candidates\":{},\"accepted\":{},\"capacity\":{},\"elites\":{},\"coverage_bits\":\"0x{:016x}\",\"qd_score_bits\":\"0x{:016x}\",\"seed\":\"0x{STUDY_SEED:016x}\",\"scope\":\"same-process finite complete-public-state fixture\"}}",
        original.output_digest,
        original.record.steps.len(),
        original
            .record
            .steps
            .iter()
            .filter(|step| step.accepted)
            .count(),
        original.record.final_state.capacity,
        original.record.final_state.num_elites,
        original.record.final_state.coverage,
        original.record.final_state.qd_score,
    );
}
