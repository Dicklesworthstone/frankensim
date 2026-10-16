//! G0/G3 full-study replay for the production MAP-Elites archive.
//!
//! This target retains every candidate, acceptance decision, post-operation
//! incumbent, complete ordered archive frame, aggregate metric transition,
//! and final elite for one finite deterministic illumination study. A
//! separate sorted-vector oracle regenerates candidates and
//! reconstructs cell assignment, strict replacement, coverage, QD score,
//! unique best-elite choice, and final iteration order. A disclosed mutation
//! changes one final-elite fitness bit and must fail stale-payload,
//! retained-reference, independent semantic, stable-red-record, and
//! evidence-derived merge gates.
//!
//! This target makes no search-optimality, coverage/quality-superiority, CVT,
//! novelty, all-seed, all-configuration, cross-process, cross-ISA,
//! cancellation, persistence, authenticated-admission, or performance claim.

#![deny(unsafe_code)]

use fs_archive::{Elite, MapElites};
use std::panic::catch_unwind;

const SUITE: &str = "fs-archive/map-elites-study-replay-v1";
const CASE: &str = "seeded-grid-illumination-complete-frame";
const RED_CASE: &str = "seeded-final-elite-fitness-corruption";
const STUDY_SEED: u64 = 0xA7C4_1E57_0000_0040;
const CANDIDATES: usize = 32;
const EQUALITY_SOURCE: usize = 0;
const EQUALITY_CANDIDATE: usize = 1;
const DESCRIPTOR_DIMENSION: usize = 2;
const SOLUTION_DIMENSION: usize = 3;
const LOWER: [f64; DESCRIPTOR_DIMENSION] = [0.0, 0.0];
const UPPER: [f64; DESCRIPTOR_DIMENSION] = [1.0, 1.0];
const BINS: [usize; DESCRIPTOR_DIMENSION] = [4, 4];
const CAPACITY: usize = BINS[0] * BINS[1];
const MUTATION_SEED: u64 = 0xA7C4_1E57_0000_00D5;
const GENERATOR_FRAME_DIGEST: u64 = 0x718a_6ed7_ccbd_d877;
const GENERATOR_GOLDENS: [(usize, u64, u64, u64); 5] = [
    (0, 0, 0xf812_a1a5_18a4_4204, 0x3fef_0254_34a3_1488),
    (0, 1, 0xa24d_9e9e_c9c9_a924, 0x3fe4_49b3_d3d9_3935),
    (2, 0, 0x1d54_2041_32f3_3788, 0x3fbd_5420_4132_f330),
    (7, 1, 0x6a67_40d9_48ab_5619, 0x3fda_99d0_3652_2ad4),
    (31, 0, 0x8cff_479f_72d0_4ce3, 0x3fe1_9fe8_f3ee_5a09),
];

const _: () = assert!(CAPACITY == 16 && CANDIDATES > CAPACITY);
const _: () = assert!(EQUALITY_CANDIDATE == EQUALITY_SOURCE + 1);

#[derive(Debug, Clone, PartialEq, Eq)]
struct CandidateBits {
    solution: Vec<u64>,
    descriptor: Vec<u64>,
    fitness: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StepBits {
    candidate: CandidateBits,
    cell: Vec<usize>,
    accepted: bool,
    incumbent: CandidateBits,
    archive_elites: Vec<CandidateBits>,
    coverage: u64,
    qd_score: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    steps: Vec<StepBits>,
    capacity: usize,
    num_elites: usize,
    coverage: u64,
    qd_score: u64,
    best: Option<CandidateBits>,
    final_elites: Vec<CandidateBits>,
}

impl StudyRecord {
    fn canonical_bytes(&self, config_digest: u64) -> Vec<u8> {
        let mut bytes = b"fs-archive-map-elites-study-output-v1".to_vec();
        push_u64(&mut bytes, config_digest);
        push_len(&mut bytes, self.steps.len());
        for step in &self.steps {
            push_candidate(&mut bytes, &step.candidate);
            push_usize_slice(&mut bytes, &step.cell);
            bytes.push(u8::from(step.accepted));
            push_candidate(&mut bytes, &step.incumbent);
            push_len(&mut bytes, step.archive_elites.len());
            for elite in &step.archive_elites {
                push_candidate(&mut bytes, elite);
            }
            push_u64(&mut bytes, step.coverage);
            push_u64(&mut bytes, step.qd_score);
        }
        push_len(&mut bytes, self.capacity);
        push_len(&mut bytes, self.num_elites);
        push_u64(&mut bytes, self.coverage);
        push_u64(&mut bytes, self.qd_score);
        match &self.best {
            Some(best) => {
                bytes.push(1);
                push_candidate(&mut bytes, best);
            }
            None => bytes.push(0),
        }
        push_len(&mut bytes, self.final_elites.len());
        for elite in &self.final_elites {
            push_candidate(&mut bytes, elite);
        }
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
        // Deliberately use a sorted vector instead of production's BTreeMap:
        // this separately reconstructs both cell ownership and iteration order.
        let mut oracle: Vec<(Vec<usize>, CandidateBits)> = Vec::new();
        for (index, step) in self.steps.iter().enumerate() {
            let generated = oracle_candidate(index);
            let expected_candidate = candidate_bits(&generated);
            if step.candidate != expected_candidate {
                return Some(format!("steps[{index}].candidate"));
            }
            let expected_cell = independent_cell(&generated.descriptor);
            if step.cell != expected_cell {
                return Some(format!(
                    "steps[{index}].cell:{:?}!=expected-{expected_cell:?}",
                    step.cell
                ));
            }
            let position = oracle.iter().position(|(cell, _)| cell == &expected_cell);
            let accepted = position.is_none_or(|position| {
                generated.fitness > f64::from_bits(oracle[position].1.fitness)
            });
            if accepted {
                if let Some(position) = position {
                    oracle[position].1 = expected_candidate;
                } else {
                    oracle.push((expected_cell, expected_candidate));
                    oracle.sort_by(|left, right| left.0.cmp(&right.0));
                }
            }
            if step.accepted != accepted {
                return Some(format!(
                    "steps[{index}].accepted:{}!=expected-{accepted}",
                    step.accepted
                ));
            }
            let incumbent = oracle
                .iter()
                .find(|(cell, _)| cell == &step.cell)
                .map(|(_, elite)| elite)
                .expect("the candidate's cell has an incumbent after insertion");
            if &step.incumbent != incumbent {
                return Some(format!("steps[{index}].incumbent"));
            }
            if index == EQUALITY_CANDIDATE {
                let source = &self.steps[EQUALITY_SOURCE];
                if step.candidate.descriptor != source.candidate.descriptor
                    || step.candidate.fitness != source.candidate.fitness
                    || step.candidate.solution == source.candidate.solution
                    || step.accepted
                    || step.incumbent != source.incumbent
                {
                    return Some("strict-equality-replacement-control".to_string());
                }
            }
            let expected_archive: Vec<CandidateBits> =
                oracle.iter().map(|(_, elite)| elite.clone()).collect();
            if step.archive_elites.len() != expected_archive.len() {
                return Some(format!(
                    "steps[{index}].archive_elites.length:{}!=expected-{}",
                    step.archive_elites.len(),
                    expected_archive.len()
                ));
            }
            for (elite, (found, expected)) in step
                .archive_elites
                .iter()
                .zip(&expected_archive)
                .enumerate()
            {
                if found != expected {
                    return Some(format!("steps[{index}].archive_elites[{elite}]"));
                }
            }
            let coverage = oracle.len() as f64 / CAPACITY as f64;
            let qd_score = oracle
                .iter()
                .map(|(_, elite)| f64::from_bits(elite.fitness))
                .sum::<f64>();
            if step.coverage != coverage.to_bits() {
                return Some(format!(
                    "steps[{index}].coverage:0x{:016x}!=expected-0x{:016x}",
                    step.coverage,
                    coverage.to_bits()
                ));
            }
            if step.qd_score != qd_score.to_bits() {
                return Some(format!(
                    "steps[{index}].qd_score:0x{:016x}!=expected-0x{:016x}",
                    step.qd_score,
                    qd_score.to_bits()
                ));
            }
        }

        let expected_elites: Vec<CandidateBits> =
            oracle.iter().map(|(_, elite)| elite.clone()).collect();
        if self.final_elites.len() != expected_elites.len() {
            return Some(format!(
                "final_elites.length:{}!=expected-{}",
                self.final_elites.len(),
                expected_elites.len()
            ));
        }
        for (index, (found, expected)) in self.final_elites.iter().zip(&expected_elites).enumerate()
        {
            if found != expected {
                let field = if found.solution != expected.solution {
                    "solution"
                } else if found.descriptor != expected.descriptor {
                    "descriptor"
                } else {
                    "fitness"
                };
                return Some(format!("final_elites[{index}].{field}"));
            }
        }
        let mut expected_best: Option<CandidateBits> = None;
        let mut best_multiplicity = 0usize;
        for (_, elite) in &oracle {
            let fitness = f64::from_bits(elite.fitness);
            let Some(incumbent) = &expected_best else {
                expected_best = Some(elite.clone());
                best_multiplicity = 1;
                continue;
            };
            let incumbent_bits = incumbent.fitness;
            let incumbent_fitness = f64::from_bits(incumbent_bits);
            if fitness > incumbent_fitness {
                expected_best = Some(elite.clone());
                best_multiplicity = 1;
            } else if elite.fitness == incumbent_bits {
                best_multiplicity += 1;
            }
        }
        if best_multiplicity != 1 {
            return Some(format!(
                "fixture.best_multiplicity:{best_multiplicity}!=expected-1"
            ));
        }
        if self.best != expected_best {
            return Some("best".to_string());
        }
        let coverage = oracle.len() as f64 / CAPACITY as f64;
        let qd_score = oracle
            .iter()
            .map(|(_, elite)| f64::from_bits(elite.fitness))
            .sum::<f64>();
        if self.capacity != CAPACITY {
            return Some(format!("capacity:{}!=expected-{CAPACITY}", self.capacity));
        }
        if self.num_elites != oracle.len() {
            return Some(format!(
                "num_elites:{}!=expected-{}",
                self.num_elites,
                oracle.len()
            ));
        }
        if self.coverage != coverage.to_bits() {
            return Some("coverage".to_string());
        }
        if self.qd_score != qd_score.to_bits() {
            return Some("qd_score".to_string());
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
    elite: usize,
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

fn push_usize_slice(bytes: &mut Vec<u8>, values: &[usize]) {
    push_len(bytes, values.len());
    for &value in values {
        push_len(bytes, value);
    }
}

fn push_candidate(bytes: &mut Vec<u8>, candidate: &CandidateBits) {
    push_u64_slice(bytes, &candidate.solution);
    push_u64_slice(bytes, &candidate.descriptor);
    push_u64(bytes, candidate.fitness);
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
    let mut bytes = b"fs-archive-map-elites-generator-frame-v1".to_vec();
    push_len(&mut bytes, CANDIDATES);
    for candidate in 0..CANDIDATES {
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

fn raw_candidate(index: usize) -> Candidate {
    let x = unit_interval(generator_word(STUDY_SEED, index, 0));
    let y = unit_interval(generator_word(STUDY_SEED, index, 1));
    let ordinal = index as f64 / CANDIDATES as f64;
    let fitness = (0.5 * x).mul_add(1.0, (0.75 * y).mul_add(1.0, ordinal + 1.0));
    Candidate {
        solution: vec![x, y, ordinal],
        descriptor: vec![x, y],
        fitness,
    }
}

fn candidate(index: usize) -> Candidate {
    if index == EQUALITY_CANDIDATE {
        let mut equal = raw_candidate(EQUALITY_SOURCE);
        equal.solution[2] = index as f64 / CANDIDATES as f64;
        equal
    } else {
        raw_candidate(index)
    }
}

fn oracle_candidate(index: usize) -> Candidate {
    let generator_index = if index == EQUALITY_CANDIDATE {
        EQUALITY_SOURCE
    } else {
        index
    };
    let x_word = generator_word(STUDY_SEED, generator_index, 0);
    let y_word = generator_word(STUDY_SEED, generator_index, 1);
    let x = (x_word >> 11) as f64 / 9_007_199_254_740_992.0;
    let y = (y_word >> 11) as f64 / 9_007_199_254_740_992.0;
    let solution_ordinal = index as f64 / CANDIDATES as f64;
    let fitness_ordinal = generator_index as f64 / CANDIDATES as f64;
    let fitness = (0.5 * x).mul_add(1.0, (0.75 * y).mul_add(1.0, fitness_ordinal + 1.0));
    Candidate {
        solution: vec![x, y, solution_ordinal],
        descriptor: vec![x, y],
        fitness,
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

fn independent_cell(descriptor: &[f64]) -> Vec<usize> {
    descriptor
        .iter()
        .enumerate()
        .map(|(axis, value)| {
            let fraction = (value - LOWER[axis]) / (UPPER[axis] - LOWER[axis]);
            let raw = (fraction * BINS[axis] as f64).floor();
            if raw < 0.0 {
                0
            } else if raw >= BINS[axis] as f64 {
                BINS[axis] - 1
            } else {
                raw as usize
            }
        })
        .collect()
}

fn config_bytes() -> Vec<u8> {
    let mut bytes = b"fs-archive-map-elites-study-config-v1".to_vec();
    push_str(&mut bytes, CASE);
    push_str(&mut bytes, "fs_archive::MapElites");
    push_str(&mut bytes, "splitmix64-indexed-candidate-generator-v1");
    push_str(&mut bytes, "dimensionless");
    push_u64(&mut bytes, STUDY_SEED);
    push_len(&mut bytes, CANDIDATES);
    push_len(&mut bytes, EQUALITY_SOURCE);
    push_len(&mut bytes, EQUALITY_CANDIDATE);
    push_len(&mut bytes, SOLUTION_DIMENSION);
    push_len(&mut bytes, DESCRIPTOR_DIMENSION);
    for value in LOWER.into_iter().chain(UPPER) {
        push_u64(&mut bytes, value.to_bits());
    }
    push_usize_slice(&mut bytes, &BINS);
    push_len(&mut bytes, CAPACITY);
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
        "no-optimality-coverage-quality-superiority-CVT-novelty-all-seed-all-config-cross-process-cross-ISA-Cx-persistence-auth-performance-claim",
    );
    bytes
}

fn run_study(config_digest: u64) -> SealedStudy {
    let mut archive = MapElites::new(LOWER.to_vec(), UPPER.to_vec(), BINS.to_vec());
    let mut steps = Vec::with_capacity(CANDIDATES);
    for index in 0..CANDIDATES {
        let candidate = candidate(index);
        let descriptor_query = candidate.descriptor.clone();
        let bits = candidate_bits(&candidate);
        let cell = archive.cell_of(&candidate.descriptor);
        let accepted = archive.add(candidate.solution, candidate.descriptor, candidate.fitness);
        let incumbent = archive
            .elite_at(&descriptor_query)
            .map(elite_bits)
            .expect("the candidate's cell has an incumbent after insertion");
        let archive_elites = archive.elites().map(elite_bits).collect();
        steps.push(StepBits {
            candidate: bits,
            cell,
            accepted,
            incumbent,
            archive_elites,
            coverage: archive.coverage().to_bits(),
            qd_score: archive.qd_score().to_bits(),
        });
    }
    let record = StudyRecord {
        steps,
        capacity: archive.capacity(),
        num_elites: archive.num_elites(),
        coverage: archive.coverage().to_bits(),
        qd_score: archive.qd_score().to_bits(),
        best: archive.best().map(elite_bits),
        final_elites: archive.elites().map(elite_bits).collect(),
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
        if expected != found {
            return Some(format!("steps[{index}]"));
        }
    }
    if expected.capacity != found.capacity {
        return Some("capacity".to_string());
    }
    if expected.num_elites != found.num_elites {
        return Some("num_elites".to_string());
    }
    if expected.coverage != found.coverage {
        return Some("coverage".to_string());
    }
    if expected.qd_score != found.qd_score {
        return Some("qd_score".to_string());
    }
    if expected.best != found.best {
        return Some("best".to_string());
    }
    if expected.final_elites.len() != found.final_elites.len() {
        return Some("final_elites.length".to_string());
    }
    for (elite, (expected, found)) in expected
        .final_elites
        .iter()
        .zip(&found.final_elites)
        .enumerate()
    {
        if expected.solution != found.solution {
            return Some(format!("final_elites[{elite}].solution"));
        }
        if expected.descriptor != found.descriptor {
            return Some(format!("final_elites[{elite}].descriptor"));
        }
        if expected.fitness != found.fitness {
            return Some(format!("final_elites[{elite}].fitness"));
        }
    }
    None
}

fn mutate_final_elite(reference: &SealedStudy) -> CorruptionRun {
    let selector_word = generator_word(MUTATION_SEED, reference.record.final_elites.len(), 0);
    let elite = usize::try_from(selector_word % usize_u64(reference.record.final_elites.len()))
        .expect("elite index fits usize");
    let mantissa_bit = u32::try_from((selector_word >> 32) % 20).expect("mantissa bit fits u32");
    let mut stale = reference.clone();
    let before = stale.record.final_elites[elite].fitness;
    let after = before ^ (1u64 << mantissa_bit);
    stale.record.final_elites[elite].fitness = after;
    let stale_error = stale
        .validate_payload()
        .expect_err("unsealed final-elite mutation must refuse");
    let mutant = SealedStudy::seal(reference.config_digest, stale.record);
    let reference_error = mutant
        .admit_against(reference.output_digest)
        .expect_err("resealed final-elite mutation must miss retained reference");
    let first_mismatch = first_record_mismatch(&reference.record, &mutant.record)
        .expect("one final-elite fitness bit changes the record");
    let semantic_mismatch = mutant
        .record
        .semantic_mismatch()
        .expect("mutated elite must fail independent reconstruction");
    let mutation = Mutation {
        seed: MUTATION_SEED,
        elite,
        mantissa_bit,
        selector_word,
        before,
        after,
    };
    let red_line = format!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"{RED_CASE}\",\"pass\":false,\"config\":\"0x{:016x}\",\"reference\":\"0x{:016x}\",\"mutant\":\"0x{:016x}\",\"corruption_seed\":\"0x{:016x}\",\"selector_word\":\"0x{:016x}\",\"target\":\"final_elites[{}].fitness\",\"mantissa_bit\":{},\"before\":\"0x{:016x}\",\"after\":\"0x{:016x}\",\"stale_gate\":\"PayloadIdentityMismatch\",\"reference_gate\":\"ReferenceIdentityMismatch\",\"first_mismatch\":\"{}\",\"semantic_mismatch\":\"{}\"}}",
        reference.config_digest,
        reference.output_digest,
        mutant.output_digest,
        mutation.seed,
        mutation.selector_word,
        mutation.elite,
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

fn exact_elite_bit_delta(
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
    if expected.final_elites[mutation.elite].fitness != mutation.before {
        return false;
    }
    expected.final_elites[mutation.elite].fitness = mutation.after;
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
#[allow(clippy::too_many_lines)]
fn map_elites_full_study_replays_and_seeded_failure_is_refused() {
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
        "complete candidate, transition, aggregate, and elite frames must repeat bit-for-bit"
    );
    assert_eq!(original, replay);
    assert_mergeable(&original, original.output_digest);

    let first = mutate_final_elite(&original);
    let second = mutate_final_elite(&replay);
    assert_eq!(first, second, "seeded red state must repeat exactly");
    assert!(exact_elite_bit_delta(
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
    let expected_path = format!("final_elites[{}].fitness", first.mutation.elite);
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
        "{{\"suite\":\"{SUITE}\",\"case\":\"{CASE}\",\"pass\":true,\"config\":\"0x{config_digest:016x}\",\"output\":\"0x{:016x}\",\"candidates\":{},\"accepted\":{},\"capacity\":{},\"elites\":{},\"coverage_bits\":\"0x{:016x}\",\"qd_score_bits\":\"0x{:016x}\",\"seed\":\"0x{STUDY_SEED:016x}\",\"scope\":\"same-process finite fixture\"}}",
        original.output_digest,
        original.record.steps.len(),
        original
            .record
            .steps
            .iter()
            .filter(|step| step.accepted)
            .count(),
        original.record.capacity,
        original.record.num_elites,
        original.record.coverage,
        original.record.qd_score,
    );
}
