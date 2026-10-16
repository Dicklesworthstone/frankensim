//! G3/G5 full-study replay for the production NSGA-II path (7tv.21.34).
//!
//! The fixture retains every ordered ZDT1 objective callback plus every bit of
//! the returned first front. It independently reconstructs the seeded initial
//! population, recomputes objective values, checks exact evaluation accounting,
//! and requires a second run to reproduce the complete frame. A disclosed
//! mutation changes one returned-objective mantissa bit and must fail stale
//! payload, retained-reference, semantic, and test-local merge gates.
//!
//! This target does not claim convergence, hypervolume, coverage, diversity,
//! optimizer superiority, all-seed behavior, cross-ISA equality, cancellation,
//! persistence, authenticated admission, internal selection/variation history,
//! external-oracle parity, or performance.

#![deny(unsafe_code)]

use fs_dfo::{NsgaParams, nsga2};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/nsga2-study-replay";
const CASE: &str = "short-zdt1-complete-callback-and-front";
const RED_CASE: &str = "seeded-returned-objective-corruption";

const DIMENSION: usize = 4;
const OBJECTIVES: usize = 2;
const POPULATION: usize = 8;
const GENERATIONS: usize = 4;
const ETA_C: f64 = 15.0;
const ETA_M: f64 = 20.0;
const MUTATION_PROBABILITY: f64 = 0.25;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND: f64 = 1.0;
const STUDY_SEED: u64 = 41;
const OPTIMIZER_STREAM_KERNEL: u32 = 0x05A2;
const OPTIMIZER_STREAM_TILE: u32 = 0;
const EXPECTED_EVALUATIONS: usize = POPULATION * (GENERATIONS + 1);

const MUTATION_SEED: u64 = 0xD0F0_7E1D_0000_0034;
const MUTATION_KERNEL: u32 = 0xD034;
const MUTATION_TILE: u32 = 0;

const _: () = assert!(DIMENSION > 1);
const _: () = assert!(OBJECTIVES == 2);
const _: () = assert!(POPULATION % 2 == 0);
const _: () = assert!(EXPECTED_EVALUATIONS == 40);

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

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    record: StudyRecord,
    result: ReplayIdentity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
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
struct SeededCorruption {
    run: StudyRun,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    first_mismatch: String,
    semantic_mismatch: String,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
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

fn zdt1(decision: &[f64]) -> Vec<f64> {
    let f1 = decision[0];
    let tail_mean = decision[1..].iter().sum::<f64>() / (decision.len() - 1) as f64;
    let g = 9.0f64.mul_add(tail_mean, 1.0);
    vec![f1, g * (1.0 - fs_math::det::sqrt(f1 / g))]
}

fn evaluation_bits(decision: &[f64], objectives: &[f64]) -> EvaluationBits {
    EvaluationBits {
        decision: decision.iter().map(|value| value.to_bits()).collect(),
        objectives: objectives.iter().map(|value| value.to_bits()).collect(),
    }
}

fn fixture_identity() -> ReplayIdentity {
    IdentityBuilder::new("fs-dfo-nsga2-study-fixture-v1")
        .str("algorithm", "fs_dfo::nsga2")
        .str("objective", "zdt1")
        .str("units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .u64("objectives", usize_u64(OBJECTIVES))
        .u64("population", usize_u64(POPULATION))
        .u64("generations", usize_u64(GENERATIONS))
        .f64_bits("lower-bound", LOWER_BOUND)
        .f64_bits("upper-bound", UPPER_BOUND)
        .f64_bits("eta-c", ETA_C)
        .f64_bits("eta-m", ETA_M)
        .f64_bits("mutation-probability", MUTATION_PROBABILITY)
        .u64("input-seed", STUDY_SEED)
        .u64(
            "optimizer-stream-kernel",
            u64::from(OPTIMIZER_STREAM_KERNEL),
        )
        .u64("optimizer-stream-tile", u64::from(OPTIMIZER_STREAM_TILE))
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .finish()
}

fn result_identity(fixture: &ReplayIdentity, record: &StudyRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-nsga2-study-result-v1")
        .child("fixture", fixture)
        .u64("evaluation-count", usize_u64(record.evaluations.len()));
    for (index, evaluation) in record.evaluations.iter().enumerate() {
        builder = builder.u64("evaluation-index", usize_u64(index));
        for &bits in &evaluation.decision {
            builder = builder.f64_bits("evaluation-decision", f64::from_bits(bits));
        }
        for &bits in &evaluation.objectives {
            builder = builder.f64_bits("evaluation-objective", f64::from_bits(bits));
        }
    }
    builder = builder.u64("front-count", usize_u64(record.front.len()));
    for (index, individual) in record.front.iter().enumerate() {
        builder = builder.u64("front-index", usize_u64(index));
        for &bits in &individual.decision {
            builder = builder.f64_bits("front-decision", f64::from_bits(bits));
        }
        for &bits in &individual.objectives {
            builder = builder.f64_bits("front-objective", f64::from_bits(bits));
        }
    }
    builder.finish()
}

fn run_study() -> StudyRun {
    let mut evaluations = Vec::with_capacity(EXPECTED_EVALUATIONS);
    let front = {
        let mut objective = |decision: &[f64]| {
            let objectives = zdt1(decision);
            evaluations.push(evaluation_bits(decision, &objectives));
            objectives
        };
        nsga2(
            &mut objective,
            DIMENSION,
            (LOWER_BOUND, UPPER_BOUND),
            &parameters(),
        )
    };
    let front = front
        .iter()
        .map(|individual| evaluation_bits(&individual.x, &individual.f))
        .collect();
    let fixture = fixture_identity();
    let record = StudyRecord { evaluations, front };
    let result = result_identity(&fixture, &record);
    StudyRun {
        fixture,
        record,
        result,
    }
}

fn decode(bits: &[u64]) -> Vec<f64> {
    bits.iter().copied().map(f64::from_bits).collect()
}

fn independent_dominates(left: &[u64], right: &[u64]) -> bool {
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

fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    if record.evaluations.len() != EXPECTED_EVALUATIONS {
        return Some(format!(
            "evaluation-count:{}!=expected-{EXPECTED_EVALUATIONS}",
            record.evaluations.len()
        ));
    }
    let mut initializer = StreamKey {
        seed: STUDY_SEED,
        kernel: OPTIMIZER_STREAM_KERNEL,
        tile: OPTIMIZER_STREAM_TILE,
    }
    .stream();
    for (index, evaluation) in record.evaluations.iter().enumerate() {
        if evaluation.decision.len() != DIMENSION {
            return Some(format!(
                "evaluation[{index}]-decision-length:{}!=expected-{DIMENSION}",
                evaluation.decision.len()
            ));
        }
        if evaluation.objectives.len() != OBJECTIVES {
            return Some(format!(
                "evaluation[{index}]-objective-length:{}!=expected-{OBJECTIVES}",
                evaluation.objectives.len()
            ));
        }
        let decision = decode(&evaluation.decision);
        if decision
            .iter()
            .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
        {
            return Some(format!("evaluation[{index}]-decision-domain"));
        }
        let expected = zdt1(&decision);
        let expected_bits: Vec<u64> = expected.iter().map(|value| value.to_bits()).collect();
        if evaluation.objectives != expected_bits {
            return Some(format!("evaluation[{index}]-objectives"));
        }
        if index < POPULATION {
            for coordinate in 0..DIMENSION {
                let expected = (UPPER_BOUND - LOWER_BOUND)
                    .mul_add(initializer.next_f64(), LOWER_BOUND)
                    .to_bits();
                if evaluation.decision[coordinate] != expected {
                    return Some(format!(
                        "initial-population[{index}].x[{coordinate}]:0x{:016x}!=0x{expected:016x}",
                        evaluation.decision[coordinate]
                    ));
                }
            }
        }
    }
    if record.front.is_empty() || record.front.len() > POPULATION {
        return Some(format!(
            "front-count:{} outside 1..={POPULATION}",
            record.front.len()
        ));
    }
    for (individual, member) in record.front.iter().enumerate() {
        if member.decision.len() != DIMENSION || member.objectives.len() != OBJECTIVES {
            return Some(format!("front[{individual}]-dimensions"));
        }
        let decision = decode(&member.decision);
        if decision
            .iter()
            .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
        {
            return Some(format!("front[{individual}]-decision-domain"));
        }
        let expected = zdt1(&decision);
        let expected_bits: Vec<u64> = expected.iter().map(|value| value.to_bits()).collect();
        if member.objectives != expected_bits {
            return Some(format!("front[{individual}]-objectives"));
        }
        if !record
            .evaluations
            .iter()
            .any(|evaluation| evaluation == member)
        {
            return Some(format!("front[{individual}]-missing-callback"));
        }
        for (other, candidate_dominator) in record.front.iter().enumerate() {
            if individual != other
                && independent_dominates(&candidate_dominator.objectives, &member.objectives)
            {
                return Some(format!("front[{individual}]-dominated-by-front[{other}]"));
            }
        }
    }
    None
}

fn first_record_mismatch(left: &StudyRecord, right: &StudyRecord) -> Option<String> {
    if left.evaluations.len() != right.evaluations.len() {
        return Some(format!(
            "evaluation-count:{}!={}",
            left.evaluations.len(),
            right.evaluations.len()
        ));
    }
    for (index, (left, right)) in left.evaluations.iter().zip(&right.evaluations).enumerate() {
        if left != right {
            return Some(format!("evaluation[{index}]"));
        }
    }
    if left.front.len() != right.front.len() {
        return Some(format!(
            "front-count:{}!={}",
            left.front.len(),
            right.front.len()
        ));
    }
    for (individual, (left, right)) in left.front.iter().zip(&right.front).enumerate() {
        if left.decision != right.decision {
            return Some(format!("front[{individual}].x"));
        }
        if left.objectives.len() != right.objectives.len() {
            return Some(format!(
                "front[{individual}].f.length:{}!={}",
                left.objectives.len(),
                right.objectives.len()
            ));
        }
        if let Some((objective, (left, right))) = left
            .objectives
            .iter()
            .zip(&right.objectives)
            .enumerate()
            .find(|(_, (left, right))| left != right)
        {
            return Some(format!(
                "front[{individual}].f[{objective}]:0x{left:016x}!=0x{right:016x}"
            ));
        }
    }
    None
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let computed = result_identity(&run.fixture, &run.record);
    if computed.canonical_bytes() == run.result.canonical_bytes() {
        Ok(())
    } else {
        Err(AdmissionError::PayloadIdentityMismatch {
            declared: run.result.root(),
            computed: computed.root(),
        })
    }
}

fn admit_against(run: &StudyRun, reference: &ReplayIdentity) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result.canonical_bytes() == reference.canonical_bytes() {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: reference.root(),
            found: run.result.root(),
        })
    }
}

fn exact_returned_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let mut expected = reference.record.clone();
    let Some(individual) = expected.front.get_mut(mutation.individual) else {
        return false;
    };
    let Some(bits) = individual.objectives.get_mut(mutation.objective) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || *bits != mutation.before
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }
    *bits = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
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

    let mut run = reference.clone();
    let before = run.record.front[individual].objectives[objective];
    let after = before ^ (1u64 << mantissa_bit);
    run.record.front[individual].objectives[objective] = after;
    let stale_error = validate_payload(&run).expect_err("unsealed result mutation must refuse");
    run.result = result_identity(&run.fixture, &run.record);
    let reference_error = admit_against(&run, &reference.result)
        .expect_err("resealed mutation must not match retained reference");
    let first_mismatch = first_record_mismatch(&reference.record, &run.record)
        .expect("seeded mutation changes the returned record");
    let semantic_mismatch =
        semantic_mismatch(&run.record).expect("mutated objective must fail recomputation");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            individual,
            objective,
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

fn corruption_detail(reference: &StudyRun, corruption: &SeededCorruption) -> String {
    format!(
        "fixture={}; reference={}; mutant={}; optimizer_seed=0x{STUDY_SEED:016x}; corruption_seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=front[{}].f[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale_gate={:?}; reference_gate={:?}; first_mismatch={}; semantic_mismatch={}",
        reference.fixture.hex(),
        reference.result.hex(),
        corruption.run.result.hex(),
        corruption.mutation.seed,
        corruption.mutation.kernel,
        corruption.mutation.tile,
        corruption.mutation.selector_draws,
        corruption.mutation.individual,
        corruption.mutation.objective,
        corruption.mutation.mantissa_bit,
        corruption.mutation.before,
        corruption.mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.first_mismatch,
        corruption.semantic_mismatch,
    )
}

fn failure_event(detail: &str, mutation: Mutation) -> Event {
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail: detail.to_string(),
            seed: mutation.seed,
        },
        None,
    )
}

fn assert_mergeable(event: &Event) {
    let EventKind::ConformanceCase {
        case, pass, detail, ..
    } = &event.kind
    else {
        panic!("merge gate accepts only ConformanceCase evidence");
    };
    assert!(*pass, "merge gate refused {case}: {detail}");
}

fn emit_green_receipt(run: &StudyRun) {
    let mut emitter = Emitter::new(SUITE, CASE);
    let event = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "nsga2-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"result_identity\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::nsga2\",\"objective\":\"zdt1\",",
                    "\"input_seed\":{},\"dimension\":{},\"objectives\":{},",
                    "\"population\":{},\"generations\":{},",
                    "\"expected_evaluations\":{},\"actual_evaluations\":{},",
                    "\"returned_front\":{},\"stream_semantics_version\":{},",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_math\":\"{}\",",
                    "\"fs_obs\":\"{}\",\"fs_rand\":\"{}\"}},",
                    "\"no_claims\":[\"convergence\",\"hypervolume\",",
                    "\"coverage\",\"diversity\",\"optimizer-superiority\",",
                    "\"all-seeds\",\"cross-ISA\",\"cancellation\",",
                    "\"persistence\",\"authenticated-admission\",",
                    "\"internal-selection-variation-history\",",
                    "\"external-oracle\",\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.result.hex(),
                STUDY_SEED,
                DIMENSION,
                OBJECTIVES,
                POPULATION,
                GENERATIONS,
                EXPECTED_EVALUATIONS,
                run.record.evaluations.len(),
                run.record.front.len(),
                fs_rand::STREAM_SEMANTICS_VERSION,
                fs_dfo::VERSION,
                fs_math::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
            ),
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("NSGA-II receipt must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("NSGA-II receipt must use the fs-obs wire schema");
    let receipt = event.content_identity_receipt();
    event
        .admit_content_identity(&receipt)
        .expect("fresh receipt identity must admit exactly");
    println!("{line}");
}

fn emit_green_verdict(run: &StudyRun) -> Event {
    let detail = format!(
        "fixture={}; result={}; callbacks={}; returned_front={}; complete_callback_and_front=bit-exact",
        run.fixture.hex(),
        run.result.hex(),
        run.record.evaluations.len(),
        run.record.front.len(),
    );
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    let event = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail,
            seed: STUDY_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("NSGA-II verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("NSGA-II verdict must use the fs-obs wire schema");
    println!("{line}");
    event
}

fn exercise_seeded_corruption(original: &StudyRun, replay: &StudyRun) {
    let first = seeded_corruption(original);
    let second = seeded_corruption(replay);
    assert_eq!(first, second, "seeded red state must replay exactly");
    assert!(
        exact_returned_bit_delta(original, &first.run, first.mutation),
        "corruption must change exactly one returned objective bit"
    );
    assert!(matches!(
        first.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if declared == original.result.root() && computed == first.run.result.root()
    ));
    assert!(matches!(
        first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.result.root() && found == first.run.result.root()
    ));
    assert!(first.first_mismatch.starts_with(&format!(
        "front[{}].f[{}]",
        first.mutation.individual, first.mutation.objective
    )));
    assert_eq!(
        first.semantic_mismatch,
        format!("front[{}]-objectives", first.mutation.individual)
    );

    let first_detail = corruption_detail(original, &first);
    let second_detail = corruption_detail(replay, &second);
    assert_eq!(first_detail, second_detail);
    let first_event = failure_event(&first_detail, first.mutation);
    let second_event = failure_event(&second_detail, second.mutation);
    for event in [&first_event, &second_event] {
        fs_obs::lint_failure_record(event)
            .expect("seeded NSGA-II corruption must retain replay inputs");
        fs_obs::validate_line(&event.to_jsonl())
            .expect("seeded NSGA-II corruption must remain wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("red event identity must admit its exact content");
    }
    assert_eq!(first_event, second_event);
    assert_eq!(first_event.to_jsonl(), second_event.to_jsonl());
    assert_eq!(
        first_event.content_identity().canonical_bytes(),
        second_event.content_identity().canonical_bytes()
    );
    println!("{}", first_event.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_event))
        .expect_err("merge gate must refuse seeded returned-bit corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains(&format!(
        "front[{}].f[{}]",
        first.mutation.individual, first.mutation.objective
    )));
    assert!(message.contains(&first.semantic_mismatch));
}

#[test]
fn nsga2_full_study_replays_and_seeded_failure_is_refused() {
    let original = run_study();
    let replay = run_study();
    assert_eq!(semantic_mismatch(&original.record), None);
    assert_eq!(semantic_mismatch(&replay.record), None);
    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));
    assert_eq!(admit_against(&original, &original.result), Ok(()));
    assert_eq!(admit_against(&replay, &original.result), Ok(()));
    assert_eq!(
        first_record_mismatch(&original.record, &replay.record),
        None,
        "full NSGA-II study replay drifted"
    );
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.result, replay.result);
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete callback and returned-front frames must replay byte-for-byte"
    );

    emit_green_receipt(&original);
    let green = emit_green_verdict(&original);
    assert_mergeable(&green);
    exercise_seeded_corruption(&original, &replay);
}
