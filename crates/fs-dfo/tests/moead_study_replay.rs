//! G5 full-study replay for the production MOEA/D path (7tv.21.32).
//!
//! The fixture runs a short four-variable ZDT1 study over eight ordered
//! Das-Dennis weights. It retains every decision/objective callback and every
//! ordered public-front bit, reconstructs the initializer and exact evaluation
//! count, recomputes every objective, and audits returned-front membership and
//! nondominance. Independent runs must reproduce the complete canonical frame.
//! A disclosed `StreamKey` mutation changes one returned-objective mantissa bit,
//! is rejected while stale, is resealed into a self-consistent payload, emits
//! reproducible red fs-obs evidence, and is refused by the retained-reference
//! and test-local merge gates.
//!
//! This fixture does not claim convergence, Pareto-front quality, hypervolume,
//! coverage, diversity, optimizer superiority, broad-input behavior, internal
//! replacement-history visibility, cross-ISA equality, cancellation,
//! persistence, authenticated admission, external-oracle parity, or
//! performance.

use fs_dfo::{Individual, MoeadParams, das_dennis, moead};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/moead-study-replay";
const CASE: &str = "short-zdt1-full-callback-and-front";
const RED_CASE: &str = "seeded-returned-objective-corruption";

const DIMENSION: usize = 4;
const OBJECTIVES: usize = 2;
const REFERENCE_DIVISIONS: usize = 7;
const EXPECTED_WEIGHTS: usize = 8;
const NEIGHBORS: usize = 3;
const MAX_REPLACE: usize = 2;
const GENERATIONS: usize = 4;
const ETA_C: f64 = 20.0;
const ETA_M: f64 = 20.0;
const MUTATION_PROBABILITY: f64 = 0.25;
const LOWER_BOUND: f64 = 0.0;
const UPPER_BOUND: f64 = 1.0;
const STUDY_SEED: u64 = 29;
const OPTIMIZER_STREAM_KERNEL: u32 = 0x0D0E;
const OPTIMIZER_STREAM_TILE: u32 = 0;
const EXPECTED_EVALUATIONS: usize = EXPECTED_WEIGHTS * (GENERATIONS + 1);

const MUTATION_SEED: u64 = 0xD0F0_7E1D_0000_0032;
const MUTATION_KERNEL: u32 = 0xD032;
const MUTATION_TILE: u32 = 0;

const _: () = assert!(DIMENSION == 4);
const _: () = assert!(OBJECTIVES == 2);
const _: () = assert!(EXPECTED_WEIGHTS == 8);
const _: () = assert!(EXPECTED_EVALUATIONS == 40);

#[derive(Debug, Clone, PartialEq, Eq)]
struct EvaluationBits {
    decision: Vec<u64>,
    objectives: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct IndividualBits {
    decision: Vec<u64>,
    objectives: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    evaluations: Vec<EvaluationBits>,
    front: Vec<IndividualBits>,
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
    accounting_mismatch: String,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn parameters() -> MoeadParams {
    MoeadParams {
        neighbors: NEIGHBORS,
        max_replace: MAX_REPLACE,
        generations: GENERATIONS,
        eta_c: ETA_C,
        eta_m: ETA_M,
        p_mut: MUTATION_PROBABILITY,
        seed: STUDY_SEED,
    }
}

fn weights() -> Vec<Vec<f64>> {
    das_dennis(OBJECTIVES, REFERENCE_DIVISIONS)
}

fn zdt1(decision: &[f64]) -> Vec<f64> {
    let f1 = decision[0];
    let g = 1.0 + 9.0 * decision[1..].iter().sum::<f64>() / (decision.len() - 1) as f64;
    vec![f1, g * (1.0 - fs_math::det::sqrt(f1 / g))]
}

fn evaluation_bits(decision: &[f64], objectives: &[f64]) -> EvaluationBits {
    EvaluationBits {
        decision: decision.iter().map(|value| value.to_bits()).collect(),
        objectives: objectives.iter().map(|value| value.to_bits()).collect(),
    }
}

fn individual_bits(individual: &Individual) -> IndividualBits {
    IndividualBits {
        decision: individual.x.iter().map(|value| value.to_bits()).collect(),
        objectives: individual.f.iter().map(|value| value.to_bits()).collect(),
    }
}

fn fixture_identity() -> ReplayIdentity {
    let params = parameters();
    let weights = weights();
    let mut builder = IdentityBuilder::new("fs-dfo-moead-study-fixture-v1")
        .str("algorithm", "fs_dfo::moead-tchebycheff-v1")
        .str("objective", "zdt1-v1")
        .str("coordinate-units", "dimensionless")
        .str("objective-units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .u64("objectives", usize_u64(OBJECTIVES))
        .f64_bits("lower-bound", LOWER_BOUND)
        .f64_bits("upper-bound", UPPER_BOUND)
        .u64("reference-divisions", usize_u64(REFERENCE_DIVISIONS))
        .u64("weight-count", usize_u64(weights.len()))
        .u64("neighbors", usize_u64(params.neighbors))
        .u64("max-replace", usize_u64(params.max_replace))
        .u64("generations", usize_u64(params.generations))
        .f64_bits("eta-c", params.eta_c)
        .f64_bits("eta-m", params.eta_m)
        .f64_bits("per-variable-mutation-probability", params.p_mut)
        .u64("optimizer-input-seed", params.seed)
        .u64(
            "expected-objective-callbacks",
            usize_u64(EXPECTED_EVALUATIONS),
        )
        .u64(
            "optimizer-stream-kernel",
            u64::from(OPTIMIZER_STREAM_KERNEL),
        )
        .u64("optimizer-stream-tile", u64::from(OPTIMIZER_STREAM_TILE))
        .str(
            "evaluation-accounting",
            "weight-count-initializers+generations-times-weight-count-children",
        )
        .str(
            "returned-output",
            "ordered-rank-zero-subset-of-final-population",
        )
        .str("execution-context", "single-threaded-direct-test-no-Cx")
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .u64(
            "fs-rand-stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str(
            "fs-rand-stream-position-domain",
            fs_rand::STREAM_POSITION_IDENTITY_DOMAIN,
        )
        .str(
            "no-claims",
            "convergence;front-quality;hypervolume;coverage;diversity;superiority;all-objectives;all-dimensions;all-configurations;all-seeds;internal-population;ideal-history;neighborhood-history;replacement-counts;cross-ISA;Cx;checkpoint;parallelism;authenticated-ledger;external-oracle;performance",
        );
    for (weight_index, weight) in weights.iter().enumerate() {
        builder = builder
            .u64("weight-index", usize_u64(weight_index))
            .u64("weight-length", usize_u64(weight.len()));
        for (objective, &value) in weight.iter().enumerate() {
            builder = builder
                .u64("weight-objective-index", usize_u64(objective))
                .f64_bits("weight-value", value);
        }
    }
    builder.finish()
}

fn result_identity(fixture: &ReplayIdentity, record: &StudyRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-moead-study-result-v1")
        .child("fixture", fixture)
        .u64(
            "objective-callback-count",
            usize_u64(record.evaluations.len()),
        )
        .u64("returned-front-count", usize_u64(record.front.len()));
    for (evaluation_index, evaluation) in record.evaluations.iter().enumerate() {
        builder = builder
            .u64("evaluation-index", usize_u64(evaluation_index))
            .u64(
                "evaluation-decision-length",
                usize_u64(evaluation.decision.len()),
            )
            .u64(
                "evaluation-objective-length",
                usize_u64(evaluation.objectives.len()),
            );
        for (coordinate, &value) in evaluation.decision.iter().enumerate() {
            builder = builder
                .u64("evaluation-coordinate-index", usize_u64(coordinate))
                .f64_bits("evaluation-coordinate", f64::from_bits(value));
        }
        for (objective, &value) in evaluation.objectives.iter().enumerate() {
            builder = builder
                .u64("evaluation-objective-index", usize_u64(objective))
                .f64_bits("evaluation-objective", f64::from_bits(value));
        }
    }
    for (individual_index, individual) in record.front.iter().enumerate() {
        builder = builder
            .u64("front-individual-index", usize_u64(individual_index))
            .u64(
                "front-decision-length",
                usize_u64(individual.decision.len()),
            )
            .u64(
                "front-objective-length",
                usize_u64(individual.objectives.len()),
            );
        for (coordinate, &value) in individual.decision.iter().enumerate() {
            builder = builder
                .u64("front-coordinate-index", usize_u64(coordinate))
                .f64_bits("front-coordinate", f64::from_bits(value));
        }
        for (objective, &value) in individual.objectives.iter().enumerate() {
            builder = builder
                .u64("front-objective-index", usize_u64(objective))
                .f64_bits("front-objective", f64::from_bits(value));
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
        moead(
            &mut objective,
            DIMENSION,
            (LOWER_BOUND, UPPER_BOUND),
            &weights(),
            &parameters(),
        )
    };
    let record = StudyRecord {
        evaluations,
        front: front.iter().map(individual_bits).collect(),
    };
    let fixture = fixture_identity();
    let result = result_identity(&fixture, &record);
    StudyRun {
        fixture,
        record,
        result,
    }
}

fn expected_initial_decisions() -> Vec<Vec<u64>> {
    let mut stream = StreamKey {
        seed: STUDY_SEED,
        kernel: OPTIMIZER_STREAM_KERNEL,
        tile: OPTIMIZER_STREAM_TILE,
    }
    .stream();
    let span = UPPER_BOUND - LOWER_BOUND;
    (0..EXPECTED_WEIGHTS)
        .map(|_| {
            (0..DIMENSION)
                .map(|_| span.mul_add(stream.next_f64(), LOWER_BOUND).to_bits())
                .collect()
        })
        .collect()
}

fn weights_mismatch() -> Option<String> {
    let weights = weights();
    if weights.len() != EXPECTED_WEIGHTS {
        return Some(format!(
            "weight-count:{}!=expected-{EXPECTED_WEIGHTS}",
            weights.len()
        ));
    }
    for (index, weight) in weights.iter().enumerate() {
        if weight.len() != OBJECTIVES {
            return Some(format!(
                "weight[{index}]-dimension:{}!=expected-{OBJECTIVES}",
                weight.len()
            ));
        }
        let expected_first = index as f64 / REFERENCE_DIVISIONS as f64;
        let expected_second = (REFERENCE_DIVISIONS - index) as f64 / REFERENCE_DIVISIONS as f64;
        if weight[0].to_bits() != expected_first.to_bits()
            || weight[1].to_bits() != expected_second.to_bits()
        {
            return Some(format!(
                "weight[{index}]:actual={:016x?};expected=[0x{:016x},0x{:016x}]",
                weight,
                expected_first.to_bits(),
                expected_second.to_bits()
            ));
        }
    }
    None
}

fn dominates(left: &[u64], right: &[u64]) -> bool {
    let mut strict = false;
    for (&left, &right) in left.iter().zip(right) {
        let left = f64::from_bits(left);
        let right = f64::from_bits(right);
        if left > right {
            return false;
        }
        if left < right {
            strict = true;
        }
    }
    strict
}

#[allow(clippy::too_many_lines)] // Full callback/front accounting is the receipt.
fn accounting_mismatch(record: &StudyRecord) -> Option<String> {
    if let Some(mismatch) = weights_mismatch() {
        return Some(mismatch);
    }
    if record.evaluations.len() != EXPECTED_EVALUATIONS {
        return Some(format!(
            "callback-count:{}!=expected-{EXPECTED_EVALUATIONS}",
            record.evaluations.len()
        ));
    }
    if record.front.is_empty() {
        return Some("returned-front-is-empty".to_string());
    }
    if record.front.len() > EXPECTED_WEIGHTS {
        return Some(format!(
            "returned-front-count:{}>population-{EXPECTED_WEIGHTS}",
            record.front.len()
        ));
    }

    let expected_initial = expected_initial_decisions();
    for (index, expected) in expected_initial.iter().enumerate() {
        if record.evaluations[index].decision != *expected {
            return Some(format!(
                "initializer[{index}]:actual={:016x?};expected={expected:016x?}",
                record.evaluations[index].decision
            ));
        }
    }

    for (index, evaluation) in record.evaluations.iter().enumerate() {
        if evaluation.decision.len() != DIMENSION {
            return Some(format!(
                "evaluation[{index}]-decision-dimension:{}!=expected-{DIMENSION}",
                evaluation.decision.len()
            ));
        }
        if evaluation.objectives.len() != OBJECTIVES {
            return Some(format!(
                "evaluation[{index}]-objective-dimension:{}!=expected-{OBJECTIVES}",
                evaluation.objectives.len()
            ));
        }
        let decision: Vec<f64> = evaluation
            .decision
            .iter()
            .copied()
            .map(f64::from_bits)
            .collect();
        if decision
            .iter()
            .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
        {
            return Some(format!(
                "evaluation[{index}]-decision-outside-box:{:016x?}",
                evaluation.decision
            ));
        }
        let recomputed: Vec<u64> = zdt1(&decision).into_iter().map(f64::to_bits).collect();
        if recomputed != evaluation.objectives {
            return Some(format!(
                "evaluation[{index}]-objectives:actual={:016x?};recomputed={recomputed:016x?}",
                evaluation.objectives
            ));
        }
        if evaluation
            .objectives
            .iter()
            .any(|value| !f64::from_bits(*value).is_finite())
        {
            return Some(format!(
                "evaluation[{index}]-non-finite-objectives:{:016x?}",
                evaluation.objectives
            ));
        }
    }

    for (index, individual) in record.front.iter().enumerate() {
        if individual.decision.len() != DIMENSION {
            return Some(format!(
                "front[{index}]-decision-dimension:{}!=expected-{DIMENSION}",
                individual.decision.len()
            ));
        }
        if individual.objectives.len() != OBJECTIVES {
            return Some(format!(
                "front[{index}]-objective-dimension:{}!=expected-{OBJECTIVES}",
                individual.objectives.len()
            ));
        }
        let decision: Vec<f64> = individual
            .decision
            .iter()
            .copied()
            .map(f64::from_bits)
            .collect();
        if decision
            .iter()
            .any(|value| !value.is_finite() || !(LOWER_BOUND..=UPPER_BOUND).contains(value))
        {
            return Some(format!(
                "front[{index}]-decision-outside-box:{:016x?}",
                individual.decision
            ));
        }
        let recomputed: Vec<u64> = zdt1(&decision).into_iter().map(f64::to_bits).collect();
        if recomputed != individual.objectives {
            return Some(format!(
                "front[{index}]-objectives:actual={:016x?};recomputed={recomputed:016x?}",
                individual.objectives
            ));
        }
        if !record.evaluations.iter().any(|evaluation| {
            evaluation.decision == individual.decision
                && evaluation.objectives == individual.objectives
        }) {
            return Some(format!("front[{index}]-individual-was-never-evaluated"));
        }
    }

    for (left, left_individual) in record.front.iter().enumerate() {
        for (right, right_individual) in record.front.iter().enumerate() {
            if left != right && dominates(&left_individual.objectives, &right_individual.objectives)
            {
                return Some(format!("front[{left}]-dominates-front[{right}]"));
            }
        }
    }
    None
}

fn first_record_mismatch(left: &StudyRecord, right: &StudyRecord) -> Option<String> {
    if left.evaluations.len() != right.evaluations.len() {
        return Some(format!(
            "evaluations.length:{}!={}",
            left.evaluations.len(),
            right.evaluations.len()
        ));
    }
    for (index, (a, b)) in left.evaluations.iter().zip(&right.evaluations).enumerate() {
        if a != b {
            return Some(format!("evaluations[{index}]:left={a:?};right={b:?}"));
        }
    }
    if left.front.len() != right.front.len() {
        return Some(format!(
            "front.length:{}!={}",
            left.front.len(),
            right.front.len()
        ));
    }
    for (individual, (a, b)) in left.front.iter().zip(&right.front).enumerate() {
        if a.decision.len() != b.decision.len() {
            return Some(format!(
                "front[{individual}].x.length:{}!={}",
                a.decision.len(),
                b.decision.len()
            ));
        }
        if let Some((coordinate, (left, right))) = a
            .decision
            .iter()
            .zip(&b.decision)
            .enumerate()
            .find(|(_, (left, right))| left != right)
        {
            return Some(format!(
                "front[{individual}].x[{coordinate}]:0x{left:016x}!=0x{right:016x}"
            ));
        }
        if a.objectives.len() != b.objectives.len() {
            return Some(format!(
                "front[{individual}].f.length:{}!={}",
                a.objectives.len(),
                b.objectives.len()
            ));
        }
        if let Some((objective, (left, right))) = a
            .objectives
            .iter()
            .zip(&b.objectives)
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
    let Some(reference_individual) = reference.record.front.get(mutation.individual) else {
        return false;
    };
    let Some(&reference_bits) = reference_individual.objectives.get(mutation.objective) else {
        return false;
    };
    let Some(mutant_individual) = mutant.record.front.get(mutation.individual) else {
        return false;
    };
    let Some(&mutant_bits) = mutant_individual.objectives.get(mutation.objective) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || reference_bits != mutation.before
        || mutant_bits != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.front[mutation.individual].objectives[mutation.objective] = mutation.after;
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
        .expect_err("resealed mutation must not match the retained reference");
    let first_mismatch = first_record_mismatch(&reference.record, &run.record)
        .expect("seeded mutation changes the returned record");
    let accounting_mismatch = accounting_mismatch(&run.record)
        .expect("seeded returned-objective mutation must fail accounting");
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
        accounting_mismatch,
    }
}

fn corruption_detail(reference: &StudyRun, corruption: &SeededCorruption) -> String {
    format!(
        "fixture={}; reference={}; mutant={}; optimizer_seed=0x{STUDY_SEED:016x}; corruption_seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=front[{}].f[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale_gate={:?}; reference_gate={:?}; first_mismatch={}; accounting_mismatch={}",
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
        corruption.accounting_mismatch,
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
            name: "moead-full-study-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"result_identity\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::moead\",\"objective\":\"zdt1\",",
                    "\"input_seed\":{},\"dimension\":{},\"objectives\":{},\"weights\":{},",
                    "\"generations\":{},\"expected_evaluations\":{},",
                    "\"actual_evaluations\":{},\"returned_front\":{},",
                    "\"stream_semantics_version\":{},",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_math\":\"{}\",",
                    "\"fs_obs\":\"{}\",\"fs_rand\":\"{}\"}},",
                    "\"no_claims\":[\"convergence\",\"front-quality\",",
                    "\"hypervolume\",\"coverage\",\"diversity\",",
                    "\"optimizer-superiority\",\"internal-replacement-history\",",
                    "\"cross-ISA\",\"cancellation\",\"checkpointing\",",
                    "\"authenticated-ledger\",\"external-oracle\",\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.result.hex(),
                STUDY_SEED,
                DIMENSION,
                OBJECTIVES,
                EXPECTED_WEIGHTS,
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
    fs_obs::lint_failure_record(&event).expect("MOEA/D receipt must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("MOEA/D receipt must use the fs-obs wire schema");
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
    fs_obs::lint_failure_record(&event).expect("MOEA/D verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("MOEA/D verdict must use the fs-obs wire schema");
    println!("{line}");
    event
}

fn exercise_seeded_corruption(original: &StudyRun, replay: &StudyRun) {
    let first = seeded_corruption(original);
    let second = seeded_corruption(replay);
    assert_eq!(first, second, "seeded red state must replay exactly");
    assert!(
        exact_returned_bit_delta(original, &first.run, first.mutation),
        "the corruption must change exactly one returned objective bit"
    );
    assert!(
        exact_returned_bit_delta(replay, &second.run, second.mutation),
        "the replay corruption must change exactly one returned objective bit"
    );
    assert!(first.mutation.individual < original.record.front.len());
    assert!(first.mutation.objective < OBJECTIVES);
    assert!((0..20).contains(&first.mutation.mantissa_bit));
    assert!(f64::from_bits(first.mutation.before).is_finite());
    assert!(f64::from_bits(first.mutation.after).is_finite());
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
    let expected_mismatch = format!(
        "front[{}].f[{}]",
        first.mutation.individual, first.mutation.objective
    );
    assert!(first.first_mismatch.starts_with(&expected_mismatch));
    assert!(
        first
            .accounting_mismatch
            .starts_with(&format!("front[{}]-objectives", first.mutation.individual))
    );

    let first_detail = corruption_detail(original, &first);
    let second_detail = corruption_detail(replay, &second);
    assert_eq!(first_detail, second_detail);
    let first_event = failure_event(&first_detail, first.mutation);
    let second_event = failure_event(&second_detail, second.mutation);
    for event in [&first_event, &second_event] {
        fs_obs::lint_failure_record(event)
            .expect("seeded MOEA/D corruption must retain replay inputs");
        fs_obs::validate_line(&event.to_jsonl())
            .expect("seeded MOEA/D corruption must remain wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("red event identity must admit its exact content");
    }
    assert_eq!(first_event, second_event);
    assert_eq!(
        first_event.content_identity().canonical_bytes(),
        second_event.content_identity().canonical_bytes()
    );
    println!("{}", first_event.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_event))
        .expect_err("the merge gate must refuse seeded returned-objective corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&expected_mismatch));
    assert!(message.contains("ReferenceIdentityMismatch"));
}

#[test]
fn moead_full_study_replays_and_seeded_failure_is_refused() {
    let original = run_study();
    let replay = run_study();
    let original_accounting = accounting_mismatch(&original.record);
    let replay_accounting = accounting_mismatch(&replay.record);
    assert_eq!(original_accounting, None, "original accounting failed");
    assert_eq!(replay_accounting, None, "replay accounting failed");
    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));

    let mismatch = first_record_mismatch(&original.record, &replay.record);
    assert_eq!(mismatch, None, "full study replay drifted");
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.result, replay.result);
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete result frames must replay byte-for-byte"
    );

    emit_green_receipt(&original);
    let green = emit_green_verdict(&original);
    assert_mergeable(&green);
    exercise_seeded_corruption(&original, &replay);
}
