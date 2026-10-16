//! G0/G3/G5 full-public-report replay for the discrete Wasserstein-DRO
//! inner supremum (7tv.21.44).
//!
//! A fixed three-sample/four-support problem has a unique public worst-case
//! distribution and an independently specified primal coupling that exactly
//! exhausts the radius. The fixture retains every public `DroReport` bit and
//! every coupling cell under canonical and domain-separated BLAKE3 identities.
//! An independent tiny-LP enumeration oracle, dual evaluation, and coupling
//! accounting verify feasibility, row and column marginals, transport cost,
//! expected loss, and the active-radius certificate. A disclosed mutation
//! stream flips one mantissa bit in one positive reported mass; stale,
//! identity-consistently resealed, and semantically invalid forms all refuse,
//! and the stable red receipt cannot pass the test-local merge gate.
//!
//! This test covers one finite, feasible, well-conditioned, same-build
//! fixture. It does not claim arbitrary losses/costs/supports, optimizer
//! quality, internal iteration or evaluation accounting (the public report
//! exposes none), cross-ISA equality, cancellation/checkpointing, persistence,
//! authenticated-ledger authority, or performance.

#![deny(unsafe_code)]

use fs_blake3::{ContentHash, hash_domain};
use fs_dfo::{DroReport, wasserstein_worst_case};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/dro-study-replay";
const CASE: &str = "three-sample-four-support-active-radius";
const RED_CASE: &str = "seeded-public-mass-corruption";

const FIXTURE_IDENTITY_KIND: &str = "fs-dfo-dro-study-fixture-v1";
const RESULT_IDENTITY_KIND: &str = "fs-dfo-dro-study-result-v1";
const FIXTURE_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.dro-study-fixture.v1";
const RESULT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.dro-study-result.v1";
const EVENT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.dro-study-event.v1";

const EMPIRICAL_SAMPLES: usize = 3;
const SUPPORT_POINTS: usize = 4;
const COST_CELLS: usize = EMPIRICAL_SAMPLES * SUPPORT_POINTS;
const LOSSES: [f64; SUPPORT_POINTS] = [0.0, 1.0, 3.0, 5.0];
const COSTS: [f64; COST_CELLS] = [
    0.0, 1.0, 2.0, 4.0, // empirical sample 0
    1.0, 0.0, 1.0, 3.0, // empirical sample 1
    2.0, 1.0, 0.0, 2.0, // empirical sample 2
];
const RHO: f64 = 0.5;
const EMPIRICAL_MASS: f64 = 1.0 / 3.0;

// Independent primal witness. Row 1 moves completely to support 2; the
// remaining radius moves one quarter of row 0's mass to support 2.
const PRIMAL_PLAN: [f64; COST_CELLS] = [
    0.25,
    0.0,
    1.0 / 12.0,
    0.0,
    0.0,
    0.0,
    1.0 / 3.0,
    0.0,
    0.0,
    0.0,
    1.0 / 3.0,
    0.0,
];
const EXPECTED_Q: [f64; SUPPORT_POINTS] = [0.25, 0.0, 0.75, 0.0];
const EXPECTED_WORST_CASE: f64 = 2.25;
const EXPECTED_LAMBDA: f64 = 1.5;
const ORACLE_TOLERANCE: f64 = 2.0e-10;
const MARGINAL_TOLERANCE: f64 = 2.0e-12;

const MUTATION_SEED: u64 = 0xD20A_FA11_0000_0044;
const MUTATION_KERNEL: u32 = 0xD044;
const MUTATION_TILE: u32 = 0;
const POSITIVE_Q_CELLS: [usize; 2] = [0, 2];
const MUTATION_BIT_BASE: u32 = 32;
const MUTATION_BIT_COUNT: u64 = 8;

const _: () = assert!(EMPIRICAL_SAMPLES > 0);
const _: () = assert!(SUPPORT_POINTS > 0);
const _: () = assert!(COST_CELLS == 12);
const _: () = assert!(LOSSES.len() == SUPPORT_POINTS);
const _: () = assert!(COSTS.len() == COST_CELLS);
const _: () = assert!(PRIMAL_PLAN.len() == COST_CELLS);
const _: () = assert!(EXPECTED_Q.len() == SUPPORT_POINTS);

#[derive(Debug, Clone, PartialEq, Eq)]
struct DroReportBits {
    worst_case: u64,
    lambda: u64,
    q: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    report: DroReportBits,
    primal_plan: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRun {
    fixture: ReplayIdentity,
    fixture_digest: ContentHash,
    record: StudyRecord,
    result: ReplayIdentity,
    result_digest: ContentHash,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch {
        declared: [u8; 32],
        computed: [u8; 32],
    },
    ReferenceIdentityMismatch {
        expected: [u8; 32],
        found: [u8; 32],
    },
    SemanticInconsistency(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    q_index: usize,
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
    semantic_error: AdmissionError,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixed DRO fixture cardinality fits u64")
}

fn digest_bytes(digest: ContentHash) -> [u8; 32] {
    *digest.as_bytes()
}

fn report_bits(report: DroReport) -> DroReportBits {
    DroReportBits {
        worst_case: report.worst_case.to_bits(),
        lambda: report.lambda.to_bits(),
        q: report.q.into_iter().map(f64::to_bits).collect(),
    }
}

fn fixture_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(FIXTURE_IDENTITY_KIND)
        .str("algorithm", "fs_dfo::wasserstein_worst_case")
        .str("algorithm-randomness", "none")
        .str(
            "problem",
            "discrete-support-Wasserstein-inner-supremum",
        )
        .str(
            "dual-formula",
            "min-lambda>=0 lambda*rho+(1/n)*sum_i max_j(loss_j-lambda*cost_ij)",
        )
        .str("cost-layout", "row-major-empirical-sample-by-support")
        .str("probability-units", "unit-mass")
        .str("loss-units", "dimensionless")
        .str("cost-and-radius-units", "dimensionless-transport-cost")
        .u64("empirical-samples", usize_u64(EMPIRICAL_SAMPLES))
        .u64("support-points", usize_u64(SUPPORT_POINTS))
        .u64("loss-count", usize_u64(LOSSES.len()))
        .u64("cost-cells", usize_u64(COST_CELLS))
        .u64("retained-primal-plan-cells", usize_u64(PRIMAL_PLAN.len()))
        .f64_bits("empirical-row-mass", EMPIRICAL_MASS)
        .f64_bits("wasserstein-radius", RHO)
        .str("radius-admission", "primal-plan-cost<=rho")
        .str("fixture-radius-state", "active-to-tolerance")
        .str("public-report-fields", "worst_case;lambda;q")
        .str(
            "public-work-accounting",
            "unavailable-no-iteration-or-evaluation-counter-claimed",
        )
        .str(
            "independent-oracle",
            "tiny-LP-basic-solution-enumeration-plus-explicit-primal-coupling",
        )
        .f64_bits("expected-worst-case", EXPECTED_WORST_CASE)
        .f64_bits("expected-dual-multiplier", EXPECTED_LAMBDA)
        .f64_bits("oracle-tolerance", ORACLE_TOLERANCE)
        .f64_bits("marginal-tolerance", MARGINAL_TOLERANCE)
        .str("fs-dfo-version", fs_dfo::VERSION)
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
        .u64("mutation-seed", MUTATION_SEED)
        .u64("mutation-kernel", u64::from(MUTATION_KERNEL))
        .u64("mutation-tile", u64::from(MUTATION_TILE))
        .str("fixture-digest-domain", FIXTURE_DIGEST_DOMAIN)
        .str("result-digest-domain", RESULT_DIGEST_DOMAIN)
        .str("event-digest-domain", EVENT_DIGEST_DOMAIN)
        .str(
            "no-claims",
            "arbitrary-instances;optimizer-quality;internal-work-counters;cross-ISA;Cx;checkpoint;persistence;authenticated-ledger;performance",
        );
    for (index, loss) in LOSSES.into_iter().enumerate() {
        builder = builder
            .u64("loss-index", usize_u64(index))
            .f64_bits("loss", loss);
    }
    for (index, cost) in COSTS.into_iter().enumerate() {
        builder = builder
            .u64("cost-index", usize_u64(index))
            .f64_bits("cost", cost);
    }
    for (index, expected_mass) in EXPECTED_Q.into_iter().enumerate() {
        builder = builder
            .u64("expected-q-index", usize_u64(index))
            .f64_bits("expected-q-mass", expected_mass);
    }
    for (index, plan_mass) in PRIMAL_PLAN.into_iter().enumerate() {
        builder = builder
            .u64("oracle-plan-index", usize_u64(index))
            .f64_bits("oracle-plan-mass", plan_mass);
    }
    builder.finish()
}

fn fixture_digest(fixture: &ReplayIdentity) -> ContentHash {
    hash_domain(FIXTURE_DIGEST_DOMAIN, fixture.canonical_bytes())
}

fn result_identity(
    fixture: &ReplayIdentity,
    strong_fixture: ContentHash,
    record: &StudyRecord,
) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(RESULT_IDENTITY_KIND)
        .child("fixture-compatibility-root", fixture)
        .bytes("fixture-canonical-bytes", fixture.canonical_bytes())
        .bytes("fixture-blake3", strong_fixture.as_bytes())
        .f64_bits(
            "reported-worst-case",
            f64::from_bits(record.report.worst_case),
        )
        .f64_bits(
            "reported-dual-multiplier",
            f64::from_bits(record.report.lambda),
        )
        .u64("reported-q-length", usize_u64(record.report.q.len()))
        .u64(
            "retained-primal-plan-length",
            usize_u64(record.primal_plan.len()),
        );
    for (index, &mass) in record.report.q.iter().enumerate() {
        builder = builder
            .u64("reported-q-index", usize_u64(index))
            .f64_bits("reported-q-mass", f64::from_bits(mass));
    }
    for (index, &mass) in record.primal_plan.iter().enumerate() {
        builder = builder
            .u64("retained-plan-index", usize_u64(index))
            .f64_bits("retained-plan-mass", f64::from_bits(mass));
    }
    builder.finish()
}

fn result_digest(result: &ReplayIdentity) -> ContentHash {
    hash_domain(RESULT_DIGEST_DOMAIN, result.canonical_bytes())
}

fn event_digest(event: &Event) -> ContentHash {
    hash_domain(
        EVENT_DIGEST_DOMAIN,
        event.content_identity().canonical_bytes(),
    )
}

fn run_study() -> StudyRun {
    let report = report_bits(wasserstein_worst_case(
        &LOSSES,
        &COSTS,
        EMPIRICAL_SAMPLES,
        RHO,
    ));
    let record = StudyRecord {
        report,
        primal_plan: PRIMAL_PLAN.into_iter().map(f64::to_bits).collect(),
    };
    let fixture = fixture_identity();
    let fixture_digest = fixture_digest(&fixture);
    let result = result_identity(&fixture, fixture_digest, &record);
    let result_digest = result_digest(&result);
    StudyRun {
        fixture,
        fixture_digest,
        record,
        result,
        result_digest,
    }
}

fn within(found: f64, expected: f64, tolerance: f64) -> bool {
    (found - expected).abs() <= tolerance
}

fn primal_plan_marginals(
    plan: &[f64],
) -> Option<([f64; EMPIRICAL_SAMPLES], [f64; SUPPORT_POINTS])> {
    if plan.len() != COST_CELLS {
        return None;
    }
    let mut rows = [0.0f64; EMPIRICAL_SAMPLES];
    let mut columns = [0.0f64; SUPPORT_POINTS];
    for (row_index, row) in plan.chunks_exact(SUPPORT_POINTS).enumerate() {
        rows[row_index] = row.iter().sum();
        for (column_index, &mass) in row.iter().enumerate() {
            columns[column_index] += mass;
        }
    }
    Some((rows, columns))
}

fn primal_plan_cost(plan: &[f64]) -> Option<f64> {
    if plan.len() != COST_CELLS {
        return None;
    }
    Some(
        plan.iter()
            .zip(&COSTS)
            .fold(0.0f64, |acc, (&mass, &cost)| mass.mul_add(cost, acc)),
    )
}

fn distribution_expectation(q: &[f64]) -> Option<f64> {
    if q.len() != SUPPORT_POINTS {
        return None;
    }
    Some(
        q.iter()
            .zip(&LOSSES)
            .fold(0.0f64, |acc, (&mass, &loss)| mass.mul_add(loss, acc)),
    )
}

fn primal_plan_expectation(plan: &[f64]) -> Option<f64> {
    let (_, columns) = primal_plan_marginals(plan)?;
    distribution_expectation(&columns)
}

fn dual_value(lambda: f64) -> f64 {
    let mut acc = 0.0f64;
    for row in COSTS.chunks_exact(SUPPORT_POINTS) {
        let best = LOSSES
            .iter()
            .zip(row)
            .fold(f64::NEG_INFINITY, |best, (&loss, &cost)| {
                best.max(lambda.mul_add(-cost, loss))
            });
        acc += best;
    }
    lambda.mul_add(RHO, acc / EMPIRICAL_SAMPLES as f64)
}

/// Independent primal LP oracle by basic-solution enumeration. Every extreme
/// point assigns each empirical row to one destination, except at most one row
/// that splits between two destinations to bind the radius.
fn primal_oracle() -> f64 {
    let mut best = f64::NEG_INFINITY;
    let mut assignment_count = 1usize;
    for _ in 0..EMPIRICAL_SAMPLES {
        assignment_count = assignment_count
            .checked_mul(SUPPORT_POINTS)
            .expect("fixed DRO oracle assignment count fits usize");
    }
    for code in 0..assignment_count {
        let mut remainder = code;
        let mut cost = 0.0f64;
        let mut value = 0.0f64;
        let mut assignment = [0usize; EMPIRICAL_SAMPLES];
        for (row_index, destination) in assignment.iter_mut().enumerate() {
            *destination = remainder % SUPPORT_POINTS;
            remainder /= SUPPORT_POINTS;
            cost += EMPIRICAL_MASS * COSTS[row_index * SUPPORT_POINTS + *destination];
            value += EMPIRICAL_MASS * LOSSES[*destination];
        }
        if cost <= RHO + MARGINAL_TOLERANCE {
            best = best.max(value);
        }
        for row_index in 0..EMPIRICAL_SAMPLES {
            let source = assignment[row_index];
            for destination in 0..SUPPORT_POINTS {
                if destination == source {
                    continue;
                }
                let delta_cost = EMPIRICAL_MASS
                    * (COSTS[row_index * SUPPORT_POINTS + destination]
                        - COSTS[row_index * SUPPORT_POINTS + source]);
                if delta_cost.abs() < f64::EPSILON {
                    continue;
                }
                let fraction = (RHO - cost) / delta_cost;
                if (0.0..=1.0).contains(&fraction) {
                    let delta_loss = EMPIRICAL_MASS * (LOSSES[destination] - LOSSES[source]);
                    best = best.max(fraction.mul_add(delta_loss, value));
                }
            }
        }
    }
    best
}

#[allow(clippy::too_many_lines)] // Public report, primal coupling, and dual/oracle gates meet here.
fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    if record.report.q.len() != SUPPORT_POINTS {
        return Some(format!(
            "reported-q-length:{}!=support-points:{SUPPORT_POINTS}",
            record.report.q.len()
        ));
    }
    if record.primal_plan.len() != COST_CELLS {
        return Some(format!(
            "primal-plan-length:{}!=cost-cells:{COST_CELLS}",
            record.primal_plan.len()
        ));
    }

    let worst_case = f64::from_bits(record.report.worst_case);
    let lambda = f64::from_bits(record.report.lambda);
    if !worst_case.is_finite() || worst_case.is_sign_negative() {
        return Some(format!(
            "invalid-worst-case:0x{:016x}",
            record.report.worst_case
        ));
    }
    if !lambda.is_finite() || lambda.is_sign_negative() {
        return Some(format!(
            "invalid-dual-multiplier:0x{:016x}",
            record.report.lambda
        ));
    }

    let q: Vec<f64> = record
        .report
        .q
        .iter()
        .copied()
        .map(f64::from_bits)
        .collect();
    for (index, (&mass, &bits)) in q.iter().zip(&record.report.q).enumerate() {
        if !mass.is_finite() || mass.is_sign_negative() {
            return Some(format!("invalid-q[{index}]:0x{bits:016x}"));
        }
    }
    let plan: Vec<f64> = record
        .primal_plan
        .iter()
        .copied()
        .map(f64::from_bits)
        .collect();
    for (index, (&mass, &bits)) in plan.iter().zip(&record.primal_plan).enumerate() {
        if !mass.is_finite() || mass.is_sign_negative() {
            return Some(format!("invalid-primal-plan[{index}]:0x{bits:016x}"));
        }
    }
    let Some((rows, columns)) = primal_plan_marginals(&plan) else {
        return Some("primal-plan-marginal-shape-refusal".to_string());
    };
    for (row_index, row_mass) in rows.into_iter().enumerate() {
        if !within(row_mass, EMPIRICAL_MASS, MARGINAL_TOLERANCE) {
            return Some(format!(
                "primal-row-marginal[{row_index}]:{row_mass:.17e}!=empirical-mass:{EMPIRICAL_MASS:.17e}"
            ));
        }
    }
    for (column_index, (&column_mass, &reported_mass)) in columns.iter().zip(&q).enumerate() {
        if !within(column_mass, reported_mass, MARGINAL_TOLERANCE) {
            return Some(format!(
                "reported-q[{column_index}]-plan-column-marginal:reported={reported_mass:.17e};plan={column_mass:.17e}"
            ));
        }
    }
    let q_sum = q.iter().sum::<f64>();
    if !within(q_sum, 1.0, MARGINAL_TOLERANCE) {
        return Some(format!("q-mass:{q_sum:.17e}!=1"));
    }

    let Some(transport_cost) = primal_plan_cost(&plan) else {
        return Some("primal-plan-cost-shape-refusal".to_string());
    };
    if transport_cost > RHO + MARGINAL_TOLERANCE {
        return Some(format!(
            "primal-plan-cost:{transport_cost:.17e}>radius:{RHO:.17e}"
        ));
    }
    if lambda > ORACLE_TOLERANCE && !within(transport_cost, RHO, MARGINAL_TOLERANCE) {
        return Some(format!(
            "active-radius-not-exhausted:lambda={lambda:.17e};cost={transport_cost:.17e};rho={RHO:.17e}"
        ));
    }

    let Some(q_expectation) = distribution_expectation(&q) else {
        return Some("reported-q-expectation-shape-refusal".to_string());
    };
    let Some(plan_expectation) = primal_plan_expectation(&plan) else {
        return Some("primal-plan-expectation-shape-refusal".to_string());
    };
    if !within(q_expectation, worst_case, ORACLE_TOLERANCE) {
        return Some(format!(
            "reported-q-expectation:{q_expectation:.17e}!=worst-case:{worst_case:.17e}"
        ));
    }
    if !within(plan_expectation, worst_case, ORACLE_TOLERANCE) {
        return Some(format!(
            "primal-plan-expectation:{plan_expectation:.17e}!=worst-case:{worst_case:.17e}"
        ));
    }

    let independent_primal = primal_oracle();
    if !within(worst_case, independent_primal, ORACLE_TOLERANCE) {
        return Some(format!(
            "worst-case:{worst_case:.17e}!=independent-primal:{independent_primal:.17e}"
        ));
    }
    let independent_dual = dual_value(lambda);
    if !within(worst_case, independent_dual, ORACLE_TOLERANCE) {
        return Some(format!(
            "worst-case:{worst_case:.17e}!=dual-at-reported-lambda:{independent_dual:.17e}"
        ));
    }
    if !within(worst_case, EXPECTED_WORST_CASE, ORACLE_TOLERANCE)
        || !within(lambda, EXPECTED_LAMBDA, ORACLE_TOLERANCE)
    {
        return Some(format!(
            "analytic-optimum:worst={worst_case:.17e};lambda={lambda:.17e};expected-worst={EXPECTED_WORST_CASE:.17e};expected-lambda={EXPECTED_LAMBDA:.17e}"
        ));
    }
    for (index, (&found, &expected)) in q.iter().zip(&EXPECTED_Q).enumerate() {
        if !within(found, expected, MARGINAL_TOLERANCE) {
            return Some(format!(
                "analytic-q[{index}]:found={found:.17e};expected={expected:.17e}"
            ));
        }
    }
    None
}

fn validate_payload(run: &StudyRun) -> Result<(), AdmissionError> {
    let expected_fixture = fixture_identity();
    let computed_fixture_digest = fixture_digest(&run.fixture);
    if run.fixture.canonical_bytes() != expected_fixture.canonical_bytes()
        || computed_fixture_digest != run.fixture_digest
    {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: digest_bytes(run.fixture_digest),
            computed: digest_bytes(computed_fixture_digest),
        });
    }
    let computed_result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    let computed_result_digest = result_digest(&computed_result);
    if run.result.canonical_bytes() != computed_result.canonical_bytes()
        || run.result_digest != computed_result_digest
    {
        return Err(AdmissionError::PayloadIdentityMismatch {
            declared: digest_bytes(run.result_digest),
            computed: digest_bytes(computed_result_digest),
        });
    }
    Ok(())
}

fn validate_semantics(run: &StudyRun) -> Result<(), AdmissionError> {
    match semantic_mismatch(&run.record) {
        Some(mismatch) => Err(AdmissionError::SemanticInconsistency(mismatch)),
        None => Ok(()),
    }
}

fn admit_reference(run: &StudyRun, reference: &StudyRun) -> Result<(), AdmissionError> {
    validate_payload(run)?;
    if run.result.canonical_bytes() == reference.result.canonical_bytes()
        && run.result_digest == reference.result_digest
    {
        Ok(())
    } else {
        Err(AdmissionError::ReferenceIdentityMismatch {
            expected: digest_bytes(reference.result_digest),
            found: digest_bytes(run.result_digest),
        })
    }
}

fn reseal(run: &mut StudyRun) {
    run.result = result_identity(&run.fixture, run.fixture_digest, &run.record);
    run.result_digest = result_digest(&run.result);
}

fn exact_q_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let Some(&reference_bits) = reference.record.report.q.get(mutation.q_index) else {
        return false;
    };
    let Some(&mutant_bits) = mutant.record.report.q.get(mutation.q_index) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || reference.fixture_digest != mutant.fixture_digest
        || reference_bits != mutation.before
        || mutant_bits != mutation.after
        || mutation.before ^ mutation.after != mask
        || reference.record.primal_plan != mutant.record.primal_plan
    {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.report.q[mutation.q_index] = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let positive_slot = usize::try_from(selector.next_below(usize_u64(POSITIVE_Q_CELLS.len())))
        .expect("selected positive q slot fits usize");
    let q_index = POSITIVE_Q_CELLS[positive_slot];
    let mantissa_bit = MUTATION_BIT_BASE
        + u32::try_from(selector.next_below(MUTATION_BIT_COUNT)).expect("selected bit fits u32");
    let selector_draws = selector.index();

    let mut run = reference.clone();
    let before = run.record.report.q[q_index];
    let after = before ^ (1u64 << mantissa_bit);
    run.record.report.q[q_index] = after;
    let stale_error = validate_payload(&run).expect_err("unsealed DRO mutation must refuse");
    reseal(&mut run);
    let reference_error = admit_reference(&run, reference)
        .expect_err("resealed DRO mutation must not match retained reference");
    let semantic_error = validate_semantics(&run)
        .expect_err("resealed DRO mass mutation must remain semantically invalid");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            q_index,
            mantissa_bit,
            selector_draws,
            before,
            after,
        },
        stale_error,
        reference_error,
        semantic_error,
    }
}

fn green_receipt(run: &StudyRun) -> Event {
    let plan: Vec<f64> = run
        .record
        .primal_plan
        .iter()
        .copied()
        .map(f64::from_bits)
        .collect();
    let transport_cost = primal_plan_cost(&plan).expect("retained DRO plan has fixed shape");
    let realized_loss = primal_plan_expectation(&plan).expect("retained DRO plan has fixed shape");
    let mut emitter = Emitter::new(SUITE, CASE);
    emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "wasserstein-dro-full-plan-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"fixture_blake3\":\"{}\",",
                    "\"result_identity\":\"{}\",\"result_blake3\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::wasserstein_worst_case\",",
                    "\"algorithm_seed\":null,\"empirical_samples\":{},",
                    "\"support_points\":{},\"loss_count\":{},\"cost_cells\":{},",
                    "\"retained_plan_cells\":{},\"rho_bits\":\"0x{:016x}\",",
                    "\"worst_case_bits\":\"0x{:016x}\",\"lambda_bits\":\"0x{:016x}\",",
                    "\"q_bits\":[\"0x{:016x}\",\"0x{:016x}\",\"0x{:016x}\",\"0x{:016x}\"],",
                    "\"transport_cost_bits\":\"0x{:016x}\",\"realized_loss_bits\":\"0x{:016x}\",",
                    "\"public_work_counters\":null,",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_obs\":\"{}\",\"fs_rand\":\"{}\",",
                    "\"stream_semantics\":{}}},",
                    "\"no_claims\":[\"arbitrary-instances\",\"optimizer-quality\",",
                    "\"internal-work-accounting\",\"cross-ISA\",\"cancellation\",",
                    "\"checkpointing\",\"persistence\",\"authenticated-ledger\",",
                    "\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.fixture_digest.to_hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                EMPIRICAL_SAMPLES,
                SUPPORT_POINTS,
                LOSSES.len(),
                COST_CELLS,
                run.record.primal_plan.len(),
                RHO.to_bits(),
                run.record.report.worst_case,
                run.record.report.lambda,
                run.record.report.q[0],
                run.record.report.q[1],
                run.record.report.q[2],
                run.record.report.q[3],
                transport_cost.to_bits(),
                realized_loss.to_bits(),
                fs_dfo::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
                fs_rand::STREAM_SEMANTICS_VERSION,
            ),
        },
        None,
    )
}

fn green_verdict(run: &StudyRun) -> Event {
    let mut emitter = Emitter::new(SUITE, format!("{CASE}/verdict"));
    emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail: format!(
                "fixture={}; result={}; blake3={}; public-q-cells={}; retained-plan-cells={}; radius-active",
                run.fixture.hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                run.record.report.q.len(),
                run.record.primal_plan.len(),
            ),
            seed: 0,
        },
        None,
    )
}

fn corruption_event(reference: &StudyRun, corruption: &SeededCorruption) -> Event {
    let mutation = corruption.mutation;
    let detail = format!(
        "reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=report.q[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale={:?}; reference_gate={:?}; semantic_gate={:?}",
        reference.result_digest.to_hex(),
        corruption.run.result_digest.to_hex(),
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.q_index,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        corruption.stale_error,
        corruption.reference_error,
        corruption.semantic_error,
    );
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail,
            seed: MUTATION_SEED,
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

fn assert_event_pair(first: &Event, second: &Event, label: &str) {
    assert_eq!(
        first.content_identity().canonical_bytes(),
        second.content_identity().canonical_bytes(),
        "{label} content must replay byte-for-byte"
    );
    assert_eq!(event_digest(first), event_digest(second));
    for event in [first, second] {
        fs_obs::lint_failure_record(event).expect("DRO evidence retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("DRO evidence is fs-obs wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("DRO evidence content identity admits exactly");
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One causal test spans replay plus all refusal gates.
fn wasserstein_dro_full_plan_replays_and_seeded_failure_is_refused() {
    let original = run_study();
    let replay = run_study();
    assert_eq!(validate_payload(&original), Ok(()));
    assert_eq!(validate_payload(&replay), Ok(()));
    assert_eq!(validate_semantics(&original), Ok(()));
    assert_eq!(validate_semantics(&replay), Ok(()));
    assert_eq!(admit_reference(&original, &replay), Ok(()));
    assert_eq!(admit_reference(&replay, &original), Ok(()));
    assert_eq!(original.record, replay.record);
    assert_eq!(original.fixture, replay.fixture);
    assert_eq!(original.fixture_digest, replay.fixture_digest);
    assert_eq!(original.result, replay.result);
    assert_eq!(original.result_digest, replay.result_digest);
    assert_eq!(
        original.result.canonical_bytes(),
        replay.result.canonical_bytes(),
        "complete DRO result frames must replay byte-for-byte"
    );

    let first_receipt = green_receipt(&original);
    let second_receipt = green_receipt(&replay);
    assert_event_pair(&first_receipt, &second_receipt, "green DRO receipt");
    println!("{}", first_receipt.to_jsonl());

    let first_green = green_verdict(&original);
    let second_green = green_verdict(&replay);
    assert_event_pair(&first_green, &second_green, "green DRO verdict");
    assert_mergeable(&first_green);
    assert_mergeable(&second_green);
    println!("{}", first_green.to_jsonl());

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    assert_eq!(first, second, "seeded DRO corruption must replay exactly");
    assert!(
        exact_q_bit_delta(&original, &first.run, first.mutation),
        "mutation must change exactly one retained public q bit"
    );
    assert_eq!(
        validate_payload(&first.run),
        Ok(()),
        "resealed DRO mutation must be internally self-consistent"
    );
    let after = f64::from_bits(first.mutation.after);
    assert!(after.is_finite() && after > 0.0);
    assert!(matches!(
        &first.stale_error,
        AdmissionError::PayloadIdentityMismatch { declared, computed }
            if declared == original.result_digest.as_bytes()
                && computed == first.run.result_digest.as_bytes()
    ));
    assert!(matches!(
        &first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { expected, found }
            if expected == original.result_digest.as_bytes()
                && found == first.run.result_digest.as_bytes()
    ));
    assert!(matches!(
        &first.semantic_error,
        AdmissionError::SemanticInconsistency(mismatch)
            if mismatch.contains("plan-column-marginal")
    ));

    let first_red = corruption_event(&original, &first);
    let second_red = corruption_event(&replay, &second);
    assert_event_pair(&first_red, &second_red, "red DRO evidence");
    println!("{}", first_red.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("merge gate must refuse seeded DRO corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&format!("report.q[{}]", first.mutation.q_index)));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains("SemanticInconsistency"));
}
