//! G0/G3/G5 full-public-report replay for the standalone Sinkhorn OT
//! surface (7tv.21.42).
//!
//! A fixed balanced two-point problem is solved at three entropic
//! regularizations. The fixture retains every public `OtReport` field and
//! every transport-plan cell as exact IEEE-754 bits under canonical and
//! domain-separated BLAKE3 identities. An analytic symmetric-plan oracle and
//! independent plan accounting check the retained semantics. A disclosed
//! mutation stream flips one mantissa bit in one off-diagonal plan cell; the
//! stale payload, resealed foreign result, and inconsistent report all refuse.
//!
//! This test covers one finite, symmetric, well-conditioned, same-build
//! fixture. It does not claim arbitrary-measure/cost/epsilon convergence,
//! unbalanced or partial OT, Sinkhorn divergence/debiasing/barycenters,
//! production admission or fallibility, dual-potential/complementarity or
//! iteration-trace replay, an independent elementary-function implementation,
//! cancellation/checkpointing, cross-ISA equality, authenticated-ledger
//! authority, optimizer quality, or performance.

use fs_blake3::{ContentHash, hash_domain};
use fs_dfo::{OtReport, sinkhorn};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-dfo/ot-study-replay";
const CASE: &str = "balanced-two-point-epsilon-ladder";
const RED_CASE: &str = "seeded-off-diagonal-plan-corruption";

const FIXTURE_IDENTITY_KIND: &str = "fs-dfo-ot-study-fixture-v1";
const RESULT_IDENTITY_KIND: &str = "fs-dfo-ot-study-result-v1";
const FIXTURE_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.ot-study-fixture.v1";
const RESULT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.ot-study-result.v1";
const EVENT_DIGEST_DOMAIN: &str = "frankensim.fs-dfo.ot-study-event.v1";

const ROWS: usize = 2;
const COLUMNS: usize = 2;
const PLAN_CELLS: usize = ROWS * COLUMNS;
const A: [f64; ROWS] = [0.5, 0.5];
const B: [f64; COLUMNS] = [0.5, 0.5];
const COSTS: [f64; PLAN_CELLS] = [0.0, 1.0, 1.0, 0.0];
const EPSILONS: [f64; 3] = [1.0, 0.5, 0.25];
const MAX_ITERATIONS: usize = 100;
const EXPECTED_ITERATIONS: usize = 10;
const ANALYTIC_TOLERANCE: f64 = 1.0e-12;
const RESIDUAL_GATE: f64 = 1.0e-12;

const MUTATION_SEED: u64 = 0x07A4_FA11_0000_0042;
const MUTATION_KERNEL: u32 = 0x0A42;
const MUTATION_TILE: u32 = 0;
const OFF_DIAGONAL_CELLS: [usize; 2] = [1, 2];
const MUTATION_BIT_BASE: u32 = 32;
const MUTATION_BIT_COUNT: u64 = 8;

const _: () = assert!(ROWS > 0);
const _: () = assert!(COLUMNS > 0);
const _: () = assert!(PLAN_CELLS == 4);
const _: () = assert!(MAX_ITERATIONS >= EXPECTED_ITERATIONS);
const _: () = assert!(MAX_ITERATIONS % 10 == 0);

#[derive(Debug, Clone, PartialEq, Eq)]
struct OtReportBits {
    cost: u64,
    plan: Vec<u64>,
    marginal_residual: u64,
    iters: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    reports: Vec<OtReportBits>,
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
    report_index: usize,
    plan_index: usize,
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
    u64::try_from(value).expect("fixed OT fixture cardinality fits u64")
}

fn digest_bytes(digest: ContentHash) -> [u8; 32] {
    *digest.as_bytes()
}

fn report_bits(report: OtReport) -> OtReportBits {
    OtReportBits {
        cost: report.cost.to_bits(),
        plan: report.plan.into_iter().map(f64::to_bits).collect(),
        marginal_residual: report.marginal_residual.to_bits(),
        iters: report.iters,
    }
}

fn fixture_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new(FIXTURE_IDENTITY_KIND)
        .str("algorithm", "fs_dfo::sinkhorn/log-domain")
        .str("algorithm-randomness", "none")
        .str("transport-semantics", "balanced-entropic-OT")
        .str("matrix-layout", "row-major")
        .str("mass-units", "unit-probability")
        .str("cost-units", "dimensionless-transport-penalty")
        .u64("rows", usize_u64(ROWS))
        .u64("columns", usize_u64(COLUMNS))
        .u64("plan-cells", usize_u64(PLAN_CELLS))
        .u64("epsilon-count", usize_u64(EPSILONS.len()))
        .u64("max-iterations", usize_u64(MAX_ITERATIONS))
        .u64("residual-check-cadence", 10)
        .str("residual-stop-comparison", "row-residual<1e-10")
        .str(
            "returned-residual",
            "max-row-then-column-absolute-marginal-error",
        )
        .str("reported-cost", "row-major-plan-cost-mul-add")
        .str(
            "analytic-oracle",
            "r=det-exp(-1/epsilon);diag=1/(2*(1+r));offdiag=r/(2*(1+r));cost=r/(1+r)",
        )
        .f64_bits("analytic-tolerance", ANALYTIC_TOLERANCE)
        .f64_bits("returned-residual-gate", RESIDUAL_GATE)
        .str("deterministic-math", "fs-math-det-exp-ln")
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
        .str("fixture-digest-domain", FIXTURE_DIGEST_DOMAIN)
        .str("result-digest-domain", RESULT_DIGEST_DOMAIN)
        .str("event-digest-domain", EVENT_DIGEST_DOMAIN)
        .str(
            "no-claims",
            "arbitrary-inputs;unbalanced-OT;partial-OT;sinkhorn-divergence;debiasing;barycenters;production-admission;dual-potentials;complementarity;iteration-trace;independent-elementary-function-oracle;Cx;checkpoint;cross-ISA;authenticated-ledger;optimizer-quality;performance",
        );
    for (index, mass) in A.into_iter().enumerate() {
        builder = builder
            .u64("row-mass-index", usize_u64(index))
            .f64_bits("row-mass", mass);
    }
    for (index, mass) in B.into_iter().enumerate() {
        builder = builder
            .u64("column-mass-index", usize_u64(index))
            .f64_bits("column-mass", mass);
    }
    for (index, cost) in COSTS.into_iter().enumerate() {
        builder = builder
            .u64("cost-index", usize_u64(index))
            .f64_bits("cost", cost);
    }
    for (index, epsilon) in EPSILONS.into_iter().enumerate() {
        builder = builder
            .u64("epsilon-index", usize_u64(index))
            .f64_bits("epsilon", epsilon);
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
        .u64("report-count", usize_u64(record.reports.len()));
    for (report_index, report) in record.reports.iter().enumerate() {
        builder = builder
            .u64("report-index", usize_u64(report_index))
            .f64_bits("reported-cost", f64::from_bits(report.cost))
            .u64("reported-plan-length", usize_u64(report.plan.len()))
            .f64_bits(
                "reported-marginal-residual",
                f64::from_bits(report.marginal_residual),
            )
            .u64("reported-iterations", usize_u64(report.iters));
        for (plan_index, &bits) in report.plan.iter().enumerate() {
            builder = builder
                .u64("plan-index", usize_u64(plan_index))
                .f64_bits("plan-mass", f64::from_bits(bits));
        }
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
    let reports = EPSILONS
        .into_iter()
        .map(|epsilon| report_bits(sinkhorn(&A, &B, &COSTS, epsilon, MAX_ITERATIONS)))
        .collect();
    let record = StudyRecord { reports };
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

fn analytic_report(epsilon: f64) -> ([f64; PLAN_CELLS], f64) {
    let r = fs_math::det::exp(-1.0 / epsilon);
    let diagonal = 1.0 / (2.0 * (1.0 + r));
    let off_diagonal = r * diagonal;
    (
        [diagonal, off_diagonal, off_diagonal, diagonal],
        r / (1.0 + r),
    )
}

fn plan_cost(plan: &[f64]) -> Option<f64> {
    if plan.len() != PLAN_CELLS {
        return None;
    }
    Some(
        plan.iter()
            .zip(&COSTS)
            .fold(0.0f64, |acc, (&mass, &cost)| mass.mul_add(cost, acc)),
    )
}

fn plan_residual(plan: &[f64]) -> Option<f64> {
    if plan.len() != PLAN_CELLS {
        return None;
    }
    let mut worst = 0.0f64;
    for (row, &target) in plan.chunks_exact(COLUMNS).zip(&A) {
        let mass: f64 = row.iter().sum();
        worst = worst.max((mass - target).abs());
    }
    for (column, &target) in B.iter().enumerate() {
        let mass: f64 = plan.chunks_exact(COLUMNS).map(|row| row[column]).sum();
        worst = worst.max((mass - target).abs());
    }
    Some(worst)
}

fn within(value: f64, expected: f64, tolerance: f64) -> bool {
    (value - expected).abs() <= tolerance
}

#[allow(clippy::too_many_lines)] // All public report fields and analytic/accounting gates meet here.
fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    if record.reports.len() != EPSILONS.len() {
        return Some(format!(
            "report-count:{}!=epsilon-count:{}",
            record.reports.len(),
            EPSILONS.len()
        ));
    }
    let row_mass = A.iter().sum::<f64>();
    let column_mass = B.iter().sum::<f64>();
    if row_mass.to_bits() != 1.0f64.to_bits()
        || column_mass.to_bits() != 1.0f64.to_bits()
        || row_mass.to_bits() != column_mass.to_bits()
    {
        return Some(format!(
            "fixture-mass:rows=0x{:016x};columns=0x{:016x}",
            row_mass.to_bits(),
            column_mass.to_bits()
        ));
    }

    let mut previous_cost = f64::INFINITY;
    let mut previous_off_diagonal = f64::INFINITY;
    for (report_index, (&epsilon, report)) in EPSILONS.iter().zip(&record.reports).enumerate() {
        if report.plan.len() != PLAN_CELLS {
            return Some(format!(
                "report[{report_index}]-plan-length:{}!=expected-{PLAN_CELLS}",
                report.plan.len()
            ));
        }
        let plan: Vec<f64> = report.plan.iter().copied().map(f64::from_bits).collect();
        if plan.iter().any(|mass| !mass.is_finite() || *mass < 0.0) {
            return Some(format!(
                "report[{report_index}]-invalid-plan:{:016x?}",
                report.plan
            ));
        }
        let cost = f64::from_bits(report.cost);
        let residual = f64::from_bits(report.marginal_residual);
        if !cost.is_finite() || cost <= 0.0 {
            return Some(format!(
                "report[{report_index}]-invalid-cost:0x{:016x}",
                report.cost
            ));
        }
        if !residual.is_finite() || residual < 0.0 {
            return Some(format!(
                "report[{report_index}]-invalid-residual:0x{:016x}",
                report.marginal_residual
            ));
        }
        if report.iters != EXPECTED_ITERATIONS {
            return Some(format!(
                "report[{report_index}]-iterations:{}!=expected-{EXPECTED_ITERATIONS}",
                report.iters
            ));
        }

        let Some(recomputed_cost) = plan_cost(&plan) else {
            return Some(format!("report[{report_index}]-cost-shape-refusal"));
        };
        if recomputed_cost.to_bits() != report.cost {
            return Some(format!(
                "report[{report_index}]-reported-cost!=plan-dot-cost:reported=0x{:016x};recomputed=0x{:016x}",
                report.cost,
                recomputed_cost.to_bits()
            ));
        }
        let Some(recomputed_residual) = plan_residual(&plan) else {
            return Some(format!("report[{report_index}]-residual-shape-refusal"));
        };
        if recomputed_residual.to_bits() != report.marginal_residual {
            return Some(format!(
                "report[{report_index}]-reported-residual!=plan-marginals:reported=0x{:016x};recomputed=0x{:016x}",
                report.marginal_residual,
                recomputed_residual.to_bits()
            ));
        }
        if residual > RESIDUAL_GATE {
            return Some(format!(
                "report[{report_index}]-residual:{}>gate:{RESIDUAL_GATE}",
                residual
            ));
        }

        if report.plan[0] != report.plan[3] || report.plan[1] != report.plan[2] {
            return Some(format!(
                "report[{report_index}]-symmetric-fixture-plan-asymmetry"
            ));
        }
        let (analytic_plan, analytic_cost) = analytic_report(epsilon);
        for (plan_index, (&found, &expected)) in plan.iter().zip(&analytic_plan).enumerate() {
            if !within(found, expected, ANALYTIC_TOLERANCE) {
                return Some(format!(
                    "report[{report_index}]-analytic-plan[{plan_index}]:found=0x{:016x};expected=0x{:016x}",
                    found.to_bits(),
                    expected.to_bits()
                ));
            }
        }
        if !within(cost, analytic_cost, ANALYTIC_TOLERANCE) {
            return Some(format!(
                "report[{report_index}]-analytic-cost:found=0x{:016x};expected=0x{:016x}",
                cost.to_bits(),
                analytic_cost.to_bits()
            ));
        }
        let off_diagonal = plan[1] + plan[2];
        if report_index > 0 && (cost >= previous_cost || off_diagonal >= previous_off_diagonal) {
            return Some(format!(
                "report[{report_index}]-epsilon-ladder-not-strictly-decreasing"
            ));
        }
        previous_cost = cost;
        previous_off_diagonal = off_diagonal;
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

fn exact_plan_bit_delta(reference: &StudyRun, mutant: &StudyRun, mutation: Mutation) -> bool {
    let Some(mask) = 1u64.checked_shl(mutation.mantissa_bit) else {
        return false;
    };
    let Some(reference_report) = reference.record.reports.get(mutation.report_index) else {
        return false;
    };
    let Some(mutant_report) = mutant.record.reports.get(mutation.report_index) else {
        return false;
    };
    let Some(&reference_bits) = reference_report.plan.get(mutation.plan_index) else {
        return false;
    };
    let Some(&mutant_bits) = mutant_report.plan.get(mutation.plan_index) else {
        return false;
    };
    if reference.fixture != mutant.fixture
        || reference.fixture_digest != mutant.fixture_digest
        || reference_bits != mutation.before
        || mutant_bits != mutation.after
        || mutation.before ^ mutation.after != mask
    {
        return false;
    }
    let mut expected = reference.record.clone();
    expected.reports[mutation.report_index].plan[mutation.plan_index] = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(reference: &StudyRun) -> SeededCorruption {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let report_index = usize::try_from(selector.next_below(usize_u64(EPSILONS.len())))
        .expect("selected report fits usize");
    let off_diagonal_index =
        usize::try_from(selector.next_below(usize_u64(OFF_DIAGONAL_CELLS.len())))
            .expect("selected off-diagonal slot fits usize");
    let plan_index = OFF_DIAGONAL_CELLS[off_diagonal_index];
    let mantissa_bit = MUTATION_BIT_BASE
        + u32::try_from(selector.next_below(MUTATION_BIT_COUNT)).expect("selected bit fits u32");
    let selector_draws = selector.index();

    let mut run = reference.clone();
    let before = run.record.reports[report_index].plan[plan_index];
    let after = before ^ (1u64 << mantissa_bit);
    run.record.reports[report_index].plan[plan_index] = after;
    let stale_error = validate_payload(&run).expect_err("unsealed OT mutation must refuse");
    reseal(&mut run);
    let reference_error = admit_reference(&run, reference)
        .expect_err("resealed OT mutation must not match retained reference");
    let semantic_error = validate_semantics(&run)
        .expect_err("resealed OT plan mutation must remain semantically invalid");
    SeededCorruption {
        run,
        mutation: Mutation {
            seed: MUTATION_SEED,
            kernel: MUTATION_KERNEL,
            tile: MUTATION_TILE,
            report_index,
            plan_index,
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
    let mut emitter = Emitter::new(SUITE, CASE);
    emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "sinkhorn-full-public-report-replay-receipt".to_string(),
            json: format!(
                concat!(
                    "{{\"fixture_identity\":\"{}\",\"fixture_blake3\":\"{}\",",
                    "\"result_identity\":\"{}\",\"result_blake3\":\"{}\",",
                    "\"algorithm\":\"fs_dfo::sinkhorn\",\"algorithm_seed\":null,",
                    "\"rows\":{},\"columns\":{},\"epsilon_count\":{},",
                    "\"max_iterations\":{},\"actual_iterations\":[{},{},{}],",
                    "\"cost_bits\":[\"0x{:016x}\",\"0x{:016x}\",\"0x{:016x}\"],",
                    "\"residual_bits\":[\"0x{:016x}\",\"0x{:016x}\",\"0x{:016x}\"],",
                    "\"versions\":{{\"fs_dfo\":\"{}\",\"fs_math\":\"{}\",",
                    "\"fs_obs\":\"{}\",\"fs_rand\":\"{}\"}},",
                    "\"no_claims\":[\"arbitrary-inputs\",\"unbalanced-or-partial-OT\",",
                    "\"dual-or-iteration-trace\",\"production-admission\",",
                    "\"cross-ISA\",\"cancellation\",\"checkpointing\",",
                    "\"authenticated-ledger\",\"performance\"]}}"
                ),
                run.fixture.hex(),
                run.fixture_digest.to_hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                ROWS,
                COLUMNS,
                EPSILONS.len(),
                MAX_ITERATIONS,
                run.record.reports[0].iters,
                run.record.reports[1].iters,
                run.record.reports[2].iters,
                run.record.reports[0].cost,
                run.record.reports[1].cost,
                run.record.reports[2].cost,
                run.record.reports[0].marginal_residual,
                run.record.reports[1].marginal_residual,
                run.record.reports[2].marginal_residual,
                fs_dfo::VERSION,
                fs_math::VERSION,
                fs_obs::VERSION,
                fs_rand::VERSION,
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
                "fixture={}; result={}; blake3={}; reports={}; all-plan-cells-bound",
                run.fixture.hex(),
                run.result.hex(),
                run.result_digest.to_hex(),
                run.record.reports.len(),
            ),
            seed: 0,
        },
        None,
    )
}

fn corruption_event(reference: &StudyRun, corruption: &SeededCorruption) -> Event {
    let mutation = corruption.mutation;
    let detail = format!(
        "reference={}; mutant={}; seed=0x{:016x}; kernel=0x{:04x}; tile={}; selector_draws={}; target=reports[{}].plan[{}]; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale={:?}; reference_gate={:?}; semantic_gate={:?}",
        reference.result_digest.to_hex(),
        corruption.run.result_digest.to_hex(),
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_draws,
        mutation.report_index,
        mutation.plan_index,
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
        fs_obs::lint_failure_record(event).expect("OT evidence retains replay inputs");
        fs_obs::validate_line(&event.to_jsonl()).expect("OT evidence is fs-obs wire-valid");
        let receipt = event.content_identity_receipt();
        event
            .admit_content_identity(&receipt)
            .expect("OT evidence content identity admits exactly");
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One causal test spans replay plus all three refusal gates.
fn sinkhorn_full_public_report_replays_and_seeded_failure_is_refused() {
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
        "complete OT result frames must replay byte-for-byte"
    );

    let first_receipt = green_receipt(&original);
    let second_receipt = green_receipt(&replay);
    assert_event_pair(&first_receipt, &second_receipt, "green OT receipt");
    println!("{}", first_receipt.to_jsonl());

    let first_green = green_verdict(&original);
    let second_green = green_verdict(&replay);
    assert_event_pair(&first_green, &second_green, "green OT verdict");
    assert_mergeable(&first_green);
    assert_mergeable(&second_green);
    println!("{}", first_green.to_jsonl());

    let first = seeded_corruption(&original);
    let second = seeded_corruption(&replay);
    assert_eq!(first, second, "seeded OT corruption must replay exactly");
    assert!(
        exact_plan_bit_delta(&original, &first.run, first.mutation),
        "mutation must change exactly one retained plan-cell bit"
    );
    assert_eq!(
        validate_payload(&first.run),
        Ok(()),
        "resealed OT mutation must be internally self-consistent"
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
            if mismatch.contains("reported-cost!=plan-dot-cost")
    ));

    let first_red = corruption_event(&original, &first);
    let second_red = corruption_event(&replay, &second);
    assert_event_pair(&first_red, &second_red, "red OT evidence");
    println!("{}", first_red.to_jsonl());

    let panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("merge gate must refuse seeded OT corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("merge-gate panic carries text");
    assert!(message.contains(RED_CASE));
    assert!(message.contains(&format!("0x{MUTATION_SEED:016x}")));
    assert!(message.contains(&format!(
        "reports[{}].plan[{}]",
        first.mutation.report_index, first.mutation.plan_index
    )));
    assert!(message.contains("ReferenceIdentityMismatch"));
    assert!(message.contains("SemanticInconsistency"));
}
