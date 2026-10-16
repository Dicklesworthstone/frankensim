//! G0/G3 coverage for the typed BIPOP restart ledger and production callback
//! trace (7tv.23.1, 7tv.23.3).
//!
//! The report binds each completed CMA restart, exact root-input identity, and
//! every objective query/result bit while keeping legacy summary fields checked
//! projections. Cancellation and pause/resume/fork remain separate scope.

#![deny(unsafe_code)]

use fs_dfo::{
    BIPOP_EVALUATION_SCHEMA_VERSION, BIPOP_REPORT_SCHEMA_VERSION, BIPOP_RESTART_SCHEMA_VERSION,
    BIPOP_ROOT_IDENTITY_KIND, BIPOP_TRACE_IDENTITY_DOMAIN, BIPOP_TRACE_IDENTITY_SCHEMA_VERSION,
    BipopLane, BipopReport, BipopRestartRecord, CmaReport, CmaStopReason, bipop_cmaes,
    try_bipop_cmaes,
};

const ROOT_SEED: u64 = 0xB1_90_00_01;
const RESTART_SEED_STRIDE: u64 = 0x9E37_79B9;

fn sphere(point: &[f64]) -> f64 {
    point.iter().map(|coordinate| coordinate * coordinate).sum()
}

fn run_sphere(total_budget: usize, seed: u64) -> BipopReport {
    let mut objective = |point: &[f64]| sphere(point);
    bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, total_budget, -1.0, seed)
}

fn assert_slice_bits(left: &[f64], right: &[f64]) {
    assert_eq!(left.len(), right.len());
    for (left, right) in left.iter().zip(right) {
        assert_eq!(left.to_bits(), right.to_bits());
    }
}

fn assert_report_bits(left: &CmaReport, right: &CmaReport) {
    assert_slice_bits(&left.x_best, &right.x_best);
    assert_eq!(left.f_best.to_bits(), right.f_best.to_bits());
    assert_eq!(left.evals, right.evals);
    assert_eq!(left.generations, right.generations);
    assert_eq!(left.converged, right.converged);
    assert_eq!(left.sigma.to_bits(), right.sigma.to_bits());
}

fn assert_record_bits(left: &BipopRestartRecord, right: &BipopRestartRecord) {
    assert_eq!(left.schema_version(), right.schema_version());
    assert_eq!(left.ordinal(), right.ordinal());
    assert_eq!(left.lane(), right.lane());
    assert_eq!(left.lambda(), right.lambda());
    assert_eq!(left.allocated_budget(), right.allocated_budget());
    assert_eq!(left.seed(), right.seed());
    assert_slice_bits(left.start(), right.start());
    assert_eq!(left.trace_start(), right.trace_start());
    assert_eq!(left.trace_end(), right.trace_end());
    assert_eq!(left.stop_reason(), right.stop_reason());
    assert_report_bits(left.report(), right.report());
}

fn assert_complete_evidence_bits(left: &BipopReport, right: &BipopReport) {
    assert_eq!(left.schema_version(), right.schema_version());
    let left_root = left.root_inputs();
    let right_root = right.root_inputs();
    assert_slice_bits(left_root.start(), right_root.start());
    assert_eq!(left_root.sigma().to_bits(), right_root.sigma().to_bits());
    assert_eq!(left_root.total_budget(), right_root.total_budget());
    assert_eq!(
        left_root.target().map(f64::to_bits),
        right_root.target().map(f64::to_bits)
    );
    assert_eq!(left_root.seed(), right_root.seed());
    assert_eq!(left_root.identity(), right_root.identity());
    assert_eq!(left.trace_identity(), right.trace_identity());

    let left_trace = left.evaluations().collect::<Vec<_>>();
    let right_trace = right.evaluations().collect::<Vec<_>>();
    assert_eq!(left_trace.len(), right_trace.len());
    for (left, right) in left_trace.iter().zip(right_trace) {
        assert_eq!(left.schema_version(), right.schema_version());
        assert_eq!(left.global_offset(), right.global_offset());
        assert_eq!(left.restart(), right.restart());
        assert_eq!(left.local_offset(), right.local_offset());
        assert_slice_bits(left.point(), right.point());
        assert_eq!(left.objective().to_bits(), right.objective().to_bits());
    }
}

/// G0: ordered intervals partition the aggregate evaluation trace, every run
/// obeys its local cap and population accounting, and the public best report is
/// a bit-exact projection of the explicitly named restart.
#[test]
fn g0_restart_records_partition_the_trace_and_name_the_best_run() {
    let report = run_sphere(20, ROOT_SEED);
    report
        .validate_ledger()
        .expect("generated ledger validates");

    assert_eq!(report.total_evals, 20);
    assert_eq!(report.schema_version(), BIPOP_REPORT_SCHEMA_VERSION);
    assert_eq!(report.records().len(), 2);
    assert_eq!(report.schedule, vec![6, 6]);
    let root = report.root_inputs();
    assert_slice_bits(root.start(), &[2.0, -1.0]);
    assert_eq!(root.sigma().to_bits(), 0.75_f64.to_bits());
    assert_eq!(root.total_budget(), 20);
    assert_eq!(root.target().map(f64::to_bits), Some((-1.0_f64).to_bits()));
    assert_eq!(root.seed(), ROOT_SEED);
    assert_eq!(root.identity().kind(), BIPOP_ROOT_IDENTITY_KIND);

    let trace_identity = report.trace_identity();
    assert_eq!(
        trace_identity.schema_version(),
        BIPOP_TRACE_IDENTITY_SCHEMA_VERSION
    );
    assert_eq!(trace_identity.rows(), report.total_evals);
    assert_eq!(trace_identity.dimension(), 2);
    assert_eq!(
        BIPOP_TRACE_IDENTITY_DOMAIN,
        "frankensim.fs-dfo.bipop-callback-trace.v3"
    );
    let trace = report.evaluations().collect::<Vec<_>>();
    assert_eq!(trace.len(), report.total_evals);
    for (global_offset, evaluation) in trace.iter().enumerate() {
        assert_eq!(evaluation.schema_version(), BIPOP_EVALUATION_SCHEMA_VERSION);
        assert_eq!(evaluation.global_offset(), global_offset);
        assert_eq!(evaluation.point().len(), 2);
        assert_eq!(
            evaluation.objective().to_bits(),
            sphere(evaluation.point()).to_bits()
        );
    }

    let first = &report.records()[0];
    assert_eq!(first.schema_version(), BIPOP_RESTART_SCHEMA_VERSION);
    assert_eq!(first.ordinal(), 0);
    assert_eq!(first.lane(), BipopLane::Large);
    assert_eq!(first.lambda(), 6);
    assert_eq!(first.allocated_budget(), 20);
    assert_eq!(first.seed(), ROOT_SEED);
    assert_slice_bits(first.start(), &[2.0, -1.0]);
    assert_eq!((first.trace_start(), first.trace_end()), (0, 19));
    assert_eq!(first.report().evals, 19);
    assert_eq!(first.report().generations, 3);
    assert_eq!(first.stop_reason(), CmaStopReason::BudgetExhausted);

    let second = &report.records()[1];
    assert_eq!(second.ordinal(), 1);
    assert_eq!(second.lane(), BipopLane::Small);
    assert_eq!(second.lambda(), 6);
    assert_eq!(second.allocated_budget(), 1);
    assert_eq!(second.seed(), ROOT_SEED.wrapping_add(RESTART_SEED_STRIDE));
    assert_eq!((second.trace_start(), second.trace_end()), (19, 20));
    assert_eq!(second.report().evals, 1);
    assert_eq!(second.report().generations, 0);
    assert_eq!(second.stop_reason(), CmaStopReason::BudgetExhausted);

    let best = report.best_record().expect("named best restart exists");
    assert!(report.best_restart() < report.records().len());
    assert_report_bits(&report.best, best.report());
}

/// G0: target and stagnation stops are retained causally rather than inferred
/// from the remaining local budget after the run has returned.
#[test]
fn g0_restart_records_retain_causal_terminal_reasons() {
    let mut at_target = |point: &[f64]| sphere(point);
    let target_report = bipop_cmaes(&mut at_target, &[0.0, 0.0], 0.5, 20, 0.0, ROOT_SEED);
    target_report
        .validate_ledger()
        .expect("initial-target ledger validates");
    let target_record = &target_report.records()[0];
    assert_eq!(target_record.report().evals, 1);
    assert_eq!(target_record.report().generations, 0);
    assert!(target_record.report().converged);
    assert_eq!(target_record.stop_reason(), CmaStopReason::TargetReached);

    let mut constant = |_point: &[f64]| 7.0;
    let stagnated = bipop_cmaes(&mut constant, &[0.5, -0.25], 0.5, 800, -1.0, ROOT_SEED);
    stagnated
        .validate_ledger()
        .expect("stagnation ledger validates");
    let first = &stagnated.records()[0];
    assert_eq!(first.stop_reason(), CmaStopReason::Stagnated);
    assert!(first.report().evals < first.allocated_budget());
    assert!(!first.report().converged);
}

/// G0: exact budget edges admit only complete populations. A final remainder
/// becomes its own restart record instead of disappearing from aggregate
/// accounting.
#[test]
fn g0_budget_edges_retain_complete_and_partial_restart_records() {
    for (budget, expected_evals, expected_records) in
        [(1, vec![1], 1), (7, vec![7], 1), (8, vec![7, 1], 2)]
    {
        let report = run_sphere(budget, ROOT_SEED);
        report
            .validate_ledger()
            .expect("budget-edge ledger validates");
        assert_eq!(report.total_evals, budget);
        assert_eq!(report.records().len(), expected_records);
        assert_eq!(
            report
                .records()
                .iter()
                .map(|record| record.report().evals)
                .collect::<Vec<_>>(),
            expected_evals
        );
    }
}

/// G3: exact objective ties select the earliest restart, and mutations to any
/// still-public compatibility projection are refused by the ledger validator.
#[test]
fn g3_earliest_tie_wins_and_projection_mutations_are_refused() {
    let mut constant = |_point: &[f64]| 7.0;
    let report = bipop_cmaes(&mut constant, &[0.5, -0.25], 0.5, 14, -1.0, ROOT_SEED);
    report
        .validate_ledger()
        .expect("constant-objective ledger validates");
    assert_eq!(report.records().len(), 2);
    assert_eq!(report.best_restart(), 0, "equal totals use earliest wins");
    assert_eq!(report.best.f_best.to_bits(), 7.0_f64.to_bits());

    // `-0.0 < +0.0` is false under ordinary comparison but is strictly less
    // under total_cmp. The eighth callback is the one-point final restart.
    let mut calls = 0usize;
    let mut signed_zero = |_point: &[f64]| {
        calls += 1;
        if calls <= 7 { 0.0 } else { -0.0 }
    };
    let signed_zero_report = bipop_cmaes(
        &mut signed_zero,
        &[0.5, -0.25],
        0.5,
        8,
        f64::NEG_INFINITY,
        ROOT_SEED,
    );
    signed_zero_report
        .validate_ledger()
        .expect("signed-zero ledger validates");
    assert_eq!(signed_zero_report.records().len(), 2);
    assert_eq!(signed_zero_report.best_restart(), 1);
    assert_eq!(
        signed_zero_report.best.f_best.to_bits(),
        (-0.0_f64).to_bits()
    );

    let mut bad_total = report.clone();
    bad_total.total_evals += 1;
    let error = bad_total
        .validate_ledger()
        .expect_err("total mutation refuses");
    assert_eq!(error.restart(), None);
    assert_eq!(error.invariant(), "total-evaluations");

    let mut bad_schedule = report.clone();
    bad_schedule.schedule[0] += 1;
    let error = bad_schedule
        .validate_ledger()
        .expect_err("schedule mutation refuses");
    assert_eq!(error.restart(), Some(0));
    assert_eq!(error.invariant(), "population-schedule");

    let mut bad_best = report.clone();
    bad_best.best.f_best = f64::from_bits(bad_best.best.f_best.to_bits() ^ 1);
    let error = bad_best
        .validate_ledger()
        .expect_err("best mutation refuses");
    assert_eq!(error.restart(), None);
    assert_eq!(error.invariant(), "best-projection");
}

/// G3: once large and small lanes have spent the same retained evaluation
/// budget, the next restart returns to the large lane and doubles lambda.
#[test]
fn g3_equal_lane_budgets_advance_the_large_population_ladder() {
    let mut constant = |_point: &[f64]| 7.0;
    let report = bipop_cmaes(&mut constant, &[0.5, -0.25], 0.5, 1_600, -1.0, ROOT_SEED);
    report
        .validate_ledger()
        .expect("large/small/large ledger validates");
    assert!(report.records().len() >= 3);
    let first = &report.records()[0];
    let second = &report.records()[1];
    let third = &report.records()[2];
    assert_eq!((first.lane(), first.lambda()), (BipopLane::Large, 6));
    assert_eq!((second.lane(), second.lambda()), (BipopLane::Small, 6));
    assert_eq!((third.lane(), third.lambda()), (BipopLane::Large, 12));
    assert_eq!(first.report().evals, second.report().evals);
    assert_eq!(third.ordinal(), 2);
    assert_eq!(
        third.seed(),
        ROOT_SEED.wrapping_add(2 * RESTART_SEED_STRIDE)
    );
    assert_eq!(second.trace_end(), third.trace_start());
}

/// G3/G5 precursor: identical root inputs replay every retained restart field
/// and both compatibility projections bit for bit.
#[test]
fn g3_same_seed_replays_the_complete_restart_ledger() {
    let first = run_sphere(20, ROOT_SEED);
    let second = run_sphere(20, ROOT_SEED);
    first.validate_ledger().expect("first ledger validates");
    second.validate_ledger().expect("second ledger validates");

    assert_eq!(first.best_restart(), second.best_restart());
    assert_eq!(first.schedule, second.schedule);
    assert_eq!(first.total_evals, second.total_evals);
    assert_report_bits(&first.best, &second.best);
    assert_eq!(first.records().len(), second.records().len());
    for (left, right) in first.records().iter().zip(second.records()) {
        assert_record_bits(left, right);
    }
    assert_complete_evidence_bits(&first, &second);
}

/// G0: objective non-finiteness is retained by exact bits without weakening
/// the finite-query invariant or inventing a finiteness claim.
#[test]
fn g0_trace_retains_nonfinite_objective_bits_as_data() {
    let payload_nan = f64::from_bits(0x7ff8_0000_0000_0042);
    let mut objective = |_point: &[f64]| payload_nan;
    let report = try_bipop_cmaes(&mut objective, &[1.0], 0.5, 1, None, ROOT_SEED)
        .expect("non-finite objective output remains retained data");
    report
        .validate_ledger()
        .expect("exact non-finite objective trace validates");
    let evaluation = report.evaluation(0).expect("one callback retained");
    assert_eq!(evaluation.objective().to_bits(), payload_nan.to_bits());
    assert_eq!(report.best.f_best.to_bits(), payload_nan.to_bits());
    assert!(evaluation.point().iter().all(|value| value.is_finite()));
}

/// G0/G3: local CMA best selection uses the same exact total order as the
/// aggregate restart selector. This exercises signed zero inside one complete
/// population; ordinary `<` would incorrectly retain the initial `+0.0`.
#[test]
fn g0_within_restart_signed_zero_uses_earliest_total_cmp_minimum() {
    let mut calls = 0usize;
    let values = [0.0, -0.0, 3.0, 2.0, 1.0];
    let mut objective = |_point: &[f64]| {
        let value = values[calls];
        calls += 1;
        value
    };
    let report = try_bipop_cmaes(&mut objective, &[1.0], 0.5, values.len(), None, ROOT_SEED)
        .expect("one-generation signed-zero fixture admits");
    report
        .validate_ledger()
        .expect("total-ordered signed-zero ledger validates");

    assert_eq!(report.records().len(), 1);
    assert_eq!(report.best.f_best.to_bits(), (-0.0_f64).to_bits());
    let expected = report
        .evaluation(1)
        .expect("negative-zero callback retained");
    assert_eq!(expected.objective().to_bits(), (-0.0_f64).to_bits());
    assert_slice_bits(&report.best.x_best, expected.point());
}

/// G0/G3: positive and negative NaN payloads participate in the exact IEEE-754
/// total order without normalization. A positive NaN initial value, a finite
/// improvement, and then a negative NaN minimum force two updates that ordinary
/// comparison would both suppress.
#[test]
fn g0_within_restart_nan_and_finite_values_use_exact_total_order() {
    let positive_nan = f64::from_bits(0x7ff8_0000_0000_0042);
    let negative_nan = f64::from_bits(0xfff8_0000_0000_0024);
    let values = [positive_nan, 4.0, negative_nan, f64::NEG_INFINITY, 1.0];
    let mut calls = 0usize;
    let mut objective = |_point: &[f64]| {
        let value = values[calls];
        calls += 1;
        value
    };
    let report = try_bipop_cmaes(&mut objective, &[1.0], 0.5, values.len(), None, ROOT_SEED)
        .expect("one-generation NaN/finite fixture admits");
    report
        .validate_ledger()
        .expect("total-ordered NaN/finite ledger validates");

    assert_eq!(report.records().len(), 1);
    assert_eq!(report.best.f_best.to_bits(), negative_nan.to_bits());
    let expected = report
        .evaluation(2)
        .expect("negative-NaN callback retained");
    assert_eq!(expected.objective().to_bits(), negative_nan.to_bits());
    assert_slice_bits(&report.best.x_best, expected.point());
}

/// G0/G3: target stopping is an existential numeric callback witness, not an
/// implication of the total-order representative. A negative NaN remains the
/// exact `f_best`, while a later finite value in the same complete generation
/// still reaches the target and terminates the study.
#[test]
fn g0_finite_target_witness_survives_a_negative_nan_total_order_best() {
    let negative_nan = f64::from_bits(0xfff8_0000_0000_0024);
    let values = [negative_nan, 4.0, 0.25, 3.0, 2.0];
    let mut calls = 0usize;
    let mut objective = |_point: &[f64]| {
        let value = values[calls];
        calls += 1;
        value
    };
    let report = try_bipop_cmaes(
        &mut objective,
        &[1.0],
        0.5,
        values.len(),
        Some(0.5),
        ROOT_SEED,
    )
    .expect("one-generation target-witness fixture admits");
    report
        .validate_ledger()
        .expect("retained finite target witness validates independently of f_best");

    assert_eq!(calls, values.len());
    assert_eq!(report.records().len(), 1);
    assert_eq!(report.best.f_best.to_bits(), negative_nan.to_bits());
    assert!(report.best.converged);
    assert_eq!(
        report.records()[0].stop_reason(),
        CmaStopReason::TargetReached
    );
    assert_eq!(
        report
            .evaluation(2)
            .expect("finite target callback retained")
            .objective()
            .to_bits(),
        0.25_f64.to_bits()
    );
}
