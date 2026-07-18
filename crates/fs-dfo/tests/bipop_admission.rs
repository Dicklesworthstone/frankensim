//! G0/G3/G4 admission coverage for fallible BIPOP construction (7tv.23.2).
//!
//! This tranche proves callback-free raw-input refusal, conservative checked
//! seed admission, exact one-callback budget behavior, and replay equivalence
//! between the typed target surface and the legacy compatibility spelling. It
//! does not claim callback fault containment, cancellation, resumable state,
//! allocation recovery, or authenticated root-input identity.

#![deny(unsafe_code)]

use fs_dfo::{
    BIPOP_ADMISSION_SCHEMA_VERSION, BipopError, BipopReport, BipopRestartRecord, CmaReport,
    CmaStopReason, admit_bipop, bipop_cmaes, try_bipop_cmaes,
};
use std::cell::{Cell, RefCell};

const ROOT_SEED: u64 = 0xB1_90_00_02;
const RESTART_SEED_STRIDE: u64 = 0x9E37_79B9;

fn sphere(point: &[f64]) -> f64 {
    point.iter().map(|coordinate| coordinate * coordinate).sum()
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

fn assert_bipop_bits(left: &BipopReport, right: &BipopReport) {
    assert_report_bits(&left.best, &right.best);
    assert_eq!(left.schedule, right.schedule);
    assert_eq!(left.total_evals, right.total_evals);
    assert_eq!(left.best_restart(), right.best_restart());
    assert_eq!(left.records().len(), right.records().len());
    for (left, right) in left.records().iter().zip(right.records()) {
        assert_record_bits(left, right);
    }
}

/// G0/G4: malformed raw inputs refuse before the objective sees any point.
#[test]
fn malformed_inputs_and_zero_budget_are_callback_free() {
    let calls = Cell::new(0usize);
    let mut objective = |_point: &[f64]| {
        calls.set(calls.get() + 1);
        0.0
    };

    let error = try_bipop_cmaes(&mut objective, &[], 0.5, 1, None, ROOT_SEED)
        .expect_err("empty start must refuse");
    assert_eq!(error, BipopError::EmptyStart);

    let payload_nan = f64::from_bits(0x7ff8_0000_0000_0042);
    let error = try_bipop_cmaes(&mut objective, &[1.0, payload_nan], 0.5, 1, None, ROOT_SEED)
        .expect_err("non-finite start must refuse");
    assert_eq!(
        error,
        BipopError::NonFiniteStart {
            component: 1,
            bits: payload_nan.to_bits(),
        }
    );

    for sigma in [-0.0, 0.0, f64::INFINITY, f64::NAN] {
        let error = try_bipop_cmaes(&mut objective, &[1.0], sigma, 1, None, ROOT_SEED)
            .expect_err("invalid sigma must refuse");
        assert_eq!(
            error,
            BipopError::InvalidSigma {
                bits: sigma.to_bits(),
            }
        );
    }

    for target in [f64::NEG_INFINITY, f64::INFINITY, payload_nan] {
        let error = try_bipop_cmaes(&mut objective, &[1.0], 0.5, 1, Some(target), ROOT_SEED)
            .expect_err("explicit non-finite target must refuse");
        assert_eq!(
            error,
            BipopError::NonFiniteTarget {
                bits: target.to_bits(),
            }
        );
    }

    let error = try_bipop_cmaes(&mut objective, &[1.0], 0.5, 0, None, ROOT_SEED)
        .expect_err("zero budget must refuse");
    assert_eq!(error, BipopError::ZeroBudget);

    let overflow_seed = u64::MAX - RESTART_SEED_STRIDE + 1;
    let error = try_bipop_cmaes(&mut objective, &[1.0], 0.5, 2, None, overflow_seed)
        .expect_err("conservative seed overflow must refuse");
    assert_eq!(
        error,
        BipopError::SeedRangeOverflow {
            seed: overflow_seed,
            max_restart_ordinal: 1,
        }
    );
    assert_eq!(calls.get(), 0, "no refused input may invoke the objective");
}

/// G0: admission exposes deterministic checked maxima and refuses the first
/// root seed whose conservative second-restart coordinate would wrap.
#[test]
fn admission_and_seed_boundary_are_exact() {
    let boundary_seed = u64::MAX - RESTART_SEED_STRIDE;
    let admitted = admit_bipop(&[2.0, -1.0], 0.75, 2, None, boundary_seed)
        .expect("the exact seed boundary admits");
    assert_eq!(admitted.schema_version(), BIPOP_ADMISSION_SCHEMA_VERSION);
    assert_eq!(admitted.dimension(), 2);
    assert_eq!(admitted.total_budget(), 2);
    assert_eq!(admitted.base_lambda(), 6);
    assert_eq!(admitted.max_large_lambda(), 12);
    assert_eq!(admitted.max_local_budget(), 3_000);
    assert_eq!(admitted.max_restart_ordinal(), 1);
    assert_eq!(admitted.max_matrix_entries(), 4);
    assert_eq!(admitted.max_population_entries(), 24);

    let overflow_seed = boundary_seed + 1;
    assert_eq!(
        admit_bipop(&[2.0, -1.0], 0.75, 2, None, overflow_seed)
            .expect_err("the next root seed must refuse"),
        BipopError::SeedRangeOverflow {
            seed: overflow_seed,
            max_restart_ordinal: 1,
        }
    );

    let one = admit_bipop(&[2.0, -1.0], 0.75, 1, None, u64::MAX)
        .expect("budget one never derives a second seed");
    assert_eq!(one.max_restart_ordinal(), 0);
    assert_eq!(one.max_large_lambda(), one.base_lambda());
}

/// G0/G4: one admitted callback is a complete final partial restart, not a
/// hidden population overshoot or an empty report.
#[test]
fn one_callback_budget_is_exact_and_replayable() {
    let x0 = [2.0, -1.0];
    let seen = RefCell::new(Vec::<Vec<u64>>::new());
    let mut objective = |point: &[f64]| {
        seen.borrow_mut()
            .push(point.iter().map(|value| value.to_bits()).collect());
        sphere(point)
    };
    let report =
        try_bipop_cmaes(&mut objective, &x0, 0.75, 1, None, u64::MAX).expect("one callback admits");
    report
        .validate_ledger()
        .expect("generated ledger validates");

    assert_eq!(seen.borrow().as_slice(), &[x0.map(f64::to_bits).to_vec()]);
    assert_eq!(report.total_evals, 1);
    assert_eq!(report.records().len(), 1);
    assert_eq!(report.best_restart(), 0);
    let record = &report.records()[0];
    assert_eq!(record.ordinal(), 0);
    assert_eq!(record.allocated_budget(), 1);
    assert_eq!((record.trace_start(), record.trace_end()), (0, 1));
    assert_eq!(record.stop_reason(), CmaStopReason::BudgetExhausted);
    assert_eq!(record.report().evals, 1);
    assert_eq!(record.report().generations, 0);
    assert!(!record.report().converged);
    assert_slice_bits(record.start(), &x0);
}

/// G3/G5-local: the checked target spelling replays every retained bit and is
/// a compatibility projection of both legacy finite and legacy no-target use.
#[test]
fn typed_and_legacy_targets_preserve_complete_report_bits() {
    let x0 = [1.25, -0.75];
    let run_typed = |target: Option<f64>| {
        let mut objective = |point: &[f64]| sphere(point);
        try_bipop_cmaes(&mut objective, &x0, 0.5, 20, target, ROOT_SEED).expect("fixture admits")
    };
    let run_legacy = |target: f64| {
        let mut objective = |point: &[f64]| sphere(point);
        bipop_cmaes(&mut objective, &x0, 0.5, 20, target, ROOT_SEED)
    };

    let finite_first = run_typed(Some(-1.0));
    let finite_replay = run_typed(Some(-1.0));
    let finite_legacy = run_legacy(-1.0);
    assert_bipop_bits(&finite_first, &finite_replay);
    assert_bipop_bits(&finite_first, &finite_legacy);

    let none_first = run_typed(None);
    let none_replay = run_typed(None);
    let none_legacy = run_legacy(f64::NEG_INFINITY);
    assert_bipop_bits(&none_first, &none_replay);
    assert_bipop_bits(&none_first, &none_legacy);
}
