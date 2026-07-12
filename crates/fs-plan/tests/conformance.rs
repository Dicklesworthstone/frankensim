//! fs-plan conformance suite (plan §13.3; the gp3.8 bead). Acceptance:
//! composed error models bound measured end-to-end errors on fixture
//! pipelines (tightness tracked); cost predictions calibrated within
//! stated quantile bands; attribution trees complete (no unattributed
//! mass); online updates demonstrably improve prediction accuracy;
//! deterministic fits from ledger snapshots; the Rep Router plans with
//! these models live.

use fs_plan::{
    Contribution, CostModel, CostObservation, ErrorLedger, ErrorSource, LedgerDefect,
    PlanCostOracle, Rigor, TimeLedger, TimeLedgerDefect, TimeStage, TuneModelError,
    cost_model_from_tune,
};
use std::collections::BTreeMap;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-plan/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

struct MeasuredEdgeRunner {
    error_abs: f64,
    cost_s: f64,
}

impl fs_geom::EdgeRunner for MeasuredEdgeRunner {
    fn run(&self, _cx: &fs_exec::Cx<'_>) -> Result<fs_geom::EdgeOutcome, String> {
        fs_geom::EdgeOutcome::estimated(
            self.error_abs,
            0.0,
            self.error_abs,
            fs_evidence::ProvenanceHash::of_bytes(b"fs-plan-conformance-edge"),
            self.cost_s,
        )
        .map_err(|error| error.to_string())
    }
}

fn with_cx<R>(f: impl FnOnce(&fs_exec::Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 1,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn record_edge_measurement(
    spec: &fs_geom::ConverterSpec,
    oracle: &mut PlanCostOracle,
    cost_s: f64,
    error_abs: f64,
    cx: &fs_exec::Cx<'_>,
) {
    let mut router = fs_geom::Router::new();
    router.register(spec.clone()).unwrap();
    let plan = router
        .plan(
            &fs_geom::RouteRequest {
                from: spec.from.clone(),
                to: spec.to.clone(),
                scale: 1.0,
                max_abs_error: 1.0,
                max_cost_s: 100.0,
            },
            oracle,
        )
        .unwrap();
    let runners: BTreeMap<String, Box<dyn fs_geom::EdgeRunner>> = BTreeMap::from([(
        spec.name.clone(),
        Box::new(MeasuredEdgeRunner { error_abs, cost_s }) as Box<dyn fs_geom::EdgeRunner>,
    )]);
    router.execute(&plan, &runners, oracle, cx).unwrap();
}

#[test]
fn pl_001_error_ledger_bounds_fixture_pipeline() {
    // A 3-stage synthetic pipeline with KNOWN per-stage error bounds:
    // geometry conversion (2e-3), discretization (5e-3), solver (1e-4).
    // The true end-to-end error is drawn adversarially under each bound;
    // the composed ledger total must bound the measured sum every time.
    let mut ledger = ErrorLedger::new();
    ledger.attribute(Contribution {
        source: ErrorSource::Geometry,
        label: "frep->sdf/sampled".to_string(),
        abs: 2e-3,
        rigor: Rigor::Certified,
    });
    ledger.attribute(Contribution {
        source: ErrorSource::Discretization,
        label: "cutfem/h=2mm".to_string(),
        abs: 5e-3,
        rigor: Rigor::Estimated,
    });
    ledger.attribute(Contribution {
        source: ErrorSource::Algebraic,
        label: "cg/tol=1e-4".to_string(),
        abs: 1e-4,
        rigor: Rigor::Certified,
    });
    ledger.lint().expect("complete attribution");
    let total = ledger.total();
    let mut seed = 0x5EED_9147_0000_0011u64;
    let mut worst_ratio = 0.0f64;
    for _ in 0..2000 {
        let measured = 2e-3 * lcg(&mut seed) + 5e-3 * lcg(&mut seed) + 1e-4 * lcg(&mut seed);
        assert!(
            total >= measured,
            "ledger total {total} under measured {measured}"
        );
        worst_ratio = worst_ratio.max(measured / total);
    }
    // Tightness tracked: the bound is conservative but useful.
    assert!(
        worst_ratio > 0.5,
        "bound uselessly loose: worst ratio {worst_ratio}"
    );
    assert_eq!(
        ledger.dominant().map(|(s, _)| s),
        Some(ErrorSource::Discretization),
        "escalation must attack the dominant source"
    );
    let json = ledger.explain();
    assert!(json.contains("\"discretization\"") && json.contains("\"dominant\""));
    assert!(json.contains("\"dominant\":\"discretization\""));
    assert!(!json.contains("Some(") && !json.contains("None"));
    verdict(
        "pl-001",
        &format!("composed bound held on 2000 draws; worst tightness ratio {worst_ratio:.3}"),
    );
}

#[test]
fn pl_002_ledger_lint_refuses_silent_error_mass() {
    let mut bad = ErrorLedger::new();
    bad.attribute(Contribution {
        source: ErrorSource::Surrogate,
        label: "fno".to_string(),
        abs: f64::NAN,
        rigor: Rigor::Estimated,
    });
    assert!(bad.lint().is_err(), "NaN error mass must refuse");
    let mut neg = ErrorLedger::new();
    neg.declared_residual = -1.0;
    assert!(neg.lint().is_err(), "negative residual must refuse");
    let mut anonymous = ErrorLedger::new();
    anonymous.attribute(Contribution {
        source: ErrorSource::Geometry,
        label: " \t".to_string(),
        abs: 0.1,
        rigor: Rigor::Certified,
    });
    assert_eq!(anonymous.lint(), Err(LedgerDefect::BlankLabel));
    verdict(
        "pl-002",
        "completeness lint refuses NaN contributions, negative residuals, and anonymous sources",
    );
}

#[test]
fn pl_003_cost_calibration_within_stated_bands() {
    // Truth: cost = 5e-9·size^1.4 with ±25% multiplicative noise. Fit on
    // 40 observations, audit on 200 held-out draws: the [p10, p90] band
    // must cover roughly 80% (loose gate: > 60% — calibration, not luck).
    let mut seed = 0x5EED_CA11_0000_0021u64;
    let draw = |seed: &mut u64| {
        let size = 1e3 * (1.0 + 999.0 * lcg(seed));
        let noise = 0.75 + 0.5 * lcg(seed);
        CostObservation {
            size,
            cost_s: 5e-9 * size.powf(1.4) * noise,
        }
    };
    let train: Vec<CostObservation> = (0..40).map(|_| draw(&mut seed)).collect();
    let test: Vec<CostObservation> = (0..200).map(|_| draw(&mut seed)).collect();
    let model = CostModel::fit(&train).expect("fit");
    let coverage = model.calibration(&test).expect("audit");
    assert!(
        (0.6..=1.0).contains(&coverage),
        "p10-p90 coverage {coverage} outside credible range"
    );
    println!(
        "{{\"suite\":\"fs-plan/conformance\",\"metric\":\"cost_band_coverage\",\
         \"value\":{coverage:.3},\"n_train\":40,\"n_test\":200}}"
    );
    verdict(
        "pl-003",
        &format!("held-out band coverage {coverage:.3} (target ~0.8)"),
    );
}

#[test]
fn pl_004_online_updates_improve_predictions() {
    // The machine drifts: costs double mid-campaign. A model refit with
    // the new observations must beat the stale model on new traffic.
    let mut seed = 0x5EED_04D1_0000_0031u64;
    let regime = |seed: &mut u64, scale: f64| -> CostObservation {
        let size = 1e4 * (1.0 + 9.0 * lcg(seed));
        let noise = 0.9 + 0.2 * lcg(seed);
        CostObservation {
            size,
            cost_s: scale * 1e-8 * size.powf(1.2) * noise,
        }
    };
    let old: Vec<CostObservation> = (0..30).map(|_| regime(&mut seed, 1.0)).collect();
    let new_obs: Vec<CostObservation> = (0..30).map(|_| regime(&mut seed, 2.0)).collect();
    let probes: Vec<CostObservation> = (0..50).map(|_| regime(&mut seed, 2.0)).collect();
    let stale = CostModel::fit(&old).unwrap();
    let mut updated = CostModel::fit(&old).unwrap();
    for &o in &new_obs {
        updated.observe(o).unwrap();
    }
    let stale_err = stale.median_rel_error(&probes).unwrap();
    let updated_err = updated.median_rel_error(&probes).unwrap();
    assert!(
        updated_err < stale_err,
        "online updates must improve accuracy: stale {stale_err:.3} vs updated {updated_err:.3}"
    );
    verdict(
        "pl-004",
        &format!("median rel error improved {stale_err:.3} -> {updated_err:.3} after refit"),
    );
}

#[test]
fn pl_005_tune_loader_refuses_legacy_schema_and_foreign_scope() {
    let ledger = fs_ledger::Ledger::open(":memory:").expect("open");
    let foreign_machine = [1_u8; 40];
    let requested_machine = [2_u8; 40];
    ledger
        .tune_put(
            "simd-axpy-f64",
            "roofline-v1",
            &foreign_machine,
            r#"{"reps":9}"#,
            r#"{"elems_per_sec":2000000000.0,"elements":4194304.0}"#,
        )
        .expect("legacy tune row");
    assert!(matches!(
        cost_model_from_tune(&ledger, "simd-axpy-f64", "roofline-v1", &requested_machine,),
        Err(TuneModelError::MissingRow)
    ));
    assert!(matches!(
        cost_model_from_tune(&ledger, "simd-axpy-f64", "roofline-v1", &foreign_machine,),
        Err(TuneModelError::InvalidReceipt { .. })
    ));
    verdict(
        "pl-005",
        "exact tune lookup ignores foreign machines and refuses legacy ad-hoc JSON",
    );
}

#[test]
fn pl_006_router_plans_with_live_cost_models() {
    use fs_geom::{ConverterSpec, ErrorModel, RouteRequest, Router};
    let mut router = Router::new();
    let coarse = ConverterSpec {
        from: "frep".to_string(),
        to: "sdf".to_string(),
        name: "frep->sdf/coarse".to_string(),
        base_cost_s: 1.0,
        error: ErrorModel::AdditiveAbs(0.01),
        certified: false,
    };
    let fine = ConverterSpec {
        base_cost_s: 4.0,
        name: "frep->sdf/fine".to_string(),
        ..coarse.clone()
    };
    router.register(coarse.clone()).unwrap();
    router.register(fine.clone()).unwrap();
    let mut oracle = PlanCostOracle::new();
    oracle.register_edge(&coarse, 1e6).unwrap();
    oracle.register_edge(&fine, 1e6).unwrap();
    let req = RouteRequest {
        from: "frep".to_string(),
        to: "sdf".to_string(),
        scale: 1.0,
        max_abs_error: 0.05,
        max_cost_s: 100.0,
    };
    // A-priori: the coarse edge (base 1.0s) wins.
    let before = router.plan(&req, &oracle).unwrap();
    assert_eq!(before.edges(), ["frep->sdf/coarse"]);
    // Measured history says coarse is actually slow on THIS machine.
    with_cx(|cx| {
        for _ in 0..5 {
            record_edge_measurement(&coarse, &mut oracle, 9.0, 0.005, cx);
            record_edge_measurement(&fine, &mut oracle, 2.0, 0.004, cx);
        }
    });
    assert_eq!(
        oracle.register_edge(&coarse, 2e6),
        Err(fs_plan::PlanOracleError::ReferenceSizeConflict)
    );
    let after = router.plan(&req, &oracle).unwrap();
    assert_eq!(
        after.edges(),
        ["frep->sdf/fine"],
        "quantile cost models must reroute the router: {after:?}"
    );
    verdict(
        "pl-006",
        "Rep Router replanned from quantile models after measured drift",
    );
}

#[test]
fn pl_007_time_ledger_attributes_and_calibrates() {
    let mut tl = TimeLedger::new();
    tl.record(TimeStage {
        op: "flux.lbm".to_string(),
        predicted: Some((10.0, 12.0, 15.0)),
        measured_s: Some(13.0),
    });
    tl.record(TimeStage {
        op: "ascent.optimize".to_string(),
        predicted: Some((1.0, 2.0, 3.0)),
        measured_s: Some(5.0), // blew its band
    });
    tl.record(TimeStage {
        op: "report".to_string(),
        predicted: None,
        measured_s: Some(0.5),
    });
    assert!((tl.total_measured_s() - 18.5).abs() < 1e-12);
    assert!((tl.total_p50_s() - 14.0).abs() < 1e-12);
    assert_eq!(
        tl.calibration(),
        Some(0.5),
        "one of two comparable stages in band"
    );
    let json = tl.explain();
    assert!(json.contains("flux.lbm") && json.contains("calibration"));
    assert!(json.contains("\"calibration\":0.5"));
    assert!(!json.contains("Some(") && !json.contains("None"));
    verdict(
        "pl-007",
        "time attribution totals, band calibration, and explain() payload",
    );
}

#[test]
fn pl_008_explain_payloads_remain_strict_json_under_hostile_state() {
    let mut valid_error = ErrorLedger::new();
    valid_error.attribute(Contribution {
        source: ErrorSource::Statistical,
        label: "hostile\"label\\line\n\u{0001}".to_string(),
        abs: 0.25,
        rigor: Rigor::Estimated,
    });
    let valid_error_json = valid_error.explain();
    assert!(valid_error_json.contains("hostile\\\"label\\\\line\\n\\u0001"));

    let empty_error = ErrorLedger::new().explain();
    assert!(empty_error.contains("\"dominant\":null"));
    assert!(empty_error.contains("\"valid\":true"));

    let mut overflow = ErrorLedger::new();
    for label in ["first", "second"] {
        overflow.attribute(Contribution {
            source: ErrorSource::Geometry,
            label: label.to_string(),
            abs: f64::MAX,
            rigor: Rigor::Certified,
        });
    }
    assert_eq!(overflow.lint(), Err(LedgerDefect::AggregateOverflow));
    let error_json = overflow.explain();
    assert!(error_json.contains("\"valid\":false"));
    assert!(!error_json.contains("inf") && !error_json.contains("NaN"));

    let mut invalid_time = TimeLedger::new();
    invalid_time.record(TimeStage {
        op: "hostile\"stage\n".to_string(),
        predicted: Some((1.0, f64::NAN, 3.0)),
        measured_s: Some(2.0),
    });
    assert!(matches!(
        invalid_time.lint(),
        Err(TimeLedgerDefect::BadStage { .. })
    ));
    let time_json = invalid_time.explain();
    assert!(time_json.contains("\"valid\":false"));
    assert!(!time_json.contains("NaN") && !time_json.contains("Some("));

    let empty_time = TimeLedger::new().explain();
    assert!(empty_time.contains("\"calibration\":null"));
    assert!(empty_time.contains("\"valid\":true"));

    let mut valid_time = TimeLedger::new();
    valid_time.record(TimeStage {
        op: "hostile\"op\\line\n\u{0001}".to_string(),
        predicted: Some((1.0, 2.0, 3.0)),
        measured_s: Some(2.5),
    });
    let valid_time_json = valid_time.explain();
    assert!(valid_time_json.contains("hostile\\\"op\\\\line\\n\\u0001"));

    let validator = fs_ledger::Ledger::open(":memory:").expect("JSON validator ledger");
    for (ordinal, payload) in [
        &valid_error_json,
        &empty_error,
        &error_json,
        &time_json,
        &empty_time,
        &valid_time_json,
    ]
    .into_iter()
    .enumerate()
    {
        validator
            .append_event(&fs_ledger::EventRow {
                session: None,
                t: i64::try_from(ordinal).expect("small ordinal"),
                kind: "fs-plan.explain-json",
                payload: Some(payload),
            })
            .expect("every explain payload must pass the production JSON validator");
    }
    assert_eq!(validator.table_count("events").unwrap(), 6);
}
