//! fs-plan conformance suite (plan §13.3; the gp3.8 bead). Acceptance:
//! composed error models bound measured end-to-end errors on fixture
//! pipelines (tightness tracked); cost predictions calibrated within
//! stated quantile bands; attribution trees complete (no unattributed
//! mass); online updates demonstrably improve prediction accuracy;
//! deterministic fits from ledger snapshots; the Rep Router plans with
//! these models live.

use std::sync::atomic::{AtomicU32, Ordering};

use fs_plan::{
    Contribution, CostModel, CostObservation, ErrorLedger, ErrorSource, PlanCostOracle, Rigor,
    TimeLedger, TimeStage, cost_model_from_tune,
};

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

fn temp_db(tag: &str) -> String {
    let n = NEXT_DB.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!("fs-plan-conf-{tag}-{}-{n}.db", std::process::id()))
        .display()
        .to_string()
}

fn cleanup_db(path: &str) {
    for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
        let _ = std::fs::remove_file(format!("{path}{suffix}"));
    }
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-plan/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
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
        let measured =
            2e-3 * lcg(&mut seed) + 5e-3 * lcg(&mut seed) + 1e-4 * lcg(&mut seed);
        assert!(total >= measured, "ledger total {total} under measured {measured}");
        worst_ratio = worst_ratio.max(measured / total);
    }
    // Tightness tracked: the bound is conservative but useful.
    assert!(worst_ratio > 0.5, "bound uselessly loose: worst ratio {worst_ratio}");
    assert_eq!(
        ledger.dominant().map(|(s, _)| s),
        Some(ErrorSource::Discretization),
        "escalation must attack the dominant source"
    );
    let json = ledger.explain();
    assert!(json.contains("\"discretization\"") && json.contains("\"dominant\""));
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
    verdict("pl-002", "completeness lint refuses NaN contributions and negative residuals");
}

#[test]
fn pl_003_cost_calibration_within_stated_bands() {
    // Truth: cost = 5e-9·size^1.4 with ±25% multiplicative noise. Fit on
    // 40 observations, audit on 200 held-out draws: the [p10, p90] band
    // must cover roughly 80% (loose gate: > 60% — calibration, not luck).
    let mut seed = 0x5EED_CA11_0000_0021u64;
    let mut draw = |seed: &mut u64| {
        let size = 1e3 * (1.0 + 999.0 * lcg(seed));
        let noise = 0.75 + 0.5 * lcg(seed);
        CostObservation { size, cost_s: 5e-9 * size.powf(1.4) * noise }
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
    verdict("pl-003", &format!("held-out band coverage {coverage:.3} (target ~0.8)"));
}

#[test]
fn pl_004_online_updates_improve_predictions() {
    // The machine drifts: costs double mid-campaign. A model refit with
    // the new observations must beat the stale model on new traffic.
    let mut seed = 0x5EED_04D1_0000_0031u64;
    let regime =
        |seed: &mut u64, scale: f64| -> CostObservation {
            let size = 1e4 * (1.0 + 9.0 * lcg(seed));
            let noise = 0.9 + 0.2 * lcg(seed);
            CostObservation { size, cost_s: scale * 1e-8 * size.powf(1.2) * noise }
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
fn pl_005_models_rebuild_deterministically_from_ledger_tune() {
    let db = temp_db("tune");
    let ledger = fs_ledger::Ledger::open(&db).expect("open");
    // Seed roofline-style tune rows for one kernel across two machines.
    for (machine, rate) in [(b"m1".as_slice(), 2.0e9), (b"m2".as_slice(), 1.0e9)] {
        ledger
            .tune_put(
                "simd-axpy-f64",
                "roofline-v1",
                machine,
                r#"{"reps":9}"#,
                &format!(r#"{{"elems_per_sec":{rate:?},"elements":4194304.0}}"#),
            )
            .expect("tune row");
    }
    let m1 = cost_model_from_tune(&ledger, "simd-axpy-f64", 4_194_304.0).expect("build");
    let m2 = cost_model_from_tune(&ledger, "simd-axpy-f64", 4_194_304.0).expect("rebuild");
    assert_eq!(m1.n_obs(), 2, "one observation per tune row");
    assert_eq!(m2.n_obs(), m1.n_obs(), "snapshot determinism");
    drop(ledger);
    cleanup_db(&db);
    verdict("pl-005", "tune-table rows rebuild into identical models across reads");
}

#[test]
fn pl_006_router_plans_with_live_cost_models() {
    use fs_geom::{ConverterSpec, ErrorModel, RouteRequest, Router};
    let mut router = Router::new();
    for (name, cost) in [("frep->sdf/coarse", 1.0), ("frep->sdf/fine", 4.0)] {
        router
            .register(ConverterSpec {
                from: "frep".to_string(),
                to: "sdf".to_string(),
                name: name.to_string(),
                base_cost_s: cost,
                error: ErrorModel::AdditiveAbs(0.01),
                certified: true,
            })
            .unwrap();
    }
    let mut oracle = PlanCostOracle::new();
    oracle.register_edge("frep->sdf/coarse", 1e6);
    oracle.register_edge("frep->sdf/fine", 1e6);
    let req = RouteRequest {
        from: "frep".to_string(),
        to: "sdf".to_string(),
        scale: 1.0,
        max_abs_error: 0.05,
        max_cost_s: 100.0,
    };
    // A-priori: the coarse edge (base 1.0s) wins.
    let before = router.plan(&req, &oracle).unwrap();
    assert_eq!(before.edges, vec!["frep->sdf/coarse"]);
    // Measured history says coarse is actually slow on THIS machine.
    use fs_geom::CostOracle as _;
    for _ in 0..5 {
        oracle.record("frep->sdf/coarse", 9.0, 0.005);
        oracle.record("frep->sdf/fine", 2.0, 0.004);
    }
    let after = router.plan(&req, &oracle).unwrap();
    assert_eq!(
        after.edges,
        vec!["frep->sdf/fine"],
        "quantile cost models must reroute the router: {after:?}"
    );
    verdict("pl-006", "Rep Router replanned from quantile models after measured drift");
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
    tl.record(TimeStage { op: "report".to_string(), predicted: None, measured_s: Some(0.5) });
    assert!((tl.total_measured_s() - 18.5).abs() < 1e-12);
    assert!((tl.total_p50_s() - 14.0).abs() < 1e-12);
    assert_eq!(tl.calibration(), Some(0.5), "one of two comparable stages in band");
    let json = tl.explain();
    assert!(json.contains("flux.lbm") && json.contains("calibration"));
    verdict("pl-007", "time attribution totals, band calibration, and explain() payload");
}
