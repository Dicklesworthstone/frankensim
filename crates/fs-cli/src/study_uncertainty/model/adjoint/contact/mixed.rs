//! Joint-state and fixed-wall uncertainty through real native child solves.
use super::*;

fn source(qmc: bool) -> String {
    let mut text = SOURCE.replace(":high 6W)))", concat!(":high 6W)\n",
        " (uniform :name \"wall\" :target fixed-temperature :entity \"cold\" :low 280K :high 300K)))"));
    assert!(text.contains(":target fixed-temperature"));
    if qmc {
        text = text.replace("monte-carlo", "quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)")
            .replace("independent", "(gaussian-copula :latent-correlation ((1 0.25 -0.2) (0.25 1 0.1) (-0.2 0.1 1)))");
    }
    text
}
fn controlled(text: &str) -> String {
    text.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)")
}

#[test]
fn mixed_calibration_uses_one_solve_and_preserves_raw_native_observations() {
    assert!(crate::SOLVE_DRIVER_VERSION >= 44, "this test must execute state-aware contact physics");
    for qmc in [false, true] {
        let f = Fixture::new();
        let raw = source(qmc);
        let text = controlled(&raw);
        let (a, ap) = f.run(&raw, "raw-mixed", crate::exit::SUCCESS);
        let (b, bp) = f.run(&text, "controlled-mixed", crate::exit::SUCCESS);
        let a = Fixture::report(&a, &ap);
        let b = Fixture::report(&b, &bp);
        assert_eq!(a.get("observations"), b.get("observations"));
        let control = b.get("mean_control").unwrap();
        assert_eq!(control.f64_field("probe_solves_planned"), Some(1.0));
        assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
        assert_eq!(control.f64_field("probe_solves_attempted"), Some(1.0));
        let units: Vec<_> = control.get("parameter_units").unwrap().as_array().unwrap().iter()
            .map(|v| v.as_str().unwrap()).collect();
        assert_eq!(units, ["Pa", "W", "K"]);
        let gradient: Vec<_> = control.get("gradient").unwrap().as_array().unwrap().iter()
            .map(|v| v.as_f64().unwrap()).collect();
        assert_eq!(gradient.len(), 3);
        assert!(gradient[0] < 0.0 && gradient[1] > 0.0 && gradient[2] > 0.0);
        let model = Model::load(&f.inputs.join("controlled-mixed.fsim")).unwrap();
        let view = model.nominal_view().unwrap();
        let outputs = view.bound.base().outputs.as_ref().unwrap();
        assert_eq!(outputs.iter().filter(|r| r.name.ends_with("-adjoint")).count(), 1);
        assert!(outputs.iter().any(|r| r.name == "temperature-max-contact-boundary-adjoint"));
        let ledger = Ledger::open(bp.to_str().unwrap()).unwrap();
        let gate = CancelGate::new_clock_free();
        let mean = [900000.0, 5.0, 290.0];
        // These are deterministic derivative probes, never probability samples.
        // The base fixture is at a different pressure AND wall temperature.
        for (direction, step) in [([1.0, 0.0, 0.0], 10.0), ([0.0, 0.0, 1.0], 0.002),
            ([2000.0, 0.2, 0.5], 0.002)] {
            let lo: Vec<_> = mean.iter().zip(direction).map(|(m, d)| m - step * d).collect();
            let hi: Vec<_> = mean.iter().zip(direction).map(|(m, d)| m + step * d).collect();
            let lo = model.sample(&ledger, &gate, &lo, 120.0).unwrap().unwrap();
            let hi = model.sample(&ledger, &gate, &hi, 120.0).unwrap().unwrap();
            let expected = (hi.value_k - lo.value_k) / (2.0 * step);
            let actual: f64 = gradient.iter().zip(direction).map(|(g, d)| g * d).sum();
            assert!((actual - expected).abs() < 5e-4 * expected.abs().max(1e-9),
                "qmc={qmc} direction={direction:?}: {actual:e} != {expected:e}");
        }
    }
}

#[test]
fn combined_calibration_resumes_from_retained_inputs_and_keeps_kink_refusals() {
    let f = Fixture::new();
    let text = controlled(&source(true));
    let (baseline, bp) = f.run(&text, "baseline-mixed", crate::exit::SUCCESS);
    let baseline = Fixture::report(&baseline, &bp);
    let path = f.inputs.join("split-mixed.fsim");
    std::fs::write(&path, &text).unwrap();
    let ledger = f.root.join("split-mixed.db");
    let prefix = crate::study::study_path(&path, &ledger, Some("1"), crate::OutputMode::Json);
    assert_eq!(prefix.exit_code, crate::exit::BUDGET, "{}", prefix.stderr);
    let prefix = J::parse(&prefix.stdout).unwrap();
    let report = Fixture::report(&prefix, &ledger);
    assert!(report.get("observations").unwrap().as_array().unwrap().is_empty());
    assert_eq!(report.get("mean_control").unwrap().f64_field("probe_solves_completed"), Some(1.0));
    let kink = text.replace("700000Pa", "1100000Pa").replace("1100000Pa)", "1300000Pa)");
    let (failed, fp) = f.run(&kink, "mixed-kink", crate::exit::REFUSED);
    let failed = Fixture::report(&failed, &fp);
    assert!(failed.str_field("failure").unwrap().contains("unequal-slope"));
    assert!(failed.get("observations").unwrap().as_array().unwrap().is_empty());
    std::fs::rename(&f.inputs, f.root.join("moved-mixed-inputs")).unwrap();
    let complete = crate::study::resume_path(prefix.str_field("run").unwrap(), &ledger, None, crate::OutputMode::Json);
    assert_eq!(complete.exit_code, crate::exit::SUCCESS, "{}", complete.stderr);
    let complete = J::parse(&complete.stdout).unwrap();
    let resumed = Fixture::report(&complete, &ledger);
    assert_eq!(resumed.get("observations"), baseline.get("observations"));
    assert_eq!(resumed.get("mean_control"), baseline.get("mean_control"));
}
