use super::*;
const MODEL: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/equilibrium-uncertainty/joint-reliability.model"));
const DESIGN: &[u8] = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/equilibrium-uncertainty/joint-reliability.fit"));
fn loaded() -> EquilibriumDesignFile { EquilibriumDesignFile::from_bytes(MODEL, DESIGN, &CancelGate::new()).unwrap() }
fn bands() -> Vec<(String, f64)> { vec![("settled-band".into(), 1.0/256.0)] }
fn rows(force: f64) -> Vec<ForwardConstraintResult> {
    let q = if force <= 1.0 { force/64.0 } else { (force+2.0)/192.0 };
    let penetration = (q-1.0/64.0).max(0.0);
    [128.0*penetration, q, penetration, -64.0*q, q].into_iter()
        .map(|value| ForwardConstraintResult { value, residual: 0.0 }).collect()
}

#[test]
fn joint_event_uses_same_draw_dependence_not_a_product_or_rounded_zero_residual() {
    let loaded = loaded();
    let mut event = Event::new(loaded.problem(), &bands()).unwrap();
    // Deliberately zero residual placeholders: event classification must use
    // physical values, not assume that a normalized underflow means feasible.
    let values: Vec<_> = [0.625, 0.875, 1.125, 1.375].map(|f| event.record(&rows(f)).unwrap()).into();
    assert_eq!(values, [1.0, 0.0, 0.0, 1.0]);
    assert_eq!(event.failures, [1, 1, 1, 0, 1]);
    assert_eq!(event.samples, 4);
    let joint = 1.0-values.iter().sum::<f64>()/4.0;
    let product = (1.0-event.failures[0] as f64/4.0)*(1.0-event.failures[1] as f64/4.0);
    assert_eq!(joint, 0.5);
    assert_eq!(product, 0.5625);
}

#[test]
fn equality_bands_are_explicit_named_physical_units_and_never_borrow_solver_tolerances() {
    let loaded = loaded(); let p = loaded.problem();
    assert!(Event::new(p, &[]).is_err());
    for band in [vec![("unknown".into(), 0.1)], vec![("force-ceiling".into(), 0.1)],
        vec![("settled-band".into(), -1.0)], vec![("settled-band".into(), f64::INFINITY)],
        vec![("settled-band".into(), 0.1), ("settled-band".into(), 0.2)]] {
        assert!(Event::new(p, &band).is_err());
    }
    let mut exact = Event::new(p, &[("settled-band".into(), 0.0)]).unwrap();
    assert_eq!(exact.record(&rows(1.0)).unwrap(), 0.0);
    assert_eq!(exact.record(&rows(1.001)).unwrap(), 1.0);
    let mut band = Event::new(p, &bands()).unwrap();
    assert_eq!(band.record(&rows(1.001)).unwrap(), 0.0);
    let bare = std::str::from_utf8(DESIGN).unwrap().split("constraint_limits").next().unwrap();
    let no_constraints = EquilibriumDesignFile::from_bytes(MODEL, bare.as_bytes(), &CancelGate::new()).unwrap();
    assert!(Event::new(no_constraints.problem(), &[]).is_err());
}

#[test]
fn bad_response_family_does_not_change_any_diagnostic_counts_or_ranges() {
    let loaded = loaded(); let mut event = Event::new(loaded.problem(), &bands()).unwrap();
    event.record(&rows(1.0)).unwrap();
    let ranges = event.ranges.clone(); let failures = event.failures.clone();
    assert!(event.record(&[]).is_err());
    let mut malformed = rows(1.0); malformed[4].value = f64::NAN;
    assert!(event.record(&malformed).is_err());
    assert_eq!(event.samples, 1); assert_eq!(event.ranges, ranges); assert_eq!(event.failures, failures);
}

#[test]
fn joint_flags_are_exclusive_and_malformed_bands_refuse_before_file_access() {
    let base = "missing.model missing.fit --method mc --samples 32 --seed 73 --independent --fixed-x force-N 0 --all-constraints";
    let parse = |tail: &str| options(&format!("{base} {tail}").split_whitespace().map(str::to_owned).collect::<Vec<_>>());
    assert!(parse("--equality-tolerance settled-band 0.01").unwrap().all_constraints);
    for tail in ["--case load-a", "--target 0", "--limit-m 1", "--all-constraints",
        "--equality-tolerance a NaN", "--equality-tolerance a -0.1", "--equality-tolerance a 0 --equality-tolerance a 0"] {
        assert!(parse(tail).is_err(), "{tail}");
    }
    let scalar = "m d --method mc --samples 32 --seed 73 --independent --fixed-x force-N 0 --case load-a --target 0 --limit-m 1 --equality-tolerance settled-band 0";
    assert!(options(&scalar.split_whitespace().map(str::to_owned).collect::<Vec<_>>()).is_err());
}
