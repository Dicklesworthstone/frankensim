//! G1/G3: the shared search sizes a real discrete pulse, not the cooled endpoint.
use super::*;
use fs_project::uncertainty::Target;

fn policy() -> DesignPolicy {
    DesignPolicy { control: DesignControl::PowerScheduleScale { tolerance: 0.001 },
        temperature_tolerance_k: 0.001 }
}
fn control() -> UniformParameter {
    UniformParameter { name: "pulse".into(), entity: "chip".into(),
        target: Target::PowerScheduleScale, low: 0.0, high: 10.0 }
}
fn trial(amplitude: f64, peak: f64, final_k: f64) -> Trial {
    Trial { sample: Sample { run: "physical-reference".into(), project_hash: String::new(),
        qoi_receipt: hash_bytes(b"qoi"), value_k: final_k, parameters: vec![amplitude] },
        temperature_k: peak, time_s: Some(0.25), state: Some(2), conduction: hash_bytes(b"conduction") }
}
fn pulse(amplitude: f64, initial: f64) -> Trial {
    let mut t = initial;
    let mut peak = initial;
    let mut start = 0.0;
    for end in [0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0] {
        let dt = end - start;
        let watts = if end <= 0.25 { 20.0 * amplitude } else { 0.0 };
        // Independent backward-Euler balance of C=10 J/K, hA=2 W/K.
        t = (10.0 * t + dt * (watts + 2.0 * 300.0)) / (10.0 + 2.0 * dt);
        peak = peak.max(t);
        start = end;
    }
    trial(amplitude, peak, t)
}

#[test]
fn g1_pulse_search_resolves_the_independent_peak_limit_not_the_final_limit() {
    let mut rows = Vec::new();
    for _ in 0..40 {
        match next(&control(), policy(), 303.0, &rows) {
            Next::Evaluate(amplitude) => rows.push(pulse(amplitude, 300.0)),
            Next::Complete("target-bracketed") => break,
            other => panic!("unexpected pulse search {other:?}"),
        }
    }
    assert_eq!(next(&control(), policy(), 303.0, &rows), Next::Complete("target-bracketed"));
    let exact = 3.0 / (10.0 * (1.0 - (10.0_f64 / 10.25).powi(2)));
    let chosen = selected(&rows, policy(), 303.0).unwrap();
    assert!(chosen.coordinate() <= exact + 1e-10);
    assert!((chosen.coordinate() - exact).abs() <= 0.001);
    assert!(chosen.temperature_k <= 303.0);
    assert!(chosen.sample.value_k < 302.5, "the cooled endpoint would admit an excessive pulse");
    assert_eq!(next(&control(), policy(), 303.0, &[pulse(0.0, 310.0)]),
        Next::Refuse("minimum-schedule-amplitude-infeasible"));
}

#[test]
fn g3_partial_operating_cases_and_repeated_continuation_do_not_relax_feasibility() {
    let mut rows = Vec::new();
    let mut resumed = Vec::new();
    for _ in 0..80 {
        let completed = cases::completed(&rows, 2);
        let decision = next(&control(), policy(), 303.0, &completed);
        let restored = cases::completed(&resumed, 2);
        assert_eq!(next(&control(), policy(), 303.0, &restored), decision);
        match decision {
            Next::Evaluate(a) => {
                let odd = rows.len() % 2 == 1;
                let t = if odd { 301.0 + 0.5 * a } else { 300.0 + 0.8 * a };
                rows.push(trial(a, t, 300.0));
                resumed = rows.clone();
                assert_eq!(cases::completed(&rows, 2).len(), rows.len() / 2);
            }
            Next::Complete("target-bracketed") => break,
            other => panic!("unexpected case search {other:?}"),
        }
    }
    let full = cases::completed(&rows, 2);
    assert_eq!(next(&control(), policy(), 303.0, &full), Next::Complete("target-bracketed"));
    let chosen = selected(&full, policy(), 303.0).unwrap();
    assert!((3.749..=3.75).contains(&chosen.coordinate()));
    assert!(300.0 + 0.8 * chosen.coordinate() <= 303.0);
    assert!(301.0 + 0.5 * chosen.coordinate() <= 303.0);
}

#[test]
fn g0_schedule_identity_and_units_cannot_alias_static_watts_or_fan_speed() {
    let row = J::parse(&trial(1.5, 303.0, 301.0).json(policy())).unwrap();
    assert_eq!(row.f64_field("power_schedule_scale"), Some(1.5));
    assert!(row.get("power_w").is_none() && row.get("fan_speed_ratio").is_none());
    assert_eq!(driver(policy()), SCHEDULE_DRIVER);
    for name in [SCHEDULE_DRIVER, SCHEDULE_CASE_DRIVER] {
        assert!(is_driver(Some(name)));
        assert_ne!(observations_kind(name), OBSERVATIONS);
        assert_ne!(observations_kind(name), FAN_OBSERVATIONS);
        assert_ne!(ir(name, hash_bytes(b"model"), 1), ir(DRIVER, hash_bytes(b"model"), 1));
    }
    let metadata = coordinate::presentation(policy());
    assert_eq!(metadata.bracket_width, "bracket_width_scale");
    assert_eq!(metadata.objective, "maximize-regional-power-schedule-scale");
    assert!(!metadata.unit.contains('W'));
}
