use super::*;
use super::super::super::tests::with_cx;

// Compacted so fixture edits are independent of the example's formatting.
static FIXTURE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| crate::json_read::compact(include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/repeated-contact-pulse.json"))));

fn controlled_text(controller: &str, repeat_prefix: &str) -> String {
    let old = "\"repeat\":{\"cycles\":10,\"max_total_steps\":10000}";
    assert!(FIXTURE.contains(old));
    FIXTURE.replace(old, &format!("\"repeat\":{{{repeat_prefix},\"max_total_steps\":10000,\"fan_controller\":{controller}}}"))
}

fn controller() -> &'static str {
    r#"{"sensor_vertex":4,"low_temperature_k":300.2,"high_temperature_k":301.0,
        "low_speed_multiplier":0.7,"high_speed_multiplier":1.3,"initial_speed_multiplier":0.7}"#
}

#[test]
fn thermostat_hysteresis_has_explicit_memory_and_thresholds() {
    let c = FanController::parse(&J::parse(controller()).unwrap()).unwrap();
    let mut state = FanControlState { multiplier: c.initial_speed_multiplier, switches: 0 };
    assert!(!c.update(300.5, &mut state).unwrap());
    assert_eq!(state.multiplier, 0.7);
    assert!(c.update(301.0, &mut state).unwrap());
    assert_eq!(state.multiplier, 1.3);
    assert!(!c.update(300.5, &mut state).unwrap());
    assert_eq!(state.multiplier, 1.3);
    assert!(c.update(300.2, &mut state).unwrap());
    assert_eq!(state.multiplier, 0.7);
    assert_eq!(state.switches, 2);
}

#[test]
fn repeated_run_changes_fan_state_from_the_simulated_sensor() {
    let text = controlled_text(controller(), "\"cycles\":4");
    let request = Request::parse(&text).unwrap();
    let doc = J::parse(&execute(&request, &CancelGate::new_clock_free()).unwrap()).unwrap();
    let repeated = doc.get("repeated_cycles").unwrap();
    let control = repeated.get("fan_controller").unwrap();
    assert_eq!(control.f64_field("sensor_vertex"), Some(4.0));
    assert!(control.f64_field("switches").unwrap() >= 1.0);
    let cycles = repeated.get("cycles").unwrap().as_array().unwrap();
    assert_eq!(cycles[0].f64_field("controller_speed_multiplier"), Some(0.7));
    assert!(cycles.iter().any(|cycle| cycle.f64_field("controller_speed_multiplier") == Some(1.3)));
    assert!(cycles.iter().all(|cycle| cycle.get("controller_sensor_start_k").unwrap().as_f64().unwrap().is_finite()));
}

#[test]
fn periodic_streak_can_be_reset_when_controller_state_would_change() {
    let periodic = Periodic { tolerance_k: 0.1, consecutive: 2 };
    let mut streak = 0;
    assert!(!periodic.observe(0.01, &mut streak));
    assert_eq!(streak, 1);
    streak = 0; // the repeat driver does this when the controller would switch
    assert_eq!(streak, 0);
    assert!(!periodic.observe(0.01, &mut streak));
    assert!(periodic.observe(0.01, &mut streak));
}

#[test]
fn malformed_sensor_and_out_of_domain_control_refuse_before_a_cycle() {
    for bad_controller in [
        r#"{"sensor_vertex":4,"low_temperature_k":301,"high_temperature_k":300,
            "low_speed_multiplier":0.7,"high_speed_multiplier":1.3,"initial_speed_multiplier":0.7}"#,
        r#"{"sensor_vertex":4,"low_temperature_k":300,"high_temperature_k":301,
            "low_speed_multiplier":1.3,"high_speed_multiplier":0.7,"initial_speed_multiplier":0.7}"#,
        r#"{"sensor_vertex":4,"low_temperature_k":300,"high_temperature_k":301,
            "low_speed_multiplier":0.7,"high_speed_multiplier":1.3,"initial_speed_multiplier":1.5}"#,
    ] {
        assert!(Config::parse(&J::parse(&format!(
            "{{\"cycles\":2,\"max_total_steps\":1000,\"fan_controller\":{bad_controller}}}"
        )).unwrap(), 75, false).is_err());
    }

    let bad_vertex = r#"{"sensor_vertex":99999,"low_temperature_k":300.2,"high_temperature_k":301,
        "low_speed_multiplier":0.7,"high_speed_multiplier":1.3,"initial_speed_multiplier":0.7}"#;
    let request = Request::parse(&controlled_text(bad_vertex, "\"cycles\":2")).unwrap();
    assert!(execute(&request, &CancelGate::new_clock_free()).is_err());

    let bad_speed = r#"{"sensor_vertex":4,"low_temperature_k":300.2,"high_temperature_k":301,
        "low_speed_multiplier":0.7,"high_speed_multiplier":2.0,"initial_speed_multiplier":0.7}"#;
    let request = Request::parse(&controlled_text(bad_speed, "\"cycles\":2")).unwrap();
    with_cx(|cx| {
        let schedule = request.transient.as_ref().unwrap();
        assert!(simulate(&request, cx, schedule, 1.0, schedule.repeat.unwrap()).is_err());
    });
}
