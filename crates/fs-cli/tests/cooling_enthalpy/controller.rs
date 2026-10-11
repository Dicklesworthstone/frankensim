//! Actual-binary cycle-rate fan feedback through sensible and latent storage.
//! Manual replays carry accepted nodal h, never an inverse of temperature.

use super::*;

fn fan_request(initial_h: f64) -> J {
    let mut request = J::parse(FIXTURE).unwrap();
    let fan_example = J::parse(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/enthalpy-fan-sizing.json"
    )))
    .unwrap();
    put(
        &mut request,
        "hydraulics",
        fan_example.get("hydraulics").unwrap().clone(),
    );
    remove(&mut request, "radiation");
    let schedule = member(&mut request, "transient");
    put(
        member(schedule, "enthalpy"),
        "initial_specific_enthalpy_j_kg",
        number(initial_h),
    );
    let J::Array(intervals) = member(schedule, "intervals") else {
        panic!()
    };
    for (interval, speed) in intervals.iter_mut().zip([1.0, 1.2]) {
        put(interval, "fan_speed_ratio", number(speed));
    }
    request
}

fn controller(low_k: f64, high_k: f64, initial: f64) -> J {
    J::parse(&format!(
        r#"{{"sensor_vertex":0,"low_temperature_k":{low_k},"high_temperature_k":{high_k},"low_speed_multiplier":0.75,"high_speed_multiplier":1.5,"initial_speed_multiplier":{initial}}}"#
    ))
    .unwrap()
}

fn repeated(request: &mut J, control: J, cycles: usize) {
    let mut policy = J::parse(&format!(r#"{{"cycles":{cycles},"max_total_steps":256}}"#)).unwrap();
    put(&mut policy, "fan_controller", control);
    put(member(request, "transient"), "repeat", policy);
}

fn periodic(request: &mut J, max_cycles: usize, tolerance_k: f64, tolerance_h: f64) {
    let policy = member(member(request, "transient"), "repeat");
    remove(policy, "cycles");
    put(
        policy,
        "until_periodic",
        J::parse(&format!(
            r#"{{"max_cycles":{max_cycles},"temperature_tolerance_k":{tolerance_k},"specific_enthalpy_tolerance_j_kg":{tolerance_h},"consecutive_cycles":2}}"#
        ))
        .unwrap(),
    );
}

fn scale_speeds(request: &mut J, multiplier: f64) {
    let J::Array(intervals) = member(member(request, "transient"), "intervals") else {
        panic!()
    };
    for interval in intervals {
        let speed = n(interval, "fan_speed_ratio");
        put(interval, "fan_speed_ratio", number(speed * multiplier));
    }
}

fn same_state(actual: &J, expected: &J) {
    for key in [
        "solid_specific_enthalpies_j_kg",
        "solid_temperatures_k",
        "solid_liquid_mass_fractions",
    ] {
        assert_eq!(actual.get(key), expected.get(key), "{key}");
    }
}

#[test]
fn cycle_start_hysteresis_matches_manual_fan_schedules_with_full_h_memory() {
    let base = fan_request(500.0); // 300 K on the sensible branch.
    let mut controlled = base.clone();
    repeated(&mut controlled, controller(310.0, 330.0, 1.25), 4);
    let result = run(&controlled);
    let report = result.get("repeated_cycles").unwrap();
    let rows = report.get("cycles").unwrap().as_array().unwrap();
    assert_eq!(rows.len(), 4);

    // The first heating cycle ends inside the hold band. The next cycle
    // crosses the high threshold, and the final cycles enter latent storage.
    assert!(n(&rows[1], "controller_sensor_start_k") > 310.0);
    assert!(n(&rows[1], "controller_sensor_start_k") < 330.0);
    assert!(n(&rows[2], "controller_sensor_start_k") > 330.0);
    let mut initial_h = J::Array(vec![number(500.0); 4]);
    let mut sensor = 300.0;
    let mut multiplier = 1.25;
    let mut switches = 0;
    let mut sums = [0.0; 3];
    let energy_keys = [
        "input_energy_j",
        "stored_energy_change_j",
        "air_energy_gain_j",
    ];
    for (index, row) in rows.iter().enumerate() {
        let next = if sensor >= 330.0 {
            1.5
        } else if sensor <= 310.0 {
            0.75
        } else {
            multiplier
        };
        let switched = next != multiplier;
        switches += usize::from(switched);
        multiplier = next;
        near(n(row, "controller_sensor_start_k"), sensor, 0.0);
        near(n(row, "controller_speed_multiplier"), multiplier, 0.0);
        assert_eq!(row.get("controller_switched"), Some(&J::Bool(switched)));
        assert_eq!(
            values(row, "applied_interval_fan_speed_ratios"),
            vec![multiplier, 1.2 * multiplier]
        );

        let mut manual = base.clone();
        scale_speeds(&mut manual, multiplier);
        let storage = member(member(&mut manual, "transient"), "enthalpy");
        remove(storage, "initial_specific_enthalpy_j_kg");
        put(storage, "initial_specific_enthalpies_j_kg", initial_h);
        let accepted = run(&manual);
        sensor = values(&accepted, "solid_temperatures_k")[0];
        near(n(row, "controller_sensor_end_k"), sensor, 0.0);
        for (sum, key) in sums.iter_mut().zip(energy_keys) {
            let energy = n(accepted.get("transient").unwrap(), key);
            near(n(row, key), energy, 1e-9);
            *sum += energy;
        }
        initial_h = accepted
            .get("solid_specific_enthalpies_j_kg")
            .unwrap()
            .clone();
        if index == rows.len() - 1 {
            same_state(&result, &accepted);
        }
    }
    assert_eq!(switches, 2);
    let control = report.get("fan_controller").unwrap();
    near(n(control, "switches"), switches as f64, 0.0);
    near(n(control, "final_speed_multiplier"), 1.5, 0.0);
    for (sum, key) in sums.into_iter().zip(energy_keys) {
        near(n(report, key), sum, 1e-9);
    }
    near(n(report, "energy_residual_j"), 0.0, 3.2e-7);
    for fraction in values(&result, "solid_liquid_mass_fractions") {
        assert!(fraction > 0.0 && fraction < 1.0);
    }
}

#[test]
fn hold_band_unit_multiplier_preserves_adaptive_repetition_exactly() {
    let mut controlled = fan_request(2000.0);
    repeated(&mut controlled, controller(340.0, 360.0, 1.0), 3);
    let schedule = member(&mut controlled, "transient");
    put(schedule, "max_steps", number(64.0));
    put(
        schedule,
        "adaptive",
        J::parse(r#"{"absolute_tolerance_k":1e-7,"absolute_specific_enthalpy_tolerance_j_kg":0.01,"relative_tolerance":0,"minimum_trial_step_s":0.001,"max_trials":128}"#).unwrap(),
    );
    let mut plain = controlled.clone();
    remove(
        member(member(&mut plain, "transient"), "repeat"),
        "fan_controller",
    );
    let actual = run(&controlled);
    let expected = run(&plain);
    same_state(&actual, &expected);
    assert_eq!(actual.get("transient"), expected.get("transient"));
    let report = actual.get("repeated_cycles").unwrap();
    let baseline = expected.get("repeated_cycles").unwrap();
    for key in [
        "total_accepted_steps",
        "total_solid_solves",
        "input_energy_j",
        "stored_energy_change_j",
        "air_energy_gain_j",
        "energy_residual_j",
    ] {
        assert_eq!(report.get(key), baseline.get(key), "{key}");
    }
    near(n(report, "total_accepted_steps"), 24.0, 0.0);
    near(
        n(report.get("fan_controller").unwrap(), "switches"),
        0.0,
        0.0,
    );
    for row in report.get("cycles").unwrap().as_array().unwrap() {
        near(n(row, "controller_speed_multiplier"), 1.0, 0.0);
        assert_eq!(row.get("controller_state_stable"), Some(&J::Bool(true)));
    }
}

#[test]
fn periodic_control_requires_h_temperature_and_stable_next_fan_state() {
    let mut equilibrium = fan_request(2000.0);
    repeated(&mut equilibrium, controller(340.0, 360.0, 1.25), 3);
    periodic(&mut equilibrium, 3, 1e-9, 1e-6);
    let drifting = equilibrium.clone();
    assert!(refuses(&drifting).contains("periodic"));
    put(
        member(&mut equilibrium, "solid"),
        "source_w_m3",
        number(0.0),
    );
    put(
        member(member(&mut equilibrium, "hydraulics"), "fan"),
        "temperature_k",
        number(350.0),
    );
    let result = run(&equilibrium);
    let report = result.get("repeated_cycles").unwrap();
    near(n(report, "cycles_completed"), 2.0, 0.0);
    let convergence = report.get("periodic").unwrap();
    near(n(convergence, "full_field_residual_k"), 0.0, 1e-9);
    near(
        n(convergence, "full_specific_enthalpy_residual_j_kg"),
        0.0,
        1e-6,
    );
    assert_eq!(
        convergence.get("controller_state_stable"),
        Some(&J::Bool(true))
    );

    // Deliberately loose physical tolerances isolate controller stability:
    // cycle one passes them but its end sensor demands a switch next cycle.
    let mut switching = fan_request(500.0);
    repeated(&mut switching, controller(305.0, 310.0, 0.75), 3);
    periodic(&mut switching, 2, 1000.0, 1_000_000.0);
    assert!(refuses(&switching).contains("periodic"));
    periodic(&mut switching, 3, 1000.0, 1_000_000.0);
    let result = run(&switching);
    let report = result.get("repeated_cycles").unwrap();
    near(n(report, "cycles_completed"), 3.0, 0.0);
    let rows = report.get("cycles").unwrap().as_array().unwrap();
    assert_eq!(
        rows[0].get("controller_state_stable"),
        Some(&J::Bool(false))
    );
    for row in &rows[1..] {
        assert_eq!(row.get("controller_state_stable"), Some(&J::Bool(true)));
    }
    near(
        n(report.get("periodic").unwrap(), "consecutive_cycles_met"),
        2.0,
        0.0,
    );
}

#[test]
fn controller_adjoint_and_invalid_sensor_or_unused_speed_refuse_without_output() {
    let mut request = fan_request(2000.0);
    repeated(&mut request, controller(340.0, 360.0, 1.0), 2);
    let mut differentiated = request.clone();
    put(
        member(&mut differentiated, "transient"),
        "adjoint",
        J::parse(r#"{"qoi":"final","max_checkpoint_bytes":1048576}"#).unwrap(),
    );
    let diagnostic = refuses(&differentiated);
    assert!(diagnostic.contains("adjoint") && diagnostic.contains("controller"));
    for (key, value) in [("sensor_vertex", 4.0), ("high_speed_multiplier", 2.0)] {
        let mut invalid = request.clone();
        put(
            member(
                member(member(&mut invalid, "transient"), "repeat"),
                "fan_controller",
            ),
            key,
            number(value),
        );
        let diagnostic = refuses(&invalid);
        assert!(
            diagnostic.contains(if key == "sensor_vertex" {
                "sensor"
            } else {
                "speed"
            }),
            "{diagnostic}"
        );
    }
}
