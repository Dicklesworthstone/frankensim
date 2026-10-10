//! Real repeated-cycle endpoints, with an independent latent-storage oracle.
//! Reuses the parent's JSON/process helpers and the actual CLI binary.

use super::*;

fn repeated(cycles: usize) -> J {
    let mut request = J::parse(FIXTURE).unwrap();
    put(
        member(&mut request, "transient"),
        "repeat",
        J::parse(&format!(
            r#"{{"cycles":{cycles},"max_total_steps":{}}}"#,
            4 * cycles
        ))
        .unwrap(),
    );
    request
}

fn periodic() -> J {
    let mut request = J::parse(FIXTURE).unwrap();
    put(
        member(&mut request, "transient"),
        "repeat",
        J::parse(
            r#"{"until_periodic":{"max_cycles":3,"temperature_tolerance_k":1e-9,"specific_enthalpy_tolerance_j_kg":1e-6,"consecutive_cycles":2},"max_total_steps":12}"#,
        )
        .unwrap(),
    );
    request
}

fn phase(result: &J) -> &J {
    result.get("repeated_cycles").unwrap()
}

fn expected_h(cycles: usize) -> [f64; 4] {
    // Each complete cycle supplies 240 J/kg and removes the same convection
    // and radiation loads while every vertex remains on the 350 K plateau.
    enthalpies(0.8, Some(300.0)).map(|h| 2000.0 + cycles as f64 * (h - 2000.0))
}

#[test]
fn fixed_cycles_carry_full_enthalpy_match_restarts_and_conserve_cumulative_energy() {
    let request = repeated(3);
    let result = run(&request);
    let report = phase(&result);
    assert_eq!(report.str_field("status"), Some("fixed-count-complete"));
    assert_eq!(report.get("periodic"), Some(&J::Null));
    assert_eq!(n(report, "cycles_completed"), 3.0);
    assert_eq!(n(report, "total_accepted_steps"), 12.0);
    near(n(report, "cycle_duration_s"), 0.8, 1e-14);
    near(n(report, "elapsed_time_s"), 2.4, 1e-14);
    near(n(report, "last_cycle_start_time_s"), 1.6, 1e-14);
    near(n(report, "sampled_peak_objective_k"), 350.0, 1e-9);
    near(n(report, "sampled_peak_time_s"), 0.0, 0.0);

    let (_, air, radiation) = heat(Some(300.0));
    let cycle_storage = 400.0 - 0.8 * (air + radiation);
    let drift = expected_h(1)
        .iter()
        .map(|h| (h - 2000.0).abs())
        .fold(0.0_f64, f64::max);
    assert!(drift > 100.0);
    for (key, expected) in [
        ("input_energy_j", 1200.0),
        ("stored_energy_change_j", 3.0 * cycle_storage),
        ("air_energy_gain_j", 2.4 * air),
        ("radiative_energy_loss_j", 2.4 * radiation),
    ] {
        near(n(report, key), expected, 3e-6);
    }
    near(n(report, "energy_residual_j"), 0.0, 2.4e-7);
    let rows = report.get("cycles").unwrap().as_array().unwrap();
    assert_eq!(rows.len(), 3);
    for (index, row) in rows.iter().enumerate() {
        assert_eq!(n(row, "cycle"), (index + 1) as f64);
        near(n(row, "start_time_s"), index as f64 * 0.8, 1e-14);
        near(n(row, "end_time_s"), (index + 1) as f64 * 0.8, 1e-14);
        near(n(row, "start_to_end_field_residual_k"), 0.0, 1e-9);
        near(
            n(row, "start_to_end_specific_enthalpy_residual_j_kg"),
            drift,
            2e-6,
        );
        near(n(row, "stored_energy_change_j"), cycle_storage, 1e-6);
        near(n(row, "input_energy_j"), 400.0, 1e-8);
        near(n(row, "air_energy_gain_j"), 0.8 * air, 1e-7);
        near(n(row, "radiative_energy_loss_j"), 0.8 * radiation, 1e-6);
        assert_eq!(n(row, "accepted_steps"), 4.0);
    }
    let final_h = values(&result, "solid_specific_enthalpies_j_kg");
    for (actual, expected) in final_h.iter().zip(expected_h(3)) {
        assert!((1000.0..3000.0).contains(&expected));
        near(*actual, expected, 4e-6);
    }
    for temperature in values(&result, "solid_temperatures_k") {
        near(temperature, 350.0, 1e-9);
    }
    // Total storage is referenced to the original cold-start h, while the
    // ordinary transient report is just the final cycle in local time.
    near(
        n(report, "stored_energy_change_j"),
        RHO / 24.0 * final_h.iter().map(|h| h - 2000.0).sum::<f64>(),
        1e-8,
    );
    let last = result.get("transient").unwrap();
    near(n(last, "time_s"), 0.8, 1e-14);
    assert_eq!(n(last, "steps"), 4.0);
    near(n(last, "stored_energy_change_j"), cycle_storage, 1e-6);
    near(
        n(last.get("enthalpy").unwrap(), "initial_total_enthalpy_j"),
        RHO / 24.0 * expected_h(2).iter().sum::<f64>(),
        3e-6,
    );
    let history = last.get("history").unwrap().as_array().unwrap();
    assert_eq!(history.len(), 5);
    near(n(&history[0], "time_s"), 0.0, 0.0);
    near(
        n(&history[0], "minimum_specific_enthalpy_j_kg"),
        expected_h(2).into_iter().fold(f64::INFINITY, f64::min),
        3e-6,
    );

    // Independent invocations restart from the actual accepted nodal h,
    // which is not recoverable from their identical plateau temperatures.
    let mut single = J::parse(FIXTURE).unwrap();
    let mut energies = [0.0; 4];
    let energy_keys = [
        "input_energy_j",
        "stored_energy_change_j",
        "air_energy_gain_j",
        "radiative_energy_loss_j",
    ];
    for cycle in 0..3 {
        let accepted = run(&single);
        for (sum, key) in energies.iter_mut().zip(energy_keys) {
            *sum += n(accepted.get("transient").unwrap(), key);
        }
        if cycle == 2 {
            for key in [
                "solid_specific_enthalpies_j_kg",
                "solid_temperatures_k",
                "solid_liquid_mass_fractions",
            ] {
                assert_eq!(result.get(key), accepted.get(key), "{key}");
            }
        }
        let chart = member(member(&mut single, "transient"), "enthalpy");
        remove(chart, "initial_specific_enthalpy_j_kg");
        put(
            chart,
            "initial_specific_enthalpies_j_kg",
            accepted
                .get("solid_specific_enthalpies_j_kg")
                .unwrap()
                .clone(),
        );
    }
    for (sum, key) in energies.into_iter().zip(energy_keys) {
        near(n(report, key), sum, 1e-8);
    }

    let mut unrolled = J::parse(FIXTURE).unwrap();
    let schedule = member(&mut unrolled, "transient");
    let J::Array(intervals) = member(schedule, "intervals") else {
        panic!()
    };
    *intervals = intervals.iter().cloned().cycle().take(6).collect();
    put(schedule, "max_steps", number(12.0));
    let explicit = run(&unrolled);
    for key in [
        "solid_specific_enthalpies_j_kg",
        "solid_temperatures_k",
        "solid_liquid_mass_fractions",
    ] {
        for (a, b) in values(&result, key).iter().zip(values(&explicit, key)) {
            near(*a, b, 2e-7);
        }
    }
    for key in energy_keys {
        near(
            n(report, key),
            n(explicit.get("transient").unwrap(), key),
            3e-6,
        );
    }
}

#[test]
fn constant_plateau_temperature_cannot_hide_nonperiodic_enthalpy_drift() {
    let request = periodic();
    let diagnostic = refuses(&request);
    assert!(diagnostic.contains("periodic"), "{diagnostic}");
    assert!(diagnostic.contains("exhaust"), "{diagnostic}");
    assert!(
        diagnostic.contains("enthalpy") || diagnostic.contains("J/kg"),
        "{diagnostic}"
    );
    let mut missing = request.clone();
    remove(
        member(
            member(member(&mut missing, "transient"), "repeat"),
            "until_periodic",
        ),
        "specific_enthalpy_tolerance_j_kg",
    );
    assert!(refuses(&missing).contains("specific_enthalpy_tolerance_j_kg"));
    let mut controlled = repeated(2);
    put(
        member(member(&mut controlled, "transient"), "repeat"),
        "fan_controller",
        J::parse(r#"{"sensor_vertex":0,"low_temperature_k":340,"high_temperature_k":360,"low_speed_multiplier":0.5,"high_speed_multiplier":1.5,"initial_speed_multiplier":1}"#).unwrap(),
    );
    let diagnostic = refuses(&controlled);
    assert!(
        diagnostic.contains("enthalpy") && diagnostic.contains("fan_controller"),
        "{diagnostic}"
    );
}

#[test]
fn periodic_stopping_requires_consecutive_full_state_passes_and_keeps_energy() {
    let mut request = periodic();
    put(member(&mut request, "solid"), "source_w_m3", number(0.0));
    let J::Array(boundaries) = member(member(&mut request, "hydraulics"), "boundaries") else {
        panic!()
    };
    put(&mut boundaries[0], "temperature_k", number(350.0));
    let J::Array(patches) = member(member(&mut request, "radiation"), "surfaces") else {
        panic!()
    };
    put(&mut patches[0], "ambient_temperature_k", number(350.0));
    let result = run(&request);
    let report = phase(&result);
    assert_eq!(
        report.str_field("status"),
        Some("periodic-field-tolerance-met")
    );
    assert_eq!(n(report, "cycles_completed"), 2.0);
    assert_eq!(n(report, "total_accepted_steps"), 8.0);
    near(n(report, "elapsed_time_s"), 1.6, 1e-14);
    let periodic = report.get("periodic").unwrap();
    near(n(periodic, "full_field_residual_k"), 0.0, 1e-9);
    near(
        n(periodic, "full_specific_enthalpy_residual_j_kg"),
        0.0,
        1e-6,
    );
    near(n(periodic, "specific_enthalpy_tolerance_j_kg"), 1e-6, 0.0);
    assert_eq!(n(periodic, "consecutive_cycles_met"), 2.0);
    assert_eq!(n(periodic, "required_consecutive_cycles"), 2.0);
    for h in values(&result, "solid_specific_enthalpies_j_kg") {
        near(h, 2000.0, 1e-6);
    }
    for t in values(&result, "solid_temperatures_k") {
        near(t, 350.0, 1e-9);
    }
    for key in [
        "input_energy_j",
        "stored_energy_change_j",
        "air_energy_gain_j",
        "radiative_energy_loss_j",
        "energy_residual_j",
    ] {
        near(n(report, key), 0.0, 1e-7);
    }
}
