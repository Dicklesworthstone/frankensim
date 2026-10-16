//! G3 complete-history comparisons: physical h carries through cycle boundaries,
//! interval controls apply in every cycle, and design must retain warm-up peaks.

use super::*;

fn repeated_request(qoi: Option<&str>) -> J {
    let mut result = fan::fan_request(true, qoi);
    let schedule = member(&mut result, "transient");
    // Binary-exact clocks isolate cycle history from differences in rounding
    // decimal dt values. The entire first cycle remains inside latent storage.
    put(schedule, "max_step_s", number(0.25));
    put(
        schedule,
        "intervals",
        parse(
            r#"[
          {"duration_s":0.5,"power_scale":1,"fan_speed_ratio":1},
          {"duration_s":0.5,"power_scale":2,"fan_speed_ratio":1.2}
        ]"#,
        ),
    );
    put(
        schedule,
        "repeat",
        parse(r#"{"cycles":3,"max_total_steps":12}"#),
    );
    result
}

fn flattened(request: &J) -> J {
    let mut result = request.clone();
    let schedule = member(&mut result, "transient");
    remove(schedule, "repeat");
    let intervals = array_mut(member(schedule, "intervals"));
    let cycle = intervals.clone();
    intervals.extend(cycle.iter().cloned());
    intervals.extend(cycle);
    put(schedule, "max_steps", number(12.0));
    result
}

fn repeated(result: &J) -> &J {
    result.get("repeated_cycles").unwrap()
}

#[test]
fn repeated_latent_history_matches_flattening_and_shared_control_derivatives() {
    let base = repeated_request(None);
    let plain = run(&base);
    let differentiated = run(&repeated_request(Some("final")));
    let flat = run(&flattened(&repeated_request(Some("final"))));
    fan::same_physical_history(&plain, &differentiated);
    assert_eq!(
        repeated(&plain).get("cycles"),
        repeated(&differentiated).get("cycles")
    );
    for key in ["solid_specific_enthalpies_j_kg", "solid_temperatures_k"] {
        assert_eq!(plain.get(key), flat.get(key));
    }
    let history = at(&flat, &["transient", "history"]).as_array().unwrap();
    assert_eq!(history.len(), 13);
    for row in &history[1..=4] {
        assert_eq!(n(row, "objective_temperature_k"), 350.0);
        assert!((1000.0..3000.0).contains(&n(row, "minimum_specific_enthalpy_j_kg")));
        assert!((1000.0..3000.0).contains(&n(row, "maximum_specific_enthalpy_j_kg")));
    }
    assert!(final_temperature(&plain) > 400.0);
    let total = repeated(&differentiated);
    assert_eq!(n(total, "cycles_completed"), 3.0);
    assert_eq!(n(total, "total_accepted_steps"), 12.0);
    assert_eq!(at(&differentiated, &["transient", "adjoint"]), &J::Null);
    let adjoint = total.get("adjoint").unwrap();
    let flat_adjoint = at(&flat, &["transient", "adjoint"]);
    assert_eq!(n(adjoint, "cycles"), 3.0);
    assert_eq!(n(adjoint, "state_index"), 12.0);
    assert_eq!(n(adjoint, "time_s"), 3.0);
    close(n(adjoint, "value_k"), final_temperature(&plain));
    let rows = adjoint.get("intervals").unwrap().as_array().unwrap();
    let flat_rows = flat_adjoint.get("intervals").unwrap().as_array().unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!(flat_rows.len(), 6);
    for control in [
        "dtemperature_dpower_multiplier_k",
        "dtemperature_dlog_fan_speed_ratio_k",
    ] {
        for (interval, row) in rows.iter().enumerate() {
            close(
                n(row, control),
                (0..3)
                    .map(|cycle| n(&flat_rows[2 * cycle + interval], control))
                    .sum(),
            );
        }
    }
    assert!(n(&flat_rows[0], "dtemperature_dpower_multiplier_k") > 0.1);
    assert!(
        n(&flat_rows[0], "dtemperature_dlog_fan_speed_ratio_k") < -1e-3,
        "an entirely latent first cycle must influence the final sensible field"
    );
    let key = "dtemperature_dinitial_specific_enthalpies_k_kg_j";
    for (actual, expected) in vector(adjoint, key).iter().zip(vector(flat_adjoint, key)) {
        close(*actual, expected);
    }
    let initial = n(
        adjoint,
        "dtemperature_duniform_initial_specific_enthalpy_k_kg_j",
    );
    assert!(initial > 0.01);
    let mut plus = base.clone();
    let mut minus = base.clone();
    shift_initial(&mut plus, None, 0.1);
    shift_initial(&mut minus, None, -0.1);
    close(
        initial,
        (final_temperature(&run(&plus)) - final_temperature(&run(&minus))) / 0.2,
    );
    let epsilon = 1e-3_f64;
    for (interval, row) in rows.iter().enumerate() {
        for fan_control in [false, true] {
            let mut plus = base.clone();
            let mut minus = base.clone();
            let key = if fan_control {
                fan::scale_speed(&mut plus, Some(interval), epsilon.exp());
                fan::scale_speed(&mut minus, Some(interval), (-epsilon).exp());
                "dtemperature_dlog_fan_speed_ratio_k"
            } else {
                scale_interval(&mut plus, interval, 1.0 + epsilon);
                scale_interval(&mut minus, interval, 1.0 - epsilon);
                "dtemperature_dpower_multiplier_k"
            };
            close(
                n(row, key),
                (final_temperature(&run(&plus)) - final_temperature(&run(&minus)))
                    / (2.0 * epsilon),
            );
        }
    }
}

fn warmup(multiplier: f64, adjoint: bool) -> J {
    let mut result = repeated_request(adjoint.then_some("sampled-peak"));
    put(
        &mut result,
        "objective",
        parse(r#"{"mean_wall_region":"first","gradient":false}"#),
    );
    let schedule = member(&mut result, "transient");
    put(
        member(schedule, "enthalpy"),
        "initial_specific_enthalpies_j_kg",
        parse("[2000,9000,1970,2010]"),
    );
    put(member(schedule, "repeat"), "max_total_steps", number(24.0));
    let rows = array_mut(member(schedule, "intervals"));
    for (row, power) in rows.iter_mut().zip([0.2, 0.4]) {
        put(row, "duration_s", number(1.0));
        put(row, "power_scale", number(power * multiplier));
    }
    result
}

#[test]
fn repeated_power_sizing_keeps_an_earlier_warmup_peak_and_its_adjoint() {
    // The initially hot vertex is outside the observed first patch. Its stored
    // energy heats that patch during cycle 1, then leaves through air/radiation.
    // Define the target with a separate full run, never an assumed steady limit.
    let target = run(&warmup(0.85, false));
    let limit = n(repeated(&target), "sampled_peak_objective_k");
    let mut request = warmup(1.0, true);
    let schedule = member(&mut request, "transient");
    put(schedule, "temperature_limit_k", number(limit));
    put(
        schedule,
        "power_design",
        parse(
            r#"{
          "min_power_multiplier":0.5,"max_power_multiplier":1.4,
          "power_multiplier_tolerance":0.0001,"temperature_tolerance_k":0.0001,
          "max_evaluations":48
        }"#,
        ),
    );
    let result = run(&request);
    let total = repeated(&result);
    let cycles = total.get("cycles").unwrap().as_array().unwrap();
    assert!(
        n(&cycles[0], "sampled_peak_objective_k") > n(&cycles[2], "sampled_peak_objective_k") + 1.0
    );
    let peak = n(total, "sampled_peak_objective_k");
    assert_eq!(peak, n(&cycles[0], "sampled_peak_objective_k"));
    assert!(peak <= limit && limit - peak <= 0.0001);
    let design = result.get("transient_power_design").unwrap();
    assert!(n(design, "newton_trials") > 0.0);
    assert!(n(design, "multiplier_bracket_width") <= 0.0001);
    let failed = design.get("failed_upper").unwrap();
    assert!(n(failed, "sampled_peak_objective_k") > limit);
    let selected = n(design, "selected_power_multiplier");
    assert!((selected - 0.85).abs() < 0.001);
    let replay = run(&warmup(selected, false));
    fan::same_physical_history(&result, &replay);
    assert_eq!(total.get("cycles"), repeated(&replay).get("cycles"));
    let adjoint = total.get("adjoint").unwrap();
    assert_eq!(at(&result, &["transient", "adjoint"]), &J::Null);
    assert!(n(adjoint, "state_index") > 0.0 && n(adjoint, "state_index") <= 8.0);
    assert!(n(adjoint, "state_index") < n(total, "total_accepted_steps"));
    assert!(n(adjoint, "time_s") <= 2.0);
    close(n(adjoint, "value_k"), peak);
    let slope = adjoint
        .get("intervals")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|row| n(row, "dtemperature_dpower_multiplier_k"))
        .sum::<f64>()
        / selected;
    let trials = design.get("history").unwrap().as_array().unwrap();
    let accepted = trials
        .iter()
        .find(|row| n(row, "power_multiplier") == selected)
        .unwrap();
    close(n(accepted, "dpeak_dmultiplier_k"), slope);
    let epsilon = 1e-4;
    let plus = run(&warmup(selected + epsilon, false));
    let minus = run(&warmup(selected - epsilon, false));
    close(
        slope,
        (n(repeated(&plus), "sampled_peak_objective_k")
            - n(repeated(&minus), "sampled_peak_objective_k"))
            / (2.0 * epsilon),
    );
}
