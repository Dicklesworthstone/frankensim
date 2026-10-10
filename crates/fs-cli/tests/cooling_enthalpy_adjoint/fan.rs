//! G3 trajectory checks for single-bank affinity controls through latent storage.
//! Reuse the real-binary enthalpy fixture and JSON helpers from the parent test.

use super::*;

fn fan_request(correlated: bool, qoi: Option<&str>) -> J {
    let mut result = request(qoi);
    // The unit-speed intersection is Q=7, pressure=2: the first and bypass
    // paths carry 5 and 2, then mix before the final exchanger. The synthetic
    // fan obeys affinity scaling; these values are not measured hardware data.
    put(
        &mut result,
        "hydraulics",
        parse(
            r#"{
      "node_count":3,
      "fan":{
        "name":"phase-history fan","source":"Synthetic affinity fixture",
        "source_id":"enthalpy-fan-test-v1","inlet":0,"outlet":2,
        "temperature_k":330,"count":1,"arrangement":"series",
        "speed_ratio":1,"min_speed_ratio":0.4,"max_speed_ratio":2,
        "min_flow_m3_s":0.001,"pressure_tolerance_rel":0.05,
        "points":[{"flow_m3_s":0,"static_pressure_pa":4},
                  {"flow_m3_s":14,"static_pressure_pa":0}],
        "max_iterations":160,"flow_tolerance_m3_s":1e-11,
        "pressure_tolerance_pa":1e-10
      },
      "branches":[
        {"name":"first","from":0,"to":1,"resistance_pa_s2_m6":0.04,"source":"synthetic first path","regions":["first"]},
        {"name":"bypass","from":0,"to":1,"resistance_pa_s2_m6":0.25,"source":"synthetic bypass","regions":[]},
        {"name":"last","from":1,"to":2,"resistance_pa_s2_m6":0.02040816326530612,"source":"synthetic last path","regions":["last"]}
      ]}"#,
        ),
    );
    let surfaces = array_mut(member(member(&mut result, "solid"), "surfaces"));
    for (i, surface) in surfaces.iter_mut().enumerate() {
        if correlated {
            remove(surface, "htc_w_m2_k");
            // Re=1000 and Pr=1 at unit speed, L/D=5. The complete sizing
            // bracket keeps both ports strictly inside the Hausen domain.
            let (area, property) = if i == 0 {
                (0.00625, 0.08)
            } else {
                (0.0175, 0.04)
            };
            put(
                surface,
                "convection",
                parse(&format!(
                    r#"{{
                  "card":"convection.circular-duct-hausen-developing",
                  "hydraulic_diameter_m":0.1,"flow_area_m2":{area},
                  "channel_length_m":0.5,"dynamic_viscosity_pa_s":{property},
                  "thermal_conductivity_w_m_k":{property},
                  "source":"Synthetic fixed-property derivative fixture"
                }}"#
                )),
            );
        } else {
            put(
                surface,
                "htc_w_m2_k",
                number(if i == 0 { 8.0 } else { 4.0 }),
            );
        }
    }
    let intervals = array_mut(member(member(&mut result, "transient"), "intervals"));
    for (interval, speed) in intervals.iter_mut().zip([1.0, 1.2]) {
        put(interval, "fan_speed_ratio", number(speed));
    }
    result
}

fn scale_speed(request: &mut J, interval: Option<usize>, multiplier: f64) {
    let intervals = array_mut(member(member(request, "transient"), "intervals"));
    for (index, row) in intervals.iter_mut().enumerate() {
        if interval.is_none_or(|i| i == index) {
            let speed = n(row, "fan_speed_ratio") * multiplier;
            put(row, "fan_speed_ratio", number(speed));
        }
    }
}

fn scaled_fan(multiplier: f64) -> J {
    let mut result = fan_request(true, None);
    scale_speed(&mut result, None, multiplier);
    result
}

fn same_physical_history(first: &J, second: &J) {
    for path in [
        &["solid_specific_enthalpies_j_kg"][..],
        &["solid_temperatures_k"][..],
        &["transient", "history"][..],
    ] {
        assert_eq!(at(first, path), at(second, path));
    }
}

#[test]
fn interval_fan_gradients_cross_a_latent_plateau_with_air_and_convection_feedback() {
    for correlated in [false, true] {
        let base = fan_request(correlated, None);
        let plain = run(&base);
        let differentiated = run(&fan_request(correlated, Some("final")));
        same_physical_history(&plain, &differentiated);
        let trajectory = plain.get("transient").unwrap();
        let history = trajectory.get("history").unwrap().as_array().unwrap();
        for row in &history[1..=2] {
            assert_eq!(n(row, "objective_temperature_k"), 350.0);
            assert!((1000.0..3000.0).contains(&n(row, "minimum_specific_enthalpy_j_kg")));
            assert!((1000.0..3000.0).contains(&n(row, "maximum_specific_enthalpy_j_kg")));
        }
        assert!(final_temperature(&plain) > 360.0);
        assert!(n(trajectory, "radiative_energy_loss_j") > 0.0);
        close(n(plain.get("fan").unwrap(), "flow_m3_s"), 7.0 * 1.2);
        let branches = plain.get("branches").unwrap().as_array().unwrap();
        close(n(&branches[0], "flow_m3_s"), 5.0 * 1.2);
        close(n(&branches[1], "flow_m3_s"), 2.0 * 1.2);
        // The bypass really mixes with the warmed first branch; this is not
        // an independent set of fixed-temperature reservoir boundaries.
        close(
            n(&branches[2], "inlet_k"),
            (5.0 * n(&branches[0], "outlet_k") + 2.0 * 330.0) / 7.0,
        );
        let adjoint = at(&differentiated, &["transient", "adjoint"]);
        let rows = adjoint.get("intervals").unwrap().as_array().unwrap();
        assert!(
            n(&rows[0], "dtemperature_dlog_fan_speed_ratio_k") < -1e-3,
            "fan cooling during T'=0 must still affect the later sensible state"
        );
        let epsilon = 1e-3_f64;
        for (index, row) in rows.iter().enumerate() {
            let mut plus = base.clone();
            let mut minus = base.clone();
            scale_speed(&mut plus, Some(index), epsilon.exp());
            scale_speed(&mut minus, Some(index), (-epsilon).exp());
            close(
                n(row, "dtemperature_dlog_fan_speed_ratio_k"),
                (final_temperature(&run(&plus)) - final_temperature(&run(&minus)))
                    / (2.0 * epsilon),
            );
        }
    }
}

#[test]
fn fan_sizing_replays_evaluated_phase_histories_and_the_nonunit_chain_rule() {
    let limit = peak(&run(&scaled_fan(0.85)));
    assert!(limit > 360.0, "the target is past latent storage");
    for with_adjoint in [true, false] {
        let mut request = fan_request(true, with_adjoint.then_some("sampled-peak"));
        let schedule = member(&mut request, "transient");
        put(schedule, "temperature_limit_k", number(limit));
        put(
            schedule,
            "fan_speed_design",
            parse(
                r#"{
              "min_speed_multiplier":0.65,"max_speed_multiplier":1.4,
              "speed_multiplier_tolerance":0.0001,
              "temperature_tolerance_k":0.0001,"max_evaluations":48
            }"#,
            ),
        );
        let result = run(&request);
        let design = result.get("transient_fan_speed_design").unwrap();
        assert_eq!(
            design.str_field("search_method"),
            Some(if with_adjoint {
                "safeguarded-adjoint-newton-bisection"
            } else {
                "bisection"
            })
        );
        let selected = n(design, "selected_speed_multiplier");
        assert!((selected - 0.85).abs() < 0.001);
        assert!(peak(&result) <= limit);
        assert!(limit - peak(&result) <= 0.0001);
        assert!(n(design, "multiplier_bracket_width") <= 0.0001);
        let failed = design.get("failed_lower").unwrap();
        assert!(n(failed, "speed_multiplier") < selected);
        assert!(n(failed, "sampled_peak_objective_k") > limit);
        same_physical_history(&result, &run(&scaled_fan(selected)));
        let trajectory = result.get("transient").unwrap();
        let trials = design.get("history").unwrap().as_array().unwrap();
        if with_adjoint {
            assert!(n(design, "newton_trials") > 0.0);
            let adjoint = trajectory.get("adjoint").unwrap();
            close(n(adjoint, "value_k"), peak(&result));
            let relative: f64 = adjoint
                .get("intervals")
                .unwrap()
                .as_array()
                .unwrap()
                .iter()
                .map(|row| n(row, "dtemperature_dlog_fan_speed_ratio_k"))
                .sum();
            let slope = relative / selected;
            let accepted = trials
                .iter()
                .find(|row| n(row, "speed_multiplier") == selected)
                .unwrap();
            close(n(accepted, "dpeak_dmultiplier_k"), slope);
            let epsilon = 1e-4;
            close(
                slope,
                (peak(&run(&scaled_fan(selected + epsilon)))
                    - peak(&run(&scaled_fan(selected - epsilon))))
                    / (2.0 * epsilon),
            );
        } else {
            assert_eq!(n(design, "newton_trials"), 0.0);
            assert_eq!(trajectory.get("adjoint"), Some(&J::Null));
            assert_eq!(
                n(trajectory, "total_solid_solves"),
                n(trajectory, "forward_solid_solves")
            );
            assert!(
                trials
                    .iter()
                    .all(|row| row.get("dpeak_dmultiplier_k") == Some(&J::Null))
            );
        }
    }
}
