//! Actual-binary return mixing, accepted enthalpy energy and fresh-air controls.
//! The unit tetrahedron and phase chart are synthetic numerical references.

use super::*;

const RETURN_FRACTION: f64 = 0.6;

fn request(fraction: f64) -> J {
    let mut input = J::parse(FIXTURE).unwrap();
    remove(&mut input, "radiation");
    put(
        &mut input,
        "recirculation",
        J::parse(&format!(
            r#"{{"model":"prescribed-adiabatic-return",
          "source":"Synthetic unit-capacity return fixture",
          "temperature_tolerance_k":1e-9,
          "links":[{{"supply_node":0,"return_node":1,"fraction":{fraction}}}]}}"#
        ))
        .unwrap(),
    );
    input
}

fn history(result: &J) -> &[J] {
    result
        .path(&["transient", "history"])
        .unwrap()
        .as_array()
        .unwrap()
}

fn rows(value: &mut J) -> &mut Vec<J> {
    let J::Array(rows) = value else {
        panic!("expected array")
    };
    rows
}

fn plateau(fraction: f64) -> (f64, f64, f64) {
    // C=1 W/K, Ts=350 K, Tf=300 K, b=exp(-hA/C).
    // Tout=(1-b)Ts+b*Tin and Tin=(1-r)Tf+r*Tout close independently.
    let area = (3.0 + 3.0_f64.sqrt()) / 2.0;
    let decay = (-2.0 * area).exp();
    let exhaust = 300.0 + 50.0 * (1.0 - decay) / (1.0 - fraction * decay);
    let mixed = (1.0 - fraction) * 300.0 + fraction * exhaust;
    (mixed, exhaust, (1.0 - fraction) * (exhaust - 300.0))
}

fn expected_h(time: f64, fraction: f64) -> [f64; 4] {
    let area = (3.0 + 3.0_f64.sqrt()) / 2.0;
    let flux = plateau(fraction).2 / area;
    let mass = RHO / 24.0;
    [3.0_f64.sqrt() / 2.0, 0.5, 0.5, 0.5].map(|opposite| {
        2000.0 + 6000.0 * time.min(0.4) / RHO - time * (area - opposite) / 3.0 * flux / mass
    })
}

fn check_accepted_energy(result: &J) {
    let trajectory = result.get("transient").unwrap();
    let samples = history(result);
    assert_eq!(samples.len(), n(trajectory, "steps") as usize + 1);
    let mut fresh = 0.0;
    for sample in &samples[1..] {
        let report = sample.get("recirculation").unwrap();
        fresh += n(sample, "dt_s") * n(report, "external_heat_gain_w");
        near(n(report, "heat_imbalance_w"), 0.0, 1e-7);
        assert!(n(report, "max_mixing_residual_k") <= 1e-9);
    }
    near(n(trajectory, "fresh_exhaust_energy_gain_j"), fresh, 1e-8);
    let residual =
        n(trajectory, "stored_energy_change_j") - n(trajectory, "input_energy_j") + fresh;
    near(
        n(trajectory, "fresh_exhaust_energy_residual_j"),
        residual,
        1e-8,
    );
    near(residual, 0.0, n(trajectory, "time_s") * 1e-7);
}

#[test]
fn zero_return_preserves_once_through_fields_and_accepted_history_exactly() {
    let zero = request(0.0);
    let mut once = zero.clone();
    remove(&mut once, "recirculation");
    let once = run(&once);
    let zero = run(&zero);
    assert_eq!(
        zero.get("recirculation").unwrap().str_field("status"),
        Some("once-through-zero-returns")
    );
    for field in [
        "solid_temperatures_k",
        "solid_specific_enthalpies_j_kg",
        "solid_liquid_mass_fractions",
        "objective",
    ] {
        assert_eq!(
            zero.get(field).unwrap(),
            once.get(field).unwrap(),
            "{field}"
        );
    }
    let actual = zero.get("transient").unwrap();
    let expected = once.get("transient").unwrap();
    for field in [
        "steps",
        "input_energy_j",
        "stored_energy_change_j",
        "air_energy_gain_j",
        "energy_residual_j",
        "sampled_peak_objective_k",
        "sampled_peak_time_s",
    ] {
        assert_eq!(actual.get(field).unwrap(), expected.get(field).unwrap());
    }
    near(
        n(actual, "fresh_exhaust_energy_gain_j"),
        n(expected, "air_energy_gain_j"),
        0.0,
    );
    near(
        n(actual, "fresh_exhaust_energy_residual_j"),
        n(expected, "energy_residual_j"),
        0.0,
    );
    assert_eq!(history(&zero).len(), history(&once).len());
    for (index, (actual, expected)) in history(&zero).iter().zip(history(&once)).enumerate() {
        if index > 0 {
            assert_eq!(
                actual.get("recirculation").unwrap().str_field("status"),
                Some("once-through-zero-returns")
            );
        }
        let mut physical = actual.clone();
        remove(&mut physical, "recirculation");
        assert_eq!(&physical, expected);
    }
}

#[test]
fn latent_return_matches_analytic_mixing_nodal_storage_and_fresh_exhaust_energy() {
    let result = run(&request(RETURN_FRACTION));
    let (mixed, exhaust, heat) = plateau(RETURN_FRACTION);
    let report = result.get("recirculation").unwrap();
    assert_eq!(report.str_field("status"), Some("solved"));
    let supplies = report.get("supplies").unwrap().as_array().unwrap();
    let streams = report.get("streams").unwrap().as_array().unwrap();
    assert_eq!(supplies.len(), 1);
    assert_eq!(streams.len(), 1);
    near(n(&supplies[0], "fresh_temperature_k"), 300.0, 0.0);
    near(n(&supplies[0], "mixed_temperature_k"), mixed, 1e-8);
    near(n(&supplies[0], "fresh_capacity_w_per_k"), 0.4, 1e-12);
    near(n(&supplies[0], "return_capacity_w_per_k"), 0.6, 1e-12);
    near(n(&streams[0], "temperature_k"), exhaust, 1e-8);
    near(n(report, "external_heat_gain_w"), heat, 1e-7);
    for sample in &history(&result)[1..] {
        let mixing = sample.get("recirculation").unwrap();
        let supplies = mixing.get("mixed_supplies").unwrap().as_array().unwrap();
        near(n(&supplies[0], "mixed_temperature_k"), mixed, 1e-8);
        near(n(sample, "air_heat_gain_w"), heat, 1e-7);
        let expected = expected_h(n(sample, "time_s"), RETURN_FRACTION);
        near(
            n(sample, "minimum_specific_enthalpy_j_kg"),
            expected.into_iter().fold(f64::INFINITY, f64::min),
            2e-6,
        );
        near(
            n(sample, "maximum_specific_enthalpy_j_kg"),
            expected.into_iter().fold(f64::NEG_INFINITY, f64::max),
            2e-6,
        );
    }
    for ((actual, liquid), expected) in values(&result, "solid_specific_enthalpies_j_kg")
        .iter()
        .zip(values(&result, "solid_liquid_mass_fractions"))
        .zip(expected_h(0.8, RETURN_FRACTION))
    {
        assert!((1000.0..3000.0).contains(&expected));
        near(*actual, expected, 2e-6);
        near(liquid, (expected - 1000.0) / 2000.0, 1e-9);
    }
    for temperature in values(&result, "solid_temperatures_k") {
        near(temperature, 350.0, 1e-9);
    }
    let trajectory = result.get("transient").unwrap();
    near(n(trajectory, "input_energy_j"), 400.0, 1e-8);
    near(
        n(trajectory, "fresh_exhaust_energy_gain_j"),
        0.8 * heat,
        1e-7,
    );
    near(
        n(trajectory, "stored_energy_change_j"),
        400.0 - 0.8 * heat,
        1e-6,
    );
    check_accepted_energy(&result);
}

#[test]
fn adaptive_repeated_cycles_sum_only_accepted_fresh_exhaust_energy() {
    let mut input = request(RETURN_FRACTION);
    let schedule = member(&mut input, "transient");
    put(schedule, "max_steps", number(32.0));
    put(
        schedule,
        "adaptive",
        J::parse(
            r#"{"absolute_tolerance_k":1000000,
          "absolute_specific_enthalpy_tolerance_j_kg":0.01,"relative_tolerance":0,
          "minimum_trial_step_s":0.000001,"max_trials":64}"#,
        )
        .unwrap(),
    );
    put(
        schedule,
        "repeat",
        J::parse(r#"{"cycles":2,"max_total_steps":32}"#).unwrap(),
    );
    let result = run(&input);
    let trajectory = result.get("transient").unwrap();
    let repeated = result.get("repeated_cycles").unwrap();
    let heat = plateau(RETURN_FRACTION).2;
    assert_eq!(n(repeated, "cycles_completed"), 2.0);
    assert_eq!(n(repeated, "total_accepted_steps"), 16.0);
    assert_eq!(n(trajectory, "steps"), 8.0);
    let accepted_work: f64 = history(&result)[1..]
        .iter()
        .map(|sample| n(sample, "coupling_iterations"))
        .sum();
    assert!(n(trajectory, "forward_solid_solves") > accepted_work);
    near(n(repeated, "fresh_exhaust_energy_gain_j"), 1.6 * heat, 2e-7);
    near(n(repeated, "fresh_exhaust_energy_residual_j"), 0.0, 1.6e-7);
    near(n(repeated, "input_energy_j"), 800.0, 1e-7);
    near(
        n(repeated, "stored_energy_change_j"),
        800.0 - 1.6 * heat,
        2e-6,
    );
    let cycles = repeated.get("cycles").unwrap().as_array().unwrap();
    assert_eq!(cycles.len(), 2);
    for cycle in cycles {
        near(n(cycle, "fresh_exhaust_energy_gain_j"), 0.8 * heat, 1e-7);
        near(n(cycle, "fresh_exhaust_energy_residual_j"), 0.0, 0.8e-7);
    }
    for (actual, once) in values(&result, "solid_specific_enthalpies_j_kg")
        .iter()
        .zip(expected_h(0.8, RETURN_FRACTION))
    {
        near(*actual, 2000.0 + 2.0 * (once - 2000.0), 4e-6);
    }
    check_accepted_energy(&result);
}

#[test]
fn mixed_phase_fresh_inlet_and_power_adjoints_include_return_feedback() {
    let mut input = request(RETURN_FRACTION);
    put(
        &mut input,
        "objective",
        J::parse(r#"{"mean_wall_region":"wall","gradient":false}"#).unwrap(),
    );
    let schedule = member(&mut input, "transient");
    put(schedule, "max_step_s", number(0.125));
    put(
        schedule,
        "intervals",
        J::parse(
            r#"[{"duration_s":0.25,"power_scale":1},
                {"duration_s":0.25,"power_scale":0.25}]"#,
        )
        .unwrap(),
    );
    let chart = member(schedule, "enthalpy");
    remove(chart, "initial_specific_enthalpy_j_kg");
    put(
        chart,
        "initial_specific_enthalpies_j_kg",
        J::parse("[500,2000,2100,2200]").unwrap(),
    );
    let plain = run(&input);
    let temperatures = values(&plain, "solid_temperatures_k");
    assert!((300.0..350.0).contains(&temperatures[0]));
    for temperature in &temperatures[1..] {
        near(*temperature, 350.0, 1e-9);
    }
    let mixed = |sample: &J| {
        n(
            &sample
                .path(&["recirculation", "mixed_supplies"])
                .unwrap()
                .as_array()
                .unwrap()[0],
            "mixed_temperature_k",
        )
    };
    assert!((mixed(history(&plain).last().unwrap()) - mixed(&history(&plain)[1])).abs() > 1e-4);
    check_accepted_energy(&plain);
    let mut differentiated = input.clone();
    put(
        member(&mut differentiated, "transient"),
        "adjoint",
        J::parse(r#"{"qoi":"final","max_checkpoint_bytes":1048576}"#).unwrap(),
    );
    let differentiated = run(&differentiated);
    assert_eq!(history(&plain), history(&differentiated));
    assert_eq!(
        plain.get("solid_specific_enthalpies_j_kg"),
        differentiated.get("solid_specific_enthalpies_j_kg")
    );
    let adjoint = differentiated.path(&["transient", "adjoint"]).unwrap();
    let objective = |result: &J| n(result.get("objective").unwrap(), "value_k");
    near(n(adjoint, "value_k"), objective(&plain), 1e-9);
    let inlets = values(adjoint, "dtemperature_dinlet_temperatures");
    assert_eq!(inlets.len(), 2);
    assert!(inlets[0] > 1e-5, "fresh inlet sensitivity must be nonzero");
    assert_eq!(inlets[1], 0.0);
    let mut plus = input.clone();
    let mut minus = input.clone();
    for (input, delta) in [(&mut plus, 0.02), (&mut minus, -0.02)] {
        let boundaries = rows(member(member(input, "hydraulics"), "boundaries"));
        put(&mut boundaries[0], "temperature_k", number(300.0 + delta));
    }
    near(
        inlets[0],
        (objective(&run(&plus)) - objective(&run(&minus))) / 0.04,
        2e-6,
    );
    let intervals = adjoint.get("intervals").unwrap().as_array().unwrap();
    assert_eq!(intervals.len(), 2);
    for (index, interval) in intervals.iter().enumerate() {
        let derivative = n(interval, "dtemperature_dpower_multiplier_k");
        assert!(
            derivative > 1e-3,
            "active interval must affect the sensible node"
        );
        let mut plus = input.clone();
        let mut minus = input.clone();
        for (input, multiplier) in [(&mut plus, 1.001), (&mut minus, 0.999)] {
            let intervals = rows(member(member(input, "transient"), "intervals"));
            let scale = n(&intervals[index], "power_scale");
            put(
                &mut intervals[index],
                "power_scale",
                number(scale * multiplier),
            );
        }
        let finite_difference = (objective(&run(&plus)) - objective(&run(&minus))) / 0.002;
        near(
            derivative,
            finite_difference,
            2e-5 * finite_difference.abs().max(1.0),
        );
    }
}
