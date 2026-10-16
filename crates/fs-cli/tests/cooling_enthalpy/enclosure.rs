//! Actual-binary closed-enclosure exchange through heterogeneous latent storage.
//! The two equal-area traces use prescribed synthetic parallel-plate factors;
//! no visibility calculation or physical validation of the tetrahedra is claimed.

use super::*;

fn rows(value: &mut J) -> &mut Vec<J> {
    let J::Array(rows) = value else {
        panic!("expected array")
    };
    rows
}

fn request() -> J {
    J::parse(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../examples/cooling-network/enthalpy-enclosure.json"
    )))
    .unwrap()
}

fn patch<'a>(result: &'a J, name: &str) -> &'a J {
    result
        .path(&["radiation", "surfaces"])
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field("surface") == Some(name))
        .unwrap()
}

fn history(result: &J) -> &[J] {
    result
        .path(&["transient", "history"])
        .unwrap()
        .as_array()
        .unwrap()
}

fn reference() -> (f64, f64, f64, f64) {
    let area = 1.0 + 3.0_f64.sqrt() / 2.0;
    let exchange =
        area * SIGMA * (350.0_f64.powi(4) - 330.0_f64.powi(4)) / (1.0 / 0.8 + 1.0 / 0.6 - 1.0);
    let effectiveness = 1.0 - (-2.0 * area).exp();
    (area, exchange, 50.0 * effectiveness, 30.0 * effectiveness)
}

fn expected_h(time: f64, powered_time: f64) -> [f64; 8] {
    let (area, radiation, hot_air, cold_air) = reference();
    let diagonal = 3.0_f64.sqrt() / 2.0;
    // Integral of each P1 node over the three radiating/convecting faces.
    let weights = [
        1.0 / 3.0,
        (diagonal + 0.5) / 3.0,
        (diagonal + 0.5) / 3.0,
        area / 3.0,
    ];
    let mut result = [0.0; 8];
    for side in 0..2 {
        let (rho, initial, surface_heat, contact_sign) = if side == 0 {
            (10.0, 2000.0, hot_air + radiation, 1.0)
        } else {
            (20.0, 3000.0, cold_air - radiation, -1.0)
        };
        for node in 0..4 {
            // The matching contact face has area 1/2; its first three nodal
            // loads are (350-330)/R times (area/3). Its fourth node is absent.
            let contact = if node < 3 {
                contact_sign * (20.0 / 0.1) / 6.0
            } else {
                0.0
            };
            result[4 * side + node] = initial + 6000.0 * powered_time / rho
                - time * (weights[node] * surface_heat / area + contact) / (rho / 24.0);
        }
    }
    result
}

fn check_accepted_energy(result: &J) {
    let trajectory = result.get("transient").unwrap();
    let samples = history(result);
    assert_eq!(samples.len(), n(trajectory, "steps") as usize + 1);
    let mut input = 0.0;
    let mut storage = 0.0;
    let mut air = 0.0;
    let mut radiation = 0.0;
    for sample in &samples[1..] {
        let dt = n(sample, "dt_s");
        input += dt * n(sample, "source_w");
        storage += n(sample, "stored_energy_change_j");
        air += dt * n(sample, "air_heat_gain_w");
        radiation += dt * n(sample, "radiative_heat_w");
    }
    for (key, expected) in [
        ("input_energy_j", input),
        ("stored_energy_change_j", storage),
        ("air_energy_gain_j", air),
        ("radiative_energy_loss_j", radiation),
    ] {
        near(n(trajectory, key), expected, 1e-8);
    }
    near(radiation, 0.0, 1e-8);
    near(
        storage - input + air + radiation,
        n(trajectory, "energy_residual_j"),
        1e-8,
    );
    near(
        n(trajectory, "energy_residual_j"),
        0.0,
        1e-7 * n(trajectory, "time_s"),
    );
}

#[test]
fn unequal_latent_temperatures_exchange_gray_heat_without_an_external_sink() {
    let input = request();
    let result = run(&input);
    let (area, exchange, hot_air, cold_air) = reference();
    let radiation = result.get("radiation").unwrap();
    assert_eq!(
        radiation.str_field("model"),
        Some("closed-gray-diffuse-enclosure")
    );
    assert_eq!(
        radiation.str_field("temporal_scope"),
        Some("final-accepted-endpoint")
    );
    assert!(exchange > 100.0);
    for (name, temperature, heat) in [
        ("first-wall", 350.0, exchange),
        ("second-wall", 330.0, -exchange),
    ] {
        let row = patch(&result, name);
        near(n(row, "area_m2"), area, 1e-12);
        near(n(row, "mean_temperature_k"), temperature, 1e-9);
        near(n(row, "outward_heat_w"), heat, 1e-7);
        near(n(row, "applied_heat_w"), heat, 1e-7);
    }
    near(n(radiation, "radiative_out_w"), 0.0, 1e-8);
    check_accepted_energy(&result);
    let trajectory = result.get("transient").unwrap();
    near(n(trajectory, "input_energy_j"), 400.0, 1e-8);
    near(
        n(trajectory, "air_energy_gain_j"),
        0.4 * (hot_air + cold_air),
        1e-7,
    );
    near(
        n(trajectory, "stored_energy_change_j"),
        400.0 - 0.4 * (hot_air + cold_air),
        1e-6,
    );
    for (actual, expected) in values(&result, "solid_specific_enthalpies_j_kg")
        .iter()
        .zip(expected_h(0.4, 0.2))
    {
        near(*actual, expected, 3e-6);
    }
    for (index, temperature) in values(&result, "solid_temperatures_k").iter().enumerate() {
        near(*temperature, if index < 4 { 350.0 } else { 330.0 }, 1e-9);
    }
    let mut absent = input;
    remove(&mut absent, "radiation");
    let absent = run(&absent);
    let actual = values(&result, "solid_specific_enthalpies_j_kg");
    let unexchanged = values(&absent, "solid_specific_enthalpies_j_kg");
    assert!(
        actual[..4]
            .iter()
            .zip(&unexchanged[..4])
            .all(|(a, b)| a < b)
    );
    assert!(
        actual[4..]
            .iter()
            .zip(&unexchanged[4..])
            .all(|(a, b)| a > b)
    );
    near(
        n(trajectory, "stored_energy_change_j"),
        n(absent.get("transient").unwrap(), "stored_energy_change_j"),
        1e-7,
    );
}

fn mixed_request() -> J {
    let mut input = request();
    put(
        &mut input,
        "objective",
        J::parse(r#"{"mean_wall_region":"first-wall","gradient":false}"#).unwrap(),
    );
    let schedule = member(&mut input, "transient");
    put(schedule, "max_step_s", number(0.03125));
    put(
        schedule,
        "intervals",
        J::parse(
            r#"[{"duration_s":0.0625,"power_scale":1},{"duration_s":0.0625,"power_scale":0.5}]"#,
        )
        .unwrap(),
    );
    let storage = member(schedule, "enthalpy");
    put(
        storage,
        "initial_specific_enthalpies_j_kg",
        J::parse("[500,2000,2100,2200,1500,1700,1900,1600]").unwrap(),
    );
    let materials = rows(member(storage, "materials"));
    put(&mut materials[0], "phase_conductivity", J::parse(
        r#"{"law":"linear-liquid-mass-fraction","solid_multiplier":0.5,"liquid_multiplier":4,"source":"Synthetic phase-conductivity enclosure fixture"}"#).unwrap());
    put(&mut materials[1], "phase", J::Str("solid".into()));
    put(
        &mut materials[1],
        "knots",
        J::parse(
            r#"[
        {"specific_enthalpy_j_kg":0,"temperature_k":230,"liquid_mass_fraction":0},
        {"specific_enthalpy_j_kg":7000,"temperature_k":530,"liquid_mass_fraction":0}
    ]"#,
        )
        .unwrap(),
    );
    input
}

#[test]
fn enclosure_total_history_adjoint_matches_initial_inlet_and_power_perturbations() {
    let input = mixed_request();
    let plain = run(&input);
    let mut differentiated = input.clone();
    put(
        member(&mut differentiated, "transient"),
        "adjoint",
        J::parse(r#"{"qoi":"final","max_checkpoint_bytes":1048576}"#).unwrap(),
    );
    let result = run(&differentiated);
    assert_eq!(history(&result), history(&plain));
    assert_eq!(
        result.get("solid_specific_enthalpies_j_kg"),
        plain.get("solid_specific_enthalpies_j_kg")
    );
    assert_eq!(result.get("radiation"), plain.get("radiation"));
    check_accepted_energy(&result);
    let temperatures = values(&result, "solid_temperatures_k");
    assert!(temperatures[0] > 300.0 && temperatures[0] < 350.0);
    for temperature in &temperatures[1..4] {
        near(*temperature, 350.0, 1e-9);
    }
    let adjoint = result.path(&["transient", "adjoint"]).unwrap();
    let value = |result: &J| n(result.get("objective").unwrap(), "value_k");
    let initial = values(adjoint, "dtemperature_dinitial_specific_enthalpies_k_kg_j");
    assert_eq!(initial.len(), 8);
    let inlets = values(adjoint, "dtemperature_dinlet_temperatures");
    assert_eq!(inlets.len(), 2);
    near(inlets[1], 0.0, 0.0);
    let power = n(
        &adjoint.get("intervals").unwrap().as_array().unwrap()[0],
        "dtemperature_dpower_multiplier_k",
    );
    // Vertex seven is not on the contact face: its response reaches the
    // objective through the physical receiver, enclosure and air coupling.
    let gradients = [initial[0], initial[7], inlets[0], power];
    for (control, gradient) in gradients.into_iter().enumerate() {
        assert!(
            gradient.abs() > 1e-7,
            "control {control} must exercise feedback"
        );
        let step = if control == 3 { 0.001 } else { 0.1 };
        let mut plus = input.clone();
        let mut minus = input.clone();
        for (trial, delta) in [(&mut plus, step), (&mut minus, -step)] {
            match control {
                0 | 1 => {
                    let vertex = if control == 0 { 0 } else { 7 };
                    let old = rows(member(
                        member(member(trial, "transient"), "enthalpy"),
                        "initial_specific_enthalpies_j_kg",
                    ));
                    old[vertex] = number(old[vertex].as_f64().unwrap() + delta);
                }
                2 => {
                    let boundaries = rows(member(member(trial, "hydraulics"), "boundaries"));
                    put(&mut boundaries[0], "temperature_k", number(300.0 + delta));
                }
                _ => {
                    let intervals = rows(member(member(trial, "transient"), "intervals"));
                    put(&mut intervals[0], "power_scale", number(1.0 + delta));
                }
            }
        }
        let finite_difference = (value(&run(&plus)) - value(&run(&minus))) / (2.0 * step);
        near(
            gradient,
            finite_difference,
            3e-6 * finite_difference.abs().max(1.0),
        );
    }
}

#[test]
fn adaptive_repetition_carries_accepted_h_and_only_counts_internal_heat_once() {
    let mut input = request();
    let schedule = member(&mut input, "transient");
    put(schedule, "max_steps", number(32.0));
    put(
        schedule,
        "adaptive",
        J::parse(
            r#"{
        "absolute_tolerance_k":1e-7,"absolute_specific_enthalpy_tolerance_j_kg":0.01,
        "relative_tolerance":0,"minimum_trial_step_s":0.0001,"max_trials":64
    }"#,
        )
        .unwrap(),
    );
    put(
        schedule,
        "repeat",
        J::parse(r#"{"cycles":2,"max_total_steps":32}"#).unwrap(),
    );
    let result = run(&input);
    let repeated = result.get("repeated_cycles").unwrap();
    near(n(repeated, "cycles_completed"), 2.0, 0.0);
    near(n(repeated, "total_accepted_steps"), 16.0, 0.0);
    check_accepted_energy(&result);
    let (_, _, hot_air, cold_air) = reference();
    near(n(repeated, "input_energy_j"), 800.0, 1e-8);
    near(
        n(repeated, "air_energy_gain_j"),
        0.8 * (hot_air + cold_air),
        2e-7,
    );
    near(n(repeated, "radiative_energy_loss_j"), 0.0, 1e-8);
    near(n(repeated, "energy_residual_j"), 0.0, 0.8e-7);
    for (actual, expected) in values(&result, "solid_specific_enthalpies_j_kg")
        .iter()
        .zip(expected_h(0.8, 0.4))
    {
        near(*actual, expected, 6e-6);
    }
    // Two separate CLI invocations must replay the same accepted physical h.
    remove(member(&mut input, "transient"), "repeat");
    let first = run(&input);
    put(
        member(member(&mut input, "transient"), "enthalpy"),
        "initial_specific_enthalpies_j_kg",
        first.get("solid_specific_enthalpies_j_kg").unwrap().clone(),
    );
    let second = run(&input);
    assert_eq!(
        result.get("solid_specific_enthalpies_j_kg"),
        second.get("solid_specific_enthalpies_j_kg")
    );
    assert_eq!(result.get("radiation"), second.get("radiation"));
    let trajectory = result.get("transient").unwrap();
    let accepted_work: f64 = history(&result)[1..]
        .iter()
        .map(|row| n(row, "coupling_iterations"))
        .sum();
    assert!(
        n(trajectory, "forward_solid_solves") > accepted_work,
        "adaptive coarse solves must be charged without advancing physical h"
    );
}

#[test]
fn invalid_or_unresolved_enclosure_refuses_without_a_partial_trajectory() {
    let mut nonreciprocal = request();
    put(
        member(member(&mut nonreciprocal, "radiation"), "enclosure"),
        "view_factors",
        J::parse("[[0,1],[0.5,0.5]]").unwrap(),
    );
    let diagnostic = refuses(&nonreciprocal);
    assert!(diagnostic.contains("recipro"), "{diagnostic}");
    let mut singular = request();
    for surface in rows(member(
        member(member(&mut singular, "radiation"), "enclosure"),
        "surfaces",
    )) {
        put(surface, "emissivity", number(1e-20));
    }
    let diagnostic = refuses(&singular);
    assert!(
        diagnostic.contains("singular") || diagnostic.contains("pivot"),
        "{diagnostic}"
    );
    let mut exhausted = mixed_request();
    put(
        member(&mut exhausted, "radiation"),
        "max_iterations",
        number(1.0),
    );
    let diagnostic = refuses(&exhausted);
    assert!(
        diagnostic.contains("enclosure") && diagnostic.contains("exhaust"),
        "{diagnostic}"
    );
}
