//! G3 complete-trajectory controls with original P1 footprints and mixed
//! sensible/latent storage. All finite differences rerun the actual binary.
use super::*;

const COMPONENTS: [&str; 3] = ["chip", "memory", "standby"];
// Small relative to the 60--250 W applied loads, but large enough to resolve
// the weaker cross-contact response above the 1e-9 K coupling tolerance.
const WATT_STEP: f64 = 0.1;
const WATT_TOLERANCE: f64 = 3e-7;

fn rows(value: &mut J) -> &mut Vec<J> {
    let J::Array(rows) = value else {
        panic!("expected array")
    };
    rows
}

fn request(qoi: Option<&str>) -> J {
    let mut input = J::parse(CONTACT_FIXTURE).unwrap();
    let solid = member(&mut input, "solid");
    remove(solid, "source_w_m3");
    put(
        solid,
        "component_power",
        J::parse(
            r#"{
      "total_w":280,"relative_tolerance":1e-12,"components":[
        {"name":"chip","watts":200,"vertices":[1]},
        {"name":"memory","watts":80,"vertices":[4,6]},
        {"name":"standby","watts":0,"vertices":[1,3]}
      ]}"#,
        )
        .unwrap(),
    );
    put(
        &mut input,
        "objective",
        J::parse(r#"{"mean_wall_region":"first-wall","gradient":false}"#).unwrap(),
    );
    let schedule = member(&mut input, "transient");
    // Binary-exact clocks make repeated and flattened schedules identical.
    // dt != 1 also exposes an extra time factor in either load contraction.
    put(schedule, "max_step_s", number(0.03125));
    put(schedule, "max_steps", number(16.0));
    put(
        schedule,
        "intervals",
        J::parse(
            r#"[
      {"duration_s":0.0625,"power_scale":1.25},
      {"duration_s":0.0625,"component_powers_w":{"chip":120,"memory":60,"standby":0}}
    ]"#,
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
    // Nonzero df/dh changes heat transfer from plateau nodes to the sensible
    // node, despite their local dT/dh=0. Radiation remains in the fixture.
    put(
        &mut materials[0],
        "phase_conductivity",
        J::parse(
            r#"{
      "law":"linear-liquid-mass-fraction","solid_multiplier":0.5,
      "liquid_multiplier":4,"source":"Synthetic phase-conductivity control fixture"
    }"#,
        )
        .unwrap(),
    );
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
    if let Some(qoi) = qoi {
        put(
            schedule,
            "adjoint",
            J::parse(&format!(
                r#"{{
          "qoi":"{qoi}","max_checkpoint_bytes":1048576,
          "component_power":true,"contact_resistance":true
        }}"#
            ))
            .unwrap(),
        );
    }
    input
}

fn primal(input: &J) -> J {
    let mut input = input.clone();
    remove(member(&mut input, "transient"), "adjoint");
    input
}

fn section(result: &J) -> &J {
    result
        .get("repeated_cycles")
        .unwrap_or_else(|| result.get("transient").unwrap())
}

fn adjoint(result: &J) -> &J {
    section(result).get("adjoint").unwrap()
}

fn value(result: &J, qoi: &str) -> f64 {
    if qoi == "final" {
        n(result.get("objective").unwrap(), "value_k")
    } else {
        n(section(result), "sampled_peak_objective_k")
    }
}

fn component<'a>(result: &'a J, interval: usize, name: &str) -> &'a J {
    adjoint(result)
        .path(&["component_power_sensitivities", "intervals"])
        .unwrap()
        .as_array()
        .unwrap()[interval]
        .get("rows")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field("component") == Some(name))
        .unwrap()
}

fn contact(result: &J) -> &J {
    let report = adjoint(result)
        .get("contact_resistance_sensitivities")
        .unwrap();
    assert_eq!(
        report.str_field("method"),
        Some("trajectory-coupled-adjoint-contact-bilinear-form")
    );
    let rows = report.get("rows").unwrap().as_array().unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].str_field("contact"), Some("bondline"));
    &rows[0]
}

fn applied(input: &J, interval: usize) -> J {
    let row = &input
        .path(&["transient", "intervals"])
        .unwrap()
        .as_array()
        .unwrap()[interval];
    if let Some(powers) = row.get("component_powers_w") {
        return powers.clone();
    }
    let scale = n(row, "power_scale");
    J::Object(
        input
            .path(&["solid", "component_power", "components"])
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .map(|row| {
                (
                    row.str_field("name").unwrap().into(),
                    number(scale * n(row, "watts")),
                )
            })
            .collect(),
    )
}

fn change_interval(input: &mut J, interval: usize, name: &str, watts: f64) {
    let mut powers = applied(input, interval);
    put(&mut powers, name, number(watts));
    let row = &mut rows(member(member(input, "transient"), "intervals"))[interval];
    remove(row, "power_scale");
    put(row, "component_powers_w", powers);
}

fn change_base(input: &mut J, name: &str, watts: f64) {
    let power = member(member(input, "solid"), "component_power");
    let components = rows(member(power, "components"));
    let row = components
        .iter_mut()
        .find(|row| row.str_field("name") == Some(name))
        .unwrap();
    put(row, "watts", number(watts));
    let total = components.iter().map(|row| n(row, "watts")).sum();
    put(power, "total_w", number(total));
}

fn interval_fd(input: &J, qoi: &str, interval: usize, name: &str) -> f64 {
    let baseline = n(&applied(input, interval), name);
    let mut plus = primal(input);
    let mut minus = primal(input);
    change_interval(&mut plus, interval, name, baseline + WATT_STEP);
    // Zero watts has an admissible one-sided derivative; the original map
    // and all other components stay unchanged in both complete trajectories.
    let denominator = if baseline >= WATT_STEP {
        change_interval(&mut minus, interval, name, baseline - WATT_STEP);
        2.0 * WATT_STEP
    } else {
        WATT_STEP
    };
    (value(&run(&plus), qoi) - value(&run(&minus), qoi)) / denominator
}

fn contact_fd(input: &J, qoi: &str) -> f64 {
    let epsilon = 1e-3_f64;
    let mut plus = primal(input);
    let mut minus = primal(input);
    for (input, scale) in [(&mut plus, epsilon.exp()), (&mut minus, (-epsilon).exp())] {
        let row = &mut rows(member(member(input, "solid"), "contacts"))[0];
        put(row, "resistance_m2_k_w", number(0.1 * scale));
    }
    (value(&run(&plus), qoi) - value(&run(&minus), qoi)) / (2.0 * epsilon)
}

fn same_history(first: &J, second: &J) {
    for path in [
        &["solid_specific_enthalpies_j_kg"][..],
        &["solid_temperatures_k"][..],
        &["transient", "history"][..],
    ] {
        assert_eq!(first.path(path), second.path(path));
    }
    assert_eq!(section(first).get("cycles"), section(second).get("cycles"));
}

#[test]
fn absolute_and_base_watts_preserve_overlapping_and_zero_power_footprints() {
    let input = request(Some("final"));
    let result = run(&input);
    same_history(&result, &run(&primal(&input)));
    let h = values(&result, "solid_specific_enthalpies_j_kg");
    assert!((0.0..1000.0).contains(&h[0]));
    assert!(h[1..4].iter().all(|h| (1000.0..3000.0).contains(h)));
    assert!(
        values(&result, "solid_liquid_mass_fractions")[4..]
            .iter()
            .all(|&f| f == 0.0)
    );
    let report = adjoint(&result)
        .get("component_power_sensitivities")
        .unwrap();
    assert_eq!(
        report.str_field("method"),
        Some("consistent-P1-source-transpose")
    );
    for interval in 0..2 {
        let mut relative_sum = 0.0;
        for name in COMPONENTS {
            let row = component(&result, interval, name);
            let slope = n(row, "dtemperature_dpower_w_k_per_w");
            near(
                n(row, "applied_power_w"),
                n(&applied(&input, interval), name),
                1e-12,
            );
            near(
                slope,
                interval_fd(&input, "final", interval, name),
                WATT_TOLERANCE,
            );
            let relative = n(row, "dtemperature_dpower_multiplier_k");
            near(relative, slope * n(row, "applied_power_w"), 1e-12);
            relative_sum += relative;
        }
        let old = &adjoint(&result)
            .get("intervals")
            .unwrap()
            .as_array()
            .unwrap()[interval];
        near(
            relative_sum,
            n(old, "dtemperature_dpower_multiplier_k"),
            1e-8,
        );
        let dormant = component(&result, interval, "standby");
        assert_eq!(n(dormant, "applied_power_w"), 0.0);
        assert_eq!(n(dormant, "dtemperature_dpower_multiplier_k"), 0.0);
        assert!(n(dormant, "dtemperature_dpower_w_k_per_w").abs() > 1e-5);
    }
    let base = report.get("base_rows").unwrap().as_array().unwrap();
    assert_eq!(base.len(), COMPONENTS.len());
    for row in base {
        let name = row.str_field("component").unwrap();
        let watts = n(row, "base_power_w");
        // Only the first interval scales base watts; the absolute second
        // interval ignores them, including the updated declared system total.
        let slope = n(row, "dtemperature_dbase_power_w_k_per_w");
        near(
            slope,
            1.25 * n(component(&result, 0, name), "dtemperature_dpower_w_k_per_w"),
            1e-12,
        );
        let mut plus = primal(&input);
        let mut minus = primal(&input);
        change_base(&mut plus, name, watts + WATT_STEP);
        let denominator = if watts >= WATT_STEP {
            change_base(&mut minus, name, watts - WATT_STEP);
            2.0 * WATT_STEP
        } else {
            WATT_STEP
        };
        near(
            slope,
            (value(&run(&plus), "final") - value(&run(&minus), "final")) / denominator,
            WATT_TOLERANCE,
        );
    }
}

#[test]
fn log_contact_resistance_differentiates_mixed_storage_and_radiative_history() {
    let input = request(Some("final"));
    let result = run(&input);
    let row = contact(&result);
    let slope = n(row, "dtemperature_dlog_resistance_k");
    assert!(slope.abs() > 1e-4);
    near(n(row, "resistance_m2_k_w"), 0.1, 0.0);
    near(n(row, "dtemperature_dresistance_w_m2"), slope / 0.1, 1e-12);
    near(slope, contact_fd(&input, "final"), 2e-5);
    assert!(n(result.get("transient").unwrap(), "radiative_energy_loss_j") > 1.0);

    let mut ordinary = input.clone();
    let options = member(member(&mut ordinary, "transient"), "adjoint");
    remove(options, "component_power");
    remove(options, "contact_resistance");
    let existing = run(&ordinary);
    same_history(&existing, &result);
    assert!(
        adjoint(&existing)
            .get("component_power_sensitivities")
            .is_none()
    );
    assert!(
        adjoint(&existing)
            .get("contact_resistance_sensitivities")
            .is_none()
    );
    for key in [
        "reconstruction_solid_solves",
        "reconstructed_solid_endpoints",
        "adjoint_sweeps",
    ] {
        assert_eq!(
            adjoint(&existing).get(key),
            adjoint(&result).get(key),
            "no solve per control"
        );
    }
    put(
        member(member(&mut ordinary, "transient"), "adjoint"),
        "contact_resistance",
        J::Bool(true),
    );
    let contact_only = run(&ordinary);
    assert_eq!(
        adjoint(&contact_only).get("component_power_sensitivities"),
        Some(&J::Null)
    );
    near(
        n(contact(&contact_only), "dtemperature_dlog_resistance_k"),
        slope,
        1e-12,
    );
}

#[test]
fn repeated_component_watts_and_persistent_contacts_match_flattened_controls() {
    let qoi = "sampled-peak";
    let mut input = request(Some(qoi));
    put(
        member(&mut input, "transient"),
        "repeat",
        J::parse(r#"{"cycles":3,"max_total_steps":12}"#).unwrap(),
    );
    let result = run(&input);
    let mut flat = input.clone();
    let schedule = member(&mut flat, "transient");
    remove(schedule, "repeat");
    let intervals = rows(member(schedule, "intervals"));
    let cycle = intervals.clone();
    intervals.extend(cycle.iter().cloned());
    intervals.extend(cycle);
    let flattened = run(&flat);
    for key in ["solid_specific_enthalpies_j_kg", "solid_temperatures_k"] {
        assert_eq!(result.get(key), flattened.get(key));
    }
    assert_eq!(result.path(&["transient", "adjoint"]), Some(&J::Null));
    assert_eq!(n(adjoint(&result), "cycles"), 3.0);
    assert_eq!(n(section(&result), "total_accepted_steps"), 12.0);
    assert!(n(adjoint(&result), "state_index") > 8.0);
    near(value(&result, qoi), value(&flattened, qoi), 1e-10);
    for interval in 0..2 {
        for name in COMPONENTS {
            let key = "dtemperature_dpower_w_k_per_w";
            near(
                n(component(&result, interval, name), key),
                (0..3)
                    .map(|cycle| n(component(&flattened, 2 * cycle + interval, name), key))
                    .sum(),
                1e-8,
            );
        }
    }
    near(
        n(contact(&result), "dtemperature_dlog_resistance_k"),
        n(contact(&flattened), "dtemperature_dlog_resistance_k"),
        1e-8,
    );
    // Each perturbation affects every occurrence of the selected base interval.
    for (interval, name) in [(0, "chip"), (1, "standby")] {
        near(
            n(
                component(&result, interval, name),
                "dtemperature_dpower_w_k_per_w",
            ),
            interval_fd(&input, qoi, interval, name),
            WATT_TOLERANCE,
        );
    }
    near(
        n(contact(&result), "dtemperature_dlog_resistance_k"),
        contact_fd(&input, qoi),
        2e-5,
    );
}

#[test]
fn initial_peaks_zero_future_controls_and_control_storage_is_admitted_before_output() {
    let mut initial_peak = request(Some("sampled-peak"));
    put(
        &mut initial_peak,
        "objective",
        J::parse(r#"{"max_solid_temperature":true,"gradient":false}"#).unwrap(),
    );
    let schedule = member(&mut initial_peak, "transient");
    put(
        member(schedule, "enthalpy"),
        "initial_specific_enthalpies_j_kg",
        J::parse("[4000,2000,2100,2200,1500,1700,1900,1600]").unwrap(),
    );
    for interval in rows(member(schedule, "intervals")) {
        remove(interval, "component_powers_w");
        put(interval, "power_scale", number(0.0));
    }
    let result = run(&initial_peak);
    assert_eq!(n(adjoint(&result), "state_index"), 0.0);
    assert_eq!(n(adjoint(&result), "reconstructed_solid_endpoints"), 0.0);
    assert_eq!(n(contact(&result), "dtemperature_dlog_resistance_k"), 0.0);
    for interval in 0..2 {
        for name in COMPONENTS {
            assert_eq!(
                n(
                    component(&result, interval, name),
                    "dtemperature_dpower_w_k_per_w"
                ),
                0.0
            );
        }
    }

    let mut ordinary = request(Some("final"));
    let options = member(member(&mut ordinary, "transient"), "adjoint");
    remove(options, "component_power");
    remove(options, "contact_resistance");
    let bytes = n(adjoint(&run(&ordinary)), "checkpoint_bytes");
    let mut exhausted = request(Some("final"));
    put(
        member(member(&mut exhausted, "transient"), "adjoint"),
        "max_checkpoint_bytes",
        number(bytes),
    );
    let failure = output(&exhausted);
    assert_eq!(failure.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(failure.stdout.is_empty());
    assert!(String::from_utf8_lossy(&failure.stderr).contains("max_checkpoint_bytes"));

    let mut missing = request(Some("final"));
    let solid = member(&mut missing, "solid");
    remove(solid, "component_power");
    put(solid, "source_w_m3", number(0.0));
    for interval in rows(member(member(&mut missing, "transient"), "intervals")) {
        remove(interval, "component_powers_w");
        put(interval, "power_scale", number(0.0));
    }
    assert!(refuses(&missing).contains("component_power"));
}
