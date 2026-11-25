//! G1/G3: an equivalent-circuit battery scene discharges into a resistor
//! branch on its own closed form.

use super::scene_context;

/// A flat OCV `E` behind `Rs` and the `Rp`/`Cp` pair, on the bus capacitor,
/// discharging into `R_L` whose far end is held at 0 V.
const E: f64 = 12.0;
const RS: f64 = 0.05;
const RP: f64 = 0.02;
const RL: f64 = 1.93;
const CP: f64 = 100.0;
const CBUS: f64 = 1.0;
const CAPACITY_C: f64 = 3600.0;
const SOC0: f64 = 0.9;
/// The heat capacities \[J/K\] of the series resistance's body, the
/// polarization resistance's body and the load's body, all from 300 K.
const BODIES: [f64; 3] = [50.0, 25.0, 1000.0];

fn close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * (1.0 + expected.abs()),
        "{what}: {actual} != {expected}"
    );
}

/// The scene at the stationary start of the discharge current `E / (Rs + Rp
/// + R_L)`, with `extra` added to the battery's fields.
fn scene(extra: &str) -> String {
    let current = E / (RS + RP + RL);
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 10.0, "step_s": 0.1,
  "elements": {{
    "cell": {{"kind": "battery", "capacity_c": {CAPACITY_C}, "soc_voltage": [[0, {E}], [1, {E}]],
      "series_resistance_ohm": {RS}, "polarization_resistance_ohm": {RP},
      "polarization_capacitance_f": {CP}, "bus_capacitance_f": {CBUS},
      "soc": {SOC0}, "polarization_v": {polarization}, "bus_v": {bus},
      "energy_tolerance_j": 1e-9{extra}}},
    "load": {{"kind": "resistor_branch", "resistance_ohm": {RL}}},
    "series_body": {{"kind": "heat_capacity", "capacity_j_k": {s}, "temperature_k": 300.0}},
    "polarization_body": {{"kind": "heat_capacity", "capacity_j_k": {p}, "temperature_k": 300.0}},
    "load_body": {{"kind": "heat_capacity", "capacity_j_k": {l}, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"common_effort": ["cell.terminal", "load.a"]}},
    {{"prescribe": "load.b", "voltage_v": 0}},
    {{"join": ["cell.heat_series", "series_body.heat"]}},
    {{"join": ["cell.heat_polarization", "polarization_body.heat"]}},
    {{"route": "load.heat", "to": [["load_body.heat", 1.0]]}}
  ]
}}"#,
        polarization = RP * current,
        bus = E - (RS + RP) * current,
        s = BODIES[0],
        p = BODIES[1],
        l = BODIES[2],
    )
}

fn probe(scene: &fs_cli::SystemScene, path: &str) -> f64 {
    scene
        .probe(path)
        .unwrap()
        .unwrap_or_else(|| panic!("{path} is null"))
}

#[test]
fn a_flat_battery_holds_its_stationary_discharge_and_heats_each_body_by_its_resistance() {
    // With a flat OCV the polarization voltage Rp i and bus voltage
    // E - (Rs + Rp) i are stationary under the load current i: only the
    // chemical charge moves, by i t, and each resistance R takes R i^2.
    let current = E / (RS + RP + RL);
    let mut scene = fs_cli::SystemScene::from_text(&scene(""), "battery.json").unwrap();
    scene_context(|cx, _| scene.run(cx).unwrap());
    let t = scene.time_s();
    close(t, 10.0, "time");
    let drawn = current * t;
    close(
        probe(&scene, "cell.remaining_charge_c"),
        SOC0 * CAPACITY_C - drawn,
        "charge",
    );
    close(probe(&scene, "cell.soc"), SOC0 - drawn / CAPACITY_C, "soc");
    // The two capacitor charges, Cp Rp i and Cbus (E - (Rs + Rp) i), are
    // the owner's second and third coordinates.
    close(
        probe(&scene, "cell.polarization_charge_c"),
        CP * RP * current,
        "polarization charge",
    );
    close(
        probe(&scene, "cell.bus_charge_c"),
        CBUS * (E - (RS + RP) * current),
        "bus charge",
    );
    close(
        probe(&scene, "cell.polarization_v"),
        RP * current,
        "polarization",
    );
    close(probe(&scene, "cell.bus_v"), E - (RS + RP) * current, "bus");
    close(
        probe(&scene, "cell.terminal_v"),
        E - (RS + RP) * current,
        "terminal",
    );
    // The terminal's current is positive into the battery.
    close(probe(&scene, "cell.current_a"), -current, "current");
    let rates = [RS, RP, RL].map(|r| r * current * current);
    close(probe(&scene, "cell.heat_series_w"), rates[0], "series heat");
    close(
        probe(&scene, "cell.heat_polarization_w"),
        rates[1],
        "polarization heat",
    );
    for (body, (rate, capacity)) in ["series_body", "polarization_body", "load_body"]
        .into_iter()
        .zip(rates.into_iter().zip(BODIES))
    {
        close(
            probe(&scene, &format!("{body}.temperature_k")),
            300.0 + rate * t / capacity,
            body,
        );
    }
    assert!(scene.ledger().energy_residual_j().abs() < 1e-9);
    let report = scene.render_json();
    assert!(
        report.contains("a battery is stepped by its fs-system owner on its authored curve"),
        "{report}"
    );
    // The controls the run used are reported as given: here fs-system's
    // defaults.
    for (control, value) in [
        ("solve_relative", 1e-13),
        ("solve_absolute", 1e-12),
        ("solve_max_iterations", 40.0),
        ("backtracking_halvings", 10.0),
        ("energy_tolerance_j", 1e-9),
    ] {
        assert_eq!(
            probe(&scene, &format!("cell.{control}")).to_bits(),
            f64::to_bits(value),
            "{control}"
        );
    }
}

/// The scene stepped adaptively at `order`, with a tolerance for every
/// account except `without`.
fn adaptive(order: u32, without: &str) -> String {
    let tolerances = [
        "cell.remaining_charge_c",
        "cell.polarization_charge_c",
        "cell.bus_charge_c",
        "series_body.temperature_k",
        "polarization_body.temperature_k",
        "load_body.temperature_k",
    ]
    .into_iter()
    .filter(|name| *name != without)
    .map(|name| format!(r#""{name}": {{"absolute": 1e-9, "relative": 1e-9}}"#))
    .collect::<Vec<_>>()
    .join(", ");
    scene("").replace(
        r#""step_s": 0.1,"#,
        &format!(
            r#""step_s": 0.1, "adaptive": {{"order": {order}, "tolerances": {{{tolerances}}}}},"#
        ),
    )
}

#[test]
fn an_adaptive_battery_accounts_its_three_charges_at_order_two() {
    let current = E / (RS + RP + RL);
    let mut scene = fs_cli::SystemScene::from_text(&adaptive(2, ""), "battery.json").unwrap();
    scene_context(|cx, _| scene.run(cx).unwrap());
    close(scene.time_s(), 10.0, "time");
    let drawn = current * 10.0;
    close(
        probe(&scene, "cell.remaining_charge_c"),
        SOC0 * CAPACITY_C - drawn,
        "charge",
    );
    close(probe(&scene, "cell.bus_v"), E - (RS + RP) * current, "bus");
    for (bad, expected) in [
        (adaptive(3, ""), "is a battery, whose step is of order 2"),
        (
            adaptive(2, "cell.bus_charge_c"),
            "cell.bus_charge_c: missing explicit tolerance",
        ),
    ] {
        let error = fs_cli::SystemScene::from_text(&bad, "battery.json")
            .err()
            .unwrap_or_else(|| panic!("admitted: {expected}"))
            .to_string();
        assert!(error.contains(expected), "{expected}: {error}");
    }
}

#[test]
fn authored_battery_controls_are_reported_and_bad_batteries_are_refused() {
    let authored = scene(
        r#", "solve": {"relative": 1e-12, "absolute": 1e-11, "max_iterations": 25},
      "backtracking_halvings": 4"#,
    );
    let admitted = fs_cli::SystemScene::from_text(&authored, "battery.json").unwrap();
    for (control, value) in [
        ("solve_relative", 1e-12),
        ("solve_absolute", 1e-11),
        ("solve_max_iterations", 25.0),
        ("backtracking_halvings", 4.0),
    ] {
        assert_eq!(
            probe(&admitted, &format!("cell.{control}")).to_bits(),
            f64::to_bits(value),
            "{control}"
        );
    }
    // Before the first step there is no last step.
    assert_eq!(admitted.probe("cell.current_a").unwrap(), None);
    let good = scene("");
    for (bad, expected) in [
        (
            good.replace("[[0, 12], [1, 12]]", "[[0, 12], [0.5, 12]]"),
            "battery owner refused",
        ),
        (
            good.replace("\"soc\": 0.9", "\"soc\": 1.5"),
            "battery owner refused",
        ),
        (
            good.replace(",\n      \"energy_tolerance_j\": 1e-9", ""),
            "missing energy_tolerance_j",
        ),
        (
            good.replace("\"soc\": 0.9", "\"soc\": 0.9, \"temperature_law\": 1"),
            "unsupported field \"temperature_law\"",
        ),
        (
            good.replace(
                "\"energy_tolerance_j\": 1e-9",
                "\"energy_tolerance_j\": 1e-9, \"backtracking_halvings\": 2.5",
            ),
            "backtracking_halvings",
        ),
    ] {
        assert_ne!(bad, good, "the replacement must change the scene");
        let error = fs_cli::SystemScene::from_text(&bad, "battery.json")
            .err()
            .unwrap_or_else(|| panic!("admitted: {expected}"))
            .to_string();
        assert!(error.contains("elements.cell"), "{error}");
        assert!(error.contains(expected), "{expected}: {error}");
    }
}
