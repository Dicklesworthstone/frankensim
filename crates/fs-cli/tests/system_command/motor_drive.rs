//! G1/G3: `motor_drive` scenes step fs-phs's battery-fed commutated drive by
//! its fs-system owner: free shafts keep their momentum through a switching
//! program, both lanes hold the DC operating point of the resistive network
//! in closed form, and bad scenes are refused.

use super::scene_context;

/// The battery: capacity, its two resistances and two capacitances, and
/// the lead.
const CAPACITY_C: f64 = 3600.0;
const RS: f64 = 0.1;
const RP: f64 = 0.1;
const CP: f64 = 2.0;
const CB: f64 = 0.5;
const RLEAD: f64 = 0.1;
/// The converter's inductance, bus capacitances and winding.
const L: f64 = 0.3;
const BUS_F: [f64; 2] = [0.2, 0.4];
const RW: f64 = 0.1;
/// The four switches, input high, input low, output high, output low: on
/// and off resistances, each different so that an order mistake shows.
const ON: [f64; 4] = [0.10, 0.11, 0.12, 0.13];
const OFF: [f64; 4] = [100.0, 120.0, 140.0, 160.0];
/// The switches' output capacitances in the eleven-state lane.
const COSS: [f64; 4] = [0.025, 0.04, 0.035, 0.05];
/// The motor: inductance, winding, the two inertias, coupling and drag.
const LM: f64 = 0.05;
const RM: f64 = 0.1;
const JR: f64 = 0.4;
const JS: f64 = 0.8;
const K: f64 = 0.2;
const B: f64 = 3e-4;
/// The flat OCV of the DC test.
const E: f64 = 12.0;
const SOC0: f64 = 0.9;
/// The heat ports in the owner's order.
const HEAT: [&str; 10] = [
    "heat_converter_winding",
    "heat_input_high",
    "heat_input_low",
    "heat_output_high",
    "heat_output_low",
    "heat_battery_series",
    "heat_battery_polarization",
    "heat_lead",
    "heat_motor_winding",
    "heat_bearing",
];
/// A synchronous buck at half duty: input high, then input low, for 2.5 ms
/// each, output high on throughout.
const BUCK: &str = r#"{"intervals": [{"on": [true, false, true, false], "duration_s": 0.0025},
      {"on": [false, true, true, false], "duration_s": 0.0025}], "substeps": 3}"#;
/// Both high switches held on for one 1/256 s interval.
const HIGHS: &str = r#"{"intervals": [{"on": [true, false, true, false], "duration_s": 0.00390625}],
      "substeps": 2}"#;

fn close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * (1.0 + expected.abs()),
        "{what}: {actual} != {expected}"
    );
}

/// The drive's report coordinate `name`.
fn probe(scene: &fs_cli::SystemScene, name: &str) -> f64 {
    scene
        .probe(&format!("d.{name}"))
        .unwrap()
        .unwrap_or_else(|| panic!("d.{name} is null"))
}

fn run(text: &str, duration_s: f64) -> fs_cli::SystemScene {
    let mut scene = fs_cli::SystemScene::from_text(text, "drive.json").unwrap();
    scene_context(|cx, _| scene.run(cx).unwrap());
    close(scene.time_s(), duration_s, "time");
    assert!(
        scene.ledger().energy_residual_j().abs() < 1e-9,
        "ledger residual {} J",
        scene.ledger().energy_residual_j()
    );
    scene
}

/// The scene: the drive `d` on the OCV `curve` with `program`, `start` and
/// the `extra` fields, run for `duration_s` in steps of `step_s`, with no
/// torque on either shaft and every heat port held at 300 K. A switch that
/// turns on across a charged node carries about 120 A at first, so the
/// channels' limit is 1000 A.
fn scene(
    duration_s: f64,
    step_s: f64,
    curve: &str,
    program: &str,
    start: &str,
    extra: &str,
) -> String {
    let switches = ON
        .iter()
        .zip(OFF)
        .map(|(on, off)| {
            format!(
                r#"{{"on_resistance_ohm": {on}, "off_resistance_ohm": {off},
          "diode_forward_v": 0.7, "diode_resistance_ohm": 0.2, "max_reverse_v": 100,
          "max_channel_a": 1000, "max_diode_a": 1000}}"#
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    let heat = HEAT
        .map(|port| format!(r#"{{"prescribe": "d.{port}", "temperature_k": 300}}"#))
        .join(",\n    ");
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": {duration_s}, "step_s": {step_s},
  "elements": {{
    "d": {{"kind": "motor_drive",
      "battery": {{"capacity_c": {CAPACITY_C}, "soc_voltage": {curve},
        "series_resistance_ohm": {RS}, "polarization_resistance_ohm": {RP},
        "polarization_capacitance_f": {CP}, "bus_capacitance_f": {CB}}},
      "lead_resistance_ohm": {RLEAD},
      "converter": {{"inductance_h": {L}, "bus_capacitance_f": [{c0}, {c1}],
        "winding_resistance_ohm": {RW}, "switch": [{switches}]}},
      "motor": {{"inductance_h": {LM}, "winding_resistance_ohm": {RM},
        "rotor_inertia_kg_m2": {JR}, "stator_inertia_kg_m2": {JS},
        "coupling_nm_per_a": {K}, "bearing_drag_nm_s_per_rad": {B}}},
      "program": {program},
      "start": {start},
      "energy_tolerance_j": 1e-9{extra}}}
  }},
  "connections": [
    {{"prescribe": "d.rotor", "torque_n_m": 0}},
    {{"prescribe": "d.stator", "torque_n_m": 0}},
    {heat}
  ]
}}"#,
        c0 = BUS_F[0],
        c1 = BUS_F[1],
    )
}

/// The eleven-state lane's switch storage, as an extra field.
fn storage() -> String {
    format!(
        r#", "output_capacitance_f": [{}, {}, {}, {}]"#,
        COSS[0], COSS[1], COSS[2], COSS[3]
    )
}

/// The switching scene of 140 steps of 5 ms on a sloped curve, starting
/// with the rotor at 30 rad/s and the stator at 2 rad/s.
fn buck() -> String {
    let start = r#"{"soc": 0.8, "polarization_v": 0.05, "bus_v": 12.2,
      "converter_current_a": 1.0, "converter_bus_v": [12, 7], "node_v": [12, 7],
      "motor_current_a": 0.8, "shaft_rad_s": [30, 2]}"#;
    scene(
        0.7,
        0.005,
        "[[0, 11], [0.5, 12], [1, 12.6]]",
        BUCK,
        start,
        &storage(),
    )
}

#[test]
fn free_shafts_keep_their_momentum_through_a_switching_program() {
    // No torque acts on either shaft: the motor's torque and the bearing's
    // act between them, so J_r w_r + J_s w_s keeps its start. 0.7 s is 140
    // periods of 5 ms whose summed time misses the last one by more than
    // the owner's 64 rounding units: the last step is still the period.
    let scene = run(&buck(), 0.7);
    let momentum = [
        probe(&scene, "rotor_momentum_n_m_s"),
        probe(&scene, "stator_momentum_n_m_s"),
    ];
    let start = JR * 30.0 + JS * 2.0;
    assert!(
        (momentum[0] + momentum[1] - start).abs() <= 1e-12 * start,
        "momentum {momentum:?} from {start}"
    );
    // The drive turned the shafts apart.
    assert!((momentum[0] - JR * 30.0).abs() > 1e-3, "{momentum:?}");
    close(probe(&scene, "rotor_rad_s"), momentum[0] / JR, "rotor");
    close(probe(&scene, "stator_rad_s"), momentum[1] / JS, "stator");
    close(
        probe(&scene, "relative_rad_s"),
        momentum[0] / JR - momentum[1] / JS,
        "relative",
    );
    let report = scene.render_json();
    assert!(
        report.contains("a motor_drive is fs-phs's battery-fed commutated drive"),
        "{report}"
    );
    // The controls the run used are reported as given: here fs-system's
    // test values.
    for (control, value) in [
        ("solve_relative", 1e-13),
        ("solve_absolute", 1e-12),
        ("solve_max_iterations", 40.0),
        ("backtracking_halvings", 10.0),
        ("interval_absolute_j", 1e-10),
        ("interval_relative", 1e-12),
        ("cycle_max_intervals", 64.0),
        ("cycle_absolute_j", 1e-9),
        ("cycle_relative", 1e-11),
        ("energy_tolerance_j", 1e-9),
        ("period_s", 0.005),
    ] {
        assert_eq!(
            probe(&scene, control).to_bits(),
            f64::to_bits(value),
            "{control}"
        );
    }
}

/// The DC operating point with both high switches on, by this reading of
/// the owner's topology: the OCV behind the series and polarization
/// resistances (the capacitors carry no current), the lead to the input
/// bus, the input high channel to node a, the converter's winding from a to
/// node b, the output high channel from b to the output bus and the motor.
/// Each low switch is off and leaks from its node to ground through its off
/// resistance. The free shafts settle where the motor's torque k i meets the
/// bearing's b w, so the back EMF k w makes the motor a resistance
/// R_m + k^2 / b.
struct Dc {
    motor_a: f64,
    output_v: f64,
    b_v: f64,
    converter_a: f64,
    a_v: f64,
    battery_a: f64,
    input_v: f64,
}

impl Dc {
    fn new() -> Self {
        // Every voltage and current is proportional to the motor's current:
        // at 1 A, then scaled so the battery's loop closes on E.
        let output_v = RM + K * K / B;
        let b_v = output_v + ON[2];
        let converter_a = b_v / OFF[3] + 1.0;
        let a_v = b_v + RW * converter_a;
        let battery_a = a_v / OFF[1] + converter_a;
        let input_v = a_v + ON[0] * battery_a;
        let scale = E / (input_v + (RS + RP + RLEAD) * battery_a);
        Self {
            motor_a: scale,
            output_v: scale * output_v,
            b_v: scale * b_v,
            converter_a: scale * converter_a,
            a_v: scale * a_v,
            battery_a: scale * battery_a,
            input_v: scale * input_v,
        }
    }

    fn relative_rad_s(&self) -> f64 {
        K * self.motor_a / B
    }

    /// Each heat port's rate, in the owner's order.
    fn heat_w(&self) -> [f64; 10] {
        let (i_b, i_l, i_m) = (self.battery_a, self.converter_a, self.motor_a);
        [
            RW * i_l * i_l,
            ON[0] * i_b * i_b,
            self.a_v * self.a_v / OFF[1],
            ON[2] * i_m * i_m,
            self.b_v * self.b_v / OFF[3],
            RS * i_b * i_b,
            RP * i_b * i_b,
            RLEAD * i_b * i_b,
            RM * i_m * i_m,
            B * self.relative_rad_s().powi(2),
        ]
    }
}

#[test]
fn both_lanes_hold_the_dc_operating_point_of_the_resistive_network() {
    let dc = Dc::new();
    let i_b = dc.battery_a;
    let (period_s, duration_s) = (1.0 / 256.0, 0.0625);
    for eleven in [false, true] {
        let node_v = if eleven {
            format!(r#", "node_v": [{}, {}]"#, dc.a_v, dc.b_v)
        } else {
            String::new()
        };
        let start = format!(
            r#"{{"soc": {SOC0}, "polarization_v": {pv}, "bus_v": {bv},
      "converter_current_a": {il}, "converter_bus_v": [{v1}, {v2}]{node_v},
      "motor_current_a": {im}, "shaft_rad_s": [{w}, 0]}}"#,
            pv = RP * i_b,
            bv = E - (RS + RP) * i_b,
            il = dc.converter_a,
            v1 = dc.input_v,
            v2 = dc.output_v,
            im = dc.motor_a,
            w = dc.relative_rad_s(),
        );
        let extra = if eleven { storage() } else { String::new() };
        let text = scene(
            duration_s,
            period_s,
            &format!("[[0, {E}], [1, {E}]]"),
            HIGHS,
            &start,
            &extra,
        );
        let scene = run(&text, duration_s);
        let lane = if eleven { "eleven" } else { "nine" };
        let drawn = i_b * duration_s;
        for (name, expected) in [
            ("remaining_charge_c", SOC0 * CAPACITY_C - drawn),
            ("soc", SOC0 - drawn / CAPACITY_C),
            ("polarization_charge_c", CP * RP * i_b),
            ("polarization_v", RP * i_b),
            ("battery_bus_v", E - (RS + RP) * i_b),
            ("converter_flux_wb", L * dc.converter_a),
            ("converter_current_a", dc.converter_a),
            ("motor_flux_wb", LM * dc.motor_a),
            ("motor_current_a", dc.motor_a),
            ("rotor_momentum_n_m_s", JR * dc.relative_rad_s()),
            ("stator_momentum_n_m_s", 0.0),
            ("rotor_rad_s", dc.relative_rad_s()),
            ("relative_rad_s", dc.relative_rad_s()),
        ] {
            close(probe(&scene, name), expected, &format!("{lane}: {name}"));
        }
        if eleven {
            assert_eq!(scene.probe("d.input_bus_v").ok().flatten(), None);
            // fs-phs's electrode charges: a bus holds its capacitor's and its
            // high switch's, a node its low switch's less its high switch's.
            let (bus, node) = ([dc.input_v, dc.output_v], [dc.a_v, dc.b_v]);
            for (leg, side) in ["input", "output"].into_iter().enumerate() {
                let high = COSS[2 * leg] * (bus[leg] - node[leg]);
                let low = COSS[2 * leg + 1] * node[leg];
                let name = format!("{side}_bus_charge_c");
                close(probe(&scene, &name), BUS_F[leg] * bus[leg] + high, &name);
                let name = format!("{side}_node_charge_c");
                close(probe(&scene, &name), low - high, &name);
            }
        } else {
            close(probe(&scene, "input_bus_v"), dc.input_v, "input bus");
            close(probe(&scene, "output_bus_v"), dc.output_v, "output bus");
            close(
                probe(&scene, "input_bus_charge_c"),
                BUS_F[0] * dc.input_v,
                "input bus charge",
            );
        }
        for (port, rate) in HEAT.into_iter().zip(dc.heat_w()) {
            let name = format!("{port}_w");
            close(probe(&scene, &name), rate, &format!("{lane}: {port}"));
        }
        // The heat rates account for the battery's whole power.
        close(dc.heat_w().iter().sum(), E * i_b, "power");
    }
}

#[test]
fn a_drive_reports_null_before_a_step_its_authored_controls_and_refuses_bad_scenes() {
    let good = buck();
    let admitted = fs_cli::SystemScene::from_text(&good, "drive.json").unwrap();
    assert_eq!(admitted.probe("d.heat_lead_w").unwrap(), None);
    close(
        probe(&admitted, "rotor_momentum_n_m_s"),
        JR * 30.0,
        "start rotor",
    );
    close(probe(&admitted, "stator_rad_s"), 2.0, "start stator");
    close(probe(&admitted, "soc"), 0.8, "start soc");
    close(probe(&admitted, "converter_current_a"), 1.0, "current");
    let authored = good.replace(
        r#""energy_tolerance_j": 1e-9"#,
        r#""energy_tolerance_j": 2e-9,
      "solve": {"relative": 2e-13, "absolute": 3e-12, "max_iterations": 30},
      "backtracking_halvings": 8, "interval_tolerance": {"absolute_j": 2e-10, "relative": 3e-12},
      "cycle": {"max_intervals": 32, "absolute_j": 2e-9, "relative": 3e-11}"#,
    );
    let admitted = fs_cli::SystemScene::from_text(&authored, "drive.json").unwrap();
    for (control, value) in [
        ("solve_relative", 2e-13),
        ("solve_absolute", 3e-12),
        ("solve_max_iterations", 30.0),
        ("backtracking_halvings", 8.0),
        ("interval_absolute_j", 2e-10),
        ("interval_relative", 3e-12),
        ("cycle_max_intervals", 32.0),
        ("cycle_absolute_j", 2e-9),
        ("cycle_relative", 3e-11),
        ("energy_tolerance_j", 2e-9),
    ] {
        assert_eq!(
            probe(&admitted, control).to_bits(),
            f64::to_bits(value),
            "{control}"
        );
    }
    let period = "step_s must be the period and duration_s a whole number of steps";
    for (bad, expected) in [
        (
            good.replace(r#""step_s": 0.005"#, r#""step_s": 0.004"#),
            period,
        ),
        (
            good.replace(r#""duration_s": 0.7"#, r#""duration_s": 0.7025"#),
            period,
        ),
        (
            good.replace(
                r#""switch": [{"#,
                r#""switch": [{"on_resistance_ohm": 1}, {"#,
            ),
            "one switch, or a list of four",
        ),
        (
            good.replace(
                r#""on": [true, false, true, false]"#,
                r#""on": [true, false, true]"#,
            ),
            "four booleans",
        ),
        (
            good.replace(r#""soc": 0.8"#, r#""soc": 1.5"#),
            "elements.d.start: ",
        ),
        (
            good.replace(&storage(), ""),
            "elements.d.start: unsupported field \"node_v\"",
        ),
        (
            good.replace(
                r#""step_s": 0.005,"#,
                r#""step_s": 0.005, "adaptive": {"order": 2, "tolerances": {}},"#,
            ),
            "is a motor_drive, whose state has no error account here",
        ),
    ] {
        assert_ne!(bad, good, "the replacement must change the scene");
        let error = fs_cli::SystemScene::from_text(&bad, "drive.json")
            .err()
            .unwrap_or_else(|| panic!("admitted: {expected}"))
            .to_string();
        assert!(error.contains(expected), "{expected}: {error}");
    }
}
