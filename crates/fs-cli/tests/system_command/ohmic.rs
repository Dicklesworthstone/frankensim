//! Card-backed resistors: nickel_wire_heater.json, a 1 m Nickel 200 wire of
//! 2e-7 m^2 cross-section on a fixed 1 V, its Joule heat into a lumped body of
//! the wire's own mass on the same card. The resistance is not typed: the
//! wire names the scene material `nickel200`, whose card is the Special
//! Metals SMC-061 table (resistivity 9.6e-8 ohm m at 293.15 K, 1.30e-7 at
//! 373.15 K, 1.85e-7 at 473.15 K) and whose specific heat is one constant,
//! 456 J/(kg K).
//!
//! Claims, written before the first run:
//! - The pack's resistivity curve is the seed source's printed table, knot
//!   for knot.
//! - After each of 120 steps of 0.5 s the body's temperature is within
//!   1e-8 K of a hand-written backward step: the end temperature bisected so
//!   that `m c (T1 - T0) = dt V^2 / R(T_lm)`, with `R` the printed table's
//!   linear interpolation times `L / A` and `T_lm` the logarithmic mean of
//!   the step's two end temperatures, which the body's heat port publishes.
//!   The run crosses the 373.15 K knot. The reported resistance is the
//!   printed table's at the reported temperature, within 1e-12 relative,
//!   and the current is the voltage over it; the ledger's residual is
//!   within 1e-9 J.
//! - The series form (`resistor`, a current loop with an inductance): the
//!   inductance's lost energy is the body's gain, within 1e-9 J, and the
//!   current decays.
//! - The binary runs the example, binds the wire to its material, lists the
//!   resistivity receipt, and the no-claim sentence names the card-backed
//!   resistor.
//! - Refusals: a resistance beside a material, a material without a
//!   resistivity, an unknown material, a geometry that is not positive, a
//!   span past the card's resistivity, and 40 V, which would carry the wire
//!   off the declared span in the first step; that scene does not move.

use super::*;
use fs_cli::SystemScene;
use fs_matdb::{NormalizedMaterialCardPack, PropertyValue};

/// The example, its pack and the seed source the pack is compiled from.
const NICKEL_EXAMPLE: &str = "nickel_wire_heater.json";
const NICKEL_PACK: &str = "nickel200/nickel-200-annealed.fsmcdpk";
const NICKEL_SOURCE: &str =
    "data/matdb/seed-v1/nickel-200-annealed-special-metals-2006/properties.tsv";
/// The example's wire [m, m^2], the body's mass [kg] (the card's 8890 kg/m^3
/// times the wire's volume, as the scene writes it), the card's specific
/// heat [J/(kg K)], the voltage [V], the step [s] and the steps.
const LENGTH_M: f64 = 1.0;
const AREA_M2: f64 = 2.0e-7;
const MASS_KG: f64 = 1.778e-3;
const SPECIFIC_HEAT_J_KG_K: f64 = 456.0;
const VOLTAGE_V: f64 = 1.0;
const STEP_S: f64 = 0.5;
const STEPS: u64 = 120;
const START_K: f64 = 293.15;

fn probe(scene: &SystemScene, path: &str) -> f64 {
    scene
        .probe(path)
        .unwrap()
        .unwrap_or_else(|| panic!("{path} is not reported"))
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64, what: &str) {
    let allowed = absolute.max(relative * expected.abs());
    assert!(
        (actual - expected).abs() <= allowed,
        "{what}: {actual} against {expected} (allowed {allowed:e})"
    );
}

/// The resistivity table the seed source prints, from its `curve` row:
/// [K, ohm m].
fn printed_resistivity() -> Vec<(f64, f64)> {
    let path = format!("{}/../../{NICKEL_SOURCE}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).unwrap();
    let row = text
        .lines()
        .find(|line| line.starts_with("curve\tresistivity\t"))
        .expect("the source has the resistivity curve");
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(fields[3], "electrical_resistivity");
    assert_eq!(
        (fields[4], fields[5], fields[6], fields[8]),
        ("temperature", "K", "Ohm*m", "linear")
    );
    fields[7]
        .split(',')
        .map(|knot| {
            let (kelvin, resistivity) = knot.split_once(':').unwrap();
            (kelvin.parse().unwrap(), resistivity.parse().unwrap())
        })
        .collect()
}

/// `rho(T) L / A` from the printed table, linear inside each segment.
fn printed_resistance(printed: &[(f64, f64)], kelvin: f64) -> f64 {
    let segment = printed
        .iter()
        .rposition(|&(knot, _)| knot <= kelvin)
        .expect("inside the table")
        .min(printed.len() - 2);
    let ((t0, r0), (t1, r1)) = (printed[segment], printed[segment + 1]);
    (r0 + (r1 - r0) * (kelvin - t0) / (t1 - t0)) * LENGTH_M / AREA_M2
}

/// The logarithmic mean of a step's two end temperatures.
fn log_mean(start: f64, end: f64) -> f64 {
    let increment = end - start;
    if increment == 0.0 {
        start
    } else {
        increment / (increment / start).ln_1p()
    }
}

/// One backward step of the body heated by `V^2 / R` at the logarithmic
/// mean: the end temperature by bisection.
fn hand_step(printed: &[(f64, f64)], kelvin: f64) -> f64 {
    let capacity = MASS_KG * SPECIFIC_HEAT_J_KG_K;
    let rise = |next: f64| {
        STEP_S * VOLTAGE_V * VOLTAGE_V
            / printed_resistance(printed, log_mean(kelvin, next))
            / capacity
    };
    // The power falls as the wire warms, so the rise at the start's power
    // bounds the end temperature from above.
    let (mut lower, mut upper) = (kelvin, kelvin + rise(kelvin));
    for _ in 0..80 {
        let trial = 0.5 * lower + 0.5 * upper;
        if trial - kelvin < rise(trial) {
            lower = trial;
        } else {
            upper = trial;
        }
    }
    0.5 * lower + 0.5 * upper
}

#[test]
fn a_nickel_wire_on_a_fixed_voltage_is_the_hand_step_of_the_printed_resistivity() {
    let printed = printed_resistivity();
    let pack =
        NormalizedMaterialCardPack::from_bytes(&std::fs::read(example(NICKEL_PACK)).unwrap())
            .unwrap();
    let claims = pack.card().claims().claims_for("electrical-resistivity");
    assert_eq!(claims.len(), 1);
    let PropertyValue::Curve {
        knots: compiled, ..
    } = &claims[0].1.value
    else {
        panic!("the resistivity is not a curve");
    };
    assert_eq!(compiled, &printed, "the pack is not the source's table");
    let text = std::fs::read_to_string(example(NICKEL_EXAMPLE)).unwrap();
    scene_context(|cx, _| {
        let mut scene = SystemScene::from_text(&text, example(NICKEL_EXAMPLE)).unwrap();
        let (mut kelvin, mut mean) = (START_K, START_K);
        while scene.advance(f64::MAX, cx).unwrap() {
            let next = hand_step(&printed, kelvin);
            mean = log_mean(kelvin, next);
            kelvin = next;
            close(
                probe(&scene, "body.temperature_k"),
                kelvin,
                0.0,
                1e-8,
                &format!("step {}", scene.steps()),
            );
        }
        assert_eq!(scene.steps(), STEPS);
        let read = probe(&scene, "wire.temperature_k");
        let resistance = probe(&scene, "wire.resistance_ohm");
        let current = probe(&scene, "wire.current_a");
        close(
            read,
            mean,
            0.0,
            1e-8,
            "the last step's heat-port temperature",
        );
        close(
            resistance,
            printed_resistance(&printed, read),
            1e-12,
            0.0,
            "the reported resistance",
        );
        close(current, VOLTAGE_V / resistance, 1e-12, 0.0, "the current");
        let ledger = scene.ledger();
        let (stored, residual) = (ledger.stored_change_j(), ledger.energy_residual_j());
        eprintln!(
            "{} steps; body {kelvin} K, wire read {read} K, {resistance} ohm, {current} A; stored change {stored} J, ledger residual {residual:e} J",
            scene.steps()
        );
        close(
            stored,
            MASS_KG * SPECIFIC_HEAT_J_KG_K * (kelvin - START_K),
            0.0,
            1e-7,
            "the stored change",
        );
        assert!(residual.abs() <= 1e-9, "ledger residual {residual:e} J");
        // The run crosses the table's 373.15 K knot and stays inside the
        // declared span.
        assert!(380.0 < kelvin && kelvin < 473.15, "{kelvin} K");
    });
}

/// The example's materials section, for scenes written around it.
fn materials() -> String {
    let text = std::fs::read_to_string(example(NICKEL_EXAMPLE)).unwrap();
    let start = text.find("  \"materials\": {").expect("the part begins");
    let end = start + text[start..].find("  \"elements\"").expect("the part ends");
    text[start..end].to_owned()
}

#[test]
fn in_a_current_loop_the_inductances_lost_energy_is_the_wires_body_gain() {
    let inductance = 0.5;
    let text = format!(
        "{{\n  \"schema\": \"frankensim.system.v0\",\n  \"duration_s\": 2.0,\n  \"step_s\": 0.05,\n{}  \"elements\": {{\n    \"coil\": {{\"kind\": \"inductor\", \"inductance_h\": {inductance}, \"current_a\": 4.0}},\n    \"wire\": {{\"kind\": \"resistor\", \"material\": \"nickel200\", \"length_m\": 1.0, \"area_m2\": 2.0e-7}},\n    \"body\": {{\"kind\": \"enthalpy_body\", \"material\": \"nickel200\", \"mass_kg\": 1.778e-3, \"specific_enthalpy_j_kg\": 0.0}}\n  }},\n  \"connections\": [\n    {{\"join\": [\"coil.p\", \"wire.a\"]}},\n    {{\"prescribe\": \"wire.b\", \"current_a\": 0}},\n    {{\"join\": [\"body.heat\", \"wire.heat\"]}}\n  ]\n}}\n",
        materials()
    );
    scene_context(|cx, _| {
        let mut scene = SystemScene::from_text(&text, example(NICKEL_EXAMPLE)).unwrap();
        scene.run(cx).unwrap();
        let current = probe(&scene, "coil.current_a");
        let lost = 0.5 * inductance * (4.0 * 4.0 - current * current);
        let gained = MASS_KG * probe(&scene, "body.specific_enthalpy_j_kg");
        eprintln!(
            "loop: {} steps, current 4 A -> {current} A, body {} K, wire {} ohm, lost {lost} J, gained {gained} J",
            scene.steps(),
            probe(&scene, "body.temperature_k"),
            probe(&scene, "wire.resistance_ohm"),
        );
        assert!(
            current > 0.0 && current < 1.0,
            "the current decays: {current}"
        );
        close(gained, lost, 0.0, 1e-9, "the body's gain");
    });
}

#[test]
fn the_wire_example_runs_through_the_binary_and_reports_its_card() {
    let path = example(NICKEL_EXAMPLE);
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(value(&json, &element("wire"), "temperature_k") > 380.0);
    for member in [
        "\"wire\":[\"nickel200\"]",
        "\"body\":[\"nickel200\"]",
        "\"property\":\"electrical-resistivity\"",
    ] {
        assert!(json.contains(member), "{member} not in {json}");
    }
    let (code, text, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(text.contains("wire <- nickel200"), "{text}");
    for said in [
        "lumped elements with constant parameters, card-backed resistances or",
        "source resolution is limited to declared pure-phase thermal laws and resistors' \
         resistivity; no contact, friction or step-size error control",
        "a card-backed resistor is fs-material's Ohmic node of fixed length and \
         cross-section",
    ] {
        assert!(text.contains(said), "{said} not in {text}");
    }
    assert!(!text.contains("no material-card resolution"), "{text}");
}

#[test]
fn card_backed_resistors_refuse_what_the_card_cannot_answer() {
    let base = std::fs::read_to_string(example(NICKEL_EXAMPLE)).unwrap();
    let wire = "\"material\": \"nickel200\", \"length_m\": 1.0,";
    let cases = [
        (
            "resistance beside a material",
            wire,
            "\"material\": \"nickel200\", \"resistance_ohm\": 0.5, \"length_m\": 1.0,",
            "elements.wire: unsupported field \"resistance_ohm\"",
        ),
        (
            "no resistivity",
            ",\n      \"resistivity_property\": \"electrical-resistivity\"",
            "",
            "names no resistivity_property",
        ),
        (
            "unknown material",
            wire,
            "\"material\": \"copper\", \"length_m\": 1.0,",
            "elements.wire.material: unknown material \"copper\"",
        ),
        (
            "no area",
            "\"area_m2\": 2.0e-7",
            "\"area_m2\": 0.0",
            "elements.wire: area_m2 must be positive and finite",
        ),
        (
            "a span past the resistivity",
            "\"temperature_max_k\": 473.15,",
            "\"temperature_max_k\": 1500.0,",
            "elements.wire.material: conductor material refusal: Support { property: \
             \"electrical-resistivity\"",
        ),
    ];
    for (case, from, to, reason) in cases {
        assert!(base.contains(from), "{case}: test fixture changed");
        let text = base.replacen(from, to, 1);
        let Err(error) = SystemScene::from_text(&text, example(NICKEL_EXAMPLE)) else {
            panic!("{case}: admitted");
        };
        eprintln!("{case}: exit {}: {}", error.exit_code, error.message);
        assert!(error.message.contains(reason), "{case}: {}", error.message);
    }
    // 40 V puts more than 1700 W into 0.81 J/K even at the span's highest
    // resistance: the first step would carry the wire past 473.15 K, it is
    // refused, and the scene does not move.
    let hot = base.replacen("\"voltage_v\": 1.0", "\"voltage_v\": 40.0", 1);
    assert_ne!(hot, base);
    let mut scene = SystemScene::from_text(&hot, example(NICKEL_EXAMPLE)).unwrap();
    let before = scene.render_json();
    scene_context(|cx, _| {
        let error = scene.advance(f64::MAX, cx).unwrap_err();
        eprintln!("40 V: exit {}: {}", error.exit_code, error.message);
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert_eq!(scene.render_json(), before);
        assert_eq!(scene.steps(), 0);
    });
}
