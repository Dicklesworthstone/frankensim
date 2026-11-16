//! Card-backed pipes: shell_s2_oil_line.json, a 1 m line of 0.5 mm radius
//! on a fixed 1e5 Pa, carrying Shell Heat Transfer Oil S2 while a 200 W
//! heater warms 0.1 kg of the same oil, the loop's well-mixed inventory, for
//! 90 s. The line's resistance is not typed: it names the scene material
//! `s2oil`, whose card is Shell's May 2011 data sheet (kinematic viscosity
//! 2.5e-5 m^2/s at 313.15 K, 4.7e-6 at 373.15 K and 1.1e-6 at 473.15 K, density
//! 850, 811 and 746 kg/m^3 there), and fs-material's laminar pipe reads its
//! viscosity `rho nu` again at the temperature the line's heat port receives.
//!
//! Claims, written before the first run:
//! - The pack's kinematic-viscosity and density curves are the seed source's
//!   printed tables, knot for knot.
//! - After each of 180 steps the reported flow is `pi r^4 dp / (8 rho nu L)`
//!   with `rho` and `nu` the printed tables' linear interpolation at the
//!   temperature the line reports, within 1e-12 relative; the reported
//!   viscosity is `rho nu` there and the heat is the flow times the drop,
//!   both within 1e-12. That temperature lies between the oil's temperatures
//!   at the step's two ends, and the flow rises at every step. The run
//!   crosses the tables' 373.15 K knot and stays inside the declared span.
//! - The heater's work is 200 W times 90 s within 1e-9 J, the oil's stored
//!   change is that work plus the line's heat summed over the steps within
//!   1e-6 J, and the ledger's residual is within 1e-6 J. The last step's
//!   Reynolds number `2 rho Q / (pi r mu)`, by hand, is below 500: the
//!   example is laminar, which the pipe assumes and does not check.
//! - The binary runs the example, binds the line and the oil to their
//!   material, lists the viscosity and density receipts, and the no-claim
//!   sentence names the card-backed pipe.
//! - Refusals: a resistance beside a material, a material without a
//!   viscosity, a material with both a dynamic and a kinematic one, a
//!   dynamic viscosity the card does not have, an unknown material, a radius
//!   that is not positive, and a heater that would carry the oil off the
//!   declared span in the first step; that scene does not move.
//!
//! Corrected after the first two runs, in the test, not the code: the test
//! first read the heater's `energy_in_j` as the work it gave, but the scene's
//! account of a source is the energy into it (-18000 J here), as the
//! heated-body tests read it. It then took the ledger's total stored change
//! for the oil's, but fs-system's heater books its account as stored energy,
//! so that total is the oil's gain less the heater's work: the line's heat
//! (3.5726 J). The oil's gain is now its mass times its specific enthalpy,
//! from the material's datum, and the ledger's total is checked against the
//! line's heat.

use super::*;
use fs_cli::SystemScene;
use fs_matdb::{NormalizedMaterialCardPack, PropertyValue};

/// The example, its pack and the seed source the pack is compiled from.
const OIL_EXAMPLE: &str = "shell_s2_oil_line.json";
const OIL_PACK: &str = "shell_s2/shell-heat-transfer-oil-s2-2011.fsmcdpk";
const OIL_SOURCE: &str = "data/matdb/seed-v1/shell-heat-transfer-oil-s2-2011/properties.tsv";
/// The example's line [m, m], drop [Pa], heater [W], step [s] and steps.
const RADIUS_M: f64 = 5.0e-4;
const LENGTH_M: f64 = 1.0;
const DROP_PA: f64 = 1.0e5;
const HEATER_W: f64 = 200.0;
const OIL_KG: f64 = 0.1;
const STEP_S: f64 = 0.5;
const STEPS: u64 = 180;

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

/// One table the seed source prints, from its `curve` row: [K, SI].
fn printed(row: &str, property: &str, unit: &str) -> Vec<(f64, f64)> {
    let path = format!("{}/../../{OIL_SOURCE}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).unwrap();
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("curve\t{row}\t")))
        .unwrap_or_else(|| panic!("the source has the {row} curve"));
    let fields: Vec<&str> = line.split('\t').collect();
    assert_eq!(fields[3], property);
    assert_eq!(
        (fields[4], fields[5], fields[6], fields[8]),
        ("temperature", "K", unit, "linear")
    );
    fields[7]
        .split(',')
        .map(|knot| {
            let (kelvin, value) = knot.split_once(':').unwrap();
            (kelvin.parse().unwrap(), value.parse().unwrap())
        })
        .collect()
}

fn printed_viscosity() -> Vec<(f64, f64)> {
    printed("nu-0-200c", "kinematic_viscosity", "m2/s")
}

fn printed_density() -> Vec<(f64, f64)> {
    printed("rho-0-200c", "density", "kg/m3")
}

/// A printed table's linear interpolation at `kelvin`.
fn interpolate(table: &[(f64, f64)], kelvin: f64) -> f64 {
    let segment = table
        .iter()
        .rposition(|&(knot, _)| knot <= kelvin)
        .expect("inside the table")
        .min(table.len() - 2);
    let ((t0, v0), (t1, v1)) = (table[segment], table[segment + 1]);
    v0 + (v1 - v0) * (kelvin - t0) / (t1 - t0)
}

#[test]
fn an_oil_line_follows_the_printed_viscosity_of_the_oil_its_heater_warms() {
    let (nu, rho) = (printed_viscosity(), printed_density());
    let pack =
        NormalizedMaterialCardPack::from_bytes(&std::fs::read(example(OIL_PACK)).unwrap()).unwrap();
    for (property, table) in [("kinematic_viscosity", &nu), ("density", &rho)] {
        let claims = pack.card().claims().claims_for(property);
        assert_eq!(claims.len(), 1, "{property}");
        let PropertyValue::Curve {
            knots: compiled, ..
        } = &claims[0].1.value
        else {
            panic!("the {property} is not a curve");
        };
        assert_eq!(
            compiled, table,
            "the pack is not the source's {property} table"
        );
    }
    let text = std::fs::read_to_string(example(OIL_EXAMPLE)).unwrap();
    scene_context(|cx, _| {
        let mut scene = SystemScene::from_text(&text, example(OIL_EXAMPLE)).unwrap();
        let (mut previous_oil, mut previous_flow, mut heat_j) = (313.15, 0.0, 0.0);
        while scene.advance(f64::MAX, cx).unwrap() {
            let step = scene.steps();
            let kelvin = probe(&scene, "line.temperature_k");
            let oil = probe(&scene, "oil.temperature_k");
            assert!(
                previous_oil <= kelvin && kelvin <= oil,
                "step {step}: the line read {kelvin} K, the oil went {previous_oil} K -> {oil} K"
            );
            let viscosity = interpolate(&rho, kelvin) * interpolate(&nu, kelvin);
            let flow =
                core::f64::consts::PI * RADIUS_M.powi(4) * DROP_PA / (8.0 * viscosity * LENGTH_M);
            let reported = probe(&scene, "line.flow_m3_s");
            close(
                reported,
                flow,
                1e-12,
                0.0,
                &format!("step {step}: the flow"),
            );
            close(
                probe(&scene, "line.viscosity_pa_s"),
                viscosity,
                1e-12,
                0.0,
                &format!("step {step}: the viscosity"),
            );
            let heat = probe(&scene, "line.heat_w");
            close(
                heat,
                reported * DROP_PA,
                1e-12,
                0.0,
                &format!("step {step}: the heat"),
            );
            assert!(reported > previous_flow, "step {step}: the flow rises");
            heat_j += STEP_S * heat;
            (previous_oil, previous_flow) = (oil, reported);
        }
        assert_eq!(scene.steps(), STEPS);
        // A source's account is the energy into it: the heater's is the
        // negative of what it gave.
        let work = -probe(&scene, "heater.energy_in_j");
        close(
            work,
            HEATER_W * STEP_S * STEPS as f64,
            0.0,
            1e-9,
            "the heater's work",
        );
        // The oil starts at the material's enthalpy datum.
        let gained = OIL_KG * probe(&scene, "oil.specific_enthalpy_j_kg");
        close(gained, work + heat_j, 0.0, 1e-6, "the oil's gain");
        // The ledger books the heater's account, the negative of its work,
        // as stored energy beside the oil's gain.
        let ledger = scene.ledger();
        let (stored, residual) = (ledger.stored_change_j(), ledger.energy_residual_j());
        close(stored, heat_j, 0.0, 1e-6, "the ledger's stored change");
        assert!(residual.abs() <= 1e-6, "ledger residual {residual:e} J");
        let kelvin = probe(&scene, "line.temperature_k");
        let viscosity = probe(&scene, "line.viscosity_pa_s");
        let reynolds = 2.0 * interpolate(&rho, kelvin) * previous_flow
            / (core::f64::consts::PI * RADIUS_M * viscosity);
        eprintln!(
            "{} steps; oil {previous_oil} K, line read {kelvin} K, viscosity {viscosity} Pa s, flow {previous_flow} m3/s, Reynolds {reynolds}; heater {work} J, line heat {heat_j} J, stored change {stored} J, ledger residual {residual:e} J",
            scene.steps()
        );
        assert!(reynolds < 500.0, "Reynolds {reynolds}");
        // The run crosses the tables' 373.15 K knot and stays inside the
        // declared span.
        assert!(
            380.0 < previous_oil && previous_oil < 473.15,
            "{previous_oil} K"
        );
    });
}

#[test]
fn the_oil_line_example_runs_through_the_binary_and_reports_its_card() {
    let path = example(OIL_EXAMPLE);
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(value(&json, &element("line"), "temperature_k") > 380.0);
    for member in [
        "\"line\":[\"s2oil\"]",
        "\"oil\":[\"s2oil\"]",
        "\"property\":\"kinematic_viscosity\"",
        "\"property\":\"density\"",
    ] {
        assert!(json.contains(member), "{member} not in {json}");
    }
    let (code, text, stderr) = run_binary(&["system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert!(text.contains("line <- s2oil"), "{text}");
    for said in [
        "lumped elements with constant parameters, card-backed pipes or",
        "source resolution is limited to declared pure-phase thermal laws and pipes' \
         viscosity; no ",
        "a card-backed pipe is fs-material's laminar circular pipe of fixed radius and \
         length",
        "laminar by assumption with no Reynolds number checked",
    ] {
        assert!(text.contains(said), "{said} not in {text}");
    }
    assert!(!text.contains("no material-card resolution"), "{text}");
}

#[test]
fn card_backed_pipes_refuse_what_the_card_cannot_answer() {
    let base = std::fs::read_to_string(example(OIL_EXAMPLE)).unwrap();
    let line = "\"material\": \"s2oil\", \"radius_m\": 5.0e-4,";
    let kinematic = "\"kinematic_viscosity_property\": \"kinematic_viscosity\"";
    let cases = [
        (
            "resistance beside a material",
            line,
            "\"material\": \"s2oil\", \"resistance_pa_s_m3\": 1.0e9, \"radius_m\": 5.0e-4,",
            "elements.line: unsupported field \"resistance_pa_s_m3\"",
        ),
        (
            "no viscosity",
            ",\n      \"kinematic_viscosity_property\": \"kinematic_viscosity\"",
            "",
            "names no dynamic_viscosity_property or kinematic_viscosity_property",
        ),
        (
            "both viscosities",
            kinematic,
            "\"dynamic_viscosity_property\": \"dynamic_viscosity\", \
             \"kinematic_viscosity_property\": \"kinematic_viscosity\"",
            "materials.s2oil: a dynamic_viscosity_property or a \
             kinematic_viscosity_property, not both",
        ),
        (
            "a dynamic viscosity the card does not have",
            kinematic,
            "\"dynamic_viscosity_property\": \"dynamic_viscosity\"",
            "materials.s2oil.dynamic_viscosity_property: ",
        ),
        (
            "unknown material",
            line,
            "\"material\": \"water\", \"radius_m\": 5.0e-4,",
            "elements.line.material: unknown material \"water\"",
        ),
        (
            "no radius",
            "\"radius_m\": 5.0e-4",
            "\"radius_m\": 0.0",
            "elements.line: radius_m must be positive and finite",
        ),
    ];
    for (case, from, to, reason) in cases {
        assert!(base.contains(from), "{case}: test fixture changed");
        let text = base.replacen(from, to, 1);
        let Err(error) = SystemScene::from_text(&text, example(OIL_EXAMPLE)) else {
            panic!("{case}: admitted");
        };
        eprintln!("{case}: exit {}: {}", error.exit_code, error.message);
        assert!(error.message.contains(reason), "{case}: {}", error.message);
    }
    // 2 MW into 0.1 kg of oil would carry it thousands of kelvin past
    // 473.15 K in the first step: the step is refused, and the scene does
    // not move.
    let hot = base.replacen("\"power_w\": 200.0", "\"power_w\": 2.0e6", 1);
    assert_ne!(hot, base);
    let mut scene = SystemScene::from_text(&hot, example(OIL_EXAMPLE)).unwrap();
    let before = scene.render_json();
    scene_context(|cx, _| {
        let error = scene.advance(f64::MAX, cx).unwrap_err();
        eprintln!("2 MW: exit {}: {}", error.exit_code, error.message);
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert_eq!(scene.render_json(), before);
        assert_eq!(scene.steps(), 0);
    });
}
