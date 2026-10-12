//! Native material-coordinate adjoints checked by fresh, independently sourced
//! conductivity and contact cards. The phase fixture exits a latent plateau
//! only at its final fine endpoint, so its temperature gradient needs history.

use super::*;

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";
const DK: &str = "dtemperature_dconductivity_multiplier_k";
const DR: &str = "dtemperature_dresistance_multiplier_k";

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn rows<'a>(row: &'a J, key: &str) -> &'a [J] {
    row.get(key).unwrap().as_array().unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

struct ThermalRun {
    receipt: J,
    temperature: Vec<f64>,
}

impl ThermalRun {
    fn transient(&self) -> &J {
        self.receipt.get("transient").unwrap()
    }

    fn adjoint(&self) -> &J {
        let adjoint = self.transient().get("adjoint").unwrap();
        assert_eq!(adjoint.str_field("schema"), Some("fs-cli-native-enthalpy-adjoint-v1"));
        adjoint
    }

    fn controls(&self) -> &J {
        let adjoint = self.adjoint();
        let controls = adjoint.get("material_controls").unwrap();
        close(n(controls, "control_steps"), n(adjoint, "state_index"), 0.0);
        controls
    }

    fn maximum(&self) -> f64 {
        self.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }
}

fn run(
    spec: &ProjectSpec,
    cards: &CardPackSet,
    bytes: Option<&[u8]>,
    output: Option<&str>,
) -> ThermalRun {
    let mut spec = spec.clone();
    if let Some(name) = output {
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: name.into(), kind: "report".into(), region: None,
        });
    }
    let decoded = decode(&spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    if let Some(bytes) = bytes {
        import_fixture(&ledger, &spec, bytes.to_vec());
    } else {
        import_multi_region_contact(&ledger, &spec);
    }
    let (_, _, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    let field = J::parse(&field).unwrap();
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        temperature: rows(&field, "temperature").iter()
            .map(|value| value.as_f64().unwrap()).collect(),
    }
}

fn phase_fixture(scale: f64) -> (ProjectSpec, CardPackSet) {
    let cards = CardPackSet::admit(vec![raw_pack(
        CardPackKind::Material,
        "fixtures/enthalpy-conductivity-control.fsmcdpk",
        material_pack_bytes_with_property(
            "AA6061", &format!("enthalpy-conductivity-control-{scale}"), 280.0, 360.0,
            fs_matdb::PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_conduction::TEMPERATURE_DIMS,
                // K(T) has nonzero slope throughout the accepted trajectory.
                // Scale every ordinate; leave all phase-law parameters fixed.
                knots: [(280.0, 1.0), (298.0, 1.0), (302.0, 3.0), (360.0, 3.0)]
                    .into_iter().map(|(t, k)| (t, scale * k)).collect(),
                dims: CONDUCTIVITY_DIMS,
            },
            fs_matdb::InterpolationPolicy::LinearInside,
        ),
    )]).unwrap();
    let mut spec = transient_enthalpy_conductivity::project(&cards, Some((1.0, 4.0)));
    let time = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap();
    time.horizon.value = 6.0;
    time.max_steps = 36;
    assert_eq!(
        time.enthalpy.as_ref().unwrap().materials[0].material_card_identity.to_hex(),
        cards.materials()[0].card().to_hex(),
    );
    (spec, cards)
}

fn check_latent_exit(solved: &ThermalRun) {
    let fine = rows(solved.transient(), "fine");
    assert_eq!(fine.len(), 24);
    assert!(n(&fine[1], "mean_liquid_mass_fraction") > 0.0);
    let before = &fine[fine.len() - 2];
    close(n(before, "final_region_max_k"), 301.0, 0.0);
    assert!(n(before, "mean_liquid_mass_fraction") > 0.24);
    assert!(n(before, "mean_liquid_mass_fraction") < 0.25);
    close(n(fine.last().unwrap(), "mean_liquid_mass_fraction"), 0.25, 1e-12);
    assert_eq!(solved.temperature.len(), 4);
    assert_eq!(solved.temperature.iter().filter(|&&t| t > 301.0).count(), 1);
    assert!(solved.maximum() > 301.05 && solved.maximum() < 301.12);
    let law = rows(solved.transient(), "materials")[0]
        .get("phase_conductivity").unwrap();
    close(n(law, "solid_multiplier"), 1.0, 0.0);
    close(n(law, "liquid_multiplier"), 4.0, 0.0);
}

#[test]
fn g1_native_enthalpy_conductivity_multiplier_carries_nonlinear_phase_history() {
    let bytes = tetra_stl();
    let (spec, cards) = phase_fixture(1.0);
    let solved = run(&spec, &cards, Some(&bytes), Some(FINAL));
    check_latent_exit(&solved);
    let adjoint = solved.adjoint();
    close(n(adjoint, "value_k"), solved.maximum(), 0.0);
    close(n(adjoint, "state_index"), 24.0, 0.0);
    assert!(n(adjoint, "maximum_dual_relative_residual") < 1e-9);
    let controls = solved.controls();
    assert!(rows(controls, "contact_interfaces").is_empty());
    let regions = rows(controls, "conductivity_regions");
    assert_eq!(regions.len(), 1);
    let control = &regions[0];
    let binding = &spec.materials.as_ref().unwrap()[0];
    assert_eq!(control.str_field("region"), Some(binding.region.as_str()));
    assert_eq!(control.str_field("material_card"), Some(binding.card.as_str()));
    assert_eq!(control.str_field("material_state"), Some(binding.state.as_str()));
    close(n(control, "elements"), 1.0, 0.0);
    close(n(control, "reference_value"), 1.0, 0.0);
    close(n(control, "log_reference_value"), 0.0, 0.0);
    assert!(n(control, DK) < -1.0, "latent history must have a nonzero final effect");

    // New immutable packs and project bindings perturb the actual physical
    // conductivity, with no adjoint requested in either comparison solve.
    let epsilon = 0.001;
    let perturb = |scale: f64| {
        let (changed, changed_cards) = phase_fixture(scale);
        assert_ne!(changed.materials.as_ref().unwrap()[0].card, binding.card);
        let trial = run(&changed, &changed_cards, Some(&bytes), None);
        check_latent_exit(&trial);
        trial.maximum()
    };
    let plus = perturb(1.0 + epsilon);
    let minus = perturb(1.0 - epsilon);
    assert!(plus < minus);
    close(n(control, DK), (plus - minus) / (2.0 * epsilon), 5e-4);
}

fn contact_fixture(scale: f64) -> (ProjectSpec, CardPackSet) {
    let cards = contact_cards_with_interface(
        interface_pack_bytes_with_resistance_and_peak(0.1 * scale, None),
    );
    let mut spec = transient_enthalpy_adjoint::contact_project(0.0);
    spec.interface_cards.as_mut().unwrap()[0].card = cards.interfaces()[0].card().to_hex();
    (spec, cards)
}

fn contact_control<'a>(solved: &'a ThermalRun, spec: &ProjectSpec) -> &'a J {
    let controls = solved.controls();
    let regions = rows(controls, "conductivity_regions");
    let bindings = spec.materials.as_ref().unwrap();
    assert_eq!(regions.len(), 2, "sharing a card must not merge regional controls");
    assert_eq!(bindings[0].card, bindings[1].card);
    for (region, binding) in regions.iter().zip(bindings) {
        assert_eq!(region.str_field("region"), Some(binding.region.as_str()));
        assert_eq!(region.str_field("material_card"), Some(binding.card.as_str()));
        assert_eq!(region.str_field("material_state"), Some(binding.state.as_str()));
        assert!(n(region, "elements") > 0.0);
        close(n(region, "reference_value"), 1.0, 0.0);
        close(n(region, "log_reference_value"), 0.0, 0.0);
    }
    let contacts = rows(controls, "contact_interfaces");
    assert_eq!(contacts.len(), 1);
    let control = &contacts[0];
    let binding = &spec.interface_cards.as_ref().unwrap()[0];
    assert_eq!(control.str_field("interface"), Some(binding.interface.as_str()));
    assert_eq!(control.str_field("interface_card"), Some(binding.card.as_str()));
    assert_eq!(control.get("mapped"), Some(&J::Bool(false)));
    assert!(n(control, "face_pairs") > 0.0);
    close(n(control, "reference_value"), 1.0, 0.0);
    close(n(control, "log_reference_value"), 0.0, 0.0);
    control
}

#[test]
fn g1_native_enthalpy_contact_resistance_multiplier_matches_fresh_two_body_resolves() {
    let (spec, cards) = contact_fixture(1.0);
    let solved = run(&spec, &cards, None, Some(FINAL));
    let control = contact_control(&solved, &spec);
    close(n(solved.controls(), "control_steps"), 8.0, 0.0);
    close(n(solved.adjoint(), "value_k"), solved.maximum(), 0.0);
    assert!(solved.maximum() > 300.0 && solved.maximum() < 350.0);
    assert!(n(control, DR) > 0.01, "more joint resistance must retain heat in the hot body");
    let flux = &solved.receipt.path(&["interfaces", "fluxes"]).unwrap()
        .as_array().unwrap()[0];
    assert!(n(flux, "heat_rate_a_to_b_w") < 0.0);
    close(n(solved.transient().get("energy").unwrap(), "stored_change_j"), 0.0, 2e-6);

    let epsilon = 0.002;
    let perturb = |scale: f64| {
        let (changed, changed_cards) = contact_fixture(scale);
        assert_ne!(
            changed.interface_cards.as_ref().unwrap()[0].card,
            spec.interface_cards.as_ref().unwrap()[0].card,
        );
        run(&changed, &changed_cards, None, None).maximum()
    };
    let plus = perturb(1.0 + epsilon);
    let minus = perturb(1.0 - epsilon);
    assert!(plus > minus);
    close(n(control, DR), (plus - minus) / (2.0 * epsilon), 5e-4);
}

#[test]
fn g1_native_enthalpy_initial_peak_retains_material_identities_with_zero_future_controls() {
    let (spec, cards) = contact_fixture(1.0);
    let solved = run(&spec, &cards, None, Some(PEAK));
    let adjoint = solved.adjoint();
    close(n(adjoint, "value_k"), 350.0, 1e-12);
    assert!(solved.maximum() < 350.0);
    close(n(adjoint, "time_s"), 0.0, 0.0);
    close(n(adjoint, "state_index"), 0.0, 0.0);
    close(n(adjoint, "revalidated_solid_endpoints"), 0.0, 0.0);
    close(n(solved.controls(), "control_steps"), 0.0, 0.0);
    close(n(contact_control(&solved, &spec), DR), 0.0, 0.0);
    for region in rows(solved.controls(), "conductivity_regions") {
        close(n(region, DK), 0.0, 0.0);
    }
}
