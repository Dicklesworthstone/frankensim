//! MR20: initially hot and cold bodies exchange heat through retained contact.

use super::*;
use fs_project::{
    TransientEnthalpy, TransientEnthalpyKnot, TransientInitialTemperature,
    TransientRegionEnthalpy, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const COLD_K: f64 = 300.0;
const HOT_K: f64 = 350.0;
const CP: f64 = 10.0;
const COLD_DENSITY: f64 = 10.0;
const HOT_DENSITY: f64 = 20.0;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} differs from {expected}, tolerance {tolerance}"
    );
}

fn n(value: &J, key: &str) -> f64 {
    value.f64_field(key).unwrap()
}

fn values<'a>(row: &'a J, key: &str) -> impl Iterator<Item = f64> + 'a {
    row.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
}

fn time(spec: &mut ProjectSpec) -> &mut fs_project::ConductionTransient {
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap()
}

fn isolated_contact_project() -> ProjectSpec {
    let mut spec = multi_region_contact_project();
    transient_product::declaration(&mut spec, COLD_DENSITY * CP, 2.0);
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-11;
    spec.requirements.as_mut().unwrap()[0].region = "hot".into();
    for source in spec.power.as_mut().unwrap() {
        source.watts.value = 0.0;
    }
    for boundary in &mut spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries
    {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
    }
    let transient = time(&mut spec);
    transient.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    transient.capacities.clear();
    transient.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: [
            ("cold", COLD_DENSITY, COLD_K),
            ("hot", HOT_DENSITY, HOT_K),
        ]
        .into_iter()
        .map(|(region, density, initial_temperature)| TransientRegionEnthalpy {
            region: region.into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(
                &fixture_card_identity().0,
            )
            .unwrap(),
            source: "synthetic Cp = 10 J/(kg K) and unequal reference densities for an isolated contact oracle; not a material-property claim"
                .into(),
            reference_density: QtyAny::new(density, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(
                CP * (initial_temperature - COLD_K),
                dims::SPECIFIC_ENTHALPY,
            ),
            phase_conductivity: None,
            knots: [250.0, COLD_K, HOT_K, 390.0]
                .into_iter()
                .map(|temperature| TransientEnthalpyKnot {
                    specific_enthalpy: QtyAny::new(
                        CP * (temperature - COLD_K),
                        dims::SPECIFIC_ENTHALPY,
                    ),
                    temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                    liquid_mass_fraction: 0.0,
                })
                .collect(),
        })
        .collect(),
    });
    spec
}

#[test]
fn g1_regional_initial_enthalpy_exchanges_contact_heat_and_replays_deterministically() {
    let spec = isolated_contact_project();
    let decoded = decode(&spec);
    let cards = contact_cards();
    let ledger = Ledger::open(":memory:").unwrap();
    import_multi_region_contact(&ledger, &spec);
    let first = run_conjugate_to_completion(&ledger, &decoded, &cards);
    let receipts = stage_receipt_hashes(&ledger, &first.0);
    let replay = run_conjugate_to_completion(&ledger, &decoded, &cards);
    assert_eq!(first, replay);
    assert_eq!(receipts, stage_receipt_hashes(&ledger, &first.0));

    let independent = Ledger::open(":memory:").unwrap();
    import_multi_region_contact(&independent, &spec);
    assert_eq!(
        first,
        run_conjugate_to_completion(&independent, &decoded, &cards),
        "fresh execution and sealed replay retain identical run, QoI, and field bytes"
    );

    let conduction = J::parse(&first.2).unwrap();
    let field = J::parse(&first.3).unwrap();
    let transient = conduction.get("transient").unwrap();
    assert_eq!(
        transient.str_field("initial_state"),
        Some("regional-enthalpy")
    );
    assert_eq!(transient.get("initial_temperature_k"), Some(&J::Null));
    close(n(transient, "final_time_s"), 2.0, 0.0);
    let peak = transient.get("sampled_peak").unwrap();
    assert_eq!(peak.str_field("region"), Some("hot"));
    close(n(peak, "value_k"), HOT_K, 1e-12);
    close(n(peak, "time_s"), 0.0, 0.0);
    close(n(peak, "state"), 0.0, 0.0);
    close(n(peak, "samples"), n(transient, "fine_steps") + 1.0, 0.0);
    assert_eq!(peak.get("includes_initial"), Some(&J::Bool(true)));
    assert!(transient.get("adjoint").is_none());
    let materials = transient.get("materials").unwrap().as_array().unwrap();
    assert_eq!(materials.len(), 2);
    for (region, density, initial_temperature) in [
        ("cold", COLD_DENSITY, COLD_K),
        ("hot", HOT_DENSITY, HOT_K),
    ] {
        let material = materials
            .iter()
            .find(|row| row.str_field("region") == Some(region))
            .unwrap();
        close(n(material, "initial_temperature_k"), initial_temperature, 1e-12);
        close(n(material, "reference_density_kg_m3"), density, 0.0);
        close(
            n(material, "initial_specific_enthalpy_j_kg"),
            CP * (initial_temperature - COLD_K),
            0.0,
        );
    }

    let endpoint = transient.get("endpoint").unwrap();
    let temperatures: Vec<_> = values(&field, "temperature").collect();
    let enthalpies: Vec<_> = values(endpoint, "specific_enthalpy_j_kg").collect();
    let masses: Vec<_> = values(endpoint, "reference_nodal_mass_kg").collect();
    assert!(!temperatures.is_empty());
    assert_eq!(temperatures.len(), enthalpies.len());
    assert_eq!(temperatures.len(), masses.len());
    for (index, &temperature) in temperatures.iter().enumerate() {
        // Every cold node started at the lower endpoint and every hot node
        // at the upper endpoint. Strict interior values therefore establish
        // warming and cooling without assuming retained vertex order.
        assert!(
            temperature > COLD_K + 1e-6 && temperature < HOT_K - 1e-6,
            "each body must move inside the initial temperature bracket: {temperature}"
        );
        close(enthalpies[index], CP * (temperature - COLD_K), 1e-9);
        assert!(masses[index] > 0.0);
    }
    for fraction in values(endpoint, "liquid_mass_fraction") {
        close(fraction, 0.0, 0.0);
    }
    // Both retained regions are unit cubes. The hotter cube alone starts
    // with h = 500 J/kg. Its 20 kg fixes an absolute 10,000 J energy oracle;
    // equal weighting would be incorrect for these unequal densities.
    close(masses.iter().sum(), COLD_DENSITY + HOT_DENSITY, 3e-12);
    close(
        masses.iter().zip(&enthalpies).map(|(m, h)| m * h).sum(),
        HOT_DENSITY * CP * (HOT_K - COLD_K),
        2e-6,
    );
    let energy = transient.get("energy").unwrap();
    close(n(energy, "stored_change_j"), 0.0, 2e-6);
    close(n(energy, "integrated_net_input_j"), 0.0, 1e-12);
    close(n(energy, "window_residual_j"), 0.0, 2e-6);
    for grid in ["coarse", "fine"] {
        for row in transient.get(grid).unwrap().as_array().unwrap() {
            close(n(row, "source_w"), 0.0, 1e-12);
            close(n(row, "net_input_w"), 0.0, 1e-12);
            close(n(row, "stored_energy_change_j"), 0.0, 1e-6);
            close(n(row, "energy_residual_j"), 0.0, 1e-6);
            assert!(n(row, "residual_j") <= n(row, "threshold_j"));
        }
    }

    let fluxes = conduction
        .path(&["interfaces", "fluxes"])
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(fluxes.len(), 1);
    let flux = &fluxes[0];
    assert_eq!(flux.str_field("interface"), Some("cold-hot-joint"));
    close(n(flux, "area_m2"), 1.0, 2e-12);
    close(n(flux, "conductance_w_per_k"), 10.0, 2e-10);
    assert!(n(flux, "mean_jump_k") < -0.01);
    assert!(
        n(flux, "heat_rate_a_to_b_w") < -0.1,
        "the retained cold-to-hot orientation must report negative heat"
    );
    close(
        n(flux, "heat_rate_a_to_b_w"),
        10.0 * n(flux, "mean_jump_k"),
        2e-10,
    );
}

#[test]
fn g0_explicit_uniform_temperature_still_rejects_inconsistent_regional_enthalpy() {
    let mut spec = isolated_contact_project();
    time(&mut spec).initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(COLD_K, dims::TEMPERATURE));
    // This intentionally invalid declaration must still print and parse,
    // then stop at admission before any import or solve stage can publish.
    let source = print_sexpr(&spec).unwrap();
    let decoded = fs_project::parse_sexpr(&source).unwrap();
    assert!(
        decoded.findings().iter().any(|finding| {
            finding.code == "project-conduction-transient-enthalpy-initial-state"
        }),
        "{:?}",
        decoded.findings()
    );
    let ledger = Ledger::open(":memory:").unwrap();
    let refusal = run_solve(
        &ledger,
        &CancelGate::new_clock_free(),
        &mut benign_clock(),
        &decoded,
        &contact_cards(),
        &mut Vec::new(),
    )
    .expect_err("Uniform(300 K) cannot silently override the hot body's 350 K enthalpy");
    assert_eq!(refusal.code, "cli-solve-project-invalid");
    assert!(refusal.run.is_none());
    assert!(refusal.stage.is_none());
}

#[test]
fn g1_regional_enthalpy_peak_uses_the_selected_regions_original_chart_temperature() {
    let mut spec = isolated_contact_project();
    spec.requirements.as_mut().unwrap()[0].region = "cold".into();
    let ledger = Ledger::open(":memory:").unwrap();
    import_multi_region_contact(&ledger, &spec);
    let (_, qoi, conduction, _) =
        run_conjugate_to_completion(&ledger, &decode(&spec), &contact_cards());
    let conduction = J::parse(&conduction).unwrap();
    let transient = conduction.get("transient").unwrap();
    let peak = transient.get("sampled_peak").unwrap();
    let rows = transient.get("fine").unwrap().as_array().unwrap();
    let mut expected = (COLD_K, 0.0, 0_usize);
    for (index, row) in rows.iter().enumerate() {
        let temperature = n(row, "final_region_max_k");
        if temperature > expected.0 {
            expected = (temperature, n(row, "time_s"), index + 1);
        }
    }
    assert_eq!(peak.str_field("region"), Some("cold"));
    close(n(peak, "value_k"), expected.0, 1e-12);
    close(n(peak, "time_s"), expected.1, 0.0);
    close(n(peak, "state"), expected.2 as f64, 0.0);
    close(n(peak, "samples"), (rows.len() + 1) as f64, 0.0);
    assert!(n(peak, "value_k") > COLD_K + 1e-6);
    assert!(n(peak, "value_k") < HOT_K - 1e-6,
        "the hotter neighboring body's initial temperature is outside this QoI region");
    close(receipt_number_field(&qoi, "value"),
        n(rows.last().unwrap(), "final_region_max_k"), 0.0);
    assert_eq!(peak.get("includes_initial"), Some(&J::Bool(true)));
    assert_eq!(peak.get("continuous_time_bound"), Some(&J::Bool(false)));
    assert!(transient.get("adjoint").is_none());
}
