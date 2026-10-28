//! Native accepted-enthalpy adjoints through imported thermal meshes.
//! Scalar balances independently fix sensible and latent sensitivities; fresh
//! perturbed solves also test the regional history and reference-density chain.

use super::*;
use fs_project::{
    ConductionTransient, TransientEnthalpy, TransientEnthalpyKnot,
    TransientInitialTemperature, TransientRegionEnthalpy, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

mod prescribed;
mod flux;

const CP: f64 = 10.0;
const VOLUME: f64 = 8.0 / 3.0;
const HORIZON: f64 = 2.0;
const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn values(row: &J, key: &str) -> Vec<f64> {
    row.get(key).unwrap().as_array().unwrap().iter()
        .map(|value| value.as_f64().unwrap()).collect()
}

fn region<'a>(row: &'a J, key: &str, name: &str) -> &'a J {
    row.get(key).unwrap().as_array().unwrap().iter()
        .find(|entry| entry.str_field("region") == Some(name)).unwrap()
}

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap()
}

fn material(region: &str, density: f64, initial_h: f64, latent: bool) -> TransientRegionEnthalpy {
    let knots = if latent {
        vec![(-500.0, 250.0, 0.0), (0.0, 300.0, 0.0),
            (20.0, 300.0, 1.0), (920.0, 390.0, 1.0)]
    } else {
        vec![(-500.0, 250.0, 0.0), (900.0, 390.0, 0.0)]
    };
    TransientRegionEnthalpy {
        region: region.into(),
        material_card_identity: fs_blake3::ContentHash::from_hex(&fixture_card_identity().0)
            .unwrap(),
        source: "synthetic Cp, equilibrium latent heat and reference density for an independent derivative oracle; not measured material data".into(),
        reference_density: QtyAny::new(density, dims::DENSITY),
        initial_specific_enthalpy: QtyAny::new(initial_h, dims::SPECIFIC_ENTHALPY),
        knots: knots.into_iter().map(|(h, temperature, fraction)| TransientEnthalpyKnot {
            specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
            temperature: QtyAny::new(temperature, dims::TEMPERATURE),
            liquid_mass_fraction: fraction,
        }).collect(),
        phase_conductivity: None,
    }
}

fn install_storage(spec: &mut ProjectSpec, materials: Vec<TransientRegionEnthalpy>) {
    transient_product::declaration(spec, 1.0, HORIZON);
    let policy = time(spec);
    policy.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    policy.max_steps = 12;
    policy.capacities.clear();
    policy.enthalpy = Some(TransientEnthalpy {
        materials,
        max_newton_iterations: 32,
        max_backtracks: 24,
    });
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
}

fn uniform_project(density: f64, initial_h: f64, power: f64, latent: bool)
    -> (ProjectSpec, Vec<u8>)
{
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conduction_fixture_project(7, &bytes);
    install_storage(&mut spec, vec![material("air", density, initial_h, latent)]);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.regions[0].seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    setup.boundaries[0].condition = ThermalBoundaryCondition::HeatFlux {
        outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
    };
    spec.power.as_mut().unwrap()[0].watts.value = power;
    (spec, bytes)
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
        let result = self.transient().get("adjoint").unwrap();
        assert_eq!(result.str_field("schema"), Some("fs-cli-native-enthalpy-adjoint-v1"));
        result
    }

    fn final_temperature(&self) -> f64 {
        self.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }
}

fn run(spec: &ProjectSpec, bytes: Option<&[u8]>, cards: &CardPackSet,
    output: Option<&str>) -> ThermalRun
{
    let mut spec = spec.clone();
    if let Some(name) = output {
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: name.into(), kind: "report".into(), region: None,
        });
    }
    let ledger = Ledger::open(":memory:").unwrap();
    if let Some(bytes) = bytes {
        import_fixture(&ledger, &spec, bytes.to_vec());
    } else {
        import_multi_region_contact(&ledger, &spec);
    }
    let (_, _, receipt, field) =
        run_conjugate_to_completion(&ledger, &decode(&spec), cards);
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        temperature: values(&J::parse(&field).unwrap(), "temperature"),
    }
}

fn initial_gradient(row: &J) -> f64 {
    n(row, "dtemperature_dinitial_specific_enthalpy_k_kg_per_j")
}

fn density_gradient(row: &J) -> f64 {
    n(row, "dtemperature_dreference_density_k_m3_per_kg")
}

fn check_final_adjoint(solved: &ThermalRun) {
    let gradient = solved.adjoint();
    assert_eq!(gradient.str_field("output"), Some(FINAL));
    close(n(gradient, "value_k"), solved.final_temperature(), 0.0);
    close(n(gradient, "time_s"), HORIZON, 0.0);
    let fine = solved.transient().get("fine").unwrap().as_array().unwrap();
    close(n(gradient, "state_index"), fine.len() as f64, 0.0);
    close(n(gradient, "revalidated_solid_endpoints"), fine.len() as f64, 0.0);
    assert!(n(gradient, "maximum_dual_relative_residual") < 1e-10);
    assert!(n(gradient, "adjoint_krylov_iterations")
        <= n(gradient, "max_adjoint_krylov_iterations"));
    let initial = values(gradient, "dtemperature_dinitial_specific_enthalpy_k_kg_per_j");
    assert_eq!(initial.len(), solved.temperature.len());
    close(initial.iter().sum(),
        n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"), 1e-12);
}

#[test]
fn g1_native_enthalpy_adjoint_sensible_controls_match_energy_and_fresh_resolves() {
    let density = 10.0;
    let power = 40.0;
    let (spec, bytes) = uniform_project(density, 0.0, power, false);
    let cards = fixture_cards();
    let solved = run(&spec, Some(&bytes), &cards, Some(FINAL));
    check_final_adjoint(&solved);
    assert_eq!(solved.temperature.len(), 4);
    let mass = density * VOLUME;
    let rise = power * HORIZON / (mass * CP);
    for &temperature in &solved.temperature {
        close(temperature, 300.0 + rise, 2e-8);
    }
    let gradient = solved.adjoint();
    let storage = region(gradient, "regional_storage", "air");
    let source = region(gradient, "regional_power", "air");
    let dh = 1.0 / CP;
    let dp = HORIZON / (mass * CP);
    let drho = -rise / density;
    close(n(storage, "reference_density_kg_m3"), density, 0.0);
    close(initial_gradient(storage), dh, 1e-9);
    close(n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"), dh, 1e-9);
    close(n(source, "dtemperature_duniform_added_watts_k_per_w"), dp, 1e-10);
    close(density_gradient(storage), drho, 1e-9);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"), -rise, 1e-8);

    for interval in gradient.get("intervals").unwrap().as_array().unwrap() {
        let dt = n(interval, "end_s") - n(interval, "start_s");
        let expected = dt / (mass * CP);
        close(n(region(interval, "regional_power", "air"),
            "dtemperature_ddelivered_watts_k_per_w"), expected, 1e-10);
        close(n(interval, "dtemperature_dsource_multiplier_k"), power * expected, 1e-9);
    }

    // Perturb all three independent declarations, and re-solve without an
    // adjoint request. This probes the physical result, not the derivative code.
    let epsilon = 0.02;
    let perturb = |delta: f64| {
        let (changed, _) = uniform_project(
            density + 0.4 * delta, 2.0 * delta, power + 3.0 * delta, false,
        );
        run(&changed, Some(&bytes), &cards, None).final_temperature()
    };
    let fd = (perturb(epsilon) - perturb(-epsilon)) / (2.0 * epsilon);
    let directional = 2.0 * initial_gradient(storage)
        + 3.0 * n(source, "dtemperature_duniform_added_watts_k_per_w")
        + 0.4 * density_gradient(storage);
    close(directional, 2.0 * dh + 3.0 * dp + 0.4 * drho, 1e-8);
    close(directional, fd, 2e-6);
}

#[test]
fn g1_native_enthalpy_adjoint_carries_latent_history_into_later_sensible_temperature() {
    let (spec, bytes) = uniform_project(1.0, 10.0, 40.0, true);
    let cards = fixture_cards();
    let solved = run(&spec, Some(&bytes), &cards, Some(FINAL));
    check_final_adjoint(&solved);
    let gradient = solved.adjoint();
    let storage = region(gradient, "regional_storage", "air");
    let final_h = 10.0 + 40.0 * HORIZON / VOLUME;
    close(final_h, 40.0, 1e-12);
    close(solved.final_temperature(), 300.0 + (final_h - 20.0) / CP, 2e-8);
    let fine = solved.transient().get("fine").unwrap().as_array().unwrap();
    assert!(n(&fine[0], "mean_liquid_mass_fraction") > 0.5);
    assert!(n(&fine[0], "mean_liquid_mass_fraction") < 1.0);
    close(n(fine.last().unwrap(), "mean_liquid_mass_fraction"), 1.0, 0.0);
    close(initial_gradient(storage), 1.0 / CP, 1e-9);
    close(n(region(gradient, "regional_power", "air"),
        "dtemperature_duniform_added_watts_k_per_w"), HORIZON / (VOLUME * CP), 1e-9);
    close(density_gradient(storage), -40.0 * HORIZON / (VOLUME * CP), 1e-8);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"),
        -40.0 * HORIZON / (VOLUME * CP), 1e-8);

    let epsilon = 0.25;
    let (plus_spec, _) = uniform_project(1.0, 10.0 + epsilon, 40.0, true);
    let (minus_spec, _) = uniform_project(1.0, 10.0 - epsilon, 40.0, true);
    let plus = run(&plus_spec, Some(&bytes), &cards, None);
    let minus = run(&minus_spec, Some(&bytes), &cards, None);
    for trial in [&solved, &plus, &minus] {
        // The initial temperature is identical: dT/dh is zero inside this
        // plateau. A temperature-only history or an inverted initial slope
        // would erase the nonzero final sensitivity measured below.
        close(n(region(trial.transient(), "materials", "air"),
            "initial_temperature_k"), 300.0, 0.0);
    }
    let fd = (plus.final_temperature() - minus.final_temperature()) / (2.0 * epsilon);
    close(fd, 1.0 / CP, 2e-7);
    close(initial_gradient(storage), fd, 2e-7);
}

#[test]
fn g1_native_enthalpy_adjoint_initial_sampled_peak_has_no_future_control_effect() {
    let (mut spec, bytes) = uniform_project(10.0, 300.0, 0.0, false);
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .boundaries[0].condition = ThermalBoundaryCondition::Convection {
            coefficient: QtyAny::new(8.0, dims::HEAT_TRANSFER_COEFFICIENT),
            reference_temperature: QtyAny::new(293.15, dims::TEMPERATURE),
        };
    let solved = run(&spec, Some(&bytes), &fixture_cards(), Some(PEAK));
    let gradient = solved.adjoint();
    assert_eq!(gradient.str_field("output"), Some(PEAK));
    close(n(gradient, "value_k"), 330.0, 1e-12);
    assert!(solved.final_temperature() < 329.0);
    close(n(gradient, "time_s"), 0.0, 0.0);
    close(n(gradient, "state_index"), 0.0, 0.0);
    close(n(gradient, "active_ties"), 4.0, 0.0);
    assert_eq!(gradient.get("unique_sampled_maximum"), Some(&J::Bool(false)));
    close(n(gradient, "revalidated_solid_endpoints"), 0.0, 0.0);
    close(n(gradient, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"), 0.0, 0.0);
    let initial = values(gradient, "dtemperature_dinitial_specific_enthalpy_k_kg_per_j");
    assert_eq!(initial.len(), 4);
    let active = n(gradient, "active_vertex") as usize;
    assert!(active < initial.len());
    for (vertex, &value) in initial.iter().enumerate() {
        close(value, if vertex == active { 1.0 / CP } else { 0.0 }, 1e-14);
    }
    close(n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"),
        1.0 / CP, 1e-14);
    let storage = region(gradient, "regional_storage", "air");
    close(initial_gradient(storage), 1.0 / CP, 1e-14);
    close(density_gradient(storage), 0.0, 0.0);
    close(n(region(gradient, "regional_power", "air"),
        "dtemperature_duniform_added_watts_k_per_w"), 0.0, 0.0);
    for interval in gradient.get("intervals").unwrap().as_array().unwrap() {
        close(n(interval, "dtemperature_dsource_multiplier_k"), 0.0, 0.0);
        for source in interval.get("regional_power").unwrap().as_array().unwrap() {
            close(n(source, "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
        }
    }
}

pub(super) fn contact_project(delta: f64) -> ProjectSpec {
    let mut spec = multi_region_contact_project();
    install_storage(&mut spec, vec![
        material("cold", 10.0 + 0.3 * delta, delta, false),
        material("hot", 20.0 - 0.4 * delta, 500.0 - 2.0 * delta, false),
    ]);
    spec.requirements.as_mut().unwrap()[0].region = "hot".into();
    for source in spec.power.as_mut().unwrap() {
        source.watts.value = 0.0;
    }
    for boundary in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
    }
    spec
}

#[test]
fn g1_native_enthalpy_adjoint_resolves_unequal_regional_history_and_density_through_contact() {
    let cards = contact_cards();
    let solved = run(&contact_project(0.0), None, &cards, Some(FINAL));
    check_final_adjoint(&solved);
    assert!(solved.final_temperature() > 300.0 && solved.final_temperature() < 350.0);
    let gradient = solved.adjoint();
    let cold = region(gradient, "regional_storage", "cold");
    let hot = region(gradient, "regional_storage", "hot");
    close(n(cold, "reference_density_kg_m3"), 10.0, 0.0);
    close(n(hot, "reference_density_kg_m3"), 20.0, 0.0);
    assert!(initial_gradient(cold) > 0.0);
    assert!(initial_gradient(hot) > 0.0);
    assert!(density_gradient(cold) < 0.0, "a heavier cold body must lower the hot maximum");
    assert!(density_gradient(hot) > 0.0, "a heavier hot body must retain more heat");
    // A common h shift raises both bodies by exactly dh/Cp. Their contact
    // exchange is unchanged, despite their different reference masses.
    close(initial_gradient(cold) + initial_gradient(hot), 1.0 / CP, 1e-8);
    close(n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"),
        1.0 / CP, 1e-8);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"),
        10.0 * density_gradient(cold) + 20.0 * density_gradient(hot), 1e-10);
    close(n(solved.transient().get("energy").unwrap(), "stored_change_j"), 0.0, 2e-6);
    let flux = &solved.receipt.path(&["interfaces", "fluxes"]).unwrap()
        .as_array().unwrap()[0];
    assert!(n(flux, "heat_rate_a_to_b_w") < 0.0);

    let epsilon = 0.01;
    let plus = run(&contact_project(epsilon), None, &cards, None).final_temperature();
    let minus = run(&contact_project(-epsilon), None, &cards, None).final_temperature();
    let fd = (plus - minus) / (2.0 * epsilon);
    let directional = initial_gradient(cold) - 2.0 * initial_gradient(hot)
        + 0.3 * density_gradient(cold) - 0.4 * density_gradient(hot);
    assert!(directional.abs() > 1e-3);
    close(directional, fd, 5e-5);
}
