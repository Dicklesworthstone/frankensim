//! Native radiative enthalpy adjoints against independent scalar heat balances.
//!
//! A regular tetrahedron has one exact uniform P1 mode. The oracle differentiates
//! its physical backward-Euler equation, including 4 epsilon sigma T^3, and
//! never calls the conduction, enthalpy-chart or radiation derivative owners.

use super::*;
use fs_project::{
    ConductionTransient, TransientEnthalpy, TransientEnthalpyKnot,
    TransientInitialTemperature, TransientRegionEnthalpy, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const SIGMA: f64 = 5.670_374_419e-8;
const EMISSIVITY: f64 = 0.85;
const VOLUME: f64 = 8.0 / 3.0;
const CP: f64 = 10.0;
const INITIAL_K: f64 = 330.0;
const AMBIENT_K: f64 = 293.15;
const HTC: f64 = 10.0;
const LATENT_H: f64 = 200.0;
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

fn region<'a>(row: &'a J, key: &str) -> &'a J {
    row.get(key).unwrap().as_array().unwrap().iter()
        .find(|entry| entry.str_field("region") == Some("air")).unwrap()
}

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap()
}

#[derive(Clone, Copy)]
struct Inputs {
    density: f64,
    initial_h: f64,
    power: f64,
    reservoir: f64,
    latent: bool,
}

impl Inputs {
    fn sensible(reservoir: f64) -> Self {
        Self { density: 10.0, initial_h: 0.0, power: 5.0, reservoir, latent: false }
    }

    fn perturb(self, delta: f64) -> Self {
        Self {
            density: self.density + 0.4 * delta,
            initial_h: self.initial_h + 2.0 * delta,
            power: self.power + 3.0 * delta,
            ..self
        }
    }
}

fn project(card: &str, input: Inputs) -> ProjectSpec {
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conduction_fixture_project(7, &bytes);
    transient_product::declaration(&mut spec, input.density * CP, 1.0);
    radiation_product::declare(&mut spec, card, "air", input.reservoir);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.regions[0].seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    setup.boundaries[0].condition = ThermalBoundaryCondition::Convection {
        coefficient: QtyAny::new(HTC, dims::HEAT_TRANSFER_COEFFICIENT),
        reference_temperature: QtyAny::new(AMBIENT_K, dims::TEMPERATURE),
    };
    let knots = if input.latent {
        vec![
            (CP * (233.15 - INITIAL_K), 233.15, 0.0),
            (0.0, INITIAL_K, 0.0),
            (LATENT_H, INITIAL_K, 1.0),
            (LATENT_H + CP * (398.15 - INITIAL_K), 398.15, 1.0),
        ]
    } else {
        vec![
            (CP * (233.15 - INITIAL_K), 233.15, 0.0),
            (CP * (398.15 - INITIAL_K), 398.15, 0.0),
        ]
    };
    let policy = time(&mut spec);
    policy.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    policy.max_steps = 6;
    policy.capacities.clear();
    policy.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            phase_conductivity: None,
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&fixture_card_identity().0)
                .unwrap(),
            source: "synthetic Cp, latent interval and reference density for an independent radiation-gradient oracle; not measured material data".into(),
            reference_density: QtyAny::new(input.density, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(input.initial_h, dims::SPECIFIC_ENTHALPY),
            knots: knots.into_iter().map(|(h, temperature, fraction)| TransientEnthalpyKnot {
                specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                liquid_mass_fraction: fraction,
            }).collect(),
        }],
    });
    spec.power.as_mut().unwrap()[0].watts.value = input.power;
    spec
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
        let gradient = self.transient().get("adjoint").unwrap();
        assert_eq!(gradient.str_field("schema"), Some("fs-cli-native-enthalpy-adjoint-v1"));
        gradient
    }

    fn final_temperature(&self) -> f64 {
        self.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }
}

fn run(input: Inputs, cards: &CardPackSet, card: &str, output: Option<&str>) -> ThermalRun {
    let mut spec = project(card, input);
    if let Some(name) = output {
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: name.into(), kind: "report".into(), region: None,
        });
    }
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, transient_radiation::regular_tetra_stl());
    let decoded = decode(&spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let (_, _, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        temperature: values(&J::parse(&field).unwrap(), "temperature"),
    }
}

fn state(h: f64, latent: bool) -> (f64, f64) {
    if !latent || h < 0.0 {
        (INITIAL_K + h / CP, 1.0 / CP)
    } else if h < LATENT_H {
        (INITIAL_K, 0.0)
    } else {
        (INITIAL_K + (h - LATENT_H) / CP, 1.0 / CP)
    }
}

struct Oracle {
    temperatures: Vec<f64>,
    enthalpies: Vec<f64>,
    initial: f64,
    density: f64,
    emissivity: f64,
    reservoir: f64,
    power_intervals: Vec<f64>,
    frozen_initial: f64,
}

fn oracle(input: Inputs) -> Oracle {
    let area = 8.0 * 3.0_f64.sqrt();
    let mass = input.density * VOLUME;
    let dt = 0.25;
    let mut old = input.initial_h;
    let mut initial = 1.0;
    let mut density = 0.0;
    let mut emissivity = 0.0;
    let mut reservoir = 0.0;
    let mut power_intervals = vec![0.0; 4];
    let mut frozen_initial = 1.0;
    let mut temperatures = Vec::new();
    let mut enthalpies = Vec::new();
    let mut final_slope = 0.0;
    for index in 0..4 {
        let balance = |h: f64| {
            let (temperature, _) = state(h, input.latent);
            let convection = HTC * area * (temperature - AMBIENT_K);
            let radiation = EMISSIVITY * SIGMA * area
                * (temperature.powi(4) - input.reservoir.powi(4));
            mass * (h - old) + dt * (convection + radiation - input.power)
        };
        let (mut low, mut high) = (
            CP * (233.15 - INITIAL_K),
            CP * (398.15 - INITIAL_K) + if input.latent { LATENT_H } else { 0.0 },
        );
        assert!(balance(low) < 0.0 && balance(high) > 0.0);
        for _ in 0..100 {
            let middle = 0.5 * (low + high);
            if balance(middle) > 0.0 { high = middle; } else { low = middle; }
        }
        let h = 0.5 * (low + high);
        let (temperature, slope) = state(h, input.latent);
        // d[M(h-h_old)+dt(Q_conv(T(h))+Q_rad(T(h))-P)]/dh.
        // This independent scalar law includes the entire radiative feedback.
        let diagonal = mass + dt * area
            * (HTC + 4.0 * EMISSIVITY * SIGMA * temperature.powi(3)) * slope;
        initial *= mass / diagonal;
        density = (mass * density - VOLUME * (h - old)) / diagonal;
        // The declarations are constant over the trajectory, so their direct
        // residual partials accumulate at every consumed time step. A plateau
        // has zero T'(h), but retains both boundary inputs in its h history.
        emissivity = (mass * emissivity - dt * SIGMA * area
            * (temperature.powi(4) - input.reservoir.powi(4))) / diagonal;
        reservoir = (mass * reservoir + dt * 4.0 * EMISSIVITY * SIGMA * area
            * input.reservoir.powi(3)) / diagonal;
        for derivative in &mut power_intervals[..index] {
            *derivative *= mass / diagonal;
        }
        power_intervals[index] = dt / diagonal;
        // Negative control: differentiate with the converged radiative secant
        // incorrectly frozen. It must be distinguishable from the actual law.
        let secant = EMISSIVITY * SIGMA * (temperature + input.reservoir)
            * (temperature.powi(2) + input.reservoir.powi(2));
        frozen_initial *= mass / (mass + dt * area * (HTC + secant) * slope);
        temperatures.push(temperature);
        enthalpies.push(h);
        old = h;
        final_slope = slope;
    }
    for derivative in &mut power_intervals { *derivative *= final_slope; }
    Oracle {
        temperatures, enthalpies,
        initial: final_slope * initial,
        density: final_slope * density,
        emissivity: final_slope * emissivity,
        reservoir: final_slope * reservoir,
        power_intervals,
        frozen_initial: final_slope * frozen_initial,
    }
}

fn initial_gradient(gradient: &J) -> f64 {
    n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j")
}

fn density_gradient(gradient: &J) -> f64 {
    n(region(gradient, "regional_storage"), "dtemperature_dreference_density_k_m3_per_kg")
}

fn power_gradient(gradient: &J) -> f64 {
    n(region(gradient, "regional_power"), "dtemperature_duniform_added_watts_k_per_w")
}

fn radiation_surface(row: &J) -> &J {
    let surfaces = row.get("radiation").unwrap().get("surfaces").unwrap()
        .as_array().unwrap();
    assert_eq!(surfaces.len(), 1);
    &surfaces[0]
}

fn check_radiation_controls(solved: &ThermalRun, expected: &Oracle, input: Inputs) {
    let surface = radiation_surface(solved.adjoint());
    assert_eq!(surface.str_field("name"), Some("gray-exterior"));
    assert_eq!(surface.str_field("target"), Some("air"));
    close(n(surface, "emissivity"), EMISSIVITY, 0.0);
    close(n(surface, "reservoir_temperature_k"), input.reservoir, 0.0);
    let emissivity = n(surface, "dtemperature_demissivity_k");
    let reservoir = n(surface, "dtemperature_dreservoir_temperature_k_per_k");
    close(emissivity, expected.emissivity, 2e-6);
    close(reservoir, expected.reservoir, 2e-7);
    assert_eq!(emissivity > 0.0, input.reservoir > INITIAL_K);
    assert!(reservoir > 0.0);

    // Edit actual native declarations and query newly admitted immutable
    // emissivity packs. Neither an accepted patch nor a derived Robin row is
    // modified. Geometry, reference mass, initial h, sources and grid stay fixed.
    let source_card = radiation_surface(&solved.receipt).str_field("card").unwrap();
    let evaluate = |delta: f64| {
        let changed_emissivity = EMISSIVITY + 0.02 * delta;
        let (cards, card) = radiation_product::emissivity_cards(
            &fixture_cards(), changed_emissivity,
        );
        assert_ne!(card.as_str(), source_card);
        let changed = Inputs { reservoir: input.reservoir + 0.5 * delta, ..input };
        let trial = run(changed, &cards, &card, None);
        let retained = radiation_surface(&trial.receipt);
        assert_eq!(retained.str_field("card"), Some(card.as_str()));
        close(n(retained, "emissivity"), changed_emissivity, 0.0);
        close(n(retained, "ambient_temperature_k"), changed.reservoir, 0.0);
        trial.final_temperature()
    };
    let step = 0.02;
    let difference = (evaluate(step) - evaluate(-step)) / (2.0 * step);
    close(0.02 * emissivity + 0.5 * reservoir, difference, 3e-6);
}

fn check_final(solved: &ThermalRun, expected: &Oracle, input: Inputs) {
    assert_eq!(solved.temperature.len(), 4, "oracle requires one regular tetra");
    let gradient = solved.adjoint();
    assert_eq!(gradient.str_field("output"), Some(FINAL));
    close(n(gradient, "state_index"), 4.0, 0.0);
    close(n(gradient, "time_s"), 1.0, 0.0);
    close(n(gradient, "value_k"), solved.final_temperature(), 0.0);
    close(n(gradient, "revalidated_solid_endpoints"), 4.0, 0.0);
    assert!(n(gradient, "maximum_dual_relative_residual") < 1e-10);
    assert!(n(gradient, "adjoint_krylov_iterations")
        <= n(gradient, "max_adjoint_krylov_iterations"));
    close(initial_gradient(gradient), expected.initial, 2e-7);
    close(density_gradient(gradient), expected.density, 2e-7);
    close(power_gradient(gradient), expected.power_intervals.iter().sum(), 2e-8);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"),
        input.density * expected.density, 2e-6);
    close(n(region(gradient, "regional_storage"),
        "dtemperature_dinitial_specific_enthalpy_k_kg_per_j"), expected.initial, 2e-7);
    close(values(gradient, "dtemperature_dinitial_specific_enthalpy_k_kg_per_j")
        .iter().sum(), expected.initial, 2e-7);
    let fine = solved.transient().get("fine").unwrap().as_array().unwrap();
    assert_eq!(fine.len(), 4);
    for (row, &temperature) in fine.iter().zip(&expected.temperatures) {
        close(n(row, "final_region_max_k"), temperature, 3e-8);
        let radiation = row.get("radiation").unwrap();
        assert!(n(radiation, "physical_residual_norm_j")
            <= n(radiation, "physical_residual_tolerance_j"));
    }
    for &temperature in &solved.temperature {
        close(temperature, *expected.temperatures.last().unwrap(), 3e-8);
    }
    let intervals = gradient.get("intervals").unwrap().as_array().unwrap();
    assert_eq!(intervals.len(), 4);
    for (interval, &derivative) in intervals.iter().zip(&expected.power_intervals) {
        close(n(interval, "end_s") - n(interval, "start_s"), 0.25, 0.0);
        close(n(region(interval, "regional_power"),
            "dtemperature_ddelivered_watts_k_per_w"), derivative, 2e-8);
        close(n(interval, "dtemperature_dsource_multiplier_k"),
            input.power * derivative, 2e-8);
    }
    check_radiation_controls(solved, expected, input);
}

#[test]
fn g1_native_radiative_enthalpy_gradients_match_physical_roots_and_fresh_resolves() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for reservoir in [270.0, 400.0] {
        let input = Inputs::sensible(reservoir);
        let expected = oracle(input);
        let solved = run(input, &cards, &card, Some(FINAL));
        check_final(&solved, &expected, input);
        let gradient = solved.adjoint();
        assert!((expected.initial - expected.frozen_initial).abs() > 0.002,
            "the fixture must distinguish the full radiation tangent from a frozen secant");
        assert_eq!(solved.final_temperature() > INITIAL_K, reservoir > INITIAL_K);
        assert_eq!(density_gradient(gradient) > 0.0, reservoir < INITIAL_K);
        let step = 0.02;
        let plus = run(input.perturb(step), &cards, &card, None).final_temperature();
        let minus = run(input.perturb(-step), &cards, &card, None).final_temperature();
        let difference = (plus - minus) / (2.0 * step);
        let derivative = 2.0 * initial_gradient(gradient)
            + 3.0 * power_gradient(gradient) + 0.4 * density_gradient(gradient);
        close(derivative, difference, 3e-6);
    }
}

#[test]
fn g1_native_radiative_enthalpy_reverse_preserves_plateau_history_before_sensible_heating() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let input = Inputs { latent: true, initial_h: 80.0, ..Inputs::sensible(400.0) };
    let expected = oracle(input);
    let solved = run(input, &cards, &card, Some(FINAL));
    check_final(&solved, &expected, input);
    let fine = solved.transient().get("fine").unwrap().as_array().unwrap();
    for (row, &h) in fine[..3].iter().zip(&expected.enthalpies) {
        assert!(h > 0.0 && h < LATENT_H);
        close(n(row, "final_region_max_k"), INITIAL_K, 0.0);
        close(n(row, "mean_liquid_mass_fraction"), h / LATENT_H, 2e-9);
    }
    assert!(*expected.enthalpies.last().unwrap() > LATENT_H);
    close(n(fine.last().unwrap(), "mean_liquid_mass_fraction"), 1.0, 0.0);
    assert!(initial_gradient(solved.adjoint()) > 0.08);
    let step = 0.25;
    let plus = run(Inputs { initial_h: input.initial_h + step, ..input },
        &cards, &card, None);
    let minus = run(Inputs { initial_h: input.initial_h - step, ..input },
        &cards, &card, None);
    for trial in [&solved, &plus, &minus] {
        close(n(region(trial.transient(), "materials"), "initial_temperature_k"),
            INITIAL_K, 0.0);
    }
    close(initial_gradient(solved.adjoint()),
        (plus.final_temperature() - minus.final_temperature()) / (2.0 * step), 2e-7);
}

#[test]
fn g1_native_radiative_initial_peak_has_no_future_source_or_storage_mass_effect() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let solved = run(Inputs::sensible(270.0), &cards, &card, Some(PEAK));
    let gradient = solved.adjoint();
    assert_eq!(gradient.str_field("output"), Some(PEAK));
    close(n(gradient, "value_k"), INITIAL_K, 1e-12);
    assert!(solved.final_temperature() < INITIAL_K - 20.0);
    close(n(gradient, "state_index"), 0.0, 0.0);
    close(n(gradient, "time_s"), 0.0, 0.0);
    close(n(gradient, "revalidated_solid_endpoints"), 0.0, 0.0);
    close(n(gradient, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(initial_gradient(gradient), 1.0 / CP, 1e-14);
    close(density_gradient(gradient), 0.0, 0.0);
    close(power_gradient(gradient), 0.0, 0.0);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"), 0.0, 0.0);
    let surface = radiation_surface(gradient);
    assert_eq!(surface.str_field("name"), Some("gray-exterior"));
    assert_eq!(surface.str_field("target"), Some("air"));
    close(n(surface, "emissivity"), EMISSIVITY, 0.0);
    close(n(surface, "reservoir_temperature_k"), 270.0, 0.0);
    close(n(surface, "dtemperature_demissivity_k"), 0.0, 0.0);
    close(n(surface, "dtemperature_dreservoir_temperature_k_per_k"), 0.0, 0.0);
    for interval in gradient.get("intervals").unwrap().as_array().unwrap() {
        close(n(interval, "dtemperature_dsource_multiplier_k"), 0.0, 0.0);
        close(n(region(interval, "regional_power"),
            "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
    }
}
