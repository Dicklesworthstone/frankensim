//! Native convective enthalpy adjoints against independent physical derivatives.
//! The scalar natural-convection oracle differentiates the film temperature
//! and both density factors in Ra. Fresh native solves exercise the streamed
//! air reference and spatial mixed-boundary chains without freezing feedback.

use super::*;
use fs_project::{
    ConductionTransient, TransientEnthalpy, TransientEnthalpyKnot, TransientInitialTemperature,
    TransientRegionEnthalpy, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const CP: f64 = 10.0;
const BASE_K: f64 = 300.0;
const INLET_K: f64 = 293.15;
const VOLUME: f64 = 8.0 / 3.0;
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
    row.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
        .collect()
}

fn region<'a>(row: &'a J, key: &str) -> &'a J {
    row.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| entry.str_field("region") == Some("air"))
        .unwrap()
}

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
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

#[derive(Clone, Copy)]
struct Inputs {
    initial_h: f64,
    density: f64,
    power: f64,
}

impl Inputs {
    fn hot() -> Self {
        Self {
            initial_h: CP * (330.0 - BASE_K),
            density: 10.0,
            power: 5.0,
        }
    }
}

fn storage(spec: &mut ProjectSpec, input: Inputs, horizon: f64) {
    transient_product::declaration(spec, input.density * CP, horizon);
    let policy = time(spec);
    policy.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    policy.max_step.value = horizon / 2.0;
    policy.max_steps = 6;
    policy.capacities.clear();
    policy.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 16,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            region: "air".into(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&fixture_card_identity().0)
                .unwrap(),
            source: "synthetic fixed Cp chart and reference density for native convective-gradient checks; not measured material data".into(),
            reference_density: QtyAny::new(input.density, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(input.initial_h, dims::SPECIFIC_ENTHALPY),
            knots: [233.15, 398.15].into_iter().map(|temperature| TransientEnthalpyKnot {
                specific_enthalpy: QtyAny::new(CP * (temperature - BASE_K), dims::SPECIFIC_ENTHALPY),
                temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                liquid_mass_fraction: 0.0,
            }).collect(),
            phase_conductivity: None,
        }],
    });
    spec.power.as_mut().unwrap()[0].watts.value = input.power;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
}

fn natural_condition(ambient: f64) -> ThermalBoundaryCondition {
    ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    }
}

fn uniform_project(input: Inputs, natural: bool) -> (ProjectSpec, Vec<u8>) {
    let bytes = transient_radiation::regular_tetra_stl();
    let mut spec = conjugate_fixture_project(7, &bytes, "convection.gnielinski");
    storage(&mut spec, input, 1.0);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.regions[0].seed = [0.0; 3].map(|value| QtyAny::new(value, dims::LENGTH));
    if natural {
        setup.boundaries[0].condition = natural_condition(INLET_K);
    }
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
        let gradient = self.transient().get("adjoint").unwrap();
        assert_eq!(
            gradient.str_field("schema"),
            Some("fs-cli-native-enthalpy-adjoint-v1")
        );
        gradient
    }

    fn final_temperature(&self) -> f64 {
        self.temperature
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

fn run(spec: &ProjectSpec, bytes: &[u8], cards: &CardPackSet, output: Option<&str>) -> ThermalRun {
    let mut spec = spec.clone();
    if let Some(name) = output {
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: name.into(),
            kind: "report".into(),
            region: None,
        });
    }
    let decoded = decode(&spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, bytes.to_vec());
    let (_, _, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        temperature: values(&J::parse(&field).unwrap(), "temperature"),
    }
}

fn gradients(gradient: &J) -> [f64; 3] {
    [
        n(
            gradient,
            "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j",
        ),
        n(
            region(gradient, "regional_power"),
            "dtemperature_duniform_added_watts_k_per_w",
        ),
        n(
            region(gradient, "regional_storage"),
            "dtemperature_dreference_density_k_m3_per_kg",
        ),
    ]
}

fn convection_control<'a>(gradient: &'a J, list: &str, key: &str, name: &str) -> &'a J {
    gradient
        .get("convection")
        .unwrap()
        .get(list)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field(key) == Some(name))
        .unwrap()
}

fn set_natural_ambient(spec: &mut ProjectSpec, target: &str, ambient: f64) {
    let boundary = spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries
        .iter_mut()
        .find(|row| row.target == target)
        .unwrap();
    let ThermalBoundaryCondition::NaturalConvection {
        ambient_temperature,
        ..
    } = &mut boundary.condition
    else {
        panic!("expected natural-convection target")
    };
    ambient_temperature.value = ambient;
}

fn set_air_inlet(spec: &mut ProjectSpec, inlet: f64) {
    let mut count = 0;
    for boundary in &mut spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries
    {
        if let ThermalBoundaryCondition::AirflowConvection {
            branch,
            inlet_temperature,
            ..
        } = &mut boundary.condition
        {
            assert_eq!(branch, "air");
            inlet_temperature.value = inlet;
            count += 1;
        }
    }
    assert!(count > 0);
}

fn set_fan_speed(spec: &mut ProjectSpec, speed: f64) {
    let system = spec.cooling.as_mut().unwrap().fan_system.as_mut().unwrap();
    assert_eq!(system.banks.len(), 1);
    assert_eq!(system.banks[0].bank_id, "fixture-bank");
    system.banks[0].speed_ratio = speed;
}

fn check_final(solved: &ThermalRun, horizon: f64) {
    let gradient = solved.adjoint();
    assert_eq!(gradient.str_field("output"), Some(FINAL));
    close(n(gradient, "value_k"), solved.final_temperature(), 0.0);
    close(n(gradient, "time_s"), horizon, 0.0);
    close(n(gradient, "state_index"), 4.0, 0.0);
    close(n(gradient, "revalidated_solid_endpoints"), 4.0, 0.0);
    assert!(n(gradient, "maximum_dual_relative_residual") < 1e-9);
    assert!(
        n(gradient, "adjoint_krylov_iterations") <= n(gradient, "max_adjoint_krylov_iterations")
    );
    close(
        values(
            gradient,
            "dtemperature_dinitial_specific_enthalpy_k_kg_per_j",
        )
        .iter()
        .sum(),
        gradients(gradient)[0],
        1e-12,
    );
    close(
        n(
            region(gradient, "regional_storage"),
            "dtemperature_dinitial_specific_enthalpy_k_kg_per_j",
        ),
        gradients(gradient)[0],
        1e-12,
    );
}

/// h, dQ/dT, the incorrect frozen-density derivative, and dQ/dTa.
/// Ra = C |T-Ta| / Tf^3: one film factor is beta and two are from nu^2.
fn natural_law(wall: f64) -> (f64, f64, f64, f64) {
    let film = 0.5 * (wall + INLET_K);
    let delta = wall - INLET_K;
    assert!(delta.abs() > 0.1);
    let density = 101325.0 / (287.05 * film);
    let nu = 1.846e-5 / density;
    let ra = 9.80665 * delta.abs() * 0.1_f64.powi(3) / (film * nu * (nu / 0.707));
    let b = 0.387 / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0);
    let r = ra.powf(1.0 / 6.0);
    let f = 0.825 + b * r;
    let h = 0.0263 / 0.1 * f * f;
    let dh_dlog_ra = 0.0263 / 0.1 * f * b * r / 3.0;
    let area = 8.0 * 3.0_f64.sqrt();
    let full = area * (h + delta * dh_dlog_ra * (1.0 / delta - 1.5 / film));
    let fixed_density = area * (h + delta * dh_dlog_ra * (1.0 / delta - 0.5 / film));
    let ambient = area * (-h + delta * dh_dlog_ra * (-1.0 / delta - 1.5 / film));
    (h, full, fixed_density, ambient)
}

struct NaturalOracle {
    temperatures: Vec<f64>,
    gradients: [f64; 3],
    intervals: Vec<f64>,
    frozen_h_initial: f64,
    frozen_density_initial: f64,
    ambient: f64,
}

fn natural_oracle(input: Inputs) -> NaturalOracle {
    let area = 8.0 * 3.0_f64.sqrt();
    let mass = input.density * VOLUME;
    let dt = 0.25;
    let mut old = input.initial_h;
    let mut derivative = [1.0, 0.0, 0.0];
    let mut intervals = vec![0.0; 4];
    let mut frozen_h_initial = 1.0;
    let mut frozen_density_initial = 1.0;
    let mut ambient = 0.0;
    let mut temperatures = Vec::new();
    for index in 0..4 {
        let residual = |h: f64| {
            let t = BASE_K + h / CP;
            mass * (h - old) + dt * (area * natural_law(t).0 * (t - INLET_K) - input.power)
        };
        // Stay on the same heated/cooled branch while bracketing the root;
        // the natural card deliberately does not admit exact equilibrium.
        let initial_t = BASE_K + input.initial_h / CP;
        let (mut low, mut high) = if initial_t > INLET_K {
            (CP * (INLET_K + 0.2 - BASE_K), CP * (390.0 - BASE_K))
        } else {
            (CP * (240.0 - BASE_K), CP * (INLET_K - 0.2 - BASE_K))
        };
        assert!(residual(low) < 0.0 && residual(high) > 0.0);
        for _ in 0..100 {
            let middle = 0.5 * (low + high);
            if residual(middle) > 0.0 {
                high = middle;
            } else {
                low = middle;
            }
        }
        let h = 0.5 * (low + high);
        let t = BASE_K + h / CP;
        let (htc, dq_dt, fixed_density, dq_dambient) = natural_law(t);
        let diagonal = mass + dt * dq_dt / CP;
        derivative[0] *= mass / diagonal;
        derivative[1] = (mass * derivative[1] + dt) / diagonal;
        derivative[2] = (mass * derivative[2] - VOLUME * (h - old)) / diagonal;
        ambient = (mass * ambient - dt * dq_dambient) / diagonal;
        for previous in &mut intervals[..index] {
            *previous *= mass / diagonal;
        }
        intervals[index] = dt / diagonal;
        frozen_h_initial *= mass / (mass + dt * area * htc / CP);
        frozen_density_initial *= mass / (mass + dt * fixed_density / CP);
        temperatures.push(t);
        old = h;
    }
    for value in &mut derivative {
        *value /= CP;
    }
    for value in &mut intervals {
        *value /= CP;
    }
    NaturalOracle {
        temperatures,
        gradients: derivative,
        intervals,
        frozen_h_initial: frozen_h_initial / CP,
        frozen_density_initial: frozen_density_initial / CP,
        ambient: ambient / CP,
    }
}

#[test]
fn g1_native_natural_enthalpy_adjoint_differentiates_heated_and_cooled_film_properties() {
    let cards = fixture_cards();
    for initial_k in [330.0, 280.0] {
        let input = Inputs {
            initial_h: CP * (initial_k - BASE_K),
            ..Inputs::hot()
        };
        let (spec, bytes) = uniform_project(input, true);
        let expected = natural_oracle(input);
        let solved = run(&spec, &bytes, &cards, Some(FINAL));
        check_final(&solved, 1.0);
        let actual = gradients(solved.adjoint());
        for (found, reference) in actual.into_iter().zip(expected.gradients) {
            close(found, reference, 3e-7);
        }
        assert!((expected.gradients[0] - expected.frozen_h_initial).abs() > 5e-4);
        assert!(
            (expected.gradients[0] - expected.frozen_density_initial).abs() > 1e-5,
            "the oracle must detect omitting density's film-temperature derivative"
        );
        let fine = solved.transient().get("fine").unwrap().as_array().unwrap();
        for (row, &temperature) in fine.iter().zip(&expected.temperatures) {
            close(n(row, "final_region_max_k"), temperature, 5e-8);
        }
        for (interval, derivative) in solved
            .adjoint()
            .get("intervals")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .zip(&expected.intervals)
        {
            close(
                n(
                    region(interval, "regional_power"),
                    "dtemperature_ddelivered_watts_k_per_w",
                ),
                *derivative,
                3e-8,
            );
        }
        let perturb = |delta: f64| {
            let input = Inputs {
                initial_h: input.initial_h + 2.0 * delta,
                power: input.power + 3.0 * delta,
                density: input.density + 0.4 * delta,
            };
            run(&uniform_project(input, true).0, &bytes, &cards, None).final_temperature()
        };
        let epsilon = 0.02;
        close(
            2.0 * actual[0] + 3.0 * actual[1] + 0.4 * actual[2],
            (perturb(epsilon) - perturb(-epsilon)) / (2.0 * epsilon),
            4e-6,
        );
        let reservoir = convection_control(solved.adjoint(), "natural_reservoirs", "target", "air");
        close(n(reservoir, "ambient_temperature_k"), INLET_K, 0.0);
        let ambient = n(reservoir, "dtemperature_dambient_temperature_k_per_k");
        close(ambient, expected.ambient, 3e-7);
        assert!(ambient > 0.0);
        let evaluate_ambient = |delta: f64| {
            let mut changed = spec.clone();
            set_natural_ambient(&mut changed, "air", INLET_K + delta);
            // A reservoir control does not alter the hydraulic operating point.
            assert_eq!(changed.envelope, spec.envelope);
            run(&changed, &bytes, &cards, None).final_temperature()
        };
        close(
            ambient,
            (evaluate_ambient(epsilon) - evaluate_ambient(-epsilon)) / (2.0 * epsilon),
            4e-6,
        );
        if initial_k > INLET_K {
            let peak = run(&spec, &bytes, &cards, Some(PEAK));
            close(n(peak.adjoint(), "state_index"), 0.0, 0.0);
            close(
                n(
                    convection_control(peak.adjoint(), "natural_reservoirs", "target", "air"),
                    "dtemperature_dambient_temperature_k_per_k",
                ),
                0.0,
                0.0,
            );
        }
    }
}

#[test]
fn g1_native_forced_air_enthalpy_controls_match_fresh_initial_power_and_density_resolves() {
    let input = Inputs::hot();
    let (mut spec, bytes) = uniform_project(input, false);
    let speed = 1.2; // Distinguish an absolute-speed derivative from a logarithmic one.
    set_fan_speed(&mut spec, speed);
    let cards = fixture_cards();
    let solved = run(&spec, &bytes, &cards, Some(FINAL));
    check_final(&solved, 1.0);
    let derivative = gradients(solved.adjoint());
    assert!(derivative[0] > 0.0 && derivative[0] < 1.0 / CP);
    assert!(derivative[1] > 0.0 && derivative[2] > 0.0);
    for control in 0..3 {
        let evaluate = |delta: f64| {
            let mut changed = input;
            match control {
                0 => changed.initial_h += delta,
                1 => changed.power += delta,
                _ => changed.density += delta,
            }
            let mut changed_spec = uniform_project(changed, false).0;
            set_fan_speed(&mut changed_spec, speed);
            // Changing initial h cannot re-gauge the constitutive chart.
            assert_eq!(
                time(&mut changed_spec.clone())
                    .enthalpy
                    .as_ref()
                    .unwrap()
                    .materials[0]
                    .knots,
                time(&mut spec.clone()).enthalpy.as_ref().unwrap().materials[0].knots
            );
            run(&changed_spec, &bytes, &cards, None).final_temperature()
        };
        let epsilon = if control == 2 { 0.005 } else { 0.05 };
        close(
            derivative[control],
            (evaluate(epsilon) - evaluate(-epsilon)) / (2.0 * epsilon),
            4e-6,
        );
    }
    let inlet = convection_control(solved.adjoint(), "air_inlets", "branch", "air");
    close(n(inlet, "inlet_temperature_k"), INLET_K, 0.0);
    let inlet_gradient = n(inlet, "dtemperature_dinlet_temperature_k_per_k");
    assert!(inlet_gradient > 0.0);
    let fan = convection_control(solved.adjoint(), "fan_banks", "bank", "fixture-bank");
    assert_eq!(fan.str_field("status"), Some("available"));
    close(n(fan, "speed_ratio"), speed, 0.0);
    let fan_gradient = n(fan, "dtemperature_dspeed_ratio_k");
    assert!(fan_gradient < 0.0);
    for control in ["inlet", "speed"] {
        let evaluate = |delta: f64| {
            let mut changed = spec.clone();
            if control == "inlet" {
                set_air_inlet(&mut changed, INLET_K + delta);
            } else {
                set_fan_speed(&mut changed, speed + delta);
            }
            assert_eq!(changed.envelope, spec.envelope);
            run(&changed, &bytes, &cards, None).final_temperature()
        };
        let (epsilon, expected) = if control == "inlet" {
            (0.02, inlet_gradient)
        } else {
            (0.0005, fan_gradient)
        };
        close(
            expected,
            (evaluate(epsilon) - evaluate(-epsilon)) / (2.0 * epsilon),
            1e-5,
        );
    }
    let peak = run(&spec, &bytes, &cards, Some(PEAK));
    close(n(peak.adjoint(), "time_s"), 0.0, 0.0);
    close(n(peak.adjoint(), "value_k"), 330.0, 1e-12);
    close(n(peak.adjoint(), "state_index"), 0.0, 0.0);
    for (actual, expected) in gradients(peak.adjoint())
        .into_iter()
        .zip([1.0 / CP, 0.0, 0.0])
    {
        close(actual, expected, 1e-12);
    }
    close(
        n(
            convection_control(peak.adjoint(), "air_inlets", "branch", "air"),
            "dtemperature_dinlet_temperature_k_per_k",
        ),
        0.0,
        0.0,
    );
    let peak_fan = convection_control(peak.adjoint(), "fan_banks", "bank", "fixture-bank");
    assert_eq!(peak_fan.str_field("status"), Some("available"));
    close(n(peak_fan, "dtemperature_dspeed_ratio_k"), 0.0, 0.0);

    // The same physical operating point at a declared domain endpoint has
    // no two-sided speed derivative. Keep its bank and explicit reason.
    let mut endpoint_spec = spec.clone();
    endpoint_spec
        .cooling
        .as_mut()
        .unwrap()
        .fan_system
        .as_mut()
        .unwrap()
        .banks[0]
        .speed_ratio_domain
        .0 = speed;
    let endpoint = run(&endpoint_spec, &bytes, &cards, Some(FINAL));
    close(
        endpoint.final_temperature(),
        solved.final_temperature(),
        1e-12,
    );
    let unavailable = convection_control(endpoint.adjoint(), "fan_banks", "bank", "fixture-bank");
    assert_eq!(unavailable.str_field("status"), Some("unsupported"));
    assert!(
        unavailable
            .str_field("reason")
            .unwrap()
            .contains("interior")
    );
    assert!(
        unavailable
            .get("dtemperature_dspeed_ratio_k")
            .is_none_or(|value| matches!(value, J::Null))
    );
}

#[test]
fn g1_native_air_feedback_adjoint_carries_plateau_history_into_later_sensible_cooling() {
    let input = Inputs::hot();
    let (mut spec, bytes) = uniform_project(input, false);
    let cards = fixture_cards();
    let sensible = run(&spec, &bytes, &cards, None);
    let exchange = sensible.receipt.get("conjugate").unwrap();
    let segment = &exchange.get("segments").unwrap().as_array().unwrap()[0];
    let capacity_rate = 1007.0 * n(exchange, "mass_flow_kg_s");
    let conductance = capacity_rate
        * -(-n(segment, "htc_w_m2_k") * n(segment, "wetted_area_m2") / capacity_rate).exp_m1();
    let mass = input.density * VOLUME;
    let dt = 0.25;
    let plateau_loss_h = dt * (conductance * (330.0 - INLET_K) - input.power) / mass;
    assert!(plateau_loss_h > 0.0);
    // Place the phase exit strictly between the first and second fine-grid
    // endpoints. The source chart stays fixed during all initial-h probes.
    let initial_h = 1.5 * plateau_loss_h;
    let latent_h = 4.0 * plateau_loss_h;
    let material = &mut time(&mut spec).enthalpy.as_mut().unwrap().materials[0];
    material.initial_specific_enthalpy.value = initial_h;
    material.knots = [
        (CP * (233.15 - 330.0), 233.15, 0.0),
        (0.0, 330.0, 0.0),
        (latent_h, 330.0, 1.0),
        (latent_h + CP * (398.15 - 330.0), 398.15, 1.0),
    ]
    .into_iter()
    .map(|(h, temperature, fraction)| TransientEnthalpyKnot {
        specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
        temperature: QtyAny::new(temperature, dims::TEMPERATURE),
        liquid_mass_fraction: fraction,
    })
    .collect();
    let solved = run(&spec, &bytes, &cards, Some(FINAL));
    check_final(&solved, 1.0);
    let fine = solved.transient().get("fine").unwrap().as_array().unwrap();
    let mut old_h = initial_h;
    let mut initial_derivative = 1.0;
    let mut inlet_derivative = 0.0;
    for (index, row) in fine.iter().enumerate() {
        let next_on_plateau = old_h - plateau_loss_h;
        let diagonal = if next_on_plateau > 0.0 {
            mass
        } else {
            mass + dt * conductance / CP
        };
        initial_derivative *= mass / diagonal;
        inlet_derivative = (mass * inlet_derivative + dt * conductance) / diagonal;
        let h = mass * next_on_plateau / diagonal;
        let temperature = if h > 0.0 { 330.0 } else { 330.0 + h / CP };
        close(n(row, "final_region_max_k"), temperature, 5e-8);
        close(
            n(row, "mean_liquid_mass_fraction"),
            (h / latent_h).clamp(0.0, 1.0),
            2e-8,
        );
        assert_eq!(h > 0.0, index == 0);
        old_h = h;
    }
    let expected = initial_derivative / CP;
    assert!(
        expected > 0.01,
        "a zero T'(h) on the plateau must not erase its stored history"
    );
    close(gradients(solved.adjoint())[0], expected, 3e-7);
    let evaluate = |delta: f64| {
        let mut changed = spec.clone();
        time(&mut changed).enthalpy.as_mut().unwrap().materials[0]
            .initial_specific_enthalpy
            .value += delta;
        let result = run(&changed, &bytes, &cards, None);
        close(
            n(
                region(result.transient(), "materials"),
                "initial_temperature_k",
            ),
            330.0,
            0.0,
        );
        result.final_temperature()
    };
    let epsilon = (0.05 * plateau_loss_h).min(0.05);
    close(
        expected,
        (evaluate(epsilon) - evaluate(-epsilon)) / (2.0 * epsilon),
        4e-6,
    );
    let inlet = n(
        convection_control(solved.adjoint(), "air_inlets", "branch", "air"),
        "dtemperature_dinlet_temperature_k_per_k",
    );
    close(inlet, inlet_derivative / CP, 3e-7);
    let evaluate_inlet = |delta: f64| {
        let mut changed = spec.clone();
        set_air_inlet(&mut changed, INLET_K + delta);
        run(&changed, &bytes, &cards, None).final_temperature()
    };
    close(
        inlet,
        (evaluate_inlet(0.02) - evaluate_inlet(-0.02)) / 0.04,
        4e-6,
    );
}

fn mixed_project(input: Inputs, card: &str) -> (ProjectSpec, Vec<u8>) {
    let bytes = tetra_stl();
    let mut spec = conjugate_fixture_project(7, &bytes, "convection.gnielinski");
    storage(&mut spec, input, 0.5);
    spec.assignments.as_mut().unwrap()[0].allow_overlap = true;
    for (name, normal) in [
        ("upstream-wall", [1.0, 0.0, 0.0]),
        ("natural-wall", [0.0, 0.0, 1.0]),
    ] {
        spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
            parent: "enclosure".into(),
            name: name.into(),
            display: name.into(),
            expect_id: None,
        });
        spec.assignments.as_mut().unwrap().push(GeometryAssignment {
            artifact: "enclosure".into(),
            target: name.into(),
            length_unit: "m".into(),
            selector: MeshSelector::HalfSpace {
                normal,
                offset: 0.0,
                side: HalfSpaceSide::AtMost,
                tolerance: 0.0,
            },
            allow_overlap: true,
        });
    }
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let first = setup.boundaries[0].condition.clone();
    if let ThermalBoundaryCondition::AirflowConvection { order, .. } =
        &mut setup.boundaries[0].condition
    {
        *order = 1;
    }
    setup.boundaries.push(ThermalBoundary {
        target: "upstream-wall".into(),
        condition: first,
    });
    setup.boundaries.push(ThermalBoundary {
        target: "natural-wall".into(),
        condition: natural_condition(350.0),
    });
    radiation_product::declare(&mut spec, card, "natural-wall", 270.0);
    let radiation = spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .radiation
        .as_mut()
        .unwrap();
    radiation.max_iterations = 24;
    radiation.relaxation = 1.0;
    time(&mut spec)
        .enthalpy
        .as_mut()
        .unwrap()
        .max_newton_iterations = 8;
    (spec, bytes)
}

#[test]
fn g3_native_mixed_enthalpy_adjoint_includes_upstream_air_and_radiative_feedback() {
    let input = Inputs::hot();
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), 0.85);
    let (spec, bytes) = mixed_project(input, &card);
    let solved = run(&spec, &bytes, &cards, Some(FINAL));
    check_final(&solved, 0.5);
    let air = solved.receipt.get("conjugate").unwrap();
    let segments = air.get("segments").unwrap().as_array().unwrap();
    assert_eq!(segments.len(), 2);
    assert_eq!(segments[0].str_field("target"), Some("upstream-wall"));
    assert_eq!(segments[1].str_field("target"), Some("air"));
    close(
        n(&segments[0], "air_out_k"),
        n(&segments[1], "air_in_k"),
        1e-12,
    );
    assert!(n(&segments[1], "air_in_k") > INLET_K + 0.1);
    let minimum = solved
        .temperature
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    assert!(
        solved.final_temperature() - minimum > 0.01,
        "distinct reservoirs must exercise a spatially nonuniform response"
    );
    assert!(
        n(
            solved.receipt.get("radiation").unwrap(),
            "nonlinear_radiative_out_w"
        ) > 0.0
    );
    let natural = solved
        .receipt
        .get("natural_convection")
        .unwrap()
        .get("laws")
        .unwrap()
        .as_array()
        .unwrap();
    assert!(n(&natural[0], "heat_rate_w") < 0.0);
    let derivative = gradients(solved.adjoint());
    let evaluate = |delta: f64| {
        let changed = Inputs {
            initial_h: input.initial_h + 2.0 * delta,
            power: input.power + 3.0 * delta,
            density: input.density + 0.4 * delta,
        };
        run(&mixed_project(changed, &card).0, &bytes, &cards, None).final_temperature()
    };
    let epsilon = 0.02;
    close(
        2.0 * derivative[0] + 3.0 * derivative[1] + 0.4 * derivative[2],
        (evaluate(epsilon) - evaluate(-epsilon)) / (2.0 * epsilon),
        8e-6,
    );
    let gradient = solved.adjoint();
    let reservoir = convection_control(gradient, "natural_reservoirs", "target", "natural-wall");
    let inlet = convection_control(gradient, "air_inlets", "branch", "air");
    let fan = convection_control(gradient, "fan_banks", "bank", "fixture-bank");
    close(n(reservoir, "ambient_temperature_k"), 350.0, 0.0);
    close(n(inlet, "inlet_temperature_k"), INLET_K, 0.0);
    assert_eq!(fan.str_field("status"), Some("available"));
    let radiation = &gradient
        .get("radiation")
        .unwrap()
        .get("surfaces")
        .unwrap()
        .as_array()
        .unwrap()[0];
    assert_eq!(radiation.str_field("target"), Some("natural-wall"));
    let control_direction = 0.5 * n(reservoir, "dtemperature_dambient_temperature_k_per_k")
        + 0.3 * n(inlet, "dtemperature_dinlet_temperature_k_per_k")
        + 0.02 * n(fan, "dtemperature_dspeed_ratio_k")
        + 0.015 * n(radiation, "dtemperature_demissivity_k")
        + 0.4 * n(radiation, "dtemperature_dreservoir_temperature_k_per_k");
    let evaluate_controls = |delta: f64| {
        let (changed_cards, changed_card) =
            radiation_product::emissivity_cards(&fixture_cards(), 0.85 + 0.015 * delta);
        let mut changed = mixed_project(input, &changed_card).0;
        set_natural_ambient(&mut changed, "natural-wall", 350.0 + 0.5 * delta);
        // Every segment of the same branch declares the same external inlet;
        // the second segment's actual inlet must still come from the march.
        set_air_inlet(&mut changed, INLET_K + 0.3 * delta);
        set_fan_speed(&mut changed, 1.0 + 0.02 * delta);
        changed
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .radiation
            .as_mut()
            .unwrap()
            .surfaces[0]
            .reservoir_temperature
            .value = 270.0 + 0.4 * delta;
        assert_eq!(changed.envelope, spec.envelope);
        run(&changed, &bytes, &changed_cards, None).final_temperature()
    };
    close(
        control_direction,
        (evaluate_controls(epsilon) - evaluate_controls(-epsilon)) / (2.0 * epsilon),
        1e-5,
    );
}


/// The admitted forward card supplies Nu; this independent derivative uses
/// Ra proportional to |Tw-Ta|/Tfilm^3 and never calls the adjoint implementation.
fn actual_natural_coefficient(wall: f64) -> (f64, f64) {
    use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity, evaluate};
    let delta = wall - INLET_K;
    let film = INLET_K + 0.5 * delta;
    let density = 101325.0 / (287.05 * film);
    let kinematic = 1.846e-5 / density;
    let rayleigh = 9.80665 * delta.abs() * 0.1_f64.powi(3)
        / (film * kinematic * (kinematic / 0.707));
    let card = evaluate(
        CorrelationId::ChurchillChuVerticalPlate,
        CorrelationInputs::natural(rayleigh, 0.707),
    ).unwrap();
    assert!(card.evidence().model.in_domain);
    let h = card.heat_transfer_coefficient(
        ThermalConductivity::new(0.0263), fs_qty::Length::new(0.1),
    ).unwrap().value.value();
    let exponent = (1.0 - 0.825 / card.evidence().value.sqrt()) / 3.0;
    (h, h * exponent * (1.0 / delta - 1.5 / film))
}

/// Uniform imported-tetra balance with a 0..300 J/kg latent interval at 330 K.
/// Derivative coordinates are initial h, power, density, emissivity, reservoir.
fn natural_radiative_phase_oracle(
    input: Inputs, emissivity: f64, reservoir: f64,
) -> (Vec<[f64; 2]>, [f64; 5], f64) {
    let sigma = 5.670_374_419e-8;
    let area = 8.0 * 3.0_f64.sqrt();
    let mass = input.density * VOLUME;
    let dt = 0.25;
    let state = |h: f64| {
        if h < 0.0 {
            (330.0 + h / CP, 1.0 / CP)
        } else if h < 300.0 {
            (330.0, 0.0)
        } else {
            (330.0 + (h - 300.0) / CP, 1.0 / CP)
        }
    };
    let mut old_h = input.initial_h;
    let mut derivatives = [1.0, 0.0, 0.0, 0.0, 0.0];
    let mut frozen_initial = 1.0;
    let mut endpoints = Vec::new();
    let mut final_slope = 0.0;
    for _ in 0..4 {
        let balance = |h: f64| {
            let (t, _) = state(h);
            let convection = area * actual_natural_coefficient(t).0 * (t - INLET_K);
            let radiation = emissivity * sigma * area * (t.powi(4) - reservoir.powi(4));
            mass * (h - old_h) + dt * (convection + radiation - input.power)
        };
        // All bracket states stay within the source chart and above ambient.
        let (mut low, mut high) = (-100.0, 900.0);
        assert!(balance(low) < 0.0 && balance(high) > 0.0);
        for _ in 0..100 {
            let middle = 0.5 * (low + high);
            if balance(middle) > 0.0 { high = middle; } else { low = middle; }
        }
        let h = 0.5 * (low + high);
        let (t, slope) = state(h);
        let (htc, wall_slope) = actual_natural_coefficient(t);
        let radiation_slope = 4.0 * emissivity * sigma * area * t.powi(3);
        let convection_slope = area * (htc + (t - INLET_K) * wall_slope);
        let diagonal = mass + dt * (convection_slope + radiation_slope) * slope;
        let forcing = [
            0.0,
            dt,
            -VOLUME * (h - old_h),
            -dt * sigma * area * (t.powi(4) - reservoir.powi(4)),
            4.0 * dt * emissivity * sigma * area * reservoir.powi(3),
        ];
        for (derivative, rhs) in derivatives.iter_mut().zip(forcing) {
            *derivative = (mass * *derivative + rhs) / diagonal;
        }
        // Omit only natural h'(T), retaining the exact radiation derivative.
        frozen_initial *= mass / (mass + dt * (area * htc + radiation_slope) * slope);
        endpoints.push([h, t]);
        old_h = h;
        final_slope = slope;
    }
    for derivative in &mut derivatives { *derivative *= final_slope; }
    (endpoints, derivatives, frozen_initial * final_slope)
}

#[test]
fn g1_native_natural_radiative_reverse_retains_latent_history_and_surface_controls() {
    let input = Inputs { initial_h: 80.0, ..Inputs::hot() };
    let run_phase = |input: Inputs, emissivity: f64, reservoir: f64, output: Option<&str>| {
        // A fresh immutable pack and binding are built for every emissivity probe.
        let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), emissivity);
        let (mut spec, bytes) = uniform_project(input, true);
        let policy = time(&mut spec).enthalpy.as_mut().unwrap();
        policy.max_newton_iterations = 8;
        policy.max_backtracks = 12;
        policy.materials[0].knots = [
            (CP * (233.15 - 330.0), 233.15, 0.0),
            (0.0, 330.0, 0.0),
            (300.0, 330.0, 1.0),
            (300.0 + CP * (398.15 - 330.0), 398.15, 1.0),
        ].into_iter().map(|(h, temperature, fraction)| TransientEnthalpyKnot {
            specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
            temperature: QtyAny::new(temperature, dims::TEMPERATURE),
            liquid_mass_fraction: fraction,
        }).collect();
        radiation_product::declare(&mut spec, &card, "air", reservoir);
        let radiation = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .radiation.as_mut().unwrap();
        radiation.surfaces[0].query_temperature.value = 330.0;
        radiation.max_iterations = 24;
        radiation.relaxation = 1.0;
        run(&spec, &bytes, &cards, output)
    };
    let surface = |receipt: &J| -> J {
        let surfaces = receipt.get("radiation").unwrap().get("surfaces").unwrap()
            .as_array().unwrap();
        assert_eq!(surfaces.len(), 1);
        surfaces[0].clone()
    };
    let (expected_states, expected, frozen_initial) =
        natural_radiative_phase_oracle(input, 0.85, 400.0);
    let solved = run_phase(input, 0.85, 400.0, Some(FINAL));
    check_final(&solved, 1.0);
    assert_eq!(solved.temperature.len(), 4, "the balance assumes one regular tetra");
    for (row, &[h, temperature]) in solved.transient().get("fine").unwrap()
        .as_array().unwrap().iter().zip(&expected_states)
    {
        close(n(row, "final_region_max_k"), temperature, 5e-8);
        close(n(row, "mean_liquid_mass_fraction"), (h / 300.0).clamp(0.0, 1.0), 2e-8);
        let natural = row.get("natural").unwrap();
        assert!(n(natural, "physical_residual_norm_j")
            <= n(natural, "physical_residual_tolerance_j"));
        close(n(natural.get("physical_coefficients_w_m2_k").unwrap(), "air"),
            actual_natural_coefficient(temperature).0, 2e-8);
    }
    assert!(expected_states[..3].iter().all(|[h, t]| *h > 0.0 && *h < 300.0 && *t == 330.0));
    assert!(expected_states[3][0] > 300.0);
    assert!(expected[0] > 0.08, "stored latent history must reach the final sensible state");
    assert!((expected[0] - frozen_initial).abs() > 1e-3);
    for (&actual, &reference) in gradients(solved.adjoint()).iter().zip(&expected[..3]) {
        close(actual, reference, 3e-7);
    }
    let controls = surface(solved.adjoint());
    assert_eq!(controls.str_field("target"), Some("air"));
    close(n(&controls, "emissivity"), 0.85, 0.0);
    close(n(&controls, "reservoir_temperature_k"), 400.0, 0.0);
    let emissivity = n(&controls, "dtemperature_demissivity_k");
    let reservoir = n(&controls, "dtemperature_dreservoir_temperature_k_per_k");
    close(emissivity, expected[3], 2e-6);
    close(reservoir, expected[4], 2e-7);
    let epsilon = 0.02;
    let material_probe = |delta: f64| {
        let changed = Inputs {
            initial_h: input.initial_h + 2.0 * delta,
            power: input.power + 3.0 * delta,
            density: input.density + 0.4 * delta,
        };
        let trial = run_phase(changed, 0.85, 400.0, None);
        close(n(region(trial.transient(), "materials"), "initial_temperature_k"), 330.0, 0.0);
        trial.final_temperature()
    };
    let actual = gradients(solved.adjoint());
    close(2.0 * actual[0] + 3.0 * actual[1] + 0.4 * actual[2],
        (material_probe(epsilon) - material_probe(-epsilon)) / (2.0 * epsilon), 4e-6);
    let source_card = surface(&solved.receipt).str_field("card").unwrap().to_string();
    let radiation_probe = |delta: f64| {
        let changed_emissivity = 0.85 + 0.02 * delta;
        let changed_reservoir = 400.0 + 0.5 * delta;
        let trial = run_phase(input, changed_emissivity, changed_reservoir, None);
        let retained = surface(&trial.receipt);
        assert_ne!(retained.str_field("card").unwrap(), source_card);
        close(n(&retained, "emissivity"), changed_emissivity, 0.0);
        close(n(&retained, "ambient_temperature_k"), changed_reservoir, 0.0);
        trial.final_temperature()
    };
    close(0.02 * emissivity + 0.5 * reservoir,
        (radiation_probe(epsilon) - radiation_probe(-epsilon)) / (2.0 * epsilon), 4e-6);
}
