//! Pointwise emission shares the native enthalpy trajectory transpose.
//! The uniform tetrahedron has an independent scalar backward-Euler oracle.
use super::*;

fn pointwise(card: &str, reservoir: f64) -> ProjectSpec {
    let mut spec = project(card, reservoir);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let radiation = setup.radiation.as_mut().unwrap();
    radiation.spatial_model = Some(fs_project::RadiationSpatialModel::Pointwise);
    radiation.max_iterations = 1;
    let policy = time(&mut spec);
    policy.initial_temperature = fs_project::TransientInitialTemperature::FromRegionalEnthalpy;
    policy.power_schedules = vec![TransientRegionPower {
        region: "air".into(), source: "synthetic constant delivered-power adjoint oracle".into(),
        steps: vec![TransientPowerStep {
            until: policy.horizon, watts: QtyAny::new(5.0, dims::POWER),
        }],
    }];
    spec
}

fn request(spec: &mut ProjectSpec, name: &str) {
    spec.outputs.as_mut().unwrap().push(OutputRequest {
        name: name.into(), kind: "report".into(), region: None,
    });
}

fn derivative_close(actual: f64, expected: f64) {
    close(actual, expected, 3e-7 + 3e-5 * expected.abs());
}

#[test]
fn g1_pointwise_ambient_native_dual_matches_independent_storage_and_radiation_recurrences() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    for reservoir in [270.0, 400.0] {
        let spec = pointwise(&card, reservoir);
        let baseline = run(&spec, &cards);
        let mut differentiated = spec.clone();
        request(&mut differentiated, "temperature-final-adjoint");
        let solved = run(&differentiated, &cards);
        assert_eq!(solved.temperature, baseline.temperature);
        assert_eq!(solved.receipt.get("energy"), baseline.receipt.get("energy"));
        assert_eq!(solved.receipt.get("radiation"), baseline.receipt.get("radiation"));
        let history = solved.transient();
        let gradient = history.get("adjoint").unwrap();
        let radiation = gradient.get("radiation").unwrap();
        assert!(radiation.str_field("method").unwrap().contains("pointwise"));
        assert!(n(gradient, "maximum_dual_relative_residual") <= 1e-9);
        assert!(n(gradient, "adjoint_krylov_iterations")
            <= n(gradient, "max_adjoint_krylov_iterations"));
        close(n(gradient, "state_index"), n(history, "fine_steps"), 0.0);
        let surface = &radiation.get("surfaces").unwrap().as_array().unwrap()[0];
        let convection = &gradient.get("prescribed_convection").unwrap()
            .get("boundaries").unwrap().as_array().unwrap()[0];
        let density = &gradient.get("regional_storage").unwrap().as_array().unwrap()[0];
        let power = &gradient.get("regional_power").unwrap().as_array().unwrap()[0];

        // Derivatives of M Cp (T-new - T-old)/dt + Ah(T-Tc)
        // + A epsilon sigma (T^4-Ta^4) - Q = 0. No FEM, radiation or
        // adjoint assembly is used to form this seven-control oracle.
        let area = area_m2();
        let mass = mass_kg(DENSITY);
        let volume = mass / DENSITY;
        let mut old = INITIAL_K;
        let mut derivatives = [0.0; 7];
        derivatives[2] = 1.0 / CP;
        for row in solved.rows("fine") {
            let dt = n(row, "dt_s");
            let temperature = uniform_step(old, dt, reservoir, 5.0);
            let storage = mass * CP / dt;
            let denominator = storage + area * 10.0
                + 4.0 * area * EMISSIVITY * SIGMA * temperature.powi(3);
            let forcing = [
                -area * SIGMA * (temperature.powi(4) - reservoir.powi(4)),
                4.0 * area * EMISSIVITY * SIGMA * reservoir.powi(3),
                0.0,
                -volume * CP * (temperature - old) / dt,
                -area * (temperature - AMBIENT_K),
                area * 10.0,
                1.0,
            ];
            for (value, forcing) in derivatives.iter_mut().zip(forcing) {
                *value = (storage * *value + forcing) / denominator;
            }
            old = temperature;
        }
        let actual = [
            n(surface, "dtemperature_demissivity_k"),
            n(surface, "dtemperature_dreservoir_temperature_k_per_k"),
            n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"),
            n(density, "dtemperature_dreference_density_k_m3_per_kg"),
            n(convection, "dtemperature_dcoefficient_k_per_w_m2_k"),
            n(convection, "dtemperature_dreference_temperature_k_per_k"),
            n(power, "dtemperature_duniform_added_watts_k_per_w"),
        ];
        for (actual, expected) in actual.into_iter().zip(derivatives) {
            derivative_close(actual, expected);
        }
        let vertex = n(gradient, "active_vertex") as usize;
        close(solved.temperature[vertex], old, 2e-8);
        for (coordinate, step) in [(0, 0.0002), (1, 0.002)] {
            let probe = |delta: f64| {
                let mut changed = spec.clone();
                if coordinate == 0 {
                    let (trial, identity) = radiation_product::emissivity_cards(
                        &cards, EMISSIVITY + delta);
                    changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .radiation.as_mut().unwrap().surfaces[0].card = identity;
                    run(&changed, &trial).temperature[vertex]
                } else {
                    changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                        .radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value += delta;
                    run(&changed, &cards).temperature[vertex]
                }
            };
            derivative_close(actual[coordinate], (probe(step) - probe(-step)) / (2.0 * step));
        }
    }
}

#[test]
fn g1_pointwise_ambient_initial_peak_keeps_future_controls_zero() {
    let (cards, card) = radiation_product::emissivity_cards(&fixture_cards(), EMISSIVITY);
    let mut spec = pointwise(&card, 270.0);
    request(&mut spec, "temperature-sampled-peak-adjoint");
    let solved = run(&spec, &cards);
    let gradient = solved.transient().get("adjoint").unwrap();
    close(n(gradient, "state_index"), 0.0, 0.0);
    close(n(gradient, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(n(gradient, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"),
        1.0 / CP, 1e-14);
    close(n(gradient, "dtemperature_dreference_density_multiplier_k"), 0.0, 0.0);
    for row in gradient.get("radiation").unwrap().get("surfaces").unwrap().as_array().unwrap() {
        close(n(row, "dtemperature_demissivity_k"), 0.0, 0.0);
        close(n(row, "dtemperature_dreservoir_temperature_k_per_k"), 0.0, 0.0);
    }
    for row in gradient.get("regional_power").unwrap().as_array().unwrap() {
        close(n(row, "dtemperature_duniform_added_watts_k_per_w"), 0.0, 0.0);
    }
}
