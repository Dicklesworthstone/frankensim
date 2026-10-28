//! Closed exchange without an artificial convection anchor. The equal-area
//! two-patch factors are a manufactured algebraic closure, not CAD view factors.
use super::*;
use fs_project::TransientInitialTemperature;

const SECOND_EPSILON: f64 = 0.6;

#[derive(Clone, Copy, Debug)]
enum Storage {
    Capacity,
    Enthalpy,
}

fn fixture(storage: Storage) -> (ProjectSpec, CardPackSet) {
    let (cards, first) = radiation_product::emissivity_cards(&fixture_cards(), EPSILON);
    let (cards, second) = radiation_product::emissivity_cards(&cards, SECOND_EPSILON);
    let mut spec = project(&first);
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.enclosure_radiation.as_mut().unwrap().surfaces[1].card = second;
    for boundary in &mut setup.boundaries {
        boundary.condition = flux(0.0);
    }
    let policy = time(&mut spec);
    match storage {
        Storage::Capacity => {
            policy.enthalpy = None;
            policy.capacities = vec![TransientRegionCapacity {
                region: "air".into(),
                volumetric_heat_capacity: QtyAny::new(DENSITY * CP, dims::VOLUMETRIC_HEAT_CAPACITY),
                source: "synthetic capacity for a closed enclosure without convection".into(),
            }];
            policy.initial_temperature = TransientInitialTemperature::Uniform(
                QtyAny::new(INITIAL_K, dims::TEMPERATURE));
        }
        Storage::Enthalpy => {
            policy.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
        }
    }
    (spec, cards)
}

fn flux(outward_w_m2: f64) -> ThermalBoundaryCondition {
    ThermalBoundaryCondition::HeatFlux {
        outward_flux: QtyAny::new(outward_w_m2, dims::HEAT_FLUX),
    }
}

fn maximum(field: &[f64]) -> f64 {
    field.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn endpoint(spec: &ProjectSpec, cards: &CardPackSet) -> f64 {
    maximum(&values(&J::parse(&run(spec, cards).field).unwrap(), "temperature"))
}

fn derivative_close(actual: f64, expected: f64) {
    close(actual, expected, 5e-6 + 3e-4 * expected.abs());
}

fn assert_exchange(receipt: &J, nonzero: bool) {
    let surfaces = receipt.get("radiation").unwrap().get("surfaces").unwrap()
        .as_array().unwrap();
    let a = n(&surfaces[0], "mean_temperature_k");
    let b = n(&surfaces[1], "mean_temperature_k");
    let exchange = SIGMA * 0.5 * (a - b) * (a + b) * (a * a + b * b)
        / (1.0 / EPSILON + 1.0 / SECOND_EPSILON - 1.0);
    if nonzero {
        assert!(exchange.abs() > 1e-5, "the fixture must exercise reflected internal heat");
    }
    close(n(&surfaces[0], "nonlinear_heat_w"), exchange, 1e-6);
    close(n(&surfaces[1], "nonlinear_heat_w"), -exchange, 1e-6);
}

/// The external Neumann account is independent of the applied radiosity
/// iteration. A mixed Robin/flux enclosure has nonzero internal exchange in
/// each subset even though the complete enclosure's net radiation is zero.
fn check_energy(solved: &Solved, prescribed_w: f64, dynamic: Option<&str>) {
    let receipt = J::parse(&solved.conduction).unwrap();
    let history = receipt.get("transient").unwrap();
    let energy = receipt.get("energy").unwrap();
    close(n(energy, "neumann_out_w"), prescribed_w, 1e-12);
    let mut fine_storage = 0.0;
    for grid in ["coarse", "fine"] {
        let mut storage = 0.0;
        let mut input = 0.0;
        for row in history.get(grid).unwrap().as_array().unwrap() {
            let radiation = row.get("radiation").unwrap();
            close(n(radiation, "prescribed_neumann_out_w"), prescribed_w, 1e-12);
            close(n(radiation, "nonlinear_radiation_out_w"), 0.0, 1e-6);
            assert!(n(radiation, "decomposition_residual_w").abs() <= 1e-6);
            assert!(n(radiation, "physical_residual_j")
                <= n(radiation, "physical_residual_threshold_j"));
            let external_convection = match dynamic {
                Some("airflow") => {
                    let air = row.get("conjugate").unwrap();
                    assert!(n(air, "physical_residual_norm_j")
                        <= n(air, "physical_residual_tolerance_j"));
                    assert!(n(air, "physical_energy_residual_j").abs()
                        <= n(history, "energy_tolerance_j"));
                    assert!(n(air, "air_heat_gain_w") > 0.0);
                    // The other surface is Neumann. Radiation on it cannot
                    // create a convection path or feed transported-air heat.
                    close(n(air, "off_path_convective_out_w"), 0.0, 1e-12);
                    let net = n(row, "source_w") - prescribed_w
                        - n(air, "air_heat_gain_w")
                        - n(air, "off_path_convective_out_w")
                        - n(air, "radiative_out_w")
                        + n(air, "physical_dirichlet_in_w");
                    close(n(row, "stored_energy_change_j") - n(row, "dt_s") * net,
                        n(air, "coupled_energy_residual_j"), 2e-8);
                    n(air, "air_heat_gain_w") + n(air, "off_path_convective_out_w")
                }
                Some("natural") => {
                    let natural = row.get("natural").unwrap();
                    assert!(n(natural, "physical_residual_norm_j")
                        <= n(natural, "physical_residual_tolerance_j"));
                    assert!(n(natural, "physical_energy_residual_j").abs()
                        <= n(history, "energy_tolerance_j"));
                    let fluxes = natural.get("physical_convective_robin_fluxes").unwrap()
                        .as_array().unwrap();
                    assert_eq!(fluxes.len(), 1);
                    assert_eq!(fluxes[0].str_field("region"), Some("wall-z"));
                    n(natural, "physical_convective_out_w")
                }
                None => n(radiation, "convective_out_w"),
                Some(_) => panic!("unknown boundary test"),
            };
            close(n(row, "net_input_w"),
                n(row, "source_w") - prescribed_w - external_convection, 2e-6);
            let dt = n(row, "dt_s");
            close(n(row, "stored_energy_change_j") - dt * n(row, "net_input_w"),
                n(row, "energy_residual_j"), 2e-9);
            assert!(n(row, "energy_residual_j").abs() <= n(history, "energy_tolerance_j"));
            storage += n(row, "stored_energy_change_j");
            input += dt * n(row, "net_input_w");
        }
        close(storage, input, 1e-5);
        if grid == "fine" { fine_storage = storage; }
    }
    let temperature = values(&J::parse(&solved.field).unwrap(), "temperature");
    assert_eq!(temperature.len(), 4);
    let independently_stored: f64 = temperature.iter()
        .map(|temperature| DENSITY * CP / 24.0 * (temperature - INITIAL_K)).sum();
    close(fine_storage, independently_stored, 2e-8);
    close(n(history.get("energy").unwrap(), "stored_change_j"), fine_storage, 2e-9);
}

#[test]
fn g1_native_vacuum_enclosure_stores_known_input_without_a_convection_anchor() {
    for storage in [Storage::Capacity, Storage::Enthalpy] {
        let (spec, cards) = fixture(storage);
        let solved = run(&spec, &cards);
        check_energy(&solved, 0.0, None);
        let receipt = J::parse(&solved.conduction).unwrap();
        assert_exchange(&receipt, false);
        let history = receipt.get("transient").unwrap();
        // Unit tetrahedron volume = 1/6, declared rho*Cp = 60; total
        // capacity is 10 J/K. The event-exact pulse delivers exactly 9 J.
        close(n(history.get("energy").unwrap(), "stored_change_j"), 9.0, 1e-5);
        close(n(receipt.get("energy").unwrap(), "robin_out_w"), 0.0, 1e-6);
        for temperature in values(&J::parse(&solved.field).unwrap(), "temperature") {
            close(temperature, INITIAL_K + 0.9, 2e-6);
        }
        for row in history.get("fine").unwrap().as_array().unwrap() {
            close(n(row.get("radiation").unwrap(), "convective_out_w"), 0.0, 0.0);
        }
    }
}

#[test]
fn g1_native_mixed_flux_enclosure_preserves_external_heat_through_dynamic_convection() {
    for storage in [Storage::Capacity, Storage::Enthalpy] {
        for kind in ["prescribed", "natural", "airflow"] {
            let (mut spec, cards) = fixture(storage);
            spec.envelope.as_mut().unwrap().ambient_lo.value = 290.0;
            spec.envelope.as_mut().unwrap().ambient_hi.value = 310.0;
            let boundaries = &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries;
            boundaries[0].condition = flux(-12.0);
            boundaries[1].condition = match kind {
                "prescribed" => ThermalBoundaryCondition::Convection {
                    coefficient: QtyAny::new(10.0, dims::HEAT_TRANSFER_COEFFICIENT),
                    reference_temperature: QtyAny::new(290.0, dims::TEMPERATURE),
                },
                "natural" => ThermalBoundaryCondition::NaturalConvection {
                    characteristic_length: QtyAny::new(0.1, dims::LENGTH),
                    ambient_temperature: QtyAny::new(290.0, dims::TEMPERATURE),
                    correlation: "convection.churchill-chu-vertical-plate".into(),
                },
                "airflow" => ThermalBoundaryCondition::AirflowConvection {
                    branch: "air".into(), order: 0,
                    inlet_temperature: QtyAny::new(290.0, dims::TEMPERATURE),
                    hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH),
                    flow_area: QtyAny::new(0.004, dims::AREA),
                    channel_length: QtyAny::new(0.3, dims::LENGTH),
                    correlation: "convection.gnielinski".into(),
                },
                _ => unreachable!(),
            };
            let solved = run(&spec, &cards);
            check_energy(&solved, -6.0, (kind != "prescribed").then_some(kind));
            assert_exchange(&J::parse(&solved.conduction).unwrap(), true);
        }
    }
}

#[test]
fn g1_native_flux_only_enclosure_history_derivatives_use_complete_reflection_and_storage() {
    for storage in [Storage::Capacity, Storage::Enthalpy] {
        let (mut spec, cards) = fixture(storage);
        let boundaries = &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries;
        boundaries[0].condition = flux(-12.0);
        boundaries[1].condition = flux(2.0);
        let mut differentiated = spec.clone();
        differentiated.outputs.as_mut().unwrap().push(OutputRequest {
            name: "temperature-final-adjoint".into(), kind: "report".into(), region: None,
        });
        let solved = run(&differentiated, &cards);
        check_energy(&solved, -5.0, None);
        let receipt = J::parse(&solved.conduction).unwrap();
        assert_exchange(&receipt, true);
        let history = receipt.get("transient").unwrap();
        // The authored Neumann input adds 5 J to the 9 J volume pulse.
        close(n(history.get("energy").unwrap(), "stored_change_j"), 14.0, 1e-5);
        let gradient = history.get("adjoint").unwrap();
        close(n(gradient, "active_ties"), 1.0, 0.0);
        assert!(n(gradient, "maximum_dual_relative_residual") <= 1e-9);
        assert!(n(gradient, "adjoint_krylov_iterations")
            <= n(gradient, "max_adjoint_krylov_iterations"));
        let radiation = gradient.get("radiation").unwrap();
        let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
        let emissivity = n(&surfaces[0], "dtemperature_demissivity_k");
        assert!(emissivity.abs() > 1e-8, "this history must exercise reflected emissivity sensitivity");
        let epsilon_step = 0.0005;
        let epsilon_probe = |delta: f64| {
            let (trial_cards, card) = radiation_product::emissivity_cards(&cards, EPSILON + delta);
            let mut changed = spec.clone();
            changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .enclosure_radiation.as_mut().unwrap().surfaces[0].card = card;
            endpoint(&changed, &trial_cards)
        };
        derivative_close(emissivity,
            (epsilon_probe(epsilon_step) - epsilon_probe(-epsilon_step)) / (2.0 * epsilon_step));

        let initial = match storage {
            Storage::Capacity => n(gradient, "dtemperature_duniform_initial_k_per_k"),
            Storage::Enthalpy => CP * n(gradient,
                "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"),
        };
        let initial_probe = |delta: f64| {
            let mut changed = spec.clone();
            match storage {
                Storage::Capacity => time(&mut changed).initial_temperature =
                    TransientInitialTemperature::Uniform(QtyAny::new(INITIAL_K + delta, dims::TEMPERATURE)),
                Storage::Enthalpy => time(&mut changed).enthalpy.as_mut().unwrap()
                    .materials[0].initial_specific_enthalpy.value += CP * delta,
            }
            endpoint(&changed, &cards)
        };
        let step = 0.002;
        derivative_close(initial, (initial_probe(step) - initial_probe(-step)) / (2.0 * step));

        let first_pulse: f64 = gradient.get("intervals").unwrap().as_array().unwrap().iter()
            .filter(|row| n(row, "end_s") <= 0.3)
            .map(|row| n(&row.get("regional_power").unwrap().as_array().unwrap()[0],
                "dtemperature_ddelivered_watts_k_per_w")).sum();
        let pulse_probe = |delta: f64| {
            let mut changed = spec.clone();
            time(&mut changed).power_schedules[0].steps[0].watts.value += delta;
            endpoint(&changed, &cards)
        };
        derivative_close(first_pulse, (pulse_probe(step) - pulse_probe(-step)) / (2.0 * step));
    }
}

