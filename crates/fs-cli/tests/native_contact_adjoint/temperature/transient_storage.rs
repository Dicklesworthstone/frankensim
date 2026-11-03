//! Actual CLI contact trajectories, storage equivalence and physical heat exchange.
use super::*;
use fs_project::{
    ConductionTransient, ThermalBoundaryCondition, TransientBoundarySchedule,
    TransientBoundaryStep, TransientBoundaryTarget, TransientEnthalpy, TransientEnthalpyKnot,
    TransientInitialTemperature, TransientRegionCapacity, TransientRegionEnthalpy, spec::dims,
};
use fs_qty::QtyAny;

const INITIAL_K: f64 = 293.15;
const CP: f64 = 10.0;
const DENSITY: f64 = 10.0;
const ENERGY_TOLERANCE_J: f64 = 1e-6;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

fn time(f: &mut Fixture) -> &mut ConductionTransient {
    f.project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap()
}

fn fixture(enthalpy: bool, scheduled: bool) -> Fixture {
    let mut f = Fixture::new(false, false);
    f.project.power.as_mut().unwrap()[0].watts.value = 10.0;
    f.project.power.as_mut().unwrap()[0].duty = 1.0;
    let bindings = f.project.materials.as_ref().unwrap().clone();
    let capacities = bindings
        .iter()
        .map(|binding| TransientRegionCapacity {
            region: binding.region.clone(),
            volumetric_heat_capacity: QtyAny::new(DENSITY * CP, dims::VOLUMETRIC_HEAT_CAPACITY),
            source:
                "synthetic Cv = 100 J/(m3 K); numerical storage equivalence, not experimental data"
                    .into(),
        })
        .collect();
    let storage = enthalpy.then(|| TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: bindings.iter().map(|binding| TransientRegionEnthalpy {
            region: binding.region.clone(),
            material_card_identity: fs_blake3::ContentHash::from_hex(&binding.card).unwrap(),
            source: "synthetic Cp = 10 J/(kg K), density = 10 kg/m3; same storage as capacity fixture".into(),
            reference_density: QtyAny::new(DENSITY, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            phase_conductivity: None,
            knots: [253.15, 393.15].into_iter().map(|temperature| TransientEnthalpyKnot {
                specific_enthalpy: QtyAny::new(CP * (temperature - INITIAL_K), dims::SPECIFIC_ENTHALPY),
                temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                liquid_mass_fraction: 0.0,
            }).collect(),
        }).collect(),
    });
    let setup = f
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    for boundary in &mut setup.boundaries {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
    }
    setup.transient = Some(ConductionTransient {
        initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(
            INITIAL_K,
            dims::TEMPERATURE,
        )),
        horizon: QtyAny::new(1.0, dims::TIME),
        max_step: QtyAny::new(0.6, dims::TIME),
        max_steps: 600,
        energy_tolerance: QtyAny::new(ENERGY_TOLERANCE_J, dims::ENERGY),
        capacities: if enthalpy { Vec::new() } else { capacities },
        enthalpy: storage,
        power_schedules: Vec::new(),
        boundary_schedules: if scheduled {
            vec![TransientBoundarySchedule {
                target: "cold".into(),
                kind: TransientBoundaryTarget::HeatFlux,
                source: "synthetic held boundary pulse with independently known zero net energy"
                    .into(),
                steps: [(0.25, -2.0), (0.75, 1.0), (1.0, 0.0)]
                    .into_iter()
                    .map(|(until, flux)| TransientBoundaryStep {
                        until: QtyAny::new(until, dims::TIME),
                        value: QtyAny::new(flux, dims::HEAT_FLUX),
                    })
                    .collect(),
            }]
        } else {
            Vec::new()
        },
    });
    f
}

fn solve(f: &Fixture, ordinal: usize) -> (J, J) {
    let (source, contact, ledger) = prepare(f, ordinal, None);
    let result = run(f, &source, &contact, &ledger);
    assert_eq!(
        result.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        result.stdout,
        result.stderr
    );
    let result = J::parse(&result.stdout).unwrap();
    let run_receipt = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = run_receipt
        .get("stages")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row.str_field("stage") == Some("conduction"))
        .unwrap();
    let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
    (receipt, field)
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

fn check_contact(row: &J, selected_claim: &str) -> f64 {
    let contact = row.get("temperature_contact").unwrap();
    assert!(matches!(
        contact.str_field("method"),
        Some("source-backed-contact-newton" | "source-backed-contact-enthalpy-newton")
    ));
    assert!(contact.f64_field("thermal_responses").unwrap() >= 1.0);
    assert!(
        contact.f64_field("physical_residual_norm_j").unwrap()
            <= contact.f64_field("physical_residual_tolerance_j").unwrap()
    );
    assert!(
        contact
            .f64_field("physical_energy_residual_j")
            .unwrap()
            .abs()
            <= ENERGY_TOLERANCE_J
    );
    let faces = contact.get("face_states").unwrap().as_array().unwrap();
    assert!(!faces.is_empty());
    for face in faces {
        let a = face.f64_field("side_a_temperature_k").unwrap();
        let b = face.f64_field("side_b_temperature_k").unwrap();
        let query = face.f64_field("query_temperature_k").unwrap();
        assert_eq!(query, 0.5 * a + 0.5 * b);
        // Independent material law supplied by temperature_pack(), not the solver.
        close(
            face.f64_field("resistance_m2_k_w").unwrap(),
            0.001 * query,
            1e-14,
        );
        assert_eq!(face.str_field("selected_claim"), Some(selected_claim));
        assert_eq!(face.str_field("usage_receipt").unwrap().len(), 64);
    }
    faces[0].f64_field("query_temperature_k").unwrap()
}

fn interface_heat(receipt: &J) -> f64 {
    receipt
        .get("interfaces")
        .unwrap()
        .get("fluxes")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|row| row.f64_field("heat_rate_a_to_b_w").unwrap())
        .sum()
}

#[test]
fn contact_capacity_and_enthalpy_share_boundary_events_and_the_same_physical_storage() {
    let claim = temperature_pack().2;
    let mut endpoints = Vec::new();
    for enthalpy in [false, true] {
        let f = fixture(enthalpy, true);
        let (receipt, field) = solve(&f, 40);
        let transient = receipt.get("transient").unwrap();
        assert_eq!(
            transient
                .get("temperature_contact")
                .unwrap()
                .str_field("method"),
            Some(if enthalpy {
                "source-backed-contact-enthalpy-newton"
            } else {
                "source-backed-contact-newton"
            })
        );
        for (grid, expected_ends) in [
            ("coarse", vec![0.25, 0.75, 1.0]),
            ("fine", vec![0.125, 0.25, 0.5, 0.75, 0.875, 1.0]),
        ] {
            let rows = transient.get(grid).unwrap().as_array().unwrap();
            assert_eq!(rows.len(), expected_ends.len());
            let mut start = 0.0;
            let mut integrated = 0.0;
            let mut queries = Vec::new();
            for (row, end) in rows.iter().zip(expected_ends) {
                close(row.f64_field("time_s").unwrap(), end, 0.0);
                close(row.f64_field("dt_s").unwrap(), end - start, 0.0);
                let outward_flux = if start < 0.25 {
                    -2.0
                } else if start < 0.75 {
                    1.0
                } else {
                    0.0
                };
                // Unit cubes have five external square-metre faces each.
                // The shared contact face is internal and supplies no net heat.
                let expected_net = 10.0 - 5.0 * outward_flux;
                close(row.f64_field("source_w").unwrap(), 10.0, 1e-10);
                close(row.f64_field("net_input_w").unwrap(), expected_net, 1e-10);
                integrated += (end - start) * expected_net;
                queries.push(check_contact(row, &claim));
                start = end;
            }
            close(integrated, 10.0, 0.0);
            assert!(
                (queries.last().unwrap() - queries[0]).abs() > 1e-8,
                "the source contact query must follow changing physical temperatures"
            );
        }
        let energy = transient.get("energy").unwrap();
        close(
            energy.f64_field("integrated_net_input_j").unwrap(),
            10.0,
            1e-9,
        );
        close(
            energy.f64_field("stored_change_j").unwrap(),
            10.0,
            6.0 * ENERGY_TOLERANCE_J,
        );
        assert!(energy.f64_field("maximum_step_residual_j").unwrap() <= ENERGY_TOLERANCE_J);
        assert!(
            interface_heat(&receipt) < 0.0,
            "hot-to-cold heat must survive final report assembly"
        );
        let temperature = values(&field, "temperature");
        assert!(temperature.iter().any(|&value| value > INITIAL_K));
        endpoints.push(temperature);
    }
    assert_eq!(endpoints[0].len(), endpoints[1].len());
    for (&capacity, &enthalpy) in endpoints[0].iter().zip(&endpoints[1]) {
        close(capacity, enthalpy, 2e-7);
    }
}

#[test]
fn unequal_initial_enthalpy_exchanges_heat_through_the_same_source_contact() {
    let mut f = fixture(true, false);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    let transient = time(&mut f);
    transient.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    for material in &mut transient.enthalpy.as_mut().unwrap().materials {
        let (initial, density) = if material.region == "cold" {
            (300.0, 10.0)
        } else {
            (350.0, 20.0)
        };
        material.reference_density.value = density;
        material.initial_specific_enthalpy.value = CP * (initial - INITIAL_K);
    }
    let (receipt, field) = solve(&f, 41);
    let transient = receipt.get("transient").unwrap();
    assert_eq!(
        transient.str_field("initial_state"),
        Some("regional-enthalpy")
    );
    let claim = temperature_pack().2;
    for grid in ["coarse", "fine"] {
        for row in transient.get(grid).unwrap().as_array().unwrap() {
            check_contact(row, &claim);
            close(row.f64_field("net_input_w").unwrap(), 0.0, 1e-12);
            close(row.f64_field("source_w").unwrap(), 0.0, 0.0);
        }
    }
    let energy = transient.get("energy").unwrap();
    close(
        energy.f64_field("stored_change_j").unwrap(),
        0.0,
        4.0 * ENERGY_TOLERANCE_J,
    );
    close(
        energy.f64_field("integrated_net_input_j").unwrap(),
        0.0,
        1e-12,
    );
    let endpoint = transient.get("endpoint").unwrap();
    let masses = values(endpoint, "reference_nodal_mass_kg");
    let enthalpies = values(endpoint, "specific_enthalpy_j_kg");
    let temperature = values(&field, "temperature");
    assert_eq!(masses.len(), enthalpies.len());
    assert_eq!(masses.len(), temperature.len());
    close(masses.iter().sum(), 30.0, 1e-10);
    let final_energy = masses
        .iter()
        .zip(&enthalpies)
        .map(|(m, h)| m * h)
        .sum::<f64>();
    let initial_energy = CP * (10.0 * (300.0 - INITIAL_K) + 20.0 * (350.0 - INITIAL_K));
    close(final_energy, initial_energy, 4.0 * ENERGY_TOLERANCE_J);
    for (&h, &t) in enthalpies.iter().zip(&temperature) {
        close(t, INITIAL_K + h / CP, 1e-11);
    }
    let final_row = transient
        .get("fine")
        .unwrap()
        .as_array()
        .unwrap()
        .last()
        .unwrap();
    for face in final_row
        .get("temperature_contact")
        .unwrap()
        .get("face_states")
        .unwrap()
        .as_array()
        .unwrap()
    {
        let cold = face.f64_field("side_a_temperature_k").unwrap();
        let hot = face.f64_field("side_b_temperature_k").unwrap();
        assert!(cold > 300.0 && hot < 350.0 && cold < hot);
    }
    assert!(interface_heat(&receipt) < 0.0);
    let peak = transient.get("sampled_peak").unwrap();
    close(peak.f64_field("value_k").unwrap(), 350.0, 1e-11);
    close(peak.f64_field("time_s").unwrap(), 0.0, 0.0);
}

#[test]
fn contact_trajectory_refuses_unimplemented_history_derivatives() {
    for (ordinal, enthalpy, output) in [
        (42, false, "temperature-final-adjoint"),
        (43, true, "temperature-sampled-peak-adjoint"),
    ] {
        let f = fixture(enthalpy, false);
        let (source, contact, ledger) = prepare(&f, ordinal, Some(output));
        let result = run(&f, &source, &contact, &ledger);
        assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
        assert!(
            format!("{}{}", result.stdout, result.stderr)
                .contains("cli-solve-conduction-contact-composition")
        );
    }
}

#[test]
fn contact_heat_changes_phase_fraction_on_a_latent_plateau_without_duplicate_storage() {
    let mut f = fixture(true, false);
    f.project.power.as_mut().unwrap()[0].watts.value = 0.0;
    let transient = time(&mut f);
    transient.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    for material in &mut transient.enthalpy.as_mut().unwrap().materials {
        if material.region == "hot" {
            material.initial_specific_enthalpy.value = 1000.0;
            material.knots = [
                (0.0, 300.0, 0.0),
                (500.0, 350.0, 0.0),
                (1500.0, 350.0, 1.0),
                (1931.5, 393.15, 1.0),
            ]
            .into_iter()
            .map(|(h, temperature, fraction)| TransientEnthalpyKnot {
                specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                liquid_mass_fraction: fraction,
            })
            .collect();
            material.source =
                "synthetic 1000 J/kg latent interval at 350 K; numerical phase-storage oracle"
                    .into();
        } else {
            material.initial_specific_enthalpy.value = CP * (300.0 - INITIAL_K);
        }
    }
    let (receipt, field) = solve(&f, 44);
    let transient = receipt.get("transient").unwrap();
    let endpoint = transient.get("endpoint").unwrap();
    let mass = values(endpoint, "reference_nodal_mass_kg");
    let enthalpy = values(endpoint, "specific_enthalpy_j_kg");
    let fraction = values(endpoint, "liquid_mass_fraction");
    let temperature = values(&field, "temperature");
    assert_eq!(mass.len(), enthalpy.len());
    assert_eq!(mass.len(), fraction.len());
    assert_eq!(mass.len(), temperature.len());
    close(mass.iter().sum(), 2.0 * DENSITY, 1e-10);
    let initial_energy = DENSITY * (CP * (300.0 - INITIAL_K) + 1000.0);
    let final_energy = mass.iter().zip(&enthalpy).map(|(m, h)| m * h).sum::<f64>();
    close(final_energy, initial_energy, 4.0 * ENERGY_TOLERANCE_J);
    let mut latent_nodes = 0;
    let mut cold_nodes = 0;
    for ((&h, &t), &liquid) in enthalpy.iter().zip(&temperature).zip(&fraction) {
        if liquid > 0.0 {
            latent_nodes += 1;
            close(t, 350.0, 1e-12);
            close(liquid, (h - 500.0) / 1000.0, 1e-13);
            assert!(
                liquid < 0.5,
                "contact heat must consume the original hot phase storage"
            );
        } else {
            cold_nodes += 1;
            assert!(t > 300.0 && t < 350.0);
            close(h, CP * (t - INITIAL_K), 1e-10);
        }
    }
    assert!(latent_nodes > 0 && cold_nodes > 0);
    let claim = temperature_pack().2;
    for grid in ["coarse", "fine"] {
        for row in transient.get(grid).unwrap().as_array().unwrap() {
            check_contact(row, &claim);
            close(row.f64_field("net_input_w").unwrap(), 0.0, 1e-12);
        }
    }
    let energy = transient.get("energy").unwrap();
    close(
        energy.f64_field("stored_change_j").unwrap(),
        0.0,
        4.0 * ENERGY_TOLERANCE_J,
    );
    assert!(interface_heat(&receipt) < 0.0);
}
