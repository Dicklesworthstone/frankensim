//! Imported nonlinear contacts must transport real latent enthalpy, not a
//! temperature-derived replacement history or a frozen interface resistance.

use super::*;
#[path = "enthalpy/convection.rs"]
mod convection;
#[path = "enthalpy/airflow.rs"]
mod airflow;
#[path = "enthalpy/surface_radiation.rs"]
mod surface_radiation;
use fs_project::{
    ConductionEnclosure, ConductionTransient, EnclosureRadiatingSurface,
    EnclosureViewFactorEvidence, RadiationSpatialModel, ThermalBoundaryCondition,
    TransientEnthalpy, TransientEnthalpyKnot, TransientInitialTemperature, TransientPowerStep,
    TransientRegionEnthalpy, TransientRegionPower, spec::dims,
};
use fs_qty::QtyAny;

const MELT_K: f64 = 295.15;
// The retained regions are unit cubes: 1 kg at h=5 and 2 kg at h=218.5.
const INITIAL_ENERGY_J: f64 = 442.0;

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
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
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn fixture(horizon: f64, max_step: f64) -> Fixture {
    let mut f = Fixture::new(true, false);
    let identity =
        fs_blake3::ContentHash::from_hex(&f.project.materials.as_ref().unwrap()[0].card).unwrap();
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
        initial_temperature: TransientInitialTemperature::FromRegionalEnthalpy,
        horizon: QtyAny::new(horizon, dims::TIME),
        max_step: QtyAny::new(max_step, dims::TIME),
        max_steps: 256,
        energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
        capacities: Vec::new(),
        boundary_schedules: Vec::new(),
        power_schedules: vec![TransientRegionPower {
            region: "hot".into(),
            source: "synthetic delivered-power pulse".into(),
            steps: [(0.375 * horizon, 12.0), (horizon, 4.0)]
                .into_iter()
                .map(|(until, watts)| TransientPowerStep {
                    until: QtyAny::new(until, dims::TIME),
                    watts: QtyAny::new(watts, dims::POWER),
                })
                .collect(),
        }],
        enthalpy: Some(TransientEnthalpy {
            max_newton_iterations: 64,
            max_backtracks: 32,
            materials: [("cold", 1.0, 5.0), ("hot", 2.0, 218.5)]
                .into_iter()
                .map(|(region, density, initial)| TransientRegionEnthalpy {
                    region: region.into(),
                    material_card_identity: identity,
                    source: "synthetic Cp=10 and latent heat=100; numerical verification only"
                        .into(),
                    reference_density: QtyAny::new(density, dims::DENSITY),
                    initial_specific_enthalpy: QtyAny::new(initial, dims::SPECIFIC_ENTHALPY),
                    phase_conductivity: None,
                    knots: [
                        (-200.0, 273.15, 0.0),
                        (20.0, MELT_K, 0.0),
                        (120.0, MELT_K, 1.0),
                        (1000.0, 383.15, 1.0),
                    ]
                    .into_iter()
                    .map(|(h, temperature, fraction)| TransientEnthalpyKnot {
                        specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                        temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                        liquid_mass_fraction: fraction,
                    })
                    .collect(),
                })
                .collect(),
        }),
    });
    f
}

fn receipt(f: &Fixture, ordinal: usize) -> (J, J) {
    receipt_with_pack(f, ordinal, temperature_pack())
}

fn receipt_with_pack(f: &Fixture, ordinal: usize, pack: (Vec<u8>, String, String)) -> (J, J) {
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, None, pack);
    let output = run_with_radiation_pack(f, &source, &contact, &ledger);
    assert_eq!(
        output.exit_code,
        fs_cli::exit::SUCCESS,
        "{}\n{}",
        output.stdout,
        output.stderr
    );
    let output = J::parse(&output.stdout).unwrap();
    let run = artifact(&ledger, output.str_field("run_receipt").unwrap());
    let stage = run
        .get("stages")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .find(|stage| stage.str_field("stage") == Some("conduction"))
        .unwrap();
    let conduction = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, conduction.str_field("solution_artifact").unwrap());
    (conduction, field)
}

fn run_with_radiation_pack(f: &Fixture, source: &Path, contact: &Path, ledger: &Path) -> fs_cli::CommandOutput {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    fs_cli::run(vec![
        "--json".into(),
        "solve".into(),
        source.to_str().unwrap().into(),
        ledger.to_str().unwrap().into(),
        "--materials".into(),
        f.dir.join("solid.fsmcdpk").to_str().unwrap().into(),
        "--materials".into(),
        root.join("data/reference-project/gray-surface.fsmcdpk")
            .to_str()
            .unwrap()
            .into(),
        "--interfaces".into(),
        contact.to_str().unwrap().into(),
    ])
}

fn audit_contact(time: &J) {
    let (_, card, selected) = temperature_pack();
    let mut responses = 0.0;
    let mut krylov = 0.0;
    let mut min_query = f64::INFINITY;
    let mut max_query = f64::NEG_INFINITY;
    for grid in ["coarse", "fine"] {
        let rows = time.get(grid).unwrap().as_array().unwrap();
        assert!(
            rows.iter()
                .any(|row| n(row, "time_s") == 0.375 * n(time, "final_time_s"))
        );
        for row in rows {
            let contact = row.get("temperature_contact").unwrap();
            assert_eq!(
                contact.str_field("method"),
                Some("source-backed-contact-enthalpy-newton")
            );
            assert!(
                n(contact, "physical_residual_norm_j")
                    <= n(contact, "physical_residual_tolerance_j")
            );
            assert!(
                n(contact, "physical_energy_residual_j").abs() <= n(time, "energy_tolerance_j")
            );
            close(
                n(row, "residual_j"),
                n(contact, "physical_residual_norm_j"),
                0.0,
            );
            close(
                n(row, "krylov_iterations"),
                n(contact, "krylov_iterations"),
                0.0,
            );
            responses += n(contact, "thermal_responses");
            krylov += n(contact, "krylov_iterations");
            for state in contact.get("face_states").unwrap().as_array().unwrap() {
                let query = n(state, "query_temperature_k");
                close(
                    query,
                    0.5 * (n(state, "side_a_temperature_k") + n(state, "side_b_temperature_k")),
                    1e-12,
                );
                close(n(state, "resistance_m2_k_w"), 0.001 * query, 1e-14);
                assert_eq!(state.str_field("interface_card"), Some(card.as_str()));
                assert_eq!(state.str_field("selected_claim"), Some(selected.as_str()));
                min_query = min_query.min(query);
                max_query = max_query.max(query);
            }
        }
    }
    let contact = time.get("temperature_contact").unwrap();
    close(n(contact, "total_thermal_responses"), responses, 0.0);
    close(
        n(contact, "total_thermal_responses"),
        n(time, "total_steps"),
        0.0,
    );
    close(n(contact, "total_krylov_iterations"), krylov, 0.0);
    assert!(
        max_query - min_query > 1e-5,
        "the original source must follow the evolving faces"
    );
    close(n(time.get("nonlinear").unwrap(), "residual_rtol"), 0.0, 0.0);
}

fn stored(time: &J) -> f64 {
    let endpoint = time.get("endpoint").unwrap();
    values(endpoint, "specific_enthalpy_j_kg")
        .iter()
        .zip(values(endpoint, "reference_nodal_mass_kg"))
        .map(|(h, mass)| h * mass)
        .sum()
}

#[test]
fn native_enthalpy_contact_crosses_latent_plateaus_with_original_history_and_live_source() {
    let f = fixture(2.0, 0.25);
    let (conduction, field) = receipt(&f, 90);
    let time = conduction.get("transient").unwrap();
    audit_contact(time);
    let endpoint = time.get("endpoint").unwrap();
    close(
        values(endpoint, "reference_nodal_mass_kg").iter().sum(),
        3.0,
        1e-12,
    );
    // 12 W for 0.75 s + 4 W for 1.25 s: contact is entirely internal.
    close(stored(time), INITIAL_ENERGY_J + 14.0, 2e-6);
    close(
        n(time.get("energy").unwrap(), "stored_change_j"),
        14.0,
        2e-6,
    );
    let h = values(endpoint, "specific_enthalpy_j_kg");
    let fraction = values(endpoint, "liquid_mass_fraction");
    let temperature = values(&field, "temperature");
    let mut latent = 0;
    for ((&h, &fraction), &temperature) in h.iter().zip(&fraction).zip(&temperature) {
        if fraction > 1e-5 && fraction < 1.0 - 1e-5 {
            latent += 1;
            close(temperature, MELT_K, 1e-12);
            close(h, 20.0 + 100.0 * fraction, 1e-10);
        }
    }
    assert!(
        latent > 0,
        "heat must evolve a resolved latent state, not just sensible temperature"
    );
    let (frozen, _, _) = f.solve(0.299325, None, 91);
    let frozen_h = values(
        frozen.get("transient").unwrap().get("endpoint").unwrap(),
        "specific_enthalpy_j_kg",
    );
    assert!(
        h.iter().zip(&frozen_h).any(|(a, b)| (a - b).abs() > 1e-4),
        "a live resistance must change the heat partition from its initial frozen value"
    );
    close(
        stored(frozen.get("transient").unwrap()),
        INITIAL_ENERGY_J + 14.0,
        2e-6,
    );
}

fn radiation(f: &mut Fixture, model: usize) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let reference = fs_project::parse_sexpr_migrating(
        &std::fs::read_to_string(root.join("data/reference-project/cooling-radiation.fsim"))
            .unwrap(),
    )
    .unwrap()
    .decoded
    .spec;
    let mut ambient = reference
        .cooling
        .unwrap()
        .conduction
        .unwrap()
        .radiation
        .unwrap();
    ambient.spatial_model = Some(RadiationSpatialModel::Pointwise);
    ambient.max_iterations = 1;
    ambient.surfaces[0].name = "cold-radiator".into();
    ambient.surfaces[0].target = "cold".into();
    ambient.surfaces[0].reservoir_temperature.value = 300.0;
    let setup = f
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    if model == 1 {
        setup.radiation = Some(ambient);
    } else {
        let view = if model == 2 { 1.0 } else { 0.35 };
        setup.enclosure_radiation = Some(ConductionEnclosure {
            spatial_model: Some(RadiationSpatialModel::Pointwise),
            ambient_temperature: (model == 3).then_some(QtyAny::new(300.0, dims::TEMPERATURE)),
            surfaces: ["cold", "hot"]
                .into_iter()
                .map(|target| EnclosureRadiatingSurface {
                    name: format!("{target}-radiator"),
                    target: target.into(),
                    card: ambient.surfaces[0].card.clone(),
                    claim: ambient.surfaces[0].claim.clone(),
                    query_temperature: QtyAny::new(300.0, dims::TEMPERATURE),
                })
                .collect(),
            view_factors: vec![vec![0.0, view], vec![view, 0.0]],
            evidence: EnclosureViewFactorEvidence::Analytic {
                geometry: "synthetic equal-area exchange; not geometric visibility evidence".into(),
            },
            closure_tolerance: 1e-12,
            reciprocity_relative_tolerance: 1e-12,
            max_iterations: 1,
            temperature_tolerance: QtyAny::new(1e-9, dims::TEMPERATURE),
            heat_tolerance: QtyAny::new(1e-6, dims::POWER),
            relaxation: 0.5,
        });
    }
}

#[test]
fn native_enthalpy_contact_composes_pointwise_ambient_closed_and_open_radiation() {
    let baseline = fixture(0.2, 0.1);
    let (base, _) = receipt(&baseline, 100);
    let base_h = values(
        base.get("transient").unwrap().get("endpoint").unwrap(),
        "specific_enthalpy_j_kg",
    );
    for model in 1..=3 {
        let mut f = fixture(0.2, 0.1);
        radiation(&mut f, model);
        let (conduction, _) = receipt(&f, 100 + model);
        let time = conduction.get("transient").unwrap();
        audit_contact(time);
        let mut radiated_j = 0.0;
        for row in time.get("fine").unwrap().as_array().unwrap() {
            let radiation = row.get("radiation").unwrap();
            let outgoing_w = n(radiation, "nonlinear_radiative_out_w");
            close(n(radiation, "convective_out_w"), 0.0, 0.0);
            close(
                n(row, "stored_energy_change_j"),
                n(row, "dt_s") * (n(row, "source_w") - outgoing_w),
                1e-6,
            );
            radiated_j += n(row, "dt_s") * outgoing_w;
        }
        // The same exact 1.4 J pulse heats every case. Contact and reflected
        // internal exchange never become a second external heat account.
        close(stored(time), INITIAL_ENERGY_J + 1.4 - radiated_j, 2e-6);
        let h = values(time.get("endpoint").unwrap(), "specific_enthalpy_j_kg");
        assert!(h.iter().zip(&base_h).any(|(a, b)| (a - b).abs() > 1e-4));
        let final_radiation = conduction.get("radiation").unwrap();
        if model == 2 {
            close(radiated_j, 0.0, 1e-6);
            let surfaces = final_radiation.get("surfaces").unwrap().as_array().unwrap();
            assert!(surfaces.iter().any(|row| n(row, "outward_heat_w") > 1.0));
            assert!(surfaces.iter().any(|row| n(row, "outward_heat_w") < -1.0));
        } else {
            assert!(
                radiated_j.abs() > 1e-3,
                "an external radiation reservoir must exchange real stored heat"
            );
            if model == 3 {
                close(
                    n(final_radiation, "ambient_absorbed_w"),
                    n(final_radiation, "radiative_out_w"),
                    1e-6,
                );
            }
        }
    }
}

#[test]
fn native_enthalpy_contact_refuses_late_source_exit_without_completed_history() {
    let mut f = fixture(1.0, 0.25);
    f.project.envelope.as_mut().unwrap().ambient_hi.value = 293.25;
    let requirement = &mut f.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 293.25;
    requirement.margin.value = 0.0;
    let time = f
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap();
    for material in &mut time.enthalpy.as_mut().unwrap().materials {
        material.initial_specific_enthalpy.value = 0.0;
    }
    time.power_schedules[0].steps = [(0.5, 0.0), (1.0, 100.0)]
        .into_iter()
        .map(|(until, watts)| TransientPowerStep {
            until: QtyAny::new(until, dims::TIME),
            watts: QtyAny::new(watts, dims::POWER),
        })
        .collect();
    let source = temperature_pack_from_knots(vec![(293.0, 0.293), (293.3, 0.2933)]);
    let (project, contact, ledger) = prepare_with_pack(&f, 110, None, source);
    let result = run(&f, &project, &contact, &ledger);
    assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
    let refusal = format!("{}{}", result.stdout, result.stderr);
    assert!(
        refusal.contains("cli-solve-conduction-transient"),
        "{refusal}"
    );
    assert!(
        !refusal.contains("run_receipt"),
        "an unsupported endpoint cannot publish a completed trajectory: {refusal}"
    );
    // Two original zero-load intervals are exactly stationary and supported.
    // The same history completes when only its original source span is widened.
    let (conduction, field) = receipt(&f, 111);
    let time = conduction.get("transient").unwrap();
    for row in time
        .get("coarse")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .filter(|row| n(row, "time_s") <= 0.5)
    {
        close(n(row, "final_region_max_k"), 293.15, 1e-12);
        close(n(row, "stored_energy_change_j"), 0.0, 1e-12);
    }
    assert!(values(&field, "temperature").iter().any(|&t| t > 293.3));
}
