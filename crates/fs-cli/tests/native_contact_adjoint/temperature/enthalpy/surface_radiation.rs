//! The actual CLI composes source surfaces, latent storage and local radiation.
use super::*;

fn surface_pack(lower: f64, upper: f64) -> (Vec<u8>, String, String) {
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: [
            [lower, 0.5e6],
            [upper, 0.5e6],
            [lower, 1.5e6],
            [upper, 1.5e6],
        ]
        .into_iter()
        .map(|p| (p, 0.001 * p[0] + 1e-8 * (p[1] - 1e6)))
        .collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

fn audit(time: &J, initial_j: f64, card: &str, claim: &str) {
    let summary = time.get("temperature_contact").unwrap();
    let radiative = time.get("radiation").unwrap();
    let mut responses = 0.0;
    let mut updates = 0.0;
    let mut krylov = 0.0;
    for grid in ["coarse", "fine"] {
        let rows = time.get(grid).unwrap().as_array().unwrap();
        let mut grid_responses = 0.0;
        let mut net_j = 0.0;
        for row in rows {
            let contact = row.get("temperature_contact").unwrap();
            let radiation = row.get("radiation").unwrap();
            assert_eq!(
                contact.str_field("method"),
                Some("source-backed-contact-enthalpy-fixed-point")
            );
            assert_eq!(radiation.str_field("model"), Some("pointwise-p1"));
            assert!(
                n(contact, "physical_residual_norm_j")
                    <= n(contact, "physical_residual_tolerance_j")
            );
            assert!(
                n(radiation, "physical_residual_j")
                    <= n(radiation, "physical_residual_threshold_j")
            );
            assert!(
                n(contact, "physical_energy_residual_j").abs() <= n(time, "energy_tolerance_j")
            );
            close(
                n(contact, "thermal_responses"),
                n(radiation, "solid_solves"),
                0.0,
            );
            close(
                n(contact, "nonlinear_updates"),
                n(radiation, "solid_iterations"),
                0.0,
            );
            close(
                n(contact, "krylov_iterations"),
                n(radiation, "krylov_iterations"),
                0.0,
            );
            assert!(
                n(contact, "thermal_responses") <= n(summary, "max_thermal_responses_per_step")
            );
            assert!(n(contact, "krylov_iterations") >= n(row, "krylov_iterations"));
            grid_responses += n(contact, "thermal_responses");
            updates += n(contact, "nonlinear_updates");
            krylov += n(contact, "krylov_iterations");
            for state in contact.get("face_states").unwrap().as_array().unwrap() {
                close(
                    n(state, "resistance_m2_k_w"),
                    0.001 * n(state, "query_temperature_k"),
                    1e-13,
                );
                assert_eq!(state.str_field("interface_card"), Some(card));
                assert_eq!(state.str_field("selected_claim"), Some(claim));
            }
            let dt = n(row, "dt_s");
            let net = n(row, "source_w") + n(radiation, "physical_dirichlet_in_w")
                - n(radiation, "prescribed_neumann_out_w")
                - n(radiation, "convective_out_w")
                - n(radiation, "nonlinear_radiative_out_w");
            close(n(row, "net_input_w"), net, 2e-6);
            close(
                n(row, "stored_energy_change_j"),
                dt * net,
                n(time, "energy_tolerance_j"),
            );
            net_j += dt * net;
            if let Some(natural) = row.get("natural") {
                close(
                    n(natural, "solid_solves"),
                    n(contact, "thermal_responses"),
                    0.0,
                );
                close(
                    n(natural, "physical_convective_out_w"),
                    n(radiation, "convective_out_w"),
                    0.0,
                );
                close(
                    n(natural, "physical_radiative_out_w"),
                    n(radiation, "nonlinear_radiative_out_w"),
                    0.0,
                );
            }
        }
        assert!(grid_responses > rows.len() as f64);
        close(
            n(summary, &format!("{grid}_thermal_responses")),
            grid_responses,
            0.0,
        );
        responses += grid_responses;
        if grid == "fine" {
            close(
                stored(time),
                initial_j + net_j,
                rows.len() as f64 * n(time, "energy_tolerance_j"),
            );
        }
    }
    for (field, value) in [
        ("total_thermal_responses", responses),
        ("total_nonlinear_updates", updates),
        ("total_krylov_iterations", krylov),
    ] {
        close(n(summary, field), value, 0.0);
    }
    close(n(radiative, "total_trials"), responses, 0.0);
    close(n(radiative, "total_newton_updates"), updates, 0.0);
    close(n(radiative, "total_krylov_iterations"), krylov, 0.0);
    let nonlinear = time.get("nonlinear").unwrap();
    assert!(updates <= n(nonlinear, "total_newton_limit"));
    assert!(krylov <= n(nonlinear, "total_krylov_limit"));
    assert_eq!(
        nonlinear.str_field("history_scope"),
        Some("final-inner-newton-solve-of-each-accepted-contact-endpoint")
    );
}

#[test]
fn source_surface_enthalpy_runs_ambient_closed_and_open_pointwise_radiation() {
    for model in 1..=3 {
        let mut f = fixture(0.2, 0.1);
        radiation(&mut f, model);
        let (reference, _) = receipt(&f, 150 + model);
        let pack = surface_pack(200.0, 450.0);
        let (card, claim) = (pack.1.clone(), pack.2.clone());
        let (actual, _) = receipt_with_pack(&f, 155 + model, pack);
        let time = actual.get("transient").unwrap();
        audit(time, INITIAL_ENERGY_J, &card, &claim);
        let actual_h = values(time.get("endpoint").unwrap(), "specific_enthalpy_j_kg");
        let reference_h = values(
            reference.get("transient").unwrap().get("endpoint").unwrap(),
            "specific_enthalpy_j_kg",
        );
        for (a, b) in actual_h.iter().zip(reference_h) {
            close(*a, b, 3e-5);
        }
        let report = actual.get("radiation").unwrap();
        if model == 2 {
            close(stored(time), INITIAL_ENERGY_J + 1.4, 4e-6);
            let patches = report.get("surfaces").unwrap().as_array().unwrap();
            assert!(patches.iter().any(|p| n(p, "outward_heat_w") > 1.0));
            assert!(patches.iter().any(|p| n(p, "outward_heat_w") < -1.0));
        } else if model == 3 {
            close(
                n(report, "ambient_absorbed_w"),
                n(report, "radiative_out_w"),
                1e-6,
            );
        }
    }
}

#[test]
fn source_surface_pointwise_enthalpy_preserves_the_natural_convection_outer_owner() {
    let mut f = fixture(0.2, 0.1);
    radiation(&mut f, 3);
    let setup = f
        .project
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    setup
        .boundaries
        .iter_mut()
        .find(|row| row.target == "hot")
        .unwrap()
        .condition = ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(293.15, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    };
    setup
        .transient
        .as_mut()
        .unwrap()
        .enthalpy
        .as_mut()
        .unwrap()
        .materials
        .iter_mut()
        .find(|row| row.region == "cold")
        .unwrap()
        .initial_specific_enthalpy
        .value = 70.0;
    let (reference, _) = receipt(&f, 160);
    let pack = surface_pack(200.0, 450.0);
    let (card, claim) = (pack.1.clone(), pack.2.clone());
    let (actual, _) = receipt_with_pack(&f, 161, pack);
    let time = actual.get("transient").unwrap();
    audit(time, 70.0 + 2.0 * 218.5, &card, &claim);
    let original = values(
        reference.get("transient").unwrap().get("endpoint").unwrap(),
        "specific_enthalpy_j_kg",
    );
    for (a, b) in values(time.get("endpoint").unwrap(), "specific_enthalpy_j_kg")
        .iter()
        .zip(original)
    {
        close(*a, b, 3e-5);
    }
    assert!(n(time.get("natural").unwrap(), "total_iterations") > n(time, "total_steps"));
    assert!(
        values(time.get("endpoint").unwrap(), "liquid_mass_fraction")
            .iter()
            .any(|f| *f > 0.0 && *f < 1.0)
    );
}

#[test]
fn source_surface_pointwise_enthalpy_refuses_late_domain_exit_and_retries_original_history() {
    let mut f = fixture(1.0, 0.25);
    radiation(&mut f, 2);
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
    let (project, contact, ledger) = prepare_with_pack(&f, 162, None, surface_pack(293.0, 293.3));
    let result = run_with_radiation_pack(&f, &project, &contact, &ledger);
    assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
    let refusal = format!("{}{}", result.stdout, result.stderr);
    assert!(
        refusal.contains("cli-solve-conduction-transient"),
        "{refusal}"
    );
    assert!(!refusal.contains("run_receipt"), "{refusal}");
    let pack = surface_pack(200.0, 450.0);
    let (actual, field) = receipt_with_pack(&f, 163, pack);
    let time = actual.get("transient").unwrap();
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
    close(stored(time), 50.0, 8e-6);
    assert!(values(&field, "temperature").iter().any(|t| *t > 293.3));
}
