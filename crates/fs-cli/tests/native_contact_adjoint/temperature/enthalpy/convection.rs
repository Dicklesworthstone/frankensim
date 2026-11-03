//! Native convection and live contact share the original latent history.

use super::*;

const AMBIENT_K: f64 = 293.15;
const LENGTH_M: f64 = 0.1;
const NATURAL_CARD: &str = "convection.churchill-chu-vertical-plate";
const INITIAL_J: f64 = 70.0 + 2.0 * 218.5;

fn natural(f: &mut Fixture, ambient: f64) {
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.iter_mut().find(|row| row.target == "hot").unwrap().condition =
        ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: QtyAny::new(LENGTH_M, dims::LENGTH),
            ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
            correlation: NATURAL_CARD.into(),
        };
}

fn natural_fixture() -> Fixture {
    let mut f = fixture(0.5, 0.25);
    natural(&mut f, AMBIENT_K);
    let policy = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap();
    policy.enthalpy.as_mut().unwrap().materials.iter_mut()
        .find(|row| row.region == "cold").unwrap().initial_specific_enthalpy.value = 70.0;
    f
}

fn surface_pack() -> (Vec<u8>, String, String) {
    // The same resistance as the curve at the original fixed pressure 1 MPa.
    temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: [[200.0, 0.5e6], [450.0, 0.5e6], [200.0, 1.5e6], [450.0, 1.5e6]]
            .into_iter().map(|p| (p, 0.001 * p[0] + 1e-8 * (p[1] - 1e6))).collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    })
}

fn solve_with_pack(f: &Fixture, ordinal: usize, pack: (Vec<u8>, String, String)) -> (J, J) {
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, None, pack);
    let output = run(f, &source, &contact, &ledger);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", output.stdout, output.stderr);
    let output = J::parse(&output.stdout).unwrap();
    let completed = artifact(&ledger, output.str_field("run_receipt").unwrap());
    let stage = completed.get("stages").unwrap().as_array().unwrap().iter()
        .find(|stage| stage.str_field("stage") == Some("conduction")).unwrap();
    let conduction = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, conduction.str_field("solution_artifact").unwrap());
    (conduction, field)
}

// Independently evaluate the declared Churchill-Chu card and native air policy.
fn coefficient(wall: f64, ambient: f64) -> (f64, f64, f64) {
    let film = 0.5 * (wall + ambient);
    let density = 101325.0 / (287.05 * film);
    let kinematic = 1.846e-5 / density;
    let diffusivity = kinematic / 0.707;
    let rayleigh = 9.80665 * (wall - ambient).abs() * LENGTH_M.powi(3)
        / (film * kinematic * diffusivity);
    let factor = 0.825 + 0.387 * rayleigh.powf(1.0 / 6.0)
        / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0);
    let nusselt = factor * factor;
    (0.0263 * nusselt / LENGTH_M, rayleigh, nusselt)
}

fn audit(time: &J, initial_j: f64, input_j: f64, card: &str, claim: &str, surface: bool) -> f64 {
    let method = if surface {
        "source-backed-contact-enthalpy-fixed-point"
    } else { "source-backed-contact-enthalpy-newton" };
    let contact_work = time.get("temperature_contact").unwrap();
    let natural_work = time.get("natural").unwrap();
    assert_eq!(contact_work.str_field("method"), Some(method));
    assert_eq!(natural_work.str_field("method"), Some("implicit-area-mean-natural-convection"));
    let mut responses = 0.0;
    let mut updates = 0.0;
    let mut krylov = 0.0;
    let mut natural_iterations = 0.0;
    let mut fine_out_j = 0.0;
    let mut minimum_query = f64::INFINITY;
    let mut maximum_query = f64::NEG_INFINITY;
    for grid in ["coarse", "fine"] {
        let rows = time.get(grid).unwrap().as_array().unwrap();
        let mut delivered_j = 0.0;
        let mut removed_j = 0.0;
        for row in rows {
            let contact = row.get("temperature_contact").unwrap();
            let natural = row.get("natural").unwrap();
            assert_eq!(contact.str_field("method"), Some(method));
            for physical in [contact, natural] {
                assert!(n(physical, "physical_residual_norm_j") <= n(physical, "physical_residual_tolerance_j"));
                assert!(n(physical, "physical_energy_residual_j").abs() <= n(time, "energy_tolerance_j"));
            }
            close(n(contact, "physical_residual_norm_j"), n(natural, "physical_residual_norm_j"), 0.0);
            close(n(contact, "physical_energy_residual_j"), n(natural, "physical_energy_residual_j"), 0.0);
            close(n(row, "residual_j"), n(natural, "physical_residual_norm_j"), 0.0);
            close(n(contact, "thermal_responses"), n(natural, "solid_solves"), 0.0);
            close(n(contact, "nonlinear_updates"), n(natural, "nonlinear_iterations"), 0.0);
            close(n(contact, "krylov_iterations"), n(natural, "krylov_iterations"), 0.0);
            assert!(n(contact, "thermal_responses") <= n(contact_work, "max_thermal_responses_per_step"));
            assert!(n(contact, "krylov_iterations") >= n(row, "krylov_iterations"));
            responses += n(contact, "thermal_responses");
            updates += n(contact, "nonlinear_updates");
            krylov += n(contact, "krylov_iterations");
            natural_iterations += n(natural, "iterations");
            let laws = natural.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap();
            assert_eq!(laws.len(), 1);
            let law = &laws[0];
            assert_eq!(law.str_field("target"), Some("hot"));
            assert_eq!(law.str_field("card"), Some(NATURAL_CARD));
            assert_eq!(law.get("in_domain"), Some(&J::Bool(true)));
            close(n(law, "characteristic_length_m"), LENGTH_M, 0.0);
            let wall = n(law, "mean_wall_k");
            let ambient = n(law, "ambient_k");
            let (h, rayleigh, nusselt) = coefficient(wall, ambient);
            assert!((0.1..=1e12).contains(&rayleigh));
            close(n(law, "rayleigh"), rayleigh, 1e-8 * rayleigh.max(1.0));
            close(n(law, "nusselt"), nusselt, 1e-8);
            close(n(law, "htc_w_m2_k"), h, 1e-8);
            let outward = n(natural, "physical_convective_out_w");
            // The hot unit cube has five exterior faces; the joint is internal.
            close(outward, 5.0 * h * (wall - ambient), 2e-6);
            close(outward, n(law, "heat_rate_w"), 2e-6);
            assert!(outward > 0.0);
            close(n(natural, "physical_radiative_out_w"), 0.0, 0.0);
            let dt = n(row, "dt_s");
            let net_w = n(row, "source_w") - outward;
            close(n(row, "net_input_w"), net_w, 2e-6);
            close(n(row, "stored_energy_change_j"), dt * net_w, n(time, "energy_tolerance_j"));
            delivered_j += dt * n(row, "source_w");
            removed_j += dt * outward;
            for state in contact.get("face_states").unwrap().as_array().unwrap() {
                let query = n(state, "query_temperature_k");
                close(query, 0.5 * n(state, "side_a_temperature_k") + 0.5 * n(state, "side_b_temperature_k"), 1e-12);
                close(n(state, "resistance_m2_k_w"), 0.001 * query, 1e-13);
                assert_eq!(state.str_field("interface_card"), Some(card));
                assert_eq!(state.str_field("selected_claim"), Some(claim));
                minimum_query = minimum_query.min(query);
                maximum_query = maximum_query.max(query);
            }
        }
        close(delivered_j, input_j, 1e-10);
        if grid == "fine" {
            fine_out_j = removed_j;
            close(stored(time), initial_j + input_j - removed_j,
                rows.len() as f64 * n(time, "energy_tolerance_j"));
        }
    }
    close(n(contact_work, "total_thermal_responses"), responses, 0.0);
    close(n(contact_work, "total_nonlinear_updates"), updates, 0.0);
    close(n(contact_work, "total_krylov_iterations"), krylov, 0.0);
    close(n(natural_work, "total_solid_solves"), responses, 0.0);
    close(n(natural_work, "total_newton_updates"), updates, 0.0);
    close(n(natural_work, "total_krylov_iterations"), krylov, 0.0);
    close(n(natural_work, "total_iterations"), natural_iterations, 0.0);
    assert!(natural_iterations > n(time, "total_steps"));
    if surface { assert!(responses > natural_iterations); }
    else { close(responses, natural_iterations, 0.0); }
    assert!(maximum_query - minimum_query > 1e-5);
    close(n(time.get("energy").unwrap(), "stored_change_j"), stored(time) - initial_j, 2e-6);
    fine_out_j
}

#[test]
fn native_enthalpy_natural_contact_keeps_live_source_and_latent_energy() {
    let mut reference_h: Option<Vec<f64>> = None;
    for surface in [false, true] {
        let f = natural_fixture();
        let source = if surface { surface_pack() } else { temperature_pack() };
        let (card, claim) = (source.1.clone(), source.2.clone());
        let (conduction, field) = solve_with_pack(&f, if surface { 121 } else { 120 }, source);
        let time = conduction.get("transient").unwrap();
        let removed_j = audit(time, INITIAL_J, 3.5, &card, &claim, surface);
        assert!(removed_j > 3.5, "natural cooling must remove real stored heat");
        let endpoint = time.get("endpoint").unwrap();
        close(values(endpoint, "reference_nodal_mass_kg").iter().sum(), 3.0, 1e-12);
        let h = values(endpoint, "specific_enthalpy_j_kg");
        let fraction = values(endpoint, "liquid_mass_fraction");
        let temperature = values(&field, "temperature");
        let mut advanced_latent = false;
        for ((&h, &fraction), &temperature) in h.iter().zip(&fraction).zip(&temperature) {
            if fraction > 1e-5 && fraction < 1.0 - 1e-5 {
                close(temperature, MELT_K, 1e-12);
                close(h, 20.0 + 100.0 * fraction, 1e-10);
                advanced_latent |= h > 70.0 + 1e-5;
            }
        }
        assert!(advanced_latent, "contact heat must change the original cold latent fraction");
        let fluxes = conduction.get("interfaces").unwrap().get("fluxes").unwrap().as_array().unwrap();
        assert!(fluxes.iter().any(|row| n(row, "heat_rate_a_to_b_w") < 0.0));
        if let Some(reference) = &reference_h {
            for (&actual, &expected) in h.iter().zip(reference) {
                close(actual, expected, 3e-5);
            }
        } else {
            reference_h = Some(h.clone());
            let (frozen, _, _) = f.solve(0.300075, None, 122);
            let frozen_h = values(frozen.get("transient").unwrap().get("endpoint").unwrap(), "specific_enthalpy_j_kg");
            assert!(h.iter().zip(frozen_h).any(|(live, frozen)| (live - frozen).abs() > 1e-5),
                "fresh source evidence must accompany a physically evolving resistance");
        }
    }
}

#[test]
fn native_enthalpy_natural_contact_refuses_unsupported_equilibrium() {
    let mut f = natural_fixture();
    natural(&mut f, MELT_K);
    let policy = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap();
    for material in &mut policy.enthalpy.as_mut().unwrap().materials {
        material.initial_specific_enthalpy.value = 70.0;
    }
    for step in &mut policy.power_schedules[0].steps { step.watts.value = 0.0; }
    let (source, contact, ledger) = prepare(&f, 130, None);
    let output = run(&f, &source, &contact, &ledger);
    assert_ne!(output.exit_code, fs_cli::exit::SUCCESS);
    let refusal = format!("{}{}", output.stdout, output.stderr);
    assert!(refusal.contains("cli-solve-conduction-natural"), "{refusal}");
    assert!(!refusal.contains("run_receipt"), "unsupported natural equilibrium cannot publish a trajectory: {refusal}");
}

#[test]
fn native_enthalpy_natural_contact_refuses_late_source_exit_and_retries_original_history() {
    let mut f = fixture(1.0, 0.25);
    natural(&mut f, 293.0);
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 293.0;
    envelope.ambient_hi.value = 293.25;
    let requirement = &mut f.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 293.25;
    requirement.margin.value = 0.0;
    let policy = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap();
    for material in &mut policy.enthalpy.as_mut().unwrap().materials {
        material.initial_specific_enthalpy.value = 0.0;
    }
    policy.power_schedules[0].steps = [(0.5, 0.0), (1.0, 100.0)].into_iter()
        .map(|(until, watts)| TransientPowerStep {
            until: QtyAny::new(until, dims::TIME), watts: QtyAny::new(watts, dims::POWER),
        }).collect();
    let limited = temperature_pack_from_knots(vec![(293.0, 0.293), (293.3, 0.2933)]);
    let (source, contact, ledger) = prepare_with_pack(&f, 140, None, limited);
    let output = run(&f, &source, &contact, &ledger);
    assert_ne!(output.exit_code, fs_cli::exit::SUCCESS);
    let refusal = format!("{}{}", output.stdout, output.stderr);
    assert!(refusal.contains("cli-solve-conduction-transient"), "{refusal}");
    assert!(!refusal.contains("run_receipt"), "an unsupported contact endpoint cannot publish history: {refusal}");
    let wide = temperature_pack();
    let (card, claim) = (wide.1.clone(), wide.2.clone());
    let (conduction, field) = solve_with_pack(&f, 141, wide);
    let time = conduction.get("transient").unwrap();
    audit(time, 0.0, 50.0, &card, &claim, false);
    for row in time.get("coarse").unwrap().as_array().unwrap().iter()
        .filter(|row| n(row, "time_s") <= 0.5)
    {
        close(n(row, "source_w"), 0.0, 0.0);
        assert!(n(row, "stored_energy_change_j") < 0.0);
    }
    assert!(values(&field, "temperature").iter().any(|&temperature| temperature > 293.3));
}
