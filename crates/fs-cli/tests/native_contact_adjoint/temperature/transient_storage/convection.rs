//! G3: native capacity/contact/natural cooling agrees with equivalent storage.

use super::*;
use fs_project::{TransientPowerStep, TransientRegionPower};

const START_K: f64 = 300.0;
const LENGTH_M: f64 = 0.1;
const CARD: &str = "convection.churchill-chu-vertical-plate";

fn n(value: &J, key: &str) -> f64 { value.f64_field(key).unwrap() }

fn natural(f: &mut Fixture, ambient: f64) {
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.iter_mut().find(|row| row.target == "cold").unwrap().condition =
        ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: QtyAny::new(LENGTH_M, dims::LENGTH),
            ambient_temperature: QtyAny::new(ambient, dims::TEMPERATURE),
            correlation: CARD.into(),
        };
}

fn natural_fixture(enthalpy: bool) -> Fixture {
    let mut f = fixture(enthalpy, true);
    natural(&mut f, INITIAL_K);
    let policy = time(&mut f);
    policy.initial_temperature = TransientInitialTemperature::Uniform(
        QtyAny::new(START_K, dims::TEMPERATURE));
    policy.boundary_schedules[0].target = "hot".into();
    f
}

fn audit(row: &J, card: &str, claim: &str, ambient: f64) -> (f64, f64) {
    let query = check_contact(row, claim);
    let contact = row.get("temperature_contact").unwrap();
    let natural = row.get("natural").unwrap();
    for state in contact.get("face_states").unwrap().as_array().unwrap() {
        assert_eq!(state.str_field("interface_card"), Some(card));
    }
    for key in ["physical_residual_norm_j", "physical_energy_residual_j"] {
        close(n(contact, key), n(natural, key), 0.0);
    }
    assert!(n(natural, "physical_residual_norm_j") <= n(natural, "physical_residual_tolerance_j"));
    close(n(contact, "thermal_responses"), n(natural, "solid_solves"), 0.0);
    close(n(contact, "nonlinear_updates"), n(natural, "nonlinear_iterations"), 0.0);
    close(n(contact, "nonlinear_backtracks"), n(natural, "nonlinear_backtracks"), 0.0);
    close(n(contact, "krylov_iterations"), n(natural, "krylov_iterations"), 0.0);
    let laws = natural.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap();
    assert_eq!(laws.len(), 1);
    let law = &laws[0];
    assert_eq!(law.str_field("target"), Some("cold"));
    assert_eq!(law.str_field("card"), Some(CARD));
    assert_eq!(law.get("in_domain"), Some(&J::Bool(true)));
    close(n(law, "ambient_k"), ambient, 0.0);
    close(n(law, "characteristic_length_m"), LENGTH_M, 0.0);
    let wall = n(law, "mean_wall_k");
    // Independent Churchill-Chu evaluation using the declared native air policy.
    let film = 0.5 * (wall + ambient);
    let kinematic = 1.846e-5 / (101325.0 / (287.05 * film));
    let rayleigh = 9.80665 * (wall - ambient).abs() * LENGTH_M.powi(3)
        / (film * kinematic * (kinematic / 0.707));
    let factor = 0.825 + 0.387 * rayleigh.powf(1.0 / 6.0)
        / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0);
    let h = 0.0263 * factor * factor / LENGTH_M;
    assert!((0.1..=1e12).contains(&rayleigh));
    close(n(law, "rayleigh"), rayleigh, 1e-8 * rayleigh.max(1.0));
    close(n(law, "htc_w_m2_k"), h, 1e-8);
    close(n(natural.get("physical_coefficients_w_m2_k").unwrap(), "cold"), h, 1e-8);
    let outward = n(natural, "physical_convective_out_w");
    // The cold unit cube exposes five square metres; its contact face is internal.
    close(outward, 5.0 * h * (wall - ambient), 2e-6);
    close(outward, n(law, "heat_rate_w"), 2e-6);
    assert!(outward > 0.0);
    (query, outward)
}

#[test]
fn native_capacity_natural_contact_preserves_live_laws_and_both_grid_storage() {
    // The existing sensible-enthalpy fixture has exactly the same Cv = rho*Cp.
    // It provides an independent storage owner and the retained mesh mass weights.
    let (reference, reference_field) = solve(&natural_fixture(true), 160);
    let reference_time = reference.get("transient").unwrap();
    let masses = values(reference_time.get("endpoint").unwrap(), "reference_nodal_mass_kg");
    close(masses.iter().sum(), 20.0, 1e-10);
    for surface in [false, true] {
        let source = if surface { temperature_surface_pack() } else { temperature_pack() };
        let (card, claim) = (source.1.clone(), source.2.clone());
        let (receipt, field) = solve_with_pack(&natural_fixture(false), 161 + usize::from(surface), source);
        let trajectory = receipt.get("transient").unwrap();
        let contact_work = trajectory.get("temperature_contact").unwrap();
        let natural_work = trajectory.get("natural").unwrap();
        let method = if surface { "source-backed-contact-fixed-point" } else { "source-backed-contact-newton" };
        assert_eq!(contact_work.str_field("method"), Some(method));
        let (mut responses, mut iterations, mut fine_energy) = (0.0, 0.0, 0.0);
        let (mut min_query, mut max_query) = (f64::INFINITY, f64::NEG_INFINITY);
        for (grid, ends) in [
            ("coarse", vec![0.25, 0.75, 1.0]),
            ("fine", vec![0.125, 0.25, 0.5, 0.75, 0.875, 1.0]),
        ] {
            let rows = trajectory.get(grid).unwrap().as_array().unwrap();
            let reference_rows = reference_time.get(grid).unwrap().as_array().unwrap();
            assert_eq!(rows.len(), ends.len());
            assert_eq!(rows.len(), reference_rows.len());
            let (mut start, mut integrated, mut grid_responses) = (0.0, 0.0, 0.0);
            for ((row, reference_row), end) in rows.iter().zip(reference_rows).zip(ends) {
                close(n(row, "time_s"), end, 0.0);
                close(n(row, "dt_s"), end - start, 0.0);
                let (query, removed_w) = audit(row, &card, &claim, INITIAL_K);
                min_query = min_query.min(query);
                max_query = max_query.max(query);
                let flux = if start < 0.25 { -2.0 } else if start < 0.75 { 1.0 } else { 0.0 };
                let net_w = 10.0 - 5.0 * flux - removed_w;
                close(n(row, "source_w"), 10.0, 1e-10);
                close(n(row, "net_input_w"), net_w, 2e-6);
                close(n(row, "stored_energy_change_j"), (end - start) * net_w, 2.0 * ENERGY_TOLERANCE_J);
                close(n(row, "stored_energy_change_j"), n(reference_row, "stored_energy_change_j"), 4e-6);
                close(n(row, "final_region_max_k"), n(reference_row, "final_region_max_k"), 2e-6);
                let contact = row.get("temperature_contact").unwrap();
                assert_eq!(contact.str_field("method"), Some(method));
                assert!(n(contact, "thermal_responses") <= n(contact_work, "max_thermal_responses_per_step"));
                grid_responses += n(contact, "thermal_responses");
                iterations += n(row.get("natural").unwrap(), "iterations");
                integrated += (end - start) * net_w;
                start = end;
            }
            close(n(contact_work, &format!("{grid}_thermal_responses")), grid_responses, 0.0);
            close(n(natural_work, &format!("{grid}_solid_solves")), grid_responses, 0.0);
            responses += grid_responses;
            if grid == "fine" { fine_energy = integrated; }
        }
        close(n(contact_work, "total_thermal_responses"), responses, 0.0);
        close(n(natural_work, "total_solid_solves"), responses, 0.0);
        close(n(natural_work, "total_iterations"), iterations, 0.0);
        assert!(iterations > n(trajectory, "total_steps"));
        if surface { assert!(responses > iterations); } else { close(responses, iterations, 0.0); }
        assert!(max_query - min_query > 1e-4, "contact must requery changing temperatures");
        let temperature = values(&field, "temperature");
        assert_eq!(temperature.len(), masses.len());
        let stored: f64 = temperature.iter().zip(&masses).map(|(t, mass)| mass * CP * (t - START_K)).sum();
        close(stored, fine_energy, 6.0 * ENERGY_TOLERANCE_J);
        let energy = trajectory.get("energy").unwrap();
        close(n(energy, "stored_change_j"), stored, 1e-8);
        close(n(energy, "integrated_net_input_j"), fine_energy, 2e-6);
        for (actual, expected) in temperature.iter().zip(values(&reference_field, "temperature")) {
            close(*actual, expected, 2e-6);
        }
        assert!(stored < 0.0, "cooling must remove stored heat despite the positive source");
        assert!(interface_heat(&receipt) < 0.0, "internal hot-to-cold exchange must remain live");
    }
}

#[test]
fn native_capacity_natural_contact_refuses_late_source_exit_and_retries_original_history() {
    let mut f = fixture(false, false);
    natural(&mut f, 293.0);
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_lo.value = 293.0;
    envelope.ambient_hi.value = 293.25;
    let requirement = &mut f.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 293.25;
    requirement.margin.value = 0.0;
    let policy = time(&mut f);
    policy.max_step.value = 0.25;
    policy.power_schedules = vec![TransientRegionPower {
        region: "hot".into(), source: "synthetic late pulse testing source-domain refusal and retry".into(),
        steps: [(0.5, 0.0), (1.0, 200.0)].into_iter().map(|(until, watts)| TransientPowerStep {
            until: QtyAny::new(until, dims::TIME), watts: QtyAny::new(watts, dims::POWER),
        }).collect(),
    }];
    let limited = temperature_pack_from_knots(vec![(293.0, 0.293), (293.3, 0.2933)]);
    let (source, contact, ledger) = prepare_with_pack(&f, 163, None, limited);
    let output = run(&f, &source, &contact, &ledger);
    assert_ne!(output.exit_code, fs_cli::exit::SUCCESS);
    let refusal = format!("{}{}", output.stdout, output.stderr);
    assert!(refusal.contains("cli-solve-conduction-transient"), "{refusal}");
    assert!(!refusal.contains("run_receipt"), "unsupported source state published history: {refusal}");
    let wide = temperature_pack();
    let (card, claim) = (wide.1.clone(), wide.2.clone());
    let (receipt, _) = solve_with_pack(&f, 164, wide);
    let trajectory = receipt.get("transient").unwrap();
    for grid in ["coarse", "fine"] {
        let (mut supplied, mut crossed) = (0.0, false);
        for row in trajectory.get(grid).unwrap().as_array().unwrap() {
            let (query, removed_w) = audit(row, &card, &claim, 293.0);
            let power = if n(row, "time_s") <= 0.5 { 0.0 } else { 200.0 };
            close(n(row, "source_w"), power, 1e-10);
            close(n(row, "stored_energy_change_j"), n(row, "dt_s") * (power - removed_w), 2.0 * ENERGY_TOLERANCE_J);
            if power == 0.0 {
                assert!(n(row, "stored_energy_change_j") < 0.0);
                assert!(query < INITIAL_K, "retry must begin with the original cooling history");
            }
            supplied += n(row, "dt_s") * power;
            crossed |= query > 293.3;
        }
        close(supplied, 100.0, 1e-10);
        assert!(crossed, "the wider card must admit the actual previously unsupported query");
    }
}
