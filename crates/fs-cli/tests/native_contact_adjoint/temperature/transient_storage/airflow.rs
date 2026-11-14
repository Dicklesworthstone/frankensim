//! G3: actual air enthalpy and off-path cooling share one capacity history.

use super::*;
use fs_project::{TransientPowerStep, TransientRegionPower};

const START_K: f64 = 300.0;
const AIR_CP: f64 = 1007.0;

fn n(value: &J, key: &str) -> f64 { value.f64_field(key).unwrap() }

fn airflow_fixture(enthalpy: bool, mixed: bool) -> Fixture {
    let mut f = fixture(enthalpy, false);
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.iter_mut().find(|row| row.target == "cold").unwrap().condition =
        ThermalBoundaryCondition::AirflowConvection {
            branch: "cold".into(), order: 0,
            inlet_temperature: QtyAny::new(INITIAL_K, dims::TEMPERATURE),
            hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH),
            flow_area: QtyAny::new(0.004, dims::AREA),
            channel_length: QtyAny::new(0.3, dims::LENGTH),
            correlation: "convection.gnielinski".into(),
        };
    if mixed {
        setup.boundaries.iter_mut().find(|row| row.target == "hot").unwrap().condition =
            ThermalBoundaryCondition::NaturalConvection {
                characteristic_length: QtyAny::new(0.1, dims::LENGTH),
                ambient_temperature: QtyAny::new(INITIAL_K, dims::TEMPERATURE),
                correlation: "convection.churchill-chu-vertical-plate".into(),
            };
    }
    let policy = setup.transient.as_mut().unwrap();
    policy.initial_temperature = TransientInitialTemperature::Uniform(
        QtyAny::new(START_K, dims::TEMPERATURE));
    // The enthalpy materials start at the same temperature: their chart's zero
    // is INITIAL_K, and the project refuses an initial enthalpy that does not
    // resolve to the declared uniform start.
    if let Some(storage) = policy.enthalpy.as_mut() {
        for material in &mut storage.materials {
            material.initial_specific_enthalpy.value = CP * (START_K - INITIAL_K);
        }
    }
    // Preserve the existing event times and 10 J input. Both exterior targets
    // now cool, so the 20/5/10 W holds are authored as a volume workload.
    policy.power_schedules = vec![TransientRegionPower {
        region: "hot".into(), source: "synthetic held power preserving the storage fixture's event grid and input".into(),
        steps: [(0.25, 20.0), (0.75, 5.0), (1.0, 10.0)].into_iter().map(|(until, watts)| TransientPowerStep {
            until: QtyAny::new(until, dims::TIME), watts: QtyAny::new(watts, dims::POWER),
        }).collect(),
    }];
    f
}

fn air_heat(row: &J) -> f64 {
    let coupled = row.get("conjugate").unwrap();
    let exchange = coupled.get("exchange").unwrap();
    assert_eq!(exchange.str_field("branch"), Some("cold"));
    assert_eq!(exchange.str_field("path"), Some("vent:cold"));
    close(n(exchange, "inlet_k"), INITIAL_K, 0.0);
    close(n(exchange.get("air_properties").unwrap(), "specific_heat_j_kg_k"), AIR_CP, 0.0);
    close(n(exchange, "mass_flow_kg_s"), n(exchange, "air_density_kg_m3")
        * n(exchange.get("flow_m3_s").unwrap(), "mid"), 2e-14);
    let capacity_rate = AIR_CP * n(exchange, "mass_flow_kg_s");
    assert!(capacity_rate > 0.0);
    let heat = capacity_rate * (n(exchange, "outlet_k") - INITIAL_K);
    assert!(heat > 0.0, "the actual air stream must remove stored heat");
    let segments = exchange.get("segments").unwrap().as_array().unwrap();
    assert_eq!(segments.len(), 1);
    let segment = &segments[0];
    assert_eq!(segment.str_field("target"), Some("cold"));
    assert_eq!(segment.str_field("card"), Some("convection.gnielinski"));
    assert_eq!(segment.get("in_domain"), Some(&J::Bool(true)));
    close(n(segment, "wetted_area_m2"), 5.0, 2e-12);
    let ntu = 5.0 * n(segment, "htc_w_m2_k") / capacity_rate;
    close(n(segment, "ntu"), ntu, 2e-12);
    close(n(segment, "effectiveness"), -(-ntu).exp_m1(), 2e-12);
    close(n(segment, "air_heat_rate_w"), heat, 3e-6);
    close(n(exchange, "air_total_w"), heat, 3e-6);
    close(n(coupled, "air_heat_gain_w"), heat, 3e-6);
    assert!(n(exchange, "max_region_imbalance_w").abs() <= n(exchange, "balance_tolerance_w"));
    assert!(n(exchange, "decomposition_residual_w").abs() <= n(exchange, "balance_tolerance_w"));
    close(n(coupled, "radiative_out_w"), 0.0, 0.0);
    close(n(coupled, "physical_dirichlet_in_w"), 0.0, 0.0);
    heat
}

fn natural_heat(row: &J) -> f64 {
    let natural = row.get("natural").unwrap();
    let laws = natural.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap();
    assert_eq!(laws.len(), 1);
    let law = &laws[0];
    assert_eq!(law.str_field("target"), Some("hot"));
    assert_eq!(law.str_field("card"), Some("convection.churchill-chu-vertical-plate"));
    assert_eq!(law.get("in_domain"), Some(&J::Bool(true)));
    close(n(law, "ambient_k"), INITIAL_K, 0.0);
    let wall = n(law, "mean_wall_k");
    let film = 0.5 * (wall + INITIAL_K);
    let kinematic = 1.846e-5 / (101325.0 / (287.05 * film));
    let rayleigh = 9.80665 * (wall - INITIAL_K).abs() * 0.1_f64.powi(3)
        / (film * kinematic * (kinematic / 0.707));
    let factor = 0.825 + 0.387 * rayleigh.powf(1.0 / 6.0)
        / (1.0 + (0.492_f64 / 0.707).powf(9.0 / 16.0)).powf(8.0 / 27.0);
    let h = 0.0263 * factor * factor / 0.1;
    close(n(law, "htc_w_m2_k"), h, 1e-8);
    let heat = 5.0 * h * (wall - INITIAL_K);
    assert!(heat > 0.0);
    close(n(law, "heat_rate_w"), heat, 3e-6);
    close(n(row.get("conjugate").unwrap(), "off_path_convective_out_w"), heat, 3e-6);
    heat
}

fn check_airflow(mixed: bool, surfaces: &[bool]) {
    let (reference, reference_field) = solve(&airflow_fixture(true, mixed), 170);
    let reference_time = reference.get("transient").unwrap();
    let masses = values(reference_time.get("endpoint").unwrap(), "reference_nodal_mass_kg");
    close(masses.iter().sum(), 20.0, 1e-10);
    for &surface in surfaces {
        let pack = if surface { temperature_surface_pack() } else { temperature_pack() };
        let (card, claim) = (pack.1.clone(), pack.2.clone());
        let (receipt, field) = solve_with_pack(&airflow_fixture(false, mixed), 171 + usize::from(surface), pack);
        let trajectory = receipt.get("transient").unwrap();
        let contact_work = trajectory.get("temperature_contact").unwrap();
        let air_work = trajectory.get("conjugate").unwrap();
        let method = if surface { "source-backed-contact-fixed-point" } else { "source-backed-contact-newton" };
        assert_eq!(contact_work.str_field("method"), Some(method));
        let (mut responses, mut air_iterations, mut fine_energy) = (0.0, 0.0, 0.0);
        let (mut min_query, mut max_query) = (f64::INFINITY, f64::NEG_INFINITY);
        for (grid, ends) in [
            ("coarse", vec![0.25, 0.75, 1.0]),
            ("fine", vec![0.125, 0.25, 0.5, 0.75, 0.875, 1.0]),
        ] {
            let rows = trajectory.get(grid).unwrap().as_array().unwrap();
            let reference_rows = reference_time.get(grid).unwrap().as_array().unwrap();
            assert_eq!(rows.len(), ends.len());
            assert_eq!(rows.len(), reference_rows.len());
            let (mut start, mut supplied, mut integrated, mut grid_responses) = (0.0, 0.0, 0.0, 0.0);
            for ((row, reference_row), end) in rows.iter().zip(reference_rows).zip(ends) {
                let dt = end - start;
                close(n(row, "time_s"), end, 0.0);
                close(n(row, "dt_s"), dt, 0.0);
                let query = check_contact(row, &claim);
                min_query = min_query.min(query); max_query = max_query.max(query);
                let contact = row.get("temperature_contact").unwrap();
                assert_eq!(contact.str_field("method"), Some(method));
                for state in contact.get("face_states").unwrap().as_array().unwrap() {
                    assert_eq!(state.str_field("interface_card"), Some(card.as_str()));
                }
                let owner = row.get(if mixed { "natural" } else { "conjugate" }).unwrap();
                assert!(n(owner, "physical_residual_norm_j") <= n(owner, "physical_residual_tolerance_j"));
                assert!(n(owner, "coupled_energy_residual_j").abs() <= ENERGY_TOLERANCE_J);
                for key in ["physical_residual_norm_j", "physical_energy_residual_j"] {
                    close(n(contact, key), n(owner, key), 0.0);
                }
                for (contact_key, owner_key) in [("thermal_responses", "solid_solves"),
                    ("nonlinear_updates", "nonlinear_iterations"), ("nonlinear_backtracks", "nonlinear_backtracks"),
                    ("krylov_iterations", "krylov_iterations")] {
                    close(n(contact, contact_key), n(owner, owner_key), 0.0);
                }
                let air_w = air_heat(row);
                let reservoir_w = if mixed { natural_heat(row) } else {
                    close(n(owner, "off_path_convective_out_w"), 0.0, 0.0); 0.0
                };
                if mixed { close(n(owner, "physical_convective_out_w"), air_w + reservoir_w, 3e-6); }
                let power = if start < 0.25 { 20.0 } else if start < 0.75 { 5.0 } else { 10.0 };
                let net_w = power - air_w - reservoir_w;
                close(n(row, "source_w"), power, 1e-10);
                close(n(row, "net_input_w"), net_w, 2.0 * ENERGY_TOLERANCE_J / dt + 1e-8);
                close(n(row, "stored_energy_change_j"), dt * net_w, ENERGY_TOLERANCE_J + 1e-8);
                close(n(row, "stored_energy_change_j"), n(reference_row, "stored_energy_change_j"), 4e-6);
                close(n(row, "final_region_max_k"), n(reference_row, "final_region_max_k"), 2e-6);
                grid_responses += n(contact, "thermal_responses");
                air_iterations += if mixed { n(owner, "airflow_iterations") }
                    else { n(owner.get("exchange").unwrap(), "iterations") };
                supplied += dt * power; integrated += dt * net_w; start = end;
            }
            close(supplied, 10.0, 0.0);
            close(n(contact_work, &format!("{grid}_thermal_responses")), grid_responses, 0.0);
            close(n(air_work, &format!("{grid}_solid_solves")), grid_responses, 0.0);
            responses += grid_responses;
            if grid == "fine" { fine_energy = integrated; }
        }
        close(n(contact_work, "total_thermal_responses"), responses, 0.0);
        close(n(air_work, "total_solid_solves"), responses, 0.0);
        close(n(air_work, "total_iterations"), air_iterations, 0.0);
        if mixed { close(n(trajectory.get("natural").unwrap(), "total_solid_solves"), responses, 0.0); }
        assert!(air_iterations > n(trajectory, "total_steps"));
        if surface { assert!(responses > air_iterations); } else { close(responses, air_iterations, 0.0); }
        assert!(max_query - min_query > 1e-4);
        let temperature = values(&field, "temperature");
        assert_eq!(temperature.len(), masses.len());
        let stored: f64 = temperature.iter().zip(&masses).map(|(t, mass)| mass * CP * (t - START_K)).sum();
        close(stored, fine_energy, 6.0 * ENERGY_TOLERANCE_J + 1e-8);
        close(n(trajectory.get("energy").unwrap(), "stored_change_j"), stored, 1e-8);
        close(n(air_work.get("energy").unwrap(), "integrated_net_input_j"), fine_energy, 2e-6);
        for (actual, expected) in temperature.iter().zip(values(&reference_field, "temperature")) {
            close(*actual, expected, 2e-6);
        }
        assert!(stored < 0.0);
        // Competing reservoirs can reverse which named body is warmer. The
        // unit-square contact is meshed as four equal triangles of 0.25 m2
        // (the receipt's interface area is 1 m2 over four face pairs); their
        // live trace temperatures and resistance determine the heat sign.
        let final_row = trajectory.get("fine").unwrap().as_array().unwrap().last().unwrap();
        let states = final_row.get("temperature_contact").unwrap()
            .get("face_states").unwrap().as_array().unwrap();
        assert_eq!(states.len(), 4);
        let expected_contact_w: f64 = states.iter().map(|state| {
                0.25 * (n(state, "side_a_temperature_k") - n(state, "side_b_temperature_k"))
                    / n(state, "resistance_m2_k_w")
            }).sum();
        assert!(expected_contact_w.abs() > 1e-8);
        close(interface_heat(&receipt), expected_contact_w, 1e-8);
    }
}

#[test]
fn native_capacity_airflow_contact_matches_actual_air_heat_and_storage() {
    check_airflow(false, &[false, true]);
}

#[test]
fn native_capacity_mixed_contact_keeps_reservoir_heat_out_of_the_air_path() {
    check_airflow(true, &[false]);
}
