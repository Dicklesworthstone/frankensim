//! Actual transported-air heat, reservoir heat and live contact share storage.

use super::*;

const INLET_K: f64 = 293.15;
const INITIAL_J: f64 = 70.0 + 2.0 * 218.5;

fn airflow_fixture() -> Fixture {
    let mut f = fixture(0.2, 0.1);
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries.iter_mut().find(|row| row.target == "cold").unwrap().condition =
        ThermalBoundaryCondition::AirflowConvection {
            branch: "cold".into(), order: 0,
            inlet_temperature: QtyAny::new(INLET_K, dims::TEMPERATURE),
            hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH),
            flow_area: QtyAny::new(0.004, dims::AREA),
            channel_length: QtyAny::new(0.3, dims::LENGTH),
            correlation: "convection.gnielinski".into(),
        };
    setup.transient.as_mut().unwrap().enthalpy.as_mut().unwrap().materials.iter_mut()
        .find(|row| row.region == "cold").unwrap().initial_specific_enthalpy.value = 70.0;
    f
}

fn surface_receipt(f: &Fixture) -> (J, J, String, String) {
    let pack = temperature_pack_from_value(PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
        samples: [[200.0, 0.5e6], [450.0, 0.5e6], [200.0, 1.5e6], [450.0, 1.5e6]]
            .into_iter().map(|p| (p, 0.001 * p[0] + 1e-8 * (p[1] - 1e6))).collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    });
    let (card, claim) = (pack.1.clone(), pack.2.clone());
    let (source, contact, ledger) = prepare_with_pack(f, 151, None, pack);
    let output = run(f, &source, &contact, &ledger);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", output.stdout, output.stderr);
    let output = J::parse(&output.stdout).unwrap();
    let completed = artifact(&ledger, output.str_field("run_receipt").unwrap());
    let stage = completed.get("stages").unwrap().as_array().unwrap().iter()
        .find(|stage| stage.str_field("stage") == Some("conduction")).unwrap();
    let conduction = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, conduction.str_field("solution_artifact").unwrap());
    (conduction, field, card, claim)
}

fn audit(conduction: &J, field: &J, source: (&str, &str), surface: bool, mixed: bool) -> Vec<f64> {
    let time = conduction.get("transient").unwrap();
    let contact_work = time.get("temperature_contact").unwrap();
    let air_work = time.get("conjugate").unwrap();
    let method = if surface {
        "source-backed-contact-enthalpy-fixed-point"
    } else { "source-backed-contact-enthalpy-newton" };
    assert_eq!(contact_work.str_field("method"), Some(method));
    assert_eq!(air_work.str_field("method"), Some("implicit-total-enthalpy-quasi-steady-airflow"));
    let mut totals = [0.0; 3];
    let mut minimum_query = f64::INFINITY;
    let mut maximum_query = f64::NEG_INFINITY;
    for grid in ["coarse", "fine"] {
        let rows = time.get(grid).unwrap().as_array().unwrap();
        assert!(rows.iter().any(|row| n(row, "time_s") == 0.375 * n(time, "final_time_s")));
        let mut input_j = 0.0;
        let mut outward_j = 0.0;
        for row in rows {
            let contact = row.get("temperature_contact").unwrap();
            let coupled = row.get("conjugate").unwrap();
            assert_eq!(contact.str_field("method"), Some(method));
            for physical in [contact, coupled] {
                assert!(n(physical, "physical_residual_norm_j") <= n(physical, "physical_residual_tolerance_j"));
                assert!(n(physical, "physical_energy_residual_j").abs() <= n(time, "energy_tolerance_j"));
            }
            close(n(contact, "physical_residual_norm_j"), n(coupled, "physical_residual_norm_j"), 0.0);
            close(n(contact, "physical_energy_residual_j"), n(coupled, "physical_energy_residual_j"), 0.0);
            assert!(n(coupled, "coupled_energy_residual_j").abs() <= n(time, "energy_tolerance_j"));
            close(n(row, "energy_residual_j"), n(coupled, "coupled_energy_residual_j"), 0.0);
            for (index, (contact_key, air_key)) in [
                ("thermal_responses", "solid_solves"),
                ("nonlinear_updates", "nonlinear_iterations"),
                ("krylov_iterations", "krylov_iterations"),
            ].into_iter().enumerate() {
                let count = n(contact, contact_key);
                close(count, n(coupled, air_key), 0.0);
                totals[index] += count;
            }
            assert!(n(contact, "thermal_responses") <= n(contact_work, "max_thermal_responses_per_step"));
            let exchange = coupled.get("exchange").unwrap();
            assert_eq!(exchange.str_field("branch"), Some("cold"));
            assert_eq!(exchange.str_field("path"), Some("vent:cold"));
            close(n(exchange, "inlet_k"), INLET_K, 0.0);
            close(n(exchange.get("air_properties").unwrap(), "specific_heat_j_kg_k"), 1007.0, 0.0);
            close(n(exchange, "mass_flow_kg_s"),
                n(exchange, "air_density_kg_m3") * n(exchange.get("flow_m3_s").unwrap(), "mid"), 2e-14);
            let capacity_rate = 1007.0 * n(exchange, "mass_flow_kg_s");
            assert!(capacity_rate > 0.0);
            let segments = exchange.get("segments").unwrap().as_array().unwrap();
            assert_eq!(segments.len(), 1);
            let segment = &segments[0];
            assert_eq!(segment.str_field("target"), Some("cold"));
            assert_eq!(segment.str_field("card"), Some("convection.gnielinski"));
            assert_eq!(segment.get("in_domain"), Some(&J::Bool(true)));
            close(n(segment, "wetted_area_m2"), 5.0, 2e-12);
            let ntu = n(segment, "htc_w_m2_k") * 5.0 / capacity_rate;
            let effectiveness = -(-ntu).exp_m1();
            close(n(segment, "ntu"), ntu, 2e-12);
            close(n(segment, "effectiveness"), effectiveness, 2e-12);
            let fluxes = coupled.get("physical_convective_robin_fluxes").unwrap().as_array().unwrap();
            let forced = fluxes.iter().find(|flux| flux.str_field("region") == Some("cold")).unwrap();
            close(n(forced, "area_m2"), 5.0, 2e-12);
            // One original branch segment obeys the exponential air law at
            // its actual area-mean wall temperature, even with latent storage.
            let air_w = capacity_rate * effectiveness
                * (n(forced, "mean_wall_temperature_k") - INLET_K);
            assert!(air_w > 0.0);
            close(capacity_rate * (n(exchange, "outlet_k") - INLET_K), air_w, 3e-6);
            close(n(exchange, "air_total_w"), air_w, 3e-6);
            close(n(segment, "air_heat_rate_w"), air_w, 3e-6);
            close(n(coupled, "air_heat_gain_w"), air_w, 3e-6);
            let balance_tolerance = n(exchange, "balance_tolerance_w");
            assert!(n(exchange, "max_region_imbalance_w").abs() <= balance_tolerance);
            assert!(n(exchange, "decomposition_residual_w").abs() <= balance_tolerance);
            close(n(forced, "heat_rate_w"), air_w, balance_tolerance.max(3e-6));
            let off_path = n(coupled, "off_path_convective_out_w");
            let radiative = n(coupled, "radiative_out_w");
            if mixed {
                let natural = row.get("natural").unwrap();
                let laws = natural.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap();
                assert_eq!(laws.len(), 1);
                assert_eq!(laws[0].str_field("target"), Some("hot"));
                close(n(&laws[0], "heat_rate_w"), off_path, 3e-6);
                assert!(off_path > 0.0);
                assert!(radiative < 0.0, "the 300 K radiation bath heats the colder latent body");
                let radiation = row.get("radiation").unwrap();
                close(n(radiation, "nonlinear_radiative_out_w"), radiative, 3e-6);
                close(n(radiation, "physical_residual_j"), n(coupled, "physical_residual_norm_j"), 0.0);
                close(n(radiation, "physical_energy_residual_j"), n(coupled, "physical_energy_residual_j"), 0.0);
                close(n(radiation, "convective_out_w"), n(coupled, "physical_convective_out_w"), 0.0);
                assert!(n(coupled, "solid_solves") >= n(natural, "solid_solves"));
            } else {
                close(off_path, 0.0, 1e-12);
                close(radiative, 0.0, 1e-12);
            }
            close(n(coupled, "physical_dirichlet_in_w"), 0.0, 0.0);
            let dt = n(row, "dt_s");
            let net_w = n(row, "source_w") - air_w - off_path - radiative;
            close(n(row, "net_input_w"), net_w, 3e-6);
            close(n(row, "stored_energy_change_j"), dt * net_w, n(time, "energy_tolerance_j"));
            input_j += dt * n(row, "source_w");
            outward_j += dt * (air_w + off_path + radiative);
            for state in contact.get("face_states").unwrap().as_array().unwrap() {
                let query = n(state, "query_temperature_k");
                close(query, 0.5 * n(state, "side_a_temperature_k") + 0.5 * n(state, "side_b_temperature_k"), 1e-12);
                close(n(state, "resistance_m2_k_w"), 0.001 * query, 1e-13);
                assert_eq!(state.str_field("interface_card"), Some(source.0));
                assert_eq!(state.str_field("selected_claim"), Some(source.1));
                minimum_query = minimum_query.min(query);
                maximum_query = maximum_query.max(query);
            }
        }
        close(input_j, 1.4, 1e-10);
        if grid == "fine" {
            let allowance = rows.len() as f64 * n(time, "energy_tolerance_j");
            close(stored(time), INITIAL_J + input_j - outward_j, allowance);
            close(n(time.get("energy").unwrap(), "integrated_net_input_j"), input_j - outward_j, allowance);
        }
    }
    for ((contact_key, air_key), total) in [
        ("total_thermal_responses", "total_solid_solves"),
        ("total_nonlinear_updates", "total_newton_updates"),
        ("total_krylov_iterations", "total_krylov_iterations"),
    ].into_iter().zip(totals) {
        close(n(contact_work, contact_key), total, 0.0);
        close(n(air_work, air_key), total, 0.0);
    }
    assert!(totals[0] > n(time, "total_steps"));
    assert!(maximum_query - minimum_query > 1e-5);
    close(n(time.get("energy").unwrap(), "stored_change_j"), stored(time) - INITIAL_J, 2e-6);
    let endpoint = time.get("endpoint").unwrap();
    close(values(endpoint, "reference_nodal_mass_kg").iter().sum(), 3.0, 1e-12);
    let h = values(endpoint, "specific_enthalpy_j_kg");
    let fraction = values(endpoint, "liquid_mass_fraction");
    let temperature = values(field, "temperature");
    let mut changed_latent = false;
    for ((&h, &fraction), &temperature) in h.iter().zip(&fraction).zip(&temperature) {
        if fraction > 1e-5 && fraction < 1.0 - 1e-5 {
            close(temperature, MELT_K, 1e-12);
            close(h, 20.0 + 100.0 * fraction, 1e-10);
            changed_latent |= (h - 70.0).abs() > 1e-5;
        }
    }
    assert!(changed_latent, "the original cold latent state must evolve once per physical step");
    h
}

#[test]
fn native_enthalpy_airflow_contact_preserves_live_source_storage_and_complete_work() {
    let mut reference: Option<Vec<f64>> = None;
    for surface in [false, true] {
        let f = airflow_fixture();
        let (conduction, field, card, claim) = if surface {
            surface_receipt(&f)
        } else {
            let (conduction, field) = receipt(&f, 150);
            let (_, card, claim) = temperature_pack();
            (conduction, field, card, claim)
        };
        let h = audit(&conduction, &field, (&card, &claim), surface, false);
        if let Some(reference) = &reference {
            for (&actual, &expected) in h.iter().zip(reference) {
                close(actual, expected, 3e-5);
            }
        } else { reference = Some(h); }
    }
}

#[test]
fn native_enthalpy_airflow_contact_keeps_natural_and_pointwise_heat_out_of_air() {
    let mut f = airflow_fixture();
    f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries
        .iter_mut().find(|row| row.target == "hot").unwrap().condition =
        ThermalBoundaryCondition::NaturalConvection {
            characteristic_length: QtyAny::new(0.1, dims::LENGTH),
            ambient_temperature: QtyAny::new(INLET_K, dims::TEMPERATURE),
            correlation: "convection.churchill-chu-vertical-plate".into(),
        };
    radiation(&mut f, 1);
    let (conduction, field) = receipt(&f, 152);
    let (_, card, claim) = temperature_pack();
    audit(&conduction, &field, (&card, &claim), false, true);
    assert_eq!(conduction.get("radiation").unwrap().str_field("spatial_model"), Some("pointwise"));
}
