//! Native capacity storage joins the same enclosure and convection owners.
//! The sensible-enthalpy comparison is a separately authored numerical oracle;
//! no density or enthalpy conversion is used by the capacity product path.

use super::*;

const CV: f64 = 60.0;

fn capacity(spec: &mut ProjectSpec) {
    let policy = time(spec);
    policy.enthalpy = None;
    policy.capacities = vec![TransientRegionCapacity {
        region: "air".into(),
        volumetric_heat_capacity: QtyAny::new(CV, dims::VOLUMETRIC_HEAT_CAPACITY),
        source: "independently declared synthetic 60 J/(m3 K) capacity".into(),
    }];
}

fn material_project(kind: &str) -> (ProjectSpec, CardPackSet) {
    let base = transient_product::transient_conductivity_cards(kind);
    let (cards, card) = radiation_product::emissivity_cards(&base, EPSILON);
    let mut spec = project(&card);
    rebind_to(&mut spec, &base);
    let binding = &mut spec.materials.as_mut().unwrap()[0];
    binding.temp_lo.value = 280.0;
    binding.temp_hi.value = 360.0;
    let identity = fs_blake3::ContentHash::from_hex(&binding.card).unwrap();
    let material = &mut time(&mut spec).enthalpy.as_mut().unwrap().materials[0];
    material.material_card_identity = identity;
    // This chart independently declares rho=6, Cp=10; their product matches
    // CV only inside the regression oracle, not inside native lowering.
    material.knots = knots(&[
        (10.0 * (280.0 - INITIAL_K), 280.0, 0.0),
        (10.0 * (360.0 - INITIAL_K), 360.0, 0.0),
    ]);
    (spec, cards)
}

fn check_capacity(solved: &Solved, single_tet: bool) {
    let receipt = J::parse(&solved.conduction).unwrap();
    let field = J::parse(&solved.field).unwrap();
    let history = receipt.get("transient").unwrap();
    assert_eq!(history.str_field("schema"), Some("fs-cli-transient-conduction-v1"));
    assert_eq!(history.str_field("authority"), Some("Estimated"));
    assert_eq!(history.get("radiation").unwrap().str_field("method"),
        Some("implicit-closed-gray-diffuse-enclosure"));
    close(n(history, "total_steps"), 9.0, 0.0);
    close(n(history.get("workload").unwrap(), "integrated_scheduled_input_j"), 9.0, 1e-12);
    let mut fine_storage = 0.0;
    for grid in ["coarse", "fine"] {
        let rows = history.get(grid).unwrap().as_array().unwrap();
        for switch in [0.3, 0.7, 1.0] {
            assert!(rows.iter().any(|row| n(row, "time_s") == switch));
        }
        let (mut start, mut storage, mut input) = (0.0, 0.0, 0.0);
        for row in rows {
            let dt = n(row, "dt_s");
            close(dt, n(row, "time_s") - start, 1e-15);
            let watts = if start < 0.3 { 10.0 } else if start < 0.7 { 0.0 } else { 20.0 };
            close(n(row, "source_w"), watts, 1e-10);
            let radiation = row.get("radiation").unwrap();
            assert_eq!(radiation.str_field("model"), Some("closed-gray-diffuse-enclosure"));
            assert_eq!(radiation.str_field("storage"), Some("volumetric-heat-capacity"));
            assert!(n(radiation, "solid_solves") >= 1.0);
            assert!(n(radiation, "physical_residual_j")
                <= n(radiation, "physical_residual_threshold_j"));
            assert!(n(radiation, "physical_energy_residual_j").abs() <= n(history, "energy_tolerance_j"));
            close(n(radiation, "applied_radiation_out_w"), 0.0, 1e-6);
            close(n(radiation, "nonlinear_radiation_out_w"), 0.0, 1e-6);
            close(n(radiation, "prescribed_neumann_out_w"), 0.0, 1e-12);
            close(n(radiation, "physical_dirichlet_in_w"), 0.0, 1e-12);
            close(n(row, "net_input_w"), watts - n(radiation, "convective_out_w"), 2e-6);
            close(n(row, "stored_energy_change_j") - dt * n(row, "net_input_w"),
                n(row, "energy_residual_j"), 2e-10);
            assert!(n(row, "energy_residual_j").abs() <= n(history, "energy_tolerance_j"));
            storage += n(row, "stored_energy_change_j");
            input += dt * n(row, "net_input_w");
            start = n(row, "time_s");
        }
        close(storage, input, 1e-5);
        if grid == "fine" { fine_storage = storage; }
    }
    let energy = history.get("energy").unwrap();
    close(n(energy, "stored_change_j"), fine_storage, 2e-9);
    close(n(energy, "integrated_net_input_j"), fine_storage, 1e-5);
    if single_tet {
        let temperature = values(&field, "temperature");
        assert_eq!(temperature.len(), 4);
        let actual: f64 = temperature.iter().map(|t| CV / 24.0 * (t - INITIAL_K)).sum();
        close(fine_storage, actual, 2e-8);
    }
    let radiation = receipt.get("radiation").unwrap();
    assert_eq!(radiation.str_field("time_integration"), Some("backward-euler-capacity"));
    let surfaces = radiation.get("surfaces").unwrap().as_array().unwrap();
    let a = n(&surfaces[0], "mean_temperature_k");
    let b = n(&surfaces[1], "mean_temperature_k");
    let exchange = SIGMA * 0.5 * (a - b) * (a + b) * (a*a + b*b) / (2.0 / EPSILON - 1.0);
    assert!(exchange.abs() > 1e-5, "fixture must exchange reflected radiation");
    close(n(&surfaces[0], "nonlinear_heat_w"), exchange, 1e-6);
    close(n(&surfaces[1], "nonlinear_heat_w"), -exchange, 1e-6);
    let final_energy = receipt.get("energy").unwrap();
    close(n(final_energy, "source_w") + n(final_energy, "dirichlet_in_w")
        - n(final_energy, "neumann_out_w") - n(final_energy, "robin_out_w")
        - n(final_energy, "storage_w"), n(final_energy, "closure_w"), 2e-10);
    close(n(final_energy, "closure_w"), 0.0, 2e-5);
    assert!(solved.html.contains("closed-gray-diffuse-enclosure"));
    assert!(solved.html.contains("Estimated"));
}

#[test]
fn g1_g3_g5_native_capacity_enclosure_matches_authored_sensible_enthalpy_and_replays() {
    let mut endpoints = Vec::new();
    for kind in ["constant-scalar", "nonlinear"] {
        let (sensible, cards) = material_project(kind);
        let enthalpy = run(&sensible, &cards);
        let mut spec = sensible.clone();
        capacity(&mut spec);
        let solved = run(&spec, &cards);
        check_capacity(&solved, true);
        let expected = values(&J::parse(&enthalpy.field).unwrap(), "temperature");
        let actual = values(&J::parse(&solved.field).unwrap(), "temperature");
        assert_eq!(actual.len(), expected.len());
        for (&a, &b) in actual.iter().zip(&expected) { close(a, b, 2e-6); }
        let receipt = J::parse(&solved.conduction).unwrap();
        let history = receipt.get("transient").unwrap();
        if kind == "nonlinear" {
            assert!(n(history.get("nonlinear").unwrap(), "total_updates") > 0.0);
            assert_eq!(solved, run(&spec, &cards), "retained artifacts replay exactly");
        }
        endpoints.push(actual);
    }
    assert!(endpoints[0].iter().zip(&endpoints[1]).any(|(a,b)| (a-b).abs() > 1e-4),
        "temperature-dependent conductivity must change this nonuniform history");
}

fn natural() -> ThermalBoundaryCondition {
    ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(290.0, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    }
}

fn air() -> ThermalBoundaryCondition {
    ThermalBoundaryCondition::AirflowConvection {
        branch: "air".into(), order: 0,
        inlet_temperature: QtyAny::new(290.0, dims::TEMPERATURE),
        hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH),
        flow_area: QtyAny::new(0.004, dims::AREA),
        channel_length: QtyAny::new(0.3, dims::LENGTH),
        correlation: "convection.gnielinski".into(),
    }
}

fn independent_natural_htc(wall: f64) -> f64 {
    let film = 0.5 * (wall + 290.0);
    let density = 101325.0 / (287.05 * film);
    let viscosity = 1.846e-5 / density;
    let rayleigh = 9.80665 * (wall - 290.0).abs() * 0.1_f64.powi(3)
        / (film * viscosity * (viscosity / 0.707));
    assert!((0.1..=1e12).contains(&rayleigh));
    let nu = (0.825 + 0.387 * rayleigh.powf(1.0/6.0)
        / (1.0 + (0.492_f64 / 0.707).powf(9.0/16.0)).powf(8.0/27.0)).powi(2);
    nu * 0.0263 / 0.1
}

#[test]
fn g1_native_capacity_enclosure_composes_actual_natural_and_airflow_with_nonlinear_storage_history() {
    for variant in ["natural", "airflow", "mixed"] {
        let (mut spec, cards) = material_project("nonlinear");
        capacity(&mut spec);
        spec.envelope.as_mut().unwrap().ambient_lo.value = 290.0;
        spec.envelope.as_mut().unwrap().ambient_hi.value = 310.0;
        let boundaries = &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries;
        boundaries[0].condition = if variant == "airflow" { air() } else { natural() };
        if variant == "mixed" { boundaries[1].condition = air(); }
        let solved = run(&spec, &cards);
        check_capacity(&solved, true);
        let receipt = J::parse(&solved.conduction).unwrap();
        let history = receipt.get("transient").unwrap();
        assert!(n(history.get("nonlinear").unwrap(), "total_updates") > 0.0);
        for row in history.get("fine").unwrap().as_array().unwrap() {
            let dt = n(row, "dt_s");
            let outer = row.get(if variant == "airflow" { "conjugate" } else { "natural" }).unwrap();
            assert!(n(outer, "physical_residual_norm_j") <= n(outer, "physical_residual_tolerance_j"));
            assert!(n(outer, "physical_energy_residual_j").abs() <= n(history, "energy_tolerance_j"));
            assert!(n(outer, "coupled_energy_residual_j").abs() <= n(history, "energy_tolerance_j"));
            if variant != "airflow" {
                let laws = outer.get("exchange").unwrap().get("laws").unwrap().as_array().unwrap();
                assert_eq!(laws.len(), 1);
                let wall = n(&laws[0], "mean_wall_k");
                let h = independent_natural_htc(wall);
                close(n(&laws[0], "htc_w_m2_k"), h, 2e-8);
                close(n(&laws[0], "heat_rate_w"), h * 0.5 * (wall - 290.0), 2e-7);
                close(n(row, "stored_energy_change_j")
                    - dt * n(outer, "coupled_net_input_w"),
                    n(outer, "coupled_energy_residual_j"), 2e-10);
            }
            if variant != "natural" {
                let air = row.get("conjugate").unwrap();
                let net = n(row, "source_w") + n(air, "physical_dirichlet_in_w")
                    - n(air, "air_heat_gain_w") - n(air, "off_path_convective_out_w")
                    - n(air, "radiative_out_w");
                close(n(row, "stored_energy_change_j") - dt * net,
                    n(air, "coupled_energy_residual_j"), 2e-7);
                assert!(n(air, "air_heat_gain_w") > 0.0);
                assert!(n(air, "radiation_trials") >= n(air, "solid_solves"));
                assert!(n(row, "krylov_iterations") >= n(air, "krylov_iterations"));
                if variant == "mixed" {
                    assert!(n(outer, "solid_solves") >= n(air, "solid_solves"));
                    assert!(n(outer, "airflow_iterations") >= 1.0);
                }
            }
        }
    }
}

#[test]
fn g3_native_capacity_enclosure_spatial_ladder_shares_time_work_and_rebinds_patches() {
    let (mut spec, cards) = material_project("nonlinear");
    capacity(&mut spec);
    spec.solver.as_mut().unwrap().fidelity = "ladder".into();
    time(&mut spec).max_steps = 18;
    let solved = run(&spec, &cards);
    check_capacity(&solved, false);
    let receipt = J::parse(&solved.conduction).unwrap();
    let ladder = receipt.get("ladder").unwrap();
    let rungs = ladder.get("rungs").unwrap().as_array().unwrap();
    assert_eq!(rungs.len(), 2);
    close(n(&rungs[1], "tets"), 8.0 * n(&rungs[0], "tets"), 0.0);
    assert_eq!(ladder.str_field("stop"), Some("time-step-budget"));
    close(n(ladder.get("space_time").unwrap(), "forward_steps_completed"), 18.0, 0.0);
    close(n(ladder.get("space_time").unwrap(), "max_forward_steps"), 18.0, 0.0);
    let radiation = receipt.get("radiation").unwrap();
    for surface in radiation.get("surfaces").unwrap().as_array().unwrap() {
        close(n(surface, "area_m2"), 0.5, 2e-12);
    }
}

