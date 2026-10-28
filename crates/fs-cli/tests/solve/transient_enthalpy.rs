//! G1/G3 native phase-changing heat storage with independent energy oracles.

use super::*;
use fs_project::{
    TransientEnthalpy, TransientEnthalpyKnot, TransientPowerStep, TransientRegionEnthalpy,
    TransientRegionPower, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} differs from {expected}, tolerance {tolerance}"
    );
}

fn n(value: &J, key: &str) -> f64 {
    value.f64_field(key).unwrap()
}

fn time(spec: &mut ProjectSpec) -> &mut fs_project::ConductionTransient {
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .transient
        .as_mut()
        .unwrap()
}

fn project(density: f64, horizon: f64, knots: &[(f64, f64, f64)]) -> ProjectSpec {
    let mut spec = conduction_fixture_project(7, &tetra_stl());
    transient_product::declaration(&mut spec, 100.0, horizon);
    let identity = fs_blake3::ContentHash::from_hex(&fixture_card_identity().0).unwrap();
    time(&mut spec).capacities.clear();
    time(&mut spec).enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            phase_conductivity: None,
            region: "air".into(),
            material_card_identity: identity,
            source: "synthetic enthalpy relation for an independent energy oracle; not a material-property claim".into(),
            reference_density: QtyAny::new(density, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            knots: knots.iter().map(|&(h, temperature, fraction)| TransientEnthalpyKnot {
                specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                liquid_mass_fraction: fraction,
            }).collect(),
        }],
    });
    spec.cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries[0]
        .condition = ThermalBoundaryCondition::HeatFlux {
        outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
    };
    spec
}

fn run(spec: &ProjectSpec) -> (String, J, J, J, Ledger) {
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, tetra_stl());
    let decoded = decode(spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let (run, qoi, conduction, field) =
        run_conjugate_to_completion(&ledger, &decoded, &fixture_cards());
    (
        run,
        J::parse(&qoi).unwrap(),
        J::parse(&conduction).unwrap(),
        J::parse(&field).unwrap(),
        ledger,
    )
}

fn values<'a>(row: &'a J, name: &str) -> impl Iterator<Item = f64> + 'a {
    row.get(name)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
}

#[test]
fn g1_native_enthalpy_constant_cp_recovers_mass_scaled_heating_and_deterministic_replay() {
    let knots = [(-90.0, 203.15, 0.0), (100.0, 393.15, 0.0)];
    let mut first = None;
    for density in [6.0, 12.0, 6.0] {
        let spec = project(density, 2.0, &knots);
        let (run, qoi, conduction, field, _) = run(&spec);
        // Unit right tetrahedron V=1/6 m³, P=5 W, Cp=1 J/(kg K).
        let expected_h = 10.0 / (density / 6.0);
        let expected_temperature = 293.15 + expected_h;
        let time = conduction.get("transient").unwrap();
        assert_eq!(
            time.str_field("schema"),
            Some("fs-cli-transient-enthalpy-v1")
        );
        assert_eq!(time.str_field("authority"), Some("Estimated"));
        let endpoint = time.get("endpoint").unwrap();
        for value in values(endpoint, "specific_enthalpy_j_kg") {
            close(value, expected_h, 2e-8);
        }
        for value in values(endpoint, "liquid_mass_fraction") {
            close(value, 0.0, 0.0);
        }
        for value in values(&field, "temperature") {
            close(value, expected_temperature, 2e-8);
        }
        close(
            values(endpoint, "reference_nodal_mass_kg").sum(),
            density / 6.0,
            2e-12,
        );
        close(
            n(time.get("energy").unwrap(), "stored_change_j"),
            10.0,
            2e-7,
        );
        let qoi_value = &qoi.get("qoi").unwrap().as_array().unwrap()[0];
        close(n(qoi_value, "value"), expected_temperature, 2e-8);
        assert!(
            conduction
                .str_field("no_claim")
                .unwrap()
                .contains("total-enthalpy")
        );
        if density == 6.0 {
            let retained = (run, qoi, conduction, field);
            if let Some(original) = &first {
                assert_eq!(original, &retained);
            } else {
                first = Some(retained);
            }
        }
    }
}

#[test]
fn g1_native_enthalpy_prescribed_bath_retains_initial_storage_and_reservoir_heat() {
    const DENSITY: f64 = 6.0;
    const CP: f64 = 10.0;
    const INITIAL_K: f64 = 293.15;
    const HORIZON: f64 = 2.0;
    let knots = [(-200.0, 273.15, 0.0), (200.0, 313.15, 0.0)];
    for bath_k in [283.15, 303.15] {
        let mut spec = project(DENSITY, HORIZON, &knots);
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0]
            .condition = ThermalBoundaryCondition::FixedTemperature {
                temperature: QtyAny::new(bath_k, dims::TEMPERATURE),
            };
        let (_, qoi, conduction, field, _) = run(&spec);
        let history = conduction.get("transient").unwrap();
        let expected_h = CP * (bath_k - INITIAL_K);
        let stored_j = (DENSITY / 6.0) * expected_h;
        close(n(conduction.get("solver").unwrap(), "free_dofs"), 0.0, 0.0);
        close(n(history, "initial_temperature_k"), INITIAL_K, 0.0);
        for temperature in values(&field, "temperature") {
            close(temperature, bath_k, 2e-12);
        }
        for h in values(history.get("endpoint").unwrap(), "specific_enthalpy_j_kg") {
            close(h, expected_h, 2e-11);
        }
        let qoi_value = &qoi.get("qoi").unwrap().as_array().unwrap()[0];
        close(n(qoi_value, "value"), bath_k, 2e-12);
        for grid in ["coarse", "fine"] {
            let rows = history.get(grid).unwrap().as_array().unwrap();
            let mut reservoir_j = 0.0;
            for (index, row) in rows.iter().enumerate() {
                let dt = n(row, "dt_s");
                let step_storage = if index == 0 { stored_j } else { 0.0 };
                let reaction_w = step_storage / dt - 5.0;
                close(n(row, "stored_energy_change_j"), step_storage, 2e-10);
                close(n(row, "dirichlet_in_w"), reaction_w, 2e-9);
                close(n(row, "physical_dirichlet_in_w"), reaction_w, 2e-9);
                close(n(row, "net_input_w"), step_storage / dt, 2e-9);
                close(n(row, "energy_residual_j"), 0.0, 2e-10);
                reservoir_j += dt * n(row, "physical_dirichlet_in_w");
            }
            // Initial history remains at 293.15 K. Its complete jump is
            // charged to the bath once, even though every endpoint is fixed.
            close(reservoir_j, stored_j - 5.0 * HORIZON, 2e-9);
        }
        let energy = history.get("energy").unwrap();
        close(n(energy, "stored_change_j"), stored_j, 2e-10);
        close(n(energy, "integrated_net_input_j"), stored_j, 2e-9);
        close(n(energy, "window_residual_j"), 0.0, 2e-9);
        let endpoint = conduction.get("energy").unwrap();
        close(n(endpoint, "dirichlet_in_w"), -5.0, 2e-10);
        close(n(endpoint, "storage_w"), 0.0, 2e-10);
        close(n(endpoint, "closure_w"), 0.0, 2e-10);
    }
}

#[test]
fn g1_native_enthalpy_constant_cp_matches_two_body_contact_and_conserves_total_energy() {
    const DENSITY: f64 = 10.0;
    const CP: f64 = 10.0;
    const INITIAL_T: f64 = 293.15;
    const HORIZON: f64 = 2.0;
    const DELIVERED_J: f64 = 5.0 * HORIZON;

    let mut capacity_spec = multi_region_contact_project();
    transient_product::declaration(&mut capacity_spec, DENSITY * CP, HORIZON);
    capacity_spec.requirements.as_mut().unwrap()[0].region = "hot".into();
    for boundary in &mut capacity_spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap()
        .boundaries
    {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX),
        };
    }
    let mut enthalpy_spec = capacity_spec.clone();
    time(&mut enthalpy_spec).capacities.clear();
    time(&mut enthalpy_spec).enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: ["cold", "hot"]
            .into_iter()
            .map(|region| TransientRegionEnthalpy {
                phase_conductivity: None,
                region: region.into(),
                material_card_identity: fs_blake3::ContentHash::from_hex(
                    &fixture_card_identity().0,
                )
                .unwrap(),
                source:
                    "synthetic Cp = 10 J/(kg K) contact comparison; not a material-property claim"
                        .into(),
                reference_density: QtyAny::new(DENSITY, dims::DENSITY),
                initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
                knots: [-50.0, 0.0, 50.0]
                    .into_iter()
                    .map(|delta_t| TransientEnthalpyKnot {
                        specific_enthalpy: QtyAny::new(CP * delta_t, dims::SPECIFIC_ENTHALPY),
                        temperature: QtyAny::new(INITIAL_T + delta_t, dims::TEMPERATURE),
                        liquid_mass_fraction: 0.0,
                    })
                    .collect(),
            })
            .collect(),
    });
    let solve = |spec: &ProjectSpec| {
        let ledger = Ledger::open(":memory:").unwrap();
        import_multi_region_contact(&ledger, spec);
        let (_, _, conduction, field) =
            run_conjugate_to_completion(&ledger, &decode(spec), &contact_cards());
        (J::parse(&conduction).unwrap(), J::parse(&field).unwrap())
    };
    let (capacity, capacity_field) = solve(&capacity_spec);
    let (enthalpy, enthalpy_field) = solve(&enthalpy_spec);
    let capacity_time = capacity.get("transient").unwrap();
    let enthalpy_time = enthalpy.get("transient").unwrap();
    assert_eq!(
        capacity_time.str_field("schema"),
        Some("fs-cli-transient-conduction-v1")
    );
    assert_eq!(
        enthalpy_time.str_field("schema"),
        Some("fs-cli-transient-enthalpy-v1")
    );

    // Both P1 formulations use the same mesh and lumped storage. Comparing
    // every retained node exercises spatial transport and the old-state march,
    // while h = Cp (T - T0) independently fixes the enthalpy gauge and units.
    let temperatures: Vec<_> = values(&enthalpy_field, "temperature").collect();
    let capacity_temperatures: Vec<_> = values(&capacity_field, "temperature").collect();
    let endpoint = enthalpy_time.get("endpoint").unwrap();
    let enthalpies: Vec<_> = values(endpoint, "specific_enthalpy_j_kg").collect();
    let masses: Vec<_> = values(endpoint, "reference_nodal_mass_kg").collect();
    assert!(!temperatures.is_empty());
    assert_eq!(temperatures.len(), capacity_temperatures.len());
    assert_eq!(temperatures.len(), enthalpies.len());
    assert_eq!(temperatures.len(), masses.len());
    for (index, &temperature) in temperatures.iter().enumerate() {
        close(temperature, capacity_temperatures[index], 2e-8);
        close(enthalpies[index], CP * (temperature - INITIAL_T), 2e-9);
        close(
            enthalpies[index],
            CP * (capacity_temperatures[index] - INITIAL_T),
            2e-7,
        );
        assert!(masses[index] > 0.0);
    }
    assert!(
        temperatures
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
            - temperatures.iter().copied().fold(f64::INFINITY, f64::min)
            > 0.01,
        "heating only one body must produce a resolved nonuniform temperature"
    );
    for fraction in values(endpoint, "liquid_mass_fraction") {
        close(fraction, 0.0, 0.0);
    }
    // The two unit cubes contain 20 kg. Contact redistributes heat internally:
    // summing physical nodal storage must still recover exactly 5 W * 2 s.
    close(masses.iter().sum(), 2.0 * DENSITY, 2e-12);
    close(
        masses.iter().zip(&enthalpies).map(|(m, h)| m * h).sum(),
        DELIVERED_J,
        2e-7,
    );
    close(
        masses
            .iter()
            .zip(&capacity_temperatures)
            .map(|(m, t)| m * CP * (t - INITIAL_T))
            .sum(),
        DELIVERED_J,
        2e-7,
    );
    for history in [capacity_time, enthalpy_time] {
        let energy = history.get("energy").unwrap();
        close(n(energy, "stored_change_j"), DELIVERED_J, 2e-7);
        close(n(energy, "integrated_net_input_j"), DELIVERED_J, 2e-12);
        close(n(energy, "window_residual_j"), 0.0, 2e-7);
    }
    for grid in ["coarse", "fine"] {
        let capacity_rows = capacity_time.get(grid).unwrap().as_array().unwrap();
        let enthalpy_rows = enthalpy_time.get(grid).unwrap().as_array().unwrap();
        assert_eq!(capacity_rows.len(), enthalpy_rows.len());
        for (a, b) in capacity_rows.iter().zip(enthalpy_rows) {
            close(n(a, "time_s"), n(b, "time_s"), 0.0);
            close(n(a, "final_region_max_k"), n(b, "final_region_max_k"), 2e-8);
            for row in [a, b] {
                close(n(row, "source_w"), 5.0, 2e-12);
                close(n(row, "stored_energy_change_j"), 5.0 * n(row, "dt_s"), 2e-7);
            }
        }
    }

    let capacity_fluxes = capacity
        .path(&["interfaces", "fluxes"])
        .unwrap()
        .as_array()
        .unwrap();
    let enthalpy_fluxes = enthalpy
        .path(&["interfaces", "fluxes"])
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(capacity_fluxes.len(), 1);
    assert_eq!(enthalpy_fluxes.len(), 1);
    let a = &capacity_fluxes[0];
    let b = &enthalpy_fluxes[0];
    assert_eq!(a.str_field("interface"), Some("cold-hot-joint"));
    assert_eq!(a.str_field("interface"), b.str_field("interface"));
    assert_eq!(a.str_field("card"), b.str_field("card"));
    for key in [
        "area_m2",
        "conductance_w_per_k",
        "mean_jump_k",
        "heat_rate_a_to_b_w",
    ] {
        close(n(a, key), n(b, key), 2e-7);
    }
    for flux in [a, b] {
        close(n(flux, "area_m2"), 1.0, 2e-12);
        close(n(flux, "conductance_w_per_k"), 10.0, 2e-10);
        close(
            n(flux, "heat_rate_a_to_b_w"),
            10.0 * n(flux, "mean_jump_k"),
            2e-12,
        );
        // The interface is oriented cold -> hot, so positive hot -> cold heat
        // has a negative recorded sign. The 5 W source is still partly stored.
        let hot_to_cold_w = -n(flux, "heat_rate_a_to_b_w");
        assert!(hot_to_cold_w > 0.01 && hot_to_cold_w < 5.0);
    }
}

#[test]
fn g1_native_enthalpy_resolves_latent_storage_instead_of_heating_through_the_plateau() {
    let knots = [
        (-10.0, 283.15, 0.0),
        (5.3, 298.45, 0.0),
        (25.3, 298.45, 1.0),
        (125.0, 398.15, 1.0),
    ];
    let spec = project(6.0, 4.0, &knots);
    let (run, _, conduction, field, ledger) = run(&spec);
    let time = conduction.get("transient").unwrap();
    let endpoint = time.get("endpoint").unwrap();
    // The reference mass is 1 kg: 5 W * 4 s -> 20 J/kg.
    for value in values(endpoint, "specific_enthalpy_j_kg") {
        close(value, 20.0, 3e-8);
    }
    for value in values(endpoint, "liquid_mass_fraction") {
        close(value, (20.0 - 5.3) / 20.0, 2e-9);
    }
    for value in values(&field, "temperature") {
        close(value, 298.45, 2e-9);
    }
    close(
        n(time.get("energy").unwrap(), "stored_change_j"),
        20.0,
        2e-7,
    );
    close(
        n(time.get("energy").unwrap(), "integrated_net_input_j"),
        20.0,
        2e-10,
    );
    assert_eq!(
        time.get("temporal_error")
            .unwrap()
            .get("estimated_half_width_k"),
        Some(&J::Null),
        "identical plateau temperatures must not claim zero temporal uncertainty"
    );
    let receipts = stage_receipt_hashes(&ledger, &run);
    let report = String::from_utf8(artifact_bytes(&ledger, &receipts[6])).unwrap();
    let html = String::from_utf8(artifact_bytes(
        &ledger,
        &receipt_str_field(&report, "report_html"),
    ))
    .unwrap();
    assert!(
        html.contains("total-enthalpy"),
        "the native report must carry the real storage model"
    );
    assert!(
        html.contains("final time 4"),
        "the reported QoI must remain a final-time observation"
    );
}

#[test]
fn g1_native_enthalpy_workload_switches_deliver_exact_pulse_energy() {
    let knots = [
        (-10.0, 283.15, 0.0),
        (1.0, 294.15, 0.0),
        (21.0, 294.15, 1.0),
        (121.0, 394.15, 1.0),
    ];
    let mut spec = project(6.0, 2.0, &knots);
    time(&mut spec).power_schedules = vec![TransientRegionPower {
        region: "air".into(),
        source: "declared 0.3-second pulse at 10 W".into(),
        steps: vec![(0.3, 10.0), (2.0, 0.0)]
            .into_iter()
            .map(|(until, watts)| TransientPowerStep {
                until: QtyAny::new(until, dims::TIME),
                watts: QtyAny::new(watts, dims::POWER),
            })
            .collect(),
    }];
    let (_, _, conduction, field, _) = run(&spec);
    let time = conduction.get("transient").unwrap();
    for grid in ["coarse", "fine"] {
        let rows = time.get(grid).unwrap().as_array().unwrap();
        assert!(rows.iter().any(|row| n(row, "time_s") == 0.3));
        let energy: f64 = rows
            .iter()
            .map(|row| n(row, "dt_s") * n(row, "source_w"))
            .sum();
        close(energy, 3.0, 2e-12);
    }
    let endpoint = time.get("endpoint").unwrap();
    for value in values(endpoint, "specific_enthalpy_j_kg") {
        close(value, 3.0, 2e-8);
    }
    for value in values(endpoint, "liquid_mass_fraction") {
        close(value, 0.1, 2e-9);
    }
    for value in values(&field, "temperature") {
        close(value, 294.15, 2e-9);
    }
    close(n(time.get("energy").unwrap(), "stored_change_j"), 3.0, 2e-7);
}

#[test]
fn g1_native_enthalpy_shipped_example_runs_through_the_existing_report_stage() {
    let decoded = fs_project::parse_sexpr(include_str!(
        "../../../../data/reference-project/cooling-enthalpy.fsim"
    ))
    .unwrap();
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(
        &ledger,
        &decoded.spec,
        include_bytes!("../../../../data/reference-project/plate.stl").to_vec(),
    );
    let (_, _, conduction, _) = run_conjugate_to_completion(&ledger, &decoded, &fixture_cards());
    let parsed = J::parse(&conduction).unwrap();
    let time = parsed.get("transient").unwrap();
    close(
        n(time.get("energy").unwrap(), "stored_change_j"),
        10.0,
        2e-7,
    );
    for fraction in values(time.get("endpoint").unwrap(), "liquid_mass_fraction") {
        close(fraction, 2.0 / 3.0, 2e-8);
    }
}

#[test]
fn g0_native_enthalpy_refusals_leave_no_published_conduction_field() {
    for variant in ["domain", "latent-temperature-boundary", "outside-temperature-boundary"] {
        let mut spec = if variant == "latent-temperature-boundary" {
            project(6.0, 2.0, &[
                (-10.0, 283.15, 0.0), (0.0, 293.15, 0.0),
                (10.0, 293.15, 1.0), (20.0, 303.15, 1.0),
            ])
        } else {
            project(6.0, 2.0, &[(-1.0, 292.15, 0.0), (2.0, 295.15, 0.0)])
        };
        if variant != "domain" {
            spec.cooling
                .as_mut()
                .unwrap()
                .conduction
                .as_mut()
                .unwrap()
                .boundaries[0]
                .condition = ThermalBoundaryCondition::FixedTemperature {
                temperature: QtyAny::new(
                    if variant == "latent-temperature-boundary" { 293.15 } else { 300.0 },
                    dims::TEMPERATURE,
                ),
            };
        }
        let ledger = Ledger::open(":memory:").unwrap();
        import_fixture(&ledger, &spec, tetra_stl());
        let error = run_solve(
            &ledger,
            &CancelGate::new_clock_free(),
            &mut benign_clock(),
            &decode(&spec),
            &fixture_cards(),
            &mut Vec::new(),
        )
        .unwrap_err();
        assert_eq!(
            error.code, "cli-solve-conduction-transient",
            "{variant}: {error:?}"
        );
        assert_eq!(
            stage_receipt_hashes(&ledger, error.run.as_ref().unwrap()).len(),
            4,
            "{variant}: an incomplete physical history must publish no conduction stage"
        );
    }
}


#[test]
fn g1_native_enthalpy_spatial_ladder_preserves_mass_energy_and_shared_time_allowance() {
    let mut spec = project(6.0, 2.0, &[(-90.0, 203.15, 0.0), (100.0, 393.15, 0.0)]);
    spec.solver.as_mut().unwrap().fidelity = "ladder".into();
    time(&mut spec).max_steps = 36;
    let (_, qoi, conduction, field, _) = run(&spec);
    let ladder = conduction.get("ladder").unwrap();
    assert_eq!(ladder.get("rungs").unwrap().as_array().unwrap().len(), 3);
    close(n(ladder.get("space_time").unwrap(), "forward_steps_completed"), 36.0, 0.0);
    let history = conduction.get("transient").unwrap();
    assert_eq!(history.str_field("schema"), Some("fs-cli-transient-enthalpy-v1"));
    close(n(history, "total_steps"), 12.0, 0.0);
    close(n(history.get("energy").unwrap(), "stored_change_j"), 10.0, 2e-6);
    let endpoint = history.get("endpoint").unwrap();
    close(values(endpoint, "reference_nodal_mass_kg").sum(), 1.0, 2e-12);
    for h in values(endpoint, "specific_enthalpy_j_kg") { close(h, 10.0, 2e-7); }
    for temperature in values(&field, "temperature") { close(temperature, 303.15, 2e-7); }
    close(n(&qoi.get("qoi").unwrap().as_array().unwrap()[0], "value"), 303.15, 2e-7);
}
