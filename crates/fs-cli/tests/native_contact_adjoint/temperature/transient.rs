//! Native imported contacts evolve their original source law in physical time.

use super::*;
use fs_project::{
    ConductionTransient, ThermalBoundaryCondition, TransientInitialTemperature, TransientPowerStep,
    TransientRegionCapacity, TransientRegionPower, spec::dims,
};
use fs_qty::QtyAny;

#[path = "transient/adjoint.rs"]
mod adjoint;

fn fixture(nonlinear: bool, max_step: f64) -> Fixture {
    let mut f = Fixture::new(nonlinear, false);
    let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient = Some(ConductionTransient {
        initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(293.15, dims::TEMPERATURE)),
        horizon: QtyAny::new(2.0, dims::TIME),
        max_step: QtyAny::new(max_step, dims::TIME),
        max_steps: 256,
        energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
        boundary_schedules: Vec::new(),
        power_schedules: vec![TransientRegionPower {
            region: "hot".into(),
            source: "synthetic delivered-power pulse; no empirical claim".into(),
            steps: vec![
                TransientPowerStep { until: QtyAny::new(0.75, dims::TIME), watts: QtyAny::new(12.0, dims::POWER) },
                TransientPowerStep { until: QtyAny::new(2.0, dims::TIME), watts: QtyAny::new(4.0, dims::POWER) },
            ],
        }],
        enthalpy: None,
        capacities: setup.regions.iter().map(|row| TransientRegionCapacity {
            region: row.region.clone(),
            volumetric_heat_capacity: QtyAny::new(10.0, dims::VOLUMETRIC_HEAT_CAPACITY),
            source: "synthetic capacity for numerical contact verification".into(),
        }).collect(),
    });
    f
}

fn receipt(f: &Fixture, ordinal: usize) -> (J, J) {
    receipt_with_pack(f, ordinal, temperature_pack())
}

fn receipt_with_pack(f: &Fixture, ordinal: usize, pack: (Vec<u8>, String, String)) -> (J, J) {
    let (source, contact, ledger) = prepare_with_pack(f, ordinal, None, pack);
    let result = run(f, &source, &contact, &ledger);
    assert_eq!(result.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", result.stdout, result.stderr);
    let result = J::parse(&result.stdout).unwrap();
    let completed = artifact(&ledger, result.str_field("run_receipt").unwrap());
    let stage = completed.get("stages").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("stage") == Some("conduction")).unwrap();
    let conduction = artifact(&ledger, stage.str_field("receipt").unwrap());
    let field = artifact(&ledger, conduction.str_field("solution_artifact").unwrap());
    (conduction, field)
}

fn n(row: &J, key: &str) -> f64 { row.f64_field(key).unwrap() }

fn maximum(field: &J) -> f64 {
    field.get("temperature").unwrap().as_array().unwrap().iter()
        .map(|value| value.as_f64().unwrap()).fold(f64::NEG_INFINITY, f64::max)
}

#[test]
fn native_transient_contact_requeries_endpoints_and_preserves_physical_history() {
    for nonlinear in [false, true] {
        let f = fixture(nonlinear, 0.5);
        let (conduction, field) = receipt(&f, 50);
        let time = conduction.get("transient").unwrap();
        let contact = time.get("temperature_contact").unwrap();
        let numerical = time.get("nonlinear").unwrap();
        assert_eq!(contact.str_field("method"), Some("source-backed-contact-newton"));
        assert_eq!(n(contact, "max_thermal_responses_per_step"), 1.0);
        assert_eq!(numerical.str_field("method"), Some("newton-fgmres"));
        assert_eq!(numerical.str_field("conductivity"), Some(if nonlinear {
            "endpoint-k(T)-with-k-prime-tangent"
        } else { "temperature-independent" }));
        assert_eq!(numerical.str_field("temperature_contact"),
            Some("endpoint-source-law-with-generalized-primal-tangent"));
        let (_, card, selected) = temperature_pack();
        let mut total_responses = 0.0;
        let mut total_updates = 0.0;
        let mut total_backtracks = 0.0;
        for grid in ["coarse", "fine"] {
            let rows = time.get(grid).unwrap().as_array().unwrap();
            let mut delivered_j = 0.0;
            let mut grid_updates = 0.0;
            let mut min_query = f64::INFINITY;
            let mut max_query = f64::NEG_INFINITY;
            assert!(rows.iter().any(|row| n(row, "time_s") == 0.75), "each grid must land on the authored power switch");
            for row in rows {
                let physical = row.get("temperature_contact").unwrap();
                assert_eq!(physical.str_field("method"), Some("source-backed-contact-newton"));
                assert_eq!(n(physical, "thermal_responses"), 1.0);
                assert!(n(physical, "physical_residual_norm_j") <= n(physical, "physical_residual_tolerance_j"));
                assert!(n(physical, "physical_energy_residual_j").abs() <= n(time, "energy_tolerance_j"));
                assert_eq!(n(row, "energy_residual_j"), n(physical, "physical_energy_residual_j"));
                assert_eq!(n(row, "krylov_iterations"), n(physical, "krylov_iterations"));
                assert!(n(physical, "krylov_iterations") <= n(contact, "max_krylov_iterations_per_step"));
                assert!(n(physical, "nonlinear_updates") <= n(numerical, "max_updates_per_step"));
                assert!(n(physical, "nonlinear_backtracks") <= n(physical, "nonlinear_updates")
                    * n(numerical, "max_backtracks_per_update"));
                total_responses += n(physical, "thermal_responses");
                grid_updates += n(physical, "nonlinear_updates");
                total_backtracks += n(physical, "nonlinear_backtracks");
                delivered_j += n(row, "dt_s") * n(row, "source_w");
                for state in physical.get("face_states").unwrap().as_array().unwrap() {
                    let query = n(state, "query_temperature_k");
                    assert_eq!(query, 0.5 * n(state, "side_a_temperature_k") + 0.5 * n(state, "side_b_temperature_k"));
                    assert!((n(state, "resistance_m2_k_w") - 0.001 * query).abs() < 1e-14);
                    assert_eq!(state.str_field("interface_card"), Some(card.as_str()));
                    assert_eq!(state.str_field("selected_claim"), Some(selected.as_str()));
                    assert_eq!(state.str_field("usage_receipt").unwrap().len(), 64);
                    min_query = min_query.min(query);
                    max_query = max_query.max(query);
                }
            }
            assert!((delivered_j - 14.0).abs() < 1e-10, "contact trials must not replay workload energy: {delivered_j}");
            assert!(max_query - min_query > 1e-4, "the source must follow the evolving interface");
            assert_eq!(grid_updates, n(numerical, &format!("{grid}_updates")));
            total_updates += grid_updates;
        }
        assert_eq!(n(contact, "total_thermal_responses"), total_responses);
        assert_eq!(total_responses, n(time, "total_steps"));
        assert_eq!(n(contact, "coarse_thermal_responses"), n(time, "coarse_steps"));
        assert_eq!(n(contact, "fine_thermal_responses"), n(time, "fine_steps"));
        assert!(total_updates > 0.0);
        assert_eq!(total_updates, n(numerical, "total_updates"));
        assert_eq!(total_updates, n(conduction.get("solver").unwrap(), "iterations"));
        assert_eq!(total_backtracks, n(numerical, "total_backtracks"));
        let energy = time.get("energy").unwrap();
        assert!(n(energy, "window_residual_j").abs() <= n(time, "fine_steps") * n(time, "energy_tolerance_j"));
        assert!(n(energy, "maximum_step_residual_j") <= n(time, "energy_tolerance_j"));
        assert_eq!(maximum(&field), n(time.get("temporal_error").unwrap(), "fine_final_k"));
        let fluxes = conduction.get("interfaces").unwrap().get("fluxes").unwrap().as_array().unwrap();
        assert!(fluxes.iter().all(|row| n(row, "heat_rate_a_to_b_w") < 0.0), "heat must flow from hot B to cold A");
        let (_, frozen, _) = f.solve(0.29315, None, 51);
        assert!((maximum(&field) - maximum(&frozen)).abs() > 1e-6,
            "a live source law must not reduce to the initial frozen resistance");
    }
}

#[test]
fn native_transient_contact_newton_starts_at_the_original_resistance_knot() {
    for nonlinear in [false, true] {
        let mut f = fixture(nonlinear, 0.5);
        let setup = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
        // 300/3 is exact: both original face means start at the source corner.
        setup.transient.as_mut().unwrap().initial_temperature =
            TransientInitialTemperature::Uniform(QtyAny::new(300.0, dims::TEMPERATURE));
        setup.boundaries[0].condition = ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(300.0, dims::TEMPERATURE),
        };
        let corner = temperature_pack_from_knots(vec![(200.0, 0.2), (300.0, 0.3), (450.0, 0.6)]);
        let (conduction, field) = receipt_with_pack(&f, 55, corner);
        let time = conduction.get("transient").unwrap();
        assert_eq!(n(time, "initial_temperature_k"), 300.0);
        assert_eq!(time.get("temperature_contact").unwrap().str_field("method"),
            Some("source-backed-contact-newton"));
        assert!(n(time.get("nonlinear").unwrap(), "total_updates") > 0.0);
        let mut moved_above_corner = false;
        for grid in ["coarse", "fine"] {
            for row in time.get(grid).unwrap().as_array().unwrap() {
                let physical = row.get("temperature_contact").unwrap();
                assert_eq!(n(physical, "thermal_responses"), 1.0);
                assert!(n(physical, "physical_residual_norm_j") <= n(physical, "physical_residual_tolerance_j"));
                assert!(n(physical, "physical_energy_residual_j").abs() <= n(time, "energy_tolerance_j"));
                for state in physical.get("face_states").unwrap().as_array().unwrap() {
                    let query = n(state, "query_temperature_k");
                    let slope = if query <= 300.0 { 0.001 } else { 0.002 };
                    let expected = 0.3 + slope * (query - 300.0);
                    assert!((n(state, "resistance_m2_k_w") - expected).abs() < 1e-13,
                        "a primal knot tangent must preserve the original resistance law");
                    moved_above_corner |= query > 300.0 + 1e-8;
                }
            }
        }
        assert!(moved_above_corner);
        assert!(maximum(&field) > 300.0);
    }
}

#[test]
fn native_transient_contact_converges_on_the_existing_nested_time_grids() {
    let mut values = Vec::new();
    for max_step in [0.5, 0.25, 0.125] {
        let mut f = fixture(false, max_step);
        f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .transient.as_mut().unwrap().power_schedules.clear();
        let (conduction, field) = receipt(&f, 60);
        let time = conduction.get("transient").unwrap();
        assert_eq!(time.get("temporal_error").unwrap().str_field("status"), Some("estimated-order-one"));
        values.push(maximum(&field));
    }
    let coarse = (values[1] - values[0]).abs();
    let fine = (values[2] - values[1]).abs();
    assert!(coarse > 1e-7 && fine > 0.0 && fine < 0.8 * coarse,
        "the same physical contact trajectory must converge as dt shrinks: {values:?}");
}

#[test]
fn native_transient_contact_refuses_late_source_exit_without_a_completed_field() {
    let mut f = fixture(false, 0.25);
    let envelope = f.project.envelope.as_mut().unwrap();
    envelope.ambient_hi.value = 293.25;
    let requirement = &mut f.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 293.25;
    requirement.margin.value = 0.0;
    let time = f.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap();
    time.horizon.value = 1.0;
    time.power_schedules[0].steps = vec![
        TransientPowerStep { until: QtyAny::new(0.5, dims::TIME), watts: QtyAny::new(0.0, dims::POWER) },
        TransientPowerStep { until: QtyAny::new(1.0, dims::TIME), watts: QtyAny::new(100.0, dims::POWER) },
    ];
    let limited = temperature_pack_from_knots(vec![(293.0, 0.293), (293.3, 0.2933)]);
    let (source, contact, ledger) = prepare_with_pack(&f, 70, None, limited);
    let result = run(&f, &source, &contact, &ledger);
    assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
    let text = format!("{}{}", result.stdout, result.stderr);
    assert!(text.contains("cli-solve-conduction-transient"), "{text}");
    assert!(text.contains("contact"), "{text}");
    assert!(!text.contains("run_receipt"), "an unsupported source endpoint cannot publish a completed run: {text}");
    // The original two zero-power intervals were admissible; widening only
    // the source support lets the unchanged driver complete from its prefix.
    let (completed, _) = receipt(&f, 71);
    assert_eq!(completed.get("transient").unwrap().str_field("status"), Some("completed"));
}

#[test]
fn native_transient_contact_keeps_material_derivatives_explicit() {
    let f = fixture(false, 0.5);
    for (ordinal, output) in [
        "temperature-final-material-adjoint",
        "temperature-sampled-peak-material-adjoint",
    ]
    .into_iter()
    .enumerate()
    {
        let (source, contact, ledger) = prepare(&f, 80 + ordinal, Some(output));
        let result = run(&f, &source, &contact, &ledger);
        assert_ne!(result.exit_code, fs_cli::exit::SUCCESS);
        assert!(
            format!("{}{}", result.stdout, result.stderr)
                .contains("cli-solve-conduction-contact-composition")
        );
    }
}
