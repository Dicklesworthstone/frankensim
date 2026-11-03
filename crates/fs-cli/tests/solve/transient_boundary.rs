//! Real imported .fsim thermal cycles, with independent storage/flux oracles.
use super::*;
use fs_project::{
    ConductionTransient, TransientBoundarySchedule, TransientBoundaryStep,
    TransientBoundaryTarget, TransientInitialTemperature, TransientPowerStep,
    TransientRegionPower, spec::dims,
};
#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

fn time(spec: &mut ProjectSpec) -> &mut ConductionTransient {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap()
}

fn project(enthalpy: bool) -> ProjectSpec {
    let mut spec = if enthalpy {
        // rho=6 kg/m³, Cp=10 J/(kg K), V=1/6 m³: total capacity 10 J/K.
        transient_enthalpy::project(6.0, 1.0, &[(-400.0, 253.15, 0.0), (1000.0, 393.15, 0.0)])
    } else {
        let mut spec = conduction_fixture_project(7, &tetra_stl());
        transient_product::declaration(&mut spec, 60.0, 1.0);
        spec
    };
    time(&mut spec).max_step.value = 0.6;
    for source in spec.power.as_mut().unwrap() { source.watts.value = 0.0; }
    spec.envelope.as_mut().unwrap().ambient_lo.value = 280.0;
    spec.envelope.as_mut().unwrap().ambient_hi.value = 315.0;
    spec
}

fn schedule(spec: &mut ProjectSpec, kind: TransientBoundaryTarget, steps: &[(f64, f64)]) {
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient.as_mut().unwrap().boundary_schedules = vec![TransientBoundarySchedule {
        target: setup.boundaries[0].target.clone(), kind,
        source: "manufactured thermal cycle for an independent physical oracle".into(),
        steps: steps.iter().map(|&(until, value)| TransientBoundaryStep {
            until: QtyAny::new(until, dims::TIME), value: QtyAny::new(value, kind.dims()),
        }).collect(),
    }];
}

fn run(spec: &ProjectSpec) -> (J, J) {
    let decoded = decode(spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, spec, tetra_stl());
    let (_, _, conduction, field) =
        run_conjugate_to_completion(&ledger, &decoded, &fixture_cards());
    (J::parse(&conduction).unwrap(), J::parse(&field).unwrap())
}

fn n(row: &J, key: &str) -> f64 { row.f64_field(key).unwrap() }
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}");
}

#[test]
fn g1_boundary_cycle_rebuilds_dirichlet_lifts_and_retains_original_storage_and_peak() {
    for initial in [293.15, 320.0] {
        let mut spec = project(true);
        time(&mut spec).initial_temperature =
            TransientInitialTemperature::Uniform(QtyAny::new(initial, dims::TEMPERATURE));
        if let Some(storage) = &mut time(&mut spec).enthalpy {
            storage.materials[0].initial_specific_enthalpy.value = 10.0 * (initial - 293.15);
        }
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::FixedTemperature {
                temperature: QtyAny::new(280.0, dims::TEMPERATURE),
            };
        schedule(&mut spec, TransientBoundaryTarget::FixedTemperature,
            &[(0.3, 290.0), (0.7, 310.0), (1.0, 300.0)]);
        let (conduction, field) = run(&spec);
        let history = conduction.get("transient").unwrap();
        close(n(history, "initial_temperature_k"), initial, 0.0);
        for grid in ["coarse", "fine"] {
            let rows = history.get(grid).unwrap().as_array().unwrap();
            assert!(rows.iter().any(|row| n(row, "time_s") == 0.3));
            assert!(rows.iter().any(|row| n(row, "time_s") == 0.7));
            let mut previous = initial;
            for row in rows {
                let end = n(row, "time_s");
                let expected = if end <= 0.3 { 290.0 } else if end <= 0.7 { 310.0 } else { 300.0 };
                close(n(row, "final_region_max_k"), expected, 2e-9);
                close(n(row, "stored_energy_change_j"), 10.0 * (expected - previous), 2e-7);
                close(n(row, "net_input_w") * n(row, "dt_s"), 10.0 * (expected - previous), 2e-7);
                previous = expected;
            }
        }
        let peak = history.get("sampled_peak").unwrap();
        close(n(peak, "value_k"), initial.max(310.0), 2e-9);
        close(n(peak, "time_s"), if initial > 310.0 { 0.0 } else { 0.5 }, 2e-15);
        let energy = history.get("energy").unwrap();
        close(n(energy, "stored_change_j"), 10.0 * (300.0 - initial), 2e-7);
        assert!(n(energy, "maximum_step_residual_j") <= 1e-6);
        for value in field.get("temperature").unwrap().as_array().unwrap() {
            close(value.as_f64().unwrap(), 300.0, 2e-9);
        }
        assert_eq!(history.get("boundary_history").unwrap().str_field("sampling"),
            Some("held-interval-start"));
    }
}

#[test]
fn g1_signed_boundary_heat_pulses_deliver_the_original_surface_area_energy() {
    // Unit right tetrahedron: three area-1/2 faces and one sqrt(3)/2 face.
    let area = (3.0 + 3.0_f64.sqrt()) / 2.0;
    for enthalpy in [false, true] {
        let mut spec = project(enthalpy);
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::HeatFlux { outward_flux: QtyAny::new(99.0, dims::HEAT_FLUX) };
        schedule(&mut spec, TransientBoundaryTarget::HeatFlux,
            &[(0.25, 2.0), (0.75, -1.0), (1.0, 3.0)]);
        let (conduction, _) = run(&spec);
        let history = conduction.get("transient").unwrap();
        // Boundary-only switches must reach the fine grid even without a
        // regional source schedule. A uniform six-step grid misses both.
        for (grid, expected_ends) in [
            ("coarse", &[0.25, 0.75, 1.0][..]),
            ("fine", &[0.125, 0.25, 0.5, 0.75, 0.875, 1.0][..]),
        ] {
            let rows = history.get(grid).unwrap().as_array().unwrap();
            assert_eq!(rows.len(), expected_ends.len(), "{grid}");
            for (row, &expected_end) in rows.iter().zip(expected_ends) {
                let end = n(row, "time_s");
                close(end, expected_end, 0.0);
                let flux = if end <= 0.25 { 2.0 } else if end <= 0.75 { -1.0 } else { 3.0 };
                close(n(row, "net_input_w"), -area * flux, 2e-8);
                close(n(row, "stored_energy_change_j"), -area * flux * n(row, "dt_s"), 2e-7);
            }
        }
        let energy = history.get("energy").unwrap();
        close(n(energy, "stored_change_j"), -0.75 * area, 2e-6);
        close(n(energy, "integrated_net_input_j"), -0.75 * area, 2e-7);
        assert!(n(energy, "maximum_step_residual_j") <= 1e-6);
    }
}

#[test]
fn g1_source_and_boundary_pulses_share_one_nested_physical_time_grid() {
    let area = (3.0 + 3.0_f64.sqrt()) / 2.0;
    for enthalpy in [false, true] {
        let mut spec = project(enthalpy);
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::HeatFlux { outward_flux: QtyAny::new(99.0, dims::HEAT_FLUX) };
        schedule(&mut spec, TransientBoundaryTarget::HeatFlux,
            &[(0.25, 2.0), (0.75, -1.0), (1.0, 3.0)]);
        time(&mut spec).power_schedules = vec![TransientRegionPower {
            region: "air".into(),
            source: "manufactured regional heat pulse with a coincident boundary switch".into(),
            steps: [(0.5, 4.0), (0.75, 0.0), (1.0, 2.0)].iter()
                .map(|&(until, watts)| TransientPowerStep {
                    until: QtyAny::new(until, dims::TIME),
                    watts: QtyAny::new(watts, dims::POWER),
                }).collect(),
        }];
        let (conduction, field) = run(&spec);
        let history = conduction.get("transient").unwrap();
        for (grid, expected_ends) in [
            ("coarse", &[0.25, 0.5, 0.75, 1.0][..]),
            ("fine", &[0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0][..]),
        ] {
            let rows = history.get(grid).unwrap().as_array().unwrap();
            assert_eq!(rows.len(), expected_ends.len(), "{grid}");
            for (row, &end) in rows.iter().zip(expected_ends) {
                close(n(row, "time_s"), end, 0.0);
                let watts = if end <= 0.5 { 4.0 } else if end <= 0.75 { 0.0 } else { 2.0 };
                let flux = if end <= 0.25 { 2.0 } else if end <= 0.75 { -1.0 } else { 3.0 };
                let net_w = watts - area * flux;
                close(n(row, "source_w"), watts, 2e-8);
                close(n(row, "net_input_w"), net_w, 2e-8);
                close(n(row, "stored_energy_change_j"), net_w * n(row, "dt_s"), 2e-7);
            }
        }
        // Source energy is 4 W * 0.5 s + 2 W * 0.25 s; prescribed outward
        // boundary energy is A * (2 * 0.25 - 1 * 0.5 + 3 * 0.25) J/m².
        let expected_j = 2.5 - 0.75 * area;
        let energy = history.get("energy").unwrap();
        close(n(energy, "stored_change_j"), expected_j, 2e-6);
        close(n(energy, "integrated_net_input_j"), expected_j, 2e-7);
        assert!(n(energy, "maximum_step_residual_j") <= 1e-6);
        let temperatures = field.get("temperature").unwrap().as_array().unwrap();
        // Both storage declarations give this tetrahedron 10 J/K total,
        // distributed equally over its four lumped nodal masses.
        assert_eq!(temperatures.len(), 4);
        let mean_k = temperatures.iter().map(|value| value.as_f64().unwrap()).sum::<f64>()
            / temperatures.len() as f64;
        close(10.0 * (mean_k - 293.15), expected_j, 2e-6);
    }
}

#[test]
fn g1_reference_history_retains_the_original_robin_coefficient_in_both_storage_models() {
    for enthalpy in [false, true] {
        let mut spec = project(enthalpy);
        time(&mut spec).max_step.value = 0.25;
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::Convection {
                coefficient: QtyAny::new(2.0, dims::HEAT_TRANSFER_COEFFICIENT),
                reference_temperature: QtyAny::new(285.0, dims::TEMPERATURE),
            };
        let (_, original) = run(&spec);
        schedule(&mut spec, TransientBoundaryTarget::ConvectionTemperature,
            &[(0.25, 285.0), (0.75, 285.0), (1.0, 285.0)]);
        let (_, held) = run(&spec);
        assert_eq!(original.get("temperature"), held.get("temperature"),
            "constant history must preserve the same native Robin operator and grid");
        time(&mut spec).boundary_schedules[0].steps[1].value.value = 310.0;
        let (conduction, heated) = run(&spec);
        let original = held.get("temperature").unwrap().as_array().unwrap();
        let heated = heated.get("temperature").unwrap().as_array().unwrap();
        assert!(heated.iter().zip(original).all(|(hot, cold)|
            hot.as_f64().unwrap() > cold.as_f64().unwrap() + 0.01));
        assert!(n(conduction.get("transient").unwrap().get("energy").unwrap(),
            "maximum_step_residual_j") <= 1e-6);
    }
}

#[test]
fn g1_scheduled_capacity_adjoints_retain_the_original_forward_field() {
    for output in ["temperature-final-adjoint", "temperature-sampled-peak-adjoint"] {
        let mut spec = project(false);
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
            ThermalBoundaryCondition::HeatFlux { outward_flux: QtyAny::new(0.0, dims::HEAT_FLUX) };
        schedule(&mut spec, TransientBoundaryTarget::HeatFlux, &[(0.3, 2.0), (1.0, -1.0)]);
        let (original, original_field) = run(&spec);
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: output.into(), kind: "report".into(), region: None,
        });
        let (differentiated, field) = run(&spec);
        assert_eq!(field.get("temperature"), original_field.get("temperature"),
            "a history-aware reverse must not alter the physical trajectory");
        let ordinary = original.get("transient").unwrap();
        let history = differentiated.get("transient").unwrap();
        for key in ["coarse", "fine", "energy", "sampled_peak", "boundary_history"] {
            assert_eq!(history.get(key), ordinary.get(key), "{key}");
        }
        let adjoint = history.get("adjoint").unwrap();
        assert_eq!(adjoint.str_field("output"), Some(output));
        assert!(n(adjoint, "maximum_dual_relative_residual") <= 1e-8);
        assert!(adjoint.get("prescribed_heat_flux").is_none(),
            "the unused static flux is not an independent schedule coordinate");
    }
}

#[test]
fn g1_capacity_boundary_cycle_matches_one_free_node_and_full_prescribed_reaction() {
    let bytes = tetra_stl();
    let mut spec = transient_product::driven_face_project(&bytes);
    time(&mut spec).horizon.value = 1.0;
    time(&mut spec).max_step.value = 0.6;
    time(&mut spec).initial_temperature =
        TransientInitialTemperature::Uniform(QtyAny::new(320.0, dims::TEMPERATURE));
    schedule(&mut spec, TransientBoundaryTarget::FixedTemperature,
        &[(0.3, 310.0), (0.7, 300.0), (1.0, 305.0)]);
    let cards = transient_product::transient_conductivity_cards("constant-scalar");
    rebind_to(&mut spec, &cards);
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, bytes);
    let (_, _, conduction, field) =
        run_conjugate_to_completion(&ledger, &decode(&spec), &cards);
    let conduction = J::parse(&conduction).unwrap();
    let history = conduction.get("transient").unwrap();
    close(n(conduction.get("solver").unwrap(), "free_dofs"), 1.0, 0.0);
    let mut fine_final = 320.0;
    for grid in ["coarse", "fine"] {
        let mut free = 320.0;
        let mut bath = 320.0;
        for row in history.get(grid).unwrap().as_array().unwrap() {
            let end = n(row, "time_s");
            let next_bath = if end <= 0.3 { 310.0 } else if end <= 0.7 { 300.0 } else { 305.0 };
            let dt = n(row, "dt_s");
            // The original unit tetra has lumped row-sum C=15/24 J/K.
            // k=2 W/(m K) and the opposite face fixes three vertices, leaving
            // C*(T-new-T-old)+dt*(k/2)*(T-new-bath)=0 at the fourth.
            let next_free = (0.625 * free + dt * next_bath) / (0.625 + dt);
            let stored = 0.625 * (next_free - free + 3.0 * (next_bath - bath));
            close(n(row, "final_region_max_k"), next_free.max(next_bath), 2e-7);
            close(n(row, "stored_energy_change_j"), stored, 2e-7);
            close(dt * n(row, "net_input_w"), stored, 2e-7);
            free = next_free;
            bath = next_bath;
        }
        if grid == "fine" { fine_final = free; }
    }
    close(n(history.get("sampled_peak").unwrap(), "value_k"), 320.0, 0.0);
    close(n(history.get("sampled_peak").unwrap(), "time_s"), 0.0, 0.0);
    let field = J::parse(&field).unwrap();
    let mut values: Vec<_> = field.get("temperature").unwrap().as_array().unwrap()
        .iter().map(|value| value.as_f64().unwrap()).collect();
    values.sort_by(f64::total_cmp);
    for &fixed in &values[..3] { close(fixed, 305.0, 0.0); }
    close(values[3], fine_final, 2e-7);
    close(n(history.get("energy").unwrap(), "stored_change_j"),
        0.625 * (values.iter().sum::<f64>() - 4.0 * 320.0), 2e-7);
}
