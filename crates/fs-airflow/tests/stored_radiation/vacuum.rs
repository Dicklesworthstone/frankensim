//! Radiation-only traces use the actual FEM/secant owners, not a tiny air HTC.
use super::*;
use fs_conduction::{ConductionProblem, InitialGuess, SolveConfig};

fn vacuum_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("exposed", |f| f.outward_normal[2] < -0.99, ThermalBc::adiabatic()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}

fn mixed_boundary(mesh: &ConductionMesh, rows: &[CorrelationRobinBoundary])
    -> std::result::Result<ThermalBoundary, ConductionError>
{
    ThermalBoundaryBuilder::new(mesh)
        .region("exposed", |f| f.outward_normal[2] < -0.99, ThermalBc::adiabatic())?
        .region("body", |f| f.outward_normal[2] > 0.99, rows[0].boundary_condition().clone())?
        .adiabatic_remainder().finish()
}

#[test]
fn pure_radiation_heats_cools_and_preserves_equilibrium_without_air() {
    for ambient in [250.0, 300.0, 400.0] {
        fixture(|cx, _, _, _, state, mesh, material| {
            let boundary = vacuum_boundary(mesh);
            let original = boundary.clone();
            let heat = ScalarField::uniform("no volumetric heat", 0.0).unwrap();
            let old = state.solid_specific_enthalpy_j_kg().to_vec();
            let result = state.storage().advance_with_ambient_radiation(cx,
                ConductionProblem { mesh, boundary: &boundary, material, source: &heat, element_materials: None },
                None, &old, 1.0, controls().thermal.solid,
                &[patch(0.9, ambient, (250.0,450.0))], controls().radiation).unwrap();
            assert!(!boundary.has_robin());
            assert!(result.combined_boundary.has_robin());
            assert_eq!(boundary, original);
            assert_eq!(old.as_slice(), state.solid_specific_enthalpy_j_kg());
            assert!(result.convective_robin_fluxes.is_empty());
            assert_eq!(result.convective_out_w, 0.0);
            assert_eq!(result.conduction.robin_fluxes.len(), 1);
            let radiation = result.radiation.nonlinear_radiation_out_w;
            assert!((result.conduction.stored_energy_change_j + radiation).abs() <= controls().thermal.solid.energy_tolerance_j);
            assert!(result.physical_residual_norm_j <= result.physical_residual_tolerance_j);
            if ambient < 300.0 {
                assert!(radiation > 0.0 && result.conduction.stored_energy_change_j < 0.0);
            } else if ambient > 300.0 {
                assert!(radiation < 0.0 && result.conduction.stored_energy_change_j > 0.0);
            } else {
                assert!(radiation.abs() < 1e-9);
                assert!(result.conduction.temperature.iter().all(|t| (t-300.0).abs() < 1e-8));
            }
        });
    }
}

#[test]
fn radiation_supplies_the_only_steady_thermal_anchor() {
    fixture(|cx, _, _, _, _, mesh, material| {
        let boundary = vacuum_boundary(mesh);
        let heat = ScalarField::uniform("no volumetric heat", 0.0).unwrap();
        let result = fs_conduction::radiation::solve_with_ambient_radiation(cx,
            ConductionProblem { mesh, boundary: &boundary, material, source: &heat, element_materials: None },
            None, &[patch(0.9,350.0,(250.0,450.0))],
            SolveConfig { initial: InitialGuess::Uniform(300.0), ..Default::default() },
            controls().radiation).unwrap();
        assert!(boundary.dirichlet().is_empty() && !boundary.has_robin());
        assert!(result.convective_robin_fluxes.is_empty());
        assert_eq!(result.convective_out_w, 0.0);
        assert!(result.conduction.temperature.iter().all(|t| (t-350.0).abs() < 1e-6));
        assert!(result.radiation.nonlinear_radiation_out_w.abs() < 1e-7);
    });
}

#[test]
fn vacuum_finish_changes_both_histories_without_creating_a_coolant_port() {
    let run = |emissivity| fixture(|cx, _, network, pump, state, mesh, material| {
        let heat = ScalarField::uniform("20 W body heat", 20.0/0.00006).unwrap();
        let result = network.advance_stored_solid_with_pump_and_radiation(cx, state,
            1.0, 0.0, &[300.0], pump,
            StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
            &[patch(emissivity,250.0,(250.0,450.0))], controls(),
            |_, rows| mixed_boundary(mesh, rows)).unwrap();
        let response = &result.step.liquid.circuit.state.solid.state;
        assert_eq!(response.external_convection_out_w, 0.0);
        assert_eq!(response.thermal.convective_robin_fluxes.len(), 1);
        assert_eq!(response.thermal.convective_robin_fluxes[0].region, "body");
        assert_eq!(response.coolant_robin_fluxes.len(), 1);
        assert_eq!(response.thermal.conduction.robin_fluxes.len(), 2);
        assert!((response.coolant_out_w - result.step.liquid.returned_heat_j).abs() < 2e-6);
        assert!(result.step.combined_energy_residual_j.abs() < 2e-6);
        assert!(result.step.interface_energy_residual_j.abs() < 2e-6);
        assert!((result.step.solid_external_heat_j - (20.0-response.thermal.radiation.nonlinear_radiation_out_w)).abs() < 1e-8);
        (state.reservoir().state().temperature_k(), response.thermal.conduction.temperature.clone())
    });
    let low = run(0.1);
    let high = run(0.9);
    assert!(high.0 < low.0);
    assert!(high.1.iter().zip(&low.1).all(|(a,b)| a < b));
}

#[test]
fn radiation_cannot_erase_prescribed_flux_or_temperature() {
    fixture(|cx, _, _, _, state, mesh, material| {
        let heat = ScalarField::uniform("body", 0.0).unwrap();
        let zero = vacuum_boundary(mesh);
        // The new admission is private to radiation: ordinary Robin invariants
        // and the public convection replacement contract are unchanged.
        assert!(ThermalBc::robin(0.0,300.0).is_err());
        assert!(zero.with_uniform_robin_replacements(&[(0,5.0,300.0)]).is_err());
        for condition in [ThermalBc::neumann(1.0).unwrap(), ThermalBc::dirichlet(300.0).unwrap(),
            ThermalBc::Neumann { outward_flux: ScalarField::Nodal(vec![0.0; mesh.vertex_count()]) }]
        {
            let boundary = ThermalBoundaryBuilder::new(mesh)
                .region("exposed", |f| f.outward_normal[2] < -0.99, condition).unwrap()
                .adiabatic_remainder().finish().unwrap();
            let saved = boundary.clone();
            assert!(state.storage().advance_with_ambient_radiation(cx,
                ConductionProblem { mesh, boundary: &boundary, material, source: &heat, element_materials: None },
                None, state.solid_specific_enthalpy_j_kg(), 1.0, controls().thermal.solid,
                &[patch(0.9,250.0,(250.0,450.0))], controls().radiation).is_err());
            assert_eq!(saved, boundary);
        }
    });
}

#[test]
fn vacuum_source_work_and_cancel_failures_preserve_both_histories() {
    for case in 0..3 {
        fixture(|cx, gate, network, pump, state, mesh, material| {
            let old = state.solid_specific_enthalpy_j_kg().to_vec();
            let tank = state.reservoir().clone();
            let heat = ScalarField::uniform("20 W body heat", 20.0/0.00006).unwrap();
            let span = if case == 0 { (300.0,300.0) } else { (250.0,450.0) };
            let mut limits = controls();
            if case == 1 { limits.max_total_enthalpy_solves = 1; }
            let result = network.advance_stored_solid_with_pump_and_radiation(cx, state,
                1.0, 0.0, &[300.0], pump,
                StoredSolidProblem { mesh, material, source: &heat, element_materials: None, interfaces: None },
                &[patch(0.9,250.0,span)], limits, |_, rows| {
                    if case == 2 { gate.request(); }
                    mixed_boundary(mesh, rows)
                });
            assert!(result.is_err());
            assert_eq!(state.reservoir(), &tank);
            assert_eq!(state.solid_specific_enthalpy_j_kg(), old.as_slice());
        });
    }
}

#[test]
fn existing_transient_adjoint_differentiates_the_radiation_only_endpoint() {
    fixture(|cx, _, _, _, state, mesh, material| {
        let boundary = vacuum_boundary(mesh);
        let heat = ScalarField::uniform("20 W body heat", 20.0/0.00006).unwrap();
        let problem = ConductionProblem { mesh, boundary: &boundary, material,
            source: &heat, element_materials: None };
        let old = state.solid_specific_enthalpy_j_kg();
        let storage = state.storage();
        let linearization = storage.linearize_step_with_ambient_radiation(cx, problem, None,
            old, 1.0, controls().thermal.solid, &[patch(0.5,250.0,(250.0,450.0))],
            controls().radiation, 2*mesh.vertex_count()).unwrap();
        let weights = vec![1.0/mesh.vertex_count() as f64; mesh.vertex_count()];
        let seed = linearization.temperature_pullback(cx, &weights).unwrap();
        let gradient = linearization.pullback(cx, &seed, fs_conduction::LinearConfig {
            tolerance: 1e-11, max_iterations: 64, restart: 24 }).unwrap();
        let value = |epsilon, ambient, dt| {
            let result = storage.advance_with_ambient_radiation(cx, problem, None, old,
                dt, controls().thermal.solid, &[patch(epsilon,ambient,(250.0,450.0))],
                controls().radiation).unwrap();
            result.conduction.temperature.iter().sum::<f64>()/mesh.vertex_count() as f64
        };
        for (analytic, difference) in [
            (gradient.emissivities[0], (value(0.5001,250.0,1.0)-value(0.4999,250.0,1.0))/0.0002),
            (gradient.reservoir_temperatures[0], (value(0.5,250.01,1.0)-value(0.5,249.99,1.0))/0.02),
            (gradient.time_step_s, (value(0.5,250.0,1.001)-value(0.5,250.0,0.999))/0.002),
        ] {
            assert!((analytic-difference).abs() <= 2e-6*analytic.abs().max(difference.abs()).max(1.0),
                "analytic={analytic}, finite difference={difference}");
        }
        assert!(gradient.emissivities[0] < 0.0);
        assert!(gradient.reservoir_temperatures[0] > 0.0);
    });
}
