//! Finite solid walls exchange the original measured boundary heat with
//! both native thermal lattices; they are not prescribed-temperature baths.

use fs_lbm::ThermalLbm;
use fs_lbm::d3q19::ThermalLbm3;
use fs_lbm::thermal::boundary::{FiniteThermalWall, ThermalWall, ThermalWallField};
use fs_lbm::thermal::sources::{ThermalSourceError, ThermalSourceLimits};

fn limits() -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: 256,
        temperature_range: [-1.0, 2.0],
        balance_tolerance: 2e-11,
    }
}

fn solid(temperature: f64, capacity: f64, resistance: f64, heat: f64) -> ThermalWall {
    ThermalWall::FiniteSolid(FiniteThermalWall {
        temperature,
        heat_capacity: capacity,
        resistance,
        heat_per_step: heat,
    })
}

fn solid_energy(fields: [ThermalWallField<'_>; 2]) -> f64 {
    fields
        .into_iter()
        .flat_map(ThermalWallField::laws)
        .map(|wall| match wall {
            ThermalWall::FiniteSolid(solid) => solid.heat_capacity * solid.temperature,
            _ => 0.0,
        })
        .sum()
}

fn temperatures(field: ThermalWallField<'_>) -> Vec<f64> {
    field
        .laws()
        .filter_map(|wall| match wall {
            ThermalWall::FiniteSolid(solid) => Some(solid.temperature),
            _ => None,
        })
        .collect()
}

// Exercise identical public operations on both actual native solvers.
macro_rules! for_slabs {
    ($state:ident, $plane:ident, $forces:ident, $body:block) => {{
        {
            let mut $state = ThermalLbm::slab(3, 4, 0.9, 0.9, 0.0);
            let $plane = 3;
            let $forces = vec![[0.0; 2]; $state.gpop.len()];
            $body
        }
        {
            let mut $state = ThermalLbm3::slab(3, 2, 4, 0.9, 0.9, 0.0);
            let $plane = 6;
            let $forces = vec![[0.0; 3]; $state.gpop.len()];
            $body
        }
    }};
}

#[test]
fn closed_slabs_conserve_total_energy_and_reach_capacity_weighted_equilibrium() {
    for_slabs!(state, plane, forces, {
        let weights = state.gpop[0];
        state.gpop.fill(weights.map(|w| 0.25 * w));
        state
            .set_thermal_walls([solid(1.0, 2.0, 0.5, 0.0), solid(0.5, 5.0, 1.0, 0.0)])
            .unwrap();
        assert!(
            state
                .thermal_wall_fields()
                .into_iter()
                .all(|field| matches!(field, ThermalWallField::Spatial(_)))
        );
        let n = state.gpop.len();
        let zero_heat = vec![0.0; n];
        let initial = state.gpop[plane..n - plane].iter().flatten().sum::<f64>()
            + solid_energy(state.thermal_wall_fields());
        for _ in 0..2000 {
            let fluid = state
                .step_with_sources(&zero_heat, &forces, limits())
                .unwrap();
            for (heat, wall) in [fluid.bottom_wall, fluid.top_wall]
                .into_iter()
                .zip(state.last_wall_exchange().unwrap())
            {
                let solid = wall.solid.unwrap();
                assert_eq!(solid.heat_in, heat);
                assert_eq!(solid.source, 0.0);
                assert!(solid.defect().abs() < 2e-12);
            }
        }
        let final_energy = state.gpop[plane..n - plane].iter().flatten().sum::<f64>()
            + solid_energy(state.thermal_wall_fields());
        assert!((final_energy - initial).abs() < 2e-10 * plane as f64);
        // Per normal column: (4*.25 + 2*1 + 5*.5)/(4+2+5) = .5.
        for populations in &state.gpop[plane..n - plane] {
            assert!((populations.iter().sum::<f64>() - 0.5).abs() < 2e-8);
        }
        for temperature in state
            .thermal_wall_fields()
            .into_iter()
            .flat_map(temperatures)
        {
            assert!((temperature - 0.5).abs() < 2e-8);
        }
    });
}

#[test]
fn wall_source_is_deposited_once_and_partitioned_between_solid_and_fluid() {
    for_slabs!(state, plane, forces, {
        let weights = state.gpop[0];
        state.gpop.fill(weights.map(|w| 0.25 * w));
        state
            .set_thermal_walls([solid(0.25, 2.0, 3.0, 0.002), ThermalWall::Adiabatic])
            .unwrap();
        let n = state.gpop.len();
        let zero_heat = vec![0.0; n];
        for _ in 0..30 {
            let before = state.gpop[plane..n - plane].iter().flatten().sum::<f64>()
                + solid_energy(state.thermal_wall_fields());
            let fluid = state
                .step_with_sources(&zero_heat, &forces, limits())
                .unwrap();
            let report = state.last_wall_exchange().unwrap();
            let solid = report[0].solid.unwrap();
            assert_eq!(fluid.source, 0.0);
            assert_eq!(fluid.top_wall, 0.0);
            assert_eq!(fluid.bottom_wall, solid.heat_in);
            assert!(solid.stored_change > 0.0);
            assert!(solid.heat_in > 0.0 && solid.heat_in < solid.source);
            assert!((solid.source - 0.002 * plane as f64).abs() < 1e-16);
            let after = state.gpop[plane..n - plane].iter().flatten().sum::<f64>()
                + solid_energy(state.thermal_wall_fields());
            assert!((after - before - solid.source).abs() < 2e-13);
        }
    });
}

#[test]
fn uniform_declarations_evolve_independent_surface_cells() {
    for_slabs!(state, plane, forces, {
        let weights = state.gpop[0];
        for (cell, populations) in state.gpop.iter_mut().enumerate() {
            let temperature = 0.1 + 0.1 * (cell % plane) as f64;
            *populations = weights.map(|w| temperature * w);
        }
        state
            .set_thermal_walls([solid(1.0, 2.0, 0.5, 0.0), ThermalWall::Adiabatic])
            .unwrap();
        let zero_heat = vec![0.0; state.gpop.len()];
        state
            .step_with_sources(&zero_heat, &forces, limits())
            .unwrap();
        let wall = temperatures(state.thermal_wall_fields()[0]);
        assert_eq!(wall.len(), plane);
        assert!(wall.windows(2).all(|pair| pair[0] < pair[1]));
    });
}

#[test]
fn a_hot_solid_endpoint_refuses_without_publishing_fluid_or_wall_state() {
    for_slabs!(state, plane, forces, {
        let weights = state.gpop[0];
        state.gpop.fill(weights.map(|w| 0.25 * w));
        // The large wall resistance keeps the fluid-facing surface in range,
        // while this source pushes the finite bulk beyond the admitted 2.0.
        state
            .set_thermal_walls([solid(0.25, 0.1, 1000.0, 1.0), ThermalWall::Adiabatic])
            .unwrap();
        let zero_heat = vec![0.0; state.gpop.len()];
        let before = state.gpop.clone();
        let flow = state.grid.f.clone();
        let applied_forces = state.grid.fext.clone();
        let wall = state
            .thermal_wall_fields()
            .map(|field| field.laws().collect::<Vec<_>>());
        let report = state.last_wall_exchange();
        assert!(matches!(
            state.step_with_sources(&zero_heat, &forces, limits()),
            Err(ThermalSourceError::WallTemperature { wall: 0, .. }),
        ));
        assert_eq!(state.gpop, before);
        assert_eq!(state.grid.f, flow);
        assert_eq!(state.grid.fext, applied_forces);
        assert_eq!(state.last_wall_exchange(), report);
        assert_eq!(
            state
                .thermal_wall_fields()
                .map(|field| field.laws().collect::<Vec<_>>()),
            wall
        );
        assert_eq!(temperatures(state.thermal_wall_fields()[0]).len(), plane);
    });
}
