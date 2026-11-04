//! G3: original modal wall and mouth histories share one exact physical clock.
use super::*;
use fs_duct::modal::time_domain::{ModalViscothermalBudget, ModalViscothermalResolution};

#[test]
fn wall_and_loaded_mouth_keep_all_states_heat_and_clock_subdivision() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let radius = 0.002;
    let profile = RadiusProfile::exponential(radius, 0.0, 0.006).unwrap();
    let omega = 100.0_f64.powi(2) * gas.dynamic_viscosity / (gas.density * radius * radius);
    let resolution = ModalViscothermalResolution {
        poles: 12,
        band_rad_s: [0.99 * omega, 1.01 * omega],
        maximum_diffusive_relative_error: 0.02,
        maximum_boundary_layer_ratio: 0.1,
    };
    let wall_budget = ModalViscothermalBudget {
        maximum_states: 124,
        maximum_payload_bytes: 32 * 1024 * 1024,
        maximum_preparation_work: 200_000_000,
    };
    let duration = 1.0 / 48_000.0;
    for loaded in [false, true] {
        let mut results = Vec::new();
        for subdivision in [1, 2, 4] {
            let mut line = if loaded {
                ContinuousModalLine::new_viscothermal_loaded(
                    &profile,
                    &gas,
                    ModeCount(2),
                    1,
                    coupled_load(radius, gas, &[1.0, 0.6]),
                    resolution,
                    duration / subdivision as f64,
                    wall_budget,
                    budget(),
                    || false,
                )
            } else {
                ContinuousModalLine::new_viscothermal(
                    &profile,
                    &gas,
                    ModeCount(2),
                    1,
                    ModalEnd::Ports,
                    resolution,
                    duration / subdivision as f64,
                    wall_budget,
                    budget(),
                    || false,
                )
            }
            .unwrap();
            assert_eq!(line.duct_state().len(), 44);
            assert_eq!(line.load_state().len(), if loaded { 2 } else { 0 });
            assert_eq!(line.source().viscothermal_state_range().unwrap(), 8..44);
            let input = if loaded {
                vec![75.0, -20.0]
            } else {
                vec![75.0, -20.0, 0.0, 0.0]
            };
            let (mut incident, mut outgoing, mut heat) = (0.0, 0.0, 0.0);
            for _ in 0..subdivision {
                let step = line.push(&input).unwrap();
                incident += step.incident_energy;
                outgoing += step.outgoing_energy;
                heat += step.physical_heat;
            }
            let wall = line
                .source()
                .viscothermal_heat_w(line.duct_state())
                .unwrap();
            assert!(wall.viscous_w > 0.0 && wall.thermal_w > 0.0);
            assert!(heat > 0.0);
            assert!(line.duct_state()[8..].iter().any(|v| v.abs() > 1e-16));
            if loaded {
                assert!(line.load_state().iter().any(|v| v.abs() > 1e-16));
            }
            close(
                line.stored_energy_j() + heat + outgoing,
                incident,
                3e-8,
                1e-15,
            );
            results.push((
                line.state().to_vec(),
                heat,
                outgoing,
                line.stored_energy_j(),
            ));
        }
        for result in &results[1..] {
            for (actual, expected) in result.0.iter().zip(&results[0].0) {
                close(*actual, *expected, 2e-7, 1e-12);
            }
            close(result.1, results[0].1, 2e-7, 1e-15);
            close(result.2, results[0].2, 2e-7, 1e-15);
            close(result.3, results[0].3, 2e-7, 1e-15);
        }
    }
}
