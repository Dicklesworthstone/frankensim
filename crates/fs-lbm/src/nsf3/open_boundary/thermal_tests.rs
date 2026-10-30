//! Incoming thermal data uses the original material and conservative fluxes.
use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{NsfSgsConfig, NsfSpec, WallThermal};
use crate::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_material::gas::{GasSpec, GasState};

fn grid(axis: usize, inlet_high: bool, workers: usize, temperature: f64) -> NsfGrid3<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&gas, 300.0, 100_000.0).unwrap();
    let velocity_unit = (3.0 * reference.specific_gas_constant * 300.0).sqrt();
    let mut dims = [4; 3];
    dims[axis] = 8;
    let mut periodic = [true; 3];
    periodic[axis] = false;
    let mut velocity = [0.0; 3];
    velocity[axis] = if inlet_high { -2.0 } else { 2.0 } / velocity_unit;
    NsfGrid3::new(NsfSpec {
        dims, periodic, gas, temperature_k: 300.0, pressure_pa: 100_000.0,
        dx_m: 5e-5, temperature_window_k: [250.0, 600.0],
        maximum_mach: 0.2, maximum_temperature_ratio: 1.0,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "manufactured dry-air boundary fixture".to_owned(),
        maximum_cells: dims.iter().product(), wall: WallThermal::Adiabatic,
        workers, transport_points: 17,
    }, |_| false, |_| (300.0 / temperature, velocity, temperature)).unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis,
            faces: [NsfCharacteristicBoundary {
                reference_pressure_pa: 100_000.0, reference_temperature_k: 300.0,
                reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
            }; 2],
            maximum_relative_perturbation: 0.05,
        }).unwrap()
}

fn thermal(pressure: f64, temperature: f64) -> NsfBoundaryDrive {
    NsfBoundaryDrive::PressureTemperature {
        pressure_offset_pa: pressure, inflow_temperature_k: temperature,
    }
}

#[test]
fn g1_inflow_temperature_and_outgoing_entropy_use_each_original_face_and_gas() {
    for axis in 0..3 {
        for high in [false, true] {
            let mut grid = grid(axis, high, 1, 300.4);
            let original_energy = grid.total_energy_j();
            let before: Vec<_> = (0..grid.n).map(|i| grid.cell(grid.at(i))).collect();
            let mut drives = [thermal(-3.0, 309.0); 2];
            drives[usize::from(high)] = thermal(8.0, 304.0);
            grid.set_boundary_drives(drives).unwrap();
            assert_eq!(grid.boundary_drives(), Some(drives));
            assert_eq!(grid.total_energy_j(), original_energy);
            for (i, &(density, velocity, temperature)) in before.iter().enumerate() {
                let actual = grid.cell(grid.at(i));
                assert!((actual.0 - density).abs() < 2e-13);
                assert_eq!(actual.2, temperature);
                for a in 0..3 { assert!((actual.1[a] - velocity[a]).abs() < 2e-13); }
            }
            let reference = grid.characteristic_reference_gas(high).unwrap();
            let alpha = (reference.gamma - 1.0) / reference.gamma
                * reference.temperature / reference.pressure;
            for face_high in [false, true] {
                for face in grid.open_boundary_traces(face_high).unwrap() {
                    if face_high == high {
                        assert!(face.outward_normal_velocity_m_s < 0.0);
                        assert!((face.temperature_k - 304.0).abs() < 1e-12);
                        assert_eq!(face.pressure_pa, 100_008.0);
                    } else {
                        assert!(face.outward_normal_velocity_m_s > 0.0);
                        let (rho, _, temperature) = grid.cell(face.cell);
                        let pressure = rho * reference.density * reference.specific_gas_constant * temperature;
                        let expected = reference.temperature + alpha * (face.pressure_pa - reference.pressure)
                            + (temperature - reference.temperature - alpha * (pressure - reference.pressure));
                        assert!((face.temperature_k - expected).abs() < 1e-10);
                        assert!((face.temperature_k - 309.0).abs() > 1.0);
                    }
                    let reconstructed_pressure = face.density_kg_m3
                        * reference.specific_gas_constant * face.temperature_k;
                    assert!((face.pressure_pa - reconstructed_pressure).abs() < 1e-9);
                }
            }
        }
    }
}

#[test]
fn g5_matching_inflow_and_ignored_outflow_temperature_preserve_legacy_physical_bits() {
    let mut original = grid(0, false, 1, 300.0);
    let mut controlled = grid(0, false, 1, 300.0);
    original.set_boundary_drives([NsfBoundaryDrive::PressureOffset(0.0); 2]).unwrap();
    controlled.set_boundary_drives([thermal(0.0, 300.0), thermal(0.0, 309.0)]).unwrap();
    for _ in 0..8 {
        assert_eq!(original.f, controlled.f);
        assert_eq!(original.e, controlled.e);
        assert_eq!(original.open_boundary_traces(false), controlled.open_boundary_traces(false));
        assert_eq!(original.open_boundary_traces(true), controlled.open_boundary_traces(true));
        assert_eq!(original.step().unwrap(), controlled.step().unwrap());
    }
    // A dormant thermal control still changes the accepted source identity.
    let hash = controlled.state_hash();
    controlled.set_boundary_drives([thermal(0.0, 300.0), thermal(0.0, 308.0)]).unwrap();
    assert_ne!(controlled.state_hash(), hash);
    controlled.set_boundary_drives([NsfBoundaryDrive::PressureOffset(0.0); 2]).unwrap();
    assert_eq!(original.boundary_drives(), controlled.boundary_drives());
    assert_eq!(controlled.characteristic_inflow_temperature_k, [None; 2]);
}

#[test]
fn g3_g4_warm_source_trials_conserve_mass_energy_and_replay_with_modeled_transport() {
    for modeled in [false, true] {
        let build = |workers| {
            let grid = grid(0, false, workers, 300.0);
            if modeled {
                grid.with_sgs(NsfSgsConfig {
                    closure: SgsConfig {
                        model: SgsModel::Smagorinsky { cs: 0.17 }, filter_width_lu: 1.0,
                        maximum_eddy_to_molecular: 100.0, usage: SgsUse::Exploration,
                    }, turbulent_prandtl: 0.9,
                }).unwrap()
            } else { grid }
        };
        let mut grid = build(1);
        let mut other_workers = build(3);
        let controls = [thermal(2.0, 301.0), NsfBoundaryDrive::PressureOffset(0.0)];
        let before = (grid.state_hash(), grid.boundary_drives(), grid.total_energy_j(), grid.steps());
        let work = 3 * grid.cell_count();
        let trial = grid.preview_boundary_interval(controls, 3, work, || false).unwrap();
        let mass_out: f64 = trial.boundary_exchange().iter().map(|face| face.mass_out_kg).sum();
        let energy_out: f64 = trial.boundary_exchange().iter().map(|face| face.total_energy_out_j).sum();
        assert!((trial.grid().total_mass_kg() - trial.previous_mass_kg() + mass_out).abs()
            < 5e-12 * trial.previous_mass_kg());
        assert!((trial.grid().total_energy_j() - trial.previous_energy_j() + energy_out).abs()
            < 5e-12 * trial.previous_energy_j());
        let expected = (trial.grid().state_hash(), *trial.energy_ledger(), *trial.boundary_exchange());
        drop(trial);
        assert_eq!((grid.state_hash(), grid.boundary_drives(), grid.total_energy_j(), grid.steps()), before);
        let calls = std::cell::Cell::new(0);
        assert!(grid.preview_boundary_interval(controls, 3, work, || {
            calls.set(calls.get() + 1); calls.get() == 3
        }).is_err());
        assert_eq!((grid.state_hash(), grid.boundary_drives(), grid.total_energy_j(), grid.steps()), before);
        let retry = grid.preview_boundary_interval(controls, 3, work, || false).unwrap();
        assert_eq!((retry.grid().state_hash(), *retry.energy_ledger(), *retry.boundary_exchange()), expected);
        retry.accept();
        other_workers.preview_boundary_interval(controls, 3, work, || false).unwrap().accept();
        assert_eq!(grid.state_hash(), other_workers.state_hash());
        assert_eq!(grid.cumulative_ledger(), other_workers.cumulative_ledger());
        assert!(grid.cell([0; 3]).2 > 300.0, "original gas must receive incoming heat");
    }
}

#[test]
fn g0_patch_thermal_controls_are_independent_and_invalid_sources_are_atomic() {
    let patch = NsfCharacteristicPatch { start: [1, 1], end: [3, 3] };
    let mut grid = grid(0, false, 1, 300.0)
        .with_characteristic_patch(false, patch, thermal(2.0, 304.0)).unwrap();
    let controls = [thermal(2.0, 302.0), thermal(0.0, 309.0)];
    grid.set_boundary_drives(controls).unwrap();
    for face in grid.open_boundary_traces(false).unwrap() {
        let expected = if patch.contains(face.cell, 0) { 304.0 } else { 302.0 };
        assert!((face.temperature_k - expected).abs() < 1e-12);
    }
    let before = (grid.state_hash(), grid.boundary_drives(), grid.characteristic_patch_drive(),
        grid.open_boundary_traces(false), grid.open_boundary_traces(true), grid.total_energy_j());
    for temperature in [f64::NAN, f64::INFINITY, 0.0, 249.0, 320.0, 700.0] {
        // This face is outgoing: invalid future inflow still cannot be retained.
        assert!(grid.set_boundary_drives([controls[0], thermal(0.0, temperature)]).is_err());
        assert!(grid.set_boundary_drives_with_patch(controls, thermal(2.0, temperature)).is_err());
        assert_eq!((grid.state_hash(), grid.boundary_drives(), grid.characteristic_patch_drive(),
            grid.open_boundary_traces(false), grid.open_boundary_traces(true), grid.total_energy_j()), before);
    }
    assert!(grid.set_boundary_drives([thermal(-100_000.0, 300.0), controls[1]]).is_err());
    assert_eq!(grid.state_hash(), before.0);
}
