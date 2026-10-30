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
fn g5_matching_initial_entropy_and_dormant_outflow_temperature_keep_their_physical_meaning() {
    let mut original = grid(0, false, 1, 300.0);
    let mut controlled = grid(0, false, 1, 300.0);
    let mut dormant = grid(0, false, 1, 300.0);
    original.set_boundary_drives([NsfBoundaryDrive::PressureOffset(0.0); 2]).unwrap();
    controlled.set_boundary_drives([thermal(0.0, 300.0), thermal(0.0, 309.0)]).unwrap();
    dormant.set_boundary_drives([thermal(0.0, 300.0), thermal(0.0, 308.0)]).unwrap();
    // Matching initial entropy produces the same original boundary material.
    // Later acoustic entropy and a held thermal inlet are distinct conditions.
    assert_eq!(original.f, controlled.f);
    assert_eq!(original.e, controlled.e);
    assert_eq!(original.open_boundary_traces(false), controlled.open_boundary_traces(false));
    assert_eq!(original.open_boundary_traces(true), controlled.open_boundary_traces(true));
    for _ in 0..8 {
        assert_eq!(dormant.f, controlled.f);
        assert_eq!(dormant.e, controlled.e);
        assert_eq!(dormant.open_boundary_traces(false), controlled.open_boundary_traces(false));
        assert_eq!(dormant.open_boundary_traces(true), controlled.open_boundary_traces(true));
        for face in controlled.open_boundary_traces(false).unwrap() {
            assert!(face.outward_normal_velocity_m_s < 0.0);
            assert_eq!(face.temperature_k, 300.0);
        }
        assert!(controlled.open_boundary_traces(true).unwrap().iter().all(|face|
            face.outward_normal_velocity_m_s > 0.0 && (face.temperature_k - 309.0).abs() > 1.0));
        assert_eq!(dormant.step().unwrap(), controlled.step().unwrap());
    }
    // A dormant thermal control still changes the accepted source identity.
    let hash = controlled.state_hash();
    controlled.set_boundary_drives([thermal(0.0, 300.0), thermal(0.0, 308.0)]).unwrap();
    assert_ne!(controlled.state_hash(), hash);
    controlled.set_boundary_drives([NsfBoundaryDrive::PressureOffset(0.0); 2]).unwrap();
    assert_eq!(original.boundary_drives(), controlled.boundary_drives());
    assert_eq!(controlled.characteristic_inflow_temperature_k, [None; 2]);

    let mut warm = grid(0, false, 1, 300.4);
    warm.set_boundary_drives([NsfBoundaryDrive::PressureOffset(0.0); 2]).unwrap();
    assert_convected_donor(&warm, false, true);
    let acoustic = warm.open_boundary_traces(false).unwrap();
    let energy = warm.total_energy_j();
    warm.set_boundary_drives([thermal(0.0, 300.0), thermal(0.0, 309.0)]).unwrap();
    assert_eq!(warm.total_energy_j(), energy);
    for (face, retained) in warm.open_boundary_traces(false).unwrap().iter().zip(acoustic) {
        assert_eq!(face.temperature_k, 300.0);
        assert!((retained.temperature_k - face.temperature_k).abs() > 0.3);
    }
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

fn convected_grid(axis: usize, high: bool, reference_normal: f64) -> NsfGrid3<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&gas, 300.0, 100_000.0).unwrap();
    let unit = (3.0 * reference.specific_gas_constant * reference.temperature).sqrt();
    let mut dims = [4; 3];
    dims[axis] = 8;
    let mut periodic = [true; 3];
    periodic[axis] = false;
    let mut velocity = [0.0; 3];
    velocity[(axis + 1) % 3] = 0.2 / unit;
    velocity[(axis + 2) % 3] = -0.1 / unit;
    let mut faces = [NsfCharacteristicBoundary {
        reference_pressure_pa: reference.pressure, reference_temperature_k: reference.temperature,
        reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
    }; 2];
    faces[usize::from(high)].reference_velocity_m_s[axis] =
        if high { reference_normal } else { -reference_normal };
    NsfGrid3::new(NsfSpec {
        dims, periodic, gas, temperature_k: 300.0, pressure_pa: 100_000.0,
        dx_m: 5e-5, temperature_window_k: [250.0, 600.0],
        maximum_mach: 0.2, maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "manufactured entropy/shear reversal fixture".to_owned(),
        maximum_cells: dims.iter().product(), wall: WallThermal::Adiabatic,
        workers: 1, transport_points: 17,
    }, |_| false, |_| (300.0 / 300.4, velocity, 300.4)).unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis, faces, maximum_relative_perturbation: 0.05,
        }).unwrap()
}

// Invert the acoustic characteristic identities to put the selected face on
// either side of zero normal speed without changing its interior material.
fn convected_drives(
    grid: &NsfGrid3<D3q27>, high: bool, normal: f64, pressure: bool, temperature: Option<f64>,
) -> [NsfBoundaryDrive; 2] {
    let open = grid.open.unwrap();
    let side = usize::from(high);
    let sign = if high { 1.0 } else { -1.0 };
    let mut cell = [0; 3];
    if high { cell[open.axis] = grid.spec.dims[open.axis] - 1; }
    let (rho, u, t) = grid.cell(cell);
    let reference = grid.characteristic_reference_gas(high).unwrap();
    let impedance = reference.density * reference.sound_speed;
    let normal_reference = sign * open.faces[side].reference_velocity_m_s[open.axis];
    let outgoing = sign * u[open.axis] * grid.law.dx_m / grid.law.dt_s - normal_reference
        + (grid.law.pressure_pa(rho, grid.law.theta_of_kelvin(t)) - reference.pressure) / impedance;
    let offset = impedance * (normal_reference + outgoing - normal);
    let mut drives = [NsfBoundaryDrive::IncomingCharacteristic(0.0); 2];
    drives[side] = if let Some(temperature) = temperature {
        thermal(offset, temperature)
    } else if pressure {
        NsfBoundaryDrive::PressureOffset(offset)
    } else {
        NsfBoundaryDrive::IncomingCharacteristic(2.0 * (normal - normal_reference) - outgoing)
    };
    drives
}

fn assert_convected_donor(grid: &NsfGrid3<D3q27>, high: bool, interior: bool) {
    let side = usize::from(high);
    let open = grid.open.unwrap();
    let reference = grid.characteristic_reference_gas(high).unwrap();
    let alpha = (reference.gamma - 1.0) / reference.gamma * reference.temperature / reference.pressure;
    for trace in grid.open_boundary_traces(high).unwrap() {
        let (rho, u, temperature) = grid.cell(trace.cell);
        let entropy = temperature - reference.temperature
            - alpha * (grid.law.pressure_pa(rho, grid.law.theta_of_kelvin(temperature)) - reference.pressure);
        assert!(entropy > 0.3, "fixture must retain a nonzero entropy characteristic");
        let expected = reference.temperature + alpha * (trace.pressure_pa - reference.pressure)
            + if interior { entropy } else { 0.0 };
        assert!((trace.temperature_k - expected).abs() < 1e-10);
        let index = plane_index(trace.cell, grid.spec.dims, open.axis);
        let face = grid.open_faces.as_ref().unwrap().values[side][index];
        for axis in 0..3 {
            if axis == open.axis { continue; }
            let donor = if interior { u[axis] } else {
                open.faces[side].reference_velocity_m_s[axis] * grid.law.dt_s / grid.law.dx_m
            };
            assert!((face.u[axis] - donor).abs() < 1e-15);
        }
        assert!((grid.law.pressure_pa(face.rho, face.theta) - trace.pressure_pa).abs() < 1e-8);
    }
}

#[test]
fn g1_stationary_entropy_and_shear_pass_acoustic_reversal_continuously() {
    for axis in 0..3 {
        for high in [false, true] {
            for pressure in [false, true] {
                let mut grid = convected_grid(axis, high, 0.0);
                let mut pressure_jumps = [0.0; 2];
                for (rung, speed) in [1e-4, 1e-8].into_iter().enumerate() {
                    let mut means = [0.0; 2];
                    for (side, normal) in [-speed, speed].into_iter().enumerate() {
                        let drives = convected_drives(&grid, high, normal, pressure, None);
                        grid.set_boundary_drives(drives).unwrap();
                        assert_convected_donor(&grid, high, true);
                        for trace in grid.open_boundary_traces(high).unwrap() {
                            assert!((trace.outward_normal_velocity_m_s - normal).abs() < 1e-12);
                        }
                        let before = grid.state_hash();
                        let work = 2 * grid.cell_count();
                        let trial = grid.preview_boundary_interval(drives, 2, work, || false).unwrap();
                        let exchange = trial.boundary_exchange();
                        means[side] = exchange[usize::from(high)].mean_acoustic_pressure_pa;
                        let mass_out: f64 = exchange.iter().map(|face| face.mass_out_kg).sum();
                        let energy_out: f64 = exchange.iter().map(|face| face.total_energy_out_j).sum();
                        assert!((trial.grid().total_mass_kg() - trial.previous_mass_kg() + mass_out).abs()
                            < 5e-12 * trial.previous_mass_kg());
                        assert!((trial.grid().total_energy_j() - trial.previous_energy_j() + energy_out).abs()
                            < 5e-12 * trial.previous_energy_j());
                        drop(trial);
                        assert_eq!(grid.state_hash(), before);
                    }
                    pressure_jumps[rung] = (means[1] - means[0]).abs();
                }
                // A finite donor jump would survive shrinking the acoustic
                // perturbation by 10,000. Original time-averaged pressure must not.
                assert!(pressure_jumps[0] > 1e-4);
                assert!(pressure_jumps[1] < 1e-3 * pressure_jumps[0] + 1e-8,
                    "axis={axis} high={high} pressure={pressure}: {pressure_jumps:?}");
            }
        }
    }
}

#[test]
fn g1_nonzero_reference_speed_owns_convected_donor_through_acoustic_reversal() {
    for axis in 0..3 {
        for high in [false, true] {
            for reference_normal in [-0.01, 0.01] {
                for pressure in [false, true] {
                    let mut grid = convected_grid(axis, high, reference_normal);
                    let normal = -reference_normal;
                    let drives = convected_drives(&grid, high, normal, pressure, None);
                    grid.set_boundary_drives(drives).unwrap();
                    assert_convected_donor(&grid, high, reference_normal > 0.0);
                    assert!(grid.open_boundary_traces(high).unwrap().iter().all(|face|
                        (face.outward_normal_velocity_m_s - normal).abs() < 1e-12));
                }
            }
        }
    }
}

#[test]
fn g1_explicit_thermal_source_keeps_actual_inflow_even_at_stationary_reference() {
    for axis in 0..3 {
        for high in [false, true] {
            let mut grid = convected_grid(axis, high, 0.0);
            for normal in [-0.01, 0.01] {
                let drives = convected_drives(&grid, high, normal, true, Some(304.0));
                grid.set_boundary_drives(drives).unwrap();
                if normal > 0.0 {
                    assert_convected_donor(&grid, high, true);
                } else {
                    for trace in grid.open_boundary_traces(high).unwrap() {
                        assert_eq!(trace.temperature_k, 304.0);
                        assert!(trace.outward_normal_velocity_m_s < 0.0);
                    }
                    for face in &grid.open_faces.as_ref().unwrap().values[usize::from(high)] {
                        for tangent in [(axis + 1) % 3, (axis + 2) % 3] {
                            assert_eq!(face.u[tangent], 0.0);
                        }
                    }
                }
            }
        }
    }
}
