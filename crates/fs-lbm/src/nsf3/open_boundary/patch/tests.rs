//! Two apertures exercise the original pressure reconstruction, streamed
//! mass/momentum, FV energy and transactional trial on the same NSF owner.
use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{NsfOpenState, NsfReservoirAperture, NsfSpec, WallThermal};
use fs_material::gas::{GasSpec, GasState};
use std::cell::Cell;

fn grid(axis: usize, workers: usize) -> NsfGrid3<D3q27> {
    NsfGrid3::new(NsfSpec {
        dims: [6; 3], periodic: [false; 3], gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0, pressure_pa: 100_000.0, dx_m: 5e-5,
        temperature_window_k: [250.0, 600.0], maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "manufactured two-aperture conservation fixture".to_owned(),
        maximum_cells: 216, wall: WallThermal::Adiabatic, workers, transport_points: 17,
    }, |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis, faces: [NsfCharacteristicBoundary {
                reference_pressure_pa: 100_000.0, reference_temperature_k: 300.0,
                reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
            }; 2], maximum_relative_perturbation: 0.05,
        }).unwrap()
}

fn shapes() -> [NsfCharacteristicPatch; 2] {
    [NsfCharacteristicPatch { start: [1, 1], end: [3, 3] },
        NsfCharacteristicPatch { start: [0, 1], end: [4, 3] }]
}

fn paired(axis: usize, workers: usize) -> NsfGrid3<D3q27> {
    grid(axis, workers)
        .with_characteristic_patch(false, shapes()[0], NsfBoundaryDrive::PressureOffset(0.0)).unwrap()
        .with_characteristic_patch(true, shapes()[1], NsfBoundaryDrive::PressureOffset(0.0)).unwrap()
}

fn pressure(value: f64) -> NsfBoundaryDrive { NsfBoundaryDrive::PressureOffset(value) }

#[test]
fn g1_two_faces_keep_distinct_aperture_pressure_temperature_and_complements() {
    for axis in 0..3 {
        let mut owner = paired(axis, 1);
        let initial = (owner.e.clone(), owner.total_energy_j());
        let drives = [pressure(1.0), pressure(-2.0)];
        let patch_drives = [Some(NsfBoundaryDrive::PressureTemperature {
            pressure_offset_pa: 9.0, inflow_temperature_k: 301.0,
        }), Some(pressure(-7.0))];
        owner.set_boundary_drives_with_patches(drives, patch_drives).unwrap();
        assert_eq!(owner.e, initial.0);
        assert_eq!(owner.total_energy_j(), initial.1);
        assert!(owner.characteristic_patch().is_none());
        assert!(owner.characteristic_patch_drive().is_none());
        assert!(owner.open_patch_traces().is_none());
        for side in 0..2 {
            let high = side == 1;
            assert_eq!(owner.characteristic_patch_on(high), Some(shapes()[side]));
            assert_eq!(owner.characteristic_patch_drive_on(high), patch_drives[side]);
            assert_eq!(owner.open_patch_traces_on(high).unwrap().len(), [4, 8][side]);
            for trace in owner.open_boundary_traces(high).unwrap() {
                let selected = shapes()[side].contains(trace.cell, axis);
                let offset = if selected { [9.0, -7.0][side] } else { [1.0, -2.0][side] };
                assert_eq!(trace.pressure_pa, 100_000.0 + offset);
                if side == 0 && selected {
                    assert!(trace.outward_normal_velocity_m_s < 0.0);
                    assert!((trace.temperature_k - 301.0).abs() < 1e-12);
                }
            }
        }
        owner.set_boundary_drives_with_patches(drives, [Some(pressure(12.0)), None]).unwrap();
        assert_eq!(owner.characteristic_patch_drive_on(false), Some(pressure(12.0)));
        assert_eq!(owner.characteristic_patch_drive_on(true), Some(pressure(-7.0)));
    }
}

#[test]
fn g3_full_face_apertures_reproduce_original_fields_and_both_complete_accounts() {
    let full = NsfCharacteristicPatch { start: [0, 0], end: [6, 6] };
    for axis in 0..3 {
        let mut direct = grid(axis, 1);
        let mut patches = grid(axis, 3)
            .with_characteristic_patch(false, full, pressure(0.0)).unwrap()
            .with_characteristic_patch(true, full, pressure(0.0)).unwrap();
        let drives = [pressure(3.0), pressure(-1.0)];
        direct.set_boundary_drives(drives).unwrap();
        patches.set_boundary_drives_with_patches([pressure(0.0); 2], drives.map(Some)).unwrap();
        for _ in 0..3 {
            assert_eq!(patches.f, direct.f);
            assert_eq!(patches.e, direct.e);
            assert_eq!(patches.step().unwrap(), direct.step().unwrap());
            assert_eq!(patches.last_boundary_exchange(), direct.last_boundary_exchange());
            for side in 0..2 {
                assert_eq!(patches.last_patch_exchange_on(side == 1),
                    Some(&direct.last_boundary_exchange().unwrap()[side]));
            }
        }
    }
}

#[test]
fn g3_two_patch_trials_accumulate_original_work_and_close_full_mass_energy() {
    let mut direct = paired(0, 1);
    let mut previewed = paired(0, 3);
    let drives = [pressure(0.0); 2];
    let patch_drives = [Some(pressure(3.0)), Some(pressure(-1.0))];
    direct.set_boundary_drives_with_patches(drives, patch_drives).unwrap();
    let mut mass = [0.0; 2];
    let mut work = [0.0; 2];
    let mut mean_pressure = [0.0; 2];
    let mut flow = [0.0; 2];
    for _ in 0..3 {
        direct.step().unwrap();
        for side in 0..2 {
            let patch = direct.last_patch_exchange_on(side == 1).unwrap();
            mass[side] += patch.mass_out_kg;
            work[side] += patch.acoustic_pressure_work_out_j;
            mean_pressure[side] += patch.mean_acoustic_pressure_pa;
            flow[side] += patch.normal_volume_flow_m3_s;
        }
    }
    let trial = previewed.preview_boundary_interval_with_patches(
        drives, patch_drives, 3, 648, || false).unwrap();
    assert!(trial.patch_exchange().is_none());
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    for side in 0..2 {
        let patch = trial.patch_exchange_on(side == 1).unwrap();
        assert_eq!(patch.mass_out_kg, mass[side]);
        assert_eq!(patch.acoustic_pressure_work_out_j, work[side]);
        assert_eq!(patch.mean_acoustic_pressure_pa, mean_pressure[side] / 3.0);
        assert_eq!(patch.normal_volume_flow_m3_s, flow[side] / 3.0);
        assert!(patch.acoustic_pressure_work_out_j.abs() > 0.0);
    }
    let mass_residual = trial.grid().total_mass_kg() - trial.previous_mass_kg()
        + trial.boundary_exchange().iter().map(|face| face.mass_out_kg).sum::<f64>();
    assert!(mass_residual.abs() <= 1e-12 * trial.previous_mass_kg());
    assert!(trial.energy_ledger().residual_rel < 1e-12);
    trial.accept();
    assert_eq!(previewed.cumulative_ledger(), direct.cumulative_ledger());
}

#[test]
fn g4_two_patch_refusal_and_late_cancellation_preserve_both_histories() {
    let mut owner = paired(0, 1);
    let drives = [pressure(0.0); 2];
    let controls = [Some(pressure(3.0)), Some(pressure(-1.0))];
    owner.set_boundary_drives_with_patches(drives, controls).unwrap();
    owner.step().unwrap();
    let before = (owner.state_hash(), owner.steps(), owner.cumulative_ledger(),
        owner.boundary_drives(), owner.last_patch_exchange, owner.last_boundary_exchange,
        owner.f.clone(), owner.e.clone());
    drop(owner.preview_boundary_interval_with_patches(drives, controls, 3, 648, || false).unwrap());
    let calls = Cell::new(0);
    assert!(owner.preview_boundary_interval_with_patches(drives, controls, 3, 648, || {
        let n = calls.get(); calls.set(n + 1); n >= 2
    }).is_err());
    assert!(calls.get() >= 3);
    assert!(owner.set_boundary_drives_with_patch(drives, pressure(1.0)).is_err());
    assert!(owner.set_boundary_drives_with_patches(drives,
        [Some(pressure(2.0)), Some(pressure(f64::NAN))]).is_err());
    assert!(owner.set_boundary_drives_with_patches(drives,
        [Some(pressure(2.0)), Some(NsfBoundaryDrive::PressureTemperature {
            pressure_offset_pa: 1.0, inflow_temperature_k: 400.0,
        })]).is_err());
    assert_eq!(before, (owner.state_hash(), owner.steps(), owner.cumulative_ledger(),
        owner.boundary_drives(), owner.last_patch_exchange, owner.last_boundary_exchange,
        owner.f.clone(), owner.e.clone()));
    assert_eq!(owner.characteristic_patch_drive_on(false), controls[0]);
    assert_eq!(owner.characteristic_patch_drive_on(true), controls[1]);
    let mut missing = grid(0, 1);
    assert!(missing.set_boundary_drives_with_patches(drives, controls).is_err());
    assert!(paired(0, 1).with_characteristic_patch(false, shapes()[0], pressure(0.0)).is_err());
}

#[test]
fn g3_both_apertures_and_warm_reservoir_share_one_original_candidate() {
    let reference = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let state = NsfOpenState {
        density_kg_m3: reference.density, velocity_m_s: [0.0; 3], temperature_k: 300.0,
    };
    let mut owner = paired(0, 1).with_reservoir_aperture(NsfReservoirAperture {
        axis: 1, high: false, cells: vec![[2, 0, 2], [3, 0, 2]], state,
    }).unwrap();
    let warm = owner.reservoir_aperture_state_at(100_003.0, 301.0).unwrap();
    let drives = [pressure(0.0); 2];
    let controls = [Some(pressure(3.0)), Some(pressure(-1.0))];
    let hash = owner.state_hash();
    let trial = owner.preview_boundary_interval_with_patches_and_reservoir(
        drives, controls, warm, 2, 432, || false).unwrap();
    assert_eq!(trial.grid().reservoir_aperture().unwrap().state, warm);
    assert!(trial.patch_exchange_on(false).is_some());
    assert!(trial.patch_exchange_on(true).is_some());
    assert!(trial.reservoir_exchange().unwrap().pressure_work_out_j.is_finite());
    let residual = trial.grid().total_mass_kg() - trial.previous_mass_kg()
        + trial.boundary_exchange().iter().map(|face| face.mass_out_kg).sum::<f64>()
        + trial.reservoir_mass_out_kg();
    assert!(residual.abs() < 1e-12 * trial.previous_mass_kg());
    assert!(trial.energy_ledger().residual_rel < 1e-12);
    drop(trial);
    assert_eq!(owner.state_hash(), hash);
    assert_eq!(owner.reservoir_aperture().unwrap().state, state);
    let trial = owner.preview_boundary_interval_with_patches_and_reservoir(
        drives, controls, warm, 2, 432, || false).unwrap();
    trial.accept();
    assert_eq!(owner.reservoir_aperture().unwrap().state, warm);
    assert_eq!(owner.characteristic_patch_drive_on(false), controls[0]);
    assert_eq!(owner.characteristic_patch_drive_on(true), controls[1]);
}
