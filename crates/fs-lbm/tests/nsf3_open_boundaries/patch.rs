//! Partial-aperture controls and their original conservative quadrature.
use super::*;
use fs_lbm::nsf3::{NsfBoundaryDrive, NsfCharacteristicPatch, NsfOpenState, NsfReservoirAperture};
use std::cell::Cell;

const PATCH: NsfCharacteristicPatch = NsfCharacteristicPatch {
    start: [1, 0],
    end: [3, 4],
};
const AMBIENT: NsfBoundaryDrive = NsfBoundaryDrive::PressureOffset(0.0);

fn grid(axis: usize, speed: f64, walls: bool) -> NsfGrid3<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let state = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let velocity_unit = (3.0 * state.specific_gas_constant * state.temperature).sqrt();
    let mut dims = [4; 3];
    dims[axis] = 8;
    let mut source = spec(gas, dims, axis, 1);
    if walls {
        source.periodic = [false; 3];
    }
    let mut velocity = [0.0; 3];
    velocity[axis] = speed / velocity_unit;
    NsfGrid3::new(source, |_| false, |_| (1.0, velocity, 300.0))
        .unwrap()
        .with_characteristic_axis(boundary(axis, speed))
        .unwrap()
}

fn selected(cell: [usize; 3], axis: usize) -> bool {
    let tangent: Vec<_> = (0..3).filter(|&a| a != axis).collect();
    (PATCH.start[0]..PATCH.end[0]).contains(&cell[tangent[0]])
        && (PATCH.start[1]..PATCH.end[1]).contains(&cell[tangent[1]])
}

#[test]
fn aperture_pressure_and_complement_are_distinct_on_every_axis_and_side() {
    for axis in 0..3 {
        for high in [false, true] {
            let mut source = grid(axis, 0.0, false);
            source.set_boundary_drives([AMBIENT; 2]).unwrap();
            let mass = source.total_mass_kg();
            let energy = source.total_energy_j();
            let mut source = source
                .with_characteristic_patch(high, PATCH, NsfBoundaryDrive::PressureOffset(3.0))
                .unwrap();
            assert_eq!(source.characteristic_patch(), Some((high, PATCH)));
            assert_eq!(source.open_patch_traces().unwrap().len(), 8);
            assert_eq!(source.total_mass_kg(), mass);
            assert_eq!(source.total_energy_j(), energy);
            for _ in 0..3 {
                for trace in source.open_boundary_traces(high).unwrap() {
                    assert_eq!(
                        trace.pressure_pa,
                        100000.0 + if selected(trace.cell, axis) { 3.0 } else { 0.0 }
                    );
                }
                assert!(
                    source
                        .open_boundary_traces(!high)
                        .unwrap()
                        .iter()
                        .all(|trace| trace.pressure_pa == 100000.0)
                );
                source.step().unwrap();
            }
            let patch = source.last_patch_exchange().unwrap();
            let complete = source.last_boundary_exchange().unwrap()[usize::from(high)];
            assert_eq!(patch.mean_acoustic_pressure_pa, 3.0);
            assert_eq!(complete.mean_acoustic_pressure_pa, 1.5);
            let work = 3.0 * patch.normal_volume_flow_m3_s * source.law().dt_s;
            assert!(
                (patch.acoustic_pressure_work_out_j - work).abs() < 1e-13 * work.abs().max(1e-30)
            );
        }
    }
}

#[test]
fn uniform_patch_quadrature_has_its_actual_mass_enthalpy_and_pressure_traction() {
    let state = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    for axis in 0..3 {
        for high in [false, true] {
            let mut source = grid(axis, 5.0, false)
                .with_characteristic_patch(
                    high,
                    PATCH,
                    NsfBoundaryDrive::IncomingCharacteristic(0.0),
                )
                .unwrap();
            let area = 8.0 * source.law().dx_m.powi(2);
            let dt = source.law().dt_s;
            let sign = if high { 1.0 } else { -1.0 };
            source.step().unwrap();
            let patch = source.last_patch_exchange().unwrap();
            let expected_mass = sign * state.density * 5.0 * area * dt;
            let expected_energy = expected_mass * (state.specific_heat_cp * 300.0 + 12.5);
            let expected_momentum = sign * (state.pressure + state.density * 25.0) * area * dt;
            assert!((patch.mass_out_kg - expected_mass).abs() < 2e-11 * expected_mass.abs());
            assert!(
                (patch.total_energy_out_j - expected_energy).abs() < 2e-10 * expected_energy.abs()
            );
            assert!(
                (patch.momentum_out_kg_m_s[axis] - expected_momentum).abs()
                    < 2e-11 * expected_momentum.abs()
            );
            assert!((patch.normal_volume_flow_m3_s - sign * 5.0 * area).abs() < 2e-11 * 5.0 * area);
        }
    }
}

#[test]
fn identical_aperture_drive_retains_unpatched_native_motion_and_complete_accounts() {
    let drives = [NsfBoundaryDrive::PressureOffset(2.0), AMBIENT];
    let mut plain = grid(0, 0.0, false);
    plain.set_boundary_drives(drives).unwrap();
    let mut selected = grid(0, 0.0, false)
        .with_characteristic_patch(true, PATCH, AMBIENT)
        .unwrap();
    selected.set_boundary_drives(drives).unwrap();
    for _ in 0..5 {
        assert_eq!(plain.step().unwrap(), selected.step().unwrap());
        assert_eq!(
            plain.last_boundary_exchange(),
            selected.last_boundary_exchange()
        );
        for x in 0..8 {
            for y in 0..4 {
                for z in 0..4 {
                    assert_eq!(plain.cell([x, y, z]), selected.cell([x, y, z]));
                }
            }
        }
    }
    assert!(plain.characteristic_patch().is_none());
    assert!(plain.last_patch_exchange().is_none());
    assert!(plain.open_patch_traces().is_none());
}

#[test]
fn patch_and_escape_preserve_complete_mass_energy_with_transverse_wall_corners() {
    for axis in 0..3 {
        let mut source = grid(axis, 0.0, true)
            .with_characteristic_patch(true, PATCH, NsfBoundaryDrive::PressureOffset(-1.0))
            .unwrap();
        source
            .set_boundary_drives([NsfBoundaryDrive::PressureOffset(2.0), AMBIENT])
            .unwrap();
        for _ in 0..6 {
            let mass = source.total_mass_kg();
            let energy = source.total_energy_j();
            let ledger = source.step().unwrap();
            let faces = source.last_boundary_exchange().unwrap();
            let mass_out: f64 = faces.iter().map(|face| face.mass_out_kg).sum();
            let energy_out: f64 = faces.iter().map(|face| face.total_energy_out_j).sum();
            assert!((source.total_mass_kg() - mass + mass_out).abs() < 2e-12 * mass);
            assert!((source.total_energy_j() - energy + energy_out).abs() < 2e-12 * energy);
            assert!((ledger.open_boundary_energy_j + energy_out).abs() < 2e-12 * energy);
            let patch = source.last_patch_exchange().unwrap();
            assert!(patch.mass_out_kg.is_finite() && patch.total_energy_out_j.is_finite());
            assert_eq!(patch.mean_acoustic_pressure_pa, -1.0);
        }
    }
}

#[test]
fn complementary_apertures_partition_original_nonuniform_exchange() {
    for axis in 0..3 {
        let left = NsfCharacteristicPatch {
            start: [0, 0],
            end: [2, 4],
        };
        let right = NsfCharacteristicPatch {
            start: [2, 0],
            end: [4, 4],
        };
        let low = NsfBoundaryDrive::PressureOffset(2.0);
        let positive = NsfBoundaryDrive::PressureOffset(3.0);
        let negative = NsfBoundaryDrive::PressureOffset(-1.0);
        let mut a = grid(axis, 0.0, true)
            .with_characteristic_patch(true, left, AMBIENT)
            .unwrap();
        let mut b = grid(axis, 0.0, true)
            .with_characteristic_patch(true, right, AMBIENT)
            .unwrap();
        a.set_boundary_drives_with_patch([low, negative], positive)
            .unwrap();
        b.set_boundary_drives_with_patch([low, positive], negative)
            .unwrap();
        for _ in 0..6 {
            assert_eq!(a.step().unwrap(), b.step().unwrap());
            assert_eq!(a.last_boundary_exchange(), b.last_boundary_exchange());
            let first = a.last_patch_exchange().unwrap();
            let second = b.last_patch_exchange().unwrap();
            let full = &a.last_boundary_exchange().unwrap()[1];
            // Link accounts subtract equilibrium populations, so their
            // roundoff scale is the original population sum, not a nearly
            // zero net mass or tangential momentum.
            let face_mass = 16.0 * a.law().rho_ref_kg_m3 * a.law().dx_m.powi(3);
            let face_momentum = face_mass * a.law().dx_m / a.law().dt_s;
            let partition = |x: f64, y: f64, total: f64, roundoff_scale: f64| {
                let scale = x
                    .abs()
                    .max(y.abs())
                    .max(total.abs())
                    .max(roundoff_scale)
                    .max(1e-30);
                assert!(
                    (x + y - total).abs() <= 2e-12 * scale,
                    "axis {axis}: {x} + {y} != {total}"
                );
            };
            partition(
                first.mass_out_kg,
                second.mass_out_kg,
                full.mass_out_kg,
                face_mass,
            );
            partition(
                first.total_energy_out_j,
                second.total_energy_out_j,
                full.total_energy_out_j,
                0.0,
            );
            partition(
                first.pressure_work_out_j,
                second.pressure_work_out_j,
                full.pressure_work_out_j,
                0.0,
            );
            partition(
                first.acoustic_pressure_work_out_j,
                second.acoustic_pressure_work_out_j,
                full.acoustic_pressure_work_out_j,
                0.0,
            );
            partition(
                first.normal_volume_flow_m3_s,
                second.normal_volume_flow_m3_s,
                full.normal_volume_flow_m3_s,
                0.0,
            );
            partition(
                0.5 * first.mean_pressure_pa,
                0.5 * second.mean_pressure_pa,
                full.mean_pressure_pa,
                0.0,
            );
            partition(
                0.5 * first.mean_acoustic_pressure_pa,
                0.5 * second.mean_acoustic_pressure_pa,
                full.mean_acoustic_pressure_pa,
                0.0,
            );
            for component in 0..3 {
                partition(
                    first.momentum_out_kg_m_s[component],
                    second.momentum_out_kg_m_s[component],
                    full.momentum_out_kg_m_s[component],
                    face_momentum,
                );
            }
        }
    }
}

#[test]
fn patch_trials_retain_controls_accounts_and_cancellation_atomicity() {
    let make = || {
        grid(0, 0.0, false)
            .with_characteristic_patch(true, PATCH, AMBIENT)
            .unwrap()
    };
    let mut direct = make();
    let mut preview = make();
    let drives = [NsfBoundaryDrive::PressureOffset(2.0), AMBIENT];
    let patch_drive = NsfBoundaryDrive::IncomingCharacteristic(0.01);
    direct
        .set_boundary_drives_with_patch(drives, patch_drive)
        .unwrap();
    let mut sum = [0.0; 4];
    for _ in 0..3 {
        direct.step().unwrap();
        let patch = direct.last_patch_exchange().unwrap();
        sum[0] += patch.mean_acoustic_pressure_pa;
        sum[1] += patch.normal_volume_flow_m3_s;
        sum[2] += patch.acoustic_pressure_work_out_j;
        sum[3] += patch.total_energy_out_j;
    }
    let trial = preview
        .preview_boundary_interval_with_patch(drives, patch_drive, 3, 384, || false)
        .unwrap();
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    let patch = trial.patch_exchange().unwrap();
    assert_eq!(patch.mean_acoustic_pressure_pa, sum[0] / 3.0);
    assert_eq!(patch.normal_volume_flow_m3_s, sum[1] / 3.0);
    assert_eq!(patch.acoustic_pressure_work_out_j, sum[2]);
    assert_eq!(patch.total_energy_out_j, sum[3]);
    trial.accept();
    assert_eq!(preview.last_patch_exchange(), direct.last_patch_exchange());
    let before = (
        preview.state_hash(),
        preview.steps(),
        preview.cumulative_ledger(),
        preview.boundary_drives(),
        preview.characteristic_patch_drive(),
        preview.last_patch_exchange().copied(),
        preview.open_boundary_traces(true),
    );
    drop(
        preview
            .preview_boundary_interval_with_patch(drives, AMBIENT, 2, 256, || false)
            .unwrap(),
    );
    let calls = Cell::new(0);
    assert!(
        preview
            .preview_boundary_interval_with_patch(drives, AMBIENT, 3, 384, || {
                let n = calls.get();
                calls.set(n + 1);
                n >= 2
            })
            .is_err()
    );
    assert!(calls.get() >= 3);
    assert!(
        preview
            .set_boundary_drives_with_patch(drives, NsfBoundaryDrive::PressureOffset(f64::NAN))
            .is_err()
    );
    assert_eq!(
        before,
        (
            preview.state_hash(),
            preview.steps(),
            preview.cumulative_ledger(),
            preview.boundary_drives(),
            preview.characteristic_patch_drive(),
            preview.last_patch_exchange().copied(),
            preview.open_boundary_traces(true)
        )
    );
    // The original interval entry point preserves, rather than dropping, the patch.
    let trial = preview
        .preview_boundary_interval(drives, 1, 128, || false)
        .unwrap();
    assert_eq!(trial.grid().characteristic_patch_drive(), Some(patch_drive));
    assert!(trial.patch_exchange().is_some());
}

#[test]
fn characteristic_patch_retains_disjoint_reservoir_escape_and_its_trial_accounts() {
    let reservoir_gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.02, 100000.0).unwrap();
    let escape = NsfReservoirAperture {
        axis: 1,
        high: false,
        cells: vec![[3, 0, 1], [3, 0, 2], [4, 0, 1], [4, 0, 2]],
        state: NsfOpenState {
            density_kg_m3: reservoir_gas.density,
            velocity_m_s: [0.0; 3],
            temperature_k: 300.02,
        },
    };
    let make = || {
        grid(0, 0.0, true)
            .with_characteristic_patch(true, PATCH, AMBIENT)
            .unwrap()
            .with_reservoir_aperture(escape.clone())
            .unwrap()
    };
    let mut direct = make();
    let mut preview = make();
    let drives = [NsfBoundaryDrive::PressureOffset(2.0), AMBIENT];
    let drive = NsfBoundaryDrive::PressureOffset(-1.0);
    direct
        .set_boundary_drives_with_patch(drives, drive)
        .unwrap();
    let mut escape_mass = 0.0;
    for _ in 0..3 {
        direct.step().unwrap();
        escape_mass += direct.last_reservoir_mass_out_kg().unwrap();
    }
    let before = (
        preview.state_hash(),
        preview.steps(),
        preview.cumulative_ledger(),
        preview.last_reservoir_mass_out_kg(),
        preview.last_patch_exchange().copied(),
    );
    let trial = preview
        .preview_boundary_interval_with_patch(drives, drive, 3, 384, || false)
        .unwrap();
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    assert_eq!(trial.grid().reservoir_aperture(), Some(&escape));
    assert_eq!(trial.grid().characteristic_patch(), Some((true, PATCH)));
    assert_eq!(trial.reservoir_mass_out_kg(), escape_mass);
    assert!(escape_mass.abs() > 1e-12 * trial.previous_mass_kg());
    assert!(trial.patch_exchange().is_some());
    let mass_out = trial
        .boundary_exchange()
        .iter()
        .map(|face| face.mass_out_kg)
        .sum::<f64>()
        + trial.reservoir_mass_out_kg();
    let mass_error = trial.grid().total_mass_kg() - trial.previous_mass_kg() + mass_out;
    assert!(mass_error.abs() < 2e-12 * trial.previous_mass_kg());
    let ledger = trial.energy_ledger();
    assert_eq!(ledger.viscous_heating_j, 0.0);
    assert!(ledger.open_energy_j > 0.0);
    let energy_error = trial.grid().total_energy_j()
        - trial.previous_energy_j()
        - ledger.open_boundary_energy_j
        - ledger.open_energy_j;
    assert!(energy_error.abs() < 2e-12 * trial.previous_energy_j());
    drop(trial);
    assert_eq!(
        before,
        (
            preview.state_hash(),
            preview.steps(),
            preview.cumulative_ledger(),
            preview.last_reservoir_mass_out_kg(),
            preview.last_patch_exchange().copied()
        )
    );
    assert_eq!(preview.reservoir_aperture(), Some(&escape));
    preview
        .preview_boundary_interval_with_patch(drives, drive, 3, 384, || false)
        .unwrap()
        .accept();
    assert_eq!(preview.state_hash(), direct.state_hash());
    assert_eq!(
        preview.last_reservoir_mass_out_kg(),
        direct.last_reservoir_mass_out_kg()
    );
    assert_eq!(preview.last_patch_exchange(), direct.last_patch_exchange());
    // Either cold attachment order preserves both independent boundary owners.
    let reversed = grid(0, 0.0, true)
        .with_reservoir_aperture(escape.clone())
        .unwrap()
        .with_characteristic_patch(true, PATCH, AMBIENT)
        .unwrap();
    assert_eq!(reversed.reservoir_aperture(), Some(&escape));
    assert_eq!(reversed.characteristic_patch(), Some((true, PATCH)));
}

#[test]
fn invalid_patch_topology_and_late_attachment_refuse() {
    for patch in [
        NsfCharacteristicPatch {
            start: [1, 0],
            end: [1, 2],
        },
        NsfCharacteristicPatch {
            start: [0, 0],
            end: [5, 2],
        },
        NsfCharacteristicPatch {
            start: [3, 2],
            end: [2, 3],
        },
    ] {
        assert!(
            grid(0, 0.0, false)
                .with_characteristic_patch(true, patch, AMBIENT)
                .is_err()
        );
    }
    let mut late = grid(0, 0.0, false);
    late.step().unwrap();
    assert!(
        late.with_characteristic_patch(true, PATCH, AMBIENT)
            .is_err()
    );
    let duplicate = grid(0, 0.0, false)
        .with_characteristic_patch(true, PATCH, AMBIENT)
        .unwrap();
    assert!(
        duplicate
            .with_characteristic_patch(false, PATCH, AMBIENT)
            .is_err()
    );
    let mut plain = grid(0, 0.0, false);
    assert!(
        plain
            .preview_boundary_interval_with_patch([AMBIENT; 2], AMBIENT, 1, 128, || false)
            .is_err()
    );
}
