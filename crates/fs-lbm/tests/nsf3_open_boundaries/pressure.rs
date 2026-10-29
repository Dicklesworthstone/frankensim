//! Prescribed pressure uses the original characteristic and conservative owners.
use super::*;
use fs_lbm::nsf3::NsfBoundaryDrive;

fn pressure_grid(axis: usize, temperature: f64, pressure: f64) -> NsfGrid3<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&gas, temperature, pressure).unwrap();
    let velocity_unit = (3.0 * reference.specific_gas_constant * temperature).sqrt();
    let mut dims = [2; 3];
    dims[axis] = 8;
    let mut description = spec(gas, dims, axis, 1);
    description.temperature_k = temperature;
    description.pressure_pa = pressure;
    let mut open = boundary(axis, 0.0);
    for face in &mut open.faces {
        face.reference_pressure_pa = pressure;
        face.reference_temperature_k = temperature;
    }
    NsfGrid3::new(
        description,
        |_| false,
        |_| {
            let mut velocity = [0.0; 3];
            velocity[axis] = 0.4 / velocity_unit;
            (1.0001, velocity, temperature)
        },
    )
    .unwrap()
    .with_characteristic_axis(open)
    .unwrap()
}

#[test]
fn g1_pressure_face_retains_outgoing_wave_and_material_impedance_on_each_axis() {
    for axis in 0..3 {
        for (temperature, pressure) in [(300.0, 100000.0), (310.0, 100000.0), (300.0, 200000.0)] {
            let mut grid = pressure_grid(axis, temperature, pressure);
            let reference = grid.characteristic_reference_gas(false).unwrap();
            let velocity_unit = grid.law().dx_m / grid.law().dt_s;
            let before = grid.cell([0; 3]);
            let energy = grid.total_energy_j();
            let drives = [
                NsfBoundaryDrive::PressureOffset(8.0),
                NsfBoundaryDrive::PressureOffset(-3.0),
            ];
            grid.set_boundary_drives(drives).unwrap();
            assert_eq!(grid.boundary_drives(), Some(drives));
            assert_eq!(grid.total_energy_j(), energy);
            let after = grid.cell([0; 3]);
            assert!((before.0 - after.0).abs() < 1e-13);
            assert!((before.2 - after.2).abs() < 1e-11);
            for (a, b) in before.1.into_iter().zip(after.1) {
                assert!((a - b).abs() < 1e-13);
            }
            for (side, offset) in [8.0, -3.0].into_iter().enumerate() {
                let sign = if side == 0 { -1.0 } else { 1.0 };
                for face in grid.open_boundary_traces(side == 1).unwrap() {
                    let (rho, velocity, t) = grid.cell(face.cell);
                    let p_cell = rho * reference.density * reference.specific_gas_constant * t;
                    let outgoing = sign * velocity[axis] * velocity_unit
                        + (p_cell - pressure) / (reference.density * reference.sound_speed);
                    let expected = outgoing - offset / (reference.density * reference.sound_speed);
                    assert_eq!(face.pressure_pa, pressure + offset);
                    assert!((face.outward_normal_velocity_m_s - expected).abs() < 1e-10);
                }
            }
            // The same returning field under a held incoming wave is a
            // different boundary problem, not a pressure-calibration shortcut.
            let pressure_flow =
                grid.open_boundary_traces(false).unwrap()[0].outward_normal_velocity_m_s;
            grid.set_incoming_characteristics([0.0; 2]).unwrap();
            assert_eq!(
                grid.boundary_drives(),
                Some([NsfBoundaryDrive::IncomingCharacteristic(0.0); 2])
            );
            let wave_flow =
                grid.open_boundary_traces(false).unwrap()[0].outward_normal_velocity_m_s;
            assert!((pressure_flow - wave_flow).abs() > 0.01);
        }
    }
}

#[test]
fn g1_fixed_pressure_has_original_integrated_work_and_mass_energy_balance() {
    let mut grid = pressure_grid(0, 300.0, 100000.0);
    grid.set_boundary_drives([
        NsfBoundaryDrive::PressureOffset(20.0),
        NsfBoundaryDrive::PressureOffset(0.0),
    ])
    .unwrap();
    for _ in 0..8 {
        let mass = grid.total_mass_kg();
        let energy = grid.total_energy_j();
        let ledger = grid.step().unwrap();
        let exchange = grid.last_boundary_exchange().unwrap();
        let mass_out = exchange.iter().map(|f| f.mass_out_kg).sum::<f64>();
        let energy_out = exchange.iter().map(|f| f.total_energy_out_j).sum::<f64>();
        assert!((grid.total_mass_kg() - mass + mass_out).abs() < 5e-12 * mass);
        assert!((grid.total_energy_j() - energy + energy_out).abs() < 5e-12 * energy);
        assert!((ledger.open_boundary_energy_j + energy_out).abs() < 5e-12 * energy);
        for (face, pressure) in exchange.iter().zip([20.0, 0.0]) {
            assert!((face.mean_acoustic_pressure_pa - pressure).abs() < 1e-9);
            let expected = pressure * face.normal_volume_flow_m3_s * grid.law().dt_s;
            assert!(
                (face.acoustic_pressure_work_out_j - expected).abs()
                    <= 1e-12 * expected.abs().max(1e-24)
            );
        }
    }
}

#[test]
fn g0_pressure_refusal_and_cancelled_trial_publish_no_control_or_state() {
    let mut grid = pressure_grid(0, 300.0, 100000.0);
    grid.set_boundary_drives([
        NsfBoundaryDrive::PressureOffset(10.0),
        NsfBoundaryDrive::IncomingCharacteristic(0.0),
    ])
    .unwrap();
    let before = (
        grid.state_hash(),
        grid.boundary_drives(),
        grid.open_boundary_traces(false),
        grid.open_boundary_traces(true),
        grid.total_energy_j(),
        grid.steps(),
    );
    for pressure in [f64::NAN, f64::INFINITY, -100000.0, 10000.0] {
        assert!(
            grid.set_boundary_drives([
                NsfBoundaryDrive::PressureOffset(pressure),
                NsfBoundaryDrive::IncomingCharacteristic(0.0)
            ])
            .is_err()
        );
        assert_eq!(
            before,
            (
                grid.state_hash(),
                grid.boundary_drives(),
                grid.open_boundary_traces(false),
                grid.open_boundary_traces(true),
                grid.total_energy_j(),
                grid.steps()
            )
        );
    }
    let calls = std::cell::Cell::new(0);
    let result = grid.preview_boundary_interval(
        [
            NsfBoundaryDrive::PressureOffset(20.0),
            NsfBoundaryDrive::IncomingCharacteristic(0.0),
        ],
        4,
        4 * grid.cell_count(),
        || {
            calls.set(calls.get() + 1);
            calls.get() >= 3
        },
    );
    assert!(result.is_err());
    drop(result);
    assert!(calls.get() >= 3);
    assert_eq!(
        before,
        (
            grid.state_hash(),
            grid.boundary_drives(),
            grid.open_boundary_traces(false),
            grid.open_boundary_traces(true),
            grid.total_energy_j(),
            grid.steps()
        )
    );
}

#[test]
fn g5_typed_characteristic_and_existing_preview_keep_identical_physical_bits() {
    let mut old = pressure_grid(0, 300.0, 100000.0);
    let mut typed = pressure_grid(0, 300.0, 100000.0);
    for value in [-0.02, 0.03, 0.0] {
        let count = old.cell_count();
        old.preview_characteristic_interval([value, 0.0], 2, 2 * count, || false)
            .unwrap()
            .accept();
        typed
            .preview_boundary_interval(
                [
                    NsfBoundaryDrive::IncomingCharacteristic(value),
                    NsfBoundaryDrive::IncomingCharacteristic(0.0),
                ],
                2,
                2 * count,
                || false,
            )
            .unwrap()
            .accept();
        assert_eq!(old.state_hash(), typed.state_hash());
        assert_eq!(old.boundary_drives(), typed.boundary_drives());
        assert_eq!(old.last_boundary_exchange(), typed.last_boundary_exchange());
    }
}
