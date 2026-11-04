//! G1/G3/G4: actual NSF reflection, source impedance and conservative fluxes.
use super::*;
use fs_lbm::nsf3::{NsfBoundaryDrive, NsfCharacteristicPatch};

fn load(resistance: f64, pressure: f64) -> NsfBoundaryDrive {
    NsfBoundaryDrive::ResistiveImpedance {
        specific_impedance_pa_s_m: resistance,
        pressure_offset_pa: pressure,
    }
}

fn grid(axis: usize, temperature: f64, pressure: f64) -> NsfGrid3<D3q27> {
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
        face.reference_temperature_k = temperature;
        face.reference_pressure_pa = pressure;
        face.reference_velocity_m_s[axis] = 0.1;
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
fn g1_impedance_uses_live_outgoing_trace_original_medium_and_outward_normal() {
    for axis in 0..3 {
        for (temperature, pressure) in [(300.0, 100000.0), (310.0, 100000.0), (300.0, 200000.0)] {
            let mut grid = grid(axis, temperature, pressure);
            let reference = grid.characteristic_reference_gas(false).unwrap();
            let z0 = reference.density * reference.sound_speed;
            for ratio in [0.0, 0.25, 1.0, 4.0, 1e12] {
                let resistance = ratio * z0;
                let before = grid.cell([0; 3]);
                let energy = grid.total_energy_j();
                let drives = [load(resistance, 8.0), load(resistance, -3.0)];
                grid.set_boundary_drives(drives).unwrap();
                assert_eq!(grid.boundary_drives(), Some(drives));
                assert_eq!(grid.total_energy_j(), energy);
                let after = grid.cell([0; 3]);
                assert!((after.0 - before.0).abs() < 1e-13);
                assert!((after.2 - before.2).abs() < 1e-11);
                assert!((after.1[axis] - before.1[axis]).abs() < 1e-13);
                for (side, source) in [8.0, -3.0].into_iter().enumerate() {
                    let sign = if side == 0 { -1.0 } else { 1.0 };
                    for face in grid.open_boundary_traces(side != 0).unwrap() {
                        let (rho, u, t) = grid.cell(face.cell);
                        let pi = rho * reference.density * reference.specific_gas_constant * t;
                        let outgoing = sign * (u[axis] * grid.law().dx_m / grid.law().dt_s - 0.1)
                            + (pi - pressure) / z0;
                        let du = (z0 * outgoing - source) / (resistance + z0);
                        let dp = (resistance * outgoing + source) / (ratio + 1.0);
                        assert!(
                            (face.outward_normal_velocity_m_s - sign * 0.1 - du).abs()
                                < 1e-11 * du.abs().max(1.0)
                        );
                        assert!((face.pressure_pa - pressure - dp).abs() < 2e-8);
                    }
                }
            }
        }
    }
}

#[test]
fn g1_passive_load_keeps_native_mass_energy_and_pressure_work_accounts() {
    let mut grid = NsfGrid3::<D3q27>::new(
        spec(GasSpec::dry_air_ussa1976(), [12, 2, 2], 0, 1),
        |_| false,
        |_| (1.0001, [0.0; 3], 300.0),
    )
    .unwrap()
    .with_characteristic_axis(boundary(0, 0.0))
    .unwrap();
    let gas = grid.characteristic_reference_gas(false).unwrap();
    let resistance = 2.0 * gas.density * gas.sound_speed;
    grid.set_boundary_drives([load(resistance, 0.0); 2])
        .unwrap();
    for _ in 0..12 {
        let mass = grid.total_mass_kg();
        let energy = grid.total_energy_j();
        let ledger = grid.step().unwrap();
        let exchange = grid.last_boundary_exchange().unwrap();
        let mass_out = exchange.iter().map(|f| f.mass_out_kg).sum::<f64>();
        let energy_out = exchange.iter().map(|f| f.total_energy_out_j).sum::<f64>();
        assert!((grid.total_mass_kg() - mass + mass_out).abs() < 5e-12 * mass);
        assert!((grid.total_energy_j() - energy + energy_out).abs() < 5e-12 * energy);
        assert!((ledger.open_boundary_energy_j + energy_out).abs() < 5e-12 * energy);
        for (side, account) in exchange.iter().enumerate() {
            assert!(account.acoustic_pressure_work_out_j >= 0.0);
            for face in grid.open_boundary_traces(side != 0).unwrap() {
                let dp = face.pressure_pa - gas.pressure;
                assert!((dp - resistance * face.outward_normal_velocity_m_s).abs() < 2e-8);
            }
        }
    }
}

// Integrate the signed returning pressure packet after reflection and before
// it reaches the opposite matched face. Molecular attenuation is retained.
fn returned_packet(n: usize, ratio: f64, high: bool) -> f64 {
    let gas = GasSpec::dry_air_ussa1976();
    let g = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let c = (g.gamma / 3.0).sqrt();
    let sign = if high { 1.0 } else { -1.0 };
    let mut grid = NsfGrid3::<D3q27>::new(
        spec(gas, [n, 1, 1], 0, 1),
        |_| false,
        |x| {
            let distance = if high { x[0] } else { n - 1 - x[0] } as f64;
            let z = (distance - 0.3 * n as f64) / (0.07 * n as f64);
            let amplitude = 1e-5 * fs_math::det::exp(-z * z);
            (
                1.0 + amplitude / g.gamma,
                [sign * c * amplitude / g.gamma, 0.0, 0.0],
                300.0 * (1.0 + (g.gamma - 1.0) / g.gamma * amplitude),
            )
        },
    )
    .unwrap()
    .with_characteristic_axis(boundary(0, 0.0))
    .unwrap();
    let mut drives = [NsfBoundaryDrive::IncomingCharacteristic(0.0); 2];
    drives[usize::from(high)] = load(ratio * g.density * g.sound_speed, 0.0);
    grid.set_boundary_drives(drives).unwrap();
    let integrated_pressure = |grid: &NsfGrid3<D3q27>| -> f64 {
        (0..n)
            .map(|x| {
                let (rho, _, t) = grid.cell([x, 0, 0]);
                rho * t / 300.0 - 1.0
            })
            .sum()
    };
    let incident = integrated_pressure(&grid);
    for _ in 0..(1.15 * n as f64 / c).round() as usize {
        grid.step().unwrap();
    }
    integrated_pressure(&grid) / incident
}

#[test]
fn g2_pressure_packet_reflects_with_the_analytic_signed_impedance_ratio() {
    for (ratio, expected) in [
        (0.0, -1.0),
        (0.25, -0.6),
        (1.0, 0.0),
        (4.0, 0.6),
        (1e12, 1.0),
    ] {
        let coarse = returned_packet(64, ratio, true);
        let fine = returned_packet(128, ratio, true);
        assert!(
            (fine - expected).abs() < 0.04,
            "R/Z0={ratio}: {coarse} -> {fine}; exact={expected}"
        );
        assert!(
            (fine - expected).abs() <= (coarse - expected).abs() + 0.003,
            "refinement R/Z0={ratio}: {coarse} -> {fine}"
        );
    }
    let low = returned_packet(96, 4.0, false);
    let high = returned_packet(96, 4.0, true);
    assert!(
        (low - high).abs() < 1e-7,
        "outward-normal reflection {low} versus {high}"
    );
}

#[test]
fn g4_impedance_patches_and_trials_preserve_accepted_state_on_refusal_or_drop() {
    let mut grid = grid(0, 300.0, 100000.0)
        .with_characteristic_patch(
            true,
            NsfCharacteristicPatch {
                start: [0, 0],
                end: [1, 1],
            },
            load(800.0, 1.0),
        )
        .unwrap();
    let drives = [load(300.0, 0.0), load(1200.0, -2.0)];
    grid.set_boundary_drives(drives).unwrap();
    let before = (
        grid.state_hash(),
        grid.boundary_drives(),
        grid.open_boundary_traces(true),
        grid.steps(),
    );
    for resistance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            grid.set_boundary_drives([load(resistance, 0.0); 2])
                .is_err()
        );
        assert!(
            grid.set_boundary_drives_with_patch(drives, load(resistance, 0.0))
                .is_err()
        );
        assert_eq!(
            before,
            (
                grid.state_hash(),
                grid.boundary_drives(),
                grid.open_boundary_traces(true),
                grid.steps()
            )
        );
    }
    let count = grid.cell_count();
    let trial = grid
        .preview_boundary_interval(drives, 3, 3 * count, || false)
        .unwrap();
    assert_eq!(trial.grid().boundary_drives(), Some(drives));
    assert_eq!(
        trial.grid().characteristic_patch_drive(),
        Some(load(800.0, 1.0))
    );
    for face in trial.grid().open_boundary_traces(true).unwrap() {
        let (resistance, source) = if face.cell[1] == 0 && face.cell[2] == 0 {
            (800.0, 1.0)
        } else {
            (1200.0, -2.0)
        };
        let expected = resistance * (face.outward_normal_velocity_m_s - 0.1) + source;
        assert!((face.pressure_pa - 100000.0 - expected).abs() < 2e-8);
    }
    drop(trial);
    assert_eq!(
        before,
        (
            grid.state_hash(),
            grid.boundary_drives(),
            grid.open_boundary_traces(true),
            grid.steps()
        )
    );
    let cancelled = grid.preview_boundary_interval(drives, 3, 3 * count, || true);
    assert!(cancelled.is_err());
    drop(cancelled);
    assert_eq!(
        before,
        (
            grid.state_hash(),
            grid.boundary_drives(),
            grid.open_boundary_traces(true),
            grid.steps()
        )
    );
    grid.preview_boundary_interval(drives, 3, 3 * count, || false)
        .unwrap()
        .accept();
    assert_eq!(grid.steps(), 3);
    assert_eq!(grid.boundary_drives(), Some(drives));
    grid.set_incoming_characteristics([0.0; 2]).unwrap();
    assert_eq!(
        grid.boundary_drives(),
        Some([NsfBoundaryDrive::IncomingCharacteristic(0.0); 2])
    );
}
