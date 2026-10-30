use super::*;
use fs_material::gas::GasState;

fn characteristic(axis: usize, velocity_m_s: [f64; 3]) -> NsfCharacteristicAxis {
    NsfCharacteristicAxis {
        axis,
        faces: [NsfCharacteristicBoundary {
            reference_pressure_pa: 101_325.0,
            reference_temperature_k: 293.15,
            reference_velocity_m_s: velocity_m_s,
            incoming_characteristic_m_s: 0.0,
        }; 2],
        maximum_relative_perturbation: 0.05,
    }
}

#[test]
fn uniform_throughflow_has_no_spurious_wall_shear_at_characteristic_faces() {
    for axis in 0..3 {
        let mut periodic = [true; 3];
        periodic[axis] = false;
        let velocity = [0.003, -0.002, 0.004];
        let make = || {
            NsfGrid3::<D3q27>::new(
                spec([6; 3], periodic, WallThermal::Adiabatic),
                |_| false,
                |_| (1.0, velocity, 293.15),
            )
            .unwrap()
        };
        let original = make();
        let reference = velocity.map(|v| v * original.law.dx_m / original.law.dt_s);
        let open = characteristic(axis, reference);
        let config = config(SgsModel::Smagorinsky { cs: 0.17 });
        let mut first = original
            .with_sgs(config)
            .unwrap()
            .with_characteristic_axis(open)
            .unwrap();
        let mut last = make()
            .with_characteristic_axis(open)
            .unwrap()
            .with_sgs(config)
            .unwrap();
        assert_eq!(first.sgs_transport, last.sgs_transport);
        for _ in 0..4 {
            for transport in &first.sgs_transport {
                assert!(transport.eddy_viscosity_lu < 1e-16);
            }
            assert_eq!(first.step().unwrap(), last.step().unwrap());
            assert_eq!(first.state_hash(), last.state_hash());
            for i in 0..first.n {
                let (rho, u, temperature) = first.cell(first.at(i));
                assert!((rho - 1.0).abs() < 1e-12);
                assert!((temperature - 293.15).abs() < 1e-9);
                for a in 0..3 {
                    assert!((u[a] - velocity[a]).abs() < 1e-12);
                }
            }
        }
    }
}

#[test]
fn zero_open_sgs_is_exact_and_sgs_rebind_refusal_keeps_all_fields() {
    let make = || {
        NsfGrid3::<D3q27>::new(
            spec([6, 4, 4], [false, true, true], WallThermal::Adiabatic),
            |_| false,
            rest,
        )
        .unwrap()
        .with_characteristic_axis(characteristic(0, [0.0; 3]))
        .unwrap()
    };
    let mut molecular = make();
    let mut zero = make()
        .with_sgs(config(SgsModel::Vreman { cv: 0.0 }))
        .unwrap();
    let drives = [
        NsfBoundaryDrive::PressureOffset(1.0),
        NsfBoundaryDrive::PressureOffset(0.0),
    ];
    molecular.set_boundary_drives(drives).unwrap();
    zero.set_boundary_drives(drives).unwrap();
    for _ in 0..4 {
        assert_eq!(molecular.step().unwrap(), zero.step().unwrap());
        assert_eq!(molecular.state_hash(), zero.state_hash());
    }
    let mut bounded = config(SgsModel::Smagorinsky { cs: 0.17 });
    bounded.closure.maximum_eddy_to_molecular = 1e-8;
    let mut grid = make().with_sgs(bounded).unwrap();
    let before = (
        grid.state_hash(),
        grid.sgs_transport.clone(),
        grid.div.clone(),
        grid.force.clone(),
        grid.boundary_drives(),
        grid.cumulative_ledger(),
    );
    assert!(matches!(
        grid.set_boundary_drives(drives),
        Err(NsfError::Sgs(SgsError::ViscosityRatioExceeded { .. }))
    ));
    assert_eq!(
        (
            grid.state_hash(),
            grid.sgs_transport.clone(),
            grid.div.clone(),
            grid.force.clone(),
            grid.boundary_drives(),
            grid.cumulative_ledger()
        ),
        before
    );
}

fn aperture_grid(workers: usize, sgs_first: bool) -> NsfGrid3<D3q27> {
    let mut s = spec([8; 3], [false; 3], WallThermal::Adiabatic);
    s.workers = workers;
    let mut grid = NsfGrid3::<D3q27>::new(
        s,
        |_| false,
        |x| {
            let y = core::f64::consts::TAU * (x[1] as f64 + 0.5) / 8.0;
            (1.0, [0.002 * y.sin(), 0.0, 0.0], 293.15 + 0.1 * y.cos())
        },
    )
    .unwrap()
    .with_characteristic_axis(characteristic(0, [0.0; 3]))
    .unwrap();
    let config = config(SgsModel::Smagorinsky { cs: 0.17 });
    if sgs_first {
        grid = grid.with_sgs(config).unwrap();
    }
    grid = grid
        .with_characteristic_patch(
            true,
            NsfCharacteristicPatch {
                start: [2, 2],
                end: [6, 6],
            },
            NsfBoundaryDrive::PressureOffset(0.0),
        )
        .unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.35, 101_325.0).unwrap();
    grid = grid
        .with_reservoir_aperture(NsfReservoirAperture {
            axis: 1,
            high: false,
            cells: vec![[3, 0, 3], [3, 0, 4], [4, 0, 3], [4, 0, 4]],
            state: NsfOpenState {
                density_kg_m3: gas.density,
                velocity_m_s: [0.0; 3],
                temperature_k: gas.temperature,
            },
        })
        .unwrap();
    if !sgs_first {
        grid = grid.with_sgs(config).unwrap();
    }
    grid
}

#[test]
fn open_sgs_heat_flux_is_the_same_flux_in_storage_and_boundary_accounts() {
    let grid = aperture_grid(1, false);
    let (rho, mom) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &mom, &grid.force);
    let stress = vec![[0.0; 6]; grid.n];
    let conductivity = 0.004;
    let transport = vec![
        NsfSgsTransport {
            eddy_viscosity_lu: 0.001,
            turbulent_conductivity_lu: conductivity,
        };
        grid.n
    ];
    let molecular = grid
        .flux_divergence(
            &grid.e,
            &rho,
            &velocity,
            &grid.theta,
            &stress,
            &[],
            grid.open_faces.as_ref(),
        )
        .unwrap();
    let modeled = grid
        .flux_divergence(
            &grid.e,
            &rho,
            &velocity,
            &grid.theta,
            &stress,
            &transport,
            grid.open_faces.as_ref(),
        )
        .unwrap();
    let mut expected_characteristic = [0.0; 2];
    let mut expected_reservoir = 0.0;
    for i in 0..grid.n {
        for axis in 0..3 {
            for (side, direction) in [(0, -1), (1, 1)] {
                if let Some(face) =
                    grid.characteristic_face(grid.open_faces.as_ref(), i, axis, direction)
                {
                    expected_characteristic[side] -=
                        2.0 * conductivity * (face.theta - grid.theta[i]);
                } else if let Some(ghost) = grid.open_at(i, axis, direction) {
                    expected_reservoir -= conductivity * (ghost.theta - grid.theta[i]);
                }
            }
        }
    }
    for side in 0..2 {
        let delta = modeled.2.total_out[side] - molecular.2.total_out[side];
        assert!((delta - expected_characteristic[side]).abs() < 1e-14);
        let internal = modeled.2.internal_out[side] - molecular.2.internal_out[side];
        assert!((internal - expected_characteristic[side]).abs() < 1e-14);
    }
    assert!(expected_reservoir.abs() > 1e-8);
    assert!(
        (modeled.2.prescribed_out - molecular.2.prescribed_out - expected_reservoir).abs() < 1e-14
    );
    assert_eq!(modeled.1.to_bits(), molecular.1.to_bits());
    let divergence_delta: f64 = modeled.0.iter().zip(&molecular.0).map(|(a, b)| a - b).sum();
    let expected = expected_reservoir + expected_characteristic.iter().sum::<f64>();
    assert!((divergence_delta - expected).abs() < 1e-13);
}

#[test]
fn aperture_ghost_gradient_uses_its_full_cell_distance() {
    let mut grid = aperture_grid(1, false);
    // A quadratic transverse velocity has an exact centred derivative at the
    // adjacent centre when the real reservoir ghost is used one spacing away.
    let velocity: Vec<_> = (0..grid.n)
        .map(|i| {
            let y = grid.at(i)[1] as f64;
            [0.001 * (y + 1.0).powi(2), 0.0, 0.0]
        })
        .collect();
    // The retained reservoir velocity is zero at y=-1 for this polynomial.
    let i = grid.index([3, 0, 3]);
    let gradient = grid.sgs_gradient(i, &velocity, grid.open_faces.as_ref());
    assert!((gradient[0][1] - 0.002).abs() < 1e-15);
    // A neighboring impermeable tile uses its halfway no-slip value instead.
    let wall = grid.index([2, 0, 3]);
    let wall_gradient = grid.sgs_gradient(wall, &velocity, grid.open_faces.as_ref());
    assert!((wall_gradient[0][1] - (0.001 + 0.004 / 3.0)).abs() < 1e-15);
    // Updating the acoustic controls recomputes the modeled transport too.
    let before = grid.sgs_transport.clone();
    grid.set_boundary_drives([NsfBoundaryDrive::PressureOffset(1.0); 2])
        .unwrap();
    assert!(grid.sgs_transport.iter().zip(before).any(|(a, b)| *a != b));
}

#[test]
fn coupled_open_sgs_trials_conserve_mass_energy_and_keep_atomic_histories() {
    let mut serial = aperture_grid(1, true);
    let mut parallel = aperture_grid(3, true);
    assert!(serial.sgs_transport == parallel.sgs_transport);
    // Attaching SGS on the other side of the force-chart rebind performs its
    // final moment reconstruction in a different order. That equivalence is
    // physical to roundoff; worker identity below keeps identical operations.
    let other_order = aperture_grid(1, false);
    for (left, right) in serial.sgs_transport.iter().zip(&other_order.sgs_transport) {
        for (a, b) in [
            (left.eddy_viscosity_lu, right.eddy_viscosity_lu),
            (
                left.turbulent_conductivity_lu,
                right.turbulent_conductivity_lu,
            ),
        ] {
            assert!((a - b).abs() < 1e-12 * a.abs().max(b.abs()).max(1e-12));
        }
    }
    let before = (
        serial.state_hash(),
        serial.sgs_transport.clone(),
        serial.cumulative_ledger(),
    );
    assert!(
        serial
            .set_boundary_drives([NsfBoundaryDrive::PressureOffset(100_000.0); 2])
            .is_err()
    );
    assert_eq!(
        (
            serial.state_hash(),
            serial.sgs_transport.clone(),
            serial.cumulative_ledger()
        ),
        before
    );
    let drives = [
        NsfBoundaryDrive::PressureOffset(1.0),
        NsfBoundaryDrive::PressureOffset(0.0),
    ];
    let patch = NsfBoundaryDrive::PressureOffset(-0.5);
    assert!(
        serial
            .preview_boundary_interval_with_patch(drives, patch, 2, 1024, || true)
            .is_err()
    );
    {
        let trial = serial
            .preview_boundary_interval_with_patch(drives, patch, 2, 1024, || false)
            .unwrap();
        assert!(
            trial
                .grid()
                .sgs_transport
                .iter()
                .any(|t| t.eddy_viscosity_lu > 1e-8)
        );
    }
    assert_eq!(
        (
            serial.state_hash(),
            serial.sgs_transport.clone(),
            serial.cumulative_ledger()
        ),
        before
    );
    for _ in 0..4 {
        let a = serial
            .preview_boundary_interval_with_patch(drives, patch, 2, 1024, || false)
            .unwrap();
        let b = parallel
            .preview_boundary_interval_with_patch(drives, patch, 2, 1024, || false)
            .unwrap();
        assert_eq!(a.grid().state_hash(), b.grid().state_hash());
        assert!(a.grid().sgs_transport == b.grid().sgs_transport);
        assert_eq!(a.energy_ledger(), b.energy_ledger());
        let mass_out = a
            .boundary_exchange()
            .iter()
            .map(|e| e.mass_out_kg)
            .sum::<f64>()
            + a.reservoir_mass_out_kg();
        let residual = a.grid().total_mass_kg() - a.previous_mass_kg() + mass_out;
        assert!(residual.abs() < 2e-12 * a.previous_mass_kg());
        let ledger = a.energy_ledger();
        let energy_in = ledger.open_boundary_energy_j + ledger.open_energy_j + ledger.wall_heat_j;
        let residual = a.grid().total_energy_j() - a.previous_energy_j() - energy_in;
        assert!(residual.abs() < 2e-12 * a.previous_energy_j());
        assert!(ledger.residual_rel < 2e-12);
        assert_eq!(a.patch_exchange(), b.patch_exchange());
        a.accept();
        b.accept();
    }
}
