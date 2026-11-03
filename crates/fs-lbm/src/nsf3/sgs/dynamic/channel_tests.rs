use super::*;
use crate::d3q27::D3q27;
use crate::sgs3::SgsUse;
use fs_material::thermochemical_gas::energy::GasEnergySolve;

#[path = "../../../../tests/nsf3_thermochemical/source.rs"]
mod source;

const SHAPE: [usize; 3] = [5, 4, 5];

fn spec(workers: usize, wall: WallThermal) -> NsfSpec {
    NsfSpec {
        dims: SHAPE,
        periodic: [true, true, false],
        wall,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored dynamic-channel coupling fixture",
        maximum_cells: 4096,
        workers,
        transport_points: 9,
    }
}

fn config() -> NsfChannelDynamicThermalSgsConfig {
    NsfChannelDynamicThermalSgsConfig {
        momentum: ChannelDynamicSgsConfig {
            maximum_coefficient_squared: 16.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: ChannelDynamicScalarConfig {
            maximum_coefficient: 16.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn noise(i: usize, salt: u64) -> f64 {
    let mut value = (i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^= value >> 31;
    ((value >> 11) as f64 / (1_u64 << 53) as f64) * 2.0 - 1.0
}

fn mixed(p: [usize; 3]) -> (f64, [f64; 3], f64) {
    let i = (p[2] * SHAPE[1] + p[1]) * SHAPE[0] + p[0];
    let height = p[2] as f64 + 0.5;
    let envelope = 4.0 * height * (SHAPE[2] as f64 - height) / (SHAPE[2] * SHAPE[2]) as f64;
    (
        1.0,
        core::array::from_fn(|a| 0.002 * envelope * noise(i, 31 * a as u64 + 7)),
        293.15 + 0.1 * (1.0 - height / SHAPE[2] as f64 + 0.05 * envelope * noise(i, 101)),
    )
}

fn rest(_: [usize; 3]) -> (f64, [f64; 3], f64) {
    (1.0, [0.0; 3], 293.15)
}

fn workspace(grid: &NsfGrid3<D3q27>) -> std::sync::MutexGuard<'_, ChannelWorkspace> {
    let DynamicState::Channel { workspace, .. } = &grid.dynamic_sgs.as_ref().unwrap().state else {
        panic!("expected channel workspace");
    };
    workspace.lock().unwrap()
}

#[test]
fn g1_channel_nsf_uses_current_source_heat_capacity_force_and_moving_walls() {
    let gas = ThermochemicalNsfGas::new(
        source::model(-4000.0, 250000.0, 1.0),
        source::point(300.0, 100000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 64,
        },
    )
    .unwrap();
    let mut settings = spec(1, WallThermal::Isothermal(293.15)).map_gas(|_| gas);
    settings.temperature_k = 300.0;
    settings.pressure_pa = 100000.0;
    settings.temperature_window_k = [250.0, 600.0];
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        settings,
        |_| false,
        mixed,
        Some(Box::new(|p, _| {
            (0.0, [0.0001 * (p[1] as f64 + 0.5).sin(), 0.0, 0.0], 0.0)
        })),
    )
    .unwrap()
    .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
    .unwrap();
    let mut flow = ChannelDynamicSgs::new(SHAPE, 1.0, 2, usize::MAX).unwrap();
    let mut heat = ChannelDynamicScalar::new(SHAPE, 2, usize::MAX).unwrap();
    for wall_speed in [0.0, 0.0003] {
        grid.set_edge_velocity(2, true, [wall_speed, 0.0, 0.0])
            .unwrap();
        let (rho, mom) = grid.macro_fields(&grid.f);
        let velocity = grid.velocity_field(&rho, &mom, &grid.force);
        let walls = [ScalarWall::Value(grid.law.theta_of_kelvin(293.15)); 2];
        let reports = flow
            .prepare_with_scalar(
                config().momentum,
                config().heat,
                &rho,
                &velocity,
                grid.edge_velocity[2],
                &grid.theta,
                walls,
                &mut heat,
                || false,
            )
            .unwrap();
        assert!(reports.0.maximum_eddy_viscosity > 0.0);
        assert!(reports.1.maximum_diffusivity > 0.0);
        for i in 0..grid.n {
            let local = grid.material[i].as_ref().unwrap();
            let cp_over_r = local.state.specific_heat_cp / grid.law.gas_constant;
            assert!(
                (cp_over_r - (3.0 + 0.001 * grid.law.kelvin_of_theta(grid.theta[i]))).abs() < 1e-12
            );
            assert_eq!(
                grid.sgs_transport[i].eddy_viscosity_lu,
                flow.viscosities()[i]
            );
            assert_eq!(
                grid.sgs_transport[i].turbulent_conductivity_lu,
                rho[i] * cp_over_r * heat.diffusivities()[i]
            );
            assert_eq!(
                grid.sgs_rates(
                    rho[i],
                    grid.theta[i],
                    Some(local),
                    Some(&grid.sgs_transport[i])
                )
                .1,
                grid.rates(rho[i], grid.theta[i], Some(local)).1
            );
        }
        let physical = flow.viscosities().to_vec();
        let raw: Vec<_> = mom
            .iter()
            .zip(&rho)
            .map(|(m, r)| m.map(|v| v / r))
            .collect();
        flow.prepare_with_scalar(
            config().momentum,
            config().heat,
            &rho,
            &raw,
            grid.edge_velocity[2],
            &grid.theta,
            walls,
            &mut heat,
            || false,
        )
        .unwrap();
        assert!(
            flow.viscosities()
                .iter()
                .zip(physical)
                .any(|(a, b)| (a - b).abs() > 1e-10)
        );
    }
}

#[test]
fn g1_channel_nsf_heat_only_planes_keep_molecular_wall_flux_and_shared_face_limit() {
    let grid = NsfGrid3::<D3q27>::new(spec(1, WallThermal::Isothermal(293.15)), |_| false, mixed)
        .unwrap()
        .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
        .unwrap();
    assert!(
        grid.sgs_transport
            .iter()
            .any(|t| t.eddy_viscosity_lu == 0.0 && t.turbulent_conductivity_lu > 0.0)
    );
    let (rho, mom) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &mom, &grid.force);
    let stress = vec![[0.0; 6]; grid.n];
    let molecular = grid
        .flux_divergence(&grid.e, &rho, &velocity, &grid.theta, &stress, &[], None)
        .unwrap();
    let modeled = grid
        .flux_divergence(
            &grid.e,
            &rho,
            &velocity,
            &grid.theta,
            &stress,
            &grid.sgs_transport,
            None,
        )
        .unwrap();
    assert_eq!(modeled.1, molecular.1);
    assert!(
        modeled
            .0
            .iter()
            .zip(&molecular.0)
            .any(|(a, b)| (a - b).abs() > 1e-14)
    );
    assert!((modeled.0.iter().sum::<f64>() - molecular.0.iter().sum::<f64>()).abs() < 1e-14);
    // A heat-only fit must still be admitted against the actual shared-face
    // molecular source; there is no momentum-based shortcut past this bound.
    let cell = grid
        .sgs_transport
        .iter()
        .position(|t| t.eddy_viscosity_lu == 0.0 && t.turbulent_conductivity_lu > 0.0)
        .unwrap();
    let modeled = grid.sgs_transport[cell].turbulent_conductivity_lu;
    assert!(matches!(
        grid.admit_dynamic_heat_ratio(cell, modeled, modeled / 2000.0),
        Err(NsfError::SgsHeatRatio { .. })
    ));
}

#[test]
fn g1_g5_channel_nsf_reuses_buffers_and_balances_moving_wall_energy_across_workers() {
    for wall in [WallThermal::Adiabatic, WallThermal::Isothermal(293.15)] {
        let make = |workers| {
            let mut grid = NsfGrid3::<D3q27>::new(spec(workers, wall), |_| false, mixed).unwrap();
            grid.set_edge_velocity(2, true, [0.0003, 0.0, 0.0]).unwrap();
            grid.with_channel_dynamic_thermal_sgs(config(), usize::MAX)
                .unwrap()
        };
        let mut serial = make(1);
        let mut parallel = make(3);
        let allocation = {
            let w = workspace(&serial);
            (
                w.momentum.viscosities().as_ptr(),
                w.heat.diffusivities().as_ptr(),
            )
        };
        for _ in 0..4 {
            let energy = serial.total_energy_j();
            let ledger = serial.step().unwrap();
            assert_eq!(ledger, parallel.step().unwrap());
            assert_eq!(serial.state_hash(), parallel.state_hash());
            assert_eq!(serial.sgs_transport, parallel.sgs_transport);
            assert!(ledger.residual_rel < 1e-12);
            assert!(
                (serial.total_energy_j() - energy - ledger.wall_heat_j - ledger.wall_work_j).abs()
                    < energy.abs() * 2e-12
            );
            if matches!(wall, WallThermal::Adiabatic) {
                assert_eq!(ledger.wall_heat_j, 0.0);
            }
            let w = workspace(&serial);
            assert_eq!(w.momentum.viscosities().as_ptr(), allocation.0);
            assert_eq!(w.heat.diffusivities().as_ptr(), allocation.1);
        }
        // A freely decaying channel need not retain a positive scalar fit.
        // Establish the refusal trigger on the authored initial field: an
        // uncapped twin must produce nonzero heat after both requested actions.
        let mut serial = make(1);
        let mut admitted = make(1);
        admitted
            .set_edge_velocity(2, true, [0.0004, 0.0, 0.0])
            .unwrap();
        assert!(
            admitted
                .sgs_transport
                .iter()
                .any(|t| t.turbulent_conductivity_lu > 0.0),
            "wall-change fixture must have positive modeled heat: {wall:?}"
        );
        admitted
            .set_edge_velocity(2, true, [0.0003, 0.0, 0.0])
            .unwrap();
        admitted.step().unwrap();
        assert!(
            admitted
                .sgs_transport
                .iter()
                .any(|t| t.turbulent_conductivity_lu > 0.0),
            "step fixture must have positive modeled heat: {wall:?}"
        );
        let previous = (
            serial.state_hash(),
            serial.sgs_transport.clone(),
            serial.cumulative_ledger(),
            serial.steps(),
            serial.edge_velocity,
            serial.div.clone(),
            serial.wall_out,
        );
        let DynamicState::Channel { config, .. } = &mut serial.dynamic_sgs.as_mut().unwrap().state
        else {
            unreachable!();
        };
        config.heat.maximum_diffusivity_to_molecular = 0.0;
        let wall_result = serial.set_edge_velocity(2, true, [0.0004, 0.0, 0.0]);
        assert!(
            matches!(&wall_result, Err(NsfError::SgsHeatRatio { .. })),
            "expected heat admission refusal for {wall:?}, got {wall_result:?}"
        );
        let step_result = serial.step();
        assert!(
            matches!(&step_result, Err(NsfError::SgsHeatRatio { .. })),
            "expected heat admission refusal for {wall:?}, got {step_result:?}"
        );
        assert_eq!(
            (
                serial.state_hash(),
                serial.sgs_transport.clone(),
                serial.cumulative_ledger(),
                serial.steps(),
                serial.edge_velocity,
                serial.div.clone(),
                serial.wall_out
            ),
            previous
        );
    }
}

#[test]
fn g0_g4_channel_nsf_zero_control_topology_and_workspace_admissions() {
    for axis in 0..3 {
        let mut settings = spec(1, WallThermal::Adiabatic);
        settings.periodic = [true; 3];
        settings.periodic[axis] = false;
        let make = || NsfGrid3::<D3q27>::new(settings.clone(), |_| false, rest).unwrap();
        let bytes = 200 * SHAPE.iter().product::<usize>();
        assert_eq!(
            make().required_dynamic_sgs_workspace_bytes().unwrap(),
            bytes
        );
        assert!(matches!(
            make().with_channel_dynamic_thermal_sgs(config(), bytes - 1),
            Err(NsfError::DynamicSgs(DynamicSgsError::WorkspaceBytes { .. }))
        ));
        let mut plain = make();
        let mut fitted = make()
            .with_channel_dynamic_thermal_sgs(config(), bytes)
            .unwrap();
        assert_eq!(fitted.channel_dynamic_thermal_sgs_config(), Some(config()));
        assert_eq!(fitted.dynamic_thermal_sgs_config(), None);
        assert_eq!(fitted.dynamic_sgs_workspace_bytes(), Some(bytes));
        for _ in 0..2 {
            assert_eq!(plain.step().unwrap(), fitted.step().unwrap());
            assert_eq!(plain.state_hash(), fitted.state_hash());
        }
        let hash = fitted.state_hash();
        assert!(fitted.set_open_face(axis, false, None).is_err());
        assert_eq!(fitted.state_hash(), hash);
        let characteristic = NsfCharacteristicAxis {
            axis,
            faces: [NsfCharacteristicBoundary {
                reference_pressure_pa: 101_325.0,
                reference_temperature_k: 293.15,
                reference_velocity_m_s: [0.0; 3],
                incoming_characteristic_m_s: 0.0,
            }; 2],
            maximum_relative_perturbation: 0.05,
        };
        assert!(
            make()
                .with_channel_dynamic_thermal_sgs(config(), bytes)
                .unwrap()
                .with_characteristic_axis(characteristic)
                .is_err()
        );
        assert!(
            make()
                .with_characteristic_axis(characteristic)
                .unwrap()
                .with_channel_dynamic_thermal_sgs(config(), bytes)
                .is_err()
        );
    }
    for periodic in [[true; 3], [false; 3], [true, false, false]] {
        let mut settings = spec(1, WallThermal::Adiabatic);
        settings.periodic = periodic;
        assert!(
            NsfGrid3::<D3q27>::new(settings, |_| false, rest)
                .unwrap()
                .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
                .is_err()
        );
    }
    for axis in 0..3 {
        let mut settings = spec(1, WallThermal::Adiabatic);
        settings.dims[axis] = 2;
        assert!(
            NsfGrid3::<D3q27>::new(settings, |_| false, rest)
                .unwrap()
                .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
                .is_err()
        );
    }
    assert!(
        NsfGrid3::<D3q27>::new(spec(1, WallThermal::Adiabatic), |p| p == [1, 1, 1], rest)
            .unwrap()
            .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
            .is_err()
    );
    let mut acceptance = config();
    acceptance.heat.usage = SgsUse::AcceptanceRender;
    assert!(
        NsfGrid3::<D3q27>::new(spec(1, WallThermal::Adiabatic), |_| false, rest)
            .unwrap()
            .with_channel_dynamic_thermal_sgs(acceptance, usize::MAX)
            .is_err()
    );
}
