use super::*;
use crate::d3q27::D3q27;
use crate::sgs3::SgsUse;
use crate::sgs3::dynamic::DynamicAveraging;

const SHAPE: [usize; 3] = [5, 4, 5];

fn spec(normal: Option<usize>) -> NsfSpec {
    NsfSpec {
        dims: SHAPE,
        periodic: core::array::from_fn(|axis| normal != Some(axis)),
        wall: WallThermal::Adiabatic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored sensor-protected dynamic WALE fixture",
        maximum_cells: 4096,
        workers: 1,
        transport_points: 9,
    }
}

fn config() -> NsfDynamicThermalSgsConfig {
    NsfDynamicThermalSgsConfig {
        momentum: DynamicSgsConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient_squared: 16.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: DynamicScalarConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: 16.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn channel_config() -> NsfChannelDynamicThermalSgsConfig {
    NsfChannelDynamicThermalSgsConfig {
        momentum: ChannelDynamicSgsConfig {
            maximum_coefficient_squared: config().momentum.maximum_coefficient_squared,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: ChannelDynamicScalarConfig {
            maximum_coefficient: config().heat.maximum_coefficient,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn policy() -> DynamicWaleConfig {
    DynamicWaleConfig {
        sensor_threshold: 0.09,
        near_wall_coefficient: 0.5,
    }
}

fn initial(p: [usize; 3]) -> (f64, [f64; 3], f64) {
    let i = (p[2] * SHAPE[1] + p[1]) * SHAPE[0] + p[0];
    (
        1.0 + (i % 11) as f64 / 1024.0,
        core::array::from_fn(|a| {
            (((i * (17 + 8 * a) + i * i * (a + 1) + 13 * a) % 127) as f64 - 63.0) / 65_536.0
        }),
        293.15 + ((i * i * 3 + 7 * i) % 41) as f64 / 512.0,
    )
}

#[test]
fn g1_wale_nsf_uses_same_force_corrected_momentum_and_independent_heat() {
    let grid = NsfGrid3::<D3q27>::new_with_sources(
        spec(None),
        |_| false,
        initial,
        Some(Box::new(|p, _| {
            (0.0, [0.0001 * (p[1] as f64 + 0.5).sin(), 0.0, 0.0], 0.0)
        })),
    )
    .unwrap()
    .with_dynamic_wale_thermal_sgs(config(), policy(), usize::MAX)
    .unwrap();
    assert_eq!(grid.dynamic_wale_config().unwrap(), Some(policy()));
    let (rho, mom) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &mom, &grid.force);
    let mut flow = PeriodicDynamicSgs::new(SHAPE, 1.0, usize::MAX)
        .unwrap()
        .with_wale(policy(), usize::MAX)
        .unwrap();
    let mut heat = PeriodicDynamicScalar::new(SHAPE, usize::MAX).unwrap();
    flow.prepare_with_scalar(
        config().momentum,
        config().heat,
        &rho,
        &velocity,
        &grid.theta,
        &mut heat,
        || false,
    )
    .unwrap();
    assert!(flow.viscosities().iter().any(|nu| *nu > 0.0));
    for cell in 0..grid.n {
        let local = grid.material.get(cell).and_then(Option::as_ref);
        let cp_over_r = local.map_or(grid.law.gamma / (grid.law.gamma - 1.0), |m| {
            m.state.specific_heat_cp / grid.law.gas_constant
        });
        assert_eq!(
            grid.sgs_transport[cell].eddy_viscosity_lu,
            flow.viscosities()[cell]
        );
        assert_eq!(
            grid.sgs_transport[cell].turbulent_conductivity_lu,
            rho[cell] * cp_over_r * heat.diffusivities()[cell]
        );
        assert_eq!(
            grid.sgs_rates(
                rho[cell],
                grid.theta[cell],
                local,
                Some(&grid.sgs_transport[cell])
            )
            .1,
            grid.rates(rho[cell], grid.theta[cell], local).1
        );
    }
    let physical = flow.viscosities().to_vec();
    let raw: Vec<_> = mom
        .iter()
        .zip(&rho)
        .map(|(m, rho)| m.map(|v| v / rho))
        .collect();
    flow.prepare(config().momentum, &rho, &raw, || false)
        .unwrap();
    assert!(
        flow.viscosities()
            .iter()
            .zip(physical)
            .any(|(raw, physical)| (raw - physical).abs() > 1e-10)
    );
}

#[test]
fn g1_wale_channel_owner_retains_current_moving_walls_on_every_axis() {
    for normal in 0..3 {
        let mut grid = NsfGrid3::<D3q27>::new(spec(Some(normal)), |_| false, initial)
            .unwrap()
            .with_channel_dynamic_wale_thermal_sgs(channel_config(), policy(), usize::MAX)
            .unwrap();
        let tangential = (normal + 1) % 3;
        let mut flow = ChannelDynamicSgs::new(SHAPE, 1.0, normal, usize::MAX)
            .unwrap()
            .with_wale(policy(), usize::MAX)
            .unwrap();
        let mut heat = ChannelDynamicScalar::new(SHAPE, normal, usize::MAX).unwrap();
        for speed in [0.0, 0.0004] {
            let mut wall = [0.0; 3];
            wall[tangential] = speed;
            grid.set_edge_velocity(normal, true, wall).unwrap();
            assert_eq!(grid.dynamic_wale_config().unwrap(), Some(policy()));
            let (rho, mom) = grid.macro_fields(&grid.f);
            let velocity = grid.velocity_field(&rho, &mom, &grid.force);
            flow.prepare_with_scalar(
                channel_config().momentum,
                channel_config().heat,
                &rho,
                &velocity,
                grid.edge_velocity[normal],
                &grid.theta,
                [ScalarWall::ZeroGradient; 2],
                &mut heat,
                || false,
            )
            .unwrap();
            assert!(flow.viscosities().iter().any(|nu| *nu > 0.0));
            for cell in 0..grid.n {
                let local = grid.material.get(cell).and_then(Option::as_ref);
                let cp_over_r = local.map_or(grid.law.gamma / (grid.law.gamma - 1.0), |m| {
                    m.state.specific_heat_cp / grid.law.gas_constant
                });
                assert_eq!(
                    grid.sgs_transport[cell].eddy_viscosity_lu,
                    flow.viscosities()[cell]
                );
                assert_eq!(
                    grid.sgs_transport[cell].turbulent_conductivity_lu,
                    rho[cell] * cp_over_r * heat.diffusivities()[cell]
                );
            }
        }
    }
}

#[test]
fn g3_wale_native_trials_retain_model_state_budgets_and_publication_rules() {
    let n = SHAPE.iter().product::<usize>();
    let molecular = || NsfGrid3::<D3q27>::new(spec(None), |_| false, initial).unwrap();
    assert_eq!(
        molecular().required_dynamic_sgs_workspace_bytes().unwrap(),
        200 * n
    );
    assert_eq!(
        molecular()
            .required_dynamic_wale_sgs_workspace_bytes()
            .unwrap(),
        216 * n
    );
    assert!(matches!(
        molecular().with_dynamic_wale_thermal_sgs(config(), policy(), 216 * n - 1),
        Err(NsfError::DynamicSgs(DynamicSgsError::WorkspaceBytes { .. }))
    ));
    let mut too_small = config();
    too_small.momentum.maximum_coefficient_squared = 0.24;
    assert!(matches!(
        molecular().with_dynamic_wale_thermal_sgs(too_small, policy(), usize::MAX),
        Err(NsfError::DynamicSgs(DynamicSgsError::Invalid(_)))
    ));
    assert!(
        NsfGrid3::<D3q27>::new(spec(Some(0)), |_| false, initial)
            .unwrap()
            .with_dynamic_wale_thermal_sgs(config(), policy(), usize::MAX)
            .is_err()
    );
    let make = || {
        molecular()
            .with_dynamic_wale_thermal_sgs(config(), policy(), 216 * n)
            .unwrap()
    };
    let mut grid = make();
    let mut direct = make();
    let before = grid.state_hash();
    let transport = grid.sgs_transport.clone();
    assert_eq!(grid.dynamic_sgs_workspace_bytes(), Some(216 * n));
    assert!(
        grid.preview_thermal_wall_step(Vec::new(), n, || true)
            .is_err()
    );
    assert_eq!(grid.state_hash(), before);
    assert_eq!(grid.sgs_transport, transport);
    let expected_ledger = direct.step().unwrap();
    let trial = grid
        .preview_thermal_wall_step(Vec::new(), n, || false)
        .unwrap();
    assert_eq!(trial.grid().dynamic_wale_config().unwrap(), Some(policy()));
    assert_eq!(*trial.energy_ledger(), expected_ledger);
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    assert_eq!(trial.grid().sgs_transport, direct.sgs_transport);
    drop(trial);
    assert_eq!(grid.state_hash(), before);
    assert_eq!(grid.sgs_transport, transport);
    grid.preview_thermal_wall_step(Vec::new(), n, || false)
        .unwrap()
        .accept();
    assert_eq!(grid.state_hash(), direct.state_hash());
    assert_eq!(grid.dynamic_wale_config().unwrap(), Some(policy()));
}
