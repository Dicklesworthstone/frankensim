use super::*;
use crate::nsf3::{NsfThermalWallFace, NsfWallThermal};

fn grid(workers: usize) -> NsfGrid3<D3q27> {
    let mut grid = NsfGrid3::<D3q27>::new(spec(workers, WallThermal::Adiabatic), |_| false, mixed)
        .unwrap()
        .with_channel_dynamic_thermal_sgs(config(), usize::MAX)
        .unwrap();
    grid.set_edge_velocity(2, true, [0.0003, 0.0, 0.0]).unwrap();
    grid
}

fn walls() -> Vec<NsfThermalWallFace> {
    vec![
        NsfThermalWallFace {
            cell: [1, 1, 0],
            axis: 2,
            high: false,
            thermal: NsfWallThermal::Isothermal(293.25),
        },
        NsfThermalWallFace {
            cell: [3, 2, 0],
            axis: 2,
            high: false,
            thermal: NsfWallThermal::HeatFlux(20.0),
        },
        NsfThermalWallFace {
            cell: [2, 1, 4],
            axis: 2,
            high: true,
            thermal: NsfWallThermal::Film {
                coefficient_w_m2_k: 25.0,
                ambient_temperature_k: 293.05,
            },
        },
    ]
}

fn buffer_addresses(grid: &NsfGrid3<D3q27>) -> (*const f64, *const f64) {
    let work = workspace(grid);
    (
        work.momentum.viscosities().as_ptr(),
        work.heat.diffusivities().as_ptr(),
    )
}

fn accepted(grid: &NsfGrid3<D3q27>) -> (u64, EnergyLedger, Vec<NsfSgsTransport>, Vec<f64>, usize) {
    (
        grid.state_hash(),
        grid.cumulative_ledger(),
        grid.sgs_transport.clone(),
        grid.div.clone(),
        grid.steps(),
    )
}

#[test]
fn dynamic_wall_interval_retains_original_solver_heat_and_shared_scratch_on_drop_and_accept() {
    let steps = 3;
    let mut fluid = grid(1);
    let mut direct = grid(3);
    assert!(
        fluid
            .sgs_transport
            .iter()
            .any(|x| x.turbulent_conductivity_lu > 0.0)
    );
    let original = accepted(&fluid);
    let allocation = buffer_addresses(&fluid);
    let bytes = fluid.dynamic_sgs_workspace_bytes().unwrap();
    direct
        .set_thermal_walls_recorded(walls(), || false)
        .unwrap();
    let native: Vec<_> = (0..steps).map(|_| direct.step().unwrap()).collect();
    let mut expected_exchange: Option<Vec<crate::nsf3::NsfThermalWallExchange>> = None;
    for keep in [false, true] {
        let trial = fluid
            .preview_thermal_wall_interval(walls(), steps, 100, 100 * steps, || false)
            .unwrap();
        assert_eq!(trial.grid().state_hash(), direct.state_hash());
        assert_eq!(trial.grid().sgs_transport, direct.sgs_transport);
        assert_eq!(
            trial.grid().last_coupling_report(),
            direct.last_coupling_report()
        );
        assert_eq!(
            trial.grid().channel_dynamic_thermal_sgs_config(),
            Some(config())
        );
        assert_eq!(trial.grid().dynamic_sgs_workspace_bytes(), Some(bytes));
        assert_eq!(buffer_addresses(trial.grid()), allocation);
        assert_eq!(
            trial.energy_ledger().wall_heat_j,
            native.iter().map(|x| x.wall_heat_j).sum::<f64>()
        );
        assert_eq!(
            trial.energy_ledger().wall_work_j,
            native.iter().map(|x| x.wall_work_j).sum::<f64>()
        );
        let local_heat: f64 = trial
            .wall_exchange()
            .iter()
            .map(|x| x.heat_into_fluid_j)
            .sum();
        assert!((local_heat - trial.energy_ledger().wall_heat_j).abs() < 1e-22);
        assert!(trial.energy_ledger().residual_rel < 3e-12);
        if let Some(expected) = &expected_exchange {
            assert_eq!(trial.wall_exchange(), expected.as_slice());
        } else {
            expected_exchange = Some(trial.wall_exchange().to_vec());
        }
        if keep {
            trial.accept();
        } else {
            drop(trial);
            assert_eq!(accepted(&fluid), original);
            assert_eq!(buffer_addresses(&fluid), allocation);
        }
    }
    assert_eq!(accepted(&fluid), accepted(&direct));
    assert_eq!(buffer_addresses(&fluid), allocation);
}

#[test]
fn dynamic_interval_late_cancellation_and_fit_refusal_preserve_accepted_coefficients() {
    let mut fluid = grid(1);
    let before = accepted(&fluid);
    let allocation = buffer_addresses(&fluid);
    let mut one_step_polls = 0;
    let trial = fluid
        .preview_thermal_wall_interval(walls(), 1, 100, 100, || {
            one_step_polls += 1;
            false
        })
        .unwrap();
    drop(trial);
    let mut polls = 0;
    assert!(
        fluid
            .preview_thermal_wall_interval(walls(), 3, 100, 300, || {
                polls += 1;
                polls > one_step_polls
            })
            .is_err()
    );
    assert!(polls > one_step_polls);
    assert_eq!(accepted(&fluid), before);
    assert_eq!(buffer_addresses(&fluid), allocation);
    assert_eq!(fluid.channel_dynamic_thermal_sgs_config(), Some(config()));
    // The initial original fit is known to have nonzero modeled heat. Its cap
    // refusal must not overwrite accepted coefficients or consume fluid time.
    let DynamicState::Channel { config, .. } = &mut fluid.dynamic_sgs.as_mut().unwrap().state
    else {
        unreachable!();
    };
    config.heat.maximum_diffusivity_to_molecular = 0.0;
    assert!(matches!(
        fluid.preview_thermal_wall_interval(Vec::new(), 2, 100, 200, || false),
        Err(NsfError::SgsHeatRatio { .. })
    ));
    assert_eq!(accepted(&fluid), before);
    assert_eq!(buffer_addresses(&fluid), allocation);
}

#[test]
fn periodic_dynamic_preview_replays_the_original_empty_wall_interval() {
    let mut settings = spec(1, WallThermal::Adiabatic);
    settings.periodic = [true; 3];
    let dynamic = NsfDynamicThermalSgsConfig {
        momentum: DynamicSgsConfig {
            averaging: crate::sgs3::dynamic::DynamicAveraging::Domain,
            maximum_coefficient_squared: 16.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: DynamicScalarConfig {
            averaging: crate::sgs3::dynamic::DynamicAveraging::Domain,
            maximum_coefficient: 16.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    };
    let make = || {
        NsfGrid3::<D3q27>::new(settings.clone(), |_| false, mixed)
            .unwrap()
            .with_dynamic_thermal_sgs(dynamic, usize::MAX)
            .unwrap()
    };
    let mut fluid = make();
    let mut direct = make();
    direct
        .set_thermal_walls_recorded(Vec::new(), || false)
        .unwrap();
    for _ in 0..2 {
        direct.step().unwrap();
    }
    let trial = fluid
        .preview_thermal_wall_interval(Vec::new(), 2, 100, 200, || false)
        .unwrap();
    assert_eq!(trial.grid().dynamic_thermal_sgs_config(), Some(dynamic));
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    assert_eq!(trial.grid().sgs_transport, direct.sgs_transport);
    assert!(trial.wall_exchange().is_empty());
    assert_eq!(trial.energy_ledger().wall_heat_j, 0.0);
    assert_eq!(trial.energy_ledger().wall_work_j, 0.0);
    trial.accept();
    assert_eq!(accepted(&fluid), accepted(&direct));
}
