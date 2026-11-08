use super::*;
use crate::d3q27::D3q27;
use crate::sgs3::SgsUse;

const SHAPE: [usize; 3] = [5, 6, 7];
const CELLS: usize = 210;

fn spec(periodic: [bool; 3]) -> NsfSpec {
    NsfSpec {
        dims: SHAPE,
        periodic,
        wall: WallThermal::Adiabatic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored Lagrangian NSF fixture",
        maximum_cells: CELLS,
        workers: 1,
        transport_points: 9,
    }
}

fn config() -> NsfLagrangianThermalSgsConfig {
    NsfLagrangianThermalSgsConfig {
        momentum: LagrangianSgsConfig {
            memory_time: 4.0,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: LagrangianScalarConfig {
            memory_time: 1.5,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        adaptive_momentum_memory: Some(AdaptiveLagrangianMemory {
            coefficient: 0.15,
            minimum_time: 0.25,
        }),
    }
}

fn initial(p: [usize; 3]) -> (f64, [f64; 3], f64) {
    let [x, y, z] =
        core::array::from_fn(|a| core::f64::consts::TAU * (p[a] as f64 + 0.5) / SHAPE[a] as f64);
    (
        1.0 + 0.01 * x.sin(),
        [
            0.006 * x.sin() + 0.004 * (y + z).sin(),
            0.005 * (x - z).cos() + 0.003 * y.sin(),
            0.004 * (x + y).sin() + 0.002 * z.cos(),
        ],
        293.15 + 0.4 * (x + y).sin() + 0.2 * z.cos(),
    )
}

fn make(periodic: [bool; 3]) -> NsfGrid3<D3q27> {
    NsfGrid3::new(spec(periodic), |_| false, initial)
        .unwrap()
        .with_lagrangian_thermal_sgs(config(), 264 * CELLS)
        .unwrap()
}

#[test]
fn g1_native_endpoint_matches_one_joint_history_update_and_conserves_energy() {
    let mut grid = make([true; 3]);
    assert_eq!(
        grid.required_lagrangian_sgs_workspace_bytes().unwrap(),
        264 * CELLS
    );
    assert_eq!(grid.dynamic_sgs_workspace_bytes(), Some(264 * CELLS));
    let mut reference = LagrangianDynamicScalar::new(SHAPE, [1.0; 3], 264 * CELLS)
        .unwrap()
        .with_adaptive_momentum_memory(config().adaptive_momentum_memory.unwrap())
        .unwrap();
    let (rho, mom) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &mom, &grid.force);
    reference
        .prepare(
            config().momentum,
            config().heat,
            0.0,
            &rho,
            &velocity,
            &grid.theta,
            || false,
        )
        .unwrap()
        .commit();
    assert_eq!(grid.lagrangian_sgs_accepted_time().unwrap(), Some(0.0));
    for step in 1..=3 {
        let ledger = grid.step().unwrap();
        assert!(ledger.residual_rel < 2e-13, "{ledger:?}");
        assert_eq!(
            grid.lagrangian_sgs_accepted_time().unwrap(),
            Some(step as f64)
        );
        let (rho, mom) = grid.macro_fields(&grid.f);
        let velocity = grid.velocity_field(&rho, &mom, &grid.force);
        let trial = reference
            .prepare(
                config().momentum,
                config().heat,
                step as f64,
                &rho,
                &velocity,
                &grid.theta,
                || false,
            )
            .unwrap();
        for i in 0..CELLS {
            assert_eq!(
                grid.sgs_transport[i].eddy_viscosity_lu,
                trial.viscosities()[i]
            );
            let expected =
                rho[i] * (grid.law.gamma / (grid.law.gamma - 1.0)) * trial.diffusivities()[i];
            assert_eq!(grid.sgs_transport[i].turbulent_conductivity_lu, expected);
        }
        assert_eq!(
            grid.lagrangian_sgs_memory_time_range().unwrap(),
            Some(trial.momentum_memory_time_range())
        );
        trial.commit();
    }
}

#[test]
fn g4_final_cancellation_discards_both_histories_and_retry_replays() {
    let mut direct = make([true; 3]);
    let mut checkpoints = 0;
    direct
        .step_body_acceleration(None, &mut || {
            checkpoints += 1;
            false
        })
        .unwrap();
    assert!(checkpoints > 2);
    let mut grid = make([true; 3]);
    let before = grid.state_hash();
    let transport = grid.sgs_transport.clone();
    let mut reached = 0;
    assert!(
        grid.step_body_acceleration(None, &mut || {
            reached += 1;
            reached == checkpoints
        })
        .is_err()
    );
    assert_eq!(reached, checkpoints);
    assert_eq!(grid.state_hash(), before);
    assert_eq!(grid.sgs_transport, transport);
    assert_eq!(grid.lagrangian_sgs_accepted_time().unwrap(), Some(0.0));
    grid.step().unwrap();
    assert_eq!(grid.state_hash(), direct.state_hash());
    assert_eq!(grid.sgs_transport, direct.sgs_transport);
    assert_eq!(grid.lagrangian_sgs_accepted_time().unwrap(), Some(1.0));
}

#[test]
fn g3_channel_wall_refresh_and_interval_history_are_independent() {
    let mut grid = make([true, true, false]);
    let mut direct = make([true, true, false]);
    for owner in [&mut grid, &mut direct] {
        owner.set_edge_velocity(2, true, [0.004, 0.0, 0.0]).unwrap();
        assert_eq!(owner.lagrangian_sgs_accepted_time().unwrap(), Some(0.0));
    }
    direct.set_thermal_walls(Vec::new()).unwrap();
    direct.step().unwrap();
    direct.step().unwrap();
    let before = grid.state_hash();
    let before_transport = grid.sgs_transport.clone();
    let trial = grid
        .preview_thermal_wall_interval(Vec::new(), 2, CELLS, 2 * CELLS, || false)
        .unwrap();
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    assert_eq!(trial.grid().sgs_transport, direct.sgs_transport);
    assert_eq!(
        trial.grid().lagrangian_sgs_accepted_time().unwrap(),
        Some(2.0)
    );
    drop(trial);
    assert_eq!(grid.state_hash(), before);
    assert_eq!(grid.sgs_transport, before_transport);
    assert_eq!(grid.lagrangian_sgs_accepted_time().unwrap(), Some(0.0));
    grid.preview_thermal_wall_interval(Vec::new(), 2, CELLS, 2 * CELLS, || false)
        .unwrap()
        .accept();
    assert_eq!(grid.state_hash(), direct.state_hash());
    assert_eq!(grid.lagrangian_sgs_accepted_time().unwrap(), Some(2.0));
}

#[test]
fn g0_workspace_usage_and_geometry_refuse_before_attachment() {
    let grid = NsfGrid3::<D3q27>::new(spec([true; 3]), |_| false, initial).unwrap();
    assert!(
        grid.with_lagrangian_thermal_sgs(config(), 264 * CELLS - 1)
            .is_err()
    );
    let mut invalid = config();
    invalid.momentum.usage = SgsUse::AcceptanceRender;
    assert!(invalid.validate().is_err());
    let grid = NsfGrid3::<D3q27>::new(spec([true, false, false]), |_| false, initial).unwrap();
    assert!(
        grid.with_lagrangian_thermal_sgs(config(), 264 * CELLS)
            .is_err()
    );
    let mut grid = make([true; 3]);
    grid.steps = 1usize << (usize::BITS.min(54) - 1);
    if usize::BITS >= 54 {
        let before = grid.state_hash();
        assert!(grid.step().is_err());
        assert_eq!(grid.state_hash(), before);
    }
}
