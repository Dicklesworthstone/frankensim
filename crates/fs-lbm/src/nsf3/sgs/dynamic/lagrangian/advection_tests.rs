use super::*;
use crate::d3q27::D3q27;
use crate::sgs3::SgsUse;

const SHAPE: [usize; 3] = [5, 4, 5];

fn spec(channel: bool) -> NsfSpec {
    NsfSpec {
        dims: SHAPE,
        periodic: [true, true, !channel],
        wall: WallThermal::Adiabatic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored midpoint SGS history fixture",
        maximum_cells: 100,
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
            0.004 * (x.sin() + (y + z).sin()),
            0.003 * ((x - z).cos() + y.sin()),
            0.002 * ((x + y).sin() + z.cos()),
        ],
        293.15 + 0.2 * (x + y).sin() + 0.1 * z.cos(),
    )
}

#[test]
fn g3_native_midpoint_history_forks_preserve_policy_epoch_and_conservative_endpoint() {
    let policy = LagrangianAdvection::Midpoint { substeps: 3 };
    for channel in [false, true] {
        let molecular = || NsfGrid3::<D3q27>::new(spec(channel), |_| false, initial).unwrap();
        assert!(
            molecular()
                .with_lagrangian_thermal_sgs_advection(
                    config(),
                    LagrangianAdvection::Midpoint { substeps: 0 },
                    264 * 100
                )
                .is_err()
        );
        let make = || {
            molecular()
                .with_lagrangian_thermal_sgs_advection(config(), policy, 264 * 100)
                .unwrap()
        };
        let mut tested = make();
        let mut direct = make();
        if channel {
            for grid in [&mut tested, &mut direct] {
                grid.set_edge_velocity(2, true, [0.0004, 0.0, 0.0]).unwrap();
                assert_eq!(grid.lagrangian_sgs_accepted_time().unwrap(), Some(0.0));
            }
        }
        assert_eq!(tested.lagrangian_sgs_advection().unwrap(), Some(policy));
        assert_eq!(tested.dynamic_sgs_workspace_bytes(), Some(264 * 100));
        let before = tested.state_hash();
        assert!(
            tested
                .preview_thermal_wall_step(Vec::new(), 100, || true)
                .is_err()
        );
        assert_eq!(tested.state_hash(), before);
        assert_eq!(tested.lagrangian_sgs_accepted_time().unwrap(), Some(0.0));
        let ledger = direct.step().unwrap();
        let trial = tested
            .preview_thermal_wall_step(Vec::new(), 100, || false)
            .unwrap();
        assert_eq!(
            trial.grid().lagrangian_sgs_advection().unwrap(),
            Some(policy)
        );
        assert_eq!(
            trial.grid().lagrangian_sgs_accepted_time().unwrap(),
            Some(1.0)
        );
        assert_eq!(*trial.energy_ledger(), ledger);
        assert_eq!(trial.grid().state_hash(), direct.state_hash());
        assert_eq!(trial.grid().sgs_transport, direct.sgs_transport);
        drop(trial);
        assert_eq!(tested.state_hash(), before);
        assert_eq!(tested.lagrangian_sgs_accepted_time().unwrap(), Some(0.0));
        tested
            .preview_thermal_wall_step(Vec::new(), 100, || false)
            .unwrap()
            .accept();
        assert_eq!(tested.state_hash(), direct.state_hash());
        assert_eq!(tested.lagrangian_sgs_advection().unwrap(), Some(policy));
        assert_eq!(tested.lagrangian_sgs_accepted_time().unwrap(), Some(1.0));
    }
}
