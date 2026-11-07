use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{NsfCouplingControl, NsfSpec, NsfWallThermal, WallThermal};
use fs_material::gas::GasSpec;

fn spec() -> NsfSpec {
    NsfSpec {
        dims: [4, 4, 4],
        periodic: [false; 3],
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0,
        pressure_pa: 100000.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer N2 ratio; native finite-wall trial fixture",
        maximum_cells: 64,
        wall: WallThermal::Adiabatic,
        workers: 1,
        transport_points: 17,
    }
}

fn grid() -> NsfGrid3<D3q27> {
    NsfGrid3::new(spec(), |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap()
}

fn face(axis: usize, high: bool, thermal: NsfWallThermal) -> NsfThermalWallFace {
    let mut cell = [1; 3];
    cell[axis] = if high { 3 } else { 0 };
    NsfThermalWallFace {
        cell,
        axis,
        high,
        thermal,
    }
}

fn near(a: f64, b: f64, allowance: f64) {
    assert!(
        (a - b).abs() <= allowance,
        "{a:e} != {b:e} (allowance {allowance:e})"
    );
}

#[test]
fn face_joules_keep_physical_orientation_and_native_single_owned_heat() {
    for axis in 0..3 {
        for high in [false, true] {
            for flux in [-250.0, 250.0] {
                let mut fluid = grid();
                let before = fluid.state_hash();
                let expected = flux * fluid.law().dx_m.powi(2) * fluid.law().dt_s;
                let wall = face(axis, high, NsfWallThermal::HeatFlux(flux));
                let trial = fluid
                    .preview_thermal_wall_step(vec![wall], 64, || false)
                    .unwrap();
                let exchange = trial.wall_exchange()[0];
                assert_eq!(
                    (exchange.cell, exchange.axis, exchange.high),
                    (wall.cell, axis, high)
                );
                near(exchange.heat_into_fluid_j, expected, 5e-14 * expected.abs());
                near(
                    trial.energy_ledger().wall_heat_j,
                    expected,
                    5e-14 * expected.abs(),
                );
                near(
                    trial.grid().total_energy_j() - trial.previous_energy_j(),
                    expected,
                    1e-18,
                );
                assert_eq!(trial.grid().steps(), 1);
                assert_eq!(trial.energy_ledger().wall_work_j, 0.0);
                drop(trial);
                assert_eq!(fluid.state_hash(), before);
                assert_eq!(fluid.steps(), 0);
            }
        }
    }
}

#[test]
fn thermal_trial_keeps_the_exact_native_corrector_and_its_optional_residual_probe() {
    let walls = vec![
        face(1, false, NsfWallThermal::Isothermal(310.0)),
        face(
            1,
            true,
            NsfWallThermal::Film {
                coefficient_w_m2_k: 500.0,
                ambient_temperature_k: 290.0,
            },
        ),
        face(2, true, NsfWallThermal::HeatFlux(125.0)),
    ];
    for controlled in [false, true] {
        let mut fluid = grid();
        let mut direct = grid();
        if controlled {
            let control = NsfCouplingControl {
                maximum_endpoint_evaluations: 32,
                force_relative_tolerance: 1e-8,
                velocity_relative_tolerance: 1e-8,
                caloric_relative_tolerance: 1e-8,
                energy_relative_tolerance: 1e-8,
            };
            fluid.set_coupling_control(Some(control)).unwrap();
            direct.set_coupling_control(Some(control)).unwrap();
        }
        direct.set_thermal_walls(walls.clone()).unwrap();
        let expected = direct.step().unwrap();
        let trial = fluid
            .preview_thermal_wall_step(walls.clone(), 64, || false)
            .unwrap();
        assert_eq!(trial.grid().state_hash(), direct.state_hash());
        assert_eq!(trial.energy_ledger().wall_heat_j, expected.wall_heat_j);
        assert_eq!(trial.energy_ledger().internal_j, expected.internal_j);
        let heat: f64 = trial
            .wall_exchange()
            .iter()
            .map(|face| face.heat_into_fluid_j)
            .sum();
        near(heat, expected.wall_heat_j, 1e-22);
        trial.accept();
        assert_eq!(fluid.state_hash(), direct.state_hash());
        // Reusing held controls still gives a complete old/new native heat
        // pair: the recorded trial explicitly refreshes its starting history.
        let trial = fluid
            .preview_thermal_wall_step(walls.clone(), 64, || false)
            .unwrap();
        assert_eq!(trial.wall_exchange().len(), walls.len());
        near(
            trial
                .wall_exchange()
                .iter()
                .map(|face| face.heat_into_fluid_j)
                .sum(),
            trial.energy_ledger().wall_heat_j,
            1e-22,
        );
        trial.accept();
        assert_eq!(fluid.steps(), 2);
    }
}

#[test]
fn thermal_trials_refuse_bad_faces_work_and_noncopyable_sources_without_publishing() {
    let mut fluid = grid();
    let before = fluid.state_hash();
    let wall = face(1, false, NsfWallThermal::HeatFlux(250.0));
    assert!(
        fluid
            .preview_thermal_wall_step(vec![wall], 63, || false)
            .is_err()
    );
    let mut bad = wall;
    bad.cell = [1; 3];
    assert!(
        fluid
            .preview_thermal_wall_step(vec![bad], 64, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        fluid
            .preview_thermal_wall_step(vec![wall], 64, || {
                polls += 1;
                polls >= 12
            })
            .is_err()
    );
    assert!(polls >= 12);
    assert_eq!(fluid.state_hash(), before);
    assert!(fluid.thermal_walls().is_empty());
    assert_eq!(fluid.steps(), 0);
    let mut source = NsfGrid3::<D3q27>::new_with_sources(
        spec(),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(|_, _| (0.0, [0.0; 3], 0.0))),
    )
    .unwrap();
    let before = source.state_hash();
    assert!(
        source
            .preview_thermal_wall_step(vec![wall], 64, || false)
            .is_err()
    );
    assert_eq!(source.state_hash(), before);
}
