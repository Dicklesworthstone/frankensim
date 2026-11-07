use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfSpec, THETA0, WallThermal};
use fs_material::gas::GasSpec;

fn spec(periodic: [bool; 3], workers: usize) -> NsfSpec {
    NsfSpec {
        dims: [8, 4, 4],
        periodic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0,
        pressure_pa: 100_000.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 ratio; spatial thermal wall fixture",
        maximum_cells: 128,
        wall: WallThermal::Adiabatic,
        workers,
        transport_points: 17,
    }
}

fn grid(periodic: [bool; 3], workers: usize) -> NsfGrid3<D3q27> {
    NsfGrid3::new(
        spec(periodic, workers),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap()
}

fn face(cell: [usize; 3], axis: usize, high: bool, thermal: NsfWallThermal) -> NsfThermalWallFace {
    NsfThermalWallFace {
        cell,
        axis,
        high,
        thermal,
    }
}

fn resting_flux(grid: &NsfGrid3<D3q27>) -> (Vec<f64>, f64, super::super::BoundaryFlux) {
    grid.flux_divergence(
        &grid.e,
        &vec![1.0; grid.n],
        &vec![[0.0; 3]; grid.n],
        &vec![THETA0; grid.n],
        &vec![[0.0; 6]; grid.n],
        &[],
        grid.open_faces.as_ref(),
    )
    .unwrap()
}

#[test]
fn prescribed_heat_enters_each_orientation_and_interior_solid_face_once() {
    for axis in 0..3 {
        for high in [false, true] {
            for q in [-250.0, 250.0] {
                let mut grid = grid([false; 3], 1);
                let mut cell = [1; 3];
                cell[axis] = if high { grid.spec.dims[axis] - 1 } else { 0 };
                grid.set_thermal_walls(vec![face(cell, axis, high, NsfWallThermal::HeatFlux(q))])
                    .unwrap();
                let (div, total, boundary) = resting_flux(&grid);
                let expected =
                    -q * grid.law.dt_s * grid.law.dx_m.powi(2) / grid.law.energy_unit_j();
                assert!((div[grid.index(cell)] - expected).abs() < 1e-20);
                assert!((total - expected).abs() < 1e-20);
                assert_eq!(total, boundary.wall_heat_out);
                assert_eq!(boundary.wall_work_out, 0.0);
                assert_eq!(div.iter().sum::<f64>(), total);
            }
        }
    }
    let mut grid = NsfGrid3::<D3q27>::new(
        spec([true; 3], 1),
        |x| x == [2, 2, 2],
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap();
    grid.set_thermal_walls(vec![face(
        [1, 2, 2],
        0,
        true,
        NsfWallThermal::HeatFlux(400.0),
    )])
    .unwrap();
    let (div, total, boundary) = resting_flux(&grid);
    let expected = -400.0 / grid.law.heat_flux_unit_w_m2();
    assert_eq!(div[grid.index([2, 2, 2])], 0.0);
    assert_eq!(div[grid.index([1, 2, 2])], expected);
    assert_eq!(total, expected);
    assert_eq!(boundary.wall_heat_out, expected);
}

#[test]
fn film_matches_independent_series_resistance_and_actual_interface_conductivity() {
    let mut constant = grid([true, false, true], 1);
    let k_lu = constant.law.kappa_lattice(THETA0);
    constant.law.kappa.fill(k_lu);
    let k_si =
        k_lu * constant.law.gas_constant * constant.law.rho_ref_kg_m3 * constant.law.dx_m.powi(2)
            / constant.law.dt_s;
    let i = constant.index([2, 0, 2]);
    for ambient in [280.0, 320.0] {
        for h in [0.0, 1e-12, 50.0, 1e18] {
            let flux = constant
                .thermal_wall_flux(
                    i,
                    1,
                    -1,
                    NsfWallThermal::Film {
                        coefficient_w_m2_k: h,
                        ambient_temperature_k: ambient,
                    },
                    1.0,
                    THETA0,
                )
                .unwrap()
                * constant.law.heat_flux_unit_w_m2();
            let expected = if h == 0.0 {
                0.0
            } else {
                (ambient - 300.0) / (constant.law.dx_m / (2.0 * k_si) + 1.0 / h)
            };
            assert!(
                (flux - expected).abs() <= 2e-13 * expected.abs() + 1e-22,
                "h={h}, ambient={ambient}, flux={flux}, expected={expected}"
            );
        }
    }
    let variable = grid([true, false, true], 1);
    let h = 50.0;
    let ambient = 320.0;
    let q = variable
        .thermal_wall_flux(
            i,
            1,
            -1,
            NsfWallThermal::Film {
                coefficient_w_m2_k: h,
                ambient_temperature_k: ambient,
            },
            1.0,
            THETA0,
        )
        .unwrap()
        * variable.law.heat_flux_unit_w_m2();
    // Recover Tw from the independent ambient-film equation, then test the
    // native molecular half-cell equation at its actual fluid/interface mean.
    let interface = ambient - q / h;
    assert!(interface > 300.0 && interface < ambient);
    let conversion =
        variable.law.gas_constant * variable.law.rho_ref_kg_m3 * variable.law.dx_m.powi(2)
            / variable.law.dt_s;
    let k = variable
        .law
        .kappa_lattice(variable.law.theta_of_kelvin(0.5 * (300.0 + interface)))
        * conversion;
    let molecular_q = 2.0 * k * (interface - 300.0) / variable.law.dx_m;
    assert!((q - molecular_q).abs() < 2e-10 * q.abs());
    let wrong_k = variable
        .law
        .kappa_lattice(variable.law.theta_of_kelvin(310.0))
        * conversion;
    let bulk_ambient_guess = (ambient - 300.0) / (variable.law.dx_m / (2.0 * wrong_k) + 1.0 / h);
    assert!((q - bulk_ambient_guess).abs() > 1e-5 * q.abs());
}

#[test]
fn heterogeneous_walls_heat_and_cool_original_flow_and_close_energy() {
    let walls: Vec<_> = (0..8)
        .flat_map(|x| {
            (0..4).flat_map(move |z| {
                [
                    face(
                        [x, 0, z],
                        1,
                        false,
                        if x < 4 {
                            NsfWallThermal::Isothermal(310.0)
                        } else {
                            NsfWallThermal::HeatFlux(500.0)
                        },
                    ),
                    face(
                        [x, 3, z],
                        1,
                        true,
                        NsfWallThermal::Film {
                            coefficient_w_m2_k: 50.0,
                            ambient_temperature_k: 290.0,
                        },
                    ),
                ]
            })
        })
        .collect();
    let mut serial = grid([true, false, true], 1)
        .with_thermal_walls(walls.clone())
        .unwrap();
    let mut pooled = grid([true, false, true], 3)
        .with_thermal_walls(walls)
        .unwrap();
    // A tangential lid supplies genuine mechanical work alongside heterogeneous
    // heat, on the existing flow/stress/energy owner.
    for grid in [&mut serial, &mut pooled] {
        grid.set_edge_velocity(1, true, [0.01, 0.0, 0.0]).unwrap();
    }
    let initial = serial.total_energy_j();
    let mut supplied = 0.0;
    for _ in 0..8 {
        let ledger = serial.step().unwrap();
        assert_eq!(ledger, pooled.step().unwrap());
        assert_eq!(serial.state_hash(), pooled.state_hash());
        supplied += ledger.wall_heat_j + ledger.wall_work_j;
        assert!(ledger.residual_rel < 3e-12);
        assert!(ledger.wall_heat_j > 0.0);
    }
    assert!((serial.total_energy_j() - initial - supplied).abs() < 3e-12 * initial);
    assert!(serial.theta[serial.index([1, 0, 1])] > THETA0);
    assert!(serial.theta[serial.index([1, 3, 1])] < THETA0);
    assert!(serial.cumulative_ledger().wall_work_j > 0.0);
}

#[test]
fn empty_and_uniform_overrides_keep_original_wall_arithmetic() {
    let mut input = spec([true, false, true], 1);
    input.wall = WallThermal::Isothermal(305.0);
    let make =
        || NsfGrid3::<D3q27>::new(input.clone(), |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap();
    let mut original = make();
    let mut empty = make().with_thermal_walls(Vec::new()).unwrap();
    let walls = (0..8)
        .flat_map(|x| {
            (0..4).flat_map(move |z| {
                [
                    face([x, 0, z], 1, false, NsfWallThermal::Isothermal(305.0)),
                    face([x, 3, z], 1, true, NsfWallThermal::Isothermal(305.0)),
                ]
            })
        })
        .collect();
    let mut local = make().with_thermal_walls(walls).unwrap();
    for _ in 0..3 {
        let expected = original.step().unwrap();
        assert_eq!(empty.step().unwrap(), expected);
        assert_eq!(local.step().unwrap(), expected);
        assert_eq!(empty.state_hash(), original.state_hash());
        assert_eq!(local.f, original.f);
        assert_eq!(local.e, original.e);
        assert_eq!(local.div, original.div);
    }
    local.set_thermal_walls(Vec::new()).unwrap();
    assert_eq!(local.state_hash(), original.state_hash());
    // Rebinding at an evolved accepted state refreshes the native history,
    // even though a fixed-pass predictor history can differ from that endpoint.
    assert_eq!(local.f, original.f);
    assert_eq!(local.e, original.e);
}

#[test]
fn malformed_faces_and_final_cancellation_preserve_controls_fields_and_history() {
    let mut grid = grid([true, false, true], 1);
    let accepted = face([1, 0, 1], 1, false, NsfWallThermal::HeatFlux(75.0));
    grid.set_thermal_walls(vec![accepted]).unwrap();
    let before = (
        grid.state_hash(),
        grid.div.clone(),
        grid.wall_out,
        grid.open_out.wall_heat_out,
        grid.e.clone(),
        grid.theta.clone(),
        grid.force.clone(),
        grid.cumulative_ledger(),
    );
    for bad in [
        face([1, 1, 1], 1, false, NsfWallThermal::Adiabatic),
        face([0, 1, 1], 0, false, NsfWallThermal::Adiabatic),
        face([1, 0, 1], 3, false, NsfWallThermal::Adiabatic),
        face([8, 0, 1], 1, false, NsfWallThermal::Adiabatic),
        face([1, 0, 1], 1, false, NsfWallThermal::HeatFlux(f64::NAN)),
        face([1, 0, 1], 1, false, NsfWallThermal::Isothermal(270.0)),
        face(
            [1, 0, 1],
            1,
            false,
            NsfWallThermal::Film {
                coefficient_w_m2_k: -1.0,
                ambient_temperature_k: 300.0,
            },
        ),
    ] {
        assert!(grid.set_thermal_walls(vec![bad]).is_err());
        assert_eq!(grid.thermal_walls(), &[accepted]);
        assert_eq!(grid.state_hash(), before.0);
        assert_eq!(grid.div, before.1);
    }
    assert!(grid.set_thermal_walls(vec![accepted, accepted]).is_err());
    let mut polls = 0;
    assert_eq!(
        grid.set_thermal_walls_cancellable(
            vec![face(
                [2, 3, 2],
                1,
                true,
                NsfWallThermal::Film {
                    coefficient_w_m2_k: 100.0,
                    ambient_temperature_k: 280.0
                }
            )],
            || {
                polls += 1;
                polls == 3
            }
        ),
        Err(NsfError::Input("thermal wall update cancelled"))
    );
    assert_eq!(polls, 3);
    assert_eq!(
        (
            grid.state_hash(),
            grid.div.clone(),
            grid.wall_out,
            grid.open_out.wall_heat_out,
            grid.e.clone(),
            grid.theta.clone(),
            grid.force.clone(),
            grid.cumulative_ledger()
        ),
        before
    );
    assert_eq!(grid.steps(), 0);
}

fn characteristic() -> NsfCharacteristicAxis {
    NsfCharacteristicAxis {
        axis: 0,
        maximum_relative_perturbation: 0.1,
        faces: [NsfCharacteristicBoundary {
            reference_pressure_pa: 100_000.0,
            reference_temperature_k: 300.0,
            reference_velocity_m_s: [0.0; 3],
            incoming_characteristic_m_s: 0.0,
        }; 2],
    }
}

#[test]
fn wall_opening_ownership_and_characteristic_trials_retain_local_heat() {
    let conflicting = vec![face([0, 1, 1], 0, false, NsfWallThermal::HeatFlux(100.0))];
    assert!(
        grid([false, true, true], 1)
            .with_thermal_walls(conflicting.clone())
            .unwrap()
            .with_characteristic_axis(characteristic())
            .is_err()
    );
    let mut reservoir = grid([false, true, true], 1)
        .with_thermal_walls(conflicting.clone())
        .unwrap();
    let ghost = crate::nsf3::NsfOpenState {
        density_kg_m3: reservoir.law.rho_ref_kg_m3,
        velocity_m_s: [0.0; 3],
        temperature_k: 300.0,
    };
    let reservoir_hash = reservoir.state_hash();
    assert!(reservoir.set_open_face(0, false, Some(ghost)).is_err());
    assert_eq!(reservoir.state_hash(), reservoir_hash);
    let mut open = grid([false, false, true], 1)
        .with_characteristic_axis(characteristic())
        .unwrap();
    assert!(open.set_thermal_walls(conflicting).is_err());
    let wall = face([3, 0, 2], 1, false, NsfWallThermal::HeatFlux(250.0));
    open.set_thermal_walls(vec![wall]).unwrap();
    let old_hash = open.state_hash();
    let trial = open
        .preview_characteristic_interval([0.0; 2], 2, 256, || false)
        .unwrap();
    assert_eq!(trial.grid().thermal_walls(), &[wall]);
    assert!(trial.energy_ledger().wall_heat_j > 0.0);
    drop(trial);
    assert_eq!(open.state_hash(), old_hash);
    assert_eq!(open.steps(), 0);
    open.preview_characteristic_interval([0.0; 2], 2, 256, || false)
        .unwrap()
        .accept();
    assert_eq!(open.thermal_walls(), &[wall]);
    assert_eq!(open.steps(), 2);
}

#[test]
fn film_uses_actual_pressure_and_retains_named_source_domain_refusal() {
    use crate::nsf3::NamedMixtureNsfGas;
    use fs_material::gas_species::GasSpecies;
    use fs_qty::semantic::{Composition, CompositionBasis};
    let gas = NamedMixtureNsfGas::new(
        vec![GasSpecies::Nitrogen, GasSpecies::Oxygen],
        Composition::new(CompositionBasis::MoleFraction, vec![0.79, 0.21]).unwrap(),
        300.0,
        100_000.0,
    )
    .unwrap();
    let mut input = spec([true, false, true], 1).map_gas(|_| gas);
    input.temperature_window_k = [280.0, 320.0];
    let grid = NsfGrid3::<D3q27>::new(input, |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap();
    let film = NsfWallThermal::Film {
        coefficient_w_m2_k: 50.0,
        ambient_temperature_k: 310.0,
    };
    assert!(
        grid.thermal_wall_flux(grid.index([2, 0, 2]), 1, -1, film, 1.0, THETA0)
            .is_ok()
    );
    assert!(matches!(
        grid.thermal_wall_flux(grid.index([2, 0, 2]), 1, -1, film, 3.0, THETA0),
        Err(NsfError::Gas(_))
    ));
}
