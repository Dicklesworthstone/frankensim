//! G0/G1 physical wall-work and heat checks on the original NSF operator.

use super::*;
use crate::d3q27::D3q27;

fn spec(axis: usize, wall: WallThermal) -> NsfSpec {
    let mut periodic = [true; 3];
    periodic[axis] = false;
    NsfSpec {
        dims: [4; 3],
        periodic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: MAX_QUALIFIED_MACH,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012, N2 at 293 K",
        maximum_cells: 64,
        wall,
        workers: 1,
        transport_points: 17,
    }
}

#[test]
fn g0_uniform_pressure_work_telescopes_across_wall_and_solid_stencils() {
    for axis in 0..3 {
        for partitioned in [false, true] {
            let mut setup = spec(axis, WallThermal::Adiabatic);
            setup.dims[axis] = 11;
            setup.maximum_cells = setup.dims.iter().product();
            let mut grid = NsfGrid3::<D3q27>::new(
                setup,
                |x| partitioned && x[axis] == 5,
                |_| (1.0, [0.0; 3], 293.15),
            )
            .unwrap();
            let rho = vec![1.0; grid.n];
            let theta = vec![THETA0; grid.n];
            let u: Vec<_> = (0..grid.n)
                .map(|i| {
                    let x = (grid.at(i)[axis] + 1) as f64;
                    let mut velocity = [0.0; 3];
                    velocity[axis] = 1e-5 * x * x * x;
                    velocity
                })
                .collect();
            let work = grid.pressure_work_density(&rho, &u, &theta, None);
            let regions = if partitioned { 2 } else { 1 };
            for region in 0..regions {
                let cells: Vec<_> = (0..grid.n)
                    .filter(|&i| {
                        !grid.solid[i]
                            && (!partitioned || usize::from(grid.at(i)[axis] > 5) == region)
                    })
                    .collect();
                let total: f64 = cells.iter().map(|&i| work[i]).sum();
                let magnitude: f64 = cells.iter().map(|&i| work[i].abs()).sum();
                assert!(magnitude > 0.0);
                assert!(
                    total.abs() < 64.0 * f64::EPSILON * magnitude,
                    "axis={axis}, partitioned={partitioned}, region={region}: {total:e}"
                );
            }
            assert!((0..grid.n).all(|i| !grid.solid[i] || work[i] == 0.0));
            // Cell 1 shares its left two-cell face and its right four-cell
            // face with their neighbors. For u=1e-5*(x+1)^3, those fluxes
            // are 4.5e-5 and 15e-5. A per-cell switch gives 13e-5 instead.
            let mut cell = [0; 3];
            cell[axis] = 1;
            assert!((work[grid.index(cell)] - THETA0 * 10.5e-5).abs() < 1e-18);
            grid.spec.workers = 3;
            assert_eq!(work, grid.pressure_work_density(&rho, &u, &theta, None));
        }
    }
}

#[test]
fn g1_isolated_boundary_cells_retain_reservoir_compression_work() {
    for axis in 0..3 {
        let mut setup = spec(axis, WallThermal::Adiabatic);
        setup.dims = [2; 3];
        setup.dims[axis] = 3;
        let mut grid = NsfGrid3::<D3q27>::new(
            setup,
            |x| x[axis] == 1,
            |_| (1.0, [0.0; 3], 293.15),
        )
        .unwrap();
        // Two admitted reservoirs separated by a solid layer: each wet
        // boundary cell has no fluid neighbor along its open axis.
        let ghost_speeds = [0.014, -0.006];
        for (side, &speed) in ghost_speeds.iter().enumerate() {
            let mut velocity_m_s = [0.0; 3];
            velocity_m_s[axis] = speed * grid.law.dx_m / grid.law.dt_s;
            grid.set_open_face(
                axis,
                side == 1,
                Some(NsfOpenState {
                    density_kg_m3: grid.law.rho_ref_kg_m3,
                    velocity_m_s,
                    temperature_k: 293.15,
                }),
            )
            .unwrap();
        }
        let rho = vec![1.0; grid.n];
        let theta = vec![THETA0; grid.n];
        let mut velocity = [0.0; 3];
        velocity[axis] = 0.004;
        let work = grid.pressure_work_density(&rho, &vec![velocity; grid.n], &theta, None);
        for (i, &actual) in work.iter().enumerate() {
            let coordinate = grid.at(i)[axis];
            let expected = match coordinate {
                0 => -THETA0 * 0.5 * (velocity[axis] + ghost_speeds[0]),
                2 => THETA0 * 0.5 * (velocity[axis] + ghost_speeds[1]),
                _ => 0.0,
            };
            assert!(
                (actual - expected).abs() < 1e-17,
                "axis={axis}, cell={i}"
            );
        }
        assert!(work.iter().sum::<f64>() < 0.0, "both chambers compress");
    }
}

#[test]
fn g1_pressure_work_integrates_the_actual_characteristic_face_velocity() {
    for axis in 0..3 {
        let mut setup = spec(axis, WallThermal::Adiabatic);
        setup.dims = [2; 3];
        setup.dims[axis] = 7;
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: setup.pressure_pa,
            reference_temperature_k: setup.temperature_k,
            reference_velocity_m_s: [0.0; 3],
            incoming_characteristic_m_s: 1.0,
        };
        let grid = NsfGrid3::<D3q27>::new(
            setup,
            |_| false,
            |x| {
                let mut velocity = [0.0; 3];
                velocity[axis] = 0.0001 * ((x[axis] + 1) as f64).powi(2);
                (1.0, velocity, 293.15)
            },
        )
        .unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis,
            faces: [face; 2],
            maximum_relative_perturbation: 0.05,
        })
        .unwrap();
        let (rho, momentum) = grid.macro_fields(&grid.f);
        let velocity = grid.velocity_field(&rho, &momentum, &grid.force);
        let work = grid.pressure_work_density(
            &vec![1.0; grid.n],
            &velocity,
            &vec![THETA0; grid.n],
            grid.open_faces.as_ref(),
        );
        // Independent public SI face traces, including their outward signs.
        // A characteristic trace already lives on the face: it must not be
        // averaged with the cell as if it were a prescribed ghost-center state.
        let expected = THETA0
            * [false, true]
                .into_iter()
                .flat_map(|high| grid.open_boundary_traces(high).unwrap())
                .map(|trace| {
                    trace.outward_normal_velocity_m_s * grid.law.dt_s / grid.law.dx_m
                })
                .sum::<f64>();
        let actual: f64 = work.iter().sum();
        let scale: f64 = work.iter().map(|value| value.abs()).sum();
        assert!((actual - expected).abs() < 64.0 * f64::EPSILON * scale);
    }
}

#[test]
fn g1_periodic_compression_work_retains_fourth_order_spatial_convergence() {
    for axis in 0..3 {
        let mut errors = Vec::new();
        for n in [16, 32, 64] {
            let mut setup = spec(axis, WallThermal::Adiabatic);
            setup.periodic = [true; 3];
            setup.dims = [1; 3];
            setup.dims[axis] = n;
            let grid = NsfGrid3::<D3q27>::new(
                setup,
                |_| false,
                |_| (1.0, [0.0; 3], 293.15),
            )
            .unwrap();
            let wave_number = std::f64::consts::TAU / n as f64;
            let amplitude = 0.01;
            let phase: Vec<_> = (0..grid.n)
                .map(|i| wave_number * (grid.at(i)[axis] as f64 + 0.5))
                .collect();
            let rho: Vec<_> = phase.iter().map(|&x| 1.0 + 0.02 * x.cos()).collect();
            let velocity: Vec<_> = phase
                .iter()
                .map(|&x| {
                    let mut velocity = [0.0; 3];
                    velocity[axis] = amplitude * x.sin();
                    velocity
                })
                .collect();
            let work =
                grid.pressure_work_density(&rho, &velocity, &vec![THETA0; grid.n], None);
            let scale = THETA0 * amplitude * wave_number;
            let error = work
                .iter()
                .enumerate()
                .map(|(i, &actual)| {
                    let exact = rho[i] * scale * phase[i].cos();
                    (actual - exact).abs() / scale
                })
                .fold(0.0_f64, f64::max);
            errors.push(error);
            assert!(work.iter().any(|&value| value < 0.0));
            assert!(work.iter().any(|&value| value > 0.0));
            assert!(work.iter().sum::<f64>() > 0.0);
        }
        assert!(errors[2] < 4e-6, "axis={axis}: {errors:?}");
        assert!(
            errors
                .windows(2)
                .all(|e| (15.0..16.1).contains(&(e[0] / e[1]))),
            "axis={axis}: {errors:?}"
        );
    }
}

#[test]
fn g1_wall_traction_work_and_fourier_heat_have_independent_signs() {
    let stress_tensor = [0.003, -0.002, 0.004, 0.005, -0.006, 0.007];
    for axis in 0..3 {
        for high in [false, true] {
            for direction in [-1.0, 1.0] {
                for thermal in [WallThermal::Adiabatic, WallThermal::Isothermal(291.15)] {
                    let mut grid = NsfGrid3::<D3q27>::new(
                        spec(axis, thermal),
                        |_| false,
                        |_| (1.0, [0.0; 3], 293.15),
                    )
                    .unwrap();
                    let mut velocity = [0.013 * direction, -0.017 * direction, 0.019 * direction];
                    velocity[axis] = 0.0;
                    grid.set_edge_velocity(axis, high, velocity).unwrap();
                    let stress = vec![stress_tensor; grid.n];
                    let (div, total, boundary) = grid
                        .flux_divergence(
                            &grid.e,
                            &vec![1.0; grid.n],
                            &vec![[0.0; 3]; grid.n],
                            &grid.theta,
                            &stress,
                            &[],
                            None,
                        )
                        .unwrap();
                    // Independent tensor traction n.sigma paired with wall
                    // velocity. The outward total-energy flux has minus this
                    // work; reversing the normal or velocity reverses its sign.
                    let row = match axis {
                        0 => [0.003, 0.005, -0.006],
                        1 => [0.005, -0.002, 0.007],
                        _ => [-0.006, 0.007, 0.004],
                    };
                    let normal = if high { 1.0 } else { -1.0 };
                    let work = -16.0
                        * normal
                        * row
                            .iter()
                            .zip(velocity)
                            .map(|(traction, speed)| traction * speed)
                            .sum::<f64>();
                    let heat = match thermal {
                        WallThermal::Adiabatic => 0.0,
                        WallThermal::Isothermal(kelvin) => {
                            let wall_theta = grid.law.theta_of_kelvin(kelvin);
                            let conductivity = grid
                                .spec
                                .gas
                                .face_conductivity(
                                    &grid.law,
                                    0.5 * (THETA0 + wall_theta),
                                    grid.law.pressure_pa(1.0, THETA0),
                                )
                                .unwrap();
                            32.0 * conductivity * (THETA0 - wall_theta) / 0.5
                        }
                    };
                    assert!((boundary.wall_work_out - work).abs() < 1e-16);
                    assert!((boundary.wall_heat_out - heat).abs() < 1e-16);
                    assert!((total - heat - work).abs() < 1e-16);
                    assert!((div.iter().sum::<f64>() - total).abs() < 1e-16);
                }
            }
        }
    }
}

#[test]
fn g0_equal_opposite_heat_and_work_do_not_erase_either_transfer() {
    let mut grid = NsfGrid3::<D3q27>::new(
        spec(1, WallThermal::Isothermal(291.15)),
        |_| false,
        |_| (1.0, [0.0; 3], 293.15),
    )
    .unwrap();
    let speed = 0.01;
    grid.set_edge_velocity(1, true, [speed, 0.0, 0.0]).unwrap();
    let wall_theta = grid.law.theta_of_kelvin(291.15);
    let conductivity = grid
        .spec
        .gas
        .face_conductivity(
            &grid.law,
            0.5 * (THETA0 + wall_theta),
            grid.law.pressure_pa(1.0, THETA0),
        )
        .unwrap();
    let heat_per_face = conductivity * (THETA0 - wall_theta) / 0.5;
    let stress = vec![[0.0, 0.0, 0.0, 2.0 * heat_per_face / speed, 0.0, 0.0]; grid.n];
    let (_, total, boundary) = grid
        .flux_divergence(
            &grid.e,
            &vec![1.0; grid.n],
            &vec![[0.0; 3]; grid.n],
            &grid.theta,
            &stress,
            &[],
            None,
        )
        .unwrap();
    let expected_heat = 32.0 * heat_per_face;
    assert!(expected_heat > 0.0);
    assert!(total.abs() < 1e-15 * expected_heat);
    assert!((boundary.wall_heat_out - expected_heat).abs() < 1e-15 * expected_heat);
    assert!((boundary.wall_work_out + expected_heat).abs() < 1e-15 * expected_heat);
}

#[test]
fn g1_adiabatic_couette_wall_work_enters_energy_and_viscous_heat() {
    let mut setup = spec(1, WallThermal::Adiabatic);
    setup.dims = [2, 8, 2];
    let speed = 0.01;
    let mut grid = NsfGrid3::<D3q27>::new(
        setup,
        |_| false,
        |x| (1.0, [speed * (x[1] as f64 + 0.5) / 8.0, 0.0, 0.0], 293.15),
    )
    .unwrap();
    grid.set_edge_velocity(1, true, [speed, 0.0, 0.0]).unwrap();
    let initial_total = grid.total_energy_j();
    let initial_internal = initial_total - grid.kinetic_energy_j();
    let mut summed_work = 0.0;
    for _ in 0..64 {
        let before = grid.total_energy_j();
        let ledger = grid.step().unwrap();
        assert_eq!(
            ledger.wall_heat_j, 0.0,
            "an adiabatic wall supplies no heat"
        );
        assert_eq!(ledger.unresolved_internal_exchange_j, 0.0);
        assert!((grid.total_energy_j() - before - ledger.wall_work_j).abs() < 3e-14 * before);
        summed_work += ledger.wall_work_j;
    }
    let ledger = grid.cumulative_ledger();
    assert!(summed_work > 0.0);
    assert_eq!(summed_work, ledger.wall_work_j);
    assert_eq!(ledger.wall_heat_j, 0.0);
    assert!(ledger.internal_j > initial_internal);
    assert!(ledger.viscous_heating_j > 0.0);
    assert!((grid.total_energy_j() - initial_total - summed_work).abs() < 3e-13 * initial_total);
    assert!(
        (ledger.internal_j - initial_internal - ledger.viscous_heating_j - ledger.pressure_work_j)
            .abs()
            < 3e-13 * initial_total
    );

    // Refusing a new wall law neither publishes a transfer nor corrupts the
    // next accepted quadrature history or the last accepted physical state.
    let before = grid.state_hash();
    assert!(grid.set_edge_velocity(1, true, [0.0, 0.01, 0.0]).is_err());
    assert_eq!(grid.state_hash(), before);
    assert_eq!(grid.cumulative_ledger(), ledger);
}

#[test]
fn g0_characteristic_intervals_retain_transverse_wall_work() {
    let make = || {
        let mut setup = spec(1, WallThermal::Adiabatic);
        setup.periodic[0] = false;
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: setup.pressure_pa,
            reference_temperature_k: setup.temperature_k,
            reference_velocity_m_s: [0.0; 3],
            incoming_characteristic_m_s: 0.0,
        };
        let mut grid = NsfGrid3::<D3q27>::new(
            setup,
            |_| false,
            |x| (1.0, [0.0, 0.0, 0.01 * (x[1] as f64 + 0.5) / 4.0], 293.15),
        )
        .unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0,
            faces: [face; 2],
            maximum_relative_perturbation: 0.05,
        })
        .unwrap();
        grid.set_edge_velocity(1, true, [0.0, 0.0, 0.01]).unwrap();
        grid
    };
    let mut direct = make();
    let mut previewed = make();
    direct.set_incoming_characteristics([0.0; 2]).unwrap();
    let before_hash = previewed.state_hash();
    let before_ledger = previewed.cumulative_ledger();
    let mut work = 0.0;
    for _ in 0..3 {
        work += direct.step().unwrap().wall_work_j;
    }
    assert!(work > 0.0);
    {
        let trial = previewed
            .preview_characteristic_interval([0.0; 2], 3, 192, || false)
            .unwrap();
        let ledger = trial.energy_ledger();
        assert_eq!(ledger.wall_work_j, work);
        assert_eq!(ledger.wall_heat_j, 0.0);
        assert_eq!(trial.grid().state_hash(), direct.state_hash());
        assert!(ledger.residual_rel < 3e-14);
    }
    assert_eq!(previewed.state_hash(), before_hash);
    assert_eq!(previewed.cumulative_ledger(), before_ledger);
    previewed
        .preview_characteristic_interval([0.0; 2], 3, 192, || false)
        .unwrap()
        .accept();
    assert_eq!(previewed.cumulative_ledger(), direct.cumulative_ledger());
}
