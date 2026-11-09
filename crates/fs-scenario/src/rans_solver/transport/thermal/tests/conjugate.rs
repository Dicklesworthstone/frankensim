//! Physical-unit coupled-system oracles. The reference assembles and pivots
//! the entire solid/fluid matrix without production elimination or impedance.

use super::*;

fn layers() -> [RansWallLayer; 3] {
    [
        RansWallLayer {
            thickness_m: 0.05,
            volumetric_heat_capacity_j_m3_k: 15.0,
            conductivity_w_m_k: 0.02,
        },
        RansWallLayer {
            thickness_m: 0.2,
            volumetric_heat_capacity_j_m3_k: 8.0,
            conductivity_w_m_k: 0.9,
        },
        RansWallLayer {
            thickness_m: 0.1,
            volumetric_heat_capacity_j_m3_k: 40.0,
            conductivity_w_m_k: 0.08,
        },
    ]
}

fn wall(temperatures: &[f64; 3]) -> RansConjugateWall {
    RansConjugateWall::try_new(&layers(), temperatures, 3).unwrap()
}

fn solid_energy(wall: &RansConjugateWall) -> f64 {
    wall.layers()
        .iter()
        .zip(wall.temperatures_k())
        .map(|(p, t)| p.thickness_m * p.volumetric_heat_capacity_j_m3_k * t)
        .sum()
}

#[allow(clippy::too_many_arguments)]
fn dense_joint(
    solver: &RansSolver,
    grid: &RansChannelGrid,
    state: &RansFieldState,
    wall: &RansConjugateWall,
    nut: &[f64],
    fluid_sources: &[f64],
    wall_sources: &[f64],
    exterior_flux: f64,
    dt: f64,
) -> Vec<f64> {
    let m = wall.layers().len();
    let n = state.temp_k.len();
    let mut a = vec![vec![0.0; n + m]; n + m];
    let mut b = vec![0.0; n + m];
    for (i, layer) in wall.layers().iter().enumerate() {
        let capacity = layer.thickness_m * layer.volumetric_heat_capacity_j_m3_k;
        a[i][i] += capacity / dt;
        b[i] = capacity * wall.temperatures_k()[i] / dt + layer.thickness_m * wall_sources[i];
        let resistance = layer.thickness_m / (2.0 * layer.conductivity_w_m_k)
            + wall.layers().get(i + 1).map_or(0.0, |next| {
                next.thickness_m / (2.0 * next.conductivity_w_m_k)
            });
        let g = 1.0 / resistance;
        a[i][i] += g;
        a[i + 1][i + 1] += g;
        a[i][i + 1] -= g;
        a[i + 1][i] -= g;
    }
    b[0] += exterior_flux;
    let mut volume = vec![0.0; n];
    for i in 0..n - 1 {
        let width = grid.y_nodes[i + 1] - grid.y_nodes[i];
        volume[i] += width / 2.0;
        volume[i + 1] += width / 2.0;
        let g = (solver.k_thermal
            + solver.rho * solver.cp * (nut[i] + nut[i + 1]) / (2.0 * solver.pr_turbulent))
            / width;
        a[m + i][m + i] += g;
        a[m + i + 1][m + i + 1] += g;
        a[m + i][m + i + 1] -= g;
        a[m + i + 1][m + i] -= g;
    }
    for i in 0..n {
        let storage = solver.rho * solver.cp * volume[i] / dt;
        a[m + i][m + i] += storage;
        b[m + i] = storage * state.temp_k[i] + volume[i] * fluid_sources[i];
    }
    for column in 0..n + m {
        let pivot = (column..n + m)
            .max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))
            .unwrap();
        a.swap(column, pivot);
        b.swap(column, pivot);
        let diagonal = a[column][column];
        for value in &mut a[column][column..] {
            *value /= diagonal;
        }
        b[column] /= diagonal;
        for i in 0..n + m {
            if i == column {
                continue;
            }
            let factor = a[i][column];
            for j in column..n + m {
                a[i][j] -= factor * a[column][j];
            }
            b[i] -= factor * b[column];
        }
    }
    b
}

#[test]
fn g1_joint_wall_fluid_endpoints_match_the_independent_full_system() {
    let grid = grid();
    let solver = solver();
    let initial = seeded(&grid);
    let initial_wall = wall(&[350.0, 320.0, 270.0]);
    let fluid_heat = [0.2, -0.1, 0.3, 0.0, -0.05, 0.1, 0.25];
    let wall_heat = [0.6, -0.2, 0.1];
    for coupling in [RansHeatCoupling::Passive, RansHeatCoupling::Dissipative] {
        for dt in [0.001, 0.07, 1.3, 12.0] {
            let mut mechanical = initial.clone();
            let mut native_work = RansTransportWorkspace::try_new(&grid, 7).unwrap();
            solver
                .advance_launder_sharma(&grid, &mut mechanical, -0.8, dt, &mut native_work, || {
                    false
                })
                .unwrap();
            let mut sources = fluid_heat.to_vec();
            if coupling == RansHeatCoupling::Dissipative {
                let physical =
                    dissipative::physical_heating(&solver, &grid, &initial, &mechanical).0;
                for (q, loss) in sources.iter_mut().zip(physical) {
                    *q += loss;
                }
            }
            let expected = dense_joint(
                &solver,
                &grid,
                &initial,
                &initial_wall,
                &mechanical.nu_t,
                &sources,
                &wall_heat,
                0.7,
                dt,
            );
            let mut state = initial.clone();
            let mut wall = initial_wall.clone();
            let mut work = RansConjugateWorkspace::try_new(&grid, &wall, 7, 3).unwrap();
            let report = solver
                .advance_launder_sharma_with_conjugate_wall(
                    &grid,
                    &mut state,
                    &mut wall,
                    -0.8,
                    dt,
                    RansConjugateHeatLoad {
                        exterior_heat_flux_w_m2: 0.7,
                        wall_volumetric_heating_w_m3: Some(&wall_heat),
                        fluid_volumetric_heating_w_m3: Some(&fluid_heat),
                        coupling,
                        ..RansConjugateHeatLoad::default()
                    },
                    &mut work,
                    || false,
                )
                .unwrap();
            for (&a, &b) in wall
                .temperatures_k()
                .iter()
                .chain(&state.temp_k)
                .zip(&expected)
            {
                close(a, b, 2e-12);
            }
            assert_eq!(&bits(&state)[..4], &bits(&mechanical)[..4]);
            let last = layers()[2];
            let flux = 2.0 * last.conductivity_w_m_k / last.thickness_m
                * (wall.temperatures_k()[2] - state.temp_k[0]);
            close(report.interface_heat_flux_w_m2, flux, 3e-15);
            let allowance = 1e-9 * (1.0 + 1.0 / dt);
            assert!(
                report.combined_thermal_balance_defect_w_m2.abs() < allowance,
                "{report:?}"
            );
            assert!(report.wall_balance_defect_w_m2.abs() < allowance);
            assert!(report.fluid.thermal_balance_defect_w_m2.abs() < allowance);
            assert!(report.interface_equation_defect_w_m2.abs() < allowance);
            close(
                report.wall_energy_before_j_m2,
                solid_energy(&initial_wall),
                3e-15,
            );
            close(report.wall_energy_after_j_m2, solid_energy(&wall), 3e-15);
            assert_eq!(work.numerical_bytes(), (9 * 7 + 4 * 3) * 8);
        }
    }
}

#[test]
fn g1_layered_steady_conduction_and_uniform_equilibrium_are_retained() {
    let grid = grid();
    let solver = solver();
    for flux in [0.0, 2.0, -1.0] {
        let mut initial = laminar(&grid, 350.0);
        for (t, &y) in initial.temp_k.iter_mut().zip(&grid.y_nodes) {
            *t -= flux / solver.k_thermal * (y - y * y / (2.0 * grid.half_height_m));
        }
        let layers = layers();
        let mut solid = [0.0; 3];
        solid[2] =
            initial.temp_k[0] + flux * layers[2].thickness_m / (2.0 * layers[2].conductivity_w_m_k);
        for i in (0..2).rev() {
            solid[i] = solid[i + 1]
                + flux
                    * (layers[i].thickness_m / (2.0 * layers[i].conductivity_w_m_k)
                        + layers[i + 1].thickness_m / (2.0 * layers[i + 1].conductivity_w_m_k));
        }
        let original_wall = wall(&solid);
        let source = vec![-flux / grid.half_height_m; 7];
        for dt in [0.0001, 0.1, 1000.0] {
            let mut state = initial.clone();
            let mut wall = original_wall.clone();
            let mut work = RansConjugateWorkspace::try_new(&grid, &wall, 7, 3).unwrap();
            let r = solver
                .advance_launder_sharma_with_conjugate_wall(
                    &grid,
                    &mut state,
                    &mut wall,
                    0.0,
                    dt,
                    RansConjugateHeatLoad {
                        exterior_heat_flux_w_m2: flux,
                        fluid_volumetric_heating_w_m3: Some(&source),
                        ..Default::default()
                    },
                    &mut work,
                    || false,
                )
                .unwrap();
            for (&a, &b) in state.temp_k.iter().zip(&initial.temp_k) {
                close(a, b, 2e-13);
            }
            for (&a, &b) in wall.temperatures_k().iter().zip(&solid) {
                close(a, b, 2e-13);
            }
            assert!((r.interface_heat_flux_w_m2 - flux).abs() < 3e-11);
            if flux == 0.0 {
                assert_eq!(bits(&state), bits(&initial));
                assert_eq!(wall, original_wall);
                assert_eq!(r.combined_thermal_balance_defect_w_m2, 0.0);
            }
        }
    }
}

#[test]
fn g1_finite_wall_equilibrates_with_the_actual_fluid_and_solid_capacities() {
    let grid = grid();
    let solver = solver();
    let mut state = laminar(&grid, 300.0);
    let mut wall = wall(&[600.0; 3]);
    let capacity_wall: f64 = layers()
        .iter()
        .map(|p| p.thickness_m * p.volumetric_heat_capacity_j_m3_k)
        .sum();
    let capacity_fluid = solver.rho * solver.cp * grid.half_height_m;
    let energy = capacity_fluid * 300.0 + capacity_wall * 600.0;
    let equilibrium = energy / (capacity_fluid + capacity_wall);
    let mut work = RansConjugateWorkspace::try_new(&grid, &wall, 7, 3).unwrap();
    for _ in 0..1000 {
        let r = solver
            .advance_launder_sharma_with_conjugate_wall(
                &grid,
                &mut state,
                &mut wall,
                0.0,
                2.0,
                RansConjugateHeatLoad::default(),
                &mut work,
                || false,
            )
            .unwrap();
        assert!(r.combined_thermal_balance_defect_w_m2.abs() < 2e-10);
        assert!(
            state
                .temp_k
                .iter()
                .chain(wall.temperatures_k())
                .all(|&t| (300.0..=600.0).contains(&t))
        );
    }
    for &t in state.temp_k.iter().chain(wall.temperatures_k()) {
        assert!((t - equilibrium).abs() < 0.002, "{t} vs {equilibrium}");
    }
    let final_energy =
        solver.rho * solver.cp * integral(&grid, &state.temp_k) + solid_energy(&wall);
    close(final_energy, energy, 5e-13);
}

#[test]
fn g4_every_coupled_cancellation_keeps_both_states_and_replays_exactly() {
    let grid = grid();
    let solver = solver();
    let initial = seeded(&grid);
    let initial_wall = wall(&[350.0, 320.0, 270.0]);
    let heat = RansConjugateHeatLoad {
        exterior_heat_flux_w_m2: 0.3,
        coupling: RansHeatCoupling::Dissipative,
        ..Default::default()
    };
    let mut work = RansConjugateWorkspace::try_new(&grid, &initial_wall, 7, 3).unwrap();
    let bytes = work.numerical_bytes();
    let mut expected = initial.clone();
    let mut expected_wall = initial_wall.clone();
    let mut polls = 0;
    let expected_report = solver
        .advance_launder_sharma_with_conjugate_wall(
            &grid,
            &mut expected,
            &mut expected_wall,
            -0.8,
            0.2,
            heat,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in 0..polls {
        let mut state = initial.clone();
        let mut wall = initial_wall.clone();
        let mut calls = 0;
        let error = solver
            .advance_launder_sharma_with_conjugate_wall(
                &grid,
                &mut state,
                &mut wall,
                -0.8,
                0.2,
                heat,
                &mut work,
                || {
                    let cancel = calls == stop;
                    calls += 1;
                    cancel
                },
            )
            .unwrap_err();
        assert_eq!(error, RansTransportError::Cancelled);
        assert_eq!(bits(&state), bits(&initial));
        assert_eq!(wall, initial_wall);
        let r = solver
            .advance_launder_sharma_with_conjugate_wall(
                &grid,
                &mut state,
                &mut wall,
                -0.8,
                0.2,
                heat,
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(r, expected_report);
        assert_eq!(bits(&state), bits(&expected));
        assert_eq!(wall, expected_wall);
        assert_eq!(work.numerical_bytes(), bytes);
    }
}

#[test]
fn g0_material_load_and_nonpositive_endpoint_refusals_are_atomic() {
    let grid = grid();
    let solver = solver();
    let initial = seeded(&grid);
    let initial_wall = wall(&[350.0, 320.0, 270.0]);
    assert!(RansConjugateWall::try_new(&[], &[], 3).is_err());
    assert!(RansConjugateWall::try_new(&layers(), &[300.0; 3], 2).is_err());
    assert!(RansConjugateWorkspace::try_new(&grid, &initial_wall, 7, 2).is_err());
    assert!(RansConjugateWorkspace::try_new(&grid, &initial_wall, 6, 3).is_err());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut material = layers();
        material[1].conductivity_w_m_k = bad;
        assert!(RansConjugateWall::try_new(&material, &[300.0; 3], 3).is_err());
    }
    let invalid_source = [0.0, f64::NAN, 0.0];
    let short_source = [0.0];
    let loads = [
        RansConjugateHeatLoad {
            exterior_heat_flux_w_m2: f64::INFINITY,
            ..Default::default()
        },
        RansConjugateHeatLoad {
            exterior_heat_flux_w_m2: -1e10,
            ..Default::default()
        },
        RansConjugateHeatLoad {
            wall_volumetric_heating_w_m3: Some(&invalid_source),
            ..Default::default()
        },
        RansConjugateHeatLoad {
            wall_volumetric_heating_w_m3: Some(&short_source),
            ..Default::default()
        },
        RansConjugateHeatLoad {
            fluid_volumetric_heating_w_m3: Some(&short_source),
            ..Default::default()
        },
    ];
    let mut work = RansConjugateWorkspace::try_new(&grid, &initial_wall, 7, 3).unwrap();
    for heat in loads {
        let mut state = initial.clone();
        let mut wall = initial_wall.clone();
        assert!(
            solver
                .advance_launder_sharma_with_conjugate_wall(
                    &grid,
                    &mut state,
                    &mut wall,
                    -0.8,
                    0.2,
                    heat,
                    &mut work,
                    || false,
                )
                .is_err()
        );
        assert_eq!(bits(&state), bits(&initial));
        assert_eq!(wall, initial_wall);
    }
}
