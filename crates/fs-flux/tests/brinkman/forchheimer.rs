//! Native nonlinear porous momentum, actual endpoint work and boundary reactions.

use super::*;
use fs_flux::{ForchheimerCell, ForchheimerControls, StokesBoundaryError};

fn constant_dofs(mesh: &TriMesh, velocity: [f64; 2]) -> Vec<f64> {
    mesh.edges
        .iter()
        .flat_map(|edge| {
            [
                velocity[0] * edge.normal[0] + velocity[1] * edge.normal[1],
                0.0,
            ]
        })
        .collect()
}

#[test]
fn g1_scalar_pressure_free_resistance_matches_independent_quadratic_root() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.2,
        ..FluxParams::default()
    };
    let (k, beta, force_x) = (0.1, 50.0, 100.0);
    let cells = vec![
        ForchheimerCell::try_new(
            BrinkmanCell::try_new(IntrinsicPermeability::isotropic(k).unwrap(), 1.0).unwrap(),
            beta,
        )
        .unwrap();
        mesh.tris.len()
    ];
    let traction = |_: [f64; 2], n: [f64; 2]| [-3.5 * n[0], -3.5 * n[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let solution = system
        .solve_brinkman_forchheimer_mixed(
            params,
            &cells,
            &|_| [force_x, 0.0],
            &boundaries,
            &vec![0.0; system.n_u],
            ForchheimerControls::default(),
        )
        .unwrap();
    let rate = params.nu / k;
    let exact = 2.0 * force_x / (rate + (rate * rate + 4.0 * beta * force_x).sqrt());
    let (error, divergence) = system.velocity_error(&solution.flow.field.x, &|_| [exact, 0.0]);
    assert!(
        error < 2e-9 && divergence < 2e-9,
        "error={error:e}, divergence={divergence:e}"
    );
    assert!(solution.iterations > 1 && solution.line_search_reductions > 0);
    assert!(solution.flow.field.converged(1e-9));
    assert_eq!(
        solution.flow.pressure_reference,
        StokesPressureReference::TractionDatum
    );
    for pressure in &solution.flow.pressure {
        near(*pressure, 3.5, 2e-8);
    }
    for ((inertial, darcy), area) in solution
        .cell_forchheimer_exchange
        .iter()
        .zip(&solution.cell_darcy_exchange)
        .zip(&mesh.areas)
    {
        near(
            inertial.skeleton_force[0],
            area * beta * exact * exact,
            2e-8,
        );
        near(darcy.skeleton_force[0], area * rate * exact, 2e-8);
        near(inertial.skeleton_force[1], 0.0, 2e-8);
        near(inertial.dissipation, area * beta * exact.powi(3), 2e-8);
        near(inertial.fluid_power, -inertial.dissipation, 2e-12);
        assert_eq!(inertial.skeleton_power, 0.0);
        near(
            inertial.skeleton_force[0] + darcy.skeleton_force[0],
            area * force_x,
            2e-8,
        );
    }
}

#[test]
fn g1_rotated_darcy_and_moving_forchheimer_preserve_pressure_and_frame_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let (s, c) = 0.53_f64.sin_cos();
    let rxx = params.nu * (c * c / 0.002 + s * s / 0.0005);
    let rxy = params.nu * s * c * (1.0 / 0.002 - 1.0 / 0.0005);
    let ryy = params.nu * (s * s / 0.002 + c * c / 0.0005);
    let (beta, relative) = (11.0, [-0.6_f64, -0.3_f64]);
    let speed = relative[0].hypot(relative[1]);
    let darcy = [
        rxx * relative[0] + rxy * relative[1],
        rxy * relative[0] + ryy * relative[1],
    ];
    let inertial = relative.map(|w| beta * speed * w);
    let gradient = [-2.0, 1.0];
    let force = |_| {
        [
            darcy[0] + inertial[0] + gradient[0],
            darcy[1] + inertial[1] + gradient[1],
        ]
    };
    let traction = |x: [f64; 2], n: [f64; 2]| {
        let p = 7.0 + gradient[0] * x[0] + gradient[1] * x[1];
        [-p * n[0], -p * n[1]]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    for boost in [[0.0, 0.0], [1.3, -0.8]] {
        let skeleton = [0.9 + boost[0], 0.2 + boost[1]];
        let exact = [skeleton[0] + relative[0], skeleton[1] + relative[1]];
        let cells = vec![
            ForchheimerCell::try_new(
                BrinkmanCell::try_new(tensor(0.53), 1.7)
                    .unwrap()
                    .with_skeleton_velocity(skeleton)
                    .unwrap(),
                beta,
            )
            .unwrap();
            mesh.tris.len()
        ];
        let solution = system
            .solve_brinkman_forchheimer_mixed(
                params,
                &cells,
                &force,
                &boundaries,
                &vec![0.0; system.n_u],
                ForchheimerControls::default(),
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&solution.flow.field.x, &|_| exact);
        assert!(error < 2e-9 && divergence < 2e-9);
        for (pressure, point) in solution.flow.pressure.iter().zip(&mesh.centroids) {
            near(
                *pressure,
                7.0 + gradient[0] * point[0] + gradient[1] * point[1],
                3e-8,
            );
        }
        for ((exchange, darcy_exchange), area) in solution
            .cell_forchheimer_exchange
            .iter()
            .zip(&solution.cell_darcy_exchange)
            .zip(&mesh.areas)
        {
            for axis in 0..2 {
                near(exchange.skeleton_force[axis], area * inertial[axis], 3e-9);
                near(
                    darcy_exchange.skeleton_force[axis],
                    area * darcy[axis],
                    3e-9,
                );
            }
            near(exchange.dissipation, area * beta * speed.powi(3), 3e-9);
            near(
                exchange.fluid_power,
                -area * (inertial[0] * exact[0] + inertial[1] * exact[1]),
                3e-9,
            );
            near(
                exchange.skeleton_power,
                area * (inertial[0] * skeleton[0] + inertial[1] * skeleton[1]),
                3e-9,
            );
            near(
                exchange.fluid_power + exchange.skeleton_power + exchange.dissipation,
                0.0,
                3e-13,
            );
            if boost == [0.0, 0.0] {
                assert!(exchange.fluid_power > 0.0 && exchange.skeleton_power < 0.0);
            }
        }
    }
}

#[test]
fn g1_material_heterogeneous_endpoint_retains_actual_constrained_reaction_rows() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let (density, mu) = (7.0, 0.2);
    let material = CountedMaterial {
        state: FluidTransportState::try_new(density, 1000.0, 0.5, mu).unwrap(),
        reads: Cell::new(0),
    };
    let velocity = |_: [f64; 2]| [0.3, -0.1];
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| {
            if point[0] > 0.7 {
                ForchheimerCell::try_new(BrinkmanCell::free_fluid(), 0.0).unwrap()
            } else {
                ForchheimerCell::try_new(
                    BrinkmanCell::try_new(tensor(0.53), 1.6)
                        .unwrap()
                        .with_skeleton_velocity([0.8 - 0.3 * point[1], 0.1 + 0.2 * point[0]])
                        .unwrap(),
                    12.0 + 5.0 * point[1],
                )
                .unwrap()
            }
        })
        .collect();
    let (s, c) = 0.53_f64.sin_cos();
    let rxx = c * c / 0.002 + s * s / 0.0005;
    let rxy = s * c * (1.0 / 0.002 - 1.0 / 0.0005);
    let ryy = s * s / 0.002 + c * c / 0.0005;
    let force = |point: [f64; 2]| {
        let cell = cells[mesh.containing_cell(point).unwrap()];
        if cell.brinkman().permeability().is_none() {
            return [0.0; 2];
        }
        let skeleton = cell.brinkman().skeleton_velocity();
        let u = velocity(point);
        let w = [u[0] - skeleton[0], u[1] - skeleton[1]];
        let rate = density * cell.beta_per_m() * w[0].hypot(w[1]);
        [
            mu * (rxx * w[0] + rxy * w[1]) + rate * w[0],
            mu * (rxy * w[0] + ryy * w[1]) + rate * w[1],
        ]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| MaterialStokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    let binding =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let solution = binding
        .solve_brinkman_forchheimer_mixed(
            &cells,
            &force,
            &boundaries,
            &vec![0.0; system.n_u],
            ForchheimerControls::default(),
        )
        .unwrap();
    assert_eq!(material.reads.get(), 1);
    let (error, divergence) = system.velocity_error(&solution.flow().discrete().x, &velocity);
    assert!(
        error < 2e-9 && divergence < 2e-9,
        "error={error:e}, divergence={divergence:e}"
    );
    for pressure in solution.flow().pressure_pa() {
        near(*pressure, 0.0, 3e-7);
    }
    for (cell, exchange) in solution.cell_forchheimer_exchange().iter().enumerate() {
        let coefficient = cells[cell];
        if coefficient.beta_per_m() == 0.0 {
            assert_eq!(*exchange, fs_flux::MaterialBrinkmanDragExchange::default());
            continue;
        }
        let skeleton = coefficient.brinkman().skeleton_velocity();
        let u = velocity(mesh.centroids[cell]);
        let w = [u[0] - skeleton[0], u[1] - skeleton[1]];
        let rate = mesh.areas[cell] * density * coefficient.beta_per_m() * w[0].hypot(w[1]);
        for axis in 0..2 {
            near(exchange.skeleton_force_n_per_m[axis], rate * w[axis], 3e-8);
        }
        near(
            exchange.dissipation_w_per_m,
            rate * (w[0] * w[0] + w[1] * w[1]),
            3e-8,
        );
        near(
            exchange.fluid_power_w_per_m
                + exchange.skeleton_power_w_per_m
                + exchange.dissipation_w_per_m,
            0.0,
            3e-12,
        );
    }
    near(
        solution.forchheimer_power_w_per_m(),
        solution
            .cell_forchheimer_exchange()
            .iter()
            .map(|exchange| exchange.dissipation_w_per_m)
            .sum(),
        3e-12,
    );
    let selected: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| edge)
        .collect();
    // The exact state has zero stress. A frozen final Newton resistance or a
    // missing J*u-F load on constrained rows would produce spurious wall loads.
    let load = solution
        .flow()
        .boundary_loads(&selected, [0.0; 2], &[Some(0.0)])
        .unwrap();
    near(load.stress.force_n_per_m[0], 0.0, 3e-7);
    near(load.stress.force_n_per_m[1], 0.0, 3e-7);
    near(load.stress.moment_n, 0.0, 3e-7);
}

#[test]
fn g1_implicit_quadratic_relaxation_matches_root_and_native_energy_balance() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.4,
        ..FluxParams::default()
    };
    let (k, beta, dt, skeleton, old) = (0.1, 9.0, 0.7, 0.2, 1.2);
    let cells = vec![
        ForchheimerCell::try_new(
            BrinkmanCell::try_new(IntrinsicPermeability::isotropic(k).unwrap(), 1.0)
                .unwrap()
                .with_skeleton_velocity([skeleton, 0.0])
                .unwrap(),
            beta,
        )
        .unwrap();
        mesh.tris.len()
    ];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let previous = constant_dofs(&mesh, [old, 0.0]);
    let saved = previous.clone();
    let solution = system
        .bdf1_brinkman_forchheimer_mixed_step(
            params,
            &cells,
            &|_| [0.0; 2],
            &boundaries,
            &previous,
            dt,
            ForchheimerControls::default(),
        )
        .unwrap();
    assert_eq!(previous, saved);
    let rate = params.nu / k;
    let a = 1.0 + dt * rate;
    let w = 2.0 * (old - skeleton) / (a + (a * a + 4.0 * dt * beta * (old - skeleton)).sqrt());
    let exact = skeleton + w;
    let (error, divergence) = system.velocity_error(&solution.flow.field.x, &|_| [exact, 0.0]);
    assert!(error < 2e-9 && divergence < 2e-9);
    assert!(dt * (rate + 2.0 * beta * (old - skeleton)) > 2.0);
    let fluid_power: f64 = solution
        .cell_darcy_exchange
        .iter()
        .chain(&solution.cell_forchheimer_exchange)
        .map(|exchange| exchange.fluid_power)
        .sum();
    let physical_loss: f64 = solution
        .cell_darcy_exchange
        .iter()
        .chain(&solution.cell_forchheimer_exchange)
        .map(|exchange| exchange.dissipation)
        .sum();
    let skeleton_power: f64 = solution
        .cell_darcy_exchange
        .iter()
        .chain(&solution.cell_forchheimer_exchange)
        .map(|exchange| exchange.skeleton_power)
        .sum();
    near(
        0.5 * (exact * exact - old * old) + 0.5 * (exact - old).powi(2),
        dt * fluid_power,
        3e-9,
    );
    near(fluid_power + skeleton_power + physical_loss, 0.0, 3e-12);
    near(physical_loss, rate * w * w + beta * w.powi(3), 3e-9);
    assert!(physical_loss > 0.0);
}

#[test]
fn g0_zero_beta_retains_exact_steady_and_bdf1_native_fields_and_darcy_accounts() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let base: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| {
            if point[0] > 0.7 {
                BrinkmanCell::free_fluid()
            } else {
                BrinkmanCell::try_new(tensor(0.53), 1.6)
                    .unwrap()
                    .with_skeleton_velocity([0.2, -0.1])
                    .unwrap()
            }
        })
        .collect();
    let cells: Vec<_> = base
        .iter()
        .map(|cell| ForchheimerCell::try_new(*cell, 0.0).unwrap())
        .collect();
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::Velocity {
            edge,
            value: &no_slip,
        })
        .collect();
    let force = |_| [-0.2, 0.8];
    let previous = constant_dofs(&mesh, [0.1, -0.05]);
    for dt in [None, Some(0.3)] {
        let original = if let Some(dt) = dt {
            system.assemble_brinkman_mixed_bdf1(params, &base, &force, &boundaries, &previous, dt)
        } else {
            system.assemble_brinkman_mixed(params, &base, &force, &boundaries)
        }
        .unwrap();
        let expected = original.solve().unwrap();
        let actual = if let Some(dt) = dt {
            system.bdf1_brinkman_forchheimer_mixed_step(
                params,
                &cells,
                &force,
                &boundaries,
                &previous,
                dt,
                ForchheimerControls::default(),
            )
        } else {
            system.solve_brinkman_forchheimer_mixed(
                params,
                &cells,
                &force,
                &boundaries,
                &previous,
                ForchheimerControls::default(),
            )
        }
        .unwrap();
        assert_eq!(actual.flow.field.x, expected.field.x);
        assert_eq!(actual.flow.pressure, expected.pressure);
        assert_eq!(actual.flow.field.rel_residual, expected.field.rel_residual);
        assert_eq!(
            actual.flow.field.scaled_rel_residual,
            expected.field.scaled_rel_residual
        );
        assert_eq!(
            actual.cell_darcy_exchange,
            original
                .cell_darcy_exchange_per_density(&expected.field.x[..system.n_u])
                .unwrap()
        );
        assert_eq!(
            actual.cell_forchheimer_exchange,
            vec![fs_flux::BrinkmanDragExchange::default(); cells.len()]
        );
        assert_eq!(actual.iterations, 0);
        assert_eq!(actual.relative_update, 0.0);
        assert_eq!(actual.velocity_update_norm_m_s, 0.0);
    }
}

#[test]
fn g0_invalid_initial_history_and_controls_refuse_before_callbacks() {
    assert!(ForchheimerCell::try_new(BrinkmanCell::free_fluid(), 1.0).is_err());
    let porous = BrinkmanCell::try_new(tensor(0.0), 1.0).unwrap();
    for beta in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(ForchheimerCell::try_new(porous, beta).is_err());
    }
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let cells = vec![ForchheimerCell::try_new(porous, 1.0).unwrap(); mesh.tris.len()];
    let calls = Cell::new(0);
    let counted = |_: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::Velocity {
            edge,
            value: &counted,
        })
        .collect();
    let zero = vec![0.0; system.n_u];
    for initial in [&[][..], &vec![f64::NAN; system.n_u]] {
        assert!(
            system
                .solve_brinkman_forchheimer_mixed(
                    FluxParams::default(),
                    &cells,
                    &counted,
                    &boundaries,
                    initial,
                    ForchheimerControls::default()
                )
                .is_err()
        );
    }
    for controls in [
        ForchheimerControls {
            max_iterations: 0,
            ..ForchheimerControls::default()
        },
        ForchheimerControls {
            max_backtracks: 61,
            ..ForchheimerControls::default()
        },
        ForchheimerControls {
            relative_tolerance: f64::NAN,
            ..ForchheimerControls::default()
        },
        ForchheimerControls {
            absolute_velocity_tolerance_m_s: -1.0,
            ..ForchheimerControls::default()
        },
        ForchheimerControls {
            absolute_velocity_tolerance_m_s: f64::NAN,
            ..ForchheimerControls::default()
        },
        ForchheimerControls {
            absolute_velocity_tolerance_m_s: f64::INFINITY,
            ..ForchheimerControls::default()
        },
    ] {
        assert!(
            system
                .solve_brinkman_forchheimer_mixed(
                    FluxParams::default(),
                    &cells,
                    &counted,
                    &boundaries,
                    &zero,
                    controls
                )
                .is_err()
        );
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(
            system
                .bdf1_brinkman_forchheimer_mixed_step(
                    FluxParams::default(),
                    &cells,
                    &counted,
                    &boundaries,
                    &zero,
                    dt,
                    ForchheimerControls::default()
                )
                .is_err()
        );
    }
    assert_eq!(calls.get(), 0);
}

#[test]
fn g0_exhausted_newton_or_backtracking_returns_no_accepted_state() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.2,
        ..FluxParams::default()
    };
    let cells = vec![
        ForchheimerCell::try_new(
            BrinkmanCell::try_new(IntrinsicPermeability::isotropic(0.1).unwrap(), 1.0).unwrap(),
            50.0,
        )
        .unwrap();
        mesh.tris.len()
    ];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let initial = vec![0.0; system.n_u];
    let saved = initial.clone();
    for controls in [
        ForchheimerControls {
            max_iterations: 1,
            ..ForchheimerControls::default()
        },
        ForchheimerControls {
            max_backtracks: 0,
            ..ForchheimerControls::default()
        },
    ] {
        assert!(matches!(
            system.solve_brinkman_forchheimer_mixed(
                params,
                &cells,
                &|_| [100.0, 0.0],
                &boundaries,
                &initial,
                controls
            ),
            Err(StokesBoundaryError::Rheology(
                fs_flux::FluxRheologyError::NonlinearConvergence { .. }
            ))
        ));
    }
    assert_eq!(initial, saved);
}

#[test]
fn g0_homogeneous_steady_flow_converges_to_rest_from_nonzero_velocity() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.2,
        ..FluxParams::default()
    };
    let cells = vec![
        ForchheimerCell::try_new(
            BrinkmanCell::try_new(IntrinsicPermeability::isotropic(0.1).unwrap(), 1.0).unwrap(),
            50.0,
        )
        .unwrap();
        mesh.tris.len()
    ];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let initial = constant_dofs(&mesh, [0.8, -0.4]);
    let solution = system
        .solve_brinkman_forchheimer_mixed(
            params,
            &cells,
            &|_| [0.0; 2],
            &boundaries,
            &initial,
            ForchheimerControls::default(),
        )
        .unwrap();
    let (error, divergence) = system.velocity_error(&solution.flow.field.x, &|_| [0.0; 2]);
    assert!(
        error < 2e-9 && divergence < 2e-9,
        "error={error:e}, divergence={divergence:e}"
    );
    assert!(solution.iterations > 1 && solution.flow.field.converged(1e-9));
    assert!(solution.relative_update <= 1e-9);
    assert!(
        solution
            .cell_forchheimer_exchange
            .iter()
            .all(|exchange| exchange.dissipation < 1e-20)
    );
}

#[test]
fn g0_initial_normal_moment_mismatch_is_projected_from_original_boundary_owner() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let base = vec![
        BrinkmanCell::try_new(IntrinsicPermeability::isotropic(0.1).unwrap(), 1.0)
            .unwrap();
        mesh.tris.len()
    ];
    let cells: Vec<_> = base
        .iter()
        .map(|cell| ForchheimerCell::try_new(*cell, 60.0).unwrap())
        .collect();
    let velocity = |x: [f64; 2]| [0.2 + 0.1 * x[1], -0.1];
    let force = |x: [f64; 2]| {
        let u = velocity(x);
        let rate = params.nu / 0.1 + 60.0 * u[0].hypot(u[1]);
        [rate * u[0], rate * u[1]]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    let native = system
        .assemble_brinkman_mixed(params, &base, &force, &boundaries)
        .unwrap();
    let initial = constant_dofs(&mesh, [2.0, -3.0]);
    let solution = system
        .solve_brinkman_forchheimer_mixed(
            params,
            &cells,
            &force,
            &boundaries,
            &initial,
            ForchheimerControls::default(),
        )
        .unwrap();
    for edge in 0..mesh.edges.len() {
        if mesh.is_boundary(edge) {
            assert_eq!(solution.flow.field.x[2 * edge], native.rhs()[2 * edge]);
            assert_eq!(
                solution.flow.field.x[2 * edge + 1],
                native.rhs()[2 * edge + 1]
            );
        }
    }
    let (error, divergence) = system.velocity_error(&solution.flow.field.x, &velocity);
    assert!(error < 2e-9 && divergence < 2e-9);
}

#[test]
fn g1_material_bdf1_keeps_one_density_for_inertia_force_and_cubic_loss() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (density, mu, beta, dt, old) = (7.0, 0.4, 9.0, 0.7, 1.2);
    let material = CountedMaterial {
        state: FluidTransportState::try_new(density, 1000.0, 0.5, mu).unwrap(),
        reads: Cell::new(0),
    };
    let cells = vec![
        ForchheimerCell::try_new(
            BrinkmanCell::try_new(IntrinsicPermeability::isotropic(0.1).unwrap(), 1.0).unwrap(),
            beta,
        )
        .unwrap();
        mesh.tris.len()
    ];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| MaterialStokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let previous = constant_dofs(&mesh, [old, 0.0]);
    let binding =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let solution = binding
        .bdf1_brinkman_forchheimer_mixed_step(
            &cells,
            &|_| [0.0; 2],
            &boundaries,
            &previous,
            dt,
            ForchheimerControls::default(),
        )
        .unwrap();
    assert_eq!(material.reads.get(), 1);
    let a = 1.0 + dt * mu / (density * 0.1);
    let exact = 2.0 * old / (a + (a * a + 4.0 * dt * beta * old).sqrt());
    let (error, divergence) =
        system.velocity_error(&solution.flow().discrete().x, &|_| [exact, 0.0]);
    assert!(error < 2e-9 && divergence < 2e-9);
    near(
        solution.forchheimer_exchange().skeleton_force_n_per_m[0],
        density * beta * exact * exact,
        3e-8,
    );
    near(
        solution.forchheimer_power_w_per_m(),
        density * beta * exact.powi(3),
        3e-8,
    );
    near(
        solution.flow().darcy_power_w_per_m(),
        mu / 0.1 * exact * exact,
        3e-8,
    );
}

#[test]
fn g1_pressure_balanced_rest_accepts_an_explicit_absolute_update_in_steady_and_bdf1() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.3, 0.9, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.023,
        ..FluxParams::default()
    };
    // Synthetic large beta makes roundoff velocity change the Newton tangent.
    // This is a numerical equilibrium regression, not an empirical calibration.
    let cells =
        vec![
            ForchheimerCell::try_new(BrinkmanCell::try_new(tensor(0.37), 1.4).unwrap(), 1e18,)
                .unwrap();
            mesh.tris.len()
        ];
    let gradient = [0.37, -1.13];
    let pressure = |x: [f64; 2]| 9.7 + gradient[0] * x[0] + gradient[1] * x[1];
    let traction = |x: [f64; 2], normal: [f64; 2]| {
        let p = pressure(x);
        [-p * normal[0], -p * normal[1]]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let zero = vec![0.0; system.n_u];
    let controls = ForchheimerControls {
        absolute_velocity_tolerance_m_s: 1e-12,
        ..ForchheimerControls::default()
    };
    for dt in [None, Some(0.37)] {
        let solution = if let Some(dt) = dt {
            system.bdf1_brinkman_forchheimer_mixed_step(
                params,
                &cells,
                &|_| gradient,
                &boundaries,
                &zero,
                dt,
                controls,
            )
        } else {
            system.solve_brinkman_forchheimer_mixed(
                params,
                &cells,
                &|_| gradient,
                &boundaries,
                &zero,
                controls,
            )
        }
        .unwrap();
        let (error, divergence) = system.velocity_error(&solution.flow.field.x, &|_| [0.0; 2]);
        assert!(
            error < 2e-12 && divergence < 2e-12,
            "error={error:e}, divergence={divergence:e}"
        );
        assert!(solution.flow.field.converged(controls.relative_tolerance));
        assert!(solution.velocity_update_norm_m_s.is_finite());
        assert!(solution.velocity_update_norm_m_s <= controls.absolute_velocity_tolerance_m_s);
        for (p, centroid) in solution.flow.pressure.iter().zip(&mesh.centroids) {
            near(*p, pressure(*centroid), 3e-8);
        }
        for exchange in &solution.cell_forchheimer_exchange {
            assert!(
                exchange
                    .skeleton_force
                    .iter()
                    .all(|force| force.abs() < 1e-8)
            );
            assert!(exchange.dissipation >= 0.0 && exchange.dissipation < 1e-18);
        }
    }
    assert!(zero.iter().all(|value| *value == 0.0));
}

#[test]
fn g1_explicit_absolute_update_does_not_hide_weak_nonlinear_flow_under_pressure_offset() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let (beta, driving) = (1e6, 1e-6);
    let cells = vec![
        ForchheimerCell::try_new(
            BrinkmanCell::try_new(IntrinsicPermeability::isotropic(params.nu).unwrap(), 1.0)
                .unwrap(),
            beta,
        )
        .unwrap();
        mesh.tris.len()
    ];
    let exact = 2.0 * driving / (1.0 + (1.0 + 4.0 * beta * driving).sqrt());
    let zero = vec![0.0; system.n_u];
    let controls = ForchheimerControls {
        absolute_velocity_tolerance_m_s: 1e-9,
        ..ForchheimerControls::default()
    };
    for datum in [0.0, 1e5] {
        let traction = |x: [f64; 2], normal: [f64; 2]| {
            let p = datum - driving * x[0];
            [-p * normal[0], -p * normal[1]]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, _)| StokesBoundary::TotalTraction {
                edge,
                value: &traction,
            })
            .collect();
        // The first Newton velocity is driving/alpha = 1e-6, about62% above
        // the nonlinear root. Pressure datum must not relax the update gate.
        assert!(matches!(
            system.solve_brinkman_forchheimer_mixed(
                params,
                &cells,
                &|_| [0.0; 2],
                &boundaries,
                &zero,
                ForchheimerControls {
                    max_iterations: 1,
                    ..controls
                }
            ),
            Err(StokesBoundaryError::Rheology(
                fs_flux::FluxRheologyError::NonlinearConvergence { .. }
            ))
        ));
        let solution = system
            .solve_brinkman_forchheimer_mixed(
                params,
                &cells,
                &|_| [0.0; 2],
                &boundaries,
                &zero,
                controls,
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&solution.flow.field.x, &|_| [exact, 0.0]);
        assert!(
            error < 5e-9 && divergence < 5e-9,
            "datum={datum:e}, error={error:e}, divergence={divergence:e}"
        );
        assert!(
            solution.iterations > 1 && solution.flow.field.converged(controls.relative_tolerance)
        );
        assert!(
            solution.relative_update <= controls.relative_tolerance
                || solution.velocity_update_norm_m_s <= controls.absolute_velocity_tolerance_m_s
        );
        if datum == 0.0 {
            // Even a deliberately coarse absolute update allowance cannot
            // replace the mandatory nonlinear force/continuity residuals.
            assert!(matches!(
                system.solve_brinkman_forchheimer_mixed(
                    params,
                    &cells,
                    &|_| [0.0; 2],
                    &boundaries,
                    &zero,
                    ForchheimerControls {
                        max_iterations: 1,
                        absolute_velocity_tolerance_m_s: 1.0,
                        ..controls
                    }
                ),
                Err(StokesBoundaryError::Rheology(
                    fs_flux::FluxRheologyError::NonlinearConvergence { .. }
                ))
            ));
        }
    }
}

#[test]
fn g1_material_pressure_rest_reports_the_actual_full_absolute_update() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.3, 0.9, 3, 3));
    let system = FluxSystem::new(&mesh);
    let density = 7.0;
    let material = CountedMaterial {
        state: FluidTransportState::try_new(density, 1000.0, 0.5, density * 0.023).unwrap(),
        reads: Cell::new(0),
    };
    let cells =
        vec![
            ForchheimerCell::try_new(BrinkmanCell::try_new(tensor(0.37), 1.4).unwrap(), 1e18,)
                .unwrap();
            mesh.tris.len()
        ];
    let gradient = [0.37, -1.13];
    let pressure = |x: [f64; 2]| density * (9.7 + gradient[0] * x[0] + gradient[1] * x[1]);
    let traction = |x: [f64; 2], normal: [f64; 2]| {
        let p = pressure(x);
        [-p * normal[0], -p * normal[1]]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| MaterialStokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let binding =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let controls = ForchheimerControls {
        absolute_velocity_tolerance_m_s: 1e-12,
        ..ForchheimerControls::default()
    };
    let zero = vec![0.0; system.n_u];
    let solution = binding
        .bdf1_brinkman_forchheimer_mixed_step(
            &cells,
            &|_| gradient.map(|value| density * value),
            &boundaries,
            &zero,
            0.37,
            controls,
        )
        .unwrap();
    assert_eq!(material.reads.get(), 1);
    assert!(
        solution
            .flow()
            .discrete()
            .converged(controls.relative_tolerance)
    );
    assert!(solution.velocity_update_norm_m_s().is_finite());
    assert!(solution.velocity_update_norm_m_s() <= controls.absolute_velocity_tolerance_m_s);
    for (p, centroid) in solution.flow().pressure_pa().iter().zip(&mesh.centroids) {
        near(*p, pressure(*centroid), 2e-7);
    }
}
