//! G1/G3: endpoint-advector density/flow coupling and atomic refusal.

use super::*;
use fs_flux::ns::boundary::StokesRheologyControls;

#[path = "coupled/rheology.rs"]
mod rheology;

#[path = "coupled/windkessel.rs"]
mod windkessel;

fn nonlinear() -> StokesRheologyControls {
    StokesRheologyControls {
        maximum_iterations: 80,
        relative_tolerance: 2e-10,
        relaxation: 0.8,
    }
}

#[test]
fn accelerating_endpoint_transports_density_and_every_inventory_once() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let endpoint = |_: [f64; 2]| [0.4, 0.0];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| Boundary::Velocity {
            edge,
            value: &endpoint,
        })
        .collect();
    let initial_velocity = vec![0.0; system.n_u];
    let initial: Vec<Cell<2>> = mesh
        .centroids
        .iter()
        .map(|point| Cell {
            density_kg_m3: if point[0] < 0.5 { 2.0 } else { 80.0 },
            specific: if point[0] < 0.5 {
                [100.0, 0.2]
            } else {
                [300.0, 0.8]
            },
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| {
            mesh.is_boundary(edge).then_some(Cell {
                density_kg_m3: 8.0,
                specific: [220.0, 0.4],
            })
        })
        .collect();
    let dt = 0.04;
    let viscosity = |_, _, _: &Cell<2>| 0.2;
    // For a uniform endpoint velocity, the shared mass equation makes
    // rho_new*u/dt + div(j*u) = rho_old*u/dt exactly in real arithmetic.
    let force =
        |cell: usize, _: [f64; 2], _: &Cell<2>| [initial[cell].density_kg_m3 * 0.4 / dt, 0.0];
    for reference in [1.0, 80.0] {
        let controls = controls(reference, 1.7);
        let result = system
            .density_momentum_coupled_bdf1_step(
                numerics(),
                &viscosity,
                &force,
                &boundaries,
                &initial_velocity,
                &initial,
                &inflow,
                dt,
                controls,
                nonlinear(),
                || false,
            )
            .unwrap();
        let step = &result.step;
        assert!(result.iterations >= 2);
        assert!(result.relative_update <= nonlinear().relative_tolerance);
        assert!(result.nonlinear_residual <= nonlinear().relative_tolerance);
        let (error, divergence) = system.velocity_error(&step.flow.field.x, &endpoint);
        assert!(
            error < 2e-10 && divergence < 2e-10,
            "{error:e}, {divergence:e}"
        );
        let plan = system
            .prepare_density_transport(
                &step.flow.field.x[..system.n_u],
                controls.depth_m,
                controls.relative_volume_divergence_tolerance,
                controls.maximum_transport_storage_bytes,
                || false,
            )
            .unwrap();
        let mut workspace = plan
            .workspace(controls.maximum_transport_storage_bytes)
            .unwrap();
        let mut expected = vec![Cell::<2>::default(); initial.len()];
        let report = plan
            .advance_into(
                &initial,
                &inflow,
                dt,
                controls.maximum_courant,
                &mut expected,
                &mut workspace,
                || false,
            )
            .unwrap();
        assert_eq!(step.cells, expected);
        assert_eq!(step.transport, report);
        assert!(
            expected
                .iter()
                .zip(&initial)
                .any(|(a, b)| a.density_kg_m3 != b.density_kg_m3 && a.specific != b.specific)
        );
        near(report.mass_balance_residual_kg, 0.0, 2e-12);
        for residual in report.carried_balance_residual {
            near(residual, 0.0, 4e-10);
        }
        // The earlier frozen-old-advector API is a distinct physical map:
        // starting from rest, it performs no material transport this interval.
        let lagged = system
            .density_momentum_bdf1_step(
                numerics(),
                &viscosity,
                &force,
                &boundaries,
                &initial_velocity,
                &initial_velocity,
                &initial,
                &inflow,
                dt,
                controls,
                || false,
            )
            .unwrap();
        assert_eq!(lagged.cells, initial);
        assert_eq!(lagged.transport.inlet_mass_kg, 0.0);
        assert!(step.transport.inlet_mass_kg > 0.0);
    }
    assert!(initial_velocity.iter().all(|value| *value == 0.0));
    assert_eq!(initial[0].specific, [100.0, 0.2]);
}

#[test]
fn coupled_closed_counterflow_preserves_inventory_and_native_kinetic_balance() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.2 / (1.5 * mesh.edges[internal].len);
    let original_velocity = previous.clone();
    let cells = [
        Cell {
            density_kg_m3: 2.0,
            specific: [100.0, 0.2],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [300.0, 0.8],
        },
    ];
    let zero = |_: [f64; 2]| [0.0; 2];
    let selected: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect();
    let boundaries: Vec<_> = selected
        .iter()
        .map(|&edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let native: Vec<_> = selected
        .iter()
        .map(|&edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let reference = 4.0;
    let depth = 1.7;
    let dt = 0.05;
    let mu = 0.03;
    let viscous = system
        .assemble_stokes_mixed(
            FluxParams {
                nu: mu / reference,
                sigma: numerics().penalty,
                tol: numerics().linear_tolerance,
                max_iters: numerics().max_linear_iterations,
            },
            &zero,
            &native,
        )
        .unwrap();
    let result = system
        .density_momentum_coupled_bdf1_step(
            numerics(),
            &|_, _, _| mu,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &cells,
            &vec![None; mesh.edges.len()],
            dt,
            controls(reference, depth),
            nonlinear(),
            || false,
        )
        .unwrap();
    let step = result.step;
    let mut applied = vec![0.0; step.flow.field.x.len()];
    viscous.matrix().spmv(&step.flow.field.x, &mut applied);
    let power: f64 = step.flow.field.x[..system.n_u]
        .iter()
        .zip(&applied[..system.n_u])
        .map(|(a, b)| a * b)
        .sum::<f64>()
        * reference
        * depth;
    let observations = &step.observations;
    near(
        observations.kinetic_energy_after_j - observations.kinetic_energy_before_j
            + observations.temporal_dissipation_j
            + observations.interior_upwind_dissipation_j
            + dt * power,
        0.0,
        2e-10,
    );
    assert!(observations.kinetic_energy_after_j < observations.kinetic_energy_before_j);
    near(step.transport.mass_balance_residual_kg, 0.0, 2e-13);
    for residual in step.transport.carried_balance_residual {
        near(residual, 0.0, 3e-11);
    }
    assert_eq!(previous, original_velocity);
    assert_eq!(cells[0].density_kg_m3, 2.0);
}

#[test]
fn coupled_trial_refusals_do_not_publish_or_advance_old_state() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let endpoint = |_: [f64; 2]| [0.4, 0.0];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| Boundary::Velocity {
            edge,
            value: &endpoint,
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    let cells = [
        Cell::<1> {
            density_kg_m3: 2.0,
            specific: [100.0],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [300.0],
        },
    ];
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| mesh.is_boundary(edge).then_some(cells[0]))
        .collect();
    let force = |cell: usize, _: [f64; 2], _: &Cell<1>| [cells[cell].density_kg_m3 * 10.0, 0.0];
    let run = |control, solve, exterior: &[Option<Cell<1>>], cancel: &mut dyn FnMut() -> bool| {
        system.density_momentum_coupled_bdf1_step(
            numerics(),
            &|_, _, _| 0.1,
            &force,
            &boundaries,
            &previous,
            &cells,
            exterior,
            0.04,
            control,
            solve,
            cancel,
        )
    };
    assert!(matches!(
        run(controls(4.0, 1.0), nonlinear(), &inflow, &mut || true),
        Err(Error::Cancelled)
    ));
    assert!(
        run(
            controls(4.0, 1.0),
            StokesRheologyControls {
                maximum_iterations: 1,
                ..nonlinear()
            },
            &inflow,
            &mut || false,
        )
        .is_err()
    );
    assert!(
        run(
            controls(4.0, 1.0),
            nonlinear(),
            &vec![None; mesh.edges.len()],
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        run(
            Controls {
                maximum_courant: 1e-6,
                ..controls(4.0, 1.0)
            },
            nonlinear(),
            &inflow,
            &mut || false,
        ),
        Err(Error::Transport(DensityTransportError::Courant { .. }))
    ));
    for solve in [
        StokesRheologyControls {
            maximum_iterations: 0,
            ..nonlinear()
        },
        StokesRheologyControls {
            relaxation: 0.0,
            ..nonlinear()
        },
        StokesRheologyControls {
            relative_tolerance: f64::NAN,
            ..nonlinear()
        },
    ] {
        assert!(run(controls(4.0, 1.0), solve, &inflow, &mut || false).is_err());
    }
    let mut polls = 0;
    run(controls(4.0, 1.0), nonlinear(), &inflow, &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        run(controls(4.0, 1.0), nonlinear(), &inflow, &mut || {
            seen += 1;
            seen == polls
        }),
        Err(Error::Cancelled)
    ));
    assert!(previous.iter().all(|u| *u == 0.0));
    assert_eq!(cells[0].specific, [100.0]);
    assert_eq!(cells[1].density_kg_m3, 8.0);
}

#[test]
fn controlled_port_couples_heterogeneous_density_pressure_and_physical_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let endpoint = |_: [f64; 2]| [0.4, 0.0];
    let zero = |_: [f64; 2]| [0.0; 2];
    let base_pressure = 7.0;
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -base_pressure * n);
    let west: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].mid[0] == 0.0)
        .collect();
    let open: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| {
            mesh.is_boundary(edge)
                && (mesh.edges[edge].mid[0] == 0.0 || mesh.edges[edge].mid[0] == 1.0)
        })
        .collect();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if open.contains(&edge) {
                Boundary::OpenReservoir {
                    edge,
                    traction: &traction,
                    velocity: &endpoint,
                }
            } else {
                Boundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_pa_s_per_m: 0.0,
                }
            }
        })
        .collect();
    let initial_velocity = vec![0.0; system.n_u];
    let initial: Vec<Cell<2>> = mesh
        .centroids
        .iter()
        .map(|point| Cell {
            density_kg_m3: if point[0] < 0.5 { 2.0 } else { 80.0 },
            specific: if point[0] < 0.5 {
                [100.0, 0.2]
            } else {
                [300.0, 0.8]
            },
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| {
            open.contains(&edge).then_some(Cell {
                density_kg_m3: 8.0,
                specific: [220.0, 0.4],
            })
        })
        .collect();
    let dt = 0.04;
    let increment = 3.0;
    // The SAME mass equation cancels the new-density/convection terms for
    // this constant endpoint. The remaining body load and port pressure
    // recover rho_old*U/dt, even while the transported density changes.
    let force = |cell: usize, _: [f64; 2], _: &Cell<2>| {
        [initial[cell].density_kg_m3 * 0.4 / dt - increment, 0.0]
    };
    for (reference, depth) in [(1.0, 1.7), (80.0, 1.7), (80.0, 0.5)] {
        let controls = controls(reference, depth);
        let result = system
            .density_momentum_coupled_bdf1_flow_control_step(
                numerics(),
                &|_, _, _| 0.2,
                &force,
                &boundaries,
                &[&west],
                &[-0.4],
                2e-11,
                &initial_velocity,
                &initial,
                &inflow,
                dt,
                controls,
                nonlinear(),
                || false,
            )
            .unwrap();
        let coupled = &result.step;
        let step = &coupled.step;
        assert!(coupled.iterations >= 2);
        assert!(coupled.relative_update <= nonlinear().relative_tolerance);
        assert!(coupled.nonlinear_residual <= nonlinear().relative_tolerance);
        assert!(result.maximum_flow_residual_m2_s <= 2e-11);
        near(result.outward_volume_rates_m2_s[0], -0.4, 2e-11);
        near(result.pressure_offsets_pa[0], increment, 5e-8);
        let (error, divergence) = system.velocity_error(&step.flow.field.x, &endpoint);
        assert!(
            error < 2e-10 && divergence < 2e-10,
            "{error:e}, {divergence:e}"
        );
        for (&pressure, centroid) in step.pressure_pa.iter().zip(&mesh.centroids) {
            near(
                pressure,
                base_pressure + increment * (1.0 - centroid[0]),
                5e-8,
            );
        }
        let plan = system
            .prepare_density_transport(
                &step.flow.field.x[..system.n_u],
                depth,
                controls.relative_volume_divergence_tolerance,
                controls.maximum_transport_storage_bytes,
                || false,
            )
            .unwrap();
        let mut workspace = plan
            .workspace(controls.maximum_transport_storage_bytes)
            .unwrap();
        let mut expected = vec![Cell::<2>::default(); initial.len()];
        let report = plan
            .advance_into(
                &initial,
                &inflow,
                dt,
                controls.maximum_courant,
                &mut expected,
                &mut workspace,
                || false,
            )
            .unwrap();
        // The reported material state is exactly the endpoint's transport,
        // including the compatible flux used by its pressure-controlled solve.
        assert_eq!(step.cells, expected);
        assert_eq!(step.transport, report);
        assert!(expected.iter().zip(&initial).any(|(a, b)| a != b));
        near(report.inlet_mass_kg, 8.0 * 0.4 * depth * dt, 2e-11);
        near(report.outlet_mass_kg, 80.0 * 0.4 * depth * dt, 2e-11);
        near(report.mass_balance_residual_kg, 0.0, 2e-12);
        for residual in report.carried_balance_residual {
            near(residual, 0.0, 4e-10);
        }
        let loads = step
            .flow
            .reactions
            .resultant(&open, [0.0; 2], &[None])
            .unwrap();
        // Native reactions act ON the exterior. Restore physical stress by
        // the common reference and depth; both base traction and increments
        // contribute, and the common base pressure cancels over the two ports.
        let delivered_power =
            -loads.imposed.rigid_power(endpoint([0.0; 2]), 0.0).unwrap() * reference * depth;
        near(delivered_power, increment * 0.4 * depth, 5e-8);
        near(
            loads.stress.force[0] * reference * depth,
            -increment * depth,
            5e-8,
        );
    }
    assert!(initial_velocity.iter().all(|value| *value == 0.0));
    assert_eq!(initial[0].density_kg_m3, 2.0);
}

#[test]
fn controlled_volume_flow_requires_actual_density_inertia_and_reversed_donors() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let west: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].mid[0] == 0.0)
        .collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let previous = vec![0.0; system.n_u];
    let dt = 0.04;
    let depth = 2.5;
    for density in [1.2, 998.0] {
        for speed in [-0.3, 0.3] {
            let endpoint = |_: [f64; 2]| [speed, 0.0];
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|&edge| mesh.is_boundary(edge))
                .map(|edge| {
                    if mesh.edges[edge].mid[0] == 0.0 || mesh.edges[edge].mid[0] == 1.0 {
                        Boundary::OpenReservoir {
                            edge,
                            traction: &traction,
                            velocity: &endpoint,
                        }
                    } else {
                        Boundary::NavierSlip {
                            edge,
                            value: &zero,
                            tangential_resistance_pa_s_per_m: 0.0,
                        }
                    }
                })
                .collect();
            let cells = vec![
                Cell::<1> {
                    density_kg_m3: density,
                    specific: [12.0]
                };
                mesh.tris.len()
            ];
            let inflow: Vec<_> = (0..mesh.edges.len())
                .map(|edge| mesh.is_boundary(edge).then_some(cells[0]))
                .collect();
            let result = system
                .density_momentum_coupled_bdf1_flow_control_step(
                    numerics(),
                    &|_, _, _| 0.2,
                    &|_, _, _| [0.0; 2],
                    &boundaries,
                    &[&west],
                    &[-speed],
                    2e-11,
                    &previous,
                    &cells,
                    &inflow,
                    dt,
                    controls(17.0, depth),
                    nonlinear(),
                    || false,
                )
                .unwrap();
            let step = result.step.step;
            let pressure = density * speed / dt;
            near(result.pressure_offsets_pa[0], pressure, 2e-7);
            let (error, divergence) = system.velocity_error(&step.flow.field.x, &endpoint);
            assert!(
                error < 2e-10 && divergence < 2e-10,
                "{error:e}, {divergence:e}"
            );
            near(
                step.transport.inlet_mass_kg,
                density * speed.abs() * depth * dt,
                2e-10,
            );
            near(
                step.transport.outlet_mass_kg,
                density * speed.abs() * depth * dt,
                2e-10,
            );
            near(step.transport.mass_balance_residual_kg, 0.0, 2e-11);
            for cell in &step.cells {
                near(cell.density_kg_m3, density, 2e-10);
                near(cell.specific[0], 12.0, 2e-12);
            }
            // With reversed target the uncontrolled east datum supplies
            // incoming material. A donor only at the nominal inlet is invalid.
            if speed < 0.0 {
                let only_west: Vec<_> = (0..mesh.edges.len())
                    .map(|edge| west.contains(&edge).then_some(cells[0]))
                    .collect();
                assert!(matches!(
                    system.density_momentum_coupled_bdf1_flow_control_step(
                        numerics(),
                        &|_, _, _| 0.2,
                        &|_, _, _| [0.0; 2],
                        &boundaries,
                        &[&west],
                        &[-speed],
                        2e-11,
                        &previous,
                        &cells,
                        &only_west,
                        dt,
                        controls(17.0, depth),
                        nonlinear(),
                        || false,
                    ),
                    Err(Error::Transport(_))
                ));
            }
        }
    }
    assert!(previous.iter().all(|value| *value == 0.0));
}

#[test]
fn controlled_density_trials_preserve_atomic_refusal_and_native_port_admission() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let endpoint = |_: [f64; 2]| [0.3, 0.0];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let zero = |_: [f64; 2]| [0.0; 2];
    let west: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].mid[0] == 0.0)
        .collect();
    let east: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge) && mesh.edges[edge].mid[0] == 1.0)
        .collect();
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if west.contains(&edge) || east.contains(&edge) {
                Boundary::OpenReservoir {
                    edge,
                    traction: &traction,
                    velocity: &endpoint,
                }
            } else {
                Boundary::NavierSlip {
                    edge,
                    value: &zero,
                    tangential_resistance_pa_s_per_m: 0.0,
                }
            }
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    let cells = vec![
        Cell::<1> {
            density_kg_m3: 8.0,
            specific: [31.0]
        };
        mesh.tris.len()
    ];
    let original = cells.clone();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| mesh.is_boundary(edge).then_some(cells[0]))
        .collect();
    let run =
        |ports: &[&[usize]], targets: &[f64], control, solve, cancel: &mut dyn FnMut() -> bool| {
            system.density_momentum_coupled_bdf1_flow_control_step(
                numerics(),
                &|_, _, _| 0.2,
                &|_, _, _| [0.0; 2],
                &boundaries,
                ports,
                targets,
                2e-11,
                &previous,
                &cells,
                &inflow,
                0.04,
                control,
                solve,
                cancel,
            )
        };
    let control = controls(4.0, 1.0);
    assert!(
        run(
            &[&west, &east],
            &[-0.3, 0.3],
            control,
            nonlinear(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        run(
            &[&west, &west],
            &[-0.3, -0.3],
            control,
            nonlinear(),
            &mut || false
        )
        .is_err()
    );
    assert!(run(&[&west], &[], control, nonlinear(), &mut || false).is_err());
    assert!(run(&[&west], &[f64::NAN], control, nonlinear(), &mut || false).is_err());
    assert!(
        run(
            &[&west],
            &[-0.3],
            control,
            StokesRheologyControls {
                maximum_iterations: 1,
                ..nonlinear()
            },
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        run(
            &[&west],
            &[-0.3],
            Controls {
                maximum_courant: 1e-6,
                ..control
            },
            nonlinear(),
            &mut || false
        ),
        Err(Error::Transport(DensityTransportError::Courant { .. }))
    ));
    assert!(matches!(
        run(&[&west], &[-0.3], control, nonlinear(), &mut || true),
        Err(Error::Cancelled)
    ));
    let mut polls = 0;
    run(&[&west], &[-0.3], control, nonlinear(), &mut || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        run(&[&west], &[-0.3], control, nonlinear(), &mut || {
            seen += 1;
            seen == polls
        }),
        Err(Error::Cancelled)
    ));
    assert!(previous.iter().all(|value| *value == 0.0));
    assert_eq!(cells, original);
}
