//! Transient drag regressions run through the original mixed momentum owner.

use super::*;

fn uniform_moments(mesh: &TriMesh, velocity: [f64; 2]) -> Vec<f64> {
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

fn relaxed_velocity(previous: [f64; 2], skeleton: [f64; 2], dt: f64) -> [f64; 2] {
    let (s, c) = 0.53_f64.sin_cos();
    let relative = [previous[0] - skeleton[0], previous[1] - skeleton[1]];
    // Independent diagonal solution of the uniform physical ODE in material axes.
    let a = (c * relative[0] + s * relative[1]) / (1.0 + dt * 10.0);
    let b = (-s * relative[0] + c * relative[1]) / (1.0 + dt * 40.0);
    [skeleton[0] + c * a - s * b, skeleton[1] + s * a + c * b]
}

#[test]
fn g1_implicit_rotated_relaxation_retains_history_and_signed_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let traction = |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0], -7.0 * n[1]];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let mut reference = Vec::new();
    for (frame, boost) in [[0.0, 0.0], [0.4, -0.1]].into_iter().enumerate() {
        let mut expected = boost;
        let mut previous = uniform_moments(&mesh, expected);
        // dt*lambda_max=8: a forward-Euler drag split would be unstable.
        for (step, skeleton) in [[0.8, -0.4], [-0.3, 0.6]].into_iter().enumerate() {
            let skeleton = [skeleton[0] + boost[0], skeleton[1] + boost[1]];
            let cells = vec![
                BrinkmanCell::try_new(tensor(0.53), 1.7)
                    .unwrap()
                    .with_skeleton_velocity(skeleton)
                    .unwrap();
                mesh.tris.len()
            ];
            let history = previous.clone();
            let assembly = system
                .assemble_brinkman_mixed_bdf1(
                    params,
                    &cells,
                    &|_| [0.0; 2],
                    &boundaries,
                    &previous,
                    0.2,
                )
                .unwrap();
            let flow = assembly.solve().unwrap();
            expected = relaxed_velocity(expected, skeleton, 0.2);
            let (error, divergence) = system.velocity_error(&flow.field.x, &|_| expected);
            assert!(error < 2e-10 && divergence < 2e-10);
            assert_eq!(
                flow.pressure_reference,
                StokesPressureReference::TractionDatum
            );
            for pressure in &flow.pressure {
                near(*pressure, 7.0, 2e-9);
            }
            assert_eq!(previous, history);
            let next = &flow.field.x[..system.n_u];
            let change: Vec<_> = next.iter().zip(&previous).map(|(u, old)| u - old).collect();
            let exchanges = assembly.cell_darcy_exchange_per_density(next).unwrap();
            let loss: f64 = exchanges.iter().map(|e| e.dissipation).sum();
            let solid_work: f64 = exchanges.iter().map(|e| e.skeleton_power).sum();
            let fluid_work: f64 = exchanges.iter().map(|e| e.fluid_power).sum();
            let force: [f64; 2] =
                core::array::from_fn(|a| exchanges.iter().map(|e| e.skeleton_force[a]).sum());
            near(fluid_work + solid_work + loss, 0.0, 2e-12);
            // The true mass quadrature closes the BE energy identity. Its
            // positive numerical loss must not be reported as physical heating.
            near(
                system.kinetic_energy(next) - system.kinetic_energy(&previous)
                    + system.kinetic_energy(&change)
                    + 0.2 * (loss + solid_work),
                0.0,
                2e-10,
            );
            if frame == 0 {
                reference.push((force, loss));
                if step == 0 {
                    assert!(fluid_work > 0.0 && solid_work < 0.0);
                }
            } else {
                for a in 0..2 {
                    near(force[a], reference[step].0[a], 2e-9);
                }
                near(loss, reference[step].1, 2e-9);
            }
            previous = next.to_vec();
        }
    }
}

#[test]
fn g2_uniform_drag_relaxation_converges_to_continuous_time_at_first_order() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.02,
        ..FluxParams::default()
    };
    let skeleton = [0.3, -0.2];
    let (s, c) = 0.53_f64.sin_cos();
    let initial = [skeleton[0] + c - 0.5 * s, skeleton[1] + s + 0.5 * c];
    let exact = [
        skeleton[0] + c * (-2.0_f64).exp() - 0.5 * s * (-8.0_f64).exp(),
        skeleton[1] + s * (-2.0_f64).exp() + 0.5 * c * (-8.0_f64).exp(),
    ];
    let cells = vec![
        BrinkmanCell::try_new(tensor(0.53), 1.0)
            .unwrap()
            .with_skeleton_velocity(skeleton)
            .unwrap();
        mesh.tris.len()
    ];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut previous = uniform_moments(&mesh, initial);
        for _ in 0..steps {
            let flow = system
                .assemble_brinkman_mixed_bdf1(
                    params,
                    &cells,
                    &|_| [0.0; 2],
                    &boundaries,
                    &previous,
                    0.2 / f64::from(steps),
                )
                .unwrap()
                .solve()
                .unwrap();
            previous = flow.field.x[..system.n_u].to_vec();
        }
        let (error, divergence) = system.velocity_error(&previous, &|_| exact);
        assert!(divergence < 2e-10);
        errors.push(error);
    }
    for pair in errors.windows(2) {
        assert!((1.7..2.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn g3_transient_cellwise_motion_preserves_source_and_constrained_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (density, mu, dt) = (7.0, 0.2, 0.03);
    let material = CountedMaterial {
        state: FluidTransportState::try_new(density, 1000.0, 0.5, mu).unwrap(),
        reads: Cell::new(0),
    };
    let initial = [-0.1, 0.2];
    let endpoint = [0.4, -0.1];
    let previous = uniform_moments(&mesh, initial);
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| {
            if point[0] > 0.5 {
                BrinkmanCell::free_fluid()
            } else {
                BrinkmanCell::try_new(tensor(0.53), 1.4)
                    .unwrap()
                    .with_skeleton_velocity([0.8 - 0.3 * point[1], 0.1 + 0.2 * point[0]])
                    .unwrap()
            }
        })
        .collect();
    let (s, c) = 0.53_f64.sin_cos();
    let rxx = c * c / 0.002 + s * s / 0.0005;
    let rxy = s * c * (1.0 / 0.002 - 1.0 / 0.0005);
    let ryy = s * s / 0.002 + c * c / 0.0005;
    let resistance = |cell: BrinkmanCell| {
        if cell.permeability().is_none() {
            return [0.0; 2];
        }
        let solid = cell.skeleton_velocity();
        let r = [endpoint[0] - solid[0], endpoint[1] - solid[1]];
        [
            mu * (rxx * r[0] + rxy * r[1]),
            mu * (rxy * r[0] + ryy * r[1]),
        ]
    };
    let force = |point: [f64; 2]| {
        let drag = resistance(cells[mesh.containing_cell(point).unwrap()]);
        core::array::from_fn(|a| density * (endpoint[a] - initial[a]) / dt + drag[a])
    };
    let velocity = |_: [f64; 2]| endpoint;
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| MaterialStokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    let binding =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let solution = binding
        .bdf1_brinkman_mixed_step(&cells, &force, &boundaries, &previous, dt)
        .unwrap();
    assert_eq!(material.reads.get(), 1);
    assert_eq!(solution.transport(), &material.state);
    assert_eq!(solution.cells(), cells.as_slice());
    assert_eq!(previous, uniform_moments(&mesh, initial));
    let (error, divergence) = system.velocity_error(&solution.discrete().x, &velocity);
    assert!(error < 2e-10 && divergence < 2e-10);
    for pressure in solution.pressure_pa() {
        near(*pressure, 0.0, 2e-8);
    }
    for (i, exchange) in solution.cell_darcy_exchange().iter().enumerate() {
        let force = resistance(cells[i]).map(|value| mesh.areas[i] * value);
        let solid = cells[i].skeleton_velocity();
        for a in 0..2 {
            near(exchange.skeleton_force_n_per_m[a], force[a], 2e-9);
        }
        near(
            exchange.dissipation_w_per_m,
            force[0] * (endpoint[0] - solid[0]) + force[1] * (endpoint[1] - solid[1]),
            2e-9,
        );
        near(
            exchange.fluid_power_w_per_m,
            -force[0] * endpoint[0] - force[1] * endpoint[1],
            2e-9,
        );
        near(
            exchange.skeleton_power_w_per_m,
            force[0] * solid[0] + force[1] * solid[1],
            2e-9,
        );
    }
    // With uniform endpoint velocity and zero pressure, true boundary stress
    // vanishes. Missing mass/history/solid loads in reaction rows would not.
    let selected: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect();
    let load = solution
        .boundary_loads(&selected, [0.0; 2], &[Some(0.0)])
        .unwrap();
    near(load.stress.force_n_per_m[0], 0.0, 2e-8);
    near(load.stress.force_n_per_m[1], 0.0, 2e-8);
    near(load.stress.moment_n, 0.0, 2e-8);
}

#[test]
fn g0_transient_free_fluid_is_native_stokes_and_bad_history_never_samples_loads() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let cells = vec![BrinkmanCell::free_fluid(); mesh.tris.len()];
    let previous = uniform_moments(&mesh, [0.7, -0.4]);
    let sampled = Cell::new(0);
    let force = |_: [f64; 2]| {
        sampled.set(sampled.get() + 1);
        [0.0; 2]
    };
    let traction = |_: [f64; 2], n: [f64; 2]| {
        sampled.set(sampled.get() + 1);
        [-2.0 * n[0], -2.0 * n[1]]
    };
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| StokesBoundary::TotalTraction {
            edge,
            value: &traction,
        })
        .collect();
    let stokes = system
        .assemble_stokes_mixed_bdf1(params, &force, &boundaries, &previous, 0.1)
        .unwrap();
    let brinkman = system
        .assemble_brinkman_mixed_bdf1(params, &cells, &force, &boundaries, &previous, 0.1)
        .unwrap();
    assert_eq!(stokes.matrix().to_dense(), brinkman.matrix().to_dense());
    assert_eq!(stokes.rhs(), brinkman.rhs());
    let flow = brinkman.solve().unwrap();
    let (error, divergence) = system.velocity_error(&flow.field.x, &|_| [0.7, -0.4]);
    assert!(error < 2e-10 && divergence < 2e-10);
    assert_eq!(
        brinkman
            .cell_darcy_power_per_density(&flow.field.x[..system.n_u])
            .unwrap(),
        vec![0.0; mesh.tris.len()]
    );
    sampled.set(0);
    for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(
            system
                .assemble_brinkman_mixed_bdf1(params, &cells, &force, &boundaries, &previous, dt)
                .is_err()
        );
    }
    for bad in [
        vec![],
        vec![f64::NAN; system.n_u],
        vec![f64::INFINITY; system.n_u],
    ] {
        let before: Vec<_> = bad.iter().map(|value| value.to_bits()).collect();
        assert!(
            system
                .assemble_brinkman_mixed_bdf1(params, &cells, &force, &boundaries, &bad, 0.1)
                .is_err()
        );
        assert_eq!(
            bad.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            before
        );
    }
    assert!(
        system
            .assemble_brinkman_mixed_bdf1(params, &cells[1..], &force, &boundaries, &previous, 0.1)
            .is_err()
    );
    assert_eq!(sampled.get(), 0);
    assert_eq!(previous, uniform_moments(&mesh, [0.7, -0.4]));
}
