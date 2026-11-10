//! G1/G3: gravity must use transported mass and the native pressure balance.
use super::*;

#[test]
fn transported_density_has_one_physical_free_fall_acceleration() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 0.6, 3, 2));
    let system = FluxSystem::new(&mesh);
    let density: Vec<_> = (0..mesh.tris.len())
        .map(|cell| Cell {
            density_kg_m3: 2.0 + 10.0 * cell as f64,
            specific: [],
        })
        .collect();
    let u = [0.7, -0.2];
    let velocity = constant(&mesh, u);
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .map(|face| {
            (face.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 321.0,
                specific: [],
            })
        })
        .collect();
    let exterior: Vec<_> = mesh
        .edges
        .iter()
        .map(|face| (face.tris.1 == usize::MAX).then_some([u; 2]))
        .collect();
    let traction = |_, _| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::TotalTraction {
                edge,
                value: &traction,
            })
        })
        .collect();
    let gravity = [1.7, -3.2];
    let mut reference_velocity: Option<Vec<f64>> = None;
    for reference in [0.5, 7.0, 1000.0] {
        let setting = controls(reference);
        let accepted = system
            .assemble_conservative_momentum_bdf1_with_gravity(
                History {
                    density: &density,
                    velocity: &velocity,
                    inflow_density: &inlet,
                    inflow_velocity: &exterior,
                },
                setting,
                MaterialFluxNumerics::default(),
                &|_, _| 0.03,
                &|_| [0.0; 2],
                gravity,
                &boundaries,
                || false,
            )
            .unwrap()
            .solve(|| false)
            .unwrap();
        assert!(accepted.density().iter().zip(&density).any(|(a, b)| a != b));
        let expected_u = [
            u[0] + setting.timestep_s * gravity[0],
            u[1] + setting.timestep_s * gravity[1],
        ];
        let expected = constant(&mesh, expected_u);
        for (&actual, expected) in accepted.velocity_m_s().iter().zip(expected) {
            near(actual, expected, 3e-10);
        }
        for &pressure in accepted.pressure_pa() {
            near(pressure, 0.0, 3e-8);
        }
        let mass = accepted.mass_transport().mass_after_kg / setting.depth_m;
        let actual_momentum = momentum(&system, accepted.density(), accepted.velocity_m_s());
        for axis in 0..2 {
            near(actual_momentum[axis], mass * expected_u[axis], 3e-9);
        }
        if let Some(previous) = &reference_velocity {
            for (&a, &b) in accepted.velocity_m_s().iter().zip(previous) {
                near(a, b, 3e-10);
            }
        } else {
            reference_velocity = Some(accepted.velocity_m_s().to_vec());
        }
    }
}

#[test]
fn layered_hydrostatic_pressure_balances_gravity_without_spurious_flow() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 4));
    let system = FluxSystem::new(&mesh);
    let density: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| Cell {
            density_kg_m3: if point[1] < 0.5 { 1100.0 } else { 900.0 },
            specific: [],
        })
        .collect();
    let velocity = vec![0.0; system.n_u];
    let inlet = vec![None; mesh.edges.len()];
    let exterior = vec![None; mesh.edges.len()];
    let wall = |_| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::Velocity {
                edge,
                value: &wall,
            })
        })
        .collect();
    let gravity = [0.0, -9.81];
    let accepted = system
        .assemble_conservative_momentum_bdf1_with_gravity(
            History {
                density: &density,
                velocity: &velocity,
                inflow_density: &inlet,
                inflow_velocity: &exterior,
            },
            controls(1000.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.0; 2],
            gravity,
            &boundaries,
            || false,
        )
        .unwrap()
        .solve(|| false)
        .unwrap();
    assert_eq!(accepted.density(), density);
    for &value in accepted.velocity_m_s() {
        near(value, 0.0, 2e-9);
    }
    // Continuous, piecewise linear pressure, projected onto cell-average P0.
    // The interface lies exactly on a mesh row. The sealed component removes
    // only its area-weighted constant pressure gauge.
    let expected: Vec<_> = mesh
        .centroids
        .iter()
        .zip(&density)
        .map(|(point, rho)| rho.density_kg_m3 * gravity[1] * (point[1] - 0.5))
        .collect();
    let area: f64 = mesh.areas.iter().sum();
    let mean = expected
        .iter()
        .zip(&mesh.areas)
        .map(|(p, a)| p * a)
        .sum::<f64>()
        / area;
    for (&actual, &expected) in accepted.pressure_pa().iter().zip(&expected) {
        near(actual, expected - mean, 3e-7);
    }
    assert!(accepted.kinetic_energy_j_per_m() < 1e-12);
}

#[test]
fn lateral_density_contrast_sinks_heavy_fluid_and_rises_light_fluid() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let density: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| Cell {
            density_kg_m3: if point[0] < 0.5 { 1200.0 } else { 800.0 },
            specific: [],
        })
        .collect();
    let velocity = vec![0.0; system.n_u];
    let inlet = vec![None; mesh.edges.len()];
    let exterior = vec![None; mesh.edges.len()];
    let wall = |_| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::Velocity {
                edge,
                value: &wall,
            })
        })
        .collect();
    let solve = |gravity| {
        system
            .assemble_conservative_momentum_bdf1_with_gravity(
                History {
                    density: &density,
                    velocity: &velocity,
                    inflow_density: &inlet,
                    inflow_velocity: &exterior,
                },
                controls(1000.0),
                MaterialFluxNumerics::default(),
                &|_, _| 0.1,
                &|_| [0.0; 2],
                gravity,
                &boundaries,
                || false,
            )
            .unwrap()
            .solve(|| false)
            .unwrap()
    };
    let down = solve([0.0, -9.81]);
    let up = solve([0.0, 9.81]);
    let mut vertical = [0.0; 2];
    for (cell, nodes) in mesh.tris.iter().enumerate() {
        let side = usize::from(mesh.centroids[cell][0] >= 0.5);
        for (point, weight) in tri_quad(nodes.map(|node| mesh.verts[node]), mesh.areas[cell]) {
            for local in 0..6 {
                let basis = eval_basis(&system.bases[cell], local, point);
                let row = 2 * mesh.tri_edges[cell][local / 2].0 + local % 2;
                vertical[side] += weight * basis[1] * down.velocity_m_s()[row];
            }
        }
    }
    assert!(vertical[0] < -1e-6, "heavy fluid must sink: {vertical:?}");
    assert!(vertical[1] > 1e-6, "light fluid must rise: {vertical:?}");
    near(vertical[0] + vertical[1], 0.0, 2e-10);
    for (&a, &b) in down.velocity_m_s().iter().zip(up.velocity_m_s()) {
        near(a, -b, 3e-10);
    }
}

#[test]
fn zero_gravity_is_bitwise_unchanged_and_invalid_gravity_is_atomic() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let density = vec![
        Cell {
            density_kg_m3: 2.0,
            specific: [],
        },
        Cell {
            density_kg_m3: 9.0,
            specific: [],
        },
    ];
    let velocity = vec![0.0; system.n_u];
    let inlet = vec![None; mesh.edges.len()];
    let exterior = vec![None; mesh.edges.len()];
    let wall = |_| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::Velocity {
                edge,
                value: &wall,
            })
        })
        .collect();
    let history = || History {
        density: &density,
        velocity: &velocity,
        inflow_density: &inlet,
        inflow_velocity: &exterior,
    };
    let original = system
        .assemble_conservative_momentum_bdf1(
            history(),
            controls(1.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.2, -0.3],
            &boundaries,
            || false,
        )
        .unwrap();
    let zero = system
        .assemble_conservative_momentum_bdf1_with_gravity(
            history(),
            controls(1.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.2, -0.3],
            [0.0; 2],
            &boundaries,
            || false,
        )
        .unwrap();
    assert_eq!(original.rhs(), zero.rhs());
    let a = original.solve(|| false).unwrap();
    let b = zero.solve(|| false).unwrap();
    assert_eq!(a.field().x, b.field().x);
    assert_eq!(a.density(), b.density());
    for invalid in [[f64::NAN, 0.0], [0.0, f64::INFINITY]] {
        assert!(system
            .assemble_conservative_momentum_bdf1_with_gravity(
                history(),
                controls(1.0),
                MaterialFluxNumerics::default(),
                &|_, _| panic!("invalid gravity must precede viscosity callbacks"),
                &|_| panic!("invalid gravity must precede force callbacks"),
                invalid,
                &boundaries,
                || false,
            )
            .is_err());
    }
    assert!(matches!(
        system.assemble_conservative_momentum_bdf1_with_gravity(
            history(),
            controls(1.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.0; 2],
            [0.0, -9.81],
            &boundaries,
            || true,
        ),
        Err(Error::Transport(DensityTransportError::Cancelled))
    ));
    assert_eq!(density[0].density_kg_m3, 2.0);
    assert_eq!(density[1].density_kg_m3, 9.0);
    assert!(velocity.iter().all(|value| *value == 0.0));
}

#[test]
fn spatial_acceleration_uses_transported_mass_and_resolves_angular_motion() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 0.6, 3, 2));
    let system = FluxSystem::new(&mesh);
    let density: Vec<_> = (0..mesh.tris.len())
        .map(|cell| Cell {
            density_kg_m3: 2.0 + 10.0 * cell as f64,
            specific: [],
        })
        .collect();
    let u = [0.7, -0.2];
    let velocity = constant(&mesh, u);
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .map(|face| {
            (face.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 321.0,
                specific: [],
            })
        })
        .collect();
    let exterior: Vec<_> = mesh
        .edges
        .iter()
        .map(|face| (face.tris.1 == usize::MAX).then_some([u; 2]))
        .collect();
    let traction = |_, _| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::TotalTraction {
                edge,
                value: &traction,
            })
        })
        .collect();
    let history = || History {
        density: &density,
        velocity: &velocity,
        inflow_density: &inlet,
        inflow_velocity: &exterior,
    };
    // An affine angular acceleration has zero divergence and zero symmetric
    // gradient. The exact native endpoint is a rigid rotation superposed on
    // the transported free stream, even though rho_new differs in every cell.
    let acceleration = |p: [f64; 2]| [1.3 * (p[1] + 0.2), -1.3 * (p[0] - 0.4)];
    let mut reference_velocity: Option<Vec<f64>> = None;
    for reference in [0.5, 7.0, 1000.0] {
        let samples = std::cell::Cell::new(0);
        let source = |_: usize, point| {
            samples.set(samples.get() + 1);
            acceleration(point)
        };
        let unloaded = system
            .assemble_conservative_momentum_bdf1(
                history(),
                controls(reference),
                MaterialFluxNumerics::default(),
                &|_, _| 0.03,
                &|_| [0.0; 2],
                &boundaries,
                || false,
            )
            .unwrap();
        let loaded = system
            .assemble_conservative_momentum_bdf1_with_acceleration(
                history(),
                controls(reference),
                MaterialFluxNumerics::default(),
                &|_, _| 0.03,
                &|_| [0.0; 2],
                &source,
                &boundaries,
                || false,
            )
            .unwrap();
        assert_eq!(samples.get(), 6 * mesh.tris.len());
        assert!(loaded.density().iter().zip(&density).any(|(a, b)| a != b));
        assert_eq!(loaded.density(), unloaded.density());
        // Independent degree-two edge-midpoint triangle rule. The product
        // of the affine acceleration and BDM test is quadratic: a centroid
        // source sample would miss its covariance, including the odd modes.
        let mut expected_load = vec![0.0; system.n];
        for (cell, nodes) in mesh.tris.iter().enumerate() {
            let vertices = nodes.map(|node| mesh.verts[node]);
            for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                let point = core::array::from_fn(|axis| {
                    f64::midpoint(vertices[a][axis], vertices[b][axis])
                });
                let g = acceleration(point);
                let mass = mesh.areas[cell] * loaded.density()[cell].density_kg_m3 / 3.0;
                for local in 0..6 {
                    let test = eval_basis(&system.bases[cell], local, point);
                    let row = 2 * mesh.tri_edges[cell][local / 2].0 + local % 2;
                    expected_load[row] += mass * (g[0] * test[0] + g[1] * test[1]);
                }
            }
        }
        for ((&with, &without), &expected) in
            loaded.rhs().iter().zip(unloaded.rhs()).zip(&expected_load)
        {
            near(reference * (with - without), expected, 3e-9);
        }
        let accepted = loaded.solve(|| false).unwrap();
        for (cell, nodes) in mesh.tris.iter().enumerate() {
            for (point, _) in tri_quad(nodes.map(|node| mesh.verts[node]), mesh.areas[cell]) {
                let mut actual = [0.0; 2];
                for local in 0..6 {
                    let test = eval_basis(&system.bases[cell], local, point);
                    let row = 2 * mesh.tri_edges[cell][local / 2].0 + local % 2;
                    for axis in 0..2 {
                        actual[axis] += test[axis] * accepted.velocity_m_s()[row];
                    }
                }
                for axis in 0..2 {
                    near(
                        actual[axis],
                        u[axis] + controls(reference).timestep_s * acceleration(point)[axis],
                        3e-9,
                    );
                }
            }
        }
        for &pressure in accepted.pressure_pa() {
            near(pressure, 0.0, 3e-7);
        }
        if let Some(previous) = &reference_velocity {
            for (&a, &b) in accepted.velocity_m_s().iter().zip(previous) {
                near(a, b, 3e-9);
            }
        } else {
            reference_velocity = Some(accepted.velocity_m_s().to_vec());
        }
    }
}

#[test]
fn centrifugal_acceleration_balances_quadratic_pressure_at_rest() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 0.6, 3, 2));
    let system = FluxSystem::new(&mesh);
    let rho = 17.0;
    let omega_squared = 6.25;
    let center = [-0.2, 0.4];
    let density = vec![
        Cell {
            density_kg_m3: rho,
            specific: []
        };
        mesh.tris.len()
    ];
    let velocity = vec![0.0; system.n_u];
    let inlet = vec![None; mesh.edges.len()];
    let exterior = vec![None; mesh.edges.len()];
    let pressure = |point: [f64; 2]| {
        30.0 + 0.5
            * rho
            * omega_squared
            * ((point[0] - center[0]).powi(2) + (point[1] - center[1]).powi(2))
    };
    let traction = |point, normal: [f64; 2]| normal.map(|value| -pressure(point) * value);
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::TotalTraction {
                edge,
                value: &traction,
            })
        })
        .collect();
    let accepted = system
        .assemble_conservative_momentum_bdf1_with_acceleration(
            History {
                density: &density,
                velocity: &velocity,
                inflow_density: &inlet,
                inflow_velocity: &exterior,
            },
            controls(1000.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.03,
            &|_| [0.0; 2],
            &|_, point| core::array::from_fn(|axis| omega_squared * (point[axis] - center[axis])),
            &boundaries,
            || false,
        )
        .unwrap()
        .solve(|| false)
        .unwrap();
    for cell in accepted.density() {
        near(cell.density_kg_m3, rho, 3e-13);
    }
    for &value in accepted.velocity_m_s() {
        near(value, 0.0, 3e-10);
    }
    for (cell, nodes) in mesh.tris.iter().enumerate() {
        let vertices = nodes.map(|node| mesh.verts[node]);
        let expected = [(0, 1), (1, 2), (2, 0)]
            .into_iter()
            .map(|(a, b)| {
                pressure(core::array::from_fn(|axis| {
                    f64::midpoint(vertices[a][axis], vertices[b][axis])
                })) / 3.0
            })
            .sum::<f64>();
        near(accepted.pressure_pa()[cell], expected, 3e-8);
    }
}

#[test]
fn spatial_acceleration_preserves_uniform_paths_and_refuses_late_sources_atomically() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let density = vec![
        Cell {
            density_kg_m3: 2.0,
            specific: [],
        },
        Cell {
            density_kg_m3: 9.0,
            specific: [],
        },
    ];
    let velocity = vec![0.0; system.n_u];
    let inlet = vec![None; mesh.edges.len()];
    let exterior = vec![None; mesh.edges.len()];
    let traction = |_, _| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::TotalTraction {
                edge,
                value: &traction,
            })
        })
        .collect();
    let history = || History {
        density: &density,
        velocity: &velocity,
        inflow_density: &inlet,
        inflow_velocity: &exterior,
    };
    let assemble = |source: &dyn Fn(usize, [f64; 2]) -> [f64; 2]| {
        system.assemble_conservative_momentum_bdf1_with_acceleration(
            history(),
            controls(7.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.2, -0.3],
            source,
            &boundaries,
            || false,
        )
    };
    for gravity in [[0.0; 2], [1.7, -3.2]] {
        let constant = system
            .assemble_conservative_momentum_bdf1_with_gravity(
                history(),
                controls(7.0),
                MaterialFluxNumerics::default(),
                &|_, _| 0.1,
                &|_| [0.2, -0.3],
                gravity,
                &boundaries,
                || false,
            )
            .unwrap();
        let field = assemble(&|_, _| gravity).unwrap();
        assert_eq!(constant.rhs(), field.rhs());
        let a = constant.solve(|| false).unwrap();
        let b = field.solve(|| false).unwrap();
        assert_eq!(a.field().x, b.field().x);
        assert_eq!(a.density(), b.density());
    }
    let samples = std::cell::Cell::new(0);
    let refused = assemble(&|_, _| {
        let count = samples.get() + 1;
        samples.set(count);
        if count == 12 {
            [f64::NAN, 0.0]
        } else {
            [0.4, -0.7]
        }
    });
    let expected_point = tri_quad(mesh.tris[1].map(|node| mesh.verts[node]), mesh.areas[1])[5].0;
    assert!(
        matches!(refused, Err(Error::NonFiniteAcceleration { cell: 1, point_m }) if point_m == expected_point)
    );
    assert_eq!(samples.get(), 12);
    samples.set(0);
    assert!(matches!(
        system.assemble_conservative_momentum_bdf1_with_acceleration(
            history(),
            controls(7.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.2, -0.3],
            &|_, _| {
                samples.set(samples.get() + 1);
                [0.4, -0.7]
            },
            &boundaries,
            || samples.get() == 6,
        ),
        Err(Error::Transport(DensityTransportError::Cancelled))
    ));
    let first = assemble(&|_, _| [0.4, -0.7])
        .unwrap()
        .solve(|| false)
        .unwrap();
    let retry = assemble(&|_, _| [0.4, -0.7])
        .unwrap()
        .solve(|| false)
        .unwrap();
    assert_eq!(first.field().x, retry.field().x);
    assert_eq!(first.density(), retry.density());
    assert_eq!(density[0].density_kg_m3, 2.0);
    assert_eq!(density[1].density_kg_m3, 9.0);
    assert!(velocity.iter().all(|value| *value == 0.0));
}
