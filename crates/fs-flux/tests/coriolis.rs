//! G0/G1/G3 planar Coriolis physics on the original BDM momentum owner.
#![cfg(feature = "continuum")]

use fs_flux::ale::rigid::{RigidAleError, RigidAleFlow, RigidAleMotion};
use fs_flux::ale::translating::{TranslatingAleError, TranslatingAleFlow};
use fs_flux::{
    FluxParams, FluxSystem, MaterialFluxNumerics, MaterialFluxSystem, MaterialStokesBoundary,
    StokesBoundary, StokesBoundaryError, StokesPressureReference, TriMesh,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_solid::Mesh2;

fn mesh() -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2))
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (tolerance {tolerance:.3e})"
    );
}

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

fn boundary_edges(mesh: &TriMesh) -> Vec<usize> {
    (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect()
}

fn no_stress(_: [f64; 2], _: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}

fn open_boundaries(mesh: &TriMesh) -> Vec<StokesBoundary<'_>> {
    boundary_edges(mesh)
        .into_iter()
        .map(|edge| StokesBoundary::TotalTraction {
            edge,
            value: &no_stress,
        })
        .collect()
}

fn backward_euler(velocity: [f64; 2], parameter: f64, dt: f64) -> [f64; 2] {
    let angle = parameter * dt;
    let denominator = 1.0 + angle * angle;
    [
        (velocity[0] + angle * velocity[1]) / denominator,
        (velocity[1] - angle * velocity[0]) / denominator,
    ]
}

/// The unconstrained uniform mode has u' = -f J u. Its backward-Euler
/// amplitude falls by 1/sqrt(1+(f dt)^2); that numerical damping is not heating.
#[test]
fn g1_uniform_inertial_oscillation_has_correct_sign_phase_and_bdf1_damping() {
    let mesh = mesh();
    let boundaries = open_boundaries(&mesh);
    for parameter in [-3.5, 3.5] {
        let system = FluxSystem::new(&mesh)
            .with_coriolis_parameter(parameter)
            .unwrap();
        let mut exact = [0.7, -1.1];
        let mut previous = uniform_moments(&mesh, exact);
        for dt in [0.03, 0.11, 0.06, 0.08] {
            let old_energy = 0.5 * (exact[0] * exact[0] + exact[1] * exact[1]);
            exact = backward_euler(exact, parameter, dt);
            let snapshot = previous.clone();
            let solved = system
                .bdf1_stokes_mixed_step(
                    FluxParams::default(),
                    &|_| [0.0; 2],
                    &boundaries,
                    &previous,
                    dt,
                )
                .unwrap();
            let (error, divergence) = system.velocity_error(&solved.field.x, &|_| exact);
            assert!(error < 3e-11 && divergence < 3e-11, "{error}, {divergence}");
            for &pressure in &solved.pressure {
                near(pressure, 0.0, 3e-10);
            }
            near(
                system.kinetic_energy(&solved.field.x),
                old_energy / (1.0 + (parameter * dt).powi(2)),
                3e-11,
            );
            assert_eq!(previous, snapshot);
            previous.copy_from_slice(&solved.field.x[..system.n_u]);
        }
    }
}

/// r=1+2x+3y has total mass 7/2 on the unit square. For u=(1,2),
/// v=(x,y), the independent exact integral of r*v.Ju is -11/6.
#[test]
fn g0_g1_density_metric_preserves_skew_work_and_uniform_angular_response() {
    let mesh = mesh();
    let plain = FluxSystem::new(&mesh);
    let parameter = 1.7;
    let rotating = FluxSystem::new(&mesh)
        .with_coriolis_parameter(parameter)
        .unwrap();
    let params = FluxParams {
        nu: 0.2,
        ..FluxParams::default()
    };
    let boundaries = open_boundaries(&mesh);
    let density = |_: usize, point: [f64; 2]| 1.0 + 2.0 * point[0] + 3.0 * point[1];
    let previous = uniform_moments(&mesh, [1.0, 2.0]);
    let dt = 0.2;
    let base = plain
        .assemble_stokes_mixed_bdf1_variable_density(
            params,
            &density,
            &|_, _| params.nu,
            &|_| [0.0; 2],
            &boundaries,
            &previous,
            dt,
        )
        .unwrap();
    let source = rotating
        .assemble_stokes_mixed_bdf1_variable_density(
            params,
            &density,
            &|_, _| params.nu,
            &|_| [0.0; 2],
            &boundaries,
            &previous,
            dt,
        )
        .unwrap();
    assert_eq!(base.rhs(), source.rhs());
    let n = rotating.n_u;
    let mut coupling = vec![vec![0.0; n]; n];
    for (row, entries) in coupling.iter_mut().enumerate() {
        for (matrix, sign) in [(source.matrix(), 1.0), (base.matrix(), -1.0)] {
            let (columns, values) = matrix.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                if column < n {
                    entries[column] += sign * value;
                }
            }
        }
    }
    for i in 0..n {
        for j in 0..n {
            near(coupling[i][j] + coupling[j][i], 0.0, 2e-12);
        }
    }
    let arbitrary: Vec<_> = (0..n)
        .map(|i| ((17 * i % 23) as f64 - 11.0) / 7.0)
        .collect();
    let work: f64 = coupling
        .iter()
        .enumerate()
        .map(|(i, row)| arbitrary[i] * row.iter().zip(&arbitrary).map(|(c, u)| c * u).sum::<f64>())
        .sum();
    near(work, 0.0, 3e-11);
    // (x,y).n is constant along each straight edge; its odd moment is zero.
    let radial: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|edge| {
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
            [
                0.5 * ((a[0] + b[0]) * edge.normal[0] + (a[1] + b[1]) * edge.normal[1]),
                0.0,
            ]
        })
        .collect();
    let mixed_work: f64 = coupling
        .iter()
        .enumerate()
        .map(|(i, row)| radial[i] * row.iter().zip(&previous).map(|(c, u)| c * u).sum::<f64>())
        .sum();
    near(mixed_work, -parameter * 11.0 / 6.0, 3e-11);

    let exact = backward_euler([1.0, 2.0], parameter, dt);
    let solved = source.solve().unwrap();
    let (error, divergence) = rotating.velocity_error(&solved.field.x, &|_| exact);
    assert!(error < 3e-11 && divergence < 3e-11, "{error}, {divergence}");
    near(
        source.kinetic_energy(&solved.field.x).unwrap(),
        1.75 * (exact[0] * exact[0] + exact[1] * exact[1]),
        3e-11,
    );
    let sparse = rotating.solve_sparse(source.matrix(), source.rhs(), params);
    assert!(sparse.converged(params.tol));
    let (error, divergence) = rotating.velocity_error(&sparse.x, &|_| exact);
    assert!(error < 3e-8 && divergence < 3e-8, "{error}, {divergence}");
}

/// Constant u balances grad(p)=-f J u. The P0 pressure is the exact cell
/// average of this linear geostrophic pressure; boundary reaction uses the
/// original momentum rows, including Coriolis on constrained normal moments.
#[test]
fn g1_geostrophic_pressure_and_boundary_reaction_balance_without_body_force() {
    let mesh = mesh();
    let parameter = 6.0;
    let system = FluxSystem::new(&mesh)
        .with_coriolis_parameter(parameter)
        .unwrap();
    let velocity = [0.4, -0.7];
    let prescribed = |_| velocity;
    let edges = boundary_edges(&mesh);
    let boundaries: Vec<_> = edges
        .iter()
        .map(|&edge| StokesBoundary::Velocity {
            edge,
            value: &prescribed,
        })
        .collect();
    let solved = system
        .assemble_stokes_mixed(FluxParams::default(), &|_| [0.0; 2], &boundaries)
        .unwrap()
        .solve()
        .unwrap();
    let (error, divergence) = system.velocity_error(&solved.field.x, &|_| velocity);
    assert!(error < 3e-11 && divergence < 3e-11, "{error}, {divergence}");
    assert_eq!(solved.pressure_reference, StokesPressureReference::MeanZero);
    for (pressure, point) in solved.pressure.iter().zip(&mesh.centroids) {
        let exact = parameter * (velocity[1] * (point[0] - 0.5) - velocity[0] * (point[1] - 0.5));
        near(*pressure, exact, 3e-10);
    }
    let load = solved
        .reactions
        .resultant(&edges, [0.5, 0.5], &[Some(0.0)])
        .unwrap();
    near(load.stress.force[0], parameter * velocity[1], 3e-10);
    near(load.stress.force[1], -parameter * velocity[0], 3e-10);
    near(load.stress.rigid_power(velocity, 0.0).unwrap(), 0.0, 3e-10);
}

#[test]
fn g3_material_density_changes_physical_energy_without_changing_inertial_frequency() {
    let mesh = mesh();
    let parameter = -2.3;
    let system = FluxSystem::new(&mesh)
        .with_coriolis_parameter(parameter)
        .unwrap();
    let previous = uniform_moments(&mesh, [0.3, 0.8]);
    let dt = 0.12;
    let exact = backward_euler([0.3, 0.8], parameter, dt);
    let boundaries: Vec<_> = boundary_edges(&mesh)
        .into_iter()
        .map(|edge| MaterialStokesBoundary::TotalTraction {
            edge,
            value: &no_stress,
        })
        .collect();
    for density in [0.7, 5.0] {
        let material = FluidTransportState::try_new(density, 1000.0, 0.5, 0.2 * density).unwrap();
        let flow =
            MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
        let solved = flow
            .bdf1_stokes_mixed_step(&|_| [0.0; 2], &boundaries, &previous, dt)
            .unwrap();
        let (error, divergence) = system.velocity_error(&solved.discrete().x, &|_| exact);
        assert!(error < 3e-11 && divergence < 3e-11, "{error}, {divergence}");
        near(
            solved.kinetic_energy_j_per_m(),
            0.5 * density * (exact[0] * exact[0] + exact[1] * exact[1]),
            3e-11,
        );
    }
}

#[test]
fn g0_zero_configuration_preserves_original_assembly_and_nonfinite_values_refuse() {
    let mesh = mesh();
    let plain = FluxSystem::new(&mesh);
    assert_eq!(plain.coriolis_parameter(), 0.0);
    let params = FluxParams::default();
    let velocity = |_| [0.4, -0.2];
    let force = |point: [f64; 2]| [point[0] * point[1], -0.7 * point[0]];
    let previous = uniform_moments(&mesh, velocity([0.0; 2]));
    let (base, rhs) = plain.assemble_transient(
        params,
        &force,
        &velocity,
        Some(&previous),
        4.0,
        Some(&previous),
    );
    for parameter in [0.0, -0.0] {
        let configured = FluxSystem::new(&mesh)
            .with_coriolis_parameter(parameter)
            .unwrap();
        let (matrix, load) = configured.assemble_transient(
            params,
            &force,
            &velocity,
            Some(&previous),
            4.0,
            Some(&previous),
        );
        assert_eq!(configured.coriolis_parameter(), parameter);
        for row in 0..plain.n {
            let (columns, values) = base.row(row);
            let (actual_columns, actual_values) = matrix.row(row);
            assert_eq!(columns, actual_columns);
            assert!(
                values
                    .iter()
                    .zip(actual_values)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
        }
        assert!(
            rhs.iter()
                .zip(load)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            FluxSystem::new(&mesh).with_coriolis_parameter(invalid),
            Err(StokesBoundaryError::InvalidInput(_))
        ));
    }
}

#[test]
fn g0_world_ale_owners_refuse_relative_coriolis_without_changing_zero_behavior() {
    let mesh = mesh();
    let material = FluidTransportState::try_new(1.0, 1000.0, 0.5, 0.2).unwrap();
    let previous = uniform_moments(&mesh, [0.3, -0.1]);
    for parameter in [0.0, -0.0, 2.0, -2.0] {
        let system = FluxSystem::new(&mesh)
            .with_coriolis_parameter(parameter)
            .unwrap();
        let translating = TranslatingAleFlow::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            [0.2, -0.4],
            &previous,
        );
        let rigid = RigidAleFlow::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            [0.5, 0.5],
            RigidAleMotion::ZERO,
            &previous,
        );
        if parameter == 0.0 {
            assert!(translating.is_ok() && rigid.is_ok());
        } else {
            assert!(matches!(translating, Err(TranslatingAleError::Invalid(_))));
            assert!(matches!(rigid, Err(RigidAleError::Invalid(_))));
        }
    }
}
