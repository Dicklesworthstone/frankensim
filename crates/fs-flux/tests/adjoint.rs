//! Discrete viscosity sensitivities through the actual transposed flow system.

use fs_flux::bdm::{edge_gauss_pub, eval_basis, tri_quad};
use fs_flux::{FluxAdjointError, FluxParams, FluxSystem, TriMesh};
use fs_solid::Mesh2;
use fs_sparse::{Coo, Csr};

fn mesh(n: usize) -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n))
}

fn moments(mesh: &TriMesh, field: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut values = vec![0.0; 2 * mesh.edges.len()];
    for (e, edge) in mesh.edges.iter().enumerate() {
        let (a, b) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
        for (point, weight) in edge_gauss_pub(a, b) {
            let u = field(point);
            let normal = u[0] * edge.normal[0] + u[1] * edge.normal[1];
            let s = ((point[0] - a[0]) * (b[0] - a[0]) + (point[1] - a[1]) * (b[1] - a[1]))
                / edge.len.powi(2)
                - 0.5;
            values[2 * e] += weight * normal / edge.len;
            values[2 * e + 1] += weight * normal * s / edge.len;
        }
    }
    values
}

fn objective(system: &FluxSystem<'_>, weight: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut gradient = vec![0.0; system.n];
    for cell in 0..system.mesh.tris.len() {
        let triangle = core::array::from_fn(|k| system.mesh.verts[system.mesh.tris[cell][k]]);
        for (point, area) in tri_quad(triangle, system.mesh.areas[cell]) {
            let w = weight(point);
            for i in 0..6 {
                let basis = eval_basis(&system.bases[cell], i, point);
                let dof = 2 * system.mesh.tri_edges[cell][i / 2].0 + i % 2;
                gradient[dof] += area * (w[0] * basis[0] + w[1] * basis[1]);
            }
        }
    }
    gradient
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn scaled_viscosity_operator(
    system: &FluxSystem<'_>,
    params: FluxParams,
    scale: f64,
    advect: &[f64],
) -> (Csr, Vec<f64>) {
    system.assemble_viscous_transient(
        params,
        &|cell, x| {
            scale
                * (0.2 + x[0] * x[0])
                * if system.mesh.centroids[cell][1] < 0.5 {
                    1.0
                } else {
                    20.0
                }
        },
        &|x| [1.0 + x[1] * x[1], 0.2 * x[0]],
        &|_| [0.0; 2],
        Some(advect),
        7.0,
        Some(advect),
    )
}

#[test]
fn adjoint_transposes_constraints_and_nonsymmetric_convection() {
    let mesh = mesh(4);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-11,
        ..FluxParams::default()
    };
    let advect = moments(&mesh, &|x| [0.4 * x[1], -0.2 * x[0]]);
    let (matrix, _) = scaled_viscosity_operator(&system, params, 0.7, &advect);
    let exact: Vec<_> = (0..system.n)
        .map(|i| ((17 * i % 31) as f64 - 15.0) / 15.0)
        .collect();
    // Independent scatter from ORIGINAL rows; do not use the production
    // transpose routine to manufacture its own test right-hand side.
    let mut gradient = vec![0.0; system.n];
    for (row, &multiplier) in exact.iter().enumerate() {
        let (columns, coefficients) = matrix.row(row);
        for (&column, &coefficient) in columns.iter().zip(coefficients) {
            gradient[column] += coefficient * multiplier;
        }
    }
    let solution = system
        .solve_adjoint_checked(&matrix, &gradient, params)
        .unwrap();
    let error = solution
        .multipliers
        .iter()
        .zip(&exact)
        .fold(0.0_f64, |error, (a, b)| error.max((a - b).abs()));
    assert!(error < 2e-9, "adjoint maximum error {error:e}");
    let mut wrong = vec![0.0; system.n];
    matrix.spmv(&exact, &mut wrong);
    assert!(
        wrong
            .iter()
            .zip(&gradient)
            .any(|(a, b)| (a - b).abs() > 1.0)
    );
}

#[test]
fn sparse_adjoint_matches_viscosity_gradient_above_dense_cutoff() {
    let mesh = mesh(16);
    let system = FluxSystem::new(&mesh);
    assert!(system.n > 1500);
    let params = FluxParams {
        tol: 1e-11,
        max_iters: 5000,
        ..FluxParams::default()
    };
    let scale = 0.7;
    let advect = moments(&mesh, &|x| [0.4 * x[1], -0.2 * x[0]]);
    let gradient = objective(&system, &|x| [1.0 + x[1] * x[1], 0.2 * x[0]]);
    let solve = |scale| {
        let (matrix, rhs) = scaled_viscosity_operator(&system, params, scale, &advect);
        let solution = system.solve_linear(&matrix, &rhs, params);
        assert!(solution.converged(params.tol));
        (dot(&gradient, &solution.x), solution.x, matrix, rhs)
    };
    let (_, primal, matrix, rhs) = solve(scale);
    let adjoint = system
        .solve_adjoint_checked(&matrix, &gradient, params)
        .unwrap();
    assert!(adjoint.iters > 0 && adjoint.iters <= params.max_iters);
    let (twice, twice_rhs) = scaled_viscosity_operator(&system, params, 2.0 * scale, &advect);
    let mut applied = vec![0.0; system.n];
    let mut applied_twice = vec![0.0; system.n];
    matrix.spmv(&primal, &mut applied);
    twice.spmv(&primal, &mut applied_twice);
    // The harmonic interface weights are unchanged by a common viscosity
    // multiplier. All viscous terms scale linearly, while mass/convection,
    // pressure and strong constraints stay fixed: this derivative is exact.
    let sensitivity: f64 = (0..system.n)
        .map(|i| {
            adjoint.multipliers[i] * (twice_rhs[i] - rhs[i] - applied_twice[i] + applied[i]) / scale
        })
        .sum();
    let h = 1e-3 * scale;
    let finite_difference = (solve(scale + h).0 - solve(scale - h).0) / (2.0 * h);
    let relative = (sensitivity - finite_difference).abs() / finite_difference.abs();
    eprintln!(
        "sparse viscosity adjoint: n={}, iterations={}, residual={:e}, scaled={:e}, dJ={sensitivity:e}, FD={finite_difference:e}, relative={relative:e}",
        system.n, adjoint.iters, adjoint.rel_residual, adjoint.scaled_rel_residual
    );
    assert!(finite_difference.abs() > 1e-7);
    assert!(relative < 2e-5);
}

#[test]
fn adjoint_accounts_for_viscosity_dependent_boundary_loads() {
    let mesh = mesh(4);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-11,
        ..FluxParams::default()
    };
    let scale = 0.7;
    let assemble = |scale| {
        system.assemble_variable_viscosity(
            params,
            &|_, _| scale,
            &|_| [0.0; 2],
            &|x| [scale * x[1], 0.0],
            None,
        )
    };
    let gradient = objective(&system, &|_| [1.0, 0.0]);
    let (matrix, rhs) = assemble(scale);
    let primal = system.solve_linear(&matrix, &rhs, params);
    let (error, divergence) = system.velocity_error(&primal.x, &|x| [scale * x[1], 0.0]);
    assert!(error < 1e-11 && divergence < 1e-11);
    assert!((dot(&gradient, &primal.x) - scale / 2.0).abs() < 1e-11);
    let adjoint = system
        .solve_adjoint_checked(&matrix, &gradient, params)
        .unwrap();
    let (twice, twice_rhs) = assemble(2.0 * scale);
    let mut applied = vec![0.0; system.n];
    let mut applied_twice = vec![0.0; system.n];
    matrix.spmv(&primal.x, &mut applied);
    twice.spmv(&primal.x, &mut applied_twice);
    // b(theta)=theta*b_normal + theta^2*b_SIP, hence b'(theta)=b(2theta)/(2theta).
    let sensitivity: f64 = (0..system.n)
        .map(|i| {
            adjoint.multipliers[i] * (0.5 * twice_rhs[i] - applied_twice[i] + applied[i]) / scale
        })
        .sum();
    assert!(
        (sensitivity - 0.5).abs() < 2e-10,
        "dJ/dtheta={sensitivity:e}"
    );
}

#[test]
fn adjoint_refuses_invalid_inputs_and_honors_sparse_budgets() {
    let mesh = mesh(16);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-11,
        max_iters: 3,
        ..FluxParams::default()
    };
    let (matrix, _) = system.assemble(params, &|_| [0.0; 2], &|_| [0.0; 2], None);
    let gradient: Vec<_> = (0..system.n)
        .map(|i| ((i * 13 % 29) as f64 - 14.0) / 14.0)
        .collect();
    assert!(matches!(
        system.solve_adjoint_checked(&matrix, &gradient, params),
        Err(FluxAdjointError::NonConvergence { iterations: 3, .. })
    ));
    let zero = system
        .solve_adjoint_checked(
            &matrix,
            &vec![0.0; system.n],
            FluxParams {
                max_iters: 0,
                ..params
            },
        )
        .unwrap();
    assert_eq!(zero.iters, 0);
    assert!(zero.multipliers.iter().all(|x| *x == 0.0));
    assert!(matches!(
        system.solve_adjoint_checked(&matrix, &gradient[..system.n - 1], params),
        Err(FluxAdjointError::Invalid(_))
    ));
    assert!(matches!(
        system.solve_adjoint_checked(
            &matrix,
            &gradient,
            FluxParams {
                tol: f64::NAN,
                ..params
            }
        ),
        Err(FluxAdjointError::Invalid(_))
    ));
    let mut bad = Coo::new(system.n, system.n);
    bad.push(0, 0, f64::NAN);
    assert!(matches!(
        system.solve_adjoint_checked(&bad.assemble(), &gradient, params),
        Err(FluxAdjointError::Invalid(_))
    ));
    let small_mesh = self::mesh(2);
    let small = FluxSystem::new(&small_mesh);
    assert!(matches!(
        small.solve_adjoint_checked(
            &Coo::new(small.n, small.n).assemble(),
            &vec![1.0; small.n],
            params
        ),
        Err(FluxAdjointError::SingularSystem)
    ));
}
