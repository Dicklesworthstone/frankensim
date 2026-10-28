//! G1/G3: independent residual differences and nonlinear re-solve gradients.
#![cfg(feature = "continuum")]

use fs_flux::bdm::{edge_gauss_pub, eval_basis, tri_quad};
use fs_flux::ns::linearization::{FluxLinearizationError, FluxViscosityTangent};
use fs_flux::{FluxParams, FluxSystem, TriMesh};
use fs_solid::Mesh2;
use std::cell::Cell;

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

fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |n, v| n.hypot(*v))
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn objective(system: &FluxSystem<'_>) -> Vec<f64> {
    let mut gradient = vec![0.0; system.n];
    for cell in 0..system.mesh.tris.len() {
        let triangle = core::array::from_fn(|k| system.mesh.verts[system.mesh.tris[cell][k]]);
        for (point, area) in tri_quad(triangle, system.mesh.areas[cell]) {
            let weight = [1.0 + point[1] * point[1], 0.3 * point[0]];
            for i in 0..6 {
                let basis = eval_basis(&system.bases[cell], i, point);
                let dof = 2 * system.mesh.tri_edges[cell][i / 2].0 + i % 2;
                gradient[dof] += area * (weight[0] * basis[0] + weight[1] * basis[1]);
            }
        }
    }
    gradient
}

#[test]
fn constitutive_and_convective_jacobian_matches_full_residual_differences() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        sigma: 20.0,
        ..FluxParams::default()
    };
    let mut state = moments(&mesh, &|x| [0.8 + 0.3 * x[1], 0.21 + 0.1 * x[0]]);
    for (i, value) in state.iter_mut().enumerate() {
        *value += 0.003 * ((i * 17 % 29) as f64 - 14.0);
    }
    state.extend((0..mesh.tris.len()).map(|i| (i as f64 * 0.31).sin()));
    let direction: Vec<_> = (0..system.n)
        .map(|i| ((i * 13 % 43) as f64 - 21.0) / 23.0)
        .collect();
    let previous = moments(&mesh, &|x| [0.3 * x[1], 0.1 * x[0]]);
    let force = |x: [f64; 2]| [1.0 + x[1] * x[1], 0.2 * x[0]];
    let boundary = |x: [f64; 2]| [0.6 + 0.2 * x[1], -0.3 + 0.1 * x[0]];
    for exponent in [0.4, 1.7] {
        // Point variation AND a 100:1 cell contrast exercise the derivatives
        // of one-sided harmonic weights, rather than only a common scale.
        let law = |cell: usize, point: [f64; 2], rate: f64| {
            let scale = (0.5 + point[0]) * if cell % 3 == 0 { 100.0 } else { 1.0 };
            let nu = scale * (1.0 + rate * rate).powf(0.5 * (exponent - 1.0));
            FluxViscosityTangent {
                nu,
                derivative: nu * (exponent - 1.0) * rate / (1.0 + rate * rate),
            }
        };
        for convection in [false, true] {
            let linear = system
                .linearize_rheology(
                    params,
                    &law,
                    &force,
                    &boundary,
                    &state,
                    convection,
                    3.5,
                    Some(&previous),
                )
                .unwrap();
            assert_eq!(linear.nonsmooth_upwind_points, 0);
            let mut applied = vec![0.0; system.n];
            linear.jacobian.spmv(&direction, &mut applied);
            // This oracle uses the original primal assembler at perturbed
            // states; it never evaluates the new derivative callback there.
            let residual = |x: &[f64]| {
                let rates = system.shear_rates(x);
                let (operator, rhs) = system.assemble_viscous_transient(
                    params,
                    &|cell, point| law(cell, point, rates[cell]).nu,
                    &force,
                    &boundary,
                    convection.then_some(&x[..system.n_u]),
                    3.5,
                    Some(&previous),
                );
                let mut value = vec![0.0; system.n];
                operator.spmv(x, &mut value);
                for (v, b) in value.iter_mut().zip(rhs) {
                    *v -= b;
                }
                (value, operator)
            };
            let (original, frozen) = residual(&state);
            assert_eq!(original, linear.residual);
            let mut frozen_applied = vec![0.0; system.n];
            frozen.spmv(&direction, &mut frozen_applied);
            let missing: Vec<_> = applied
                .iter()
                .zip(&frozen_applied)
                .map(|(a, b)| a - b)
                .collect();
            assert!(
                norm(&missing) / norm(&applied) > 1e-3,
                "negative control must expose missing derivatives"
            );
            let mut errors = Vec::new();
            for h in [2e-5, 1e-5] {
                let plus: Vec<_> = state
                    .iter()
                    .zip(&direction)
                    .map(|(x, d)| x + h * d)
                    .collect();
                let minus: Vec<_> = state
                    .iter()
                    .zip(&direction)
                    .map(|(x, d)| x - h * d)
                    .collect();
                let rp = residual(&plus).0;
                let rm = residual(&minus).0;
                let error: Vec<_> = rp
                    .iter()
                    .zip(&rm)
                    .zip(&applied)
                    .map(|((p, m), a)| (p - m) / (2.0 * h) - a)
                    .collect();
                errors.push(norm(&error) / norm(&applied));
            }
            eprintln!(
                "rheology Jacobian n={exponent}, convection={convection}, FD errors={errors:?}, frozen defect={:e}",
                norm(&missing) / norm(&applied)
            );
            assert!(errors[1] < 2e-7, "{errors:?}");
            assert!(
                errors[1] < 0.4 * errors[0] + 2e-10,
                "centered differences must converge: {errors:?}"
            );
        }
    }
}

fn carreau(scale: f64, rate: f64) -> FluxViscosityTangent {
    let nu = scale * 0.2 * (1.0 + 9.0 * rate * rate).powf(-0.25);
    FluxViscosityTangent {
        nu,
        derivative: -4.5 * rate * nu / (1.0 + 9.0 * rate * rate),
    }
}

#[test]
fn nonlinear_bdf1_viscosity_adjoint_matches_independent_resolves() {
    let mesh = mesh(4);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-12,
        ..FluxParams::default()
    };
    let boundary = |x: [f64; 2]| [0.4 + 0.3 * x[1], 0.15 + 0.1 * x[0]];
    let force = |x: [f64; 2]| [1.0 + x[1], 0.2 * x[0]];
    let previous = moments(&mesh, &boundary);
    let gradient = objective(&system);
    let dt = 0.2;
    let solve = |scale| {
        system
            .bdf1_rheology_step(
                params,
                &|_, _, rate| carreau(scale, rate).nu,
                &force,
                &boundary,
                &previous,
                dt,
                100,
                5e-12,
            )
            .unwrap()
    };
    let endpoint = solve(1.0);
    let linear = system
        .linearize_rheology(
            params,
            &|_, _, rate| carreau(1.0, rate),
            &force,
            &boundary,
            &endpoint.solution.x,
            true,
            1.0 / dt,
            Some(&previous),
        )
        .unwrap();
    assert_eq!(linear.nonsmooth_upwind_points, 0);
    let adjoint = system
        .solve_adjoint_checked(&linear.jacobian, &gradient, params)
        .unwrap();
    let rates = system.shear_rates(&endpoint.solution.x);
    let frozen = |scale| {
        system.assemble_viscous_transient(
            params,
            &|cell, _| carreau(scale, rates[cell]).nu,
            &force,
            &boundary,
            Some(&endpoint.solution.x[..system.n_u]),
            1.0 / dt,
            Some(&previous),
        )
    };
    let (base, load) = frozen(1.0);
    let (twice, twice_load) = frozen(2.0);
    let mut base_applied = vec![0.0; system.n];
    let mut twice_applied = vec![0.0; system.n];
    base.spmv(&endpoint.solution.x, &mut base_applied);
    twice.spmv(&endpoint.solution.x, &mut twice_applied);
    // At fixed state, all viscous terms AND their boundary loads are linear
    // in this common material scale; all other physical inputs stay fixed.
    let parameter_rhs: Vec<_> = (0..system.n)
        .map(|i| twice_load[i] - load[i] - twice_applied[i] + base_applied[i])
        .collect();
    let sensitivity = dot(&adjoint.multipliers, &parameter_rhs);
    let h = 1e-4;
    let difference = (dot(&gradient, &solve(1.0 + h).solution.x)
        - dot(&gradient, &solve(1.0 - h).solution.x))
        / (2.0 * h);
    let relative = (sensitivity - difference).abs() / difference.abs();
    let frozen_adjoint = system
        .solve_adjoint_checked(&base, &gradient, params)
        .unwrap();
    let frozen_sensitivity = dot(&frozen_adjoint.multipliers, &parameter_rhs);
    let frozen_error = (frozen_sensitivity - difference).abs() / difference.abs();
    eprintln!(
        "nonlinear BDF1 adjoint: iterations={}, dJ={sensitivity:e}, FD={difference:e}, relative={relative:e}, frozen relative={frozen_error:e}",
        endpoint.iterations
    );
    assert!(difference.abs() > 1e-6);
    assert!(relative < 2e-6);
    assert!(frozen_error > 1e-3);
}

#[test]
fn zero_strain_and_newtonian_reduction_preserve_actual_operator() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let zero = vec![0.0; system.n];
    let law = |cell: usize, point: [f64; 2], _: f64| FluxViscosityTangent {
        nu: 0.3 + point[0] + cell as f64 * 0.1,
        derivative: 0.0,
    };
    let linear = system
        .linearize_rheology(
            params,
            &law,
            &|_| [0.0; 2],
            &|_| [0.0; 2],
            &zero,
            false,
            0.0,
            None,
        )
        .unwrap();
    let (base, _) = system.assemble_variable_viscosity(
        params,
        &|cell, point| law(cell, point, 0.0).nu,
        &|_| [0.0; 2],
        &|_| [0.0; 2],
        None,
    );
    for row in 0..system.n {
        let (cols, values) = linear.jacobian.row(row);
        let (base_cols, base_values) = base.row(row);
        assert_eq!(cols, base_cols);
        assert_eq!(
            values.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            base_values.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
    }
    assert!(linear.residual.iter().all(|x| *x == 0.0));
    let switching = system
        .linearize_rheology(
            params,
            &law,
            &|_| [0.0; 2],
            &|_| [1.0, 0.0],
            &zero,
            true,
            0.0,
            None,
        )
        .unwrap();
    assert!(
        switching.nonsmooth_upwind_points > 0,
        "zero-flux upwind ambiguity must remain visible"
    );
}

#[test]
fn invalid_states_and_constitutive_tangents_are_refused_without_result() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let zero = vec![0.0; system.n];
    let calls = Cell::new(0);
    let law = |_, _, _| {
        calls.set(calls.get() + 1);
        FluxViscosityTangent {
            nu: 1.0,
            derivative: 0.0,
        }
    };
    assert!(
        system
            .linearize_rheology(
                params,
                &law,
                &|_| [0.0; 2],
                &|_| [0.0; 2],
                &[],
                false,
                0.0,
                None
            )
            .is_err()
    );
    assert_eq!(calls.get(), 0);
    assert!(
        system
            .linearize_rheology(
                params,
                &law,
                &|_| [0.0; 2],
                &|_| [0.0; 2],
                &zero,
                false,
                1.0,
                None
            )
            .is_err()
    );
    assert_eq!(calls.get(), 0);
    for tangent in [
        FluxViscosityTangent {
            nu: 0.0,
            derivative: 0.0,
        },
        FluxViscosityTangent {
            nu: f64::INFINITY,
            derivative: 0.0,
        },
        FluxViscosityTangent {
            nu: 1.0,
            derivative: f64::NAN,
        },
        FluxViscosityTangent {
            nu: 1.0,
            derivative: 0.1,
        },
    ] {
        let result = system.linearize_rheology(
            params,
            &|_, _, _| tangent,
            &|_| [0.0; 2],
            &|_| [0.0; 2],
            &zero,
            false,
            0.0,
            None,
        );
        assert!(matches!(
            result,
            Err(FluxLinearizationError::Constitutive { .. })
        ));
    }
    assert!(
        system
            .linearize_rheology(
                params,
                &law,
                &|_| [f64::NAN, 0.0],
                &|_| [0.0; 2],
                &zero,
                false,
                0.0,
                None
            )
            .is_err()
    );
    assert!(
        system
            .linearize_rheology(
                params,
                &law,
                &|_| [0.0; 2],
                &|_| [f64::INFINITY, 0.0],
                &zero,
                false,
                0.0,
                None
            )
            .is_err()
    );
}
