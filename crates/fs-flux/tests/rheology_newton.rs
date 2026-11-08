//! G1/G3: stiff constitutive solves, direct residual proof and budget refusal.
#![cfg(feature = "continuum")]

use fs_flux::bdm::{cell_basis, eval_basis, tri_quad};
use fs_flux::ns::linearization::FluxViscosityTangent;
use fs_flux::ns::newton::{FluxNewtonControls, FluxNewtonError};
use fs_flux::{FluxParams, FluxSystem, TriMesh};
use fs_solid::Mesh2;

fn mesh(n: usize) -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n))
}

fn thickening(_: usize, _: [f64; 2], rate: f64) -> FluxViscosityTangent {
    let nu = 0.1 * (1.0 + rate * rate).powi(2);
    FluxViscosityTangent {
        nu,
        derivative: 4.0 * rate * nu / (1.0 + rate * rate),
    }
}

// Exact BDM moments of u = [shear*y, 0]. For an affine normal trace,
// its mean is the endpoint average and its centered first moment is
// the endpoint difference divided by twelve; no flow solve initializes it.
fn couette_moments(mesh: &TriMesh, shear: f64) -> Vec<f64> {
    mesh.edges
        .iter()
        .flat_map(|edge| {
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
            [
                shear * edge.normal[0] * (a[1] + b[1]) / 2.0,
                shear * edge.normal[0] * (b[1] - a[1]) / 12.0,
            ]
        })
        .collect()
}

#[test]
fn bdf2_newton_quadratic_shear_matches_exact_flow_and_independent_mass_residual() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-12,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls {
        residual_tolerance: 1e-10,
        ..FluxNewtonControls::default()
    };
    let amplitude = |t: f64| 1.0 + t + t * t;
    let law = |_: usize, _: [f64; 2], rate: f64| FluxViscosityTangent {
        nu: 0.2 * (1.0 + rate * rate),
        derivative: 0.4 * rate,
    };
    let previous_time = 0.4;
    let previous_dt = 0.3;
    let older = couette_moments(&mesh, amplitude(previous_time - previous_dt));
    let previous = couette_moments(&mesh, amplitude(previous_time));
    let saved_histories = (older.clone(), previous.clone());
    for timestep in [0.15, 0.3, 0.6] {
        let time = previous_time + timestep;
        let shear = amplitude(time);
        let velocity = |p: [f64; 2]| [shear * p[1], 0.0];
        // u_t = [(1+2t)y,0]; convection and the divergence of the uniform
        // shear-dependent viscous stress vanish. BDF2 differentiates this
        // quadratic time dependence exactly even when the intervals differ.
        let force = |p: [f64; 2]| [(1.0 + 2.0 * time) * p[1], 0.0];
        let result = system
            .bdf2_rheology_newton_variable_step(
                params,
                &law,
                &force,
                &velocity,
                &previous,
                Some((&older, previous_dt)),
                timestep,
                controls,
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&result.solution.x, &velocity);
        assert!(error < 2e-9 && divergence < 2e-9);
        let rates = system.shear_rates(&result.solution.x);
        assert!(rates.iter().all(|rate| (rate - shear).abs() < 2e-8));
        assert!(
            system
                .pressure(&result.solution.x)
                .iter()
                .all(|p| p.abs() < 1e-7)
        );

        // Independently integrate M * du/dt using the three raw BDF2
        // coefficients. This never calls the production history helper or
        // substitutes its normalized history into a transient assembler.
        let ratio = timestep / previous_dt;
        let a0 = (1.0 + 2.0 * ratio) / (1.0 + ratio);
        let a1 = -(1.0 + ratio);
        let a2 = ratio * ratio / (1.0 + ratio);
        let derivative: Vec<_> = (0..system.n_u)
            .map(|i| (a0 * result.solution.x[i] + a1 * previous[i] + a2 * older[i]) / timestep)
            .collect();
        let (operator, rhs) = system.assemble_variable_viscosity(
            params,
            &|cell, point| law(cell, point, rates[cell]).nu,
            &force,
            &velocity,
            Some(&result.solution.x[..system.n_u]),
        );
        let mut residual = vec![0.0; system.n];
        operator.spmv(&result.solution.x, &mut residual);
        for (value, load) in residual.iter_mut().zip(&rhs) {
            *value -= load;
        }
        for (cell, triangle) in mesh.tris.iter().enumerate() {
            let basis = cell_basis(&mesh, cell);
            let vertices = triangle.map(|vertex| mesh.verts[vertex]);
            for (point, weight) in tri_quad(vertices, mesh.areas[cell]) {
                let values: [_; 6] = core::array::from_fn(|i| eval_basis(&basis, i, point));
                let mut acceleration = [0.0; 2];
                for (i, value) in values.iter().enumerate() {
                    let dof = 2 * mesh.tri_edges[cell][i / 2].0 + i % 2;
                    acceleration[0] += derivative[dof] * value[0];
                    acceleration[1] += derivative[dof] * value[1];
                }
                for (i, value) in values.iter().enumerate() {
                    let edge = mesh.tri_edges[cell][i / 2].0;
                    // Prescribed normal-moment rows are constraints, not
                    // momentum equations; all interior rows retain inertia.
                    if !mesh.is_boundary(edge) {
                        residual[2 * edge + i % 2] +=
                            weight * (value[0] * acceleration[0] + value[1] * acceleration[1]);
                    }
                }
            }
        }
        let defect = residual.iter().fold(0.0_f64, |v, r| v.hypot(*r));
        let load = rhs.iter().fold(0.0_f64, |v, r| v.hypot(*r));
        assert!(defect / load.max(1.0) < 2e-9);
    }
    assert_eq!((older, previous), saved_histories);
}

#[test]
fn bdf2_newton_startup_and_invalid_history_are_explicit() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-12,
        ..FluxParams::default()
    };
    let previous = vec![0.0; system.n_u];
    let force = |p: [f64; 2]| [2.0 * (p[1] - 0.5), -2.0 * (p[0] - 0.5)];
    let controls = FluxNewtonControls::default();
    let bdf1 = system
        .bdf1_rheology_newton_step(
            params,
            &thickening,
            &force,
            &|_| [0.0; 2],
            &previous,
            0.1,
            controls,
        )
        .unwrap();
    for startup in [
        system.bdf2_rheology_newton_step(
            params,
            &thickening,
            &force,
            &|_| [0.0; 2],
            &previous,
            None,
            0.1,
            controls,
        ),
        system.bdf2_rheology_newton_variable_step(
            params,
            &thickening,
            &force,
            &|_| [0.0; 2],
            &previous,
            None,
            0.1,
            controls,
        ),
    ] {
        let startup = startup.unwrap();
        assert_eq!(startup.solution.x, bdf1.solution.x);
        assert_eq!(startup.iterations, bdf1.iterations);
        assert_eq!(startup.backtracks, bdf1.backtracks);
        assert_eq!(startup.nonlinear_residual, bdf1.nonlinear_residual);
    }
    let mut nonfinite = previous.clone();
    nonfinite[0] = f64::NAN;
    for (older, previous_dt, timestep) in [
        (&previous[..system.n_u - 1], 0.1, 0.1),
        (&nonfinite[..], 0.1, 0.1),
        (&previous[..], 0.0, 0.1),
        (&previous[..], 0.1, 0.21),
    ] {
        let result = system.bdf2_rheology_newton_variable_step(
            params,
            &|_, _, _| panic!("invalid history queried constitutive data"),
            &|_| panic!("invalid history queried force"),
            &|_| panic!("invalid history queried boundary"),
            &previous,
            Some((older, previous_dt)),
            timestep,
            controls,
        );
        assert!(matches!(result, Err(FluxNewtonError::History(_))));
    }
    assert!(previous.iter().all(|value| *value == 0.0));
}

#[test]
fn steady_newton_solves_stiff_viscosity_without_a_transient_mass_term() {
    let mesh = mesh(4);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-11,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let initial = vec![0.0; system.n_u];
    let force = |x: [f64; 2]| [100.0 * (x[1] - 0.5), -100.0 * (x[0] - 0.5)];
    let controls = FluxNewtonControls {
        maximum_iterations: 30,
        maximum_backtracks: 16,
        residual_tolerance: 1e-9,
    };
    for convection in [false, true] {
        let result = system
            .steady_rheology_newton(
                params,
                &thickening,
                &force,
                &|_| [0.0; 2],
                &initial,
                convection,
                controls,
            )
            .unwrap();
        assert!(result.iterations < controls.maximum_iterations);
        assert!(result.backtracks > 0, "the stiff solve requires damping");
        assert!(result.nonlinear_residual <= controls.residual_tolerance);
        assert!(system.kinetic_energy(&result.solution.x) > 1e-5);
        let rates = system.shear_rates(&result.solution.x);
        // Assemble the original steady equations independently, with no mass
        // term, time history or large finite timestep standing in for steady.
        let (operator, rhs) = system.assemble_variable_viscosity(
            params,
            &|cell, point| thickening(cell, point, rates[cell]).nu,
            &force,
            &|_| [0.0; 2],
            convection.then_some(&result.solution.x[..system.n_u]),
        );
        let mut action = vec![0.0; system.n];
        operator.spmv(&result.solution.x, &mut action);
        let defect = action
            .iter()
            .zip(&rhs)
            .fold(0.0_f64, |v, (a, b)| v.hypot(a - b));
        let load = rhs.iter().fold(0.0_f64, |v, b| v.hypot(*b));
        assert!(defect / load < 1e-9);
        let (_, divergence) = system.velocity_error(&result.solution.x, &|_| [0.0; 2]);
        assert!(divergence < 1e-9);
    }
    assert!(initial.iter().all(|value| *value == 0.0));
}

#[test]
fn steady_shear_thinning_couette_matches_the_exact_boundary_driven_solution() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-12,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let initial = vec![0.0; system.n_u];
    let law = |_: usize, _: [f64; 2], rate: f64| {
        let excess = 0.7 * (1.0 + 9.0 * rate * rate).powf(-0.3);
        FluxViscosityTangent {
            nu: 0.1 + excess,
            derivative: -5.4 * rate * excess / (1.0 + 9.0 * rate * rate),
        }
    };
    let velocity = |x: [f64; 2]| [2.0 * x[1], 0.0];
    for convection in [false, true] {
        let result = system
            .steady_rheology_newton(
                params,
                &law,
                &|_| [0.0; 2],
                &velocity,
                &initial,
                convection,
                FluxNewtonControls {
                    residual_tolerance: 1e-10,
                    ..FluxNewtonControls::default()
                },
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&result.solution.x, &velocity);
        assert!(error < 1e-9 && divergence < 1e-9);
        for rate in system.shear_rates(&result.solution.x) {
            assert!((rate - 2.0).abs() < 1e-8);
        }
        for pressure in system.pressure(&result.solution.x) {
            assert!(pressure.abs() < 1e-8);
        }
    }
}

#[test]
fn damped_newton_solves_stiff_thickening_that_exhausts_picard() {
    let mesh = mesh(4);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-11,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let previous = vec![0.0; system.n_u];
    let original = previous.clone();
    let force = |x: [f64; 2]| [100.0 * (x[1] - 0.5), -100.0 * (x[0] - 0.5)];
    let controls = FluxNewtonControls {
        maximum_iterations: 30,
        maximum_backtracks: 16,
        residual_tolerance: 1e-9,
    };
    let result = system
        .bdf1_rheology_newton_step(
            params,
            &thickening,
            &force,
            &|_| [0.0; 2],
            &previous,
            0.3,
            controls,
        )
        .unwrap();
    assert_eq!(previous, original);
    assert!(result.iterations < 20);
    assert!(
        result.backtracks > 0,
        "fixture must actually require damping"
    );
    assert!(result.nonlinear_residual <= controls.residual_tolerance);
    assert!(system.kinetic_energy(&result.solution.x) > 1e-5);
    let rates = system.shear_rates(&result.solution.x);
    let (operator, rhs) = system.assemble_viscous_transient(
        params,
        &|cell, point| thickening(cell, point, rates[cell]).nu,
        &force,
        &|_| [0.0; 2],
        Some(&result.solution.x[..system.n_u]),
        1.0 / 0.3,
        Some(&previous),
    );
    let mut defect = vec![0.0; system.n];
    operator.spmv(&result.solution.x, &mut defect);
    let norm = defect
        .iter()
        .zip(&rhs)
        .fold(0.0_f64, |v, (a, b)| v.hypot(a - b));
    let load = rhs.iter().fold(0.0_f64, |v, b| v.hypot(*b));
    assert!(norm / load < 1e-9);
    let (_, divergence) = system.velocity_error(&result.solution.x, &|_| [0.0; 2]);
    assert!(divergence < 1e-9);
    // Independent existing Picard implementation, identical physical inputs
    // and solve budget. Its failure is not manufactured by an invalid law.
    let picard = system.bdf1_rheology_step(
        params,
        &|cell, point, rate| thickening(cell, point, rate).nu,
        &force,
        &|_| [0.0; 2],
        &previous,
        0.3,
        30,
        1e-9,
    );
    assert!(
        picard.is_err(),
        "this stiff case should expose Picard's convergence gap"
    );
    eprintln!(
        "stiff BDF1 Newton: iterations={}, rejected={}, residual={:e}, div={divergence:e}, Picard refused={}",
        result.iterations,
        result.backtracks,
        result.nonlinear_residual,
        picard.is_err()
    );
}

#[test]
fn newton_and_picard_agree_for_resolved_carreau_step() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        tol: 1e-12,
        ..FluxParams::default()
    };
    let previous = vec![0.0; system.n_u];
    let force = |x: [f64; 2]| [2.0 * (x[1] - 0.5), -2.0 * (x[0] - 0.5)];
    let law = |_: usize, _: [f64; 2], rate: f64| {
        let excess = 0.7 * (1.0 + 9.0 * rate * rate).powf(-0.3);
        FluxViscosityTangent {
            nu: 0.1 + excess,
            derivative: -5.4 * rate * excess / (1.0 + 9.0 * rate * rate),
        }
    };
    let controls = FluxNewtonControls {
        residual_tolerance: 1e-10,
        ..FluxNewtonControls::default()
    };
    let newton = system
        .bdf1_rheology_newton_step(
            params,
            &law,
            &force,
            &|_| [0.0; 2],
            &previous,
            0.1,
            controls,
        )
        .unwrap();
    let picard = system
        .bdf1_rheology_step(
            params,
            &|cell, point, rate| law(cell, point, rate).nu,
            &force,
            &|_| [0.0; 2],
            &previous,
            0.1,
            100,
            1e-10,
        )
        .unwrap();
    let difference = newton
        .solution
        .x
        .iter()
        .zip(&picard.solution.x)
        .fold(0.0_f64, |v, (a, b)| v.max((a - b).abs()));
    assert!(
        difference < 1e-9,
        "same nonlinear equations: difference={difference:e}"
    );
    eprintln!(
        "Carreau BDF1: Newton={} Picard={} maximum difference={difference:e}",
        newton.iterations, picard.iterations
    );
}

#[test]
fn newton_budgets_boundary_compatibility_and_zero_state_are_checked() {
    let mesh = mesh(3);
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let previous = vec![0.0; system.n_u];
    let controls = FluxNewtonControls::default();
    let zero = system
        .bdf1_rheology_newton_step(
            params,
            &thickening,
            &|_| [0.0; 2],
            &|_| [0.0; 2],
            &previous,
            0.1,
            controls,
        )
        .unwrap();
    assert_eq!(zero.iterations, 0);
    assert_eq!(zero.backtracks, 0);
    assert!(zero.solution.x.iter().all(|v| *v == 0.0));
    let force = |x: [f64; 2]| [100.0 * (x[1] - 0.5), -100.0 * (x[0] - 0.5)];
    let limited = FluxNewtonControls {
        maximum_iterations: 1,
        ..controls
    };
    assert!(matches!(
        system.bdf1_rheology_newton_step(
            params,
            &thickening,
            &force,
            &|_| [0.0; 2],
            &previous,
            0.3,
            limited
        ),
        Err(FluxNewtonError::NonConvergence { iterations: 1, .. })
    ));
    let undamped = FluxNewtonControls {
        maximum_backtracks: 0,
        ..controls
    };
    assert!(matches!(
        system.bdf1_rheology_newton_step(
            params,
            &thickening,
            &force,
            &|_| [0.0; 2],
            &previous,
            0.3,
            undamped
        ),
        Err(FluxNewtonError::LineSearch {
            attempted_steps: 1,
            ..
        })
    ));
    let incompatible = system.bdf1_rheology_newton_step(
        params,
        &thickening,
        &|_| [0.0; 2],
        &|x| [x[0], 0.0],
        &previous,
        0.1,
        controls,
    );
    assert!(matches!(incompatible, Err(FluxNewtonError::Boundary(_))));
    let invalid = FluxNewtonControls {
        maximum_backtracks: 53,
        ..controls
    };
    assert!(matches!(
        system.bdf1_rheology_newton_step(
            params,
            &thickening,
            &force,
            &|_| [0.0; 2],
            &previous,
            0.1,
            invalid
        ),
        Err(FluxNewtonError::Invalid(_))
    ));
    assert!(previous.iter().all(|v| *v == 0.0));
}
