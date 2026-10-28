//! G1/G3: stiff constitutive solves, direct residual proof and budget refusal.
#![cfg(feature = "continuum")]

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
