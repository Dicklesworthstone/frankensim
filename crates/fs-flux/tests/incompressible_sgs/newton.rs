//! G1/G3: original nonlinear SGS residual, stiff solves and real time history.
use super::*;
use fs_flux::ns::linearization::FluxViscosityTangent;
use fs_flux::ns::newton::{FluxNewtonControls, FluxNewtonError};

#[path = "mixed_newton.rs"]
mod mixed;

#[path = "gradient_newton.rs"]
mod gradient;

#[test]
fn consistent_sgs_tangent_resolves_stiff_flow_within_a_bounded_newton_budget() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.7; mesh.tris.len()];
    let molecular = vec![1e-4; mesh.tris.len()];
    let model = PlanarSgs::new(&system, 0.5, &widths, &molecular).unwrap();
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-11,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls {
        maximum_iterations: 30,
        maximum_backtracks: 20,
        residual_tolerance: 1e-9,
    };
    let zero = |_: [f64; 2]| [0.0; 2];
    let force = |p: [f64; 2]| [100.0 * (p[1] - 0.5), -100.0 * (p[0] - 0.5)];
    let previous = vec![0.0; system.n_u];
    let original = previous.clone();
    let dt = 0.3;
    let result = model
        .bdf1_newton_step(params, &force, &zero, &previous, dt, controls)
        .unwrap();
    assert!(result.step.iterations < controls.maximum_iterations);
    assert!(result.step.nonlinear_residual <= controls.residual_tolerance);
    assert!(result.sgs.integrated_dissipation_m4_s3 > 0.0);

    // Independently reassemble the ORIGINAL nonlinear momentum equations.
    // The control with SGS held at its zero-strain initial value must fail.
    let x = &result.step.solution.x;
    let rates = system.shear_rates(x);
    for active in [true, false] {
        let (operator, rhs) = system.assemble_viscous_transient(
            params,
            &|cell, _| {
                molecular[cell]
                    + if active {
                        (0.5 * widths[cell]).powi(2) * rates[cell]
                    } else {
                        0.0
                    }
            },
            &force,
            &zero,
            Some(&x[..system.n_u]),
            dt.recip(),
            Some(&previous),
        );
        let mut action = vec![0.0; system.n];
        operator.spmv(x, &mut action);
        let defect = action
            .iter()
            .zip(&rhs)
            .fold(0.0_f64, |v, (a, b)| v.hypot(a - b));
        let load = rhs.iter().fold(0.0_f64, |v, b| v.hypot(*b));
        if active {
            assert!(
                defect / load < 1e-9,
                "nonlinear SGS defect {}",
                defect / load
            );
        } else {
            assert!(defect / load > 1e-3, "frozen model control was ineffective");
        }
    }
    let boundaries = velocity_boundaries(&mesh, &zero);
    let picard = model.bdf1_step(
        params,
        &force,
        &boundaries,
        &previous,
        dt,
        NavierStokesControls {
            maximum_iterations: controls.maximum_iterations,
            relative_tolerance: controls.residual_tolerance,
            relaxation: 1.0,
        },
    );
    eprintln!(
        "SGS Newton: {} iterations, {} backtracks, residual {:e}; Picard: {:?}",
        result.step.iterations,
        result.step.backtracks,
        result.step.nonlinear_residual,
        picard
            .as_ref()
            .map(|step| (step.step.iterations, step.step.nonlinear_residual))
    );
    assert!(
        picard.is_err(),
        "the bounded fixture must expose slow SGS fixed-point convergence"
    );
    assert_eq!(previous, original);
}

#[test]
fn smagorinsky_newton_respects_unequal_physical_intervals_and_endpoint_transfer() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.5; mesh.tris.len()];
    let molecular = vec![0.02; mesh.tris.len()];
    let model = PlanarSgs::new(&system, 0.3, &widths, &molecular).unwrap();
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-12,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls {
        residual_tolerance: 1e-10,
        ..FluxNewtonControls::default()
    };
    let amplitude = |time: f64| 1.0 + time + time * time;
    let previous_dt = 0.3;
    let time = 0.4;
    let older = moments(&mesh, &|p| [amplitude(time - previous_dt) * p[1], 0.0]);
    let previous = moments(&mesh, &|p| [amplitude(time) * p[1], 0.0]);
    let snapshots = (older.clone(), previous.clone());
    for dt in [0.15, 0.3, 0.6] {
        let endpoint = time + dt;
        let shear = amplitude(endpoint);
        let velocity = |p: [f64; 2]| [shear * p[1], 0.0];
        // BDF2 differentiates this quadratic time dependence exactly, including
        // unequal intervals. Uniform stress and convection have zero divergence.
        let force = |p: [f64; 2]| [(1.0 + 2.0 * endpoint) * p[1], 0.0];
        let result = model
            .bdf2_newton_variable_step(
                params,
                &force,
                &velocity,
                &previous,
                Some((&older, previous_dt)),
                dt,
                controls,
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&result.step.solution.x, &velocity);
        assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
        let eddy = (0.3_f64 * 0.5).powi(2) * shear;
        for cell in &result.sgs.cells {
            close(cell.strain_rate_per_s, shear, 4e-8);
            close(cell.eddy_viscosity_m2_s, eddy, 3e-9);
            close(cell.dissipation_m2_s3, eddy * shear * shear, 3e-8);
        }
        close(
            result.sgs.integrated_dissipation_m4_s3,
            eddy * shear * shear,
            3e-8,
        );
    }
    assert_eq!((older, previous), snapshots);
}

#[test]
fn zero_newton_sgs_is_the_native_molecular_steady_bdf1_and_bdf2_owner() {
    for closure in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        zero_newton_sgs(closure);
    }
}

fn zero_newton_sgs(closure: SgsModel) {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let molecular: Vec<_> = (0..mesh.tris.len())
        .map(|cell| 0.2 + 0.01 * cell as f64)
        .collect();
    let model =
        PlanarSgs::with_model(&system, closure, &vec![0.5; mesh.tris.len()], &molecular).unwrap();
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-12,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls::default();
    let zero = |_: [f64; 2]| [0.0; 2];
    let force = |p: [f64; 2]| [p[1] - 0.5, 0.5 - p[0]];
    let initial = vec![0.0; system.n_u];
    let old = moments(&mesh, &|p| [0.1 * p[1], 0.05 * p[0]]);
    let law = |cell: usize, _: [f64; 2], _: f64| FluxViscosityTangent {
        nu: molecular[cell],
        derivative: 0.0,
    };
    let reference = system
        .steady_rheology_newton(params, &law, &force, &zero, &initial, true, controls)
        .unwrap();
    let result = model
        .solve_steady_newton(params, &force, &zero, &initial, controls)
        .unwrap();
    assert_eq!(result.step.solution.x, reference.solution.x);
    assert_eq!(result.step.nonlinear_residual, reference.nonlinear_residual);
    for history in [None, Some((old.as_slice(), 0.04))] {
        let reference = system
            .bdf2_rheology_newton_variable_step(
                params, &law, &force, &zero, &initial, history, 0.02, controls,
            )
            .unwrap();
        let result = model
            .bdf2_newton_variable_step(params, &force, &zero, &initial, history, 0.02, controls)
            .unwrap();
        assert_eq!(result.step.solution.x, reference.solution.x);
        assert_eq!(result.step.iterations, reference.iterations);
        assert_eq!(result.step.backtracks, reference.backtracks);
        assert_eq!(result.step.nonlinear_residual, reference.nonlinear_residual);
        assert!(
            result
                .sgs
                .cells
                .iter()
                .all(|cell| cell.eddy_viscosity_m2_s == 0.0 && cell.dissipation_m2_s3 == 0.0)
        );
    }
    let startup = model
        .bdf2_newton_variable_step(params, &force, &zero, &initial, None, 0.02, controls)
        .unwrap();
    let first = model
        .bdf1_newton_step(params, &force, &zero, &initial, 0.02, controls)
        .unwrap();
    assert_eq!(startup.step.solution.x, first.step.solution.x);
}

#[test]
fn sgs_newton_admits_zero_gradient_and_refuses_bad_history_or_exhausted_work_atomically() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let molecular = vec![0.01; mesh.tris.len()];
    let widths = vec![0.5; mesh.tris.len()];
    let params = FluxParams::default();
    let controls = FluxNewtonControls::default();
    let previous = vec![0.0; system.n_u];
    let snapshot = previous.clone();
    let zero = |_: [f64; 2]| [0.0; 2];
    let never = |_: [f64; 2]| panic!("refusal must precede external data queries");
    for closure in [SgsModel::Wale { cw: 0.5 }, SgsModel::Vreman { cv: 0.07 }] {
        let model = PlanarSgs::with_model(&system, closure, &widths, &molecular).unwrap();
        let solved = model
            .solve_steady_newton(params, &zero, &zero, &previous, controls)
            .unwrap();
        assert_eq!(solved.step.iterations, 0);
        assert_eq!(solved.step.solution.x, vec![0.0; system.n]);
        assert!(
            solved
                .sgs
                .cells
                .iter()
                .all(|cell| cell.eddy_viscosity_m2_s == 0.0)
        );
    }
    // The law value can be finite at tiny strain even when its Newton tangent
    // is not representable. Refuse the tangent instead of silently freezing it.
    for width in [1e160, 1e-200] {
        let model =
            PlanarSgs::new(&system, 1.0, &vec![width; mesh.tris.len()], &molecular).unwrap();
        assert!(matches!(
            model.solve_steady_newton(params, &never, &never, &previous, controls),
            Err(PlanarSgsError::UnrepresentableTransport { .. })
        ));
    }
    let model = PlanarSgs::new(&system, 0.3, &widths, &molecular).unwrap();
    assert!(matches!(
        model.bdf2_newton_variable_step(
            params,
            &never,
            &never,
            &previous,
            Some((&previous, 0.1)),
            0.21,
            controls
        ),
        Err(PlanarSgsError::Newton(FluxNewtonError::History(_)))
    ));
    let force = |p: [f64; 2]| [10.0 * (p[1] - 0.5), 10.0 * (0.5 - p[0])];
    assert!(matches!(
        model.bdf1_newton_step(
            params,
            &force,
            &zero,
            &previous,
            0.1,
            FluxNewtonControls {
                maximum_iterations: 1,
                ..controls
            }
        ),
        Err(PlanarSgsError::Newton(
            FluxNewtonError::NonConvergence { .. }
        ))
    ));
    assert_eq!(previous, snapshot);
}
