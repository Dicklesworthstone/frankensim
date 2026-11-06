//! G0/G1/G3: true steady SGS flow, original mixed boundaries and typed refusal.
use super::*;
use fs_flux::ns::boundary::{StokesBoundaryError, StokesPressureReference};
use std::cell::Cell;

#[test]
fn steady_sgs_closes_full_gradient_traction_and_convection_without_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-12,
        ..FluxParams::default()
    };
    let controls = NavierStokesControls {
        maximum_iterations: 100,
        relative_tolerance: 1e-10,
        relaxation: 0.8,
    };
    let widths = vec![0.7; mesh.tris.len()];
    let molecular = vec![0.5; mesh.tris.len()];
    // Equal symmetric strain but different rotation: a scalar-rate substitute
    // cannot reproduce both WALE/Vreman solutions and their physical tractions.
    let gradients = [[[0.3, 0.2], [0.0, -0.3]], [[0.3, 0.6], [-0.4, -0.3]]];
    let initial_guesses = [
        vec![0.0; system.n_u],
        affine_moments(&mesh, [[0.1, 0.15], [0.0, -0.1]], [0.2, 0.03]),
    ];
    let snapshots = initial_guesses.clone();
    for closure in [
        SgsModel::Smagorinsky { cs: 0.5 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let model = PlanarSgs::with_model(&system, closure, &widths, &molecular).unwrap();
        let mut accepted_eddy = Vec::new();
        for (case, g) in gradients.into_iter().enumerate() {
            let target = |p: [f64; 2]| {
                [
                    0.4 + g[0][0] * p[0] + g[0][1] * p[1],
                    g[1][0] * p[0] + g[1][1] * p[1],
                ]
            };
            let eddy = trace_free_planar_reference(closure, 0.7, g);
            let nu = 0.5 + eddy;
            assert!(eddy > 0.0);
            let policy = if case == 0 {
                OpenBoundaryPolicy::NaturalTraction
            } else {
                OpenBoundaryPolicy::DissipativeBackflow
            };
            // Independent affine derivatives: div(2 nu D)=0 and div(u)=0.
            let force = |p: [f64; 2]| {
                let u = target(p);
                [
                    g[0][0] * u[0] + g[0][1] * u[1] + 2.0,
                    g[1][0] * u[0] + g[1][1] * u[1] - 3.0,
                ]
            };
            let traction = |p: [f64; 2], n: [f64; 2]| {
                let pressure = 5.0 + 2.0 * p[0] - 3.0 * p[1];
                let shear = nu * (g[0][1] + g[1][0]);
                let mut t = [
                    (2.0 * nu * g[0][0] - pressure) * n[0] + shear * n[1],
                    shear * n[0] + (2.0 * nu * g[1][1] - pressure) * n[1],
                ];
                if policy == OpenBoundaryPolicy::DissipativeBackflow {
                    let u = target(p);
                    let inward = (u[0] * n[0] + u[1] * n[1]).min(0.0);
                    for axis in 0..2 {
                        t[axis] -= inward * u[axis];
                    }
                }
                t
            };
            // For the first field, y=0 is an impermeable Navier-slip wall:
            // tau_xn=-.2*nu=-beta*(u_x-wall_x), with beta=.8.
            let wall = |p: [f64; 2]| [target(p)[0] - 0.25 * nu, 0.0];
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, face)| {
                    if face.normal[0].abs() > 0.5 {
                        NavierStokesBoundary::Open {
                            edge,
                            traction: &traction,
                            backflow: policy,
                        }
                    } else if case == 0 && face.normal[1] < 0.0 {
                        NavierStokesBoundary::NavierSlip {
                            edge,
                            value: &wall,
                            tangential_resistance_m_s: 0.8,
                        }
                    } else {
                        NavierStokesBoundary::Velocity {
                            edge,
                            value: &target,
                        }
                    }
                })
                .collect();
            let mut solutions = Vec::new();
            for initial in &initial_guesses {
                let result = model
                    .solve_steady(params, &force, &boundaries, initial, controls)
                    .unwrap();
                let x = &result.step.solution.field.x;
                let (error, divergence) = system.velocity_error(x, &target);
                assert!(error < 5e-9 && divergence < 5e-9, "{error}, {divergence}");
                assert!(result.step.relative_update <= controls.relative_tolerance);
                assert!(result.step.nonlinear_residual <= controls.relative_tolerance);
                assert!(result.step.iterations > 1);
                assert_eq!(
                    result.step.solution.pressure_reference,
                    StokesPressureReference::TractionDatum
                );
                assert!(
                    result
                        .step
                        .open_boundary_policies
                        .iter()
                        .all(|(_, retained)| *retained == policy)
                );
                for (actual, p) in result.step.solution.pressure.iter().zip(&mesh.centroids) {
                    close(*actual, 5.0 + 2.0 * p[0] - 3.0 * p[1], 6e-8);
                }
                for cell in &result.sgs.cells {
                    close(cell.strain_rate_per_s, 0.4_f64.sqrt(), 3e-8);
                    close(cell.eddy_viscosity_m2_s, eddy, 5e-9);
                    close(cell.effective_viscosity_m2_s, nu, 5e-9);
                }
                close(result.sgs.integrated_dissipation_m4_s3, 0.4 * eddy, 5e-9);

                // Independent original-operator residual. This diagnostic BDF
                // assembly uses IDENTICAL solution/history velocities, cancelling
                // its mass terms; the steady solve above receives no timestep.
                // Freezing SGS at its zero-strain initial value must fail.
                for include_sgs in [true, false] {
                    let assembled = system
                        .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                            params,
                            &|cell, _| {
                                if include_sgs {
                                    result.sgs.cells[cell].effective_viscosity_m2_s
                                } else {
                                    molecular[cell]
                                }
                            },
                            &force,
                            &boundaries,
                            &x[..system.n_u],
                            &x[..system.n_u],
                            0.125,
                        )
                        .unwrap();
                    let mut action = vec![0.0; system.n];
                    assembled.matrix().spmv(x, &mut action);
                    let defect = action
                        .iter()
                        .zip(assembled.rhs())
                        .fold(0.0_f64, |norm, (a, b)| norm.hypot(a - b));
                    if include_sgs {
                        assert!(defect < 5e-8, "accepted steady defect {defect}");
                    } else {
                        assert!(
                            defect > 1e-5,
                            "frozen molecular stress was not detected: {defect}"
                        );
                    }
                }
                accepted_eddy.push(result.sgs.cells[0].eddy_viscosity_m2_s);
                solutions.push(x.clone());
            }
            // Initial guesses affect only nonlinear work, never a mass load.
            for (a, b) in solutions[0].iter().zip(&solutions[1]) {
                close(*a, *b, 8e-8);
            }
        }
        if !matches!(closure, SgsModel::Smagorinsky { .. }) {
            assert!((accepted_eddy[0] - accepted_eddy[2]).abs() > 1e-5);
        }
    }
    assert_eq!(initial_guesses, snapshots);
}

#[test]
fn steady_zero_closures_match_the_original_variable_molecular_flow_bitwise() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.5; mesh.tris.len()];
    let molecular: Vec<_> = (0..mesh.tris.len())
        .map(|cell| 0.3 + 0.01 * cell as f64)
        .collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    let force = |p: [f64; 2]| [p[1] - 0.5, 0.5 - p[0]];
    let initial = moments(&mesh, &|p| [0.1 * p[1], 0.05 * p[0]]);
    let snapshot = initial.clone();
    let params = FluxParams {
        nu: f64::NAN,
        ..FluxParams::default()
    };
    let controls = NavierStokesControls::default();
    let reference = system
        .solve_navier_stokes_mixed_rheology(
            params,
            &|cell, _, _| molecular[cell],
            &force,
            &boundaries,
            &initial,
            controls,
        )
        .unwrap();
    for closure in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let model = PlanarSgs::with_model(&system, closure, &widths, &molecular).unwrap();
        let result = model
            .solve_steady(params, &force, &boundaries, &initial, controls)
            .unwrap();
        assert_eq!(result.step.solution.field.x, reference.solution.field.x);
        assert_eq!(result.step.solution.pressure, reference.solution.pressure);
        assert_eq!(result.step.iterations, reference.iterations);
        assert_eq!(result.step.relative_update, reference.relative_update);
        assert_eq!(result.step.nonlinear_residual, reference.nonlinear_residual);
        assert!(
            result
                .sgs
                .cells
                .iter()
                .all(|cell| cell.eddy_viscosity_m2_s == 0.0 && cell.dissipation_m2_s3 == 0.0)
        );
    }
    assert_eq!(initial, snapshot);
}

#[test]
fn steady_sgs_refuses_bad_guesses_unanchored_flow_overflow_and_exhausted_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let molecular = vec![0.1; mesh.tris.len()];
    let initial = affine_moments(&mesh, [[0.3, 0.2], [0.0, -0.3]], [0.1; 2]);
    let snapshot = initial.clone();
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    let calls = Cell::new(0);
    let observed_force = |_: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let observed_traction = |_: [f64; 2], _: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let open: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| NavierStokesBoundary::Open {
            edge,
            traction: &observed_traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    for closure in [
        SgsModel::Smagorinsky { cs: 0.5 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let model =
            PlanarSgs::with_model(&system, closure, &vec![0.5; mesh.tris.len()], &molecular)
                .unwrap();
        let params = FluxParams::default();
        let controls = NavierStokesControls::default();
        for invalid in [vec![], vec![f64::NAN; system.n_u]] {
            assert!(matches!(
                model.solve_steady(params, &observed_force, &boundaries, &invalid, controls),
                Err(PlanarSgsError::Flow(StokesBoundaryError::InvalidInput(_)))
            ));
        }
        assert!(matches!(
            model.solve_steady(params, &observed_force, &open, &initial, controls),
            Err(PlanarSgsError::Flow(
                StokesBoundaryError::NoVelocityBoundary
            ))
        ));
        assert_eq!(calls.get(), 0);
        let oversized =
            PlanarSgs::with_model(&system, closure, &vec![1e160; mesh.tris.len()], &molecular)
                .unwrap();
        for _ in 0..2 {
            assert!(matches!(
                oversized.solve_steady(params, &zero, &boundaries, &initial, controls),
                Err(PlanarSgsError::UnrepresentableTransport { .. })
            ));
        }
        assert!(matches!(
            model.solve_steady(
                params,
                &|p| [p[1] - 0.5, 0.5 - p[0]],
                &boundaries,
                &initial,
                NavierStokesControls {
                    maximum_iterations: 1,
                    ..controls
                },
            ),
            Err(PlanarSgsError::Flow(StokesBoundaryError::Rheology(
                fs_flux::FluxRheologyError::NonlinearConvergence { .. }
            )))
        ));
        assert_eq!(initial, snapshot);
    }
}
