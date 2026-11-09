//! Exact mixed convection, physical traction/work and implicit SGS endpoints.
use super::*;
use fs_flux::ns::boundary::{StokesBoundaryError, StokesPressureReference};

fn norm(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |value, next| value.hypot(*next))
}

#[test]
fn mixed_newton_jacobian_matches_original_fluxes_under_reversal() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: f64::NAN,
        sigma: 15.0,
        ..FluxParams::default()
    };
    // Both vertical edges reverse normal flow inside an original face.
    // The horizontal slip traces also have nonzero normal transport.
    let target = |p: [f64; 2]| [0.61 - 1.4 * p[1], 0.2 * p[0] - 0.11];
    let reservoir = |p: [f64; 2]| [-0.37 + 0.1 * p[1], 0.29];
    let traction = |_: [f64; 2], n: [f64; 2]| [-1.3 * n[0], -1.3 * n[1]];
    let force = |_: [f64; 2]| [0.13, -0.07];
    let previous = moments(&mesh, &|_| [0.12, -0.03]);
    let mut state = moments(&mesh, &target);
    state.resize(system.n, 0.23);
    let direction: Vec<_> = (0..system.n)
        .map(|i| ((i * 13 % 29) as f64 - 14.0) / 31.0)
        .collect();
    let dt = 0.4;
    for slip in [false, true] {
        for policy in 0..3 {
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, face)| {
                    if face.normal[0].abs() > 0.5 {
                        if policy == 2 {
                            NavierStokesBoundary::OpenReservoir {
                                edge,
                                traction: &traction,
                                velocity: &reservoir,
                            }
                        } else {
                            NavierStokesBoundary::Open {
                                edge,
                                traction: &traction,
                                backflow: if policy == 0 {
                                    OpenBoundaryPolicy::NaturalTraction
                                } else {
                                    OpenBoundaryPolicy::DissipativeBackflow
                                },
                            }
                        }
                    } else if slip {
                        NavierStokesBoundary::NavierSlip {
                            edge,
                            value: &target,
                            tangential_resistance_m_s: 0.7,
                        }
                    } else {
                        NavierStokesBoundary::Velocity {
                            edge,
                            value: &target,
                        }
                    }
                })
                .collect();
            for tangent in [0.0, 0.07] {
                let law = |_: usize, _: [f64; 2], rate: f64| FluxViscosityTangent {
                    nu: 0.03 + tangent * rate,
                    derivative: if rate == 0.0 { 0.0 } else { tangent },
                };
                let linearization = system
                    .linearize_navier_stokes_mixed_rheology(
                        params,
                        &law,
                        &force,
                        &boundaries,
                        &state,
                        1.0 / dt,
                        Some(&previous),
                    )
                    .unwrap();
                let mut exact = vec![0.0; system.n];
                linearization.jacobian.spmv(&direction, &mut exact);
                // Independent centered differences of the ORIGINAL public
                // primal assembly, not differences of Newton's residual helper.
                let epsilon = 1e-6;
                let residual = |sign: f64| {
                    let x: Vec<_> = state
                        .iter()
                        .zip(&direction)
                        .map(|(value, step)| value + sign * epsilon * step)
                        .collect();
                    let rates = system.shear_rates(&x);
                    let assembled = system
                        .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                            params,
                            &|cell, _| 0.03 + tangent * rates[cell],
                            &force,
                            &boundaries,
                            &x[..system.n_u],
                            &previous,
                            dt,
                        )
                        .unwrap();
                    let mut result = vec![0.0; system.n];
                    assembled.matrix().spmv(&x, &mut result);
                    for (value, load) in result.iter_mut().zip(assembled.rhs()) {
                        *value -= load;
                    }
                    result
                };
                let plus = residual(1.0);
                let minus = residual(-1.0);
                let difference: Vec<_> = exact
                    .iter()
                    .zip(plus.iter().zip(&minus))
                    .map(|(a, (p, m))| a - (p - m) / (2.0 * epsilon))
                    .collect();
                let relative = norm(&difference) / norm(&exact);
                assert!(
                    relative < 2e-7,
                    "slip={slip}, policy={policy}, tangent={tangent}: {relative:e}"
                );
            }
        }
    }
}

#[test]
fn mixed_sgs_newton_preserves_analytic_pressure_traction_and_actual_time_history() {
    use std::cell::Cell;
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let model = PlanarSgs::new(
        &system,
        0.5,
        &vec![0.7; mesh.tris.len()],
        &vec![0.3; mesh.tris.len()],
    )
    .unwrap();
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
    let time = 0.4;
    let previous_dt = 0.3;
    for policy in 0..3 {
        for (scheme, dt) in [(0, 0.0), (1, 0.2), (2, 0.15), (2, 0.3), (2, 0.6)] {
            let amplitude = |t: f64| {
                if scheme == 0 {
                    1.0
                } else if scheme == 1 {
                    1.0 + t
                } else {
                    1.0 + t + t * t
                }
            };
            let a = amplitude(time + dt);
            let rate = if scheme == 0 {
                0.0
            } else if scheme == 1 {
                1.0
            } else {
                1.0 + 2.0 * (time + dt)
            };
            let base = |p: [f64; 2]| [0.4 + 0.3 * p[0] + 0.2 * p[1], -0.3 * p[1]];
            let target = |p: [f64; 2]| base(p).map(|u| a * u);
            let force = |p: [f64; 2]| {
                let u = base(p);
                [
                    rate * u[0] + a * a * (0.3 * u[0] + 0.2 * u[1]) + 2.0,
                    rate * u[1] - a * a * 0.3 * u[1] - 3.0,
                ]
            };
            let gamma = a * 0.4_f64.sqrt();
            let eddy = (0.5_f64 * 0.7).powi(2) * gamma;
            let nu = 0.3 + eddy;
            let physical_stress = |p: [f64; 2], n: [f64; 2]| {
                let pressure = 5.0 + 2.0 * p[0] - 3.0 * p[1];
                [
                    (0.6 * a * nu - pressure) * n[0] + 0.2 * a * nu * n[1],
                    0.2 * a * nu * n[0] + (-0.6 * a * nu - pressure) * n[1],
                ]
            };
            let traction = |p: [f64; 2], n: [f64; 2]| {
                let mut t = physical_stress(p, n);
                if policy == 1 {
                    let u = target(p);
                    let incoming = (u[0] * n[0] + u[1] * n[1]).min(0.0);
                    for axis in 0..2 {
                        t[axis] -= incoming * u[axis];
                    }
                }
                t
            };
            let wall = |p: [f64; 2]| [target(p)[0] - 0.25 * a * nu, 0.0];
            let reservoir_queries = Cell::new(0);
            let reservoir = |p: [f64; 2]| {
                reservoir_queries.set(reservoir_queries.get() + 1);
                target(p)
            };
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, face)| {
                    if face.normal[0].abs() > 0.5 {
                        if policy == 2 {
                            NavierStokesBoundary::OpenReservoir {
                                edge,
                                traction: &traction,
                                velocity: &reservoir,
                            }
                        } else {
                            NavierStokesBoundary::Open {
                                edge,
                                traction: &traction,
                                backflow: if policy == 0 {
                                    OpenBoundaryPolicy::NaturalTraction
                                } else {
                                    OpenBoundaryPolicy::DissipativeBackflow
                                },
                            }
                        }
                    } else if face.normal[1] < 0.0 {
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
            let previous = moments(&mesh, &|p| base(p).map(|u| amplitude(time) * u));
            let older = moments(&mesh, &|p| {
                base(p).map(|u| amplitude(time - previous_dt) * u)
            });
            let original = (previous.clone(), older.clone());
            let result = match scheme {
                0 => model.solve_steady_mixed_newton(
                    params,
                    &force,
                    &boundaries,
                    &vec![0.0; system.n_u],
                    controls,
                ),
                1 => model.bdf1_mixed_newton_step(
                    params,
                    &force,
                    &boundaries,
                    &previous,
                    dt,
                    controls,
                ),
                _ => model.bdf2_mixed_newton_variable_step(
                    params,
                    &force,
                    &boundaries,
                    &previous,
                    Some((&older, previous_dt)),
                    dt,
                    controls,
                ),
            }
            .unwrap();
            assert_eq!((previous, older), original);
            let (error, divergence) = system.velocity_error(&result.step.solution.field.x, &target);
            assert!(
                error < 2e-8 && divergence < 2e-8,
                "scheme={scheme}, policy={policy}, dt={dt}: {error:e}, {divergence:e}"
            );
            assert_eq!(
                result.step.solution.pressure_reference,
                StokesPressureReference::TractionDatum
            );
            for (actual, p) in result.step.solution.pressure.iter().zip(&mesh.centroids) {
                close(*actual, 5.0 + 2.0 * p[0] - 3.0 * p[1], 8e-8);
            }
            if policy == 2 {
                assert_eq!(
                    reservoir_queries.get(),
                    2 * boundaries
                        .iter()
                        .filter(|b| matches!(b, NavierStokesBoundary::OpenReservoir { .. }))
                        .count()
                );
            }
            for diagnostic in &result.sgs.cells {
                close(diagnostic.strain_rate_per_s, gamma, 5e-8);
                close(diagnostic.eddy_viscosity_m2_s, eddy, 1e-8);
            }
            close(
                result.sgs.integrated_dissipation_m4_s3,
                eddy * gamma * gamma,
                8e-8,
            );
            let mut slip_loss = 0.0;
            for face in result.step.solution.reactions.edges() {
                let geometry = &mesh.edges[face.edge];
                let load = face.load([0.0; 2], None).unwrap();
                slip_loss += load.slip_dissipation;
                if geometry.normal[0].abs() > 0.5 {
                    // Independent rigid virtual work of the actual Cauchy
                    // traction; added backflow stress must retain its sign.
                    let expected: f64 =
                        edge_gauss_pub(mesh.verts[geometry.verts.0], mesh.verts[geometry.verts.1])
                            .iter()
                            .map(|(p, w)| {
                                let t = physical_stress(*p, geometry.normal);
                                -w * (t[0] * (0.2 - 0.17 * p[1]) + t[1] * (-0.1 + 0.17 * p[0]))
                            })
                            .sum();
                    close(
                        load.stress.rigid_power([0.2, -0.1], 0.17).unwrap(),
                        expected,
                        1e-7,
                    );
                    if policy == 1 && geometry.normal[0] < 0.0 {
                        assert!(load.backflow.force[0] > 0.0);
                    }
                }
            }
            close(slip_loss, 0.8 * (0.25 * a * nu).powi(2), 1e-8);
        }
    }
}

#[test]
fn mixed_newton_zero_sgs_matches_native_molecular_work_and_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let molecular: Vec<_> = (0..mesh.tris.len())
        .map(|i| 0.2 + 0.01 * i as f64)
        .collect();
    let model = PlanarSgs::new(&system, 0.0, &vec![0.5; mesh.tris.len()], &molecular).unwrap();
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-12,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls::default();
    let zero = |_: [f64; 2]| [0.0; 2];
    let force = |p: [f64; 2]| [0.3 - 0.2 * p[1], 0.1 * p[0]];
    let traction = |p: [f64; 2], n: [f64; 2]| n.map(|x| -(2.0 - p[0]) * x);
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
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                NavierStokesBoundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let previous = vec![0.0; system.n_u];
    let older = moments(&mesh, &|p| [0.02 * p[1], 0.01 * p[0]]);
    let law = |cell: usize, _: [f64; 2], _: f64| FluxViscosityTangent {
        nu: molecular[cell],
        derivative: 0.0,
    };
    for scheme in 0..3 {
        let native = match scheme {
            0 => system.solve_navier_stokes_mixed_rheology_newton(
                params,
                &law,
                &force,
                &boundaries,
                &previous,
                controls,
            ),
            1 => system.bdf1_navier_stokes_mixed_rheology_newton_step(
                params,
                &law,
                &force,
                &boundaries,
                &previous,
                0.02,
                controls,
            ),
            _ => system.bdf2_navier_stokes_mixed_rheology_newton_variable_step(
                params,
                &law,
                &force,
                &boundaries,
                &previous,
                Some((&older, 0.04)),
                0.02,
                controls,
            ),
        }
        .unwrap();
        let result = match scheme {
            0 => model.solve_steady_mixed_newton(params, &force, &boundaries, &previous, controls),
            1 => {
                model.bdf1_mixed_newton_step(params, &force, &boundaries, &previous, 0.02, controls)
            }
            _ => model.bdf2_mixed_newton_variable_step(
                params,
                &force,
                &boundaries,
                &previous,
                Some((&older, 0.04)),
                0.02,
                controls,
            ),
        }
        .unwrap();
        assert_eq!(result.step.solution.field.x, native.solution.field.x);
        assert_eq!(result.step.solution.pressure, native.solution.pressure);
        assert_eq!(result.step.iterations, native.iterations);
        assert_eq!(result.step.backtracks, native.backtracks);
        assert_eq!(
            result.step.open_boundary_policies,
            native.open_boundary_policies
        );
        for (a, b) in result
            .step
            .solution
            .reactions
            .edges()
            .iter()
            .zip(native.solution.reactions.edges())
        {
            assert_eq!(
                a.load([0.0; 2], None).unwrap(),
                b.load([0.0; 2], None).unwrap()
            );
        }
        assert!(
            result
                .sgs
                .cells
                .iter()
                .all(|cell| cell.eddy_viscosity_m2_s == 0.0 && cell.dissipation_m2_s3 == 0.0)
        );
    }
}

#[test]
fn pressure_driven_sgs_newton_converges_with_original_open_stress() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 4, 3));
    let system = FluxSystem::new(&mesh);
    let molecular = vec![0.005; mesh.tris.len()];
    let model = PlanarSgs::new(&system, 0.5, &vec![0.7; mesh.tris.len()], &molecular).unwrap();
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
    let traction = |p: [f64; 2], n: [f64; 2]| n.map(|v| if p[0] < 1.0 { -3.0 * v } else { 0.0 });
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
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            } else {
                NavierStokesBoundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let initial = vec![0.0; system.n_u];
    let result = model
        .solve_steady_mixed_newton(params, &zero, &boundaries, &initial, controls)
        .unwrap();
    assert!(result.step.iterations < controls.maximum_iterations);
    assert!(result.sgs.integrated_dissipation_m4_s3 > 0.0);
    let outlet: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(i, e)| mesh.is_boundary(*i) && e.normal[0] > 0.5)
        .map(|(i, _)| i)
        .collect();
    assert!(
        result
            .step
            .solution
            .reactions
            .volume_flux_moments(&outlet, [0.0; 2])
            .unwrap()
            .0
            > 0.0
    );
    let x = &result.step.solution.field.x;
    for active in [true, false] {
        // Equal diagnostic history cancels mass without adding any physical
        // timestep to the steady solve. Reuse the original public assembly.
        let assembled = system
            .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                params,
                &|cell, _| {
                    if active {
                        result.sgs.cells[cell].effective_viscosity_m2_s
                    } else {
                        molecular[cell]
                    }
                },
                &zero,
                &boundaries,
                &x[..system.n_u],
                &x[..system.n_u],
                0.125,
            )
            .unwrap();
        let mut action = vec![0.0; system.n];
        assembled.matrix().spmv(x, &mut action);
        for (value, load) in action.iter_mut().zip(assembled.rhs()) {
            *value -= load;
        }
        let relative = norm(&action) / norm(assembled.rhs());
        if active {
            assert!(relative < 1e-9, "{relative:e}");
        } else {
            assert!(relative > 1e-3, "frozen SGS failed to expose stress error");
        }
    }
    let picard = model.solve_steady(
        params,
        &zero,
        &boundaries,
        &initial,
        NavierStokesControls {
            maximum_iterations: 30,
            relative_tolerance: 1e-9,
            relaxation: 1.0,
        },
    );
    eprintln!(
        "Open SGS Newton: {} iterations, {} backtracks, residual {:e}; Picard: {:?}",
        result.step.iterations,
        result.step.backtracks,
        result.step.nonlinear_residual,
        picard
            .as_ref()
            .map(|r| (r.step.iterations, r.step.nonlinear_residual))
    );
    match picard {
        Ok(reference) => {
            for (a, b) in x.iter().zip(reference.step.solution.field.x) {
                close(*a, b, 1e-7);
            }
        }
        Err(PlanarSgsError::Flow(StokesBoundaryError::Rheology(
            fs_flux::FluxRheologyError::NonlinearConvergence { iterations: 30, .. },
        ))) => {}
        Err(other) => panic!("unrelated Picard refusal: {other:?}"),
    }
    assert_eq!(initial, vec![0.0; system.n_u]);
}

#[test]
fn mixed_newton_mass_anchor_and_refusals_preserve_original_histories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.7; mesh.tris.len()];
    let molecular = vec![0.03; mesh.tris.len()];
    let model = PlanarSgs::new(&system, 0.5, &widths, &molecular).unwrap();
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-12,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls::default();
    let previous = vec![0.0; system.n_u];
    let snapshot = previous.clone();
    let no_query = |_: [f64; 2], _: [f64; 2]| panic!("boundary refusal must precede samples");
    let zero = |_: [f64; 2]| [0.0; 2];
    let refused: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| NavierStokesBoundary::Open {
            edge,
            traction: &no_query,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    assert!(matches!(
        model.solve_steady_mixed_newton(params, &zero, &refused, &previous, controls),
        Err(PlanarSgsError::Newton(FluxNewtonError::Boundary(
            StokesBoundaryError::NoVelocityBoundary
        )))
    ));
    assert!(matches!(
        model.bdf2_mixed_newton_variable_step(
            params,
            &zero,
            &refused,
            &previous,
            Some((&previous, 0.1)),
            0.21,
            controls
        ),
        Err(PlanarSgsError::Newton(FluxNewtonError::History(_)))
    ));
    for closure in [SgsModel::Wale { cw: 0.5 }, SgsModel::Vreman { cv: 0.07 }] {
        let other = PlanarSgs::with_model(&system, closure, &widths, &molecular).unwrap();
        assert!(matches!(
            other.bdf1_mixed_newton_step(params, &zero, &refused, &previous, 0.1, controls),
            Err(PlanarSgsError::InvalidInput(_))
        ));
    }
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| NavierStokesBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    // Finite physical inertia supplies the absent steady velocity anchor.
    // Uniform acceleration has zero strain and conservative convective divergence.
    let result = model
        .bdf1_mixed_newton_step(
            params,
            &|_| [2.0, -1.0],
            &boundaries,
            &previous,
            0.1,
            controls,
        )
        .unwrap();
    let (error, divergence) =
        system.velocity_error(&result.step.solution.field.x, &|_| [0.2, -0.1]);
    assert!(error < 1e-10 && divergence < 1e-10, "{error}, {divergence}");
    assert!(result.sgs.integrated_dissipation_m4_s3 < 1e-25);
    assert_eq!(previous, snapshot);
}
