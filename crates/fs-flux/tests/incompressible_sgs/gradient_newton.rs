//! Full-gradient constitutive derivatives against the independent primal owner.
use super::*;
use fs_flux::ns::linearization::{FluxLinearizationError, FluxViscosityGradientTangent};
use fs_lbm::sgs3::eddy_viscosity_tangent;
use std::cell::Cell;

fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |v, x| v.hypot(*x))
}

fn response(
    model: SgsModel,
    width: f64,
    molecular: f64,
    gradient: [[f64; 2]; 2],
) -> FluxViscosityGradientTangent {
    let value = eddy_viscosity_tangent(
        model,
        width,
        [
            [gradient[0][0], gradient[0][1], 0.0],
            [gradient[1][0], gradient[1][1], 0.0],
            [0.0; 3],
        ],
    )
    .unwrap();
    FluxViscosityGradientTangent {
        nu: molecular + value.viscosity,
        derivative: core::array::from_fn(|i| core::array::from_fn(|j| value.derivative[i][j])),
    }
}

#[test]
fn full_gradient_jacobian_matches_primal_sip_boundary_and_reversing_convection() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let widths: Vec<_> = (0..mesh.tris.len())
        .map(|i| 0.4 + 0.1 * (i % 3) as f64)
        .collect();
    let molecular: Vec<_> = (0..mesh.tris.len())
        .map(|i| 0.03 + 0.01 * (i % 4) as f64)
        .collect();
    let params = FluxParams {
        nu: f64::NAN,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let target = |p: [f64; 2]| [0.61 - 1.4 * p[1], 0.2 * p[0] - 0.11];
    let reservoir = |p: [f64; 2]| [-0.37 + 0.1 * p[1], 0.29];
    let traction = |_: [f64; 2], n: [f64; 2]| n.map(|v| -1.3 * v);
    let force = |_: [f64; 2]| [0.13, -0.07];
    let previous = moments(&mesh, &|_| [0.12, -0.03]);
    let mut state = moments(&mesh, &target);
    // Retain reversed open flux while making DG jumps and cell gradients live.
    for (i, value) in state.iter_mut().enumerate() {
        *value += 0.015 * ((i % 7) as f64 - 3.0);
    }
    state.resize(system.n, 0.23);
    let general: Vec<_> = (0..system.n)
        .map(|i| ((i * 13 % 29) as f64 - 14.0) / 31.0)
        .collect();
    // A strain-only derivative would miss the constitutive response to this
    // pure rotational gradient variation. Boundary data remain fixed.
    let mut rotation = moments(&mesh, &|p| [-p[1], p[0]]);
    rotation.resize(system.n, 0.0);
    let dt = 0.4;
    for model in [SgsModel::Wale { cw: 0.7 }, SgsModel::Vreman { cv: 0.4 }] {
        let closure = PlanarSgs::with_model(&system, model, &widths, &molecular).unwrap();
        for kind in 0usize..8 {
            let policy = (kind.saturating_sub(2)) % 3;
            let slip = kind >= 5;
            let convection = kind != 0;
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, face)| {
                    if kind < 2 || face.normal[0].abs() < 0.5 && !slip {
                        NavierStokesBoundary::Velocity {
                            edge,
                            value: &target,
                        }
                    } else if face.normal[0].abs() < 0.5 {
                        NavierStokesBoundary::NavierSlip {
                            edge,
                            value: &target,
                            tangential_resistance_m_s: 0.7,
                        }
                    } else if policy == 2 {
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
                })
                .collect();
            let samples = Cell::new(0usize);
            let law = |cell: usize, gradient| {
                samples.set(samples.get() + 1);
                response(model, widths[cell], molecular[cell], gradient)
            };
            let linearize = |law: &dyn Fn(usize, [[f64; 2]; 2]) -> FluxViscosityGradientTangent| {
                if kind < 2 {
                    system
                        .linearize_gradient_rheology(
                            params,
                            law,
                            &force,
                            &target,
                            &state,
                            convection,
                            1.0 / dt,
                            Some(&previous),
                        )
                        .unwrap()
                } else {
                    system
                        .linearize_navier_stokes_mixed_gradient_rheology(
                            params,
                            law,
                            &force,
                            &boundaries,
                            &state,
                            1.0 / dt,
                            Some(&previous),
                        )
                        .unwrap()
                }
            };
            let linearization = linearize(&law);
            assert_eq!(
                samples.get(),
                mesh.tris.len(),
                "one closure per cell per assembly"
            );
            let frozen = linearize(&|cell, gradient| FluxViscosityGradientTangent {
                derivative: [[0.0; 2]; 2],
                ..response(model, widths[cell], molecular[cell], gradient)
            });
            for direction in [&general, &rotation] {
                let mut exact = vec![0.0; system.n];
                let mut incomplete = vec![0.0; system.n];
                linearization.jacobian.spmv(direction, &mut exact);
                frozen.jacobian.spmv(direction, &mut incomplete);
                let epsilon = 1e-6;
                let residual = |sign: f64| {
                    let x: Vec<_> = state
                        .iter()
                        .zip(direction)
                        .map(|(x, d)| x + sign * epsilon * d)
                        .collect();
                    // Values come from the ORIGINAL primal closure/diagnostic
                    // path, independently of the new derivative callback.
                    let field = closure.field(&x[..system.n_u]).unwrap();
                    let viscosity =
                        |cell: usize, _: [f64; 2]| field.cells[cell].effective_viscosity_m2_s;
                    let mut result = vec![0.0; system.n];
                    if kind < 2 {
                        let (matrix, rhs) = system.assemble_viscous_transient(
                            params,
                            &viscosity,
                            &force,
                            &target,
                            convection.then_some(&x[..system.n_u]),
                            1.0 / dt,
                            Some(&previous),
                        );
                        matrix.spmv(&x, &mut result);
                        for (r, b) in result.iter_mut().zip(rhs) {
                            *r -= b;
                        }
                    } else {
                        let assembly = system
                            .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                                params,
                                &viscosity,
                                &force,
                                &boundaries,
                                &x[..system.n_u],
                                &previous,
                                dt,
                            )
                            .unwrap();
                        assembly.matrix().spmv(&x, &mut result);
                        for (r, b) in result.iter_mut().zip(assembly.rhs()) {
                            *r -= b;
                        }
                    }
                    result
                };
                let plus = residual(1.0);
                let minus = residual(-1.0);
                let defect: Vec<_> = exact
                    .iter()
                    .zip(plus.iter().zip(&minus))
                    .map(|(a, (p, m))| a - (p - m) / (2.0 * epsilon))
                    .collect();
                let relative = norm(&defect) / norm(&exact);
                assert!(
                    relative < 3e-7,
                    "{model:?}, boundary kind {kind}: {relative:e}"
                );
                let omitted: Vec<_> = exact.iter().zip(&incomplete).map(|(a, b)| a - b).collect();
                assert!(
                    norm(&omitted) > 1e-5,
                    "full-gradient tangent negative control must be live"
                );
            }
        }
    }
}

#[test]
fn gradient_newton_resolves_active_stress_from_zero_and_matches_original_endpoint() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let widths = vec![0.7; mesh.tris.len()];
    let molecular = vec![0.01; mesh.tris.len()];
    let params = FluxParams {
        nu: f64::NAN,
        tol: 1e-11,
        sigma: 15.0,
        ..FluxParams::default()
    };
    let controls = FluxNewtonControls {
        maximum_iterations: 40,
        maximum_backtracks: 20,
        residual_tolerance: 1e-9,
    };
    let force = |p: [f64; 2]| [20.0 * (p[1] - 0.5), -20.0 * (p[0] - 0.5)];
    let zero = |_: [f64; 2]| [0.0; 2];
    let previous = vec![0.0; system.n_u];
    let snapshot = previous.clone();
    let dt = 0.2;
    for model in [SgsModel::Wale { cw: 0.5 }, SgsModel::Vreman { cv: 0.2 }] {
        let closure = PlanarSgs::with_model(&system, model, &widths, &molecular).unwrap();
        let result = closure
            .bdf1_newton_step(params, &force, &zero, &previous, dt, controls)
            .unwrap();
        let x = &result.step.solution.x;
        assert!(result.sgs.integrated_dissipation_m4_s3 > 0.0);
        for active in [true, false] {
            let (matrix, rhs) = system.assemble_viscous_transient(
                params,
                &|cell, _| {
                    if active {
                        result.sgs.cells[cell].effective_viscosity_m2_s
                    } else {
                        molecular[cell]
                    }
                },
                &force,
                &zero,
                Some(&x[..system.n_u]),
                dt.recip(),
                Some(&previous),
            );
            let mut residual = vec![0.0; system.n];
            matrix.spmv(x, &mut residual);
            for (r, b) in residual.iter_mut().zip(&rhs) {
                *r -= b;
            }
            let defect = norm(&residual) / norm(&rhs);
            if active {
                assert!(defect < 1e-9, "{model:?}: {defect:e}");
            } else {
                assert!(
                    defect > 1e-3,
                    "zero-SGS initial law must fail at the endpoint"
                );
            }
        }
        eprintln!(
            "{model:?} Newton: {} iterations, {} backtracks, residual {:e}",
            result.step.iterations, result.step.backtracks, result.step.nonlinear_residual
        );
        assert_eq!(previous, snapshot);
    }
}

#[test]
fn nonfinite_full_gradient_tangents_refuse_before_primal_callbacks() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let state = vec![0.0; system.n];
    let never = |_: [f64; 2]| panic!("invalid tangent must refuse before external loads");
    for tangent in [
        FluxViscosityGradientTangent {
            nu: 0.0,
            derivative: [[0.0; 2]; 2],
        },
        FluxViscosityGradientTangent {
            nu: 0.03,
            derivative: [[f64::NAN; 2]; 2],
        },
        FluxViscosityGradientTangent {
            nu: f64::INFINITY,
            derivative: [[0.0; 2]; 2],
        },
    ] {
        assert!(matches!(
            system.linearize_gradient_rheology(
                FluxParams::default(),
                &|_, _| tangent,
                &never,
                &never,
                &state,
                true,
                0.0,
                None
            ),
            Err(FluxLinearizationError::GradientConstitutive { cell: 0, .. })
        ));
    }
    assert_eq!(state, vec![0.0; system.n]);
}
