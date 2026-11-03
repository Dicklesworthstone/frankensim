//! Nonlinear pressure/traction flow and independently integrated boundary work.

use super::*;
use fs_flux::bdm::edge_gauss_pub;
use fs_flux::{NavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy};

#[path = "reservoir.rs"]
mod reservoir;
#[path = "navier_stokes_rheology.rs"]
mod rheology;

fn moments(mesh: &TriMesh, velocity: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut values = Vec::new();
    for edge in &mesh.edges {
        let a = mesh.verts[edge.verts.0];
        let b = mesh.verts[edge.verts.1];
        let mut pair = [0.0; 2];
        for (point, weight) in edge_gauss_pub(a, b) {
            let u = velocity(point);
            let normal = u[0] * edge.normal[0] + u[1] * edge.normal[1];
            let s = ((point[0] - a[0]) * (b[0] - a[0]) + (point[1] - a[1]) * (b[1] - a[1]))
                / edge.len.powi(2)
                - 0.5;
            pair[0] += weight * normal / edge.len;
            pair[1] += weight * normal * s / edge.len;
        }
        values.extend(pair);
    }
    values
}

fn open<'a>(
    mesh: &TriMesh,
    traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
    backflow: OpenBoundaryPolicy,
) -> Vec<NavierStokesBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| NavierStokesBoundary::Open {
            edge,
            traction,
            backflow,
        })
        .collect()
}

#[test]
fn nonlinear_pressure_flow_has_live_convection_and_keeps_its_datum() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let dt = 0.02;
    let old = |p: [f64; 2]| [0.1 * p[0] + 0.1 * p[1] + 0.2, -0.1 * p[1] + 0.05];
    let target = |p: [f64; 2]| [0.3 * p[0] + 0.2 * p[1] + 0.4, -0.3 * p[1] + 0.1];
    let previous = moments(&mesh, &old);
    let snapshot = previous.clone();
    let viscosity = |_: usize, p: [f64; 2]| 0.5 + 0.1 * p[0];
    // grad(u) = [[.3,.2],[0,-.3]], so u.grad(u) =
    // [.09*x+.14, .09*y-.03], and div(2 nu D(u))=[.06,.02].
    let force = |p: [f64; 2]| {
        let u = target(p);
        let v = old(p);
        [
            (u[0] - v[0]) / dt + 0.09 * p[0] + 0.14 - 0.06 + 2.0,
            (u[1] - v[1]) / dt + 0.09 * p[1] - 0.03 - 0.02 - 3.0,
        ]
    };
    let mut velocities = Vec::new();
    for shift in [0.0, 17.0] {
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let pressure = 5.0 + shift + 2.0 * p[0] - 3.0 * p[1];
            let nu = viscosity(0, p);
            [
                (0.6 * nu - pressure) * n[0] + 0.2 * nu * n[1],
                0.2 * nu * n[0] + (-0.6 * nu - pressure) * n[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, e)| {
                if e.normal[0].abs() > 0.5 {
                    NavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                } else {
                    NavierStokesBoundary::Velocity {
                        edge,
                        value: &target,
                    }
                }
            })
            .collect();
        let solved = system
            .bdf1_navier_stokes_mixed_variable_viscosity_step(
                params,
                &viscosity,
                &force,
                &boundaries,
                &previous,
                dt,
                NavierStokesControls {
                    relative_tolerance: 1e-10,
                    ..NavierStokesControls::default()
                },
            )
            .unwrap();
        assert_eq!(previous, snapshot);
        let (error, divergence) = system.velocity_error(&solved.solution.field.x, &target);
        assert!(error < 3e-10 && divergence < 3e-10, "{error}, {divergence}");
        assert!(solved.iterations > 1);
        assert!(solved.nonlinear_residual < 1e-10);
        assert_eq!(
            solved.solution.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        assert!(
            solved
                .open_boundary_policies
                .iter()
                .all(|(_, p)| *p == OpenBoundaryPolicy::NaturalTraction)
        );
        for (actual, point) in solved.solution.pressure.iter().zip(&mesh.centroids) {
            assert!((actual - (5.0 + shift + 2.0 * point[0] - 3.0 * point[1])).abs() < 2e-8);
        }
        // A separate endpoint assembly must accept this velocity, not the
        // previous (lagged) advection field or just the last linear system.
        let endpoint = system
            .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
                params,
                &viscosity,
                &force,
                &boundaries,
                &solved.solution.field.x[..system.n_u],
                &previous,
                dt,
            )
            .unwrap();
        let mut residual = vec![0.0; system.n];
        endpoint
            .matrix()
            .spmv(&solved.solution.field.x, &mut residual);
        let norm = residual
            .iter()
            .zip(endpoint.rhs())
            .fold(0.0_f64, |r, (ax, b)| r.hypot(ax - b));
        assert!(norm < 1e-8);
        if shift == 0.0 {
            let stokes: Vec<_> = boundaries
                .iter()
                .map(|b| match *b {
                    NavierStokesBoundary::Velocity { edge, value } => {
                        StokesBoundary::Velocity { edge, value }
                    }
                    NavierStokesBoundary::NavierSlip {
                        edge,
                        value,
                        tangential_resistance_m_s,
                    } => StokesBoundary::NavierSlip {
                        edge,
                        value,
                        tangential_resistance_m_s,
                    },
                    NavierStokesBoundary::Open { edge, traction, .. }
                    | NavierStokesBoundary::OpenReservoir { edge, traction, .. } => {
                        StokesBoundary::TotalTraction {
                            edge,
                            value: traction,
                        }
                    }
                })
                .collect();
            let without_convection = system
                .assemble_stokes_mixed_bdf1_variable_viscosity(
                    params, &viscosity, &force, &stokes, &previous, dt,
                )
                .unwrap()
                .solve()
                .unwrap();
            assert!(
                system
                    .velocity_error(&without_convection.field.x, &target)
                    .0
                    > 1e-4
            );
        }
        velocities.push(solved.solution.field.x[..system.n_u].to_vec());
    }
    for (a, b) in velocities[0].iter().zip(&velocities[1]) {
        assert!((a - b).abs() < 3e-10);
    }
}

#[test]
fn open_convection_has_the_declared_natural_or_dissipative_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let zero = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let stokes: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| StokesBoundary::TotalTraction { edge, value: &zero })
        .collect();
    let reversing_work = (0.37_f64.powi(4) + 0.63_f64.powi(4)) / 4.0;
    for (reversing, natural_work, dissipative_work) in
        [(false, -2.0, 20.0 / 3.0), (true, 0.0, reversing_work)]
    {
        let velocity = |p: [f64; 2]| {
            if reversing {
                [p[1] - 0.37, 0.0]
            } else {
                [2.0 - p[0], p[1]]
            }
        };
        let history = moments(&mesh, &velocity);
        let mut state = history.clone();
        state.resize(system.n, 0.0);
        let base = system
            .assemble_stokes_mixed_bdf1(params, &|_| [0.0; 2], &stokes, &history, 0.1)
            .unwrap();
        let mut base_action = vec![0.0; system.n];
        base.matrix().spmv(&state, &mut base_action);
        for (policy, expected) in [
            (OpenBoundaryPolicy::NaturalTraction, natural_work),
            (OpenBoundaryPolicy::DissipativeBackflow, dissipative_work),
        ] {
            let boundaries = open(&mesh, &zero, policy);
            let full = system
                .assemble_navier_stokes_mixed_bdf1(
                    params,
                    &|_| [0.0; 2],
                    &boundaries,
                    &history,
                    &history,
                    0.1,
                )
                .unwrap();
            let mut action = vec![0.0; system.n];
            full.matrix().spmv(&state, &mut action);
            let measured: f64 = state
                .iter()
                .zip(action.iter().zip(&base_action))
                .map(|(u, (a, b))| u * (a - b))
                .sum();
            // Analytic surface integrals for the two continuous affine fields:
            // natural: int (u.n)|u|²/2; dissipative: int |u.n||u|²/2.
            // The second field reverses INSIDE a vertical edge at y=.37.
            assert!(
                (measured - expected).abs() < 2e-11,
                "{policy:?}: {measured} vs {expected}"
            );
            assert_eq!(full.rhs(), base.rhs());
        }
    }
}

#[test]
fn closed_oseen_assembly_reuses_the_original_operator_exactly() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let g = |p: [f64; 2]| [p[1] - 0.5, 0.0];
    let f = |p: [f64; 2]| [p[0], -p[1]];
    let viscosity = |_: usize, p: [f64; 2]| 0.5 + p[0];
    let previous = moments(&mesh, &g);
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| NavierStokesBoundary::Velocity { edge, value: &g })
        .collect();
    let mixed = system
        .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
            params,
            &viscosity,
            &f,
            &boundaries,
            &previous,
            &previous,
            0.1,
        )
        .unwrap();
    let (original, rhs) = system.assemble_viscous_transient(
        params,
        &viscosity,
        &f,
        &g,
        Some(&previous),
        10.0,
        Some(&previous),
    );
    assert_eq!(mixed.matrix().to_dense(), original.to_dense());
    assert_eq!(mixed.rhs(), rhs.as_slice());
}

#[test]
fn invalid_or_unconverged_nonlinear_steps_never_publish_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let calls = Cell::new(0);
    let traction = |p: [f64; 2], n: [f64; 2]| {
        calls.set(calls.get() + 1);
        [-(3.0 - p[0]) * n[0], -(3.0 - p[0]) * n[1]]
    };
    let boundaries = open(&mesh, &traction, OpenBoundaryPolicy::NaturalTraction);
    let previous = vec![0.0; system.n_u];
    for dt in [0.0, -1.0, f64::NAN, f64::from_bits(1)] {
        assert!(
            system
                .bdf1_navier_stokes_mixed_step(
                    params,
                    &|_| [0.0; 2],
                    &boundaries,
                    &previous,
                    dt,
                    NavierStokesControls::default()
                )
                .is_err()
        );
    }
    for controls in [
        NavierStokesControls {
            maximum_iterations: 0,
            ..NavierStokesControls::default()
        },
        NavierStokesControls {
            relative_tolerance: f64::NAN,
            ..NavierStokesControls::default()
        },
        NavierStokesControls {
            relaxation: 0.0,
            ..NavierStokesControls::default()
        },
        NavierStokesControls {
            relaxation: 1.01,
            ..NavierStokesControls::default()
        },
    ] {
        assert!(
            system
                .bdf1_navier_stokes_mixed_step(
                    params,
                    &|_| [0.0; 2],
                    &boundaries,
                    &previous,
                    0.1,
                    controls
                )
                .is_err()
        );
    }
    assert_eq!(calls.get(), 0);
    assert!(matches!(
        system.bdf1_navier_stokes_mixed_step(
            params,
            &|_| [0.0; 2],
            &boundaries,
            &previous,
            0.1,
            NavierStokesControls {
                maximum_iterations: 1,
                ..NavierStokesControls::default()
            }
        ),
        Err(StokesBoundaryError::Rheology(
            fs_flux::FluxRheologyError::NonlinearConvergence { .. }
        ))
    ));
    assert!(calls.get() > 0);
    assert!(previous.iter().all(|v| *v == 0.0));
}
