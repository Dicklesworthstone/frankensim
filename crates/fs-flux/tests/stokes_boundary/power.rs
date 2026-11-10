//! Independent continuum and polynomial oracles for native source-work traces.

use super::*;
use fs_flux::bdm::eval_basis;
use fs_flux::{NavierStokesBoundary, OpenBoundaryPolicy};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}"
    );
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[test]
fn nonrigid_prescribed_work_keeps_pressure_covariance_and_explicit_gauge() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let nu = 1.3;
    let reads = Cell::new(0);
    let velocity = |p: [f64; 2]| {
        reads.set(reads.get() + 1);
        [
            0.2 * p[0] + 0.3 * p[1] + 0.4,
            -0.1 * p[0] - 0.2 * p[1] + 0.1,
        ]
    };
    let boundaries = velocity_boundaries(&mesh, &velocity);
    let edges: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    for gradient in [[0.0; 2], [0.7, -0.4]] {
        let solved = system
            .assemble_stokes_mixed(
                FluxParams {
                    nu,
                    tol: 1e-12,
                    ..FluxParams::default()
                },
                &|_| gradient,
                &boundaries,
            )
            .unwrap()
            .solve()
            .unwrap();
        let (error, divergence) = system.velocity_error(&solved.field.x, &velocity);
        assert!(error < 2e-9 && divergence < 2e-9);
        let before = reads.get();
        let mut missed_covariance = 0.0_f64;
        for offset in [0.0, 7.0] {
            for face in solved.reactions.edges() {
                let e = &mesh.edges[face.edge];
                let a = mesh.verts[e.verts.0];
                let b = mesh.verts[e.verts.1];
                let delta = [b[0] - a[0], b[1] - a[1]];
                let u = [
                    0.2 * e.mid[0] + 0.3 * e.mid[1] + 0.4,
                    -0.1 * e.mid[0] - 0.2 * e.mid[1] + 0.1,
                ];
                let du = [
                    0.2 * delta[0] + 0.3 * delta[1],
                    -0.1 * delta[0] - 0.2 * delta[1],
                ];
                let pressure = offset + dot(gradient, [e.mid[0] - 0.5, e.mid[1] - 0.5]);
                let traction = [
                    (0.4 * nu - pressure) * e.normal[0] + 0.2 * nu * e.normal[1],
                    0.2 * nu * e.normal[0] + (-0.4 * nu - pressure) * e.normal[1],
                ];
                // Exact integral of two affine factors, independently of the
                // rigid wrench and of the native row/edge quadrature routines.
                let covariance = -e.len * dot(gradient, delta) * dot(e.normal, du) / 12.0;
                let expected = e.len * dot(traction, u) + covariance;
                missed_covariance = missed_covariance.max(covariance.abs());
                let power = face.power(Some(offset)).unwrap();
                near(power.prescribed_velocity, expected, 3e-9);
                assert_eq!(power.imposed_traction, 0.0);
                assert_eq!(power.backflow_resistance, 0.0);
                assert_eq!(power.pressure_offsets, vec![(0, Some(offset))]);
                let base = face.power(Some(0.0)).unwrap();
                near(
                    power.prescribed_velocity - base.prescribed_velocity,
                    -offset * e.len * dot(e.normal, u),
                    3e-9,
                );
                assert!(face.power(None).is_err());
                assert!(face.power(Some(f64::NAN)).is_err());
            }
            let total = solved.reactions.power(&edges, &[Some(offset)]).unwrap();
            let body = dot(gradient, [0.65, -0.05]);
            near(total.prescribed_velocity + body, nu * 0.20, 3e-9);
            if gradient == [0.0; 2] {
                assert!(total.prescribed_velocity > 0.0);
            }
        }
        if gradient != [0.0; 2] {
            assert!(missed_covariance > 1e-4);
        }
        assert_eq!(
            reads.get(),
            before,
            "work must use original admitted samples"
        );
        assert!(solved.reactions.power(&[], &[Some(0.0)]).is_err());
        assert!(
            solved
                .reactions
                .power(&[edges[0], edges[0]], &[Some(0.0)])
                .is_err()
        );
        assert!(solved.reactions.power(&edges, &[]).is_err());
    }
}

#[test]
fn finite_slip_wall_input_equals_bulk_plus_relative_slip_loss_not_twice_heat() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let (nu, beta, offset) = (0.7, 1.1, 0.2);
    for shear in [-0.4, 0.4] {
        let velocity = |p: [f64; 2]| [offset + shear * p[1], 0.0];
        let wall = |_| [offset + shear + nu * shear / beta, 0.0];
        let traction = |_: [f64; 2], n: [f64; 2]| {
            [
                nu * shear * n[1] - 3.0 * n[0],
                nu * shear * n[0] - 3.0 * n[1],
            ]
        };
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| {
                let n = mesh.edges[edge].normal;
                if n[1] > 0.5 {
                    StokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_m_s: beta,
                    }
                } else if n[1] < -0.5 {
                    StokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                } else {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                }
            })
            .collect();
        let solved = system
            .assemble_stokes_mixed(
                FluxParams {
                    nu,
                    tol: 1e-12,
                    ..FluxParams::default()
                },
                &|_| [0.0; 2],
                &boundaries,
            )
            .unwrap()
            .solve()
            .unwrap();
        assert!(system.velocity_error(&solved.field.x, &velocity).0 < 2e-9);
        let edges: Vec<_> = solved
            .reactions
            .edges()
            .iter()
            .map(|face| face.edge)
            .collect();
        let power = solved.reactions.power(&edges, &[None]).unwrap();
        let slip: f64 = solved
            .reactions
            .edges()
            .iter()
            .map(|face| face.slip_dissipation())
            .sum();
        let expected_slip = nu * nu * shear * shear / beta;
        near(slip, expected_slip, 2e-10);
        near(
            power.prescribed_velocity,
            nu * shear * shear + expected_slip,
            3e-9,
        );
        near(power.imposed_traction, 0.0, 2e-9);
        assert_eq!(power.backflow_resistance, 0.0);
        for face in solved.reactions.edges() {
            let e = &mesh.edges[face.edge];
            if e.normal[1] > 0.5 {
                near(
                    face.power(None).unwrap().prescribed_velocity,
                    e.len * nu * shear * wall(e.mid)[0],
                    2e-9,
                );
            }
            assert!(face.power(Some(0.0)).is_err());
        }
    }
}

#[test]
fn reversing_frozen_advector_backflow_uses_actual_split_cubic_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let w = |p: [f64; 2]| [p[1] - 0.5, 0.0];
    let advect: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|e| {
            let a = w(mesh.verts[e.verts.0]);
            let b = w(mesh.verts[e.verts.1]);
            [
                0.5 * dot([a[0] + b[0], a[1] + b[1]], e.normal),
                dot([b[0] - a[0], b[1] - a[1]], e.normal) / 12.0,
            ]
        })
        .collect();
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    for policy in [
        OpenBoundaryPolicy::NaturalTraction,
        OpenBoundaryPolicy::DissipativeBackflow,
    ] {
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&e| mesh.is_boundary(e))
            .map(|edge| {
                if mesh.edges[edge].normal[0].abs() > 0.5 {
                    NavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: policy,
                    }
                } else {
                    NavierStokesBoundary::Velocity { edge, value: &zero }
                }
            })
            .collect();
        let solved = system
            .assemble_navier_stokes_mixed_bdf1(
                FluxParams {
                    nu: 0.4,
                    tol: 1e-12,
                    ..FluxParams::default()
                },
                &|_| [1.0, 0.0],
                &boundaries,
                &advect,
                &vec![0.0; system.n_u],
                0.04,
            )
            .unwrap()
            .solve()
            .unwrap();
        let at = |cell: usize, point: [f64; 2]| {
            let mut u = [0.0; 2];
            for local in 0..6 {
                let basis = eval_basis(&system.bases[cell], local, point);
                let dof = 2 * mesh.tri_edges[cell][local / 2].0 + local % 2;
                for axis in 0..2 {
                    u[axis] += solved.field.x[dof] * basis[axis];
                }
            }
            u
        };
        let mut expected = 0.0;
        for face in solved.reactions.edges() {
            let e = &mesh.edges[face.edge];
            if e.normal[0].abs() < 0.5 {
                continue;
            }
            let a = mesh.verts[e.verts.0];
            let b = mesh.verts[e.verts.1];
            let u = at(e.tris.0, a);
            let end = at(e.tris.0, b);
            let du = [end[0] - u[0], end[1] - u[1]];
            let n0 = dot(w(a), e.normal);
            let dn = dot(w(b), e.normal) - n0;
            assert!(n0 * (n0 + dn) < 0.0);
            let split = -n0 / dn;
            let (lower, upper) = if n0 < 0.0 { (0.0, split) } else { (split, 1.0) };
            // Analytically integrate the cubic (w.n)*|u_h|² on its negative
            // half-face, without calling the native split quadrature.
            let c = [dot(u, u), 2.0 * dot(u, du), dot(du, du)];
            let coefficients = [
                n0 * c[0],
                n0 * c[1] + dn * c[0],
                n0 * c[2] + dn * c[1],
                dn * c[2],
            ];
            let integral: f64 = coefficients
                .iter()
                .enumerate()
                .map(|(i, &value)| {
                    let degree = (i + 1) as i32;
                    value * (upper.powi(degree) - lower.powi(degree)) / f64::from(degree)
                })
                .sum();
            expected += e.len * integral;
        }
        assert!(expected < -1e-8);
        let edges: Vec<_> = solved
            .reactions
            .edges()
            .iter()
            .map(|face| face.edge)
            .collect();
        let power = solved.reactions.power(&edges, &[None]).unwrap();
        near(power.prescribed_velocity, 0.0, 1e-12);
        assert_eq!(power.imposed_traction, 0.0);
        near(
            power.backflow_resistance,
            if policy == OpenBoundaryPolicy::DissipativeBackflow {
                expected
            } else {
                0.0
            },
            2e-12,
        );
    }
}
