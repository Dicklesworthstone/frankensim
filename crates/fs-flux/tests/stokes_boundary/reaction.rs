//! Physical reaction oracles independent of the recovered normal rows.

use super::*;
use fs_flux::bdm::{edge_gauss_pub, eval_basis, tri_quad};
use fs_flux::{
    BoundaryReactionOwner, NavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}"
    );
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

#[test]
fn static_box_retains_pressure_datum_in_wall_reactions() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], n: [f64; 2]| [-7.0 * n[0], -7.0 * n[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, geometry)| {
            if geometry.normal[0] > 0.5 {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            } else {
                StokesBoundary::Velocity { edge, value: &zero }
            }
        })
        .collect();
    let solved = system
        .assemble_stokes_mixed(FluxParams::default(), &zero, &boundaries)
        .unwrap()
        .solve()
        .unwrap();
    for item in solved.reactions.edges() {
        let edge = &mesh.edges[item.edge];
        let load = item.load([0.2, 0.3], None).unwrap();
        let exact = edge.normal.map(|n| 7.0 * edge.len * n);
        near(load.stress.force[0], exact[0], 3e-10);
        near(load.stress.force[1], exact[1], 3e-10);
        near(
            load.stress.moment,
            cross([edge.mid[0] - 0.2, edge.mid[1] - 0.3], exact),
            3e-10,
        );
        assert_eq!(load.pressure_offsets, vec![(0, None)]);
        assert_eq!(load.outward_momentum.force, [0.0; 2]);
        assert!(
            item.load([0.0; 2], Some(7.0)).is_err(),
            "traction datum cannot be silently shifted"
        );
    }
}

#[test]
fn affine_strain_restores_normal_and_shear_stress_without_reporting_queries() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let nu = 1.3;
    let stress = |n: [f64; 2]| {
        [
            (0.4 * nu - 7.0) * n[0] + 0.3 * nu * n[1],
            0.3 * nu * n[0] + (-0.4 * nu - 7.0) * n[1],
        ]
    };
    let traction = |_: [f64; 2], n: [f64; 2]| stress(n);
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, geometry)| {
            if geometry.normal[0].abs() > 0.5 {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            } else {
                StokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            }
        })
        .collect();
    let queries = Cell::new(0);
    let viscosity = |_: usize, _: [f64; 2]| {
        queries.set(queries.get() + 1);
        nu
    };
    let assembled = system
        .assemble_stokes_mixed_variable_viscosity(
            FluxParams::default(),
            &viscosity,
            &|_| [0.0; 2],
            &boundaries,
        )
        .unwrap();
    let before = queries.get();
    let solved = assembled.solve().unwrap();
    assert!(system.velocity_error(&solved.field.x, &velocity).0 < 2e-10);
    let origin = [0.2, -0.3];
    for item in solved.reactions.edges() {
        let edge = &mesh.edges[item.edge];
        let load = item.load(origin, None).unwrap();
        let exact = stress(edge.normal).map(|t| -edge.len * t);
        near(load.stress.force[0], exact[0], 3e-9);
        near(load.stress.force[1], exact[1], 3e-9);
        let moment = cross([edge.mid[0] - origin[0], edge.mid[1] - origin[1]], exact);
        near(load.stress.moment, moment, 3e-9);
        near(
            load.stress.rigid_power([0.7, -0.2], 0.4).unwrap(),
            0.7 * exact[0] - 0.2 * exact[1] + 0.4 * moment,
            3e-9,
        );
        if edge.normal[1].abs() > 0.5 {
            assert_eq!(item.owner, BoundaryReactionOwner::VelocityConstraint);
        }
    }
    assert_eq!(
        queries.get(),
        before,
        "solving/reporting must reuse admitted viscosity samples"
    );
}

#[test]
fn hydrostatic_edge_moments_use_explicit_independent_pressure_references() {
    let mesh = separate_domains();
    let system = FluxSystem::new(&mesh);
    let force = |p: [f64; 2]| if p[0] < 2.0 { [2.0, -3.0] } else { [-4.0, 1.0] };
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries = velocity_boundaries(&mesh, &zero);
    let solved = system
        .assemble_stokes_mixed(FluxParams::default(), &force, &boundaries)
        .unwrap()
        .solve()
        .unwrap();
    let offsets = [Some(5.0), Some(11.0)];
    let origin = [0.2, 0.3];
    let mut naive_error = 0.0_f64;
    for item in solved.reactions.edges() {
        let edge = &mesh.edges[item.edge];
        let a = mesh.verts[edge.verts.0];
        let b = mesh.verts[edge.verts.1];
        let delta = [b[0] - a[0], b[1] - a[1]];
        let gradient = force(edge.mid);
        let center = if item.component == 0 {
            [0.5, 0.5]
        } else {
            [3.5, 0.5]
        };
        let p_mid = offsets[item.component].unwrap()
            + gradient[0] * (edge.mid[0] - center[0])
            + gradient[1] * (edge.mid[1] - center[1]);
        let exact = edge.normal.map(|n| edge.len * p_mid * n);
        let moment = cross([edge.mid[0] - origin[0], edge.mid[1] - origin[1]], exact)
            + edge.len
                * (gradient[0] * delta[0] + gradient[1] * delta[1])
                * cross(delta, edge.normal)
                / 12.0;
        let load = item.load(origin, offsets[item.component]).unwrap();
        near(load.stress.force[0], exact[0], 3e-9);
        near(load.stress.force[1], exact[1], 3e-9);
        near(load.stress.moment, moment, 3e-9);
        assert!(item.load(origin, None).is_err());
        let naive = solved.pressure[edge.tris.0] + offsets[item.component].unwrap();
        naive_error = naive_error.max(edge.len * (naive - p_mid).abs());
    }
    assert!(
        naive_error > 0.1,
        "cell-pressure trace must expose the incorrect force extraction"
    );
    let all: Vec<_> = solved
        .reactions
        .edges()
        .iter()
        .map(|edge| edge.edge)
        .collect();
    let total = solved.reactions.resultant(&all, origin, &offsets).unwrap();
    near(total.stress.force[0], -2.0, 3e-9);
    near(total.stress.force[1], -2.0, 3e-9);
    assert_eq!(
        total.pressure_offsets,
        vec![(0, Some(5.0)), (1, Some(11.0))]
    );
    assert!(
        solved
            .reactions
            .resultant(&[all[0], all[0]], origin, &offsets)
            .is_err()
    );
}

#[test]
fn bdf_momentum_and_moment_close_with_sip_reactions_and_separate_open_flux() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let previous: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|edge| [0.03 * edge.normal[0], 0.0])
        .collect();
    let dt = 0.04;
    let zero = |_: [f64; 2]| [0.0; 2];
    let force = |p: [f64; 2]| [0.2 + 0.1 * p[1], 0.0];
    let traction = |p: [f64; 2], n: [f64; 2]| n.map(|v| -(5.0 - p[0]) * v);
    let law = |_: usize, _: [f64; 2], rate: f64| 0.3 / (1.0 + 0.15 * rate);
    for policy in [
        OpenBoundaryPolicy::NaturalTraction,
        OpenBoundaryPolicy::DissipativeBackflow,
    ] {
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
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
        let step = system
            .bdf1_navier_stokes_mixed_rheology_step(
                FluxParams::default(),
                &law,
                &force,
                &boundaries,
                &previous,
                dt,
                NavierStokesControls {
                    relative_tolerance: 1e-10,
                    relaxation: 0.8,
                    ..NavierStokesControls::default()
                },
            )
            .unwrap();
        let all: Vec<_> = step
            .solution
            .reactions
            .edges()
            .iter()
            .map(|edge| edge.edge)
            .collect();
        let load = step
            .solution
            .reactions
            .resultant(&all, [0.0; 2], &[None])
            .unwrap();
        // DG tangential traces belong to the actual integration cell; a
        // nearest-centroid point lookup can select a different polynomial.
        let velocity_in_cell = |cell: usize, point: [f64; 2]| {
            let mut velocity = [0.0; 2];
            for i in 0..6 {
                let basis = eval_basis(&system.bases[cell], i, point);
                let global = 2 * mesh.tri_edges[cell][i / 2].0 + i % 2;
                for axis in 0..2 {
                    velocity[axis] += step.solution.field.x[global] * basis[axis];
                }
            }
            velocity
        };
        let mut force_minus_inertia = [0.0; 2];
        let mut moment_minus_inertia = 0.0;
        for (cell, vertices) in mesh.tris.iter().enumerate() {
            let points = vertices.map(|vertex| mesh.verts[vertex]);
            for (point, weight) in tri_quad(points, mesh.areas[cell]) {
                let u = velocity_in_cell(cell, point);
                let f = force(point);
                let residual = [f[0] - (u[0] - 0.03) / dt, f[1] - u[1] / dt];
                for axis in 0..2 {
                    force_minus_inertia[axis] += weight * residual[axis];
                }
                moment_minus_inertia += weight * cross(point, residual);
            }
        }
        let mut flux = [0.0; 2];
        let mut flux_moment = 0.0;
        let mut slip = 0.0_f64;
        for item in step.solution.reactions.edges() {
            let edge = &mesh.edges[item.edge];
            for (point, weight) in
                edge_gauss_pub(mesh.verts[edge.verts.0], mesh.verts[edge.verts.1])
            {
                let u = velocity_in_cell(edge.tris.0, point);
                let wn = u[0] * edge.normal[0] + u[1] * edge.normal[1];
                let momentum = u.map(|v| weight * wn * v);
                for axis in 0..2 {
                    flux[axis] += momentum[axis];
                }
                flux_moment += cross(point, momentum);
                if item.owner == BoundaryReactionOwner::VelocityConstraint {
                    slip = slip.max(u[0].abs());
                }
            }
        }
        assert!(
            slip > 1e-8,
            "non-affine flow must exercise weak SIP boundary consistency"
        );
        for axis in 0..2 {
            near(load.outward_momentum.force[axis], flux[axis], 2e-9);
            near(
                load.stress.force[axis] + flux[axis],
                force_minus_inertia[axis],
                2e-8,
            );
        }
        near(load.outward_momentum.moment, flux_moment, 2e-9);
        near(load.stress.moment + flux_moment, moment_minus_inertia, 2e-8);
        if policy == OpenBoundaryPolicy::NaturalTraction {
            assert_eq!(load.backflow.force, [0.0; 2]);
        } else {
            assert!(load.backflow.force[0] > 0.0);
        }
    }
}
