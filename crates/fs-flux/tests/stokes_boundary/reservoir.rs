//! Incoming reservoir momentum reaches the original Oseen and BDF1 owners.

use super::*;
use fs_flux::bdm::{eval_basis, tri_quad};
use fs_flux::{BoundaryReactionOwner, FluxRheologyError};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}"
    );
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn velocity_in_cell(
    system: &FluxSystem<'_>,
    values: &[f64],
    cell: usize,
    point: [f64; 2],
) -> [f64; 2] {
    let mut velocity = [0.0; 2];
    for i in 0..6 {
        let basis = eval_basis(&system.bases[cell], i, point);
        let global = 2 * system.mesh.tri_edges[cell][i / 2].0 + i % 2;
        for axis in 0..2 {
            velocity[axis] += values[global] * basis[axis];
        }
    }
    velocity
}

// Exact polynomial antiderivatives, independent of the production Gauss rule.
fn product_integral(factors: &[[f64; 2]], lo: f64, hi: f64) -> f64 {
    let mut coefficients = vec![1.0];
    for &[constant, slope] in factors {
        let mut next = vec![0.0; coefficients.len() + 1];
        for (power, value) in coefficients.into_iter().enumerate() {
            next[power] += value * constant;
            next[power + 1] += value * slope;
        }
        coefficients = next;
    }
    coefficients
        .iter()
        .enumerate()
        .map(|(power, coefficient)| {
            let n = (power + 1) as i32;
            coefficient * (hi.powi(n) - lo.powi(n)) / f64::from(n)
        })
        .sum()
}

// Integrate the affine upwind flux and its cubic moment using endpoint values.
fn flux(
    mesh: &TriMesh,
    edge: usize,
    advector: &dyn Fn([f64; 2]) -> [f64; 2],
    interior: &dyn Fn([f64; 2]) -> [f64; 2],
    reservoir: &dyn Fn([f64; 2]) -> [f64; 2],
    origin: [f64; 2],
) -> ([f64; 2], f64) {
    let geometry = &mesh.edges[edge];
    let a = mesh.verts[geometry.verts.0];
    let b = mesh.verts[geometry.verts.1];
    let wa = advector(a);
    let wb = advector(b);
    let normal = [
        wa[0] * geometry.normal[0] + wa[1] * geometry.normal[1],
        (wb[0] - wa[0]) * geometry.normal[0] + (wb[1] - wa[1]) * geometry.normal[1],
    ];
    let mut cuts = vec![0.0, 1.0];
    if normal[1] != 0.0 {
        let zero = -normal[0] / normal[1];
        if zero > 0.0 && zero < 1.0 {
            cuts.insert(1, zero);
        }
    }
    let position: [[f64; 2]; 2] =
        core::array::from_fn(|axis| [a[axis] - origin[axis], b[axis] - a[axis]]);
    let mut force = [0.0; 2];
    let mut moment = 0.0;
    for interval in cuts.windows(2) {
        let (lo, hi) = (interval[0], interval[1]);
        let trace = if normal[0] + normal[1] * (lo + hi) * 0.5 < 0.0 {
            reservoir
        } else {
            interior
        };
        let ua = trace(a);
        let ub = trace(b);
        let u: [[f64; 2]; 2] = core::array::from_fn(|axis| [ua[axis], ub[axis] - ua[axis]]);
        for axis in 0..2 {
            force[axis] += geometry.len * product_integral(&[normal, u[axis]], lo, hi);
        }
        moment += geometry.len
            * (product_integral(&[position[0], normal, u[1]], lo, hi)
                - product_integral(&[position[1], normal, u[0]], lo, hi));
    }
    (force, moment)
}

fn selected_reservoir<'a>(
    mesh: &TriMesh,
    selected: usize,
    traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
    velocity: &'a dyn Fn([f64; 2]) -> [f64; 2],
) -> Vec<NavierStokesBoundary<'a>> {
    (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .map(|edge| {
            if edge == selected {
                NavierStokesBoundary::OpenReservoir {
                    edge,
                    traction,
                    velocity,
                }
            } else {
                NavierStokesBoundary::Open {
                    edge,
                    traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            }
        })
        .collect()
}

#[test]
fn reversing_face_imports_exact_affine_momentum_and_torque_without_fixing_flow_rate() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let left = mesh
        .edges
        .iter()
        .enumerate()
        .find(|(edge, e)| mesh.is_boundary(*edge) && e.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .unwrap();
    let advector = |p: [f64; 2]| [p[1] - 0.3, 0.0];
    let incoming = |p: [f64; 2]| [0.7 + 0.2 * p[1], -0.4 + 0.9 * p[1]];
    let zero = |_: [f64; 2]| [0.0; 2];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let history = moments(&mesh, &|_| [0.04, -0.02]);
    let w = moments(&mesh, &advector);
    let params = FluxParams {
        nu: 0.3,
        ..FluxParams::default()
    };
    let boundaries = selected_reservoir(&mesh, left, &traction, &incoming);
    let empty = selected_reservoir(&mesh, left, &traction, &zero);
    let assembled = system
        .assemble_navier_stokes_mixed_bdf1(params, &zero, &boundaries, &w, &history, 0.03)
        .unwrap();
    let baseline = system
        .assemble_navier_stokes_mixed_bdf1(params, &zero, &empty, &w, &history, 0.03)
        .unwrap();
    for row in 0..system.n {
        assert_eq!(assembled.matrix().row(row), baseline.matrix().row(row));
    }
    let (force, torque) = flux(&mesh, left, &advector, &zero, &incoming, [0.0; 2]);
    let rigid: [&dyn Fn([f64; 2]) -> [f64; 2]; 3] =
        [&|_| [1.0, 0.0], &|_| [0.0, 1.0], &|p| [-p[1], p[0]]];
    for (test, expected) in rigid.into_iter().zip([-force[0], -force[1], -torque]) {
        let work: f64 = moments(&mesh, test)
            .iter()
            .enumerate()
            .map(|(i, v)| v * (assembled.rhs()[i] - baseline.rhs()[i]))
            .sum();
        near(work, expected, 3e-13);
    }
    assert!(force[0].abs() > 0.1);
    let solved = assembled.solve().unwrap();
    let cell = mesh.edges[left].tris.0;
    let interior = |point| velocity_in_cell(&system, &solved.field.x, cell, point);
    let origin = [0.2, -0.4];
    let (expected, moment) = flux(&mesh, left, &advector, &interior, &incoming, origin);
    let reaction = solved
        .reactions
        .edges()
        .iter()
        .find(|r| r.edge == left)
        .unwrap();
    assert_eq!(reaction.owner, BoundaryReactionOwner::OpenReservoir);
    let load = reaction.load(origin, None).unwrap();
    for axis in 0..2 {
        near(load.outward_momentum.force[axis], expected[axis], 3e-12);
    }
    near(load.outward_momentum.moment, moment, 3e-12);
    assert_eq!(load.stress.force, [0.0; 2]);
    assert_eq!(load.backflow.force, [0.0; 2]);
    // The original frozen advector owns volume flow, not u_reservoir.n.
    let (volume, _) = solved
        .reactions
        .volume_flux_moments(&[left], origin)
        .unwrap();
    near(volume, -0.2, 3e-13);
}

#[test]
fn reservoir_momentum_drives_bdf_evolution_with_independent_force_and_torque_balance() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let history = moments(&mesh, &|_| [0.08, -0.02]);
    let snapshot = history.clone();
    let force = |p: [f64; 2]| [0.2 + 0.1 * p[1], -0.05 * p[0]];
    let traction = |p: [f64; 2], n: [f64; 2]| n.map(|x| -(4.0 - 0.3 * p[0]) * x);
    let dt = 0.04;
    let origin = [0.2, -0.15];
    let mut endpoints = Vec::new();
    for transverse in [0.0, 0.4] {
        let reservoir = |p: [f64; 2]| [0.1, transverse + 0.1 * p[1]];
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| NavierStokesBoundary::OpenReservoir {
                edge,
                traction: &traction,
                velocity: &reservoir,
            })
            .collect();
        let step = system
            .bdf1_navier_stokes_mixed_step(
                FluxParams {
                    nu: 0.3,
                    tol: 1e-12,
                    ..FluxParams::default()
                },
                &force,
                &boundaries,
                &history,
                dt,
                NavierStokesControls {
                    relative_tolerance: 1e-11,
                    ..NavierStokesControls::default()
                },
            )
            .unwrap();
        assert!(step.iterations > 1);
        let all: Vec<_> = step
            .solution
            .reactions
            .edges()
            .iter()
            .map(|r| r.edge)
            .collect();
        let load = step
            .solution
            .reactions
            .resultant(&all, origin, &[None])
            .unwrap();
        let mut body_minus_inertia = [0.0; 2];
        let mut moment_minus_inertia = 0.0;
        for (cell, vertices) in mesh.tris.iter().enumerate() {
            for (point, weight) in tri_quad(vertices.map(|v| mesh.verts[v]), mesh.areas[cell]) {
                let u = velocity_in_cell(&system, &step.solution.field.x, cell, point);
                let f = force(point);
                let balance = [f[0] - (u[0] - 0.08) / dt, f[1] - (u[1] + 0.02) / dt];
                for axis in 0..2 {
                    body_minus_inertia[axis] += weight * balance[axis];
                }
                moment_minus_inertia +=
                    weight * cross([point[0] - origin[0], point[1] - origin[1]], balance);
            }
        }
        let mut outward = [0.0; 2];
        let mut torque = 0.0;
        for &edge in &all {
            let cell = mesh.edges[edge].tris.0;
            let u = |point| velocity_in_cell(&system, &step.solution.field.x, cell, point);
            let (momentum, moment) = flux(&mesh, edge, &u, &u, &reservoir, origin);
            for axis in 0..2 {
                outward[axis] += momentum[axis];
            }
            torque += moment;
        }
        for axis in 0..2 {
            near(load.outward_momentum.force[axis], outward[axis], 2e-9);
            near(
                load.stress.force[axis] + outward[axis],
                body_minus_inertia[axis],
                2e-8,
            );
        }
        near(load.outward_momentum.moment, torque, 2e-9);
        near(load.stress.moment + torque, moment_minus_inertia, 2e-8);
        assert_eq!(load.backflow.force, [0.0; 2]);
        endpoints.push(step.solution.field.x);
    }
    let changed = endpoints[0][..system.n_u]
        .iter()
        .zip(&endpoints[1])
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        changed > 1e-5,
        "reservoir momentum must reach the actual evolving fluid"
    );
    assert_eq!(history, snapshot);
}

#[test]
fn reservoir_admission_and_cached_samples_preserve_history_on_refusal() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let history = moments(&mesh, &|_| [0.08, 0.0]);
    let snapshot = history.clone();
    let queries = Cell::new(0);
    let force_queries = Cell::new(0);
    let reservoir = |_: [f64; 2]| {
        queries.set(queries.get() + 1);
        [0.1, 0.2]
    };
    let traction = |_: [f64; 2], n: [f64; 2]| n.map(|v| -4.0 * v);
    let force = |_: [f64; 2]| {
        force_queries.set(force_queries.get() + 1);
        [0.0; 2]
    };
    let edges: Vec<_> = (0..mesh.edges.len())
        .filter(|&edge| mesh.is_boundary(edge))
        .collect();
    let boundaries: Vec<_> = edges
        .iter()
        .map(|&edge| NavierStokesBoundary::OpenReservoir {
            edge,
            traction: &traction,
            velocity: &reservoir,
        })
        .collect();
    let params = FluxParams {
        nu: 0.3,
        ..FluxParams::default()
    };
    let controls = NavierStokesControls {
        relative_tolerance: 1e-10,
        ..NavierStokesControls::default()
    };
    let mut duplicate = boundaries.clone();
    duplicate.push(boundaries[0]);
    assert!(
        system
            .bdf1_navier_stokes_mixed_step(params, &force, &duplicate, &history, 0.02, controls,)
            .is_err()
    );
    assert_eq!(queries.get(), 0);
    assert_eq!(force_queries.get(), 0);
    let invalid = |_: [f64; 2]| [f64::NAN, 0.0];
    let bad: Vec<_> = edges
        .iter()
        .map(|&edge| NavierStokesBoundary::OpenReservoir {
            edge,
            traction: &traction,
            velocity: &invalid,
        })
        .collect();
    assert!(matches!(
        system.bdf1_navier_stokes_mixed_step(params, &force, &bad, &history, 0.02, controls,),
        Err(StokesBoundaryError::NonFiniteBoundary { .. })
    ));
    assert_eq!(force_queries.get(), 0);
    let solved = system
        .bdf1_navier_stokes_mixed_step(params, &force, &boundaries, &history, 0.02, controls)
        .unwrap();
    assert_eq!(queries.get(), 2 * edges.len());
    for reaction in solved.solution.reactions.edges() {
        reaction.load([0.2, 0.3], None).unwrap();
    }
    assert_eq!(
        queries.get(),
        2 * edges.len(),
        "reporting must not re-query reservoir data"
    );
    assert!(matches!(
        system.bdf1_navier_stokes_mixed_step(
            params,
            &force,
            &boundaries,
            &history,
            0.02,
            NavierStokesControls {
                maximum_iterations: 1,
                ..controls
            },
        ),
        Err(StokesBoundaryError::Rheology(
            FluxRheologyError::NonlinearConvergence { .. }
        ))
    ));
    assert_eq!(history, snapshot);
}
