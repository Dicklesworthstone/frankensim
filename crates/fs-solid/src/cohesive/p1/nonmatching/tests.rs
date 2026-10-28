use super::*;
use crate::cohesive::{Budget, CohesiveSystem, Config};
use fs_solver::{LinearOp, NewtonKrylovConfig};

fn mesh() -> (
    Vec<[f64; 3]>,
    Vec<[usize; 4]>,
    [[usize; 3]; 2],
    [[usize; 3]; 2],
) {
    // A and B divide the unit square along opposite diagonals. Each face is
    // exterior to a distinct tetrahedron, with shared vertices within a side.
    let nodes = vec![
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [0.5, 0.5, -1.],
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [0.5, 0.5, 1.],
    ];
    let a = [[0, 1, 2], [0, 2, 3]];
    let b = [[5, 6, 8], [6, 7, 8]];
    (
        nodes,
        vec![[0, 1, 2, 4], [0, 2, 3, 4], [5, 6, 8, 9], [6, 7, 8, 9]],
        a,
        b,
    )
}
fn budget() -> NonmatchingBudget {
    NonmatchingBudget {
        mesh: TraceBudget {
            maximum_nodes: 100,
            maximum_tetrahedra: 100,
            maximum_interfaces: 10,
            minimum_face_quality: 1e-10,
        },
        maximum_surface_faces: 40,
        maximum_pair_checks: 100,
        maximum_quadrature_points: 100,
        relative_geometry_tolerance: 1e-11,
    }
}
fn law() -> BilinearCohesive {
    BilinearCohesive::new(1e9, 2e5, 100., 2e9, 0.01).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn build<'a>(
    nodes: &[[f64; 3]],
    tets: &[[usize; 4]],
    free: &[usize],
    a: &'a [[usize; 3]],
    b: &'a [[usize; 3]],
    law: &'a BilinearCohesive,
) -> Result<P1Interfaces<'a>, Error> {
    P1Interfaces::from_nonmatching_controlled(
        nodes,
        tets,
        free,
        &[NonmatchingSurface {
            id: 71,
            side_a: a,
            side_b: b,
            law,
        }],
        budget(),
        &mut || false,
    )
}

#[test]
fn crossed_diagonals_integrate_area_and_linear_bond_energy_exactly() {
    let (nodes, tets, a, b) = mesh();
    let law = law();
    let faces = build(&nodes, &tets, &(0..30).collect::<Vec<_>>(), &a, &b, &law).unwrap();
    assert_eq!(faces.points().len(), 12);
    close(
        faces.points().iter().map(|p| p.area_weight_m2).sum(),
        1.,
        3e-16,
    );
    let mut u = vec![0.; 30];
    for i in 5..9 {
        u[3 * i + 2] = 1e-5 * (1. + nodes[i][0] + nodes[i][1]);
    }
    let jumps = faces.jumps(&u, &mut || false).unwrap();
    let integral: f64 = jumps
        .iter()
        .zip(faces.points())
        .map(|(jump, point)| point.area_weight_m2 * jump[0] * jump[0])
        .sum();
    // Integral_0^1 Integral_0^1 (1+x+y)^2 dx dy = 25/6.
    close(integral, 1e-10 * 25. / 6., 1e-24);
    let traction: Vec<_> = jumps.iter().map(|j| j.map(|x| 1e9 * x)).collect();
    let forces = faces.full_forces(&traction, &mut || false).unwrap();
    close(
        0.5 * forces.iter().zip(&u).map(|(f, u)| f * u).sum::<f64>(),
        0.5e9 * integral,
        2e-16,
    );
    assert!(
        faces
            .points()
            .iter()
            .any(|p| p.shape_values != p.side_b_shape_values)
    );
}

#[test]
fn overlap_map_preserves_virtual_work_resultants_moments_and_rigid_motion() {
    let (nodes, tets, a, b) = mesh();
    let law = law();
    let faces = build(&nodes, &tets, &(0..30).collect::<Vec<_>>(), &a, &b, &law).unwrap();
    let u: Vec<_> = (0..30).map(|i| (i as f64 * 0.7).sin() * 0.01).collect();
    let tractions: Vec<_> = (0..faces.points().len())
        .map(|i| [17. + i as f64, -31., 12.])
        .collect();
    let jumps = faces.jumps(&u, &mut || false).unwrap();
    let force = faces.full_forces(&tractions, &mut || false).unwrap();
    close(
        force.iter().zip(&u).map(|(f, u)| f * u).sum(),
        jumps
            .iter()
            .zip(&tractions)
            .zip(faces.points())
            .map(|((j, t), p)| p.area_weight_m2 * dot(*j, *t))
            .sum(),
        1e-15,
    );
    let mut resultant = [0.; 3];
    let mut moment = [0.; 3];
    let mut rigid = Vec::new();
    for (node, f) in nodes.iter().zip(force.chunks_exact(3)) {
        let torque = super::super::cross(*node, [f[0], f[1], f[2]]);
        for c in 0..3 {
            resultant[c] += f[c];
            moment[c] += torque[c];
        }
        rigid.extend(super::super::cross([0.003, 0.001, -0.002], *node).map(|v| v + 0.01));
    }
    for v in resultant.into_iter().chain(moment) {
        close(v, 0., 1e-14);
    }
    for jump in faces.jumps(&rigid, &mut || false).unwrap() {
        for v in jump {
            close(v, 0., 1e-17);
        }
    }
}

#[test]
fn identical_triangulations_recover_matching_consistent_stiffness() {
    let (nodes, _, a, _) = mesh();
    let law = law();
    let b = [[5, 6, 7], [5, 7, 8]];
    let tets = [[0, 1, 2, 4], [0, 2, 3, 4], [5, 6, 7, 9], [5, 7, 8, 9]];
    let free: Vec<_> = (0..30).collect();
    let overlap = build(&nodes, &tets, &free, &a, &b, &law).unwrap();
    let pairs: Vec<_> = (0..2)
        .map(|i| super::super::FacePair {
            id: i as u64,
            side_a: a[i],
            side_b: b[i],
            law: &law,
        })
        .collect();
    let matching =
        P1Interfaces::new_controlled(&nodes, &tets, &free, &pairs, budget().mesh, &mut || false)
            .unwrap();
    for dof in 0..30 {
        let mut u = vec![0.; 30];
        u[dof] = 1.;
        let overlap_jumps = overlap.jumps(&u, &mut || false).unwrap();
        let matching_jumps = matching.jumps(&u, &mut || false).unwrap();
        let f1 = overlap.full_forces(&overlap_jumps, &mut || false).unwrap();
        let f2 = matching
            .full_forces(&matching_jumps, &mut || false)
            .unwrap();
        for (x, y) in f1.into_iter().zip(f2) {
            close(x, y, 1e-16);
        }
    }
}

#[test]
fn nonmatching_force_controlled_pull_softens_fractures_and_recloses() {
    // B's independent z DOFs couple to a positive bulk stiffness. The exact
    // uniform solution has F/A = k_bulk*u + t(u), including the softening rung.
    struct Bulk([f64; 4]);
    impl LinearOp for Bulk {
        fn n(&self) -> usize {
            4
        }
        fn apply(&self, x: &[f64], y: &mut [f64]) {
            for i in 0..4 {
                y[i] = 1e9 * self.0[i] * x[i];
            }
        }
    }
    let (nodes, tets, a, b) = mesh();
    let law = law();
    let faces = build(&nodes, &tets, &[17, 20, 23, 26], &a, &b, &law).unwrap();
    let areas = [1. / 6., 1. / 3., 1. / 6., 1. / 3.];
    let bulk = Bulk(areas);
    let system = CohesiveSystem::new(
        &bulk,
        faces.patches().to_vec(),
        Budget {
            max_dofs: 4,
            max_patches: 20,
            max_jump_terms: 36,
            max_krylov_basis_entries: 4096,
        },
    )
    .unwrap();
    let config = Config {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-13,
            linear_restart: 4,
            forcing_maximum: 0.05,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 80,
        max_operator_evaluations: 2000,
        force_tolerance_n: 1e-7,
        local_energy_tolerance_j_m2: 1e-7,
        total_cohesive_energy_tolerance_j: 1e-7,
        maximum_displacement_m: 0.01,
    };
    let initial = system.bonded_state();
    let mut state = initial.clone();
    let offsets = vec![[0.; 3]; faces.points().len()];
    let mut dissipation = 0.;
    for (displacement, traction) in [(1e-4, 1e5), (3e-4, 1.75e5), (7e-4, 7.5e4), (1.2e-3, 0.)] {
        let force = areas.map(|a| a * (1e9 * displacement + traction));
        let result = state
            .solve_controlled(&force, &offsets, config, &mut || false)
            .unwrap();
        for &u in result.state.displacements_m() {
            close(u, displacement, 2e-15);
        }
        dissipation += result.cohesive_energy_j.dissipated;
        state = result.state;
    }
    assert!(initial.points().iter().all(|p| p.damage() == 0.));
    assert!(state.points().iter().all(|p| p.failed()));
    close(dissipation, 100., 1e-8);
    let force = areas.map(|a| -3e4 * a);
    let closed = state
        .solve_controlled(&force, &offsets, config, &mut || false)
        .unwrap();
    for &u in closed.state.displacements_m() {
        close(u, -1e-5, 2e-15);
    }
    for p in closed.state.points() {
        assert!(p.failed());
        close(p.traction_pa()[0], -2e4, 1e-6);
    }
}

#[test]
fn sweep_spends_pair_budget_only_on_spatial_candidates() {
    let law = law();
    let mut nodes = Vec::new();
    let mut tets = Vec::new();
    let mut a = Vec::new();
    let mut b = Vec::new();
    for i in 0..32 {
        let first = nodes.len();
        let x = 3. * i as f64;
        nodes.extend([
            [x, 0., 0.],
            [x + 1., 0., 0.],
            [x, 1., 0.],
            [x, 0., -1.],
            [x, 0., 0.],
            [x + 1., 0., 0.],
            [x, 1., 0.],
            [x, 0., 1.],
        ]);
        a.push([first, first + 1, first + 2]);
        b.push([first + 4, first + 5, first + 6]);
        tets.extend([
            [first, first + 1, first + 2, first + 3],
            [first + 4, first + 5, first + 6, first + 7],
        ]);
    }
    let mut cap = budget();
    cap.mesh.maximum_nodes = 256;
    cap.maximum_surface_faces = 64;
    cap.maximum_pair_checks = 32;
    cap.maximum_quadrature_points = 96;
    let faces = P1Interfaces::from_nonmatching_controlled(
        &nodes,
        &tets,
        &[],
        &[NonmatchingSurface {
            id: 9,
            side_a: &a,
            side_b: &b,
            law: &law,
        }],
        cap,
        &mut || false,
    )
    .unwrap();
    assert_eq!(faces.points().len(), 96);
    close(
        faces.points().iter().map(|p| p.area_weight_m2).sum(),
        16.,
        3e-14,
    );
}

#[test]
fn missing_coverage_same_side_overlap_nonplanarity_budgets_and_cancel_refuse() {
    let (nodes, tets, a, b) = mesh();
    let law = law();
    let free: Vec<_> = (0..30).collect();
    assert!(build(&nodes, &tets, &free, &a, &b[..1], &law).is_err());
    let mut shifted = nodes.clone();
    shifted[7][2] = 1e-4;
    assert!(build(&shifted, &tets, &free, &a, &b, &law).is_err());
    let mut same_side = nodes.clone();
    same_side[9][2] = -1.;
    assert!(build(&same_side, &tets, &free, &a, &b, &law).is_err());
    let mut overlap_nodes = nodes.clone();
    overlap_nodes.extend([nodes[0], nodes[1], nodes[2], [0.5, 0.5, -2.]]);
    let mut overlap_tets = tets.clone();
    overlap_tets.push([10, 11, 12, 13]);
    let overlap_faces = [a[0], a[1], [10, 11, 12]];
    assert!(
        build(
            &overlap_nodes,
            &overlap_tets,
            &free,
            &overlap_faces,
            &b,
            &law
        )
        .is_err()
    );
    let surfaces = [NonmatchingSurface {
        id: 7,
        side_a: &a,
        side_b: &b,
        law: &law,
    }];
    let mut tiny = budget();
    tiny.maximum_pair_checks = 1;
    assert!(matches!(
        P1Interfaces::from_nonmatching_controlled(
            &nodes,
            &tets,
            &free,
            &surfaces,
            tiny,
            &mut || false
        ),
        Err(Error::Budget(_))
    ));
    tiny = budget();
    tiny.maximum_quadrature_points = 11;
    assert!(matches!(
        P1Interfaces::from_nonmatching_controlled(
            &nodes,
            &tets,
            &free,
            &surfaces,
            tiny,
            &mut || false
        ),
        Err(Error::Budget(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        P1Interfaces::from_nonmatching_controlled(
            &nodes,
            &tets,
            &free,
            &surfaces,
            budget(),
            &mut || {
                polls += 1;
                polls == 12
            }
        ),
        Err(Error::Cancelled)
    ));
}
