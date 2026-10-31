//! Independent enclosed-volume identities and native solid pressure consumption.
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::pressure::cavity::ClosedTetCavity;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_solid::hyper3::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 59,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        )
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn pairing(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn rotate(a: [f64; 3]) -> [f64; 3] {
    [-a[1], a[0], a[2]]
}

/// A conforming hollow tetrahedron. Each radial triangular frustum is split
/// into three tets with the same diagonal rule on its shared quadrilaterals.
fn shell() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let nodes = vec![
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 3.0],
        [0.5, 0.5, 0.5],
        [1.5, 0.5, 0.5],
        [0.5, 1.5, 0.5],
        [0.5, 0.5, 1.5],
    ];
    let mut tets = Vec::new();
    for [a, b, c] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
        tets.push([a + 4, b + 4, c + 4, c]);
        tets.push([a + 4, b + 4, b, c]);
        tets.push([a + 4, a, b, c]);
    }
    (nodes, tets)
}

fn surface(nodes: &[[f64; 3]], tets: &[[usize; 4]], cx: &Cx<'_>) -> Arc<ClosedTetPressureSurface> {
    Arc::new(
        ClosedTetPressureSurface::new(
            nodes,
            tets,
            TetAssemblyBudget::standard(),
            4 * tets.len(),
            cx,
        )
        .unwrap(),
    )
}

fn inner_faces(surface: &ClosedTetPressureSurface) -> Vec<usize> {
    surface
        .faces()
        .iter()
        .enumerate()
        .filter_map(|(index, face)| face.iter().all(|&node| node >= 4).then_some(index))
        .collect()
}

fn cavity(nodes: &[[f64; 3]], tets: &[[usize; 4]], cx: &Cx<'_>) -> ClosedTetCavity {
    let surface = surface(nodes, tets, cx);
    let inner = inner_faces(&surface);
    assert_eq!(inner.len(), 4);
    assert_eq!(surface.faces().len(), 8);
    ClosedTetCavity::new(surface, &inner, cx).unwrap()
}

/// Scalar triple product and its cofactor derivatives. This oracle has no
/// triangle incidence, surface divergence formula or pressure-load kernel.
fn tetra_volume_gradient(x: [[f64; 3]; 4]) -> (f64, Vec<f64>) {
    let a = sub(x[1], x[0]);
    let b = sub(x[2], x[0]);
    let c = sub(x[3], x[0]);
    let mut g = vec![0.0; 12];
    for (i, cofactor) in [cross(b, c), cross(c, a), cross(a, b)]
        .into_iter()
        .enumerate()
    {
        for axis in 0..3 {
            g[3 * (i + 1) + axis] = cofactor[axis] / 6.0;
            g[axis] -= cofactor[axis] / 6.0;
        }
    }
    (dot(a, cross(b, c)) / 6.0, g)
}

fn affine_displacement(nodes: &[[f64; 3]]) -> Vec<f64> {
    let f = [[0.84, 0.06, -0.02], [0.02, 0.92, 0.04], [0.01, -0.03, 0.90]];
    let translation = [0.2, -0.3, 0.1];
    nodes
        .iter()
        .flat_map(|&x| (0..3).map(move |axis| dot(f[axis], x) + translation[axis] - x[axis]))
        .collect()
}

fn current_inner(nodes: &[[f64; 3]], u: &[f64]) -> [[f64; 3]; 4] {
    std::array::from_fn(|i| std::array::from_fn(|axis| nodes[4 + i][axis] + u[3 * (4 + i) + axis]))
}

#[test]
fn g0_gas_volume_derivatives_and_pressure_sign_match_an_independent_tetrahedron() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        assert_eq!(cavity.nodes(), &[4, 5, 6, 7]);
        near(cavity.reference_volume_m3(), 1.0 / 6.0, 2.0e-15);
        // Material and gas volumes are different even in the reference state.
        near(
            cavity.surface().volume_m3(&[0.0; 24], cx).unwrap(),
            13.0 / 3.0,
            2.0e-14,
        );

        let u = affine_displacement(&nodes);
        let state = cavity.evaluate(&u, cx).unwrap();
        assert_eq!(state.nodes(), cavity.nodes());
        let (expected_volume, expected_gradient) = tetra_volume_gradient(current_inner(&nodes, &u));
        near(state.volume_m3(), expected_volume, 2.0e-15);
        for (&actual, &expected) in state.gradient_m2().iter().zip(&expected_gradient) {
            near(actual, expected, 2.0e-14);
        }
        let pressure = 2.7;
        let face_pressures = cavity.face_pressures_pa(pressure, cx).unwrap();
        let forces = cavity
            .surface()
            .force_n(&u, TetPressure::PerFace(&face_pressures), cx)
            .unwrap();
        assert_eq!(&forces[..12], &[0.0; 12]);
        for i in 0..12 {
            near(forces[12 + i], pressure * expected_gradient[i], 3.0e-14);
        }
        let h = 1.0e-6;
        for (local, &node) in cavity.nodes().iter().enumerate() {
            for axis in 0..3 {
                let mut plus = u.clone();
                let mut minus = u.clone();
                plus[3 * node + axis] += h;
                minus[3 * node + axis] -= h;
                let derivative = (cavity.evaluate(&plus, cx).unwrap().volume_m3()
                    - cavity.evaluate(&minus, cx).unwrap().volume_m3())
                    / (2.0 * h);
                near(state.gradient_m2()[3 * local + axis], derivative, 2.0e-10);
            }
        }
        let work = cavity
            .held_pressure_work_j(&[0.0; 24], &u, pressure, cx)
            .unwrap();
        near(work, pressure * (expected_volume - 1.0 / 6.0), 2.0e-14);
        assert!(
            work < 0.0,
            "compression transfers positive work into the gas"
        );
        near(
            cavity
                .surface()
                .work_j(&[0.0; 24], &u, TetPressure::PerFace(&face_pressures), cx)
                .unwrap(),
            work,
            2.0e-14,
        );
        near(
            cavity
                .held_pressure_work_j(&u, &[0.0; 24], pressure, cx)
                .unwrap(),
            -work,
            2.0e-14,
        );
        near(
            cavity
                .held_pressure_work_j(&[0.0; 24], &u, -pressure, cx)
                .unwrap(),
            -work,
            2.0e-14,
        );
    });
}

#[test]
fn g0_matrix_free_volume_hessian_is_the_true_gradient_derivative_and_is_symmetric() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let u = affine_displacement(&nodes);
        let state = cavity.evaluate(&u, cx).unwrap();
        let d = [
            -0.2, 0.11, 0.07, 0.31, -0.08, 0.15, 0.05, 0.16, -0.18, -0.12, 0.04, 0.23,
        ];
        let e = [
            0.12, 0.04, -0.09, -0.13, 0.07, 0.14, -0.18, -0.19, 0.11, 0.05, -0.03, 0.17,
        ];
        let action = state.hessian_action(&d, cx).unwrap();
        let h = 1.0e-6;
        let mut plus = u.clone();
        let mut minus = u.clone();
        for (i, &node) in cavity.nodes().iter().enumerate() {
            for axis in 0..3 {
                plus[3 * node + axis] += h * d[3 * i + axis];
                minus[3 * node + axis] -= h * d[3 * i + axis];
            }
        }
        // The independent tetrahedron cofactors provide both perturbed gradients.
        let (_, gp) = tetra_volume_gradient(current_inner(&nodes, &plus));
        let (_, gm) = tetra_volume_gradient(current_inner(&nodes, &minus));
        for i in 0..d.len() {
            near(action[i], (gp[i] - gm[i]) / (2.0 * h), 1.0e-10);
        }
        let he = state.hessian_action(&e, cx).unwrap();
        near(pairing(&e, &action), pairing(&d, &he), 2.0e-14);
        let translation: Vec<f64> = (0..4).flat_map(|_| [0.2, -0.3, 0.4]).collect();
        assert_eq!(
            state.hessian_action(&translation, cx).unwrap(),
            vec![0.0; 12]
        );

        // A uniform closed cavity is conservative even though it is represented
        // as a sparse field on the full solid boundary.
        let pressures = cavity.face_pressures_pa(1.0, cx).unwrap();
        let jacobian = cavity
            .surface()
            .force_tangent(&u, TetPressure::PerFace(&pressures), cx)
            .unwrap();
        for row in 0..12 {
            let applied: f64 = (0..12)
                .map(|column| jacobian.get(12 + row, 12 + column) * d[column])
                .sum();
            near(applied, action[row], 2.0e-14);
        }
    });
}

#[test]
fn g3_cavity_volume_and_derivatives_are_objective_and_keep_large_rigid_translation() {
    with_cx(|_, cx| {
        let (nodes, mut tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let u = affine_displacement(&nodes);
        let state = cavity.evaluate(&u, cx).unwrap();
        let rotated: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(node, &x)| {
                let current = std::array::from_fn(|axis| x[axis] + u[3 * node + axis]);
                let qx = rotate(current);
                (0..3).map(move |axis| qx[axis] - x[axis])
            })
            .collect();
        let d: Vec<_> = (0..12).map(|i| 0.02 * (i as f64 - 4.0)).collect();
        let qd: Vec<_> = d
            .chunks_exact(3)
            .flat_map(|v| rotate([v[0], v[1], v[2]]))
            .collect();
        let before = state.hessian_action(&d, cx).unwrap();
        let after = cavity.evaluate(&rotated, cx).unwrap();
        let rotated_action = after.hessian_action(&qd, cx).unwrap();
        near(after.volume_m3(), state.volume_m3(), 2.0e-15);
        for node in 0..4 {
            let expected = rotate(std::array::from_fn(|axis| {
                state.gradient_m2()[3 * node + axis]
            }));
            let expected_action = rotate(std::array::from_fn(|axis| before[3 * node + axis]));
            for axis in 0..3 {
                near(
                    after.gradient_m2()[3 * node + axis],
                    expected[axis],
                    3.0e-14,
                );
                near(
                    rotated_action[3 * node + axis],
                    expected_action[axis],
                    3.0e-14,
                );
            }
        }
        let translation: Vec<_> = (0..nodes.len())
            .flat_map(|_| [1.0e12, -1.0e12, 1.0e12])
            .collect();
        near(
            cavity.evaluate(&translation, cx).unwrap().volume_m3(),
            1.0 / 6.0,
            2.0e-15,
        );

        // Either material-tet winding and either selection order give the same
        // cavity orientation and compact coordinate chart.
        for tet in &mut tets {
            tet.swap(1, 2);
        }
        let reversed_surface = surface(&nodes, &tets, cx);
        let mut selection = inner_faces(&reversed_surface);
        selection.reverse();
        let reversed = ClosedTetCavity::new(reversed_surface, &selection, cx).unwrap();
        let other = reversed.evaluate(&u, cx).unwrap();
        near(other.volume_m3(), state.volume_m3(), 2.0e-15);
        for (&a, &b) in other.gradient_m2().iter().zip(state.gradient_m2()) {
            near(a, b, 2.0e-14);
        }
    });
}

#[test]
fn g1_native_solid_solver_consumes_cavity_pressure_and_reaches_analytic_compression() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let (mu, lambda, stretch) = (2.0, 3.0, 0.92_f64);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let piola = mu * (stretch - 1.0 / stretch) + 3.0 * lambda * stretch.ln() / stretch;
        let pressure = -piola / stretch.powi(2);
        let pressures = cavity.face_pressures_pa(pressure, cx).unwrap();
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let (_, outer_gradient) = tetra_volume_gradient([nodes[0], nodes[1], nodes[2], nodes[3]]);
        // Reference dead forces on the OUTER shell equal its current pressure
        // traction at the analytic target. The inner shell uses real follower
        // pressure from the cavity mapping. No produced pressure force is used
        // to manufacture these independent outer forces.
        let dead: Vec<_> = (0..12)
            .filter(|dof| !fixed.iter().any(|&(id, _)| id == *dof))
            .map(|dof| (dof, piola * outer_gradient[dof]))
            .collect();
        let problem = HyperTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &fixed,
            nodal_forces_n: &dead,
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings {
                load_steps: 8,
                ..HyperTetSettings::default()
            },
        };
        let solution = problem
            .solve_with_face_pressures(cavity.surface(), &pressures, cx)
            .unwrap();
        let u: Vec<_> = solution.displacement_m.iter().flatten().copied().collect();
        for (&actual, &reference) in u.iter().zip(nodes.iter().flatten()) {
            near(actual, (stretch - 1.0) * reference, 3.0e-8);
        }
        let state = cavity.evaluate(&u, cx).unwrap();
        near(state.volume_m3(), stretch.powi(3) / 6.0, 2.0e-9);
        near(solution.min_det_f, stretch.powi(3), 3.0e-8);
        let log_j = 3.0 * stretch.ln();
        let density =
            0.5 * mu * (3.0 * stretch * stretch - 3.0) - mu * log_j + 0.5 * lambda * log_j * log_j;
        near(solution.strain_energy_j, (13.0 / 3.0) * density, 5.0e-8);
        for &(dof, reaction) in &solution.reactions_n {
            near(reaction, piola * outer_gradient[dof], 5.0e-8);
        }
        let work = cavity
            .held_pressure_work_j(&[0.0; 24], &u, pressure, cx)
            .unwrap();
        near(work, pressure * (stretch.powi(3) - 1.0) / 6.0, 5.0e-9);
        assert!(work < 0.0);
    });
}

#[test]
fn g4_cavity_refuses_bad_selections_invalid_geometry_and_cancellation() {
    with_cx(|gate, cx| {
        let (nodes, tets) = shell();
        let surface = surface(&nodes, &tets, cx);
        let inner = inner_faces(&surface);
        let outer: Vec<_> = (0..surface.faces().len())
            .filter(|i| !inner.contains(i))
            .collect();
        let all: Vec<_> = (0..surface.faces().len()).collect();
        for selection in [
            vec![],
            inner[..3].to_vec(),
            vec![inner[0], inner[0]],
            vec![surface.faces().len()],
            outer,
            all,
        ] {
            assert!(matches!(
                ClosedTetCavity::new(Arc::clone(&surface), &selection, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
        }
        let cavity = ClosedTetCavity::new(Arc::clone(&surface), &inner, cx).unwrap();
        for displacement in [
            vec![0.0; 23],
            vec![f64::NAN; 24],
            nodes.iter().flatten().map(|x| -x).collect(),
        ] {
            assert!(matches!(
                cavity.evaluate(&displacement, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
        }
        // An inverted element away from the selected cavity must also refuse.
        let mut bad = vec![0.0; 24];
        bad[3] = -4.0;
        assert!(cavity.evaluate(&bad, cx).is_err());
        let state = cavity.evaluate(&[0.0; 24], cx).unwrap();
        assert!(state.hessian_action(&[0.0; 11], cx).is_err());
        assert!(state.hessian_action(&[f64::INFINITY; 12], cx).is_err());
        for p in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(cavity.face_pressures_pa(p, cx).is_err());
            assert!(
                cavity
                    .held_pressure_work_j(&[0.0; 24], &[0.0; 24], p, cx)
                    .is_err()
            );
        }
        assert_eq!(
            cavity.face_pressures_pa(0.0, cx).unwrap(),
            vec![0.0; surface.faces().len()]
        );
        gate.request();
        assert!(matches!(
            ClosedTetCavity::new(Arc::clone(&surface), &inner, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            cavity.evaluate(&[0.0; 24], cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            cavity.face_pressures_pa(1.0, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            state.hessian_action(&[0.0; 12], cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            cavity.held_pressure_work_j(&[0.0; 24], &[0.0; 24], 1.0, cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}
