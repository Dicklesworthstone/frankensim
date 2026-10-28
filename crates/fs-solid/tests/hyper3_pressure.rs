//! Closed-pressure forces checked against independent tetrahedron volume
//! derivatives, not the surface integration used by the production owner.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use fs_solid::linear3::TetAssemblyBudget;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 19,
                    kernel_id: 58,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        )
    })
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
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}

/// d det[x1-x0,x2-x0,x3-x0]/dx, derived directly from the scalar triple product.
fn tet_oracle(x: [[f64; 3]; 4], pressure: f64) -> (f64, [[f64; 3]; 4]) {
    let a = sub(x[1], x[0]);
    let b = sub(x[2], x[0]);
    let c = sub(x[3], x[0]);
    let volume = dot(a, cross(b, c)) / 6.0;
    let mut force = [[0.0; 3]; 4];
    for (node, gradient) in [cross(b, c), cross(c, a), cross(a, b)]
        .into_iter()
        .enumerate()
    {
        for i in 0..3 {
            force[node + 1][i] = -pressure * gradient[i] / 6.0;
            force[0][i] -= force[node + 1][i];
        }
    }
    (volume, force)
}

#[test]
fn g0_current_area_force_tangent_and_exact_work_match_volume_derivatives() {
    with_cx(|_, cx| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        assert_eq!(surface.faces().len(), 4);
        assert!(surface.matches_reference(&NODES, &TETS));
        let x = [
            [0.1, -0.2, 0.3],
            [1.4, -0.1, 0.35],
            [0.25, 0.6, 0.25],
            [0.15, -0.05, 1.4],
        ];
        let u: Vec<f64> = x
            .iter()
            .zip(NODES)
            .flat_map(|(x, reference)| (0..3).map(move |i| x[i] - reference[i]))
            .collect();
        let pressure = 3.7;
        let (volume, expected) = tet_oracle(x, pressure);
        near(surface.volume_m3(&u, cx).unwrap(), volume, 2.0e-15);
        let force = surface.force_n(&u, pressure, cx).unwrap();
        for (actual, expected) in force.iter().zip(expected.iter().flatten()) {
            near(*actual, *expected, 2.0e-14);
        }
        let tangent = surface.force_tangent(&u, pressure, cx).unwrap();
        let epsilon = 1.0e-6;
        for column in 0..12 {
            let mut plus = u.clone();
            plus[column] += epsilon;
            let mut minus = u.clone();
            minus[column] -= epsilon;
            let forward = surface.force_n(&plus, pressure, cx).unwrap();
            let backward = surface.force_n(&minus, pressure, cx).unwrap();
            for row in 0..12 {
                near(
                    tangent.get(row, column),
                    (forward[row] - backward[row]) / (2.0 * epsilon),
                    5.0e-9,
                );
                near(tangent.get(row, column), tangent.get(column, row), 2.0e-14);
            }
        }
        let before = vec![0.0; 12];
        let work = surface.work_j(&before, &u, pressure, cx).unwrap();
        near(work, -pressure * (volume - 1.0 / 6.0), 2.0e-14);
        near(
            surface.work_j(&u, &before, pressure, cx).unwrap(),
            -work,
            2.0e-14,
        );
        near(
            surface.work_j(&before, &u, -pressure, cx).unwrap(),
            -work,
            2.0e-14,
        );
        assert_eq!(surface.force_n(&u, 0.0, cx).unwrap(), vec![0.0; 12]);
    });
}

#[test]
fn g3_pressure_rotates_with_the_body_and_has_zero_resultant_and_moment() {
    with_cx(|_, cx| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let pressure = 2.0;
        let stretch = [1.3, 0.9, 1.1];
        let x: [[f64; 3]; 4] = NODES.map(|point| std::array::from_fn(|i| stretch[i] * point[i]));
        let rotate = |point: [f64; 3]| [-point[1], point[0], point[2]];
        let shift = [8.0, -4.0, 2.0];
        let transformed = x.map(|point| {
            let r = rotate(point);
            std::array::from_fn::<_, 3, _>(|i| r[i] + shift[i])
        });
        let displacement = |points: [[f64; 3]; 4]| -> Vec<f64> {
            points
                .into_iter()
                .zip(NODES)
                .flat_map(|(x, reference)| sub(x, reference))
                .collect()
        };
        let u = displacement(x);
        let moved_u = displacement(transformed);
        let force = surface.force_n(&u, pressure, cx).unwrap();
        let moved_force = surface.force_n(&moved_u, pressure, cx).unwrap();
        near(
            surface.volume_m3(&u, cx).unwrap(),
            surface.volume_m3(&moved_u, cx).unwrap(),
            1.0e-14,
        );
        let mut resultant = [0.0; 3];
        let mut moment = [0.0; 3];
        for node in 0..4 {
            let original = [force[3 * node], force[3 * node + 1], force[3 * node + 2]];
            let transformed_force = [
                moved_force[3 * node],
                moved_force[3 * node + 1],
                moved_force[3 * node + 2],
            ];
            let expected = rotate(original);
            let torque = cross(transformed[node], transformed_force);
            for i in 0..3 {
                near(transformed_force[i], expected[i], 2.0e-14);
                resultant[i] += transformed_force[i];
                moment[i] += torque[i];
            }
        }
        for value in resultant.into_iter().chain(moment) {
            near(value, 0.0, 5.0e-14);
        }
    });
}

#[test]
fn g1_shared_interior_face_cancels_and_boundary_loads_match_volume_union() {
    with_cx(|_, cx| {
        let nodes = [NODES[0], NODES[1], NODES[2], NODES[3], [0.0, 0.0, -1.0]];
        let tets = [[0, 1, 2, 3], [0, 2, 1, 4]];
        let surface =
            ClosedTetPressureSurface::new(&nodes, &tets, TetAssemblyBudget::standard(), 8, cx)
                .unwrap();
        assert_eq!(surface.faces().len(), 6);
        for face in surface.faces() {
            let mut sorted = *face;
            sorted.sort_unstable();
            assert_ne!(sorted, [0, 1, 2], "interior face must not receive pressure");
        }
        let u = vec![0.0; 15];
        near(surface.volume_m3(&u, cx).unwrap(), 1.0 / 3.0, 1.0e-15);
        let pressure = 1.2;
        let mut expected = [0.0; 15];
        for tet in tets {
            let (_, local) = tet_oracle(tet.map(|node| nodes[node]), pressure);
            for (a, node) in tet.into_iter().enumerate() {
                for i in 0..3 {
                    expected[3 * node + i] += local[a][i];
                }
            }
        }
        for (actual, expected) in surface
            .force_n(&u, pressure, cx)
            .unwrap()
            .into_iter()
            .zip(expected)
        {
            near(actual, expected, 1.0e-14);
        }
        assert!(
            ClosedTetPressureSurface::new(
                &nodes,
                &[tets[0], tets[0]],
                TetAssemblyBudget::standard(),
                8,
                cx,
            )
            .is_err(),
            "duplicate cells do not make an admissible closed solid"
        );
    });
}

#[test]
fn g4_invalid_geometry_budget_and_cancel_refuse_before_a_load_is_returned() {
    with_cx(|gate, cx| {
        assert!(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 3, cx,)
                .is_err()
        );
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let mut inverted = vec![0.0; 12];
        inverted[3] = -2.0;
        assert!(surface.force_n(&inverted, 1.0, cx).is_err());
        assert!(surface.volume_m3(&[0.0; 11], cx).is_err());
        assert!(surface.force_n(&[0.0; 12], f64::NAN, cx).is_err());
        gate.request();
        assert!(surface.force_n(&[0.0; 12], 1.0, cx).is_err());
        assert!(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx,)
                .is_err()
        );
    });
}

#[test]
fn g3_pressure_geometry_is_invariant_to_large_common_displacement() {
    with_cx(|_, cx| {
        // The reference edges are smaller than one ULP of the translation.
        // They remain known exactly through X and the common nodal motion.
        let scale = 2.0_f64.powi(-20);
        let nodes = NODES.map(|point| point.map(|value| scale * value));
        let surface =
            ClosedTetPressureSurface::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let zero = [0.0; 12];
        let translation = [2.0_f64.powi(40), -2.0_f64.powi(39), 2.0_f64.powi(38)];
        let translated: Vec<f64> = (0..4).flat_map(|_| translation).collect();
        let pressure = 3.7;
        let volume = surface.volume_m3(&zero, cx).unwrap();
        assert!(volume > 0.0);
        assert_eq!(surface.volume_m3(&translated, cx).unwrap(), volume);
        assert_eq!(
            surface.force_n(&translated, pressure, cx).unwrap(),
            surface.force_n(&zero, pressure, cx).unwrap(),
        );
        let reference_tangent = surface.force_tangent(&zero, pressure, cx).unwrap();
        let translated_tangent = surface.force_tangent(&translated, pressure, cx).unwrap();
        for row in 0..12 {
            assert_eq!(reference_tangent.row(row), translated_tangent.row(row));
        }
        assert_eq!(surface.work_j(&zero, &translated, pressure, cx).unwrap(), 0.0);
    });
}
