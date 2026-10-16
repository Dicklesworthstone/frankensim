//! G0/G1 outline contracts; all geometries are declared analytic fixtures.
use fs_exec::{CancelGate, Cx};
use fs_mesh::planar::*;
use std::collections::{BTreeMap, BTreeSet};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 55,
                kernel_id: 55,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn options(h: f64, angle: f64) -> PlanarOptions<'static> {
    PlanarOptions {
        maximum_edge_m: h,
        minimum_angle_deg: angle,
        size_field: None,
        budget: PlanarBudget {
            max_vertices: 20000,
            max_triangles: 40000,
            max_refinement_steps: 20000,
        },
    }
}
struct Ellipse {
    a: f64,
    b: f64,
}
fn point_in_loop(p: [f64; 2], ring: &[[f64; 2]]) -> bool {
    let mut crossings = 0;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if (a[1] > p[1]) != (b[1] > p[1])
            && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}
impl ParametricCurve for Ellipse {
    fn position(&self, t: f64) -> [f64; 2] {
        if t == 1.0 {
            return [self.a, 0.0];
        }
        let angle = 2.0 * core::f64::consts::PI * t;
        [
            self.a * fs_math::det::cos(angle),
            self.b * fs_math::det::sin(angle),
        ]
    }
    fn second_derivative_bound_m(&self) -> f64 {
        4.0 * core::f64::consts::PI * core::f64::consts::PI * self.a.max(self.b)
    }
}
fn audit(mesh: &PlanarMesh, holes: usize) {
    let mut edges = BTreeSet::new();
    let mut incidence = BTreeMap::new();
    for t in &mesh.triangles {
        assert_eq!(
            fs_ivl::orient2d(
                mesh.vertices[t[0]],
                mesh.vertices[t[1]],
                mesh.vertices[t[2]]
            ),
            fs_ivl::Sign::Positive
        );
        for k in 0..3 {
            let a = t[k];
            let b = t[(k + 1) % 3];
            let edge = (a.min(b), a.max(b));
            edges.insert(edge);
            *incidence.entry(edge).or_insert(0) += 1;
        }
    }
    for chain in &mesh.receipt.source_segment_chains {
        for pair in chain.windows(2) {
            let edge = (pair[0].min(pair[1]), pair[0].max(pair[1]));
            assert!(edges.contains(&edge));
            assert_eq!(incidence[&edge], 1);
        }
        let a = mesh.vertices[chain[0]];
        let b = mesh.vertices[*chain.last().unwrap()];
        for &i in chain {
            assert_eq!(fs_ivl::orient2d(a, b, mesh.vertices[i]), fs_ivl::Sign::Zero);
        }
    }
    assert_eq!(
        mesh.vertices.len() as isize - edges.len() as isize + mesh.triangles.len() as isize,
        1 - holes as isize
    );
    assert!((mesh.receipt.mesh_area_m2 / mesh.receipt.polygon_area_m2 - 1.0).abs() < 1e-12);
}
#[test]
fn g0_cdt_conforms_to_outline_and_holes() {
    with_cx(|cx| {
        let outer = [[-2.0, -2.0], [2.0, -2.0], [2.0, 2.0], [-2.0, 2.0]];
        let circle = Ellipse { a: 0.5, b: 0.5 };
        let sound = mesh_outline(
            &PlanarOutline {
                outer: Loop::Polyline(&outer),
                holes: vec![Loop::Parametric {
                    curve: &circle,
                    chord_tolerance_m: 0.01,
                }],
            },
            &options(10.0, 0.0),
            cx,
        )
        .unwrap();
        audit(&sound, 1);
        // No triangle centroid lies inside the hole's inscribed polygon:
        // every radial centroid must lie outside its apothem.
        let apothem = 0.5
            * fs_math::det::cos(
                core::f64::consts::PI / sound.receipt.boundary_loops[1].len() as f64,
            );
        for t in &sound.triangles {
            let x = t.iter().map(|&i| sound.vertices[i][0] / 3.0).sum::<f64>();
            let y = t.iter().map(|&i| sound.vertices[i][1] / 3.0).sum::<f64>();
            assert!(x * x + y * y >= apothem * apothem * (1.0 - 1e-12));
        }
        let left = [
            [-1.0, -1.0],
            [-0.5, -0.5],
            [-0.75, 0.0],
            [-0.5, 1.0],
            [-1.0, 0.5],
            [-0.875, 0.0],
        ];
        let right: Vec<_> = left.iter().map(|p| [-p[0], p[1]]).collect();
        let slots = mesh_outline(
            &PlanarOutline {
                outer: Loop::Polyline(&outer),
                holes: vec![Loop::Polyline(&left), Loop::Polyline(&right)],
            },
            &options(0.5, 20.0),
            cx,
        )
        .unwrap();
        audit(&slots, 2);
        for t in &slots.triangles {
            let p = core::array::from_fn(|k| {
                t.iter().map(|&i| slots.vertices[i][k] / 3.0).sum::<f64>()
            });
            assert!(!point_in_loop(p, &left) && !point_in_loop(p, &right));
        }
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W6\",\"test\":\"g0_cdt_conforms_to_outline_and_holes\",\"tier\":\"G0\",\"inputs\":{{\"outline\":\"fixture:two-nonconvex-slots\",\"holes\":2}},\"vertices\":{},\"triangles\":{},\"area_relative_error\":{},\"verdict\":true}}",
            slots.vertices.len(),
            slots.triangles.len(),
            (slots.receipt.mesh_area_m2 / slots.receipt.polygon_area_m2 - 1.0).abs()
        );
    });
}
#[test]
fn g0_quality_bound() {
    with_cx(|cx| {
        let outer = [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]];
        let hole = [[0.75, 0.25], [1.25, 0.25], [1.25, 0.75], [0.75, 0.75]];
        let mesh = mesh_outline(
            &PlanarOutline {
                outer: Loop::Polyline(&outer),
                holes: vec![Loop::Polyline(&hole)],
            },
            &options(0.25, 20.0),
            cx,
        )
        .unwrap();
        audit(&mesh, 1);
        for (i, t) in mesh.triangles.iter().enumerate() {
            if mesh.receipt.exempt_triangles.contains(&i) {
                assert!(
                    t.iter()
                        .any(|v| mesh.receipt.acute_input_vertices.contains(v))
                );
            } else {
                assert!(minimum_angle(t.map(|i| mesh.vertices[i])) >= 20.0 * (1.0 - 1e-12));
            }
        }
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W6\",\"test\":\"g0_quality_bound\",\"tier\":\"G0\",\"min_angle_deg\":{},\"exempt\":{:?},\"tolerance\":20.0,\"verdict\":true}}",
            mesh.receipt.minimum_angle_deg, mesh.receipt.exempt_triangles
        );
    });
}
#[test]
fn g1_curved_outline_chord_ladder() {
    with_cx(|cx| {
        for curve in [Ellipse { a: 1.0, b: 1.0 }, Ellipse { a: 1.0, b: 0.5 }] {
            let mut hs = Vec::new();
            let mut errors = Vec::new();
            for n in [16, 32, 64] {
                let eps = curve.second_derivative_bound_m() / (8.0 * (n * n) as f64);
                let mesh = mesh_outline(
                    &PlanarOutline {
                        outer: Loop::Parametric {
                            curve: &curve,
                            chord_tolerance_m: eps,
                        },
                        holes: vec![],
                    },
                    &options(10.0, 0.0),
                    cx,
                )
                .unwrap();
                audit(&mesh, 0);
                let exact = core::f64::consts::PI * curve.a * curve.b;
                let error = (mesh.receipt.mesh_area_m2 / exact - 1.0).abs();
                hs.push(1.0 / mesh.receipt.boundary_loops[0].len() as f64);
                errors.push(error);
                println!(
                    "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W6\",\"test\":\"g1_curved_outline_chord_ladder\",\"tier\":\"G1\",\"inputs\":{{\"ellipse_axes_m\":[{},{}],\"chord_eps_m\":{eps}}},\"vertices\":{},\"area_error\":{error},\"provenance\":\"analytic ellipse area\"}}",
                    curve.a,
                    curve.b,
                    mesh.vertices.len()
                );
            }
            let ladder = fs_mms::RefinementLadder::new(hs, errors).unwrap();
            let verdict = fs_mms::OrderGate { theoretical: 2.0 }
                .check("outline-chord-area", fs_mms::LadderSide::Primal, &ladder)
                .unwrap();
            println!("{}", verdict.json_line(true));
        }
    });
}
