use super::*;
use crate::shell::ShellMesh;

// Kinematic fixture only. Production ShellReduction construction still verifies
// the complete physical pencil/mass basis; fs-couple tests use that constructor.
fn source() -> ShellReduction {
    let mesh = ShellMesh::new(
        vec![[0., 0., 0.], [1., 0., 0.2], [1., 1., 0.4], [0., 1., -0.1]],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();
    let translations = (0..2)
        .flat_map(|mode| mesh.nodes.iter().map(move |&p| motion(mode, p)))
        .collect();
    let area_normals = (0..2)
        .map(|f| {
            let g = mesh.facet(f).unwrap();
            g.frame[2].map(|v| v * g.area_m2)
        })
        .collect();
    ShellReduction {
        omegas: vec![1., 2.],
        remainder: vec![1., 0., 0., 4.],
        facets: vec![],
        translations,
        rotations: vec![[12., -7., 31.]; 8],
        reference_positions: mesh.nodes,
        section_thicknesses: vec![0.01, 0.02],
        nodes: 4,
        triangles: mesh.tris,
        area_normals,
        rigid_motion: None,
    }
}
fn motion(mode: usize, p: [f64; 3]) -> [f64; 3] {
    if mode == 0 {
        [0.4, -0.3, 0.9]
    } else {
        [
            0.2 * p[0] + 0.1,
            -0.4 * p[1],
            0.6 * p[0] - 0.5 * p[1] + 0.2 * p[2],
        ]
    }
}
fn budget() -> OpenScreenRadiationBudget {
    OpenScreenRadiationBudget {
        maximum_vertices: 256,
        maximum_facets: 128,
        maximum_corner_weights: 768,
    }
}
fn pressure(p: [f64; 3]) -> f64 {
    0.7 * p[0] - 0.2 * p[1] + 0.3 * p[2] + 0.9
}
fn work(chart: &ShellOpenScreenRadiation, mode: usize) -> f64 {
    let mesh = ShellMesh::new(chart.vertices.clone(), chart.triangles.clone()).unwrap();
    chart
        .triangles
        .iter()
        .enumerate()
        .map(|(f, tri)| {
            let area = mesh.facet(f).unwrap().area_m2;
            let values = tri.map(|node| pressure(chart.vertices[node]));
            let b = chart.weights[mode][f];
            // Independent degree-two exact midpoint rule, not a lumped nodal mass.
            [(0, 1), (1, 2), (2, 0)]
                .iter()
                .map(|&(i, j)| {
                    (area / 3.) * (0.5 * b[i] + 0.5 * b[j]) * (0.5 * values[i] + 0.5 * values[j])
                })
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn g1_folded_midsurface_keeps_each_facet_trace_and_exact_p1_work_under_three_refinements() {
    let shell = source();
    let coarse = shell.open_screen_radiation(0, budget(), || false).unwrap();
    assert_ne!(
        coarse.weights[0][0][0], coarse.weights[0][1][0],
        "one shared vertex has two physical normal velocities"
    );
    assert!(
        coarse.weights[0].iter().flatten().all(|v| *v != 0.),
        "free-edge drive is physical, not a pressure boundary value"
    );
    for level in 0..=2 {
        let chart = shell
            .open_screen_radiation(level, budget(), || false)
            .unwrap();
        assert_eq!(chart.acoustic_refinements(), level);
        assert_eq!(chart.triangles.len(), 2 * 4usize.pow(level));
        assert_eq!(chart.weights.len(), 2);
        for (f, tri) in chart.triangles.iter().enumerate() {
            let parent = f / 4usize.pow(level);
            let normal = shell.area_normals[parent];
            let normal = normal.map(|v| v / dot(normal, normal).sqrt());
            for mode in 0..2 {
                for a in 0..3 {
                    let expected = dot(normal, motion(mode, chart.vertices[tri[a]]));
                    assert!((chart.weights[mode][f][a] - expected).abs() < 2e-15);
                }
            }
        }
        for mode in 0..2 {
            assert!((work(&chart, mode) - work(&coarse, mode)).abs() < 2e-14);
        }
        if level == 1 {
            assert_eq!(
                chart.vertices.len(),
                9,
                "shared diagonal midpoint must be welded by topology"
            );
        }
    }
}

#[test]
fn g0_complete_source_budget_and_cancellation_do_not_publish_partial_refinement() {
    let shell = source();
    let expected = shell.open_screen_radiation(2, budget(), || false).unwrap();
    let mut small = budget();
    small.maximum_corner_weights = 191;
    assert!(shell.open_screen_radiation(2, small, || false).is_err());
    let mut small = budget();
    small.maximum_vertices = 8;
    assert!(shell.open_screen_radiation(1, small, || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        shell.open_screen_radiation(2, budget(), || {
            calls += 1;
            calls > 8
        }),
        Err(PlateError::Cancelled)
    ));
    let replay = shell.open_screen_radiation(2, budget(), || false).unwrap();
    assert_eq!(expected.vertices, replay.vertices);
    assert_eq!(expected.triangles, replay.triangles);
    assert_eq!(expected.weights, replay.weights);
}
