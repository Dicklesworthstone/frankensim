use super::*;

fn polynomial(x: f64, y: f64) -> f64 {
    0.1 + 0.2 * x - 0.3 * y + 0.4 * x * x * y - 0.1 * x * y * y * y + 0.2 * x * x * x * y * y
}
fn grid(f: impl Fn(f64, f64) -> f64) -> BicubicGrid {
    let x = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
    let y = x.clone();
    let values = y
        .iter()
        .flat_map(|&b| x.iter().map(|&a| f(a, b)).collect::<Vec<_>>())
        .collect();
    BicubicGrid::new(x, y, values).unwrap()
}
fn budget(n: usize, t: usize) -> ProfileBudget {
    ProfileBudget {
        max_nodes: n,
        max_triangles: t,
        max_feature_evaluations: 100 * n,
    }
}

#[test]
fn g0_bicubic_reproduction() {
    let g = grid(polynomial);
    for i in 0..=20 {
        for j in 0..=20 {
            let x = -1.0 + i as f64 / 10.0;
            let y = -1.0 + j as f64 / 10.0;
            let (z, [dx, dy]) = g.sample([x, y]).unwrap();
            assert!((z - polynomial(x, y)).abs() < 1e-12);
            assert!(
                (dx - (0.2 + 0.8 * x * y - 0.1 * y * y * y + 0.6 * x * x * y * y)).abs() < 1e-12
            );
            assert!(
                (dy - (-0.3 + 0.4 * x * x - 0.3 * x * y * y + 0.4 * x * x * x * y)).abs() < 1e-12
            );
        }
    }
}
#[test]
fn g0_section_interpolation_exact() {
    let ys = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
    let long = Curve2::new(
        ys.clone(),
        ys.iter().map(|y| 0.1 + 0.02 * y * y).collect(),
        ys.iter().map(|y| 0.04 * y).collect(),
    )
    .unwrap();
    let cross = ys
        .iter()
        .map(|&y| {
            (
                y,
                Curve2::new(
                    vec![-1.0, 0.0, 1.0],
                    vec![0.13 + 0.02 * y * y, 0.1 + 0.02 * y * y, 0.13 + 0.02 * y * y],
                    vec![-0.06, 0.0, 0.06],
                )
                .unwrap(),
            )
        })
        .collect();
    let surface = ArchingSurface::new(
        ArchingParameters::LoftedSections {
            long: long.clone(),
            cross,
        },
        ReferenceLoadState::Unloaded,
    )
    .unwrap();
    for &y in &ys {
        for i in 0..21 {
            let x = -1.0 + i as f64 / 10.0;
            assert!(
                (surface.sample([x, y]).unwrap().0 - (0.1 + 0.02 * y * y + 0.03 * x * x)).abs()
                    < 1e-12
            );
        }
    }
    for j in 0..21 {
        let y = -1.0 + j as f64 / 10.0;
        assert!((surface.sample([0.0, y]).unwrap().0 - long.sample(y).unwrap().0).abs() < 1e-12);
    }
}
#[test]
fn g0_missing_parameter_refuses() {
    assert!(Curve2::new(vec![0.0, 1.0], vec![0.0, 0.0], vec![]).is_err());
    assert!(Curve2::new(vec![0.0, 0.0], vec![0.0, 0.0], vec![0.0, 0.0]).is_err());
    assert!(BicubicGrid::new(vec![0.0, 1.0], vec![0.0, 1.0], vec![0.0; 4]).is_err());
    let long = Curve2::new(vec![-1.0, 1.0], vec![0.0, 0.0], vec![0.0, 0.0]).unwrap();
    assert!(
        ArchingSurface::new(
            ArchingParameters::LoftedSections {
                long,
                cross: vec![]
            },
            ReferenceLoadState::Unloaded
        )
        .is_err()
    );
    assert!(
        ArchingSurface::new(
            ArchingParameters::HeightMap(grid(|_, _| 0.0)),
            ReferenceLoadState::Loaded {
                description: String::new()
            }
        )
        .is_err()
    );
    assert!(grid(polynomial).sample([2.0, 0.0]).is_err());
    let mut values = vec![0.0; 16];
    values[5] = f64::NAN;
    assert!(BicubicGrid::new(vec![0.0, 1.0, 2.0, 3.0], vec![0.0, 1.0, 2.0, 3.0], values).is_err());
}
#[test]
fn g0_lifted_directors_unit_and_consistent() {
    let surface = ArchingSurface::new(
        ArchingParameters::HeightMap(grid(|x, y| 0.02 * (x * x + y * y))),
        ReferenceLoadState::Loaded {
            description: "declared 10 N fixture preload".into(),
        },
    )
    .unwrap();
    let mesh = PlateMesh::rectangle(1.0, 1.0, 4, 4);
    let lifted = surface
        .lift(&mesh, budget(mesh.nodes.len(), mesh.tris.len()), || false)
        .unwrap();
    assert_eq!(lifted.mesh.tris, mesh.tris);
    let work = 5 * mesh.nodes.len();
    let exact = ProfileBudget {
        max_nodes: mesh.nodes.len(),
        max_triangles: mesh.tris.len(),
        max_feature_evaluations: work,
    };
    assert!(surface.lift(&mesh, exact, || false).is_ok());
    assert!(
        surface
            .lift(
                &mesh,
                ProfileBudget {
                    max_feature_evaluations: work - 1,
                    ..exact
                },
                || false
            )
            .is_err()
    );
    assert_eq!(lifted.reference_state, surface.reference_state);
    for ((&(x, y), p), n) in mesh
        .nodes
        .iter()
        .zip(&lifted.mesh.nodes)
        .zip(&lifted.directors)
    {
        assert!((p[2] - 0.02 * (x * x + y * y)).abs() < 1e-12);
        assert!((n.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-14);
        assert!((n[0] + 0.04 * x * n[2]).abs() < 1e-14);
        assert!((n[1] + 0.04 * y * n[2]).abs() < 1e-14);
    }
    assert!(
        surface
            .lift(&mesh, budget(mesh.nodes.len() - 1, mesh.tris.len()), || {
                false
            })
            .is_err()
    );
    assert!(
        surface
            .lift(&mesh, budget(mesh.nodes.len(), mesh.tris.len() - 1), || {
                false
            })
            .is_err()
    );
    assert!(
        surface
            .lift(&mesh, budget(mesh.nodes.len(), mesh.tris.len()), || true)
            .is_err()
    );
    let thickness = ThicknessField::new(ArchingParameters::HeightMap(grid(|x, y| {
        0.003 + 0.0001 * (x + y)
    })))
    .unwrap();
    assert_eq!(
        thickness
            .sample_mesh(&mesh, mesh.nodes.len(), || false)
            .unwrap()
            .len(),
        mesh.nodes.len()
    );
    assert!(
        thickness
            .sample_mesh(&mesh, mesh.nodes.len() - 1, || false)
            .is_err()
    );
    assert!(
        ThicknessField::new(ArchingParameters::HeightMap(grid(|_, _| -0.01)))
            .unwrap()
            .sample([0.0, 0.0])
            .is_err()
    );
}

#[test]
fn g0_nurbs_seam_director_is_the_mean_of_one_sided_slopes() {
    use fs_rep_nurbs::KnotVector;
    // z = |x| + 0.25 y + 0.1 x y on [-1, 1]^2: a degree-one graph whose cell
    // line x = 0 is a C0 seam with one-sided slopes dz/dx = -1 + 0.1 y and
    // 1 + 0.1 y. The ordinary first partial is undefined there.
    let z = |x: f64, y: f64| x.abs() + 0.25 * y + 0.1 * x * y;
    let chart = |knots_x: Vec<f64>, degree: usize, xs: &[f64]| {
        let points: Vec<Vec<[f64; 3]>> = xs
            .iter()
            .map(|&x| [-1.0, 1.0].iter().map(|&y| [x, y, z(x, y)]).collect())
            .collect();
        Arc::new(
            NurbsSurface::new(
                KnotVector::new(knots_x, degree).unwrap(),
                KnotVector::new(vec![-1.0, -1.0, 1.0, 1.0], 1).unwrap(),
                &points,
                &vec![vec![1.0; 2]; xs.len()],
            )
            .unwrap(),
        )
    };
    let seam = chart(vec![-1.0, -1.0, 0.0, 1.0, 1.0], 1, &[-1.0, 0.0, 1.0]);
    let arch = ArchingSurface::new(
        ArchingParameters::NurbsHeightMap(seam.clone()),
        ReferenceLoadState::Unloaded,
    )
    .unwrap();
    assert!(
        seam.partials(0.0, 0.3).is_err(),
        "the ordinary partial refuses"
    );
    for y in [-1.0, 0.3, 0.75] {
        // Off the seam: the ordinary partial.
        let (h, [dx, dy]) = arch.sample([0.5, y]).unwrap();
        assert!((h - z(0.5, y)).abs() < 1e-15);
        assert!((dx - (1.0 + 0.1 * y)).abs() < 1e-14 && (dy - 0.3).abs() < 1e-14);
        // On the seam: the exact height and the mean of the one-sided slopes.
        let (h, [dx, dy]) = arch.sample([0.0, y]).unwrap();
        assert!((h - z(0.0, y)).abs() < 1e-15, "{h}");
        assert!((dx - 0.1 * y).abs() < 1e-14, "{dx}");
        assert!((dy - 0.25).abs() < 1e-12, "{dy}");
    }
    // A C1 interior knot (degree two, multiplicity one) is not a seam: the
    // sample is the ordinary partial, bit for bit.
    let smooth = chart(
        vec![-1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 1.0],
        2,
        &[-1.0, -0.5, 0.5, 1.0],
    );
    let arch = ArchingSurface::new(
        ArchingParameters::NurbsHeightMap(smooth.clone()),
        ReferenceLoadState::Unloaded,
    )
    .unwrap();
    let (point, du, dv) = smooth.partials(0.0, 0.3).unwrap();
    let (h, [dx, dy]) = arch.sample([0.0, 0.3]).unwrap();
    assert_eq!(
        [h, dx, dy].map(f64::to_bits),
        [point[2], du[2], dv[2]].map(f64::to_bits)
    );
}
