//! Battery for scientific visualization (fs-viz). Each test checks a primitive
//! against ANALYTIC ground truth: rotation streamlines are circles, saddle
//! streamlines conserve xy, Hessian classification recovers the known Morse
//! type, and a circle-SDF isocontour lies on the circle.

use fs_viz::{
    CriticalKind, Grid2, Grid3, Grid3Error, IsoSurfaceError, classify_hessian, streamline,
};

fn radius(p: [f64; 2]) -> f64 {
    (p[0] * p[0] + p[1] * p[1]).sqrt()
}

#[test]
fn a_rotation_field_streams_along_a_circle() {
    // u = (-y, x): rigid rotation, so the radius is conserved.
    let line = streamline(|p| [-p[1], p[0]], [1.0, 0.0], 0.01, 400);
    for p in &line {
        assert!(
            (radius(*p) - 1.0).abs() < 1e-3,
            "radius {} drifted",
            radius(*p)
        );
    }
    // it actually goes somewhere (not a fixed point).
    assert!((line.last().unwrap()[1]).abs() > 0.1);
}

#[test]
fn a_saddle_field_conserves_the_hyperbola_invariant() {
    // u = (x, -y): flow x·y is invariant along a streamline.
    let line = streamline(|p| [p[0], -p[1]], [1.0, 1.0], 0.01, 50);
    for p in &line {
        assert!(
            (p[0] * p[1] - 1.0).abs() < 1e-4,
            "xy = {} drifted",
            p[0] * p[1]
        );
    }
    // x grows, y shrinks (the saddle's unstable/stable manifolds).
    assert!(line.last().unwrap()[0] > 1.4 && line.last().unwrap()[1] < 0.7);
}

#[test]
fn hessian_classification_recovers_the_morse_type() {
    let t = 1e-9;
    // f = x² + y²  -> minimum, index 0.
    assert_eq!(
        classify_hessian([[2.0, 0.0], [0.0, 2.0]], t).kind,
        CriticalKind::Minimum
    );
    // f = x² - y²  -> saddle, index 1.
    let s = classify_hessian([[2.0, 0.0], [0.0, -2.0]], t);
    assert_eq!(s.kind, CriticalKind::Saddle);
    assert_eq!(s.morse_index, 1);
    // f = -(x² + y²) -> maximum, index 2.
    assert_eq!(
        classify_hessian([[-2.0, 0.0], [0.0, -2.0]], t).morse_index,
        2
    );
    // f = xy -> saddle (off-diagonal Hessian, eigenvalues ±1).
    assert_eq!(
        classify_hessian([[0.0, 1.0], [1.0, 0.0]], t).kind,
        CriticalKind::Saddle
    );
    // a zero eigenvalue is degenerate.
    assert_eq!(
        classify_hessian([[2.0, 0.0], [0.0, 0.0]], t).kind,
        CriticalKind::Degenerate
    );
}

#[test]
fn a_circle_sdf_isocontour_lies_on_the_circle() {
    // f(x,y) = sqrt(x²+y²) - 1, zero level set is the unit circle.
    let grid = Grid2::from_fn(41, 41, [-2.0, -2.0], [2.0, 2.0], |p| radius(p) - 1.0);
    let crossings = grid.isocontour_crossings(0.0);
    assert!(!crossings.is_empty());
    for c in &crossings {
        assert!(
            (radius(*c) - 1.0).abs() < 0.02,
            "crossing radius {}",
            radius(*c)
        );
    }
    // a level set outside the field's range has no crossings.
    assert!(grid.isocontour_crossings(100.0).is_empty());
}

#[test]
fn the_grid_samples_and_addresses_correctly() {
    let grid = Grid2::from_fn(3, 3, [0.0, 0.0], [2.0, 2.0], |p| p[0] + p[1]);
    let (p00, p22) = (grid.point(0, 0), grid.point(2, 2));
    assert!(p00[0].abs() < 1e-12 && p00[1].abs() < 1e-12);
    assert!((p22[0] - 2.0).abs() < 1e-12 && (p22[1] - 2.0).abs() < 1e-12);
    assert!((grid.at(1, 1) - 2.0).abs() < 1e-12); // (1,1) -> value 1+1
}

#[test]
fn visualization_is_deterministic() {
    let a = streamline(|p| [-p[1], p[0]], [1.0, 0.0], 0.01, 100);
    let b = streamline(|p| [-p[1], p[0]], [1.0, 0.0], 0.01, 100);
    assert_eq!(a.len(), b.len());
    assert_eq!(
        a.last().unwrap()[0].to_bits(),
        b.last().unwrap()[0].to_bits()
    );
}

#[test]
fn marching_tetrahedra_extracts_an_exact_oriented_plane() {
    let dimensions = [9, 10, 11];
    let node_limit = dimensions.into_iter().product();
    let grid = Grid3::from_fn(dimensions, [-1.0; 3], [1.0; 3], node_limit, |point| {
        point[0] - 0.13
    })
    .expect("bounded finite plane grid");
    assert_eq!(grid.dimensions(), dimensions);
    assert!((grid.at(0, 0, 0).expect("in bounds") + 1.13).abs() < 1e-15);
    let upper = grid.point(8, 9, 10).expect("upper node is in bounds");
    assert!(upper.into_iter().all(|coordinate| (coordinate - 1.0).abs() < 1e-15));
    assert_eq!(grid.point(9, 0, 0), None);

    let mesh = grid.isosurface(0.0, 10_000).expect("plane isosurface");
    assert!(!mesh.triangles().is_empty());
    assert!(mesh.vertices().len() < mesh.triangles().len() * 3);
    assert!((mesh.surface_area() - 4.0).abs() < 1e-12);
    for vertex in mesh.vertices() {
        assert!((vertex[0] - 0.13).abs() < 1e-12);
    }
    for triangle in mesh.triangles() {
        let a = mesh.vertices()[triangle[0] as usize];
        let b = mesh.vertices()[triangle[1] as usize];
        let c = mesh.vertices()[triangle[2] as usize];
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let normal_x = ab[1] * ac[2] - ab[2] * ac[1];
        assert!(
            normal_x > 0.0,
            "plane triangle must point toward increasing field"
        );
    }
    assert!(matches!(
        grid.isosurface(0.0, 1),
        Err(IsoSurfaceError::TriangleBudgetExceeded { limit: 1 })
    ));
}

#[test]
fn sphere_isosurface_area_converges_under_refinement() {
    let radius = 0.7;
    let sphere = |resolution: usize| {
        let dimensions = [resolution; 3];
        let node_limit = dimensions.into_iter().product();
        Grid3::from_fn(dimensions, [-1.2; 3], [1.2; 3], node_limit, |point| {
            point[0]
                .mul_add(point[0], point[1].mul_add(point[1], point[2] * point[2]))
                .sqrt()
                - radius
        })
        .expect("bounded finite sphere grid")
        .isosurface(0.0, 200_000)
        .expect("sphere isosurface")
    };
    let coarse = sphere(17);
    let fine = sphere(33);
    let exact_area = 4.0 * std::f64::consts::PI * radius * radius;
    let coarse_error = (coarse.surface_area() - exact_area).abs();
    let fine_error = (fine.surface_area() - exact_area).abs();
    assert!(
        fine_error < coarse_error,
        "sphere area must converge: coarse {coarse_error:.3e}, fine {fine_error:.3e}"
    );
    assert!(fine_error / exact_area < 0.03);

    // Negative values are inside the sphere, so outward winding must align
    // each nondegenerate face normal with its centroid radius.
    for triangle in fine.triangles() {
        let a = fine.vertices()[triangle[0] as usize];
        let b = fine.vertices()[triangle[1] as usize];
        let c = fine.vertices()[triangle[2] as usize];
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let normal = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let centroid = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ];
        let orientation = normal[0].mul_add(
            centroid[0],
            normal[1].mul_add(centroid[1], normal[2] * centroid[2]),
        );
        assert!(orientation > 0.0);
    }
}

#[test]
fn gyroid_extraction_is_indexed_symmetric_and_deterministic() {
    let dimensions = [19; 3];
    let node_limit = dimensions.into_iter().product();
    let bound = std::f64::consts::PI;
    let grid = Grid3::from_fn(dimensions, [-bound; 3], [bound; 3], node_limit, |point| {
        point[0].sin() * point[1].cos()
            + point[1].sin() * point[2].cos()
            + point[2].sin() * point[0].cos()
    })
    .expect("bounded finite gyroid grid");
    let first = grid.isosurface(0.0, 100_000).expect("gyroid surface");
    let replay = grid.isosurface(0.0, 100_000).expect("gyroid replay");
    assert_eq!(first, replay);
    assert!(!first.triangles().is_empty());
    assert!(first.vertices().len() < first.triangles().len() * 3);

    let mut lower = [f64::INFINITY; 3];
    let mut upper = [f64::NEG_INFINITY; 3];
    for vertex in first.vertices() {
        for axis in 0..3 {
            lower[axis] = lower[axis].min(vertex[axis]);
            upper[axis] = upper[axis].max(vertex[axis]);
        }
    }
    for axis in 0..3 {
        assert!((lower[axis] + upper[axis]).abs() < 1e-12);
    }
}

#[test]
fn grid3_admission_fails_before_unbounded_or_nonfinite_work() {
    let calls = std::cell::Cell::new(0usize);
    let over_budget = Grid3::from_fn([100, 100, 100], [-1.0; 3], [1.0; 3], 1_000, |_| {
        calls.set(calls.get() + 1);
        0.0
    });
    assert!(matches!(
        over_budget,
        Err(Grid3Error::NodeBudgetExceeded {
            required: 1_000_000,
            limit: 1_000
        })
    ));
    assert_eq!(calls.get(), 0);
    assert!(matches!(
        Grid3::from_values([2, 2, 2], [-1.0; 3], [1.0; 3], 8, vec![0.0; 7]),
        Err(Grid3Error::ValueCountMismatch {
            expected: 8,
            actual: 7
        })
    ));
    assert!(matches!(
        Grid3::from_fn([2, 2, 2], [-1.0; 3], [1.0; 3], 8, |_| f64::NAN),
        Err(Grid3Error::NonFiniteValue { index: 0, .. })
    ));
    let grid = Grid3::from_fn([2, 2, 2], [-1.0; 3], [1.0; 3], 8, |point| point[0])
        .expect("small admitted grid");
    assert!(matches!(
        grid.isosurface(f64::INFINITY, 10),
        Err(IsoSurfaceError::NonFiniteIso { .. })
    ));
    assert!(matches!(
        grid.isosurface(0.0, 0),
        Err(IsoSurfaceError::ZeroTriangleLimit)
    ));
}
