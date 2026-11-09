use super::*;
use fs_material::thermal::transport::expansion::surface::MAX_ORIENTED_EXPANSION_PATCHES;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    a.into_iter().map(|v| v * v).sum::<f64>().sqrt()
}

#[test]
fn oblique_face_matches_the_actual_deformed_triangle_and_its_area_tangent() {
    let alpha = [0.006, -0.002, 0.003];
    let values = alpha.map(|a| vec![(300.0, 0.0), (400.0, 100.0 * a)]);
    let path = source(300.0, &values);
    let q = [2.0, 1.0, -1.0];
    let r = [-1.0, 2.0, 3.0];
    let a0 = cross(q, r).map(|v| 0.5 * v);
    let aligned = [1.0, 2.0, 3.0];
    let surface = PreparedExpansionSurface::try_new(path, aligned)
        .unwrap()
        .with_oriented_patches(&[a0])
        .unwrap();
    near(surface.reference_total_area_m2(), 6.0 + norm(a0));
    assert_eq!(surface.reference_oriented_area_vectors_m2(), &[a0]);
    let triangle_area = |t: f64| {
        let s = alpha.map(|a| 1.0 + a * (t - 300.0));
        let q = core::array::from_fn(|i| q[i] * s[i]);
        let r = core::array::from_fn(|i| r[i] * s[i]);
        0.5 * norm(cross(q, r))
            + (0..3)
                .map(|i| aligned[i] * s[(i + 1) % 3] * s[(i + 2) % 3])
                .sum::<f64>()
    };
    for t in [300.0, 337.0, 400.0] {
        let sample = surface.at_temperature(t).unwrap();
        near(sample.total_area_m2(), triangle_area(t));
        let numerical = (triangle_area(t + 0.001) - triangle_area(t - 0.001)) / 0.002;
        assert!((sample.total_area_derivative_m2_k() - numerical).abs() < 5e-9);
    }
}

#[test]
fn oriented_area_preserves_isotropy_normal_sign_and_principal_axis_permutation() {
    let vectors = [[1.0, 2.0, -3.0], [-2.0, 1.0, 4.0]];
    let isotropic = PreparedExpansionSurface::try_new(
        source(300.0, &[vec![(300.0, 0.0), (400.0, 0.2)]]),
        [0.0; 3],
    )
    .unwrap()
    .with_oriented_patches(&vectors)
    .unwrap();
    let original_area = vectors.into_iter().map(norm).sum::<f64>();
    let sample = isotropic.at_temperature(350.0).unwrap();
    near(sample.total_area_m2(), original_area * 1.1 * 1.1);
    near(
        sample.total_area_derivative_m2_k(),
        original_area * 2.0 * 1.1 * 0.002,
    );

    let values = [
        vec![(300.0, 0.0), (400.0, 0.3)],
        vec![(300.0, 0.0), (400.0, -0.2)],
        vec![(300.0, 0.0), (400.0, 0.1)],
    ];
    let a = PreparedExpansionSurface::try_new(source(300.0, &values), [0.0; 3])
        .unwrap()
        .with_oriented_patches(&vectors)
        .unwrap()
        .at_temperature(365.0)
        .unwrap();
    let rotated = vectors.map(|v| [-v[1], -v[2], -v[0]]);
    let reordered = [values[1].clone(), values[2].clone(), values[0].clone()];
    let b = PreparedExpansionSurface::try_new(source(300.0, &reordered), [0.0; 3])
        .unwrap()
        .with_oriented_patches(&rotated)
        .unwrap()
        .at_temperature(365.0)
        .unwrap();
    near(a.total_area_m2(), b.total_area_m2());
    near(a.total_area_derivative_m2_k(), b.total_area_derivative_m2_k());
}

#[test]
fn oriented_patches_refuse_invalid_geometry_and_keep_representable_small_areas() {
    let path = source(300.0, &[vec![(300.0, 0.0), (400.0, 0.1)]]);
    let surface = PreparedExpansionSurface::try_new(path, [0.0; 3]).unwrap();
    for vector in [
        [0.0; 3],
        [f64::NAN, 1.0, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        [f64::MAX; 3],
    ] {
        assert!(surface.clone().with_oriented_patches(&[vector]).is_err());
    }
    assert!(
        surface
            .clone()
            .with_oriented_patches(&vec![[1.0, 0.0, 0.0]; MAX_ORIENTED_EXPANSION_PATCHES + 1])
            .is_err()
    );
    let tiny = surface
        .clone()
        .with_oriented_patches(&[[1e-300, 0.0, 0.0]])
        .unwrap();
    let area = tiny.at_temperature(350.0).unwrap().total_area_m2();
    assert!((area / 1e-300 - 1.05 * 1.05).abs() < 1e-14);
    let large = surface
        .with_oriented_patches(&[[f64::MAX * 0.9, 0.0, 0.0]])
        .unwrap();
    assert!(large.at_temperature(400.0).is_err());
}
