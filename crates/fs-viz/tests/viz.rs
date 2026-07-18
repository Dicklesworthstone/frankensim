//! Battery for scientific visualization (fs-viz). Each test checks a primitive
//! against ANALYTIC ground truth: rotation streamlines are circles, saddle
//! streamlines conserve xy, Hessian classification recovers the known Morse
//! type, and a circle-SDF isocontour lies on the circle.

use fs_viz::{
    CriticalKind, Grid2, Grid2Error, Grid3, Grid3Error, IsoContourError, IsoSurfaceError,
    SCALAR_FIELD3_ARTIFACT_KIND, SCALAR_FIELD3_SCHEMA_VERSION, ScalarField3, ScalarField3Error,
    ScalarFieldSemantics, ScalarLayout3, classify_hessian, streamline,
};
use std::cell::Cell;

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
    // Positive scaling preserves inertia even when directly squaring the
    // finite off-diagonal entry would overflow.
    let normalized = classify_hessian([[1.0, 1e-108], [1e-108, 1.0]], t);
    let large = classify_hessian([[1e308, 1e200], [1e200, 1e308]], t);
    assert_eq!(large.kind, CriticalKind::Minimum);
    assert_eq!(large.morse_index, 0);
    assert_eq!(large, normalized);
    // Invalid numerics cannot manufacture a confident Morse type.
    let invalid = classify_hessian([[1.0, f64::INFINITY], [f64::INFINITY, 1.0]], t);
    assert_eq!(invalid.kind, CriticalKind::Degenerate);
    assert_eq!(invalid.morse_index, 0);
    assert_eq!(
        classify_hessian([[1.0, 0.0], [0.0, 1.0]], f64::NAN).kind,
        CriticalKind::Degenerate
    );
}

#[test]
fn a_circle_sdf_isocontour_lies_on_the_circle() {
    // f(x,y) = sqrt(x²+y²) - 1, zero level set is the unit circle.
    let grid = Grid2::from_fn(41, 41, [-2.0, -2.0], [2.0, 2.0], 41 * 41, |p| {
        radius(p) - 1.0
    })
    .expect("finite circle grid within its exact node budget");
    let crossing_limit = 2 * 41 * 40;
    let crossings = grid
        .isocontour_crossings(0.0, crossing_limit)
        .expect("finite non-coincident circle crossings within edge budget");
    assert!(!crossings.is_empty());
    for c in &crossings {
        assert!(
            (radius(*c) - 1.0).abs() < 0.02,
            "crossing radius {}",
            radius(*c)
        );
    }
    // a level set outside the field's range has no crossings.
    assert!(
        grid.isocontour_crossings(100.0, crossing_limit)
            .expect("finite out-of-range levels are valid")
            .is_empty()
    );
}

#[test]
fn the_grid_samples_and_addresses_correctly() {
    let grid = Grid2::from_fn(3, 3, [-0.0, 0.0], [2.0, 2.0], 9, |p| p[0] + p[1])
        .expect("finite 3x3 grid within its exact node budget");
    let (p00, p22) = (grid.point(0, 0), grid.point(2, 2));
    assert_eq!(p00[0].to_bits(), (-0.0_f64).to_bits());
    assert!(p00[0].abs() < 1e-12 && p00[1].abs() < 1e-12);
    assert!((p22[0] - 2.0).abs() < 1e-12 && (p22[1] - 2.0).abs() < 1e-12);
    assert!((grid.at(1, 1) - 2.0).abs() < 1e-12); // (1,1) -> value 1+1
}

#[test]
fn grid2_layout_admission_precedes_sampling() {
    let calls = Cell::new(0usize);
    let mut sample = |_| {
        calls.set(calls.get() + 1);
        0.0
    };
    assert!(matches!(
        Grid2::from_fn(1, 2, [0.0; 2], [1.0; 2], 2, &mut sample),
        Err(Grid2Error::InvalidDimensions { dimensions: [1, 2] })
    ));
    assert!(matches!(
        Grid2::from_fn(2, 1, [0.0; 2], [1.0; 2], 2, &mut sample),
        Err(Grid2Error::InvalidDimensions { dimensions: [2, 1] })
    ));
    assert!(matches!(
        Grid2::from_fn(usize::MAX, 2, [0.0; 2], [1.0; 2], usize::MAX, &mut sample),
        Err(Grid2Error::NodeCountOverflow { .. })
    ));
    assert!(matches!(
        Grid2::from_fn(2, 2, [0.0; 2], [1.0; 2], 3, &mut sample),
        Err(Grid2Error::NodeBudgetExceeded {
            required: 4,
            limit: 3
        })
    ));

    let invalid_bounds = [
        ([f64::NAN, 0.0], [1.0, 1.0]),
        ([0.0, 0.0], [f64::INFINITY, 1.0]),
        ([0.0, 0.0], [0.0, 1.0]),
        ([1.0, 0.0], [0.0, 1.0]),
        ([-f64::MAX, 0.0], [f64::MAX, 1.0]),
    ];
    for (lo, hi) in invalid_bounds {
        assert!(matches!(
            Grid2::from_fn(2, 2, lo, hi, 4, &mut sample),
            Err(Grid2Error::InvalidBounds { axis: 0, .. })
        ));
    }
    let adjacent_to_one = 1.0_f64.next_up();
    assert!(matches!(
        Grid2::from_fn(3, 2, [1.0, 0.0], [adjacent_to_one, 1.0], 6, &mut sample),
        Err(Grid2Error::UnrepresentableCoordinates {
            axis: 0,
            first_index: 0,
            first,
            second_index: 1,
            second
        }) if first.to_bits() == 1.0_f64.to_bits()
            && second.to_bits() == 1.0_f64.to_bits()
    ));
    assert_eq!(calls.get(), 0, "invalid layouts must not invoke the field");
}

#[test]
fn grid2_rejects_the_first_nonfinite_sample() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let index = Cell::new(0usize);
        let result = Grid2::from_fn(3, 2, [0.0; 2], [1.0; 2], 6, |_| {
            let current = index.get();
            index.set(current + 1);
            if current == 4 { bad } else { current as f64 }
        });
        assert!(matches!(
            result,
            Err(Grid2Error::NonFiniteValue {
                index: 4,
                value
            }) if value.to_bits() == bad.to_bits()
        ));
        assert_eq!(index.get(), 5, "sampling stops at the first bad value");
    }
}

#[test]
fn isocontour_admission_distinguishes_invalid_from_empty() {
    let grid = Grid2::from_fn(2, 2, [-1.0; 2], [1.0; 2], 4, |p| p[0]).expect("finite affine grid");
    for iso in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            grid.isocontour_crossings(iso, 4),
            Err(IsoContourError::NonFiniteIso { iso: rejected })
                if rejected.to_bits() == iso.to_bits()
        ));
    }
    assert_eq!(
        grid.isocontour_crossings(0.0, 0),
        Err(IsoContourError::ZeroCrossingLimit)
    );
    assert_eq!(
        grid.isocontour_crossings(0.0, 1),
        Err(IsoContourError::CrossingBudgetExceeded { limit: 1 })
    );
    assert!(
        grid.isocontour_crossings(2.0, 1)
            .expect("a finite absent level is valid")
            .is_empty()
    );
}

#[test]
fn exact_level_nodes_are_unique_and_coincident_edges_are_refused() {
    let grid = Grid2::from_fn(3, 3, [-1.0; 2], [1.0; 2], 9, |p| p[0] + 2.0 * p[1])
        .expect("finite affine grid");
    let crossings = grid
        .isocontour_crossings(0.0, 16)
        .expect("isolated exact nodes have unique point intersections");
    assert_eq!(
        crossings
            .iter()
            .filter(|point| point[0] == 0.0 && point[1] == 0.0)
            .count(),
        1,
        "the exact center is shared by incident edges but emitted once"
    );

    let signed_zero = Grid2::from_fn(3, 2, [0.0; 2], [2.0, 1.0], 6, |p| {
        if p[0] == 1.0 && p[1] == 0.0 {
            -0.0
        } else {
            p[0] + p[1] - 1.0
        }
    })
    .expect("finite signed-zero grid");
    let signed_crossings = signed_zero
        .isocontour_crossings(0.0, 8)
        .expect("signed zero is one exact level node");
    assert_eq!(
        signed_crossings
            .iter()
            .filter(|point| point[0] == 1.0 && point[1] == 0.0)
            .count(),
        1
    );

    let plateau = Grid2::from_fn(3, 2, [-1.0, 0.0], [1.0, 1.0], 6, |p| p[0])
        .expect("finite grid containing a level-coincident vertical edge");
    assert_eq!(
        plateau.isocontour_crossings(0.0, 8),
        Err(IsoContourError::CoincidentLevelEdge {
            first: [1, 0],
            second: [1, 1]
        })
    );
}

#[test]
fn isocontour_interpolation_handles_extreme_finite_values() {
    for magnitude in [f64::MAX, f64::from_bits(1)] {
        let grid = Grid2::from_fn(2, 2, [0.0; 2], [1.0; 2], 4, |p| {
            if p[0] == 0.0 { -magnitude } else { magnitude }
        })
        .expect("finite extreme samples are admissible");
        let crossings = grid
            .isocontour_crossings(0.0, 2)
            .expect("scaled interpolation remains finite");
        assert_eq!(crossings, vec![[0.5, 0.0], [0.5, 1.0]]);
    }
}

#[test]
fn g3_isocontour_value_transformations_preserve_crossings() {
    let base = Grid2::from_fn(4, 4, [-1.0; 2], [1.0; 2], 16, |p| p[0] + 0.5 * p[1])
        .expect("finite affine base grid");
    let inverted = Grid2::from_fn(4, 4, [-1.0; 2], [1.0; 2], 16, |p| -(p[0] + 0.5 * p[1]))
        .expect("finite sign-inverted grid");
    let scaled = Grid2::from_fn(4, 4, [-1.0; 2], [1.0; 2], 16, |p| 8.0 * (p[0] + 0.5 * p[1]))
        .expect("finite power-of-two-scaled grid");

    let expected = base
        .isocontour_crossings(0.1, 24)
        .expect("base affine contour");
    assert_eq!(
        inverted
            .isocontour_crossings(-0.1, 24)
            .expect("sign inversion contour"),
        expected
    );
    assert_eq!(
        scaled
            .isocontour_crossings(0.8, 24)
            .expect("positive scaling contour"),
        expected
    );
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
    assert!(
        upper
            .into_iter()
            .all(|coordinate| (coordinate - 1.0).abs() < 1e-15)
    );
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

fn density_semantics() -> ScalarFieldSemantics {
    ScalarFieldSemantics {
        quantity: "density".to_string(),
        coordinate_unit: "m".to_string(),
        value_unit: "kg/m^3".to_string(),
    }
}

#[test]
fn scalar_field_artifact_round_trips_bit_exactly_into_node_viz() {
    assert_eq!(SCALAR_FIELD3_ARTIFACT_KIND, "frankensim.scalar-field3");
    assert_eq!(SCALAR_FIELD3_SCHEMA_VERSION, 1);
    let values: Vec<f64> = (0..12).map(|index| f64::from(index) * 0.125).collect();
    let field = ScalarField3::new(
        ScalarLayout3::NodeCentered,
        [3, 2, 2],
        [-1.0, -2.0, 0.0],
        [0.5, 2.0, 1.0],
        density_semantics(),
        12,
        values,
    )
    .expect("valid node field");
    assert_eq!(field.world_bounds(), [[-1.0, -2.0, 0.0], [0.0, 0.0, 1.0]]);
    let encoded = field.encode(4096).expect("bounded encode");
    assert_eq!(encoded, field.encode(4096).expect("replay encode"));
    let decoded = ScalarField3::decode(&encoded, 12, encoded.len()).expect("bounded decode");
    assert_eq!(decoded, field);
    assert_eq!(decoded.encode(encoded.len()).expect("re-encode"), encoded);
    let grid = decoded.into_node_grid(12).expect("node-grid conversion");
    assert_eq!(grid.dimensions(), [3, 2, 2]);
    assert_eq!(grid.bounds(), [[-1.0, -2.0, 0.0], [0.0, 0.0, 1.0]]);
    assert_eq!(grid.at(2, 1, 1), Some(11.0 * 0.125));
}

#[test]
fn scalar_field_artifact_keeps_one_cell_thick_lbm_layout_honest() {
    let values: Vec<f64> = (0..12).map(|index| f64::from(index) / 11.0).collect();
    let field = ScalarField3::new(
        ScalarLayout3::CellCentered,
        [4, 3, 1],
        [0.0; 3],
        [1.0, 1.0, 24.0],
        ScalarFieldSemantics {
            quantity: "liquid_mass_fraction".to_string(),
            coordinate_unit: "cell".to_string(),
            value_unit: "1".to_string(),
        },
        12,
        values,
    )
    .expect("valid one-cell-thick field");
    assert_eq!(field.world_bounds(), [[0.0; 3], [4.0, 3.0, 24.0]]);
    let encoded = field.encode(4096).expect("bounded encode");
    let decoded = ScalarField3::decode(&encoded, 12, 4096).expect("bounded decode");
    assert_eq!(decoded.layout(), ScalarLayout3::CellCentered);
    assert_eq!(decoded.dimensions(), [4, 3, 1]);
    assert_eq!(decoded.origin(), [0.0; 3]);
    assert_eq!(decoded.spacing(), [1.0, 1.0, 24.0]);
    assert_eq!(decoded.semantics().value_unit, "1");
    assert!(matches!(
        decoded.into_node_grid(12),
        Err(ScalarField3Error::NotNodeCentered)
    ));
}

#[test]
fn scalar_field_codec_refuses_before_unbounded_or_ambiguous_work() {
    assert!(matches!(
        ScalarField3::new(
            ScalarLayout3::NodeCentered,
            [2, 2, 2],
            [0.0; 3],
            [1.0; 3],
            density_semantics(),
            7,
            vec![0.0; 8],
        ),
        Err(ScalarField3Error::SampleBudgetExceeded {
            required: 8,
            limit: 7,
        })
    ));
    assert!(matches!(
        ScalarField3::new(
            ScalarLayout3::CellCentered,
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            ScalarFieldSemantics {
                quantity: String::new(),
                coordinate_unit: "m".to_string(),
                value_unit: "1".to_string(),
            },
            1,
            vec![0.0],
        ),
        Err(ScalarField3Error::InvalidSemantic { field: "quantity" })
    ));
    let field = ScalarField3::new(
        ScalarLayout3::CellCentered,
        [1, 1, 1],
        [0.0; 3],
        [1.0; 3],
        density_semantics(),
        1,
        vec![0.5],
    )
    .expect("small valid field");
    let encoded = field.encode(1024).expect("bounded encode");
    assert!(matches!(
        field.encode(encoded.len() - 1),
        Err(ScalarField3Error::ByteBudgetExceeded { .. })
    ));
    assert!(matches!(
        ScalarField3::decode(&encoded, 1, encoded.len() - 1),
        Err(ScalarField3Error::ByteBudgetExceeded { .. })
    ));
    assert!(matches!(
        ScalarField3::decode(&encoded[..encoded.len() - 1], 1, encoded.len()),
        Err(ScalarField3Error::Malformed { .. })
    ));

    let mut bad_magic = encoded.clone();
    bad_magic[0] ^= 1;
    assert!(matches!(
        ScalarField3::decode(&bad_magic, 1, bad_magic.len()),
        Err(ScalarField3Error::Malformed { what: "bad magic" })
    ));
    let mut future = encoded.clone();
    future[8..12].copy_from_slice(&(SCALAR_FIELD3_SCHEMA_VERSION + 1).to_le_bytes());
    assert!(matches!(
        ScalarField3::decode(&future, 1, future.len()),
        Err(ScalarField3Error::UnsupportedSchema { found: 2 })
    ));
    let mut nonfinite = encoded;
    let tail = nonfinite.len() - 8;
    nonfinite[tail..].copy_from_slice(&f64::NAN.to_bits().to_le_bytes());
    assert!(matches!(
        ScalarField3::decode(&nonfinite, 1, nonfinite.len()),
        Err(ScalarField3Error::NonFiniteValue { index: 0 })
    ));
}
