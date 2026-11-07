//! Original source slices, pressure edges and inward temperature endpoints.
use fs_matdb::{PreparedSurfaceResponse, PropertyValue};
use fs_qty::Dims;

fn surface(values: [f64; 4], triangles: Vec<[u32; 3]>) -> PreparedSurfaceResponse {
    PropertyValue::Surface2D {
        axes: ["temperature".into(), "pressure".into()],
        axis_dims: [Dims::NONE; 2],
        samples: [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [2.0, 2.0]]
            .into_iter()
            .zip(values)
            .collect(),
        triangles,
        dims: Dims::NONE,
    }
    .prepare_surface_response()
    .unwrap()
}

#[test]
fn fixed_pressure_slices_keep_inward_temperature_support_without_extrapolation() {
    for triangles in [vec![[0, 1, 2], [1, 3, 2]], vec![[0, 2, 1], [1, 2, 3]]] {
        let affine = surface([4.0, 8.0, -2.0, 2.0], triangles);
        for t in [0.0, 0.5, 1.0, 2.0] {
            for p in [0.0, 1.0, 2.0] {
                assert_eq!(
                    affine.inward_axis_derivative([t, p], 0, [true; 2]).unwrap(),
                    2.0
                );
                assert_eq!(
                    affine.inward_axis_derivative([t, p], 1, [true; 2]).unwrap(),
                    -3.0
                );
            }
        }
        // The existing classical derivative still requires both thermal sides.
        assert!(affine.gradient_for_axes([0.0, 1.0], [true, false]).is_err());
        assert!(
            affine
                .inward_axis_derivative([0.0, 1.0], 0, [true, false])
                .is_err()
        );
        assert!(
            affine
                .inward_axis_derivative([2.0, 1.0], 0, [false, true])
                .is_err()
        );
        assert!(
            affine
                .inward_axis_derivative([1.0, 1.0], 0, [false; 2])
                .is_err()
        );
        assert!(
            affine
                .inward_axis_derivative([1.0, 1.0], 2, [true; 2])
                .is_err()
        );
    }
}

#[test]
fn rounding_agreement_does_not_hide_genuine_source_corners_or_missing_cells() {
    // A decimal affine plane can produce adjacent binary slopes a few ulps
    // apart; it needs the same relative rounding policy as original 1D curves.
    let affine = surface(
        [0.03, 0.0302, 0.0300002, 0.0302002],
        vec![[0, 1, 2], [1, 3, 2]],
    );
    let slope = affine
        .inward_axis_derivative([1.0, 1.0], 0, [true; 2])
        .unwrap();
    assert!((slope - 0.0001).abs() < 128.0 * f64::EPSILON * slope.abs());

    let corner = surface([0.0, 2.0, 4.0, 10.0], vec![[0, 1, 2], [1, 3, 2]]);
    assert!(
        corner
            .inward_axis_derivative([1.0, 1.0], 0, [true; 2])
            .is_err()
    );
    // Original claim-validity endpoints can select just their admitted side.
    assert_eq!(
        corner
            .inward_axis_derivative([1.0, 1.0], 0, [true, false])
            .unwrap(),
        1.0
    );
    assert_eq!(
        corner
            .inward_axis_derivative([1.0, 1.0], 0, [false, true])
            .unwrap(),
        3.0
    );

    let hole = surface([0.0, 2.0, 4.0, 10.0], vec![[0, 1, 2]]);
    assert!(
        hole.inward_axis_derivative([1.5, 1.5], 0, [true; 2])
            .is_err()
    );
    assert!(
        hole.inward_axis_derivative([2.0, 2.0], 0, [true; 2])
            .is_err()
    );
    assert!(
        hole.inward_axis_derivative([0.0, 2.0], 0, [false, true])
            .is_err()
    );
}
