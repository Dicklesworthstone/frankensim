//! G0/G1/G3: canonical values and source-cell analytical tangents.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, EvaluationDecision, InterpolationPolicy, MatDbError, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, QuantitySpec};

fn source(values: [f64; 4], triangles: Vec<[u32; 3]>) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("response", Dims::NONE),
        value: PropertyValue::Surface2D {
            axes: ["x".into(), "y".into()],
            axis_dims: [Dims::NONE; 2],
            samples: [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [2.0, 2.0]]
                .into_iter()
                .zip(values)
                .collect(),
            triangles,
            dims: Dims::NONE,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("x", QuantitySpec::dimensional(Dims::NONE), 0.0, 2.0)
            .with_quantity("y", QuantitySpec::dimensional(Dims::NONE), 0.0, 2.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured source-cell response, not measured material data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

#[test]
fn affine_gradients_preserve_canonical_values_and_original_cell_orientation() {
    for triangles in [vec![[0, 1, 2], [1, 3, 2]], vec![[0, 2, 1], [1, 2, 3]]] {
        let claim = source([4.0, 8.0, -2.0, 2.0], triangles);
        let prepared = claim.value.prepare_surface_response().unwrap();
        let mut claims = ClaimSet::new();
        let id = claims.insert_claim(claim).unwrap();
        for x in [0.0, 0.25, 1.0, 1.75, 2.0] {
            for y in [0.0, 0.25, 1.0, 1.75, 2.0] {
                let query = QueryPoint::new()
                    .with_quantity("x", QuantitySpec::dimensional(Dims::NONE), x)
                    .unwrap()
                    .with_quantity("y", QuantitySpec::dimensional(Dims::NONE), y)
                    .unwrap();
                let canonical = claims.query_pinned("response", &query, id).unwrap();
                let state = prepared.evaluate([x, y]).unwrap();
                assert_eq!(
                    state.value.to_bits(),
                    canonical.evidence.value.value.to_bits()
                );
                assert_eq!(state.gradient, [2.0, -3.0]);
                assert!((state.value - (4.0 + 2.0 * x - 3.0 * y)).abs() < 1e-13);
            }
        }
    }
}

#[test]
fn shared_edges_and_samples_have_explicit_original_one_sided_tangents() {
    let prepared = source([0.0, 2.0, 4.0, 10.0], vec![[0, 1, 2], [1, 3, 2]])
        .value
        .prepare_surface_response()
        .unwrap();
    for point in [[0.25, 0.25], [1.0, 1.0], [2.0, 0.0], [0.0, 2.0]] {
        let sample = prepared.evaluate(point).unwrap();
        assert_eq!(sample.triangle, 0);
        assert_eq!(sample.gradient, [1.0, 2.0]);
    }
    for point in [[1.75, 1.75], [2.0, 2.0]] {
        let sample = prepared.evaluate(point).unwrap();
        assert_eq!(sample.triangle, 1);
        assert_eq!(sample.gradient, [3.0, 4.0]);
    }
    assert_eq!(
        prepared.evaluate([2.0, 2.0]).unwrap().evaluation,
        EvaluationDecision::ExactSurfaceSample { sample: 3 }
    );
}

#[test]
fn missing_cells_and_isolated_samples_do_not_gain_derivative_support() {
    let prepared = source([4.0, 8.0, -2.0, 2.0], vec![[0, 1, 2]])
        .value
        .prepare_surface_response()
        .unwrap();
    assert!(matches!(
        prepared.evaluate([1.5, 1.5]),
        Err(MatDbError::OutsideSurfaceSupport { .. })
    ));
    assert!(matches!(
        prepared.evaluate([2.0, 2.0]),
        Err(MatDbError::UnsupportedEvaluation { .. })
    ));
    assert_eq!(prepared.gradient_for_axes([2.0, 2.0], [false; 2]).unwrap(), [0.0; 2]);
    assert!(prepared.evaluate([f64::NAN, 0.0]).is_err());
    assert!(
        source([1.0; 4], vec![])
            .value
            .prepare_surface_response()
            .is_err()
    );
    let constant = source([f64::MAX; 4], vec![[0, 1, 2], [1, 3, 2]])
        .value
        .prepare_surface_response()
        .unwrap();
    assert_eq!(constant.evaluate([1.0, 1.0]).unwrap().gradient, [0.0; 2]);
}

fn response(samples: Vec<([f64; 2], f64)>, triangles: Vec<[u32; 3]>) -> fs_matdb::PreparedSurfaceResponse {
    PropertyValue::Surface2D { axes: ["x".into(), "y".into()], axis_dims: [Dims::NONE; 2],
        samples, triangles, dims: Dims::NONE }.prepare_surface_response().unwrap()
}

#[test]
fn classical_gradients_follow_the_requested_source_subspace() {
    for triangles in [vec![[0,1,2],[1,3,2]], vec![[0,2,1],[1,2,3]]] {
        let affine = source([4.0,8.0,-2.0,2.0], triangles).value.prepare_surface_response().unwrap();
        assert_eq!(affine.gradient_for_axes([1.0,1.0], [true;2]).unwrap(), [2.0,-3.0]);
        assert!(affine.gradient_for_axes([0.0,1.0], [true;2]).is_err());
        assert!(affine.gradient_for_axes([0.0,1.0], [true,false]).is_err());
        assert_eq!(affine.gradient_for_axes([0.0,1.0], [false,true]).unwrap(), [0.0,-3.0]);
    }
    // The pressure-like y slope changes at y=0, while the thermal x slice
    // is the same affine line on both source cells.
    let crease = response(vec![([-1.0,0.0],1.0),([1.0,0.0],3.0),
        ([0.0,-1.0],1.0),([0.0,1.0],5.0)], vec![[0,1,2],[0,3,1]]);
    assert_eq!(crease.gradient_for_axes([0.0,0.0], [true,false]).unwrap(), [1.0,0.0]);
    assert!(crease.gradient_for_axes([0.0,0.0], [false,true]).is_err());
    assert!(crease.gradient_for_axes([0.0,0.0], [true;2]).is_err());
}

#[test]
fn incident_wedges_cannot_hide_behind_coordinate_neighbours() {
    let mut samples: Vec<_> = [[1.0,1.0],[2.0,1.0],[2.0,2.0],[1.0,2.0],[0.0,2.0],
        [0.0,1.0],[0.0,0.0],[1.0,0.0],[2.0,0.0]].into_iter().map(|p| (p,1.0)).collect();
    let triangles: Vec<_> = (1..=8).map(|i| [0,i,if i == 8 {1} else {i+1}]).collect();
    let flat = response(samples.clone(), triangles.clone());
    assert_eq!(flat.gradient_for_axes([1.0,1.0], [true;2]).unwrap(), [0.0;2]);
    // Every coordinate ray remains supported and constant even though a
    // missing wedge leaves no open two-coordinate source neighbourhood.
    let hole = response(samples.clone(), triangles[1..].to_vec());
    for p in [[0.75,1.0],[1.25,1.0],[1.0,0.75],[1.0,1.25]] {
        assert_eq!(hole.evaluate(p).unwrap().value, 1.0);
    }
    assert!(hole.gradient_for_axes([1.0,1.0], [true;2]).is_err());
    assert_eq!(hole.gradient_for_axes([1.0,1.0], [true,false]).unwrap(), [0.0;2]);
    // Restoring the wedge with a distinct ordinate creates a genuine corner
    // that is still invisible to axis-only finite differences.
    samples[2].1 = 2.0;
    let corner = response(samples, triangles);
    for active in [[true,false],[false,true]] {
        assert_eq!(corner.gradient_for_axes([1.0,1.0], active).unwrap(), [0.0;2]);
    }
    assert!(corner.gradient_for_axes([1.0,1.0], [true;2]).is_err());
}
