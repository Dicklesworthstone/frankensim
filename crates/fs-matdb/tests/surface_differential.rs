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
