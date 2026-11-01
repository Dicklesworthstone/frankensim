//! G0/G3: simultaneous physical coordinates retain original scalar source support.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, EvaluationDecision, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Dims, QuantitySpec};

fn point(x: f64, y: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("x", QuantitySpec::dimensional(Dims::NONE), x)
        .unwrap()
        .with_quantity("y", QuantitySpec::dimensional(Dims::NONE), y)
        .unwrap()
}
fn claim(value: PropertyValue, interpolation: InterpolationPolicy) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("response", Dims::NONE),
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity("x", QuantitySpec::dimensional(Dims::NONE), 0.0, 2.0)
            .with_quantity("y", QuantitySpec::dimensional(Dims::NONE), 0.0, 2.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured affine path".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn surface(points: &[[f64; 2]], values: &[f64], cells: Vec<[u32; 3]>) -> PropertyValue {
    PropertyValue::Surface2D {
        axes: ["x".into(), "y".into()],
        axis_dims: [Dims::NONE; 2],
        samples: points.iter().copied().zip(values.iter().copied()).collect(),
        triangles: cells,
        dims: Dims::NONE,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-13, "{a} != {b}");
}

#[test]
fn original_cell_crossing_and_point_receipts_survive_reversed_paths_and_axes() {
    for reverse_axes in [false, true] {
        let mut value = surface(
            &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            &[2.0, 4.0, 10.0, 5.0],
            vec![[0, 1, 2], [0, 2, 3]],
        );
        if reverse_axes {
            if let PropertyValue::Surface2D {
                axes,
                samples,
                triangles,
                ..
            } = &mut value
            {
                axes.swap(0, 1);
                for (point, _) in samples {
                    point.swap(0, 1);
                }
                for cell in triangles {
                    cell.swap(0, 2);
                }
            }
        }
        let mut set = ClaimSet::new();
        let pin = set
            .insert_claim(claim(value, InterpolationPolicy::LinearInside))
            .unwrap();
        let key = PropertyKey::new("response", Dims::NONE);
        let a = point(0.0, 0.25);
        let b = point(1.0, 0.75);
        for reverse in [false, true] {
            let (start, end) = if reverse { (&b, &a) } else { (&a, &b) };
            let path = set.query_path_pinned_typed(&key, start, end, pin).unwrap();
            assert_eq!(path.start, *start);
            assert_eq!(path.end, *end);
            assert_eq!(
                path.knots.iter().map(|k| k.0).collect::<Vec<_>>(),
                vec![0.0, 0.5, 1.0]
            );
            let expected = if reverse {
                [8.5, 6.0, 2.75]
            } else {
                [2.75, 6.0, 8.5]
            };
            for ((_, answer), expected) in path.knots.iter().zip(expected) {
                near(answer.evidence.value.value, expected);
                assert_eq!(answer.receipt.selected, pin);
                set.verify_receipt(&answer.receipt).unwrap();
            }
        }
    }
}

#[test]
fn shared_vertex_connectivity_covers_a_path_but_disconnected_cells_and_isolated_samples_do_not() {
    let connected = surface(
        &[[0.0, 0.0], [2.0, 0.0], [1.0, 1.0], [2.0, 2.0], [0.0, 2.0]],
        &[1.0, 3.0, 3.0, 5.0, 3.0],
        vec![[0, 1, 2], [2, 3, 4]],
    );
    let mut set = ClaimSet::new();
    let pin = set
        .insert_claim(claim(connected, InterpolationPolicy::LinearInside))
        .unwrap();
    let key = PropertyKey::new("response", Dims::NONE);
    let path = set
        .query_path_pinned_typed(&key, &point(0.25, 0.25), &point(1.75, 1.75), pin)
        .unwrap();
    assert_eq!(
        path.knots.iter().map(|k| k.0).collect::<Vec<_>>(),
        vec![0.0, 0.5, 1.0]
    );
    assert!(matches!(
        path.knots[1].1.receipt.decision,
        EvaluationDecision::ExactSurfaceSample { sample: 2 }
    ));
    let hole = surface(
        &[
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.125, 1.125],
            [2.0, 1.125],
            [1.125, 2.0],
        ],
        &[1.0; 6],
        vec![[0, 1, 2], [3, 4, 5]],
    );
    let isolated = surface(&[[0.25, 0.25], [1.25, 1.25]], &[1.0; 2], vec![]);
    for value in [hole, isolated] {
        let mut set = ClaimSet::new();
        let pin = set
            .insert_claim(claim(value, InterpolationPolicy::LinearInside))
            .unwrap();
        assert!(
            set.query_pinned_typed(&key, &point(0.25, 0.25), pin)
                .is_ok()
        );
        assert!(
            set.query_pinned_typed(&key, &point(1.25, 1.25), pin)
                .is_ok()
        );
        assert!(
            set.query_path_pinned_typed(&key, &point(0.25, 0.25), &point(1.25, 1.25), pin)
                .is_err()
        );
    }
}

#[test]
fn curved_source_knots_remain_exact_while_other_coordinates_vary() {
    let value = PropertyValue::Curve {
        abscissa: "x".into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(0.1, 2.0), (0.2, 5.0), (0.9, 4.0)],
        dims: Dims::NONE,
    };
    let mut set = ClaimSet::new();
    let pin = set
        .insert_claim(claim(value.clone(), InterpolationPolicy::LinearInside))
        .unwrap();
    let key = PropertyKey::new("response", Dims::NONE);
    let path = set
        .query_path_pinned_typed(&key, &point(0.1, 0.25), &point(0.9, 1.75), pin)
        .unwrap();
    assert_eq!(path.knots.len(), 3);
    assert_eq!(path.knots[1].1.evidence.value.value, 5.0);
    assert_eq!(
        path.knots[1]
            .1
            .receipt
            .query_point
            .iter()
            .find(|(axis, _)| axis == "x")
            .unwrap()
            .1
            .to_bits(),
        0.2_f64.to_bits()
    );
    assert!(matches!(
        path.knots[1].1.receipt.decision,
        EvaluationDecision::ExactTabulated { at: 0.2 }
    ));
    set.verify_receipt(&path.knots[1].1.receipt).unwrap();
    let mut discrete = ClaimSet::new();
    let pin = discrete
        .insert_claim(claim(value, InterpolationPolicy::TabulatedOnly))
        .unwrap();
    assert!(
        discrete
            .query_path_pinned_typed(&key, &point(0.1, 0.25), &point(0.9, 1.75), pin)
            .is_err()
    );
}

#[test]
fn fixed_response_and_single_axis_paths_share_existing_owners_and_validate_full_context() {
    let mut set = ClaimSet::new();
    let pin = set
        .insert_claim(claim(
            PropertyValue::Scalar {
                value: 7.0,
                dims: Dims::NONE,
            },
            InterpolationPolicy::ConstantWithinValidity,
        ))
        .unwrap();
    let key = PropertyKey::new("response", Dims::NONE);
    for (start, end) in [
        (point(0.25, 0.5), point(0.25, 0.5)),
        (point(1.5, 0.5), point(0.25, 0.5)),
    ] {
        let path = set
            .query_path_pinned_typed(&key, &start, &end, pin)
            .unwrap();
        assert_eq!(path.knots.len(), 2);
        assert_eq!(path.knots[0].0, 0.0);
        assert_eq!(path.knots[1].0, 1.0);
        assert!(
            path.knots
                .iter()
                .all(|(_, answer)| answer.evidence.value.value == 7.0)
        );
    }
    assert!(
        set.query_path_pinned_typed(&key, &point(0.0, 0.0), &point(1.0, 3.0), pin)
            .is_err()
    );
    let incomplete = QueryPoint::new()
        .with_quantity("x", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap();
    assert!(
        set.query_path_pinned_typed(&key, &point(0.0, 0.0), &incomplete, pin)
            .is_err()
    );
}
