//! G0/G1/G3: simultaneous state dependence, explicit sampled support, and
//! portable source identity through the canonical material evaluator.
#![allow(clippy::float_cmp)]

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, EvaluationDecision, InterpolationPolicy, JointStatistics,
    MATDB_PACK_TARGET_BASIS, MATDB_SURFACE_PACK_SCHEMA_VERSION, MatDbError, NormalizationReceipt,
    NormalizationTarget, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyUsageReceipt, PropertyValue, Provenance, QueryPoint, StatisticMember, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);
const VISCOSITY: Dims = Dims([-1, 1, -1, 0, 0, 0]);

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn source() -> Provenance {
    Provenance {
        source: "synthetic manufactured viscosity(T, shear_rate); not measured material data"
            .into(),
        license: "CC0-1.0".into(),
        artifact: Some(hash_domain(
            "surface-test-source",
            b"manufactured-affine-law",
        )),
    }
}

fn law(t: f64, shear: f64) -> f64 {
    4.0 - 0.01 * (t - 300.0) - 0.005 * shear
}

fn claim() -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("dynamic_viscosity", VISCOSITY),
        value: PropertyValue::Surface2D {
            axes: ["T".into(), "shear_rate".into()],
            axis_dims: [TEMPERATURE, RATE],
            samples: [[300.0, 0.0], [400.0, 0.0], [300.0, 100.0], [400.0, 100.0]]
                .into_iter()
                .map(|point| (point, law(point[0], point[1])))
                .collect(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: VISCOSITY,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", temperature(), 300.0, 400.0)
            .with_quantity("shear_rate", QuantitySpec::dimensional(RATE), 0.0, 100.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: source(),
    }
}

fn point(t: f64, shear: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
        .with_quantity("shear_rate", QuantitySpec::dimensional(RATE), shear)
        .unwrap()
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-14 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_one_source_law_changes_with_both_temperature_and_shear() {
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(claim()).unwrap();
    for t in [300.0, 312.5, 350.0, 387.5, 400.0] {
        for shear in [0.0, 12.5, 50.0, 87.5, 100.0] {
            let answer = claims
                .query_pinned("dynamic_viscosity", &point(t, shear), id)
                .unwrap();
            assert_close(answer.evidence.value.value, law(t, shear));
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert_eq!(answer.receipt.selected, id);
            claims.verify_receipt(&answer.receipt).unwrap();
        }
    }
    assert_eq!(
        claims.claim(id).unwrap().value.axes().collect::<Vec<_>>(),
        vec![("T", TEMPERATURE), ("shear_rate", RATE)]
    );
    let exact = claims
        .query_pinned("dynamic_viscosity", &point(400.0, 100.0), id)
        .unwrap();
    assert_eq!(
        exact.receipt.decision,
        EvaluationDecision::ExactSurfaceSample { sample: 3 }
    );
    let inside = claims
        .query_pinned("dynamic_viscosity", &point(350.0, 25.0), id)
        .unwrap();
    assert_eq!(
        inside.receipt.decision,
        EvaluationDecision::LinearInsideTriangle { triangle: 0 }
    );
}

#[test]
fn g1_fixed_context_line_retains_cell_crossings_and_canonical_receipts() {
    let mut source = claim();
    if let PropertyValue::Surface2D { samples, .. } = &mut source.value {
        samples[3].1 += 1.0;
    }
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(source).unwrap();
    let key = PropertyKey::new("dynamic_viscosity", VISCOSITY);
    let slice = claims
        .query_line_typed(
            &key,
            "T",
            &point(300.0, 50.0),
            &point(400.0, 50.0),
            ClaimSelection::Pinned(id),
        )
        .unwrap();
    assert_eq!(slice.axis, "T");
    assert_eq!(
        slice.knots.iter().map(|(x, _)| *x).collect::<Vec<_>>(),
        [300.0, 350.0, 400.0]
    );
    for ((t, answer), expected) in slice.knots.iter().zip([3.75, 3.25, 3.25]) {
        assert_close(answer.evidence.value.value, expected);
        assert_eq!(
            *answer,
            claims
                .query_pinned_typed(&key, &point(*t, 50.0), id)
                .unwrap()
        );
        claims.verify_receipt(&answer.receipt).unwrap();
    }
    // Only the path has support authority; diagonal box corners are unchanged.
    assert!(
        claims
            .query_envelope_typed(
                &key,
                &point(300.0, 0.0),
                &point(400.0, 100.0),
                ClaimSelection::Pinned(id)
            )
            .is_err()
    );
    assert!(
        claims
            .query_line_typed(
                &key,
                "T",
                &point(300.0, 0.0),
                &point(400.0, 100.0),
                ClaimSelection::Pinned(id)
            )
            .is_err()
    );
    assert!(
        claims
            .query_line_typed(
                &key,
                "T",
                &point(400.0, 50.0),
                &point(300.0, 50.0),
                ClaimSelection::Pinned(id)
            )
            .is_err()
    );
}

#[test]
fn g0_line_support_does_not_bridge_even_a_sub_tolerance_source_gap() {
    let mut source = claim();
    if let PropertyValue::Surface2D {
        samples, triangles, ..
    } = &mut source.value
    {
        let right_start = f64::from_bits(350_f64.to_bits() + 1);
        *samples = [
            [300.0, 0.0],
            [350.0, 0.0],
            [300.0, 100.0],
            [350.0, 100.0],
            [right_start, 0.0],
            [400.0, 0.0],
            [right_start, 100.0],
            [400.0, 100.0],
        ]
        .map(|p| (p, law(p[0], p[1])))
        .to_vec();
        *triangles = vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]];
    }
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(source).unwrap();
    let key = PropertyKey::new("dynamic_viscosity", VISCOSITY);
    for t in [300.0, 350.0, f64::from_bits(350_f64.to_bits() + 1), 400.0] {
        claims
            .query_pinned_typed(&key, &point(t, 50.0), id)
            .unwrap();
    }
    assert!(
        claims
            .query_line_typed(
                &key,
                "T",
                &point(300.0, 50.0),
                &point(400.0, 50.0),
                ClaimSelection::Pinned(id)
            )
            .is_err()
    );
    claims
        .query_line_typed(
            &key,
            "T",
            &point(300.0, 50.0),
            &point(350.0, 50.0),
            ClaimSelection::Pinned(id),
        )
        .unwrap();
}

#[test]
fn g3_line_support_joins_real_shared_vertices_and_respects_discrete_policy() {
    let mut source = claim();
    if let PropertyValue::Surface2D {
        samples, triangles, ..
    } = &mut source.value
    {
        *samples = [
            [300.0, 0.0],
            [300.0, 100.0],
            [350.0, 50.0],
            [400.0, 0.0],
            [400.0, 100.0],
        ]
        .map(|p| (p, law(p[0], p[1])))
        .to_vec();
        *triangles = vec![[0, 2, 1], [2, 3, 4]];
    }
    let key = PropertyKey::new("dynamic_viscosity", VISCOSITY);
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(source.clone()).unwrap();
    let line = claims
        .query_line_typed(
            &key,
            "T",
            &point(300.0, 50.0),
            &point(400.0, 50.0),
            ClaimSelection::Pinned(id),
        )
        .unwrap();
    assert_eq!(
        line.knots.iter().map(|(x, _)| *x).collect::<Vec<_>>(),
        [300.0, 350.0, 400.0]
    );
    source.interpolation = InterpolationPolicy::TabulatedOnly;
    let mut discrete = ClaimSet::new();
    let id = discrete.insert_claim(source).unwrap();
    let vertex = point(350.0, 50.0);
    assert_eq!(
        discrete
            .query_line_typed(&key, "T", &vertex, &vertex, ClaimSelection::Pinned(id))
            .unwrap()
            .knots
            .len(),
        1
    );
    assert!(
        discrete
            .query_line_typed(
                &key,
                "T",
                &point(300.0, 0.0),
                &point(400.0, 0.0),
                ClaimSelection::Pinned(id)
            )
            .is_err()
    );
}

#[test]
fn g0_line_selection_cannot_cross_a_competing_claim_unless_explicitly_pinned() {
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(claim()).unwrap();
    let mut competitor = claim();
    competitor.value = PropertyValue::Scalar {
        value: 9.0,
        dims: VISCOSITY,
    };
    competitor.interpolation = InterpolationPolicy::ConstantWithinValidity;
    competitor.validity = competitor
        .validity
        .with_quantity("T", temperature(), 340.0, 360.0);
    claims.insert_claim(competitor).unwrap();
    let key = PropertyKey::new("dynamic_viscosity", VISCOSITY);
    assert!(
        claims
            .query_line_typed(
                &key,
                "T",
                &point(300.0, 50.0),
                &point(400.0, 50.0),
                ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly)
            )
            .is_err()
    );
    claims
        .query_line_typed(
            &key,
            "T",
            &point(300.0, 50.0),
            &point(400.0, 50.0),
            ClaimSelection::Pinned(id),
        )
        .unwrap();
    let wrong = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(RATE), 300.0)
        .unwrap()
        .with_quantity("shear_rate", QuantitySpec::dimensional(RATE), 50.0)
        .unwrap();
    assert!(
        claims
            .query_line_typed(&key, "T", &wrong, &wrong, ClaimSelection::Pinned(id))
            .is_err()
    );
}

#[test]
fn g0_box_and_convex_hull_do_not_fill_source_support_holes() {
    let mut one_triangle = claim();
    if let PropertyValue::Surface2D { triangles, .. } = &mut one_triangle.value {
        triangles.pop();
    }
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(one_triangle).unwrap();
    for query in [
        point(375.0, 75.0),
        point(350.0, f64::from_bits(50_f64.to_bits() + 1)),
    ] {
        assert!(matches!(
            claims.query_pinned("dynamic_viscosity", &query, id),
            Err(MatDbError::OutsideSurfaceSupport { .. })
        ));
    }
    claims
        .query_pinned("dynamic_viscosity", &point(350.0, 50.0), id)
        .unwrap();
    claims
        .query_pinned("dynamic_viscosity", &point(400.0, 100.0), id)
        .unwrap();
    assert!(
        claims
            .query_envelope(
                "dynamic_viscosity",
                &point(300.0, 0.0),
                &point(400.0, 100.0),
                ClaimSelection::Pinned(id)
            )
            .is_err()
    );

    let mut disconnected = claim();
    if let PropertyValue::Surface2D {
        samples, triangles, ..
    } = &mut disconnected.value
    {
        *samples = [
            [300.0, 0.0],
            [320.0, 0.0],
            [300.0, 100.0],
            [380.0, 0.0],
            [400.0, 0.0],
            [400.0, 100.0],
        ]
        .into_iter()
        .map(|at| (at, law(at[0], at[1])))
        .collect();
        *triangles = vec![[0, 1, 2], [3, 4, 5]];
    }
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(disconnected).unwrap();
    for (t, shear) in [(305.0, 25.0), (395.0, 25.0)] {
        let answer = claims
            .query_pinned("dynamic_viscosity", &point(t, shear), id)
            .unwrap();
        assert_close(answer.evidence.value.value, law(t, shear));
    }
    assert!(matches!(
        claims.query_pinned("dynamic_viscosity", &point(350.0, 25.0), id),
        Err(MatDbError::OutsideSurfaceSupport { .. })
    ));
}

#[test]
fn g0_source_topology_and_finite_enclosing_validity_are_required() {
    let mut malformed = Vec::new();
    for cells in [
        vec![[0, 1, 2], [0, 1, 2]],
        vec![[0, 1, 2], [0, 1, 3]],
        vec![[0, 1, 8]],
        vec![[0, 0, 2]],
    ] {
        let mut source = claim();
        if let PropertyValue::Surface2D { triangles, .. } = &mut source.value {
            *triangles = cells;
        }
        malformed.push(source);
    }
    for at in [[350.0, 0.0], [400.0, 0.0], [f64::NAN, 50.0], [1e100, 50.0]] {
        let mut source = claim();
        if let PropertyValue::Surface2D { samples, .. } = &mut source.value {
            samples[2].0 = at;
        }
        malformed.push(source);
    }
    let mut missing = claim();
    missing.validity = ValidityDomain::unconstrained();
    malformed.push(missing);
    let mut narrow = claim();
    narrow.validity = ValidityDomain::unconstrained()
        .with_quantity("T", temperature(), 300.0, 350.0)
        .with_quantity("shear_rate", QuantitySpec::dimensional(RATE), 0.0, 100.0);
    malformed.push(narrow);
    let mut duplicate_axes = claim();
    if let PropertyValue::Surface2D { axes, .. } = &mut duplicate_axes.value {
        axes[1] = "T".into();
    }
    malformed.push(duplicate_axes);
    for source in malformed {
        assert!(matches!(
            ClaimSet::new().insert_claim(source),
            Err(MatDbError::MalformedSurface { .. })
        ));
    }

    let mut wrong_dims = claim();
    if let PropertyValue::Surface2D { axis_dims, .. } = &mut wrong_dims.value {
        axis_dims[1] = TEMPERATURE;
    }
    assert!(matches!(
        ClaimSet::new().insert_claim(wrong_dims),
        Err(MatDbError::AxisQuantityMismatch { .. })
    ));

    let mut crossing = claim();
    if let PropertyValue::Surface2D {
        samples, triangles, ..
    } = &mut crossing.value
    {
        *samples = [
            [300.0, 25.0],
            [400.0, 25.0],
            [350.0, 100.0],
            [300.0, 75.0],
            [400.0, 75.0],
            [350.0, 0.0],
        ]
        .into_iter()
        .map(|at| (at, law(at[0], at[1])))
        .collect();
        *triangles = vec![[0, 1, 2], [3, 4, 5]];
    }
    assert!(matches!(
        ClaimSet::new().insert_claim(crossing),
        Err(MatDbError::MalformedSurface {
            reason: "surface triangle edges cross"
        })
    ));

    let mut hidden = claim();
    if let PropertyValue::Surface2D { samples, .. } = &mut hidden.value {
        samples.push(([325.0, 25.0], 9.0));
    }
    assert!(matches!(
        ClaimSet::new().insert_claim(hidden),
        Err(MatDbError::MalformedSurface { .. })
    ));
}

#[test]
fn g0_declared_interpolation_and_semantic_axes_cannot_be_bypassed() {
    let mut tabulated = claim();
    tabulated.interpolation = InterpolationPolicy::TabulatedOnly;
    if let PropertyValue::Surface2D { triangles, .. } = &mut tabulated.value {
        triangles.clear();
    }
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(tabulated).unwrap();
    claims
        .query_pinned("dynamic_viscosity", &point(300.0, 0.0), id)
        .unwrap();
    assert!(matches!(
        claims.query_pinned("dynamic_viscosity", &point(350.0, 25.0), id),
        Err(MatDbError::OutsideSurfaceSupport { .. })
    ));
    let difference = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::TemperatureDifference,
        ValueForm::Static,
    ));
    let wrong = QueryPoint::new()
        .with_quantity("T", difference, 300.0)
        .unwrap()
        .with_quantity("shear_rate", QuantitySpec::dimensional(RATE), 0.0)
        .unwrap();
    assert!(
        claims
            .query_pinned("dynamic_viscosity", &wrong, id)
            .is_err()
    );

    let mut constant = claim();
    constant.interpolation = InterpolationPolicy::ConstantWithinValidity;
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(constant).unwrap();
    assert!(matches!(
        claims.query_pinned("dynamic_viscosity", &point(300.0, 0.0), id),
        Err(MatDbError::UnsupportedEvaluation { .. })
    ));
}

#[test]
fn g3_axis_permutation_and_triangle_winding_preserve_the_same_law() {
    let mut first = ClaimSet::new();
    let first_id = first.insert_claim(claim()).unwrap();
    let mut reversed = claim();
    if let PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples,
        triangles,
        ..
    } = &mut reversed.value
    {
        axes.swap(0, 1);
        axis_dims.swap(0, 1);
        for (at, _) in samples {
            at.swap(0, 1);
        }
        for cell in triangles {
            cell.swap(0, 1);
        }
    }
    let mut second = ClaimSet::new();
    let second_id = second.insert_claim(reversed).unwrap();
    assert_ne!(first_id, second_id);
    for (t, shear) in [(325.0, 25.0), (375.0, 75.0), (350.0, 50.0)] {
        let a = first
            .query_pinned("dynamic_viscosity", &point(t, shear), first_id)
            .unwrap();
        let b = second
            .query_pinned("dynamic_viscosity", &point(t, shear), second_id)
            .unwrap();
        assert_close(a.evidence.value.value, b.evidence.value.value);
    }
}

#[test]
fn g0_surface_receipts_bind_source_cell_and_preserve_legacy_versions() {
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(claim()).unwrap();
    let answer = claims
        .query_pinned("dynamic_viscosity", &point(350.0, 25.0), id)
        .unwrap();
    let receipt = answer.receipt;
    assert_eq!(receipt.schema_version, 4);
    assert_eq!(receipt.evaluator_version, 3);
    let bytes = receipt.to_bytes().unwrap();
    assert!(bytes.starts_with(b"FSMATU4\0"));
    let loaded = PropertyUsageReceipt::from_bytes_verified(&bytes, receipt.content_hash()).unwrap();
    assert_eq!(loaded, receipt);
    claims.verify_receipt(&loaded).unwrap();
    let mut wrong_cell = receipt.clone();
    wrong_cell.decision = EvaluationDecision::LinearInsideTriangle { triangle: 1 };
    assert_ne!(wrong_cell.content_hash(), receipt.content_hash());
    assert!(claims.verify_receipt(&wrong_cell).is_err());
    let exact = claims
        .query_pinned("dynamic_viscosity", &point(300.0, 0.0), id)
        .unwrap()
        .receipt;
    let mut wrong_sample = exact.clone();
    wrong_sample.decision = EvaluationDecision::ExactSurfaceSample { sample: 1 };
    assert_ne!(wrong_sample.content_hash(), exact.content_hash());
    assert!(claims.verify_receipt(&wrong_sample).is_err());
    let mut wrong_schema = receipt.clone();
    wrong_schema.schema_version = 3;
    assert!(wrong_schema.to_bytes().is_err());
    let mut wrong_evaluator = receipt;
    wrong_evaluator.evaluator_version = 2;
    assert!(claims.verify_receipt(&wrong_evaluator).is_err());

    let mut untyped = claim();
    untyped.validity = ValidityDomain::unconstrained()
        .with("T", 300.0, 400.0)
        .with("shear_rate", 0.0, 100.0);
    let mut claims = ClaimSet::new();
    let id = claims.insert_claim(untyped).unwrap();
    let query = QueryPoint::new()
        .with("T", 350.0)
        .unwrap()
        .with("shear_rate", 25.0)
        .unwrap();
    let receipt = claims
        .query_pinned("dynamic_viscosity", &query, id)
        .unwrap()
        .receipt;
    assert_eq!(receipt.schema_version, 4);
    assert!(receipt.axis_quantities.is_empty());
    let loaded = PropertyUsageReceipt::from_bytes(&receipt.to_bytes().unwrap()).unwrap();
    claims.verify_receipt(&loaded).unwrap();

    for typed in [false, true] {
        let mut scalar = claim();
        scalar.value = PropertyValue::Scalar {
            value: 1.0,
            dims: VISCOSITY,
        };
        scalar.interpolation = InterpolationPolicy::ConstantWithinValidity;
        let query = if typed {
            point(350.0, 25.0)
        } else {
            scalar.validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
            QueryPoint::new().with("T", 350.0).unwrap()
        };
        let mut claims = ClaimSet::new();
        let id = claims.insert_claim(scalar).unwrap();
        let receipt = claims
            .query_pinned("dynamic_viscosity", &query, id)
            .unwrap()
            .receipt;
        assert_eq!(receipt.schema_version, if typed { 3 } else { 2 });
        assert_eq!(receipt.evaluator_version, 2);
        let loaded = PropertyUsageReceipt::from_bytes(&receipt.to_bytes().unwrap()).unwrap();
        assert_eq!(loaded, receipt);
        claims.verify_receipt(&loaded).unwrap();
    }
}

#[test]
fn g3_surface_pack_retains_samples_cells_uncertainty_and_normalization_targets() {
    let source_hash = source().artifact.unwrap();
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "manufactured two-coordinate support".into(),
            method: "analytic fixture".into(),
            artifact: source_hash,
            caveats: "not experimental evidence".into(),
            provenance: source(),
        })
        .unwrap();
    let mut measured = claim();
    measured.observations.push(observation);
    measured.uncertainty = UncertaintyModel::HalfWidth {
        half_width: 0.1,
        confidence: 0.95,
    };
    let id = claims.insert_claim(measured.clone()).unwrap();
    let members = vec![
        StatisticMember::surface_coordinate(id, 0, 0),
        StatisticMember::surface_coordinate(id, 0, 1),
        StatisticMember::surface_ordinate(id, 0),
    ];
    let statistics = JointStatistics::new(
        observation,
        "source-components",
        members.clone(),
        vec![1.0, 0.0, 4.0, 0.0, 0.0, 0.01],
        None,
    );
    let receipts = members
        .into_iter()
        .zip([TEMPERATURE, RATE, VISCOSITY])
        .map(|(member, dims)| {
            NormalizationReceipt::new(
                NormalizationTarget::ClaimValue(member),
                source_hash,
                dims,
                1.0,
                0.0,
                "declared SI source values",
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            )
        })
        .collect();
    let pack = NormalizedPack::new(
        "surface-fixture",
        "surface-test-v1",
        source_hash,
        "CC0-1.0; redistribution permitted",
        claims,
        vec![statistics],
        receipts,
    )
    .unwrap();
    assert_eq!(pack.schema_version(), MATDB_SURFACE_PACK_SCHEMA_VERSION);
    let bytes = pack.to_bytes();
    let loaded = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(loaded.to_bytes(), bytes);
    assert_eq!(loaded.claims().claim(id), Some(&measured));
    let answer = loaded
        .claims()
        .query_pinned("dynamic_viscosity", &point(375.0, 75.0), id)
        .unwrap();
    assert_close(answer.evidence.value.value, law(375.0, 75.0));
    loaded.claims().verify_receipt(&answer.receipt).unwrap();
    for old_version in [1_u32, 3, 5, 7] {
        let mut relabeled = bytes.clone();
        relabeled[8..12].copy_from_slice(&old_version.to_le_bytes());
        assert!(NormalizedPack::from_bytes(&relabeled).is_err());
    }
}
