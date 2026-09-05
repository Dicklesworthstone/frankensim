//! fs-matdb PR-4 conformance: every answer is Evidence + receipt;
//! extrapolation refuses; fusion is explicit; unstated uncertainty
//! never launders into a certificate.
#![allow(clippy::float_cmp)] // bit-exact propagation of stored values is the contract under test

use fs_blake3::{hash_bytes, hash_domain};
use fs_evidence::NumericalKind;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, EvaluationDecision, InterpolationPolicy, MATDB_EVALUATOR_VERSION,
    MAX_PROPERTY_USAGE_PROPERTY_BYTES, MAX_PROPERTY_USAGE_QUERY_AXES, MatDbError,
    ObservationDataset, PINNED_CLAIM_POLICY_TAG, PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN,
    PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION, PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION, PropertyClaim,
    PropertyKey, PropertyUsageReceipt, PropertyUsageReceiptError, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::Dims;

const DENSITY_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, 0, -1, 0]);

fn provenance(source: &str) -> Provenance {
    Provenance {
        source: source.to_string(),
        license: "internal-use".to_string(),
        artifact: None,
    }
}

fn density(value: f64, source: &str, uncertainty: UncertaintyModel) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("density", DENSITY_DIMS),
        value: PropertyValue::Scalar {
            value,
            dims: DENSITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 400.0),
        uncertainty,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: provenance(source),
    }
}

fn stated() -> UncertaintyModel {
    UncertaintyModel::HalfWidth {
        half_width: 15.0,
        confidence: 0.95,
    }
}

fn room() -> QueryPoint {
    QueryPoint::new().with("T", 293.15).expect("finite point")
}

#[test]
fn g0_typed_axis_queries_and_receipts_preserve_equal_dimension_conventions() {
    use fs_qty::{
        QuantitySpec,
        semantic::{FrequencyConvention as F, QuantityKind as K, SemanticType, ValueForm},
    };
    let spec = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (axis, actual, wrong, coordinate) in [
        (
            "frequency",
            K::Frequency(F::Cyclic),
            K::Frequency(F::Angular),
            50.0,
        ),
        ("T", K::AbsoluteTemperature, K::TemperatureDifference, 300.0),
    ] {
        let mut claim = density(
            2700.0,
            "synthetic typed-axis fixture",
            UncertaintyModel::Unstated,
        );
        claim.validity = ValidityDomain::unconstrained().with_quantity(
            axis,
            spec(actual),
            coordinate,
            coordinate,
        );
        let mut wrong_claim = claim.clone();
        wrong_claim.validity = ValidityDomain::unconstrained().with_quantity(
            axis,
            spec(wrong),
            coordinate,
            coordinate,
        );
        assert_ne!(claim.content_hash(), wrong_claim.content_hash());
        let mut claims = ClaimSet::new();
        let id = claims.insert_claim(claim).unwrap();
        let point = QueryPoint::new()
            .with_quantity(axis, spec(actual), coordinate)
            .unwrap();
        let answer = claims
            .query("density", &point, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert_eq!(point.clone().with(axis, coordinate).unwrap(), point);
        assert_eq!(answer.evidence.value.value, 2700.0);
        assert_eq!(answer.receipt.schema_version, 3);
        let bytes = answer.receipt.to_bytes().unwrap();
        let receipt =
            PropertyUsageReceipt::from_bytes_verified(&bytes, answer.receipt.content_hash())
                .unwrap();
        assert_eq!(receipt.to_bytes().unwrap(), bytes);
        claims.verify_receipt(&receipt).unwrap();
        claims
            .verify_receipt(&claims.query_pinned("density", &point, id).unwrap().receipt)
            .unwrap();
        for wrong_point in [
            QueryPoint::new().with(axis, coordinate).unwrap(),
            QueryPoint::new()
                .with_quantity(axis, spec(wrong), coordinate)
                .unwrap(),
            QueryPoint::new()
                .with_quantity(
                    axis,
                    QuantitySpec::dimensional(actual.expected_dims()),
                    coordinate,
                )
                .unwrap(),
        ] {
            assert!(matches!(
                claims.query("density", &wrong_point, SelectionPolicy::SingleClaimOnly),
                Err(MatDbError::AxisQuantityMismatch { .. })
            ));
            assert!(claims.query_pinned("density", &wrong_point, id).is_err());
        }
        let mut tampered = receipt.clone();
        tampered.axis_quantities.insert(axis.into(), spec(wrong));
        assert_ne!(tampered.content_hash(), receipt.content_hash());
        assert!(claims.verify_receipt(&tampered).is_err());
        tampered = receipt.clone();
        tampered.schema_version = 2;
        assert!(tampered.to_bytes().is_err());
        tampered = receipt.clone();
        tampered.axis_quantities.clear();
        assert!(tampered.to_bytes().is_err());
        let outside = QueryPoint::new()
            .with_quantity(axis, spec(actual), coordinate + 1.0)
            .unwrap();
        assert!(matches!(
            claims.query("density", &outside, SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
}

#[test]
fn g0_material_quantity_queries_refuse_equal_dimension_convention_mismatches() {
    use fs_qty::{
        QuantitySpec,
        semantic::{
            AttenuationConvention as A, FrequencyConvention as F, HardnessScale as H,
            MoistureBasis as M, QuantityKind as K, SemanticType, StrainBasis, StrainComponent,
            ValueForm,
        },
    };
    let pairs = [
        (K::AbsoluteTemperature, K::TemperatureDifference),
        (K::Frequency(F::Cyclic), K::Frequency(F::Angular)),
        (K::Energy, K::Torque),
        (K::Hardness(H::RockwellC), K::Hardness(H::Vickers)),
        (K::Moisture(M::DryMass), K::Moisture(M::WetMass)),
        (K::Moisture(M::RelativeHumidity), K::Moisture(M::WetMass)),
        (
            K::Attenuation(A::AmplitudeNepersPerMetre),
            K::Attenuation(A::PowerNepersPerMetre),
        ),
        (
            K::Attenuation(A::AmplitudeNepersPerMetre),
            K::Attenuation(A::DecibelsPerMetre),
        ),
        (
            K::Strain {
                basis: StrainBasis::Tensor,
                component: StrainComponent::Shear,
            },
            K::Strain {
                basis: StrainBasis::Engineering,
                component: StrainComponent::Shear,
            },
        ),
        (K::OpticalExtinctionCoefficient, K::ThermalConductivity),
    ];
    for (actual, wrong) in pairs {
        let spec = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let key = PropertyKey::with_quantity("quantity-fixture", spec(actual));
        let mut claim = density(
            1.0,
            "explicit synthetic convention test",
            UncertaintyModel::Unstated,
        );
        claim.key = key.clone();
        claim.value = PropertyValue::Scalar {
            value: 1.0,
            dims: actual.expected_dims(),
        };
        let mut set = ClaimSet::new();
        let id = set.insert_claim(claim.clone()).unwrap();
        let wrong_key = PropertyKey::with_quantity(key.name(), spec(wrong));
        assert!(matches!(
            set.query_typed(&wrong_key, &room(), SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::QuantityMismatch { .. })
        ));
        assert!(matches!(
            set.query_pinned_typed(&wrong_key, &room(), id),
            Err(MatDbError::QuantityMismatch { .. })
        ));
        assert!(matches!(
            set.query(key.name(), &room(), SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::QuantityMismatch { .. })
        ));
        assert!(matches!(
            set.query_pinned(key.name(), &room(), id),
            Err(MatDbError::QuantityMismatch { .. })
        ));
        let answer = set
            .query_typed(&key, &room(), SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert_eq!(answer.evidence.value.quantity, spec(actual));
        assert_eq!(answer.evidence.value.value, 1.0);
        let decoded =
            PropertyUsageReceipt::from_bytes(&answer.receipt.to_bytes().unwrap()).unwrap();
        set.verify_receipt(&decoded).unwrap();
        assert_eq!(decoded.selected, id);

        let mut conflicting = claim.clone();
        conflicting.key = wrong_key;
        conflicting.value = PropertyValue::Scalar {
            value: 1.0,
            dims: wrong.expected_dims(),
        };
        assert!(set.insert_claim(conflicting.clone()).is_err());
        assert_eq!(set.claim_count(), 1, "refusal is transactional");
        let mut other = ClaimSet::new();
        let other_id = other.insert_claim(conflicting).unwrap();
        assert_ne!(id, other_id, "kind changes claim and receipt identity");
        assert!(other.verify_receipt(&decoded).is_err());

        claim.key = PropertyKey::new(key.name(), actual.expected_dims());
        let mut legacy = ClaimSet::new();
        legacy.insert_claim(claim).unwrap();
        assert!(matches!(
            legacy.query_typed(&key, &room(), SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::QuantityMismatch { .. })
        ));
    }
}

#[test]
fn g0_typed_property_admission_uses_shared_scalar_domain_and_form_rules() {
    use fs_qty::{
        QuantitySpec,
        semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
    };
    for (kind, form, value) in [
        (QuantityKind::AbsoluteTemperature, ValueForm::Static, -1.0),
        (QuantityKind::AbsoluteTemperature, ValueForm::Rms, 300.0),
        (
            QuantityKind::Moisture(MoistureBasis::WetMass),
            ValueForm::Static,
            1.5,
        ),
    ] {
        let mut claim = density(value, "synthetic invalid value", UncertaintyModel::Unstated);
        claim.key = PropertyKey::with_quantity(
            "typed",
            QuantitySpec::semantic(SemanticType::new(kind, form)),
        );
        claim.value = PropertyValue::Scalar {
            value,
            dims: kind.expected_dims(),
        };
        let mut set = ClaimSet::new();
        assert!(matches!(
            set.insert_claim(claim),
            Err(MatDbError::InvalidQuantityValue { .. })
        ));
        assert_eq!(set.claim_count(), 0);
    }
}

fn portable_receipt_fixture() -> (ClaimSet, PropertyUsageReceipt) {
    let mut set = ClaimSet::new();
    let observation = set
        .register_observation(ObservationDataset {
            specimen: "AA6061-T6 plate".to_string(),
            method: "ASTM B311".to_string(),
            artifact: hash_bytes(b"portable receipt observation"),
            caveats: "none".to_string(),
            provenance: provenance("portable receipt lab report"),
        })
        .expect("observation registers");
    let mut claim = density(2699.0, "portable receipt source", stated());
    claim.observations.push(observation);
    set.insert_claim(claim).expect("claim inserts");
    let receipt = set
        .query("density", &room(), SelectionPolicy::SingleClaimOnly)
        .expect("fixture query answers")
        .receipt;
    (set, receipt)
}

fn assert_receipt_identity_moves(
    baseline: &PropertyUsageReceipt,
    changed: &PropertyUsageReceipt,
    field: &str,
) {
    assert_ne!(
        baseline.content_hash(),
        changed.content_hash(),
        "receipt identity must bind {field}"
    );
}

fn push_receipt_wire_string(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

#[test]
fn answers_carry_honest_evidence_and_complete_receipts() {
    let mut set = ClaimSet::new();
    let obs = set
        .register_observation(ObservationDataset {
            specimen: "AA6061-T6 plate".to_string(),
            method: "ASTM B311".to_string(),
            artifact: hash_bytes(b"raw table"),
            caveats: "none".to_string(),
            provenance: provenance("lab report 9"),
        })
        .expect("observation registers");
    let mut claim = density(2700.0, "MMPDS", stated());
    claim.observations.push(obs);
    let id = set.insert_claim(claim).expect("claim inserts");

    let answer = set
        .query("density", &room(), SelectionPolicy::SingleClaimOnly)
        .expect("in-domain query answers");
    assert_eq!(answer.evidence.value.value, 2700.0);
    assert_eq!(answer.evidence.qoi, 2700.0);
    assert_eq!(answer.evidence.numerical.kind, NumericalKind::Estimate);
    assert_eq!(answer.evidence.numerical.lo, 2685.0);
    assert_eq!(answer.evidence.numerical.hi, 2715.0);
    assert!(answer.evidence.model.in_domain);
    assert_eq!(
        answer.evidence.model.validity.bound("T"),
        Some((250.0, 400.0))
    );

    let receipt = &answer.receipt;
    assert_eq!(
        receipt.schema_version,
        PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION
    );
    assert_eq!(receipt.property, "density");
    assert_eq!(receipt.query_point, vec![("T".to_string(), 293.15)]);
    assert_eq!(receipt.considered, vec![id]);
    assert_eq!(receipt.in_domain, vec![id]);
    assert_eq!(receipt.selected, id);
    assert_eq!(receipt.policy, "single-claim-only");
    assert_eq!(receipt.decision, EvaluationDecision::ConstantWithinValidity);
    assert!(receipt.observation_backed);
    assert_eq!(receipt.evaluator_version, MATDB_EVALUATOR_VERSION);
    assert_eq!(receipt.source_hashes.len(), 2, "claim + one observation");

    let other_point = QueryPoint::new().with("T", 300.0).expect("finite");
    let other = set
        .query("density", &other_point, SelectionPolicy::SingleClaimOnly)
        .expect("second query");
    assert_ne!(
        receipt.content_hash(),
        other.receipt.content_hash(),
        "the receipt identity binds the query point"
    );
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"query-receipt\",\"verdict\":\"pass\",\
         \"detail\":\"evidence slices honest; receipt complete and point-sensitive\"}}"
    );
}

#[test]
fn extrapolation_and_unknown_property_refuse() {
    let mut set = ClaimSet::new();
    set.insert_claim(density(2700.0, "MMPDS", stated()))
        .expect("claim inserts");

    let cold = QueryPoint::new().with("T", 150.0).expect("finite");
    assert!(matches!(
        set.query("density", &cold, SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::NoClaimInDomain { considered: 1, .. })
    ));
    assert!(matches!(
        set.query("viscosity", &room(), SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::UnknownProperty { .. })
    ));
    assert!(matches!(
        QueryPoint::new().with("T", f64::INFINITY),
        Err(MatDbError::NonFiniteQueryPoint { .. })
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"extrapolation-refusal\",\"verdict\":\"pass\",\
         \"detail\":\"out-of-validity, unknown property, and non-finite points refuse typed\"}}"
    );
}

#[test]
fn fusion_is_explicit_and_ambiguity_refuses() {
    let mut set = ClaimSet::new();
    set.insert_claim(density(2700.0, "MMPDS", stated()))
        .expect("first claim");
    let obs = set
        .register_observation(ObservationDataset {
            specimen: "AA6061-T6 bar".to_string(),
            method: "ASTM B311".to_string(),
            artifact: hash_bytes(b"bar table"),
            caveats: "none".to_string(),
            provenance: provenance("lab report 12"),
        })
        .expect("observation registers");
    let mut backed = density(2698.5, "internal lab", stated());
    backed.observations.push(obs);
    let backed_id = set.insert_claim(backed).expect("second claim");

    assert!(matches!(
        set.query("density", &room(), SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::AmbiguousSelection { candidates, .. }) if candidates.len() == 2
    ));

    let preferred = set
        .query("density", &room(), SelectionPolicy::PreferObservationBacked)
        .expect("observation-backed claim wins");
    assert_eq!(preferred.receipt.selected, backed_id);
    assert!(preferred.receipt.observation_backed);
    assert_eq!(preferred.receipt.in_domain.len(), 2);
    assert_eq!(preferred.evidence.value.value, 2698.5);
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"explicit-fusion\",\"verdict\":\"pass\",\
         \"detail\":\"ambiguity refuses under single-claim; observation-backed preference is a \
         named policy in the receipt\"}}"
    );
}

#[test]
fn curves_interpolate_inside_and_refuse_beyond_data() {
    let mut set = ClaimSet::new();
    set.insert_claim(PropertyClaim {
        key: PropertyKey::new("electrical-conductivity", CONDUCTIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".to_string(),
            abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
            knots: vec![(256.0, 3.8e7), (320.0, 3.4e7)],
            dims: CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 400.0),
        uncertainty: stated(),
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: provenance("handbook"),
    })
    .expect("curve inserts");

    // 288 is the exact midpoint of the [256, 320] span in binary, so
    // the interpolated value is bit-exact.
    let mid = QueryPoint::new().with("T", 288.0).expect("finite");
    let answer = set
        .query(
            "electrical-conductivity",
            &mid,
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("interpolates inside");
    assert_eq!(answer.evidence.value.value, 3.6e7);
    assert_eq!(
        answer.receipt.decision,
        EvaluationDecision::LinearInside {
            x_lo: 256.0,
            x_hi: 320.0
        }
    );

    let knot = QueryPoint::new().with("T", 256.0).expect("finite");
    let hit = set
        .query(
            "electrical-conductivity",
            &knot,
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("exact knot answers");
    assert_eq!(
        hit.receipt.decision,
        EvaluationDecision::ExactTabulated { at: 256.0 }
    );

    // Inside VALIDITY but beyond the knot span: data ends, so the
    // answer refuses rather than extrapolating the last segment.
    let beyond = QueryPoint::new().with("T", 380.0).expect("finite");
    assert!(matches!(
        set.query(
            "electrical-conductivity",
            &beyond,
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::OutsideKnotSpan { .. })
    ));

    // An empty point is NOT contained by a T-constrained validity, so
    // the validity gate refuses FIRST (fail-closed ordering).
    let axisless = QueryPoint::new();
    assert!(matches!(
        set.query(
            "electrical-conductivity",
            &axisless,
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::NoClaimInDomain { .. })
    ));

    // MissingQueryAxis is reachable only through an UNCONSTRAINED
    // validity: the claim admits any point, but the curve still needs
    // its abscissa coordinate.
    set.insert_claim(PropertyClaim {
        key: PropertyKey::new("thermal-conductivity", Dims([1, 1, -3, -1, 0, 0])),
        value: PropertyValue::Curve {
            abscissa: "T".to_string(),
            abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
            knots: vec![(256.0, 200.0), (320.0, 180.0)],
            dims: Dims([1, 1, -3, -1, 0, 0]),
        },
        validity: ValidityDomain::unconstrained(),
        uncertainty: stated(),
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: provenance("handbook"),
    })
    .expect("unconstrained curve inserts");
    assert!(matches!(
        set.query(
            "thermal-conductivity",
            &QueryPoint::new(),
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::MissingQueryAxis { .. })
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"curve-evaluation\",\"verdict\":\"pass\",\
         \"detail\":\"linear inside knots, exact hits tagged, beyond-data and axis-less refuse\"}}"
    );
}

#[test]
fn receipt_replay_compares_decision_floats_by_exact_bits() {
    let mut set = ClaimSet::new();
    set.insert_claim(PropertyClaim {
        key: PropertyKey::new("signed-zero-probe", CONDUCTIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: "phase".to_string(),
            abscissa_dims: Dims([0, 0, 0, 0, 0, 0]),
            knots: vec![(-0.0, 1.0), (1.0, 2.0)],
            dims: CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("phase", -1.0, 1.0),
        uncertainty: stated(),
        interpolation: InterpolationPolicy::TabulatedOnly,
        observations: Vec::new(),
        provenance: provenance("signed-zero fixture"),
    })
    .expect("signed-zero curve inserts");
    let point = QueryPoint::new()
        .with("phase", -0.0)
        .expect("signed zero is finite");
    let authentic = set
        .query(
            "signed-zero-probe",
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("exact signed-zero knot answers")
        .receipt;
    assert!(matches!(
        &authentic.decision,
        EvaluationDecision::ExactTabulated { at } if at.to_bits() == (-0.0_f64).to_bits()
    ));

    let mut sign_tampered = authentic.clone();
    sign_tampered.decision = EvaluationDecision::ExactTabulated { at: 0.0 };
    sign_tampered
        .try_content_hash()
        .expect("the tampered shape is still portable");
    assert_ne!(
        authentic.content_hash(),
        sign_tampered.content_hash(),
        "exact-bit receipt identity must bind the zero sign"
    );
    assert!(matches!(
        set.verify_receipt(&sign_tampered),
        Err(MatDbError::ReceiptMismatch { field: "decision" })
    ));
}

#[test]
fn unstated_uncertainty_is_marked_and_never_certifies() {
    let mut set = ClaimSet::new();
    set.insert_claim(density(
        2700.0,
        "vendor datasheet",
        UncertaintyModel::Unstated,
    ))
    .expect("unstated claim inserts");
    let answer = set
        .query("density", &room(), SelectionPolicy::SingleClaimOnly)
        .expect("unstated claims still answer");
    assert_eq!(answer.evidence.numerical.kind, NumericalKind::NoClaim);
    assert!(
        answer.evidence.clone().certified().is_err(),
        "an unstated-uncertainty answer must never certify"
    );
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"no-laundering\",\"verdict\":\"pass\",\
         \"detail\":\"Unstated maps to an explicit numerical no-claim and certification refuses\"}}"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // one mutation per receipt field, in field order
fn receipt_completeness_mutation_battery() {
    // PR-5: a receipt with ANY deleted, substituted, or stale field
    // fails verification with a typed refusal. Fixture: two claims so
    // considered/in_domain/selected are all nontrivial.
    let mut set = ClaimSet::new();
    set.insert_claim(density(2700.0, "MMPDS", stated()))
        .expect("citation claim");
    let obs = set
        .register_observation(ObservationDataset {
            specimen: "AA6061-T6 plate".to_string(),
            method: "ASTM B311".to_string(),
            artifact: hash_bytes(b"plate table"),
            caveats: "none".to_string(),
            provenance: provenance("lab report 9"),
        })
        .expect("observation registers");
    let mut backed = density(2698.5, "internal lab", stated());
    backed.observations.push(obs);
    let other_id = set.insert_claim(backed).expect("backed claim");
    let citation_id = set.claims_for("density")[0].0;

    let answer = set
        .query("density", &room(), SelectionPolicy::PreferObservationBacked)
        .expect("query answers");
    let good = answer.receipt.clone();
    set.verify_receipt(&good)
        .expect("authentic receipt verifies");

    let mutations: Vec<(&str, PropertyUsageReceipt)> = vec![
        ("property", {
            let mut r = good.clone();
            r.property = "viscosity".to_string();
            r
        }),
        ("query_point", {
            let mut r = good.clone();
            r.query_point = vec![("T".to_string(), 150.0)];
            r
        }),
        ("considered", {
            let mut r = good.clone();
            r.considered = vec![other_id];
            r
        }),
        ("in_domain", {
            let mut r = good.clone();
            r.in_domain = vec![citation_id];
            r
        }),
        ("selected", {
            let mut r = good.clone();
            r.selected = citation_id;
            r
        }),
        ("policy", {
            let mut r = good.clone();
            r.policy = "single-claim-only";
            r
        }),
        ("foreign-policy", {
            let mut r = good.clone();
            r.policy = "trust-me";
            r
        }),
        ("decision", {
            let mut r = good.clone();
            r.decision = EvaluationDecision::ExactScalar;
            r
        }),
        ("observation_backed", {
            let mut r = good.clone();
            r.observation_backed = false;
            r
        }),
        ("evaluator_version", {
            let mut r = good.clone();
            r.evaluator_version = 999;
            r
        }),
        ("source_hashes", {
            let mut r = good.clone();
            r.source_hashes.pop();
            r
        }),
        ("schema_version", {
            let mut r = good.clone();
            r.schema_version += 1;
            r
        }),
    ];
    for (label, mutated) in &mutations {
        let refused = set.verify_receipt(mutated);
        assert!(
            refused.is_err(),
            "mutated receipt field '{label}' must fail verification"
        );
        assert_ne!(
            mutated.content_hash(),
            good.content_hash(),
            "mutated receipt field '{label}' must move the receipt identity"
        );
    }
    // The refusals are TYPED per class, not one blanket error.
    assert!(matches!(
        set.verify_receipt(&mutations[0].1),
        Err(MatDbError::UnknownProperty { .. })
    ));
    assert!(matches!(
        set.verify_receipt(&mutations[1].1),
        Err(MatDbError::NoClaimInDomain { .. })
    ));
    assert!(matches!(
        set.verify_receipt(&mutations[4].1),
        Err(MatDbError::PropertyUsageReceiptNotPortable {
            error: PropertyUsageReceiptError::InvalidField {
                field: "source-hashes",
                ..
            }
        })
    ));
    assert!(matches!(
        set.verify_receipt(&mutations[5].1),
        Err(MatDbError::AmbiguousSelection { .. })
    ));
    assert!(matches!(
        set.verify_receipt(&mutations[6].1),
        Err(MatDbError::UnknownPolicyTag { .. })
    ));
    assert!(matches!(
        set.verify_receipt(&mutations[7].1),
        Err(MatDbError::ReceiptMismatch { field: "decision" })
    ));
    assert!(matches!(
        set.verify_receipt(&mutations[9].1),
        Err(MatDbError::EvaluatorVersionDrift {
            receipt: 999,
            current: 1
        })
    ));
    assert!(matches!(
        set.verify_receipt(&mutations[11].1),
        Err(MatDbError::ReceiptSchemaVersionDrift {
            receipt,
            current: PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION,
        }) if receipt == PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION + 1
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"receipt-mutation-battery\",\"verdict\":\"pass\",\
         \"detail\":\"12 field mutations all refuse typed and all move the receipt identity\"}}"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // the sealed-registry version guard walks every field once
fn property_usage_receipt_v2_round_trips_and_replays_exactly() {
    let (set, receipt) = portable_receipt_fixture();
    assert_eq!(PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION, 2);
    assert_eq!(
        PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN,
        "org.frankensim.fs-matdb.property-usage-receipt.v2"
    );
    assert_eq!(
        receipt.try_content_hash().expect("portable identity"),
        receipt.content_hash()
    );

    let bytes = receipt.to_bytes().expect("receipt encodes");
    let mut manual = Vec::new();
    manual.extend_from_slice(b"FSMATUR\0");
    manual.extend_from_slice(&PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION.to_le_bytes());
    manual.push(1);
    push_receipt_wire_string(&mut manual, &receipt.property);
    manual.push(2);
    manual.extend_from_slice(&(receipt.query_point.len() as u64).to_le_bytes());
    for (axis, value) in &receipt.query_point {
        push_receipt_wire_string(&mut manual, axis);
        manual.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    manual.push(3);
    manual.extend_from_slice(&(receipt.considered.len() as u64).to_le_bytes());
    for id in &receipt.considered {
        manual.extend_from_slice(id.0.as_bytes());
    }
    manual.push(4);
    manual.extend_from_slice(&(receipt.in_domain.len() as u64).to_le_bytes());
    for id in &receipt.in_domain {
        manual.extend_from_slice(id.0.as_bytes());
    }
    manual.push(5);
    manual.extend_from_slice(receipt.selected.0.as_bytes());
    manual.push(6);
    push_receipt_wire_string(&mut manual, receipt.policy);
    manual.extend_from_slice(&[7, 1, 8, 1, 9]);
    manual.extend_from_slice(&receipt.evaluator_version.to_le_bytes());
    manual.push(10);
    manual.extend_from_slice(&(receipt.source_hashes.len() as u64).to_le_bytes());
    for hash in &receipt.source_hashes {
        manual.extend_from_slice(hash.as_bytes());
    }
    assert_eq!(
        &bytes[..bytes.len() - 32],
        manual,
        "independent manual preimage locks magic, tags, widths, order, counts, and exact bits"
    );
    assert_eq!(
        &bytes[bytes.len() - 32..],
        receipt.content_hash().as_bytes(),
        "transport retains the exact v2 content identity"
    );
    let mut identity_preimage = PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION
        .to_le_bytes()
        .to_vec();
    identity_preimage.extend_from_slice(&manual);
    assert_eq!(
        receipt.content_hash(),
        hash_domain(PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN, &identity_preimage),
        "the independent identity preimage binds version plus the exact wire field stream"
    );
    let mut alternate_version = identity_preimage.clone();
    alternate_version[..4]
        .copy_from_slice(&(PROPERTY_USAGE_RECEIPT_IDENTITY_VERSION + 1).to_le_bytes());
    assert_ne!(
        receipt.content_hash(),
        hash_domain(PROPERTY_USAGE_RECEIPT_IDENTITY_DOMAIN, &alternate_version),
        "identity-version rotation must move the digest"
    );
    assert_ne!(
        receipt.content_hash(),
        hash_domain(
            "org.frankensim.fs-matdb.property-usage-receipt.foreign",
            &identity_preimage,
        ),
        "identity-domain rotation must move the digest"
    );
    let decoded = PropertyUsageReceipt::from_bytes(&bytes).expect("receipt decodes");
    assert_eq!(decoded, receipt);
    assert_eq!(decoded.to_bytes().expect("fixed point"), bytes);
    set.verify_receipt(&decoded)
        .expect("decoded receipt replays against its claim set");

    let pinned = PropertyUsageReceipt::from_bytes_verified(&bytes, receipt.content_hash())
        .expect("externally pinned identity admits");
    assert_eq!(pinned, receipt);

    let mut stale = receipt.clone();
    stale.evaluator_version += 1;
    let stale = PropertyUsageReceipt::from_bytes(
        &stale
            .to_bytes()
            .expect("structurally valid stale evaluator receipt encodes"),
    )
    .expect("stale evaluator receipt remains structurally decodable");
    assert!(matches!(
        set.verify_receipt(&stale),
        Err(MatDbError::EvaluatorVersionDrift {
            receipt,
            current: MATDB_EVALUATOR_VERSION,
        }) if receipt == MATDB_EVALUATOR_VERSION + 1
    ));

    for decision in [
        EvaluationDecision::ConstantWithinValidity,
        EvaluationDecision::ExactScalar,
        EvaluationDecision::ExactTabulated { at: -0.0 },
        EvaluationDecision::LinearInside {
            x_lo: -0.0,
            x_hi: 1.0,
        },
    ] {
        let mut variant = receipt.clone();
        variant.decision = decision;
        let encoded = variant.to_bytes().expect("decision variant encodes");
        let recovered = PropertyUsageReceipt::from_bytes(&encoded).expect("variant decodes");
        match (&variant.decision, &recovered.decision) {
            (
                EvaluationDecision::ExactTabulated { at: expected },
                EvaluationDecision::ExactTabulated { at: actual },
            ) => assert_eq!(actual.to_bits(), expected.to_bits()),
            (
                EvaluationDecision::LinearInside {
                    x_lo: expected_lo,
                    x_hi: expected_hi,
                },
                EvaluationDecision::LinearInside {
                    x_lo: actual_lo,
                    x_hi: actual_hi,
                },
            ) => {
                assert_eq!(actual_lo.to_bits(), expected_lo.to_bits());
                assert_eq!(actual_hi.to_bits(), expected_hi.to_bits());
            }
            (expected, actual) => assert_eq!(actual, expected),
        }
    }
}

#[test]
fn property_usage_receipt_v2_binds_collection_boundaries() {
    let (_, receipt) = portable_receipt_fixture();
    let a = fs_matdb::ClaimId(hash_bytes(b"boundary-a"));
    let b = fs_matdb::ClaimId(hash_bytes(b"boundary-b"));

    let mut left = receipt.clone();
    left.considered = vec![a, b];
    left.in_domain = vec![b];
    left.selected = b;
    left.source_hashes[0] = b.0;

    let mut right = left.clone();
    right.considered = vec![a];
    right.in_domain = vec![b, b];

    let legacy_unframed_ids = |candidate: &PropertyUsageReceipt| {
        candidate
            .considered
            .iter()
            .chain(&candidate.in_domain)
            .flat_map(|id| id.0.0)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        legacy_unframed_ids(&left),
        legacy_unframed_ids(&right),
        "v1 could not recover this collection boundary"
    );
    assert_ne!(
        left.content_hash(),
        right.content_hash(),
        "v2 count framing must bind the boundary"
    );
    left.try_content_hash()
        .expect("the query-shaped side remains portable");
    assert!(matches!(
        right.try_content_hash(),
        Err(PropertyUsageReceiptError::InvalidField {
            field: "considered" | "in-domain",
            ..
        })
    ));
}

#[test]
fn property_usage_receipt_identity_fields_move_independently() {
    let (_, baseline) = portable_receipt_fixture();

    let mut changed = baseline.clone();
    changed.schema_version += 1;
    assert_receipt_identity_moves(&baseline, &changed, "wire schema version");

    let mut changed = baseline.clone();
    changed.property.push_str("-alternate");
    assert_receipt_identity_moves(&baseline, &changed, "property bytes and length");

    let mut changed = baseline.clone();
    changed.query_point[0].0 = "temperature".to_string();
    assert_receipt_identity_moves(&baseline, &changed, "query axis bytes and length");

    let mut changed = baseline.clone();
    changed.query_point[0].1 = f64::from_bits(changed.query_point[0].1.to_bits() + 1);
    assert_receipt_identity_moves(&baseline, &changed, "query coordinate exact bits");

    let mut changed = baseline.clone();
    changed
        .query_point
        .push(("pressure".to_string(), 101_325.0));
    assert_receipt_identity_moves(&baseline, &changed, "query-point count and order");

    let alternate = fs_matdb::ClaimId(hash_bytes(b"alternate claim identity"));
    let mut changed = baseline.clone();
    changed.considered.push(alternate);
    assert_receipt_identity_moves(&baseline, &changed, "considered count/order/id");

    let mut changed = baseline.clone();
    changed.considered.push(alternate);
    changed.in_domain.push(alternate);
    assert_receipt_identity_moves(&baseline, &changed, "in-domain count/order/id");

    let mut changed = baseline.clone();
    changed.selected = alternate;
    changed.considered.push(alternate);
    changed.in_domain.push(alternate);
    changed.source_hashes[0] = alternate.0;
    assert_receipt_identity_moves(&baseline, &changed, "selected claim id");

    let mut changed = baseline.clone();
    changed.policy = "prefer-observation-backed";
    assert_receipt_identity_moves(&baseline, &changed, "policy bytes and length");

    let mut changed = baseline.clone();
    changed.decision = EvaluationDecision::ExactScalar;
    assert_receipt_identity_moves(&baseline, &changed, "decision tag");

    let mut exact_a = baseline.clone();
    exact_a.decision = EvaluationDecision::ExactTabulated { at: 1.0 };
    let mut exact_b = exact_a.clone();
    exact_b.decision = EvaluationDecision::ExactTabulated {
        at: f64::from_bits(1.0_f64.to_bits() + 1),
    };
    assert_receipt_identity_moves(&exact_a, &exact_b, "exact-tabulated at bits");

    let mut linear_a = baseline.clone();
    linear_a.decision = EvaluationDecision::LinearInside {
        x_lo: 1.0,
        x_hi: 3.0,
    };
    let mut linear_b = linear_a.clone();
    linear_b.decision = EvaluationDecision::LinearInside {
        x_lo: 2.0,
        x_hi: 3.0,
    };
    assert_receipt_identity_moves(&linear_a, &linear_b, "linear x-lo bits");
    linear_b.decision = EvaluationDecision::LinearInside {
        x_lo: 1.0,
        x_hi: 4.0,
    };
    assert_receipt_identity_moves(&linear_a, &linear_b, "linear x-hi bits");

    let mut changed = baseline.clone();
    changed.observation_backed = false;
    assert_receipt_identity_moves(&baseline, &changed, "observation-backed flag");

    let mut changed = baseline.clone();
    changed.evaluator_version += 1;
    assert_receipt_identity_moves(&baseline, &changed, "evaluator version");

    let mut changed = baseline.clone();
    changed.source_hashes[1] = hash_bytes(b"alternate observation source");
    assert_receipt_identity_moves(&baseline, &changed, "source hash");

    let mut changed = baseline.clone();
    changed
        .source_hashes
        .push(hash_bytes(b"additional observation source"));
    assert_receipt_identity_moves(&baseline, &changed, "source count and order");
}

#[test]
#[allow(clippy::too_many_lines)] // one closed-decoder drill per malformation class
fn property_usage_receipt_v2_decoder_fails_closed() {
    let (_, receipt) = portable_receipt_fixture();
    let bytes = receipt.to_bytes().expect("fixture encodes");

    for length in 0..bytes.len() {
        assert!(
            PropertyUsageReceipt::from_bytes(&bytes[..length]).is_err(),
            "truncated prefix of {length} bytes must refuse"
        );
    }

    let mut wrong_schema = bytes.clone();
    wrong_schema[8..12].copy_from_slice(&1_u32.to_le_bytes());
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&wrong_schema),
        Err(PropertyUsageReceiptError::UnsupportedSchemaVersion {
            found: 1,
            supported: PROPERTY_USAGE_RECEIPT_SCHEMA_VERSION,
        })
    ));

    let mut wrong_magic = bytes.clone();
    wrong_magic[0] ^= 1;
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&wrong_magic),
        Err(PropertyUsageReceiptError::Malformed { .. })
    ));

    let mut wrong_field_tag = bytes.clone();
    wrong_field_tag[12] = u8::MAX;
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&wrong_field_tag),
        Err(PropertyUsageReceiptError::Malformed { .. })
    ));

    let mut excessive_axes = bytes.clone();
    let query_count_offset = 8 + 4 + 1 + 8 + receipt.property.len() + 1;
    excessive_axes[query_count_offset..query_count_offset + 8]
        .copy_from_slice(&((MAX_PROPERTY_USAGE_QUERY_AXES as u64) + 1).to_le_bytes());
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&excessive_axes),
        Err(PropertyUsageReceiptError::ResourceLimit {
            resource: "query-axes",
            ..
        })
    ));

    let mut hostile_property_length = bytes.clone();
    hostile_property_length[13..21].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&hostile_property_length),
        Err(PropertyUsageReceiptError::ResourceLimit {
            resource: "property-bytes",
            observed: u64::MAX,
            ..
        })
    ));

    let source_count_offset = bytes.len() - 32 - 32 * receipt.source_hashes.len() - 8;
    let mut hostile_source_count = bytes.clone();
    hostile_source_count[source_count_offset..source_count_offset + 8]
        .copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&hostile_source_count),
        Err(PropertyUsageReceiptError::ResourceLimit {
            resource: "source-hashes",
            observed: u64::MAX,
            ..
        })
    ));

    let policy_offset = bytes
        .windows(receipt.policy.len())
        .position(|window| window == receipt.policy.as_bytes())
        .expect("policy bytes are present");
    let mut unknown_policy = bytes.clone();
    unknown_policy[policy_offset] = b'x';
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&unknown_policy),
        Err(PropertyUsageReceiptError::UnknownPolicyTag { .. })
    ));

    let source_section_len = 1 + 8 + 32 * receipt.source_hashes.len();
    let observation_tag_offset = bytes.len() - 32 - source_section_len - (1 + 4) - (1 + 1);
    let decision_variant_offset = observation_tag_offset - 1;
    let mut unknown_decision = bytes.clone();
    unknown_decision[decision_variant_offset] = u8::MAX;
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&unknown_decision),
        Err(PropertyUsageReceiptError::UnknownDecisionTag { tag: u8::MAX, .. })
    ));

    let mut invalid_boolean = bytes.clone();
    invalid_boolean[observation_tag_offset + 1] = 2;
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&invalid_boolean),
        Err(PropertyUsageReceiptError::Malformed { .. })
    ));

    let mut identity_tamper = bytes.clone();
    *identity_tamper.last_mut().expect("identity byte") ^= 1;
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&identity_tamper),
        Err(PropertyUsageReceiptError::IdentityMismatch { .. })
    ));

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(matches!(
        PropertyUsageReceipt::from_bytes(&trailing),
        Err(PropertyUsageReceiptError::Malformed { .. })
    ));

    assert!(matches!(
        PropertyUsageReceipt::from_bytes_verified(&bytes, hash_bytes(b"wrong external identity")),
        Err(PropertyUsageReceiptError::ExternalIdentityMismatch { .. })
    ));
}

#[test]
fn property_usage_receipt_v2_encoder_enforces_canonical_caps_and_relations() {
    let (set, receipt) = portable_receipt_fixture();

    let mut maximum_property = receipt.clone();
    maximum_property.property = "p".repeat(MAX_PROPERTY_USAGE_PROPERTY_BYTES);
    maximum_property
        .to_bytes()
        .expect("the exact public property-byte cap is admitted");

    let mut oversized_property = receipt.clone();
    oversized_property.property = "p".repeat(MAX_PROPERTY_USAGE_PROPERTY_BYTES + 1);
    assert!(matches!(
        oversized_property.to_bytes(),
        Err(PropertyUsageReceiptError::ResourceLimit {
            resource: "property-bytes",
            ..
        })
    ));

    let mut unordered_axes = receipt.clone();
    unordered_axes.query_point = vec![("z".to_string(), 1.0), ("a".to_string(), 2.0)];
    assert!(matches!(
        unordered_axes.to_bytes(),
        Err(PropertyUsageReceiptError::InvalidField {
            field: "query-point",
            ..
        })
    ));

    let mut duplicate_axes = receipt.clone();
    duplicate_axes.query_point = vec![("T".to_string(), 0.0), ("T".to_string(), 293.15)];
    assert!(matches!(
        duplicate_axes.to_bytes(),
        Err(PropertyUsageReceiptError::InvalidField {
            field: "query-point",
            ..
        })
    ));
    assert!(matches!(
        set.verify_receipt(&duplicate_axes),
        Err(MatDbError::PropertyUsageReceiptNotPortable {
            error: PropertyUsageReceiptError::InvalidField {
                field: "query-point",
                ..
            }
        })
    ));

    let mut foreign_policy = receipt.clone();
    foreign_policy.policy = "trust-me";
    assert!(matches!(
        foreign_policy.to_bytes(),
        Err(PropertyUsageReceiptError::UnknownPolicyTag { .. })
    ));

    let mut stale_schema = receipt.clone();
    stale_schema.schema_version -= 1;
    assert!(matches!(
        stale_schema.to_bytes(),
        Err(PropertyUsageReceiptError::UnsupportedSchemaVersion { .. })
    ));

    let mut detached_source = receipt.clone();
    detached_source.source_hashes[0] = hash_bytes(b"not the selected claim");
    assert!(matches!(
        detached_source.to_bytes(),
        Err(PropertyUsageReceiptError::InvalidField {
            field: "source-hashes",
            ..
        })
    ));
}

// ---------------------------------------------------------------------------
// Pinned-claim selection (bead f85xj.6.4): conflicting claims resolve ONLY
// by an explicit recorded pin, and a pin never bypasses the extrapolation
// refusal.
// ---------------------------------------------------------------------------

#[test]
fn conflicting_claims_resolve_only_by_an_explicit_pin() {
    let mut set = ClaimSet::new();
    let handbook = set
        .insert_claim(density(2699.0, "handbook", stated()))
        .expect("handbook claim inserts");
    let vendor = set
        .insert_claim(density(2712.0, "vendor datasheet", stated()))
        .expect("vendor claim inserts");

    // No policy auto-picks between coexisting conflicting sources.
    assert!(matches!(
        set.query("density", &room(), SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::AmbiguousSelection { ref candidates, .. }) if candidates.len() == 2
    ));
    assert!(matches!(
        set.query("density", &room(), SelectionPolicy::PreferObservationBacked),
        Err(MatDbError::AmbiguousSelection { .. })
    ));

    // The pin is the explicit selection, and the receipt records it as such.
    let answer = set
        .query_pinned("density", &room(), vendor)
        .expect("pinned query answers");
    assert_eq!(answer.evidence.value.value, 2712.0);
    let receipt = &answer.receipt;
    assert_eq!(receipt.policy, PINNED_CLAIM_POLICY_TAG);
    assert_eq!(receipt.selected, vendor);
    assert_eq!(receipt.considered.len(), 2);
    assert_eq!(receipt.in_domain.len(), 2);
    assert!(receipt.considered.contains(&handbook));
    assert_eq!(receipt.source_hashes[0], vendor.0);

    // The pinned receipt replays against the live set.
    set.verify_receipt(receipt).expect("pinned receipt replays");

    // And the other pin selects the other claim: no hidden preference.
    let other = set
        .query_pinned("density", &room(), handbook)
        .expect("other pin answers");
    assert_eq!(other.evidence.value.value, 2699.0);
}

#[test]
fn a_pin_never_bypasses_domain_and_a_foreign_pin_refuses() {
    let mut set = ClaimSet::new();
    let wide = density(2700.0, "wide-range source", stated());
    let mut narrow = density(2705.0, "narrow-range source", stated());
    narrow.validity = ValidityDomain::unconstrained().with("T", 250.0, 300.0);
    set.insert_claim(wide).expect("wide claim inserts");
    let narrow_id = set.insert_claim(narrow).expect("narrow claim inserts");

    // 350 K is inside the wide claim but outside the pinned narrow one:
    // the pin refuses instead of silently extrapolating or substituting.
    let hot = QueryPoint::new().with("T", 350.0).expect("finite point");
    assert!(matches!(
        set.query_pinned("density", &hot, narrow_id),
        Err(MatDbError::PinnedClaimOutOfDomain { ref property, pinned })
            if property == "density" && pinned == narrow_id
    ));

    // A pin that names no claim under the property refuses by name.
    let foreign = fs_matdb::ClaimId(hash_bytes(b"no such claim"));
    assert!(matches!(
        set.query_pinned("density", &room(), foreign),
        Err(MatDbError::PinnedClaimUnknown { ref property, pinned })
            if property == "density" && pinned == foreign
    ));

    // Outside EVERY claim the ordinary extrapolation refusal still fires
    // first: a pin is a selection, not a domain waiver.
    let molten = QueryPoint::new().with("T", 500.0).expect("finite point");
    assert!(matches!(
        set.query_pinned("density", &molten, narrow_id),
        Err(MatDbError::NoClaimInDomain { .. })
    ));
}

#[test]
fn pinned_receipts_survive_the_portable_round_trip() {
    let mut set = ClaimSet::new();
    set.insert_claim(density(2699.0, "handbook", stated()))
        .expect("handbook claim inserts");
    let vendor = set
        .insert_claim(density(2712.0, "vendor datasheet", stated()))
        .expect("vendor claim inserts");

    let receipt = set
        .query_pinned("density", &room(), vendor)
        .expect("pinned query answers")
        .receipt;
    let identity = receipt.try_content_hash().expect("portable identity");
    let bytes = receipt.to_bytes().expect("portable encoding");
    let decoded =
        PropertyUsageReceipt::from_bytes_verified(&bytes, identity).expect("verified decode");
    assert_eq!(decoded, receipt);
    assert_eq!(decoded.policy, PINNED_CLAIM_POLICY_TAG);
    set.verify_receipt(&decoded)
        .expect("decoded receipt replays");
}
