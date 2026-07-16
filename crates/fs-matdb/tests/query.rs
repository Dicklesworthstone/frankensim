//! fs-matdb PR-4 conformance: every answer is Evidence + receipt;
//! extrapolation refuses; fusion is explicit; unstated uncertainty
//! never launders into a certificate.

use fs_blake3::hash_bytes;
use fs_evidence::NumericalKind;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, EvaluationDecision, InterpolationPolicy, MATDB_EVALUATOR_VERSION, MatDbError,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyUsageReceipt, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
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
        Err(MatDbError::ReceiptMismatch { field: "selected" })
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
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"receipt-mutation-battery\",\"verdict\":\"pass\",\
         \"detail\":\"11 field mutations all refuse typed and all move the receipt identity\"}}"
    );
}
