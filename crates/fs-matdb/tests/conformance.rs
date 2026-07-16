//! fs-matdb PR-1 conformance: the fail-closed, append-only insertion
//! boundary. Storing a claim asserts nothing about its truth — these
//! tests lock the schema discipline only.

use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MatDbError, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::Dims;

const DENSITY_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, 0, -1, 0]);

fn provenance() -> Provenance {
    Provenance {
        source: "MMPDS-2023 Table 3.2.1".to_string(),
        license: "internal-use".to_string(),
        artifact: None,
    }
}

fn density_claim(value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("density", DENSITY_DIMS),
        value: PropertyValue::Scalar {
            value,
            dims: DENSITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 273.15, 373.15),
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.01,
            confidence: 0.95,
        },
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: provenance(),
    }
}

fn observation() -> ObservationDataset {
    ObservationDataset {
        specimen: "AA6061-T6 plate, longitudinal".to_string(),
        method: "ASTM B311 (density by displacement)".to_string(),
        artifact: hash_bytes(b"raw density table v1"),
        caveats: "single lab, three specimens, no censoring".to_string(),
        provenance: provenance(),
    }
}

#[test]
fn scalar_and_curve_round_trip_with_dims_checked_at_insertion() {
    let mut set = ClaimSet::new();
    let obs = set
        .register_observation(observation())
        .expect("observation registers");
    let mut claim = density_claim(2700.0);
    claim.observations.push(obs);
    let id = set.insert_claim(claim.clone()).expect("density inserts");
    assert_eq!(set.claim(id), Some(&claim));
    assert_eq!(set.registered_dims("density"), Some(DENSITY_DIMS));

    let curve = PropertyClaim {
        key: PropertyKey::new("electrical-conductivity", CONDUCTIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".to_string(),
            abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
            knots: vec![(273.15, 3.8e7), (323.15, 3.5e7), (373.15, 3.2e7)],
            dims: CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 273.15, 373.15),
        uncertainty: UncertaintyModel::HalfWidth {
            half_width: 5.0e5,
            confidence: 0.9,
        },
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![obs],
        provenance: provenance(),
    };
    let curve_id = set.insert_claim(curve.clone()).expect("curve inserts");
    assert_eq!(set.claim(curve_id), Some(&curve));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"round-trip\",\"verdict\":\"pass\",\
         \"detail\":\"scalar+curve claims round-trip; dims registered at insertion\"}}"
    );
}

#[test]
fn dims_gates_refuse_at_the_door() {
    let mut set = ClaimSet::new();
    // Payload dims disagree with the key's dims.
    let mut wrong_payload = density_claim(2700.0);
    wrong_payload.value = PropertyValue::Scalar {
        value: 2700.0,
        dims: CONDUCTIVITY_DIMS,
    };
    assert!(matches!(
        set.insert_claim(wrong_payload),
        Err(MatDbError::DimsMismatch { .. })
    ));

    // Same name re-registered with different dims refuses even when the
    // payload matches its own (wrong) key.
    set.insert_claim(density_claim(2700.0))
        .expect("first density");
    let alias = PropertyClaim {
        key: PropertyKey::new("density", CONDUCTIVITY_DIMS),
        value: PropertyValue::Scalar {
            value: 1.0,
            dims: CONDUCTIVITY_DIMS,
        },
        ..density_claim(1.0)
    };
    assert!(matches!(
        set.insert_claim(alias),
        Err(MatDbError::DimsMismatch { .. })
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"dims-gate\",\"verdict\":\"pass\",\
         \"detail\":\"payload/key and key/registry dims mismatches refuse\"}}"
    );
}

#[test]
fn provenance_is_load_bearing() {
    let mut set = ClaimSet::new();
    let mut unlicensed = density_claim(2700.0);
    unlicensed.provenance.license = "  ".to_string();
    assert!(matches!(
        set.insert_claim(unlicensed),
        Err(MatDbError::MissingLicense { .. })
    ));

    let mut sourceless = density_claim(2700.0);
    sourceless.provenance.source = String::new();
    assert!(matches!(
        set.insert_claim(sourceless),
        Err(MatDbError::MissingSource)
    ));

    let mut bad_observation = observation();
    bad_observation.provenance.license = String::new();
    assert!(matches!(
        set.register_observation(bad_observation),
        Err(MatDbError::MissingLicense { .. })
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"provenance-gate\",\"verdict\":\"pass\",\
         \"detail\":\"missing license/source refuse for claims and observations\"}}"
    );
}

#[test]
fn payload_and_validity_pathologies_refuse() {
    let mut set = ClaimSet::new();
    assert!(matches!(
        set.insert_claim(density_claim(f64::NAN)),
        Err(MatDbError::NonFinite {
            field: "scalar value",
            ..
        })
    ));

    let mut nan_axis = density_claim(2700.0);
    nan_axis.validity = ValidityDomain::unconstrained().with("T", f64::NAN, 400.0);
    assert!(matches!(
        set.insert_claim(nan_axis),
        Err(MatDbError::UnusableValidity { .. })
    ));

    let mut bad_confidence = density_claim(2700.0);
    bad_confidence.uncertainty = UncertaintyModel::HalfWidth {
        half_width: 1.0,
        confidence: 1.0,
    };
    assert!(matches!(
        set.insert_claim(bad_confidence),
        Err(MatDbError::InvalidUncertainty { .. })
    ));

    let mut negative_width = density_claim(2700.0);
    negative_width.uncertainty = UncertaintyModel::HalfWidth {
        half_width: -1.0,
        confidence: 0.9,
    };
    assert!(matches!(
        set.insert_claim(negative_width),
        Err(MatDbError::InvalidUncertainty { .. })
    ));

    let mut short_curve = density_claim(2700.0);
    short_curve.value = PropertyValue::Curve {
        abscissa: "T".to_string(),
        abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
        knots: vec![(300.0, 2700.0)],
        dims: DENSITY_DIMS,
    };
    assert!(matches!(
        set.insert_claim(short_curve),
        Err(MatDbError::MalformedCurve { .. })
    ));

    let mut unordered = density_claim(2700.0);
    unordered.value = PropertyValue::Curve {
        abscissa: "T".to_string(),
        abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
        knots: vec![(300.0, 2700.0), (300.0, 2690.0)],
        dims: DENSITY_DIMS,
    };
    assert!(matches!(
        set.insert_claim(unordered),
        Err(MatDbError::MalformedCurve { .. })
    ));

    let mut dangling = density_claim(2700.0);
    dangling.observations = vec![fs_matdb::ObservationId(hash_bytes(b"never registered"))];
    assert!(matches!(
        set.insert_claim(dangling),
        Err(MatDbError::UnknownObservation { .. })
    ));
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"payload-gates\",\"verdict\":\"pass\",\
         \"detail\":\"NaN/validity/uncertainty/curve/dangling-ref refusals all typed\"}}"
    );
}

#[test]
fn conflicting_claims_coexist_and_nothing_overwrites() {
    let mut set = ClaimSet::new();
    let a = set.insert_claim(density_claim(2700.0)).expect("claim a");
    let b = set.insert_claim(density_claim(2698.5)).expect("claim b");
    assert_ne!(a, b, "different content gets a different id");
    let all = set.claims_for("density");
    assert_eq!(all.len(), 2, "conflicting claims BOTH survive");
    assert_eq!(all[0].0, a, "insertion order preserved");
    assert_eq!(all[1].0, b);
    assert_eq!(set.claim_count(), 2);

    // Idempotent re-insertion: same content, same id, no duplicate.
    let a_again = set.insert_claim(density_claim(2700.0)).expect("re-insert");
    assert_eq!(a, a_again);
    assert_eq!(set.claims_for("density").len(), 2);
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"no-overwrite\",\"verdict\":\"pass\",\
         \"detail\":\"conflicting claims coexist; re-insertion idempotent by content id\"}}"
    );
}

#[test]
fn content_identity_is_stable_and_field_sensitive() {
    let base = density_claim(2700.0);
    let same = density_claim(2700.0);
    assert_eq!(base.content_hash(), same.content_hash());

    let mut moved_value = density_claim(2700.0);
    moved_value.value = PropertyValue::Scalar {
        value: 2700.0000001,
        dims: DENSITY_DIMS,
    };
    assert_ne!(base.content_hash(), moved_value.content_hash());

    let mut moved_validity = density_claim(2700.0);
    moved_validity.validity = ValidityDomain::unconstrained().with("T", 273.15, 374.15);
    assert_ne!(base.content_hash(), moved_validity.content_hash());

    let mut moved_uncertainty = density_claim(2700.0);
    moved_uncertainty.uncertainty = UncertaintyModel::Unstated;
    assert_ne!(base.content_hash(), moved_uncertainty.content_hash());

    let mut moved_policy = density_claim(2700.0);
    moved_policy.interpolation = InterpolationPolicy::TabulatedOnly;
    assert_ne!(base.content_hash(), moved_policy.content_hash());

    let mut moved_license = density_claim(2700.0);
    moved_license.provenance.license = "CC-BY-4.0".to_string();
    assert_ne!(base.content_hash(), moved_license.content_hash());

    let mut set = ClaimSet::new();
    let obs_a = set
        .register_observation(observation())
        .expect("observation registers");
    let obs_b = set.register_observation(observation()).expect("idempotent");
    assert_eq!(obs_a, obs_b, "observation registration idempotent");
    println!(
        "{{\"suite\":\"fs-matdb\",\"case\":\"content-identity\",\"verdict\":\"pass\",\
         \"detail\":\"hash stable on equal content, moves on every semantic field\"}}"
    );
}
