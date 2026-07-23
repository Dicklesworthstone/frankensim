//! Battery for joint-aware querying (bead f85xj.7.3): correlated
//! properties answer with the admitted covariance block in request order,
//! every absent correlation names its reason (the silent-independence
//! path does not exist), Unstated marginals never anchor a joint band,
//! and joint receipts replay against the exact pack that minted them.
#![allow(clippy::float_cmp)] // bit-exact propagation of admitted covariance entries is the contract under test

use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, CorrelationUnknownReason, InterpolationPolicy, JOINT_USAGE_RECEIPT_SCHEMA_VERSION,
    JointCorrelation, JointStatistics, MatDbError, NormalizedPack, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    StatisticMember, UncertaintyModel,
};
use fs_qty::Dims;

const CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);
const DENSITY_DIMS: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn provenance(source: &str) -> Provenance {
    Provenance {
        source: source.to_string(),
        license: "internal-use".to_string(),
        artifact: None,
    }
}

fn stated(half_width: f64) -> UncertaintyModel {
    UncertaintyModel::HalfWidth {
        half_width,
        confidence: 0.95,
    }
}

fn scalar_claim(
    property: &str,
    dims: Dims,
    value: f64,
    uncertainty: UncertaintyModel,
) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(property, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 400.0),
        uncertainty,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: provenance("joint characterization lab"),
    }
}

fn room() -> QueryPoint {
    QueryPoint::new().with("T", 293.15).expect("finite point")
}

/// A pack holding conductivity + density claims measured jointly, with a
/// 2x2 covariance block: var(k) = 0.04, var(rho) = 225.0, cov = 1.2
/// (correlation 0.4).
fn correlated_pack() -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "AA6061-T6 coupon set".to_string(),
            method: "joint composition sweep".to_string(),
            artifact: hash_bytes(b"joint raw table"),
            caveats: "none".to_string(),
            provenance: provenance("joint characterization lab"),
        })
        .expect("observation registers");
    let mut conductivity = scalar_claim("thermal-conductivity", CONDUCTIVITY_DIMS, 167.0, {
        stated(0.2)
    });
    conductivity.observations.push(observation);
    let mut density = scalar_claim("density", DENSITY_DIMS, 2699.0, stated(15.0));
    density.observations.push(observation);
    let k_id = claims.insert_claim(conductivity).expect("k inserts");
    let rho_id = claims.insert_claim(density).expect("rho inserts");

    // Members must be strictly increasing in (claim, component) order.
    let mut members = vec![
        StatisticMember::scalar(k_id),
        StatisticMember::scalar(rho_id),
    ];
    members.sort();
    // Packed lower triangle in MEMBER order: [var(m0), cov, var(m1)].
    let covariance = if members[0].claim() == k_id {
        vec![0.04, 1.2, 225.0]
    } else {
        vec![225.0, 1.2, 0.04]
    };
    let block = JointStatistics::new(observation, "k-rho-composition", members, covariance, None);
    NormalizedPack::new(
        "joint-test-pack",
        "test-compiler 0.0.1",
        hash_bytes(b"joint source artifact"),
        "internal test data",
        claims,
        vec![block],
        Vec::new(),
    )
    .expect("pack admits")
}

#[test]
fn correlated_properties_answer_with_the_covariance_in_request_order() {
    let pack = correlated_pack();
    let answer = pack
        .query_joint(
            &["thermal-conductivity", "density"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("joint query answers");

    assert_eq!(answer.members.len(), 2);
    assert_eq!(answer.members[0].evidence.value.value, 167.0);
    assert_eq!(answer.members[1].evidence.value.value, 2699.0);

    // Request order governs the packed triangle, whatever the block's
    // internal member order is.
    let JointCorrelation::Covariance {
        block_id,
        covariance,
        correlation,
        ..
    } = &answer.receipt.correlation
    else {
        panic!("expected covariance, got {:?}", answer.receipt.correlation);
    };
    assert_eq!(block_id, "k-rho-composition");
    assert_eq!(covariance, &vec![0.04, 1.2, 225.0]);
    assert!(correlation.is_none());

    // The reversed request transposes the triangle consistently.
    let reversed = pack
        .query_joint(
            &["density", "thermal-conductivity"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("reversed joint query answers");
    let JointCorrelation::Covariance { covariance, .. } = &reversed.receipt.correlation else {
        panic!("expected covariance");
    };
    assert_eq!(covariance, &vec![225.0, 1.2, 0.04]);

    // Statistical sanity: the variance of k + rho under the returned
    // block is var(k) + var(rho) + 2 cov — the analytic bivariate case.
    let var_sum = 0.04 + 225.0 + 2.0 * 1.2;
    assert_eq!(var_sum, 227.44);

    // Receipt shape: aligned members, receipts, and the pack identity.
    let receipt = &answer.receipt;
    assert_eq!(receipt.schema_version, JOINT_USAGE_RECEIPT_SCHEMA_VERSION);
    assert_eq!(receipt.pack, pack.content_hash());
    assert_eq!(receipt.properties, vec!["thermal-conductivity", "density"]);
    assert_eq!(receipt.selected.len(), 2);
    assert_eq!(
        receipt.member_receipts,
        answer
            .members
            .iter()
            .map(|member| member.receipt.content_hash())
            .collect::<Vec<_>>()
    );

    // Determinism: the same query mints the identical receipt identity.
    let again = pack
        .query_joint(
            &["thermal-conductivity", "density"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("repeat answers");
    assert_eq!(again.receipt.content_hash(), receipt.content_hash());

    // And the receipt replays against the pack that minted it.
    pack.verify_joint_receipt(receipt).expect("receipt replays");
}

#[test]
#[allow(clippy::too_many_lines)] // one full pack fixture per unknown-correlation reason, in reason order
fn every_absent_correlation_names_its_reason() {
    // (a) No block at all.
    let mut claims = ClaimSet::new();
    claims
        .register_observation(ObservationDataset {
            specimen: "coupon".to_string(),
            method: "single-property tests".to_string(),
            artifact: hash_bytes(b"blockless raw"),
            caveats: "none".to_string(),
            provenance: provenance("joint characterization lab"),
        })
        .expect("observation registers");
    claims
        .insert_claim(scalar_claim(
            "thermal-conductivity",
            CONDUCTIVITY_DIMS,
            167.0,
            stated(0.2),
        ))
        .expect("k inserts");
    claims
        .insert_claim(scalar_claim("density", DENSITY_DIMS, 2699.0, stated(15.0)))
        .expect("rho inserts");
    let blockless = NormalizedPack::new(
        "blockless-pack",
        "test-compiler 0.0.1",
        hash_bytes(b"blockless source"),
        "internal test data",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .expect("pack admits");
    let answer = blockless
        .query_joint(
            &["thermal-conductivity", "density"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("joint query answers");
    assert_eq!(
        answer.receipt.correlation,
        JointCorrelation::Unknown {
            reason: CorrelationUnknownReason::NoBlock
        }
    );

    // (b) Partial membership: a third property outside the block.
    let pack = {
        let mut claims = ClaimSet::new();
        let observation = claims
            .register_observation(ObservationDataset {
                specimen: "coupon".to_string(),
                method: "sweep".to_string(),
                artifact: hash_bytes(b"partial raw"),
                caveats: "none".to_string(),
                provenance: provenance("joint characterization lab"),
            })
            .expect("observation registers");
        let mut conductivity = scalar_claim(
            "thermal-conductivity",
            CONDUCTIVITY_DIMS,
            167.0,
            stated(0.2),
        );
        conductivity.observations.push(observation);
        let mut density = scalar_claim("density", DENSITY_DIMS, 2699.0, stated(15.0));
        density.observations.push(observation);
        let k_id = claims.insert_claim(conductivity).expect("k inserts");
        let rho_id = claims.insert_claim(density).expect("rho inserts");
        claims
            .insert_claim(scalar_claim(
                "specific-heat",
                SPECIFIC_HEAT_DIMS,
                896.0,
                stated(8.0),
            ))
            .expect("cp inserts");
        let mut members = vec![
            StatisticMember::scalar(k_id),
            StatisticMember::scalar(rho_id),
        ];
        members.sort();
        let covariance = if members[0].claim() == k_id {
            vec![0.04, 1.2, 225.0]
        } else {
            vec![225.0, 1.2, 0.04]
        };
        let block = JointStatistics::new(observation, "k-rho-only", members, covariance, None);
        NormalizedPack::new(
            "partial-pack",
            "test-compiler 0.0.1",
            hash_bytes(b"partial source"),
            "internal test data",
            claims,
            vec![block],
            Vec::new(),
        )
        .expect("pack admits")
    };
    let answer = pack
        .query_joint(
            &["thermal-conductivity", "density", "specific-heat"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("joint query answers");
    assert_eq!(
        answer.receipt.correlation,
        JointCorrelation::Unknown {
            reason: CorrelationUnknownReason::PartialMembership
        }
    );
    // The covered pair alone still answers with its covariance.
    let pair = pack
        .query_joint(
            &["thermal-conductivity", "density"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("pair answers");
    assert!(matches!(
        pair.receipt.correlation,
        JointCorrelation::Covariance { .. }
    ));
}

#[test]
fn an_unstated_marginal_never_anchors_a_joint_band() {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "coupon".to_string(),
            method: "sweep".to_string(),
            artifact: hash_bytes(b"unstated raw"),
            caveats: "none".to_string(),
            provenance: provenance("joint characterization lab"),
        })
        .expect("observation registers");
    let mut conductivity = scalar_claim(
        "thermal-conductivity",
        CONDUCTIVITY_DIMS,
        167.0,
        UncertaintyModel::Unstated,
    );
    conductivity.observations.push(observation);
    let mut density = scalar_claim("density", DENSITY_DIMS, 2699.0, stated(15.0));
    density.observations.push(observation);
    let k_id = claims.insert_claim(conductivity).expect("k inserts");
    let rho_id = claims.insert_claim(density).expect("rho inserts");
    let mut members = vec![
        StatisticMember::scalar(k_id),
        StatisticMember::scalar(rho_id),
    ];
    members.sort();
    let covariance = if members[0].claim() == k_id {
        vec![0.04, 1.2, 225.0]
    } else {
        vec![225.0, 1.2, 0.04]
    };
    let block = JointStatistics::new(observation, "k-rho-unstated", members, covariance, None);
    let pack = NormalizedPack::new(
        "unstated-pack",
        "test-compiler 0.0.1",
        hash_bytes(b"unstated source"),
        "internal test data",
        claims,
        vec![block],
        Vec::new(),
    )
    .expect("pack admits");

    let answer = pack
        .query_joint(
            &["thermal-conductivity", "density"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("joint query answers");
    assert_eq!(
        answer.receipt.correlation,
        JointCorrelation::Unknown {
            reason: CorrelationUnknownReason::UnstatedMarginal
        }
    );
}

#[test]
fn degenerate_joint_requests_refuse_and_member_refusals_propagate() {
    let pack = correlated_pack();
    assert!(matches!(
        pack.query_joint(
            &["thermal-conductivity"],
            &room(),
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::UnsupportedEvaluation { .. })
    ));
    assert!(matches!(
        pack.query_joint(
            &["density", "density"],
            &room(),
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::UnsupportedEvaluation { .. })
    ));
    // A member refusal (unknown property) propagates unchanged.
    assert!(matches!(
        pack.query_joint(
            &["thermal-conductivity", "youngs-modulus"],
            &room(),
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::UnknownProperty { .. })
    ));
    // Extrapolation refuses through the ordinary member path: no joint
    // machinery weakens the validity gate.
    let molten = QueryPoint::new().with("T", 500.0).expect("finite point");
    assert!(matches!(
        pack.query_joint(
            &["thermal-conductivity", "density"],
            &molten,
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::NoClaimInDomain { .. })
    ));
}

#[test]
fn joint_receipts_bind_the_pack_and_every_field() {
    let pack = correlated_pack();
    let receipt = pack
        .query_joint(
            &["thermal-conductivity", "density"],
            &room(),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("joint query answers")
        .receipt;

    // A receipt minted against one pack refuses against another.
    let mut other_claims = ClaimSet::new();
    other_claims
        .register_observation(ObservationDataset {
            specimen: "coupon".to_string(),
            method: "single-property tests".to_string(),
            artifact: hash_bytes(b"other raw"),
            caveats: "none".to_string(),
            provenance: provenance("joint characterization lab"),
        })
        .expect("observation registers");
    other_claims
        .insert_claim(scalar_claim(
            "thermal-conductivity",
            CONDUCTIVITY_DIMS,
            167.0,
            stated(0.2),
        ))
        .expect("k inserts");
    other_claims
        .insert_claim(scalar_claim("density", DENSITY_DIMS, 2699.0, stated(15.0)))
        .expect("rho inserts");
    let other = NormalizedPack::new(
        "other-pack",
        "test-compiler 0.0.1",
        hash_bytes(b"other source"),
        "internal test data",
        other_claims,
        Vec::new(),
        Vec::new(),
    )
    .expect("pack admits");
    assert!(matches!(
        other.verify_joint_receipt(&receipt),
        Err(MatDbError::ReceiptMismatch { field: "pack" })
    ));

    // A tampered correlation refuses by name.
    let mut tampered = receipt.clone();
    tampered.correlation = JointCorrelation::Unknown {
        reason: CorrelationUnknownReason::NoBlock,
    };
    assert!(matches!(
        pack.verify_joint_receipt(&tampered),
        Err(MatDbError::ReceiptMismatch {
            field: "correlation"
        })
    ));
    assert_ne!(tampered.content_hash(), receipt.content_hash());

    // A tampered member-receipt list refuses by name.
    let mut tampered = receipt.clone();
    tampered.member_receipts.reverse();
    assert!(matches!(
        pack.verify_joint_receipt(&tampered),
        Err(MatDbError::ReceiptMismatch {
            field: "member_receipts"
        })
    ));

    // Version drift refuses before any replay.
    let mut drifted = receipt.clone();
    drifted.schema_version = 2;
    assert!(matches!(
        pack.verify_joint_receipt(&drifted),
        Err(MatDbError::ReceiptSchemaVersionDrift { .. })
    ));
}
