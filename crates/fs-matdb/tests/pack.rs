//! G0/G3 conformance for the normalized material-pack boundary.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSet, InterpolationPolicy, JointStatistics, MATDB_PACK_SCHEMA_VERSION,
    MATDB_PACK_TARGET_BASIS, NormalizationReceipt, NormalizationTarget, NormalizedPack,
    ObservationDataset, ObservationId, PackError, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, StatisticMember, UncertaintyModel,
};
use fs_qty::Dims;

const SOURCE_DOMAIN: &str = "org.frankensim.tests.matdb-pack.source.v1";

fn push_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_string(bytes: &mut Vec<u8>, value: &str) {
    push_u32(
        bytes,
        u32::try_from(value.len()).expect("short test string"),
    );
    bytes.extend_from_slice(value.as_bytes());
}

fn bare_pack_prefix() -> Vec<u8> {
    let mut bytes = b"FSMATPK\0".to_vec();
    push_u32(&mut bytes, MATDB_PACK_SCHEMA_VERSION);
    push_string(&mut bytes, "preflight-fixture");
    push_string(&mut bytes, "compiler-v1");
    bytes.extend_from_slice(&[0; 32]);
    push_string(&mut bytes, "redistribution permitted");
    bytes
}

fn provenance() -> Provenance {
    Provenance {
        source: "fixture handbook table 7".to_string(),
        license: "CC-BY-4.0; redistribution permitted with attribution".to_string(),
        artifact: Some(hash_domain(SOURCE_DOMAIN, b"fixture-table")),
    }
}

fn sample_claims() -> (ClaimSet, ObservationId, Vec<StatisticMember>) {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "alloy-X; solution-treated".to_string(),
            method: "ASTM fixture method".to_string(),
            artifact: hash_domain(SOURCE_DOMAIN, b"raw-observation"),
            caveats: "joint density/modulus coupon series".to_string(),
            provenance: provenance(),
        })
        .expect("licensed fixture observation");

    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let density = claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", density_dims),
            value: PropertyValue::Scalar {
                value: 7_850.0,
                dims: density_dims,
            },
            validity: ValidityDomain::unconstrained().with("temperature", 273.15, 373.15),
            uncertainty: UncertaintyModel::HalfWidth {
                half_width: 5.0,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .expect("density fixture claim");

    let pressure_dims = Dims([-1, 1, -2, 0, 0, 0]);
    let modulus = claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("young_modulus", pressure_dims),
            value: PropertyValue::Curve {
                abscissa: "temperature".to_string(),
                abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                knots: vec![(273.15, 210.0e9), (373.15, 202.0e9)],
                dims: pressure_dims,
            },
            validity: ValidityDomain::unconstrained().with("temperature", 273.15, 373.15),
            uncertainty: UncertaintyModel::RelativeHalfWidth {
                fraction: 0.02,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: provenance(),
        })
        .expect("modulus fixture claim");

    let mut members = vec![
        StatisticMember::scalar(density),
        StatisticMember::curve_ordinate(modulus, 0),
    ];
    members.sort_unstable();
    (claims, observation, members)
}

#[test]
fn g3_strain_context_pack_queries_and_receipts_preserve_all_coordinates() {
    use fs_blake3::ContentHash;
    use fs_matdb::{
        ElasticTensorOrder, MatDbError, QueryPoint, SelectionPolicy, StrainTensorBasis,
        StrainTensorComponent, StrainTensorNotation,
    };
    let basis = StrainTensorBasis {
        notation: StrainTensorNotation::Engineering,
        order: ElasticTensorOrder::XxYyZzXyYzZx,
        frame: ContentHash([1; 32]),
    };
    let component = StrainTensorComponent::new(basis, ContentHash([2; 32]), 3).unwrap();
    let contexts = [
        component,
        StrainTensorComponent::new(
            StrainTensorBasis {
                notation: StrainTensorNotation::Tensor,
                ..basis
            },
            ContentHash([2; 32]),
            3,
        )
        .unwrap(),
        StrainTensorComponent::new(
            StrainTensorBasis {
                order: ElasticTensorOrder::XxYyZzYzZxXy,
                ..basis
            },
            ContentHash([2; 32]),
            3,
        )
        .unwrap(),
        StrainTensorComponent::new(
            StrainTensorBasis {
                frame: ContentHash([3; 32]),
                ..basis
            },
            ContentHash([2; 32]),
            3,
        )
        .unwrap(),
        StrainTensorComponent::new(basis, ContentHash([3; 32]), 3).unwrap(),
        StrainTensorComponent::new(basis, ContentHash([2; 32]), 4).unwrap(),
    ];
    assert!(StrainTensorComponent::new(basis, ContentHash([0; 32]), 0).is_err());
    assert!(StrainTensorComponent::new(basis, ContentHash([2; 32]), 6).is_err());
    assert!(
        StrainTensorComponent::new(
            StrainTensorBasis {
                frame: ContentHash([0; 32]),
                ..basis
            },
            ContentHash([2; 32]),
            0
        )
        .is_err()
    );
    for dims in [fs_qty::Pressure::DIMS, Dims([0, 0, 0, -1, 0, 0])] {
        assert!(
            PropertyKey::new("bad", dims)
                .with_strain_component(component)
                .is_err()
        );
    }
    let (mut claims, observation, _) = sample_claims();
    let bare = PropertyKey::new("strain", Dims::NONE);
    let mut keys = Vec::new();
    let mut ids = std::collections::BTreeSet::new();
    for context in contexts {
        assert_eq!(
            StrainTensorComponent::from_canonical_bytes(&context.canonical_bytes()).unwrap(),
            context
        );
        let key = bare.clone().with_strain_component(context).unwrap();
        let id = claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value: -0.004,
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 300.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::TabulatedOnly,
                observations: vec![observation],
                provenance: provenance(),
            })
            .unwrap();
        assert!(ids.insert(id), "each source context field changes identity");
        keys.push((key, id));
    }
    let pack = NormalizedPack::new(
        "synthetic-strain",
        "fixture-v6",
        ContentHash([4; 32]),
        "synthetic redistribution permitted",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    assert_eq!(pack.schema_version(), 6);
    let bytes = pack.to_bytes();
    let decoded = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(decoded.to_bytes(), bytes);
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    assert!(matches!(
        decoded
            .claims()
            .query("strain", &point, SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::TensorContextMismatch { .. })
    ));
    assert!(matches!(
        decoded
            .claims()
            .query_typed(&bare, &point, SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::TensorContextMismatch { .. })
    ));
    for (key, id) in &keys {
        for answer in [
            decoded
                .claims()
                .query_typed(key, &point, SelectionPolicy::SingleClaimOnly)
                .unwrap(),
            decoded
                .claims()
                .query_pinned_typed(key, &point, *id)
                .unwrap(),
        ] {
            assert_eq!(answer.evidence.value.value, -0.004);
            assert_eq!(answer.receipt.selected, *id);
            assert_eq!(
                decoded.claims().claim(*id).unwrap().key.strain_component(),
                key.strain_component()
            );
            decoded.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    assert!(
        decoded
            .claims()
            .query_pinned_typed(&keys[0].0, &point, keys[1].1)
            .is_err()
    );
    assert!(
        decoded
            .claims()
            .query_typed(
                &keys[0].0,
                &QueryPoint::new().with("T", 301.0).unwrap(),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    let encoded = component.canonical_bytes();
    let offset = bytes
        .windows(encoded.len())
        .position(|part| part == encoded)
        .unwrap();
    for (index, value) in [(0, 255), (1, 255), (2, 6), (3, 9)] {
        let mut tampered = bytes.clone();
        tampered[offset + index] = value;
        assert!(NormalizedPack::from_bytes(&tampered).is_err());
    }
    for version in 1u32..=5 {
        let mut downgraded = bytes.clone();
        downgraded[8..12].copy_from_slice(&version.to_le_bytes());
        assert!(NormalizedPack::from_bytes(&downgraded).is_err());
    }
    assert!(StrainTensorComponent::from_canonical_bytes(&encoded[..66]).is_err());
}

fn sample_pack() -> NormalizedPack {
    let (claims, observation, members) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let modulus = claims.claims_for("young_modulus")[0].0;
    NormalizedPack::new(
        "fixture-alloy-x",
        "frankensim-matdb-pack-compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source-envelope"),
        "CC-BY-4.0: redistribution permitted with attribution",
        claims,
        vec![JointStatistics::new(
            observation,
            "coupon-joint-density-modulus",
            members,
            vec![4.0, 1.0, 9.0],
            Some(vec![1.0, 1.0 / 6.0, 1.0]),
        )],
        vec![
            NormalizationReceipt::new(
                NormalizationTarget::ClaimUncertainty { claim: density },
                hash_domain(SOURCE_DOMAIN, b"plus-or-minus 0.005 g/cm3"),
                Dims([-3, 1, 0, 0, 0, 0]),
                1_000.0,
                0.0,
                "g/cm3",
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ),
            NormalizationReceipt::new(
                NormalizationTarget::ClaimValue(StatisticMember::curve_ordinate(modulus, 0)),
                hash_domain(SOURCE_DOMAIN, b"210 GPa"),
                Dims([-1, 1, -2, 0, 0, 0]),
                1.0e9,
                0.0,
                "GPa",
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ),
            NormalizationReceipt::new(
                NormalizationTarget::ClaimValue(StatisticMember::scalar(density)),
                hash_domain(SOURCE_DOMAIN, b"7.85 g/cm3"),
                Dims([-3, 1, 0, 0, 0, 0]),
                1_000.0,
                0.0,
                "g/cm3",
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ),
            NormalizationReceipt::new(
                NormalizationTarget::ValidityBound {
                    claim: density,
                    axis: "temperature".to_string(),
                    side: fs_matdb::ValidityBoundSide::Lower,
                },
                hash_domain(SOURCE_DOMAIN, b"0 degC"),
                Dims([0, 0, 0, 1, 0, 0]),
                1.0,
                273.15,
                "degC",
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ),
            NormalizationReceipt::new(
                NormalizationTarget::JointCovariance {
                    observation,
                    block_id: "coupon-joint-density-modulus".to_string(),
                    row: 1,
                    column: 0,
                },
                hash_domain(SOURCE_DOMAIN, b"joint covariance table entry"),
                Dims([-4, 2, -2, 0, 0, 0]),
                1.0,
                0.0,
                "coherent source covariance",
                MATDB_PACK_TARGET_BASIS,
                None,
                None,
            ),
        ],
    )
    .expect("valid fixture pack")
}

#[test]
fn normalized_pack_round_trips_exact_bytes_and_semantics() {
    let pack = sample_pack();
    let bytes = pack.to_bytes();
    // Independently reconstructed from the v1 field grammar and BLAKE3
    // derive-key domains; this is deliberately not derived from the decoder.
    assert_eq!(bytes.len(), 1_951);
    assert_eq!(
        pack.content_hash().to_hex(),
        "96ec60562b2219447093f90bfd8bbbc9c87f059cfb93976c825f5caa5f6a82d2"
    );
    let decoded = NormalizedPack::from_bytes(&bytes).expect("canonical pack decodes");
    let density = decoded.claims().claims_for("density")[0].0;
    let modulus = decoded.claims().claims_for("young_modulus")[0].0;

    assert_eq!(decoded, pack);
    assert_eq!(decoded.to_bytes(), bytes);
    assert_eq!(decoded.content_hash(), pack.content_hash());
    assert_eq!(
        NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes)
            .expect("externally pinned canonical pack decodes"),
        pack
    );
    assert_eq!(decoded.claims().claim_count(), 2);
    assert_eq!(decoded.joint_statistics()[0].covariance(), &[4.0, 1.0, 9.0]);
    assert_eq!(
        decoded.joint_statistics()[0].correlation(),
        Some(&[1.0, 1.0 / 6.0, 1.0][..])
    );
    assert_eq!(decoded.normalizations().len(), 5);
    assert!(decoded.normalizations().iter().any(|receipt| {
        receipt.target() == &NormalizationTarget::ClaimValue(StatisticMember::scalar(density))
    }));
    assert!(decoded.normalizations().iter().any(|receipt| {
        receipt.target()
            == &NormalizationTarget::ClaimValue(StatisticMember::curve_ordinate(modulus, 0))
    }));
    assert!(decoded.normalizations().iter().any(|receipt| {
        matches!(
            receipt.target(),
            NormalizationTarget::ClaimUncertainty { claim } if *claim == density
        )
    }));
    assert!(decoded.normalizations().iter().any(|receipt| {
        matches!(
            receipt.target(),
            NormalizationTarget::ValidityBound { claim, axis, .. }
                if *claim == density && axis == "temperature"
        )
    }));
    assert!(decoded.normalizations().iter().any(|receipt| {
        matches!(
            receipt.target(),
            NormalizationTarget::JointCovariance {
                observation,
                block_id,
                row: 1,
                column: 0,
            } if *observation == decoded.joint_statistics()[0].observation()
                && block_id == decoded.joint_statistics()[0].block_id()
        )
    }));
}

#[test]
fn g3_typed_pack_preserves_quantity_identity_and_refuses_version_reinterpretation() {
    use fs_matdb::{MATDB_TYPED_PACK_SCHEMA_VERSION, MatDbError, QueryPoint, SelectionPolicy};
    use fs_qty::{
        QuantitySpec,
        semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
    };
    let kind = QuantityKind::Frequency(FrequencyConvention::Cyclic);
    let quantity = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let key = PropertyKey::with_quantity("reference_frequency", quantity);
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic frequency fixture".to_string(),
            method: "authored codec roundtrip case".to_string(),
            artifact: hash_domain(SOURCE_DOMAIN, b"synthetic observation: 50 Hz"),
            caveats: "synthetic fixture; no physical measurement or validation".to_string(),
            provenance: provenance(),
        })
        .unwrap();
    let id = claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Scalar {
                value: 50.0,
                dims: kind.expected_dims(),
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::TabulatedOnly,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let pack = NormalizedPack::new(
        "frequency-fixture",
        "typed-compiler-v2",
        hash_domain(SOURCE_DOMAIN, b"50 Hz"),
        "synthetic fixture permitted",
        claims,
        Vec::new(),
        vec![NormalizationReceipt::new(
            NormalizationTarget::ClaimValue(StatisticMember::scalar(id)),
            hash_domain(SOURCE_DOMAIN, b"50 Hz"),
            kind.expected_dims(),
            1.0,
            0.0,
            "Hz",
            MATDB_PACK_TARGET_BASIS,
            None,
            None,
        )],
    )
    .unwrap();
    assert_eq!(pack.schema_version(), MATDB_TYPED_PACK_SCHEMA_VERSION);
    let bytes = pack.to_bytes();
    assert_eq!(&bytes[8..12], &2_u32.to_le_bytes());
    let decoded = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(decoded.to_bytes(), bytes);
    assert_eq!(decoded.claims().claim(id).unwrap().key.quantity(), quantity);
    assert_eq!(decoded.normalizations()[0].source_basis(), "Hz");
    let point = QueryPoint::new();
    let answer = decoded
        .claims()
        .query_typed(&key, &point, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert_eq!(answer.evidence.value.value.to_bits(), 50.0_f64.to_bits());
    decoded.claims().verify_receipt(&answer.receipt).unwrap();
    let state = fs_matdb::MaterialStateId {
        chemistry: "synthetic quantity fixture".to_string(),
        phase: "test-only".to_string(),
        process: "declared fixture".to_string(),
        revision: 0,
    };
    let material =
        fs_matdb::NormalizedMaterialCardPack::new(state.clone(), decoded.clone()).unwrap();
    let material = fs_matdb::NormalizedMaterialCardPack::from_bytes_verified(
        material.content_hash(),
        &material.to_bytes(),
    )
    .unwrap();
    let material_answer = material
        .card()
        .claims()
        .query_typed(&key, &point, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert_eq!(material_answer.evidence.value.quantity, quantity);
    material
        .card()
        .claims()
        .verify_receipt(&material_answer.receipt)
        .unwrap();
    let surface = fs_matdb::SurfaceSpec {
        material: state,
        texture_frame: "declared fixture frame".to_string(),
    };
    let interface = fs_matdb::NormalizedInterfacePack::new(
        surface.clone(),
        surface,
        fs_matdb::SystemContext {
            medium: "synthetic interface medium".to_string(),
            third_body: None,
            environment: "synthetic environment".to_string(),
            history: "synthetic frequency response".to_string(),
        },
        decoded.clone(),
    )
    .unwrap();
    let interface = fs_matdb::NormalizedInterfacePack::from_bytes_verified(
        interface.content_hash(),
        &interface.to_bytes(),
    )
    .unwrap();
    let interface_answer = interface
        .card()
        .claims()
        .query_typed(&key, &point, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert_eq!(interface_answer.evidence.value.quantity, quantity);
    interface
        .card()
        .claims()
        .verify_receipt(&interface_answer.receipt)
        .unwrap();
    let wrong = PropertyKey::with_quantity(
        key.name(),
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Frequency(FrequencyConvention::Angular),
            ValueForm::Static,
        )),
    );
    assert!(matches!(
        decoded
            .claims()
            .query_typed(&wrong, &point, SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::QuantityMismatch { .. })
    ));
    let mut downgraded = bytes.clone();
    downgraded[8..12].copy_from_slice(&1_u32.to_le_bytes());
    assert!(NormalizedPack::from_bytes(&downgraded).is_err());
    let mut relabelled = bytes;
    let token = quantity.canonical_bytes();
    let position = relabelled
        .windows(token.len())
        .position(|window| window == token)
        .unwrap();
    relabelled[position + 9] = 2; // cyclic -> angular, without changing value or claim id
    assert!(NormalizedPack::from_bytes(&relabelled).is_err());

    let legacy = sample_pack();
    assert_eq!(legacy.schema_version(), 1);
    let mut upgraded = legacy.to_bytes();
    upgraded[8..12].copy_from_slice(&2_u32.to_le_bytes());
    assert!(NormalizedPack::from_bytes(&upgraded).is_err());
}

#[test]
fn g3_typed_axis_pack_roundtrips_material_interface_and_joint_queries() {
    use fs_matdb::{MatDbError, QueryPoint, SelectionPolicy};
    use fs_qty::{
        QuantitySpec,
        semantic::{FrequencyConvention as F, QuantityKind as K, SemanticType, ValueForm},
    };
    let cyclic = QuantitySpec::semantic(SemanticType::new(
        K::Frequency(F::Cyclic),
        ValueForm::Static,
    ));
    let angular = QuantitySpec::semantic(SemanticType::new(
        K::Frequency(F::Angular),
        ValueForm::Static,
    ));
    let (original, observation, _) = sample_claims();
    let mut claims = ClaimSet::new();
    claims
        .register_observation(original.observation(observation).unwrap().clone())
        .unwrap();
    let mut claim = original.claims_for("density")[0].1.clone();
    claim.validity = ValidityDomain::unconstrained().with_quantity("frequency", cyclic, 40.0, 60.0);
    let second = PropertyClaim {
        key: PropertyKey::new("young_modulus", Dims([-1, 1, -2, 0, 0, 0])),
        value: PropertyValue::Scalar {
            value: 210.0e9,
            dims: Dims([-1, 1, -2, 0, 0, 0]),
        },
        validity: claim.validity.clone(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![observation],
        provenance: provenance(),
    };
    let id = claims.insert_claim(claim).unwrap();
    claims.insert_claim(second).unwrap();
    let pack = NormalizedPack::new(
        "typed-axis-fixture",
        "fixture-v3",
        hash_domain(SOURCE_DOMAIN, b"synthetic axes"),
        "synthetic fixture only",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    assert_eq!(pack.schema_version(), 3);
    let bytes = pack.to_bytes();
    let decoded = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(decoded.to_bytes(), bytes);
    assert_eq!(
        decoded
            .claims()
            .claim(id)
            .unwrap()
            .validity
            .axis_quantities()
            .get("frequency"),
        Some(&cyclic)
    );
    let state = fs_matdb::MaterialStateId {
        chemistry: "fixture".into(),
        phase: "fixture".into(),
        process: "fixture".into(),
        revision: 0,
    };
    let material =
        fs_matdb::NormalizedMaterialCardPack::new(state.clone(), decoded.clone()).unwrap();
    let material = fs_matdb::NormalizedMaterialCardPack::from_bytes_verified(
        material.content_hash(),
        &material.to_bytes(),
    )
    .unwrap();
    let surface = fs_matdb::SurfaceSpec {
        material: state,
        texture_frame: "fixture".into(),
    };
    let interface = fs_matdb::NormalizedInterfacePack::new(
        surface.clone(),
        surface,
        fs_matdb::SystemContext {
            medium: "fixture".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "fixture".into(),
        },
        decoded.clone(),
    )
    .unwrap();
    let interface = fs_matdb::NormalizedInterfacePack::from_bytes_verified(
        interface.content_hash(),
        &interface.to_bytes(),
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("frequency", cyclic, 50.0)
        .unwrap();
    let wrong = QueryPoint::new()
        .with_quantity("frequency", angular, 50.0)
        .unwrap();
    for claims in [material.card().claims(), interface.card().claims()] {
        let answer = claims
            .query("density", &point, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert_eq!(answer.evidence.value.value, 7850.0);
        claims.verify_receipt(&answer.receipt).unwrap();
        assert!(matches!(
            claims.query("density", &wrong, SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::AxisQuantityMismatch { .. })
        ));
    }
    let joint = decoded
        .query_joint(
            &["density", "young_modulus"],
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(joint.receipt.schema_version, 2);
    decoded.verify_joint_receipt(&joint.receipt).unwrap();
    let mut tampered = joint.receipt.clone();
    tampered.axis_quantities.insert("frequency".into(), angular);
    assert_ne!(tampered.content_hash(), joint.receipt.content_hash());
    assert!(decoded.verify_joint_receipt(&tampered).is_err());
    tampered = joint.receipt.clone();
    tampered.axis_quantities.insert("orphan".into(), cyclic);
    assert!(decoded.verify_joint_receipt(&tampered).is_err());
    for version in [1_u32, 2] {
        let mut downgraded = bytes.clone();
        downgraded[8..12].copy_from_slice(&version.to_le_bytes());
        assert!(NormalizedPack::from_bytes(&downgraded).is_err());
    }
    let mut relabelled = bytes;
    let token = cyclic.canonical_bytes();
    let position = relabelled
        .windows(token.len())
        .position(|window| window == token)
        .unwrap();
    relabelled[position..position + token.len()].copy_from_slice(&angular.canonical_bytes());
    assert!(NormalizedPack::from_bytes(&relabelled).is_err());
}

#[test]
fn construction_and_decoding_canonicalize_permutable_collections() {
    let first = sample_pack();
    let (claims, observation, members) = sample_claims();
    let second = NormalizedPack::new(
        "fixture-alloy-x",
        "frankensim-matdb-pack-compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source-envelope"),
        "CC-BY-4.0: redistribution permitted with attribution",
        claims,
        vec![JointStatistics::new(
            observation,
            "coupon-joint-density-modulus",
            members,
            vec![4.0, 1.0, 9.0],
            Some(vec![1.0, 1.0 / 6.0, 1.0]),
        )],
        first.normalizations().iter().rev().cloned().collect(),
    )
    .expect("permuted receipts canonicalize");

    assert_eq!(first.to_bytes(), second.to_bytes());
    assert_eq!(first.content_hash(), second.content_hash());
}

#[test]
fn g3_elastic_components_preserve_context_through_pack_query_and_replay() {
    use fs_blake3::ContentHash;
    use fs_matdb::{
        ElasticTensorBasis, ElasticTensorComponent, ElasticTensorNotation, ElasticTensorOrder,
        ElasticTensorSymmetry, MatDbError, QueryPoint, SelectionPolicy,
    };
    let basis = ElasticTensorBasis {
        notation: ElasticTensorNotation::Engineering,
        order: ElasticTensorOrder::XxYyZzXyYzZx,
        frame: ContentHash([1; 32]),
    };
    let descriptor = |basis, tensor, row, column| {
        ElasticTensorComponent::new(
            basis,
            ElasticTensorSymmetry::MajorMinor,
            ContentHash([tensor; 32]),
            row,
            column,
        )
        .unwrap()
    };
    let base = descriptor(basis, 2, 0, 3);
    let contexts = [
        base,
        descriptor(
            ElasticTensorBasis {
                frame: ContentHash([3; 32]),
                ..basis
            },
            2,
            0,
            3,
        ),
        descriptor(
            ElasticTensorBasis {
                notation: ElasticTensorNotation::Tensor,
                ..basis
            },
            2,
            0,
            3,
        ),
        descriptor(
            ElasticTensorBasis {
                order: ElasticTensorOrder::XxYyZzYzZxXy,
                ..basis
            },
            2,
            0,
            3,
        ),
        descriptor(basis, 3, 0, 3),
        descriptor(basis, 2, 0, 4),
    ];
    let bare = PropertyKey::new("coupling", fs_qty::Pressure::DIMS);
    assert!(
        PropertyKey::new("bad", Dims::NONE)
            .with_elastic_component(base)
            .is_err()
    );
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic tensor specimen".into(),
            method: "synthetic software fixture".into(),
            artifact: ContentHash([8; 32]),
            caveats: "no empirical or calibration claim".into(),
            provenance: provenance(),
        })
        .unwrap();
    let mut keys = Vec::new();
    let mut ids = std::collections::BTreeSet::new();
    for context in contexts {
        let key = bare.clone().with_elastic_component(context).unwrap();
        let id = claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value: -4.0,
                    dims: key.dims(),
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 300.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::TabulatedOnly,
                observations: vec![observation],
                provenance: provenance(),
            })
            .unwrap();
        assert!(
            ids.insert(id),
            "every context field binds the coefficient identity"
        );
        keys.push((key, id));
    }
    let pack = NormalizedPack::new(
        "synthetic-tensor",
        "fixture-v5",
        ContentHash([4; 32]),
        "synthetic redistribution permitted",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    assert_eq!(pack.schema_version(), 5);
    let bytes = pack.to_bytes();
    let decoded = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(decoded.to_bytes(), bytes);
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    assert!(matches!(
        decoded
            .claims()
            .query("coupling", &point, SelectionPolicy::SingleClaimOnly),
        Err(MatDbError::TensorContextMismatch { .. })
    ));
    assert!(
        decoded
            .claims()
            .query_typed(&bare, &point, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    for (key, id) in &keys {
        for answer in [
            decoded
                .claims()
                .query_typed(key, &point, SelectionPolicy::SingleClaimOnly)
                .unwrap(),
            decoded
                .claims()
                .query_pinned_typed(key, &point, *id)
                .unwrap(),
        ] {
            assert_eq!(answer.receipt.selected, *id);
            assert_eq!(answer.evidence.value.value, -4.0);
            assert_eq!(
                decoded.claims().claim(*id).unwrap().key.elastic_component(),
                key.elastic_component()
            );
            decoded.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    assert!(
        decoded
            .claims()
            .query_pinned_typed(&keys[0].0, &point, keys[1].1)
            .is_err()
    );
    assert!(
        decoded
            .claims()
            .query_typed(
                &keys[0].0,
                &QueryPoint::new().with("T", 301.0).unwrap(),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    let mut downgraded = bytes.clone();
    downgraded[8..12].copy_from_slice(&4u32.to_le_bytes());
    assert!(NormalizedPack::from_bytes(&downgraded).is_err());
    let encoded = base.canonical_bytes();
    let offset = bytes
        .windows(encoded.len())
        .position(|part| part == encoded)
        .unwrap();
    let mut tampered = bytes.clone();
    tampered[offset + 5] ^= 1;
    assert!(NormalizedPack::from_bytes(&tampered).is_err());
    for (index, value) in [(0, 3), (1, 2), (2, 1), (3, 6), (4, 6)] {
        let mut malformed = encoded;
        malformed[index] = value;
        assert!(ElasticTensorComponent::from_canonical_bytes(&malformed).is_err());
    }
    let mut zero_frame = encoded;
    zero_frame[5..37].fill(0);
    assert!(ElasticTensorComponent::from_canonical_bytes(&zero_frame).is_err());
}

#[test]
fn g3_hardness_tests_select_exact_apparatus_through_material_interface_and_joint_replay() {
    use fs_matdb::{
        HardnessLoadStep, HardnessTestContext, MatDbError, QueryPoint, SelectionPolicy,
    };
    use fs_qty::{
        QtyAny, QuantitySpec,
        semantic::{HardnessScale, QuantityKind, SemanticType, ValueForm},
    };
    let scale = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Hardness(HardnessScale::Vickers),
        ValueForm::Static,
    ));
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let force = |v| QtyAny::new(v, Dims([1, 1, -2, 0, 0, 0]));
    let time = |v| QtyAny::new(v, Dims([0, 0, 1, 0, 0, 0]));
    let (source, observation, _) = sample_claims();
    let mut claims = ClaimSet::new();
    claims
        .register_observation(source.observation(observation).unwrap().clone())
        .unwrap();
    let mut other_specimen = source.observation(observation).unwrap().clone();
    other_specimen.specimen = "different synthetic specimen/process".into();
    let other_observation = claims.register_observation(other_specimen).unwrap();
    let context = |indenter: &str, load, dwell, protocol: &str, observation| {
        HardnessTestContext::new(
            indenter,
            vec![HardnessLoadStep::new(force(load), time(dwell)).unwrap()],
            protocol,
            observation,
        )
        .unwrap()
    };
    let requested = context(
        "synthetic diamond pyramid geometry A",
        20.0,
        10.0,
        "synthetic protocol revision 1",
        observation,
    );
    let bare_key = PropertyKey::with_quantity("hardness", scale);
    let key = bare_key
        .clone()
        .with_hardness_test(requested.clone())
        .unwrap();
    let domain = ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 300.0);
    let claim = PropertyClaim {
        key: key.clone(),
        value: PropertyValue::Scalar {
            value: 210.0,
            dims: Dims::NONE,
        },
        validity: domain.clone(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![observation],
        provenance: provenance(),
    };
    let mut unlinked = claim.clone();
    unlinked.observations.clear();
    let before = claims.clone();
    assert!(matches!(
        claims.insert_claim(unlinked),
        Err(MatDbError::InvalidHardnessContext { .. })
    ));
    assert_eq!(claims, before);
    let selected = claims.insert_claim(claim.clone()).unwrap();
    let alternatives = [
        context(
            "synthetic diamond pyramid geometry B",
            20.0,
            10.0,
            "synthetic protocol revision 1",
            observation,
        ),
        context(
            "synthetic diamond pyramid geometry A",
            40.0,
            10.0,
            "synthetic protocol revision 1",
            observation,
        ),
        context(
            "synthetic diamond pyramid geometry A",
            20.0,
            15.0,
            "synthetic protocol revision 1",
            observation,
        ),
        context(
            "synthetic diamond pyramid geometry A",
            20.0,
            10.0,
            "synthetic protocol revision 2",
            observation,
        ),
        context(
            "synthetic diamond pyramid geometry A",
            20.0,
            10.0,
            "synthetic protocol revision 1",
            other_observation,
        ),
        HardnessTestContext::new(
            requested.indenter(),
            vec![
                HardnessLoadStep::new(force(2.0), time(5.0)).unwrap(),
                requested.loading()[0].clone(),
            ],
            requested.protocol(),
            observation,
        )
        .unwrap(),
    ];
    let mut alternate_ids = Vec::new();
    for test in &alternatives {
        let mut different = claim.clone();
        different.key = bare_key.clone().with_hardness_test(test.clone()).unwrap();
        different.observations = vec![test.observation()];
        let id = claims.insert_claim(different).unwrap();
        assert_ne!(id, selected);
        assert!(!alternate_ids.contains(&id));
        alternate_ids.push(id);
    }
    // Incomplete old data remains representable, but never supplies the
    // conditions omitted by its source or by a caller.
    let mut legacy = claim.clone();
    legacy.key = bare_key.clone();
    claims.insert_claim(legacy).unwrap();
    let mut density = source.claims_for("density")[0].1.clone();
    density.validity = domain;
    let density_key = density.key.clone();
    claims.insert_claim(density).unwrap();
    let pack = NormalizedPack::new(
        "synthetic-hardness",
        "fixture-v4",
        hash_domain(SOURCE_DOMAIN, b"synthetic hardness only"),
        "synthetic fixture only",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    assert_eq!(
        pack.schema_version(),
        fs_matdb::MATDB_HARDNESS_PACK_SCHEMA_VERSION
    );
    let bytes = pack.to_bytes();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claim(selected).unwrap().key, key);
    let state = fs_matdb::MaterialStateId {
        chemistry: "fixture".into(),
        phase: "fixture".into(),
        process: "fixture".into(),
        revision: 0,
    };
    let material = fs_matdb::NormalizedMaterialCardPack::new(state.clone(), pack.clone()).unwrap();
    let material = fs_matdb::NormalizedMaterialCardPack::from_bytes_verified(
        material.content_hash(),
        &material.to_bytes(),
    )
    .unwrap();
    let surface = fs_matdb::SurfaceSpec {
        material: state,
        texture_frame: "fixture".into(),
    };
    let interface = fs_matdb::NormalizedInterfacePack::new(
        surface.clone(),
        surface,
        fs_matdb::SystemContext {
            medium: "fixture".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "fixture".into(),
        },
        pack.clone(),
    )
    .unwrap();
    let interface = fs_matdb::NormalizedInterfacePack::from_bytes_verified(
        interface.content_hash(),
        &interface.to_bytes(),
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", temperature, 300.0)
        .unwrap();
    let absent = bare_key
        .clone()
        .with_hardness_test(context(
            requested.indenter(),
            21.0,
            10.0,
            requested.protocol(),
            observation,
        ))
        .unwrap();
    for claims in [material.card().claims(), interface.card().claims()] {
        for policy in [
            SelectionPolicy::SingleClaimOnly,
            SelectionPolicy::PreferObservationBacked,
        ] {
            let answer = claims.query_typed(&key, &point, policy).unwrap();
            assert_eq!(answer.evidence.value.value, 210.0);
            assert_eq!(answer.receipt.selected, selected);
            assert_eq!(answer.receipt.in_domain, vec![selected]);
            assert_eq!(answer.receipt.considered.len(), 8);
            let receipt = fs_matdb::PropertyUsageReceipt::from_bytes_verified(
                &answer.receipt.to_bytes().unwrap(),
                answer.receipt.content_hash(),
            )
            .unwrap();
            claims.verify_receipt(&receipt).unwrap();
            assert!(matches!(
                claims.query_typed(&bare_key, &point, policy),
                Err(MatDbError::MissingHardnessContext { .. })
            ));
            assert!(matches!(
                claims.query_typed(&absent, &point, policy),
                Err(MatDbError::HardnessContextMismatch { .. })
            ));
            let mut changed = receipt.clone();
            changed.selected = alternate_ids[0];
            assert!(claims.verify_receipt(&changed).is_err());
        }
        assert!(claims.query_pinned_typed(&key, &point, selected).is_ok());
        assert!(
            claims
                .query_pinned_typed(&key, &point, alternate_ids[0])
                .is_err()
        );
        assert!(
            claims
                .query_pinned_typed(&bare_key, &point, selected)
                .is_err()
        );
        for (test, id) in alternatives.iter().zip(&alternate_ids) {
            let different = bare_key.clone().with_hardness_test(test.clone()).unwrap();
            assert_eq!(
                claims
                    .query_typed(&different, &point, SelectionPolicy::SingleClaimOnly)
                    .unwrap()
                    .receipt
                    .selected,
                *id
            );
        }
        let outside = QueryPoint::new()
            .with_quantity("T", temperature, 301.0)
            .unwrap();
        assert!(
            claims
                .query_typed(&key, &outside, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
    }
    let joint = pack
        .query_joint_typed(
            &[key, density_key],
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    pack.verify_joint_receipt(&joint.receipt).unwrap();
    let mut tampered = joint.receipt.clone();
    tampered.selected[0] = alternate_ids[0];
    assert!(pack.verify_joint_receipt(&tampered).is_err());
    for version in [1_u32, 2, 3] {
        let mut relabelled = bytes.clone();
        relabelled[8..12].copy_from_slice(&version.to_le_bytes());
        assert!(NormalizedPack::from_bytes(&relabelled).is_err());
    }
    let mut changed = bytes;
    let position = changed
        .windows(requested.indenter().len())
        .position(|w| w == requested.indenter().as_bytes())
        .unwrap();
    changed[position] = b'X';
    assert!(NormalizedPack::from_bytes(&changed).is_err());
}

#[test]
fn joint_blocks_name_curve_components_and_allow_disjoint_groups() {
    let (claims, observation, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let modulus = claims.claims_for("young_modulus")[0].0;
    let scalar = vec![StatisticMember::scalar(density)];
    let mut curve = vec![
        StatisticMember::curve_ordinate(modulus, 0),
        StatisticMember::curve_ordinate(modulus, 1),
    ];
    curve.sort_unstable();
    let pack = NormalizedPack::new(
        "component-fixture",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"component-source"),
        "redistribution permitted",
        claims,
        vec![
            JointStatistics::new(
                observation,
                "modulus-knots",
                curve,
                vec![9.0, 1.0, 16.0],
                Some(vec![1.0, 1.0 / 12.0, 1.0]),
            ),
            JointStatistics::new(
                observation,
                "density-only",
                scalar,
                vec![4.0],
                Some(vec![1.0]),
            ),
        ],
        Vec::new(),
    )
    .expect("disjoint named blocks for one observation are admissible");

    assert_eq!(pack.joint_statistics().len(), 2);
    assert_eq!(pack.joint_statistics()[0].block_id(), "density-only");
    assert_eq!(pack.joint_statistics()[1].members().len(), 2);
    assert_eq!(
        NormalizedPack::from_bytes(&pack.to_bytes())
            .expect("component pack round-trip")
            .to_bytes(),
        pack.to_bytes()
    );
}

#[test]
fn joint_block_and_member_duplicates_refuse() {
    let (claims, observation, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let member = StatisticMember::scalar(density);
    let block = JointStatistics::new(observation, "duplicate", vec![member], vec![1.0], None);
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims.clone(),
            vec![block.clone(), block],
            Vec::new(),
        ),
        Err(PackError::InvalidField {
            field: "joint_statistics",
            ..
        })
    ));

    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            vec![JointStatistics::new(
                observation,
                "duplicate-members",
                vec![member, member],
                vec![1.0, 0.0, 1.0],
                None,
            )],
            Vec::new(),
        ),
        Err(PackError::InvalidField {
            field: "joint_statistics.members",
            ..
        })
    ));

    let (claims, observation, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let member = StatisticMember::scalar(density);
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            vec![
                JointStatistics::new(observation, "first", vec![member], vec![1.0], None,),
                JointStatistics::new(observation, "second", vec![member], vec![2.0], None,),
            ],
            Vec::new(),
        ),
        Err(PackError::InvalidField {
            field: "joint_statistics.members",
            ..
        })
    ));
}

#[test]
fn covariance_and_correlation_refuse_malformed_or_non_psd_blocks() {
    let cases = [
        (vec![1.0, 0.0], None, "lower triangle"),
        (vec![-1.0, 0.0, 1.0], None, "negative"),
        (vec![1.0, 2.0, 1.0], None, "outside [-1,1]"),
        (vec![0.0, 1.0e-15, 0.0], None, "zero variance"),
        (vec![0.0, 1.0e-15, 1.0], None, "zero variance"),
        (vec![1.0e-20, 2.0e-20, 1.0e-20], None, "outside [-1,1]"),
        (
            vec![1.0, 0.0, 1.0],
            Some(vec![0.9, 0.0, 1.0]),
            "exactly 1.0",
        ),
        (
            vec![1.0, 0.0, 1.0],
            Some(vec![1.0, 1.1, 1.0]),
            "outside [-1,1]",
        ),
        (
            vec![1.0, 0.0, 1.0],
            Some(vec![1.0, 1.0e-15, 1.0]),
            "inconsistent",
        ),
    ];

    for (covariance, correlation, expected) in cases {
        let (claims, observation, members) = sample_claims();
        let error = NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            vec![JointStatistics::new(
                observation,
                "matrix-refusal",
                members,
                covariance,
                correlation,
            )],
            Vec::new(),
        )
        .expect_err("invalid matrix must refuse");
        assert!(error.to_string().contains(expected), "{error}");
    }
}

#[test]
fn interval_psd_gate_admits_a_well_conditioned_three_member_block() {
    let (claims, observation, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let modulus = claims.claims_for("young_modulus")[0].0;
    let mut members = vec![
        StatisticMember::scalar(density),
        StatisticMember::curve_ordinate(modulus, 0),
        StatisticMember::curve_ordinate(modulus, 1),
    ];
    members.sort_unstable();

    NormalizedPack::new(
        "positive-definite-regression",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"positive-definite-source"),
        "redistribution permitted",
        claims,
        vec![JointStatistics::new(
            observation,
            "well-conditioned",
            members,
            vec![4.0, 1.0, 9.0, 0.5, 1.5, 16.0],
            None,
        )],
        Vec::new(),
    )
    .expect("strictly diagonally dominant covariance must pass interval admission");
}

#[test]
fn interval_psd_gate_refuses_a_rounded_zero_false_certificate() {
    let (claims, observation, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let modulus = claims.claims_for("young_modulus")[0].0;
    let mut members = vec![
        StatisticMember::scalar(density),
        StatisticMember::curve_ordinate(modulus, 0),
        StatisticMember::curve_ordinate(modulus, 1),
    ];
    members.sort_unstable();

    let error = NormalizedPack::new(
        "false-certificate-regression",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"false-certificate-source"),
        "redistribution permitted",
        claims,
        vec![JointStatistics::new(
            observation,
            "rounded-zero-indefinite",
            members,
            vec![
                1.0,
                -0.320_643_490_748_954_5,
                1.0,
                0.934_029_952_841_983,
                0.038_844_170_036_952_126,
                1.0,
            ],
            None,
        )],
        Vec::new(),
    )
    .expect_err("an exact-negative determinant hidden by a rounded zero must refuse");
    assert!(error.to_string().contains("rounding-ambiguous"), "{error}");
}

#[test]
fn aggregate_psd_work_budget_refuses_cpu_amplification() {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "work-budget fixture".to_string(),
            method: "synthetic covariance workload".to_string(),
            artifact: hash_domain(SOURCE_DOMAIN, b"work-budget-observation"),
            caveats: String::new(),
            provenance: provenance(),
        })
        .expect("work-budget observation");
    let knots: Vec<_> = (0..256)
        .map(|index| (f64::from(index), f64::from(index + 1)))
        .collect();
    let curve = claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("work_budget_curve", Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: "sample".to_string(),
                abscissa_dims: Dims::NONE,
                knots,
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::TabulatedOnly,
            observations: vec![observation],
            provenance: provenance(),
        })
        .expect("work-budget curve");
    let members: Vec<_> = (0..256)
        .map(|knot| StatisticMember::curve_ordinate(curve, knot))
        .collect();
    let mut identity = Vec::with_capacity(256 * 257 / 2);
    for row in 0..256 {
        for column in 0..=row {
            identity.push(if row == column { 1.0 } else { 0.0 });
        }
    }
    let blocks = (0..9)
        .map(|index| {
            JointStatistics::new(
                observation,
                format!("block-{index}"),
                members.clone(),
                identity.clone(),
                None,
            )
        })
        .collect();

    assert!(matches!(
        NormalizedPack::new(
            "work-budget",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"work-budget-source"),
            "redistribution permitted",
            claims,
            blocks,
            Vec::new(),
        ),
        Err(PackError::ResourceLimit {
            resource: "psd_cubic_work",
            ..
        })
    ));
}

#[test]
fn joint_statistics_must_reference_cited_claims_and_observations() {
    let (mut claims, observation, members) = sample_claims();
    let uncited_observation = claims
        .register_observation(ObservationDataset {
            specimen: "different coupon".to_string(),
            method: "independent method".to_string(),
            artifact: hash_domain(SOURCE_DOMAIN, b"uncited-observation"),
            caveats: String::new(),
            provenance: provenance(),
        })
        .expect("second observation");
    let unknown_observation = ObservationId(hash_domain(SOURCE_DOMAIN, b"unknown-observation"));
    let error = NormalizedPack::new(
        "fixture",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source"),
        "redistribution permitted",
        claims.clone(),
        vec![JointStatistics::new(
            unknown_observation,
            "unknown-observation",
            members.clone(),
            vec![1.0, 0.0, 1.0],
            None,
        )],
        Vec::new(),
    )
    .expect_err("unknown observation must refuse");
    assert!(error.to_string().contains("unknown observation"), "{error}");

    let error = NormalizedPack::new(
        "fixture",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source"),
        "redistribution permitted",
        claims.clone(),
        vec![JointStatistics::new(
            uncited_observation,
            "uncited",
            members.clone(),
            vec![1.0, 0.0, 1.0],
            None,
        )],
        Vec::new(),
    )
    .expect_err("member must cite the joint observation");
    assert!(error.to_string().contains("does not cite"), "{error}");

    let unknown_claim = ClaimId(hash_domain(SOURCE_DOMAIN, b"unknown-claim"));
    let mut bad_members = vec![members[0], StatisticMember::scalar(unknown_claim)];
    bad_members.sort_unstable();
    let error = NormalizedPack::new(
        "fixture",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source"),
        "redistribution permitted",
        claims,
        vec![JointStatistics::new(
            observation,
            "unknown-claim",
            bad_members,
            vec![1.0, 0.0, 1.0],
            None,
        )],
        Vec::new(),
    )
    .expect_err("unknown claim must refuse");
    assert!(error.to_string().contains("unknown claim"), "{error}");
}

#[test]
fn canonical_profile_refuses_negative_zero_and_partial_frame_receipts() {
    let (mut claims, observation, _) = sample_claims();
    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("signed_zero_probe", density_dims),
            value: PropertyValue::Scalar {
                value: -0.0,
                dims: density_dims,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .expect("core layer preserves exact signed zero");
    let error = NormalizedPack::new(
        "fixture",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source"),
        "redistribution permitted",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .expect_err("portable profile must refuse negative zero");
    assert!(error.to_string().contains("negative zero"), "{error}");

    let (claims, _, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let error = NormalizedPack::new(
        "fixture",
        "compiler-v1",
        hash_domain(SOURCE_DOMAIN, b"source"),
        "redistribution permitted",
        claims,
        Vec::new(),
        vec![NormalizationReceipt::new(
            NormalizationTarget::ClaimValue(StatisticMember::scalar(density)),
            hash_domain(SOURCE_DOMAIN, b"7.85 g/cm3"),
            Dims([-3, 1, 0, 0, 0, 0]),
            1_000.0,
            0.0,
            "g/cm3",
            MATDB_PACK_TARGET_BASIS,
            Some("rolling-frame".to_string()),
            None,
        )],
    )
    .expect_err("a half-declared frame transform must refuse");
    assert!(error.to_string().contains("both be present"), "{error}");
}

#[test]
#[allow(clippy::too_many_lines)] // one refusal drill per linkage rule, in admission order
fn normalization_targets_are_linked_dimension_checked_and_unique() {
    let (claims, _, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let target = NormalizationTarget::ClaimValue(StatisticMember::scalar(density));
    let receipt = NormalizationReceipt::new(
        target.clone(),
        hash_domain(SOURCE_DOMAIN, b"7.85 g/cm3"),
        Dims::NONE,
        1_000.0,
        0.0,
        "g/cm3",
        MATDB_PACK_TARGET_BASIS,
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims.clone(),
            Vec::new(),
            vec![receipt],
        ),
        Err(PackError::InvalidField {
            field: "normalization.dims",
            ..
        })
    ));

    let valid = NormalizationReceipt::new(
        target,
        hash_domain(SOURCE_DOMAIN, b"7.85 g/cm3"),
        Dims([-3, 1, 0, 0, 0, 0]),
        1_000.0,
        0.0,
        "g/cm3",
        MATDB_PACK_TARGET_BASIS,
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            Vec::new(),
            vec![valid.clone(), valid],
        ),
        Err(PackError::InvalidField {
            field: "normalizations",
            ..
        })
    ));

    let pack = sample_pack();
    let density = pack.claims().claims_for("density")[0].0;
    let translated_uncertainty = NormalizationReceipt::new(
        NormalizationTarget::ClaimUncertainty { claim: density },
        hash_domain(SOURCE_DOMAIN, b"translated uncertainty"),
        Dims([-3, 1, 0, 0, 0, 0]),
        1.0,
        1.0,
        "kg/m3",
        MATDB_PACK_TARGET_BASIS,
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            pack.claims().clone(),
            Vec::new(),
            vec![translated_uncertainty],
        ),
        Err(PackError::InvalidField {
            field: "normalization.offset",
            ..
        })
    ));

    let joint = pack.joint_statistics()[0].clone();
    let translated_covariance = NormalizationReceipt::new(
        NormalizationTarget::JointCovariance {
            observation: joint.observation(),
            block_id: joint.block_id().to_string(),
            row: 1,
            column: 0,
        },
        hash_domain(SOURCE_DOMAIN, b"translated covariance"),
        Dims([-4, 2, -2, 0, 0, 0]),
        1.0,
        1.0,
        "coherent source covariance",
        MATDB_PACK_TARGET_BASIS,
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            pack.claims().clone(),
            vec![joint],
            vec![translated_covariance],
        ),
        Err(PackError::InvalidField {
            field: "normalization.offset",
            ..
        })
    ));

    let negative_uncertainty_scale = NormalizationReceipt::new(
        NormalizationTarget::ClaimUncertainty { claim: density },
        hash_domain(SOURCE_DOMAIN, b"negative uncertainty scale"),
        Dims([-3, 1, 0, 0, 0, 0]),
        -1.0,
        0.0,
        "kg/m3",
        MATDB_PACK_TARGET_BASIS,
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            pack.claims().clone(),
            Vec::new(),
            vec![negative_uncertainty_scale],
        ),
        Err(PackError::InvalidField {
            field: "normalization.scale",
            ..
        })
    ));

    let (claims, observation, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let variance_block = JointStatistics::new(
        observation,
        "density-variance",
        vec![StatisticMember::scalar(density)],
        vec![4.0],
        None,
    );
    let negative_variance_scale = NormalizationReceipt::new(
        NormalizationTarget::JointCovariance {
            observation,
            block_id: "density-variance".to_string(),
            row: 0,
            column: 0,
        },
        hash_domain(SOURCE_DOMAIN, b"negative variance scale"),
        Dims([-6, 2, 0, 0, 0, 0]),
        -1.0,
        0.0,
        "source variance basis",
        MATDB_PACK_TARGET_BASIS,
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            vec![variance_block],
            vec![negative_variance_scale],
        ),
        Err(PackError::InvalidField {
            field: "normalization.scale",
            ..
        })
    ));

    let (claims, _, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let validity_receipts = [
        (fs_matdb::ValidityBoundSide::Lower, Dims([0, 0, 0, 1, 0, 0])),
        (fs_matdb::ValidityBoundSide::Upper, Dims::NONE),
    ]
    .map(|(side, dims)| {
        NormalizationReceipt::new(
            NormalizationTarget::ValidityBound {
                claim: density,
                axis: "temperature".to_string(),
                side,
            },
            hash_domain(SOURCE_DOMAIN, b"contradictory validity units"),
            dims,
            1.0,
            0.0,
            "source temperature",
            MATDB_PACK_TARGET_BASIS,
            None,
            None,
        )
    });
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            Vec::new(),
            validity_receipts.into(),
        ),
        Err(PackError::InvalidField {
            field: "normalization.dims",
            ..
        })
    ));
}

#[test]
fn pack_metadata_and_transform_policy_refuse_structured_errors() {
    let (claims, _, _) = sample_claims();
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "   ",
            claims.clone(),
            Vec::new(),
            Vec::new(),
        ),
        Err(PackError::InvalidField {
            field: "redistribution_terms",
            ..
        })
    ));

    let density = claims.claims_for("density")[0].0;
    let receipt = NormalizationReceipt::new(
        NormalizationTarget::ClaimValue(StatisticMember::scalar(density)),
        hash_domain(SOURCE_DOMAIN, b"7.85 g/cm3"),
        Dims([-3, 1, 0, 0, 0, 0]),
        0.0,
        0.0,
        "g/cm3",
        MATDB_PACK_TARGET_BASIS,
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            Vec::new(),
            vec![receipt],
        ),
        Err(PackError::InvalidField {
            field: "normalization.scale",
            ..
        })
    ));

    let (claims, _, _) = sample_claims();
    let density = claims.claims_for("density")[0].0;
    let receipt = NormalizationReceipt::new(
        NormalizationTarget::ClaimValue(StatisticMember::scalar(density)),
        hash_domain(SOURCE_DOMAIN, b"7.85 g/cm3"),
        Dims([-3, 1, 0, 0, 0, 0]),
        1_000.0,
        0.0,
        "g/cm3",
        "unspecified target basis",
        None,
        None,
    );
    assert!(matches!(
        NormalizedPack::new(
            "fixture",
            "compiler-v1",
            hash_domain(SOURCE_DOMAIN, b"source"),
            "redistribution permitted",
            claims,
            Vec::new(),
            vec![receipt],
        ),
        Err(PackError::InvalidField {
            field: "normalization.target_basis",
            ..
        })
    ));
}

#[test]
fn wire_refuses_truncation_trailing_bytes_and_semantic_id_tamper() {
    let pack = sample_pack();
    let expected_pack_hash = pack.content_hash();
    let bytes = pack.to_bytes();

    let truncated = &bytes[..bytes.len() - 1];
    assert!(
        NormalizedPack::from_bytes(truncated)
            .expect_err("truncation must refuse")
            .to_string()
            .contains("truncated")
    );

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(
        NormalizedPack::from_bytes(&trailing)
            .expect_err("trailing byte must refuse")
            .to_string()
            .contains("trailing")
    );

    let mut top_level_tamper = bytes.clone();
    top_level_tamper[16] ^= 0x01; // first UTF-8 byte of the length-framed pack id
    assert!(NormalizedPack::from_bytes(&top_level_tamper).is_ok());
    assert!(matches!(
        NormalizedPack::from_bytes_verified(expected_pack_hash, &top_level_tamper),
        Err(PackError::IdentityMismatch { kind: "pack", .. })
    ));

    let mut tampered = bytes;
    // Header + schema + four length-framed envelope fields.  The next u32 is
    // observation count and the following 32 bytes are its semantic id.
    let mut cursor = 8 + 4;
    for field in 0..4 {
        if field == 2 {
            cursor += 32; // source_artifact is fixed-width, not a string
        } else {
            let length = u32::from_le_bytes(
                tampered[cursor..cursor + 4]
                    .try_into()
                    .expect("fixture string length"),
            ) as usize;
            cursor += 4 + length;
        }
    }
    cursor += 4; // observation count
    tampered[cursor] ^= 0x01;
    assert!(matches!(
        NormalizedPack::from_bytes(&tampered),
        Err(PackError::IdentityMismatch {
            kind: "observation",
            ..
        })
    ));
}

#[test]
fn decoder_preflights_untrusted_counts_before_semantic_allocation() {
    let mut blocks = bare_pack_prefix();
    push_u32(&mut blocks, 0); // observations
    push_u32(&mut blocks, 0); // claims
    push_u32(&mut blocks, 100_000); // joint blocks, no payload
    let error = NormalizedPack::from_bytes(&blocks).expect_err("missing block payload must refuse");
    assert!(
        error
            .to_string()
            .contains("truncated joint-statistics blocks"),
        "{error}"
    );

    let mut members = bare_pack_prefix();
    push_u32(&mut members, 0); // observations
    push_u32(&mut members, 0); // claims
    push_u32(&mut members, 1); // joint block
    members.extend_from_slice(&[0; 32]); // observation id
    push_string(&mut members, ""); // block id
    push_u32(&mut members, 256); // members, no payload
    members.push(0); // satisfy the outer block's conservative one-byte tail
    let error =
        NormalizedPack::from_bytes(&members).expect_err("missing member payload must refuse");
    assert!(
        error
            .to_string()
            .contains("truncated joint-statistics members"),
        "{error}"
    );

    let mut receipts = bare_pack_prefix();
    push_u32(&mut receipts, 0); // observations
    push_u32(&mut receipts, 0); // claims
    push_u32(&mut receipts, 0); // joint blocks
    push_u32(&mut receipts, 100_000); // receipts, no payload
    let error =
        NormalizedPack::from_bytes(&receipts).expect_err("missing receipt payload must refuse");
    assert!(
        error
            .to_string()
            .contains("truncated normalization receipts"),
        "{error}"
    );

    let mut observations = bare_pack_prefix();
    push_u32(&mut observations, 0); // observations
    push_u32(&mut observations, 1); // claims
    observations.extend_from_slice(&[0; 32]); // encoded claim id
    push_string(&mut observations, "preflight-property");
    observations.extend_from_slice(&[0; 6]); // key dims
    observations.push(0); // scalar value
    observations.extend_from_slice(&0.0f64.to_bits().to_le_bytes());
    observations.extend_from_slice(&[0; 6]); // scalar dims
    push_u32(&mut observations, 0); // validity bounds
    observations.push(0); // Unstated uncertainty
    observations.push(1); // ConstantWithinValidity
    push_u32(&mut observations, 100_000); // observation ids, no payload
    let error = NormalizedPack::from_bytes(&observations)
        .expect_err("missing claim-observation payload must refuse");
    assert!(
        error
            .to_string()
            .contains("truncated claim observation ids"),
        "{error}"
    );

    let mut curve = bare_pack_prefix();
    push_u32(&mut curve, 0); // observations
    push_u32(&mut curve, 1); // claims
    curve.extend_from_slice(&[0; 32]); // encoded claim id
    push_string(&mut curve, "preflight-curve");
    curve.extend_from_slice(&[0; 6]); // key dims
    curve.push(1); // curve value
    push_string(&mut curve, "sample");
    curve.extend_from_slice(&[0; 6]); // abscissa dims
    push_u32(&mut curve, 4_000_000); // knots, no payload
    let error =
        NormalizedPack::from_bytes(&curve).expect_err("missing curve-knot payload must refuse");
    assert!(
        error.to_string().contains("truncated curve knots"),
        "{error}"
    );
}
