//! G0/G3 coverage for the typed scenario-payload algebra and canonical V1
//! codec.  These tests intentionally exercise the closed payload/source/kind
//! tags rather than depending on the higher-level scenario IR integration.

use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::{
    AngleDomain, CompositionBasis, PhasorAmplitude, PhasorQty, QuantityKind, SemanticType,
    StrainBasis, StrainComponent, ValueForm,
};
use fs_qty::{Dims, QtyAny};
use fs_scenario::FrameId;
use fs_scenario::payload::{
    CharacteristicComponent, CharacteristicDirection, CharacteristicState, ComplexPhasorPayload,
    DistributionFamily, FieldTraceRef, MAX_PAYLOAD_ID_BYTES, MAX_PAYLOAD_ITEMS,
    MAX_PAYLOAD_WIRE_BYTES, OrientationParity, OutsideDomainPolicy, PAYLOAD_WIRE_VERSION, Payload,
    PayloadDecodeLimits, PayloadError, PayloadId, PayloadKind, PayloadMeta, PortRef,
    QuantityContract, ReferenceSemantics, SampleSource, ScalarPayload, SpeciesBundle, SpeciesValue,
    TableInterpolation, TensorPayload, VectorPayload, canonical_payload_byte_len,
    canonical_payload_bytes, decode_payload, decode_payload_with_limits,
    try_canonical_payload_bytes,
};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const AMOUNT: Dims = Dims([0, 0, 0, 0, 0, 1]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);

fn id(value: &str) -> PayloadId {
    PayloadId::new(value).expect("canonical fixture identifier")
}

fn dims_meta(dims: Dims, frame: u32) -> PayloadMeta {
    PayloadMeta::new(
        QuantityContract::Dimensions(dims),
        id("basis/world-cartesian"),
        FrameId(frame),
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .expect("valid dimensional metadata")
}

fn semantic_meta(kind: QuantityKind, form: ValueForm) -> PayloadMeta {
    PayloadMeta::new(
        QuantityContract::Semantic(SemanticType::new(kind, form)),
        id("basis/semantic-scalar"),
        FrameId(7),
        OrientationParity::Odd,
        ReferenceSemantics::ResetAtEvent(id("events/impact-1")),
    )
    .expect("valid semantic metadata")
}

fn heterogeneous_meta(frame: u32) -> PayloadMeta {
    PayloadMeta::new(
        QuantityContract::Heterogeneous,
        id("basis/characteristic"),
        FrameId(frame),
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .expect("valid heterogeneous metadata")
}

fn all_payload_variants() -> Vec<Payload> {
    let scalar = Payload::Scalar(
        ScalarPayload::new(
            semantic_meta(QuantityKind::Pressure, ValueForm::Static),
            SampleSource::fixed(QtyAny::new(101_325.0, PRESSURE)),
        )
        .expect("valid scalar"),
    );

    let vector = Payload::Vector(
        VectorPayload::new(
            dims_meta(VELOCITY, 2),
            SampleSource::table(
                vec![QtyAny::new(0.0, TIME), QtyAny::new(0.25, TIME)],
                vec![
                    vec![
                        QtyAny::new(1.0, VELOCITY),
                        QtyAny::new(2.0, VELOCITY),
                        QtyAny::new(3.0, VELOCITY),
                    ],
                    vec![
                        QtyAny::new(1.5, VELOCITY),
                        QtyAny::new(2.5, VELOCITY),
                        QtyAny::new(3.5, VELOCITY),
                    ],
                ],
                TableInterpolation::Linear,
                OutsideDomainPolicy::Refuse,
            )
            .expect("valid table"),
        )
        .expect("valid vector"),
    );

    let strain = SemanticType::new(
        QuantityKind::Strain {
            basis: StrainBasis::Engineering,
            component: StrainComponent::Shear,
        },
        ValueForm::Static,
    );
    let tensor = Payload::Tensor(
        TensorPayload::new(
            PayloadMeta::new(
                QuantityContract::Semantic(strain),
                id("basis/material"),
                FrameId(4),
                OrientationParity::Even,
                ReferenceSemantics::Continuous,
            )
            .expect("valid tensor metadata"),
            2,
            2,
            SampleSource::distribution(
                DistributionFamily::Normal,
                vec![
                    vec![
                        QtyAny::dimensionless(0.1),
                        QtyAny::dimensionless(0.2),
                        QtyAny::dimensionless(0.3),
                        QtyAny::dimensionless(0.4),
                    ],
                    vec![
                        QtyAny::dimensionless(0.01),
                        QtyAny::dimensionless(0.02),
                        QtyAny::dimensionless(0.03),
                        QtyAny::dimensionless(0.04),
                    ],
                ],
            )
            .expect("valid normal parameters"),
        )
        .expect("valid tensor"),
    );

    let phasor = PhasorQty::new(
        QtyAny::new(-12.5, PRESSURE),
        QtyAny::new(3.25, PRESSURE),
        QuantityKind::Pressure,
        PhasorAmplitude::Rms,
    )
    .expect("valid pressure phasor");
    let complex = Payload::ComplexPhasor(
        ComplexPhasorPayload::new(
            semantic_meta(QuantityKind::Pressure, ValueForm::Rms),
            SampleSource::fixed(phasor),
        )
        .expect("valid complex payload"),
    );

    let species_sample = vec![
        SpeciesValue::new(
            SpeciesId::new("CO2").expect("species id"),
            QtyAny::new(2.0, AMOUNT),
        ),
        SpeciesValue::new(
            SpeciesId::new("H2O").expect("species id"),
            QtyAny::new(3.0, AMOUNT),
        ),
    ];
    let species = Payload::SpeciesBundle(
        SpeciesBundle::new(dims_meta(AMOUNT, 0), SampleSource::fixed(species_sample))
            .expect("valid species bundle"),
    );

    let characteristic = Payload::CharacteristicState(
        CharacteristicState::new(
            heterogeneous_meta(3),
            vec![
                CharacteristicComponent::new(
                    id("pressure"),
                    CharacteristicDirection::Incoming,
                    QuantityContract::Dimensions(PRESSURE),
                )
                .expect("pressure characteristic"),
                CharacteristicComponent::new(
                    id("temperature"),
                    CharacteristicDirection::Stationary,
                    QuantityContract::Dimensions(TEMPERATURE),
                )
                .expect("temperature characteristic"),
                CharacteristicComponent::new(
                    id("velocity"),
                    CharacteristicDirection::Outgoing,
                    QuantityContract::Dimensions(VELOCITY),
                )
                .expect("velocity characteristic"),
            ],
            SampleSource::distribution(
                DistributionFamily::Empirical,
                vec![vec![
                    QtyAny::new(1.0, PRESSURE),
                    QtyAny::new(300.0, TEMPERATURE),
                    QtyAny::new(3.0, VELOCITY),
                ]],
            )
            .expect("valid empirical support"),
        )
        .expect("valid characteristic state"),
    );

    let field = Payload::FieldTraceRef(
        FieldTraceRef::new(
            dims_meta(PRESSURE, 8),
            id("blake3:abc123"),
            id("fields/wall-pressure"),
        )
        .expect("valid field ref"),
    );
    let port = Payload::PortRef(
        PortRef::new(
            dims_meta(VELOCITY, 9),
            id("components/pump-1"),
            id("ports/outlet"),
        )
        .expect("valid port ref"),
    );

    vec![
        scalar,
        vector,
        tensor,
        complex,
        species,
        characteristic,
        field,
        port,
    ]
}

#[test]
fn g0_all_payload_variants_round_trip_to_exact_canonical_bytes() {
    assert_eq!(PAYLOAD_WIRE_VERSION, 1);
    for payload in all_payload_variants() {
        let first = canonical_payload_bytes(&payload);
        let second = canonical_payload_bytes(&payload);
        assert_eq!(first, second, "encoder must be byte deterministic");
        assert_eq!(
            canonical_payload_byte_len(&payload).expect("exact byte plan"),
            first.len(),
            "counter and materializing encoders must traverse one grammar"
        );
        assert_eq!(
            try_canonical_payload_bytes(&payload).expect("fallible exact-reservation encode"),
            first,
            "fallible and compatibility encoders must be byte-identical"
        );

        let decoded = decode_payload(&first).expect("canonical bytes decode");
        assert_eq!(decoded, payload);
        assert_eq!(canonical_payload_bytes(&decoded), first);
    }
}

#[test]
fn g0_every_closed_semantic_kind_survives_the_codec() {
    let kinds = [
        QuantityKind::AbsoluteTemperature,
        QuantityKind::TemperatureDifference,
        QuantityKind::Angle(AngleDomain::Mechanical),
        QuantityKind::Angle(AngleDomain::Electrical),
        QuantityKind::AngularVelocity(AngleDomain::Mechanical),
        QuantityKind::AngularVelocity(AngleDomain::Electrical),
        QuantityKind::Torque,
        QuantityKind::Energy,
        QuantityKind::Pressure,
        QuantityKind::Stress,
        QuantityKind::Strain {
            basis: StrainBasis::Tensor,
            component: StrainComponent::Normal,
        },
        QuantityKind::Strain {
            basis: StrainBasis::Engineering,
            component: StrainComponent::Shear,
        },
        QuantityKind::Composition(CompositionBasis::MassFraction),
        QuantityKind::Composition(CompositionBasis::MoleFraction),
        QuantityKind::Composition(CompositionBasis::VolumeFraction),
        QuantityKind::Mass,
        QuantityKind::Amount,
        QuantityKind::MolarMass,
        QuantityKind::MassConcentration,
        QuantityKind::AmountConcentration,
        QuantityKind::Entropy,
        QuantityKind::HeatCapacity,
        QuantityKind::AcousticPressure,
        QuantityKind::AcousticPower,
    ];

    for kind in kinds {
        let payload = Payload::PortRef(
            PortRef::new(
                semantic_meta(kind, ValueForm::Static),
                id("component"),
                id("port"),
            )
            .expect("valid semantic port"),
        );
        let bytes = canonical_payload_bytes(&payload);
        assert_eq!(decode_payload(&bytes).expect("kind decodes"), payload);
    }
}

#[test]
fn g0_v1_payload_refuses_material_conventions_before_infallible_encoding() {
    use fs_qty::semantic::{
        AttenuationConvention, FrequencyConvention, HardnessScale, MoistureBasis,
    };
    for kind in [
        QuantityKind::Frequency(FrequencyConvention::Cyclic),
        QuantityKind::Frequency(FrequencyConvention::Angular),
        QuantityKind::Hardness(HardnessScale::RockwellC),
        QuantityKind::Moisture(MoistureBasis::DryMass),
        QuantityKind::Attenuation(AttenuationConvention::DecibelsPerMetre),
        QuantityKind::ThermalConductivity,
        QuantityKind::OpticalExtinctionCoefficient,
    ] {
        let semantic = SemanticType::new(kind, ValueForm::Static);
        assert!(matches!(
            PayloadMeta::new(
                QuantityContract::Semantic(semantic),
                id("basis/material"),
                FrameId(7),
                OrientationParity::Even,
                ReferenceSemantics::Continuous
            ),
            Err(PayloadError::InvalidSemanticContract { .. })
        ));
    }
}

#[test]
fn phasors_require_one_semantic_kind_and_peak_or_rms_convention() {
    for (form, amplitude) in [
        (ValueForm::Peak, PhasorAmplitude::Peak),
        (ValueForm::Rms, PhasorAmplitude::Rms),
    ] {
        let phasor = PhasorQty::new(
            QtyAny::new(2.0, PRESSURE),
            QtyAny::new(0.5, PRESSURE),
            QuantityKind::Pressure,
            amplitude,
        )
        .expect("phasor-capable pressure kind");
        let payload = Payload::ComplexPhasor(
            ComplexPhasorPayload::new(
                semantic_meta(QuantityKind::Pressure, form),
                SampleSource::fixed(phasor),
            )
            .expect("matching semantic phasor contract"),
        );
        let bytes = canonical_payload_bytes(&payload);
        assert_eq!(decode_payload(&bytes).expect("phasor decodes"), payload);
    }

    let pressure_peak = PhasorQty::new(
        QtyAny::new(2.0, PRESSURE),
        QtyAny::new(0.5, PRESSURE),
        QuantityKind::Pressure,
        PhasorAmplitude::Peak,
    )
    .expect("pressure peak phasor");
    assert_eq!(
        ComplexPhasorPayload::new(dims_meta(PRESSURE, 0), SampleSource::fixed(pressure_peak),),
        Err(PayloadError::PhasorSemanticContractRequired)
    );

    let pressure_rms = PhasorQty::new(
        QtyAny::new(2.0, PRESSURE),
        QtyAny::new(0.5, PRESSURE),
        QuantityKind::Pressure,
        PhasorAmplitude::Rms,
    )
    .expect("pressure RMS phasor");
    assert!(matches!(
        ComplexPhasorPayload::new(
            semantic_meta(QuantityKind::Pressure, ValueForm::Peak),
            SampleSource::fixed(pressure_rms),
        ),
        Err(PayloadError::PhasorContractMismatch { .. })
    ));
}

#[test]
fn source_policy_tags_round_trip_without_implied_evaluation_claims() {
    let periodic = Payload::Scalar(
        ScalarPayload::new(
            dims_meta(PRESSURE, 1),
            SampleSource::table(
                vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
                vec![QtyAny::new(10.0, PRESSURE), QtyAny::new(20.0, PRESSURE)],
                TableInterpolation::StepLeft,
                OutsideDomainPolicy::Periodic,
            )
            .expect("periodic table shape"),
        )
        .expect("periodic scalar payload"),
    );
    let uniform = Payload::Scalar(
        ScalarPayload::new(
            dims_meta(PRESSURE, 1),
            SampleSource::distribution(
                DistributionFamily::Uniform,
                vec![QtyAny::new(10.0, PRESSURE), QtyAny::new(20.0, PRESSURE)],
            )
            .expect("uniform arity"),
        )
        .expect("ordered uniform payload"),
    );

    for payload in [periodic, uniform] {
        let bytes = canonical_payload_bytes(&payload);
        assert_eq!(decode_payload(&bytes).expect("policy tags decode"), payload);
    }
}

#[test]
fn species_composition_validates_the_whole_six_base_bundle() {
    let meta = semantic_meta(
        QuantityKind::Composition(CompositionBasis::MassFraction),
        ValueForm::Static,
    );
    let sample = |water: f64| {
        vec![
            SpeciesValue::new(
                SpeciesId::new("CO2").expect("species"),
                QtyAny::dimensionless(0.4),
            ),
            SpeciesValue::new(
                SpeciesId::new("H2O").expect("species"),
                QtyAny::dimensionless(water),
            ),
        ]
    };
    let payload = Payload::SpeciesBundle(
        SpeciesBundle::new(meta.clone(), SampleSource::fixed(sample(0.6)))
            .expect("unit-sum composition"),
    );
    let bytes = canonical_payload_bytes(&payload);
    assert_eq!(
        decode_payload(&bytes).expect("composition decodes"),
        payload
    );

    assert!(matches!(
        SpeciesBundle::new(meta, SampleSource::fixed(sample(0.5))),
        Err(PayloadError::Semantic(_))
    ));
}

#[test]
fn constructors_refuse_noncanonical_or_incoherent_values() {
    assert!(matches!(
        PayloadId::new("bad id"),
        Err(PayloadError::InvalidIdentifier { .. })
    ));

    let wrong_dims = ScalarPayload::new(
        semantic_meta(QuantityKind::Pressure, ValueForm::Static),
        SampleSource::fixed(QtyAny::new(1.0, VELOCITY)),
    );
    assert!(matches!(
        wrong_dims,
        Err(PayloadError::DimensionMismatch { .. })
    ));

    let nonfinite = ScalarPayload::new(
        dims_meta(PRESSURE, 0),
        SampleSource::fixed(QtyAny::new(f64::INFINITY, PRESSURE)),
    );
    assert!(matches!(nonfinite, Err(PayloadError::NonFinite { .. })));

    assert!(matches!(
        SampleSource::<QtyAny>::table(
            vec![QtyAny::new(1.0, TIME), QtyAny::new(1.0, TIME)],
            vec![QtyAny::new(1.0, PRESSURE), QtyAny::new(2.0, PRESSURE)],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        ),
        Err(PayloadError::InvalidTableTime { index: 1, .. })
    ));

    let shape_drift = VectorPayload::new(
        dims_meta(VELOCITY, 0),
        SampleSource::table(
            vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
            vec![
                vec![QtyAny::new(1.0, VELOCITY)],
                vec![QtyAny::new(1.0, VELOCITY), QtyAny::new(2.0, VELOCITY)],
            ],
            TableInterpolation::StepLeft,
            OutsideDomainPolicy::Clamp,
        )
        .expect("time/value counts agree"),
    );
    assert!(matches!(
        shape_drift,
        Err(PayloadError::ShapeMismatch {
            context: "vector sample",
            ..
        })
    ));

    let component_drift = VectorPayload::new(
        dims_meta(VELOCITY, 0),
        SampleSource::table(
            vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
            vec![
                vec![QtyAny::new(1.0, VELOCITY), QtyAny::new(2.0, VELOCITY)],
                vec![QtyAny::new(3.0, VELOCITY), QtyAny::new(4.0, PRESSURE)],
            ],
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        )
        .expect("coherent table shape"),
    );
    assert!(matches!(
        component_drift,
        Err(PayloadError::DimensionMismatch { index: 3, .. })
    ));

    let unsorted = SpeciesBundle::new(
        dims_meta(AMOUNT, 0),
        SampleSource::fixed(vec![
            SpeciesValue::new(
                SpeciesId::new("H2O").expect("species"),
                QtyAny::new(1.0, AMOUNT),
            ),
            SpeciesValue::new(
                SpeciesId::new("CO2").expect("species"),
                QtyAny::new(1.0, AMOUNT),
            ),
        ]),
    );
    assert!(matches!(
        unsorted,
        Err(PayloadError::NonCanonicalSpeciesAxis { .. })
    ));

    assert!(matches!(
        PayloadMeta::new(
            QuantityContract::Semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Rms,
            )),
            id("basis/temperature"),
            FrameId(0),
            OrientationParity::Even,
            ReferenceSemantics::Continuous,
        ),
        Err(PayloadError::InvalidSemanticContract { .. })
    ));

    let negative_scale = ScalarPayload::new(
        dims_meta(PRESSURE, 0),
        SampleSource::distribution(
            DistributionFamily::Normal,
            vec![QtyAny::new(10.0, PRESSURE), QtyAny::new(-1.0, PRESSURE)],
        )
        .expect("normal has two parameters"),
    );
    assert!(matches!(
        negative_scale,
        Err(PayloadError::DistributionOrder { .. })
    ));
}

#[test]
fn aggregate_budget_counts_outer_and_nested_items_together() {
    let source = SampleSource::table(
        vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
        vec![
            vec![QtyAny::new(1.0, VELOCITY), QtyAny::new(2.0, VELOCITY)],
            vec![QtyAny::new(3.0, VELOCITY), QtyAny::new(4.0, VELOCITY)],
        ],
        TableInterpolation::Linear,
        OutsideDomainPolicy::Refuse,
    )
    .expect("small nested table");
    assert!(matches!(
        VectorPayload::new_with_item_limit(dims_meta(VELOCITY, 0), source, 8),
        Err(PayloadError::TooManyItems {
            actual: 9,
            limit: 8,
            ..
        })
    ));

    let species_source = SampleSource::fixed(vec![
        SpeciesValue::new(
            SpeciesId::new("CO2").expect("species"),
            QtyAny::new(0.4, Dims::NONE),
        ),
        SpeciesValue::new(
            SpeciesId::new("H2O").expect("species"),
            QtyAny::new(0.6, Dims::NONE),
        ),
    ]);
    assert!(matches!(
        SpeciesBundle::new_with_item_limit(dims_meta(Dims::NONE, 0), species_source, 8),
        Err(PayloadError::TooManyItems {
            actual: 9,
            limit: 8,
            ..
        })
    ));
}

#[test]
fn payload_accessors_are_checked_and_preflight_ready() {
    let payloads = all_payload_variants();
    assert_eq!(
        payloads[1]
            .bounded_dynamic_scalar_count()
            .expect("vector table scalar count"),
        8,
        "two table coordinates and six vector components are retained"
    );
    assert_eq!(
        payloads[4]
            .identity_stats()
            .expect("species identities are bounded"),
        (
            "basis/world-cartesian".len() + 2 * ("CO2".len() + "H2O".len()),
            "basis/world-cartesian".len(),
        ),
        "the source species ids and retained canonical-axis cache are both owned"
    );

    for payload in payloads {
        assert_eq!(
            payload.homogeneous_dims().is_none(),
            payload.kind() == PayloadKind::CharacteristicState
        );
        assert!(payload.bounded_dynamic_scalar_count().is_ok());
        assert!(payload.identity_bytes().expect("checked identity sum") > 0);
        assert!(payload.max_identity_component_bytes().expect("checked max") > 0);
        assert!(payload.meta().basis_key().as_str().starts_with("basis/"));
    }
}

#[test]
fn bounded_decoder_refuses_truncation_versions_extensions_and_limits() {
    assert_eq!(
        PayloadDecodeLimits::DEFAULT.max_bytes,
        MAX_PAYLOAD_WIRE_BYTES
    );
    assert!(MAX_PAYLOAD_WIRE_BYTES > MAX_PAYLOAD_ITEMS.saturating_mul(MAX_PAYLOAD_ID_BYTES));
    let payload = all_payload_variants().remove(1);
    let bytes = canonical_payload_bytes(&payload);

    assert!(matches!(
        decode_payload(&bytes[..bytes.len() - 1]),
        Err(PayloadError::Truncated { .. })
    ));

    let mut wrong_version = bytes.clone();
    wrong_version[8..10].copy_from_slice(&2_u16.to_le_bytes());
    assert_eq!(
        decode_payload(&wrong_version),
        Err(PayloadError::UnsupportedVersion { found: 2 })
    );

    let mut extended = bytes.clone();
    extended.push(0);
    assert!(matches!(
        decode_payload(&extended),
        Err(PayloadError::TrailingBytes { .. })
    ));

    let byte_limited = PayloadDecodeLimits {
        max_bytes: bytes.len() - 1,
        ..PayloadDecodeLimits::DEFAULT
    };
    assert!(matches!(
        decode_payload_with_limits(&bytes, byte_limited),
        Err(PayloadError::ByteLimit { .. })
    ));

    let item_limited = PayloadDecodeLimits {
        max_items: 1,
        ..PayloadDecodeLimits::DEFAULT
    };
    assert!(matches!(
        decode_payload_with_limits(&bytes, item_limited),
        Err(PayloadError::ItemLimit { .. })
    ));

    let identifier_limited = PayloadDecodeLimits {
        max_identifier_bytes: 5,
        ..PayloadDecodeLimits::DEFAULT
    };
    assert!(matches!(
        decode_payload_with_limits(&bytes, identifier_limited),
        Err(PayloadError::IdentifierByteLimit { .. })
    ));

    let mut invalid_magic = bytes.clone();
    invalid_magic[0] ^= 0xff;
    assert_eq!(
        decode_payload(&invalid_magic),
        Err(PayloadError::InvalidMagic)
    );

    let mut invalid_variant = bytes.clone();
    invalid_variant[10] = 0xff;
    assert!(matches!(
        decode_payload(&invalid_variant),
        Err(PayloadError::InvalidTag {
            context: "payload variant",
            ..
        })
    ));

    let mut invalid_utf8 = bytes.clone();
    let basis_at = invalid_utf8
        .windows(b"basis/world-cartesian".len())
        .position(|window| window == b"basis/world-cartesian")
        .expect("basis bytes present");
    invalid_utf8[basis_at] = 0xff;
    assert!(matches!(
        decode_payload(&invalid_utf8),
        Err(PayloadError::InvalidUtf8 { .. })
    ));
}

#[test]
fn canonical_float_encoding_collapses_signed_zero_only() {
    let positive = Payload::Scalar(
        ScalarPayload::new(
            dims_meta(Dims::NONE, 0),
            SampleSource::fixed(QtyAny::dimensionless(0.0)),
        )
        .expect("positive zero"),
    );
    let negative = Payload::Scalar(
        ScalarPayload::new(
            dims_meta(Dims::NONE, 0),
            SampleSource::fixed(QtyAny::dimensionless(-0.0)),
        )
        .expect("negative zero"),
    );
    assert_eq!(
        canonical_payload_bytes(&positive),
        canonical_payload_bytes(&negative)
    );

    let mut noncanonical_wire = canonical_payload_bytes(&positive);
    let value_at = noncanonical_wire.len() - 14;
    noncanonical_wire[value_at..value_at + 8].copy_from_slice(&(-0.0_f64).to_bits().to_le_bytes());
    assert!(matches!(
        decode_payload(&noncanonical_wire),
        Err(PayloadError::NonCanonicalFloat { .. })
    ));

    let aliased_support = ScalarPayload::new(
        dims_meta(Dims::NONE, 0),
        SampleSource::distribution(
            DistributionFamily::Empirical,
            vec![QtyAny::dimensionless(-0.0), QtyAny::dimensionless(0.0)],
        )
        .expect("empirical support arity"),
    );
    assert!(matches!(
        aliased_support,
        Err(PayloadError::DistributionOrder { index: 1, .. })
    ));
}
