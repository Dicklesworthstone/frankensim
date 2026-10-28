//! G0/G3: the same canonical material artifacts drive strings and their body.
use super::*;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MODEL_PACK_TARGET_BASIS,
    ModelNormalizationReceipt, ModelNormalizationTarget, NormalizedMaterialCardPack,
    NormalizedModelPack,
};
use fs_qty::{
    QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

fn state(name: &str) -> MaterialStateId {
    MaterialStateId {
        chemistry: format!("manufactured {name}"),
        phase: "solid".into(),
        process: "declared source specimen".into(),
        revision: 0,
    }
}
pub(super) fn canonical(name: &str, density: f64) -> NormalizedMaterialCardPack {
    NormalizedMaterialCardPack::new(state(name), pack(name, density, 200_000.0, None)).unwrap()
}
fn wire_source(
    source: &MaterialStringSource,
) -> &fs_couple::string_specimen::ResolvedStringSpecimen {
    match source {
        MaterialStringSource::Homogeneous(source) => source,
        MaterialStringSource::WithTorsion(torsion) => wire_source(torsion.bending()),
        MaterialStringSource::WithPolarizations(planes) => wire_source(planes.source()),
        _ => panic!("fixture is a homogeneous string"),
    }
}

#[test]
fn canonical_material_cards_preserve_identity_and_legacy_physics_through_pressure() {
    scoped(|cx, gate| {
        let mut legacy = compile(cx, gate, spec(), 7800.0);
        let mut source = spec()
            .compile(cx, gate, |name| Ok(canonical(name, 7800.0)))
            .unwrap();
        for (name, entry) in &source.catalog {
            let pack = canonical(name, 7800.0);
            assert_eq!(entry.card, *pack.card());
            assert_eq!(entry.pack.content_hash(), pack.content_hash());
            assert_eq!(entry.pack.to_bytes(), pack.to_bytes());
            assert_ne!(entry.card.id(), legacy.catalog[name].card.id());
        }
        let bank = source.performance.bank().source().bank();
        assert_eq!(
            wire_source(&bank.strings()[0].source).material().material(),
            &state("wire.fsmatpk")
        );
        assert_eq!(
            bank.material_plate()
                .unwrap()
                .specimen
                .material()
                .material(),
            &state("body.fsmatpk")
        );
        assert_eq!(
            source.performance.bank().system().components()[0].modes(),
            legacy.performance.bank().system().components()[0].modes()
        );
        assert_eq!(pressure(&mut source, 37), pressure(&mut legacy, 128));
    });
}

#[test]
fn canonical_material_data_swaps_reach_string_body_and_joint_sound() {
    scoped(|cx, gate| {
        let make = |changed: &str| {
            spec()
                .compile(cx, gate, |name| {
                    Ok(canonical(
                        name,
                        if name == changed { 15600.0 } else { 7800.0 },
                    ))
                })
                .unwrap()
        };
        let mut original = make("");
        let mut string = make("wire.fsmatpk");
        let mut body = make("body.fsmatpk");
        assert_eq!(original.spec.bytes, string.spec.bytes);
        assert_eq!(original.spec.bytes, body.spec.bytes);
        assert!(
            string.performance.bank().cards()[0].fundamental_hz()
                < original.performance.bank().cards()[0].fundamental_hz()
        );
        assert_eq!(
            string.performance.bank().system().components()[1].modes(),
            original.performance.bank().system().components()[1].modes()
        );
        assert_ne!(
            body.performance.bank().system().components()[1].modes(),
            original.performance.bank().system().components()[1].modes()
        );
        assert_eq!(
            body.performance.bank().cards()[0].fundamental_hz(),
            original.performance.bank().cards()[0].fundamental_hz()
        );
        let reference = pressure(&mut original, 37);
        for ready in [&mut string, &mut body] {
            assert!(
                pressure(ready, 128)
                    .iter()
                    .zip(&reference)
                    .any(|(a, b)| (a - b).abs() > 1e-12)
            );
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

fn model_source() -> NormalizedModelPack {
    let artifact = fs_blake3::ContentHash([17; 32]);
    let model = ConstitutiveModelCard {
        law: LawId("manufactured-associated-elasticity".into()),
        law_version: 1,
        parameters: BTreeMap::from([(
            "young_modulus".into(),
            LawParameter {
                value: 200e9,
                dims: Pressure::DIMS,
            },
        )]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained(),
        sources: vec![artifact],
        provenance: Provenance {
            source: "synthetic source association, not a new solver".into(),
            license: "CC0-1.0".into(),
            artifact: Some(artifact),
        },
    };
    let receipt = ModelNormalizationReceipt::new(
        ModelNormalizationTarget::Parameter {
            model: model.content_hash(),
            parameter: "young_modulus".into(),
        },
        artifact,
        Pressure::DIMS,
        1.0,
        0.0,
        "authored SI",
        MODEL_PACK_TARGET_BASIS,
        None,
        None,
    );
    NormalizedModelPack::new(
        "associated-model",
        "synthetic-v1",
        artifact,
        "CC0-1.0",
        vec![model],
        vec![receipt],
    )
    .unwrap()
}

#[test]
fn canonical_material_decoder_retains_model_association_and_refuses_damaged_cards() {
    let claims = pack("wire.fsmatpk", 7800.0, 200_000.0, None);
    let models = model_source();
    let source =
        NormalizedMaterialCardPack::new_with_models(state("wire.fsmatpk"), claims, models.clone())
            .unwrap();
    let bytes = source.to_bytes();
    let decoded = MaterialPack::from_bytes(&bytes).unwrap();
    assert_eq!(decoded.to_bytes(), bytes);
    assert_eq!(decoded.content_hash(), source.content_hash());
    assert_eq!(decoded.card().unwrap().models(), models.models());
    assert!(decoded.to_bytes().len() > source.claims_pack().to_bytes().len());
    scoped(|cx, gate| {
        let ready = spec()
            .compile(cx, gate, |name| {
                Ok(if name == "wire.fsmatpk" {
                    MaterialPack::from_bytes(&bytes).unwrap()
                } else {
                    canonical(name, 7800.0).into()
                })
            })
            .unwrap();
        assert_eq!(ready.catalog["wire.fsmatpk"].card.models(), models.models());
        assert_eq!(
            wire_source(&ready.performance.bank().source().bank().strings()[0].source)
                .material()
                .card_identity(),
            source.card().content_hash()
        );
    });
    for length in [8, 12, bytes.len() - 1] {
        assert!(MaterialPack::from_bytes(&bytes[..length]).is_err());
    }
    let mut corrupt = bytes.clone();
    corrupt[16] ^= 1;
    assert!(MaterialPack::from_bytes(&corrupt).is_err());
}

#[test]
fn canonical_material_byte_budget_counts_associated_models() {
    let base = model_source();
    let mut models = Vec::new();
    let mut receipts = Vec::new();
    for i in 0..9 {
        let mut model = base.models()[0].clone();
        model.law = LawId(format!("manufactured-associated-elasticity-{i}"));
        model.provenance.source = "x".repeat(1_000_000);
        receipts.push(ModelNormalizationReceipt::new(
            ModelNormalizationTarget::Parameter {
                model: model.content_hash(),
                parameter: "young_modulus".into(),
            },
            model.sources[0],
            Pressure::DIMS,
            1.0,
            0.0,
            "authored SI",
            MODEL_PACK_TARGET_BASIS,
            None,
            None,
        ));
        models.push(model);
    }
    let models = NormalizedModelPack::new(
        "large-model-source",
        "synthetic-v1",
        fs_blake3::ContentHash([17; 32]),
        "CC0-1.0",
        models,
        receipts,
    )
    .unwrap();
    let source = NormalizedMaterialCardPack::new_with_models(
        state("wire.fsmatpk"),
        pack("wire.fsmatpk", 7800.0, 200_000.0, None),
        models,
    )
    .unwrap();
    assert!(source.claims_pack().to_bytes().len() < input::MAX_PACK_BYTES);
    assert!(source.to_bytes().len() > input::MAX_PACK_BYTES);
    scoped(|cx, gate| {
        let error = spec()
            .compile(cx, gate, |_| Ok(source.clone()))
            .err()
            .unwrap();
        assert!(
            error
                .to_string()
                .contains("combined material and interface pack byte budget"),
            "{error}"
        );
    });
}

fn typed(
    name: &str,
    temperature_kind: QuantityKind,
    frequency_kind: QuantityKind,
) -> NormalizedMaterialCardPack {
    let base = pack(name, 7800.0, 200_000.0, None);
    let mut claims = ClaimSet::new();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.observations.clear();
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", quantity(temperature_kind), 280.0, 320.0)
            .with_quantity("omega", quantity(frequency_kind), 0.0, 200_000.0);
        claims.insert_claim(claim).unwrap();
    }
    let claims = NormalizedPack::new(
        name,
        "synthetic-typed-v1",
        base.source_artifact(),
        "CC0-1.0",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedMaterialCardPack::new(state(name), claims).unwrap()
}

#[test]
fn canonical_material_axes_preserve_kelvin_and_angular_frequency_without_relabelling() {
    scoped(|cx, gate| {
        let temperature = QuantityKind::AbsoluteTemperature;
        let frequency = QuantityKind::Frequency(FrequencyConvention::Angular);
        let mut reference = compile(cx, gate, spec(), 7800.0);
        let mut source = spec()
            .compile(cx, gate, |name| Ok(typed(name, temperature, frequency)))
            .unwrap();
        let material =
            wire_source(&source.performance.bank().source().bank().strings()[0].source).material();
        assert_eq!(
            material.axis_quantities()["T"],
            QuantitySpec::semantic(SemanticType::new(temperature, ValueForm::Static))
        );
        assert_eq!(
            material.axis_quantities()["omega"],
            QuantitySpec::semantic(SemanticType::new(frequency, ValueForm::Static))
        );
        assert_eq!(pressure(&mut source, 37), pressure(&mut reference, 128));
        for (t, omega) in [
            (QuantityKind::TemperatureDifference, frequency),
            (
                temperature,
                QuantityKind::Frequency(FrequencyConvention::Cyclic),
            ),
        ] {
            let error = spec()
                .compile(cx, gate, |name| Ok(typed(name, t, omega)))
                .err()
                .unwrap();
            assert!(
                error.to_string().contains("kelvin/radian-per-second"),
                "{error}"
            );
        }
    });
}
