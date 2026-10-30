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


fn held_pack(name: &str, typed: bool, conflict: bool) -> NormalizedMaterialCardPack {
    use fs_qty::semantic::{CompositionBasis, MoistureBasis};
    let base = pack(name, 7800.0, 200_000.0, None);
    let mut claims = ClaimSet::new();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        claim.observations.clear();
        let (value, dims) = match &claim.value {
            PropertyValue::Scalar { value, dims } => (*value, *dims),
            _ => panic!("manufactured base has scalar coefficients"),
        };
        let factor = match claim.key.name() {
            "density" => 1.2,
            "young_modulus" | "shear_modulus" => 0.6,
            ETA | SHEAR_ETA => 2.0,
            _ => 1.0,
        };
        let kind = if conflict && claim.key.name() == ETA {
            QuantityKind::Composition(CompositionBasis::MoleFraction)
        } else { QuantityKind::Moisture(MoistureBasis::DryMass) };
        claim.validity = ValidityDomain::unconstrained()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
            .with("omega", 0.0, 200_000.0);
        claim.validity = if typed {
            claim.validity.with_quantity("moisture_content", quantity(kind), 0.0, 0.2)
        } else { claim.validity.with("moisture_content", 0.0, 0.2) };
        // A source-declared T/moisture surface, not a humidity-to-moisture fit.
        // At 300 K the two moisture edges give value and factor*value.
        claim.value = PropertyValue::Surface2D {
            axes: ["T".into(), "moisture_content".into()],
            axis_dims: [fs_qty::Temperature::DIMS, Dims::NONE],
            samples: vec![([280.0, 0.0], 0.98*value), ([320.0, 0.0], 1.02*value),
                ([280.0, 0.2], (factor-0.02)*value), ([320.0, 0.2], (factor+0.02)*value)],
            triangles: vec![[0, 1, 2], [1, 3, 2]], dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(claim).unwrap();
    }
    let claims = NormalizedPack::new(name, "manufactured-held-environment-v1",
        base.source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(state(name), claims).unwrap()
}

fn held_text(wire: f64, body: f64) -> String {
    include_str!("../pluck.material-strings")
        .replace("48000 4096", "48000 64")
        .replace("state 293.15 1000", &format!(
            "state 300 1000\nmaterial-states 2\nmaterial-state wire.fsmatpk 1 moisture_content {wire}\nmaterial-state body.fsmatpk 1 moisture_content {body}"))
        .replace("initial\tbow\t0e0\t0e0", "initial\tbow\t1e-2\t1e-3")
}

#[test]
fn held_source_state_changes_string_material_and_body_in_the_same_bowed_pressure_path() {
    scoped(|cx, gate| {
        let make = |wire, body| input::Spec::parse(held_text(wire, body).as_bytes()).unwrap()
            .compile(cx, gate, |name| Ok(held_pack(name, true, false))).unwrap();
        let mut reference = make(0.0, 0.0);
        let mut string = make(0.2, 0.0);
        let mut body = make(0.0, 0.2);
        let extension = |moisture| {
            let text = held_text(moisture, 0.0)
                .replace("frankensim-material-strings-v1", "frankensim-material-strings-v2")
                .replace("string wire.fsmatpk 0.5 0.0003 20 4 3 0.8 3000 0.1",
                    "string wire.fsmatpk 0.5 0.0003 extension 0.4998 0.01 4 3 0.8 3000 0.1\nsection homogeneous");
            input::Spec::parse(text.as_bytes()).unwrap()
                .compile(cx, gate, |name| Ok(held_pack(name, true, false))).unwrap()
        };
        let mut extension_dry = extension(0.0);
        let mut extension_changed = extension(0.2);
        for name in ["wire.fsmatpk", "body.fsmatpk"] {
            assert_eq!(reference.catalog[name].pack.content_hash(), string.catalog[name].pack.content_hash());
            assert_eq!(reference.catalog[name].pack.content_hash(), body.catalog[name].pack.content_hash());
        }
        let original = wire_source(&reference.performance.bank().source().bank().strings()[0].source);
        let changed = wire_source(&string.performance.bank().source().bank().strings()[0].source);
        assert_eq!(original.material().card_identity(), changed.material().card_identity());
        assert_ne!(original.material().identity(), changed.material().identity());
        let ratio = |a: f64, b: f64, expected: f64| assert!((a/b-expected).abs() < 1e-12);
        ratio(changed.string().lin_density_kg_m, original.string().lin_density_kg_m, 1.2);
        ratio(changed.string().bending_stiffness_n_m2, original.string().bending_stiffness_n_m2, 0.6);
        ratio(changed.string().axial_stiffness_n, original.string().axial_stiffness_n, 0.6);
        ratio(changed.string().kelvin_voigt_bending.unwrap().viscous_stiffness_n_m2_s,
            original.string().kelvin_voigt_bending.unwrap().viscous_stiffness_n_m2_s, 2.0);
        assert_eq!(changed.string().tension_n, original.string().tension_n, "fixed force is not fixed pitch");
        let supported = wire_source(&extension_dry.performance.bank().source().bank().strings()[0].source);
        let supported_changed = wire_source(&extension_changed.performance.bank().source().bank().strings()[0].source);
        ratio(supported_changed.string().tension_n, supported.string().tension_n, 0.6);
        assert_eq!(supported_changed.string().length_m, supported.string().length_m);
        assert!(extension_changed.performance.bank().cards()[0].fundamental_hz()
            < extension_dry.performance.bank().cards()[0].fundamental_hz());
        assert_eq!(changed.material().property("young_modulus").unwrap().source_axes().len(), 2);
        assert_eq!(changed.material().query_point().iter().find(|(axis, _)| axis == "moisture_content").unwrap().1, 0.2);
        assert_eq!(string.performance.bank().system().components()[1].modes(),
            reference.performance.bank().system().components()[1].modes());
        assert_ne!(body.performance.bank().system().components()[1].modes(),
            reference.performance.bank().system().components()[1].modes());
        assert_eq!(body.performance.bank().cards()[0].fundamental_hz(),
            reference.performance.bank().cards()[0].fundamental_hz());
        let expected = pressure(&mut reference, 16);
        assert!(expected.iter().any(|value| value.abs() > 1e-12));
        let extension_pressure = pressure(&mut extension_dry, 16);
        assert!(pressure(&mut extension_changed, 32).iter().zip(extension_pressure)
            .any(|(a,b)| (a-b).abs() > 1e-12));
        for run in [&mut string, &mut body] {
            assert!(pressure(run, 32).iter().zip(&expected).any(|(a,b)| (a-b).abs() > 1e-12));
            let frame = run.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn held_material_axes_refuse_ambiguous_semantics_unowned_conditions_and_unused_packs() {
    let card = held_pack("wire.fsmatpk", true, false);
    for axes in [
        vec![("T".into(), 300.0)],
        vec![("omega".into(), 1000.0)],
        vec![("moisture_content".into(), 0.1), ("moisture_content".into(), 0.2)],
        vec![("undeclared".into(), 0.1)],
    ] {
        assert!(material_inputs::query_with_environment(card.card(), 300.0, 1000.0, &axes).is_err());
    }
    let axes = vec![("moisture_content".into(), 0.1)];
    for source in [held_pack("wire.fsmatpk", false, false), held_pack("wire.fsmatpk", true, true)] {
        assert!(material_inputs::query_with_environment(source.card(), 300.0, 1000.0, &axes).is_err());
    }
    // Re-querying temperature for contact heat must keep the same held context.
    let a = material_inputs::query_with_environment(card.card(), 290.0, 1000.0, &axes).unwrap();
    let b = material_inputs::query_with_environment(card.card(), 310.0, 1000.0, &axes).unwrap();
    assert_eq!(a.axes()["moisture_content"], b.axes()["moisture_content"]);
    assert_eq!(a.axis_quantities(), b.axis_quantities());
    scoped(|cx, gate| {
        assert!(input::Spec::parse(held_text(0.3, 0.0).as_bytes()).unwrap()
            .compile(cx, gate, |name| Ok(held_pack(name, true, false))).is_err());
        let text = include_str!("../pluck.material-strings").replace("state 293.15 1000",
            "state 293.15 1000\nmaterial-states 1\nmaterial-state unused.fsmatpk 1 moisture_content 0.1");
        let error = input::Spec::parse(text.as_bytes()).unwrap()
            .compile(cx, gate, |name| Ok(canonical(name, 7800.0))).err().unwrap();
        assert!(error.to_string().contains("unused"), "{error}");
    });
    let measured = include_str!("../measured-cello-g.material-strings");
    let parsed = input::Spec::parse(measured.as_bytes()).unwrap();
    let text = measured.replacen("\ngas ", &format!(
        "\nmaterial-states 1\nmaterial-state {} 1 moisture_content 0.1\ngas ",
        parsed.strings[0].pack), 1);
    assert!(input::Spec::parse(text.as_bytes()).is_err(),
        "a source's unknown measurement temperature cannot be overwritten by a bulk-state binding");
}
