//! Real normalized-source decoding through the native material stop consumer.
use super::*;
use fs_contact::interface_binding::{ADHESION_ENERGY_DIMS, ADHESION_ENERGY_PROPERTY};
use fs_matdb::{MODEL_PACK_TARGET_BASIS, ModelNormalizationReceipt, ModelNormalizationTarget,
    NormalizedModelPack, SurfaceSpec, SystemContext, ValidityBoundSide};
#[path = "../../../../fs-contact/tests/source_point/common.rs"]
mod source_fixture;
#[path = "normal_sources/shell.rs"]
mod shell;

fn provenance() -> Provenance {
    Provenance { source: "synthetic native normal-contact test".into(),
        license: "CC0-1.0".into(), artifact: Some(fs_blake3::ContentHash([29; 32])) }
}
fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim { key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, 200_000.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: provenance() }
}
fn pack_claims(name: &str, claims: ClaimSet) -> NormalizedPack {
    let pack = NormalizedPack::new(name, "synthetic-native-contact-v1",
        fs_blake3::ContentHash([29; 32]), "CC0-1.0 synthetic data", claims, vec![], vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}
fn bulk(name: &str, receiver_young: f64) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [("density", 7800.0, Density::DIMS),
        ("young_modulus", if name == "receiver.fsmatpk" { receiver_young } else { 200e9 }, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE), ("yield_stress", 1e10, Pressure::DIMS),
        (ETA, 10.0, DynViscosity::DIMS), (SHEAR, 80e9, Pressure::DIMS),
        (SHEAR_ETA, 1000.0, DynViscosity::DIMS)] {
        claims.insert_claim(claim(key, value, dims)).unwrap();
    }
    pack_claims(name, claims)
}
fn interface(outer: &str, receiver: &str) -> NormalizedInterfacePack {
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim(ADHESION_ENERGY_PROPERTY, 0.0, ADHESION_ENERGY_DIMS)).unwrap();
    let surface = |path: &str| SurfaceSpec {
        material: MaterialStateId { chemistry: format!("source-pack:{path}"),
            phase: "declared solid model".into(), process: "frozen source-pack realization".into(), revision: 0 },
        texture_frame: "explicit local point patch".into(),
    };
    NormalizedInterfacePack::new(surface(outer), surface(receiver), SystemContext {
        medium: "dry".into(), third_body: None, environment: "synthetic fixed state".into(),
        history: "initial".into(),
    }, pack_claims("normal-interface", claims)).unwrap()
}
fn models() -> NormalizedModelPack {
    let a = source_fixture::solid("A", 7800.0, 200e9);
    let b = source_fixture::solid("B", 7800.0, 2e7);
    let mut card = source_fixture::bind(&a, &b, Some(0.2), 0.0).model_card().clone();
    let literal = fs_blake3::ContentHash([31; 32]);
    card.sources = vec![literal];
    card.provenance.artifact = Some(literal);
    let model = card.content_hash();
    let mut receipts = Vec::new();
    for (name, value) in &card.parameters {
        receipts.push(ModelNormalizationReceipt::new(ModelNormalizationTarget::Parameter {
            model, parameter: name.clone() }, literal, value.dims, 1.0, 0.0,
            "synthetic SI literal", MODEL_PACK_TARGET_BASIS, None, None));
    }
    for axis in card.validity.bounds().keys() {
        for side in [ValidityBoundSide::Lower, ValidityBoundSide::Upper] {
            receipts.push(ModelNormalizationReceipt::new(ModelNormalizationTarget::ValidityBound {
                model, axis: axis.clone(), side }, literal, Dims([0,0,0,1,0,0]),
                1.0, 0.0, "K", MODEL_PACK_TARGET_BASIS, None, None));
        }
    }
    NormalizedModelPack::new("synthetic-point-law", "synthetic-native-v1", literal,
        "CC0-1.0 synthetic data", vec![card], receipts).unwrap()
}
const PROFILE: &str = "frankensim-string-stops-v1
material-stop finger 0 0.25 positive fixed -1e-9 receiver.fsmatpk normal.fsintpk normal.fsmodelpk 5000 5000 0.1 0.00005 0.00005 0.01 0.02 0.03 10 1e-8
";
fn profile(gate: &CancelGate, outer: &str, receiver: &str) -> stops::Profile {
    let mut profile = stops::Profile::parse(PROFILE.as_bytes()).unwrap();
    profile.load_sources(gate, |path| Ok(match path {
        "normal.fsintpk" => interface(outer, receiver).to_bytes(),
        "normal.fsmodelpk" => models().to_bytes(),
        _ => panic!("unexpected source {path}"),
    })).unwrap();
    profile
}
fn resting_spec() -> input::Spec {
    let mut s = spec();
    s.forces = None;
    s
}

#[test]
fn receiver_material_swap_changes_reaction_and_pressure_on_the_same_physical_bank() {
    scoped(|cx, gate| {
        let build = |young| resting_spec().compile_with_profiles(cx, gate, None,
            Some(profile(gate, "wire.fsmatpk", "receiver.fsmatpk")),
            |name| Ok(bulk(name, young))).unwrap();
        let mut a = build(2e7);
        let mut b = build(4e7);
        let sa = a.performance.bank().system().obstacle_source(0).unwrap();
        let sb = b.performance.bank().system().obstacle_source(0).unwrap();
        assert_ne!(sa.model().identity(), sb.model().identity());
        assert!(sb.stiffness_n_per_m_3_2() > sa.stiffness_n_per_m_3_2());
        assert_eq!(a.spec.bytes, b.spec.bytes);
        for path in ["wire.fsmatpk", "body.fsmatpk"] {
            assert_eq!(a.catalog[path].pack.content_hash(), b.catalog[path].pack.content_hash());
        }
        for (x, y) in a.performance.bank().system().components().iter()
            .zip(b.performance.bank().system().components()) { assert_eq!(x.modes(), y.modes()); }
        let x = pressure(&mut a, 37);
        let y = pressure(&mut b, 128);
        assert_ne!(x, y);
        assert!(x.iter().any(|p| *p != 0.0));
        assert_ne!(a.performance.bank().system().last_obstacle_frames(),
            b.performance.bank().system().last_obstacle_frames());
        let mut partitioned = build(2e7);
        assert_eq!(x, pressure(&mut partitioned, 1));
        assert!(a.performance.bank().system().obstacle_source(0).is_some());
    });
}

#[test]
fn canonical_material_stop_sources_preserve_ordered_identity_and_physics() {
    scoped(|cx, gate| {
        let card = |name: &str, young| fs_matdb::NormalizedMaterialCardPack::new(
            MaterialStateId { chemistry: format!("manufactured {name}"), phase: "solid".into(),
                process: "canonical point-contact specimen".into(), revision: 0 },
            bulk(name, young)).unwrap();
        let original = interface("wire.fsmatpk", "receiver.fsmatpk");
        let mut a = original.card().surface_a().clone();
        a.material = card("wire.fsmatpk", 2e7).card().id().clone();
        let mut b = original.card().surface_b().clone();
        b.material = card("receiver.fsmatpk", 2e7).card().id().clone();
        let ordered = NormalizedInterfacePack::new(a, b, original.card().context().clone(),
            original.claims_pack().clone()).unwrap();
        let mut sourced = profile(gate, "wire.fsmatpk", "receiver.fsmatpk");
        sourced.load_sources(gate, |path| Ok(match path {
            "normal.fsintpk" => ordered.to_bytes(),
            "normal.fsmodelpk" => models().to_bytes(),
            _ => panic!("unexpected source {path}"),
        })).unwrap();
        let mut legacy = resting_spec().compile_with_profiles(cx, gate, None,
            Some(profile(gate, "wire.fsmatpk", "receiver.fsmatpk")),
            |name| Ok(bulk(name, 2e7))).unwrap();
        let mut canonical = resting_spec().compile_with_profiles(cx, gate, None,
            Some(sourced.clone()), |name| Ok(card(name, 2e7))).unwrap();
        let mut replacement = resting_spec().compile_with_profiles(cx, gate, None,
            Some(sourced.clone()), |name| Ok(card(name, 4e7))).unwrap();
        for name in ["wire.fsmatpk", "receiver.fsmatpk"] {
            assert_eq!(canonical.catalog[name].card, *card(name, 2e7).card());
        }
        let before = canonical.performance.bank().system().obstacle_source(0).unwrap();
        let after = replacement.performance.bank().system().obstacle_source(0).unwrap();
        assert!(after.stiffness_n_per_m_3_2() > before.stiffness_n_per_m_3_2());
        let expected = pressure(&mut canonical, 37);
        assert_eq!(expected, pressure(&mut legacy, 128));
        assert_ne!(expected, pressure(&mut replacement, 128));

        // The original contact owner requires matching bulk/interface schemas.
        // A canonical envelope must not erase or retag its receiver's axes.
        use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm}};
        let mut claims = ClaimSet::new();
        for (_, original) in bulk("receiver.fsmatpk", 2e7).claims().claims_ordered() {
            let mut claim = original.clone();
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("T", QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature, ValueForm::Static)), 280.0, 320.0)
                .with_quantity("omega", QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static)), 0.0, 200_000.0);
            claims.insert_claim(claim).unwrap();
        }
        let receiver = fs_matdb::NormalizedMaterialCardPack::new(
            card("receiver.fsmatpk", 2e7).card().id().clone(),
            pack_claims("typed-receiver", claims)).unwrap();
        let error = resting_spec().compile_with_profiles(cx, gate, None,
            Some(sourced), |name| Ok(if name == "receiver.fsmatpk" {
                receiver.clone()
            } else { card(name, 2e7) })).err().unwrap();
        assert!(error.to_string().contains("AxisQuantityMismatch"), "{error}");
    });
}

#[test]
fn ordered_outer_material_and_source_domains_are_mandatory_before_runtime_binding() {
    scoped(|cx, gate| {
        let missing = stops::Profile::parse(PROFILE.as_bytes()).unwrap();
        let mut bulk_reads = 0;
        assert!(resting_spec().compile_with_profiles(cx, gate, None, Some(missing), |name| {
            bulk_reads += 1; Ok(bulk(name, 2e7))
        }).is_err());
        assert_eq!(bulk_reads, 0);
        for (outer, receiver) in [("other", "receiver.fsmatpk"), ("wire.fsmatpk", "other")] {
            assert!(resting_spec().compile_with_profiles(cx, gate, None, Some(profile(gate, outer, receiver)),
                |name| Ok(bulk(name, 2e7))).is_err());
        }
        for (surface, accepted) in [("outer.fsmatpk", true), ("wire.fsmatpk", false)] {
            let mut s = resting_spec();
            s.strings[0].section = input::WireSection::Bonded { core_radius_m: 0.0002,
                sleeves: vec![input::WireLayer { pack: "outer.fsmatpk".into(), outer_radius_m: 0.0003 }] };
            assert_eq!(s.compile_with_profiles(cx, gate, None,
                Some(profile(gate, surface, "receiver.fsmatpk")), |name| Ok(bulk(name, 2e7))).is_ok(), accepted);
        }
        let mut hot = resting_spec();
        hot.temperature = 400.0;
        assert!(hot.compile_with_profiles(cx, gate, None, Some(profile(gate, "wire.fsmatpk", "receiver.fsmatpk")),
            |name| Ok(bulk(name, 2e7))).is_err());
        let mut thin = resting_spec();
        thin.strings[0].radius = 1e-5;
        assert!(thin.compile_with_profiles(cx, gate, None, Some(profile(gate, "wire.fsmatpk", "receiver.fsmatpk")),
            |name| Ok(bulk(name, 2e7))).is_err());
    });
}

#[test]
fn source_record_parser_and_decoder_refuse_unsupported_inputs_without_fallback() {
    let gate = CancelGate::new_clock_free();
    for text in [PROFILE.replace("positive fixed", "positive bridge:0"),
        PROFILE.replace("receiver.fsmatpk", "../receiver.fsmatpk"),
        PROFILE.replace("5000 5000", "5000 0"),
        PROFILE.replace("0.1 0.00005", "-0.1 0.00005"),
        PROFILE.replace("0.01 0.02 0.03", "-0.01 0.02 0.03")] {
        assert!(stops::Profile::parse(text.as_bytes()).is_err());
    }
    let mut p = profile(&gate, "wire.fsmatpk", "receiver.fsmatpk");
    let bytes = p.source_pack_bytes().unwrap();
    assert!(p.load_sources(&gate, |_| Ok(vec![0; 16])).is_err());
    assert_eq!(p.source_pack_bytes().unwrap(), bytes, "failed loading replaced retained source packs");
}
