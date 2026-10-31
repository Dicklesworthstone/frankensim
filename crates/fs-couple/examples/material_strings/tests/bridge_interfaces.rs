//! Source-defined joints reach the original coupled mechanics and pressure.
use super::*;
use fs_matdb::{SurfaceSpec, SystemContext};
use fs_material::state_point::{interface_stiffness_per_area_dims, interface_viscosity_per_area_dims};

pub(crate) fn interface(k: f64, c: f64, outer: &str, body: &str, band: f64) -> NormalizedInterfacePack {
    let artifact = fs_blake3::ContentHash([17; 32]);
    let provenance = Provenance { source: "synthetic bridge interface".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [("normal_stiffness", k, interface_stiffness_per_area_dims()),
        ("normal_viscosity", c, interface_viscosity_per_area_dims())] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, band),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: provenance.clone() }).unwrap();
    }
    let claims = NormalizedPack::new("synthetic-bridge-interface", "synthetic-v1", artifact,
        "CC0-1.0 synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    let surface = |pack: &str| SurfaceSpec {
        material: MaterialStateId { chemistry: format!("source-pack:{pack}"),
            phase: "declared solid model".into(), process: "frozen source-pack realization".into(), revision: 0 },
        texture_frame: "declared scalar joint direction".into(),
    };
    let pack = NormalizedInterfacePack::new(surface(outer), surface(body), SystemContext {
        medium: "bonded".into(), third_body: Some("synthetic joint layer".into()),
        environment: "declared frozen state".into(), history: "synthetic fixed history".into(),
    }, claims).unwrap();
    NormalizedInterfacePack::from_bytes(&pack.to_bytes()).unwrap()
}

fn sourced_spec() -> input::Spec {
    let text = String::from_utf8(spec().bytes).unwrap().replace("0.8 3000 0.1",
        "0.8 interface bridge.fsintpk 0.0001 normal_stiffness normal_viscosity");
    input::Spec::parse(text.as_bytes()).unwrap()
}
fn load(s: &mut input::Spec, gate: &CancelGate, source: NormalizedInterfacePack) {
    s.bridge_interfaces.load(s.temperature, s.omega, gate, |_| Ok(source.clone())).unwrap();
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, k: f64, c: f64) -> Ready {
    let mut s = sourced_spec();
    load(&mut s, gate, interface(k, c, "wire.fsmatpk", "body.fsmatpk", 200_000.0));
    compile(cx, gate, s, 7800.0)
}

#[test]
fn source_law_keeps_existing_physics_and_its_receipts_through_preload_and_rendering() {
    scoped(|cx, gate| {
        let mut sourced = ready(cx, gate, 30e6, 1000.0);
        let mut original = compile(cx, gate, spec(), 7800.0);
        assert_eq!(sourced.performance.preload_energy_j(), original.performance.preload_energy_j());
        let source = sourced.performance.bank().source().bank();
        let law = &source.material_bridge_interfaces()[0].interface;
        assert_eq!(law.area_m2(), 0.0001);
        assert_eq!(law.stiffness_n_m(), source.strings()[0].connections[0].stiffness_n_m);
        assert_eq!(law.damping_n_s_m(), source.strings()[0].connections[0].damping_n_s_m);
        assert_eq!(law.resolved().surface_materials(), &source.bridge_surface_materials(0).unwrap());
        let expected = pressure(&mut original, 128);
        assert_eq!(pressure(&mut sourced, 37), expected);
        assert_eq!(sourced.performance.bank().source().bank().material_bridge_interfaces().len(), 1);
        let frame = sourced.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn changing_only_interface_source_changes_reciprocal_pressure_without_rewriting_bulk_or_geometry() {
    scoped(|cx, gate| {
        let mut original = ready(cx, gate, 30e6, 1000.0);
        let expected = pressure(&mut original, 37);
        for (k, c) in [(60e6, 1000.0), (30e6, 4000.0)] {
            let mut changed = ready(cx, gate, k, c);
            assert_eq!(original.spec.bytes, changed.spec.bytes);
            for path in ["wire.fsmatpk", "body.fsmatpk"] {
                assert_eq!(original.catalog[path].pack.content_hash(), changed.catalog[path].pack.content_hash());
            }
            assert_ne!(expected, pressure(&mut changed, 128));
            let frame = changed.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let mut partitioned = ready(cx, gate, 30e6, 1000.0);
        assert_eq!(expected, pressure(&mut partitioned, 1));
    });
}

#[test]
fn missing_stale_mismatched_and_out_of_band_sources_refuse_before_output() {
    scoped(|cx, gate| {
        let mut calls = 0;
        assert!(sourced_spec().compile(cx, gate, |name| {
            calls += 1; Ok(pack(name, 7800.0, 200_000.0, None))
        }).is_err());
        assert_eq!(calls, 0);
        let mut stale = sourced_spec();
        load(&mut stale, gate, interface(30e6, 1000.0, "wire.fsmatpk", "body.fsmatpk", 200_000.0));
        stale.temperature = 300.0;
        assert!(stale.compile(cx, gate, |name| {
            calls += 1; Ok(pack(name, 7800.0, 200_000.0, None))
        }).is_err());
        assert_eq!(calls, 0);
        for (outer, body, band) in [("other", "body.fsmatpk", 200_000.0),
            ("wire.fsmatpk", "other", 200_000.0), ("wire.fsmatpk", "body.fsmatpk", 1000.0)] {
            let mut s = sourced_spec();
            load(&mut s, gate, interface(30e6, 1000.0, outer, body, band));
            assert!(s.compile(cx, gate, |name| Ok(pack(name, 7800.0, 200_000.0, None))).is_err());
        }
    });
}

#[test]
fn layered_two_plane_strings_bind_the_physical_outer_sleeve_not_the_core() {
    scoped(|cx, gate| {
        let text = String::from_utf8(sourced_spec().bytes).unwrap()
            .replace("frankensim-material-strings-v1", "frankensim-material-strings-v3")
            .replace("0.0003 20 4", "0.0003 tension 20 4")
            .replace("normal_viscosity\ncavity", "normal_viscosity\nsection bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003\npolarizations two 1 0 1 0\ncavity");
        for (surface, accepted) in [("sleeve.fsmatpk", true), ("wire.fsmatpk", false)] {
            let mut s = input::Spec::parse(text.as_bytes()).unwrap();
            load(&mut s, gate, interface(30e6, 1000.0, surface, "body.fsmatpk", 200_000.0));
            let result = s.compile(cx, gate, |name| Ok(pack(name, 7800.0, 200_000.0, None)));
            assert_eq!(result.is_ok(), accepted);
        }
    });
}


mod held_environment {
    use super::*;
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, MoistureBasis, QuantityKind, SemanticType, ValueForm}};

    fn source() -> NormalizedInterfacePack {
        let base = interface(30e6, 1000.0, "wire.fsmatpk", "body.fsmatpk", 200_000.0);
        let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let mut claims = ClaimSet::new();
        for id in base.claims_pack().claims().observation_ids() {
            claims.register_observation(base.claims_pack().claims().observation(id).unwrap().clone()).unwrap();
        }
        for (_, original) in base.claims_pack().claims().claims_ordered() {
            let mut claim = original.clone();
            let PropertyValue::Scalar { value, dims } = claim.value else { unreachable!() };
            let factor = if claim.key.name() == "normal_stiffness" { 2.0 } else { 4.0 };
            claim.value = PropertyValue::Curve { abscissa: "moisture_content".into(),
                abscissa_dims: Dims::NONE, knots: vec![(0.0, value), (0.2, factor*value)], dims };
            claim.validity = ValidityDomain::unconstrained()
                .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 280.0, 320.0)
                .with_quantity("omega", quantity(QuantityKind::Frequency(FrequencyConvention::Angular)), 0.0, 200_000.0)
                .with_quantity("moisture_content", quantity(QuantityKind::Moisture(MoistureBasis::DryMass)), 0.0, 0.2);
            claim.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(claim).unwrap();
        }
        let claims = NormalizedPack::new("synthetic-joint-environment", "synthetic-v1",
            base.claims_pack().source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
        NormalizedInterfacePack::new(base.card().surface_a().clone(), base.card().surface_b().clone(),
            base.card().context().clone(), claims).unwrap()
    }
    fn environment(moisture: f64) -> BTreeMap<String, Vec<(String, f64)>> {
        BTreeMap::from([("bridge.fsintpk".into(), vec![("moisture_content".into(), moisture)])])
    }

    #[test]
    fn held_interface_state_reaches_bridge_reactions_and_pressure_with_the_original_artifact() {
        scoped(|cx, gate| {
            let source = source();
            let prepare = |moisture| {
                let text = String::from_utf8(sourced_spec().bytes).unwrap().replace("state 293.15 1000",
                    &format!("state 293.15 1000\nmaterial-states 1\nmaterial-state bridge.fsintpk 1 moisture_content {moisture}"));
                let mut s = input::Spec::parse(text.as_bytes()).unwrap();
                s.bridge_interfaces.load_with_environment(s.temperature, s.omega, &s.material_states,
                    gate, |_| Ok(source.clone())).unwrap();
                s
            };
            let mut dry = compile(cx, gate, prepare(0.0), 7800.0);
            let mut moist = compile(cx, gate, prepare(0.2), 7800.0);
            let original = &dry.performance.bank().source().bank().material_bridge_interfaces()[0].interface;
            let changed = &moist.performance.bank().source().bank().material_bridge_interfaces()[0].interface;
            assert_eq!(changed.resolved().card_identity(), source.card().content_hash());
            assert_eq!(original.resolved().card_identity(), changed.resolved().card_identity());
            assert_ne!(original.resolved().identity(), changed.resolved().identity());
            assert_eq!(changed.stiffness_n_m(), 2.0*original.stiffness_n_m());
            assert_eq!(changed.damping_n_s_m(), 4.0*original.damping_n_s_m());
            for name in ["wire.fsmatpk", "body.fsmatpk"] {
                assert_eq!(dry.catalog[name].pack.content_hash(), moist.catalog[name].pack.content_hash());
            }
            let expected = pressure(&mut dry, 37);
            let changed = pressure(&mut moist, 128);
            assert!(expected.iter().zip(&changed).any(|(a,b)| (a-b).abs() > 1e-15));
            assert_eq!(changed, pressure(&mut compile(cx, gate, prepare(0.2), 7800.0), 1));
            let frame = moist.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            let mut stale = prepare(0.0);
            stale.material_states = environment(0.2);
            let mut reads = 0;
            assert!(stale.compile(cx, gate, |name| {
                reads += 1; Ok(pack(name, 7800.0, 200_000.0, None))
            }).is_err());
            assert_eq!(reads, 0);
        });
    }

    #[test]
    fn axial_bearings_consume_the_same_held_joint_law_and_reject_stale_context() {
        use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::AxialBearingGeometry;
        scoped(|_, gate| {
            let mut inputs = axial_bearings::Inputs::new(vec![axial_bearings::Declaration {
                geometry: AxialBearingGeometry { string: 0, station_fraction: 0.8, shell_component: 0, node: 0,
                    string_axis: [1.0, 0.0, 0.0], contact_arm_m: [0.0; 3], rest_extension_m: 0.0 },
                source: bridge_interfaces::Declaration { string: 0, pack: "bridge.fsintpk".into(), area_m2: 0.0001,
                    stiffness_property: "normal_stiffness".into(), viscosity_property: "normal_viscosity".into() },
            }]).unwrap();
            let source = source();
            let held = environment(0.2);
            inputs.load_with_environment(293.15, 1000.0, &held, gate, |_| Ok(source.clone())).unwrap();
            inputs.require_environment(&held).unwrap();
            let law = inputs.bindings(293.15, 1000.0).unwrap()[0].interface.clone();
            assert_eq!(law.stiffness_n_m(), 6000.0);
            assert_eq!(law.damping_n_s_m(), 0.4);
            assert_eq!(law.resolved().query().axes()["moisture_content"], 0.2);
            assert_eq!(law.resolved().card_identity(), source.card().content_hash());
            assert!(inputs.require_environment(&environment(0.0)).is_err());
            assert!(inputs.load_with_environment(293.15, 1000.0, &environment(0.3),
                gate, |_| Ok(source.clone())).is_err());
            assert_eq!(inputs.bindings(293.15, 1000.0).unwrap()[0].interface, law);
        });
    }
}
