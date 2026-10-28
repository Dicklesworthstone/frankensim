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
