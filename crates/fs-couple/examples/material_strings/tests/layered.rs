//! Exercise the public V2 input through the original material and sound owners.
use super::*;
use fs_material::section::ConcentricKelvinVoigtSection;

const WIRE: &str = "string wire.fsmatpk 0.5 0.0003 20 4 3 0.8 3000 0.1";
fn v2(section: &str, support: &str) -> input::Spec {
    let text = String::from_utf8(spec().bytes).unwrap()
        .replace("frankensim-material-strings-v1", "frankensim-material-strings-v2")
        .replace(WIRE, &format!("string wire.fsmatpk 0.5 0.0003 {support} 4 3 0.8 3000 0.1\n{section}"));
    input::Spec::parse(text.as_bytes()).unwrap()
}
fn bonded(support: &str) -> input::Spec {
    v2("section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003", support)
}
fn bowed() -> input::Spec {
    let old = bonded("tension 20");
    let schedule = GestureSchedule::try_new(700, vec![GestureTrack {
        id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
        initial: GestureValue::Bow { velocity_m_per_s: 0.0, normal_force_n: 0.0, station: 0.2 },
        events: vec![GestureEvent { time_s: 0.0, transition_s: 0.002,
            value: GestureValue::Bow { velocity_m_per_s: 0.01, normal_force_n: 0.01, station: 0.2 } }],
    }]).unwrap();
    let text = String::from_utf8(old.bytes).unwrap();
    let mut bytes = format!("{}\nschedule\n", text.split_once("\nschedule\n").unwrap().0).into_bytes();
    bytes.extend(schedule.to_canonical_bytes());
    input::Spec::parse(&bytes).unwrap()
}
// Synthetic external source mutation: only the named sleeve property changes.
fn sleeve(key: &str, value: f64) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let provenance = Provenance {
        source: "synthetic independent sleeve property".into(),
        license: "CC0-1.0".into(),
        artifact: Some(fs_blake3::ContentHash([9;32])),
    };
    let obs = claims.register_observation(ObservationDataset {
        specimen: "sleeve.fsmatpk".into(),
        method: "synthetic independent sleeve property".into(),
        artifact: fs_blake3::ContentHash([9;32]),
        caveats: "synthetic test fixture".into(),
        provenance: provenance.clone(),
    }).unwrap();
    for (name, dims, original) in [("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, 200e9), (ETA, DynViscosity::DIMS, 10.0),
        (SHEAR, Pressure::DIMS, 80e9), (SHEAR_ETA, DynViscosity::DIMS, 1000.0)] {
        let value = if name == key { value } else { original };
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, 200_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![obs.clone()], provenance: provenance.clone() }).unwrap();
    }
    let pack = NormalizedPack::new("sleeve.fsmatpk", "synthetic-test-v1", fs_blake3::ContentHash([9;32]),
        "synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}
fn layered_ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, key: &str, value: f64) -> Ready {
    spec.compile(cx, gate, |name| Ok(if name == "sleeve.fsmatpk" { sleeve(key, value) }
        else { pack(name, 7800.0, 200_000.0, None) })).unwrap()
}
fn source(ready: &Ready) -> &MaterialStringSource {
    &ready.performance.bank().source().bank().strings()[0].source
}
fn section(ready: &Ready) -> &ConcentricKelvinVoigtSection {
    let s = match source(ready) { MaterialStringSource::WithTorsion(t) => t.bending(), s => s };
    let MaterialStringSource::Concentric { section, .. } = s else { panic!("bonded source retained") };
    section
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-30), "{a} != {b}");
}

fn explicit_sources(pair: &str) -> input::Spec {
    let fields = match pair {
        "young-poisson" => "E-source nu-source",
        "young-shear" => "E-source G-source",
        "bulk-shear" => "K-source G-source",
        _ => unreachable!(),
    };
    let text = String::from_utf8(bowed().bytes).unwrap().replace(
        "layer sleeve.fsmatpk 0.0003",
        &format!(
            "layer sleeve.fsmatpk 0.0003\nelastic-sources 2\n\
            elastic-source rho-source young-poisson E-source nu-source\n\
            elastic-source rho-source {pair} {fields}"
        ),
    );
    input::Spec::parse(text.as_bytes()).unwrap()
}

// Equivalent independent source pairs: E=180 GPa, nu=1/8, K=G=80 GPa.
// Only the selected pair exists in each pack; neither canonical density nor
// fictitious measured Young/shear claims can rescue an unwired consumer.
fn independent_pack(name: &str, pair: &str, mutation: Option<(&str, f64)>) -> NormalizedPack {
    let mut values = vec![
        ("rho-source", Density::DIMS, 7800.0),
        (ETA, DynViscosity::DIMS, 10.0),
        (SHEAR_ETA, DynViscosity::DIMS, 1000.0),
    ];
    values.extend(match pair {
        "young-poisson" => [
            ("E-source", Pressure::DIMS, 180e9),
            ("nu-source", Dims::NONE, 0.125),
        ],
        "young-shear" => [
            ("E-source", Pressure::DIMS, 180e9),
            ("G-source", Pressure::DIMS, 80e9),
        ],
        "bulk-shear" => [
            ("K-source", Pressure::DIMS, 80e9),
            ("G-source", Pressure::DIMS, 80e9),
        ],
        _ => unreachable!(),
    });
    let mut claims = ClaimSet::new();
    let provenance = Provenance {
        source: "manufactured independent elastic sources".into(),
        license: "CC0-1.0".into(),
        artifact: Some(fs_blake3::ContentHash([5; 32])),
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: name.into(),
            method: "manufactured algebraic equivalence".into(),
            artifact: fs_blake3::ContentHash([5; 32]),
            caveats: "synthetic; not measured metal data".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (key, dims, original) in values {
        let value = mutation
            .filter(|(selected, _)| *selected == key)
            .map_or(original, |(_, value)| value);
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, 200_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation.clone()],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let pack = NormalizedPack::new(
        name,
        "synthetic-independent-pair-v1",
        fs_blake3::ContentHash([5; 32]),
        "synthetic redistribution permitted",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

fn independent_ready(
    cx: &Cx<'_>,
    gate: &CancelGate,
    pair: &str,
    mutation: Option<(&str, f64)>,
) -> Ready {
    explicit_sources(pair)
        .compile(cx, gate, |name| {
            Ok(match name {
                "wire.fsmatpk" => independent_pack(name, "young-poisson", None),
                "sleeve.fsmatpk" => independent_pack(name, pair, mutation),
                _ => pack(name, 7800.0, 200_000.0, None),
            })
        })
        .unwrap()
}

#[test]
fn independent_elastic_source_rows_are_bounded_explicit_and_complete() {
    for pair in ["young-poisson", "young-shear", "bulk-shear"] {
        let spec = explicit_sources(pair);
        assert_eq!(spec.strings[0].elastic_sources.as_ref().unwrap().len(), 2);
    }
    let text = String::from_utf8(explicit_sources("bulk-shear").bytes).unwrap();
    for (from, to) in [
        ("elastic-sources 2", "elastic-sources 0"),
        ("elastic-sources 2", "elastic-sources 1"),
        ("elastic-sources 2", "elastic-sources 17"),
        (
            "bulk-shear K-source G-source",
            "automatic K-source G-source",
        ),
        (
            "bulk-shear K-source G-source",
            "bulk-shear K-source K-source",
        ),
        (
            "bulk-shear K-source G-source",
            "bulk-shear rho-source G-source",
        ),
        (
            "bulk-shear K-source G-source",
            "bulk-shear K-source G-source extra",
        ),
        (
            "section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003",
            "section homogeneous",
        ),
    ] {
        assert!(
            input::Spec::parse(text.replace(from, to).as_bytes()).is_err(),
            "accepted {to}"
        );
    }
    scoped(|cx, gate| {
        for mutation in [None, Some(("K-source", -1.0))] {
            assert!(
                explicit_sources("bulk-shear")
                    .compile(cx, gate, |name| Ok(match name {
                        "wire.fsmatpk" => independent_pack(name, "young-poisson", None),
                        // Missing selected K must not switch to the available E/G pair.
                        "sleeve.fsmatpk" => independent_pack(
                            name,
                            if mutation.is_none() {
                                "young-shear"
                            } else {
                                "bulk-shear"
                            },
                            mutation
                        ),
                        _ => pack(name, 7800.0, 200_000.0, None),
                    }))
                    .is_err()
            );
        }
    });
}

#[test]
fn data_only_elastic_pairs_reach_bending_torsion_and_coupled_pressure() {
    scoped(|cx, gate| {
        let mut baseline = independent_ready(cx, gate, "bulk-shear", None);
        let expected = pressure(&mut baseline, 37);
        assert!(expected.iter().any(|value| *value != 0.0));
        for pair in ["young-poisson", "young-shear", "bulk-shear"] {
            let mut ready = independent_ready(cx, gate, pair, None);
            let elastic = section(&ready).elastic();
            assert_eq!(
                elastic.isotropic_parameters(),
                ready.spec.strings[0].elastic_sources.as_deref()
            );
            for layer in elastic.layers() {
                for absent in ["density", "young_modulus", "shear_modulus"] {
                    assert!(layer.material.property(absent).is_none());
                }
            }
            let MaterialStringSource::WithTorsion(twist) = source(&ready) else {
                panic!("retained torsion")
            };
            assert_eq!(twist.section().shear().shear_property(), None);
            assert_eq!(
                pressure(&mut ready, 1),
                expected,
                "equivalent {pair} changed sound"
            );
        }
        let core_a = core::f64::consts::PI * 0.0002_f64.powi(2);
        let sleeve_a = core::f64::consts::PI * (0.0003_f64.powi(2) - 0.0002_f64.powi(2));
        let core_i = core::f64::consts::PI * 0.0002_f64.powi(4) / 4.0;
        let sleeve_i = core::f64::consts::PI * (0.0003_f64.powi(4) - 0.0002_f64.powi(4)) / 4.0;
        for (key, value) in [
            ("rho-source", 15600.0),
            ("K-source", 40e9),
            ("G-source", 40e9),
        ] {
            let mut changed = independent_ready(cx, gate, "bulk-shear", Some((key, value)));
            assert_eq!(changed.spec.bytes, baseline.spec.bytes);
            assert_eq!(
                changed.catalog["wire.fsmatpk"].pack.content_hash(),
                baseline.catalog["wire.fsmatpk"].pack.content_hash()
            );
            assert_eq!(
                changed.catalog["body.fsmatpk"].pack.content_hash(),
                baseline.catalog["body.fsmatpk"].pack.content_hash()
            );
            let density = if key == "rho-source" { value } else { 7800.0 };
            let bulk = if key == "K-source" { value } else { 80e9 };
            let shear = if key == "G-source" { value } else { 80e9 };
            let young = 9.0 * bulk * shear / (3.0 * bulk + shear);
            let elastic = section(&changed).elastic();
            close(
                elastic.line_density_kg_m(),
                7800.0 * core_a + density * sleeve_a,
            );
            close(
                elastic.axial_stiffness_n(),
                180e9 * core_a + young * sleeve_a,
            );
            close(
                elastic.bending_stiffness_n_m2(),
                180e9 * core_i + young * sleeve_i,
            );
            assert_eq!(
                elastic.layers()[0].material,
                section(&baseline).elastic().layers()[0].material
            );
            let MaterialStringSource::WithTorsion(twist) = source(&changed) else {
                panic!("retained torsion")
            };
            close(
                twist.section().rotary_inertia_kg_m(),
                2.0 * (7800.0 * core_i + density * sleeve_i),
            );
            close(
                twist.section().stiffness_n_m2(),
                2.0 * (80e9 * core_i + shear * sleeve_i),
            );
            assert_eq!(changed.performance.bank().cards()[0].tension_n, 20.0);
            let actual = pressure(&mut changed, 1);
            assert!(
                actual
                    .iter()
                    .zip(&expected)
                    .any(|(a, b)| (a - b).abs() > 1e-12),
                "{key} did not reach sound"
            );
            let frame = changed.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn explicit_v2_homogeneous_support_keeps_v1_mechanics_and_pcm_bit_for_bit() {
    scoped(|cx, gate| {
        let mut old = compile(cx, gate, spec(), 7800.0);
        let mut new = compile(cx, gate, v2("section homogeneous", "tension 20"), 7800.0);
        assert_eq!(old.performance.bank().cards()[0].fundamental_hz(), new.performance.bank().cards()[0].fundamental_hz());
        assert_eq!(old.performance.preload_energy_j(), new.performance.preload_energy_j());
        assert_eq!(old.catalog.len(), new.catalog.len());
        let (a, sa) = old.write(Cursor::new(Vec::new()), gate).unwrap();
        let (b, sb) = new.write(Cursor::new(Vec::new()), gate).unwrap();
        assert_eq!(a.into_inner(), b.into_inner());
        assert_eq!(sa.samples, sb.samples);
    });
}

#[test]
fn independent_sleeve_properties_reach_section_torsion_and_coupled_pressure() {
    scoped(|cx, gate| {
        let mut original = layered_ready(cx, gate, bowed(), "density", 7800.0);
        let expected = pressure(&mut original, 37);
        for (key, value) in [("density", 15600.0), ("young_modulus", 100e9),
            (ETA, 1e7), (SHEAR, 40e9), (SHEAR_ETA, 1e6)] {
            let mut changed = layered_ready(cx, gate, bowed(), key, value);
            assert_eq!(changed.spec.bytes, original.spec.bytes);
            assert_eq!(changed.catalog.len(), 3);
            assert_eq!(changed.catalog["wire.fsmatpk"].pack.content_hash(), original.catalog["wire.fsmatpk"].pack.content_hash());
            assert_eq!(changed.catalog["body.fsmatpk"].pack.content_hash(), original.catalog["body.fsmatpk"].pack.content_hash());
            let e = section(&changed).elastic();
            let core_a = core::f64::consts::PI * 0.0002_f64.powi(2);
            let sleeve_a = core::f64::consts::PI * (0.0003_f64.powi(2) - 0.0002_f64.powi(2));
            let core_i = core::f64::consts::PI * 0.0002_f64.powi(4) / 4.0;
            let sleeve_i = core::f64::consts::PI * (0.0003_f64.powi(4) - 0.0002_f64.powi(4)) / 4.0;
            close(e.line_density_kg_m(), 7800.0 * core_a + (if key == "density" {value} else {7800.0}) * sleeve_a);
            let young = if key == "young_modulus" {value} else {200e9};
            close(e.axial_stiffness_n(), 200e9 * core_a + young * sleeve_a);
            close(e.bending_stiffness_n_m2(), 200e9 * core_i + young * sleeve_i);
            close(section(&changed).viscous_stiffness_n_m2_s(), 10.0 * core_i + (if key == ETA {value} else {10.0}) * sleeve_i);
            assert_eq!(e.layers().len(), 2);
            assert_eq!(e.layers()[0].material, section(&original).elastic().layers()[0].material);
            let MaterialStringSource::WithTorsion(twist) = source(&changed) else { panic!("original torsion") };
            close(twist.section().rotary_inertia_kg_m(), 2.0 * (7800.0 * core_i + (if key == "density" {value} else {7800.0}) * sleeve_i));
            close(twist.section().stiffness_n_m2(), 2.0 * (80e9 * core_i + (if key == SHEAR {value} else {80e9}) * sleeve_i));
            assert_eq!(changed.performance.bank().cards()[0].tension_n, 20.0);
            let actual = pressure(&mut changed, 1);
            assert!(actual.iter().zip(&expected).any(|(a,b)| (a-b).abs() > 1e-12), "{key} did not reach actual output");
            let frame = changed.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn fixed_extension_uses_the_same_layer_axial_stiffness_instead_of_retaining_old_tension() {
    scoped(|cx, gate| {
        let support = "extension 0.4998 0.001";
        let a = layered_ready(cx, gate, bonded(support), "young_modulus", 200e9);
        let b = layered_ready(cx, gate, bonded(support), "young_modulus", 100e9);
        let strain = (0.5 - 0.4998) / 0.4998;
        close(a.performance.bank().cards()[0].tension_n, section(&a).elastic().axial_stiffness_n() * strain);
        close(b.performance.bank().cards()[0].tension_n, section(&b).elastic().axial_stiffness_n() * strain);
        assert!(b.performance.bank().cards()[0].tension_n < a.performance.bank().cards()[0].tension_n);
        let homogeneous = compile(cx, gate, v2("section homogeneous", support), 7800.0);
        close(homogeneous.performance.bank().cards()[0].tension_n, a.performance.bank().cards()[0].tension_n);
        for section in ["section homogeneous", "section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003"] {
            for bad in ["extension 0.5 0.001", "extension 0.6 0.001", "extension 0.4998 0.00001", "extension 0.4998 1"] {
                assert!(v2(section, bad).compile(cx, gate, |name| Ok(pack(name,7800.0,200_000.0,None))).is_err());
            }
        }
    });
}

#[test]
fn layer_receipts_and_original_surface_pluck_survive_callbacks_and_streamed_pcm() {
    scoped(|cx, gate| {
        let mut direct = layered_ready(cx, gate, bonded("tension 20"), "density", 15600.0);
        let expected = pressure(&mut direct, 1);
        for block in [37, 128] {
            let mut p = layered_ready(cx, gate, bonded("tension 20"), "density", 15600.0);
            assert_eq!(pressure(&mut p, block), expected);
            assert_eq!(section(&p), section(&direct));
        }
        let mut p = layered_ready(cx, gate, bonded("tension 20"), "density", 15600.0);
        let (wav, summary) = p.write(Cursor::new(Vec::new()), gate).unwrap();
        let (bytes, clips) = fs_couple::pcm_wav::encode_pcm16_wav(&expected, RATE, 0.01).unwrap();
        assert_eq!(wav.into_inner(), bytes); assert_eq!(summary.clipped_samples, clips as u64);
        assert_eq!(summary.samples, 256);
        assert_eq!(section(&p).elastic().layers()[1].material.property("density").unwrap().value_si(), 15600.0);
        // Two spatial regions may refer to ONE immutable pack without duplicate I/O.
        let same = v2("section bonded 0.0002 1\nlayer wire.fsmatpk 0.0003", "tension 20");
        let mut loaded = Vec::new();
        let p = same.compile(cx, gate, |name| { loaded.push(name.to_string()); Ok(pack(name,7800.0,200_000.0,None)) }).unwrap();
        assert_eq!(loaded, vec!["wire.fsmatpk".to_string(), "body.fsmatpk".to_string()]);
        assert_eq!(section(&p).elastic().layers().len(), 2);
    });
}

#[test]
fn sleeve_geometry_missing_loss_and_source_band_refuse_without_a_core_fallback() {
    let text = String::from_utf8(bonded("tension 20").bytes).unwrap();
    for (from, to) in [("bonded 0.0002 1", "bonded 0.0002 0"),
        ("bonded 0.0002 1", "bonded 0.0002 999999999"),
        ("bonded 0.0002 1", "bonded 0.0003 1"),
        ("layer sleeve.fsmatpk 0.0003", "layer sleeve.fsmatpk 0.0002"),
        ("layer sleeve.fsmatpk 0.0003", "layer sleeve.fsmatpk 0.00025"),
        ("layer sleeve.fsmatpk 0.0003", "layer sleeve.fsmatpk 0.0004"),
        ("sleeve.fsmatpk", "../sleeve.fsmatpk"), ("tension 20", "pitch 440"),
        ("frankensim-material-strings-v2", "frankensim-material-strings-v1")] {
        assert!(input::Spec::parse(text.replace(from,to).as_bytes()).is_err(), "accepted {to}");
    }
    scoped(|cx, gate| {
        for (band, omit) in [(200_000.0,Some(ETA)), (200_000.0,Some(SHEAR)),
            (200_000.0,Some(SHEAR_ETA)), (2000.0,None)] {
            assert!(bonded("tension 20").compile(cx,gate,|name| Ok(pack(name,7800.0,
                if name=="sleeve.fsmatpk" {band} else {200_000.0},
                if name=="sleeve.fsmatpk" {omit} else {None}))).is_err());
        }
        // The actual shipped V2 source is also consumed by this test target.
        let fixture = input::Spec::parse(include_bytes!("../bonded-pluck.material-strings")).unwrap();
        let mut ready = fixture.compile(cx,gate,|name| Ok(pack(name,7800.0,200_000.0,None))).unwrap();
        ready.performance.block(&mut [0.0; 128]).unwrap();
        assert_eq!(section(&ready).elastic().layers().len(), 2);
    });
}

const RATE: u32 = 48_000;
