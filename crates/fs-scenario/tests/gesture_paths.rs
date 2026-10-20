//! G0 bitwise path admission and legacy compatibility; G5 path ordering.
#[path = "support/gesture_description.rs"]
mod fixture;
use fs_blake3::DomainHasher;
use fs_qty::Dims;
use fs_scenario::assembly::{AssemblyDescription, DEFAULT_ASSEMBLY_BUDGET, ParameterKey};
use fs_scenario::gesture::*;
use fs_scenario::{Environment, Scenario};

fn track(id: &str, target: GestureTarget, value: GestureValue) -> GestureTrack {
    GestureTrack {
        id: id.into(),
        target,
        initial: value.clone(),
        events: vec![GestureEvent {
            time_s: 0.125,
            transition_s: 0.0,
            value,
        }],
    }
}
#[test]
fn gestures_bind_entity_id_paths() {
    let (d, cases) = fixture::fixture();
    let tracks = cases
        .iter()
        .enumerate()
        .map(|(i, (_, p, v))| {
            track(
                &format!("input-{i}"),
                GestureTarget::Parameter(p.clone()),
                v.clone(),
            )
        })
        .collect();
    let s = GestureSchedule::try_new_v2(1000, tracks).unwrap();
    let resolved = resolve_aliases(&s, &d).unwrap();
    assert_eq!(
        resolved.description_hash(),
        d.content_hash(DEFAULT_ASSEMBLY_BUDGET).unwrap()
    );
    for (_, p, _) in &cases {
        assert_eq!(
            &resolved.track(p).unwrap().target,
            &GestureTarget::Parameter(p.clone())
        );
    }
    assert_eq!(
        s.to_canonical_bytes(),
        resolved.schedule().to_canonical_bytes()
    );
}
#[test]
fn v2_canonical_bytes_round_trip_and_hash_is_stable() {
    let (d, cases) = fixture::fixture();
    let (_, p, v) = &cases[0];
    let s = GestureSchedule::try_new_v2(
        700,
        vec![track(
            "actuator",
            GestureTarget::Parameter(p.clone()),
            v.clone(),
        )],
    )
    .unwrap();
    let bytes = s.to_canonical_bytes();
    let decoded = GestureSchedule::from_canonical_bytes(&bytes).unwrap();
    assert_eq!(decoded.to_canonical_bytes(), bytes);
    assert_eq!(decoded.content_hash(), s.content_hash());
    assert_eq!(decoded.schema_version(), 2);
    let mut oracle = DomainHasher::new("org.frankensim.fs-scenario.gesture-schedule.v2");
    oracle.update(&bytes);
    assert_eq!(s.content_hash(), oracle.finalize());
    let mut old = DomainHasher::new(GESTURE_SCHEDULE_HASH_DOMAIN);
    old.update(&bytes);
    assert_ne!(s.content_hash(), old.finalize());
    assert!(resolve_aliases(&decoded, &d).is_ok());
    let mut extra = bytes.clone();
    extra.extend_from_slice(b"ignored\n");
    assert!(GestureSchedule::from_canonical_bytes(&extra).is_err());
    let empty = GestureSchedule::try_new_v2(1, vec![]).unwrap();
    assert_eq!(
        GestureSchedule::from_canonical_bytes(&empty.to_canonical_bytes()).unwrap(),
        empty
    );
    let duplicate = vec![
        track("a", GestureTarget::Parameter(p.clone()), v.clone()),
        track("b", GestureTarget::Parameter(p.clone()), v.clone()),
    ];
    assert!(matches!(
        GestureSchedule::try_new_v2(1, duplicate),
        Err(GestureError::DuplicateParameter { .. })
    ));
    let signed_zero = GestureSchedule::try_new_v2(
        1,
        vec![track(
            "zero",
            GestureTarget::Parameter(p.clone()),
            GestureValue::LengthM(-0.0),
        )],
    )
    .unwrap();
    let zero_replay =
        GestureSchedule::from_canonical_bytes(&signed_zero.to_canonical_bytes()).unwrap();
    assert_eq!(
        signed_zero.to_canonical_bytes(),
        zero_replay.to_canonical_bytes()
    );
    assert!(
        matches!(zero_replay.tracks()[0].initial,GestureValue::LengthM(v) if v.to_bits()==(-0.0f64).to_bits())
    );
}

#[test]
fn force_v2_round_trip_retains_units_and_interrupted_ramps() {
    let (_, cases) = fixture::fixture();
    let mut path = cases[0].1.clone();
    path.parameter = ParameterKey::new("actuator_force", path::FORCE);
    let schedule = GestureSchedule::try_new_v2(
        8,
        vec![GestureTrack {
            id: "force".into(),
            target: GestureTarget::Parameter(path),
            initial: GestureValue::ForceN(2.0),
            events: vec![
                GestureEvent {
                    time_s: 0.25,
                    transition_s: 0.5,
                    value: GestureValue::ForceN(10.0),
                },
                GestureEvent {
                    time_s: 0.5,
                    transition_s: 0.25,
                    value: GestureValue::ForceN(-2.0),
                },
            ],
        }],
    )
    .unwrap();
    let bytes = schedule.to_canonical_bytes();
    let decoded = GestureSchedule::from_canonical_bytes(&bytes).unwrap();
    assert_eq!(decoded.schema_version(), 2);
    assert_eq!(decoded.to_canonical_bytes(), bytes);
    assert_eq!(decoded.content_hash(), schedule.content_hash());
    // The second ramp starts at 6 N, where it interrupted the first one.
    // Generic mechanical force is signed; a closing-force consumer applies
    // its own nonnegative-magnitude constraint after schedule admission.
    for (tick, expected) in [
        (0, 2.0),
        (2, 2.0),
        (3, 4.0),
        (4, 6.0),
        (5, 2.0),
        (6, -2.0),
        (8, -2.0),
    ] {
        assert_eq!(decoded.sample("force", tick).unwrap(), expected);
        assert_eq!(
            decoded.sample_value("force", tick).unwrap(),
            GestureValue::ForceN(expected)
        );
    }
}

#[test]
fn closing_force_refuses_wrong_dimensions_and_tension_payloads() {
    let (_, cases) = fixture::fixture();
    let mut path = cases[0].1.clone();
    path.parameter = ParameterKey::new("closing_force", path::FORCE);
    assert!(
        GestureSchedule::try_new_v2(
            8,
            vec![track(
                "key",
                GestureTarget::Parameter(path.clone()),
                GestureValue::ForceN(0.5)
            )],
        )
        .is_ok()
    );
    for value in [GestureValue::PressurePa(0.5), GestureValue::Fraction(0.5)] {
        assert!(matches!(
            GestureSchedule::try_new_v2(8, vec![track("key", GestureTarget::Parameter(path.clone()), value)]),
            Err(GestureError::DimensionMismatch { path, .. }) if path.ends_with("closing_force")
        ));
    }
    // Tension has the same SI dimensions but addresses a different physical
    // quantity; a key force must not silently become a string-state command.
    assert!(matches!(
        GestureSchedule::try_new_v2(8, vec![track("key", GestureTarget::Parameter(path.clone()), GestureValue::TensionN(0.5))]),
        Err(GestureError::ParameterValueMismatch { path }) if path.ends_with("closing_force")
    ));
    path.parameter.dims = path::PRESSURE;
    assert!(matches!(
        GestureSchedule::try_new_v2(8, vec![track("key", GestureTarget::Parameter(path), GestureValue::ForceN(0.5))]),
        Err(GestureError::DimensionMismatch { path, .. }) if path.ends_with("closing_force")
    ));
}

#[test]
fn legacy_targets_resolve_as_aliases() {
    let (d, cases) = fixture::fixture();
    // All sixteen variants (and both declared hole indices) are tested alone:
    // all alias addresses and values are independent fixture declarations.
    for (alias, path, value) in &cases {
        let v1 =
            GestureSchedule::try_new(1000, vec![track("legacy", alias.clone(), value.clone())])
                .unwrap();
        let bytes = v1.to_canonical_bytes();
        let hash = v1.content_hash();
        let decoded = GestureSchedule::from_canonical_bytes(&bytes).unwrap();
        assert_eq!(decoded.content_hash(), hash);
        assert_eq!(decoded.to_canonical_bytes(), bytes);
        assert_eq!(decoded.schema_version(), 1);
        let resolved = resolve_aliases(&decoded, &d).unwrap();
        assert_eq!(resolved.source_hash(), hash);
        assert_eq!(
            resolved.track(path).unwrap().target,
            GestureTarget::Parameter(path.clone())
        );
    }
    // Two hole declarations: index count-1 is valid; count is not inferred.
    let last = GestureSchedule::try_new(
        1,
        vec![track(
            "last",
            GestureTarget::ToneHoleSigma { hole: 1 },
            GestureValue::Fraction(1.0),
        )],
    )
    .unwrap();
    assert!(resolve_aliases(&last, &d).is_ok());
    let past = GestureSchedule::try_new(
        1,
        vec![track(
            "past",
            GestureTarget::ToneHoleSigma { hole: 2 },
            GestureValue::Fraction(1.0),
        )],
    )
    .unwrap();
    assert!(
        matches!(resolve_aliases(&past,&d),Err(GestureError::AliasUnresolved {target,index:2}) if target=="open_fraction")
    );
}
#[test]
fn dimension_mismatch_refuses_by_name() {
    let (d, cases) = fixture::fixture();
    let (_, path, _) = &cases[0];
    assert!(
        matches!(GestureSchedule::try_new_v2(1,vec![track("wrong-unit",GestureTarget::Parameter(path.clone()),GestureValue::PressurePa(10.0))]),Err(GestureError::DimensionMismatch {path:p,..}) if p.contains("inserted_length"))
    );
    let mut input = d.input().clone();
    input.configurations[0].assignments[0].0.parameter.name = "custom_length".into();
    let scenario = Scenario::new(
        "test-fixture:typed-gesture-addresses",
        1,
        Environment::earth_lab(),
    );
    let d = AssemblyDescription::new(&scenario, input, DEFAULT_ASSEMBLY_BUDGET).unwrap();
    let mut p = path.clone();
    p.parameter.name = "custom_length".into();
    p.parameter.dims = Dims::NONE;
    let s = GestureSchedule::try_new_v2(
        1,
        vec![track(
            "wrong-declaration",
            GestureTarget::Parameter(p),
            GestureValue::Fraction(0.5),
        )],
    )
    .unwrap();
    assert!(
        matches!(resolve_aliases(&s,&d),Err(GestureError::DimensionMismatch {path:p,..}) if p.contains("custom_length"))
    );
    let (_, angle, _) = cases
        .iter()
        .find(|(_, p, _)| p.parameter.name == "jet_angle")
        .unwrap();
    assert!(
        matches!(GestureSchedule::try_new_v2(1,vec![track("fraction-is-not-angle",GestureTarget::Parameter(angle.clone()),GestureValue::Fraction(0.5))]),Err(GestureError::ParameterValueMismatch {path}) if path.ends_with("jet_angle"))
    );
    let custom = fs_scenario::assembly::ParameterPath {
        entity: angle.entity,
        parameter: ParameterKey::new("custom_dimensionless", Dims::NONE),
    };
    let mut mixed = track(
        "mixed-kinds",
        GestureTarget::Parameter(custom),
        GestureValue::Fraction(0.5),
    );
    mixed.events[0].value = GestureValue::AngleRad(0.1);
    assert!(
        matches!(GestureSchedule::try_new_v2(1,vec![mixed]),Err(GestureError::ParameterValueMismatch {path}) if path.ends_with("custom_dimensionless"))
    );
}
#[test]
fn unknown_parameter_refuses() {
    let (d, cases) = fixture::fixture();
    let mut p = cases[0].1.clone();
    p.parameter = ParameterKey::new("undeclared", Dims::NONE);
    let s = GestureSchedule::try_new_v2(
        1,
        vec![track(
            "unknown",
            GestureTarget::Parameter(p),
            GestureValue::Fraction(0.5),
        )],
    )
    .unwrap();
    assert!(
        matches!(resolve_aliases(&s,&d),Err(GestureError::UnknownParameter {path}) if path.contains("undeclared"))
    );
}
#[test]
fn path_order_is_deterministic() {
    let (d, cases) = fixture::fixture();
    let mut tracks: Vec<_> = cases
        .iter()
        .enumerate()
        .map(|(i, (_, p, v))| {
            track(
                &format!("input-{i}"),
                GestureTarget::Parameter(p.clone()),
                v.clone(),
            )
        })
        .collect();
    let a = GestureSchedule::try_new_v2(1000, tracks.clone()).unwrap();
    tracks.reverse();
    let b = GestureSchedule::try_new_v2(1000, tracks).unwrap();
    assert_eq!(a.to_canonical_bytes(), b.to_canonical_bytes());
    assert_eq!(a.content_hash(), b.content_hash());
    let paths: Vec<_> = a
        .tracks()
        .iter()
        .map(|t| match &t.target {
            GestureTarget::Parameter(p) => (p.entity.token(), p.parameter.name.clone()),
            _ => panic!("v2 alias"),
        })
        .collect();
    assert!(paths.windows(2).all(|w| w[0] < w[1]));
    // A display rename cannot retarget a persistent identity.
    let mut input = d.input().clone();
    input
        .catalog
        .rename(cases[0].1.entity, "renamed-valve")
        .unwrap();
    let scenario = Scenario::new(
        "test-fixture:typed-gesture-addresses",
        1,
        Environment::earth_lab(),
    );
    let renamed = AssemblyDescription::new(&scenario, input, DEFAULT_ASSEMBLY_BUDGET).unwrap();
    assert_eq!(
        resolve_aliases(&a, &renamed)
            .unwrap()
            .schedule()
            .content_hash(),
        a.content_hash()
    );
}

#[test]
fn v1_example_corpus_hashes_are_unchanged() {
    // Recorded with independent BLAKE3 derive-key hashing of the exact v1
    // schedule payloads at origin/main 3333d2b9f4aa25170d5f58346995eb84ef980786.
    // Its example tree equals the .171 starting HEAD db1e6f7df. Containers
    // without a GestureSchedule use their own force/event grammars.
    let pins = [
        (
            "bow-phrases.gesture",
            "71dd520e9d42b51a45de2af819a1c134c9ab21776db562e3ceb4dea9c155dabb",
        ),
        (
            "bowed-string.performance",
            "c88a215290ae49510a2d365492807952a2a9cb8c0deb487687e35454a256f4fd",
        ),
        (
            "flue_contraction.gesture",
            "47b72bf15170a7f381c1df1f7d591689ac924045a682b29e8c261e6797503d3c",
        ),
        (
            "material_strings/bonded-pluck.material-strings",
            "070ceaa6ba96dd450cf3ccab85ad8bc41d8023127dbb363cca5711d3830afc92",
        ),
        (
            "material_strings/pluck-and-bow.material-strings",
            "5bebc509405350b8c316b529b96d3faa7fc2253d32c48f871fb1e5d28c2c003b",
        ),
        (
            "material_strings/pluck.material-strings",
            "070ceaa6ba96dd450cf3ccab85ad8bc41d8023127dbb363cca5711d3830afc92",
        ),
        (
            "plate-valve-branched.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-exterior.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-forced.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-mixture-co2.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-mixture-low-pressure.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-mixture.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-radiating.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-regional-gas.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve-viscothermal.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "plate-valve.performance",
            "5a95a9374103399d1fa60b572fa94f370ffebd77a626078aa0dcfb497ffc97b0",
        ),
        (
            "reed-baffled-exterior.performance",
            "11226b35ceddc68fda72176ca8d01fc2ef9de3329b36e0d899f52ab020488933",
        ),
        (
            "reed-baffled-mixture.performance",
            "11226b35ceddc68fda72176ca8d01fc2ef9de3329b36e0d899f52ab020488933",
        ),
        (
            "reed-duct.performance",
            "11226b35ceddc68fda72176ca8d01fc2ef9de3329b36e0d899f52ab020488933",
        ),
        (
            "reed-flanged-baffled.performance",
            "11226b35ceddc68fda72176ca8d01fc2ef9de3329b36e0d899f52ab020488933",
        ),
        (
            "reed-mixture-co2.performance",
            "d28225d56b1b1d26f3e99aabdb7003341f40f858f0305d53a7ad31b696a125e4",
        ),
        (
            "reed-mixture-low-pressure.performance",
            "d28225d56b1b1d26f3e99aabdb7003341f40f858f0305d53a7ad31b696a125e4",
        ),
        (
            "reed-mixture.performance",
            "d28225d56b1b1d26f3e99aabdb7003341f40f858f0305d53a7ad31b696a125e4",
        ),
        (
            "reed-pressure.gesture",
            "0b645dcc15decdabef1505b11124bdd3aa75a076a89f80072ffb1f7d7c6cba4c",
        ),
        (
            "reed-regional-mixture.performance",
            "11226b35ceddc68fda72176ca8d01fc2ef9de3329b36e0d899f52ab020488933",
        ),
        (
            "two_sided_reed.gesture",
            "c49ec276f79af8dc840e9635083fd730fd1c3e0c4ed6f7c6316f38cf5042f756",
        ),
    ];
    fn visit(dir: &std::path::Path, files: &mut Vec<std::path::PathBuf>) {
        for e in std::fs::read_dir(dir).unwrap() {
            let e = e.unwrap();
            if e.file_type().unwrap().is_dir() {
                visit(&e.path(), files);
            } else if matches!(
                e.path().extension().and_then(|x| x.to_str()),
                Some("gesture" | "performance" | "material-strings")
            ) {
                files.push(e.path());
            }
        }
    }
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../fs-couple/examples");
    let mut files = Vec::new();
    visit(&root, &mut files);
    files.sort();
    let mut decoded = 0;
    for path in files {
        let bytes = std::fs::read(&path).unwrap();
        let magic = b"frankensim-gesture-schedule-v1\n";
        let Some(offset) = bytes.windows(magic.len()).position(|w| w == magic) else {
            assert_ne!(
                path.extension().unwrap(),
                "gesture",
                "standalone gesture must decode: {}",
                path.display()
            );
            continue;
        };
        let name = path
            .strip_prefix(&root)
            .unwrap()
            .to_str()
            .unwrap()
            .replace('\\', "/");
        let expected = pins
            .iter()
            .find(|(n, _)| *n == name)
            .unwrap_or_else(|| panic!("new v1 corpus entry needs a pin: {name}"));
        let s = GestureSchedule::from_canonical_bytes(&bytes[offset..]).unwrap();
        assert_eq!(s.schema_version(), 1);
        assert_eq!(s.to_canonical_bytes(), bytes[offset..], "{name}");
        assert_eq!(s.content_hash().to_hex(), expected.1, "{name}");
        decoded += 1;
    }
    assert_eq!(decoded, pins.len());
}
