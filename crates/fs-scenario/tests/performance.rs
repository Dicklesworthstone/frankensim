//! G0 exact score interpretation/refusals; G5 canonical replay. No recordings.
use fs_blake3::hash_bytes;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::*;
use fs_scenario::gesture::{GestureTarget, GestureValue};
use fs_scenario::performance::*;
use fs_scenario::{EntityCatalog, EntityDeclaration, Environment, Scenario};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
fn limits() -> SmfLimits {
    SmfLimits {
        max_input_bytes: 4096,
        max_tracks: 4,
        max_events: 32,
        max_control_ticks: 100_000,
    }
}
fn fixture() -> (AssemblyDescription, ImportPolicy, DynamicsTable) {
    let s = Scenario::new("test-fixture:authored-midi", 1, Environment::earth_lab());
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("test-fixture:player"))
        .unwrap();
    let target = ParameterPath {
        entity: root,
        parameter: ParameterKey::new("player_pressure", PA),
    };
    let desc = AssemblyDescription::new(
        &s,
        AssemblyInput {
            root,
            catalog,
            parts: vec![],
            interfaces: vec![],
            configurations: vec![ConfigurationState {
                name: "open".into(),
                assignments: vec![],
                provenance: StateProvenance::Declared {
                    rationale: "owner-authored test fixture, no real-instrument claim".into(),
                },
            }],
            observers: vec![],
            reference_state: Some(ReferenceState::Assembled),
        },
        DEFAULT_ASSEMBLY_BUDGET,
    )
    .unwrap();
    let source = hash_bytes(b"owner-authored test fixture v1");
    let policy = ImportPolicy {
        chart: NoteChart {
            source,
            entries: vec![NoteChartEntry {
                channel: 0,
                note: 60,
                configuration: "open".into(),
            }],
        },
        velocities: VelocityBins {
            source,
            bins: vec![VelocityBin {
                first: 1,
                last: 127,
                marking: "mf".into(),
            }],
        },
        bindings: vec![CommandBinding {
            path: target.clone(),
            track_id: "blow".into(),
            target: GestureTarget::BlowingPressure,
            initial: GestureValue::PressurePa(0.0),
        }],
        controllers: vec![],
        control_rate_hz: 1000,
        initial_tempo_us: None,
        observation_prior: hash_bytes(b"prior:fixed listener test fixture"),
    };
    let table = DynamicsTable {
        source,
        rows: vec![DynamicsRow {
            marking: "mf".into(),
            controls: vec![DynamicControl {
                target,
                domain: (QtyAny::new(1000.0, PA), QtyAny::new(3000.0, PA)),
                choice: ControlChoice::Fixed {
                    value: QtyAny::new(2000.0, PA),
                    source,
                },
            }],
        }],
    };
    (desc, policy, table)
}
fn file(format: u16, division: u16, tracks: &[Vec<u8>]) -> Vec<u8> {
    let mut b = b"MThd\0\0\0\x06".to_vec();
    b.extend(format.to_be_bytes());
    b.extend((tracks.len() as u16).to_be_bytes());
    b.extend(division.to_be_bytes());
    for track in tracks {
        b.extend(b"MTrk");
        b.extend((track.len() as u32).to_be_bytes());
        b.extend(track);
    }
    b
}
fn one() -> Vec<u8> {
    file(
        0,
        480,
        &[vec![
            0, 0xff, 0x51, 3, 7, 0xa1, 0x20, 0, 0x90, 60, 64, 0x83, 0x60, 60, 0, 0, 0xff, 0x2f, 0,
        ]],
    )
}
fn import(bytes: &[u8]) -> Result<Performance, PerformanceError> {
    let (d, p, t) = fixture();
    import_smf(bytes, &d, &p, &t, limits())
}

#[test]
fn pf_001_tempo_map_ticks_exact() {
    // PPQ=3: (1*500001+2*250001)/3 microseconds. At 1 MHz:
    // ceil(1000003/3)=333335, rounded ONCE across the tempo boundary.
    let conductor = vec![
        0, 0xff, 0x51, 3, 7, 0xa1, 0x21, 1, 0xff, 0x51, 3, 3, 0xd0, 0x91, 2, 0xff, 0x2f, 0,
    ];
    let notes = vec![0, 0x90, 60, 80, 3, 0x80, 60, 0, 0, 0xff, 0x2f, 0];
    let bytes = file(1, 3, &[conductor, notes]);
    let (d, mut p, t) = fixture();
    p.control_rate_hz = 1_000_000;
    let a = import_smf(
        &bytes,
        &d,
        &p,
        &t,
        SmfLimits {
            max_control_ticks: 333335,
            ..limits()
        },
    )
    .unwrap();
    assert_eq!(a.notes()[0].on_tick, 0);
    assert_eq!(a.notes()[0].off_tick, 333335);
    assert!(matches!(
        import(&file(
            0,
            480,
            &[vec![0, 0x90, 60, 64, 1, 60, 0, 0, 0xff, 0x2f, 0]]
        )),
        Err(PerformanceError::MissingTempo)
    ));
}
#[test]
fn pf_002_running_status_and_velocity_zero_note_off() {
    let p = import(&one()).unwrap();
    assert_eq!(p.notes().len(), 1);
    assert_eq!(p.notes()[0].off_tick, 500);
    assert_eq!(p.notes()[0].velocity, 64);
}
#[test]
fn pf_003_overlap_is_declared_legato_or_refused() {
    let b = file(
        0,
        480,
        &[vec![
            0, 0xff, 0x51, 3, 7, 0xa1, 0x20, 0, 0x90, 60, 64, 1, 61, 64, 1, 60, 0, 1, 61, 0, 0,
            0xff, 0x2f, 0,
        ]],
    );
    assert!(matches!(
        import(&b),
        Err(PerformanceError::OverlappingNote {
            channel: 0,
            note: 61
        })
    ));
}
#[test]
fn pf_004_smpte_division_handled_by_name() {
    let b = file(0, 0xe728, &[vec![0, 0xff, 0x2f, 0]]);
    assert!(matches!(
        import(&b),
        Err(PerformanceError::SmpteDivision { division: 0xe728 })
    ));
}
#[test]
fn pf_005_smf_limits_cap_and_cap_plus_one() {
    let (d, p, t) = fixture();
    let b = one();
    let exact = SmfLimits {
        max_input_bytes: b.len(),
        max_tracks: 1,
        max_events: 4,
        max_control_ticks: 500,
    };
    assert!(import_smf(&b, &d, &p, &t, exact).is_ok());
    for (lim, name) in [
        (
            SmfLimits {
                max_input_bytes: b.len() - 1,
                ..exact
            },
            "MIDI input bytes",
        ),
        (
            SmfLimits {
                max_tracks: 0,
                ..exact
            },
            "MIDI tracks",
        ),
        (
            SmfLimits {
                max_events: 3,
                ..exact
            },
            "MIDI events",
        ),
        (
            SmfLimits {
                max_control_ticks: 499,
                ..exact
            },
            "control ticks",
        ),
    ] {
        assert!(
            matches!(import_smf(&b,&d,&p,&t,lim),Err(PerformanceError::Limit { resource,.. }) if resource==name)
        );
    }
    // Explicit format 1 cap (two complete tracks), not merely a header count.
    let two = file(
        1,
        480,
        &[
            vec![0, 0xff, 0x51, 3, 7, 0xa1, 0x20, 0, 0xff, 0x2f, 0],
            vec![0, 0x90, 60, 64, 1, 60, 0, 0, 0xff, 0x2f, 0],
        ],
    );
    assert!(
        import_smf(
            &two,
            &d,
            &p,
            &t,
            SmfLimits {
                max_tracks: 2,
                ..limits()
            }
        )
        .is_ok()
    );
    assert!(matches!(
        import_smf(
            &two,
            &d,
            &p,
            &t,
            SmfLimits {
                max_tracks: 1,
                ..limits()
            }
        ),
        Err(PerformanceError::Limit {
            resource: "MIDI tracks",
            ..
        })
    ));
}
#[test]
fn pf_006_unknown_note_refuses_no_frequency_fallback() {
    let (d, mut p, t) = fixture();
    p.chart.entries.clear();
    assert!(matches!(
        import_smf(&one(), &d, &p, &t, limits()),
        Err(PerformanceError::UnknownNote {
            channel: 0,
            note: 60
        })
    ));
}
#[test]
fn pf_007_missing_dynamic_row_refuses() {
    let (d, p, mut t) = fixture();
    t.rows.clear();
    assert!(
        matches!(import_smf(&one(),&d,&p,&t,limits()),Err(PerformanceError::MissingDynamicRow { marking }) if marking=="mf")
    );
}
#[test]
fn pf_008_lowered_schedule_canonical_bytes_stable() {
    let (d, p, t) = fixture();
    let a = import_smf(&one(), &d, &p, &t, limits()).unwrap();
    let wire = a.to_bytes();
    let b = Performance::from_bytes(&wire, &d, &p, &t, limits()).unwrap();
    assert_eq!(a.content_hash(), b.content_hash());
    assert_eq!(wire, b.to_bytes());
    let a = lower(&a, &d, &p, &t).unwrap();
    let b = lower(&b, &d, &p, &t).unwrap();
    assert_eq!(
        a.schedule.to_canonical_bytes(),
        b.schedule.to_canonical_bytes()
    );
    assert_eq!(a.schedule.sample("blow", 0).unwrap(), 2000.0);
    assert_eq!(a.schedule.sample("blow", 499).unwrap(), 2000.0);
    assert_eq!(a.schedule.sample("blow", 500).unwrap(), 0.0);
    let mut wrong = p.clone();
    wrong.observation_prior = hash_bytes(b"different prior");
    assert!(matches!(
        Performance::from_bytes(&wire, &d, &wrong, &t, limits()),
        Err(PerformanceError::HashMismatch { input: "policy" })
    ));
    for prefix in 0..wire.len() {
        assert!(Performance::from_bytes(&wire[..prefix], &d, &p, &t, limits()).is_err());
    }
}
#[test]
fn pf_009_commands_cannot_target_physical_state() {
    let (d, mut p, t) = fixture();
    let mut catalog = d.input().catalog.clone();
    let airway = catalog
        .declare(EntityDeclaration::part(d.input().root, "airway"))
        .unwrap();
    p.bindings[0].path.entity = airway;
    assert!(
        matches!(import_smf(&one(),&d,&p,&t,limits()),Err(PerformanceError::PhysicalStateCommand { path }) if path.entity==airway)
    );
    p.bindings[0].path.entity = d.input().root;
    p.bindings[0].target = GestureTarget::PreStress;
    assert!(matches!(
        import_smf(&one(), &d, &p, &t, limits()),
        Err(PerformanceError::PhysicalStateCommand { .. })
    ));
}
#[test]
fn pf_010_controllers_require_named_rules_and_units() {
    let b = file(
        0,
        480,
        &[vec![
            0, 0xff, 0x51, 3, 7, 0xa1, 0x20, 0, 0xb0, 2, 127, 0, 0xff, 0x2f, 0,
        ]],
    );
    let (d, mut p, t) = fixture();
    assert!(matches!(
        import_smf(&b, &d, &p, &t, limits()),
        Err(PerformanceError::UnknownController {
            channel: 0,
            controller: 2
        })
    ));
    let mut target = p.bindings[0].path.clone();
    target.parameter = ParameterKey::new("pedal", Dims::NONE);
    p.bindings.push(CommandBinding {
        path: target.clone(),
        track_id: "pedal".into(),
        target: GestureTarget::SustainPedal,
        initial: GestureValue::Fraction(0.0),
    });
    p.controllers.push(ControllerRule {
        channel: 0,
        controller: 2,
        target,
        range: (QtyAny::dimensionless(0.0), QtyAny::dimensionless(1.0)),
    });
    let perf = import_smf(&b, &d, &p, &t, limits()).unwrap();
    let out = lower(&perf, &d, &p, &t).unwrap();
    assert_eq!(out.schedule.sample("pedal", 0).unwrap(), 1.0);
}
#[test]
fn pf_011_configuration_assignments_and_free_controls_are_explicit() {
    let (d, p, mut t) = fixture();
    t.rows[0].controls[0].choice = ControlChoice::PerNote {
        objective: DeclaredObjective::MinimumThresholdPressure,
    };
    let mut perf = import_smf(&one(), &d, &p, &t, limits()).unwrap();
    assert!(matches!(
        lower(&perf, &d, &p, &t),
        Err(PerformanceError::UnresolvedControl { index: 0 })
    ));
    assert!(
        perf.resolve(0, QtyAny::new(2000.0, Dims::NONE), hash_bytes(b"bad units"))
            .is_err()
    );
    assert!(
        perf.resolve(0, QtyAny::new(999.0, PA), hash_bytes(b"outside"))
            .is_err()
    );
    perf.resolve(
        0,
        QtyAny::new(2000.0, PA),
        hash_bytes(b"explicit test evaluation"),
    )
    .unwrap();
    let out = lower(&perf, &d, &p, &t).unwrap();
    assert_eq!(out.provenance.per_note[0], perf.per_note()[0]);
    let mut input = d.input().clone();
    let key = ParameterPath {
        entity: input.root,
        parameter: ParameterKey::new("hole_open_fraction", Dims::NONE),
    };
    input.configurations[0]
        .assignments
        .push((key.clone(), QtyAny::dimensionless(0.75)));
    let scenario = Scenario::new("test-fixture:authored-midi", 1, Environment::earth_lab());
    let d = AssemblyDescription::new(&scenario, input, DEFAULT_ASSEMBLY_BUDGET).unwrap();
    let mut p = p;
    p.bindings.push(CommandBinding {
        path: key,
        track_id: "hole".into(),
        target: GestureTarget::ToneHoleSigma { hole: 0 },
        initial: GestureValue::Fraction(0.0),
    });
    t.rows[0].controls[0].choice = ControlChoice::Fixed {
        value: QtyAny::new(2000.0, PA),
        source: hash_bytes(b"authored"),
    };
    let perf = import_smf(&one(), &d, &p, &t, limits()).unwrap();
    assert_eq!(
        lower(&perf, &d, &p, &t)
            .unwrap()
            .schedule
            .sample("hole", 0)
            .unwrap(),
        0.75
    );
}
