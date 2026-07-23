//! Battery for the versioned `.fsim` project schema (bead f85xj.6.1): the
//! reference cooling project parses and is admissible, canonical bytes are
//! stable across both spellings, every mandatory-field omission is a named
//! violation, unknown fields are refused, the only default is receipted, and
//! the version-bump machinery is proven with the synthetic migration. The
//! broken-project corpus doubles as documentation of the error-message
//! quality bar: every row logs its violation and fix.

use fs_project::{
    Budgets, ConsequenceClass, Cooling, DecisionGate, EntityDecl, Envelope, FSIM_VERSION, Fan,
    GeometryArtifact, GeometryAssignment, HalfSpaceSide, InterfaceCardBinding, MaterialBinding,
    MeshSelector, Metadata, OutputRequest, PowerDissipation, ProjectSpec, RequirementDirection,
    RequirementSeverity, RequirementSource, RequirementSourceKind, SafetyFactorPolicy, Seeds,
    SolverSettings, ThermalLimit, UnitsDoctrine, Vent, Versions, canonical_hash, migrate_envelope,
    parse_json, parse_sexpr, parse_sexpr_lenient, print_json, print_sexpr,
    requirement_source_reviews,
};
use fs_qty::QtyAny;
use fs_scenario::EntityDeclaration;

fn kelvin(value: f64) -> QtyAny {
    QtyAny::new(value, fs_project::spec::dims::TEMPERATURE)
}

fn watts(value: f64) -> QtyAny {
    QtyAny::new(value, fs_project::spec::dims::POWER)
}

fn reference_assembly() -> Vec<EntityDecl> {
    vec![
        EntityDecl::Assembly {
            name: "enclosure-asm".to_string(),
            display: "Enclosure".to_string(),
            expect_id: None,
        },
        EntityDecl::Part {
            parent: "enclosure-asm".to_string(),
            name: "board".to_string(),
            display: "Main board".to_string(),
            expect_id: None,
        },
        EntityDecl::Region {
            parent: "board".to_string(),
            name: "cpu".to_string(),
            display: "CPU".to_string(),
            expect_id: None,
        },
        EntityDecl::Region {
            parent: "board".to_string(),
            name: "sink-base".to_string(),
            display: "Heat sink base".to_string(),
            expect_id: None,
        },
        EntityDecl::Interface {
            parent: "enclosure-asm".to_string(),
            name: "cpu-sink-tim".to_string(),
            display: "CPU to sink TIM".to_string(),
            from: "cpu".to_string(),
            to: "sink-base".to_string(),
            expect_id: None,
        },
    ]
}

/// The reference cooling project: a small forced-convection enclosure with
/// one board, one hot region, a heat sink interface, and a junction limit.
fn reference_project() -> ProjectSpec {
    ProjectSpec {
        metadata: Some(Metadata {
            name: "reference-cooling-v1".to_string(),
            created: "2026-07-22".to_string(),
            context_of_use: "design-stage screening of a fanned telecom enclosure".to_string(),
            intended_decision: "select the heat sink and fan operating point that hold the \
                                junction limit with margin"
                .to_string(),
            decision_gate: DecisionGate::DesignSelection,
            consequence: ConsequenceClass::Reliability,
        }),
        versions: Some(Versions {
            schema: FSIM_VERSION,
            constellation: "0".repeat(64),
            workspace: "e5c8061f4faed986b831b8978d0c8d1812e960fb".to_string(),
        }),
        seeds: Some(Seeds {
            master: 0x5EED_0001,
        }),
        budgets: Some(Budgets {
            solve_time: QtyAny::new(3600.0, fs_project::spec::dims::TIME),
            memory_bytes: 8 * 1024 * 1024 * 1024,
            accuracy_rel: 0.02,
        }),
        capabilities: Some(vec!["thermal.conduction-solve".to_string()]),
        units: Some(UnitsDoctrine {
            storage: "si-base".to_string(),
            display: "engineering".to_string(),
        }),
        geometry: Some(vec![GeometryArtifact {
            role: "enclosure".to_string(),
            format: "stl".to_string(),
            source_hash: 0x00ab_cdef_0123_4567,
            parser_version: "0.0.1".to_string(),
        }]),
        assignments: Some(vec![
            GeometryAssignment {
                artifact: "enclosure".to_string(),
                target: "cpu".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "CPU".to_string(),
                },
                allow_overlap: false,
            },
            GeometryAssignment {
                artifact: "enclosure".to_string(),
                target: "sink-base".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "SINK_BASE".to_string(),
                },
                allow_overlap: false,
            },
            GeometryAssignment {
                artifact: "enclosure".to_string(),
                target: "cpu-sink-tim".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "CPU_SINK_TIM".to_string(),
                },
                allow_overlap: false,
            },
        ]),
        assembly: Some(reference_assembly()),
        materials: Some(vec![MaterialBinding {
            region: "board".to_string(),
            card: "ab".repeat(32),
            claim: None,
            state: "fr4/nominal".to_string(),
            temp_lo: kelvin(233.15),
            temp_hi: kelvin(398.15),
            source: "matdb".to_string(),
        }]),
        interface_cards: Some(vec![InterfaceCardBinding {
            interface: "cpu-sink-tim".to_string(),
            card: "cd".repeat(32),
            claim: None,
            source: "matdb".to_string(),
        }]),
        power: Some(vec![PowerDissipation {
            region: "cpu".to_string(),
            watts: watts(35.0),
            duty: 1.0,
        }]),
        cooling: Some(Cooling {
            fans: vec![Fan {
                name: "intake-1".to_string(),
                flow: QtyAny::new(0.012, fs_project::spec::dims::VOLUMETRIC_FLOW),
                static_pressure: QtyAny::new(45.0, fs_project::spec::dims::PRESSURE),
            }],
            vents: vec![Vent {
                region: "sink-base".to_string(),
                area: QtyAny::new(0.004, fs_project::spec::dims::AREA),
            }],
            leakage: watts(2.5),
        }),
        envelope: Some(Envelope {
            ambient_lo: kelvin(273.15),
            ambient_hi: kelvin(318.15),
            pressure: QtyAny::new(101_325.0, fs_project::spec::dims::PRESSURE),
        }),
        requirements: Some(vec![ThermalLimit {
            qoi: "t-junction-max".to_string(),
            class: "junction".to_string(),
            region: "cpu".to_string(),
            direction: RequirementDirection::AtMost,
            limit: kelvin(378.15),
            margin: kelvin(10.0),
            source: RequirementSource {
                kind: RequirementSourceKind::Datasheet,
                document: "cpu-thermal-specification".to_string(),
                version: "rev-7".to_string(),
                locator: "table-5:tj-max".to_string(),
            },
            safety_factor: SafetyFactorPolicy {
                factor: 1.1,
                source: RequirementSource {
                    kind: RequirementSourceKind::InternalPolicy,
                    document: "thermal-derating-policy".to_string(),
                    version: "2026.1".to_string(),
                    locator: "section-4.2".to_string(),
                },
            },
            severity: RequirementSeverity::ReliabilityDerating,
        }]),
        solver: Some(SolverSettings {
            fidelity: "auto".to_string(),
            tolerance_rel: 1e-6,
        }),
        outputs: Some(vec![OutputRequest {
            name: "t-junction-max".to_string(),
            kind: "scalar".to_string(),
        }]),
    }
}

#[test]
fn requirement_authority_and_context_gate_are_mandatory_decision_inputs() {
    let reference = reference_project();
    let metadata = reference.metadata.as_ref().expect("metadata");
    assert!(!metadata.permits_indeterminate());

    let mut scoping = reference.clone();
    let metadata = scoping.metadata.as_mut().expect("metadata");
    metadata.decision_gate = DecisionGate::ScopingEstimate;
    metadata.consequence = ConsequenceClass::Advisory;
    assert!(metadata.permits_indeterminate());
    metadata.consequence = ConsequenceClass::SafetyCritical;
    assert!(
        !metadata.permits_indeterminate(),
        "safety consequence closes the exploratory escape hatch"
    );

    let mut sourceless = reference.clone();
    sourceless.requirements.as_mut().expect("requirements")[0]
        .source
        .document
        .clear();
    assert!(
        sourceless
            .validate()
            .iter()
            .any(|finding| finding.code == "project-requirement-source-invalid")
    );

    let mut invalid_factor = reference;
    invalid_factor.requirements.as_mut().expect("requirements")[0]
        .safety_factor
        .factor = 0.99;
    assert!(
        invalid_factor
            .validate()
            .iter()
            .any(|finding| finding.code == "project-requirement-safety-factor")
    );
}

#[test]
fn source_version_bumps_flag_only_the_matching_requirement_authorities() {
    let previous = reference_project().requirements.expect("requirements");
    let mut current = previous.clone();
    current[0].source.version = "rev-8".to_string();
    current[0].safety_factor.source.version = "2026.2".to_string();

    let reviews = requirement_source_reviews(&previous, &current);
    assert_eq!(reviews.len(), 2);
    assert_eq!(reviews[0].role, "requirement");
    assert_eq!(reviews[0].previous_version, "rev-7");
    assert_eq!(reviews[0].current_version, "rev-8");
    assert_eq!(reviews[1].role, "safety-factor");
    assert_eq!(reviews[1].previous_version, "2026.1");
    assert_eq!(reviews[1].current_version, "2026.2");

    let mut replacement = current;
    replacement[0].source.document = "replacement-datasheet".to_string();
    assert_eq!(
        requirement_source_reviews(&previous, &replacement),
        reviews[1..],
        "an entirely different authority is a project diff, not a version bump"
    );
}

#[test]
fn the_reference_cooling_project_is_admissible_and_hash_stable() {
    let spec = reference_project();
    let rendered = print_sexpr(&spec).expect("reference renders");
    let decoded = parse_sexpr(&rendered).expect("canonical bytes parse strictly");
    assert_eq!(decoded.spec, spec);
    assert!(
        decoded.findings().is_empty(),
        "reference project must be admissible: {:?}",
        decoded.findings()
    );
    assert!(decoded.defaults.is_empty());
    assert!(decoded.canonicalization.is_none());
    // Deterministic: two renders, one hash.
    let again = print_sexpr(&spec).expect("still renders");
    assert_eq!(rendered, again, "canonical bytes must be deterministic");
    assert_eq!(decoded.hash(), canonical_hash(rendered.as_bytes()));
}

#[test]
fn canonical_bytes_are_stable_across_spellings() {
    let spec = reference_project();
    let sexpr = print_sexpr(&spec).expect("s-expression renders");
    let json = print_json(&spec).expect("json renders");
    assert_ne!(sexpr, json, "the two spellings are distinct byte strings");

    let from_sexpr = parse_sexpr(&sexpr).expect("sexpr parses");
    let from_json = parse_json(&json).expect("json parses");
    assert_eq!(from_sexpr.spec, from_json.spec, "one semantic model");
    assert_eq!(
        from_sexpr.canonical, from_json.canonical,
        "both spellings canonicalize to the same bytes"
    );
    assert_eq!(from_sexpr.hash(), from_json.hash(), "one canonical hash");
}

#[test]
fn parse_render_parse_is_idempotent_across_project_variants() {
    let mut variants = vec![reference_project()];
    let mut no_fans = reference_project();
    if let Some(cooling) = &mut no_fans.cooling {
        cooling.fans.clear();
        cooling.vents.clear();
    }
    variants.push(no_fans);
    let mut pinned = reference_project();
    if let Some(assembly) = &mut pinned.assembly {
        let expected = EntityDeclaration::assembly("enclosure-asm")
            .with_display_name("Enclosure")
            .identity()
            .token();
        assembly[0] = EntityDecl::Assembly {
            name: "enclosure-asm".to_string(),
            display: "Enclosure".to_string(),
            expect_id: Some(expected),
        };
    }
    variants.push(pinned);

    for (index, spec) in variants.iter().enumerate() {
        let first = print_sexpr(spec).expect("renders");
        let decoded = parse_sexpr(&first).expect("parses");
        let second = print_sexpr(&decoded.spec).expect("re-renders");
        assert_eq!(first, second, "variant {index}: parse-render-parse drifted");
        assert!(
            decoded.findings().is_empty(),
            "variant {index}: {:?}",
            decoded.findings()
        );
    }
}

#[test]
fn every_assignment_selector_round_trips_through_both_spellings() {
    let selectors = [
        MeshSelector::NamedGroup {
            name: "CPU".to_string(),
        },
        MeshSelector::HalfSpace {
            normal: [0.0, 0.0, 1.0],
            offset: 1.0,
            side: HalfSpaceSide::AtLeast,
            tolerance: 1e-9,
        },
        MeshSelector::Box {
            min: [-1.0, -2.0, -3.0],
            max: [1.0, 2.0, 3.0],
            tolerance: 1e-8,
        },
        MeshSelector::Cylinder {
            origin: [0.0, 0.0, 0.0],
            axis: [0.0, 0.0, 1.0],
            radius: 0.5,
            axial_min: -1.0,
            axial_max: 1.0,
            tolerance: 1e-8,
        },
        MeshSelector::NearestDatum {
            point: [1.0, 2.0, 3.0],
            max_distance: 0.25,
            tolerance: 1e-9,
        },
        MeshSelector::ExplicitFaceSet {
            faces: vec![0, 2, u32::MAX],
            fragility_acknowledged: true,
        },
    ];

    for selector in selectors {
        let mut spec = reference_project();
        spec.assignments.as_mut().expect("assignments")[0].selector = selector;
        let sexpr = print_sexpr(&spec).expect("selector renders");
        let from_sexpr = parse_sexpr(&sexpr).expect("selector parses");
        assert_eq!(from_sexpr.spec, spec);
        let json = print_json(&spec).expect("selector JSON renders");
        let from_json = parse_json(&json).expect("selector JSON parses");
        assert_eq!(from_json.spec, spec);
        assert_eq!(from_json.hash(), from_sexpr.hash());
        assert!(from_sexpr.findings().is_empty());
    }
}

#[test]
fn every_mandatory_section_omission_is_a_named_violation() {
    let cases: [OmissionCase; 17] = [
        ("project-metadata-missing", |s| s.metadata = None),
        ("project-versions-missing", |s| s.versions = None),
        ("project-seeds-missing", |s| s.seeds = None),
        ("project-budgets-missing", |s| s.budgets = None),
        ("project-capabilities-missing", |s| s.capabilities = None),
        ("project-units-missing", |s| s.units = None),
        ("project-geometry-missing", |s| s.geometry = None),
        ("project-assignments-missing", |s| s.assignments = None),
        ("project-assembly-missing", |s| s.assembly = None),
        ("project-materials-missing", |s| s.materials = None),
        ("project-interface-cards-missing", |s| {
            s.interface_cards = None;
        }),
        ("project-power-missing", |s| s.power = None),
        ("project-cooling-missing", |s| s.cooling = None),
        ("project-envelope-missing", |s| s.envelope = None),
        ("project-requirements-missing", |s| s.requirements = None),
        ("project-solver-missing", |s| s.solver = None),
        ("project-outputs-missing", |s| s.outputs = None),
    ];
    for (code, mutate) in cases {
        assert_omission(code, mutate);
    }
}

/// One omission case: expected code plus the section-removing mutation.
type OmissionCase = (&'static str, fn(&mut ProjectSpec));

fn assert_omission(code: &str, mutate: fn(&mut ProjectSpec)) {
    let mut spec = reference_project();
    mutate(&mut spec);
    let findings = spec.validate();
    assert!(
        findings.iter().any(|v| v.code == code),
        "omission did not surface `{code}`: {findings:?}"
    );
    for finding in &findings {
        assert!(!finding.fix.trim().is_empty(), "{code}: fix hint is empty");
    }
}

#[test]
fn unknown_fields_are_refused_by_name() {
    let rendered = print_sexpr(&reference_project()).expect("renders");
    // An unknown top-level section.
    let with_section = rendered.replacen("(metadata", "(warp-drive :q 1) (metadata", 1);
    let decoded = parse_sexpr_lenient(&with_section).expect("still recognizable");
    assert!(
        decoded
            .recognition
            .iter()
            .any(|v| v.code == "project-unknown-field" && v.what.contains("warp-drive")),
        "unknown section not named: {:?}",
        decoded.recognition
    );
    // An unknown keyword inside a known section.
    let with_field = rendered.replacen(":fidelity", ":frobnicate 3 :fidelity", 1);
    let decoded = parse_sexpr_lenient(&with_field).expect("still recognizable");
    assert!(
        decoded
            .recognition
            .iter()
            .any(|v| v.code == "project-unknown-field" && v.what.contains("frobnicate")),
        "unknown keyword not named: {:?}",
        decoded.recognition
    );
}

#[test]
fn the_duty_default_is_receipted_never_silent() {
    let rendered = print_sexpr(&reference_project()).expect("renders");
    let without_duty = rendered.replacen(" :duty 1.0", "", 1);
    assert_ne!(rendered, without_duty, "fixture edit must apply");

    // Lenient: applied default carries a receipt, and the re-emission
    // carries a canonicalization receipt binding both byte strings.
    let decoded = parse_sexpr_lenient(&without_duty).expect("lenient parse");
    assert_eq!(decoded.defaults.len(), 1);
    let receipt = &decoded.defaults[0];
    assert_eq!(receipt.field, "power.dissipation[cpu].duty");
    assert_eq!(receipt.value, "1.0");
    assert!(!receipt.rationale.trim().is_empty());
    let canonicalization = decoded
        .canonicalization
        .expect("re-emission must be receipted");
    assert!(canonicalization.verifies(without_duty.as_bytes(), decoded.canonical.as_bytes()));
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());

    // Strict: the same bytes refuse.
    let refusal = parse_sexpr(&without_duty).expect_err("strict parse must refuse");
    assert!(
        refusal.code == "fsim-non-canonical" || refusal.code == "fsim-default-in-strict-mode",
        "unexpected refusal {refusal:?}"
    );
}

#[test]
fn noncanonical_whitespace_is_refused_strictly_and_receipted_leniently() {
    let rendered = print_sexpr(&reference_project()).expect("renders");
    let padded = format!("{rendered} ");
    let refusal = parse_sexpr(&padded).expect_err("strict refuses padding");
    assert_eq!(refusal.code, "fsim-non-canonical");

    let decoded = parse_sexpr_lenient(&padded).expect("lenient accepts with receipt");
    let receipt = decoded.canonicalization.expect("receipt owed");
    assert!(receipt.verifies(padded.as_bytes(), decoded.canonical.as_bytes()));
    assert_eq!(decoded.canonical, rendered);
}

#[test]
fn the_version_bump_machinery_is_proven_with_the_synthetic_migration() {
    let rendered = print_sexpr(&reference_project()).expect("renders");
    let v0 = rendered.replacen("(fsim-project :version 1", "(fsim-project :version 0", 1);

    // The reader refuses the old envelope outright.
    let refusal = parse_sexpr(&v0).expect_err("v0 must not parse directly");
    assert_eq!(refusal.code, "fsim-unsupported-version");

    // The explicit migration path re-emits canonical current-version bytes
    // and binds both byte strings into a verifying receipt.
    let migrated = migrate_envelope(&v0, 0).expect("registered rule migrates");
    assert_eq!(migrated.decoded.canonical, rendered);
    assert!(
        migrated
            .receipt
            .verifies(v0.as_bytes(), migrated.decoded.canonical.as_bytes())
    );
    assert_eq!(migrated.receipt.source_version, 0);
    assert_eq!(migrated.receipt.target_version, FSIM_VERSION);

    // Round trip: the migrated document parses strictly to the same spec.
    let reparsed = parse_sexpr(&migrated.decoded.canonical).expect("migrated bytes parse");
    assert_eq!(reparsed.spec, reference_project());

    // Migration refuses what it does not cover.
    assert_eq!(
        migrate_envelope(&rendered, FSIM_VERSION)
            .expect_err("no-op refused")
            .code,
        "fsim-migration-not-needed"
    );
    assert_eq!(
        migrate_envelope(&v0, 7)
            .expect_err("unknown version refused")
            .code,
        "fsim-migration-unknown-version"
    );
}

#[test]
fn entity_identity_pins_catch_drift() {
    // A correct pin is silent.
    let expected = EntityDeclaration::assembly("enclosure-asm")
        .with_display_name("Enclosure")
        .identity()
        .token();
    let mut pinned = reference_project();
    pinned.assembly.as_mut().expect("assembly")[0] = EntityDecl::Assembly {
        name: "enclosure-asm".to_string(),
        display: "Enclosure".to_string(),
        expect_id: Some(expected),
    };
    assert!(pinned.validate().is_empty(), "{:?}", pinned.validate());

    // A stale pin is a named violation, not a silent rebind.
    let mut drifted = reference_project();
    drifted.assembly.as_mut().expect("assembly")[0] = EntityDecl::Assembly {
        name: "enclosure-asm".to_string(),
        display: "Enclosure".to_string(),
        expect_id: Some("assembly:0000000000000000".to_string()),
    };
    let findings = drifted.validate();
    assert!(
        findings
            .iter()
            .any(|v| v.code == "project-entity-id-mismatch"),
        "{findings:?}"
    );
}

/// One corpus row: label, deliberately broken project, expected code.
type CorpusRow = (&'static str, ProjectSpec, &'static str);

/// Each corpus row is one deliberately broken project; the log doubles as
/// documentation of the error-message quality bar.
fn broken_corpus() -> Vec<CorpusRow> {
    vec![
        (
            "empty-capabilities",
            {
                let mut s = reference_project();
                s.capabilities = Some(Vec::new());
                s
            },
            "project-capabilities-empty",
        ),
        (
            "empty-power-map",
            {
                let mut s = reference_project();
                s.power = Some(Vec::new());
                s
            },
            "project-power-empty",
        ),
        (
            "wrong-power-dims",
            {
                let mut s = reference_project();
                s.power = Some(vec![PowerDissipation {
                    region: "cpu".to_string(),
                    watts: kelvin(35.0),
                    duty: 1.0,
                }]);
                s
            },
            "project-power-dims",
        ),
        (
            "duty-out-of-range",
            {
                let mut s = reference_project();
                s.power.as_mut().expect("power")[0].duty = 1.5;
                s
            },
            "project-duty-range",
        ),
        (
            "inverted-ambient-range",
            {
                let mut s = reference_project();
                let envelope = s.envelope.as_mut().expect("envelope");
                envelope.ambient_lo = kelvin(350.0);
                envelope.ambient_hi = kelvin(300.0);
                s
            },
            "project-envelope-range",
        ),
    ]
}

fn broken_corpus_references() -> Vec<CorpusRow> {
    vec![
        (
            "orphan-region-reference",
            {
                let mut s = reference_project();
                s.power.as_mut().expect("power")[0].region = "ghost".to_string();
                s
            },
            "project-ref-unknown",
        ),
        (
            "parent-declared-after-child",
            {
                let mut s = reference_project();
                let assembly = s.assembly.as_mut().expect("assembly");
                assembly.swap(0, 1);
                s
            },
            "project-entity-parent-unknown",
        ),
        (
            "duplicate-entity-name",
            {
                let mut s = reference_project();
                let assembly = s.assembly.as_mut().expect("assembly");
                let duplicate = assembly[2].clone();
                assembly.push(duplicate);
                s
            },
            "project-entity-duplicate",
        ),
        (
            "interface-card-on-unknown-interface",
            {
                let mut s = reference_project();
                s.interface_cards.as_mut().expect("cards")[0].interface = "ghost-tim".to_string();
                s
            },
            "project-ref-unknown",
        ),
        (
            "malformed-card-hash",
            {
                let mut s = reference_project();
                s.materials.as_mut().expect("materials")[0].card = "not-a-hash".to_string();
                s
            },
            "project-material-card",
        ),
        (
            "wrong-schema-version",
            {
                let mut s = reference_project();
                s.versions.as_mut().expect("versions").schema = 99;
                s
            },
            "project-schema-version-mismatch",
        ),
        (
            "non-si-storage",
            {
                let mut s = reference_project();
                s.units.as_mut().expect("units").storage = "imperial".to_string();
                s
            },
            "project-units-storage",
        ),
        (
            "bad-output-kind",
            {
                let mut s = reference_project();
                s.outputs.as_mut().expect("outputs")[0].kind = "hologram".to_string();
                s
            },
            "project-output-kind",
        ),
        (
            "vacuous-tolerance",
            {
                let mut s = reference_project();
                s.solver.as_mut().expect("solver").tolerance_rel = 0.0;
                s
            },
            "project-solver-tolerance",
        ),
    ]
}

#[test]
fn the_broken_project_corpus_logs_every_violation_with_its_fix() {
    eprintln!("RESULT\tCASE\tCODE\tFIX");
    let mut corpus = broken_corpus();
    corpus.extend(broken_corpus_references());
    for (label, spec, expected_code) in &corpus {
        let findings = spec.validate();
        if let Some(violation) = findings.iter().find(|v| v.code == *expected_code) {
            eprintln!("PASS\t{label}\t{}\t{}", violation.code, violation.fix);
            assert!(
                !violation.fix.trim().is_empty() && !violation.what.trim().is_empty(),
                "{label}: violation must carry what+fix"
            );
        } else {
            eprintln!("FAIL\t{label}\t{expected_code}\t-");
            panic!("{label}: expected `{expected_code}`, got {findings:?}");
        }
    }
}

#[test]
fn claim_pins_travel_through_both_spellings_and_validate_as_hex() {
    // A pinned binding round-trips canonically in both spellings and
    // reaches the same canonical hash.
    let mut spec = reference_project();
    spec.materials.as_mut().expect("materials")[0].claim = Some("ef".repeat(32));
    spec.interface_cards.as_mut().expect("interface cards")[0].claim = Some("1234".repeat(16));

    let sexpr = print_sexpr(&spec).expect("pinned project renders");
    assert!(sexpr.contains(":claim"));
    let decoded = parse_sexpr(&sexpr).expect("pinned project parses strictly");
    assert_eq!(decoded.spec, spec);
    assert!(decoded.spec.validate().is_empty());

    let json = print_json(&spec).expect("json renders");
    let from_json = parse_json(&json).expect("json parses");
    assert_eq!(
        from_json.hash(),
        decoded.hash(),
        "one hash across spellings"
    );

    // An absent pin renders no `:claim` field at all (absence is the
    // canonical spelling of "no pin"), so pre-pin documents are untouched.
    let unpinned = print_sexpr(&reference_project()).expect("renders");
    assert!(!unpinned.contains(":claim"));

    // A malformed pin is a named structural violation on each binding kind.
    let mut bad = reference_project();
    bad.materials.as_mut().expect("materials")[0].claim = Some("zz".repeat(32));
    bad.interface_cards.as_mut().expect("interface cards")[0].claim = Some("abc".to_string());
    let violations = bad.validate();
    assert!(
        violations
            .iter()
            .any(|v| v.code == "project-material-claim")
    );
    assert!(
        violations
            .iter()
            .any(|v| v.code == "project-interface-claim")
    );
}
