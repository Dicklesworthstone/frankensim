//! M8-W1 G0 schema/refusal and G3 permutation/cross-reference checks.
//! Tolerance: exactly zero violations. No source is fetched by these tests.

use fs_vvreg::acoustic::{acoustic_absences, acoustic_cases};
use fs_vvreg::instrument_data::*;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

const BEAD: &str = "frankensim-irp-v5-instrument-realism-pexfp.190";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_owned()
}
fn committed_text() -> String {
    fs::read_to_string(root().join(INSTRUMENT_DATA_LEDGER_LOCATOR)).expect("read instrument ledger")
}
fn committed() -> InstrumentDataLedger {
    parse_instrument_data_ledger(&committed_text()).expect("committed ledger must parse")
}
fn fixture_fields() -> Vec<String> {
    [
        "test-source",
        "M8-W1",
        "material",
        "test source",
        "owner-measurement:test-source",
        "-",
        "FrankenSim-Original-Work",
        "permitted",
        "unstated",
        "authored",
        "fixture-only",
        "-",
        "2026-10-16",
        "authored schema fixture",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect()
}
fn document(rows: &[String]) -> String {
    format!(
        "{INSTRUMENT_DATA_LEDGER_HEADER}\n{}\n{}",
        LEDGER_FIELDS.join("\t"),
        rows.join("\n")
    )
}
fn parse_fields(fields: &[String]) -> Result<InstrumentDataLedger, LedgerError> {
    parse_instrument_data_ledger(&document(&[fields.join("\t")]))
}
fn json(value: &str) -> String {
    let mut out = String::from("\"");
    for c in value.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                use std::fmt::Write as _;
                write!(out, "\\u{:04x}", u32::from(c)).unwrap();
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
fn event(
    test: &str,
    source_ids: &[&str],
    path: Option<&Path>,
    rows: usize,
    refusals: usize,
    start: Instant,
    verdict: &str,
) {
    let artifacts = path.map_or_else(String::new, |path| {
        let bytes = fs::read(root().join(path)).expect("hash checked artifact");
        format!(
            "{{\"path\":{},\"blake3\":{}}}",
            json(&path.to_string_lossy()),
            json(&fs_blake3::hash_bytes(&bytes).to_hex())
        )
    });
    println!(
        "IRPDATA {{\"bead\":{},\"test\":{},\"source_ids\":[{}],\"artifacts\":[{}],\"rows_checked\":{},\"refusals_exercised\":{},\"ladder\":[],\"tolerance\":0,\"observed\":0,\"verdict\":{},\"elapsed_ms\":{}}}",
        json(BEAD),
        json(test),
        source_ids
            .iter()
            .map(|s| json(s))
            .collect::<Vec<_>>()
            .join(","),
        artifacts,
        rows,
        refusals,
        json(verdict),
        start.elapsed().as_millis()
    );
}
fn schema_event(test: &str, refusals: usize, start: Instant) {
    event(test, &["test-source"], None, 1, refusals, start, "pass");
}

#[test]
fn ledger_header_and_column_count_are_exact() {
    let start = Instant::now();
    let fields = fixture_fields();
    assert_eq!(parse_fields(&fields).unwrap().rows().len(), 1);
    for found in [13, 15] {
        let mut wrong = fields.clone();
        if found == 13 {
            wrong.pop();
        } else {
            wrong.push("extra".to_owned());
        }
        let error = parse_fields(&wrong).unwrap_err();
        assert_eq!(error.line, 3);
        assert_eq!(
            error.kind,
            LedgerErrorKind::ColumnCount {
                expected: 14,
                found
            }
        );
    }
    for input in [
        String::new(),
        document(&[]).replacen(INSTRUMENT_DATA_LEDGER_HEADER, "wrong", 1),
        document(&[]).replacen("source_id", "source", 1),
    ] {
        assert_eq!(
            parse_instrument_data_ledger(&input).unwrap_err().kind,
            LedgerErrorKind::Header
        );
    }
    assert!(
        parse_instrument_data_ledger(&document(&[]))
            .unwrap()
            .rows()
            .is_empty()
    );
    schema_event("ledger_header_and_column_count_are_exact", 5, start);
}

#[test]
fn every_required_field_refuses_when_blank() {
    let start = Instant::now();
    for (index, name) in LEDGER_FIELDS.iter().enumerate() {
        for blank in ["", "   "] {
            let mut fields = fixture_fields();
            fields[index] = blank.to_owned();
            assert_eq!(
                parse_fields(&fields).unwrap_err(),
                LedgerError {
                    line: 3,
                    kind: LedgerErrorKind::MissingField(name)
                }
            );
        }
    }
    schema_event("every_required_field_refuses_when_blank", 28, start);
}

#[test]
fn unknown_enumerated_values_refuse() {
    let start = Instant::now();
    for index in [2, 7, 8, 9, 10, 11] {
        let mut fields = fixture_fields();
        fields[index] = "not-a-known-value".to_owned();
        assert!(
            matches!(parse_fields(&fields).unwrap_err().kind, LedgerErrorKind::UnknownValue { field, .. } if field == LEDGER_FIELDS[index])
        );
    }
    let mut fields = fixture_fields();
    fields[11] = "foreign:pack".to_owned();
    assert!(matches!(
        parse_fields(&fields).unwrap_err().kind,
        LedgerErrorKind::UnknownValue {
            field: "artifacts",
            ..
        }
    ));
    schema_event("unknown_enumerated_values_refuse", 7, start);
}

#[test]
fn duplicate_source_id_refuses() {
    let start = Instant::now();
    let row = fixture_fields().join("\t");
    assert_eq!(
        parse_instrument_data_ledger(&document(&[row.clone(), row])).unwrap_err(),
        LedgerError {
            line: 4,
            kind: LedgerErrorKind::DuplicateSourceId("test-source".to_owned())
        }
    );
    schema_event("duplicate_source_id_refuses", 1, start);
}

#[test]
fn field_byte_cap_and_cap_plus_one() {
    let start = Instant::now();
    let mut fields = fixture_fields();
    // Multi-byte input proves that the cap counts bytes, not characters.
    fields[13] = "é".repeat(LEDGER_FIELD_BYTE_CAP / 2);
    assert_eq!(fields[13].len(), 4096);
    assert!(parse_fields(&fields).is_ok());
    fields[13].push('x');
    assert_eq!(
        parse_fields(&fields).unwrap_err().kind,
        LedgerErrorKind::FieldTooLong {
            field: "notes",
            cap: 4096
        }
    );
    schema_event("field_byte_cap_and_cap_plus_one", 1, start);
}

#[test]
fn decision_consistency_rules_refuse_each_violation() {
    let start = Instant::now();
    let mut cases = vec![
        vec![(10, "admit"), (7, "permission-required")],
        vec![(10, "admit"), (6, "unknown")],
        vec![(10, "admit-facts-only"), (7, "refused")],
        vec![(9, "figure-digitized"), (7, "facts-only"), (8, "stated")],
        vec![(9, "figure-digitized"), (8, "unstated")],
        vec![(9, "table-transcribed"), (13, "-")],
    ];
    for decision in ["refuse", "pending-hunt", "pending-permission"] {
        cases.push(vec![(10, decision), (11, "matdb:test-pack")]);
    }
    for kind in ["geometry", "anatomy", "passive-response"] {
        cases.push(vec![(2, kind)]);
    }
    for edits in &cases {
        let mut fields = fixture_fields();
        for &(index, value) in edits {
            fields[index] = value.to_owned();
        }
        assert!(
            matches!(
                parse_fields(&fields).unwrap_err().kind,
                LedgerErrorKind::Inconsistent(_)
            ),
            "edits={edits:?}"
        );
    }
    // The corresponding permitted boundaries do not refuse.
    for decision in ["admit", "admit-facts-only", "fixture-only"] {
        let mut fields = fixture_fields();
        fields[10] = decision.to_owned();
        assert!(parse_fields(&fields).is_ok());
    }
    for uncertainty in ["stated", "bounded"] {
        let mut fields = fixture_fields();
        fields[9] = "figure-digitized".to_owned();
        fields[8] = uncertainty.to_owned();
        assert!(parse_fields(&fields).is_ok());
    }
    for decision in ["pending-hunt", "pending-permission"] {
        let mut fields = fixture_fields();
        fields[6] = "unknown".to_owned();
        fields[10] = decision.to_owned();
        fields[11] = "vvreg-absence:acoustic-absent-felt-coupon".to_owned();
        assert!(parse_fields(&fields).is_ok());
    }
    schema_event(
        "decision_consistency_rules_refuse_each_violation",
        cases.len(),
        start,
    );
}

#[test]
fn identifiers_dates_and_artifact_paths_are_checked() {
    let start = Instant::now();
    let bad = [
        (0, "Bad_ID"),
        (5, "../specimen"),
        (6, "terms:"),
        (12, "2025-02-29"),
        (12, "2026-13-01"),
        (12, "2026-04-31"),
        (12, "2026-1-01"),
        (11, "table:/etc/passwd"),
        (11, "table:data/../escape"),
        (11, "matdb:../escape"),
        (11, "reference:"),
        (11, "table:data/file;"),
        (13, "hidden\u{0000}control"),
    ];
    for (index, value) in bad {
        let mut fields = fixture_fields();
        fields[index] = value.to_owned();
        assert!(
            matches!(parse_fields(&fields).unwrap_err().kind, LedgerErrorKind::UnknownValue { field, .. } if field == LEDGER_FIELDS[index])
        );
    }
    let mut fields = fixture_fields();
    fields[12] = "2024-02-29".to_owned();
    fields[11] = "matdb:test-pack;vvreg:test-case;vvreg-absence:test-absence;geometry:test-specimen;reference:test-reference;table:data/test.tsv".to_owned();
    let ledger = parse_fields(&fields).unwrap();
    assert_eq!(ledger.rows()[0].artifacts.len(), 6);
    assert!(ledger.by_source("missing").is_none());
    assert!(ledger.for_item("missing").is_empty());
    assert!(
        ledger
            .for_artifact(&ArtifactRef::Matdb("missing".to_owned()))
            .is_empty()
    );
    // All vocabulary tags round-trip, including currently unseeded methods.
    for digitization in [
        Digitization::None,
        Digitization::TextTranscribed,
        Digitization::TableTranscribed,
        Digitization::FigureDigitized,
        Digitization::OwnerMeasured,
        Digitization::ServiceMeasured,
        Digitization::Scan,
        Digitization::Photogrammetry,
        Digitization::Authored,
    ] {
        fields[8] = UncertaintyStatement::Bounded.tag().to_owned();
        fields[9] = digitization.tag().to_owned();
        assert_eq!(
            parse_fields(&fields).unwrap().rows()[0].digitization,
            digitization
        );
    }
    for kind in [
        DataKind::Material,
        DataKind::Geometry,
        DataKind::ReferenceAudio,
        DataKind::PlayerData,
        DataKind::MachinePlayerData,
        DataKind::PassiveResponse,
        DataKind::Anatomy,
        DataKind::Method,
    ] {
        fields[2] = kind.tag().to_owned();
        fields[5] = "test-specimen".to_owned();
        assert_eq!(parse_fields(&fields).unwrap().rows()[0].kind, kind);
    }
    schema_event(
        "identifiers_dates_and_artifact_paths_are_checked",
        bad.len(),
        start,
    );
}

#[test]
fn row_permutation_is_lookup_invariant() {
    let start = Instant::now();
    let text = committed_text();
    let original = parse_instrument_data_ledger(&text).unwrap();
    let mut rows: Vec<_> = text.lines().skip(2).map(str::to_owned).collect();
    rows.reverse();
    let permuted = parse_instrument_data_ledger(&document(&rows)).unwrap();
    assert_eq!(original, permuted);
    for row in original.rows() {
        assert_eq!(
            original.by_source(&row.source_id),
            permuted.by_source(&row.source_id)
        );
        assert_eq!(original.for_item(&row.item), permuted.for_item(&row.item));
        for artifact in &row.artifacts {
            assert_eq!(
                original.for_artifact(artifact),
                permuted.for_artifact(artifact)
            );
        }
    }
    event(
        "row_permutation_is_lookup_invariant",
        &[],
        Some(Path::new(INSTRUMENT_DATA_LEDGER_LOCATOR)),
        original.rows().len(),
        0,
        start,
        "pass",
    );
}

fn manifest_field<'a>(text: &'a str, key: &str) -> Result<&'a str, String> {
    let values: Vec<_> = text
        .lines()
        .filter_map(|line| {
            let mut columns = line.split('\t');
            (columns.next() == Some(key))
                .then(|| columns.next())
                .flatten()
        })
        .collect();
    if values.len() == 1 {
        Ok(values[0])
    } else {
        Err(format!(
            "manifest field {key}: expected exactly one, found {}",
            values.len()
        ))
    }
}

fn validate_references(ledger: &InstrumentDataLedger, emit: bool) -> Result<usize, String> {
    let start = Instant::now();
    let mut checked = 0;
    for row in ledger.rows() {
        for artifact in &row.artifacts {
            let mut verdict = "pass";
            let path = match artifact {
                ArtifactRef::Matdb(slug) => {
                    let path = PathBuf::from(format!("data/matdb/seed-v1/{slug}/manifest.tsv"));
                    let text = fs::read_to_string(root().join(&path))
                        .map_err(|e| format!("{}: {}: {e}", row.source_id, path.display()))?;
                    for (key, expected) in [
                        ("license", row.licence.as_str()),
                        ("redistribution", row.redistribution.tag()),
                    ] {
                        let actual = manifest_field(&text, key)
                            .map_err(|e| format!("{}: {e}", row.source_id))?;
                        if actual != expected {
                            return Err(format!(
                                "{}: {key} mismatch for {slug}: ledger={expected}, manifest={actual}",
                                row.source_id
                            ));
                        }
                    }
                    Some(path)
                }
                ArtifactRef::Vvreg(id) => {
                    if !acoustic_cases().iter().any(|case| case.id == id.as_str()) {
                        return Err(format!("{}: unresolved vvreg:{id}", row.source_id));
                    }
                    Some(PathBuf::from("crates/fs-vvreg/src/acoustic.rs"))
                }
                ArtifactRef::VvregAbsence(id) => {
                    if !acoustic_absences()
                        .iter()
                        .any(|case| case.id == id.as_str())
                    {
                        return Err(format!("{}: unresolved vvreg-absence:{id}", row.source_id));
                    }
                    Some(PathBuf::from("crates/fs-vvreg/src/acoustic.rs"))
                }
                ArtifactRef::Geometry(id) => Some(PathBuf::from(format!(
                    "data/instruments/geometry/{id}/specimen.tsv"
                ))),
                ArtifactRef::Table(path) => Some(PathBuf::from(path)),
                ArtifactRef::Reference(_) => {
                    verdict = "deferred-to-M8-W6a";
                    None
                }
            };
            if let Some(path) = &path {
                if !root().join(path).is_file() {
                    return Err(format!("{}: missing {}", row.source_id, path.display()));
                }
            }
            if emit {
                event(
                    "committed_ledger_parses_and_every_artifact_resolves",
                    &[&row.source_id],
                    path.as_deref(),
                    1,
                    0,
                    start,
                    verdict,
                );
            }
            checked += 1;
        }
    }
    Ok(checked)
}

#[test]
fn committed_ledger_parses_and_every_artifact_resolves() {
    let start = Instant::now();
    let ledger = committed();
    assert!(
        validate_references(&ledger, true).unwrap() >= 50,
        "seed checks must be non-vacuous"
    );
    event(
        "committed_ledger_parses_and_every_artifact_resolves",
        &[],
        Some(Path::new(INSTRUMENT_DATA_LEDGER_LOCATOR)),
        ledger.rows().len(),
        0,
        start,
        "pass",
    );
}

fn retained_path(ledger: &InstrumentDataLedger, path: &Path) -> bool {
    ledger
        .rows()
        .iter()
        .filter(|row| row.decision.retains_artifacts())
        .any(|row| {
            row.artifacts.iter().any(|artifact| match artifact {
                ArtifactRef::Table(table) => path == Path::new(table),
                ArtifactRef::Geometry(id) => {
                    path.starts_with(PathBuf::from(format!("data/instruments/geometry/{id}")))
                }
                _ => false,
            })
        })
}
fn files_under(dir: &Path, files: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).expect("read instrument directory") {
        let entry = entry.unwrap();
        let kind = entry.file_type().unwrap();
        assert!(
            !kind.is_symlink(),
            "instrument data symlink needs an explicit retention policy: {}",
            entry.path().display()
        );
        if kind.is_dir() {
            files_under(&entry.path(), files);
        } else {
            files.push(entry.path());
        }
    }
}

#[test]
fn committed_ledger_has_no_orphan_instrument_data_files() {
    let start = Instant::now();
    let ledger = committed();
    let mut paths = Vec::new();
    files_under(&root().join("data/instruments"), &mut paths);
    paths.sort();
    for absolute in paths {
        let base = root();
        let path = absolute.strip_prefix(&base).unwrap();
        if path == Path::new(INSTRUMENT_DATA_LEDGER_LOCATOR) {
            continue;
        }
        assert!(
            retained_path(&ledger, path),
            "orphan instrument file: {}",
            path.display()
        );
        event(
            "committed_ledger_has_no_orphan_instrument_data_files",
            &[],
            Some(path),
            1,
            0,
            start,
            "pass",
        );
    }
    assert!(!retained_path(
        &ledger,
        Path::new("data/instruments/not-registered.tsv")
    ));
    event(
        "committed_ledger_has_no_orphan_instrument_data_files",
        &[],
        Some(Path::new(INSTRUMENT_DATA_LEDGER_LOCATOR)),
        ledger.rows().len(),
        1,
        start,
        "pass",
    );
}

const SEED_PACKS: &[&str] = &[
    "cane-arundo-damping-scielo-mr20170795",
    "cane-arundo-reedblank-mdpi-ma18122759",
    "cane-arundo-reedheel-mdpi-ma13204566",
    "lip-perioral-invivo-mdpi-ma17153654",
    "vocalfold-felid-rheometry-plos-pone0027029",
    "vocalfold-human-indentation-jvoice-pmc12180296",
    "vocalfold-porcine-aspiration-mdpi-s21092923",
    "tissue-density-stated-plos-pcbi1004907",
    "twomass-if72-standard-plos-pone0187486",
    "spruce-sitka-loss-qiu2026",
    "spruce-norway-danihelova2022",
    "maple-sycamore-danihelova2022",
    "ash-white-fpl-gtr282",
    "baldcypress-fpl-gtr282",
    "balsa-fpl-gtr282",
    "basswood-american-fpl-gtr282",
    "birch-yellow-fpl-gtr282",
    "cedar-western-red-fpl-gtr282",
    "cherry-black-fpl-gtr282",
    "cottonwood-eastern-fpl-gtr282",
    "douglas-fir-coast-fpl-gtr282",
    "hemlock-western-fpl-gtr282",
    "larch-western-fpl-gtr282",
    "mahogany-african-fpl-gtr282",
    "mahogany-honduras-fpl-gtr282",
    "maple-red-fpl-gtr282",
    "maple-sugar-fpl-gtr282",
    "redwood-old-growth-fpl-gtr282",
    "rosewood-brazilian-fpl-gtr282",
    "rosewood-indian-fpl-gtr282",
    "spruce-engelmann-fpl-gtr282",
    "spruce-sitka-fpl-gtr282",
    "sweetgum-fpl-gtr282",
    "walnut-black-fpl-gtr282",
    "yellow-poplar-fpl-gtr282",
    "aluminum-2024-t4-damping-nasa-tn-d2893",
    "aluminum-2024-t3-nasa-tn-d6448",
    "brass-yellow-half-hard-damping-nasa-tn-d1467",
    "brass-cartridge-c26000-mil-hdbk-698a",
    "phosphor-bronze-c51000-nist-mono177",
    "phosphor-bronze-5a-mil-hdbk-698a",
    "music-wire-nbs-c447",
    "music-wire-spring-moduli-fuchs-1968",
    "bronze-gunmetal-damping-rsic508",
    "contact-cane-reed-on-mouthpiece-lay",
    "contact-string-on-fret-nickel",
    "contact-jawari-bone-bridge",
];

#[test]
fn committed_ledger_seed_covers_existing_music_sources() {
    let start = Instant::now();
    let ledger = committed();
    for slug in SEED_PACKS {
        assert!(
            !ledger
                .for_artifact(&ArtifactRef::Matdb((*slug).to_owned()))
                .is_empty(),
            "missing seed matdb:{slug}"
        );
    }
    for prefix in [
        "acoustic-ernoult-2021-",
        "acoustic-carcagno-2018-",
        "acoustic-assaneo-2011-",
    ] {
        let cases: Vec<_> = acoustic_cases()
            .iter()
            .filter(|case| case.id.starts_with(prefix))
            .collect();
        assert!(!cases.is_empty());
        for case in cases {
            assert!(
                !ledger
                    .for_artifact(&ArtifactRef::Vvreg(case.id.to_owned()))
                    .is_empty(),
                "missing {}",
                case.id
            );
        }
    }
    for id in [
        "acoustic-refused-openwind-curves",
        "acoustic-absent-felt-coupon",
        "acoustic-absent-trumpet-impedance",
        "acoustic-absent-measured-frf-trace",
    ] {
        assert!(
            !ledger
                .for_artifact(&ArtifactRef::VvregAbsence(id.to_owned()))
                .is_empty(),
            "missing {id}"
        );
    }
    assert!(
        !ledger
            .for_artifact(&ArtifactRef::Table(
                "data/felt/steinway-us5125310-strip-taper.tsv".to_owned()
            ))
            .is_empty()
    );
    for slug in [
        "contact-cane-reed-on-mouthpiece-lay",
        "contact-string-on-fret-nickel",
        "contact-jawari-bone-bridge",
        "twomass-if72-standard-plos-pone0187486",
    ] {
        assert!(
            ledger
                .for_artifact(&ArtifactRef::Matdb(slug.to_owned()))
                .iter()
                .all(|row| row.decision == Decision::FixtureOnly),
            "fixture promoted: {slug}"
        );
    }
    for id in [
        "strad3d-ct",
        "brandeis-messiah-ct",
        "usc-morphdb",
        "usc-75-speaker-mri",
        "dvtd",
        "goettingen-brass-mri",
        "fraunhofer-musices",
        "maugeais-dalmont-2024",
        "davies-pinho-english-2006",
        "dequand-2003",
        "fletcher-2005",
        "met-1995-541-8",
        "horniman-m32-1990",
        "swift-bosworth-2021",
        "uiowa-mis",
    ] {
        let row = ledger
            .by_source(id)
            .unwrap_or_else(|| panic!("missing acquisition lead {id}"));
        assert!(matches!(
            row.decision,
            Decision::PendingHunt | Decision::PendingPermission
        ));
        assert!(
            row.artifacts.is_empty(),
            "unacquired lead must not provide artifacts"
        );
    }
    let jawari = ledger.for_artifact(&ArtifactRef::Matdb("contact-jawari-bone-bridge".to_owned()));
    assert!(
        jawari[0]
            .notes
            .contains("calibrated to a fixture's output, BR-03")
    );
    event(
        "committed_ledger_seed_covers_existing_music_sources",
        &[],
        Some(Path::new(INSTRUMENT_DATA_LEDGER_LOCATOR)),
        ledger.rows().len(),
        0,
        start,
        "pass",
    );
}

#[test]
fn falsifier_licence_mismatch_is_caught() {
    let start = Instant::now();
    let text = committed_text();
    let mut rows: Vec<_> = text.lines().skip(2).map(str::to_owned).collect();
    let index = rows
        .iter()
        .position(|r| {
            r.split('\t')
                .nth(11)
                .unwrap()
                .contains("matdb:cane-arundo-damping-")
        })
        .expect("cane seed");
    let mut fields: Vec<_> = rows[index].split('\t').map(str::to_owned).collect();
    let source_id = fields[0].clone();
    fields[6] = "MIT".to_owned();
    rows[index] = fields.join("\t");
    let mutant = parse_instrument_data_ledger(&document(&rows)).unwrap();
    let error = validate_references(&mutant, false).unwrap_err();
    assert!(
        error.contains(&source_id) && error.contains("license mismatch"),
        "{error}"
    );
    fields[6] = "CC-BY-4.0".to_owned();
    fields[7] = "facts-only".to_owned();
    rows[index] = fields.join("\t");
    let mutant = parse_instrument_data_ledger(&document(&rows)).unwrap();
    let error = validate_references(&mutant, false).unwrap_err();
    assert!(
        error.contains(&source_id) && error.contains("redistribution mismatch"),
        "{error}"
    );
    event(
        "falsifier_licence_mismatch_is_caught",
        &[&source_id],
        Some(Path::new(INSTRUMENT_DATA_LEDGER_LOCATOR)),
        1,
        2,
        start,
        "pass",
    );
}
