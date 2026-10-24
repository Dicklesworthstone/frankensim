//! Piano string strength packs through the real offline compiler and the typed evaluator (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.369`, U2 follow-up): the music-wire minimum tensile
//! strength of QQ-W-470B Table II is answered per wire diameter exactly as declared: a tabulated row
//! exactly, a diameter between rows by the declared linear rule, and a diameter outside the transcribed
//! rows not at all. This establishes source transport, not experimental validation.
//!
//! The engine-piano string tests load compiled packs committed under
//! `apps/instrument-lab/engine-piano/tests/fixtures/packs` (that crate has its own workspace and cannot
//! run this compiler). Each fixture must be byte-for-byte this compiler's output for its committed seed,
//! so a seed edit without a refreshed fixture fails here. The grand piano board's loss packs
//! (`crates/fs-couple/examples/grand_piano/packs`, pexfp.373) are guarded the same way; they are plain
//! claim packs, compiled without a material state, as the board source loads them. So are the plain
//! packs the fs-couple literature tests resolve their materials from
//! (`crates/fs-couple/tests/fixtures/matdb`, pexfp.32): seeds, and labelled authored fixtures whose
//! sources live beside them.

#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;

use fs_matdb::{
    EvaluationDecision, MatDbError, NormalizedMaterialCardPack, QueryPoint, SelectionPolicy,
};
use fs_qty::{
    Length, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const QQ_W_470B: &str = "data/matdb/seed-v1/steel-music-wire-min-strength-qq-w-470b/manifest.tsv";
const MPA_PER_KSI: f64 = 6.894757293168;
const MM_PER_IN: f64 = 25.4;

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a workspace parent")
        .to_path_buf()
}

/// Compiled packs the engine-piano string tests load, by seed directory and material state.
const PIANO_FIXTURES: [(&str, [&str; 4]); 6] = [
    ("gold-pure-nbs-c447", ["gold", "solid", "nbs-c447", "0"]),
    (
        "steel-music-wire-string-heptonstall-ligo-p1400229",
        ["steel", "solid", "music-wire-untreated", "0"],
    ),
    (
        "steel-music-wire-min-strength-qq-w-470b",
        ["steel", "solid", "qq-w-470b-minimum", "0"],
    ),
    (
        "copper-pure-nbs-c447-string-wrap",
        ["copper", "solid", "nbs-c447", "0"],
    ),
    ("silver-pure-nbs-c447", ["silver", "solid", "nbs-c447", "0"]),
    (
        "silver-drawn-hard-uts-nbs-c101",
        ["silver", "solid", "drawn-hard-nbs-c101", "0"],
    ),
];

/// Loss packs the grand piano's material-damping board loads (plain claim packs).
const BOARD_LOSS_FIXTURES: [&str; 6] = [
    "spruce-sitka-loss-moisture-qiu2026",
    "spruce-sitka-loss-ono1996",
    "maple-sugar-loss-ono1996",
    "spruce-norway-loss-frequency-spycher2008",
    "spruce-norway-loss-frequency-moisture-merhar2025",
    "pine-loss-longitudinal-bremaud2012",
];

/// Plain packs the fs-couple literature tests load, by source manifest (pexfp.32): seed packs, then
/// the authored test fixtures that carry values no source pack does.
const COUPLE_TEST_FIXTURES: [&str; 13] = [
    "data/matdb/seed-v1/spruce-sitka-fpl-gtr282/manifest.tsv",
    "data/matdb/seed-v1/spruce-sitka-loss-qiu2026/manifest.tsv",
    "data/matdb/seed-v1/rosewood-brazilian-fpl-gtr282/manifest.tsv",
    "data/matdb/seed-v1/rosewood-indian-fpl-gtr282/manifest.tsv",
    "data/matdb/seed-v1/aluminum-2024-t3-nasa-tn-d6448/manifest.tsv",
    "data/matdb/seed-v1/music-wire-nbs-c447/manifest.tsv",
    "data/matdb/seed-v1/music-wire-spring-moduli-fuchs-1968/manifest.tsv",
    "data/matdb/seed-v1/cane-arundo-reedheel-mdpi-ma13204566/manifest.tsv",
    "data/matdb/seed-v1/cane-arundo-reedblank-mdpi-ma18122759/manifest.tsv",
    "data/matdb/seed-v1/cane-arundo-damping-scielo-mr20170795/manifest.tsv",
    "crates/fs-couple/tests/fixtures/matdb/rosewood-brazilian-back-authored-fixture/manifest.tsv",
    "crates/fs-couple/tests/fixtures/matdb/rosewood-indian-bar-authored-fixture/manifest.tsv",
    "crates/fs-couple/tests/fixtures/matdb/cane-reed-authored-fixture/manifest.tsv",
];

/// Compile one committed seed manifest as a material-state pack, exactly as the piano tools do.
fn compile(manifest: &str, state: [&str; 4]) -> NormalizedMaterialCardPack {
    let (_, bytes) = compile_bytes(manifest, Some(state), &state.join("-"));
    NormalizedMaterialCardPack::from_bytes(&bytes).expect("decode compiled pack")
}

/// The compiler's bytes for one seed (a material-state pack when `state` is
/// given, else a plain claim pack), written to `<temp>/<name>.fsmatpk`, and that path.
fn compile_bytes(manifest: &str, state: Option<[&str; 4]>, name: &str) -> (PathBuf, Vec<u8>) {
    let dir = std::env::temp_dir().join(format!(
        "frankensim-piano-string-packs-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).expect("fixture directory");
    let out = dir.join(format!("{name}.fsmatpk"));
    // The compiler refuses to overwrite; a reused process id may have left this temp file.
    let _ = std::fs::remove_file(&out);
    let mut command = Command::new(env!("CARGO_BIN_EXE_xtask"));
    command
        .arg("matdb-pack")
        .arg("--manifest")
        .arg(workspace().join(manifest))
        .arg("--out")
        .arg(&out);
    if let Some(state) = state {
        command.arg("--material-state").args(state);
    }
    let run = command
        .env("CARGO_WORKSPACE_DIR", workspace())
        .output()
        .expect("run xtask matdb-pack");
    assert!(
        run.status.success(),
        "{manifest}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(&out).expect("compiled pack");
    (out, bytes)
}

#[test]
fn engine_piano_pack_fixtures_are_the_compiled_committed_seeds() {
    let fixtures = workspace().join("apps/instrument-lab/engine-piano/tests/fixtures/packs");
    let mut stale = Vec::new();
    for (seed, state) in PIANO_FIXTURES {
        let manifest = format!("data/matdb/seed-v1/{seed}/manifest.tsv");
        let (fresh, bytes) = compile_bytes(&manifest, Some(state), seed);
        let committed = fixtures.join(format!("{seed}.fsmatpk"));
        if std::fs::read(&committed).ok().as_deref() != Some(bytes.as_slice()) {
            stale.push(format!(
                "{} (fresh: {})",
                committed.display(),
                fresh.display()
            ));
        }
    }
    assert!(
        stale.is_empty(),
        "stale or missing engine-piano pack fixtures; commit the fresh compiler output:\n{}",
        stale.join("\n")
    );
}

#[test]
fn grand_piano_loss_pack_fixtures_are_the_compiled_committed_seeds() {
    let fixtures = workspace().join("crates/fs-couple/examples/grand_piano/packs");
    let mut stale = Vec::new();
    for seed in BOARD_LOSS_FIXTURES {
        let manifest = format!("data/matdb/seed-v1/{seed}/manifest.tsv");
        let (fresh, bytes) = compile_bytes(&manifest, None, &format!("board-{seed}"));
        let committed = fixtures.join(format!("{seed}.fsmatpk"));
        if std::fs::read(&committed).ok().as_deref() != Some(bytes.as_slice()) {
            stale.push(format!(
                "{} (fresh: {})",
                committed.display(),
                fresh.display()
            ));
        }
    }
    assert!(
        stale.is_empty(),
        "stale or missing grand piano loss pack fixtures; commit the fresh compiler output:\n{}",
        stale.join("\n")
    );
}

#[test]
fn couple_test_pack_fixtures_are_the_compiled_committed_sources() {
    let fixtures = workspace().join("crates/fs-couple/tests/fixtures/matdb");
    let mut stale = Vec::new();
    for manifest in COUPLE_TEST_FIXTURES {
        let id = Path::new(manifest)
            .parent()
            .and_then(Path::file_name)
            .and_then(|name| name.to_str())
            .expect("a pack directory");
        let (fresh, bytes) = compile_bytes(manifest, None, &format!("couple-{id}"));
        let committed = fixtures.join(format!("{id}.fsmatpk"));
        if std::fs::read(&committed).ok().as_deref() != Some(bytes.as_slice()) {
            stale.push(format!(
                "{} (fresh: {})",
                committed.display(),
                fresh.display()
            ));
        }
    }
    assert!(
        stale.is_empty(),
        "stale or missing fs-couple test pack fixtures; commit the fresh compiler output:\n{}",
        stale.join("\n")
    );
}

/// Room temperature on the typed `T` axis and one wire diameter [m].
fn at(diameter_m: f64) -> QueryPoint {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    QueryPoint::new()
        .with_quantity("T", temperature, 293.15)
        .and_then(|point| {
            point.with_quantity(
                "wire_diameter",
                QuantitySpec::dimensional(Length::DIMS),
                diameter_m,
            )
        })
        .expect("query point")
}

/// A printed row in MPa as the pack stores it: to four decimals.
fn stored_mpa(ksi: f64) -> f64 {
    (ksi * MPA_PER_KSI * 1e4).round() / 1e4
}

#[test]
fn qq_w_470b_strength_is_answered_per_diameter_inside_its_rows_only() {
    let pack = compile(QQ_W_470B, ["steel", "solid", "qq-w-470b-minimum", "0"]);
    let claims = pack.card().claims();
    let query = |diameter_mm: f64| {
        claims.query(
            "minimum_tensile_strength",
            &at(diameter_mm * 1e-3),
            SelectionPolicy::SingleClaimOnly,
        )
    };

    // A tabulated row: 0.040 in, 315 thousand p.s.i.
    let row = query(0.040 * MM_PER_IN).expect("tabulated row");
    assert_eq!(row.evidence.value.value, stored_mpa(315.0) * 1e6);
    assert!(
        matches!(
            row.receipt.decision,
            EvaluationDecision::ExactTabulated { .. }
        ),
        "{:?}",
        row.receipt.decision
    );

    // The Model D's thickest plain core (key 41, 1.117 mm) lies between 0.042 in (313) and
    // 0.045 in (309): the declared linear rule, between exactly those two rows.
    let (x0, x1) = (0.042 * MM_PER_IN, 0.045 * MM_PER_IN);
    let (y0, y1) = (stored_mpa(313.0), stored_mpa(309.0));
    let diameter = 1.117;
    let expected = (y0 + (y1 - y0) * (diameter - x0) / (x1 - x0)) * 1e6;
    let inside = query(diameter).expect("between rows");
    assert!(
        (inside.evidence.value.value / expected - 1.0).abs() < 1e-12,
        "{} vs {expected}",
        inside.evidence.value.value
    );
    match inside.receipt.decision {
        EvaluationDecision::LinearInside { x_lo, x_hi } => {
            assert!((x_lo - x0 * 1e-3).abs() < 1e-15 && (x_hi - x1 * 1e-3).abs() < 1e-15);
        }
        other => panic!("expected a linear interpolation, got {other:?}"),
    }

    // Outside the transcribed rows (0.508-2.413 mm) nothing is claimed: refused, never extrapolated.
    for diameter in [0.4, 2.6] {
        let refused = query(diameter).expect_err("outside the transcribed rows");
        assert!(
            matches!(
                refused,
                MatDbError::NoClaimInDomain { .. } | MatDbError::OutsideKnotSpan { .. }
            ),
            "{diameter} mm: {refused:?}"
        );
    }
}
