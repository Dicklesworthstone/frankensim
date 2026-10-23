//! Piano string strength packs through the real offline compiler and the typed evaluator (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.369`, U2 follow-up): the music-wire minimum tensile
//! strength of QQ-W-470B Table II is answered per wire diameter exactly as declared: a tabulated row
//! exactly, a diameter between rows by the declared linear rule, and a diameter outside the transcribed
//! rows not at all. This establishes source transport, not experimental validation.

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

/// Compile one committed seed manifest as a material-state pack, exactly as the piano tools do.
fn compile(manifest: &str, state: [&str; 4]) -> NormalizedMaterialCardPack {
    let dir = std::env::temp_dir().join(format!(
        "frankensim-piano-string-packs-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).expect("fixture directory");
    let out = dir.join(format!("{}.fsmatpk", state.join("-")));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .arg("matdb-pack")
        .arg("--manifest")
        .arg(workspace().join(manifest))
        .arg("--out")
        .arg(&out)
        .arg("--material-state")
        .args(state)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .output()
        .expect("run xtask matdb-pack");
    assert!(
        run.status.success(),
        "{manifest}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    NormalizedMaterialCardPack::from_bytes(&std::fs::read(&out).expect("compiled pack"))
        .expect("decode compiled pack")
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
