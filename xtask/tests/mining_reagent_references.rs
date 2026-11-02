//! Source-to-pack tests for pure anhydrous sodium carbonate.
//! Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C497198&Mask=2
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    ClaimSelection, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const MOLAR_MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

const PACK_ID: &str = "sodium-carbonate-solid-nist-shomate";
const MASS_KG_PER_MOL: f64 = 0.1059884;
// Independent transcription of the first solid Shomate column (298–723 K).
const SHOMATE: [f64; 5] = [175.2010, -348.0580, 743.0720, -305.5510, -1.634221];

fn workspace() -> PathBuf {
    std::env::var_os("CARGO_WORKSPACE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .into()
        })
}

fn compile() -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!(
        "mining-reagent-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{PACK_ID}/manifest.tsv")))
        .arg("--out")
        .arg(&out)
        .current_dir(workspace())
        .env("CARGO_WORKSPACE_DIR", workspace())
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        PACK_ID,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), PACK_ID);
    assert_eq!(pack.claims().claims_ordered().count(), 4);
    assert_eq!(pack.claims().observation_ids().count(), 4);
    pack
}

fn point(temperature: f64) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", semantic(QuantityKind::Pressure), 100000.0)
        .unwrap();
    for (axis, value) in [
        ("source-standard-state", 1.0),
        ("source-elemental-reference-janaf", 1.0),
        ("source-phase-solid", 1.0),
        ("source-pure-species", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn shomate_cp(temperature: f64) -> f64 {
    let [a, b, c, d, e] = SHOMATE;
    let t = temperature / 1000.0;
    (a + b * t + c * t * t + d * t.powi(3) + e / (t * t)) / MASS_KG_PER_MOL
}

fn shomate_h(temperature: f64) -> f64 {
    let [a, b, c, d, e] = SHOMATE;
    let t = temperature / 1000.0;
    // Integration constant cancels in differences; not a formation enthalpy.
    1000.0 * (a * t + b * t * t / 2.0 + c * t.powi(3) / 3.0 + d * t.powi(4) / 4.0 - e / t)
        / MASS_KG_PER_MOL
}

#[test]
fn sodium_carbonate_compiles_and_preserves_source_molar_units() {
    let pack = compile();
    for (name, dims, expected) in [
        ("molar-mass", MOLAR_MASS, MASS_KG_PER_MOL),
        (
            "standard-molar-formation-enthalpy",
            MOLAR_ENERGY,
            -1130770.0,
        ),
        ("standard-molar-entropy", MOLAR_ENTROPY, 138.79),
    ] {
        let answer = pack
            .claims()
            .query_typed(
                &PropertyKey::new(name, dims),
                &point(298.15),
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        close(answer.evidence.value.value, expected, 2e-14);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        assert!(answer.receipt.source_hashes.len() >= 3);
        pack.claims().verify_receipt(&answer.receipt).unwrap();
    }
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    let line = pack
        .claims()
        .query_line_typed(&cp, "temperature", &point(298.15), &point(700.0), SINGLE)
        .unwrap();
    assert_eq!(line.knots.len(), 42);
    for (t, answer) in &line.knots {
        close(answer.evidence.value.value, shomate_cp(*t), 2e-12);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        pack.claims().verify_receipt(&answer.receipt).unwrap();
    }
    // Independent rounded source table entries (J/mol/K), not mass-specific
    // values or additional uncertainty claims. Allow half the printed 0.1 unit.
    // https://webbook.nist.gov/cgi/cbook.cgi?ID=C497198&Mask=2&Table=on&Type=JANAFS
    for (t, source_cp_molar) in [(300.0, 111.3), (400.0, 125.1), (700.0, 187.5)] {
        let answer = pack
            .claims()
            .query_typed(&cp, &point(t), SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert!((answer.evidence.value.value * MASS_KG_PER_MOL - source_cp_molar).abs() <= 0.05);
    }
}

#[test]
fn sodium_carbonate_refuses_unsupported_states_and_reference_temperatures() {
    let pack = compile();
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    let bad_points = [
        QueryPoint::new().with("temperature", 300.0).unwrap(),
        point(298.0),
        point(701.0),
        point(300.0).with("source-pure-species", 0.0).unwrap(),
        point(300.0).with("source-phase-solid", 0.0).unwrap(),
        point(300.0)
            .with("source-elemental-reference-janaf", 0.0)
            .unwrap(),
        point(300.0).with("source-specimen-qualified", 1.0).unwrap(),
        point(300.0)
            .with_quantity("pressure", semantic(QuantityKind::Pressure), 101325.0)
            .unwrap(),
    ];
    for p in bad_points {
        assert!(
            pack.claims()
                .query_typed(&cp, &p, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
    }
    for (name, dims) in [
        ("standard-molar-formation-enthalpy", MOLAR_ENERGY),
        ("standard-molar-entropy", MOLAR_ENTROPY),
        ("density", Dims([-3, 1, 0, 0, 0, 0])),
    ] {
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &point(300.0),
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn sodium_carbonate_sensible_heat_matches_source_enthalpy_differences() {
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    let pack = compile();
    let line = pack
        .claims()
        .query_line_typed(&cp, "temperature", &point(300.0), &point(400.0), SINGLE)
        .unwrap();
    // Authored fixed-mass 2 kg sensible-heating example at standard pressure,
    // without reaction or phase change. Integrate the compiled linear Cp curve.
    let heat_j: f64 = line
        .knots
        .windows(2)
        .map(|w| (w[1].0 - w[0].0) * (w[0].1.evidence.value.value + w[1].1.evidence.value.value))
        .sum();
    let expected = 2.0 * (shomate_h(400.0) - shomate_h(300.0));
    close(heat_j, expected, 5e-5);
    assert!(heat_j > 0.0);
    // Independently printed source enthalpy increments are 0.20 and 12.01
    // kJ/mol at 300 and 400 K. Their difference has at most 0.01 kJ/mol
    // rounding ambiguity; the extra 0.001 allows the authored interpolation.
    let heat_kj_per_mol = heat_j * MASS_KG_PER_MOL / 2.0 / 1000.0;
    assert!((heat_kj_per_mol - (12.01 - 0.20)).abs() < 0.011);
    // Finite regression stencil, not a continuous interpolation or measurement
    // error certificate. Keep the shipped uncertainty explicitly unstated.
    for w in line.knots.windows(2) {
        let t = 0.5 * (w[0].0 + w[1].0);
        let answer = pack
            .claims()
            .query_typed(&cp, &point(t), SelectionPolicy::SingleClaimOnly)
            .unwrap();
        close(answer.evidence.value.value, shomate_cp(t), 5e-5);
    }
}
