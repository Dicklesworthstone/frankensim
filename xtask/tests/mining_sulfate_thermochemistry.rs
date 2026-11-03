//! Gypsum and anhydrite from the scanned USGS B2131 primary tables.
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
const USGS: &str = "source-elemental-reference-usgs-b2131";

struct Sulfate {
    pack: &'static str,
    phase: &'static str,
    mass_kg_per_mol: f64,
    formation_j_per_mol: f64,
    entropy_j_per_mol_k: f64,
    claims: usize,
}

// Literal reference rows: USGS Bulletin 2131, printed page 28.
// https://digital.library.unt.edu/iiif/ark:/67531/metadc1114954/m1/34/full/max/0/default.jpg
const GYPSUM: Sulfate = Sulfate {
    pack: "gypsum-solid-usgs-b2131-thermochemistry",
    phase: "source-phase-gypsum",
    mass_kg_per_mol: 0.172172,
    formation_j_per_mol: -2023000.0,
    entropy_j_per_mol_k: 193.8,
    claims: 3,
};
const ANHYDRITE: Sulfate = Sulfate {
    pack: "anhydrite-solid-usgs-b2131-thermochemistry",
    phase: "source-phase-anhydrite",
    mass_kg_per_mol: 0.136142,
    formation_j_per_mol: -1434400.0,
    entropy_j_per_mol_k: 107.4,
    claims: 4,
};

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

fn compile(sulfate: &Sulfate) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!(
        "sulfate-thermochemistry-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{}/manifest.tsv", sulfate.pack)))
        .arg("--out")
        .arg(&out)
        .current_dir(workspace())
        .env("CARGO_WORKSPACE_DIR", workspace())
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        sulfate.pack,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), sulfate.pack);
    assert_eq!(pack.claims().claims_ordered().count(), sulfate.claims);
    pack
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(sulfate: &Sulfate, temperature: f64) -> QueryPoint {
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
        (USGS, 1.0),
        ("source-phase-solid", 1.0),
        (sulfate.phase, 1.0),
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

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn sulfate_reference_claims_compile_and_retain_literal_usgs_values() {
    for sulfate in [&GYPSUM, &ANHYDRITE] {
        let pack = compile(sulfate);
        for (name, dims, expected) in [
            ("molar-mass", MOLAR_MASS, sulfate.mass_kg_per_mol),
            (
                "standard-molar-formation-enthalpy",
                MOLAR_ENERGY,
                sulfate.formation_j_per_mol,
            ),
            (
                "standard-molar-entropy",
                MOLAR_ENTROPY,
                sulfate.entropy_j_per_mol_k,
            ),
        ] {
            let answer = pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &point(sulfate, 298.15),
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
            let claim = pack.claims().claim(answer.receipt.selected).unwrap();
            assert_eq!(claim.validity.bound(USGS), Some((1.0, 1.0)));
        }
    }
}

#[test]
fn sulfate_claims_refuse_other_hydration_states_and_reference_conditions() {
    let hf = PropertyKey::new("standard-molar-formation-enthalpy", MOLAR_ENERGY);
    let entropy = PropertyKey::new("standard-molar-entropy", MOLAR_ENTROPY);
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    let mass = PropertyKey::new("molar-mass", MOLAR_MASS);
    for (sulfate, other_phase) in [(&GYPSUM, &ANHYDRITE), (&ANHYDRITE, &GYPSUM)] {
        let pack = compile(sulfate);
        for bad in [
            QueryPoint::new().with("temperature", 298.15).unwrap(),
            point(other_phase, 298.15),
            point(sulfate, 300.0),
            point(sulfate, 298.15)
                .with("source-pure-species", 0.0)
                .unwrap(),
            point(sulfate, 298.15)
                .with("source-specimen-qualified", 1.0)
                .unwrap(),
            point(sulfate, 298.15)
                .with_quantity("pressure", semantic(QuantityKind::Pressure), 101325.0)
                .unwrap(),
            point(sulfate, 298.15)
                .with(USGS, 0.0)
                .unwrap()
                .with("source-elemental-reference-janaf", 1.0)
                .unwrap(),
        ] {
            for property in [&hf, &entropy] {
                assert!(
                    pack.claims()
                        .query_typed(property, &bad, SelectionPolicy::SingleClaimOnly)
                        .is_err(),
                    "{}",
                    sulfate.pack
                );
            }
        }
        for temperature in [298.0, 1000.01] {
            for property in [&cp, &mass] {
                assert!(
                    pack.claims()
                        .query_typed(
                            property,
                            &point(sulfate, temperature),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err()
                );
            }
        }
        if sulfate.pack == GYPSUM.pack {
            // The selected gypsum reference table does not license borrowing
            // anhydrite or wallboard Cp, even at the reference temperature.
            assert!(
                pack.claims()
                    .query_typed(
                        &cp,
                        &point(sulfate, 298.15),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
}

#[test]
fn anhydrite_cp_and_sensible_heat_match_independent_printed_columns() {
    let pack = compile(&ANHYDRITE);
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let line = pack
        .claims()
        .query_line_typed(
            &cp,
            "temperature",
            &point(&ANHYDRITE, 298.15),
            &point(&ANHYDRITE, 1000.0),
            selection,
        )
        .unwrap();
    assert_eq!(line.knots.len(), 9);
    // USGS B2131 printed page 285: original molar Cp entries.
    // https://digital.library.unt.edu/iiif/ark:/67531/metadc1114954/m1/291/full/max/0/default.jpg
    for ((t, answer), (source_t, source_cp_molar)) in line.knots.iter().zip([
        (298.15, 101.23),
        (300.0, 101.57),
        (400.0, 116.68),
        (500.0, 127.18),
        (600.0, 135.03),
        (700.0, 141.54),
        (800.0, 147.55),
        (900.0, 153.59),
        (1000.0, 160.04),
    ]) {
        assert_eq!(*t, source_t);
        close(
            answer.evidence.value.value * ANHYDRITE.mass_kg_per_mol,
            source_cp_molar,
            2e-13,
        );
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        pack.claims().verify_receipt(&answer.receipt).unwrap();
    }
    let mass_at_upper_endpoint = pack
        .claims()
        .query_typed(
            &PropertyKey::new("molar-mass", MOLAR_MASS),
            &point(&ANHYDRITE, 1000.0),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    close(
        mass_at_upper_endpoint.evidence.value.value,
        ANHYDRITE.mass_kg_per_mol,
        2e-14,
    );
    let heating = pack
        .claims()
        .query_line_typed(
            &cp,
            "temperature",
            &point(&ANHYDRITE, 300.0),
            &point(&ANHYDRITE, 400.0),
            selection,
        )
        .unwrap();
    // Authored 2 kg, constant-pressure heating without a chemical/phase change.
    let mass_kg = 2.0;
    let heat_j: f64 = mass_kg
        * heating
            .knots
            .windows(2)
            .map(|w| {
                0.5 * (w[1].0 - w[0].0)
                    * (w[0].1.evidence.value.value + w[1].1.evidence.value.value)
            })
            .sum::<f64>();
    // The independent printed column is (H(T)-H(298.15))/T [J/mol/K].
    // Multiply each entry by its own temperature before taking the difference.
    let source_heat_j_per_mol = 400.0 * 27.87 - 300.0 * 0.63;
    let expected = mass_kg / ANHYDRITE.mass_kg_per_mol * source_heat_j_per_mol;
    close(heat_j, expected, 0.005);
    assert!(heat_j > 0.0);
    // This finite regression tolerates the coarse linearization; it does not
    // establish a global interpolation bound or source uncertainty model.
}
