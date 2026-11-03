//! Actual-source checks for the dilute-water oxygen Henry reference.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    AmountConcentration, Dims, Pressure, QtyAny, QuantitySpec,
    parse::parse_qty,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const SLUG: &str = "oxygen-water-henry-sander2023-rettich2000";
const PROPERTY: &str = "henry_molar_solubility";
const HENRY_DIMS: Dims = Dims([-2, -1, 2, 0, 0, 1]);
const FLAGS: [(&str, f64); 5] = [
    ("source-solvent-pure-water", 1.0),
    ("source-henry-infinite-dilution", 1.0),
    ("source-solute-molecular-oxygen", 1.0),
    ("engineering-reference-model", 1.0),
    ("source-specimen-qualified", 0.0),
];

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile() -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-oxygen-henry-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed),
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{SLUG}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{SLUG}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), SLUG);
    assert_eq!(pack.claims().claim_count(), 1);
    pack
}

fn temperature_type() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(omit: Option<&str>) -> QueryPoint {
    let mut point = QueryPoint::new();
    if omit != Some("temperature") {
        point = point
            .with_quantity("temperature", temperature_type(), 298.15)
            .unwrap();
    }
    for (axis, value) in FLAGS {
        if omit != Some(axis) {
            point = point
                .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
                .unwrap();
        }
    }
    point
}

#[test]
fn actual_henry_source_preserves_literal_units_receipts_and_reference_scope() {
    let pack = compile();
    let key = PropertyKey::new(PROPERTY, HENRY_DIMS);
    let query = point(None);
    let answer = pack
        .claims()
        .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    // Sander's normalized Rettich et al. row gives two significant digits;
    // neither more precise measurements nor other rows are substituted.
    assert_eq!(answer.evidence.value.value.to_bits(), 1.3e-5_f64.to_bits());
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();

    for (axis, expected) in FLAGS {
        for wrong in [
            point(Some(axis)),
            query
                .clone()
                .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0 - expected)
                .unwrap(),
        ] {
            assert!(
                pack.claims()
                    .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err(),
                "missing or wrong {axis} must refuse"
            );
        }
    }
    for wrong in [
        point(Some("temperature")),
        query
            .clone()
            .with_quantity("temperature", temperature_type(), 310.0)
            .unwrap(),
        query
            .clone()
            .with_quantity(
                "temperature",
                QuantitySpec::dimensional(Dims([0, 0, 0, 1, 0, 0])),
                298.15,
            )
            .unwrap(),
    ] {
        assert!(
            pack.claims()
                .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
    }
    // The coefficient cannot be queried as its reciprocal volatility
    // convention, as concentration alone, or as pressure alone.
    for wrong_dims in [
        Dims([2, 1, -2, 0, 0, -1]),
        AmountConcentration::DIMS,
        Pressure::DIMS,
    ] {
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new(PROPERTY, wrong_dims),
                    &query,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn dilute_water_oxygen_partial_pressure_gives_dimensioned_concentration() {
    let pack = compile();
    let key = PropertyKey::new(PROPERTY, HENRY_DIMS);
    let answer = pack
        .claims()
        .query_typed(&key, &point(None), SelectionPolicy::SingleClaimOnly)
        .unwrap();
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    let coefficient = QtyAny::new(answer.evidence.value.value, HENRY_DIMS);
    assert_eq!(coefficient, parse_qty("1.3e-5 mol/m^3/Pa").unwrap());

    // Explicit modeling assumptions: dilute O2 in pure water at 298.15 K,
    // and an ideal gas with the specified O2 partial pressure. This does
    // not infer a finite pressure-validity range or a transfer rate.
    let oxygen_partial_pressure = parse_qty("20 kPa").unwrap();
    assert_eq!(oxygen_partial_pressure.dims, Pressure::DIMS);
    let concentration = coefficient.try_mul(oxygen_partial_pressure).unwrap();
    assert_eq!(concentration.dims, AmountConcentration::DIMS);
    let typed: AmountConcentration = concentration.to_typed().unwrap();
    assert!((typed.value() - 0.26).abs() < 1e-15);
    assert!((concentration.value - parse_qty("0.00026 mol/L").unwrap().value).abs() < 1e-15);

    // Litre and cubic-metre coefficients share dimensions but differ by
    // a factor of 1000. A dimensional check alone cannot catch that error;
    // the literal value and unit conversion above must also agree.
    let mislabeled_coefficient = parse_qty("1.3e-5 mol/L/Pa").unwrap();
    assert_eq!(mislabeled_coefficient.dims, coefficient.dims);
    assert!((mislabeled_coefficient.value / coefficient.value - 1000.0).abs() < 1e-10);
    assert_ne!(mislabeled_coefficient.value, coefficient.value);

    // Total gas pressure and species partial pressure also share dimensions.
    // The caller must supply p_O2 explicitly; the source is not a gas-mixture
    // model and dimensional arithmetic cannot choose a component fraction.
}
