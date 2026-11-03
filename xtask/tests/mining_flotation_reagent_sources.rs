//! G0/G3: source-backed frother properties and explicit dosing/phase bases.
#![deny(unsafe_code)]

use fs_matdb::{
    MatDbError, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const ID: &str = "mibc-arkema-flotation-reference";
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);

fn compile() -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!("{ID}-{}-{nonce}.fsmatpk", std::process::id()));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{ID}/manifest.tsv")))
        .arg("--out")
        .arg(&out)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), ID);
    assert_eq!(pack.claims().claims_ordered().count(), 5);
    pack
}

fn product() -> QueryPoint {
    [
        ("source-arkema-mibc-product", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-lot-qualified", 0.0),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |p, (name, value)| {
        p.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

fn liquid(temperature: f64) -> QueryPoint {
    product()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "source-test-pressure-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
}

fn value(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) -> f64 {
    let result = pack
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(
        result.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(result.receipt.observation_backed);
    assert!(!result.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&result.receipt).unwrap();
    result.evidence.value.value
}

#[test]
fn source_density_purity_and_directional_solubilities_keep_original_bases() {
    let pack = compile();
    let p = liquid(293.15);
    for (name, dims, expected) in [
        ("density-reported-range-minimum", DENSITY, 806.0),
        ("density-reported-range-maximum", DENSITY, 808.0),
        (
            "water-rich-phase-mibc-saturation-mass-fraction",
            Dims::NONE,
            0.017,
        ),
        (
            "mibc-rich-phase-water-saturation-mass-fraction",
            Dims::NONE,
            0.058,
        ),
    ] {
        let actual = value(&pack, name, dims, &p);
        assert!((actual - expected).abs() <= 1e-13 * expected.abs().max(1.0));
    }
    assert_eq!(
        value(
            &pack,
            "mibc-product-purity-mass-fraction-minimum",
            Dims::NONE,
            &product()
        ),
        0.985
    );
    // An authored volume flow of 1 L/min and 60 dry tonnes/h: product kg/t,
    // not a recommended frother dose or a forecast of mineral recovery.
    let volume_m3_s = 0.001 / 60.0;
    let dry_kg_s = 60000.0 / 3600.0;
    for (name, expected_kg_t) in [
        ("density-reported-range-minimum", 0.806),
        ("density-reported-range-maximum", 0.808),
    ] {
        let product_kg_t = volume_m3_s * value(&pack, name, DENSITY, &p) / dry_kg_s * 1000.0;
        assert!((product_kg_t - expected_kg_t).abs() < 1e-14);
    }
}

#[test]
fn frother_reference_refuses_invented_temperature_lot_and_recovery_properties() {
    let pack = compile();
    let key = PropertyKey::new("density-reported-range-minimum", DENSITY);
    for temperature in [293.14, 293.16, 298.15] {
        assert!(matches!(
            pack.claims()
                .query_typed(&key, &liquid(temperature), SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    for (name, wrong) in [
        ("source-arkema-mibc-product", 0.0),
        ("engineering-reference-model", 0.0),
        ("source-lot-qualified", 1.0),
        ("source-test-pressure-known", 1.0),
    ] {
        let p = liquid(293.15)
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), wrong)
            .unwrap();
        assert!(matches!(
            pack.claims()
                .query_typed(&key, &p, SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    assert!(
        pack.claims()
            .query_typed(&key, &product(), SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    for name in [
        "density",
        "flotation-rate-constant",
        "surface-tension",
        "standard-molar-formation-enthalpy",
    ] {
        assert!(matches!(
            pack.claims().query_typed(
                &PropertyKey::new(name, Dims::NONE),
                &liquid(293.15),
                SelectionPolicy::SingleClaimOnly
            ),
            Err(MatDbError::UnknownProperty { .. })
        ));
    }
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::new("water-rich-phase-mibc-saturation-mass-fraction", DENSITY),
                &liquid(293.15),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}

#[test]
fn frother_pack_is_deterministic_and_discoverable_by_ordinary_catalog() {
    assert_eq!(compile().to_bytes(), compile().to_bytes());
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let args = [
        "matdb-catalog",
        "--pack",
        ID,
        "--property",
        "density-reported-range-minimum",
        "--at",
        "temperature=293.15",
        "--at",
        "source-arkema-mibc-product=1",
        "--at",
        "engineering-reference-model=1",
        "--at",
        "source-lot-qualified=0",
        "--at",
        "source-test-pressure-known=0",
    ];
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(args)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let output = String::from_utf8(run.stdout).unwrap();
    assert!(output.contains(ID));
    assert!(output.contains("density-reported-range-minimum"));
}
