//! G0/G3/G5: source mineral densities survive the real compiler and query path.
//! Mineral reference values do not qualify bulk ores, work indices or recovery.
#![deny(unsafe_code)]

use fs_matdb::{
    MatDbError, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Dims, QuantitySpec};
use std::{path::Path, process::Command};

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SOURCE_FLAGS: [(&str, f64); 4] = [
    ("mineral-reference-model", 1.0),
    ("source-test-temperature-known", 0.0),
    ("source-test-pressure-known", 0.0),
    ("bulk-ore-qualified", 0.0),
];

fn compile(mineral: &str, claims: usize) -> NormalizedPack {
    let slug = format!("{mineral}-handbook-mineralogy-reference");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "frankensim-mining-minerals-{}-{nonce}-{slug}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), claims);
    pack
}

fn reference_point() -> QueryPoint {
    SOURCE_FLAGS
        .iter()
        .fold(QueryPoint::new(), |point, &(axis, value)| {
            point
                .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
                .unwrap()
        })
}

fn density_key(property: &str) -> PropertyKey {
    PropertyKey::with_quantity(property, QuantitySpec::dimensional(DENSITY))
}

fn assert_density(pack: &NormalizedPack, property: &str, expected_si: f64) {
    let answer = pack
        .claims()
        .query_typed(
            &density_key(property),
            &reference_point(),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!((answer.evidence.value.value - expected_si).abs() <= 1e-12 * expected_si);
    assert_eq!(answer.evidence.value.dims, DENSITY);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
}

fn assert_reference_boundaries(pack: &NormalizedPack, property: &str) {
    // Source conditions are required, and a reference cannot become a qualified ore.
    let missing = pack.claims().query_typed(
        &density_key(property),
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    );
    assert!(
        matches!(
            &missing,
            Err(MatDbError::AxisQuantityMismatch {
                axis,
                expected: Some(quantity),
                found: None,
                ..
            }) if axis == "bulk-ore-qualified"
                && *quantity == QuantitySpec::dimensional(Dims::NONE)
        ),
        "missing source coordinates must retain their typed-axis refusal: {missing:?}"
    );
    for (axis, value) in SOURCE_FLAGS {
        let point = reference_point()
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0 - value)
            .unwrap();
        assert!(matches!(
            pack.claims().query_typed(
                &density_key(property),
                &point,
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    assert!(matches!(
        pack.claims().query_typed(
            &PropertyKey::with_quantity(
                "bond-ball-mill-work-index",
                QuantitySpec::dimensional(Dims([2, 0, -2, 0, 0, 0])),
            ),
            &reference_point(),
            SelectionPolicy::SingleClaimOnly,
        ),
        Err(MatDbError::UnknownProperty { .. })
    ));
}

#[test]
fn g0_g3_reported_mineral_densities_keep_identity_and_unknown_conditions() {
    // Handbook Physical Properties D(meas.), converted from g/cm3 to kg/m3.
    for (mineral, expected_si) in [
        ("pyrite", 5018.0),
        ("dolomite", 2860.0),
        ("magnetite", 5175.0),
    ] {
        let pack = compile(mineral, 1);
        assert_density(&pack, "density", expected_si);
        assert_reference_boundaries(&pack, "density");
    }
}

#[test]
fn g0_g3_reported_ranges_remain_endpoints_without_a_nominal_density() {
    for (mineral, lower, upper) in [
        ("chalcopyrite", 4100.0, 4300.0),
        ("chalcocite", 5500.0, 5800.0),
        ("cassiterite", 6980.0, 7010.0),
    ] {
        let pack = compile(mineral, 2);
        assert_density(&pack, "density-reported-range-minimum", lower);
        assert_density(&pack, "density-reported-range-maximum", upper);
        assert_reference_boundaries(&pack, "density-reported-range-minimum");
        assert!(matches!(
            pack.claims().query_typed(
                &density_key("density"),
                &reference_point(),
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(MatDbError::UnknownProperty { .. })
        ));
    }
}

#[test]
fn g5_source_pack_compilation_is_byte_identical() {
    let first = compile("cassiterite", 2);
    let second = compile("cassiterite", 2);
    assert_eq!(first.content_hash(), second.content_hash());
    assert_eq!(first.to_bytes(), second.to_bytes());
}

#[test]
fn catalog_exposes_range_endpoints_and_refuses_an_invented_test_temperature() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let mut args = vec![
        "matdb-catalog".to_owned(),
        "--pack".to_owned(),
        "chalcocite-handbook-mineralogy-reference".to_owned(),
        "--property".to_owned(),
        "density-reported-range-minimum".to_owned(),
        "--property".to_owned(),
        "density-reported-range-maximum".to_owned(),
    ];
    for (axis, value) in SOURCE_FLAGS {
        args.extend(["--at".to_owned(), format!("{axis}={value}")]);
    }
    let run = |args: &[String]| {
        Command::new(env!("CARGO_BIN_EXE_xtask"))
            .args(args)
            .env("CARGO_WORKSPACE_DIR", root)
            .current_dir(root)
            .output()
            .unwrap()
    };
    let result = run(&args);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let output = String::from_utf8(result.stdout).unwrap();
    assert!(output.contains("density-reported-range-minimum"));
    assert!(output.contains("density-reported-range-maximum"));
    args.extend(["--at".to_owned(), "temperature=298.15".to_owned()]);
    let result = run(&args);
    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
}
