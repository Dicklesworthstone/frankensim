//! Integration regression for `.11.7` Stage 1: pins the thermal heatsink
//! run identity, report content hash, and evidence package Merkle root across
//! the factored durable stage driver (`crate::stage_driver`).

use std::path::Path;

use fs_cli::{exit, run};

const PINNED_HEATSINK_RUN_ID: &str =
    "5aa11590e72d49ec218c9d5049c71ae19c1ce5e647efa083a08d20d9e21b0a91";
const PINNED_HEATSINK_REPORT_HASH: &str =
    "5fe2c82237d11ea78a526a348f074c9e69311cbe9167e8e8fe92f4113c250399";
const PINNED_HEATSINK_PACKAGE_ROOT: &str =
    "1b42e3d70134c85c021e3e0442c75dd17e5d45a2c11b7ea5b750c533747a80b4";

fn json_field<'a>(stdout: &'a str, key: &str) -> &'a str {
    let needle = format!("\"{key}\":\"");
    stdout
        .split(&needle)
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .unwrap_or_else(|| panic!("missing JSON field `{key}` in stdout: {stdout}"))
}

#[test]
fn heatsink_run_id_report_hash_and_package_root_are_bit_for_bit_unchanged() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("examples/heatsink-fan/heatsink-fan.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");

    let dir = std::env::temp_dir().join(format!("fs-cli-stage-driver-pin-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch directory");
    let ledger = dir.join("heatsink.db");

    let imported = run([
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]
    .into_iter()
    .map(String::from));
    assert_eq!(imported.exit_code, exit::SUCCESS, "{}", imported.stderr);

    let output = run([
        "--json",
        "run",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]
    .into_iter()
    .map(String::from));
    assert_eq!(output.exit_code, exit::SUCCESS, "{}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"completed\""));
    assert!(output.stdout.contains("\"stages_completed\":7"));
    assert!(output.stdout.contains("\"checker\":\"pass\""));

    assert_eq!(json_field(&output.stdout, "run"), PINNED_HEATSINK_RUN_ID);
    assert_eq!(
        json_field(&output.stdout, "report_content_hash"),
        PINNED_HEATSINK_REPORT_HASH
    );
    assert_eq!(
        json_field(&output.stdout, "package_root"),
        PINNED_HEATSINK_PACKAGE_ROOT
    );

    let pkg = std::fs::read_to_string(json_field(&output.stdout, "package")).expect("package");
    let parsed = fs_package::EvidencePackage::from_json(&pkg).expect("format-9 package");
    assert!(fs_checker::check(&parsed).passed());
}
