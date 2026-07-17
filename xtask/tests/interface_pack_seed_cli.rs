//! G3 CLI conformance for the committed source-bound interface seed packs.

#![deny(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedInterfacePack, PropertyValue, UncertaintyModel};

const INTERFACE_COMPILER_ID: &str = "frankensim-matdb-interface-pack-compiler-v1";
const DRY_52100_MANIFEST: &str = "data/matdb/seed-v1/nasa-52100-dry-air-interface/manifest.tsv";
const GREASED_52100_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-52100-gxl320a-vacuum-interface/manifest.tsv";
const JOURNAL_4340_BRONZE_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-tn-d-2223-4340-high-lead-bronze-journal/manifest.tsv";
const NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

fn workspace_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a workspace parent")
        .join(relative)
}

fn fixture_dir() -> PathBuf {
    loop {
        let sequence = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "frankensim-interface-seed-cli-test-{}-{sequence}",
            std::process::id()
        ));
        match fs::create_dir(&path) {
            Ok(()) => return path,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => panic!("unique fixture directory: {error}"),
        }
    }
}

fn run_compiler(manifest: &Path, output: &Path) -> Output {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a workspace parent");
    Command::new(env!("CARGO_BIN_EXE_xtask"))
        .arg("matdb-pack")
        .arg("--manifest")
        .arg(manifest)
        .arg("--out")
        .arg(output)
        .env("CARGO_WORKSPACE_DIR", workspace)
        .output()
        .expect("run xtask matdb-pack")
}

fn compile_twice(manifest: &str) -> (NormalizedInterfacePack, String) {
    let directory = fixture_dir();
    let first_path = directory.join("first.fsintpk");
    let second_path = directory.join("second.fsintpk");
    let manifest = workspace_path(manifest);

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first interface-seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second interface-seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "decision stream moved");

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    let expected_prefix =
        format!("{{\"check\":\"matdb-pack\",\"compiler\":\"{INTERFACE_COMPILER_ID}\",");
    assert!(!decisions.is_empty(), "compiler emitted no decisions");
    assert!(
        decisions
            .lines()
            .all(|row| row.starts_with(&expected_prefix)),
        "decision row used the wrong compiler identity:\n{decisions}"
    );
    assert!(decisions.contains("\"reason_code\":\"runtime_interface_pack_self_verified\""));

    let first_bytes = fs::read(first_path).expect("read first interface seed pack");
    let second_bytes = fs::read(second_path).expect("read second interface seed pack");
    assert_eq!(first_bytes, second_bytes, "published interface bytes moved");
    assert_eq!(&first_bytes[..8], b"FSINTPK\0");

    let decoded = NormalizedInterfacePack::from_bytes(&first_bytes)
        .expect("compiler output re-admits at runtime");
    NormalizedInterfacePack::from_bytes_verified(decoded.content_hash(), &first_bytes)
        .expect("externally pinned interface bytes re-admit");
    assert_eq!(decoded.compiler(), INTERFACE_COMPILER_ID);
    assert!(decoded.card().models().is_empty(), "v1 carries no models");
    (decoded, decisions)
}

fn scalar(pack: &NormalizedInterfacePack, property: &str) -> f64 {
    let claims = pack.card().claims_for(property);
    assert_eq!(claims.len(), 1, "expected one {property} claim");
    let claim = claims[0].1;
    let PropertyValue::Scalar { value, .. } = &claim.value else {
        panic!("{property} was not scalar");
    };
    assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
    *value
}

#[test]
fn g3_cli_compiles_committed_nasa_52100_dry_air_interface() {
    let (pack, decisions) = compile_twice(DRY_52100_MANIFEST);

    assert_eq!(
        pack.pack_id(),
        "nasa-20150020881-aisi-52100-dry-air-sliding-interface"
    );
    assert_eq!(
        pack.card().surface_a().material.chemistry,
        "SAE 52100 source-Table-1 composition"
    );
    assert_eq!(
        pack.card().surface_b().material.chemistry,
        "SAE 52100 source-Table-1 composition"
    );
    assert_ne!(
        pack.card().surface_a().texture_frame,
        pack.card().surface_b().texture_frame,
        "ordered rider and disk surfaces collapsed"
    );
    assert_eq!(pack.card().medium(), "dry-unlubricated-contact");
    assert_eq!(pack.card().third_body(), None);
    assert_eq!(pack.card().environment(), "air-at-760-mmHg");
    assert_eq!(pack.card().history(), "fresh-specimens-one-hour-run");
    assert_eq!(pack.claims_pack().claims().claim_count(), 1);
    assert_eq!(scalar(&pack, "kinetic-friction-coefficient"), 0.45);

    let claim = pack.card().claims_for("kinetic-friction-coefficient")[0].1;
    for (axis, value) in [
        ("temperature", 297.038_888_888_888_9),
        ("ambient_pressure", 101_325.0),
        ("sliding_speed", 1.981_2),
        ("normal_force", 9.806_65),
        ("test_duration", 3_600.0),
        ("rider_radius", 0.004_762_5),
        ("disk_diameter", 0.050_8),
        ("source_atmosphere_air", 1.0),
    ] {
        assert_eq!(claim.validity.bound(axis), Some((value, value)));
    }
    assert_eq!(claim.validity.bounds().len(), 8);
    for refused_property in [
        "static-friction-coefficient",
        "wear-rate",
        "transferable-friction-law",
    ] {
        assert!(
            pack.card().claims_for(refused_property).is_empty(),
            "dry source crossed the {refused_property} no-claim boundary"
        );
    }
    assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
}

#[test]
fn g3_cli_compiles_committed_nasa_52100_gxl320a_vacuum_interface() {
    let (pack, decisions) = compile_twice(GREASED_52100_MANIFEST);

    assert_eq!(
        pack.pack_id(),
        "nasa-tm-1999-209064-aisi-52100-gxl320a-vacuum-interface"
    );
    assert_eq!(
        pack.card().surface_a().material.chemistry,
        "AISI 52100 chrome bearing steel"
    );
    assert_eq!(
        pack.card().surface_b().material.chemistry,
        "AISI 52100 chrome bearing steel"
    );
    assert_ne!(
        pack.card().surface_a().texture_frame,
        pack.card().surface_b().texture_frame,
        "rotating and stationary ball roles collapsed"
    );
    assert_eq!(pack.card().medium(), "boundary-lubricated-contact");
    assert_eq!(pack.card().third_body(), Some("GXL-320A-base-grease"));
    assert_eq!(pack.card().environment(), "vacuum-below-6.7e-4-Pa");
    assert_eq!(
        pack.card().history(),
        "four-hour-test-with-hourly-wear-measurements"
    );
    assert_eq!(pack.claims_pack().claims().claim_count(), 3);
    assert_eq!(scalar(&pack, "kinetic-friction-coefficient"), 0.11);
    assert_eq!(
        scalar(&pack, "kinetic-friction-coefficient-observed-minimum"),
        0.09
    );
    assert_eq!(
        scalar(&pack, "kinetic-friction-coefficient-observed-maximum"),
        0.23
    );

    for property in [
        "kinetic-friction-coefficient",
        "kinetic-friction-coefficient-observed-minimum",
        "kinetic-friction-coefficient-observed-maximum",
    ] {
        let claim = pack.card().claims_for(property)[0].1;
        for (axis, value) in [
            ("normal_force", 200.0),
            ("sliding_speed", 0.028_8),
            ("initial_hertz_mean_stress", 3_500_000_000.0),
            ("test_duration", 14_400.0),
            ("number_of_runs", 4.0),
            ("ball_diameter", 0.009_5),
            ("source_temperature_approximately_23c", 1.0),
            ("source_pressure_less_than_6.7e-4_pa", 1.0),
        ] {
            assert_eq!(claim.validity.bound(axis), Some((value, value)));
        }
        assert_eq!(claim.validity.bounds().len(), 8);
    }
    for refused_property in ["wear-rate", "friction-law", "bearing-life"] {
        assert!(
            pack.card().claims_for(refused_property).is_empty(),
            "greased source crossed the {refused_property} no-claim boundary"
        );
    }
    assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
}

#[test]
fn g3_cli_compiles_committed_nasa_4340_high_lead_bronze_journal_interface() {
    let (pack, decisions) = compile_twice(JOURNAL_4340_BRONZE_MANIFEST);

    assert_eq!(
        pack.pack_id(),
        "nasa-tn-d-2223-sae-4340-high-lead-bronze-hexane-journal-interface"
    );
    assert_eq!(pack.card().surface_a().material.chemistry, "SAE 4340 steel");
    assert_eq!(
        pack.card().surface_b().material.chemistry,
        "high-lead bearing bronze 70wt%Cu-26wt%Pb-4wt%Sn"
    );
    assert_eq!(pack.card().medium(), "hexane");
    assert_eq!(pack.card().third_body(), None);
    assert_eq!(
        pack.card().environment(),
        "source-surrounding-environment-unstated"
    );
    assert_eq!(
        pack.card().history(),
        "table-I-220-psi-screening-run-test-objective-attained"
    );
    assert_eq!(pack.claims_pack().claims().claim_count(), 1);
    assert_eq!(
        scalar(&pack, "maximum-demonstrated-unit-bearing-load"),
        1_516_846.604_497_039_5
    );

    let claim = pack
        .card()
        .claims_for("maximum-demonstrated-unit-bearing-load")[0]
        .1;
    for (axis, bounds) in [
        ("journal_rotation_frequency", (250.0, 250.0)),
        ("journal_surface_speed", (32.613_6, 32.613_6)),
        ("time_at_maximum_speed_and_load", (1_800.0, 1_800.0)),
        ("total_test_time", (5_400.0, 5_400.0)),
        ("bearing_inside_diameter", (0.038_1, 0.038_1)),
        ("bearing_length", (0.038_1, 0.038_1)),
        (
            "room_temperature_diametral_clearance",
            (0.000_038_1, 0.000_038_1),
        ),
        ("journal_surface_finish_rms", (0.000_000_127, 0.000_000_254)),
        (
            "hexane_inlet_gauge_pressure",
            (68_947.572_931_683_61, 68_947.572_931_683_61),
        ),
        ("hexane_inlet_hole_diameter", (0.003_175, 0.003_175)),
        (
            "hexane_dynamic_viscosity_at_source_reference_temperature",
            (0.000_296_474_563_606_239, 0.000_296_474_563_606_239),
        ),
        (
            "hexane_viscosity_reference_temperature",
            (297.038_888_888_888_9, 297.038_888_888_888_9),
        ),
        ("source_groove_type_figure_2a", (1.0, 1.0)),
    ] {
        assert_eq!(claim.validity.bound(axis), Some(bounds));
    }
    assert_eq!(claim.validity.bounds().len(), 13);
    for refused_property in [
        "kinetic-friction-coefficient",
        "wear-rate",
        "design-allowable-bearing-load",
        "transferable-journal-bearing-law",
    ] {
        assert!(
            pack.card().claims_for(refused_property).is_empty(),
            "journal source crossed the {refused_property} no-claim boundary"
        );
    }
    assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
}
