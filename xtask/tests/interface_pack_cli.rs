//! G3 CLI conformance for the bounded `interface-tsv-v1` compiler profile.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::NormalizedInterfacePack;

const INTERFACE_COMPILER_ID: &str = "frankensim-matdb-interface-pack-compiler-v1";
const NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

const MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tfixture-steel-bronze-journal-interface\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture pin-on-disk campaign POD-19\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tinterface.tsv\tinterface-tsv-v1\n",
);

const SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "surface_a\tAISI-52100\ttempered-martensite\thardened-60HRC\t2\tjournal-ground-frame-3\n",
    "surface_b\tC93200\tcast-bearing-bronze\tmachined-bore\t1\tbore-honed-frame-8\n",
    "context\toil-film\tnamed-reference-oil-lot-4\tlaboratory-air\trun-in-1000-cycles\n",
    "observation\tpod19\tordered steel journal on bronze bearing coupon\tPOD-19 run-in campaign\tfixture value; not a seed-dataset authority\n",
    "scalar\tfriction\tpod19\tkinetic_friction_coefficient\t0.08\t1\tconstant\n",
    "uncertainty\tfriction\tunstated\t-\t-\t-\t-\n",
    "validity\tfriction\ttemperature\t293.15\t313.15\tK\n",
    "validity\tfriction\tnormal_pressure\t100000\t5000000\tPa\n",
);

fn fixture_dir() -> PathBuf {
    loop {
        let sequence = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "frankensim-interface-pack-cli-test-{}-{sequence}",
            std::process::id()
        ));
        match fs::create_dir(&path) {
            Ok(()) => return path,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => panic!("unique fixture directory: {error}"),
        }
    }
}

fn write_fixture(source: &str) -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, MANIFEST).expect("write interface manifest fixture");
    fs::write(directory.join("interface.tsv"), source).expect("write interface source fixture");
    (directory, manifest)
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

fn decision_text(output: &Output) -> &str {
    std::str::from_utf8(&output.stdout).expect("decision stream is UTF-8")
}

#[test]
fn g3_cli_compiles_identical_verified_ordered_interface_packs() {
    let (directory, manifest) = write_fixture(SOURCE);
    let first_path = directory.join("first.fsintpk");
    let second_path = directory.join("second.fsintpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first compiler run failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second compiler run failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "decision stream moved");
    let decisions = decision_text(&first);
    let expected_prefix =
        format!("{{\"check\":\"matdb-pack\",\"compiler\":\"{INTERFACE_COMPILER_ID}\",");
    assert!(
        decisions
            .lines()
            .all(|row| row.starts_with(&expected_prefix)),
        "decision row used the wrong compiler identity:\n{decisions}"
    );
    assert!(decisions.contains("\"reason_code\":\"ordered_surface_identity_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_interface_pack_self_verified\""));

    let first_bytes = fs::read(&first_path).expect("read first interface pack");
    let second_bytes = fs::read(&second_path).expect("read second interface pack");
    assert_eq!(first_bytes, second_bytes, "published interface bytes moved");
    assert_eq!(&first_bytes[..8], b"FSINTPK\0");

    let decoded = NormalizedInterfacePack::from_bytes(&first_bytes)
        .expect("compiler output re-admits at runtime");
    NormalizedInterfacePack::from_bytes_verified(decoded.content_hash(), &first_bytes)
        .expect("externally pinned interface bytes re-admit");
    assert_eq!(decoded.pack_id(), "fixture-steel-bronze-journal-interface");
    assert_eq!(decoded.compiler(), INTERFACE_COMPILER_ID);
    assert_eq!(decoded.card().surface_a().material.chemistry, "AISI-52100");
    assert_eq!(decoded.card().surface_b().material.chemistry, "C93200");
    assert_eq!(decoded.card().medium(), "oil-film");
    assert_eq!(
        decoded.card().third_body(),
        Some("named-reference-oil-lot-4")
    );
    assert_eq!(decoded.card().environment(), "laboratory-air");
    assert_eq!(decoded.card().history(), "run-in-1000-cycles");
    assert_eq!(
        decoded
            .card()
            .claims_for("kinetic_friction_coefficient")
            .len(),
        1
    );
    assert!(decoded.card().models().is_empty(), "v1 carries no models");
}

#[test]
fn g3_cli_refuses_an_interface_without_complete_system_identity() {
    let incomplete = SOURCE.replace(
        "context\toil-film\tnamed-reference-oil-lot-4\tlaboratory-air\trun-in-1000-cycles\n",
        "",
    );
    let (directory, manifest) = write_fixture(&incomplete);
    let output_path = directory.join("refused.fsintpk");
    let output = run_compiler(&manifest, &output_path);

    assert!(
        !output.status.success(),
        "incomplete identity was published"
    );
    assert!(!output_path.exists(), "refusal left a partial artifact");
    let decisions = decision_text(&output);
    assert!(decisions.contains(&format!("\"compiler\":\"{INTERFACE_COMPILER_ID}\"")));
    assert!(decisions.contains("\"reason_code\":\"missing_interface_context\""));
}

#[test]
fn g3_cli_refuses_noncanonical_material_revisions() {
    let noncanonical = SOURCE.replacen(
        "\t2\tjournal-ground-frame-3",
        "\t02\tjournal-ground-frame-3",
        1,
    );
    let (directory, manifest) = write_fixture(&noncanonical);
    let output_path = directory.join("refused-revision.fsintpk");
    let output = run_compiler(&manifest, &output_path);

    assert!(
        !output.status.success(),
        "noncanonical revision was published"
    );
    assert!(!output_path.exists(), "refusal left a partial artifact");
    let decisions = decision_text(&output);
    assert!(decisions.contains("\"reason_code\":\"noncanonical_unsigned_integer\""));
}
