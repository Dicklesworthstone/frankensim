#![deny(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::NormalizedPack;

static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

const PACK_BYTES_GOLDEN: usize = 3_177;
const PACK_HASH_GOLDEN: &str = "c1fb2f443708d297423179f4ac6024ee26b1d0c940a229d1d9084726ccbd2bc5";

const MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tfixture-alloy-x\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture handbook table 7\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tsource.tsv\tmaterial-tsv-v1\n",
);

const SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tcoupon\talloy-X-solution-treated\tASTM-fixture\tjoint coupon series\n",
    "scalar\tdensity\tcoupon\tdensity\t7.85\tg/cm3\tconstant\n",
    "uncertainty\tdensity\tabsolute\t0.005\tg/cm3\t0.95\t1\n",
    "validity\tdensity\ttemperature\t0\t100\tdegC\n",
    "curve\tmodulus\tcoupon\tyoung_modulus\ttemperature\tdegC\tGPa\t0:210,100:202\tlinear\n",
    "uncertainty\tmodulus\trelative\t2\t%\t0.95\t1\n",
    "validity\tmodulus\ttemperature\t0\t100\tdegC\n",
    "frame\tmodulus\tspecimen\tlab\n",
    "joint\tcoupon\tdensity-modulus\tdensity:scalar,modulus:y:0\t0.000025,0,0.000009\t1,0,1\t1\n",
);

fn fixture_dir() -> PathBuf {
    loop {
        let sequence = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "frankensim-matdb-pack-cli-test-{}-{sequence}",
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
    fs::write(&manifest, MANIFEST).expect("write manifest fixture");
    fs::write(directory.join("source.tsv"), source).expect("write source fixture");
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

#[test]
fn g3_cli_compiles_two_identical_pinned_packs() {
    let (directory, manifest) = write_fixture(SOURCE);
    let first_path = directory.join("first.fsmatpk");
    let second_path = directory.join("second.fsmatpk");

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

    let first_bytes = fs::read(&first_path).expect("read first normalized pack");
    let second_bytes = fs::read(&second_path).expect("read second normalized pack");
    assert_eq!(first_bytes, second_bytes, "published pack bytes moved");
    assert_eq!(first_bytes.len(), PACK_BYTES_GOLDEN);
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode compiler output");
    assert_eq!(decoded.content_hash().to_string(), PACK_HASH_GOLDEN);
    NormalizedPack::from_bytes_verified(decoded.content_hash(), &first_bytes)
        .expect("externally pinned bytes re-admit");

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    let rows: Vec<_> = decisions.lines().collect();
    assert!(!rows.is_empty());
    assert!(
        rows.iter()
            .all(|row| row.starts_with("{\"check\":\"matdb-pack\""))
    );
    assert!(
        rows.iter()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{PACK_HASH_GOLDEN}\"")))
    );
    assert!(decisions.contains("\"reason_code\":\"published_new_verified_artifact\""));
    assert!(decisions.contains("\"reason_code\":\"joint_statistics_normalized\""));
}

#[test]
fn g3_cli_retains_prior_admissions_when_later_claim_refuses() {
    let rejected_source = SOURCE.replacen(
        "uncertainty\tmodulus\trelative\t2\t%\t0.95\t1\n",
        "uncertainty\tmodulus\tabsolute\t2\tkg\t0.95\t1\n",
        1,
    );
    assert_ne!(rejected_source, SOURCE);
    let (directory, manifest) = write_fixture(&rejected_source);
    let first_path = directory.join("refused-first.fsmatpk");
    let second_path = directory.join("refused-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(!first.status.success());
    assert!(!second.status.success());
    assert_eq!(first.stdout, second.stdout, "refusal transcript moved");
    assert!(!first_path.exists(), "refused compilation published output");
    assert!(
        !second_path.exists(),
        "refused compilation published output"
    );

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"subject\":\"source:primary\""));
    assert!(decisions.contains("\"subject\":\"claim:density\""));
    assert!(decisions.contains("\"reason_code\":\"claim_normalized\""));
    assert_eq!(decisions.matches("\"verdict\":\"refuse\"").count(), 1);
    let refusal = decisions
        .lines()
        .find(|row| row.contains("\"verdict\":\"refuse\""))
        .expect("one terminal refusal row");
    assert!(refusal.contains("\"reason_code\":\"uncertainty_dims_mismatch\""));
    assert!(!refusal.contains("\"source_hash\":\"\""));
    assert!(refusal.contains("\"pack_hash\":\"\""));
    assert!(
        decisions
            .lines()
            .filter(|row| row.contains("\"verdict\":\"admit\""))
            .all(|row| row.contains("\"pack_hash\":\"\""))
    );
    assert!(
        String::from_utf8_lossy(&first.stderr)
            .contains("error: matdb pack refused [uncertainty_dims_mismatch]")
    );
}
