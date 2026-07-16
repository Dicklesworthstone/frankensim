#![deny(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedModelPack, NormalizedPack, PropertyValue};

static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

const PACK_BYTES_GOLDEN: usize = 3_177;
const PACK_HASH_GOLDEN: &str = "c1fb2f443708d297423179f4ac6024ee26b1d0c940a229d1d9084726ccbd2bc5";
const NASA9_PACK_BYTES_GOLDEN: usize = 4_940;
const NASA9_PACK_HASH_GOLDEN: &str =
    "006177a7cc6f7b4ae10a9eb4a5bf49faaf21911ef9473190a29ecfc3a818a162";

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

const MATERIAL_FAMILIES_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tfixture-material-families\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture material-family extracts\n",
    "license\tCC-BY-4.0\n",
    "source\thandbook\thandbook.tsv\tmaterial-tsv-v1\n",
    "source\tbh-curve\tbh.tsv\tmaterial-tsv-v1\n",
    "source\tsn-curve\tsn.tsv\tmaterial-tsv-v1\n",
    "source\tlubricant\tlubricant.tsv\tmaterial-tsv-v1\n",
);

const HANDBOOK_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\thandbook-coupon\talloy-X-solution-treated\thandbook-table\tdensity extract\n",
    "scalar\thandbook-density\thandbook-coupon\tdensity\t7.85\tg/cm3\tconstant\n",
    "uncertainty\thandbook-density\trelative\t0.5\t%\t0.95\t1\n",
    "validity\thandbook-density\ttemperature\t0\t100\tdegC\n",
);

const BH_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tbh-loop\talloy-X-ring\tquasistatic-hysteresis\tdemagnetized branch\n",
    "curve\tbh-curve\tbh-loop\tmagnetic_flux_density\tmagnetic_field_strength\tA/m\tT\t0:0,100:0.2,1000:1.5\tlinear\n",
    "uncertainty\tbh-curve\trelative\t1\t%\t0.95\t1\n",
    "validity\tbh-curve\tmagnetic_field_strength\t0\t1000\tA/m\n",
    "validity\tbh-curve\ttemperature\t20\t25\tdegC\n",
);

const SN_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tsn-coupons\talloy-X-polished\tconstant-amplitude-fatigue\tfully reversed\n",
    "curve\tsn-curve\tsn-coupons\tfatigue_life\tstress_amplitude\tMPa\t1\t100:10000000,250:500000,400:20000\ttabulated\n",
    "uncertainty\tsn-curve\trelative\t5\t%\t0.90\t1\n",
    "validity\tsn-curve\tstress_amplitude\t100\t400\tMPa\n",
);

const LUBRICANT_SOURCE: &str = concat!(
    "frankensim.matdb-source.v1\n",
    "observation\tlubricant-batch\tPAO-4-batch-A\trotational-rheometer\tnew fluid\n",
    "curve\tlubricant-viscosity\tlubricant-batch\tdynamic_viscosity\ttemperature\tdegC\tPa*s\t-20:0.12,40:0.018,100:0.005\tlinear\n",
    "uncertainty\tlubricant-viscosity\trelative\t3\t%\t0.95\t1\n",
    "validity\tlubricant-viscosity\ttemperature\t-20\t100\tdegC\n",
);

const NASA9_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tN2\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture NASA-9 species table\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tnasa9.tsv\tnasa9-v1\n",
);

const NASA9_SOURCE: &str = concat!(
    "frankensim.nasa9-source.v1\n",
    "region\tN2\tlow\t-73.15\t700\tdegC\t100\tkPa\n",
    "coefficient\tN2\tlow\ta0\t0\tK^2\n",
    "coefficient\tN2\tlow\ta1\t0\tK\n",
    "coefficient\tN2\tlow\ta2\t3.5\t1\n",
    "coefficient\tN2\tlow\ta3\t0.001\tK^-1\n",
    "coefficient\tN2\tlow\ta4\t0\tK^-2\n",
    "coefficient\tN2\tlow\ta5\t0\tK^-3\n",
    "coefficient\tN2\tlow\ta6\t0\tK^-4\n",
    "coefficient\tN2\tlow\ta7\t100\tK\n",
    "coefficient\tN2\tlow\ta8\t1\t1\n",
    "region\tN2\thigh\t1000\t6000\tK\t100000\tPa\n",
    "coefficient\tN2\thigh\ta0\t0\tK^2\n",
    "coefficient\tN2\thigh\ta1\t0\tK\n",
    "coefficient\tN2\thigh\ta2\t4\t1\n",
    "coefficient\tN2\thigh\ta3\t0.0001\tK^-1\n",
    "coefficient\tN2\thigh\ta4\t0\tK^-2\n",
    "coefficient\tN2\thigh\ta5\t0\tK^-3\n",
    "coefficient\tN2\thigh\ta6\t0\tK^-4\n",
    "coefficient\tN2\thigh\ta7\t200\tK\n",
    "coefficient\tN2\thigh\ta8\t2\t1\n",
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

fn write_material_families_fixture() -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, MATERIAL_FAMILIES_MANIFEST).expect("write family manifest fixture");
    for (name, source) in [
        ("handbook.tsv", HANDBOOK_SOURCE),
        ("bh.tsv", BH_SOURCE),
        ("sn.tsv", SN_SOURCE),
        ("lubricant.tsv", LUBRICANT_SOURCE),
    ] {
        fs::write(directory.join(name), source).expect("write family source fixture");
    }
    (directory, manifest)
}

fn write_nasa9_fixture() -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, NASA9_MANIFEST).expect("write NASA-9 manifest fixture");
    fs::write(directory.join("nasa9.tsv"), NASA9_SOURCE).expect("write NASA-9 source fixture");
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
fn g3_cli_compiles_handbook_bh_sn_and_lubricant_material_claims() {
    let (directory, manifest) = write_material_families_fixture();
    let first_path = directory.join("families-first.fsmatpk");
    let second_path = directory.join("families-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first material-family compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second material-family compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "decision stream moved");

    let first_bytes = fs::read(first_path).expect("read first family pack");
    let second_bytes = fs::read(second_path).expect("read second family pack");
    assert_eq!(first_bytes, second_bytes, "material-family pack moved");
    let decoded = NormalizedPack::from_bytes_verified(
        NormalizedPack::from_bytes(&first_bytes)
            .expect("decode material-family pack")
            .content_hash(),
        &first_bytes,
    )
    .expect("verified material-family pack");
    assert_eq!(decoded.claims().claim_count(), 4);

    for (property, source, expects_curve) in [
        ("density", "handbook", false),
        ("magnetic_flux_density", "bh-curve", true),
        ("fatigue_life", "sn-curve", true),
        ("dynamic_viscosity", "lubricant", true),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique {property} claim");
        let claim = claims[0].1;
        assert_eq!(
            matches!(&claim.value, PropertyValue::Curve { .. }),
            expects_curve,
            "unexpected payload kind for {property}"
        );
        assert!(
            claim
                .provenance
                .source
                .contains(&format!("[source:{source}]")),
            "{property} lost source-local provenance: {:?}",
            claim.provenance
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    for source in ["handbook", "bh-curve", "sn-curve", "lubricant"] {
        assert!(
            decisions.contains(&format!("\"subject\":\"source:{source}\"")),
            "missing admission row for {source}"
        );
    }
}

#[test]
fn g3_cli_compiles_nasa9_regions_into_identical_verified_model_packs() {
    let (directory, manifest) = write_nasa9_fixture();
    let first_path = directory.join("nasa9-first.fsmodpk");
    let second_path = directory.join("nasa9-second.fsmodpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-9 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-9 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "NASA-9 decision stream moved");

    let first_bytes = fs::read(first_path).expect("read first NASA-9 pack");
    let second_bytes = fs::read(second_path).expect("read second NASA-9 pack");
    assert_eq!(first_bytes, second_bytes, "NASA-9 pack bytes moved");
    assert_eq!(first_bytes.len(), NASA9_PACK_BYTES_GOLDEN);
    let decoded = NormalizedModelPack::from_bytes(&first_bytes).expect("decode NASA-9 model pack");
    assert_eq!(decoded.content_hash().to_string(), NASA9_PACK_HASH_GOLDEN);
    let decoded = NormalizedModelPack::from_bytes_verified(decoded.content_hash(), &first_bytes)
        .expect("verified NASA-9 model pack");
    assert_eq!(decoded.pack_id(), "N2");
    assert_eq!(
        decoded.compiler(),
        "frankensim-matdb-nasa9-model-pack-compiler-v1"
    );
    assert_eq!(decoded.models().len(), 2);
    assert_eq!(decoded.normalizations().len(), 24);
    assert!(decoded.models().iter().all(|card| {
        card.law.0 == "nasa9-standard-state"
            && card.law_version == 1
            && card.parameters.len() == 10
            && card.validity.bound("T").is_some()
            && card.provenance.source.contains("[species:N2]")
    }));

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"species_pack_id_bound\""));
    assert!(decisions.contains("\"reason_code\":\"nasa9_region_normalized\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_model_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{NASA9_PACK_HASH_GOLDEN}\"")))
    );
}

#[test]
fn g3_cli_refuses_standalone_species_and_kinetics_without_runtime_codecs() {
    for profile in ["species-v1", "kinetics-v1"] {
        let manifest_text = MANIFEST.replacen("material-tsv-v1", profile, 1);
        assert_ne!(manifest_text, MANIFEST);
        let directory = fixture_dir();
        let manifest = directory.join("manifest.tsv");
        let output = directory.join(format!("{profile}.fsmatpk"));
        fs::write(&manifest, manifest_text).expect("write unsupported-profile manifest");
        fs::write(directory.join("source.tsv"), SOURCE).expect("write source fixture");

        let refused = run_compiler(&manifest, &output);
        assert!(!refused.status.success(), "{profile} unexpectedly compiled");
        assert!(!output.exists(), "{profile} refusal published an output");
        let decisions = String::from_utf8(refused.stdout).expect("decision stream is UTF-8");
        assert_eq!(decisions.matches("\"verdict\":\"refuse\"").count(), 1);
        assert!(decisions.contains("\"reason_code\":\"unsupported_source_profile\""));
        assert!(decisions.contains("\"subject\":\"source:primary\""));
        assert!(
            String::from_utf8_lossy(&refused.stderr)
                .contains("error: matdb pack refused [unsupported_source_profile]")
        );
    }
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
