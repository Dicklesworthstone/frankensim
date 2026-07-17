#![deny(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    InterpolationPolicy, NormalizedModelPack, NormalizedPack, NormalizedSpeciesPack, PropertyValue,
    SPECIES_MOLAR_MASS_DIMS, SPECIES_PACK_TARGET_BASIS, SPECIES_REFERENCE_PRESSURE_DIMS,
    SpeciesNormalizationTarget, UncertaintyModel,
};
use fs_qty::Dims;

static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

const PACK_BYTES_GOLDEN: usize = 3_177;
const PACK_HASH_GOLDEN: &str = "c1fb2f443708d297423179f4ac6024ee26b1d0c940a229d1d9084726ccbd2bc5";
const NASA9_PACK_BYTES_GOLDEN: usize = 4_940;
const NASA9_PACK_HASH_GOLDEN: &str =
    "006177a7cc6f7b4ae10a9eb4a5bf49faaf21911ef9473190a29ecfc3a818a162";
const MATERIAL_COMPILER_ID: &str = "frankensim-matdb-pack-compiler-v1";
const NASA9_COMPILER_ID: &str = "frankensim-matdb-nasa9-model-pack-compiler-v1";
const KINETICS_COMPILER_ID: &str = "frankensim-matdb-kinetics-model-pack-compiler-v1";
const SPECIES_COMPILER_ID: &str = "frankensim-matdb-species-pack-compiler-v1";
const METHANE_SEED_MANIFEST: &str = "data/matdb/seed-v1/methane/manifest.tsv";
const ALUMINUM_6061_T6_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aluminum-6061-t6-cryogenic/manifest.tsv";
const OFHC_COPPER_SEED_MANIFEST: &str = "data/matdb/seed-v1/ofhc-copper-rrr100/manifest.tsv";
const AISI_4140_RC33_SEED_MANIFEST: &str = "data/matdb/seed-v1/aisi-4140-rc33/manifest.tsv";
const AISI_1045_COLD_DRAWN_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-1045-cold-drawn/manifest.tsv";
const AISI_52100_CVM_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-52100-cvm-hot-hardness/manifest.tsv";
const GRAY_CAST_IRON_S2_S_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/gray-cast-iron-s2-s/manifest.tsv";
const NASA_SEED_LICENSE: &str = "Work-of-the-US-Government-Public-Use-Permitted";
const CC_BY_4_0_LICENSE: &str = "CC-BY-4.0";
const NIST_PUBLIC_INFORMATION_LICENSE: &str = "NIST-Public-Information-Attribution-Requested";
const NASA_METHANE_MOLAR_MASS_G_PER_MOL: f64 = 16.042_46;
const NIST_SRD69_METHANE_MOLAR_MASS_KG_PER_MOL: f64 = 0.016_042_5;
const NIST_SRD69_DISPLAY_ROUNDING_HALF_WIDTH_KG_PER_MOL: f64 = 0.000_000_05;

#[derive(Clone, Copy)]
struct CommittedSpeciesSeed {
    manifest: &'static str,
    species: &'static str,
    nasa_molar_mass_g_per_mol: f64,
    nist_molar_mass_g_per_mol: f64,
    nist_display_rounding_half_width_g_per_mol: f64,
}

const AIR_EXHAUST_SPECIES_SEEDS: [CommittedSpeciesSeed; 6] = [
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/nitrogen/manifest.tsv",
        species: "N2",
        nasa_molar_mass_g_per_mol: 28.013_40,
        nist_molar_mass_g_per_mol: 28.013_4,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/oxygen/manifest.tsv",
        species: "O2",
        nasa_molar_mass_g_per_mol: 31.998_80,
        nist_molar_mass_g_per_mol: 31.998_8,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/argon/manifest.tsv",
        species: "Ar",
        nasa_molar_mass_g_per_mol: 39.948_00,
        nist_molar_mass_g_per_mol: 39.948,
        nist_display_rounding_half_width_g_per_mol: 0.000_5,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/carbon-dioxide/manifest.tsv",
        species: "CO2",
        nasa_molar_mass_g_per_mol: 44.009_50,
        nist_molar_mass_g_per_mol: 44.009_5,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/water-vapor/manifest.tsv",
        species: "H2O",
        nasa_molar_mass_g_per_mol: 18.015_28,
        nist_molar_mass_g_per_mol: 18.015_3,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
    CommittedSpeciesSeed {
        manifest: "data/matdb/seed-v1/carbon-monoxide/manifest.tsv",
        species: "CO",
        nasa_molar_mass_g_per_mol: 28.010_10,
        nist_molar_mass_g_per_mol: 28.010_1,
        nist_display_rounding_half_width_g_per_mol: 0.000_05,
    },
];

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

const KINETICS_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\twater-formation\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture first-order kinetics table\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tkinetics.tsv\tkinetics-v1\n",
);

const KINETICS_SOURCE: &str = concat!(
    "frankensim.kinetics-source.v1\n",
    "reaction\twater-formation\tfirst-order\t300\t2500\tK\n",
    "parameter\twater-formation\tactivation_temperature\t12000\tK\n",
    "parameter\twater-formation\tpre_exponential\t2.5e7\ts^-1\n",
);

const SPECIES_MANIFEST: &str = concat!(
    "frankensim.matdb-manifest.v1\n",
    "pack_id\tN2\n",
    "redistribution\tpermitted\tCC-BY-4.0 redistribution with attribution\n",
    "citation\tfixture licensed species metadata\n",
    "license\tCC-BY-4.0\n",
    "source\tprimary\tspecies.tsv\tspecies-v1\n",
);

const SPECIES_SOURCE: &str = concat!(
    "frankensim.species-source.v1\n",
    "species\tN2\t28.0134\tg/mol\tgas\tideal-gas\t100\tkPa\tNASA-TP-2002-211556\n",
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

fn write_kinetics_fixture() -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, KINETICS_MANIFEST).expect("write kinetics manifest fixture");
    fs::write(directory.join("kinetics.tsv"), KINETICS_SOURCE)
        .expect("write kinetics source fixture");
    (directory, manifest)
}

fn write_species_fixture(source: &str) -> (PathBuf, PathBuf) {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    fs::write(&manifest, SPECIES_MANIFEST).expect("write species manifest fixture");
    fs::write(directory.join("species.tsv"), source).expect("write species source fixture");
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

fn workspace_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a workspace parent")
        .join(relative)
}

fn assert_decision_compiler(output: &Output, expected: &str) {
    let stdout = std::str::from_utf8(&output.stdout).expect("decision stream is UTF-8");
    assert!(!stdout.is_empty(), "compiler emitted no decision rows");
    let expected_prefix = format!("{{\"check\":\"matdb-pack\",\"compiler\":\"{expected}\",");
    assert!(
        stdout.lines().all(|row| row.starts_with(&expected_prefix)),
        "decision row used the wrong compiler identity:\n{stdout}"
    );
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
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

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
fn g3_cli_compiles_committed_aluminum_6061_t6_exact_point_seed() {
    let manifest = workspace_path(ALUMINUM_6061_T6_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Aluminum 6061-T6 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aluminum-6061-t6-first.fsmatpk");
    let second_path = directory.join("aluminum-6061-t6-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Aluminum 6061-T6 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Aluminum 6061-T6 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Aluminum 6061-T6 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Aluminum 6061-T6 pack");
    let second_bytes = fs::read(second_path).expect("read second Aluminum 6061-T6 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Aluminum 6061-T6 pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode Aluminum 6061-T6 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Aluminum 6061-T6 pack identity");

    assert_eq!(decoded.pack_id(), "aluminum-6061-t6-cryogenic");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public information")
    );
    assert_eq!(decoded.claims().claim_count(), 6);
    assert!(decoded.joint_statistics().is_empty());

    let expected = [
        (
            "thermal_conductivity",
            77.0,
            83.531_441_947_072_3,
            Dims([1, 1, -3, -1, 0, 0]),
            "0.5 percent curve-fit error",
            "a..i=0.07918,1.0957",
        ),
        (
            "thermal_conductivity",
            293.0,
            154.345_205_650_720,
            Dims([1, 1, -3, -1, 0, 0]),
            "0.5 percent curve-fit error",
            "a..i=0.07918,1.0957",
        ),
        (
            "specific_heat_capacity",
            77.0,
            348.127_924_911_502,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=46.6467,-314.292",
        ),
        (
            "specific_heat_capacity",
            293.0,
            942.911_235_990_911,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=46.6467,-314.292",
        ),
        (
            "young_modulus",
            77.0,
            77.145_050_657_273_1e9,
            Dims([-1, 1, -2, 0, 0, 0]),
            "1 percent curve-fit error",
            "a..e=77.71221,0.01030646",
        ),
        (
            "young_modulus",
            293.0,
            70.358_592_182_729_1e9,
            Dims([-1, 1, -2, 0, 0, 0]),
            "1 percent curve-fit error",
            "a..e=77.71221,0.01030646",
        ),
    ];

    for (property, temperature, expected_value, expected_dims, fit_error_note, coefficient_note) in
        expected
    {
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing {property} claim at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("{property} at {temperature} K was not an exact-point scalar");
        };
        assert_eq!(*dims, expected_dims, "{property} dimensions moved");
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "{property} at {temperature} K moved by {relative_error:e} relative"
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NIST_PUBLIC_INFORMATION_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("National Institute of Standards and Technology")
        );
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nist-cryogenic-fit]")
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "aluminum-6061-t6-uns-aa96061-temper-t6"
        );
        assert!(observation.method.contains("NIST"));
        assert!(
            observation
                .method
                .contains("exact-temperature derived scalars")
        );
        assert!(observation.caveats.contains(fit_error_note));
        assert!(observation.caveats.contains(coefficient_note));
        assert!(
            observation
                .caveats
                .contains("without a confidence level or degrees of freedom")
        );
    }

    // G3 independent-source evidence only: NASA's 1966 compilation reports
    // 82 W/(m K) at 75 K and 155 W/(m K) at 300 K. These nearby-temperature
    // checks do not replace the exact NIST-derived values stored above.
    for (nist_temperature, nasa_temperature, nasa_value) in
        [(77.0, 75.0, 82.0), (293.0, 300.0, 155.0)]
    {
        let (_, claim) = decoded
            .claims()
            .claims_for("thermal_conductivity")
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((nist_temperature, nist_temperature))
            })
            .expect("thermal-conductivity comparison point");
        let PropertyValue::Scalar { value, .. } = &claim.value else {
            panic!("thermal-conductivity comparison point is scalar");
        };
        assert!((nist_temperature - nasa_temperature).abs() <= 7.0);
        let relative_difference = (*value - nasa_value).abs() / nasa_value;
        assert!(
            relative_difference <= 0.03,
            "NIST-derived {nist_temperature} K conductivity and NASA {nasa_temperature} K comparison differ by {relative_difference:e}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_ofhc_copper_exact_point_seed() {
    let manifest = workspace_path(OFHC_COPPER_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed OFHC Copper seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("ofhc-copper-first.fsmatpk");
    let second_path = directory.join("ofhc-copper-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first OFHC Copper seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second OFHC Copper seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "OFHC Copper decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first OFHC Copper pack");
    let second_bytes = fs::read(second_path).expect("read second OFHC Copper pack");
    assert_eq!(first_bytes, second_bytes, "OFHC Copper pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode OFHC Copper pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify OFHC Copper pack identity");

    assert_eq!(decoded.pack_id(), "ofhc-copper-cryogenic");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public information")
    );
    assert_eq!(decoded.claims().claim_count(), 4);
    assert!(decoded.joint_statistics().is_empty());

    let expected = [
        (
            "thermal_conductivity",
            77.0,
            547.199_698_079_367,
            Dims([1, 1, -3, -1, 0, 0]),
            "1 percent curve-fit error",
            "a..i=2.2154,-0.47461",
            "ofhc-copper-uns-c10100-c10200-rrr100",
        ),
        (
            "thermal_conductivity",
            293.0,
            396.908_547_137_121,
            Dims([1, 1, -3, -1, 0, 0]),
            "1 percent curve-fit error",
            "a..i=2.2154,-0.47461",
            "ofhc-copper-uns-c10100-c10200-rrr100",
        ),
        (
            "specific_heat_capacity",
            77.0,
            195.920_875_203_329,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=-1.91844,-0.15973",
            "ofhc-copper-uns-c10100-c10200-source-rrr-unspecified",
        ),
        (
            "specific_heat_capacity",
            293.0,
            389.085_653_150_371,
            Dims([2, 0, -2, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=-1.91844,-0.15973",
            "ofhc-copper-uns-c10100-c10200-source-rrr-unspecified",
        ),
    ];

    for (
        property,
        temperature,
        expected_value,
        expected_dims,
        fit_error_note,
        coefficient_note,
        expected_specimen,
    ) in expected
    {
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing OFHC {property} claim at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("OFHC {property} at {temperature} K was not an exact-point scalar");
        };
        assert_eq!(*dims, expected_dims, "OFHC {property} dimensions moved");
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "OFHC {property} at {temperature} K moved by {relative_error:e} relative"
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NIST_PUBLIC_INFORMATION_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("Material Properties: OFHC Copper")
        );
        assert!(claim.provenance.source.contains("[source:nist-ofhc-fit]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("OFHC claim observation remains linked");
        assert_eq!(observation.specimen, expected_specimen);
        assert!(observation.method.contains("NIST"));
        assert!(
            observation
                .method
                .contains("exact-temperature derived scalars")
        );
        assert!(observation.caveats.contains(fit_error_note));
        assert!(observation.caveats.contains(coefficient_note));
        assert!(
            observation
                .caveats
                .contains("without a confidence level or degrees of freedom")
        );
    }

    // G3 independent-source evidence only: NASA-CR-134806 reports typical
    // room-temperature OFHC Copper values of 390 W/(m K) and 386 J/(kg K).
    // They remain comparisons and do not replace the NIST-derived claims.
    for (property, nasa_value) in [
        ("thermal_conductivity", 390.0),
        ("specific_heat_capacity", 386.0),
    ] {
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| claim.validity.bound("temperature") == Some((293.0, 293.0)))
            .unwrap_or_else(|| panic!("missing OFHC room-temperature {property}"));
        let PropertyValue::Scalar { value, .. } = &claim.value else {
            panic!("OFHC room-temperature comparison point was not scalar");
        };
        let relative_difference = (*value - nasa_value).abs() / nasa_value;
        assert!(
            relative_difference <= 0.02,
            "NIST-derived 293 K OFHC {property} and NASA room-temperature comparison differ by {relative_difference:e}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_aisi_4140_rc33_exact_condition_seed() {
    let manifest = workspace_path(AISI_4140_RC33_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 4140 Rockwell C33 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-4140-rc33-first.fsmatpk");
    let second_path = directory.join("aisi-4140-rc33-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 4140 Rockwell C33 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 4140 Rockwell C33 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 4140 Rockwell C33 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 4140 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 4140 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 4140 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 4140 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 4140 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-4140-qq-s-624-rc33");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 14);
    assert!(decoded.joint_statistics().is_empty());

    let pressure_dims = Dims([-1, 1, -2, 0, 0, 0]);
    let energy_dims = Dims([2, 1, -2, 0, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let expected = [
        (
            "ultimate_tensile_strength",
            26.7,
            1.074 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "yield_strength_0p2_offset",
            26.7,
            0.985 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_elongation_2in",
            26.7,
            19.4 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_reduction_of_area",
            26.7,
            62.5 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "charpy_v_notch_impact_energy",
            26.7,
            95.2,
            energy_dims,
            "MIL-STD-151 Charpy V-notched impact",
            "four impact tests",
        ),
        (
            "double_shear_ultimate_strength",
            26.7,
            0.66 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
        (
            "double_shear_yield_strength",
            26.7,
            0.56 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
        (
            "ultimate_tensile_strength",
            -73.0,
            1.158 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "yield_strength_0p2_offset",
            -73.0,
            1.060 * 1.0e9,
            pressure_dims,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_elongation_2in",
            -73.0,
            20.0 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "tensile_reduction_of_area",
            -73.0,
            61.0 * 0.01,
            dimensionless,
            "longitudinal round smooth tensile",
            "five smooth and five notched tensile specimens",
        ),
        (
            "charpy_v_notch_impact_energy",
            -73.0,
            84.6,
            energy_dims,
            "MIL-STD-151 Charpy V-notched impact",
            "four impact tests",
        ),
        (
            "double_shear_ultimate_strength",
            -73.0,
            0.73 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
        (
            "double_shear_yield_strength",
            -73.0,
            0.60 * 1.0e9,
            pressure_dims,
            "double-shear specimens",
            "four shear specimens",
        ),
    ];

    for (property, temperature_c, expected_value, expected_dims, method_note, sample_note) in
        expected
    {
        let expected_value: f64 = expected_value;
        let temperature_k = temperature_c + 273.15;
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
            })
            .unwrap_or_else(|| {
                panic!("missing AISI 4140 {property} claim at {temperature_c} degC")
            });
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 4140 {property} at {temperature_c} degC was not scalar");
        };
        assert_eq!(
            *dims, expected_dims,
            "AISI 4140 {property} dimensions moved"
        );
        let scale = f64::abs(expected_value).max(1.0);
        let relative_error = (*value - expected_value).abs() / scale;
        assert!(
            relative_error <= 2.0e-15,
            "AISI 4140 {property} at {temperature_c} degC moved by {relative_error:e} relative"
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TM-X-64791"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 4140 claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-4140-QQ-S-624-heat-137M186-1in-bar-Rockwell-C33"
        );
        assert!(observation.method.contains(method_note));
        assert!(observation.caveats.contains(sample_note));
        assert!(observation.caveats.contains("oil quenched"));
        assert!(observation.caveats.contains("tempered 566 degC"));
    }

    // NASA Table IV prints both ksi and GN/m2. This checks transcription and
    // unit normalization against the redundant source columns; it is not an
    // independent-source agreement claim.
    for (property, temperature_c, source_ksi) in [
        ("ultimate_tensile_strength", 26.7, 155.8),
        ("yield_strength_0p2_offset", 26.7, 142.9),
        ("ultimate_tensile_strength", -73.0, 168.0),
        ("yield_strength_0p2_offset", -73.0, 153.7),
    ] {
        let temperature_k = temperature_c + 273.15;
        let (_, claim) = decoded
            .claims()
            .claims_for(property)
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
            })
            .expect("AISI 4140 redundant-unit comparison point");
        let PropertyValue::Scalar { value, .. } = &claim.value else {
            panic!("AISI 4140 redundant-unit comparison point was not scalar");
        };
        let source_ksi_in_pa = source_ksi * 6_894_757.293_168_361;
        let relative_rounding_difference = (*value - source_ksi_in_pa).abs() / *value;
        assert!(
            relative_rounding_difference <= 5.0e-4,
            "AISI 4140 {property} at {temperature_c} degC disagrees with the source ksi column by {relative_rounding_difference:e}"
        );
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_aisi_1045_cold_drawn_tensile_seed() {
    let manifest = workspace_path(AISI_1045_COLD_DRAWN_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 1045 cold-drawn seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-1045-cold-drawn-first.fsmatpk");
    let second_path = directory.join("aisi-1045-cold-drawn-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 1045 cold-drawn seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 1045 cold-drawn seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 1045 cold-drawn decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 1045 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 1045 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 1045 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 1045 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 1045 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-1045-cold-drawn-tensile");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("Attribution 4.0 International")
    );
    assert_eq!(decoded.claims().claim_count(), 3);
    assert!(
        decoded.joint_statistics().is_empty(),
        "paired source rows do not authorize an inferred covariance block"
    );

    let pressure_dims = Dims([-1, 1, -2, 0, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let expected = [
        (
            "yield_strength",
            550.51,
            0.005,
            [540.73, 557.59, 553.20],
            1.0e6,
            pressure_dims,
        ),
        (
            "ultimate_tensile_strength",
            695.31,
            0.005,
            [684.58, 707.75, 693.60],
            1.0e6,
            pressure_dims,
        ),
        (
            "tensile_elongation_50mm",
            14.1,
            0.05,
            [14.42, 14.20, 13.68],
            0.01,
            dimensionless,
        ),
    ];
    let student_t_0p975_df2 = 4.302_652_729_911_275;
    let crosshead_speed_m_per_s = 10.0 * 1.0e-3 / 60.0;

    for (
        property,
        reported_mean,
        reported_rounding_half_width,
        samples,
        source_unit_scale,
        expected_dims,
    ) in expected
    {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(
            claims.len(),
            1,
            "expected exactly one AISI 1045 {property} claim"
        );
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 1045 {property} was not scalar");
        };
        assert_eq!(
            *dims, expected_dims,
            "AISI 1045 {property} dimensions moved"
        );
        let expected_value = reported_mean * source_unit_scale;
        let relative_value_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
        assert!(
            relative_value_error <= 2.0e-15,
            "AISI 1045 {property} moved by {relative_value_error:e} relative"
        );

        let sample_mean = samples.iter().copied().sum::<f64>() / 3.0;
        assert!(
            (sample_mean - reported_mean).abs() <= reported_rounding_half_width,
            "AISI 1045 {property} source mean is inconsistent with its printed replicates"
        );
        let sample_variance = samples
            .iter()
            .map(|sample| (sample - sample_mean).powi(2))
            .sum::<f64>()
            / 2.0;
        let expected_half_width =
            student_t_0p975_df2 * sample_variance.sqrt() / 3.0_f64.sqrt() * source_unit_scale;
        let UncertaintyModel::HalfWidth {
            half_width,
            confidence,
        } = &claim.uncertainty
        else {
            panic!("AISI 1045 {property} lost its derived Student-t half-width");
        };
        let relative_half_width_error =
            (*half_width - expected_half_width).abs() / expected_half_width;
        assert!(
            relative_half_width_error <= 2.0e-14,
            "AISI 1045 {property} Student-t half-width moved by {relative_half_width_error:e} relative"
        );
        assert_eq!(*confidence, 0.95);

        let Some((speed_lo, speed_hi)) = claim.validity.bound("crosshead_speed") else {
            panic!("AISI 1045 {property} lost its crosshead-speed validity point");
        };
        for speed in [speed_lo, speed_hi] {
            let relative_speed_error =
                (speed - crosshead_speed_m_per_s).abs() / crosshead_speed_m_per_s;
            assert!(
                relative_speed_error <= 2.0e-15,
                "AISI 1045 {property} crosshead speed moved by {relative_speed_error:e} relative"
            );
        }
        assert_eq!(
            claim.validity.bound("source_test_temperature_known"),
            Some((0.0, 0.0)),
            "AISI 1045 {property} must require explicit acknowledgement of missing temperature"
        );
        assert_eq!(claim.validity.bounds().len(), 2);
        assert_eq!(claim.validity.bound("temperature"), None);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(claim.provenance.source.contains("doi:10.3390/pr12061171"));
        assert!(claim.provenance.source.contains("[source:primary]"));

        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 1045 claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-1045-cold-drawn-bar-37mm-OD-102mm-length-test-temperature-not-reported"
        );
        assert!(observation.method.contains("ASTM E8"));
        assert!(observation.method.contains("50 mm gauge length"));
        assert!(observation.method.contains("10 mm/min crosshead speed"));
        assert!(observation.caveats.contains("three specimens"));
        assert!(observation.caveats.contains("t(0.975, df=2)"));
        assert!(
            observation
                .caveats
                .contains("source does not report test temperature")
        );
        assert!(
            observation
                .caveats
                .contains("no joint covariance is inferred")
        );
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_aisi_52100_cvm_heat_treatment_states() {
    let manifest = workspace_path(AISI_52100_CVM_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 52100 CVM seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-52100-cvm-first.fsmatpk");
    let second_path = directory.join("aisi-52100-cvm-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 52100 CVM seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 52100 CVM seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 52100 CVM decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 52100 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 52100 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 52100 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 52100 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 52100 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-52100-cvm-nasa-tn-d-6632");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 15);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    for (property, source_percent) in [
        ("carbon_mass_fraction", 0.96),
        ("silicon_mass_fraction", 0.22),
        ("manganese_mass_fraction", 0.36),
        ("sulfur_mass_fraction", 0.012),
        ("phosphorus_mass_fraction", 0.007),
        ("chromium_mass_fraction", 1.36),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one AISI 52100 {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 52100 {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "AISI 52100 {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TN-D-6632"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 52100 chemistry observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-52100-consumable-vacuum-melted-single-ingot-NASA-TN-D-6632"
        );
        assert!(observation.method.contains("Table I actual composition"));
        assert!(observation.caveats.contains("Balance iron"));
        assert!(observation.caveats.contains("no heat identifier"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let hardness_claims = decoded.claims().claims_for("rockwell_c_scale_reading");
    let austenite_claims = decoded
        .claims()
        .claims_for("retained_austenite_volume_fraction");
    assert_eq!(hardness_claims.len(), 5);
    assert_eq!(austenite_claims.len(), 4);
    let states = [
        (Some(505.0), "second-temper-505K", 59.7, None),
        (Some(450.0), "second-temper-450K", 62.3, Some(12.8)),
        (Some(433.0), "second-temper-433K", 63.4, Some(15.6)),
        (Some(394.0), "second-temper-394K", 64.6, Some(18.4)),
        (None, "no-second-temper", 65.1, Some(11.8)),
    ];

    for (second_temper_k, specimen_state, expected_hardness, expected_austenite_percent) in states {
        let matches_state = |claim: &fs_matdb::PropertyClaim| match second_temper_k {
            Some(second_temper_k) => {
                claim.validity.bound("second_temper_temperature")
                    == Some((second_temper_k, second_temper_k))
                    && claim.validity.bound("second_temper_applied").is_none()
            }
            None => {
                claim.validity.bound("second_temper_applied") == Some((0.0, 0.0))
                    && claim.validity.bound("second_temper_temperature").is_none()
            }
        };
        let (hardness_id, hardness_claim) = hardness_claims
            .iter()
            .copied()
            .find(|(_, claim)| matches_state(claim))
            .unwrap_or_else(|| panic!("missing AISI 52100 hardness state {specimen_state}"));
        let PropertyValue::Scalar { value, dims } = &hardness_claim.value else {
            panic!("AISI 52100 hardness state {specimen_state} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        assert_eq!(*value, expected_hardness);
        assert_eq!(
            hardness_claim.validity.bound("temperature"),
            Some((294.0, 294.0))
        );
        assert_eq!(hardness_claim.validity.bounds().len(), 2);
        assert_eq!(hardness_claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            hardness_claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(hardness_claim.observations.len(), 1);
        assert_eq!(hardness_claim.provenance.license, NASA_SEED_LICENSE);
        let hardness_observation = decoded
            .claims()
            .observation(hardness_claim.observations[0])
            .expect("AISI 52100 hardness observation remains linked");
        assert!(hardness_observation.specimen.contains(specimen_state));
        assert!(
            hardness_observation
                .specimen
                .contains("austenitize-1116-to-1144K-30min")
        );
        assert!(hardness_observation.specimen.contains("oil-quench-325K"));
        assert!(
            hardness_observation
                .specimen
                .contains("first-temper-394K-60min")
        );
        assert!(hardness_observation.method.contains("150 kg load"));
        assert!(hardness_observation.method.contains("Rockwell C"));
        assert!(hardness_observation.method.contains("294 K reading"));
        assert!(
            hardness_observation
                .caveats
                .contains("Minimum two hardness measurements")
        );
        assert!(
            hardness_observation
                .caveats
                .contains("dispersion not reported")
        );
        assert!(hardness_observation.caveats.contains("ASTM grain size 12"));
        assert!(
            hardness_observation
                .caveats
                .contains("predictive equation is not measurement uncertainty")
        );
        if second_temper_k == Some(505.0) {
            assert!(
                hardness_observation
                    .caveats
                    .contains("censored as less than 2 volume percent")
            );
        }
        assert_eq!(
            hardness_claim.observations[0].0,
            hardness_observation.content_hash()
        );
        assert_eq!(hardness_id.0, hardness_claim.content_hash());

        let matching_austenite = austenite_claims
            .iter()
            .copied()
            .find(|(_, claim)| matches_state(claim));
        match (matching_austenite, expected_austenite_percent) {
            (Some((austenite_id, austenite_claim)), Some(expected_percent)) => {
                let PropertyValue::Scalar { value, dims } = &austenite_claim.value else {
                    panic!("AISI 52100 retained-austenite state {specimen_state} was not scalar");
                };
                assert_eq!(*dims, dimensionless);
                let expected_value = expected_percent * 0.01;
                let relative_error = (*value - expected_value).abs() / expected_value;
                assert!(
                    relative_error <= 2.0e-15,
                    "AISI 52100 retained austenite {specimen_state} moved by {relative_error:e} relative"
                );
                assert_eq!(
                    austenite_claim.validity.bound("temperature"),
                    Some((294.0, 294.0))
                );
                assert_eq!(austenite_claim.validity.bounds().len(), 2);
                assert_eq!(austenite_claim.uncertainty, UncertaintyModel::Unstated);
                assert_eq!(austenite_claim.observations.len(), 1);
                let observation = decoded
                    .claims()
                    .observation(austenite_claim.observations[0])
                    .expect("AISI 52100 retained-austenite observation remains linked");
                assert!(observation.specimen.contains(specimen_state));
                assert!(observation.method.contains("X-ray diffraction"));
                assert!(
                    observation
                        .caveats
                        .contains("uncertainty and replicate count not reported")
                );
                assert!(
                    observation
                        .caveats
                        .contains("no covariance with hardness is inferred")
                );
                assert_ne!(
                    austenite_claim.observations[0],
                    hardness_claim.observations[0]
                );
                assert_eq!(
                    austenite_claim.observations[0].0,
                    observation.content_hash()
                );
                assert_eq!(austenite_id.0, austenite_claim.content_hash());
            }
            (None, None) => {}
            (Some(_), None) => panic!("censored AISI 52100 austenite became an exact scalar"),
            (None, Some(_)) => panic!("missing exact AISI 52100 austenite state {specimen_state}"),
        }
    }

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_s2_s_gray_cast_iron_seed() {
    let manifest = workspace_path(GRAY_CAST_IRON_S2_S_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed S2-S gray-cast-iron seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("gray-cast-iron-s2-s-first.fsmatpk");
    let second_path = directory.join("gray-cast-iron-s2-s-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first S2-S gray-cast-iron seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second S2-S gray-cast-iron seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "S2-S gray-cast-iron decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first S2-S gray-iron pack");
    let second_bytes = fs::read(second_path).expect("read second S2-S gray-iron pack");
    assert_eq!(first_bytes, second_bytes, "S2-S gray-iron pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode S2-S gray-iron pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify S2-S gray-iron pack identity");

    assert_eq!(decoded.pack_id(), "pearlitic-gray-cast-iron-s2-s-sr-fesi");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("Attribution 4.0 International")
    );
    assert_eq!(decoded.claims().claim_count(), 15);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    for (property, source_percent) in [
        ("carbon_mass_fraction", 3.54),
        ("silicon_mass_fraction", 1.62),
        ("manganese_mass_fraction", 0.51),
        ("phosphorus_mass_fraction", 0.025),
        ("sulfur_mass_fraction", 0.028),
        ("molybdenum_mass_fraction", 0.35),
        ("copper_mass_fraction", 0.58),
        ("tin_mass_fraction", 0.060),
        ("carbon_equivalent_ce", 4.05),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one S2-S {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("S2-S {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "S2-S {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(claim.provenance.source.contains("doi:10.3390/ma11101876"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("S2-S composition observation remains linked");
        assert!(
            observation
                .specimen
                .contains("S2-S-pearlitic-gray-cast-iron")
        );
        assert!(observation.specimen.contains("0p4wtpct-SrFeSi-Ino2"));
        assert!(observation.method.contains("Table 1"));
        assert!(observation.caveats.contains("2.0 wt% Sr"));
        assert!(observation.caveats.contains("no balance-iron scalar"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let carbon_equivalent_from_printed_composition: f64 = 3.54 + 0.31 * 1.62 + 0.33 * 0.025;
    assert!(
        (carbon_equivalent_from_printed_composition - 4.05).abs() <= 0.005,
        "S2-S carbon-equivalent transcription exceeds the source's printed rounding"
    );

    for (property, expected_value, expected_dims, caveat_fragment) in [
        (
            "graphite_area_fraction",
            9.0 * 0.01,
            dimensionless,
            "graphite area 9.0 +/- 0.2 percent",
        ),
        (
            "maximum_graphite_flake_length",
            273.0 * 1.0e-6,
            Dims([1, 0, 0, 0, 0, 0]),
            "maximum graphite length 273 +/- 19 um",
        ),
        (
            "primary_dendrite_area_fraction",
            15.6 * 0.01,
            dimensionless,
            "primary-dendrite area 15.6 +/- 0.9 percent",
        ),
        (
            "eutectic_colony_areal_density",
            371.0 * 1.0e4,
            Dims([-2, 0, 0, 0, 0, 0]),
            "eutectic-colony count 371 +/- 19 per cm2",
        ),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one S2-S {property} claim");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("S2-S {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "S2-S {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("S2-S microstructure observation remains linked");
        assert!(observation.method.contains("eight cross-section fields"));
        assert!(observation.caveats.contains("type-A graphite"));
        assert!(observation.caveats.contains(caveat_fragment));
        assert!(observation.caveats.contains("one standard deviation"));
        assert!(
            observation
                .caveats
                .contains("runtime uncertainty remains Unstated")
        );
    }

    let (_, tensile) = decoded
        .claims()
        .claims_for("ultimate_tensile_strength")
        .into_iter()
        .next()
        .expect("S2-S ultimate tensile strength claim");
    let PropertyValue::Scalar {
        value: tensile_value,
        dims: tensile_dims,
    } = &tensile.value
    else {
        panic!("S2-S ultimate tensile strength was not scalar");
    };
    assert_eq!(*tensile_dims, Dims([-1, 1, -2, 0, 0, 0]));
    assert_eq!(*tensile_value, 326.0e6);
    assert_eq!(
        tensile.validity.bound("source_test_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(tensile.validity.bounds().len(), 1);
    assert_eq!(tensile.uncertainty, UncertaintyModel::Unstated);
    let tensile_observation = decoded
        .claims()
        .observation(tensile.observations[0])
        .expect("S2-S tensile observation remains linked");
    assert!(tensile_observation.method.contains("GB/T T228.1-2010"));
    assert!(tensile_observation.method.contains("three tests averaged"));
    assert!(tensile_observation.caveats.contains("nearest 1 MPa"));
    assert!(tensile_observation.caveats.contains("approximately 8 MPa"));
    assert!(
        tensile_observation
            .caveats
            .contains("exact test temperature is not reported")
    );

    let (_, conductivity) = decoded
        .claims()
        .claims_for("thermal_conductivity")
        .into_iter()
        .next()
        .expect("S2-S thermal conductivity claim");
    let PropertyValue::Scalar {
        value: conductivity_value,
        dims: conductivity_dims,
    } = &conductivity.value
    else {
        panic!("S2-S thermal conductivity was not scalar");
    };
    assert_eq!(*conductivity_dims, Dims([1, 1, -3, -1, 0, 0]));
    assert_eq!(*conductivity_value, 58.8);
    assert_eq!(
        conductivity.validity.bound("source_test_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(conductivity.validity.bounds().len(), 1);
    assert_eq!(conductivity.uncertainty, UncertaintyModel::Unstated);
    let conductivity_observation = decoded
        .claims()
        .observation(conductivity.observations[0])
        .expect("S2-S thermal observation remains linked");
    assert!(conductivity_observation.method.contains("NETZSCH LFA 457"));
    assert!(
        conductivity_observation
            .method
            .contains("Archimedes density")
    );
    assert!(
        conductivity_observation
            .caveats
            .contains("nearest 0.1 W/(m K)")
    );
    assert!(
        conductivity_observation
            .caveats
            .contains("approximately 0.3 W/(m K)")
    );
    assert!(
        conductivity_observation
            .caveats
            .contains("exact room temperature is not reported")
    );

    // G3 independent-source plausibility evidence only: ORNL/TM-2012/506
    // Appendix C gives a broad 42..62 W/(m K) range for generic gray cast
    // iron. It neither identifies S2-S nor overwrites the primary claim.
    assert!((42.0..=62.0).contains(conductivity_value));

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"uncertainty_policy_admitted\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
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
    assert_decision_compiler(&first, NASA9_COMPILER_ID);

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
fn g3_cli_compiles_first_order_kinetics_into_an_identical_verified_model_pack() {
    let (directory, manifest) = write_kinetics_fixture();
    let first_path = directory.join("kinetics-first.fsmodpk");
    let second_path = directory.join("kinetics-second.fsmodpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first kinetics compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second kinetics compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "kinetics decision stream moved"
    );
    assert_decision_compiler(&first, KINETICS_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first kinetics pack");
    let second_bytes = fs::read(second_path).expect("read second kinetics pack");
    assert_eq!(first_bytes, second_bytes, "kinetics pack bytes moved");
    let decoded =
        NormalizedModelPack::from_bytes(&first_bytes).expect("decode kinetics model pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedModelPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verified kinetics model pack");
    assert_eq!(decoded.pack_id(), "water-formation");
    assert_eq!(
        decoded.compiler(),
        "frankensim-matdb-kinetics-model-pack-compiler-v1"
    );
    assert_eq!(decoded.models().len(), 1);
    assert_eq!(decoded.normalizations().len(), 4);
    let card = &decoded.models()[0];
    assert_eq!(card.law.0, "arrhenius-first-order-rate");
    assert_eq!(card.law_version, 1);
    assert_eq!(
        card.parameters["activation_temperature"].dims,
        Dims([0, 0, 0, 1, 0, 0])
    );
    assert_eq!(
        card.parameters["pre_exponential"].dims,
        Dims([0, 0, -1, 0, 0, 0])
    );
    assert!(
        card.provenance
            .source
            .contains("[reaction:water-formation]")
    );
    assert!(card.provenance.source.contains("[rate-basis:first-order]"));

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"reaction_pack_id_bound\""));
    assert!(decisions.contains("\"reason_code\":\"kinetics_reaction_normalized\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_model_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
#[allow(clippy::too_many_lines)] // The end-to-end receipt audit is clearer as one ordered assertion path.
fn g3_cli_compiles_species_association_into_identical_verified_species_packs() {
    let (directory, manifest) = write_species_fixture(SPECIES_SOURCE);
    let first_path = directory.join("species-first.fsspcpk");
    let second_path = directory.join("species-second.fsspcpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first species compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second species compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "species decision stream moved");
    assert_decision_compiler(&first, SPECIES_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first species pack");
    let second_bytes = fs::read(second_path).expect("read second species pack");
    assert_eq!(first_bytes, second_bytes, "species pack bytes moved");
    let decoded = NormalizedSpeciesPack::from_bytes(&first_bytes).expect("decode species pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedSpeciesPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verified species pack");
    assert_eq!(decoded.pack_id(), "N2");
    assert_eq!(
        decoded.compiler(),
        "frankensim-matdb-species-pack-compiler-v1"
    );
    let association = decoded.association();
    assert_eq!(association.species().as_str(), "N2");
    assert_eq!(association.molar_mass().to_bits(), 0.028_013_4f64.to_bits());
    assert_eq!(association.standard_state_phase(), "gas");
    assert_eq!(association.reference_eos(), "ideal-gas");
    assert_eq!(
        association.reference_pressure().to_bits(),
        100_000.0f64.to_bits()
    );
    assert_eq!(association.elemental_reference(), "NASA-TP-2002-211556");
    assert_eq!(association.sources().len(), 2);
    assert_eq!(association.provenance().license.as_str(), "CC-BY-4.0");
    assert!(
        association
            .provenance()
            .artifact
            .is_some_and(|artifact| association.sources().contains(&artifact))
    );
    assert_eq!(decoded.normalizations().len(), 2);
    assert_eq!(
        decoded.normalizations()[0].target(),
        SpeciesNormalizationTarget::MolarMass
    );
    assert_eq!(decoded.normalizations()[0].dims(), SPECIES_MOLAR_MASS_DIMS);
    assert_eq!(
        decoded.normalizations()[0].scale().to_bits(),
        0.001f64.to_bits()
    );
    assert_eq!(
        decoded.normalizations()[0].offset().to_bits(),
        0.0f64.to_bits()
    );
    assert_eq!(decoded.normalizations()[0].source_basis(), "g/mol");
    assert_eq!(
        decoded.normalizations()[0].target_basis(),
        SPECIES_PACK_TARGET_BASIS
    );
    assert_eq!(
        decoded.normalizations()[1].target(),
        SpeciesNormalizationTarget::ReferencePressure
    );
    assert_eq!(
        decoded.normalizations()[1].dims(),
        SPECIES_REFERENCE_PRESSURE_DIMS
    );
    assert_eq!(
        decoded.normalizations()[1].scale().to_bits(),
        1_000.0f64.to_bits()
    );

    let decisions = String::from_utf8(first.stdout).expect("decision stream is UTF-8");
    assert!(decisions.contains("\"reason_code\":\"species_pack_id_bound\""));
    assert!(decisions.contains("\"reason_code\":\"species_association_normalized\""));
    assert!(decisions.contains("\"reason_code\":\"runtime_species_pack_self_verified\""));
    assert!(
        decisions
            .lines()
            .all(|row| row.contains(&format!("\"pack_hash\":\"{pack_hash}\"")))
    );
}

#[test]
fn g3_cli_compiles_committed_methane_seed_and_records_independent_agreement() {
    let manifest = workspace_path(METHANE_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed methane seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("methane-first.fsspcpk");
    let second_path = directory.join("methane-second.fsspcpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first methane seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second methane seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(first.stdout, second.stdout, "methane decision stream moved");
    assert_decision_compiler(&first, SPECIES_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first methane pack");
    let second_bytes = fs::read(second_path).expect("read second methane pack");
    assert_eq!(first_bytes, second_bytes, "methane pack bytes moved");
    let decoded = NormalizedSpeciesPack::from_bytes(&first_bytes).expect("decode methane pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedSpeciesPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify methane pack identity");

    assert_eq!(decoded.pack_id(), "CH4");
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    let association = decoded.association();
    assert_eq!(association.species().as_str(), "CH4");
    assert_eq!(
        association.molar_mass().to_bits(),
        (NASA_METHANE_MOLAR_MASS_G_PER_MOL * 0.001).to_bits()
    );
    assert_eq!(association.standard_state_phase(), "gas");
    assert_eq!(association.reference_eos(), "ideal-gas");
    assert_eq!(
        association.reference_pressure().to_bits(),
        100_000.0f64.to_bits()
    );
    assert_eq!(
        association.elemental_reference(),
        "NASA-TP-2002-211556-reference-elements-298.15K-1bar"
    );
    assert_eq!(association.provenance().license.as_str(), NASA_SEED_LICENSE);
    assert!(
        association
            .provenance()
            .source
            .contains("NASA/TP-2002-211556")
    );
    assert!(association.provenance().source.contains("[source:primary]"));

    let independent_difference =
        (association.molar_mass() - NIST_SRD69_METHANE_MOLAR_MASS_KG_PER_MOL).abs();
    assert!(
        independent_difference <= NIST_SRD69_DISPLAY_ROUNDING_HALF_WIDTH_KG_PER_MOL,
        "NASA seed and NIST SRD 69 display disagree beyond the recorded rounding band: {independent_difference:e} kg/mol"
    );
}

#[test]
fn g3_cli_compiles_committed_air_exhaust_constituents_without_inventing_a_mixture() {
    for seed in AIR_EXHAUST_SPECIES_SEEDS {
        let manifest = workspace_path(seed.manifest);
        assert!(
            manifest.is_file(),
            "committed {} seed manifest is missing",
            seed.species
        );
        let directory = fixture_dir();
        let first_path = directory.join(format!("{}-first.fsspcpk", seed.species));
        let second_path = directory.join(format!("{}-second.fsspcpk", seed.species));

        let first = run_compiler(&manifest, &first_path);
        let second = run_compiler(&manifest, &second_path);
        assert!(
            first.status.success(),
            "first {} seed compilation failed: {}",
            seed.species,
            String::from_utf8_lossy(&first.stderr)
        );
        assert!(
            second.status.success(),
            "second {} seed compilation failed: {}",
            seed.species,
            String::from_utf8_lossy(&second.stderr)
        );
        assert_eq!(
            first.stdout, second.stdout,
            "{} decision stream moved",
            seed.species
        );
        assert_decision_compiler(&first, SPECIES_COMPILER_ID);

        let first_bytes = fs::read(first_path).expect("read first constituent pack");
        let second_bytes = fs::read(second_path).expect("read second constituent pack");
        assert_eq!(
            first_bytes, second_bytes,
            "{} pack bytes moved",
            seed.species
        );
        let decoded =
            NormalizedSpeciesPack::from_bytes(&first_bytes).expect("decode constituent pack");
        let pack_hash = decoded.content_hash();
        let decoded = NormalizedSpeciesPack::from_bytes_verified(pack_hash, &first_bytes)
            .expect("verify constituent pack identity");

        assert_eq!(decoded.pack_id(), seed.species);
        assert!(
            decoded
                .redistribution_terms()
                .contains("public use permitted")
        );
        let association = decoded.association();
        assert_eq!(association.species().as_str(), seed.species);
        assert_eq!(
            association.molar_mass().to_bits(),
            (seed.nasa_molar_mass_g_per_mol * 0.001).to_bits()
        );
        assert_eq!(association.standard_state_phase(), "gas");
        assert_eq!(association.reference_eos(), "ideal-gas");
        assert_eq!(
            association.reference_pressure().to_bits(),
            100_000.0f64.to_bits()
        );
        assert_eq!(
            association.elemental_reference(),
            "NASA-TP-2002-211556-reference-elements-298.15K-1bar"
        );
        assert_eq!(association.provenance().license.as_str(), NASA_SEED_LICENSE);
        assert!(
            association
                .provenance()
                .source
                .contains("NASA/TP-2002-211556")
        );

        let independent_difference_g_per_mol =
            (association.molar_mass() * 1_000.0 - seed.nist_molar_mass_g_per_mol).abs();
        assert!(
            independent_difference_g_per_mol <= seed.nist_display_rounding_half_width_g_per_mol,
            "NASA {} seed and NIST SRD 69 display disagree beyond the recorded rounding band: {independent_difference_g_per_mol:e} g/mol",
            seed.species
        );
    }
}

#[test]
fn g3_cli_refuses_malformed_species_without_publishing() {
    let malformed = SPECIES_SOURCE.replacen("28.0134\tg/mol", "28.0134\tkg", 1);
    assert_ne!(malformed, SPECIES_SOURCE);
    let (directory, manifest) = write_species_fixture(&malformed);
    let output = directory.join("refused-species.fsspcpk");

    let refused = run_compiler(&manifest, &output);
    assert!(
        !refused.status.success(),
        "invalid species unexpectedly compiled"
    );
    assert!(!output.exists(), "species refusal published an output");
    assert_decision_compiler(&refused, SPECIES_COMPILER_ID);
    let decisions = String::from_utf8(refused.stdout).expect("decision stream is UTF-8");
    assert_eq!(decisions.matches("\"verdict\":\"refuse\"").count(), 1);
    assert!(decisions.contains("\"reason_code\":\"species_molar_mass_dims_mismatch\""));
    assert!(decisions.contains("\"subject\":\"species:N2\""));
    assert!(
        String::from_utf8_lossy(&refused.stderr)
            .contains("error: matdb pack refused [species_molar_mass_dims_mismatch]")
    );
}

#[test]
fn g3_cli_uses_generic_driver_identity_for_an_unknown_profile() {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    let unsupported = MANIFEST.replace("material-tsv-v1", "future-profile-v1");
    fs::write(&manifest, unsupported).expect("write unsupported-profile manifest");
    fs::write(directory.join("source.tsv"), SOURCE).expect("write source fixture");
    let output = directory.join("unsupported.fsmatpk");

    let refused = run_compiler(&manifest, &output);
    assert!(!refused.status.success());
    assert!(!output.exists(), "unsupported profile published an output");
    assert_decision_compiler(&refused, MATERIAL_COMPILER_ID);
    assert!(
        String::from_utf8_lossy(&refused.stdout)
            .contains("\"reason_code\":\"unsupported_source_profile\"")
    );
}

#[test]
fn g3_cli_uses_generic_driver_identity_before_profile_selection() {
    let directory = fixture_dir();
    let manifest = directory.join("manifest.tsv");
    let incomplete = MANIFEST.replace("license\tCC-BY-4.0\n", "");
    fs::write(&manifest, incomplete).expect("write incomplete manifest");
    fs::write(directory.join("source.tsv"), SOURCE).expect("write source fixture");
    let output = directory.join("incomplete.fsmatpk");

    let refused = run_compiler(&manifest, &output);
    assert!(!refused.status.success());
    assert!(!output.exists(), "incomplete manifest published an output");
    assert_decision_compiler(&refused, MATERIAL_COMPILER_ID);
    assert!(
        String::from_utf8_lossy(&refused.stdout).contains("\"reason_code\":\"missing_license\"")
    );
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
