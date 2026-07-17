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
const NASA_SEED_LICENSE: &str = "Work-of-the-US-Government-Public-Use-Permitted";
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
