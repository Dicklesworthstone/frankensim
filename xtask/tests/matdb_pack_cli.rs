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
const PTFE_TEFLON_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/ptfe-teflon-nist-cryogenic/manifest.tsv";
const PEEK_THERMIC_SEED_MANIFEST: &str = "data/matdb/seed-v1/peek-nasa-thermic-plate/manifest.tsv";
const NASA_CR_115153_WATER_ETHYLENE_GLYCOL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-cr-115153-water-ethylene-glycol/manifest.tsv";
const N0602_001_NITRILE_JP8_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/n0602-001-nitrile-jp8-compatibility/manifest.tsv";
const NASA_TN_D_8184_M19_MATERIAL_DECK_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-tn-d-8184-m19-material-deck/manifest.tsv";
const NGYC_N42_SINTERED_NICKEL_COATED_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/ngyc-n42-sintered-nickel-coated/manifest.tsv";
const JINSHAN_N42_PRISTINE_TEMPERATURE_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/jinshan-n42-pristine-temperature/manifest.tsv";
const NACA_TN_2680_ISOOCTANE_FLAME_SPEED_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/naca-tn-2680-isooctane-flame-speed/manifest.tsv";
const FACE_G_CDTRF_G_2023_V1_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/face-g-cdtrf-g-2023-v1/manifest.tsv";
const WO2018_125520_FORMULATION_8_5W30_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/wo2018-125520-formulation-8-5w30/manifest.tsv";
const NASA_UAM_MW16C_POLYIMIDE_WIRE_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-uam-mw16c-polyimide-magnet-wire/manifest.tsv";
const NASA_UAM_NOMEX_410_SLOT_LINER_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-uam-nomex-410-slot-liner/manifest.tsv";
const NASA_UAM_COOLTHERM_EP2000_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-uam-cooltherm-ep2000-180c-cure/manifest.tsv";
const AISI_4140_RC33_SEED_MANIFEST: &str = "data/matdb/seed-v1/aisi-4140-rc33/manifest.tsv";
const AISI_1045_COLD_DRAWN_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-1045-cold-drawn/manifest.tsv";
const AISI_52100_CVM_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-52100-cvm-hot-hardness/manifest.tsv";
const AISI_9310_CVM_CARBURIZED_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/aisi-9310-cvm-carburized/manifest.tsv";
const NAPC_PE_5_L_1274_GEAR_OIL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/napc-pe-5-l-1274-gear-oil/manifest.tsv";
const NAPC_PE_5_L_1307_1553_GEAR_OIL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/napc-pe-5-l-1307-1553-gear-oil/manifest.tsv";
const RHEOLUBE_2000_PENNZANE_GREASE_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/rheolube-2000-pennzane-grease/manifest.tsv";
const PENNZANE_SHF_X_2000_BEARING_OIL_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/pennzane-shf-x-2000-bearing-oil/manifest.tsv";
const GRAY_CAST_IRON_S2_S_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/gray-cast-iron-s2-s/manifest.tsv";
const NASA_CR_195445_OMC_PS200_ROTARY_COATING_SEED_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-cr-195445-omc-ps200-rotary-coating/manifest.tsv";
const NASA_SEED_LICENSE: &str = "Work-of-the-US-Government-Public-Use-Permitted";
const PUBLIC_USE_PERMITTED_LICENSE: &str = "Public-Use-Permitted";
const PUBLIC_USE_AND_PATENT_PUBLICATION_LICENSE: &str =
    "Public-Use-Permitted-and-US-Patent-Publication";
const CC_BY_4_0_LICENSE: &str = "CC-BY-4.0";
const NIST_PUBLIC_INFORMATION_LICENSE: &str = "NIST-Public-Information-Attribution-Requested";
const USPTO_PATENT_TEXT_LICENSE: &str = "USPTO-Patent-Text-Typically-No-Copyright-Restrictions";
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
fn g3_cli_compiles_committed_ptfe_teflon_cryogenic_seed() {
    let manifest = workspace_path(PTFE_TEFLON_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed PTFE/Teflon seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("ptfe-teflon-first.fsmatpk");
    let second_path = directory.join("ptfe-teflon-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first PTFE/Teflon seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second PTFE/Teflon seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "PTFE/Teflon decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first PTFE/Teflon pack");
    let second_bytes = fs::read(second_path).expect("read second PTFE/Teflon pack");
    assert_eq!(first_bytes, second_bytes, "PTFE/Teflon pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode PTFE/Teflon pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify PTFE/Teflon pack identity");

    assert_eq!(decoded.pack_id(), "ptfe-teflon-nist-cryogenic");
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
            0.232_391_801_023_681,
            Dims([1, 1, -3, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=2.7380,-30.677",
        ),
        (
            "thermal_conductivity",
            293.0,
            0.272_587_209_362_470,
            Dims([1, 1, -3, -1, 0, 0]),
            "5 percent curve-fit error",
            "a..i=2.7380,-30.677",
        ),
        (
            "specific_heat_capacity",
            77.0,
            301.115_701_345_352,
            Dims([2, 0, -2, -1, 0, 0]),
            "1.5 percent curve-fit error",
            "a..i=31.88256,-166.51949",
        ),
        (
            "specific_heat_capacity",
            293.0,
            1_015.656_817_896_37,
            Dims([2, 0, -2, -1, 0, 0]),
            "1.5 percent curve-fit error",
            "a..i=31.88256,-166.51949",
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
            .unwrap_or_else(|| panic!("missing PTFE/Teflon {property} claim at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PTFE/Teflon {property} at {temperature} K was not an exact-point scalar");
        };
        assert_eq!(
            *dims, expected_dims,
            "PTFE/Teflon {property} dimensions moved"
        );
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "PTFE/Teflon {property} at {temperature} K moved by {relative_error:e} relative"
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
                .contains("Material Properties: Teflon")
        );
        assert!(claim.provenance.source.contains("[source:nist-ptfe-fit]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PTFE/Teflon claim observation remains linked");
        assert_eq!(
            observation.specimen,
            "ptfe-teflon-source-grade-and-process-unspecified"
        );
        assert!(observation.method.contains("NIST Teflon"));
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
                .contains("data and equation range 4-300 K")
        );
        assert!(
            observation
                .caveats
                .contains("without a confidence level or degrees of freedom")
        );
        assert!(
            observation
                .caveats
                .contains("does not identify resin grade")
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
fn g3_cli_compiles_committed_peek_thermic_plate_seed() {
    let manifest = workspace_path(PEEK_THERMIC_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed PEEK THERMIC seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("peek-thermic-first.fsmatpk");
    let second_path = directory.join("peek-thermic-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first PEEK THERMIC seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second PEEK THERMIC seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "PEEK THERMIC decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first PEEK THERMIC pack");
    let second_bytes = fs::read(second_path).expect("read second PEEK THERMIC pack");
    assert_eq!(first_bytes, second_bytes, "PEEK THERMIC pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode PEEK THERMIC pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify PEEK THERMIC pack identity");

    assert_eq!(decoded.pack_id(), "peek-nasa-thermic-plate-2021");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 9);
    assert!(decoded.joint_statistics().is_empty());

    let conductivity = [
        (300.0, 0.224_458_9),
        (400.0, 0.243_077_8),
        (500.0, 0.265_855_5),
        (525.0, 0.274_943_823_437_5),
    ];
    for (temperature, expected_value) in conductivity {
        let (_, claim) = decoded
            .claims()
            .claims_for("thermal_conductivity")
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing PEEK conductivity at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PEEK conductivity at {temperature} K was not scalar");
        };
        assert_eq!(*dims, Dims([1, 1, -3, -1, 0, 0]));
        assert_eq!(*value, expected_value);
        assert_eq!(
            claim.validity.bound("source_pressure_atmospheric"),
            Some((1.0, 1.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA/TM-20210014330"));
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nasa-thermic-peek]")
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PEEK conductivity observation remains linked");
        assert_eq!(
            observation.specimen,
            "nasa-larc-thermic-peek-plate-grade-and-process-unspecified"
        );
        assert!(observation.method.contains("Continuous Genetic Algorithm"));
        assert!(observation.caveats.contains("c0..c3=-4.0607e-2"));
        assert!(
            observation
                .caveats
                .contains("narrower repeated range governs")
        );
        assert!(observation.caveats.contains("differed by about 3 percent"));
        assert!(observation.caveats.contains("does not identify PEEK grade"));
    }

    let specific_heat = [
        (300.0, 1_058.931),
        (400.0, 1_347.916),
        (500.0, 1_765.685),
        (525.0, 1_897.616_390_625),
    ];
    for (temperature, expected_value) in specific_heat {
        let (_, claim) = decoded
            .claims()
            .claims_for("specific_heat_capacity")
            .into_iter()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing PEEK specific heat at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PEEK specific heat at {temperature} K was not scalar");
        };
        assert_eq!(*dims, Dims([2, 0, -2, -1, 0, 0]));
        assert_eq!(*value, expected_value);
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PEEK specific-heat observation remains linked");
        assert!(
            observation
                .method
                .contains("differential-scanning-calorimeter")
        );
        assert!(observation.caveats.contains("1.0477e-5*T^3"));
        assert!(
            observation
                .caveats
                .contains("Equation 1's dimensional balance")
        );
        assert!(observation.caveats.contains("no residual, dispersion"));
    }

    let density_claims = decoded.claims().claims_for("density");
    assert_eq!(density_claims.len(), 1);
    let (_, density) = density_claims[0];
    let PropertyValue::Scalar { value, dims } = &density.value else {
        panic!("PEEK density was not scalar");
    };
    assert_eq!(*value, 1_264.0);
    assert_eq!(*dims, Dims([-3, 1, 0, 0, 0, 0]));
    assert_eq!(
        density.validity.bound("source_test_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(density.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(density.provenance.license, NASA_SEED_LICENSE);
    let density_observation = decoded
        .claims()
        .observation(density.observations[0])
        .expect("PEEK density observation remains linked");
    assert!(density_observation.method.contains("Commercial-laboratory"));
    assert!(
        density_observation
            .caveats
            .contains("Netzsch report 621004797")
    );
    assert!(density_observation.caveats.contains("test temperature"));

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
fn g3_cli_compiles_committed_nasa_cr_115153_water_ethylene_glycol_seed() {
    let manifest = workspace_path(NASA_CR_115153_WATER_ETHYLENE_GLYCOL_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NASA-CR-115153 water/ethylene-glycol seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nasa-cr-115153-water-glycol-first.fsmatpk");
    let second_path = directory.join("nasa-cr-115153-water-glycol-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-CR-115153 water/glycol seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-CR-115153 water/glycol seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NASA-CR-115153 water/glycol decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NASA water/glycol pack");
    let second_bytes = fs::read(second_path).expect("read second NASA water/glycol pack");
    assert_eq!(
        first_bytes, second_bytes,
        "NASA-CR-115153 water/glycol pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode NASA water/glycol pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NASA water/glycol pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nasa-cr-115153-inhibited-water-ethylene-glycol-coolant"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 11);
    assert!(decoded.joint_statistics().is_empty());

    let expected_bounds = [
        ("sodium_nitrite_mass_fraction_lower_bound", 0.10),
        ("sodium_nitrite_mass_fraction_upper_bound", 0.25),
        ("sodium_benzoate_mass_fraction_lower_bound", 1.33),
        ("sodium_benzoate_mass_fraction_upper_bound", 1.57),
        ("water_mass_fraction_lower_bound", 36.0),
        ("water_mass_fraction_upper_bound", 38.5),
    ];
    for (property, source_percent) in expected_bounds {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique {property} claim");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("{property} was not scalar");
        };
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "{property} moved by {relative_error:e} relative"
        );
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-CR-115153"));
        assert!(
            claim
                .provenance
                .source
                .contains("[source:nasa-cr-115153-table-5]")
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("composition-bound observation remains linked");
        assert_eq!(
            observation.specimen,
            "nasa-cr-115153-water-ethylene-glycol-sodium-nitrite-sodium-benzoate-solution"
        );
        assert!(observation.method.contains("formulation specification"));
        assert!(observation.caveats.contains("formulation bounds"));
        assert!(observation.caveats.contains("without inventing a midpoint"));
    }
    assert!(
        decoded
            .claims()
            .claims_for("ethylene_glycol_mass_fraction")
            .is_empty(),
        "an unreported ethylene-glycol balance must not be inferred"
    );

    let bulk_properties = [
        (
            "density",
            1_081.246_277_742_31,
            Dims([-3, 1, 0, 0, 0, 0]),
            false,
        ),
        (
            "thermal_conductivity",
            0.380_761_626_601_706,
            Dims([1, 1, -3, -1, 0, 0]),
            true,
        ),
    ];
    for (property, expected_value, expected_dims, btu_conversion) in bulk_properties {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique water/glycol {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("water/glycol {property} was not scalar");
        };
        assert_eq!(*value, expected_value);
        assert_eq!(*dims, expected_dims);
        assert_eq!(
            claim.validity.bound("source_test_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_test_pressure_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("bulk-property observation remains linked");
        assert!(observation.caveats.contains("extra SI digits"));
        assert!(observation.caveats.contains("not source precision"));
        if btu_conversion {
            assert_eq!(
                claim.validity.bound("source_btu_convention_known"),
                Some((0.0, 0.0))
            );
            assert!(observation.caveats.contains("Btu_IT=1055.05585262 J"));
        }
    }

    let specific_heat_points = [
        (255.372_222_222_222, 2_805.156),
        (283.15, 4_479.876),
        (310.927_777_777_778, 6_154.596),
    ];
    let specific_heat_claims = decoded.claims().claims_for("specific_heat_capacity");
    assert_eq!(specific_heat_claims.len(), specific_heat_points.len());
    for (temperature, expected_value) in specific_heat_points {
        let (_, claim) = specific_heat_claims
            .iter()
            .copied()
            .find(|(_, claim)| {
                claim.validity.bound("temperature") == Some((temperature, temperature))
            })
            .unwrap_or_else(|| panic!("missing NASA water/glycol cp at {temperature} K"));
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("NASA water/glycol cp at {temperature} K was not scalar");
        };
        assert_eq!(*value, expected_value);
        assert_eq!(*dims, Dims([2, 0, -2, -1, 0, 0]));
        assert_eq!(
            claim.validity.bound("source_test_pressure_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_btu_convention_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("specific-heat observation remains linked");
        assert!(observation.method.contains("three exact temperatures"));
        assert!(
            observation
                .caveats
                .contains("approximately cp=(0.67 + 0.008*T_degF)")
        );
        assert!(
            observation
                .caveats
                .contains("do not expose a continuous law")
        );
    }

    // G3 comparison evidence only: NASA/TM-2019-220019 Table VIII lists a
    // separately sourced, composition-basis-unspecified 50-50 water/ethylene-
    // glycol fluid at 1082 kg/m3 and 0.402 W/(m K). Those condition-mismatched
    // values do not overwrite this NASA-CR-115153 formulation; they only bound
    // a coarse transcription plausibility check.
    let (_, density) = decoded.claims().claims_for("density")[0];
    let PropertyValue::Scalar {
        value: density_value,
        ..
    } = &density.value
    else {
        panic!("NASA water/glycol density was not scalar");
    };
    let density_relative_difference = (*density_value - 1_082.0_f64).abs() / 1_082.0;
    assert!(density_relative_difference <= 0.001);

    let (_, conductivity) = decoded.claims().claims_for("thermal_conductivity")[0];
    let PropertyValue::Scalar {
        value: conductivity_value,
        ..
    } = &conductivity.value
    else {
        panic!("NASA water/glycol conductivity was not scalar");
    };
    let conductivity_relative_difference = (*conductivity_value - 0.402_f64).abs() / 0.402;
    assert!(conductivity_relative_difference <= 0.06);

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
fn g3_cli_compiles_committed_n0602_001_nitrile_jp8_compatibility_seed() {
    let manifest = workspace_path(N0602_001_NITRILE_JP8_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed N0602-001 nitrile seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("n0602-001-nitrile-first.fsmatpk");
    let second_path = directory.join("n0602-001-nitrile-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first N0602-001 nitrile seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second N0602-001 nitrile seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "N0602-001 nitrile decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first N0602-001 nitrile pack");
    let second_bytes = fs::read(second_path).expect("read second N0602-001 nitrile pack");
    assert_eq!(
        first_bytes, second_bytes,
        "N0602-001 nitrile pack bytes moved"
    );
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode N0602-001 nitrile pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify N0602-001 nitrile pack identity");

    assert_eq!(
        decoded.pack_id(),
        "n0602-001-nitrile-o-ring-jp8-compatibility"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 10);
    assert!(decoded.joint_statistics().is_empty());

    let tga_claims = decoded
        .claims()
        .claims_for("tga_semivolatile_mass_fraction");
    assert_eq!(tga_claims.len(), 1);
    let (_, tga) = tga_claims[0];
    let PropertyValue::Scalar { value, dims } = &tga.value else {
        panic!("N0602-001 TGA semi-volatiles were not scalar");
    };
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    let tga_expected = 10.1_f64 * 0.01;
    assert!((*value - tga_expected).abs() / tga_expected <= 2.0e-15);
    assert_eq!(
        tga.validity.bound("source_tga_temperature_program_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(tga.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(tga.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
    assert!(tga.provenance.source.contains("NTRS 20080003822"));
    assert!(tga.provenance.source.contains("[source:primary]"));
    let tga_observation = decoded
        .claims()
        .observation(tga.observations[0])
        .expect("N0602-001 TGA observation remains linked");
    assert_eq!(
        tga_observation.specimen,
        "n0602-001-nitrile-rubber-o-ring-source-formulation-and-lot-unspecified"
    );
    assert!(tga_observation.method.contains("Thermogravimetric"));
    assert!(tga_observation.caveats.contains("propensity to shrink"));
    assert!(tga_observation.caveats.contains("compound formulation"));

    let absorbed_claims = decoded.claims().claims_for("absorbed_fuel_volume_fraction");
    assert_eq!(absorbed_claims.len(), 2);
    for (source_aromatic_percent, source_absorbed_percent) in [(0.0, 8.7), (25.0, 27.9)] {
        let aromatic_fraction = source_aromatic_percent * 0.01;
        let (_, claim) = absorbed_claims
            .iter()
            .copied()
            .find(|(_, claim)| {
                claim.validity.bound("fuel_aromatic_volume_fraction")
                    == Some((aromatic_fraction, aromatic_fraction))
            })
            .unwrap_or_else(|| {
                panic!("missing N0602-001 fuel absorption at {source_aromatic_percent}% aromatic")
            });
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("N0602-001 fuel absorption was not scalar");
        };
        let expected_value = source_absorbed_percent * 0.01;
        assert!((*value - expected_value).abs() / expected_value <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(
            claim.validity.bound("source_exposure_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_duration_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("N0602-001 fuel-absorption observation remains linked");
        assert!(observation.method.contains("thermal-desorption GC-MS"));
        assert!(
            observation
                .caveats
                .contains("does not define service compatibility")
        );
    }

    for (property, expected_value) in [
        ("jp8_alkane_fuel_polymer_partition_coefficient", 0.120),
        ("jp8_aromatic_fuel_polymer_partition_coefficient", 0.412),
        ("jp8_aromatic_to_alkane_partition_ratio", 3.4),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique N0602-001 {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("N0602-001 {property} was not scalar");
        };
        assert_eq!(*value, expected_value);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        assert_eq!(
            claim.validity.bound("fuel_aromatic_volume_fraction"),
            Some((0.0, 25.0_f64 * 0.01))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_duration_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
    }

    let slope_claims = decoded
        .claims()
        .claims_for("jp8_volume_swell_per_aromatic_volume_fraction");
    assert_eq!(slope_claims.len(), 2);
    let mut slopes = slope_claims
        .iter()
        .map(|(_, claim)| match &claim.value {
            PropertyValue::Scalar { value, dims } => {
                assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
                *value
            }
            PropertyValue::Curve { .. } => panic!("N0602-001 slope was not scalar"),
        })
        .collect::<Vec<_>>();
    slopes.sort_by(f64::total_cmp);
    assert_eq!(slopes, vec![0.451, 0.463]);
    assert_ne!(
        slope_claims[0].1.observations[0], slope_claims[1].1.observations[0],
        "conflicting printed slopes must retain distinct observations"
    );
    for (_, claim) in slope_claims {
        assert_eq!(
            claim.validity.bound("fuel_aromatic_volume_fraction"),
            Some((0.0, 25.0_f64 * 0.01))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_temperature_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_exposure_duration_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("N0602-001 slope observation remains linked");
        assert!(observation.caveats.contains("conflict"));
    }

    let r_squared = decoded
        .claims()
        .claims_for("jp8_volume_swell_aromatic_fraction_r_squared");
    assert_eq!(r_squared.len(), 1);
    let PropertyValue::Scalar { value, dims } = &r_squared[0].1.value else {
        panic!("N0602-001 R-squared was not scalar");
    };
    assert_eq!(*value, 0.948);
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    assert_eq!(
        r_squared[0]
            .1
            .validity
            .bound("source_exposure_temperature_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        r_squared[0].1.provenance.license,
        PUBLIC_USE_PERMITTED_LICENSE
    );

    let intercepts = decoded
        .claims()
        .claims_for("jp8_volume_swell_zero_aromatic_intercept");
    assert_eq!(intercepts.len(), 1);
    let PropertyValue::Scalar { value, dims } = &intercepts[0].1.value else {
        panic!("N0602-001 regression intercept was not scalar");
    };
    let expected_intercept = -1.167_f64 * 0.01;
    assert!((*value - expected_intercept).abs() / expected_intercept.abs() <= 2.0e-15);
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    assert_eq!(
        intercepts[0]
            .1
            .validity
            .bound("source_exposure_duration_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        intercepts[0].1.provenance.license,
        PUBLIC_USE_PERMITTED_LICENSE
    );
    let intercept_observation = decoded
        .claims()
        .observation(intercepts[0].1.observations[0])
        .expect("N0602-001 intercept observation remains linked");
    assert!(
        intercept_observation
            .caveats
            .contains("not a certified shrinkage value")
    );
    assert!(
        decoded
            .claims()
            .claims_for("jp8_prediction_interval_overlap")
            .is_empty(),
        "the source's approximate 57% overlap must remain observation-only"
    );

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
fn g3_cli_compiles_committed_nasa_tn_d_8184_m19_material_deck_without_inventing_process_state() {
    let manifest = workspace_path(NASA_TN_D_8184_M19_MATERIAL_DECK_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NASA-TN-D-8184 M-19 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nasa-tn-d-8184-m19-first.fsmatpk");
    let second_path = directory.join("nasa-tn-d-8184-m19-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-TN-D-8184 M-19 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-TN-D-8184 M-19 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NASA-TN-D-8184 M-19 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NASA-TN-D-8184 M-19 pack");
    let second_bytes = fs::read(second_path).expect("read second NASA-TN-D-8184 M-19 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "NASA-TN-D-8184 M-19 pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode NASA-TN-D-8184 M-19 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NASA-TN-D-8184 M-19 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nasa-tn-d-8184-m19-silicon-steel-material-deck"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("public use"));
    assert_eq!(decoded.claims().claim_count(), 6);
    assert!(decoded.joint_statistics().is_empty());

    let magnetization = decoded.claims().claims_for("magnetic_flux_density");
    assert_eq!(magnetization.len(), 1);
    let magnetization = magnetization[0].1;
    let PropertyValue::Curve {
        abscissa,
        abscissa_dims,
        knots,
        dims,
    } = &magnetization.value
    else {
        panic!("NASA-TN-D-8184 M-19 magnetization data was not a curve");
    };
    assert_eq!(abscissa, "magnetic_field_strength");
    assert_eq!(*abscissa_dims, Dims([-1, 0, 0, 0, 1, 0]));
    assert_eq!(*dims, Dims([0, 1, -2, 0, -1, 0]));
    assert_eq!(knots.len(), 14);
    assert_eq!(
        magnetization.interpolation,
        InterpolationPolicy::TabulatedOnly
    );
    assert_eq!(magnetization.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(magnetization.provenance.license, NASA_SEED_LICENSE);
    assert!(magnetization.provenance.source.contains("NASA-TN-D-8184"));
    assert!(magnetization.provenance.source.contains("[source:primary]"));
    for missing_axis in [
        "source_manufacturer_known",
        "source_processing_anneal_state_known",
        "source_lamination_thickness_known",
        "source_magnetic_test_method_known",
        "source_test_frequency_known",
        "source_test_temperature_known",
        "source_chemistry_known",
        "source_waveform_known",
        "source_direction_known",
    ] {
        assert_eq!(
            magnetization.validity.bound(missing_axis),
            Some((0.0, 0.0)),
            "M-19 B-H curve must retain missing identity axis {missing_axis}"
        );
    }
    assert_eq!(
        magnetization
            .validity
            .bound("source_curve_points_printed_not_digitized"),
        Some((1.0, 1.0))
    );

    let source_b_kilolines_per_square_inch = [
        26.0_f64, 30.0, 40.0, 50.0, 60.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 110.0, 116.0,
    ];
    let source_h_ampere_turns_per_inch = [
        1.30_f64, 1.45, 1.95, 2.55, 3.50, 5.1, 6.5, 8.8, 13.0, 21.0, 37.0, 60.0, 130.0, 185.0,
    ];
    for ((actual_h, actual_b), (source_h, source_b)) in knots.iter().zip(
        source_h_ampere_turns_per_inch
            .iter()
            .zip(source_b_kilolines_per_square_inch.iter()),
    ) {
        let expected_h = source_h / 0.0254_f64;
        let expected_b = source_b * 1.0e-5_f64 / 0.0254_f64.powi(2);
        assert!((actual_h - expected_h).abs() / expected_h <= 2.0e-14);
        assert!((actual_b - expected_b).abs() / expected_b <= 2.0e-14);
    }

    let scalar_expectations = [
        (
            "specific_core_loss",
            9.4_f64 / 0.453_592_37_f64,
            Dims([2, 0, -3, 0, 0, 0]),
        ),
        (
            "core_loss_frequency_power_law_exponent",
            1.47_f64,
            Dims([0, 0, 0, 0, 0, 0]),
        ),
        (
            "lamination_thickness",
            0.014_f64 * 0.0254_f64,
            Dims([1, 0, 0, 0, 0, 0]),
        ),
        (
            "core_loss_reference_frequency",
            400.0_f64,
            Dims([0, 0, -1, 0, 0, 0]),
        ),
        (
            "core_loss_reference_flux_density",
            64.5_f64 * 1.0e-5_f64 / 0.0254_f64.powi(2),
            Dims([0, 1, -2, 0, -1, 0]),
        ),
    ];
    let mut scalar_observations = Vec::new();
    for (property, expected_value, expected_dims) in scalar_expectations {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing NASA M-19 {property}");
        let claim = claims[0].1;
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("NASA M-19 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let scale = expected_value.abs().max(1.0e-12_f64);
        assert!((*value - expected_value).abs() / scale <= 2.0e-15);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TN-D-8184"));
        assert_eq!(
            claim.validity.bound("source_material_grade_m19"),
            Some((1.0, 1.0))
        );
        for missing_axis in [
            "source_manufacturer_known",
            "source_chemistry_known",
            "source_processing_anneal_state_known",
            "source_magnetic_test_method_known",
            "source_test_temperature_known",
            "source_waveform_known",
            "source_direction_known",
        ] {
            assert_eq!(
                claim.validity.bound(missing_axis),
                Some((0.0, 0.0)),
                "NASA M-19 {property} must retain missing identity axis {missing_axis}"
            );
        }
        scalar_observations.push(claim.observations[0]);
    }
    assert!(
        scalar_observations
            .windows(2)
            .all(|pair| pair[0] == pair[1]),
        "NASA M-19 frequency-loss parameters must share one source observation"
    );

    let curve_observation = decoded
        .claims()
        .observation(magnetization.observations[0])
        .expect("NASA M-19 curve observation remains linked");
    assert!(curve_observation.method.contains("Figure 10"));
    assert!(curve_observation.caveats.contains("tabulated-only"));
    assert!(curve_observation.caveats.contains("anneal"));
    let loss_observation = decoded
        .claims()
        .observation(scalar_observations[0])
        .expect("NASA M-19 frequency-loss observation remains linked");
    assert!(loss_observation.method.contains("WCORE=9.4 W/lb"));
    assert!(
        loss_observation
            .caveats
            .contains("not a complete Steinmetz law")
    );
    assert!(loss_observation.caveats.contains("test method"));

    for refused_property in [
        "steinmetz_coefficient",
        "core_loss_flux_density_power_law_exponent",
        "recoil_relative_permeability",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent NASA M-19 property must remain refused: {refused_property}"
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
fn g3_cli_compiles_committed_ngyc_n42_sintered_magnet_seed() {
    let manifest = workspace_path(NGYC_N42_SINTERED_NICKEL_COATED_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NGYC N42 magnet seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("ngyc-n42-first.fsmatpk");
    let second_path = directory.join("ngyc-n42-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NGYC N42 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NGYC N42 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NGYC N42 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NGYC N42 pack");
    let second_bytes = fs::read(second_path).expect("read second NGYC N42 pack");
    assert_eq!(first_bytes, second_bytes, "NGYC N42 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode NGYC N42 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NGYC N42 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "ngyc-n42-sintered-ndfeb-nickel-coated-cubes"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 4);
    assert!(decoded.joint_statistics().is_empty());

    let remanence = decoded.claims().claims_for("remanent_flux_density");
    assert_eq!(remanence.len(), 1);
    let PropertyValue::Scalar { value, dims } = &remanence[0].1.value else {
        panic!("NGYC N42 remanence was not scalar");
    };
    let expected_remanence = 1350.0_f64 * 1.0e-3;
    assert!((*value - expected_remanence).abs() / expected_remanence <= 2.0e-15);
    assert_eq!(*dims, Dims([0, 1, -2, 0, -1, 0]));

    let coercivity = decoded.claims().claims_for("coercive_field_strength");
    assert_eq!(coercivity.len(), 1);
    let PropertyValue::Scalar { value, dims } = &coercivity[0].1.value else {
        panic!("NGYC N42 coercivity was not scalar");
    };
    assert_eq!(*value, 923.0_f64 * 1.0e3);
    assert_eq!(*dims, Dims([-1, 0, 0, 0, 1, 0]));

    let energy_products = decoded
        .claims()
        .claims_for("maximum_magnetic_energy_product");
    assert_eq!(energy_products.len(), 2);
    assert_ne!(
        energy_products[0].1.observations[0], energy_products[1].1.observations[0],
        "conflicting printed energy products must retain distinct observations"
    );
    let mut energy_values = energy_products
        .iter()
        .map(|(_, claim)| match &claim.value {
            PropertyValue::Scalar { value, dims } => {
                assert_eq!(*dims, Dims([-1, 1, -2, 0, 0, 0]));
                *value
            }
            PropertyValue::Curve { .. } => panic!("NGYC N42 energy product was not scalar"),
        })
        .collect::<Vec<_>>();
    energy_values.sort_by(f64::total_cmp);
    let printed_si = 318.3_f64 * 1.0e3;
    let normalized_42_mgoe = 42.0_f64 * (100_000.0 / (4.0 * std::f64::consts::PI));
    assert!((energy_values[0] - printed_si).abs() / printed_si <= 2.0e-15);
    assert!((energy_values[1] - normalized_42_mgoe).abs() / normalized_42_mgoe <= 2.0e-15);
    assert!(energy_values[1] > energy_values[0]);

    for property in [
        "remanent_flux_density",
        "coercive_field_strength",
        "maximum_magnetic_energy_product",
    ] {
        for (_, claim) in decoded.claims().claims_for(property) {
            assert_eq!(
                claim
                    .validity
                    .bound("source_magnetic_test_temperature_known"),
                Some((0.0, 0.0))
            );
            assert_eq!(
                claim.validity.bound("source_magnetic_test_method_known"),
                Some((0.0, 0.0))
            );
            assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
            assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
            assert!(
                claim
                    .provenance
                    .source
                    .contains("10.1038/s41598-023-47689-2")
            );
            assert!(claim.provenance.source.contains("[source:primary]"));
            let observation = decoded
                .claims()
                .observation(claim.observations[0])
                .expect("NGYC N42 observation remains linked");
            assert_eq!(
                observation.specimen,
                "ngyc-yinxian-ningbo-n42-sintered-ndfeb-nickel-coated-cubes-paper-lot-unspecified"
            );
        }
    }

    let si_observation = decoded
        .claims()
        .observation(remanence[0].1.observations[0])
        .expect("NGYC N42 SI observation remains linked");
    assert!(si_observation.method.contains("Telfah et al. 2023"));
    assert!(si_observation.caveats.contains("supplier nominal values"));
    assert!(si_observation.caveats.contains("not SI-equivalent"));
    let cgs_observation = energy_products
        .iter()
        .find_map(|(_, claim)| {
            let observation = decoded.claims().observation(claim.observations[0])?;
            observation
                .method
                .contains("Exact unit normalization")
                .then_some(observation)
        })
        .expect("NGYC N42 CGS observation remains linked");
    assert!(cgs_observation.method.contains("1 Oe=1000/(4*pi) A/m"));
    assert!(cgs_observation.caveats.contains("source conflict"));

    for refused_property in [
        "intrinsic_coercive_field_strength",
        "recoil_relative_permeability",
        "remanence_temperature_coefficient",
        "coercivity_temperature_coefficient",
        "demagnetization_curve",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent NGYC N42 property must remain refused: {refused_property}"
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
fn g3_cli_compiles_committed_jinshan_n42_pristine_temperature_endpoints() {
    let manifest = workspace_path(JINSHAN_N42_PRISTINE_TEMPERATURE_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Jinshan N42 temperature seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("jinshan-n42-temperature-first.fsmatpk");
    let second_path = directory.join("jinshan-n42-temperature-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Jinshan N42 temperature seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Jinshan N42 temperature seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Jinshan N42 temperature decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Jinshan N42 temperature pack");
    let second_bytes = fs::read(second_path).expect("read second Jinshan N42 temperature pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Jinshan N42 temperature pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode Jinshan N42 temperature pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Jinshan N42 temperature pack identity");

    assert_eq!(
        decoded.pack_id(),
        "jinshan-n42-pristine-sintered-temperature-endpoints"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 8);
    assert!(decoded.joint_statistics().is_empty());

    let endpoint_expectations = [
        (
            "remanent_flux_density",
            25.0_f64,
            12.75_f64 * 1.0e3 * 1.0e-4,
            Dims([0, 1, -2, 0, -1, 0]),
        ),
        (
            "remanent_flux_density",
            120.0_f64,
            11.18_f64 * 1.0e3 * 1.0e-4,
            Dims([0, 1, -2, 0, -1, 0]),
        ),
        (
            "intrinsic_coercive_field_strength",
            25.0_f64,
            12.07_f64 * 1.0e3 * (1.0e3 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 0, 0, 0, 1, 0]),
        ),
        (
            "intrinsic_coercive_field_strength",
            120.0_f64,
            5.17_f64 * 1.0e3 * (1.0e3 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 0, 0, 0, 1, 0]),
        ),
        (
            "maximum_magnetic_energy_product",
            25.0_f64,
            40.14_f64 * (100_000.0 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 1, -2, 0, 0, 0]),
        ),
        (
            "maximum_magnetic_energy_product",
            120.0_f64,
            29.29_f64 * (100_000.0 / (4.0 * std::f64::consts::PI)),
            Dims([-1, 1, -2, 0, 0, 0]),
        ),
    ];
    let endpoint_observation = endpoint_expectations
        .iter()
        .map(
            |(property, source_temperature_c, expected_value, expected_dims)| {
                let temperature_k = source_temperature_c + 273.15;
                let claims = decoded.claims().claims_for(property);
                let (_, claim) = claims
                    .into_iter()
                    .find(|(_, claim)| {
                        claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
                    })
                    .unwrap_or_else(|| {
                        panic!("missing Jinshan N42 {property} at {source_temperature_c} degC")
                    });
                let PropertyValue::Scalar { value, dims } = &claim.value else {
                    panic!("Jinshan N42 {property} endpoint was not scalar");
                };
                assert_eq!(dims, expected_dims);
                let scale = expected_value.abs().max(1.0e-12);
                assert!((*value - expected_value).abs() / scale <= 2.0e-15);
                assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
                assert!(
                    claim
                        .provenance
                        .source
                        .contains("10.1016/j.jmrt.2024.12.235")
                );
                assert!(claim.provenance.source.contains("[source:primary]"));
                assert_eq!(
                    claim.validity.bound("source_instrument_nim_6500c"),
                    Some((1.0, 1.0))
                );
                assert_eq!(
                    claim
                        .validity
                        .bound("source_authors_heat_treatment_applied"),
                    Some((0.0, 0.0))
                );
                for missing_axis in [
                    "source_supplier_production_lot_known",
                    "source_composition_known",
                    "source_temperature_control_method_known",
                ] {
                    assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                }
                claim.observations[0]
            },
        )
        .collect::<Vec<_>>();
    assert!(
        endpoint_observation
            .windows(2)
            .all(|pair| pair[0] == pair[1]),
        "all Table 1 endpoints must retain one shared observation"
    );
    let endpoint_observation = decoded
        .claims()
        .observation(endpoint_observation[0])
        .expect("Jinshan N42 endpoint observation remains linked");
    assert_eq!(
        endpoint_observation.specimen,
        "jinshan-magnetic-materials-commercial-n42-sintered-ndfeb-pristine-wire-cut-10x10x6-mm-lot-unspecified"
    );
    assert!(endpoint_observation.method.contains("NIM 6500C"));
    assert!(
        endpoint_observation
            .caveats
            .contains("calls the 10 mm x 10 mm x 6 mm pieces cubes")
    );
    assert!(
        endpoint_observation
            .caveats
            .contains("no curve points are digitized")
    );

    for (property, expected_value) in [
        ("remanence_temperature_coefficient", -0.00129_f64),
        ("coercivity_temperature_coefficient", -0.00602_f64),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing Jinshan N42 {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Jinshan N42 {property} was not scalar");
        };
        assert_eq!(*dims, Dims([0, 0, 0, -1, 0, 0]));
        assert_eq!(*value, expected_value);
        assert_eq!(
            claim.validity.bound("temperature"),
            Some((25.0 + 273.15, 120.0 + 273.15))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_coefficient_uses_25c_and_120c_endpoints"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_rounded_endpoints_reproduce_printed_coefficient_exactly"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("Jinshan N42 coefficient observation remains linked");
        assert!(
            observation
                .caveats
                .contains("do not reproduce both printed coefficients exactly")
        );
        assert!(
            observation
                .caveats
                .contains("not continuous constitutive laws")
        );
    }

    for refused_property in [
        "coercive_field_strength",
        "recoil_relative_permeability",
        "demagnetization_curve",
        "irreversible_demagnetization_loss_boundary",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent Jinshan N42 property must remain refused: {refused_property}"
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
fn g3_cli_compiles_committed_naca_tn_2680_isooctane_flame_speed_seed() {
    let manifest = workspace_path(NACA_TN_2680_ISOOCTANE_FLAME_SPEED_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NACA TN 2680 iso-octane flame-speed seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("naca-tn-2680-isooctane-first.fsmatpk");
    let second_path = directory.join("naca-tn-2680-isooctane-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NACA TN 2680 iso-octane seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NACA TN 2680 iso-octane seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NACA TN 2680 iso-octane decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first NACA TN 2680 iso-octane pack");
    let second_bytes = fs::read(second_path).expect("read second NACA TN 2680 iso-octane pack");
    assert_eq!(
        first_bytes, second_bytes,
        "NACA TN 2680 iso-octane pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode NACA TN 2680 iso-octane pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify NACA TN 2680 iso-octane pack identity");

    assert_eq!(
        decoded.pack_id(),
        "naca-tn-2680-2-2-4-trimethylpentane-flame-speed"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("Work of the US Gov. Public Use Permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 16);
    assert!(decoded.joint_statistics().is_empty());

    let purity = decoded
        .claims()
        .claims_for("minimum_reported_fuel_mole_fraction_purity");
    assert_eq!(purity.len(), 1);
    let (_, purity_claim) = purity[0];
    let PropertyValue::Scalar { value, dims } = &purity_claim.value else {
        panic!("NACA TN 2680 fuel minimum purity was not scalar");
    };
    let expected_purity = 99.6_f64 * 0.01;
    assert!((*value - expected_purity).abs() / expected_purity <= 2.0e-15);
    assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
    assert_eq!(
        purity_claim
            .validity
            .bound("source_fuel_supplier_identity_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        purity_claim.validity.bound("source_fuel_lot_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(
        purity_claim.validity.bound("source_exact_assay_known"),
        Some((0.0, 0.0))
    );
    assert_eq!(purity_claim.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(purity_claim.provenance.license, NASA_SEED_LICENSE);
    assert!(
        purity_claim
            .provenance
            .source
            .contains("NACA Technical Note 2680")
    );
    let purity_observation = decoded
        .claims()
        .observation(purity_claim.observations[0])
        .expect("NACA TN 2680 purity observation remains linked");
    assert!(
        purity_observation
            .method
            .contains("minimum-purity statement")
    );
    assert!(purity_observation.caveats.contains("lower-bound statement"));

    let flame_claims = decoded.claims().claims_for("maximum_laminar_flame_speed");
    assert_eq!(flame_claims.len(), 15);
    let expected_rows = [
        (311.0, 0.210, 1000.0, 1.256, 34.6),
        (311.0, 0.250, 1000.0, 0.838, 52.1),
        (311.0, 0.294, 1600.0, 0.838, 72.2),
        (311.0, 0.294, 900.0, 0.617, 67.2),
        (311.0, 0.347, 1200.0, 0.617, 89.1),
        (311.0, 0.496, 1800.0, 0.297, 152.2),
        (367.0, 0.210, 1000.0, 1.256, 44.8),
        (422.0, 0.210, 1000.0, 1.256, 56.1),
        (422.0, 0.210, 1000.0, 1.256, 59.0),
        (422.0, 0.210, 700.0, 0.838, 57.4),
        (422.0, 0.250, 1400.0, 0.838, 83.1),
        (422.0, 0.294, 900.0, 0.617, 108.0),
        (422.0, 0.294, 900.0, 0.617, 102.1),
        (422.0, 0.347, 1400.0, 0.617, 138.0),
        (422.0, 0.496, 1800.0, 0.297, 229.9),
    ];

    for (temperature, oxygen_fraction, reynolds_number, diameter_cm, speed_cm_per_s) in
        expected_rows
    {
        let diameter_m = diameter_cm * 0.01;
        let speed_m_per_s = speed_cm_per_s * 0.01;
        let matching = flame_claims
            .iter()
            .filter(|(_, claim)| {
                claim.validity.bound("initial_mixture_temperature")
                    == Some((temperature, temperature))
                    && claim
                        .validity
                        .bound("oxygen_mole_fraction_in_oxygen_nitrogen")
                        == Some((oxygen_fraction, oxygen_fraction))
                    && claim.validity.bound("stream_flow_reynolds_number")
                        == Some((reynolds_number, reynolds_number))
                    && claim.validity.bound("burner_inside_diameter")
                        == Some((diameter_m, diameter_m))
                    && matches!(
                        &claim.value,
                        PropertyValue::Scalar { value, .. } if *value == speed_m_per_s
                    )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            matching.len(),
            1,
            "expected one NACA Table I row at T={temperature}, O2={oxygen_fraction}, Re={reynolds_number}, diameter={diameter_cm} cm, speed={speed_cm_per_s} cm/s"
        );
        let (_, claim) = matching[0];
        let PropertyValue::Scalar { dims, .. } = &claim.value else {
            unreachable!("matching predicate admitted only scalar flame speeds");
        };
        assert_eq!(*dims, Dims([1, 0, -1, 0, 0, 0]));
        assert_eq!(
            claim.validity.bound("average_atmospheric_pressure"),
            Some((99.2_f64 * 1.0e3, 99.2_f64 * 1.0e3))
        );
        assert_eq!(
            claim.validity.bound("source_pressure_per_row_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_equivalence_ratio_at_maximum_exact_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_maximum_equivalence_ratio_lower_bound"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_maximum_equivalence_ratio_upper_bound"),
            Some((1.1, 1.1))
        );
        assert_eq!(
            claim.validity.bound("source_fuel_supplier_identity_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_fuel_lot_known"),
            Some((0.0, 0.0))
        );
        if oxygen_fraction == 0.210 {
            assert_eq!(
                claim
                    .validity
                    .bound("source_oxidizer_analysis_half_width_known"),
                Some((0.0, 0.0))
            );
            assert_eq!(
                claim
                    .validity
                    .bound("source_oxidizer_analysis_absolute_half_width"),
                None
            );
        } else {
            assert_eq!(
                claim
                    .validity
                    .bound("source_oxidizer_analysis_absolute_half_width"),
                Some((0.1_f64 * 0.01, 0.1_f64 * 0.01))
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NTRS 19930083861"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("NACA TN 2680 flame-speed observation remains linked");
        assert!(observation.method.contains("total-area method"));
        assert!(
            observation
                .caveats
                .contains("not geometry-free bulk-material constants")
        );
        assert!(observation.caveats.contains("Repeated rows"));
    }

    for refused_property in [
        "density",
        "dynamic_viscosity",
        "surface_tension",
        "specific_heat_capacity",
        "heat_of_vaporization",
        "vapor_pressure",
        "research_octane_number",
        "empirical_maximum_flame_speed_fit",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent or model-only NACA TN 2680 property must remain refused: {refused_property}"
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
fn g3_cli_compiles_committed_face_g_cdtrf_g_2023_v1_surrogate_seed() {
    let manifest = workspace_path(FACE_G_CDTRF_G_2023_V1_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed FACE G CDTRF-G 2023 v1 surrogate seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("face-g-cdtrf-g-2023-v1-first.fsmatpk");
    let second_path = directory.join("face-g-cdtrf-g-2023-v1-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first FACE G CDTRF-G seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second FACE G CDTRF-G seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "FACE G CDTRF-G decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first FACE G CDTRF-G pack");
    let second_bytes = fs::read(second_path).expect("read second FACE G CDTRF-G pack");
    assert_eq!(first_bytes, second_bytes, "FACE G CDTRF-G pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode FACE G CDTRF-G pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify FACE G CDTRF-G pack identity");

    assert_eq!(decoded.pack_id(), "face-g-cdtrf-g-2023-v1");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(decoded.redistribution_terms().contains("Creative Commons"));
    assert_eq!(decoded.claims().claim_count(), 7);
    assert!(decoded.joint_statistics().is_empty());

    let expected_components = [
        ("isooctane_component_volume_fraction", 23.75),
        ("n_heptane_component_volume_fraction", 19.0),
        ("toluene_component_volume_fraction", 42.75),
        ("diisobutylene_component_volume_fraction", 9.5),
        ("cyclohexane_component_volume_fraction", 5.0),
    ];
    let mut fraction_sum = 0.0;
    for (property, source_percent) in expected_components {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique CDTRF-G {property}");
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("CDTRF-G {property} was not scalar");
        };
        let expected_fraction = source_percent * 0.01;
        assert!((*value - expected_fraction).abs() / expected_fraction <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        fraction_sum += *value;
        assert_eq!(
            claim
                .validity
                .bound("source_composition_basis_is_volume_fraction"),
            Some((1.0, 1.0))
        );
        for missing_axis in [
            "source_component_supplier_known",
            "source_component_lot_known",
            "source_component_purity_known",
            "source_mixing_temperature_known",
            "source_mixing_pressure_known",
            "source_volume_contraction_treatment_known",
        ] {
            assert_eq!(
                claim.validity.bound(missing_axis),
                Some((0.0, 0.0)),
                "CDTRF-G {property} must retain missing axis {missing_axis}"
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        assert!(
            claim
                .provenance
                .source
                .contains("10.3390/molecules28114273")
        );
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("CDTRF-G composition observation remains linked");
        assert_eq!(
            observation.specimen,
            "face-g-targeted-cdtrf-g-2023-v1-source-component-lots-and-purities-unspecified"
        );
        assert!(observation.method.contains("Table 2 CDTRF-G"));
        assert!(observation.caveats.contains("sum to exactly 100 percent"));
        assert!(observation.caveats.contains("molar basis"));
    }
    assert!((fraction_sum - 1.0_f64).abs() <= 2.0e-15);

    let ron_claims = decoded
        .claims()
        .claims_for("reported_calculated_research_octane_number");
    assert_eq!(ron_claims.len(), 2);
    assert_ne!(
        ron_claims[0].1.observations[0], ron_claims[1].1.observations[0],
        "conflicting CDTRF-G RON prints must retain distinct observations"
    );
    let mut ron_values = ron_claims
        .iter()
        .map(|(_, claim)| match &claim.value {
            PropertyValue::Scalar { value, dims } => {
                assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
                *value
            }
            PropertyValue::Curve { .. } => panic!("CDTRF-G RON was not scalar"),
        })
        .collect::<Vec<_>>();
    ron_values.sort_by(f64::total_cmp);
    assert_eq!(ron_values, vec![93.9, 94.0]);
    for (_, claim) in ron_claims {
        assert_eq!(
            claim.validity.bound("source_octane_measurement_performed"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_octane_calculation_basis_unambiguous"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_component_purity_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, CC_BY_4_0_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("CDTRF-G RON observation remains linked");
        assert!(observation.caveats.contains("Table 2"));
        assert!(observation.caveats.contains("Table 7"));
        assert!(observation.caveats.contains("separate"));
    }

    for refused_property in [
        "density",
        "dynamic_viscosity",
        "surface_tension",
        "specific_heat_capacity",
        "heat_of_vaporization",
        "vapor_pressure",
        "motor_octane_number",
        "laminar_flame_speed",
        "ignition_delay_time",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent or model-only CDTRF-G property must remain refused: {refused_property}"
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
fn g3_cli_compiles_committed_wo2018_formulation_8_5w30_seed() {
    let manifest = workspace_path(WO2018_125520_FORMULATION_8_5W30_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed WO 2018/125520 Formulation 8 5W30 seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("wo2018-125520-formulation-8-5w30-first.fsmatpk");
    let second_path = directory.join("wo2018-125520-formulation-8-5w30-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first WO 2018/125520 Formulation 8 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second WO 2018/125520 Formulation 8 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "WO 2018/125520 Formulation 8 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first WO 2018 Formulation 8 pack");
    let second_bytes = fs::read(second_path).expect("read second WO 2018 Formulation 8 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "WO 2018/125520 Formulation 8 pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode WO 2018/125520 Formulation 8 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify WO 2018/125520 Formulation 8 pack identity");

    assert_eq!(decoded.pack_id(), "wo2018-125520-formulation-8-5w30");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("not a patent-practice or trademark license")
    );
    assert_eq!(decoded.claims().claim_count(), 12);
    assert!(decoded.joint_statistics().is_empty());

    let expected_components = [
        ("spectrasyn_4_component_mass_fraction", 60.0),
        ("synesstic_5_component_mass_fraction", 10.0),
        ("spectrasyn_elite_150_component_mass_fraction", 18.0),
        ("infineum_p6003_component_mass_fraction", 12.0),
    ];
    let mut fraction_sum = 0.0;
    let mut composition_observation = None;
    for (property, source_percent) in expected_components {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique Formulation 8 {property}");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Formulation 8 {property} was not scalar");
        };
        let expected_fraction = source_percent * 0.01;
        assert!((*value - expected_fraction).abs() <= 2.0e-15);
        assert_eq!(*dims, Dims([0, 0, 0, 0, 0, 0]));
        fraction_sum += *value;
        assert_eq!(
            claim
                .validity
                .bound("source_composition_basis_is_mass_fraction"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim.validity.bound("source_formulation_number"),
            Some((8.0, 8.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_component_commercial_identifier_known"),
            Some((1.0, 1.0))
        );
        for missing_axis in [
            "source_component_lot_known",
            "source_component_detailed_chemistry_known",
            "source_final_blend_protocol_known",
            "source_patent_practice_license_granted",
        ] {
            assert_eq!(
                claim.validity.bound(missing_axis),
                Some((0.0, 0.0)),
                "Formulation 8 {property} must retain missing axis {missing_axis}"
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, USPTO_PATENT_TEXT_LICENSE);
        assert!(claim.provenance.source.contains("WO 2018/125520 A1"));
        assert!(claim.provenance.source.contains("US 2018/0179462 A1"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        assert_eq!(id.0, claim.content_hash());
        match composition_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => composition_observation = Some(claim.observations[0]),
        }
    }
    assert!((fraction_sum - 1.0_f64).abs() <= 2.0e-15);

    let composition = decoded
        .claims()
        .observation(composition_observation.expect("composition observation exists"))
        .expect("Formulation 8 composition observation remains linked");
    assert_eq!(
        composition.specimen,
        "wo2018-125520-table-ix-formulation-8-source-products-lots-unspecified"
    );
    assert!(composition.method.contains("Table IX Formulation 8"));
    assert!(composition.caveats.contains("sum to exactly 100.00 wt%"));
    assert!(composition.caveats.contains("not present-day fungible"));
    assert!(composition.caveats.contains("without implying endorsement"));

    let kinematic_viscosity_dims = Dims([2, 0, -1, 0, 0, 0]);
    let viscosity_claims = decoded.claims().claims_for("kinematic_viscosity");
    assert_eq!(viscosity_claims.len(), 2);
    for (source_temperature_c, source_mm2_per_s) in [(40.0, 61.49), (100.0, 10.62)] {
        let temperature_k = source_temperature_c + 273.15;
        let mut matches = viscosity_claims.iter().copied().filter(|(_, claim)| {
            claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
        });
        let (_, claim) = matches.next().unwrap_or_else(|| {
            panic!("missing Formulation 8 viscosity at {source_temperature_c} degC")
        });
        assert!(matches.next().is_none());
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Formulation 8 viscosity was not scalar");
        };
        assert_eq!(*dims, kinematic_viscosity_dims);
        let expected_m2_per_s = source_mm2_per_s * 1.0e-6;
        assert!((*value - expected_m2_per_s).abs() / expected_m2_per_s <= 2.0e-15);
    }

    let dynamic_viscosity_dims = Dims([-1, 1, -1, 0, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let temperature_dims = Dims([0, 0, 0, 1, 0, 0]);
    let expected_unique = [
        ("viscosity_index_scale_reading", 164.0, dimensionless, None),
        (
            "pour_point_temperature",
            -66.0 + 273.15,
            temperature_dims,
            None,
        ),
        (
            "cold_cranking_simulator_dynamic_viscosity",
            4_886.0 * 1.0e-3,
            dynamic_viscosity_dims,
            Some(-30.0 + 273.15),
        ),
        (
            "mini_rotary_viscometer_dynamic_viscosity",
            10_782.0 * 1.0e-3,
            dynamic_viscosity_dims,
            Some(-35.0 + 273.15),
        ),
        (
            "high_temperature_high_shear_dynamic_viscosity",
            3.395 * 1.0e-3,
            dynamic_viscosity_dims,
            Some(150.0 + 273.15),
        ),
        (
            "noack_mass_loss_fraction",
            9.2 * 0.01,
            dimensionless,
            Some(250.0 + 273.15),
        ),
    ];
    let mut performance_observation = None;
    for (property, expected_value, expected_dims, validity_temperature) in expected_unique {
        let expected_value: f64 = expected_value;
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "missing unique Formulation 8 {property}");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Formulation 8 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let scale = expected_value.abs().max(1.0e-12);
        assert!((*value - expected_value).abs() / scale <= 2.0e-15);
        match validity_temperature {
            Some(temperature_k) => assert_eq!(
                claim.validity.bound("temperature"),
                Some((temperature_k, temperature_k))
            ),
            None => assert_eq!(claim.validity.bound("temperature"), None),
        }
        assert_eq!(
            claim.validity.bound("source_formulation_number"),
            Some((8.0, 8.0))
        );
        assert_eq!(
            claim.validity.bound("source_viscosity_grade_is_5w30"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_patent_practice_license_granted"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_test_method_edition_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(
            claim.validity.bound("source_repeat_count_known"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, USPTO_PATENT_TEXT_LICENSE);
        assert_eq!(id.0, claim.content_hash());
        match performance_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => performance_observation = Some(claim.observations[0]),
        }
    }
    assert_eq!(
        decoded.claims().claims_for("noack_mass_loss_fraction")[0]
            .1
            .validity
            .bound("source_test_duration_known"),
        Some((0.0, 0.0))
    );

    for (_, claim) in viscosity_claims {
        assert_eq!(
            claim.validity.bound("source_formulation_number"),
            Some((8.0, 8.0))
        );
        assert_eq!(
            claim.validity.bound("source_viscosity_grade_is_5w30"),
            Some((1.0, 1.0))
        );
        assert_eq!(
            claim
                .validity
                .bound("source_patent_practice_license_granted"),
            Some((0.0, 0.0))
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, USPTO_PATENT_TEXT_LICENSE);
        match performance_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => performance_observation = Some(claim.observations[0]),
        }
    }

    let performance = decoded
        .claims()
        .observation(performance_observation.expect("performance observation exists"))
        .expect("Formulation 8 performance observation remains linked");
    assert!(performance.method.contains("5W30 property-row"));
    assert!(performance.caveats.contains("absences, not zero-valued"));
    assert!(performance.caveats.contains("do not generalize"));

    assert!(decoded.claims().claims_for("dynamic_viscosity").is_empty());
    assert!(decoded.claims().claims_for("density").is_empty());
    assert!(decoded.claims().claims_for("total_base_number").is_empty());
    assert!(
        decoded
            .claims()
            .claims_for("flash_point_temperature")
            .is_empty()
    );

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
fn g3_cli_compiles_committed_nasa_uam_insulation_stack_constituents() {
    let seeds = [
        (
            NASA_UAM_MW16C_POLYIMIDE_WIRE_SEED_MANIFEST,
            "nasa-uam-mw16c-polyimide-magnet-wire",
            2_usize,
        ),
        (
            NASA_UAM_NOMEX_410_SLOT_LINER_SEED_MANIFEST,
            "nasa-uam-nomex-410-slot-liner",
            1_usize,
        ),
        (
            NASA_UAM_COOLTHERM_EP2000_SEED_MANIFEST,
            "nasa-uam-cooltherm-ep2000-180c-cure",
            2_usize,
        ),
    ];
    let directory = fixture_dir();

    for (manifest_relative, expected_pack_id, expected_claim_count) in seeds {
        let manifest = workspace_path(manifest_relative);
        assert!(
            manifest.is_file(),
            "committed NASA UAM insulation seed manifest is missing: {manifest_relative}"
        );
        let first_path = directory.join(format!("{expected_pack_id}-first.fsmatpk"));
        let second_path = directory.join(format!("{expected_pack_id}-second.fsmatpk"));

        let first = run_compiler(&manifest, &first_path);
        let second = run_compiler(&manifest, &second_path);
        assert!(
            first.status.success(),
            "first {expected_pack_id} compilation failed: {}",
            String::from_utf8_lossy(&first.stderr)
        );
        assert!(
            second.status.success(),
            "second {expected_pack_id} compilation failed: {}",
            String::from_utf8_lossy(&second.stderr)
        );
        assert_eq!(
            first.stdout, second.stdout,
            "{expected_pack_id} decision stream moved"
        );
        assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

        let first_bytes = fs::read(first_path).expect("read first NASA insulation pack");
        let second_bytes = fs::read(second_path).expect("read second NASA insulation pack");
        assert_eq!(
            first_bytes, second_bytes,
            "{expected_pack_id} pack bytes moved"
        );
        let decoded =
            NormalizedPack::from_bytes(&first_bytes).expect("decode NASA insulation pack");
        let pack_hash = decoded.content_hash();
        let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
            .expect("verify NASA insulation pack identity");

        assert_eq!(decoded.pack_id(), expected_pack_id);
        assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
        assert!(
            decoded
                .redistribution_terms()
                .contains("government public use permitted")
        );
        assert_eq!(decoded.claims().claim_count(), expected_claim_count);
        assert!(decoded.joint_statistics().is_empty());

        match expected_pack_id {
            "nasa-uam-mw16c-polyimide-magnet-wire" => {
                let temperature_claims = decoded
                    .claims()
                    .claims_for("thermal_endurance_reference_temperature");
                let duration_claims = decoded
                    .claims()
                    .claims_for("thermal_endurance_reference_duration");
                assert_eq!(temperature_claims.len(), 1);
                assert_eq!(duration_claims.len(), 1);

                let (temperature_id, temperature_claim) = temperature_claims[0];
                let PropertyValue::Scalar {
                    value: temperature,
                    dims: temperature_dims,
                } = &temperature_claim.value
                else {
                    panic!("MW-16C thermal-endurance temperature was not scalar");
                };
                assert_eq!(*temperature_dims, Dims([0, 0, 0, 1, 0, 0]));
                assert!((*temperature - (240.0 + 273.15)).abs() <= 1.0e-12);
                assert_eq!(
                    temperature_claim.validity.bound("reference_duration"),
                    Some((20_000.0 * 3_600.0, 20_000.0 * 3_600.0))
                );

                let (duration_id, duration_claim) = duration_claims[0];
                let PropertyValue::Scalar {
                    value: duration,
                    dims: duration_dims,
                } = &duration_claim.value
                else {
                    panic!("MW-16C thermal-endurance duration was not scalar");
                };
                assert_eq!(*duration_dims, Dims([0, 0, 1, 0, 0, 0]));
                assert_eq!(*duration, 20_000.0 * 3_600.0);
                assert_eq!(
                    duration_claim.validity.bound("reference_temperature"),
                    Some((240.0 + 273.15, 240.0 + 273.15))
                );

                for (id, claim) in [
                    (temperature_id, temperature_claim),
                    (duration_id, duration_claim),
                ] {
                    for required_axis in [
                        "source_wire_spec_nema_mw16c",
                        "source_nema_mw1000_2003",
                        "source_test_standard_astm_d2307_2013",
                    ] {
                        assert_eq!(claim.validity.bound(required_axis), Some((1.0, 1.0)));
                    }
                    for missing_axis in [
                        "source_wire_vendor_known",
                        "source_wire_lot_known",
                        "source_thermal_endurance_raw_data_available",
                    ] {
                        assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                    }
                    assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                    assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
                    assert!(claim.provenance.source.contains("NTRS 20240007451"));
                    assert!(claim.provenance.source.contains("[source:primary]"));
                    assert_eq!(id.0, claim.content_hash());
                }
                assert_eq!(temperature_claim.observations, duration_claim.observations);
                let observation = decoded
                    .claims()
                    .observation(temperature_claim.observations[0])
                    .expect("MW-16C observation remains linked");
                assert!(observation.method.contains("ASTM D2307-2013"));
                assert!(
                    observation
                        .caveats
                        .contains("cross-bound classification basis")
                );
                assert!(observation.caveats.contains("not an Arrhenius law"));
                assert!(observation.caveats.contains("different unspecified spool"));
            }
            "nasa-uam-nomex-410-slot-liner" => {
                let claims = decoded.claims().claims_for("selected_slot_liner_thickness");
                assert_eq!(claims.len(), 1);
                let (id, claim) = claims[0];
                let PropertyValue::Scalar { value, dims } = &claim.value else {
                    panic!("Nomex 410 selected thickness was not scalar");
                };
                assert_eq!(*dims, Dims([1, 0, 0, 0, 0, 0]));
                assert!((*value - 0.08e-3).abs() <= 1.0e-15);
                assert_eq!(
                    claim.validity.bound("source_product_is_dupont_nomex_410"),
                    Some((1.0, 1.0))
                );
                assert_eq!(
                    claim.validity.bound("source_material_is_aramid_paper"),
                    Some((1.0, 1.0))
                );
                for missing_axis in [
                    "source_product_lot_known",
                    "source_thickness_measurement_method_known",
                    "source_moisture_condition_known",
                ] {
                    assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                }
                assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
                assert_eq!(id.0, claim.content_hash());
                let observation = decoded
                    .claims()
                    .observation(claim.observations[0])
                    .expect("Nomex 410 observation remains linked");
                assert!(observation.method.contains("component-selection statement"));
                assert!(
                    observation
                        .caveats
                        .contains("not a generic Nomex 410 design allowable")
                );
            }
            "nasa-uam-cooltherm-ep2000-180c-cure" => {
                let completed_claims = decoded
                    .claims()
                    .claims_for("highest_completed_post_cure_temperature");
                let omitted_claims = decoded
                    .claims()
                    .claims_for("manufacturer_recommended_final_cure_temperature");
                assert_eq!(completed_claims.len(), 1);
                assert_eq!(omitted_claims.len(), 1);

                for (claims, source_temperature_c, step_completed) in [
                    (&completed_claims, 180.0, 1.0),
                    (&omitted_claims, 210.0, 0.0),
                ] {
                    let (id, claim) = claims[0];
                    let PropertyValue::Scalar { value, dims } = &claim.value else {
                        panic!("CoolTherm EP-2000 cure temperature was not scalar");
                    };
                    assert_eq!(*dims, Dims([0, 0, 0, 1, 0, 0]));
                    assert!((*value - (source_temperature_c + 273.15)).abs() <= 1.0e-12);
                    assert_eq!(
                        claim
                            .validity
                            .bound("source_product_is_parker_lord_cooltherm_ep2000"),
                        Some((1.0, 1.0))
                    );
                    assert_eq!(
                        claim.validity.bound("source_step_completed"),
                        Some((step_completed, step_completed))
                    );
                    for missing_axis in [
                        "source_epoxy_lot_known",
                        "source_cure_hold_duration_known",
                        "source_degree_of_cure_known",
                    ] {
                        assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
                    }
                    assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
                    assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
                    assert_eq!(id.0, claim.content_hash());
                }
                assert_eq!(
                    completed_claims[0].1.observations,
                    omitted_claims[0].1.observations
                );
                let observation = decoded
                    .claims()
                    .observation(completed_claims[0].1.observations[0])
                    .expect("CoolTherm EP-2000 observation remains linked");
                assert!(observation.caveats.contains("intentionally not completed"));
                assert!(
                    observation
                        .caveats
                        .contains("deliberately incomplete source process")
                );
            }
            unexpected => panic!("unexpected NASA insulation pack {unexpected}"),
        }

        for refused_property in [
            "partial_discharge_inception_voltage",
            "dielectric_strength",
            "thermal_conductivity",
            "service_life",
            "arrhenius_activation_energy",
        ] {
            assert!(
                decoded.claims().claims_for(refused_property).is_empty(),
                "assembly- or source-absent property must stay refused for {expected_pack_id}: {refused_property}"
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
fn g3_cli_compiles_committed_aisi_9310_cvm_carburized_gear_seed() {
    let manifest = workspace_path(AISI_9310_CVM_CARBURIZED_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed AISI 9310 CVM carburized seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("aisi-9310-cvm-carburized-first.fsmatpk");
    let second_path = directory.join("aisi-9310-cvm-carburized-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first AISI 9310 CVM carburized seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second AISI 9310 CVM carburized seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "AISI 9310 CVM carburized decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first AISI 9310 pack");
    let second_bytes = fs::read(second_path).expect("read second AISI 9310 pack");
    assert_eq!(first_bytes, second_bytes, "AISI 9310 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode AISI 9310 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify AISI 9310 pack identity");

    assert_eq!(decoded.pack_id(), "aisi-9310-cvm-carburized-nasa-tm-104352");
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use is permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 13);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    for (property, source_percent) in [
        ("carbon_mass_fraction", 0.10),
        ("manganese_mass_fraction", 0.63),
        ("silicon_mass_fraction", 0.27),
        ("nickel_mass_fraction", 3.22),
        ("chromium_mass_fraction", 1.21),
        ("molybdenum_mass_fraction", 0.12),
        ("copper_mass_fraction", 0.13),
        ("phosphorus_mass_fraction", 0.005),
        ("sulfur_mass_fraction", 0.005),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one AISI 9310 {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 9310 {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        let relative_error = (*value - expected_value).abs() / expected_value;
        assert!(
            relative_error <= 2.0e-15,
            "AISI 9310 {property} moved by {relative_error:e} relative"
        );
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-TM-104352"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 9310 chemistry observation remains linked");
        assert_eq!(
            observation.specimen,
            "AISI-9310-CVM-single-lot-single-heat-28-tooth-spur-gear-NASA-TM-104352"
        );
        assert!(observation.method.contains("Table I nominal composition"));
        assert!(observation.caveats.contains("Nominal grade chemistry"));
        assert!(observation.caveats.contains("balance iron"));
        assert!(observation.caveats.contains("not inferred"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    let case_claims = decoded.claims().claims_for("case_rockwell_c_scale_reading");
    assert_eq!(
        case_claims.len(),
        2,
        "the report's conflicting C58 and C60 case statements must both survive"
    );
    let mut case_values = Vec::with_capacity(2);
    for (id, claim) in case_claims {
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("AISI 9310 case hardness was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        assert!(claim.validity.bounds().is_empty());
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("AISI 9310 case-hardness observation remains linked");
        match *value {
            58.0 => {
                assert!(observation.method.contains("Test Materials detailed"));
                assert!(observation.caveats.contains("carburize 1172 K for 8 h"));
                assert!(observation.caveats.contains("austenitize 1117 K for 2.5 h"));
                assert!(
                    observation
                        .caveats
                        .contains("subzero treat 180 K for 3.5 h")
                );
                assert!(observation.caveats.contains("double temper 450 K"));
                assert!(observation.caveats.contains("stress relieve 450 K for 2 h"));
                assert!(observation.caveats.contains("conflicts"));
            }
            60.0 => {
                assert!(observation.method.contains("abstract and summary"));
                assert!(observation.caveats.contains("conflicts"));
                assert!(observation.caveats.contains("not averaged or selected"));
            }
            other => panic!("unexpected AISI 9310 case-hardness value {other}"),
        }
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
        case_values.push(*value);
    }
    case_values.sort_by(f64::total_cmp);
    assert_eq!(case_values, [58.0, 60.0]);

    let core_claims = decoded.claims().claims_for("core_rockwell_c_scale_reading");
    assert_eq!(core_claims.len(), 1);
    let (core_id, core_claim) = core_claims[0];
    let PropertyValue::Scalar {
        value: core_value,
        dims: core_dims,
    } = &core_claim.value
    else {
        panic!("AISI 9310 core hardness was not scalar");
    };
    assert_eq!(*core_dims, dimensionless);
    assert_eq!(*core_value, 40.0);
    assert!(core_claim.validity.bounds().is_empty());
    assert_eq!(core_claim.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(core_claim.observations.len(), 1);

    let depth_claims = decoded.claims().claims_for("carburized_case_depth");
    assert_eq!(depth_claims.len(), 1);
    let (depth_id, depth_claim) = depth_claims[0];
    let PropertyValue::Scalar {
        value: depth_value,
        dims: depth_dims,
    } = &depth_claim.value
    else {
        panic!("AISI 9310 carburized case depth was not scalar");
    };
    assert_eq!(*depth_dims, Dims([1, 0, 0, 0, 0, 0]));
    let expected_depth_m = 0.97e-3;
    let relative_depth_error = (*depth_value - expected_depth_m).abs() / expected_depth_m;
    assert!(relative_depth_error <= 2.0e-15);
    assert!(depth_claim.validity.bounds().is_empty());
    assert_eq!(depth_claim.uncertainty, UncertaintyModel::Unstated);
    assert_eq!(depth_claim.observations.len(), 1);
    assert_eq!(core_claim.observations, depth_claim.observations);

    let detailed_observation = decoded
        .claims()
        .observation(core_claim.observations[0])
        .expect("AISI 9310 detailed gear observation remains linked");
    assert!(detailed_observation.method.contains("case/core hardness"));
    assert!(detailed_observation.method.contains("case-depth"));
    assert!(
        detailed_observation
            .caveats
            .contains("One lot from one CVM heat")
    );
    assert!(detailed_observation.caveats.contains("replicate count"));
    assert_eq!(
        core_claim.observations[0].0,
        detailed_observation.content_hash()
    );
    assert_eq!(core_id.0, core_claim.content_hash());
    assert_eq!(depth_id.0, depth_claim.content_hash());

    // G3 plausibility only: NASA SP-410 (NTRS 19750018303) reports a different
    // VAR AISI 9310 gear lot at nominal C62 case, C45 core, and 1 mm case depth.
    // These checks bound transcription-scale agreement without fusing the lots.
    let independent_case_hardness: f64 = 62.0;
    let independent_core_hardness: f64 = 45.0;
    let independent_case_depth_m: f64 = 1.0e-3;
    assert!((58.0 - independent_case_hardness).abs() <= 4.0);
    assert!((*core_value - independent_core_hardness).abs() <= 5.0);
    assert!((*depth_value - independent_case_depth_m).abs() / independent_case_depth_m <= 0.031);

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
fn g3_cli_compiles_committed_napc_gear_oils_without_fusing_batches() {
    let temperature_dims = Dims([0, 0, 0, 1, 0, 0]);
    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    type ExpectedClaim = (&'static str, &'static str, f64, Dims, Option<f64>);
    let seeds: [(&str, &str, &str, &[ExpectedClaim]); 2] = [
        (
            NAPC_PE_5_L_1274_GEAR_OIL_SEED_MANIFEST,
            "napc-pe-5-l-1274",
            "napc-pe-5-l-1274-polyol-ester-gear-oil",
            &[
                (
                    "flash_point_temperature",
                    "PE-5-L-1274",
                    516.0,
                    temperature_dims,
                    None,
                ),
                (
                    "reported_specific_gravity",
                    "PE-5-L-1274",
                    0.998,
                    dimensionless,
                    Some(298.0),
                ),
                (
                    "total_acid_number_as_koh_mass_per_oil_mass",
                    "PE-5-L-1274",
                    0.07e-3,
                    dimensionless,
                    None,
                ),
            ],
        ),
        (
            NAPC_PE_5_L_1307_1553_GEAR_OIL_SEED_MANIFEST,
            "napc-pe-5-l-1307-1553",
            "napc-pe-5-l-1307-1553-mil-l-23699-gear-oil",
            &[
                (
                    "flash_point_temperature",
                    "PE-5-L-1307-NASA",
                    539.0,
                    temperature_dims,
                    None,
                ),
                (
                    "pour_point_temperature",
                    "PE-5-L-1307-NASA",
                    220.0,
                    temperature_dims,
                    None,
                ),
                (
                    "flash_point_temperature",
                    "PE-5-L-1553-NASA",
                    539.0,
                    temperature_dims,
                    None,
                ),
                (
                    "pour_point_temperature",
                    "PE-5-L-1553-NASA",
                    213.0,
                    temperature_dims,
                    None,
                ),
                (
                    "reported_specific_gravity",
                    "PE-5-L-1307-and-PE-5-L-1553",
                    1.0,
                    dimensionless,
                    Some(289.0),
                ),
                (
                    "total_acid_number_as_koh_mass_per_oil_mass",
                    "PE-5-L-1307-and-PE-5-L-1553",
                    0.03e-3,
                    dimensionless,
                    None,
                ),
            ],
        ),
    ];
    let directory = fixture_dir();

    for (manifest_relative, stem, expected_pack_id, expected_rows) in seeds {
        let manifest = workspace_path(manifest_relative);
        assert!(
            manifest.is_file(),
            "committed NASA/NAPC gear-oil seed manifest is missing: {manifest_relative}"
        );
        let first_path = directory.join(format!("{stem}-first.fsmatpk"));
        let second_path = directory.join(format!("{stem}-second.fsmatpk"));

        let first = run_compiler(&manifest, &first_path);
        let second = run_compiler(&manifest, &second_path);
        assert!(
            first.status.success(),
            "first {stem} seed compilation failed: {}",
            String::from_utf8_lossy(&first.stderr)
        );
        assert!(
            second.status.success(),
            "second {stem} seed compilation failed: {}",
            String::from_utf8_lossy(&second.stderr)
        );
        assert_eq!(first.stdout, second.stdout, "{stem} decision stream moved");
        assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

        let first_bytes = fs::read(first_path).expect("read first NASA/NAPC gear-oil pack");
        let second_bytes = fs::read(second_path).expect("read second NASA/NAPC gear-oil pack");
        assert_eq!(first_bytes, second_bytes, "{stem} pack bytes moved");
        let decoded =
            NormalizedPack::from_bytes(&first_bytes).expect("decode NASA/NAPC gear-oil pack");
        let pack_hash = decoded.content_hash();
        let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
            .expect("verify NASA/NAPC gear-oil pack identity");

        assert_eq!(decoded.pack_id(), expected_pack_id);
        assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
        assert!(
            decoded
                .redistribution_terms()
                .contains("public use is permitted")
        );
        assert_eq!(decoded.claims().claim_count(), expected_rows.len());
        assert!(decoded.joint_statistics().is_empty());

        assert!(
            decoded
                .claims()
                .claims_for("kinematic_viscosity")
                .is_empty(),
            "Table IV omits the viscosity unit, so no viscosity claim is admissible"
        );
        for &(property, observation_token, expected_value, expected_dims, validity_temperature) in
            expected_rows
        {
            let claims = decoded.claims().claims_for(property);
            let mut matches = claims.iter().copied().filter(|(_, claim)| {
                if claim.observations.len() != 1 {
                    return false;
                }
                decoded
                    .claims()
                    .observation(claim.observations[0])
                    .is_some_and(|observation| observation.specimen.contains(observation_token))
            });
            let (id, claim) = matches
                .next()
                .unwrap_or_else(|| panic!("missing {property} for {observation_token}"));
            assert!(
                matches.next().is_none(),
                "duplicate {property} for {observation_token}"
            );
            let PropertyValue::Scalar { value, dims } = &claim.value else {
                panic!("{property} for {observation_token} was not scalar");
            };
            assert_eq!(*dims, expected_dims);
            let relative_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
            assert!(
                relative_error <= 2.0e-15,
                "{property} for {observation_token} moved by {relative_error:e} relative"
            );
            match validity_temperature {
                Some(temperature_k) => {
                    assert_eq!(
                        claim.validity.bound("temperature"),
                        Some((temperature_k, temperature_k))
                    );
                    assert_eq!(claim.validity.bounds().len(), 1);
                }
                None => assert!(claim.validity.bounds().is_empty()),
            }
            assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
            assert_eq!(
                claim.interpolation,
                InterpolationPolicy::ConstantWithinValidity
            );
            assert_eq!(claim.provenance.license, NASA_SEED_LICENSE);
            assert!(claim.provenance.source.contains("NASA-TM-104352"));
            assert!(claim.provenance.source.contains("[source:primary]"));

            let observation = decoded
                .claims()
                .observation(claim.observations[0])
                .expect("NASA/NAPC gear-oil observation remains linked");
            assert!(observation.method.contains("Table IV"));
            assert!(observation.caveats.contains("proprietary"));
            if observation_token == "PE-5-L-1274" {
                assert!(observation.caveats.contains("NASA reference lubricant"));
                assert!(observation.caveats.contains("without stating their unit"));
                assert!(observation.caveats.contains("less than 200 K"));
            } else if observation_token.contains("-NASA") {
                assert!(
                    observation
                        .caveats
                        .contains("two batches of the same lubricant")
                );
                assert!(observation.caveats.contains("MIL-L-23699"));
                assert!(
                    observation
                        .caveats
                        .contains("batch-specific values remain separate")
                );
                assert!(observation.caveats.contains("without stating their unit"));
            } else {
                assert!(observation.caveats.contains("MIL-L-23699"));
                assert!(observation.caveats.contains("same MIL-L-23699 lubricant"));
            }
            assert_eq!(claim.observations[0].0, observation.content_hash());
            assert_eq!(id.0, claim.content_hash());
        }

        if expected_rows.len() == 6 {
            for property in ["flash_point_temperature", "pour_point_temperature"] {
                let batch_claims = decoded.claims().claims_for(property);
                assert_eq!(batch_claims.len(), 2);
                assert_ne!(
                    batch_claims[0].1.observations, batch_claims[1].1.observations,
                    "the two NASA/NAPC batches were fused for {property}"
                );
            }
        } else {
            assert!(
                decoded
                    .claims()
                    .claims_for("pour_point_temperature")
                    .is_empty(),
                "the censored PE-5-L-1274 pour point became an exact scalar"
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
}

#[test]
fn g3_cli_compiles_committed_rheolube_2000_bearing_grease_seed() {
    let manifest = workspace_path(RHEOLUBE_2000_PENNZANE_GREASE_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Rheolube 2000 bearing-grease seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("rheolube-2000-first.fsmatpk");
    let second_path = directory.join("rheolube-2000-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Rheolube 2000 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Rheolube 2000 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Rheolube 2000 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Rheolube 2000 pack");
    let second_bytes = fs::read(second_path).expect("read second Rheolube 2000 pack");
    assert_eq!(first_bytes, second_bytes, "Rheolube 2000 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode Rheolube 2000 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Rheolube 2000 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "rheolube-2000-pennzane-shf-x-2000-bearing-grease"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 3);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let expected = [
        ("nlgi_consistency_grade", 2.0, dimensionless, None, None),
        ("density", 890.0, density_dims, Some(298.15), None),
        (
            "oil_separation_mass_fraction",
            0.033,
            dimensionless,
            Some(373.15),
            Some(86_400.0),
        ),
    ];
    let mut shared_observation = None;

    for (property, expected_value, expected_dims, temperature_k, duration_s) in expected {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(
            claims.len(),
            1,
            "expected one Rheolube 2000 {property} claim"
        );
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Rheolube 2000 {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let relative_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
        assert!(
            relative_error <= 2.0e-15,
            "Rheolube 2000 {property} moved by {relative_error:e} relative"
        );
        match temperature_k {
            Some(temperature_k) => assert_eq!(
                claim.validity.bound("temperature"),
                Some((temperature_k, temperature_k))
            ),
            None => assert_eq!(claim.validity.bound("temperature"), None),
        }
        match duration_s {
            Some(duration_s) => assert_eq!(
                claim.validity.bound("duration"),
                Some((duration_s, duration_s))
            ),
            None => assert_eq!(claim.validity.bound("duration"), None),
        }
        assert_eq!(
            claim.validity.bounds().len(),
            temperature_k.is_some() as usize + duration_s.is_some() as usize
        );
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-CP-3350"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        assert_eq!(claim.observations.len(), 1);
        match shared_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => shared_observation = Some(claim.observations[0]),
        }

        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("Rheolube 2000 observation remains linked");
        assert_eq!(
            observation.specimen,
            "Rheolube-2000-Pennzane-SHF-X-2000-sodium-octadecylterephthalamate-bearing-grease"
        );
        assert!(observation.method.contains("Bessette Table 7"));
        assert!(observation.method.contains("typical grease properties"));
        assert!(observation.caveats.contains("approximately 20 percent"));
        assert!(observation.caveats.contains("approximately 260 degC"));
        assert!(observation.caveats.contains("labels results typical"));
        assert!(observation.caveats.contains("vacuum-hardening state"));
        assert!(observation.caveats.contains("no printed unit or method"));
        assert!(observation.caveats.contains("not admitted as bulk claims"));
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }

    for refused_property in [
        "dropping_point_temperature",
        "penetration_scale_reading",
        "wear_scar_diameter",
        "oxidation_pressure_drop",
        "vapor_pressure",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "Rheolube 2000 {refused_property} crossed the no-claim boundary"
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
fn g3_cli_compiles_committed_pennzane_shf_x_2000_bearing_oil_seed() {
    let manifest = workspace_path(PENNZANE_SHF_X_2000_BEARING_OIL_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed Pennzane SHF X-2000 bearing-oil seed manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("pennzane-shf-x-2000-first.fsmatpk");
    let second_path = directory.join("pennzane-shf-x-2000-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first Pennzane SHF X-2000 seed compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second Pennzane SHF X-2000 seed compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "Pennzane SHF X-2000 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first Pennzane SHF X-2000 pack");
    let second_bytes = fs::read(second_path).expect("read second Pennzane SHF X-2000 pack");
    assert_eq!(
        first_bytes, second_bytes,
        "Pennzane SHF X-2000 pack bytes moved"
    );
    let decoded =
        NormalizedPack::from_bytes(&first_bytes).expect("decode Pennzane SHF X-2000 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify Pennzane SHF X-2000 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "pennzane-shf-x-2000-mac-aerospace-bearing-oil"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert_eq!(decoded.claims().claim_count(), 7);
    assert!(decoded.joint_statistics().is_empty());

    let kinematic_viscosity_dims = Dims([2, 0, -1, 0, 0, 0]);
    let viscosity_claims = decoded.claims().claims_for("kinematic_viscosity");
    assert_eq!(viscosity_claims.len(), 3);
    let mut shared_observation = None;
    for (source_temperature_c, source_mm2_per_s) in
        [(100.0, 14.3), (40.0, 107.0), (-40.0, 80_500.0)]
    {
        let temperature_k = source_temperature_c + 273.15;
        let mut matches = viscosity_claims.iter().copied().filter(|(_, claim)| {
            claim.validity.bound("temperature") == Some((temperature_k, temperature_k))
        });
        let (id, claim) = matches
            .next()
            .unwrap_or_else(|| panic!("missing Pennzane viscosity at {source_temperature_c} degC"));
        assert!(matches.next().is_none());
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Pennzane viscosity at {source_temperature_c} degC was not scalar");
        };
        assert_eq!(*dims, kinematic_viscosity_dims);
        let expected_m2_per_s = source_mm2_per_s * 1.0e-6;
        let relative_error = (*value - expected_m2_per_s).abs() / expected_m2_per_s;
        assert!(relative_error <= 2.0e-15);
        assert_eq!(claim.validity.bounds().len(), 1);
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.interpolation,
            InterpolationPolicy::ConstantWithinValidity
        );
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        assert!(claim.provenance.source.contains("NASA-CP-3350"));
        assert_eq!(claim.observations.len(), 1);
        match shared_observation {
            Some(observation) => assert_eq!(claim.observations[0], observation),
            None => shared_observation = Some(claim.observations[0]),
        }
        assert_eq!(id.0, claim.content_hash());
    }

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let temperature_dims = Dims([0, 0, 0, 1, 0, 0]);
    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let expected = [
        ("viscosity_index_scale_reading", 137.0, dimensionless, None),
        (
            "flash_point_temperature",
            300.0 + 273.15,
            temperature_dims,
            None,
        ),
        (
            "pour_point_temperature",
            -55.0 + 273.15,
            temperature_dims,
            None,
        ),
        ("density", 0.84 * 1_000.0, density_dims, Some(298.15)),
    ];

    for (property, expected_value, expected_dims, validity_temperature) in expected {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one Pennzane {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("Pennzane {property} was not scalar");
        };
        assert_eq!(*dims, expected_dims);
        let relative_error = (*value - expected_value).abs() / expected_value.abs().max(1.0);
        assert!(relative_error <= 2.0e-15);
        match validity_temperature {
            Some(temperature_k) => {
                assert_eq!(
                    claim.validity.bound("temperature"),
                    Some((temperature_k, temperature_k))
                );
                assert_eq!(claim.validity.bounds().len(), 1);
            }
            None => assert!(claim.validity.bounds().is_empty()),
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(claim.provenance.license, PUBLIC_USE_PERMITTED_LICENSE);
        assert_eq!(claim.observations.len(), 1);
        assert_eq!(
            claim.observations[0],
            shared_observation.expect("shared observation")
        );
        assert_eq!(id.0, claim.content_hash());
    }

    let observation_id = shared_observation.expect("Pennzane observation id");
    let observation = decoded
        .claims()
        .observation(observation_id)
        .expect("Pennzane observation remains linked");
    assert_eq!(
        observation.specimen,
        "Pennzane-SHF-X-2000-multiply-alkylated-cyclopentane-aerospace-bearing-oil"
    );
    assert!(observation.method.contains("Bessette Table 6"));
    assert!(observation.method.contains("typical Pennzane properties"));
    assert!(
        observation
            .caveats
            .contains("Tris(2-octyldodecyl) cyclopentane")
    );
    assert!(
        observation
            .caveats
            .contains("approximate molecular weight 910 g/mol")
    );
    assert!(observation.caveats.contains("labels results typical"));
    assert!(
        observation
            .caveats
            .contains("no temperature-interval degree-Celsius token")
    );
    assert!(observation.caveats.contains("tribometer-system result"));
    assert!(observation.caveats.contains("not unambiguous"));
    assert_eq!(observation_id.0, observation.content_hash());

    for refused_property in [
        "volumetric_thermal_expansion_coefficient",
        "wear_scar_diameter",
        "vapor_pressure",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "Pennzane {refused_property} crossed the no-claim boundary"
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
fn g3_cli_compiles_committed_nasa_cr_195445_omc_ps200_rotary_coating_system() {
    let manifest = workspace_path(NASA_CR_195445_OMC_PS200_ROTARY_COATING_SEED_MANIFEST);
    assert!(
        manifest.is_file(),
        "committed NASA-CR-195445 OMC PS-200 coating manifest is missing"
    );
    let directory = fixture_dir();
    let first_path = directory.join("nasa-cr-195445-omc-ps200-first.fsmatpk");
    let second_path = directory.join("nasa-cr-195445-omc-ps200-second.fsmatpk");

    let first = run_compiler(&manifest, &first_path);
    let second = run_compiler(&manifest, &second_path);
    assert!(
        first.status.success(),
        "first NASA-CR-195445 OMC PS-200 compilation failed: {}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(
        second.status.success(),
        "second NASA-CR-195445 OMC PS-200 compilation failed: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    assert_eq!(
        first.stdout, second.stdout,
        "NASA-CR-195445 OMC PS-200 decision stream moved"
    );
    assert_decision_compiler(&first, MATERIAL_COMPILER_ID);

    let first_bytes = fs::read(first_path).expect("read first OMC PS-200 pack");
    let second_bytes = fs::read(second_path).expect("read second OMC PS-200 pack");
    assert_eq!(first_bytes, second_bytes, "OMC PS-200 pack bytes moved");
    let decoded = NormalizedPack::from_bytes(&first_bytes).expect("decode OMC PS-200 pack");
    let pack_hash = decoded.content_hash();
    let decoded = NormalizedPack::from_bytes_verified(pack_hash, &first_bytes)
        .expect("verify OMC PS-200 pack identity");

    assert_eq!(
        decoded.pack_id(),
        "nasa-cr-195445-omc-ps200-rotary-coating-system"
    );
    assert_eq!(decoded.compiler(), MATERIAL_COMPILER_ID);
    assert!(
        decoded
            .redistribution_terms()
            .contains("public use permitted")
    );
    assert!(decoded.redistribution_terms().contains("US4728448A"));
    assert_eq!(decoded.claims().claim_count(), 7);
    assert!(decoded.joint_statistics().is_empty());

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let mut composition_sum = 0.0_f64;
    let mut composition_observation = None;
    for (property, source_percent) in [
        (
            "ps200_bonded_chromium_carbide_feedstock_mass_fraction",
            80.0_f64,
        ),
        ("ps200_silver_feedstock_mass_fraction", 10.0_f64),
        ("ps200_baf2_caf2_eutectic_feedstock_mass_fraction", 10.0_f64),
    ] {
        let claims = decoded.claims().claims_for(property);
        assert_eq!(claims.len(), 1, "expected one PS-200 {property} claim");
        let (id, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PS-200 {property} was not scalar");
        };
        assert_eq!(*dims, dimensionless);
        let expected_value = source_percent * 0.01;
        assert_eq!(value.to_bits(), expected_value.to_bits());
        composition_sum += *value;
        for required_axis in [
            "source_feedstock_composition_basis_is_mass_fraction",
            "source_ps200_plasma_sprayed_in_engine_report",
        ] {
            assert_eq!(claim.validity.bound(required_axis), Some((1.0, 1.0)));
        }
        for missing_or_refused_axis in [
            "source_post_spray_phase_fractions_known",
            "source_powder_lot_known",
            "source_patent_practice_license_granted",
        ] {
            assert_eq!(
                claim.validity.bound(missing_or_refused_axis),
                Some((0.0, 0.0))
            );
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.provenance.license,
            PUBLIC_USE_AND_PATENT_PUBLICATION_LICENSE
        );
        assert!(claim.provenance.source.contains("NASA-CR-195445"));
        assert!(claim.provenance.source.contains("US Patent 4,728,448"));
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("PS-200 composition observation remains linked");
        assert!(observation.specimen.contains("80wtpct-nickel-bonded-Cr3C2"));
        assert!(observation.method.contains("US4728448A Table II"));
        assert!(observation.caveats.contains("pre-spray feedstock"));
        assert!(observation.caveats.contains("not a patent-practice"));
        match composition_observation {
            Some(expected) => assert_eq!(claim.observations[0], expected),
            None => composition_observation = Some(claim.observations[0]),
        }
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }
    assert!((composition_sum - 1.0).abs() <= f64::EPSILON);

    let roughness_claims = decoded.claims().claims_for("surface_roughness_rms");
    assert_eq!(roughness_claims.len(), 4);
    let mut observed_conditions = Vec::new();
    for (id, claim) in roughness_claims {
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("OMC PS-200 RMS finish was not scalar");
        };
        assert_eq!(*dims, Dims([1, 0, 0, 0, 0, 0]));
        let (test_number, _) = claim
            .validity
            .bound("source_engine_test_number")
            .expect("engine test number remains pinned");
        let (stage_index, _) = claim
            .validity
            .bound("source_surface_stage_index")
            .expect("surface stage remains pinned");
        let condition = (test_number as u8, stage_index as u8);
        let expected_from_microinch: f64 = match condition {
            (3, 0) => 21.0 * 25.4e-9,
            (3, 1) => 7.0 * 25.4e-9,
            (6, 0) => 24.0 * 25.4e-9,
            (6, 1) => 17.0 * 25.4e-9,
            other => panic!("unexpected OMC PS-200 test/stage condition: {other:?}"),
        };
        assert!(
            (*value - expected_from_microinch).abs() <= 1.0e-21,
            "OMC PS-200 microinch-to-metre transcription moved for {condition:?}"
        );
        observed_conditions.push(condition);
        for required_axis in [
            "source_ps200_over_zirconia_or_sx331",
            "source_substrate_is_aluminum_alloy",
        ] {
            assert_eq!(claim.validity.bound(required_axis), Some((1.0, 1.0)));
        }
        for missing_axis in [
            "source_aluminum_alloy_grade_known",
            "source_coating_thickness_known",
            "source_surface_finish_method_known",
        ] {
            assert_eq!(claim.validity.bound(missing_axis), Some((0.0, 0.0)));
        }
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert_eq!(
            claim.provenance.license,
            PUBLIC_USE_AND_PATENT_PUBLICATION_LICENSE
        );
        let observation = decoded
            .claims()
            .observation(claim.observations[0])
            .expect("OMC PS-200 finish observation remains linked");
        match test_number as u8 {
            3 => {
                assert_eq!(
                    claim.validity.bound("source_narrative_run_duration"),
                    Some((2.5 * 3_600.0, 2.5 * 3_600.0))
                );
                assert!(observation.specimen.contains("Test3-air-cooled-OMC"));
                assert!(observation.caveats.contains("TBC crack"));
                assert!(observation.caveats.contains("scrap that housing"));
                assert!(observation.caveats.contains("no zero-wear"));
            }
            6 => {
                assert_eq!(
                    claim.validity.bound("source_actual_run_duration"),
                    Some((1.5 * 3_600.0, 1.5 * 3_600.0))
                );
                assert_eq!(
                    claim
                        .validity
                        .bound("source_local_ps200_breakthrough_present"),
                    Some((1.0, 1.0))
                );
                assert!(observation.specimen.contains("Test6-air-cooled-OMC"));
                assert!(observation.caveats.contains("500 degF"));
                assert!(observation.caveats.contains("local PS-200 breakthrough"));
                assert!(
                    observation
                        .caveats
                        .contains("does not establish a wear rate")
                );
            }
            other => panic!("unexpected OMC PS-200 engine test number: {other}"),
        }
        assert_eq!(claim.observations[0].0, observation.content_hash());
        assert_eq!(id.0, claim.content_hash());
    }
    observed_conditions.sort_unstable();
    assert_eq!(observed_conditions, [(3, 0), (3, 1), (6, 0), (6, 1)]);

    for refused_property in [
        "coefficient_of_friction",
        "wear_rate",
        "coating_thickness",
        "thermal_conductivity",
        "specific_fuel_consumption",
        "service_life",
    ] {
        assert!(
            decoded.claims().claims_for(refused_property).is_empty(),
            "source-absent or system-level coating property must remain refused: {refused_property}"
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
