//! G3 CLI conformance for the committed source-bound interface seed packs.

#![deny(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedInterfacePack, PropertyValue, UncertaintyModel};
use fs_qty::Dims;

const INTERFACE_COMPILER_ID: &str = "frankensim-matdb-interface-pack-compiler-v1";
const DRY_52100_MANIFEST: &str = "data/matdb/seed-v1/nasa-52100-dry-air-interface/manifest.tsv";
const GREASED_52100_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-52100-gxl320a-vacuum-interface/manifest.tsv";
const JOURNAL_4340_BRONZE_MANIFEST: &str =
    "data/matdb/seed-v1/nasa-tn-d-2223-4340-high-lead-bronze-journal/manifest.tsv";
const CARBON_PTFE_CR_ROD_MANIFEST: &str =
    "data/matdb/seed-v1/zhang-2021-carbon-ptfe-cr-piston-rod/manifest.tsv";
const PTFE_CF10_STEEL_BORE_MANIFEST: &str =
    "data/matdb/seed-v1/deaconescu-2020-ptfe-cf10-steel-bore-interface/manifest.tsv";
const MAHLE_PTFE_NICKEL_SIC_COATED_BORE_MANIFEST: &str =
    "data/matdb/seed-v1/mahle-wo2019072721-ptfe-nickel-sic-coated-bore-interface/manifest.tsv";
const A2017_LLC_RA005_MANIFEST: &str =
    "data/matdb/seed-v1/yilmaz-2026-a2017-seiken-llc-ra005-wetting/manifest.tsv";
const A2017_LLC_RA3_MANIFEST: &str =
    "data/matdb/seed-v1/yilmaz-2026-a2017-seiken-llc-ra3-wetting/manifest.tsv";
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

#[test]
fn g3_cli_compiles_committed_carbon_ptfe_chrome_rod_interface() {
    let (pack, decisions) = compile_twice(CARBON_PTFE_CR_ROD_MANIFEST);

    assert_eq!(
        pack.pack_id(),
        "zhang-2021-carbon-fiber-ptfe-cr-piston-rod-kunlun15-interface"
    );
    assert_eq!(
        pack.card().surface_a().material.chemistry,
        "PTFE filled with about 15 percent carbon fiber"
    );
    assert_eq!(
        pack.card().surface_b().material.chemistry,
        "stainless steel with about 50 micrometer electroplated chromium coating"
    );
    assert_ne!(
        pack.card().surface_a().texture_frame,
        pack.card().surface_b().texture_frame,
        "ordered seal and coated-rod roles collapsed"
    );
    assert_eq!(pack.card().medium(), "Kunlun-15-aviation-hydraulic-oil");
    assert_eq!(pack.card().third_body(), None);
    assert_eq!(pack.card().environment(), "laboratory-at-293.15-K");
    assert_eq!(
        pack.card().history(),
        "new-Sterling-seals-300-reciprocations-two-experiment-average"
    );
    assert_eq!(pack.claims_pack().claims().claim_count(), 4);

    for (property, expected_force, pressure, speed) in [
        (
            "single-seal-instroke-friction-force-at-source-observed-minimum",
            97.8955,
            10_000_000.0,
            0.5,
        ),
        (
            "single-seal-outstroke-friction-force-at-source-observed-minimum",
            75.0238,
            10_000_000.0,
            0.5,
        ),
        (
            "single-seal-instroke-friction-force-at-source-observed-maximum",
            516.9906,
            35_000_000.0,
            0.1,
        ),
        (
            "single-seal-outstroke-friction-force-at-source-observed-maximum",
            404.4382,
            35_000_000.0,
            0.1,
        ),
    ] {
        assert_eq!(scalar(&pack, property), expected_force);
        let claim = pack.card().claims_for(property)[0].1;
        for (axis, value) in [
            ("experimental_temperature", 293.15),
            ("hydraulic_oil_pressure", pressure),
            ("reciprocating_speed", speed),
            ("reciprocation_count", 300.0),
            ("sensor_sampling_frequency", 2_400.0),
            ("experiments_averaged", 2.0),
            ("seals_per_experiment", 2.0),
            ("source_fill_fraction_approximately_15_percent", 1.0),
            (
                "source_chromium_thickness_approximately_50_micrometers",
                1.0,
            ),
        ] {
            assert_eq!(claim.validity.bound(axis), Some((value, value)));
        }
        assert_eq!(claim.validity.bounds().len(), 9);
    }

    for refused_property in [
        "coated-bore-friction-force",
        "kinetic-friction-coefficient",
        "wear-rate",
        "seal-life",
        "leakage-rate",
        "transferable-friction-law",
    ] {
        assert!(
            pack.card().claims_for(refused_property).is_empty(),
            "seal source crossed the {refused_property} no-claim boundary"
        );
    }
    assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
}

#[test]
fn g3_cli_compiles_committed_ptfe_cf10_steel_bore_iso_vg32_interface() {
    let (pack, decisions) = compile_twice(PTFE_CF10_STEEL_BORE_MANIFEST);

    assert_eq!(
        pack.pack_id(),
        "deaconescu-2020-ptfe-cf10-steel-bore-iso-vg32-interface"
    );
    assert_eq!(
        pack.card().surface_a().material.chemistry,
        "90 wt% PTFE plus 10 wt% carbon fiber"
    );
    assert_eq!(
        pack.card().surface_b().material.chemistry,
        "hydraulic-cylinder steel inner bore"
    );
    assert_eq!(
        pack.card().surface_b().texture_frame,
        "cylinder-inner-bore-coating-and-finish-unstated"
    );
    assert_ne!(
        pack.card().surface_a().texture_frame,
        pack.card().surface_b().texture_frame,
        "ordered piston-seal and cylinder-bore roles collapsed"
    );
    assert_eq!(
        pack.card().medium(),
        "source-linked-CITGO-Anti-Wear-Hydraulic-Oil-Grade-32-lot-and-formulation-unstated"
    );
    assert_eq!(pack.card().third_body(), None);
    assert_eq!(
        pack.card().environment(),
        "hydraulic-cylinder-opposed-chamber-pressure-and-atmosphere-unstated"
    );
    assert_eq!(
        pack.card().history(),
        "mounted-coaxial-piston-seal-15pct-o-ring-radial-deformation-table-2-resistive-film-measurement"
    );
    assert_eq!(pack.claims_pack().claims().claim_count(), 7);

    let claims = pack.card().claims_for("fluid-film-thickness");
    assert_eq!(claims.len(), 7);
    let mut observed = Vec::new();
    for (_, claim) in claims {
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("PTFE CF10 fluid-film thickness was not scalar");
        };
        assert_eq!(*dims, Dims([1, 0, 0, 0, 0, 0]));
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert!(
            claim
                .provenance
                .license
                .contains("Creative Commons Attribution 4.0")
        );
        assert!(claim.provenance.source.contains("10.3390/polym12010155"));
        assert!(claim.provenance.source.contains("[source:primary]"));

        let pressure = claim
            .validity
            .bound("working_pressure")
            .expect("pressure remains exact")
            .0;
        let temperature = claim
            .validity
            .bound("oil_temperature")
            .expect("temperature remains exact")
            .0;
        let speed = claim
            .validity
            .bound("piston_speed")
            .expect("speed remains exact")
            .0;
        let voltage = claim
            .validity
            .bound("source_measured_voltage_drop")
            .expect("voltage remains exact")
            .0;
        for (axis, expected) in [
            ("source_o_ring_shore_a_scale_reading", (70.0, 70.0)),
            (
                "source_o_ring_initial_radial_deformation_fraction",
                (0.15, 0.15),
            ),
            (
                "source_film_thickness_derived_by_resistive_method",
                (1.0, 1.0),
            ),
            ("source_series_resistance", (10_000_000.0, 10_000_000.0)),
            ("source_cylinder_material_steel", (1.0, 1.0)),
            ("source_cylinder_bore_grade_known", (0.0, 0.0)),
            ("source_cylinder_bore_coating_state_known", (0.0, 0.0)),
            ("source_cylinder_bore_measured_finish_known", (0.0, 0.0)),
            ("source_ptfe_cf10_supplier_lot_process_known", (0.0, 0.0)),
            ("source_oil_lot_and_formulation_known", (0.0, 0.0)),
            ("source_replicates_dispersion_confidence_known", (0.0, 0.0)),
        ] {
            assert_eq!(claim.validity.bound(axis), Some(expected));
        }
        assert_eq!(claim.validity.bounds().len(), 15);

        let observation = pack
            .claims_pack()
            .claims()
            .observation(claim.observations[0])
            .expect("PTFE CF10 observation remains linked");
        let block = if observation.method.contains("velocity block") {
            "velocity"
        } else if observation.method.contains("pressure block") {
            "pressure"
        } else if observation.method.contains("temperature block") {
            "temperature"
        } else {
            panic!("unrecognized Table 2 block: {}", observation.method);
        };
        assert!(observation.specimen.contains("steel cylinder bore"));
        observed.push((block, pressure, temperature, speed, voltage, *value));
    }

    let expected_rows = [
        ("velocity", 10_000_000.0, 333.0, 0.2, 0.186, 8.3e-6),
        ("velocity", 10_000_000.0, 333.0, 0.5, 0.143, 10.8e-6),
        ("pressure", 5_000_000.0, 333.0, 0.2, 0.173, 8.9e-6),
        ("pressure", 10_000_000.0, 333.0, 0.2, 0.186, 8.3e-6),
        ("pressure", 15_000_000.0, 333.0, 0.2, 0.222, 6.9e-6),
        ("pressure", 20_000_000.0, 333.0, 0.2, 0.591, 2.5e-6),
        ("temperature", 10_000_000.0, 333.0, 0.2, 0.19, 8.1e-6),
    ];
    for expected in expected_rows {
        assert_eq!(
            observed.iter().filter(|row| **row == expected).count(),
            1,
            "missing or duplicate Table 2 row {expected:?}"
        );
    }

    let mut same_condition_values = observed
        .iter()
        .filter_map(|(_, pressure, temperature, speed, _, thickness)| {
            (*pressure == 10_000_000.0 && *temperature == 333.0 && *speed == 0.2)
                .then_some(*thickness)
        })
        .collect::<Vec<_>>();
    same_condition_values.sort_by(f64::total_cmp);
    assert_eq!(same_condition_values, [8.1e-6, 8.3e-6, 8.3e-6]);

    for refused_property in [
        "coated-bore-fluid-film-thickness",
        "kinetic-friction-coefficient",
        "leakage-rate",
        "wear-rate",
        "seal-life",
        "transferable-lubrication-law",
    ] {
        assert!(
            pack.card().claims_for(refused_property).is_empty(),
            "steel-bore source crossed the {refused_property} no-claim boundary"
        );
    }
    assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
}

#[test]
fn g3_cli_compiles_committed_mahle_ptfe_nickel_sic_coated_bore_patent_interface() {
    let (pack, decisions) = compile_twice(MAHLE_PTFE_NICKEL_SIC_COATED_BORE_MANIFEST);

    assert_eq!(
        pack.pack_id(),
        "mahle-wo2019072721-ptfe-nickel-sic-coated-bore-interface"
    );
    assert_eq!(
        pack.card().surface_a().material.chemistry,
        "polytetrafluoroethylene PTFE piston ring"
    );
    assert_eq!(
        pack.card().surface_b().material.chemistry,
        "cylinder running surface with nickel-main-component liner layer containing silicon carbide"
    );
    assert_eq!(
        pack.card().surface_b().texture_frame,
        "honed-nickel-sic-liner-patent-claimed-ra-approx-0.1um-not-measured"
    );
    assert_eq!(
        pack.card().medium(),
        "oil-free-compressor-working-gas-composition-pressure-humidity-and-temperature-unstated"
    );
    assert_eq!(pack.card().third_body(), None);
    assert_eq!(
        pack.card().history(),
        "patent-disclosed-ptfe-piston-ring-against-nickel-sic-coated-cylinder-running-surface"
    );
    assert_eq!(pack.claims_pack().claims().claim_count(), 11);

    let dimensionless = Dims([0, 0, 0, 0, 0, 0]);
    let length = Dims([1, 0, 0, 0, 0, 0]);
    for (property, expected_value, expected_dims, scope, semantic_axis) in [
        (
            "patent_claim_liner_silicon_carbide_mass_fraction_approximate_lower_bound",
            2.0 * 0.01,
            dimensionless,
            "claims",
            "source_value_is_approximate",
        ),
        (
            "patent_claim_liner_silicon_carbide_mass_fraction_approximate_upper_bound",
            6.0 * 0.01,
            dimensionless,
            "claims",
            "source_value_is_approximate",
        ),
        (
            "patent_claim_silicon_carbide_particle_diameter_lower_exclusive_bound",
            1.0e-6,
            length,
            "claims",
            "source_bound_is_strict_exclusive",
        ),
        (
            "patent_claim_silicon_carbide_particle_diameter_preferred",
            2.5e-6,
            length,
            "claims",
            "source_value_is_preferred",
        ),
        (
            "patent_claim_silicon_carbide_particle_diameter_upper_exclusive_bound",
            4.0e-6,
            length,
            "claims",
            "source_bound_is_strict_exclusive",
        ),
        (
            "patent_claim_liner_layer_thickness_lower_exclusive_bound",
            40.0e-6,
            length,
            "claims",
            "source_bound_is_strict_exclusive",
        ),
        (
            "patent_claim_liner_layer_thickness_preferred",
            50.0e-6,
            length,
            "claims",
            "source_value_is_preferred",
        ),
        (
            "patent_claim_liner_layer_thickness_upper_exclusive_bound",
            100.0e-6,
            length,
            "claims",
            "source_bound_is_strict_exclusive",
        ),
        (
            "patent_claim_honed_liner_arithmetic_mean_roughness_approximate",
            0.1e-6,
            length,
            "claims",
            "source_value_is_approximate",
        ),
        (
            "patent_embodiment_liner_silicon_carbide_mass_fraction_lower_bound",
            2.0 * 0.01,
            dimensionless,
            "embodiment",
            "source_statement_scope_embodiment_description",
        ),
        (
            "patent_embodiment_liner_silicon_carbide_mass_fraction_upper_bound",
            4.0 * 0.01,
            dimensionless,
            "embodiment",
            "source_statement_scope_embodiment_description",
        ),
    ] {
        let claims = pack.card().claims_for(property);
        assert_eq!(
            claims.len(),
            1,
            "missing or duplicate patent claim {property}"
        );
        let (_, claim) = claims[0];
        let PropertyValue::Scalar { value, dims } = &claim.value else {
            panic!("MAHLE patent claim {property} was not scalar");
        };
        assert_eq!(*value, expected_value);
        assert_eq!(*dims, expected_dims);
        assert_eq!(claim.uncertainty, UncertaintyModel::Unstated);
        assert!(claim.provenance.source.contains("WO 2019/072721 A1"));
        assert!(claim.provenance.source.contains("[source:primary]"));
        assert!(
            claim.provenance.license.contains(
                "German-UrhG-Section-5-2-Published-Patent-Document-No-Copyright-Protection"
            )
        );

        for (axis, expected) in [
            (
                "source_is_patent_design_statement_not_measured_specimen",
                (1.0, 1.0),
            ),
            ("source_lubricant_present", (0.0, 0.0)),
            (
                "source_ptfe_grade_filler_manufacturing_finish_known",
                (0.0, 0.0),
            ),
            (
                "source_cylinder_substrate_liner_alloy_application_process_known",
                (0.0, 0.0),
            ),
            (
                "source_operating_pressure_speed_temperature_gas_known",
                (0.0, 0.0),
            ),
            (
                "source_quantitative_friction_wear_leakage_lifetime_measured",
                (0.0, 0.0),
            ),
            ("source_patent_practice_license_granted", (0.0, 0.0)),
        ] {
            assert_eq!(claim.validity.bound(axis), Some(expected));
        }
        assert_eq!(claim.validity.bound(semantic_axis), Some((1.0, 1.0)));

        let observation = pack
            .claims_pack()
            .claims()
            .observation(claim.observations[0])
            .expect("MAHLE patent observation remains linked");
        match scope {
            "claims" => {
                assert!(observation.method.contains("claims 1 and 4 through 10"));
                assert_eq!(
                    claim.validity.bound("source_statement_scope_patent_claims"),
                    Some((1.0, 1.0))
                );
                assert_eq!(
                    claim
                        .validity
                        .bound("source_statement_scope_embodiment_description"),
                    Some((0.0, 0.0))
                );
            }
            "embodiment" => {
                assert!(observation.method.contains("disclosed embodiment"));
                assert_eq!(
                    claim.validity.bound("source_statement_scope_patent_claims"),
                    Some((0.0, 0.0))
                );
            }
            _ => panic!("unexpected patent-statement scope {scope}"),
        }
    }

    let roughness = pack
        .card()
        .claims_for("patent_claim_honed_liner_arithmetic_mean_roughness_approximate");
    assert_eq!(
        roughness[0].1.validity.bound("source_surface_is_honed"),
        Some((1.0, 1.0))
    );

    for refused_property in [
        "kinetic-friction-coefficient",
        "wear-rate",
        "leakage-rate",
        "seal-life",
        "lubricated-coated-bore-performance",
        "transferable-friction-or-durability-law",
    ] {
        assert!(
            pack.card().claims_for(refused_property).is_empty(),
            "patent design disclosure crossed the {refused_property} no-claim boundary"
        );
    }
    assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
}

#[test]
fn g3_cli_compiles_committed_a2017_seiken_llc_air_wetting_interfaces() {
    for (manifest, pack_id, texture_frame, roughness, angle) in [
        (
            A2017_LLC_RA005_MANIFEST,
            "yilmaz-2026-a2017-seiken-llc-ra005-air-wetting-interface",
            "pre-boiling-a2017-ra005um",
            0.000_000_05,
            1.094_321_441_000_444_7,
        ),
        (
            A2017_LLC_RA3_MANIFEST,
            "yilmaz-2026-a2017-seiken-llc-ra3-air-wetting-interface",
            "pre-boiling-a2017-ra3um-edm",
            0.000_003,
            1.795_769_267_376_965_6,
        ),
    ] {
        let (pack, decisions) = compile_twice(manifest);

        assert_eq!(pack.pack_id(), pack_id);
        assert_eq!(
            pack.card().surface_a().material.chemistry,
            "A2017 aluminum alloy; source wt% Si 0.2..0.8 Fe max 0.7 Cu 3.5..4.5 Mn 0.4..1.0 Mg max 0.1 Cr max 0.25 Zn max 0.15 Ti max 0.15 others max 0.2"
        );
        assert_eq!(pack.card().surface_a().texture_frame, texture_frame);
        assert_eq!(
            pack.card().surface_b().material.chemistry,
            "Seiken Chemical Industry long-life coolant source sample"
        );
        assert_eq!(
            pack.card().surface_b().texture_frame,
            "approximately-1-microliter-free-surface"
        );
        assert_eq!(
            pack.card().medium(),
            "Seiken-Chemical-Industry-long-life-coolant-source-sample"
        );
        assert_eq!(pack.card().third_body(), None);
        assert_eq!(
            pack.card().environment(),
            "ambient-air-297.65-K-60pct-RH-pressure-and-composition-unstated"
        );
        assert_eq!(
            pack.card().history(),
            "pre-boiling-IPA-cleaned-five-static-angle-measurements"
        );
        assert_eq!(pack.claims_pack().claims().claim_count(), 1);
        assert_eq!(scalar(&pack, "static-contact-angle"), angle);

        let claim = pack.card().claims_for("static-contact-angle")[0].1;
        for (axis, value) in [
            ("ambient_temperature", 297.65),
            ("relative_humidity", 0.6),
            ("surface_roughness_ra", roughness),
            ("nominal_droplet_volume", 0.000_000_001),
            ("source_droplet_volume_approximately_1_microliter", 1.0),
            ("measurements_averaged", 5.0),
            ("source_theta_over_2_method", 1.0),
            ("source_measurement_before_boiling", 1.0),
            ("source_surface_cleaned_with_isopropyl_alcohol", 1.0),
            ("source_air_pressure_known", 0.0),
            ("source_air_composition_known", 0.0),
            ("source_llc_product_code_known", 0.0),
            ("source_llc_formulation_known", 0.0),
        ] {
            assert_eq!(claim.validity.bound(axis), Some((value, value)));
        }
        assert_eq!(claim.validity.bounds().len(), 13);

        for refused_property in [
            "advancing-contact-angle",
            "receding-contact-angle",
            "contact-angle-hysteresis",
            "post-boiling-contact-angle",
            "surface-energy",
            "contact-angle-temperature-law",
            "transferable-wetting-law",
        ] {
            assert!(
                pack.card().claims_for(refused_property).is_empty(),
                "wetting source crossed the {refused_property} no-claim boundary"
            );
        }
        assert!(decisions.contains("\"reason_code\":\"interface_context_admitted\""));
    }
}
