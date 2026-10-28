//! Exact 28-day water-cured concrete/mortar reference profiles through the real compiler.
//! The source's derived Cp and ultrasonic modulus retain their measurement roles.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const ALPHA: Dims = Dims([2, 0, -1, 0, 0, 0]);
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

struct Profile {
    slug: &'static str,
    dry_density: f64,
    wet_density: f64,
    conductivity: f64,
    specific_heat: f64,
    measured_diffusivity: f64,
    compressive_strength: f64,
    flexural_strength: f64,
    dynamic_modulus: f64,
    wave_speed: f64,
    porosity: f64,
}

const PROFILES: [Profile; 4] = [
    Profile {
        slug: "concrete-nac-wc045-moreno2019",
        dry_density: 2240.0,
        wet_density: 2390.0,
        conductivity: 1.88,
        specific_heat: 816.4,
        measured_diffusivity: 1.21e-6,
        compressive_strength: 54.9e6,
        flexural_strength: 9.9e6,
        dynamic_modulus: 37.9e9,
        wave_speed: 4290.0,
        porosity: 0.146,
    },
    Profile {
        slug: "concrete-argex-wc045-moreno2019",
        dry_density: 1660.0,
        wet_density: 1870.0,
        conductivity: 1.02,
        specific_heat: 922.1,
        measured_diffusivity: 0.70e-6,
        compressive_strength: 17.7e6,
        flexural_strength: 7.2e6,
        dynamic_modulus: 28.7e9,
        wave_speed: 4220.0,
        porosity: 0.286,
    },
    Profile {
        slug: "concrete-ca25-wc045-moreno2019",
        dry_density: 1990.0,
        wet_density: 2150.0,
        conductivity: 1.39,
        specific_heat: 743.5,
        measured_diffusivity: 0.98e-6,
        compressive_strength: 54.7e6,
        flexural_strength: 10.0e6,
        dynamic_modulus: 30.9e9,
        wave_speed: 4080.0,
        porosity: 0.170,
    },
    Profile {
        slug: "cement-mortar-wc045-moreno2019",
        dry_density: 2090.0,
        wet_density: 2290.0,
        conductivity: 1.67,
        specific_heat: 702.9,
        measured_diffusivity: 1.17e-6,
        compressive_strength: 54.2e6,
        flexural_strength: 9.3e6,
        dynamic_modulus: 33.2e9,
        wave_speed: 4100.0,
        porosity: 0.191,
    },
];

fn compile(slug: &str) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-cementitious-thermal-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&path)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(path).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), 10);
    pack
}

fn flags(point: QueryPoint, pairs: &[(&str, f64)]) -> QueryPoint {
    pairs.iter().fold(point, |p, &(name, value)| {
        p.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

fn reference(moisture_protocol: f64) -> QueryPoint {
    flags(
        QueryPoint::new(),
        &[
            ("source-reference-state", 1.0),
            ("source-temperature-known", 0.0),
            ("water-cement-mass-ratio", 0.45),
            ("source-moisture-protocol", moisture_protocol),
        ],
    )
    .with_quantity(
        "water-curing-duration",
        QuantitySpec::dimensional(TIME),
        28.0 * 86400.0,
    )
    .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let spec = if dims == K {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, spec)
}

fn value(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) -> f64 {
    let result = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    assert_eq!(result.evidence.value.dims, dims);
    assert!(matches!(
        result.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    ));
    assert!(result.receipt.observation_backed);
    assert!(!result.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&result.receipt).unwrap();
    result.evidence.value.value
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} admitted {name} at {point:?}",
        pack.pack_id()
    );
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g3_dry_density_conductivity_and_source_derived_cp_retain_table_values_and_limits() {
    for profile in &PROFILES {
        let pack = compile(profile.slug);
        let dry = reference(0.0);
        close(value(&pack, "density", RHO, &dry), profile.dry_density);
        close(
            value(&pack, "thermal-conductivity", K, &dry),
            profile.conductivity,
        );
        // The caller must acknowledge both source derivation and unreconciled closure.
        refuses(&pack, "specific-heat-capacity", CP, &dry);
        let thermal = flags(
            dry.clone(),
            &[
                ("source-derived-heat-capacity", 1.0),
                ("source-thermal-closure-verified", 0.0),
            ],
        );
        close(
            value(&pack, "specific-heat-capacity", CP, &thermal),
            profile.specific_heat,
        );
        close(
            value(&pack, "measured-thermal-diffusivity", ALPHA, &thermal),
            profile.measured_diffusivity,
        );
        let closure = profile.conductivity / (profile.dry_density * profile.specific_heat);
        assert!((closure - profile.measured_diffusivity).abs() > 1e-8);
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &flags(thermal.clone(), &[("source-thermal-closure-verified", 1.0)]),
        );
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &flags(thermal.clone(), &[("source-derived-heat-capacity", 0.0)]),
        );
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &flags(thermal, &[("source-moisture-protocol", 1.0)]),
        );
        refuses(&pack, "thermal-conductivity", K, &reference(1.0));
        refuses(&pack, "thermal-diffusivity", ALPHA, &dry);
        for axis in ["source-reference-state", "source-temperature-known"] {
            let unsupported = if axis == "source-reference-state" {
                0.0
            } else {
                1.0
            };
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &flags(dry.clone(), &[(axis, unsupported)]),
            );
        }
    }
}

#[test]
fn g0_g3_water_cured_mechanical_observations_are_distinct_from_dry_thermal_state() {
    for profile in &PROFILES {
        let pack = compile(profile.slug);
        let wet = reference(1.0);
        close(value(&pack, "density", RHO, &wet), profile.wet_density);
        assert!(profile.wet_density > profile.dry_density);
        for (name, dims, expected) in [
            ("compressive-strength", PA, profile.compressive_strength),
            ("flexural-strength", PA, profile.flexural_strength),
            ("mean-p-wave-speed", SPEED, profile.wave_speed),
        ] {
            close(value(&pack, name, dims, &wet), expected);
            refuses(&pack, name, dims, &reference(0.0));
        }
        refuses(&pack, "dynamic-youngs-modulus", PA, &wet);
        let dynamic = flags(wet.clone(), &[("source-ultrasonic-dynamic-modulus", 1.0)]);
        close(
            value(&pack, "dynamic-youngs-modulus", PA, &dynamic),
            profile.dynamic_modulus,
        );
        // Wave-derived dynamic stiffness is not an independently measured static law.
        for name in [
            "young-modulus",
            "youngs-modulus",
            "static-compressive-youngs-modulus",
        ] {
            refuses(&pack, name, PA, &dynamic);
        }
        refuses(&pack, "poisson-ratio", Dims::NONE, &dynamic);
        refuses(
            &pack,
            "compressive-strength",
            PA,
            &flags(wet.clone(), &[("water-cement-mass-ratio", 0.55)]),
        );
        let immature = wet
            .with_quantity(
                "water-curing-duration",
                QuantitySpec::dimensional(TIME),
                7.0 * 86400.0,
            )
            .unwrap();
        refuses(&pack, "compressive-strength", PA, &immature);
    }
}

#[test]
fn g0_g3_connected_porosity_keeps_its_vacuum_saturation_protocol() {
    for profile in &PROFILES {
        let pack = compile(profile.slug);
        close(
            value(&pack, "connected-porosity", Dims::NONE, &reference(2.0)),
            profile.porosity,
        );
        refuses(&pack, "connected-porosity", Dims::NONE, &reference(0.0));
        refuses(&pack, "connected-porosity", Dims::NONE, &reference(1.0));
        refuses(&pack, "density", RHO, &reference(2.0));
        refuses(&pack, "density", RHO, &reference(0.5));
        refuses(&pack, "total-porosity", Dims::NONE, &reference(2.0));
    }
}
