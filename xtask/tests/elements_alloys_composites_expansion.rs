//! G0/G3 source -> real compiler -> typed-query checks for the material expansion.
//! Numerical/source transport evidence only; no experimental validation claim.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::path::Path;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const HEAT_CAPACITY: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CONDUCTIVITY: Dims = Dims([1, 1, -3, -1, 0, 0]);

fn compile(slug: &str, expected_claims: usize) -> NormalizedPack {
    static SERIAL: AtomicU64 = AtomicU64::new(0);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-material-conditions-{}-{time}-{}",
        std::process::id(),
        SERIAL.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    let output = dir.join("pack.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(
        pack.claims().claims_ordered().count(),
        expected_claims,
        "{slug}"
    );
    pack
}

fn point(kelvin: f64, selectors: &[(&str, f64)]) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap();
    for &(name, value) in selectors {
        point = point
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn reference(kelvin: f64, temperatures_matched: f64) -> QueryPoint {
    point(
        kelvin,
        &[
            ("engineering-reference-model", 1.0),
            ("source-temperatures-matched", temperatures_matched),
        ],
    )
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let quantity = if dims == CONDUCTIVITY {
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, quantity)
}

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    assert!(
        (answer.evidence.value.value - expected).abs() <= expected.abs().max(1.0) * 1e-11,
        "{} {name}: {} != {expected}",
        pack.pack_id(),
        answer.evidence.value.value
    );
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} {name} must reject {point:?}",
        pack.pack_id()
    );
}

#[test]
fn mineral_and_element_reference_profiles_retain_values_and_model_boundaries() {
    for (slug, count, matched, rho, cp, conductivity) in [
        (
            "barium-fluoride-baf2-crystran-reference",
            3,
            0.0,
            4890.0,
            410.0,
            11.72,
        ),
        (
            "zinc-selenide-crystran-reference",
            5,
            0.0,
            5270.0,
            339.0,
            18.0,
        ),
        ("chromium-espi-reference", 4, 1.0, 7190.0, 459.8, 67.0),
        ("niobium-espi-reference", 4, 0.0, 8570.0, 267.776, 53.7),
    ] {
        let pack = compile(slug, count);
        for t in [293.15, 293.65, 294.15] {
            let query = reference(t, matched);
            sample(&pack, "density", DENSITY, &query, rho);
            sample(&pack, "specific-heat-capacity", HEAT_CAPACITY, &query, cp);
            sample(
                &pack,
                "thermal-conductivity",
                CONDUCTIVITY,
                &query,
                conductivity,
            );
        }
        refuses(&pack, "density", DENSITY, &point(293.15, &[]));
        refuses(&pack, "density", DENSITY, &reference(293.15, 1.0 - matched));
        refuses(
            &pack,
            "thermal-conductivity",
            CONDUCTIVITY,
            &reference(294.16, matched),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            CONDUCTIVITY,
            &reference(293.14, matched),
        );
    }
    let sapphire = compile("sapphire-al2o3-crystran-reference", 3);
    let query = reference(293.15, 0.0);
    sample(&sapphire, "density", DENSITY, &query, 3970.0);
    sample(
        &sapphire,
        "specific-heat-capacity",
        HEAT_CAPACITY,
        &query,
        763.0,
    );
    sample(
        &sapphire,
        "thermal-conductivity-direction-unspecified",
        CONDUCTIVITY,
        &query,
        27.21,
    );
    refuses(&sapphire, "thermal-conductivity", CONDUCTIVITY, &query);
}

#[test]
fn crystal_heat_conduction_requires_the_actual_principal_direction() {
    for (slug, rho, cp, parallel, perpendicular) in [
        (
            "magnesium-fluoride-mgf2-crystran-reference",
            3180.0,
            1003.0,
            21.0,
            33.6,
        ),
        (
            "quartz-crystal-crystran-reference",
            2649.0,
            710.0,
            10.7,
            6.2,
        ),
        (
            "calcite-caco3-crystran-reference",
            2710.0,
            852.0,
            5.526,
            4.646,
        ),
    ] {
        let pack = compile(slug, 4);
        let common = reference(293.65, 0.0);
        sample(&pack, "density", DENSITY, &common, rho);
        sample(&pack, "specific-heat-capacity", HEAT_CAPACITY, &common, cp);
        refuses(&pack, "thermal-conductivity", CONDUCTIVITY, &common);
        for (direction, expected) in [(1.0, parallel), (0.0, perpendicular)] {
            let query = common
                .clone()
                .with_quantity(
                    "heat-flow-parallel-to-c",
                    QuantitySpec::dimensional(Dims::NONE),
                    direction,
                )
                .unwrap();
            sample(
                &pack,
                "thermal-conductivity",
                CONDUCTIVITY,
                &query,
                expected,
            );
        }
        let oblique = common
            .with_quantity(
                "heat-flow-parallel-to-c",
                QuantitySpec::dimensional(Dims::NONE),
                0.5,
            )
            .unwrap();
        refuses(&pack, "thermal-conductivity", CONDUCTIVITY, &oblique);
    }
}

#[test]
fn g10cr_keeps_directional_equation_ranges_and_unstated_uncertainty() {
    let pack = compile("g10cr-nist-cryogenic", 3);
    let source = [("source-pressure-known", 0.0)];
    sample(
        &pack,
        "specific-heat-capacity",
        HEAT_CAPACITY,
        &point(4.0, &source),
        2.01588856638622,
    );
    sample(
        &pack,
        "specific-heat-capacity",
        HEAT_CAPACITY,
        &point(77.0, &source),
        239.352663464469,
    );
    sample(
        &pack,
        "specific-heat-capacity",
        HEAT_CAPACITY,
        &point(300.0, &source),
        998.742717161209,
    );
    for (direction, lower, at_77, at_300) in [
        (1.0, 10.0, 0.279965413170793, 0.607982754731163),
        (2.0, 12.0, 0.385822232550798, 0.863637841736190),
    ] {
        let selectors = [
            ("source-pressure-known", 0.0),
            ("conductivity-direction", direction),
        ];
        sample(
            &pack,
            "thermal-conductivity",
            CONDUCTIVITY,
            &point(77.0, &selectors),
            at_77,
        );
        sample(
            &pack,
            "thermal-conductivity",
            CONDUCTIVITY,
            &point(300.0, &selectors),
            at_300,
        );
        refuses(
            &pack,
            "thermal-conductivity",
            CONDUCTIVITY,
            &point(lower - 0.01, &selectors),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            CONDUCTIVITY,
            &point(300.01, &selectors),
        );
    }
    refuses(
        &pack,
        "thermal-conductivity",
        CONDUCTIVITY,
        &point(77.0, &source),
    );
    refuses(
        &pack,
        "thermal-conductivity",
        CONDUCTIVITY,
        &point(
            77.0,
            &[
                ("source-pressure-known", 0.0),
                ("conductivity-direction", 3.0),
            ],
        ),
    );
    refuses(&pack, "density", DENSITY, &point(77.0, &source));
    refuses(
        &pack,
        "specific-heat-capacity",
        HEAT_CAPACITY,
        &point(3.99, &source),
    );
    // The stored samples come from high-precision evaluation. An independent
    // ordinary-f64 Horner calculation checks the source polynomial at 77 K.
    let x = 77.0_f64.log10();
    let exponent = [
        -2.4083, 7.6006, -8.2982, 7.3301, -4.2386, 1.4294, -0.24396, 0.015236,
    ]
    .iter()
    .rev()
    .fold(0.0, |value, &coefficient| value * x + coefficient);
    assert!((10.0_f64.powf(exponent) / 239.352663464469 - 1.0).abs() < 1e-10);
}
