//! Producer wood-panel thermal sets, declared layer resistance and humidity references.
//! Exercise real pack compilation, typed quantities, evidence receipts and domain refusal.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const LAYER_R: Dims = Dims([0, -1, 3, 1, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const ALPHA: Dims = Dims([2, 0, -1, 0, 0, 0]);

struct Profile {
    slug: &'static str,
    claims: usize,
    density: f64,
    specific_heat: f64,
    conductivity: f64,
    conductivity_key: &'static str,
    volumetric_heat_capacity: f64,
    dried: bool,
}

const PROFILES: [Profile; 5] = [
    Profile {
        slug: "kronobuild-osb-superfinish-eco-dry",
        claims: 4,
        density: 550.0,
        specific_heat: 1221.0,
        conductivity: 0.098,
        conductivity_key: "thermal-conductivity",
        volumetric_heat_capacity: 671_550.0,
        dried: true,
    },
    Profile {
        slug: "kronobuild-dfp-thermal-reference",
        claims: 5,
        density: 550.0,
        specific_heat: 1600.0,
        conductivity: 0.10,
        conductivity_key: "reference-thermal-conductivity",
        volumetric_heat_capacity: 880_000.0,
        dried: false,
    },
    Profile {
        slug: "kronobuild-betonyp-thermal-reference",
        claims: 5,
        density: 1350.0,
        specific_heat: 1600.0,
        conductivity: 0.26,
        conductivity_key: "reference-thermal-conductivity",
        volumetric_heat_capacity: 2_160_000.0,
        dried: false,
    },
    Profile {
        slug: "gutex-thermowall-thermal-reference",
        claims: 12,
        density: 160.0,
        specific_heat: 2100.0,
        conductivity: 0.040,
        conductivity_key: "declared-thermal-conductivity",
        volumetric_heat_capacity: 336_000.0,
        dried: false,
    },
    Profile {
        slug: "steico-therm-dry-v2-2026",
        claims: 10,
        density: 110.0,
        specific_heat: 2100.0,
        conductivity: 0.037,
        conductivity_key: "declared-thermal-conductivity",
        volumetric_heat_capacity: 231_000.0,
        dried: false,
    },
];

fn compile(profile: &Profile) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-wood-panel-thermal-{}-{nonce}-{}",
        std::process::id(),
        profile.slug
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{}/manifest.tsv", profile.slug)))
        .arg("--out")
        .arg(&path)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        profile.slug,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(path).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), profile.slug);
    assert_eq!(pack.claims().claims_ordered().count(), profile.claims);
    pack
}

fn flags(point: QueryPoint, pairs: &[(&str, f64)]) -> QueryPoint {
    pairs.iter().fold(point, |p, &(name, value)| {
        p.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

fn reference(profile: &Profile) -> QueryPoint {
    let product = format!("source-product-{}", profile.slug);
    let point = flags(
        QueryPoint::new(),
        &[
            (&product, 1.0),
            ("source-reference-state", 1.0),
            ("source-specimens-matched", 0.0),
            ("source-temperature-known", 0.0),
            (
                "source-moisture-condition-known",
                if profile.dried { 1.0 } else { 0.0 },
            ),
        ],
    );
    if profile.dried {
        flags(point, &[("source-dried-state", 1.0)])
    } else {
        point
    }
}

fn transport(profile: &Profile) -> QueryPoint {
    let point = flags(
        reference(profile),
        &[
            ("source-thermal-direction-known", 0.0),
            ("source-test-thickness-known", 0.0),
        ],
    );
    if profile.conductivity_key == "declared-thermal-conductivity" {
        flags(point, &[("source-is-declared-value", 1.0)])
    } else {
        point
    }
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
fn thermal_sets_support_capacity_and_layer_arithmetic_without_invented_conditions() {
    for profile in &PROFILES {
        let pack = compile(profile);
        let point = transport(profile);
        let rho = value(&pack, "density", RHO, &point);
        let cp = value(&pack, "specific-heat-capacity", CP, &point);
        let conductivity = value(&pack, profile.conductivity_key, K, &point);
        close(rho, profile.density);
        close(cp, profile.specific_heat);
        close(conductivity, profile.conductivity);
        close(rho * cp, profile.volumetric_heat_capacity);

        // Caller-chosen 100 mm, 1 m2 reference layer and 1 K temperature difference.
        // This checks SI thermal inputs, not a measured temperature or thickness law.
        close(rho * cp * 0.1, profile.volumetric_heat_capacity / 10.0);
        close(0.1 / conductivity, 0.1 / profile.conductivity);
        for (axis, wrong) in [
            ("source-reference-state", 0.0),
            ("source-specimens-matched", 1.0),
            ("source-temperature-known", 1.0),
            (
                "source-moisture-condition-known",
                if profile.dried { 0.0 } else { 1.0 },
            ),
        ] {
            let unsupported = flags(point.clone(), &[(axis, wrong)]);
            refuses(&pack, "density", RHO, &unsupported);
            refuses(&pack, "specific-heat-capacity", CP, &unsupported);
            refuses(&pack, profile.conductivity_key, K, &unsupported);
        }
        refuses(
            &pack,
            profile.conductivity_key,
            K,
            &flags(point.clone(), &[("source-thermal-direction-known", 1.0)]),
        );
        refuses(
            &pack,
            profile.conductivity_key,
            K,
            &flags(point.clone(), &[("source-test-thickness-known", 1.0)]),
        );
        if profile.dried {
            let alpha = value(&pack, "source-tabulated-thermal-diffusivity", ALPHA, &point);
            close(alpha, 1.46e-7);
            // The producer table rounds alpha; this is a rounding comparison, not a CI.
            assert!((alpha / (conductivity / (rho * cp)) - 1.0).abs() < 0.001);
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &flags(point, &[("source-dried-state", 0.0)]),
            );
        } else {
            // Reference and declared values are not silently relabelled as measurements.
            refuses(&pack, "thermal-conductivity", K, &point);
        }
    }
}

fn vapour_point(profile: &Profile, humidity: f64, temperature: f64) -> QueryPoint {
    flags(
        reference(profile),
        &[
            ("source-temperature-known", 1.0),
            ("source-moisture-condition-known", 1.0),
            ("test-relative-humidity", humidity),
        ],
    )
    .with_quantity(
        "temperature",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        temperature,
    )
    .unwrap()
}

#[test]
fn mdf_and_cement_bonded_board_vapour_references_preserve_discrete_humidity_states() {
    for (profile, dry_max, wet_min) in [(&PROFILES[1], 10.0, 8.0), (&PROFILES[2], 50.0, 20.0)] {
        let pack = compile(profile);
        for (humidity, expected) in [(0.25, dry_max), (0.72, wet_min)] {
            let point = vapour_point(profile, humidity, 296.15);
            close(
                value(
                    &pack,
                    "water-vapor-diffusion-resistance-factor",
                    Dims::NONE,
                    &point,
                ),
                expected,
            );
            refuses(&pack, "density", RHO, &point);
            refuses(&pack, "specific-heat-capacity", CP, &point);
            refuses(
                &pack,
                "water-vapor-diffusion-resistance-factor",
                Dims::NONE,
                &vapour_point(profile, humidity, 303.15),
            );
            refuses(
                &pack,
                "water-vapor-diffusion-resistance-factor",
                Dims::NONE,
                &flags(point, &[("source-moisture-condition-known", 0.0)]),
            );
        }
        refuses(
            &pack,
            "water-vapor-diffusion-resistance-factor",
            Dims::NONE,
            &vapour_point(profile, 0.50, 296.15),
        );
    }
}

#[test]
fn wood_fibre_declarations_preserve_thickness_and_mechanical_roles() {
    for (profile, layers, unsupported_thickness) in [
        (
            &PROFILES[3],
            &[
                (0.08, 2.0),
                (0.10, 2.5),
                (0.12, 3.0),
                (0.14, 3.5),
                (0.16, 4.0),
            ][..],
            0.09,
        ),
        (
            &PROFILES[4],
            &[(0.04, 1.05), (0.06, 1.60), (0.08, 2.15), (0.10, 2.70)][..],
            0.05,
        ),
    ] {
        let pack = compile(profile);
        let layer = |thickness| {
            flags(reference(profile), &[("source-is-declared-value", 1.0)])
                .with_quantity(
                    "nominal-board-thickness",
                    QuantitySpec::dimensional(LENGTH),
                    thickness,
                )
                .unwrap()
        };
        for &(thickness, expected) in layers {
            let point = layer(thickness);
            close(
                value(&pack, "declared-layer-thermal-resistance", LAYER_R, &point),
                expected,
            );
            // R_D is rounded downward in the STEICO source. No uncertainty is inferred.
            let exact_reference_r = thickness / profile.conductivity;
            assert!(expected <= exact_reference_r + 1e-12);
            assert!(exact_reference_r - expected < 0.05);
            refuses(
                &pack,
                "declared-layer-thermal-resistance",
                LAYER_R,
                &flags(point, &[("source-is-declared-value", 0.0)]),
            );
        }
        refuses(
            &pack,
            "declared-layer-thermal-resistance",
            LAYER_R,
            &layer(unsupported_thickness),
        );
        refuses(
            &pack,
            "declared-thermal-conductivity",
            K,
            &flags(transport(profile), &[("source-is-declared-value", 0.0)]),
        );
        let load = flags(
            reference(profile),
            &[
                ("loading-normal-to-board", 1.0),
                ("source-test-rate-known", 0.0),
            ],
        );
        if profile.slug.starts_with("gutex-") {
            let design = flags(
                reference(profile),
                &[
                    ("source-thermal-direction-known", 0.0),
                    ("source-test-thickness-known", 0.0),
                    ("source-is-design-value", 1.0),
                ],
            );
            close(
                value(&pack, "design-thermal-conductivity", K, &design),
                0.042,
            );
            refuses(&pack, "design-thermal-conductivity", K, &transport(profile));
            let declared_load = flags(load.clone(), &[("source-is-declared-value", 1.0)]);
            close(
                value(
                    &pack,
                    "compressive-stress-or-strength-lower-limit",
                    PA,
                    &declared_load,
                ),
                100_000.0,
            );
            close(
                value(
                    &pack,
                    "tensile-strength-perpendicular-lower-limit",
                    PA,
                    &declared_load,
                ),
                10_000.0,
            );
            refuses(
                &pack,
                "compressive-stress-or-strength-lower-limit",
                PA,
                &load,
            );
            // CS(10/Y) does not establish the stress at an exact strain before/after yield.
            refuses(&pack, "compression-stress-10-percent", PA, &declared_load);
            refuses(
                &pack,
                "compressive-stress-or-strength-lower-limit",
                PA,
                &flags(declared_load, &[("loading-normal-to-board", 0.0)]),
            );
        } else {
            let compression = flags(load.clone(), &[("engineering-compression-strain", 0.10)]);
            close(
                value(&pack, "compression-stress-10-percent", PA, &compression),
                50_000.0,
            );
            close(
                value(&pack, "tensile-strength-perpendicular-reference", PA, &load),
                10_000.0,
            );
            refuses(
                &pack,
                "compression-stress-10-percent",
                PA,
                &flags(
                    compression.clone(),
                    &[("engineering-compression-strain", 0.02)],
                ),
            );
            refuses(
                &pack,
                "compression-stress-10-percent",
                PA,
                &flags(compression, &[("source-test-rate-known", 1.0)]),
            );
        }
        refuses(&pack, "youngs-modulus", PA, &load);
    }
}
