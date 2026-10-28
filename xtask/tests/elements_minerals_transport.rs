//! G0/G3 source compiler and typed-query checks for elemental transport data.
//! Source transport and bounded interpolation evidence, not physical qualification.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::path::Path;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const RESISTIVITY: Dims = Dims([3, 1, -3, 0, -2, 0]);
const INVERSE_TEMP: Dims = Dims([0, 0, 0, -1, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    static SERIAL: AtomicU64 = AtomicU64::new(0);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "frankensim-elements-{}-{nanos}-{}.fsmatpk",
        std::process::id(),
        SERIAL.fetch_add(1, Ordering::Relaxed)
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
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
    assert_eq!(pack.claims().claims_ordered().count(), count, "{slug}");
    pack
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values
        .iter()
        .fold(QueryPoint::new(), |point, (name, value)| {
            point
                .with_quantity(*name, QuantitySpec::dimensional(Dims::NONE), *value)
                .unwrap()
        })
}

fn temperature(point: QueryPoint, name: &str, kelvin: f64) -> QueryPoint {
    point
        .with_quantity(
            name,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let quantity = if name == "thermal-conductivity" {
        assert_eq!(dims, K);
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
        .unwrap_or_else(|error| panic!("{} {name}: {error:?}", pack.pack_id()));
    let actual = answer.evidence.value.value;
    assert!(
        (actual - expected).abs() <= (expected.abs() * 1e-11).max(1e-30),
        "{} {name}: {actual} != {expected}",
        pack.pack_id()
    );
    assert_eq!(answer.evidence.value.dims, dims);
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
        "{} {name} unexpectedly admits {point:?}",
        pack.pack_id()
    );
}

fn alloy_reference() -> QueryPoint {
    flags(&[
        ("engineering-reference-model", 1.0),
        ("source-specimens-matched", 0.0),
    ])
}

#[test]
fn copper_alloys_preserve_units_source_temperatures_and_interval_expansion() {
    for (slug, pounds_per_inch3, k_imperial, e_ksi, g_ksi, alpha_micro_per_f) in [
        (
            "brass-c36000-cda-reference",
            0.307,
            67.0,
            14000.0,
            5300.0,
            11.4,
        ),
        (
            "brass-c46400-cda-reference",
            0.304,
            67.0,
            15000.0,
            5600.0,
            11.8,
        ),
        (
            "bronze-c65500-cda-reference",
            0.308,
            21.0,
            15000.0,
            5600.0,
            10.0,
        ),
    ] {
        let pack = compile(slug, 6);
        let at = temperature(alloy_reference(), "temperature", 293.15);
        // Independent conversions from source imperial values and SI definitions.
        let psi = 0.45359237 * 9.80665 / 0.0254_f64.powi(2);
        let rho = pounds_per_inch3 * 0.45359237 / 0.0254_f64.powi(3);
        let cp = 0.09 * 1055.05585262 / (0.45359237 * (5.0 / 9.0));
        let k = k_imperial * 1055.05585262 / (3600.0 * 0.3048 * (5.0 / 9.0));
        for (name, dims, value) in [
            ("density", RHO, rho),
            ("specific-heat-capacity", CP, cp),
            ("thermal-conductivity", K, k),
        ] {
            sample(&pack, name, dims, &at, value);
            refuses(
                &pack,
                name,
                dims,
                &temperature(alloy_reference(), "temperature", 294.15),
            );
            refuses(
                &pack,
                name,
                dims,
                &at.clone().with("source-specimens-matched", 1.0).unwrap(),
            );
            refuses(&pack, name, dims, &QueryPoint::new());
        }
        refuses(&pack, "young-modulus", PA, &at);
        let elastic = at
            .clone()
            .with_quantity(
                "source-elastic-temperature-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        sample(&pack, "young-modulus", PA, &elastic, e_ksi * 1000.0 * psi);
        sample(&pack, "shear-modulus", PA, &elastic, g_ksi * 1000.0 * psi);
        refuses(
            &pack,
            "young-modulus",
            PA,
            &elastic
                .clone()
                .with("source-elastic-temperature-known", 1.0)
                .unwrap(),
        );
        refuses(&pack, "poisson-ratio", Dims::NONE, &elastic);
        let endpoints = temperature(
            temperature(alloy_reference(), "source-reference-temperature", 293.15),
            "source-range-end-temperature",
            573.15,
        );
        sample(
            &pack,
            "mean-linear-thermal-expansion-coefficient",
            INVERSE_TEMP,
            &endpoints,
            alpha_micro_per_f * 1e-6 / (5.0 / 9.0),
        );
        refuses(
            &pack,
            "mean-linear-thermal-expansion-coefficient",
            INVERSE_TEMP,
            &temperature(endpoints.clone(), "source-range-end-temperature", 473.15),
        );
        refuses(
            &pack,
            "linear-thermal-expansion-coefficient",
            INVERSE_TEMP,
            &endpoints,
        );
        // A dimension-only request must not erase conductivity's semantic kind.
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::with_quantity(
                        "thermal-conductivity",
                        QuantitySpec::dimensional(K)
                    ),
                    &at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn nickel_grades_keep_annealed_transport_curves_distinct_and_bounded() {
    for (slug, count, k_knots, rho_knots) in [
        (
            "nickel-200-annealed-special-metals-2006",
            7,
            [
                75.5, 70.3, 66.5, 61.6, 56.8, 55.4, 57.6, 59.7, 61.8, 64.0, 66.1, 68.2,
            ],
            [
                0.050, 0.096, 0.130, 0.185, 0.260, 0.330, 0.365, 0.400, 0.430, 0.460, 0.485, 0.510,
            ],
        ),
        (
            "nickel-201-annealed-special-metals-2006",
            5,
            [
                88.3, 79.3, 73.4, 66.3, 59.9, 56.1, 58.2, 60.6, 62.8, 65.1, 67.7, 69.9,
            ],
            [
                0.040, 0.085, 0.125, 0.175, 0.250, 0.330, 0.375, 0.405, 0.435, 0.465, 0.490, 0.515,
            ],
        ),
    ] {
        let pack = compile(slug, count);
        let base = alloy_reference()
            .with_quantity(
                "annealed-material",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
        // Exact declared Kelvin endpoints avoid subtraction-rounding below support.
        let kelvin = [
            173.15, 293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15, 1073.15,
            1173.15, 1273.15,
        ];
        for i in 0..kelvin.len() {
            let at = temperature(base.clone(), "temperature", kelvin[i]);
            sample(&pack, "thermal-conductivity", K, &at, k_knots[i]);
            sample(
                &pack,
                "electrical-resistivity",
                RESISTIVITY,
                &at,
                rho_knots[i] * 1e-6,
            );
            if i + 1 < kelvin.len() {
                let middle = temperature(
                    base.clone(),
                    "temperature",
                    (kelvin[i] + kelvin[i + 1]) / 2.0,
                );
                sample(
                    &pack,
                    "thermal-conductivity",
                    K,
                    &middle,
                    (k_knots[i] + k_knots[i + 1]) / 2.0,
                );
                sample(
                    &pack,
                    "electrical-resistivity",
                    RESISTIVITY,
                    &middle,
                    (rho_knots[i] + rho_knots[i + 1]) * 0.5e-6,
                );
            }
        }
        let at = temperature(base.clone(), "temperature", 293.15);
        for name_dims in [
            ("thermal-conductivity", K),
            ("electrical-resistivity", RESISTIVITY),
        ] {
            refuses(
                &pack,
                name_dims.0,
                name_dims.1,
                &temperature(base.clone(), "temperature", 173.14),
            );
            refuses(
                &pack,
                name_dims.0,
                name_dims.1,
                &at.clone().with("annealed-material", 0.0).unwrap(),
            );
            refuses(
                &pack,
                name_dims.0,
                name_dims.1,
                &temperature(alloy_reference(), "temperature", 293.15),
            );
        }
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(base.clone(), "temperature", 1273.16),
        );
        let hot = temperature(base, "temperature", 1373.15);
        if count == 7 {
            sample(&pack, "electrical-resistivity", RESISTIVITY, &hot, 0.540e-6);
            refuses(
                &pack,
                "electrical-resistivity",
                RESISTIVITY,
                &temperature(hot.clone(), "temperature", 1373.16),
            );
        } else {
            refuses(&pack, "electrical-resistivity", RESISTIVITY, &hot);
        }
        let constants = alloy_reference()
            .with_quantity(
                "source-constant-temperature-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap()
            .with_quantity(
                "source-constant-heat-treatment-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        sample(&pack, "density", RHO, &constants, 8890.0);
        sample(&pack, "specific-heat-capacity", CP, &constants, 456.0);
        refuses(&pack, "specific-heat-capacity", CP, &at);
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &constants
                .clone()
                .with("source-constant-temperature-known", 1.0)
                .unwrap(),
        );
    }
}

#[test]
fn nickel_dynamic_elasticity_is_not_a_static_modulus_or_a_known_frequency() {
    let pack = compile("nickel-200-annealed-special-metals-2006", 7);
    let base = alloy_reference()
        .with_quantity(
            "annealed-material",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap()
        .with_quantity(
            "source-dynamic-method",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap()
        .with_quantity(
            "source-frequency-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    for (celsius, e, g, nu) in [
        (26.0, 205.0, 79.6, 0.29),
        (100.0, 200.0, 77.9, 0.28),
        (200.0, 195.0, 75.8, 0.29),
        (300.0, 190.0, 73.8, 0.29),
        (400.0, 183.0, 71.4, 0.28),
        (500.0, 177.0, 69.0, 0.28),
    ] {
        let at = temperature(base.clone(), "temperature", celsius + 273.15);
        sample(&pack, "dynamic-young-modulus", PA, &at, e * 1e9);
        sample(&pack, "dynamic-shear-modulus", PA, &at, g * 1e9);
        sample(&pack, "dynamic-poisson-ratio", Dims::NONE, &at, nu);
        refuses(&pack, "young-modulus", PA, &at);
        refuses(
            &pack,
            "dynamic-young-modulus",
            PA,
            &at.clone().with("source-dynamic-method", 0.0).unwrap(),
        );
        refuses(
            &pack,
            "dynamic-young-modulus",
            PA,
            &at.with("source-frequency-known", 1.0).unwrap(),
        );
    }
    sample(
        &pack,
        "dynamic-young-modulus",
        PA,
        &temperature(base.clone(), "temperature", 423.15),
        197.5e9,
    );
    refuses(
        &pack,
        "dynamic-young-modulus",
        PA,
        &temperature(base.clone(), "temperature", 299.14),
    );
    refuses(
        &pack,
        "dynamic-young-modulus",
        PA,
        &temperature(base, "temperature", 773.16),
    );
    let nickel201 = compile("nickel-201-annealed-special-metals-2006", 5);
    let constants = alloy_reference()
        .with_quantity(
            "source-constant-temperature-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
        .with_quantity(
            "source-constant-heat-treatment-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
        .with_quantity(
            "source-elastic-method-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    sample(&nickel201, "tensile-modulus", PA, &constants, 207e9);
    refuses(&nickel201, "dynamic-young-modulus", PA, &constants);
    refuses(
        &nickel201,
        "tensile-modulus",
        PA,
        &constants.with("source-elastic-method-known", 1.0).unwrap(),
    );
}

#[test]
fn halide_crystals_preserve_thermal_mismatches_and_dielectric_frequency() {
    for (slug, rho, cp, k, k_kelvin, eps, rho_kelvin) in [
        (
            "potassium-bromide-crystran-reference",
            2753.0,
            435.0,
            4.816,
            319.0,
            4.9,
            None,
        ),
        (
            "sodium-fluoride-crystran-reference",
            2790.0,
            1088.0,
            3.746,
            273.0,
            6.0,
            Some(293.15),
        ),
        (
            "caesium-iodide-crystran-reference",
            4510.0,
            201.0,
            1.1,
            298.0,
            5.65,
            None,
        ),
    ] {
        let pack = compile(slug, 4);
        let base = flags(&[
            ("source-optical-crystal", 1.0),
            ("source-specimens-matched", 0.0),
        ]);
        let unknown = base
            .clone()
            .with_quantity(
                "source-test-temperature-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        let density_point = match rho_kelvin {
            Some(kelvin) => temperature(base.clone(), "temperature", kelvin),
            None => unknown.clone(),
        };
        sample(&pack, "density", RHO, &density_point, rho);
        sample(&pack, "specific-heat-capacity", CP, &unknown, cp);
        refuses(&pack, "specific-heat-capacity", CP, &base);
        let thermal = temperature(base.clone(), "temperature", k_kelvin);
        sample(&pack, "thermal-conductivity", K, &thermal, k);
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(base.clone(), "temperature", k_kelvin + 0.15),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &thermal.with("source-optical-crystal", 0.0).unwrap(),
        );
        let dielectric = unknown
            .with_quantity("frequency", QuantitySpec::dimensional(FREQUENCY), 1e6)
            .unwrap();
        sample(&pack, "relative-permittivity", Dims::NONE, &dielectric, eps);
        refuses(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &dielectric
                .clone()
                .with_quantity("frequency", QuantitySpec::dimensional(FREQUENCY), 1e3)
                .unwrap(),
        );
        refuses(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &dielectric
                .clone()
                .with("source-test-temperature-known", 1.0)
                .unwrap(),
        );
        refuses(&pack, "young-modulus", PA, &dielectric);
        if rho_kelvin.is_some() {
            refuses(
                &pack,
                "density",
                RHO,
                &temperature(base, "temperature", k_kelvin),
            );
        }
    }
}

#[test]
fn cobalt_keeps_twenty_and_twenty_five_degrees_distinct() {
    let pack = compile("cobalt-goodfellow-reference", 4);
    let base = flags(&[
        ("source-process-known", 0.0),
        ("source-specimens-matched", 0.0),
    ]);
    let twenty = temperature(base.clone(), "temperature", 293.15);
    let twenty_five = temperature(base.clone(), "temperature", 298.15);
    sample(&pack, "density", RHO, &twenty, 8900.0);
    sample(
        &pack,
        "electrical-resistivity",
        RESISTIVITY,
        &twenty,
        6.34e-8,
    );
    sample(&pack, "specific-heat-capacity", CP, &twenty_five, 456.0);
    refuses(&pack, "specific-heat-capacity", CP, &twenty);
    refuses(&pack, "density", RHO, &twenty_five);
    refuses(&pack, "electrical-resistivity", RESISTIVITY, &twenty_five);
    refuses(&pack, "thermal-conductivity", K, &twenty);
    let interval = base
        .with_quantity(
            "source-interval-reference",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    for kelvin in [273.15, 293.15, 298.15, 373.15] {
        sample(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(interval.clone(), "temperature", kelvin),
            100.0,
        );
    }
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &temperature(interval.clone(), "temperature", 373.16),
    );
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &temperature(interval, "temperature", 273.14),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &twenty.with("source-process-known", 1.0).unwrap(),
    );
}
