//! Named polyolefin measurements and directly measured volumetric thermal inputs.
//! The source's distinct temperatures and dependent quantities remain explicit.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, InterpolationPolicy, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const DIFFUSIVITY: Dims = Dims([2, 0, -1, 0, 0, 0]);
const CAPACITY: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 1, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const MTPS: &str = "polypropylene-mtps-thermtest-26c-reference";

fn compile(slug: &str, expected_claims: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "polymer-thermal-completion-{}-{nonce}-{slug}.fsmatpk",
        std::process::id()
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), expected_claims);
    pack
}

fn axis(p: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    p.with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn flags(mut p: QueryPoint, values: &[(&str, f64)]) -> QueryPoint {
    for &(name, value) in values {
        p = axis(p, name, Dims::NONE, value);
    }
    p
}

fn temperature(p: QueryPoint, name: &str, value: f64) -> QueryPoint {
    p.with_quantity(
        name,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        value,
    )
    .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    PropertyKey::with_quantity(
        name,
        if matches!(
            name,
            "thermal-conductivity" | "source-derived-thermal-conductivity"
        ) {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        },
    )
}

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint, expected: f64) -> f64 {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|error| panic!("{} {name}: {error:?}", pack.pack_id()));
    assert!((answer.evidence.value.value - expected).abs() <= expected.abs() * 1e-12);
    assert_eq!(answer.evidence.value.dims, dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer.evidence.value.value
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} must refuse {name} at {p:?}",
        pack.pack_id()
    );
}

fn grade(slug: &str, temperature_known: f64) -> QueryPoint {
    flags(
        QueryPoint::new(),
        &[
            (&format!("source-grade-{slug}"), 1.0),
            ("source-test-temperature-known", temperature_known),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-specimens-matched", 0.0),
            ("source-properties-share-temperature", 0.0),
        ],
    )
}

fn film(slug: &str, press_temperature: f64) -> QueryPoint {
    let p = flags(grade(slug, 0.0), &[("source-is-hotpressed-film", 1.0)]);
    let p = temperature(p, "source-film-hotpress-temperature", press_temperature);
    let p = axis(p, "source-film-hotpress-pressure", PRESSURE, 1e7);
    axis(p, "source-film-hotpress-duration", TIME, 180.0)
}

#[test]
fn named_polyolefins_keep_density_calorimetry_and_tls_coordinates_distinct() {
    for (slug, rho, cp, d, k, hotpress) in [
        (
            "hdpe-dr720-abbasgholi2022",
            958.0,
            1760.0,
            3.20e-7,
            0.54,
            453.15,
        ),
        (
            "lldpe-22b02-abbasgholi2022",
            922.0,
            2030.0,
            2.19e-7,
            0.41,
            453.15,
        ),
        (
            "ldpe-2420e02-abbasgholi2022",
            921.0,
            2090.0,
            1.92e-7,
            0.37,
            453.15,
        ),
        (
            "pp-zr340r-abbasgholi2022",
            903.0,
            1670.0,
            1.46e-7,
            0.22,
            468.15,
        ),
    ] {
        let pack = compile(slug, 4);
        let density_point = temperature(grade(slug, 1.0), "temperature", 296.15);
        let cp_point = axis(
            temperature(grade(slug, 1.0), "temperature", 298.15),
            "source-cp-heating-rate",
            RATE,
            1.0 / 3.0,
        );
        let diffusivity_point = film(slug, hotpress);
        sample(&pack, "density", RHO, &density_point, rho);
        sample(&pack, "specific-heat-capacity", CP, &cp_point, cp);
        sample(
            &pack,
            "thermal-diffusivity",
            DIFFUSIVITY,
            &diffusivity_point,
            d,
        );
        refuses(&pack, "density", RHO, &cp_point);
        refuses(&pack, "specific-heat-capacity", CP, &density_point);
        refuses(
            &pack,
            "thermal-diffusivity",
            DIFFUSIVITY,
            &film(slug, hotpress + 1.0),
        );
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &axis(cp_point.clone(), "source-cp-heating-rate", RATE, 1.0 / 6.0),
        );

        let derived_point = temperature(
            temperature(
                flags(
                    diffusivity_point.clone(),
                    &[
                        ("engineering-reference-model", 1.0),
                        ("source-conductivity-independently-measured", 0.0),
                    ],
                ),
                "density-reference-temperature",
                296.15,
            ),
            "specific-heat-reference-temperature",
            298.15,
        );
        sample(
            &pack,
            "source-derived-thermal-conductivity",
            K,
            &derived_point,
            k,
        );
        // Equation 4 reproduces the printed two-decimal k within rounding;
        // this does not manufacture a common-temperature experimental state.
        assert!((rho * cp * d - k).abs() < 0.005);
        let derived = pack
            .claims()
            .query_typed(
                &key("source-derived-thermal-conductivity", K),
                &derived_point,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(
            pack.claims()
                .claim(derived.receipt.selected)
                .unwrap()
                .observations
                .len(),
            5
        );
        refuses(&pack, "thermal-conductivity", K, &derived_point);
        refuses(
            &pack,
            "source-derived-thermal-conductivity",
            K,
            &axis(
                derived_point.clone(),
                "source-conductivity-independently-measured",
                Dims::NONE,
                1.0,
            ),
        );
        refuses(
            &pack,
            "source-derived-thermal-conductivity",
            K,
            &axis(
                derived_point.clone(),
                "engineering-reference-model",
                Dims::NONE,
                0.0,
            ),
        );
        refuses(
            &pack,
            "source-derived-thermal-conductivity",
            K,
            &temperature(derived_point, "density-reference-temperature", 298.15),
        );
        for (name, dims, p) in [
            ("density", RHO, density_point),
            ("specific-heat-capacity", CP, cp_point),
            ("thermal-diffusivity", DIFFUSIVITY, diffusivity_point),
        ] {
            refuses(
                &pack,
                name,
                dims,
                &axis(
                    p.clone(),
                    "source-properties-share-temperature",
                    Dims::NONE,
                    1.0,
                ),
            );
            refuses(
                &pack,
                name,
                dims,
                &axis(p, &format!("source-grade-{slug}"), Dims::NONE, 0.0),
            );
        }
    }
}

fn mtps(mean: bool, repeat: usize) -> QueryPoint {
    let p = flags(
        QueryPoint::new(),
        &[
            ("engineering-reference-model", 1.0),
            ("source-sample-thermtest-pp", 1.0),
            ("source-grade-known", 0.0),
            ("source-processing-history-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-anisotropy-direction-known", 0.0),
            ("source-same-specimen", 1.0),
            ("source-method-mtps-asymmetric", 1.0),
            ("source-statistic-is-mean", if mean { 1.0 } else { 0.0 }),
            ("source-repeat-number", repeat as f64),
        ],
    );
    let p = axis(p, "source-test-duration", TIME, 80.0);
    let p = axis(p, "source-sensor-radius", LENGTH, 0.0128);
    let p = axis(p, "source-sensor-power", POWER, 0.030);
    temperature(p, "temperature", 299.15)
}

#[test]
fn same_specimen_volumetric_measurements_preserve_all_runs_and_printed_statistics() {
    let pack = compile(MTPS, 21);
    assert!(
        pack.claims()
            .claims_ordered()
            .all(|(_, c)| c.interpolation == InterpolationPolicy::TabulatedOnly)
    );
    let runs = [
        (0.270, 1.72e-7, 1570000.0),
        (0.280, 1.81e-7, 1544000.0),
        (0.282, 1.83e-7, 1547000.0),
        (0.286, 1.90e-7, 1510000.0),
        (0.283, 1.87e-7, 1510000.0),
    ];
    for (i, &(k, d, cv)) in runs.iter().enumerate() {
        let p = mtps(false, i + 1);
        sample(&pack, "thermal-conductivity", K, &p, k);
        sample(&pack, "thermal-diffusivity", DIFFUSIVITY, &p, d);
        sample(&pack, "volumetric-heat-capacity", CAPACITY, &p, cv);
        // The printed instrument outputs are correlated and rounded. Retain
        // all three instead of replacing a source value with an exact identity.
        assert!((k / cv - d).abs() / d < 0.01);
    }
    let p = mtps(true, 0);
    for (name, dims, mean, sd) in [
        ("thermal-conductivity", K, 0.280, 0.006),
        ("thermal-diffusivity", DIFFUSIVITY, 1.83e-7, 7e-9),
        ("volumetric-heat-capacity", CAPACITY, 1536000.0, 26000.0),
    ] {
        sample(&pack, name, dims, &p, mean);
        sample(
            &pack,
            &format!("{name}-reported-standard-deviation"),
            dims,
            &p,
            sd,
        );
        refuses(&pack, name, dims, &mtps(true, 1));
        refuses(&pack, name, dims, &mtps(false, 0));
    }
    let printed_mean_k = sample(&pack, "thermal-conductivity", K, &p, 0.280);
    assert!((runs.iter().map(|r| r.0).sum::<f64>() / 5.0 - printed_mean_k).abs() < 0.0005);
    assert!(
        pack.claims()
            .query_envelope_typed(
                &key("thermal-conductivity", K),
                &mtps(false, 1),
                &mtps(false, 5),
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .is_err()
    );
}

#[test]
fn volumetric_source_requires_its_actual_experiment_and_does_not_invent_mass_properties() {
    let pack = compile(MTPS, 21);
    let p = mtps(true, 0);
    for (name, dims) in [
        ("thermal-conductivity", K),
        ("volumetric-heat-capacity", CAPACITY),
    ] {
        for (axis_name, axis_dims, changed) in [
            ("engineering-reference-model", Dims::NONE, 0.0),
            ("source-grade-known", Dims::NONE, 1.0),
            ("source-processing-history-known", Dims::NONE, 1.0),
            ("source-same-specimen", Dims::NONE, 0.0),
            ("source-test-duration", TIME, 40.0),
            ("source-sensor-power", POWER, 0.060),
        ] {
            refuses(
                &pack,
                name,
                dims,
                &axis(p.clone(), axis_name, axis_dims, changed),
            );
        }
        for t in [298.15, 299.14, 299.16, 300.15] {
            refuses(&pack, name, dims, &temperature(p.clone(), "temperature", t));
        }
        refuses(
            &pack,
            name,
            dims,
            &axis(p.clone(), "temperature", TEMP, 299.15),
        );
    }
    refuses(&pack, "density", RHO, &p);
    refuses(&pack, "specific-heat-capacity", CP, &p);
    refuses(
        &pack,
        "linear-thermal-expansion",
        Dims([0, 0, 0, -1, 0, 0]),
        &p,
    );
}
