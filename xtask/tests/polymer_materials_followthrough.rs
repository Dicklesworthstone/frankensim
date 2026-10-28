//! Source-backed cured silicones, extruded PP and closed-cell EPDM insulation.
//! Bounds, conversion assumptions and authored local-conductivity models stay explicit.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const HZ: Dims = Dims([0, 0, -1, 0, 0, 0]);
const SECOND: Dims = Dims([0, 0, 1, 0, 0, 0]);
const METRE: Dims = Dims([1, 0, 0, 0, 0, 0]);
const OHM_M: Dims = Dims([3, 1, -3, 0, -2, 0]);
const V_M: Dims = Dims([1, 1, -3, 0, -1, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-polymer-followthrough-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let output = dir.join("material.fsmatpk");
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
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count, "{slug}");
    pack
}
fn axis(p: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    p.with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}
fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (name, value)| {
        axis(p, name, Dims::NONE, *value)
    })
}
fn temp(p: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    p.with_quantity(
        axis,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        value,
    )
    .unwrap()
}
fn key(name: &str, dims: Dims) -> PropertyKey {
    let spec = if name.contains("thermal-conductivity") {
        assert_eq!(dims, K);
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, spec)
}
fn sample(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    assert!((answer.evidence.value.value - expected).abs() <= 1e-11 * expected.abs().max(1e-6));
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
        "{} must refuse {name}",
        pack.pack_id()
    );
}
fn silicone(grade: &str) -> QueryPoint {
    let p = flags(&[
        (&format!("source-grade-momentive-rtv{grade}"), 1.0),
        ("mix-resin-to-hardener-mass-ratio", 10.0),
        ("matched-kit-components", 1.0),
        ("source-test-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-test-rate-known", 0.0),
        ("source-optical-wavelength-known", 0.0),
        ("thermochemical-calorie-conversion", 1.0),
    ]);
    temp(
        axis(p, "cure-time", SECOND, 3600.0),
        "cure-temperature",
        373.15,
    )
}

#[test]
fn g0_g3_cured_silicone_thermal_and_dielectric_data_keep_their_actual_cure() {
    for (grade, shore, index, epsilon, loss, resistivity) in [
        ("615", 44.0, 1.406, 2.7, 0.0006, 1.8e13),
        ("655", 45.0, 1.430, 2.69, 0.0004, 1.2e13),
    ] {
        let pack = compile(&format!("momentive-rtv{grade}-100c-cure-2021"), 11);
        let p = silicone(grade);
        sample(&pack, "tensile-strength", PA, &p, 65.0 * 9.80665 / 1e-4);
        sample(&pack, "tensile-elongation", Dims::NONE, &p, 1.2);
        sample(&pack, "hardness-shore-a", Dims::NONE, &p, shore);
        sample(&pack, "reported-cure-shrinkage", Dims::NONE, &p, 0.002);
        sample(&pack, "refractive-index", Dims::NONE, &p, index);
        sample(&pack, "thermal-conductivity", K, &p, 0.00045 * 4.184 / 0.01);
        sample(&pack, "specific-heat-capacity", CP, &p, 0.3 * 4.184 / 0.001);
        sample(&pack, "electrical-resistivity", OHM_M, &p, resistivity);
        let electrical = axis(
            axis(p.clone(), "frequency", HZ, 1000.0),
            "test-thickness",
            METRE,
            0.0019,
        );
        sample(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &electrical,
            epsilon,
        );
        sample(
            &pack,
            "dielectric-loss-tangent",
            Dims::NONE,
            &electrical,
            loss,
        );
        sample(&pack, "dielectric-strength", V_M, &electrical, 19.7e6);
        refuses(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &axis(electrical.clone(), "frequency", HZ, 100.0),
        );
        refuses(
            &pack,
            "dielectric-strength",
            V_M,
            &axis(electrical, "test-thickness", METRE, 0.001),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &axis(
                p.clone(),
                "thermochemical-calorie-conversion",
                Dims::NONE,
                0.0,
            ),
        );
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &temp(p.clone(), "cure-temperature", 298.15),
        );
        refuses(
            &pack,
            "tensile-strength",
            PA,
            &axis(p.clone(), "cure-time", SECOND, 7.0 * 86400.0),
        );
        refuses(
            &pack,
            "tensile-strength",
            PA,
            &axis(
                p.clone(),
                "mix-resin-to-hardener-mass-ratio",
                Dims::NONE,
                1.0,
            ),
        );
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &axis(p.clone(), "source-test-temperature-known", Dims::NONE, 1.0),
        );
        // Component-A specific gravity cannot silently provide cured density.
        refuses(&pack, "density", RHO, &p);
        refuses(&pack, "youngs-modulus", PA, &p);
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &silicone(if grade == "615" { "655" } else { "615" }),
        );
    }
}

#[test]
fn g0_g3_polypropylene_heat_capacity_does_not_promote_guideline_bounds_to_measurements() {
    let pack = compile("roechling-polystone-p-h-grey-extruded-2025", 9);
    let p = flags(&[
        ("source-grade-polystone-p-h-grey-extruded-2025", 1.0),
        ("source-extruded-product", 1.0),
        ("source-test-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-test-rate-known", 0.0),
    ]);
    sample(&pack, "specific-heat-capacity", CP, &p, 1700.0);
    sample(&pack, "thermal-conductivity", K, &p, 0.20);
    sample(
        &pack,
        "density-strict-lower-guideline-bound",
        RHO,
        &p,
        900.0,
    );
    sample(
        &pack,
        "tensile-modulus-strict-lower-guideline-bound",
        PA,
        &p,
        1.5e9,
    );
    sample(
        &pack,
        "tensile-yield-strength-strict-lower-guideline-bound",
        PA,
        &p,
        30e6,
    );
    sample(
        &pack,
        "tensile-elongation-strict-lower-guideline-bound",
        Dims::NONE,
        &p,
        0.5,
    );
    sample(
        &pack,
        "notched-charpy-impact-strict-lower-guideline-bound",
        Dims([0, 1, -2, 0, 0, 0]),
        &p,
        6000.0,
    );
    sample(
        &pack,
        "water-absorption-strict-upper-guideline-bound",
        Dims::NONE,
        &p,
        0.001,
    );
    sample(
        &pack,
        "dielectric-loss-tangent",
        Dims::NONE,
        &axis(p.clone(), "frequency", HZ, 1e6),
        0.00019,
    );
    refuses(
        &pack,
        "dielectric-loss-tangent",
        Dims::NONE,
        &axis(p.clone(), "frequency", HZ, 1000.0),
    );
    refuses(
        &pack,
        "specific-heat-capacity",
        CP,
        &axis(p.clone(), "source-extruded-product", Dims::NONE, 0.0),
    );
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &axis(p.clone(), "source-test-temperature-known", Dims::NONE, 1.0),
    );
    for (name, dims) in [
        ("density", RHO),
        ("tensile-modulus", PA),
        ("tensile-yield-strength", PA),
    ] {
        refuses(&pack, name, dims, &p);
    }
}

#[test]
fn g0_g3_epdm_foam_forms_keep_declarations_separate_from_local_transport_references() {
    for (form, offset, limit, mu) in [
        ("tube", 0.0, 0.042, 4000.0),
        ("sheet", 0.003, 0.045, 3000.0),
    ] {
        let pack = compile(&format!("armacell-ht-armaflex-{form}-reference-2026"), 3);
        let p = flags(&[
            (&format!("source-grade-ht-armaflex-{form}"), 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-gas-pressure-known", 0.0),
            ("source-test-temperature-difference-known", 0.0),
            ("source-test-temperature-known", 0.0),
        ]);
        sample(
            &pack,
            "declared-thermal-conductivity-upper-limit",
            K,
            &temp(p.clone(), "mean-test-temperature", 313.15),
            limit,
        );
        sample(
            &pack,
            "vapour-diffusion-resistance-factor-lower-limit",
            Dims::NONE,
            &p,
            mu,
        );
        let reference = axis(
            axis(p.clone(), "engineering-reference-model", Dims::NONE, 1.0),
            "declared-conductivity-limit-as-local-model",
            Dims::NONE,
            1.0,
        );
        for celsius in [0.0, 20.0, 40.0, 60.0, 80.0, 100.0] {
            let polynomial =
                (36.92 + 0.125 * celsius + 0.0008 * (celsius - 30.0_f64).powi(2)) / 1000.0;
            sample(
                &pack,
                "thermal-conductivity",
                K,
                &temp(reference.clone(), "temperature", celsius + 273.15),
                polynomial + offset,
            );
        }
        // The stored local law is piecewise-linear, not the producer quadratic.
        sample(
            &pack,
            "thermal-conductivity",
            K,
            &temp(reference.clone(), "temperature", 303.15),
            (0.03950 + 0.04200) / 2.0 + offset,
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temp(p.clone(), "temperature", 313.15),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temp(reference.clone(), "temperature", 373.16),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temp(reference, "temperature", 273.14),
        );
        refuses(
            &pack,
            "declared-thermal-conductivity-upper-limit",
            K,
            &temp(p.clone(), "mean-test-temperature", 303.15),
        );
        refuses(&pack, "density", RHO, &p);
        refuses(&pack, "specific-heat-capacity", CP, &p);
    }
}
