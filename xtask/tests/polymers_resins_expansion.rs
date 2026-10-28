//! G0/G3: real source compiler, SI values, directional/rate/cure domains and refusals.
//! Manufacturer typical values are neither design allowables nor validated constitutive laws.
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
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const HZ: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-polymers-{}-{now}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    let output = dir.join("material.fsmatpk");
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
    assert_eq!(pack.claims().claims_ordered().count(), count, "{slug}");
    pack
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (axis, value)| {
        axis_value(p, axis, Dims::NONE, *value)
    })
}

fn axis_value(p: QueryPoint, axis: &str, dims: Dims, value: f64) -> QueryPoint {
    p.with_quantity(axis, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn temperature(p: QueryPoint, axis: &str, kelvin: f64) -> QueryPoint {
    p.with_quantity(axis, semantic(QuantityKind::AbsoluteTemperature), kelvin)
        .unwrap()
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let quantity = if name.contains("thermal-conductivity") {
        assert_eq!(dims, K);
        semantic(QuantityKind::ThermalConductivity)
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, quantity)
}

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    assert!(
        (answer.evidence.value.value - expected).abs() <= 1e-11 * expected.abs().max(1.0),
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

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} {name}: unsupported point {p:?}",
        pack.pack_id()
    );
}

#[test]
fn g0_g3_pa66_pomh_pei_keep_thermal_values_rates_and_discrete_compression() {
    for (grade, rho, cp, k, tensile, flexural, compression, stress) in [
        (
            "tecamid-66",
            1150.0,
            1500.0,
            0.36,
            3500e6,
            3100e6,
            2700e6,
            81e6,
        ),
        (
            "tecaform-ad",
            1430.0,
            1300.0,
            0.43,
            3400e6,
            3600e6,
            2700e6,
            69e6,
        ),
        (
            "tecapei", 1280.0, 1200.0, 0.21, 3200e6, 3300e6, 2800e6, 92e6,
        ),
    ] {
        let pack = compile(&format!("ensinger-{grade}-natural"), 8);
        let condition = format!("source-grade-{grade}-natural");
        let p = flags(&[
            (&condition, 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
        ]);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "specific-heat-capacity", CP, &p, cp);
        sample(&pack, "thermal-conductivity", K, &p, k);
        refuses(
            &pack,
            "density",
            RHO,
            &p.clone().with(&condition, 0.0).unwrap(),
        );
        refuses(
            &pack,
            "density",
            RHO,
            &p.clone()
                .with("source-moisture-condition-known", 1.0)
                .unwrap(),
        );
        sample(
            &pack,
            "tensile-modulus",
            PA,
            &axis_value(p.clone(), "test-speed", SPEED, 1.0 / 60000.0),
            tensile,
        );
        sample(
            &pack,
            "flexural-modulus",
            PA,
            &axis_value(p.clone(), "test-speed", SPEED, 2.0 / 60000.0),
            flexural,
        );
        let compression_p = axis_value(p.clone(), "test-speed", SPEED, 5.0 / 60000.0);
        sample(
            &pack,
            "compression-modulus",
            PA,
            &compression_p,
            compression,
        );
        sample(
            &pack,
            "compressive-stress",
            PA,
            &axis_value(
                compression_p.clone(),
                "compressive-strain",
                Dims::NONE,
                0.05,
            ),
            stress,
        );
        refuses(
            &pack,
            "compressive-stress",
            PA,
            &axis_value(compression_p, "compressive-strain", Dims::NONE, 0.03),
        );
        refuses(&pack, "tensile-modulus", PA, &p);
        refuses(&pack, "young-modulus", PA, &p);
    }
}

#[test]
fn g0_g3_cast_pmma_has_real_thermal_profile_and_tabulated_temperature_frequency() {
    let pack = compile("plexiglas-gs-clear-reference", 8);
    let p = temperature(
        flags(&[
            ("source-grade-plexiglas-gs-0f00-233", 1.0),
            ("relative-humidity", 0.5),
            ("source-pressure-known", 0.0),
        ]),
        "temperature",
        296.15,
    );
    sample(&pack, "density", RHO, &p, 1190.0);
    sample(&pack, "specific-heat-capacity", CP, &p, 1470.0);
    sample(&pack, "thermal-conductivity", K, &p, 0.19);
    sample(&pack, "flexural-strength", PA, &p, 115e6);
    sample(&pack, "compressive-yield-stress", PA, &p, 110e6);
    sample(
        &pack,
        "tensile-modulus",
        PA,
        &axis_value(p.clone(), "test-speed", SPEED, 1.0 / 60000.0),
        3300e6,
    );
    let strength_p = axis_value(p.clone(), "test-speed", SPEED, 5.0 / 60000.0);
    for (t, strength) in [(233.15, 110e6), (296.15, 80e6), (343.15, 40e6)] {
        sample(
            &pack,
            "tensile-strength",
            PA,
            &strength_p.clone().with("temperature", t).unwrap(),
            strength,
        );
    }
    refuses(
        &pack,
        "tensile-strength",
        PA,
        &strength_p.with("temperature", 310.0).unwrap(),
    );
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &p.clone().with("temperature", 300.0).unwrap(),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &p.clone().with("relative-humidity", 0.8).unwrap(),
    );
    for (frequency, eps) in [(50.0, 3.6), (100000.0, 2.7)] {
        sample(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &axis_value(p.clone(), "frequency", HZ, frequency),
            eps,
        );
    }
    refuses(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &axis_value(p, "frequency", HZ, 1000.0),
    );
}

#[test]
fn g0_g3_mylar_preserves_md_td_and_separates_1000a_conductivity() {
    let pack = compile("mylar-92a-dupont-physical-thermal", 8);
    let p = flags(&[
        ("source-grade-mylar-92a", 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-test-rate-known", 0.0),
    ]);
    sample(&pack, "density", RHO, &p, 1390.0);
    sample(&pack, "specific-heat-capacity", CP, &p, 0.28 * 4184.0);
    for (direction, modulus, strength, elongation) in [
        (0.0, 490.0 * 9.80665e6, 20.0 * 9.80665e6, 1.16),
        (1.0, 510.0 * 9.80665e6, 24.0 * 9.80665e6, 0.91),
    ] {
        let directional = axis_value(p.clone(), "film-material-axis", Dims::NONE, direction);
        sample(&pack, "tensile-modulus", PA, &directional, modulus);
        sample(&pack, "tensile-strength", PA, &directional, strength);
        sample(
            &pack,
            "tensile-elongation",
            Dims::NONE,
            &directional,
            elongation,
        );
    }
    refuses(&pack, "tensile-modulus", PA, &p);
    refuses(
        &pack,
        "tensile-modulus",
        PA,
        &axis_value(p.clone(), "film-material-axis", Dims::NONE, 2.0),
    );
    refuses(&pack, "thermal-conductivity", K, &p);
    let thick = compile("mylar-1000a-dupont-conductivity", 1);
    let interval = temperature(
        temperature(
            flags(&[
                ("source-grade-mylar-1000a", 1.0),
                ("source-moisture-condition-known", 0.0),
                ("source-direction-known", 0.0),
            ]),
            "interval-temperature-lower",
            298.15,
        ),
        "interval-temperature-upper",
        348.15,
    );
    sample(
        &thick,
        "mean-thermal-conductivity",
        K,
        &interval,
        0.00037 * 4.184 * 100.0,
    );
    refuses(&thick, "thermal-conductivity", K, &interval);
    refuses(
        &thick,
        "mean-thermal-conductivity",
        K,
        &interval.with("interval-temperature-upper", 373.15).unwrap(),
    );
}

#[test]
fn g0_g3_epoxy_retains_mix_cure_rate_and_thermal_pressure() {
    let pack = compile("3m-dp420-black-cured", 4);
    let base = flags(&[
        ("source-grade-dp420-black", 1.0),
        ("base-to-accelerator-volume-ratio", 2.0),
        ("source-cure-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
    ]);
    let p = axis_value(
        axis_value(
            temperature(base.clone(), "temperature", 296.15),
            "test-speed",
            SPEED,
            10.0 / 60000.0,
        ),
        "cure-age",
        TIME,
        604800.0,
    );
    sample(&pack, "tensile-strength", PA, &p, 35.4e6);
    refuses(
        &pack,
        "tensile-strength",
        PA,
        &p.clone().with("cure-age", 86400.0).unwrap(),
    );
    refuses(
        &pack,
        "tensile-strength",
        PA,
        &p.with("base-to-accelerator-volume-ratio", 1.0).unwrap(),
    );
    let thermal = temperature(
        axis_value(base, "source-cure-age-known", Dims::NONE, 0.0),
        "temperature",
        323.15,
    )
    .with_quantity(
        "pressure",
        semantic(QuantityKind::Pressure),
        25.0 * 6894.757293168,
    )
    .unwrap();
    sample(&pack, "thermal-conductivity", K, &thermal, 0.22);
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &thermal.clone().with("pressure", 101325.0).unwrap(),
    );
    refuses(&pack, "density", RHO, &thermal);
}

#[test]
fn g0_g3_silicone_sealant_cure_and_discrete_dielectric_data() {
    let pack = compile("dowsil-732-seven-day-cure", 6);
    let p = temperature(
        axis_value(
            flags(&[
                ("source-grade-dowsil-732", 1.0),
                ("cure-relative-humidity", 0.5),
                ("source-test-temperature-known", 0.0),
                ("source-test-rate-known", 0.0),
            ]),
            "cure-age",
            TIME,
            604800.0,
        ),
        "cure-temperature",
        298.15,
    );
    sample(&pack, "specific-gravity", Dims::NONE, &p, 1.04);
    sample(&pack, "tensile-strength", PA, &p, 2.3e6);
    sample(&pack, "tensile-elongation", Dims::NONE, &p, 5.4);
    sample(
        &pack,
        "electrical-resistivity",
        Dims([3, 1, -3, 0, -2, 0]),
        &p,
        1.5e13,
    );
    for frequency in [100.0, 100000.0] {
        let electrical = axis_value(p.clone(), "frequency", HZ, frequency);
        sample(&pack, "relative-permittivity", Dims::NONE, &electrical, 2.8);
        sample(
            &pack,
            "dielectric-loss-tangent",
            Dims::NONE,
            &electrical,
            0.0015,
        );
    }
    refuses(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &axis_value(p.clone(), "frequency", HZ, 1000.0),
    );
    refuses(
        &pack,
        "tensile-strength",
        PA,
        &p.clone().with("cure-age", 86400.0).unwrap(),
    );
    refuses(&pack, "density", RHO, &p);
}

#[test]
fn g0_g3_melamine_foam_limits_and_ldpe_resin_are_not_promoted_to_laws() {
    let foam = compile("basotect-gplus-melamine-foam", 3);
    let p = flags(&[
        ("source-grade-basotect-gplus", 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-direction-known", 0.0),
    ]);
    sample(&foam, "density", RHO, &p, 9.0);
    sample(&foam, "thermal-conductivity-upper-limit", K, &p, 0.035);
    sample(&foam, "compression-hardness-lower-limit", PA, &p, 9000.0);
    refuses(&foam, "thermal-conductivity", K, &p);
    refuses(&foam, "specific-heat-capacity", CP, &p);
    let resin = compile("sabic-ldpe-2102tx00", 5);
    let p = flags(&[
        ("source-grade-ldpe-2102tx00", 1.0),
        ("unfoamed-resin", 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-test-rate-known", 0.0),
        ("source-dsc-rate-known", 0.0),
    ]);
    sample(&resin, "density", RHO, &p, 921.0);
    sample(&resin, "elastic-modulus-reference", PA, &p, 250e6);
    sample(
        &resin,
        "dsc-heat-of-fusion",
        Dims([2, 0, -2, 0, 0, 0]),
        &p,
        110000.0,
    );
    let transition_key = PropertyKey::with_quantity(
        "dsc-melting-temperature",
        semantic(QuantityKind::AbsoluteTemperature),
    );
    let transition = resin
        .claims()
        .query_typed(&transition_key, &p, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert!((transition.evidence.value.value - 381.15).abs() < 1e-12);
    refuses(
        &resin,
        "density",
        RHO,
        &p.clone().with("unfoamed-resin", 0.0).unwrap(),
    );
    refuses(&resin, "young-modulus", PA, &p);
}

#[test]
fn g0_g3_wood_glue_liquid_and_polyurethane_cure_remain_distinct_states() {
    let glue = compile("titebond-original-liquid-emulsion", 4);
    let liquid = flags(&[
        ("source-grade-titebond-original", 1.0),
        ("uncured-liquid-emulsion", 1.0),
        ("source-temperature-known", 0.0),
        ("source-shear-rate-known", 0.0),
        ("source-solids-basis-known", 0.0),
    ]);
    sample(
        &glue,
        "density",
        RHO,
        &liquid,
        9.2 * 0.45359237 / 0.003785411784,
    );
    sample(
        &glue,
        "apparent-dynamic-viscosity",
        Dims([-1, 1, -1, 0, 0, 0]),
        &liquid,
        3.2,
    );
    sample(&glue, "reported-solids-fraction", Dims::NONE, &liquid, 0.46);
    sample(&glue, "ph", Dims::NONE, &liquid, 4.0);
    refuses(
        &glue,
        "density",
        RHO,
        &liquid.clone().with("uncured-liquid-emulsion", 0.0).unwrap(),
    );
    refuses(&glue, "tensile-strength", PA, &liquid);

    let resin = compile("smooth-on-smooth-cast-300-cured", 8);
    let cured = temperature(
        axis_value(
            flags(&[
                ("source-grade-smooth-cast-300", 1.0),
                ("part-a-to-b-volume-ratio", 1.0),
                ("source-moisture-condition-known", 0.0),
                ("source-test-rate-known", 0.0),
                ("source-processing-history-known", 0.0),
            ]),
            "cure-age",
            TIME,
            7.0 * 86400.0,
        ),
        "temperature",
        296.15,
    );
    for (property, psi) in [
        ("tensile-strength", 3000.0),
        ("tensile-modulus", 139500.0),
        ("flexural-strength", 4510.0),
        ("flexural-modulus", 128000.0),
        ("compressive-strength", 4000.0),
        ("compression-modulus", 45800.0),
    ] {
        sample(&resin, property, PA, &cured, psi * 6894.757293168);
    }
    sample(&resin, "tensile-elongation", Dims::NONE, &cured, 0.05);
    sample(&resin, "linear-molding-shrinkage", Dims::NONE, &cured, 0.01);
    refuses(
        &resin,
        "tensile-strength",
        PA,
        &cured.clone().with("cure-age", 600.0).unwrap(),
    );
    refuses(
        &resin,
        "tensile-strength",
        PA,
        &cured
            .clone()
            .with("part-a-to-b-volume-ratio", 100.0 / 90.0)
            .unwrap(),
    );
    refuses(&resin, "young-modulus", PA, &cured);
}

#[test]
fn g0_g3_cured_potting_resins_keep_thermal_points_and_reported_diffusivity() {
    let diffusivity_dims = Dims([2, 0, -1, 0, 0, 0]);
    for (grade, count, ratio, rho, cp, k, alpha) in [
        ("834htc", 7, 10.0, 1700.0, 1300.0, 0.9, 4e-7),
        ("832tc", 7, 1.1, 1700.0, 1400.0, 0.7, 3e-7),
        ("832c", 7, 2.3, 1100.0, 2100.0, 0.3, 2e-7),
        ("832fx", 6, 1.22, 1100.0, 2700.0, 0.3, 1e-7),
        ("8810", 7, 1.7, 1200.0, 1700.0, 0.3, 1e-7),
    ] {
        let pack = compile(&format!("mg-chemicals-{grade}-cured-2026"), count);
        let grade_axis = format!("source-grade-mg-{grade}");
        let without_cure = temperature(
            flags(&[
                (&grade_axis, 1.0),
                ("cured-material", 1.0),
                ("specified-resin-to-hardener-mass-ratio", ratio),
                ("source-moisture-condition-known", 0.0),
                ("source-pressure-known", 0.0),
                ("source-specimens-matched", 0.0),
                ("source-measurement-direction-known", 0.0),
            ]),
            "temperature",
            298.15,
        );
        let point = axis_value(
            without_cure.clone(),
            "source-cure-schedule-known",
            Dims::NONE,
            0.0,
        );
        for (property, dims, expected) in [
            ("thermal-conductivity", K, k),
            ("specific-heat-capacity", CP, cp),
            ("thermal-diffusivity", diffusivity_dims, alpha),
        ] {
            sample(&pack, property, dims, &point, expected);
            refuses(&pack, property, dims, &without_cure);
            refuses(
                &pack,
                property,
                dims,
                &point.clone().with(&grade_axis, 0.0).unwrap(),
            );
            for outside in [297.15, 299.15] {
                refuses(
                    &pack,
                    property,
                    dims,
                    &point.clone().with("temperature", outside).unwrap(),
                );
            }
        }
        // Density carries an explicitly unknown test temperature; it must not
        // be accepted under a claim that its measurement temperature is known.
        refuses(&pack, "density", RHO, &point);
        let density_point = axis_value(
            point.clone(),
            "source-density-temperature-known",
            Dims::NONE,
            0.0,
        );
        sample(&pack, "density", RHO, &density_point, rho);
        refuses(
            &pack,
            "density",
            RHO,
            &density_point
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        if grade == "832c" {
            let reported = pack
                .claims()
                .query_typed(
                    &key("thermal-diffusivity", diffusivity_dims),
                    &point,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            assert_eq!(reported.evidence.value.value.to_bits(), 2e-7_f64.to_bits());
            assert!((reported.evidence.value.value - k / (rho * cp)).abs() > 5e-8);
            pack.claims().verify_receipt(&reported.receipt).unwrap();
        }
    }
}
