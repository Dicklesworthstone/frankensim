//! G0/G3 source regressions for missing polymer and cured-adhesive families.
//! These are producer comparison facts; no complete constitutive law is implied.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::sync::atomic::{AtomicU64, Ordering};
use std::{path::Path, process::Command};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const HZ: Dims = Dims([0, 0, -1, 0, 0, 0]);
const OHM_M: Dims = Dims([3, 1, -3, 0, -2, 0]);
const V_M: Dims = Dims([1, 1, -3, 0, -1, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-polymer-families-{}-{stamp}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    let out = dir.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&out)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count, "{slug}");
    pack
}

fn axis(p: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    p.with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn temp(p: QueryPoint, name: &str, kelvin: f64) -> QueryPoint {
    p.with_quantity(
        name,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        kelvin,
    )
    .unwrap()
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (name, value)| {
        axis(p, name, Dims::NONE, *value)
    })
}

fn stock(grade: &str) -> QueryPoint {
    flags(&[
        (grade, 1.0),
        ("source-test-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
    ])
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let qty = if name == "thermal-conductivity" {
        assert_eq!(dims, K);
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, qty)
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
        "{} {name}: accepted unsupported {p:?}",
        pack.pack_id()
    );
}

#[test]
fn g0_g3_pa12_retains_rate_strain_and_mean_expansion_conditions() {
    let pack = compile("ensinger-tecamid-12-natural", 13);
    let p = axis(
        stock("source-grade-tecamid-12-natural"),
        "immediately-after-machining",
        Dims::NONE,
        1.0,
    );
    sample(&pack, "density", RHO, &p, 1020.0);
    sample(&pack, "specific-heat-capacity", CP, &p, 1800.0);
    sample(&pack, "thermal-conductivity", K, &p, 0.30);
    let e = axis(p.clone(), "test-speed", SPEED, 1.6666666666666667e-5);
    sample(&pack, "tensile-modulus", PA, &e, 1800e6);
    refuses(&pack, "tensile-strength", PA, &e);
    let tension = axis(p.clone(), "test-speed", SPEED, 8.333333333333334e-4);
    sample(&pack, "tensile-strength", PA, &tension, 53e6);
    sample(&pack, "tensile-elongation", Dims::NONE, &tension, 2.0);
    sample(&pack, "tensile-yield-strain", Dims::NONE, &tension, 0.09);
    let bend = axis(p.clone(), "test-speed", SPEED, 3.3333333333333335e-5);
    sample(&pack, "flexural-modulus", PA, &bend, 1700e6);
    sample(&pack, "flexural-strength", PA, &bend, 68e6);
    let compression = axis(p.clone(), "test-speed", SPEED, 8.333333333333333e-5);
    sample(&pack, "compression-modulus", PA, &compression, 1600e6);
    for (strain, stress) in [(0.01, 13e6), (0.02, 24e6), (0.05, 55e6)] {
        sample(
            &pack,
            "compressive-stress",
            PA,
            &axis(
                compression.clone(),
                "compressive-strain",
                Dims::NONE,
                strain,
            ),
            stress,
        );
    }
    refuses(
        &pack,
        "compressive-stress",
        PA,
        &axis(compression, "compressive-strain", Dims::NONE, 0.03),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &axis(p.clone(), "immediately-after-machining", Dims::NONE, 0.0),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &axis(
            p.clone(),
            "source-moisture-condition-known",
            Dims::NONE,
            1.0,
        ),
    );
    for (end, alpha) in [(333.15, 15e-5), (373.15, 16e-5)] {
        let mean = temp(
            temp(p.clone(), "source-reference-temperature", 296.15),
            "source-range-end-temperature",
            end,
        );
        sample(
            &pack,
            "mean-longitudinal-thermal-expansion-coefficient-from-23c",
            Dims([0, 0, 0, -1, 0, 0]),
            &mean,
            alpha,
        );
    }
    refuses(
        &pack,
        "linear-thermal-expansion-coefficient",
        Dims([0, 0, 0, -1, 0, 0]),
        &p,
    );
}

#[test]
fn g0_g3_uhmwpe_retains_yield_and_indentation_roles_without_thermal_guesses() {
    let pack = compile("ensinger-tecafine-pe1000-natural-2017", 6);
    let p = stock("source-grade-tecafine-pe1000-natural");
    sample(&pack, "density", RHO, &p, 930.0);
    let mech = axis(p.clone(), "source-test-rate-known", Dims::NONE, 0.0);
    sample(&pack, "tensile-modulus", PA, &mech, 700e6);
    sample(&pack, "tensile-yield-strength", PA, &mech, 19e6);
    sample(&pack, "tensile-yield-strain", Dims::NONE, &mech, 0.11);
    sample(&pack, "ball-indentation-hardness", PA, &mech, 30e6);
    sample(
        &pack,
        "dielectric-breakdown-strength",
        V_M,
        &axis(p.clone(), "source-test-thickness-known", Dims::NONE, 0.0),
        44e6,
    );
    refuses(&pack, "tensile-strength", PA, &mech);
    refuses(&pack, "specific-heat-capacity", CP, &p);
    refuses(&pack, "thermal-conductivity", K, &p);
    refuses(
        &pack,
        "tensile-modulus",
        PA,
        &axis(mech, "source-test-rate-known", Dims::NONE, 1.0),
    );
}

#[test]
fn g0_g3_ptfe_preserves_compression_point_and_dielectric_film() {
    let pack = compile("ensinger-tecaflon-ptfe-natural", 8);
    let p = stock("source-grade-tecaflon-ptfe-natural");
    sample(&pack, "density", RHO, &p, 2150.0);
    sample(&pack, "thermal-conductivity", K, &p, 0.20);
    sample(&pack, "electrical-resistivity", OHM_M, &p, 1e15);
    let mech = axis(p.clone(), "source-test-rate-known", Dims::NONE, 0.0);
    sample(&pack, "tensile-strength", PA, &mech, 22e6);
    sample(&pack, "tensile-elongation", Dims::NONE, &mech, 2.2);
    sample(
        &pack,
        "compressive-stress",
        PA,
        &axis(mech.clone(), "compressive-strain", Dims::NONE, 0.01),
        5e6,
    );
    refuses(
        &pack,
        "compressive-stress",
        PA,
        &axis(mech, "compressive-strain", Dims::NONE, 0.02),
    );
    let film = axis(
        axis(p.clone(), "film-thickness", LENGTH, 0.000125),
        "test-medium-air",
        Dims::NONE,
        1.0,
    );
    sample(&pack, "dielectric-breakdown-strength", V_M, &film, 80e6);
    refuses(
        &pack,
        "dielectric-breakdown-strength",
        V_M,
        &axis(film, "film-thickness", LENGTH, 0.001),
    );
    refuses(&pack, "specific-heat-capacity", CP, &p);
    let mean = temp(
        temp(p, "source-reference-temperature", 296.15),
        "source-range-end-temperature",
        373.15,
    );
    sample(
        &pack,
        "mean-longitudinal-thermal-expansion-coefficient-from-23c",
        Dims([0, 0, 0, -1, 0, 0]),
        &mean,
        13e-5,
    );
}

#[test]
fn g0_g3_polyimide_separates_thermal_temperatures_and_frequency_samples() {
    let pack = compile("ensinger-tecasint-2011-natural-2025", 13);
    let p = stock("source-grade-tecasint-2011-natural");
    sample(&pack, "density", RHO, &p, 1380.0);
    sample(&pack, "specific-heat-capacity", CP, &p, 925.0);
    sample(
        &pack,
        "thermal-conductivity",
        K,
        &temp(p.clone(), "temperature", 313.15),
        0.22,
    );
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &temp(p.clone(), "temperature", 296.15),
    );
    sample(
        &pack,
        "electrical-resistivity",
        OHM_M,
        &temp(p.clone(), "temperature", 296.15),
        1e13,
    );
    refuses(
        &pack,
        "electrical-resistivity",
        OHM_M,
        &temp(p.clone(), "temperature", 313.15),
    );
    let e = axis(p.clone(), "test-speed", SPEED, 1.6666666666666667e-5);
    sample(&pack, "tensile-modulus", PA, &e, 3600e6);
    sample(&pack, "compression-modulus", PA, &e, 3430e6);
    let tensile = axis(p.clone(), "test-speed", SPEED, 8.333333333333334e-4);
    sample(&pack, "tensile-strength", PA, &tensile, 130e6);
    sample(&pack, "tensile-elongation", Dims::NONE, &tensile, 0.08);
    sample(
        &pack,
        "flexural-modulus",
        PA,
        &axis(p.clone(), "test-speed", SPEED, 3.3333333333333335e-5),
        3600e6,
    );
    let fast = axis(p.clone(), "test-speed", SPEED, 1.6666666666666666e-4);
    sample(&pack, "flexural-strength", PA, &fast, 177e6);
    sample(&pack, "compressive-strength", PA, &fast, 470e6);
    sample(
        &pack,
        "compressive-stress",
        PA,
        &axis(fast.clone(), "compressive-strain", Dims::NONE, 0.1),
        170e6,
    );
    refuses(
        &pack,
        "compressive-stress",
        PA,
        &axis(fast, "compressive-strain", Dims::NONE, 0.2),
    );
    for (hz, epsilon) in [(100.0, 3.5), (1000.0, 3.5), (10000.0, 3.4), (100000.0, 3.4)] {
        sample(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &axis(p.clone(), "frequency", HZ, hz),
            epsilon,
        );
    }
    refuses(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &axis(p.clone(), "frequency", HZ, 500.0),
    );
    refuses(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &axis(p, "frequency", HZ, 1e6),
    );
}

#[test]
fn g0_g3_acrylic_cure_and_frequency_conditions_exclude_liquid_and_joint_data() {
    let pack = compile("henkel-loctite-aa326-seven-day-cure-2020", 9);
    let p = temp(
        axis(
            flags(&[
                ("source-grade-loctite-aa326-2020", 1.0),
                ("cured-material", 1.0),
                ("source-cure-activation-known", 0.0),
                ("source-test-temperature-known", 0.0),
                ("source-moisture-condition-known", 0.0),
                ("source-pressure-known", 0.0),
            ]),
            "cure-age",
            Dims([0, 0, 1, 0, 0, 0]),
            604800.0,
        ),
        "cure-temperature",
        296.15,
    );
    sample(&pack, "specific-heat-capacity", CP, &p, 300.0);
    sample(&pack, "thermal-conductivity", K, &p, 0.1);
    sample(&pack, "electrical-resistivity", OHM_M, &p, 2e11);
    let mech = axis(p.clone(), "source-test-rate-known", Dims::NONE, 0.0);
    sample(&pack, "tensile-modulus", PA, &mech, 300e6);
    sample(&pack, "tensile-strength", PA, &mech, 34e6);
    sample(&pack, "tensile-elongation", Dims::NONE, &mech, 1.35);
    sample(
        &pack,
        "dielectric-breakdown-strength",
        V_M,
        &axis(p.clone(), "source-test-thickness-known", Dims::NONE, 0.0),
        30e6,
    );
    for (hz, epsilon, loss) in [(100.0, 5.6, 0.03), (1000.0, 5.3, 0.03), (1e6, 4.6, 0.04)] {
        let freq = axis(p.clone(), "frequency", HZ, hz);
        sample(&pack, "relative-permittivity", Dims::NONE, &freq, epsilon);
        sample(
            &pack,
            "dielectric-dissipation-factor",
            Dims::NONE,
            &freq,
            loss,
        );
    }
    refuses(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &axis(p.clone(), "frequency", HZ, 500.0),
    );
    refuses(
        &pack,
        "specific-heat-capacity",
        CP,
        &axis(p.clone(), "cure-age", Dims([0, 0, 1, 0, 0, 0]), 86400.0),
    );
    refuses(
        &pack,
        "specific-heat-capacity",
        CP,
        &temp(p.clone(), "cure-temperature", 298.15),
    );
    refuses(
        &pack,
        "specific-heat-capacity",
        CP,
        &axis(p.clone(), "cured-material", Dims::NONE, 0.0),
    );
    refuses(&pack, "density", RHO, &p);
    refuses(&pack, "lap-shear-strength", PA, &mech);
}

#[test]
fn g0_g3_neat_epoxies_keep_recipe_cure_and_specific_gravity_separate() {
    // Independent exact SI conversion: pound-force per square inch.
    let psi = 0.45359237 * 9.80665 / (0.0254_f64 * 0.0254);
    for (grade, ratio, sg, strain, values) in [
        (
            "105-205",
            5.19,
            1.18,
            0.034,
            [11418.0, 7846.0, 408000.0, 14112.0, 461000.0],
        ),
        (
            "105-207",
            3.64,
            1.15,
            0.038,
            [11043.0, 6748.0, 440000.0, 11324.0, 412000.0],
        ),
        (
            "gflex-650",
            1.2,
            1.11,
            0.327,
            [5268.0, 3440.0, 144000.0, 5192.0, 156000.0],
        ),
    ] {
        let slug = format!("west-system-{grade}-neat-cured-reference");
        let pack = compile(&slug, 7);
        let p = axis(
            flags(&[
                (&format!("source-grade-west-system-{grade}"), 1.0),
                ("cured-material", 1.0),
                ("specified-resin-hardener-mass-ratio", ratio),
                ("added-fillers", 0.0),
                ("added-fiber-reinforcement", 0.0),
                ("source-cure-temperature-known", 0.0),
                ("source-test-temperature-known", 0.0),
                ("source-moisture-condition-known", 0.0),
                ("source-pressure-known", 0.0),
            ]),
            "cure-age",
            Dims([0, 0, 1, 0, 0, 0]),
            1209600.0,
        );
        let specific = axis(
            p.clone(),
            "specific-gravity-reference-temperature-known",
            Dims::NONE,
            0.0,
        );
        sample(&pack, "specific-gravity", Dims::NONE, &specific, sg);
        refuses(&pack, "density", RHO, &specific);
        let mech = axis(p.clone(), "source-test-rate-known", Dims::NONE, 0.0);
        for (name, value) in [
            "compressive-yield-strength",
            "tensile-strength",
            "tensile-modulus",
            "flexural-strength",
            "flexural-modulus",
        ]
        .into_iter()
        .zip(values)
        {
            sample(&pack, name, PA, &mech, value * psi);
        }
        sample(&pack, "tensile-elongation", Dims::NONE, &mech, strain);
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &axis(mech.clone(), "cure-age", Dims([0, 0, 1, 0, 0, 0]), 604800.0),
        );
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &axis(mech.clone(), "added-fiber-reinforcement", Dims::NONE, 1.0),
        );
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &axis(mech.clone(), "added-fillers", Dims::NONE, 1.0),
        );
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &axis(
                mech.clone(),
                "specified-resin-hardener-mass-ratio",
                Dims::NONE,
                ratio + 0.1,
            ),
        );
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &axis(mech, "cured-material", Dims::NONE, 0.0),
        );
        refuses(&pack, "thermal-conductivity", K, &p);
    }
}
