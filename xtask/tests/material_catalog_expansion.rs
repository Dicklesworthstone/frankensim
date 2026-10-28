//! G0/G3 source-to-pack/query checks for the October 2026 material additions.
//! These establish source transport and bounded evaluation, not experimental validation.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    MaterialAnswer, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(slug: &str, count: usize) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-material-expansion-{}-{time}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    let output = dir.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
        .output()
        .expect("run the actual material compiler");
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes, "{slug}: canonical round trip");
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count, "{slug}");
    pack
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (name, value)| {
        p.with_quantity(*name, QuantitySpec::dimensional(Dims::NONE), *value)
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

fn sample(
    pack: &NormalizedPack,
    name: &str,
    dims: Dims,
    point: &QueryPoint,
    expected: f64,
) -> MaterialAnswer {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    let actual = answer.evidence.value.value;
    assert!(
        (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
        "{} {name}: {actual} != {expected}",
        pack.pack_id()
    );
    assert_eq!(answer.evidence.value.dims, dims);
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    assert!(!answer.receipt.axis_quantities.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer
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
fn g0_g3_clear_wood_density_stiffness_and_moisture_are_preserved() {
    for (species, count, sg, bending_mpa, rupture_kpa) in [
        ("pine-loblolly", 19, 0.51, 12300.0, 88000.0),
        ("pine-longleaf", 19, 0.59, 13700.0, 100000.0),
        ("pine-slash", 17, 0.59, 13700.0, 112000.0),
        ("pine-ponderosa", 19, 0.40, 8900.0, 65000.0),
        ("oak-northern-red", 10, 0.63, 12500.0, 99000.0),
        ("oak-white", 10, 0.68, 12300.0, 105000.0),
        ("cedar-northern-white", 19, 0.31, 5500.0, 45000.0),
        ("fir-subalpine", 18, 0.32, 8900.0, 59000.0),
    ] {
        let pack = compile(&format!("{species}-fpl-gtr282"), count);
        let p = flags(&[
            ("moisture-content", 0.12),
            ("clear-specimen", 1.0),
            ("source-test-temperature-known", 0.0),
        ]);
        sample(&pack, "density", RHO, &p, sg * 1000.0 * 1.12);
        sample(
            &pack,
            "modulus_of_elasticity_bending",
            PA,
            &p,
            bending_mpa * 1e6,
        );
        sample(
            &pack,
            "young_modulus_longitudinal",
            PA,
            &p,
            bending_mpa * 1.1e6,
        );
        sample(&pack, "bending_mor", PA, &p, rupture_kpa * 1000.0);
        for bad in [
            p.clone().with("moisture-content", 0.20).unwrap(),
            p.clone().with("clear-specimen", 0.0).unwrap(),
            QueryPoint::new(),
        ] {
            refuses(&pack, "density", RHO, &bad);
        }
        refuses(&pack, "young_modulus", PA, &p); // no isotropic substitution
        if species.starts_with("oak-") {
            refuses(&pack, "er_over_el", Dims::NONE, &p);
        }
    }
}

#[test]
fn g0_g3_elements_minerals_and_alloys_keep_reference_model_guards() {
    for (slug, count, t, rho, cp, k) in [
        (
            "diamond-crystran-thermal-reference",
            3,
            273.0,
            3510.0,
            502.0,
            2600.0,
        ),
        (
            "fluorite-caf2-crystran-reference",
            3,
            293.0,
            3180.0,
            854.0,
            9.71,
        ),
        (
            "fused-silica-crystran-reference",
            5,
            298.15,
            2203.0,
            703.0,
            1.38,
        ),
        (
            "germanium-crystran-thermal-reference",
            3,
            293.0,
            5330.0,
            310.0,
            58.62,
        ),
        (
            "halite-nacl-crystran-reference",
            3,
            273.0,
            2170.0,
            854.0,
            6.5,
        ),
        (
            "periclase-mgo-crystran-reference",
            3,
            273.0,
            3580.0,
            877.0,
            42.0,
        ),
        (
            "aluminum-5005-austral-wright-reference",
            7,
            293.15,
            2700.0,
            900.0,
            201.0,
        ),
        (
            "aluminum-5052-austral-wright-reference",
            7,
            293.15,
            2680.0,
            900.0,
            138.0,
        ),
        (
            "bronze-c93200-cda-reference",
            4,
            293.15,
            8912.929316685405,
            376.812,
            58.15268479007874,
        ),
        (
            "bronze-c95400-cda-reference",
            4,
            293.15,
            7445.89436704464,
            418.68,
            58.67190518999016,
        ),
        (
            "cupronickel-c70600-cda-reference",
            6,
            293.15,
            8940.609221395608,
            376.812,
            44.99910132565617,
        ),
        (
            "magnesium-espi-reference",
            4,
            298.15,
            1738.0,
            1016.712,
            156.0,
        ),
        (
            "molybdenum-attl-reference",
            5,
            293.15,
            10220.0,
            240.0,
            140.0,
        ),
        ("tantalum-efinea-reference", 4, 293.15, 16650.0, 140.0, 57.5),
    ] {
        let pack = compile(slug, count);
        let p = temperature(
            flags(&[
                ("engineering-reference-model", 1.0),
                ("source-temperatures-matched", 0.0),
            ]),
            "temperature",
            t,
        );
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "specific-heat-capacity", CP, &p, cp);
        sample(&pack, "thermal-conductivity", K, &p, k);
        for bad in [
            p.clone().with("engineering-reference-model", 0.0).unwrap(),
            p.clone().with("source-temperatures-matched", 1.0).unwrap(),
            p.clone().with("temperature", t - 0.1).unwrap(),
            QueryPoint::new(),
        ] {
            refuses(&pack, "density", RHO, &bad);
        }
        if count == 3 {
            refuses(&pack, "young-modulus", PA, &p);
        }
        if slug == "aluminum-5052-austral-wright-reference" {
            sample(&pack, "young-modulus", PA, &p, 69.3e9);
            sample(
                &pack,
                "electrical-resistivity",
                Dims([3, 1, -3, 0, -2, 0]),
                &p,
                5e-8,
            );
        }
    }
}

fn plastic_point(grade: &str) -> QueryPoint {
    let mut p = flags(&[
        (grade, 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
    ]);
    if grade.contains("tecamid") {
        p = p
            .with_quantity(
                "immediately-after-machining",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
    }
    p
}

#[test]
fn g0_g3_plastics_keep_grade_test_speed_and_discrete_compression() {
    for (grade, rho, cp, k, modulus) in [
        ("tecaform-ah-natural", 1410.0, 1400.0, 0.39, 2.8e9),
        ("tecamid-6-natural", 1140.0, 1600.0, 0.37, 3.3e9),
        ("tecapet-white", 1360.0, 1200.0, 0.31, 3.1e9),
        ("tecaflon-pvdf-natural", 1780.0, 1300.0, 0.25, 2.2e9),
    ] {
        let pack = compile(&format!("ensinger-{grade}"), 8);
        let grade_axis = format!("source-grade-{grade}");
        let p = plastic_point(&grade_axis);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "specific-heat-capacity", CP, &p, cp);
        sample(&pack, "thermal-conductivity", K, &p, k);
        let tension = p
            .clone()
            .with_quantity(
                "test-speed",
                QuantitySpec::dimensional(SPEED),
                1.6666666666666667e-5,
            )
            .unwrap();
        sample(&pack, "tensile-modulus", PA, &tension, modulus);
        refuses(&pack, "tensile-modulus", PA, &p);
        refuses(
            &pack,
            "density",
            RHO,
            &p.clone().with(&grade_axis, 0.0).unwrap(),
        );
        if grade == "tecaform-ah-natural" {
            let compression = p
                .with_quantity(
                    "test-speed",
                    QuantitySpec::dimensional(SPEED),
                    8.333333333333333e-5,
                )
                .unwrap();
            for (strain, stress) in [(0.01, 20e6), (0.02, 35e6), (0.05, 68e6)] {
                let at = compression
                    .clone()
                    .with_quantity(
                        "compressive-strain",
                        QuantitySpec::dimensional(Dims::NONE),
                        strain,
                    )
                    .unwrap();
                sample(&pack, "compressive-stress", PA, &at, stress);
            }
            let between = compression
                .with_quantity(
                    "compressive-strain",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.03,
                )
                .unwrap();
            refuses(&pack, "compressive-stress", PA, &between);
        }
    }
}

#[test]
fn g0_g3_foam_grade_and_cured_resin_conditions_survive_compilation() {
    for (grade, rho, modulus) in [("51", 52.0, 43e6), ("110", 110.0, 120e6)] {
        let pack = compile(&format!("evonik-rohacell-{grade}-igf-2022"), 7);
        let grade_axis = format!("source-grade-rohacell-{grade}-igf-2022");
        let p = flags(&[
            (&grade_axis, 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-test-rate-known", 0.0),
            ("full-sheet-at-least-10mm", 1.0),
        ]);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "compression-modulus", PA, &p, modulus);
        refuses(
            &pack,
            "density",
            RHO,
            &p.clone().with("full-sheet-at-least-10mm", 0.0).unwrap(),
        );
    }
    for (slug, grade, modulus) in [
        (
            "west-system-105-205-neat-cured",
            "source-grade-west-105-205",
            2813060975.612544,
        ),
        (
            "west-system-gflex-650-neat-cured",
            "source-grade-west-gflex-650",
            992845050.216192,
        ),
    ] {
        let pack = compile(slug, 7);
        let p = flags(&[
            (grade, 1.0),
            ("room-temperature-cure", 1.0),
            ("neat-unreinforced-resin", 1.0),
            ("source-temperature-known", 0.0),
            ("source-test-rate-known", 0.0),
        ])
        .with_quantity(
            "cure-age",
            QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])),
            1209600.0,
        )
        .unwrap();
        sample(&pack, "tensile-modulus", PA, &p, modulus);
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &p.clone().with("cure-age", 604800.0).unwrap(),
        );
        refuses(&pack, "density", RHO, &p); // specific gravity alone has no declared water reference
    }
}

#[test]
fn g0_g3_insulation_curves_and_density_temperature_surface_keep_typed_axes() {
    let rock = compile("rockwool-prorox-sl960-en-2020", 3);
    let p = temperature(
        flags(&[
            ("source-grade-prorox-sl960-en-2020", 1.0),
            ("source-temperature-known", 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
        ]),
        "temperature",
        348.15,
    );
    sample(&rock, "thermal-conductivity", K, &p, 0.0425);
    refuses(
        &rock,
        "thermal-conductivity",
        K,
        &p.clone().with("temperature", 300.0).unwrap(),
    );
    refuses(
        &rock,
        "thermal-conductivity",
        K,
        &p.with_quantity("temperature", QuantitySpec::dimensional(TEMP), 348.15)
            .unwrap(),
    );

    let eps = compile("eps-nist-srm1453-2013", 1);
    let p = temperature(
        flags(&[
            ("source-lot-srm1453", 1.0),
            ("source-conditioning-protocol-satisfied", 1.0),
            ("unstacked-unsliced-uncontaminated-board", 1.0),
        ]),
        "mean-temperature",
        297.0,
    )
    .with_quantity("bulk-density", QuantitySpec::dimensional(RHO), 41.5)
    .unwrap()
    .with_quantity(
        "temperature-difference",
        QuantitySpec::dimensional(TEMP),
        20.0,
    )
    .unwrap();
    let answer = sample(
        &eps,
        "thermal-conductivity",
        K,
        &p,
        0.00111 - 0.0000424 * 41.5 + 0.000115 * 297.0,
    );
    assert!(matches!(answer.evidence.value.uncertainty,
        UncertaintyModel::RelativeHalfWidth { fraction, confidence } if fraction == 0.015 && confidence == 0.95));
    for bad in [
        p.clone().with("bulk-density", 36.0).unwrap(),
        p.clone().with("mean-temperature", 280.0).unwrap(),
        p.clone().with("temperature-difference", 10.0).unwrap(),
        p.with("source-lot-srm1453", 0.0).unwrap(),
    ] {
        refuses(&eps, "thermal-conductivity", K, &bad);
    }
    let repeated = compile("eps-nist-srm1453-2013", 1);
    assert_eq!(
        eps.to_bytes(),
        repeated.to_bytes(),
        "source surface compilation is deterministic"
    );
}

#[test]
fn g0_g3_xps_keeps_declared_conductivity_thickness_gaps_and_grade_limits() {
    let cases = [
        (
            "30",
            9,
            30.0,
            12e6,
            300e3,
            &[
                (0.03, 0.033),
                (0.04, 0.032),
                (0.05, 0.032),
                (0.06, 0.033),
                (0.07, 0.035),
                (0.12, 0.035),
                (0.14, 0.036),
                (0.16, 0.036),
            ] as &[(f64, f64)],
            [0.035, 0.13, 0.20],
        ),
        (
            "70",
            7,
            39.0,
            25e6,
            700e3,
            &[
                (0.05, 0.033),
                (0.06, 0.033),
                (0.08, 0.035),
                (0.14, 0.035),
                (0.16, 0.036),
            ] as &[(f64, f64)],
            [0.07, 0.15, 0.20],
        ),
    ];
    for (grade, count, rho, modulus, minimum_strength, thicknesses, gaps) in cases {
        let slug = format!("austrotherm-xps-top{grade}-sf-2026");
        let pack = compile(&slug, count);
        let grade_axis = format!("source-grade-{slug}");
        let p = flags(&[
            (&grade_axis, 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
        ]);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "closed-cell-fraction-minimum", Dims::NONE, &p, 0.95);
        let mechanical = p
            .clone()
            .with_quantity(
                "source-test-rate-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        sample(&pack, "compression-modulus", PA, &mechanical, modulus);
        let strength = mechanical
            .with_quantity(
                "source-cs10y-class",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
        sample(
            &pack,
            "compression-strength-minimum",
            PA,
            &strength,
            minimum_strength,
        );

        let declared = p
            .clone()
            .with_quantity(
                "source-is-declared-value",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
        let at_thickness = |metres| {
            declared
                .clone()
                .with_quantity(
                    "nominal-thickness",
                    QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                    metres,
                )
                .unwrap()
        };
        for &(metres, expected) in thicknesses {
            sample(
                &pack,
                "thermal-conductivity",
                K,
                &at_thickness(metres),
                expected,
            );
        }
        for unsupported in gaps {
            refuses(&pack, "thermal-conductivity", K, &at_thickness(unsupported));
        }
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &at_thickness(0.05)
                .with("source-is-declared-value", 0.0)
                .unwrap(),
        );
        refuses(&pack, "density", RHO, &p.with(&grade_axis, 0.0).unwrap());
    }
}
