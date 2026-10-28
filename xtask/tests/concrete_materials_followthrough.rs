//! G0/G3: concrete product, strength role and exposure class survive real compilation.
//! These are declared AAC block properties, not complete concrete constitutive laws.
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
const GRADES: [(&str, f64, f64, f64, [f64; 4]); 4] = [
    ("solar", 460.0, 2.9e6, 2.3e6, [0.10, 0.11, 0.12, 0.17]),
    ("standard", 600.0, 3.6e6, 2.9e6, [0.13, 0.15, 0.16, 0.24]),
    (
        "high-strength",
        730.0,
        7.3e6,
        5.8e6,
        [0.16, 0.18, 0.19, 0.29],
    ),
    (
        "super-strength",
        730.0,
        8.7e6,
        7.0e6,
        [0.16, 0.18, 0.19, 0.29],
    ),
];

fn compile(grade: &str) -> NormalizedPack {
    let slug = format!("hh-celcon-{grade}-bba-2026");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-concrete-products-{}-{nonce}-{grade}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("concrete.fsmatpk");
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
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), 12);
    pack
}

fn flags(point: QueryPoint, pairs: &[(&str, f64)]) -> QueryPoint {
    pairs.iter().fold(point, |p, &(axis, value)| {
        p.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

fn base(grade: &str) -> QueryPoint {
    flags(
        QueryPoint::new(),
        &[
            (&format!("source-grade-hh-celcon-{grade}-bba-2026"), 1.0),
            ("source-cure-age-known", 0.0),
            ("source-specimens-matched", 0.0),
        ],
    )
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

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) -> f64 {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    assert_eq!(answer.evidence.value.dims, dims);
    assert!(matches!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    ));
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer.evidence.value.value
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} unexpectedly admits {name} at {point:?}",
        pack.pack_id()
    );
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g3_four_aac_grades_keep_dry_density_and_strength_requirements_distinct() {
    for (grade, rho, mean, individual, _) in GRADES {
        let pack = compile(grade);
        let p = flags(
            base(grade),
            &[
                ("source-temperature-known", 0.0),
                ("dry-condition", 1.0),
                ("source-moisture-condition-known", 1.0),
            ],
        );
        close(sample(&pack, "density", RHO, &p), rho);
        close(
            sample(&pack, "density-declared-range-lower", RHO, &p),
            rho - 50.0,
        );
        close(
            sample(&pack, "density-declared-range-upper", RHO, &p),
            rho + 50.0,
        );
        refuses(
            &pack,
            "density",
            RHO,
            &flags(p.clone(), &[("dry-condition", 0.0)]),
        );
        refuses(
            &pack,
            "density",
            RHO,
            &flags(p.clone(), &[("source-cure-age-known", 1.0)]),
        );
        refuses(
            &pack,
            "density",
            RHO,
            &flags(
                p,
                &[(&format!("source-grade-hh-celcon-{grade}-bba-2026"), 0.0)],
            ),
        );
        let p = flags(
            base(grade),
            &[
                ("source-temperature-known", 0.0),
                ("source-moisture-condition-known", 0.0),
                ("intended-block-loading-orientation", 1.0),
                ("source-test-rate-known", 0.0),
            ],
        );
        close(
            sample(&pack, "mean-compressive-strength-requirement", PA, &p),
            mean,
        );
        close(
            sample(
                &pack,
                "individual-block-compressive-strength-requirement",
                PA,
                &p,
            ),
            individual,
        );
        refuses(
            &pack,
            "mean-compressive-strength-requirement",
            PA,
            &flags(p.clone(), &[("intended-block-loading-orientation", 0.0)]),
        );
        // A passed specification is not a measured stress-strain curve or wall allowable.
        for property in [
            "compressive-strength",
            "young-modulus",
            "masonry-compressive-strength",
        ] {
            refuses(&pack, property, PA, &p);
        }
    }
}

#[test]
fn g0_g3_dry_thermal_declarations_and_service_exposure_values_are_not_interchangeable() {
    for (grade, _, _, _, ks) in GRADES {
        let pack = compile(grade);
        let dry = flags(
            base(grade),
            &[
                ("source-temperature-known", 1.0),
                ("dry-condition", 1.0),
                ("source-moisture-condition-known", 1.0),
                ("source-thermal-direction-known", 0.0),
                ("source-test-temperature-difference-known", 0.0),
            ],
        )
        .with_quantity(
            "mean-test-temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            283.15,
        )
        .unwrap();
        close(
            sample(&pack, "declared-dry-thermal-conductivity", K, &dry),
            ks[0],
        );
        let other_temperature = dry
            .clone()
            .with_quantity(
                "mean-test-temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                293.15,
            )
            .unwrap();
        refuses(
            &pack,
            "declared-dry-thermal-conductivity",
            K,
            &other_temperature,
        );
        refuses(&pack, "design-thermal-conductivity", K, &dry);
        for exposure in 0..3 {
            let p = flags(
                base(grade),
                &[
                    ("source-temperature-known", 0.0),
                    ("source-moisture-content-known", 0.0),
                    ("source-thermal-direction-known", 0.0),
                    ("blockwork-exposure-class", exposure as f64),
                ],
            );
            close(
                sample(&pack, "design-thermal-conductivity", K, &p),
                ks[exposure + 1],
            );
            refuses(
                &pack,
                "design-thermal-conductivity",
                K,
                &flags(p.clone(), &[("source-moisture-content-known", 1.0)]),
            );
            refuses(
                &pack,
                "design-thermal-conductivity",
                K,
                &flags(p.clone(), &[("blockwork-exposure-class", 0.5)]),
            );
            // Declared/design conductivity never silently becomes a local k(T) law.
            refuses(&pack, "thermal-conductivity", K, &p);
            refuses(&pack, "declared-dry-thermal-conductivity", K, &p);
        }
        assert!(ks[3] > ks[1]);
    }
}

#[test]
fn g0_g3_moisture_movement_and_directional_vapor_factors_preserve_their_roles() {
    for (grade, _, _, _, _) in GRADES {
        let pack = compile(grade);
        let p = flags(
            base(grade),
            &[
                ("source-temperature-known", 0.0),
                ("source-moisture-endpoints-known", 0.0),
            ],
        );
        close(
            sample(&pack, "maximum-declared-moisture-movement", Dims::NONE, &p),
            0.0004,
        );
        refuses(
            &pack,
            "maximum-declared-moisture-movement",
            Dims::NONE,
            &flags(p.clone(), &[("source-moisture-endpoints-known", 1.0)]),
        );
        refuses(&pack, "thermal-expansion", Dims([0, 0, 0, -1, 0, 0]), &p);
        for (direction, mu) in [(1.0, 5.0), (-1.0, 10.0)] {
            let p = flags(
                base(grade),
                &[
                    ("source-temperature-known", 0.0),
                    ("source-moisture-condition-known", 0.0),
                    ("vapor-diffusion-direction", direction),
                ],
            );
            close(
                sample(
                    &pack,
                    "water-vapor-diffusion-resistance-factor",
                    Dims::NONE,
                    &p,
                ),
                mu,
            );
            refuses(
                &pack,
                "water-vapor-diffusion-resistance-factor",
                Dims::NONE,
                &flags(p, &[("vapor-diffusion-direction", 0.0)]),
            );
        }
    }
}
