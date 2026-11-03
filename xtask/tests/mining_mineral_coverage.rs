//! Additional ore phases beyond the original mining mineral-density tranche.
//! Checks source transport, not ore characterization or process calibration.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{Dims, QuantitySpec};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);

struct Mineral {
    name: &'static str,
    // Independent transcription of D(meas.), in original g/cm3. A source
    // range or parenthetical magnitude never becomes a confidence interval.
    density: &'static [(&'static str, f64)],
    mohs: [f64; 2],
}

const MINERALS: &[Mineral] = &[
    Mineral {
        name: "bornite",
        density: &[
            ("density-reported-range-minimum", 5.06),
            ("density-reported-range-maximum", 5.08),
        ],
        mohs: [3.0, 3.25],
    },
    Mineral {
        name: "sphalerite",
        density: &[
            ("density-reported-range-minimum", 3.9),
            ("density-reported-range-maximum", 4.1),
        ],
        mohs: [3.5, 4.0],
    },
    Mineral {
        name: "galena",
        density: &[("density", 7.58)],
        mohs: [2.5, 2.75],
    },
    Mineral {
        name: "arsenopyrite",
        density: &[
            ("density", 6.07),
            ("density-source-parenthetical-spread", 0.15),
        ],
        mohs: [5.5, 6.0],
    },
    Mineral {
        name: "malachite",
        density: &[
            ("density", 4.05),
            ("density-source-parenthetical-spread", 0.02),
        ],
        mohs: [3.5, 4.0],
    },
    Mineral {
        name: "stannite",
        density: &[
            ("density-reported-range-minimum", 4.3),
            ("density-reported-range-maximum", 4.5),
        ],
        mohs: [4.0, 4.0],
    },
    Mineral {
        name: "cuprite",
        density: &[("density", 6.14)],
        mohs: [3.5, 4.0],
    },
    Mineral {
        name: "tenorite",
        density: &[("density", 6.45)],
        mohs: [3.5, 3.5],
    },
    Mineral {
        name: "kaolinite",
        density: &[
            ("density-reported-range-minimum", 2.61),
            ("density-reported-range-maximum", 2.68),
        ],
        mohs: [2.0, 2.5],
    },
];

fn workspace() -> PathBuf {
    std::env::var_os("CARGO_WORKSPACE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .into()
        })
}

fn packs() -> &'static [NormalizedPack] {
    static PACKS: OnceLock<Vec<NormalizedPack>> = OnceLock::new();
    PACKS.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let output_dir = std::env::temp_dir().join(format!(
            "frankensim-mining-mineral-coverage-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&output_dir).unwrap();
        MINERALS
            .iter()
            .map(|mineral| {
                let slug = format!("{}-handbook-mineralogy", mineral.name);
                let output = output_dir.join(format!("{slug}.fsmatpk"));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
                    .arg("--out")
                    .arg(&output)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{slug}: {}",
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(output).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(pack.pack_id(), slug);
                assert_eq!(
                    pack.claims().claims_ordered().count(),
                    mineral.density.len() + 2
                );
                pack
            })
            .collect()
    })
}

fn context(mineral: &Mineral) -> QueryPoint {
    let point = base_context(mineral);
    if mineral.name == "tenorite" {
        point
            .with_quantity(
                "source-density-synthetic",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap()
    } else {
        point
    }
}

fn base_context(mineral: &Mineral) -> QueryPoint {
    [
        (format!("source-mineral-{}", mineral.name), 1.0),
        ("mineral-reference-model".into(), 1.0),
        ("source-test-temperature-known".into(), 0.0),
        ("source-test-pressure-known".into(), 0.0),
        ("bulk-ore-qualified".into(), 0.0),
        ("source-specimen-qualified".into(), 0.0),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |point, (axis, value)| {
        point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

#[test]
fn additional_ore_phase_values_survive_compilation_and_typed_queries() {
    let mut claims = 0;
    for (mineral, pack) in MINERALS.iter().zip(packs()) {
        let mut expected: Vec<(&str, Dims, f64)> = mineral
            .density
            .iter()
            .map(|&(property, value)| (property, RHO, value * 1000.0))
            .collect();
        expected.extend([
            ("hardness-mohs-source-minimum", Dims::NONE, mineral.mohs[0]),
            ("hardness-mohs-source-maximum", Dims::NONE, mineral.mohs[1]),
        ]);
        for (property, dims, value) in expected {
            let answer = pack
                .claims()
                .query_typed(
                    &PropertyKey::new(property, dims),
                    &context(mineral),
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            assert!(
                (answer.evidence.value.value - value).abs() <= 1e-12 * value.max(1.0),
                "{} {property}: {} != {value}",
                mineral.name,
                answer.evidence.value.value
            );
            assert_eq!(answer.evidence.value.dims, dims);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            claims += 1;
        }
    }
    assert_eq!(claims, 33);
}

#[test]
fn tenorite_density_preserves_its_synthetic_source_qualification() {
    let index = MINERALS
        .iter()
        .position(|mineral| mineral.name == "tenorite")
        .unwrap();
    let mineral = &MINERALS[index];
    let pack = &packs()[index];
    let density = PropertyKey::new("density", RHO);
    for unsupported in [
        base_context(mineral),
        context(mineral)
            .with("source-density-synthetic", 0.0)
            .unwrap(),
    ] {
        assert!(
            pack.claims()
                .query_typed(&density, &unsupported, SelectionPolicy::SingleClaimOnly)
                .is_err(),
            "synthetic measurement must not silently become natural-ore density"
        );
    }

    // The source does not identify the hardness specimen as synthetic.
    let hardness = pack
        .claims()
        .query_typed(
            &PropertyKey::new("hardness-mohs-source-minimum", Dims::NONE),
            &base_context(mineral),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(hardness.evidence.value.value, 3.5);
    let claim = pack.claims().claim(hardness.receipt.selected).unwrap();
    assert_eq!(claim.validity.bound("source-density-synthetic"), None);
}

#[test]
fn mineral_ranges_and_ordinal_hardness_do_not_invent_process_properties() {
    for (mineral, pack) in MINERALS.iter().zip(packs()) {
        let property = PropertyKey::new(mineral.density[0].0, RHO);
        let valid = context(mineral);
        assert!(
            pack.claims()
                .query_typed(
                    &property,
                    &QueryPoint::new(),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        for (axis, value) in [
            (format!("source-mineral-{}", mineral.name), 0.0),
            ("mineral-reference-model".into(), 0.0),
            ("source-test-temperature-known".into(), 1.0),
            ("source-test-pressure-known".into(), 1.0),
            ("bulk-ore-qualified".into(), 1.0),
            ("source-specimen-qualified".into(), 1.0),
        ] {
            let point = valid
                .clone()
                .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
                .unwrap();
            assert!(
                pack.claims()
                    .query_typed(&property, &point, SelectionPolicy::SingleClaimOnly)
                    .is_err(),
                "{} {point:?}",
                mineral.name
            );
        }
        for (name, dims) in [
            ("bulk-ore-density", RHO),
            ("bond-ball-mill-work-index", Dims([2, 0, -2, 0, 0, 0])),
            ("metallurgical-recovery", Dims::NONE),
            ("hardness-mohs-source-minimum", RHO),
        ] {
            assert!(
                pack.claims()
                    .query_typed(
                        &PropertyKey::new(name, dims),
                        &valid,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err(),
                "{} cannot supply {name}",
                mineral.name
            );
        }
        if mineral.density[0].0 == "density-reported-range-minimum" {
            assert!(
                pack.claims()
                    .query_typed(
                        &PropertyKey::new("density", RHO),
                        &valid,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err(),
                "range must not invent a point density"
            );
        }
    }
}
