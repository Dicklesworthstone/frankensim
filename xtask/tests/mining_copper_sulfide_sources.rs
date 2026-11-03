//! Copper-sulfide and reagent reference data from USGS Bulletin 2131.
//! These verify source transport and applicability, not measured leach recovery.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const USGS: &str = "source-elemental-reference-usgs-b2131";

struct Source {
    pack: &'static str,
    solid_phase: Option<&'static str>,
    // Independent transcription of the printed row, in original units:
    // g/mol, kJ/mol, kJ/mol and J/(mol K).
    values: [f64; 4],
}

// Sulfides: printed p11 / scan m1/17; Cu2+: p6 / m1/12; sulfur: p8 / m1/14.
// https://digital.library.unt.edu/ark:/67531/metadc1114954/
const SOURCES: &[Source] = &[
    Source {
        pack: "chalcopyrite-solid-usgs-b2131-thermochemistry",
        solid_phase: Some("source-phase-chalcopyrite"),
        values: [183.525, -194.9, -195.1, 124.9],
    },
    Source {
        pack: "chalcocite-solid-usgs-b2131-thermochemistry",
        solid_phase: Some("source-phase-chalcocite"),
        values: [159.158, -83.9, -89.2, 116.2],
    },
    Source {
        pack: "covellite-solid-usgs-b2131-thermochemistry",
        solid_phase: Some("source-phase-covellite"),
        values: [95.612, -54.6, -55.3, 67.4],
    },
    Source {
        pack: "bornite-solid-usgs-b2131-thermochemistry",
        solid_phase: Some("source-phase-bornite"),
        values: [501.841, -371.6, -396.7, 398.5],
    },
    Source {
        pack: "copper-ii-ion-aqueous-usgs-b2131",
        solid_phase: None,
        values: [63.546, 64.9, 65.1, -98.0],
    },
    Source {
        pack: "sulfur-orthorhombic-solid-usgs-b2131",
        solid_phase: Some("source-phase-sulfur-orthorhombic"),
        values: [32.066, 0.0, 0.0, 32.05],
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
        let out = std::env::temp_dir().join(format!(
            "frankensim-copper-sulfide-sources-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        SOURCES
            .iter()
            .map(|source| {
                let target = out.join(format!("{}.fsmatpk", source.pack));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(
                        workspace()
                            .join(format!("data/matdb/seed-v1/{}/manifest.tsv", source.pack)),
                    )
                    .arg("--out")
                    .arg(&target)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{}: {}",
                    source.pack,
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(target).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), source.pack);
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(
                    pack.claims().claims_ordered().count(),
                    if source.solid_phase.is_none() { 6 } else { 4 }
                );
                pack
            })
            .collect()
    })
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn context(source: &Source) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            298.15,
        )
        .unwrap()
        .with_quantity("pressure", semantic(QuantityKind::Pressure), 100000.0)
        .unwrap();
    let mut axes = vec![
        (USGS, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if let Some(phase) = source.solid_phase {
        axes.extend([
            ("source-phase-solid", 1.0),
            (phase, 1.0),
            ("source-pure-species", 1.0),
        ]);
    } else {
        axes.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]);
    }
    for (axis, value) in axes {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn properties(source: &Source) -> [(&'static str, Dims, f64); 4] {
    [
        ("molar-mass", MASS, source.values[0] * 0.001),
        (
            "standard-molar-formation-enthalpy",
            ENERGY,
            source.values[1] * 1000.0,
        ),
        (
            "standard-molar-formation-gibbs-energy",
            ENERGY,
            source.values[2] * 1000.0,
        ),
        (
            if source.solid_phase.is_some() {
                "standard-molar-entropy"
            } else {
                "conventional-standard-molar-entropy"
            },
            ENTROPY,
            source.values[3],
        ),
    ]
}

#[test]
fn copper_sulfide_reference_rows_compile_and_replay_in_typed_si() {
    let mut checked = 0;
    for (source, pack) in SOURCES.iter().zip(packs()) {
        for (name, dims, expected) in properties(source) {
            let answer = pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &context(source),
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            assert!(
                (answer.evidence.value.value - expected).abs() <= 2e-14 * expected.abs().max(1.0),
                "{} {name}",
                source.pack
            );
            assert_eq!(answer.evidence.value.dims, dims);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            let claim = pack.claims().claim(answer.receipt.selected).unwrap();
            assert_eq!(claim.validity.bound(USGS), Some((1.0, 1.0)));
            checked += 1;
        }
    }
    assert_eq!(checked, 24);
}

#[test]
fn sources_refuse_temperature_pressure_phase_and_reference_substitution() {
    for (source, pack) in SOURCES.iter().zip(packs()) {
        let phase_axis = source
            .solid_phase
            .unwrap_or("source-ionic-reference-hydrogen-zero");
        let points = [
            QueryPoint::new(),
            context(source)
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::AbsoluteTemperature),
                    300.0,
                )
                .unwrap(),
            context(source)
                .with_quantity("pressure", semantic(QuantityKind::Pressure), 101325.0)
                .unwrap(),
            context(source)
                .with(USGS, 0.0)
                .unwrap()
                .with("source-elemental-reference-janaf", 1.0)
                .unwrap(),
            context(source).with(phase_axis, 0.0).unwrap(),
            context(source)
                .with("source-specimen-qualified", 1.0)
                .unwrap(),
        ];
        for point in points {
            for (name, dims, _) in properties(source) {
                assert!(
                    pack.claims()
                        .query_typed(
                            &PropertyKey::new(name, dims),
                            &point,
                            SelectionPolicy::SingleClaimOnly
                        )
                        .is_err(),
                    "{} {name}",
                    source.pack
                );
            }
        }
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new("specific-heat-capacity", Dims([2, 0, -2, -1, 0, 0])),
                    &context(source),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn copper_ionic_entropy_and_sulfur_atom_basis_remain_distinct() {
    let copper = &SOURCES[4];
    let pack = &packs()[4];
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::new("standard-molar-entropy", ENTROPY),
                &context(copper),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    let conventional = pack
        .claims()
        .query_typed(
            &PropertyKey::new("conventional-standard-molar-entropy", ENTROPY),
            &context(copper),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(conventional.evidence.value.value, -98.0);
    assert_eq!(SOURCES[5].values[0], 32.066);
    let sulfur = packs()[5]
        .claims()
        .query_typed(
            &PropertyKey::new("molar-mass", MASS),
            &context(&SOURCES[5]),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert!((sulfur.evidence.value.value - 0.032066).abs() < 1e-16);
    assert!(
        packs()[5]
            .claims()
            .query_typed(
                &PropertyKey::new("molar-mass", MASS),
                &context(&SOURCES[5])
                    .with("source-phase-sulfur-orthorhombic", 0.0)
                    .unwrap()
                    .with("source-phase-sulfur-monoclinic", 1.0)
                    .unwrap(),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}
