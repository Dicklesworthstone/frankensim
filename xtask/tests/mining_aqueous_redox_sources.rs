//! Source-coherent reference data for sulfide oxidation and acid neutralization.
//! Primary numbers are transcribed from original USGS B2131 page images.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{FormationEnthalpySource, minerals::MineralReaction},
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const HF: &str = "standard-molar-formation-enthalpy";
const GF: &str = "standard-molar-formation-gibbs-energy";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const PYRITE: &str = "pyrite-solid-usgs-b2131-thermochemistry";
const OXYGEN: &str = "oxygen-gas-usgs-b2131";
const WATER: &str = "water-liquid-usgs-b2131";
const FERROUS: &str = "iron-ii-ion-aqueous-usgs-b2131";
const FERRIC: &str = "iron-iii-ion-aqueous-usgs-b2131";
const SULFATE: &str = "sulfate-ion-aqueous-usgs-b2131";
const HYDROGEN: &str = "hydrogen-ion-aqueous-usgs-b2131";

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(slug: &str) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-aqueous-redox-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("material.fsmatpk");
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
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
    pack
}

fn point(phase: &str) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![
        ("source-standard-state", 1.0),
        (REFERENCE, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if phase == "aqueous" {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
            ("source-individual-ion", 1.0),
        ]);
    } else {
        flags.push(("source-pure-species", 1.0));
        flags.push((
            match phase {
                "gas" => "source-phase-gas",
                "liquid" => "source-phase-liquid",
                "pyrite" | "calcite" => "source-phase-solid",
                _ => panic!("unknown test phase"),
            },
            1.0,
        ));
        if phase == "pyrite" {
            flags.push(("source-phase-pyrite", 1.0));
        } else if phase == "calcite" {
            flags.push(("source-phase-calcite", 1.0));
        }
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn assert_value(pack: &NormalizedPack, at: &QueryPoint, key: PropertyKey, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&key, at, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert!(
        (answer.evidence.value.value - expected).abs() <= expected.abs().max(1.0) * 1e-12,
        "{} {} != {expected}",
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
    let claim = pack.claims().claim(answer.receipt.selected).unwrap();
    for (axis, value) in at.axes() {
        assert_eq!(claim.validity.bound(axis), Some((*value, *value)));
        assert_eq!(
            claim.validity.axis_quantities().get(axis),
            at.axis_quantities().get(axis)
        );
    }
    for wrong in [
        QueryPoint::new(),
        at.clone().with("temperature", 300.0).unwrap(),
        at.clone().with("pressure", 101_325.0).unwrap(),
        at.clone().with(REFERENCE, 0.0).unwrap(),
    ] {
        assert!(
            pack.claims()
                .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
    }
}

#[test]
fn g0_g3_redox_profiles_preserve_printed_quantities_and_reference_states() {
    // USGS B2131 p7 (Fe2+/Fe3+), p18 (SO4), p8 (O2), p12 (pyrite), p16 (water).
    // Gibbs and enthalpy columns remain distinct, including sulfate's -744.0.
    for (slug, phase, mass, hf, gf, entropy) in [
        (FERROUS, "aqueous", 0.055847, -91_100.0, -90_000.0, -107.1),
        (FERRIC, "aqueous", 0.055847, -49_900.0, -16_700.0, -280.0),
        (SULFATE, "aqueous", 0.096064, -909_300.0, -744_000.0, 18.5),
        (OXYGEN, "gas", 0.031999, 0.0, 0.0, 205.15),
        (PYRITE, "pyrite", 0.119979, -171_500.0, -160_100.0, 52.9),
        (WATER, "liquid", 0.018015, -285_800.0, -237_100.0, 70.0),
    ] {
        let pack = compile(slug);
        let count = match slug {
            FERRIC => 4,
            OXYGEN => 5,
            _ => 6,
        };
        assert_eq!(pack.claims().claim_count(), count);
        let at = point(phase);
        let entropy_key = if phase == "aqueous" {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        for (name, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, mass),
            (HF, MOLAR_ENERGY, hf),
            (GF, MOLAR_ENERGY, gf),
            (entropy_key, MOLAR_ENTROPY, entropy),
        ] {
            assert_value(&pack, &at, PropertyKey::new(name, dims), expected);
        }
        if phase == "aqueous" {
            assert!(
                pack.claims()
                    .query_typed(
                        &PropertyKey::new("standard-molar-entropy", MOLAR_ENTROPY),
                        &at,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
}

#[test]
fn g1_existing_carbonate_profiles_gain_separate_gibbs_formation_claims() {
    // Original source rows: p268 (calcite), p7 (H+), p6 (Ca2+), p15 (HCO3-).
    for (slug, phase, expected) in [
        (
            "calcite-solid-usgs-b2131-thermochemistry",
            "calcite",
            -1_128_500.0,
        ),
        (HYDROGEN, "aqueous", 0.0),
        ("calcium-ion-aqueous-usgs-b2131", "aqueous", -553_600.0),
        ("bicarbonate-ion-aqueous-usgs-b2131", "aqueous", -586_800.0),
    ] {
        assert_value(
            &compile(slug),
            &point(phase),
            PropertyKey::new(GF, MOLAR_ENERGY),
            expected,
        );
    }
}

#[test]
fn g1_real_sources_close_pyrite_oxidation_heat_with_limited_oxygen() {
    // 2 FeS2(s) + 7 O2(g) + 2 H2O(l) -> 2 Fe2+ + 4 SO4(2-) + 4 H+.
    // O2 is gas-referenced; this does not assume a dissolved-oxygen activity.
    let packs = [PYRITE, OXYGEN, WATER, FERROUS, SULFATE, HYDROGEN].map(compile);
    let ids = [
        "FeS2(s)",
        "O2(g)",
        "H2O(l)",
        "Fe[2+](aq)",
        "SO4[2-](aq)",
        "H+(aq)",
    ]
    .map(|id| SpeciesId::new(id).unwrap());
    let points = ["pyrite", "gas", "liquid", "aqueous", "aqueous", "aqueous"].map(point);
    let sources: [_; 6] = std::array::from_fn(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let reaction = MineralReaction::PyriteOxidationFerrous.build().unwrap();
    let heat = reaction.standard_enthalpy_change(&sources).unwrap();
    // Independently summed printed source entries, per integer equation.
    assert_eq!(heat.molar().value(), -2_904_800.0);
    let initial: [_; 6] = std::array::from_fn(|i| {
        (
            ids[i].clone(),
            Amount::new([10.0, 7.0, 10.0, 0.0, 0.0, 0.0][i]),
        )
    });
    let extent = reaction.limiting_extent(&initial[..3]).unwrap();
    assert_eq!(extent.value(), 1.0);
    let after = reaction.advance_inventory(&initial, extent).unwrap();
    for (i, expected) in [8.0, 0.0, 8.0, 2.0, 4.0, 4.0].into_iter().enumerate() {
        assert_eq!(
            after
                .iter()
                .find(|(id, _)| id == &ids[i])
                .unwrap()
                .1
                .value(),
            expected
        );
    }
    assert_eq!(heat.at_extent(extent).unwrap().value(), -2_904_800.0);
    assert!(
        reaction
            .advance_inventory(&initial, Amount::new(1.1))
            .is_err()
    );
    assert!(reaction.standard_enthalpy_change(&sources[..5]).is_err());
    for term in heat.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
    }
}
