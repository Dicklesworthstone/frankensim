//! Complete historical NBS/1989-errata reference set for the Elsner equation.
//! Literal NIST rows are checked through the real pack compiler and original
//! reaction owners. No plant recovery or finite-concentration activity is fitted.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, Mass, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        FormationEnthalpySource, ReactionBalanceError, minerals::MineralReaction,
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);

struct Row {
    species: &'static str,
    pack: &'static str,
    phase: &'static str,
    // Independent literal transcription: g/mol, kJ/mol, kJ/mol, J/(mol K).
    values: [f64; 4],
}

// NIST DOI 10.18434/M32124, NBS_Tables_with_ErratumValues.csv.
// Original printed pages 2-164 (Au), 2-89 (CN-), 2-37 (O2), 2-38 (H2O/OH-).
// The gold table was prepared in 1967; this is not the later USBM evaluation.
const ROWS: &[Row] = &[
    Row {
        species: "Au(s)",
        pack: "gold-solid-nbs1989-thermochemistry",
        phase: "source-phase-solid",
        values: [196.9670, 0.0, 0.0, 47.40],
    },
    Row {
        species: "CN[-](aq)",
        pack: "cyanide-ion-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        values: [26.0179, 150.6, 172.4, 94.1],
    },
    Row {
        species: "O2(g)",
        pack: "oxygen-gas-nbs1989-thermochemistry",
        phase: "source-phase-gas",
        values: [31.9988, 0.0, 0.0, 205.138],
    },
    Row {
        species: "H2O(l)",
        pack: "water-liquid-nbs1989-thermochemistry",
        phase: "source-phase-liquid",
        values: [18.0154, -285.830, -237.129, 69.91],
    },
    Row {
        species: "Au(CN)2[-](aq)",
        pack: "dicyanoaurate-ion-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        values: [249.0028, 242.3, 285.8, 172.0],
    },
    Row {
        species: "OH[-](aq)",
        pack: "hydroxide-ion-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        values: [17.0074, -229.994, -157.244, -10.75],
    },
];

struct Source {
    row: &'static Row,
    species: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
}

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

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(phase: &str) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            298.15,
        )
        .unwrap()
        .with_quantity("pressure", semantic(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        (phase, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if phase == "source-phase-aqueous" {
        flags.extend([
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]);
    } else {
        flags.push(("source-pure-species", 1.0));
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}

fn sources() -> &'static [Source] {
    static SOURCES: OnceLock<Vec<Source>> = OnceLock::new();
    SOURCES.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-gold-cyanidation-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        ROWS.iter()
            .map(|row| {
                let target = out.join(format!("{}.fsmatpk", row.pack));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{}/manifest.tsv", row.pack)))
                    .arg("--out")
                    .arg(&target)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{}: {}",
                    row.pack,
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(target).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), row.pack);
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(pack.claims().claims_ordered().count(), 4);
                Source {
                    row,
                    species: SpeciesId::new(row.species).unwrap(),
                    pack,
                    at: point(row.phase),
                }
            })
            .collect()
    })
}

fn heat_sources() -> Vec<FormationEnthalpySource<'static>> {
    sources()
        .iter()
        .map(|source| FormationEnthalpySource {
            species: &source.species,
            pack: &source.pack,
            at: &source.at,
            reference_axis: REFERENCE,
        })
        .collect()
}

fn gibbs_sources() -> Vec<FormationGibbsSource<'static>> {
    sources()
        .iter()
        .map(|source| FormationGibbsSource {
            species: &source.species,
            pack: &source.pack,
            at: &source.at,
            reference_axis: REFERENCE,
        })
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g1_all_24_nbs_values_replay_in_typed_si_with_original_uncertainty() {
    let mut checked = 0;
    for source in sources() {
        let [mass, h, g, s] = source.row.values;
        let entropy = if source.row.phase == "source-phase-aqueous" {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        for (name, dims, expected) in [
            ("molar-mass", MASS, mass * 0.001),
            ("standard-molar-formation-enthalpy", ENERGY, h * 1000.0),
            ("standard-molar-formation-gibbs-energy", ENERGY, g * 1000.0),
            (entropy, ENTROPY, s),
        ] {
            let answer = source
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &source.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected);
            assert_eq!(answer.evidence.value.dims, dims);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            source
                .pack
                .claims()
                .verify_receipt(&answer.receipt)
                .unwrap();
            let claim = source.pack.claims().claim(answer.receipt.selected).unwrap();
            assert_eq!(claim.validity.bound(REFERENCE), Some((1.0, 1.0)));
            checked += 1;
        }
    }
    assert_eq!(checked, 24);
}

#[test]
fn g1_elsner_heat_gibbs_entropy_and_activities_keep_the_nbs_evaluation_distinct() {
    let reaction = MineralReaction::GoldCyanidation.build().unwrap();
    let caloric = reaction.standard_equilibrium(&heat_sources()).unwrap();
    let direct = reaction.standard_gibbs_change(&gibbs_sources()).unwrap();
    close(caloric.enthalpy().molar().value(), -583_916.0);
    close(caloric.molar_entropy().value(), -642.358);
    close(direct.molar().value(), -390_718.0);
    close(caloric.molar_gibbs().value(), -392_396.9623);
    // Preserve the actual discrepancy; it is not explained away as rounding,
    // reconciled by changing a datum, or replaced by the later USBM Gf value.
    close(
        caloric.molar_gibbs().value() - direct.molar().value(),
        -1_678.9623,
    );
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    close(direct.log_equilibrium_constant().value(), 390_718.0 / rt);
    for term in caloric.entropy_terms() {
        let source = sources()
            .iter()
            .find(|s| s.species == term.term.species)
            .unwrap();
        source
            .pack
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
    for term in direct.terms() {
        let source = sources()
            .iter()
            .find(|s| s.species == term.term.species)
            .unwrap();
        source
            .pack
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
    let unit = reaction
        .terms()
        .iter()
        .map(|term| (term.species.clone(), Dimensionless::new(1.0)))
        .collect::<Vec<_>>();
    close(
        direct.at_activities(&unit).unwrap().molar.value(),
        -390_718.0,
    );
    for (species, expected_shift) in [
        ("CN[-](aq)", 8.0 * rt * 10.0_f64.ln()),
        ("OH[-](aq)", -4.0 * rt * 10.0_f64.ln()),
    ] {
        let mut activities = unit.clone();
        activities
            .iter_mut()
            .find(|(id, _)| id.as_str() == species)
            .unwrap()
            .1 = Dimensionless::new(0.1);
        close(
            direct.at_activities(&activities).unwrap().molar.value(),
            -390_718.0 + expected_shift,
        );
    }
}

#[test]
fn g1_finite_cyanide_and_oxygen_stocks_bound_gold_dissolution_and_heat() {
    let reaction = MineralReaction::GoldCyanidation.build().unwrap();
    let heat = reaction.standard_enthalpy_change(&heat_sources()).unwrap();
    // Four reactants in source order; products are explicitly absent initially.
    for (stocks, extent, expected) in [
        ([4.0, 4.0, 1.0, 2.0], 0.5, [2.0, 0.0, 0.5, 1.0, 2.0, 2.0]),
        ([4.0, 8.0, 0.25, 2.0], 0.25, [3.0, 6.0, 0.0, 1.5, 1.0, 1.0]),
    ] {
        let initial = sources()
            .iter()
            .enumerate()
            .map(|(i, source)| {
                (
                    source.species.clone(),
                    Amount::new(stocks.get(i).copied().unwrap_or(0.0)),
                )
            })
            .collect::<Vec<_>>();
        let bound = reaction.limiting_extent(&initial[..4]).unwrap();
        close(bound.value(), extent);
        let final_inventory = reaction.advance_inventory(&initial, bound).unwrap();
        for (source, expected) in sources().iter().zip(expected) {
            let actual = final_inventory
                .iter()
                .find(|(id, _)| *id == source.species)
                .unwrap()
                .1;
            close(actual.value(), expected);
        }
        close(heat.at_extent(bound).unwrap().value(), -583_916.0 * extent);
        assert!(matches!(
            reaction.advance_inventory(&initial, Amount::new(extent + 0.125)),
            Err(ReactionBalanceError::ExtentExceedsInventory { .. })
        ));
    }
    let gold = &sources()[0];
    let mass = gold
        .pack
        .claims()
        .query_typed(
            &PropertyKey::new("molar-mass", MASS),
            &gold.at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    // One mole of actually reacted gold is one quarter of the written equation.
    let extent = reaction
        .extent_from_reacted_mass(
            &gold.species,
            Mass::new(0.196967),
            MolarMass::new(mass.evidence.value.value),
        )
        .unwrap();
    close(extent.value(), 0.25);
    close(heat.at_extent(extent).unwrap().value(), -145_979.0);
}

#[test]
fn g0_g3_reference_phase_temperature_and_ionic_entropy_substitution_refuse() {
    for source in sources() {
        let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
        let changed = [
            QueryPoint::new(),
            source
                .at
                .clone()
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::AbsoluteTemperature),
                    300.0,
                )
                .unwrap(),
            source
                .at
                .clone()
                .with_quantity("pressure", semantic(QuantityKind::Pressure), 101_325.0)
                .unwrap(),
            source
                .at
                .clone()
                .with(REFERENCE, 0.0)
                .unwrap()
                .with("source-elemental-reference-usgs-b2131", 1.0)
                .unwrap(),
            source.at.clone().with(source.row.phase, 0.0).unwrap(),
            source
                .at
                .clone()
                .with("source-specimen-qualified", 1.0)
                .unwrap(),
        ];
        for at in changed {
            assert!(
                source
                    .pack
                    .claims()
                    .query_typed(&key, &at, SelectionPolicy::SingleClaimOnly,)
                    .is_err(),
                "unexpected admission for {}",
                source.row.pack
            );
        }
        if source.row.phase == "source-phase-aqueous" {
            assert!(
                source
                    .pack
                    .claims()
                    .query_typed(
                        &PropertyKey::new("standard-molar-entropy", ENTROPY),
                        &source.at,
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .is_err()
            );
            let no_ionic_reference = source
                .at
                .clone()
                .with("source-ionic-reference-hydrogen-zero", 0.0)
                .unwrap();
            assert!(
                source
                    .pack
                    .claims()
                    .query_typed(&key, &no_ionic_reference, SelectionPolicy::SingleClaimOnly,)
                    .is_err()
            );
        }
    }
}
