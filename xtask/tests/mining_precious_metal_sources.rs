//! Actual NBS packs close native-silver leaching and finite shared-zinc recovery.
//! Source facts, equation amounts and stage limits remain distinct from kinetics.
#![deny(unsafe_code)]

use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const ROWS: &[(&str, &str, &str, bool)] = &[
    (
        "Ag(s)",
        "silver-solid-nbs1989-thermochemistry",
        "source-phase-solid",
        false,
    ),
    (
        "Ag(CN)2[-](aq)",
        "dicyanoargentate-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "Zn(s)",
        "zinc-solid-nbs1989-thermochemistry",
        "source-phase-solid",
        false,
    ),
    (
        "Zn(CN)4[2-](aq)",
        "tetracyanozincate-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "Au(s)",
        "gold-solid-nbs1989-thermochemistry",
        "source-phase-solid",
        false,
    ),
    (
        "Au(CN)2[-](aq)",
        "dicyanoaurate-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "CN[-](aq)",
        "cyanide-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "OH[-](aq)",
        "hydroxide-ion-aqueous-nbs1989",
        "source-phase-aqueous",
        true,
    ),
    (
        "O2(aq)",
        "oxygen-aqueous-nbs1989-thermochemistry",
        "source-phase-aqueous",
        false,
    ),
    (
        "O2(g)",
        "oxygen-gas-nbs1989-thermochemistry",
        "source-phase-gas",
        false,
    ),
    (
        "H2O(l)",
        "water-liquid-nbs1989-thermochemistry",
        "source-phase-liquid",
        false,
    ),
];

struct Source {
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

fn point(phase: &str, ion: bool) -> QueryPoint {
    let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
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
        flags.push(("source-aqueous-unit-activity-infinite-dilution", 1.0));
        if ion {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
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
            "frankensim-precious-metals-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        ROWS.iter()
            .map(|(species, slug, phase, ion)| {
                let target = out.join(format!("{slug}.fsmatpk"));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
                    .arg("--out")
                    .arg(&target)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{slug}: {}",
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(target).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), *slug);
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(pack.claims().claim_count(), 4);
                Source {
                    species: SpeciesId::new(*species).unwrap(),
                    pack,
                    at: point(phase, *ion),
                }
            })
            .collect()
    })
}

fn source(species: &SpeciesId) -> &'static Source {
    sources()
        .iter()
        .find(|source| source.species == *species)
        .unwrap()
}
fn heat_sources(
    species: impl IntoIterator<Item = SpeciesId>,
) -> Vec<FormationEnthalpySource<'static>> {
    species
        .into_iter()
        .map(|species| {
            let s = source(&species);
            FormationEnthalpySource {
                species: &s.species,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn gibbs_sources(reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'static>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let s = source(&term.species);
            FormationGibbsSource {
                species: &s.species,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn inventory(network: &ReactionNetwork, stocks: &[(&str, f64)]) -> Vec<(SpeciesId, Amount)> {
    network
        .species()
        .iter()
        .map(|species| {
            let value = stocks
                .iter()
                .find(|(name, _)| *name == species.as_str())
                .map_or(0.0, |(_, n)| *n);
            (species.clone(), Amount::new(value))
        })
        .collect()
}
fn amount(stocks: &[(SpeciesId, Amount)], species: &str) -> f64 {
    stocks
        .iter()
        .find(|(name, _)| name.as_str() == species)
        .unwrap()
        .1
        .value()
}
fn mass(stocks: &[(SpeciesId, Amount)]) -> f64 {
    stocks
        .iter()
        .map(|(species, n)| {
            let s = source(species);
            let m = s
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new("molar-mass", MolarMass::DIMS),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            n.value() * m.evidence.value.value
        })
        .sum()
}

#[test]
fn literal_silver_and_zinc_rows_retain_all_sixteen_typed_facts_and_receipts() {
    // Independent literal NBS Tables33:Zn and37:Ag: g/mol, kJ/mol, kJ/mol, J/(mol K).
    for (s, values) in sources()[..4].iter().zip([
        [107.8700, 0.0, 0.0, 42.55],
        [159.9058, 270.3, 305.5, 192.0],
        [65.3700, 0.0, 0.0, 41.63],
        [169.4416, 342.3, 446.9, 226.0],
    ]) {
        let entropy = if s.species.as_str().ends_with("(aq)") {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        for (name, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, values[0] * 0.001),
            (
                "standard-molar-formation-enthalpy",
                ENERGY,
                values[1] * 1000.0,
            ),
            (
                "standard-molar-formation-gibbs-energy",
                ENERGY,
                values[2] * 1000.0,
            ),
            (entropy, ENTROPY, values[3]),
        ] {
            let answer = s
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            s.pack.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
}

#[test]
fn actual_sources_close_four_reaction_heats_and_keep_direct_gibbs_independent() {
    for (recipe, h, g, entropy) in [
        (
            MineralReaction::SilverCyanidation,
            -471_916.0,
            -311_918.0,
            -542.958,
        ),
        (
            MineralReaction::SilverCyanidationDissolvedOxygen,
            -460_216.0,
            -328_318.0,
            -448.720,
        ),
        (
            MineralReaction::GoldCyanideCementationWithZinc,
            -142_300.0,
            -124_700.0,
            -64.83,
        ),
        (
            MineralReaction::SilverCyanideCementationWithZinc,
            -198_300.0,
            -164_100.0,
            -114.53,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let heats = heat_sources(reaction.terms().iter().map(|term| term.species.clone()));
        let caloric = reaction.standard_equilibrium(&heats).unwrap();
        let direct = reaction
            .standard_gibbs_change(&gibbs_sources(&reaction))
            .unwrap();
        close(caloric.enthalpy().molar().value(), h);
        close(caloric.molar_entropy().value(), entropy);
        close(caloric.molar_gibbs().value(), h - 298.15 * entropy);
        close(direct.molar().value(), g);
        // The independent source evaluations are not reconciled by editing data.
        assert!((caloric.molar_gibbs().value() - direct.molar().value()).abs() > 1.0);
        for term in caloric.enthalpy().terms() {
            source(&term.term.species)
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        for term in direct.terms() {
            source(&term.term.species)
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
    }
}

#[test]
fn mixed_gold_silver_cementation_shares_finite_zinc_and_conserves_source_mass() {
    let gold = MineralReaction::GoldCyanideCementationWithZinc
        .build()
        .unwrap();
    let silver = MineralReaction::SilverCyanideCementationWithZinc
        .build()
        .unwrap();
    let network = ReactionNetwork::new(vec![silver.clone(), gold.clone()]).unwrap();
    let initial = inventory(
        &network,
        &[
            ("Au(CN)2[-](aq)", 10.0),
            ("Ag(CN)2[-](aq)", 20.0),
            ("Zn(s)", 8.0),
        ],
    );
    let requested = vec![
        (gold.id().clone(), Amount::new(3.0)),
        (silver.id().clone(), Amount::new(5.0)),
    ];
    let h = heat_sources(network.species().iter().cloned());
    let result = network
        .advance_with_standard_enthalpy(&initial, &requested, &h)
        .unwrap();
    for (species, expected) in [
        ("Au(s)", 6.0),
        ("Ag(s)", 10.0),
        ("Zn(s)", 0.0),
        ("Au(CN)2[-](aq)", 4.0),
        ("Ag(CN)2[-](aq)", 10.0),
        ("Zn(CN)4[2-](aq)", 8.0),
    ] {
        close(amount(&result.inventory.final_inventory, species), expected);
    }
    close(result.enthalpy.value(), -1_418_400.0);
    close(mass(&initial), mass(&result.inventory.final_inventory));
    close(amount(&result.inventory.consumed, "Zn(s)"), 8.0);
    let excessive = vec![
        (gold.id().clone(), Amount::new(3.0)),
        (silver.id().clone(), Amount::new(6.0)),
    ];
    assert!(
        matches!(network.advance_with_standard_enthalpy(&initial, &excessive, &h),
        Err(ReactionNetworkError::InsufficientInventory { species, available, required })
        if species.as_str() == "Zn(s)" && available.value() == 8.0 && required.value() == 9.0)
    );
    // The rejected stage leaves the caller's complete original stocks intact.
    assert_eq!(amount(&initial, "Zn(s)"), 8.0);
    let mut reversed = initial.clone();
    reversed.reverse();
    let mut requests = requested.clone();
    requests.reverse();
    let replay = network
        .advance_with_standard_enthalpy(&reversed, &requests, &h)
        .unwrap();
    assert_eq!(
        result.inventory.final_inventory,
        replay.inventory.final_inventory
    );
    assert_eq!(
        result.enthalpy.value().to_bits(),
        replay.enthalpy.value().to_bits()
    );
}

#[test]
fn silver_leaching_products_fund_only_a_subsequent_explicit_cementation_stage() {
    let leach = MineralReaction::SilverCyanidationDissolvedOxygen
        .build()
        .unwrap();
    let cement = MineralReaction::SilverCyanideCementationWithZinc
        .build()
        .unwrap();
    let network = ReactionNetwork::new(vec![leach.clone(), cement.clone()]).unwrap();
    let initial = inventory(
        &network,
        &[
            ("Ag(s)", 4.0),
            ("CN[-](aq)", 8.0),
            ("O2(aq)", 1.0),
            ("H2O(l)", 2.0),
            ("Zn(s)", 2.0),
        ],
    );
    let h = heat_sources(network.species().iter().cloned());
    let both = [
        (leach.id().clone(), Amount::new(1.0)),
        (cement.id().clone(), Amount::new(2.0)),
    ];
    assert!(
        matches!(network.advance_with_standard_enthalpy(&initial, &both, &h),
        Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "Ag(CN)2[-](aq)")
    );
    let leached = network
        .advance_with_standard_enthalpy(
            &initial,
            &[
                (leach.id().clone(), Amount::new(1.0)),
                (cement.id().clone(), Amount::new(0.0)),
            ],
            &h,
        )
        .unwrap();
    close(amount(&leached.inventory.final_inventory, "O2(aq)"), 0.0);
    close(
        amount(&leached.inventory.final_inventory, "Ag(CN)2[-](aq)"),
        4.0,
    );
    let recovered = network
        .advance_with_standard_enthalpy(
            &leached.inventory.final_inventory,
            &[
                (leach.id().clone(), Amount::new(0.0)),
                (cement.id().clone(), Amount::new(2.0)),
            ],
            &h,
        )
        .unwrap();
    close(amount(&recovered.inventory.final_inventory, "Ag(s)"), 4.0);
    close(
        amount(&recovered.inventory.final_inventory, "Zn(CN)4[2-](aq)"),
        2.0,
    );
    close(
        amount(&recovered.inventory.final_inventory, "OH[-](aq)"),
        4.0,
    );
    close(
        amount(&recovered.inventory.final_inventory, "CN[-](aq)"),
        0.0,
    );
    close(mass(&initial), mass(&recovered.inventory.final_inventory));
    close(
        leached.enthalpy.value() + recovered.enthalpy.value(),
        -856_816.0,
    );
}

#[test]
fn point_phase_and_ionic_reference_cannot_be_relabelled_into_source_authority() {
    for (index, s) in sources()[..4].iter().enumerate() {
        let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
        for at in [
            QueryPoint::new(),
            s.at.clone()
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    300.0,
                )
                .unwrap(),
            s.at.clone()
                .with(REFERENCE, 0.0)
                .unwrap()
                .with("source-elemental-reference-usgs-b2131", 1.0)
                .unwrap(),
            s.at.clone().with(ROWS[index].2, 0.0).unwrap(),
        ] {
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &at, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        if ROWS[index].3 {
            let changed =
                s.at.clone()
                    .with("source-ionic-reference-hydrogen-zero", 0.0)
                    .unwrap();
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &changed, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
            assert!(
                s.pack
                    .claims()
                    .query_typed(
                        &PropertyKey::new("standard-molar-entropy", ENTROPY),
                        &s.at,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
}
