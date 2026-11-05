//! Actual NBS packs drive silver-sulfide dissolution and sulfur cyanide demand.
//! G0/G3: conserved source identities, finite stocks, direct heat/Gibbs and replay.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, Mass, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    ReactionId, SpeciesId,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::FormationGibbsSource,
    },
};
use std::{path::PathBuf, process::Command, sync::OnceLock};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const AG2S: &str = "Ag2S(s-nbs-alpha-orthorhombic)";
const ROWS: &[(&str, &str, &str)] = &[
    (AG2S, "silver-sulfide-alpha-solid-nbs1989", "solid"),
    ("S2O3[2-](aq)", "thiosulfate-ion-aqueous-nbs1989", "ion"),
    ("SCN[-](aq)", "thiocyanate-ion-aqueous-nbs1989", "ion"),
    ("S(s)", "sulfur-orthorhombic-solid-nbs1989", "solid"),
    ("SO3[2-](aq)", "sulfite-ion-aqueous-nbs1989", "ion"),
    (
        "Ag(CN)2[-](aq)",
        "dicyanoargentate-ion-aqueous-nbs1989",
        "ion",
    ),
    ("Au(s)", "gold-solid-nbs1989-thermochemistry", "solid"),
    ("Au(CN)2[-](aq)", "dicyanoaurate-ion-aqueous-nbs1989", "ion"),
    ("CN[-](aq)", "cyanide-ion-aqueous-nbs1989", "ion"),
    ("OH[-](aq)", "hydroxide-ion-aqueous-nbs1989", "ion"),
    (
        "O2(aq)",
        "oxygen-aqueous-nbs1989-thermochemistry",
        "aqueous",
    ),
    ("H2O(l)", "water-liquid-nbs1989-thermochemistry", "liquid"),
];

struct Source {
    id: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
}

fn workspace() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .into()
}

fn point(species: &str, phase: &str) -> QueryPoint {
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
    let aqueous = phase == "ion" || phase == "aqueous";
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (
            match phase {
                "ion" | "aqueous" => "source-phase-aqueous",
                "solid" => "source-phase-solid",
                "liquid" => "source-phase-liquid",
                _ => panic!("unknown source phase"),
            },
            1.0,
        ),
        (
            if aqueous {
                "source-aqueous-unit-activity-infinite-dilution"
            } else {
                "source-pure-species"
            },
            1.0,
        ),
    ];
    if phase == "ion" {
        flags.extend([
            ("source-individual-ion", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]);
    }
    match species {
        AG2S => flags.push(("source-phase-silver-sulfide-nbs-alpha-orthorhombic", 1.0)),
        "S(s)" => flags.push(("source-phase-sulfur-orthorhombic", 1.0)),
        _ => (),
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
            "frankensim-sulfur-cyanide-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        ROWS.iter()
            .map(|(id, slug, phase)| {
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
                let bytes = std::fs::read(&target).unwrap();
                let raw = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(pack.pack_id(), *slug);
                Source {
                    id: SpeciesId::new(*id).unwrap(),
                    pack,
                    at: point(id, phase),
                }
            })
            .collect()
    })
}

fn source(id: &SpeciesId) -> &'static Source {
    sources().iter().find(|s| s.id == *id).unwrap()
}

fn value(s: &Source, name: &str, dims: Dims) -> f64 {
    let answer = s
        .pack
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            &s.at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    s.pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer.evidence.value.value
}

fn heat(ids: impl IntoIterator<Item = SpeciesId>) -> Vec<FormationEnthalpySource<'static>> {
    ids.into_iter()
        .map(|id| {
            let s = source(&id);
            FormationEnthalpySource {
                species: &s.id,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}

fn gibbs(reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'static>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let s = source(&term.species);
            FormationGibbsSource {
                species: &s.id,
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

#[test]
fn fourteen_literal_facts_preserve_phase_erratum_and_negative_conventional_heat_capacity() {
    let expected = [
        [0.2478040, -32_590.0, -40_670.0, 144.01, 76.53],
        [0.1121262, -652_300.0, -522_500.0, 67.0, f64::NAN],
        [0.0580819, 76_440.0, 92_710.0, 144.3, -40.2],
    ];
    for (i, (s, row)) in sources()[..3].iter().zip(expected).enumerate() {
        assert_eq!(s.pack.claims().claim_count(), if i == 1 { 4 } else { 5 });
        for ((key, dims), expected) in [
            ("molar-mass", MolarMass::DIMS),
            ("standard-molar-formation-enthalpy", ENERGY),
            ("standard-molar-formation-gibbs-energy", ENERGY),
            (
                if i == 0 {
                    "standard-molar-entropy"
                } else {
                    "conventional-standard-molar-entropy"
                },
                ENTROPY,
            ),
        ]
        .into_iter()
        .zip(row[..4].iter().copied())
        {
            close(value(s, key, dims), expected);
        }
        if i != 1 {
            close(
                value(
                    s,
                    if i == 0 {
                        "standard-molar-isobaric-heat-capacity"
                    } else {
                        "conventional-standard-molar-isobaric-heat-capacity"
                    },
                    ENTROPY,
                ),
                row[4],
            );
        }
        if i != 0 {
            assert!(
                s.pack
                    .claims()
                    .query_typed(
                        &PropertyKey::new("standard-molar-isobaric-heat-capacity", ENTROPY),
                        &s.at,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err(),
                "no bulk heat capacity is supplied for an aqueous ion"
            );
        }
    }
    assert!(
        sources()[1]
            .pack
            .claims()
            .query_typed(
                &PropertyKey::new(
                    "conventional-standard-molar-isobaric-heat-capacity",
                    ENTROPY
                ),
                &sources()[1].at,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err(),
        "missing thiosulfate Cp remains missing"
    );
}

#[test]
fn new_sources_refuse_other_phases_reference_families_and_extrapolated_conditions() {
    let key = PropertyKey::new("standard-molar-formation-enthalpy", ENERGY);
    for (s, qualifier) in sources()[..3].iter().zip([
        "source-phase-silver-sulfide-nbs-alpha-orthorhombic",
        "source-individual-ion",
        "source-ionic-reference-hydrogen-zero",
    ]) {
        for wrong in [
            QueryPoint::new(),
            s.at.clone().with("temperature", 298.16).unwrap(),
            s.at.clone().with("pressure", 101_325.0).unwrap(),
            s.at.clone().with(qualifier, 0.0).unwrap(),
            s.at.clone().with(REFERENCE, 0.0).unwrap(),
            s.at.clone().with("source-specimen-qualified", 1.0).unwrap(),
        ] {
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        assert!(
            s.pack
                .claims()
                .query_typed(
                    &PropertyKey::new("cyanide-consumption-rate", Dims::NONE),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn matched_sources_drive_three_independent_mass_enthalpy_and_gibbs_balances() {
    for (recipe, reactant, kg, extent, dh, dg, ds) in [
        (
            MineralReaction::SilverSulfideCyanidationThiosulfate,
            AG2S,
            0.2478040,
            0.5,
            -861_478.0,
            -708_519.0,
            -519.03,
        ),
        (
            MineralReaction::ThiosulfateCyanideSulfurTransfer,
            "S2O3[2-](aq)",
            0.1121262,
            1.0,
            -57_360.0,
            -43_690.0,
            -45.8,
        ),
        (
            MineralReaction::ElementalSulfurCyanideConsumption,
            "S(s)",
            0.0320640,
            1.0,
            -74_160.0,
            -79_690.0,
            18.4,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let id = SpeciesId::new(reactant).unwrap();
        let molar_mass = value(source(&id), "molar-mass", MolarMass::DIMS);
        close(
            reaction
                .extent_from_reacted_mass(&id, Mass::new(kg), MolarMass::new(molar_mass))
                .unwrap()
                .value(),
            extent,
        );
        let h_sources = heat(reaction.terms().iter().map(|t| t.species.clone()));
        let h = reaction.standard_enthalpy_change(&h_sources).unwrap();
        close(h.molar().value(), dh);
        let g_sources = gibbs(&reaction);
        let g = reaction.standard_gibbs_change(&g_sources).unwrap();
        close(g.molar().value(), dg);
        let hs = reaction.standard_equilibrium(&h_sources).unwrap();
        close(hs.molar_gibbs().value(), dh - 298.15 * ds);
        assert!(
            (g.molar().value() - hs.molar_gibbs().value()).abs() > 1.0,
            "retain independent tabulated-G discrepancy"
        );
        let mass_residual: f64 = reaction
            .terms()
            .iter()
            .map(|term| {
                term.coefficient as f64
                    * value(source(&term.species), "molar-mass", MolarMass::DIMS)
            })
            .sum();
        close(mass_residual, 0.0);
        assert!(
            reaction
                .standard_gibbs_change(&g_sources[..g_sources.len() - 1])
                .is_err()
        );
        let mut wrong = g_sources;
        wrong[0].reference_axis = "source-elemental-reference-usgs-b2131";
        assert!(reaction.standard_gibbs_change(&wrong).is_err());
    }
}

fn network() -> ReactionNetwork {
    ReactionNetwork::new(
        [
            MineralReaction::GoldCyanidationDissolvedOxygen,
            MineralReaction::SilverSulfideCyanidationThiosulfate,
            MineralReaction::ElementalSulfurCyanideConsumption,
            MineralReaction::ThiosulfateCyanideSulfurTransfer,
        ]
        .into_iter()
        .map(|r| r.build().unwrap())
        .collect(),
    )
    .unwrap()
}

fn requests(gold: f64, silver: f64, sulfur: f64, thiosulfate: f64) -> Vec<(ReactionId, Amount)> {
    [
        ("gold-cyanidation-elsner-dissolved-oxygen", gold),
        ("silver-sulfide-cyanidation-thiosulfate", silver),
        ("elemental-sulfur-cyanide-consumption", sulfur),
        ("thiosulfate-cyanide-sulfur-transfer", thiosulfate),
    ]
    .into_iter()
    .map(|(id, n)| (ReactionId::new(id).unwrap(), Amount::new(n)))
    .collect()
}

fn initial(network: &ReactionNetwork, cyanide: f64, oxygen: f64) -> Vec<(SpeciesId, Amount)> {
    network
        .species()
        .iter()
        .map(|id| {
            let n = match id.as_str() {
                "Au(s)" => 4.0,
                AG2S => 2.0,
                "S(s)" => 1.0,
                "CN[-](aq)" => cyanide,
                "O2(aq)" => oxygen,
                "H2O(l)" => 10.0,
                _ => 0.0, // Every initial product is explicitly zero in this authored fixture.
            };
            (id.clone(), Amount::new(n))
        })
        .collect()
}

fn amount(stock: &[(SpeciesId, Amount)], name: &str) -> f64 {
    stock
        .iter()
        .find(|(id, _)| id.as_str() == name)
        .unwrap()
        .1
        .value()
}

#[test]
fn finite_gold_silver_and_sulfur_stages_preserve_all_cyanide_units_and_source_mass() {
    let network = network();
    let initial = initial(&network, 18.0, 3.0);
    let original = initial.clone();
    let heat = heat(network.species().iter().cloned());
    let first = network
        .advance_with_standard_enthalpy(&initial, &requests(1.0, 1.0, 1.0, 0.0), &heat)
        .unwrap();
    let after_first = &first.inventory.final_inventory;
    close(first.enthalpy.value(), -1_507_854.0);
    for (name, expected) in [
        ("CN[-](aq)", 1.0),
        ("O2(aq)", 0.0),
        ("S2O3[2-](aq)", 1.0),
        ("SCN[-](aq)", 1.0),
    ] {
        close(amount(after_first, name), expected);
    }
    let second = network
        .advance_with_standard_enthalpy(after_first, &requests(0.0, 0.0, 0.0, 1.0), &heat)
        .unwrap();
    close(second.enthalpy.value(), -57_360.0);
    close(
        first.enthalpy.value() + second.enthalpy.value(),
        -1_565_214.0,
    );
    let final_stock = &second.inventory.final_inventory;
    for (name, expected) in [
        ("Au(s)", 0.0),
        (AG2S, 0.0),
        ("S(s)", 0.0),
        ("CN[-](aq)", 0.0),
        ("O2(aq)", 0.0),
        ("S2O3[2-](aq)", 0.0),
        ("Au(CN)2[-](aq)", 4.0),
        ("Ag(CN)2[-](aq)", 4.0),
        ("SCN[-](aq)", 2.0),
        ("SO3[2-](aq)", 1.0),
        ("OH[-](aq)", 6.0),
        ("H2O(l)", 7.0),
    ] {
        close(amount(final_stock, name), expected);
    }
    for stock in [&initial, after_first, final_stock] {
        let n = |name| amount(stock, name);
        close(n("Au(s)") + n("Au(CN)2[-](aq)"), 4.0);
        close(2.0 * n(AG2S) + n("Ag(CN)2[-](aq)"), 4.0);
        close(
            n(AG2S) + n("S(s)") + 2.0 * n("S2O3[2-](aq)") + n("SCN[-](aq)") + n("SO3[2-](aq)"),
            3.0,
        );
        close(
            n("CN[-](aq)")
                + 2.0 * n("Au(CN)2[-](aq)")
                + 2.0 * n("Ag(CN)2[-](aq)")
                + n("SCN[-](aq)"),
            18.0,
        );
        close(
            2.0 * n("O2(aq)")
                + n("H2O(l)")
                + n("OH[-](aq)")
                + 3.0 * n("S2O3[2-](aq)")
                + 3.0 * n("SO3[2-](aq)"),
            16.0,
        );
        close(2.0 * n("H2O(l)") + n("OH[-](aq)"), 20.0);
        close(
            -n("CN[-](aq)")
                - n("Au(CN)2[-](aq)")
                - n("Ag(CN)2[-](aq)")
                - n("OH[-](aq)")
                - n("SCN[-](aq)")
                - 2.0 * n("S2O3[2-](aq)")
                - 2.0 * n("SO3[2-](aq)"),
            -18.0,
        );
        let mass: f64 = stock
            .iter()
            .map(|(id, n)| n.value() * value(source(id), "molar-mass", MolarMass::DIMS))
            .sum();
        close(mass, 2.0600126);
    }
    assert_eq!(initial, original);
    let replay = network
        .advance_with_standard_enthalpy(&initial, &requests(1.0, 1.0, 1.0, 0.0), &heat)
        .unwrap();
    assert_eq!(
        first.inventory.final_inventory,
        replay.inventory.final_inventory
    );
    assert_eq!(first.enthalpy, replay.enthalpy);
}

#[test]
fn shared_reagents_and_new_thiosulfate_cannot_fund_unadmitted_withdrawals() {
    let network = network();
    let heat = heat(network.species().iter().cloned());
    for (cn, oxygen, extents, missing) in [
        (16.0, 3.0, [1.0, 1.0, 1.0, 0.0], "CN[-](aq)"),
        (18.0, 2.0, [1.0, 1.0, 0.0, 0.0], "O2(aq)"),
        (18.0, 3.0, [0.0, 1.0, 0.0, 1.0], "S2O3[2-](aq)"),
    ] {
        let stock = initial(&network, cn, oxygen);
        let before = stock.clone();
        let [au, ag, sulfur, thio] = extents;
        assert!(matches!(
            network.advance_with_standard_enthalpy(&stock, &requests(au, ag, sulfur, thio), &heat),
            Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == missing
        ));
        assert_eq!(stock, before);
    }
    // Each metal's oxygen demand is feasible on its own, but their sum is not.
    let limited_oxygen = initial(&network, 18.0, 2.0);
    network
        .advance_inventory(&limited_oxygen, &requests(1.0, 0.0, 0.0, 0.0))
        .unwrap();
    network
        .advance_inventory(&limited_oxygen, &requests(0.0, 1.0, 0.0, 0.0))
        .unwrap();
    // Thiocyanate retains C/N atoms; it does not replenish free CN in this model.
    let exhausted = network
        .advance_inventory(&initial(&network, 17.0, 3.0), &requests(1.0, 1.0, 1.0, 0.0))
        .unwrap();
    close(amount(&exhausted.final_inventory, "SCN[-](aq)"), 1.0);
    assert!(matches!(
        network.advance_inventory(&exhausted.final_inventory, &requests(0.0, 0.0, 0.0, 1.0)),
        Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "CN[-](aq)"
    ));
}
