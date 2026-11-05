//! Actual NBS packs drive finite copper/clay acid demand and staged Al hydrolysis.
//! Known-answer extents are constructed inputs, not a clay-reactivity prediction.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, Mass, MolarMass, QuantitySpec,
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
const KAOLINITE: &str = "Al2Si2O5(OH)4(s-kaolinite)";
const GIBBSITE: &str = "Al2O3.3H2O(s-gibbsite)";
const ROWS: &[(&str, &str, &str)] = &[
    (
        KAOLINITE,
        "kaolinite-solid-nbs1989-thermochemistry",
        "solid",
    ),
    ("Al[3+](aq)", "aluminum-iii-ion-aqueous-nbs1989", "ion"),
    (
        "H4SiO4(aq)",
        "silicic-acid-component-aqueous-nbs1989",
        "aqueous",
    ),
    (GIBBSITE, "gibbsite-solid-nbs1989-thermochemistry", "solid"),
    (
        "CuO(s)",
        "copper-ii-oxide-solid-nbs1989-thermochemistry",
        "solid",
    ),
    ("Cu[2+](aq)", "copper-ii-ion-aqueous-nbs1989", "ion"),
    ("H+(aq)", "hydrogen-ion-aqueous-nbs1989", "ion"),
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
                "liquid" => "source-phase-liquid",
                "solid" => "source-phase-solid",
                _ => panic!("unknown phase"),
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
        KAOLINITE => flags.push(("source-phase-kaolinite", 1.0)),
        GIBBSITE => flags.extend([
            ("source-phase-gibbsite", 1.0),
            ("source-formula-aluminum-atoms", 2.0),
        ]),
        "H4SiO4(aq)" => flags.push(("source-silicic-component-h2sio3-plus-water", 1.0)),
        "CuO(s)" => flags.push(("source-phase-copper-ii-oxide-crystalline", 1.0)),
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
            "frankensim-clay-acid-{}-{nonce}",
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

fn heat(species: impl IntoIterator<Item = SpeciesId>) -> Vec<FormationEnthalpySource<'static>> {
    species
        .into_iter()
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
fn twenty_three_nbs_facts_preserve_phase_component_and_formula_unit_basis() {
    let expected = [
        [0.2581616, -4_119_600.0, -3_799_700.0, 205.0, 246.14],
        [0.0269815, -531_000.0, -485_000.0, -321.7, f64::NAN],
        [0.0961156, -1_468_600.0, -1_316_600.0, 180.0, f64::NAN],
        [0.1560074, -2_586_670.0, -2_310_210.0, 136.90, 183.47],
        [0.0795394, -157_300.0, -129_700.0, 42.63, 42.30],
    ];
    for ((s, row), index) in sources()[..5].iter().zip(expected).zip(0..) {
        assert_eq!(
            s.pack.claims().claim_count(),
            if row[4].is_nan() { 4 } else { 5 }
        );
        for ((key, dims), expected) in [
            ("molar-mass", MolarMass::DIMS),
            ("standard-molar-formation-enthalpy", ENERGY),
            ("standard-molar-formation-gibbs-energy", ENERGY),
            (
                if index == 1 {
                    "conventional-standard-molar-entropy"
                } else {
                    "standard-molar-entropy"
                },
                ENTROPY,
            ),
        ]
        .into_iter()
        .zip(row[..4].iter().copied())
        {
            close(value(s, key, dims), expected);
        }
        if row[4].is_nan() {
            assert!(
                s.pack
                    .claims()
                    .query_typed(
                        &PropertyKey::new("standard-molar-isobaric-heat-capacity", ENTROPY),
                        &s.at,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        } else {
            close(
                value(s, "standard-molar-isobaric-heat-capacity", ENTROPY),
                row[4],
            );
        }
    }
}

#[test]
fn source_domains_refuse_wrong_hydration_basis_phase_reference_and_inferred_rates() {
    let key = PropertyKey::new("standard-molar-formation-enthalpy", ENERGY);
    for (s, qualifier) in sources()[..5].iter().zip([
        "source-phase-kaolinite",
        "source-individual-ion",
        "source-silicic-component-h2sio3-plus-water",
        "source-phase-gibbsite",
        "source-phase-copper-ii-oxide-crystalline",
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
        for key in [
            "acid-consumption-rate",
            "ore-leach-recovery",
            "precipitation-rate-constant",
        ] {
            assert!(
                s.pack
                    .claims()
                    .query_typed(
                        &PropertyKey::new(key, Dims::NONE),
                        &s.at,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
    let gibbsite = &sources()[3];
    let per_aloh3 = gibbsite
        .at
        .clone()
        .with("source-formula-aluminum-atoms", 1.0)
        .unwrap();
    assert!(
        gibbsite
            .pack
            .claims()
            .query_typed(&key, &per_aloh3, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
}

#[test]
fn source_bound_mass_heat_and_gibbs_match_three_independent_written_equations() {
    for (recipe, first, kg, extent, dh, dg) in [
        (
            MineralReaction::KaoliniteAcidLeach,
            KAOLINITE,
            0.2581616,
            1.0,
            -165_430.0,
            -40_629.0,
        ),
        (
            MineralReaction::AluminumHydrolysisGibbsite,
            "Al[3+](aq)",
            0.0269815,
            0.5,
            190_310.0,
            82_564.0,
        ),
        (
            MineralReaction::CupricOxideAcidLeach,
            "CuO(s)",
            0.0795394,
            1.0,
            -63_760.0,
            -41_939.0,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let id = SpeciesId::new(first).unwrap();
        let molar_mass = value(source(&id), "molar-mass", MolarMass::DIMS);
        close(
            reaction
                .extent_from_reacted_mass(&id, Mass::new(kg), MolarMass::new(molar_mass))
                .unwrap()
                .value(),
            extent,
        );
        let h = reaction
            .standard_enthalpy_change(&heat(reaction.terms().iter().map(|t| t.species.clone())))
            .unwrap();
        close(h.molar().value(), dh);
        let g_sources = gibbs(&reaction);
        let g = reaction.standard_gibbs_change(&g_sources).unwrap();
        close(g.molar().value(), dg);
        let activities: Vec<_> = reaction
            .terms()
            .iter()
            .map(|t| (t.species.clone(), Dimensionless::new(1.0)))
            .collect();
        close(g.at_activities(&activities).unwrap().molar.value(), dg);
        assert!(
            reaction
                .standard_gibbs_change(&g_sources[..g_sources.len() - 1])
                .is_err()
        );
        let mut mismatch = g_sources;
        mismatch[0].reference_axis = "source-elemental-reference-usgs-b2131";
        assert!(reaction.standard_gibbs_change(&mismatch).is_err());
    }
}

fn network() -> ReactionNetwork {
    ReactionNetwork::new(
        [
            MineralReaction::CupricOxideAcidLeach,
            MineralReaction::KaoliniteAcidLeach,
            MineralReaction::AluminumHydrolysisGibbsite,
        ]
        .into_iter()
        .map(|r| r.build().unwrap())
        .collect(),
    )
    .unwrap()
}

fn requests(copper: f64, clay: f64, gibbsite: f64) -> Vec<(ReactionId, Amount)> {
    [
        ("cupric-oxide-acid-leach", copper),
        ("kaolinite-acid-leach", clay),
        ("aluminum-hydrolysis-gibbsite", gibbsite),
    ]
    .into_iter()
    .map(|(id, n)| (ReactionId::new(id).unwrap(), Amount::new(n)))
    .collect()
}

fn inventory(network: &ReactionNetwork, stocks: &[(&str, f64)]) -> Vec<(SpeciesId, Amount)> {
    network
        .species()
        .iter()
        .map(|id| {
            let n = stocks
                .iter()
                .find(|(name, _)| *name == id.as_str())
                .map_or(0.0, |(_, n)| *n);
            (id.clone(), Amount::new(n))
        })
        .collect()
}

fn amount(inventory: &[(SpeciesId, Amount)], name: &str) -> f64 {
    inventory
        .iter()
        .find(|(id, _)| id.as_str() == name)
        .unwrap()
        .1
        .value()
}

#[test]
fn clay_and_copper_compete_for_acid_and_only_later_stages_receive_hydrolysis_credit() {
    let network = network();
    let initial = inventory(
        &network,
        &[
            ("CuO(s)", 2.0),
            (KAOLINITE, 1.0),
            ("H+(aq)", 5.0),
            ("H2O(l)", 100.0),
        ],
    );
    let original = initial.clone();
    let heat = heat(network.species().iter().cloned());
    // Individually feasible extents do not authorize their seven-mole combined acid demand.
    network
        .advance_inventory(&initial, &requests(2.0, 0.0, 0.0))
        .unwrap();
    network
        .advance_inventory(&initial, &requests(0.0, 0.5, 0.0))
        .unwrap();
    assert!(
        matches!(network.advance_with_standard_enthalpy(&initial, &requests(2.0, 0.5, 0.0), &heat),
        Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "H+(aq)")
    );
    assert_eq!(initial, original);
    let leached = network
        .advance_with_standard_enthalpy(&initial, &requests(1.0, 0.5, 0.0), &heat)
        .unwrap();
    close(leached.enthalpy.value(), -146_475.0);
    let after_leach = &leached.inventory.final_inventory;
    for (id, expected) in [
        ("H+(aq)", 0.0),
        ("Al[3+](aq)", 1.0),
        ("H4SiO4(aq)", 1.0),
        ("Cu[2+](aq)", 1.0),
        ("H2O(l)", 101.5),
    ] {
        close(amount(after_leach, id), expected);
    }
    // Net acid would be positive, but spending the three protons formed later
    // is forbidden in this simultaneous stage's original zero-proton stock.
    assert!(
        matches!(network.advance_inventory(after_leach, &requests(1.0, 0.0, 0.5)),
        Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "H+(aq)")
    );
    let precipitated = network
        .advance_with_standard_enthalpy(after_leach, &requests(0.0, 0.0, 0.5), &heat)
        .unwrap();
    close(precipitated.enthalpy.value(), 95_155.0);
    close(
        amount(&precipitated.inventory.final_inventory, "H+(aq)"),
        3.0,
    );
    let later = network
        .advance_with_standard_enthalpy(
            &precipitated.inventory.final_inventory,
            &requests(1.0, 0.0, 0.0),
            &heat,
        )
        .unwrap();
    for (id, expected) in [
        ("CuO(s)", 0.0),
        (KAOLINITE, 0.5),
        ("H+(aq)", 1.0),
        ("H2O(l)", 99.5),
        ("Cu[2+](aq)", 2.0),
        ("Al[3+](aq)", 0.0),
        ("H4SiO4(aq)", 1.0),
        (GIBBSITE, 0.5),
    ] {
        close(amount(&later.inventory.final_inventory, id), expected);
    }
    close(
        leached.enthalpy.value() + precipitated.enthalpy.value() + later.enthalpy.value(),
        -115_080.0,
    );
    let mass = |stock: &[(SpeciesId, Amount)]| {
        stock
            .iter()
            .map(|(id, n)| n.value() * value(source(id), "molar-mass", MolarMass::DIMS))
            .sum::<f64>()
    };
    close(mass(&initial), 2.2238204);
    close(mass(&later.inventory.final_inventory), 2.2238204);
    // The explicit final elemental inventory retains the unreacted clay and
    // selected precipitate; no aluminum/silicon mass becomes a loss or plug.
    let final_stock = &later.inventory.final_inventory;
    close(
        2.0 * amount(final_stock, KAOLINITE)
            + amount(final_stock, "Al[3+](aq)")
            + 2.0 * amount(final_stock, GIBBSITE),
        2.0,
    );
    close(
        2.0 * amount(final_stock, KAOLINITE) + amount(final_stock, "H4SiO4(aq)"),
        2.0,
    );
}

#[test]
fn selected_gibbsite_extent_cannot_borrow_aluminum_or_water() {
    let network = network();
    for (stocks, expected) in [
        (
            vec![(KAOLINITE, 1.0), ("H+(aq)", 6.0), ("H2O(l)", 100.0)],
            "Al[3+](aq)",
        ),
        (vec![("Al[3+](aq)", 1.0), ("H2O(l)", 2.0)], "H2O(l)"),
    ] {
        let initial = inventory(&network, &stocks);
        let before = initial.clone();
        // In the first case clay could release Al, but the hydrolysis request
        // cannot consume products from another reaction in the same stage.
        let clay = if expected == "Al[3+](aq)" { 0.5 } else { 0.0 };
        assert!(
            matches!(network.advance_inventory(&initial, &requests(0.0, clay, 0.5)),
            Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == expected)
        );
        assert_eq!(initial, before);
    }
}
