//! Actual NBS packs drive conserved pyrrhotite reagent, heat and Gibbs balances.
//! Endpoint selection and reacted extents are explicit inputs, not kinetics.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, Mass, MolarMass, QuantitySpec,
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
use std::{path::PathBuf, process::Command, sync::OnceLock};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const FES: &str = "FeS(s-alpha-pyrrhotite)";
const FE7S8: &str = "Fe7S8(s-pyrrhotite)";
const ROWS: &[(&str, &str, &str)] = &[
    (FES, "pyrrhotite-fes-alpha-solid-nbs1989", "solid"),
    (FE7S8, "pyrrhotite-fe7s8-solid-nbs1989", "solid"),
    ("S(s)", "sulfur-orthorhombic-solid-nbs1989", "solid"),
    ("Fe[2+](aq)", "iron-ii-ion-aqueous-nbs1989", "ion"),
    ("Fe[3+](aq)", "iron-iii-ion-aqueous-nbs1989", "ion"),
    ("SO4[2-](aq)", "sulfate-ion-aqueous-nbs1989", "ion"),
    ("H+(aq)", "hydrogen-ion-aqueous-nbs1989", "ion"),
    ("H2O(l)", "water-liquid-nbs1989-thermochemistry", "liquid"),
    (
        "O2(aq)",
        "oxygen-aqueous-nbs1989-thermochemistry",
        "aqueous",
    ),
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
        FES => flags.extend([
            ("source-phase-pyrrhotite-fes-alpha", 1.0),
            ("source-iron-atoms-per-sulfur-atom", 1.0),
        ]),
        FE7S8 => flags.extend([
            ("source-phase-pyrrhotite-fe7s8", 1.0),
            ("source-iron-atoms-per-sulfur-atom", 0.875),
        ]),
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
            "frankensim-pyrrhotite-{}-{nonce}",
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

#[test]
fn fifteen_literal_nbs_facts_keep_formula_units_and_reference_phases() {
    for (s, expected) in sources()[..3].iter().zip([
        [0.087911, -100_000.0, -100_400.0, 60.29, 50.54],
        [0.647441, -736_400.0, -748_500.0, 485.8, 398.57],
        [0.032064, 0.0, 0.0, 31.8, 22.64],
    ]) {
        assert_eq!(s.pack.claims().claim_count(), 5);
        for ((property, dims), expected) in [
            ("molar-mass", MolarMass::DIMS),
            ("standard-molar-formation-enthalpy", ENERGY),
            ("standard-molar-formation-gibbs-energy", ENERGY),
            ("standard-molar-entropy", ENTROPY),
            ("standard-molar-isobaric-heat-capacity", ENTROPY),
        ]
        .into_iter()
        .zip(expected)
        {
            close(value(s, property, dims), expected);
        }
    }
}

#[test]
fn new_sources_refuse_wrong_composition_phase_temperature_and_missing_physics() {
    let key = PropertyKey::new("standard-molar-formation-enthalpy", ENERGY);
    for (s, phase) in sources()[..3].iter().zip([
        "source-phase-pyrrhotite-fes-alpha",
        "source-phase-pyrrhotite-fe7s8",
        "source-phase-sulfur-orthorhombic",
    ]) {
        for wrong in [
            QueryPoint::new(),
            s.at.clone().with("temperature", 298.16).unwrap(),
            s.at.clone().with("pressure", 101_325.0).unwrap(),
            s.at.clone().with(phase, 0.0).unwrap(),
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
        for property in [
            "oxidation-rate-constant",
            "bulk-ore-acid-potential",
            "flotation-recovery",
        ] {
            assert!(
                s.pack
                    .claims()
                    .query_typed(
                        &PropertyKey::new(property, Dims::NONE),
                        &s.at,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
    for s in &sources()[..2] {
        for ratio in [0.9, 0.99] {
            let wrong =
                s.at.clone()
                    .with("source-iron-atoms-per-sulfur-atom", ratio)
                    .unwrap();
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
    }
}

#[test]
fn four_endpoints_use_complete_sources_for_mass_heat_and_direct_gibbs() {
    use MineralReaction::*;
    for (recipe, mineral, reacted_kg, extent, expected_h, expected_g) in [
        (
            IronRichPyrrhotiteDissolvedOxygenOxidation,
            FES,
            87.911,
            1000.0,
            -874_970.0,
            -755_830.0,
        ),
        (
            SulfurRichPyrrhotiteDissolvedOxygenOxidation,
            FE7S8,
            647.441,
            500.0,
            -13_388_560.0,
            -11_554_222.0,
        ),
        (
            FerricIronRichPyrrhotiteLeach,
            FES,
            87.911,
            1000.0,
            -70_300.0,
            -126_900.0,
        ),
        (
            FerricSulfurRichPyrrhotiteLeach,
            FE7S8,
            647.441,
            1000.0,
            -455_700.0,
            -842_600.0,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let id = SpeciesId::new(mineral).unwrap();
        let molar_mass = value(source(&id), "molar-mass", MolarMass::DIMS);
        close(
            reaction
                .extent_from_reacted_mass(&id, Mass::new(reacted_kg), MolarMass::new(molar_mass))
                .unwrap()
                .value(),
            extent,
        );
        let h = reaction
            .standard_enthalpy_change(&heat(reaction.terms().iter().map(|t| t.species.clone())))
            .unwrap();
        close(h.molar().value(), expected_h);
        let g_sources = gibbs(&reaction);
        let g = reaction.standard_gibbs_change(&g_sources).unwrap();
        close(g.molar().value(), expected_g);
        let unit_activities: Vec<_> = reaction
            .terms()
            .iter()
            .map(|t| (t.species.clone(), Dimensionless::new(1.0)))
            .collect();
        close(
            g.at_activities(&unit_activities).unwrap().molar.value(),
            expected_g,
        );
        // Even the sulfur reference zero must have its original source.
        assert!(
            reaction
                .standard_gibbs_change(&g_sources[..g_sources.len() - 1])
                .is_err()
        );
    }
}

#[test]
fn both_pyrrhotite_compositions_compete_for_one_finite_dissolved_oxygen_stock() {
    let network = ReactionNetwork::new(vec![
        MineralReaction::IronRichPyrrhotiteDissolvedOxygenOxidation
            .build()
            .unwrap(),
        MineralReaction::SulfurRichPyrrhotiteDissolvedOxygenOxidation
            .build()
            .unwrap(),
    ])
    .unwrap();
    let initial: Vec<_> = network
        .species()
        .iter()
        .map(|id| {
            (
                id.clone(),
                Amount::new(match id.as_str() {
                    FES => 3.0,
                    FE7S8 => 2.0,
                    "O2(aq)" => 20.0,
                    "H2O(l)" => 100.0,
                    _ => 0.0,
                }),
            )
        })
        .collect();
    let extents = |fes: f64, fe7s8: f64| {
        network
            .reactions()
            .iter()
            .map(|r| {
                (
                    r.id().clone(),
                    Amount::new(
                        if r.id().as_str() == "pyrrhotite-fes-oxidation-dissolved-oxygen" {
                            fes
                        } else {
                            fe7s8
                        },
                    ),
                )
            })
            .collect::<Vec<_>>()
    };
    let heat_sources = heat(network.species().iter().cloned());
    let accepted = network
        .advance_with_standard_enthalpy(&initial, &extents(2.0, 0.5), &heat_sources)
        .unwrap();
    close(accepted.enthalpy.value(), -8_444_220.0);
    for (id, expected) in [
        (FES, 1.0),
        (FE7S8, 1.0),
        ("O2(aq)", 0.5),
        ("H2O(l)", 99.0),
        ("Fe[2+](aq)", 9.0),
        ("SO4[2-](aq)", 10.0),
        ("H+(aq)", 2.0),
    ] {
        let actual = accepted
            .inventory
            .final_inventory
            .iter()
            .find(|(s, _)| s.as_str() == id)
            .unwrap()
            .1
            .value();
        close(actual, expected);
    }
    // Each requested conversion is separately feasible, but their shared
    // oxygen demand is 21.5 mol. No air reservoir or rate is invented.
    assert!(
        matches!(network.advance_inventory(&initial, &extents(3.0, 0.5)),
        Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "O2(aq)")
    );
    assert_eq!(
        initial
            .iter()
            .find(|(s, _)| s.as_str() == "O2(aq)")
            .unwrap()
            .1
            .value(),
        20.0
    );
}
