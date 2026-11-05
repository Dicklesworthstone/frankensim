//! G0/G1/G3: original NBS sulfite sources, conserved SO2/air endpoints and
//! finite competing demands. These checks do not establish process kinetics.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpeciationOptions, SpectatorIon,
                WaterAutoprotolysis,
            },
        },
    },
};
use std::{path::PathBuf, process::Command, sync::OnceLock};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
// Original source phases: aqueous ion, aqueous neutral, liquid or solid.
const ROWS: &[(&str, &str, &str)] = &[
    ("HSO3[-](aq)", "hydrogen-sulfite-ion-aqueous-nbs1989", "ion"),
    ("SO3[2-](aq)", "sulfite-ion-aqueous-nbs1989", "ion"),
    (
        "Na2S2O5(s)",
        "sodium-metabisulfite-solid-nbs1989-thermochemistry",
        "solid",
    ),
    (
        "Na2SO3(s)",
        "sodium-sulfite-solid-nbs1989-thermochemistry",
        "solid",
    ),
    ("CN[-](aq)", "cyanide-ion-aqueous-nbs1989", "ion"),
    ("OCN[-](aq)", "cyanate-ion-aqueous-nbs1989", "ion"),
    ("SO4[2-](aq)", "sulfate-ion-aqueous-nbs1989", "ion"),
    (
        "O2(aq)",
        "oxygen-aqueous-nbs1989-thermochemistry",
        "aqueous",
    ),
    ("Na[+](aq)", "sodium-ion-aqueous-nbs1989", "ion"),
    ("H+(aq)", "hydrogen-ion-aqueous-nbs1989", "ion"),
    ("H2O(l)", "water-liquid-nbs1989-thermochemistry", "liquid"),
    ("OH[-](aq)", "hydroxide-ion-aqueous-nbs1989", "ion"),
];

struct Source {
    id: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
}

fn workspace() -> PathBuf {
    std::env::var_os("CARGO_WORKSPACE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .into()
        })
}

fn point(phase: &str) -> QueryPoint {
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
    let phase_axis = match phase {
        "ion" | "aqueous" => "source-phase-aqueous",
        "solid" => "source-phase-solid",
        "liquid" => "source-phase-liquid",
        _ => panic!("unexpected test phase"),
    };
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        (phase_axis, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
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
        let out =
            std::env::temp_dir().join(format!("frankensim-sulfite-{}-{nonce}", std::process::id()));
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
                let bytes = std::fs::read(target).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), *slug);
                assert_eq!(pack.to_bytes(), bytes);
                Source {
                    id: SpeciesId::new(*id).unwrap(),
                    pack,
                    at: point(phase),
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

#[test]
fn fifteen_source_facts_retain_literal_values_units_conventions_and_receipts() {
    for (index, literals) in [
        (0, [0.0810702, -626_220.0, -527_730.0, 139.7]),
        (1, [0.0800622, -635_500.0, -486_500.0, -29.0]),
        (3, [0.1260418, -1_100_800.0, -1_012_500.0, 145.94]),
    ] {
        let s = &sources()[index];
        assert_eq!(
            s.pack.claims().claim_count(),
            if index == 3 { 5 } else { 4 }
        );
        for (name, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, literals[0]),
            ("standard-molar-formation-enthalpy", ENERGY, literals[1]),
            ("standard-molar-formation-gibbs-energy", ENERGY, literals[2]),
            (
                if index == 3 {
                    "standard-molar-entropy"
                } else {
                    "conventional-standard-molar-entropy"
                },
                ENTROPY,
                literals[3],
            ),
        ] {
            check_value(s, name, dims, expected);
        }
    }
    let meta = &sources()[2];
    assert_eq!(meta.pack.claims().claim_count(), 2);
    check_value(meta, "molar-mass", MolarMass::DIMS, 0.1901046);
    check_value(
        meta,
        "standard-molar-formation-enthalpy",
        ENERGY,
        -1_478_200.0,
    );
    check_value(
        &sources()[3],
        "standard-molar-isobaric-heat-capacity",
        ENTROPY,
        120.25,
    );
}

fn check_value(s: &Source, name: &str, dims: Dims, expected: f64) {
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

#[test]
fn absent_thermodynamics_and_unsupported_reference_states_do_not_become_numbers() {
    for (i, s) in sources()[..4].iter().enumerate() {
        let key = PropertyKey::new("standard-molar-formation-enthalpy", ENERGY);
        for wrong in [
            QueryPoint::new(),
            s.at.clone().with("temperature", 310.0).unwrap(),
            s.at.clone().with(REFERENCE, 0.0).unwrap(),
            s.at.clone()
                .with(
                    if i < 2 {
                        "source-phase-aqueous"
                    } else {
                        "source-phase-solid"
                    },
                    0.0,
                )
                .unwrap(),
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
                    &PropertyKey::new(key.name(), Dims::NONE),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
    let meta = &sources()[2];
    for (name, dims) in [
        ("standard-molar-formation-gibbs-energy", ENERGY),
        ("standard-molar-entropy", ENTROPY),
        ("standard-molar-isobaric-heat-capacity", ENTROPY),
        ("cyanide-removal-efficiency", Dims::NONE),
    ] {
        assert!(
            meta.pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &meta.at,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
    let reaction = MineralReaction::SodiumMetabisulfiteDissolution
        .build()
        .unwrap();
    assert!(reaction.standard_gibbs_change(&gibbs(&reaction)).is_err());
}

#[test]
fn source_bound_heat_and_equilibrium_preserve_preparation_and_protonation() {
    use MineralReaction::*;
    // Independent sums of the original NBS Hf literals, J per written equation.
    for (recipe, expected) in [
        (SodiumMetabisulfiteDissolution, 31_350.0),
        (SodiumSulfiteDissolution, -14_940.0),
        (HydrogenSulfiteDissociation, -9_280.0),
        (CyanideBisulfiteOxidation, -567_950.0),
        (CyanideSulfiteOxidation, -558_670.0),
        (SulfiteOxygenOxidation, -535_840.0),
        (BisulfiteOxygenOxidation, -554_400.0),
    ] {
        let reaction = recipe.build().unwrap();
        let h = reaction
            .standard_enthalpy_change(&heat(reaction.terms().iter().map(|t| t.species.clone())))
            .unwrap();
        close(h.molar().value(), expected);
    }
    let reaction = HydrogenSulfiteDissociation.build().unwrap();
    let step = AcidDissociation::from_sources(
        reaction.clone(),
        SpeciesId::new("HSO3[-](aq)").unwrap(),
        SpeciesId::new("SO3[2-](aq)").unwrap(),
        SpeciesId::new("H+(aq)").unwrap(),
        &gibbs(&reaction),
    )
    .unwrap();
    close(step.standard().molar().value(), 41_230.0);
    close(
        step.standard().log_equilibrium_constant().value(),
        -41_230.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15),
    );
}

fn network() -> ReactionNetwork {
    use MineralReaction::*;
    ReactionNetwork::new(
        [
            SodiumMetabisulfiteDissolution,
            CyanideBisulfiteOxidation,
            BisulfiteOxygenOxidation,
        ]
        .into_iter()
        .map(|r| r.build().unwrap())
        .collect(),
    )
    .unwrap()
}
fn inventory(network: &ReactionNetwork, values: &[(&str, f64)]) -> Vec<(SpeciesId, Amount)> {
    network
        .species()
        .iter()
        .map(|id| {
            (
                id.clone(),
                Amount::new(
                    values
                        .iter()
                        .find(|(s, _)| *s == id.as_str())
                        .map_or(0.0, |(_, n)| *n),
                ),
            )
        })
        .collect()
}
fn extents(
    network: &ReactionNetwork,
    prepare: f64,
    productive: f64,
    side: f64,
) -> Vec<(fs_thermochem::ReactionId, Amount)> {
    network
        .reactions()
        .iter()
        .map(|r| {
            (
                r.id().clone(),
                Amount::new(match r.id().as_str() {
                    "sodium-metabisulfite-dissolution" => prepare,
                    "cyanide-bisulfite-oxidation" => productive,
                    "bisulfite-oxygen-oxidation" => side,
                    _ => panic!("unexpected test pathway"),
                }),
            )
        })
        .collect()
}
fn amount(values: &[(SpeciesId, Amount)], name: &str) -> f64 {
    values
        .iter()
        .find(|(id, _)| id.as_str() == name)
        .unwrap()
        .1
        .value()
}

#[test]
fn staged_metabisulfite_treatment_accounts_for_finite_shared_oxygen_and_acid() {
    let network = network();
    let initial = inventory(
        &network,
        &[
            ("Na2S2O5(s)", 4.0),
            ("H2O(l)", 100.0),
            ("CN[-](aq)", 10.0),
            ("O2(aq)", 7.0),
        ],
    );
    let h = heat(network.species().iter().cloned());
    let prepared = network
        .advance_with_standard_enthalpy(&initial, &extents(&network, 4.0, 0.0, 0.0), &h)
        .unwrap();
    close(prepared.enthalpy.value(), 125_400.0);
    let treated = network
        .advance_with_standard_enthalpy(
            &prepared.inventory.final_inventory,
            &extents(&network, 0.0, 6.0, 1.0),
            &h,
        )
        .unwrap();
    for (species, expected) in [
        ("CN[-](aq)", 4.0),
        ("OCN[-](aq)", 6.0),
        ("SO4[2-](aq)", 8.0),
        ("H+(aq)", 8.0),
        ("Na[+](aq)", 8.0),
        ("O2(aq)", 0.0),
        ("HSO3[-](aq)", 0.0),
        ("H2O(l)", 96.0),
    ] {
        close(
            amount(&treated.inventory.final_inventory, species),
            expected,
        );
    }
    close(treated.enthalpy.value(), -3_962_100.0);
    // C and N stay in the explicit cyanide/cyanate pool; oxidation is not removal.
    close(
        amount(&treated.inventory.final_inventory, "CN[-](aq)")
            + amount(&treated.inventory.final_inventory, "OCN[-](aq)"),
        10.0,
    );
    let mut reversed = prepared.inventory.final_inventory.clone();
    reversed.reverse();
    let replay = network
        .advance_with_standard_enthalpy(&reversed, &extents(&network, 0.0, 6.0, 1.0), &h)
        .unwrap();
    assert_eq!(
        treated.inventory.final_inventory,
        replay.inventory.final_inventory
    );
    assert_eq!(
        treated.enthalpy.value().to_bits(),
        replay.enthalpy.value().to_bits()
    );
}

#[test]
fn same_stage_reagent_borrowing_and_shared_oxygen_overdraw_refuse_atomically() {
    let network = network();
    let initial = inventory(
        &network,
        &[
            ("Na2S2O5(s)", 4.0),
            ("H2O(l)", 100.0),
            ("CN[-](aq)", 10.0),
            ("O2(aq)", 7.0),
        ],
    );
    let saved = initial.clone();
    assert!(
        matches!(network.advance_inventory(&initial, &extents(&network, 4.0, 6.0, 1.0)), Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "HSO3[-](aq)")
    );
    let prepared = network
        .advance_inventory(&initial, &extents(&network, 4.0, 0.0, 0.0))
        .unwrap();
    assert!(
        matches!(network.advance_inventory(&prepared.final_inventory, &extents(&network, 0.0, 8.0, 0.0)), Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "O2(aq)")
    );
    assert!(matches!(
        network.advance_inventory(&prepared.final_inventory, &extents(&network, 0.0, 7.0, 1.0)),
        Err(ReactionNetworkError::InsufficientInventory { .. })
    ));
    let dry = inventory(&network, &[("Na2S2O5(s)", 4.0), ("H2O(l)", 3.0)]);
    assert!(
        matches!(network.advance_inventory(&dry, &extents(&network, 4.0, 0.0, 0.0)), Err(ReactionNetworkError::InsufficientInventory { species, .. }) if species.as_str() == "H2O(l)")
    );
    assert_eq!(initial, saved);
}

#[test]
fn source_sulfite_family_partitions_at_one_charge_balanced_ph() {
    let id = |s| SpeciesId::new(s).unwrap();
    let acid = MineralReaction::HydrogenSulfiteDissociation
        .build()
        .unwrap();
    let step = AcidDissociation::from_sources(
        acid.clone(),
        id("HSO3[-](aq)"),
        id("SO3[2-](aq)"),
        id("H+(aq)"),
        &gibbs(&acid),
    )
    .unwrap();
    let water = MineralReaction::WaterDissociation.build().unwrap();
    let water = WaterAutoprotolysis::from_sources(
        water.clone(),
        id("H2O(l)"),
        id("H+(aq)"),
        id("OH[-](aq)"),
        &gibbs(&water),
    )
    .unwrap();
    let model = IdealAcidBaseSystem::new(
        water,
        vec![AcidFamily::new("dissolved-sulfur-iv", vec![step]).unwrap()],
    )
    .unwrap();
    let answer = model
        .speciate(
            &[("dissolved-sulfur-iv".to_owned(), Amount::new(2.0))],
            &[SpectatorIon {
                species: id("Na[+](aq)"),
                charge: 1,
                amount: Amount::new(3.0),
            }],
            Amount::new(1000.0),
            SpeciationOptions {
                ph_bracket: [Dimensionless::new(6.0), Dimensionless::new(8.0)],
                ph_tolerance: Dimensionless::new(1e-11),
                ..SpeciationOptions::default()
            },
            || false,
        )
        .unwrap();
    let ka = (-41_230.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15)).exp();
    let proton = answer.proton_activity.value();
    for form in &answer.families[0].forms {
        let expected = if form.species.as_str() == "HSO3[-](aq)" {
            2.0 * proton / (proton + ka)
        } else {
            2.0 * ka / (proton + ka)
        };
        close(form.amount.value(), expected);
    }
    close(
        answer.families[0]
            .forms
            .iter()
            .map(|f| f.amount.value())
            .sum(),
        2.0,
    );
    assert!(answer.charge_balance_residual.value().abs() < 1e-9);
    // Only the named HSO3-/SO3(2-) family is admitted, not low-pH SO2 speciation.
    assert_eq!(answer.families[0].forms.len(), 2);
}
