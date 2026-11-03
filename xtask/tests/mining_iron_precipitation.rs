//! G0/G1/G3: phase-specific ferric precipitation closes the compatible lime cycle.
//! Hydrolysis endpoint: https://pubs.usgs.gov/of/2005/1148/of2005-1148.pdf.
#![deny(unsafe_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource, ReactionBalanceError, ReactionEnthalpy,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SLUGS: [&str; 6] = [
    "goethite-solid-usgs-b2131-thermochemistry",
    "portlandite-solid-usgs-b2131-thermochemistry",
    "iron-iii-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
    "hydrogen-ion-aqueous-usgs-b2131",
    "calcium-ion-aqueous-usgs-b2131",
];
const SPECIES: [&str; 6] = [
    "FeO(OH)(s)",
    "Ca(OH)2(s)",
    "Fe[3+](aq)",
    "H2O(l)",
    "H+(aq)",
    "Ca[2+](aq)",
];

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
        "frankensim-iron-precipitation-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("material.fsmatpk");
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
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    pack
}

fn point(index: usize) -> QueryPoint {
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
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    match index {
        0 | 1 => flags.extend([
            ("source-phase-solid", 1.0),
            (
                if index == 0 {
                    "source-phase-goethite"
                } else {
                    "source-phase-portlandite"
                },
                1.0,
            ),
            ("source-pure-species", 1.0),
        ]),
        3 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn heat(
    reaction: &ConservedReaction,
    packs: &[NormalizedPack; 6],
    points: &[QueryPoint; 6],
) -> ReactionEnthalpy {
    let sources = reaction
        .terms()
        .iter()
        .map(|term| {
            let index = SPECIES
                .iter()
                .position(|name| *name == term.species.as_str())
                .unwrap();
            FormationEnthalpySource {
                species: &term.species,
                pack: &packs[index],
                at: &points[index],
                reference_axis: REFERENCE,
            }
        })
        .collect::<Vec<_>>();
    reaction.standard_enthalpy_change(&sources).unwrap()
}

fn gibbs_sources<'a>(
    reaction: &'a ConservedReaction,
    packs: &'a [NormalizedPack; 6],
    points: &'a [QueryPoint; 6],
) -> Vec<FormationGibbsSource<'a>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let index = SPECIES
                .iter()
                .position(|name| *name == term.species.as_str())
                .unwrap();
            FormationGibbsSource {
                species: &term.species,
                pack: &packs[index],
                at: &points[index],
                reference_axis: REFERENCE,
            }
        })
        .collect()
}

#[test]
fn g0_g3_source_values_require_the_named_solid_phase_and_reference() {
    for (index, values) in [
        (0, [0.088854, -562_600.0, -491_800.0, 60.4, 2100.0, 0.6]),
        (1, [0.074093, -986_100.0, -898_000.0, 83.4, 1300.0, 0.4]),
    ] {
        let pack = compile(SLUGS[index]);
        let query = point(index);
        assert_eq!(pack.claims().claim_count(), 6);
        let properties = [
            ("molar-mass", MolarMass::DIMS),
            ("standard-molar-formation-enthalpy", MOLAR_ENERGY),
            ("standard-molar-formation-gibbs-energy", MOLAR_ENERGY),
            ("standard-molar-entropy", MOLAR_ENTROPY),
            (
                "reported-formation-enthalpy-uncertainty-magnitude",
                MOLAR_ENERGY,
            ),
            (
                "reported-standard-entropy-uncertainty-magnitude",
                MOLAR_ENTROPY,
            ),
        ];
        for ((property, dims), expected) in properties.into_iter().zip(values) {
            let key = PropertyKey::new(property, dims);
            let answer = pack
                .claims()
                .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            assert!(
                (answer.evidence.value.value - expected).abs() <= 1e-12 * expected.abs().max(1.0)
            );
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            for wrong in [
                point(1 - index),
                point(2),
                query.clone().with("source-phase-solid", 0.0).unwrap(),
                query.clone().with("temperature", 300.0).unwrap(),
                query.clone().with(REFERENCE, 0.0).unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err()
                );
            }
        }
    }
}

#[test]
fn g1_actual_sources_close_the_goethite_and_lime_hess_cycle() {
    let packs = SLUGS.map(compile);
    let points = std::array::from_fn(point);
    let hydrolysis = MineralReaction::FerricHydrolysisGoethite.build().unwrap();
    let neutralization = MineralReaction::HydratedLimeNeutralization.build().unwrap();
    let combined = MineralReaction::FerricGoethitePrecipitationWithHydratedLime
        .build()
        .unwrap();
    let hydrolysis_heat = heat(&hydrolysis, &packs, &points);
    let neutralization_heat = heat(&neutralization, &packs, &points);
    let combined_heat = heat(&combined, &packs, &points);
    // Literal B2131 rows give an endothermic hydrolysis and exothermic lime treatment.
    assert_eq!(hydrolysis_heat.molar().value(), 58_900.0);
    assert_eq!(neutralization_heat.molar().value(), -128_500.0);
    assert_eq!(combined_heat.molar().value(), -267_700.0);
    assert_eq!(
        hydrolysis_heat.at_extent(Amount::new(2.0)).unwrap().value()
            + neutralization_heat
                .at_extent(Amount::new(3.0))
                .unwrap()
                .value(),
        combined_heat.at_extent(Amount::new(1.0)).unwrap().value(),
    );
    for heat in [&hydrolysis_heat, &neutralization_heat, &combined_heat] {
        for term in heat.terms() {
            let pack = packs
                .iter()
                .find(|pack| pack.pack_id() == term.pack_id)
                .unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
            assert!(!term.answer.receipt.source_hashes.is_empty());
        }
    }
    // The intermediate six protons and four water molecules cancel exactly.
    let mut sum = BTreeMap::<SpeciesId, i128>::new();
    for (reaction, multiplier) in [(&hydrolysis, 2), (&neutralization, 3)] {
        for term in reaction.terms() {
            *sum.entry(term.species.clone()).or_default() += multiplier * term.coefficient;
        }
    }
    sum.retain(|_, coefficient| *coefficient != 0);
    assert_eq!(
        sum,
        combined
            .terms()
            .iter()
            .map(|term| (term.species.clone(), term.coefficient))
            .collect()
    );
    assert!(
        combined
            .terms()
            .iter()
            .all(|term| term.species.as_str() != "Fe(OH)3(s)")
    );
}

#[test]
fn g1_g3_direct_gibbs_and_explicit_activities_close_the_precipitation_cycle() {
    let packs = SLUGS.map(compile);
    let points = std::array::from_fn(point);
    let reactions = [
        MineralReaction::FerricHydrolysisGoethite,
        MineralReaction::HydratedLimeNeutralization,
        MineralReaction::FerricGoethitePrecipitationWithHydratedLime,
    ]
    .map(|reaction| reaction.build().unwrap());
    let standard = std::array::from_fn::<_, 3, _>(|i| {
        reactions[i]
            .standard_gibbs_change(&gibbs_sources(&reactions[i], &packs, &points))
            .unwrap()
    });
    // Independent literal Gf columns in B2131 p16, not fitted H - T*S values.
    for (result, expected) in standard.iter().zip([-900.0, -129_800.0, -391_200.0]) {
        assert_eq!(result.molar().value(), expected);
        for term in result.terms() {
            let pack = packs.iter().find(|p| p.pack_id() == term.pack_id).unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
            assert!(term.answer.receipt.observation_backed);
            assert!(!term.answer.receipt.source_hashes.is_empty());
            assert_eq!(
                term.answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
    }
    assert_eq!(
        2.0 * standard[0].molar().value() + 3.0 * standard[1].molar().value(),
        standard[2].molar().value()
    );

    // Shared declared activities make the proton cancellation observable in
    // both the reaction quotient and the local Gibbs energy. Solids and water
    // are explicitly present; the evaluator does not invent unit activities.
    let activities = reactions.each_ref().map(|reaction| {
        reaction
            .terms()
            .iter()
            .map(|term| {
                let activity = match term.species.as_str() {
                    "Fe[3+](aq)" | "Ca[2+](aq)" => 1e-3,
                    _ => 1.0,
                };
                (term.species.clone(), Dimensionless::new(activity))
            })
            .collect::<Vec<_>>()
    });
    let states =
        std::array::from_fn::<_, 3, _>(|i| standard[i].at_activities(&activities[i]).unwrap());
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    for (i, log_q) in [1e3_f64.ln(), 1e-3_f64.ln(), 1e-3_f64.ln()]
        .into_iter()
        .enumerate()
    {
        assert!((states[i].log_quotient.value() - log_q).abs() < 1e-12);
        assert!(
            (states[i].molar.value() - (standard[i].molar().value() + rt * log_q)).abs() < 1e-8
        );
        assert!(std::ptr::eq(states[i].standard, &standard[i]));
        assert_eq!(states[i].activities, activities[i]);
        assert!(standard[i].at_activities(&activities[i][1..]).is_err());
    }
    assert!(
        (2.0 * states[0].molar.value() + 3.0 * states[1].molar.value() - states[2].molar.value())
            .abs()
            < 1e-8
    );
    assert!(states[0].molar.value() > 0.0);
    let mut less_acidic = activities[0].clone();
    let proton = less_acidic
        .iter_mut()
        .find(|(species, _)| species.as_str() == "H+(aq)")
        .unwrap();
    proton.1 = Dimensionless::new(1e-3);
    let lower_proton_activity = standard[0].at_activities(&less_acidic).unwrap();
    assert!((lower_proton_activity.log_quotient.value() - 1e-6_f64.ln()).abs() < 1e-12);
    assert!(lower_proton_activity.molar.value() < 0.0);

    for (reaction_index, solid_index, phase) in [
        (0, 0, "source-phase-goethite"),
        (2, 1, "source-phase-portlandite"),
    ] {
        let reaction = &reactions[reaction_index];
        let mut sources = gibbs_sources(reaction, &packs, &points);
        assert!(reaction.standard_gibbs_change(&sources[1..]).is_err());
        let wrong_phase = points[solid_index].clone().with(phase, 0.0).unwrap();
        sources
            .iter_mut()
            .find(|source| source.species.as_str() == SPECIES[solid_index])
            .unwrap()
            .at = &wrong_phase;
        assert!(reaction.standard_gibbs_change(&sources).is_err());
    }
}

#[test]
fn g1_g3_finite_lime_bounds_precipitation_and_hydrolysis_protons_fund_only_later_stages() {
    let combined = MineralReaction::FerricGoethitePrecipitationWithHydratedLime
        .build()
        .unwrap();
    let initial = combined
        .terms()
        .iter()
        .map(|term| {
            let amount = match term.species.as_str() {
                "Fe[3+](aq)" => 10.0,
                "Ca(OH)2(s)" => 6.0,
                _ => 0.0,
            };
            (term.species.clone(), Amount::new(amount))
        })
        .collect::<Vec<_>>();
    let reactants = initial
        .iter()
        .filter(|(id, _)| matches!(id.as_str(), "Fe[3+](aq)" | "Ca(OH)2(s)"))
        .cloned()
        .collect::<Vec<_>>();
    assert_eq!(combined.limiting_extent(&reactants).unwrap().value(), 2.0);
    assert!(matches!(
        combined.advance_inventory(&initial, Amount::new(2.1)),
        Err(ReactionBalanceError::ExtentExceedsInventory { .. })
    ));
    let final_inventory = combined
        .advance_inventory(&initial, Amount::new(2.0))
        .unwrap();
    for (name, expected) in [
        ("Fe[3+](aq)", 6.0),
        ("Ca(OH)2(s)", 0.0),
        ("FeO(OH)(s)", 4.0),
        ("Ca[2+](aq)", 6.0),
        ("H2O(l)", 4.0),
    ] {
        assert_eq!(
            final_inventory
                .iter()
                .find(|(id, _)| id.as_str() == name)
                .unwrap()
                .1
                .value(),
            expected
        );
    }

    let hydrolysis = MineralReaction::FerricHydrolysisGoethite.build().unwrap();
    let neutralization = MineralReaction::HydratedLimeNeutralization.build().unwrap();
    let network = ReactionNetwork::new(vec![hydrolysis.clone(), neutralization.clone()]).unwrap();
    let staged_initial = network
        .species()
        .iter()
        .map(|id| {
            let amount = match id.as_str() {
                "Fe[3+](aq)" => 2.0,
                "Ca(OH)2(s)" => 3.0,
                "H2O(l)" => 4.0,
                _ => 0.0,
            };
            (id.clone(), Amount::new(amount))
        })
        .collect::<Vec<_>>();
    let request = |hydrolysis_extent, neutralization_extent| {
        [
            (hydrolysis.id().clone(), Amount::new(hydrolysis_extent)),
            (
                neutralization.id().clone(),
                Amount::new(neutralization_extent),
            ),
        ]
    };
    assert!(matches!(
        network.advance_inventory(&staged_initial, &request(2.0, 3.0)),
        Err(ReactionNetworkError::InsufficientInventory { .. })
    ));
    let hydrolyzed = network
        .advance_inventory(&staged_initial, &request(2.0, 0.0))
        .unwrap();
    let treated = network
        .advance_inventory(&hydrolyzed.final_inventory, &request(0.0, 3.0))
        .unwrap();
    // Four initial water moles remain part of the inventory: the net reaction adds two.
    for (name, expected) in SPECIES.into_iter().zip([2.0, 0.0, 0.0, 6.0, 0.0, 3.0]) {
        assert_eq!(
            treated
                .final_inventory
                .iter()
                .find(|(id, _)| id.as_str() == name)
                .unwrap()
                .1
                .value(),
            expected
        );
    }
}
