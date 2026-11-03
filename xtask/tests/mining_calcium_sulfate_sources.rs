//! Actual NBS calcium-sulfate phases through conserved reactions and finite stocks.
//! G0 source/stoichiometry checks and G1 independent ideal-solubility answers.
//! Authored inventories are numerical scenarios, not measured treatment results.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, Mass, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        network::ReactionNetwork,
        thermodynamics::{
            FormationGibbsSource,
            equilibrium::{
                EquilibriumExtentStatus, EquilibriumInventory, EquilibriumSolveOptions,
                IdealActivityModel, SingleReactionEquilibrium,
            },
        },
    },
};
use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const PACKS: [(&str, &str); 7] = [
    (
        "gypsum-selenite-solid-nbs1989-thermochemistry",
        "CaSO4.2H2O(s)",
    ),
    ("anhydrite-solid-nbs1989-thermochemistry", "CaSO4(s)"),
    ("calcium-ion-aqueous-nbs1989", "Ca[2+](aq)"),
    ("sulfate-ion-aqueous-nbs1989", "SO4[2-](aq)"),
    ("water-liquid-nbs1989-thermochemistry", "H2O(l)"),
    ("hydrogen-ion-aqueous-nbs1989", "H+(aq)"),
    (
        "calcium-hydroxide-solid-nbs1989-thermochemistry",
        "Ca(OH)2(s)",
    ),
];

fn point(index: usize) -> QueryPoint {
    let mut at = QueryPoint::new()
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
        0 | 1 | 6 => flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]),
        4 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
    }
    if index == 0 {
        flags.extend([
            ("source-phase-gypsum", 1.0),
            ("source-variety-selenite", 1.0),
        ]);
    } else if index == 1 {
        flags.push(("source-phase-anhydrite", 1.0));
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
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
struct Materials {
    packs: Vec<NormalizedPack>,
    points: Vec<QueryPoint>,
    ids: Vec<SpeciesId>,
}
fn materials() -> &'static Materials {
    static DATA: OnceLock<Materials> = OnceLock::new();
    DATA.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-calcium-sulfate-{}-{nonce}",
            std::process::id(),
        ));
        std::fs::create_dir(&out).unwrap();
        let packs = PACKS
            .iter()
            .map(|(name, _)| {
                let target = out.join(format!("{name}.fsmatpk"));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{name}/manifest.tsv")))
                    .arg("--out")
                    .arg(&target)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{name}: {}",
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(target).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(pack.pack_id(), *name);
                pack
            })
            .collect();
        Materials {
            packs,
            points: (0..PACKS.len()).map(point).collect(),
            ids: PACKS
                .iter()
                .map(|(_, id)| SpeciesId::new(*id).unwrap())
                .collect(),
        }
    })
}
fn index(id: &SpeciesId) -> usize {
    PACKS.iter().position(|(_, s)| *s == id.as_str()).unwrap()
}
fn sources(reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'static>> {
    let m = materials();
    reaction
        .terms()
        .iter()
        .map(|term| {
            let i = index(&term.species);
            FormationGibbsSource {
                species: &m.ids[i],
                pack: &m.packs[i],
                at: &m.points[i],
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
fn value(i: usize, name: &str, dims: Dims) -> f64 {
    let m = materials();
    let answer = m.packs[i]
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            &m.points[i],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    m.packs[i].claims().verify_receipt(&answer.receipt).unwrap();
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    answer.evidence.value.value
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} +/- {tolerance}"
    );
}
fn amount(result: &SingleReactionEquilibrium, species: &str) -> f64 {
    result
        .final_inventory
        .iter()
        .find(|(id, _)| id.as_str() == species)
        .unwrap()
        .1
        .value()
}
fn options() -> EquilibriumSolveOptions {
    EquilibriumSolveOptions {
        extent_tolerance: Amount::new(1e-12),
        log_residual_tolerance: Dimensionless::new(1e-9),
        max_iterations: 256,
    }
}
fn inventory(stocks: &[(usize, f64)], basis: f64) -> Vec<EquilibriumInventory> {
    stocks
        .iter()
        .map(|&(i, n)| EquilibriumInventory {
            species: materials().ids[i].clone(),
            amount: Amount::new(n),
            activity: if matches!(i, 0 | 1 | 4 | 6) {
                IdealActivityModel::PureCondensed
            } else {
                IdealActivityModel::IdealDilute {
                    amount_at_unit_activity: Amount::new(basis),
                }
            },
        })
        .collect()
}
fn solve(
    recipe: MineralReaction,
    stocks: &[(usize, f64)],
    basis: f64,
) -> SingleReactionEquilibrium {
    let reaction = recipe.build().unwrap();
    reaction
        .equilibrate_ideal(&sources(&reaction), &inventory(stocks, basis), options())
        .unwrap()
}

#[test]
fn literal_anhydrite_and_selenite_rows_preserve_hydration_phase_and_source_authority() {
    let m = materials();
    for (i, literal) in [
        (0, [0.1721724, -2_022_630.0, -1_797_280.0, 194.1]),
        (1, [0.1361416, -1_434_110.0, -1_321_790.0, 106.7]),
    ] {
        assert_eq!(m.packs[i].claims().claim_count(), 4);
        for (name, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, literal[0]),
            ("standard-molar-formation-enthalpy", ENERGY, literal[1]),
            ("standard-molar-formation-gibbs-energy", ENERGY, literal[2]),
            ("standard-molar-entropy", ENTROPY, literal[3]),
        ] {
            close(
                value(i, name, dims),
                expected,
                expected.abs().max(1.0) * 1e-12,
            );
        }
        let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
        // Missing the actual source phase/variety cannot be repaired by caller intent.
        for axis in if i == 0 {
            vec!["source-phase-gypsum", "source-variety-selenite"]
        } else {
            vec!["source-phase-anhydrite"]
        } {
            let mut missing = QueryPoint::new();
            for (name, &v) in m.points[i].axes() {
                if name != axis {
                    missing = missing
                        .with_quantity(name, m.points[i].axis_quantities()[name], v)
                        .unwrap();
                }
            }
            assert!(
                m.packs[i]
                    .claims()
                    .query_typed(&key, &missing, SelectionPolicy::SingleClaimOnly,)
                    .is_err()
            );
        }
        for wrong in [
            point(1 - i),
            point(i)
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    310.0,
                )
                .unwrap(),
            point(i)
                .with_quantity(REFERENCE, QuantitySpec::dimensional(Dims::NONE), 0.0)
                .unwrap(),
        ] {
            assert!(
                m.packs[i]
                    .claims()
                    .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly,)
                    .is_err()
            );
        }
    }
    close(
        value(0, "molar-mass", MolarMass::DIMS) - value(1, "molar-mass", MolarMass::DIMS),
        2.0 * value(4, "molar-mass", MolarMass::DIMS),
        1e-16,
    );
}

#[test]
fn sulfate_endpoints_close_atoms_hydrate_water_reagent_stocks_heat_and_direct_gibbs() {
    use MineralReaction::*;
    let m = materials();
    for (recipe, h, g, s) in [
        (GypsumPrecipitation, 1_130.0, -24_912.0, 87.28),
        (AnhydritePrecipitation, 17_990.0, -23_680.0, 139.7),
        (AnhydriteHydration, -16_860.0, -1_232.0, -52.42),
        (
            GypsumPrecipitationWithHydratedLime,
            -127_270.0,
            -154_260.0,
            90.61,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let enthalpy_sources = reaction
            .terms()
            .iter()
            .map(|term| {
                let i = index(&term.species);
                FormationEnthalpySource {
                    species: &m.ids[i],
                    pack: &m.packs[i],
                    at: &m.points[i],
                    reference_axis: REFERENCE,
                }
            })
            .collect::<Vec<_>>();
        let caloric = reaction.standard_equilibrium(&enthalpy_sources).unwrap();
        let direct = reaction.standard_gibbs_change(&sources(&reaction)).unwrap();
        close(caloric.enthalpy().molar().value(), h, 1e-8);
        close(caloric.molar_entropy().value(), s, 1e-10);
        close(caloric.molar_gibbs().value(), h - 298.15 * s, 1e-8);
        close(direct.molar().value(), g, 1e-8);
        let mass_residual: f64 = reaction
            .terms()
            .iter()
            .map(|term| {
                term.coefficient as f64 * value(index(&term.species), "molar-mass", MolarMass::DIMS)
            })
            .sum();
        close(mass_residual, 0.0, 2e-15);
        for term in direct.terms() {
            m.packs[index(&term.term.species)]
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        let mut mismatched = sources(&reaction);
        mismatched[0].reference_axis = "source-elemental-reference-usgs-b2131";
        assert!(reaction.standard_gibbs_change(&mismatched).is_err());
    }
    ReactionNetwork::new(vec![
        HydratedLimeNeutralization.build().unwrap(),
        GypsumPrecipitation.build().unwrap(),
        GypsumPrecipitationWithHydratedLime.build().unwrap(),
        AnhydritePrecipitation.build().unwrap(),
        AnhydriteHydration.build().unwrap(),
    ])
    .unwrap();

    let reaction = GypsumPrecipitationWithHydratedLime.build().unwrap();
    let extent = reaction
        .extent_from_reacted_mass(
            &m.ids[6],
            Mass::new(0.00370474),
            MolarMass::new(value(6, "molar-mass", MolarMass::DIMS)),
        )
        .unwrap();
    close(extent.value(), 0.05, 1e-15);
    let stocks = vec![
        (m.ids[6].clone(), Amount::new(0.1)),
        (m.ids[5].clone(), Amount::new(0.12)),
        (m.ids[3].clone(), Amount::new(0.1)),
        (m.ids[0].clone(), Amount::new(0.0)),
    ];
    close(
        reaction.limiting_extent(&stocks[..3]).unwrap().value(),
        0.06,
        1e-15,
    );
    let final_stocks = reaction.advance_inventory(&stocks, extent).unwrap();
    for (i, expected) in [(6, 0.05), (5, 0.02), (3, 0.05), (0, 0.05)] {
        close(
            final_stocks
                .iter()
                .find(|(id, _)| id == &m.ids[i])
                .unwrap()
                .1
                .value(),
            expected,
            1e-15,
        );
    }
    assert!(
        reaction
            .advance_inventory(&stocks, Amount::new(0.07))
            .is_err()
    );
}

#[test]
fn finite_saturation_common_ion_and_phase_exhaustion_match_independent_solutions() {
    use MineralReaction::*;
    // Independent source arithmetic for gypsum DISSOLUTION:
    // -553580 -744530 +2*(-237129) -(-1797280) = +24912 J/mol.
    let ksp = (-24_912.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15)).exp();
    let saturation = ksp.sqrt();
    let gypsum_cases = [
        (
            0.0,
            0.0,
            0.02,
            1.0,
            saturation,
            EquilibriumExtentStatus::Interior,
        ),
        (
            0.05,
            0.05,
            0.0,
            1.0,
            saturation,
            EquilibriumExtentStatus::Interior,
        ),
        (
            0.0,
            0.0,
            0.001,
            1.0,
            0.001,
            EquilibriumExtentStatus::ReverseInventoryLimit,
        ),
        (
            0.0,
            0.0,
            0.03,
            2.0,
            2.0 * saturation,
            EquilibriumExtentStatus::Interior,
        ),
    ];
    for (ca, sulfate, solid, basis, expected_ca, status) in gypsum_cases {
        let result = solve(
            GypsumPrecipitation,
            &[(2, ca), (3, sulfate), (4, 55.5), (0, solid)],
            basis,
        );
        assert_eq!(result.status, status);
        let calcium = amount(&result, "Ca[2+](aq)");
        let sulfur = amount(&result, "SO4[2-](aq)");
        let gypsum = amount(&result, "CaSO4.2H2O(s)");
        let water = amount(&result, "H2O(l)");
        close(calcium, expected_ca, 2e-11);
        close(calcium + gypsum, ca + solid, 2e-15);
        close(sulfur + gypsum, sulfate + solid, 2e-15);
        close(water + 2.0 * gypsum, 55.5 + 2.0 * solid, 2e-14);
        close(
            2.0 * calcium - 2.0 * sulfur,
            2.0 * ca - 2.0 * sulfate,
            2e-15,
        );
        if status == EquilibriumExtentStatus::Interior {
            assert!(result.log_quotient_over_equilibrium.unwrap().value().abs() <= 1e-9);
            close(calcium * sulfur / (basis * basis), ksp, 1e-13);
        } else {
            close(gypsum, 0.0, 1e-16);
            assert!(result.log_quotient_over_equilibrium.unwrap().value() > 0.0);
        }
    }
    // Declared excess sulfate has fixed spectator Na+ outside this one reaction.
    // Analytical common-ion answer s(s+0.05)=Ksp, evaluated without cancellation.
    let excess = 0.05;
    let expected = 2.0 * ksp / (excess + (excess * excess + 4.0 * ksp).sqrt());
    let common = solve(
        GypsumPrecipitation,
        &[(2, 0.0), (3, excess), (4, 55.5), (0, 0.02)],
        1.0,
    );
    close(amount(&common, "Ca[2+](aq)"), expected, 2e-11);
    close(amount(&common, "SO4[2-](aq)"), excess + expected, 2e-11);
    assert!(expected < saturation);
    let anhydrite = solve(
        AnhydritePrecipitation,
        &[(2, 0.0), (3, 0.0), (1, 0.02)],
        1.0,
    );
    let anhydrite_solubility =
        (-23_680.0 / (2.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15)).exp();
    close(
        amount(&anhydrite, "Ca[2+](aq)"),
        anhydrite_solubility,
        2e-11,
    );
    assert!(anhydrite_solubility > saturation);

    let hydration = solve(AnhydriteHydration, &[(1, 0.05), (4, 0.06), (0, 0.0)], 1.0);
    assert_eq!(
        hydration.status,
        EquilibriumExtentStatus::ForwardInventoryLimit
    );
    close(hydration.extent.value(), 0.03, 1e-15);
    close(amount(&hydration, "CaSO4(s)"), 0.02, 1e-15);
    close(amount(&hydration, "H2O(l)"), 0.0, 1e-15);

    let reaction = GypsumPrecipitation.build().unwrap();
    let stocks = [(2, 0.0), (3, 0.0), (4, 55.5), (0, 0.02)];
    assert!(
        reaction
            .equilibrate_ideal(&sources(&reaction), &inventory(&stocks, 0.0), options())
            .is_err()
    );
    let mut wrong_phase = inventory(&stocks, 1.0);
    wrong_phase[3].activity = IdealActivityModel::IdealDilute {
        amount_at_unit_activity: Amount::new(1.0),
    };
    assert!(
        reaction
            .equilibrate_ideal(&sources(&reaction), &wrong_phase, options())
            .is_err()
    );
    let mut exhausted = options();
    exhausted.max_iterations = 1;
    assert!(
        reaction
            .equilibrate_ideal(&sources(&reaction), &inventory(&stocks, 1.0), exhausted)
            .is_err()
    );
}
