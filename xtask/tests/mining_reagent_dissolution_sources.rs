//! Actual NBS reagent sources through conserved dissolution and shared pH.
//! Dissolved extents are explicit inputs; no dissolution completion is inferred.
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
        FormationEnthalpySource,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidBaseSpeciation, AcidDissociation, AcidFamily, IdealAcidBaseSystem,
                SpeciationOptions, SpectatorIon, WaterAutoprotolysis,
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const CUBIC: &str = "source-phase-sodium-cyanide-cubic";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const FAMILY: &str = "uncomplexed-cyanide";
const PACKS: [&str; 8] = [
    "sodium-hydroxide-solid-nbs1989-thermochemistry",
    "sodium-cyanide-cubic-solid-nbs1989-thermochemistry",
    "sodium-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "cyanide-ion-aqueous-nbs1989",
    "hydrogen-cyanide-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
];
const SPECIES: [&str; 8] = [
    "NaOH(s)",
    "NaCN(s)",
    "Na[+](aq)",
    "OH[-](aq)",
    "CN[-](aq)",
    "HCN(aq)",
    "H+(aq)",
    "H2O(l)",
];

struct Materials {
    packs: Vec<NormalizedPack>,
    points: Vec<QueryPoint>,
    ids: Vec<SpeciesId>,
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
    if index <= 1 {
        flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]);
        if index == 1 {
            flags.push((CUBIC, 1.0));
        }
    } else if index == 7 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if index != 5 {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}

fn materials() -> &'static Materials {
    static MATERIALS: OnceLock<Materials> = OnceLock::new();
    MATERIALS.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-reagent-dissolution-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        let packs = PACKS
            .iter()
            .map(|name| {
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
            points: (0..8).map(point).collect(),
            ids: SPECIES
                .iter()
                .map(|s| SpeciesId::new(*s).unwrap())
                .collect(),
        }
    })
}

fn value(index: usize, name: &str, dims: Dims) -> f64 {
    let m = materials();
    let answer = m.packs[index]
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            &m.points[index],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    m.packs[index]
        .claims()
        .verify_receipt(&answer.receipt)
        .unwrap();
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    answer.evidence.value.value
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn species_index(species: &SpeciesId) -> usize {
    SPECIES.iter().position(|s| *s == species.as_str()).unwrap()
}
fn stock(inventory: &[(SpeciesId, Amount)], species: &str) -> f64 {
    inventory
        .iter()
        .find(|(s, _)| s.as_str() == species)
        .unwrap()
        .1
        .value()
}

#[test]
fn original_reagent_rows_retain_phase_reference_units_and_uncertainty() {
    // Independent transcription of Table 99:Na, Sodium (Prepared 1980),
    // retaining NaCN cr I, cubic rather than the adjacent orthorhombic H row.
    let m = materials();
    for (i, literals) in [
        (0, [0.0399972, -425_609.0, -379_494.0, 64.455]),
        (1, [0.0490077, -87_490.0, -76_430.0, 115.60]),
        (2, [0.0229898, -240_120.0, -261_905.0, 59.0]),
    ] {
        assert_eq!(m.packs[i].claims().claim_count(), 4);
        let entropy = if i == 2 {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        for (name, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, literals[0]),
            ("standard-molar-formation-enthalpy", ENERGY, literals[1]),
            ("standard-molar-formation-gibbs-energy", ENERGY, literals[2]),
            (entropy, ENTROPY, literals[3]),
        ] {
            close(value(i, name, dims), expected);
            let key = PropertyKey::new(name, dims);
            for wrong in [
                point(if i == 2 { 0 } else { 2 }),
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
                    .with_quantity(
                        "pressure",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::Pressure,
                            ValueForm::Static,
                        )),
                        101_325.0,
                    )
                    .unwrap(),
                point(i)
                    .with_quantity(REFERENCE, QuantitySpec::dimensional(Dims::NONE), 0.0)
                    .unwrap(),
                point(i)
                    .with_quantity(
                        "source-specimen-qualified",
                        QuantitySpec::dimensional(Dims::NONE),
                        1.0,
                    )
                    .unwrap(),
            ] {
                assert!(
                    m.packs[i]
                        .claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err()
                );
            }
            if i == 1 {
                for wrong in [
                    point(0), // Generic solid support does not establish cubic NaCN.
                    point(1)
                        .with_quantity(CUBIC, QuantitySpec::dimensional(Dims::NONE), 0.0)
                        .unwrap(),
                ] {
                    assert!(
                        m.packs[i]
                            .claims()
                            .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                            .is_err()
                    );
                }
            }
        }
        let wrong_entropy = if i == 2 {
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        assert!(
            m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new(wrong_entropy, ENTROPY),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn actual_sources_close_dissolution_energies_mass_and_finite_reagent_accounts() {
    let m = materials();
    for (recipe, solid, anion, expected_h, expected_g, expected_s) in [
        (
            MineralReaction::SodiumHydroxideDissolution,
            0,
            3,
            -44_505.0,
            -39_655.0,
            -16.205,
        ),
        (
            MineralReaction::SodiumCyanideDissolution,
            1,
            4,
            -2030.0,
            -13_075.0,
            37.5,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let enthalpies = reaction
            .terms()
            .iter()
            .map(|term| {
                let i = species_index(&term.species);
                FormationEnthalpySource {
                    species: &term.species,
                    pack: &m.packs[i],
                    at: &m.points[i],
                    reference_axis: REFERENCE,
                }
            })
            .collect::<Vec<_>>();
        let gibbs = reaction
            .terms()
            .iter()
            .map(|term| {
                let i = species_index(&term.species);
                FormationGibbsSource {
                    species: &term.species,
                    pack: &m.packs[i],
                    at: &m.points[i],
                    reference_axis: REFERENCE,
                }
            })
            .collect::<Vec<_>>();
        let caloric = reaction.standard_equilibrium(&enthalpies).unwrap();
        let direct = reaction.standard_gibbs_change(&gibbs).unwrap();
        close(caloric.enthalpy().molar().value(), expected_h);
        close(caloric.molar_entropy().value(), expected_s);
        close(
            caloric.molar_gibbs().value(),
            expected_h - 298.15 * expected_s,
        );
        close(direct.molar().value(), expected_g);
        assert_eq!(direct.reference_axis(), REFERENCE);
        for term in direct.terms() {
            m.packs[species_index(&term.term.species)]
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        let molar_mass = value(solid, "molar-mass", MolarMass::DIMS);
        close(
            molar_mass,
            value(2, "molar-mass", MolarMass::DIMS) + value(anion, "molar-mass", MolarMass::DIMS),
        );
        let extent = reaction
            .extent_from_reacted_mass(
                &m.ids[solid],
                Mass::new(0.1 * molar_mass),
                MolarMass::new(molar_mass),
            )
            .unwrap();
        close(extent.value(), 0.1);
        let inventory = reaction
            .terms()
            .iter()
            .map(|term| {
                (
                    term.species.clone(),
                    Amount::new(if term.species == m.ids[solid] {
                        0.25
                    } else {
                        0.0
                    }),
                )
            })
            .collect::<Vec<_>>();
        let dissolved = reaction.advance_inventory(&inventory, extent).unwrap();
        close(stock(&dissolved, SPECIES[solid]), 0.15);
        close(stock(&dissolved, SPECIES[2]), 0.1);
        close(stock(&dissolved, SPECIES[anion]), 0.1);
        assert!(
            reaction
                .advance_inventory(&inventory, Amount::new(0.3))
                .is_err()
        );
    }
}

fn gibbs_sources(indices: &[usize]) -> Vec<FormationGibbsSource<'static>> {
    let m = materials();
    indices
        .iter()
        .map(|&i| FormationGibbsSource {
            species: &m.ids[i],
            pack: &m.packs[i],
            at: &m.points[i],
            reference_axis: REFERENCE,
        })
        .collect()
}

fn acid_base_system() -> IdealAcidBaseSystem {
    let ids = &materials().ids;
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build().unwrap(),
        ids[7].clone(),
        ids[6].clone(),
        ids[3].clone(),
        &gibbs_sources(&[7, 6, 3]),
    )
    .unwrap();
    let hcn = AcidDissociation::from_sources(
        MineralReaction::HydrogenCyanideDissociation
            .build()
            .unwrap(),
        ids[5].clone(),
        ids[4].clone(),
        ids[6].clone(),
        &gibbs_sources(&[5, 6, 4]),
    )
    .unwrap();
    IdealAcidBaseSystem::new(water, vec![AcidFamily::new(FAMILY, vec![hcn]).unwrap()]).unwrap()
}

fn declared_dissolved_mass(
    recipe: MineralReaction,
    solid: usize,
    kg: f64,
) -> Vec<(SpeciesId, Amount)> {
    let m = materials();
    let reaction = recipe.build().unwrap();
    let extent = reaction
        .extent_from_reacted_mass(
            &m.ids[solid],
            Mass::new(kg),
            MolarMass::new(value(solid, "molar-mass", MolarMass::DIMS)),
        )
        .unwrap();
    let initial = reaction
        .terms()
        .iter()
        .map(|term| {
            (
                term.species.clone(),
                if term.species == m.ids[solid] {
                    extent
                } else {
                    Amount::new(0.0)
                },
            )
        })
        .collect::<Vec<_>>();
    reaction.advance_inventory(&initial, extent).unwrap()
}

fn form(answer: &AcidBaseSpeciation, name: &str) -> f64 {
    answer.families[0]
        .forms
        .iter()
        .find(|s| s.species.as_str() == name)
        .unwrap()
        .amount
        .value()
}

#[test]
fn declared_dissolved_reagent_masses_feed_shared_ph_with_their_counterions() {
    // These masses are declared actually dissolved, not assumed to dissolve
    // merely because the material packs exist. Both are pure constituents.
    let salt = declared_dissolved_mass(MineralReaction::SodiumCyanideDissolution, 1, 0.00490077);
    let caustic =
        declared_dissolved_mass(MineralReaction::SodiumHydroxideDissolution, 0, 0.000799944);
    close(stock(&salt, "CN[-](aq)"), 0.1);
    close(stock(&salt, "Na[+](aq)"), 0.1);
    close(stock(&caustic, "OH[-](aq)"), 0.02);
    close(stock(&caustic, "Na[+](aq)"), 0.02);
    let system = acid_base_system();
    let total = stock(&salt, "CN[-](aq)");
    let solve = |sodium| {
        system
            .speciate(
                &[(FAMILY.to_owned(), Amount::new(total))],
                &[SpectatorIon {
                    species: materials().ids[2].clone(),
                    charge: 1,
                    amount: Amount::new(sodium),
                }],
                // Explicit 10 kg solvent basis for NBS unit-molal activity.
                Amount::new(10.0),
                SpeciationOptions {
                    ph_tolerance: Dimensionless::new(1e-11),
                    absolute_charge_tolerance: Amount::new(1e-14),
                    relative_charge_tolerance: Dimensionless::new(1e-12),
                    ..SpeciationOptions::default()
                },
                || false,
            )
            .unwrap()
    };
    let without_caustic = solve(stock(&salt, "Na[+](aq)"));
    let with_caustic = solve(stock(&salt, "Na[+](aq)") + stock(&caustic, "Na[+](aq)"));
    assert!(with_caustic.ph.value() > without_caustic.ph.value());
    assert!(form(&with_caustic, "CN[-](aq)") > form(&without_caustic, "CN[-](aq)"));
    // After dissolution, OH is solved by the shared water/protonation owner;
    // appending it as another fixed spectator would double-count the base.
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    let ka = (-52_700.0 / rt).exp();
    let kw = (-79_885.0 / rt).exp();
    for answer in [&without_caustic, &with_caustic] {
        let cyanide = form(answer, "CN[-](aq)");
        let hcn = form(answer, "HCN(aq)");
        close(cyanide + hcn, total);
        close(answer.charge_balance_residual.value(), 0.0);
        assert!(answer.max_log_mass_action_residual.value() < 1e-10);
        assert!((answer.proton_activity.value() * cyanide / hcn / ka - 1.0).abs() < 1e-9);
        assert!(
            (answer.proton_activity.value() * answer.hydroxide_activity.value() / kw - 1.0).abs()
                < 1e-9
        );
        assert_eq!(answer.spectators.len(), 1);
        assert_eq!(answer.spectators[0].species.as_str(), "Na[+](aq)");
    }
}
