//! G0/G1 source and finite-stock checks for ferrite-bearing zinc calcine.
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
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SLUGS: [&str; 8] = [
    "franklinite-solid-usgs-b2131-thermochemistry",
    "hematite-solid-usgs-b2131-thermochemistry",
    "zincite-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
    "zinc-ion-aqueous-usgs-b2131",
    "iron-iii-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
    "goethite-solid-usgs-b2131-thermochemistry",
];
const SPECIES: [&str; 8] = [
    "ZnFe2O4(s)",
    "Fe2O3(s)",
    "ZnO(s)",
    "H+(aq)",
    "Zn[2+](aq)",
    "Fe[3+](aq)",
    "H2O(l)",
    "FeO(OH)(s)",
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
        "frankensim-zinc-ferrite-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed),
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
        3..=5 => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
        6 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]),
    }
    match index {
        0 => flags.push(("source-phase-franklinite", 1.0)),
        1 => flags.push(("source-phase-hematite", 1.0)),
        2 => flags.push(("source-phase-zincite", 1.0)),
        7 => flags.push(("source-phase-goethite", 1.0)),
        _ => {}
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn index(species: &SpeciesId) -> usize {
    SPECIES
        .iter()
        .position(|name| *name == species.as_str())
        .unwrap()
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
fn new_phase_sources_retain_literal_values_and_refuse_other_material_states() {
    // Independent transcription of USGS B2131 printed p22 (ZnFe2O4) and
    // p16 (hematite); in particular these are not JANAF hematite values.
    for (i, values, magnitudes) in [
        (
            0,
            [0.241082, -1_188_100.0, -1_082_100.0, 150.7],
            [4000.0, 4000.0, 0.3],
        ),
        (
            1,
            [0.159692, -826_200.0, -744_400.0, 87.4],
            [1300.0, 1300.0, 0.2],
        ),
    ] {
        let pack = compile(SLUGS[i]);
        let query = point(i);
        assert_eq!(pack.claims().claim_count(), 7);
        for (property, dims, expected) in [
            ("molar-mass", MolarMass::DIMS, values[0]),
            ("standard-molar-formation-enthalpy", ENERGY, values[1]),
            ("standard-molar-formation-gibbs-energy", ENERGY, values[2]),
            ("standard-molar-entropy", ENTROPY, values[3]),
            (
                "reported-formation-enthalpy-uncertainty-magnitude",
                ENERGY,
                magnitudes[0],
            ),
            (
                "reported-formation-gibbs-energy-uncertainty-magnitude",
                ENERGY,
                magnitudes[1],
            ),
            (
                "reported-standard-entropy-uncertainty-magnitude",
                ENTROPY,
                magnitudes[2],
            ),
        ] {
            let key = PropertyKey::new(property, dims);
            let answer = pack
                .claims()
                .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            close(answer.evidence.value.value, expected);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            for wrong in [
                point(1 - i),
                query
                    .clone()
                    .with_quantity(
                        "temperature",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        900.0,
                    )
                    .unwrap(),
                query
                    .clone()
                    .with_quantity(
                        "pressure",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::Pressure,
                            ValueForm::Static,
                        )),
                        101_325.0,
                    )
                    .unwrap(),
                query
                    .clone()
                    .with_quantity(REFERENCE, QuantitySpec::dimensional(Dims::NONE), 0.0)
                    .unwrap(),
                query
                    .clone()
                    .with_quantity(
                        "source-specimen-qualified",
                        QuantitySpec::dimensional(Dims::NONE),
                        1.0,
                    )
                    .unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly,)
                        .is_err()
                );
            }
        }
    }
}

#[test]
fn actual_ferrite_sources_close_dissolution_hydrolysis_and_net_hess_accounts() {
    let packs = SLUGS.map(compile);
    let points: [QueryPoint; 8] = std::array::from_fn(point);
    let mut heats = Vec::new();
    let mut energies = Vec::new();
    for (recipe, expected_heat, expected_gibbs, expected_entropy) in [
        (
            MineralReaction::ZincFerriteFormation,
            -11_400.0,
            -17_300.0,
            20.1,
        ),
        (
            MineralReaction::ZincFerriteAcidLeach,
            -208_300.0,
            -47_000.0,
            -540.5,
        ),
        (
            MineralReaction::FerricHydrolysisGoethite,
            58_900.0,
            -900.0,
            200.4,
        ),
        (
            MineralReaction::ZincFerriteToGoethiteLeach,
            -90_500.0,
            -48_800.0,
            -139.7,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let enthalpy_sources = reaction
            .terms()
            .iter()
            .map(|term| FormationEnthalpySource {
                species: &term.species,
                pack: &packs[index(&term.species)],
                at: &points[index(&term.species)],
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>();
        let gibbs_sources = reaction
            .terms()
            .iter()
            .map(|term| FormationGibbsSource {
                species: &term.species,
                pack: &packs[index(&term.species)],
                at: &points[index(&term.species)],
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>();
        let caloric = reaction.standard_equilibrium(&enthalpy_sources).unwrap();
        let gibbs = reaction.standard_gibbs_change(&gibbs_sources).unwrap();
        close(caloric.enthalpy().molar().value(), expected_heat);
        close(caloric.molar_entropy().value(), expected_entropy);
        close(
            caloric.molar_gibbs().value(),
            expected_heat - 298.15 * expected_entropy,
        );
        close(gibbs.molar().value(), expected_gibbs);
        for term in gibbs.terms() {
            let pack = packs
                .iter()
                .find(|pack| pack.pack_id() == term.pack_id)
                .unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        }
        heats.push(caloric.enthalpy().molar().value());
        energies.push(gibbs.molar().value());
    }
    // Dissolution followed by two hydrolysis extents regenerates six protons;
    // the resulting net endpoint still consumes two per zinc.
    close(heats[1] + 2.0 * heats[2], heats[3]);
    close(energies[1] + 2.0 * energies[2], energies[3]);
}

#[test]
fn ferrite_and_acid_stocks_are_admitted_before_later_iron_precipitation() {
    let formation = MineralReaction::ZincFerriteFormation.build().unwrap();
    let initial = formation
        .terms()
        .iter()
        .map(|term| {
            (
                term.species.clone(),
                Amount::new(match term.species.as_str() {
                    "ZnO(s)" => 1.0,
                    "Fe2O3(s)" => 0.5,
                    _ => 0.0,
                }),
            )
        })
        .collect::<Vec<_>>();
    let partial = formation
        .advance_inventory(&initial, Amount::new(0.5))
        .unwrap();
    close(stock(&partial, "ZnFe2O4(s)"), 0.5);
    close(stock(&partial, "ZnO(s)"), 0.5);
    assert!(
        formation
            .advance_inventory(&initial, Amount::new(1.0))
            .is_err()
    );

    let dissolve = MineralReaction::ZincFerriteAcidLeach.build().unwrap();
    let hydrolyze = MineralReaction::FerricHydrolysisGoethite.build().unwrap();
    let network = ReactionNetwork::new(vec![dissolve.clone(), hydrolyze.clone()]).unwrap();
    let inventory = |acid| {
        network
            .species()
            .iter()
            .map(|species| {
                (
                    species.clone(),
                    Amount::new(match species.as_str() {
                        "ZnFe2O4(s)" => 1.0,
                        "H+(aq)" => acid,
                        _ => 0.0,
                    }),
                )
            })
            .collect::<Vec<_>>()
    };
    let extents = |leach, precipitation| {
        [
            (dissolve.id().clone(), Amount::new(leach)),
            (hydrolyze.id().clone(), Amount::new(precipitation)),
        ]
    };
    assert!(matches!(
        network.advance_inventory(&inventory(2.0), &extents(1.0, 2.0)),
        Err(ReactionNetworkError::InsufficientInventory { .. })
    ));
    assert!(
        network
            .advance_inventory(&inventory(2.0), &extents(1.0, 0.0))
            .is_err()
    );
    let leached = network
        .advance_inventory(&inventory(8.0), &extents(1.0, 0.0))
        .unwrap();
    close(stock(&leached.final_inventory, "H+(aq)"), 0.0);
    close(stock(&leached.final_inventory, "Fe[3+](aq)"), 2.0);
    close(stock(&leached.final_inventory, "H2O(l)"), 4.0);
    let finished = network
        .advance_inventory(&leached.final_inventory, &extents(0.0, 2.0))
        .unwrap();
    for (species, value) in [
        ("ZnFe2O4(s)", 0.0),
        ("H+(aq)", 6.0),
        ("Zn[2+](aq)", 1.0),
        ("Fe[3+](aq)", 0.0),
        ("H2O(l)", 0.0),
        ("FeO(OH)(s)", 2.0),
    ] {
        close(stock(&finished.final_inventory, species), value);
    }

    // The declared aggregate equation reproduces the final account, while
    // leaving the actual staged admission above intact.
    let aggregate: ConservedReaction = MineralReaction::ZincFerriteToGoethiteLeach.build().unwrap();
    let combined_initial = aggregate
        .terms()
        .iter()
        .map(|term| {
            (
                term.species.clone(),
                Amount::new(match term.species.as_str() {
                    "ZnFe2O4(s)" => 1.0,
                    "H+(aq)" => 8.0,
                    _ => 0.0,
                }),
            )
        })
        .collect::<Vec<_>>();
    let combined = aggregate
        .advance_inventory(&combined_initial, Amount::new(1.0))
        .unwrap();
    for (species, amount) in combined {
        close(
            amount.value(),
            stock(&finished.final_inventory, species.as_str()),
        );
    }
}
