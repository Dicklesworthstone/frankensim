//! G0/G1/G3: the iron redox stages retain stocks and close the source heat cycle.
//! Equations 3-4: https://pubs.usgs.gov/of/1999/0348/report.pdf.
#![deny(unsafe_code)]

use std::collections::BTreeMap;
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
        ConservedReaction, FormationEnthalpySource, ReactionEnthalpy,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SLUGS: [&str; 7] = [
    "pyrite-solid-usgs-b2131-thermochemistry",
    "oxygen-gas-usgs-b2131",
    "water-liquid-usgs-b2131",
    "iron-ii-ion-aqueous-usgs-b2131",
    "sulfate-ion-aqueous-usgs-b2131",
    "hydrogen-ion-aqueous-usgs-b2131",
    "iron-iii-ion-aqueous-usgs-b2131",
];
const SPECIES: [&str; 7] = [
    "FeS2(s)",
    "O2(g)",
    "H2O(l)",
    "Fe[2+](aq)",
    "SO4[2-](aq)",
    "H+(aq)",
    "Fe[3+](aq)",
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
        "frankensim-iron-redox-{}-{nonce}-{}",
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
        0 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-pyrite", 1.0),
            ("source-pure-species", 1.0),
        ]),
        1 => flags.extend([("source-phase-gas", 1.0), ("source-pure-species", 1.0)]),
        2 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
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
    packs: &[NormalizedPack; 7],
    points: &[QueryPoint; 7],
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

#[test]
fn g0_g3_ferric_source_keeps_its_conventional_entropy_and_uncertainty_magnitudes() {
    let pack = compile(SLUGS[6]);
    let query = point(6);
    // Keep the directly tabulated Gibbs claim alongside Hf/S and uncertainty metadata.
    assert_eq!(pack.claims().claim_count(), 6);
    for (property, dims, expected) in [
        ("molar-mass", MolarMass::DIMS, 0.055847),
        ("standard-molar-formation-enthalpy", MOLAR_ENERGY, -49_900.0),
        ("standard-molar-formation-gibbs-energy", MOLAR_ENERGY, -16_700.0),
        ("conventional-standard-molar-entropy", MOLAR_ENTROPY, -280.0),
        (
            "reported-formation-enthalpy-uncertainty-magnitude",
            MOLAR_ENERGY,
            5000.0,
        ),
        (
            "reported-conventional-entropy-uncertainty-magnitude",
            MOLAR_ENTROPY,
            13.0,
        ),
    ] {
        let key = PropertyKey::new(property, dims);
        let answer = pack
            .claims()
            .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert!((answer.evidence.value.value - expected).abs() <= 1e-12 * expected.abs().max(1.0));
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        assert!(answer.receipt.observation_backed);
        pack.claims().verify_receipt(&answer.receipt).unwrap();
        for wrong in [
            point(0),
            query.clone().with("temperature", 300.0).unwrap(),
            query.clone().with(REFERENCE, 0.0).unwrap(),
            query
                .clone()
                .with("source-ionic-reference-hydrogen-zero", 0.0)
                .unwrap(),
            query.clone().with("source-individual-ion", 0.0).unwrap(),
        ] {
            assert!(
                pack.claims()
                    .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
    }
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::new("standard-molar-entropy", MOLAR_ENTROPY),
                &query,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}

#[test]
fn g1_real_sources_and_exact_coefficients_close_the_iron_redox_hess_cycle() {
    let packs = SLUGS.map(compile);
    let points = std::array::from_fn(point);
    let oxidation = MineralReaction::FerrousOxidationFerric.build().unwrap();
    let attack = MineralReaction::FerricPyriteOxidation.build().unwrap();
    let net = MineralReaction::PyriteOxidationFerrous.build().unwrap();
    let oxidizing_heat = heat(&oxidation, &packs, &points);
    let attacking_heat = heat(&attack, &packs, &points);
    let net_heat = heat(&net, &packs, &points);
    // Literal USGS rows: Fe2+ -91.1, Fe3+ -49.9, H2O -285.8,
    // pyrite -171.5, SO4(2-) -909.3 kJ/mol. O2 and H+ are explicit zeros.
    assert_eq!(oxidizing_heat.molar().value(), -406_800.0);
    assert_eq!(attacking_heat.molar().value(), -28_600.0);
    assert_eq!(net_heat.molar().value(), -2_904_800.0);
    assert_eq!(
        oxidizing_heat.at_extent(Amount::new(7.0)).unwrap().value()
            + attacking_heat.at_extent(Amount::new(2.0)).unwrap().value(),
        net_heat.at_extent(Amount::new(1.0)).unwrap().value(),
    );
    for heat in [&oxidizing_heat, &attacking_heat, &net_heat] {
        for term in heat.terms() {
            let pack = packs
                .iter()
                .find(|pack| pack.pack_id() == term.pack_id)
                .unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
            assert!(!term.answer.receipt.source_hashes.is_empty());
        }
    }
    // Check exact stoichiometric cancellation independently of the heat sum:
    // 7 integer ferrous oxidations + 2 ferric attacks = 1 two-pyrite endpoint.
    let mut combined = BTreeMap::<SpeciesId, i128>::new();
    for (reaction, multiplier) in [(&oxidation, 7), (&attack, 2)] {
        for term in reaction.terms() {
            *combined.entry(term.species.clone()).or_default() += multiplier * term.coefficient;
        }
    }
    combined.retain(|_, coefficient| *coefficient != 0);
    assert_eq!(
        combined,
        net.terms()
            .iter()
            .map(|term| (term.species.clone(), term.coefficient))
            .collect()
    );
}

#[test]
fn g1_g3_explicit_stages_restore_ferric_inventory_without_borrowing_future_products() {
    let oxidation = MineralReaction::FerrousOxidationFerric.build().unwrap();
    let attack = MineralReaction::FerricPyriteOxidation.build().unwrap();
    let network = ReactionNetwork::new(vec![oxidation.clone(), attack.clone()]).unwrap();
    let initial = network
        .species()
        .iter()
        .map(|id| {
            let amount = match id.as_str() {
                "FeS2(s)" => 2.0,
                "O2(g)" => 7.0,
                "H2O(l)" => 16.0,
                "Fe[3+](aq)" => 28.0,
                _ => 0.0,
            };
            (id.clone(), Amount::new(amount))
        })
        .collect::<Vec<_>>();
    let request = |oxidation_extent, attack_extent| {
        [
            (oxidation.id().clone(), Amount::new(oxidation_extent)),
            (attack.id().clone(), Amount::new(attack_extent)),
        ]
    };
    // Fe2+ and H+ are initially absent. Ferric attack must produce them in an
    // earlier declared stage; a net-zero algebraic cycle cannot grant stocks.
    assert!(matches!(
        network.advance_inventory(&initial, &request(7.0, 2.0)),
        Err(ReactionNetworkError::InsufficientInventory { .. })
    ));
    let first = network
        .advance_inventory(&initial, &request(0.0, 2.0))
        .unwrap();
    let final_state = network
        .advance_inventory(&first.final_inventory, &request(7.0, 0.0))
        .unwrap();
    let final_values = [0.0, 0.0, 14.0, 2.0, 4.0, 4.0, 28.0];
    for (name, expected) in SPECIES.into_iter().zip(final_values) {
        let actual = final_state
            .final_inventory
            .iter()
            .find(|(id, _)| id.as_str() == name)
            .unwrap()
            .1;
        assert_eq!(actual.value(), expected, "{name}");
    }
    let ferric_limited = [
        (SpeciesId::new("FeS2(s)").unwrap(), Amount::new(10.0)),
        (SpeciesId::new("Fe[3+](aq)").unwrap(), Amount::new(28.0)),
        (SpeciesId::new("H2O(l)").unwrap(), Amount::new(100.0)),
    ];
    assert_eq!(
        attack.limiting_extent(&ferric_limited).unwrap().value(),
        2.0
    );
}
