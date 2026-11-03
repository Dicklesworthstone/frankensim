//! G0/G1 source, Hess-cycle and finite-stock checks for cassiterite reduction.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, Mass, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        FormationEnthalpySource,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SLUGS: [&str; 5] = [
    "cassiterite-solid-usgs-b2131-thermochemistry",
    "tin-white-solid-usgs-b2131",
    "graphite-solid-usgs-b2131",
    "carbon-monoxide-gas-usgs-b2131",
    "carbon-dioxide-gas-usgs-b2131",
];
const SPECIES: [&str; 5] = ["SnO2(s)", "Sn(s)", "C(s)", "CO(g)", "CO2(g)"];
const PHASES: [&str; 3] = [
    "source-phase-cassiterite",
    "source-phase-tin-white-tetragonal",
    "source-phase-graphite",
];

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn packs() -> &'static [NormalizedPack; 5] {
    static PACKS: OnceLock<[NormalizedPack; 5]> = OnceLock::new();
    PACKS.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory =
            std::env::temp_dir().join(format!("frankensim-tin-{}-{nonce}", std::process::id(),));
        std::fs::create_dir(&directory).unwrap();
        SLUGS.map(|slug| {
            let output = directory.join(format!("{slug}.fsmatpk"));
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
        })
    })
}

fn point(index: usize) -> QueryPoint {
    let mut query = QueryPoint::new();
    for (name, kind, value) in [
        ("temperature", QuantityKind::AbsoluteTemperature, 298.15),
        ("pressure", QuantityKind::Pressure, 100_000.0),
    ] {
        query = query
            .with_quantity(
                name,
                QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static)),
                value,
            )
            .unwrap();
    }
    for (axis, value) in [
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("source-pure-species", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (
            if index < 3 {
                "source-phase-solid"
            } else {
                "source-phase-gas"
            },
            1.0,
        ),
    ] {
        query = query
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if index < 3 {
        query = query
            .with_quantity(PHASES[index], QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    query
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-10 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn index(species: &SpeciesId) -> usize {
    SPECIES
        .iter()
        .position(|name| *name == species.as_str())
        .unwrap()
}

#[test]
fn literal_usgs_values_keep_their_source_phase_precision_and_limits() {
    // Independent transcription from original page images, not candidate output:
    // USGS B2131 pp19 (cassiterite), 9 (white tin), 5 (graphite), 15 (CO/CO2).
    let rows: [([f64; 4], &[(&str, Dims, f64)]); 5] = [
        (
            [0.150709, -577_600.0, -515_800.0, 49.0],
            &[
                (
                    "reported-formation-enthalpy-uncertainty-magnitude",
                    ENERGY,
                    200.0,
                ),
                (
                    "reported-formation-gibbs-energy-uncertainty-magnitude",
                    ENERGY,
                    200.0,
                ),
                (
                    "reported-standard-entropy-uncertainty-magnitude",
                    ENTROPY,
                    0.2,
                ),
            ],
        ),
        (
            [0.118710, 0.0, 0.0, 51.18],
            &[(
                "reported-standard-entropy-uncertainty-magnitude",
                ENTROPY,
                0.08,
            )],
        ),
        (
            [0.012011, 0.0, 0.0, 5.74],
            &[(
                "reported-standard-entropy-uncertainty-magnitude",
                ENTROPY,
                0.10,
            )],
        ),
        (
            [0.028010, -110_500.0, -137_100.0, 197.3],
            &[
                (
                    "reported-formation-enthalpy-uncertainty-magnitude",
                    ENERGY,
                    200.0,
                ),
                (
                    "reported-formation-gibbs-energy-uncertainty-magnitude",
                    ENERGY,
                    200.0,
                ),
                (
                    "reported-standard-entropy-uncertainty-magnitude",
                    ENTROPY,
                    0.0,
                ),
            ],
        ),
        (
            [0.044010, -393_500.0, -394_400.0, 213.8],
            &[
                (
                    "reported-formation-enthalpy-uncertainty-magnitude",
                    ENERGY,
                    100.0,
                ),
                (
                    "reported-formation-gibbs-energy-uncertainty-magnitude",
                    ENERGY,
                    200.0,
                ),
                (
                    "reported-standard-entropy-uncertainty-magnitude",
                    ENTROPY,
                    0.0,
                ),
            ],
        ),
    ];
    for (i, (values, magnitudes)) in rows.iter().enumerate() {
        let pack = &packs()[i];
        assert_eq!(pack.claims().claim_count(), 4 + magnitudes.len());
        let quantities = [
            ("molar-mass", MolarMass::DIMS, values[0]),
            ("standard-molar-formation-enthalpy", ENERGY, values[1]),
            ("standard-molar-formation-gibbs-energy", ENERGY, values[2]),
            ("standard-molar-entropy", ENTROPY, values[3]),
        ];
        for (property, dims, expected) in quantities.iter().chain(magnitudes.iter()) {
            let key = PropertyKey::new(*property, *dims);
            let query = point(i);
            let answer = pack
                .claims()
                .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            close(answer.evidence.value.value, *expected);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            let wrongs = [
                QueryPoint::new(),
                query
                    .clone()
                    .with_quantity(
                        "temperature",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        1300.0,
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
                        101325.0,
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
                query
                    .clone()
                    .with_quantity(
                        if i < 3 { PHASES[i] } else { "source-phase-gas" },
                        QuantitySpec::dimensional(Dims::NONE),
                        0.0,
                    )
                    .unwrap(),
            ];
            for wrong in wrongs {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err(),
                    "{property}: unsupported source state accepted"
                );
            }
        }
    }
}

#[test]
fn tin_reduction_and_boudouard_sources_close_the_independent_hess_cycle() {
    let points: [QueryPoint; 5] = std::array::from_fn(point);
    let mut heats = Vec::new();
    let mut gibbs = Vec::new();
    let mut entropies = Vec::new();
    for (recipe, expected_heat, expected_gibbs, expected_entropy, residual) in [
        (
            MineralReaction::CassiteriteCarbonMonoxideReduction,
            11_600.0,
            1_200.0,
            35.18,
            1e-6,
        ),
        (
            MineralReaction::GraphiteCarbonDioxideGasification,
            172_500.0,
            120_200.0,
            175.06,
            -1e-6,
        ),
        (
            MineralReaction::CassiteriteGraphiteReduction,
            356_600.0,
            241_600.0,
            385.30,
            -1e-6,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let caloric_sources = reaction
            .terms()
            .iter()
            .map(|term| FormationEnthalpySource {
                species: &term.species,
                pack: &packs()[index(&term.species)],
                at: &points[index(&term.species)],
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>();
        let gibbs_sources = reaction
            .terms()
            .iter()
            .map(|term| FormationGibbsSource {
                species: &term.species,
                pack: &packs()[index(&term.species)],
                at: &points[index(&term.species)],
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>();
        let thermal = reaction.standard_equilibrium(&caloric_sources).unwrap();
        let direct = reaction.standard_gibbs_change(&gibbs_sources).unwrap();
        close(thermal.enthalpy().molar().value(), expected_heat);
        close(thermal.molar_entropy().value(), expected_entropy);
        close(direct.molar().value(), expected_gibbs);
        close(
            thermal.molar_gibbs().value(),
            expected_heat - 298.15 * expected_entropy,
        );
        // Keep the source's independently rounded H/S and G evaluations distinct.
        assert!((direct.molar().value() - thermal.molar_gibbs().value()).abs() > 1.0);
        for term in direct.terms() {
            packs()
                .iter()
                .find(|pack| pack.pack_id() == term.pack_id)
                .unwrap()
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        let mut source_mass_residual = 0.0;
        for term in reaction.terms() {
            let i = index(&term.species);
            let mass = packs()[i]
                .claims()
                .query_typed(
                    &PropertyKey::new("molar-mass", MolarMass::DIMS),
                    &points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap()
                .evidence
                .value
                .value;
            source_mass_residual += term.coefficient as f64 * mass;
        }
        close(source_mass_residual, residual);
        heats.push(thermal.enthalpy().molar().value());
        gibbs.push(direct.molar().value());
        entropies.push(thermal.molar_entropy().value());
    }
    close(heats[0] + 2.0 * heats[1], heats[2]);
    close(gibbs[0] + 2.0 * gibbs[1], gibbs[2]);
    close(entropies[0] + 2.0 * entropies[1], entropies[2]);
    let reaction = MineralReaction::CassiteriteGraphiteReduction
        .build()
        .unwrap();
    close(
        reaction
            .extent_from_reacted_mass(
                &SpeciesId::new("SnO2(s)").unwrap(),
                Mass::new(0.150709),
                MolarMass::new(0.150709),
            )
            .unwrap()
            .value(),
        1.0,
    );
}

#[test]
fn reducing_gas_must_exist_before_the_cassiterite_stage_can_consume_it() {
    let gasify = MineralReaction::GraphiteCarbonDioxideGasification
        .build()
        .unwrap();
    let reduce = MineralReaction::CassiteriteCarbonMonoxideReduction
        .build()
        .unwrap();
    let network = ReactionNetwork::new(vec![gasify.clone(), reduce.clone()]).unwrap();
    let inventory = network
        .species()
        .iter()
        .map(|species| {
            (
                species.clone(),
                Amount::new(match species.as_str() {
                    "SnO2(s)" => 1.0,
                    "C(s)" | "CO2(g)" => 2.0,
                    _ => 0.0,
                }),
            )
        })
        .collect::<Vec<_>>();
    let extents = |g, r| {
        [
            (gasify.id().clone(), Amount::new(g)),
            (reduce.id().clone(), Amount::new(r)),
        ]
    };
    assert!(matches!(
        network.advance_inventory(&inventory, &extents(2.0, 1.0)),
        Err(ReactionNetworkError::InsufficientInventory { .. })
    ));
    assert!(
        network
            .advance_inventory(&inventory, &extents(2.1, 0.0))
            .is_err()
    );
    let gas = network
        .advance_inventory(&inventory, &extents(2.0, 0.0))
        .unwrap();
    let metal = network
        .advance_inventory(&gas.final_inventory, &extents(0.0, 1.0))
        .unwrap();
    let expected = [
        ("SnO2(s)", 0.0),
        ("Sn(s)", 1.0),
        ("C(s)", 0.0),
        ("CO(g)", 2.0),
        ("CO2(g)", 2.0),
    ];
    for (species, amount) in expected {
        close(
            metal
                .final_inventory
                .iter()
                .find(|(s, _)| s.as_str() == species)
                .unwrap()
                .1
                .value(),
            amount,
        );
    }
    assert!(
        network
            .advance_inventory(&metal.final_inventory, &extents(0.0, 0.1))
            .is_err()
    );
    // The aggregate carbon endpoint consumes the original finite graphite,
    // and cannot silently substitute CO or another carbon-bearing species.
    let aggregate = MineralReaction::CassiteriteGraphiteReduction
        .build()
        .unwrap();
    let initial = aggregate
        .terms()
        .iter()
        .map(|term| {
            (
                term.species.clone(),
                Amount::new(match term.species.as_str() {
                    "SnO2(s)" => 1.0,
                    "C(s)" => 1.0,
                    _ => 0.0,
                }),
            )
        })
        .collect::<Vec<_>>();
    assert!(
        aggregate
            .advance_inventory(&initial, Amount::new(1.0))
            .is_err()
    );
    let half = aggregate
        .advance_inventory(&initial, Amount::new(0.5))
        .unwrap();
    close(
        half.iter()
            .find(|(s, _)| s.as_str() == "Sn(s)")
            .unwrap()
            .1
            .value(),
        0.5,
    );
}
