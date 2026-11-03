//! Actual USGS packs close zinc roast, acid leach, electrowinning and ferric attack.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, ElectricCharge, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction, FARADAY_C_PER_MOL, FormationEnthalpySource,
        electrolysis::ElectrolysisReaction,
        minerals::MineralReaction,
        network::{ReactionNetwork, ReactionNetworkError},
        thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SLUGS: [&str; 13] = [
    "sphalerite-solid-usgs-b2131-thermochemistry",
    "zincite-solid-usgs-b2131-thermochemistry",
    "zinc-ion-aqueous-usgs-b2131",
    "zinc-solid-usgs-b2131",
    "sulfur-dioxide-gas-usgs-b2131",
    "oxygen-gas-usgs-b2131",
    "hydrogen-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
    "iron-iii-ion-aqueous-usgs-b2131",
    "iron-ii-ion-aqueous-usgs-b2131",
    "sulfur-orthorhombic-solid-usgs-b2131",
    "copper-ii-ion-aqueous-usgs-b2131",
    "copper-solid-usgs-b2131",
];
const SPECIES: [&str; 13] = [
    "ZnS(s)",
    "ZnO(s)",
    "Zn[2+](aq)",
    "Zn(s)",
    "SO2(g)",
    "O2(g)",
    "H+(aq)",
    "H2O(l)",
    "Fe[3+](aq)",
    "Fe[2+](aq)",
    "S(s)",
    "Cu[2+](aq)",
    "Cu(s)",
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
        "frankensim-zinc-{}-{nonce}-{}",
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
        2 | 6 | 8 | 9 | 11 => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
        4 | 5 => flags.extend([("source-phase-gas", 1.0), ("source-pure-species", 1.0)]),
        7 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]),
    }
    match index {
        0 => flags.push(("source-phase-sphalerite", 1.0)),
        1 => flags.push(("source-phase-zincite", 1.0)),
        10 => flags.push(("source-phase-sulfur-orthorhombic", 1.0)),
        _ => {}
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn index(species: &SpeciesId) -> usize {
    SPECIES
        .iter()
        .position(|name| *name == species.as_str())
        .unwrap()
}

fn gibbs_sources<'a>(
    reaction: &'a ConservedReaction,
    packs: &'a [NormalizedPack; 13],
    points: &'a [QueryPoint; 13],
) -> Vec<FormationGibbsSource<'a>> {
    reaction
        .terms()
        .iter()
        .map(|term| FormationGibbsSource {
            species: &term.species,
            pack: &packs[index(&term.species)],
            at: &points[index(&term.species)],
            reference_axis: REFERENCE,
        })
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn literal_zinc_sources_retain_phase_reference_and_reported_uncertainties() {
    // Independent transcription of B2131 printed pp10,13,18,20, including
    // orthorhombic-S sphalerite rather than its S2(g)-reference high-T table.
    for (i, values, magnitudes) in [
        (
            0,
            [0.097456, -204_100.0, -199_600.0, 58.7],
            [1500.0, 1500.0, 0.2],
        ),
        (
            1,
            [0.081389, -350_500.0, -320_400.0, 43.2],
            [300.0, 300.0, 0.1],
        ),
        (
            2,
            [0.06539, -153_400.0, -147_300.0, -109.8],
            [200.0, 200.0, 0.5],
        ),
        (3, [0.06539, 0.0, 0.0, 41.63], [0.0, 0.0, 0.15]),
        (
            4,
            [0.064065, -296_800.0, -300_100.0, 248.2],
            [200.0, 200.0, 0.1],
        ),
    ] {
        let pack = compile(SLUGS[i]);
        let query = point(i);
        assert_eq!(pack.claims().claim_count(), if i == 3 { 5 } else { 7 });
        let entropy = if i == 2 {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        let entropy_uncertainty = if i == 2 {
            "reported-conventional-entropy-uncertainty-magnitude"
        } else {
            "reported-standard-entropy-uncertainty-magnitude"
        };
        let mut properties = vec![
            ("molar-mass", MolarMass::DIMS, values[0]),
            ("standard-molar-formation-enthalpy", ENERGY, values[1]),
            ("standard-molar-formation-gibbs-energy", ENERGY, values[2]),
            (entropy, ENTROPY, values[3]),
            (entropy_uncertainty, ENTROPY, magnitudes[2]),
        ];
        if i != 3 {
            properties.extend([
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
            ]);
        }
        for (property, dims, expected) in properties {
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
                QueryPoint::new(),
                query.clone().with("temperature", 900.0).unwrap(),
                query.clone().with("pressure", 101_325.0).unwrap(),
                query.clone().with(REFERENCE, 0.0).unwrap(),
                query
                    .clone()
                    .with("source-specimen-qualified", 1.0)
                    .unwrap(),
                point(if i == 0 {
                    1
                } else if i == 3 {
                    2
                } else {
                    0
                }),
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
fn actual_sources_close_roast_leach_electrowinning_hess_cycle_and_ferric_attack() {
    let packs = SLUGS.map(compile);
    let points = std::array::from_fn(point);
    let mut heats = Vec::new();
    let mut energies = Vec::new();
    for (recipe, expected_heat, expected_gibbs) in [
        (
            MineralReaction::SphaleriteRoastingToZincite,
            -886_400.0,
            -841_800.0,
        ),
        (MineralReaction::ZinciteAcidLeach, -88_700.0, -64_000.0),
        (MineralReaction::ZincElectrowinning, 878_400.0, 768_800.0),
        (MineralReaction::FerricSphaleriteLeach, -31_700.0, -94_300.0),
        (
            MineralReaction::CopperCementationWithZinc,
            -218_300.0,
            -212_400.0,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let sources = reaction
            .terms()
            .iter()
            .map(|term| FormationEnthalpySource {
                species: &term.species,
                pack: &packs[index(&term.species)],
                at: &points[index(&term.species)],
                reference_axis: REFERENCE,
            })
            .collect::<Vec<_>>();
        let heat = reaction.standard_enthalpy_change(&sources).unwrap();
        let gibbs = reaction
            .standard_gibbs_change(&gibbs_sources(&reaction, &packs, &points))
            .unwrap();
        close(heat.molar().value(), expected_heat);
        close(gibbs.molar().value(), expected_gibbs);
        for term in gibbs.terms() {
            let pack = packs.iter().find(|p| p.pack_id() == term.pack_id).unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        }
        heats.push(heat.molar().value());
        energies.push(gibbs.molar().value());
    }
    // Net 2 ZnS + 2 O2 -> 2 Zn + 2 SO2 after the intermediate ZnO/Zn2+
    // and the acid/water carrier cancel. Source-rounded masses stay unchanged.
    close(heats[0] + 2.0 * heats[1] + heats[2], -185_400.0);
    close(energies[0] + 2.0 * energies[1] + energies[2], -201_000.0);
}

#[test]
fn zinc_stages_cannot_borrow_future_calcine_dissolved_zinc_or_regenerated_acid() {
    let roast = MineralReaction::SphaleriteRoastingToZincite
        .build()
        .unwrap();
    let leach = MineralReaction::ZinciteAcidLeach.build().unwrap();
    let winning = MineralReaction::ZincElectrowinning.build().unwrap();
    let network =
        ReactionNetwork::new(vec![roast.clone(), leach.clone(), winning.clone()]).unwrap();
    let initial = network
        .species()
        .iter()
        .map(|id| {
            (
                id.clone(),
                Amount::new(match id.as_str() {
                    "ZnS(s)" => 2.0,
                    "O2(g)" => 3.0,
                    "H+(aq)" => 4.0,
                    _ => 0.0,
                }),
            )
        })
        .collect::<Vec<_>>();
    let extents = |r, l, w| {
        [
            (roast.id().clone(), Amount::new(r)),
            (leach.id().clone(), Amount::new(l)),
            (winning.id().clone(), Amount::new(w)),
        ]
    };
    assert!(matches!(
        network.advance_inventory(&initial, &extents(1.0, 2.0, 1.0)),
        Err(ReactionNetworkError::InsufficientInventory { .. })
    ));
    let roasted = network
        .advance_inventory(&initial, &extents(1.0, 0.0, 0.0))
        .unwrap();
    let mut without_acid = roasted.final_inventory.clone();
    without_acid
        .iter_mut()
        .find(|(id, _)| id.as_str() == "H+(aq)")
        .unwrap()
        .1 = Amount::new(0.0);
    assert!(
        network
            .advance_inventory(&without_acid, &extents(0.0, 2.0, 1.0))
            .is_err()
    );
    let leached = network
        .advance_inventory(&roasted.final_inventory, &extents(0.0, 2.0, 0.0))
        .unwrap();
    let finished = network
        .advance_inventory(&leached.final_inventory, &extents(0.0, 0.0, 1.0))
        .unwrap();
    for (id, amount) in finished.final_inventory {
        let expected = match id.as_str() {
            "Zn(s)" | "SO2(g)" => 2.0,
            "O2(g)" => 1.0,
            "H+(aq)" => 4.0,
            _ => 0.0,
        };
        close(amount.value(), expected);
    }
}

#[test]
fn copper_cementation_spends_finite_zinc_and_credits_dissolved_zinc() {
    let reaction = MineralReaction::CopperCementationWithZinc.build().unwrap();
    let initial = [
        (SpeciesId::new("Cu[2+](aq)").unwrap(), Amount::new(2.0)),
        (SpeciesId::new("Zn(s)").unwrap(), Amount::new(1.5)),
        (SpeciesId::new("Cu(s)").unwrap(), Amount::new(0.0)),
        (SpeciesId::new("Zn[2+](aq)").unwrap(), Amount::new(10.0)),
    ];
    let limit = reaction.limiting_extent(&initial[..2]).unwrap();
    close(limit.value(), 1.5);
    assert!(
        reaction
            .advance_inventory(&initial, Amount::new(1.5001))
            .is_err()
    );
    let final_inventory = reaction.advance_inventory(&initial, limit).unwrap();
    for ((id, _), expected) in initial.iter().zip([0.5, 0.0, 1.5, 11.5]) {
        close(
            final_inventory
                .iter()
                .find(|(s, _)| s == id)
                .unwrap()
                .1
                .value(),
            expected,
        );
    }
    // Existing dissolved zinc cannot replace the metallic reductant. Zero
    // zinc dust leaves copper in solution even when Zn(II) is abundant.
    let mut no_reductant = initial.clone();
    no_reductant[1].1 = Amount::new(0.0);
    close(
        reaction
            .limiting_extent(&no_reductant[..2])
            .unwrap()
            .value(),
        0.0,
    );
    assert!(
        reaction
            .advance_inventory(&no_reductant, Amount::new(0.1))
            .is_err()
    );
    assert_eq!(initial[1].1.value(), 1.5);
}

#[test]
fn source_zinc_electrowinning_uses_two_electrons_and_preserves_finite_stocks() {
    let packs = SLUGS.map(compile);
    let points = std::array::from_fn(point);
    let reaction = MineralReaction::ZincElectrowinning.build().unwrap();
    let mass = packs[3]
        .claims()
        .query_typed(
            &PropertyKey::new("molar-mass", MolarMass::DIMS),
            &points[3],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    let model = ElectrolysisReaction::new(
        reaction.clone(),
        SpeciesId::new("Zn(s)").unwrap(),
        MolarMass::new(mass.evidence.value.value),
        2,
    )
    .unwrap();
    assert_eq!(model.electrons_per_equation(), 4);
    let potential = model
        .standard_potential(&gibbs_sources(&reaction, &packs, &points))
        .unwrap();
    close(potential.potential().value(), -1.992_012_627_933_59);
    let initial = reaction
        .terms()
        .iter()
        .map(|t| {
            (
                t.species.clone(),
                Amount::new(if t.coefficient < 0 { 2.0 } else { 0.0 }),
            )
        })
        .collect::<Vec<_>>();
    let charge = ElectricCharge::new(4.0 * FARADAY_C_PER_MOL);
    let result = model
        .advance_charge(&initial, charge, Dimensionless::new(1.0))
        .unwrap();
    close(result.product_yield.amount().value(), 2.0);
    close(result.product_yield.mass().value(), 0.13078);
    assert!(
        model
            .advance_charge(
                &initial,
                ElectricCharge::new(4.01 * FARADAY_C_PER_MOL),
                Dimensionless::new(1.0)
            )
            .is_err()
    );
    let partial = model
        .advance_charge(&initial, charge, Dimensionless::new(0.5))
        .unwrap();
    close(partial.product_yield.amount().value(), 1.0);
    close(partial.unallocated_charge.value(), 2.0 * FARADAY_C_PER_MOL);
    close(
        partial
            .modeled_final_inventory
            .iter()
            .find(|(s, _)| s.as_str() == "H+(aq)")
            .unwrap()
            .1
            .value(),
        2.0,
    );
}
