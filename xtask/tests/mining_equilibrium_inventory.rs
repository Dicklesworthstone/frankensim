//! Finite iron hydrolysis equilibrium from the actual USGS source packs.
//! Fixed aqueous amount scale and gamma=1 are declared model assumptions.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedPack, QueryPoint, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            equilibrium::{
                EquilibriumExtentStatus, EquilibriumInventory, EquilibriumSolveOptions,
                IdealActivityModel, SingleReactionEquilibrium,
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const SLUGS: [&str; 4] = [
    "iron-iii-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
    "goethite-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
];
const SPECIES: [&str; 4] = ["Fe[3+](aq)", "H2O(l)", "FeO(OH)(s)", "H+(aq)"];

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
        "frankensim-mineral-equilibrium-{}-{nonce}-{}",
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
    assert_eq!(pack.pack_id(), slug);
    pack
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
    match index {
        1 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        2 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-goethite", 1.0),
            ("source-pure-species", 1.0),
        ]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}

fn solve(packs: &[NormalizedPack; 4], amounts: [f64; 4]) -> SingleReactionEquilibrium {
    let ids = SPECIES.map(|name| SpeciesId::new(name).unwrap());
    let points = std::array::from_fn::<_, 4, _>(point);
    let sources = std::array::from_fn::<_, 4, _>(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let inventory = std::array::from_fn::<_, 4, _>(|i| EquilibriumInventory {
        species: ids[i].clone(),
        amount: Amount::new(amounts[i]),
        activity: if i == 0 || i == 3 {
            IdealActivityModel::IdealDilute {
                amount_at_unit_activity: Amount::new(1.0),
            }
        } else {
            IdealActivityModel::PureCondensed
        },
    });
    MineralReaction::FerricHydrolysisGoethite
        .build()
        .unwrap()
        .equilibrate_ideal(
            &sources,
            &inventory,
            EquilibriumSolveOptions {
                extent_tolerance: Amount::new(1e-12),
                log_residual_tolerance: Dimensionless::new(1e-9),
                max_iterations: 200,
            },
        )
        .unwrap()
}

fn amounts(result: &SingleReactionEquilibrium) -> [f64; 4] {
    SPECIES.map(|name| {
        result
            .final_inventory
            .iter()
            .find(|(species, _)| species.as_str() == name)
            .unwrap()
            .1
            .value()
    })
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
}

#[test]
fn g1_source_driven_precipitation_solves_extent_and_preserves_iron_acid_and_water() {
    let packs = SLUGS.map(compile);
    let initial = [0.001, 55.5, 0.0, 0.001];
    let result = solve(&packs, initial);
    assert_eq!(result.status, EquilibriumExtentStatus::Interior);
    assert!(result.extent.value() > 0.0 && result.extent.value() < initial[0]);
    assert_eq!(result.standard.molar().value(), -900.0);
    assert!(result.extent_bracket[1].value() - result.extent_bracket[0].value() <= 1e-12);
    assert!(result.log_quotient_over_equilibrium.unwrap().value().abs() <= 1e-9);
    assert!(result.iterations > 0 && result.iterations <= 200);
    let final_amounts = amounts(&result);
    for ((start, finish), coefficient) in initial
        .into_iter()
        .zip(final_amounts)
        .zip([-1.0, -2.0, 1.0, 3.0])
    {
        assert!(finish >= 0.0);
        close(finish, coefficient * result.extent.value() + start);
    }
    // Independent conserved totals: iron, oxygen, hydrogen, then charge.
    for stocks in [initial, final_amounts] {
        close(stocks[0] + stocks[2], initial[0] + initial[2]);
        close(stocks[1] + 2.0 * stocks[2], initial[1] + 2.0 * initial[2]);
        close(
            2.0 * stocks[1] + stocks[2] + stocks[3],
            2.0 * initial[1] + initial[2] + initial[3],
        );
        close(3.0 * stocks[0] + stocks[3], 3.0 * initial[0] + initial[3]);
    }
    // Independent mass-action check for the declared one-mole aqueous scale.
    let q = final_amounts[3].powi(3) / final_amounts[0];
    let log_k = 900.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15);
    assert!((q.ln() - log_k).abs() <= 1e-9);
    for term in result.standard.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        assert!(!term.answer.receipt.source_hashes.is_empty());
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    for entry in &result.initial {
        let index = SPECIES
            .iter()
            .position(|name| *name == entry.species.as_str())
            .unwrap();
        assert_eq!(entry.amount.value(), initial[index]);
        assert_eq!(
            entry.activity,
            if index == 0 || index == 3 {
                IdealActivityModel::IdealDilute {
                    amount_at_unit_activity: Amount::new(1.0),
                }
            } else {
                IdealActivityModel::PureCondensed
            }
        );
    }
}

#[test]
fn g1_g3_finite_goethite_dissolves_to_a_phase_limit_and_absence_preserves_inequality() {
    let packs = SLUGS.map(compile);
    let initial = [1e-8, 55.5, 1e-7, 0.01];
    let dissolved = solve(&packs, initial);
    assert_eq!(
        dissolved.status,
        EquilibriumExtentStatus::ReverseInventoryLimit
    );
    assert_eq!(dissolved.extent.value(), -initial[2]);
    assert_eq!(dissolved.extent_bracket, [dissolved.extent; 2]);
    assert_eq!(dissolved.iterations, 0);
    let final_amounts = amounts(&dissolved);
    assert_eq!(final_amounts[2], 0.0);
    close(final_amounts[0], initial[0] + initial[2]);
    close(final_amounts[1], initial[1] + 2.0 * initial[2]);
    close(final_amounts[3], initial[3] - 3.0 * initial[2]);
    let residual = (final_amounts[3].powi(3) / final_amounts[0]).ln()
        - 900.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15);
    assert!(residual > 1.0);
    close(
        dissolved.log_quotient_over_equilibrium.unwrap().value(),
        residual,
    );

    // No solid remains to fund further dissolution. Q=K must not be invented.
    let absent = solve(&packs, final_amounts);
    assert_eq!(
        absent.status,
        EquilibriumExtentStatus::ReverseInventoryLimit
    );
    assert_eq!(absent.extent.value(), 0.0);
    assert_eq!(amounts(&absent), final_amounts);
    assert!(absent.log_quotient_over_equilibrium.unwrap().value() > 1.0);
}
