//! Matched NBS aqueous species consumed by conserved reactions and ideal closure.
//! Numerical rows are literal primary-source values; no pH or gas-transfer fit.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        FormationEnthalpySource, ReactionBalanceError,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            equilibrium::{
                EquilibriumExtentStatus, EquilibriumInventory, EquilibriumSolveOptions,
                IdealActivityModel,
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const ACID: &[&str] = &["HCN(aq)", "H+(aq)", "CN[-](aq)"];
const WATER: &[&str] = &["H2O(l)", "H+(aq)", "OH[-](aq)"];
const GOLD: &[&str] = &[
    "Au(s)",
    "CN[-](aq)",
    "O2(aq)",
    "H2O(l)",
    "Au(CN)2[-](aq)",
    "OH[-](aq)",
];

struct Row {
    species: &'static str,
    pack: &'static str,
    phase: &'static str,
    ionic: bool,
    // Literal g/mol, kJ/mol, kJ/mol, J/(mol K), not a fitted reconciliation.
    values: [f64; 4],
}

// NIST corrected CSV, DOI 10.18434/M32124; original pp. 2-37, 2-38,
// 2-89 and 2-164. HCN ao is neutral molecular solute; HCN ai is NOT that row.
const ROWS: &[Row] = &[
    Row {
        species: "HCN(aq)",
        pack: "hydrogen-cyanide-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        ionic: false,
        values: [27.0259, 107.1, 119.7, 124.7],
    },
    Row {
        species: "H+(aq)",
        pack: "hydrogen-ion-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        ionic: true,
        values: [1.0080, 0.0, 0.0, 0.0],
    },
    Row {
        species: "O2(aq)",
        pack: "oxygen-aqueous-nbs1989-thermochemistry",
        phase: "source-phase-aqueous",
        ionic: false,
        values: [31.9988, -11.7, 16.4, 110.9],
    },
    Row {
        species: "CN[-](aq)",
        pack: "cyanide-ion-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        ionic: true,
        values: [26.0179, 150.6, 172.4, 94.1],
    },
    Row {
        species: "H2O(l)",
        pack: "water-liquid-nbs1989-thermochemistry",
        phase: "source-phase-liquid",
        ionic: false,
        values: [18.0154, -285.830, -237.129, 69.91],
    },
    Row {
        species: "OH[-](aq)",
        pack: "hydroxide-ion-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        ionic: true,
        values: [17.0074, -229.994, -157.244, -10.75],
    },
    Row {
        species: "Au(s)",
        pack: "gold-solid-nbs1989-thermochemistry",
        phase: "source-phase-solid",
        ionic: false,
        values: [196.9670, 0.0, 0.0, 47.40],
    },
    Row {
        species: "Au(CN)2[-](aq)",
        pack: "dicyanoaurate-ion-aqueous-nbs1989",
        phase: "source-phase-aqueous",
        ionic: true,
        values: [249.0028, 242.3, 285.8, 172.0],
    },
    Row {
        species: "O2(g)",
        pack: "oxygen-gas-nbs1989-thermochemistry",
        phase: "source-phase-gas",
        ionic: false,
        values: [31.9988, 0.0, 0.0, 205.138],
    },
];

struct Source {
    row: &'static Row,
    id: SpeciesId,
    pack: NormalizedPack,
    at: QueryPoint,
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

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(row: &Row) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            298.15,
        )
        .unwrap()
        .with_quantity("pressure", semantic(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        (row.phase, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if row.phase == "source-phase-aqueous" {
        flags.push(("source-aqueous-unit-activity-infinite-dilution", 1.0));
        if row.ionic {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    } else {
        flags.push(("source-pure-species", 1.0));
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
            "frankensim-aqueous-cyanide-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        ROWS.iter()
            .map(|row| {
                let target = out.join(format!("{}.fsmatpk", row.pack));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{}/manifest.tsv", row.pack)))
                    .arg("--out")
                    .arg(&target)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{}: {}",
                    row.pack,
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(target).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), row.pack);
                assert_eq!(pack.to_bytes(), bytes);
                assert_eq!(pack.claims().claims_ordered().count(), 4);
                Source {
                    row,
                    id: SpeciesId::new(row.species).unwrap(),
                    pack,
                    at: point(row),
                }
            })
            .collect()
    })
}

fn source(id: &str) -> &'static Source {
    sources().iter().find(|s| s.row.species == id).unwrap()
}

fn heat_sources(ids: &[&str]) -> Vec<FormationEnthalpySource<'static>> {
    ids.iter()
        .map(|id| {
            let s = source(id);
            FormationEnthalpySource {
                species: &s.id,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}

fn gibbs_sources(ids: &[&str]) -> Vec<FormationGibbsSource<'static>> {
    ids.iter()
        .map(|id| {
            let s = source(id);
            FormationGibbsSource {
                species: &s.id,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}

fn inventory(ids: &[&str], amounts: &[f64], scale: f64) -> Vec<EquilibriumInventory> {
    ids.iter()
        .zip(amounts)
        .map(|(id, amount)| {
            let s = source(id);
            EquilibriumInventory {
                species: s.id.clone(),
                amount: Amount::new(*amount),
                activity: if s.row.phase == "source-phase-aqueous" {
                    IdealActivityModel::IdealDilute {
                        amount_at_unit_activity: Amount::new(scale),
                    }
                } else {
                    IdealActivityModel::PureCondensed
                },
            }
        })
        .collect()
}

fn options() -> EquilibriumSolveOptions {
    EquilibriumSolveOptions {
        extent_tolerance: Amount::new(1e-15),
        log_residual_tolerance: Dimensionless::new(1e-10),
        max_iterations: 256,
    }
}

fn amount(values: &[(SpeciesId, Amount)], name: &str) -> f64 {
    values
        .iter()
        .find(|(id, _)| id.as_str() == name)
        .unwrap()
        .1
        .value()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn small_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_literal_rows_keep_neutral_ionic_and_aqueous_reference_states_distinct() {
    for s in sources() {
        let entropy = if s.row.ionic {
            "conventional-standard-molar-entropy"
        } else {
            "standard-molar-entropy"
        };
        let [mass, h, g, entropy_value] = s.row.values;
        for (key, dims, expected) in [
            ("molar-mass", MASS, mass * 0.001),
            ("standard-molar-formation-enthalpy", ENERGY, h * 1000.0),
            ("standard-molar-formation-gibbs-energy", ENERGY, g * 1000.0),
            (entropy, ENTROPY, entropy_value),
        ] {
            let answer = s
                .pack
                .claims()
                .query_typed(
                    &PropertyKey::new(key, dims),
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
    }
    for s in &sources()[..3] {
        let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
        for at in [
            s.at.clone().with("source-phase-aqueous", 0.0).unwrap(),
            s.at.clone()
                .with("source-aqueous-unit-activity-infinite-dilution", 0.0)
                .unwrap(),
            s.at.clone().with(REFERENCE, 0.0).unwrap(),
            s.at.clone()
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::AbsoluteTemperature),
                    300.0,
                )
                .unwrap(),
        ] {
            assert!(
                s.pack
                    .claims()
                    .query_typed(&key, &at, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        let wrong_entropy = if s.row.ionic {
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        assert!(
            s.pack
                .claims()
                .query_typed(
                    &PropertyKey::new(wrong_entropy, ENTROPY),
                    &s.at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
    // An aqueous query cannot substitute the gas pack's elemental-reference zero.
    assert!(
        source("O2(g)")
            .pack
            .claims()
            .query_typed(
                &PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY),
                &source("O2(aq)").at,
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
}

#[test]
fn g1_hcn_dissociation_uses_the_molecular_row_and_matches_a_finite_stock_root() {
    let reaction = MineralReaction::HydrogenCyanideDissociation
        .build()
        .unwrap();
    let standard = reaction
        .standard_gibbs_change(&gibbs_sources(ACID))
        .unwrap();
    let caloric = reaction.standard_equilibrium(&heat_sources(ACID)).unwrap();
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    close(standard.molar().value(), 52_700.0);
    close(standard.log_equilibrium_constant().value(), -52_700.0 / rt);
    close(
        -standard.log_equilibrium_constant().value() / 10.0_f64.ln(),
        9.232640551245206,
    );
    close(caloric.enthalpy().molar().value(), 43_500.0);
    close(caloric.molar_entropy().value(), -30.6);
    close(caloric.molar_gibbs().value(), 52_623.39);
    let ka = (-52_700.0 / rt).exp();
    let mut free_cyanide = Vec::new();
    for proton in [0.0, 1e-4] {
        let total = 0.01;
        let result = reaction
            .equilibrate_ideal(
                &gibbs_sources(ACID),
                &inventory(ACID, &[total, proton, 0.0], 1.0),
                options(),
            )
            .unwrap();
        assert_eq!(result.status, EquilibriumExtentStatus::Interior);
        assert!(result.log_quotient_over_equilibrium.unwrap().value().abs() <= 1e-10);
        // Independent stable quadratic: (H0+x)*x/(HCN0-x) = Ka for n_ref=1.
        let b = proton + ka;
        let expected = 2.0 * ka * total / (b + (b * b + 4.0 * ka * total).sqrt());
        small_close(result.extent.value(), expected);
        let hcn = amount(&result.final_inventory, "HCN(aq)");
        let h = amount(&result.final_inventory, "H+(aq)");
        let cn = amount(&result.final_inventory, "CN[-](aq)");
        close(hcn + cn, total);
        close(h - cn, proton);
        small_close(h * cn / hcn, ka);
        free_cyanide.push(cn);
        // Scaling solvent reference AND stocks preserves activities and doubles moles.
        let scaled = reaction
            .equilibrate_ideal(
                &gibbs_sources(ACID),
                &inventory(ACID, &[2.0 * total, 2.0 * proton, 0.0], 2.0),
                options(),
            )
            .unwrap();
        small_close(scaled.extent.value(), 2.0 * result.extent.value());
    }
    assert!(free_cyanide[1] < free_cyanide[0]);
}

#[test]
fn g1_water_dissociation_has_a_source_derived_kw_and_conserves_charge_and_solvent() {
    let reaction = MineralReaction::WaterDissociation.build().unwrap();
    let standard = reaction
        .standard_gibbs_change(&gibbs_sources(WATER))
        .unwrap();
    let caloric = reaction.standard_equilibrium(&heat_sources(WATER)).unwrap();
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    close(standard.molar().value(), 79_885.0);
    close(standard.log_equilibrium_constant().value(), -79_885.0 / rt);
    close(caloric.enthalpy().molar().value(), 55_836.0);
    close(caloric.molar_entropy().value(), -80.66);
    close(caloric.molar_gibbs().value(), 79_884.779);
    let root_kw = (-79_885.0 / (2.0 * rt)).exp();
    for scale in [1.0, 2.0] {
        let water = 55.5 * scale;
        let result = reaction
            .equilibrate_ideal(
                &gibbs_sources(WATER),
                &inventory(WATER, &[water, 0.0, 0.0], scale),
                options(),
            )
            .unwrap();
        assert_eq!(result.status, EquilibriumExtentStatus::Interior);
        assert!(result.log_quotient_over_equilibrium.unwrap().value().abs() <= 1e-10);
        let h = amount(&result.final_inventory, "H+(aq)");
        let oh = amount(&result.final_inventory, "OH[-](aq)");
        small_close(h, scale * root_kw);
        small_close(oh, h);
        close(amount(&result.final_inventory, "H2O(l)") + h, water);
    }
    // No solvent means no dissociation; zero products cannot synthesize a stock.
    let dry = reaction
        .equilibrate_ideal(
            &gibbs_sources(WATER),
            &inventory(WATER, &[0.0, 0.0, 0.0], 1.0),
            options(),
        )
        .unwrap();
    assert_eq!(dry.status, EquilibriumExtentStatus::InventoryBlocked);
    assert_eq!(dry.extent.value(), 0.0);
}

#[test]
fn g1_dissolved_oxygen_changes_gold_thermodynamics_and_bounds_available_gold() {
    let reaction = MineralReaction::GoldCyanidationDissolvedOxygen
        .build()
        .unwrap();
    let standard = reaction
        .standard_gibbs_change(&gibbs_sources(GOLD))
        .unwrap();
    let caloric = reaction.standard_equilibrium(&heat_sources(GOLD)).unwrap();
    close(standard.molar().value(), -407_118.0);
    close(caloric.enthalpy().molar().value(), -572_216.0);
    close(caloric.molar_entropy().value(), -548.12);
    close(caloric.molar_gibbs().value(), -408_794.022);
    let gas_ids = [
        "Au(s)",
        "CN[-](aq)",
        "O2(g)",
        "H2O(l)",
        "Au(CN)2[-](aq)",
        "OH[-](aq)",
    ];
    let gas = MineralReaction::GoldCyanidation
        .build()
        .unwrap()
        .standard_gibbs_change(&gibbs_sources(&gas_ids))
        .unwrap();
    // Subtracting the two conserved equations is O2(g) -> O2(aq).
    let dissolution_g = gas.molar().value() - standard.molar().value();
    close(dissolution_g, 16_400.0);
    let dissolution_k = (-dissolution_g / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15)).exp();
    close(dissolution_k, 0.0013391965031285008);
    // K relates aqueous activity to gas fugacity/p0; no molarity or transfer rate.
    for dissolved in [0.0, 0.00025] {
        let amounts = [4.0, 8.0, dissolved, 55.5, 0.0, 0.0];
        let initial = GOLD
            .iter()
            .zip(amounts)
            .map(|(id, n)| (source(id).id.clone(), Amount::new(n)))
            .collect::<Vec<_>>();
        let extent = reaction.limiting_extent(&initial[..4]).unwrap();
        close(extent.value(), dissolved);
        let result = reaction.advance_inventory(&initial, extent).unwrap();
        close(amount(&result, "Au(s)"), 4.0 - 4.0 * dissolved);
        close(amount(&result, "CN[-](aq)"), 8.0 - 8.0 * dissolved);
        assert_eq!(amount(&result, "O2(aq)"), 0.0);
        close(amount(&result, "Au(CN)2[-](aq)"), 4.0 * dissolved);
        close(amount(&result, "OH[-](aq)"), 4.0 * dissolved);
        close(
            caloric.enthalpy().at_extent(extent).unwrap().value(),
            -572_216.0 * dissolved,
        );
        assert!(matches!(
            reaction.advance_inventory(&initial, Amount::new(dissolved + 0.0001)),
            Err(ReactionBalanceError::ExtentExceedsInventory { .. })
        ));
    }
}
