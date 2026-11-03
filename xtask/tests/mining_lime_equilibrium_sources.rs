//! Actual NBS lime sources through conserved endpoints and finite solid/shared pH.
//! All finite inventories below are authored scenarios, not measured plant doses.
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
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        FormationEnthalpySource,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon,
                WaterAutoprotolysis,
                precipitation::{
                    HydrolysisPrecipitation, PrecipitationOptions, PrecipitationStatus,
                },
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const FAMILY: &str = "uncomplexed-cyanide";
const PACKS: [&str; 8] = [
    "calcium-oxide-solid-nbs1989-thermochemistry",
    "calcium-hydroxide-solid-nbs1989-thermochemistry",
    "calcium-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
    "hydrogen-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "hydrogen-cyanide-aqueous-nbs1989",
    "cyanide-ion-aqueous-nbs1989",
];
const SPECIES: [&str; 8] = [
    "CaO(s)",
    "Ca(OH)2(s)",
    "Ca[2+](aq)",
    "H2O(l)",
    "H+(aq)",
    "OH[-](aq)",
    "HCN(aq)",
    "CN[-](aq)",
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
    match index {
        0 | 1 => flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]),
        3 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => {
            flags.extend([
                ("source-phase-aqueous", 1.0),
                ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ]);
            if index != 6 {
                flags.extend([
                    ("source-individual-ion", 1.0),
                    ("source-ionic-reference-hydrogen-zero", 1.0),
                ]);
            }
        }
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}

fn without_axis(at: &QueryPoint, omitted: &str) -> QueryPoint {
    let mut result = QueryPoint::new();
    for (axis, value) in at.axes() {
        if axis != omitted {
            result = result
                .with_quantity(axis.clone(), at.axis_quantities()[axis], *value)
                .unwrap();
        }
    }
    result
}

fn materials() -> &'static Materials {
    static MATERIALS: OnceLock<Materials> = OnceLock::new();
    MATERIALS.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-lime-equilibrium-{}-{nonce}",
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
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}
fn species_index(id: &SpeciesId) -> usize {
    SPECIES.iter().position(|s| *s == id.as_str()).unwrap()
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

#[test]
fn original_calcium_rows_keep_units_ionic_entropy_phase_and_source_scope() {
    let m = materials();
    // Independent literal transcription: NBS Table 94:Ca, Calcium (Prepared 1968).
    for (i, literals) in [
        (0, [0.0560794, -635_090.0, -604_030.0, 39.75]),
        (1, [0.0740948, -986_090.0, -898_490.0, 83.39]),
        (2, [0.0400800, -542_830.0, -553_580.0, -53.1]),
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
            close(
                value(i, name, dims),
                expected,
                1e-12 * expected.abs().max(1.0),
            );
        }
        let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
        for axis in m.points[i].axes().keys() {
            assert!(
                m.packs[i]
                    .claims()
                    .query_typed(
                        &key,
                        &without_axis(&m.points[i], axis),
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .is_err(),
                "missing {axis} must refuse for {}",
                PACKS[i]
            );
        }
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
        assert!(
            m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new("molar-mass", Dims::NONE),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
        // No source row names a solid polymorph, and the admitted solid query needs none.
        assert!(!m.points[i].axes().contains_key("source-phase-portlandite"));
    }
}

#[test]
fn actual_lime_sources_close_hydration_neutralization_and_reverse_precipitation() {
    let m = materials();
    for (recipe, expected_h, expected_g, expected_s) in [
        (MineralReaction::LimeHydration, -65_170.0, -57_331.0, -26.27),
        (
            MineralReaction::HydratedLimeNeutralization,
            -128_400.0,
            -129_348.0,
            3.33,
        ),
        (
            MineralReaction::CalciumHydroxidePrecipitation,
            128_400.0,
            129_348.0,
            -3.33,
        ),
    ] {
        let reaction = recipe.build().unwrap();
        let enthalpies = reaction
            .terms()
            .iter()
            .map(|term| {
                let i = species_index(&term.species);
                FormationEnthalpySource {
                    species: &m.ids[i],
                    pack: &m.packs[i],
                    at: &m.points[i],
                    reference_axis: REFERENCE,
                }
            })
            .collect::<Vec<_>>();
        let indices = reaction
            .terms()
            .iter()
            .map(|term| species_index(&term.species))
            .collect::<Vec<_>>();
        let caloric = reaction.standard_equilibrium(&enthalpies).unwrap();
        let direct = reaction
            .standard_gibbs_change(&gibbs_sources(&indices))
            .unwrap();
        close(caloric.enthalpy().molar().value(), expected_h, 1e-8);
        close(caloric.molar_entropy().value(), expected_s, 1e-10);
        close(
            caloric.molar_gibbs().value(),
            expected_h - 298.15 * expected_s,
            1e-8,
        );
        close(direct.molar().value(), expected_g, 1e-8);
        assert_eq!(direct.reference_axis(), REFERENCE);
        for term in direct.terms() {
            m.packs[species_index(&term.term.species)]
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        // A caller cannot relabel a retained NBS source as the USGS reference.
        let mut incompatible = gibbs_sources(&indices);
        incompatible[0].reference_axis = "source-elemental-reference-usgs-b2131";
        assert!(reaction.standard_gibbs_change(&incompatible).is_err());
    }
    let reaction = MineralReaction::LimeHydration.build().unwrap();
    let extent = reaction
        .extent_from_reacted_mass(
            &m.ids[0],
            Mass::new(0.00560794),
            MolarMass::new(value(0, "molar-mass", MolarMass::DIMS)),
        )
        .unwrap();
    close(extent.value(), 0.1, 1e-15);
    let inventory = [
        (m.ids[0].clone(), Amount::new(0.25)),
        (m.ids[3].clone(), Amount::new(0.20)),
        (m.ids[1].clone(), Amount::new(0.0)),
    ];
    let result = reaction.advance_inventory(&inventory, extent).unwrap();
    for (id, expected) in [(&m.ids[0], 0.15), (&m.ids[3], 0.10), (&m.ids[1], 0.10)] {
        close(
            result.iter().find(|(s, _)| s == id).unwrap().1.value(),
            expected,
            1e-15,
        );
    }
    assert!(
        reaction
            .advance_inventory(&inventory, Amount::new(0.3))
            .is_err()
    );
}

fn lime_model() -> HydrolysisPrecipitation {
    let ids = &materials().ids;
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build().unwrap(),
        ids[3].clone(),
        ids[4].clone(),
        ids[5].clone(),
        &gibbs_sources(&[3, 4, 5]),
    )
    .unwrap();
    let hcn = AcidDissociation::from_sources(
        MineralReaction::HydrogenCyanideDissociation
            .build()
            .unwrap(),
        ids[6].clone(),
        ids[7].clone(),
        ids[4].clone(),
        &gibbs_sources(&[6, 4, 7]),
    )
    .unwrap();
    let aqueous =
        IdealAcidBaseSystem::new(water, vec![AcidFamily::new(FAMILY, vec![hcn]).unwrap()]).unwrap();
    HydrolysisPrecipitation::from_sources(
        aqueous,
        MineralReaction::CalciumHydroxidePrecipitation
            .build()
            .unwrap(),
        ids[2].clone(),
        ids[1].clone(),
        &gibbs_sources(&[2, 3, 1, 4]),
    )
    .unwrap()
}

#[test]
fn finite_lime_and_cyanide_share_ph_and_preserve_calcium_charge_and_saturation() {
    let model = lime_model();
    let totals = [(FAMILY.to_owned(), Amount::new(0.01))];
    let sodium = [SpectatorIon {
        species: SpeciesId::new("Na[+](aq)").unwrap(),
        charge: 1,
        amount: Amount::new(0.01),
    }];
    // The source molal reference uses n_ref=1 mol for the explicit 1 kg solvent.
    // These totals are authored scenarios, not an assertion about a mine reagent.
    let mut results = Vec::new();
    for total in [0.0, 0.001, 0.1] {
        let result = model
            .speciate(
                Amount::new(total),
                &totals,
                &sodium,
                Amount::new(1.0),
                PrecipitationOptions::default(),
                || false,
            )
            .unwrap();
        close(
            result.dissolved_amount.value() + result.solid_amount.value(),
            total,
            1e-14,
        );
        close(result.component_balance_residual.value(), 0.0, 1e-14);
        let cyanide = &result.aqueous.families[0];
        close(
            cyanide.forms.iter().map(|f| f.amount.value()).sum(),
            0.01,
            1e-14,
        );
        let cn = cyanide
            .forms
            .iter()
            .find(|f| f.species.as_str() == "CN[-](aq)")
            .unwrap()
            .amount
            .value();
        let charge =
            2.0 * result.dissolved_amount.value() + 0.01 + result.aqueous.proton_amount.value()
                - cn
                - result.aqueous.hydroxide_amount.value();
        close(charge, 0.0, 2e-13);
        close(
            result.aqueous.charge_balance_residual.value(),
            charge,
            2e-15,
        );
        assert_eq!(result.model().standard().reference_axis(), REFERENCE);
        for term in result.model().standard().terms() {
            materials().packs[species_index(&term.term.species)]
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
        results.push((result, cn));
    }
    assert_eq!(results[0].0.status, PrecipitationStatus::EmptyComponent);
    assert_eq!(results[1].0.status, PrecipitationStatus::SolidAbsent);
    assert_eq!(results[1].0.solid_amount, Amount::new(0.0));
    close(results[1].0.dissolved_amount.value(), 0.001, 1e-15);
    assert!(results[1].0.log_saturation_ratio.unwrap().value() < 0.0);
    assert_eq!(results[2].0.status, PrecipitationStatus::SolidPresent);
    assert!(results[2].0.solid_amount.value() > 0.0);
    let present = &results[2].0;
    // Independent dissolution direction: Ca(OH)2 -> Ca2+ + 2 OH-.
    // Literal NBS formation G values give +30422 J/mol, with no fitted Ksp.
    let ksp = (-30_422.0 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15)).exp();
    let iap = present.dissolved_amount.value() * present.aqueous.hydroxide_amount.value().powi(2);
    close(iap / ksp, 1.0, 2e-9);
    close(present.log_saturation_ratio.unwrap().value(), 0.0, 1e-9);
    assert!(results[0].0.aqueous.ph.value() < results[1].0.aqueous.ph.value());
    assert!(results[1].0.aqueous.ph.value() < results[2].0.aqueous.ph.value());
    assert!(results[0].1 < results[1].1 && results[1].1 < results[2].1);
    // No second ownership of calcium, and a cancelled solve returns no result.
    let duplicate = [SpectatorIon {
        species: materials().ids[2].clone(),
        charge: 2,
        amount: Amount::new(0.001),
    }];
    assert!(
        model
            .speciate(
                Amount::new(0.1),
                &totals,
                &duplicate,
                Amount::new(1.0),
                PrecipitationOptions::default(),
                || false,
            )
            .is_err()
    );
    assert!(
        model
            .speciate(
                Amount::new(0.1),
                &totals,
                &sodium,
                Amount::new(1.0),
                PrecipitationOptions::default(),
                || true,
            )
            .is_err()
    );
}
