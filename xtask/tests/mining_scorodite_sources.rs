//! Consistent source chemistry through finite Fe/As/solid and shared-pH closure.
//! Analytical fixtures test the declared ideal selected-phase model, not a
//! measured removal yield, stable-phase assemblage or mine-water calibration.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon,
                WaterAutoprotolysis,
                precipitation::{
                    PrecipitationOptions, PrecipitationStatus, salt::SaltPrecipitation,
                },
            },
        },
    },
};
use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

const REFERENCE: &str = "source-elemental-reference-nmk2014";
const NBS_REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-arsenate";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const PACKS: [&str; 10] = [
    "arsenic-acid-aqueous-nmk2014",
    "dihydrogen-arsenate-ion-aqueous-nmk2014",
    "hydrogen-arsenate-ion-aqueous-nmk2014",
    "arsenate-ion-aqueous-nmk2014",
    "iron-iii-ion-aqueous-nmk2014",
    "scorodite-solid-nmk2014",
    "water-liquid-nmk2014",
    "hydrogen-ion-aqueous-nmk2014",
    "hydroxide-ion-aqueous-nmk2014",
    "iron-iii-ion-aqueous-nbs1989",
];
const SPECIES: [&str; 10] = [
    "H3AsO4(aq)",
    "H2AsO4[-](aq)",
    "HAsO4[2-](aq)",
    "AsO4[3-](aq)",
    "Fe[3+](aq)",
    "FeAsO4.2H2O(s)",
    "H2O(l)",
    "H+(aq)",
    "OH[-](aq)",
    "Fe[3+](aq)",
];
struct Materials {
    packs: Vec<NormalizedPack>,
    ids: Vec<SpeciesId>,
    points: Vec<QueryPoint>,
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
        (if index == 9 { NBS_REFERENCE } else { REFERENCE }, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    match index {
        5 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-pure-species", 1.0),
            ("source-phase-scorodite", 1.0),
        ]),
        6 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => {
            flags.extend([
                ("source-phase-aqueous", 1.0),
                ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ]);
            if index != 0 {
                flags.extend([
                    ("source-individual-ion", 1.0),
                    ("source-ionic-reference-hydrogen-zero", 1.0),
                ]);
            }
        }
    }
    if index < 4 || index == 5 {
        flags.push(("source-arsenic-oxidation-state-v", 1.0));
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}

fn materials() -> &'static Materials {
    static DATA: OnceLock<Materials> = OnceLock::new();
    DATA.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-scorodite-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        let packs = PACKS
            .iter()
            .map(|name| {
                let output = out.join(format!("{name}.fsmatpk"));
                let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                    .args(["matdb-pack", "--manifest"])
                    .arg(workspace().join(format!("data/matdb/seed-v1/{name}/manifest.tsv")))
                    .arg("--out")
                    .arg(&output)
                    .env("CARGO_WORKSPACE_DIR", workspace())
                    .current_dir(workspace())
                    .output()
                    .unwrap();
                assert!(
                    run.status.success(),
                    "{name}: {}",
                    String::from_utf8_lossy(&run.stderr)
                );
                let bytes = std::fs::read(output).unwrap();
                let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
                let pack =
                    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
                assert_eq!(pack.pack_id(), *name);
                assert_eq!(pack.to_bytes(), bytes);
                pack
            })
            .collect();
        Materials {
            packs,
            ids: SPECIES
                .iter()
                .map(|s| SpeciesId::new(*s).unwrap())
                .collect(),
            points: (0..PACKS.len()).map(point).collect(),
        }
    })
}

fn sources(indices: &[usize]) -> Vec<FormationGibbsSource<'static>> {
    let m = materials();
    indices
        .iter()
        .map(|&i| FormationGibbsSource {
            species: &m.ids[i],
            pack: &m.packs[i],
            at: &m.points[i],
            reference_axis: if i == 9 { NBS_REFERENCE } else { REFERENCE },
        })
        .collect()
}

fn aqueous() -> IdealAcidBaseSystem {
    let m = materials();
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build().unwrap(),
        m.ids[6].clone(),
        m.ids[7].clone(),
        m.ids[8].clone(),
        &sources(&[6, 7, 8]),
    )
    .unwrap();
    let steps = [
        MineralReaction::ArsenicAcidDissociation,
        MineralReaction::DihydrogenArsenateDissociation,
        MineralReaction::HydrogenArsenateDissociation,
    ]
    .into_iter()
    .enumerate()
    .map(|(i, recipe)| {
        AcidDissociation::from_sources(
            recipe.build().unwrap(),
            m.ids[i].clone(),
            m.ids[i + 1].clone(),
            m.ids[7].clone(),
            &sources(&[i, 7, i + 1]),
        )
        .unwrap()
    })
    .collect();
    IdealAcidBaseSystem::new(water, vec![AcidFamily::new(FAMILY, steps).unwrap()]).unwrap()
}

fn model() -> SaltPrecipitation {
    let m = materials();
    SaltPrecipitation::from_sources(
        aqueous(),
        MineralReaction::ScoroditePrecipitation.build().unwrap(),
        m.ids[4].clone(),
        m.ids[5].clone(),
        FAMILY,
        &sources(&[4, 3, 6, 5]),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, atol: f64, rtol: f64) {
    assert!(
        (actual - expected).abs() <= atol + rtol * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_consistent_evaluated_gibbs_rows_and_quoted_uncertainties_survive_the_compiler() {
    // NMK2014 Tables8/9 and their CODATA-based auxiliaries, not legacy NBS iron.
    let gibbs = [
        -766420.0, -753575.0, -713727.0, -647618.0, -16280.0, -1284800.0, -237140.0, 0.0, -157220.0,
    ];
    let quoted = [
        Some(350.0),
        Some(350.0),
        Some(350.0),
        Some(1500.0),
        Some(1100.0),
        Some(2900.0),
        Some(41.0),
        None,
        Some(72.0),
    ];
    let m = materials();
    for (i, expected) in gibbs.into_iter().enumerate() {
        let answer = m.packs[i]
            .claims()
            .query_typed(
                &PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY),
                &m.points[i],
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        close(answer.evidence.value.value, expected, 1e-9, 0.0);
        let uncertainty = match i {
            6 => UncertaintyModel::HalfWidth {
                half_width: 41.0,
                confidence: 0.95,
            },
            8 => UncertaintyModel::HalfWidth {
                half_width: 72.0,
                confidence: 0.95,
            },
            _ => UncertaintyModel::Unstated,
        };
        assert_eq!(answer.evidence.value.uncertainty, uncertainty);
        assert!(answer.receipt.observation_backed);
        m.packs[i].claims().verify_receipt(&answer.receipt).unwrap();
        if let Some(error) = quoted[i] {
            let reported = m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new("source-reported-formation-gibbs-uncertainty", ENERGY),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(reported.evidence.value.value, error, 1e-10, 0.0);
            assert_eq!(
                reported.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            m.packs[i]
                .claims()
                .verify_receipt(&reported.receipt)
                .unwrap();
        }
    }
}

#[test]
fn g1_source_dissolution_constant_and_proton_steps_match_literal_thermodynamics() {
    let reaction = MineralReaction::ScoroditePrecipitation.build().unwrap();
    let standard = reaction
        .standard_gibbs_change(&sources(&[4, 3, 6, 5]))
        .unwrap();
    close(standard.molar().value(), -146622.0, 1e-8, 0.0);
    let log_ksp = -standard.log_equilibrium_constant().value() / core::f64::consts::LN_10;
    close(log_ksp, -25.68706305321963, 1e-12, 0.0);
    // The published table rounds separately; do not alter formation values to fit it.
    close(log_ksp, -25.68, 0.02, 0.0);
    let aqueous = aqueous();
    for (step, expected) in aqueous.families()[0].steps().iter().zip([
        2.2503466391033142,
        6.9810675651996,
        11.58179571541308,
    ]) {
        close(
            -step.standard().log_equilibrium_constant().value() / core::f64::consts::LN_10,
            expected,
            1e-12,
            0.0,
        );
    }
}

#[test]
fn g1_g3_real_sources_recover_independent_fe_as_saturation_and_charge_roots() {
    let model = model();
    let rt = 8.31446261815324 * 298.15;
    let k = [12845.0_f64, 39848.0, 66109.0].map(|dg| (-dg / rt).exp());
    let kw = (-79920.0_f64 / rt).exp();
    let ksp = (-146622.0_f64 / rt).exp();
    // Independent polynomial alpha and quadratic activity-product solution.
    // The production owner receives component/strong-ion totals, not target pH.
    for ph in [2.0_f64, 3.0, 4.0] {
        let h = 10.0_f64.powf(-ph);
        let weights = [h * h * h, k[0] * h * h, k[0] * k[1] * h, k[0] * k[1] * k[2]];
        let sum: f64 = weights.iter().sum();
        let alpha = weights.map(|w| w / sum);
        for ratio in [0.1_f64, 1.0, 10.0] {
            for scale in [1.0, 1000.0] {
                let basis = 10.0 * scale;
                let geometric_mean = basis * (ksp / alpha[3]).sqrt();
                let fe = geometric_mean * ratio.sqrt();
                let arsenic = geometric_mean / ratio.sqrt();
                let solid = 0.002 * scale;
                let fixed_charge = -(3.0 * fe
                    - arsenic * (alpha[1] + 2.0 * alpha[2] + 3.0 * alpha[3])
                    + basis * (h - kw / h));
                let spectator = SpectatorIon {
                    species: SpeciesId::new(if fixed_charge >= 0.0 {
                        "Na[+](aq)"
                    } else {
                        "Cl[-](aq)"
                    })
                    .unwrap(),
                    charge: if fixed_charge >= 0.0 { 1 } else { -1 },
                    amount: Amount::new(fixed_charge.abs()),
                };
                let answer = model
                    .speciate(
                        Amount::new(fe + solid),
                        &[(FAMILY.into(), Amount::new(arsenic + solid))],
                        &[spectator],
                        Amount::new(basis),
                        PrecipitationOptions::default(),
                        || false,
                    )
                    .unwrap();
                assert_eq!(answer.status, PrecipitationStatus::SolidPresent);
                close(answer.aqueous.ph.value(), ph, 2e-7, 0.0);
                close(answer.dissolved_cation.value(), fe, 2e-14 * scale, 2e-7);
                close(
                    answer.dissolved_family.value(),
                    arsenic,
                    2e-14 * scale,
                    2e-7,
                );
                close(answer.solid_amount.value(), solid, 2e-14 * scale, 2e-7);
                let forms = &answer.aqueous.families[0].forms;
                for (form, fraction) in forms.iter().zip(alpha) {
                    close(form.fraction.value(), fraction, 1e-15, 2e-7);
                }
                let dissolved_as: f64 = forms.iter().map(|f| f.amount.value()).sum();
                close(
                    dissolved_as + answer.solid_amount.value(),
                    arsenic + solid,
                    2e-14 * scale,
                    0.0,
                );
                close(
                    answer.dissolved_cation.value() + answer.solid_amount.value(),
                    fe + solid,
                    2e-14 * scale,
                    0.0,
                );
                let iap =
                    (answer.dissolved_cation.value() / basis) * (forms[3].amount.value() / basis);
                close((iap / ksp).ln(), 0.0, 2e-9, 0.0);
                let charge = fixed_charge + 3.0 * answer.dissolved_cation.value()
                    - forms
                        .iter()
                        .enumerate()
                        .map(|(i, f)| i as f64 * f.amount.value())
                        .sum::<f64>()
                    + answer.aqueous.proton_amount.value()
                    - answer.aqueous.hydroxide_amount.value();
                close(charge, 0.0, 2e-12 * scale, 0.0);
            }
        }
    }
}

#[test]
fn g0_g1_undersaturation_and_incompatible_legacy_reference_do_not_create_a_solid() {
    let model = model();
    let answer = model
        .speciate(
            Amount::new(1e-12),
            &[(FAMILY.into(), Amount::new(1e-12))],
            &[SpectatorIon {
                species: SpeciesId::new("Cl[-](aq)").unwrap(),
                charge: -1,
                amount: Amount::new(0.1),
            }],
            Amount::new(10.0),
            PrecipitationOptions::default(),
            || false,
        )
        .unwrap();
    assert_eq!(answer.status, PrecipitationStatus::SolidAbsent);
    assert_eq!(answer.solid_amount.value(), 0.0);
    assert_eq!(answer.dissolved_cation.value(), 1e-12);
    assert_eq!(answer.dissolved_family.value(), 1e-12);
    assert!(answer.log_saturation_ratio.unwrap().value() < 0.0);
    let m = materials();
    let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
    // Both original Fe3+ sources are individually queryable but are different evaluations.
    close(
        m.packs[9]
            .claims()
            .query_typed(&key, &m.points[9], SelectionPolicy::SingleClaimOnly)
            .unwrap()
            .evidence
            .value
            .value,
        -4700.0,
        1e-10,
        0.0,
    );
    assert!(
        SaltPrecipitation::from_sources(
            aqueous(),
            MineralReaction::ScoroditePrecipitation.build().unwrap(),
            m.ids[4].clone(),
            m.ids[5].clone(),
            FAMILY,
            &sources(&[9, 3, 6, 5])
        )
        .is_err()
    );
    let wrong_phase = m.points[5]
        .clone()
        .with_quantity(
            "source-phase-scorodite",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    assert!(
        m.packs[5]
            .claims()
            .query_typed(&key, &wrong_phase, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    let warm = m.points[5]
        .clone()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            350.0,
        )
        .unwrap();
    assert!(
        m.packs[5]
            .claims()
            .query_typed(&key, &warm, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
}
