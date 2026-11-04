//! Actual NBS arsenic source records through conserved chemistry and shared pH.
//! No measured mine-water pH, reaction rate, removal or fixation is claimed.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        FormationEnthalpySource, ReactionBalanceError,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                AcidBaseError, AcidDissociation, AcidFamily, IdealAcidBaseSystem,
                SpeciationOptions, SpectatorIon, WaterAutoprotolysis,
            },
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-arsenate";
const ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const PACKS: [&str; 10] = [
    "arsenic-acid-aqueous-nbs1989",
    "dihydrogen-arsenate-ion-aqueous-nbs1989",
    "hydrogen-arsenate-ion-aqueous-nbs1989",
    "arsenate-ion-aqueous-nbs1989",
    "arsenous-acid-aqueous-nbs1989",
    "dihydrogen-arsenite-ion-aqueous-nbs1989",
    "hydrogen-ion-aqueous-nbs1989",
    "hydroxide-ion-aqueous-nbs1989",
    "water-liquid-nbs1989-thermochemistry",
    "hydrogen-peroxide-aqueous-nbs1989",
];
const SPECIES: [&str; 10] = [
    "H3AsO4(aq)",
    "H2AsO4[-](aq)",
    "HAsO4[2-](aq)",
    "AsO4[3-](aq)",
    "H3AsO3(aq)",
    "H2AsO3[-](aq)",
    "H+(aq)",
    "OH[-](aq)",
    "H2O(l)",
    "H2O2(aq)",
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

fn point(index: usize, arsenic_context: bool) -> QueryPoint {
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
    if index == 8 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
        ]);
        if ![0, 4, 9].contains(&index) {
            flags.extend([
                ("source-individual-ion", 1.0),
                ("source-ionic-reference-hydrogen-zero", 1.0),
            ]);
        }
    }
    if arsenic_context {
        if index < 4 {
            flags.push(("source-arsenic-oxidation-state-v", 1.0));
        } else if index < 6 {
            flags.extend([
                ("source-arsenic-oxidation-state-iii", 1.0),
                ("source-arsenite-component-plus-water", 1.0),
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
    static DATA: OnceLock<Materials> = OnceLock::new();
    DATA.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out =
            std::env::temp_dir().join(format!("frankensim-arsenic-{}-{nonce}", std::process::id()));
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
            points: (0..PACKS.len()).map(|i| point(i, true)).collect(),
            ids: SPECIES
                .iter()
                .map(|s| SpeciesId::new(*s).unwrap())
                .collect(),
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
            reference_axis: REFERENCE,
        })
        .collect()
}

fn step(recipe: MineralReaction, acid: usize, base: usize) -> AcidDissociation {
    let m = materials();
    AcidDissociation::from_sources(
        recipe.build().unwrap(),
        m.ids[acid].clone(),
        m.ids[base].clone(),
        m.ids[6].clone(),
        &sources(&[acid, 6, base]),
    )
    .unwrap()
}

fn system() -> IdealAcidBaseSystem {
    let m = materials();
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build().unwrap(),
        m.ids[8].clone(),
        m.ids[6].clone(),
        m.ids[7].clone(),
        &sources(&[8, 6, 7]),
    )
    .unwrap();
    let family = AcidFamily::new(
        FAMILY,
        vec![
            step(MineralReaction::ArsenicAcidDissociation, 0, 1),
            step(MineralReaction::DihydrogenArsenateDissociation, 1, 2),
            step(MineralReaction::HydrogenArsenateDissociation, 2, 3),
        ],
    )
    .unwrap();
    IdealAcidBaseSystem::new(water, vec![family]).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; tolerance {tolerance}"
    );
}

#[test]
fn g0_original_arsenic_rows_retain_units_states_uncertainty_and_receipts() {
    // Literal ao rows in corrected NBS Table 20:As, not the different aq/cr rows.
    let expected = [
        [0.1419432, -902500.0, -766000.0, 184.0],
        [0.1409352, -909560.0, -753170.0, 117.0],
        [0.1399272, -906340.0, -714600.0, -1.7],
        [0.1389192, -888140.0, -648410.0, -162.8],
        [0.1259438, -742200.0, -639800.0, 195.0],
        [0.1249358, -714790.0, -587130.0, 110.5],
    ];
    let m = materials();
    for (i, values) in expected.into_iter().enumerate() {
        assert_eq!(m.packs[i].claims().claim_count(), 4);
        let entropy = if [0, 4].contains(&i) {
            "standard-molar-entropy"
        } else {
            "conventional-standard-molar-entropy"
        };
        for ((property, dims), expected) in [
            ("molar-mass", MolarMass::DIMS),
            ("standard-molar-formation-enthalpy", ENERGY),
            ("standard-molar-formation-gibbs-energy", ENERGY),
            (entropy, ENTROPY),
        ]
        .into_iter()
        .zip(values)
        {
            let answer = m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new(property, dims),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected, 1e-9);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            m.packs[i].claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
}

#[test]
fn g1_four_proton_steps_use_original_gibbs_and_enthalpy_without_refitting() {
    for (recipe, acid, base, dg, dh, pka) in [
        (
            MineralReaction::ArsenicAcidDissociation,
            0,
            1,
            12830.0,
            -7060.0,
            2.2477187527984053,
        ),
        (
            MineralReaction::DihydrogenArsenateDissociation,
            1,
            2,
            38570.0,
            3220.0,
            6.757171652021395,
        ),
        (
            MineralReaction::HydrogenArsenateDissociation,
            2,
            3,
            66190.0,
            18200.0,
            11.595986301459584,
        ),
        (
            MineralReaction::ArsenousAcidDissociation,
            4,
            5,
            52670.0,
            27410.0,
            9.227384778635386,
        ),
    ] {
        let acid_step = step(recipe, acid, base);
        close(acid_step.standard().molar().value(), dg, 1e-8);
        close(
            -acid_step.standard().log_equilibrium_constant().value() / core::f64::consts::LN_10,
            pka,
            1e-11,
        );
        let sources = sources(&[acid, 6, base]);
        let heat = sources
            .iter()
            .map(|s| FormationEnthalpySource {
                species: s.species,
                pack: s.pack,
                at: s.at,
                reference_axis: s.reference_axis,
            })
            .collect::<Vec<_>>();
        close(
            acid_step
                .reaction()
                .standard_enthalpy_change(&heat)
                .unwrap()
                .molar()
                .value(),
            dh,
            1e-8,
        );
    }
}

#[test]
fn g1_g3_triprotic_partition_matches_independent_polynomial_and_charge_balance() {
    let system = system();
    // Closed-form triprotic fractions, evaluated from literal source differences.
    // Choose spectator inventory independently to make each prescribed h a root;
    // the production owner must still solve its charge equation to recover h.
    let rt = 8.31446261815324 * 298.15;
    let k = [12830.0_f64, 38570.0, 66190.0].map(|dg| (-dg / rt).exp());
    let kw = (-79885.0_f64 / rt).exp();
    for ph in [2.0_f64, 4.0, 7.0, 9.0, 12.0] {
        let h = 10.0_f64.powf(-ph);
        let weights = [h * h * h, k[0] * h * h, k[0] * k[1] * h, k[0] * k[1] * k[2]];
        let denominator: f64 = weights.iter().sum();
        let alpha = weights.map(|w| w / denominator);
        for scale in [1.0, 10.0] {
            let total = 0.02 * scale;
            let basis = 10.0 * scale;
            let fixed_charge =
                total * (alpha[1] + 2.0 * alpha[2] + 3.0 * alpha[3]) + basis * (kw / h - h);
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
            let answer = system
                .speciate(
                    &[(FAMILY.to_owned(), Amount::new(total))],
                    &[spectator],
                    Amount::new(basis),
                    SpeciationOptions {
                        ph_tolerance: Dimensionless::new(1e-11),
                        absolute_charge_tolerance: Amount::new(1e-14),
                        relative_charge_tolerance: Dimensionless::new(1e-12),
                        ..SpeciationOptions::default()
                    },
                    || false,
                )
                .unwrap();
            close(answer.ph.value(), ph, 2e-10);
            let mut charge =
                fixed_charge + answer.proton_amount.value() - answer.hydroxide_amount.value();
            let mut sum = 0.0;
            for (i, fraction) in alpha.into_iter().enumerate() {
                let form = answer.families[0]
                    .forms
                    .iter()
                    .find(|f| f.species.as_str() == SPECIES[i])
                    .unwrap();
                close(form.fraction.value(), fraction, 2e-10);
                close(form.amount.value(), total * fraction, 2e-10 * total);
                charge -= i as f64 * form.amount.value();
                sum += form.amount.value();
            }
            close(sum, total, 2e-15 * scale);
            close(charge, 0.0, 2e-12 * scale);
            assert!(answer.max_log_mass_action_residual.value() < 1e-9);
        }
    }
}

#[test]
fn g1_peroxide_oxidation_has_sourced_heat_and_cannot_borrow_oxidant() {
    let reaction = MineralReaction::ArsenousAcidPeroxideOxidation
        .build()
        .unwrap();
    let sources = sources(&[4, 9, 0, 8]);
    let gibbs = reaction.standard_gibbs_change(&sources).unwrap();
    close(gibbs.molar().value(), -229299.0, 1e-8);
    let heat = sources
        .iter()
        .map(|s| FormationEnthalpySource {
            species: s.species,
            pack: s.pack,
            at: s.at,
            reference_axis: s.reference_axis,
        })
        .collect::<Vec<_>>();
    close(
        reaction
            .standard_enthalpy_change(&heat)
            .unwrap()
            .molar()
            .value(),
        -254960.0,
        1e-8,
    );
    let m = materials();
    let inventory = vec![
        (m.ids[4].clone(), Amount::new(0.2)),
        (m.ids[9].clone(), Amount::new(0.05)),
        (m.ids[0].clone(), Amount::new(0.0)),
        (m.ids[8].clone(), Amount::new(100.0)),
    ];
    close(
        reaction.limiting_extent(&inventory[..2]).unwrap().value(),
        0.05,
        1e-16,
    );
    let out = reaction
        .advance_inventory(&inventory, Amount::new(0.05))
        .unwrap();
    let amount = |index: usize| {
        out.iter()
            .find(|(id, _)| *id == m.ids[index])
            .unwrap()
            .1
            .value()
    };
    close(amount(9), 0.0, 1e-16);
    close(amount(4), 0.15, 1e-16);
    close(amount(0), 0.05, 1e-16);
    close(amount(8), 100.05, 1e-13);
    close(amount(4) + amount(0), 0.2, 1e-16);
    let mass = |stock: &[(SpeciesId, Amount)]| {
        stock
            .iter()
            .map(|(id, amount)| {
                let i = m.ids.iter().position(|known| known == id).unwrap();
                let molar_mass = m.packs[i]
                    .claims()
                    .query_typed(
                        &PropertyKey::new("molar-mass", MolarMass::DIMS),
                        &m.points[i],
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .unwrap()
                    .evidence
                    .value
                    .value;
                amount.value() * molar_mass
            })
            .sum::<f64>()
    };
    close(mass(&out), mass(&inventory), 1e-13);
    assert!(matches!(
        reaction.advance_inventory(&inventory, Amount::new(0.06)),
        Err(ReactionBalanceError::ExtentExceedsInventory { .. })
    ));
    assert_eq!(inventory[1].1.value(), 0.05);
}

#[test]
fn g0_missing_arsenic_context_temperature_extrapolation_and_invented_properties_refuse() {
    let m = materials();
    for i in 0..6 {
        let key = PropertyKey::new("standard-molar-formation-gibbs-energy", ENERGY);
        assert!(
            m.packs[i]
                .claims()
                .query_typed(&key, &point(i, false), SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
        if i >= 4 {
            // Each convention is independently required. Keeping just one must
            // not accidentally admit the source component as another species.
            for retained_axis in [
                "source-arsenic-oxidation-state-iii",
                "source-arsenite-component-plus-water",
            ] {
                let incomplete = point(i, false)
                    .with_quantity(retained_axis, QuantitySpec::dimensional(Dims::NONE), 1.0)
                    .unwrap();
                assert!(
                    m.packs[i]
                        .claims()
                        .query_typed(&key, &incomplete, SelectionPolicy::SingleClaimOnly,)
                        .is_err()
                );
            }
        }
        let warm = m.points[i]
            .clone()
            .with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                310.0,
            )
            .unwrap();
        assert!(
            m.packs[i]
                .claims()
                .query_typed(&key, &warm, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
        assert!(
            m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new("standard-molar-heat-capacity", ENTROPY),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
        assert!(
            m.packs[i]
                .claims()
                .query_typed(
                    &PropertyKey::new("standard-molar-formation-gibbs-energy", Dims::NONE),
                    &m.points[i],
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
    assert!(matches!(
        system().speciate(
            &[(FAMILY.to_owned(), Amount::new(0.01))],
            &[],
            Amount::new(10.0),
            SpeciationOptions::default(),
            || true,
        ),
        Err(AcidBaseError::Cancelled)
    ));
}
