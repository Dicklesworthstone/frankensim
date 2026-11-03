//! Actual USGS sources: simultaneous pH and finite goethite partition.
//! Authored solution inventories are model scenarios, not field calibration.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, QueryPoint};
use fs_qty::{
    Amount, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{
                IdealAcidBaseSystem, SpectatorIon, WaterAutoprotolysis,
                precipitation::{
                    HydrolysisPrecipitation, PrecipitationOptions, PrecipitationSpeciation,
                    PrecipitationStatus,
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

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const PACKS: [&str; 5] = [
    "iron-iii-ion-aqueous-usgs-b2131",
    "goethite-solid-usgs-b2131-thermochemistry",
    "hydrogen-ion-aqueous-usgs-b2131",
    "hydroxide-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
];
const SPECIES: [&str; 5] = ["Fe[3+](aq)", "FeO(OH)(s)", "H+(aq)", "OH[-](aq)", "H2O(l)"];
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
        1 => flags.extend([
            ("source-phase-solid", 1.0),
            ("source-phase-goethite", 1.0),
            ("source-pure-species", 1.0),
        ]),
        4 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-individual-ion", 1.0),
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
fn materials() -> &'static Materials {
    static MATERIALS: OnceLock<Materials> = OnceLock::new();
    MATERIALS.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-precipitation-{}-{nonce}",
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
                assert_eq!(pack.pack_id(), *name);
                pack
            })
            .collect();
        Materials {
            packs,
            points: (0..5).map(point).collect(),
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
fn model() -> HydrolysisPrecipitation {
    let m = materials();
    let water = WaterAutoprotolysis::from_sources(
        MineralReaction::WaterDissociation.build().unwrap(),
        m.ids[4].clone(),
        m.ids[2].clone(),
        m.ids[3].clone(),
        &sources(&[4, 2, 3]),
    )
    .unwrap();
    let aqueous = IdealAcidBaseSystem::new(water, vec![]).unwrap();
    HydrolysisPrecipitation::from_sources(
        aqueous,
        MineralReaction::FerricHydrolysisGoethite.build().unwrap(),
        m.ids[0].clone(),
        m.ids[1].clone(),
        &sources(&[0, 4, 1, 2]),
    )
    .unwrap()
}
fn fixed_ions(chloride: f64, sodium: f64) -> Vec<SpectatorIon> {
    [("Cl[-](aq)", -1, chloride), ("Na[+](aq)", 1, sodium)]
        .into_iter()
        .map(|(name, charge, amount)| SpectatorIon {
            species: SpeciesId::new(name).unwrap(),
            charge,
            amount: Amount::new(amount),
        })
        .collect()
}
fn solve(
    model: &HydrolysisPrecipitation,
    iron: f64,
    chloride: f64,
    sodium: f64,
) -> PrecipitationSpeciation {
    model
        .speciate(
            Amount::new(iron),
            &[],
            &fixed_ions(chloride, sodium),
            Amount::new(10.0),
            PrecipitationOptions::default(),
            || false,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}
fn relative(actual: f64, expected: f64) {
    close(actual, expected, expected.abs() * 3e-8);
}
fn check_closure(result: &PrecipitationSpeciation, total: f64, chloride: f64, sodium: f64) {
    assert!(result.dissolved_amount.value() >= 0.0 && result.solid_amount.value() >= 0.0);
    close(
        result.dissolved_amount.value() + result.solid_amount.value(),
        total,
        3e-18,
    );
    close(result.component_balance_residual.value(), 0.0, 3e-18);
    let aqueous = &result.aqueous;
    // Independent complete aqueous charge: the neutral solid contributes zero.
    let charge = aqueous.proton_amount.value() - aqueous.hydroxide_amount.value()
        + 3.0 * result.dissolved_amount.value()
        + sodium
        - chloride;
    let unsigned_charge = aqueous.proton_amount.value()
        + aqueous.hydroxide_amount.value()
        + 3.0 * result.dissolved_amount.value()
        + sodium
        + chloride;
    close(charge, 0.0, 2e-15 + 2e-12 * unsigned_charge);
    close(
        aqueous.charge_balance_residual.value(),
        charge,
        8.0 * f64::EPSILON * unsigned_charge,
    );
    let owned_iron = aqueous
        .spectators
        .iter()
        .find(|s| s.species.as_str() == SPECIES[0])
        .unwrap();
    assert_eq!(owned_iron.charge, 3);
    assert_eq!(owned_iron.amount, result.dissolved_amount);
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 298.15;
    relative(
        aqueous.proton_activity.value() * aqueous.hydroxide_activity.value(),
        (-79_800.0 / rt).exp(),
    );
    if total > 0.0 {
        let log_saturation = (result.dissolved_amount.value() / 10.0).ln()
            - 3.0 * aqueous.proton_activity.value().ln()
            + 900.0 / rt;
        close(
            result.log_saturation_ratio.unwrap().value(),
            log_saturation,
            3e-14,
        );
        if result.status == PrecipitationStatus::SolidPresent {
            assert!(log_saturation.abs() <= 1e-9);
        } else {
            assert!(log_saturation <= 0.0);
        }
    }
}

#[test]
fn g1_actual_goethite_sources_close_iron_chloride_charge_and_phase_saturation_together() {
    let model = model();
    let answer = solve(&model, 0.001, 0.003, 0.0);
    assert_eq!(answer.status, PrecipitationStatus::SolidPresent);
    check_closure(&answer, 0.001, 0.003, 0.0);
    // Independent 70-digit Decimal solution of the scalar saturated balance:
    // h + 3*h^3/Kprecip - Kw/h = (n_Cl-n_Na)/n_ref.
    relative(
        answer.aqueous.proton_activity.value(),
        0.000_299_999_978_536_417_3,
    );
    relative(answer.dissolved_amount.value(), 1.877_977_966_836_884e-10);
    close(answer.aqueous.ph.value(), 3.522_878_776_352_057, 1e-9);
    assert!(answer.aqueous_solves > 1 && answer.aqueous_iterations > 0);
    assert!(answer.aqueous.families.is_empty());
    assert_eq!(answer.model().standard().molar().value(), -900.0);
    assert_eq!(answer.model().standard().reference_axis(), REFERENCE);
    for standard in [
        answer.model().standard(),
        answer.model().aqueous_model().water().standard(),
    ] {
        for term in standard.terms() {
            let pack = materials()
                .packs
                .iter()
                .find(|p| p.pack_id() == term.pack_id)
                .unwrap();
            pack.claims().verify_receipt(&term.answer.receipt).unwrap();
            assert!(term.answer.receipt.observation_backed);
        }
    }
}

#[test]
fn g1_g3_acid_suppresses_the_solid_and_declared_base_retains_trace_dissolved_iron() {
    let model = model();
    let acidic = solve(&model, 0.001, 1.003, 0.0);
    assert_eq!(acidic.status, PrecipitationStatus::SolidAbsent);
    assert_eq!(acidic.dissolved_amount.value(), 0.001);
    assert_eq!(acidic.solid_amount.value(), 0.0);
    assert!(acidic.log_saturation_ratio.unwrap().value() < 0.0);
    check_closure(&acidic, 0.001, 1.003, 0.0);

    let neutralized = solve(&model, 0.001, 0.003, 0.003);
    let excess_base = solve(&model, 0.001, 0.003, 0.004);
    for (state, sodium, expected) in [
        (&neutralized, 0.003, 7.443_789_931_007_482e-21),
        (&excess_base, 0.004, 7.966_361_905_894_26e-30),
    ] {
        assert_eq!(state.status, PrecipitationStatus::SolidPresent);
        assert!(state.dissolved_amount.value() > 0.0);
        relative(state.dissolved_amount.value(), expected);
        check_closure(state, 0.001, 0.003, sodium);
    }
    assert!(neutralized.aqueous.ph.value() > acidic.aqueous.ph.value());
    assert!(excess_base.aqueous.ph.value() > neutralized.aqueous.ph.value());
    assert!(excess_base.dissolved_amount.value() < neutralized.dissolved_amount.value());
    // Trace Fe is retained as its own log-solved state even when subtraction
    // from the much larger solid inventory cannot resolve it.
    assert_eq!(0.001 - excess_base.solid_amount.value(), 0.0);
    assert!(excess_base.dissolved_amount.value() < f64::EPSILON * 0.001);
}

#[test]
fn g3_empty_component_and_a_caller_floor_do_not_fabricate_precipitation() {
    let model = model();
    let empty = solve(&model, 0.0, 0.0, 0.0);
    assert_eq!(empty.status, PrecipitationStatus::EmptyComponent);
    assert_eq!(empty.dissolved_amount.value(), 0.0);
    assert_eq!(empty.solid_amount.value(), 0.0);
    assert_eq!(empty.log_saturation_ratio, None);
    check_closure(&empty, 0.0, 0.0, 0.0);
    let ions = fixed_ions(0.003, 0.004);
    let original = ions.clone();
    assert!(
        model
            .speciate(
                Amount::new(0.001),
                &[],
                &ions,
                Amount::new(10.0),
                PrecipitationOptions {
                    minimum_dissolved_amount: Amount::new(1e-20),
                    ..PrecipitationOptions::default()
                },
                || false
            )
            .is_err()
    );
    assert_eq!(ions, original);
}
