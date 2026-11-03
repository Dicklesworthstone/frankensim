//! The existing matched NBS oxygen rows drive an ideal gas/aqueous boundary.
//! Primary rows: https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv
//! Table 1:O, printed p. 2-37; source conventions are unit molality and 1 bar.
//! These checks supply no measured gas-transfer rate or salinity correction.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, QueryPoint, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, Pressure, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{FormationGibbsSource, gas_liquid::IdealGasLiquidPartition},
    },
};

const REFERENCE: &str = "source-elemental-reference-nbs1989";
const PACKS: [(&str, &str, &str); 2] = [
    (
        "O2(g)",
        "oxygen-gas-nbs1989-thermochemistry",
        "source-phase-gas",
    ),
    (
        "O2(aq)",
        "oxygen-aqueous-nbs1989-thermochemistry",
        "source-phase-aqueous",
    ),
];

struct Source {
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

fn point(phase: &str) -> QueryPoint {
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
    for (axis, value) in [
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        (phase, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (
            if phase == "source-phase-gas" {
                "source-pure-species"
            } else {
                "source-aqueous-unit-activity-infinite-dilution"
            },
            1.0,
        ),
    ] {
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
            "frankensim-gas-solubility-{}-{nonce}",
            std::process::id(),
        ));
        std::fs::create_dir(&out).unwrap();
        PACKS
            .iter()
            .map(|(species, name, phase)| {
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
                assert_eq!(pack.to_bytes(), bytes);
                Source {
                    id: SpeciesId::new(*species).unwrap(),
                    pack,
                    at: point(phase),
                }
            })
            .collect()
    })
}

fn gibbs_sources() -> Vec<FormationGibbsSource<'static>> {
    sources()
        .iter()
        .map(|s| FormationGibbsSource {
            species: &s.id,
            pack: &s.pack,
            at: &s.at,
            reference_axis: REFERENCE,
        })
        .collect()
}

fn model() -> IdealGasLiquidPartition {
    IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build().unwrap(),
        sources()[0].id.clone(),
        sources()[1].id.clone(),
        &gibbs_sources(),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1e-12),
        "{actual} != {expected}",
    );
}

#[test]
fn nbs_rows_give_an_oxygen_amount_on_the_explicit_solvent_basis() {
    let model = model();
    close(model.standard().molar().value(), 16_400.0);
    close(model.standard().temperature().value(), 298.15);
    close(model.standard().reference_pressure().value(), 100_000.0);
    // Independent 70-digit Decimal evaluation of exp(-16400/(R*298.15)),
    // R=8.31446261815324; the standard constant is dimensionless.
    close(
        model.standard().equilibrium_constant().unwrap().value(),
        0.001_339_196_503_128_502_6,
    );
    let amount = model
        .equilibrium_aqueous_amount(Pressure::new(21_000.0), Amount::new(10.0))
        .unwrap();
    close(amount.value(), 0.002_812_312_656_569_855_4);
    // n_ref=10 mol corresponds to 10 kg solvent for the NBS molal state.
    // An independently supplied liquid volume is needed to report mol/m3.
    let liquid_volume_m3 = 0.010;
    close(amount.value() / liquid_volume_m3, 0.281_231_265_656_985_54);
    assert_ne!(amount.value(), amount.value() / liquid_volume_m3);
    for term in model.standard().terms() {
        let source = sources()
            .iter()
            .find(|s| s.id == term.term.species)
            .unwrap();
        assert_eq!(term.pack_id, source.pack.pack_id());
        close(
            term.answer.evidence.value.value,
            if term.term.species == *model.gas() {
                0.0
            } else {
                16_400.0
            },
        );
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        assert!(term.answer.receipt.observation_backed);
        source
            .pack
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
}

#[test]
fn equilibrium_activity_and_pressure_scaling_use_the_same_retained_sources() {
    let model = model();
    let amount = model
        .equilibrium_aqueous_amount(Pressure::new(21_000.0), Amount::new(10.0))
        .unwrap();
    let activities = [
        (model.gas().clone(), Dimensionless::new(0.21)),
        (
            model.aqueous().clone(),
            Dimensionless::new(amount.value() / 10.0),
        ),
    ];
    assert!(
        model
            .standard()
            .at_activities(&activities)
            .unwrap()
            .molar
            .value()
            .abs()
            < 1e-9
    );
    close(
        model
            .equilibrium_aqueous_amount(Pressure::new(42_000.0), Amount::new(10.0))
            .unwrap()
            .value(),
        2.0 * amount.value(),
    );
    close(
        model
            .equilibrium_aqueous_amount(Pressure::new(21_000.0), Amount::new(20.0))
            .unwrap()
            .value(),
        2.0 * amount.value(),
    );
    assert_eq!(
        model
            .equilibrium_aqueous_amount(Pressure::new(0.0), Amount::new(10.0))
            .unwrap()
            .value(),
        0.0,
    );
    let mut reversed = gibbs_sources();
    reversed.reverse();
    let reordered = IdealGasLiquidPartition::from_sources(
        model.reaction().clone(),
        model.gas().clone(),
        model.aqueous().clone(),
        &reversed,
    )
    .unwrap();
    for (a, b) in model
        .standard()
        .terms()
        .iter()
        .zip(reordered.standard().terms())
    {
        assert_eq!(a.term, b.term);
        assert_eq!(a.answer, b.answer);
    }
    // Dissolution plus dissolved-O2 Elsner exactly equals the gas-O2 Elsner
    // endpoint. This is an atom/charge identity, not a transfer-rate model.
    let mut combined = std::collections::BTreeMap::new();
    for reaction in [
        model.reaction().clone(),
        MineralReaction::GoldCyanidationDissolvedOxygen
            .build()
            .unwrap(),
    ] {
        for term in reaction.terms() {
            *combined.entry(term.species.clone()).or_insert(0_i128) += term.coefficient;
        }
    }
    combined.retain(|_, coefficient| *coefficient != 0);
    let gas_elsner = MineralReaction::GoldCyanidation.build().unwrap();
    let expected: std::collections::BTreeMap<_, _> = gas_elsner
        .terms()
        .iter()
        .map(|term| (term.species.clone(), term.coefficient))
        .collect();
    assert_eq!(combined, expected);
}

#[test]
fn real_phase_sources_cannot_be_swapped_or_extrapolated() {
    let loaded = sources();
    let reaction = MineralReaction::OxygenDissolution.build().unwrap();
    let swapped = [
        FormationGibbsSource {
            species: &loaded[0].id,
            pack: &loaded[1].pack,
            at: &loaded[1].at,
            reference_axis: REFERENCE,
        },
        FormationGibbsSource {
            species: &loaded[1].id,
            pack: &loaded[0].pack,
            at: &loaded[0].at,
            reference_axis: REFERENCE,
        },
    ];
    // General standard-G admission cannot infer chemical identity from labels.
    // The gas/liquid owner additionally requires the original phase roles.
    assert!(reaction.standard_gibbs_change(&swapped).is_ok());
    assert!(
        IdealGasLiquidPartition::from_sources(
            reaction.clone(),
            loaded[0].id.clone(),
            loaded[1].id.clone(),
            &swapped,
        )
        .is_err()
    );
    let off_points: Vec<_> = loaded
        .iter()
        .map(|source| {
            source
                .at
                .clone()
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    300.0,
                )
                .unwrap()
        })
        .collect();
    let off_sources: Vec<_> = loaded
        .iter()
        .zip(&off_points)
        .map(|(source, at)| FormationGibbsSource {
            species: &source.id,
            pack: &source.pack,
            at,
            reference_axis: REFERENCE,
        })
        .collect();
    assert!(
        IdealGasLiquidPartition::from_sources(
            reaction,
            loaded[0].id.clone(),
            loaded[1].id.clone(),
            &off_sources,
        )
        .is_err()
    );
}
