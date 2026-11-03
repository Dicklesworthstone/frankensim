//! Actual JANAF and USGS packs drive standard-state mineral thermodynamics.
//! Known answers use independent arithmetic from the published source values;
//! they are not measured solution equilibria or propagated uncertainty claims.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, QueryPoint, UncertaintyModel};
use fs_qty::{
    Dimensionless, Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        FormationEnthalpySource, ReactionBalanceError, minerals::MineralReaction,
        thermodynamics::caloric::EntropyConvention,
    },
};
use std::{
    path::Path,
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
};

const JANAF: &str = "source-elemental-reference-janaf";
const USGS: &str = "source-elemental-reference-usgs-b2131";
type Participant = (&'static str, &'static str, &'static str);
const SULFUR: [Participant; 3] = [
    (
        "SO2(g)",
        "sulfur-dioxide-gas-janaf-mining",
        "source-phase-gas",
    ),
    ("O2(g)", "oxygen-gas-janaf-mining", "source-phase-gas"),
    (
        "SO3(g)",
        "sulfur-trioxide-gas-janaf-mining",
        "source-phase-gas",
    ),
];
const CALCITE: [Participant; 4] = [
    (
        "CaCO3(s)",
        "calcite-solid-usgs-b2131-thermochemistry",
        "source-phase-solid",
    ),
    (
        "H+(aq)",
        "hydrogen-ion-aqueous-usgs-b2131",
        "source-phase-aqueous",
    ),
    (
        "Ca[2+](aq)",
        "calcium-ion-aqueous-usgs-b2131",
        "source-phase-aqueous",
    ),
    (
        "HCO3[-](aq)",
        "bicarbonate-ion-aqueous-usgs-b2131",
        "source-phase-aqueous",
    ),
];
const PYRITE: [Participant; 6] = [
    (
        "FeS2(s)",
        "pyrite-solid-usgs-b2131-thermochemistry",
        "source-phase-solid",
    ),
    ("O2(g)", "oxygen-gas-usgs-b2131", "source-phase-gas"),
    ("H2O(l)", "water-liquid-usgs-b2131", "source-phase-liquid"),
    (
        "Fe[2+](aq)",
        "iron-ii-ion-aqueous-usgs-b2131",
        "source-phase-aqueous",
    ),
    (
        "SO4[2-](aq)",
        "sulfate-ion-aqueous-usgs-b2131",
        "source-phase-aqueous",
    ),
    (
        "H+(aq)",
        "hydrogen-ion-aqueous-usgs-b2131",
        "source-phase-aqueous",
    ),
];

fn compile(slug: &str, omit: Option<(&[&str], &str)>) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!(
        "frankensim-equilibrium-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&out).unwrap();
    let source = root.join(format!("data/matdb/seed-v1/{slug}"));
    let manifest = if let Some((properties, axis)) = omit {
        std::fs::copy(source.join("manifest.tsv"), out.join("manifest.tsv")).unwrap();
        let original = std::fs::read_to_string(source.join("properties.tsv")).unwrap();
        // Remove source authority while preserving every tabulated value.
        let changed = original
            .lines()
            .filter(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                !(fields.first() == Some(&"validity")
                    && fields.get(1).is_some_and(|p| properties.contains(p))
                    && fields.get(2) == Some(&axis))
            })
            .collect::<Vec<_>>()
            .join("\n")
            + "\n";
        assert_ne!(original, changed, "mutation must affect the actual source");
        std::fs::write(out.join("properties.tsv"), changed).unwrap();
        out.join("manifest.tsv")
    } else {
        source.join("manifest.tsv")
    };
    let output = out.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(manifest)
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap()
}

fn point(reference: &str, participant: Participant) -> QueryPoint {
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
        ("source-standard-state", 1.0),
        (reference, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (participant.2, 1.0),
    ];
    if participant.2 == "source-phase-aqueous" {
        flags.extend([
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]);
    } else {
        flags.push(("source-pure-species", 1.0));
    }
    match participant.0 {
        "CaCO3(s)" => flags.push(("source-phase-calcite", 1.0)),
        "FeS2(s)" => flags.push(("source-phase-pyrite", 1.0)),
        _ => {}
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

struct SourceSet {
    packs: Vec<NormalizedPack>,
    ids: Vec<SpeciesId>,
    points: Vec<QueryPoint>,
    reference: &'static str,
}
impl SourceSet {
    fn new(reference: &'static str, participants: &[Participant]) -> Self {
        Self {
            packs: participants.iter().map(|p| compile(p.1, None)).collect(),
            ids: participants
                .iter()
                .map(|p| SpeciesId::new(p.0).unwrap())
                .collect(),
            points: participants.iter().map(|p| point(reference, *p)).collect(),
            reference,
        }
    }
    fn offered(&self) -> Vec<FormationEnthalpySource<'_>> {
        (0..self.packs.len())
            .map(|i| FormationEnthalpySource {
                species: &self.ids[i],
                pack: &self.packs[i],
                at: &self.points[i],
                reference_axis: self.reference,
            })
            .collect()
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs().max(1.0) * 2e-12,
        "{actual} != {expected}"
    );
}

#[test]
fn janaf_sulfur_gibbs_log_k_and_activity_driving_force_use_actual_gas_sources() {
    let sources = SourceSet::new(JANAF, &SULFUR);
    let reaction = MineralReaction::SulfurDioxideOxidation.build().unwrap();
    let mut offered = sources.offered();
    let state = reaction.standard_equilibrium(&offered).unwrap();
    // Literal JANAF values: H=2*(-395765)-2*(-296842); S=2*256.769-2*248.212-205.147.
    close(state.enthalpy().molar().value(), -197_846.0);
    close(state.molar_entropy().value(), -188.033);
    close(state.molar_gibbs().value(), -141_783.96105);
    close(
        state.log_equilibrium_constant().value(),
        57.195_005_543_416_536,
    );
    close(
        state.equilibrium_constant().unwrap().value(),
        6.909_956_285_957_111e24,
    );
    assert_eq!(state.temperature().value(), 298.15);
    assert_eq!(state.pressure().value(), 100_000.0);
    for term in state.entropy_terms() {
        assert_eq!(term.convention, EntropyConvention::Absolute);
        assert_eq!(
            term.answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        sources
            .packs
            .iter()
            .find(|p| p.pack_id() == term.pack_id)
            .unwrap()
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
    offered.reverse();
    let reversed = reaction.standard_equilibrium(&offered).unwrap();
    assert_eq!(
        state.molar_gibbs().value().to_bits(),
        reversed.molar_gibbs().value().to_bits()
    );
    assert_eq!(
        state.log_equilibrium_constant().value().to_bits(),
        reversed.log_equilibrium_constant().value().to_bits()
    );
    let logs = |product_log| {
        sources
            .ids
            .iter()
            .map(|id| {
                (
                    id.clone(),
                    Dimensionless::new(if id.as_str() == "SO3(g)" {
                        product_log
                    } else {
                        0.0
                    }),
                )
            })
            .collect::<Vec<_>>()
    };
    let unit = state.at_log_activities(&logs(0.0)).unwrap();
    assert_eq!(unit.log_reaction_quotient.value(), 0.0);
    close(unit.molar_gibbs.value(), state.molar_gibbs().value());
    let equilibrium_log = state.log_equilibrium_constant().value() / 2.0;
    let at_equilibrium = state.at_log_activities(&logs(equilibrium_log)).unwrap();
    assert!(at_equilibrium.molar_gibbs.value().abs() < 1e-9);
    assert!(
        state
            .at_log_activities(&logs(equilibrium_log + 1.0))
            .unwrap()
            .molar_gibbs
            .value()
            > 0.0
    );
    assert!(
        state
            .at_log_activities(&logs(equilibrium_log - 1.0))
            .unwrap()
            .molar_gibbs
            .value()
            < 0.0
    );
}

#[test]
fn usgs_calcite_retains_negative_conventional_ionic_entropy_and_validates_activity_axes() {
    let sources = SourceSet::new(USGS, &CALCITE);
    let state = MineralReaction::CalciteBicarbonateNeutralization
        .build()
        .unwrap()
        .standard_equilibrium(&sources.offered())
        .unwrap();
    // B2131: Ca=-56.2, HCO3=98.4, calcite=91.71, H+=0 J/(mol K).
    close(state.enthalpy().molar().value(), -25_500.0);
    close(state.molar_entropy().value(), -49.51);
    close(state.molar_gibbs().value(), -10_738.5935);
    close(
        state.log_equilibrium_constant().value(),
        4.331_899_815_836_022,
    );
    close(
        state.equilibrium_constant().unwrap().value(),
        76.088_703_867_164_63,
    );
    let calcium = state
        .entropy_terms()
        .iter()
        .find(|t| t.term.species.as_str() == "Ca[2+](aq)")
        .unwrap();
    assert_eq!(calcium.convention, EntropyConvention::AqueousHydrogenZero);
    assert_eq!(calcium.answer.evidence.value.value, -56.2);
    assert_eq!(
        state
            .entropy_terms()
            .iter()
            .filter(|t| t.convention == EntropyConvention::AqueousHydrogenZero)
            .count(),
        3
    );
    let logs: Vec<_> = sources
        .ids
        .iter()
        .map(|id| (id.clone(), Dimensionless::new(0.0)))
        .collect();
    assert!(state.at_log_activities(&logs[..3]).is_err());
    let mut duplicate = logs.clone();
    duplicate[1] = duplicate[0].clone();
    assert!(state.at_log_activities(&duplicate).is_err());
    let mut foreign = logs.clone();
    foreign[1].0 = SpeciesId::new("CO2(aq)").unwrap();
    assert!(state.at_log_activities(&foreign).is_err());
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut bad = logs.clone();
        bad[1].1 = Dimensionless::new(invalid);
        assert!(state.at_log_activities(&bad).is_err());
    }
}

#[test]
fn pyrite_ferrous_endpoint_keeps_finite_log_k_when_k_cannot_be_represented() {
    let sources = SourceSet::new(USGS, &PYRITE);
    let state = MineralReaction::PyriteOxidationFerrous
        .build()
        .unwrap()
        .standard_equilibrium(&sources.offered())
        .unwrap();
    // Per 2 FeS2 +7 O2 +2 H2O ->2 Fe2+ +4 SO4(2-) +4 H+; literal B2131 rows.
    close(state.enthalpy().molar().value(), -2_904_800.0);
    close(state.molar_entropy().value(), -1822.05);
    close(state.molar_gibbs().value(), -2_361_555.7925);
    close(
        state.log_equilibrium_constant().value(),
        952.640_874_488_567_1,
    );
    assert!(matches!(
        state.equilibrium_constant(),
        Err(ReactionBalanceError::Arithmetic(_))
    ));
    let logs: Vec<_> = sources
        .ids
        .iter()
        .map(|id| {
            (
                id.clone(),
                Dimensionless::new(if id.as_str() == "H+(aq)" {
                    state.log_equilibrium_constant().value() / 4.0
                } else {
                    0.0
                }),
            )
        })
        .collect();
    assert!(
        state
            .at_log_activities(&logs)
            .unwrap()
            .molar_gibbs
            .value()
            .abs()
            < 1e-8
    );
}

#[test]
fn missing_ionic_source_reference_cannot_be_repaired_by_query_flags() {
    let mut sources = SourceSet::new(USGS, &CALCITE);
    let reaction = MineralReaction::CalciteBicarbonateNeutralization
        .build()
        .unwrap();
    for properties in [
        &["conventional-entropy"][..],
        &["formation-enthalpy", "conventional-entropy"][..],
    ] {
        sources.packs[2] = compile(CALCITE[2].1, Some((properties, "ionic")));
        let offered = sources.offered();
        assert!(reaction.standard_enthalpy_change(&offered).is_ok());
        assert!(matches!(
            reaction.standard_equilibrium(&offered),
            Err(ReactionBalanceError::Invalid(_))
        ));
    }
}

#[test]
fn enthalpy_entropy_phase_and_polymorph_must_remain_compatible() {
    let mut calcite = SourceSet::new(USGS, &CALCITE);
    calcite.packs[0] = compile(CALCITE[0].1, Some((&["standard-entropy"], "polymorph")));
    let reaction = MineralReaction::CalciteBicarbonateNeutralization
        .build()
        .unwrap();
    assert!(
        reaction
            .standard_enthalpy_change(&calcite.offered())
            .is_ok()
    );
    assert!(matches!(
        reaction.standard_equilibrium(&calcite.offered()),
        Err(ReactionBalanceError::Invalid(_))
    ));
    let mut sulfur = SourceSet::new(JANAF, &SULFUR);
    sulfur.packs[0] = compile(SULFUR[0].1, Some((&["entropy298"], "phase")));
    let reaction = MineralReaction::SulfurDioxideOxidation.build().unwrap();
    assert!(reaction.standard_enthalpy_change(&sulfur.offered()).is_ok());
    assert!(matches!(
        reaction.standard_equilibrium(&sulfur.offered()),
        Err(ReactionBalanceError::Invalid(_))
    ));
}
