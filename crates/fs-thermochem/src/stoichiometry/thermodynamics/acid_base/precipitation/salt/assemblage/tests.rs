//! G0/G1/G3/G4 constructed equilibria. These are not measured mineral data.
use super::*;
use crate::{
    MolarEnergyQuantityV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{AcidDissociation, AcidFamily, SpeciationOptions, WaterAutoprotolysis},
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const REFERENCE: &str = "source-elemental-reference-authored-assemblage-test";
const TEMPERATURE: f64 = 298.15;
const CALCITE: &str = "CaCO3(s)";
const GYPSUM: &str = "CaSO4.2H2O(s)";
const CALCIUM: &str = "Ca[2+](aq)";

fn id(value: &str) -> SpeciesId {
    SpeciesId::new(value).unwrap()
}
fn point(phase: &str, charge: i32) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![REFERENCE, "source-standard-state", phase];
    if phase == "source-phase-aqueous" {
        flags.push("source-aqueous-unit-activity-infinite-dilution");
        if charge != 0 {
            flags.extend([
                "source-individual-ion",
                "source-ionic-reference-hydrogen-zero",
            ]);
        }
    } else {
        flags.push("source-pure-species");
    }
    for flag in flags {
        at = at
            .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    at
}
fn pack(name: &str, gibbs: f64, at: &QueryPoint) -> NormalizedPack {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in at.axes() {
        validity = validity.with_quantity(axis, at.axis_quantities()[axis], *value, *value);
    }
    let dims = MolarEnergyQuantityV1::DIMS;
    let provenance = Provenance {
        source: format!("authored salt assemblage fixture {name}"),
        license: "project fixture".into(),
        artifact: None,
    };
    let mut claim = PropertyClaim {
        key: PropertyKey::new("standard-molar-formation-gibbs-energy", dims),
        value: PropertyValue::Scalar { value: gibbs, dims },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: provenance.clone(),
    };
    let hash = claim.content_hash();
    let mut claims = ClaimSet::new();
    claim.observations.push(
        claims
            .register_observation(ObservationDataset {
                specimen: name.into(),
                method: "independent manufactured equilibrium".into(),
                artifact: hash,
                caveats: "not measured mineral solubility".into(),
                provenance,
            })
            .unwrap(),
    );
    claims.insert_claim(claim).unwrap();
    NormalizedPack::new(
        name,
        "fixture-v1",
        hash,
        "project fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}
struct Record {
    name: &'static str,
    at: QueryPoint,
    pack: NormalizedPack,
}
struct Fixture {
    records: Vec<Record>,
}
impl Fixture {
    fn new(k_calcite: f64, k_gypsum: f64) -> Self {
        let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE;
        let acid = -rt * 1e-7_f64.ln();
        let rows = [
            ("H+(aq)", 0.0, "source-phase-aqueous", 1),
            (
                "OH[-](aq)",
                -rt * 1e-14_f64.ln(),
                "source-phase-aqueous",
                -1,
            ),
            ("H2O(l)", 0.0, "source-phase-liquid", 0),
            (CALCIUM, 0.0, "source-phase-aqueous", 2),
            ("H2CO3(aq)", 0.0, "source-phase-aqueous", 0),
            ("HCO3[-](aq)", acid, "source-phase-aqueous", -1),
            ("CO3[2-](aq)", 2.0 * acid, "source-phase-aqueous", -2),
            ("HSO4[-](aq)", 0.0, "source-phase-aqueous", -1),
            ("SO4[2-](aq)", acid, "source-phase-aqueous", -2),
            (
                CALCITE,
                2.0 * acid + rt * k_calcite.ln(),
                "source-phase-solid",
                0,
            ),
            (GYPSUM, acid + rt * k_gypsum.ln(), "source-phase-solid", 0),
        ];
        Self {
            records: rows
                .into_iter()
                .map(|(name, g, phase, charge)| {
                    let at = point(phase, charge);
                    Record {
                        name,
                        pack: pack(name, g, &at),
                        at,
                    }
                })
                .collect(),
        }
    }
    fn sources<'a>(&'a self, reaction: &'a ConservedReaction) -> Vec<FormationGibbsSource<'a>> {
        reaction
            .terms()
            .iter()
            .map(|term| {
                let row = self
                    .records
                    .iter()
                    .find(|row| row.name == term.species.as_str())
                    .unwrap();
                FormationGibbsSource {
                    species: &term.species,
                    pack: &row.pack,
                    at: &row.at,
                    reference_axis: REFERENCE,
                }
            })
            .collect()
    }
    fn aqueous(&self) -> IdealAcidBaseSystem {
        let water_reaction = MineralReaction::WaterDissociation.build().unwrap();
        let water = WaterAutoprotolysis::from_sources(
            water_reaction.clone(),
            id("H2O(l)"),
            id("H+(aq)"),
            id("OH[-](aq)"),
            &self.sources(&water_reaction),
        )
        .unwrap();
        let step = |kind: MineralReaction, acid: &str, base: &str| {
            let reaction = kind.build().unwrap();
            AcidDissociation::from_sources(
                reaction.clone(),
                id(acid),
                id(base),
                id("H+(aq)"),
                &self.sources(&reaction),
            )
            .unwrap()
        };
        IdealAcidBaseSystem::new(
            water,
            vec![
                AcidFamily::new(
                    "carbon",
                    vec![
                        step(
                            MineralReaction::CarbonicComponentDissociation,
                            "H2CO3(aq)",
                            "HCO3[-](aq)",
                        ),
                        step(
                            MineralReaction::BicarbonateDissociation,
                            "HCO3[-](aq)",
                            "CO3[2-](aq)",
                        ),
                    ],
                )
                .unwrap(),
                AcidFamily::new(
                    "sulfate",
                    vec![step(
                        MineralReaction::HydrogenSulfateDissociation,
                        "HSO4[-](aq)",
                        "SO4[2-](aq)",
                    )],
                )
                .unwrap(),
            ],
        )
        .unwrap()
    }
    fn phases(&self) -> Vec<SaltPrecipitation> {
        [
            (MineralReaction::CalcitePrecipitation, CALCITE, "carbon"),
            (MineralReaction::GypsumPrecipitation, GYPSUM, "sulfate"),
        ]
        .into_iter()
        .map(|(kind, solid, family)| {
            let reaction = kind.build().unwrap();
            SaltPrecipitation::from_sources(
                self.aqueous(),
                reaction.clone(),
                id(CALCIUM),
                id(solid),
                family,
                &self.sources(&reaction),
            )
            .unwrap()
        })
        .collect()
    }
    fn model(&self) -> CommonCationSaltSystem {
        CommonCationSaltSystem::new(self.phases()).unwrap()
    }
}
fn options(scale: f64) -> SaltAssemblageOptions {
    SaltAssemblageOptions {
        phase: PrecipitationOptions {
            aqueous: SpeciationOptions {
                ph_tolerance: Dimensionless::new(1e-13),
                absolute_charge_tolerance: Amount::new(1e-20 * scale),
                relative_charge_tolerance: Dimensionless::new(1e-13),
                max_iterations: 256,
                ..SpeciationOptions::default()
            },
            log_activity_tolerance: Dimensionless::new(1e-11),
            log_saturation_tolerance: Dimensionless::new(1e-10),
            ..PrecipitationOptions::default()
        },
        absolute_balance_tolerance: Amount::new(1e-14 * scale),
        ..SaltAssemblageOptions::default()
    }
}
fn totals(carbon: f64, sulfate: f64) -> Vec<(String, Amount)> {
    vec![
        ("carbon".into(), Amount::new(carbon)),
        ("sulfate".into(), Amount::new(sulfate)),
    ]
}
fn sodium(amount: f64) -> Vec<SpectatorIon> {
    vec![SpectatorIon {
        species: id("Na[+](aq)"),
        charge: 1,
        amount: Amount::new(amount),
    }]
}
fn close(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:.16e}, expected={expected:.16e}, allowance={allowance:.3e}"
    );
}
fn phase<'a>(answer: &'a SaltAssemblageSpeciation, solid: &str) -> &'a SaltPhaseEquilibrium {
    answer
        .phases
        .iter()
        .find(|p| p.solid.as_str() == solid)
        .unwrap()
}
fn ledger(answer: &SaltAssemblageSpeciation, scale: f64) {
    let ca = answer
        .aqueous
        .spectators
        .iter()
        .find(|ion| ion.species == id(CALCIUM))
        .unwrap()
        .amount
        .value();
    let calcium_total = ca + answer.phases.iter().map(|p| p.amount.value()).sum::<f64>();
    close(calcium_total, answer.cation_total.value(), 1e-11 * scale);
    for f in &answer.families {
        let aqueous = answer
            .aqueous
            .families
            .iter()
            .find(|q| q.id == f.family)
            .unwrap();
        let dissolved: f64 = aqueous.forms.iter().map(|q| q.amount.value()).sum();
        let held: f64 = answer
            .phases
            .iter()
            .filter(|p| p.family == f.family)
            .map(|p| p.amount.value())
            .sum();
        close(dissolved + held, f.total.value(), 1e-11 * scale);
    }
}

#[test]
fn g1_independent_two_phase_state_has_one_ph_and_closes_all_components() {
    // At pH7 both carbonate Ka values equal aH: fractions=(1,1,1)/3.
    // Sulfate Ka=aH: fractions=(1,1)/2. Ca=.1, C=.3, S=.4 give IAPs .01,.02.
    // Na=.7 balances +.2 Ca charge, -.3 C and -.6 S; water ions cancel.
    let model = Fixture::new(0.01, 0.02).model();
    let a = model
        .speciate(
            Amount::new(0.6),
            &totals(0.5, 0.7),
            &sodium(0.7),
            Amount::new(1.0),
            options(1.0),
            || false,
        )
        .unwrap();
    close(a.aqueous.ph.value(), 7.0, 2e-9);
    close(a.dissolved_cation.value(), 0.1, 2e-9);
    close(phase(&a, CALCITE).amount.value(), 0.2, 2e-9);
    close(phase(&a, GYPSUM).amount.value(), 0.3, 2e-9);
    for p in &a.phases {
        assert_eq!(p.status, PrecipitationStatus::SolidPresent);
        assert!(
            p.log_saturation_ratio.unwrap().value().abs()
                + p.log_saturation_numerical_spread.unwrap().value()
                <= 1e-10
        );
    }
    ledger(&a, 1.0);
    assert!(a.sweeps > 1);
    assert!(a.aqueous_solves <= a.options.max_aqueous_solves);
}

#[test]
fn g1_one_phase_reduces_to_the_original_coordinate_owner() {
    let fixture = Fixture::new(0.01, 0.02);
    let original = fixture.phases().remove(0);
    let expected = original
        .speciate(
            Amount::new(0.6),
            &totals(0.5, 0.7),
            &sodium(0.7),
            Amount::new(1.0),
            options(1.0).phase,
            || false,
        )
        .unwrap();
    let model = CommonCationSaltSystem::new(vec![original]).unwrap();
    let a = model
        .speciate(
            Amount::new(0.6),
            &totals(0.5, 0.7),
            &sodium(0.7),
            Amount::new(1.0),
            options(1.0),
            || false,
        )
        .unwrap();
    assert_eq!(a.dissolved_cation, expected.dissolved_cation);
    assert_eq!(a.phases[0].amount, expected.solid_amount);
    assert_eq!(a.aqueous.ph, expected.aqueous.ph);
    assert_eq!(a.sweeps, 1);
}

#[test]
fn g1_absent_and_empty_phases_do_not_require_false_saturation_equalities() {
    let model = Fixture::new(0.01, 0.02).model();
    // The calcite IAP stays .01; gypsum IAP=.1*(.1/2)=.005 < .02.
    let a = model
        .speciate(
            Amount::new(0.3),
            &totals(0.5, 0.1),
            &sodium(0.25),
            Amount::new(1.0),
            options(1.0),
            || false,
        )
        .unwrap();
    close(a.aqueous.ph.value(), 7.0, 2e-9);
    close(a.dissolved_cation.value(), 0.1, 2e-9);
    close(phase(&a, CALCITE).amount.value(), 0.2, 2e-9);
    assert_eq!(phase(&a, GYPSUM).status, PrecipitationStatus::SolidAbsent);
    assert_eq!(phase(&a, GYPSUM).amount.value(), 0.0);
    close(
        phase(&a, GYPSUM).log_saturation_ratio.unwrap().value(),
        0.25_f64.ln(),
        2e-9,
    );
    let empty = model
        .speciate(
            Amount::new(0.0),
            &totals(0.0, 0.0),
            &[],
            Amount::new(1.0),
            options(1.0),
            || false,
        )
        .unwrap();
    assert_eq!(empty.sweeps, 0);
    assert!(
        empty
            .phases
            .iter()
            .all(|p| p.status == PrecipitationStatus::EmptyComponent
                && p.amount.value() == 0.0
                && p.log_saturation_ratio.is_none())
    );
}

#[test]
fn g3_scaled_stocks_and_canonical_phase_order_retain_the_same_equilibrium() {
    let fixture = Fixture::new(0.01, 0.02);
    let model = fixture.model();
    let mut reversed = fixture.phases();
    reversed.reverse();
    let reversed = CommonCationSaltSystem::new(reversed).unwrap();
    let reference = model
        .speciate(
            Amount::new(0.6),
            &totals(0.5, 0.7),
            &sodium(0.7),
            Amount::new(1.0),
            options(1.0),
            || false,
        )
        .unwrap();
    let permuted = reversed
        .speciate(
            Amount::new(0.6),
            &totals(0.5, 0.7),
            &sodium(0.7),
            Amount::new(1.0),
            options(1.0),
            || false,
        )
        .unwrap();
    assert_eq!(reference.dissolved_cation, permuted.dissolved_cation);
    assert_eq!(reference.aqueous.ph, permuted.aqueous.ph);
    for scale in [1e-8, 1e8] {
        let a = model
            .speciate(
                Amount::new(0.6 * scale),
                &totals(0.5 * scale, 0.7 * scale),
                &sodium(0.7 * scale),
                Amount::new(scale),
                options(scale),
                || false,
            )
            .unwrap();
        close(a.dissolved_cation.value() / scale, 0.1, 2e-9);
        close(a.aqueous.ph.value(), 7.0, 2e-9);
        close(phase(&a, CALCITE).amount.value() / scale, 0.2, 2e-9);
        close(phase(&a, GYPSUM).amount.value() / scale, 0.3, 2e-9);
        ledger(&a, scale);
    }
}

#[test]
fn g3_positive_dissolved_trace_survives_a_solid_that_rounds_to_the_total() {
    let model = Fixture::new(1e-40, 0.02).model();
    let a = model
        .speciate(
            Amount::new(1.0),
            &totals(1.0, 0.0),
            &[],
            Amount::new(1.0),
            options(1.0),
            || false,
        )
        .unwrap();
    assert_eq!(phase(&a, CALCITE).amount.value(), 1.0);
    assert!(a.dissolved_cation.value() > 0.0);
    close(a.dissolved_cation.value() / 1e-20, 3.0_f64.sqrt(), 1e-8);
    assert_eq!(
        phase(&a, GYPSUM).status,
        PrecipitationStatus::EmptyComponent
    );
}

#[test]
fn g0_source_disagreement_and_ambiguous_phase_family_are_refused() {
    let fixture = Fixture::new(0.01, 0.02);
    let phases = fixture.phases();
    assert!(CommonCationSaltSystem::new(vec![phases[0].clone(), phases[0].clone()]).is_err());
    let mut changed = Fixture::new(0.01, 0.02);
    let calcium = changed
        .records
        .iter_mut()
        .find(|r| r.name == CALCIUM)
        .unwrap();
    calcium.pack = pack(CALCIUM, 100.0, &calcium.at);
    assert!(
        CommonCationSaltSystem::new(vec![phases[0].clone(), changed.phases().remove(1)]).is_err()
    );
    assert!(CommonCationSaltSystem::new(vec![]).is_err());
    assert!(CommonCationSaltSystem::new(vec![phases[0].clone(); MAX_SALT_PHASES + 1]).is_err());
}

#[test]
fn g4_joint_convergence_budget_identity_and_cancellation_refuse_atomically() {
    let model = Fixture::new(0.01, 0.02).model();
    let inputs = totals(0.5, 0.7);
    let original = inputs.clone();
    for control in [
        SaltAssemblageOptions {
            max_sweeps: 1,
            ..options(1.0)
        },
        SaltAssemblageOptions {
            max_aqueous_solves: 1,
            ..options(1.0)
        },
    ] {
        assert!(matches!(
            model.speciate(
                Amount::new(0.6),
                &inputs,
                &sodium(0.7),
                Amount::new(1.0),
                control,
                || false
            ),
            Err(AcidBaseError::Budget)
        ));
    }
    assert!(matches!(
        model.speciate(
            Amount::new(0.6),
            &inputs,
            &sodium(0.7),
            Amount::new(1.0),
            options(1.0),
            || true
        ),
        Err(AcidBaseError::Cancelled)
    ));
    assert_eq!(inputs, original);
    assert!(
        model
            .speciate(
                Amount::new(0.6),
                &inputs[..1],
                &sodium(0.7),
                Amount::new(1.0),
                options(1.0),
                || false
            )
            .is_err()
    );
    let duplicate = [SpectatorIon {
        species: id(CALCIUM),
        charge: 2,
        amount: Amount::new(0.1),
    }];
    assert!(
        model
            .speciate(
                Amount::new(0.6),
                &inputs,
                &duplicate,
                Amount::new(1.0),
                options(1.0),
                || false
            )
            .is_err()
    );
}
