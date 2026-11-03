//! Real source/stoichiometry owners with explicitly authored Gibbs fixtures.
//! Ka=1e-7 and Kw=1e-14 are manufactured mathematical inputs, not HCN measurements.

use super::*;

#[path = "tests/gas_transfer.rs"]
mod gas_transfer;

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, LawId, LawParameter, NormalizedPack, ObservationDataset,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{
    Frequency,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    MolarEnergyQuantityV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    stoichiometry::thermodynamics::{
        FormationGibbsSource,
        acid_base::{AcidDissociation, AcidFamily, WaterAutoprotolysis},
    },
};

pub(super) const REFERENCE: &str = "source-elemental-reference-authored-cyanidation-test";
const TEMPERATURE: f64 = 298.15;
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
pub(super) fn point(name: &str) -> QueryPoint {
    let aqueous = name != "H2O(l)";
    let ionic = !["H2O(l)", HCN, "H2CO3(aq)", OXYGEN].contains(&name);
    let mut point = QueryPoint::new()
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
    let mut flags = vec![
        REFERENCE,
        "source-standard-state",
        if aqueous {
            "source-phase-aqueous"
        } else {
            "source-phase-liquid"
        },
    ];
    if aqueous {
        flags.push("source-aqueous-unit-activity-infinite-dilution");
    } else {
        flags.push("source-pure-species");
    }
    if ionic {
        flags.extend([
            "source-ionic-reference-hydrogen-zero",
            "source-individual-ion",
        ]);
    }
    for flag in flags {
        point = point
            .with_quantity(flag, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    point
}
pub(super) fn pack(name: &str, gibbs: f64, point: &QueryPoint) -> NormalizedPack {
    let artifact = hash_domain("cyanidation-gibbs-fixture-v1", name.as_bytes());
    let provenance = Provenance {
        source: "authored Gibbs fixture; not measured HCN/carbonate thermochemistry".into(),
        license: "project test fixture".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: name.into(),
            method: "manufactured acid/base constants".into(),
            artifact,
            caveats: "synthetic values for an independent kinetic known answer".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in point.axes() {
        validity = validity.with_quantity(axis, point.axis_quantities()[axis], *value, *value);
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                "standard-molar-formation-gibbs-energy",
                MolarEnergyQuantityV1::DIMS,
            ),
            value: PropertyValue::Scalar {
                value: gibbs,
                dims: MolarEnergyQuantityV1::DIMS,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    NormalizedPack::new(
        name,
        "authored-v1",
        artifact,
        "project test fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}
fn sources<'a>(
    reaction: &'a ConservedReaction,
    records: &'a [(SpeciesId, QueryPoint, NormalizedPack)],
) -> Vec<FormationGibbsSource<'a>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let record = records
                .iter()
                .find(|(id, _, _)| id == &term.species)
                .unwrap();
            FormationGibbsSource {
                species: &term.species,
                pack: &record.2,
                at: &record.1,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}
pub(super) fn system(buffer: bool) -> IdealAcidBaseSystem {
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE;
    let records = [
        ("H+(aq)", 0.0),
        ("H2O(l)", 0.0),
        (HYDROXIDE, -rt * 1e-14_f64.ln()),
        (HCN, 0.0),
        (CN, -rt * 1e-7_f64.ln()),
        ("H2CO3(aq)", 0.0),
        ("HCO3[-](aq)", -rt * 1e-7_f64.ln()),
        ("CO3[2-](aq)", -rt * 1e-14_f64.ln()),
    ]
    .into_iter()
    .map(|(name, g)| {
        let query = point(name);
        let source = pack(name, g, &query);
        (id(name), query, source)
    })
    .collect::<Vec<_>>();
    let water = MineralReaction::WaterDissociation.build().unwrap();
    let water = WaterAutoprotolysis::from_sources(
        water.clone(),
        id("H2O(l)"),
        id("H+(aq)"),
        id(HYDROXIDE),
        &sources(&water, &records),
    )
    .unwrap();
    let step = |kind: MineralReaction, acid: &str, base: &str| {
        let reaction = kind.build().unwrap();
        AcidDissociation::from_sources(
            reaction.clone(),
            id(acid),
            id(base),
            id("H+(aq)"),
            &sources(&reaction, &records),
        )
        .unwrap()
    };
    let mut families = vec![
        AcidFamily::new(
            "uncomplexed-cyanide",
            vec![step(MineralReaction::HydrogenCyanideDissociation, HCN, CN)],
        )
        .unwrap(),
    ];
    if buffer {
        families.push(
            AcidFamily::new(
                "carbonate-buffer",
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
        );
    }
    IdealAcidBaseSystem::new(water, families).unwrap()
}
pub(super) fn card(rate: f64) -> ConstitutiveModelCard {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (PH_AXIS, Dims::NONE, -1.0, 15.0),
        (FREE_CYANIDE_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (GOLD_AXIS, Amount::DIMS, 0.0, 2.0),
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 100.0),
        (
            TEMPERATURE_AXIS,
            Temperature::DIMS,
            TEMPERATURE,
            TEMPERATURE,
        ),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 10.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    let artifact = hash_domain(
        "cyanidation-rate-fixture-v1",
        b"independent empirical first/first/first rate, authored only",
    );
    ConstitutiveModelCard {
        law: LawId(CYANIDATION_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(
            RATE_PARAMETER.into(),
            LawParameter {
                value: rate,
                dims: CyanidationRateConstant::DIMS,
            },
        )]
        .into(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: "authored kinetic fixture, not a fitted gold ore".into(),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    }
}
pub(super) fn options() -> SpeciationOptions {
    SpeciationOptions {
        ph_tolerance: Dimensionless::new(1e-12),
        absolute_charge_tolerance: Amount::new(1e-16),
        relative_charge_tolerance: Dimensionless::new(1e-12),
        max_iterations: 256,
        ..SpeciationOptions::default()
    }
}
pub(super) fn charge(sodium: f64) -> CyanidationCharge {
    CyanidationCharge {
        accessible_gold: Amount::new(1.0),
        uncomplexed_cyanide: Amount::new(2.0),
        dissolved_oxygen: Amount::new(0.25),
        gold_complex: Amount::new(0.0),
        amount_at_unit_activity: Amount::new(1.0),
        solution_volume: Volume::new(1.0),
        temperature: Temperature::new(TEMPERATURE),
        spectators: vec![SpectatorIon {
            species: id("Na[+](aq)"),
            charge: 1,
            amount: Amount::new(sodium),
        }],
        other_families: vec![],
    }
}
fn batch(charge: CyanidationCharge) -> GoldCyanidationBatch {
    GoldCyanidationBatch::new(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", card(4.0)).unwrap(),
        charge,
        QueryPoint::new(),
        options(),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}
fn unchanged(bath: &GoldCyanidationBatch, before: &GoldCyanidationBatch) {
    assert_eq!(bath.elapsed(), before.elapsed());
    assert_eq!(bath.state().gold_dissolved, before.state().gold_dissolved);
    assert_eq!(bath.state().gold_complex, before.state().gold_complex);
    assert_eq!(
        bath.state().speciation.families,
        before.state().speciation.families
    );
    assert_eq!(bath.state().conditions, before.state().conditions);
    assert_eq!(bath.accumulated_error, before.accumulated_error);
}
fn balances(answer: &CyanidationAdvance) {
    for residual in [
        answer.conservation.gold,
        answer.conservation.cyanide,
        answer.conservation.dissolved_oxygen,
    ] {
        close(residual.value(), 0.0, 2e-14);
    }
    assert!(answer.state.accessible_gold.value() >= 0.0);
    assert!(answer.state.uncomplexed_cyanide.value() >= 0.0);
    assert!(answer.state.dissolved_oxygen.value() >= 0.0);
}

#[test]
fn independent_coupled_known_answer_uses_total_cyanide_and_actual_equilibrium_oh() {
    // At Ka=sqrt(Kw)=1e-7, n_ref=V=1 and Na=1, charge fixes pH=7.
    // HCN=CN=1-x, O2=(1-x)/4 and r=k(1-x)^3/4, independently giving
    // Au_remaining=1/sqrt(1+k*t/2). k=4,t=1.5 therefore gives x=0.5.
    let mut bath = batch(charge(1.0));
    let initial_free_cn = bath.state().free_cyanide.value();
    close(bath.state().speciation.ph.value(), 7.0, 2e-12);
    let answer = bath
        .advance(Time::new(1.5), CyanidationNumerics::default(), || false)
        .unwrap();
    close(answer.state.gold_dissolved.value(), 0.5, 2e-7);
    close(answer.state.accessible_gold.value(), 0.5, 2e-7);
    close(answer.state.free_cyanide.value(), 0.5, 2e-7);
    close(answer.state.uncomplexed_cyanide.value(), 1.0, 4e-7);
    close(answer.state.dissolved_oxygen.value(), 0.125, 5e-8);
    close(answer.state.gold_complex.value(), 0.5, 2e-7);
    close(answer.reaction_extent.value(), 0.125, 5e-8);
    close(amount(&answer.stoichiometric_amounts, HYDROXIDE), 0.5, 2e-7);
    close(
        answer.state.speciation.hydroxide_amount.value(),
        1e-7,
        1e-18,
    );
    assert_eq!(
        answer.model().reaction().id().as_str(),
        "gold-cyanidation-elsner-dissolved-oxygen"
    );
    assert_eq!(answer.model().source(), &card(4.0));
    assert!(
        !answer
            .state
            .speciation
            .model()
            .water()
            .standard()
            .terms()
            .is_empty()
    );
    balances(&answer);
    let continued = bath
        .advance(Time::new(4.5), CyanidationNumerics::default(), || false)
        .unwrap();
    close(
        continued.state.accessible_gold.value(),
        1.0 / 13.0_f64.sqrt(),
        3e-7,
    );
    assert!(amount(&continued.stoichiometric_amounts, CN) > initial_free_cn);
    assert_eq!(continued.previous_elapsed, Time::new(1.5));
    assert_eq!(continued.elapsed, Time::new(6.0));
    balances(&continued);
}

#[test]
fn actual_ph_changes_free_cyanide_rate_and_volume_is_not_the_activity_basis() {
    let mut low = batch(charge(0.8));
    let mut high = batch(charge(1.2));
    assert!(high.state().speciation.ph.value() > low.state().speciation.ph.value());
    assert!(high.state().free_cyanide.value() > low.state().free_cyanide.value());
    assert!(high.state().rate.value() > low.state().rate.value());
    low.advance(Time::new(0.1), CyanidationNumerics::default(), || false)
        .unwrap();
    high.advance(Time::new(0.1), CyanidationNumerics::default(), || false)
        .unwrap();
    assert!(high.state().gold_dissolved.value() > low.state().gold_dissolved.value());
    let ordinary = batch(charge(1.0));
    let mut larger = charge(1.0);
    larger.solution_volume = Volume::new(2.0);
    let larger = batch(larger);
    assert_eq!(
        larger.state().speciation.amount_at_unit_activity,
        Amount::new(1.0)
    );
    close(
        larger.state().speciation.ph.value(),
        ordinary.state().speciation.ph.value(),
        1e-13,
    );
    close(
        larger.state().rate.value(),
        ordinary.state().rate.value() / 4.0,
        1e-13,
    );
    assert_eq!(
        larger.state().conditions.axes()[ACTIVITY_REFERENCE_AXIS],
        1.0
    );
    assert_eq!(larger.state().conditions.axes()[VOLUME_AXIS], 2.0);
}

#[test]
fn finite_oxygen_and_cyanide_limits_remain_positive_without_clipping() {
    let mut oxygen_limited = charge(1.0);
    oxygen_limited.dissolved_oxygen = Amount::new(0.025);
    let mut oxygen_limited = batch(oxygen_limited);
    let mut cyanide_limited = charge(0.1);
    cyanide_limited.uncomplexed_cyanide = Amount::new(0.2);
    let mut cyanide_limited = batch(cyanide_limited);
    for bath in [&mut oxygen_limited, &mut cyanide_limited] {
        close(bath.maximum_gold_dissolution().value(), 0.1, 1e-15);
        let answer = bath
            .advance(Time::new(10.0), CyanidationNumerics::default(), || false)
            .unwrap();
        assert!(
            answer.state.gold_dissolved.value() > 0.099
                && answer.state.gold_dissolved.value() < 0.1
        );
        assert!(
            answer.state.dissolved_oxygen.value() > 0.0
                && answer.state.uncomplexed_cyanide.value() > 0.0
        );
        balances(&answer);
    }
    // These two manufactured limiting-stock cases have the same scalar ODE.
    close(
        oxygen_limited.state().gold_dissolved.value(),
        cyanide_limited.state().gold_dissolved.value(),
        3e-7,
    );
    let mut no_oxygen = charge(1.0);
    no_oxygen.dissolved_oxygen = Amount::new(0.0);
    let mut no_oxygen = batch(no_oxygen);
    let answer = no_oxygen
        .advance(Time::new(4.0), CyanidationNumerics::default(), || false)
        .unwrap();
    assert_eq!(answer.state.gold_dissolved.value(), 0.0);
    assert_eq!(answer.state.accessible_gold.value(), 1.0);
    assert_eq!(answer.state.uncomplexed_cyanide.value(), 2.0);
    assert_eq!(answer.elapsed, Time::new(4.0));
    balances(&answer);
}

#[test]
fn other_family_totals_and_initial_gold_complex_are_retained_in_charge_balance() {
    let mut charge = charge(1.7); // 1 mol cyanide anion + .3 carbonate charge + .4 initial complex.
    charge.gold_complex = Amount::new(0.4);
    charge
        .other_families
        .push(("carbonate-buffer".into(), Amount::new(0.3)));
    let model = GoldCyanidation::from_card(system(true), "uncomplexed-cyanide", card(4.0)).unwrap();
    let mut bath = GoldCyanidationBatch::new(model, charge, QueryPoint::new(), options()).unwrap();
    let result = bath
        .advance(Time::new(1.5), CyanidationNumerics::default(), || false)
        .unwrap();
    close(result.state.gold_dissolved.value(), 0.5, 2e-7);
    close(result.state.gold_complex.value(), 0.9, 2e-7);
    let buffer = result
        .state
        .speciation
        .families
        .iter()
        .find(|f| f.id == "carbonate-buffer")
        .unwrap();
    assert_eq!(buffer.total, Amount::new(0.3));
    close(
        buffer.forms.iter().map(|f| f.amount.value()).sum(),
        0.3,
        1e-14,
    );
    assert_eq!(
        result
            .state
            .speciation
            .spectators
            .iter()
            .find(|s| s.species.as_str() == COMPLEX)
            .unwrap()
            .amount,
        result.state.gold_complex
    );
    close(result.state.speciation.ph.value(), 7.0, 2e-12);
    balances(&result);
}

#[test]
fn refinement_and_continuation_match_the_independent_solution() {
    let mut coarse = batch(charge(1.0));
    let mut fine = coarse.clone();
    let mut continued = coarse.clone();
    let coarse_options = CyanidationNumerics {
        absolute_gold_error: Amount::new(1e-3),
        relative_gold_error: Dimensionless::new(1e-3),
        ..CyanidationNumerics::default()
    };
    coarse
        .advance(Time::new(1.5), coarse_options, || false)
        .unwrap();
    let tight = CyanidationNumerics {
        absolute_gold_error: Amount::new(1e-11),
        relative_gold_error: Dimensionless::new(1e-9),
        ..CyanidationNumerics::default()
    };
    let result = fine.advance(Time::new(1.5), tight, || false).unwrap();
    let fine_error = (fine.state().gold_dissolved.value() - 0.5).abs();
    let coarse_error = (coarse.state().gold_dissolved.value() - 0.5).abs();
    assert!(fine_error < coarse_error / 10.0);
    close(fine.state().gold_dissolved.value(), 0.5, 1e-8);
    assert!(result.rejected_steps > 0 && result.accepted_steps > 1);
    assert!(result.accumulated_error_estimate.value() <= 1.011e-9);
    continued.advance(Time::new(0.75), tight, || false).unwrap();
    let second = continued.advance(Time::new(0.75), tight, || false).unwrap();
    close(
        continued.state().gold_dissolved.value(),
        fine.state().gold_dissolved.value(),
        2e-8,
    );
    assert_eq!(second.previous_elapsed, Time::new(0.75));
    assert_eq!(second.elapsed, Time::new(1.5));
    balances(&second);
}

#[test]
fn source_temperature_family_identity_spectators_and_owned_axes_are_enforced() {
    let chemistry = system(false);
    assert!(GoldCyanidation::from_card(chemistry.clone(), "not-a-family", card(4.0)).is_err());
    let mut wrong = card(4.0);
    wrong.parameters.get_mut(RATE_PARAMETER).unwrap().dims = Frequency::DIMS;
    assert!(GoldCyanidation::from_card(chemistry.clone(), "uncomplexed-cyanide", wrong).is_err());
    for rate in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            GoldCyanidation::from_card(chemistry.clone(), "uncomplexed-cyanide", card(rate))
                .is_err()
        );
    }
    let model = GoldCyanidation::from_card(chemistry, "uncomplexed-cyanide", card(4.0)).unwrap();
    let mut bad = charge(1.0);
    bad.temperature = Temperature::new(300.0);
    assert!(GoldCyanidationBatch::new(model.clone(), bad, QueryPoint::new(), options()).is_err());
    let mut bad = charge(1.0);
    bad.spectators.push(SpectatorIon {
        species: id(COMPLEX),
        charge: -1,
        amount: Amount::new(0.1),
    });
    assert!(GoldCyanidationBatch::new(model.clone(), bad, QueryPoint::new(), options()).is_err());
    let mut bad = charge(1.0);
    bad.other_families
        .push(("uncomplexed-cyanide".into(), Amount::new(2.0)));
    assert!(GoldCyanidationBatch::new(model.clone(), bad, QueryPoint::new(), options()).is_err());
    let owned = QueryPoint::new()
        .with_quantity(PH_AXIS, QuantitySpec::dimensional(Dims::NONE), 7.0)
        .unwrap();
    assert!(GoldCyanidationBatch::new(model, charge(1.0), owned, options()).is_err());
    let mut zero = GoldCyanidationBatch::new(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", card(0.0)).unwrap(),
        charge(1.0),
        QueryPoint::new(),
        options(),
    )
    .unwrap();
    let before = zero.clone();
    zero.advance(Time::new(0.0), CyanidationNumerics::default(), || false)
        .unwrap();
    unchanged(&zero, &before);
    zero.advance(Time::new(1.0), CyanidationNumerics::default(), || false)
        .unwrap();
    assert_eq!(zero.state().gold_dissolved.value(), 0.0);
}

#[test]
fn actual_ph_support_cancellation_and_budgets_leave_no_partial_state() {
    let mut source = card(4.0);
    source.validity =
        source
            .validity
            .with_quantity(PH_AXIS, QuantitySpec::dimensional(Dims::NONE), 6.8, 7.2);
    let mut bath = GoldCyanidationBatch::new(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", source).unwrap(),
        charge(1.2),
        QueryPoint::new(),
        options(),
    )
    .unwrap();
    let before = bath.clone();
    assert!(
        matches!(bath.advance(Time::new(1.5), CyanidationNumerics::default(), || false), Err(CyanidationError::OutsideValidity { axis, .. }) if axis == PH_AXIS)
    );
    unchanged(&bath, &before);
    let mut bath = batch(charge(1.0));
    let before = bath.clone();
    let mut polls = 0;
    assert_eq!(
        bath.advance(Time::new(1.5), CyanidationNumerics::default(), || {
            polls += 1;
            polls > 1000
        })
        .unwrap_err(),
        CyanidationError::Cancelled
    );
    unchanged(&bath, &before);
    assert_eq!(
        bath.advance(
            Time::new(1.5),
            CyanidationNumerics {
                maximum_work: 1,
                ..CyanidationNumerics::default()
            },
            || false
        )
        .unwrap_err(),
        CyanidationError::Budget
    );
    unchanged(&bath, &before);
    assert_eq!(
        bath.advance(
            Time::new(1.5),
            CyanidationNumerics {
                maximum_attempts: 1,
                ..CyanidationNumerics::default()
            },
            || false
        )
        .unwrap_err(),
        CyanidationError::Budget
    );
    unchanged(&bath, &before);
    assert!(
        matches!(bath.advance(Time::new(101.0), CyanidationNumerics::default(), || false), Err(CyanidationError::OutsideValidity { axis, .. }) if axis == ELAPSED_TIME_AXIS)
    );
    unchanged(&bath, &before);
}

#[test]
fn other_acid_families_cannot_also_own_dissolved_oxygen() {
    use fs_thermochem::{
        ChargeVector, ElementId, ElementalMatrix, ReactionId, StoichiometricMatrix,
    };
    // This artificial family is valid for the acid/base owner, but O2 must not
    // simultaneously be a buffered analytical pool and the leach's DO stock.
    let acid = id("HO2[+](aq)");
    let species = vec![acid.clone(), id("H+(aq)"), id(OXYGEN)];
    let reaction_id = ReactionId::new("authored-protonated-oxygen").unwrap();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("H").unwrap(), ElementId::new("O").unwrap()],
        species.clone(),
        vec![vec![1, 1, 0], vec![2, 0, 2]],
    )
    .unwrap();
    let stoichiometry = StoichiometricMatrix::new(
        species.clone(),
        vec![reaction_id.clone()],
        vec![vec![-1], vec![1], vec![1]],
    )
    .unwrap();
    let charges = ChargeVector::new(species, vec![1, 1, 0]).unwrap();
    let reaction =
        ConservedReaction::new(&elements, &stoichiometry, &charges, reaction_id).unwrap();
    let acid_gibbs = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE * 1e-7_f64.ln();
    let records = [("HO2[+](aq)", acid_gibbs), ("H+(aq)", 0.0), (OXYGEN, 0.0)]
        .into_iter()
        .map(|(name, gibbs)| {
            let query = point(name);
            let pack = pack(name, gibbs, &query);
            (id(name), query, pack)
        })
        .collect::<Vec<_>>();
    let step = AcidDissociation::from_sources(
        reaction.clone(),
        acid,
        id(OXYGEN),
        id("H+(aq)"),
        &sources(&reaction, &records),
    )
    .unwrap();
    let original = system(false);
    let mut families = original.families().to_vec();
    families.push(AcidFamily::new("oxygen-family", vec![step]).unwrap());
    let overlapping = IdealAcidBaseSystem::new(original.water().clone(), families).unwrap();
    assert!(matches!(
        GoldCyanidation::from_card(overlapping, "uncomplexed-cyanide", card(4.0)),
        Err(CyanidationError::Invalid(
            "other acid families cannot own a cyanidation reactant or product"
        ))
    ));
}
