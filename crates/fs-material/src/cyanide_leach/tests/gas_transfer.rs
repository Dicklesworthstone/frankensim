//! G0/G1/G3/G4 finite-headspace and coupled-leach behavior with authored sources.

use super::*;
use crate::cyanide_leach::gas_transfer::{
    AeratedGoldCyanidation, AerationAdvance, AerationNumerics, GAS_VOLUME_AXIS, GasTransferError,
    HENRY_SOLUBILITY_PROPERTY, HenrySolubilityCp, KLA_PARAMETER, OXYGEN_TRANSFER_LAW,
    OxygenHenrySource, OxygenTransfer, PARTIAL_PRESSURE_AXIS,
};
use fs_matdb::{ClaimSelection, SelectionPolicy};
use fs_qty::Pressure;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, StoichiometricMatrix,
    stoichiometry::thermodynamics::gas_liquid::IdealGasLiquidPartition,
};

const P0: f64 = 100_000.0;

fn henry_point(temperature: f64) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            TEMPERATURE_AXIS,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap();
    for (axis, value) in [
        ("source-solvent-pure-water", 1.0),
        ("source-henry-infinite-dilution", 1.0),
        ("source-solute-molecular-oxygen", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}
fn henry_pack(point: &QueryPoint, dims: Dims, observation_backed: bool) -> NormalizedPack {
    let artifact = hash_domain(
        "authored-Hcp-fixture-v1",
        b"Hcp=1/P0 mathematical fixture, not O2 measurements",
    );
    let provenance = Provenance {
        source: "authored Hcp analytic test; not measured oxygen solubility".into(),
        license: "project test fixture".into(),
        artifact: Some(artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "authored pure-water reference fixture".into(),
            method: "manufactured independent finite-headspace exponential".into(),
            artifact,
            caveats: "synthetic known-answer inputs, no physical calibration".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let mut validity = ValidityDomain::unconstrained();
    for (axis, value) in point.axes() {
        validity = validity.with_quantity(axis, point.axis_quantities()[axis], *value, *value);
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(HENRY_SOLUBILITY_PROPERTY, dims),
            value: PropertyValue::Scalar {
                value: 1.0 / P0,
                dims,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: if observation_backed {
                vec![observation]
            } else {
                vec![]
            },
            provenance,
        })
        .unwrap();
    NormalizedPack::new(
        "authored-Hcp-fixture",
        "authored-v1",
        artifact,
        "project test fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}
fn measured_source(temperature: f64) -> OxygenHenrySource {
    let point = henry_point(temperature);
    OxygenHenrySource::from_pack(
        henry_pack(&point, HenrySolubilityCp::DIMS, true),
        point,
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap()
}
fn henry_context(source: &OxygenHenrySource) -> QueryPoint {
    let mut point = QueryPoint::new();
    for (axis, value) in source.point().axes() {
        if axis != TEMPERATURE_AXIS {
            point = point
                .with_quantity(axis, source.point().axis_quantities()[axis], *value)
                .unwrap();
        }
    }
    point
}
fn attach_henry(
    c: CyanidationCharge,
    rate: f64,
    gas: f64,
    kla: f64,
    source: OxygenHenrySource,
    context: QueryPoint,
) -> std::result::Result<AeratedGoldCyanidation, GasTransferError> {
    let batch = GoldCyanidationBatch::new(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", card(rate)).unwrap(),
        c,
        QueryPoint::new(),
        options(),
    )
    .unwrap();
    AeratedGoldCyanidation::from_batch(
        batch,
        OxygenTransfer::from_henry_source(source, transfer_card(kla)).unwrap(),
        Amount::new(gas),
        gas_volume(),
        context,
    )
}
pub(super) fn gas_volume() -> Volume {
    // K=1 and n_ref=1 make D_equilibrium=G; this is a manufactured partition.
    Volume::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE / P0)
}

#[test]
fn measured_henry_uses_physical_volume_and_matches_finite_exchange_in_both_directions() {
    for (liquid, gas) in [(0.01, 0.25), (0.24, 0.01)] {
        let source = measured_source(TEMPERATURE);
        let context = henry_context(&source);
        let source_hash = source.pack_hash();
        let mut c = charge(1.0);
        c.dissolved_oxygen = Amount::new(liquid);
        // Deliberately unequal: Hcp uses physical volume, never activity amount.
        c.solution_volume = Volume::new(2.0);
        c.amount_at_unit_activity = Amount::new(7.0);
        let mut coupled = attach_henry(c, 0.0, gas, 0.7, source, context).unwrap();
        close(
            coupled.state().equilibrium_dissolved_oxygen.value(),
            2.0 * gas,
            2e-15,
        );
        assert!(coupled.transfer_model().equilibrium().is_none());
        let answer = coupled
            .advance(Time::new(0.4), true, tight(), || false)
            .unwrap();
        // Hcp*R*T*V_liquid/V_gas = 2; lambda=kLa*(1+2).
        let total = liquid + gas;
        let equilibrium = total * 2.0 / 3.0;
        let expected = equilibrium + (liquid - equilibrium) * (-0.7_f64 * 3.0 * 0.4).exp();
        close(
            answer.state.solution.dissolved_oxygen.value(),
            expected,
            3e-10,
        );
        close(answer.state.gas_oxygen.value(), total - expected, 3e-10);
        assert_eq!(
            answer.interval_oxygen_transfer.value().is_sign_positive(),
            expected > liquid
        );
        let retained = answer.transfer_model().henry_source().unwrap();
        assert_eq!(retained.pack_hash(), source_hash);
        let actual = answer.state.henry_answer.as_ref().unwrap();
        retained
            .pack()
            .claims()
            .verify_receipt(&actual.receipt)
            .unwrap();
        assert_eq!(actual.receipt.selected, retained.answer().receipt.selected);
        assert_eq!(
            actual.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        assert_eq!(
            actual.receipt.axis_quantities[TEMPERATURE_AXIS],
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static
            ))
        );
        conserved(&answer);
    }
}

#[test]
fn direct_henry_couples_reaction_and_retains_partitioned_clock_and_doses() {
    let source = measured_source(TEMPERATURE);
    let context = henry_context(&source);
    let mut c = charge(1.0);
    c.dissolved_oxygen = Amount::new(0.0001);
    let mut measured = attach_henry(c, 4.0, 0.25, 1.0, source, context).unwrap();
    let mut gibbs = bath(0.0001, 0.25, 4.0, 1.0);
    let mut partitioned = measured.clone();
    let answer = measured
        .advance(Time::new(0.5), true, tight(), || false)
        .unwrap();
    let reference = gibbs
        .advance(Time::new(0.5), true, tight(), || false)
        .unwrap();
    assert!(answer.state.solution.gold_dissolved.value() > 4.0 * 0.0001);
    close(
        answer.state.solution.gold_dissolved.value(),
        reference.state.solution.gold_dissolved.value(),
        3e-10,
    );
    close(
        answer.state.solution.free_cyanide.value(),
        reference.state.solution.free_cyanide.value(),
        3e-10,
    );
    partitioned
        .advance(Time::new(0.2), true, tight(), || false)
        .unwrap();
    let split = partitioned
        .advance(Time::new(0.3), true, tight(), || false)
        .unwrap();
    close(
        split.state.solution.gold_dissolved.value(),
        answer.state.solution.gold_dissolved.value(),
        3e-10,
    );
    assert_eq!(split.state.solution.elapsed, Time::new(0.5));
    let before_dose = measured.state().clone();
    measured.dose_oxygen(Amount::new(0.1), || false).unwrap();
    close(
        measured.state().gas_oxygen.value(),
        before_dose.gas_oxygen.value() + 0.1,
        1e-15,
    );
    assert_eq!(
        measured.state().solution.elapsed,
        before_dose.solution.elapsed
    );
    assert_eq!(
        measured.state().solution.dissolved_oxygen,
        before_dose.solution.dissolved_oxygen
    );
    let stopped = measured
        .advance(Time::new(0.1), false, tight(), || false)
        .unwrap();
    assert_eq!(stopped.interval_oxygen_transfer.value(), 0.0);
    conserved(&stopped);
    let before = measured.clone();
    assert!(matches!(
        measured.advance(Time::new(0.1), true, tight(), || true),
        Err(GasTransferError::Cyanidation(CyanidationError::Cancelled))
    ));
    unchanged_aeration(&measured, &before);
    conserved(&answer);
}

#[test]
fn measured_henry_admits_actual_temperature_and_all_source_flags() {
    let wrong_temperature = measured_source(TEMPERATURE + 1.0);
    let context = henry_context(&wrong_temperature);
    assert!(attach_henry(charge(1.0), 0.0, 0.25, 1.0, wrong_temperature, context).is_err());
    let source = measured_source(TEMPERATURE);
    assert!(
        attach_henry(
            charge(1.0),
            0.0,
            0.25,
            1.0,
            source.clone(),
            QueryPoint::new()
        )
        .is_err()
    );
    for (axis, value) in [
        ("source-solvent-pure-water", 0.0),
        ("source-henry-infinite-dilution", 0.0),
        ("source-solute-molecular-oxygen", 0.0),
        ("engineering-reference-model", 0.0),
        ("source-specimen-qualified", 1.0),
    ] {
        let context = henry_context(&source)
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
        assert!(
            attach_henry(charge(1.0), 0.0, 0.25, 1.0, source.clone(), context).is_err(),
            "accepted altered {axis}"
        );
    }
}

#[test]
fn measured_henry_refuses_wrong_conventions_unowned_flags_and_dynamic_source_axes() {
    let point = henry_point(TEMPERATURE);
    for (dims, observed) in [(Dims::NONE, true), (HenrySolubilityCp::DIMS, false)] {
        assert!(
            OxygenHenrySource::from_pack(
                henry_pack(&point, dims, observed),
                point.clone(),
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)
            )
            .is_err()
        );
    }
    // Query flags cannot manufacture source-owned identity or qualification.
    let mut incomplete = QueryPoint::new();
    for (axis, value) in point.axes() {
        if axis != "source-specimen-qualified" {
            incomplete = incomplete
                .with_quantity(axis, point.axis_quantities()[axis], *value)
                .unwrap();
        }
    }
    assert!(
        OxygenHenrySource::from_pack(
            henry_pack(&incomplete, HenrySolubilityCp::DIMS, true),
            point.clone(),
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)
        )
        .is_err()
    );
    for (axis, dims) in [
        (PARTIAL_PRESSURE_AXIS, Pressure::DIMS),
        (OXYGEN_AXIS, AmountConcentration::DIMS),
        ("ionic-strength", AmountConcentration::DIMS),
    ] {
        let dynamic = point
            .clone()
            .with_quantity(axis, QuantitySpec::dimensional(dims), 1.0)
            .unwrap();
        assert!(
            OxygenHenrySource::from_pack(
                henry_pack(&dynamic, HenrySolubilityCp::DIMS, true),
                dynamic,
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly)
            )
            .is_err(),
            "accepted unsupported dynamic axis {axis}"
        );
    }
}

fn partition() -> IdealGasLiquidPartition {
    partition_with_atoms(2)
}
fn partition_with_atoms(atoms: i128) -> IdealGasLiquidPartition {
    let gas = id("O2(g)");
    let aqueous = id(OXYGEN);
    let species = vec![gas.clone(), aqueous.clone()];
    let rid = ReactionId::new("authored-oxygen-phase-transfer").unwrap();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("O").unwrap()],
        species.clone(),
        vec![vec![atoms, atoms]],
    )
    .unwrap();
    let stoichiometry =
        StoichiometricMatrix::new(species.clone(), vec![rid.clone()], vec![vec![-1], vec![1]])
            .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let reaction = ConservedReaction::new(&elements, &stoichiometry, &charge, rid).unwrap();
    let aqueous_point = point(OXYGEN);
    let mut gas_point = QueryPoint::new();
    for (axis, value) in aqueous_point.axes() {
        if ![
            "source-phase-aqueous",
            "source-aqueous-unit-activity-infinite-dilution",
        ]
        .contains(&axis.as_str())
        {
            gas_point = gas_point
                .with_quantity(axis, aqueous_point.axis_quantities()[axis], *value)
                .unwrap();
        }
    }
    for axis in ["source-phase-gas", "source-pure-species"] {
        gas_point = gas_point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    let gas_pack = pack("authored-oxygen-gas", 0.0, &gas_point);
    let aqueous_pack = pack("authored-oxygen-aqueous", 0.0, &aqueous_point);
    IdealGasLiquidPartition::from_sources(
        reaction,
        gas.clone(),
        aqueous.clone(),
        &[
            FormationGibbsSource {
                species: &gas,
                pack: &gas_pack,
                at: &gas_point,
                reference_axis: REFERENCE,
            },
            FormationGibbsSource {
                species: &aqueous,
                pack: &aqueous_pack,
                at: &aqueous_point,
                reference_axis: REFERENCE,
            },
        ],
    )
    .unwrap()
}
fn transfer_card(kla: f64) -> ConstitutiveModelCard {
    let mut source = card(0.0);
    source.law = LawId(OXYGEN_TRANSFER_LAW.into());
    source.parameters = [(
        KLA_PARAMETER.into(),
        LawParameter {
            value: kla,
            dims: Frequency::DIMS,
        },
    )]
    .into();
    source.validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (PARTIAL_PRESSURE_AXIS, Pressure::DIMS, 0.0, 2_000_000.0),
        (
            TEMPERATURE_AXIS,
            Temperature::DIMS,
            TEMPERATURE,
            TEMPERATURE,
        ),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 10.0),
        (GAS_VOLUME_AXIS, Volume::DIMS, 0.001, 10.0),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 0.001, 10.0),
    ] {
        source.validity =
            source
                .validity
                .with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    source.provenance.source =
        "authored independent liquid-film coefficient, not measured aeration data".into();
    source
}
pub(super) fn transfer(kla: f64) -> OxygenTransfer {
    OxygenTransfer::from_card(partition(), transfer_card(kla)).unwrap()
}
fn configured(
    charge: CyanidationCharge,
    leach_rate: f64,
    gas: f64,
    kla: f64,
    volume: Volume,
) -> AeratedGoldCyanidation {
    let batch = GoldCyanidationBatch::new(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", card(leach_rate)).unwrap(),
        charge,
        QueryPoint::new(),
        options(),
    )
    .unwrap();
    AeratedGoldCyanidation::from_batch(
        batch,
        transfer(kla),
        Amount::new(gas),
        volume,
        QueryPoint::new(),
    )
    .unwrap()
}
fn bath(dissolved: f64, gas: f64, leach_rate: f64, kla: f64) -> AeratedGoldCyanidation {
    let mut c = charge(1.0);
    c.dissolved_oxygen = Amount::new(dissolved);
    configured(c, leach_rate, gas, kla, gas_volume())
}
fn tight() -> AerationNumerics {
    AerationNumerics {
        absolute_amount_error: Amount::new(1e-12),
        relative_amount_error: Dimensionless::new(1e-10),
        ..AerationNumerics::default()
    }
}
fn conserved(answer: &AerationAdvance) {
    for residual in [
        answer.conservation.gold,
        answer.conservation.cyanide,
        answer.conservation.liquid_oxygen,
        answer.conservation.gas_oxygen,
        answer.conservation.total_oxygen,
    ] {
        close(residual.value(), 0.0, 3e-13);
    }
    assert!(answer.state.solution.accessible_gold.value() >= 0.0);
    assert!(answer.state.solution.dissolved_oxygen.value() >= 0.0);
    assert!(answer.state.gas_oxygen.value() >= 0.0);
    close(
        answer
            .state
            .solution
            .speciation
            .charge_balance_residual
            .value(),
        0.0,
        2e-12,
    );
}
fn unchanged_aeration(current: &AeratedGoldCyanidation, before: &AeratedGoldCyanidation) {
    assert_eq!(
        format!("{:?}", current.state()),
        format!("{:?}", before.state())
    );
    assert_eq!(
        current.cumulative_oxygen_dose(),
        before.cumulative_oxygen_dose()
    );
    assert_eq!(
        current.maximum_gold_dissolution(),
        before.maximum_gold_dissolution()
    );
}

#[test]
fn independent_closed_gas_relaxation_has_both_transfer_signs_and_no_invented_reaction() {
    let kla = 0.75;
    let duration = 1.3;
    for (dissolved, gas) in [(0.1, 0.9), (0.9, 0.1)] {
        let mut tank = bath(dissolved, gas, 0.0, kla);
        let answer = tank
            .advance(Time::new(duration), true, tight(), || false)
            .unwrap();
        // Independently solve D' = k[(N-D)-D] for the manufactured K*n_ref*RT/(p0*Vg)=1.
        let expected = 0.5 + (dissolved - 0.5) * (-2.0 * kla * duration).exp();
        close(
            answer.state.solution.dissolved_oxygen.value(),
            expected,
            3e-10,
        );
        close(answer.state.gas_oxygen.value(), 1.0 - expected, 3e-10);
        close(
            answer.interval_oxygen_transfer.value(),
            expected - dissolved,
            3e-10,
        );
        assert_eq!(
            answer.interval_oxygen_transfer.value().is_sign_positive(),
            gas > dissolved
        );
        assert_eq!(answer.state.solution.gold_dissolved.value(), 0.0);
        assert_eq!(answer.state.solution.uncomplexed_cyanide.value(), 2.0);
        conserved(&answer);
    }
    let equilibrium = transfer(1.0)
        .equilibrium()
        .unwrap()
        .equilibrium_aqueous_amount(Pressure::new(P0 * 0.5), Amount::new(1.0))
        .unwrap();
    let mut stationary = bath(equilibrium.value(), 0.5, 0.0, 1.0);
    // Use the actual admitted pressure to eliminate a second authored pressure-rounding path.
    let target = stationary.state().equilibrium_dissolved_oxygen;
    stationary = bath(target.value(), 0.5, 0.0, 1.0);
    let held = stationary
        .advance(Time::new(1.0), true, tight(), || false)
        .unwrap();
    assert_eq!(held.interval_oxygen_transfer.value(), 0.0);
    let mut empty = bath(0.0, 0.0, 4.0, 1.0);
    let empty = empty
        .advance(Time::new(1.0), true, tight(), || false)
        .unwrap();
    assert_eq!(empty.interval_oxygen_transfer.value(), 0.0);
    assert_eq!(empty.interval_gold_dissolved.value(), 0.0);
    conserved(&empty);
}

#[test]
fn finite_headspace_feeds_zero_do_and_exceeds_the_original_closed_oxygen_ceiling() {
    for initial_do in [0.0, 0.025] {
        let mut c = charge(1.0);
        c.dissolved_oxygen = Amount::new(initial_do);
        let closed = batch(c.clone());
        let initial_ceiling = closed.maximum_gold_dissolution().value();
        let mut aerated = AeratedGoldCyanidation::from_batch(
            closed,
            transfer(1.0),
            Amount::new(0.5),
            gas_volume(),
            QueryPoint::new(),
        )
        .unwrap();
        close(
            aerated.maximum_gold_dissolution().unwrap().value(),
            1.0,
            0.0,
        );
        let answer = aerated
            .advance(Time::new(3.0), true, tight(), || false)
            .unwrap();
        assert!(answer.state.solution.gold_dissolved.value() > initial_ceiling + 0.1);
        assert!(answer.state.gas_oxygen.value() < 0.5);
        assert!(answer.state.solution.dissolved_oxygen.value() > 0.0);
        assert!(answer.interval_oxygen_transfer.value() > 0.0);
        close(
            answer.state.gas_oxygen.value()
                + answer.state.solution.dissolved_oxygen.value()
                + answer.state.solution.gold_dissolved.value() / 4.0,
            initial_do + 0.5,
            2e-13,
        );
        close(answer.state.solution.speciation.ph.value(), 7.0, 2e-12);
        assert_eq!(
            answer.state.solution.speciation.amount_at_unit_activity,
            c.amount_at_unit_activity
        );
        conserved(&answer);
    }
}

#[test]
fn zero_transfer_retains_original_closed_leach_and_attachment_keeps_prior_progress() {
    let mut c = charge(1.0);
    c.dissolved_oxygen = Amount::new(0.025);
    let mut closed = batch(c);
    let mut aerated = AeratedGoldCyanidation::from_batch(
        closed.clone(),
        transfer(0.0),
        Amount::new(0.5),
        gas_volume(),
        QueryPoint::new(),
    )
    .unwrap();
    let original = closed
        .advance(Time::new(3.0), CyanidationNumerics::default(), || false)
        .unwrap();
    let answer = aerated
        .advance(Time::new(3.0), true, tight(), || false)
        .unwrap();
    close(
        answer.state.solution.gold_dissolved.value(),
        original.state.gold_dissolved.value(),
        3e-8,
    );
    assert!(answer.state.solution.gold_dissolved.value() < 0.1);
    assert_eq!(answer.state.gas_oxygen.value(), 0.5);
    assert_eq!(answer.interval_oxygen_transfer.value(), 0.0);
    conserved(&answer);
    let progress = closed.state().gold_dissolved;
    let mut attached = AeratedGoldCyanidation::from_batch(
        closed,
        transfer(1.0),
        Amount::new(0.25),
        gas_volume(),
        QueryPoint::new(),
    )
    .unwrap();
    assert_eq!(attached.state().solution.elapsed, Time::new(3.0));
    assert_eq!(attached.state().solution.gold_dissolved, progress);
    assert_eq!(attached.state().cumulative_oxygen_transfer.value(), 0.0);
    let result = attached
        .advance(Time::new(1.0), true, tight(), || false)
        .unwrap();
    assert_eq!(result.previous_elapsed, Time::new(3.0));
    assert!(result.state.solution.gold_dissolved.value() > progress.value());
    conserved(&result);
}

#[test]
fn metered_doses_stops_and_restart_preserve_absolute_history_and_independent_accounts() {
    let mut tank = bath(0.0, 0.0, 4.0, 0.7);
    tank.advance(Time::new(1.0), true, tight(), || false)
        .unwrap();
    tank.dose_oxygen(Amount::new(0.03), || false).unwrap();
    assert_eq!(tank.state().solution.elapsed, Time::new(1.0));
    assert_eq!(tank.state().solution.dissolved_oxygen.value(), 0.0);
    close(
        tank.maximum_gold_dissolution().unwrap().value(),
        0.12,
        1e-15,
    );
    let paused = tank
        .advance(Time::new(1.0), false, tight(), || false)
        .unwrap();
    assert_eq!(paused.interval_gold_dissolved.value(), 0.0);
    assert_eq!(paused.interval_oxygen_transfer.value(), 0.0);
    assert_eq!(paused.state.gas_oxygen.value(), 0.03);
    let restarted = tank
        .advance(Time::new(2.0), true, tight(), || false)
        .unwrap();
    assert!(restarted.interval_gold_dissolved.value() > 0.0);
    assert!(restarted.interval_oxygen_transfer.value() > 0.0);
    let gas_before = tank.state().gas_oxygen;
    let do_before = tank.state().solution.dissolved_oxygen;
    let stopped = tank
        .advance(Time::new(1.0), false, tight(), || false)
        .unwrap();
    assert_eq!(stopped.state.gas_oxygen, gas_before);
    assert_eq!(stopped.interval_oxygen_transfer.value(), 0.0);
    assert!(stopped.state.solution.dissolved_oxygen.value() < do_before.value());
    assert!(stopped.interval_gold_dissolved.value() > 0.0);
    let before_dose = tank.state().solution.clone();
    tank.dose_oxygen(Amount::new(0.08), || false).unwrap();
    assert_eq!(
        tank.state().solution.gold_dissolved,
        before_dose.gold_dissolved
    );
    assert_eq!(
        tank.state().solution.dissolved_oxygen,
        before_dose.dissolved_oxygen
    );
    assert_eq!(tank.state().solution.elapsed, before_dose.elapsed);
    let continued = tank
        .advance(Time::new(10.0), true, tight(), || false)
        .unwrap();
    close(continued.cumulative_oxygen_dose.value(), 0.11, 1e-15);
    assert!(continued.state.solution.gold_dissolved.value() > 0.12);
    close(
        tank.maximum_gold_dissolution().unwrap().value(),
        0.44,
        1e-15,
    );
    for result in [&paused, &restarted, &stopped, &continued] {
        conserved(result);
    }
}

#[test]
fn extensive_scaling_continuation_and_refinement_respect_physics() {
    let mut original = bath(0.025, 0.3, 4.0, 0.7);
    let mut continued = original.clone();
    let mut scaled_charge = charge(2.0);
    scaled_charge.accessible_gold = Amount::new(2.0);
    scaled_charge.uncomplexed_cyanide = Amount::new(4.0);
    scaled_charge.dissolved_oxygen = Amount::new(0.05);
    scaled_charge.amount_at_unit_activity = Amount::new(2.0);
    scaled_charge.solution_volume = Volume::new(2.0);
    let mut scaled = configured(
        scaled_charge,
        4.0,
        0.6,
        0.7,
        Volume::new(2.0 * gas_volume().value()),
    );
    let a = original
        .advance(Time::new(1.3), true, tight(), || false)
        .unwrap();
    let b = scaled
        .advance(Time::new(1.3), true, tight(), || false)
        .unwrap();
    close(
        b.state.solution.gold_dissolved.value(),
        2.0 * a.state.solution.gold_dissolved.value(),
        2e-9,
    );
    close(
        b.state.solution.dissolved_oxygen.value(),
        2.0 * a.state.solution.dissolved_oxygen.value(),
        2e-9,
    );
    close(
        b.state.gas_oxygen.value(),
        2.0 * a.state.gas_oxygen.value(),
        2e-9,
    );
    close(
        b.state.solution.speciation.ph.value(),
        a.state.solution.speciation.ph.value(),
        3e-12,
    );
    continued
        .advance(Time::new(0.6), true, tight(), || false)
        .unwrap();
    let c = continued
        .advance(Time::new(0.7), true, tight(), || false)
        .unwrap();
    close(
        c.state.solution.gold_dissolved.value(),
        a.state.solution.gold_dissolved.value(),
        2e-9,
    );
    close(c.state.gas_oxygen.value(), a.state.gas_oxygen.value(), 2e-9);
    conserved(&a);
    conserved(&b);
    conserved(&c);
    let mut coarse = bath(0.1, 0.9, 0.0, 0.75);
    let mut fine = coarse.clone();
    let loose = AerationNumerics {
        absolute_amount_error: Amount::new(1e-3),
        relative_amount_error: Dimensionless::new(1e-3),
        ..tight()
    };
    let coarse = coarse
        .advance(Time::new(1.3), true, loose, || false)
        .unwrap();
    let fine = fine
        .advance(Time::new(1.3), true, tight(), || false)
        .unwrap();
    let expected = 0.5 - 0.4 * (-1.95_f64).exp();
    let coarse_error = (coarse.state.solution.dissolved_oxygen.value() - expected).abs();
    let fine_error = (fine.state.solution.dissolved_oxygen.value() - expected).abs();
    assert!(
        fine_error < coarse_error / 100.0,
        "coarse={coarse_error:e}, fine={fine_error:e}"
    );
    assert!(fine.accepted_steps > coarse.accepted_steps);
}

#[test]
fn transfer_source_identity_dimensions_locality_and_actual_coordinates_are_enforced() {
    for bad_rate in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(OxygenTransfer::from_card(partition(), transfer_card(bad_rate)).is_err());
    }
    let mut wrong = transfer_card(1.0);
    wrong.parameters.get_mut(KLA_PARAMETER).unwrap().dims = Amount::DIMS;
    assert!(OxygenTransfer::from_card(partition(), wrong).is_err());
    let mut wrong = transfer_card(1.0);
    wrong.law = LawId("residence-time-oxygen-correlation".into());
    assert!(OxygenTransfer::from_card(partition(), wrong).is_err());
    let mut aged = transfer_card(1.0);
    aged.validity = aged.validity.with_quantity(
        ELAPSED_TIME_AXIS,
        QuantitySpec::dimensional(Time::DIMS),
        0.0,
        100.0,
    );
    assert!(OxygenTransfer::from_card(partition(), aged).is_err());
    let mut untyped = transfer_card(1.0);
    untyped.validity = untyped.validity.with_quantity(
        PARTIAL_PRESSURE_AXIS,
        QuantitySpec::dimensional(Volume::DIMS),
        0.0,
        1e6,
    );
    assert!(OxygenTransfer::from_card(partition(), untyped).is_err());
    let mut unsupported = transfer_card(1.0);
    unsupported.sources.clear();
    assert!(OxygenTransfer::from_card(partition(), unsupported).is_err());
    let original = batch(charge(1.0));
    let wrong_formula =
        OxygenTransfer::from_card(partition_with_atoms(3), transfer_card(1.0)).unwrap();
    assert!(matches!(
        AeratedGoldCyanidation::from_batch(
            original.clone(),
            wrong_formula,
            Amount::new(0.5),
            gas_volume(),
            QueryPoint::new()
        ),
        Err(GasTransferError::Cyanidation(CyanidationError::Network(_)))
    ));
    let spoofed_ph = QueryPoint::new()
        .with_quantity(PH_AXIS, QuantitySpec::dimensional(Dims::NONE), 7.0)
        .unwrap();
    assert!(
        AeratedGoldCyanidation::from_batch(
            original.clone(),
            transfer(1.0),
            Amount::new(0.5),
            gas_volume(),
            spoofed_ph
        )
        .is_err()
    );
    for gas in [-1.0, f64::NAN] {
        assert!(
            AeratedGoldCyanidation::from_batch(
                original.clone(),
                transfer(1.0),
                Amount::new(gas),
                gas_volume(),
                QueryPoint::new()
            )
            .is_err()
        );
    }
    assert!(
        AeratedGoldCyanidation::from_batch(
            original.clone(),
            transfer(1.0),
            Amount::new(0.5),
            Volume::new(0.0),
            QueryPoint::new()
        )
        .is_err()
    );
    let owned = QueryPoint::new()
        .with_quantity(
            PARTIAL_PRESSURE_AXIS,
            QuantitySpec::dimensional(Pressure::DIMS),
            1.0,
        )
        .unwrap();
    assert!(
        AeratedGoldCyanidation::from_batch(
            original,
            transfer(1.0),
            Amount::new(0.5),
            gas_volume(),
            owned
        )
        .is_err()
    );
    let mut ph_limited = transfer_card(1.0);
    ph_limited.validity =
        ph_limited
            .validity
            .with_quantity(PH_AXIS, QuantitySpec::dimensional(Dims::NONE), 6.8, 7.2);
    let mut actual_ph = AeratedGoldCyanidation::from_batch(
        batch(charge(1.2)),
        OxygenTransfer::from_card(partition(), ph_limited).unwrap(),
        Amount::new(0.5),
        gas_volume(),
        QueryPoint::new(),
    )
    .unwrap();
    assert_eq!(
        actual_ph.state().transfer_conditions.axes()[PH_AXIS],
        actual_ph.state().solution.speciation.ph.value()
    );
    let before = actual_ph.clone();
    assert!(
        matches!(actual_ph.advance(Time::new(1.5), true, tight(), || false), Err(GasTransferError::Cyanidation(CyanidationError::OutsideValidity { axis, .. })) if axis == PH_AXIS)
    );
    unchanged_aeration(&actual_ph, &before);
}

#[test]
fn cancellation_work_source_departure_and_bad_doses_leave_the_accepted_state_unchanged() {
    let mut tank = bath(0.0, 0.5, 4.0, 1.0);
    let before = tank.clone();
    let mut polls = 0;
    let error = tank
        .advance(Time::new(3.0), true, tight(), || {
            polls += 1;
            polls > 500
        })
        .unwrap_err();
    assert_eq!(
        error,
        GasTransferError::Cyanidation(CyanidationError::Cancelled)
    );
    unchanged_aeration(&tank, &before);
    for budget in [
        AerationNumerics {
            maximum_work: 1,
            ..tight()
        },
        AerationNumerics {
            maximum_attempts: 1,
            maximum_step: Time::new(0.1),
            ..tight()
        },
    ] {
        assert_eq!(
            tank.advance(Time::new(3.0), true, budget, || false)
                .unwrap_err(),
            GasTransferError::Cyanidation(CyanidationError::Budget)
        );
        unchanged_aeration(&tank, &before);
    }
    for duration in [-1.0, f64::NAN] {
        assert!(
            tank.advance(Time::new(duration), true, tight(), || false)
                .is_err()
        );
        unchanged_aeration(&tank, &before);
    }
    for dose in [-1.0, f64::NAN] {
        assert!(tank.dose_oxygen(Amount::new(dose), || false).is_err());
        unchanged_aeration(&tank, &before);
    }
    assert_eq!(
        tank.dose_oxygen(Amount::new(0.1), || true).unwrap_err(),
        GasTransferError::Cyanidation(CyanidationError::Cancelled)
    );
    unchanged_aeration(&tank, &before);
    let mut narrow = transfer_card(1.0);
    narrow.validity = narrow.validity.with_quantity(
        OXYGEN_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        0.0,
        0.05,
    );
    let mut c = charge(1.0);
    c.dissolved_oxygen = Amount::new(0.0);
    let mut unsupported = AeratedGoldCyanidation::from_batch(
        batch(c),
        OxygenTransfer::from_card(partition(), narrow).unwrap(),
        Amount::new(0.5),
        gas_volume(),
        QueryPoint::new(),
    )
    .unwrap();
    let before = unsupported.clone();
    assert!(
        matches!(unsupported.advance(Time::new(1.0), true, tight(), || false), Err(GasTransferError::Cyanidation(CyanidationError::OutsideValidity { axis, .. })) if axis == OXYGEN_AXIS)
    );
    unchanged_aeration(&unsupported, &before);
    let mut pressure_limited = transfer_card(1.0);
    pressure_limited.validity = pressure_limited.validity.with_quantity(
        PARTIAL_PRESSURE_AXIS,
        QuantitySpec::dimensional(Pressure::DIMS),
        0.0,
        20_000.0,
    );
    let mut limited = AeratedGoldCyanidation::from_batch(
        batch(charge(1.0)),
        OxygenTransfer::from_card(partition(), pressure_limited).unwrap(),
        Amount::new(0.1),
        gas_volume(),
        QueryPoint::new(),
    )
    .unwrap();
    let before = limited.clone();
    assert!(
        matches!(limited.dose_oxygen(Amount::new(0.2), || false), Err(GasTransferError::Cyanidation(CyanidationError::OutsideValidity { axis, .. })) if axis == PARTIAL_PRESSURE_AXIS)
    );
    unchanged_aeration(&limited, &before);
}
