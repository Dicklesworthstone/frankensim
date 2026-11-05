//! G0/G1/G3: finite stocks, independent exponential/event oracles and splitting.

use super::*;
use crate::SpeciesId;
use fs_qty::AmountConcentration;

fn model(capacity: f64, conductance: f64, policy: PrecipitationPolicyV1) -> SolubilityKineticsV1 {
    SolubilityKineticsV1::new(
        IdealSolubilityV1::new(
            SpeciesId::new("solvent").unwrap(),
            SpeciesId::new("solute").unwrap(),
            Temperature::new(300.0),
            AmountConcentration::new(capacity),
            0.0,
            [Temperature::new(280.0), Temperature::new(350.0)],
        )
        .unwrap(),
        VolumetricFlowRate::new(conductance),
        policy,
    )
    .unwrap()
}

fn state(total: f64, dissolved: f64) -> SolubilityKineticStateV1 {
    SolubilityKineticStateV1::new(Amount::new(total), Amount::new(dissolved)).unwrap()
}

fn step(
    model: &SolubilityKineticsV1,
    state: &SolubilityKineticStateV1,
    dt: f64,
) -> SolubilityKineticStepV1 {
    model
        .advance(
            Temperature::new(300.0),
            Volume::new(1.0),
            state,
            Time::new(dt),
        )
        .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-14 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn exact_exponential_dissolution_and_precipitation_are_opposite_transfers() {
    let model = model(2.0, 0.5, PrecipitationPolicyV1::ExistingSolidOnly);
    for (initial, expected, direction) in [
        (state(5.0, 0.5), 2.0 - 1.5 * (-0.75_f64).exp(), 1.0),
        (state(5.0, 4.0), 2.0 + 2.0 * (-0.75_f64).exp(), -1.0),
    ] {
        let outcome = step(&model, &initial, 1.5);
        near(outcome.state().dissolved().value(), expected);
        near(outcome.state().solid().value(), 5.0 - expected);
        assert_eq!(
            outcome.transferred_to_solution().value().signum(),
            direction
        );
        assert_eq!(outcome.state().total_solute(), initial.total_solute());
        assert!(outcome.state().balance_residual().value().abs() < 1e-14);
        assert_eq!(outcome.solid_exhaustion_time(), None);
    }
}

#[test]
fn declared_half_life_halves_active_deficits_and_composes_across_steps() {
    let active = model(2.0, 1.5, PrecipitationPolicyV1::ExistingSolidOnly);
    let volume = Volume::new(3.0);
    let half_life = active.saturation_deficit_half_life(volume).unwrap().unwrap();
    near(half_life.value(), 2.0 * core::f64::consts::LN_2);
    let advance = |state: &SolubilityKineticStateV1, duration| {
        active
            .advance(Temperature::new(300.0), volume, state, duration)
            .unwrap()
    };

    for initial in [state(15.0, 1.5), state(15.0, 12.0)] {
        let deficit = 6.0 - initial.dissolved().value();
        let direct = advance(&initial, half_life);
        near(6.0 - direct.state().dissolved().value(), 0.5 * deficit);
        assert_eq!(direct.solid_exhaustion_time(), None);
        assert_eq!(direct.state().total_solute(), initial.total_solute());

        let first = advance(&initial, Time::new(half_life.value() / 3.0));
        let split = advance(first.state(), Time::new(2.0 * half_life.value() / 3.0));
        near(
            split.state().dissolved().value(),
            direct.state().dissolved().value(),
        );
        near(split.state().solid().value(), direct.state().solid().value());
        let second_half = advance(direct.state(), half_life);
        near(
            6.0 - second_half.state().dissolved().value(),
            0.25 * deficit,
        );
    }
}

#[test]
fn half_life_preserves_solid_exhaustion_and_blocked_precipitation() {
    let exhausting = model(4.0, 1.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let half_life = exhausting
        .saturation_deficit_half_life(Volume::new(1.0))
        .unwrap()
        .unwrap();
    let end = step(&exhausting, &state(1.0, 0.0), half_life.value());
    let event = end.solid_exhaustion_time().unwrap().value();
    near(event, (4.0_f64 / 3.0).ln());
    assert!(event < half_life.value());
    assert_eq!(end.state().dissolved().value(), 1.0);
    assert_eq!(end.state().solid().value(), 0.0);
    // Saturation remains unreachable after the finite donor has run out.
    assert_eq!(4.0 - end.state().dissolved().value(), 3.0);

    let unseeded = model(2.0, 1.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let initial = state(5.0, 5.0);
    let blocked = step(&unseeded, &initial, half_life.value());
    assert!(blocked.precipitation_blocked());
    assert_eq!(blocked.state(), &initial);
}

#[test]
fn half_life_distinguishes_closed_contact_from_unrepresentable_times() {
    let closed = model(2.0, 0.0, PrecipitationPolicyV1::SeededSurface);
    assert_eq!(
        closed
            .saturation_deficit_half_life(Volume::new(1.0))
            .unwrap(),
        None,
    );
    let active = model(2.0, 1.0, PrecipitationPolicyV1::SeededSurface);
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        for owner in [&active, &closed] {
            assert!(matches!(
                owner.saturation_deficit_half_life(Volume::new(volume)),
                Err(SolubilityErrorV1::Invalid(_))
            ));
        }
    }

    // Scaling both extensive quantities preserves V/G, including subnormal
    // inputs for which forming ln(2)/G first would overflow.
    for scale in [1e-320, 1.0, 1e300] {
        let scaled = model(2.0, scale, PrecipitationPolicyV1::SeededSurface);
        near(
            scaled
                .saturation_deficit_half_life(Volume::new(scale))
                .unwrap()
                .unwrap()
                .value(),
            core::f64::consts::LN_2,
        );
    }
    for (volume, conductance) in [(1e300, 1e-300), (1e-300, 1e300)] {
        let extreme = model(2.0, conductance, PrecipitationPolicyV1::SeededSurface);
        assert!(matches!(
            extreme.saturation_deficit_half_life(Volume::new(volume)),
            Err(SolubilityErrorV1::Unrepresentable(_))
        ));
    }
}

#[test]
fn finite_solid_exhaustion_occurs_at_the_analytic_event() {
    let model = model(4.0, 0.25, PrecipitationPolicyV1::ExistingSolidOnly);
    let initial = state(2.0, 1.0);
    let event = 4.0 * (3.0_f64 / 2.0).ln();
    let before = step(&model, &initial, 0.8 * event);
    assert!(before.state().solid().value() > 0.0);
    assert_eq!(before.solid_exhaustion_time(), None);
    let after = step(&model, &initial, 2.0 * event);
    near(after.solid_exhaustion_time().unwrap().value(), event);
    assert_eq!(after.state().dissolved().value(), 2.0);
    assert_eq!(after.state().solid().value(), 0.0);
    assert_eq!(after.transferred_to_solution().value(), 1.0);
    assert_eq!(step(&model, after.state(), event).state(), after.state());
}

#[test]
fn a_small_positive_capacity_excess_keeps_its_finite_exhaustion_event() {
    let capacity = f64::from_bits(1.0_f64.to_bits() + 1);
    let active = model(capacity, 1.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let initial = state(1.0, 0.25);
    let event = (1.0 + 0.75 / (capacity - 1.0)).ln();
    let before = step(&active, &initial, event - 0.25);
    assert!(before.state().solid().value() > 0.0);
    assert_eq!(before.solid_exhaustion_time(), None);
    let after = step(&active, before.state(), 0.5);
    near(after.solid_exhaustion_time().unwrap().value(), 0.25);
    assert_eq!(after.state().solid().value(), 0.0);
    assert_eq!(after.state().dissolved().value(), 1.0);
}

#[test]
fn overflowing_exposure_uses_a_finite_phase_or_event_limit() {
    let saturated = model(2.0, 1e300, PrecipitationPolicyV1::ExistingSolidOnly);
    let end = step(&saturated, &state(5.0, 0.0), 1e300);
    assert_eq!(end.state().dissolved().value(), 2.0);
    assert_eq!(end.state().solid().value(), 3.0);
    assert_eq!(end.solid_exhaustion_time(), None);
    let exhausting = model(4.0, 1e300, PrecipitationPolicyV1::ExistingSolidOnly);
    let end = step(&exhausting, &state(2.0, 1.0), 1e300);
    near(
        end.solid_exhaustion_time().unwrap().value(),
        (1.5_f64).ln() / 1e300,
    );
    assert_eq!(end.state().dissolved().value(), 2.0);
    assert_eq!(end.state().solid().value(), 0.0);
}

#[test]
fn absence_of_solid_does_not_invent_nuclei() {
    let initial = state(5.0, 5.0);
    let unseeded = model(2.0, 1.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let blocked = step(&unseeded, &initial, 10.0);
    assert_eq!(blocked.state(), &initial);
    assert!(blocked.precipitation_blocked());
    let seeded = model(2.0, 1.0, PrecipitationPolicyV1::SeededSurface);
    let growth = step(&seeded, &initial, 1.0);
    near(
        growth.state().dissolved().value(),
        2.0 + 3.0 * (-1.0_f64).exp(),
    );
    assert!(growth.state().solid().value() > 0.0);
    assert!(!growth.precipitation_blocked());
}

#[test]
fn short_times_and_late_positive_donors_are_retained_across_steps() {
    let active = model(1.0, 1.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let initial = state(1.0, 0.0);
    let early = step(&active, &initial, 1e-18);
    near(early.state().dissolved().value(), 1e-18);
    let late = step(&active, &initial, 100.0);
    assert_eq!(late.state().dissolved().value(), 1.0);
    near(late.state().solid().value(), (-100.0_f64).exp());
    let later = step(&active, late.state(), 5.0);
    near(later.state().solid().value(), (-105.0_f64).exp());
    assert!(later.state().solid().value() < late.state().solid().value());

    // Subnormal/underflowed exponential factors must not lose the accuracy or
    // existence of a representable scaled donor.
    let huge = model(1e100, 1.0, PrecipitationPolicyV1::ExistingSolidOnly);
    for time in [720.0, 745.0, 800.0] {
        let tail = step(&huge, &state(1e100, 0.0), time);
        near(
            tail.state().solid().value(),
            (100.0 * 10.0_f64.ln() - time).exp(),
        );
    }
}

#[test]
fn semigroup_and_extensive_rescaling_hold_across_exhaustion() {
    for (capacity, initial) in [
        (2.0, state(8.0, 0.0)),
        (4.0, state(2.0, 1.0)),
        (1.0, state(6.0, 5.0)),
    ] {
        let base = model(capacity, 0.7, PrecipitationPolicyV1::ExistingSolidOnly);
        let direct = step(&base, &initial, 3.0);
        let first = step(&base, &initial, 1.1);
        let second = step(&base, first.state(), 1.9);
        near(
            second.state().dissolved().value(),
            direct.state().dissolved().value(),
        );
        near(
            second.state().solid().value(),
            direct.state().solid().value(),
        );
        for scale in [1e-12, 1e12] {
            let scaled = model(
                capacity,
                0.7 * scale,
                PrecipitationPolicyV1::ExistingSolidOnly,
            )
            .advance(
                Temperature::new(300.0),
                Volume::new(scale),
                &state(
                    initial.total_solute().value() * scale,
                    initial.dissolved().value() * scale,
                ),
                Time::new(3.0),
            )
            .unwrap();
            near(
                scaled.state().dissolved().value() / scale,
                direct.state().dissolved().value(),
            );
            near(
                scaled.state().solid().value() / scale,
                direct.state().solid().value(),
            );
        }
    }
}

#[test]
fn zero_time_closed_contact_and_equilibrium_are_exact() {
    let active = model(2.0, 1.0, PrecipitationPolicyV1::SeededSurface);
    for initial in [state(0.0, 0.0), state(5.0, 2.0), state(1.0, 1.0)] {
        assert_eq!(step(&active, &initial, 20.0).state(), &initial);
    }
    let unequilibrated = state(5.0, 1.0);
    assert_eq!(step(&active, &unequilibrated, 0.0).state(), &unequilibrated);
    let closed = model(2.0, 0.0, PrecipitationPolicyV1::SeededSurface);
    assert_eq!(
        step(&closed, &unequilibrated, 20.0).state(),
        &unequilibrated
    );
}

#[test]
fn explicit_phase_stocks_preserve_trace_donors_and_refuse_invalid_totals() {
    let trace = 1e-100;
    for (dissolved, solid) in [(1.0, trace), (trace, 1.0)] {
        let initial =
            SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))
                .unwrap();
        assert_eq!(initial.dissolved().value(), dissolved);
        assert_eq!(initial.solid().value(), solid);
        assert_eq!(initial.total_solute().value(), 1.0);
        assert_eq!(initial.balance_residual().value(), trace);
    }

    // Replacing the dissolved inventory after transport must not reconstruct
    // a rounded-away solid stock from total minus dissolved. The retained
    // trace remains an active donor on the next physical phase step.
    let initial =
        SolubilityKineticStateV1::from_phases(Amount::new(1.0), Amount::new(trace)).unwrap();
    let active = model(1.0, 1.0, PrecipitationPolicyV1::ExistingSolidOnly);
    let evolved = step(&active, &initial, 5.0);
    near(evolved.state().solid().value(), trace * (-5.0_f64).exp());
    assert_eq!(evolved.state().total_solute(), initial.total_solute());
    assert_eq!(evolved.solid_exhaustion_time(), None);

    for (dissolved, solid) in [
        (-1.0, 2.0),
        (1.0, -1.0),
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
        (f64::MAX, f64::MAX),
    ] {
        assert!(
            SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))
                .is_err()
        );
    }
}

#[test]
fn invalid_states_and_conditions_refuse() {
    for (total, dissolved) in [
        (-1.0, 0.0),
        (1.0, -1.0),
        (1.0, 2.0),
        (f64::INFINITY, 0.0),
        (1.0, f64::NAN),
    ] {
        assert!(SolubilityKineticStateV1::new(Amount::new(total), Amount::new(dissolved)).is_err());
    }
    let model = model(2.0, 1.0, PrecipitationPolicyV1::SeededSurface);
    let initial = state(5.0, 0.0);
    for (temperature, volume, duration) in [
        (279.0, 1.0, 1.0),
        (300.0, 0.0, 1.0),
        (300.0, 1.0, -1.0),
        (300.0, 1.0, f64::NAN),
    ] {
        assert!(
            model
                .advance(
                    Temperature::new(temperature),
                    Volume::new(volume),
                    &initial,
                    Time::new(duration)
                )
                .is_err()
        );
        assert_eq!(initial, state(5.0, 0.0));
    }
    assert!(
        SolubilityKineticsV1::new(
            model.law().clone(),
            VolumetricFlowRate::new(-1.0),
            PrecipitationPolicyV1::SeededSurface
        )
        .is_err()
    );
}
