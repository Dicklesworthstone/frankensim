//! G0/G1 particle-area conservation, independent trajectories and event oracles.

use super::*;
use crate::SpeciesId;
use fs_qty::AmountConcentration;

fn model(capacity: f64, conductance: f64, reference: f64) -> SolubilityKineticsV1 {
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
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
    .unwrap()
    .with_particle_surface(Amount::new(reference))
    .unwrap()
}

fn phases(total: f64, dissolved: f64) -> SolubilityKineticStateV1 {
    SolubilityKineticStateV1::new(Amount::new(total), Amount::new(dissolved)).unwrap()
}

fn advance(
    model: &SolubilityKineticsV1,
    initial: &SolubilityKineticStateV1,
    duration: f64,
) -> SolubilityKineticStepV1 {
    model
        .advance(
            Temperature::new(300.0),
            Volume::new(1.0),
            initial,
            Time::new(duration),
        )
        .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

#[test]
fn particle_area_uses_the_original_reference_stock() {
    let model = model(4.0, 3.0, 2.0);
    assert_eq!(model.conductance().value(), 3.0);
    assert_eq!(model.particle_reference_solid(), Some(Amount::new(2.0)));
    assert_eq!(model.conductance_at(Amount::new(2.0)).unwrap().value(), 3.0);
    assert_eq!(model.conductance_at(Amount::new(0.0)).unwrap().value(), 0.0);
    close(
        model.conductance_at(Amount::new(16.0)).unwrap().value(),
        12.0,
        1e-14,
    );
    close(
        model.conductance_at(Amount::new(0.25)).unwrap().value(),
        0.75,
        1e-14,
    );
    assert!(model.conductance_at(Amount::new(-1.0)).is_err());
    assert!(
        model
            .saturation_deficit_half_life(Volume::new(1.0))
            .is_err()
    );
}

#[test]
fn capacity_equal_to_total_has_an_algebraic_dissolution_tail() {
    let model = model(8.0, 3.0, 8.0);
    let initial = phases(8.0, 0.0);
    for duration in [1e-5, 0.1, 1.0, 10.0, 1e6] {
        let end = advance(&model, &initial, duration);
        // Exact solution of ds/dt = -(3/8^(2/3))*s^(5/3).
        let expected_solid = 8.0 / (1.0 + 2.0 * duration).powf(1.5);
        close(end.state().solid().value(), expected_solid, 2e-13);
        assert_eq!(end.solid_exhaustion_time(), None);
        assert_eq!(end.state().total_solute(), initial.total_solute());
        let first = advance(&model, &initial, duration * 0.375);
        let split = advance(&model, first.state(), duration * 0.625);
        close(split.state().solid().value(), expected_solid, 3e-13);
    }
}

// Independent RK4 integration in the physical radius coordinate. This oracle
// never evaluates the implementation's separated integral or scalar root.
fn radius_oracle(
    initial_solid: f64,
    capacity_minus_total: f64,
    g: f64,
    reference: f64,
    t: f64,
) -> f64 {
    let coefficient = g / (3.0 * reference.cbrt().powi(2));
    let rhs = |radius: f64| -coefficient * (radius.powi(3) + capacity_minus_total);
    let count = 16_384;
    let h = t / f64::from(count);
    let mut radius = initial_solid.cbrt();
    for _ in 0..count {
        let k1 = rhs(radius);
        let k2 = rhs(radius + 0.5 * h * k1);
        let k3 = rhs(radius + 0.5 * h * k2);
        let k4 = rhs(radius + h * k3);
        radius += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    radius.powi(3)
}

#[test]
fn finite_bath_dissolution_and_growth_match_independent_radius_evolution() {
    for (capacity, total, dissolved, conductance, reference, duration) in [
        (4.0, 10.0, 1.0, 0.7, 3.0, 1.5),
        (4.0, 10.0, 7.0, 0.7, 3.0, 1.5),
        (12.0, 10.0, 1.0, 0.3, 2.0, 1.0),
        (0.25, 6.0, 2.0, 1.2, 7.0, 0.75),
    ] {
        let model = model(capacity, conductance, reference);
        let initial = phases(total, dissolved);
        let end = advance(&model, &initial, duration);
        let expected = radius_oracle(
            total - dissolved,
            capacity - total,
            conductance,
            reference,
            duration,
        );
        close(end.state().solid().value(), expected, 2e-11);
        close(end.state().dissolved().value(), total - expected, 2e-11);
        assert!(end.state().balance_residual().value().abs() <= 8.0 * f64::EPSILON * total);
        let first = advance(&model, &initial, 0.4 * duration);
        let split = advance(&model, first.state(), 0.6 * duration);
        close(
            split.state().solid().value(),
            end.state().solid().value(),
            5e-12,
        );
    }
}

#[test]
fn exhaustion_event_matches_an_independent_radius_integral() {
    let model = model(4.0, 1.0, 1.0);
    let initial = phases(1.0, 0.0);
    // t* = 3 integral_0^1 dr/(r^3+3), independently Simpson-integrated.
    let intervals = 8192;
    let h = 1.0 / f64::from(intervals);
    let mut weighted = 1.0 / 3.0 + 1.0 / 4.0;
    for i in 1..intervals {
        let radius = f64::from(i) * h;
        weighted += if i % 2 == 0 { 2.0 } else { 4.0 } / (radius.powi(3) + 3.0);
    }
    let exhaustion = weighted * h;
    let before = advance(&model, &initial, exhaustion * (1.0 - 1e-5));
    assert!(before.state().solid().value() > 0.0);
    assert_eq!(before.solid_exhaustion_time(), None);
    let after = advance(&model, &initial, exhaustion * (1.0 + 1e-5));
    assert_eq!(after.state().solid().value(), 0.0);
    assert_eq!(after.state().dissolved().value(), 1.0);
    close(
        after.solid_exhaustion_time().unwrap().value(),
        exhaustion,
        2e-12,
    );
    assert_eq!(advance(&model, after.state(), 10.0).state(), after.state());
}

#[test]
fn changing_area_slows_dissolution_and_accelerates_seeded_growth() {
    for initial in [phases(5.0, 0.5), phases(5.0, 4.0)] {
        let particle = model(2.0, 0.5, initial.solid().value());
        let fixed = SolubilityKineticsV1::new(
            particle.law().clone(),
            particle.conductance(),
            PrecipitationPolicyV1::ExistingSolidOnly,
        )
        .unwrap();
        let moving = advance(&particle, &initial, 1.5);
        let constant = advance(&fixed, &initial, 1.5);
        assert!(moving.state().dissolved().value() < constant.state().dissolved().value());
        assert_eq!(moving.state().total_solute(), initial.total_solute());
    }
}

#[test]
fn particle_surface_cannot_nucleate_or_replace_an_inert_seed() {
    let owner = model(2.0, 1.0, 1.0);
    let clear = phases(5.0, 5.0);
    let blocked = advance(&owner, &clear, 100.0);
    assert!(blocked.precipitation_blocked());
    assert_eq!(blocked.state(), &clear);
    assert_eq!(
        advance(&owner, &phases(5.0, 1.0), 0.0).state(),
        &phases(5.0, 1.0)
    );
    for reference in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            owner
                .clone()
                .with_particle_surface(Amount::new(reference))
                .is_err()
        );
    }
    let seed = SolubilityKineticsV1::new(
        owner.law().clone(),
        owner.conductance(),
        PrecipitationPolicyV1::SeededSurface,
    )
    .unwrap();
    assert!(seed.with_particle_surface(Amount::new(1.0)).is_err());
    let shut = model(2.0, 0.0, 1.0);
    assert_eq!(
        advance(&shut, &phases(5.0, 1.0), 100.0).state(),
        &phases(5.0, 1.0)
    );
    assert_eq!(
        shut.saturation_deficit_half_life(Volume::new(1.0)).unwrap(),
        None
    );
}

#[test]
fn amount_and_volume_rescaling_preserve_the_particle_trajectory() {
    let original = model(4.0, 0.7, 3.0);
    let initial = phases(10.0, 1.0);
    let end = advance(&original, &initial, 1.5);
    for scale in [1e-12, 1e-3, 1e3, 1e12] {
        let scaled = model(4.0, 0.7 * scale, 3.0 * scale);
        let result = scaled
            .advance(
                Temperature::new(300.0),
                Volume::new(scale),
                &phases(10.0 * scale, scale),
                Time::new(1.5),
            )
            .unwrap();
        close(
            result.state().solid().value() / scale,
            end.state().solid().value(),
            3e-12,
        );
        close(
            result.state().dissolved().value() / scale,
            end.state().dissolved().value(),
            3e-12,
        );
    }
}

#[test]
fn trace_solid_survives_a_dominant_dissolved_stock() {
    let model = model(1e20, 3.0, 1.0);
    let initial =
        SolubilityKineticStateV1::from_phases(Amount::new(1e20), Amount::new(1.0)).unwrap();
    let end = advance(&model, &initial, 0.5);
    let expected = 1.0 / 2.0_f64.powf(1.5);
    close(end.state().solid().value(), expected, 2e-13);
    close(end.transferred_to_solution().value(), 1.0 - expected, 2e-13);
    assert_eq!(end.state().dissolved().value(), 1e20);
    let split = advance(&model, end.state(), 0.5);
    close(
        split.state().solid().value(),
        1.0 / 3.0_f64.powf(1.5),
        3e-13,
    );
}

#[test]
fn tiny_intervals_resolve_transfer_and_long_growth_retains_dissolved_saturation() {
    let initial = phases(1.0, 0.0);
    let owner = model(2.0, 1.0, 1.0);
    let tiny = advance(&owner, &initial, 1e-20);
    close(tiny.transferred_to_solution().value(), 2e-20, 2e-13);
    close(tiny.state().dissolved().value(), 2e-20, 2e-13);
    assert_eq!(tiny.state().solid().value(), 1.0);

    let grow = model(1e-120, 1.0, 1.0);
    let saturated = advance(&grow, &phases(8.0, 5.0), 2000.0);
    close(saturated.state().dissolved().value(), 1e-120, 1e-13);
    assert_eq!(saturated.state().solid().value(), 8.0);
    assert_eq!(saturated.solid_exhaustion_time(), None);
}
