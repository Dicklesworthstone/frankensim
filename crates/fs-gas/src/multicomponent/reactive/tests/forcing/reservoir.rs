//! Interior acceleration with original characteristic and molecular bath work.

use super::*;
use crate::multicomponent::{Boundary, Cell};

fn open_field(mechanism: &IdealGasKineticsV1) -> MixtureEuler1d {
    let mut flow = field(mechanism, 8, 0.2, |x| [0.3 + 0.2 * x, 0.7 - 0.2 * x]);
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
    };
    flow
}

#[test]
fn inert_open_reactor_matches_native_force_maps_and_original_characteristic_flux() {
    let mechanism = mechanism(0.0, 0.0, 1e4);
    let mut native = open_field(&mechanism);
    let initial = native.duplicate(&mut || false).unwrap();
    let bath = native.gas().conserved(&[0.9, 0.1], 0.0, 1.5).unwrap();
    let baths = [Some(bath.as_cell()), None];
    let force = BodyForce::UniformAcceleration(12.0);
    let h = Time::new(0.0004);
    let controls = ReactiveStepControl::default();
    // Independent public native source/transport owners retain their own
    // physical clock and receipts; the prescribed bath never receives a kick.
    let first = native
        .apply_body_force(Time::new(h.value() / 2.0), force, FORCE_CONTROL, || false)
        .unwrap();
    let transport = native
        .step_with_reservoirs(h, controls.transport, baths, || false)
        .unwrap();
    assert_eq!(transport.elapsed, h);
    let second = native
        .apply_body_force(Time::new(h.value() / 2.0), force, FORCE_CONTROL, || false)
        .unwrap();
    let mut reacting = ReactingEuler1d::new(initial, mechanism, chemistry(), || false).unwrap();
    let report = reacting
        .step_with_body_force_and_reservoirs(h, controls, force, FORCE_CONTROL, baths, || false)
        .unwrap();
    assert_eq!(report.elapsed, h);
    assert_eq!(reacting.flow().time(), native.time());
    assert!(report.transport.left_boundary_change.mass() > 0.0);
    near(
        report.transport.left_boundary_change.species_masses[0],
        0.9 * report.transport.left_boundary_change.mass(),
        3e-12,
    );
    near(
        report.transport.source_change.impulse,
        first.impulse + second.impulse,
        3e-12,
    );
    near(
        report.transport.source_change.work,
        first.work + second.work,
        3e-12,
    );
    near(report.absolute_energy_residual_j, 0.0, 3e-9);
    for (&actual, &expected) in reacting.accepted.states.iter().zip(&native.states) {
        near(actual, expected, 3e-9);
    }
    assert_eq!(bath.momentum, 0.0);
}

#[test]
fn reacting_molecular_baths_keep_interior_impulse_and_all_energy_sources_distinct() {
    let mechanism = mechanism(10.0, 30.0, 1e4);
    let initial = open_field(&mechanism);
    let left = initial.gas().conserved(&[0.9, 0.1], 0.0, 1.5).unwrap();
    let right = initial.gas().conserved(&[0.1, 0.9], 0.0, 1.0).unwrap();
    let mut operator = molecular::molecular(&initial, 0.05, 0.01, 0.1);
    let mut simulation = ReactingEuler1d::new(initial, mechanism, chemistry(), || false).unwrap();
    let acceleration = 100.0;
    let report = simulation
        .step_with_dissipation_body_force_and_reservoirs(
            &mut operator,
            Time::new(0.001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            BodyForce::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            [Some(left.as_cell()), Some(right.as_cell())],
            || false,
        )
        .unwrap();
    let halves = report.dissipation.unwrap();
    // Reaction and barycentric molecular transport conserve total mass.
    // The two force maps therefore use these independent endpoint masses,
    // including the actual mass that entered through the advective ports.
    let expected_impulse =
        0.5 * acceleration * report.elapsed.value() * (report.before.mass() + report.after.mass());
    near(
        report.transport.source_change.impulse,
        expected_impulse,
        3e-12,
    );
    assert!(report.transport.source_change.work > 0.0);
    assert!(report.chemistry.referenced_energy_change_j.abs() > 0.0);
    assert!(halves.first.boundary_change.energy.abs() > 0.0);
    for k in 0..report.before.species_count + 2 {
        near(report.conservation_residual.get(k), 0.0, 3e-10);
    }
    near(report.mass_residual_kg, 0.0, 3e-12);
    near(report.absolute_energy_residual_j, 0.0, 3e-8);
}

#[test]
fn zero_force_preserves_bath_composition_and_late_refusal_preserves_the_field() {
    let mechanism = mechanism(10.0, 30.0, 0.0);
    let make = || {
        ReactingEuler1d::new(
            open_field(&mechanism),
            mechanism.clone(),
            chemistry(),
            || false,
        )
        .unwrap()
    };
    let mut expected = make();
    let mut actual = make();
    let bath = expected
        .flow()
        .gas()
        .conserved(&[0.9, 0.1], 0.0, 1.5)
        .unwrap();
    let baths = [Some(bath.as_cell()), None];
    let h = Time::new(0.0004);
    let controls = ReactiveStepControl::default();
    let a = expected
        .step_with_reservoirs(h, controls, baths, || false)
        .unwrap();
    let b = actual
        .step_with_body_force_and_reservoirs(
            h,
            controls,
            BodyForce::UniformAcceleration(0.0),
            FORCE_CONTROL,
            baths,
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(expected.accepted.states, actual.accepted.states);

    let force = BodyForce::UniformAcceleration(12.0);
    let mut expected = make();
    let mut polls = 0;
    let reference = expected
        .step_with_body_force_and_reservoirs(h, controls, force, FORCE_CONTROL, baths, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut actual = make();
    let before = actual.accepted.states.clone();
    let mut seen = 0;
    let refusal = actual
        .step_with_body_force_and_reservoirs(h, controls, force, FORCE_CONTROL, baths, || {
            seen += 1;
            seen >= polls
        })
        .unwrap_err();
    assert_eq!(refusal, ReactiveError::Cancelled);
    assert_eq!(actual.accepted.states, before);
    assert_eq!(actual.flow().time().value(), 0.0);

    let moving = Cell {
        momentum: 0.1,
        ..bath.as_cell()
    };
    assert!(
        actual
            .step_with_body_force_and_reservoirs(
                h,
                controls,
                force,
                FORCE_CONTROL,
                [Some(moving), None],
                || false,
            )
            .is_err()
    );
    assert_eq!(actual.accepted.states, before);
    let replay = actual
        .step_with_body_force_and_reservoirs(h, controls, force, FORCE_CONTROL, baths, || false)
        .unwrap();
    assert_eq!(replay.elapsed, reference.elapsed);
    near(
        replay.transport.source_change.impulse,
        reference.transport.source_change.impulse,
        3e-12,
    );
    for (&actual, &expected) in actual.accepted.states.iter().zip(&expected.accepted.states) {
        near(actual, expected, 3e-11);
    }
}
