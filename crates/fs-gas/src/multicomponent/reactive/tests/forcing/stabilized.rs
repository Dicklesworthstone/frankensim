//! Original RKL2 work and complete source balances under physical acceleration.

use super::*;
use crate::multicomponent::Boundary;

fn accelerated_mode(cells: usize) -> f64 {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let flow = field(&mechanism, cells, 0.7, |x| {
        let a = 0.6 + 0.2 * (core::f64::consts::TAU * x).sin();
        [a, 1.0 - a]
    });
    let mut operator = molecular::molecular(&flow, 0.3, 0.0, 0.0);
    let mut simulation = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let end = 0.04;
    let acceleration = 3.0;
    let report = simulation
        .advance_stabilized_with_body_force_to(
            &mut operator,
            Time::new(end),
            StabilizedReactiveIntervalControl::default(),
            BodyForce::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            [None; 2],
            || false,
        )
        .unwrap();
    assert_eq!(simulation.flow().time().value(), end);
    assert_eq!(report.final_time.value(), end);
    assert!(report.coupled_steps > 1);
    assert!(report.maximum_stages_used > 2);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    let mass = report.before.mass();
    near(
        report.source_change.impulse,
        mass * acceleration * end,
        3e-12,
    );
    near(
        report.source_change.work,
        mass * acceleration * end * (0.7 + 0.5 * acceleration * end),
        3e-11,
    );
    near(report.conservation_residual.momentum, 0.0, 3e-12);
    near(report.absolute_energy_residual_j, 0.0, 3e-11);
    near(report.mass_residual_kg, 0.0, 3e-13);
    // Isomer reaction and discrete periodic diffusion share this eigenmode.
    // Exact acceleration translates it by integral(v dt), not final v * t.
    let k = core::f64::consts::TAU;
    let laplacian =
        4.0 * (cells * cells) as f64 * (core::f64::consts::PI / cells as f64).sin().powi(2);
    let mean = 0.5 + 0.1 * (-4.0_f64 * end).exp();
    let amplitude = 0.2 * (-(4.0 + 0.3 * laplacian) * end).exp();
    let displacement = 0.7 * end + 0.5 * acceleration * end * end;
    let mut error = 0.0;
    for i in 0..cells {
        let x = (i as f64 + 0.5) / cells as f64;
        let q = simulation.flow().cell(i).unwrap();
        error += (q.partial_densities[0] - (mean + amplitude * (k * (x - displacement)).sin()))
            .abs()
            / cells as f64;
        let p = simulation.flow().gas().primitive(q, || false).unwrap();
        near(p.velocity, 0.7 + acceleration * end, 3e-9);
        near(p.temperature, 1.0, 3e-9);
    }
    error
}

#[test]
fn accelerating_rkl2_window_matches_independent_reaction_diffusion_translation() {
    let coarse = accelerated_mode(40);
    let fine = accelerated_mode(80);
    assert!(coarse < 1e-3, "coarse={coarse:e}");
    assert!(fine < 0.45 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn zero_source_reuses_original_stabilized_window_and_noop_still_admits_source() {
    let mechanism = mechanism(10.0, 30.0, 0.0);
    let mut original = uniform(&mechanism);
    let mut forced = uniform(&mechanism);
    let mut a = molecular::molecular(original.flow(), 0.2, 0.01, 0.1);
    let mut b = molecular::molecular(forced.flow(), 0.2, 0.01, 0.1);
    let controls = StabilizedReactiveIntervalControl::default();
    let expected = original
        .advance_stabilized_for(&mut a, Time::new(0.01), controls, || false)
        .unwrap();
    let actual = forced
        .advance_stabilized_with_body_force_for(
            &mut b,
            Time::new(0.01),
            controls,
            BodyForce::Acceleration(&[0.0; 4]),
            FORCE_CONTROL,
            [None; 2],
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(forced.accepted.states, original.accepted.states);
    assert_eq!(actual.source_change, ForceChange::default());
    let before = forced.accepted.states.clone();
    let time = forced.flow().time();
    for force in [
        BodyForce::UniformAcceleration(f64::NAN),
        BodyForce::Acceleration(&[0.0; 2]),
    ] {
        assert!(
            forced
                .advance_stabilized_with_body_force_to(
                    &mut b,
                    time,
                    controls,
                    force,
                    FORCE_CONTROL,
                    [None; 2],
                    || false,
                )
                .is_err()
        );
        assert_eq!(forced.accepted.states, before);
        assert_eq!(forced.flow().time(), time);
    }
}

#[test]
fn accelerated_bath_window_rolls_back_a_completed_prefix_on_cancel_or_work_refusal() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let make = || {
        let mut flow = field(&mechanism, 16, 0.2, |_| [0.5, 0.5]);
        flow.boundaries = Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        };
        let operator = molecular::molecular(&flow, 0.2, 0.0, 0.0);
        (
            ReactingEuler1d::new(flow, mechanism.clone(), chemistry(), || false).unwrap(),
            operator,
        )
    };
    let (mut expected, mut a) = make();
    let left = expected
        .flow()
        .gas()
        .conserved(&[0.8, 0.2], 0.0, 1.0)
        .unwrap();
    let right = expected
        .flow()
        .gas()
        .conserved(&[0.2, 0.8], 0.0, 1.0)
        .unwrap();
    let baths = [Some(left.as_cell()), Some(right.as_cell())];
    let force = BodyForce::UniformAcceleration(3.0);
    let target = Time::new(0.025);
    let controls = StabilizedReactiveIntervalControl::default();
    let mut polls = 0;
    let report = expected
        .advance_stabilized_with_body_force_to(
            &mut a,
            target,
            controls,
            force,
            FORCE_CONTROL,
            baths,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(report.coupled_steps > 1 && report.operator_evaluations > 1);
    assert!(report.source_change.impulse > 0.0);
    assert!(report.source_change.work > 0.0);
    near(report.absolute_energy_residual_j, 0.0, 5e-11);
    near(report.conservation_residual.momentum, 0.0, 5e-12);
    near(report.mass_residual_kg, 0.0, 5e-12);
    let (mut actual, mut b) = make();
    let before = actual.accepted.states.clone();
    let mut seen = 0;
    let refusal = actual
        .advance_stabilized_with_body_force_to(
            &mut b,
            target,
            controls,
            force,
            FORCE_CONTROL,
            baths,
            || {
                seen += 1;
                seen >= polls
            },
        )
        .unwrap_err();
    assert_eq!(refusal, ReactiveError::Cancelled);
    assert_eq!(actual.accepted.states, before);
    assert_eq!(actual.flow().time().value(), 0.0);
    assert!(
        actual
            .advance_stabilized_with_body_force_to(
                &mut b,
                target,
                StabilizedReactiveIntervalControl {
                    maximum_coupled_attempts: 1,
                    ..controls
                },
                force,
                FORCE_CONTROL,
                baths,
                || false,
            )
            .is_err()
    );
    assert_eq!(actual.accepted.states, before);
    assert_eq!(actual.flow().time().value(), 0.0);
    let replay = actual
        .advance_stabilized_with_body_force_to(
            &mut b,
            target,
            controls,
            force,
            FORCE_CONTROL,
            baths,
            || false,
        )
        .unwrap();
    assert_eq!(replay.final_time, report.final_time);
    near(
        replay.source_change.impulse,
        report.source_change.impulse,
        3e-12,
    );
    for (&actual, &expected) in actual.accepted.states.iter().zip(&expected.accepted.states) {
        near(actual, expected, 3e-11);
    }
}

#[test]
fn ramped_stabilized_windows_share_absolute_time_and_preserve_cold_isomer_chemistry() {
    let cells = 40;
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let flow = field(&mechanism, cells, 0.7, |x| {
        let a = 0.6 + 0.2 * (core::f64::consts::TAU * x).sin();
        [a, 1.0 - a]
    });
    let mut operator = molecular::molecular(&flow, 0.3, 0.0, 0.0);
    let mut simulation = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let initial = simulation.flow().totals();
    let (reference_time, acceleration, jerk) = (-0.025, 2.0, 120.0);
    let at_zero = acceleration - jerk * reference_time;
    let velocity_at = |time: f64| 0.7 + at_zero * time + 0.5 * jerk * time * time;
    let source = BodyForce::UniformAccelerationRamp {
        reference_time,
        acceleration,
        jerk,
    };
    let (mut impulse, mut work) = (0.0, 0.0);
    for end in [0.015, 0.04] {
        let start = simulation.flow().time().value();
        let report = simulation
            .advance_stabilized_with_body_force_to(
                &mut operator,
                Time::new(end),
                StabilizedReactiveIntervalControl::default(),
                source,
                FORCE_CONTROL,
                [None; 2],
                || false,
            )
            .unwrap();
        assert_eq!(simulation.flow().time().value(), end);
        assert!(report.coupled_steps > 1);
        assert!(report.maximum_stages_used > 2);
        assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
        let before_velocity = velocity_at(start);
        let after_velocity = velocity_at(end);
        near(
            report.source_change.impulse,
            initial.mass() * (after_velocity - before_velocity),
            3e-11,
        );
        near(
            report.source_change.work,
            0.5 * initial.mass() * (after_velocity.powi(2) - before_velocity.powi(2)),
            3e-10,
        );
        near(report.conservation_residual.momentum, 0.0, 3e-11);
        near(report.absolute_energy_residual_j, 0.0, 3e-10);
        near(report.mass_residual_kg, 0.0, 3e-12);
        impulse += report.source_change.impulse;
        work += report.source_change.work;
        // Equal-caloric isomers react without changing temperature; diffusion
        // redistributes their mode while preserving its global reaction mean.
        let expected_mean = 0.5 + 0.1 * (-4.0 * end).exp();
        near(simulation.flow().totals().species_masses[0], initial.mass() * expected_mean, 3e-7);
        for i in 0..cells {
            let q = simulation.flow().cell(i).unwrap();
            let p = simulation.flow().gas().primitive(q, || false).unwrap();
            near(p.velocity, after_velocity, 3e-9);
            near(p.temperature, 1.0, 3e-9);
        }
    }
    let end = simulation.flow().time().value();
    near(impulse, initial.mass() * (velocity_at(end) - 0.7), 3e-11);
    near(work, 0.5 * initial.mass() * (velocity_at(end).powi(2) - 0.49), 3e-10);
    let laplacian =
        4.0 * (cells * cells) as f64 * (core::f64::consts::PI / cells as f64).sin().powi(2);
    let mean = 0.5 + 0.1 * (-4.0 * end).exp();
    let amplitude = 0.2 * (-(4.0 + 0.3 * laplacian) * end).exp();
    let displacement = 0.7 * end + 0.5 * at_zero * end * end + jerk * end * end * end / 6.0;
    let mut error = 0.0;
    for i in 0..cells {
        let x = (i as f64 + 0.5) / cells as f64;
        let expected = mean + amplitude * (core::f64::consts::TAU * (x - displacement)).sin();
        error += (simulation.flow().cell(i).unwrap().partial_densities[0] - expected).abs()
            / cells as f64;
    }
    assert!(error < 1e-3, "accelerated reaction/diffusion/advection L1={error:e}");
}
