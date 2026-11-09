//! G0/G1/G3/G4: acceleration, chemical energy and molecular transport together.

use super::*;

mod reservoir;

const FORCE_CONTROL: ForceControl = ForceControl { maximum_mach: 80.0 };

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g1_accelerated_heat_release_uses_one_retried_duration_and_original_chemistry() {
    let mechanism = mechanism(40.0, 500.0, 1e6);
    let mut simulation = uniform(&mechanism);
    let before = simulation.flow().totals();
    let acceleration = 2000.0;
    let report = simulation
        .step_with_body_force(
            Time::new(0.2),
            ReactiveStepControl::default(),
            BodyForce::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    let h = report.elapsed.value();
    assert!(report.coupling_retries > 0);
    assert_eq!(report.elapsed, report.transport.elapsed);
    assert_eq!(simulation.flow().time(), report.elapsed);
    let mass = before.species_masses[..2].iter().sum::<f64>();
    let impulse = mass * acceleration * h;
    let work = impulse * (0.7 + 0.5 * acceleration * h);
    near(report.transport.source_change.impulse, impulse, 3e-12);
    near(report.transport.source_change.work, work, 3e-11);
    near(report.after.momentum - before.momentum, impulse, 3e-12);
    near(report.conservation_residual.momentum, 0.0, 3e-12);
    near(report.conservation_residual.energy, 0.0, 3e-10);
    near(report.absolute_energy_residual_j, 0.0, 3e-8);
    near(report.mass_residual_kg, 0.0, 3e-13);

    // Independent homogeneous reactor: mechanical acceleration must not change
    // reaction progress or temperature through spurious source heating.
    let volume = simulation.flow().mesh().volumes()[0];
    let reference = advance_adiabatic(
        &mechanism,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        report.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for i in 0..simulation.flow().cell_count() {
        let q = simulation.flow().cell(i).unwrap();
        let p = simulation.flow().gas().primitive(q, || false).unwrap();
        near(p.velocity, 0.7 + acceleration * h, 3e-11);
        near(p.temperature, reference.temperature().value(), 2e-4);
        near(
            q.partial_densities[1],
            reference.amounts()[1].value() * R / volume,
            2e-5,
        );
    }
}

fn forced_reacting_mode(cells: usize) -> f64 {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let wave_number = 2.0 * core::f64::consts::PI;
    let flow = field(&mechanism, cells, 0.7, |x| {
        let perturbation = 0.2 * (wave_number * x).sin();
        [0.5 + perturbation, 0.5 - perturbation]
    });
    let mut operator = molecular::molecular(&flow, 0.05, 0.01, 0.1);
    let mut simulation = ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap();
    let end = 0.1;
    let acceleration = 2.0;
    let mut steps = 0;
    while simulation.flow().time().value() < end {
        let report = simulation
            .step_with_dissipation_and_body_force(
                &mut operator,
                Time::new(end - simulation.flow().time().value()),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                BodyForce::UniformAcceleration(acceleration),
                FORCE_CONTROL,
                || false,
            )
            .unwrap();
        let halves = report.dissipation.unwrap();
        assert_eq!(halves.first.elapsed.value(), 0.5 * report.elapsed.value());
        assert_eq!(halves.second.elapsed.value(), 0.5 * report.elapsed.value());
        near(report.conservation_residual.momentum, 0.0, 2e-12);
        near(report.absolute_energy_residual_j, 0.0, 3e-11);
        near(report.mass_residual_kg, 0.0, 3e-13);
        steps += 1;
        assert!(steps < 1000);
    }
    // Isomers have equal caloric laws: the independent solution is a reacting
    // diffusing Fourier mode translated by integral(v dt), including a*t²/2.
    let amplitude = 0.2 * (-(4.0 + 0.05 * wave_number * wave_number) * end).exp();
    let displacement = 0.7 * end + 0.5 * acceleration * end * end;
    let mut error = 0.0;
    for i in 0..cells {
        let x = (i as f64 + 0.5) / cells as f64;
        let expected = 0.5 + amplitude * (wave_number * (x - displacement)).sin();
        let q = simulation.flow().cell(i).unwrap();
        let p = simulation.flow().gas().primitive(q, || false).unwrap();
        error += (q.partial_densities[0] - expected).abs() / cells as f64;
        near(p.velocity, 0.7 + acceleration * end, 3e-9);
        near(p.temperature, 1.0, 3e-8);
    }
    error
}

#[test]
fn g1_accelerating_advection_diffusion_reaction_matches_analytic_motion_and_refines() {
    let coarse = forced_reacting_mode(40);
    let fine = forced_reacting_mode(80);
    assert!(coarse < 1e-3, "coarse={coarse:e}");
    assert!(fine < 0.45 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn g0_zero_force_preserves_unforced_chemistry_and_molecular_paths_bitwise() {
    let mechanism = mechanism(10.0, 30.0, 0.0);
    for molecular in [false, true] {
        let mut reference = uniform(&mechanism);
        let mut forced = uniform(&mechanism);
        let duration = Time::new(0.001);
        let control = ReactiveStepControl::default();
        let (expected, actual) = if molecular {
            let mut a = molecular::molecular(reference.flow(), 0.05, 0.01, 0.1);
            let mut b = molecular::molecular(forced.flow(), 0.05, 0.01, 0.1);
            (
                reference
                    .step_with_dissipation(
                        &mut a,
                        duration,
                        control,
                        DissipationControl::default(),
                        || false,
                    )
                    .unwrap(),
                forced
                    .step_with_dissipation_and_body_force(
                        &mut b,
                        duration,
                        control,
                        DissipationControl::default(),
                        BodyForce::Acceleration(&[0.0; 4]),
                        FORCE_CONTROL,
                        || false,
                    )
                    .unwrap(),
            )
        } else {
            (
                reference.step(duration, control, || false).unwrap(),
                forced
                    .step_with_body_force(
                        duration,
                        control,
                        BodyForce::UniformAcceleration(0.0),
                        FORCE_CONTROL,
                        || false,
                    )
                    .unwrap(),
            )
        };
        assert_eq!(actual, expected);
        assert_eq!(forced.accepted.states, reference.accepted.states);
    }
}

#[test]
fn g4_forced_chemistry_late_cancel_and_malformed_field_preserve_accepted_state() {
    let mechanism = mechanism(10.0, 30.0, 0.0);
    let duration = Time::new(0.001);
    let control = ReactiveStepControl::default();
    let force = BodyForce::Acceleration(&[2.0, -1.0, 3.0, 0.5]);
    let mut reference = uniform(&mechanism);
    let mut polls = 0;
    let expected = reference
        .step_with_body_force(duration, control, force, FORCE_CONTROL, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut simulation = uniform(&mechanism);
    let before = simulation.accepted.states.clone();
    let mut attempted = 0;
    let refusal = simulation
        .step_with_body_force(duration, control, force, FORCE_CONTROL, || {
            attempted += 1;
            attempted >= polls - 2
        })
        .unwrap_err();
    assert_eq!(refusal, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, before);
    assert_eq!(simulation.flow().time().value(), 0.0);
    assert!(
        simulation
            .step_with_body_force(
                duration,
                control,
                BodyForce::Acceleration(&[1.0, 2.0]),
                FORCE_CONTROL,
                || false
            )
            .is_err()
    );
    assert_eq!(simulation.accepted.states, before);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let replay = simulation
        .step_with_body_force(duration, control, force, FORCE_CONTROL, || false)
        .unwrap();
    assert_eq!(replay.elapsed, expected.elapsed);
    for (actual, expected) in simulation
        .accepted
        .states
        .iter()
        .zip(&reference.accepted.states)
    {
        near(*actual, *expected, 2e-11);
    }
}
