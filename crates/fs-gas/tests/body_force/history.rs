//! Time-resolved physical forcing through the original native evolution paths.

use super::*;

fn cells(flow: &MixtureEuler1d) -> Vec<(Vec<f64>, f64, f64)> {
    (0..flow.cell_count())
        .map(|i| {
            let q = flow.cell(i).unwrap();
            (q.partial_densities.to_vec(), q.momentum, q.energy)
        })
        .collect()
}

#[test]
fn linear_acceleration_uses_absolute_native_time_and_preserves_calorics() {
    let (reference_time, acceleration, jerk) = (-0.25, -10.0, 300.0);
    let source = BodyForce::UniformAccelerationRamp {
        reference_time,
        acceleration,
        jerk,
    };
    let mut frozen = frozen(-0.3);
    let mut mixture = mixture(-0.3);
    let initial_frozen = frozen.totals();
    let initial_mixture = mixture.totals();
    let (mut frozen_impulse, mut mixture_impulse) = (0.0, 0.0);
    let (mut frozen_work, mut mixture_work) = (0.0, 0.0);
    let h = 0.0005;
    for _ in 0..80 {
        let start = frozen.time().value();
        assert_eq!(mixture.time(), frozen.time());
        let delta = h * (acceleration + jerk * (start - reference_time) + 0.5 * jerk * h);
        let a = frozen
            .step_with_body_force(Time::new(h), StepControl::default(), source, FORCE, || false)
            .unwrap();
        let b = mixture
            .step_with_body_force(Time::new(h), StepControl::default(), source, FORCE, || false)
            .unwrap();
        assert_eq!(a.elapsed, Time::new(h));
        assert_eq!(b.elapsed, a.elapsed);
        close(a.source_change.impulse, initial_frozen.mass * delta, 2e-13);
        close(b.source_change.impulse, initial_mixture.mass() * delta, 2e-13);
        close(a.conservation_residual.energy, 0.0, 2e-12);
        close(b.conservation_residual.energy, 0.0, 2e-12);
        frozen_impulse += a.source_change.impulse;
        mixture_impulse += b.source_change.impulse;
        frozen_work += a.source_change.work;
        mixture_work += b.source_change.work;
    }
    let time = frozen.time().value();
    close(time, 0.04, 1e-15);
    let velocity = -0.3 + (acceleration - jerk * reference_time) * time
        + 0.5 * jerk * time * time;
    let work_per_mass = 0.5 * (velocity * velocity - 0.09);
    close(frozen_impulse, initial_frozen.mass * (velocity + 0.3), 5e-12);
    close(mixture_impulse, initial_mixture.mass() * (velocity + 0.3), 5e-12);
    close(frozen_work, initial_frozen.mass * work_per_mass, 5e-11);
    close(mixture_work, initial_mixture.mass() * work_per_mass, 5e-11);
    close(frozen.totals().energy - initial_frozen.energy, frozen_work, 5e-11);
    close(mixture.totals().energy - initial_mixture.energy, mixture_work, 5e-11);
    assert_eq!(mixture.totals().species_masses, initial_mixture.species_masses);
    for i in 0..8 {
        let q = frozen.states()[i];
        let p = frozen.gas().primitive(q, || false).unwrap();
        assert_eq!(q.density, 1.0);
        close(p.velocity, velocity, 5e-12);
        close(p.temperature, 1.7, 3e-9);
        let q = mixture.cell(i).unwrap();
        assert_eq!(q.partial_densities, &[0.8, 0.2, 0.0]);
        let p = mixture.gas().primitive(q, || false).unwrap();
        close(p.velocity, velocity, 5e-12);
        close(p.temperature, 3.0, 3e-9);
    }
}

#[test]
fn cell_force_density_histories_use_each_actual_mass_and_reference_time() {
    let densities = [1.0, 2.0];
    let partials = [[0.8, 0.2, 0.0], [0.4, 1.6, 0.0]];
    let velocities = [-0.6, 0.7];
    let gas = frozen_gas();
    let states = (0..2)
        .map(|i| gas.conserved(densities[i], velocities[i], 3.0 * densities[i]).unwrap())
        .collect();
    let mut frozen = Euler1d::new(
        gas, Length::new(0.4), 0.3, states, Boundaries::Periodic, || false,
    )
    .unwrap();
    let gas = mixture_gas();
    let states = (0..2)
        .map(|i| {
            let d = partials[i];
            gas.conserved(&d, velocities[i], 3.0 * (d[0] + 0.5 * d[1])).unwrap()
        })
        .collect();
    let mut mixture = MixtureEuler1d::new(
        gas, Length::new(0.4), 0.3, states, Boundaries::Periodic, || false,
    )
    .unwrap();
    let before_frozen = frozen.totals();
    let before_mixture = mixture.totals();
    let force_density = [4.0, -3.0];
    let rate = [120.0, -80.0];
    let reference_time = -0.03;
    let h = 0.04;
    let source = BodyForce::ForceDensityRamp {
        reference_time,
        force_density: &force_density,
        rate: &rate,
    };
    let a = frozen.apply_body_force(Time::new(h), source, FORCE, || false).unwrap();
    let b = mixture.apply_body_force(Time::new(h), source, FORCE, || false).unwrap();
    let (mut impulse, mut work) = (0.0, 0.0);
    for i in 0..2 {
        let integral = h * (force_density[i] - rate[i] * reference_time + 0.5 * rate[i] * h);
        let velocity = velocities[i] + integral / densities[i];
        impulse += 0.12 * integral;
        work += 0.06 * densities[i] * (velocity * velocity - velocities[i] * velocities[i]);
        let q = frozen.states()[i];
        assert_eq!(q.density, densities[i]);
        let p = frozen.gas().primitive(q, || false).unwrap();
        close(p.velocity, velocity, 2e-13);
        close(p.temperature, 3.0, 3e-10);
        let q = mixture.cell(i).unwrap();
        assert_eq!(q.partial_densities, &partials[i]);
        let p = mixture.gas().primitive(q, || false).unwrap();
        close(p.velocity, velocity, 2e-13);
        close(p.temperature, 3.0, 3e-10);
    }
    for exchange in [a, b] {
        close(exchange.impulse, impulse, 2e-14);
        close(exchange.work, work, 2e-13);
    }
    close(frozen.totals().energy - before_frozen.energy - a.work, a.energy_residual, 2e-13);
    close(mixture.totals().energy - before_mixture.energy - b.work, b.energy_residual, 2e-13);
    assert_eq!(mixture.totals().species_masses, before_mixture.species_masses);
    assert_eq!(frozen.time(), Time::new(0.0));
    assert_eq!(mixture.time(), Time::new(0.0));
}

#[test]
fn all_zero_rate_histories_replay_the_held_native_paths_bit_for_bit() {
    let values = [1.5, -0.2, 0.7, -0.1, 0.4, 0.3, -0.9, 0.8];
    let zeros = [0.0; 8];
    for (held, history) in [
        (
            BodyForce::UniformAcceleration(1.5),
            BodyForce::UniformAccelerationRamp {
                reference_time: -1e200,
                acceleration: 1.5,
                jerk: 0.0,
            },
        ),
        (
            BodyForce::Acceleration(&values),
            BodyForce::AccelerationRamp {
                reference_time: 1e200,
                acceleration: &values,
                jerk: &zeros,
            },
        ),
        (
            BodyForce::UniformForceDensity(-0.7),
            BodyForce::UniformForceDensityRamp {
                reference_time: -1e200,
                force_density: -0.7,
                rate: 0.0,
            },
        ),
        (
            BodyForce::ForceDensity(&values),
            BodyForce::ForceDensityRamp {
                reference_time: 1e200,
                force_density: &values,
                rate: &zeros,
            },
        ),
    ] {
        let mut original_frozen = frozen(0.3);
        let mut history_frozen = frozen(0.3);
        let mut original_mixture = mixture(0.3);
        let mut history_mixture = mixture(0.3);
        for _ in 0..3 {
            let expected = original_frozen
                .step_with_body_force(Time::new(0.001), StepControl::default(), held, FORCE, || false)
                .unwrap();
            let actual = history_frozen
                .step_with_body_force(Time::new(0.001), StepControl::default(), history, FORCE, || false)
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(history_frozen.states(), original_frozen.states());
            assert_eq!(history_frozen.time(), original_frozen.time());
            let expected = original_mixture
                .step_with_body_force(Time::new(0.001), StepControl::default(), held, FORCE, || false)
                .unwrap();
            let actual = history_mixture
                .step_with_body_force(Time::new(0.001), StepControl::default(), history, FORCE, || false)
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(cells(&history_mixture), cells(&original_mixture));
            assert_eq!(history_mixture.time(), original_mixture.time());
        }
    }
}

#[test]
fn zero_net_impulse_does_not_hide_an_interior_mach_violation_or_publish_scratch() {
    let source = BodyForce::UniformAccelerationRamp {
        reference_time: 0.0,
        acceleration: 10_000.0,
        jerk: -1_000_000.0,
    };
    let limit = ForceControl { maximum_mach: 5.0 };
    let h = Time::new(0.02);
    // The endpoint velocity returns to zero, but u(0.01) = 50 m/s.
    let mut frozen_flow = frozen(0.0);
    let initial_frozen = frozen_flow.states().to_vec();
    let mut mixture_flow = mixture(0.0);
    let initial_mixture = cells(&mixture_flow);
    assert!(frozen_flow.apply_body_force(h, source, limit, || false).is_err());
    assert!(mixture_flow.apply_body_force(h, source, limit, || false).is_err());
    assert_eq!(frozen_flow.states(), initial_frozen);
    assert_eq!(cells(&mixture_flow), initial_mixture);
    assert_eq!(frozen_flow.time(), Time::new(0.0));
    assert_eq!(mixture_flow.time(), Time::new(0.0));

    let mut reference_frozen = frozen(0.0);
    let mut reference_mixture = mixture(0.0);
    let a = frozen_flow.apply_body_force(h, source, FORCE, || false).unwrap();
    let expected = reference_frozen.apply_body_force(h, source, FORCE, || false).unwrap();
    assert_eq!(a, expected);
    assert_eq!(frozen_flow.states(), reference_frozen.states());
    let b = mixture_flow.apply_body_force(h, source, FORCE, || false).unwrap();
    let expected = reference_mixture.apply_body_force(h, source, FORCE, || false).unwrap();
    assert_eq!(b, expected);
    assert_eq!(cells(&mixture_flow), cells(&reference_mixture));
    for exchange in [a, b] {
        close(exchange.impulse, 0.0, 2e-13);
        close(exchange.work, 0.0, 2e-13);
    }
}

#[test]
fn interior_source_cfl_shortens_the_whole_step_and_replays_the_accepted_clock() {
    let source = BodyForce::UniformAccelerationRamp {
        reference_time: 0.0,
        acceleration: 20_000.0,
        jerk: -40_000_000.0,
    };
    let mut flow = frozen(20.0);
    let initial = flow.totals();
    let sound_speed = flow.gas().primitive(flow.states()[0], || false).unwrap().sound_speed;
    let controls = StepControl::default();
    let requested = 0.002;
    // Both endpoints of the first half kick are 20 m/s; its interior
    // maximum is 25 m/s. Only that interior state violates the native CFL.
    assert!(requested * (20.0 + sound_speed) / 0.125 < controls.cfl);
    assert!(requested * (25.0 + sound_speed) / 0.125 > controls.cfl);
    let report = flow
        .step_with_body_force(Time::new(requested), controls, source, FORCE, || false)
        .unwrap();
    assert!(report.halvings > 0);
    assert_eq!(report.elapsed, Time::new(0.001));
    assert!(!report.first_order_retry);
    assert!(report.maximum_cfl <= controls.cfl * (1.0 + 64.0 * f64::EPSILON));
    let h = report.elapsed.value();
    let velocity = 20.0 + 20_000.0 * h - 20_000_000.0 * h * h;
    close(report.source_change.impulse, initial.mass * (velocity - 20.0), 5e-12);
    close(report.source_change.work, 0.5 * initial.mass * (velocity * velocity - 400.0), 5e-11);
    close(report.conservation_residual.energy, 0.0, 5e-11);
    for &q in flow.states() {
        let p = flow.gas().primitive(q, || false).unwrap();
        close(p.velocity, velocity, 2e-12);
        close(p.temperature, 1.7, 5e-9);
    }
    let mut reference = frozen(20.0);
    let replay = reference
        .step_with_body_force(report.elapsed, controls, source, FORCE, || false)
        .unwrap();
    assert_eq!(replay.halvings, 0);
    assert_eq!(replay.source_change, report.source_change);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(flow.time(), reference.time());
}
