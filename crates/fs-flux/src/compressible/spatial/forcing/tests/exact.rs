use super::*;

#[test]
fn g1_exact_coriolis_resolves_many_inertial_turns_on_the_native_euler_path() {
    let velocity = [0.23, -0.17, 0.11];
    let mut flow = field([3; 3], [0.3; 3], 0.3, |_| {
        [1.7, velocity[0], velocity[1], velocity[2], 3.4]
    });
    let initial = flow.states()[0].primitive(flow.gas(), || false).unwrap();
    let energy = flow.states()[0].energy;
    let mut workspace = flow.forcing_workspace(usize::MAX).unwrap();
    let pointer = workspace.original.as_ptr();
    let omega = [[1000.0, 2000.0, 2000.0]];
    let source = BodyForce3d::CoriolisExact {
        angular_velocity: &omega,
        acceleration: &[[0.0; 3]],
    };
    for _ in 0..64 {
        let report = flow
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                source,
                cap(),
                &mut workspace,
                || false,
            )
            .unwrap();
        assert_eq!(report.source_change.work, 0.0);
        balance(&report);
    }
    assert_eq!(workspace.original.as_ptr(), pointer);
    let duration = flow.time().value();
    // A full source interval turns through 6000*dt, independent of the two
    // source half steps and the original advective RK stages between them.
    assert!(6000.0 * duration > 100.0 * std::f64::consts::TAU);
    let (sine, cosine) = (-6000.0 * duration).sin_cos();
    let axis = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let projection = (velocity[0] + 2.0 * velocity[1] + 2.0 * velocity[2]) / 3.0;
    let cross = [
        (2.0 * velocity[2] - 2.0 * velocity[1]) / 3.0,
        (2.0 * velocity[0] - velocity[2]) / 3.0,
        (velocity[1] - 2.0 * velocity[0]) / 3.0,
    ];
    for q in flow.states() {
        let p = q.primitive(flow.gas(), || false).unwrap();
        assert_eq!(q.energy.to_bits(), energy.to_bits());
        close(p.temperature, initial.temperature, 5e-11);
        close(p.pressure, initial.pressure, 5e-11);
        for a in 0..3 {
            let expected =
                cosine * velocity[a] + sine * cross[a] + (1.0 - cosine) * projection * axis[a];
            close(p.velocity[a], expected, 4e-11);
        }
    }
}

#[test]
fn g1_exact_coriolis_integrated_drive_supplies_native_impulse_and_kinetic_work() {
    let rho = 1.3;
    let velocity = [0.23, -0.17, 0.11];
    let acceleration = [[0.4, -0.6, 0.9]];
    let omega = [[0.0, 0.0, 200.0]];
    let mut flow = field([3; 3], [0.3; 3], 0.3, |_| {
        [rho, velocity[0], velocity[1], velocity[2], 2.0 * rho]
    });
    let initial = flow.states()[0].primitive(flow.gas(), || false).unwrap();
    let before = flow.totals();
    let mut workspace = flow.forcing_workspace(usize::MAX).unwrap();
    let mut source_work = 0.0;
    let mut source_impulse = [0.0; 3];
    for _ in 0..16 {
        let report = flow
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                BodyForce3d::CoriolisExact {
                    angular_velocity: &omega,
                    acceleration: &acceleration,
                },
                cap(),
                &mut workspace,
                || false,
            )
            .unwrap();
        source_work += report.source_change.work;
        for a in 0..3 {
            source_impulse[a] += report.source_change.impulse[a];
        }
        balance(&report);
    }
    let duration = flow.time().value();
    let (sine, cosine) = (400.0 * duration).sin_cos();
    let acceleration = acceleration[0];
    let expected = [
        cosine * velocity[0]
            + sine * velocity[1]
            + (sine * acceleration[0] + (1.0 - cosine) * acceleration[1]) / 400.0,
        -sine * velocity[0]
            + cosine * velocity[1]
            + (-(1.0 - cosine) * acceleration[0] + sine * acceleration[1]) / 400.0,
        velocity[2] + duration * acceleration[2],
    ];
    let work = 0.5
        * before.mass
        * (0..3)
            .map(|a| expected[a] * expected[a] - velocity[a] * velocity[a])
            .sum::<f64>();
    close(source_work, work, 2e-12);
    close(flow.totals().energy - before.energy, work, 2e-12);
    for a in 0..3 {
        close(
            source_impulse[a],
            before.mass * (expected[a] - velocity[a]),
            2e-12,
        );
    }
    for q in flow.states() {
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(p.temperature, initial.temperature, 5e-11);
        close(p.pressure, initial.pressure, 5e-11);
        for a in 0..3 {
            close(p.velocity[a], expected[a], 2e-12);
        }
    }
}
