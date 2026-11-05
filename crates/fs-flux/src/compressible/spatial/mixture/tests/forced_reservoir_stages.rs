use super::*;

const CAP: ForceControl3d = ForceControl3d { maximum_mach: 50.0 };

fn force_work(flow: &MixtureEuler3d) -> ForcingWorkspace3d {
    flow.forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap()
}

fn history(time: f64, face: usize, rho: &mut [f64]) -> ReservoirFaceState {
    rho.copy_from_slice(&[0.65 + 0.02 * face as f64 + time, 0.35 - 0.4 * time, 0.0]);
    ReservoirFaceState {
        velocity: [0.08 + 2.0 * time, 0.03 - time, -0.02 + 0.5 * time],
        pressure_pa: 1.02 + time,
    }
}

fn fixture() -> MixtureEuler3d {
    prepared(
        CartesianGrid::new(2, 2, 2, 0.3, 0.4, 0.5).unwrap(),
        open_axis(0),
        &[BoundarySide::XMin],
        [0.08, 0.03, -0.02],
    )
}

/// Independent cellwise momentum impulse and kinetic-energy difference. The
/// density is recovered again after advection, including exchanged species.
fn analytic_kick(
    flow: &mut MixtureEuler3d,
    field: &[[f64; 3]],
    force_density: bool,
    duration: f64,
) -> ForceChange3d {
    let mut source = ForceChange3d::default();
    let volume = flow.grid.cell_volume();
    for (row, force) in flow.states.chunks_exact_mut(7).zip(field) {
        let rho: f64 = row[..3].iter().sum();
        let momentum = [row[3], row[4], row[5]];
        let impulse = force.map(|value| value * duration * if force_density { 1.0 } else { rho });
        let mut work = 0.0;
        for axis in 0..3 {
            work += impulse[axis] * (momentum[axis] + 0.5 * impulse[axis]) / rho;
            row[3 + axis] += impulse[axis];
            source.impulse[axis] += volume * impulse[axis];
        }
        row[6] += work;
        source.work += volume * work;
    }
    source
}

#[test]
fn g1_forced_profiles_match_independent_kinetic_work_and_canonical_species() {
    for force_density in [false, true] {
        let mut actual = fixture();
        let mut expected = fixture();
        let mut workspace = force_work(&actual);
        let forcing: Vec<_> = (0..actual.cell_count())
            .map(|i| [2.0 + i as f64 * 0.1, -1.5, 0.7])
            .collect();
        let force = if force_density {
            BodyForce3d::ForceDensity(&forcing)
        } else {
            BodyForce3d::Acceleration(&forcing)
        };
        let duration = 0.002;
        let first = analytic_kick(&mut expected, &forcing, force_density, 0.5 * duration);
        let reference = expected
            .step_with_reservoir_profiles(
                Time::new(duration),
                StepControl::default(),
                8,
                |time, _, face, _, rho| Ok(history(time, face, rho)),
                || false,
            )
            .unwrap();
        assert_eq!(reference.elapsed.value(), duration);
        let second = analytic_kick(&mut expected, &forcing, force_density, 0.5 * duration);
        let report = actual
            .step_with_forcing_and_reservoir_profiles(
                Time::new(duration),
                StepControl::default(),
                force,
                CAP,
                &mut workspace,
                8,
                |time, _, face, _, rho| Ok(history(time, face, rho)),
                || false,
            )
            .unwrap();
        assert_eq!(report.elapsed.value(), duration);
        assert_eq!(report.halvings, 0);
        balance(&report);
        for (&actual, &expected) in actual.states().iter().zip(expected.states()) {
            close(actual, expected, 3e-12);
        }
        for axis in 0..3 {
            close(
                report.source_change.impulse[axis],
                first.impulse[axis] + second.impulse[axis],
                3e-13,
            );
        }
        close(report.source_change.work, first.work + second.work, 3e-13);
        for row in actual.states().chunks_exact(7) {
            assert_eq!(row[2], 0.0);
        }
        for (&actual, &expected) in actual
            .boundary_impulses
            .iter()
            .zip(&expected.boundary_impulses)
        {
            close(actual, expected, 3e-12);
        }
        assert_eq!(snapshot(&actual).profiles, snapshot(&expected).profiles);
    }
}

#[test]
fn g0_forced_profiles_preserve_zero_force_and_constant_history_bits() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        let mut zero = fixture();
        let mut unforced = fixture();
        let mut workspace = force_work(&zero);
        let bytes = workspace.workspace_bytes();
        let expected = unforced
            .step_with_reservoir_profiles(
                Time::new(0.002),
                controls,
                8,
                |time, _, face, _, rho| Ok(history(time, face, rho)),
                || false,
            )
            .unwrap();
        let actual = zero
            .step_with_forcing_and_reservoir_profiles(
                Time::new(0.002),
                controls,
                BodyForce3d::UniformAcceleration([0.0; 3]),
                CAP,
                &mut workspace,
                8,
                |time, _, face, _, rho| Ok(history(time, face, rho)),
                || false,
            )
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(snapshot(&zero), snapshot(&unforced));
        assert_eq!(workspace.workspace_bytes(), bytes);

        let mut held = fixture();
        let mut timed = fixture();
        let mut held_work = force_work(&held);
        let mut timed_work = force_work(&timed);
        let force = BodyForce3d::UniformAcceleration([0.2, -0.1, 0.05]);
        let reference = held
            .step_with_forcing(
                Time::new(0.002),
                controls,
                force,
                CAP,
                &mut held_work,
                || false,
            )
            .unwrap();
        let report = timed
            .step_with_forcing_and_reservoir_profiles(
                Time::new(0.002),
                controls,
                force,
                CAP,
                &mut timed_work,
                8,
                |_, _, _, _, rho| Ok(input(rho, [0.08, 0.03, -0.02], 1.0)),
                || false,
            )
            .unwrap();
        assert_eq!(report, reference);
        assert_eq!(snapshot(&timed), snapshot(&held));
    }
}

#[test]
fn g4_forced_profile_shortening_replays_both_kicks_with_one_provider_budget() {
    let force = BodyForce3d::UniformAcceleration([0.3, 0.2, -0.1]);
    let mut flow = retry_fixture();
    let mut workspace = force_work(&flow);
    let start = flow.time().value();
    let mut calls = Vec::new();
    let report = flow
        .step_with_forcing_and_reservoir_profiles(
            Time::new(0.02),
            StepControl::default(),
            force,
            CAP,
            &mut workspace,
            4,
            |time, _, _, _, rho| {
                calls.push(time);
                Ok(input(rho, [4.0 + 1000.0 * (time - start), 0.0, 0.0], 1.0))
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.halvings, 1);
    assert_eq!(report.elapsed.value(), 0.01);
    assert_eq!(calls, [start, start + 0.02, start, start + 0.01]);
    balance(&report);
    let mut replay = retry_fixture();
    let mut replay_work = force_work(&replay);
    let repeated = replay
        .step_with_forcing_and_reservoir_profiles(
            report.elapsed,
            StepControl::default(),
            force,
            CAP,
            &mut replay_work,
            2,
            |time, _, _, _, rho| Ok(input(rho, [4.0 + 1000.0 * (time - start), 0.0, 0.0], 1.0)),
            || false,
        )
        .unwrap();
    assert_eq!(flow.states(), replay.states());
    assert_eq!(flow.boundary_impulses, replay.boundary_impulses);
    assert_eq!(report.source_change, repeated.source_change);
    assert_eq!(snapshot(&flow).profiles, snapshot(&replay).profiles);

    let mut limited = retry_fixture();
    let before = snapshot(&limited);
    let mut limited_work = force_work(&limited);
    let mut calls = 0;
    let error = limited
        .step_with_forcing_and_reservoir_profiles(
            Time::new(0.02),
            StepControl::default(),
            force,
            CAP,
            &mut limited_work,
            3,
            |time, _, _, _, rho| {
                calls += 1;
                Ok(input(rho, [4.0 + 1000.0 * (time - start), 0.0, 0.0], 1.0))
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error,
        Error::Invalid("mixture reservoir profile face evaluation budget")
    ));
    assert_eq!(calls, 3);
    assert_eq!(snapshot(&limited), before);
}
