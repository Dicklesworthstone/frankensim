use super::*;
use crate::compressible::spatial::mixture::molecular::MolecularStepReport3d;

type Owner = MolecularEuler3d<ConstantTransport>;
const CAP: ForceControl3d = ForceControl3d { maximum_mach: 50.0 };
const ACCELERATION: [f64; 3] = [0.3, -0.2, 0.1];

fn fixture() -> (Owner, ForcingWorkspace3d) {
    let (flow, _) = prepared(
        CartesianGrid::new(2, 2, 1, 0.4, 0.5, 0.6).unwrap(),
        standard(0.2),
        0.6,
        0.02,
        0.03,
    );
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.02, 0.04, 0.01, 0.03),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source).unwrap();
    let molecular = MolecularTransport3d::new(&flow, source, bytes, || false).unwrap();
    let work = flow
        .forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap();
    let bytes = Owner::required_stabilized_workspace_bytes(&flow, &molecular).unwrap()
        + work.workspace_bytes();
    (Owner::new(flow, molecular, bytes, || false).unwrap(), work)
}

fn history(time: f64, side: BoundarySide, face: usize, rho: &mut [f64]) -> ReservoirFaceState {
    let mut result = input(rho, 0.55 + time + 0.01 * face as f64, 0.98 + 2.0 * time);
    result.velocity = [0.08 + time, side.sign() * (0.02 + time), -0.01 + 0.2 * time];
    result
}

#[allow(clippy::too_many_arguments)]
fn run(
    owner: &mut Owner,
    work: &mut ForcingWorkspace3d,
    stabilized: bool,
    duration: Time,
    force: BodyForce3d<'_>,
    maximum_face_evaluations: usize,
    sample: impl FnMut(
        f64,
        BoundarySide,
        usize,
        [f64; 3],
        &mut [f64],
    ) -> Result<ReservoirFaceState, Error>,
    cancelled: impl FnMut() -> bool,
) -> Result<MolecularStepReport3d, DissipationError> {
    if stabilized {
        owner
            .step_stabilized_with_forcing_and_reservoir_profiles(
                duration,
                StabilizedMolecularStepControl::default(),
                force,
                CAP,
                work,
                maximum_face_evaluations,
                sample,
                cancelled,
            )
            .map(|report| report.transport)
    } else {
        owner.step_with_forcing_and_reservoir_profiles(
            duration,
            MolecularStepControl::default(),
            force,
            CAP,
            work,
            maximum_face_evaluations,
            sample,
            cancelled,
        )
    }
}

#[test]
fn g1_forced_molecular_profiles_conserve_species_work_and_split_time_with_zero_force_parity() {
    for stabilized in [false, true] {
        let (mut owner, mut work) = fixture();
        owner
            .step(Time::new(0.002), MolecularStepControl::default(), || false)
            .unwrap();
        let start = owner.flow().time().value();
        let dt = 0.001;
        let mut times = Vec::new();
        let report = run(
            &mut owner,
            &mut work,
            stabilized,
            Time::new(dt),
            BodyForce3d::UniformAcceleration(ACCELERATION),
            if stabilized { 36 } else { 28 },
            |time, side, face, _, rho| {
                if side == BoundarySide::XMin && face == 0 {
                    times.push(time);
                }
                Ok(history(time, side, face, rho))
            },
            || false,
        )
        .unwrap();
        assert_eq!(report.elapsed.value(), dt);
        assert_eq!(report.coupling_retries, 0);
        let offsets: &[f64] = if stabilized {
            &[0.0, 0.0, 1.0 / 6.0, 0.5, 0.0, 1.0, 0.5, 2.0 / 3.0, 1.0]
        } else {
            &[0.0, 0.0, 0.5, 0.0, 1.0, 0.5, 1.0]
        };
        assert_eq!(times.len(), offsets.len());
        for (time, offset) in times.iter().zip(offsets) {
            close(*time, start + offset * dt, 1e-15);
        }
        split_balance(&report);
        for axis in 0..3 {
            close(report.conservation_residual.momentum[axis], 0.0, 2e-11);
            // Each kick uses its own actual density, including species that
            // crossed a prescribed reservoir during the intervening advection.
            close(
                report.source_change.impulse[axis],
                ACCELERATION[axis]
                    * 0.5
                    * dt
                    * (report.advection.before.mass() + report.advection.after.mass()),
                2e-12,
            );
        }
        assert_eq!(report.source_change, report.advection.source_change);
        assert!(report.source_change.work.abs() > 1e-12);
        assert!(report.first.boundary_change.species_masses[0].abs() > 1e-10);
        assert!(report.first.boundary_change.momentum[0].abs() > 1e-10);
        for row in owner.flow().states().chunks_exact(7) {
            assert_eq!(row[2], 0.0);
        }
        for side in SIDES {
            for face in 0..owner.flow().grid().face_count(side) {
                let mut rho = [0.0; 3];
                let expected = history(start + dt, side, face, &mut rho);
                let profile = owner.flow().reservoir_face(side, face).unwrap();
                assert_eq!(profile.partial_densities, rho);
                assert_eq!(profile.velocity, expected.velocity);
                assert_eq!(profile.pressure_pa, expected.pressure_pa);
            }
        }

        let (mut zero, mut zero_work) = fixture();
        let (mut original, _) = fixture();
        let sample = |time, side, face, _, rho: &mut [f64]| Ok(history(time, side, face, rho));
        let expected = if stabilized {
            original
                .step_stabilized_with_reservoir_profiles(
                    Time::new(dt),
                    StabilizedMolecularStepControl::default(),
                    36,
                    sample,
                    || false,
                )
                .unwrap()
                .transport
        } else {
            original
                .step_with_reservoir_profiles(
                    Time::new(dt),
                    MolecularStepControl::default(),
                    28,
                    sample,
                    || false,
                )
                .unwrap()
        };
        let actual = run(
            &mut zero,
            &mut zero_work,
            stabilized,
            Time::new(dt),
            BodyForce3d::UniformAcceleration([0.0; 3]),
            36,
            sample,
            || false,
        )
        .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(snapshot(&zero), snapshot(&original));
    }
}

#[test]
fn g4_forced_molecular_retries_keep_one_budget_and_replay_the_accepted_interval() {
    for stabilized in [false, true] {
        let history = |time: f64, rho: &mut [f64]| {
            let mut result = input(rho, 0.6, 1.0);
            result.velocity[0] = 1e6 * time;
            result
        };
        let force = BodyForce3d::UniformAcceleration(ACCELERATION);
        let (mut owner, mut work) = fixture();
        let mut calls = 0;
        let report = run(
            &mut owner,
            &mut work,
            stabilized,
            Time::new(0.001),
            force,
            10000,
            |time, _, _, _, rho| {
                calls += 1;
                Ok(history(time, rho))
            },
            || false,
        )
        .unwrap();
        assert!(report.coupling_retries > 0);
        assert!(report.elapsed.value() < 0.001);
        split_balance(&report);
        let (mut replay, mut replay_work) = fixture();
        let repeated = run(
            &mut replay,
            &mut replay_work,
            stabilized,
            report.elapsed,
            force,
            10000,
            |time, _, _, _, rho| Ok(history(time, rho)),
            || false,
        )
        .unwrap();
        assert_eq!(owner.flow().states(), replay.flow().states());
        assert_eq!(snapshot(&owner).profiles, snapshot(&replay).profiles);
        assert_eq!(snapshot(&owner).loads, snapshot(&replay).loads);
        assert_eq!(report.source_change, repeated.source_change);
        let (mut limited, mut limited_work) = fixture();
        let before = snapshot(&limited);
        let mut limited_calls = 0;
        let error = run(
            &mut limited,
            &mut limited_work,
            stabilized,
            Time::new(0.001),
            force,
            calls - 1,
            |time, _, _, _, rho| {
                limited_calls += 1;
                Ok(history(time, rho))
            },
            || false,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            DissipationError::Flow(Error::Invalid(
                "mixture reservoir profile face evaluation budget"
            ))
        ));
        assert_eq!(limited_calls, calls - 1);
        assert_eq!(snapshot(&limited), before);
    }
}

#[test]
fn g4_forced_molecular_late_provider_source_budget_and_cancel_refuse_atomically() {
    for stabilized in [false, true] {
        let force = BodyForce3d::UniformAcceleration(ACCELERATION);
        let (mut reference, mut reference_work) = fixture();
        let mut calls = 0;
        let mut polls = 0;
        let expected = run(
            &mut reference,
            &mut reference_work,
            stabilized,
            Time::new(0.001),
            force,
            1000,
            |time, side, face, _, rho| {
                calls += 1;
                Ok(history(time, side, face, rho))
            },
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
        for failure in 0..4 {
            let (mut owner, mut work) = fixture();
            let before = snapshot(&owner);
            let mut observed_calls = 0;
            let mut observed_polls = 0;
            let error = run(
                &mut owner,
                &mut work,
                stabilized,
                Time::new(0.001),
                force,
                if failure == 0 { calls - 1 } else { 1000 },
                |time, side, face, _, rho| {
                    observed_calls += 1;
                    let result = history(time, side, face, rho);
                    if observed_calls == calls {
                        if failure == 1 {
                            return Err(Error::Invalid("late forced bath refusal"));
                        }
                        if failure == 2 {
                            rho[0] = f64::NAN;
                        }
                    }
                    Ok(result)
                },
                || {
                    observed_polls += 1;
                    failure == 3 && observed_polls == polls
                },
            )
            .unwrap_err();
            if failure == 3 {
                assert!(matches!(error, DissipationError::Cancelled), "{error:?}");
            } else {
                assert!(matches!(error, DissipationError::Flow(_)), "{error:?}");
                assert_eq!(observed_calls, if failure == 0 { calls - 1 } else { calls });
            }
            assert_eq!(snapshot(&owner), before);
            // Reusing trial buffers after terminal refusal must reproduce a
            // fresh forced timed transaction and publish only its own loads.
            let actual = run(
                &mut owner,
                &mut work,
                stabilized,
                Time::new(0.001),
                force,
                1000,
                |time, side, face, _, rho| Ok(history(time, side, face, rho)),
                || false,
            )
            .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(snapshot(&owner), snapshot(&reference));
        }
    }
}
