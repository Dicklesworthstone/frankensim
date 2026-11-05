//! Original force work and time-dependent open-face transport in one transaction.

use super::*;

const SIDES: [BoundarySide; 2] = [BoundarySide::XMin, BoundarySide::XMax];
const FACES: usize = 8;

fn bath(time: f64, side: BoundarySide, face: usize) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.05 + 0.01 * face as f64 + 0.02 * (4.0 * time).sin(),
        velocity: [0.18 + 0.2 * time, 0.03 * side.sign(), -0.01],
        pressure_pa: 2.1 + 0.3 * time,
    }
}

fn profiled() -> Euler3d {
    let source = gas(0.2);
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let states = (0..grid.len())
        .map(|i| {
            State::from_primitive(
                &source,
                1.0 + 0.03 * i as f64,
                [0.15 + 0.01 * i as f64, -0.04, 0.02],
                2.0 + 0.02 * i as f64,
            )
            .unwrap()
        })
        .collect();
    let mut flow = Euler3d::new(
        source,
        grid,
        states,
        Boundaries3d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Transmissive,
                upper: BoundaryCondition::Transmissive,
            },
            ..PERIODIC
        },
        || false,
    )
    .unwrap();
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&SIDES)
        .unwrap();
    flow.prepare_reservoir_profiles(
        &SIDES,
        bytes,
        |side, face, _| Ok(bath(0.0, side, face)),
        || false,
    )
    .unwrap();
    flow
}

// Independent kinetic-energy identity, without calling the source kernel.
fn kick_by_kinetic_energy(
    flow: &mut Euler3d,
    values: &[[f64; 3]],
    acceleration: bool,
    duration: f64,
) -> ForceChange3d {
    let mut change = ForceChange3d::default();
    let volume = flow.grid.cell_volume();
    for (q, load) in flow.states.iter_mut().zip(values) {
        let old = *q;
        for axis in 0..3 {
            q.momentum[axis] += duration * load[axis] * if acceleration { q.density } else { 1.0 };
            change.impulse[axis] += (q.momentum[axis] - old.momentum[axis]) * volume;
        }
        let work = (0..3)
            .map(|axis| (q.momentum[axis].powi(2) - old.momentum[axis].powi(2)) / (2.0 * q.density))
            .sum::<f64>();
        q.energy += work;
        change.work += work * volume;
    }
    change
}

#[test]
fn g1_timed_forcing_matches_independent_kicks_and_original_open_face_quadrature() {
    for acceleration in [false, true] {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut flow = profiled();
            let mut reference = profiled();
            let loads: Vec<_> = (0..flow.states.len())
                .map(|i| [0.7 + 0.02 * i as f64, -0.2, 0.15])
                .collect();
            let force = if acceleration {
                BodyForce3d::Acceleration(&loads)
            } else {
                BodyForce3d::ForceDensity(&loads)
            };
            let controls = StepControl {
                reconstruction,
                ..StepControl::default()
            };
            let mut work = flow.forcing_workspace(usize::MAX).unwrap();
            let storage = work.original.as_ptr();
            let report = flow
                .step_with_forcing_and_reservoir_profiles(
                    Time::new(0.002),
                    controls,
                    force,
                    cap(),
                    &mut work,
                    2 * FACES,
                    |time, side, face, _| Ok(bath(time, side, face)),
                    || false,
                )
                .unwrap();
            assert_eq!(report.elapsed.value(), 0.002);
            assert_eq!(report.halvings, 0);
            let first = kick_by_kinetic_energy(&mut reference, &loads, acceleration, 0.001);
            let central = reference
                .step_with_reservoir_profiles(
                    report.elapsed,
                    controls,
                    2 * FACES,
                    |time, side, face, _| Ok(bath(time, side, face)),
                    || false,
                )
                .unwrap();
            let second = kick_by_kinetic_energy(&mut reference, &loads, acceleration, 0.001);
            for (actual, expected) in flow.states().iter().zip(reference.states()) {
                for (a, b) in actual.array().into_iter().zip(expected.array()) {
                    close(a, b, 3e-12);
                }
            }
            for axis in 0..3 {
                close(
                    report.source_change.impulse[axis],
                    first.impulse[axis] + second.impulse[axis],
                    3e-13,
                );
            }
            close(report.source_change.work, first.work + second.work, 3e-13);
            assert!(report.source_change.work > 0.0);
            for (actual, expected) in report
                .boundary_change
                .array()
                .into_iter()
                .zip(central.boundary_change.array())
            {
                close(actual, expected, 3e-13);
            }
            for side in SIDES {
                for face in 0..flow.grid.face_count(side) {
                    assert_eq!(
                        flow.reservoir_face(side, face),
                        Some(bath(0.002, side, face))
                    );
                }
                for (actual, expected) in flow
                    .boundary_impulses(side)
                    .iter()
                    .zip(reference.boundary_impulses(side))
                {
                    for (a, b) in actual.array().into_iter().zip(expected.array()) {
                        close(a, b, 3e-13);
                    }
                }
            }
            assert_eq!(storage, work.original.as_ptr());
            balance(&report);
        }
    }
}

#[test]
fn g0_zero_timed_force_retains_cells_profiles_loads_and_report_bits() {
    for force in [
        BodyForce3d::UniformAcceleration([0.0; 3]),
        BodyForce3d::UniformForceDensity([-0.0; 3]),
    ] {
        let mut flow = profiled();
        let mut reference = profiled();
        let mut work = flow.forcing_workspace(usize::MAX).unwrap();
        for _ in 0..2 {
            let expected = reference
                .step_with_reservoir_profiles(
                    Time::new(0.002),
                    StepControl::default(),
                    2 * FACES,
                    |time, side, face, _| Ok(bath(time, side, face)),
                    || false,
                )
                .unwrap();
            let actual = flow
                .step_with_forcing_and_reservoir_profiles(
                    Time::new(0.002),
                    StepControl::default(),
                    force,
                    cap(),
                    &mut work,
                    2 * FACES,
                    |time, side, face, _| Ok(bath(time, side, face)),
                    || false,
                )
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(flow.states(), reference.states());
            assert_eq!(flow.time(), reference.time());
            assert_eq!(impulses(&flow), impulses(&reference));
            for side in SIDES {
                for face in 0..flow.grid.face_count(side) {
                    assert_eq!(
                        flow.reservoir_face(side, face),
                        reference.reservoir_face(side, face)
                    );
                }
            }
        }
    }
}

#[test]
fn g2_timed_force_mach_retry_replays_accepted_cells_with_one_source_budget() {
    let controls = StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    };
    let force = BodyForce3d::UniformAcceleration([80.0, 0.0, 0.0]);
    let cap = ForceControl3d { maximum_mach: 0.25 };
    let mut flow = profiled();
    let mut work = flow.forcing_workspace(usize::MAX).unwrap();
    let mut times = Vec::new();
    let mut calls = 0;
    let report = flow
        .step_with_forcing_and_reservoir_profiles(
            Time::new(0.01),
            controls,
            force,
            cap,
            &mut work,
            1000,
            |time, side, face, _| {
                calls += 1;
                if side == BoundarySide::XMin && face == 0 {
                    times.push(time);
                }
                Ok(bath(time, side, face))
            },
            || false,
        )
        .unwrap();
    assert!(report.halvings > 0);
    assert!(times.iter().filter(|&&t| t == 0.0).count() > 1);
    assert_eq!(times.last(), Some(&report.elapsed.value()));
    for q in flow.states() {
        let p = q.primitive(flow.gas(), || false).unwrap();
        let mach = p
            .velocity
            .iter()
            .map(|u| (u / p.sound_speed).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(mach <= cap.maximum_mach);
    }
    let mut direct = profiled();
    let mut direct_work = direct.forcing_workspace(usize::MAX).unwrap();
    let shorter = direct
        .step_with_forcing_and_reservoir_profiles(
            report.elapsed,
            controls,
            force,
            cap,
            &mut direct_work,
            1000,
            |time, side, face, _| Ok(bath(time, side, face)),
            || false,
        )
        .unwrap();
    assert_eq!(shorter.halvings, 0);
    assert_eq!(flow.states(), direct.states());
    assert_eq!(report.source_change, shorter.source_change);
    assert_eq!(impulses(&flow), impulses(&direct));
    balance(&report);
    let mut refused = profiled();
    let original = refused.states().to_vec();
    let mut refused_work = refused.forcing_workspace(usize::MAX).unwrap();
    let mut used = 0;
    assert_eq!(
        refused.step_with_forcing_and_reservoir_profiles(
            Time::new(0.01),
            controls,
            force,
            cap,
            &mut refused_work,
            calls - FACES,
            |time, side, face, _| {
                used += 1;
                Ok(bath(time, side, face))
            },
            || false,
        ),
        Err(Error::Invalid(
            "timed reservoir profile face evaluation budget exhausted"
        ))
    );
    assert_eq!(used, calls - FACES);
    assert_eq!(refused.states(), original);
    assert_eq!(refused.time().value(), 0.0);
    assert!(refused.last_step().is_none());
}

#[test]
fn g4_timed_force_late_refusal_preserves_previous_publication_and_held_recovery() {
    let force = BodyForce3d::UniformForceDensity([0.7, -0.2, 0.15]);
    let make = || {
        let mut flow = profiled();
        let mut work = flow.forcing_workspace(usize::MAX).unwrap();
        flow.step_with_forcing_and_reservoir_profiles(
            Time::new(0.002),
            StepControl::default(),
            force,
            cap(),
            &mut work,
            2 * FACES,
            |time, side, face, _| Ok(bath(time, side, face)),
            || false,
        )
        .unwrap();
        (flow, work)
    };
    let (mut probe, mut probe_work) = make();
    let mut total_polls = 0;
    probe
        .step_with_forcing_and_reservoir_profiles(
            Time::new(0.002),
            StepControl::default(),
            force,
            cap(),
            &mut probe_work,
            2 * FACES,
            |time, side, face, _| Ok(bath(time, side, face)),
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    for failure in 0..3 {
        let (mut flow, mut work) = make();
        let (mut reference, mut reference_work) = make();
        let states = flow.states().to_vec();
        let previous = flow.last_step().copied();
        let loads = impulses(&flow);
        let mut calls = 0;
        let mut polls = 0;
        let result = flow.step_with_forcing_and_reservoir_profiles(
            Time::new(0.002),
            StepControl::default(),
            force,
            cap(),
            &mut work,
            if failure == 0 { FACES } else { 2 * FACES },
            |time, side, face, _| {
                calls += 1;
                if failure == 1 && calls == FACES + 2 {
                    Err(Error::Invalid("late forced trace unavailable"))
                } else {
                    Ok(bath(time, side, face))
                }
            },
            || {
                polls += 1;
                failure == 2 && polls == total_polls
            },
        );
        assert_eq!(
            result,
            Err(match failure {
                0 => Error::Invalid("timed reservoir profile face evaluation budget exhausted"),
                1 => Error::Invalid("late forced trace unavailable"),
                _ => Error::Cancelled,
            })
        );
        assert_eq!(flow.states(), states);
        assert_eq!(flow.time().value(), 0.002);
        assert_eq!(flow.last_step().copied(), previous);
        assert_eq!(impulses(&flow), loads);
        for side in SIDES {
            for face in 0..flow.grid.face_count(side) {
                assert_eq!(
                    flow.reservoir_face(side, face),
                    Some(bath(0.002, side, face))
                );
            }
        }
        let actual = flow
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        let expected = reference
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                force,
                cap(),
                &mut reference_work,
                || false,
            )
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(flow.states(), reference.states());
        assert_eq!(impulses(&flow), impulses(&reference));
    }
}
