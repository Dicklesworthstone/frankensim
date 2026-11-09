use super::*;
use crate::compressible::spatial::{BodyForce3d, ForceControl3d, MeshMotion};

fn cap() -> ForceControl3d {
    ForceControl3d { maximum_mach: 10.0 }
}
fn motion() -> MeshMotion {
    MeshMotion {
        lower_velocity: [-0.08, 0.0, 0.0],
        upper_velocity: [0.12, 0.0, 0.0],
    }
}
fn initial() -> Euler3d {
    field(
        gas(0.0),
        CartesianGrid::new(3, 2, 2, 0.3, 0.4, 0.5).unwrap(),
        Boundaries3d {
            x: WALL,
            ..PERIODIC
        },
        |x| {
            [
                1.0 + 0.03 * x[0],
                0.08 + 0.04 * x[0],
                0.1 - 0.03 * x[0],
                -0.06 + 0.02 * x[0],
                2.0 + 0.2 * x[0],
            ]
        },
    )
}
fn wall(temperature: f64) -> MolecularWall3d {
    MolecularWall3d::NoSlip {
        velocity: [0.0; 3],
        thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
    }
}
fn owner() -> MolecularEuler3d<ConstantTransport> {
    let flow = initial();
    let source = law(flow.gas(), 0.03, 0.01, 0.05);
    MolecularEuler3d::new(flow, source, || false)
        .unwrap()
        .with_wall(BoundarySide::XMin, wall(3.0))
        .unwrap()
        .with_wall(BoundarySide::XMax, wall(1.0))
        .unwrap()
}
fn check(report: MovingMolecularStepReport3d) {
    let r = report.transport;
    assert_eq!(r.first.elapsed, r.second.elapsed);
    assert_eq!(2.0 * r.first.elapsed.value(), r.elapsed.value());
    assert_eq!(r.advection.elapsed, r.elapsed);
    assert_eq!(r.source_change, r.advection.source_change);
    assert!(r.first.maximum_parabolic_number <= 0.40000000000001);
    assert!(r.second.maximum_parabolic_number <= 0.40000000000001);
    for x in r.conservation_residual.array() {
        close(x, 0.0, 5e-12);
    }
}

#[test]
fn g1_moving_forced_molecular_halves_use_actual_geometry_and_separate_source_work_heat() {
    let force = BodyForce3d::UniformForceDensity([0.7, -0.4, 0.2]);
    let f = [0.7, -0.4, 0.2];
    let mut coupled = owner();
    let mut work = coupled.flow().forcing_workspace(usize::MAX).unwrap();
    let report = coupled
        .step_moving_with_forcing(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion(),
            force,
            cap(),
            &mut work,
            || false,
        )
        .unwrap();
    check(report);
    let dt = report.transport.elapsed.value();
    assert_eq!(dt, 0.001);

    // Independently compose the original operators on their actual incoming
    // fields. A misplaced kick, stale geometry or cached wall law cannot match.
    let mut reference = initial();
    let v0 = reference.grid().cell_volume();
    let mut operator =
        MolecularTransport3d::new(&reference, law(reference.gas(), 0.03, 0.01, 0.05), || false)
            .unwrap()
            .with_wall(BoundarySide::XMin, wall(3.0))
            .unwrap()
            .with_wall(BoundarySide::XMax, wall(1.0))
            .unwrap();
    let (first, first_transfer) = operator
        .advance_at_geometry(
            &mut reference,
            Time::new(0.5 * dt),
            MolecularControl::default(),
            timed::Clock::split(0.0, dt, false),
            None,
            Some(motion()),
            false,
            &mut || false,
        )
        .unwrap();
    let first_work: f64 = reference
        .states()
        .iter()
        .map(|q| {
            (0..3)
                .map(|d| {
                    let j = 0.5 * dt * f[d];
                    j * (q.momentum[d] + 0.5 * j) / q.density * v0
                })
                .sum::<f64>()
        })
        .sum();
    let mut ref_work = reference.forcing_workspace(usize::MAX).unwrap();
    let advection = reference
        .step_moving_with_forcing(
            Time::new(dt),
            StepControl::default(),
            motion(),
            force,
            cap(),
            &mut ref_work,
            || false,
        )
        .unwrap();
    let v1 = reference.grid().cell_volume();
    assert_ne!(v0, v1);
    let second_work: f64 = reference
        .states()
        .iter()
        .map(|q| {
            (0..3)
                .map(|d| {
                    let j = 0.5 * dt * f[d];
                    j * (q.momentum[d] - 0.5 * j) / q.density * v1
                })
                .sum::<f64>()
        })
        .sum();
    let (second, second_transfer) = operator
        .advance_at_geometry(
            &mut reference,
            Time::new(0.5 * dt),
            MolecularControl::default(),
            timed::Clock::split(0.0, dt, true),
            None,
            Some(motion()),
            false,
            &mut || false,
        )
        .unwrap();
    assert_eq!(report.transport.first, first);
    assert_eq!(report.transport.advection, advection);
    assert_eq!(report.transport.second, second);
    assert_eq!(coupled.flow().states(), reference.states());
    close(
        report.transport.source_change.work,
        first_work + second_work,
        5e-13,
    );
    for d in 0..3 {
        close(
            report.transport.source_change.impulse[d],
            reference.states().len() as f64 * 0.5 * dt * f[d] * (v0 + v1),
            5e-13,
        );
    }
    assert!(report.transport.first.boundary_change.energy.abs() > 1e-7);
    assert!(report.transport.second.boundary_change.energy.abs() > 1e-7);
    for side in 0..6 {
        close(
            report.molecular_wall_work_j[side],
            first_transfer.work[side] + second_transfer.work[side],
            3e-14,
        );
        close(
            report.wall_heat_j[side],
            first_transfer.heat[side] + second_transfer.heat[side],
            3e-14,
        );
    }
    assert!(report.wall_heat_j[0] > 0.0 && report.wall_heat_j[1] < 0.0);
    assert!(
        report.molecular_wall_work_j[..2]
            .iter()
            .any(|q| q.abs() > 1e-8)
    );
}

#[test]
fn g0_moving_forced_molecular_zero_cases_preserve_native_fields_and_reports() {
    for stationary in [false, true] {
        let mut baseline = owner();
        let mut test = owner();
        let mut a = baseline.flow().forcing_workspace(usize::MAX).unwrap();
        let mut b = test.flow().forcing_workspace(usize::MAX).unwrap();
        let m = if stationary {
            MeshMotion::default()
        } else {
            motion()
        };
        let force = BodyForce3d::UniformAcceleration(if stationary {
            [0.3, -0.2, 0.1]
        } else {
            [0.0; 3]
        });
        let actual = test
            .step_moving_with_forcing(
                Time::new(0.001),
                MolecularStepControl::default(),
                m,
                force,
                cap(),
                &mut b,
                || false,
            )
            .unwrap();
        if stationary {
            let expected = baseline
                .step_with_forcing(
                    Time::new(0.001),
                    MolecularStepControl::default(),
                    force,
                    cap(),
                    &mut a,
                    || false,
                )
                .unwrap();
            assert_eq!(actual.transport, expected);
        } else {
            let expected = baseline
                .step_moving(Time::new(0.001), MolecularStepControl::default(), m, || {
                    false
                })
                .unwrap();
            assert_eq!(actual, expected);
        }
        assert_eq!(baseline.flow().states(), test.flow().states());
        assert_eq!(baseline.flow().grid(), test.flow().grid());
        assert_eq!(baseline.flow().time(), test.flow().time());
        assert_eq!(baseline.last_step(), test.last_step());
        for side in BoundarySide::ALL {
            assert_eq!(
                baseline.boundary_impulses(side),
                test.boundary_impulses(side)
            );
        }
    }
}

#[test]
fn g1_moving_forced_molecular_final_cooling_mach_refuses_and_retries_the_whole_interval() {
    let make = || {
        let flow = field(
            gas(0.0),
            CartesianGrid::new(1, 1, 1, 0.2, 0.2, 0.2).unwrap(),
            Boundaries3d {
                z: WALL,
                ..PERIODIC
            },
            |_| [1.0, 0.4, 0.0, 0.02, 2.0],
        );
        let source = law(flow.gas(), 0.0, 0.0, 0.5);
        MolecularEuler3d::new(flow, source, || false)
            .unwrap()
            .with_wall(BoundarySide::ZMin, wall(0.02))
            .unwrap()
            .with_wall(BoundarySide::ZMax, wall(0.02))
            .unwrap()
    };
    let m = MeshMotion {
        lower_velocity: [0.0, 0.0, 0.02],
        upper_velocity: [0.0, 0.0, 0.02],
    };
    let force = BodyForce3d::UniformAcceleration([0.02, 0.0, 0.0]);
    let mut probe = make();
    let mut work = probe.flow().forcing_workspace(usize::MAX).unwrap();
    let unrestricted = probe
        .step_moving_with_forcing(
            Time::new(0.005),
            MolecularStepControl::default(),
            m,
            force,
            cap(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(unrestricted.transport.elapsed, Time::new(0.005));
    let volume = probe.flow().grid().cell_volume();
    let central = State::from_array(
        unrestricted
            .transport
            .advection
            .after
            .array()
            .map(|x| x / volume),
    )
    .primitive(probe.flow().gas(), || false)
    .unwrap();
    let end = probe.flow().states()[0]
        .primitive(probe.flow().gas(), || false)
        .unwrap();
    let mach = |p: crate::compressible::spatial::Primitive| {
        p.velocity.into_iter().fold(0.0_f64, |v, u| v.hypot(u)) / p.sound_speed
    };
    assert!(mach(end) > mach(central) + 0.005);
    let cap = ForceControl3d {
        maximum_mach: 0.5 * (mach(end) + mach(central)),
    };
    let mut flow = make();
    let original = flow.flow().states().to_vec();
    let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
    let controls = MolecularStepControl {
        maximum_coupling_retries: 0,
        ..MolecularStepControl::default()
    };
    assert!(
        flow.step_moving_with_forcing(Time::new(0.005), controls, m, force, cap, &mut work, || {
            false
        },)
            .is_err()
    );
    assert_eq!(flow.flow().states(), original);
    assert_eq!(flow.flow().time(), Time::new(0.0));
    assert!(flow.last_step().is_none());
    let report = flow
        .step_moving_with_forcing(
            Time::new(0.005),
            MolecularStepControl::default(),
            m,
            force,
            cap,
            &mut work,
            || false,
        )
        .unwrap();
    check(report);
    assert!(report.transport.coupling_retries > 0);
    assert!(report.transport.elapsed.value() < 0.005);
    let mut replay = make();
    let mut replay_work = replay.flow().forcing_workspace(usize::MAX).unwrap();
    let direct = replay
        .step_moving_with_forcing(
            report.transport.elapsed,
            MolecularStepControl::default(),
            m,
            force,
            cap,
            &mut replay_work,
            || false,
        )
        .unwrap();
    assert_eq!(direct.transport.first, report.transport.first);
    assert_eq!(direct.transport.advection, report.transport.advection);
    assert_eq!(direct.transport.second, report.transport.second);
    assert_eq!(replay.flow().states(), flow.flow().states());
}

#[test]
fn g4_moving_forced_molecular_final_poll_preserves_accepted_binding_and_replays() {
    let force = BodyForce3d::UniformForceDensity([0.3, -0.2, 0.1]);
    let make = || {
        let mut flow = owner();
        let mut work = flow.flow().forcing_workspace(usize::MAX).unwrap();
        flow.step_moving_with_forcing(
            Time::new(0.001),
            MolecularStepControl::default(),
            motion(),
            force,
            cap(),
            &mut work,
            || false,
        )
        .unwrap();
        (flow, work)
    };
    let (mut probe, mut work) = make();
    let mut polls = 0;
    let expected = probe
        .step_moving_with_forcing(
            Time::new(0.002),
            MolecularStepControl::default(),
            motion(),
            force,
            cap(),
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [polls / 2, polls] {
        let (mut flow, mut work) = make();
        let original = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let previous = flow.last_step().copied();
        let time = flow.flow().time();
        let loads = BoundarySide::ALL.map(|s| flow.boundary_impulses(s).to_vec());
        let mut calls = 0;
        assert_eq!(
            flow.step_moving_with_forcing(
                Time::new(0.002),
                MolecularStepControl::default(),
                motion(),
                force,
                cap(),
                &mut work,
                || {
                    calls += 1;
                    calls == cutoff
                },
            ),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), original);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.last_step().copied(), previous);
        assert_eq!(flow.flow().time(), time);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(flow.boundary_impulses(side), loads[i]);
        }
        // The accepted molecular binding still describes the old box.
        flow.stress_field_moving(motion(), || false).unwrap();
        let accepted = flow
            .step_moving_with_forcing(
                Time::new(0.002),
                MolecularStepControl::default(),
                motion(),
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(accepted, expected);
        assert_eq!(flow.flow().states(), probe.flow().states());
    }
}
