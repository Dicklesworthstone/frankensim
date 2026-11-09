use super::*;

fn work(flow: &ReactingEuler3d) -> ForcingWorkspace3d {
    flow.flow()
        .forcing_workspace(flow.flow().required_forcing_workspace_bytes().unwrap())
        .unwrap()
}

fn force_control() -> ForceControl3d {
    ForceControl3d { maximum_mach: 4.0 }
}

#[allow(clippy::too_many_arguments)]
fn forced(
    flow: &mut ReactingEuler3d,
    op: &mut MolecularTransport3d<ConstantTransport>,
    work: &mut ForcingWorkspace3d,
    molecular: bool,
    duration: f64,
    motion: MeshMotion,
    force: BodyForce3d<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<MovingReactiveStepReport3d, ReactiveError> {
    if molecular {
        let bytes = flow.workspace_bytes() + op.workspace_bytes() + work.workspace_bytes();
        flow.step_moving_with_dissipation_and_forcing(
            op,
            Time::new(duration),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            motion,
            force,
            force_control(),
            work,
            cancelled,
        )
    } else {
        flow.step_moving_with_forcing(
            Time::new(duration),
            ReactiveStepControl::default(),
            motion,
            force,
            force_control(),
            work,
            cancelled,
        )
    }
}

#[test]
fn g1_accelerating_translating_reactor_retains_chemical_state_and_exact_kinetic_work() {
    let mut flow = uniform();
    let mut reference = uniform();
    let mut scratch = work(&flow);
    let bytes = flow.workspace_bytes() + scratch.workspace_bytes();
    let initial = flow
        .flow()
        .cell(0)
        .unwrap()
        .primitive(flow.flow().gas(), || false)
        .unwrap();
    let acceleration = [0.4, -0.3, 0.2];
    let motion = MeshMotion {
        lower_velocity: [0.2, -0.1, 0.05],
        upper_velocity: [0.2, -0.1, 0.05],
    };
    let report = flow
        .step_moving_with_forcing(
            Time::new(0.0005),
            ReactiveStepControl::default(),
            motion,
            BodyForce3d::UniformAcceleration(acceleration),
            force_control(),
            &mut scratch,
            || false,
        )
        .unwrap();
    let unforced = reference
        .step_moving(
            report.reactive.elapsed,
            ReactiveStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
    balance(&report.reactive);
    let dt = report.reactive.elapsed.value();
    let mass = report.reactive.before.mass();
    let expected_work: f64 = (0..3)
        .map(|axis| {
            mass * acceleration[axis]
                * dt
                * (initial.velocity[axis] + 0.5 * acceleration[axis] * dt)
        })
        .sum();
    close(report.reactive.source_change.work, expected_work, 3e-12);
    for axis in 0..3 {
        close(
            report.reactive.source_change.impulse[axis],
            mass * acceleration[axis] * dt,
            3e-12,
        );
        close(
            report.reactive.after.momentum[axis] - report.reactive.before.momentum[axis],
            mass * acceleration[axis] * dt,
            3e-12,
        );
    }
    for index in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(index).unwrap();
        let r = reference.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        let pr = r.primitive(reference.flow().gas(), || false).unwrap();
        for k in 0..2 {
            close(q.partial_densities[k], r.partial_densities[k], 3e-11);
        }
        close(p.temperature, pr.temperature, 3e-11);
        close(p.pressure, pr.pressure, 3e-11);
        for axis in 0..3 {
            close(
                p.velocity[axis],
                initial.velocity[axis] + acceleration[axis] * dt,
                3e-12,
            );
        }
    }
    assert!(report.reactive.chemistry.species_change_kg[1] > 0.0);
    close(
        report.reactive.chemistry.referenced_energy_change_j,
        unforced.reactive.chemistry.referenced_energy_change_j,
        3e-11,
    );
    assert_eq!(flow.workspace_bytes() + scratch.workspace_bytes(), bytes);
}

#[test]
fn g1_exact_coriolis_rotates_reacting_momentum_without_inventing_heat_on_a_moving_grid() {
    let mut flow = uniform();
    let mut reference = uniform();
    let mut scratch = work(&flow);
    let initial = flow
        .flow()
        .cell(0)
        .unwrap()
        .primitive(flow.flow().gas(), || false)
        .unwrap();
    let omega = [[0.0, 0.0, 17.0]];
    let acceleration = [[0.0; 3]];
    let motion = MeshMotion {
        lower_velocity: [0.2, -0.1, 0.05],
        upper_velocity: [0.2, -0.1, 0.05],
    };
    let report = flow
        .step_moving_with_forcing(
            Time::new(0.0005),
            ReactiveStepControl::default(),
            motion,
            BodyForce3d::CoriolisExact {
                angular_velocity: &omega,
                acceleration: &acceleration,
            },
            force_control(),
            &mut scratch,
            || false,
        )
        .unwrap();
    reference
        .step_moving(
            report.reactive.elapsed,
            ReactiveStepControl::default(),
            motion,
            || false,
        )
        .unwrap();
    let angle = -2.0 * omega[0][2] * report.reactive.elapsed.value();
    let (sine, cosine) = angle.sin_cos();
    let velocity = [
        initial.velocity[0] * cosine - initial.velocity[1] * sine,
        initial.velocity[0] * sine + initial.velocity[1] * cosine,
        initial.velocity[2],
    ];
    assert_eq!(report.reactive.source_change.work, 0.0);
    balance(&report.reactive);
    for axis in 0..3 {
        close(
            report.reactive.source_change.impulse[axis],
            report.reactive.before.mass() * (velocity[axis] - initial.velocity[axis]),
            3e-12,
        );
    }
    for index in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(index).unwrap();
        let qr = reference.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        let pr = qr.primitive(reference.flow().gas(), || false).unwrap();
        for axis in 0..3 {
            close(p.velocity[axis], velocity[axis], 3e-12);
        }
        for k in 0..2 {
            close(q.partial_densities[k], qr.partial_densities[k], 3e-11);
        }
        close(p.temperature, pr.temperature, 3e-11);
        close(p.pressure, pr.pressure, 3e-11);
    }
}

#[test]
fn g2_reacting_molecular_force_density_integrates_actual_old_and_new_volume() {
    let force_density = [0.8, -0.4, 0.3];
    for axis in 0..3 {
        let (mut flow, mut op) = case(0.0);
        let mut scratch = work(&flow);
        let volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
        let bytes = flow.workspace_bytes() + op.workspace_bytes() + scratch.workspace_bytes();
        let mut motion = MeshMotion::default();
        motion.lower_velocity[axis] = -0.04;
        motion.upper_velocity[axis] = 0.08;
        let report = forced(
            &mut flow,
            &mut op,
            &mut scratch,
            true,
            0.0005,
            motion,
            BodyForce3d::UniformForceDensity(force_density),
            || false,
        )
        .unwrap();
        moving_balance(&flow, report, motion);
        let next_volume = flow.flow().grid().cell_volume() * flow.flow().cell_count() as f64;
        for a in 0..3 {
            close(
                report.reactive.source_change.impulse[a],
                force_density[a] * 0.5 * (volume + next_volume) * report.reactive.elapsed.value(),
                3e-13,
            );
        }
        assert!(report.reactive.source_change.work.abs() > 1e-8);
        assert!(report.reactive.chemistry.species_change_kg[1] > 0.0);
        assert!(report.wall_heat_j.iter().any(|v| v.abs() > 1e-8));
        assert_eq!(
            flow.workspace_bytes() + op.workspace_bytes() + scratch.workspace_bytes(),
            bytes
        );
    }
}

#[test]
fn g0_zero_force_retains_both_original_moving_reacting_paths_and_complete_budget() {
    let motion = MeshMotion {
        lower_velocity: [-0.04, 0.0, 0.0],
        upper_velocity: [0.08, 0.0, 0.0],
    };
    for molecular in [false, true] {
        let (mut a, mut ao) = case(0.0);
        let (mut b, mut bo) = case(0.0);
        let mut scratch = work(&b);
        let expected = step(&mut a, &mut ao, molecular, 0.0005, motion, || false).unwrap();
        let actual = forced(
            &mut b,
            &mut bo,
            &mut scratch,
            molecular,
            0.0005,
            motion,
            BodyForce3d::UniformAcceleration([0.0; 3]),
            || false,
        )
        .unwrap();
        assert_eq!(expected, actual);
        assert_eq!(a.flow().states(), b.flow().states());
        assert_eq!(loads(&a), loads(&b));
        assert_eq!(molecular_loads(&ao), molecular_loads(&bo));
    }
    let (mut flow, mut op) = case(0.0);
    let mut scratch = work(&flow);
    let state = flow.flow().states().to_vec();
    let grid = *flow.flow().grid();
    let bytes = flow.workspace_bytes() + op.workspace_bytes() + scratch.workspace_bytes();
    assert!(
        flow.step_moving_with_dissipation_and_forcing(
            &mut op,
            Time::new(0.0005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes - 1,
            motion,
            BodyForce3d::UniformAcceleration([0.4, -0.3, 0.2]),
            force_control(),
            &mut scratch,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.flow().states(), state);
    assert_eq!(*flow.flow().grid(), grid);
    assert_eq!(flow.last_step(), None);
    assert_eq!(op.last_step(), None);
}

#[test]
fn g4_forced_moving_reacting_last_checkpoint_and_mach_refusal_preserve_accepted_state() {
    let motion = MeshMotion {
        lower_velocity: [-0.04, 0.0, 0.0],
        upper_velocity: [0.08, 0.0, 0.0],
    };
    let force = BodyForce3d::UniformAcceleration([0.4, -0.3, 0.2]);
    let make = || {
        let (mut flow, mut op) = case(0.0);
        let mut scratch = work(&flow);
        forced(
            &mut flow,
            &mut op,
            &mut scratch,
            true,
            0.0002,
            motion,
            force,
            || false,
        )
        .unwrap();
        (flow, op, scratch)
    };
    let (mut expected, mut expected_op, mut expected_work) = make();
    let mut calls = 0;
    let reference = forced(
        &mut expected,
        &mut expected_op,
        &mut expected_work,
        true,
        0.0005,
        motion,
        force,
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    for cutoff in [calls / 2, calls] {
        let (mut flow, mut op, mut scratch) = make();
        let state = flow.flow().states().to_vec();
        let grid = *flow.flow().grid();
        let time = flow.flow().time();
        let report = flow.last_step().copied();
        let op_report = op.last_step().copied();
        let boundary = loads(&flow);
        let op_boundary = molecular_loads(&op);
        let mut seen = 0;
        assert_eq!(
            forced(
                &mut flow,
                &mut op,
                &mut scratch,
                true,
                0.0005,
                motion,
                force,
                || {
                    seen += 1;
                    seen == cutoff
                }
            ),
            Err(ReactiveError::Cancelled)
        );
        assert_eq!(flow.flow().states(), state);
        assert_eq!(*flow.flow().grid(), grid);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), report);
        assert_eq!(op.last_step().copied(), op_report);
        assert_eq!(loads(&flow), boundary);
        assert_eq!(molecular_loads(&op), op_boundary);
        let recovered = forced(
            &mut flow,
            &mut op,
            &mut scratch,
            true,
            0.0005,
            motion,
            force,
            || false,
        )
        .unwrap();
        assert_eq!(recovered, reference);
        assert_eq!(flow.flow().states(), expected.flow().states());
    }
    let (mut flow, mut op, mut scratch) = make();
    let state = flow.flow().states().to_vec();
    let before = flow.last_step().copied();
    let op_before = op.last_step().copied();
    let bytes = flow.workspace_bytes() + op.workspace_bytes() + scratch.workspace_bytes();
    assert!(
        flow.step_moving_with_dissipation_and_forcing(
            &mut op,
            Time::new(0.0005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            motion,
            force,
            ForceControl3d {
                maximum_mach: 0.001
            },
            &mut scratch,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.flow().states(), state);
    assert_eq!(flow.last_step().copied(), before);
    assert_eq!(op.last_step().copied(), op_before);
}
