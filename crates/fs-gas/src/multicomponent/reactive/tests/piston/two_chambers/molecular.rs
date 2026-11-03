use super::super::molecular::molecular;
use super::*;
use crate::dynamics::piston::PistonStep;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationOperator,
};
use crate::multicomponent::reactive::ReactiveSharedPistonGasReport;

fn molecular_step(
    left: &mut ReactingEuler1d,
    right: &mut ReactingEuler1d,
    left_operator: &mut DissipationOperator<ConstantTransport>,
    right_operator: &mut DissipationOperator<ConstantTransport>,
    piston: &mut Piston,
    duration: Time,
    cancelled: impl FnMut() -> bool,
) -> Result<PistonStep<ReactiveSharedPistonGasReport>, ReactiveError> {
    left.step_with_shared_piston_and_dissipation(
        right,
        left_operator,
        right_operator,
        piston,
        duration,
        ReactiveStepControl::default(),
        [DissipationControl::default(); 2],
        0.01,
        piston_controls(),
        cancelled,
    )
}

#[test]
fn g1_shared_piston_pairs_both_viscous_tractions_with_work_and_physical_energy() {
    let mut left = chamber(&mechanism(8.0, 30.0, 100.0), 0.0, 0.8, 1, chemistry());
    let mut right = chamber(&mechanism(2.0, 8.0, -30.0), 0.8, 1.2, 1, chemistry());
    let mut left_operator = molecular(left.flow());
    let mut right_operator = molecular(right.flow());
    let mut piston = Piston::new(0.2, 0.2).unwrap();
    let initial = [left.flow().totals(), right.flow().totals()];
    let initial_velocity = piston.velocity_m_s();
    let mut support = 0.0;
    let mut external_work = 0.0;
    let mut external_impulse = 0.0;
    let mut stress_impulse = 0.0;
    for step in 0..12 {
        let report = molecular_step(
            &mut left,
            &mut right,
            &mut left_operator,
            &mut right_operator,
            &mut piston,
            Time::new(0.001),
            || false,
        )
        .unwrap();
        let speed = report.midpoint_velocity_m_s;
        let dt = report.gas.left.elapsed.value();
        assert_eq!(report.gas.left.elapsed, report.gas.right.elapsed);
        assert_eq!(left.flow().time(), right.flow().time());
        assert_eq!(left_operator.mesh(), left.flow().mesh());
        assert_eq!(right_operator.mesh(), right.flow().mesh());
        assert_eq!(
            left.flow().mesh().positions()[1],
            right.flow().mesh().positions()[0],
        );
        let mut shared_stress = 0.0;
        for (side, flow) in [&report.gas.left, &report.gas.right]
            .into_iter()
            .enumerate()
        {
            let halves = flow.dissipation.unwrap();
            assert_eq!(flow.chemistry.cell_solves, 2);
            assert!(flow.chemistry.referenced_energy_change_j > 0.0);
            close(flow.mass_residual_kg, 0.0, 3e-12);
            close(flow.absolute_energy_residual_j, 0.0, 2e-8);
            for half in [halves.first, halves.second] {
                let ends = [half.left_boundary_change, half.right_boundary_change];
                let moving = 1 - side;
                assert_eq!(ends[side].energy, 0.0);
                close(ends[moving].energy, speed * ends[moving].momentum, 3e-14);
                assert_eq!(half.boundary_change.species_masses[..2], [0.0; 2]);
                shared_stress += ends[moving].momentum;
            }
            if step == 0 {
                // Independent SSP-RK2 solution of rho*du/dt =
                // 2*eta*(v_left+v_right-2*u)/L^2, initially rho=1,u=0.
                // The first D stage precedes chemistry, so neither reaction
                // heat nor the moved-volume A stage enters this oracle.
                let h = 0.5 * dt;
                let length = [0.8_f64, 1.2][side];
                let eta = (4.0 / 3.0) * 0.02;
                let stage_velocity = h * 2.0 * eta * speed / length.powi(2);
                let average_velocity = 0.5 * stage_velocity;
                let end_speeds = if side == 0 {
                    [0.0, speed]
                } else {
                    [speed, 0.0]
                };
                let expected = end_speeds
                    .map(|wall| h * AREA * (2.0 * eta / length) * (wall - average_velocity));
                assert_eq!(halves.first.accepted_steps, 1);
                close(
                    halves.first.left_boundary_change.momentum,
                    expected[0],
                    3e-14,
                );
                close(
                    halves.first.right_boundary_change.momentum,
                    expected[1],
                    3e-14,
                );
            }
        }
        let pressure_impulse = -report.gas.left.transport.right_boundary_change.momentum
            - report.gas.right.transport.left_boundary_change.momentum;
        close(
            report.gas_impulse_n_s,
            pressure_impulse - shared_stress,
            3e-14,
        );
        stress_impulse += shared_stress.abs();
        close(
            report.total_energy_residual_j,
            0.0,
            piston_controls().energy_tolerance_j,
        );
        close(
            report.combined_momentum_residual_n_s,
            0.0,
            piston_controls().momentum_tolerance_n_s,
        );
        support += report.fixed_wall_impulse_n_s;
        external_work += report.external_force_work_j;
        external_impulse += 0.01 * dt;
    }
    assert!(
        stress_impulse > 1e-7,
        "both molecular tractions must drive the root"
    );
    let kinetic_change =
        0.5 * piston.mass_kg() * (piston.velocity_m_s().powi(2) - initial_velocity.powi(2));
    close(
        physical_change(left.flow(), initial[0])
            + physical_change(right.flow(), initial[1])
            + kinetic_change
            - external_work,
        0.0,
        3e-7,
    );
    close(
        left.flow().totals().momentum - initial[0].momentum + right.flow().totals().momentum
            - initial[1].momentum
            + piston.mass_kg() * (piston.velocity_m_s() - initial_velocity)
            + support
            - external_impulse,
        0.0,
        3e-10,
    );
    assert_eq!(left.flow().mesh().positions()[0], 0.0);
    assert_eq!(right.flow().mesh().positions()[1], 2.0);
}

#[test]
fn g1_shortened_shared_molecular_step_replays_both_original_stage_sequences() {
    let make = || {
        (
            chamber(&mechanism(40.0, 500.0, 0.0), 0.0, 1.0, 4, chemistry()),
            chamber(&mechanism(0.0, 2.0, 30.0), 1.0, 1.0, 4, chemistry()),
        )
    };
    let (mut left, mut right) = make();
    let mut left_operator = molecular(left.flow());
    let mut right_operator = molecular(right.flow());
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let report = molecular_step(
        &mut left,
        &mut right,
        &mut left_operator,
        &mut right_operator,
        &mut piston,
        Time::new(0.2),
        || false,
    )
    .unwrap();
    assert!(report.gas.left.coupling_retries > 0);
    assert_eq!(
        report.gas.left.coupling_retries,
        report.gas.right.coupling_retries
    );
    let dt = report.gas.left.elapsed.value();
    let speed = report.midpoint_velocity_m_s;
    let (mut replay_left, mut replay_right) = make();
    for (side, (actual, replay, accepted_report)) in [
        (&left, &mut replay_left, &report.gas.left),
        (&right, &mut replay_right, &report.gas.right),
    ]
    .into_iter()
    .enumerate()
    {
        let mut operator = molecular(replay.flow());
        operator
            .admit_motion(&replay.accepted, &mut || false)
            .unwrap();
        let velocities = if side == 0 {
            [0.0, speed]
        } else {
            [speed, 0.0]
        };
        let original_volumes = replay.flow().mesh().volumes().to_vec();
        let first = operator
            .advance_moving_trial(
                &mut replay.accepted,
                Time::new(0.5 * dt),
                DissipationControl::default(),
                velocities,
                &mut || false,
            )
            .unwrap();
        let mut chemistry = ChemistryReport::default();
        react_cells(
            &mut replay.accepted,
            &mut replay.reactor,
            &mut replay.amounts,
            0.5 * dt,
            &mut chemistry,
            &mut || false,
        )
        .unwrap();
        let transport = replay
            .accepted
            .step_with_moving_walls(Time::new(dt), StepControl::default(), velocities, || false)
            .unwrap();
        assert_eq!(transport.elapsed.value(), dt);
        assert_ne!(replay.flow().mesh().volumes(), original_volumes);
        react_cells(
            &mut replay.accepted,
            &mut replay.reactor,
            &mut replay.amounts,
            0.5 * dt,
            &mut chemistry,
            &mut || false,
        )
        .unwrap();
        let second = operator
            .advance_moving_trial(
                &mut replay.accepted,
                Time::new(0.5 * dt),
                DissipationControl::default(),
                velocities,
                &mut || false,
            )
            .unwrap();
        assert_eq!(actual.flow().mesh(), replay.flow().mesh());
        assert_eq!(actual.flow().time(), replay.flow().time());
        for (&a, &b) in actual.accepted.states.iter().zip(&replay.accepted.states) {
            close(a, b, 3e-12);
        }
        assert_eq!(accepted_report.chemistry, chemistry);
        assert_eq!(
            accepted_report.dissipation,
            Some(SplitDissipationReport { first, second })
        );
    }
}

#[test]
fn g4_shared_molecular_final_cancellation_preserves_both_bindings_and_exact_retry() {
    let (mut left, mut right) = pair();
    let (mut reference_left, mut reference_right) = pair();
    let mut left_operator = molecular(left.flow());
    let mut right_operator = molecular(right.flow());
    let mut reference_left_operator = molecular(reference_left.flow());
    let mut reference_right_operator = molecular(reference_right.flow());
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let mut reference_piston = piston;
    let dt = Time::new(0.002);
    molecular_step(
        &mut left,
        &mut right,
        &mut left_operator,
        &mut right_operator,
        &mut piston,
        dt,
        || false,
    )
    .unwrap();
    molecular_step(
        &mut reference_left,
        &mut reference_right,
        &mut reference_left_operator,
        &mut reference_right_operator,
        &mut reference_piston,
        dt,
        || false,
    )
    .unwrap();
    let states = [left.accepted.states.clone(), right.accepted.states.clone()];
    let meshes = [left.flow().mesh().clone(), right.flow().mesh().clone()];
    let time = left.flow().time();
    let original_piston = piston;
    let mut polls = 0;
    let expected = molecular_step(
        &mut reference_left,
        &mut reference_right,
        &mut reference_left_operator,
        &mut reference_right_operator,
        &mut reference_piston,
        dt,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    let error = molecular_step(
        &mut left,
        &mut right,
        &mut left_operator,
        &mut right_operator,
        &mut piston,
        dt,
        || {
            seen += 1;
            seen == polls
        },
    )
    .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(seen, polls);
    for (side, flow) in [&left, &right].into_iter().enumerate() {
        assert_eq!(flow.accepted.states, states[side]);
        assert_eq!(flow.flow().mesh(), &meshes[side]);
        assert_eq!(flow.flow().time(), time);
    }
    assert_eq!(left_operator.mesh(), &meshes[0]);
    assert_eq!(right_operator.mesh(), &meshes[1]);
    assert_eq!(piston, original_piston);
    let recovered = molecular_step(
        &mut left,
        &mut right,
        &mut left_operator,
        &mut right_operator,
        &mut piston,
        dt,
        || false,
    )
    .unwrap();
    assert_eq!(recovered, expected);
    assert_eq!(left.accepted.states, reference_left.accepted.states);
    assert_eq!(right.accepted.states, reference_right.accepted.states);
    assert_eq!(left_operator.mesh(), reference_left_operator.mesh());
    assert_eq!(right_operator.mesh(), reference_right_operator.mesh());
    assert_eq!(piston, reference_piston);
}

#[test]
fn g4_right_source_refusal_keeps_both_molecular_chambers_and_piston_accepted() {
    let source = mechanism(0.0, 2.0, 0.0);
    let mut left = chamber(&source, 0.0, 1.0, 8, chemistry());
    let mut right = chamber(
        &source,
        1.0,
        1.0,
        8,
        IsochoricKineticsSolve {
            maximum_temperature_k: 1.0000001,
            ..chemistry()
        },
    );
    let mut left_operator = molecular(left.flow());
    let mut right_operator = molecular(right.flow());
    let states = [left.accepted.states.clone(), right.accepted.states.clone()];
    let meshes = [left.flow().mesh().clone(), right.flow().mesh().clone()];
    let mut piston = Piston::new(0.1, 0.0).unwrap();
    let original_piston = piston;
    let result = left.step_with_shared_piston_and_dissipation(
        &mut right,
        &mut left_operator,
        &mut right_operator,
        &mut piston,
        Time::new(0.001),
        ReactiveStepControl {
            maximum_coupling_retries: 0,
            ..ReactiveStepControl::default()
        },
        [DissipationControl::default(); 2],
        1.0,
        piston_controls(),
        || false,
    );
    assert!(matches!(result, Err(ReactiveError::RetryLimit { .. })));
    for (side, flow) in [&left, &right].into_iter().enumerate() {
        assert_eq!(flow.accepted.states, states[side]);
        assert_eq!(flow.flow().mesh(), &meshes[side]);
        assert_eq!(flow.flow().time(), Time::new(0.0));
    }
    assert_eq!(left_operator.mesh(), &meshes[0]);
    assert_eq!(right_operator.mesh(), &meshes[1]);
    assert_eq!(piston, original_piston);
}
