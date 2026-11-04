use super::*;
use crate::dynamics::piston::{PistonSuspension, SuspendedPistonStep};

fn small_pair() -> (ReactingEuler1d, ReactingEuler1d) {
    (
        chamber(&mechanism(8.0, 30.0, 100.0), 0.0, 0.8, 2, chemistry()),
        chamber(&mechanism(2.0, 8.0, -30.0), 0.8, 1.2, 2, chemistry()),
    )
}

#[allow(clippy::too_many_arguments)]
fn suspended_step(
    left: &mut ReactingEuler1d,
    right: &mut ReactingEuler1d,
    left_operator: &mut DissipationOperator<ConstantTransport>,
    right_operator: &mut DissipationOperator<ConstantTransport>,
    piston: &mut Piston,
    suspension: PistonSuspension,
    cancelled: impl FnMut() -> bool,
) -> Result<SuspendedPistonStep<ReactiveSharedPistonGasReport>, ReactiveError> {
    left.step_with_shared_suspended_piston_and_dissipation(
        right,
        left_operator,
        right_operator,
        piston,
        Time::new(0.001),
        ReactiveStepControl::default(),
        [DissipationControl::default(); 2],
        0.01,
        suspension,
        piston_controls(),
        cancelled,
    )
}

#[test]
fn g1_shared_suspension_separates_viscous_gas_heating_chemistry_and_damper_loss() {
    let (mut left, mut right) = small_pair();
    let mut left_operator = molecular(left.flow());
    let mut right_operator = molecular(right.flow());
    let mut piston = Piston::new(0.2, 0.2).unwrap();
    let suspension = PistonSuspension::new(1.5, 0.3, 0.75).unwrap();
    let before = [left.flow().totals(), right.flow().totals()];
    let mut loss = 0.0;
    let mut work = 0.0;
    let mut impulse = 0.0;
    let mut support = 0.0;
    let mut viscous_impulse = 0.0;
    for iteration in 0..12 {
        let x0 = left.flow().mesh().positions()[2];
        let v0 = piston.velocity_m_s();
        let report = suspended_step(
            &mut left,
            &mut right,
            &mut left_operator,
            &mut right_operator,
            &mut piston,
            suspension,
            || false,
        )
        .unwrap();
        let mechanical = &report.piston;
        let dt = mechanical.gas.left.elapsed.value();
        let w = mechanical.midpoint_velocity_m_s;
        let x1 = left.flow().mesh().positions()[2];
        let expected =
            (2.0 * 0.2 * v0 + mechanical.gas_impulse_n_s + dt * (0.01 - 1.5 * (x0 - 0.75)))
                / (2.0 * 0.2 + 0.3 * dt + 0.5 * 1.5 * dt * dt);
        close(w, expected, 2e-10);
        let mut shared_stress = 0.0;
        for (side, flow) in [&mechanical.gas.left, &mechanical.gas.right]
            .into_iter()
            .enumerate()
        {
            let halves = flow.dissipation.unwrap();
            assert!(flow.chemistry.referenced_energy_change_j > 0.0);
            assert!(flow.chemistry.species_change_kg[1] > 0.0);
            close(flow.mass_residual_kg, 0.0, 3e-12);
            for half in [halves.first, halves.second] {
                let ends = [half.left_boundary_change, half.right_boundary_change];
                assert_eq!(ends[side].energy, 0.0);
                close(ends[1 - side].energy, w * ends[1 - side].momentum, 3e-14);
                assert_eq!(half.boundary_change.species_masses[..2], [0.0; 2]);
                shared_stress += ends[1 - side].momentum;
            }
            if iteration == 0 {
                // Independent two-cell SSP-RK2 shear solution, before chemistry.
                // The moving end accelerates its cell in stage one; the inner
                // face transfers momentum to the other cell in stage two.
                let h = 0.5 * dt;
                let length = [0.8_f64, 1.2][side];
                let dx = 0.5 * length;
                let eta = (4.0 / 3.0) * 0.02;
                let u_stage = h * 2.0 * eta * w / dx.powi(2);
                let u_moving = 0.5 * (u_stage + h * eta * (2.0 * w - 3.0 * u_stage) / dx.powi(2));
                let u_fixed = 0.5 * h * eta * u_stage / dx.powi(2);
                let moving_impulse = h * AREA * 2.0 * eta / dx * (w - 0.5 * u_stage);
                let expected_heat =
                    w * moving_impulse - 0.5 * AREA * dx * (u_moving.powi(2) + u_fixed.powi(2));
                let first = halves.first;
                let (mut replay_left, mut replay_right) = small_pair();
                let replay = if side == 0 {
                    &mut replay_left
                } else {
                    &mut replay_right
                };
                let mut operator = molecular(replay.flow());
                operator
                    .admit_motion(&replay.accepted, &mut || false)
                    .unwrap();
                let replayed = operator
                    .advance_moving_trial(
                        &mut replay.accepted,
                        Time::new(h),
                        DissipationControl::default(),
                        if side == 0 { [0.0, w] } else { [w, 0.0] },
                        &mut || false,
                    )
                    .unwrap();
                let mut kinetic = 0.0;
                for i in 0..2 {
                    let p = replay
                        .flow()
                        .gas()
                        .primitive(replay.flow().cell(i).unwrap(), || false)
                        .unwrap();
                    close(
                        p.velocity,
                        if i == 1 - side { u_moving } else { u_fixed },
                        3e-14,
                    );
                    kinetic +=
                        0.5 * p.density * p.velocity.powi(2) * replay.flow().mesh().volumes()[i];
                }
                let actual_heat = replayed.after.energy - replayed.before.energy - kinetic;
                assert_eq!(first.accepted_steps, 1);
                assert!(expected_heat > 0.0);
                close(
                    first.after.energy - first.before.energy,
                    w * moving_impulse,
                    3e-14,
                );
                close(
                    first.after.momentum - first.before.momentum,
                    moving_impulse,
                    3e-14,
                );
                close(actual_heat, expected_heat, 3e-14);
            }
        }
        let pressure = -mechanical.gas.left.transport.right_boundary_change.momentum
            - mechanical.gas.right.transport.left_boundary_change.momentum;
        close(mechanical.gas_impulse_n_s, pressure - shared_stress, 3e-14);
        close(report.damping_energy_j, 0.3 * w * w * dt, 3e-16);
        assert!(report.total_energy_residual_j.abs() <= piston_controls().energy_tolerance_j);
        assert_eq!(x1, right.flow().mesh().positions()[0]);
        assert_eq!(left_operator.mesh(), left.flow().mesh());
        assert_eq!(right_operator.mesh(), right.flow().mesh());
        loss += report.damping_energy_j;
        work += 0.01 * (x1 - x0);
        impulse += 0.01 * dt + report.suspension_impulse_n_s;
        support += mechanical.fixed_wall_impulse_n_s;
        viscous_impulse += shared_stress.abs();
    }
    let x = left.flow().mesh().positions()[2];
    let kinetic = 0.5 * 0.2 * (piston.velocity_m_s().powi(2) - 0.2_f64.powi(2));
    let stored = 0.5 * 1.5 * ((x - 0.75).powi(2) - 0.05_f64.powi(2));
    close(
        physical_change(left.flow(), before[0])
            + physical_change(right.flow(), before[1])
            + kinetic
            + stored
            + loss,
        work,
        3e-7,
    );
    close(
        left.flow().totals().momentum - before[0].momentum + right.flow().totals().momentum
            - before[1].momentum
            + 0.2 * (piston.velocity_m_s() - 0.2)
            + support,
        impulse,
        3e-10,
    );
    assert!(viscous_impulse > 1e-7);
    assert!(loss > 1e-5);
}

#[test]
fn g4_suspended_molecular_zero_parity_final_cancel_and_both_bindings_replay() {
    let (mut left, mut right) = small_pair();
    let (mut reference_left, mut reference_right) = small_pair();
    let mut left_operator = molecular(left.flow());
    let mut right_operator = molecular(right.flow());
    let mut reference_left_operator = molecular(reference_left.flow());
    let mut reference_right_operator = molecular(reference_right.flow());
    let mut piston = Piston::new(0.2, 0.2).unwrap();
    let mut reference_piston = piston;
    let original = molecular_step(
        &mut reference_left,
        &mut reference_right,
        &mut reference_left_operator,
        &mut reference_right_operator,
        &mut reference_piston,
        Time::new(0.001),
        || false,
    )
    .unwrap();
    let zero = suspended_step(
        &mut left,
        &mut right,
        &mut left_operator,
        &mut right_operator,
        &mut piston,
        PistonSuspension::new(0.0, 0.0, 10.0).unwrap(),
        || false,
    )
    .unwrap();
    assert_eq!(zero.piston, original);
    assert_eq!(left.accepted.states, reference_left.accepted.states);
    assert_eq!(right.accepted.states, reference_right.accepted.states);
    let states = [left.accepted.states.clone(), right.accepted.states.clone()];
    let meshes = [left.flow().mesh().clone(), right.flow().mesh().clone()];
    let time = left.flow().time();
    let before_piston = piston;
    let suspension = PistonSuspension::new(1.5, 0.3, 0.75).unwrap();
    let mut polls = 0;
    let expected = suspended_step(
        &mut reference_left,
        &mut reference_right,
        &mut reference_left_operator,
        &mut reference_right_operator,
        &mut reference_piston,
        suspension,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        suspended_step(
            &mut left,
            &mut right,
            &mut left_operator,
            &mut right_operator,
            &mut piston,
            suspension,
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(seen, polls);
    for (j, flow) in [&left, &right].into_iter().enumerate() {
        assert_eq!(flow.accepted.states, states[j]);
        assert_eq!(flow.flow().mesh(), &meshes[j]);
        assert_eq!(flow.flow().time(), time);
    }
    assert_eq!(left_operator.mesh(), &meshes[0]);
    assert_eq!(right_operator.mesh(), &meshes[1]);
    assert_eq!(piston, before_piston);
    let recovered = suspended_step(
        &mut left,
        &mut right,
        &mut left_operator,
        &mut right_operator,
        &mut piston,
        suspension,
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
