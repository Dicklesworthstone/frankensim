use super::*;
use crate::dynamics::piston::{PistonSuspension, SuspendedPistonStep};
use crate::multicomponent::SharedPistonGasReport;

fn suspended_step(
    left: &mut ReactingEuler1d,
    right: &mut ReactingEuler1d,
    piston: &mut Piston,
    suspension: PistonSuspension,
    cancelled: impl FnMut() -> bool,
) -> Result<ReactiveSharedPistonStep<SuspendedPistonStep<SharedPistonGasReport>>, ReactiveError> {
    left.step_with_shared_suspended_piston(
        right,
        piston,
        Time::new(0.002),
        ReactiveStepControl::default(),
        -0.001,
        suspension,
        piston_controls(),
        cancelled,
    )
}

#[test]
fn g1_both_chemical_sources_drive_spring_storage_and_independent_damper_loss() {
    let (mut left, mut right) = pair();
    let initial = [left.flow().totals(), right.flow().totals()];
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let suspension = PistonSuspension::new(1.6, 0.2, 0.8).unwrap();
    let mut loss = 0.0;
    let mut work = 0.0;
    let mut impulse = 0.0;
    let mut support = 0.0;
    let mut species = [[0.0; 2]; 2];
    for _ in 0..30 {
        let before = [left.flow().totals(), right.flow().totals()];
        let x0 = left.flow().mesh().positions()[left.flow().cell_count()];
        let v0 = piston.velocity_m_s();
        let report =
            suspended_step(&mut left, &mut right, &mut piston, suspension, || false).unwrap();
        let spring = &report.piston;
        let mechanical = &spring.piston;
        let dt = report.left.elapsed.value();
        let x1 = left.flow().mesh().positions()[left.flow().cell_count()];
        let w = 0.5 * (v0 + piston.velocity_m_s());
        assert_eq!(report.left.elapsed, report.right.elapsed);
        assert_eq!(report.left.elapsed, mechanical.gas.left.elapsed);
        assert_eq!(x1, right.flow().mesh().positions()[0]);
        for (j, (flow, side)) in [(&left, &report.left), (&right, &report.right)]
            .into_iter()
            .enumerate()
        {
            close(
                report.gas_energy_change_j[j],
                physical_change(flow.flow(), before[j]),
                3e-12,
            );
            assert!(side.chemistry.referenced_energy_change_j > 0.0);
            close(side.mass_residual_kg, 0.0, 3e-12);
            for k in 0..2 {
                assert_eq!(side.transport.boundary_change.species_masses[k], 0.0);
                species[j][k] += side.chemistry.species_change_kg[k];
            }
        }
        close(
            spring.spring_energy_change_j,
            0.5 * 1.6 * ((x1 - 0.8).powi(2) - (x0 - 0.8).powi(2)),
            2e-16,
        );
        close(spring.damping_energy_j, 0.2 * w * w * dt, 2e-16);
        assert!(report.total_energy_residual_j.abs() <= piston_controls().energy_tolerance_j);
        assert!(
            report.combined_momentum_residual_n_s.abs() <= piston_controls().momentum_tolerance_n_s
        );
        loss += spring.damping_energy_j;
        work += -0.001 * (x1 - x0);
        impulse += -0.001 * dt + spring.suspension_impulse_n_s;
        support += mechanical.fixed_wall_impulse_n_s;
    }
    let x = left.flow().mesh().positions()[left.flow().cell_count()];
    let kinetic = 0.5 * 0.2 * piston.velocity_m_s().powi(2);
    close(
        physical_change(left.flow(), initial[0])
            + physical_change(right.flow(), initial[1])
            + kinetic
            + 0.5 * 1.6 * (x - 0.8).powi(2)
            + loss,
        work,
        3e-7,
    );
    close(
        left.flow().totals().momentum - initial[0].momentum + right.flow().totals().momentum
            - initial[1].momentum
            + 0.2 * piston.velocity_m_s()
            + support,
        impulse,
        3e-10,
    );
    for (j, flow) in [&left, &right].into_iter().enumerate() {
        let after = flow.flow().totals();
        for k in 0..2 {
            close(
                after.species_masses[k] - initial[j].species_masses[k],
                species[j][k],
                3e-12,
            );
        }
        assert!(after.species_masses[1] > 0.0);
    }
    assert!(piston.velocity_m_s() > 0.001);
    assert!(x > 0.8);
    assert!(loss > 0.0);
}

#[test]
fn g4_zero_law_parity_and_final_cancel_keep_both_reacting_owners_reusable() {
    let (mut left, mut right) = pair();
    let (mut reference_left, mut reference_right) = pair();
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let mut reference_piston = piston;
    let free = PistonSuspension::new(0.0, 0.0, 10.0).unwrap();
    let original = reference_left
        .step_with_shared_piston(
            &mut reference_right,
            &mut reference_piston,
            Time::new(0.002),
            ReactiveStepControl::default(),
            -0.001,
            piston_controls(),
            || false,
        )
        .unwrap();
    let zero = suspended_step(&mut left, &mut right, &mut piston, free, || false).unwrap();
    assert_eq!(zero.left, original.left);
    assert_eq!(zero.right, original.right);
    assert_eq!(zero.piston.piston, original.piston);
    assert_eq!(
        zero.total_energy_residual_j,
        original.total_energy_residual_j
    );
    assert_eq!(
        zero.combined_momentum_residual_n_s,
        original.combined_momentum_residual_n_s
    );
    assert_eq!(left.accepted.states, reference_left.accepted.states);
    assert_eq!(right.accepted.states, reference_right.accepted.states);
    let states = [left.accepted.states.clone(), right.accepted.states.clone()];
    let meshes = [left.flow().mesh().clone(), right.flow().mesh().clone()];
    let time = left.flow().time();
    let before_piston = piston;
    let suspension = PistonSuspension::new(1.6, 0.2, 0.8).unwrap();
    let mut polls = 0;
    let expected = suspended_step(
        &mut reference_left,
        &mut reference_right,
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
        suspended_step(&mut left, &mut right, &mut piston, suspension, || {
            seen += 1;
            seen == polls
        },),
        Err(ReactiveError::Cancelled)
    );
    assert_eq!(seen, polls);
    for (j, flow) in [&left, &right].into_iter().enumerate() {
        assert_eq!(flow.accepted.states, states[j]);
        assert_eq!(flow.flow().mesh(), &meshes[j]);
        assert_eq!(flow.flow().time(), time);
    }
    assert_eq!(piston, before_piston);
    let recovered =
        suspended_step(&mut left, &mut right, &mut piston, suspension, || false).unwrap();
    assert_eq!(recovered, expected);
    assert_eq!(left.accepted.states, reference_left.accepted.states);
    assert_eq!(right.accepted.states, reference_right.accepted.states);
    assert_eq!(piston, reference_piston);
}
