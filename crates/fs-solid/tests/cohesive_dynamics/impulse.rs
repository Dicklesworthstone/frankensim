//! G0: original consistent inertia, irreversible history and signed event work.
use super::*;
use fs_solid::cohesive::dynamics::impulse::{ImpulseControls, ImpulseError};
use fs_solver::CheckedCgConfig;
use fs_sparse::precond::IdentityPrecond;

fn impulse_controls() -> ImpulseControls {
    ImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 32,
            max_corrections: 2,
        },
        max_workspace_entries: 1024,
        absolute_momentum_tolerance_ns: 1e-13,
        relative_momentum_tolerance: 1e-10,
        absolute_work_tolerance_j: 1e-13,
        relative_work_tolerance: 1e-10,
    }
}

#[test]
fn an_impulse_changes_original_mass_momentum_without_restarting_a_damaged_motion() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law();
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let u0 = 0.002;
        let traction = (0.004 - u0) / 0.003;
        let accepted = system
            .bonded_state()
            .solve(
                &[10.0 * u0 + 0.01 * traction],
                &[[0.0; 3]],
                config().solve,
                cx,
            )
            .unwrap()
            .state;
        assert!((0.0..1.0).contains(&accepted.points()[0].damage()));
        let motion = accepted
            .into_dynamics(&mass, vec![-0.001], 3.0, config())
            .unwrap()
            .advance(&[0.0], &[[0.0; 3]], 0.001, cx)
            .unwrap()
            .motion;
        let old_point = motion.state().points()[0];
        let old_u = motion.state().displacements_m()[0];
        let old_v = motion.velocity_m_s()[0];
        assert_eq!(motion.steps(), 1);

        // Oppose the current velocity without reversing it: impulse work is
        // negative, not spurious fracture dissipation or a positive-only loss.
        let impulse = 5e-5;
        let expected_v = old_v + impulse / 0.1;
        assert!(old_v < expected_v && expected_v < 0.0);
        let event = motion
            .impulse_candidate(cx, &[impulse], &IdentityPrecond, impulse_controls())
            .unwrap();
        close(event.motion.velocity_m_s()[0], expected_v, 1e-14);
        close(
            0.1 * (event.motion.velocity_m_s()[0] - old_v),
            impulse,
            1e-15,
        );
        assert!(event.report.applied_impulse_work_j < 0.0);
        close(
            event.report.applied_impulse_work_j,
            0.5 * (old_v + expected_v) * impulse,
            1e-16,
        );
        close(
            event.report.kinetic_energy_before_j,
            0.05 * old_v * old_v,
            1e-16,
        );
        close(
            event.report.kinetic_energy_after_j,
            0.05 * expected_v * expected_v,
            1e-16,
        );
        close(
            event.report.kinetic_energy_change_j,
            event.report.applied_impulse_work_j,
            1e-16,
        );
        close(event.report.work_residual_j, 0.0, 1e-16);
        assert!(event.report.solve.converged_euclidean());
        assert_eq!(event.motion.state().displacements_m(), &[old_u]);
        assert_eq!(event.motion.state().points(), &[old_point]);
        assert!(core::ptr::eq(event.motion.state().system(), &system));
        assert_eq!(event.motion.time_s(), motion.time_s());
        assert_eq!(event.motion.steps(), motion.steps());

        // The following physical interval uses the same softened stiffness,
        // native mass, material history and accepted clock/step budget.
        let dt = 0.001;
        let stiffness = 10.0 + (1.0 - old_point.damage()) * 10.0;
        let expected_u = ((2.0 * 0.1 / (dt * dt) - 0.5 * stiffness) * old_u
            + 2.0 * 0.1 * expected_v / dt)
            / (2.0 * 0.1 / (dt * dt) + 0.5 * stiffness);
        let next = event.motion.advance(&[0.0], &[[0.0; 3]], dt, cx).unwrap();
        close(next.motion.state().displacements_m()[0], expected_u, 1e-12);
        close(
            next.kinetic_before_j,
            event.report.kinetic_energy_after_j,
            1e-16,
        );
        assert_eq!(next.motion.state().points()[0].damage(), old_point.damage());
        assert_eq!(
            next.motion.state().points()[0].dissipated_j_m2(),
            old_point.dissipated_j_m2()
        );
        assert_eq!(next.cohesive_energy_j.dissipated, 0.0);
        assert_eq!(next.viscous_dissipated_j, 0.0);
        assert_eq!(next.motion.steps(), motion.steps() + 1);
        close(next.energy_residual_j, 0.0, 1e-12);

        // Neither the event nor its continuation mutates the accepted input.
        assert_eq!(motion.state().points(), &[old_point]);
        assert_eq!(motion.state().displacements_m(), &[old_u]);
        assert_eq!(motion.velocity_m_s(), &[old_v]);
        assert_eq!(motion.steps(), 1);
    });
}

#[test]
fn impulse_refusals_preserve_the_accepted_motion_and_allow_exact_retry() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law();
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let motion = system
            .bonded_state()
            .into_dynamics(&mass, vec![0.01], 2.0, config())
            .unwrap()
            .advance(&[0.0], &[[0.0; 3]], 0.001, cx)
            .unwrap()
            .motion;
        let saved_point = motion.state().points()[0];
        let saved_u = motion.state().displacements_m().to_vec();
        let saved_v = motion.velocity_m_s().to_vec();
        let mut complete_polls = 0;
        let expected = motion
            .impulse_candidate_controlled(
                &[0.001],
                &IdentityPrecond,
                impulse_controls(),
                &mut || {
                    complete_polls += 1;
                    false
                },
            )
            .unwrap();
        let mut polls = 0;
        assert!(matches!(
            motion.impulse_candidate_controlled(
                &[0.001],
                &IdentityPrecond,
                impulse_controls(),
                &mut || {
                    polls += 1;
                    polls == complete_polls
                },
            ),
            Err(ImpulseError::Cancelled)
        ));
        let mut limited = impulse_controls();
        limited.max_workspace_entries = 0;
        assert!(
            motion
                .impulse_candidate(cx, &[0.001], &IdentityPrecond, limited)
                .is_err()
        );
        assert!(
            motion
                .impulse_candidate(cx, &[1.0], &IdentityPrecond, impulse_controls())
                .is_err()
        ); // Original velocity applicability limit remains active.
        assert_eq!(motion.state().points(), &[saved_point]);
        assert_eq!(motion.state().displacements_m(), saved_u);
        assert_eq!(motion.velocity_m_s(), saved_v);
        assert_eq!(motion.time_s(), 2.001);
        assert_eq!(motion.steps(), 1);
        let retry = motion
            .impulse_candidate(cx, &[0.001], &IdentityPrecond, impulse_controls())
            .unwrap();
        assert_eq!(retry.motion.velocity_m_s(), expected.motion.velocity_m_s());
        assert_eq!(
            retry.motion.state().points(),
            expected.motion.state().points()
        );
        assert_eq!(
            retry.report.applied_impulse_work_j,
            expected.report.applied_impulse_work_j
        );
    });
}
