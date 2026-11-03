//! G0/G1: accepted offsets and exact clocks retain the original damaged owner.

use super::*;

#[test]
fn accepted_offsets_and_common_clock_keep_the_original_cohesive_evolution() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law();
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let prescribed = [[0.0002, 0.00003, -0.00004]];
        let state = system
            .bonded_state()
            .solve(&[0.026], &prescribed, config().solve, cx)
            .unwrap()
            .state;
        let old_point = state.points()[0];
        assert!(old_point.damage() > 0.0 && old_point.damage() < 1.0);
        let old_u = state.displacements_m().to_vec();
        let start = 0.2 / 64.0;
        let target = 0.9 / 64.0;
        let h = target - start;
        // Re-adding this rounded remainder misses the target by one ULP.
        assert_ne!(start + h, target);
        let motion = state
            .into_dynamics(&mass, vec![-0.001], start, config())
            .unwrap();
        let offsets = motion.prescribed_jumps_m_controlled(&mut || false).unwrap();
        assert_eq!(offsets, motion.prescribed_jumps_m(cx).unwrap());
        for axis in 0..3 {
            close(offsets[0][axis], prescribed[0][axis], 1e-15);
        }
        let controls = motion.dynamics_config();
        assert_eq!(controls.maximum_step_s, config().maximum_step_s);
        assert_eq!(controls.maximum_steps, config().maximum_steps);
        assert_eq!(
            controls.solve.max_operator_evaluations,
            config().solve.max_operator_evaluations
        );

        let original = motion
            .advance_controlled(&[0.0], &prescribed, h, &mut || false)
            .unwrap();
        let timed = motion
            .advance_to_time_controlled(&[0.0], &prescribed, target, &mut || false)
            .unwrap();
        let scoped = motion
            .advance_to_time(&[0.0], &prescribed, target, cx)
            .unwrap();
        assert_ne!(original.motion.time_s().to_bits(), target.to_bits());
        assert_eq!(timed.motion.time_s().to_bits(), target.to_bits());
        assert_eq!(scoped.motion.time_s().to_bits(), target.to_bits());
        for candidate in [&timed, &scoped] {
            assert_eq!(
                candidate.motion.state().displacements_m(),
                original.motion.state().displacements_m()
            );
            assert_eq!(
                candidate.motion.velocity_m_s(),
                original.motion.velocity_m_s()
            );
            assert_eq!(
                candidate.motion.state().points(),
                original.motion.state().points()
            );
            assert_eq!(candidate.motion.steps(), 1);
            assert!(core::ptr::eq(candidate.motion.state().system(), &system));
            assert!(core::ptr::eq(
                candidate.motion.mass_operator(),
                motion.mass_operator()
            ));
            assert_eq!(
                candidate.motion.state().points()[0].damage(),
                old_point.damage()
            );
            assert_eq!(candidate.cohesive_energy_j.dissipated, 0.0);
            close(candidate.energy_residual_j, 0.0, 1e-12);
            let retained = candidate.motion.prescribed_jumps_m(cx).unwrap();
            for axis in 0..3 {
                close(retained[0][axis], prescribed[0][axis], 1e-15);
            }
        }
        assert_eq!(motion.state().points(), &[old_point]);
        assert_eq!(motion.state().displacements_m(), old_u);
        assert_eq!(motion.velocity_m_s(), &[-0.001]);
        assert_eq!(motion.time_s(), start);
        assert_eq!(motion.steps(), 0);
    });
}

#[test]
fn offset_and_target_time_refusals_leave_the_accepted_body_unchanged() {
    with_cx(|cx| {
        let assembled = scalar_tet(cx);
        let mass = CsrOp::symmetric(assembled.mass);
        let elastic = CsrOp::symmetric(assembled.stiffness);
        let law = law();
        let system = CohesiveSystem::new(&elastic, vec![patch(&law)], budget()).unwrap();
        let motion = system
            .bonded_state()
            .into_dynamics(&mass, vec![0.01], 2.0, config())
            .unwrap();
        let old_points = motion.state().points().to_vec();
        let mut total_polls = 0;
        motion
            .prescribed_jumps_m_controlled(&mut || {
                total_polls += 1;
                false
            })
            .unwrap();
        let mut polls = 0;
        assert!(matches!(
            motion.prescribed_jumps_m_controlled(&mut || {
                polls += 1;
                polls == total_polls
            }),
            Err(DynamicsError::Model(Error::Cancelled))
        ));
        assert!(matches!(
            motion.advance_to_time_controlled(&[0.0], &[[0.0; 3]], 2.01, &mut || true),
            Err(DynamicsError::Model(Error::Cancelled))
        ));
        for target in [2.0, 1.0, f64::NAN, f64::INFINITY] {
            assert!(
                motion
                    .advance_to_time(&[0.0], &[[0.0; 3]], target, cx)
                    .is_err()
            );
        }
        assert!(matches!(
            motion.advance_to_time(&[0.0], &[[0.0; 3]], 2.1, cx),
            Err(DynamicsError::Model(Error::Budget("time increment")))
        ));
        assert_eq!(motion.state().points(), old_points);
        assert_eq!(motion.state().displacements_m(), &[0.0]);
        assert_eq!(motion.velocity_m_s(), &[0.01]);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.steps(), 0);
        assert_eq!(motion.prescribed_jumps_m(cx).unwrap(), vec![[0.0; 3]]);
    });
}
