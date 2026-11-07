use super::*;
use fs_solid::dynamics::sustained::events::{
    CollinearCoulombContact, CoulombPhaseMode, SolidCoulombEventStep,
};

fn law(h: f64, surface: [f64; 3]) -> CollinearCoulombContact<'static> {
    CollinearCoulombContact {
        contact: contact(h, 0.6, surface),
        kinetic_coefficient: 0.3,
        tangent_axis: [1.0, 0.0],
    }
}

fn balance<P>(r: &SolidCoulombEventStep<P>) {
    near(r.native.energy.balance_error_j, 0.0, 1e-9);
    near(r.native.energy.cumulative_balance_error_j, 0.0, 1e-9);
    near(r.energy_balance_error_j, 0.0, 1e-9);
    near(
        r.native.energy.supplied_work_j,
        r.external_work_j
            + r.prescribed_surface_work_j
            + r.normal_contact_work_j
            + r.frictional_work_j,
        1e-9,
    );
    assert!(r.frictional_work_j <= 1e-10);
}

#[test]
fn event_stop_recaptures_static_friction_without_a_velocity_reset() {
    with_cx(|cx, _| {
        let mut s = point(1.0, [1.0, 0.0, 0.0], cx);
        let r = s
            .advance_collinear_coulomb_events(
                &law(1.0, [0.0; 3]),
                [20.0, -6.0, 0.0],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                |p| {
                    assert_eq!(p.before.time_s(), 0.0);
                    assert_eq!(p.phases.len(), 2);
                    Ok::<_, ()>("fixed planar point")
                },
            )
            .unwrap();
        assert_eq!(r.phases[0].mode, CoulombPhaseMode::Sliding);
        assert_eq!(r.phases[1].mode, CoulombPhaseMode::Sticking);
        near(r.stopping_time_offset_s.unwrap(), 1.0 / 3.0, 1e-15);
        assert!(r.stopping_speed_residual_mps.unwrap() <= controls().velocity_tolerance_mps);
        near(s.displacement_m()[0], 1.0 / 6.0, 1e-10);
        near(s.velocity_mps()[0], 0.0, 1e-9);
        near(r.frictional_dissipation_j, 1.0, 1e-9);
        near(r.sliding_distance_m, 1.0 / 6.0, 1e-10);
        assert_eq!(r.native.step, 1);
        assert_eq!(r.native.time_s, 1.0);
        assert_eq!(r.native.energy.dissipated_j, 0.0);
        balance(&r);
        let next = s
            .advance_collinear_coulomb_events(
                &law(1.0, [0.0; 3]),
                [20.0, 0.0, 0.0],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        near(s.displacement_m()[0], 1.0 / 6.0, 1e-10);
        near(s.velocity_mps()[0], 0.0, 1e-9);
        assert_eq!(next.native.step, 2);
        balance(&next);
    });
}

#[test]
fn event_capture_on_moving_surface_separates_transport_work_and_sliding_loss() {
    with_cx(|cx, _| {
        let mut s = point(1.0, [3.0, 0.0, 0.0], cx);
        let r = s
            .advance_collinear_coulomb_events(
                &law(1.0, [0.0, 1.0, 0.0]),
                [0.0; 3],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        near(r.stopping_time_offset_s.unwrap(), 2.0 / 3.0, 1e-15);
        near(s.velocity_mps()[0], 1.0, 1e-9);
        near(s.displacement_m()[0], 5.0 / 3.0, 1e-10);
        near(r.prescribed_surface_work_j, -4.0, 1e-9);
        near(r.frictional_dissipation_j, 4.0, 1e-9);
        near(s.kinetic_energy_j(), 1.0, 1e-9);
        balance(&r);
    });
}

#[test]
fn event_reversal_counts_both_sliding_distances_and_discards_static_trial() {
    with_cx(|cx, _| {
        let h = 0.25;
        let mut s = point(h, [1.0, 0.0, 0.0], cx);
        let r = s
            .advance_collinear_coulomb_events(
                &law(h, [0.0; 3]),
                [0.0; 3],
                &[-20.0, -20.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        let stop = 1.0 / 13.0;
        let remainder = h - stop;
        let forward = 0.5 * stop;
        let backward = 3.5 * remainder * remainder;
        near(r.stopping_time_offset_s.unwrap(), stop, 1e-15);
        assert_eq!(r.phases.len(), 2);
        assert!(r.phases.iter().all(|p| p.mode == CoulombPhaseMode::Sliding));
        near(r.phases[0].force_n[1], -6.0, 1e-8);
        near(r.phases[1].force_n[1], 6.0, 1e-8);
        near(s.velocity_mps()[0], -7.0 * remainder, 1e-9);
        near(s.displacement_m()[0], forward - backward, 1e-10);
        near(r.sliding_distance_m, forward + backward, 1e-10);
        near(r.frictional_dissipation_j, 6.0 * (forward + backward), 1e-9);
        assert!(r.frictional_dissipation_j > 6.0 * s.displacement_m()[0].abs());
        assert_eq!(r.static_capture_trials, 1);
        assert!(
            r.response_evaluations > r.phases.iter().map(|p| p.solve.response_evaluations).sum()
        );
        balance(&r);
    });
}

#[test]
fn event_static_threshold_and_hysteretic_recapture_use_distinct_coefficients() {
    with_cx(|cx, _| {
        for (load, acceleration) in [(9.0, 0.0), (14.0, 4.0)] {
            let mut s = point(0.5, [0.0; 3], cx);
            let r = s
                .advance_collinear_coulomb_events(
                    &law(0.5, [0.0; 3]),
                    [0.0; 3],
                    &[load, -20.0, 0.0],
                    controls(),
                    cx,
                    |_| Ok::<_, ()>(()),
                )
                .unwrap();
            assert_eq!(r.phases.len(), 1);
            near(s.velocity_mps()[0], 0.5 * acceleration, 1e-9);
            near(s.displacement_m()[0], 0.125 * acceleration, 1e-10);
            near(
                r.phases[0].friction_coefficient,
                if load == 9.0 { 0.6 } else { 0.3 },
                0.0,
            );
            balance(&r);
        }
        // Kinetic friction plus the positive source stops negative incoming
        // motion; the stronger static cone can then hold that same source.
        let mut s = point(0.5, [-1.0, 0.0, 0.0], cx);
        let r = s
            .advance_collinear_coulomb_events(
                &law(0.5, [0.0; 3]),
                [0.0; 3],
                &[9.0, -20.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        near(r.stopping_time_offset_s.unwrap(), 1.0 / 7.5, 1e-15);
        near(r.phases[0].force_n[1], 6.0, 1e-8);
        near(r.phases[1].force_n[1], -9.0, 1e-8);
        assert_eq!(r.phases[1].mode, CoulombPhaseMode::Sticking);
        near(s.velocity_mps()[0], 0.0, 1e-9);
        balance(&r);
    });
}

#[test]
fn event_contact_keeps_original_off_diagonal_consistent_inertia() {
    with_cx(|cx, _| {
        let mass = [2.0, 0.5, 0.0, 0.5, 2.0, 0.0, 0.0, 0.0, 2.0];
        let mut s = ExactLinearSolid::prepare(
            &assembly(&[0.0; 9], &mass, 3),
            1.0,
            &[0.0; 3],
            &[1.0, 0.0, 0.0],
            native_controls(),
            cx,
        )
        .unwrap();
        let r = s
            .advance_collinear_coulomb_events(
                &law(1.0, [0.0; 3]),
                [0.0; 3],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        let deceleration = 6.0 / 2.15;
        near(r.stopping_time_offset_s.unwrap(), 1.0 / deceleration, 1e-14);
        near(r.phases[0].force_n[0], 20.0 - 0.5 * deceleration, 1e-8);
        near(r.phases[1].force_n[0], 20.0, 1e-8);
        near(s.displacement_m()[0], 0.5 / deceleration, 1e-10);
        near(s.velocity_mps()[0], 0.0, 1e-9);
        near(s.velocity_mps()[1], 0.0, 1e-9);
        near(r.frictional_dissipation_j, 1.0, 1e-9);
        balance(&r);
    });
}

#[test]
fn event_retains_complete_elastic_viscous_spectator_motion_and_loss() {
    with_cx(|cx, _| {
        let mut m = [0.0; 16];
        for i in 0..4 {
            m[5 * i] = if i < 3 { 2.0 } else { 1.0 };
        }
        let mut k = [0.0; 16];
        k[15] = 4.0;
        let assembly = assembly(&k, &m, 4);
        let mut c = Coo::new(4, 4);
        c.push(3, 3, 0.5);
        let c = c.assemble();
        let prepare = || {
            ExactLinearSolid::prepare_viscous(
                &assembly,
                &c,
                1.0,
                &[0.0, 0.0, 0.0, 0.1],
                &[1.0, 0.0, 0.0, 0.2],
                native_controls(),
                cx,
            )
            .unwrap()
        };
        let (mut s, mut free) = (prepare(), prepare());
        let weights = [
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        let mut law = law(1.0, [0.0; 3]);
        law.contact.velocity_weights = [&weights[0], &weights[1], &weights[2]];
        let load = [0.0, -20.0, 0.0, 0.3];
        let reference = free.advance(&load, cx).unwrap();
        let r = s
            .advance_collinear_coulomb_events(&law, [0.0; 3], &load, controls(), cx, |_| {
                Ok::<_, ()>(())
            })
            .unwrap();
        assert_eq!(
            s.displacement_m()[3].to_bits(),
            free.displacement_m()[3].to_bits()
        );
        assert_eq!(
            s.velocity_mps()[3].to_bits(),
            free.velocity_mps()[3].to_bits()
        );
        assert_eq!(
            r.native.energy.dissipated_j.to_bits(),
            reference.energy.dissipated_j.to_bits()
        );
        assert!(r.native.energy.dissipated_j > 0.0);
        near(r.frictional_dissipation_j, 1.0, 1e-9);
        balance(&r);
        // Continue through the original native flow to expose stale or clipped
        // hidden state, not only a coincident first published endpoint.
        s.advance(&load, cx).unwrap();
        free.advance(&load, cx).unwrap();
        assert_eq!(
            s.displacement_m()[3].to_bits(),
            free.displacement_m()[3].to_bits()
        );
        assert_eq!(
            s.velocity_mps()[3].to_bits(),
            free.velocity_mps()[3].to_bits()
        );
    });
}

#[test]
fn event_nonrigid_or_rotating_contact_refuses_before_path_admission_atomically() {
    with_cx(|cx, _| {
        for viscous in [false, true] {
            let mut matrix = [0.0; 9];
            matrix[0] = 0.5;
            let m = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
            let mut s = if viscous {
                let c = assembly(&matrix, &m, 3).stiffness;
                ExactLinearSolid::prepare_viscous(
                    &assembly(&[0.0; 9], &m, 3),
                    &c,
                    1.0,
                    &[0.0; 3],
                    &[1.0, 0.0, 0.0],
                    native_controls(),
                    cx,
                )
                .unwrap()
            } else {
                ExactLinearSolid::prepare(
                    &assembly(&matrix, &m, 3),
                    1.0,
                    &[0.0; 3],
                    &[1.0, 0.0, 0.0],
                    native_controls(),
                    cx,
                )
                .unwrap()
            };
            let before = bits(&s);
            let r = s.advance_collinear_coulomb_events(
                &law(1.0, [0.0; 3]),
                [0.0; 3],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                |_| -> Result<(), ()> { panic!("unsupported trace") },
            );
            assert!(matches!(r, Err(SolidSustainedError::Invalid(_))));
            assert_eq!(bits(&s), before);
        }
        let mut s = point(1.0, [1.0, 0.0, 0.1], cx);
        let before = bits(&s);
        let r = s.advance_collinear_coulomb_events(
            &law(1.0, [0.0; 3]),
            [0.0; 3],
            &[0.0, -20.0, 0.0],
            controls(),
            cx,
            |_| -> Result<(), ()> { panic!("rotating trace") },
        );
        assert!(matches!(r, Err(SolidSustainedError::ContactMode { .. })));
        assert_eq!(bits(&s), before);
    });
}

#[test]
fn event_global_budget_path_refusal_and_late_cancellation_leave_no_partial_stop() {
    with_cx(|cx, gate| {
        let law = law(1.0, [0.0; 3]);
        let load = [0.0, -20.0, 0.0];
        let mut reference = point(1.0, [1.0, 0.0, 0.0], cx);
        let accepted = reference
            .advance_collinear_coulomb_events(&law, [0.0; 3], &load, controls(), cx, |_| {
                Ok::<_, ()>(())
            })
            .unwrap();
        let mut controls = controls();
        controls.solve.maximum_response_evaluations =
            accepted.phases[0].solve.response_evaluations + 1;
        let mut s = point(1.0, [1.0, 0.0, 0.0], cx);
        let before = bits(&s);
        let r = s.advance_collinear_coulomb_events(
            &law,
            [0.0; 3],
            &load,
            controls,
            cx,
            |_| -> Result<(), ()> { panic!("spent event budget") },
        );
        assert!(matches!(r, Err(SolidSustainedError::Budget { .. })));
        assert_eq!(bits(&s), before);
        let r = s.advance_collinear_coulomb_events(
            &law,
            [0.0; 3],
            &load,
            super::controls(),
            cx,
            |_| Err::<(), _>("changed contact feature"),
        );
        assert!(matches!(
            r,
            Err(SolidSustainedError::Path("changed contact feature"))
        ));
        assert_eq!(bits(&s), before);
        let r = s.advance_collinear_coulomb_events(
            &law,
            [0.0; 3],
            &load,
            super::controls(),
            cx,
            |_| {
                gate.request();
                Ok::<_, ()>(())
            },
        );
        assert!(matches!(r, Err(SolidSustainedError::Cancelled)));
        assert_eq!(bits(&s), before);
    });
}
