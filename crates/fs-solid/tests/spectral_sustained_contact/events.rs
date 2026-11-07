use super::*;
use fs_solid::dynamics::sustained::events::{CollinearCoulombContact, CoulombPhaseMode};

#[test]
fn spectral_stop_preserves_elastic_spectator_and_actual_modal_continuation() {
    with_cx(|cx, _| {
        let mut mass = [0.0; 16];
        for i in 0..4 {
            mass[5 * i] = 1.0;
        }
        let mut k = [0.0; 16];
        k[15] = 4.0;
        let assembly = assembly(&k, &mass, 4);
        let prepare = || {
            CompleteSpectralSolid::prepare_with_basis(
                &assembly,
                &[0.0, 0.0, 0.0, 4.0],
                &mass,
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
        let mut c = contact(1.0, 0.75, [0.0; 3]);
        c.velocity_weights = [&weights[0], &weights[1], &weights[2]];
        let law = CollinearCoulombContact {
            contact: c,
            kinetic_coefficient: 0.25,
            tangent_axis: [1.0, 0.0],
        };
        let load = [0.0, -8.0, 0.0, 0.3];
        free.advance(&load, cx).unwrap();
        let r = s
            .advance_collinear_coulomb_events(&law, [0.0; 3], &load, controls(), cx, |_| {
                Ok::<_, ()>(())
            })
            .unwrap();
        near(r.stopping_time_offset_s.unwrap(), 0.5, 1e-15);
        assert_eq!(r.phases[1].mode, CoulombPhaseMode::Sticking);
        near(s.displacement_m()[0], 0.25, 1e-10);
        near(s.velocity_mps()[0], 0.0, 1e-9);
        near(r.frictional_dissipation_j, 0.5, 1e-9);
        near(r.native.energy.balance_error_j, 0.0, 1e-10);
        near(r.native.energy.cumulative_balance_error_j, 0.0, 1e-10);
        near(r.energy_balance_error_j, 0.0, 1e-9);
        assert_eq!(
            s.displacement_m()[3].to_bits(),
            free.displacement_m()[3].to_bits()
        );
        assert_eq!(
            s.velocity_mps()[3].to_bits(),
            free.velocity_mps()[3].to_bits()
        );
        // Re-enter the ordinary modal flow. Its rigid coordinates must start
        // from the event result rather than the uncorrected free modal state.
        let next = s.advance(&[0.0, 0.0, 0.0, 0.3], cx).unwrap();
        free.advance(&load, cx).unwrap();
        near(s.displacement_m()[0], 0.25, 1e-10);
        near(s.velocity_mps()[0], 0.0, 1e-9);
        assert_eq!(
            s.displacement_m()[3].to_bits(),
            free.displacement_m()[3].to_bits()
        );
        assert_eq!(
            s.velocity_mps()[3].to_bits(),
            free.velocity_mps()[3].to_bits()
        );
        near(next.energy.cumulative_balance_error_j, 0.0, 1e-10);
    });
}

#[test]
fn spectral_reversal_retains_both_work_segments_and_one_native_clock_advance() {
    with_cx(|cx, _| {
        let mut s = point(0.5, [1.0, 0.0, 0.0], cx);
        let law = CollinearCoulombContact {
            contact: contact(0.5, 0.4, [0.0; 3]),
            kinetic_coefficient: 0.25,
            tangent_axis: [7.0, 0.0],
        };
        let r = s
            .advance_collinear_coulomb_events(
                &law,
                [0.0; 3],
                &[-4.0, -8.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        near(r.stopping_time_offset_s.unwrap(), 1.0 / 6.0, 1e-15);
        near(s.velocity_mps()[0], -2.0 / 3.0, 1e-9);
        near(s.displacement_m()[0], -1.0 / 36.0, 1e-10);
        near(r.sliding_distance_m, 7.0 / 36.0, 1e-10);
        near(r.frictional_dissipation_j, 7.0 / 18.0, 1e-9);
        near(r.external_work_j, 1.0 / 9.0, 1e-9);
        near(r.energy_balance_error_j, 0.0, 1e-9);
        near(r.native.energy.cumulative_balance_error_j, 0.0, 1e-10);
        assert_eq!(r.native.step, 1);
        assert_eq!(s.time_s(), 0.5);
        assert_eq!(r.phases.len(), 2);
        assert!(r.phases.iter().all(|p| p.mode == CoulombPhaseMode::Sliding));
    });
}

#[test]
fn spectral_nonrigid_trace_and_event_path_refusal_leave_modal_state_unchanged() {
    with_cx(|cx, _| {
        let k = [4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut s = CompleteSpectralSolid::prepare_with_basis(
            &assembly(&k, &I, 3),
            &[4.0, 0.0, 0.0],
            &I,
            1.0,
            &[0.0; 3],
            &[1.0, 0.0, 0.0],
            native_controls(),
            cx,
        )
        .unwrap();
        let law = CollinearCoulombContact {
            contact: contact(1.0, 0.75, [0.0; 3]),
            kinetic_coefficient: 0.25,
            tangent_axis: [1.0, 0.0],
        };
        let before = bits(&s);
        let r = s.advance_collinear_coulomb_events(
            &law,
            [0.0; 3],
            &[0.0, -8.0, 0.0],
            controls(),
            cx,
            |_| -> Result<(), ()> { panic!("nonrigid trace") },
        );
        assert!(matches!(r, Err(SolidSustainedError::Native(_))));
        assert_eq!(bits(&s), before);
        let mut s = point(1.0, [1.0, 0.0, 0.0], cx);
        let before = bits(&s);
        let r = s.advance_collinear_coulomb_events(
            &law,
            [0.0; 3],
            &[0.0, -8.0, 0.0],
            controls(),
            cx,
            |_| Err::<(), _>("uncovered feature"),
        );
        assert!(matches!(
            r,
            Err(SolidSustainedError::Path("uncovered feature"))
        ));
        assert_eq!(bits(&s), before);
        let mut reference = point(1.0, [1.0, 0.0, 0.0], cx);
        s.advance(&[0.0, -8.0, 0.0], cx).unwrap();
        reference.advance(&[0.0, -8.0, 0.0], cx).unwrap();
        assert_eq!(bits(&s), bits(&reference));
    });
}
