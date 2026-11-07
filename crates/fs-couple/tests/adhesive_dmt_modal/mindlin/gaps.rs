//! G1/G3 signed native release, attractive flight and reattachment.
use super::*;

#[test]
fn complete_adhesive_release_flight_and_recontact_keep_actual_motion_and_work() {
    let source = source_fixture::source(0.01, Some(MU), "initial").unwrap();
    let zero = source
        .normal_source()
        .evaluate_contact_load(0.0)
        .unwrap()
        .indentation_m;
    for cross in [false, true] {
        let mut system = make(&source, 0.6 * zero, cross);
        let initial = system.total_energy_j().unwrap();
        let h = system.sample_period_s();
        let (s, v) = response(cross, h);
        let mut supplied = 0.0;
        let mut losses = 0.0;
        let mut detached = false;
        let mut reattached = false;
        let mut attractive_flight = false;
        for (target, slip) in [
            (0.8 * zero, 0.003 * Z),
            (-0.8 * Z, 0.004 * Z),
            (-1.4 * Z, -0.002 * Z),
            (-0.5 * Z, 0.007 * Z),
            (0.9 * zero, 0.003 * Z),
            (zero, -0.002 * Z),
        ] {
            let old = system.state().indentation_m();
            let old_slip = system.state().slip_m();
            let original = source
                .trial_coupled_path(
                    system.state(),
                    &SourceMindlinPathInput {
                        indentation_m: target,
                        free_slip_increment_m: slip,
                        dt_s: h,
                        compliance_m_per_n: 0.0,
                        solve: config().shear,
                    },
                    &mut || false,
                )
                .unwrap();
            let f = [
                original.mean_normal_force_n,
                original.tangential.average_force_n,
            ];
            let (q0, v0) = relative(&system);
            let (free, free_v) = released(q0, v0, cross, h);
            let desired = [target - old, slip];
            let rate = core::array::from_fn(|i| {
                (desired[i] - free[i] + s[i][0] * f[0] + s[i][1] * f[1]) / h
            });
            let frame = system.step(&[0.0; 6], rate, [-1.5 * Z, Z]).unwrap().clone();
            close(frame.position_m[0], target, 5e-23, 5e-9);
            close(frame.position_m[1], old_slip + slip, 5e-23, 5e-9);
            close(
                frame.microslip_dissipation_j,
                original.tangential.dissipated_j,
                1e-28,
                3e-8,
            );
            let fraction = if old > 0.0 && target > 0.0 {
                1.0
            } else if old > 0.0 {
                old / (old - target)
            } else if target > 0.0 {
                target / (target - old)
            } else {
                0.0
            };
            close(frame.contact_fraction, fraction, 2e-10, 2e-8);
            for i in 0..2 {
                close(
                    frame.contact_interval_fraction[i],
                    original.contact_interval_fraction[i],
                    2e-10,
                    2e-8,
                );
                close(frame.mean_force_n[i], f[i], 5e-18, 5e-8);
            }
            detached |= old > 0.0 && target < 0.0;
            reattached |= old < 0.0 && target > 0.0;
            if target < 0.0 {
                assert_eq!(system.state().endpoint_force_n(), 0.0);
                assert_eq!(system.state().tangential_energy_j(), 0.0);
                assert_eq!(system.state().contact_radius_m(), 0.0);
            }
            if old < 0.0 && target < 0.0 {
                assert_eq!(frame.mean_force_n[1], 0.0);
                assert_eq!(frame.microslip_dissipation_j, 0.0);
                attractive_flight |= frame.mean_force_n[0] < 0.0;
            }
            let (q1, v1) = relative(&system);
            for i in 0..3 {
                close(
                    q1[i] - q0[i],
                    free[i] - s[i][0] * f[0] - s[i][1] * f[1],
                    8e-23,
                    8e-8,
                );
                close(
                    v1[i],
                    free_v[i] - v[i][0] * f[0] - v[i][1] * f[1],
                    8e-16,
                    8e-8,
                );
                close(momentum(&system)[i], 0.0, 1e-24, 0.0);
            }
            supplied += frame.prescribed_work_j;
            losses += frame.microslip_dissipation_j
                + frame.network_dissipation_j
                + frame.numerical_dissipation_j;
            close(
                system.total_energy_j().unwrap() + losses,
                initial + supplied,
                2e-25,
                8e-8,
            );
            assert!(frame.contact_work_residual_j.abs() <= frame.contact_work_tolerance_j);
        }
        assert!(detached && reattached && attractive_flight);
        assert_eq!(system.samples_rendered(), 6);
    }
}

#[test]
fn initially_separated_native_bodies_recontact_only_after_receiver_admission() {
    let source = source_fixture::source(0.01, Some(MU), "initial").unwrap();
    let initial = -1.2 * Z;
    let target = 0.5
        * source
            .normal_source()
            .evaluate_contact_load(0.0)
            .unwrap()
            .indentation_m;
    let mut system = make(&source, initial, true);
    let mut replay = make(&source, initial, true);
    let native = all_states(&system);
    let state = system.state().clone();
    let h = system.sample_period_s();
    let slip = 0.003 * Z;
    let direct = source
        .trial_coupled_path(
            &state,
            &SourceMindlinPathInput {
                indentation_m: target,
                free_slip_increment_m: slip,
                dt_s: h,
                compliance_m_per_n: 0.0,
                solve: config().shear,
            },
            &mut || false,
        )
        .unwrap();
    let c = system.compliance_m_per_n();
    let f = [
        direct.mean_normal_force_n,
        direct.tangential.average_force_n,
    ];
    let rate = [
        (target - initial + c[0][0] * f[0] + c[0][1] * f[1]) / h,
        (slip + c[1][0] * f[0] + c[1][1] * f[1]) / h,
    ];
    let gate = CancelGate::new_clock_free();
    let rejected: Result<(), AdhesiveDmtMindlinError> =
        system.step_with_admission(&[0.0; 6], rate, [-1.5 * Z, Z], &gate, |frame, _| {
            assert!(frame.contact_fraction > 0.0 && frame.contact_fraction < 1.0);
            assert!(frame.endpoint_shear_force_n > 0.0);
            gate.request();
            Ok(())
        });
    assert!(rejected.is_err());
    assert_eq!(all_states(&system), native);
    assert_eq!(system.state(), &state);
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
    let accepted = system.step(&[0.0; 6], rate, [-1.5 * Z, Z]).unwrap().clone();
    assert_eq!(
        replay.step(&[0.0; 6], rate, [-1.5 * Z, Z]).unwrap(),
        &accepted
    );
    assert_eq!(all_states(&system), all_states(&replay));
    assert_eq!(system.state(), replay.state());
    close(system.state().slip_m(), slip, 3e-23, 3e-9);
}
