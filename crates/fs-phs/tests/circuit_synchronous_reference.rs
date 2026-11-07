//! Source-owned torque allocation, including nonlinear multiple branches.
use fs_phs::PreparedStepError;
use fs_phs::circuit::inductance::InductorFluxCurve;
use fs_phs::circuit::synchronous::{SynchronousMotor, SynchronousMotorParameters};

fn parameters() -> SynchronousMotorParameters {
    SynchronousMotorParameters {
        d_inductance_h: 0.02,
        q_inductance_h: 0.04,
        winding_resistance_ohm: 0.1,
        pm_flux_linkage_wb: 0.2,
        pole_pairs: 3,
        rotor_inertia_kg_m2: 0.1,
        stator_inertia_kg_m2: 0.2,
        bearing_drag_nm_s_per_rad: 0.01,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn nonlinear() -> SynchronousMotor {
    let d = InductorFluxCurve::try_new(vec![(-2.0, -2.0), (0.0, 0.0), (2.0, 2.0)]).unwrap();
    // Strictly monotone flux/current source, but torque at fixed id=1 has
    // multiple affine branches. These are manufactured constitutive data.
    let q = InductorFluxCurve::try_new(vec![
        (-4.0, -4.0),
        (0.0, 0.0),
        (0.2, 1.0),
        (2.2, 2.0),
        (2.4, 3.0),
        (4.4, 4.0),
    ])
    .unwrap();
    SynchronousMotor::with_flux_curves(
        SynchronousMotorParameters {
            pole_pairs: 2,
            pm_flux_linkage_wb: 0.0,
            ..parameters()
        },
        d,
        q,
    )
    .unwrap()
}

#[test]
fn linear_torque_limits_use_the_current_circle_and_both_torque_signs() {
    let motor = SynchronousMotor::new(parameters()).unwrap();
    for sign in [-1.0, 1.0] {
        let r = motor
            .current_reference_for_torque_controlled(sign * 1.44, -2.0, 5.0, || false)
            .unwrap();
        assert!(!r.torque_limited);
        close(r.current_dq_a[0], -2.0, 0.0);
        close(r.current_dq_a[1], sign * 2.0, 2e-15);
        close(r.achievable_torque_nm, sign * 1.44, 2e-15);
        let limited = motor
            .current_reference_for_torque_controlled(sign * 100.0, -2.0, 5.0, || false)
            .unwrap();
        assert!(limited.torque_limited);
        close(limited.current_dq_a[1], sign * 21.0_f64.sqrt(), 2e-15);
        close(
            limited.current_dq_a[0].hypot(limited.current_dq_a[1]),
            5.0,
            2e-15,
        );
        let x = motor.state(limited.current_dq_a, 0.0, [0.0; 2]).unwrap();
        close(
            motor.operating_point(&x).unwrap().electromagnetic_torque_nm,
            limited.achievable_torque_nm,
            2e-15,
        );
    }
    let zero = motor
        .current_reference_for_torque_controlled(10.0, 0.0, 0.0, || false)
        .unwrap();
    assert_eq!(zero.current_dq_a, [0.0; 2]);
    assert!(zero.torque_limited);
    let edge = motor
        .current_reference_for_torque_controlled(10.0, 5.0, 5.0, || false)
        .unwrap();
    assert_eq!(edge.current_dq_a, [5.0, 0.0]);
    assert!(edge.torque_limited);
}

#[test]
fn nonlinear_source_search_selects_the_least_current_of_all_torque_branches() {
    let motor = nonlinear();
    let result = motor
        .current_reference_for_torque_controlled(0.8, 1.0, 5.0, || false)
        .unwrap();
    assert!(!result.torque_limited);
    // Other roots are iq=1.4,2.75,3.2. The original curve's first branch gives
    // lambda_q=.2*iq, hence tau=2*(iq-.2*iq) and the minimum-|iq| root is .5.
    close(result.current_dq_a[1], 0.5, 5e-16);
    close(result.achievable_torque_nm, 0.8, 5e-16);
    let x = motor.state(result.current_dq_a, 0.0, [0.0; 2]).unwrap();
    close(
        motor.operating_point(&x).unwrap().electromagnetic_torque_nm,
        0.8,
        5e-16,
    );

    let upper = motor
        .current_reference_for_torque_controlled(1.8, 1.0, 20.0, || false)
        .unwrap();
    assert!(upper.torque_limited);
    close(upper.current_dq_a[1], 1.0, 0.0);
    close(upper.achievable_torque_nm, 1.6, 1e-15);
    let lower = motor
        .current_reference_for_torque_controlled(-1.0, 1.0, 20.0, || false)
        .unwrap();
    assert!(lower.torque_limited);
    close(lower.current_dq_a[1], 4.0, 0.0);
    close(lower.achievable_torque_nm, -0.8, 2e-15);
    let circle = motor
        .current_reference_for_torque_controlled(0.8, 1.0, 1.02, || false)
        .unwrap();
    assert!(circle.torque_limited);
    close(
        circle.current_dq_a[1],
        (1.02_f64.powi(2) - 1.0).sqrt(),
        6e-16,
    );
    close(
        circle.achievable_torque_nm,
        1.6 * circle.current_dq_a[1],
        1e-15,
    );
}

#[test]
fn pure_reluctance_and_zero_authority_do_not_invent_pm_torque() {
    let p = SynchronousMotorParameters {
        pm_flux_linkage_wb: 0.0,
        ..parameters()
    };
    let reluctance = SynchronousMotor::new(p).unwrap();
    let r = reluctance
        .current_reference_for_torque_controlled(0.24, 2.0, 5.0, || false)
        .unwrap();
    close(r.current_dq_a[1], -2.0, 2e-15);
    assert!(!r.torque_limited);
    let no_authority = SynchronousMotor::new(SynchronousMotorParameters {
        q_inductance_h: p.d_inductance_h,
        ..p
    })
    .unwrap();
    for requested in [-1.0, 0.0, 1.0] {
        let r = no_authority
            .current_reference_for_torque_controlled(requested, 1.0, 5.0, || false)
            .unwrap();
        assert_eq!(r.current_dq_a, [1.0, 0.0]);
        assert_eq!(r.achievable_torque_nm, 0.0);
        assert_eq!(r.torque_limited, requested != 0.0);
    }
}

#[test]
fn unsupported_d_current_and_cancelled_search_leave_the_original_law_reusable() {
    let motor = nonlinear();
    for (torque, d, limit) in [
        (f64::NAN, 0.0, 1.0),
        (0.0, 2.1, 5.0),
        (0.0, 1.0, 0.5),
        (0.0, 0.0, -1.0),
    ] {
        assert!(
            motor
                .current_reference_for_torque_controlled(torque, d, limit, || false)
                .is_err()
        );
    }
    let mut calls = 0;
    let original = motor
        .current_reference_for_torque_controlled(0.8, 1.0, 5.0, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls] {
        let mut count = 0;
        assert!(matches!(
            motor.current_reference_for_torque_controlled(0.8, 1.0, 5.0, || {
                count += 1;
                count == stop
            }),
            Err(PreparedStepError::Cancelled)
        ));
    }
    assert_eq!(
        motor
            .current_reference_for_torque_controlled(0.8, 1.0, 5.0, || false)
            .unwrap(),
        original
    );
}
