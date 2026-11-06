use fs_contact::tangential::impact::multiple::relaxation::{
    CoupledFrictionFlowControls, CoupledFrictionFlowError, CoupledFrictionFlowInput,
    CoupledFrictionFlowReport, advance_coupled_friction,
};
use fs_contact::tangential::impact::multiple::{
    CoupledFrictionImpulseControls, CoupledFrictionImpulseError, MAX_FRICTION_IMPULSE_PORTS,
};

fn controls() -> CoupledFrictionFlowControls {
    CoupledFrictionFlowControls {
        projection: CoupledFrictionImpulseControls {
            maximum_iterations: 128,
            impulse_tolerance_ns: 2e-11,
            velocity_tolerance_m_s: 2e-11,
            energy_tolerance_j: 2e-10,
            relative_matrix_tolerance: 1e-12,
            maximum_impulse_ns: 100.0,
        },
        maximum_segments: 64,
        acceleration_tolerance_m_s2: 2e-11,
        maximum_force_n: 100.0,
    }
}

fn flow(v: &[f64], w: &[f64], c: &[f64], d: &[f64], t: f64) -> CoupledFrictionFlowReport {
    advance_coupled_friction(
        &CoupledFrictionFlowInput {
            velocity_before_m_s: v,
            mobility_per_kg: w,
            capacity_n: c,
            drive_n: d,
            duration_s: t,
        },
        controls(),
        || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}

fn audit(out: &CoupledFrictionFlowReport) {
    near(out.energy_residual_j, 0.0, 2e-11);
    near(
        out.kinetic_energy_change_j,
        out.total_drive_work_j - out.dissipated_j,
        2e-11,
    );
    assert!(out.dissipated_j >= 0.0);
    for i in 0..out.count {
        assert!(out.friction_work_j[i] <= 0.0);
        assert!(out.distance_m[i] >= out.displacement_m[i].abs());
        near(out.velocity_residual_m_s[i], 0.0, 2e-11);
        near(out.impulse_residual_ns[i], 0.0, 2e-11);
        assert!(out.impulse_capacity_residual_ns[i] <= 2e-11);
    }
}

#[test]
fn g1_scalar_motion_has_exact_finite_arrest_without_a_capture_speed() {
    // m=1/2, F=1, v0=3: stop at 3/2 s after 9/4 m, losing 9/4 J.
    for duration in [0.0_f64, 0.25, 1.5, 4.0] {
        let out = flow(&[3.0], &[2.0], &[1.0], &[0.0], duration);
        let moving = duration.min(1.5);
        let endpoint = 3.0 - 2.0 * moving;
        let distance = 3.0 * moving - moving * moving;
        near(out.velocity_after_m_s[0], endpoint, 2e-14);
        near(out.friction_impulse_ns[0], -moving, 2e-14);
        near(out.displacement_m[0], distance, 2e-14);
        near(out.distance_m[0], distance, 2e-14);
        near(out.dissipated_j, distance, 2e-14);
        near(
            out.kinetic_energy_change_j,
            0.25 * (endpoint * endpoint - 9.0),
            2e-14,
        );
        audit(&out);
    }
    let tiny = flow(&[1e-15], &[1.0], &[1.0], &[0.0], 1.0);
    near(tiny.friction_impulse_ns[0], -1e-15, 1e-30);
    near(tiny.distance_m[0], 0.5e-30, 1e-44);
    assert!(tiny.dissipated_j > 0.0);
    assert_eq!(tiny.flow_velocity_after_m_s[0], 0.0);
    assert_eq!(tiny.stops, 1);
}

#[test]
fn g1_held_drives_stop_hold_and_reverse_with_physical_source_work() {
    // Opposing initial motion stops at t=1/3. The drive then exceeds the
    // same friction bound, so v=t-1/3 on the second part of the path.
    let reversed = flow(&[-1.0], &[1.0], &[1.0], &[2.0], 2.0);
    near(reversed.velocity_after_m_s[0], 5.0 / 3.0, 2e-14);
    near(reversed.friction_impulse_ns[0], -4.0 / 3.0, 2e-14);
    near(reversed.drive_impulse_ns[0], 4.0, 2e-14);
    near(reversed.displacement_m[0], 11.0 / 9.0, 2e-14);
    near(reversed.distance_m[0], 14.0 / 9.0, 2e-14);
    near(reversed.total_drive_work_j, 22.0 / 9.0, 2e-14);
    near(reversed.dissipated_j, 14.0 / 9.0, 2e-14);
    near(reversed.kinetic_energy_change_j, 8.0 / 9.0, 2e-14);
    assert_eq!(reversed.stops, 1);
    assert_eq!(reversed.segments, 2);
    audit(&reversed);

    // A smaller drive holds after the original negative slip reaches zero.
    let held = flow(&[-2.0], &[1.0], &[1.0], &[0.25], 3.0);
    near(held.velocity_after_m_s[0], 0.0, 2e-14);
    near(held.friction_impulse_ns[0], 1.25, 2e-14);
    near(held.drive_impulse_ns[0], 0.75, 2e-14);
    near(held.displacement_m[0], -1.6, 2e-14);
    near(held.total_drive_work_j, -0.4, 2e-14);
    near(held.dissipated_j, 1.6, 2e-14);
    audit(&held);

    let rest = flow(&[0.0], &[3.0], &[1.0], &[0.75], 2.0);
    near(rest.friction_impulse_ns[0], -1.5, 2e-14);
    assert_eq!(rest.flow_velocity_after_m_s[0], 0.0);
    assert_eq!(rest.distance_m[0], 0.0);
    assert_eq!(rest.dissipated_j, 0.0);
    assert_eq!(rest.total_drive_work_j, 0.0);
    audit(&rest);

    // An explicit subsequent load segment starts again from the held state.
    let restart = flow(&rest.velocity_after_m_s[..1], &[3.0], &[0.5], &[-1.0], 2.0);
    near(restart.velocity_after_m_s[0], -3.0, 2e-14);
    near(restart.distance_m[0], 3.0, 2e-14);
    near(restart.dissipated_j, 1.5, 2e-14);
    audit(&restart);
}

#[test]
fn g1_coupled_arrest_retains_true_port_work_and_path_after_first_stop() {
    // W=[[2,1],[1,2]]. First phase a=[-3,-3], duration 1/3.
    // Then f=[1/2,-1] holds port zero while port one slows at -3/2
    // for another 4/3 s. Both stop at 5/3 s.
    let out = flow(
        &[1.0, 3.0],
        &[2.0, 1.0, 1.0, 2.0],
        &[1.0, 1.0],
        &[0.0, 0.0],
        3.0,
    );
    near(out.friction_impulse_ns[0], 1.0 / 3.0, 3e-14);
    near(out.friction_impulse_ns[1], -5.0 / 3.0, 3e-14);
    near(out.distance_m[0], 1.0 / 6.0, 3e-14);
    near(out.distance_m[1], 13.0 / 6.0, 3e-14);
    near(out.friction_work_j[0], -1.0 / 6.0, 3e-14);
    near(out.friction_work_j[1], -13.0 / 6.0, 3e-14);
    near(out.dissipated_j, 7.0 / 3.0, 3e-14);
    assert_eq!(out.flow_velocity_after_m_s[..2], [0.0, 0.0]);
    assert_eq!(out.stops, 2);
    assert_eq!(out.segments, 3);
    // Its integrated first impulse is positive, while its continuous
    // friction work is negative. Endpoint port pairing is not local heat.
    assert!(0.5 * out.friction_impulse_ns[0] * (1.0 + out.velocity_after_m_s[0]) > 0.0);
    audit(&out);
}

#[test]
fn g1_coupling_can_release_an_initially_resting_port_without_external_drive() {
    // At v=[0,3], holding port zero would require 1/2 N, exceeding 1/4 N.
    // It leaves with a=[-1/2,-7/4]. Port one stops at 12/7 s;
    // f=[1/4,-1/8] then holds it while port zero stops after another 16/7 s.
    let out = flow(
        &[0.0, 3.0],
        &[2.0, 1.0, 1.0, 2.0],
        &[0.25, 1.0],
        &[0.0, 0.0],
        5.0,
    );
    near(out.friction_impulse_ns[0], 1.0, 3e-14);
    near(out.friction_impulse_ns[1], -2.0, 3e-14);
    near(out.displacement_m[0], -12.0 / 7.0, 3e-14);
    near(out.displacement_m[1], 18.0 / 7.0, 3e-14);
    near(out.distance_m[0], 12.0 / 7.0, 3e-14);
    near(out.friction_work_j[0], -3.0 / 7.0, 3e-14);
    near(out.friction_work_j[1], -18.0 / 7.0, 3e-14);
    near(out.dissipated_j, 3.0, 3e-14);
    assert_eq!(out.flow_velocity_after_m_s[..2], [0.0, 0.0]);
    assert_eq!(out.stops, 2);
    audit(&out);
}

#[test]
fn g3_exact_path_is_consistent_under_time_partition_order_and_work_conjugate_units() {
    let w = [2.0, 1.0, 1.0, 2.0];
    let full = flow(&[1.0, 3.0], &w, &[1.0, 1.0], &[0.0, 0.0], 3.0);
    let first = flow(&[1.0, 3.0], &w, &[1.0, 1.0], &[0.0, 0.0], 0.2);
    let second = flow(
        &first.velocity_after_m_s[..2],
        &w,
        &[1.0, 1.0],
        &[0.0, 0.0],
        2.8,
    );
    for i in 0..2 {
        near(
            full.friction_impulse_ns[i],
            first.friction_impulse_ns[i] + second.friction_impulse_ns[i],
            5e-14,
        );
        near(
            full.displacement_m[i],
            first.displacement_m[i] + second.displacement_m[i],
            5e-14,
        );
        near(
            full.distance_m[i],
            first.distance_m[i] + second.distance_m[i],
            5e-14,
        );
        near(
            full.friction_work_j[i],
            first.friction_work_j[i] + second.friction_work_j[i],
            5e-14,
        );
    }
    let reversed = flow(&[3.0, 1.0], &w, &[1.0, 1.0], &[0.0, 0.0], 3.0);
    for i in 0..2 {
        near(
            reversed.friction_impulse_ns[1 - i],
            full.friction_impulse_ns[i],
            5e-14,
        );
        near(reversed.distance_m[1 - i], full.distance_m[i], 5e-14);
    }
    for scale in [[3.0, 0.25], [-1.0, 1.0]] {
        let matrix = [
            2.0 * scale[0] * scale[0],
            scale[0] * scale[1],
            scale[0] * scale[1],
            2.0 * scale[1] * scale[1],
        ];
        let out = flow(
            &[scale[0], 3.0 * scale[1]],
            &matrix,
            &[1.0 / scale[0].abs(), 1.0 / scale[1].abs()],
            &[0.0, 0.0],
            3.0,
        );
        for i in 0..2 {
            near(
                out.friction_impulse_ns[i] * scale[i],
                full.friction_impulse_ns[i],
                5e-14,
            );
            near(
                out.displacement_m[i] / scale[i],
                full.displacement_m[i],
                5e-14,
            );
            near(
                out.distance_m[i] / scale[i].abs(),
                full.distance_m[i],
                5e-14,
            );
            near(out.friction_work_j[i], full.friction_work_j[i], 5e-14);
        }
        near(out.dissipated_j, full.dissipated_j, 5e-14);
        audit(&out);
    }
}

#[test]
fn g1_empty_unforced_and_maximum_port_cases_preserve_declared_physics() {
    let empty = flow(&[], &[], &[], &[], 2.0);
    assert_eq!(empty.count, 0);
    assert_eq!(empty.segments, 0);
    assert_eq!(empty.factorizations, 0);
    assert_eq!(empty.dissipated_j, 0.0);

    let free = flow(&[2.0], &[0.5], &[0.0], &[4.0], 2.0);
    near(free.velocity_after_m_s[0], 6.0, 2e-14);
    near(free.distance_m[0], 8.0, 2e-14);
    near(free.total_drive_work_j, 32.0, 2e-14);
    assert_eq!(free.dissipated_j, 0.0);
    audit(&free);

    let n = MAX_FRICTION_IMPULSE_PORTS;
    let mut matrix = [0.0; MAX_FRICTION_IMPULSE_PORTS * MAX_FRICTION_IMPULSE_PORTS];
    let velocity = core::array::from_fn::<_, MAX_FRICTION_IMPULSE_PORTS, _>(|i| 2.0 + i as f64);
    for i in 0..n {
        matrix[i * n + i] = 1.0;
    }
    let out = flow(
        &velocity,
        &matrix,
        &[0.5; MAX_FRICTION_IMPULSE_PORTS],
        &[0.0; MAX_FRICTION_IMPULSE_PORTS],
        0.5,
    );
    for (i, &v) in velocity.iter().enumerate() {
        near(out.velocity_after_m_s[i], v - 0.25, 2e-14);
        near(out.friction_impulse_ns[i], -0.25, 2e-14);
        near(out.distance_m[i], 0.5 * v - 0.0625, 2e-14);
    }
    assert_eq!(out.factorizations, 1);
    assert_eq!(out.projection_iterations, 0);
    audit(&out);
}

#[test]
fn g4_budgets_refuse_without_replacing_physical_force_capacities() {
    let input = CoupledFrictionFlowInput {
        velocity_before_m_s: &[1.0, 3.0],
        mobility_per_kg: &[2.0, 1.0, 1.0, 2.0],
        capacity_n: &[1.0, 1.0],
        drive_n: &[0.0, 0.0],
        duration_s: 3.0,
    };
    let mut limited = controls();
    limited.maximum_segments = 1;
    assert!(matches!(
        advance_coupled_friction(&input, limited, || false),
        Err(CoupledFrictionFlowError::SegmentLimit { segments: 1, .. })
    ));
    limited = controls();
    limited.maximum_force_n = 0.5;
    assert!(matches!(
        advance_coupled_friction(&input, limited, || false),
        Err(CoupledFrictionFlowError::ForceBudget { .. })
    ));
    limited = controls();
    limited.projection.maximum_impulse_ns = 0.1;
    assert!(matches!(
        advance_coupled_friction(&input, limited, || false),
        Err(CoupledFrictionFlowError::ImpulseBudget { .. })
    ));

    // The physical set may be much wider than the numerical force budget.
    // A genuinely small held reaction must still succeed without shrinking it.
    let held = CoupledFrictionFlowInput {
        velocity_before_m_s: &[0.0],
        mobility_per_kg: &[1.0],
        capacity_n: &[1e6],
        drive_n: &[0.5],
        duration_s: 2.0,
    };
    limited = controls();
    limited.maximum_force_n = 1.0;
    let out = advance_coupled_friction(&held, limited, || false).unwrap();
    near(out.friction_impulse_ns[0], -1.0, 1e-14);
    assert_eq!(out.distance_m[0], 0.0);
    audit(&out);

    let singular = CoupledFrictionFlowInput {
        mobility_per_kg: &[1.0, 1.0, 1.0, 1.0],
        ..input
    };
    assert!(matches!(
        advance_coupled_friction(&singular, controls(), || false),
        Err(CoupledFrictionFlowError::Projection(
            CoupledFrictionImpulseError::Mobility { .. }
        ))
    ));
    let invalid = CoupledFrictionFlowInput {
        capacity_n: &[-1.0, 1.0],
        ..input
    };
    assert!(matches!(
        advance_coupled_friction(&invalid, controls(), || false),
        Err(CoupledFrictionFlowError::Invalid(_))
    ));
}

#[test]
fn g4_every_cancellation_boundary_returns_no_partial_path() {
    let velocity = [0.0, 3.0];
    let input = CoupledFrictionFlowInput {
        velocity_before_m_s: &velocity,
        mobility_per_kg: &[2.0, 1.0, 1.0, 2.0],
        capacity_n: &[0.25, 1.0],
        drive_n: &[0.0, 0.0],
        duration_s: 5.0,
    };
    let mut polls = 0;
    let expected = advance_coupled_friction(&input, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for at in 1..=polls {
        let mut seen = 0;
        assert!(
            matches!(
                advance_coupled_friction(&input, controls(), || {
                    seen += 1;
                    seen == at
                }),
                Err(CoupledFrictionFlowError::Cancelled)
            ),
            "poll {at}/{polls}"
        );
        assert_eq!(velocity, [0.0, 3.0]);
    }
    assert_eq!(
        expected,
        advance_coupled_friction(&input, controls(), || false).unwrap()
    );
}
