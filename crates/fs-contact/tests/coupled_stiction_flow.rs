use fs_contact::tangential::impact::multiple::relaxation::stiction::{
    CoupledStictionFlowControls, CoupledStictionFlowError, CoupledStictionFlowInput,
    CoupledStictionFlowReport, MAX_STICTION_FLOW_PORTS, MAX_STICTION_MODES_PER_PHASE,
    advance_coupled_stiction,
};
use fs_contact::tangential::impact::multiple::relaxation::{
    CoupledFrictionFlowControls, CoupledFrictionFlowError, CoupledFrictionFlowInput,
    advance_coupled_friction,
};
use fs_contact::tangential::impact::multiple::{
    CoupledFrictionImpulseControls, CoupledFrictionImpulseError,
};

fn controls() -> CoupledStictionFlowControls {
    CoupledStictionFlowControls {
        flow: CoupledFrictionFlowControls {
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
        },
        maximum_mode_evaluations: MAX_STICTION_MODES_PER_PHASE * 64,
    }
}

fn flow(
    v: &[f64],
    w: &[f64],
    cs: &[f64],
    ck: &[f64],
    d: &[f64],
    t: f64,
) -> CoupledStictionFlowReport {
    let out = advance_coupled_stiction(
        &CoupledStictionFlowInput {
            velocity_before_m_s: v,
            mobility_per_kg: w,
            static_capacity_n: cs,
            kinetic_capacity_n: ck,
            drive_n: d,
            duration_s: t,
        },
        controls(),
        || false,
    )
    .unwrap();
    audit(&out);
    out
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}

fn audit(out: &CoupledStictionFlowReport) {
    let path = &out.flow;
    near(path.energy_residual_j, 0.0, 2e-11);
    near(
        path.kinetic_energy_change_j,
        path.total_drive_work_j - path.dissipated_j,
        2e-11,
    );
    assert!(path.dissipated_j >= 0.0);
    assert!(out.diagnostics.mode_evaluations <= MAX_STICTION_MODES_PER_PHASE * path.segments);
    assert!(out.diagnostics.maximum_holding_acceleration_residual_m_s2 <= 2e-11);
    for i in 0..path.count {
        assert!(path.friction_work_j[i] <= 0.0);
        assert!(path.distance_m[i] >= path.displacement_m[i].abs());
        assert!(path.minimum_speed_m_s[i] >= 0.0);
        assert!(path.minimum_speed_m_s[i] <= path.flow_velocity_after_m_s[i].abs());
        assert!(path.maximum_speed_m_s[i] >= path.flow_velocity_after_m_s[i].abs());
        near(path.velocity_residual_m_s[i], 0.0, 2e-11);
        near(path.impulse_residual_ns[i], 0.0, 2e-11);
        assert!(path.impulse_capacity_residual_ns[i] <= 2e-11);
    }
}

#[test]
fn static_threshold_holds_but_identical_load_on_a_moving_port_uses_kinetic_force() {
    let held = flow(&[0.0], &[1.0], &[2.0], &[1.0], &[1.5], 2.0);
    assert_eq!(held.flow.flow_velocity_after_m_s[0], 0.0);
    assert_eq!(held.flow.friction_impulse_ns[0], -3.0);
    assert_eq!(held.flow.dissipated_j, 0.0);
    assert_eq!(held.flow.distance_m[0], 0.0);
    assert_eq!(held.diagnostics.last_phase_holding_mask, 1);

    // There is no capture band: an arbitrarily small positive slip already
    // uses Ck, so this same sub-breakaway load now accelerates the port.
    let moving = flow(&[1e-15], &[1.0], &[2.0], &[1.0], &[1.5], 2.0);
    near(moving.flow.velocity_after_m_s[0], 1.0 + 1e-15, 1e-15);
    near(moving.flow.distance_m[0], 1.0 + 2e-15, 1e-15);
    assert_eq!(moving.flow.friction_impulse_ns[0], -2.0);
    assert_eq!(moving.flow.minimum_speed_m_s[0], 1e-15);
    assert_eq!(moving.flow.stops, 0);
    assert_eq!(moving.diagnostics.last_phase_holding_mask, 0);
    assert!(moving.flow.dissipated_j > 0.0);

    let threshold = flow(&[0.0], &[3.0], &[2.0], &[1.0], &[2.0], 1.0);
    assert_eq!(threshold.flow.flow_velocity_after_m_s[0], 0.0);
    assert_eq!(threshold.flow.friction_impulse_ns[0], -2.0);
    assert_eq!(threshold.diagnostics.breakaways, 0);
    let above = f64::from_bits(2.0_f64.to_bits() + 1);
    let release = flow(&[0.0], &[3.0], &[2.0], &[1.0], &[above], 1.0);
    near(
        release.flow.velocity_after_m_s[0],
        3.0 * (above - 1.0),
        2e-14,
    );
    assert_eq!(release.flow.friction_impulse_ns[0], -1.0);
    assert_eq!(release.diagnostics.breakaways, 1);
}

#[test]
fn arrest_recaptures_below_static_limit_and_reverses_above_it() {
    let caught = flow(&[-1.0], &[1.0], &[2.0], &[1.0], &[1.5], 1.0);
    assert_eq!(caught.flow.flow_velocity_after_m_s[0], 0.0);
    near(caught.flow.friction_impulse_ns[0], -0.5, 2e-14);
    near(caught.flow.distance_m[0], 0.2, 2e-14);
    near(caught.flow.displacement_m[0], -0.2, 2e-14);
    near(caught.flow.dissipated_j, 0.2, 2e-14);
    near(caught.flow.total_drive_work_j, -0.3, 2e-14);
    assert_eq!(caught.flow.stops, 1);
    assert_eq!(caught.flow.segments, 2);
    assert_eq!(caught.diagnostics.breakaways, 0);
    assert_eq!(caught.diagnostics.last_phase_holding_mask, 1);

    let reversed = flow(&[-1.0], &[1.0], &[2.0], &[1.0], &[3.0], 1.0);
    near(reversed.flow.velocity_after_m_s[0], 1.5, 2e-14);
    near(reversed.flow.friction_impulse_ns[0], -0.5, 2e-14);
    near(reversed.flow.displacement_m[0], 7.0 / 16.0, 2e-14);
    near(reversed.flow.distance_m[0], 11.0 / 16.0, 2e-14);
    near(reversed.flow.dissipated_j, 11.0 / 16.0, 2e-14);
    near(reversed.flow.total_drive_work_j, 21.0 / 16.0, 2e-14);
    assert_eq!(reversed.flow.minimum_speed_m_s[0], 0.0);
    near(reversed.flow.maximum_speed_m_s[0], 1.5, 2e-14);
    assert_eq!(reversed.diagnostics.breakaways, 1);
    assert_eq!(reversed.diagnostics.last_phase_holding_mask, 0);
    assert_eq!(reversed.diagnostics.endpoint_resting_mask, 0);

    // An event at the requested endpoint is reported without inventing a
    // subsequent release phase or calling its not-yet-needed mode selection.
    let terminal_stop = flow(&[-1.0], &[1.0], &[2.0], &[1.0], &[3.0], 0.25);
    assert_eq!(terminal_stop.flow.flow_velocity_after_m_s[0], 0.0);
    assert_eq!(terminal_stop.flow.segments, 1);
    assert_eq!(terminal_stop.diagnostics.last_phase_holding_mask, 0);
    assert_eq!(terminal_stop.diagnostics.endpoint_resting_mask, 1);
    assert_eq!(terminal_stop.diagnostics.breakaways, 0);
}

#[test]
fn shared_mobility_supplies_static_reaction_to_another_ports_sliding_force() {
    let w = [2.0, 1.0, 1.0, 2.0];
    // f=[1/2,-1] gives a=[0,-3/2]. Independent diagonal sliders cannot
    // reproduce either the nonzero static reaction or the moving acceleration.
    let out = flow(&[0.0, 2.0], &w, &[2.0, 2.0], &[0.5, 1.0], &[0.0, 0.0], 1.0);
    near(out.flow.velocity_after_m_s[0], 0.0, 2e-14);
    near(out.flow.velocity_after_m_s[1], 0.5, 2e-14);
    near(out.flow.friction_impulse_ns[0], 0.5, 2e-14);
    near(out.flow.friction_impulse_ns[1], -1.0, 2e-14);
    assert_eq!(out.flow.distance_m[0], 0.0);
    assert_eq!(out.flow.friction_work_j[0], 0.0);
    near(out.flow.distance_m[1], 1.25, 2e-14);
    near(out.flow.dissipated_j, 1.25, 2e-14);
    assert_eq!(out.diagnostics.last_phase_holding_mask, 1);
    assert!(out.diagnostics.mode_factorizations > 0);

    let stopped = flow(&[0.0, 2.0], &w, &[2.0, 2.0], &[0.5, 1.0], &[0.0, 0.0], 2.0);
    assert_eq!(stopped.flow.flow_velocity_after_m_s[..2], [0.0, 0.0]);
    near(stopped.flow.dissipated_j, 4.0 / 3.0, 3e-14);
    assert_eq!(stopped.diagnostics.last_phase_holding_mask, 3);
    assert_eq!(stopped.diagnostics.breakaways, 0);
}

#[test]
fn coupled_breakaway_can_preserve_a_different_resting_coordinate() {
    // Full hold needs f=[-2,0] and is impossible. Releasing port zero
    // leaves f=[-1/2,-3/4], a=[9/4,0], and port one still at rest.
    let out = flow(
        &[0.0, 0.0],
        &[2.0, 1.0, 1.0, 2.0],
        &[1.0, 1.0],
        &[0.5, 0.5],
        &[2.0, 0.0],
        1.0,
    );
    near(out.flow.velocity_after_m_s[0], 2.25, 2e-14);
    near(out.flow.velocity_after_m_s[1], 0.0, 2e-14);
    near(out.flow.friction_impulse_ns[0], -0.5, 2e-14);
    near(out.flow.friction_impulse_ns[1], -0.75, 2e-14);
    near(out.flow.displacement_m[0], 1.125, 2e-14);
    near(out.flow.dissipated_j, 0.5625, 2e-14);
    near(out.flow.total_drive_work_j, 2.25, 2e-14);
    assert_eq!(out.diagnostics.last_phase_holding_mask, 2);
    assert_eq!(out.diagnostics.breakaways, 1);
}

#[test]
fn whole_path_speed_extrema_and_physics_survive_permutation_and_conjugate_units() {
    // Coupling releases port zero at a=-1/2. Port one stops at 4/7 s;
    // port zero then returns to rest at 4/3 s. Both outer velocities of
    // port zero are zero, but its actual maximum speed is 2/7 m/s.
    let w = [2.0, 1.0, 1.0, 2.0];
    let out = flow(&[0.0, 1.0], &w, &[0.3, 1.2], &[0.25, 1.0], &[0.0, 0.0], 2.0);
    assert_eq!(out.flow.flow_velocity_after_m_s[..2], [0.0, 0.0]);
    assert_eq!(out.flow.minimum_speed_m_s[..2], [0.0, 0.0]);
    near(out.flow.maximum_speed_m_s[0], 2.0 / 7.0, 3e-14);
    near(out.flow.maximum_speed_m_s[1], 1.0, 3e-14);
    assert_eq!(out.flow.stops, 2);
    let permuted = flow(&[1.0, 0.0], &w, &[1.2, 0.3], &[1.0, 0.25], &[0.0, 0.0], 2.0);
    for i in 0..2 {
        near(
            permuted.flow.friction_impulse_ns[1 - i],
            out.flow.friction_impulse_ns[i],
            5e-14,
        );
        near(
            permuted.flow.distance_m[1 - i],
            out.flow.distance_m[i],
            5e-14,
        );
        near(
            permuted.flow.maximum_speed_m_s[1 - i],
            out.flow.maximum_speed_m_s[i],
            5e-14,
        );
    }
    for scale in [[4.0_f64, 0.25], [-1.0, 1.0]] {
        let matrix = [
            2.0 * scale[0] * scale[0],
            scale[0] * scale[1],
            scale[0] * scale[1],
            2.0 * scale[1] * scale[1],
        ];
        let transformed = flow(
            &[0.0, scale[1]],
            &matrix,
            &[0.3 / scale[0].abs(), 1.2 / scale[1].abs()],
            &[0.25 / scale[0].abs(), 1.0 / scale[1].abs()],
            &[0.0, 0.0],
            2.0,
        );
        near(transformed.flow.dissipated_j, out.flow.dissipated_j, 5e-14);
        for i in 0..2 {
            near(
                transformed.flow.friction_impulse_ns[i] * scale[i],
                out.flow.friction_impulse_ns[i],
                5e-14,
            );
            near(
                transformed.flow.displacement_m[i] / scale[i],
                out.flow.displacement_m[i],
                5e-14,
            );
            near(
                transformed.flow.maximum_speed_m_s[i] / scale[i].abs(),
                out.flow.maximum_speed_m_s[i],
                5e-14,
            );
        }
    }
}

#[test]
fn incomparable_static_continuations_refuse_in_both_port_orders() {
    // Exact binary data has two maximal holding sets: hold zero with
    // f=[-7/8,1/2], or hold one with f=[-1/2,7/8]. Both-release is
    // dominated. Neither index order nor a warm start may choose a branch.
    let velocity = [0.0, 0.0];
    let w = [1.0, 0.5, 0.5, 1.0];
    for drive in [[1.25, -1.25], [-1.25, 1.25]] {
        let error = advance_coupled_stiction(
            &CoupledStictionFlowInput {
                velocity_before_m_s: &velocity,
                mobility_per_kg: &w,
                static_capacity_n: &[1.0, 1.0],
                kinetic_capacity_n: &[0.5, 0.5],
                drive_n: &drive,
                duration_s: 1.0,
            },
            controls(),
            || false,
        )
        .unwrap_err();
        match error {
            CoupledStictionFlowError::Ambiguous {
                resting_mask: 3,
                maximal_modes: 2,
                first_holding_mask,
                second_holding_mask,
            } => {
                assert_ne!(first_holding_mask, second_holding_mask);
                assert_eq!(first_holding_mask | second_holding_mask, 3);
                assert_eq!(first_holding_mask & second_holding_mask, 0);
            }
            other => panic!("expected unresolved physical branching, got {other:?}"),
        }
        assert_eq!(velocity, [0.0, 0.0]);
        assert_eq!(w, [1.0, 0.5, 0.5, 1.0]);
    }
}

#[test]
fn equal_static_and_kinetic_limits_delegate_to_the_original_complete_receipt() {
    let velocity = [1.0, 3.0];
    let w = [2.0, 1.0, 1.0, 2.0];
    let capacity = [1.0, 1.0];
    let drive = [0.0, 0.0];
    let legacy = advance_coupled_friction(
        &CoupledFrictionFlowInput {
            velocity_before_m_s: &velocity,
            mobility_per_kg: &w,
            capacity_n: &capacity,
            drive_n: &drive,
            duration_s: 3.0,
        },
        controls().flow,
        || false,
    )
    .unwrap();
    let out = flow(&velocity, &w, &capacity, &capacity, &drive, 3.0);
    assert_eq!(out.flow, legacy);
    assert!(out.diagnostics.used_equal_bound_solver);
    assert_eq!(out.diagnostics.mode_evaluations, 0);
    assert_eq!(out.diagnostics.endpoint_resting_mask, 3);
    assert_eq!(out.flow.minimum_speed_m_s[..2], [0.0, 0.0]);
    assert_eq!(out.flow.maximum_speed_m_s[..2], [1.0, 3.0]);
}

#[test]
fn all_six_port_modes_are_bounded_and_exhaustion_or_cancellation_returns_no_path() {
    let mut w = [0.0; MAX_STICTION_FLOW_PORTS * MAX_STICTION_FLOW_PORTS];
    for i in 0..MAX_STICTION_FLOW_PORTS {
        w[i * MAX_STICTION_FLOW_PORTS + i] = 1.0;
    }
    let velocity = [0.0; MAX_STICTION_FLOW_PORTS];
    let drive = [2.0; MAX_STICTION_FLOW_PORTS];
    let input = CoupledStictionFlowInput {
        velocity_before_m_s: &velocity,
        mobility_per_kg: &w,
        static_capacity_n: &[1.0; MAX_STICTION_FLOW_PORTS],
        kinetic_capacity_n: &[0.5; MAX_STICTION_FLOW_PORTS],
        drive_n: &drive,
        duration_s: 1.0,
    };
    let mut bounded = controls();
    bounded.maximum_mode_evaluations = MAX_STICTION_MODES_PER_PHASE - 1;
    assert_eq!(
        advance_coupled_stiction(&input, bounded, || false).unwrap_err(),
        CoupledStictionFlowError::ModeBudget {
            evaluations: MAX_STICTION_MODES_PER_PHASE - 1,
            maximum: MAX_STICTION_MODES_PER_PHASE - 1,
        }
    );
    bounded.maximum_mode_evaluations = MAX_STICTION_MODES_PER_PHASE;
    let resolved = advance_coupled_stiction(&input, bounded, || false).unwrap();
    audit(&resolved);
    assert_eq!(
        resolved.diagnostics.mode_evaluations,
        MAX_STICTION_MODES_PER_PHASE
    );
    assert_eq!(resolved.diagnostics.breakaways, MAX_STICTION_FLOW_PORTS);
    assert_eq!(resolved.flow.flow_velocity_after_m_s[..6], [1.5; 6]);
    near(resolved.flow.dissipated_j, 2.25, 3e-14);

    let mut polls = 0;
    assert_eq!(
        advance_coupled_stiction(&input, bounded, || {
            polls += 1;
            polls > 75
        })
        .unwrap_err(),
        CoupledStictionFlowError::Flow(CoupledFrictionFlowError::Cancelled)
    );
    assert!(polls > 75);
    assert_eq!(velocity, [0.0; MAX_STICTION_FLOW_PORTS]);
    assert_eq!(drive, [2.0; MAX_STICTION_FLOW_PORTS]);
}

#[test]
fn identity_paths_and_original_matrix_and_capacity_admission_remain_explicit() {
    let identity = flow(&[-2.0], &[3.0], &[2.0], &[1.0], &[1.5], 0.0);
    assert_eq!(identity.flow.velocity_after_m_s[0], -2.0);
    assert_eq!(identity.flow.minimum_speed_m_s[0], 2.0);
    assert_eq!(identity.flow.maximum_speed_m_s[0], 2.0);
    assert_eq!(identity.flow.dissipated_j, 0.0);
    assert_eq!(identity.diagnostics.mode_evaluations, 0);
    let empty = flow(&[], &[], &[], &[], &[], 1.0);
    assert_eq!(empty.flow.count, 0);
    assert_eq!(empty.flow.segments, 0);

    let input = CoupledStictionFlowInput {
        velocity_before_m_s: &[0.0],
        mobility_per_kg: &[1.0],
        static_capacity_n: &[1.0],
        kinetic_capacity_n: &[2.0],
        drive_n: &[0.0],
        duration_s: 1.0,
    };
    assert!(matches!(
        advance_coupled_stiction(&input, controls(), || false),
        Err(CoupledStictionFlowError::Invalid(_))
    ));
    let indefinite = CoupledStictionFlowInput {
        velocity_before_m_s: &[0.0, 0.0],
        mobility_per_kg: &[1.0, 2.0, 2.0, 1.0],
        static_capacity_n: &[2.0, 2.0],
        kinetic_capacity_n: &[1.0, 1.0],
        drive_n: &[0.0, 0.0],
        duration_s: 1.0,
    };
    assert!(matches!(
        advance_coupled_stiction(&indefinite, controls(), || false),
        Err(CoupledStictionFlowError::Flow(
            CoupledFrictionFlowError::Projection(CoupledFrictionImpulseError::Mobility { .. })
        ))
    ));
    let too_many = CoupledStictionFlowInput {
        velocity_before_m_s: &[0.0; 7],
        mobility_per_kg: &[0.0; 49],
        static_capacity_n: &[2.0; 7],
        kinetic_capacity_n: &[1.0; 7],
        drive_n: &[0.0; 7],
        duration_s: 1.0,
    };
    assert!(matches!(
        advance_coupled_stiction(&too_many, controls(), || false),
        Err(CoupledStictionFlowError::Invalid(_))
    ));
}
