//! Independent G1 trajectories and G3 phase/partition checks of the full pencil.
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::*;
use fs_exec::CancelGate;

fn config() -> LinearNetworkEventConfig {
    LinearNetworkEventConfig {
        maximum_modes: 16,
        maximum_connections: 16,
        maximum_workspace_bytes: 1_000_000,
        maximum_cells: 8192,
        maximum_events: 64,
        maximum_operator_evaluations: 1_000_000,
        maximum_root_evaluations: 1_000_000,
        maximum_root_iterations: 80,
        maximum_step_s: 10.0,
        event_time_tolerance_s: 2e-14,
        sticking_velocity_tolerance_m_s: 1e-11,
        maximum_abs_modal_coordinate: 100.0,
        maximum_abs_modal_velocity: 100.0,
        energy_absolute_tolerance_j: 1e-11,
        energy_relative_tolerance: 1e-9,
        maximum_numerical_work_uncertainty_j: 1e-7,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "actual={a:.17e}, reference={b:.17e}, difference={:.4e}",
        a - b
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn link(column: Vec<f64>, k: f64, c: f64, rest: f64) -> LinearNetworkLink {
    LinearNetworkLink {
        column,
        stiffness_n_per_m: k,
        damping_n_s_per_m: c,
        rest_extension_m: rest,
    }
}

#[test]
fn complete_chain_modes_and_original_rest_forces_match_independent_spectrum() {
    let links = vec![
        link(vec![1.0, -1.0, 0.0], 3.0, 0.0, 0.02),
        link(vec![0.0, 1.0, -1.0], 3.0, 0.0, 0.03),
    ];
    let mut flow = LinearNetworkEventFlow::new(
        vec![1.0; 3],
        vec![0.0; 3],
        links,
        vec![1.0, 0.5, -0.25],
        config(),
    )
    .unwrap();
    let q0 = [0.2, -0.1, 0.15];
    let v0 = [0.3, 0.2, -0.4];
    let force = [0.02, -0.03, 0.01];
    let dt = 2.3;
    let s = flow
        .advance(
            &q0,
            &v0,
            &force,
            LinearNetworkPhase::Fresh,
            dt,
            0.0,
            0.0,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let forcing = [force[0] + 0.06, force[1] - 0.06 + 0.09, force[2] - 0.09];
    let basis: [([f64; 3], f64); 3] = [
        ([1.0 / 3.0_f64.sqrt(); 3], 1.0),
        ([1.0 / 2.0_f64.sqrt(), 0.0, -1.0 / 2.0_f64.sqrt()], 4.0),
        (
            [
                1.0 / 6.0_f64.sqrt(),
                -2.0 / 6.0_f64.sqrt(),
                1.0 / 6.0_f64.sqrt(),
            ],
            10.0,
        ),
    ];
    let mut q = [0.0; 3];
    let mut v = [0.0; 3];
    for (axis, lambda) in basis {
        let w = lambda.sqrt();
        let x = dot(&axis, &q0);
        let u = dot(&axis, &v0);
        let equilibrium = dot(&axis, &forcing) / lambda;
        let next = equilibrium + (x - equilibrium) * (w * dt).cos() + u / w * (w * dt).sin();
        let speed = -w * (x - equilibrium) * (w * dt).sin() + u * (w * dt).cos();
        for i in 0..3 {
            q[i] += axis[i] * next;
            v[i] += axis[i] * speed;
        }
    }
    for i in 0..3 {
        close(flow.positions()[i], q[i], 3e-14);
        close(flow.velocities()[i], v[i], 3e-14);
    }
    assert!(s.numerical_work_uncertainty_j > 0.0);
    assert!(s.energy_residual_j.abs() < 2e-14);
    assert_eq!(s.modal_damping_dissipation_j, 0.0);
    assert_eq!(s.connection_damping_dissipation_j, 0.0);
    assert_eq!(s.friction_dissipation_j, 0.0);
    close(
        s.external_work_j,
        dot(&force, &[q[0] - q0[0], q[1] - q0[1], q[2] - q0[2]]),
        2e-15,
    );
}

fn rhs(x: [f64; 9], force: [f64; 3], links: &[LinearNetworkLink]) -> [f64; 9] {
    let mut r = [0.0; 9];
    let k = [1.0, 4.0, 9.0];
    let c = [0.1, 0.2, 0.4];
    for i in 0..3 {
        r[i] = x[3 + i];
        r[3 + i] = force[i] - k[i] * x[i] - c[i] * x[3 + i];
        r[6] += c[i] * x[3 + i] * x[3 + i];
    }
    for (j, l) in links.iter().enumerate() {
        let z = dot(&l.column, &x[..3]) - l.rest_extension_m;
        let v = dot(&l.column, &x[3..6]);
        let f = l.stiffness_n_per_m * z + l.damping_n_s_per_m * v;
        for i in 0..3 {
            r[3 + i] -= l.column[i] * f;
        }
        r[7 + j] = l.damping_n_s_per_m * v * v;
    }
    r
}
fn add(x: [f64; 9], r: [f64; 9], h: f64) -> [f64; 9] {
    core::array::from_fn(|i| x[i] + h * r[i])
}

#[test]
fn noncommuting_modal_and_link_damping_match_independent_runge_kutta_work() {
    let links = vec![
        link(vec![1.0, -1.0, 0.2], 2.0, 0.15, 0.03),
        link(vec![0.1, 1.0, -1.0], 0.7, 0.2, -0.02),
    ];
    let mut flow = LinearNetworkEventFlow::new(
        vec![1.0, 4.0, 9.0],
        vec![0.1, 0.2, 0.4],
        links.clone(),
        vec![1.0, -0.5, 0.3],
        config(),
    )
    .unwrap();
    let q0 = [0.1, -0.2, 0.05];
    let v0 = [0.7, -0.1, 0.3];
    let force = [0.1, -0.05, 0.02];
    let dt = 0.7;
    let actual = flow
        .advance(
            &q0,
            &v0,
            &force,
            LinearNetworkPhase::Fresh,
            dt,
            0.0,
            0.0,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut x = [0.1, -0.2, 0.05, 0.7, -0.1, 0.3, 0.0, 0.0, 0.0];
    let h = dt / 8192.0;
    for _ in 0..8192 {
        let a = rhs(x, force, &links);
        let b = rhs(add(x, a, h / 2.0), force, &links);
        let c = rhs(add(x, b, h / 2.0), force, &links);
        let d = rhs(add(x, c, h), force, &links);
        for i in 0..9 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    for i in 0..3 {
        close(flow.positions()[i], x[i], 2e-13);
        close(flow.velocities()[i], x[3 + i], 2e-13);
    }
    close(actual.modal_damping_dissipation_j, x[6], 2e-13);
    for i in 0..2 {
        close(flow.connection_dissipation_j()[i], x[7 + i], 2e-13);
        assert!(flow.connection_dissipation_j()[i] > 0.0);
        assert!(flow.connection_force_bounds_n()[i] >= flow.connection_forces_n()[i].abs());
    }
    close(actual.connection_damping_dissipation_j, x[7] + x[8], 2e-13);
    assert!(actual.energy_residual_j.abs() < 2e-14);
}

#[test]
fn elastic_common_motion_changes_static_reaction_and_resolves_release() {
    let mut flow = LinearNetworkEventFlow::new(
        vec![1.0, 4.0],
        vec![0.0; 2],
        vec![],
        vec![1.0, -1.0],
        config(),
    )
    .unwrap();
    let q0 = [0.2; 2];
    let v0 = [2.0; 2];
    let dt = 0.6;
    let w = 2.5_f64.sqrt();
    let radius = 0.2_f64.hypot(2.0 / w);
    let phase = 0.2_f64.atan2(2.0 / w);
    let release = ((2.0 / 3.0 / radius).asin() - phase) / w;
    let x = 0.2 * (w * release).cos() + 2.0 / w * (w * release).sin();
    let v = -0.2 * w * (w * release).sin() + 2.0 * (w * release).cos();
    let s = flow
        .advance(
            &q0,
            &v0,
            &[0.0; 2],
            LinearNetworkPhase::Fresh,
            dt,
            1.0,
            0.6,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(
        s.events[0].unwrap().kind,
        LinearNetworkEventKind::SlipStarted { direction: 1 }
    );
    close(s.events[0].unwrap().time_s, release, 5e-13);
    let t = dt - release;
    for (i, (omega, force)) in [(1.0_f64, -0.6), (2.0, 0.6)].into_iter().enumerate() {
        let equilibrium = force / (omega * omega);
        let q = equilibrium + (x - equilibrium) * (omega * t).cos() + v / omega * (omega * t).sin();
        let velocity = -omega * (x - equilibrium) * (omega * t).sin() + v * (omega * t).cos();
        close(flow.positions()[i], q, 2e-12);
        close(flow.velocities()[i], velocity, 2e-12);
    }
    assert!(s.friction_dissipation_j > 0.0);
    assert!(s.energy_residual_j.abs() < 2e-13);
}

#[test]
fn stopping_static_hold_restart_and_arbitrarily_small_initial_motion() {
    let create = || {
        LinearNetworkEventFlow::new(
            vec![0.0; 2],
            vec![0.0; 2],
            vec![],
            vec![1.0, -1.0],
            config(),
        )
        .unwrap()
    };
    let mut flow = create();
    let s = flow
        .advance(
            &[0.0; 2],
            &[0.25, -0.25],
            &[0.0; 2],
            LinearNetworkPhase::Fresh,
            1.0,
            1.0,
            0.5,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(s.events[0].unwrap().time_s, 0.5, 2e-14);
    assert_eq!(s.phase, LinearNetworkPhase::Sticking);
    assert_eq!(flow.velocities(), &[0.0; 2]);
    close(flow.positions()[0], 0.0625, 2e-15);
    close(flow.positions()[1], -0.0625, 2e-15);
    close(s.friction_dissipation_j, 0.0625, 2e-15);
    let q = flow.positions().to_vec();
    let v = flow.velocities().to_vec();
    let restart = flow
        .advance(
            &q,
            &v,
            &[2.0, -2.0],
            s.phase,
            0.2,
            1.0,
            0.5,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(
        restart.events[0].unwrap().kind,
        LinearNetworkEventKind::SlipStarted { direction: 1 }
    );
    close(flow.velocities()[0], 0.3, 2e-15);
    close(flow.positions()[0], 0.0925, 2e-15);
    close(restart.external_work_j, 0.12, 2e-15);
    close(restart.friction_dissipation_j, 0.03, 2e-15);
    let mut tiny = create();
    let s = tiny
        .advance(
            &[0.0; 2],
            &[1e-13, -1e-13],
            &[0.0; 2],
            LinearNetworkPhase::Fresh,
            1e-14,
            1.0,
            0.5,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(s.event_count, 0);
    assert!(tiny.velocities()[0] > 0.0);
    assert!(s.accumulated_slip_m > 0.0);
    // If a force command changes exactly at zero-velocity release, the new
    // subthreshold load can be held. Occupancy alone must not launch motion.
    let mut changed = create();
    let held = changed
        .advance(
            &[0.0; 2],
            &[0.0; 2],
            &[0.75, -0.75],
            LinearNetworkPhase::Sliding { direction: 1 },
            0.1,
            1.0,
            0.5,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(held.phase, LinearNetworkPhase::Sticking);
    assert_eq!(changed.velocities(), &[0.0; 2]);
    assert_eq!(held.friction_dissipation_j, 0.0);
}

#[test]
fn exactly_stationary_high_frequency_network_retains_static_impulse_in_one_cell() {
    let mut budget = config();
    budget.maximum_cells = 1;
    let mut flow =
        LinearNetworkEventFlow::new(vec![1e8; 2], vec![1.0; 2], vec![], vec![1.0, -1.0], budget)
            .unwrap();
    let q = [1e-8, -1e-8];
    let s = flow
        .advance(
            &q,
            &[0.0; 2],
            &[0.0; 2],
            LinearNetworkPhase::Fresh,
            1.0,
            1.1,
            0.5,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(flow.positions(), &q);
    assert_eq!(flow.velocities(), &[0.0; 2]);
    assert_eq!(s.phase, LinearNetworkPhase::Sticking);
    assert_eq!(s.cells, 1);
    assert_eq!(s.event_count, 0);
    close(s.contact_impulse_n_s, -1.0, 2e-15);
    close(s.contact_time_moment_n_s2, -0.5, 2e-15);
    assert_eq!(s.friction_dissipation_j, 0.0);
    assert_eq!(s.modal_damping_dissipation_j, 0.0);
}

#[test]
fn two_internal_stops_are_found_when_cell_endpoint_velocities_have_the_same_sign() {
    // v_contact = cos(t-centre)-cos(half_width) has two roots in this ONE
    // admitted Taylor cell; checking just its endpoints misses both.
    let centre = 0.06_f64;
    let half_width = 0.01_f64;
    let q0 = [0.0, -centre.sin()];
    let v0 = [-half_width.cos(), centre.cos()];
    let mut flow =
        LinearNetworkEventFlow::new(vec![0.0, 1.0], vec![0.0; 2], vec![], vec![1.0; 2], config())
            .unwrap();
    let s = flow
        .advance(
            &q0,
            &v0,
            &[0.0; 2],
            LinearNetworkPhase::Fresh,
            0.12,
            1e-4,
            0.0,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let stops: Vec<_> = s
        .events
        .iter()
        .flatten()
        .filter(|e| e.kind == LinearNetworkEventKind::Stopped)
        .collect();
    assert_eq!(stops.len(), 2);
    close(stops[0].time_s, 0.05, 1e-12);
    close(stops[1].time_s, 0.07, 1e-12);
    close(flow.positions()[0], v0[0] * 0.12, 2e-13);
    close(flow.positions()[1], (0.12 - centre).sin(), 2e-13);
    assert_eq!(s.friction_dissipation_j, 0.0);
    assert!(s.numerical_work_uncertainty_j > 0.0);
}

#[test]
fn multiple_modes_events_and_losses_are_stable_under_time_partition() {
    let create = || {
        LinearNetworkEventFlow::new(
            vec![1.0, 9.0, 25.0],
            vec![0.01, 0.02, 0.03],
            vec![link(vec![1.0, -1.0, 0.0], 0.4, 0.02, 0.01)],
            vec![1.0, -0.6, 0.25],
            config(),
        )
        .unwrap()
    };
    let q0 = [0.1, -0.05, 0.02];
    let v0 = [0.3, -0.2, 0.1];
    let force = [0.02, -0.01, 0.015];
    let dt = 2.0;
    let mut whole = create();
    let reference = whole
        .advance(
            &q0,
            &v0,
            &force,
            LinearNetworkPhase::Fresh,
            dt,
            0.02,
            0.01,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(reference.event_count >= 2);
    for parts in [7, 23, 101] {
        let mut flow = create();
        let mut q = q0.to_vec();
        let mut v = v0.to_vec();
        let mut phase = LinearNetworkPhase::Fresh;
        let mut friction = 0.0;
        let mut modal = 0.0;
        let mut links = 0.0;
        let mut impulse = 0.0;
        let mut moment = 0.0;
        for j in 0..parts {
            let s = flow
                .advance(
                    &q,
                    &v,
                    &force,
                    phase,
                    dt / parts as f64,
                    0.02,
                    0.01,
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            q.copy_from_slice(flow.positions());
            v.copy_from_slice(flow.velocities());
            phase = s.phase;
            friction += s.friction_dissipation_j;
            modal += s.modal_damping_dissipation_j;
            links += s.connection_damping_dissipation_j;
            impulse += s.contact_impulse_n_s;
            moment +=
                s.contact_time_moment_n_s2 + (j as f64 * dt / parts as f64) * s.contact_impulse_n_s;
        }
        for i in 0..3 {
            close(q[i], whole.positions()[i], 3e-11);
            close(v[i], whole.velocities()[i], 3e-11);
        }
        close(friction, reference.friction_dissipation_j, 3e-12);
        close(modal, reference.modal_damping_dissipation_j, 3e-12);
        close(links, reference.connection_damping_dissipation_j, 3e-12);
        close(impulse, reference.contact_impulse_n_s, 3e-12);
        close(moment, reference.contact_time_moment_n_s2, 3e-12);
    }
}

#[test]
fn cancellation_work_limits_and_late_accuracy_refusal_preserve_identical_retry() {
    let create = |c| {
        LinearNetworkEventFlow::new(vec![1.0, 4.0], vec![0.1, 0.2], vec![], vec![1.0, -1.0], c)
            .unwrap()
    };
    let mut flow = create(config());
    let mut independent = create(config());
    let q = [0.1, 0.2];
    let v = [0.3, -0.1];
    let f = [0.0; 2];
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert_eq!(
        flow.advance(
            &q,
            &v,
            &f,
            LinearNetworkPhase::Fresh,
            0.2,
            0.1,
            0.05,
            &cancel
        ),
        Err(LinearNetworkEventError::Cancelled)
    );
    assert!(
        flow.advance(
            &q,
            &v,
            &[1e8, -1e8],
            LinearNetworkPhase::Fresh,
            0.2,
            0.1,
            0.05,
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let a = flow
        .advance(
            &q,
            &v,
            &f,
            LinearNetworkPhase::Fresh,
            0.2,
            0.1,
            0.05,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let b = independent
        .advance(
            &q,
            &v,
            &f,
            LinearNetworkPhase::Fresh,
            0.2,
            0.1,
            0.05,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(flow.positions(), independent.positions());
    assert_eq!(flow.velocities(), independent.velocities());
    let mut strict = config();
    strict.maximum_numerical_work_uncertainty_j = a.numerical_work_uncertainty_j / 2.0;
    assert!(matches!(
        create(strict).advance(
            &q,
            &v,
            &f,
            LinearNetworkPhase::Fresh,
            0.2,
            0.1,
            0.05,
            &CancelGate::new_clock_free()
        ),
        Err(LinearNetworkEventError::EnergyAccuracy { .. })
    ));
    let mut exhausted = config();
    exhausted.maximum_root_evaluations = 1;
    assert!(matches!(
        create(exhausted).advance(
            &q,
            &v,
            &f,
            LinearNetworkPhase::Fresh,
            0.2,
            0.1,
            0.05,
            &CancelGate::new_clock_free()
        ),
        Err(LinearNetworkEventError::Budget("dense event evaluations"))
    ));
    let mut workspace = config();
    workspace.maximum_workspace_bytes = 1;
    assert!(
        LinearNetworkEventFlow::new(vec![1.0], vec![0.0], vec![], vec![1.0], workspace).is_err()
    );
}
