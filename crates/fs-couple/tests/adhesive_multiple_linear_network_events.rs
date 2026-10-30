//! G1 analytic joint events and independent ODE work; G3 covariance and partition; G4 retry.
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::*;
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::multiple::*;
use fs_exec::CancelGate;

fn config() -> MultipleLinearNetworkEventConfig {
    MultipleLinearNetworkEventConfig {
        network: LinearNetworkEventConfig {
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
        },
        maximum_contacts: 8,
        minimum_column_pivot: 1e-12,
        maximum_phase_trials: 10000,
    }
}
fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "actual={a:.17e}, expected={b:.17e}, defect={:.3e}",
        a - b
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn chain(config: MultipleLinearNetworkEventConfig) -> MultipleLinearNetworkEventFlow {
    MultipleLinearNetworkEventFlow::new(
        vec![0.0; 3],
        vec![0.0; 3],
        vec![],
        vec![vec![1.0, -1.0, 0.0], vec![0.0, 1.0, -1.0]],
        config,
    )
    .unwrap()
}

#[test]
fn nonorthogonal_junctions_stop_jointly_with_the_analytic_reaction_impulse_and_first_moment() {
    let mut flow = chain(config());
    let s = flow
        .advance(
            &[0.0; 3],
            &[2.0, 1.0, -3.0],
            &[0.0; 3],
            &[LinearNetworkPhase::Fresh; 2],
            4.0,
            &[2.0; 2],
            &[1.0; 2],
            &gate(),
        )
        .unwrap();
    assert_eq!(s.contact_count, 2);
    assert_eq!(&s.phases[..2], &[LinearNetworkPhase::Sticking; 2]);
    assert_eq!(s.event_count, 2);
    for (i, t) in [1.0, 3.0].into_iter().enumerate() {
        let e = s.events[i].unwrap();
        assert_eq!(e.contact_index, i);
        assert_eq!(e.kind, LinearNetworkEventKind::Stopped);
        close(e.time_s, t, 3e-13);
    }
    for (a, b) in flow.positions().iter().zip([2.5, 2.0, -4.5]) {
        close(*a, b, 5e-13);
    }
    for v in flow.velocities() {
        close(*v, 0.0, 2e-14);
    }
    for i in 0..2 {
        close(s.contact_impulse_n_s[i], [2.0, 3.0][i], 2e-13);
        close(s.contact_time_moment_n_s2[i], [2.5, 4.5][i], 4e-13);
        close(s.contact_friction_dissipation_j[i], [0.5, 6.5][i], 4e-13);
        close(s.accumulated_slip_m[i], [0.5, 6.5][i], 4e-13);
    }
    close(s.friction_dissipation_j, 7.0, 5e-13);
    close(flow.positions().iter().sum(), 0.0, 5e-13);
    assert!(s.energy_residual_j.abs() < 3e-13);
}

#[test]
fn simultaneous_stops_are_one_physical_time_without_ordering_the_contacts() {
    let mut flow = MultipleLinearNetworkEventFlow::new(
        vec![0.0; 3],
        vec![0.0; 3],
        vec![],
        vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]],
        config(),
    )
    .unwrap();
    let s = flow
        .advance(
            &[0.0; 3],
            &[1.0, -2.0, 3.0],
            &[0.0; 3],
            &[LinearNetworkPhase::Fresh; 2],
            1.0,
            &[3.0, 6.0],
            &[2.0, 4.0],
            &gate(),
        )
        .unwrap();
    assert_eq!(s.event_count, 2);
    for e in s.events.iter().flatten() {
        close(e.time_s, 0.5, 2e-14);
    }
    for (a, b) in flow.positions().iter().zip([0.25, -0.5, 3.0]) {
        close(*a, b, 3e-14);
    }
    assert_eq!(flow.velocities(), &[0.0, 0.0, 3.0]);
    close(s.contact_impulse_n_s[0], 1.0, 2e-14);
    close(s.contact_impulse_n_s[1], -2.0, 2e-14);
    close(s.contact_time_moment_n_s2[0], 0.25, 2e-14);
    close(s.contact_time_moment_n_s2[1], -0.5, 2e-14);
    close(s.friction_dissipation_j, 2.5, 3e-14);
}

#[test]
fn elastic_internal_motion_changes_both_reactions_and_releases_only_the_overloaded_junction() {
    let mut flow = MultipleLinearNetworkEventFlow::new(
        vec![0.0, 0.0, 12.0],
        vec![0.0; 3],
        vec![],
        vec![vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 1.0]],
        config(),
    )
    .unwrap();
    let s = flow
        .advance(
            &[0.0; 3],
            &[-0.6, -0.6, 0.6],
            &[0.0; 3],
            &[LinearNetworkPhase::Fresh; 2],
            0.3,
            &[0.4, 10.0],
            &[0.2, 0.5],
            &gate(),
        )
        .unwrap();
    let release = (1.0_f64 / 3.0).asin() / 2.0;
    assert_eq!(s.event_count, 1);
    let e = s.events[0].unwrap();
    assert_eq!(e.contact_index, 0);
    assert_eq!(
        e.kind,
        LinearNetworkEventKind::SlipStarted { direction: -1 }
    );
    close(e.time_s, release, 3e-13);
    assert_eq!(s.phases[0], LinearNetworkPhase::Sliding { direction: -1 });
    assert_eq!(s.phases[1], LinearNetworkPhase::Sticking);
    close(flow.positions()[1] + flow.positions()[2], 0.0, 5e-15);
    close(flow.velocities()[1] + flow.velocities()[2], 0.0, 5e-15);
    assert!(s.contact_friction_dissipation_j[0] > 0.0);
    assert_eq!(s.contact_friction_dissipation_j[1], 0.0);
    // After release: q0''=0.2, q1=-q2, q2''=-6*q2+0.1.
    let q_release = 0.3 * (2.0 * release).sin();
    let v_release = 0.6 * (2.0 * release).cos();
    let t = 0.3 - release;
    let w = 6.0_f64.sqrt();
    let equilibrium = 0.1 / 6.0;
    let q2 =
        equilibrium + (q_release - equilibrium) * (w * t).cos() + v_release / w * (w * t).sin();
    let v2 = -(q_release - equilibrium) * w * (w * t).sin() + v_release * (w * t).cos();
    close(
        flow.positions()[0],
        -q_release - v_release * t + 0.1 * t * t,
        2e-13,
    );
    close(flow.velocities()[0], -v_release + 0.2 * t, 2e-13);
    close(flow.positions()[2], q2, 2e-13);
    close(flow.velocities()[2], v2, 2e-13);
}

#[test]
fn hidden_contact_stops_are_found_while_another_junction_remains_stuck() {
    let centre = 0.06_f64;
    let half = 0.01_f64;
    let mut flow = MultipleLinearNetworkEventFlow::new(
        vec![0.0, 1.0, 3.0],
        vec![0.0; 3],
        vec![],
        vec![vec![1.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]],
        config(),
    )
    .unwrap();
    let s = flow
        .advance(
            &[0.0, -centre.sin(), 0.1],
            &[-half.cos(), centre.cos(), 0.0],
            &[0.0; 3],
            &[LinearNetworkPhase::Fresh; 2],
            0.12,
            &[1e-4, 1.0],
            &[0.0, 0.5],
            &gate(),
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
    assert!(stops.iter().all(|e| e.contact_index == 0));
    close(flow.positions()[2], 0.1, 1e-15);
    close(flow.velocities()[2], 0.0, 1e-15);
    close(s.contact_impulse_n_s[1], -0.036, 1e-15);
    assert_eq!(s.friction_dissipation_j, 0.0);
}

#[test]
fn full_rank_stationary_constraints_retain_all_static_forces_in_one_high_frequency_cell() {
    let mut cfg = config();
    cfg.network.maximum_cells = 1;
    let mut flow = MultipleLinearNetworkEventFlow::new(
        vec![1e8; 2],
        vec![1.0; 2],
        vec![],
        vec![vec![1.0, 0.0], vec![0.5, 1.0]],
        cfg,
    )
    .unwrap();
    let q = [1e-8, -1e-8];
    let s = flow
        .advance(
            &q,
            &[0.0; 2],
            &[0.0; 2],
            &[LinearNetworkPhase::Fresh; 2],
            1.0,
            &[2.0; 2],
            &[0.5; 2],
            &gate(),
        )
        .unwrap();
    assert_eq!(flow.positions(), &q);
    assert_eq!(flow.velocities(), &[0.0; 2]);
    assert_eq!(s.cells, 1);
    assert_eq!(s.event_count, 0);
    close(s.contact_impulse_n_s[0], -1.5, 3e-15);
    close(s.contact_impulse_n_s[1], 1.0, 3e-15);
    close(s.contact_time_moment_n_s2[0], -0.75, 3e-15);
    close(s.contact_time_moment_n_s2[1], 0.5, 3e-15);
}

#[test]
fn partition_and_contact_unit_scaling_preserve_joint_motion_and_separate_work() {
    let mut whole = chain(config());
    let s = whole
        .advance(
            &[0.0; 3],
            &[2.0, 1.0, -3.0],
            &[0.0; 3],
            &[LinearNetworkPhase::Fresh; 2],
            4.0,
            &[2.0; 2],
            &[1.0; 2],
            &gate(),
        )
        .unwrap();
    for parts in [7, 23, 101] {
        let mut f = chain(config());
        let mut q = [0.0; 3];
        let mut v = [2.0, 1.0, -3.0];
        let mut phase = [LinearNetworkPhase::Fresh; 2];
        let mut loss = [0.0; 2];
        let mut impulse = [0.0; 2];
        let mut moment = [0.0; 2];
        for j in 0..parts {
            let x = f
                .advance(
                    &q,
                    &v,
                    &[0.0; 3],
                    &phase,
                    4.0 / parts as f64,
                    &[2.0; 2],
                    &[1.0; 2],
                    &gate(),
                )
                .unwrap();
            q.copy_from_slice(f.positions());
            v.copy_from_slice(f.velocities());
            phase.copy_from_slice(&x.phases[..2]);
            for i in 0..2 {
                loss[i] += x.contact_friction_dissipation_j[i];
                impulse[i] += x.contact_impulse_n_s[i];
                moment[i] += x.contact_time_moment_n_s2[i]
                    + j as f64 * 4.0 / parts as f64 * x.contact_impulse_n_s[i];
            }
        }
        for i in 0..3 {
            close(q[i], whole.positions()[i], 2e-12);
            close(v[i], whole.velocities()[i], 2e-12);
        }
        for i in 0..2 {
            close(loss[i], s.contact_friction_dissipation_j[i], 2e-12);
            close(impulse[i], s.contact_impulse_n_s[i], 2e-12);
            close(moment[i], s.contact_time_moment_n_s2[i], 3e-12);
        }
    }
    let scales = [1e4, 1e-3];
    let mut scaled = MultipleLinearNetworkEventFlow::new(
        vec![0.0; 3],
        vec![0.0; 3],
        vec![],
        vec![
            vec![scales[0], -scales[0], 0.0],
            vec![0.0, scales[1], -scales[1]],
        ],
        config(),
    )
    .unwrap();
    let t = scaled
        .advance(
            &[0.0; 3],
            &[2.0, 1.0, -3.0],
            &[0.0; 3],
            &[LinearNetworkPhase::Fresh; 2],
            4.0,
            &scales.map(|v| 2.0 / v),
            &scales.map(|v| 1.0 / v),
            &gate(),
        )
        .unwrap();
    for i in 0..3 {
        close(scaled.positions()[i], whole.positions()[i], 3e-12);
    }
    for i in 0..2 {
        close(
            t.contact_friction_dissipation_j[i],
            s.contact_friction_dissipation_j[i],
            2e-12,
        );
        close(
            t.accumulated_slip_m[i] / scales[i],
            s.accumulated_slip_m[i],
            2e-12,
        );
        close(
            t.contact_impulse_n_s[i] * scales[i],
            s.contact_impulse_n_s[i],
            2e-12,
        );
    }
}

#[test]
fn ambiguous_static_release_and_dependent_columns_refuse_without_corrupting_retry() {
    let ambiguous_columns = vec![vec![1.0, 0.0], vec![-0.5, 0.75_f64.sqrt()]];
    let make = || {
        MultipleLinearNetworkEventFlow::new(
            vec![0.0; 2],
            vec![0.0; 2],
            vec![],
            ambiguous_columns.clone(),
            config(),
        )
        .unwrap()
    };
    let mut flow = make();
    let phase = [LinearNetworkPhase::Fresh; 2];
    let force = [0.6, 0.9 / 0.75_f64.sqrt()];
    assert_eq!(
        flow.advance(
            &[0.0; 2],
            &[0.0; 2],
            &force,
            &phase,
            0.1,
            &[1.0; 2],
            &[0.2; 2],
            &gate()
        ),
        Err(LinearNetworkEventError::EventUnresolved)
    );
    let cancel = gate();
    cancel.request();
    assert_eq!(
        flow.advance(
            &[0.0; 2], &[0.0; 2], &[0.0; 2], &phase, 0.1, &[1.0; 2], &[0.2; 2], &cancel
        ),
        Err(LinearNetworkEventError::Cancelled)
    );
    let a = flow
        .advance(
            &[0.0; 2],
            &[0.0; 2],
            &[0.1; 2],
            &phase,
            0.1,
            &[1.0; 2],
            &[0.2; 2],
            &gate(),
        )
        .unwrap();
    let mut fresh = make();
    let b = fresh
        .advance(
            &[0.0; 2],
            &[0.0; 2],
            &[0.1; 2],
            &phase,
            0.1,
            &[1.0; 2],
            &[0.2; 2],
            &gate(),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(flow.positions(), fresh.positions());
    assert_eq!(flow.velocities(), fresh.velocities());
    assert!(
        MultipleLinearNetworkEventFlow::new(
            vec![0.0; 2],
            vec![0.0; 2],
            vec![],
            vec![vec![1.0, 1.0], vec![2.0, 2.0]],
            config()
        )
        .is_err()
    );
    assert!(
        MultipleLinearNetworkEventFlow::new(
            vec![0.0; 2],
            vec![0.0; 2],
            vec![],
            vec![vec![1.0, 0.0], vec![1.0, 1e-8]],
            config()
        )
        .is_err()
    );
}

#[test]
fn tiny_fresh_sliding_and_aggregate_phase_and_accuracy_budgets_are_preserved() {
    let columns = vec![vec![1.0, 0.0], vec![1.0, 1.0]];
    let make = |cfg| {
        MultipleLinearNetworkEventFlow::new(
            vec![0.0; 2],
            vec![0.0; 2],
            vec![],
            columns.clone(),
            cfg,
        )
        .unwrap()
    };
    let mut flow = make(config());
    let s = flow
        .advance(
            &[0.0; 2],
            &[1e-13, 0.0],
            &[0.0; 2],
            &[LinearNetworkPhase::Fresh; 2],
            1e-14,
            &[1.0; 2],
            &[0.5; 2],
            &gate(),
        )
        .unwrap();
    assert_eq!(s.event_count, 0);
    assert!(
        s.phases[..2]
            .iter()
            .all(|p| *p == LinearNetworkPhase::Sliding { direction: 1 })
    );
    assert!(s.accumulated_slip_m[..2].iter().all(|v| *v > 0.0));
    let mut limited = config();
    limited.maximum_phase_trials = 1;
    assert!(matches!(
        make(limited).advance(
            &[0.0; 2],
            &[0.0; 2],
            &[3.0, 4.0],
            &[LinearNetworkPhase::Fresh; 2],
            0.1,
            &[1.0; 2],
            &[0.5; 2],
            &gate()
        ),
        Err(LinearNetworkEventError::Budget(_))
    ));
    let mut strict = config();
    strict.network.maximum_numerical_work_uncertainty_j = s.numerical_work_uncertainty_j / 2.0;
    assert!(matches!(
        make(strict).advance(
            &[0.0; 2],
            &[1e-13, 0.0],
            &[0.0; 2],
            &[LinearNetworkPhase::Fresh; 2],
            1e-14,
            &[1.0; 2],
            &[0.5; 2],
            &gate()
        ),
        Err(LinearNetworkEventError::EnergyAccuracy { .. })
    ));
    let mut near_parallel = MultipleLinearNetworkEventFlow::new(
        vec![1.0, 2.0, 3.0],
        vec![0.0; 3],
        vec![],
        vec![vec![1.0, 0.0, 1.0], vec![1.0, 1e-4, 1.0]],
        config(),
    )
    .unwrap();
    assert!(matches!(
        near_parallel.advance(
            &[0.1, 0.0, -0.1],
            &[0.0; 3],
            &[0.0; 3],
            &[LinearNetworkPhase::Fresh; 2],
            0.1,
            &[1e6; 2],
            &[0.1; 2],
            &gate()
        ),
        Err(LinearNetworkEventError::EnergyAccuracy { .. })
    ));
}

#[test]
fn full_elastic_pencil_and_each_physical_loss_match_independent_runge_kutta() {
    let k = [2.0, 3.0, 5.0, 7.0];
    let c = [0.12, 0.13, 0.2, 0.23];
    let links = vec![
        LinearNetworkLink {
            column: vec![1.0, 0.4, -0.2, 0.1],
            stiffness_n_per_m: 4.0,
            damping_n_s_per_m: 0.17,
            rest_extension_m: 0.2,
        },
        LinearNetworkLink {
            column: vec![-0.1, 0.8, 0.3, -0.6],
            stiffness_n_per_m: 2.0,
            damping_n_s_per_m: 0.09,
            rest_extension_m: -0.1,
        },
    ];
    let columns = vec![vec![1.0, -1.0, 0.0, 0.0], vec![0.3, 0.2, -0.4, 1.0]];
    let f = [0.3, -0.1, 0.2, -0.1];
    let capacities = [0.2, 0.15];
    let mut flow = MultipleLinearNetworkEventFlow::new(
        k.to_vec(),
        c.to_vec(),
        links.clone(),
        columns.clone(),
        config(),
    )
    .unwrap();
    let q = [0.2, -0.1, 0.04, 0.3];
    let v = [1.0, -0.5, 0.2, 0.4];
    let dt = 0.08;
    let s = flow
        .advance(
            &q,
            &v,
            &f,
            &[LinearNetworkPhase::Fresh; 2],
            dt,
            &[0.3, 0.25],
            &capacities,
            &gate(),
        )
        .unwrap();
    assert_eq!(s.event_count, 0);
    // Assemble independent dense K/C once. The production owner applies each
    // physical rank-one connection directly and uses Taylor propagation.
    let mut kk = [[0.0; 4]; 4];
    let mut cc = [[0.0; 4]; 4];
    let mut ff = f;
    for i in 0..4 {
        kk[i][i] = k[i];
        cc[i][i] = c[i];
    }
    for link in &links {
        for i in 0..4 {
            ff[i] += link.stiffness_n_per_m * link.rest_extension_m * link.column[i];
            for j in 0..4 {
                kk[i][j] += link.stiffness_n_per_m * link.column[i] * link.column[j];
                cc[i][j] += link.damping_n_s_per_m * link.column[i] * link.column[j];
            }
        }
    }
    for p in 0..2 {
        for i in 0..4 {
            ff[i] -= columns[p][i] * capacities[p];
        }
    }
    let rhs = |x: [f64; 13]| {
        let mut r = [0.0; 13];
        for i in 0..4 {
            r[i] = x[4 + i];
            r[4 + i] = ff[i];
            for j in 0..4 {
                r[4 + i] -= kk[i][j] * x[j] + cc[i][j] * x[4 + j];
            }
            r[8] += c[i] * x[4 + i] * x[4 + i];
        }
        for p in 0..2 {
            let speed = dot(&links[p].column, &x[4..8]);
            r[9 + p] = links[p].damping_n_s_per_m * speed * speed;
            let slip = dot(&columns[p], &x[4..8]);
            assert!(slip > 0.1);
            r[11 + p] = capacities[p] * slip;
        }
        r
    };
    let add = |x: [f64; 13], r: [f64; 13], h: f64| core::array::from_fn(|i| x[i] + h * r[i]);
    let mut x = [0.0; 13];
    x[..4].copy_from_slice(&q);
    x[4..8].copy_from_slice(&v);
    let h = dt / 4000.0;
    for _ in 0..4000 {
        let a = rhs(x);
        let b = rhs(add(x, a, h / 2.0));
        let d = rhs(add(x, b, h / 2.0));
        let e = rhs(add(x, d, h));
        for i in 0..13 {
            x[i] += h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * d[i] + e[i]);
        }
    }
    for i in 0..4 {
        close(flow.positions()[i], x[i], 3e-13);
        close(flow.velocities()[i], x[4 + i], 3e-13);
    }
    close(s.modal_damping_dissipation_j, x[8], 2e-14);
    for p in 0..2 {
        close(flow.connection_dissipation_j()[p], x[9 + p], 2e-14);
        close(s.contact_friction_dissipation_j[p], x[11 + p], 2e-14);
    }
    close(
        s.external_work_j,
        dot(&f, &core::array::from_fn::<f64, 4, _>(|i| x[i] - q[i])),
        3e-14,
    );
    assert!(s.energy_residual_j.abs() < 3e-14);
}
