//! G1 disk mechanics and original ODE refinement; G3 frame/partition; G4 retry.
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::{
    LinearNetworkLink, vector::*,
};
use fs_exec::CancelGate;

fn config() -> VectorLinearNetworkConfig {
    VectorLinearNetworkConfig {
        maximum_modes: 32,
        maximum_connections: 16,
        maximum_contacts: 8,
        maximum_workspace_bytes: 1_000_000,
        maximum_cells: 100_000,
        maximum_nonlinear_iterations: 1_000_000,
        maximum_backtracks: 24,
        maximum_step_s: 10.0,
        maximum_cell_s: 0.002,
        minimum_column_pivot: 1e-10,
        force_absolute_tolerance_n: 2e-12,
        force_relative_tolerance: 2e-12,
        velocity_absolute_tolerance_m_s: 1e-12,
        maximum_abs_modal_coordinate: 100.0,
        maximum_abs_modal_velocity: 100.0,
        energy_absolute_tolerance_j: 1e-9,
        energy_relative_tolerance: 1e-8,
        maximum_numerical_work_uncertainty_j: 1e-7,
        maximum_numerical_dissipation_j: 100.0,
    }
}
fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "actual={a:.16e}, expected={b:.16e}, error={:.3e}",
        a - b
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn disk(c: VectorLinearNetworkConfig) -> VectorLinearNetworkFlow {
    VectorLinearNetworkFlow::new(
        vec![0.0; 2],
        vec![0.0; 2],
        vec![],
        vec![[vec![1.0, 0.0], vec![0.0, 1.0]]],
        c,
    )
    .unwrap()
}

#[test]
fn oblique_velocity_stops_on_one_disk_without_bounce_and_with_explicit_numerical_loss() {
    let mut c = config();
    c.maximum_cell_s = 0.125;
    let mut flow = disk(c);
    let s = flow
        .advance(&[0.0; 2], &[3.0, 4.0], &[0.0; 2], 2.0, &[2.5], &gate())
        .unwrap();
    close(flow.velocities()[0], 0.0, 1e-13);
    close(flow.velocities()[1], 0.0, 1e-13);
    // The exact discrete endpoint sum, distinct from the continuum distance 5.
    close(flow.positions()[0], 0.6 * 4.6875, 2e-13);
    close(flow.positions()[1], 0.8 * 4.6875, 2e-13);
    close(s.friction_dissipation_j, 11.71875, 1e-12);
    close(s.numerical_dissipation_j, 0.78125, 1e-12);
    close(s.contact_impulse_n_s[0][0], 3.0, 3e-13);
    close(s.contact_impulse_n_s[0][1], 4.0, 3e-13);
    close(s.contact_time_moment_n_s2[0][0], 3.0, 3e-13);
    close(s.contact_time_moment_n_s2[0][1], 4.0, 3e-13);
    close(
        s.initial_energy_j,
        s.friction_dissipation_j + s.numerical_dissipation_j,
        2e-12,
    );
    assert_eq!(s.regimes[0], VectorFrictionRegime::Sticking);
    let q = flow.positions().to_vec();
    let v = flow.velocities().to_vec();
    let hold = flow
        .advance(&q, &v, &[0.7, 0.2], 1.0, &[2.5], &gate())
        .unwrap();
    assert_eq!(hold.regimes[0], VectorFrictionRegime::Sticking);
    close(flow.positions()[0], q[0], 5e-13);
    close(flow.positions()[1], q[1], 5e-13);
    close(hold.contact_impulse_n_s[0][0], 0.7, 1e-12);
    close(hold.contact_impulse_n_s[0][1], 0.2, 1e-12);
    assert!(hold.friction_dissipation_j < 1e-12);
}

#[test]
fn force_diagonal_outside_disk_slides_even_when_both_scalar_components_fit() {
    let mut c = config();
    c.maximum_cell_s = 0.1;
    let mut flow = disk(c);
    let s = flow
        .advance(&[0.0; 2], &[0.0; 2], &[0.8, 0.8], 0.1, &[1.0], &gate())
        .unwrap();
    let expected = 0.1 * (0.8 - 1.0 / 2.0_f64.sqrt());
    close(flow.velocities()[0], expected, 2e-14);
    close(flow.velocities()[1], expected, 2e-14);
    close(s.maximum_contact_force_n[0], 1.0, 2e-13);
    assert_eq!(s.regimes[0], VectorFrictionRegime::Sliding);
    assert!(s.friction_dissipation_j > 0.0);
}

fn chain_columns() -> Vec<[Vec<f64>; 2]> {
    vec![
        [
            vec![1.0, 0.0, -1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0],
        ],
        [
            vec![0.0, 0.0, 1.0, 0.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 0.0, -1.0],
        ],
    ]
}
fn chain(c: VectorLinearNetworkConfig, cols: Vec<[Vec<f64>; 2]>) -> VectorLinearNetworkFlow {
    VectorLinearNetworkFlow::new(vec![0.0; 6], vec![0.0; 6], vec![], cols, c).unwrap()
}

#[test]
fn nonorthogonal_contacts_stop_jointly_preserving_both_momentum_components() {
    let mut c = config();
    c.maximum_cell_s = 0.005;
    let mut flow = chain(c, chain_columns());
    let v0 = [2.0, 1.0, -1.0, 2.0, -1.0, -3.0];
    let s = flow
        .advance(&[0.0; 6], &v0, &[0.0; 6], 4.0, &[2.0, 2.5], &gate())
        .unwrap();
    for v in flow.velocities() {
        close(*v, 0.0, 5e-12);
    }
    close(flow.positions().iter().step_by(2).sum(), 0.0, 2e-12);
    close(flow.positions().iter().skip(1).step_by(2).sum(), 0.0, 2e-12);
    for j in 0..2 {
        assert_eq!(s.regimes[j], VectorFrictionRegime::Sticking);
    }
    close(s.contact_impulse_n_s[0][0], 2.0, 4e-12);
    close(s.contact_impulse_n_s[0][1], 1.0, 4e-12);
    close(s.contact_impulse_n_s[1][0], 1.0, 4e-12);
    close(s.contact_impulse_n_s[1][1], 3.0, 4e-12);
    close(
        s.initial_energy_j,
        s.friction_dissipation_j + s.numerical_dissipation_j,
        5e-11,
    );
}

#[test]
fn tangent_basis_rotation_and_contact_order_preserve_the_coupled_trajectory() {
    let cols = chain_columns();
    let (c, s) = (0.8, 0.6);
    let rotated: Vec<_> = cols
        .iter()
        .rev()
        .map(|pair| {
            [
                pair[0]
                    .iter()
                    .zip(&pair[1])
                    .map(|(x, y)| c * x + s * y)
                    .collect(),
                pair[0]
                    .iter()
                    .zip(&pair[1])
                    .map(|(x, y)| -s * x + c * y)
                    .collect(),
            ]
        })
        .collect();
    let mut original = chain(config(), cols);
    let mut changed = chain(config(), rotated);
    let v = [2.0, 1.0, -1.0, 2.0, -1.0, -3.0];
    let a = original
        .advance(&[0.0; 6], &v, &[0.0; 6], 0.5, &[0.7, 1.3], &gate())
        .unwrap();
    let b = changed
        .advance(&[0.0; 6], &v, &[0.0; 6], 0.5, &[1.3, 0.7], &gate())
        .unwrap();
    for (x, y) in original.positions().iter().zip(changed.positions()) {
        close(*x, *y, 2e-12);
    }
    for (x, y) in original.velocities().iter().zip(changed.velocities()) {
        close(*x, *y, 2e-12);
    }
    close(a.friction_dissipation_j, b.friction_dissipation_j, 1e-11);
    close(a.numerical_dissipation_j, b.numerical_dissipation_j, 1e-11);
    for j in 0..2 {
        let t = a.contact_impulse_n_s[j];
        let r = b.contact_impulse_n_s[1 - j];
        close(r[0], c * t[0] + s * t[1], 3e-12);
        close(r[1], -s * t[0] + c * t[1], 3e-12);
    }
}

fn full_data() -> (
    Vec<f64>,
    Vec<f64>,
    Vec<LinearNetworkLink>,
    Vec<[Vec<f64>; 2]>,
) {
    (
        vec![4.0, 6.0, 2.0, 3.0, 1.0, 2.0],
        vec![0.1, 0.2, 0.15, 0.05, 0.12, 0.08],
        vec![
            LinearNetworkLink {
                column: vec![1.0, 0.3, -0.8, 0.2, 0.4, -0.1],
                stiffness_n_per_m: 3.0,
                damping_n_s_per_m: 0.4,
                rest_extension_m: 0.02,
            },
            LinearNetworkLink {
                column: vec![0.1, -0.4, 0.2, 0.8, -0.2, 0.3],
                stiffness_n_per_m: 1.5,
                damping_n_s_per_m: 0.15,
                rest_extension_m: -0.03,
            },
        ],
        chain_columns(),
    )
}
fn ode(
    y: &[f64; 17],
    k: &[f64],
    c: &[f64],
    links: &[LinearNetworkLink],
    cols: &[[Vec<f64>; 2]],
    f: &[f64],
) -> [f64; 17] {
    let mut d = [0.0; 17];
    for i in 0..6 {
        d[i] = y[6 + i];
        d[6 + i] = f[i] - k[i] * y[i] - c[i] * y[6 + i];
        d[12] += c[i] * y[6 + i] * y[6 + i];
    }
    for (j, l) in links.iter().enumerate() {
        let speed = dot(&l.column, &y[6..12]);
        let extension = dot(&l.column, &y[..6]) - l.rest_extension_m;
        let force = l.stiffness_n_per_m * extension + l.damping_n_s_per_m * speed;
        for i in 0..6 {
            d[6 + i] -= l.column[i] * force;
        }
        d[13 + j] = l.damping_n_s_per_m * speed * speed;
    }
    for (j, pair) in cols.iter().enumerate() {
        let u = [dot(&pair[0], &y[6..12]), dot(&pair[1], &y[6..12])];
        let speed = u[0].hypot(u[1]);
        assert!(speed > 0.1, "reference never reaches a nonsmooth stop");
        let cap = [0.2, 0.35][j];
        for i in 0..6 {
            d[6 + i] -= cap * (pair[0][i] * u[0] + pair[1][i] * u[1]) / speed;
        }
        d[15 + j] = cap * speed;
    }
    d
}
#[test]
fn full_original_elastic_pencil_and_each_physical_loss_converge_to_independent_rk4() {
    let (k, c, links, cols) = full_data();
    let q0 = [0.02, -0.01, 0.03, 0.01, -0.02, 0.04];
    let v0 = [2.0, 0.2, -0.3, 1.2, -0.4, -0.6];
    let f = [0.1, -0.2, 0.3, 0.05, -0.1, 0.15];
    let duration = 0.2;
    let mut reference = [0.0; 17];
    reference[..6].copy_from_slice(&q0);
    reference[6..12].copy_from_slice(&v0);
    let h = duration / 8192.0;
    for _ in 0..8192 {
        let a = ode(&reference, &k, &c, &links, &cols, &f);
        let b = ode(
            &core::array::from_fn(|i| reference[i] + 0.5 * h * a[i]),
            &k,
            &c,
            &links,
            &cols,
            &f,
        );
        let cc = ode(
            &core::array::from_fn(|i| reference[i] + 0.5 * h * b[i]),
            &k,
            &c,
            &links,
            &cols,
            &f,
        );
        let d = ode(
            &core::array::from_fn(|i| reference[i] + h * cc[i]),
            &k,
            &c,
            &links,
            &cols,
            &f,
        );
        for i in 0..17 {
            reference[i] += h * (a[i] + 2.0 * b[i] + 2.0 * cc[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    let mut numerical = Vec::new();
    for cell in [0.004, 0.002, 0.001] {
        let mut cfg = config();
        cfg.maximum_cell_s = cell;
        let mut flow =
            VectorLinearNetworkFlow::new(k.clone(), c.clone(), links.clone(), cols.clone(), cfg)
                .unwrap();
        let s = flow
            .advance(&q0, &v0, &f, duration, &[0.2, 0.35], &gate())
            .unwrap();
        let error = flow
            .positions()
            .iter()
            .chain(flow.velocities())
            .zip(&reference[..12])
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(error);
        numerical.push(s.numerical_dissipation_j);
        close(s.modal_damping_dissipation_j, reference[12], 0.005);
        for j in 0..2 {
            close(flow.connection_dissipation_j()[j], reference[13 + j], 0.005);
            close(
                s.contact_friction_dissipation_j[j],
                reference[15 + j],
                0.005,
            );
            assert!(flow.connection_numerical_dissipation_j()[j] > 0.0);
        }
        assert!(s.energy_residual_j.abs() < 1e-10);
    }
    assert!(
        errors[0] / errors[1] > 1.9 && errors[1] / errors[2] > 1.9,
        "{errors:?}"
    );
    assert!(
        numerical[0] / numerical[1] > 1.9 && numerical[1] / numerical[2] > 1.9,
        "{numerical:?}"
    );
}

#[test]
fn static_disk_preserves_internal_motion_and_original_connection_work() {
    // The contact holds both endpoint coordinates while an uncontacted original
    // oscillator moves. Its spring coupling changes the required disk reaction.
    let mut cfg = config();
    cfg.maximum_cell_s = 0.001;
    let mut flow = VectorLinearNetworkFlow::new(
        vec![0.0, 0.0, 4.0],
        vec![0.0, 0.0, 0.2],
        vec![LinearNetworkLink {
            column: vec![1.0, 0.0, -1.0],
            stiffness_n_per_m: 2.0,
            damping_n_s_per_m: 0.3,
            rest_extension_m: 0.0,
        }],
        vec![[vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]]],
        cfg,
    )
    .unwrap();
    let s = flow
        .advance(
            &[0.0; 3],
            &[0.0, 0.0, 1.0],
            &[0.0, 0.1, 0.0],
            0.3,
            &[3.0],
            &gate(),
        )
        .unwrap();
    assert_eq!(s.regimes[0], VectorFrictionRegime::Sticking);
    close(flow.positions()[0], 0.0, 2e-13);
    close(flow.positions()[1], 0.0, 2e-13);
    assert!(flow.positions()[2] > 0.2);
    assert!(flow.velocities()[2] > 0.5);
    assert!(s.modal_damping_dissipation_j > 0.0);
    assert!(flow.connection_dissipation_j()[0] > 0.0);
    assert!(flow.connection_numerical_dissipation_j()[0] > 0.0);
    assert!(s.friction_dissipation_j < 1e-12);
    // Backward Euler's endpoint spring work differs from stored-energy change
    // by exactly the separate connection numerical account.
    let connection_energy = flow.connection_forces_n()[0];
    assert!(connection_energy > 0.0);
    close(
        s.initial_energy_j + s.external_work_j,
        s.modal_energy_j
            + s.connection_energy_j
            + s.modal_damping_dissipation_j
            + s.connection_damping_dissipation_j
            + s.friction_dissipation_j
            + s.numerical_dissipation_j,
        2e-12,
    );
}

#[test]
fn equal_partition_retains_discrete_trajectory_and_vector_impulse_moment() {
    let mut whole = chain(config(), chain_columns());
    let mut split = chain(config(), chain_columns());
    let q0 = [0.0; 6];
    let v0 = [2.0, 1.0, -1.0, 2.0, -1.0, -3.0];
    let cap = [0.7, 1.3];
    let a = whole
        .advance(&q0, &v0, &[0.0; 6], 0.4, &cap, &gate())
        .unwrap();
    let b = split
        .advance(&q0, &v0, &[0.0; 6], 0.2, &cap, &gate())
        .unwrap();
    let q = split.positions().to_vec();
    let v = split.velocities().to_vec();
    let c = split
        .advance(&q, &v, &[0.0; 6], 0.2, &cap, &gate())
        .unwrap();
    for (x, y) in whole.positions().iter().zip(split.positions()) {
        close(*x, *y, 2e-12);
    }
    for (x, y) in whole.velocities().iter().zip(split.velocities()) {
        close(*x, *y, 2e-12);
    }
    close(
        a.numerical_dissipation_j,
        b.numerical_dissipation_j + c.numerical_dissipation_j,
        2e-12,
    );
    for j in 0..2 {
        for axis in 0..2 {
            close(
                a.contact_impulse_n_s[j][axis],
                b.contact_impulse_n_s[j][axis] + c.contact_impulse_n_s[j][axis],
                2e-12,
            );
            close(
                a.contact_time_moment_n_s2[j][axis],
                b.contact_time_moment_n_s2[j][axis]
                    + c.contact_time_moment_n_s2[j][axis]
                    + 0.2 * c.contact_impulse_n_s[j][axis],
                2e-12,
            );
        }
    }
}

#[test]
fn cancellation_late_loss_and_work_budget_refusals_leave_identical_retry() {
    let mut strict = config();
    strict.maximum_numerical_dissipation_j = 1e-8;
    let mut flow = disk(strict);
    assert!(matches!(
        flow.advance(&[0.0; 2], &[3.0, 4.0], &[0.0; 2], 0.1, &[2.5], &gate()),
        Err(VectorLinearNetworkError::Budget(
            "numerical time-discretization dissipation"
        ))
    ));
    let expected = flow
        .advance(&[0.0; 2], &[0.0; 2], &[0.3, 0.4], 0.1, &[2.5], &gate())
        .unwrap();
    let q = flow.positions().to_vec();
    let v = flow.velocities().to_vec();
    let cancelled = gate();
    cancelled.request();
    assert_eq!(
        flow.advance(&[0.0; 2], &[0.0; 2], &[0.3, 0.4], 0.1, &[2.5], &cancelled),
        Err(VectorLinearNetworkError::Cancelled)
    );
    let retry = flow
        .advance(&[0.0; 2], &[0.0; 2], &[0.3, 0.4], 0.1, &[2.5], &gate())
        .unwrap();
    assert_eq!(expected, retry);
    assert_eq!(q, flow.positions());
    assert_eq!(v, flow.velocities());
    let mut tiny = config();
    tiny.maximum_cells = 1;
    assert!(matches!(
        disk(tiny).advance(&[0.0; 2], &[1.0; 2], &[0.0; 2], 0.1, &[1.0], &gate()),
        Err(VectorLinearNetworkError::Budget(_))
    ));
    let mut finite = config();
    finite.maximum_nonlinear_iterations = 1;
    let mut coupled = chain(finite, chain_columns());
    assert!(matches!(
        coupled.advance(
            &[0.0; 6],
            &[2.0, 1.0, -1.0, 2.0, -1.0, -3.0],
            &[0.0; 6],
            0.1,
            &[0.7, 1.3],
            &gate()
        ),
        Err(VectorLinearNetworkError::Budget(_))
    ));
}

#[test]
fn dependent_columns_invalid_domains_and_strict_numerical_work_accuracy_refuse() {
    let bad = VectorLinearNetworkFlow::new(
        vec![0.0; 2],
        vec![0.0; 2],
        vec![],
        vec![[vec![1.0, 0.0], vec![2.0, 0.0]]],
        config(),
    );
    assert!(matches!(
        bad,
        Err(VectorLinearNetworkError::Invalid(
            "independent conditioned tangent columns"
        ))
    ));
    let mut flow = disk(config());
    assert!(matches!(
        flow.advance(&[0.0; 2], &[1.0; 2], &[0.0; 2], 0.1, &[-1.0], &gate()),
        Err(VectorLinearNetworkError::Invalid(_))
    ));
    let mut c = config();
    c.maximum_numerical_work_uncertainty_j = 1e-30;
    assert!(matches!(
        disk(c).advance(&[0.0; 2], &[3.0, 4.0], &[0.0; 2], 0.1, &[2.5], &gate()),
        Err(VectorLinearNetworkError::EnergyAccuracy { .. })
    ));
}
