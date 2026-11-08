//! Physical modal viscosity, including rigid motion and coupled loss.
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{
    PhsError, modal_bank_ports, modal_bank_ports_with_damping_matrix,
    modal_bank_ports_with_damping_rates, modal_bank_with_damping_matrix,
    modal_bank_with_damping_rates,
};

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget { max_modes: n, max_omega_dt: 0.04 }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:.17e} versus {expected:.17e}, tolerance {tolerance:e}");
}

#[test]
fn physical_rates_preserve_a_dragged_rigid_coordinate_beside_an_oscillator() {
    let system = modal_bank_with_damping_rates(&[0.0, 80.0], &[3.0, 0.0], &[1.0, 0.0])
        .unwrap();
    let duration = 0.4;
    let flow = PreparedExactZoh::new(&system, duration, budget(4)).unwrap();
    let initial = [0.2, -0.4, 0.01, 0.3];
    let force = 0.9;
    let record = flow.step(&initial, &[force]).unwrap();
    let decay = (-3.0_f64 * duration).exp();
    let terminal = force / 3.0;
    close(record.x[0], initial[0] + terminal * duration
        + (initial[1] - terminal) * (1.0 - decay) / 3.0, 3e-13);
    close(record.x[1], terminal + (initial[1] - terminal) * decay, 3e-13);
    // More than five cycles in this exact interval: neither a midpoint phase
    // approximation nor a stiffness floor can reproduce these two limits.
    let phase = 80.0 * duration;
    close(record.x[2], initial[2] * phase.cos() + initial[3] * phase.sin() / 80.0, 3e-13);
    close(record.x[3], initial[3] * phase.cos() - 80.0 * initial[2] * phase.sin(), 3e-12);
    let offset = initial[1] - terminal;
    let heat = 3.0 * (terminal * terminal * duration
        + 2.0 * terminal * offset * (1.0 - decay) / 3.0
        + offset * offset * (1.0 - decay * decay) / 6.0);
    close(record.dissipated, heat, 3e-13);
    close(record.supplied, force * (record.x[0] - initial[0]), 3e-13);
    close(record.balance_residual(), 0.0, 3e-12);
}

fn rk4<const N: usize>(
    mut x: [f64; N], duration: f64, steps: usize,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let dt = duration / steps as f64;
    for _ in 0..steps {
        let a = rhs(&x);
        let b = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]));
        let c = rhs(&core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]));
        let d = rhs(&core::array::from_fn(|i| x[i] + dt * c[i]));
        for i in 0..N { x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0; }
    }
    x
}

#[test]
fn full_viscosity_matches_unequal_physical_masses_and_independent_heat_work() {
    // m=[4,9], K=diag(4*17^2,9*29^2), and three physical dashpots:
    // 40 and 60 to ground, 120 between the two moving coordinates.
    // Phi=diag(1/2,1/3), so Phi^T C_physical Phi=[[40,-20],[-20,20]].
    let damping = [40.0, -20.0, -20.0, 20.0];
    let drives: [&[f64]; 2] = [&[0.5, 0.0], &[0.0, 1.0 / 3.0]];
    let system = modal_bank_ports_with_damping_matrix(&[17.0, 29.0], &damping, &drives)
        .unwrap();
    let initial = [0.02, 0.06, -0.021, -0.03];
    let duration = 0.2;
    let forces = [0.7, -0.4];
    let reference = rk4([0.01, 0.03, -0.007, -0.01, 0.0, 0.0], duration, 4096, |x| {
        let [q0, v0, q1, v1, _, _] = *x;
        [v0, (forces[0] - 1156.0 * q0 - 160.0 * v0 + 120.0 * v1) / 4.0,
         v1, (forces[1] - 7569.0 * q1 + 120.0 * v0 - 180.0 * v1) / 9.0,
         40.0 * v0 * v0 + 60.0 * v1 * v1 + 120.0 * (v0 - v1).powi(2),
         forces[0] * v0 + forces[1] * v1]
    });
    let physical_energy = 0.5 * (4.0 * 0.03_f64.powi(2) + 9.0 * 0.01_f64.powi(2)
        + 1156.0 * 0.01_f64.powi(2) + 7569.0 * 0.007_f64.powi(2));
    close(system.hamiltonian(&initial), physical_energy, 2e-15);
    for (actual, expected) in system.output(&initial).iter().zip([0.03, -0.01]) {
        close(*actual, expected, 2e-16);
    }
    let flow = PreparedExactZoh::new(&system, duration, budget(4)).unwrap();
    let record = flow.step(&initial, &forces).unwrap();
    for i in 0..4 {
        close(record.x[i] / if i < 2 { 2.0 } else { 3.0 }, reference[i], 3e-11);
    }
    assert!(record.dissipated > 0.01);
    close(record.dissipated, reference[4], 3e-11);
    close(record.supplied, reference[5], 3e-12);
    close(record.balance_residual(), 0.0, 3e-12);
    let (_, r, _) = system.structure();
    assert_eq!([r[5], r[7], r[13], r[15]], damping);
}

#[test]
fn relative_velocity_drag_retains_its_lossless_rigid_nullspace() {
    let system = modal_bank_with_damping_matrix(
        &[0.0, 0.0], &[2.0, -2.0, -2.0, 2.0], &[1.0, 1.0],
    ).unwrap();
    let duration = 0.5;
    let flow = PreparedExactZoh::new(&system, duration, budget(4)).unwrap();
    let common = flow.step(&[0.3, 0.7, -0.2, 0.7], &[0.2]).unwrap();
    for (actual, expected) in common.x.iter().zip([0.675, 0.8, 0.175, 0.8]) {
        close(*actual, expected, 3e-14);
    }
    assert!(common.dissipated >= 0.0 && common.dissipated < 1e-28);
    close(common.supplied, 0.15, 3e-14);
    // Relative motion decays at 4/s, while both zero-stiffness positions
    // integrate that same velocity without a restoring force.
    let relative = flow.step(&[0.0, 0.5, 0.0, -0.5], &[0.0]).unwrap();
    let decay = (-4.0_f64 * duration).exp();
    let displacement = 0.125 * (1.0 - decay);
    for (actual, expected) in relative.x.iter().zip([
        displacement, 0.5 * decay, -displacement, -0.5 * decay,
    ]) { close(*actual, expected, 3e-14); }
    close(relative.dissipated, 0.25 * (1.0 - decay * decay), 3e-14);
    close(relative.balance_residual(), 0.0, 3e-14);
}

#[test]
fn explicit_rates_and_full_matrix_keep_the_legacy_positive_mode_physics() {
    let omega = [2.0, 10.0, 30.0, 80.0];
    let zeta = [0.0, 0.02, 1.0, 2.0];
    let rates: Vec<_> = omega.iter().zip(zeta).map(|(w, z)| 2.0 * z * w).collect();
    let drives: [&[f64]; 2] = [&[1.0, 0.3, -0.2, 0.1], &[0.0, -0.4, 0.8, 1.0]];
    let legacy = modal_bank_ports(&omega, &zeta, &drives).unwrap();
    let explicit = modal_bank_ports_with_damping_rates(&omega, &rates, &drives).unwrap();
    let mut matrix = [0.0; 16];
    for i in 0..4 { matrix[5 * i] = rates[i]; }
    let complete = modal_bank_ports_with_damping_matrix(&omega, &matrix, &drives).unwrap();
    for system in [&explicit, &complete] {
        assert_eq!(system.structure(), legacy.structure());
        assert_eq!(system.storage().exact_quadratic_matrix(),
            legacy.storage().exact_quadratic_matrix());
        assert!(system.resistance_witness().is_some());
    }
}

#[test]
fn malformed_or_nonpassive_physical_modal_banks_refuse() {
    for rates in [&[0.0][..], &[-1.0, 0.0], &[0.0, f64::NAN], &[0.0, f64::INFINITY]] {
        assert!(modal_bank_with_damping_rates(&[0.0, 2.0], rates, &[1.0, 1.0]).is_err());
    }
    for omega in [-1.0, f64::NAN, f64::INFINITY, f64::MAX, f64::MIN_POSITIVE] {
        assert!(modal_bank_with_damping_rates(&[omega], &[0.0], &[1.0]).is_err());
    }
    assert!(modal_bank_ports_with_damping_rates(&[0.0], &[0.0], &[]).is_err());
    assert!(modal_bank_with_damping_rates(&[], &[], &[]).is_err());
    assert!(modal_bank_with_damping_rates(&[0.0], &[0.0], &[f64::NAN]).is_err());
    for matrix in [
        &[1.0, 0.0][..], // incomplete matrix
        &[1.0, 2.0, 2.0, 1.0], // indefinite, despite nonnegative diagonal
        &[1.0, 0.0, 0.0, -1e-300], // a strong mode cannot hide active weak drag
        &[1.0, f64::NAN, f64::NAN, 1.0],
        &[1.0, 0.0, 0.0, f64::INFINITY],
    ] {
        assert!(modal_bank_with_damping_matrix(&[0.0, 2.0], matrix, &[1.0, 1.0]).is_err());
    }
    assert!(matches!(modal_bank_with_damping_matrix(
        &[0.0, 2.0], &[1.0, 0.0, 1e-15, 1.0], &[1.0, 1.0],
    ), Err(PhsError::NotSymmetric { .. })));
}
