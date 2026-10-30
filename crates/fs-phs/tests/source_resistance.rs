//! Physical source squares retain low-rank loss without changing the trajectory.
use fs_math::det;
use fs_phs::linear_exact::{
    LinearSplitBudget, PreparedExactZoh, integrated_resistance_loss_with_witness, split_quadratic,
    step_render,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, ResistanceWitness};

fn budget(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 0.04,
    }
}

fn identity(n: usize) -> Vec<f64> {
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        matrix[i * n + i] = 1.0;
    }
    matrix
}

fn source_matrix(n: usize, rows: &[f64]) -> Vec<f64> {
    let mut matrix = vec![0.0; n * n];
    for row in rows.chunks_exact(n) {
        for i in 0..n {
            for j in 0..n {
                matrix[i * n + j] += row[i] * row[j];
            }
        }
    }
    matrix
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * scale,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn rank_one_source_retains_original_matrix_nullspace_and_analytic_heat() {
    let (n, damping) = (4, 0.1);
    let column = [0.5, -2.0, 4.0, 0.125];
    // This is the original source assembly order c*b_i*b_j, not a repaired R.
    let matrix: Vec<_> = column
        .iter()
        .flat_map(|a| column.iter().map(move |b| damping * a * b))
        .collect();
    let rows: Vec<_> = column
        .iter()
        .map(|value| det::sqrt(damping) * value)
        .collect();
    let witness = ResistanceWitness::try_from_source_rows(n, matrix.clone(), rows).unwrap();
    let system = PortHamiltonian::new(
        n,
        0,
        vec![0.0; n * n],
        matrix.clone(),
        Vec::new(),
        Box::new(QuadraticStorage::new(identity(n), n).unwrap()),
    )
    .unwrap()
    .with_resistance_witness(witness.clone())
    .unwrap();
    assert_eq!(system.structure().1, matrix);
    let split = split_quadratic(&system, budget(n)).unwrap();
    let dt = 0.7;
    let initial = [0.7, -0.4, 0.12, 0.3];
    let norm: f64 = column.iter().map(|value| value * value).sum();
    let projection: f64 = initial.iter().zip(column).map(|(x, b)| x * b).sum();
    let decay = det::exp(-damping * norm * dt);
    let expected_heat = 0.5 * projection.powi(2) / norm * (1.0 - decay.powi(2));
    let step = step_render(&system, &split, &initial, &[], dt).unwrap();
    for i in 0..n {
        close(
            step.x[i],
            initial[i] + column[i] * projection / norm * (decay - 1.0),
            1.0,
        );
    }
    close(step.dissipated, expected_heat, expected_heat);
    assert!(step.dissipated > 0.0);
    close(step.balance_residual(), 0.0, system.hamiltonian(&initial));
    let channel =
        integrated_resistance_loss_with_witness(&system, &split, &initial, &[], dt, &witness)
            .unwrap();
    assert_eq!(channel.to_bits(), step.dissipated.to_bits());
    let prepared = PreparedExactZoh::new(&system, dt, budget(n)).unwrap();
    close(
        prepared.step(&initial, &[]).unwrap().dissipated,
        expected_heat,
        expected_heat,
    );

    // Powers-of-two column ratios make this physical source-null effort exact.
    let null = [4.0, 1.0, 0.0, 0.0];
    let silent = step_render(&system, &split, &null, &[], dt).unwrap();
    assert_eq!(silent.dissipated, 0.0);
    for i in 0..n {
        close(silent.x[i], null[i], 1.0);
    }
    assert_eq!(system.structure().1, matrix);
}

fn rk4<const N: usize>(
    mut state: [f64; N],
    dt: f64,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let h = dt / 2048.0;
    for _ in 0..2048 {
        let a = rhs(&state);
        let b = rhs(&core::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = rhs(&core::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = rhs(&core::array::from_fn(|i| state[i] + h * c[i]));
        for i in 0..N {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn rectangular_source_channels_share_full_forced_dynamics_and_material_composition() {
    let n = 3;
    let a = vec![0.4, -0.1, 0.2, -0.8, 0.2, -0.4];
    let b = vec![0.0, 0.3, -0.2, 0.0, -0.6, 0.4];
    let ra = source_matrix(n, &a);
    let rb = source_matrix(n, &b);
    let matrix: Vec<_> = ra.iter().zip(&rb).map(|(a, b)| a + b).collect();
    let source_a = ResistanceWitness::try_from_source_rows(n, ra, a.clone()).unwrap();
    let source_b = ResistanceWitness::try_from_source_rows(n, rb, b.clone()).unwrap();
    let rows: Vec<_> = a.iter().chain(&b).copied().collect();
    // Four physical squares in three states exercise the existing QR owner.
    let total = ResistanceWitness::try_from_source_rows(n, matrix.clone(), rows).unwrap();
    let j = vec![0.0, 1.0, -0.2, -1.0, 0.0, 0.5, 0.2, -0.5, 0.0];
    let g = vec![0.3, -0.2, 0.7];
    let system = PortHamiltonian::new(
        n,
        1,
        j.clone(),
        matrix.clone(),
        g.clone(),
        Box::new(QuadraticStorage::new(identity(n), n).unwrap()),
    )
    .unwrap()
    .with_resistance_witness(total)
    .unwrap();
    let split = split_quadratic(&system, budget(n)).unwrap();
    let initial = [0.4, -0.1, 0.2];
    let input = 0.25;
    let dt = 0.4;
    let reference = rk4(
        [initial[0], initial[1], initial[2], 0.0, 0.0, 0.0],
        dt,
        |x| {
            let mut out = [0.0; 6];
            for i in 0..n {
                out[i] = g[i] * input
                    + (0..n)
                        .map(|k| (j[i * n + k] - matrix[i * n + k]) * x[k])
                        .sum::<f64>();
                out[5] += input * g[i] * x[i];
            }
            for (heat, rows) in [(3, &a), (4, &b)] {
                for row in rows.chunks_exact(n) {
                    out[heat] += row.iter().zip(x).map(|(r, e)| r * e).sum::<f64>().powi(2);
                }
            }
            out
        },
    );
    let step = step_render(&system, &split, &initial, &[input], dt).unwrap();
    let heat_a =
        integrated_resistance_loss_with_witness(&system, &split, &initial, &[input], dt, &source_a)
            .unwrap();
    let heat_b =
        integrated_resistance_loss_with_witness(&system, &split, &initial, &[input], dt, &source_b)
            .unwrap();
    assert!(heat_a > 0.0 && heat_b > 0.0);
    close(heat_a, reference[3], reference[3]);
    close(heat_b, reference[4], reference[4]);
    close(heat_a + heat_b, step.dissipated, step.dissipated);
    close(step.supplied, reference[5], 1.0);
    for i in 0..n {
        close(step.x[i], reference[i], 1.0);
    }
    close(step.balance_residual(), 0.0, system.hamiltonian(&initial));
    assert_eq!(system.structure().1, matrix);

    // Genuine Maxwell history extends, rather than replaces, the same source R.
    let extended = system
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0, -0.25],
            stiffness: 2.0,
            relaxation_time_s: 0.05,
        }])
        .unwrap();
    for i in 0..n {
        for k in 0..n {
            assert_eq!(extended.structure().1[i * 4 + k], matrix[i * n + k]);
        }
    }
    let state = [initial[0], initial[1], initial[2], 0.1];
    let record = PreparedExactZoh::new(&extended, dt, budget(4))
        .unwrap()
        .step(&state, &[input])
        .unwrap();
    assert!(record.dissipated > 0.0);
    close(record.balance_residual(), 0.0, extended.hamiltonian(&state));
}

#[test]
fn source_witness_refuses_unsupported_entries_mismatched_laws_and_foreign_splits() {
    let witness = || {
        ResistanceWitness::try_from_source_rows(2, vec![1.0, -2.0, -2.0, 4.0], vec![1.0, -2.0])
            .unwrap()
    };
    let system = |damping| {
        PortHamiltonian::new(
            2,
            0,
            vec![0.0; 4],
            vec![damping, -2.0 * damping, -2.0 * damping, 4.0 * damping],
            Vec::new(),
            Box::new(QuadraticStorage::new(identity(2), 2).unwrap()),
        )
        .unwrap()
    };
    assert!(system(2.0).with_resistance_witness(witness()).is_err());
    assert!(ResistanceWitness::try_from_source_rows(2, vec![1.0; 4], vec![1.0]).is_err());
    assert!(ResistanceWitness::try_from_source_rows(2, identity(2), vec![f64::NAN, 0.0]).is_err());
    // A large unrelated entry cannot conceal unsupported negative tiny heat.
    assert!(
        ResistanceWitness::try_from_source_rows(2, vec![1.0, 0.0, 0.0, -1e-200], vec![1.0, 0.0])
            .is_err()
    );
    assert!(
        ResistanceWitness::try_from_source_rows(2, vec![1.0, 0.0, 1e-12, 1.0], identity(2))
            .is_err()
    );
    let first = system(1.0).with_resistance_witness(witness()).unwrap();
    let other = system(2.0);
    let split = split_quadratic(&first, budget(2)).unwrap();
    assert!(
        integrated_resistance_loss_with_witness(&other, &split, &[0.1, 0.2], &[], 0.1, &witness())
            .is_err()
    );
    let foreign = ResistanceWitness::try_from_source_rows(1, vec![1.0], vec![1.0]).unwrap();
    assert!(
        integrated_resistance_loss_with_witness(&first, &split, &[0.1, 0.2], &[], 0.1, &foreign)
            .is_err()
    );
}
