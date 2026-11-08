//! G1/G3/G4: exact endpoint-control tangents of original hereditary histories.
use fs_material::visco_tensor::stress_control::{
    PreparedTensorMixedControl as Mixed, PreparedTensorStressControl as Stress,
    TensorStressControlBudget as Budget,
};
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch, TensorPronyError};
use fs_material::{Tangent6, Voigt};

fn diagonal(values: Voigt) -> Tangent6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e}, expected={b:e}, tolerance={tolerance:e}"
    );
}
fn coupled() -> TensorProny {
    let mut eq = diagonal([1200.0, 1700.0, 1000.0, 900.0, 700.0, 1300.0]);
    eq[0][2] = -130.0;
    eq[2][0] = -130.0;
    eq[1][4] = 90.0;
    eq[4][1] = 90.0;
    let mut factors = [[[0.0; 6]; 6]; 2];
    factors[0][0] = [17.0, -3.0, 14.0, 7.0, 1.0, -2.0];
    factors[1][0] = [1.0, 12.0, -2.0, 6.0, 9.0, -4.0];
    factors[1][1] = [-3.0, 2.0, 8.0, -4.0, 0.0, 11.0];
    let branches = factors
        .into_iter()
        .zip([0.23, 1.1])
        .map(|(b, tau)| TensorPronyBranch::from_gram_rows(b, tau).unwrap())
        .collect();
    TensorProny::new_mandel(eq, branches, 0.1).unwrap()
}
fn diagonal_sls(eq: Voigt, branch: Voigt, tau: f64) -> TensorProny {
    TensorProny::new_mandel(
        diagonal(eq),
        vec![TensorPronyBranch::from_gram_rows(diagonal(branch.map(f64::sqrt)), tau).unwrap()],
        0.1,
    )
    .unwrap()
}

#[test]
fn g1_mixed_tangent_matches_actual_preloaded_history_differences_including_shears() {
    let mut initial = coupled().history();
    initial
        .advance(
            [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005],
            0.4,
            || false,
        )
        .unwrap();
    initial
        .advance(
            [0.0005, 0.0002, -0.0003, -0.0001, 0.0004, 0.0002],
            0.17,
            || false,
        )
        .unwrap();
    let before = initial.clone();
    let mask = [true, false, true, false, false, true];
    let control = Mixed::new(&initial, mask, Budget::default(), || false).unwrap();
    let duration = 0.47;
    let tangent = control.ramp_tangent(duration, || false).unwrap();
    let target = [-0.5, 0.0007, 0.4, -0.0003, 0.0002, -0.6];
    for column in 0..6 {
        let h = if mask[column] { 5e-5 } else { 5e-8 };
        let mut plus = target;
        let mut minus = target;
        plus[column] += h;
        minus[column] -= h;
        let a = control
            .advance(&mut initial.clone(), plus, duration, || false)
            .unwrap();
        let b = control
            .advance(&mut initial.clone(), minus, duration, || false)
            .unwrap();
        for row in 0..6 {
            let strain = (a.strain[row] - b.strain[row]) / (2.0 * h);
            let stress = (a.stress[row] - b.stress[row]) / (2.0 * h);
            close(
                tangent.strain[row][column],
                strain,
                3e-11 + 3e-9 * strain.abs(),
            );
            close(
                tangent.stress[row][column],
                stress,
                1e-8 + 3e-9 * stress.abs(),
            );
        }
    }
    assert_eq!(initial, before);
}

#[test]
fn g1_pure_control_limits_match_existing_strain_and_analytic_creep_tangents() {
    let law = coupled();
    let history = law.history();
    let all_strain = Mixed::new(&history, [false; 6], Budget::default(), || false).unwrap();
    for duration in [0.0, 0.37, 3.5] {
        let tangent = all_strain.ramp_tangent(duration, || false).unwrap();
        // tangent() is explicitly tensor-shear: normal/shear cross entries
        // differ from their Mandel values by the appropriate sqrt(2) ratios.
        let expected = law.tangent(duration).unwrap();
        for (row, expected_row) in expected.iter().enumerate() {
            for (column, &expected_value) in expected_row.iter().enumerate() {
                close(
                    tangent.strain[row][column],
                    if row == column { 1.0 } else { 0.0 },
                    2e-13,
                );
                close(
                    tangent.stress[row][column],
                    expected_value,
                    1e-10 + 2e-11 * expected_value.abs(),
                );
            }
        }
    }
    let eq = [400.0, 700.0, 1000.0, 1400.0, 900.0, 600.0];
    let branch = [1600.0, 900.0, 300.0, 500.0, 700.0, 2400.0];
    let tau = 0.6;
    let history = diagonal_sls(eq, branch, tau).history();
    let all_stress = Stress::new(&history, Budget::default(), || false).unwrap();
    for duration in [0.0, 0.37, 3.5] {
        let tangent = all_stress.ramp_tangent(duration, || false).unwrap();
        for row in 0..6 {
            let c0 = eq[row] + branch[row];
            let r = duration * eq[row] / (c0 * tau);
            let phi = if r == 0.0 { 1.0 } else { -(-r).exp_m1() / r };
            // Average of the SLS step-creep compliance over the prescribed
            // stress ramp. It is not the inverse of a linear-strain tangent.
            let compliance = 1.0 / c0 + branch[row] / (eq[row] * c0) * (1.0 - phi);
            for column in 0..6 {
                close(
                    tangent.strain[row][column],
                    if row == column { compliance } else { 0.0 },
                    2e-14,
                );
                close(
                    tangent.stress[row][column],
                    if row == column { 1.0 } else { 0.0 },
                    2e-12,
                );
            }
        }
    }
}

// Independent elimination of the six physical elastic constraints. The matrix
// has stress rows from C0 and identity rows for prescribed strain components.
fn solve(mut a: Tangent6, mut rhs: Voigt) -> Voigt {
    for column in 0..6 {
        let pivot = (column..6)
            .max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))
            .unwrap();
        a.swap(column, pivot);
        rhs.swap(column, pivot);
        let pivot_row = a[column];
        for row in column + 1..6 {
            let factor = a[row][column] / pivot_row[column];
            for k in column..6 {
                a[row][k] -= factor * pivot_row[k];
            }
            rhs[row] -= factor * rhs[column];
        }
    }
    let mut x = [0.0; 6];
    for row in (0..6).rev() {
        x[row] = (rhs[row] - (row + 1..6).map(|k| a[row][k] * x[k]).sum::<f64>()) / a[row][row];
    }
    x
}

#[test]
fn g1_zero_time_partition_and_stiff_weak_equilibrium_keep_physical_tangent_limits() {
    let law = coupled();
    let history = law.history();
    let mask = [true, false, true, false, true, false];
    let control = Mixed::new(&history, mask, Budget::default(), || false).unwrap();
    let tangent = control.ramp_tangent(0.0, || false).unwrap();
    let c0 = law.tangent(0.0).unwrap();
    let constraints = core::array::from_fn(|row| {
        if mask[row] {
            c0[row]
        } else {
            core::array::from_fn(|column| if row == column { 1.0 } else { 0.0 })
        }
    });
    for column in 0..6 {
        let strain = solve(
            constraints,
            core::array::from_fn(|row| if row == column { 1.0 } else { 0.0 }),
        );
        for row in 0..6 {
            let stress: f64 = (0..6).map(|k| c0[row][k] * strain[k]).sum();
            close(tangent.strain[row][column], strain[row], 2e-13);
            close(tangent.stress[row][column], stress, 3e-10);
        }
    }

    // The equilibrium stiffness is twenty orders below the branch stiffness.
    // In strain-controlled directions a subtraction C0 - relaxed_branch would
    // lose it completely. Stress-controlled directions have long-time creep.
    let (eq, branch, tau, duration) = (1e-20, 1.0, 1e-6, 1e20);
    let history = diagonal_sls([eq; 6], [branch; 6], tau).history();
    let control = Mixed::new(&history, mask, Budget::default(), || false).unwrap();
    let tangent = control.ramp_tangent(duration, || false).unwrap();
    let r = duration * (eq / (eq + branch)) / tau;
    let creep = 1.0 / (eq + branch) + branch / (eq * (eq + branch)) * (1.0 + (-r).exp_m1() / r);
    let strain_stiffness = eq + branch * (tau / duration) * -(-duration / tau).exp_m1();
    for (row, &stress_controlled) in mask.iter().enumerate() {
        let expected_strain = if stress_controlled { creep } else { 1.0 };
        let expected_stress = if stress_controlled { 1.0 } else { strain_stiffness };
        close(tangent.strain[row][row] / expected_strain, 1.0, 3e-9);
        close(tangent.stress[row][row] / expected_stress, 1.0, 3e-9);
        assert!(tangent.stress[row][row] > 0.0);
        for column in 0..6 {
            if row != column {
                assert_eq!(tangent.strain[row][column], 0.0);
                assert_eq!(tangent.stress[row][column], 0.0);
            }
        }
    }
}

#[test]
fn g4_tangent_invalid_duration_and_cancellation_leave_history_unchanged() {
    let mut history = coupled().history();
    history
        .advance(
            [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005],
            0.4,
            || false,
        )
        .unwrap();
    let before = history.clone();
    let mixed = Mixed::new(
        &history,
        [true, false, true, false, false, true],
        Budget::default(),
        || false,
    )
    .unwrap();
    let full = Stress::new(&history, Budget::default(), || false).unwrap();
    for duration in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(mixed.ramp_tangent(duration, || false).is_err());
        assert!(full.ramp_tangent(duration, || false).is_err());
    }
    let mut polls = 0_usize;
    let baseline = mixed
        .ramp_tangent(0.47, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls.div_ceil(2), polls] {
        let mut count = 0;
        assert!(matches!(
            mixed.ramp_tangent(0.47, || {
                count += 1;
                count >= stop
            }),
            Err(TensorPronyError::Cancelled)
        ));
    }
    assert!(matches!(
        full.ramp_tangent(0.47, || true),
        Err(TensorPronyError::Cancelled)
    ));
    let retry = mixed.ramp_tangent(0.47, || false).unwrap();
    assert_eq!(baseline.strain, retry.strain);
    assert_eq!(baseline.stress, retry.stress);
    assert_eq!(history, before);
}
