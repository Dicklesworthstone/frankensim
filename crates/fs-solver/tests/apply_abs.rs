//! `LinearOp::apply_abs`: the absolute-term action behind a residual gate that knows the rounding floor of `fl(b - A x)`.

use fs_solver::{CsrOp, LinearOp};
use fs_sparse::Coo;

const U: f64 = f64::EPSILON / 2.0;

/// An operator that does not override `apply_abs`.
struct Unaware;

impl LinearOp for Unaware {
    fn n(&self) -> usize {
        3
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.copy_from_slice(x);
    }
}

#[test]
fn the_default_makes_no_claim_and_writes_nothing() {
    let mut y = vec![7.0; 3];
    assert_eq!(Unaware.apply_abs(&[1.0, -2.0, 3.0], &mut y), None);
    assert_eq!(y, vec![7.0; 3]);
}

fn mixed_sign_matrix() -> CsrOp {
    // Rows of 1, 2, 3, 4 and 5 stored entries with mixed signs; no entry is pushed twice, so no stored value is a sum.
    let n = 5;
    let mut coo = Coo::new(n, n);
    for row in 0..n {
        for column in 0..=row {
            let sign = if (row + column) % 2 == 0 { 1.0 } else { -1.0 };
            coo.push(row, column, sign * (1.0 + 0.37 * row as f64 + 0.11 * column as f64));
        }
    }
    CsrOp::general(coo.assemble())
}

fn probe(op: &impl LinearOp) -> Vec<Vec<f64>> {
    let n = op.n();
    (0..n)
        .map(|j| {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            let mut column = vec![0.0; n];
            op.apply(&e, &mut column);
            column
        })
        .collect()
}

#[test]
fn csr_apply_abs_equals_the_dense_absolute_action_and_counts_the_longest_row() {
    let op = mixed_sign_matrix();
    let n = op.n();
    let columns = probe(&op);
    let x = [1.5, -2.0, 0.25, -3.0, 4.0];
    let mut y = vec![0.0; n];
    let m = op.apply_abs(&x, &mut y).expect("CsrOp states its rounding chain");
    let mut probed_longest = 0;
    let mut ax = vec![0.0; n];
    op.apply(&x, &mut ax);
    for i in 0..n {
        let dense: f64 = (0..n).map(|j| columns[j][i].abs() * x[j].abs()).sum();
        probed_longest = probed_longest.max((0..n).filter(|&j| columns[j][i] != 0.0).count());
        assert!((y[i] - dense).abs() <= (m as f64 + 1.0) * U * dense, "row {i}: {} against {dense}", y[i]);
        // The triangle inequality, up to the rounding the count itself bounds.
        assert!(ax[i].abs() <= y[i] * (1.0 + (m as f64 + 1.0) * U) , "row {i}");
    }
    assert_eq!(m, probed_longest, "a CSR row is one fused chain over its stored entries");
    assert_eq!(m, 5);
}

/// The real-number value of `sum a_k b_k` for the given pairs, to about u^2: exact products by fused multiply-add, then Neumaier summation.
fn nearly_exact(pairs: &[(f64, f64)]) -> f64 {
    let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
    for &(a, b) in pairs {
        let product = a * b;
        let error = a.mul_add(b, -product);
        for term in [product, error] {
            let next = sum + term;
            compensation += if sum.abs() >= term.abs() { (sum - next) + term } else { (term - next) + sum };
            sum = next;
        }
    }
    sum + compensation
}

#[test]
fn csr_apply_abs_bounds_the_rounding_of_apply_by_gamma_m() {
    // Row 0 sums terms of size 3e7 that cancel to about 1e-1: the computed value is far below the absolute-term sum, and its error against the real-number value stays inside gamma_m * y.
    let n = 4;
    let mut coo = Coo::new(n, n);
    for column in 0..n {
        coo.push(0, column, if column % 2 == 0 { 1.0e8 / 3.0 } else { -1.0e8 / 3.0 });
        coo.push(column, column, 2.0);
    }
    let op = CsrOp::symmetric(coo.assemble());
    let columns = probe(&op);
    let x = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 7.0, 1.0 / 7.0];
    let mut y = vec![0.0; n];
    let m = op.apply_abs(&x, &mut y).unwrap();
    let mut ax = vec![0.0; n];
    op.apply(&x, &mut ax);
    let gamma = (m as f64) * U / (1.0 - (m as f64) * U);
    for i in 0..n {
        let pairs: Vec<(f64, f64)> = (0..n).filter(|&j| columns[j][i] != 0.0).map(|j| (columns[j][i], x[j])).collect();
        let exact = nearly_exact(&pairs);
        assert!((ax[i] - exact).abs() <= gamma * y[i], "row {i}: apply {} exact {exact} gamma*y {}", ax[i], gamma * y[i]);
    }
    assert!(y[0] > 1.0e7 && ax[0].abs() < 1.0e-3 * y[0], "row 0 cancels: |apply| {} against y {}", ax[0].abs(), y[0]);
}
