//! The true-residual gate of `checked_solve_preconditioned`: a residual inside the rounding floor of `fl(b - A x)` is admitted, one above tol plus that floor is refused, and a solve that passes
//! the plain test, or whose operator states no floor, behaves exactly as before.
//!
//! The data are dyadic and the sums fused, so `apply` is exact and every residual below is the residual of the exact arithmetic: the zones do not depend on rounding noise.

use super::{admits_true_residual, checked_solve_preconditioned, evaluation_floor};
use crate::control::{SolveBudget, SolveControl};
use fs_solver::op::LinearOp;
use fs_sparse::precond::IdentityPrecond;
use std::ops::ControlFlow;

const TOL: f64 = 1e-11;
const U: f64 = f64::EPSILON / 2.0;
const K: f64 = 131072.0; // 2^17

/// A dense operator. `stated` decides whether it answers `apply_abs`; when it does, the sums run in the order of `apply`, so the rounding chain is the row length `n`.
struct Dense {
    n: usize,
    a: Vec<f64>,
    stated: bool,
}

impl LinearOp for Dense {
    fn n(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..self.n {
            let mut sum = 0.0;
            for j in 0..self.n {
                sum = self.a[i * self.n + j].mul_add(x[j], sum);
            }
            y[i] = sum;
        }
    }
    fn apply_abs(&self, x: &[f64], y: &mut [f64]) -> Option<usize> {
        if !self.stated {
            return None;
        }
        for i in 0..self.n {
            let mut sum = 0.0;
            for j in 0..self.n {
                sum = self.a[i * self.n + j].abs().mul_add(x[j].abs(), sum);
            }
            y[i] = sum;
        }
        Some(self.n)
    }
}

/// Springs of stiffness K = 2^17 in a chain, the last one grounded by 1: x = (1, 1, 1) solves it with b = (0, 1, 1), and the rows of |A||x| are of size K while b is of size 1, as in a stiff problem with a soft load.
fn stiff_chain(stated: bool) -> (Dense, Vec<f64>, Vec<f64>) {
    let a = vec![K, -K, 0.0, -K, 2.0 * K + 1.0, -K, 0.0, -K, K + 1.0];
    (Dense { n: 3, a, stated }, vec![0.0, 1.0, 1.0], vec![1.0, 1.0, 1.0])
}

/// The relative residual ||b - A x|| / ||b|| and the gate's decision for x = x* + 2^-p e_3 (p at most 52, below which 1 + 2^-p is not representable), or for x = x* itself. The residual is exactly
/// 2^-p (0, K, -(K+1)) up to the exact sums, so it is about 2^(17-p).
fn decide(op: &Dense, b: &[f64], x_star: &[f64], p: Option<i32>) -> (f64, bool) {
    let mut x = x_star.to_vec();
    if let Some(p) = p {
        assert!(p <= 52, "1 + 2^-{p} is not representable");
        x[2] += 2.0_f64.powi(-p);
    }
    let mut applied = vec![0.0; 3];
    op.apply(&x, &mut applied);
    let r: f64 = b.iter().zip(&applied).map(|(b, a)| (b - a) * (b - a)).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|b| b * b).sum::<f64>().sqrt();
    (r / nb, admits_true_residual(op, b, &x, &applied, 1.0, TOL))
}

/// The closed form of the floor in units of ||b||: gamma_{m+1} (1 + || |A||x| || / ||b||) with m = n = 3 and |A||x| recomputed densely here.
fn closed_form(op: &Dense, b: &[f64], x: &[f64]) -> f64 {
    let n = op.n;
    let abs: Vec<f64> = (0..n).map(|i| (0..n).map(|j| op.a[i * n + j].abs() * x[j].abs()).sum()).collect();
    let nb: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    let na: f64 = abs.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mu = (n as f64 + 1.0) * U;
    mu / (1.0 - mu) * (1.0 + na / nb)
}

#[test]
fn the_floor_is_the_closed_form_and_is_absent_when_the_operator_states_none() {
    let (op, b, x) = stiff_chain(true);
    let nb: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    let floor = evaluation_floor(&op, &b, &x, 1.0).expect("a stating operator has a floor");
    let expected = closed_form(&op, &b, &x) * nb;
    assert!((floor - expected).abs() <= 1.0e-12 * expected, "{floor:e} against {expected:e}");
    assert!(floor / nb > 10.0 * TOL, "the fixture has a floor well above tol: {:e}", floor / nb);
    assert_eq!(evaluation_floor(&stiff_chain(false).0, &b, &x, 1.0), None);
}

#[test]
fn a_residual_inside_the_rounding_floor_is_admitted_and_one_above_it_is_refused() {
    let (op, b, x) = stiff_chain(true);
    let allowance = TOL + closed_form(&op, &b, &x);
    // The exact solution: the plain test passes.
    let (rel, admitted) = decide(&op, &b, &x, None);
    assert!(rel == 0.0 && admitted);
    // Between tol and tol + floor (p = 51: 2^-34 = 5.8e-11): the plain test refuses, the floor admits.
    let (rel, admitted) = decide(&op, &b, &x, Some(51));
    assert!(rel > 5.0 * TOL && rel < allowance && admitted, "{rel:e} against {allowance:e}");
    // Above tol + floor (p = 48: 2^-31 = 4.7e-10): refused.
    let (rel, admitted) = decide(&op, &b, &x, Some(48));
    assert!(rel > allowance && !admitted, "{rel:e} against {allowance:e}");
}

#[test]
fn the_decision_flips_at_the_closed_form() {
    let (op, b, x) = stiff_chain(true);
    let allowance = TOL + closed_form(&op, &b, &x);
    // The residual halves with each p: the smallest admitted p and the one before it straddle the closed form within a factor two.
    let admitted: Vec<i32> = (40..=52).filter(|&p| decide(&op, &b, &x, Some(p)).1).collect();
    let first = *admitted.first().expect("some p is admitted");
    assert!(admitted.iter().all(|&p| p >= first), "admission is monotone in the residual");
    assert!(first > 40, "the largest residuals are refused");
    let (rel_in, _) = decide(&op, &b, &x, Some(first));
    let (rel_out, _) = decide(&op, &b, &x, Some(first - 1));
    assert!(rel_in <= allowance && allowance < rel_out, "{rel_in:e} <= {allowance:e} < {rel_out:e}");
}

#[test]
fn a_solve_stopped_at_a_thousand_times_tol_still_refuses_beside_its_derived_floor() {
    let (op, b, x) = stiff_chain(true);
    let floor = closed_form(&op, &b, &x);
    // p = 43: 2^-26 = 1.5e-8, about 1.5e3 tol.
    let (rel, admitted) = decide(&op, &b, &x, Some(43));
    assert!(rel > 1.0e3 * TOL && !admitted, "{rel:e}");
    assert!(floor < 0.05 * rel, "the floor {floor:e} sits far under {rel:e}");
}

#[test]
fn an_operator_that_states_no_floor_keeps_the_plain_gate() {
    let (op, b, x) = stiff_chain(false);
    // Between tol and tol + floor the stating operator admits (above); the unaware one refuses, exactly as before this change.
    let (rel, admitted) = decide(&op, &b, &x, Some(51));
    assert!(rel > 5.0 * TOL && !admitted, "{rel:e}");
    let (rel, admitted) = decide(&op, &b, &x, None);
    assert!(rel == 0.0 && admitted);
}

fn tridiagonal(n: usize, stated: bool) -> Dense {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = 4.0;
        if i + 1 < n {
            a[i * n + i + 1] = -1.0;
            a[(i + 1) * n + i] = -1.0;
        }
    }
    Dense { n, a, stated }
}

#[test]
fn a_solve_that_passes_the_plain_gate_returns_the_same_bits_whether_or_not_the_operator_states_a_floor() {
    let rhs: Vec<f64> = (0..8).map(|i| 1.0 + 0.25 * i as f64).collect();
    let mut solved = Vec::new();
    for stated in [false, true] {
        let op = tridiagonal(8, stated);
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget { total_iterations: 10_000, per_solve_iterations: 1_000 }, &mut poll);
        let x = checked_solve_preconditioned(&op, &IdentityPrecond, &rhs, TOL, "true-residual-test", &mut control).unwrap();
        solved.push(x.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
    }
    assert_eq!(solved[0], solved[1]);
}
