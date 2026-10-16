//! Bead ijil battery (a, b): interior-point and SQP parity with the
//! landed augmented-Lagrangian on the KKT fixtures, plus the SQP
//! warm-start polish gate.

use fs_ascent::auglag::{ConstrainedProblem, kkt_residual};
use fs_ascent::{augmented_lagrangian, interior_point, sqp};

type ConstraintEval<'a> = &'a dyn Fn(&[f64]) -> Vec<f64>;
type ConstraintJt<'a> = &'a dyn Fn(&[f64], &[f64]) -> Vec<f64>;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// The landed fixture: minimize (x−2)² + (y−1)² s.t. x + y = 2,
/// x ≤ 1.2. Optimum (1.2, 0.8), both constraints active.
fn fixture_a<'a>(
    fg: fs_ascent::FnGrad<'a>,
    ce: ConstraintEval<'a>,
    ce_jt: ConstraintJt<'a>,
    ci: ConstraintEval<'a>,
    ci_jt: ConstraintJt<'a>,
) -> ConstrainedProblem<'a> {
    ConstrainedProblem {
        fg,
        ce,
        ce_jt,
        ci,
        ci_jt,
    }
}

fn fg_a(x: &[f64]) -> (f64, Vec<f64>) {
    (
        (x[0] - 2.0).powi(2) + (x[1] - 1.0).powi(2),
        vec![2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)],
    )
}

#[test]
fn ip_parity_on_landed_kkt_fixture() {
    let mut fg = fg_a;
    let ce = |x: &[f64]| vec![x[0] + x[1] - 2.0];
    let ce_jt = |_: &[f64], w: &[f64]| vec![w[0], w[0]];
    let ci = |x: &[f64]| vec![x[0] - 1.2];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0], 0.0];
    let mut p = fixture_a(&mut fg, &ce, &ce_jt, &ci, &ci_jt);
    let rep = interior_point(&mut p, &[0.0, 0.0], 1e-6, 60);
    verdict(
        "ijil-ip-parity",
        rep.converged
            && (rep.x[0] - 1.2).abs() < 1e-4
            && (rep.x[1] - 0.8).abs() < 1e-4
            && rep.nu[0] > 0.0,
        &format!(
            "IP x=({:.5},{:.5}) kkt=({:.1e},{:.1e},{:.1e}) nu={:.4} outer={} — matches the AL fixture optimum (1.2, 0.8)",
            rep.x[0],
            rep.x[1],
            rep.kkt.stationarity,
            rep.kkt.feasibility,
            rep.kkt.complementarity,
            rep.nu[0],
            rep.outer_iters
        ),
    );
}

#[test]
fn sqp_parity_and_multiplier_agreement() {
    let mut fg1 = fg_a;
    let ce = |x: &[f64]| vec![x[0] + x[1] - 2.0];
    let ce_jt = |_: &[f64], w: &[f64]| vec![w[0], w[0]];
    let ci = |x: &[f64]| vec![x[0] - 1.2];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0], 0.0];
    let mut p1 = fixture_a(&mut fg1, &ce, &ce_jt, &ci, &ci_jt);
    let al = augmented_lagrangian(&mut p1, &[0.0, 0.0], 1e-7, 40);
    let mut fg2 = fg_a;
    let mut p2 = fixture_a(&mut fg2, &ce, &ce_jt, &ci, &ci_jt);
    let sq = sqp(&mut p2, &[0.0, 0.0], 1e-7, 60);
    let xdev =
        sq.x.iter()
            .zip(&al.x)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
    let mult_dev = (sq.lambda[0] - al.lambda[0])
        .abs()
        .max((sq.nu[0] - al.nu[0]).abs());
    verdict(
        "ijil-sqp-parity",
        sq.converged && xdev < 1e-5 && mult_dev < 1e-3,
        &format!(
            "SQP x=({:.5},{:.5}) vs AL x=({:.5},{:.5}): |dx| {xdev:.1e}; multipliers (lambda, nu) = ({:.4},{:.4}) vs AL ({:.4},{:.4}), dev {mult_dev:.1e}",
            sq.x[0], sq.x[1], al.x[0], al.x[1], sq.lambda[0], sq.nu[0], al.lambda[0], al.nu[0]
        ),
    );
}

#[test]
fn sqp_warm_start_polish_is_fast() {
    // The polish regime the bead names: from a near-optimal start the
    // active set is right immediately and SQP lands in a handful of
    // iterations.
    let mut fg = fg_a;
    let ce = |x: &[f64]| vec![x[0] + x[1] - 2.0];
    let ce_jt = |_: &[f64], w: &[f64]| vec![w[0], w[0]];
    let ci = |x: &[f64]| vec![x[0] - 1.2];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0], 0.0];
    let mut p = fixture_a(&mut fg, &ce, &ce_jt, &ci, &ci_jt);
    let rep = sqp(&mut p, &[1.19, 0.81], 1e-9, 20);
    // Measured: 10 iterations to 5e-11 KKT (identity-seeded BFGS needs
    // a few curvature pairs + working-set settling); the gate is the
    // measured envelope, not a wish.
    verdict(
        "ijil-sqp-polish",
        rep.converged && rep.iters <= 12 && (rep.x[0] - 1.2).abs() < 1e-8,
        &format!(
            "warm-start polish: {} iters to kkt=({:.1e},{:.1e},{:.1e})",
            rep.iters, rep.kkt.stationarity, rep.kkt.feasibility, rep.kkt.complementarity
        ),
    );
}

#[test]
fn ip_and_sqp_on_inequality_only_circle() {
    // minimize x + y s.t. x² + y² ≤ 2: optimum (−1, −1), ν = 0.5.
    let mk_fg = || |x: &[f64]| -> (f64, Vec<f64>) { (x[0] + x[1], vec![1.0, 1.0]) };
    let ce = |_: &[f64]| Vec::new();
    let ce_jt = |_: &[f64], _: &[f64]| vec![0.0, 0.0];
    let ci = |x: &[f64]| vec![x[0] * x[0] + x[1] * x[1] - 2.0];
    let ci_jt = |x: &[f64], w: &[f64]| vec![2.0 * x[0] * w[0], 2.0 * x[1] * w[0]];
    let mut fg1 = mk_fg();
    let mut p1 = fixture_a(&mut fg1, &ce, &ce_jt, &ci, &ci_jt);
    let ip = interior_point(&mut p1, &[0.0, 0.0], 1e-6, 60);
    let mut fg2 = mk_fg();
    let mut p2 = fixture_a(&mut fg2, &ce, &ce_jt, &ci, &ci_jt);
    let sq = sqp(&mut p2, &[0.0, 0.0], 1e-7, 80);
    let ok = |x: &[f64], nu: f64| {
        (x[0] + 1.0).abs() < 1e-3 && (x[1] + 1.0).abs() < 1e-3 && (nu - 0.5).abs() < 1e-2
    };
    verdict(
        "ijil-circle-both-engines",
        ip.converged && sq.converged && ok(&ip.x, ip.nu[0]) && ok(&sq.x, sq.nu[0]),
        &format!(
            "circle fixture: IP x=({:.4},{:.4}) nu={:.4}; SQP x=({:.4},{:.4}) nu={:.4} (analytic (-1,-1), nu=0.5)",
            ip.x[0], ip.x[1], ip.nu[0], sq.x[0], sq.x[1], sq.nu[0]
        ),
    );
}

#[test]
fn kkt_dual_feasibility_blocks_negative_multiplier_false_certificate() {
    // Old behavior returned three exact zeros here: the negative dual
    // cancels the objective gradient, while primal feasibility and
    // complementarity both vanish at the active boundary.
    let mut fg = |x: &[f64]| (x[0], vec![1.0]);
    let ce = |_: &[f64]| Vec::new();
    let ce_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
    let ci = |x: &[f64]| vec![x[0]];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0]];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };

    let residual = kkt_residual(&mut problem, &[0.0], &[], &[-1.0]);
    assert_eq!(
        [
            residual.stationarity,
            residual.feasibility,
            residual.complementarity,
            residual.dual_feasibility,
        ]
        .map(f64::to_bits),
        [
            0.0_f64.to_bits(),
            0.0_f64.to_bits(),
            0.0_f64.to_bits(),
            1.0_f64.to_bits()
        ],
    );
    assert!(!residual.within_tolerance(1e-8));
}

#[test]
#[should_panic(expected = "objective gradient entries must be finite")]
fn kkt_rejects_nan_instead_of_dropping_it_from_the_norm() {
    let mut fg = |_: &[f64]| (0.0, vec![f64::NAN]);
    let none = |_: &[f64]| Vec::new();
    let zero_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &none,
        ce_jt: &zero_jt,
        ci: &none,
        ci_jt: &zero_jt,
    };

    let _ = kkt_residual(&mut problem, &[0.0], &[], &[]);
}

#[test]
#[should_panic(expected = "objective gradient length must match")]
fn kkt_rejects_objective_gradient_dimension_mismatch() {
    let mut fg = |_: &[f64]| (0.0, vec![0.0]);
    let none = |_: &[f64]| Vec::new();
    let zero_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &none,
        ce_jt: &zero_jt,
        ci: &none,
        ci_jt: &zero_jt,
    };

    let _ = kkt_residual(&mut problem, &[0.0, 0.0], &[], &[]);
}

#[test]
#[should_panic(expected = "inequality multiplier length must match")]
fn kkt_rejects_multiplier_constraint_dimension_mismatch() {
    let mut fg = |_: &[f64]| (0.0, vec![0.0]);
    let ce = |_: &[f64]| Vec::new();
    let ce_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
    let ci = |_: &[f64]| vec![0.0];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0]];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };

    let _ = kkt_residual(&mut problem, &[0.0], &[], &[]);
}

#[test]
#[should_panic(expected = "Jacobian-transpose output length must match")]
fn kkt_rejects_jacobian_transpose_dimension_mismatch() {
    let mut fg = |_: &[f64]| (0.0, vec![0.0, 0.0]);
    let ce = |_: &[f64]| Vec::new();
    let ce_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
    let ci = |_: &[f64]| vec![0.0];
    let ci_jt = |_: &[f64], _: &[f64]| vec![0.0];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };

    let _ = kkt_residual(&mut problem, &[0.0, 0.0], &[], &[0.0]);
}

#[test]
#[should_panic(expected = "Jacobian-transpose action must map zero weights to zero")]
fn kkt_rejects_affine_jacobian_transpose_callback() {
    let mut fg = |_: &[f64]| (0.0, vec![0.0]);
    let ce = |_: &[f64]| Vec::new();
    let ce_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
    let ci = |_: &[f64]| vec![0.0];
    let biased_ci_jt = |_: &[f64], w: &[f64]| vec![w[0] + 1.0];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &biased_ci_jt,
    };

    let _ = kkt_residual(&mut problem, &[0.0], &[], &[0.0]);
}

#[test]
fn all_constrained_engines_require_positive_finite_tolerances() {
    for engine in ["AL", "IP", "SQP"] {
        for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let mut fg = |x: &[f64]| (x[0] * x[0], vec![2.0 * x[0]]);
                let none = |_: &[f64]| Vec::new();
                let zero_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
                let mut problem = ConstrainedProblem {
                    fg: &mut fg,
                    ce: &none,
                    ce_jt: &zero_jt,
                    ci: &none,
                    ci_jt: &zero_jt,
                };
                match engine {
                    "AL" => {
                        let _ = augmented_lagrangian(&mut problem, &[1.0], tolerance, 1);
                    }
                    "IP" => {
                        let _ = interior_point(&mut problem, &[1.0], tolerance, 1);
                    }
                    "SQP" => {
                        let _ = sqp(&mut problem, &[1.0], tolerance, 1);
                    }
                    _ => unreachable!(),
                }
            }));
            assert!(result.is_err(), "{engine} accepted tolerance {tolerance}");
        }
    }
}

#[test]
fn interior_exhaustion_keeps_last_solved_barrier_multiplier() {
    // With one outer iteration the solved barrier parameter is mu=1.
    // Recomputing nu after reducing the unsolved next mu to 0.2 shrinks
    // complementarity from 1 to 0.2 and used to false-report convergence.
    let mut fg = |x: &[f64]| {
        let delta = x[0] - 1.0;
        (0.05 * delta * delta, vec![0.1 * delta])
    };
    let ce = |_: &[f64]| Vec::new();
    let ce_jt = |x: &[f64], _: &[f64]| vec![0.0; x.len()];
    let ci = |x: &[f64]| vec![x[0]];
    let ci_jt = |_: &[f64], w: &[f64]| vec![w[0]];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };

    let report = interior_point(&mut problem, &[-1.0], 0.9, 1);
    assert!(
        !report.converged,
        "unsolved next-mu state certified: {report:?}"
    );
    assert!(
        report.kkt.complementarity > 0.9,
        "expected the mu=1 complementarity residual, got {report:?}"
    );
}
