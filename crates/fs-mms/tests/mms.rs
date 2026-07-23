//! fs-mms G1 harness conformance (bead 6nb.2): exact slope recovery, the
//! 0.2 gate in both directions, adjoint treatment, ladder refusals, a
//! REAL 1-D finite-difference Poisson MMS ladder end-to-end, and the
//! declared battery matrix with lintable gaps.

use fs_mms::{
    Coverage, LadderSide, MmsMatrix, MmsMatrixRow, ORDER_GATE_TOLERANCE, OrderFit, OrderGate,
    OrderVerdict, RefinementLadder, fit_order,
};

fn ladder(hs: &[f64], errors: &[f64]) -> RefinementLadder {
    RefinementLadder::new(hs.to_vec(), errors.to_vec()).expect("test ladder admits")
}

#[test]
fn synthetic_power_law_ladders_recover_their_exact_order() {
    for p in [1.0f64, 1.5, 2.0, 3.0, 4.0] {
        let hs: Vec<f64> = (0..5).map(|k| fs_math::det::powi(0.5, k)).collect();
        let errors: Vec<f64> = hs.iter().map(|h| 3.7 * fs_math::det::pow(*h, p)).collect();
        let fit = fit_order(&ladder(&hs, &errors));
        assert!(
            (fit.observed - p).abs() < 1e-10,
            "pure h^{p} ladder must fit exactly, got {}",
            fit.observed
        );
        assert!(
            fit.rms_residual < 1e-10,
            "pure power law is a straight line"
        );
    }
}

#[test]
fn the_gate_passes_within_0_2_and_build_fails_beyond_it() {
    let hs: Vec<f64> = (0..4).map(|k| fs_math::det::powi(0.5, k)).collect();
    let inside: Vec<f64> = hs.iter().map(|h| fs_math::det::pow(*h, 1.85)).collect();
    let outside: Vec<f64> = hs.iter().map(|h| fs_math::det::pow(*h, 1.75)).collect();
    let gate = OrderGate { theoretical: 2.0 };

    let verdict = gate
        .check("gate/inside", LadderSide::Primal, &ladder(&hs, &inside))
        .expect("1.85 vs 2.0 sits inside the 0.2 gate");
    assert!(verdict.deviation < ORDER_GATE_TOLERANCE);
    let line = verdict.json_line(true);
    assert!(line.contains("\"side\":\"primal\"") && line.contains("\"pass\":true"));

    let refused = gate
        .check("gate/outside", LadderSide::Primal, &ladder(&hs, &outside))
        .expect_err("1.75 vs 2.0 deviates by 0.25 > 0.2: the quiet death of accuracy");
    assert_eq!(refused.rule(), "mms-order-gate");
    assert!(refused.detail().contains("\"pass\":false"));

    let superconvergent: Vec<f64> = hs.iter().map(|h| fs_math::det::pow(*h, 2.31)).collect();
    assert!(
        gate.check(
            "gate/super",
            LadderSide::Primal,
            &ladder(&hs, &superconvergent)
        )
        .is_err(),
        "too-good orders fail too: superconvergence claims need their own theory"
    );
}

#[test]
fn json_lines_escape_every_caller_controlled_string() {
    let verdict = OrderGate { theoretical: 2.0 }
        .check(
            "quote\"\\slash\nsecond\u{0001}",
            LadderSide::Primal,
            &ladder(&[1.0, 0.5, 0.25], &[1.0, 0.25, 0.0625]),
        )
        .expect("exact second-order ladder");
    let order_line = verdict.json_line(true);
    assert_eq!(order_line.lines().count(), 1);
    assert!(!order_line.chars().any(|ch| ch < ' '));
    assert!(
        order_line.contains("\"case\":\"quote\\\"\\\\slash\\nsecond\\u0001\""),
        "{order_line}"
    );

    let matrix_line = MmsMatrix {
        rows: vec![MmsMatrixRow {
            frontend: "front\"end".into(),
            family: "fam\\ily".into(),
            bc: "bc\nline".into(),
            coverage: Coverage::Gap {
                reason: "gap\u{0008}\u{000c}\r\t\u{0002}".into(),
            },
        }],
    }
    .json_lines()
    .pop()
    .expect("one matrix row");
    assert_eq!(matrix_line.lines().count(), 1);
    assert!(!matrix_line.chars().any(|ch| ch < ' '));
    assert!(matrix_line.contains("\"frontend\":\"front\\\"end\""));
    assert!(matrix_line.contains("\"family\":\"fam\\\\ily\""));
    assert!(matrix_line.contains("\"bc\":\"bc\\nline\""));
    assert!(matrix_line.contains("\"detail\":\"gap\\b\\f\\r\\t\\u0002\""));

    let invalid_numeric_line = OrderVerdict {
        case: "invalid-numerics".into(),
        side: LadderSide::Adjoint,
        fit: OrderFit {
            observed: f64::NAN,
            intercept: f64::INFINITY,
            rms_residual: f64::NEG_INFINITY,
        },
        theoretical: f64::INFINITY,
        deviation: f64::NAN,
    }
    .json_line(true);
    assert!(invalid_numeric_line.contains("\"observed\":null"));
    assert!(invalid_numeric_line.contains("\"theoretical\":null"));
    assert!(invalid_numeric_line.contains("\"deviation\":null"));
    assert!(invalid_numeric_line.contains("\"rms_residual\":null"));
    assert!(invalid_numeric_line.contains("\"pass\":false"));
    assert!(!invalid_numeric_line.contains("NaN"));
    assert!(!invalid_numeric_line.contains("inf"));
}

#[test]
fn adjoint_ladders_get_the_identical_gate() {
    let hs: Vec<f64> = (0..4).map(|k| fs_math::det::powi(0.5, k)).collect();
    let dual_errors: Vec<f64> = hs
        .iter()
        .map(|h| 0.9 * fs_math::det::pow(*h, 3.92))
        .collect();
    let verdict = OrderGate { theoretical: 4.0 }
        .check(
            "adjoint/goal",
            LadderSide::Adjoint,
            &ladder(&hs, &dual_errors),
        )
        .expect("dual consistency verified, not assumed");
    assert!(verdict.json_line(true).contains("\"side\":\"adjoint\""));
}

#[test]
fn malformed_ladders_refuse_by_named_rule() {
    let refusal = |hs: &[f64], es: &[f64]| {
        RefinementLadder::new(hs.to_vec(), es.to_vec())
            .expect_err("must refuse")
            .rule()
            .to_owned()
    };
    assert_eq!(refusal(&[1.0, 0.5], &[1.0, 0.5]), "mms-ladder-shape");
    assert_eq!(refusal(&[1.0, 0.5, 0.25], &[1.0, 0.5]), "mms-ladder-shape");
    assert_eq!(
        refusal(&[1.0, 1.0, 0.5], &[1.0, 0.5, 0.2]),
        "mms-ladder-order"
    );
    assert_eq!(
        refusal(&[1.0, 0.5, 0.25], &[1.0, 0.0, 0.1]),
        "mms-ladder-domain"
    );
    assert_eq!(
        refusal(&[1.0, 0.5, 0.25], &[1.0, f64::NAN, 0.1]),
        "mms-ladder-domain"
    );
    assert_eq!(
        OrderGate { theoretical: 0.0 }
            .check(
                "gate/bad",
                LadderSide::Primal,
                &ladder(&[1.0, 0.5, 0.25], &[1.0, 0.25, 0.0625]),
            )
            .expect_err("zero theoretical order")
            .rule(),
        "mms-gate-domain"
    );
}

/// A REAL manufactured solution end-to-end: −u″ = f on (0,1) with
/// u* = sin(πx), f = π²·sin(πx), homogeneous Dirichlet, second-order
/// central differences on a refinement ladder — the observed L2 order
/// must gate green against the theoretical 2.0.
#[test]
fn real_fd_poisson_mms_ladder_gates_at_second_order() {
    let pi = std::f64::consts::PI;
    let mut hs = Vec::new();
    let mut errors = Vec::new();
    for n in [16usize, 32, 64, 128] {
        let h = 1.0 / n as f64;
        // Assemble and solve the tridiagonal system by the Thomas
        // algorithm (deterministic, direct).
        let m = n - 1;
        let mut sub = vec![-1.0f64; m];
        let mut diag = vec![2.0f64; m];
        let mut sup = vec![-1.0f64; m];
        let mut rhs: Vec<f64> = (1..n)
            .map(|i| {
                let x = i as f64 * h;
                h * h * pi * pi * fs_math::det::sin(pi * x)
            })
            .collect();
        for i in 1..m {
            let w = sub[i] / diag[i - 1];
            diag[i] -= w * sup[i - 1];
            rhs[i] -= w * rhs[i - 1];
        }
        let mut u = vec![0.0f64; m];
        u[m - 1] = rhs[m - 1] / diag[m - 1];
        for i in (0..m - 1).rev() {
            u[i] = (rhs[i] - sup[i] * u[i + 1]) / diag[i];
        }
        // Discrete L2 error against the manufactured solution.
        let mut ss = 0.0f64;
        for (i, ui) in u.iter().enumerate() {
            let x = (i + 1) as f64 * h;
            let d = ui - fs_math::det::sin(pi * x);
            ss += d * d;
        }
        hs.push(h);
        errors.push((h * ss).sqrt());
    }
    let verdict = OrderGate { theoretical: 2.0 }
        .check(
            "fd-poisson-1d/dirichlet",
            LadderSide::Primal,
            &RefinementLadder::new(hs, errors).expect("real ladder admits"),
        )
        .expect("central differences must gate at second order");
    println!("{}", verdict.json_line(true));
    assert!(
        verdict.deviation < 0.05,
        "FD Poisson sits well inside the gate"
    );
}

/// The declared battery matrix: current coverage in data, gaps explicit
/// and lintable — the coverage map the bead demands.
#[test]
fn battery_matrix_declares_coverage_and_lintable_gaps() {
    let matrix = MmsMatrix {
        rows: vec![
            MmsMatrixRow {
                frontend: "feec-body-fitted".into(),
                family: "p1-simplicial".into(),
                bc: "dirichlet".into(),
                coverage: Coverage::Covered {
                    test: "fs-feec/tests/feec_battery.rs::mms_poisson_primal_converges_at_second_order".into(),
                },
            },
            MmsMatrixRow {
                frontend: "cutfem-sdf".into(),
                family: "cut-p1".into(),
                bc: "sliver-cut".into(),
                coverage: Coverage::Covered {
                    test: "fs-cutfem/tests/elasticity.rs (log-log slope battery)".into(),
                },
            },
            MmsMatrixRow {
                frontend: "iga-patch".into(),
                family: "nurbs-p2".into(),
                bc: "mortar-seam".into(),
                coverage: Coverage::Gap {
                    reason: "multi-patch mortar MMS awaits the fs-iga mortar battery migration".into(),
                },
            },
            MmsMatrixRow {
                frontend: "opdsl-generated".into(),
                family: "any".into(),
                bc: "any".into(),
                coverage: Coverage::Gap {
                    reason: "forcing-term generation from operator definitions lands with fs-opdsl (tfz.4)".into(),
                },
            },
        ],
    };
    let gaps = matrix.gaps();
    assert_eq!(gaps.len(), 2, "every hole is explicit, none silent");
    let lines = matrix.json_lines();
    assert_eq!(lines.len(), 4);
    assert!(lines[0].contains("\"status\":\"covered\""));
    assert!(lines[3].contains("fs-opdsl"));
    for line in lines {
        println!("{line}");
    }
}
