//! Gradient-certificate conformance (the bk0o.3 bead; runs under the
//! `gradient-certs` feature). Acceptance: every emitted gradient
//! carries a color, verified ones carry an interval residual bound; FD
//! spot checks run via the falsifier pairing and agree within
//! conditioning-aware tolerance; a SEEDED TRANSPOSE BUG is caught and
//! BLOCKS merge; color assignment is correct across
//! differentiable/flagged/anchored paths; the no-falsifier-no-ship gate
//! covers the adjoint-gradient class.
#![cfg(feature = "gradient-certs")]

use std::collections::BTreeSet;

use fs_adjoint::certs::{
    Anchor, SparseLinear, adjoint_residual_bound, certify, fd_spot_checks, merge_gate,
};
use fs_adjoint::mitigate::grade_ops;
use fs_evidence::falsify::FalsifierRegistry;
use fs_evidence::{Color, ValidityDomain};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-adjoint/certs\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn smoothing_op(m: usize, n: usize) -> SparseLinear {
    let rows = (0..m)
        .map(|r| {
            (0..n)
                .map(|c| {
                    #[allow(clippy::cast_precision_loss)]
                    let w = 0.3 / (1.0 + (c as f64 - r as f64).abs());
                    (c, w)
                })
                .collect()
        })
        .collect();
    SparseLinear { rows, ncols: n }
}

#[test]
fn gc_001_interval_residual_bound_is_sound_and_tight() {
    let op = smoothing_op(10, 7);
    let bound = adjoint_residual_bound(&op, 24);
    // The registered pair IS a true transpose: the mathematical residual
    // is zero, so the verified enclosure only carries fp rounding —
    // tiny, and NEVER negative.
    assert!(bound >= 0.0);
    assert!(
        bound < 1e-12,
        "a true transpose pair has a rounding-level residual bound: {bound}"
    );
    verdict(
        "gc-001",
        "outward-rounded residual enclosure on a true transpose pair: bound in \
         [0, 1e-12] over 24 probes",
    );
}

#[test]
fn gc_002_seeded_transpose_bug_is_caught_and_blocks_merge() {
    // The objective J(x) = g·(Ax): its TRUE gradient is Aᵀg. The buggy
    // adjoint flips a sign mid-operator — exactly the class of bug the
    // FD falsifier exists to trip.
    let op = smoothing_op(10, 7);
    let g: Vec<f64> = (0..10).map(|i| 1.0 + 0.2 * f64::from(i as u8)).collect();
    let x = vec![0.4; 7];
    let f = |xs: &[f64]| -> f64 {
        op.rows
            .iter()
            .zip(&g)
            .map(|(row, gi)| gi * row.iter().map(|&(c, w)| w * xs[c]).sum::<f64>())
            .sum()
    };
    // Correct gradient: Aᵀ g.
    let mut grad_ok = vec![0.0f64; 7];
    for (row, gi) in op.rows.iter().zip(&g) {
        for &(c, w) in row {
            grad_ok[c] += w * gi;
        }
    }
    // Buggy gradient: a sign flip on column 3 (the seeded transpose bug).
    let mut grad_bad = grad_ok.clone();
    grad_bad[3] = -grad_bad[3];
    let checks_ok = fd_spot_checks(&f, &x, &grad_ok, 4, 0xabcd);
    let checks_bad = fd_spot_checks(&f, &x, &grad_bad, 4, 0xabcd);
    let cert_ok = certify(
        &grade_ops(&["convert/restrict"], &BTreeSet::new()),
        Some(adjoint_residual_bound(&op, 8)),
        checks_ok,
        None,
    );
    let cert_bad = certify(
        &grade_ops(&["convert/restrict"], &BTreeSet::new()),
        Some(adjoint_residual_bound(&op, 8)),
        checks_bad,
        None,
    );
    merge_gate(&cert_ok).expect("the correct gradient merges");
    let refusal = merge_gate(&cert_bad).expect_err("the sign bug must block merge");
    assert!(
        format!("{refusal}").contains("transpose or sign bug"),
        "teaches the cause: {refusal}"
    );
    // And a certificate with NO checks at all is refused too.
    let cert_none = certify(
        &grade_ops(&["convert/restrict"], &BTreeSet::new()),
        None,
        Vec::new(),
        None,
    );
    let refusal2 = merge_gate(&cert_none).expect_err("no checks, no merge");
    assert!(format!("{refusal2}").contains("mandatory"), "{refusal2}");
    // The COLOR must honor the falsifier, not just the merge gate (bead 9sf6):
    // the refuted gradient cannot be laundered into a Verified ledger color even
    // though its structural transpose-residual bound is tiny; the correct
    // gradient stays Verified.
    assert!(
        matches!(cert_ok.color, Color::Verified { .. }),
        "correct gradient stays Verified, got {:?}",
        cert_ok.color
    );
    assert!(
        matches!(cert_bad.color, Color::Estimated { .. }),
        "a gradient refuted by its own FD falsifier must NOT wear Verified, got {:?}",
        cert_bad.color
    );
    verdict(
        "gc-002",
        "seeded sign flip caught by 4 FD spot checks and blocked at the merge gate; \
         missing checks also refuse; the correct gradient merges",
    );
}

#[test]
fn gc_003_color_assignment_property() {
    let smooth = grade_ops(&["convert/restrict", "solver/spd"], &BTreeSet::new());
    let flagged = grade_ops(&["mesh/remesh"], &{
        let mut s = BTreeSet::new();
        s.insert("mesh/remesh".to_string());
        s
    });
    // Differentiable + interval-bounded → VERIFIED with the bound.
    let c1 = certify(&smooth, Some(3e-15), Vec::new(), None);
    match &c1.color {
        Color::Verified { lo, hi } => {
            assert!((*lo).abs() < f64::EPSILON && *hi > 0.0);
        }
        other => panic!("must be verified: {other:?}"),
    }
    assert!(c1.discontinuity.is_none());
    // Flagged remesh → ESTIMATED (inherited, never upgraded — that
    // would be laundering).
    let c2 = certify(&flagged, Some(1e-15), Vec::new(), None);
    assert!(
        matches!(c2.color, Color::Estimated { .. }),
        "flagged stays estimated even with a residual bound: {:?}",
        c2.color
    );
    assert!(c2.discontinuity.is_some());
    // Anchored to experimental data → VALIDATED in the regime.
    let c3 = certify(
        &smooth,
        Some(1e-15),
        Vec::new(),
        Some(&Anchor {
            dataset: "wind-tunnel-2026-03".to_string(),
            regime: ValidityDomain::unconstrained(),
        }),
    );
    match &c3.color {
        Color::Validated { dataset, .. } => assert_eq!(dataset, "wind-tunnel-2026-03"),
        other => panic!("anchored must be validated: {other:?}"),
    }
    // Smooth with NO evidence at all → ESTIMATED (folklore rule).
    let c4 = certify(&smooth, None, Vec::new(), None);
    assert!(
        matches!(&c4.color, Color::Estimated { estimator, .. } if estimator.contains("folklore") || estimator.contains("without")),
        "no certificate, no verified color: {:?}",
        c4.color
    );
    verdict(
        "gc-003",
        "colors: smooth+bounded=verified, flagged=estimated (no laundering), \
         anchored=validated, evidence-free=estimated",
    );
}

#[test]
fn gc_004_no_falsifier_no_ship_covers_gradients() {
    // The standard registry already pairs adjoint-gradient with the FD
    // spot check (qmao.4): gradients ship.
    let registry = FalsifierRegistry::standard();
    assert!(
        registry.ship_gate(&["adjoint-gradient"]).is_empty(),
        "the pairing exists in the standard registry"
    );
    // A hypothetical unpaired gradient class is blocked by name.
    let blocked = registry.ship_gate(&["adjoint-gradient", "hyper-gradient-v2"]);
    assert_eq!(blocked.len(), 1);
    assert!(blocked[0].contains("hyper-gradient-v2"));
    verdict(
        "gc-004",
        "no-falsifier-no-ship applies to gradients: adjoint-gradient ships via its FD \
         pairing, an unpaired class is blocked by name",
    );
}

#[test]
fn non_finite_residual_bound_cannot_mint_verified() {
    // Bead 9sf6 F5 regression: certify(Some(INFINITY)) used to wear
    // Verified{0, inf} — a vacuous certificate in the strongest color.
    use fs_adjoint::certs::certify;
    use fs_adjoint::mitigate::GradientGrade;
    let grade = GradientGrade::Smooth { route: vec![] };
    let cert = certify(&grade, Some(f64::INFINITY), vec![], None);
    assert!(
        !matches!(cert.color, fs_evidence::Color::Verified { .. }),
        "a non-finite bound must not certify: {:?}",
        cert.color
    );
    let cert_nan = certify(&grade, Some(f64::NAN), vec![], None);
    assert!(!matches!(
        cert_nan.color,
        fs_evidence::Color::Verified { .. }
    ));
    // A finite bound still verifies.
    let cert_ok = certify(&grade, Some(1e-12), vec![], None);
    assert!(matches!(
        cert_ok.color,
        fs_evidence::Color::Verified { .. }
    ));
}

#[test]
#[should_panic(expected = "ZERO probes")]
fn zero_probe_residual_bound_panics() {
    use fs_adjoint::certs::{SparseLinear, adjoint_residual_bound};
    let op = SparseLinear {
        rows: vec![vec![(0usize, 1.0f64)]],
        ncols: 1,
    };
    let _ = adjoint_residual_bound(&op, 0);
}
