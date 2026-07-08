//! Sheaf-merge conformance (the lmp4.12 crown jewel; runs under the
//! `sheaf-merge` feature). Acceptance: coboundary-only mismatches
//! auto-reconcile to a RE-VERIFIED certificate-passing state; harmonic
//! mismatches report structural conflicts localized to the right cells
//! with both provenances; type-level collisions are caught before any
//! decomposition; degraded-gap merges are flagged low-confidence; the
//! Sev-0 adversarial case (reconciliation cannot reach a passing state)
//! ESCALATES rather than certifying falsely; trivial merges take the
//! fast paths; the kill-criterion harness measures the harmonic rate.
#![cfg(feature = "sheaf-merge")]

use std::collections::BTreeMap;

use fs_geom::sheaf_merge::{
    BranchState, Confidence, MergeOutcome, harmonic_conflict_rate, spectral_gap,
    three_way_merge,
};
use fs_geom::sheaf_repair::SheafSkeleton;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/sheaf-merge\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// 3-patch triangle (contractible: no harmonic space).
fn triangle() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2), (0, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

/// 4-patch ring (one cycle: harmonic space is 1-dimensional).
fn ring() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (1, 2), (2, 3), (0, 3)],
        triangles: vec![],
    }
}

fn branch(name: &str, mismatch: Vec<f64>) -> BranchState {
    BranchState {
        provenance: name.to_string(),
        mismatch,
        assignments: BTreeMap::new(),
    }
}

#[test]
fn sm_001_coboundary_auto_reconciles_with_reverified_certificate() {
    let sk = triangle();
    let base = vec![0.0; 3];
    // X re-gauges patch 1, Y re-gauges patch 2: pure coboundary edits.
    let x = branch("agent-x@c1", sk.d0(&[0.0, 0.02, 0.0]));
    let y = branch("agent-y@c2", sk.d0(&[0.0, 0.0, -0.015]));
    let out = three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6);
    match out {
        MergeOutcome::Resolved {
            merged,
            gauge,
            certificate,
            confidence,
        } => {
            assert!(certificate.post_norm <= certificate.tol, "re-verified");
            assert!(
                merged.iter().all(|v| v.abs() < 1e-9),
                "certificate-passing state: {merged:?}"
            );
            // The canonical gauge recovered both branches' offsets.
            assert!((gauge[1] - 0.02).abs() < 1e-9, "{gauge:?}");
            assert!((gauge[2] + 0.015).abs() < 1e-9, "{gauge:?}");
            assert!(matches!(confidence, Confidence::Normal { .. }));
        }
        other => panic!("coboundary edits must auto-resolve: {other:?}"),
    }
    verdict(
        "sm-001",
        "two gauge edits auto-reconcile; the certificate is re-verified (post-norm \
         under tol) and the recovered gauge matches both branches",
    );
}

#[test]
fn sm_002_harmonic_reports_structural_conflict_localized() {
    let sk = ring();
    let base = vec![0.0; 4];
    // X and Y push OPPOSING circulations around the cycle — their union
    // carries a net harmonic class no gauge can remove.
    let x = branch("agent-x@c7", vec![0.03, 0.03, 0.03, -0.03]);
    let y = branch("agent-y@c9", vec![0.01, 0.01, 0.01, -0.01]);
    let out = three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6);
    match out {
        MergeOutcome::Conflicted {
            structural,
            type_conflicts,
            ..
        } => {
            assert!(type_conflicts.is_empty());
            assert_eq!(structural.len(), 1);
            let c = &structural[0];
            assert_eq!(c.cells.len(), 4, "the whole cycle is the support");
            assert_eq!(
                c.parents,
                ("agent-x@c7".to_string(), "agent-y@c9".to_string()),
                "both provenances attached"
            );
        }
        other => panic!("a harmonic class must conflict: {other:?}"),
    }
    verdict(
        "sm-002",
        "opposing cycle edits surface ONE structural conflict supported on the full \
         cycle, carrying both parents",
    );
}

#[test]
fn sm_003_auto_resolution_iff_coboundary() {
    // Property sweep: seeded gauge-only edits ALWAYS resolve; any edit
    // with a cycle component NEVER auto-resolves on the ring.
    let sk = ring();
    let base = vec![0.0; 4];
    let mut state = 0xd00d_u64;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    for trial in 0..20 {
        let gx: Vec<f64> = (0..4).map(|_| 0.1 * lcg()).collect();
        let gy: Vec<f64> = (0..4).map(|_| 0.1 * lcg()).collect();
        let x = branch("x", sk.d0(&gx));
        let y = branch("y", sk.d0(&gy));
        let out = three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6);
        assert!(
            matches!(out, MergeOutcome::Resolved { .. }),
            "gauge-only trial {trial} must resolve"
        );
        // Now inject a cycle component into X.
        let mut mx = sk.d0(&gx);
        let eps = 0.02 + 0.05 * lcg().abs();
        for (k, v) in mx.iter_mut().enumerate() {
            *v += if k == 3 { -eps } else { eps };
        }
        let x2 = branch("x", mx);
        let y2 = branch("y", sk.d0(&gy));
        let out2 = three_way_merge(&sk, &base, &x2, &y2, None, 1e-9, 1e-6);
        assert!(
            matches!(out2, MergeOutcome::Conflicted { ref structural, .. } if !structural.is_empty()),
            "cycle-tainted trial {trial} must conflict"
        );
    }
    verdict(
        "sm-003",
        "20 seeded trials: auto-resolution applied IFF the union mismatch is a \
         coboundary — never on a harmonic class",
    );
}

#[test]
fn sm_004_sev0_escalates_instead_of_false_certificate() {
    // A COEXACT residue on the triangle: gauge reconciliation cannot
    // reach a passing state (circulation around the triple junction is
    // not a coboundary), and there is NO harmonic space here — the
    // naive path would "resolve" and attach a certificate over a state
    // that fails watertightness. The Sev-0 guard must escalate.
    let sk = triangle();
    let base = vec![0.0; 3];
    let circulation = sk.d1t(&[0.05]);
    let x = branch("x", circulation.clone());
    let y = branch("y", vec![0.0; 3]);
    // Wait: Y unchanged from base triggers the trivial path — perturb Y
    // slightly so the merge genuinely runs.
    let y = BranchState {
        mismatch: sk.d0(&[0.0, 1e-3, 0.0]),
        ..y
    };
    let out = three_way_merge(&sk, &base, &x, &y, None, 1e-6, 1e-6);
    match out {
        MergeOutcome::EscalatedUnresolved {
            post_norm,
            tol,
            fractions,
        } => {
            assert!(post_norm > tol, "the failure is real: {post_norm} > {tol}");
            assert!(
                fractions.1 > 0.5,
                "the residue is coexact (converter-side): {fractions:?}"
            );
        }
        other => panic!("Sev-0: must escalate, never certify falsely: {other:?}"),
    }
    // And the trivial fast paths themselves.
    let same = branch("s", sk.d0(&[0.0, 0.01, 0.0]));
    let t1 = three_way_merge(&sk, &base, &same, &same.clone(), None, 1e-9, 1e-6);
    assert!(
        matches!(t1, MergeOutcome::Trivial { reason, .. } if reason == "branches identical")
    );
    let unchanged = branch("u", base.clone());
    let t2 = three_way_merge(&sk, &base, &unchanged, &same, None, 1e-9, 1e-6);
    assert!(
        matches!(t2, MergeOutcome::Trivial { reason, .. } if reason == "X unchanged from base")
    );
    verdict(
        "sm-004",
        "a coexact residue escalates unresolved (never a false certificate); trivial \
         fast paths fire without decomposition",
    );
}

#[test]
fn sm_005_type_conflicts_and_degraded_gap() {
    let sk = triangle();
    let base = vec![0.0; 3];
    // Both branches edit the same load case differently: caught BEFORE
    // decomposition, even though the geometry would resolve.
    let mut x = branch("x", sk.d0(&[0.0, 0.01, 0.0]));
    x.assignments
        .insert("loadcase/cruise".to_string(), "2.5g".to_string());
    let mut y = branch("y", sk.d0(&[0.0, 0.0, 0.01]));
    y.assignments
        .insert("loadcase/cruise".to_string(), "3.0g".to_string());
    let out = three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6);
    match out {
        MergeOutcome::Conflicted {
            structural,
            type_conflicts,
            ..
        } => {
            assert!(structural.is_empty(), "no geometric conflict");
            assert_eq!(type_conflicts.len(), 1);
            assert_eq!(type_conflicts[0].key, "loadcase/cruise");
            assert_eq!(type_conflicts[0].x_value, "2.5g");
            assert_eq!(type_conflicts[0].y_value, "3.0g");
        }
        other => panic!("type collision must conflict: {other:?}"),
    }
    // Degraded gap: two clusters joined by ONE weak interface — the
    // weighted algebraic connectivity collapses and the merge is
    // flagged low-confidence (R5).
    let barbell = SheafSkeleton {
        n_patches: 6,
        edges: vec![(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)],
        triangles: vec![(0, 1, 2), (3, 4, 5)],
    };
    let weights = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1e-4];
    let gap = spectral_gap(&barbell, Some(&weights));
    assert!(gap < 1e-3, "weak-link gap is tiny: {gap}");
    let xb = branch("x", barbell.d0(&[0.0, 0.01, 0.0, 0.0, 0.0, 0.0]));
    let yb = branch("y", barbell.d0(&[0.0, 0.0, 0.0, 0.0, -0.01, 0.0]));
    let base_b = vec![0.0; 7];
    let out_b = three_way_merge(&barbell, &base_b, &xb, &yb, Some(&weights), 1e-9, 1e-3);
    match out_b {
        MergeOutcome::Resolved { confidence, .. } => {
            assert!(
                matches!(confidence, Confidence::LowGap { gap, threshold }
                    if gap < threshold),
                "degraded-gap merge must be flagged: {confidence:?}"
            );
        }
        other => panic!("the barbell gauge merge itself resolves: {other:?}"),
    }
    // A healthy complex at the same threshold is Normal.
    let healthy_gap = spectral_gap(&sk, None);
    assert!(healthy_gap > 1e-3, "triangle is well-coupled: {healthy_gap}");
    verdict(
        "sm-005",
        "same-key load-case collision caught before decomposition; the weak-link \
         barbell merge resolves but is FLAGGED LowGap (R5)",
    );
}

#[test]
fn sm_006_kill_criterion_harness() {
    // The measurement Proposal 10's kill criterion needs: harmonic
    // conflict rate over seeded realistic (gauge-dominated) edits.
    let ring4 = ring();
    let rate_ring = harmonic_conflict_rate(&ring4, 60, 0.1, 0xfeed);
    // Gauge-dominated edits on a cycle-bearing complex: the small
    // interface noise projects a little into the harmonic space, but
    // the rate must sit FAR below the 25% kill line.
    assert!(
        rate_ring < 0.25,
        "gauge-dominated merges must not structurally collide: {rate_ring}"
    );
    let tri = triangle();
    let rate_tri = harmonic_conflict_rate(&tri, 60, 0.1, 0xfeed);
    assert!(
        rate_tri.abs() < f64::EPSILON,
        "no harmonic space, no structural conflicts: {rate_tri}"
    );
    println!(
        "{{\"metric\":\"harmonic-conflict-rate\",\"ring\":{rate_ring:.3},\"triangle\":{rate_tri:.3},\
         \"kill_line\":0.25}}"
    );
    verdict(
        "sm-006",
        "the kill-criterion harness runs: triangle rate 0, ring rate below the 25% line, \
         both ledgered",
    );
}
