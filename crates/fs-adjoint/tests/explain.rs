//! Explanation-object conformance (the knh1.5 bead; runs under
//! `explanation-objects`). Acceptance: channels + residual = observed
//! ΔQoI within certified bounds — THE PERMANENT INVARIANT (the
//! Proposal-B kill criterion, measured over a case battery with 0%
//! failures allowed above the 10% line); the honesty gate refuses on
//! high-residual fixtures; every node re-derivable (stable
//! fingerprints, G5); the NL rendering is marked non-authoritative;
//! the flagship — far-field drag decomposition reconciling against the
//! analytic lifting-line envelope.
#![cfg(feature = "explanation-objects")]

use fs_adjoint::explain::{
    Elliptic1d, Explanation, ExplanationNode, LiftingLine, adjoint_attribution, drag_decomposition,
    finalize, provenance_attribution,
};
use fs_evidence::Color;
use std::panic::{AssertUnwindSafe, catch_unwind};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-adjoint/explain\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

#[test]
fn xp_001_permanent_invariant_kill_battery() {
    // 20 seeded conductivity edits; the adjoint engine must reconcile
    // on EVERY one (the kill line is 10% — we allow zero).
    let fixture = Elliptic1d { n: 120 };
    let mut lcg = 0xb00c_u64;
    let mut rnd = move || {
        lcg = lcg
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((lcg >> 11) as f64) / (1u64 << 53) as f64
    };
    let mut failures = 0usize;
    for case in 0..20 {
        let a0: Vec<f64> = (0..=120).map(|_| 1.0 + 0.5 * rnd()).collect();
        let a1: Vec<f64> = a0.iter().map(|a| a * (1.0 + 0.3 * (rnd() - 0.5))).collect();
        let observed =
            fixture.compliance(&fixture.solve(&a1)) - fixture.compliance(&fixture.solve(&a0));
        let channels = [
            ("left-half", (0..60).collect::<Vec<_>>()),
            ("right-half", (60..=120).collect::<Vec<_>>()),
        ];
        let nodes = adjoint_attribution(&fixture, &a0, &a1, &channels);
        // F2 (bead 9sf6): the Verified interval MUST equal the node's declared
        // (relative) bound. A fixed absolute ±1e-12 half-width is not a rigorous
        // enclosure once |value| is large (it would exclude the exact value).
        for node in &nodes {
            if let Color::Verified { lo, hi } = node.color {
                let half = 0.5 * (hi - lo);
                assert!(
                    (half - node.bound).abs() <= 1e-9 * node.bound.max(1e-30),
                    "Verified half-width {half:.3e} != declared bound {:.3e}",
                    node.bound
                );
            }
        }
        let explanation = finalize(nodes, observed, 1e-10);
        assert!(
            matches!(explanation, Explanation::Explained { .. }),
            "case {case}: the exact identity attributes fully"
        );
        if !explanation.reconciles() {
            failures += 1;
        }
    }
    println!(
        "{{\"metric\":\"kill-battery\",\"cases\":20,\"reconciliation_failures\":{failures},\
         \"kill_line\":\"10%\"}}"
    );
    assert_eq!(failures, 0, "the attribution engine never lies");
    verdict(
        "xp-001",
        "20-case battery: the exact bilinear adjoint identity attributes the full ΔJ to \
         channel masks with zero reconciliation failures (kill line 10%, measured 0%)",
    );
}

#[test]
fn xp_002_honesty_gate_refuses_hidden_channels() {
    // Declare only the left half, but edit BOTH halves: the residual
    // (the right half's true effect) exceeds the threshold — the gate
    // REFUSES rather than smearing it into the declared channel.
    let fixture = Elliptic1d { n: 120 };
    let a0 = vec![1.0f64; 121];
    let mut a1 = a0.clone();
    for (e, ae) in a1.iter_mut().enumerate() {
        *ae = if e < 60 { 1.2 } else { 0.7 }; // both halves edited
    }
    let observed =
        fixture.compliance(&fixture.solve(&a1)) - fixture.compliance(&fixture.solve(&a0));
    let declared_only = [("left-half", (0..60).collect::<Vec<_>>())];
    let nodes = adjoint_attribution(&fixture, &a0, &a1, &declared_only);
    let explanation = finalize(nodes, observed, 1e-6);
    assert!(
        matches!(explanation, Explanation::Refused { .. }),
        "the gate must refuse: {explanation:?}"
    );
    let Explanation::Refused {
        residual,
        threshold,
        ref partial,
    } = explanation
    else {
        return;
    };
    assert!(residual.abs() > threshold);
    assert_eq!(
        partial.len(),
        1,
        "the partial tree is forensics, not a claim"
    );
    // The narrative SAYS it refused.
    assert!(explanation.render_narrative().contains("REFUSED"));
    // Declaring the full mask set explains the same edit completely.
    let full = [
        ("left-half", (0..60).collect::<Vec<_>>()),
        ("right-half", (60..=120).collect::<Vec<_>>()),
    ];
    let ok = finalize(
        adjoint_attribution(&fixture, &a0, &a1, &full),
        observed,
        1e-6,
    );
    assert!(matches!(ok, Explanation::Explained { .. }) && ok.reconciles());
    verdict(
        "xp-002",
        "an undeclared channel's effect lands in the residual and the gate refuses; \
         declaring the full mask set explains the same edit completely",
    );
}

#[test]
fn xp_003_provenance_attribution_and_rederivability() {
    // Telescoping edit attribution is exact, and node fingerprints are
    // bit-stable under replay (G5 — the re-derivability witness).
    let edits = vec![
        ("thicken-spar".to_string(), 10.00, 10.40),
        ("trim-flange".to_string(), 10.40, 10.15),
        ("re-route-duct".to_string(), 10.15, 10.90),
    ];
    let nodes = provenance_attribution(&edits);
    let observed = 10.90 - 10.00;
    let explanation = finalize(nodes.clone(), observed, 1e-12);
    assert!(matches!(explanation, Explanation::Explained { .. }));
    assert!(explanation.reconciles(), "telescoping is exact");
    // Replay: identical fingerprints, node for node.
    let replay = provenance_attribution(&edits);
    for (a, b) in nodes.iter().zip(&replay) {
        assert_eq!(a.fingerprint, b.fingerprint, "re-derivable: {}", a.channel);
    }
    // Every node carries evidence links.
    assert!(nodes.iter().all(|n| !n.evidence.is_empty()));
    verdict(
        "xp-003",
        "edit attribution telescopes exactly to the observed change; fingerprints replay \
         bit-stable and every node carries its ledger evidence links",
    );
}

#[test]
fn xp_004_flagship_farfield_drag_decomposition() {
    // The lifting-line flagship: elliptic circulation at CL ~ 0.5,
    // AR = 8. The Trefftz wake integral must land on the analytic
    // envelope CDi = CL²/(π·AR) within its certified discretization
    // bound, and the three-channel decomposition must reconcile with
    // the near-field total.
    let (b, v, s_ref) = (8.0f64, 1.0f64, 8.0f64); // AR = 8
    // Γ0 chosen for CL ≈ 0.5: CL = π Γ0 b / (4 · ½ v S) ⇒ Γ0 = 2 CL v S/(π b).
    let cl_target = 0.5;
    let gamma0 = 2.0 * cl_target * v * s_ref / (std::f64::consts::PI * b);
    let wing = LiftingLine::elliptic(gamma0, b, v, s_ref, 400);
    let cl = wing.cl();
    assert!((cl - cl_target).abs() < 5e-3, "CL calibrated: {cl}");
    let cdi_analytic = cl * cl / (std::f64::consts::PI * wing.aspect_ratio());
    let cdi = wing.induced_drag_coefficient();
    let rel = (cdi - cdi_analytic).abs() / cdi_analytic;
    println!(
        "{{\"metric\":\"trefftz\",\"cdi\":{cdi:.6},\"analytic\":{cdi_analytic:.6},\
         \"rel\":{rel:.4}}}"
    );
    assert!(
        rel < 0.02,
        "the wake integral lands on CL^2/(pi AR): {cdi} vs {cdi_analytic}"
    );
    // Near-field 'observed' total = analytic induced + viscous strip.
    let (cf, wetted) = (0.006, 2.05);
    let cd_total = cdi_analytic + cf * wetted;
    let explanation = drag_decomposition(&wing, cf, wetted, cd_total, 2e-3);
    assert!(
        matches!(explanation, Explanation::Explained { .. }),
        "the flagship reconciles: {explanation:?}"
    );
    let Explanation::Explained {
        ref nodes,
        residual,
        ..
    } = explanation
    else {
        return;
    };
    assert!(explanation.reconciles(), "the permanent invariant holds");
    assert_eq!(nodes.len(), 3, "induced + viscous + declared-zero wave");
    assert!(
        nodes[2].channel.contains("declared zero"),
        "the wave channel is DECLARED, not omitted"
    );
    println!(
        "{{\"metric\":\"drag-decomposition\",\"induced\":{:.6},\"viscous\":{:.6},\
         \"wave\":0.0,\"residual\":{residual:.2e}}}",
        nodes[0].contribution, nodes[1].contribution
    );
    verdict(
        "xp-004",
        "flagship: the Trefftz wake integral matches CL^2/(pi AR) within 2%; the \
         induced/viscous/wave tree reconciles with the near-field total and the wave \
         channel is declared zero-subsonic rather than silently missing",
    );
}

#[test]
fn xp_005_narrative_is_non_authoritative() {
    let edits = vec![("polish".to_string(), 1.0, 1.1)];
    let explanation = finalize(provenance_attribution(&edits), 0.1, 1e-9);
    let text = explanation.render_narrative();
    assert!(
        text.starts_with("NON-AUTHORITATIVE RENDERING"),
        "the rendering leads with its own demotion"
    );
    assert!(text.contains("the explanation tree is the artifact"));
    verdict(
        "xp-005",
        "the natural-language rendering opens by declaring itself non-authoritative — \
         the tree is the artifact",
    );
}

#[test]
fn xp_006_edge_contracts_fail_fast_and_one_node_solves() {
    // The one-interior-node fixture is the smallest real SPD problem;
    // the tridiagonal solve must handle its empty off-diagonal.
    let fixture = Elliptic1d { n: 1 };
    let a0 = vec![1.0, 1.0];
    let a1 = vec![1.2, 0.9];
    let observed =
        fixture.compliance(&fixture.solve(&a1)) - fixture.compliance(&fixture.solve(&a0));
    let explanation = finalize(
        adjoint_attribution(&fixture, &a0, &a1, &[("all", vec![0, 1])]),
        observed,
        1e-12,
    );
    assert!(
        matches!(explanation, Explanation::Explained { .. }) && explanation.reconciles(),
        "one-node elliptic fixture must reconcile"
    );

    let short_conductivity = catch_unwind(AssertUnwindSafe(|| fixture.solve(&[1.0])));
    assert!(
        short_conductivity.is_err(),
        "conductivity length mismatches must fail fast"
    );
    let bad_channel = catch_unwind(AssertUnwindSafe(|| {
        adjoint_attribution(&fixture, &a0, &a1, &[("bad", vec![2])])
    }));
    assert!(
        bad_channel.is_err(),
        "out-of-range channel masks must fail fast"
    );
    let bad_node = catch_unwind(AssertUnwindSafe(|| {
        ExplanationNode::new(
            "bad",
            1.0,
            -1.0,
            Color::Estimated {
                estimator: "fixture".to_string(),
                dispersion: 0.1,
            },
            vec!["evidence".to_string()],
        )
    }));
    assert!(
        bad_node.is_err(),
        "negative contribution bounds must fail fast"
    );
    let bad_finalize = catch_unwind(AssertUnwindSafe(|| {
        finalize(provenance_attribution(&[]), f64::NAN, 1e-9)
    }));
    assert!(
        bad_finalize.is_err(),
        "non-finite observations must not produce explanation objects"
    );
    let bad_wing = catch_unwind(AssertUnwindSafe(|| {
        let _ = LiftingLine::elliptic(1.0, 1.0, 1.0, 1.0, 0);
    }));
    assert!(
        bad_wing.is_err(),
        "zero-station lifting-line fixtures must fail fast"
    );

    // Color is part of the re-derivation payload. Same numeric term
    // with different epistemic color must not collide.
    let verified = ExplanationNode::new(
        "same",
        1.0,
        0.1,
        Color::Verified { lo: 0.9, hi: 1.1 },
        vec!["evidence".to_string()],
    );
    let estimated = ExplanationNode::new(
        "same",
        1.0,
        0.1,
        Color::Estimated {
            estimator: "surrogate".to_string(),
            dispersion: 0.1,
        },
        vec!["evidence".to_string()],
    );
    assert_ne!(
        verified.fingerprint, estimated.fingerprint,
        "fingerprints must include evidence color"
    );
    verdict(
        "xp-006",
        "edge contracts fail fast; the one-node elliptic fixture reconciles; fingerprints include color",
    );
}
