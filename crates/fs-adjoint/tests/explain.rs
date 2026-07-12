//! Explanation-object conformance (the knh1.5 bead; runs under
//! `explanation-objects`). Acceptance: channels + residual = observed
//! ΔQoI within honestly admitted bounds — THE PERMANENT INVARIANT (the
//! Proposal-B kill criterion, measured over a case battery with 0%
//! failures allowed above the 10% line); the honesty gate refuses on
//! high-residual fixtures; built-in nodes input-bound and replay-stable
//! (payload-integrity fingerprints, G5); the NL rendering is non-authoritative;
//! the flagship — far-field drag decomposition reconciling against the
//! analytic lifting-line envelope.
#![cfg(feature = "explanation-objects")]

use fs_adjoint::explain::{
    Elliptic1d, Explanation, ExplanationError, ExplanationNode, LiftingLine, adjoint_attribution,
    drag_decomposition, finalize, provenance_attribution,
};
use fs_evidence::{Color, ValidityDomain};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-adjoint/explain\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn elliptic_fixture(n: usize) -> Elliptic1d {
    Elliptic1d::new(n).expect("bounded test fixture")
}

#[test]
fn xp_001_unproved_rounding_allowance_never_certifies_the_kill_battery() {
    // Twenty seeded conductivity edits exercise the exact discrete identity,
    // but the floating-point solve/accumulation allowance is not proved. The
    // honesty gate must refuse rather than use that heuristic as coverage.
    let fixture = elliptic_fixture(120);
    let mut lcg = 0xb00c_u64;
    let mut rnd = move || {
        lcg = lcg
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((lcg >> 11) as f64) / (1u64 << 53) as f64
    };
    let mut honest_refusals = 0usize;
    for case in 0..20 {
        let a0: Vec<f64> = (0..=120).map(|_| 1.0 + 0.5 * rnd()).collect();
        let a1: Vec<f64> = a0.iter().map(|a| a * (1.0 + 0.3 * (rnd() - 0.5))).collect();
        let observed = fixture
            .compliance(&fixture.solve(&a1).expect("positive conductivity solves"))
            .expect("finite compliance")
            - fixture
                .compliance(&fixture.solve(&a0).expect("positive conductivity solves"))
                .expect("finite compliance");
        let channels = [
            ("left-half", (0..60).collect::<Vec<_>>()),
            ("right-half", (60..=120).collect::<Vec<_>>()),
        ];
        let nodes =
            adjoint_attribution(&fixture, &a0, &a1, &channels).expect("valid attribution inputs");
        assert!(
            nodes
                .iter()
                .all(|node| matches!(node.color(), Color::Estimated { .. })),
            "heuristic solve/accumulation allowances must not mint Verified"
        );
        let explanation = finalize(nodes, observed, 1e-10).expect("valid finalization inputs");
        assert!(explanation.is_structurally_valid(), "case {case}");
        assert!(matches!(
            explanation.receipt().aggregate_color,
            Color::Estimated { .. }
        ));
        if matches!(explanation, Explanation::Refused { .. }) {
            honest_refusals += 1;
        }
        assert!(!explanation.reconciles(), "case {case}: no proof, no claim");
    }
    println!(
        "{{\"metric\":\"kill-battery\",\"cases\":20,\"honest_refusals\":{honest_refusals},\
         \"certified_rounding_proof\":false}}"
    );
    assert_eq!(honest_refusals, 20, "the attribution engine must not lie");
    verdict(
        "xp-001",
        "20-case battery: exact discrete attribution remains inspectable, while every \
         unproved floating-point allowance refuses certification",
    );
}

#[test]
fn xp_002_honesty_gate_refuses_hidden_channels() {
    // Declare only the left half, but edit BOTH halves: the residual
    // (the right half's true effect) exceeds the threshold — the gate
    // REFUSES rather than smearing it into the declared channel.
    let fixture = elliptic_fixture(120);
    let a0 = vec![1.0f64; 121];
    let mut a1 = a0.clone();
    for (e, ae) in a1.iter_mut().enumerate() {
        *ae = if e < 60 { 1.2 } else { 0.7 }; // both halves edited
    }
    let observed = fixture
        .compliance(&fixture.solve(&a1).expect("positive conductivity solves"))
        .expect("finite compliance")
        - fixture
            .compliance(&fixture.solve(&a0).expect("positive conductivity solves"))
            .expect("finite compliance");
    let declared_only = [("left-half", (0..60).collect::<Vec<_>>())];
    let nodes =
        adjoint_attribution(&fixture, &a0, &a1, &declared_only).expect("valid partial attribution");
    let explanation = finalize(nodes, observed, 1e-6).expect("valid finalization inputs");
    assert!(
        matches!(explanation, Explanation::Refused { .. }),
        "the gate must refuse: {explanation:?}"
    );
    let Explanation::Refused {
        residual,
        ref receipt,
        ref partial,
        ..
    } = explanation
    else {
        return;
    };
    assert!(residual.abs() > receipt.effective_limit);
    assert!(explanation.is_structurally_valid());
    assert!(!explanation.reconciles(), "a refusal is never reconciled");
    assert_eq!(
        partial.len(),
        1,
        "the partial tree is forensics, not a claim"
    );
    // The narrative SAYS it refused.
    assert!(explanation.render_narrative().contains("REFUSED"));
    // Declaring the full mask set removes the hidden-channel residual, but the
    // unproved floating-point allowance still cannot certify an explanation.
    let full = [
        ("left-half", (0..60).collect::<Vec<_>>()),
        ("right-half", (60..=120).collect::<Vec<_>>()),
    ];
    let ok = finalize(
        adjoint_attribution(&fixture, &a0, &a1, &full).expect("valid full attribution"),
        observed,
        1e-6,
    )
    .expect("valid finalization inputs");
    assert!(matches!(ok, Explanation::Refused { .. }));
    assert!(ok.is_structurally_valid());
    assert!(!ok.reconciles());
    verdict(
        "xp-002",
        "an undeclared channel's effect lands in the residual and the gate refuses; \
         a full mask removes the modeling omission but remains an honest no-claim until \
         floating-point attribution is certified",
    );
}

#[test]
fn xp_003_provenance_attribution_and_rederivability() {
    // Caller-supplied edit attribution is deterministic and input-bound, but
    // remains Estimated until an authenticated ledger replay receipt exists.
    let edits = vec![
        ("thicken-spar".to_string(), 10.00, 10.40),
        ("trim-flange".to_string(), 10.40, 10.15),
        ("re-route-duct".to_string(), 10.15, 10.90),
    ];
    let nodes = provenance_attribution(&edits).expect("valid telescoping history");
    let observed = 10.90 - 10.00;
    let explanation = finalize(nodes.clone(), observed, 1e-12).expect("valid finalization inputs");
    assert!(matches!(explanation, Explanation::Explained { .. }));
    assert!(explanation.reconciles(), "telescoping is exact");
    assert!(
        nodes
            .iter()
            .all(|node| { node.bound() > 0.0 && matches!(node.color(), Color::Estimated { .. }) })
    );
    assert!(matches!(
        explanation.receipt().aggregate_color,
        Color::Estimated { .. }
    ));
    // Replay: identical fingerprints, node for node.
    let replay = provenance_attribution(&edits).expect("valid replay history");
    for (a, b) in nodes.iter().zip(&replay) {
        assert_eq!(
            a.fingerprint(),
            b.fingerprint(),
            "re-derivable: {}",
            a.channel()
        );
    }
    // Every node carries evidence links.
    assert!(nodes.iter().all(|n| !n.evidence().is_empty()));
    verdict(
        "xp-003",
        "caller edit attribution is input-bound and ULP-enclosed but Estimated; fingerprints \
         replay bit-stable and every node carries its declared evidence links",
    );
}

#[test]
fn xp_004_flagship_farfield_drag_decomposition() {
    // The lifting-line flagship: elliptic circulation at CL ~ 0.5,
    // AR = 8. The Trefftz wake integral must land near the analytic
    // envelope CDi = CL²/(π·AR), but that measured conformance result is
    // not a certified discretization enclosure.
    let (b, v, s_ref) = (8.0f64, 1.0f64, 8.0f64); // AR = 8
    // Γ0 chosen for CL ≈ 0.5: CL = π Γ0 b / (4 · ½ v S) ⇒ Γ0 = 2 CL v S/(π b).
    let cl_target = 0.5;
    let gamma0 = 2.0 * cl_target * v * s_ref / (std::f64::consts::PI * b);
    let wing = LiftingLine::elliptic(gamma0, b, v, s_ref, 400).expect("valid wing fixture");
    let cl = wing.cl().expect("finite lift coefficient");
    assert!((cl - cl_target).abs() < 5e-3, "CL calibrated: {cl}");
    let cdi_analytic =
        cl * cl / (std::f64::consts::PI * wing.aspect_ratio().expect("finite aspect ratio"));
    let cdi = wing
        .induced_drag_coefficient()
        .expect("finite induced drag");
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
    let explanation = drag_decomposition(
        &wing,
        cf,
        wetted,
        0.3,
        "mach-probe:flagship",
        cd_total,
        2e-3,
    )
    .expect("valid drag decomposition inputs");
    let Explanation::Refused {
        ref partial,
        residual,
        ref receipt,
        ..
    } = explanation
    else {
        panic!("heuristic channel bounds must not certify the analytic/model discrepancy");
    };
    assert!(explanation.is_structurally_valid());
    assert!(!explanation.reconciles());
    let nodes = partial;
    assert_eq!(nodes.len(), 3, "induced + viscous + declared-zero wave");
    assert!(
        nodes
            .iter()
            .all(|node| matches!(node.color(), Color::Estimated { .. })),
        "no heuristic drag channel may carry built-in Verified authority"
    );
    assert_eq!(
        receipt.certified_coverage.to_bits(),
        receipt.aggregation_roundoff.to_bits(),
        "the Trefftz O(1/n) heuristic contributes no certified coverage"
    );
    assert!(residual.abs() > receipt.effective_limit);
    assert!(
        nodes[2].channel().contains("declared zero"),
        "the wave channel is DECLARED, not omitted"
    );
    assert!(
        matches!(nodes[2].color(), Color::Estimated { dispersion, .. } if dispersion.is_infinite())
    );
    assert!(matches!(
        explanation.receipt().aggregate_color,
        Color::Estimated { .. }
    ));
    println!(
        "{{\"metric\":\"drag-decomposition\",\"induced\":{:.6},\"viscous\":{:.6},\
         \"wave\":0.0,\"residual\":{residual:.2e}}}",
        nodes[0].contribution(),
        nodes[1].contribution()
    );
    verdict(
        "xp-004",
        "flagship: the Trefftz wake integral matches CL^2/(pi AR) within 2%; the \
         measured/heuristic induced, viscous, and wave channels remain Estimated and \
         honestly refuse to certify the analytic/model discrepancy",
    );
}

#[test]
fn xp_005_narrative_is_non_authoritative() {
    let edits = vec![("polish".to_string(), 1.0, 1.1)];
    let explanation = finalize(
        provenance_attribution(&edits).expect("valid history"),
        0.1,
        1e-9,
    )
    .expect("valid finalization inputs");
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

fn assert_xp_006_one_node_fixture() {
    let fixture = elliptic_fixture(1);
    let a0 = vec![1.0, 1.0];
    let a1 = vec![1.2, 0.9];
    let observed = fixture
        .compliance(&fixture.solve(&a1).expect("positive conductivity solves"))
        .expect("finite compliance")
        - fixture
            .compliance(&fixture.solve(&a0).expect("positive conductivity solves"))
            .expect("finite compliance");
    let explanation = finalize(
        adjoint_attribution(&fixture, &a0, &a1, &[("all", vec![0, 1])])
            .expect("valid one-node attribution"),
        observed,
        1e-12,
    )
    .expect("valid finalization inputs");
    assert!(matches!(explanation, Explanation::Refused { .. }));
    assert!(explanation.is_structurally_valid());
    assert!(!explanation.reconciles());
}

fn assert_xp_006_fail_fast_contracts() {
    let fixture = elliptic_fixture(1);
    let a0 = vec![1.0, 1.0];
    let a1 = vec![1.2, 0.9];
    assert!(matches!(
        fixture.solve(&[1.0]),
        Err(ExplanationError::LengthMismatch {
            field: "conductivity",
            expected: 2,
            actual: 1,
        })
    ));
    assert!(matches!(
        adjoint_attribution(&fixture, &a0, &a1, &[("bad", vec![2])]),
        Err(ExplanationError::InvalidIndex {
            field: "adjoint channel element",
            index: 2,
            upper_bound: 2,
        })
    ));
    assert!(matches!(
        ExplanationNode::new(
            "bad",
            1.0,
            -1.0,
            Color::Estimated {
                estimator: "fixture".to_string(),
                dispersion: 0.1,
            },
            vec!["evidence".to_string()],
        ),
        Err(ExplanationError::InvalidNumber { field: "bound", .. })
    ));
    assert!(matches!(
        ExplanationNode::new(
            "bad-color-envelope",
            2.0,
            0.1,
            Color::Verified { lo: 1.0, hi: 1.0 },
            vec!["evidence".to_string()],
        ),
        Err(ExplanationError::InvalidColor { .. })
    ));
    let valid_nodes =
        provenance_attribution(&[("valid".to_string(), 0.0, 1.0)]).expect("valid history");
    assert!(matches!(
        finalize(valid_nodes, f64::NAN, 1e-9),
        Err(ExplanationError::InvalidNumber {
            field: "observed change",
            ..
        })
    ));
    assert!(matches!(
        LiftingLine::elliptic(1.0, 1.0, 1.0, 1.0, 0),
        Err(ExplanationError::InvalidCount {
            field: "lifting-line stations",
            value: 0,
            ..
        })
    ));
}

fn assert_xp_006_color_bound_fingerprint() {
    let verified = ExplanationNode::new(
        "same",
        1.0,
        0.1,
        Color::Verified { lo: 0.9, hi: 1.1 },
        vec!["evidence".to_string()],
    )
    .expect("valid verified node");
    let estimated = ExplanationNode::new(
        "same",
        1.0,
        0.1,
        Color::Estimated {
            estimator: "surrogate".to_string(),
            dispersion: 0.1,
        },
        vec!["evidence".to_string()],
    )
    .expect("valid estimated node");
    assert_ne!(
        verified.fingerprint(),
        estimated.fingerprint(),
        "fingerprints must include evidence color"
    );
}

#[test]
fn xp_006_edge_contracts_fail_fast_and_one_node_solves() {
    assert_xp_006_one_node_fixture();
    assert_xp_006_fail_fast_contracts();
    assert_xp_006_color_bound_fingerprint();
    verdict(
        "xp-006",
        "edge contracts fail fast; the one-node heuristic attribution refuses certification; fingerprints include color",
    );
}

#[test]
fn estimated_bounds_do_not_certify_residual_coverage() {
    let node = ExplanationNode::new(
        "partial-channel",
        0.0,
        1e6,
        Color::Estimated {
            estimator: "coverage-regression".to_string(),
            dispersion: 0.0,
        },
        vec!["test-evidence".to_string()],
    )
    .expect("valid estimated node");
    let outcome = finalize(vec![node], 1e6, 1e6).expect("valid finalization inputs");
    let Explanation::Refused {
        residual,
        ref receipt,
        ..
    } = outcome
    else {
        panic!("an Estimated bound must not discharge certified coverage");
    };
    assert_eq!(residual.to_bits(), 1e6_f64.to_bits());
    assert_eq!(receipt.certified_coverage.to_bits(), 0.0_f64.to_bits());
    assert!(matches!(receipt.aggregate_color, Color::Estimated { .. }));
    assert!(outcome.is_structurally_valid());
    assert!(!outcome.reconciles());
}

#[test]
fn zero_width_verified_node_does_not_receive_a_hidden_tolerance() {
    let node = ExplanationNode::new(
        "exact-channel",
        1.0,
        0.0,
        Color::Verified { lo: 1.0, hi: 1.0 },
        vec!["test-evidence".to_string()],
    )
    .expect("valid exact node");
    let outcome = finalize(vec![node], 1.0_f64.next_up(), 1.0).expect("valid finalization inputs");
    assert!(matches!(outcome, Explanation::Refused { .. }));
    assert_eq!(
        outcome.receipt().aggregation_roundoff.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        outcome.receipt().certified_coverage.to_bits(),
        0.0_f64.to_bits()
    );
    assert!(outcome.is_structurally_valid());
}

#[test]
fn exact_color_bytes_detect_same_rendered_payload_tamper() {
    let first_color = Color::Estimated {
        estimator: "same-render".to_string(),
        dispersion: 0.100_000_01,
    };
    let second_color = Color::Estimated {
        estimator: "same-render".to_string(),
        dispersion: 0.100_000_02,
    };
    assert_eq!(
        first_color.payload_json(),
        second_color.payload_json(),
        "the legacy display encoding demonstrates the collision"
    );
    let first = ExplanationNode::new(
        "color-tamper",
        1.0,
        0.1,
        first_color,
        vec!["test-evidence".to_string()],
    )
    .expect("valid first node");
    let second = ExplanationNode::new(
        "color-tamper",
        1.0,
        0.1,
        second_color,
        vec!["test-evidence".to_string()],
    )
    .expect("valid second node");
    assert_ne!(first.fingerprint(), second.fingerprint());
}

#[test]
fn provenance_fingerprint_binds_the_exact_history() {
    let first = provenance_attribution(&[("same-edit".to_string(), 1.0, 2.0)])
        .expect("valid first history");
    let second = provenance_attribution(&[("same-edit".to_string(), 10.0, 11.0)])
        .expect("valid second history");
    assert_eq!(
        first[0].contribution().to_bits(),
        second[0].contribution().to_bits()
    );
    assert_ne!(first[0].derivation_digest(), second[0].derivation_digest());
    assert_ne!(first[0].fingerprint(), second[0].fingerprint());

    let fixture = elliptic_fixture(1);
    let channels = [("all", vec![0, 1])];
    let zero_at_one = adjoint_attribution(&fixture, &[1.0, 1.0], &[1.0, 1.0], &channels)
        .expect("valid first attribution");
    let zero_at_two = adjoint_attribution(&fixture, &[2.0, 2.0], &[2.0, 2.0], &channels)
        .expect("valid second attribution");
    assert_eq!(
        zero_at_one[0].contribution().to_bits(),
        zero_at_two[0].contribution().to_bits()
    );
    assert_ne!(
        zero_at_one[0].derivation_digest(),
        zero_at_two[0].derivation_digest()
    );
    assert_ne!(zero_at_one[0].fingerprint(), zero_at_two[0].fingerprint());
}

#[test]
fn finalized_receipt_binds_policy_color_and_top_level_payload() {
    let mut input_bound = provenance_attribution(&[("receipt-channel".to_string(), 0.0, 1.0)])
        .expect("valid receipt history");
    let node = input_bound.pop().expect("one input-bound receipt channel");
    let explanation = finalize(vec![node], 1.0, 10.0).expect("valid finalization inputs");
    assert!(explanation.reconciles());
    assert!(explanation.is_structurally_valid());
    assert_eq!(
        explanation.receipt().requested_threshold.to_bits(),
        10.0_f64.to_bits()
    );
    assert_eq!(
        explanation.receipt().aggregation_roundoff.to_bits(),
        0.0_f64.to_bits()
    );
    assert!(matches!(
        explanation.receipt().aggregate_color,
        Color::Estimated { .. }
    ));

    let mut policy_tamper = explanation.clone();
    let Explanation::Explained { receipt, .. } = &mut policy_tamper else {
        unreachable!();
    };
    receipt.requested_threshold = 11.0;
    assert!(!policy_tamper.is_structurally_valid());
    assert!(!policy_tamper.reconciles());

    let mut observed_tamper = explanation;
    let Explanation::Explained { observed, .. } = &mut observed_tamper else {
        unreachable!();
    };
    *observed = observed.next_up();
    assert!(!observed_tamper.is_structurally_valid());

    let ordered_nodes = vec![
        ExplanationNode::new(
            "first",
            0.25,
            0.0,
            Color::Verified { lo: 0.25, hi: 0.25 },
            vec!["first-evidence".to_string()],
        )
        .expect("valid first ordered node"),
        ExplanationNode::new(
            "second",
            0.75,
            0.0,
            Color::Verified { lo: 0.75, hi: 0.75 },
            vec!["second-evidence".to_string()],
        )
        .expect("valid second ordered node"),
    ];
    let mut order_tamper = finalize(ordered_nodes, 1.0, 1.0).expect("valid ordered finalization");
    assert!(order_tamper.is_structurally_valid());
    let Explanation::Explained { nodes, .. } = &mut order_tamper else {
        unreachable!();
    };
    nodes.swap(0, 1);
    assert!(!order_tamper.is_structurally_valid());
}

#[test]
fn validated_color_without_regime_witness_is_rejected() {
    let attempt = ExplanationNode::new(
        "validated-without-witness",
        1.0,
        0.1,
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("mach", 0.0, 0.8),
            dataset: "wind-tunnel-anchor".to_string(),
        },
        vec!["test-evidence".to_string()],
    );
    assert!(matches!(
        attempt,
        Err(ExplanationError::InvalidColor {
            reason: "Validated requires a retained regime-membership witness"
        })
    ));
}

#[test]
fn refused_receipt_is_valid_but_never_reconciled_and_detects_tamper() {
    let node = ExplanationNode::new(
        "partial-channel",
        1.0,
        0.0,
        Color::Verified { lo: 1.0, hi: 1.0 },
        vec!["test-evidence".to_string()],
    )
    .expect("valid partial node");
    let outcome = finalize(vec![node], 2.0, 1.0).expect("valid finalization inputs");
    assert!(matches!(outcome, Explanation::Refused { .. }));
    assert!(outcome.is_structurally_valid());
    assert!(!outcome.reconciles());

    let mut forged = outcome;
    let Explanation::Refused { receipt, .. } = &mut forged else {
        unreachable!();
    };
    receipt.root.push('0');
    assert!(!forged.is_structurally_valid());
}

#[test]
fn node_and_text_identity_tampering_fail_closed() {
    let unretained_verified = ExplanationNode::new(
        "unretained-verified",
        0.0,
        100.0,
        Color::Verified {
            lo: -100.0,
            hi: 100.0,
        },
        vec!["test-evidence".to_string()],
    )
    .expect("valid unretained node");
    let refused =
        finalize(vec![unretained_verified], 100.0, 100.0).expect("valid finalization inputs");
    assert!(matches!(refused, Explanation::Refused { .. }));
    assert_eq!(
        refused.receipt().certified_coverage.to_bits(),
        0.0_f64.to_bits()
    );
    assert!(matches!(
        refused.receipt().aggregate_color,
        Color::Estimated { .. }
    ));

    let fixture = elliptic_fixture(1);
    let mut trusted =
        adjoint_attribution(&fixture, &[1.0, 1.0], &[1.0, 1.0], &[("all", vec![0, 1])])
            .expect("valid trusted attribution");
    let trusted_node = trusted.pop().expect("one trusted channel");
    let duplicate = trusted_node.clone();
    assert!(matches!(
        finalize(vec![trusted_node, duplicate], 0.0, 1.0),
        Err(ExplanationError::DuplicateIdentity { .. })
    ));
    let mut first_batch = adjoint_attribution(
        &fixture,
        &[1.0, 1.0],
        &[1.0, 1.0],
        &[("first-batch", vec![0])],
    )
    .expect("valid first attribution batch");
    let mut second_batch = adjoint_attribution(
        &fixture,
        &[1.0, 1.0],
        &[1.0, 1.0],
        &[("second-batch", vec![0])],
    )
    .expect("valid second attribution batch");
    let first_node = first_batch.pop().expect("one first-batch node");
    let second_node = second_batch.pop().expect("one second-batch node");
    assert!(matches!(
        finalize(vec![first_node, second_node], 0.0, 1.0),
        Err(ExplanationError::IntegrityMismatch {
            field: "built-in attribution batch",
            ..
        })
    ));

    for bad_channel in ["", " ", "line\nbreak"] {
        assert!(matches!(
            ExplanationNode::new(
                bad_channel,
                0.0,
                0.0,
                Color::Verified { lo: 0.0, hi: 0.0 },
                vec!["test-evidence".to_string()],
            ),
            Err(ExplanationError::InvalidText {
                field: "channel",
                ..
            })
        ));
    }
    assert!(matches!(
        ExplanationNode::new(
            "valid-channel",
            0.0,
            0.0,
            Color::Verified { lo: 0.0, hi: 0.0 },
            vec!["bad\nevidence".to_string()],
        ),
        Err(ExplanationError::InvalidText {
            field: "evidence link",
            index: Some(0),
            ..
        })
    ));
    assert!(matches!(
        provenance_attribution(&[(" ".to_string(), 0.0, 1.0)]),
        Err(ExplanationError::InvalidText {
            field: "provenance edit",
            index: Some(0),
            ..
        })
    ));
}

fn assert_provenance_collection_contracts() {
    assert!(matches!(
        provenance_attribution(&[]),
        Err(ExplanationError::InvalidCount {
            field: "provenance edits",
            value: 0,
            ..
        })
    ));
    assert!(matches!(
        provenance_attribution(&[
            ("duplicate".to_string(), 0.0, 1.0),
            ("duplicate".to_string(), 1.0, 2.0),
        ]),
        Err(ExplanationError::DuplicateIdentity {
            field: "provenance edit name"
        })
    ));
    assert!(matches!(
        provenance_attribution(&[("nonfinite".to_string(), 0.0, f64::INFINITY)]),
        Err(ExplanationError::InvalidNumber {
            field: "provenance after state",
            index: Some(0),
            ..
        })
    ));
    assert!(matches!(
        provenance_attribution(&[
            ("first".to_string(), 0.0, 1.0),
            ("second".to_string(), 2.0, 3.0),
        ]),
        Err(ExplanationError::DisconnectedHistory { edit_index: 1 })
    ));
    assert!(matches!(
        provenance_attribution(&[("overflow".to_string(), -f64::MAX, f64::MAX,)]),
        Err(ExplanationError::InvalidNumber {
            field: "provenance contribution",
            index: Some(0),
            ..
        })
    ));
    assert!(matches!(
        provenance_attribution(&[("one-sided-extreme".to_string(), 0.0, f64::MAX)]),
        Err(ExplanationError::InvalidNumber {
            field: "provenance rounding envelope",
            index: Some(0),
            ..
        })
    ));
    let oversized_history = (0..1_025)
        .map(|index| (format!("edit-{index}"), 0.0, 0.0))
        .collect::<Vec<_>>();
    assert!(matches!(
        provenance_attribution(&oversized_history),
        Err(ExplanationError::InvalidCount {
            field: "provenance edits",
            value: 1_025,
            ..
        })
    ));
    let signed_zero_chain = provenance_attribution(&[
        ("first".to_string(), 1.0, -0.0),
        ("second".to_string(), 0.0, 2.0),
    ])
    .expect("signed-zero history telescopes");
    assert_eq!(signed_zero_chain.len(), 2);
}

fn assert_adjoint_and_elliptic_collection_contracts() {
    let fixture = elliptic_fixture(1);
    assert!(matches!(
        Elliptic1d::new(0),
        Err(ExplanationError::InvalidCount {
            field: "Elliptic1d interior nodes",
            value: 0,
            ..
        })
    ));
    assert!(matches!(
        adjoint_attribution(
            &fixture,
            &[1.0, 1.0],
            &[1.0, 1.0],
            &[("left", vec![0]), ("right", vec![0, 1])],
        ),
        Err(ExplanationError::OverlappingChannelElement { element: 0 })
    ));
    assert!(matches!(
        adjoint_attribution(
            &fixture,
            &[1.0, 1.0],
            &[1.0, 1.0],
            &[("same", vec![0]), ("same", vec![1])],
        ),
        Err(ExplanationError::DuplicateIdentity {
            field: "adjoint channel name"
        })
    ));
    assert!(matches!(
        adjoint_attribution(&fixture, &[1.0, 1.0], &[1.0, 1.0], &[("empty", Vec::new())],),
        Err(ExplanationError::InvalidCount {
            field: "adjoint channel mask elements",
            value: 0,
            ..
        })
    ));
    assert!(matches!(
        fixture.solve(&[f64::MAX, f64::MAX]),
        Err(ExplanationError::InvalidNumber {
            field: "assembled stiffness weight",
            ..
        })
    ));
    assert!(matches!(
        fixture.solve(&[1.0, f64::NAN]),
        Err(ExplanationError::InvalidNumber {
            field: "conductivity",
            index: Some(1),
            ..
        })
    ));
    assert!(matches!(
        fixture.compliance(&[]),
        Err(ExplanationError::LengthMismatch {
            field: "elliptic state",
            expected: 1,
            actual: 0,
        })
    ));
    assert!(matches!(
        fixture.compliance(&[f64::INFINITY]),
        Err(ExplanationError::InvalidNumber {
            field: "elliptic state",
            index: Some(0),
            ..
        })
    ));
    assert!(matches!(
        LiftingLine::elliptic(0.1, 8.0, 1.0, 8.0, 4_097),
        Err(ExplanationError::InvalidCount {
            field: "lifting-line stations",
            value: 4_097,
            ..
        })
    ));
    assert!(matches!(
        Elliptic1d::new(65_537),
        Err(ExplanationError::InvalidCount {
            field: "Elliptic1d interior nodes",
            value: 65_537,
            ..
        })
    ));
}

fn assert_drag_and_finalize_collection_contracts() {
    let finite_wing = LiftingLine::elliptic(0.1, 8.0, 1.0, 8.0, 8).expect("valid finite wing");
    assert!(matches!(
        LiftingLine::elliptic(f64::NAN, 8.0, 1.0, 8.0, 8),
        Err(ExplanationError::InvalidNumber {
            field: "elliptic circulation amplitude",
            ..
        })
    ));
    assert!(matches!(
        drag_decomposition(
            &finite_wing,
            f64::MAX,
            2.0,
            0.3,
            "mach-probe:overflow",
            0.0,
            1.0,
        ),
        Err(ExplanationError::InvalidNumber {
            field: "drag diagnostic contribution/dispersion",
            ..
        })
    ));
    assert!(matches!(
        drag_decomposition(&finite_wing, 0.01, 2.0, 0.3, "bad\nevidence", 0.1, 1.0),
        Err(ExplanationError::InvalidText {
            field: "subsonic evidence identity",
            ..
        })
    ));
    let extreme_wing =
        LiftingLine::elliptic(0.1, 1.0, f64::MAX, 1.0, 8).expect("finite inputs construct");
    assert!(matches!(
        extreme_wing.cl(),
        Err(ExplanationError::InvalidNumber {
            field: "lift coefficient",
            ..
        })
    ));

    let oversized_nodes = (0..1_025)
        .map(|index| {
            ExplanationNode::new(
                &format!("channel-{index}"),
                0.0,
                0.0,
                Color::Estimated {
                    estimator: format!("estimate-{index}"),
                    dispersion: 0.0,
                },
                vec![format!("evidence-{index}")],
            )
            .expect("valid bounded node")
        })
        .collect();
    assert!(matches!(
        finalize(oversized_nodes, 0.0, 0.0),
        Err(ExplanationError::InvalidCount {
            field: "explanation nodes",
            value: 1_025,
            ..
        })
    ));
}

#[test]
fn exact_engine_collection_contracts_fail_closed() {
    assert_provenance_collection_contracts();
    assert_adjoint_and_elliptic_collection_contracts();
    assert_drag_and_finalize_collection_contracts();
}

#[test]
fn drag_wave_claim_demotes_outside_declared_subsonic_regime() {
    let wing = LiftingLine::elliptic(0.1, 8.0, 1.0, 8.0, 100).expect("valid wing fixture");
    let cf = 0.006;
    let wetted = 2.0;
    let observed = wing
        .induced_drag_coefficient()
        .expect("finite induced drag")
        + cf * wetted;
    let outcome = drag_decomposition(
        &wing,
        cf,
        wetted,
        0.9,
        "mach-probe:outside-regime",
        observed,
        1.0,
    )
    .expect("valid drag decomposition inputs");
    let Explanation::Explained {
        ref nodes,
        ref receipt,
        ..
    } = outcome
    else {
        panic!("an exactly reconciled but demoted decomposition remains inspectable");
    };
    assert!(matches!(nodes[2].color(), Color::Estimated { .. }));
    assert!(matches!(receipt.aggregate_color, Color::Estimated { .. }));
    assert!(outcome.is_structurally_valid());
}
