//! Gradient-certificate conformance (the bk0o.3 bead; runs under the
//! `gradient-certs` feature). Acceptance: every emitted gradient
//! carries an honest color and retains sampled transpose diagnostics; FD
//! spot checks run via the falsifier pairing and agree within
//! conditioning-aware tolerance; a SEEDED TRANSPOSE BUG is caught and
//! BLOCKS merge; color assignment is correct across
//! differentiable/flagged/anchored paths without converting raw inputs into
//! authority; bounded input validation precedes allocation/indexing; the
//! no-falsifier-no-ship gate covers the adjoint-gradient class.
#![cfg(feature = "gradient-certs")]

use std::collections::BTreeSet;

use fs_adjoint::certs::{
    Anchor, GradientCertError, MAX_ADJOINT_RESIDUAL_ENTRY_VISITS, MAX_ADJOINT_RESIDUAL_PROBES,
    MAX_FD_DIRECTIONS, MAX_PROBE_COMPONENTS, SparseLinear, adjoint_residual_bound, certify,
    fd_spot_checks, merge_gate,
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
    let bound = adjoint_residual_bound(&op, 24).expect("bounded valid transpose fixture");
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
    let checks_ok = fd_spot_checks("objective:linear-smoothing", &f, &x, &grad_ok, 4, 0xabcd)
        .expect("finite fixture");
    let checks_bad = fd_spot_checks("objective:linear-smoothing", &f, &x, &grad_bad, 4, 0xabcd)
        .expect("finite fixture");
    let ok_context = checks_ok.context_digest().to_string();
    let bad_context = checks_bad.context_digest().to_string();
    let cert_ok = certify(
        &grade_ops(&["convert/restrict"], &BTreeSet::new()),
        Some(adjoint_residual_bound(&op, 8).expect("bounded valid transpose fixture")),
        Some(checks_ok),
        None,
    );
    let cert_bad = certify(
        &grade_ops(&["convert/restrict"], &BTreeSet::new()),
        Some(adjoint_residual_bound(&op, 8).expect("bounded valid transpose fixture")),
        Some(checks_bad),
        None,
    );
    merge_gate(&cert_ok, &ok_context).expect("the correct gradient merges");
    let refusal = merge_gate(&cert_bad, &bad_context).expect_err("the sign bug must block merge");
    assert!(
        format!("{refusal}").contains("transpose or sign bug"),
        "teaches the cause: {refusal}"
    );
    // And a certificate with NO checks at all is refused too.
    let cert_none = certify(
        &grade_ops(&["convert/restrict"], &BTreeSet::new()),
        None,
        None,
        None,
    );
    let refusal2 = merge_gate(&cert_none, &ok_context).expect_err("no checks, no merge");
    assert!(format!("{refusal2}").contains("mandatory"), "{refusal2}");
    // The COLOR must honor the falsifier, not just the merge gate (bead 9sf6):
    // the refuted gradient cannot be laundered into a Verified ledger color even
    // though its structural transpose-residual diagnostic is tiny. The correct
    // gradient also remains Estimated: the sampled transpose identity does not
    // prove the objective's derivative.
    assert!(
        matches!(cert_ok.color(), Color::Estimated { .. }),
        "a sampled transpose diagnostic cannot mint Verified, got {:?}",
        cert_ok.color()
    );
    assert!(
        matches!(cert_bad.color(), Color::Estimated { .. }),
        "a gradient refuted by its own FD falsifier must NOT wear Verified, got {:?}",
        cert_bad.color()
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
    // A tiny caller-provided residual remains a diagnostic, not a gradient
    // certificate. It is retained without upgrading the color.
    let c1 = certify(&smooth, Some(3e-15), None, None);
    assert!(matches!(c1.color(), Color::Estimated { .. }));
    assert_eq!(c1.residual_bound(), Some(3e-15));
    assert!(c1.discontinuity().is_none());
    // Flagged remesh → ESTIMATED (inherited, never upgraded — that
    // would be laundering).
    let c2 = certify(&flagged, Some(1e-15), None, None);
    assert!(
        matches!(c2.color(), Color::Estimated { .. }),
        "flagged stays estimated even with a residual bound: {:?}",
        c2.color()
    );
    assert!(c2.discontinuity().is_some());
    // Raw caller anchor metadata is retained but cannot authenticate itself.
    let anchor = Anchor {
        dataset: "wind-tunnel-2026-03".to_string(),
        regime: ValidityDomain::unconstrained().with("mach", 0.0, 0.8),
    };
    let c3 = certify(&smooth, Some(1e-15), None, Some(&anchor));
    assert!(matches!(c3.color(), Color::Estimated { .. }));
    assert_eq!(c3.anchor(), Some(&anchor));
    // Smooth with NO evidence at all → ESTIMATED (folklore rule).
    let c4 = certify(&smooth, None, None, None);
    assert!(
        matches!(c4.color(), Color::Estimated { estimator, .. } if estimator.contains("without")),
        "no certificate, no verified color: {:?}",
        c4.color()
    );
    verdict(
        "gc-003",
        "raw residuals and anchors remain Estimated diagnostics; flagged and evidence-free \
         gradients also remain Estimated",
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
    let cert = certify(&grade, Some(f64::INFINITY), None, None);
    assert!(
        !matches!(cert.color(), fs_evidence::Color::Verified { .. }),
        "a non-finite bound must not certify: {:?}",
        cert.color()
    );
    let cert_nan = certify(&grade, Some(f64::NAN), None, None);
    assert!(!matches!(
        cert_nan.color(),
        fs_evidence::Color::Verified { .. }
    ));
    assert!(cert.residual_bound().is_none());
    assert!(cert_nan.residual_bound().is_none());
    // A finite caller value is retained but still cannot verify a gradient.
    let cert_ok = certify(&grade, Some(1e-12), None, None);
    assert!(matches!(
        cert_ok.color(),
        fs_evidence::Color::Estimated { .. }
    ));
    assert_eq!(cert_ok.residual_bound(), Some(1e-12));
}

#[test]
fn malformed_and_unbounded_probe_inputs_return_structured_errors() {
    let op = SparseLinear {
        rows: vec![vec![(0usize, 1.0f64)]],
        ncols: 1,
    };
    assert!(matches!(
        adjoint_residual_bound(&op, 0),
        Err(GradientCertError::InvalidCount { value: 0, .. })
    ));
    assert!(matches!(
        adjoint_residual_bound(&op, MAX_ADJOINT_RESIDUAL_PROBES + 1),
        Err(GradientCertError::InvalidCount { .. })
    ));

    let invalid_column = SparseLinear {
        rows: vec![vec![(1, 1.0)]],
        ncols: 1,
    };
    assert!(matches!(
        adjoint_residual_bound(&invalid_column, 1),
        Err(GradientCertError::InvalidSparseEntry { .. })
    ));
    let non_finite_weight = SparseLinear {
        rows: vec![vec![(0, f64::NAN)]],
        ncols: 1,
    };
    assert!(matches!(
        adjoint_residual_bound(&non_finite_weight, 1),
        Err(GradientCertError::InvalidSparseEntry { .. })
    ));

    // Entry and probe counts can each satisfy their individual caps while
    // their product would still trigger billions of interval operations. The
    // aggregate two-apply visit budget must refuse before probe allocation or
    // execution.
    let entries = MAX_ADJOINT_RESIDUAL_ENTRY_VISITS / (2 * MAX_ADJOINT_RESIDUAL_PROBES) + 1;
    let dense_single_row = SparseLinear {
        rows: vec![vec![(0, 1.0); entries]],
        ncols: 1,
    };
    let expected_visits = entries * 2 * MAX_ADJOINT_RESIDUAL_PROBES;
    assert!(matches!(
        adjoint_residual_bound(&dense_single_row, MAX_ADJOINT_RESIDUAL_PROBES),
        Err(GradientCertError::WorkBudgetExceeded {
            operation: "adjoint residual sparse-entry visits",
            requested,
            max: MAX_ADJOINT_RESIDUAL_ENTRY_VISITS,
        }) if requested == expected_visits
    ));

    let f = |values: &[f64]| values.iter().sum::<f64>();
    assert!(matches!(
        fd_spot_checks("objective:resource-check", &f, &[1.0], &[1.0], 0, 1),
        Err(GradientCertError::InvalidCount { value: 0, .. })
    ));
    assert!(matches!(
        fd_spot_checks(
            "objective:resource-check",
            &f,
            &[1.0],
            &[1.0],
            MAX_FD_DIRECTIONS + 1,
            1,
        ),
        Err(GradientCertError::InvalidCount { .. })
    ));
    assert!(matches!(
        fd_spot_checks("objective:resource-check", &f, &[1.0, 2.0], &[1.0], 1, 1,),
        Err(GradientCertError::DimensionMismatch { .. })
    ));
    let dimension = 1_024;
    let directions = MAX_PROBE_COMPONENTS / dimension + 1;
    assert!(matches!(
        fd_spot_checks(
            "objective:resource-check",
            &f,
            &vec![0.0; dimension],
            &vec![0.0; dimension],
            directions,
            1,
        ),
        Err(GradientCertError::WorkBudgetExceeded { .. })
    ));
}

#[test]
fn non_finite_objective_response_is_an_actionable_fd_error() {
    let f = |_values: &[f64]| f64::NAN;
    assert!(matches!(
        fd_spot_checks("objective:nonfinite", &f, &[1.0], &[1.0], 1, 1),
        Err(GradientCertError::InvalidFdVerdict { direction: 0, .. })
    ));
}

#[test]
fn fd_steps_that_round_away_at_large_coordinates_are_refused() {
    let constant = |_values: &[f64]| 1.0;
    let error = fd_spot_checks(
        "objective:large-coordinate",
        &constant,
        &[1e100],
        &[0.0],
        1,
        2,
    )
    .expect_err("a zero-difference batch at an unrepresentable step is not evidence");
    assert!(matches!(
        error,
        GradientCertError::UnrepresentablePerturbation {
            direction: 0,
            component: 0,
            reason: "perturbation rounded back to the unperturbed coordinate",
            ..
        }
    ));

    assert!(matches!(
        fd_spot_checks(
            "objective:nonfinite-coordinate",
            &constant,
            &[f64::INFINITY],
            &[0.0],
            1,
            2,
        ),
        Err(GradientCertError::InvalidVectorValue {
            field: "point",
            index: 0,
            ..
        })
    ));
}

#[test]
fn sealed_fd_context_detects_mutation_and_cross_context_replay() {
    use fs_adjoint::mitigate::GradientGrade;

    let f = |values: &[f64]| values.iter().map(|value| value * value).sum::<f64>();
    let first = fd_spot_checks(
        "objective:quadratic-energy",
        &f,
        &[1.0, 2.0],
        &[2.0, 4.0],
        4,
        7,
    )
    .expect("first bounded context");
    let second = fd_spot_checks(
        "objective:quadratic-energy",
        &f,
        &[1.0, 3.0],
        &[2.0, 6.0],
        4,
        7,
    )
    .expect("mutated bounded context");
    assert_eq!(first.objective_identity(), "objective:quadratic-energy");
    assert_ne!(
        first.context_digest(),
        second.context_digest(),
        "exact x/gradient bits are bound into the sealed batch"
    );

    let expected_first = first.context_digest().to_string();
    let expected_second = second.context_digest().to_string();
    let cert = certify(
        &GradientGrade::Smooth { route: vec![] },
        None,
        Some(first),
        None,
    );
    merge_gate(&cert, &expected_first).expect("matching orchestrator context passes");
    let replay = merge_gate(&cert, &expected_second).expect_err("cross-context replay refuses");
    assert!(replay.to_string().contains("context mismatch"));

    assert!(matches!(
        fd_spot_checks("", &f, &[1.0], &[2.0], 1, 7),
        Err(GradientCertError::InvalidObjectiveIdentity { .. })
    ));
    assert!(matches!(
        fd_spot_checks("derived:v2:forged-objective", &f, &[1.0], &[2.0], 1, 7,),
        Err(GradientCertError::InvalidObjectiveIdentity { .. })
    ));
}
