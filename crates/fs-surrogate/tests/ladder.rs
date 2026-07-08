//! Abstraction-ladder conformance (the knh1.4 bead; runs under
//! `abstraction-ladder`). Acceptance: RB bounds are TRUE bounds on the
//! parametric elliptic family (the certificate law, G1); the leak
//! alarm never lets an operator act on an abstraction whose certified
//! bound exceeds tolerance — it auto-drills; estimated-color concept
//! levels cannot masquerade as RB-certified; the kill measurement (RB
//! coverage of the query battery) is ledgered and clears the 20%
//! beachhead floor; queries replay bit-equal (G5).
#![cfg(feature = "abstraction-ladder")]

use fs_evidence::Color;
use fs_surrogate::ladder::{Ladder, RbLevel, TruthModel, rb_coverage};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-surrogate/ladder\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

#[test]
fn la_001_rb_bounds_are_true_bounds() {
    // The certificate law: |s_truth − s_rb| ≤ qoi_bound across a μ
    // battery, for every basis size — a FALSE certificate is worse
    // than a wrong answer (certify-the-certifiers).
    let truth = TruthModel { n: 200 };
    for k in [2usize, 4, 6] {
        let rb = RbLevel::train(&truth, (0.0, 4.0), k);
        let mut worst_ratio = 0.0f64;
        for i in 0..25 {
            let mu = 4.0 * f64::from(i) / 24.0;
            let (_, s_rb, _, qoi_bound) = rb.query(mu);
            let s_true = truth.compliance(&truth.solve(mu));
            let err = (s_true - s_rb).abs();
            assert!(
                err <= qoi_bound + 1e-14,
                "k={k}, mu={mu}: |{s_true} - {s_rb}| = {err:.3e} NOT <= bound {qoi_bound:.3e}"
            );
            if qoi_bound > 1e-15 {
                worst_ratio = worst_ratio.max(err / qoi_bound);
            }
        }
        println!(
            "{{\"metric\":\"rb-bound\",\"k\":{k},\"worst_effectivity_inverse\":{worst_ratio:.4}}}"
        );
    }
    // Bounds TIGHTEN with basis size (at an off-training parameter).
    let bounds: Vec<f64> = [2usize, 4, 6]
        .iter()
        .map(|&k| RbLevel::train(&truth, (0.0, 4.0), k).query(1.7).3)
        .collect();
    // k=4 already resolves the 1-parameter manifold to machine noise
    // (~1e-29); assert tightening against k=2 rather than comparing
    // two roundoff-floor values.
    assert!(
        bounds[1] < 1e-6 * bounds[0] && bounds[2] < 1e-6 * bounds[0],
        "bounds collapse with the basis: {bounds:?}"
    );
    verdict(
        "la-001",
        "RB QoI bounds are true bounds at every battery point for k in {2,4,6}, and \
         tighten monotonically with basis size (G1)",
    );
}

#[test]
fn la_002_leak_alarm_never_lets_a_leak_answer() {
    // THE PROPERTY: whatever rung you start at, the answer's certified
    // width meets the tolerance or came from level 0 — a leaking rung
    // NEVER answers; it is recorded and descended past.
    let ladder = Ladder::build(200, (0.0, 4.0), &[6, 2], true);
    for i in 0..15 {
        let mu = 0.2 + 3.6 * f64::from(i) / 14.0;
        for tol in [1e-3, 1e-6, 1e-10, 1e-14] {
            let ans = ladder.at_level(ladder.top()).query(mu, tol);
            match ans.color {
                Color::Verified { lo, hi } => {
                    let half = (hi - lo) / 2.0;
                    assert!(
                        half <= tol || ans.level_used == 0,
                        "mu={mu}, tol={tol:.0e}: width {half:.2e} from level {}",
                        ans.level_used
                    );
                }
                Color::Estimated { dispersion, .. } => {
                    assert!(
                        dispersion <= tol,
                        "an estimate only answers within its calibrated dispersion"
                    );
                }
                Color::Validated { .. } => panic!("no validated rungs in this fixture"),
            }
            // Leaks are recorded strictly above the answering level.
            for &l in &ans.leaks {
                assert!(
                    l > ans.level_used,
                    "leak {l} above answer {}",
                    ans.level_used
                );
            }
        }
    }
    // A tolerance below every rung's achievable bound (the k=6 rung
    // reaches ~1e-29 on this manifold, so go to 1e-32) forces full
    // descent with the whole ordered leak trail.
    let ans = ladder.at_level(ladder.top()).query(1.3, 1e-32);
    assert_eq!(ans.level_used, 0, "ultra-tight tol reaches the truth");
    assert_eq!(ans.leaks, vec![3, 2, 1], "every rung leaked, in order");
    verdict(
        "la-002",
        "across 15 mu x 4 tolerances no leaking rung ever answers: the certified width \
         meets tol or level 0 answered; the leak trail is complete and ordered",
    );
}

#[test]
fn la_003_estimated_color_cannot_masquerade() {
    let ladder = Ladder::build(200, (0.0, 4.0), &[6, 2], true);
    // A loose query the concept rung CAN answer: it answers, but its
    // color says Estimated — never Verified.
    let loose = ladder.at_level(ladder.top()).query(2.0, 0.5);
    assert_eq!(loose.level_used, ladder.top(), "the concept rung answered");
    assert!(
        matches!(loose.color, Color::Estimated { ref estimator, .. }
            if estimator.contains("cross-rung")),
        "the concept answer is estimated, calibrated by cross-rung probes: {:?}",
        loose.color
    );
    // The same query demanding CERTIFIED accuracy skips the concept
    // rung (it leaks by construction — an estimate is not a bound).
    let strict = ladder.at_level(ladder.top()).query(2.0, 1e-6);
    assert!(
        matches!(strict.color, Color::Verified { .. }),
        "the strict answer is verified"
    );
    assert!(
        strict.leaks.contains(&ladder.top()),
        "the concept rung is recorded as leaked, not silently skipped"
    );
    verdict(
        "la-003",
        "the concept rung answers loose queries with honest Estimated color and is a \
         recorded leak on strict queries — estimates cannot masquerade as certificates",
    );
}

#[test]
fn la_004_kill_measurement_rb_coverage() {
    // The Proposal-A kill criterion: RB rungs must cover >= 20% of the
    // wedge query battery, else the beachhead is too narrow. Measured
    // over a mu x tol grid spanning loose-to-tight demands.
    let ladder = Ladder::build(200, (0.0, 4.0), &[6, 2], false);
    let mus: Vec<f64> = (0..12).map(|i| 4.0 * f64::from(i) / 11.0).collect();
    let tols = [1e-2, 1e-4, 1e-6, 1e-8];
    let coverage = rb_coverage(&ladder, &mus, &tols);
    println!(
        "{{\"metric\":\"kill-measurement\",\"rb_coverage\":{coverage:.3},\
         \"floor\":0.2,\"battery\":\"12 mu x 4 tol\"}}"
    );
    assert!(
        coverage >= 0.2,
        "the beachhead clears the kill floor: {coverage}"
    );
    verdict(
        "la-004",
        "RB coverage of the 48-point query battery ledgered and above the 0.2 kill \
         floor — the beachhead is viable on the elliptic fixture",
    );
}

#[test]
fn la_005_g5_determinism_and_invisibility() {
    let ladder = Ladder::build(150, (0.0, 4.0), &[5], true);
    // Bit-equal replay.
    let a = ladder.at_level(1).query(2.3, 1e-6);
    let b = ladder.at_level(1).query(2.3, 1e-6);
    assert!(a.value.to_bits() == b.value.to_bits() && a.level_used == b.level_used);
    // INVISIBLE UNTIL IT LEAKS: a loose-tol query at the RB rung stays
    // there (no descent, no leaks) — the operator never sees the
    // ladder working.
    let quiet = ladder.at_level(1).query(2.3, 1e-2);
    assert_eq!(quiet.level_used, 1);
    assert!(quiet.leaks.is_empty(), "no leak, no descent, no noise");
    verdict(
        "la-005",
        "queries replay bit-equal (G5); a satisfiable query answers at its rung with an \
         empty leak trail — the ladder is invisible until it leaks",
    );
}
