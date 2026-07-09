//! End-to-end battery: an illuminated credibility map. Run to steady state,
//! EVERY point matches the analytic Poiseuille solution — so the credibility
//! differentiation is the REGIME certificate (near-τ=0.5 points are flagged as
//! risky even where they happen to be accurate), not a step-budget artifact.

use fs_evidence::Color;
use fs_flowcert_e2e::{certify_point, default_sweep, run_campaign};

#[test]
fn the_credibility_map_is_illuminated_and_certified() {
    let (re, ny) = default_sweep();
    let report = run_campaign(&re, &ny, 100_000, 0.03);
    // ILLUMINATION: a real (Reynolds x resolution) atlas.
    assert!(
        report.num_niches >= 5,
        "too few niches: {}",
        report.num_niches
    );
    assert!(report.coverage > 0.0 && report.qd_score > 0.0);
    // ACCURACY IS NOT THE ISSUE: every point reaches steady state and matches the
    // analytic profile — so inaccuracy is never mistaken for a regime problem.
    assert!(
        report.points.iter().all(|p| p.converged),
        "a point did not reach steady state"
    );
    assert!(
        report.all_accurate,
        "a converged point missed the analytic solution"
    );
    assert!(report.best_error < 0.03, "best error {}", report.best_error);
    // CREDIBILITY MAP: the differentiation is the REGIME certificate — some
    // operating points sit in a comfortably-stable regime, others near τ=0.5 do
    // not, even though all are accurate here.
    let fully_credible = report
        .points
        .iter()
        .filter(|p| p.accurate && p.regime_stable)
        .count();
    assert!(fully_credible > 0, "no fully-credible operating point");
    assert!(
        fully_credible < report.points.len(),
        "the atlas should flag some points by regime"
    );
    assert!(
        report.stable_fraction > 0.0 && report.stable_fraction < 1.0,
        "no regime boundary: {}",
        report.stable_fraction
    );
    // every point carries a resolution and a positive viscosity.
    assert!(
        report
            .points
            .iter()
            .all(|p| p.viscosity > 0.0 && p.tau > 0.5)
    );
    println!(
        "{{\"campaign\":\"flowcert\",\"niches\":{},\"coverage\":{:.3},\"best_error\":{:.4},\
         \"all_accurate\":{},\"stable_fraction\":{:.3},\"points\":{:?}}}",
        report.num_niches,
        report.coverage,
        report.best_error,
        report.all_accurate,
        report.stable_fraction,
        report
            .points
            .iter()
            .map(|p| (p.reynolds, p.ny, p.profile_error, p.regime_stable))
            .collect::<Vec<_>>(),
    );
}

#[test]
fn a_single_low_reynolds_point_is_fully_verified() {
    // low Re, decent resolution: converges, accurate, AND in a stable regime.
    let p = certify_point(20.0, 24, 0.05, 100_000, 0.03);
    assert!(p.converged, "did not reach steady state");
    assert!(p.accurate, "error {}", p.profile_error);
    assert!(p.regime_stable, "regime not stable (tau {})", p.tau);
}

#[test]
fn a_near_tau_half_point_is_flagged_by_regime_even_when_accurate() {
    // high Re => τ near 0.5: once converged it IS accurate, yet the regime
    // certificate flags it as risky (small stability margin).
    let p = certify_point(90.0, 24, 0.05, 100_000, 0.03);
    assert!(p.converged && p.accurate, "error {}", p.profile_error);
    assert!(
        !p.regime_stable,
        "expected a flagged regime (tau {})",
        p.tau
    );
}

#[test]
fn the_campaign_is_deterministic() {
    // A low cap keeps this fast; the run is deterministic at any cap.
    let (re, ny) = default_sweep();
    let a = run_campaign(&re, &ny, 10_000, 0.05);
    let b = run_campaign(&re, &ny, 10_000, 0.05);
    assert_eq!(a.best_error.to_bits(), b.best_error.to_bits());
    assert_eq!(a.num_niches, b.num_niches);
    assert_eq!(
        matches!(a.credibility_color, Color::Verified { .. }),
        matches!(b.credibility_color, Color::Verified { .. })
    );
}
