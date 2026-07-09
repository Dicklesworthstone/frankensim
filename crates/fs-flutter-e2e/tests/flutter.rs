//! End-to-end battery: two independent methods agree on the proven flutter
//! boundary μ*=2, and only the Aitken partitioned solver reaches it.

use fs_evidence::Color;
use fs_flutter_e2e::run_campaign;

#[test]
fn the_flutter_boundary_is_proven_and_cross_checked() {
    let report = run_campaign(0.55, 2.45, 20);
    // the Lyapunov certificate recovers the exact boundary μ* = 2 (last stable
    // sample sits just below 2).
    assert!(
        report.lyapunov_boundary < 2.0 && report.lyapunov_boundary > 1.7,
        "lyapunov boundary {}",
        report.lyapunov_boundary
    );
    // an INDEPENDENT spectral test agrees — certifying the certifier.
    assert!(
        report.boundaries_agree,
        "lyapunov {} vs spectral {}",
        report.lyapunov_boundary, report.spectral_boundary
    );
    // the naive partitioned solver quits early; only Aitken reaches the boundary.
    assert!(
        report.naive_boundary < 1.05,
        "naive boundary {}",
        report.naive_boundary
    );
    assert!(
        report.aitken_beats_naive,
        "aitken {} vs naive {}",
        report.aitken_boundary, report.naive_boundary
    );
    // a witness in the certified-stable range that ONLY Aitken can compute,
    // carrying a Verified Lyapunov certificate: past the naive solver's reach
    // yet strictly below the proven flutter boundary.
    let mu = report
        .witness_mu
        .expect("a stable-but-naive-fails witness exists");
    assert!(
        mu > report.naive_boundary && mu < report.lyapunov_boundary,
        "witness μ {mu} not between naive {} and lyapunov {}",
        report.naive_boundary,
        report.lyapunov_boundary
    );
    assert!(matches!(report.witness_color, Some(Color::Verified { .. })));
    println!(
        "{{\"campaign\":\"fluttercert\",\"lyapunov_boundary\":{:.3},\"spectral_boundary\":{:.3},\
         \"naive_boundary\":{:.3},\"aitken_boundary\":{:.3},\"witness_mu\":{:.3}}}",
        report.lyapunov_boundary,
        report.spectral_boundary,
        report.naive_boundary,
        report.aitken_boundary,
        mu,
    );
}

#[test]
fn beyond_the_boundary_is_not_certified() {
    let report = run_campaign(2.1, 3.0, 10);
    // every sample is past μ* = 2 → no Lyapunov certificate anywhere.
    assert!(report.samples.iter().all(|s| !s.lyapunov_stable));
    assert!(report.samples.iter().all(|s| s.abscissa >= 0.0));
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(0.55, 2.45, 20);
    let b = run_campaign(0.55, 2.45, 20);
    assert_eq!(a.lyapunov_boundary.to_bits(), b.lyapunov_boundary.to_bits());
    assert_eq!(
        a.witness_mu.map(f64::to_bits),
        b.witness_mu.map(f64::to_bits)
    );
}
