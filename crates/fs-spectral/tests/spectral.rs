//! Battery for spectral health monitoring (addendum Proposal 5). Covers the
//! Jacobi eigensolver (known spectra + error paths), the λ-gap ratio, the
//! hysteresis health monitor, mandatory low-confidence propagation (demote,
//! never promote), end-to-end conditioning composition, and the router's
//! conditioning-aware path choice.

use fs_evidence::{Color, ColorRank};
use fs_spectral::{
    GapHealthMonitor, Health, RouterPath, SpectralError, compose_conditioning, propagate, route,
    spectral_gap, symmetric_eigenvalues,
};

fn approx_sorted(got: &[f64], want: &[f64], tol: f64) {
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() < tol, "got {got:?}, want {want:?}");
    }
}

#[test]
fn jacobi_recovers_known_spectra() {
    // 2x2: [[2,1],[1,2]] -> {1, 3}.
    approx_sorted(
        &symmetric_eigenvalues(&[vec![2.0, 1.0], vec![1.0, 2.0]]).unwrap(),
        &[1.0, 3.0],
        1e-9,
    );
    // [[0,1],[1,0]] -> {-1, 1}.
    approx_sorted(
        &symmetric_eigenvalues(&[vec![0.0, 1.0], vec![1.0, 0.0]]).unwrap(),
        &[-1.0, 1.0],
        1e-9,
    );
    // diagonal is its own spectrum (sorted).
    approx_sorted(
        &symmetric_eigenvalues(&[
            vec![3.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ])
        .unwrap(),
        &[1.0, 2.0, 3.0],
        1e-12,
    );
    // tridiagonal Toeplitz [[4,1,0],[1,4,1],[0,1,4]] -> 4 + √2·{-1,0,1}.
    let s = 2.0_f64.sqrt();
    approx_sorted(
        &symmetric_eigenvalues(&[
            vec![4.0, 1.0, 0.0],
            vec![1.0, 4.0, 1.0],
            vec![0.0, 1.0, 4.0],
        ])
        .unwrap(),
        &[4.0 - s, 4.0, 4.0 + s],
        1e-9,
    );
}

#[test]
fn eigensolver_rejects_malformed_matrices() {
    assert_eq!(symmetric_eigenvalues(&[]), Err(SpectralError::Empty));
    assert_eq!(
        symmetric_eigenvalues(&[vec![1.0, 2.0]]),
        Err(SpectralError::NotSquare)
    );
    assert_eq!(
        symmetric_eigenvalues(&[vec![1.0, 2.0], vec![3.0, 4.0]]),
        Err(SpectralError::NotSymmetric)
    );
}

#[test]
fn the_gap_ratio_reflects_separation() {
    // well separated: eigenvalues 1, 5, 9 -> gap 4, spread 8, ratio 0.5.
    let g = spectral_gap(&[1.0, 5.0, 9.0]).unwrap();
    assert!((g.gap - 4.0).abs() < 1e-12);
    assert!((g.ratio - 0.5).abs() < 1e-12);
    // near-degenerate lowest pair -> tiny ratio (a collapsing gap).
    let d = spectral_gap(&[1.0, 1.001, 9.0]).unwrap();
    assert!(d.ratio < 0.001);
    // fewer than two eigenvalues -> no gap.
    assert!(spectral_gap(&[1.0]).is_none());
}

#[test]
fn the_health_monitor_has_hysteresis() {
    let mut mon = GapHealthMonitor::new(0.1, 0.3);
    assert_eq!(mon.health(), Health::Healthy);
    // a collapse below the lower threshold degrades.
    assert_eq!(mon.update(0.05), Health::Degraded);
    // a partial recovery INTO the band does not flap back.
    assert_eq!(mon.update(0.2), Health::Degraded);
    // recovery above the upper threshold restores.
    assert_eq!(mon.update(0.35), Health::Healthy);
}

#[test]
#[should_panic(expected = "inverted hysteresis band")]
fn the_monitor_rejects_an_inverted_band() {
    let _ = GapHealthMonitor::new(0.3, 0.1);
}

#[test]
fn a_degraded_gap_propagates_low_confidence() {
    let verified = Color::Verified { lo: -1.0, hi: 1.0 };
    // healthy leaves the color untouched.
    assert_eq!(propagate(verified.clone(), Health::Healthy), verified);
    // degraded DEMOTES verified to estimated (never silently trusted).
    assert_eq!(
        propagate(verified, Health::Degraded).rank(),
        ColorRank::Estimated
    );
    // and it never PROMOTES: an estimated claim stays estimated.
    let est = Color::Estimated {
        estimator: "s".into(),
        dispersion: 1.0,
    };
    assert_eq!(
        propagate(est, Health::Degraded).rank(),
        ColorRank::Estimated
    );
}

#[test]
fn conditioning_composes_multiplicatively() {
    assert!((compose_conditioning(&[2.0, 3.0, 1.5]).unwrap() - 9.0).abs() < 1e-12);
    // an empty pipeline is perfectly conditioned.
    assert!((compose_conditioning(&[]).unwrap() - 1.0).abs() < 1e-12);
    // a bad amplification factor is rejected.
    assert!(compose_conditioning(&[2.0, -1.0]).is_err());
    assert!(compose_conditioning(&[f64::INFINITY]).is_err());
}

#[test]
fn the_router_trades_cheapness_for_conditioning() {
    let paths = vec![
        RouterPath {
            label: "cheap-illcond".into(),
            base_cost: 1.0,
            conditioning: 1e6,
        },
        RouterPath {
            label: "costly-wellcond".into(),
            base_cost: 3.0,
            conditioning: 1.0,
        },
    ];
    // no conditioning weight -> pick the cheapest.
    assert_eq!(route(&paths, 0.0).unwrap().label, "cheap-illcond");
    // a strong conditioning weight -> prefer the well-posed path.
    assert_eq!(route(&paths, 2.0).unwrap().label, "costly-wellcond");
    assert!(route(&[], 1.0).is_none());
}

#[test]
fn the_eigensolver_is_deterministic() {
    let a = vec![
        vec![4.0, 1.0, 0.0],
        vec![1.0, 4.0, 1.0],
        vec![0.0, 1.0, 4.0],
    ];
    assert_eq!(symmetric_eigenvalues(&a), symmetric_eigenvalues(&a));
}
