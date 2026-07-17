//! Battery for surrogates-with-guarantees (fs-surrogate). Covers the POD ROM
//! (exact low-rank reproduction, orthonormal modes, energy-based rank + reduced
//! error), distribution-free conformal bands with empirical coverage, and the
//! certify-or-escalate policy (including cost reduction vs all-high-fidelity).

use fs_surrogate::{
    ConformalBand, Decision, SurrogateError, certify_or_escalate, conformal_band,
    empirical_coverage, pod,
};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn pod_reproduces_a_low_rank_snapshot_set_exactly() {
    // snapshots spanning a 2-D subspace of R^4.
    let snaps = vec![
        vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0],
        vec![1.0, 1.0, 0.0, 0.0],
        vec![2.0, 1.0, 0.0, 0.0],
    ];
    let rom = pod(&snaps, 0.999).unwrap();
    assert!(rom.rank() <= 2);
    assert!(rom.energy_captured() >= 0.999);
    // Every snapshot lies in the reduced space, so POD reconstructs it TO
    // ROUNDOFF, not merely "closely": worst error measured 2.2e-16 over the four
    // snapshots (the Jacobi eigensolver is IEEE-only — sqrt/±/×/÷, no libm — so
    // this is cross-ISA bit-stable). Gate at 1e-12 (~4500× the roundoff floor),
    // 1000× tighter than the old 1e-9, which sat far above any real projection
    // defect while the test still promised the set is reproduced "exactly".
    for s in &snaps {
        assert!(
            rom.reconstruction_error(s) < 1e-12,
            "err {}",
            rom.reconstruction_error(s)
        );
    }
}

#[test]
fn pod_modes_are_orthonormal() {
    let snaps = vec![
        vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0],
        vec![1.0, 1.0, 0.0, 0.0],
        vec![-1.0, 2.0, 0.0, 0.0],
    ];
    let rom = pod(&snaps, 0.999).unwrap();
    assert_eq!(rom.rank(), 2);
    // reconstruct the modes by projecting the unit-energy directions is awkward;
    // instead check orthonormality via the reduced-coordinate identity: the
    // reduced coordinates of a mode-reconstructed vector are canonical.
    let e0 = rom.reconstruct(&[1.0, 0.0]);
    let e1 = rom.reconstruct(&[0.0, 1.0]);
    let mode0: Vec<f64> = e0
        .iter()
        .zip(rom.reconstruct(&[0.0, 0.0]))
        .map(|(a, m)| a - m)
        .collect();
    let mode1: Vec<f64> = e1
        .iter()
        .zip(rom.reconstruct(&[0.0, 0.0]))
        .map(|(a, m)| a - m)
        .collect();
    assert!((dot(&mode0, &mode0) - 1.0).abs() < 1e-9);
    assert!((dot(&mode1, &mode1) - 1.0).abs() < 1e-9);
    assert!(dot(&mode0, &mode1).abs() < 1e-9);
}

#[test]
fn pod_rank_captures_the_dominant_energy() {
    // dominant variation in dim 0 (±10), some in dim 1 (±1), tiny in dim 2.
    let snaps = vec![
        vec![10.0, 1.0, 0.01, 0.0],
        vec![-10.0, 1.0, -0.01, 0.0],
        vec![10.0, -1.0, 0.01, 0.0],
        vec![-10.0, -1.0, -0.01, 0.0],
    ];
    let rom = pod(&snaps, 0.999).unwrap();
    assert!(rom.rank() <= 2, "rank {}", rom.rank());
    assert!(rom.energy_captured() >= 0.999);
    // the tiny dim-2 tail is the only reduced-vs-full error.
    assert!(rom.reconstruction_error(&snaps[0]) < 0.1);
}

#[test]
fn pod_rejects_bad_input() {
    assert_eq!(pod(&[], 0.9), Err(SurrogateError::NoSnapshots));
    assert!(matches!(
        pod(&[vec![1.0, 2.0], vec![1.0]], 0.9),
        Err(SurrogateError::DimMismatch { .. })
    ));
    assert_eq!(pod(&[vec![1.0]], 1.5), Err(SurrogateError::BadThreshold));
}

#[test]
fn the_conformal_band_achieves_its_nominal_coverage() {
    // calibration residuals uniformly on (0, 1].
    let calib: Vec<f64> = (1..=20).map(|i| f64::from(i) / 20.0).collect();
    let band = conformal_band(&calib, 0.1);
    // held-out pairs from the same distribution (residual = |pred - truth|).
    let held_out: Vec<(f64, f64)> = (1..=20).map(|i| (0.0, f64::from(i) / 20.0)).collect();
    let coverage = empirical_coverage(&band, &held_out);
    assert!(coverage >= 0.9, "coverage {coverage} below nominal 0.9");
    assert!(band.covers(0.0, 0.5) && !band.covers(0.0, 100.0));
}

#[test]
fn conformal_band_fails_closed_when_alpha_below_one_over_n_plus_one() {
    // Distribution-free (1-alpha) coverage needs rank ceil((1-alpha)(n+1)) <= n,
    // i.e. alpha >= 1/(n+1). Below that the honest band is unbounded: the old
    // `.clamp(1, n)` instead returned the MAX residual (rank n), whose true
    // coverage is only n/(n+1) < 1-alpha -- a silent under-coverage on a crate
    // that advertises "at least (1-alpha)" coverage. It must fail closed to +inf.
    let residuals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]; // n = 8
    // alpha = 0.1 < 1/9: uncertifiable -> infinite band -> escalate.
    let band = conformal_band(&residuals, 0.1);
    assert!(
        band.half_width.is_infinite(),
        "n=8, alpha=0.1 is below 1/(n+1); the honest band is infinite"
    );
    assert!(
        matches!(
            certify_or_escalate(&band, true, 0.5),
            Decision::Escalate { .. }
        ),
        "an uncertifiable (infinite) band must escalate, never use the surrogate"
    );
    // alpha = 0.15 >= 1/9: certifiable; rank = ceil(0.85*9) = 8 -> 8th residual.
    let ok = conformal_band(&residuals, 0.15);
    assert!(
        (ok.half_width - 0.8).abs() < 1e-12,
        "n=8, alpha=0.15 must return the 8th-smallest residual (0.8), got {}",
        ok.half_width
    );
}

#[test]
fn certify_or_escalate_uses_the_surrogate_only_when_trustworthy() {
    // n = 20 (not 8): at alpha = 0.1 a split-conformal band needs n >= 9
    // (alpha >= 1/(n+1)) to certify (1-alpha) coverage at all; with n = 8 the
    // honest band is +inf. These fixtures exercise the DECISION logic, so use a
    // certifiable calibration size and let the residual magnitude drive it.
    let narrow = conformal_band(&[0.001; 20], 0.1);
    let wide = conformal_band(&[0.5; 20], 0.1);
    // narrow band, in domain, tight enough -> use the surrogate.
    assert!(matches!(
        certify_or_escalate(&narrow, true, 0.1),
        Decision::UseSurrogate { .. }
    ));
    // band too wide for the decision -> escalate.
    assert!(matches!(
        certify_or_escalate(&wide, true, 0.1),
        Decision::Escalate { .. }
    ));
    // outside the validity domain -> escalate regardless.
    assert!(matches!(
        certify_or_escalate(&narrow, false, 0.1),
        Decision::Escalate { .. }
    ));
}

#[test]
fn certify_or_escalate_rejects_unbounded_and_malformed_numeric_policy() {
    // Below 1/(n+1), conformal_band deliberately returns the only honest band:
    // an unbounded one. IEEE-754 makes +inf <= +inf true, so validating the
    // decision boundary is essential when a caller supplies an infinite
    // tolerance. The old comparison-only policy incorrectly used the surrogate.
    let unbounded = conformal_band(&[0.1; 8], 0.1);
    assert_eq!(unbounded.half_width, f64::INFINITY);
    assert!(matches!(
        certify_or_escalate(&unbounded, true, f64::INFINITY),
        Decision::Escalate { .. }
    ));

    // ConformalBand fields are public, so the policy must reject malformed
    // caller-constructed state rather than relying only on conformal_band.
    let malformed_bands = [
        (
            ConformalBand {
                half_width: f64::NAN,
                alpha: 0.1,
            },
            "NaN half-width",
        ),
        (
            ConformalBand {
                half_width: f64::NEG_INFINITY,
                alpha: 0.1,
            },
            "negative-infinite half-width",
        ),
        (
            ConformalBand {
                half_width: -0.1,
                alpha: 0.1,
            },
            "negative half-width",
        ),
        (
            ConformalBand {
                half_width: 0.1,
                alpha: f64::NAN,
            },
            "NaN alpha",
        ),
        (
            ConformalBand {
                half_width: 0.1,
                alpha: 0.0,
            },
            "zero alpha",
        ),
        (
            ConformalBand {
                half_width: 0.1,
                alpha: 1.0,
            },
            "unit alpha",
        ),
    ];
    for (band, label) in malformed_bands {
        assert!(
            matches!(
                certify_or_escalate(&band, true, 0.1),
                Decision::Escalate { .. }
            ),
            "{label} must fail closed"
        );
    }

    let valid_band = ConformalBand {
        half_width: 0.1,
        alpha: 0.1,
    };
    for invalid_tolerance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -f64::EPSILON] {
        assert!(matches!(
            certify_or_escalate(&valid_band, true, invalid_tolerance),
            Decision::Escalate { .. }
        ));
    }

    // The finite inclusive boundary remains intentional: a band exactly as
    // tight as the decision tolerance is admissible.
    assert!(matches!(
        certify_or_escalate(&valid_band, true, 0.1),
        Decision::UseSurrogate {
            band_half_width: 0.1
        }
    ));
    let exact_band = ConformalBand {
        half_width: 0.0,
        alpha: 0.1,
    };
    assert!(matches!(
        certify_or_escalate(&exact_band, true, 0.0),
        Decision::UseSurrogate {
            band_half_width: 0.0
        }
    ));
}

#[test]
fn the_policy_reduces_cost_versus_all_high_fidelity() {
    // n = 20 (not 8): at alpha = 0.1 a split-conformal band needs n >= 9
    // (alpha >= 1/(n+1)) to certify (1-alpha) coverage at all; with n = 8 the
    // honest band is +inf. These fixtures exercise the DECISION logic, so use a
    // certifiable calibration size and let the residual magnitude drive it.
    let narrow = conformal_band(&[0.001; 20], 0.1);
    let wide = conformal_band(&[0.5; 20], 0.1);
    let (surrogate_cost, full_cost) = (1.0, 100.0);
    // a fleet of queries: most trustworthy, a couple must escalate.
    let queries = [
        (narrow, true),
        (narrow, true),
        (narrow, true),
        (wide, true),
        (narrow, false),
    ];
    let total: f64 = queries
        .iter()
        .map(
            |(band, in_dom)| match certify_or_escalate(band, *in_dom, 0.1) {
                Decision::UseSurrogate { .. } => surrogate_cost,
                Decision::Escalate { .. } => full_cost,
            },
        )
        .sum();
    // 3 surrogate + 2 escalate = 203, far below 5 full solves = 500.
    assert!(total < queries.len() as f64 * full_cost);
    assert!((total - 203.0).abs() < 1e-9);
}

#[test]
fn pod_is_deterministic() {
    let snaps = vec![
        vec![1.0, 0.5, 0.0],
        vec![0.0, 1.0, 0.2],
        vec![2.0, 0.3, 0.1],
    ];
    assert_eq!(pod(&snaps, 0.99), pod(&snaps, 0.99));
}
