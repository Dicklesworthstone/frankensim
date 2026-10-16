//! End-to-end battery: two independent methods agree on the proven flutter
//! boundary μ*=2, and only the Aitken partitioned solver reaches it.

use fs_evidence::Color;
use fs_flutter_e2e::{run_campaign, spectral_abscissa, spectral_abscissa_interval};
use fs_math::dd::Dd;

/// Absolute slop allowed against the double-double oracle. `fs_math::dd`
/// documents ≤ 2⁻¹⁰³ relative error for sqrt on finite operands, and every
/// quantity here is O(1), so the oracle is accurate to ~1e-31. Each gap this
/// file asserts is an f64 ulp (~1e-17), fourteen decades larger.
const ORACLE_SLOP: f64 = 1e-25;

/// The EXACT largest eigenvalue real part of `A(μ)` at ~106-bit precision:
/// `−1 + √(max(μ−1, 0))`. For `μ < 1` the eigenvalues are a complex conjugate
/// pair whose real part is exactly `−1`.
fn oracle_largest_real_part(mu: f64) -> Dd {
    let radicand = Dd::from_f64(mu) - Dd::ONE;
    if radicand.hi < 0.0 {
        return -Dd::ONE;
    }
    radicand.sqrt() - Dd::ONE
}

/// `oracle − x` as an f64 (the residual is far above the oracle's own error).
fn oracle_minus(oracle: Dd, x: f64) -> f64 {
    (oracle - Dd::from_f64(x)).to_f64()
}

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
    // the INDEPENDENT eigenvalue criterion (necessary+sufficient) reaches the
    // same boundary — the P=I certificate is tight; and the fs-sos / fs-spectral
    // implementations of the sufficient condition agree at every sample.
    assert!(
        report.boundaries_agree,
        "lyapunov {} vs eigen {}",
        report.lyapunov_boundary, report.eigen_boundary
    );
    assert!(report.stability_classifications_agree);
    let bracket = report
        .boundary_bracket
        .expect("the sweep must witness a shared stability transition");
    assert!(bracket[0] < 2.0 && bracket[1] >= 2.0, "{bracket:?}");
    assert!(report.impl_consistent, "fs-sos and fs-spectral disagree");
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
    assert!(matches!(
        report.witness_decay_rate_color,
        Some(Color::Verified { .. })
    ));
    println!(
        "{{\"campaign\":\"fluttercert\",\"lyapunov_boundary\":{:.3},\"eigen_boundary\":{:.3},\
         \"naive_boundary\":{:.3},\"aitken_boundary\":{:.3},\"witness_mu\":{:.3}}}",
        report.lyapunov_boundary,
        report.eigen_boundary,
        report.naive_boundary,
        report.aitken_boundary,
        mu,
    );
}

#[test]
fn beyond_the_boundary_is_not_certified() {
    let report = run_campaign(2.1, 3.0, 10);
    // every sample is past μ* = 2 → no Lyapunov certificate anywhere, and BOTH
    // abscissae (symmetric-part and actual-eigenvalue) are non-negative.
    assert!(report.samples.iter().all(|s| !s.lyapunov_stable));
    assert!(report.samples.iter().all(|s| s.numerical_abscissa >= 0.0));
    assert!(report.samples.iter().all(|s| s.spectral_abscissa >= 0.0));
    assert!(!report.boundaries_agree);
    assert_eq!(report.boundary_bracket, None);
}

#[test]
fn co_truncated_sweeps_do_not_mint_boundary_agreement() {
    let below = run_campaign(0.2, 1.8, 9);
    assert_eq!(
        below.lyapunov_boundary.to_bits(),
        below.eigen_boundary.to_bits(),
        "both criteria legitimately share the same largest stable sample"
    );
    assert!(below.stability_classifications_agree);
    assert!(
        !below.boundaries_agree,
        "equal co-truncated maxima do not locate a transition"
    );
    assert_eq!(below.boundary_bracket, None);

    let descending = run_campaign(2.45, 0.55, 20);
    assert!(descending.stability_classifications_agree);
    assert!(!descending.boundaries_agree);
    assert_eq!(descending.boundary_bracket, None);
}

/// Regression for bead `frankensim-extreal-program-f85xj.2.34`.
///
/// The witness color used to be `Verified{lo: spectral_abscissa(μ), hi: 0.0}`
/// over a quantity the public surface never named. Two things were wrong:
/// the round-to-nearest `−1 + √(μ−1)` can land ABOVE the exact abscissa, so the
/// "certified lower bound" excluded the true value, and reading `[lo, 0]` as an
/// enclosure of *the eigenvalues' real parts* is falsified for `μ > 1` by the
/// second eigenvalue `−1 − √(μ−1)`.
///
/// `μ = 1.3` is the first sample of this sweep and its witness. It is one of the
/// reachable `μ` where the nearest-rounded evaluation overshoots — the default
/// sweep's witness (`μ = 1.05`) happens to round the other way, which is what
/// kept the defect invisible.
#[test]
fn the_witness_decay_rate_enclosure_is_outward_rounded_and_names_one_eigenvalue() {
    let report = run_campaign(1.3, 1.9, 7);
    let mu = report.witness_mu.expect("this sweep has a witness");
    assert_eq!(mu.to_bits(), 1.3_f64.to_bits(), "witness μ {mu}");
    let Some(Color::Verified { lo, hi }) = report.witness_decay_rate_color else {
        panic!(
            "expected a Verified decay-rate enclosure, got {:?}",
            report.witness_decay_rate_color
        );
    };

    let oracle = oracle_largest_real_part(mu);
    // (1) The published interval really encloses the exact largest real part.
    assert!(
        oracle_minus(oracle, lo) >= -ORACLE_SLOP,
        "lo {lo} is above the exact abscissa by {}",
        -oracle_minus(oracle, lo)
    );
    assert!(
        oracle_minus(oracle, hi) <= ORACLE_SLOP,
        "hi {hi} is below the exact abscissa by {}",
        oracle_minus(oracle, hi)
    );

    // (2) The OLD endpoint would have failed exactly that check here: the
    // nearest-rounded evaluation sits strictly ABOVE the exact value, so a
    // "certified lower bound" published from it excluded the truth.
    let nearest = spectral_abscissa(mu);
    assert!(
        oracle_minus(oracle, nearest) < -ORACLE_SLOP,
        "fixture is stale: the nearest-rounded abscissa {nearest} no longer overshoots"
    );
    assert!(lo < nearest, "outward rounding must move lo down");

    // (3) The claim names the LARGEST real part only. The operator's second
    // eigenvalue is far below the enclosure and is deliberately not covered —
    // this is why the field is not called `witness_color`.
    let second_real_part = -1.0 - (mu - 1.0).sqrt();
    assert!(
        second_real_part < lo,
        "second eigenvalue {second_real_part} vs enclosure [{lo}, {hi}]"
    );
    assert!(hi < 0.0, "the witness must still be a decaying sample");
}

/// The enclosure is sound on both branches of the abscissa and refuses to
/// fabricate a bound for a non-finite μ.
#[test]
fn the_decay_rate_enclosure_is_sound_across_the_branch_and_fails_open_to_no_claim() {
    // μ < 1: complex pair, real part exactly −1 — the branch whose exact
    // `sqrt(0)` used to make the defect invisible.
    for mu in [0.05, 0.55, 0.95, 1.0] {
        let iv = spectral_abscissa_interval(mu);
        let oracle = oracle_largest_real_part(mu);
        assert!(
            oracle_minus(oracle, iv.lo()) >= -ORACLE_SLOP
                && oracle_minus(oracle, iv.hi()) <= ORACLE_SLOP,
            "μ={mu}: [{}, {}] does not enclose {oracle:?}",
            iv.lo(),
            iv.hi()
        );
    }
    // μ ≥ 1: the real branch, over the whole reachable sweep range.
    for k in 0..=200 {
        let mu = 1.0 + 2.0 * f64::from(k) / 200.0;
        let iv = spectral_abscissa_interval(mu);
        let oracle = oracle_largest_real_part(mu);
        assert!(
            oracle_minus(oracle, iv.lo()) >= -ORACLE_SLOP
                && oracle_minus(oracle, iv.hi()) <= ORACLE_SLOP,
            "μ={mu}: [{}, {}] does not enclose {oracle:?}",
            iv.lo(),
            iv.hi()
        );
        assert!(iv.lo() <= iv.hi());
    }
    for mu in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let iv = spectral_abscissa_interval(mu);
        assert_eq!(
            iv.lo().to_bits(),
            f64::NEG_INFINITY.to_bits(),
            "μ={mu} must mint no bound"
        );
        assert_eq!(
            iv.hi().to_bits(),
            f64::INFINITY.to_bits(),
            "μ={mu} must mint no bound"
        );
    }
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
