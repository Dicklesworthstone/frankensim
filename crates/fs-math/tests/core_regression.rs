//! Core-only and worst-case-point regression tests for fs-math's strict
//! trig. Born as the debugging probe that isolated two real bugs during
//! implementation: (1) decimal-rendered Cody–Waite constants losing their
//! exact trailing-zero bit patterns, and (2) Taylor cores two terms short
//! (166 ULP at the |r| = π/4 interval edge). These assertions keep both
//! mistakes dead.

#[test]
fn sin_core_only_stays_within_one_ulp() {
    // |x| ≤ 0.78 needs no argument reduction: isolates the polynomial core.
    let mut worst = 0u64;
    let mut s = 1u64;
    for _ in 0..200_000 {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let x = (((s >> 11) as f64) / (1u64 << 53) as f64 - 0.5) * 1.56;
        worst = worst.max(fs_math::ulp_distance(fs_math::det::sin(x), x.sin()));
    }
    assert!(
        worst <= 1,
        "core-only sin max ULP {worst} (was 166 with the short core)"
    );
}

#[test]
fn known_worst_case_points_stay_tight() {
    // The two reduction-stress points found by the ULP sweep, plus the
    // near-zero quadrant boundaries that demanded the 4th pi/2 chunk.
    for (x, budget) in [
        (172_940.748_752_312_73_f64, 3u64),
        (953_301.368_295_506_6, 3),
        (std::f64::consts::PI, 3),
        (3.0 * std::f64::consts::FRAC_PI_2, 3),
    ] {
        let ds = fs_math::ulp_distance(fs_math::det::sin(x), x.sin());
        let dc = fs_math::ulp_distance(fs_math::det::cos(x), x.cos());
        assert!(ds <= budget, "sin({x}): {ds} ULP");
        assert!(dc <= budget, "cos({x}): {dc} ULP");
    }
}

#[test]
fn pio2_split_parts_reconstruct_pi_over_two() {
    // The 3+1-part split must sum to pi/2 well beyond f64 precision; checked
    // here in f64 to the extent representable (guards typo'd bit patterns).
    let p1 = f64::from_bits(0x3FF9_21FB_5440_0000);
    let p2 = f64::from_bits(0x3DD0_B461_1A60_0000);
    let p3 = f64::from_bits(0x3BA3_198A_2E00_0000);
    // Sterbenz-exact tail of the ROUNDED pi/2:
    let tail = std::f64::consts::FRAC_PI_2 - p1;
    // p2 + p3 must agree with that tail to within pi/2's own rounding error.
    assert!(
        ((tail - p2) - p3).abs() < 1.2e-16,
        "split drifted: {}",
        ((tail - p2) - p3).abs()
    );
    // Trailing-zero exactness property: products with |j| <= 2^20 are exact.
    for p in [p1, p2, p3] {
        assert_eq!(
            p.to_bits() & 0xF_FFFF,
            0,
            "part {p:e} lost its trailing zeros"
        );
    }
}
