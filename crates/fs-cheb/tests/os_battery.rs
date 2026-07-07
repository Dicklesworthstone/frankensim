//! Orr–Sommerfeld battery (urvw item 2): stability verdicts on both
//! sides of criticality and THE acceptance test — the neutral crossing
//! at α = 1.02056 reproducing the published Re_c ≈ 5772.22 — plus the
//! cross-ISA golden hash over the growth-rate query.

use fs_cheb::orr_sommerfeld::{critical_reynolds, growth_rates, max_growth};

const N: usize = 48; // collocation order (47 interior nodes)
const ALPHA_C: f64 = 1.02056;

#[test]
fn stability_verdicts_bracket_criticality() {
    // Well below critical: stable (all growth rates negative).
    let stable = max_growth(2000.0, ALPHA_C, N).expect("solve");
    assert!(stable < 0.0, "Re=2000 must be stable: {stable}");
    // Well above critical: the classic instability appears.
    let unstable = max_growth(10_000.0, ALPHA_C, N).expect("solve");
    assert!(unstable > 0.0, "Re=10000 must be unstable: {unstable}");
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"os-verdicts\",\"verdict\":\"pass\",\"detail\":\"growth(2000)={stable:.3e} growth(10000)={unstable:.3e}\"}}"
    );
}

#[test]
fn critical_reynolds_matches_published_value() {
    // THE acceptance test (plan §15.3 / bead 6ys.15 criterion): published
    // Re_c = 5772.22 at α = 1.02056 (Orszag 1971).
    let rc = critical_reynolds(ALPHA_C, N, 4000.0, 8000.0).expect("bisection");
    assert!(
        (rc - 5772.22).abs() < 5.0,
        "critical Reynolds {rc} vs published 5772.22"
    );
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"os-critical\",\"verdict\":\"pass\",\"detail\":\"Re_c = {rc:.2} vs published 5772.22 (N={N})\"}}"
    );
}

#[test]
fn growth_rate_query_shape() {
    // σ₁..σ₈ as a first-class query: descending real parts, finite.
    let sigmas = growth_rates(5000.0, ALPHA_C, N, 8).expect("solve");
    assert_eq!(sigmas.len(), 8);
    for w in sigmas.windows(2) {
        assert!(w[0].re >= w[1].re, "growth rates must descend: {sigmas:?}");
    }
    assert!(sigmas.iter().all(|z| z.re.is_finite() && z.im.is_finite()));
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x0; // placeholder: set from first run

#[test]
fn os_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    for &(re, al) in &[(2000.0f64, 1.0f64), (5772.0, ALPHA_C), (9000.0, 0.9)] {
        for z in growth_rates(re, al, 40, 4).expect("solve") {
            feed(z.re);
            feed(z.im);
        }
    }
    println!(
        "{{\"suite\":\"fs-cheb\",\"case\":\"os-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "OS spectrum bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
