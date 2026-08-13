//! Exact-predicate conformance suite (plan §13.3; the 6ys.14 bead): any
//! reimplementation must pass. The deep batteries live as unit tests beside
//! the implementation; this suite exercises the PUBLIC surface and emits
//! the JSONL verdicts the bead's acceptance criteria name.

use fs_ivl::{
    Sign, Stage, incircle, insphere, orient2d, orient2d_sos, orient2d_with_stage, orient3d,
};
use std::panic::{AssertUnwindSafe, catch_unwind};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ivl/predicates\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

#[test]
fn pd_001_adversarial_degeneracies_are_exact_zero() {
    // Cocircular lattice points (x² + y² = 25) in every rotation.
    let ring = [
        [5.0, 0.0],
        [3.0, 4.0],
        [-3.0, 4.0],
        [0.0, -5.0],
        [4.0, -3.0],
    ];
    let mut zeros = 0;
    for i in 0..ring.len() {
        for j in 0..ring.len() {
            for k in 0..ring.len() {
                for l in 0..ring.len() {
                    if [i, j, k, l].windows(2).all(|w| w[0] != w[1]) && i != k && i != l && j != l {
                        assert_eq!(
                            incircle(ring[i], ring[j], ring[k], ring[l]),
                            Sign::Zero,
                            "cocircular points must be exactly degenerate ({i},{j},{k},{l})"
                        );
                        zeros += 1;
                    }
                }
            }
        }
    }
    // Cospherical lattice points (x² + y² + z² = 9).
    let s = [
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 3.0],
        [-3.0, 0.0, 0.0],
        [0.0, 0.0, -3.0],
    ];
    assert_eq!(insphere(s[0], s[1], s[2], s[3], s[4]), Sign::Zero);
    assert_eq!(
        orient3d(
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],
            [3.0, 6.0, 9.0],
            [4.0, 8.0, 12.0]
        ),
        Sign::Zero,
        "collinear 3D grid must be coplanar-degenerate"
    );
    verdict(
        "pd-001",
        &format!("{zeros} cocircular orderings + cospherical + collinear all Zero"),
    );
}

#[test]
fn pd_002_one_ulp_class_perturbations_resolve_correctly() {
    // Exact dyadic radial scaling: truth known analytically (outside).
    let grow = 1.0 + fs_math::det::powi(2.0, -40);
    let (a, b, c) = ([5.0, 0.0], [3.0, 4.0], [-3.0, 4.0]);
    assert_eq!(incircle(a, b, c, [0.0, -5.0]), Sign::Zero);
    assert_eq!(incircle(a, b, c, [0.0, -5.0 * grow]), Sign::Negative);
    // The same points 1 ulp INSIDE via next-toward-zero.
    let inward = f64::from_bits((5.0f64).to_bits() - 1);
    assert_eq!(incircle(a, b, c, [0.0, -inward]), Sign::Positive);
    verdict(
        "pd-002",
        "dyadic-scaled and 1-ulp radial perturbations classified correctly",
    );
}

#[test]
fn pd_003_sos_ties_break_deterministically() {
    // Coincident + collinear worst case, all index assignments.
    let p = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]];
    let mut seen = Vec::new();
    for perm in [[0u64, 1, 2], [1, 0, 2], [2, 0, 1]] {
        let s = orient2d_sos(p[0], p[1], p[2], perm[0], perm[1], perm[2]);
        assert_ne!(s, Sign::Zero, "SoS is total even on coincident points");
        seen.push((perm, s));
        // Re-ask: identical.
        assert_eq!(orient2d_sos(p[0], p[1], p[2], perm[0], perm[1], perm[2]), s);
    }
    verdict(
        "pd-003",
        &format!("coincident-point ties resolved deterministically: {seen:?}"),
    );
}

#[test]
fn pd_004_filter_rate_measured_and_logged() {
    let mut seed = 0x5EED_F117_0000_0BBBu64;
    let mut rnd = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64) / (1u64 << 53) as f64 * 1000.0 - 500.0
    };
    let (mut filtered, mut total) = (0usize, 0usize);
    for _ in 0..5000 {
        let (a, b, c) = ([rnd(), rnd()], [rnd(), rnd()], [rnd(), rnd()]);
        total += 1;
        if orient2d_with_stage(a, b, c).1 == Stage::Filtered {
            filtered += 1;
        }
        // Exactness spot-check: antisymmetry on the same sample.
        assert_eq!(orient2d(a, b, c), orient2d(b, a, c).flip());
    }
    let rate = filtered as f64 / total as f64;
    assert!(rate > 0.99, "stage-A filter rate collapsed: {rate}");
    verdict(
        "pd-004",
        &format!("orient2d stage-A filter rate {rate:.4} over {total} general-position samples"),
    );
}

#[test]
fn pd_005_invalid_numeric_domains_fail_closed() {
    let fails_closed = |case: &str, f: &mut dyn FnMut()| {
        assert!(
            catch_unwind(AssertUnwindSafe(f)).is_err(),
            "{case} returned a predicate sign outside the certified numeric domain"
        );
    };

    fails_closed("orient2d non-finite", &mut || {
        let _ = orient2d([f64::NAN, 0.0], [1.0, 0.0], [0.0, 1.0]);
    });
    fails_closed("orient3d non-finite", &mut || {
        let _ = orient3d(
            [f64::INFINITY, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0; 3],
        );
    });
    fails_closed("incircle non-finite", &mut || {
        let _ = incircle([0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [f64::NAN, 0.0]);
    });
    fails_closed("insphere non-finite", &mut || {
        let _ = insphere(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, f64::INFINITY, 0.0],
        );
    });

    let m = f64::MAX;
    fails_closed("orient2d finite overflow", &mut || {
        let _ = orient2d([m, 0.0], [0.0, m], [-m, 0.0]);
    });
    fails_closed("orient3d finite overflow", &mut || {
        let _ = orient3d([m, 0.0, 0.0], [0.0, m, 0.0], [0.0, 0.0, m], [-m, 0.0, 0.0]);
    });
    fails_closed("incircle finite overflow", &mut || {
        let _ = incircle([m, 0.0], [0.0, m], [0.0, -m], [-m, 0.0]);
    });
    fails_closed("insphere finite overflow", &mut || {
        let _ = insphere(
            [m, 0.0, 0.0],
            [0.0, m, 0.0],
            [0.0, 0.0, m],
            [0.0, -m, 0.0],
            [-m, 0.0, 0.0],
        );
    });

    verdict(
        "pd-005",
        "all four public predicates rejected non-finite coordinates and finite-coordinate overflow",
    );
}

/// Bead frankensim-d8yax (mutation-campaign survivor): the stage-A filter
/// threshold was never load-bearing on any committed input, because every
/// prior fixture's naive determinant is EXACT (integer lattices, small
/// dyadics) — a 64× overconfident CCWERRBOUND_A survived the whole battery.
/// This fixture family makes the threshold load-bearing: two points exactly
/// on `y = x` at ~1e7 with a third displaced by k tiny steps produce a
/// catastrophically cancelled naive determinant whose float value is
/// rounding garbage BELOW the true stage-A bound. The battery asserts, for
/// every case that lands inside a would-be-weakened acceptance band, that
/// (a) the ladder does NOT resolve at the filtered stage, and (b) the
/// returned sign equals the sign known by construction. A filter accepting
/// inside its proven bound fails (a) immediately and (b) whenever the
/// cancellation garbage disagrees — this is the test that kills the
/// `weaken-orient-filter` mutant in scripts/ci/e2e_extreal_arithmetic_audit.sh.
#[test]
fn pd_006_stage_a_filter_threshold_is_load_bearing_on_inexact_determinants() {
    // Shewchuk's bound, recomputed independently of the crate's private
    // constant (half-ulp epsilon).
    let eps = f64::EPSILON / 2.0;
    let ccwerrbound_a = (3.0 + 16.0 * eps) * eps;
    // The mutant under audit scales the bound by 2^-6; the band between the
    // two is where an overconfident filter resolves and the real one must
    // not. Checking a FAMILY of weakenings (2^-1 .. 2^-6) kills any
    // wrong-side scaling, not just the campaign's exact operator.
    let weakest = ccwerrbound_a * 0.015_625;

    let base1 = 10_000_001.0f64;
    let base2 = 10_000_002.0f64; // b2 − b1 = 1 exactly
    let c = 0.5f64;
    // Displacement quantum = ULP at the ~1e7 magnitude (2^-29, exact bit
    // construction): the ladder's pc-relative subtractions then preserve δ
    // exactly, while the ~1e14 products still round — cancellation garbage
    // with a known true sign. A smaller δ silently rounds away in the
    // subtraction and every determinant is exactly zero (this test's own
    // vacuity guard caught that on its first run).
    let step = f64::from_bits((1023u64 - 29) << 52);

    let mut in_band = 0u64;
    let mut checked = 0u64;
    for k in 1..=512i64 {
        for sign in [1.0f64, -1.0] {
            #[allow(clippy::cast_precision_loss)]
            let delta = sign * (k as f64) * step;
            let pa = [base1, base1];
            let pb = [base2, base2];
            let pc = [c, c + delta];
            // Recompute the naive stage-A quantities exactly as the ladder
            // arranges them, to classify this case's band membership.
            let acx = pa[0] - pc[0];
            let bcx = pb[0] - pc[0];
            let acy = pa[1] - pc[1];
            let bcy = pb[1] - pc[1];
            let det = acx * bcy - acy * bcx;
            let detsum = (acx * bcy).abs() + (acy * bcx).abs();
            if det == 0.0 {
                continue;
            }
            let inside_true_bound = det.abs() <= ccwerrbound_a * detsum;
            let above_weakest = det.abs() > weakest * detsum;
            if !(inside_true_bound && above_weakest) {
                continue;
            }
            in_band += 1;
            // Exact truth by construction: det_exact = δ·(b2−b1) = δ.
            let expected = if delta > 0.0 {
                Sign::Positive
            } else {
                Sign::Negative
            };
            let (got, stage) = orient2d_with_stage(pa, pb, pc);
            assert_eq!(
                got, expected,
                "k={k} sign={sign}: cancellation garbage must not decide the sign"
            );
            assert_ne!(
                stage,
                Stage::Filtered,
                "k={k} sign={sign}: |det| is inside the proven stage-A error \
                 bound; a filter resolving here is running a weakened bound"
            );
            checked += 1;
        }
    }
    assert!(
        in_band >= 64,
        "the fixture family must actually populate the weakened acceptance \
         band, or this battery is vacuous: {in_band}"
    );
    assert_eq!(in_band, checked);
    verdict(
        "pd-006",
        &format!("{in_band} inexact-determinant cases pinned the stage-A threshold"),
    );
}
