//! Exact-predicate conformance suite (plan §13.3; the 6ys.14 bead): any
//! reimplementation must pass. The deep batteries live as unit tests beside
//! the implementation; this suite exercises the PUBLIC surface and emits
//! the JSONL verdicts the bead's acceptance criteria name.

use fs_ivl::{
    Sign, Stage, incircle, insphere, orient2d, orient2d_sos, orient2d_with_stage, orient3d,
};

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
    let grow = 1.0 + (2.0f64).powi(-40);
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
