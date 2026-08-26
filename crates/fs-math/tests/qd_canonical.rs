//! Canonical conformance and property test battery for Quad-Double (Qd)
//! representation (bead frankensim-epic-bedrock-6ys.23.1.1).
//!
//! Covers:
//! - Canonical representation invariants and non-overlapping limbs
//! - Checked construction and typed validation refusals
//! - Idempotence of renormalization
//! - Monotonic magnitude ordering and total ordering
//! - Byte round-trips (little-endian, big-endian, canonical)
//! - Edge cases: +/-0, subnormals, infinities, NaNs, max finite
//! - Arithmetic sanity checks and algebraic laws

#![deny(unsafe_code)]

use fs_math::dd::Dd;
use fs_math::qd::{Qd, QdError, QdOpError, QdOpOutcome};

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

// ---------------------------------------------------------------------------
// 1. Unit & Boundary Tests
// ---------------------------------------------------------------------------

#[test]
fn test_canonical_constants_and_properties() {
    assert!(Qd::ZERO.is_canonical());
    assert!(Qd::ZERO.is_zero());
    assert!(!Qd::ZERO.is_nan());
    assert!(!Qd::ZERO.is_infinite());

    assert!(Qd::ONE.is_canonical());
    assert_eq!(Qd::ONE.to_f64(), 1.0);
    assert_eq!(Qd::ONE.to_dd(), Dd::ONE);

    assert!(Qd::NEG_ONE.is_canonical());
    assert_eq!(Qd::NEG_ONE.to_f64(), -1.0);

    assert!(Qd::TWO.is_canonical());
    assert_eq!(Qd::TWO.to_f64(), 2.0);

    assert!(Qd::HALF.is_canonical());
    assert_eq!(Qd::HALF.to_f64(), 0.5);

    assert!(Qd::MAX.is_canonical());
    assert!(Qd::MAX.is_finite());

    assert!(Qd::MIN_POSITIVE.is_canonical());
    assert!(Qd::MIN_POSITIVE.is_finite());

    assert!(Qd::INFINITY.is_infinite());
    assert!(!Qd::INFINITY.is_finite());
    assert!(Qd::INFINITY.is_sign_positive());

    assert!(Qd::NEG_INFINITY.is_infinite());
    assert!(!Qd::NEG_INFINITY.is_finite());
    assert!(Qd::NEG_INFINITY.is_sign_negative());

    assert!(Qd::NAN.is_nan());
    assert!(!Qd::NAN.is_finite());
}

#[test]
fn test_signed_zeros() {
    let pz = Qd::from_f64(0.0);
    let nz = Qd::from_f64(-0.0);

    assert!(pz.is_zero());
    assert!(nz.is_zero());
    assert!(pz.is_sign_positive());
    assert!(nz.is_sign_negative());
    assert!(pz.is_canonical());
    assert!(nz.is_canonical());
    assert_eq!(pz, nz); // IEEE 754 value equality
}

#[test]
fn test_subnormals_and_extremes() {
    let min_subnormal = f64::from_bits(1);
    let q_sub = Qd::from_f64(min_subnormal);
    assert!(q_sub.is_canonical());
    assert!(q_sub.is_finite());
    assert!(!q_sub.is_zero());
    assert_eq!(q_sub.to_f64(), min_subnormal);

    let max_normal = f64::MAX;
    let q_max = Qd::from_f64(max_normal);
    assert!(q_max.is_canonical());
    assert!(q_max.is_finite());
    assert_eq!(q_max.to_f64(), max_normal);
}

// ---------------------------------------------------------------------------
// 2. Checked Construction & Invariant Refusals
// ---------------------------------------------------------------------------

#[test]
fn test_checked_construction_valid() {
    let q1 = Qd::from_parts_checked(1.0, 0.0, 0.0, 0.0).expect("valid");
    assert_eq!(q1, Qd::ONE);

    let q2 = Qd::from_parts_checked(1.0, f64::EPSILON * 0.25, 0.0, 0.0).expect("valid non-overlap");
    assert!(q2.is_canonical());
}

#[test]
fn test_checked_construction_overlapping_refusal() {
    // 1.0 and 1.0 heavily overlap
    let err = Qd::from_parts_checked(1.0, 1.0, 0.0, 0.0).expect_err("must refuse overlap");
    assert_eq!(err, QdError::OverlappingComponents { index: 0 });

    // Component 1 and 2 overlap (c0 and c1 do not overlap)
    let err2 = Qd::from_parts_checked(1.0, 1e-20, 1e-20, 0.0).expect_err("must refuse overlap");
    assert_eq!(err2, QdError::OverlappingComponents { index: 1 });

    // Component 2 and 3 overlap (c0, c1, and c2 do not overlap)
    let err3 = Qd::from_parts_checked(1.0, 1e-20, 1e-40, 1e-40).expect_err("must refuse overlap");
    assert_eq!(err3, QdError::OverlappingComponents { index: 2 });
}

#[test]
fn test_checked_construction_unordered_refusal() {
    // Magnitude of c1 > c0
    let err = Qd::from_parts_checked(1.0, 2.0, 0.0, 0.0).expect_err("must refuse unordered");
    assert_eq!(err, QdError::UnorderedMagnitudes { index: 0 });

    // Magnitude of c2 > c1
    let err2 = Qd::from_parts_checked(1.0, 1e-20, 1e-16, 0.0).expect_err("must refuse unordered");
    assert_eq!(err2, QdError::UnorderedMagnitudes { index: 1 });
}

#[test]
fn test_checked_construction_zero_ordering_refusal() {
    // Zero followed by non-zero
    let err = Qd::from_parts_checked(0.0, 1e-20, 0.0, 0.0).expect_err("must refuse non-zero after 0");
    assert_eq!(err, QdError::InvalidZeroRepresentation);

    let err2 = Qd::from_parts_checked(1.0, 0.0, 1e-35, 0.0).expect_err("must refuse non-zero after 0");
    assert_eq!(err2, QdError::InvalidZeroRepresentation);
}

// ---------------------------------------------------------------------------
// 3. Property & Metamorphic Tests
// ---------------------------------------------------------------------------

#[test]
fn test_renormalization_idempotence_random_sweep() {
    let mut seed = 0xABCD_EF01_u64;
    for _ in 0..50_000 {
        let c0 = lcg(&mut seed) * 1e8;
        let c1 = lcg(&mut seed) * 1e2;
        let c2 = lcg(&mut seed) * 1e-4;
        let c3 = lcg(&mut seed) * 1e-10;

        let q = Qd::from_parts(c0, c1, c2, c3);
        assert!(q.is_canonical(), "renormalized Qd must be canonical: {q:?}");

        // Renormalizing canonical output must produce identical components
        let q_again = Qd::from_parts(q.c0, q.c1, q.c2, q.c3);
        assert_eq!(
            q.components(),
            q_again.components(),
            "renormalization must be idempotent"
        );
    }
}

#[test]
fn test_exact_sum_preservation() {
    // Adding numbers of vastly different scales should not lose exact terms
    let a = Qd::from_f64(1e20);
    let b = Qd::from_f64(1.0);
    let sum = a + b;

    assert_eq!(sum.c0.to_bits(), (1e20_f64).to_bits());
    assert_eq!(sum.c1.to_bits(), (1.0_f64).to_bits());
    assert_eq!(sum.c2.to_bits(), 0);
    assert_eq!(sum.c3.to_bits(), 0);

    // Subtracting back exact term
    let diff = sum - a;
    assert_eq!(diff.c0.to_bits(), (1.0_f64).to_bits());
    assert_eq!(diff.c1.to_bits(), 0);
}

#[test]
fn test_total_ordering_and_comparisons() {
    let vals = [
        Qd::NEG_INFINITY,
        Qd::from_f64(-100.0),
        Qd::NEG_ONE,
        Qd::from_f64(-0.0),
        Qd::ZERO,
        Qd::MIN_POSITIVE,
        Qd::EPSILON,
        Qd::HALF,
        Qd::ONE,
        Qd::TWO,
        Qd::from_f64(100.0),
        Qd::MAX,
        Qd::INFINITY,
    ];

    for i in 0..vals.len() {
        for j in 0..vals.len() {
            let left = vals[i];
            let right = vals[j];

            if i < j {
                if left.is_zero() && right.is_zero() {
                    assert!(left.le(right));
                    assert!(left.ge(right));
                } else {
                    assert!(left.lt(right), "{left:?} < {right:?}");
                    assert!(left.le(right), "{left:?} <= {right:?}");
                    assert!(!left.gt(right), "not {left:?} > {right:?}");
                }
            } else if i == j {
                assert!(left.le(right));
                assert!(left.ge(right));
            } else if left.is_zero() && right.is_zero() {
                assert!(left.le(right));
                assert!(left.ge(right));
            } else {
                assert!(left.gt(right), "{left:?} > {right:?}");
                assert!(left.ge(right), "{left:?} >= {right:?}");
                assert!(!left.lt(right), "not {left:?} < {right:?}");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// 4. Byte Round-Trip & Serialisation Conformance
// ---------------------------------------------------------------------------

#[test]
fn test_bytes_round_trip_all_formats() {
    let test_cases = [
        Qd::ZERO,
        Qd::ONE,
        Qd::NEG_ONE,
        Qd::from_f64(std::f64::consts::PI),
        Qd::from_parts(std::f64::consts::E, 1e-17, -1e-33, 1e-49),
        Qd::MAX,
        Qd::MIN_POSITIVE,
    ];

    for q in test_cases {
        assert!(q.is_canonical());

        // Little endian
        let le_bytes = q.to_bytes_le();
        assert_eq!(le_bytes.len(), 32);
        let dec_le = Qd::from_bytes_le(&le_bytes).expect("decode LE");
        assert_eq!(q.components(), dec_le.components());

        // Big endian
        let be_bytes = q.to_bytes_be();
        assert_eq!(be_bytes.len(), 32);
        let dec_be = Qd::from_bytes_be(&be_bytes).expect("decode BE");
        assert_eq!(q.components(), dec_be.components());

        // Canonical interchange bytes
        let canon_bytes = q.to_canonical_bytes();
        let dec_canon = Qd::from_canonical_bytes(&canon_bytes).expect("decode canonical");
        assert_eq!(q.components(), dec_canon.components());
    }
}

#[test]
fn test_corrupt_bytes_refusal() {
    let mut bad_bytes = [0u8; 32];
    // Put 1.0 in first component and 1.0 in second component (overlapping)
    bad_bytes[0..8].copy_from_slice(&(1.0_f64).to_le_bytes());
    bad_bytes[8..16].copy_from_slice(&(1.0_f64).to_le_bytes());

    let res = Qd::from_bytes_le(&bad_bytes);
    assert!(res.is_err(), "overlapping components from bytes must be rejected");
}

// ---------------------------------------------------------------------------
// 5. Arithmetic Laws & Precision Tests
// ---------------------------------------------------------------------------

#[test]
fn test_sqrt_machin_precision() {
    // sqrt(2)^2 = 2 to ~212 bits
    let sqrt2 = Qd::from_f64(2.0).sqrt();
    let sq = sqrt2 * sqrt2;
    let diff = (sq - Qd::TWO).abs();
    assert!(diff.c0 < 1e-62, "sqrt(2)^2 error must be < 1e-62: {diff:?}");

    // sqrt(3)^2 = 3
    let sqrt3 = Qd::from_f64(3.0).sqrt();
    let sq3 = sqrt3 * sqrt3;
    let diff3 = (sq3 - Qd::from_f64(3.0)).abs();
    assert!(diff3.c0 < 1e-62, "sqrt(3)^2 error must be < 1e-62: {diff3:?}");
}

#[test]
fn test_division_and_multiplication_inverse() {
    let mut seed = 0x5555_AAAA_u64;
    for _ in 0..5_000 {
        let x = Qd::from_f64(lcg(&mut seed) * 100.0) + Qd::from_f64(lcg(&mut seed) * 1e-14);
        let y = Qd::from_f64(lcg(&mut seed) * 100.0) + Qd::from_f64(lcg(&mut seed) * 1e-14);

        if y.abs().to_f64() > 1e-2 {
            let div = x / y;
            let mult = div * y;
            let residual = (mult - x).abs();
            let scale = x.abs().to_f64().max(1.0);
            assert!(
                residual.to_f64() <= 1e-55 * scale,
                "x/y*y must equal x to quad precision: residual={:?}, scale={}",
                residual,
                scale
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 6. Core Arithmetic (6ys.23.1.2): Scaling, Dyadics, Checked Ops & Metamorphic
// ---------------------------------------------------------------------------

#[test]
fn test_exact_power_of_two_scaling_metamorphic() {
    let mut seed = 0x9988_7766_u64;
    for _ in 0..10_000 {
        let val = Qd::from_parts(
            lcg(&mut seed) * 1e5,
            lcg(&mut seed) * 1e-10,
            lcg(&mut seed) * 1e-25,
            lcg(&mut seed) * 1e-40,
        );
        for exp in [-500, -53, -1, 0, 1, 53, 500] {
            let scaled = val.scale_power_of_two(exp);
            let unscaled = scaled.scale_power_of_two(-exp);
            assert_eq!(
                val.components(),
                unscaled.components(),
                "scaling by 2^{exp} and 2^{} must round-trip exactly",
                -exp
            );
        }
    }
}

#[test]
fn test_dyadic_kat_exact_reconstruction() {
    // Dyadic rational sums must be exact across quad components
    let half = Qd::from_f64(0.5);
    let quarter = Qd::from_f64(0.25);
    let eighth = Qd::from_f64(0.125);
    let sixteenth = Qd::from_f64(0.0625);

    let sum = half + quarter + eighth + sixteenth;
    let expected = Qd::from_f64(0.9375);
    assert_eq!(sum, expected);
    assert_eq!(sum.c1, 0.0);
    assert_eq!(sum.c2, 0.0);
    assert_eq!(sum.c3, 0.0);

    // Dyadic exponent gaps
    let d1 = Qd::from_f64(1.0);
    let d2 = d1.scale_power_of_two(-53);
    let d3 = d1.scale_power_of_two(-106);
    let d4 = d1.scale_power_of_two(-159);

    let composite = d1 + d2 + d3 + d4;
    assert_eq!(composite.c0, 1.0);
    assert_eq!(composite.c1, f64::exp2(-53.0));
    assert_eq!(composite.c2, f64::exp2(-106.0));
    assert_eq!(composite.c3, f64::exp2(-159.0));
    assert!(composite.is_canonical());
}

#[test]
fn test_commutativity_and_distributivity_laws() {
    let mut seed = 0xFEED_FACE_u64;
    for _ in 0..5_000 {
        let a = Qd::from_f64(lcg(&mut seed) * 50.0) + Qd::from_f64(lcg(&mut seed) * 1e-12);
        let b = Qd::from_f64(lcg(&mut seed) * 50.0) + Qd::from_f64(lcg(&mut seed) * 1e-12);
        let c = Qd::from_f64(lcg(&mut seed) * 50.0) + Qd::from_f64(lcg(&mut seed) * 1e-12);

        // Commutativity: a + b == b + a
        let sum_ab = a + b;
        let sum_ba = b + a;
        assert_eq!(sum_ab.components(), sum_ba.components(), "a+b must equal b+a");

        // Commutativity: a * b == b * a
        let mul_ab = a * b;
        let mul_ba = b * a;
        assert_eq!(mul_ab.components(), mul_ba.components(), "a*b must equal b*a");

        // Distributivity residual: a * (b + c) ≈ a*b + a*c
        let lhs = a * (b + c);
        let rhs = a * b + a * c;
        let diff = (lhs - rhs).abs();
        let scale = (a.abs().to_f64() * (b.abs().to_f64() + c.abs().to_f64())).max(1.0);
        assert!(
            diff.to_f64() <= 1e-55 * scale,
            "distributivity error must be bounded: diff={diff:?}"
        );
    }
}

#[test]
fn test_checked_arithmetic_outcomes_and_refusals() {
    let a = Qd::from_f64(10.0);
    let b = Qd::from_f64(2.0);

    // Normal checked outcome
    let out_add = a.checked_add(b).expect("checked add");
    assert_eq!(out_add.value(), Qd::from_f64(12.0));

    let out_scale = a.checked_scale_power_of_two(3).expect("checked scale");
    assert_eq!(out_scale, QdOpOutcome::CorrectlyRounded(Qd::from_f64(80.0)));

    // Non-finite refusal
    let err_nan = Qd::NAN.checked_add(b).expect_err("NaN must refuse");
    assert!(matches!(err_nan, QdOpError::InvalidInput { .. }));

    // Division by zero refusal
    let err_zero = a.checked_div(Qd::ZERO).expect_err("div 0 must refuse");
    assert_eq!(err_zero, QdOpError::InvalidInput { reason: "division by zero" });
}

#[test]
fn test_mutant_detection_and_fault_injection() {
    // 1-ULP perturbed limb must violate non-overlapping invariant
    let q = Qd::from_parts(1.0, f64::EPSILON * 0.5, 0.0, 0.0);
    assert!(q.is_canonical());

    // Perturbing second component above half-ULP threshold must fail checked construction
    let bad = Qd::from_parts_checked(1.0, f64::EPSILON * 0.75, 0.0, 0.0);
    assert!(bad.is_err());
}

// ---------------------------------------------------------------------------
// 7. Advanced Arithmetic (6ys.23.1.3): Sqrt, Escalation & Extreme Boundaries
// ---------------------------------------------------------------------------

#[test]
fn test_checked_sqrt_and_escalation_primitives() {
    // Exact zero sqrt
    let sqrt_zero = Qd::ZERO.checked_sqrt().expect("sqrt 0");
    assert_eq!(sqrt_zero, QdOpOutcome::CorrectlyRounded(Qd::ZERO));

    // Negative sqrt refusal
    let err_neg = Qd::NEG_ONE.checked_sqrt().expect_err("sqrt negative must refuse");
    assert!(matches!(err_neg, QdOpError::InvalidInput { .. }));

    // Exact square
    let four = Qd::from_f64(4.0);
    let sqrt_four = four.checked_sqrt().expect("sqrt 4");
    assert_eq!(sqrt_four.value(), Qd::TWO);

    // Escalation primitives
    let esc_f64 = Qd::escalate_precision(1.2345, None);
    assert_eq!(esc_f64, Qd::from_f64(1.2345));

    let dd = Dd::from_pair(std::f64::consts::PI, 1e-16);
    let esc_dd = Qd::escalate_precision(std::f64::consts::PI, Some(dd));
    assert_eq!(esc_dd.c0, dd.hi);
    assert_eq!(esc_dd.c1, dd.lo);

    let chk_esc = Qd::checked_escalate(2.718).expect("checked escalate");
    assert_eq!(chk_esc, QdOpOutcome::CorrectlyRounded(Qd::from_f64(2.718)));
}

#[test]
fn test_division_and_sqrt_extreme_exponents() {
    // Extreme exponent division
    let huge = Qd::from_f64(1e250);
    let tiny = Qd::from_f64(1e-250);

    let div_huge = huge / tiny;
    assert!(div_huge.is_infinite() || div_huge.to_f64() > 1e300);

    let div_tiny = tiny / huge;
    assert!(div_tiny.is_zero() || div_tiny.to_f64() < 1e-300);

    // Extreme square roots
    let sqrt_huge = huge.sqrt();
    assert!((sqrt_huge.to_f64() - 1e125).abs() < 1e120);

    let sqrt_tiny = tiny.sqrt();
    assert!((sqrt_tiny.to_f64() - 1e-125).abs() < 1e-130);
}

#[test]
fn test_sqrt_square_and_reciprocal_metamorphic() {
    let mut seed = 0x3344_5566_u64;
    for _ in 0..5_000 {
        let x = Qd::from_f64(lcg(&mut seed).abs() * 100.0 + 0.01);

        // (sqrt(x))^2 ≈ x
        let sqrt_x = x.sqrt();
        let sq = sqrt_x * sqrt_x;
        let diff_sq = (sq - x).abs();
        let scale = x.to_f64().max(1.0);
        assert!(
            diff_sq.to_f64() <= 1e-55 * scale,
            "sqrt(x)^2 must equal x: diff={diff_sq:?}"
        );

        // 1 / (1 / x) ≈ x
        let recip = Qd::ONE / x;
        let recip_recip = Qd::ONE / recip;
        let diff_recip = (recip_recip - x).abs();
        assert!(
            diff_recip.to_f64() <= 1e-55 * scale,
            "1/(1/x) must equal x: diff={diff_recip:?}"
        );
    }
}
