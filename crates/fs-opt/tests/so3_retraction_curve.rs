//! G0/G3 retraction-curve derivative conformance for the fs-opt SO(3) authority.
//!
//! The fixture stays away from the canonical-representative sign seam. It
//! independently derives the right-composed quaternion curve velocity, checks
//! it against two central-difference scales, and requires antipodal bases to
//! produce the same canonical curve bits.
//!
//! This target makes no claim about sign-seam differentiability, arbitrary
//! points or increments, projection, vector transport, Stiefel manifolds,
//! solver convergence, budgets or cancellation, fs-ascent consumer migration,
//! cross-ISA equality, or performance.

#![deny(unsafe_code)]

use fs_opt::Manifold;

const BASE: [f64; 4] = [0.5, 0.5, 0.5, 0.5];
const BODY_VELOCITY: [f64; 3] = [0.7, -0.4, 0.9];
const COARSE_H: f64 = 1.0 / 512.0;
const FINE_H: f64 = COARSE_H / 2.0;

fn dot<const N: usize>(left: &[f64; N], right: &[f64; N]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn norm<const N: usize>(values: &[f64; N]) -> f64 {
    dot(values, values).sqrt()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn scale_body(scale: f64) -> [f64; 3] {
    BODY_VELOCITY.map(|component| scale * component)
}

/// Independent derivative of `q * exp(h omega / 2)` at `h = 0`.
fn right_curve_velocity(q: &[f64; 4], omega: &[f64; 3]) -> [f64; 4] {
    let [w, x, y, z] = *q;
    let [a, b, c] = *omega;
    [
        -0.5 * (x * a + y * b + z * c),
        0.5 * (w * a + y * c - z * b),
        0.5 * (w * b - x * c + z * a),
        0.5 * (w * c + x * b - y * a),
    ]
}

/// Deliberately wrong-side derivative used to prove the fixture distinguishes
/// body-frame right composition from space-frame left composition.
fn left_curve_velocity(q: &[f64; 4], omega: &[f64; 3]) -> [f64; 4] {
    let [w, x, y, z] = *q;
    let [a, b, c] = *omega;
    [
        -0.5 * (a * x + b * y + c * z),
        0.5 * (a * w + b * z - c * y),
        0.5 * (b * w - a * z + c * x),
        0.5 * (c * w + a * y - b * x),
    ]
}

fn central_difference(base: &[f64; 4], h: f64) -> ([f64; 4], Vec<u64>, Vec<u64>) {
    let manifold = Manifold::So3;
    let plus = manifold
        .retract(base, &scale_body(h))
        .expect("positive curve sample");
    let minus = manifold
        .retract(base, &scale_body(-h))
        .expect("negative curve sample");
    let derivative = core::array::from_fn(|index| (plus[index] - minus[index]) / (2.0 * h));
    (derivative, bits(&plus), bits(&minus))
}

fn max_error(actual: &[f64; 4], expected: &[f64; 4]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max)
}

#[test]
fn g0_declared_right_curve_velocity_is_tangent_and_metric_scaled() {
    let derivative = right_curve_velocity(&BASE, &BODY_VELOCITY);
    let wrong_side = left_curve_velocity(&BASE, &BODY_VELOCITY);

    assert_eq!(norm(&BASE).to_bits(), 1.0f64.to_bits());
    assert!(
        dot(&BASE, &derivative).abs() <= 8.0e-17,
        "right-composed velocity must be tangent: q_dot={derivative:?}"
    );
    assert!(
        (norm(&derivative) - 0.5 * norm(&BODY_VELOCITY)).abs() <= 2.0e-16,
        "unit-quaternion metric must scale body velocity by one half"
    );
    assert!(
        max_error(&derivative, &wrong_side) >= 0.5,
        "fixture must distinguish right from left quaternion composition"
    );
}

#[test]
fn g3_central_difference_refines_toward_independent_curve_velocity() {
    let expected = right_curve_velocity(&BASE, &BODY_VELOCITY);
    let (coarse, _, _) = central_difference(&BASE, COARSE_H);
    let (fine, _, _) = central_difference(&BASE, FINE_H);
    let coarse_error = max_error(&coarse, &expected);
    let fine_error = max_error(&fine, &expected);

    assert!(coarse_error.is_finite() && fine_error.is_finite());
    assert!(
        fine_error < coarse_error * 0.26,
        "central refinement must exhibit its second-order error reduction: coarse={coarse_error:.17e}; fine={fine_error:.17e}"
    );
    assert!(
        fine_error <= 4.0e-8,
        "fine central derivative missed the independent velocity: error={fine_error:.17e}; actual={fine:?}; expected={expected:?}"
    );
    assert!(
        dot(&BASE, &fine).abs() <= 2.0e-13,
        "finite-difference velocity must remain tangent: q_dot={fine:?}"
    );
}

#[test]
fn g3_antipodal_bases_share_canonical_curve_and_derivative_bits() {
    let antipode = BASE.map(|component| -component);
    let (canonical, plus_bits, minus_bits) = central_difference(&BASE, FINE_H);
    let (from_antipode, antipodal_plus_bits, antipodal_minus_bits) =
        central_difference(&antipode, FINE_H);

    assert_eq!(plus_bits, antipodal_plus_bits);
    assert_eq!(minus_bits, antipodal_minus_bits);
    assert_eq!(bits(&canonical), bits(&from_antipode));
}

#[test]
fn g5_fixed_input_curve_samples_replay_bit_for_bit() {
    let first = central_difference(&BASE, FINE_H);
    let second = central_difference(&BASE, FINE_H);

    assert_eq!(bits(&first.0), bits(&second.0));
    assert_eq!(first.1, second.1);
    assert_eq!(first.2, second.2);
}
