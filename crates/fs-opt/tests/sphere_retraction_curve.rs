//! G0/G3 normalization-retraction conformance for the fs-opt sphere authority.
//!
//! A fixed unit point and non-tangent ambient increment check the independently
//! derived curve velocity `v - x (x dot v)`, signed-permutation covariance,
//! exact zero-step identity, deterministic replay, and typed singular refusal.
//!
//! This target makes no claim about geodesic/exponential retractions, arbitrary
//! ambient dimensions, points, or steps, vector transport, curve acceleration,
//! solver convergence, budgets or cancellation, fs-ascent migration, cross-ISA
//! equality, or performance.

#![deny(unsafe_code)]

use fs_opt::{Manifold, OptError};

const AMBIENT: usize = 4;
const MANIFOLD: Manifold = Manifold::Sphere { ambient: 4 };
const BASE: [f64; AMBIENT] = [0.5; AMBIENT];
const AMBIENT_VELOCITY: [f64; AMBIENT] = [0.75, -0.25, 0.5, -0.5];
const COARSE_H: f64 = 1.0 / 256.0;
const FINE_H: f64 = COARSE_H / 2.0;

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn norm(values: &[f64]) -> f64 {
    dot(values, values).sqrt()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn max_error(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

/// Independent differential of `(x + h v) / ||x + h v||` at `h = 0`
/// for unit `x`.
fn normalization_curve_velocity(
    base: &[f64; AMBIENT],
    ambient_velocity: &[f64; AMBIENT],
) -> [f64; AMBIENT] {
    let radial = dot(base, ambient_velocity);
    core::array::from_fn(|index| ambient_velocity[index] - radial * base[index])
}

fn scale_velocity(scale: f64) -> [f64; AMBIENT] {
    AMBIENT_VELOCITY.map(|component| scale * component)
}

fn central_difference(h: f64) -> ([f64; AMBIENT], Vec<u64>, Vec<u64>) {
    let plus = MANIFOLD
        .retract(&BASE, &scale_velocity(h))
        .expect("positive sphere curve sample");
    let minus = MANIFOLD
        .retract(&BASE, &scale_velocity(-h))
        .expect("negative sphere curve sample");
    let derivative = core::array::from_fn(|index| (plus[index] - minus[index]) / (2.0 * h));
    (derivative, bits(&plus), bits(&minus))
}

fn signed_permutation(values: &[f64; AMBIENT]) -> [f64; AMBIENT] {
    [values[2], -values[0], values[3], -values[1]]
}

#[test]
fn g0_zero_step_is_identity_and_radial_cancellation_is_refused() {
    assert_eq!(norm(&BASE).to_bits(), 1.0f64.to_bits());
    let landed = MANIFOLD
        .retract(&BASE, &[0.0; AMBIENT])
        .expect("zero-step sphere retraction");
    assert_eq!(bits(&landed), bits(&BASE));

    let radial_cancellation = BASE.map(|component| -component);
    assert!(matches!(
        MANIFOLD.retract(&BASE, &radial_cancellation),
        Err(OptError::RetractionDomain {
            manifold: "Sphere",
            what: "candidate norm squared must be finite and nonsingular",
            location: None,
            measurement_bits,
        }) if measurement_bits == 0.0f64.to_bits()
    ));
}

#[test]
fn g0_independent_curve_velocity_removes_exact_radial_component() {
    let derivative = normalization_curve_velocity(&BASE, &AMBIENT_VELOCITY);
    let radial = dot(&BASE, &AMBIENT_VELOCITY);

    assert_eq!(radial.to_bits(), 0.25f64.to_bits());
    assert_eq!(bits(&derivative), bits(&[0.625, -0.375, 0.375, -0.625]));
    assert_eq!(dot(&BASE, &derivative).to_bits(), 0.0f64.to_bits());
    assert!(
        max_error(&derivative, &AMBIENT_VELOCITY) >= 0.125,
        "fixture must distinguish normalized-curve projection from a raw ambient step"
    );
}

#[test]
fn g3_central_difference_refines_toward_projected_curve_velocity() {
    let expected = normalization_curve_velocity(&BASE, &AMBIENT_VELOCITY);
    let (coarse, _, _) = central_difference(COARSE_H);
    let (fine, _, _) = central_difference(FINE_H);
    let coarse_error = max_error(&coarse, &expected);
    let fine_error = max_error(&fine, &expected);
    let coarse_tangent_residual = dot(&BASE, &coarse).abs();
    let fine_tangent_residual = dot(&BASE, &fine).abs();

    assert!(coarse_error.is_finite() && fine_error.is_finite());
    assert!(
        fine_error < coarse_error * 0.26,
        "central refinement must expose second-order normalization-curve convergence: coarse={coarse_error:.17e}; fine={fine_error:.17e}"
    );
    assert!(fine_error <= 2.0e-6);
    assert!(fine_tangent_residual < coarse_tangent_residual * 0.26);
    assert!(fine_tangent_residual <= 1.2e-6);
}

#[test]
fn g3_signed_permutation_commutes_with_normalization_retraction() {
    let step = scale_velocity(0.125);
    let landed = MANIFOLD
        .retract(&BASE, &step)
        .expect("reference sphere landing");
    let transformed_base = signed_permutation(&BASE);
    let transformed_step = signed_permutation(&step);
    let transformed_landing = MANIFOLD
        .retract(&transformed_base, &transformed_step)
        .expect("signed-permuted sphere landing");
    let landed: [f64; AMBIENT] = landed
        .try_into()
        .expect("sphere landing retains four components");
    let expected = signed_permutation(&landed);

    assert!((norm(&transformed_landing) - 1.0).abs() <= 2.0e-15);
    assert!(max_error(&transformed_landing, &expected) <= 2.0e-15);
}

#[test]
fn g3_fixed_curve_samples_replay_bit_for_bit() {
    let first = central_difference(FINE_H);
    let second = central_difference(FINE_H);
    assert_eq!(bits(&first.0), bits(&second.0));
    assert_eq!(first.1, second.1);
    assert_eq!(first.2, second.2);
}
