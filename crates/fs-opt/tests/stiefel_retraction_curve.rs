//! G0/G3 QR-retraction conformance for the fs-opt Stiefel authority.
//!
//! A fixed `St(4,2)` frame checks zero-step identity, orthonormal landing,
//! positive-diagonal QR sign choice, signed-permutation row covariance, and the
//! differential of the public retraction along an independently constructed
//! tangent direction.
//!
//! This target makes no claim about general QR stability, arbitrary `n`/`p`,
//! polar retractions, right-column covariance, vector transport, metric
//! isometry, solver convergence, budgets or cancellation, fs-ascent migration,
//! cross-ISA equality, or performance.

#![deny(unsafe_code)]

use fs_opt::Manifold;

const N: usize = 4;
const P: usize = 2;
const STORAGE: usize = N * P;
const MANIFOLD: Manifold = Manifold::Stiefel { n: 4, p: 2 };

const BASE: [f64; STORAGE] = [
    0.5, 0.5, 0.5, 0.5, // first column
    0.5, -0.5, 0.5, -0.5, // second column
];
const QR_STEP: [f64; STORAGE] = [
    0.0, 0.0, 0.25, -0.25, // first candidate column increment
    0.125, 0.25, -0.125, 0.0, // second candidate column increment
];
const TANGENT: [f64; STORAGE] = [
    0.4375, 0.0625, -0.0625, -0.4375, // 3/8 x1 + 1/2 z0
    -0.0625, -0.3125, -0.3125, -0.0625, // -3/8 x0 + 1/4 z1
];
const COARSE_H: f64 = 1.0 / 256.0;
const FINE_H: f64 = COARSE_H / 2.0;

fn column(values: &[f64], index: usize) -> &[f64] {
    &values[index * N..(index + 1) * N]
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
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

fn orthonormal_residual(frame: &[f64]) -> f64 {
    let mut residual = 0.0f64;
    for current in 0..P {
        for against in 0..=current {
            let expected = if current == against { 1.0 } else { 0.0 };
            residual = residual
                .max((dot(column(frame, current), column(frame, against)) - expected).abs());
        }
    }
    residual
}

fn tangent_residual(base: &[f64], tangent: &[f64]) -> f64 {
    let mut residual = 0.0f64;
    for row in 0..P {
        for col in 0..P {
            let symmetric = dot(column(base, row), column(tangent, col))
                + dot(column(tangent, row), column(base, col));
            residual = residual.max(symmetric.abs());
        }
    }
    residual
}

fn candidate(base: &[f64], step: &[f64]) -> Vec<f64> {
    base.iter().zip(step).map(|(x, t)| x + t).collect()
}

fn qr_factor(frame: &[f64], source: &[f64]) -> [[f64; P]; P] {
    core::array::from_fn(|row| {
        core::array::from_fn(|col| dot(column(frame, row), column(source, col)))
    })
}

fn qr_certificate(frame: &[f64], source: &[f64]) -> bool {
    let factor = qr_factor(frame, source);
    let triangular = factor[1][0].abs() <= 2.0e-15;
    let positive_diagonal = factor[0][0] > 0.0 && factor[1][1] > 0.0;
    let mut reconstructed = [0.0; STORAGE];
    for col in 0..P {
        for row in 0..N {
            reconstructed[col * N + row] = (0..P)
                .map(|basis| frame[basis * N + row] * factor[basis][col])
                .sum();
        }
    }
    triangular && positive_diagonal && max_error(&reconstructed, source) <= 3.0e-15
}

fn signed_row_permutation(frame: &[f64]) -> Vec<f64> {
    let mut transformed = Vec::with_capacity(STORAGE);
    for col in 0..P {
        let source = column(frame, col);
        transformed.extend([source[2], -source[0], source[3], -source[1]]);
    }
    transformed
}

fn scale_tangent(scale: f64) -> Vec<f64> {
    TANGENT.iter().map(|component| scale * component).collect()
}

fn central_difference(h: f64) -> Vec<f64> {
    let plus = MANIFOLD
        .retract(&BASE, &scale_tangent(h))
        .expect("positive Stiefel curve sample");
    let minus = MANIFOLD
        .retract(&BASE, &scale_tangent(-h))
        .expect("negative Stiefel curve sample");
    plus.iter()
        .zip(&minus)
        .map(|(positive, negative)| (positive - negative) / (2.0 * h))
        .collect()
}

#[test]
fn g0_zero_step_is_identity_for_exact_orthonormal_frame() {
    assert_eq!(orthonormal_residual(&BASE).to_bits(), 0.0f64.to_bits());
    let landed = MANIFOLD
        .retract(&BASE, &[0.0; STORAGE])
        .expect("zero-step Stiefel retraction");
    assert_eq!(bits(&landed), bits(&BASE));
}

#[test]
fn g0_qr_landing_has_positive_diagonal_and_independent_reconstruction() {
    let source = candidate(&BASE, &QR_STEP);
    let landed = MANIFOLD
        .retract(&BASE, &QR_STEP)
        .expect("nonorthogonal Stiefel candidate");
    let replay = MANIFOLD
        .retract(&BASE, &QR_STEP)
        .expect("deterministic Stiefel QR replay");

    assert!(orthonormal_residual(&landed) <= 2.0e-15);
    assert!(qr_certificate(&landed, &source));
    assert_eq!(bits(&landed), bits(&replay));

    let mut wrong_sign = landed.clone();
    for value in &mut wrong_sign[..N] {
        *value = -*value;
    }
    assert!(
        !qr_certificate(&wrong_sign, &source),
        "negative-column mutation must fail the positive-diagonal QR certificate"
    );
}

#[test]
fn g3_signed_row_permutation_commutes_with_qr_retraction() {
    let landed = MANIFOLD
        .retract(&BASE, &QR_STEP)
        .expect("reference Stiefel QR landing");
    let transformed_base = signed_row_permutation(&BASE);
    let transformed_step = signed_row_permutation(&QR_STEP);
    let transformed_landing = MANIFOLD
        .retract(&transformed_base, &transformed_step)
        .expect("signed-permuted Stiefel QR landing");
    let expected = signed_row_permutation(&landed);

    assert!(orthonormal_residual(&transformed_landing) <= 2.0e-15);
    assert!(max_error(&transformed_landing, &expected) <= 3.0e-15);
    assert!(qr_certificate(
        &transformed_landing,
        &candidate(&transformed_base, &transformed_step)
    ));
}

#[test]
fn g3_central_difference_refines_toward_independent_tangent() {
    assert_eq!(
        tangent_residual(&BASE, &TANGENT).to_bits(),
        0.0f64.to_bits()
    );
    let coarse = central_difference(COARSE_H);
    let fine = central_difference(FINE_H);
    let coarse_error = max_error(&coarse, &TANGENT);
    let fine_error = max_error(&fine, &TANGENT);
    let coarse_tangent_residual = tangent_residual(&BASE, &coarse);
    let fine_tangent_residual = tangent_residual(&BASE, &fine);

    assert!(coarse_error.is_finite() && fine_error.is_finite());
    assert!(
        fine_error < coarse_error * 0.27,
        "central refinement must expose second-order QR-curve convergence: coarse={coarse_error:.17e}; fine={fine_error:.17e}"
    );
    assert!(
        fine_error <= 1.0e-6,
        "fine QR-curve derivative missed the independent tangent: error={fine_error:.17e}; actual={fine:?}; expected={TANGENT:?}"
    );
    assert!(fine_tangent_residual < coarse_tangent_residual * 0.27);
    assert!(
        fine_tangent_residual <= 2.0e-7,
        "finite-difference tangent residual stayed too large: coarse={coarse_tangent_residual:.17e}; fine={fine_tangent_residual:.17e}; derivative={fine:?}"
    );
}
