//! Continuous Fourier calculus and nonlinear generation of retained harmonics.

use fs_cheb::fourier::FourierSeries;
use fs_math::c64::C64;

fn series(n: usize, mean: f64, modes: &[(usize, f64, f64)]) -> FourierSeries {
    let mut coeffs = vec![C64::ZERO; n / 2 + 1];
    coeffs[0] = C64::from_re(mean);
    for &(k, cosine, sine) in modes {
        let weight = if k == n / 2 { 1.0 } else { 0.5 };
        coeffs[k] = C64::new(weight * cosine, -weight * sine);
    }
    FourierSeries { coeffs, n }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-12 * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn nyquist_cosine_keeps_off_grid_derivative_and_laplacian() {
    let f = FourierSeries::build(&|theta| (8.0 * theta).cos(), 16);
    let derivative = f.differentiate();
    let laplacian = derivative.differentiate();
    for j in 0..97 {
        let theta = std::f64::consts::TAU * (j as f64 + 0.37) / 97.0;
        close(derivative.eval(theta), -8.0 * (8.0 * theta).sin());
        close(laplacian.eval(theta), -64.0 * (8.0 * theta).cos());
    }
    close(derivative.integral(), 0.0);
}

#[test]
fn imaginary_highest_coefficient_is_a_continuous_sine() {
    let f = series(16, 0.0, &[(8, 0.0, 1.25)]);
    let derivative = f.differentiate();
    for j in 0..97 {
        let theta = std::f64::consts::TAU * (j as f64 + 0.23) / 97.0;
        close(f.eval(theta), 1.25 * (8.0 * theta).sin());
        close(derivative.eval(theta), 10.0 * (8.0 * theta).cos());
    }
}

#[test]
fn nyquist_square_generates_second_harmonic_without_folding_into_mean() {
    let f = series(16, 0.0, &[(8, 1.0, 0.0)]);
    let square = f.mul(&f);
    let aliased = FourierSeries::build(&|theta| (8.0 * theta).cos().powi(2), 16);
    close(aliased.integral() / std::f64::consts::TAU, 1.0);
    close(square.integral() / std::f64::consts::TAU, 0.5);
    close(square.coeffs[16].re, 0.25);
    close(square.coeffs[16].im, 0.0);
    for j in 0..97 {
        let theta = std::f64::consts::TAU * (j as f64 + 0.19) / 97.0;
        close(square.eval(theta), 0.5 + 0.5 * (16.0 * theta).cos());
    }
}

#[test]
fn differentiated_nyquist_products_retain_highest_sine_and_cosine() {
    let cosine = series(16, 0.0, &[(8, 1.0, 0.0)]);
    let derivative = cosine.differentiate();
    let mixed = cosine.mul(&derivative);
    let square = derivative.mul(&derivative);
    close(mixed.coeffs[16].im, 2.0);
    close(square.integral() / std::f64::consts::TAU, 32.0);
    for j in 0..97 {
        let theta = std::f64::consts::TAU * (j as f64 + 0.41) / 97.0;
        close(mixed.eval(theta), -4.0 * (16.0 * theta).sin());
        close(square.eval(theta), 32.0 - 32.0 * (16.0 * theta).cos());
    }
}

#[test]
fn unequal_grids_retain_sum_and_difference_harmonics_and_phase() {
    let f = series(8, 0.5, &[(1, 0.75, 0.0), (4, 1.5, -0.25)]);
    let g = series(32, -0.75, &[(3, 0.0, 0.125), (16, 0.5, 1.0)]);
    let product = f.mul(&g);
    // Product-to-sum identities for the two highest cosine/sine pairs.
    close(product.coeffs[20].re, 0.25);
    close(product.coeffs[20].im, -0.34375);
    close(product.coeffs[12].re, 0.125);
    close(product.coeffs[12].im, -0.40625);
    for j in 0..113 {
        let theta = std::f64::consts::TAU * (j as f64 + 0.29) / 113.0;
        let left =
            0.5 + 0.75 * theta.cos() + 1.5 * (4.0 * theta).cos() - 0.25 * (4.0 * theta).sin();
        let right =
            -0.75 + 0.125 * (3.0 * theta).sin() + 0.5 * (16.0 * theta).cos() + (16.0 * theta).sin();
        close(product.eval(theta), left * right);
    }
}

#[test]
fn nonlinear_harmonics_obey_continuous_product_rule() {
    let f = series(16, 0.25, &[(3, 0.5, -0.75), (8, 1.0, 0.5)]);
    let g = series(8, -0.5, &[(1, 0.25, 0.0), (4, 0.75, 1.0)]);
    let derivative = f.mul(&g).differentiate();
    for j in 0..131 {
        let t = std::f64::consts::TAU * (j as f64 + 0.31) / 131.0;
        let f = 0.25 + 0.5 * (3.0 * t).cos() - 0.75 * (3.0 * t).sin()
            + (8.0 * t).cos()
            + 0.5 * (8.0 * t).sin();
        let df = -1.5 * (3.0 * t).sin() - 2.25 * (3.0 * t).cos() - 8.0 * (8.0 * t).sin()
            + 4.0 * (8.0 * t).cos();
        let g = -0.5 + 0.25 * t.cos() + 0.75 * (4.0 * t).cos() + (4.0 * t).sin();
        let dg = -0.25 * t.sin() - 3.0 * (4.0 * t).sin() + 4.0 * (4.0 * t).cos();
        close(derivative.eval(t), df * g + f * dg);
    }
}

#[test]
fn nonlinear_mean_agrees_with_continuous_parseval_energy() {
    let modes = [(1, 0.2, -0.4), (3, 0.75, 1.5), (8, 1.25, -0.5)];
    let f = series(16, 0.7, &modes);
    let energy = 0.7 * 0.7 + 0.5 * modes.iter().map(|(_, a, b)| a * a + b * b).sum::<f64>();
    close(f.mul(&f).integral() / std::f64::consts::TAU, energy);
}

#[test]
fn exact_padding_does_not_inflate_repeated_products_and_constants_preserve_small_modes() {
    let padded = series(1024, 0.0, &[(2, 1.0, 0.0)]);
    let square = padded.mul(&padded);
    let cube = square.mul(&padded);
    assert_eq!(square.n, 16);
    assert_eq!(cube.n, 16);
    for j in 0..53 {
        let t = std::f64::consts::TAU * (j as f64 + 0.47) / 53.0;
        close(
            cube.eval(t),
            0.75 * (2.0 * t).cos() + 0.25 * (6.0 * t).cos(),
        );
    }

    let weak = series(8, 1.0, &[(4, 1.0e-100, -1.0e-120)]);
    let constant = series(2, -2.0, &[]);
    let scaled = weak.mul(&constant);
    assert_eq!(scaled.coeffs[4], C64::new(-2.0e-100, -2.0e-120));
    assert_eq!(scaled.n, weak.n);
    close(scaled.integral() / std::f64::consts::TAU, -2.0);
}

#[test]
fn unrepresentable_derivative_refuses_instead_of_publishing_infinity() {
    let f = series(8, 0.0, &[(4, f64::MAX, 0.0)]);
    assert!(std::panic::catch_unwind(|| f.differentiate()).is_err());
}
