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
fn sparse_nonlinear_products_preserve_weak_generated_harmonics_and_exact_zeros() {
    // A small retained partial is independent information, even when adding it
    // to the fundamental on a sample grid would round it away. Product-to-sum:
    // (cos(t) + e sin(5t)) cos(2t)
    //   = (cos(t) + cos(3t) + e sin(3t) + e sin(7t))/2.
    let epsilon = 1.0e-100;
    let f = series(16, 0.0, &[(1, 1.0, 0.0), (5, 0.0, epsilon)]);
    let g = series(16, 0.0, &[(2, 1.0, 0.0)]);
    let product = f.mul(&g);
    assert_eq!(product.coeffs[1], C64::from_re(0.25));
    assert_eq!(product.coeffs[3], C64::new(0.25, -0.25 * epsilon));
    assert_eq!(product.coeffs[7], C64::new(0.0, -0.25 * epsilon));
    for (k, c) in product.coeffs.iter().enumerate() {
        if ![1, 3, 7].contains(&k) {
            assert_eq!(*c, C64::ZERO, "spurious harmonic {k}");
        }
    }
    assert_eq!(product.coeffs, g.mul(&f).coeffs);
}

#[test]
fn sparse_cubic_nonlinearity_does_not_create_forbidden_even_harmonics() {
    // Odd polynomial nonlinearities driven by an odd-harmonic signal retain
    // half-wave symmetry. No FFT roundoff floor should invent even partials.
    let f = series(32, 0.0, &[(1, 1.0, 0.0), (5, 0.125, 0.0)]);
    let cube = f.mul(&f).mul(&f);
    for (k, c) in cube.coeffs.iter().enumerate() {
        if k % 2 == 0 || k > 15 {
            assert_eq!(*c, C64::ZERO, "forbidden harmonic {k}");
        }
    }
    for j in 0..79 {
        let t = std::f64::consts::TAU * (j as f64 + 0.17) / 79.0;
        close(cube.eval(t), (t.cos() + 0.125 * (5.0 * t).cos()).powi(3));
    }
}

#[test]
fn dense_fourier_products_still_retain_all_sum_and_difference_modes() {
    // 32 nonzero modes per input exceed the sparse K*L <= N work bound.
    let modes: Vec<_> = (1..=32)
        .map(|k| (k, 1.0 / (k as f64 + 1.0), 0.25 / k as f64))
        .collect();
    let f = series(64, 0.0, &modes);
    let product = f.mul(&f);
    for j in 0..263 {
        let t = std::f64::consts::TAU * (j as f64 + 0.31) / 263.0;
        close(product.eval(t), f.eval(t).powi(2));
    }
    let (_, a, b) = modes[31];
    close(product.coeffs[64].re, 0.25 * (a * a - b * b));
    close(product.coeffs[64].im, -0.5 * a * b);
}

#[test]
fn unrepresentable_derivative_refuses_instead_of_publishing_infinity() {
    let f = series(8, 0.0, &[(4, f64::MAX, 0.0)]);
    assert!(std::panic::catch_unwind(|| f.differentiate()).is_err());
}

#[test]
fn construction_preserves_constants_across_the_finite_amplitude_range() {
    for amplitude in [f64::MAX, -f64::MAX, f64::MIN_POSITIVE, f64::from_bits(1)] {
        for n in [2, 8, 64] {
            let f = FourierSeries::build(&|_| amplitude, n);
            assert_eq!(f.coeffs[0].re, amplitude, "mean, n={n}");
            assert!(
                f.coeffs
                    .iter()
                    .all(|c| c.re.is_finite() && c.im.is_finite())
            );
            for c in &f.coeffs[1..] {
                assert_eq!(c.re, 0.0);
                assert_eq!(c.im, 0.0);
            }
        }
    }
}

#[test]
fn construction_retains_phase_when_the_unnormalized_spectrum_would_overflow() {
    let amplitude = f64::from_bits(2043_u64 << 52); // 2^1020
    let f = FourierSeries::build(
        &|t| amplitude * (0.5 + 0.75 * (3.0 * t).cos() - 0.25 * (5.0 * t).sin()),
        128,
    );
    close(f.coeffs[0].re / amplitude, 0.5);
    close(f.coeffs[3].re / amplitude, 0.375);
    close(f.coeffs[5].im / amplitude, 0.125);
    for j in 0..79 {
        let t = std::f64::consts::TAU * (j as f64 + 0.17) / 79.0;
        close(
            f.eval(t) / amplitude,
            0.5 + 0.75 * (3.0 * t).cos() - 0.25 * (5.0 * t).sin(),
        );
    }
}

#[test]
fn nonlinear_products_are_invariant_under_reciprocal_amplitude_units() {
    for exponent in [1000, -1000] {
        let large = f64::from_bits(((1023 + exponent) as u64) << 52);
        let small = f64::from_bits(((1023 - exponent) as u64) << 52);
        let f = series(8, 0.0, &[(4, large, -0.5 * large)]);
        let g = series(16, 0.25 * small, &[(3, 0.75 * small, small)]);
        let product = f.mul(&g);
        for j in 0..97 {
            let t = std::f64::consts::TAU * (j as f64 + 0.43) / 97.0;
            close(
                product.eval(t),
                ((4.0 * t).cos() - 0.5 * (4.0 * t).sin())
                    * (0.25 + 0.75 * (3.0 * t).cos() + (3.0 * t).sin()),
            );
        }
        close(product.coeffs[7].re, 0.3125);
        close(product.coeffs[7].im, -0.15625);
    }
}

#[test]
fn nonlinear_product_rescales_small_coefficients_without_losing_subnormals() {
    // (2^-536 cos(4t))² = 2^-1073 (1 + cos(8t)). The mean is
    // two minimum subnormals; the one-sided 8th coefficient is one.
    let amplitude = f64::from_bits((1023 - 536) << 52);
    let f = series(8, 0.0, &[(4, amplitude, 0.0)]);
    let product = f.mul(&f);
    assert_eq!(product.coeffs[0].re, f64::from_bits(2));
    assert_eq!(product.coeffs[8].re, f64::from_bits(1));
}

#[test]
fn unrepresentable_product_coefficient_refuses_after_rescaling() {
    let f = series(8, f64::MAX, &[(4, f64::MAX, 0.0)]);
    let g = series(8, 2.0, &[(4, 1.0, 0.0)]);
    assert!(std::panic::catch_unwind(|| f.mul(&g)).is_err());
}
