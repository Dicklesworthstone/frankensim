//! Fourier-periodic function objects (bead kw89): trigonometric
//! interpolants for closed curves and angular profiles. Coefficients
//! come from fs-fft's real transform on uniform samples of [0, 2π);
//! evaluation sums the (half-)spectrum directly, differentiation is
//! the continuous ik multiply, and the mean/integral read off c₀.
//! Nonlinear products pad beyond the sum of the input degrees so new
//! harmonics cannot alias into lower modes.

use fs_fft::{C64 as FftComplex, RealFft};
use fs_math::{c64::C64, det};

/// A real trigonometric polynomial on [0, 2π), of degree at most n/2.
/// Construction interpolates n uniform samples; subsequent calculus
/// retains the continuous polynomial, including its highest sine mode.
pub struct FourierSeries {
    /// One-sided coefficients c_0..c_{n/2}. The constant is c_0.re;
    /// interior terms are 2 Re(c_k exp(ikθ)), and the highest term is
    /// Re(c_{n/2} exp(inθ/2)). Construction produces the RealFft layout
    /// scaled by 1/n. Calculus may produce an imaginary highest coefficient:
    /// its sine is zero at the original grid points but nonzero between them.
    pub coeffs: Vec<C64>,
    /// Power-of-two resolution parameter (the construction sample count).
    pub n: usize,
}

impl FourierSeries {
    /// Interpolate `f` from `n` uniform samples (n even, ≥ 2).
    ///
    /// # Panics
    /// If `n` is not a power of two ≥ 2, a sample is non-finite, or a
    /// retained coefficient is not representable as finite f64.
    #[must_use]
    pub fn build<F: Fn(f64) -> f64>(f: &F, n: usize) -> FourierSeries {
        assert!(
            n >= 2 && n.is_power_of_two(),
            "power-of-two sample count >= 2 required"
        );
        let mut samples: Vec<f64> = (0..n)
            .map(|k| {
                let y = f(std::f64::consts::TAU * k as f64 / n as f64);
                assert!(y.is_finite(), "Fourier samples must be finite");
                y
            })
            .collect();
        let coeffs = coefficients_from_samples(&mut samples)
            .expect("Fourier coefficient is not representable as finite f64");
        FourierSeries { coeffs, n }
    }

    fn assert_valid(&self) {
        assert!(
            self.n >= 2 && self.n.is_power_of_two(),
            "FourierSeries sample count must be a power of two >= 2"
        );
        assert_eq!(
            self.coeffs.len(),
            self.n / 2 + 1,
            "FourierSeries coefficient count must be n/2 + 1"
        );
        assert!(
            self.coeffs
                .iter()
                .all(|c| c.re.is_finite() && c.im.is_finite()),
            "FourierSeries coefficients must be finite"
        );
    }

    /// Evaluate at θ (real part of the half-spectrum sum with the
    /// conjugate-symmetric double counting; endpoints k = 0 and
    /// k = n/2 count once).
    #[must_use]
    pub fn eval(&self, theta: f64) -> f64 {
        self.assert_valid();
        assert!(theta.is_finite(), "Fourier evaluation point must be finite");
        let half = self.n / 2;
        let mut s = self.coeffs[0].re;
        for (k, c) in self.coeffs.iter().enumerate().skip(1) {
            let angle = k as f64 * theta;
            let (sin, cos) = (det::sin(angle), det::cos(angle));
            let term = c.re.mul_add(cos, -(c.im * sin));
            s += if k == half { term } else { 2.0 * term };
        }
        s
    }

    /// Differentiate the continuous polynomial by multiplying every mode
    /// k by ik, including the highest mode. In particular, the derivative
    /// of a Nyquist cosine is a sine that vanishes at the construction
    /// grid, but is retained for off-grid evaluation and later derivatives.
    ///
    /// # Panics
    /// If a derivative coefficient is not representable as finite f64.
    #[must_use]
    pub fn differentiate(&self) -> FourierSeries {
        self.assert_valid();
        let coeffs: Vec<C64> = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(k, c)| {
                let kf = k as f64;
                // ik·(re + i·im) = −k·im + i·k·re.
                let derivative = C64::new(-kf * c.im, kf * c.re);
                assert!(
                    derivative.re.is_finite() && derivative.im.is_finite(),
                    "Fourier derivative coefficient is not representable as finite f64"
                );
                derivative
            })
            .collect();
        FourierSeries { coeffs, n: self.n }
    }

    /// Multiply two retained trigonometric polynomials without aliasing.
    ///
    /// Sparse inputs are convolved directly in coefficient space. This keeps
    /// weak generated partials that would disappear when added to stronger
    /// components on a sampling grid, and leaves absent harmonics exactly zero.
    /// Dense inputs use a shared real-FFT plan on the smallest
    /// power-of-two grid strictly above twice the sum of their degrees,
    /// then transforms their pointwise product. Strict oversampling also
    /// resolves a sine at the product's highest frequency: placing that
    /// frequency at the grid's Nyquist would erase it. Input endpoint
    /// coefficients receive their half weight when embedded as interior
    /// modes. All sum/difference harmonics are retained, with no threshold
    /// truncation. Modes above the known product degree are set to zero.
    ///
    /// Work is O(N + K L) for K L <= N, otherwise O(N log N), with O(N)
    /// temporary storage. K and L count nonzero input coefficients, and N is
    /// the output resolution. Constant multiplication uses coefficient scaling.
    /// This is algebra of the retained interpolants; it cannot recover
    /// frequencies already aliased when the inputs were constructed.
    ///
    /// The inputs are independently scaled by exact powers of two during
    /// convolution or transform arithmetic. This keeps the result
    /// invariant under reciprocal changes in the inputs' physical units,
    /// even when their unscaled FFT buffers would overflow or underflow.
    ///
    /// # Panics
    /// If the required grid size overflows or a retained output coefficient
    /// is not representable as finite f64.
    #[must_use]
    pub fn mul(&self, other: &FourierSeries) -> FourierSeries {
        self.assert_valid();
        other.assert_valid();
        let left_degree = self.degree();
        let right_degree = other.degree();
        if left_degree == 0 {
            return other.scaled(self.coeffs[0].re);
        }
        if right_degree == 0 {
            return self.scaled(other.coeffs[0].re);
        }
        let degree = left_degree
            .checked_add(right_degree)
            .expect("Fourier product degree overflows usize");
        let n = degree
            .checked_mul(2)
            .and_then(|size| size.checked_add(1))
            .and_then(usize::checked_next_power_of_two)
            .expect("Fourier product grid size overflows usize");
        let left_exponent = self.amplitude_exponent();
        let right_exponent = other.amplitude_exponent();
        let count = |series: &Self| {
            series
                .coeffs
                .iter()
                .filter(|c| c.re != 0.0 || c.im != 0.0)
                .count()
        };
        if (count(self) as u128) * (count(other) as u128) <= n as u128 {
            return self.sparse_product(other, n, left_exponent, right_exponent);
        }
        let fft = RealFft::new(n);
        let mut products = self.padded_samples(&fft, n, left_degree, left_exponent);
        let right = other.padded_samples(&fft, n, right_degree, right_exponent);
        for (left, right) in products.iter_mut().zip(right) {
            *left *= right;
            assert!(left.is_finite(), "Fourier sample product is not finite");
        }
        let exponent = left_exponent + right_exponent - n.ilog2() as i32;
        let coeffs = fft
            .forward(&products)
            .into_iter()
            .enumerate()
            .map(|(k, c)| {
                assert!(
                    c.re.is_finite() && c.im.is_finite(),
                    "Fourier product transform is not finite"
                );
                if k > degree {
                    C64::ZERO
                } else {
                    let coefficient = C64::new(
                        scale_pow2(c.re, exponent),
                        if k == 0 {
                            0.0
                        } else {
                            scale_pow2(c.im, exponent)
                        },
                    );
                    assert!(
                        coefficient.re.is_finite() && coefficient.im.is_finite(),
                        "Fourier product coefficient is not representable as finite f64"
                    );
                    coefficient
                }
            })
            .collect();
        FourierSeries { coeffs, n }
    }

    // Convolve the conjugate-symmetric full spectrum without constructing its
    // negative half. Interior input coefficients already have full-spectrum
    // weight; the original Nyquist coefficient needs a factor of one half.
    // Normalizing before that factor preserves representable subnormal results.
    fn sparse_product(
        &self,
        other: &Self,
        n: usize,
        left_exponent: i32,
        right_exponent: i32,
    ) -> FourierSeries {
        let modes = |series: &Self, exponent: i32| -> Vec<(usize, C64)> {
            series
                .coeffs
                .iter()
                .enumerate()
                .filter(|(_, c)| c.re != 0.0 || c.im != 0.0)
                .map(|(k, c)| {
                    let exponent = -exponent - i32::from(k == series.n / 2);
                    (
                        k,
                        C64::new(
                            scale_pow2(c.re, exponent),
                            if k == 0 {
                                0.0
                            } else {
                                scale_pow2(c.im, exponent)
                            },
                        ),
                    )
                })
                .collect()
        };
        let left = modes(self, left_exponent);
        let right = modes(other, right_exponent);
        let mut coeffs = vec![C64::ZERO; n / 2 + 1];
        let mut errors = vec![C64::ZERO; n / 2 + 1];
        for &(j, a) in &left {
            for &(k, b) in &right {
                let sum = a * b;
                accumulate(&mut coeffs[j + k], &mut errors[j + k], sum);
                if j != 0 && k != 0 {
                    let difference = a * b.conj();
                    let difference = if j == k {
                        C64::from_re(2.0 * difference.re)
                    } else if j < k {
                        difference.conj()
                    } else {
                        difference
                    };
                    let index = j.abs_diff(k);
                    accumulate(&mut coeffs[index], &mut errors[index], difference);
                }
            }
        }
        let exponent = left_exponent + right_exponent;
        for (coefficient, error) in coeffs.iter_mut().zip(errors) {
            coefficient.re = scale_pow2(coefficient.re + error.re, exponent);
            coefficient.im = scale_pow2(coefficient.im + error.im, exponent);
            assert!(
                coefficient.re.is_finite() && coefficient.im.is_finite(),
                "Fourier product coefficient is not representable as finite f64"
            );
        }
        FourierSeries { coeffs, n }
    }

    // Only exact zeros reduce the retained degree. A numerical plateau is
    // not permission to discard a small harmonic during nonlinear algebra.
    fn degree(&self) -> usize {
        self.coeffs
            .iter()
            .rposition(|c| c.re != 0.0 || c.im != 0.0)
            .unwrap_or(0)
    }

    fn scaled(&self, scale: f64) -> FourierSeries {
        let coeffs = self
            .coeffs
            .iter()
            .map(|c| {
                let c = c.scale(scale);
                assert!(
                    c.re.is_finite() && c.im.is_finite(),
                    "Fourier product coefficient is not representable as finite f64"
                );
                c
            })
            .collect();
        FourierSeries { coeffs, n: self.n }
    }

    fn amplitude_exponent(&self) -> i32 {
        binary_exponent(self.coeffs.iter().fold(0.0_f64, |largest, c| {
            largest.max(c.re.abs()).max(c.im.abs())
        }))
    }

    fn padded_samples(
        &self,
        fft: &RealFft,
        n: usize,
        degree: usize,
        amplitude_exponent: i32,
    ) -> Vec<f64> {
        let mut spectrum = vec![FftComplex::default(); n / 2 + 1];
        for (k, c) in self.coeffs.iter().take(degree + 1).enumerate() {
            let exponent = n.ilog2() as i32 - amplitude_exponent - i32::from(k == self.n / 2);
            spectrum[k] = FftComplex::new(
                scale_pow2(c.re, exponent),
                if k == 0 {
                    0.0
                } else {
                    scale_pow2(c.im, exponent)
                },
            );
            assert!(
                spectrum[k].re.is_finite() && spectrum[k].im.is_finite(),
                "Fourier padded coefficient is not representable as finite f64"
            );
        }
        let samples = fft.inverse(&spectrum);
        assert!(
            samples.iter().all(|x| x.is_finite()),
            "Fourier interpolant samples are not finite"
        );
        samples
    }

    /// ∫₀^{2π} f dθ = 2π·c₀.
    #[must_use]
    pub fn integral(&self) -> f64 {
        self.assert_valid();
        std::f64::consts::TAU * self.coeffs[0].re
    }

    /// Largest |c_k| for k ≥ `from` (the spectral-decay ledger probe).
    #[must_use]
    pub fn tail_magnitude(&self, from: usize) -> f64 {
        self.assert_valid();
        self.coeffs
            .iter()
            .skip(from)
            .fold(0.0f64, |m, c| m.max(c.abs()))
    }
}

// Componentwise Neumaier accumulation retains small intermodulation terms
// through cancellation of stronger contributions to the same harmonic.
fn accumulate(sum: &mut C64, error: &mut C64, value: C64) {
    fn scalar(sum: &mut f64, error: &mut f64, value: f64) {
        let next = *sum + value;
        *error += if sum.abs() >= value.abs() {
            (*sum - next) + value
        } else {
            (value - next) + *sum
        };
        *sum = next;
    }
    scalar(&mut sum.re, &mut error.re, value.re);
    scalar(&mut sum.im, &mut error.im, value.im);
}

// Shared with the budgeted constructor: scale in the existing sample buffer,
// not in an extra allocation. All FFT arithmetic then operates on magnitudes
// of order one; restore the physical amplitude together with the 1/n factor.
// Power-of-two scaling preserves the rounding of ordinary normal-range values.
pub(crate) fn coefficients_from_samples(samples: &mut [f64]) -> Option<Vec<C64>> {
    let n = samples.len();
    let exponent = binary_exponent(samples.iter().fold(0.0_f64, |m, x| m.max(x.abs())));
    if exponent != 0 {
        for sample in samples.iter_mut() {
            *sample = scale_pow2(*sample, -exponent);
        }
    }
    let fft = RealFft::new(n);
    let exponent = exponent - n.ilog2() as i32;
    fft.forward(samples)
        .into_iter()
        .enumerate()
        .map(|(k, c)| {
            // Real samples constrain both endpoints to be real. In particular,
            // roundoff in sin(pi) must not invent a continuous Nyquist sine:
            // only subsequent calculus or an explicit coefficient supplies it.
            let im = if k == 0 || k == n / 2 {
                0.0
            } else {
                scale_pow2(c.im, exponent)
            };
            let c = C64::new(scale_pow2(c.re, exponent), im);
            (c.re.is_finite() && c.im.is_finite()).then_some(c)
        })
        .collect()
}

// floor(log2(x)) without a logarithm or reciprocal, including subnormals.
// Callers supply a nonnegative finite magnitude; zero needs no scaling.
fn binary_exponent(x: f64) -> i32 {
    let bits = x.to_bits();
    let encoded = (bits >> 52) as i32;
    if encoded != 0 {
        encoded - 1023
    } else if bits != 0 {
        63 - bits.leading_zeros() as i32 - 1074
    } else {
        0
    }
}

// Multiply by 2^exponent without first materializing an infinite/zero factor.
// On a downscale, leave the largest negative step until last so a nonzero
// subnormal result is rounded only once, at the final multiplication.
fn scale_pow2(mut value: f64, mut exponent: i32) -> f64 {
    while exponent > 1023 {
        value *= f64::from_bits(2046_u64 << 52);
        exponent -= 1023;
    }
    while exponent < -1022 {
        let step = (exponent + 1022).max(-1022);
        value *= f64::from_bits(((1023 + step) as u64) << 52);
        exponent -= step;
    }
    value * f64::from_bits(((1023 + exponent) as u64) << 52)
}
