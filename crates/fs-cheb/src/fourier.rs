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
    /// If `n` is not a power of two ≥ 2 (caller contract — the current
    /// real transform backend is radix-2).
    #[must_use]
    pub fn build<F: Fn(f64) -> f64>(f: &F, n: usize) -> FourierSeries {
        assert!(
            n >= 2 && n.is_power_of_two(),
            "power-of-two sample count >= 2 required"
        );
        let samples: Vec<f64> = (0..n)
            .map(|k| {
                let y = f(std::f64::consts::TAU * k as f64 / n as f64);
                assert!(y.is_finite(), "Fourier samples must be finite");
                y
            })
            .collect();
        let fft = RealFft::new(n);
        let spec = fft.forward(&samples);
        let scale = 1.0 / n as f64;
        let coeffs: Vec<C64> = spec
            .into_iter()
            .map(|c| C64::new(c.re * scale, c.im * scale))
            .collect();
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
    /// A shared real-FFT plan evaluates both inputs on the smallest
    /// power-of-two grid strictly above twice the sum of their degrees,
    /// then transforms their pointwise product. Strict oversampling also
    /// resolves a sine at the product's highest frequency: placing that
    /// frequency at the grid's Nyquist would erase it. Input endpoint
    /// coefficients receive their half weight when embedded as interior
    /// modes. All sum/difference harmonics are retained, with no threshold
    /// truncation. Modes above the known product degree are set to zero.
    ///
    /// Work is O(N log N), temporary storage O(N), where N is the output
    /// resolution. Constant multiplication uses coefficient scaling.
    /// This is algebra of the retained interpolants; it cannot recover
    /// frequencies already aliased when the inputs were constructed.
    ///
    /// # Panics
    /// If the required grid size overflows, or coefficients, transform
    /// intermediates, or sample products are not representable as finite f64.
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
        let fft = RealFft::new(n);
        let mut products = self.padded_samples(&fft, n, left_degree);
        let right = other.padded_samples(&fft, n, right_degree);
        for (left, right) in products.iter_mut().zip(right) {
            *left *= right;
            assert!(left.is_finite(), "Fourier sample product is not finite");
        }
        let scale = 1.0 / n as f64;
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
                    C64::new(c.re * scale, if k == 0 { 0.0 } else { c.im * scale })
                }
            })
            .collect();
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

    fn padded_samples(&self, fft: &RealFft, n: usize, degree: usize) -> Vec<f64> {
        let mut spectrum = vec![FftComplex::default(); n / 2 + 1];
        for (k, c) in self.coeffs.iter().take(degree + 1).enumerate() {
            let scale = if k == self.n / 2 {
                0.5 * n as f64
            } else {
                n as f64
            };
            spectrum[k] = FftComplex::new(c.re * scale, if k == 0 { 0.0 } else { c.im * scale });
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
