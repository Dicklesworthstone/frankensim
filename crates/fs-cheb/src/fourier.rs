//! Fourier-periodic function objects (bead kw89): trigonometric
//! interpolants for closed curves and angular profiles. Coefficients
//! come from fs-fft's real transform on uniform samples of [0, 2π);
//! evaluation sums the (half-)spectrum directly, differentiation is
//! the ik multiply (with the Nyquist mode zeroed for odd derivatives
//! of real signals — the standard spectral convention, documented),
//! and the mean/integral read off c₀.

use fs_fft::RealFft;
use fs_math::{c64::C64, det};

/// A real trigonometric interpolant on [0, 2π), degree n/2 from n
/// uniform samples (n even).
pub struct FourierSeries {
    /// Half-spectrum c_0..c_{n/2} (RealFft layout, scaled by 1/n).
    pub coeffs: Vec<C64>,
    /// Sample count.
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

    /// The spectral derivative: multiply mode k by ik (Nyquist mode
    /// zeroed — the real-signal convention for odd derivatives).
    #[must_use]
    pub fn differentiate(&self) -> FourierSeries {
        self.assert_valid();
        let half = self.n / 2;
        let coeffs: Vec<C64> = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(k, c)| {
                if k == half {
                    C64::ZERO
                } else {
                    let kf = k as f64;
                    // ik·(re + i·im) = −k·im + i·k·re.
                    C64::new(-kf * c.im, kf * c.re)
                }
            })
            .collect();
        FourierSeries { coeffs, n: self.n }
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
