//! Prepared exact-length transforms for high-order AFT waveforms.
//!
//! The original direct, compensated path remains useful at small orders. The
//! FFT route preserves each actual even or odd sampling grid: Bluestein pads
//! only its internal convolution, never the physical quadrature or waveform.

use super::{AftSum, OrbitError, PhaseBasis};
use fs_fft::{C64 as FftComplex, GeneralFft, GeneralFftBudget};
use fs_math::c64::C64;

pub(super) fn enabled(samples: usize, harmonics: usize) -> bool {
    // The odd-grid Bluestein convolution needs two larger FFTs per call.
    // Native crossover measurements therefore select it only at higher order.
    harmonics >= if samples.is_power_of_two() { 32 } else { 256 }
}

/// Additional numeric payload beyond the existing phase and sample arrays.
/// Includes the FFT owner's preparation peak, plan and scratch, one transform
/// record and a complete retained projection. No per-state plans are retained.
pub(super) fn payload(
    samples: usize,
    harmonics: usize,
    dimension: usize,
) -> Result<usize, OrbitError> {
    if !enabled(samples, harmonics) {
        return Ok(0);
    }
    let plan = GeneralFft::required_workspace_bytes(samples).map_err(|_| invalid())?;
    let coefficients = harmonics
        .checked_add(1)
        .and_then(|n| n.checked_mul(dimension));
    coefficients
        .and_then(|n| n.checked_add(samples))
        .and_then(|n| n.checked_mul(core::mem::size_of::<C64>()))
        .and_then(|n| n.checked_add(plan))
        .filter(|&n| n <= isize::MAX as usize)
        .ok_or_else(invalid)
}

fn invalid() -> OrbitError {
    OrbitError::BadParameter {
        what: "exact-length AFT transform payload is unrepresentable",
    }
}

pub(super) struct Transform {
    plan: GeneralFft,
    scratch: Vec<FftComplex>,
    data: Vec<FftComplex>,
    pub(super) coefficients: Vec<C64>,
    pub(super) projected: bool,
}

impl Transform {
    pub(super) fn new(
        samples: usize,
        harmonics: usize,
        dimension: usize,
    ) -> Result<Self, OrbitError> {
        let bytes = payload(samples, harmonics, dimension)?;
        let plan = GeneralFft::new(
            samples,
            GeneralFftBudget {
                max_samples: samples,
                max_workspace_bytes: bytes,
            },
        )
        .map_err(|_| invalid())?;
        Ok(Self {
            scratch: vec![FftComplex::default(); plan.scratch_len()],
            data: vec![FftComplex::default(); samples],
            coefficients: vec![C64::ZERO; (harmonics + 1) * dimension],
            projected: false,
            plan,
        })
    }

    /// An unnormalized forward DFT of conjugated Fourier bins reconstructs
    /// the real AC waveform. Neither multiplication nor division by M is
    /// needed; a representable physical coefficient is never first scaled by M.
    fn reconstruct(&mut self, harmonics: usize, coefficient: impl Fn(usize) -> C64) -> bool {
        self.data.fill(FftComplex::default());
        let samples = self.data.len();
        for n in 1..=harmonics {
            let c = coefficient(n);
            self.data[n] = FftComplex::new(c.re, -c.im);
            self.data[samples - n] = FftComplex::new(c.re, c.im);
        }
        self.plan.forward(&mut self.data, &mut self.scratch).is_ok()
    }

    /// Samples 1..M contain (g_j-g_0)/M; sample zero is the unscaled reference.
    /// Transform only deviations, so a constant physical source has exact-zero
    /// AC even when its preload is many orders larger than its harmonics.
    pub(super) fn project(&mut self, centered: &[f64], dimension: usize, harmonics: usize) {
        self.projected = false;
        for k in 0..dimension {
            self.data[0] = FftComplex::default();
            let mut dc = AftSum::new(centered[k]);
            for j in 1..self.data.len() {
                let value = centered[j * dimension + k];
                self.data[j] = FftComplex::new(value, 0.0);
                dc.add(value);
            }
            if self
                .plan
                .forward(&mut self.data, &mut self.scratch)
                .is_err()
            {
                // An FFT's intermediate dynamic range can be narrower than
                // the final compensated projection. The original direct owner
                // remains available; no source callback or equation is changed.
                return;
            }
            self.coefficients[k] = C64::from_re(dc.total());
            for n in 1..=harmonics {
                let value = self.data[n];
                self.coefficients[n * dimension + k] = C64::new(value.re, value.im);
            }
        }
        self.projected = true;
    }
}

/// Visit reconstructed AC samples, with identical original coefficient units.
/// Callers add/cancel DC separately to preserve small physical residuals under
/// a large static offset. The visitor must not depend on sample traversal order.
pub(super) fn reconstruct(
    basis: &PhaseBasis,
    mut transform: Option<&mut Transform>,
    samples: usize,
    dimension: usize,
    harmonics: usize,
    coefficient: impl Fn(usize, usize) -> C64,
    mut visit: impl FnMut(usize, usize, f64),
) {
    // Preserve the original sample-major direct traversal at small orders.
    if !enabled(samples, harmonics) || transform.is_none() {
        for j in 0..samples {
            for k in 0..dimension {
                visit(
                    j,
                    k,
                    direct_sample(basis, j, harmonics, |n| coefficient(k, n)),
                );
            }
        }
        return;
    }
    for k in 0..dimension {
        if let Some(fft) = transform.as_mut()
            && fft.reconstruct(harmonics, |n| coefficient(k, n))
        {
            for j in 0..samples {
                visit(j, k, fft.data[j].re);
            }
            continue;
        }
        for j in 0..samples {
            visit(
                j,
                k,
                direct_sample(basis, j, harmonics, |n| coefficient(k, n)),
            );
        }
    }
}

fn direct_sample(
    basis: &PhaseBasis,
    sample: usize,
    harmonics: usize,
    coefficient: impl Fn(usize) -> C64,
) -> f64 {
    let mut ac = AftSum::new(0.0);
    for (n, phase) in (1..=harmonics).zip(basis.cycle(sample).skip(1)) {
        let c = coefficient(n);
        ac.add(2.0 * (c.re * phase.re - c.im * phase.im));
    }
    ac.total()
}

#[cfg(test)]
mod tests;
