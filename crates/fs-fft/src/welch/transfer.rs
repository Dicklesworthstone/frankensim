//! Frequency-dependent coherent channel transfers of a joint Welch estimate.
//!
//! Operate on retained factors, B_out(f)=H(f) B_in(f), before any power
//! contraction. This preserves S_out=H S_in H* and the original record clock,
//! rank, bin widths and one-sided convention without constructing a dense CSD.
//! It supports calibrated channels, linear propagation and delay-and-sum or
//! differential readouts. The caller supplies the physical transfer and units.
//! A spectral-bin transfer is not a causal finite-record time-domain filter:
//! taper leakage and edge transients are not reconstructed or deconvolved.

use super::{C64, Dd, WelchCrossSpectrum, WelchError, allocate, poll};

/// Admission for a complete frequency-dependent channel transfer.
///
/// Work counts complex H*B products over every bin, segment and output/input
/// pair. Numerical storage includes all new factors and one reused frequency
/// matrix. The borrowed input, callback-owned data and allocator metadata are
/// excluded. No FFT, channel covariance matrix or eigendecomposition is used.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WelchChannelTransferBudget {
    pub maximum_output_channels: usize,
    pub maximum_complex_products: usize,
    pub maximum_workspace_bytes: usize,
}

/// Exact logical payload and algebraic work of a channel transfer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WelchChannelTransferRequirements {
    /// Total coefficient callbacks; each H(f)[output,input] is requested once.
    pub coefficient_values: usize,
    pub complex_products: usize,
    pub output_factors: usize,
    pub workspace_bytes: usize,
}

impl WelchCrossSpectrum {
    /// Inspect bounds without allocating, reading coefficients or changing data.
    pub fn channel_transfer_requirements(
        &self,
        output_channels: usize,
    ) -> Result<WelchChannelTransferRequirements, WelchError> {
        if output_channels == 0 {
            return Err(WelchError::Invalid("channel transfer needs an output"));
        }
        let overflow = || WelchError::BudgetExceeded("channel transfer shape overflow");
        let matrix_values = output_channels
            .checked_mul(self.channels)
            .ok_or_else(overflow)?;
        let coefficient_values = self
            .bin_count()
            .checked_mul(matrix_values)
            .ok_or_else(overflow)?;
        let complex_products = coefficient_values
            .checked_mul(self.segments)
            .ok_or_else(overflow)?;
        let output_factors = self
            .bin_count()
            .checked_mul(self.segments)
            .and_then(|n| n.checked_mul(output_channels))
            .ok_or_else(overflow)?;
        let workspace_bytes = output_factors
            .checked_add(matrix_values)
            .and_then(|n| n.checked_mul(core::mem::size_of::<C64>()))
            .ok_or_else(overflow)?;
        Ok(WelchChannelTransferRequirements {
            coefficient_values,
            complex_products,
            output_factors,
            workspace_bytes,
        })
    }

    /// Apply the physical complex map y(f)=H(f)x(f) to EVERY retained factor.
    ///
    /// coefficient(bin, frequency_hz, output, input) supplies H in the original
    /// exp(-i omega t) convention. There is no implicit conjugation. For a
    /// beamformer written y=w* x, supply conjugated weights. A time delay tau
    /// has H=exp(+i omega tau) in this convention. Identity/order, calibration,
    /// dimensions and the validity of the physical transfer belong to callers.
    ///
    /// Each coefficient is requested once and cached across all segments of
    /// its frequency. Endpoint gains must be real at DC and EVEN-N Nyquist,
    /// preserving the conjugate symmetry of real recorded channels. Odd-N's
    /// highest positive bin remains complex. Output peak-bin views retain the
    /// same frequency and normalization, ready for the original physical
    /// propagation owners without another sqrt(df) or peak/RMS conversion.
    ///
    /// Complex products and channel sums use double-double arithmetic before
    /// each output factor is rounded once. This protects cancellation but is
    /// not an arbitrary-null or overflow certificate. Unrepresentable factors,
    /// invalid coefficients, callback errors and cancellation refuse the whole
    /// result; the input is immutable and no partial output is published.
    ///
    /// Budget admission precedes callbacks and allocation. Cancellation is
    /// checked at most every 1024 coefficient evaluations or complex products,
    /// as well as before publication. Callback-internal work remains its owner.
    pub fn transfer_channels(
        &self,
        output_channels: usize,
        mut coefficient: impl FnMut(usize, f64, usize, usize) -> Result<C64, WelchError>,
        budget: WelchChannelTransferBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WelchError> {
        poll(&mut cancelled)?;
        let needs = self.channel_transfer_requirements(output_channels)?;
        if output_channels > budget.maximum_output_channels {
            return Err(WelchError::BudgetExceeded("transferred output channels"));
        }
        if needs.complex_products > budget.maximum_complex_products {
            return Err(WelchError::BudgetExceeded("channel transfer products"));
        }
        if needs.workspace_bytes > budget.maximum_workspace_bytes {
            return Err(WelchError::BudgetExceeded("channel transfer workspace"));
        }
        let matrix_values = output_channels * self.channels;
        let mut matrix = allocate(matrix_values, C64::default())?;
        let mut factors = allocate(needs.output_factors, C64::default())?;
        let mut work_since_poll = 0usize;
        for bin in 0..self.bin_count() {
            poll(&mut cancelled)?;
            let frequency_hz = self.frequency_hz(bin)?;
            let real_endpoint = bin == 0
                || (self.config.segment_samples % 2 == 0 && bin == self.config.segment_samples / 2);
            for output in 0..output_channels {
                for input in 0..self.channels {
                    if work_since_poll == 1024 {
                        poll(&mut cancelled)?;
                        work_since_poll = 0;
                    }
                    let value = coefficient(bin, frequency_hz, output, input)?;
                    if !value.re.is_finite() || !value.im.is_finite() {
                        return Err(WelchError::NonFinite);
                    }
                    if real_endpoint && value.im != 0.0 {
                        return Err(WelchError::Invalid(
                            "real DC/Nyquist channel transfers require real coefficients",
                        ));
                    }
                    matrix[output * self.channels + input] = value;
                    work_since_poll += 1;
                }
            }
            for segment in 0..self.segments {
                let input_start = (bin * self.segments + segment) * self.channels;
                let input = &self.factors[input_start..input_start + self.channels];
                for output in 0..output_channels {
                    let mut re = Dd::ZERO;
                    let mut im = Dd::ZERO;
                    let row = &matrix[output * self.channels..(output + 1) * self.channels];
                    for (gain, sample) in row.iter().zip(input) {
                        if work_since_poll == 1024 {
                            poll(&mut cancelled)?;
                            work_since_poll = 0;
                        }
                        let hr = Dd::from_f64(gain.re);
                        let hi = Dd::from_f64(gain.im);
                        let xr = Dd::from_f64(sample.re);
                        let xi = Dd::from_f64(sample.im);
                        re = re + hr * xr - hi * xi;
                        im = im + hr * xi + hi * xr;
                        work_since_poll += 1;
                    }
                    let value = C64::new(re.to_f64(), im.to_f64());
                    if !re.hi.is_finite()
                        || !re.lo.is_finite()
                        || !im.hi.is_finite()
                        || !im.lo.is_finite()
                        || !value.re.is_finite()
                        || !value.im.is_finite()
                    {
                        return Err(WelchError::Unrepresentable);
                    }
                    factors[(bin * self.segments + segment) * output_channels + output] = value;
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            config: self.config,
            timestep_s: self.timestep_s,
            bin_width_hz: self.bin_width_hz,
            channels: output_channels,
            segments: self.segments,
            input_samples: self.input_samples,
            used_samples: self.used_samples,
            factors,
        })
    }
}

#[cfg(test)]
#[path = "transfer/tests.rs"]
mod tests;
