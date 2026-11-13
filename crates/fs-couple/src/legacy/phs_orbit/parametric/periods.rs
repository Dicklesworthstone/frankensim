//! Multiple-period response charts of the unchanged original pumped source.

use super::{
    MAX_PARAMETRIC_HARMONICS, ParametricError, ParametricHamiltonian, Result, Scratch, positive,
};
use crate::phs_orbit::stability::{gcd, repeat_rms};
use core::cell::{Cell, RefCell};
use fs_math::c64::C64;
use fs_orbit::HbOrbit;

impl ParametricHamiltonian<'_> {
    /// Fundamental of the RESPONSE chart. Original input and modulation retain
    /// their authored physical frequencies, including on a longer-period view.
    pub fn angular_frequency(&self) -> f64 {
        self.modulation.omega / self.periods as f64
    }

    /// Number of original authored source periods in this response chart.
    pub fn period_multiplier(&self) -> usize {
        self.periods
    }

    /// Minimum retained order needed to represent every original pump and input
    /// harmonic on this chart. This does not bound generated response harmonics.
    pub fn required_harmonics(&self) -> usize {
        self.bandwidth()
    }

    /// Borrow the SAME original source on a chart spanning `periods` of this
    /// view's cycles. Nested views multiply their total original-period count.
    ///
    /// Only the response clock and required numerical bandwidth change. Q,
    /// Qdot, drive amplitudes/phases, H0, J/R/G, every material-memory state and
    /// the exact native/authored H0 tangent remain bound to their original
    /// owners. No Fourier coefficients, matrices or physical frequencies are
    /// rebuilt or retuned: evaluating either view at the same physical time and
    /// state gives the same original physical field, tangent and pump work.
    ///
    /// This makes period-doubled and higher pumped responses available to the
    /// existing seeded/adaptive HB, full-cycle variational shooting, sideband,
    /// sensitivity and power/work APIs. Seed all original-state harmonics to
    /// select a branch. A longer chart also admits shorter-period solutions;
    /// it neither finds a subharmonic branch automatically nor proves its
    /// primitive period. Use `require_subharmonic_separation` independently.
    ///
    /// The total multiplier and remapped source bandwidth are bounded by 256.
    /// Unrepresentable clock/period and unsupported bandwidth refuse before
    /// scratch allocation or any original source/tangent evaluation. The view
    /// owns only ordinary source scratch; source Fourier data and tangent are
    /// borrowed unchanged. Existing per-solver resource limits still apply.
    pub fn over_periods(&self, periods: usize) -> Result<ParametricHamiltonian<'_>> {
        let total = self
            .periods
            .checked_mul(periods)
            .filter(|n| (1..=MAX_PARAMETRIC_HARMONICS).contains(n))
            .ok_or(ParametricError::Invalid(
                "response period multiplier must remain in 1..=256",
            ))?;
        let source_bandwidth = self
            .modulation
            .harmonics
            .max(self.input.highest_nonzero_harmonic());
        if source_bandwidth
            .checked_mul(total)
            .is_none_or(|n| n > MAX_PARAMETRIC_HARMONICS)
        {
            return Err(ParametricError::Invalid(
                "multiple-period pump/input bandwidth exceeds the harmonic limit",
            ));
        }
        let omega = self.modulation.omega / total as f64;
        if !positive(omega) || !positive(core::f64::consts::TAU / omega) {
            return Err(ParametricError::Invalid(
                "multiple-period response clock is not representable",
            ));
        }
        let (n, m) = (self.system.state_dim(), self.system.port_dim());
        Ok(ParametricHamiltonian {
            system: self.system,
            input: self.input,
            modulation: self.modulation,
            periods: total,
            tangent: self.tangent.borrowed(n),
            dissipation: self.dissipation,
            // These phases belong to the ORIGINAL physical source, not a
            // remapped Fourier array. Only the numerical chart has longer period.
            scratch: RefCell::new({
                let mut scratch = Scratch::new(n, m, self.modulation.harmonics);
                if self.dissipation.is_some() {
                    scratch.dissipation = Some(super::dissipation::Scratch::new(n));
                }
                scratch
            }),
            source_refused: Cell::new(false),
            failure: RefCell::new(None),
        })
    }

    /// Require the complete retained state to distinguish this response from
    /// every shorter source-synchronized repeat whose period divides the chart.
    /// The response chart must span at least two authored source periods.
    /// Each original state has its own positive physical RMS tolerance; DC and
    /// preload do not inflate it. Explicit memory can resolve a longer cycle
    /// even when a displacement observation cannot.
    ///
    /// The primitive authored source clock is determined jointly from nonzero
    /// PUMP AND INPUT harmonics. A pump-only source with a DC/zero port drive is
    /// therefore admitted. Constant sources or artificially padded source
    /// clocks refuse. The existing exact-integer Parseval shift owner performs
    /// the waveform comparison without sampling, phase alignment or small-tail
    /// thresholds hidden in an aggregate norm.
    ///
    /// This is retained-waveform separation, not an existence, stability,
    /// omitted-harmonic or general mathematical primitive-period certificate.
    /// It does not replace the original HB residual or time-domain checks.
    pub fn require_subharmonic_separation(
        &self,
        candidate: &HbOrbit,
        state_rms_tolerances: &[f64],
    ) -> Result<()> {
        let n = self.system.state_dim();
        let width = candidate.coeffs.first().map_or(0, Vec::len);
        if self.periods < 2
            || candidate.omega.to_bits() != self.angular_frequency().to_bits()
            || candidate.coeffs.len() != n
            || width < 2
            || candidate.coeffs.iter().any(|row| {
                row.len() != width
                    || row[0].im != 0.0
                    || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
            || state_rms_tolerances.len() != n
            || state_rms_tolerances.iter().any(|v| !positive(*v))
        {
            return Err(ParametricError::Invalid(
                "subharmonic separation needs the complete original waveform, response clock and positive state RMS tolerances",
            ));
        }
        self.admit_order(width - 1)?;
        let mut common = 0;
        for coefficients in self.input.coefficients.iter().chain(
            self.modulation
                .entries
                .iter()
                .map(|entry| &entry.coefficients),
        ) {
            for (harmonic, value) in coefficients.iter().enumerate().skip(1) {
                if *value != C64::ZERO {
                    common = gcd(common, harmonic);
                }
            }
        }
        if common != 1 {
            return Err(ParametricError::Invalid(
                "subharmonic order requires a nonconstant primitive authored pump/input clock",
            ));
        }
        for repeat in 1..self.periods {
            if self.periods % repeat != 0 {
                continue;
            }
            let mut separated = false;
            for (row, &tolerance) in candidate.coeffs.iter().zip(state_rms_tolerances) {
                let difference = repeat_rms(row, self.periods, repeat).map_err(|_| {
                    ParametricError::Invalid("subharmonic waveform separation is not representable")
                })?;
                if difference > tolerance {
                    separated = true;
                    break;
                }
            }
            if !separated {
                return Err(ParametricError::Invalid(
                    "candidate is unresolved from a shorter pump/input-synchronized cycle",
                ));
            }
        }
        Ok(())
    }
}
