//! Measurement seams on the existing reed stepper; no alternate renderer.
use super::{
    AcousticRealizeError, BeatingReed, Duct, GasState, PlateBank, ReedBoreVoice, ReedSolverMode,
    Termination,
};

/// Actual filter dimensions and solver stopping rule of a constructed voice.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RealizedResolution {
    /// Reciprocal of the actual step interval, in Hz.
    pub sample_rate_hz: f64,
    /// Actual DFT realization length (including any clamp).
    pub n_fft: usize,
    /// Actual stored FIR/history length.
    pub ir_len: usize,
    /// Actual solver implementation selected for this voice.
    pub solver_mode: ReedSolverMode,
    /// Relative pressure-step tolerance; None means strict adjacent-float bisection.
    pub step_tolerance: Option<f64>,
}

impl ReedBoreVoice {
    /// Construct the existing reed voice with an explicit reflectance FFT request.
    /// Admission, mechanical initialization and priming share `with_line` with
    /// the ordinary constructor; only the existing dense-line override differs.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_fft_length(
        physics: &Duct,
        gas: &GasState,
        reed: BeatingReed,
        termination: Termination,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: Option<&fs_phs::WallPin>,
        n_fft: usize,
    ) -> Result<Self, AcousticRealizeError> {
        Self::with_line(
            physics,
            gas,
            reed,
            plates,
            listener_m,
            sample_rate_hz,
            |zc| {
                crate::driving_point::characteristic_line_realized(
                    physics,
                    gas,
                    termination,
                    sample_rate_hz,
                    line_samples,
                    zc,
                    wall,
                    Some(n_fft),
                )
                .map_err(crate::acoustic_realize::map_drive)
            },
        )
    }

    /// Read the constructed filter, not the requested FFT length.
    #[must_use]
    pub fn realized_resolution(&self) -> RealizedResolution {
        let ir_len = self.line.history().len();
        RealizedResolution {
            sample_rate_hz: 1.0 / self.dt,
            n_fft: ir_len,
            ir_len,
            solver_mode: self.solver_mode,
            step_tolerance: match self.solver_mode {
                ReedSolverMode::Strict => None,
                ReedSolverMode::FastNewton => Some(self.junction_step_tolerance),
            },
        }
    }

    /// Set an explicit fast-Newton relative pressure-step tolerance, without
    /// changing controls or the default strict mode. Massive reeds currently
    /// refuse: their coupled structural solve has a different stopping rule.
    pub fn set_junction_step_tolerance(
        &mut self,
        tolerance: f64,
    ) -> Result<(), AcousticRealizeError> {
        if self.reed.mass_kg > 0.0 || !tolerance.is_finite() || tolerance <= 0.0 || tolerance > 1.0
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "massless junction tolerance must be finite and in (0, 1]",
            });
        }
        self.junction_step_tolerance = tolerance;
        Ok(())
    }

    /// Maximum measured massless-junction volume-flow defect since construction
    /// (m³/s). No samples, or a massive reed, returns None. This is a residual
    /// diagnostic, not an error bound on the observer's sound.
    #[must_use]
    pub const fn junction_residual_max(&self) -> Option<f64> {
        self.junction_residual_max
    }

    /// Recover the actual initialized FIR from the constructor's 5 Pa priming
    /// impulse, on a clone, before rendering. The last-quarter squared-tap sum
    /// is a diagnostic of wrapped/truncated energy, not a bound on the missing
    /// infinite response or on any audio QoI. No runtime state is changed.
    pub fn initial_ir_tail_energy(&self) -> Result<(f64, f64), AcousticRealizeError> {
        if self.sample_index != 0 {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "IR inspection requires an unadvanced voice",
            });
        }
        let n = self.line.history().len();
        let mut line = self.line.clone();
        let mut total = 0.0;
        let mut tail = 0.0;
        for k in 0..n {
            let h = if k == 0 {
                line.incoming() / 5.0
            } else {
                line.push(0.0).map_err(|_| AcousticRealizeError::Reed {
                    what: "nonfinite IR inspection",
                })? / 5.0
            };
            total += h * h;
            if k >= 3 * n / 4 {
                tail += h * h;
            }
        }
        Ok((tail, total))
    }
}
