//! Sample the same power-conjugate branch outputs as the time owner.
use super::{RelaxationImpedance, RelaxationImpedanceSpec, RelaxationSampling};
use crate::waveguide::WaveguideError;
use fs_material::visco::PronyRamp;
use fs_math::{c64::C64, det};

impl RelaxationImpedanceSpec {
    /// Impedance of the selected discrete load at `exp(+i*2*pi*f*n*dt)`.
    /// The output is the runtime port pressure: midpoint base R-L-C plus
    /// interval-mean exact branch pressure when that rule is selected. It is
    /// independent of initial history and characteristic port impedance.
    ///
    /// This is a sampled transfer, not the continuous material impedance or an
    /// endpoint-pressure transfer. No passivity or error grid is inferred.
    ///
    /// # Errors
    /// Invalid clock/frequency, a base compliance pole at DC, a base inertance
    /// pole at Nyquist, or unrepresentable coefficients/response.
    pub fn sampled_impedance(&self, frequency_hz: f64, dt: f64) -> Result<C64, WaveguideError> {
        let normalized = frequency_hz * dt;
        if !dt.is_finite()
            || dt <= 0.0
            || !frequency_hz.is_finite()
            || frequency_hz < 0.0
            || !normalized.is_finite()
            || normalized > 0.5
            || (frequency_hz > 0.0 && normalized == 0.0)
        {
            return Err(WaveguideError(
                "sampled relaxation needs a positive clock and finite DC-to-Nyquist frequency",
            ));
        }
        if normalized == 0.0 && self.base.compliance_m3_pa.is_some() {
            return Err(WaveguideError("sampled series compliance has a pole at DC"));
        }
        if normalized == 0.5 && self.base.inertance_pa_s2_m3 > 0.0 {
            return Err(WaveguideError(
                "sampled series inertance has a pole at Nyquist",
            ));
        }
        // Exact endpoint trigonometry preserves the static limits, including
        // arbitrarily weak RC loss at Nyquist. cos(pi/2)'s residual must not
        // become a fabricated small conductance.
        let (sin, cos) = if normalized == 0.0 {
            (0.0, 1.0)
        } else if normalized == 0.5 {
            (1.0, 0.0)
        } else {
            let angle = core::f64::consts::PI * normalized;
            (det::sin(angle), det::cos(angle))
        };
        let mut result = C64::from_re(self.base.resistance_pa_s_m3);
        if self.base.inertance_pa_s2_m3 > 0.0 {
            result.im += (2.0 * (self.base.inertance_pa_s2_m3 / dt)) * (sin / cos);
        }
        if let Some(compliance) = self.base.compliance_m3_pa {
            result.im -= (0.5 * (dt / compliance)) * (cos / sin);
        }
        for term in self.terms() {
            result = result
                + branch(
                    term.resistance_pa_s_m3,
                    term.rate_per_s * dt,
                    sin,
                    cos,
                    self.sampling,
                    false,
                )?;
        }
        for term in self.compliance_terms() {
            result = result
                + branch(
                    term.resistance_pa_s_m3,
                    term.rate_per_s * dt,
                    sin,
                    cos,
                    self.sampling,
                    true,
                )?;
        }
        if !result.re.is_finite() || !result.im.is_finite() {
            return Err(WaveguideError(
                "sampled relaxation response is not representable",
            ));
        }
        Ok(result)
    }
}

impl RelaxationImpedance {
    /// The selected load's sampled port response at its original fixed clock.
    /// This observation does not advance or change any accepted history.
    /// See [`RelaxationImpedanceSpec::sampled_impedance`] for conventions/refusals.
    pub fn sampled_impedance(&self, frequency_hz: f64) -> Result<C64, WaveguideError> {
        self.spec.sampled_impedance(frequency_hz, self.dt)
    }
}

fn branch(
    resistance: f64,
    r: f64,
    sin: f64,
    cos: f64,
    sampling: RelaxationSampling,
    compliance: bool,
) -> Result<C64, WaveguideError> {
    if !r.is_finite() || r <= 0.0 {
        return Err(WaveguideError(
            "sampled relaxation rate interval is not representable",
        ));
    }
    let (ratio, high_gain, residual) = match sampling {
        RelaxationSampling::BilinearMidpoint => (0.5 * r, resistance, 0.0),
        RelaxationSampling::ExactHeldFlow => {
            let ramp = PronyRamp::new(r);
            let (decay, mean, _) = ramp.apply(1.0, 0.0);
            let h = mean / (1.0 + decay);
            let ratio = -det::expm1(-r) / (1.0 + decay); // tanh(r/2)
            // The positive Prony heat remainder equals (1-2h)/r. Reading
            // that square-completed response avoids subtracting near-equal
            // high/low-pass transfers and preserves weak Nyquist resistance.
            let remainder = ramp.apply(-h, 1.0).2;
            let residual = if r < 1.0 {
                (resistance * r) * remainder
            } else {
                resistance * (r * remainder)
            };
            (ratio, resistance * (2.0 * h), residual)
        }
    };
    // Rotate z-a by half a sample: normalized denominator ratio*cos+i*sin.
    // Its real powers are nonnegative squares, so no passivity clamp or
    // subtractive R-Z_high formula can erase a weak physical loss.
    let x = ratio * cos;
    let scale = x.abs().max(sin.abs());
    if !scale.is_finite() || scale == 0.0 || (cos > 0.0 && x == 0.0) {
        return Err(WaveguideError(
            "sampled relaxation denominator is not representable",
        ));
    }
    let x = x / scale;
    let y = sin / scale;
    let denominator = x * x + y * y;
    let imaginary = ((high_gain * x) * y) / denominator;
    Ok(if compliance {
        C64::new(
            ((resistance * x) * x + (residual * y) * y) / denominator,
            -imaginary,
        )
    } else {
        C64::new(((high_gain * y) * y) / denominator, imaginary)
    })
}
