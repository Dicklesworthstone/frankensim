//! Fixed, lossless dispersion on the original characteristic line.
//!
//! Each authored section is `H(z)=(a+z^-1)/(1+a*z^-1)`, `|a|<1`.
//! Negative a delays low frequencies more than high ones, the direction of
//! stiff-string dispersion. Coefficients may differ across sections. The
//! caller owns their physical fit and valid band: this supplies no material
//! identification, continuum error bound or attenuation inferred from phase.
//!
//! Normalized scattering `[y,s_next]=[a,b;b,-a][x,s]`, `b=sqrt(1-a²)`,
//! retains one actual energy coordinate per section and direction. Its
//! quadratic identity closes wave power without integrating supplied work.
//! See J. O. Smith, Physical Audio Signal Processing, "Allpass Filters" and
//! "Dispersive 1D Wave Equation":
//! <https://www.dsprelated.com/freebooks/pasp/Dispersive_1D_Wave_Equation.html>.

use super::{PassiveWaveguide, WaveguideError, WaveguideSpec, fractional::FractionalTransit};
use fs_math::det;

/// Immutable admitted section; neither coefficients nor memories can be reset
/// during propagation. Zero a is an exact extra sample of transit, not bypass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AllpassSection {
    a: f64,
    b: f64,
}

impl AllpassSection {
    /// Admit a stable allpass coefficient without allocating. An endpoint
    /// coefficient has an unobservable/cancelled pole and is refused.
    pub fn new(a: f64) -> Result<Self, WaveguideError> {
        if !a.is_finite() || a.abs() >= 1.0 {
            return Err(WaveguideError("allpass coefficient must be finite and strictly inside (-1,1)"));
        }
        let b = det::sqrt((1.0 - a) * (1.0 + a));
        if !b.is_finite() || b <= 0.0 {
            return Err(WaveguideError("allpass scattering coefficient is not representable"));
        }
        Ok(Self { a: if a == 0.0 { 0.0 } else { a }, b })
    }

    #[must_use]
    pub const fn coefficient(&self) -> f64 { self.a }

    /// Unwrapped positive phase lag in radians under `exp(+i omega n)`.
    /// `omega` is normalized angular frequency [rad/sample], below Nyquist.
    pub fn phase_lag_rad(&self, omega: f64) -> Result<f64, WaveguideError> {
        frequency(omega)?;
        let (sin, cos) = det::sin_cos(0.5 * omega);
        Ok(2.0 * det::atan2((1.0 - self.a) * sin, (1.0 + self.a) * cos))
    }

    /// Analytic derivative of the unwrapped lag [samples]. Half-angle squares
    /// avoid cancellation at DC for a near -1 and near Nyquist for a near +1.
    pub fn group_delay_samples(&self, omega: f64) -> Result<f64, WaveguideError> {
        frequency(omega)?;
        allpass_group_delay(self.a, omega)
    }
}

pub(super) fn frequency(omega: f64) -> Result<(), WaveguideError> {
    if !omega.is_finite() || !(0.0..core::f64::consts::PI).contains(&omega) {
        return Err(WaveguideError("wave transit frequency must lie in [0, Nyquist)"));
    }
    Ok(())
}

pub(super) fn allpass_group_delay(a: f64, omega: f64) -> Result<f64, WaveguideError> {
    let (sin, cos) = det::sin_cos(0.5 * omega);
    let real = (1.0 + a) * cos;
    let imag = (1.0 - a) * sin;
    let delay = ((1.0 - a) * (1.0 + a)) / (real * real + imag * imag);
    if !delay.is_finite() || delay <= 0.0 {
        return Err(WaveguideError("allpass group delay is not representable"));
    }
    Ok(delay)
}

/// Read-only state of one bidirectional section, in original cascade order.
#[derive(Clone, Debug, PartialEq)]
pub struct DispersionStage {
    section: AllpassSection,
    state: [f64; 2],
}

impl DispersionStage {
    #[must_use]
    pub const fn section(&self) -> AllpassSection { self.section }

    /// Normalized memories [Pa], in forward/reverse order. Their energy is
    /// `dt/Z * (state[0]² + state[1]²)`; they are not loss or gain correction.
    #[must_use]
    pub const fn state_pa(&self) -> [f64; 2] { self.state }

    fn trial(&self, input: [f64; 2]) -> ([f64; 2], [f64; 2]) {
        let AllpassSection { a, b } = self.section;
        (std::array::from_fn(|i| a * input[i] + b * self.state[i]),
         std::array::from_fn(|i| b * input[i] - a * self.state[i]))
    }

    pub(super) fn commit(&mut self, input: [f64; 2]) -> [f64; 2] {
        let (waves, state) = self.trial(input);
        self.state = state;
        waves
    }
}

pub(super) fn stages(sections: &[AllpassSection]) -> Result<Vec<DispersionStage>, WaveguideError> {
    let mut result = Vec::new();
    result.try_reserve_exact(sections.len())
        .map_err(|_| WaveguideError("dispersion state allocation failed"))?;
    result.extend(sections.iter().map(|&section| DispersionStage { section, state: [0.0; 2] }));
    Ok(result)
}

/// One-way response of the actually retained integer, fractional and dispersive
/// propagation. This is a discrete model diagnostic, not a material certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransitResponse {
    /// Unwrapped lag [radians], including the original integer transit.
    pub phase_lag_rad: f64,
    /// Derivative of lag with respect to normalized frequency [samples].
    pub group_delay_samples: f64,
}

impl PassiveWaveguide {
    /// Construct an initially quiescent line with the same fixed cascade in
    /// each direction, after the optional fractional stage and before the
    /// original integer buffer. Every state participates in preview/commit and
    /// storage. An empty cascade preserves the fractional/integer arithmetic.
    ///
    /// Allocation is budgeted before construction; stepping allocates nothing
    /// and costs O(stages + log(integer delay)). Coefficients are fixed: live
    /// retuning would need an independent state-transfer and work law.
    pub fn with_dispersive_transit(
        spec: WaveguideSpec, fractional: FractionalTransit, sections: &[AllpassSection],
    ) -> Result<Self, WaveguideError> {
        if Self::required_memory_bytes(spec.one_way_samples, sections.len())? > spec.max_memory_bytes {
            return Err(WaveguideError("waveguide dispersion payload exceeds memory budget"));
        }
        let mut line = Self::with_fractional_transit(spec, fractional)?;
        line.dispersion = stages(sections)?;
        Ok(line)
    }

    #[must_use]
    pub fn dispersion_stages(&self) -> &[DispersionStage] { &self.dispersion }

    /// Original complete one-way phase and group delay at the same normalized
    /// frequency. Fractional interpolation dispersion remains separate in the
    /// retained model and is included once; no ideal fractional shift is used.
    pub fn propagation_response(&self, omega: f64) -> Result<TransitResponse, WaveguideError> {
        frequency(omega)?;
        let mut response = TransitResponse {
            phase_lag_rad: self.spec.one_way_samples as f64 * omega,
            group_delay_samples: self.spec.one_way_samples as f64,
        };
        if let Some(f) = self.fractional_transit() {
            response.phase_lag_rad += f.phase_lag_rad(omega)?;
            response.group_delay_samples += f.group_delay_samples(omega)?;
        }
        for stage in &self.dispersion {
            response.phase_lag_rad += stage.section.phase_lag_rad(omega)?;
            response.group_delay_samples += stage.section.group_delay_samples(omega)?;
        }
        if !response.phase_lag_rad.is_finite() || !response.group_delay_samples.is_finite() {
            return Err(WaveguideError("complete wave transit response overflowed"));
        }
        Ok(response)
    }

    pub(super) fn dispersion_energy_j(&self) -> f64 {
        self.dispersion.iter().fold(0.0, |energy, stage| {
            energy + (self.wave_energy(stage.state[0]) + self.wave_energy(stage.state[1]))
        })
    }

    // No candidate-state vector is allocated. The deterministic second pass in
    // commit uses exactly the same old-state arithmetic after all owners admit.
    pub(super) fn preview_dispersion(&self, mut waves: [f64; 2])
        -> Result<([f64; 2], f64), WaveguideError>
    {
        let mut energy = 0.0;
        for stage in &self.dispersion {
            let (next, state) = stage.trial(waves);
            energy += self.wave_energy(state[0]) + self.wave_energy(state[1]);
            if !next.iter().chain(&state).all(|v| v.is_finite()) || !energy.is_finite() {
                return Err(WaveguideError("dispersive wave transit or storage left the finite set"));
            }
            waves = next;
        }
        Ok((waves, energy))
    }
}
