//! Fixed, lossless first-order allpass transit with explicitly retained storage.
//!
//! H(z)=(a+z^-1)/(1+a z^-1), a=(1-d)/(1+d), has unit magnitude and
//! group delay d at zero frequency. It is NOT an exact broadband time shift:
//! its phase bends near Nyquist, and a small d has a long alternating tail.
//! See J. O. Smith, Physical Audio Signal Processing, First-Order Allpass
//! Interpolation: https://www.dsprelated.com/freebooks/pasp/First_Order_Allpass_Interpolation.html
//!
//! The normalized scattering realization is [y,s_next]=[a,b;b,-a][x,s],
//! b=sqrt(1-a^2). Thus x^2+s^2=y^2+s_next^2 in exact arithmetic. Unlike a
//! direct-form history or an integrated port-work estimate, s^2 is an actual
//! nonnegative storage coordinate. There is no gain correction or hidden loss.

use super::{PassiveWaveguide, WaveguideError, WaveguideSpec};

/// Immutable sub-sample transit. Zero selects the original exact integer path.
/// A nonzero delay is a first-order allpass approximation, not a retuning law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FractionalTransit {
    requested: f64,
    a: f64,
    b: f64,
}
impl FractionalTransit {
    /// Admit a residual delay in [0,1) samples without allocating or stepping.
    /// Coefficients that round onto the unit circle refuse instead of silently
    /// replacing a small physical transit with zero.
    pub fn new(samples: f64) -> Result<Self, WaveguideError> {
        if !samples.is_finite() || !(0.0..1.0).contains(&samples) {
            return Err(WaveguideError("fractional transit must be finite and in [0,1) samples"));
        }
        if samples == 0.0 { return Ok(Self { requested: 0.0, a: 1.0, b: 0.0 }); }
        let a = (1.0-samples)/(1.0+samples);
        let b = fs_math::det::sqrt((1.0-a)*(1.0+a));
        if !(0.0..1.0).contains(&a) || !b.is_finite() || b <= 0.0 {
            return Err(WaveguideError("fractional transit coefficients are not representable"));
        }
        Ok(Self { requested: samples, a, b })
    }
    /// Original requested delay at zero frequency [samples].
    #[must_use]
    pub const fn samples(&self) -> f64 { self.requested }
    /// Actual zero-frequency group delay of the rounded coefficient [samples].
    #[must_use]
    pub fn realized_dc_samples(&self) -> f64 {
        if self.requested == 0.0 { 0.0 } else { (1.0-self.a)/(1.0+self.a) }
    }
    /// Maximum absolute one-way phase discrepancy from the requested ideal
    /// fractional delay on [0,omega_max], in radians per sample.
    ///
    /// The allpass phase error relative to its own DC delay is monotone here.
    /// Add the coefficient-rounding DC error by the triangle inequality. This
    /// is an analytic-model calculation in binary64, not an interval enclosure
    /// or a bound on omitted physics, aliased inputs or a whole network path.
    pub fn phase_error_bound_rad(&self, omega_max: f64) -> Result<f64, WaveguideError> {
        if !omega_max.is_finite() || !(0.0..core::f64::consts::PI).contains(&omega_max) {
            return Err(WaveguideError("fractional transit band must lie below Nyquist"));
        }
        if self.requested == 0.0 || omega_max == 0.0 { return Ok(0.0); }
        let d = self.realized_dc_samples();
        let phase = self.phase_lag_rad(omega_max)?;
        let error = (phase-d*omega_max).abs() + (d-self.requested).abs()*omega_max;
        if !error.is_finite() { return Err(WaveguideError("fractional transit phase calculation overflowed")); }
        Ok(error)
    }
    /// Positive unwrapped phase lag of the ACTUAL rounded allpass coefficient.
    /// Under exp(+i omega n), the response is exp(-i * phase_lag_rad(omega)).
    /// This is the same fixed filter used by propagation, not the ideal d*omega.
    /// It lets physical response checks include interpolation dispersion without
    /// copying the filter law or claiming broadband-exact travel.
    ///
    /// # Errors
    /// Nonfinite frequency or a frequency outside [0, Nyquist).
    pub fn phase_lag_rad(&self, omega: f64) -> Result<f64, WaveguideError> {
        if !omega.is_finite() || !(0.0..core::f64::consts::PI).contains(&omega) {
            return Err(WaveguideError("fractional transit band must lie below Nyquist"));
        }
        if self.requested == 0.0 || omega == 0.0 { return Ok(0.0); }
        let angle = 0.5*omega;
        let tangent = fs_math::det::sin(angle)/fs_math::det::cos(angle);
        let phase = 2.0*(self.realized_dc_samples()*tangent).atan();
        if !phase.is_finite() {
            return Err(WaveguideError("fractional transit phase calculation overflowed"));
        }
        Ok(phase)
    }

}

#[derive(Debug, Clone)]
pub(super) struct FractionalPair {
    pub(super) transit: FractionalTransit,
    pub(super) state: [f64; 2],
}
impl FractionalPair {
    pub(super) fn new(transit: FractionalTransit) -> Option<Self> {
        if transit.samples() == 0.0 { None } else { Some(Self { transit, state: [0.0; 2] }) }
    }
    // Pure trial; callers validate waves AND storage before publishing any state.
    pub(super) fn trial(&self, input: [f64; 2]) -> ([f64; 2], [f64; 2]) {
        let a=self.transit.a; let b=self.transit.b;
        (std::array::from_fn(|i| a*input[i]+b*self.state[i]),
         std::array::from_fn(|i| b*input[i]-a*self.state[i]))
    }
}

impl PassiveWaveguide {
    /// Create a line with N exact integer samples followed by a fixed fractional
    /// allpass in EACH direction. They commute because coefficients are fixed.
    /// Internal normalized states store energy and participate in all previews.
    /// Zero uses precisely the original arithmetic and has no allpass state.
    ///
    /// No live delay change or reset is exposed: moving boundaries need an
    /// independent physical work/state-transfer law. The original positive
    /// integer transit is retained, so there is no algebraic junction loop.
    pub fn with_fractional_transit(spec: WaveguideSpec, transit: FractionalTransit)
        -> Result<Self, WaveguideError>
    {
        let mut line=Self::new(spec)?;
        line.fractional=FractionalPair::new(transit);
        Ok(line)
    }
    /// Actual optional allpass specification; absent means exact integer travel.
    #[must_use]
    pub fn fractional_transit(&self) -> Option<FractionalTransit> {
        self.fractional.as_ref().map(|f|f.transit)
    }
    /// Normalized allpass memory in pressure units, one per propagation direction.
    /// The energy is dt/Z times the sum of their squares. Read-only observation
    /// cannot change a wave state behind a coupled solver's accepted sample.
    #[must_use]
    pub fn fractional_state_pa(&self) -> [f64; 2] {
        self.fractional.as_ref().map_or([0.0;2],|f|f.state)
    }
    pub(super) fn preview_pair_storage(&self, outgoing: f64, reflected: f64)
        -> Result<f64, WaveguideError>
    {
        let input=[outgoing,reflected];
        let (waves,state)=self.fractional.as_ref().map_or((input,[0.0;2]),|f|f.trial(input));
        let mut stored=self.replaced_root(0,self.wave_energy(waves[0]))
            +self.replaced_root(2*self.leaves,self.wave_energy(waves[1]));
        if self.fractional.is_some() {
            stored+=self.wave_energy(state[0])+self.wave_energy(state[1]);
        }
        if !input.iter().chain(&waves).chain(&state).all(|x|x.is_finite()) || !stored.is_finite() {
            return Err(WaveguideError("fractional wave transit or storage left the finite set"));
        }
        Ok(stored)
    }
}
