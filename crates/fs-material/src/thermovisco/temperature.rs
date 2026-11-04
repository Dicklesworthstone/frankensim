//! Exact WLF/Arrhenius kinetics in the existing thermoviscoelastic history owner.
//!
//! The declared shift multiplies the original constant reference time scale.
//! It never changes spring energies, free thermal strain, branch histories or
//! the selected history integrator. In particular, changing a relaxation time
//! is not reversible thermal work and is not an independent source of heat.
//! Existing caloric feedback deposits only the original history update's
//! physical dashpot dissipation and retains its reversible-work accounting.

use super::{Error, ThermalKnot, ThermoMaxwell, finite, positive};
use crate::visco::temperature::TemperatureShift;

impl ThermoMaxwell {
    /// Attach the SAME calibrated shift used by scalar and tensor viscoelastic
    /// laws: `tau_j(T) = tau_j * reference_time_scale * a_T(T)`.
    ///
    /// The original time-scale table must be constant, so a tabulated thermal
    /// shift cannot accidentally be applied twice. Its explicit non-unit value
    /// is retained at the calibration reference temperature. Every modulus and
    /// free-strain knot, branch coefficient, limit and integration choice stays
    /// with this law. The calibration must cover the whole mechanical interval,
    /// and its reference temperature must lie inside that interval.
    ///
    /// Mechanical tangents, accepted history proposals and caloric fixed-point
    /// trials evaluate the analytic factor at the actual proposed temperature.
    /// Both history integrators retain their endpoint-temperature convention:
    /// exact strain ramps are not exact for continuously varying temperature.
    /// Refine thermal time steps to resolve that variation.
    ///
    /// # Errors
    /// Refuses a second shift, varying reference time scaling, incompatible
    /// thermal support, or unrepresentable scaled relaxation times anywhere in
    /// the admitted temperature interval. Nothing resets a material-point state.
    pub fn with_temperature_shift(mut self, shift: TemperatureShift) -> Result<Self, Error> {
        if self.temperature_shift.is_some() {
            return Err(Error::Invalid("a temperature shift is already attached"));
        }
        let reference_scale = self.knots[0].relaxation_scale;
        if self
            .knots
            .iter()
            .any(|k| k.relaxation_scale != reference_scale)
        {
            return Err(Error::Invalid(
                "analytic temperature kinetics require a constant original relaxation multiplier",
            ));
        }
        let low = self.knots[0].temperature_k;
        let high = self.knots[self.knots.len() - 1].temperature_k;
        let [valid_low, valid_high] = shift.validity_k();
        if valid_low > low
            || valid_high < high
            || !(low..=high).contains(&shift.reference_temperature_k())
        {
            return Err(Error::Invalid(
                "shift support must cover the mechanical interval and its reference must be inside",
            ));
        }
        // Both admitted shift laws are monotone. A constant positive reference
        // multiplier therefore attains every scaled-time extremum at an endpoint.
        for temperature_k in [low, high] {
            let state = shift.evaluate(temperature_k).map_err(|_| {
                Error::Invalid("temperature-shift calibration is not representable")
            })?;
            let scale = finite(reference_scale * state.factor())?;
            positive(scale)?;
            for &(_, tau) in &self.branches {
                positive(finite(tau * scale)?)?;
            }
        }
        self.temperature_shift = Some(shift);
        Ok(self)
    }

    /// Original analytic calibration, without synthesizing a tabulated curve.
    #[must_use]
    pub const fn temperature_shift(&self) -> Option<&TemperatureShift> {
        self.temperature_shift.as_ref()
    }

    // Only kinetics consume the exponential shift. Storage, thermal eigenstrain
    // and stress evaluation use the original table without redundant exponentials.
    pub(super) fn kinetics_at(&self, temperature_k: f64) -> Result<ThermalKnot, Error> {
        let mut thermal = self.at(temperature_k)?;
        if let Some(shift) = self.temperature_shift {
            let state = shift.evaluate(temperature_k).map_err(|_| {
                Error::Invalid("temperature-shift calibration is not representable")
            })?;
            // This constant was checked at attachment; avoid interpolating equal
            // reference scales and introducing a temperature-dependent roundoff.
            thermal.relaxation_scale = finite(self.knots[0].relaxation_scale * state.factor())?;
            positive(thermal.relaxation_scale)?;
        }
        Ok(thermal)
    }
}
