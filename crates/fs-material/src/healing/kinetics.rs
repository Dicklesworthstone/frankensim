//! Reduced-time welding kinetics on an explicitly supplied linear thermal path.
//!
//! `r(T) = r_ref / a_T(T)` uses the SAME immutable [`TemperatureShift`] as the
//! viscoelastic owner: heating accelerates the molecular clock, not its inverse.
//! The accumulated dimensionless dose is `r_ref * integral dt/a_T`. It enters
//! the existing contact history and separately supplied strength map unchanged.
//! See Seppala et al., Soft Matter 13 (2017), 6761–6769, equation 5:
//! <https://pmc.ncbi.nlm.nih.gov/articles/PMC5684701/>.
//!
//! An optional, explicit cutoff sets mobility to zero at/below that temperature.
//! It neither supplies a glass-transition temperature nor extends the calibrated
//! shift domain. The cutoff is split exactly before quadrature. Adaptive Simpson
//! integration has declared dose tolerances and hard work/stack bounds; its
//! error is a numerical estimate, not a certified bound or material uncertainty.
//! No fitted coefficients, intimate-contact area, crystallization, degradation,
//! strength, or simultaneous healing of mechanically damaged bonds is inferred.

use super::{HealingError, HealingLaw, nonnegative, positive, validate_identity};
use crate::visco::temperature::TemperatureShift;

/// Accuracy and bounded work for ONE actual contact-temperature ramp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MobilityIntegration {
    /// Absolute estimated quadrature error in dimensionless dose, strictly positive.
    pub absolute_dose_tolerance: f64,
    /// Relative estimated quadrature error, in `[0, 1)`.
    pub relative_dose_tolerance: f64,
    /// Maximum shift evaluations, in `5..=16385`; endpoints count toward it.
    pub maximum_evaluations: usize,
    /// Maximum bisection depth, in `1..=32`; no heap is allocated during integration.
    pub maximum_depth: usize,
}

impl MobilityIntegration {
    fn validate(self) -> Result<(), HealingError> {
        positive(self.absolute_dose_tolerance)?;
        if !self.relative_dose_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_dose_tolerance)
            || !(5..=16_385).contains(&self.maximum_evaluations)
            || !(1..=32).contains(&self.maximum_depth)
        {
            return Err(HealingError::Invalid(
                "invalid welding integration controls",
            ));
        }
        Ok(())
    }
}

/// Numerical observations for the accepted dimensionless exposure increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MobilityDose {
    pub dose: f64,
    /// Accumulated Simpson error estimate in dose units; excludes model error.
    pub estimated_absolute_error: f64,
    /// Actual canonical shift evaluations in this operation.
    pub evaluations: usize,
}

/// One declared welding clock and its original WLF or Arrhenius calibration.
#[derive(Clone, Debug, PartialEq)]
pub struct ShiftedMobility {
    identity: String,
    shift: TemperatureShift,
    reference_rate_per_s: f64,
    cutoff_temperature_k: Option<f64>,
    integration: MobilityIntegration,
}

impl ShiftedMobility {
    /// Preserve the caller's canonical shift, rate and law/version identity.
    /// The optional cutoff must lie inside its calibration band and below its
    /// reference temperature. No default material or kinetic coefficients exist.
    pub fn new(
        identity: String,
        shift: TemperatureShift,
        reference_rate_per_s: f64,
        cutoff_temperature_k: Option<f64>,
        integration: MobilityIntegration,
    ) -> Result<Self, HealingError> {
        validate_identity(&identity)?;
        positive(reference_rate_per_s)?;
        integration.validate()?;
        if let Some(cutoff) = cutoff_temperature_k {
            if !cutoff.is_finite()
                || cutoff < shift.validity_k()[0]
                || cutoff >= shift.reference_temperature_k()
            {
                return Err(HealingError::Invalid(
                    "welding cutoff must be in the shift band below its reference temperature",
                ));
            }
        }
        let law = Self {
            identity,
            shift,
            reference_rate_per_s,
            cutoff_temperature_k,
            integration,
        };
        // The active rate increases with T. Admit both extrema, including the
        // active-side limit at a cutoff, rather than hiding underflow behind zero.
        law.active_rate_per_s(cutoff_temperature_k.unwrap_or(shift.validity_k()[0]))?;
        law.active_rate_per_s(shift.validity_k()[1])?;
        Ok(law)
    }

    pub const fn shift(&self) -> &TemperatureShift {
        &self.shift
    }
    pub const fn reference_rate_per_s(&self) -> f64 {
        self.reference_rate_per_s
    }
    pub const fn cutoff_temperature_k(&self) -> Option<f64> {
        self.cutoff_temperature_k
    }
    pub const fn integration(&self) -> MobilityIntegration {
        self.integration
    }

    /// The same calibrated molecular rate used in exposure integration, in s^-1.
    /// Even zero-mobility temperatures must remain inside the original domain.
    pub fn rate_per_s(&self, temperature_k: f64) -> Result<f64, HealingError> {
        self.require_temperature(temperature_k)?;
        if self
            .cutoff_temperature_k
            .is_some_and(|cutoff| temperature_k <= cutoff)
        {
            Ok(0.0)
        } else {
            self.active_rate_per_s(temperature_k)
        }
    }

    /// Integrate one actual straight thermal path without mutating a history.
    /// Isothermal exposure is exact rate-times-time arithmetic. Reversing a
    /// linear ramp gives the identical arithmetic order and dose. Work/depth
    /// exhaustion refuses the complete increment instead of returning an
    /// under-resolved dose that could silently create a stronger bond.
    pub fn integrate(
        &self,
        start_k: f64,
        end_k: f64,
        duration_s: f64,
    ) -> Result<MobilityDose, HealingError> {
        nonnegative(duration_s)?;
        self.require_temperature(start_k)?;
        self.require_temperature(end_k)?;
        let zero = MobilityDose {
            dose: 0.0,
            estimated_absolute_error: 0.0,
            evaluations: 0,
        };
        if duration_s == 0.0 {
            return Ok(zero);
        }
        let low = start_k.min(end_k);
        let high = start_k.max(end_k);
        if self
            .cutoff_temperature_k
            .is_some_and(|cutoff| high <= cutoff)
        {
            return Ok(zero);
        }
        if low == high {
            return Ok(MobilityDose {
                dose: nonnegative(self.active_rate_per_s(low)? * duration_s)?,
                estimated_absolute_error: 0.0,
                evaluations: 1,
            });
        }
        let active_low = low.max(self.cutoff_temperature_k.unwrap_or(low));
        let active_duration = duration_s * ((high - active_low) / (high - low));
        if active_duration == 0.0 {
            return Ok(zero);
        }
        let high_rate = self.active_rate_per_s(high)?;
        let first = self.active_rate_per_s(active_low)? / high_rate;
        let mut ramp = Ramp {
            law: self,
            low: active_low,
            high,
            high_rate,
            evaluations: 2,
        };
        let middle = ramp.value(0.5)?;
        let coarse = simpson(1.0, first, middle, 1.0);
        // Scaling to the largest active rate keeps the integrand in [0, 1].
        // An infinite normalized absolute tolerance is harmless: its physical
        // tolerance exceeds any possible integral on this normalized interval.
        let absolute = (self.integration.absolute_dose_tolerance / high_rate) / active_duration;
        let (integral, error) = ramp.refine(
            Panel {
                low: 0.0,
                high: 1.0,
                values: [first, middle, 1.0],
                estimate: coarse,
            },
            absolute,
            0,
        )?;
        Ok(MobilityDose {
            dose: scaled_dose(integral, high_rate, active_duration)?,
            estimated_absolute_error: scaled_dose(error, high_rate, active_duration)?,
            evaluations: ramp.evaluations,
        })
    }

    fn require_temperature(&self, temperature_k: f64) -> Result<(), HealingError> {
        if !temperature_k.is_finite()
            || !(self.shift.validity_k()[0]..=self.shift.validity_k()[1]).contains(&temperature_k)
        {
            Err(HealingError::Temperature(temperature_k))
        } else {
            Ok(())
        }
    }

    fn active_rate_per_s(&self, temperature_k: f64) -> Result<f64, HealingError> {
        let state = self
            .shift
            .evaluate(temperature_k)
            .map_err(|_| HealingError::Invalid("original welding temperature shift refused"))?;
        positive(self.reference_rate_per_s / state.factor())
    }
}

impl HealingLaw for ShiftedMobility {
    fn identity(&self) -> &str {
        &self.identity
    }
    fn dose_increment(
        &self,
        start_k: f64,
        end_k: f64,
        duration_s: f64,
    ) -> Result<f64, HealingError> {
        self.integrate(start_k, end_k, duration_s)
            .map(|result| result.dose)
    }
}

#[derive(Clone, Copy)]
struct Panel {
    low: f64,
    high: f64,
    values: [f64; 3],
    estimate: f64,
}
struct Ramp<'a> {
    law: &'a ShiftedMobility,
    low: f64,
    high: f64,
    high_rate: f64,
    evaluations: usize,
}
impl Ramp<'_> {
    fn value(&mut self, coordinate: f64) -> Result<f64, HealingError> {
        if self.evaluations >= self.law.integration.maximum_evaluations {
            return Err(HealingError::IntegrationBudget);
        }
        self.evaluations += 1;
        let temperature = self.low + (self.high - self.low) * coordinate;
        Ok(self.law.active_rate_per_s(temperature)? / self.high_rate)
    }

    fn refine(
        &mut self,
        panel: Panel,
        absolute: f64,
        depth: usize,
    ) -> Result<(f64, f64), HealingError> {
        let width = panel.high - panel.low;
        let middle = panel.low + 0.5 * width;
        let left_mid = self.value(panel.low + 0.25 * width)?;
        let right_mid = self.value(panel.low + 0.75 * width)?;
        let [left, centre, right] = panel.values;
        let first = simpson(0.5 * width, left, left_mid, centre);
        let second = simpson(0.5 * width, centre, right_mid, right);
        let fine = first + second;
        let correction = (fine - panel.estimate) / 15.0;
        let error = correction.abs();
        let estimate = nonnegative(fine + correction)?;
        let tolerance = absolute * width + self.law.integration.relative_dose_tolerance * estimate;
        if error <= tolerance {
            return Ok((estimate, error));
        }
        if depth >= self.law.integration.maximum_depth {
            return Err(HealingError::IntegrationBudget);
        }
        let a = self.refine(
            Panel {
                low: panel.low,
                high: middle,
                values: [left, left_mid, centre],
                estimate: first,
            },
            absolute,
            depth + 1,
        )?;
        let b = self.refine(
            Panel {
                low: middle,
                high: panel.high,
                values: [centre, right_mid, right],
                estimate: second,
            },
            absolute,
            depth + 1,
        )?;
        Ok((a.0 + b.0, a.1 + b.1))
    }
}

fn simpson(width: f64, left: f64, middle: f64, right: f64) -> f64 {
    width * ((left + right) / 6.0 + (2.0 / 3.0) * middle)
}

fn scaled_dose(integral: f64, rate: f64, duration: f64) -> Result<f64, HealingError> {
    let first = integral * rate;
    if first == 0.0 && integral != 0.0 {
        nonnegative((integral * duration) * rate)
    } else {
        nonnegative(first * duration)
    }
}

#[cfg(test)]
mod tests;
