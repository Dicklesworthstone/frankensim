//! Contact-time/temperature exposure for polymer welding (issue #1, R3).
//!
//! The state is the dimensionless mobility dose D = integral r(T) dt, not a
//! strength or thermal conductance. A caller-supplied mobility table integrates
//! exactly on its piecewise-linear rate curve and a linear temperature path.
//! Cooling pauses development where the declared rate is zero; reheating keeps
//! all earlier exposure. No material-specific mobility or strength is invented.
//!
//! Equivalent-time motivation: Seppala et al., Soft Matter 13 (2017), 6761-6769,
//! doi:10.1039/C7SM00950J. A scalar dose does not model chain orientation,
//! disentanglement, crystallization, degradation, pressure or intimate-contact
//! creation. Its validity and relation to strength need independent material data.

use core::fmt;

/// Explicit, separately supplied mapping into a frozen cohesive law.
pub mod strength;

/// Refusals do not mutate accepted contact history.
#[derive(Debug, Clone, PartialEq)]
pub enum HealingError {
    /// Invalid table, time, identity or nonrepresentable arithmetic.
    Invalid(&'static str),
    /// Temperature outside the supplied mobility domain; never extrapolated.
    Temperature(f64),
    /// Integrated dose outside the supplied strength-map domain.
    Dose(f64),
}

impl fmt::Display for HealingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "bond healing refused: {self:?}")
    }
}
impl std::error::Error for HealingError {}

/// A replaceable, immutable temperature-history law. Implementations must be
/// pure and bounded; cloning a consumer must not share mutable physical history.
/// The supplied temperature path is linear on this interval, not an assertion
/// about an unobserved physical path. The increment must be finite/nonnegative.
pub trait HealingLaw: fmt::Debug {
    /// Application-owned material/law/version identity; not a validation claim.
    fn identity(&self) -> &str;
    /// Dimensionless dose acquired during continuous actual contact.
    fn dose_increment(&self, start_k: f64, end_k: f64, duration_s: f64)
        -> Result<f64, HealingError>;
}

/// A nonnegative mobility rate in inverse seconds at an absolute temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobilityKnot {
    pub temperature_k: f64,
    pub rate_per_s: f64,
}

/// Piecewise-linear mobility, including explicitly supplied zero-rate cold
/// intervals. At most 1024 knots; no WLF/Arrhenius parameters or cutoff defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct MobilityTable {
    identity: String,
    knots: Vec<MobilityKnot>,
}

impl MobilityTable {
    pub fn new(identity: String, knots: Vec<MobilityKnot>) -> Result<Self, HealingError> {
        validate_identity(&identity)?;
        if !(2..=1024).contains(&knots.len()) {
            return Err(HealingError::Invalid("mobility table needs 2..1024 knots"));
        }
        for knot in &knots {
            positive(knot.temperature_k)?;
            nonnegative(knot.rate_per_s)?;
        }
        if knots.windows(2).any(|pair| pair[0].temperature_k >= pair[1].temperature_k) {
            return Err(HealingError::Invalid("temperature knots must strictly increase"));
        }
        Ok(Self { identity, knots })
    }

    pub fn knots(&self) -> &[MobilityKnot] { &self.knots }

    /// Query the same rate used by the integral; inclusive endpoints admitted.
    pub fn rate_per_s(&self, temperature_k: f64) -> Result<f64, HealingError> {
        if !temperature_k.is_finite() || temperature_k < self.knots[0].temperature_k
            || temperature_k > self.knots[self.knots.len()-1].temperature_k
        { return Err(HealingError::Temperature(temperature_k)); }
        let upper = self.knots.partition_point(|k| k.temperature_k < temperature_k);
        if upper == 0 { return Ok(self.knots[0].rate_per_s); }
        let a = self.knots[upper-1];
        let b = self.knots[upper];
        let w = (temperature_k-a.temperature_k)/(b.temperature_k-a.temperature_k);
        nonnegative((1.0-w)*a.rate_per_s + w*b.rate_per_s)
    }
}

impl HealingLaw for MobilityTable {
    fn identity(&self) -> &str { &self.identity }

    fn dose_increment(&self, start_k: f64, end_k: f64, duration_s: f64)
        -> Result<f64, HealingError> {
        nonnegative(duration_s)?;
        let first_rate = self.rate_per_s(start_k)?;
        self.rate_per_s(end_k)?;
        if start_k == end_k { return nonnegative(first_rate*duration_s); }
        // Change variables from time to temperature. Integrating in ascending
        // temperature gives exactly the same dose on the reversed linear path.
        let low = start_k.min(end_k);
        let high = start_k.max(end_k);
        let span = high-low;
        let mut last = low;
        let mut rate = self.rate_per_s(low)?;
        let mut average = 0.0;
        for knot in &self.knots {
            if knot.temperature_k <= low || knot.temperature_k >= high { continue; }
            let weight = (knot.temperature_k-last)/span;
            average = nonnegative(average + weight*(0.5*rate+0.5*knot.rate_per_s))?;
            last = knot.temperature_k;
            rate = knot.rate_per_s;
        }
        average = nonnegative(average + ((high-last)/span)
            *(0.5*rate+0.5*self.rate_per_s(high)?))?;
        nonnegative(average*duration_s)
    }
}

/// Accepted history of ONE material contact patch. Area changes must be
/// represented by separate patches, not by applying an old dose to new area.
#[derive(Debug, Clone)]
pub struct ContactHistory<'a> {
    law: &'a dyn HealingLaw,
    dose: f64,
    contact_time_s: f64,
}

impl<'a> ContactHistory<'a> {
    pub fn new(law: &'a dyn HealingLaw) -> Result<Self, HealingError> {
        validate_identity(law.identity())?;
        Ok(Self { law, dose: 0.0, contact_time_s: 0.0 })
    }
    pub fn law_identity(&self) -> &str { self.law.identity() }
    pub fn dose(&self) -> f64 { self.dose }
    pub fn contact_time_s(&self) -> f64 { self.contact_time_s }

    /// Propose continuous-contact exposure. Only the returned candidate ages;
    /// disconnected intervals must not call this method. A zero duration is an
    /// admitted no-op only for a valid endpoint path with zero law increment.
    pub fn advance(&self, start_k: f64, end_k: f64, duration_s: f64)
        -> Result<Self, HealingError> {
        nonnegative(duration_s)?;
        positive(start_k)?;
        positive(end_k)?;
        let increment = nonnegative(self.law.dose_increment(start_k, end_k, duration_s)?)?;
        if duration_s == 0.0 && increment != 0.0 {
            return Err(HealingError::Invalid("zero-time contact cannot acquire dose"));
        }
        Ok(Self { law: self.law, dose: nonnegative(self.dose+increment)?,
            contact_time_s: nonnegative(self.contact_time_s+duration_s)? })
    }
}

pub(super) fn validate_identity(id: &str) -> Result<(), HealingError> {
    if id.trim().is_empty() || id.len() > 1024 || id.chars().any(char::is_control) {
        Err(HealingError::Invalid("supply a nonempty bounded law/version identity"))
    } else { Ok(()) }
}
pub(super) fn nonnegative(value: f64) -> Result<f64, HealingError> {
    if value.is_finite() && value >= 0.0 { Ok(value) }
    else { Err(HealingError::Invalid("expected a finite nonnegative value")) }
}
pub(super) fn positive(value: f64) -> Result<f64, HealingError> {
    if value.is_finite() && value > 0.0 { Ok(value) }
    else { Err(HealingError::Invalid("expected a finite positive value")) }
}

#[cfg(test)]
mod tests;
