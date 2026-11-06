//! Exact sensible calorics retained from the original affine Cp(T) intervals.

use super::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_blake3::DomainHasher;

impl EquilibriumEnthalpyPhaseCurve {
    /// Original endpoint heat capacities [J/(kg K)] on one sensible interval.
    ///
    /// `Some` means the interval integrates this affine Cp(T) law exactly;
    /// its h/T endpoint chord is only a boundary description. `None` means an
    /// explicitly authored affine T(h) interval, a latent interval, or an index
    /// outside the chart. Endpoints are `knots()[index]` and `knots()[index+1]`.
    /// These values carry the same source context and identity as the curve.
    #[must_use]
    pub fn sensible_heat_capacities_on_interval(&self, index: usize) -> Option<[f64; 2]> {
        self.sensible_heat_capacities.get(index).copied().flatten()
    }

    /// Minimum and maximum specific heat capacity [J/(kg K)] on one sensible
    /// interval. Sourced spans retain their original affine Cp(T) endpoint
    /// extrema; authored sensible spans return their constant h/T chord
    /// capacity twice. No average replaces a source extremum.
    ///
    /// Refuses an invalid interval, an isothermal latent span, or a capacity
    /// that is nonpositive or unrepresentable. The interval retains the
    /// original source support and thermodynamic context of this curve.
    pub fn sensible_heat_capacity_bounds_j_kg_k(
        &self,
        index: usize,
    ) -> Result<[f64; 2], PhaseStateError> {
        let upper = index.checked_add(1).and_then(|index| self.knots.get(index));
        let (Some(lower), Some(upper)) = (self.knots.get(index), upper) else {
            return Err(invalid("sensible interval is outside the original chart"));
        };
        let width = upper.temperature_k - lower.temperature_k;
        if !width.is_finite() || width <= 0.0 {
            return Err(invalid(
                "an isothermal interval has no finite sensible heat capacity",
            ));
        }
        let cp = self
            .sensible_heat_capacities_on_interval(index)
            .unwrap_or_else(|| {
                [(upper.specific_enthalpy_j_kg - lower.specific_enthalpy_j_kg) / width; 2]
            });
        if cp.iter().any(|value| !value.is_finite() || *value <= 0.0) {
            return Err(invalid(
                "sensible heat capacity must be finite and positive",
            ));
        }
        Ok([cp[0].min(cp[1]), cp[0].max(cp[1])])
    }

    /// Resolve h(T) on one explicitly selected sensible interval.
    ///
    /// Sourced intervals integrate their original affine Cp(T); authored
    /// sensible intervals retain their affine h/T law. Exact endpoints return
    /// the stored enthalpy bits. Interior integration is relative to the nearer
    /// temperature endpoint to avoid subtracting large caloric primitives.
    /// An isothermal plateau has no unique h(T) and is refused, as are source
    /// extrapolation and an interior enthalpy lost at the chosen datum.
    pub fn specific_enthalpy_on_sensible_interval(
        &self,
        index: usize,
        temperature_k: f64,
    ) -> Result<f64, PhaseStateError> {
        let upper = index.checked_add(1).and_then(|index| self.knots.get(index));
        let (Some(lower), Some(upper)) = (self.knots.get(index), upper) else {
            return Err(invalid("sensible interval is outside the original chart"));
        };
        let width = upper.temperature_k - lower.temperature_k;
        if width <= 0.0 || !width.is_finite() {
            return Err(invalid(
                "an isothermal interval has no unique temperature inverse",
            ));
        }
        if !temperature_k.is_finite() {
            return Err(PhaseStateError::NonFiniteTemperature);
        }
        if temperature_k < lower.temperature_k || temperature_k > upper.temperature_k {
            return Err(PhaseStateError::OutsideTemperatureDomain {
                temperature_k,
                lower_k: lower.temperature_k,
                upper_k: upper.temperature_k,
            });
        }
        if temperature_k == lower.temperature_k {
            return Ok(lower.specific_enthalpy_j_kg);
        }
        if temperature_k == upper.temperature_k {
            return Ok(upper.specific_enthalpy_j_kg);
        }
        let fraction = (temperature_k - lower.temperature_k) / width;
        let enthalpy = if let Some(cp) = self.sensible_heat_capacities_on_interval(index) {
            let capacity = capacity_at_fraction(cp, fraction);
            if fraction <= 0.5 {
                lower.specific_enthalpy_j_kg
                    + (0.5 * cp[0] + 0.5 * capacity) * (temperature_k - lower.temperature_k)
            } else {
                upper.specific_enthalpy_j_kg
                    - (0.5 * cp[1] + 0.5 * capacity) * (upper.temperature_k - temperature_k)
            }
        } else {
            interpolate(
                lower.specific_enthalpy_j_kg,
                upper.specific_enthalpy_j_kg,
                fraction,
            )
        };
        if !enthalpy.is_finite()
            || enthalpy <= lower.specific_enthalpy_j_kg
            || enthalpy >= upper.specific_enthalpy_j_kg
        {
            return Err(invalid(
                "interior sensible temperature has no representable enthalpy",
            ));
        }
        Ok(enthalpy)
    }

    pub(super) fn retain_sensible_heat_capacities(
        &mut self,
        capacities: Vec<Option<[f64; 2]>>,
    ) -> Result<(), PhaseStateError> {
        if capacities.len() != self.knots.len() - 1 {
            return Err(invalid(
                "retained heat capacities must cover every original interval",
            ));
        }
        for (index, cp) in capacities.iter().enumerate() {
            if let Some(cp) = cp {
                let lower = self.knots[index];
                let upper = self.knots[index + 1];
                if cp.iter().any(|value| !value.is_finite() || *value <= 0.0)
                    || upper.temperature_k <= lower.temperature_k
                    || lower.liquid_mass_fraction != upper.liquid_mass_fraction
                    || (lower.liquid_mass_fraction != 0.0 && lower.liquid_mass_fraction != 1.0)
                {
                    return Err(invalid(
                        "retained heat capacity needs a positive pure-phase sensible law",
                    ));
                }
            }
        }
        let mut identity = DomainHasher::new(
            "org.frankensim.fs-material.sensible-heat-capacity-affine-temperature.v1",
        );
        identity.update(self.identity.as_bytes());
        identity.update(&(capacities.len() as u64).to_le_bytes());
        for cp in &capacities {
            identity.update(&[u8::from(cp.is_some())]);
            if let Some(cp) = cp {
                for value in cp {
                    identity.update(&value.to_bits().to_le_bytes());
                }
            }
        }
        self.identity = identity.finalize();
        self.sensible_heat_capacities = capacities;
        Ok(())
    }
}

/// Positive endpoint interpolation; using the nearer endpoint avoids losing a
/// small positive capacity next to a much larger one.
pub(super) fn capacity_at_fraction(cp: [f64; 2], fraction: f64) -> f64 {
    interpolate(cp[0], cp[1], fraction)
}

pub(super) fn interpolate(lower: f64, upper: f64, fraction: f64) -> f64 {
    if fraction <= 0.5 {
        (upper - lower).mul_add(fraction, lower)
    } else {
        (lower - upper).mul_add(1.0 - fraction, upper)
    }
}

/// Original temperature coordinate on one admitted affine-Cp interval.
pub(super) fn temperature_fraction(
    lower: EnthalpyPhaseKnot,
    upper: EnthalpyPhaseKnot,
    cp: [f64; 2],
    h: f64,
) -> Result<f64, PhaseStateError> {
    temperature_fraction_from_bounds(
        [lower.specific_enthalpy_j_kg, upper.specific_enthalpy_j_kg],
        upper.temperature_k - lower.temperature_k,
        cp,
        h,
    )
}

pub(super) fn temperature_fraction_from_bounds(
    enthalpies: [f64; 2],
    width: f64,
    cp: [f64; 2],
    h: f64,
) -> Result<f64, PhaseStateError> {
    if h == enthalpies[0] {
        return Ok(0.0);
    }
    if h == enthalpies[1] {
        return Ok(1.0);
    }
    let from_lower = h - enthalpies[0];
    let from_upper = enthalpies[1] - h;
    if from_lower <= from_upper {
        invert_capacity_fraction(from_lower, width, cp)
    } else {
        invert_capacity_fraction(from_upper, width, [cp[1], cp[0]]).map(|fraction| 1.0 - fraction)
    }
}

// Solve dh = width * (c0*u + (c1-c0)*u^2/2), without forming c0^2,
// subtracting quadratic roots or multiplying width by a large capacity.
// Choosing the nearer enthalpy endpoint keeps a decreasing-Cp discriminant
// away from cancellation even when the opposite capacity is very small.
fn invert_capacity_fraction(energy: f64, width: f64, cp: [f64; 2]) -> Result<f64, PhaseStateError> {
    let fraction = if cp[0] == cp[1] {
        (energy / cp[0]) / width
    } else {
        let scale = cp[0].max(cp[1]);
        let a = cp[0] / scale;
        let b = cp[1] / scale;
        let energy = (energy / scale) / width;
        let discriminant = (2.0 * (b - a)).mul_add(energy, a * a);
        (2.0 * energy) / (a + discriminant.sqrt())
    };
    if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
        return Err(invalid(
            "affine heat-capacity inverse is not representable inside its source interval",
        ));
    }
    Ok(fraction)
}

fn invalid(what: &'static str) -> PhaseStateError {
    PhaseStateError::InvalidCurve { what }
}

#[cfg(test)]
mod tests;
