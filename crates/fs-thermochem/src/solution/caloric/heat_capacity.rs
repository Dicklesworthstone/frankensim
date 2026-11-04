//! Shared numerical integration of original positive molar heat capacities.

use super::{BinaryCaloricErrorV1, Sum, checked, interpolate, positive};
use crate::{MolarThermalQuantityV1, SpeciesId};
use fs_qty::Temperature;

/// Maximum source points per molar heat-capacity curve.
pub const MAX_CALORIC_KNOTS_V1: usize = 4096;

/// Numerical positive molar Cp(T), in J/(mol K), with a fixed sensible datum.
/// This curve has no phase qualification: its caller must bind actual source
/// evidence to a liquid, solid or dissolved species. It never turns pure-liquid
/// measurements into solid data. Linear interpolation preserves positive Cp.
#[derive(Clone, Debug, PartialEq)]
pub struct MolarHeatCapacityCurveV1 {
    species: SpeciesId,
    reference_temperature: Temperature,
    pub(super) knots: Vec<(f64, f64)>,
}
impl MolarHeatCapacityCurveV1 {
    /// Retain 2..=4096 strictly increasing positive temperature knots. All Cp
    /// values must be positive and finite. The datum must lie inside the curve.
    pub fn new(
        species: SpeciesId,
        reference_temperature: Temperature,
        knots: Vec<(Temperature, MolarThermalQuantityV1)>,
    ) -> Result<Self, BinaryCaloricErrorV1> {
        if !(2..=MAX_CALORIC_KNOTS_V1).contains(&knots.len())
            || knots
                .iter()
                .any(|(t, cp)| !positive(t.value()) || !positive(cp.value()))
            || knots.windows(2).any(|w| w[0].0.value() >= w[1].0.value())
            || !reference_temperature.value().is_finite()
            || reference_temperature.value() < knots[0].0.value()
            || reference_temperature.value() > knots[knots.len() - 1].0.value()
        {
            return Err(BinaryCaloricErrorV1::Invalid("molar Cp curve or datum"));
        }
        Ok(Self {
            species,
            reference_temperature,
            knots: knots
                .into_iter()
                .map(|(t, cp)| (t.value(), cp.value()))
                .collect(),
        })
    }

    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }
    pub const fn reference_temperature(&self) -> Temperature {
        self.reference_temperature
    }
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        [
            Temperature::new(self.knots[0].0),
            Temperature::new(self.knots[self.knots.len() - 1].0),
        ]
    }

    /// Exact integral of the declared piecewise-linear Cp interpolation [J/mol].
    /// Integrates only crossed segments, using their local temperature spans;
    /// cooling is the exact negative of the same ascending integration.
    pub fn enthalpy_change_j_per_mol(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, BinaryCaloricErrorV1> {
        self.admit(from.value())?;
        self.admit(to.value())?;
        let (lo, hi, sign) = if from.value() <= to.value() {
            (from.value(), to.value(), 1.0)
        } else {
            (to.value(), from.value(), -1.0)
        };
        if lo == hi {
            return Ok(0.0);
        }
        let start = self
            .knots
            .partition_point(|&(t, _)| t <= lo)
            .saturating_sub(1);
        let mut sum = Sum::default();
        for pair in self.knots[start..].windows(2) {
            let a = lo.max(pair[0].0);
            let b = hi.min(pair[1].0);
            if b > a {
                let cp_a = interpolate(pair, a);
                let cp_b = interpolate(pair, b);
                let increment =
                    checked((cp_a + 0.5 * (cp_b - cp_a)) * (b - a), "sensible increment")?;
                if increment <= 0.0 {
                    return Err(BinaryCaloricErrorV1::Arithmetic(
                        "positive sensible increment underflows",
                    ));
                }
                sum.add(increment)?;
            }
            if pair[1].0 >= hi {
                break;
            }
        }
        checked(sign * sum.value()?, "sensible integral")
    }

    pub(super) fn admit(&self, t: f64) -> Result<(), BinaryCaloricErrorV1> {
        if !t.is_finite() || t < self.knots[0].0 || t > self.knots[self.knots.len() - 1].0 {
            Err(BinaryCaloricErrorV1::Invalid(
                "temperature outside molar Cp source",
            ))
        } else {
            Ok(())
        }
    }
}
