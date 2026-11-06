//! Identified temperature-dependent microscopic reaction coefficients.
//!
//! Each immutable row contains already processed group cross sections [m²] at
//! one absolute temperature [K]. Between rows, each coefficient is interpolated
//! linearly in Kelvin. A rate is the deterministic dot product with the supplied
//! absolute group-integrated scalar flux [1/(m² s)]. The within-group spectrum
//! used to process the microscopic rows remains part of the supplied data.
//!
//! This evaluator does not Doppler broaden resonances, solve self-shielding,
//! change the neutron spectrum, or establish evaluated-data accuracy. Its
//! interpolation is a declared deterministic multigroup law, distinct from the
//! statistical temperature interpolation used in Monte Carlo transport; see
//! <https://docs.openmc.org/en/latest/methods/cross_sections.html>.

use std::sync::Arc;

use crate::irradiation::{GroupFlux, NeutronGroups, NuclideId};
use crate::{Accumulator, FissionError, finite, nonnegative, positive};

/// Physical reaction whose microscopic coefficients are represented.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NeutronReaction {
    Fission,
    /// Radiative neutron capture, ENDF MT102.
    Capture,
    /// Nonfission (n,p) absorption without neutron-producing breakup, MT103.
    ProtonProduction,
    /// Nonfission (n,d) absorption without neutron-producing breakup, MT104.
    DeuteronProduction,
    /// Nonfission (n,t) absorption without neutron-producing breakup, MT105.
    TritonProduction,
    /// Nonfission (n,He-3) absorption without neutron-producing breakup, MT106.
    Helium3Production,
    /// Nonfission (n,alpha) absorption without neutron-producing breakup, MT107.
    AlphaProduction,
}

impl NeutronReaction {
    /// Channels that remove the incident neutron and parent without producing
    /// another neutron in this model. Products and energy are supplied separately.
    /// Charged-particle MTs are the absorption components in ENDF-102 §0.4.3.2;
    /// neither a total particle-production tally nor a breakup model is inferred.
    /// <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>
    #[must_use]
    pub const fn is_nonfission_absorption(self) -> bool {
        match self {
            Self::Fission => false,
            Self::Capture
            | Self::ProtonProduction
            | Self::DeuteronProduction
            | Self::TritonProduction
            | Self::Helium3Production
            | Self::AlphaProduction => true,
        }
    }
}

/// Caller-supplied processed rows on one shared group grid. A source label does
/// not by itself certify the evaluation, temperature processing or collapse.
#[derive(Clone, Debug, PartialEq)]
pub struct CrossSectionData {
    pub data_id: String,
    pub nuclide: NuclideId,
    pub reaction: NeutronReaction,
    pub groups: Arc<NeutronGroups>,
    /// Strictly increasing, finite, positive absolute temperatures [K].
    pub temperatures_k: Vec<f64>,
    /// One nonnegative microscopic coefficient [m²] per temperature and group.
    pub cross_sections_m2: Vec<Vec<f64>>,
}

/// Explicit input-size and retained numerical-payload limits. Validation
/// checks these before allocating the flattened coefficient storage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CrossSectionLimits {
    pub max_temperatures: usize,
    pub max_neutron_groups: usize,
    pub max_values: usize,
}

impl Default for CrossSectionLimits {
    fn default() -> Self {
        Self {
            max_temperatures: 256,
            max_neutron_groups: 4096,
            max_values: 1_048_576,
        }
    }
}

/// Immutable, source-identified microscopic rows. Clone shares the group grid,
/// temperatures and contiguous coefficient payload. Evaluating one group uses
/// binary-search bracketing; flux collapse brackets once and allocates nothing.
#[derive(Clone, Debug, PartialEq)]
pub struct CrossSectionTable {
    data_id: Arc<str>,
    nuclide: NuclideId,
    reaction: NeutronReaction,
    groups: Arc<NeutronGroups>,
    temperatures_k: Arc<[f64]>,
    values_m2: Arc<[f64]>,
}

#[derive(Clone, Copy)]
enum Bracket {
    Knot(usize),
    Interior {
        lower: usize,
        upper_weight: f64,
        lower_weight: f64,
    },
}

impl CrossSectionTable {
    /// Admit all identities, row shapes, values and limits before publication.
    /// A single temperature is valid only at that exact temperature.
    pub fn new(data: CrossSectionData, limits: CrossSectionLimits) -> Result<Self, FissionError> {
        Self::new_with_checkpoint(data, limits, || Ok(()))
    }

    /// Cancellable construction, polling between rows and every 64 values.
    pub fn new_with_checkpoint(
        data: CrossSectionData,
        limits: CrossSectionLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, FissionError> {
        checkpoint()?;
        let rows = data.temperatures_k.len();
        let groups = data.groups.len();
        if data.data_id.trim().is_empty() || rows == 0 || rows != data.cross_sections_m2.len() {
            return Err(FissionError::InvalidInput(
                "identified nonempty temperature and cross-section rows required",
            ));
        }
        if rows > limits.max_temperatures {
            return Err(FissionError::Budget("cross-section temperatures"));
        }
        if groups > limits.max_neutron_groups {
            return Err(FissionError::Budget("neutron groups"));
        }
        let count = rows
            .checked_mul(groups)
            .ok_or(FissionError::Budget("cross-section table size"))?;
        if count > limits.max_values {
            return Err(FissionError::Budget("cross-section table values"));
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(count)
            .map_err(|_| FissionError::Budget("cross-section table memory"))?;
        let mut previous_temperature = 0.0;
        for (&temperature, row) in data.temperatures_k.iter().zip(&data.cross_sections_m2) {
            checkpoint()?;
            if !positive(temperature) || temperature <= previous_temperature || row.len() != groups
            {
                return Err(FissionError::InvalidInput(
                    "ascending positive temperatures and one coefficient per group required",
                ));
            }
            previous_temperature = temperature;
            for (index, &value) in row.iter().enumerate() {
                if index % 64 == 0 {
                    checkpoint()?;
                }
                if !nonnegative(value) {
                    return Err(FissionError::InvalidInput(
                        "finite nonnegative microscopic cross sections required",
                    ));
                }
                values.push(value);
            }
        }
        checkpoint()?;
        Ok(Self {
            data_id: data.data_id.into(),
            nuclide: data.nuclide,
            reaction: data.reaction,
            groups: data.groups,
            temperatures_k: data.temperatures_k.into(),
            values_m2: values.into(),
        })
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    #[must_use]
    pub const fn nuclide(&self) -> NuclideId {
        self.nuclide
    }
    #[must_use]
    pub const fn reaction(&self) -> NeutronReaction {
        self.reaction
    }
    #[must_use]
    pub fn groups(&self) -> &Arc<NeutronGroups> {
        &self.groups
    }
    #[must_use]
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    #[must_use]
    pub fn temperature_min_k(&self) -> f64 {
        self.temperatures_k[0]
    }
    #[must_use]
    pub fn temperature_max_k(&self) -> f64 {
        self.temperatures_k[self.temperatures_k.len() - 1]
    }

    /// No extrapolation, tolerance-based snapping or temperature clamping.
    pub fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionError> {
        if !positive(temperature_k)
            || temperature_k < self.temperature_min_k()
            || temperature_k > self.temperature_max_k()
        {
            return Err(FissionError::InvalidInput(
                "temperature is outside the cross-section table",
            ));
        }
        Ok(())
    }

    /// Evaluate a microscopic coefficient [m²] without allocating a row.
    pub fn cross_section_at_group(
        &self,
        temperature_k: f64,
        group: usize,
    ) -> Result<f64, FissionError> {
        if group >= self.groups.len() {
            return Err(FissionError::InvalidInput(
                "cross-section group index is out of range",
            ));
        }
        self.evaluate(self.bracket(temperature_k)?, group)
    }

    /// Per-parent reaction hazard [1/s] for one absolute spectrum. Full named
    /// grid equality and temperature applicability are mandatory even at zero
    /// flux. Work is O(log(temperatures) + groups), with a poll for each group.
    /// Rounded microscopic interpolation and products define the returned
    /// coefficient; this is not a data-uncertainty or forward-error enclosure.
    pub fn rate_constant(
        &self,
        flux: &GroupFlux,
        temperature_k: f64,
        max_neutron_groups: usize,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<f64, FissionError> {
        checkpoint()?;
        if self.groups.len() > max_neutron_groups {
            return Err(FissionError::Budget("neutron groups"));
        }
        if !Arc::ptr_eq(&self.groups, flux.groups()) && self.groups != *flux.groups() {
            return Err(FissionError::InvalidInput(
                "neutron flux and cross-section table groups do not match",
            ));
        }
        let bracket = self.bracket(temperature_k)?;
        let mut rate = Accumulator::new(0.0);
        for (group, &flux) in flux.values_m2_s().iter().enumerate() {
            checkpoint()?;
            // A zero-flux group has exactly zero hazard; it cannot amplify an
            // otherwise unrepresentable microscopic interpolation contribution.
            if flux == 0.0 {
                continue;
            }
            let cross_section = self.evaluate(bracket, group)?;
            let product = finite(cross_section * flux)?;
            if cross_section > 0.0 && product == 0.0 {
                return Err(FissionError::Arithmetic);
            }
            rate.add(product)?;
        }
        checkpoint()?;
        finite(rate.value())
    }

    fn bracket(&self, temperature_k: f64) -> Result<Bracket, FissionError> {
        self.validate_temperature(temperature_k)?;
        let upper = self
            .temperatures_k
            .partition_point(|&value| value < temperature_k);
        if self.temperatures_k[upper] == temperature_k {
            return Ok(Bracket::Knot(upper));
        }
        let lower = upper - 1;
        let low = self.temperatures_k[lower];
        let high = self.temperatures_k[upper];
        let width = high - low;
        // Compute both distances directly. Forming a small weight as 1-w can
        // erase a representable contribution near a temperature endpoint.
        Ok(Bracket::Interior {
            lower,
            upper_weight: (temperature_k - low) / width,
            lower_weight: (high - temperature_k) / width,
        })
    }

    fn evaluate(&self, bracket: Bracket, group: usize) -> Result<f64, FissionError> {
        let groups = self.groups.len();
        match bracket {
            Bracket::Knot(row) => Ok(self.values_m2[row * groups + group]),
            Bracket::Interior {
                lower,
                upper_weight,
                lower_weight,
            } => {
                let a = self.values_m2[lower * groups + group];
                let b = self.values_m2[(lower + 1) * groups + group];
                if a == b {
                    return Ok(a);
                }
                // Start at the nearer endpoint. This retains a small endpoint
                // coefficient beside a much larger coefficient at the other
                // knot, without adding two potentially overflowing terms.
                let (base, other, weight) = if upper_weight <= lower_weight {
                    (a, b, upper_weight)
                } else {
                    (b, a, lower_weight)
                };
                if weight == 0.0 {
                    return Err(FissionError::Arithmetic);
                }
                let result = finite((other - base).mul_add(weight, base))?;
                if !nonnegative(result) || result == 0.0 {
                    // Distinct nonnegative endpoints give a strictly positive
                    // interior value. Refuse an underflowed positive quantity.
                    return Err(FissionError::Arithmetic);
                }
                Ok(result)
            }
        }
    }
}

#[cfg(test)]
mod tests;
