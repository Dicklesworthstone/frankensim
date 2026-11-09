//! One original caloric/density pair across thermal, property and mechanical users.
//!
//! Exact source Cp integration must not be converted to a sampled enthalpy chart
//! just to use crystallinity-dependent properties or shrinkage. This borrowed
//! view retains either supported representation and its original identities.
//! It owns no kinetics, heat, pressure work or alternative density interpolation.

use super::{CrystallizationCaloricPair, CrystallizationVolumeState, invalid};
use crate::phase::PhaseStateError;
use crate::thermal::transport::ThermalTransportError;
use crate::thermal::transport::isobaric::crystallization::IsobaricCrystallizationPair;
use fs_blake3::ContentHash;
use fs_matdb::QueryPoint;
use fs_qty::QuantitySpec;

/// Original source failures, without replacing missing density by a default.
#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationSourceError {
    Phase(PhaseStateError),
    Transport(ThermalTransportError),
}
impl From<PhaseStateError> for CrystallizationSourceError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<ThermalTransportError> for CrystallizationSourceError {
    fn from(error: ThermalTransportError) -> Self {
        Self::Transport(error)
    }
}
impl core::fmt::Display for CrystallizationSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystallization source: {self:?}")
    }
}
impl std::error::Error for CrystallizationSourceError {}

/// Allocation-free access to the original phase pair; no charts are resampled.
#[derive(Clone, Copy, Debug)]
pub enum CrystallizationSourcePair<'a> {
    Chart(CrystallizationCaloricPair<'a>),
    ExactCp(&'a IsobaricCrystallizationPair<'a>),
}

impl<'a> From<CrystallizationCaloricPair<'a>> for CrystallizationSourcePair<'a> {
    fn from(pair: CrystallizationCaloricPair<'a>) -> Self {
        Self::Chart(pair)
    }
}
impl<'a> From<&CrystallizationCaloricPair<'a>> for CrystallizationSourcePair<'a> {
    fn from(pair: &CrystallizationCaloricPair<'a>) -> Self {
        Self::Chart(*pair)
    }
}
impl<'a> From<&'a IsobaricCrystallizationPair<'a>> for CrystallizationSourcePair<'a> {
    fn from(pair: &'a IsobaricCrystallizationPair<'a>) -> Self {
        Self::ExactCp(pair)
    }
}
impl<'a> From<&CrystallizationSourcePair<'a>> for CrystallizationSourcePair<'a> {
    fn from(pair: &CrystallizationSourcePair<'a>) -> Self {
        *pair
    }
}

impl CrystallizationSourcePair<'_> {
    pub fn identity(&self) -> ContentHash {
        match self {
            Self::Chart(pair) => pair.identity(),
            Self::ExactCp(pair) => pair.identity(),
        }
    }
    pub fn material_card_identity(&self) -> ContentHash {
        match self {
            Self::Chart(pair) => pair.parent().curve().material_card_identity(),
            Self::ExactCp(pair) => pair.parent().source().material_card_identity(),
        }
    }
    pub fn parent_identity(&self) -> ContentHash {
        match self {
            Self::Chart(pair) => pair.parent().identity(),
            Self::ExactCp(pair) => pair.parent().identity(),
        }
    }
    pub fn crystal_identity(&self) -> ContentHash {
        match self {
            Self::Chart(pair) => pair.crystal().identity(),
            Self::ExactCp(pair) => pair.crystal().identity(),
        }
    }
    pub fn temperature_range_k(&self) -> [f64; 2] {
        match self {
            Self::Chart(pair) => pair.temperature_range_k(),
            Self::ExactCp(pair) => pair.temperature_range_k(),
        }
    }
    /// Explicit fixed source pressure, when retained by the caloric owner.
    /// A chart pair still requires the caller's original pressure admission.
    pub fn pressure_pa(&self) -> Option<f64> {
        match self {
            Self::Chart(_) => None,
            Self::ExactCp(pair) => Some(pair.pressure_pa()),
        }
    }

    /// The exact-Cp owner explicitly admits this fixed pressure coordinate,
    /// including dimensional Pa. Preserve that role rather than inferring a
    /// pressure role from arbitrary stress-valued property coordinates.
    pub(crate) fn pressure_coordinate(&self) -> Option<(&str, QuantitySpec)> {
        let Self::ExactCp(pair) = self else {
            return None;
        };
        let source = pair.parent().source();
        let axis = source.pressure_axis();
        source.query_endpoints()[0]
            .axis_quantities()
            .get(axis)
            .copied()
            .map(|quantity| (axis, quantity))
    }

    /// Preserve both original density laws at the actual kinetic T and mass X.
    /// Calorics-only exact-Cp pairs remain valid for heat evolution but refuse
    /// here: neither free volume nor shrinkage follows from Cp and latent heat.
    pub fn volume_state(
        &self,
        temperature_k: f64,
        crystalline_mass_fraction: f64,
    ) -> Result<CrystallizationVolumeState, CrystallizationSourceError> {
        match self {
            Self::Chart(pair) => Ok(pair.volume_state(temperature_k, crystalline_mass_fraction)?),
            Self::ExactCp(pair) => {
                let [low, high] = pair.temperature_range_k();
                if !temperature_k.is_finite() || temperature_k < low || temperature_k > high {
                    return Err(invalid("temperature is outside the original exact-Cp pair").into());
                }
                let parent = pair
                    .parent()
                    .density_at_temperature(temperature_k)?
                    .ok_or_else(|| invalid("original exact-Cp parent density is unavailable"))?;
                let crystal = pair
                    .crystal()
                    .density_at_temperature(temperature_k)?
                    .ok_or_else(|| invalid("original exact-Cp crystal density is unavailable"))?;
                Ok(CrystallizationVolumeState::from_densities(
                    pair.identity(),
                    temperature_k,
                    crystalline_mass_fraction,
                    parent,
                    crystal,
                )?)
            }
        }
    }

    /// Exact-Cp sources retain their actual pressure and composition context.
    /// Additional kinetic fraction coordinates may be supplied by a consumer;
    /// none may replace a fixed coordinate already owned by the caloric source.
    pub(crate) fn admit_property_context(
        &self,
        query: &QueryPoint,
        temperature_axis: &str,
        fraction_axes: [Option<&str>; 2],
    ) -> Result<(), CrystallizationSourceError> {
        if let Self::ExactCp(pair) = self {
            let source = pair.parent().source();
            let original = &source.query_endpoints()[0];
            if temperature_axis != source.temperature_axis()
                || original.axes().iter().any(|(axis, value)| {
                    axis != temperature_axis
                        && (fraction_axes.contains(&Some(axis.as_str()))
                            || query.axes().get(axis) != Some(value)
                            || query.axis_quantities().get(axis)
                                != original.axis_quantities().get(axis))
                })
            {
                return Err(invalid(
                    "properties must retain exact-Cp temperature axis and fixed source context",
                )
                .into());
            }
        }
        Ok(())
    }
}
