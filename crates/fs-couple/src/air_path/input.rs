//! Source-bound observer-law selection shared by physical pressure performers.
//!
//! The gas record selects no instrument model. Its actual composition and
//! thermodynamic state admit the requested propagation law; the original
//! `ObserverFilter` owns attenuation, dispersion, approximation and latency.

/// Source-gas heat transfer with the same retained composition and material laws.
pub mod convection;

use super::{ObserverFilter, legacy_relative_humidity};
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::{gas::GasState, gas_species::input::GasRecord};

/// An explicit `observer_law` declaration, without an inferred gas identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObserverLawInput {
    /// Explicitly omit atmosphere absorption for a lossless comparison.
    Lossless,
    /// ISO 9613-1, requiring the original gas's terrestrial-air admission.
    Iso9613,
    /// Classical viscothermal absorption from the original gas transport.
    StokesKirchhoff,
    /// The sourced ISO-air molecular relaxation law on a thermochemical gas.
    IsoAirRelaxation,
}

impl ObserverLawInput {
    /// Parse exactly the fields after `observer_law`; unknown suffixes refuse.
    pub fn parse(fields: &[&str]) -> Result<Self, &'static str> {
        match fields {
            ["none"] => Ok(Self::Lossless),
            ["iso9613"] => Ok(Self::Iso9613),
            ["stokes-kirchhoff"] => Ok(Self::StokesKirchhoff),
            ["relaxation", "iso9613-terrestrial-air"] => Ok(Self::IsoAirRelaxation),
            _ => Err(
                "observer_law needs none, iso9613, stokes-kirchhoff or relaxation iso9613-terrestrial-air",
            ),
        }
    }
}

fn iso_admission(
    record: Option<&GasRecord>,
    gas: &GasState,
) -> Result<Option<(fs_material::iso9613::TerrestrialAir, f64)>, AcousticRealizeError> {
    use fs_material::iso9613::{TerrestrialAir, iso9613_relative_humidity};
    match record {
        Some(GasRecord::Legacy(state)) => Ok(legacy_relative_humidity(state).and_then(|rh| {
            let (rebuilt, air) =
                TerrestrialAir::moist_air(state.temperature, state.pressure, rh).ok()?;
            (rebuilt == *gas).then_some((air, rh))
        })),
        Some(GasRecord::Thermochemical(resolved)) => {
            let medium = fs_material::relaxation::RelaxingMedium::from_thermochemical(resolved)
                .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))?;
            let air =
                TerrestrialAir::admit(&medium).map_err(AcousticRealizeError::NonTerrestrialAir)?;
            let rh =
                iso9613_relative_humidity(gas.temperature, gas.pressure, gas.water_mole_fraction)
                    .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))?;
            Ok(Some((air, rh)))
        }
        Some(GasRecord::Mixture(..)) | None => Ok(None),
    }
}

/// Build the original causal observer filter for the actual exterior gas.
///
/// An explicit declaration never falls back to another law. With `None`, the
/// existing wind policy selects ISO only for a legacy air record inside its
/// window, otherwise Stokes--Kirchhoff. A named arbitrary mixture does not gain
/// air-specific coefficients from a similar density, heat capacity or name.
/// Callers own geometric receiver delay and must use the filter's receipt to
/// count low-frequency propagation once. This does not change mechanical gas.
pub fn observer_filter(
    declared: Option<ObserverLawInput>,
    gas: &GasState,
    records: &[GasRecord],
    range_m: f64,
    rate: u32,
) -> Result<ObserverFilter, AcousticRealizeError> {
    let record = records.iter().find(|record| record.state() == gas);
    match declared {
        Some(ObserverLawInput::Lossless) => ObserverFilter::lossless(gas, range_m, rate),
        Some(ObserverLawInput::StokesKirchhoff) => {
            ObserverFilter::stokes_kirchhoff(gas, range_m, rate)
        }
        Some(ObserverLawInput::Iso9613) => match iso_admission(record, gas)? {
            Some((air, rh)) => ObserverFilter::iso9613(&air, gas, rh, range_m, rate),
            None => Err(AcousticRealizeError::InvalidDescription {
                what: "observer_law iso9613 needs a legacy or thermochemical observer gas record",
            }),
        },
        Some(ObserverLawInput::IsoAirRelaxation) => match record {
            Some(GasRecord::Thermochemical(resolved)) => {
                fs_material::relaxation::RelaxingMedium::from_thermochemical(resolved)
                    .and_then(|medium| {
                        fs_material::relaxation::RelaxationAbsorption::iso9613_terrestrial_air()?
                            .evaluate(&medium)
                    })
                    .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))
                    .and_then(|state| ObserverFilter::molecular_relaxation(&state, range_m, rate))
            }
            _ => Err(AcousticRealizeError::InvalidDescription {
                what: "observer_law relaxation needs a thermochemical observer gas record",
            }),
        },
        None => match iso_admission(record, gas) {
            Ok(Some((air, rh)))
                if matches!(record, Some(GasRecord::Legacy(_)))
                    && fs_material::iso9613::iso9613_absorption(&air, gas, rh, 1.0).is_ok() =>
            {
                ObserverFilter::iso9613(&air, gas, rh, range_m, rate)
            }
            _ => ObserverFilter::stokes_kirchhoff(gas, range_m, rate),
        },
    }
}
