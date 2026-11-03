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
use fs_material::relaxation::{
    RelaxationAbsorption,
    input::{medium_from_record, read_relaxation},
};

/// An explicit `observer_law` declaration, without an inferred gas identity.
#[derive(Clone, Debug, PartialEq)]
pub enum ObserverLawInput {
    /// Explicitly omit atmosphere absorption for a lossless comparison.
    Lossless,
    /// ISO 9613-1, requiring the original gas's terrestrial-air admission.
    Iso9613,
    /// Classical viscothermal absorption from the original gas transport.
    StokesKirchhoff,
    /// The original species/reservoir/rate law, evaluated against the actual
    /// exterior's explicit mixture or thermochemical record. Declared rates
    /// remain authored constitutive data, not a measurement or calibration.
    MolecularRelaxation(RelaxationAbsorption),
}

impl ObserverLawInput {
    /// Parse exactly the fields after `observer_law`; unknown suffixes refuse.
    pub fn parse(fields: &[&str]) -> Result<Self, &'static str> {
        match fields {
            ["none"] => Ok(Self::Lossless),
            ["iso9613"] => Ok(Self::Iso9613),
            ["stokes-kirchhoff"] => Ok(Self::StokesKirchhoff),
            ["relaxation", declaration @ ..] => {
                // Preparation only: reuse the bounded material grammar. No
                // observer-specific species or empirical coefficients exist.
                let source = declaration.join(" ");
                let mut row = source.split_ascii_whitespace();
                let law = read_relaxation(&mut row)
                    .map_err(|_| "invalid observer molecular-relaxation declaration")?;
                if row.next().is_some() {
                    return Err("unexpected field after observer molecular-relaxation declaration");
                }
                Ok(Self::MolecularRelaxation(law))
            }
            _ => Err(
                "observer_law needs none, iso9613, stokes-kirchhoff or relaxation with iso9613-terrestrial-air or declared species data",
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
/// A declared molecular law uses the original material input grammar, including
/// its explicit zero-reservoir species, partner completeness and rate domains.
/// It contributes absorption and dispersion once through the original filter;
/// no additional ISO or classical filter is applied.
pub fn observer_filter(
    declared: Option<&ObserverLawInput>,
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
        Some(ObserverLawInput::MolecularRelaxation(law)) => match record {
            Some(record) => {
                medium_from_record(record)
                    .and_then(|medium| law.evaluate(&medium))
                    .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))
                    .and_then(|state| ObserverFilter::molecular_relaxation(&state, range_m, rate))
            }
            _ => Err(AcousticRealizeError::InvalidDescription {
                what: "observer_law relaxation needs the matching explicit-species observer gas record",
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::air_path::GasAbsorptionModel;
    use fs_material::gas_species::input::read_gas_record;

    // Synthetic constitutive data for adapter checks, not measured CO2 rates.
    const DECLARATION: &str = "relaxation declared 3 N2 0 O2 0 CO2 1 vibration 1 960 2 power 10 293.15 0 200 400";

    fn parse(text: &str) -> Result<ObserverLawInput, &'static str> {
        ObserverLawInput::parse(&text.split_ascii_whitespace().collect::<Vec<_>>())
    }

    fn mixture(temperature: f64) -> GasRecord {
        let text = format!("mixture {temperature} 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4");
        read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap()
    }

    /// G3: the source declaration reaches the SAME complex material law and
    /// causal observer, with no second classical or ISO absorption.
    #[test]
    fn declared_mixture_observer_is_the_original_material_transfer() {
        let declaration = parse(DECLARATION).unwrap();
        let record = mixture(293.15);
        let ObserverLawInput::MolecularRelaxation(law) = &declaration else {
            panic!("the molecular source was dropped");
        };
        let state = law.evaluate(&medium_from_record(&record).unwrap()).unwrap();
        let original = ObserverFilter::molecular_relaxation(&state, 1.0, 8000).unwrap();
        let actual = observer_filter(
            Some(&declaration), record.state(), core::slice::from_ref(&record), 1.0, 8000,
        ).unwrap();
        assert_eq!(actual.taps(), original.taps());
        assert_eq!(actual.receipt(), original.receipt());
        assert_eq!(actual.receipt().absorption.law, GasAbsorptionModel::MolecularRelaxation);
        assert!(actual.receipt().absorption.relaxation_omitted.is_none());
        let omega = core::f64::consts::TAU * 1000.0;
        assert!(state.causal_wavenumber(omega).unwrap().im
            > record.state().stokes_kirchhoff_absorption(omega));
        let classical = observer_filter(None, record.state(), &[record.clone()], 1.0, 8000).unwrap();
        assert_eq!(classical.receipt().absorption.law, GasAbsorptionModel::StokesKirchhoff);
        assert_ne!(actual.taps(), classical.taps());
    }

    /// G0: a selected law never falls back when its original source is missing,
    /// its coefficients are out of domain, or an ISO air law sees enriched CO2.
    #[test]
    fn declared_observer_retains_species_state_and_material_refusals() {
        let record = mixture(293.15);
        for source in [
            "relaxation declared 2 N2 0 O2 0",
            "relaxation declared 3 N2 0 O2 0 CO2 1 vibration 1 960 2 power 10 293.15 0 310 400",
            "relaxation declared 3 N2 0 O2 0 CO2 1 vibration 1 960 2 partners 1 CO2 2.302585092994046 0 200 400",
            "relaxation iso9613-terrestrial-air",
        ] {
            let law = parse(source).unwrap();
            assert!(observer_filter(Some(&law), record.state(), &[record.clone()], 1.0, 8000).is_err(), "{source}");
        }
        let law = parse(DECLARATION).unwrap();
        assert!(observer_filter(Some(&law), record.state(), &[], 1.0, 8000).is_err());
        let other = mixture(300.0);
        assert!(observer_filter(Some(&law), record.state(), &[other], 1.0, 8000).is_err());
        let legacy = GasRecord::Legacy(GasState::try_new_moist_air(293.15, 100000.0, 0.0).unwrap());
        assert!(observer_filter(Some(&law), legacy.state(), &[legacy.clone()], 1.0, 8000).is_err());
    }

    /// G0: malformed counts and ignored trailing source fields cannot select
    /// a different law or escape the material owner's allocation limits.
    #[test]
    fn observer_declarations_are_complete_and_bounded() {
        for source in [
            "relaxation", "relaxation declared 0", "relaxation declared 129",
            "relaxation declared 1 CO2 33", "relaxation declared 1 CO2 1 vibration 17",
            "relaxation declared 1 CO2 1 vibration 1 NaN 2 power 10 293.15 0 200 400",
            "relaxation iso9613-terrestrial-air ignored", "none ignored",
        ] {
            assert!(parse(source).is_err(), "{source}");
        }
        assert!(parse(&format!("{DECLARATION} ignored")).is_err());
        assert_eq!(parse("none").unwrap(), ObserverLawInput::Lossless);
        assert_eq!(parse("iso9613").unwrap(), ObserverLawInput::Iso9613);
        assert_eq!(parse("stokes-kirchhoff").unwrap(), ObserverLawInput::StokesKirchhoff);
    }
}
