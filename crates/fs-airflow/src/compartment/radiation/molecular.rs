//! TNF/RADCAL molecular Planck means with original ideal-gas composition.
//!
//! Source: <https://tnfworkshop.org/radiation/>, corrected 7 June 2000,
//! discussion updated 23 January 2003; Grosshandler, NIST TN 1402 (1993).
//! The published fits are retained for H2O, CO2, CO and CH4, with pressures
//! converted from Pa to atm and temperature restricted to 300..2500 K.
//! H2O/CO2 coefficients also appear in Cantera v3.2.0 Flow1D.cpp.
//!
//! These are local-temperature emission means, originally used in optically
//! thin flame calculations. Using the same mean for absorption of another
//! spectrum or a finite optical path is an explicit gray approximation, not
//! spectral transfer, a path-averaged fit, or certified optical accuracy. The
//! archived TNF model is a simple reference closure, not current spectroscopy.
//! Unrepresented gases require a caller-declared transparency assumption.

use super::{AbsorptionCoefficient, LayerAbsorption, RadiationError};
use fs_gas::multicomponent::MixtureGas;
use fs_gas::multicomponent::compartment::{Layer, TwoZone};
use fs_qty::{Pressure, Temperature};

/// Canonical contribution order in every report.
pub const MOLECULAR_SPECIES: [&str; 4] = ["H2O", "CO2", "CO", "CH4"];
const ATM_PA: f64 = 101_325.0;

/// Source-bound mapping from original gas species to four published fits.
/// Construction allocates; evaluating admitted layer coefficients does not.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularPlanck {
    gas: MixtureGas,
    indices: [Option<usize>; 4],
    other_species_assumption: String,
}

/// Actual local-state contributions; no propagated uncertainty is implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularPlanckReport<'a> {
    /// Original gas temperature, also the Planck weighting temperature.
    pub temperature: Temperature,
    /// Original common room pressure.
    pub pressure: Pressure,
    /// H2O, CO2, CO, CH4 partial pressures, in Pa.
    pub partial_pressures: [Pressure; 4],
    /// Corresponding additive absorption contributions, in 1/m.
    pub contributions: [AbsorptionCoefficient; 4],
    /// Sum of the four nonnegative original contributions, in 1/m.
    pub absorption: AbsorptionCoefficient,
    /// Explicit scope for every species without a supplied fit.
    pub other_species_assumption: &'a str,
}

impl MolecularPlanckReport<'_> {
    /// Explicitly adopt this emission mean as a gray enclosure coefficient.
    /// The assumption must explain use over the actual incident spectra/paths.
    pub fn gray_absorption(
        &self,
        gray_assumption: &str,
    ) -> Result<LayerAbsorption, RadiationError> {
        if gray_assumption.trim().is_empty() {
            return Err(RadiationError::Invalid(
                "molecular Planck mean requires an explicit gray-use assumption",
            ));
        }
        LayerAbsorption::new(
            self.absorption,
            format!(
                "TNF/RADCAL corrected molecular Planck fits; https://tnfworkshop.org/radiation/; T={} K; p={} Pa; partial_H2O_CO2_CO_CH4_Pa={:?}; other_species={}; gray_use={}",
                self.temperature.value(),
                self.pressure.value(),
                self.partial_pressures.map(|p| p.value()),
                self.other_species_assumption,
                gray_assumption,
            ),
        )
    }
}

impl MolecularPlanck {
    /// Bind exact canonical names and verify their elemental/charge identities.
    /// All unrepresented species are transparent only under the stated assumption.
    pub fn new(
        gas: &MixtureGas,
        other_species_assumption: impl Into<String>,
    ) -> Result<Self, RadiationError> {
        let assumption = other_species_assumption.into();
        if assumption.trim().is_empty() {
            return Err(RadiationError::Invalid(
                "unrepresented molecular species need an explicit optical assumption",
            ));
        }
        let bank = gas.source();
        let elements = bank.elemental_matrix();
        let mut indices = [None; 4];
        // Integer C/H/O counts follow the identities named in the source fits.
        let atoms = [[0, 2, 1], [1, 0, 2], [1, 0, 1], [1, 4, 0]];
        for (i, model) in bank.models().iter().enumerate() {
            let Some(k) = MOLECULAR_SPECIES
                .iter()
                .position(|name| *name == model.species().as_str())
            else {
                continue;
            };
            for (e, element) in elements.elements().iter().enumerate() {
                let expected = match element.as_str() {
                    "C" => atoms[k][0],
                    "H" => atoms[k][1],
                    "O" => atoms[k][2],
                    _ => 0,
                };
                if elements.get(e, i) != Some(expected) {
                    return Err(RadiationError::Invalid(
                        "molecular fit name disagrees with original elemental identity",
                    ));
                }
            }
            for (symbol, count) in ["C", "H", "O"].into_iter().zip(atoms[k]) {
                if count > 0 && !elements.elements().iter().any(|e| e.as_str() == symbol) {
                    return Err(RadiationError::Invalid(
                        "molecular fit identity is missing an original element",
                    ));
                }
            }
            if bank.charge_vector().get(i) != Some(0) {
                return Err(RadiationError::Invalid(
                    "neutral molecular absorption cannot bind an ion",
                ));
            }
            indices[k] = Some(i);
        }
        Ok(Self {
            gas: gas.clone(),
            indices,
            other_species_assumption: assumption,
        })
    }

    /// Evaluate actual layers with exact absent-layer handling and bounded polling.
    /// No absorption fit is queried for a species with exactly zero inventory.
    pub fn evaluate_room(
        &self,
        room: &TwoZone,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[Option<MolecularPlanckReport<'_>>; 2], RadiationError> {
        if cancelled() {
            return Err(RadiationError::Cancelled);
        }
        if room.gas() != &self.gas {
            return Err(RadiationError::Invalid(
                "molecular opacity and room original gas sources differ",
            ));
        }
        let mut reports = [None; 2];
        for layer in [Layer::Lower, Layer::Upper] {
            if cancelled() {
                return Err(RadiationError::Cancelled);
            }
            let state = room.layer(layer);
            if state.mass_kg() == 0.0 {
                continue;
            }
            let amounts: f64 = state
                .species_masses()
                .iter()
                .zip(self.gas.source().models())
                .map(|(mass, model)| mass / model.molar_mass().value())
                .sum();
            if !amounts.is_finite() || amounts <= 0.0 {
                return Err(RadiationError::Invalid(
                    "unrepresentable molecular inventory",
                ));
            }
            let mut partial = [Pressure::new(0.0); 4];
            let mut contributions = [AbsorptionCoefficient::new(0.0); 4];
            for (k, index) in self.indices.iter().enumerate() {
                let Some(i) = index else { continue };
                let m = state.species_masses()[*i];
                if m == 0.0 {
                    continue;
                }
                partial[k] = Pressure::new(
                    room.pressure().value()
                        * (m / self.gas.source().models()[*i].molar_mass().value() / amounts),
                );
                let kappa = fit(k, state.temperature().value())? * (partial[k].value() / ATM_PA);
                if !kappa.is_finite() || kappa < 0.0 || !partial[k].value().is_finite() {
                    return Err(RadiationError::Invalid("unrepresentable molecular opacity"));
                }
                contributions[k] = AbsorptionCoefficient::new(kappa);
            }
            let total: f64 = contributions.iter().map(|v| v.value()).sum();
            if !total.is_finite() {
                return Err(RadiationError::Invalid(
                    "unrepresentable molecular absorption sum",
                ));
            }
            reports[layer.index()] = Some(MolecularPlanckReport {
                temperature: state.temperature(),
                pressure: room.pressure(),
                partial_pressures: partial,
                contributions,
                absorption: AbsorptionCoefficient::new(total),
                other_species_assumption: &self.other_species_assumption,
            });
        }
        Ok(reports)
    }
}

// Published polynomial values are in m^-1 atm^-1, not m^-1 Pa^-1.
fn fit(species: usize, temperature_k: f64) -> Result<f64, RadiationError> {
    if !temperature_k.is_finite() || !(300.0..=2500.0).contains(&temperature_k) {
        return Err(RadiationError::MolecularTemperature { temperature_k });
    }
    let t = temperature_k;
    let result = match species {
        0 => polynomial(
            1000.0 / t,
            &[-0.23093, -1.12390, 9.41530, -2.99880, 0.51382, -1.86840e-5],
        ),
        1 => polynomial(
            1000.0 / t,
            &[18.741, -121.310, 273.500, -194.050, 56.310, -5.8169],
        ),
        2 if t <= 750.0 => polynomial(t, &[4.7869, -0.06953, 2.95775e-4, -4.25732e-7, 2.02894e-10]),
        2 => polynomial(t, &[10.09, -0.01183, 4.7753e-6, -5.87209e-10, -2.5334e-14]),
        3 => polynomial(t, &[6.6334, -0.0035686, 1.6682e-8, 2.5611e-10, -2.6558e-14]),
        _ => return Err(RadiationError::Invalid("unknown molecular Planck axis")),
    };
    if !result.is_finite() || result < 0.0 {
        return Err(RadiationError::Invalid(
            "source molecular Planck fit is negative or nonfinite",
        ));
    }
    Ok(result)
}

fn polynomial(x: f64, coefficients: &[f64]) -> f64 {
    coefficients
        .iter()
        .rev()
        .fold(0.0, |value, c| value * x + c)
}

#[cfg(test)]
#[path = "molecular/tests.rs"]
mod tests;
