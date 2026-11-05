//! Bind a bounded D–T burn to the original ideal-plasma material state.
//!
//! The bank contains exactly D+, T+, He4++ and one structural electron.
//! Isotopes are explicitly mapped, never inferred from species names. Their
//! original monatomic NASA-9 sources supply sensible energy, pressure and heat
//! capacity; nuclear rest energy and the neutron inventory belong to the
//! fusion owner. A source formation offset would count an unspecified part of
//! nuclear energy again, so this lane requires zero formation constants.
//!
//! This is a source binding and classical ideal-plasma admission, not a nuclear
//! equilibrium solver, Maxwellian relaxation calculation or fusion integrator.
//! The density/temperature screens use the same thermal wavelength, Coulomb
//! coupling and nonrelativistic scales as the existing `coulomb` owner, without
//! imposing its separate classical-scattering approximation. Maxwellian
//! distributions and rapid thermalization remain explicit consumer assumptions.

use super::{BOLTZMANN_J_PER_K, ELECTRON_MASS_KG, NuclearSpeciesV1, SPEED_OF_LIGHT_M_S};
use crate::electron_collision::AVOGADRO_PER_MOL;
use crate::equilibrium::composition::{FrozenIonizedGasStateV1, IdealGasEquilibriumV1};
use crate::{SpeciesId, ThermochemErrorV1};
use fs_math::det::{exp, ln};
use fs_qty::{Amount, Temperature, Volume};

const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;
const VACUUM_PERMITTIVITY_F_M: f64 = 8.854_187_818_8e-12;
const PLANCK_J_S: f64 = 6.626_070_15e-34;
// Numerical identity tolerances, not uncertainties of the physical model.
const RELATIVE_MASS_TOLERANCE: f64 = 64.0 * f64::EPSILON;
const RELATIVE_CHARGE_TOLERANCE: f64 = 64.0 * f64::EPSILON;

/// Refusal to bind sources or admit a classical D–T plasma point.
#[derive(Clone, Debug, PartialEq)]
pub enum DtFusionPlasmaErrorV1 {
    /// Unsupported species, source reference, population shape or control.
    Invalid(&'static str),
    /// The original standard-state model refused its actual temperature.
    Source {
        /// Original species identity.
        species: SpeciesId,
        /// Unmodified source refusal.
        source: ThermochemErrorV1,
    },
    /// An explicit physical-domain limit was exceeded.
    Domain(&'static str),
    /// An intermediate physical scale cannot be represented.
    Arithmetic(&'static str),
    /// Cancellation before publication of a complete result.
    Cancelled,
}
impl core::fmt::Display for DtFusionPlasmaErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "D-T plasma source: {self:?}")
    }
}
impl std::error::Error for DtFusionPlasmaErrorV1 {}

/// Caller-selected validity screens for classical, nonrelativistic ideal pools.
/// These limits are not measured accuracy or numerical error certificates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DtFusionPlasmaLimitsV1 {
    /// Largest species z²e²/(4*pi*eps0*a*kB*T), with
    /// a=[3/(4*pi*n_total)]^(1/3) from all charged particles.
    pub maximum_coulomb_coupling: f64,
    /// Largest n_s*[h/sqrt(2*pi*m_s*kB*T_s)]³; no spin factor is inferred.
    pub maximum_thermal_phase_space_density: f64,
    /// Largest kB*T_s/(m_s*c²) over the populated original species.
    pub maximum_thermal_energy_rest_energy_ratio: f64,
}
impl DtFusionPlasmaLimitsV1 {
    /// Every maximum must be positive, finite and strictly below one.
    pub fn validate(self) -> Result<(), DtFusionPlasmaErrorV1> {
        if [
            self.maximum_coulomb_coupling,
            self.maximum_thermal_phase_space_density,
            self.maximum_thermal_energy_rest_energy_ratio,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0 || *x >= 1.0)
        {
            return Err(DtFusionPlasmaErrorV1::Invalid(
                "physical plasma maxima must lie strictly inside (0,1)",
            ));
        }
        Ok(())
    }
}

/// Actual largest ideal-plasma scales, not uncertainty estimates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DtFusionPlasmaScalesV1 {
    /// Largest admitted species Coulomb coupling parameter.
    pub maximum_coulomb_coupling: f64,
    /// Largest admitted species thermal phase-space density.
    pub maximum_thermal_phase_space_density: f64,
    /// Largest admitted thermal/rest-energy ratio.
    pub maximum_thermal_energy_rest_energy_ratio: f64,
}

/// Immutable source bank and explicit canonical D/T/He4/electron associations.
#[derive(Clone, Debug, PartialEq)]
pub struct DtFusionPlasmaV1 {
    thermodynamics: IdealGasEquilibriumV1,
    indices: [usize; 3],
    electron_index: usize,
    log_masses_kg: [f64; 4],
}
impl DtFusionPlasmaV1 {
    /// Bind the original four sources without altering their cards or references.
    /// Nuclear/electron masses must match this fusion model's nominal masses to
    /// 64 epsilon relative tolerance; this prevents an independently inconsistent
    /// rest-mass ledger. All sources require fixed Cp/R=5/2, zero a7, fixed a8
    /// across their regions, and no transformed caloric reference.
    pub fn new(
        thermodynamics: IdealGasEquilibriumV1,
        deuterium: SpeciesId,
        tritium: SpeciesId,
        helium4: SpeciesId,
    ) -> Result<Self, DtFusionPlasmaErrorV1> {
        use DtFusionPlasmaErrorV1 as E;
        let elements = thermodynamics.elemental_matrix();
        if thermodynamics.models().len() != 4
            || elements.row_count() != 2
            || elements.elements()[0].as_str() != "H"
            || elements.elements()[1].as_str() != "He"
        {
            return Err(E::Invalid("require exactly D+, T+, He4++, electron on H/He axes"));
        }
        let mut indices = [0; 3];
        for (slot, id) in [deuterium, tritium, helium4].iter().enumerate() {
            indices[slot] = elements
                .species()
                .binary_search(id)
                .map_err(|_| E::Invalid("fusion species is absent from the original bank"))?;
        }
        if indices[0] == indices[1] || indices[0] == indices[2] || indices[1] == indices[2] {
            return Err(E::Invalid("D, T and He4 require distinct source species"));
        }
        let electron_index = (0..4)
            .find(|i| !indices.contains(i))
            .expect("three distinct indices in a four-species bank");
        let ordered = [indices[0], indices[1], indices[2], electron_index];
        let expected_elements = [[1, 0], [1, 0], [0, 1], [0, 0]];
        let expected_charges = [1, 1, 2, -1];
        let expected_masses = [
            NuclearSpeciesV1::Deuterium.mass_kg(),
            NuclearSpeciesV1::Tritium.mass_kg(),
            NuclearSpeciesV1::Helium4.mass_kg(),
            ELECTRON_MASS_KG,
        ];
        let mut log_masses_kg = [0.0; 4];
        for (slot, &index) in ordered.iter().enumerate() {
            if thermodynamics.charge_vector().get(index) != Some(expected_charges[slot])
                || (0..2).any(|row| elements.get(row, index) != Some(expected_elements[slot][row]))
            {
                return Err(E::Invalid("fusion mapping requires H+, H+, He++, structural electron"));
            }
            let model = &thermodynamics.models()[index];
            let mass = model.molar_mass().value() / AVOGADRO_PER_MOL;
            if ((mass - expected_masses[slot]) / expected_masses[slot]).abs()
                > RELATIVE_MASS_TOLERANCE
            {
                return Err(E::Invalid("original source mass disagrees with the mapped bare nucleus/electron"));
            }
            log_masses_kg[index] = ln(mass);
            if model.caloric_reference().is_some() {
                return Err(E::Invalid("fusion requires original untransformed sensible-energy references"));
            }
            let entropy_constant = model.regions()[0].coefficients()[8];
            for region in model.regions() {
                let c = region.coefficients();
                if c[2] != 2.5
                    || [0, 1, 3, 4, 5, 6, 7].iter().any(|&i| c[i] != 0.0)
                    || c[8] != entropy_constant
                {
                    return Err(E::Invalid("fusion requires fixed monatomic Cp/R=5/2, zero a7 and fixed a8"));
                }
            }
        }
        Ok(Self {
            thermodynamics,
            indices,
            electron_index,
            log_masses_kg,
        })
    }

    /// Original common source bank, including its unchanged cards and provenance.
    pub const fn thermodynamics(&self) -> &IdealGasEquilibriumV1 {
        &self.thermodynamics
    }
    /// Canonical bank indices in the physical order [D, T, He4].
    pub const fn indices(&self) -> [usize; 3] {
        self.indices
    }
    /// Canonical bank index of the original structural electron.
    pub const fn electron_index(&self) -> usize {
        self.electron_index
    }

    /// Check exact original-bank binding and source support at actual Th/Te.
    /// Potential products must be supported even when their amount is zero.
    pub fn admit_state(
        &self,
        state: &FrozenIonizedGasStateV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), DtFusionPlasmaErrorV1> {
        poll(&mut cancelled)?;
        if state.thermodynamics() != &self.thermodynamics
            || state.electron_index() != Some(self.electron_index)
        {
            return Err(DtFusionPlasmaErrorV1::Invalid("state belongs to a different original source bank"));
        }
        self.admit_sources(state.heavy_temperature(), state.electron_temperature(), &mut cancelled)
    }

    /// Admit an original state with explicit classical ideal-plasma limits.
    pub fn admit_state_with_limits(
        &self,
        state: &FrozenIonizedGasStateV1,
        limits: DtFusionPlasmaLimitsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DtFusionPlasmaScalesV1, DtFusionPlasmaErrorV1> {
        self.admit_state(state, &mut cancelled)?;
        self.admit_point(
            state.heavy_temperature(), state.electron_temperature(),
            state.amounts(), state.volume(), limits, cancelled,
        )
    }

    /// Allocation-free source and physical admission for a nonlinear trial.
    /// Amounts use the original bank's canonical order. Quasineutral charge
    /// must close to 64 epsilon relative to actual absolute carrier charge;
    /// no population, temperature, or charge correction is applied.
    pub fn admit_point(
        &self,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        amounts: &[Amount],
        volume: Volume,
        limits: DtFusionPlasmaLimitsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DtFusionPlasmaScalesV1, DtFusionPlasmaErrorV1> {
        use DtFusionPlasmaErrorV1 as E;
        poll(&mut cancelled)?;
        limits.validate()?;
        if amounts.len() != 4 || !positive(volume.value())
            || amounts.iter().any(|n| !n.value().is_finite() || n.value() < 0.0)
        {
            return Err(E::Invalid("four nonnegative finite amounts and positive finite volume required"));
        }
        let scale = amounts.iter().map(|n| n.value()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            return Err(E::Invalid("a positive quasineutral plasma inventory is required"));
        }
        let [d, t, he] = self.indices;
        let positive_charge = amounts[d].value() / scale + amounts[t].value() / scale
            + 2.0 * (amounts[he].value() / scale);
        let electron_charge = amounts[self.electron_index].value() / scale;
        if ((positive_charge - electron_charge) / (positive_charge + electron_charge)).abs()
            > RELATIVE_CHARGE_TOLERANCE
        {
            return Err(E::Invalid("D-T trial charge must match the original quasineutral electron inventory"));
        }
        self.admit_sources(heavy_temperature, electron_temperature, &mut cancelled)?;
        let log_temperatures = [ln(heavy_temperature.value()), ln(electron_temperature.value())];
        let mut log_densities = [f64::NEG_INFINITY; 4];
        let log_density_factor = ln(AVOGADRO_PER_MOL) - ln(volume.value());
        for (i, amount) in amounts.iter().enumerate() {
            if amount.value() > 0.0 {
                log_densities[i] = ln(amount.value()) + log_density_factor;
            }
        }
        let largest_log_density = log_densities.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let scaled_density: f64 = log_densities.iter().map(|n| exp(*n - largest_log_density)).sum();
        let log_total_density = largest_log_density + ln(scaled_density);
        let pi = core::f64::consts::PI;
        let log_coupling_factor = 2.0 * ln(ELEMENTARY_CHARGE_C)
            - ln(4.0 * pi * VACUUM_PERMITTIVITY_F_M * BOLTZMANN_J_PER_K)
            + (ln(4.0 * pi / 3.0) + log_total_density) / 3.0;
        let log_phase_factor = 3.0 * ln(PLANCK_J_S)
            - 1.5 * ln(2.0 * pi * BOLTZMANN_J_PER_K);
        let log_thermal_factor = ln(BOLTZMANN_J_PER_K) - 2.0 * ln(SPEED_OF_LIGHT_M_S);
        let mut maxima = [f64::NEG_INFINITY; 3];
        for (i, &log_density) in log_densities.iter().enumerate() {
            poll(&mut cancelled)?;
            if amounts[i].value() == 0.0 {
                continue;
            }
            let log_t = log_temperatures[usize::from(i == self.electron_index)];
            let log_z_squared = if i == he { ln(4.0) } else { 0.0 };
            maxima[0] = maxima[0].max(log_coupling_factor + log_z_squared - log_t);
            maxima[1] = maxima[1].max(log_density + log_phase_factor - 1.5 * (self.log_masses_kg[i] + log_t));
            maxima[2] = maxima[2].max(log_thermal_factor + log_t - self.log_masses_kg[i]);
        }
        for ((actual, limit), name) in maxima.iter().zip([
            limits.maximum_coulomb_coupling,
            limits.maximum_thermal_phase_space_density,
            limits.maximum_thermal_energy_rest_energy_ratio,
        ]).zip(["Coulomb coupling", "thermal phase-space density", "nonrelativistic thermal energy"]) {
            if !actual.is_finite() {
                return Err(E::Arithmetic(name));
            }
            if *actual > ln(limit) {
                return Err(E::Domain(name));
            }
        }
        poll(&mut cancelled)?;
        Ok(DtFusionPlasmaScalesV1 {
            maximum_coulomb_coupling: exp(maxima[0]),
            maximum_thermal_phase_space_density: exp(maxima[1]),
            maximum_thermal_energy_rest_energy_ratio: exp(maxima[2]),
        })
    }

    /// Admit a contiguous original source interval for each pool, including
    /// unpopulated possible products. Identical fixed coefficients guarantee
    /// continuity across adjacent original regions; gaps remain forbidden.
    pub fn admit_temperature_interval(
        &self,
        minimum_heavy: Temperature,
        maximum_heavy: Temperature,
        minimum_electron: Temperature,
        maximum_electron: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), DtFusionPlasmaErrorV1> {
        use DtFusionPlasmaErrorV1 as E;
        for (index, model) in self.thermodynamics.models().iter().enumerate() {
            poll(&mut cancelled)?;
            let (lower, upper) = if index == self.electron_index {
                (minimum_electron, maximum_electron)
            } else {
                (minimum_heavy, maximum_heavy)
            };
            if !positive(lower.value()) || !positive(upper.value()) || lower > upper {
                return Err(E::Invalid("positive ordered source temperature interval required"));
            }
            self.admit_species_temperature(index, lower)?;
            self.admit_species_temperature(index, upper)?;
            let first = model.region_for_temperature(lower.value()).expect("admitted lower endpoint").0;
            let last = model.region_for_temperature(upper.value()).expect("admitted upper endpoint").0;
            if model.regions()[first..=last].windows(2).any(|pair| {
                pair[0].temperature_max() != pair[1].temperature_min()
            }) {
                return Err(E::Invalid("source temperature interval crosses an original gap"));
            }
        }
        poll(&mut cancelled)
    }

    fn admit_sources(
        &self, heavy: Temperature, electron: Temperature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DtFusionPlasmaErrorV1> {
        for index in 0..4 {
            poll(cancelled)?;
            self.admit_species_temperature(index, if index == self.electron_index { electron } else { heavy })?;
        }
        poll(cancelled)
    }
    fn admit_species_temperature(&self, index: usize, temperature: Temperature) -> Result<(), DtFusionPlasmaErrorV1> {
        let model = &self.thermodynamics.models()[index];
        model.evaluate_properties(temperature).map_err(|source| DtFusionPlasmaErrorV1::Source {
            species: model.species().clone(), source,
        })?;
        Ok(())
    }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DtFusionPlasmaErrorV1> {
    if cancelled() { Err(DtFusionPlasmaErrorV1::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;
