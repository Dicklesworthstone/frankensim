//! Conservative internal-energy inversion with explicit real-fluid densities.
//!
//! NASA-9 supplies ideal-reference calorics and PR supplies the residual energy
//! and Cv at the SAME composition, temperature and molar volume. The bounded
//! root solver never changes density to select a different pressure root.
//! Callers declare a continuous mechanically/calorically stable homogeneous
//! temperature bracket. Every evaluated state must have p>0, (dp/dv)_T<0 and
//! Cv>0; these local checks do not certify global phase equilibrium or the
//! absence of two-phase coexistence between evaluated states.

/// Constant-mass storage with continuing fixed-volume real-gas heating.
pub mod parcel;
/// Entropy changes and reversible adiabatic density reconstruction.
pub mod entropy;
/// Source-bound partial molar calorics for species-selective real-gas transport.
pub mod partial;

use fs_matdb::QueryPoint;
use fs_qty::{Density, Temperature};
use fs_thermochem::eos::peng_robinson::PengRobinsonDensityCaloricEvaluationV1;
use fs_thermochem::eos::peng_robinson::caloric::PengRobinsonCaloricReferenceV1;

use super::{ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel};
use crate::thermochemical_gas::{coordinates, failure, positive};
use crate::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve, map_root, poll, root, updated_point,
};

/// Complete accepted real-fluid state, with an absolute energy residual.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasEnergySolution {
    state: ResolvedThermochemicalRealGasState,
    residual_j_kg: f64,
    iterations: usize,
}
impl RealGasEnergySolution {
    /// The same resolved state for heat, momentum and acoustic consumers.
    pub const fn state(&self) -> &ResolvedThermochemicalRealGasState { &self.state }
    /// Move the accepted state into an atomic storage update without cloning.
    pub fn into_state(self) -> ResolvedThermochemicalRealGasState { self.state }
    /// Actual selected energy or increment minus the requested value [J/kg].
    pub const fn residual_j_kg(&self) -> f64 { self.residual_j_kg }
    /// Bounded interior iterations; the solver also evaluates both endpoints.
    pub const fn iterations(&self) -> usize { self.iterations }
}

/// Original relative internal energy and its final fixed-density tangent.
///
/// This is a homogeneous caloric query, not an accepted transport/material
/// state. It checks the original EOS and caloric support at both endpoints;
/// consumers must resolve their complete accepted material state separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasCaloricIncrement {
    internal_energy_change_j_kg: f64,
    final_cv_j_kg_k: f64,
}
impl RealGasCaloricIncrement {
    /// Signed original NASA/PR Delta u [J/kg], including density departure.
    pub const fn internal_energy_change_j_kg(self) -> f64 {
        self.internal_energy_change_j_kg
    }
    /// Complete final Cv [J/(kg K)], not ideal-gas Cp-R or a frozen secant.
    /// This is d(Delta u)/dT_final while both endpoint densities remain fixed.
    pub const fn final_cv_j_kg_k(self) -> f64 { self.final_cv_j_kg_k }
}

struct CaloricSample {
    internal_energy_j_kg: f64,
    residual_energy_j_kg: f64,
    cv_j_kg_k: f64,
    pressure_pa: f64,
    pressure_temperature_derivative: f64,
}

impl ThermochemicalPengRobinsonModel {
    fn caloric_at_density(&self, density: f64, t: f64) -> Result<CaloricSample, GasEnergyError> {
        let direct = self.mechanics().evaluate_density_caloric(
            Temperature::new(t), Density::new(density),
        ).map_err(|error| failure("real-gas density mechanics", error))?;
        self.caloric_from_eos(&direct)
    }

    fn caloric_from_eos(&self, direct: &PengRobinsonDensityCaloricEvaluationV1)
        -> Result<CaloricSample, GasEnergyError> {
        let properties = self.source().thermodynamics()
            .evaluate_properties(direct.temperature(), direct.pressure())
            .map_err(|error| failure("real-gas reference calorics", error))?;
        let molar_mass = properties.molar_mass().value();
        let residual_energy_j_kg = direct.residual_internal_energy_j_mol() / molar_mass;
        let internal_energy_j_kg = properties.u().value() / molar_mass + residual_energy_j_kg;
        let cv_j_kg_k = properties.mass_specific_cv().value()
            + direct.residual_cv_j_mol_k() / molar_mass;
        let pressure_temperature_derivative = direct.pressure_derivatives().0;
        if !internal_energy_j_kg.is_finite() || !residual_energy_j_kg.is_finite()
            || !positive(cv_j_kg_k) || !pressure_temperature_derivative.is_finite() {
            return Err(GasEnergyError::Invalid("finite real-gas energy and positive complete Cv required"));
        }
        Ok(CaloricSample { internal_energy_j_kg, residual_energy_j_kg, cv_j_kg_k,
            pressure_pa: direct.pressure().value(), pressure_temperature_derivative })
    }

    fn caloric_reference(&self, initial_density: f64, initial_temperature: f64, final_density: f64)
        -> Result<PengRobinsonCaloricReferenceV1<'_>, GasEnergyError> {
        self.mechanics().prepare_density_caloric_change(Temperature::new(initial_temperature),
            Density::new(initial_density), Density::new(final_density))
            .map_err(|error| failure("real-gas caloric reference", error).into())
    }

    /// Reconstruct temperature and pressure from density and absolute internal
    /// energy in the original NASA/PR reference convention. Kinetic energy must
    /// already have been removed by the consumer. Transport resolves only at
    /// the accepted state; no cubic root solve occurs during Newton iterations.
    pub fn reconstruct_density_energy(
        &self, context: &QueryPoint, density_kg_m3: f64, internal_energy_j_kg: f64,
        controls: GasEnergySolve, cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasEnergySolution, GasEnergyError> {
        self.reconstruct_density_energy_change(context, density_kg_m3, 0.0,
            internal_energy_j_kg, controls, cancelled)
    }

    /// Reconstruct from a fixed energy reference and signed change without
    /// first adding the change to a potentially large formation-energy datum.
    /// For continuing heating, use the temperature-anchored increment method
    /// below, which also cancels the original NASA formation reference.
    pub fn reconstruct_density_energy_change(
        &self, context: &QueryPoint, density_kg_m3: f64,
        reference_internal_energy_j_kg: f64, internal_energy_change_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasEnergySolution, GasEnergyError> {
        controls.validate()?;
        coordinates(context)?;
        if !positive(density_kg_m3) || !reference_internal_energy_j_kg.is_finite()
            || !internal_energy_change_j_kg.is_finite() {
            return Err(GasEnergyError::Invalid("positive density and finite real-gas energy/reference required"));
        }
        let solved = root::solve(controls.minimum_temperature_k, controls.maximum_temperature_k,
            controls.energy_tolerance_j_kg, controls.maximum_iterations, |t| {
                let sample = self.caloric_at_density(density_kg_m3, t)?;
                Ok::<_, GasEnergyError>(root::Sample {
                    residual: (sample.internal_energy_j_kg - reference_internal_energy_j_kg)
                        - internal_energy_change_j_kg,
                    tangent: sample.cv_j_kg_k,
                    pressure_pa: sample.pressure_pa,
                })
            }, &mut cancelled).map_err(map_root)?;
        poll(&mut cancelled)?;
        let point = updated_point(context, solved.temperature_k, solved.sample.pressure_pa)?;
        let state = self.evaluate_density_temperature(&point, density_kg_m3)?;
        let residual = (state.specific_internal_energy_j_kg() - reference_internal_energy_j_kg)
            - internal_energy_change_j_kg;
        finish(state, density_kg_m3, residual, solved.iterations, controls, &mut cancelled)
    }

    /// Original-law specific internal-energy change at fixed density [J/kg].
    /// The ideal-reference term is integrated directly between temperatures;
    /// large formation energies are never subtracted. The PR departure uses
    /// the same density at both endpoints. Requires one declared homogeneous
    /// branch; this does not provide a phase-change energy or latent heat law.
    pub fn specific_internal_energy_change_at_density(
        &self, density_kg_m3: f64, initial_temperature_k: f64, final_temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        self.specific_internal_energy_change_between_densities(density_kg_m3,
            initial_temperature_k, density_kg_m3, final_temperature_k)
    }

    /// Original-law Delta u [J/kg] between two explicit density/temperature
    /// endpoints. The ideal-reference increment depends on temperature; the PR
    /// departure depends on BOTH temperature and density. Replacing the original
    /// departure with one at the final density would erase expansion energy.
    /// Neither heat nor boundary work is inferred from this constitutive query.
    pub fn specific_internal_energy_change_between_densities(
        &self, initial_density_kg_m3: f64, initial_temperature_k: f64,
        final_density_kg_m3: f64, final_temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        self.specific_internal_energy_change_and_cv_between_densities(
            initial_density_kg_m3, initial_temperature_k,
            final_density_kg_m3, final_temperature_k,
        ).map(RealGasCaloricIncrement::internal_energy_change_j_kg)
    }

    /// Original Delta u and final Cv from the same direct density-caloric query.
    /// The original finite NASA/PR energy increment remains independent of
    /// formation datums, including changes of density and very small heating.
    ///
    /// Both endpoints retain positive pressure, mechanical stability and
    /// positive complete Cv on the original caloric source. The caller declares
    /// a continuous homogeneous branch. No PT cubic inventory, component
    /// fugacity or transport/card resolution is performed: this supplies a
    /// nonlinear caloric trial, not a complete accepted material state.
    pub fn specific_internal_energy_change_and_cv_between_densities(
        &self, initial_density_kg_m3: f64, initial_temperature_k: f64,
        final_density_kg_m3: f64, final_temperature_k: f64,
    ) -> Result<RealGasCaloricIncrement, GasEnergyError> {
        let reference = self.caloric_reference(initial_density_kg_m3,
            initial_temperature_k, final_density_kg_m3)?;
        self.caloric_from_eos(reference.initial_state())?;
        let change = reference.evaluate(Temperature::new(final_temperature_k))
            .map_err(|error| failure("real-gas caloric increment", error))?;
        let final_state = self.caloric_from_eos(change.state())?;
        let internal_energy_change_j_kg = self.caloric_increment(
            initial_temperature_k, final_temperature_k,
            change.residual_internal_energy_change_j_mol(),
        )?;
        Ok(RealGasCaloricIncrement {
            internal_energy_change_j_kg,
            final_cv_j_kg_k: final_state.cv_j_kg_k,
        })
    }

    /// Original-law Delta h [J/kg] between explicit density/temperature states.
    /// The original NASA enthalpy increment and finite PR departure include
    /// pressure work once. No absolute formation/departure enthalpy or rounded
    /// endpoint pressure is subtracted. The two original EOS/caloric samples
    /// are shared with the finite difference; endpoint reversal is preserved.
    pub fn specific_enthalpy_change_between_densities(
        &self, initial_density_kg_m3: f64, initial_temperature_k: f64,
        final_density_kg_m3: f64, final_temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        let reference = self.caloric_reference(initial_density_kg_m3,
            initial_temperature_k, final_density_kg_m3)?;
        self.caloric_from_eos(reference.initial_state())?;
        let departure = reference.evaluate(Temperature::new(final_temperature_k))
            .map_err(|error| failure("real-gas enthalpy increment", error))?;
        self.caloric_from_eos(departure.state())?;
        let ideal = self.source().thermodynamics().specific_enthalpy_change(
            Temperature::new(initial_temperature_k), Temperature::new(final_temperature_k),
        ).map_err(|error| failure("real-gas reference enthalpy increment", error))?.value();
        let change = ideal + departure.residual_enthalpy_change_j_mol()
            / self.mechanics().molar_mass().value();
        if !change.is_finite() {
            return Err(GasEnergyError::Invalid("real-gas enthalpy increment is unrepresentable"));
        }
        Ok(change)
    }

    fn caloric_increment(&self, initial: f64, final_t: f64, residual_change_j_mol: f64)
        -> Result<f64, GasEnergyError> {
        let ideal = self.source().thermodynamics().specific_internal_energy_change(
            Temperature::new(initial), Temperature::new(final_t),
        ).map_err(|error| failure("real-gas reference energy increment", error))?.value();
        let change = ideal + residual_change_j_mol / self.mechanics().molar_mass().value();
        if !change.is_finite() {
            return Err(GasEnergyError::Invalid("real-gas energy increment is not representable"));
        }
        Ok(change)
    }

    /// Recover a real-fluid state after a signed fixed-density energy input.
    /// `context` supplies the initial temperature and retained source axes.
    /// Source caloric continuity and the same homogeneous PR branch are caller
    /// requirements on the bracket. The initial and final material sources
    /// are checked, including their pressure domains; no transport extrapolates.
    /// Nonzero heat cannot silently accept the unchanged temperature just
    /// because the input falls below the requested absolute tolerance.
    pub fn reconstruct_density_internal_energy_increment(
        &self, context: &QueryPoint, density_kg_m3: f64, change_j_kg: f64,
        controls: GasEnergySolve, cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasEnergySolution, GasEnergyError> {
        self.reconstruct_density_internal_energy_increment_between_densities(
            context, density_kg_m3, density_kg_m3, change_j_kg, controls, cancelled,
        )
    }

    /// Recover the state at a new density from the original state's Delta u.
    /// `context` supplies the reference temperature and source axes; explicit
    /// initial density selects the original homogeneous branch without an
    /// ambiguous PT root choice. Final density is independently conserved.
    /// Zero Delta u at changed density is a real-gas Joule expansion/compression
    /// endpoint, not an instruction to keep temperature fixed. The enclosing
    /// consumer owns heat, work and the physical path between these endpoints.
    pub fn reconstruct_density_internal_energy_increment_between_densities(
        &self, context: &QueryPoint, initial_density_kg_m3: f64,
        density_kg_m3: f64, change_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasEnergySolution, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let (initial_temperature, _) = coordinates(context)?;
        let initial_t = initial_temperature.value();
        if !positive(initial_density_kg_m3) || !positive(density_kg_m3) || !change_j_kg.is_finite()
            || initial_t < controls.minimum_temperature_k || initial_t > controls.maximum_temperature_k {
            return Err(GasEnergyError::Invalid("real-gas increment needs positive density, finite heat and bracketed reference temperature"));
        }
        let initial = self.evaluate_density_temperature(context, initial_density_kg_m3)?;
        if change_j_kg == 0.0 && density_kg_m3 == initial_density_kg_m3 {
            poll(&mut cancelled)?;
            return Ok(RealGasEnergySolution { state: initial, residual_j_kg: 0.0, iterations: 0 });
        }
        let reference = self.caloric_reference(initial_density_kg_m3, initial_t, density_kg_m3)?;
        let tolerance = if change_j_kg == 0.0 { controls.energy_tolerance_j_kg }
            else { controls.energy_tolerance_j_kg.min(0.5 * change_j_kg.abs()) };
        if !positive(tolerance) {
            return Err(GasEnergyError::Invalid("real-gas increment is below represented solve resolution"));
        }
        let (lo, hi) = if density_kg_m3 != initial_density_kg_m3 {
            (controls.minimum_temperature_k, controls.maximum_temperature_k)
        } else if change_j_kg > 0.0 { (initial_t, controls.maximum_temperature_k) }
            else { (controls.minimum_temperature_k, initial_t) };
        let solved = root::solve(lo, hi, tolerance, controls.maximum_iterations, |t| {
            let change = reference.evaluate(Temperature::new(t))
                .map_err(|error| failure("real-gas caloric increment", error))?;
            let sample = self.caloric_from_eos(change.state())?;
            Ok::<_, GasEnergyError>(root::Sample {
                residual: self.caloric_increment(initial_t, t,
                    change.residual_internal_energy_change_j_mol())? - change_j_kg,
                tangent: sample.cv_j_kg_k,
                pressure_pa: sample.pressure_pa,
            })
        }, &mut cancelled).map_err(map_root)?;
        poll(&mut cancelled)?;
        let point = updated_point(context, solved.temperature_k, solved.sample.pressure_pa)?;
        let state = self.evaluate_density_temperature(&point, density_kg_m3)?;
        let departure = reference.evaluate(Temperature::new(state.temperature_k()))
            .map_err(|error| failure("real-gas accepted caloric increment", error))?;
        let change = self.caloric_increment(initial_t, state.temperature_k(),
            departure.residual_internal_energy_change_j_mol())?;
        let residual = change - change_j_kg;
        if change_j_kg != 0.0 && (change == 0.0
            || change.is_sign_positive() != change_j_kg.is_sign_positive()
            || (density_kg_m3 == initial_density_kg_m3 && state.temperature_k() == initial_t)) {
            return Err(GasEnergyError::NotConverged { iterations: solved.iterations, residual_j_kg: residual });
        }
        finish(state, density_kg_m3, residual, solved.iterations, controls, &mut cancelled)
    }
}

fn finish(state: ResolvedThermochemicalRealGasState, density: f64, residual: f64,
    iterations: usize, controls: GasEnergySolve, cancelled: &mut impl FnMut() -> bool)
    -> Result<RealGasEnergySolution, GasEnergyError> {
    if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j_kg {
        return Err(GasEnergyError::NotConverged { iterations, residual_j_kg: residual });
    }
    if state.density_kg_m3().to_bits() != density.to_bits() {
        return Err(GasEnergyError::Invalid("real-gas reconstruction changed conserved density"));
    }
    poll(cancelled)?;
    Ok(RealGasEnergySolution { state, residual_j_kg: residual, iterations })
}

#[cfg(test)]
mod tests;
