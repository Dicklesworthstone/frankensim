//! Temperature-anchored caloric inversion for one unchanged gas composition.
//!
//! Trial residuals integrate original Cp or Cv. A prepared branch borrows its
//! immutable source model and amortizes continuous-domain/capacity admission
//! across repeated solves. Absolute formation references remain in full states.

use super::{
    GasEnergyError, GasEnergySolution, GasEnergySolve, ResolvedThermochemicalGasState,
    ThermochemicalGasModel, enthalpy, map_root, poll, positive, root, updated_point,
};
use crate::thermochemical_gas::{coordinates, failure};
use fs_matdb::QueryPoint;
use fs_qty::Temperature;

/// A complete constant-pressure state and its directly evaluated enthalpy change.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEnthalpyIncrementSolution {
    state: ResolvedThermochemicalGasState,
    enthalpy_change_j_kg: f64,
    residual_j_kg: f64,
    iterations: usize,
}
impl GasEnthalpyIncrementSolution {
    /// Complete source-resolved thermal, transport and acoustic state.
    pub fn state(&self) -> &ResolvedThermochemicalGasState { &self.state }
    /// Original-law Delta h between context and accepted temperature [J/kg].
    pub fn enthalpy_change_j_kg(&self) -> f64 { self.enthalpy_change_j_kg }
    /// Original-law Delta h minus the supplied enthalpy increment [J/kg].
    pub fn residual_j_kg(&self) -> f64 { self.residual_j_kg }
    /// Interior iterations of the shared safeguarded root solver.
    pub fn iterations(&self) -> usize { self.iterations }
    /// Move the accepted source-resolved state without rebuilding its receipts.
    pub fn into_state(self) -> ResolvedThermochemicalGasState { self.state }
}

/// One stable, continuous NASA-9 branch of an immutable frozen-composition gas.
///
/// Preparation runs the existing bounded source-coefficient analysis once per
/// active species. Subsequent solves reuse those bounds; no source model/card
/// cloning or interval subdivision occurs in their nonlinear iteration loops.
/// The borrow binds the exact model, composition and caloric policy, so a
/// prepared branch cannot silently be rebound to different source data.
/// Transport is still resolved at each accepted endpoint, never cached across
/// temperature/pressure changes. These are source-admission bounds, not a claim
/// that the ideal-gas model is experimentally accurate throughout the interval.
#[derive(Clone, Copy, Debug)]
pub struct GasCaloricBranch<'a> {
    model: &'a ThermochemicalGasModel,
    controls: GasEnergySolve,
    minimum_cp: f64,
    minimum_cv: f64,
}

#[derive(Clone, Copy)]
enum Constraint { Density(f64), Pressure }
struct SolvedIncrement {
    state: ResolvedThermochemicalGasState,
    change: f64,
    residual: f64,
    iterations: usize,
}

impl ThermochemicalGasModel {
    /// Prepare reusable continuous-source Cp/Cv bounds for repeated caloric solves.
    /// Raw multi-region sources require explicit caloric reconciliation first.
    /// Polls cancellation between bounded species admissions and before return.
    pub fn prepare_caloric_branch(
        &self, controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasCaloricBranch<'_>, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        let mut minimum_cp = f64::INFINITY;
        let mut minimum_cv = f64::INFINITY;
        for species in self.thermodynamics().components() {
            poll(&mut cancelled)?;
            let cp = enthalpy::minimum_cp(species, controls)?;
            let r = (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                / species.molar_mass().value()).next_up();
            let cv = (cp - r).next_down();
            minimum_cp = minimum_cp.min(cp);
            minimum_cv = minimum_cv.min(cv);
        }
        if !positive(minimum_cp) || !positive(minimum_cv) {
            return Err(GasEnergyError::Invalid("positive represented caloric capacity bounds required"));
        }
        poll(&mut cancelled)?;
        Ok(GasCaloricBranch { model: self, controls, minimum_cp, minimum_cv })
    }

    /// Reconstruct at the specified density from a temperature-anchored Delta u
    /// [J/kg]. Context temperature is the reference. Composition remains frozen.
    /// Use `prepare_caloric_branch` to amortize source admission across many solves.
    pub fn reconstruct_density_internal_energy_increment(
        &self, context: &QueryPoint, density_kg_m3: f64, change_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnergySolution, GasEnergyError> {
        self.prepare_caloric_branch(controls, &mut cancelled)?
            .reconstruct_density_internal_energy_increment(context, density_kg_m3, change_j_kg, cancelled)
    }

    /// Reconstruct at context pressure from original-law Delta h [J/kg], with
    /// context temperature as the explicit reference. Integrates Cp(T), not a
    /// frozen Cp times Delta T or subtraction of absolute formation enthalpies.
    /// The bracket contains the reference and a stable continuous source branch.
    pub fn reconstruct_pressure_enthalpy_increment(
        &self, context: &QueryPoint, change_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnthalpyIncrementSolution, GasEnergyError> {
        self.prepare_caloric_branch(controls, &mut cancelled)?
            .reconstruct_pressure_enthalpy_increment(context, change_j_kg, cancelled)
    }
}

impl GasCaloricBranch<'_> {
    /// Exact immutable source/composition binding.
    pub fn model(&self) -> &ThermochemicalGasModel { self.model }
    /// Retained physical tolerance, temperature bracket and work budget.
    pub fn controls(&self) -> GasEnergySolve { self.controls }
    /// Conservative lower bound on mass-specific Cp over the whole branch.
    pub fn minimum_cp_j_kg_k(&self) -> f64 { self.minimum_cp }
    /// Conservative lower bound on mass-specific Cv over the whole branch.
    pub fn minimum_cv_j_kg_k(&self) -> f64 { self.minimum_cv }

    /// Reuse this admitted branch for a frozen-composition density/Delta-u solve.
    /// The reference temperature must lie in the prepared bracket. The returned
    /// state's residual is the direct original-law increment minus the request.
    pub fn reconstruct_density_internal_energy_increment(
        &self, context: &QueryPoint, density_kg_m3: f64, change_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnergySolution, GasEnergyError> {
        let solved = solve(self, context, Constraint::Density(density_kg_m3), change_j_kg, &mut cancelled)?;
        Ok(GasEnergySolution { state: solved.state, residual_j_kg: solved.residual,
            iterations: solved.iterations })
    }

    /// Reuse this admitted branch for a frozen-composition pressure/Delta-h solve.
    /// Original source conditions survive; full transport is queried only at the
    /// accepted endpoint. No mutation or partial state is published on refusal.
    pub fn reconstruct_pressure_enthalpy_increment(
        &self, context: &QueryPoint, change_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnthalpyIncrementSolution, GasEnergyError> {
        let solved = solve(self, context, Constraint::Pressure, change_j_kg, &mut cancelled)?;
        Ok(GasEnthalpyIncrementSolution { state: solved.state,
            enthalpy_change_j_kg: solved.change, residual_j_kg: solved.residual,
            iterations: solved.iterations })
    }
}

fn solve(
    branch: &GasCaloricBranch<'_>, context: &QueryPoint, constraint: Constraint,
    requested: f64, cancelled: &mut impl FnMut() -> bool,
) -> Result<SolvedIncrement, GasEnergyError> {
    let (model, controls) = (branch.model, branch.controls);
    poll(cancelled)?;
    let (initial, pressure) = coordinates(context)?;
    if !requested.is_finite() || initial.value() < controls.minimum_temperature_k
        || initial.value() > controls.maximum_temperature_k
    {
        return Err(GasEnergyError::Invalid(
            "finite caloric increment and a bracket containing the reference temperature required",
        ));
    }
    if let Constraint::Density(density) = constraint {
        if !positive(density) {
            return Err(GasEnergyError::Invalid("finite positive conserved density required"));
        }
    }
    let minimum_capacity = match constraint {
        Constraint::Pressure => branch.minimum_cp,
        Constraint::Density(_) => branch.minimum_cv,
    };
    // Source bounds narrow tiny-heat solves without guessing constant capacities.
    let span = (requested.abs() / minimum_capacity).next_up();
    let (lo, hi) = if requested > 0.0 {
        (initial.value(), controls.maximum_temperature_k.min((initial.value() + span).next_up()))
    } else {
        (controls.minimum_temperature_k.max((initial.value() - span).next_down()), initial.value())
    };
    let change_at = |t| -> Result<f64, GasEnergyError> {
        let end = Temperature::new(t);
        let value = match constraint {
            Constraint::Density(_) => model.thermodynamics().specific_internal_energy_change(initial, end),
            Constraint::Pressure => model.thermodynamics().specific_enthalpy_change(initial, end),
        }.map_err(|error| failure("caloric increment", error))?;
        Ok(value.value())
    };
    let evaluate = |t| -> Result<root::Sample, GasEnergyError> {
        let properties = model.thermodynamics().evaluate_properties(Temperature::new(t), pressure)
            .map_err(|error| failure("caloric increment tangent", error))?;
        let (tangent, pressure_pa) = match constraint {
            Constraint::Density(density) => (properties.mass_specific_cv().value(),
                density * properties.mass_specific_gas_constant().value() * t),
            Constraint::Pressure => (properties.mass_specific_cp().value(), pressure.value()),
        };
        Ok(root::Sample { residual: change_at(t)? - requested, tangent, pressure_pa })
    };
    let solved = if requested == 0.0 {
        poll(cancelled)?;
        let sample = evaluate(initial.value())?;
        if !positive(sample.tangent) || !positive(sample.pressure_pa) {
            return Err(GasEnergyError::Invalid("zero-increment EOS state is unrepresentable"));
        }
        root::Solution { temperature_k: initial.value(), sample, iterations: 0 }
    } else {
        // Never accept the unchanged reference merely because a nonzero input
        // is smaller than the user's absolute tolerance. This only tightens the
        // solve; the original physical acceptance target remains unchanged.
        let tolerance = controls.energy_tolerance_j_kg.min(0.5 * requested.abs());
        if !positive(tolerance) {
            return Err(GasEnergyError::Invalid("nonzero caloric increment is below represented solve resolution"));
        }
        root::solve(lo, hi, tolerance,
            controls.maximum_iterations, evaluate, cancelled).map_err(map_root)?
    };
    poll(cancelled)?;
    let state = model.evaluate(&updated_point(context, solved.temperature_k, solved.sample.pressure_pa)?)?;
    let change = change_at(state.state().temperature)?;
    let residual = change - requested;
    if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j_kg
        || (requested != 0.0 && (change == 0.0 || change.is_sign_positive() != requested.is_sign_positive()))
    {
        return Err(GasEnergyError::NotConverged { iterations: solved.iterations, residual_j_kg: residual });
    }
    let relative = match constraint {
        Constraint::Density(density) => state.state().density / density - 1.0,
        Constraint::Pressure => state.state().pressure / pressure.value() - 1.0,
    };
    if !relative.is_finite() || relative.abs() > 64.0 * f64::EPSILON {
        return Err(GasEnergyError::Invalid("caloric reconstruction did not retain its mechanical constraint"));
    }
    poll(cancelled)?;
    Ok(SolvedIncrement { state, change, residual, iterations: solved.iterations })
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod prepared_tests;
