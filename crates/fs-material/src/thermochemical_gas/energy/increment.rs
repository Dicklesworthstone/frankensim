//! Temperature-anchored caloric inversion for one unchanged gas composition.
//!
//! Trial residuals integrate the original Cp or Cv between the explicit context
//! temperature and the trial temperature. Absolute formation references remain
//! in the published state, but cannot erase the requested sensible increment.
//! Composition changes must continue to use the full species energy balance.

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

#[derive(Clone, Copy)]
enum Constraint { Density(f64), Pressure }
struct SolvedIncrement {
    state: ResolvedThermochemicalGasState,
    change: f64,
    residual: f64,
    iterations: usize,
}

impl ThermochemicalGasModel {
    /// Reconstruct a frozen gas at the specified density from a temperature-
    /// anchored internal-energy increment [J/kg]. The context's temperature is
    /// the reference; its pressure and other axes remain explicit source inputs.
    ///
    /// Unlike an arbitrary numerical energy zero, this reference can be canceled
    /// analytically by the original NASA-9 owner. The bracket must contain the
    /// reference temperature and one stable continuous caloric branch. Raw
    /// multi-region sources require explicit caloric reconciliation first.
    /// Transport is fully resolved only at the accepted endpoint. No nonzero
    /// heat increment is accepted with an unchanged or oppositely moving state.
    /// This is not a composition-change or reaction-energy reconstruction.
    pub fn reconstruct_density_internal_energy_increment(
        &self, context: &QueryPoint, density_kg_m3: f64, change_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnergySolution, GasEnergyError> {
        let solved = solve(self, context, Constraint::Density(density_kg_m3),
            change_j_kg, controls, &mut cancelled)?;
        Ok(GasEnergySolution { state: solved.state, residual_j_kg: solved.residual,
            iterations: solved.iterations })
    }

    /// Reconstruct a frozen gas at context pressure from Delta h [J/kg], with
    /// context temperature as the explicit caloric reference. Integrates Cp(T),
    /// not a frozen Cp times Delta T or a subtraction of absolute enthalpies.
    ///
    /// Positive Cv and continuous source support are admitted over the whole
    /// bracket, including the reference temperature. The result retains the
    /// absolute formation enthalpies and the original transport/acoustic state.
    /// No chemical energy is discarded: composition and source laws stay fixed.
    pub fn reconstruct_pressure_enthalpy_increment(
        &self, context: &QueryPoint, change_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnthalpyIncrementSolution, GasEnergyError> {
        let solved = solve(self, context, Constraint::Pressure, change_j_kg,
            controls, &mut cancelled)?;
        Ok(GasEnthalpyIncrementSolution { state: solved.state,
            enthalpy_change_j_kg: solved.change, residual_j_kg: solved.residual,
            iterations: solved.iterations })
    }
}

fn solve(
    model: &ThermochemicalGasModel, context: &QueryPoint, constraint: Constraint,
    requested: f64, controls: GasEnergySolve, cancelled: &mut impl FnMut() -> bool,
) -> Result<SolvedIncrement, GasEnergyError> {
    controls.validate()?;
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
    // Initialization only: the existing bounded source analysis checks every
    // intervening region and positive Cv, not merely the two endpoints.
    let mut minimum_capacity = f64::INFINITY;
    for species in model.thermodynamics().components() {
        poll(cancelled)?;
        let cp = enthalpy::minimum_cp(species, controls)?;
        let lower = match constraint {
            Constraint::Pressure => cp,
            Constraint::Density(_) => (cp - (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                / species.molar_mass().value()).next_up()).next_down(),
        };
        minimum_capacity = minimum_capacity.min(lower);
    }
    if !positive(minimum_capacity) {
        return Err(GasEnergyError::Invalid("positive represented caloric capacity bound required"));
    }
    // The source-derived minimum capacity bounds the possible temperature
    // displacement. Outward rounding narrows tiny-heat solves without guessing
    // a constant Cp/Cv or changing the supplied physical residual tolerance.
    let span = (requested.abs() / minimum_capacity).next_up();
    let (lo, hi) = if requested > 0.0 {
        (initial.value(), controls.maximum_temperature_k.min((initial.value() + span).next_up()))
    } else {
        (controls.minimum_temperature_k.max((initial.value() - span).next_down()), initial.value())
    };
    let change_at = |t| -> Result<f64, GasEnergyError> {
        let end = Temperature::new(t);
        let value = match constraint {
            Constraint::Density(_) => model.thermodynamics()
                .specific_internal_energy_change(initial, end),
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
        root::solve(lo, hi,
            controls.energy_tolerance_j_kg, controls.maximum_iterations, evaluate, cancelled)
            .map_err(map_root)?
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
