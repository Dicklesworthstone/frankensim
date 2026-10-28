//! Pressure/entropy reconstruction on a prepared frozen ideal-gas branch.
//! At fixed outlet pressure ds/dT=Cp(T)/T, not Cv(T)/T. Source references
//! cancel in the original relative entropy owner; there is no constant-gamma fit.

use super::{GasCaloricBranch, ResolvedThermochemicalGasState, coordinates, failure, root, updated_point};
use super::super::entropy::{GasEntropyError, map_context};
use fs_matdb::QueryPoint;
use fs_qty::{Pressure, Temperature};

/// Complete source-resolved endpoint and the physical entropy residual.
#[derive(Clone, Debug, PartialEq)]
pub struct GasPressureEntropySolution {
    state: ResolvedThermochemicalGasState,
    entropy_change_j_kg_k: f64,
    residual_j_kg_k: f64,
    iterations: usize,
}
impl GasPressureEntropySolution {
    pub fn state(&self) -> &ResolvedThermochemicalGasState { &self.state }
    pub fn entropy_change_j_kg_k(&self) -> f64 { self.entropy_change_j_kg_k }
    pub fn residual_j_kg_k(&self) -> f64 { self.residual_j_kg_k }
    pub fn iterations(&self) -> usize { self.iterations }
    pub fn into_state(self) -> ResolvedThermochemicalGasState { self.state }
}

fn map_root(error: root::Error<GasEntropyError>) -> GasEntropyError {
    match error {
        root::Error::Evaluation(error) => error,
        root::Error::Invalid(message) => GasEntropyError::Invalid(message),
        root::Error::Cancelled => GasEntropyError::Cancelled,
        root::Error::OutsideBracket { lower, upper } => GasEntropyError::OutsideBracket {
            lower_residual_j_kg_k: lower, upper_residual_j_kg_k: upper,
        },
        root::Error::NotConverged { iterations, residual } => GasEntropyError::NotConverged {
            iterations, residual_j_kg_k: residual,
        },
    }
}

impl GasCaloricBranch<'_> {
    /// Reconstruct at a typed outlet pressure and prescribed entropy CHANGE
    /// from context [J/(kg K)]. Zero change defines the isentropic comparison.
    /// Context T/p and all other query axes remain explicit; the original model
    /// and composition are frozen. The prepared branch owns continuous caloric
    /// support, and the original material owner admits final transport.
    /// Entropy tolerance is absolute, positive and independent of energy zeros.
    pub fn reconstruct_pressure_entropy_change(
        &self, context: &QueryPoint, outlet_pressure: Pressure,
        entropy_change_j_kg_k: f64, entropy_tolerance_j_kg_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasPressureEntropySolution, GasEntropyError> {
        let solved = self.pressure_entropy_temperature(context, outlet_pressure,
            entropy_change_j_kg_k, entropy_tolerance_j_kg_k, &mut cancelled)?;
        super::poll(&mut cancelled).map_err(map_context)?;
        let point = updated_point(context, solved.temperature_k, outlet_pressure.value()).map_err(map_context)?;
        let state = self.model.evaluate(&point)?;
        let change = self.model.specific_entropy_change_j_kg_k(context, &point)?;
        let residual = change - entropy_change_j_kg_k;
        if !residual.is_finite() || residual.abs() > entropy_tolerance_j_kg_k {
            return Err(GasEntropyError::NotConverged { iterations: solved.iterations, residual_j_kg_k: residual });
        }
        let pressure_error = state.state().pressure / outlet_pressure.value() - 1.0;
        if !pressure_error.is_finite() || pressure_error.abs() > 64.0 * f64::EPSILON {
            return Err(GasEntropyError::Invalid("entropy reconstruction did not retain outlet pressure"));
        }
        super::poll(&mut cancelled).map_err(map_context)?;
        Ok(GasPressureEntropySolution { state, entropy_change_j_kg_k: change,
            residual_j_kg_k: residual, iterations: solved.iterations })
    }

    // Shared caloric-only comparison for flow devices. A hypothetical isentropic
    // reference must not require transport at a state the actual device never
    // visits. Publicly accepted physical states still resolve all source cards.
    pub(super) fn pressure_entropy_temperature(
        &self, context: &QueryPoint, outlet_pressure: Pressure,
        change: f64, tolerance: f64, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<root::Solution, GasEntropyError> {
        super::poll(cancelled).map_err(map_context)?;
        let (t0, p0) = coordinates(context)?;
        if !super::positive(outlet_pressure.value()) || !change.is_finite()
            || !super::positive(tolerance) || t0.value() < self.controls.minimum_temperature_k
            || t0.value() > self.controls.maximum_temperature_k
        {
            return Err(GasEntropyError::Invalid("positive outlet pressure/tolerance, finite entropy change and supported reference temperature required"));
        }
        let evaluate = |t| -> Result<root::Sample, GasEntropyError> {
            let properties = self.model.thermodynamics().evaluate_properties(Temperature::new(t), outlet_pressure)
                .map_err(|error| failure("pressure entropy tangent", error))?;
            let entropy = self.model.thermodynamics().specific_entropy_change(t0, p0,
                Temperature::new(t), outlet_pressure).map_err(|error| failure("pressure entropy increment", error))?.value();
            Ok(root::Sample { residual: entropy - change,
                tangent: properties.mass_specific_cp().value() / t,
                pressure_pa: outlet_pressure.value() })
        };
        // The immutable branch supplies a global lower bound on ds/dT.
        // Use the initial entropy residual to narrow the solve without fitting
        // a constant gamma or admitting only the endpoints of the source path.
        let reference = evaluate(t0.value())?;
        if !super::positive(reference.tangent) || !reference.residual.is_finite() {
            return Err(GasEntropyError::Invalid("pressure/entropy reference is unrepresentable"));
        }
        if reference.residual == 0.0 {
            super::poll(cancelled).map_err(map_context)?;
            return Ok(root::Solution { temperature_k: t0.value(), sample: reference, iterations: 0 });
        }
        let lower_slope = (self.minimum_cp / self.controls.maximum_temperature_k).next_down();
        let solve_tolerance = tolerance.min(0.5 * reference.residual.abs());
        if !super::positive(lower_slope) || !super::positive(solve_tolerance) {
            return Err(GasEntropyError::Invalid("pressure/entropy displacement is below represented resolution"));
        }
        let span = (reference.residual.abs() / lower_slope).next_up();
        let (lo, hi) = if reference.residual < 0.0 {
            (t0.value(), self.controls.maximum_temperature_k.min((t0.value() + span).next_up()))
        } else {
            (self.controls.minimum_temperature_k.max((t0.value() - span).next_down()), t0.value())
        };
        root::solve(lo, hi, solve_tolerance, self.controls.maximum_iterations,
            evaluate, cancelled).map_err(map_root)
    }
}

#[cfg(test)]
mod tests;
