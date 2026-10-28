//! Density/entropy reconstruction using the SAME species, EOS and transport
//! owners as density/energy reconstruction. At fixed density and composition,
//! ds/dT = cv(T)/T. The mixture owner supplies entropy including its pressure
//! and ideal-mixing terms; no constant-gamma or duplicate entropy law is used.
//!
//! Frozen composition, reversible adiabatic gas branch only. This does not
//! model shocks, reaction, entropy production, condensation or real-gas effects.

use super::{coordinates, poll, positive, root, updated_point, GasEnergyError};
use super::{ResolvedThermochemicalGasState, ThermochemicalGasModel};
use crate::MaterialError;
use fs_matdb::QueryPoint;
use fs_qty::{Pressure, Temperature};

/// Caller-owned temperature bracket and absolute mass-specific entropy tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasEntropySolve {
    /// Lower absolute temperature [K].
    pub minimum_temperature_k: f64,
    /// Upper absolute temperature [K].
    pub maximum_temperature_k: f64,
    /// Maximum absolute entropy residual [J/(kg K)], not relative to a reference.
    pub entropy_tolerance_j_kg_k: f64,
    /// At most 128 safeguarded Newton/bisection iterations.
    pub maximum_iterations: usize,
}
impl GasEntropySolve {
    pub(super) fn validate(self) -> Result<(), GasEntropyError> {
        // The producing source laws own the temperature domain, as in the
        // shared energy inversion; neither operation substitutes an air fit.
        if !positive(self.minimum_temperature_k) || !positive(self.maximum_temperature_k)
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || !positive(self.entropy_tolerance_j_kg_k)
            || !(1..=128).contains(&self.maximum_iterations)
        {
            return Err(GasEntropyError::Invalid("finite gas temperature bracket, positive absolute entropy tolerance and 1..=128 iterations required"));
        }
        Ok(())
    }
}

/// No failure exposes a partially updated physical medium.
#[derive(Clone, Debug, PartialEq)]
pub enum GasEntropyError {
    /// Invalid input or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original thermodynamic or material source refused a state.
    Material(MaterialError),
    /// Caller cancelled before publication.
    Cancelled,
    /// Requested entropy is not bracketed on an increasing branch.
    OutsideBracket { lower_residual_j_kg_k: f64, upper_residual_j_kg_k: f64 },
    /// Iteration budget or floating-point progress exhausted.
    NotConverged { iterations: usize, residual_j_kg_k: f64 },
}
impl core::fmt::Display for GasEntropyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermochemical gas entropy closure: {self:?}")
    }
}
impl std::error::Error for GasEntropyError {}
impl From<MaterialError> for GasEntropyError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
fn map_root(error: root::Error<GasEntropyError>) -> GasEntropyError {
    match error {
        root::Error::Evaluation(e) => e,
        root::Error::Invalid(s) => GasEntropyError::Invalid(s),
        root::Error::Cancelled => GasEntropyError::Cancelled,
        root::Error::OutsideBracket { lower, upper } => GasEntropyError::OutsideBracket {
            lower_residual_j_kg_k: lower, upper_residual_j_kg_k: upper,
        },
        root::Error::NotConverged { iterations, residual } => GasEntropyError::NotConverged {
            iterations, residual_j_kg_k: residual,
        },
    }
}
// The shared coordinate/poll helpers cannot perform an energy inversion.
pub(super) fn map_context(error: GasEnergyError) -> GasEntropyError {
    match error {
        GasEnergyError::Material(e) => GasEntropyError::Material(e),
        GasEnergyError::Cancelled => GasEntropyError::Cancelled,
        _ => GasEntropyError::Invalid("invalid entropy state coordinates"),
    }
}

/// One fully resolved state and its physical entropy residual.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEntropySolution {
    state: ResolvedThermochemicalGasState,
    residual_j_kg_k: f64,
    iterations: usize,
}
impl GasEntropySolution {
    /// Move the complete source-bearing state into a consumer without re-querying.
    #[must_use]
    pub fn into_state(self) -> ResolvedThermochemicalGasState { self.state }
    /// Same material state for momentum, heat and acoustic consumers.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalGasState { &self.state }
    /// Accepted specific entropy (or entropy change for the relative solve)
    /// minus the requested value [J/(kg K)].
    #[must_use]
    pub const fn residual_j_kg_k(&self) -> f64 { self.residual_j_kg_k }
    /// Interior iterations, excluding bracket and reference evaluations.
    #[must_use]
    pub const fn iterations(&self) -> usize { self.iterations }
}
impl ResolvedThermochemicalGasState {
    /// Mixture entropy [J/(kg K)], with the original species reference,
    /// pressure correction and ideal-mixing contribution retained.
    /// Refuses an unrepresentable conversion from molar to mass-specific entropy.
    #[must_use]
    pub fn specific_entropy_j_kg_k(&self) -> Result<f64, MaterialError> {
        let properties = self.thermodynamics().properties();
        let value = properties.s().value() / properties.molar_mass().value();
        if !value.is_finite() {
            return Err(super::super::invalid("mass-specific gas entropy is unrepresentable"));
        }
        Ok(value)
    }
}

impl ThermochemicalGasModel {
    /// Relative entropy from the SAME frozen-composition NASA-9/EOS owner.
    ///
    /// Cancels species reference entropy and ideal mixing before rounding.
    /// Both points need typed absolute T/p; this evaluates thermodynamics, not
    /// transport source support or a reacting/composition-changing path.
    pub fn specific_entropy_change_j_kg_k(
        &self, initial: &QueryPoint, final_point: &QueryPoint,
    ) -> Result<f64, MaterialError> {
        let (t0, p0) = coordinates(initial)?;
        let (t1, p1) = coordinates(final_point)?;
        self.thermodynamics().specific_entropy_change(t0, p0, t1, p1)
            .map(|value| value.value())
            .map_err(|error| super::super::failure("relative entropy thermodynamics", error))
    }

    /// Reconstruct the common gas state at conserved density and entropy.
    /// `context` must supply typed T/p and all other source conditions. Only
    /// T/p are replaced; composition and source selections remain on this model.
    /// The caller admits a continuous gas branch across the temperature bracket.
    /// Positive cv is checked at each query; endpoint/root checks do not prove
    /// continuous source support or establish that a physical path is reversible.
    ///
    /// Thermodynamics takes at most max_iterations+3 evaluations: one reference,
    /// two endpoints and the interior queries. Final material resolution adds
    /// one thermodynamic evaluation and one query per pure transport card.
    /// Cancellation is checked between bounded evaluations and before returning.
    pub fn reconstruct_density_entropy(
        &self, context: &QueryPoint, density_kg_m3: f64, entropy_j_kg_k: f64,
        controls: GasEntropySolve, cancelled: impl FnMut() -> bool,
    ) -> Result<GasEntropySolution, GasEntropyError> {
        self.solve_density_entropy(context, density_kg_m3, entropy_j_kg_k, controls, false, cancelled)
    }

    /// Reconstruct at a prescribed density and entropy CHANGE from `context`.
    ///
    /// A zero change follows the reversible adiabat without subtracting large
    /// absolute entropy references. The original T/p must also lie in the
    /// thermodynamic source domain. Frozen composition and source selections
    /// remain unchanged; controls, work bounds, cancellation and final transport
    /// admission are identical to the absolute-entropy closure.
    pub fn reconstruct_density_entropy_change(
        &self, context: &QueryPoint, density_kg_m3: f64, entropy_change_j_kg_k: f64,
        controls: GasEntropySolve, cancelled: impl FnMut() -> bool,
    ) -> Result<GasEntropySolution, GasEntropyError> {
        self.solve_density_entropy(context, density_kg_m3, entropy_change_j_kg_k, controls, true, cancelled)
    }

    fn solve_density_entropy(
        &self, context: &QueryPoint, density_kg_m3: f64, entropy_j_kg_k: f64,
        controls: GasEntropySolve, relative: bool, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEntropySolution, GasEntropyError> {
        controls.validate()?;
        let (initial_temperature, initial_pressure) = coordinates(context)?;
        if !positive(density_kg_m3) || !entropy_j_kg_k.is_finite() {
            return Err(GasEntropyError::Invalid("positive finite density and finite specific entropy required"));
        }
        poll(&mut cancelled).map_err(map_context)?;
        // R_mix is independent of T/p at fixed composition. Ask its existing
        // owner once, at the already-required lower bracket, never recompute it.
        let reference = self.thermodynamics().components()[0].convention().reference_pressure();
        let reference_state = self.thermodynamics()
            .evaluate_properties(Temperature::new(controls.minimum_temperature_k), reference)
            .map_err(|e| super::super::failure("entropy thermodynamics", e))?;
        let gas_constant = reference_state.mass_specific_gas_constant().value();
        let entropy_residual = |t: f64, p: f64, absolute_entropy: f64| -> Result<f64, GasEntropyError> {
            let entropy = if relative {
                self.thermodynamics().specific_entropy_change(
                    initial_temperature, initial_pressure, Temperature::new(t), Pressure::new(p))
                    .map_err(|error| super::super::failure("relative entropy thermodynamics", error))?
                    .value()
            } else {
                absolute_entropy
            };
            Ok(entropy - entropy_j_kg_k)
        };
        let evaluate = |t: f64| -> Result<root::Sample, GasEntropyError> {
            let pressure_pa = density_kg_m3 * gas_constant * t;
            if !positive(pressure_pa) {
                return Err(GasEntropyError::Invalid("entropy EOS pressure is unrepresentable"));
            }
            let properties = self.thermodynamics().evaluate_properties(Temperature::new(t), Pressure::new(pressure_pa))
                .map_err(|e| super::super::failure("entropy thermodynamics", e))?;
            Ok(root::Sample {
                residual: entropy_residual(t, pressure_pa,
                    properties.s().value() / properties.molar_mass().value())?,
                tangent: properties.mass_specific_cv().value() / t,
                pressure_pa,
            })
        };
        let solved = root::solve(controls.minimum_temperature_k, controls.maximum_temperature_k,
            controls.entropy_tolerance_j_kg_k, controls.maximum_iterations, evaluate, &mut cancelled)
            .map_err(map_root)?;
        poll(&mut cancelled).map_err(map_context)?;
        let point = updated_point(context, solved.temperature_k, solved.sample.pressure_pa).map_err(map_context)?;
        let state = self.evaluate(&point)?;
        let residual = entropy_residual(state.state().temperature, state.state().pressure,
            state.specific_entropy_j_kg_k()?)?;
        if !residual.is_finite() || residual.abs() > controls.entropy_tolerance_j_kg_k {
            return Err(GasEntropyError::NotConverged { iterations: solved.iterations, residual_j_kg_k: residual });
        }
        if (state.state().density / density_kg_m3 - 1.0).abs() > 64.0 * f64::EPSILON {
            return Err(GasEntropyError::Invalid("reconstructed EOS did not retain conserved density"));
        }
        poll(&mut cancelled).map_err(map_context)?;
        Ok(GasEntropySolution { state, residual_j_kg_k: residual, iterations: solved.iterations })
    }
}

#[cfg(test)]
mod tests;
