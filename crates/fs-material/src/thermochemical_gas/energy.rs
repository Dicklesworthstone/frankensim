//! Conservative gas-state reconstruction and fixed-volume heating.
//!
//! Solve the existing NASA-9 internal energy, not `delta_T = Q/(m cp)`.
//! The resulting temperature and density determine pressure through the same
//! ideal-gas species model; all transport and acoustic coefficients are then
//! re-evaluated together. Composition is frozen and the phase must remain gas.
//! [`parcel`] composes this heating path with reversible adiabatic volume work.
//! No reaction mechanism, piston dynamics or spatial thermal transport is
//! implied by these local conservative-state closures.

/// Isobaric enthalpy inversion and source-bound stationary-mixture caloric laws.
pub mod enthalpy;
/// Density/entropy closure for reversible adiabatic gas-state changes.
pub mod entropy;
/// Reference-invariant temperature-anchored enthalpy/internal-energy increments.
pub mod increment;
/// Closed gas storage with heating and reversible adiabatic boundary work.
pub mod parcel;
mod root;

use super::{ResolvedThermochemicalGasState, ThermochemicalGasModel, coordinates, positive, typed};
use crate::MaterialError;
use fs_matdb::QueryPoint;
use fs_qty::semantic::QuantityKind;
use fs_qty::Temperature;

/// Caller-owned temperature bracket, absolute energy tolerance and bounded work.
/// Tolerance is in J/kg, NOT relative to a potentially large formation energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasEnergySolve {
    /// Lower absolute temperature [K], inside the supported gas phase.
    pub minimum_temperature_k: f64,
    /// Upper absolute temperature [K], inside the supported gas phase.
    pub maximum_temperature_k: f64,
    /// Maximum absolute residual of the selected caloric quantity [J/kg].
    pub energy_tolerance_j_kg: f64,
    /// At most 128 safeguarded Newton/bisection iterations, plus two endpoints.
    pub maximum_iterations: usize,
}
impl GasEnergySolve {
    /// Check positive absolute temperatures and bounded iteration controls.
    /// NASA-9 and transport source coverage is checked by the original solve;
    /// the legacy air-fit envelope does not bound a source-backed gas model.
    pub fn validate(self) -> Result<(), GasEnergyError> {
        if !positive(self.minimum_temperature_k) || !positive(self.maximum_temperature_k)
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || !positive(self.energy_tolerance_j_kg) || !(1..=128).contains(&self.maximum_iterations) {
            return Err(GasEnergyError::Invalid("finite gas temperature bracket, positive absolute energy tolerance and 1..=128 iterations required"));
        }
        Ok(())
    }
}

/// Explicit failure; no partial medium or parcel update escapes.
#[derive(Clone, Debug, PartialEq)]
pub enum GasEnergyError {
    /// Invalid inputs or unrepresentable derived arithmetic.
    Invalid(&'static str),
    /// The existing thermodynamic or material owner refused the state.
    Material(MaterialError),
    /// Caller cancelled before committing any state change.
    Cancelled,
    /// Endpoint residuals do not bracket the requested energy on an increasing branch.
    OutsideBracket { lower_residual_j_kg: f64, upper_residual_j_kg: f64 },
    /// Work exhausted or floating-point stagnation without satisfying energy balance.
    NotConverged { iterations: usize, residual_j_kg: f64 },
}
impl core::fmt::Display for GasEnergyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermochemical gas energy closure: {self:?}")
    }
}
impl std::error::Error for GasEnergyError {}
impl From<MaterialError> for GasEnergyError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
fn map_root(error: root::Error<GasEnergyError>) -> GasEnergyError {
    match error {
        root::Error::Evaluation(e) => e,
        root::Error::Invalid(s) => GasEnergyError::Invalid(s),
        root::Error::Cancelled => GasEnergyError::Cancelled,
        root::Error::OutsideBracket { lower, upper } => GasEnergyError::OutsideBracket {
            lower_residual_j_kg: lower, upper_residual_j_kg: upper,
        },
        root::Error::NotConverged { iterations, residual } => GasEnergyError::NotConverged {
            iterations, residual_j_kg: residual,
        },
    }
}

/// Accepted complete material state and its actual nonlinear energy residual.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEnergySolution {
    state: ResolvedThermochemicalGasState,
    residual_j_kg: f64,
    iterations: usize,
}
impl GasEnergySolution {
    /// The same state for heat, fluid momentum and acoustic consumers.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalGasState { &self.state }
    /// `u(T) - reference_u - prescribed_change`, in J/kg.
    #[must_use]
    pub const fn residual_j_kg(&self) -> f64 { self.residual_j_kg }
    /// Interior iterations; endpoint evaluation is separately bounded by two.
    #[must_use]
    pub const fn iterations(&self) -> usize { self.iterations }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasEnergyError> {
    if cancelled() { Err(GasEnergyError::Cancelled) } else { Ok(()) }
}
fn updated_point(context: &QueryPoint, t: f64, p: f64) -> Result<QueryPoint, GasEnergyError> {
    context.clone().with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t)
        .and_then(|point| point.with_quantity("pressure", typed(QuantityKind::Pressure), p))
        .map_err(|error| super::failure("energy state coordinates", error).into())
}

impl ThermochemicalGasModel {
    /// Conservative upper bound on `dT/du` [kg K/J] for a frozen ideal gas.
    ///
    /// The original NASA-9 coefficients bound each active species' cp over the
    /// supplied source-supported bracket; crossing regions requires an explicit
    /// continuous caloric reference. Subtracting each species' own R/M
    /// bounds cv; the smallest active species bound also bounds every positive
    /// mass-weighted mixture. This is a stability bound, not a replacement
    /// caloric law: actual storage and inversion still use NASA-9 u(T).
    pub fn maximum_isochoric_temperature_slope(
        &self, controls: GasEnergySolve,
    ) -> Result<f64, GasEnergyError> {
        use fs_thermochem::{ReferenceEquationOfStateV1, StandardStatePhaseV1,
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};
        controls.validate()?;
        let convention = self.thermodynamics().components()[0].convention();
        if convention.phase() != StandardStatePhaseV1::Gas
            || convention.eos() != ReferenceEquationOfStateV1::IdealGas
            || self.thermodynamics().components().iter()
                .any(|species| species.convention() != convention)
        {
            return Err(GasEnergyError::Invalid("isochoric capacity requires one common ideal-gas convention"));
        }
        let mut lower: f64 = f64::INFINITY;
        for species in self.thermodynamics().components() {
            let gas_constant = (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                / species.molar_mass().value()).next_up();
            let cv = (enthalpy::minimum_cp(species, controls)? - gas_constant).next_down();
            if !positive(cv) {
                return Err(GasEnergyError::Invalid("source coefficient bound cannot establish positive cv; narrow the gas bracket"));
            }
            lower = lower.min(cv);
        }
        let slope = (1.0 / lower).next_up();
        if !positive(slope) {
            return Err(GasEnergyError::Invalid("isochoric temperature slope is unrepresentable"));
        }
        Ok(slope)
    }

    /// Reconstruct a common physical state from conserved density and specific
    /// internal energy. `context` supplies explicitly typed initial T/p and all
    /// other material conditions. T/p are replaced; other axes are preserved.
    /// Kinetic energy must already have been removed by the flow solver.
    ///
    /// The caller's bracket must describe one stable, continuous gas branch.
    /// Positive cv is checked at every iterate; source-domain failures propagate.
    /// A bracket width alone never counts as convergence, so a discontinuous
    /// source energy gap cannot be laundered into an accepted approximate root.
    /// Transport properties are resolved once, at the accepted state.
    pub fn reconstruct_density_energy(
        &self, context: &QueryPoint, density_kg_m3: f64, internal_energy_j_kg: f64,
        controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnergySolution, GasEnergyError> {
        self.solve_energy_change(context, density_kg_m3, 0.0, internal_energy_j_kg, controls, &mut cancelled)
    }

    /// Reconstruct from a fixed internal-energy reference and a signed change.
    ///
    /// Flow solvers may conserve sensible energy relative to an explicit
    /// reference for a frozen composition. This evaluates
    /// `(u(T) - reference_internal_energy_j_kg) - internal_energy_change_j_kg`
    /// directly instead of first adding a potentially large formation-energy
    /// offset to the transported energy. The reference must use this model's
    /// thermochemical convention and remain fixed for that conserved variable.
    /// Composition changes and reactions require their full species energy
    /// balance; an arbitrary reference does not remove chemical energy.
    ///
    /// All density, source-domain, convergence and cancellation checks are the
    /// same as `reconstruct_density_energy`. The accepted source state retains
    /// its original absolute internal energy and complete transport properties.
    pub fn reconstruct_density_energy_change(
        &self,
        context: &QueryPoint,
        density_kg_m3: f64,
        reference_internal_energy_j_kg: f64,
        internal_energy_change_j_kg: f64,
        controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasEnergySolution, GasEnergyError> {
        self.solve_energy_change(
            context,
            density_kg_m3,
            reference_internal_energy_j_kg,
            internal_energy_change_j_kg,
            controls,
            &mut cancelled,
        )
    }

    fn solve_energy_change(
        &self, context: &QueryPoint, density: f64, reference_u: f64, change: f64,
        controls: GasEnergySolve, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasEnergySolution, GasEnergyError> {
        controls.validate()?;
        coordinates(context)?;
        if !positive(density) || !reference_u.is_finite() || !change.is_finite() {
            return Err(GasEnergyError::Invalid("finite positive density and finite specific energies required"));
        }
        let thermal = |t: f64| -> Result<root::Sample, GasEnergyError> {
            // Ideal-gas u and cv do not depend on pressure. Use the retained
            // standard pressure only during inversion, then the actual EOS p.
            let reference = self.thermodynamics().components()[0].convention().reference_pressure();
            let properties = self.thermodynamics().evaluate_properties(Temperature::new(t), reference)
                .map_err(|error| super::failure("energy thermodynamics", error))?;
            let u = properties.u().value() / properties.molar_mass().value();
            Ok(root::Sample {
                residual: (u - reference_u) - change,
                tangent: properties.mass_specific_cv().value(),
                pressure_pa: density * properties.mass_specific_gas_constant().value() * t,
            })
        };
        let solved = root::solve(controls.minimum_temperature_k, controls.maximum_temperature_k,
            controls.energy_tolerance_j_kg, controls.maximum_iterations, thermal, cancelled)
            .map_err(map_root)?;
        poll(cancelled)?;
        let state = self.evaluate(&updated_point(context, solved.temperature_k, solved.sample.pressure_pa)?)?;
        let residual = (state.specific_internal_energy_j_kg() - reference_u) - change;
        if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j_kg {
            return Err(GasEnergyError::NotConverged { iterations: solved.iterations, residual_j_kg: residual });
        }
        if (state.state().density / density - 1.0).abs() > 64.0 * f64::EPSILON {
            return Err(GasEnergyError::Invalid("reconstructed EOS did not retain conserved density"));
        }
        poll(cancelled)?;
        Ok(GasEnergySolution { state, residual_j_kg: residual, iterations: solved.iterations })
    }
}

/// A closed, fixed-volume gas parcel using the same reusable energy closure.
/// Heat is prescribed energy, not a temperature animation or an implicit reaction.
#[derive(Clone, Debug, PartialEq)]
pub struct IsochoricGasParcel {
    model: ThermochemicalGasModel,
    point: QueryPoint,
    state: ResolvedThermochemicalGasState,
    mass_kg: f64,
    volume_m3: f64,
}
/// Energy balance of an accepted parcel update, in coherent SI joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasHeatStep {
    /// Signed prescribed heat: positive enters the parcel.
    pub supplied_heat_j: f64,
    /// Original-law `m Delta u` between accepted temperatures, with common
    /// formation references canceled before floating-point evaluation.
    pub stored_energy_change_j: f64,
    /// Stored energy change minus supplied heat. No formation-energy normalization.
    pub residual_j: f64,
    /// Nonlinear interior iterations used.
    pub iterations: usize,
}
impl IsochoricGasParcel {
    /// The same frozen-composition source model that owns this parcel's energy,
    /// pressure and transport. Shared immutably by coupled film evaluations.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalGasModel { &self.model }

    /// Initialize from a source-resolved temperature/pressure and explicit
    /// volume. Mass is rho*V from this SAME state, not an independent density.
    pub fn new(model: ThermochemicalGasModel, point: QueryPoint, volume_m3: f64) -> Result<Self, GasEnergyError> {
        if !positive(volume_m3) { return Err(GasEnergyError::Invalid("positive finite parcel volume required")); }
        let state = model.evaluate(&point)?;
        let mass_kg = state.state().density * volume_m3;
        if !positive(mass_kg) { return Err(GasEnergyError::Invalid("parcel mass is unrepresentable")); }
        Ok(Self { model, point, state, mass_kg, volume_m3 })
    }
    /// Current common medium for fluid, thermal and acoustic observations.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalGasState { &self.state }
    /// Conserved mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Fixed volume [m3].
    #[must_use]
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    /// Complete updated material conditions, including untouched extra axes.
    #[must_use]
    pub const fn point(&self) -> &QueryPoint { &self.point }

    /// Add signed heat atomically. Outside-domain states, insufficient work,
    /// unrepresentable balances and cancellation leave the parcel unchanged.
    /// The caller can bind cancellation to its existing runtime context.
    /// Nonzero updates require a stable continuous caloric bracket containing
    /// the initial temperature. Heat and its residual use direct original-law
    /// increments; full source references remain in the published gas state.
    pub fn add_heat(
        &mut self, heat_j: f64, controls: GasEnergySolve, mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasHeatStep, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() { return Err(GasEnergyError::Invalid("finite signed heat required")); }
        if heat_j == 0.0 {
            return Ok(GasHeatStep { supplied_heat_j: heat_j, stored_energy_change_j: 0.0, residual_j: 0.0, iterations: 0 });
        }
        let change = heat_j / self.mass_kg;
        let tolerance_j = self.mass_kg * controls.energy_tolerance_j_kg;
        if !change.is_finite() || change == 0.0 || !positive(tolerance_j) {
            return Err(GasEnergyError::Invalid("specific heat input or total energy tolerance is unrepresentable"));
        }
        let solution = self.model.reconstruct_density_internal_energy_increment(
            &self.point, self.mass_kg / self.volume_m3, change, controls, &mut cancelled,
        )?;
        let stored = self.mass_kg * self.model.thermodynamics().specific_internal_energy_change(
            Temperature::new(self.state.state().temperature),
            Temperature::new(solution.state.state().temperature),
        ).map_err(|error| super::failure("parcel heat increment", error))?.value();
        let residual = stored - heat_j;
        if !stored.is_finite() || !residual.is_finite() || residual.abs() > tolerance_j {
            return Err(GasEnergyError::Invalid("accepted parcel did not close the requested joule balance"));
        }
        let point = updated_point(&self.point, solution.state.state().temperature, solution.state.state().pressure)?;
        poll(&mut cancelled)?;
        let result = GasHeatStep { supplied_heat_j: heat_j, stored_energy_change_j: stored,
            residual_j: residual, iterations: solution.iterations };
        self.point = point;
        self.state = solution.state;
        Ok(result)
    }
}

#[cfg(test)]
#[path = "energy_tests.rs"]
mod tests;
