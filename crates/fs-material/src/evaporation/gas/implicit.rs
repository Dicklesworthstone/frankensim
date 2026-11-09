//! Simultaneous finite-atmosphere feedback on the original thermal network.
//!
//! One backward-Euler endpoint closes every native cell energy, each selected
//! volatile species inventory, and the original fixed-volume gas internal
//! energy. Films see the endpoint gas pressure, composition and temperature.
//! Evaporation carries surface vapor enthalpy; condensation carries upstream
//! gas enthalpy. Reciprocal sensible film heat is exchanged exactly once.
//!
//! The numerical coordinates are joules: native cell energy increments, a
//! scaled mass increment per volatile, and a scaled gas-temperature increment. These
//! scales condition the original Newton--Krylov solve; they do not replace the
//! temperature-dependent NASA-9 calorics or the original gas equation of state.
//! Each trial evaluates original caloric increments directly, with one original
//! gas-energy reconstruction only after convergence. A fresh physical recheck
//! uses that reconstructed gas before either endpoint can be published.
//!
//! These are first-order intervals with native reference geometry,
//! positive pure-liquid inventories, explicit species identities and caller-supplied
//! film coefficients. The endpoint solve supplies algebraic acceptance; the
//! adaptive window also controls local physical time error. The timeline clips
//! those same intervals around native births, contacts and boundary changes.
//! An empty interface roster advances original heat with unchanged gas; fully
//! empty geometry advances through the native event owner. None supplies
//! boiling, moving interfaces, dry-film onset or solid/liquid mixture advection.
//! The endpoint and event-free adaptive APIs leave target events pending;
//! the timeline API applies native events through the target inclusively.

use std::collections::BTreeMap;

use fs_solver::{NewtonError, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::BoundaryExchangeStepCoefficients;
use fs_thermal::enthalpy::evaporation::ThermalEvaporationError;
use fs_thermal::enthalpy::evaporation::implicit::{
    EvaporatingNetworkVaporEnthalpyAudit, EvaporationNetworkImplicitControls,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_thermochem::SpeciesId;

use super::{
    ClosedGasParcel, EvaporationLaw, GasEnergyError, GasEvaporationError, GasEvaporationSolve,
    GasLiquidBinding, GasSpeciesBoundaryCandidate, LiquidSurfaceState, ThermochemicalGasModel,
    nonnegative, positive,
};

mod endpoint;
mod solve;
pub use solve::{solve_deposition_with_gas, solve_deposition_with_gas_and_films};

/// Physical time accuracy on a shared finite-gas/thermal event-free window.
pub mod adaptive;

/// Native births and exposure changes around the same adaptive coupled solve.
pub mod timeline;

/// One native liquid cell, with explicit original sources and film data.
/// Each selected cell is a pure-liquid interface with its own volatile identity
/// and caloric offset. All other gas species act as nonabsorbing carriers at
/// that interface under the original supplied-film approximation. This does
/// not supply mixture-liquid activities or a multicomponent diffusion law.
/// Endpoint/event-free solves require active cells; the timeline also admits
/// declared future births, whose physical source is checked at activation.
#[derive(Clone, Copy)]
pub struct ImplicitGasSurface<'a> {
    pub cell_id: CellId,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub pure_liquid_enthalpy_range_j_kg: [f64; 2],
    pub mass_transfer_coefficient_m_s: f64,
    /// Reciprocal sensible film coefficient on the original exposed area.
    /// Vapor-carried sensible and latent energy are separate from this heat.
    pub heat_transfer_coefficient_w_m2_k: f64,
}

/// Numerical work, native physical budgets and explicitly admitted gas sources.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitGasControls {
    /// Applies to the whole coupled solve. Every uncached residual, predictor
    /// and both fresh final evaluations spend this one evaluation budget.
    pub thermal: EvaporationNetworkImplicitControls,
    /// Original gas bracket/reconstruction budget, gross donor limits,
    /// temperature-departure bounds and source compatibility tolerance.
    pub phase: GasEvaporationSolve,
    /// Positive numerical scale [J/kg] for the vapor inventory coordinate and
    /// its residual. This is not a material latent heat or a physical budget.
    pub mass_coordinate_scale_j_kg: f64,
}
impl ImplicitGasControls {
    pub fn new(phase: GasEvaporationSolve) -> Self {
        Self {
            thermal: EvaporationNetworkImplicitControls::default(),
            phase,
            mass_coordinate_scale_j_kg: 1e6,
        }
    }

    fn validate(self, duration_s: f64) -> Result<(), ImplicitGasError> {
        super::validate_controls(self.phase, duration_s, 0.0)?;
        self.thermal.newton.validate()?;
        if !positive(duration_s)
            || !positive(self.mass_coordinate_scale_j_kg)
            || !nonnegative(self.thermal.cell_energy_tolerance_j)
            || !nonnegative(self.thermal.total_energy_tolerance_j)
            || !nonnegative(self.thermal.mass_tolerance_kg)
            || !positive(self.thermal.difference_step_j)
            || self.thermal.maximum_iterations == 0
            || self.thermal.maximum_evaluations == 0
        {
            return Err(ImplicitGasError::Invalid(
                "invalid implicit finite-gas controls",
            ));
        }
        Ok(())
    }
}

/// No partial thermal or gas endpoint escapes any physical/numerical refusal.
#[derive(Debug)]
pub enum ImplicitGasError {
    Invalid(&'static str),
    Phase(GasEvaporationError),
    Thermal(ThermalEvaporationError),
    Numerical(NewtonError),
    NonConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    GasMassBalance { residual_kg: f64, tolerance_kg: f64 },
    GasEnergyBalance { residual_j: f64, tolerance_j: f64 },
    CoupledEnergyBalance { residual_j: f64, tolerance_j: f64 },
}
impl From<GasEvaporationError> for ImplicitGasError {
    fn from(error: GasEvaporationError) -> Self {
        Self::Phase(error)
    }
}
impl From<GasEnergyError> for ImplicitGasError {
    fn from(error: GasEnergyError) -> Self {
        Self::Phase(error.into())
    }
}
impl From<ThermalEvaporationError> for ImplicitGasError {
    fn from(error: ThermalEvaporationError) -> Self {
        Self::Thermal(error)
    }
}
impl From<DepositionError> for ImplicitGasError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error.into())
    }
}
impl From<NewtonError> for ImplicitGasError {
    fn from(error: NewtonError) -> Self {
        Self::Numerical(error)
    }
}
impl core::fmt::Display for ImplicitGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit finite-gas phase exchange: {self:?}")
    }
}
impl std::error::Error for ImplicitGasError {}

/// Fresh physical closure at the actual returned thermal and gas endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitGasReport {
    pub thermal: EvaporatingNetworkVaporEnthalpyAudit,
    /// Actual represented thermal phase transfer, grouped by original species.
    /// Positive evaporation and negative condensation retain separate species
    /// even when the total gas mass is unchanged.
    pub species_transferred_mass_kg: BTreeMap<SpeciesId, f64>,
    pub sensible_heat_into_liquid_j: BTreeMap<CellId, f64>,
    /// Maximum magnitude of volatile, whole-gas or carrier inventory defects.
    /// This independently gates the represented mass at the returned gas state.
    pub gas_mass_residual_kg: f64,
    /// Source-relative final-composition internal energy minus the final film
    /// streams and reciprocal sensible heat. Formation constants cancel first.
    pub gas_energy_residual_j: f64,
    /// Whole thermal-plus-gas inventory defect; independently bounded by the
    /// smaller native thermal and gas mass tolerance.
    pub mass_residual_kg: f64,
    /// Thermal plus gas energy defect, after cancelling their internal streams.
    pub energy_residual_j: f64,
    /// Native sources and prescribed baths, excluding reciprocal gas-film heat.
    pub external_heat_j: f64,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    /// Work in the single original gas caloric reconstruction.
    pub gas_iterations: usize,
    pub maximum_accepted_residual_norm_j: f64,
}

/// The native thermal network and finite gas are accepted/discarded together.
pub struct ImplicitGasCandidate<'m, 'catalog, T> {
    thermal: DepositionSimulation<'m>,
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    payload: T,
    report: ImplicitGasReport,
}
impl<'m, 'catalog, T> ImplicitGasCandidate<'m, 'catalog, T> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    pub const fn payload(&self) -> &T {
        &self.payload
    }
    pub const fn report(&self) -> &ImplicitGasReport {
        &self.report
    }
    /// The one original caloric reconstruction receipt. Its proposed input
    /// streams precede the fresh final film evaluation; their small consistency
    /// defect is independently bounded in `report`, not assumed to vanish.
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        T,
        ImplicitGasReport,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.thermal, gas, catalog, self.payload, self.report)
    }
}

#[cfg(test)]
mod tests;
