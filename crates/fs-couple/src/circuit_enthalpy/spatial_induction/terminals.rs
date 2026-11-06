//! Finite caloric state paired with its original physical electrode history.
//!
//! The nonlinear and hysteretic children use the existing material/enthalpy
//! residual. This shared candidate keeps field and caloric publication atomic
//! without exposing an electrode field as an isolated conductor.

use super::{MaterialNetworkSample, SpatialInductionThermalReport};
use fs_thermal::enthalpy::DepositionSimulation;

/// A jointly accepted native voltage- or current-driven field and finite heat.
pub struct SpatialTerminalThermalCandidate<'t, M, D, S = MaterialNetworkSample> {
    pub(super) thermal: DepositionSimulation<'t>,
    pub(super) field: M,
    pub(super) report: SpatialInductionThermalReport<D, S>,
}

impl<'t, M, D, S> SpatialTerminalThermalCandidate<'t, M, D, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    /// The original electrode owner, retaining its voltage/current causality.
    pub fn field(&self) -> &M {
        &self.field
    }
    pub fn report(&self) -> &SpatialInductionThermalReport<D, S> {
        &self.report
    }
    /// Native physical losses already deposited, including pinning when present.
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.coupled.deposited_joule_j
    }
    /// Publish both histories together. Physical heat has already been paid;
    /// numerical loss is separate and must not be deposited by the caller.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        M,
        SpatialInductionThermalReport<D, S>,
    ) {
        (self.thermal, self.field, self.report)
    }
}
