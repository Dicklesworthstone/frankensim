//! Source-backed reversible field storage on the original spatial circuit.
//!
//! The same native circuit/enthalpy trial owns conductivity, other material
//! resistances, charge/flux and C(T)/L(T). Field thermodynamics and its caloric
//! transfer are shared with the original lumped field-storage consumer.

use super::*;
use crate::circuit_enthalpy::field_storage::{
    self, FieldCircuitView, FieldThermalCandidate, FieldThermalError,
};

/// All original spatial, electrical and thermal histories with a physical-U
/// account. The nested native electrical record retains its psi convention.
pub type ConductionFieldThermalCandidate<'t, C = ResolvedConductorState> =
    FieldThermalCandidate<'t, ConductionCircuit<C>, ConductionCircuitStep>;

pub(super) fn field_thermal_error(error: CircuitEnthalpyError<FieldThermalError>) -> Error {
    match error {
        CircuitEnthalpyError::Source(error) => {
            CircuitEnthalpyError::Source(ConductionThermalError::Storage(error))
        }
        CircuitEnthalpyError::Invalid(error) => CircuitEnthalpyError::Invalid(error),
        CircuitEnthalpyError::Circuit(error) => CircuitEnthalpyError::Circuit(error),
        CircuitEnthalpyError::Thermal(error) => CircuitEnthalpyError::Thermal(error),
        CircuitEnthalpyError::Newton(error) => CircuitEnthalpyError::Newton(error),
        CircuitEnthalpyError::NotConverged(error) => CircuitEnthalpyError::NotConverged(error),
        CircuitEnthalpyError::Interrupted(error) => CircuitEnthalpyError::Interrupted(error),
        CircuitEnthalpyError::CircuitBudget => CircuitEnthalpyError::CircuitBudget,
        CircuitEnthalpyError::Cancelled => CircuitEnthalpyError::Cancelled,
    }
}

impl<C: PlanarConductorMaterial> ConductionCircuit<C> {
    pub(super) fn field_view(&self) -> FieldCircuitView<'_> {
        FieldCircuitView {
            graph: &self.graph,
            storage: &self.storage,
            charges_c: &self.charges_c,
            fluxes_wb: &self.fluxes_wb,
        }
    }
}

// Fixed-material public entry points keep their original report type. Their
// shared driver has no thermal field account to publish or silently discard.
impl<'t, C: PlanarConductorMaterial> ConductionCircuitEnthalpyCandidate<'t, C, FieldThermalSource> {
    pub(super) fn without_thermal_fields(self) -> ConductionCircuitEnthalpyCandidate<'t, C> {
        let report = self.report;
        ConductionCircuitEnthalpyCandidate {
            thermal: self.thermal,
            circuit: self.circuit,
            report: CircuitEnthalpyReport {
                electrical: report.electrical,
                thermal: report.thermal,
                source: report.source.network,
                nonlinear: report.nonlinear,
                maximum_cell_residual_j: report.maximum_cell_residual_j,
                circuit_trials: report.circuit_trials,
                deposited_joule_j: report.deposited_joule_j,
                non_electrical_heat_j: report.non_electrical_heat_j,
                total_energy_residual_j: report.total_energy_residual_j,
            },
        }
    }
}

impl ConductionThermalNetwork<'_> {
    /// Advance the spatial conductor, source R(T), C(T)/L(T) and finite heat.
    ///
    /// Every storage branch requires its original [`FieldThermalBinding`],
    /// including constant coefficients. The existing field owner admits the
    /// exact caloric card/context, whole fixed volume, affine source interval
    /// and positive combined heat capacity under the declared q/flux envelope.
    /// Every other resistor uses the same whole-cell material binding as
    /// [`Self::step_material_circuit_candidate`].
    ///
    /// Each trial starts from the accepted charge/flux and all caloric states.
    /// Source-pinned endpoint sigma, R, C and L drive the original descriptor;
    /// its midpoint voltage drives the actual spatial field. Local E.J and
    /// other resistor heat enter their original receiving cells exactly once.
    /// Signed field transfer -delta(psi_parameter)-delta(Ts) enters each field
    /// material's own caloric cell, separately from physical dissipation.
    /// The enclosing account uses U=psi+Ts and only external electrical/thermal
    /// supply; coefficient changes are not relabeled as external work.
    ///
    /// Geometry/topology and source ports remain fixed. This inherits the
    /// original first-order endpoint coupling, solid-density/affine field-law
    /// restrictions and native source-constraint refusals. It adds no spatial
    /// dielectric/magnetic field approximation, hysteresis or moving geometry.
    /// Publish the returned spatial circuit and thermal history together.
    #[allow(clippy::too_many_arguments)]
    pub fn step_field_circuit_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionCircuit<C>,
        ports: &[f64],
        graph_resistors: &[CellResistor],
        fields: &[FieldThermalBinding],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionFieldThermalCandidate<'t, C>, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        let destinations = graph_resistors
            .iter()
            .map(|r| r.destination)
            .collect::<Vec<_>>();
        let (thermal, circuit, coupled) = self
            .step_circuit_bound(
                thermal,
                target_s,
                circuit,
                ports,
                &destinations,
                graph_resistors,
                &[],
                Some(fields),
                assembly,
                field_controls,
                controls,
                initial_energy_increments_j,
                &mut cancelled,
            )?
            .into_parts();
        let supplied_j = coupled.electrical.native.supplied;
        let report =
            field_storage::close_report(coupled, supplied_j, controls.total_energy_tolerance_j)
                .map_err(field_thermal_error)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(FieldThermalCandidate::new(thermal, circuit, report))
    }
}
