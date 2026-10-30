//! Complete heat assemblies on the original reciprocal elastic trajectories.
use super::*;
use crate::circuit_enthalpy::{
    CircuitEnthalpyError,
    contacts::heat::{
        AssemblyHeatObserver, HeatAssembly, HeatAssemblyAdvanceError, observe_heat_assembly,
    },
};

/// Original source admission or the actual reciprocal solve/observer failure.
pub type ThermoelasticHeatAssemblyError<S, O> = HeatAssemblyAdvanceError<
    ThermoelasticHeatError<Infallible>,
    ThermoelasticFeedbackError<ObserverCouplingError<CircuitEnthalpyError<S>, O>>,
>;

impl<'s, 'm> ThermoelasticSolid<'s, 'm> {
    /// Resolve original elastic equilibrium and the COMPLETE heat assembly
    /// together. Additional prescribed reservoirs and finite coolant links
    /// share the same caloric timeline as source conduction, interfaces and
    /// radiation. Their retained source payload is available in the first
    /// observer; the prepared owner still supplies exact material identity.
    ///
    /// Reversible heat, material potentials, nodal work and closed-radiation
    /// accounting stay with the existing reciprocal solve. No source or thermal
    /// state is advanced twice. Retain the original second observer on
    /// continuation and refine the existing physical/thermal intervals.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_assembly_candidate<'t, 'p, 'scope, A: HeatAssembly + ?Sized, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        assembly: &'p A,
        inner: &O,
        target_s: f64,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: StepConfig,
        config: ThermoelasticSolve,
        cx: &'p Cx<'scope>,
    ) -> Result<(
        ThermoelasticCandidate<'s, 'm, 't,
            CoupledObservers<AssemblyHeatObserver<'p, 'scope, A>, O>>,
        SolidHeatAudit,
    ), ThermoelasticHeatAssemblyError<A::Error, O::Error>> {
        admit::<Infallible>(self.source, thermal, assembly.prepared_heat_network(), cx)
            .map_err(HeatAssemblyAdvanceError::Binding)?;
        let observer = observe_heat_assembly(
            assembly, cx, inner.clone(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Heat(e)))?;
        let candidate = self.advance(
            thermal, &observer, target_s, loads, thermal_config, config, cx,
        ).map_err(HeatAssemblyAdvanceError::Coupled)?;
        let report = &candidate.report;
        let audit = closed_account(
            report.external_heat_j, report.thermal_residual_j, report.combined_residual_j,
            candidate.observer.first.closed_radiation_net_j(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Invalid(e)))?;
        cx.checkpoint().map_err(|_| HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Cancelled))?;
        Ok((candidate, audit))
    }
}

impl<'s, 'm> DynamicThermoelasticSolid<'s, 'm> {
    /// Advance retained Newmark motion, reciprocal local heat and the complete
    /// original material/coolant assembly as one candidate.
    ///
    /// Displacement, velocity, acceleration, source potentials and the step
    /// counter remain with the native dynamic owner. Every trial restarts from
    /// the same accepted histories. Kinetic storage and exported Rayleigh loss
    /// remain in their original energy account; neither becomes unowned heat.
    ///
    /// Source refusal, competing conductance, cancellation or failed corrected
    /// energy gates publish no state. Fixed reference geometry/supports,
    /// small-strain sub-solidus source intervals and explicit thermal accuracy
    /// requirements are unchanged.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_assembly_candidate<'t, 'p, 'scope, A: HeatAssembly + ?Sized, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        assembly: &'p A,
        inner: &O,
        target_s: f64,
        loads: &TetLoadCase3<'_, '_>,
        thermal_config: StepConfig,
        config: ThermoelasticSolve,
        cx: &'p Cx<'scope>,
    ) -> Result<(
        DynamicThermoelasticCandidate<'s, 'm, 't,
            CoupledObservers<AssemblyHeatObserver<'p, 'scope, A>, O>>,
        SolidHeatAudit,
    ), ThermoelasticHeatAssemblyError<A::Error, O::Error>> {
        admit::<Infallible>(self.source, thermal, assembly.prepared_heat_network(), cx)
            .map_err(HeatAssemblyAdvanceError::Binding)?;
        let observer = observe_heat_assembly(
            assembly, cx, inner.clone(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Heat(e)))?;
        let candidate = self.advance(
            thermal, &observer, target_s, loads, thermal_config, config, cx,
        ).map_err(HeatAssemblyAdvanceError::Coupled)?;
        let report = &candidate.report;
        let audit = closed_account(
            report.external_heat_j, report.thermal_residual_j, report.combined_residual_j,
            candidate.observer.first.closed_radiation_net_j(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Invalid(e)))?;
        cx.checkpoint().map_err(|_| HeatAssemblyAdvanceError::Binding(ThermoelasticHeatError::Cancelled))?;
        Ok((candidate, audit))
    }
}
