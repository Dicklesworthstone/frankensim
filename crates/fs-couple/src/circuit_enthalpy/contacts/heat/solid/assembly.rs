//! Complete coolant/material assemblies on the original static and inertial J2 owners.
use super::*;
use super::super::{
    AssemblyHeatObserver, HeatAssembly, HeatAssemblyAdvanceError, observe_heat_assembly,
};
use crate::circuit_enthalpy::CircuitEnthalpyError;
use crate::solid_enthalpy::plastic::dynamics::{
    ThermoplasticDynamics, ThermoplasticDynamicsCandidate,
};
use std::convert::Infallible;

/// Original source-binding failure or the actual joint observer/solid failure.
pub type SolidHeatAssemblyError<S, O> = HeatAssemblyAdvanceError<
    SolidHeatError<Infallible>,
    ThermoplasticError<ObserverCouplingError<CircuitEnthalpyError<S>, O>>,
>;

impl<'m> ThermoplasticSolid<'m> {
    /// Advance the complete heat assembly and loaded plastic solid together.
    ///
    /// The assembly retains its prepared network as the original material and
    /// caloric identity authority. All additional source-resolved reservoirs and
    /// finite links use the SAME native thermal state and actual substeps.
    /// Physical plastic/reversible heat is added once by the existing solid
    /// iteration. Full assembly receipts remain in the returned first observer.
    ///
    /// Continue with the returned thermal/solid pair and the second observer's
    /// independent history. Fixed reference geometry, fully-solid source ranges,
    /// small strain/rotation and explicit thermal refinement remain required.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_assembly_candidate<'t, 'p, 'scope, A: HeatAssembly + ?Sized, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        assembly: &'p A,
        inner: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &'p Cx<'scope>,
    ) -> Result<(
        ThermoplasticCandidate<'m, 't,
            CoupledObservers<AssemblyHeatObserver<'p, 'scope, A>, O>>,
        SolidHeatAudit,
    ), SolidHeatAssemblyError<A::Error, O::Error>> {
        self.admit_heat_network::<Infallible>(
            thermal, assembly.prepared_heat_network(), cx,
        ).map_err(HeatAssemblyAdvanceError::Binding)?;
        let observer = observe_heat_assembly(
            assembly, cx, inner.clone(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(SolidHeatError::Heat(e)))?;
        let candidate = self.advance(
            thermal, &observer, target_s, prescribed_m, nodal_forces_n,
            settings, thermal_config, config, cx,
        ).map_err(HeatAssemblyAdvanceError::Coupled)?;
        let report = candidate.report();
        let audit = closed_account(
            report.external_heat_j, report.thermal_residual_j, report.combined_residual_j,
            candidate.observer().first.closed_radiation_net_j(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(SolidHeatError::Invalid(e)))?;
        cx.checkpoint().map_err(|_| HeatAssemblyAdvanceError::Binding(SolidHeatError::Cancelled))?;
        Ok((candidate, audit))
    }
}

impl<'m> ThermoplasticDynamics<'m> {
    /// Retain actual velocity, displacement, plastic history and the clock while
    /// the complete heat assembly drives the native inertial thermal iteration.
    ///
    /// Additional finite coolant/body storage stays on the original timeline.
    /// The native consistent-mass kinetic account and both numerical losses
    /// remain distinct from physical heating. The same source-identity and
    /// corrected closed-energy gates as the static entry point apply.
    /// Thermal substeps do not refine mechanical loading; refine physical
    /// intervals separately. Any late refusal leaves every input unchanged.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_heat_assembly_candidate<'t, 'p, 'scope, A: HeatAssembly + ?Sized, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        assembly: &'p A,
        inner: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: ThermoplasticConfig,
        cx: &'p Cx<'scope>,
    ) -> Result<(
        ThermoplasticDynamicsCandidate<'m, 't,
            CoupledObservers<AssemblyHeatObserver<'p, 'scope, A>, O>>,
        SolidHeatAudit,
    ), SolidHeatAssemblyError<A::Error, O::Error>> {
        self.solid().admit_heat_network::<Infallible>(
            thermal, assembly.prepared_heat_network(), cx,
        ).map_err(HeatAssemblyAdvanceError::Binding)?;
        let observer = observe_heat_assembly(
            assembly, cx, inner.clone(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(SolidHeatError::Heat(e)))?;
        let candidate = self.advance(
            thermal, &observer, target_s, prescribed_m, nodal_forces_n,
            settings, thermal_config, config, cx,
        ).map_err(HeatAssemblyAdvanceError::Coupled)?;
        let report = &candidate.report().material;
        let audit = closed_account(
            report.external_heat_j, report.thermal_residual_j, report.combined_residual_j,
            candidate.observer().first.closed_radiation_net_j(), config.energy_tolerance_j,
        ).map_err(|e| HeatAssemblyAdvanceError::Binding(SolidHeatError::Invalid(e)))?;
        cx.checkpoint().map_err(|_| HeatAssemblyAdvanceError::Binding(SolidHeatError::Cancelled))?;
        Ok((candidate, audit))
    }
}
