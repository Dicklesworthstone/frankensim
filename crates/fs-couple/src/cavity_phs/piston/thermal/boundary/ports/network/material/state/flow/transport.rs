//! All declared gas transport laws feed ONE inventory before the original coupled step.
//! No accepted gas, circuit, thermal or mechanical history changes during a trial.
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_PARCELS,
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::{diffusion::network::DiffusionNetworkLink,
        network::transport::{self as shared, GasTransportError, GasTransportPaths, GasTransportSolve, GasTransportStep}},
};

impl From<GasTransportError> for MaterialFlowError {
    fn from(error: GasTransportError) -> Self {
        match error { GasTransportError::Cancelled => Self::Cancelled, other => Self::Transport(other) }
    }
}
impl<E> From<GasTransportError> for MaterialGasLawFlowError<E> {
    fn from(error: GasTransportError) -> Self {
        match error { GasTransportError::Cancelled => Self::Cancelled, other => Self::Transport(other) }
    }
}
enum TransportStageError {
    Invalid(&'static str),
    Source(GasTransportError),
}
impl From<GasTransportError> for TransportStageError {
    fn from(error: GasTransportError) -> Self { Self::Source(error) }
}
impl From<TransportStageError> for MaterialFlowError {
    fn from(error: TransportStageError) -> Self {
        match error {
            TransportStageError::Invalid(message) => Self::Invalid(message),
            TransportStageError::Source(error) => error.into(),
        }
    }
}
impl<E> From<TransportStageError> for MaterialGasLawFlowError<E> {
    fn from(error: TransportStageError) -> Self {
        match error {
            TransportStageError::Invalid(message) => Self::Invalid(message),
            TransportStageError::Source(error) => error.into(),
        }
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Preview pressure restrictions and molecular diffusion on one accepted
    /// snapshot, then react the combined inventories and run the ORIGINAL joint
    /// electrical/thermal/moving-wall solve. All stages use ONE physical interval.
    ///
    /// Lists declare distinct physical paths, with their own source-backed laws;
    /// the adapter does not add a nozzle velocity to a molar-average diffusion
    /// velocity on the same aperture. Diffusion retains its explicit admitted
    /// temperature/pressure differences and source domains. Heating or motion
    /// can invalidate that admission on the next interval; nothing is clamped.
    ///
    /// All outgoing species across BOTH mechanisms share the same original
    /// inventory allowance. Incoming material cannot replenish that allowance
    /// or traverse a second path in this explicit stage. The combined transport
    /// state is reconstructed once per volume and chemistry runs once per bound
    /// chamber. Nonlinear heat/electrical/mechanical trials never repeat either.
    ///
    /// Partial enthalpy and formation energy are already included in the gas;
    /// neither is added again as heat. Existing whole-assembly mass/energy checks
    /// enclose the original pre-transport gas, circuit, thermal and kinetic state.
    /// Any GROSS species transfer makes the closed two-port chamber audit
    /// unavailable even at zero net mass, while retaining the enclosing account.
    ///
    /// This is explicit first-order transport then reaction and joint feedback,
    /// not implicit CFD, stiff kinetics, shared crank motion or jet thrust.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_step_candidate(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], nozzles: &[SnapshotNozzleLink],
        diffusion: &[DiffusionNetworkLink<'_>], transport_controls: GasTransportSolve,
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialFlowAdvance<GasTransportStep>), MaterialFlowError> {
        self.transport_paths_step_candidate(prepared, cx, target_s, ports, drives,
            reactions, gas_patches, GasTransportPaths { nozzles, diffusion, ..Default::default() },
            transport_controls, controls, guess_j)
    }

    /// The same combined transport with a caller-selected gas heat law.
    ///
    /// Nozzle and molecular fluxes sample one accepted inventory and reconstruct
    /// each volume once. The law receives the transported/reacted gas at each
    /// actual thermal/mechanical trial, so evolving composition changes heat
    /// transfer without capturing initial density, viscosity or conductivity.
    /// Surface, reaction, withdrawal, source-domain and whole-interval balance
    /// gates are identical to the original owners. Empty connection lists admit
    /// either mechanism alone; no second transport or heat integrator is used.
    ///
    /// Distinct nozzle and diffusion paths, the molar-average diffusion frame,
    /// prescribed surface-flow geometry and first-order splitting still apply.
    /// Source payloads remain typed; all endpoint owners publish together.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn transport_step_with_gas_law_candidate<I, E>(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface], nozzles: &[SnapshotNozzleLink],
        diffusion: &[DiffusionNetworkLink<'_>], transport_controls: GasTransportSolve,
        controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, MaterialFlowAdvance<GasTransportStep, I>), MaterialGasLawFlowError<E>> {
        self.transport_paths_step_with_gas_law_candidate(prepared, cx, target_s,
            ports, drives, reactions, surfaces,
            GasTransportPaths { nozzles, diffusion, ..Default::default() },
            transport_controls, controls, guess_j, resolve)
    }

    /// Preview nozzle, molecular, porous and selective membrane transport in
    /// the original reacting, electrical, thermal and moving-chamber assembly.
    ///
    /// Every path samples the same initial inventory and uses one combined
    /// species/enthalpy update. Reaction and coupled heat/work run once after
    /// transport. Total/partial-pressure driving-force checks are reapplied to
    /// the actual final chambers, including heating, reaction and piston work.
    /// Failure leaves every accepted owner unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_paths_step_candidate(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve, controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialFlowAdvance<GasTransportStep>), MaterialFlowError> {
        let stage = self.transport_stage(cx, target_s, drives, paths, transport_controls, controls)?;
        let (next, report) = self.finish_flow_candidate(prepared, cx, target_s, ports, drives,
            reactions, gas_patches, controls, guess_j, stage)?;
        next.admit_transport_endpoint(&report.flow, cx)?;
        Ok((next, report))
    }

    /// The complete transport-path set with the original caller-selected gas
    /// heat law. Evolving mixture properties reach the same nonlinear heat/work
    /// trials; transport and chemistry are never repeated inside those trials.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn transport_paths_step_with_gas_law_candidate<I, E>(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface], paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve, controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, MaterialFlowAdvance<GasTransportStep, I>), MaterialGasLawFlowError<E>> {
        let stage = self.transport_stage(cx, target_s, drives, paths, transport_controls, controls)?;
        let (next, report) = self.finish_flow_with_gas_law_candidate(prepared, cx, target_s,
            ports, drives, reactions, surfaces, controls, guess_j, stage, resolve)?;
        next.admit_transport_endpoint(&report.flow, cx)?;
        Ok((next, report))
    }

    fn admit_transport_endpoint(&self, report: &GasTransportStep, cx: &Cx<'_>)
        -> Result<(), GasTransportError>
    {
        let mut parcels = BTreeMap::new();
        for (&id, chamber) in &self.chambers {
            cx.checkpoint().map_err(|_| GasTransportError::Cancelled)?;
            parcels.insert(id, chamber.gas().clone());
        }
        report.admit_endpoint_driving_forces(&parcels, || cx.checkpoint().is_err())
    }

    #[allow(clippy::too_many_arguments)]
    fn transport_stage(
        &self, cx: &Cx<'_>, target_s: f64, drives: &[ChamberDrive],
        paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve, controls: NetworkSolve,
    ) -> Result<FlowStage<GasTransportStep>, TransportStageError> {
        let check = || cx.checkpoint().map_err(|_| GasTransportError::Cancelled);
        check()?;
        let dt = target_s - self.time_s();
        if !dt.is_finite() || dt <= 0.0 || self.chambers.is_empty()
            || self.chambers.len() > controls.maximum_chambers
            || self.chambers.len() > MAX_SNAPSHOT_PARCELS
            || drives.len() != self.chambers.len() {
            return Err(TransportStageError::Invalid(
                "one forward clock, bounded combined paths and one drive per live chamber required"));
        }
        let mut endpoints = BTreeMap::new();
        for drive in drives {
            check()?;
            if !self.chambers.contains_key(&drive.chamber)
                || endpoints.insert(drive.chamber, drive.piston_solve.gas).is_some() {
                return Err(TransportStageError::Invalid("duplicate drive or unknown live chamber"));
            }
        }
        let mut parcels = BTreeMap::new();
        for (&id, chamber) in &self.chambers {
            check()?;
            parcels.insert(id, chamber.gas().clone());
        }
        let (parcels, flow) = shared::step_candidate_with_paths(&parcels, paths, dt,
            &endpoints, transport_controls, || cx.checkpoint().is_err())?.into_parts();
        let moved = flow.exchange.parcels.iter()
            .filter(|p| p.incoming_mass_kg > 0.0 || p.outgoing_mass_kg > 0.0)
            .map(|p| p.id).collect();
        Ok(FlowStage { parcels, moved, report: flow,
            mass_tolerance_kg: transport_controls.exchange.mass_tolerance_kg })
    }

    /// Publish all gas, charge/flux, enthalpy, wall motion and port histories
    /// together, once, after combined transport and all coupled balances pass.
    /// Failed source, budget, nonlinear or final-cancellation checks leave the
    /// original time and every physical owner intact for a retry.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_transport(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], nozzles: &[SnapshotNozzleLink],
        diffusion: &[DiffusionNetworkLink<'_>], transport_controls: GasTransportSolve,
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance<GasTransportStep>, MaterialFlowError> {
        self.advance_with_transport_paths(prepared, cx, target_s, ports, drives, reactions,
            gas_patches, GasTransportPaths { nozzles, diffusion, ..Default::default() },
            transport_controls, controls, guess_j)
    }

    /// Atomically publish combined species transport and the selected gas heat
    /// law with the original reacting, electrical, thermal and moving owners.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_transport_and_gas_law<I, E>(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface], nozzles: &[SnapshotNozzleLink],
        diffusion: &[DiffusionNetworkLink<'_>], transport_controls: GasTransportSolve,
        controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<MaterialFlowAdvance<GasTransportStep, I>, MaterialGasLawFlowError<E>> {
        self.advance_with_transport_paths_and_gas_law(prepared, cx, target_s, ports,
            drives, reactions, surfaces,
            GasTransportPaths { nozzles, diffusion, ..Default::default() },
            transport_controls, controls, guess_j, resolve)
    }

    /// Atomically publish all transport mechanisms and the original coupled
    /// owners after final driving forces, whole-interval balances and cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_transport_paths(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve, controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<MaterialFlowAdvance<GasTransportStep>, MaterialFlowError> {
        let (next, report) = self.transport_paths_step_candidate(prepared, cx,
            target_s, ports, drives, reactions, gas_patches, paths,
            transport_controls, controls, guess_j)?;
        checkpoint(cx)?;
        *self = next;
        Ok(report)
    }

    /// Publish the same complete gas paths with a source-selected gas heat law.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_transport_paths_and_gas_law<I, E>(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface], paths: GasTransportPaths<'_>,
        transport_controls: GasTransportSolve, controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<MaterialFlowAdvance<GasTransportStep, I>, MaterialGasLawFlowError<E>> {
        let (next, report) = self.transport_paths_step_with_gas_law_candidate(prepared, cx,
            target_s, ports, drives, reactions, surfaces, paths,
            transport_controls, controls, guess_j, resolve)?;
        cx.checkpoint().map_err(|_| MaterialGasLawFlowError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
