//! Material-bound piston inertia from the existing rigid-body compiler.
//!
//! The body owns regional geometry, density queries and spatial inertia. The
//! piston owns the translational constraint and thermal/momentum time step.
//! Binding them does not create a second mass calculator or rigid-body solver.

/// Shared-state multiport flow with the original source-resolved wall inertias.
pub mod graph;

use super::{ClosedGasParcel, GasPiston, PistonError, PistonInput, PistonSolve, PistonStep};
use fs_mbd_material::ResolvedBody;
use super::thermal::{ConductanceSample, HeatExchangeSolve, ThermalPistonError};
use super::exchange::thermal::NozzleThermalPistonCandidate;
use super::exchange::network::NozzleNetworkCandidate;
use super::thermal::boundary::network::{NetworkGasDrive, NetworkGasError, NetworkGasSolve};
use super::thermal::boundary::network::material::GasConductionGeometry;
use crate::circuit_enthalpy::CircuitEnthalpyError;
use crate::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use crate::material_circuit::MaterialCircuit;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_thermal::enthalpy::DepositionSimulation;
use super::thermal::double_acting::{DoubleHeatExchangeSolve, DoubleThermalPistonStep};
use super::thermal::layers::{MaterialConductionLayers, LayerConductionError};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{GasNozzle, GasNozzleStepSolve};
use super::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve, DoublePistonStep,
};

/// A rigid material-bound moving wall and its coupled source-resolved gas.
///
/// All regions of `body` translate together at the piston velocity. Rotations
/// are constrained; centre-of-mass offsets and rotational inertia therefore do
/// not replace its total translational mass with an ad hoc effective mass.
/// Rigidity, geometry coverage and guide reactions remain model assumptions.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialGasPiston {
    body: ResolvedBody,
    piston: GasPiston,
}
impl MaterialGasPiston {
    /// Bind a source-resolved rigid specimen to the existing coupled dynamics.
    ///
    /// `body` is produced by `fs_mbd_material::resolve_body`, retaining density
    /// claims, region geometry and all original material state conditions. The
    /// total mass comes ONLY from its canonical spatial inertia. There is no
    /// competing numeric mass parameter or material-name selection here.
    ///
    /// `area_m2` is the pressure-exposed projected area of the translating wall,
    /// not its volume or the gas volume. Geometry must supply that independent
    /// boundary quantity. The gas and solid may have different temperatures.
    /// This binds a new physical model; it does not replace a live body's mass,
    /// conserve momentum across a material swap or infer thermal expansion.
    pub fn new(
        gas: ClosedGasParcel, body: ResolvedBody, area_m2: f64, velocity_m_s: f64,
    ) -> Result<Self, PistonError> {
        let piston = GasPiston::new(gas, body.inertia().mass(), area_m2, velocity_m_s)?;
        Ok(Self { body, piston })
    }

    /// Attach an ideal stationary support while retaining the source-resolved
    /// rigid specimen, its physical mass, exposed area and original gas owner.
    /// A moving wall refuses; no impulse or kinetic energy is silently removed.
    pub fn with_fixed_support(mut self) -> Result<Self, PistonError> {
        self.piston = self.piston.with_fixed_support()?;
        Ok(self)
    }

    /// Original source-bound specimen, including every density usage receipt.
    /// Immutable access prevents changing mass behind the dynamics owner's back.
    #[must_use]
    pub const fn body(&self) -> &ResolvedBody { &self.body }
    /// Existing accepted dynamics and gas state; no duplicate coordinates.
    #[must_use]
    pub const fn piston(&self) -> &GasPiston { &self.piston }
    /// Accepted common thermal/flow/acoustic gas state and source model.
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel { self.piston.gas() }

    /// Advance the original coupled gas and momentum equations unchanged.
    /// Failed input, source support, numerical balance or cancellation leaves
    /// body, gas, velocity and volume intact. The solid specimen is rigid and
    /// its source state is frozen: gas heating is not silently solid heating.
    pub fn step(
        &mut self, input: PistonInput, controls: PistonSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep, PistonError> {
        self.piston.step(input, controls, cancelled)
    }

    /// Preserve the original source-resolved body while a declared reaction
    /// changes gas composition, temperature and the actual pressure/work load.
    /// The immutable candidate contains no separately supplied piston mass.
    pub fn reactive_candidate(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        kinetics: super::reactive::GasReactionSolve, controls: PistonSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::ReactivePistonCandidate<Self>, super::reactive::GasReactionError> {
        let candidate = self.piston.reactive_candidate(reaction, input, kinetics, controls, &mut cancelled)?
            .map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(super::reactive::GasReactionError::Cancelled); }
        Ok(candidate)
    }

    /// Resolve chemical rates at the same endpoint as gas pressure and wall
    /// work, retaining the original material-derived inertia and source records.
    pub fn coupled_reactive_candidate(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        chemistry: super::reactive::GasReactionSolve,
        controls: super::reactive::coupled::CoupledReactivePistonSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::coupled::CoupledReactivePistonCandidate<Self>,
        super::reactive::coupled::CoupledReactivePistonError>
    {
        let candidate = self.piston.coupled_reactive_candidate(reaction, input, chemistry, controls,
            &mut cancelled)?.map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(super::reactive::coupled::CoupledReactivePistonError::Cancelled); }
        Ok(candidate)
    }

    /// Evolve material-wall heat, chemistry and pressure motion at one endpoint
    /// while retaining this translating specimen's original sourced inertia.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_thermal_candidate<S, E>(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        thermal_body: UniformEnthalpyStepInput<'_>, chemistry: super::reactive::GasReactionSolve,
        controls: super::reactive::coupled::CoupledReactivePistonSolve,
        exchange: HeatExchangeSolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::thermal::coupled::CoupledReactiveThermalCandidate<S, Self>,
        super::reactive::thermal::coupled::CoupledReactiveThermalError<E>>
    {
        let candidate = self.piston.coupled_reactive_thermal_candidate(reaction, input, thermal_body,
            chemistry, controls, exchange, resolve, &mut cancelled)?
            .map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(super::reactive::thermal::coupled::CoupledReactiveThermalError::Cancelled); }
        Ok(candidate)
    }

    /// Same-card endpoint conductivity/enthalpy and the original gas transport
    /// drive heat and reaction feedback without separately supplied piston mass.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_thermal_layers_candidate(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        thermal_body: UniformEnthalpyStepInput<'_>, chemistry: super::reactive::GasReactionSolve,
        controls: super::reactive::coupled::CoupledReactivePistonSolve,
        exchange: HeatExchangeSolve, layers: MaterialConductionLayers<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::thermal::coupled::CoupledReactiveThermalCandidate<ResolvedMaterialStatePoint, Self>,
        super::reactive::thermal::coupled::CoupledReactiveThermalError<LayerConductionError>>
    {
        let candidate = self.piston.coupled_reactive_thermal_layers_candidate(reaction, input,
            thermal_body, chemistry, controls, exchange, layers, &mut cancelled)?
            .map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(super::reactive::thermal::coupled::CoupledReactiveThermalError::Cancelled); }
        Ok(candidate)
    }

    /// Preserve source-resolved rigid inertia while reaction products exchange
    /// heat with a finite material body and drive the original pressure solve.
    /// The thermal store may be separate from this translating rigid specimen.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_candidate<S, E>(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        thermal_body: UniformEnthalpyStepInput<'_>, kinetics: super::reactive::GasReactionSolve,
        controls: PistonSolve, exchange: HeatExchangeSolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::thermal::ReactiveThermalPistonCandidate<S, Self>,
        super::reactive::thermal::ReactiveThermalPistonError<E>>
    {
        let candidate = self.piston.reactive_thermal_candidate(reaction, input, thermal_body,
            kinetics, controls, exchange, resolve, &mut cancelled)?
            .map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(super::reactive::thermal::ReactiveThermalPistonError::Cancelled); }
        Ok(candidate)
    }

    /// Source-bound wall conductivity and enthalpy coupled to chemical gas
    /// evolution. No second gas transport snapshot or piston mass is supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_layers_candidate(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        thermal_body: UniformEnthalpyStepInput<'_>, kinetics: super::reactive::GasReactionSolve,
        controls: PistonSolve, exchange: HeatExchangeSolve,
        layers: MaterialConductionLayers<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::thermal::ReactiveThermalPistonCandidate<ResolvedMaterialStatePoint, Self>,
        super::reactive::thermal::ReactiveThermalPistonError<LayerConductionError>>
    {
        self.reactive_thermal_candidate(reaction, input, thermal_body, kinetics, controls, exchange,
            |gas, thermal| layers.resolve(gas, thermal), cancelled)
    }

    /// Refine chemistry, material heat and pressure work together while keeping
    /// the original source-resolved rigid inertia and joint publication boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_subcycled_candidate<S, E>(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        thermal_body: UniformEnthalpyStepInput<'_>, kinetics: super::reactive::GasReactionSolve,
        controls: PistonSolve, exchange: HeatExchangeSolve,
        policy: super::reactive::thermal::subcycling::ReactiveThermalSubcycling,
        resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::thermal::subcycling::ReactiveThermalSubcyclingCandidate<S, Self>,
        super::reactive::thermal::ReactiveThermalPistonError<E>>
    {
        let candidate = self.piston.reactive_thermal_subcycled_candidate(reaction, input, thermal_body,
            kinetics, controls, exchange, policy, resolve, &mut cancelled)?
            .map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(super::reactive::thermal::ReactiveThermalPistonError::Cancelled); }
        Ok(candidate)
    }

    /// Source-bound wall heat feeds the next chemical substep within the same
    /// macro interval; the returned body heat and port account cover all steps.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_layers_subcycled_candidate(
        &self, reaction: &super::reactive::GasReaction, input: PistonInput,
        thermal_body: UniformEnthalpyStepInput<'_>, kinetics: super::reactive::GasReactionSolve,
        controls: PistonSolve, exchange: HeatExchangeSolve,
        policy: super::reactive::thermal::subcycling::ReactiveThermalSubcycling,
        layers: MaterialConductionLayers<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::thermal::subcycling::ReactiveThermalSubcyclingCandidate<ResolvedMaterialStatePoint, Self>,
        super::reactive::thermal::ReactiveThermalPistonError<LayerConductionError>>
    {
        self.reactive_thermal_subcycled_candidate(reaction, input, thermal_body, kinetics, controls,
            exchange, policy, |gas, thermal| layers.resolve(gas, thermal), cancelled)
    }

    /// Keep the original material-derived inertia while reaction products drive
    /// the finite circuit/thermal network, reciprocal wall heat and motion.
    /// The network's thermal material state is explicit; it does not silently
    /// replace the rigid translating body's frozen density or geometry.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn reactive_material_network_candidate<'m>(
        &self, reaction: &super::reactive::GasReaction, binding: &MaterialNetwork<'_>,
        thermal: &DepositionSimulation<'m>, target_s: f64, circuit: &MaterialCircuit,
        ports: &[f64], electrical_bases: &[(usize, ElectricalPropertyBasis)],
        drive: NetworkGasDrive, kinetics: super::reactive::GasReactionSolve,
        controls: NetworkGasSolve, geometry: GasConductionGeometry, guess_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<super::reactive::network::ReactiveNetworkCandidate<'m,
        MaterialNetworkSample, ResolvedMaterialStatePoint, Self>,
        super::reactive::network::ReactiveNetworkError<MaterialNetworkError, MaterialNetworkError>>
    {
        let candidate = self.piston.reactive_material_network_candidate(reaction, binding,
            thermal, target_s, circuit, ports, electrical_bases, drive, kinetics, controls,
            geometry, guess_j, &mut cancelled)?
            .map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(super::reactive::network::ReactiveNetworkError::Cancelled); }
        Ok(candidate)
    }

    /// Use the same conservative finite-reservoir exchange with the original
    /// source-resolved rigid-body inertia and retained material regions.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_nozzle(
        &mut self, reservoir: &mut ClosedGasParcel,
        nozzle: fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle,
        input: PistonInput,
        exchange: fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzleStepSolve,
        controls: PistonSolve, total_energy_tolerance_j: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<super::exchange::PistonNozzleStep, PistonError> {
        self.piston.step_with_nozzle(reservoir, nozzle, input, exchange, controls,
            total_energy_tolerance_j, cancelled)
    }

    /// Retain the canonical solid inertia while exchanging nozzle species and
    /// enthalpy and solving wall heat, phase and gas/piston motion together.
    /// The thermal store may be separate from this rigid inertial specimen;
    /// heating it does not silently deform or re-resolve the mechanical body.
    #[allow(clippy::too_many_arguments)]
    pub fn nozzle_thermal_candidate<S, E>(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        input: PistonInput, thermal_body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve, flow: GasNozzleStepSolve, exchange: HeatExchangeSolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleThermalPistonCandidate<S, Self>, ThermalPistonError<E>> {
        let candidate = self.piston.nozzle_thermal_candidate(reservoir, nozzle, input, thermal_body,
            controls, flow, exchange, resolve, &mut cancelled)?
            .map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(ThermalPistonError::Cancelled); }
        Ok(candidate)
    }

    /// Resolve wall conductivity from the same material card as its enthalpy
    /// chart, and gas conductivity from the actual post-exchange gas midpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn nozzle_thermal_layers_candidate(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        input: PistonInput, thermal_body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve, flow: GasNozzleStepSolve, exchange: HeatExchangeSolve,
        layers: MaterialConductionLayers<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleThermalPistonCandidate<fs_material::state_point::ResolvedMaterialStatePoint, Self>,
        ThermalPistonError<LayerConductionError>>
    {
        self.nozzle_thermal_candidate(reservoir, nozzle, input, thermal_body, controls, flow, exchange,
            |gas, thermal| layers.resolve(gas, thermal), cancelled)
    }

    /// Retain regional material inertia while a finite circuit/thermal network
    /// exchanges heat with the actual incoming gas mixture and moving wall.
    /// Electrical and conductive coefficients share the original cell sources;
    /// source-bound gas caloric/transport laws retain the accepted composition.
    /// This returns the thermal, electrical, gas/reservoir and mechanical
    /// endpoints together. The original rigid body's density binding is kept;
    /// heating a network cell does not silently change the piston material state.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn nozzle_material_network_candidate<'m>(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        binding: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        drive: NetworkGasDrive, flow: GasNozzleStepSolve, controls: NetworkGasSolve,
        geometry: GasConductionGeometry, guess_j: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleNetworkCandidate<'m, MaterialNetworkSample, ResolvedMaterialStatePoint, Self>,
        CircuitEnthalpyError<NetworkGasError<MaterialNetworkError, MaterialNetworkError>>>
    {
        let candidate = self.piston.nozzle_material_network_candidate(reservoir, nozzle, binding,
            thermal, target_s, circuit, ports, electrical_bases, drive, flow, controls, geometry,
            guess_j, &mut cancelled)?.map_piston(|piston| Self { body: self.body.clone(), piston });
        if cancelled() { return Err(CircuitEnthalpyError::Cancelled); }
        Ok(candidate)
    }
}

impl MaterialGasPiston {
    /// Add an opposing finite gas volume without discarding the original body,
    /// density receipts, primary material history or accepted wall momentum.
    /// No second inertial mass or independent opposing velocity is admitted.
    /// This consumes the old single-sided model to prevent advancing it behind
    /// the new shared coupling's back; construction failures return no new model.
    pub fn with_opposing_gas(
        self, opposing: ClosedGasParcel, opposing_area_m2: f64,
    ) -> Result<MaterialDoubleActingGasPiston, PistonError> {
        let piston = DoubleActingGasPiston::new(self.piston, opposing, opposing_area_m2)?;
        Ok(MaterialDoubleActingGasPiston { body: self.body, piston })
    }
}

/// A source-resolved rigid specimen driven on both faces by finite gas materials.
///
/// Body density/geometry determine the ONE canonical inertia. Both original gas
/// models determine thermal storage, pressure, viscosity, conductivity and sound.
/// There is no competing mass, nominal atmosphere or constant-gamma parameter.
/// The specimen is still rigid and its material state is frozen; gas heating
/// does not implicitly heat, melt, deform or change the mass of the wall.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialDoubleActingGasPiston {
    body: ResolvedBody,
    piston: DoubleActingGasPiston,
}

impl MaterialDoubleActingGasPiston {
    /// Bind two gases and one existing source-resolved body. Areas are the
    /// positive pressure-exposed [primary, opposing] face areas [m2], supplied
    /// by geometry. Positive velocity expands the primary gas.
    ///
    /// Changing a density/material card and resolving the body before this
    /// constructor changes the subsequent motion AND both gas states without
    /// editing physics. Replacing the mass of a live system is not modeled.
    pub fn new(
        primary: ClosedGasParcel, opposing: ClosedGasParcel, body: ResolvedBody,
        areas_m2: [f64; 2], velocity_m_s: f64,
    ) -> Result<Self, PistonError> {
        MaterialGasPiston::new(primary, body, areas_m2[0], velocity_m_s)?
            .with_opposing_gas(opposing, areas_m2[1])
    }

    /// Original regional geometry, material conditions and density receipts.
    #[must_use]
    pub const fn body(&self) -> &ResolvedBody { &self.body }
    /// One accepted shared dynamics owner; cannot independently step either gas.
    #[must_use]
    pub const fn piston(&self) -> &DoubleActingGasPiston { &self.piston }

    // Preserve the complete accepted-candidate owner, including both port clocks
    // and cumulative histories. Reconstructing from gas/velocity snapshots would
    // reset that history and detach continuation from its coupled acoustic clock.
    pub(in crate::cavity_phs::piston) fn bind_piston_candidate(
        &self, piston: DoubleActingGasPiston,
    ) -> Result<Self, PistonError> {
        if piston.primary().mass_kg().to_bits() != self.piston.primary().mass_kg().to_bits()
            || piston.primary().area_m2().to_bits() != self.piston.primary().area_m2().to_bits()
            || piston.opposing_area_m2().to_bits() != self.piston.opposing_area_m2().to_bits() {
            return Err(PistonError::Invalid("a material piston candidate must retain its original inertia and face geometry"));
        }
        Ok(Self { body: self.body.clone(), piston })
    }

    /// Source-bearing [primary, opposing] gas states for existing thermal,
    /// transport and acoustic consumers. No duplicate coefficient cache.
    #[must_use]
    pub const fn gases(&self) -> [&ClosedGasParcel; 2] {
        [self.piston.primary().gas(), self.piston.opposing()]
    }

    /// Reuse the same joint material/momentum solve and atomic balance gates.
    /// No material receipt or rigid-body source state is changed by gas heating.
    pub fn step(
        &mut self, input: DoublePistonInput, controls: DoublePistonSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DoublePistonStep, PistonError> {
        self.piston.step(input, controls, cancelled)
    }
}

/// Material-bound shared mechanics and its jointly accepted finite thermal store.
/// The mechanical body's original density and geometry records remain attached.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialDoubleThermalPistonCandidate<S> {
    piston: MaterialDoubleActingGasPiston,
    thermal: UniformEnthalpyStep<[S; 2]>,
    step: DoubleThermalPistonStep,
}
impl<S> MaterialDoubleThermalPistonCandidate<S> {
    /// Accepted candidate with the original source-resolved inertial specimen.
    pub const fn piston(&self) -> &MaterialDoubleActingGasPiston { &self.piston }
    /// One finite thermal endpoint retaining both interface source evaluations.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<[S; 2]> { &self.thermal }
    /// Original joint material, heat, momentum and energy accounting.
    pub const fn step(&self) -> &DoubleThermalPistonStep { &self.step }
    /// Publish the material-bound mechanism and thermal state together.
    pub fn into_parts(self) -> (MaterialDoubleActingGasPiston, UniformEnthalpyStep<[S; 2]>, DoubleThermalPistonStep) {
        (self.piston, self.thermal, self.step)
    }
}

impl MaterialDoubleActingGasPiston {
    /// Exchange heat with one finite material store without shedding the body's
    /// source records or replacing its canonical inertia with a numeric mass.
    /// The thermal store may be separate from the rigid moving specimen; its
    /// phase/temperature do not silently deform or re-resolve that specimen.
    /// Both gas materials, one momentum and the store publish only as a complete
    /// candidate. Conductance callbacks must be repeatable and free of live debits.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_candidate<S, E>(
        &self, input: DoublePistonInput, thermal_body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve, exchange: DoubleHeatExchangeSolve,
        resolve: impl FnMut(&[ResolvedThermochemicalGasState; 2], EquilibriumPhaseState)
            -> Result<[ConductanceSample<S>; 2], E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialDoubleThermalPistonCandidate<S>, ThermalPistonError<E>> {
        let candidate = self.piston.thermal_candidate(input, thermal_body, controls,
            exchange, resolve, &mut cancelled)?;
        let (piston, thermal, step) = candidate.into_parts();
        let candidate = MaterialDoubleThermalPistonCandidate {
            piston: Self { body: self.body.clone(), piston }, thermal, step,
        };
        if cancelled() { return Err(ThermalPistonError::Cancelled); }
        Ok(candidate)
    }

    /// Use existing source-bound gas/body conduction layers on both faces while
    /// preserving the mechanical material binding. Swapping supported cards and
    /// geometry changes conduction, gas response and inertia without solver edits.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_layers_candidate(
        &self, input: DoublePistonInput, thermal_body: UniformEnthalpyStepInput<'_>,
        controls: DoublePistonSolve, exchange: DoubleHeatExchangeSolve,
        layers: [MaterialConductionLayers<'_>; 2], cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialDoubleThermalPistonCandidate<fs_material::state_point::ResolvedMaterialStatePoint>,
        ThermalPistonError<LayerConductionError>>
    {
        self.thermal_candidate(input, thermal_body, controls, exchange,
            |gases, thermal| Ok([
                layers[0].resolve(&gases[0], thermal)?,
                layers[1].resolve(&gases[1], thermal)?,
            ]), cancelled)
    }
}

#[cfg(test)]
pub(super) mod tests;

