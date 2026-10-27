//! Finite gas plumbing, material heat and native elastic cavity motion.
//!
//! Every inlet/outlet samples one accepted gas snapshot. The original network
//! reconstructs all extensive inventories once; the original modal owner then
//! solves its pressure/heat/chemistry endpoint over the SAME native sample.
//! No flow law, mixing model, reaction law or wall integrator is duplicated.

use std::collections::BTreeMap;
use super::{GasValveMotion, ModalGasError, ModalGasSolve, ModalGasStep, modal_energy};
use super::reactive::thermal::{
    CoupledReactiveModalThermalSolve, CoupledReactiveModalThermalStep,
    ReactiveModalThermalError,
};
use crate::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
use crate::cavity_phs::piston::thermal::ConductanceSample;
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use fs_exec::Cx;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::{
    GasSnapshotSolve, MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::{
    self, SnapshotNozzleError, SnapshotNozzleLink, SnapshotNozzleStep,
};

/// Physical identities and controls for an elastic cavity and finite plumbing.
/// Links may join any two nodes, including two reservoirs. The native network
/// retains aggregate gross withdrawal, shared sources and pressure-order gates.
#[derive(Clone, Copy, Debug)]
pub struct ModalGasNetwork<'a> {
    pub chamber_id: u64,
    pub reservoirs: &'a BTreeMap<u64, ClosedGasParcel>,
    pub links: &'a [SnapshotNozzleLink],
    /// Exactly one source-inversion policy for every reservoir and the chamber.
    /// The chamber policy must equal its subsequent mechanical/chemical policy.
    pub endpoints: &'a BTreeMap<u64, GasEnergySolve>,
    pub exchange: GasSnapshotSolve,
    pub total_energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum ModalGasNetworkError<E> {
    Invalid(&'static str),
    Flow(SnapshotNozzleError),
    Evolution(E),
    Cancelled,
}
impl<E: core::fmt::Display> core::fmt::Display for ModalGasNetworkError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "elastic gas network: {message}"),
            Self::Flow(error) => write!(f, "elastic gas network flow: {error}"),
            Self::Evolution(error) => write!(f, "elastic gas network evolution: {error}"),
            Self::Cancelled => write!(f, "elastic gas network cancelled"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ModalGasNetworkError<E> {}

/// Local solver accounts and an independent enclosing physical balance.
#[derive(Clone, Debug)]
pub struct ModalGasNetworkStep<S> {
    pub flow: SnapshotNozzleStep,
    /// The original pressure/motion or reacting/thermal result. Its local
    /// storage and entropy window starts AFTER species transport.
    pub evolution: S,
    /// All final minus initial gas internal energy, retaining formation energy.
    pub gas_energy_change_j: f64,
    pub body_energy_change_j: f64,
    pub modal_energy_change_j: f64,
    /// Original physical modal loss leaves this enclosure, without heat credit.
    pub modal_dissipation_j: f64,
    pub exterior_pressure_work_j: f64,
    /// Gas plus finite-body external heat, excluding internal stream enthalpy.
    pub supplied_heat_j: f64,
    pub total_mass_residual_kg: f64,
    /// dU_all_gases+dH_body+dE_modes+D_modes+W_exterior-Q_external.
    pub total_energy_residual_j: f64,
}

/// One complete proposal. Publication never leaves a depleted reservoir beside
/// an old wall, old reacting mixture or old thermal body.
#[derive(Clone, Debug)]
pub struct ModalGasNetworkCandidate<S, T = ()> {
    motion: GasValveMotion,
    reservoirs: BTreeMap<u64, ClosedGasParcel>,
    thermal: T,
    step: ModalGasNetworkStep<S>,
}
impl<S, T> ModalGasNetworkCandidate<S, T> {
    pub const fn motion(&self) -> &GasValveMotion { &self.motion }
    pub fn reservoirs(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.reservoirs }
    pub const fn thermal(&self) -> &T { &self.thermal }
    pub const fn step(&self) -> &ModalGasNetworkStep<S> { &self.step }
    pub fn into_parts(self) -> (GasValveMotion, BTreeMap<u64, ClosedGasParcel>, T, ModalGasNetworkStep<S>) {
        (self.motion, self.reservoirs, self.thermal, self.step)
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), ModalGasNetworkError<E>> {
    cx.checkpoint().map_err(|_| ModalGasNetworkError::Cancelled)
}

struct Transferred {
    motion: GasValveMotion,
    reservoirs: BTreeMap<u64, ClosedGasParcel>,
    flow: SnapshotNozzleStep,
}

impl ModalGasNetwork<'_> {
    fn prepare<E>(&self, old: &GasValveMotion, gas: GasEnergySolve, cx: &Cx<'_>)
        -> Result<Transferred, ModalGasNetworkError<E>> {
        use ModalGasNetworkError as Error;
        poll(cx)?;
        if self.reservoirs.len() >= MAX_SNAPSHOT_PARCELS
            || self.links.len() > MAX_SNAPSHOT_TRANSFERS
            || self.reservoirs.contains_key(&self.chamber_id)
            || self.endpoints.get(&self.chamber_id) != Some(&gas)
            || !self.total_energy_tolerance_j.is_finite() || self.total_energy_tolerance_j <= 0.0
        {
            return Err(Error::Invalid("bounded disjoint gas owners, shared chamber inversion and positive energy tolerance required"));
        }
        let mut initial = self.reservoirs.clone();
        initial.insert(self.chamber_id, old.gas.clone());
        let candidate = snapshot::step_candidate(&initial, self.links,
            old.motion.model().sample_period_s(), self.endpoints, self.exchange,
            || cx.checkpoint().is_err()).map_err(|error| match error {
                SnapshotNozzleError::Cancelled => Error::Cancelled,
                other => Error::Flow(other),
            })?;
        let (mut reservoirs, flow) = candidate.into_parts();
        let gas = reservoirs.remove(&self.chamber_id)
            .ok_or(Error::Invalid("native flow omitted the elastic chamber"))?;
        let mut motion = old.clone();
        // Preserve the anchored volume chart, source-derived modes, equilibrium,
        // observation maps, accepted clock and every structural history.
        motion.gas = gas;
        poll(cx)?;
        Ok(Transferred { motion, reservoirs, flow })
    }

    #[allow(clippy::too_many_arguments)]
    fn finish<S, T, E>(&self, old: &GasValveMotion, next: Transferred, thermal: T,
        evolution: S, body_energy_j: f64, loss_j: f64, exterior_pressure_pa: f64,
        swept_volume_m3: f64, external_heat_j: f64, cx: &Cx<'_>)
        -> Result<ModalGasNetworkCandidate<S, T>, ModalGasNetworkError<E>> {
        use ModalGasNetworkError as Error;
        let mut mass = Sum::default();
        let mut storage = Sum::default();
        // Stable differences avoid subtracting large total formation energies.
        let mut add = |before: &ClosedGasParcel, after: &ClosedGasParcel|
            -> Result<(), ModalGasNetworkError<E>> {
            mass.add(after.mass_kg() - before.mass_kg())?;
            storage.add(after.mass_kg() * (after.state().specific_internal_energy_j_kg()
                - before.state().specific_internal_energy_j_kg())
                + (after.mass_kg() - before.mass_kg()) * before.state().specific_internal_energy_j_kg())
        };
        add(&old.gas, &next.motion.gas)?;
        for (&id, before) in self.reservoirs {
            poll(cx)?;
            add(before, &next.reservoirs[&id])?;
        }
        let modal_energy_j = modal_energy(&next.motion.motion) - modal_energy(&old.motion);
        let exterior_work_j = exterior_pressure_pa * swept_volume_m3;
        let total = storage.value + body_energy_j + modal_energy_j + loss_j
            + exterior_work_j - external_heat_j;
        if [body_energy_j, modal_energy_j, loss_j, exterior_work_j, external_heat_j, total]
            .iter().any(|x| !x.is_finite())
            || mass.value.abs() > self.exchange.mass_tolerance_kg
            || total.abs() > self.total_energy_tolerance_j
        {
            return Err(Error::Invalid("final gas inventories, thermal body and native modal motion do not close mass and energy"));
        }
        poll(cx)?;
        Ok(ModalGasNetworkCandidate { motion: next.motion, reservoirs: next.reservoirs, thermal,
            step: ModalGasNetworkStep { flow: next.flow, evolution,
                gas_energy_change_j: storage.value, body_energy_change_j: body_energy_j,
                modal_energy_change_j: modal_energy_j, modal_dissipation_j: loss_j,
                exterior_pressure_work_j: exterior_work_j, supplied_heat_j: external_heat_j,
                total_mass_residual_kg: mass.value, total_energy_residual_j: total } })
    }
}

#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add<E>(&mut self, term: f64) -> Result<(), ModalGasNetworkError<E>> {
        let adjusted = term - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if [term, adjusted, next, correction].iter().any(|x| !x.is_finite()) {
            return Err(ModalGasNetworkError::Invalid("extensive gas network sum is unrepresentable"));
        }
        self.value = next;
        self.correction = correction;
        Ok(())
    }
}

impl GasValveMotion {
    /// Exchange species/enthalpy through arbitrary finite plumbing and advance
    /// the original gas-loaded modes, at one native physical sample.
    ///
    /// All links see the same accepted pressures/compositions; incoming species
    /// cannot be sent straight out another outlet during this explicit sample.
    /// The original midpoint heat/work and held-pressure structural response
    /// then start from the jointly reconstructed chamber. The volume chart and
    /// structural history are never rebound or reset.
    ///
    /// This is a first-order transport/motion split. Flow direction may reverse
    /// after motion; the next sample observes that state. Refine the whole
    /// sample or split at a prescribed event; there is no flow event locator,
    /// jet thrust, implicit flow stability or full-network entropy claim.
    pub fn network_candidate(&self, network: ModalGasNetwork<'_>, gas_heat_j: f64,
        exterior_pressure_pa: f64, controls: ModalGasSolve, cx: &Cx<'_>)
        -> Result<ModalGasNetworkCandidate<ModalGasStep>, ModalGasNetworkError<ModalGasError>> {
        let mut next = network.prepare(self, controls.gas, cx)?;
        let step = next.motion.step(gas_heat_j, exterior_pressure_pa, controls, cx)
            .map_err(|error| match error { ModalGasError::Cancelled => ModalGasNetworkError::Cancelled,
                other => ModalGasNetworkError::Evolution(other) })?;
        let loss = step.motion.modal.viscous_dissipation_j;
        let swept = step.motion.swept_volume_change_m3;
        network.finish(self, next, (), step, 0.0, loss, exterior_pressure_pa, swept, gas_heat_j, cx)
    }

    /// Flow once, then solve the original simultaneous chemical/heat/pressure/
    /// elastic endpoint. Every nonlinear trial restarts the SAME transferred
    /// chamber; no reservoir is withdrawn again during heat or chemistry probes.
    ///
    /// Formation energy and incoming total enthalpy are already inside gas
    /// storage. Both explicit heat inputs are external only. The finite thermal
    /// body remains separate from the source-fixed elastic wall; no unmodelled
    /// thermal expansion, basis update, flame, shock or jet momentum is inferred.
    /// The local thermal entropy account covers the post-flow window only.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn reactive_thermal_network_candidate<S, E>(&self, network: ModalGasNetwork<'_>,
        reaction: &GasReaction, gas_heat_j: f64, exterior_pressure_pa: f64,
        body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
        controls: CoupledReactiveModalThermalSolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>, cx: &Cx<'_>)
        -> Result<ModalGasNetworkCandidate<CoupledReactiveModalThermalStep, UniformEnthalpyStep<S>>,
            ModalGasNetworkError<ReactiveModalThermalError<E>>> {
        let mut next = network.prepare(self, controls.coupled.mechanical.gas, cx)?;
        let candidate = next.motion.coupled_reactive_thermal_candidate(reaction, gas_heat_j,
            exterior_pressure_pa, body, chemistry, controls, resolve, cx)
            .map_err(|error| match error {
                ReactiveModalThermalError::Cancelled => ModalGasNetworkError::Cancelled,
                other => ModalGasNetworkError::Evolution(other),
            })?;
        let (motion, thermal, step) = candidate.into_parts();
        next.motion = motion;
        let body_energy = body.mass_kg * (thermal.state.specific_enthalpy_j_kg()
            - body.initial.specific_enthalpy_j_kg());
        let loss = step.coupled.motion.modal.viscous_dissipation_j;
        let swept = step.coupled.motion.swept_volume_change_m3;
        network.finish(self, next, thermal, step, body_energy, loss, exterior_pressure_pa,
            swept, gas_heat_j + body.internal_heat_j, cx)
    }

    /// The same proposal through the source-bound gas/body conduction adapter.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn reactive_thermal_network_layers_candidate(&self, network: ModalGasNetwork<'_>,
        reaction: &GasReaction, gas_heat_j: f64, exterior_pressure_pa: f64,
        body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
        controls: CoupledReactiveModalThermalSolve, layers: MaterialConductionLayers<'_>, cx: &Cx<'_>)
        -> Result<ModalGasNetworkCandidate<CoupledReactiveModalThermalStep,
            UniformEnthalpyStep<ResolvedMaterialStatePoint>>,
            ModalGasNetworkError<ReactiveModalThermalError<LayerConductionError>>> {
        self.reactive_thermal_network_candidate(network, reaction, gas_heat_j, exterior_pressure_pa,
            body, chemistry, controls, |gas, thermal| layers.resolve(gas, thermal), cx)
    }
}
