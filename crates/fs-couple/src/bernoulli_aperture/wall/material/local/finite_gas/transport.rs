//! One native material transport snapshot drives multiple continuing acoustic owners.
//!
//! Each duct remains a separate, uniform, finite gas volume. No common-pressure
//! axial field, acoustic continuity across a restriction, or jet force is invented.
//! Original nozzle, diffusion, porous and membrane laws run once before any endpoint
//! work feedback. Every acoustic owner retains its native state and source wall.

use super::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgSolve, FiniteGasDgStep,
    MaterialWallDg, RegionalMaterialWallDg, stored_change,
};
use fs_duct::dg::{DgDuct, DgPropagator, DgState, PortInput};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve, GasHeatStep};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    snapshot::MAX_SNAPSHOT_PARCELS,
    species::{
        diffusion::network::DiffusionNetworkLink,
        network::transport::{self as material, GasTransportError, GasTransportPaths, GasTransportSolve, GasTransportStep},
    },
};
use std::collections::BTreeMap;

/// Retain either existing source-wall layout without rebuilding a live operator.
#[derive(Clone, Debug)]
pub enum FiniteGasDgBinding {
    Uniform(MaterialWallDg),
    Regional(RegionalMaterialWallDg),
}
impl FiniteGasDgBinding {
    /// Original native operator, with its complete material and gas sources.
    pub fn duct(&self) -> &DgDuct {
        match self { Self::Uniform(wall) => wall.duct(), Self::Regional(wall) => wall.duct() }
    }
}

/// The actual source binding, vibration history and certified stepping owner.
/// Its map key is the SAME original parcel identity used by material transport;
/// the gas itself stays exclusively in the original parcel map.
#[derive(Clone, Debug)]
pub struct FiniteGasDgNode {
    pub binding: FiniteGasDgBinding,
    pub state: DgState,
    pub propagator: DgPropagator,
}

/// Exactly one source/control input per gas parcel, including finite reservoirs.
#[derive(Clone, Debug)]
pub struct FiniteGasDgTransportInput {
    /// One original caloric solver per volume, independent of flow direction.
    pub gas: GasEnergySolve,
    /// External heat during this interval [J], excluding transport and wave work.
    pub supplied_heat_j: f64,
    /// Required exactly for parcels with a bound acoustic owner.
    pub acoustic: Option<FiniteGasDgSolve>,
    /// Original held physical volume-flow ports; reservoirs require an empty list.
    pub ports: Vec<PortInput>,
}

/// Original gross-species transport admission and one enclosing energy tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasDgTransportSolve {
    pub transport: GasTransportSolve,
    /// Absolute sum-of-endpoint energy residual [J]; no relative storage allowance.
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum FiniteGasDgTransportError {
    Invalid(&'static str),
    Transport(GasTransportError),
    Acoustic { parcel: u64, error: FiniteGasDgError<GasEnergyError> },
    Gas { parcel: u64, error: GasEnergyError },
    Energy { residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for FiniteGasDgTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas acoustic transport: {self:?}")
    }
}
impl std::error::Error for FiniteGasDgTransportError {}

/// Actual original-owner reports and independently measured enclosure storage.
#[derive(Debug)]
pub struct FiniteGasDgTransportStep {
    pub duration_s: f64,
    /// The single original-snapshot species/partial-enthalpy transport.
    /// Its entropy account belongs to that stage, not subsequent heat/wave work.
    pub transport: GasTransportStep,
    /// Same-clock acoustic endpoints and local gas/work closure, by parcel ID.
    pub acoustics: BTreeMap<u64, FiniteGasDgStep<GasCompositionStep>>,
    /// Original heating reports for parcels without an acoustic owner.
    pub reservoirs: BTreeMap<u64, GasHeatStep>,
    pub gas_storage_change_j: f64,
    pub acoustic_storage_change_j: f64,
    pub supplied_heat_j: f64,
    pub port_work_j: f64,
    /// Unpartitioned native gas/wall/boundary loss leaves the modeled enclosure.
    /// No loss is guessed to heat a particular parcel or wall.
    pub acoustic_dissipation_j: f64,
    /// Internal gas-to-wave transfer, already included in the storage changes.
    pub parameter_work_j: f64,
    pub total_mass_residual_kg: f64,
    /// Gas + wave/wall storage + native loss - external heat - port work.
    /// Internal enthalpy transfer and gas-to-wave work are not added again.
    pub combined_energy_residual_j: f64,
}

/// Complete immutable preview of all gas and native acoustic endpoints.
/// No accepted source binding, acoustic history or gas parcel changes until the
/// corresponding advance API publishes this candidate's complete endpoint set.
#[derive(Debug)]
pub struct FiniteGasDgTransportCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    nodes: BTreeMap<u64, FiniteGasDgNode>,
    step: FiniteGasDgTransportStep,
}
impl FiniteGasDgTransportCandidate {
    /// Proposed original gas parcels, keyed by their existing identities.
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.parcels }
    /// Proposed original source bindings, states and propagators.
    pub fn nodes(&self) -> &BTreeMap<u64, FiniteGasDgNode> { &self.nodes }
    /// Original-owner reports and complete energy/mass account.
    pub fn step(&self) -> &FiniteGasDgTransportStep { &self.step }
    /// Consume the complete preview without advancing any external owner.
    pub fn into_parts(self) -> (
        BTreeMap<u64, ClosedGasParcel>, BTreeMap<u64, FiniteGasDgNode>, FiniteGasDgTransportStep,
    ) { (self.parcels, self.nodes, self.step) }
}

use FiniteGasDgTransportError as Error;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn transport_error(error: GasTransportError) -> Error {
    match error { GasTransportError::Cancelled => Error::Cancelled, error => Error::Transport(error) }
}
fn gas_error(parcel: u64, error: GasEnergyError) -> Error {
    match error { GasEnergyError::Cancelled => Error::Cancelled, error => Error::Gas { parcel, error } }
}
fn acoustic_error(parcel: u64, error: FiniteGasDgError<GasEnergyError>) -> Error {
    match error {
        FiniteGasDgError::Cancelled | FiniteGasDgError::Background(GasEnergyError::Cancelled) => Error::Cancelled,
        error => Error::Acoustic { parcel, error },
    }
}
#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        if !value.is_finite() || !adjusted.is_finite() || !next.is_finite() {
            return Err(Error::Invalid("extensive acoustic network account is unrepresentable"));
        }
        self.correction = (next-self.value)-adjusted;
        self.value = next;
        Ok(())
    }
}

/// Advance original nozzle/diffusion transport and every bound acoustic owner.
///
/// All rates use ONE unchanged gas snapshot over the native acoustic duration.
/// The original combined transport owner admits source domains, all outgoing
/// species jointly, partial enthalpies, physical direction and entropy. New
/// arrivals cannot replenish an outgoing allowance or leave through another
/// connection in this explicit stage. Transport is never rerun by work trials.
///
/// Each duct then uses the existing finite-gas acoustic transaction: its native
/// sample starts from the original medium, and its scalar parameter-work feedback
/// reconstructs the transported parcel with external heat minus gas-to-wave work.
/// Unbound reservoirs use the original caloric heat owner. Every active path's
/// original driving-force gate is reapplied after ALL endpoint heating/work.
/// Final-source, history, clock, work, energy or cancellation refusal publishes
/// none of the parcel map, source bindings, states or propagators.
///
/// This is a first-order split of well-mixed fixed-volume chambers and native
/// acoustic samples. Ducts may have different pressures/compositions; they do not
/// become a continuous spatial acoustic field. Restrictions carry mean species
/// and enthalpy only, with no inferred acoustic scattering, jet momentum/noise,
/// static pressure prestress, wall heating or gas-loss thermal partition.
/// Diffusion retains its original near-isothermal/isobaric approximation and
/// separately authored paths; later endpoint changes may require a smaller next
/// interval. Source and timestep-refinement limits of all original owners remain.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_dg_transport<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, FiniteGasDgNode>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    nozzles: &[SnapshotNozzleLink],
    diffusion: &[DiffusionNetworkLink<'_>],
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgTransportStep, Error> {
    advance_finite_gas_dg_transport_with_paths(parcels, nodes, inputs,
        GasTransportPaths { nozzles, diffusion, ..Default::default() },
        duration_s, controls, cancelled)
}

/// Atomically advance all declared gas mechanisms and the existing acoustic
/// owners. Porous and membrane paths use the same snapshot, species inventory,
/// partial-enthalpy exchange and endpoint gates as nozzle/diffusion paths.
#[allow(clippy::too_many_arguments)]
pub fn advance_finite_gas_dg_transport_with_paths<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, FiniteGasDgNode>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgTransportStep, Error> {
    let candidate = finite_gas_dg_transport_candidate_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled)?;
    check(cancelled)?;
    let (next_parcels, next_nodes, report) = candidate.into_parts();
    *parcels = next_parcels;
    *nodes = next_nodes;
    Ok(report)
}

/// Preview the complete transport, heat and acoustic-work transaction without
/// advancing any accepted parcel, source binding, state or propagator.
///
/// Source-backed porous and membrane rates enter the same original-snapshot
/// transport stage exactly once. Every final total/partial-pressure ordering is
/// checked again after heating and acoustic parameter work. The existing
/// first-order split, fixed-volume acoustic and source-domain limits apply.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_dg_transport_candidate_with_paths<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, FiniteGasDgNode>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgTransportCandidate, Error> {
    check(cancelled)?;
    if parcels.is_empty() || parcels.len() > MAX_SNAPSHOT_PARCELS || nodes.is_empty()
        || !parcels.keys().eq(inputs.keys()) || nodes.keys().any(|id| !parcels.contains_key(id))
        || !duration_s.is_finite() || duration_s <= 0.0
        || !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j <= 0.0
    {
        return Err(Error::Invalid("bounded gas graph, complete inputs, acoustic nodes and positive clock/tolerance required"));
    }
    let epoch = nodes.values().next().expect("nonempty acoustic map").state.step_count();
    let mut endpoints = BTreeMap::new();
    for (&id, input) in inputs {
        check(cancelled)?;
        if !input.supplied_heat_j.is_finite()
            || input.acoustic.is_some() != nodes.contains_key(&id)
            || (!nodes.contains_key(&id) && !input.ports.is_empty())
        {
            return Err(Error::Invalid("finite heat and exact per-parcel acoustic input binding required"));
        }
        input.gas.validate().map_err(|error| gas_error(id, error))?;
        endpoints.insert(id, input.gas);
    }
    for (&id, node) in nodes.iter() {
        check(cancelled)?;
        if node.propagator.dt() != duration_s || node.state.step_count() != epoch
            || node.binding.duct().gas_field().and_then(|field| field.uniform_state())
                != Some(parcels[&id].state())
        {
            return Err(Error::Invalid("all acoustic owners require one clock and their complete original finite gas"));
        }
        let mut sample = node.state.clone();
        node.propagator.sync_nodal(&mut sample)
            .map_err(|error| acoustic_error(id, super::native(error)))?;
        node.binding.duct().port_pressures(&sample)
            .map_err(|error| acoustic_error(id, super::native(error)))?;
    }
    let (transported, transport) = material::step_candidate_with_paths(
        parcels, paths, duration_s, &endpoints, controls.transport, &mut *cancelled,
    ).map_err(transport_error)?.into_parts();
    let mut next_parcels = transported.clone();
    let mut next_nodes = nodes.clone();
    let mut acoustics = BTreeMap::new();
    let mut reservoirs = BTreeMap::new();
    let mut heat = Sum::default();
    let mut port_work = Sum::default();
    let mut losses = Sum::default();
    let mut parameter_work = Sum::default();
    for (&id, input) in inputs {
        check(cancelled)?;
        heat.add(input.supplied_heat_j)?;
        if let Some(node) = next_nodes.get_mut(&id) {
            let transported_gas = &transported[&id];
            let enthalpy = transport.exchange.parcels.iter().find(|row| row.id == id)
                .ok_or(Error::Invalid("original transport lost a parcel account"))?.net_enthalpy_j;
            // Start the acoustic transaction with the ORIGINAL source snapshot.
            // Its repeatable callback reads the ONE immutable transported inventory.
            let mut gas = parcels[&id].clone();
            let evolve = |_: &ClosedGasParcel, _: f64, work: f64, cancel: &mut C| {
                let mut endpoint = transported_gas.clone();
                let source = endpoint.replace_model_with_heat_work_volume(
                    transported_gas.model().clone(), transported_gas.volume_m3(),
                    input.supplied_heat_j, work, input.gas, cancel,
                )?;
                Ok::<_, GasEnergyError>(FiniteGasDgBackground {
                    gas: endpoint, supplied_heat_j: input.supplied_heat_j,
                    transported_enthalpy_j: enthalpy, source,
                })
            };
            let solve = input.acoustic.expect("validated acoustic input");
            let step = match &mut node.binding {
                FiniteGasDgBinding::Uniform(wall) => wall.step_with_finite_gas_controlled(
                    &mut node.state, &mut node.propagator, &mut gas, &input.ports, solve, evolve, cancelled,
                ),
                FiniteGasDgBinding::Regional(wall) => wall.step_with_finite_gas_controlled(
                    &mut node.state, &mut node.propagator, &mut gas, &input.ports, solve, evolve, cancelled,
                ),
            }.map_err(|error| acoustic_error(id, error))?;
            port_work.add(step.acoustic.ledger.port_work_step_j)?;
            losses.add(step.acoustic_dissipation_j)?;
            parameter_work.add(step.refresh.parameter_work_j)?;
            next_parcels.insert(id, gas);
            acoustics.insert(id, step);
        } else {
            let gas = next_parcels.get_mut(&id).expect("original complete transport map");
            let step = gas.add_heat(input.supplied_heat_j, input.gas, &mut *cancelled)
                .map_err(|error| gas_error(id, error))?;
            reservoirs.insert(id, step);
        }
    }
    // Heating and wave work can change a pressure ordering that transport alone
    // admitted. Reuse the original rule against the complete FINAL endpoint map.
    transport.admit_endpoint_driving_forces(&next_parcels, &mut *cancelled)
        .map_err(transport_error)?;
    let mut gas_storage = Sum::default();
    let mut acoustic_storage = Sum::default();
    let mut mass = Sum::default();
    for (&id, gas) in &next_parcels {
        check(cancelled)?;
        gas_storage.add(stored_change(&parcels[&id], gas))?;
        mass.add(gas.mass_kg()-parcels[&id].mass_kg())?;
    }
    for (&id, node) in &next_nodes {
        check(cancelled)?;
        acoustic_storage.add(node.state.stored_energy_j()-nodes[&id].state.stored_energy_j())?;
    }
    let residual = ((gas_storage.value+acoustic_storage.value)+losses.value)
        - heat.value - port_work.value;
    if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j {
        return Err(Error::Energy { residual_j: residual, tolerance_j: controls.energy_tolerance_j });
    }
    if mass.value.abs() > controls.transport.exchange.mass_tolerance_kg {
        return Err(Error::Invalid("actual final gas network mass balance failed"));
    }
    let report = FiniteGasDgTransportStep {
        duration_s, transport, acoustics, reservoirs,
        gas_storage_change_j: gas_storage.value,
        acoustic_storage_change_j: acoustic_storage.value,
        supplied_heat_j: heat.value, port_work_j: port_work.value,
        acoustic_dissipation_j: losses.value, parameter_work_j: parameter_work.value,
        total_mass_residual_kg: mass.value, combined_energy_residual_j: residual,
    };
    check(cancelled)?;
    Ok(FiniteGasDgTransportCandidate { parcels: next_parcels, nodes: next_nodes, step: report })
}

