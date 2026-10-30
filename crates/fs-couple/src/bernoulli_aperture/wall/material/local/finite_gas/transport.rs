//! One native material transport snapshot drives multiple continuing acoustic owners.
//!
//! Each acoustic owner remains a separate, uniform, finite gas volume. No common-
//! pressure axial field, acoustic continuity across a restriction, or jet force
//! is invented. Complete port loads use this same transaction through their
//! original node adapter, including their represented channel/branch storage.
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
    /// A complete port load requires exactly its inlet, with flow in `effort`.
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
pub enum FiniteGasDgTransportError<E = GasEnergyError> {
    Invalid(&'static str),
    Transport(GasTransportError),
    Acoustic { parcel: u64, error: FiniteGasDgError<E> },
    Gas { parcel: u64, error: E },
    Energy { residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteGasDgTransportError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas acoustic transport: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FiniteGasDgTransportError<E> {}

/// Actual original-owner reports and independently measured enclosure storage.
#[derive(Debug)]
pub struct FiniteGasDgTransportStep<A = FiniteGasDgStep<GasCompositionStep>, R = GasHeatStep> {
    pub duration_s: f64,
    /// The single original-snapshot species/partial-enthalpy transport.
    /// Its entropy account belongs to that stage, not subsequent heat/wave work.
    pub transport: GasTransportStep,
    /// Same-clock bound-node endpoints and local gas/work closure, by parcel ID.
    /// Moving-source adapters include their original mechanical reports here.
    pub acoustics: BTreeMap<u64, A>,
    /// Original endpoint reports for parcels without an acoustic owner.
    pub reservoirs: BTreeMap<u64, R>,
    pub gas_storage_change_j: f64,
    /// Complete bound-node storage, including moving-source mechanics when present.
    pub acoustic_storage_change_j: f64,
    /// Actual heat into the gas boundary. A thermal companion's reciprocal
    /// transfer is included here and cancels only in its enclosing body account.
    pub supplied_heat_j: f64,
    /// Net external work into bound nodes, including mechanical drives and
    /// subtracting signed exterior exchange. Internal source/load work cancels.
    pub port_work_j: f64,
    /// Unpartitioned native mechanical/gas/wall/boundary loss leaves the enclosure.
    /// No loss is guessed to heat a particular parcel or wall.
    pub acoustic_dissipation_j: f64,
    /// Internal gas-to-wave transfer, already included in the storage changes.
    pub parameter_work_j: f64,
    pub total_mass_residual_kg: f64,
    /// Gas + complete bound-node storage + native loss - supplied gas heat - port work.
    /// Internal enthalpy transfer and gas-to-wave work are not added again.
    pub combined_energy_residual_j: f64,
}

/// Complete immutable preview of all gas and native acoustic endpoints.
/// No accepted source binding, acoustic history or gas parcel changes until the
/// corresponding advance API publishes this candidate's complete endpoint set.
#[derive(Debug)]
pub struct FiniteGasDgTransportCandidate<
    N = FiniteGasDgNode,
    A = FiniteGasDgStep<GasCompositionStep>,
    R = GasHeatStep,
> {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    nodes: BTreeMap<u64, N>,
    step: FiniteGasDgTransportStep<A, R>,
}
impl<N, A, R> FiniteGasDgTransportCandidate<N, A, R> {
    /// Proposed original gas parcels, keyed by their existing identities.
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.parcels }
    /// Proposed original source bindings, states and propagators.
    pub fn nodes(&self) -> &BTreeMap<u64, N> { &self.nodes }
    /// Original-owner reports and complete energy/mass account.
    pub fn step(&self) -> &FiniteGasDgTransportStep<A, R> { &self.step }
    /// Consume the complete preview without advancing any external owner.
    pub fn into_parts(self) -> (
        BTreeMap<u64, ClosedGasParcel>, BTreeMap<u64, N>, FiniteGasDgTransportStep<A, R>,
    ) { (self.parcels, self.nodes, self.step) }
}

/// Only bound-node ownership and its typed physical source inputs vary across
/// wall, complete-load and moving-source adapters.
/// Native transport, caloric feedback input, endpoint gates, global accounts
/// and atomic publication stay in the single transaction below.
pub(crate) trait FiniteGasTransportNode: Clone {
    type Step<S>;
    type Input;

    /// Admit the node-specific source controls before any transport runs, and
    /// expose the same caloric controls for bound nodes and finite reservoirs.
    fn transport_input(
        input: &Self::Input,
        bound: bool,
    ) -> Result<&FiniteGasDgTransportInput, &'static str>;

    fn dt(&self) -> f64;
    fn step_count(&self) -> u64;
    fn stored_energy_j(&self) -> f64;

    fn admit<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &Self::Input,
        cancelled: &mut C,
    ) -> Result<(), FiniteGasDgError<E>>;

    fn advance<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        input: &Self::Input,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C)
            -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<Self::Step<S>, FiniteGasDgError<E>>;

    /// Net port work, exported acoustic loss and internal gas-to-wave work [J].
    /// Net port work already accounts for any represented exterior terminal.
    fn accounts<S>(step: &Self::Step<S>) -> [f64; 3];

    /// The actual accepted background, never a side effect of a rejected trial.
    fn source<S>(step: &Self::Step<S>) -> &S;
}

impl FiniteGasTransportNode for FiniteGasDgNode {
    type Step<S> = FiniteGasDgStep<S>;
    type Input = FiniteGasDgTransportInput;

    fn transport_input(
        input: &Self::Input,
        _bound: bool,
    ) -> Result<&FiniteGasDgTransportInput, &'static str> {
        Ok(input)
    }

    fn dt(&self) -> f64 { self.propagator.dt() }
    fn step_count(&self) -> u64 { self.state.step_count() }
    fn stored_energy_j(&self) -> f64 { self.state.stored_energy_j() }

    fn admit<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        cancelled: &mut C,
    ) -> Result<(), FiniteGasDgError<E>> {
        super::check(cancelled)?;
        super::admit_controls(input.acoustic.expect("validated acoustic input"))?;
        let duct = self.binding.duct();
        if duct.gas_field().and_then(|field| field.uniform_state()) != Some(gas.state()) {
            return Err(FiniteGasDgError::Invalid(
                "acoustic owner and finite gas require one complete original source",
            ));
        }
        let radius = duct.uniform_circular_radius_m().ok_or(FiniteGasDgError::Invalid(
            "finite background requires the original circular fixed geometry",
        ))?;
        let volume = core::f64::consts::PI * radius * radius * duct.total_length();
        if !volume.is_finite() || volume <= 0.0
            || (gas.volume_m3() / volume - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(FiniteGasDgError::Invalid(
                "finite background volume must equal the actual duct volume",
            ));
        }
        let mut sample = self.state.clone();
        self.propagator.sync_nodal(&mut sample).map_err(super::native)?;
        duct.port_pressures(&sample).map_err(super::native)?;
        super::check(cancelled)
    }

    fn advance<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C)
            -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<Self::Step<S>, FiniteGasDgError<E>> {
        let solve = input.acoustic.expect("validated acoustic input");
        match &mut self.binding {
            FiniteGasDgBinding::Uniform(wall) => wall.step_with_finite_gas_controlled(
                &mut self.state, &mut self.propagator, gas, &input.ports,
                solve, evolve, cancelled,
            ),
            FiniteGasDgBinding::Regional(wall) => wall.step_with_finite_gas_controlled(
                &mut self.state, &mut self.propagator, gas, &input.ports,
                solve, evolve, cancelled,
            ),
        }
    }

    fn accounts<S>(step: &Self::Step<S>) -> [f64; 3] {
        [step.acoustic.ledger.port_work_step_j, step.acoustic_dissipation_j,
            step.refresh.parameter_work_j]
    }

    fn source<S>(step: &Self::Step<S>) -> &S { &step.background }
}

/// Only the post-transport caloric/chemical/thermal endpoint varies. The transaction
/// retains ownership of transport, native work feedback and complete accounts.
pub(crate) trait FiniteGasTransportEvolution {
    type Source;
    type Reservoir;
    type Error: From<GasEnergyError>;

    fn is_cancelled(error: &Self::Error) -> bool;

    fn admit<C: FnMut() -> bool>(
        &self, id: u64, gas: &ClosedGasParcel, input: &FiniteGasDgTransportInput,
        duration_s: f64, cancelled: &mut C,
    ) -> Result<(), Self::Error>;

    fn node_candidate<C: FnMut() -> bool>(
        &self, id: u64, gas: &ClosedGasParcel, input: &FiniteGasDgTransportInput,
        duration_s: f64, work_by_gas_j: f64, cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Source), Self::Error>;

    fn reservoir_candidate<C: FnMut() -> bool>(
        &self, id: u64, gas: &ClosedGasParcel, input: &FiniteGasDgTransportInput,
        duration_s: f64, cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Reservoir), Self::Error>;

    /// Actual heat crossing the gas boundary at this endpoint. A finite thermal
    /// companion includes its reciprocal Q here; its enclosing owner must then
    /// account for the companion's storage and cancel Q as an internal transfer.
    /// Mechanical work and transported enthalpy retain their separate accounts.
    fn node_supplied_heat_j(
        &self, input: &FiniteGasDgTransportInput, _: &Self::Source,
    ) -> f64 { input.supplied_heat_j }

    fn reservoir_supplied_heat_j(
        &self, input: &FiniteGasDgTransportInput, _: &Self::Reservoir,
    ) -> f64 { input.supplied_heat_j }

    /// Apply physical endpoint gates after the native work residual has closed.
    fn admit_node_endpoint(&self, id: u64, source: &Self::Source) -> Result<(), Self::Error>;
}

struct CaloricEvolution;
impl FiniteGasTransportEvolution for CaloricEvolution {
    type Source = GasCompositionStep;
    type Reservoir = GasHeatStep;
    type Error = GasEnergyError;

    fn is_cancelled(error: &Self::Error) -> bool { *error == GasEnergyError::Cancelled }

    fn admit<C: FnMut() -> bool>(
        &self, _: u64, _: &ClosedGasParcel, _: &FiniteGasDgTransportInput,
        _: f64, _: &mut C,
    ) -> Result<(), Self::Error> { Ok(()) }

    fn node_candidate<C: FnMut() -> bool>(
        &self, _: u64, gas: &ClosedGasParcel, input: &FiniteGasDgTransportInput,
        _: f64, work_by_gas_j: f64, cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Source), Self::Error> {
        let mut endpoint = gas.clone();
        let source = endpoint.replace_model_with_heat_work_volume(
            gas.model().clone(), gas.volume_m3(), input.supplied_heat_j,
            work_by_gas_j, input.gas, cancelled,
        )?;
        Ok((endpoint, source))
    }

    fn reservoir_candidate<C: FnMut() -> bool>(
        &self, _: u64, gas: &ClosedGasParcel, input: &FiniteGasDgTransportInput,
        _: f64, cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Reservoir), Self::Error> {
        let mut endpoint = gas.clone();
        let source = endpoint.add_heat(input.supplied_heat_j, input.gas, cancelled)?;
        Ok((endpoint, source))
    }

    fn admit_node_endpoint(&self, _: u64, _: &Self::Source) -> Result<(), Self::Error> {
        Ok(())
    }
}

use FiniteGasDgTransportError as Error;
fn check<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error<E>> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn transport_error<E>(error: GasTransportError) -> Error<E> {
    match error { GasTransportError::Cancelled => Error::Cancelled, error => Error::Transport(error) }
}
fn gas_error<T: FiniteGasTransportEvolution>(parcel: u64, error: T::Error) -> Error<T::Error> {
    if T::is_cancelled(&error) { Error::Cancelled } else { Error::Gas { parcel, error } }
}
fn acoustic_error<T: FiniteGasTransportEvolution>(
    parcel: u64, error: FiniteGasDgError<T::Error>,
) -> Error<T::Error> {
    match error {
        FiniteGasDgError::Cancelled => Error::Cancelled,
        FiniteGasDgError::Background(error) if T::is_cancelled(&error) => Error::Cancelled,
        error => Error::Acoustic { parcel, error },
    }
}
#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), Error<E>> {
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
    advance_finite_gas_transport_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled)
}

/// Publish one complete original-owner transaction for either acoustic layout.
#[allow(clippy::too_many_arguments)]
pub(crate) fn advance_finite_gas_transport_with_paths<N: FiniteGasTransportNode, C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, N>,
    inputs: &BTreeMap<u64, N::Input>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgTransportStep<N::Step<GasCompositionStep>>, Error> {
    let candidate = finite_gas_transport_candidate_with_paths(
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
    finite_gas_transport_candidate_with_paths(
        parcels, nodes, inputs, paths, duration_s, controls, cancelled)
}

/// Shared physical transport/accounting kernel. The node adapter can only
/// advance its original acoustic owner and report that owner's actual accounts.
#[allow(clippy::too_many_arguments)]
pub(crate) fn finite_gas_transport_candidate_with_paths<N: FiniteGasTransportNode, C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, N>,
    inputs: &BTreeMap<u64, N::Input>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    cancelled: &mut C,
) -> Result<FiniteGasDgTransportCandidate<N, N::Step<GasCompositionStep>>, Error> {
    finite_gas_transport_candidate_with_evolution(
        parcels, nodes, inputs, paths, duration_s, controls, &CaloricEvolution, cancelled,
    )
}

/// The single transport/work/accounting transaction with an original physical
/// endpoint owner. Chemistry and finite heat exchange change this endpoint; native
/// transport still runs once, and no trial can consume an accepted inventory.
#[allow(clippy::too_many_arguments)]
pub(crate) fn finite_gas_transport_candidate_with_evolution<
    N: FiniteGasTransportNode, T: FiniteGasTransportEvolution, C: FnMut() -> bool,
>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, N>,
    inputs: &BTreeMap<u64, N::Input>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: FiniteGasDgTransportSolve,
    evolution: &T,
    cancelled: &mut C,
) -> Result<FiniteGasDgTransportCandidate<N, N::Step<T::Source>, T::Reservoir>, Error<T::Error>> {
    check(cancelled)?;
    if parcels.is_empty() || parcels.len() > MAX_SNAPSHOT_PARCELS || nodes.is_empty()
        || !parcels.keys().eq(inputs.keys()) || nodes.keys().any(|id| !parcels.contains_key(id))
        || !duration_s.is_finite() || duration_s <= 0.0
        || !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j <= 0.0
    {
        return Err(Error::Invalid("bounded gas graph, complete inputs, acoustic nodes and positive clock/tolerance required"));
    }
    let epoch = nodes.values().next().expect("nonempty acoustic map").step_count();
    let mut endpoints = BTreeMap::new();
    for (&id, input) in inputs {
        check(cancelled)?;
        let input = N::transport_input(input, nodes.contains_key(&id)).map_err(Error::Invalid)?;
        if !input.supplied_heat_j.is_finite()
            || input.acoustic.is_some() != nodes.contains_key(&id)
            || (!nodes.contains_key(&id) && !input.ports.is_empty())
        {
            return Err(Error::Invalid("finite heat and exact per-parcel acoustic input binding required"));
        }
        input.gas.validate().map_err(|error| gas_error::<T>(id, error.into()))?;
        evolution.admit(id, &parcels[&id], input, duration_s, cancelled)
            .map_err(|error| gas_error::<T>(id, error))?;
        endpoints.insert(id, input.gas);
    }
    for (&id, node) in nodes.iter() {
        check(cancelled)?;
        if node.dt() != duration_s || node.step_count() != epoch {
            return Err(Error::Invalid("all acoustic owners require one clock"));
        }
        node.admit(&parcels[&id], &inputs[&id], cancelled)
            .map_err(|error| acoustic_error::<T>(id, error))?;
    }
    let (transported, transport) = material::step_candidate_with_paths(
        parcels, paths, duration_s, &endpoints, controls.transport, &mut *cancelled,
    ).map_err(transport_error)?.into_parts();
    let mut next_parcels = BTreeMap::new();
    let mut next_nodes = nodes.clone();
    let mut acoustics = BTreeMap::new();
    let mut reservoirs = BTreeMap::new();
    let mut heat = Sum::default();
    let mut port_work = Sum::default();
    let mut losses = Sum::default();
    let mut parameter_work = Sum::default();
    for (&id, node_input) in inputs {
        check(cancelled)?;
        let input = N::transport_input(node_input, nodes.contains_key(&id)).map_err(Error::Invalid)?;
        if let Some(node) = next_nodes.get_mut(&id) {
            let transported_gas = &transported[&id];
            let enthalpy = transport.exchange.parcels.iter().find(|row| row.id == id)
                .ok_or(Error::Invalid("original transport lost a parcel account"))?.net_enthalpy_j;
            // Start the acoustic transaction with the ORIGINAL source snapshot.
            // Its repeatable callback reads the ONE immutable transported inventory.
            let mut gas = parcels[&id].clone();
            let evolve = |_: &ClosedGasParcel, _: f64, work: f64, cancel: &mut C| {
                let (endpoint, source) = evolution.node_candidate(
                    id, transported_gas, input, duration_s, work, cancel,
                )?;
                let supplied_heat_j = evolution.node_supplied_heat_j(input, &source);
                Ok::<_, T::Error>(FiniteGasDgBackground {
                    gas: endpoint, supplied_heat_j,
                    transported_enthalpy_j: enthalpy, source,
                })
            };
            let step = node.advance(&mut gas, node_input, evolve, cancelled)
                .map_err(|error| acoustic_error::<T>(id, error))?;
            evolution.admit_node_endpoint(id, N::source::<T::Source>(&step))
                .map_err(|error| gas_error::<T>(id, error))?;
            heat.add(evolution.node_supplied_heat_j(input, N::source::<T::Source>(&step)))?;
            let [net_port_work, acoustic_loss, gas_work] = N::accounts::<T::Source>(&step);
            port_work.add(net_port_work)?;
            losses.add(acoustic_loss)?;
            parameter_work.add(gas_work)?;
            next_parcels.insert(id, gas);
            acoustics.insert(id, step);
        } else {
            let (gas, step) = evolution.reservoir_candidate(
                id, &transported[&id], input, duration_s, cancelled,
            ).map_err(|error| gas_error::<T>(id, error))?;
            heat.add(evolution.reservoir_supplied_heat_j(input, &step))?;
            next_parcels.insert(id, gas);
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
        acoustic_storage.add(node.stored_energy_j()-nodes[&id].stored_energy_j())?;
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

