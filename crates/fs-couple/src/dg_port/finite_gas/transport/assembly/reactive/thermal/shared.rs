//! One finite material body exchanges heat with several original gas volumes.
//!
//! Material transport is prepared once. The existing scalar feedback owner
//! closes one body enthalpy endpoint against every signed gas contact; each
//! contact reuses the original reacting gas and native acoustic work solve.
//! No thermal mass is copied per contact and no accepted history advances in a
//! discarded body, heat or work trial.

use super::{
    BoundaryGasReactionStep, ClosedGasParcel, EndpointError, EquilibriumPhaseState, Error,
    Evolution, FiniteGasAcousticInput, FiniteGasAcousticNode, FiniteGasAcousticReaction,
    FiniteGasAcousticReactionReservoir, FiniteGasAcousticReactionSource, FiniteGasAcousticStep,
    FiniteGasDgTransportInput, FiniteGasDgTransportStep, FiniteGasTransportEvolution,
    FiniteGasTransportNode, GasReactionError, GasTransportPaths, HeatExchangeSolve,
    MAX_SNAPSHOT_PARCELS, MaterialConductionLayers, MaterialPropertySelection,
    ReactiveFiniteGasAcousticTransportSolve, ResolvedMaterialStatePoint, Sum,
    UniformEnthalpyStepInput, feedback, feedback_error as endpoint_feedback_error, poll, positive,
    species_residual,
};
use crate::bernoulli_aperture::wall::material::prepare_finite_gas_transport;
use std::cell::Cell;
use std::collections::BTreeMap;

/// Hard work limit for all gas endpoint evaluations in one shared-body step.
pub const MAX_SHARED_THERMAL_ENDPOINT_EVALUATIONS: usize = 1_000_000;

/// One disjoint patch of the shared body's surface, keyed by its gas parcel.
/// The original material layer supplies both conductivities and source gates.
pub struct SharedFiniteGasAcousticThermalPort<'a> {
    pub layers: MaterialConductionLayers<'a>,
    /// Signed heat into this gas, independently bracketed from the body's net.
    pub exchange: HeatExchangeSolve,
}

/// Existing physical controls and bounded work for the complete nested solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticThermalSolve {
    pub reactive: ReactiveFiniteGasAcousticTransportSolve,
    /// Bracket for NET heat out of the body. Individual contacts may exchange
    /// larger opposite heats even when this coordinate is exactly zero.
    pub body_exchange: HeatExchangeSolve,
    /// Every original caloric/chemical endpoint, including discarded heat and
    /// acoustic-work probes and parcels without a thermal contact, costs one.
    pub maximum_endpoint_evaluations: usize,
}

/// The exact material evaluation and signed heat of one accepted gas contact.
/// Body storage is intentionally absent: it belongs to the shared report once.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticThermalStep {
    pub material: ResolvedMaterialStatePoint,
    pub body_temperature_k: f64,
    pub gas_temperature_k: f64,
    pub conductance_w_k: f64,
    pub solid_to_gas_heat_j: f64,
    /// Q_i - dt G_i(T_body,end - T_gas,end) [J].
    pub source_residual_j: f64,
    pub gas_energy_residual_j: f64,
    /// Contact evaluations within this accepted acoustic-work probe.
    pub exchange_trials: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticThermalSource {
    pub native: FiniteGasAcousticReactionSource,
    pub port: Option<SharedFiniteGasAcousticThermalStep>,
    /// Actual external gas heat plus this contact's Q; work stays separate.
    pub supplied_heat_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SharedFiniteGasAcousticThermalReservoir {
    pub native: FiniteGasAcousticReactionReservoir,
    pub port: Option<SharedFiniteGasAcousticThermalStep>,
    pub supplied_heat_j: f64,
}

/// One material body and all original native reports at the same endpoint.
#[derive(Debug)]
pub struct SharedThermalReactiveFiniteGasAcousticTransportStep {
    /// Gas/acoustic boundary only: supplied heat includes every internal Q_i.
    pub transport: FiniteGasDgTransportStep<
        FiniteGasAcousticStep<SharedFiniteGasAcousticThermalSource>,
        SharedFiniteGasAcousticThermalReservoir,
    >,
    pub body_state: EquilibriumPhaseState,
    pub body_storage_change_j: f64,
    pub external_body_heat_j: f64,
    pub external_gas_heat_j: f64,
    /// Compensated sum of the actual signed gas contact heats.
    pub solid_to_gas_heat_j: f64,
    /// Boundary heat into the body: minus the outer root's net-heat coordinate.
    pub body_boundary_heat_j: f64,
    /// Actual body storage + sum(Q_i) - prescribed body heat [J].
    pub body_energy_residual_j: f64,
    /// Outer net-heat coordinate minus the independently evaluated sum(Q_i).
    pub source_residual_j: f64,
    /// Gas + acoustics + ONE body + native losses - external heat/work [J].
    pub combined_energy_residual_j: f64,
    pub maximum_species_residual_mol: f64,
    pub body_trials: usize,
    /// All native endpoint evaluations, including discarded body/Q/work probes.
    pub endpoint_evaluations: usize,
}

impl SharedThermalReactiveFiniteGasAcousticTransportStep {
    pub fn reaction(&self, id: u64) -> Option<&BoundaryGasReactionStep> {
        if let Some(step) = self.transport.acoustics.get(&id) {
            return match &FiniteGasAcousticNode::source::<SharedFiniteGasAcousticThermalSource>(
                step,
            )
            .native
            {
                FiniteGasAcousticReactionSource::Caloric(_) => None,
                FiniteGasAcousticReactionSource::Reaction(step) => Some(step),
            };
        }
        match &self.transport.reservoirs.get(&id)?.native {
            FiniteGasAcousticReactionReservoir::Caloric(_) => None,
            FiniteGasAcousticReactionReservoir::Reaction(step) => Some(step),
        }
    }

    pub fn port(&self, id: u64) -> Option<&SharedFiniteGasAcousticThermalStep> {
        if let Some(step) = self.transport.acoustics.get(&id) {
            return FiniteGasAcousticNode::source::<SharedFiniteGasAcousticThermalSource>(step)
                .port
                .as_ref();
        }
        self.transport.reservoirs.get(&id)?.port.as_ref()
    }
}

/// No state is published until the one shared body and every gas owner accept.
pub struct SharedThermalReactiveFiniteGasAcousticTransportCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    nodes: BTreeMap<u64, FiniteGasAcousticNode>,
    step: SharedThermalReactiveFiniteGasAcousticTransportStep,
}
impl SharedThermalReactiveFiniteGasAcousticTransportCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn nodes(&self) -> &BTreeMap<u64, FiniteGasAcousticNode> {
        &self.nodes
    }
    pub const fn body_state(&self) -> &EquilibriumPhaseState {
        &self.step.body_state
    }
    pub const fn step(&self) -> &SharedThermalReactiveFiniteGasAcousticTransportStep {
        &self.step
    }
    #[allow(clippy::type_complexity)]
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, FiniteGasAcousticNode>,
        EquilibriumPhaseState,
        SharedThermalReactiveFiniteGasAcousticTransportStep,
    ) {
        (self.parcels, self.nodes, self.step.body_state, self.step)
    }
}

fn outer_feedback_error(error: feedback::ScalarError<Error>) -> Error {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => Error::Invalid(message),
        feedback::ScalarError::Cancelled => Error::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => EndpointError::OutsideBracket {
            lower_residual_j: lower,
            upper_residual_j: upper,
        }
        .into(),
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => EndpointError::NotConverged {
            iterations,
            residual_j: residual,
        }
        .into(),
    }
}

fn prepare_ports(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    body_state: &EquilibriumPhaseState,
    body: UniformEnthalpyStepInput<'_>,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
    duration_s: f64,
    controls: SharedFiniteGasAcousticThermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BTreeMap<u64, MaterialPropertySelection>, Error> {
    poll(cancelled)?;
    if parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || ports.is_empty()
        || ports.len() > parcels.len()
        || !positive(body.mass_kg)
        || !positive(body.volume_m3)
        || !positive(body.surface_area_m2)
        || !positive(duration_s)
        || body.duration_s != duration_s
        || !body.internal_heat_j.is_finite()
        || body.initial != *body_state
        || !(1..=MAX_SHARED_THERMAL_ENDPOINT_EVALUATIONS)
            .contains(&controls.maximum_endpoint_evaluations)
    {
        return Err(Error::Invalid(
            "bounded nonempty gas contacts, one canonical finite body/clock and endpoint work budget required",
        ));
    }
    if body
        .curve
        .state_at_specific_enthalpy(body.initial.specific_enthalpy_j_kg())
        .map_err(EndpointError::from)?
        != body.initial
    {
        return Err(Error::Invalid(
            "accepted shared body must belong to its original enthalpy chart",
        ));
    }
    feedback::admit::<core::convert::Infallible>(controls.body_exchange).map_err(|_| {
        Error::Invalid("valid net-body heat bracket and physical energy allowances required")
    })?;
    let mut area = Sum::default();
    let mut pins = BTreeMap::new();
    for (&id, port) in ports {
        poll(cancelled)?;
        let gas = parcels.get(&id).ok_or(Error::Invalid(
            "shared-body contact requires an existing gas parcel",
        ))?;
        if body.curve.material_card_identity() != port.layers.material.content_hash()
            || !port.layers.area_m2.is_finite()
            || port.layers.area_m2 < 0.0
        {
            return Err(Error::Invalid(
                "all shared-body layers require the original material card and finite nonnegative area",
            ));
        }
        area.add(port.layers.area_m2)?;
        if area.value > body.surface_area_m2 {
            return Err(Error::Invalid(
                "shared-body contacts exceed its one exposed surface area",
            ));
        }
        feedback::admit::<core::convert::Infallible>(port.exchange).map_err(|_| {
            Error::Invalid("each shared-body contact requires valid heat/work controls")
        })?;
        // Even an explicitly insulated patch retains its original source gate.
        let material = port
            .layers
            .resolve(gas.state(), body.initial)
            .map_err(EndpointError::from)?
            .source;
        let selected = material
            .property_by_key(port.layers.conductivity)
            .ok_or(Error::Invalid(
                "shared-body contact lacks its exact conductivity source",
            ))?
            .answer()
            .receipt
            .selected;
        pins.insert(
            id,
            MaterialPropertySelection::PinnedByProperty(vec![(
                port.layers.conductivity.name().to_owned(),
                selected,
            )]),
        );
    }
    poll(cancelled)?;
    Ok(pins)
}

struct SharedEvolution<'a, 'b, 'r> {
    native: Evolution<'a, 'r>,
    body_state: EquilibriumPhaseState,
    ports: &'a BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'b>>,
    pins: &'a BTreeMap<u64, MaterialPropertySelection>,
    evaluations: &'a Cell<usize>,
    maximum_evaluations: usize,
}

struct Endpoint<S> {
    gas: ClosedGasParcel,
    native: S,
    port: Option<SharedFiniteGasAcousticThermalStep>,
    supplied_heat_j: f64,
}

impl SharedEvolution<'_, '_, '_> {
    fn charge_endpoint(&self, cancelled: &mut impl FnMut() -> bool) -> Result<(), EndpointError> {
        poll(cancelled)?;
        if self.evaluations.get() >= self.maximum_evaluations {
            return Err(EndpointError::Invalid(
                "shared finite-body native endpoint evaluation budget exhausted",
            ));
        }
        self.evaluations.set(self.evaluations.get() + 1);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn candidate<S, C: FnMut() -> bool>(
        &self,
        id: u64,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        duration_s: f64,
        work_by_gas_j: f64,
        cancelled: &mut C,
        mut evolve: impl FnMut(
            &FiniteGasDgTransportInput,
            &mut C,
        ) -> Result<(ClosedGasParcel, S), GasReactionError>,
    ) -> Result<Endpoint<S>, EndpointError> {
        let Some(port) = self.ports.get(&id) else {
            self.charge_endpoint(cancelled)?;
            let (gas, native) = evolve(input, cancelled)?;
            return Ok(Endpoint {
                gas,
                native,
                port: None,
                supplied_heat_j: input.supplied_heat_j,
            });
        };
        let layers = MaterialConductionLayers {
            material: port.layers.material,
            context: port.layers.context,
            temperature_axis: port.layers.temperature_axis,
            liquid_fraction_axis: port.layers.liquid_fraction_axis,
            conductivity: port.layers.conductivity,
            selection: self.pins[&id].clone(),
            area_m2: port.layers.area_m2,
            gas_distance_m: port.layers.gas_distance_m,
            body_distance_m: port.layers.body_distance_m,
        };
        let candidate = feedback::solve_scalar(
            port.exchange.heat_bracket_j, port.exchange.source_tolerance_j,
            port.exchange.maximum_iterations,
            |q, trials, cancelled| {
                let heat_j = input.supplied_heat_j + q;
                if !heat_j.is_finite()
                    || (q != 0.0 && heat_j == input.supplied_heat_j)
                    || (input.supplied_heat_j != 0.0 && heat_j == q)
                {
                    return Err(EndpointError::Invalid("shared-body transfer or external gas heat disappears in energy coordinates"));
                }
                let mut heated = input.clone();
                heated.supplied_heat_j = heat_j;
                self.charge_endpoint(cancelled)?;
                let (next, native) = evolve(&heated, cancelled)?;
                poll(cancelled)?;
                let sample = layers.resolve(next.state(), self.body_state)?;
                let gas_temperature_k = next.state().state().temperature;
                let body_temperature_k = self.body_state.temperature_k();
                let delta = body_temperature_k - gas_temperature_k;
                let g = sample.conductance_w_k;
                let power = g * delta;
                let evaluated_heat = duration_s * power;
                let residual = q - evaluated_heat;
                let gas_storage = next.mass_kg()
                    * (next.state().specific_internal_energy_j_kg() - gas.state().specific_internal_energy_j_kg())
                    + (next.mass_kg() - gas.mass_kg()) * gas.state().specific_internal_energy_j_kg();
                let gas_residual = (gas_storage + work_by_gas_j) - heat_j;
                if ![g, delta, power, evaluated_heat, residual, gas_storage, gas_residual]
                    .iter().all(|x| x.is_finite()) || g < 0.0
                    || (g != 0.0 && delta != 0.0 && power == 0.0)
                    || (power != 0.0 && evaluated_heat == 0.0)
                    || gas_residual.abs() > port.exchange.energy_tolerance_j
                {
                    return Err(EndpointError::Invalid("shared-body contact law or gas energy exceeds its declared allowance"));
                }
                Ok((Endpoint { gas: next, native, supplied_heat_j: heat_j,
                    port: Some(SharedFiniteGasAcousticThermalStep {
                        material: sample.source, body_temperature_k, gas_temperature_k,
                        conductance_w_k: g, solid_to_gas_heat_j: q,
                        source_residual_j: residual, gas_energy_residual_j: gas_residual,
                        exchange_trials: trials,
                    }),
                }, residual))
            }, cancelled,
        ).map_err(endpoint_feedback_error)?;
        let report = candidate
            .port
            .as_ref()
            .expect("bound shared-body gas retains its contact");
        let difference = report.body_temperature_k - report.gas_temperature_k;
        if (report.solid_to_gas_heat_j > 0.0 && difference < 0.0)
            || (report.solid_to_gas_heat_j < 0.0 && difference > 0.0)
        {
            return Err(EndpointError::Invalid(
                "shared-body heat opposes its actual endpoint temperature ordering",
            ));
        }
        poll(cancelled)?;
        Ok(candidate)
    }
}

impl FiniteGasTransportEvolution for SharedEvolution<'_, '_, '_> {
    type Source = SharedFiniteGasAcousticThermalSource;
    type Reservoir = SharedFiniteGasAcousticThermalReservoir;
    type Error = EndpointError;
    fn is_cancelled(error: &Self::Error) -> bool {
        matches!(error, EndpointError::Cancelled)
    }
    fn admit<C: FnMut() -> bool>(
        &self,
        id: u64,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        duration_s: f64,
        cancelled: &mut C,
    ) -> Result<(), Self::Error> {
        self.native
            .admit(id, gas, input, duration_s, cancelled)
            .map_err(Into::into)
    }
    fn node_candidate<C: FnMut() -> bool>(
        &self,
        id: u64,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        duration_s: f64,
        work_by_gas_j: f64,
        cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Source), Self::Error> {
        let endpoint = self.candidate(
            id,
            gas,
            input,
            duration_s,
            work_by_gas_j,
            cancelled,
            |heated, cancel| {
                self.native
                    .node_candidate(id, gas, heated, duration_s, work_by_gas_j, cancel)
            },
        )?;
        Ok((
            endpoint.gas,
            SharedFiniteGasAcousticThermalSource {
                native: endpoint.native,
                port: endpoint.port,
                supplied_heat_j: endpoint.supplied_heat_j,
            },
        ))
    }
    fn reservoir_candidate<C: FnMut() -> bool>(
        &self,
        id: u64,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        duration_s: f64,
        cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Reservoir), Self::Error> {
        let endpoint = self.candidate(
            id,
            gas,
            input,
            duration_s,
            0.0,
            cancelled,
            |heated, cancel| {
                self.native
                    .reservoir_candidate(id, gas, heated, duration_s, cancel)
            },
        )?;
        Ok((
            endpoint.gas,
            SharedFiniteGasAcousticThermalReservoir {
                native: endpoint.native,
                port: endpoint.port,
                supplied_heat_j: endpoint.supplied_heat_j,
            },
        ))
    }
    fn admit_node_endpoint(&self, id: u64, source: &Self::Source) -> Result<(), Self::Error> {
        self.native
            .admit_node_endpoint(id, &source.native)
            .map_err(Into::into)
    }
    fn node_supplied_heat_j(&self, _: &FiniteGasDgTransportInput, source: &Self::Source) -> f64 {
        source.supplied_heat_j
    }
    fn reservoir_supplied_heat_j(
        &self,
        _: &FiniteGasDgTransportInput,
        source: &Self::Reservoir,
    ) -> f64 {
        source.supplied_heat_j
    }
}

/// Preview several gas/acoustic volumes sharing ONE finite material heat store.
///
/// The original transport owner runs once on the accepted parcel map. At each
/// outer net-heat trial the body's canonical enthalpy inverse supplies a single
/// endpoint. Every contact then closes its own original heat/chemistry/acoustic
/// work residual against that endpoint, restarting the same transported gas.
/// Only after ALL contacts have been evaluated can the net-body residual accept;
/// opposite hot/cold transfers survive even with zero net heat and unchanged body
/// temperature. Initial conductivity claims remain pinned through all trials.
///
/// All prescribed heat and native external work are independently summed. The
/// one body storage cancels the sum of internal gas transfers in the enclosing
/// energy gate; no post-step renormalization, duplicate heat capacity or duplicate
/// formation energy is introduced. The whole-call endpoint budget includes
/// discarded body/heat/work trials. Each original owner keeps its other budgets.
///
/// This is a uniform, fixed-geometry, negligible-body-pV-work model with explicit
/// effective conduction layers. The caller asserts nonoverlapping patches and
/// the validity of lumping/density/geometry. It is not spatial solid conduction,
/// automatic duct-wall stiffness or aperture updates, expansion, resolved gas
/// motion or acoustic-loss heating. The original first-order transport/acoustic
/// split, source domains and local chemical entropy gate remain; convergence is
/// not time accuracy, uniqueness or a whole-system entropy theorem.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_shared_thermal_reacting_acoustic_transport_candidate_with_paths<
    C: FnMut() -> bool,
>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    body_state: &EquilibriumPhaseState,
    body: UniformEnthalpyStepInput<'_>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    reactions: &BTreeMap<u64, FiniteGasAcousticReaction<'_>>,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: SharedFiniteGasAcousticThermalSolve,
    cancelled: &mut C,
) -> Result<SharedThermalReactiveFiniteGasAcousticTransportCandidate, Error> {
    poll(cancelled)?;
    if !positive(controls.reactive.species_tolerance_mol)
        || reactions.len() > parcels.len()
        || reactions.keys().any(|id| !parcels.contains_key(id))
    {
        return Err(Error::Invalid(
            "positive species tolerance and reactions bound to original gas parcels required",
        ));
    }
    let pins = prepare_ports(
        parcels, body_state, body, ports, duration_s, controls, cancelled,
    )?;
    let evaluations = Cell::new(0usize);
    let initial_evolution = SharedEvolution {
        native: Evolution { reactions },
        body_state: *body_state,
        ports,
        pins: &pins,
        evaluations: &evaluations,
        maximum_evaluations: controls.maximum_endpoint_evaluations,
    };
    let snapshot = prepare_finite_gas_transport(
        parcels,
        nodes,
        inputs,
        paths,
        duration_s,
        controls.reactive.transport,
        &initial_evolution,
        cancelled,
    )?;
    let mut external_gas = Sum::default();
    for (&id, input) in inputs {
        poll(cancelled)?;
        let input = FiniteGasAcousticNode::transport_input(input, nodes.contains_key(&id))
            .map_err(Error::Invalid)?;
        external_gas.add(input.supplied_heat_j)?;
    }
    let h0 = body.initial.specific_enthalpy_j_kg();
    let mut candidate = feedback::solve_scalar(
        controls.body_exchange.heat_bracket_j,
        controls.body_exchange.source_tolerance_j,
        controls.body_exchange.maximum_iterations,
        |net_heat, body_trials, cancelled| {
            let body_heat = body.internal_heat_j - net_heat;
            let dh = body_heat / body.mass_kg;
            let h1 = h0 + dh;
            if ![body_heat, dh, h1].iter().all(|x| x.is_finite())
                || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                || (net_heat != 0.0 && body_heat == body.internal_heat_j)
                || (body.internal_heat_j != 0.0 && body_heat == -net_heat)
            {
                return Err(Error::Invalid(
                    "shared-body net heat disappears in its conserved enthalpy coordinate",
                ));
            }
            let state = body
                .curve
                .state_at_specific_enthalpy(h1)
                .map_err(EndpointError::from)?;
            let evolution = SharedEvolution {
                native: Evolution { reactions },
                body_state: state,
                ports,
                pins: &pins,
                evaluations: &evaluations,
                maximum_evaluations: controls.maximum_endpoint_evaluations,
            };
            let (next_parcels, next_nodes, transport) =
                snapshot.candidate(&evolution, cancelled)?.into_parts();
            let stored = body.mass_kg * (state.specific_enthalpy_j_kg() - h0);
            let mut step = SharedThermalReactiveFiniteGasAcousticTransportStep {
                transport,
                body_state: state,
                body_storage_change_j: stored,
                external_body_heat_j: body.internal_heat_j,
                external_gas_heat_j: external_gas.value,
                solid_to_gas_heat_j: 0.0,
                body_boundary_heat_j: -net_heat,
                body_energy_residual_j: 0.0,
                source_residual_j: 0.0,
                combined_energy_residual_j: 0.0,
                maximum_species_residual_mol: 0.0,
                body_trials,
                endpoint_evaluations: evaluations.get(),
            };
            let mut transferred = Sum::default();
            // Evaluate and retain each signed port even when their sum is zero.
            for &id in ports.keys() {
                poll(cancelled)?;
                let port = step.port(id).ok_or(Error::Invalid(
                    "accepted shared-body endpoint lost an original gas contact",
                ))?;
                transferred.add(port.solid_to_gas_heat_j)?;
            }
            let source_residual = net_heat - transferred.value;
            let body_residual = (stored + transferred.value) - body.internal_heat_j;
            if ![stored, source_residual, body_residual]
                .iter()
                .all(|x| x.is_finite())
            {
                return Err(Error::Invalid(
                    "shared-body storage or net heat residual is unrepresentable",
                ));
            }
            step.solid_to_gas_heat_j = transferred.value;
            step.source_residual_j = source_residual;
            step.body_energy_residual_j = body_residual;
            Ok((
                SharedThermalReactiveFiniteGasAcousticTransportCandidate {
                    parcels: next_parcels,
                    nodes: next_nodes,
                    step,
                },
                source_residual,
            ))
        },
        cancelled,
    )
    .map_err(outer_feedback_error)?;
    let step = &mut candidate.step;
    if step.body_energy_residual_j.abs() > controls.body_exchange.energy_tolerance_j {
        return Err(Error::Energy {
            residual_j: step.body_energy_residual_j,
            tolerance_j: controls.body_exchange.energy_tolerance_j,
        });
    }
    let mut total = Sum::default();
    for term in [
        step.transport.gas_storage_change_j,
        step.transport.acoustic_storage_change_j,
        step.transport.acoustic_dissipation_j,
        step.body_storage_change_j,
        -step.external_gas_heat_j,
        -step.external_body_heat_j,
        -step.transport.port_work_j,
    ] {
        total.add(term)?;
    }
    if total.value.abs() > controls.reactive.transport.energy_tolerance_j {
        return Err(Error::Energy {
            residual_j: total.value,
            tolerance_j: controls.reactive.transport.energy_tolerance_j,
        });
    }
    step.combined_energy_residual_j = total.value;
    step.maximum_species_residual_mol = species_residual(
        parcels,
        &candidate.parcels,
        reactions,
        |id| step.reaction(id),
        cancelled,
    )?;
    if step.maximum_species_residual_mol > controls.reactive.species_tolerance_mol {
        return Err(Error::SpeciesBalance {
            maximum_residual_mol: step.maximum_species_residual_mol,
            tolerance_mol: controls.reactive.species_tolerance_mol,
        });
    }
    step.endpoint_evaluations = evaluations.get();
    poll(cancelled)?;
    Ok(candidate)
}

/// Atomically publish the one body and every original gas/acoustic owner.
/// All returned contact heats have already been applied exactly once.
#[allow(clippy::too_many_arguments)]
pub fn advance_shared_thermal_reacting_finite_gas_acoustic_transport_with_paths<
    C: FnMut() -> bool,
>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, FiniteGasAcousticNode>,
    body_state: &mut EquilibriumPhaseState,
    body: UniformEnthalpyStepInput<'_>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    reactions: &BTreeMap<u64, FiniteGasAcousticReaction<'_>>,
    ports: &BTreeMap<u64, SharedFiniteGasAcousticThermalPort<'_>>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: SharedFiniteGasAcousticThermalSolve,
    cancelled: &mut C,
) -> Result<SharedThermalReactiveFiniteGasAcousticTransportStep, Error> {
    let candidate = finite_gas_shared_thermal_reacting_acoustic_transport_candidate_with_paths(
        parcels, nodes, body_state, body, inputs, reactions, ports, paths, duration_s, controls,
        cancelled,
    )?;
    poll(cancelled)?;
    let (next_parcels, next_nodes, next_body, step) = candidate.into_parts();
    *parcels = next_parcels;
    *nodes = next_nodes;
    *body_state = next_body;
    Ok(step)
}
