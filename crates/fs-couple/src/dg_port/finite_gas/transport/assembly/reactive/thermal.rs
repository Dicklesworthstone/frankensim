//! Finite material heat stores on the same reacting acoustic transport step.
//!
//! Each declared body retains its original enthalpy/phase chart. One signed
//! interface heat is removed from that body and supplied to its finite gas.
//! The existing scalar feedback owner closes source-backed layer conduction at
//! the same endpoint as chemistry and native acoustic parameter work. Rejected
//! heat/work probes never advance accepted thermal or acoustic history.

use super::{
    BoundaryGasReactionStep, Evolution, FiniteGasAcousticInput, FiniteGasAcousticNode,
    FiniteGasAcousticReaction, FiniteGasAcousticReactionReservoir, FiniteGasAcousticReactionSource,
    FiniteGasAcousticStep, FiniteGasDgTransportError, FiniteGasDgTransportInput,
    FiniteGasDgTransportStep, FiniteGasTransportEvolution, FiniteGasTransportNode, GasEnergyError,
    GasReactionError, GasTransportPaths, ReactiveFiniteGasAcousticTransportError,
    ReactiveFiniteGasAcousticTransportSolve, finite_gas_transport_candidate_with_evolution,
    species_residual,
};
use crate::cavity_phs::piston::thermal::{
    HeatExchangeSolve, feedback,
    layers::{LayerConductionError, MaterialConductionLayers},
};
use fs_material::phase::{
    EquilibriumPhaseState, PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_PARCELS;
use std::collections::BTreeMap;

/// One distinct finite body exchanging heat with its map-key gas parcel.
///
/// The body is a separate thermal store, not an implicit change to the duct or
/// aperture material. Its immutable chart and layer card must agree. The same
/// body may not be represented under several keys; shared spatial thermal
/// networks require their original multi-cell coupling owner.
pub struct FiniteGasAcousticThermalBinding<'a> {
    pub body: UniformEnthalpyStepInput<'a>,
    pub layers: MaterialConductionLayers<'a>,
    pub exchange: HeatExchangeSolve,
}

/// Exact accepted finite-body endpoint, interface source and reciprocal heat.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasAcousticThermalStep {
    /// Boundary heat into this body is exactly `-solid_to_gas_heat_j`.
    pub thermal: UniformEnthalpyStep<ResolvedMaterialStatePoint>,
    pub conductance_w_k: f64,
    pub gas_temperature_k: f64,
    pub solid_to_gas_heat_j: f64,
    pub external_body_heat_j: f64,
    pub body_storage_change_j: f64,
    pub body_energy_residual_j: f64,
    /// Q - dt G(T_body,end - T_gas,end) [J].
    pub source_residual_j: f64,
    /// Body + gas storage + gas-to-wave work - both external heat inputs [J].
    pub endpoint_energy_residual_j: f64,
    /// Complete heat proposals for this accepted acoustic work probe.
    pub exchange_trials: usize,
}

/// The actual native gas endpoint and, when bound, its matching finite body.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasAcousticThermalSource {
    pub native: FiniteGasAcousticReactionSource,
    pub thermal: Option<FiniteGasAcousticThermalStep>,
    /// Heat across the gas boundary: external gas heat plus body transfer.
    /// This excludes mechanical work and formation energy.
    pub supplied_heat_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasAcousticThermalReservoir {
    pub native: FiniteGasAcousticReactionReservoir,
    pub thermal: Option<FiniteGasAcousticThermalStep>,
    pub supplied_heat_j: f64,
}

/// Native gas/acoustic accounts and the enclosing finite-body energy balance.
#[derive(Debug)]
pub struct ThermalReactiveFiniteGasAcousticTransportStep {
    /// This boundary encloses gas and acoustics only, so its supplied heat
    /// includes the finite bodies' Q. The enclosing account below cancels Q.
    pub transport: FiniteGasDgTransportStep<
        FiniteGasAcousticStep<FiniteGasAcousticThermalSource>,
        FiniteGasAcousticThermalReservoir,
    >,
    pub maximum_species_residual_mol: f64,
    pub body_storage_change_j: f64,
    pub external_body_heat_j: f64,
    pub external_gas_heat_j: f64,
    pub solid_to_gas_heat_j: f64,
    pub combined_energy_residual_j: f64,
}

impl ThermalReactiveFiniteGasAcousticTransportStep {
    pub fn reaction(&self, id: u64) -> Option<&BoundaryGasReactionStep> {
        if let Some(step) = self.transport.acoustics.get(&id) {
            return match &FiniteGasAcousticNode::source::<FiniteGasAcousticThermalSource>(step)
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

    pub fn thermal(&self, id: u64) -> Option<&FiniteGasAcousticThermalStep> {
        if let Some(step) = self.transport.acoustics.get(&id) {
            return FiniteGasAcousticNode::source::<FiniteGasAcousticThermalSource>(step)
                .thermal
                .as_ref();
        }
        self.transport.reservoirs.get(&id)?.thermal.as_ref()
    }
}

/// Complete proposed maps, published together only after every physical gate.
pub struct ThermalReactiveFiniteGasAcousticTransportCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    nodes: BTreeMap<u64, FiniteGasAcousticNode>,
    body_states: BTreeMap<u64, EquilibriumPhaseState>,
    step: ThermalReactiveFiniteGasAcousticTransportStep,
}
impl ThermalReactiveFiniteGasAcousticTransportCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub fn nodes(&self) -> &BTreeMap<u64, FiniteGasAcousticNode> {
        &self.nodes
    }
    pub fn body_states(&self) -> &BTreeMap<u64, EquilibriumPhaseState> {
        &self.body_states
    }
    pub const fn step(&self) -> &ThermalReactiveFiniteGasAcousticTransportStep {
        &self.step
    }
    #[allow(clippy::type_complexity)]
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, FiniteGasAcousticNode>,
        BTreeMap<u64, EquilibriumPhaseState>,
        ThermalReactiveFiniteGasAcousticTransportStep,
    ) {
        (self.parcels, self.nodes, self.body_states, self.step)
    }
}

/// Original endpoint refusals remain distinct from transport/acoustic failures.
#[derive(Debug)]
pub enum FiniteGasAcousticThermalError {
    Invalid(&'static str),
    Gas(GasEnergyError),
    Reaction(GasReactionError),
    Phase(PhaseStateError),
    Layer(LayerConductionError),
    OutsideBracket {
        lower_residual_j: f64,
        upper_residual_j: f64,
    },
    NotConverged {
        iterations: usize,
        residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for FiniteGasAcousticThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite acoustic gas/body endpoint: {self:?}")
    }
}
impl std::error::Error for FiniteGasAcousticThermalError {}
impl From<GasEnergyError> for FiniteGasAcousticThermalError {
    fn from(error: GasEnergyError) -> Self {
        if error == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(error)
        }
    }
}
impl From<GasReactionError> for FiniteGasAcousticThermalError {
    fn from(error: GasReactionError) -> Self {
        if Evolution::is_cancelled(&error) {
            Self::Cancelled
        } else {
            Self::Reaction(error)
        }
    }
}
impl From<PhaseStateError> for FiniteGasAcousticThermalError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<LayerConductionError> for FiniteGasAcousticThermalError {
    fn from(error: LayerConductionError) -> Self {
        Self::Layer(error)
    }
}

#[derive(Debug)]
pub enum ThermalReactiveFiniteGasAcousticTransportError {
    Invalid(&'static str),
    Endpoint(FiniteGasAcousticThermalError),
    Native(FiniteGasDgTransportError<FiniteGasAcousticThermalError>),
    SpeciesAudit(ReactiveFiniteGasAcousticTransportError),
    SpeciesBalance {
        maximum_residual_mol: f64,
        tolerance_mol: f64,
    },
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ThermalReactiveFiniteGasAcousticTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal reacting acoustic assembly: {self:?}")
    }
}
impl std::error::Error for ThermalReactiveFiniteGasAcousticTransportError {}
impl From<FiniteGasAcousticThermalError> for ThermalReactiveFiniteGasAcousticTransportError {
    fn from(error: FiniteGasAcousticThermalError) -> Self {
        match error {
            FiniteGasAcousticThermalError::Cancelled => Self::Cancelled,
            other => Self::Endpoint(other),
        }
    }
}
impl From<GasEnergyError> for ThermalReactiveFiniteGasAcousticTransportError {
    fn from(error: GasEnergyError) -> Self {
        FiniteGasAcousticThermalError::from(error).into()
    }
}
impl From<FiniteGasDgTransportError<FiniteGasAcousticThermalError>>
    for ThermalReactiveFiniteGasAcousticTransportError
{
    fn from(error: FiniteGasDgTransportError<FiniteGasAcousticThermalError>) -> Self {
        match error {
            FiniteGasDgTransportError::Cancelled => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl From<ReactiveFiniteGasAcousticTransportError>
    for ThermalReactiveFiniteGasAcousticTransportError
{
    fn from(error: ReactiveFiniteGasAcousticTransportError) -> Self {
        match error {
            ReactiveFiniteGasAcousticTransportError::Cancelled => Self::Cancelled,
            other => Self::SpeciesAudit(other),
        }
    }
}

use FiniteGasAcousticThermalError as EndpointError;
use ThermalReactiveFiniteGasAcousticTransportError as Error;
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EndpointError> {
    if cancelled() {
        Err(EndpointError::Cancelled)
    } else {
        Ok(())
    }
}
fn feedback_error(error: feedback::ScalarError<EndpointError>) -> EndpointError {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => EndpointError::Invalid(message),
        feedback::ScalarError::Cancelled => EndpointError::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => EndpointError::OutsideBracket {
            lower_residual_j: lower,
            upper_residual_j: upper,
        },
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => EndpointError::NotConverged {
            iterations,
            residual_j: residual,
        },
    }
}

struct ThermalEvolution<'a, 'b, 'r> {
    native: Evolution<'a, 'r>,
    bindings: &'a BTreeMap<u64, FiniteGasAcousticThermalBinding<'b>>,
    pins: BTreeMap<u64, MaterialPropertySelection>,
}

// Validate every body and pin its actual conductivity before material transport.
fn prepare(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    states: &BTreeMap<u64, EquilibriumPhaseState>,
    bindings: &BTreeMap<u64, FiniteGasAcousticThermalBinding<'_>>,
    duration_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BTreeMap<u64, MaterialPropertySelection>, EndpointError> {
    if !states.keys().eq(bindings.keys()) {
        return Err(EndpointError::Invalid(
            "thermal bindings must exactly name the accepted body-state map",
        ));
    }
    let mut pins = BTreeMap::new();
    for (&id, binding) in bindings {
        poll(cancelled)?;
        let gas = parcels.get(&id).ok_or(EndpointError::Invalid(
            "finite body requires an existing gas parcel",
        ))?;
        let body = binding.body;
        let h0 = body.initial.specific_enthalpy_j_kg();
        if !positive(body.mass_kg)
            || !positive(body.volume_m3)
            || !positive(body.surface_area_m2)
            || !positive(duration_s)
            || body.duration_s != duration_s
            || !body.internal_heat_j.is_finite()
            || body.initial != states[&id]
            || body.curve.state_at_specific_enthalpy(h0)? != body.initial
            || body.curve.material_card_identity() != binding.layers.material.content_hash()
            || !binding.layers.area_m2.is_finite()
            || binding.layers.area_m2 < 0.0
            || binding.layers.area_m2 > body.surface_area_m2
        {
            return Err(EndpointError::Invalid(
                "one clock, canonical finite body, original material card and exposed contact area required",
            ));
        }
        feedback::admit::<core::convert::Infallible>(binding.exchange).map_err(|_| {
            EndpointError::Invalid(
                "valid thermal exchange bracket and work/energy controls required",
            )
        })?;
        let source = binding.layers.resolve(gas.state(), body.initial)?.source;
        let selected = source
            .property_by_key(binding.layers.conductivity)
            .ok_or(EndpointError::Invalid(
                "initial body lacks the exact conductivity source",
            ))?
            .answer()
            .receipt
            .selected;
        pins.insert(
            id,
            MaterialPropertySelection::PinnedByProperty(vec![(
                binding.layers.conductivity.name().to_owned(),
                selected,
            )]),
        );
    }
    poll(cancelled)?;
    Ok(pins)
}

struct Endpoint<S> {
    gas: ClosedGasParcel,
    native: S,
    thermal: Option<FiniteGasAcousticThermalStep>,
    supplied_heat_j: f64,
}

impl ThermalEvolution<'_, '_, '_> {
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
        let Some(binding) = self.bindings.get(&id) else {
            let (gas, native) = evolve(input, cancelled)?;
            return Ok(Endpoint {
                gas,
                native,
                thermal: None,
                supplied_heat_j: input.supplied_heat_j,
            });
        };
        let body = binding.body;
        let h0 = body.initial.specific_enthalpy_j_kg();
        let layers = MaterialConductionLayers {
            material: binding.layers.material,
            context: binding.layers.context,
            temperature_axis: binding.layers.temperature_axis,
            liquid_fraction_axis: binding.layers.liquid_fraction_axis,
            conductivity: binding.layers.conductivity,
            selection: self.pins[&id].clone(),
            area_m2: binding.layers.area_m2,
            gas_distance_m: binding.layers.gas_distance_m,
            body_distance_m: binding.layers.body_distance_m,
        };
        let candidate = feedback::solve_scalar(
            binding.exchange.heat_bracket_j,
            binding.exchange.source_tolerance_j,
            binding.exchange.maximum_iterations,
            |q, trials, cancelled| {
                let body_heat = body.internal_heat_j - q;
                let dh = body_heat / body.mass_kg;
                let h1 = h0 + dh;
                let gas_heat = input.supplied_heat_j + q;
                if ![body_heat, dh, h1, gas_heat].iter().all(|x| x.is_finite())
                    || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (q != 0.0
                        && (body_heat == body.internal_heat_j || gas_heat == input.supplied_heat_j))
                    || (body.internal_heat_j != 0.0 && body_heat == -q)
                    || (input.supplied_heat_j != 0.0 && gas_heat == q)
                {
                    return Err(EndpointError::Invalid(
                        "finite-body exchange or prescribed heat disappears in energy coordinates",
                    ));
                }
                let state = body.curve.state_at_specific_enthalpy(h1)?;
                let mut heated = input.clone();
                heated.supplied_heat_j = gas_heat;
                let (next, native) = evolve(&heated, cancelled)?;
                poll(cancelled)?;
                let sample = layers.resolve(next.state(), state)?;
                let gas_temperature_k = next.state().state().temperature;
                let delta = state.temperature_k() - gas_temperature_k;
                let g = sample.conductance_w_k;
                let power = g * delta;
                let heat = duration_s * power;
                let source_residual_j = q - heat;
                let stored = body.mass_kg * (state.specific_enthalpy_j_kg() - h0);
                let body_residual = (stored + q) - body.internal_heat_j;
                let gas_storage = next.mass_kg()
                    * (next.state().specific_internal_energy_j_kg()
                        - gas.state().specific_internal_energy_j_kg())
                    + (next.mass_kg() - gas.mass_kg())
                        * gas.state().specific_internal_energy_j_kg();
                let endpoint_residual = (gas_storage + stored) + work_by_gas_j
                    - input.supplied_heat_j
                    - body.internal_heat_j;
                if ![
                    g,
                    delta,
                    power,
                    heat,
                    source_residual_j,
                    stored,
                    body_residual,
                    endpoint_residual,
                ]
                .iter()
                .all(|x| x.is_finite())
                    || g < 0.0
                    || (g != 0.0 && delta != 0.0 && power == 0.0)
                    || (power != 0.0 && heat == 0.0)
                    || body_residual.abs() > binding.exchange.energy_tolerance_j
                    || endpoint_residual.abs() > binding.exchange.energy_tolerance_j
                {
                    return Err(EndpointError::Invalid(
                        "material heat law or gas/body energy exceeds its declared allowance",
                    ));
                }
                Ok((
                    Endpoint {
                        gas: next,
                        native,
                        supplied_heat_j: gas_heat,
                        thermal: Some(FiniteGasAcousticThermalStep {
                            thermal: UniformEnthalpyStep {
                                state,
                                external_heat_j: -q,
                                energy_residual_tolerance_j: binding.exchange.energy_tolerance_j,
                                report: sample.source,
                            },
                            conductance_w_k: g,
                            gas_temperature_k,
                            solid_to_gas_heat_j: q,
                            external_body_heat_j: body.internal_heat_j,
                            body_storage_change_j: stored,
                            body_energy_residual_j: body_residual,
                            source_residual_j,
                            endpoint_energy_residual_j: endpoint_residual,
                            exchange_trials: trials,
                        }),
                    },
                    source_residual_j,
                ))
            },
            cancelled,
        )
        .map_err(feedback_error)?;
        let report = candidate
            .thermal
            .as_ref()
            .expect("bound thermal endpoint retains its body");
        let difference = report.thermal.state.temperature_k() - report.gas_temperature_k;
        if (report.solid_to_gas_heat_j > 0.0 && difference < 0.0)
            || (report.solid_to_gas_heat_j < 0.0 && difference > 0.0)
        {
            return Err(EndpointError::Invalid(
                "accepted finite-body heat opposes the endpoint temperature ordering",
            ));
        }
        poll(cancelled)?;
        Ok(candidate)
    }
}

impl FiniteGasTransportEvolution for ThermalEvolution<'_, '_, '_> {
    type Source = FiniteGasAcousticThermalSource;
    type Reservoir = FiniteGasAcousticThermalReservoir;
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
            FiniteGasAcousticThermalSource {
                native: endpoint.native,
                thermal: endpoint.thermal,
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
            FiniteGasAcousticThermalReservoir {
                native: endpoint.native,
                thermal: endpoint.thermal,
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

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), Error> {
        let adjusted = term - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if ![term, adjusted, next, correction]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(Error::Invalid(
                "enclosing thermal assembly energy is unrepresentable",
            ));
        }
        self.value = next;
        self.correction = correction;
        Ok(())
    }
}

/// Propose transport, finite-body heat, chemistry and native acoustic work.
///
/// Sparse thermal and reaction maps are independent. All selected bodies start
/// from their exact accepted states, and every Q/work probe restarts the same
/// transported gas. Conduction uses the original material-layer resistance law
/// at the final gas/body temperatures, with its initial source claim pinned.
/// Chemistry receives actual heat Q_external + Q_body and separate mechanical
/// work. Formation energy is already included in gas internal energy.
///
/// Body storage uses the original negligible-pV-work enthalpy model at fixed
/// geometry. Uniform-body validity, material/density/geometry consistency and
/// distinct physical ownership of bodies remain caller obligations. This is
/// effective layer conduction, not resolved wall gradients or convection. Phase
/// states follow their chart; they do not update aperture mechanics, wall
/// stiffness, expansion, topology or acoustic-loss heat destinations.
///
/// Transport and native frozen-medium sampling retain their first-order split;
/// thermal/chemical endpoints share the same interval. Source-root convergence
/// establishes neither uniqueness nor time accuracy. The original chemical
/// entropy gate remains local, not a whole-assembly passivity theorem. All
/// physical maps stay unchanged on any refusal, exhausted work or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_thermal_reacting_acoustic_transport_candidate_with_paths<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    body_states: &BTreeMap<u64, EquilibriumPhaseState>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    reactions: &BTreeMap<u64, FiniteGasAcousticReaction<'_>>,
    thermal_bindings: &BTreeMap<u64, FiniteGasAcousticThermalBinding<'_>>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: ReactiveFiniteGasAcousticTransportSolve,
    cancelled: &mut C,
) -> Result<ThermalReactiveFiniteGasAcousticTransportCandidate, Error> {
    poll(cancelled)?;
    if parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || !positive(controls.species_tolerance_mol)
        || reactions.len() > parcels.len()
        || reactions.keys().any(|id| !parcels.contains_key(id))
        || thermal_bindings.len() > parcels.len()
    {
        return Err(Error::Invalid(
            "bounded gas graph, positive species tolerance and existing gas/body/reaction bindings required",
        ));
    }
    let pins = prepare(
        parcels,
        body_states,
        thermal_bindings,
        duration_s,
        cancelled,
    )?;
    let evolution = ThermalEvolution {
        native: Evolution { reactions },
        bindings: thermal_bindings,
        pins,
    };
    let candidate = finite_gas_transport_candidate_with_evolution(
        parcels,
        nodes,
        inputs,
        paths,
        duration_s,
        controls.transport,
        &evolution,
        cancelled,
    )?;
    let (next_parcels, next_nodes, transport) = candidate.into_parts();
    let mut step = ThermalReactiveFiniteGasAcousticTransportStep {
        transport,
        maximum_species_residual_mol: 0.0,
        body_storage_change_j: 0.0,
        external_body_heat_j: 0.0,
        external_gas_heat_j: 0.0,
        solid_to_gas_heat_j: 0.0,
        combined_energy_residual_j: 0.0,
    };
    step.maximum_species_residual_mol = species_residual(
        parcels,
        &next_parcels,
        reactions,
        |id| step.reaction(id),
        cancelled,
    )?;
    if step.maximum_species_residual_mol > controls.species_tolerance_mol {
        return Err(Error::SpeciesBalance {
            maximum_residual_mol: step.maximum_species_residual_mol,
            tolerance_mol: controls.species_tolerance_mol,
        });
    }
    let mut body_storage = Sum::default();
    let mut body_heat = Sum::default();
    let mut gas_heat = Sum::default();
    let mut exchange = Sum::default();
    let mut next_bodies = BTreeMap::new();
    for &id in body_states.keys() {
        poll(cancelled)?;
        let report = step
            .thermal(id)
            .ok_or(Error::Invalid("accepted endpoint lost a bound finite body"))?;
        body_storage.add(report.body_storage_change_j)?;
        body_heat.add(report.external_body_heat_j)?;
        exchange.add(report.solid_to_gas_heat_j)?;
        next_bodies.insert(id, report.thermal.state);
    }
    for (&id, input) in inputs {
        poll(cancelled)?;
        let original = FiniteGasAcousticNode::transport_input(input, nodes.contains_key(&id))
            .map_err(Error::Invalid)?;
        gas_heat.add(original.supplied_heat_j)?;
    }
    let mut total = Sum::default();
    for term in [
        step.transport.gas_storage_change_j,
        step.transport.acoustic_storage_change_j,
        step.transport.acoustic_dissipation_j,
        body_storage.value,
        -gas_heat.value,
        -body_heat.value,
        -step.transport.port_work_j,
    ] {
        total.add(term)?;
    }
    if total.value.abs() > controls.transport.energy_tolerance_j {
        return Err(Error::Energy {
            residual_j: total.value,
            tolerance_j: controls.transport.energy_tolerance_j,
        });
    }
    step.body_storage_change_j = body_storage.value;
    step.external_body_heat_j = body_heat.value;
    step.external_gas_heat_j = gas_heat.value;
    step.solid_to_gas_heat_j = exchange.value;
    step.combined_energy_residual_j = total.value;
    poll(cancelled)?;
    Ok(ThermalReactiveFiniteGasAcousticTransportCandidate {
        parcels: next_parcels,
        nodes: next_nodes,
        body_states: next_bodies,
        step,
    })
}

/// Publish all three original state maps after the same complete proposal and
/// a final cancellation boundary. The returned Q has already been applied.
#[allow(clippy::too_many_arguments)]
pub fn advance_thermal_reacting_finite_gas_acoustic_transport_with_paths<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, FiniteGasAcousticNode>,
    body_states: &mut BTreeMap<u64, EquilibriumPhaseState>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    reactions: &BTreeMap<u64, FiniteGasAcousticReaction<'_>>,
    thermal_bindings: &BTreeMap<u64, FiniteGasAcousticThermalBinding<'_>>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: ReactiveFiniteGasAcousticTransportSolve,
    cancelled: &mut C,
) -> Result<ThermalReactiveFiniteGasAcousticTransportStep, Error> {
    let candidate = finite_gas_thermal_reacting_acoustic_transport_candidate_with_paths(
        parcels,
        nodes,
        body_states,
        inputs,
        reactions,
        thermal_bindings,
        paths,
        duration_s,
        controls,
        cancelled,
    )?;
    poll(cancelled)?;
    let (next_parcels, next_nodes, next_bodies, step) = candidate.into_parts();
    *parcels = next_parcels;
    *nodes = next_nodes;
    *body_states = next_bodies;
    Ok(step)
}
