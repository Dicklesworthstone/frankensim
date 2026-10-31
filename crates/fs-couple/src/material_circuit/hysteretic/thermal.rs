//! Finite core pinning heat and winding I^2 R on the native enthalpy timeline.
//!
//! The immutable core fit keeps its parameters throughout an admitted
//! temperature interval. A source-bound core requires the SAME original card,
//! complete query and explicit temperature coverage for every fitted coefficient;
//! authored cores retain their explicit caller-admitted model interval.
//! Caloric storage and winding resistivity/conductivity come from one prepared
//! MaterialNetwork with pinned sources. Its complete heat assembly can add
//! source-bound thermal contacts, integrated layers and enclosure radiation.
//! Every thermal trial restarts the same accepted circuit and magnetic memory.
//! Core heat enters its whole fixed solid volume; graph resistor heat enters
//! its original receiving cells. No core-loss map, eddy fit, numerical defect
//! or duplicate thermal source is added. Endpoint transport plus thermal BE
//! gives first-order coupled accuracy in general, requiring time refinement.

use super::{
    HystereticCircuit, HystereticCircuitCandidate, HystereticCircuitError, HystereticCircuitSolve,
    HystereticCircuitStep,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, DriverAccount, ThermalDriver, cell_driver_candidate,
    contacts::heat::{
        HeatAssembly, HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance,
    },
    material::{
        MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
        PreparedMaterialNetwork,
    },
};
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::{
    conductor::ElectricalPropertyBasis, field_storage::hysteresis::UniformHystereticCore,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

pub(super) mod sealed {
    pub trait History {}
    pub trait Step {}
}

/// Sealed native magnetic histories using the SAME finite core/coil heat owner.
/// The source fit, retained memory and original native electrical account stay
/// with the single-core or shared-core implementation; no thermal trial may
/// commit or advance that accepted history more than once.
pub trait HystereticThermalHistory: sealed::History + Clone {
    type Step: HystereticThermalStep;
    fn uniform_core(&self) -> &UniformHystereticCore;
    fn native_circuit(&self) -> &crate::material_circuit::MaterialCircuit;
    fn time_s(&self) -> f64;
    fn with_core_source_query(&self, query: &QueryPoint) -> Result<Self, HystereticCircuitError>;
    fn step_for_heat(
        &self,
        target_s: f64,
        resistors: Vec<crate::material_circuit::MaterialResistor>,
        external_inputs: &[f64],
        controls: HystereticCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Self::Step), HystereticCircuitError>;
}

/// Original native electrical, magnetic storage and ONE core-loss receipt.
/// Sealed implementations prevent arbitrary loss maps entering this account.
pub trait HystereticThermalStep: sealed::Step {
    fn electrical(&self) -> &crate::material_circuit::MaterialCircuitStep;
    fn core_heat_j(&self) -> f64;
    fn heat_j(&self) -> f64;
    fn numerical_dissipation_j(&self) -> f64;
    fn stored_change_j(&self) -> f64;
    fn supplied_work_j(&self) -> f64;
    fn circuit_trials(&self) -> usize;
}

impl sealed::History for HystereticCircuit {}
impl HystereticThermalHistory for HystereticCircuit {
    type Step = HystereticCircuitStep;
    fn uniform_core(&self) -> &UniformHystereticCore {
        self.core()
    }
    fn native_circuit(&self) -> &crate::material_circuit::MaterialCircuit {
        HystereticCircuit::native_circuit(self)
    }
    fn time_s(&self) -> f64 {
        HystereticCircuit::time_s(self)
    }
    fn with_core_source_query(&self, query: &QueryPoint) -> Result<Self, HystereticCircuitError> {
        Ok(Self {
            core: self.core.at_source_query(query)?,
            ..self.clone()
        })
    }
    fn step_for_heat(
        &self,
        target_s: f64,
        resistors: Vec<crate::material_circuit::MaterialResistor>,
        external_inputs: &[f64],
        controls: HystereticCircuitSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Self::Step), HystereticCircuitError> {
        self.step_candidate(target_s, resistors, external_inputs, controls, cancelled)
            .map(HystereticCircuitCandidate::into_parts)
    }
}
impl sealed::Step for HystereticCircuitStep {}
impl HystereticThermalStep for HystereticCircuitStep {
    fn electrical(&self) -> &crate::material_circuit::MaterialCircuitStep {
        &self.circuit
    }
    fn core_heat_j(&self) -> f64 {
        self.core.core_heat_j
    }
    fn heat_j(&self) -> f64 {
        self.heat_j
    }
    fn numerical_dissipation_j(&self) -> f64 {
        self.numerical_dissipation_j
    }
    fn stored_change_j(&self) -> f64 {
        self.stored_change_j
    }
    fn supplied_work_j(&self) -> f64 {
        self.supplied_work_j
    }
    fn circuit_trials(&self) -> usize {
        self.circuit_trials
    }
}

struct MagneticDriver<H: HystereticThermalHistory> {
    model: H,
    step: H::Step,
}
impl<H: HystereticThermalHistory> ThermalDriver for MagneticDriver<H> {
    type Model = H;
    type Step = H::Step;
    fn account(&self) -> DriverAccount<'_> {
        let step = &self.step;
        DriverAccount {
            resistors: &step.electrical().resistors,
            contacts: &step.electrical().contacts,
            heat_j: step.heat_j(),
            numerical_dissipation_j: step.numerical_dissipation_j(),
            stored_change_j: step.stored_change_j(),
            supplied_j: step.supplied_work_j(),
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.model, self.step)
    }
}

#[derive(Debug)]
pub enum HystereticThermalError<E = HeatNetworkError> {
    Material(MaterialNetworkError),
    Magnetic(HystereticCircuitError),
    Heat(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for HystereticThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hysteretic circuit/enthalpy: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for HystereticThermalError<E> {}
type Error<E = HeatNetworkError> = CircuitEnthalpyError<HystereticThermalError<E>>;
fn invalid<E>(s: &'static str) -> Error<E> {
    CircuitEnthalpyError::Invalid(s)
}
fn material<E>(e: MaterialNetworkError) -> Error<E> {
    match e {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(HystereticThermalError::Material(other)),
    }
}
fn magnetic<E>(e: HystereticCircuitError) -> Error<E> {
    match e {
        HystereticCircuitError::Cancelled => CircuitEnthalpyError::Cancelled,
        HystereticCircuitError::TrialBudget => CircuitEnthalpyError::CircuitBudget,
        other => CircuitEnthalpyError::Source(HystereticThermalError::Magnetic(other)),
    }
}
fn heat<E: From<HeatNetworkError>>(e: HeatNetworkError) -> Error<E> {
    match e {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(HystereticThermalError::Heat(other.into())),
    }
}
impl ThermalDriver for HystereticCircuitCandidate {
    type Model = HystereticCircuit;
    type Step = HystereticCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        let s = self.step();
        DriverAccount {
            resistors: &s.circuit.resistors,
            contacts: &s.circuit.contacts,
            heat_j: s.heat_j,
            numerical_dissipation_j: s.numerical_dissipation_j,
            stored_change_j: s.stored_change_j,
            supplied_j: s.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        HystereticCircuitCandidate::into_parts(self)
    }
}

/// Accepted circuit/core history retaining ONE prepared source selection.
pub struct ThermalHystereticCircuit<'n, 'a, H: HystereticThermalHistory = HystereticCircuit> {
    circuit: H,
    network: &'n MaterialNetwork<'a>,
    prepared: Arc<PreparedMaterialNetwork<'n, 'a>>,
    destinations: Arc<BTreeMap<usize, CellId>>,
    core_cell: CellId,
    core_density_kg_m3: f64,
    core_temperature_interval_k: [f64; 2],
    core_binding: Option<&'n MaterialCell<'a>>,
}
impl<'n, 'a, H: HystereticThermalHistory> ThermalHystereticCircuit<'n, 'a, H> {
    /// Bind the core to a distinct whole fixed-volume solid caloric cell.
    /// Its mass/density must reproduce UniformHystereticCore::volume_m3.
    /// Graph resistor geometry and heat routes stay with MaterialNetwork.
    /// Source-bound magnetic cores must retain this exact caloric card and
    /// complete initial query. Every fitted coefficient must explicitly support
    /// the requested temperature interval; each trial rechecks those SAME
    /// sources and publishes the accepted source query with magnetic memory.
    /// Source-constant parameters introduce no reversible thermal field work.
    /// An authored core keeps its explicit caller interval and authority.
    #[allow(clippy::too_many_arguments)]
    pub fn bind(
        circuit: H,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        core_cell: CellId,
        core_temperature_interval_k: [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::bind_with_link_conductances(
            circuit,
            thermal,
            network,
            electrical_bases,
            core_cell,
            core_temperature_interval_k,
            &BTreeMap::new(),
            cancelled,
        )
    }

    /// Bind the same core, winding and finite thermal inventory while another
    /// physical owner supplies the remaining native link conductances.
    ///
    /// The original material network and this initial map must partition every
    /// native thermal link exactly. Values are used only for initial admission;
    /// they are not cached as future cooling. Continue through
    /// `step_candidate_with_heat_assembly` with a provider that re-resolves all
    /// finite links at every actual trial temperature. Ordinary `step_candidate`
    /// still refuses unowned links rather than silently freezing their values.
    ///
    /// Core source identity, whole solid volume, temperature coverage, accepted
    /// magnetic memory and original winding source checks are identical to
    /// `bind`. The fitted magnetic coefficients remain constant over their
    /// admitted interval; cooling can change R(T), current and the resulting
    /// hysteresis trajectory, but supplies no unmeasured thermal magnetic law.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_with_link_conductances(
        circuit: H,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        core_cell: CellId,
        core_temperature_interval_k: [f64; 2],
        initial_link_conductances_w_k: &BTreeMap<u64, f64>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        let [low, high] = core_temperature_interval_k;
        if circuit.time_s().to_bits() != thermal.time_s().to_bits()
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || network
                .resistors
                .iter()
                .any(|r| r.destination.cell == core_cell)
        {
            return Err(invalid(
                "matching clocks, positive core temperature interval and distinct core/winding cells required",
            ));
        }
        let prepared = network
            .prepare_with_link_conductances(
                thermal,
                electrical_bases,
                initial_link_conductances_w_k,
                &mut cancelled,
            )
            .map_err(material)?;
        let binding = network
            .cells
            .iter()
            .find(|b| b.cell == core_cell)
            .ok_or_else(|| invalid("core lacks its actual caloric/material binding"))?;
        let first = binding
            .curve
            .knots()
            .first()
            .ok_or_else(|| invalid("empty core caloric chart"))?;
        let density = first.bulk_density_kg_m3;
        if !density.is_finite()
            || density <= 0.0
            || binding.curve.knots().iter().any(|k| {
                k.liquid_mass_fraction != 0.0
                    || (k.bulk_density_kg_m3 - density).abs() > 1e-10 * density
            })
        {
            return Err(invalid(
                "fixed hysteretic core requires a wholly solid constant-density caloric chart",
            ));
        }
        let mut initial = prepared
            .resolve_with_link_conductances(thermal, initial_link_conductances_w_k, &mut cancelled)
            .map_err(material)?;
        initial.resistors.sort_by_key(|r| r.branch);
        if initial.resistors != circuit.native_circuit().resistors() {
            return Err(invalid(
                "accepted circuit resistor sources must match the actual thermal endpoint",
            ));
        }
        let core_binding = if let Some(source) = circuit.uniform_core().source() {
            let cell = thermal
                .cells()
                .get(&core_cell)
                .ok_or_else(|| invalid("source-bound core cell is inactive"))?;
            let phase = binding
                .curve
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(MaterialNetworkError::Phase)
                .map_err(material)?;
            let query = binding
                .point
                .at(phase)
                .map_err(MaterialNetworkError::Material)
                .map_err(material)?;
            if source.material().card_identity() != binding.point.material.content_hash()
                || source.material().query() != query
            {
                return Err(invalid(
                    "magnetic fit and caloric core require the same original card and complete accepted query",
                ));
            }
            source
                .admit_temperature_interval(
                    binding.point.temperature_axis,
                    core_temperature_interval_k,
                )
                .map_err(HystereticCircuitError::from)
                .map_err(magnetic)?;
            Some(binding)
        } else {
            None
        };
        let destinations = network
            .resistors
            .iter()
            .map(|r| (r.destination.branch, r.destination.cell))
            .collect();
        let result = Self {
            circuit,
            network,
            prepared: Arc::new(prepared),
            destinations: Arc::new(destinations),
            core_cell,
            core_density_kg_m3: density,
            core_temperature_interval_k,
            core_binding,
        };
        let _ = result.core_query(thermal)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(result)
    }

    pub const fn circuit(&self) -> &H {
        &self.circuit
    }
    pub const fn core_cell(&self) -> CellId {
        self.core_cell
    }
    pub const fn core_temperature_interval_k(&self) -> [f64; 2] {
        self.core_temperature_interval_k
    }

    fn core_query<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<Option<QueryPoint>, Error<E>> {
        let cell = thermal
            .cells()
            .get(&self.core_cell)
            .ok_or_else(|| invalid("core heat cell is inactive"))?;
        let volume = self.circuit.uniform_core().volume_m3();
        let mass_volume = cell.mass_kg / self.core_density_kg_m3;
        let temperature = thermal.temperature_k(self.core_cell)?;
        if !mass_volume.is_finite()
            || (cell.geometry.volume_m3 / volume - 1.0).abs() > 1e-10
            || (mass_volume / volume - 1.0).abs() > 1e-10
            || temperature < self.core_temperature_interval_k[0]
            || temperature > self.core_temperature_interval_k[1]
        {
            return Err(invalid(
                "core must retain its actual fixed solid volume and admitted temperature interval",
            ));
        }
        self.core_binding
            .map(|binding| {
                let phase = binding
                    .curve
                    .state_at_specific_enthalpy(cell.specific_h_j_kg)
                    .map_err(MaterialNetworkError::Phase)
                    .map_err(material)?;
                binding
                    .point
                    .at(phase)
                    .map_err(MaterialNetworkError::Material)
                    .map_err(material)
            })
            .transpose()
    }

    /// Joint magnetic, native circuit and finite caloric proposal. Neither input
    /// is mutated. All original resistor heat plus core pinning heat are already
    /// deposited in the returned thermal state, once each.
    ///
    /// controls.circuit supplies the native electrical work/power tolerances;
    /// magnetic_controls supplies interface Newton work and magnetic acceptance.
    /// Both work budgets count every native evaluation, including discarded
    /// stencils and final reversal probes. Actual magnetic reversals and thermal
    /// events require interval cuts. Endpoint core/coil temperatures and all
    /// numerical stencils must stay within their explicit source/model domains.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        magnetic_controls: HystereticCircuitSolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticEnthalpyCandidate<'t, 'n, 'a, MaterialNetworkSample, H>, Error> {
        let cancel = RefCell::new(cancelled);
        self.step_resolved_candidate(
            thermal,
            target_s,
            external_inputs,
            controls,
            magnetic_controls,
            initial_energy_increments_j,
            |state| {
                self.prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            || (*cancel.borrow_mut())(),
        )
    }

    /// Couple the original magnetic memory and winding to its complete heat assembly.
    ///
    /// Thermal contacts, source-integrated layers and closed-enclosure radiation
    /// see the SAME trial core/winding temperatures as source-bound resistance
    /// and magnetic pinning heat. The heat assembly must wrap the original
    /// bulk MaterialNetwork supplied to `bind`; no caloric binding, conductor,
    /// branch or heat receiver is silently replaced. Retain the prepared heat
    /// network across steps to preserve every added contact and emissivity pin.
    ///
    /// Electrical contact branches need a magnetic circuit contact-work owner
    /// and are refused. Thermal interfaces and emissivity must be independent
    /// of the evolving electric/magnetic fields. Closed radiation stays an
    /// internal transfer: the native
    /// source report retains its roundoff-net power, with an additional total
    /// energy gate. True external heat is `non_electrical_heat_j` minus
    /// `dt * report.source.radiation_net_power_w()`.
    /// The original fixed magnetic coefficients, monotone-linkage intervals,
    /// fixed solid geometry and first-order thermal limits remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_heat_network<'t>(
        &self,
        cx: &Cx<'_>,
        prepared: &PreparedHeatNetwork<'_, '_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        magnetic_controls: HystereticCircuitSolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<HystereticEnthalpyCandidate<'t, 'n, 'a, HeatNetworkSample, H>, Error> {
        self.step_candidate_with_heat_assembly(
            cx,
            prepared,
            thermal,
            target_s,
            external_inputs,
            controls,
            magnetic_controls,
            initial_energy_increments_j,
        )
    }

    /// Advance original hysteresis memory, winding current and finite heat with
    /// any complete material heat assembly, including finite source-resolved
    /// coolant and prescribed material reservoirs.
    ///
    /// The provider must retain the exact bulk network supplied at binding.
    /// Every trial resolves the SAME proposed thermal state and starts the
    /// magnetic/circuit solve from the same accepted memory. Native winding
    /// Joule loss and core pinning loss are deposited once in their original
    /// distinct caloric receivers; finite coolant exchange is internal. The
    /// complete accepted provider sample remains in `report().source`.
    ///
    /// The original physical energy gate excludes closed-radiation roundoff
    /// from external supply. Typed provider failures remain under
    /// `HystereticThermalError::Heat`; native cancellation and work-budget
    /// failures keep their classifications. No partially advanced memory is
    /// published on refusal.
    ///
    /// Use `bind_with_link_conductances` when the thermal inventory has
    /// additional finite links. Existing constant magnetic-fit coefficients,
    /// admitted core temperature interval, monotone-linkage cuts, fixed solid
    /// geometry and first-order thermal coupling remain authoritative.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_heat_assembly<'t, A: HeatAssembly + ?Sized>(
        &self,
        cx: &Cx<'_>,
        assembly: &A,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        magnetic_controls: HystereticCircuitSolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<HystereticEnthalpyCandidate<'t, 'n, 'a, A::Sample, H>, Error<A::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let prepared = assembly.prepared_heat_network();
        let network = prepared
            .field_material_network()
            .map_err(heat::<A::Error>)?;
        if !std::ptr::eq(network, self.network) {
            return Err(invalid(
                "hysteretic heat assembly must retain the original bulk material network",
            ));
        }
        prepared
            .admit_field_independent_boundaries(cx, &[])
            .map_err(heat::<A::Error>)?;
        let candidate = self.step_resolved_candidate(
            thermal,
            target_s,
            external_inputs,
            controls,
            magnetic_controls,
            initial_energy_increments_j,
            |state| {
                let coefficients = assembly
                    .resolve_heat(cx, state)
                    .map_err(|error| error.map_source(HystereticThermalError::Heat))?;
                if !coefficients.contacts.is_empty() {
                    return Err(invalid(
                        "hysteretic circuit cannot discard electrical contact work or heat",
                    ));
                }
                Ok(coefficients.bulk)
            },
            || cx.checkpoint().is_err(),
        )?;
        closed_balance(
            target_s - thermal.time_s(),
            candidate.report.total_energy_residual_j,
            candidate.report.source.as_ref(),
            controls.total_energy_tolerance_j,
        )
        .map_err(heat::<A::Error>)?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }

    // One original magnetic driver, heat partition and enthalpy residual serve
    // both source assemblies. Every trial starts from the same accepted memory.
    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        magnetic_controls: HystereticCircuitSolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<CircuitThermalCoefficients<S>, Error<E>>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticEnthalpyCandidate<'t, 'n, 'a, S, H>, Error<E>> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() {
                Err(CircuitEnthalpyError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        if self.circuit.time_s().to_bits() != thermal.time_s().to_bits() {
            return Err(invalid(
                "accepted native circuit and thermal clocks must agree",
            ));
        }
        if let Some(query) = self.core_query(thermal)? {
            if self
                .circuit
                .uniform_core()
                .source()
                .map(|source| source.material().query())
                != Some(query)
            {
                return Err(invalid(
                    "accepted magnetic source query must match the actual thermal state",
                ));
            }
        }
        let mut initial = resolve(thermal)?;
        initial.resistors.sort_by_key(|r| r.branch);
        if initial.resistors != self.circuit.native_circuit().resistors() {
            return Err(invalid(
                "accepted native resistor sources do not match the thermal state",
            ));
        }
        let (thermal, circuit, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, remaining| {
                check()?;
                let query = self.core_query(state)?;
                let coefficients = resolve(state)?;
                let mut inner = magnetic_controls;
                inner.circuit = controls.circuit;
                inner.maximum_circuit_trials = inner.maximum_circuit_trials.min(remaining);
                if inner.maximum_circuit_trials < 4 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                // Recheck and retain the actual trial source context while all
                // magnetic memory remains at the one accepted beginning state.
                let beginning = match query {
                    Some(query) => self
                        .circuit
                        .with_core_source_query(&query)
                        .map_err(magnetic)?,
                    None => self.circuit.clone(),
                };
                let (model, step) = beginning
                    .step_for_heat(
                        target_s,
                        coefficients.resistors,
                        external_inputs,
                        inner,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(magnetic)?;
                let circuit_trials = step.circuit_trials();
                Ok(DrivenCoefficients {
                    circuit: MagneticDriver { model, step },
                    circuit_trials,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                })
            },
            |candidate, dt, checkpoint| {
                let mut powers = BTreeMap::new();
                for loss in &candidate.step.electrical().resistors {
                    checkpoint()?;
                    let cell = *self
                        .destinations
                        .get(&loss.material.branch)
                        .ok_or_else(|| {
                            invalid("native resistor heat lacks its original receiver")
                        })?;
                    let power = loss.heat_j / dt;
                    if !power.is_finite() || (loss.heat_j != 0.0 && power == 0.0) {
                        return Err(invalid(
                            "resistor heat cannot be represented as a thermal rate",
                        ));
                    }
                    powers.insert(cell, power);
                }
                checkpoint()?;
                let core_heat = candidate.step.core_heat_j();
                let core_power = core_heat / dt;
                if !core_power.is_finite() || (core_heat != 0.0 && core_power == 0.0) {
                    return Err(invalid(
                        "core pinning heat cannot be represented as a thermal rate",
                    ));
                }
                powers.insert(self.core_cell, core_power);
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check()?;
        Ok(HystereticEnthalpyCandidate {
            thermal,
            circuit: Self {
                circuit,
                network: self.network,
                prepared: Arc::clone(&self.prepared),
                destinations: Arc::clone(&self.destinations),
                core_cell: self.core_cell,
                core_density_kg_m3: self.core_density_kg_m3,
                core_temperature_interval_k: self.core_temperature_interval_k,
                core_binding: self.core_binding,
            },
            report,
        })
    }
}

/// Both publishable histories. The report's legacy deposited_joule_j field
/// counts TOTAL physical deposition here, including the separately reported
/// core pinning heat; numerical dissipation is never part of that deposition.
pub struct HystereticEnthalpyCandidate<
    't,
    'n,
    'a,
    S = MaterialNetworkSample,
    H: HystereticThermalHistory = HystereticCircuit,
> {
    thermal: DepositionSimulation<'t>,
    circuit: ThermalHystereticCircuit<'n, 'a, H>,
    report: CircuitEnthalpyReport<S, H::Step>,
}
impl<'t, 'n, 'a, S, H: HystereticThermalHistory> HystereticEnthalpyCandidate<'t, 'n, 'a, S, H> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn circuit(&self) -> &ThermalHystereticCircuit<'n, 'a, H> {
        &self.circuit
    }
    pub const fn report(&self) -> &CircuitEnthalpyReport<S, H::Step> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermalHystereticCircuit<'n, 'a, H>,
        CircuitEnthalpyReport<S, H::Step>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}
