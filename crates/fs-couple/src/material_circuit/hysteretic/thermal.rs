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
    contacts::heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance},
    material::{
        MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
        PreparedMaterialNetwork,
    },
};
use fs_exec::Cx;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

#[derive(Debug)]
pub enum HystereticThermalError {
    Material(MaterialNetworkError),
    Magnetic(HystereticCircuitError),
    Heat(HeatNetworkError),
}
impl core::fmt::Display for HystereticThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hysteretic circuit/enthalpy: {self:?}")
    }
}
impl std::error::Error for HystereticThermalError {}
type Error = CircuitEnthalpyError<HystereticThermalError>;
fn invalid(s: &'static str) -> Error {
    CircuitEnthalpyError::Invalid(s)
}
fn material(e: MaterialNetworkError) -> Error {
    match e {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(HystereticThermalError::Material(other)),
    }
}
fn magnetic(e: HystereticCircuitError) -> Error {
    match e {
        HystereticCircuitError::Cancelled => CircuitEnthalpyError::Cancelled,
        HystereticCircuitError::TrialBudget => CircuitEnthalpyError::CircuitBudget,
        other => CircuitEnthalpyError::Source(HystereticThermalError::Magnetic(other)),
    }
}
fn heat(e: HeatNetworkError) -> Error {
    match e {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(HystereticThermalError::Heat(other)),
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
pub struct ThermalHystereticCircuit<'n, 'a> {
    circuit: HystereticCircuit,
    network: &'n MaterialNetwork<'a>,
    prepared: Arc<PreparedMaterialNetwork<'n, 'a>>,
    destinations: Arc<BTreeMap<usize, CellId>>,
    core_cell: CellId,
    core_density_kg_m3: f64,
    core_temperature_interval_k: [f64; 2],
    core_binding: Option<&'n MaterialCell<'a>>,
}
impl<'n, 'a> ThermalHystereticCircuit<'n, 'a> {
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
        circuit: HystereticCircuit,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        core_cell: CellId,
        core_temperature_interval_k: [f64; 2],
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
            .prepare_with_bases(thermal, electrical_bases, &mut cancelled)
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
            .resolve(thermal, &mut cancelled)
            .map_err(material)?;
        initial.resistors.sort_by_key(|r| r.branch);
        if initial.resistors != circuit.native_circuit().resistors() {
            return Err(invalid(
                "accepted circuit resistor sources must match the actual thermal endpoint",
            ));
        }
        let core_binding = if let Some(source) = circuit.core().source() {
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
        let _ = result.admit_core(thermal)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(result)
    }

    pub const fn circuit(&self) -> &HystereticCircuit {
        &self.circuit
    }
    pub const fn core_cell(&self) -> CellId {
        self.core_cell
    }
    pub const fn core_temperature_interval_k(&self) -> [f64; 2] {
        self.core_temperature_interval_k
    }

    fn admit_core(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<Option<fs_material::field_storage::hysteresis::UniformHystereticCore>, Error> {
        let cell = thermal
            .cells()
            .get(&self.core_cell)
            .ok_or_else(|| invalid("core heat cell is inactive"))?;
        let volume = self.circuit.core().volume_m3();
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
                let query = binding
                    .point
                    .at(phase)
                    .map_err(MaterialNetworkError::Material)
                    .map_err(material)?;
                self.circuit
                    .core()
                    .at_source_query(&query)
                    .map_err(HystereticCircuitError::from)
                    .map_err(magnetic)
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
    ) -> Result<HystereticEnthalpyCandidate<'t, 'n, 'a>, Error> {
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
    ) -> Result<HystereticEnthalpyCandidate<'t, 'n, 'a, HeatNetworkSample>, Error> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let network = prepared.field_material_network().map_err(heat)?;
        if !std::ptr::eq(network, self.network) {
            return Err(invalid(
                "hysteretic heat assembly must retain the original bulk material network",
            ));
        }
        prepared
            .admit_field_independent_boundaries(cx, &[])
            .map_err(heat)?;
        let candidate = self.step_resolved_candidate(
            thermal,
            target_s,
            external_inputs,
            controls,
            magnetic_controls,
            initial_energy_increments_j,
            |state| {
                let coefficients = prepared.resolve(cx, state).map_err(heat)?;
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
            &candidate.report.source,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat)?;
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }

    // One original magnetic driver, heat partition and enthalpy residual serve
    // both source assemblies. Every trial starts from the same accepted memory.
    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, S>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        magnetic_controls: HystereticCircuitSolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<CircuitThermalCoefficients<S>, Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticEnthalpyCandidate<'t, 'n, 'a, S>, Error> {
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
        if self
            .admit_core(thermal)?
            .is_some_and(|core| core != *self.circuit.core())
        {
            return Err(invalid(
                "accepted magnetic source query must match the actual thermal state",
            ));
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
                let core = self.admit_core(state)?;
                let coefficients = resolve(state)?;
                let mut inner = magnetic_controls;
                inner.circuit = controls.circuit;
                inner.maximum_circuit_trials = inner.maximum_circuit_trials.min(remaining);
                if inner.maximum_circuit_trials < 4 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                // Recheck and retain the actual trial source context while all
                // magnetic memory remains at the one accepted beginning state.
                let mut beginning = self.circuit.clone();
                if let Some(core) = core {
                    beginning.core = core;
                }
                let circuit = beginning
                    .step_candidate(
                        target_s,
                        coefficients.resistors,
                        external_inputs,
                        inner,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(magnetic)?;
                let circuit_trials = circuit.step().circuit_trials;
                Ok(DrivenCoefficients {
                    circuit,
                    circuit_trials,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                })
            },
            |candidate, dt, checkpoint| {
                let mut powers = BTreeMap::new();
                for loss in &candidate.step().circuit.resistors {
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
                let core_heat = candidate.step().core.core_heat_j;
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
pub struct HystereticEnthalpyCandidate<'t, 'n, 'a, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    circuit: ThermalHystereticCircuit<'n, 'a>,
    report: CircuitEnthalpyReport<S, HystereticCircuitStep>,
}
impl<'t, 'n, 'a, S> HystereticEnthalpyCandidate<'t, 'n, 'a, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn circuit(&self) -> &ThermalHystereticCircuit<'n, 'a> {
        &self.circuit
    }
    pub const fn report(
        &self,
    ) -> &CircuitEnthalpyReport<S, HystereticCircuitStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermalHystereticCircuit<'n, 'a>,
        CircuitEnthalpyReport<S, HystereticCircuitStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}
