//! Polarization switching and native resistor heat on the finite enthalpy timeline.
//!
//! Every proposed temperature restarts the SAME accepted charge, polarization
//! and native circuit state. The existing material network supplies caloric
//! storage, R(T) and finite conduction; complete HeatAssembly providers can add
//! thermal contacts, layers, radiation and finite coolant. The existing cell
//! driver owns the thermal residual, nonlinear solve and physical energy gates.
//! Switching and Joule work enter their original receiving cells exactly once.
//! An explicitly authored leakage resistor may share the dielectric volume.
//!
//! Geometry and switching coefficients remain fixed. Source-bound fits require
//! the original caloric card, complete query and explicit temperature coverage
//! for every coefficient. This is not a temperature-dependent switching or
//! electrocaloric law. Endpoint transport and backward-Euler heat give generally
//! first-order coupled accuracy. Monotone-charge cuts and time refinement remain
//! necessary, including when thermal feedback changes the electrical trajectory.

use super::{
    DielectricCircuitCandidate, DielectricCircuitSolve, DielectricCircuitStep,
    DielectricHystereticCircuit, HystereticCircuitError,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, DriverAccount, ThermalDriver, cell_driver_candidate,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
    material::{
        MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
        PreparedMaterialNetwork,
    },
};
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

/// Typed original source, electrical or complete heat-provider refusal.
#[derive(Debug)]
pub enum DielectricHystereticThermalError<E = HeatNetworkError> {
    Material(MaterialNetworkError),
    Electrical(HystereticCircuitError),
    Heat(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for DielectricHystereticThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "dielectric hysteresis/enthalpy: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for DielectricHystereticThermalError<E> {}
type Error<E = HeatNetworkError> = CircuitEnthalpyError<DielectricHystereticThermalError<E>>;

fn invalid<E>(s: &'static str) -> Error<E> {
    CircuitEnthalpyError::Invalid(s)
}
fn material<E>(e: MaterialNetworkError) -> Error<E> {
    match e {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(DielectricHystereticThermalError::Material(other)),
    }
}
fn electrical<E>(e: HystereticCircuitError) -> Error<E> {
    match e {
        HystereticCircuitError::Cancelled => CircuitEnthalpyError::Cancelled,
        HystereticCircuitError::TrialBudget => CircuitEnthalpyError::CircuitBudget,
        other => CircuitEnthalpyError::Source(DielectricHystereticThermalError::Electrical(other)),
    }
}
fn heat<E: From<HeatNetworkError>>(e: HeatNetworkError) -> Error<E> {
    match e {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(DielectricHystereticThermalError::Heat(other.into())),
    }
}

impl ThermalDriver for DielectricCircuitCandidate {
    type Model = DielectricHystereticCircuit;
    type Step = DielectricCircuitStep;
    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &step.circuit.resistors,
            contacts: &step.circuit.contacts,
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            stored_change_j: step.stored_change_j,
            supplied_j: step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        DielectricCircuitCandidate::into_parts(self)
    }
}

/// One accepted electrical history with pinned original caloric/transport sources.
pub struct ThermalDielectricHystereticCircuit<'n, 'a> {
    circuit: DielectricHystereticCircuit,
    network: &'n MaterialNetwork<'a>,
    prepared: Arc<PreparedMaterialNetwork<'n, 'a>>,
    destinations: Arc<BTreeMap<usize, CellId>>,
    dielectric_cell: CellId,
    density_kg_m3: f64,
    temperature_interval_k: [f64; 2],
    source_binding: Option<&'n MaterialCell<'a>>,
}

impl<'n, 'a> ThermalDielectricHystereticCircuit<'n, 'a> {
    /// Bind the entire fixed solid capacitor volume to its actual finite cell.
    /// A native source-bound resistor can dissipate into this same cell; its
    /// Ohmic work and the separate switching work are accumulated, never replaced.
    /// Authored fits keep their caller-admitted temperature interval. Source fits
    /// additionally require complete original-card and selected-claim coverage.
    #[allow(clippy::too_many_arguments)]
    pub fn bind(
        circuit: DielectricHystereticCircuit,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        dielectric_cell: CellId,
        temperature_interval_k: [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::bind_with_link_conductances(
            circuit,
            thermal,
            network,
            electrical_bases,
            dielectric_cell,
            temperature_interval_k,
            &BTreeMap::new(),
            cancelled,
        )
    }

    /// Admit additional finite links owned by a complete heat assembly. The
    /// initial link map and original material links must partition the complete
    /// thermal graph. Continue through `step_candidate_with_heat_assembly` so
    /// those extra links are resolved at every proposed temperature; ordinary
    /// stepping refuses unowned links instead of freezing initial conductances.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_with_link_conductances(
        circuit: DielectricHystereticCircuit,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        dielectric_cell: CellId,
        temperature_interval_k: [f64; 2],
        initial_link_conductances_w_k: &BTreeMap<u64, f64>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        let [low, high] = temperature_interval_k;
        if circuit.time_s().to_bits() != thermal.time_s().to_bits()
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
        {
            return Err(invalid(
                "matching clocks and a positive dielectric temperature interval are required",
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
            .find(|binding| binding.cell == dielectric_cell)
            .ok_or_else(|| invalid("dielectric lacks its actual caloric/material binding"))?;
        let density_kg_m3 = binding
            .curve
            .knots()
            .first()
            .ok_or_else(|| invalid("empty dielectric caloric chart"))?
            .bulk_density_kg_m3;
        if !density_kg_m3.is_finite()
            || density_kg_m3 <= 0.0
            || binding.curve.knots().iter().any(|knot| {
                knot.liquid_mass_fraction != 0.0
                    || (knot.bulk_density_kg_m3 / density_kg_m3 - 1.0).abs() > 1e-10
            })
        {
            return Err(invalid(
                "fixed dielectric geometry requires a wholly solid constant-density caloric chart",
            ));
        }
        let mut initial = prepared
            .resolve_with_link_conductances(thermal, initial_link_conductances_w_k, &mut cancelled)
            .map_err(material)?;
        initial.resistors.sort_by_key(|resistor| resistor.branch);
        if initial.resistors != circuit.native_circuit().resistors() {
            return Err(invalid(
                "accepted circuit resistor sources must match the actual thermal endpoint",
            ));
        }
        let source_binding = if let Some(source) = circuit.capacitor().source() {
            let cell = thermal
                .cells()
                .get(&dielectric_cell)
                .ok_or_else(|| invalid("source-bound dielectric cell is inactive"))?;
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
                    "dielectric fit and caloric cell require the same original card and complete accepted query",
                ));
            }
            source
                .admit_temperature_interval(binding.point.temperature_axis, temperature_interval_k)
                .map_err(HystereticCircuitError::from)
                .map_err(electrical)?;
            Some(binding)
        } else {
            None
        };
        let destinations = network
            .resistors
            .iter()
            .map(|resistor| (resistor.destination.branch, resistor.destination.cell))
            .collect();
        let result = Self {
            circuit,
            network,
            prepared: Arc::new(prepared),
            destinations: Arc::new(destinations),
            dielectric_cell,
            density_kg_m3,
            temperature_interval_k,
            source_binding,
        };
        let _ = result.dielectric_query::<HeatNetworkError>(thermal)?;
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(result)
    }

    pub const fn circuit(&self) -> &DielectricHystereticCircuit {
        &self.circuit
    }
    pub const fn dielectric_cell(&self) -> CellId {
        self.dielectric_cell
    }
    pub const fn temperature_interval_k(&self) -> [f64; 2] {
        self.temperature_interval_k
    }

    fn dielectric_query<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<Option<QueryPoint>, Error<E>> {
        let cell = thermal
            .cells()
            .get(&self.dielectric_cell)
            .ok_or_else(|| invalid("dielectric heat cell is inactive"))?;
        let volume = self.circuit.capacitor().volume_m3();
        let mass_volume = cell.mass_kg / self.density_kg_m3;
        let temperature = thermal.temperature_k(self.dielectric_cell)?;
        if !mass_volume.is_finite()
            || (cell.geometry.volume_m3 / volume - 1.0).abs() > 1e-10
            || (mass_volume / volume - 1.0).abs() > 1e-10
            || temperature < self.temperature_interval_k[0]
            || temperature > self.temperature_interval_k[1]
        {
            return Err(invalid(
                "dielectric must retain its actual fixed solid volume and admitted temperature interval",
            ));
        }
        self.source_binding
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

    /// Propose native charge/polarization and finite-cell enthalpies together.
    /// All native resistor heat and switching heat are already deposited once
    /// in the returned thermal state. Neither accepted input is mutated.
    /// Every electrical evaluation counts against both nested work budgets,
    /// including discarded thermal stencils and electrical reversal probes.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        electrical_controls: DielectricCircuitSolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DielectricHystereticEnthalpyCandidate<'t, 'n, 'a>, Error> {
        let cancel = RefCell::new(cancelled);
        self.step_resolved_candidate(
            thermal,
            target_s,
            external_inputs,
            controls,
            electrical_controls,
            initial_energy_increments_j,
            |state| {
                self.prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            || (*cancel.borrow_mut())(),
        )
    }

    /// Use the original bulk network with its complete native heat assembly.
    /// Contacts, layers, radiation and finite coolant see the same proposed
    /// temperatures as R(T) and switching heat. Source receipts remain in the
    /// accepted report. Closed radiation is gated as internal exchange.
    /// Electrical contact branches and field-dependent thermal boundaries are
    /// refused because this electrical lane does not own their additional work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_heat_assembly<'t, A: HeatAssembly + ?Sized>(
        &self,
        cx: &Cx<'_>,
        assembly: &A,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        electrical_controls: DielectricCircuitSolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<DielectricHystereticEnthalpyCandidate<'t, 'n, 'a, A::Sample>, Error<A::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let prepared = assembly.prepared_heat_network();
        let network = prepared
            .field_material_network()
            .map_err(heat::<A::Error>)?;
        if !std::ptr::eq(network, self.network) {
            return Err(invalid(
                "dielectric heat assembly must retain the original bulk material network",
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
            electrical_controls,
            initial_energy_increments_j,
            |state| {
                let coefficients = assembly
                    .resolve_heat(cx, state)
                    .map_err(|error| error.map_source(DielectricHystereticThermalError::Heat))?;
                if !coefficients.contacts.is_empty() {
                    return Err(invalid(
                        "dielectric circuit cannot discard electrical contact work or heat",
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

    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        electrical_controls: DielectricCircuitSolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<CircuitThermalCoefficients<S>, Error<E>>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DielectricHystereticEnthalpyCandidate<'t, 'n, 'a, S>, Error<E>> {
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
            return Err(invalid("accepted dielectric and thermal clocks must agree"));
        }
        if let Some(query) = self.dielectric_query(thermal)? {
            if self
                .circuit
                .capacitor()
                .source()
                .map(|s| s.material().query())
                != Some(query)
            {
                return Err(invalid(
                    "accepted dielectric source query must match the actual thermal state",
                ));
            }
        }
        let mut initial = resolve(thermal)?;
        initial.resistors.sort_by_key(|resistor| resistor.branch);
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
                let query = self.dielectric_query(state)?;
                let coefficients = resolve(state)?;
                let mut inner = electrical_controls;
                inner.circuit = controls.circuit;
                inner.maximum_circuit_trials = inner.maximum_circuit_trials.min(remaining);
                if inner.maximum_circuit_trials < 4 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let beginning = match query {
                    Some(query) => DielectricHystereticCircuit {
                        capacitor: self
                            .circuit
                            .capacitor()
                            .at_source_query(&query)
                            .map_err(HystereticCircuitError::from)
                            .map_err(electrical)?,
                        ..self.circuit.clone()
                    },
                    None => self.circuit.clone(),
                };
                let circuit = beginning
                    .step_candidate(
                        target_s,
                        coefficients.resistors,
                        external_inputs,
                        inner,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(electrical)?;
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
                            invalid("native resistor heat lacks its original receiving cell")
                        })?;
                    add_heat(&mut powers, cell, loss.heat_j, dt)?;
                }
                checkpoint()?;
                add_heat(
                    &mut powers,
                    self.dielectric_cell,
                    candidate.step().dielectric.switching_heat_j,
                    dt,
                )?;
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check()?;
        Ok(DielectricHystereticEnthalpyCandidate {
            thermal,
            circuit: Self {
                circuit,
                network: self.network,
                prepared: Arc::clone(&self.prepared),
                destinations: Arc::clone(&self.destinations),
                dielectric_cell: self.dielectric_cell,
                density_kg_m3: self.density_kg_m3,
                temperature_interval_k: self.temperature_interval_k,
                source_binding: self.source_binding,
            },
            report,
        })
    }
}

fn add_heat<E>(
    powers: &mut BTreeMap<CellId, f64>,
    cell: CellId,
    heat_j: f64,
    dt_s: f64,
) -> Result<(), Error<E>> {
    let power = heat_j / dt_s;
    let previous = powers.get(&cell).copied().unwrap_or(0.0);
    let total = previous + power;
    if !heat_j.is_finite()
        || heat_j < 0.0
        || !power.is_finite()
        || (heat_j > 0.0 && power == 0.0)
        || !total.is_finite()
        || (power > 0.0 && total == previous)
    {
        return Err(invalid(
            "combined physical dielectric heat is unrepresentable",
        ));
    }
    powers.insert(cell, total);
    Ok(())
}

/// Both publishable histories and their complete physical energy account.
/// The native report's legacy `deposited_joule_j` field counts total physical
/// deposition here. Switching and resistor receipts remain separate under
/// `report().electrical`; numerical dissipation is never deposited as heat.
pub struct DielectricHystereticEnthalpyCandidate<'t, 'n, 'a, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    circuit: ThermalDielectricHystereticCircuit<'n, 'a>,
    report: CircuitEnthalpyReport<S, DielectricCircuitStep>,
}
impl<'t, 'n, 'a, S> DielectricHystereticEnthalpyCandidate<'t, 'n, 'a, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn circuit(&self) -> &ThermalDielectricHystereticCircuit<'n, 'a> {
        &self.circuit
    }
    pub const fn report(&self) -> &CircuitEnthalpyReport<S, DielectricCircuitStep> {
        &self.report
    }
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.deposited_joule_j
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermalDielectricHystereticCircuit<'n, 'a>,
        CircuitEnthalpyReport<S, DielectricCircuitStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}
