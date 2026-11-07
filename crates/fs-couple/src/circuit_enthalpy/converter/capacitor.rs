//! Finite capacitor supplies and all seven losses on one thermal timeline.
//!
//! The native five-state electrical owner retains converter flux, both bus
//! charges and both external capacitor charges. Every thermal endpoint trial
//! restarts those same accepted coordinates, resolves seven source-bound
//! resistances at the trial temperature/phase, and advances the complete
//! averaged or prescribed PWM trajectory. No ideal voltage source replaces
//! either finite capacitor, and source currents charge the external capacitors.
//!
//! Electrical work, all five stored energies and seven independently integrated
//! resistor heats participate in the original finite enthalpy residual. Heat is
//! deposited once into each authored receiver; explicitly shared bodies sum
//! their channels before deposition. Failed trials publish neither history.
//!
//! L, all four C values and geometry remain fixed. Thermal coefficients use the
//! endpoint, so coupled time accuracy is generally first order and requires
//! refinement. This adds no battery chemistry, C(T) parameter work, capacitor
//! leakage or semiconductor switching/dead-time/recovery loss.

use super::voltage::{
    VOLTAGE_DRIVEN_HEAT_CHANNELS, VoltageDrivenConverterThermalNetwork,
    VoltageDrivenConverterThermalSample, all_resistances, resolved_parameters, sum,
};
use super::*;
use fs_math::dd::Dd;
use fs_phs::circuit::converter::{
    capacitor_terminals::{CapacitorFedBuckBoost, CapacitorFedStep, CapacitorTerminals},
    switching::capacitor_fed::PreparedCapacitorFedSwitchCycle,
};
use fs_phs::linear_exact::DenseExponentialBudget;

/// Accepted original five-state electrical storage and actual material sources.
/// External capacitances remain fixed; all seven R values follow their original
/// channel ownership, including capacitor series resistance and shared bodies.
#[derive(Clone, Debug)]
pub struct MaterialCapacitorFedBuckBoost {
    time_s: f64,
    parameters: BuckBoostParameters,
    terminal_capacitance_f: [f64; 2],
    state: [f64; 5],
    channels: Arc<[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS]>,
}
impl MaterialCapacitorFedBuckBoost {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    /// Original flux, bus charges, then external capacitor charges.
    pub const fn state(&self) -> &[f64; 5] {
        &self.state
    }
    pub fn terminals(&self) -> CapacitorTerminals {
        CapacitorTerminals {
            capacitance_f: self.terminal_capacitance_f,
            resistance_ohm: [
                self.channels[5].resistance_ohm,
                self.channels[6].resistance_ohm,
            ],
        }
    }
    pub fn channels(&self) -> &[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS] {
        &self.channels
    }
}

pub struct CapacitorFedConverterEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    converter: MaterialCapacitorFedBuckBoost,
    report: CircuitEnthalpyReport<VoltageDrivenConverterThermalSample, CapacitorFedStep>,
}
impl<'t> CapacitorFedConverterEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn converter(&self) -> &MaterialCapacitorFedBuckBoost {
        &self.converter
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<VoltageDrivenConverterThermalSample, CapacitorFedStep> {
        &self.report
    }
    /// Publish converter, both finite supplies and thermal history together.
    /// All seven reported resistor heats have already entered their receivers.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialCapacitorFedBuckBoost,
        CircuitEnthalpyReport<VoltageDrivenConverterThermalSample, CapacitorFedStep>,
    ) {
        (self.thermal, self.converter, self.report)
    }
}

struct CapacitorDriver {
    converter: MaterialCapacitorFedBuckBoost,
    step: CapacitorFedStep,
    heat: [f64; VOLTAGE_DRIVEN_HEAT_CHANNELS],
    total_heat_j: f64,
    external_work_j: Dd,
}

// One call fixes its accepted state, sources, clock, drive and storage. Exact
// equality of all seven R values therefore permits reuse of the SAME complete
// electrical trajectory, even when the surrounding thermal trial changes.
// Source resolution and thermal admission still run on every trial.
#[derive(Clone)]
struct ElectricalTrial {
    resistance_bits: [u64; VOLTAGE_DRIVEN_HEAT_CHANNELS],
    state: [f64; 5],
    step: CapacitorFedStep,
}
impl ThermalDriver for CapacitorDriver {
    type Model = MaterialCapacitorFedBuckBoost;
    type Step = CapacitorFedStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.total.delta_h,
            supplied_j: self.external_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.total.delta_h) - self.external_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.converter, self.step)
    }
}

// Every constituent keeps its original independent native work/heat integral.
// Complete electrical storage includes the terminal capacitors; treating their
// outflow as external supply would count their energy a second time.
fn physical_account(
    step: &CapacitorFedStep,
    tolerance_j: f64,
) -> Result<([f64; VOLTAGE_DRIVEN_HEAT_CHANNELS], f64, Dd), Error> {
    let heat = core::array::from_fn(|i| {
        if i < CONDUCTION_CHANNELS {
            step.converter.conduction.joules[i]
        } else {
            step.terminal_resistor_heat_j[i - CONDUCTION_CHANNELS]
        }
    });
    if heat.iter().any(|v| !v.is_finite() || *v < 0.0)
        || finite(step.converter.native.dissipated)? < 0.0
        || finite(step.total.dissipated)? < 0.0
    {
        return Err(invalid(
            "finite capacitor converter requires nonnegative physical heat",
        ));
    }
    let internal_heat = sum(step.converter.conduction.joules)?;
    let total_heat = sum(heat)?;
    let external_work = sum(step.external_work_j)?;
    let bus_work = sum(step.converter.terminal_work_j)?;
    let capacitor_change = sum(step.terminal_storage_change_j)?;
    let converter_change = Dd::from_f64(finite(step.converter.native.delta_h)?);
    let total_change = Dd::from_f64(finite(step.total.delta_h)?);
    let converter_heat = Dd::from_f64(step.converter.native.dissipated);
    let converter_supply = Dd::from_f64(finite(step.converter.native.supplied)?);
    let mut residuals = [Dd::ZERO; 11];
    residuals[0] = internal_heat - converter_heat;
    residuals[1] = bus_work - converter_supply;
    residuals[2] = converter_change + converter_heat - converter_supply;
    residuals[3] = total_heat - Dd::from_f64(step.total.dissipated);
    residuals[4] = external_work - Dd::from_f64(finite(step.total.supplied)?);
    residuals[5] = total_change + total_heat - external_work;
    residuals[6] = total_change - converter_change - capacitor_change;
    for p in 0..2 {
        let outflow = Dd::from_f64(finite(step.terminal_outflow_work_j[p])?);
        residuals[7 + p] = Dd::from_f64(step.terminal_storage_change_j[p]) + outflow
            - Dd::from_f64(step.external_work_j[p]);
        residuals[9 + p] = outflow
            - Dd::from_f64(step.converter.terminal_work_j[p])
            - Dd::from_f64(step.terminal_resistor_heat_j[p]);
    }
    if residuals
        .iter()
        .any(|r| !r.hi.is_finite() || !r.lo.is_finite() || r.abs().to_f64() > tolerance_j)
    {
        return Err(invalid(
            "converter, finite supplies and all seven losses must close independently",
        ));
    }
    Ok((heat, finite(total_heat.to_f64())?, external_work))
}

impl VoltageDrivenConverterThermalNetwork<'_> {
    /// Bind actual material R to the original finite-capacitor electrical owner.
    /// This shares the existing seven-channel network: five converter losses,
    /// then the input/output capacitor series resistors. No charge, flux or
    /// temperature resets when a material source replaces numeric resistance.
    pub fn bind_capacitor(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &CapacitorFedBuckBoost,
        state: [f64; 5],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialCapacitorFedBuckBoost, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if state.iter().any(|v| !v.is_finite()) || !native.system().hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "converter must retain finite original flux and all four charges",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let sample = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let channels = Arc::new(resolve_channels(
            self.cells,
            self.channels,
            thermal,
            &sample.source,
            all_resistances(native.parameters(), native.terminals().resistance_ohm),
            None,
            false,
            &cancel,
        )?);
        let (parameters, terminal_r) = resolved_parameters(native.parameters(), &channels)?;
        let terminals = CapacitorTerminals {
            capacitance_f: native.terminals().capacitance_f,
            resistance_ohm: terminal_r,
        };
        AveragedBuckBoost::new(parameters, native.modulation())
            .map_err(model)?
            .with_capacitor_terminals(terminals)
            .map_err(model)?;
        check(&cancel)?;
        Ok(MaterialCapacitorFedBuckBoost {
            time_s: thermal.time_s(),
            parameters,
            terminal_capacitance_f: terminals.capacitance_f,
            state,
            channels,
        })
    }

    /// Advance all five electrical coordinates and finite caloric bodies on one
    /// accepted timeline. Inputs are held currents INTO the external capacitors,
    /// not bus currents or prescribed capacitor voltages. Zero input closes the
    /// electrical boundary while finite stored energy continues to transfer.
    ///
    /// Each endpoint trial resolves all seven R(T,phase) sources and starts at
    /// the same accepted five-state storage. The native averaged/PWM owner keeps
    /// each independent work and heat integral. The complete switching cycle
    /// must span the thermal interval; work/domain/cancellation refusal publishes
    /// no converter, capacitor or thermal state. Original material claims stay
    /// pinned across trials and continuation, including off-state channels.
    /// One bounded cache reuses the electrical candidate when all seven R values
    /// are bit-identical. Each thermal trial still consumes conservative native
    /// segment-budget units; no cache entry survives this candidate call.
    #[allow(clippy::too_many_arguments)]
    pub fn step_capacitor_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialCapacitorFedBuckBoost,
        drive: ConverterElectricalDrive<'_>,
        external_current_a: &[f64; 2],
        electrical_budget: LinearSplitBudget,
        observation_budget: DenseExponentialBudget,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitorFedConverterEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if converter.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || external_current_a.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite capacitor converter and heat require the same positive interval and finite currents",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let fixed = all_resistances(converter.parameters, converter.terminals().resistance_ohm);
        let initial = resolve_channels(
            self.cells,
            self.channels,
            thermal,
            &accepted.source,
            fixed,
            None,
            false,
            &cancel,
        )?;
        if initial
            .iter()
            .zip(converter.channels.iter())
            .any(|(a, b)| !same_channel(a, b))
        {
            return Err(invalid(
                "all seven finite-supply source identities and resistances must match accepted heat",
            ));
        }
        let segments = match drive {
            ConverterElectricalDrive::Averaged(_) => 1,
            ConverterElectricalDrive::Switching { intervals, .. } => intervals.len(),
        };
        if segments == 0 {
            return Err(invalid(
                "finite capacitor converter requires a nonempty electrical interval",
            ));
        }
        let electrical_trial = RefCell::new(None::<ElectricalTrial>);
        let (thermal, converter, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(&cancel)?;
                if segments > remaining {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let channels = Arc::new(resolve_channels(
                    self.cells,
                    self.channels,
                    trial,
                    &coefficients.source,
                    fixed,
                    Some(&converter.channels),
                    false,
                    &cancel,
                )?);
                let (parameters, terminal_r) =
                    resolved_parameters(converter.parameters, &channels)?;
                let terminals = CapacitorTerminals {
                    capacitance_f: converter.terminal_capacitance_f,
                    resistance_ohm: terminal_r,
                };
                let resistance_bits = all_resistances(parameters, terminal_r).map(f64::to_bits);
                let reused = electrical_trial
                    .borrow()
                    .as_ref()
                    .filter(|old| old.resistance_bits == resistance_bits)
                    .cloned();
                let cache_miss = reused.is_none();
                let (next, step) = if let Some(old) = reused {
                    (old.state, old.step)
                } else {
                    let (mut next, mut mean_voltage) = ([0.0; 5], [0.0; 4]);
                    let step = match drive {
                        ConverterElectricalDrive::Averaged(modulation) => {
                            let mut native = AveragedBuckBoost::new(parameters, modulation)
                                .map_err(model)?
                                .with_capacitor_terminals(terminals)
                                .map_err(model)?
                                .prepare_controlled(
                                    dt,
                                    electrical_budget,
                                    observation_budget,
                                    || (*cancel.borrow_mut())(),
                                )
                                .map_err(electrical)?;
                            native
                                .step_controlled(
                                    &converter.state,
                                    external_current_a,
                                    &mut next,
                                    &mut mean_voltage,
                                    || (*cancel.borrow_mut())(),
                                )
                                .map_err(electrical)?
                        }
                        ConverterElectricalDrive::Switching {
                            intervals,
                            controls: switching,
                        } => {
                            let mut native = PreparedCapacitorFedSwitchCycle::new_controlled(
                                parameters,
                                terminals,
                                intervals,
                                electrical_budget,
                                observation_budget,
                                switching,
                                || (*cancel.borrow_mut())(),
                            )
                            .map_err(electrical)?;
                            let period = native.period_s();
                            let represented_end = finite(thermal.time_s() + period)?;
                            let clock_roundoff =
                                0.5 * (clock_spacing(thermal.time_s()) + clock_spacing(target_s));
                            let duration_roundoff = 64.0 * f64::EPSILON * period.max(dt);
                            if represented_end <= thermal.time_s()
                                || (period - dt).abs() > clock_roundoff + duration_roundoff
                            {
                                return Err(invalid(
                                    "complete capacitor-fed switch cycle must span the same thermal interval",
                                ));
                            }
                            native
                                .step_controlled(
                                    &converter.state,
                                    external_current_a,
                                    &mut next,
                                    &mut mean_voltage,
                                    || (*cancel.borrow_mut())(),
                                )
                                .map_err(electrical)?
                        }
                    };
                    (next, step)
                };
                let (heat, total_heat_j, external_work_j) =
                    physical_account(&step, controls.total_energy_tolerance_j)?;
                check(&cancel)?;
                if cache_miss {
                    *electrical_trial.borrow_mut() = Some(ElectricalTrial {
                        resistance_bits,
                        state: next,
                        step: step.clone(),
                    });
                }
                Ok(DrivenCoefficients {
                    circuit: CapacitorDriver {
                        converter: MaterialCapacitorFedBuckBoost {
                            time_s: target_s,
                            parameters,
                            terminal_capacitance_f: converter.terminal_capacitance_f,
                            state: next,
                            channels: Arc::clone(&channels),
                        },
                        step,
                        heat,
                        total_heat_j,
                        external_work_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: VoltageDrivenConverterThermalSample {
                        network: coefficients.source,
                        channels,
                    },
                    circuit_trials: segments,
                })
            },
            |candidate, duration, checkpoint| {
                let mut by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                for (channel, &heat) in candidate.converter.channels.iter().zip(&candidate.heat) {
                    checkpoint()?;
                    let total = by_cell.entry(channel.cell).or_insert(Dd::ZERO);
                    *total = *total + Dd::from_f64(heat);
                }
                let mut powers = BTreeMap::new();
                for (cell, heat) in by_cell {
                    checkpoint()?;
                    let heat = finite(heat.to_f64())?;
                    let power = finite(heat / duration)?;
                    if heat > 0.0 && power == 0.0 {
                        return Err(invalid(
                            "finite-supply heat is unrepresentable as a thermal rate",
                        ));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(CapacitorFedConverterEnthalpyCandidate {
            thermal,
            converter,
            report,
        })
    }
}
