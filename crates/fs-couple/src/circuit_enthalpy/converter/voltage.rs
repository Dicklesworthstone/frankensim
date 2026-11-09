//! Voltage-driven conversion with seven explicit material heat destinations.
//!
//! Native source/load currents evolve with both bus capacitors throughout each
//! averaged or switched interval. The five converter losses and two terminal
//! resistor losses enter the original finite enthalpy residual exactly once.
//! Ideal-source work retains its sign, including regeneration. Resistance is
//! refreshed at each thermal endpoint trial; original flux/charges are retained.

/// Regulated loads and finite source impedance on the same material/heat timeline.
pub mod power_load;

use super::*;
use fs_math::dd::Dd;
use fs_phs::circuit::converter::{
    switching::voltage_driven::PreparedVoltageDrivenSwitchCycle, terminals::VoltageDrivenStep,
};
use fs_phs::linear_exact::DenseExponentialBudget;

/// Five native converter channels followed by the input/output terminal resistors.
pub const VOLTAGE_DRIVEN_HEAT_CHANNELS: usize = CONDUCTION_CHANNELS + 2;

/// Every physical resistor has one explicit receiving body and source law.
/// Several channels may share a body; their independently integrated heats are
/// added before deposition. Terminal channels 5/6 require strictly positive R.
pub struct VoltageDrivenConverterThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub channels: &'a [ConverterThermalChannel<'a>; VOLTAGE_DRIVEN_HEAT_CHANNELS],
}

/// Accepted native flux/charges, fixed L/C, and seven actual material snapshots.
#[derive(Clone, Debug)]
pub struct MaterialVoltageDrivenBuckBoost {
    time_s: f64,
    parameters: BuckBoostParameters,
    state: [f64; 3],
    channels: Arc<[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS]>,
}
impl MaterialVoltageDrivenBuckBoost {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    pub const fn state(&self) -> &[f64; 3] {
        &self.state
    }
    pub fn terminal_resistance_ohm(&self) -> [f64; 2] {
        [
            self.channels[5].resistance_ohm,
            self.channels[6].resistance_ohm,
        ]
    }
    pub fn channels(&self) -> &[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS] {
        &self.channels
    }
}

#[derive(Clone, Debug)]
pub struct VoltageDrivenConverterThermalSample {
    pub network: MaterialNetworkSample,
    pub channels: Arc<[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS]>,
}

pub struct VoltageDrivenConverterEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    converter: MaterialVoltageDrivenBuckBoost,
    report: CircuitEnthalpyReport<VoltageDrivenConverterThermalSample, VoltageDrivenStep>,
}
impl<'t> VoltageDrivenConverterEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn converter(&self) -> &MaterialVoltageDrivenBuckBoost {
        &self.converter
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<VoltageDrivenConverterThermalSample, VoltageDrivenStep> {
        &self.report
    }
    /// Publish both histories together. All seven reported resistor heats have
    /// already entered their thermal receivers and must not be deposited again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialVoltageDrivenBuckBoost,
        CircuitEnthalpyReport<VoltageDrivenConverterThermalSample, VoltageDrivenStep>,
    ) {
        (self.thermal, self.converter, self.report)
    }
}

pub(super) fn all_resistances(
    p: BuckBoostParameters,
    terminal: [f64; 2],
) -> [f64; VOLTAGE_DRIVEN_HEAT_CHANNELS] {
    let internal = resistances(p);
    [
        internal[0],
        internal[1],
        internal[2],
        internal[3],
        internal[4],
        terminal[0],
        terminal[1],
    ]
}
pub(super) fn resolved_parameters(
    mut parameters: BuckBoostParameters,
    channels: &[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS],
) -> Result<(BuckBoostParameters, [f64; 2]), Error> {
    parameters.winding_resistance_ohm = channels[0].resistance_ohm;
    parameters.input_on_resistance_ohm = [channels[1].resistance_ohm, channels[2].resistance_ohm];
    parameters.output_on_resistance_ohm = [channels[3].resistance_ohm, channels[4].resistance_ohm];
    let terminal = [channels[5].resistance_ohm, channels[6].resistance_ohm];
    if terminal
        .iter()
        .any(|r| !r.is_finite() || *r <= 0.0 || !r.recip().is_finite() || r.recip() <= 0.0)
    {
        return Err(invalid(
            "voltage-source and load resistances must be finite and positive",
        ));
    }
    Ok((parameters, terminal))
}
pub(super) fn sum(values: impl IntoIterator<Item = f64>) -> Result<Dd, Error> {
    let mut total = Dd::ZERO;
    for value in values {
        total = total + Dd::from_f64(finite(value)?);
    }
    finite(total.hi)?;
    finite(total.lo)?;
    Ok(total)
}

struct VoltageDriver {
    converter: MaterialVoltageDrivenBuckBoost,
    step: VoltageDrivenStep,
    heat: [f64; VOLTAGE_DRIVEN_HEAT_CHANNELS],
    total_heat_j: f64,
    source_work_j: Dd,
}
impl ThermalDriver for VoltageDriver {
    type Model = MaterialVoltageDrivenBuckBoost;
    type Step = VoltageDrivenStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.converter.native.delta_h,
            supplied_j: self.source_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.converter.native.delta_h) - self.source_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.converter, self.step)
    }
}

// Keep the converter boundary and complete source/load boundary independent.
// Neither thermal heat nor ideal-source work is an energy-balance remainder.
fn physical_account(
    step: &VoltageDrivenStep,
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
    {
        return Err(invalid(
            "negative or unrepresentable native converter/source/load heat",
        ));
    }
    let internal_heat = sum(step.converter.conduction.joules)?;
    let total_heat = sum(heat)?;
    let source_work = sum(step.source_work_j)?;
    let terminal_work = sum(step.converter.terminal_work_j)?;
    let delta = Dd::from_f64(finite(step.converter.native.delta_h)?);
    let dissipated = Dd::from_f64(finite(step.converter.native.dissipated)?);
    let supplied = Dd::from_f64(finite(step.converter.native.supplied)?);
    let residuals = [
        internal_heat - dissipated,
        terminal_work - supplied,
        delta + dissipated - supplied,
        delta + total_heat - source_work,
        Dd::from_f64(step.source_work_j[0])
            - Dd::from_f64(step.terminal_resistor_heat_j[0])
            - Dd::from_f64(step.converter.terminal_work_j[0]),
        Dd::from_f64(step.source_work_j[1])
            - Dd::from_f64(step.terminal_resistor_heat_j[1])
            - Dd::from_f64(step.converter.terminal_work_j[1]),
    ];
    if residuals
        .iter()
        .any(|r| !r.hi.is_finite() || !r.lo.is_finite() || r.abs().to_f64() > tolerance_j)
    {
        return Err(invalid(
            "native converter and ideal-source work must close with all seven heat paths",
        ));
    }
    Ok((heat, total_heat.to_f64(), source_work))
}

impl VoltageDrivenConverterThermalNetwork<'_> {
    pub(super) fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        }
    }

    /// Bind the original native state and material sources without resetting
    /// charge, flux or temperature. Fixed channels retain their numeric R;
    /// sourced channels replace it with their actual accepted-state property.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &AveragedBuckBoost,
        terminal_resistance_ohm: [f64; 2],
        state: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialVoltageDrivenBuckBoost, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if state.iter().any(|v| !v.is_finite()) || !native.system().hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "converter must retain a finite original flux/charge state",
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
            all_resistances(native.parameters(), terminal_resistance_ohm),
            None,
            false,
            &cancel,
        )?);
        let (parameters, _) = resolved_parameters(native.parameters(), &channels)?;
        AveragedBuckBoost::new(parameters, native.modulation()).map_err(model)?;
        check(&cancel)?;
        Ok(MaterialVoltageDrivenBuckBoost {
            time_s: thermal.time_s(),
            parameters,
            state,
            channels,
        })
    }

    /// Couple held ideal source voltages through their actual finite resistors.
    /// Zero volts gives a passive load; negative work is regeneration. Each
    /// trial starts at the same accepted flux/charges, resolves all seven
    /// R(T,phase) sources and integrates their heat on the native trajectory.
    ///
    /// The switched lane must span this thermal interval and charges every
    /// physical segment to the shared work budget. Source-domain, time, budget
    /// and cancellation refusals publish neither history. L/C and geometry
    /// remain fixed; no semiconductor transition/recovery loss is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialVoltageDrivenBuckBoost,
        drive: ConverterElectricalDrive<'_>,
        source_voltage_v: &[f64; 2],
        electrical_budget: LinearSplitBudget,
        heat_budget: DenseExponentialBudget,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<VoltageDrivenConverterEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if converter.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || source_voltage_v.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "voltage-driven converter and heat need the same positive interval and finite sources",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let fixed = all_resistances(converter.parameters, converter.terminal_resistance_ohm());
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
                "all seven converter source identities and resistances must match accepted heat state",
            ));
        }
        let segments = match drive {
            ConverterElectricalDrive::Averaged(_) => 1,
            ConverterElectricalDrive::Switching { intervals, .. } => intervals.len(),
        };
        if segments == 0 {
            return Err(invalid("converter requires a nonempty electrical interval"));
        }
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
                let (parameters, terminal) = resolved_parameters(converter.parameters, &channels)?;
                let (mut next, mut mean_voltage) = ([0.0; 3], [0.0; 2]);
                let step = match drive {
                    ConverterElectricalDrive::Averaged(modulation) => {
                        let mut native = AveragedBuckBoost::new(parameters, modulation)
                            .map_err(model)?
                            .prepare_voltage_driven_controlled(
                                terminal,
                                dt,
                                electrical_budget,
                                heat_budget,
                                || (*cancel.borrow_mut())(),
                            )
                            .map_err(electrical)?;
                        native
                            .step_controlled(
                                &converter.state,
                                source_voltage_v,
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
                        let mut native = PreparedVoltageDrivenSwitchCycle::new_controlled(
                            parameters,
                            terminal,
                            intervals,
                            electrical_budget,
                            heat_budget,
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
                                "complete voltage-driven switch cycle must span the same thermal interval",
                            ));
                        }
                        native
                            .step_controlled(
                                &converter.state,
                                source_voltage_v,
                                &mut next,
                                &mut mean_voltage,
                                || (*cancel.borrow_mut())(),
                            )
                            .map_err(electrical)?
                    }
                };
                let (heat, total_heat_j, source_work_j) =
                    physical_account(&step, controls.total_energy_tolerance_j)?;
                check(&cancel)?;
                Ok(DrivenCoefficients {
                    circuit: VoltageDriver {
                        converter: MaterialVoltageDrivenBuckBoost {
                            time_s: target_s,
                            parameters,
                            state: next,
                            channels: Arc::clone(&channels),
                        },
                        step,
                        heat,
                        total_heat_j,
                        source_work_j,
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
                            "converter/source/load heat is unrepresentable as a thermal rate",
                        ));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(VoltageDrivenConverterEnthalpyCandidate {
            thermal,
            converter,
            report,
        })
    }
}
