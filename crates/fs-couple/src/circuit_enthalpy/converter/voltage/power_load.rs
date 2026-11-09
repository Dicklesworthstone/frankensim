//! Regulated demand and finite source impedance on one material/heat timeline.
//!
//! Both bus voltages, source currents and regulated demands participate in the
//! original nonlinear converter solve. Seven physical resistors and optional
//! fully dissipative load receivers heat their declared material bodies once.
//! These are held ideal voltage supplies behind finite resistance, not finite
//! energy reservoirs. Work uses the native discrete effort; temporal accuracy
//! requires refinement. Prescribed PWM retains its physical gate intervals;
//! no semiconductor transition, recovery or dead-time loss is inferred.

use super::super::battery::RegulatedLoadHeat;
use super::*;
use fs_phs::circuit::converter::power_load::{
    PowerLoadControls, RegulatedPowerLoad, VoltagePowerLoadedStep,
};
use fs_phs::circuit::converter::switching::power_load::{
    PowerLoadedSwitchCycleStep, PowerSwitchCommand, PreparedPowerLoadedSwitchCycle,
};

#[derive(Clone, Debug)]
pub struct VoltagePowerLoadedThermalSample {
    pub network: MaterialNetworkSample,
    pub channels: Arc<[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS]>,
    /// None exports that bus's load work. Some records its sole heat receiver.
    pub load_heat: [Option<RegulatedLoadHeat>; 2],
}

/// Complete proposed electrical and caloric histories, published together.
pub struct PowerLoadedVoltageConverterEnthalpyCandidate<'t, Step = VoltagePowerLoadedStep> {
    thermal: DepositionSimulation<'t>,
    converter: MaterialVoltageDrivenBuckBoost,
    report: CircuitEnthalpyReport<VoltagePowerLoadedThermalSample, Step>,
}
impl<'t, Step> PowerLoadedVoltageConverterEnthalpyCandidate<'t, Step> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn converter(&self) -> &MaterialVoltageDrivenBuckBoost {
        &self.converter
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<VoltagePowerLoadedThermalSample, Step> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialVoltageDrivenBuckBoost,
        CircuitEnthalpyReport<VoltagePowerLoadedThermalSample, Step>,
    ) {
        (self.thermal, self.converter, self.report)
    }
}

/// A complete switched electrical/caloric candidate with per-cycle work.
pub type PowerLoadedSwitchConverterEnthalpyCandidate<'t> =
    PowerLoadedVoltageConverterEnthalpyCandidate<'t, PowerLoadedSwitchCycleStep>;

// Both native electrical owners expose the same physical work channels.
// This adapter adds no trajectory, averaging, heat reconstruction or solver.
trait PowerStep: Clone {
    fn converter(&self) -> &fs_phs::circuit::converter::ConverterStep;
    fn source_work(&self) -> [f64; 2];
    fn terminal_heat(&self) -> [f64; 2];
    fn load_work(&self) -> [f64; 2];
}
impl PowerStep for VoltagePowerLoadedStep {
    fn converter(&self) -> &fs_phs::circuit::converter::ConverterStep {
        &self.converter
    }
    fn source_work(&self) -> [f64; 2] {
        self.source_work_j
    }
    fn terminal_heat(&self) -> [f64; 2] {
        self.terminal_resistor_heat_j
    }
    fn load_work(&self) -> [f64; 2] {
        self.load_work_j
    }
}
impl PowerStep for PowerLoadedSwitchCycleStep {
    fn converter(&self) -> &fs_phs::circuit::converter::ConverterStep {
        &self.converter
    }
    fn source_work(&self) -> [f64; 2] {
        self.source_work_j
    }
    fn terminal_heat(&self) -> [f64; 2] {
        self.terminal_resistor_heat_j
    }
    fn load_work(&self) -> [f64; 2] {
        self.load_work_j
    }
}

struct PowerDriver<Step> {
    converter: MaterialVoltageDrivenBuckBoost,
    step: Step,
    heat: [f64; VOLTAGE_DRIVEN_HEAT_CHANNELS],
    load_heat: [Option<RegulatedLoadHeat>; 2],
    total_heat_j: f64,
    external_work_j: Dd,
}
impl<Step: PowerStep> ThermalDriver for PowerDriver<Step> {
    type Model = MaterialVoltageDrivenBuckBoost;
    type Step = Step;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.converter().native.delta_h,
            supplied_j: self.external_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.converter().native.delta_h) - self.external_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.converter, self.step)
    }
}

// Audit the original converter and each source/resistor/load boundary BEFORE
// routing any work internally. No heat or supply is a balance remainder.
fn physical_power_account<Step: PowerStep>(
    step: &Step,
    tolerance_j: f64,
) -> Result<([f64; VOLTAGE_DRIVEN_HEAT_CHANNELS], Dd, Dd), Error> {
    let converter = step.converter();
    let source_work = step.source_work();
    let terminal_heat = step.terminal_heat();
    let load_work = step.load_work();
    let heat = core::array::from_fn(|i| {
        if i < CONDUCTION_CHANNELS {
            converter.conduction.joules[i]
        } else {
            terminal_heat[i - CONDUCTION_CHANNELS]
        }
    });
    if heat
        .iter()
        .chain(&load_work)
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(invalid(
            "voltage-driven resistance and regulated load heats must be nonnegative",
        ));
    }
    let internal = sum(converter.conduction.joules)?;
    let total = sum(heat)?;
    let supply = sum(source_work)? - sum(load_work)?;
    let terminal = sum(converter.terminal_work_j)?;
    let delta = Dd::from_f64(finite(converter.native.delta_h)?);
    let dissipated = Dd::from_f64(finite(converter.native.dissipated)?);
    let native_supply = Dd::from_f64(finite(converter.native.supplied)?);
    let mut residuals = [
        internal - dissipated,
        terminal - native_supply,
        delta + internal - terminal,
        delta + total - supply,
        Dd::ZERO,
        Dd::ZERO,
    ];
    for port in 0..2 {
        residuals[4 + port] = Dd::from_f64(finite(source_work[port])?)
            - Dd::from_f64(terminal_heat[port])
            - Dd::from_f64(load_work[port])
            - Dd::from_f64(finite(converter.terminal_work_j[port])?);
    }
    if residuals
        .iter()
        .any(|r| !r.hi.is_finite() || !r.lo.is_finite() || r.abs().to_f64() > tolerance_j)
    {
        return Err(invalid(
            "regulated voltage source, resistor, load and converter work must close independently",
        ));
    }
    Ok((heat, total, supply))
}

impl VoltageDrivenConverterThermalNetwork<'_> {
    /// Advance averaged conversion with two simultaneous finite-impedance
    /// voltage supplies and two optional regulated loads. Bind the original
    /// three electrical coordinates with the existing network's bind method.
    ///
    /// Each optional load receiver declares a wholly dissipative load. Its
    /// actual accepted work, including foldback/current limiting, becomes heat;
    /// curtailed demand supplies none. None keeps that load external. The two
    /// loads and resistance channels may share a body; their work sums once.
    /// An absent load still requires a valid receiver if one is declared.
    ///
    /// Every thermal trial resolves all seven original material resistances
    /// and solves source/load feedback from the same accepted flux/charges.
    /// Negative ideal-source work remains regeneration. Storage C/L, commands
    /// and geometry are held; no source capacity, efficiency or device loss is
    /// inferred. The native discrete-gradient work is not an exact continuous
    /// RMS integral. Unsupported material states, failure and cancellation
    /// publish neither history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_power_loaded_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialVoltageDrivenBuckBoost,
        modulation: BridgeModulation,
        source_voltage_v: &[f64; 2],
        loads: &[Option<RegulatedPowerLoad>; 2],
        load_heat_cells: [Option<CellId>; 2],
        electrical_controls: PowerLoadControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerLoadedVoltageConverterEnthalpyCandidate<'t>, Error> {
        if source_voltage_v.iter().any(|v| !v.is_finite()) {
            return Err(invalid("regulated voltage sources must be finite"));
        }
        self.step_power_driver_candidate(
            thermal,
            target_s,
            converter,
            load_heat_cells,
            controls,
            initial_energy_increments_j,
            1,
            |parameters, terminal, dt, checkpoint| {
                let mut native = AveragedBuckBoost::new(parameters, modulation)
                    .map_err(model)?
                    .prepare_power_loaded_controlled(dt, electrical_controls, &mut *checkpoint)
                    .map_err(electrical)?;
                let (mut next, mut effort) = ([0.0; 3], [0.0; 2]);
                let step = native
                    .step_voltage_driven_controlled(
                        &converter.state,
                        source_voltage_v,
                        &terminal,
                        loads,
                        &mut next,
                        &mut effort,
                        &mut *checkpoint,
                    )
                    .map_err(electrical)?;
                Ok((next, step))
            },
            cancelled,
        )
    }

    /// Couple actual prescribed gate intervals and their correlated source/load
    /// commands to the same seven material heat owners. Every command belongs
    /// to one ORIGINAL gate interval; numerical subdivisions add no new edge.
    ///
    /// Every thermal trial restarts all original flux/charges and advances the
    /// complete native nonlinear switching cycle. Resistor heat sums accepted
    /// squared-current work across substeps; it is never squared mean current
    /// or the heat of an averaged duty command. Actual regulated-load work may
    /// heat the declared receivers through the same internal-transfer account.
    ///
    /// The full authored cycle must span the thermal interval. Its native
    /// integration steps are charged to the shared circuit-trial budget even
    /// when identical material coefficients reuse a cached result. Source and
    /// load commands, gate program and subdivisions stay fixed within a call.
    /// All seven R(T,phase) values use the thermal trial endpoint; coupled time
    /// accuracy requires refinement. Ideal complementary gates introduce no
    /// diode, dead-time, transition, gate-charge or reverse-recovery physics.
    #[allow(clippy::too_many_arguments)]
    pub fn step_power_loaded_switching_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialVoltageDrivenBuckBoost,
        intervals: &[SwitchInterval],
        commands: &[PowerSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cells: [Option<CellId>; 2],
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerLoadedSwitchConverterEnthalpyCandidate<'t>, Error> {
        if commands.len() != intervals.len() {
            return Err(invalid(
                "each physical switch interval requires its own source/load command",
            ));
        }
        let integration_steps = intervals
            .len()
            .checked_mul(substeps_per_interval)
            .filter(|steps| *steps != 0)
            .ok_or_else(|| {
                invalid("switched regulated conversion requires representable nonzero work")
            })?;
        self.step_power_driver_candidate(
            thermal,
            target_s,
            converter,
            load_heat_cells,
            controls,
            initial_energy_increments_j,
            integration_steps,
            |parameters, terminal, dt, checkpoint| {
                let mut native = PreparedPowerLoadedSwitchCycle::new_controlled(
                    parameters,
                    terminal,
                    intervals,
                    substeps_per_interval,
                    electrical_controls,
                    cycle_controls,
                    &mut *checkpoint,
                )
                .map_err(electrical)?;
                let period = native.period_s();
                let represented_end = finite(thermal.time_s() + period)?;
                let clock_roundoff =
                    0.5 * (clock_spacing(thermal.time_s()) + clock_spacing(target_s));
                let duration_roundoff = 64.0 * f64::EPSILON * period.max(dt);
                if represented_end <= thermal.time_s()
                    || represented_end.to_bits() != target_s.to_bits()
                    || (period - dt).abs() > clock_roundoff + duration_roundoff
                {
                    return Err(invalid(
                        "the complete regulated switch cycle must span the thermal interval",
                    ));
                }
                let (mut next, mut mean_voltage) = ([0.0; 3], [0.0; 2]);
                let step = native
                    .step_interval_commands_controlled(
                        &converter.state,
                        commands,
                        &mut next,
                        &mut mean_voltage,
                        &mut *checkpoint,
                    )
                    .map_err(electrical)?;
                Ok((next, step))
            },
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_power_driver_candidate<'t, Step: PowerStep>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialVoltageDrivenBuckBoost,
        load_heat_cells: [Option<CellId>; 2],
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        integration_steps: usize,
        integrate: impl Fn(
            BuckBoostParameters,
            [f64; 2],
            f64,
            &mut dyn FnMut() -> bool,
        ) -> Result<([f64; 3], Step), Error>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerLoadedVoltageConverterEnthalpyCandidate<'t, Step>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if converter.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || integration_steps == 0
        {
            return Err(invalid(
                "regulated voltage conversion and heat need one positive interval and finite sources",
            ));
        }
        for receiver in load_heat_cells.into_iter().flatten() {
            if !thermal.cells().contains_key(&receiver)
                || !self.cells.iter().any(|binding| binding.cell == receiver)
            {
                return Err(invalid(
                    "dissipative regulated load requires an active material heat receiver",
                ));
            }
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
                "all seven voltage/load resistance sources must match accepted heat state",
            ));
        }
        // Commands, dt, storage and accepted state are fixed within this call.
        // Exact resistance equality reuses the same native electrical trajectory;
        // caloric/material admission still runs on every thermal trial.
        let cache = RefCell::new(None::<([u64; VOLTAGE_DRIVEN_HEAT_CHANNELS], [f64; 3], Step)>);
        let (thermal, converter, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(&cancel)?;
                if integration_steps > remaining {
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
                let bits = all_resistances(parameters, terminal).map(f64::to_bits);
                let previous = cache
                    .borrow()
                    .as_ref()
                    .filter(|entry| entry.0 == bits)
                    .cloned();
                let reused = previous.is_some();
                let (next, step) = if let Some((_, next, step)) = previous {
                    (next, step)
                } else {
                    integrate(parameters, terminal, dt, &mut || (*cancel.borrow_mut())())?
                };
                let (heat, total, supply) =
                    physical_power_account(&step, controls.total_energy_tolerance_j)?;
                let load_heat: [Option<RegulatedLoadHeat>; 2] = core::array::from_fn(|port| {
                    load_heat_cells[port].map(|cell| RegulatedLoadHeat {
                        cell,
                        energy_j: step.load_work()[port],
                    })
                });
                let routed = sum(load_heat.iter().flatten().map(|load| load.energy_j))?;
                let total_heat_j = finite((total + routed).to_f64())?;
                let external_work_j = supply + routed;
                finite(external_work_j.hi)?;
                finite(external_work_j.lo)?;
                check(&cancel)?;
                if !reused {
                    *cache.borrow_mut() = Some((bits, next, step.clone()));
                }
                Ok(DrivenCoefficients {
                    circuit: PowerDriver {
                        converter: MaterialVoltageDrivenBuckBoost {
                            time_s: target_s,
                            parameters,
                            state: next,
                            channels: Arc::clone(&channels),
                        },
                        step,
                        heat,
                        load_heat,
                        total_heat_j,
                        external_work_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: VoltagePowerLoadedThermalSample {
                        network: coefficients.source,
                        channels,
                        load_heat,
                    },
                    circuit_trials: integration_steps,
                })
            },
            |candidate, duration, checkpoint| {
                let mut by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                for (channel, &heat) in candidate.converter.channels.iter().zip(&candidate.heat) {
                    checkpoint()?;
                    let total = by_cell.entry(channel.cell).or_insert(Dd::ZERO);
                    *total = *total + Dd::from_f64(heat);
                }
                for load in candidate.load_heat.iter().flatten() {
                    checkpoint()?;
                    let total = by_cell.entry(load.cell).or_insert(Dd::ZERO);
                    *total = *total + Dd::from_f64(load.energy_j);
                }
                let mut powers = BTreeMap::new();
                for (cell, heat) in by_cell {
                    checkpoint()?;
                    let heat = finite(heat.to_f64())?;
                    let power = finite(heat / duration)?;
                    if heat > 0.0 && power == 0.0 {
                        return Err(invalid(
                            "regulated converter heat is unrepresentable as a thermal rate",
                        ));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(PowerLoadedVoltageConverterEnthalpyCandidate {
            thermal,
            converter,
            report,
        })
    }
}
