//! Regulated demand and finite source impedance on one material/heat timeline.
//!
//! Both bus voltages, source currents and regulated demands participate in the
//! original nonlinear converter solve. Seven physical resistors and optional
//! fully dissipative load receivers heat their declared material bodies once.
//! These are held ideal voltage supplies behind finite resistance, not finite
//! energy reservoirs. Work uses the native discrete effort; temporal accuracy
//! requires refinement and no continuous RMS or PWM loss is inferred.

use super::super::battery::RegulatedLoadHeat;
use super::*;
use fs_phs::circuit::converter::power_load::{
    PowerLoadControls, RegulatedPowerLoad, VoltagePowerLoadedStep,
};

#[derive(Clone, Debug)]
pub struct VoltagePowerLoadedThermalSample {
    pub network: MaterialNetworkSample,
    pub channels: Arc<[ConverterChannelMaterial; VOLTAGE_DRIVEN_HEAT_CHANNELS]>,
    /// None exports that bus's load work. Some records its sole heat receiver.
    pub load_heat: [Option<RegulatedLoadHeat>; 2],
}

/// Complete proposed electrical and caloric histories, published together.
pub struct PowerLoadedVoltageConverterEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    converter: MaterialVoltageDrivenBuckBoost,
    report: CircuitEnthalpyReport<VoltagePowerLoadedThermalSample, VoltagePowerLoadedStep>,
}
impl<'t> PowerLoadedVoltageConverterEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn converter(&self) -> &MaterialVoltageDrivenBuckBoost {
        &self.converter
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<VoltagePowerLoadedThermalSample, VoltagePowerLoadedStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialVoltageDrivenBuckBoost,
        CircuitEnthalpyReport<VoltagePowerLoadedThermalSample, VoltagePowerLoadedStep>,
    ) {
        (self.thermal, self.converter, self.report)
    }
}

struct PowerDriver {
    converter: MaterialVoltageDrivenBuckBoost,
    step: VoltagePowerLoadedStep,
    heat: [f64; VOLTAGE_DRIVEN_HEAT_CHANNELS],
    load_heat: [Option<RegulatedLoadHeat>; 2],
    total_heat_j: f64,
    external_work_j: Dd,
}
impl ThermalDriver for PowerDriver {
    type Model = MaterialVoltageDrivenBuckBoost;
    type Step = VoltagePowerLoadedStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.converter.native.delta_h,
            supplied_j: self.external_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.converter.native.delta_h) - self.external_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.converter, self.step)
    }
}

// Audit the original converter and each source/resistor/load boundary BEFORE
// routing any work internally. No heat or supply is a balance remainder.
fn physical_power_account(
    step: &VoltagePowerLoadedStep,
    tolerance_j: f64,
) -> Result<([f64; VOLTAGE_DRIVEN_HEAT_CHANNELS], Dd, Dd), Error> {
    let heat = core::array::from_fn(|i| {
        if i < CONDUCTION_CHANNELS {
            step.converter.conduction.joules[i]
        } else {
            step.terminal_resistor_heat_j[i - CONDUCTION_CHANNELS]
        }
    });
    if heat
        .iter()
        .chain(&step.load_work_j)
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(invalid(
            "voltage-driven resistance and regulated load heats must be nonnegative",
        ));
    }
    let internal = sum(step.converter.conduction.joules)?;
    let total = sum(heat)?;
    let supply = sum(step.source_work_j)? - sum(step.load_work_j)?;
    let terminal = sum(step.converter.terminal_work_j)?;
    let delta = Dd::from_f64(finite(step.converter.native.delta_h)?);
    let dissipated = Dd::from_f64(finite(step.converter.native.dissipated)?);
    let native_supply = Dd::from_f64(finite(step.converter.native.supplied)?);
    let mut residuals = [
        internal - dissipated,
        terminal - native_supply,
        delta + internal - terminal,
        delta + total - supply,
        Dd::ZERO,
        Dd::ZERO,
    ];
    for port in 0..2 {
        residuals[4 + port] = Dd::from_f64(finite(step.source_work_j[port])?)
            - Dd::from_f64(step.terminal_resistor_heat_j[port])
            - Dd::from_f64(step.load_work_j[port])
            - Dd::from_f64(finite(step.converter.terminal_work_j[port])?);
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
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if converter.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || source_voltage_v.iter().any(|v| !v.is_finite())
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
        let cache = RefCell::new(
            None::<(
                [u64; VOLTAGE_DRIVEN_HEAT_CHANNELS],
                [f64; 3],
                VoltagePowerLoadedStep,
            )>,
        );
        let (thermal, converter, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(&cancel)?;
                if remaining == 0 {
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
                    let mut native = AveragedBuckBoost::new(parameters, modulation)
                        .map_err(model)?
                        .prepare_power_loaded_controlled(dt, electrical_controls, || {
                            (*cancel.borrow_mut())()
                        })
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
                            || (*cancel.borrow_mut())(),
                        )
                        .map_err(electrical)?;
                    (next, step)
                };
                let (heat, total, supply) =
                    physical_power_account(&step, controls.total_energy_tolerance_j)?;
                let load_heat: [Option<RegulatedLoadHeat>; 2] = core::array::from_fn(|port| {
                    load_heat_cells[port].map(|cell| RegulatedLoadHeat {
                        cell,
                        energy_j: step.load_work_j[port],
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
                    circuit_trials: 1,
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
