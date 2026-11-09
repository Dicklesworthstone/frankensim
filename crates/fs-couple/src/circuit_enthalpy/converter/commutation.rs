//! Source-bound channel, leakage and body-diode heat on finite material bodies.
//!
//! The original three-state commutation owner resolves every prescribed gate
//! interval. Its four device packages each have one original material source
//! and one junction heat receiver; channel and body-diode work remain separate
//! in the native report. Winding and supply resistors retain the existing bulk
//! or scalar resistance adapters. No switching loss is reconstructed from a
//! balance remainder, an efficiency, or the square of a mean current.

/// Finite chemical energy, both motor shafts and source-bound device/material heat.
pub mod battery_motor;

mod capacitors;
pub use capacitors::{DeviceCapacitanceBinding, DeviceCapacitanceMaterial};

mod charge;
pub use charge::{DeviceChargeBinding, DeviceChargeMaterial};

mod gate;
pub use gate::{
    GateChannelBinding, GateChannelMaterial, GateConductanceBinding, GateConductanceMaterial,
    GateCurrentSurfaceBinding, GateCurrentSurfaceMaterial, GateCurrentThermalScale, GateMaterial,
    GateMaterialBindings,
};

mod devices;
pub use devices::{DeviceConductionBinding, DeviceConductionKeys, DeviceConductionMaterial};

use super::battery::RegulatedLoadHeat;
use super::voltage::sum;
use super::*;
use devices::resolve_devices;
use fs_math::dd::Dd;
use fs_phs::circuit::converter::commutation::{
    BodyDiode, CommutatedStep, CommutationInterval, CommutationParameters, DeviceOperatingBounds,
    PreparedCommutatedSwitchCycle, SwitchDevice,
};
use fs_phs::circuit::converter::power_load::PowerLoadControls;
use fs_phs::circuit::converter::switching::power_load::PowerSwitchCommand;
use fs_qty::{Current, Voltage};

/// Three ordinary resistance owners: winding, input supply, output supply.
/// Four independently sourced device packages own their channel/diode heats.
pub struct CommutatedConverterThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub resistors: &'a [ConverterThermalChannel<'a>; 3],
    pub devices: &'a [DeviceConductionBinding<'a>; 4],
}

/// Accepted original flux/charges and complete material coefficient bundles.
#[derive(Clone, Debug)]
pub struct MaterialCommutatedConverter {
    time_s: f64,
    parameters: CommutationParameters,
    terminal_resistance_ohm: [f64; 2],
    state: [f64; 3],
    resistors: Arc<[ConverterChannelMaterial; 3]>,
    devices: Arc<[DeviceConductionMaterial; 4]>,
}
impl MaterialCommutatedConverter {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn parameters(&self) -> CommutationParameters {
        self.parameters
    }
    pub const fn terminal_resistance_ohm(&self) -> [f64; 2] {
        self.terminal_resistance_ohm
    }
    pub const fn state(&self) -> &[f64; 3] {
        &self.state
    }
    pub fn resistors(&self) -> &[ConverterChannelMaterial; 3] {
        &self.resistors
    }
    pub fn devices(&self) -> &[DeviceConductionMaterial; 4] {
        &self.devices
    }
}

#[derive(Clone, Debug)]
pub struct CommutatedThermalSample {
    pub network: MaterialNetworkSample,
    pub resistors: Arc<[ConverterChannelMaterial; 3]>,
    pub devices: Arc<[DeviceConductionMaterial; 4]>,
    /// None exports this load's actual work; Some deposits it internally once.
    pub load_heat: [Option<RegulatedLoadHeat>; 2],
}

/// A complete joint candidate. Neither accepted input history is mutated.
#[must_use]
pub struct CommutatedConverterEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    converter: MaterialCommutatedConverter,
    report: CircuitEnthalpyReport<CommutatedThermalSample, CommutatedStep>,
}
impl<'t> CommutatedConverterEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn converter(&self) -> &MaterialCommutatedConverter {
        &self.converter
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<CommutatedThermalSample, CommutatedStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialCommutatedConverter,
        CircuitEnthalpyReport<CommutatedThermalSample, CommutatedStep>,
    ) {
        (self.thermal, self.converter, self.report)
    }
}

struct CommutationDriver {
    converter: MaterialCommutatedConverter,
    step: CommutatedStep,
    load_heat: [Option<RegulatedLoadHeat>; 2],
    total_heat_j: f64,
    external_work_j: Dd,
}
impl ThermalDriver for CommutationDriver {
    type Model = MaterialCommutatedConverter;
    type Step = CommutatedStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.native.delta_h,
            supplied_j: self.external_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.native.delta_h) - self.external_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.converter, self.step)
    }
}

fn parameters_at(
    mut parameters: CommutationParameters,
    resistors: &[ConverterChannelMaterial; 3],
    devices: &[DeviceConductionMaterial; 4],
) -> Result<(CommutationParameters, [f64; 2]), Error> {
    parameters.winding_resistance_ohm = resistors[0].resistance_ohm();
    parameters.devices = core::array::from_fn(|i| devices[i].parameters());
    parameters.validate().map_err(model)?;
    let terminal = [resistors[1].resistance_ohm(), resistors[2].resistance_ohm()];
    if terminal.iter().any(|r| *r <= 0.0 || !(1.0 / r).is_finite()) {
        return Err(invalid(
            "commutated voltage supplies need finite positive series resistance",
        ));
    }
    Ok((parameters, terminal))
}

// Storage, gate/command history and accepted state are fixed within one call.
// Include every nonlinear coefficient AND limit: constant R with varying diode
// threshold or source qualification must never reuse a stale trajectory.
fn coefficient_bits(p: CommutationParameters, terminal: [f64; 2]) -> [u64; 31] {
    let mut bits = [0; 31];
    bits[0] = p.winding_resistance_ohm.to_bits();
    bits[1] = terminal[0].to_bits();
    bits[2] = terminal[1].to_bits();
    for (i, d) in p.devices.iter().enumerate() {
        let values = [
            d.on_resistance_ohm,
            d.off_resistance_ohm,
            d.body_diode.forward_voltage_v,
            d.body_diode.differential_resistance_ohm,
            d.bounds.maximum_reverse_voltage_v,
            d.bounds.maximum_channel_current_a,
            d.bounds.maximum_diode_current_a,
        ];
        for (j, value) in values.into_iter().enumerate() {
            bits[3 + 7 * i + j] = value.to_bits();
        }
    }
    bits
}

// Separate native, terminal, and complete source/storage balances precede all
// internal routing. Every positive heat term comes from its own constitutive law.
fn physical_account(step: &CommutatedStep, allowance_j: f64) -> Result<(Dd, Dd), Error> {
    let heats = [step.winding_heat_j]
        .into_iter()
        .chain(step.channel_heat_j)
        .chain(step.diode_heat_j)
        .chain(step.terminal_resistor_heat_j);
    if heats
        .clone()
        .chain(step.load_work_j)
        .any(|value| !value.is_finite() || value < 0.0)
    {
        return Err(invalid(
            "device, winding, supply and delivered-load heats must be nonnegative",
        ));
    }
    let internal = sum([step.winding_heat_j]
        .into_iter()
        .chain(step.channel_heat_j)
        .chain(step.diode_heat_j))?;
    let heat = sum(heats)?;
    let external = sum(step.source_work_j)? - sum(step.load_work_j)?;
    let terminal = sum(step.terminal_work_j)?;
    let delta = Dd::from_f64(finite(step.native.delta_h)?);
    let mut residuals = [
        internal - Dd::from_f64(finite(step.native.dissipated)?),
        terminal - Dd::from_f64(finite(step.native.supplied)?),
        delta + internal - terminal,
        delta + heat - external,
        Dd::ZERO,
        Dd::ZERO,
    ];
    for port in 0..2 {
        residuals[4 + port] = Dd::from_f64(finite(step.source_work_j[port])?)
            - Dd::from_f64(step.terminal_resistor_heat_j[port])
            - Dd::from_f64(step.load_work_j[port])
            - Dd::from_f64(finite(step.terminal_work_j[port])?);
    }
    if residuals
        .iter()
        .any(|r| !r.hi.is_finite() || !r.lo.is_finite() || r.abs().to_f64() > allowance_j)
    {
        return Err(invalid(
            "commutated device, source, load and storage work must close independently",
        ));
    }
    Ok((heat, external))
}

impl CommutatedConverterThermalNetwork<'_> {
    fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            links: self.links,
            resistors: &[],
        }
    }

    /// Bind original L/C and flux/charges to the actual material inventory.
    /// Device sources replace all seven authored package coefficients/limits;
    /// winding/supply resistance bindings retain their original interpretation.
    /// A gate program's operating points are admitted when that program steps.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        parameters: CommutationParameters,
        terminal_resistance_ohm: [f64; 2],
        state: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialCommutatedConverter, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if state.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "commutated conversion requires finite original flux and charges",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let coefficients = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let resistors = Arc::new(resolve_channels(
            self.cells,
            self.resistors,
            thermal,
            &coefficients.source,
            [
                parameters.winding_resistance_ohm,
                terminal_resistance_ohm[0],
                terminal_resistance_ohm[1],
            ],
            None,
            false,
            &cancel,
        )?);
        let devices = Arc::new(resolve_devices(
            self.cells,
            self.devices,
            thermal,
            &coefficients.source,
            None,
            &cancel,
        )?);
        let (parameters, terminal_resistance_ohm) =
            parameters_at(parameters, &resistors, &devices)?;
        check(&cancel)?;
        Ok(MaterialCommutatedConverter {
            time_s: thermal.time_s(),
            parameters,
            terminal_resistance_ohm,
            state,
            resistors,
            devices,
        })
    }

    /// Resolve prescribed independent gates, source/load commands and all finite
    /// thermal bodies at one endpoint. Sources are held ideal voltages behind
    /// finite resistances; actual load work may heat explicit receivers.
    ///
    /// Four original device bundles refresh on every thermal proposal, including
    /// off-state resistance, diode threshold/slope and every operating bound.
    /// Native initial/discrete/final checks apply on every electrical substep.
    /// Original source axes beyond the owned thermal state remain fixed; an
    /// unowned voltage/current/gate-dependent source cannot be frozen here.
    ///
    /// Every numerical substep is charged to the shared trial budget even on
    /// exact-coefficient cache reuse. All parameters use the proposed thermal
    /// endpoint over the electrical cycle: timestep refinement is required.
    /// This conduction-only owner has no junction/gate capacitance, recovery,
    /// switching-transition loss, avalanche, or inferred thermal impedance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialCommutatedConverter,
        intervals: &[CommutationInterval],
        commands: &[PowerSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cells: [Option<CellId>; 2],
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedConverterEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        let integration_steps = intervals
            .len()
            .checked_mul(substeps_per_interval)
            .filter(|steps| *steps > 0)
            .ok_or_else(|| invalid("commutated cycle requires bounded nonzero numerical work"))?;
        if converter.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || commands.len() != intervals.len()
        {
            return Err(invalid(
                "commutation and heat need the same positive interval and one command per physical gate interval",
            ));
        }
        for receiver in load_heat_cells.into_iter().flatten() {
            if !thermal.cells().contains_key(&receiver)
                || !self.cells.iter().any(|body| body.cell == receiver)
            {
                return Err(invalid(
                    "dissipative load requires an active original material receiver",
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
        let fixed = [
            converter.parameters.winding_resistance_ohm,
            converter.terminal_resistance_ohm[0],
            converter.terminal_resistance_ohm[1],
        ];
        let original_resistors = resolve_channels(
            self.cells,
            self.resistors,
            thermal,
            &accepted.source,
            fixed,
            Some(&converter.resistors),
            false,
            &cancel,
        )?;
        let original_devices = resolve_devices(
            self.cells,
            self.devices,
            thermal,
            &accepted.source,
            Some(&converter.devices),
            &cancel,
        )?;
        if original_resistors
            .iter()
            .zip(converter.resistors.iter())
            .any(|(a, b)| !same_channel(a, b))
            || original_devices != *converter.devices
        {
            return Err(invalid(
                "all commutated device and resistor sources must match the accepted heat state",
            ));
        }
        let cache = RefCell::new(None::<([u64; 31], [f64; 3], CommutatedStep)>);
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
                let resistors = Arc::new(resolve_channels(
                    self.cells,
                    self.resistors,
                    trial,
                    &coefficients.source,
                    fixed,
                    Some(&converter.resistors),
                    false,
                    &cancel,
                )?);
                let devices = Arc::new(resolve_devices(
                    self.cells,
                    self.devices,
                    trial,
                    &coefficients.source,
                    Some(&converter.devices),
                    &cancel,
                )?);
                let (parameters, terminal) =
                    parameters_at(converter.parameters, &resistors, &devices)?;
                let bits = coefficient_bits(parameters, terminal);
                let cached = cache.borrow().as_ref().filter(|c| c.0 == bits).cloned();
                let reused = cached.is_some();
                let (next, step) = if let Some((_, next, step)) = cached {
                    (next, step)
                } else {
                    let mut native = PreparedCommutatedSwitchCycle::new_controlled(
                        parameters,
                        terminal,
                        intervals,
                        substeps_per_interval,
                        electrical_controls,
                        cycle_controls,
                        || (*cancel.borrow_mut())(),
                    )
                    .map_err(electrical)?;
                    let period = native.period_s();
                    let represented_end = finite(thermal.time_s() + period)?;
                    let clock_allowance = 0.5
                        * (clock_spacing(thermal.time_s()) + clock_spacing(target_s))
                        + 64.0 * f64::EPSILON * period.max(dt);
                    if represented_end <= thermal.time_s()
                        || represented_end.to_bits() != target_s.to_bits()
                        || (period - dt).abs() > clock_allowance
                        || native.integration_steps() != integration_steps
                    {
                        return Err(invalid(
                            "complete commutation cycle must span the same represented thermal interval",
                        ));
                    }
                    let (mut next, mut mean_voltage) = ([0.0; 3], [0.0; 2]);
                    let step = native
                        .step_interval_commands_controlled(
                            &converter.state,
                            commands,
                            &mut next,
                            &mut mean_voltage,
                            || (*cancel.borrow_mut())(),
                        )
                        .map_err(electrical)?;
                    (next, step)
                };
                let (heat, supply) = physical_account(&step, controls.total_energy_tolerance_j)?;
                let load_heat = core::array::from_fn(|port| {
                    load_heat_cells[port].map(|cell| RegulatedLoadHeat {
                        cell,
                        energy_j: step.load_work_j[port],
                    })
                });
                let routed = sum(load_heat.iter().flatten().map(|load| load.energy_j))?;
                let total_heat_j = finite((heat + routed).to_f64())?;
                let external_work_j = supply + routed;
                finite(external_work_j.hi)?;
                finite(external_work_j.lo)?;
                check(&cancel)?;
                if !reused {
                    *cache.borrow_mut() = Some((bits, next, step.clone()));
                }
                Ok(DrivenCoefficients {
                    circuit: CommutationDriver {
                        converter: MaterialCommutatedConverter {
                            time_s: target_s,
                            parameters,
                            terminal_resistance_ohm: terminal,
                            state: next,
                            resistors: Arc::clone(&resistors),
                            devices: Arc::clone(&devices),
                        },
                        step,
                        load_heat,
                        total_heat_j,
                        external_work_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: CommutatedThermalSample {
                        network: coefficients.source,
                        resistors,
                        devices,
                        load_heat,
                    },
                    circuit_trials: integration_steps,
                })
            },
            |driver, duration, checkpoint| {
                let mut by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                let mut add = |cell, heat| {
                    let sum = by_cell.entry(cell).or_insert(Dd::ZERO);
                    *sum = *sum + Dd::from_f64(heat);
                };
                checkpoint()?;
                add(
                    driver.converter.resistors[0].cell(),
                    driver.step.winding_heat_j,
                );
                for i in 0..4 {
                    checkpoint()?;
                    let receiver = driver.converter.devices[i].cell();
                    add(receiver, driver.step.channel_heat_j[i]);
                    add(receiver, driver.step.diode_heat_j[i]);
                }
                for i in 0..2 {
                    checkpoint()?;
                    add(
                        driver.converter.resistors[i + 1].cell(),
                        driver.step.terminal_resistor_heat_j[i],
                    );
                }
                for load in driver.load_heat.iter().flatten() {
                    checkpoint()?;
                    add(load.cell, load.energy_j);
                }
                let mut powers = BTreeMap::new();
                for (cell, heat) in by_cell {
                    checkpoint()?;
                    let heat = finite(heat.to_f64())?;
                    let power = finite(heat / duration)?;
                    if heat > 0.0 && power == 0.0 {
                        return Err(invalid(
                            "commutated device heat is unrepresentable on this thermal interval",
                        ));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(CommutatedConverterEnthalpyCandidate {
            thermal,
            converter,
            report,
        })
    }
}
