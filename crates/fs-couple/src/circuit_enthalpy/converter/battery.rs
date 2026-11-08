//! Finite chemical storage, regulated conversion and heat on one timeline.
//!
//! Every thermal endpoint trial starts from the same accepted six electrical
//! coordinates and resolves all eight original resistance sources. Battery
//! series/polarization heat, the five converter channels and physical lead heat
//! enter their authored receivers exactly once. Regulated-load work leaves the
//! electrical boundary; it is never relabeled as converter or battery heat.
//!
//! OCV, capacity, inductance and capacitances stay fixed. Resistance uses the
//! thermal endpoint, so coupled temporal accuracy requires refinement. This
//! owner adds no OCV(T), reversible entropic battery heat, chemistry fit, aging,
//! PWM switching loss or protection-event model.

/// Finite winding, rotor and stator storage on the same battery/thermal timeline.
pub mod motor;

use super::voltage::sum;
use super::*;
use fs_math::dd::Dd;
use fs_phs::circuit::{
    battery::{BatteryOcvCurve, BatteryParameters, BatterySource},
    converter::{
        battery_fed::{BatteryFedBuckBoost, BatteryFedStep},
        power_load::{PowerLoadControls, RegulatedPowerLoad},
    },
};

/// Converter winding/high/low channels, battery series/polarization, then lead.
pub const BATTERY_FED_HEAT_CHANNELS: usize = CONDUCTION_CHANNELS + 3;

/// Eight physical resistors, with explicitly shared or distinct receiving cells.
pub struct BatteryFedConverterThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub channels: &'a [ConverterThermalChannel<'a>; BATTERY_FED_HEAT_CHANNELS],
}

/// Accepted chemical/polarization/battery-bus charge, converter flux and buses.
#[derive(Clone, Debug)]
pub struct MaterialBatteryFedBuckBoost {
    time_s: f64,
    curve: BatteryOcvCurve,
    battery_parameters: BatteryParameters,
    parameters: BuckBoostParameters,
    modulation: BridgeModulation,
    lead_resistance_ohm: f64,
    state: [f64; 6],
    channels: Arc<[ConverterChannelMaterial; BATTERY_FED_HEAT_CHANNELS]>,
}
impl MaterialBatteryFedBuckBoost {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        &self.curve
    }
    pub const fn battery_parameters(&self) -> BatteryParameters {
        self.battery_parameters
    }
    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    pub const fn modulation(&self) -> BridgeModulation {
        self.modulation
    }
    pub const fn lead_resistance_ohm(&self) -> f64 {
        self.lead_resistance_ohm
    }
    pub const fn state(&self) -> &[f64; 6] {
        &self.state
    }
    pub fn channels(&self) -> &[ConverterChannelMaterial; BATTERY_FED_HEAT_CHANNELS] {
        &self.channels
    }
}

#[derive(Clone, Debug)]
pub struct BatteryFedConverterThermalSample {
    pub network: MaterialNetworkSample,
    pub channels: Arc<[ConverterChannelMaterial; BATTERY_FED_HEAT_CHANNELS]>,
}

/// Publish the complete electrical and caloric histories together.
pub struct BatteryFedConverterEnthalpyCandidate<'t> {
    thermal: DepositionSimulation<'t>,
    converter: MaterialBatteryFedBuckBoost,
    report: CircuitEnthalpyReport<BatteryFedConverterThermalSample, BatteryFedStep>,
}
impl<'t> BatteryFedConverterEnthalpyCandidate<'t> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn converter(&self) -> &MaterialBatteryFedBuckBoost {
        &self.converter
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<BatteryFedConverterThermalSample, BatteryFedStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialBatteryFedBuckBoost,
        CircuitEnthalpyReport<BatteryFedConverterThermalSample, BatteryFedStep>,
    ) {
        (self.thermal, self.converter, self.report)
    }
}

fn all_resistances(p: BuckBoostParameters, b: BatteryParameters, lead: f64) -> [f64; 8] {
    let r = resistances(p);
    [
        r[0],
        r[1],
        r[2],
        r[3],
        r[4],
        b.series_resistance_ohm,
        b.polarization_resistance_ohm,
        lead,
    ]
}
fn parameters_at(
    mut p: BuckBoostParameters,
    mut b: BatteryParameters,
    channels: &[ConverterChannelMaterial; BATTERY_FED_HEAT_CHANNELS],
) -> (BuckBoostParameters, BatteryParameters, f64) {
    p.winding_resistance_ohm = channels[0].resistance_ohm;
    p.input_on_resistance_ohm = [channels[1].resistance_ohm, channels[2].resistance_ohm];
    p.output_on_resistance_ohm = [channels[3].resistance_ohm, channels[4].resistance_ohm];
    b.series_resistance_ohm = channels[5].resistance_ohm;
    b.polarization_resistance_ohm = channels[6].resistance_ohm;
    (p, b, channels[7].resistance_ohm)
}
fn native(
    curve: &BatteryOcvCurve,
    p: BuckBoostParameters,
    b: BatteryParameters,
    modulation: BridgeModulation,
    lead: f64,
) -> Result<BatteryFedBuckBoost, Error> {
    let battery = BatterySource::new(curve.clone(), b).map_err(model)?;
    AveragedBuckBoost::new(p, modulation)
        .map_err(model)?
        .with_battery_source(battery, lead)
        .map_err(model)
}

struct BatteryDriver {
    converter: MaterialBatteryFedBuckBoost,
    step: BatteryFedStep,
    heat: [f64; BATTERY_FED_HEAT_CHANNELS],
    total_heat_j: f64,
    net_external_work_j: Dd,
}
impl ThermalDriver for BatteryDriver {
    type Model = MaterialBatteryFedBuckBoost;
    type Step = BatteryFedStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.native.delta_h,
            supplied_j: self.net_external_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.native.delta_h) - self.net_external_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.converter, self.step)
    }
}

// Independent constituent balances prevent offsetting local work errors from
// disappearing in the whole electrical/thermal accounting window.
fn physical_account(step: &BatteryFedStep, tolerance_j: f64) -> Result<([f64; 8], Dd, Dd), Error> {
    let heat = step.heat_joules();
    if heat.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(invalid(
            "battery-fed conversion requires eight nonnegative physical heats",
        ));
    }
    let total = sum(heat)?;
    let exported = Dd::from_f64(finite(step.load_work_j)?);
    let external = Dd::from_f64(finite(step.external_work_j)?) - exported;
    let battery_change = sum([
        step.battery.chemical_energy_change_j,
        step.battery.polarization_energy_change_j,
        step.battery.bus_energy_change_j,
    ])?;
    let battery_heat = sum([step.battery.series_heat_j, step.battery.polarization_heat_j])?;
    let converter_change = Dd::from_f64(finite(step.converter_storage_change_j)?);
    let converter_heat = sum(step.converter_conduction.joules)?;
    let battery_work = Dd::from_f64(finite(step.battery.terminal_work_j)?);
    let input_work = Dd::from_f64(finite(step.converter_terminal_work_j[0])?);
    let output_work = Dd::from_f64(finite(step.converter_terminal_work_j[1])?);
    let delta = Dd::from_f64(finite(step.native.delta_h)?);
    let residuals = [
        battery_change + battery_heat - battery_work,
        converter_change + converter_heat - input_work - output_work,
        battery_work + Dd::from_f64(finite(step.lead_heat_j)?) + input_work,
        output_work - external,
        battery_change + converter_change - delta,
        total - Dd::from_f64(finite(step.native.dissipated)?),
        Dd::from_f64(finite(step.native.supplied)?) - external,
        delta + total - external,
    ];
    if exported.to_f64() < 0.0
        || residuals
            .iter()
            .any(|r| !r.hi.is_finite() || !r.lo.is_finite() || r.abs().to_f64() > tolerance_j)
    {
        return Err(invalid(
            "battery, lead, converter and exported load fail independent energy balances",
        ));
    }
    Ok((heat, total, external))
}

impl BatteryFedConverterThermalNetwork<'_> {
    fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            links: self.links,
            resistors: &[],
        }
    }

    /// Bind each resistance to its original material state without changing any
    /// accepted charge, flux, OCV law or capacitance. Shared heat receivers sum.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        electrical: &BatteryFedBuckBoost,
        state: [f64; 6],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialBatteryFedBuckBoost, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        if state.iter().any(|v| !v.is_finite())
            || !electrical.system().hamiltonian(&state).is_finite()
        {
            return Err(invalid(
                "battery-fed conversion requires six finite original storage coordinates",
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
            all_resistances(
                electrical.parameters(),
                electrical.battery_parameters(),
                electrical.lead_resistance_ohm(),
            ),
            None,
            false,
            &cancel,
        )?);
        let (parameters, battery_parameters, lead) = parameters_at(
            electrical.parameters(),
            electrical.battery_parameters(),
            &channels,
        );
        native(
            electrical.battery_curve(),
            parameters,
            battery_parameters,
            electrical.modulation(),
            lead,
        )?;
        check(&cancel)?;
        Ok(MaterialBatteryFedBuckBoost {
            time_s: thermal.time_s(),
            curve: electrical.battery_curve().clone(),
            battery_parameters,
            parameters,
            modulation: electrical.modulation(),
            lead_resistance_ohm: lead,
            state,
            channels,
        })
    }

    /// Resolve all eight R(T,phase) sources in the original thermal residual.
    /// Commands are held over this interval; modulation may change without
    /// resetting storage. Output current is positive INTO the output bus.
    /// Load demand is exported work, and its nonlinear voltage response is
    /// solved at the native discrete effort in every thermal trial.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        converter: &MaterialBatteryFedBuckBoost,
        modulation: BridgeModulation,
        output_current_a: f64,
        load: Option<RegulatedPowerLoad>,
        electrical_controls: PowerLoadControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryFedConverterEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if converter.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || !output_current_a.is_finite()
        {
            return Err(invalid(
                "battery-fed conversion and heat need the same positive interval and finite input",
            ));
        }
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let fixed = all_resistances(
            converter.parameters,
            converter.battery_parameters,
            converter.lead_resistance_ohm,
        );
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
                "battery-fed conversion must retain all eight accepted material sources",
            ));
        }
        let cache = RefCell::new(None::<([u64; 8], [f64; 6], BatteryFedStep)>);
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
                let (parameters, battery_parameters, lead) = parameters_at(
                    converter.parameters,
                    converter.battery_parameters,
                    &channels,
                );
                let bits = all_resistances(parameters, battery_parameters, lead).map(f64::to_bits);
                let previous = cache
                    .borrow()
                    .as_ref()
                    .filter(|entry| entry.0 == bits)
                    .cloned();
                let reused = previous.is_some();
                let (next, step) = if let Some((_, next, step)) = previous {
                    (next, step)
                } else {
                    let mut original = native(
                        &converter.curve,
                        parameters,
                        battery_parameters,
                        modulation,
                        lead,
                    )?
                    .prepare_controlled(dt, electrical_controls, || (*cancel.borrow_mut())())
                    .map_err(electrical)?;
                    let mut next = [0.0; 6];
                    let step = original
                        .step_controlled(
                            &converter.state,
                            output_current_a,
                            load,
                            &mut next,
                            || (*cancel.borrow_mut())(),
                        )
                        .map_err(electrical)?;
                    (next, step)
                };
                let (heat, total, net_external_work_j) =
                    physical_account(&step, controls.total_energy_tolerance_j)?;
                check(&cancel)?;
                if !reused {
                    *cache.borrow_mut() = Some((bits, next, step.clone()));
                }
                Ok(DrivenCoefficients {
                    circuit: BatteryDriver {
                        converter: MaterialBatteryFedBuckBoost {
                            time_s: target_s,
                            curve: converter.curve.clone(),
                            battery_parameters,
                            parameters,
                            modulation,
                            lead_resistance_ohm: lead,
                            state: next,
                            channels: Arc::clone(&channels),
                        },
                        step,
                        heat,
                        total_heat_j: total.to_f64(),
                        net_external_work_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: BatteryFedConverterThermalSample {
                        network: coefficients.source,
                        channels,
                    },
                    circuit_trials: 1,
                })
            },
            |candidate, duration, checkpoint| {
                let mut by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                for (channel, heat) in candidate.converter.channels.iter().zip(candidate.heat) {
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
                            "battery-fed conversion heat cannot be represented as a rate",
                        ));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(BatteryFedConverterEnthalpyCandidate {
            thermal,
            converter,
            report,
        })
    }
}
