//! A finite SOC/OCV battery, physical leads and a power-loaded buck/boost stage.
//!
//! Both battery and converter buses take current and return voltage. They are
//! connected by the actual lead law I=(V_battery-V_input)/R_lead, not by treating
//! two voltage outputs as a lossless effort/flow pair. All six original states
//! evolve in the same native nonlinear solve, including load current limiting
//! and foldback. No sampled source voltage or lagged load current drives a
//! separate subsystem. Regeneration reverses the same physical current path.
//!
//! Converter modulation is frozen at preparation; output current and regulated
//! demand may change every step without changing stored charge or flux. The
//! source OCV/C and resistances are frozen over that interval. This averaged
//! lane adds no PWM ripple, semiconductor switching loss or battery protection.

use super::{
    AveragedBuckBoost, BridgeModulation, BuckBoostParameters, ConductionHeat, invalid, poll,
    positive,
    power_load::{PowerLoadControls, PowerLoadPoint, RegulatedPowerLoad, load_point},
};
use crate::circuit::battery::{BatteryInterval, BatteryOcvCurve, BatteryParameters, BatterySource};
use crate::{
    InterconnectionBudget, NonlinearPortPower, PhsError, PortHamiltonian, PreparedStepError,
    PreparedStepRecord, ResistanceWitness, StepWorkspace, interconnect_mapped,
    storage_energy_change,
};
use fs_math::{dd::Dd, det};

/// The same finite electrical source driving a reciprocal DC motor/generator.
pub mod motor;

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, v| s + Dd::from_f64(v))
        .to_f64()
}
fn product(a: f64, b: f64, c: f64) -> f64 {
    (Dd::from_f64(a) * Dd::from_f64(b) * Dd::from_f64(c)).to_f64()
}

/// An external terminal's current at the accepted discrete output voltage.
///
/// `I_into_bus = constant_a + slope_a_per_v * V_output`. Both coefficients are
/// held over one interval; positive current supplies the electrical subsystem.
/// This can represent an analytically condensed dynamic load. It is NOT a
/// resistor heat owner: its complete signed work belongs to the external
/// component, whose own storage/work balance must be checked by that caller.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AffineOutputCurrent {
    pub constant_a: f64,
    pub slope_a_per_v: f64,
}

impl AffineOutputCurrent {
    /// Evaluate the same current law used by the original implicit equation.
    pub fn current_into_bus_a(self, voltage_v: f64) -> Result<f64, PhsError> {
        let current = (Dd::from_f64(self.constant_a)
            + Dd::from_f64(self.slope_a_per_v) * Dd::from_f64(voltage_v))
        .to_f64();
        if [self.constant_a, self.slope_a_per_v, voltage_v, current]
            .into_iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid("affine output current is unrepresentable"));
        }
        Ok(current)
    }
}

fn admit_balances(
    controls: PowerLoadControls,
    balances: &[(f64, f64)],
    what: &'static str,
) -> Result<(), PhsError> {
    if balances.iter().any(|&(residual, scale)| {
        let allowance =
            controls.absolute_energy_tolerance_j + controls.relative_energy_tolerance * scale;
        !scale.is_finite()
            || !allowance.is_finite()
            || !residual.is_finite()
            || residual.abs() > allowance
    }) {
        return Err(invalid(what));
    }
    Ok(())
}

/// Physical accounts of the six electrical states on one accepted interval.
/// This may be part of a larger monolithic solve: it carries no independent
/// solver receipt or convergence claim for a separately stepped subsystem.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatteryFedInterval {
    pub battery: BatteryInterval,
    pub converter_storage_change_j: f64,
    pub converter_terminal_work_j: [f64; 2],
    pub converter_conduction: ConductionHeat,
    pub lead_current_a: f64,
    pub lead_heat_j: f64,
    /// Signed work entering the six-state electrical subsystem at its output.
    /// A connected motor makes this an internal drive exchange, not a new
    /// external energy source for the complete electromechanical assembly.
    pub external_work_j: f64,
    pub load: PowerLoadPoint,
    pub load_work_j: f64,
    pub curtailed_load_energy_j: f64,
    pub discrete_bus_voltage_v: [f64; 3],
}

impl BatteryFedInterval {
    /// Converter five channels, battery series/polarization, then physical lead.
    pub fn heat_joules(&self) -> [f64; 8] {
        let c = self.converter_conduction.joules;
        [
            c[0],
            c[1],
            c[2],
            c[3],
            c[4],
            self.battery.series_heat_j,
            self.battery.polarization_heat_j,
            self.lead_heat_j,
        ]
    }
    pub fn storage_change_j(&self) -> f64 {
        sum([
            self.battery.storage_change_j(),
            self.converter_storage_change_j,
        ])
    }
    pub fn balance_residual_j(&self) -> f64 {
        sum(self.heat_joules().into_iter().chain([
            self.battery.storage_change_j(),
            self.converter_storage_change_j,
            -self.external_work_j,
            self.load_work_j,
        ]))
    }
    fn scales(&self) -> [f64; 4] {
        let b = self.battery;
        [
            sum([
                b.chemical_energy_change_j.abs(),
                b.polarization_energy_change_j.abs(),
                b.bus_energy_change_j.abs(),
                b.terminal_work_j.abs(),
                b.heat_j(),
            ]),
            sum([
                self.converter_storage_change_j.abs(),
                sum(self.converter_conduction.joules),
                self.converter_terminal_work_j[0].abs(),
                self.converter_terminal_work_j[1].abs(),
            ]),
            sum([
                b.terminal_work_j.abs(),
                self.converter_terminal_work_j[0].abs(),
                self.lead_heat_j,
            ]),
            sum([self.external_work_j.abs(), self.load_work_j]),
        ]
    }
    fn admit(&self, controls: PowerLoadControls) -> Result<(), PhsError> {
        let heat = sum(self.heat_joules());
        let [source_scale, converter_scale, lead_scale, external_scale] = self.scales();
        let values = [
            self.converter_storage_change_j,
            self.converter_terminal_work_j[0],
            self.converter_terminal_work_j[1],
            self.lead_current_a,
            self.lead_heat_j,
            self.external_work_j,
            self.load_work_j,
            self.curtailed_load_energy_j,
            heat,
            self.storage_change_j(),
        ];
        if values
            .iter()
            .chain(&self.heat_joules())
            .chain(&self.discrete_bus_voltage_v)
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "battery-converter interval power is unrepresentable",
            ));
        }
        admit_balances(
            controls,
            &[
                (self.battery.balance_residual_j(), source_scale),
                (
                    sum([
                        self.converter_storage_change_j,
                        sum(self.converter_conduction.joules),
                        -self.converter_terminal_work_j[0],
                        -self.converter_terminal_work_j[1],
                    ]),
                    converter_scale,
                ),
                (
                    sum([
                        -self.battery.terminal_work_j,
                        -self.converter_terminal_work_j[0],
                        -self.lead_heat_j,
                    ]),
                    lead_scale,
                ),
                (
                    self.balance_residual_j(),
                    source_scale
                        .max(converter_scale)
                        .max(external_scale)
                        .max(heat),
                ),
            ],
            "battery, lead, converter and exported load work do not close",
        )
    }
}

/// Complete native solve, with independent accounts at each physical boundary.
/// Load energy leaves the electrical system; it is never converter or lead heat.
#[derive(Clone, Debug)]
pub struct BatteryFedStep {
    pub native: PreparedStepRecord,
    pub battery: BatteryInterval,
    pub converter_storage_change_j: f64,
    /// Signed work entering the converter's input/output buses, in that order.
    pub converter_terminal_work_j: [f64; 2],
    pub converter_conduction: ConductionHeat,
    /// Positive from battery bus to converter input, negative for regeneration.
    pub lead_current_a: f64,
    pub lead_heat_j: f64,
    /// External terminal's signed work entering the output bus.
    pub external_work_j: f64,
    pub load: PowerLoadPoint,
    pub load_work_j: f64,
    pub curtailed_load_energy_j: f64,
    /// Native discrete voltages: battery bus, converter input, converter output.
    /// These are work-conjugate discrete efforts, not exact time averages.
    pub discrete_bus_voltage_v: [f64; 3],
}

impl BatteryFedStep {
    /// Physical subaccount, without relabeling this solve's certificate as a
    /// certificate for a subsystem of another assembled trajectory.
    pub fn interval(&self) -> BatteryFedInterval {
        BatteryFedInterval {
            battery: self.battery,
            converter_storage_change_j: self.converter_storage_change_j,
            converter_terminal_work_j: self.converter_terminal_work_j,
            converter_conduction: self.converter_conduction,
            lead_current_a: self.lead_current_a,
            lead_heat_j: self.lead_heat_j,
            external_work_j: self.external_work_j,
            load: self.load,
            load_work_j: self.load_work_j,
            curtailed_load_energy_j: self.curtailed_load_energy_j,
            discrete_bus_voltage_v: self.discrete_bus_voltage_v,
        }
    }
    fn from_interval(native: PreparedStepRecord, interval: BatteryFedInterval) -> Self {
        Self {
            native,
            battery: interval.battery,
            converter_storage_change_j: interval.converter_storage_change_j,
            converter_terminal_work_j: interval.converter_terminal_work_j,
            converter_conduction: interval.converter_conduction,
            lead_current_a: interval.lead_current_a,
            lead_heat_j: interval.lead_heat_j,
            external_work_j: interval.external_work_j,
            load: interval.load,
            load_work_j: interval.load_work_j,
            curtailed_load_energy_j: interval.curtailed_load_energy_j,
            discrete_bus_voltage_v: interval.discrete_bus_voltage_v,
        }
    }
    /// Eight independent positive heat owners: the converter's original five
    /// channels, battery series resistance, battery polarization, then the lead.
    pub fn heat_joules(&self) -> [f64; 8] {
        self.interval().heat_joules()
    }

    pub fn balance_residual_j(&self) -> f64 {
        self.interval().balance_residual_j()
    }
}

/// State order: battery chemical/polarization/bus charges, then converter flux,
/// input charge and output charge. Positive lead resistance is physical data.
pub struct BatteryFedBuckBoost {
    battery: BatterySource,
    converter: AveragedBuckBoost,
    lead_resistance_ohm: f64,
    system: PortHamiltonian,
}

impl AveragedBuckBoost {
    /// Join the two current-input buses through a retained positive resistor.
    /// Original immutable storage and all physical loss witnesses are retained.
    /// Zero lead resistance needs a common-voltage constraint and is refused.
    pub fn with_battery_source(
        mut self,
        mut battery: BatterySource,
        lead_resistance_ohm: f64,
    ) -> Result<BatteryFedBuckBoost, PhsError> {
        if !positive(lead_resistance_ohm) || !positive(lead_resistance_ohm.recip()) {
            return Err(invalid(
                "battery-converter lead resistance must be positive and representable",
            ));
        }
        let direct_sum = interconnect_mapped(
            battery.share_system(),
            self.system.share_storage(),
            &[],
            &[],
            &[],
            InterconnectionBudget {
                max_states: 6,
                max_total_ports: 3,
                max_dense_bytes: 1 << 18,
                max_coupling_multiply_terms: 128,
            },
        )?;
        let g = lead_resistance_ohm.recip();
        let root = det::sqrt(g);
        let lead =
            ResistanceWitness::try_from_source_rows(2, vec![g, -g, -g, g], vec![root, -root])?;
        let system = direct_sum.with_port_conductance_witness(&[0, 1], lead)?;
        Ok(BatteryFedBuckBoost {
            battery,
            converter: self,
            lead_resistance_ohm,
            system,
        })
    }
}

impl BatteryFedBuckBoost {
    /// Native additive current ports are battery bus, converter input, output.
    /// The prepared finite-source path holds the first two at zero; its only
    /// external electrical source/load boundary is the converter output bus.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    pub fn into_system(self) -> PortHamiltonian {
        self.system
    }
    pub fn battery_parameters(&self) -> BatteryParameters {
        self.battery.parameters()
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        self.battery.curve()
    }
    pub fn parameters(&self) -> BuckBoostParameters {
        self.converter.parameters()
    }
    pub fn modulation(&self) -> BridgeModulation {
        self.converter.modulation()
    }
    pub fn lead_resistance_ohm(&self) -> f64 {
        self.lead_resistance_ohm
    }

    /// Retain both bus capacitors separately, including physical inrush if their
    /// initial voltages differ. No equal-voltage projection removes their energy.
    pub fn state(
        &self,
        soc: f64,
        polarization_v: f64,
        battery_bus_v: f64,
        converter_current_a: f64,
        converter_bus_voltage_v: [f64; 2],
    ) -> Result<[f64; 6], PhsError> {
        let b = self.battery.state(soc, polarization_v, battery_bus_v)?;
        let c = self
            .converter
            .state(converter_current_a, converter_bus_voltage_v)?;
        let state = [b[0], b[1], b[2], c[0], c[1], c[2]];
        if !self.system.hamiltonian(&state).is_finite() {
            return Err(invalid(
                "battery-converter total storage is unrepresentable",
            ));
        }
        Ok(state)
    }

    // Shared by the original six-state solve and larger monolithic drives.
    // Efforts are the accepted native discrete-gradient prefix, not a sampled
    // voltage command or an independently propagated component trajectory.
    pub(super) fn account_interval(
        &self,
        initial: &[f64; 6],
        next: &[f64; 6],
        effort: &[f64; 6],
        output_current_a: f64,
        load: Option<RegulatedPowerLoad>,
        duration_s: f64,
    ) -> Result<BatteryFedInterval, PhsError> {
        let p = self.parameters();
        load_point(load, initial[5] / p.bus_capacitance_f[1])?;
        load_point(load, next[5] / p.bus_capacitance_f[1])?;
        let point = load_point(load, effort[5])?;
        let lead_current_a = ((Dd::from_f64(effort[2]) - Dd::from_f64(effort[4]))
            / Dd::from_f64(self.lead_resistance_ohm))
        .to_f64();
        let battery = self.battery.interval(
            &[initial[0], initial[1], initial[2]],
            &[next[0], next[1], next[2]],
            -lead_current_a,
            duration_s,
        )?;
        let converter_storage_change_j = storage_energy_change(
            self.converter.system().storage(),
            &initial[3..],
            &next[3..],
            &mut [0.0; 3],
            &mut [0.0; 3],
        );
        let joules = self.converter.effective_resistance_ohm.map(|r| {
            (Dd::from_f64(r)
                * Dd::from_f64(effort[3])
                * Dd::from_f64(effort[3])
                * Dd::from_f64(duration_s))
            .to_f64()
        });
        let lead_heat_j = (Dd::from_f64(self.lead_resistance_ohm)
            * Dd::from_f64(lead_current_a)
            * Dd::from_f64(lead_current_a)
            * Dd::from_f64(duration_s))
        .to_f64();
        let external_work_j = product(duration_s, effort[5], output_current_a);
        let load_work_j = duration_s * point.delivered_power_w;
        Ok(BatteryFedInterval {
            battery,
            converter_storage_change_j,
            converter_terminal_work_j: [
                product(duration_s, effort[4], lead_current_a),
                sum([external_work_j, -load_work_j]),
            ],
            converter_conduction: ConductionHeat { joules },
            lead_current_a,
            lead_heat_j,
            external_work_j,
            load: point,
            load_work_j,
            curtailed_load_energy_j: duration_s * point.curtailed_power_w,
            discrete_bus_voltage_v: [effort[2], effort[4], effort[5]],
        })
    }

    pub fn prepare_controlled(
        self,
        duration_s: f64,
        controls: PowerLoadControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedBatteryFedBuckBoost, PreparedStepError> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(invalid(
                "battery-fed converter requires positive duration and finite energy allowances",
            )
            .into());
        }
        let mut workspace = StepWorkspace::new(&self.system)?;
        workspace.set_solve_budget(controls.solve)?;
        workspace.set_backtracking_limit(controls.backtracking_halvings)?;
        poll(&mut cancelled)?;
        Ok(PreparedBatteryFedBuckBoost {
            model: self,
            workspace,
            duration_s,
            controls,
        })
    }
}

/// Fixed six-state native Newton scratch. Source/load refusal or cancellation
/// leaves the caller's output unchanged; successful warm steps allocate nothing.
pub struct PreparedBatteryFedBuckBoost {
    model: BatteryFedBuckBoost,
    workspace: StepWorkspace,
    duration_s: f64,
    controls: PowerLoadControls,
}

impl PreparedBatteryFedBuckBoost {
    pub fn system(&self) -> &PortHamiltonian {
        self.model.system()
    }
    pub fn parameters(&self) -> BuckBoostParameters {
        self.model.parameters()
    }
    pub fn battery_parameters(&self) -> BatteryParameters {
        self.model.battery_parameters()
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        self.model.battery_curve()
    }
    pub fn modulation(&self) -> BridgeModulation {
        self.model.modulation()
    }
    pub fn lead_resistance_ohm(&self) -> f64 {
        self.model.lead_resistance_ohm()
    }
    pub fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Positive output current supplies energy, negative current draws it.
    /// Optional regulated demand is solved at the same discrete output voltage.
    /// Voltage support is checked at initial/final and discrete states; this is
    /// not a certificate against unresolved between-endpoint transients.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 6],
        output_current_a: f64,
        load: Option<RegulatedPowerLoad>,
        next: &mut [f64; 6],
        cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryFedStep, PreparedStepError> {
        self.step_affine_output_current_controlled(
            initial,
            AffineOutputCurrent {
                constant_a: output_current_a,
                slope_a_per_v: 0.0,
            },
            load,
            next,
            cancelled,
        )
    }

    /// Solve an affine external terminal and optional demand at the SAME native
    /// discrete output voltage. Its exact analytic tangent joins the existing
    /// Newton equation; no lagged current, extra integrator or energy correction
    /// is used. All canonical battery, lead, converter and load accounts remain
    /// unchanged. A failed solve or account leaves `next` unchanged.
    ///
    /// An active coefficient is permitted because this is an externally
    /// supplied port. The caller must account for that external component's
    /// energy and admit its operating limits before publishing either state.
    pub fn step_affine_output_current_controlled(
        &mut self,
        initial: &[f64; 6],
        terminal: AffineOutputCurrent,
        load: Option<RegulatedPowerLoad>,
        next: &mut [f64; 6],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryFedStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let p = self.model.parameters();
        load_point(load, initial[5] / p.bus_capacitance_f[1])?;
        terminal.current_into_bus_a(initial[5] / p.bus_capacitance_f[1])?;
        let port = |_: &[f64], effort: &[f64], out: &mut [f64]| {
            if effort.len() != 6 || out.len() != 6 {
                return None;
            }
            let point = load_point(load, effort[5]).ok()?;
            out.fill(0.0);
            let dependent_current = terminal.slope_a_per_v * effort[5];
            out[5] = sum([point.current_a, -dependent_current]);
            Some(NonlinearPortPower {
                supplied_power_w: sum([
                    -point.delivered_power_w,
                    product(terminal.slope_a_per_v, effort[5], effort[5]),
                ]),
                dissipated_power_w: 0.0,
            })
        };
        let tangent = |_: &[f64], effort: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            if effort.len() != 6 || de.len() != 6 || out.len() != 6 {
                return false;
            }
            let Ok(point) = load_point(load, effort[5]) else {
                return false;
            };
            out.fill(0.0);
            out[5] = (point.incremental_conductance_s - terminal.slope_a_per_v) * de[5];
            out[5].is_finite()
        };
        let (mut candidate, mut voltage) = ([0.0; 6], [0.0; 3]);
        let dt = self.duration_s;
        let native = self.workspace.step_into_supplied_controlled(
            &self.model.system,
            initial,
            &[0.0, 0.0, terminal.constant_a],
            dt,
            &mut candidate,
            &mut voltage,
            &port,
            None,
            Some(&tangent),
            &mut cancelled,
        )?;
        let e = self
            .workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("battery-fed converter requires accepted native efforts"))?;
        let effort = core::array::from_fn(|i| e[i]);
        let output_current_a = terminal.current_into_bus_a(effort[5])?;
        let interval = self.model.account_interval(
            initial,
            &candidate,
            &effort,
            output_current_a,
            load,
            dt,
        )?;
        interval.admit(self.controls)?;
        let heat = sum(interval.heat_joules());
        let [source_scale, converter_scale, _, external_scale] = interval.scales();
        admit_balances(
            self.controls,
            &[
                (
                    sum([interval.storage_change_j(), -native.delta_h]),
                    source_scale.max(converter_scale),
                ),
                (sum([heat, -native.dissipated]), heat),
                (
                    sum([
                        native.supplied,
                        -interval.external_work_j,
                        interval.load_work_j,
                    ]),
                    external_scale,
                ),
            ],
            "battery-converter component accounts disagree with the native solve",
        )?;
        poll(&mut cancelled)?;
        *next = candidate;
        Ok(BatteryFedStep::from_interval(native, interval))
    }
}
