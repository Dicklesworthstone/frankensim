//! Constant or original voltage-dependent device charge in the finite battery/motor chain.
//!
//! The eleven-state lanes retain battery3, converter5, motor3. The gate-driven
//! fifteen-state lane retains battery3, converter9, motor3: converter flux,
//! bus/bus/GROUPED-node/GROUPED-node charges, then four gate charges at indices
//! 8..12; motor flux and rotor/stator momenta are at12,13,14. Grouped node
//! charge includes that leg's high-side gate electrode. Initialize converter
//! coordinates through the retained parameters with authored physical voltages;
//! old bus/node charges cannot be copied into a different capacitor topology.
//!
//! Every gate command shares one original native workspace, finite battery
//! storage and reciprocal motor law. Actual gate-source voltage owns driven
//! conductance; gate supply work and four driver resistor heats are explicit.
//! Intrinsic drain, gate/source and Miller energies are reversible storage.
//! No extra Coss loss map, node reset, sampled supply, split solve, inferred
//! efficiency, carrier recovery or finite gate-rail energy is introduced.

use super::*;
use crate::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use crate::circuit::capacitance::CapacitorChargeCurve;
use crate::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use crate::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use crate::circuit::motor::{DcMotor, DcMotorParameters};
use crate::{InterconnectionBudget, ResistanceWitness, interconnect_mapped};
use fs_math::det;

#[derive(Clone, Debug)]
pub struct CapacitiveBatteryMotorCycleStep {
    pub native: PreparedStepRecord,
    pub converter: CapacitiveConverterAccount,
    /// Chemical, polarization, battery bus, COMPLETE five- or nine-state converter,
    /// motor magnetic, rotor kinetic, stator kinetic energy changes.
    pub storage_change_j: [f64; 7],
    pub battery_terminal_work_j: f64,
    pub motor_electrical_work_j: f64,
    pub electromagnetic_work_j: f64,
    pub shaft_work_j: [f64; 2],
    pub load_work_j: f64,
    pub curtailed_load_energy_j: f64,
    pub mean_bus_voltage_v: [f64; 3],
    pub mean_lead_current_a: f64,
    pub mean_motor_current_a: f64,
    pub mean_shaft_speed_rad_s: [f64; 2],
    pub last_auxiliary_load: PowerLoadPoint,
    pub angular_momentum_residual_nms: f64,
    pub integration_steps: usize,
    /// Original series and polarization resistors, respectively.
    pub battery_heat_j: [f64; 2],
    pub lead_heat_j: f64,
    pub copper_heat_j: f64,
    pub bearing_heat_j: f64,
}
impl CapacitiveBatteryMotorCycleStep {
    /// Converter winding/four channels, battery series/polarization, lead,
    /// motor copper, bearing, then four body diodes. Device C is storage only.
    pub fn heat_joules(&self) -> [f64; 14] {
        let c = self.converter.channel_heat_j;
        let d = self.converter.diode_heat_j;
        [
            self.converter.winding_heat_j,
            c[0],
            c[1],
            c[2],
            c[3],
            self.battery_heat_j[0],
            self.battery_heat_j[1],
            self.lead_heat_j,
            self.copper_heat_j,
            self.bearing_heat_j,
            d[0],
            d[1],
            d[2],
            d[3],
        ]
    }
    /// Complete original fourteen heats plus four separately owned gate drivers.
    pub fn all_heat_joules(&self) -> [f64; 18] {
        let mut out = [0.0; 18];
        out[..14].copy_from_slice(&self.heat_joules());
        if let Some(gate) = &self.converter.gate {
            out[14..].copy_from_slice(&gate.resistor_heat_j);
        }
        out
    }
    pub fn total_storage_change_j(&self) -> f64 {
        total(self.storage_change_j)
    }
    pub fn balance_residual_j(&self) -> f64 {
        total(
            self.storage_change_j
                .into_iter()
                .chain(self.all_heat_joules())
                .chain([
                    self.load_work_j,
                    -self.shaft_work_j[0],
                    -self.shaft_work_j[1],
                ])
                .chain(self.converter.gate_source_work_j().map(|w| -w)),
        )
    }
    fn work(&self) -> [f64; 12] {
        let g = self.converter.gate_source_work_j();
        [
            self.battery_terminal_work_j,
            self.converter.terminal_work_j[0],
            self.converter.terminal_work_j[1],
            self.motor_electrical_work_j,
            self.electromagnetic_work_j,
            self.shaft_work_j[0],
            self.shaft_work_j[1],
            self.load_work_j,
            g[0],
            g[1],
            g[2],
            g[3],
        ]
    }
    fn means(&self) -> [f64; 7] {
        [
            self.mean_bus_voltage_v[0],
            self.mean_bus_voltage_v[1],
            self.mean_bus_voltage_v[2],
            self.mean_lead_current_a,
            self.mean_motor_current_a,
            self.mean_shaft_speed_rad_s[0],
            self.mean_shaft_speed_rad_s[1],
        ]
    }
    fn admit(
        &self,
        changes: [f64; 7],
        work: [f64; 12],
        absolute: f64,
        relative: f64,
    ) -> Result<(), PhsError> {
        let s = self.storage_change_j;
        let w = self.work();
        let h = self.all_heat_joules();
        let converter_heat = total(self.converter.all_heat_joules());
        let heat = total(h);
        let source_scale = total(changes[..3].iter().copied().chain([work[0], h[5], h[6]]));
        let converter_scale = total(
            [changes[3], converter_heat, work[1], work[2]]
                .into_iter()
                .chain(work[8..].iter().copied()),
        );
        let lead_scale = total([work[0], work[1], h[7]]);
        let output_scale = total([work[2], work[3], work[7]]);
        let magnetic_scale = total([changes[4], h[8], work[3], work[4]]);
        let mechanical_scale = total([changes[5], changes[6], h[9], work[4], work[5], work[6]]);
        let scale = source_scale
            .max(converter_scale)
            .max(lead_scale)
            .max(output_scale)
            .max(magnetic_scale)
            .max(mechanical_scale)
            .max(heat);
        if s.iter()
            .chain(&w)
            .chain(&self.means())
            .any(|v| !v.is_finite())
            || h.iter()
                .chain(&[self.load_work_j, self.curtailed_load_energy_j])
                .any(|v| !v.is_finite() || *v < 0.0)
            || !self.angular_momentum_residual_nms.is_finite()
            || !self.native.solver_residual.is_finite()
        {
            return Err(invalid(
                "capacitive finite-source drive work, heat or effort is unrepresentable",
            ));
        }
        for (error, magnitude) in [
            (total([s[0], s[1], s[2], h[5], h[6], -w[0]]), source_scale),
            (
                total(
                    [s[3], converter_heat, -w[1], -w[2]]
                        .into_iter()
                        .chain(w[8..].iter().map(|v| -v)),
                ),
                converter_scale,
            ),
            (total([-w[0], -w[1], -h[7]]), lead_scale),
            (total([w[2], w[3], w[7]]), output_scale),
            (total([s[4], h[8], w[4], -w[3]]), magnetic_scale),
            (
                total([s[5], s[6], h[9], -w[4], -w[5], -w[6]]),
                mechanical_scale,
            ),
            (
                total([s[3], -self.converter.storage_change_j]),
                converter_scale,
            ),
            (total(s.into_iter().chain([-self.native.delta_h])), scale),
            (total([heat, -self.native.dissipated]), heat),
            (
                total(
                    [self.native.supplied, -w[5], -w[6], w[7]]
                        .into_iter()
                        .chain(w[8..].iter().map(|v| -v)),
                ),
                total(
                    [work[5], work[6], work[7]]
                        .into_iter()
                        .chain(work[8..].iter().copied()),
                ),
            ),
            (self.balance_residual_j(), scale),
        ] {
            admit(
                error,
                magnitude,
                absolute,
                relative,
                "battery, capacitive converter, lead, motor and shaft work do not close on one native trajectory",
            )?;
        }
        Ok(())
    }
}
pub struct PreparedCapacitiveBatteryMotorCycle<
    P = CapacitiveCommutationParameters,
    const C: usize = 5,
    const N: usize = 11,
> {
    parameters: P,
    battery: BatterySource,
    core: ConverterCore<C>,
    lead_resistance_ohm: f64,
    motor_parameters: DcMotorParameters,
    system: PortHamiltonian,
    workspace: StepWorkspace,
    plan: GatePlan,
}
fn core_state<const C: usize, const N: usize>(state: &[f64; N]) -> [f64; C] {
    core::array::from_fn(|i| state[3 + i])
}

/// Original device charge curves in the SAME eleven-state native drive.
pub type PreparedChargeBatteryMotorCycle =
    PreparedCapacitiveBatteryMotorCycle<charge::ChargeCommutationParameters>;

fn converter_parameters(
    model: &BatteryMotorDrive,
    devices: [SwitchDevice; 4],
) -> CommutationParameters {
    let p = model.electrical().parameters();
    CommutationParameters {
        inductance_h: p.inductance_h,
        bus_capacitance_f: p.bus_capacitance_f,
        winding_resistance_ohm: p.winding_resistance_ohm,
        devices,
    }
}
impl PreparedCapacitiveBatteryMotorCycle {
    /// The original model supplies source/lead, converter L/Cb/winding and motor
    /// laws. Explicit devices replace averaged bridge conduction. Four authored
    /// CONSTANT device capacitances add physical node charge/storage; no old
    /// stored state is silently promoted or reset.
    #[allow(clippy::too_many_arguments)]
    pub fn new_controlled(
        model: BatteryMotorDrive,
        devices: [SwitchDevice; 4],
        output_capacitance_f: [f64; 4],
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let parameters = CapacitiveCommutationParameters {
            converter: converter_parameters(&model, devices),
            output_capacitance_f,
        };
        let core = ConverterCore::new(parameters)?;
        Self::prepare(
            model,
            parameters,
            core,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            cancelled,
        )
    }
}
impl PreparedCapacitiveBatteryMotorCycle<charge::ChargeCommutationParameters> {
    /// Four complete original V(q) laws retain reversible device energy and
    /// their finite support. Explicit total-electrode initialization is required;
    /// no sampled differential capacitance or switching-loss map is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn new_charge_controlled(
        model: BatteryMotorDrive,
        devices: [SwitchDevice; 4],
        output_charge_curves: [CapacitorChargeCurve; 4],
        maximum_cells_per_chord: usize,
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        let parameters = charge::ChargeCommutationParameters {
            converter: converter_parameters(&model, devices),
            output_charge_curves,
            maximum_cells_per_chord,
        };
        let core = ConverterCore::new_charge(&parameters)?;
        Self::prepare(
            model,
            parameters,
            core,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            cancelled,
        )
    }
}
/// Retained gate/Miller charge in the same finite battery and two-shaft drive.
/// The default gate-only conductance and complete signed I(Vds,Vgs) channel
/// responses use the same fifteen-state workspace and original heat accounts.
pub type PreparedGateBatteryMotorCycle<C = GateConductanceCurve> =
    PreparedCapacitiveBatteryMotorCycle<gate::GateCommutationParameters<C>, 9, 15>;
impl<C: GateChannelResponse + Clone + 'static>
    PreparedCapacitiveBatteryMotorCycle<gate::GateCommutationParameters<C>, 9, 15>
{
    pub fn new_gate_controlled(
        model: BatteryMotorDrive,
        parameters: gate::GateCommutationParameters<C>,
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        if converter_parameters(&model, parameters.converter.devices) != parameters.converter {
            return Err(invalid(
                "gate drive must retain original converter L, Cbus and winding resistance",
            )
            .into());
        }
        let core = ConverterCore::<9>::new_gate(&parameters)?;
        Self::prepare(
            model,
            parameters,
            core,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            cancelled,
        )
    }
}
impl<P: Clone, const C: usize, const N: usize> PreparedCapacitiveBatteryMotorCycle<P, C, N> {
    pub fn parameters(&self) -> P {
        self.parameters.clone()
    }
}
impl<P, const C: usize, const N: usize> PreparedCapacitiveBatteryMotorCycle<P, C, N> {
    #[allow(clippy::too_many_arguments)]
    fn prepare(
        model: BatteryMotorDrive,
        parameters: P,
        mut core: ConverterCore<C>,
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        if C.checked_add(6) != Some(N) {
            return Err(invalid("complete gate/converter drive dimension mismatch").into());
        }
        let plan = GatePlan::new(
            core.parameters.devices,
            intervals,
            substeps_per_interval,
            electrical,
            controls,
            core.gate_parameters().map(|p| &p.gate_drivers),
            &mut cancelled,
        )?;
        core.admit_work(&plan)?;
        let mut battery = BatterySource::new(
            model.electrical().battery_curve().clone(),
            model.electrical().battery_parameters(),
        )?;
        let lead_resistance_ohm = model.electrical().lead_resistance_ohm();
        let motor_parameters = model.motor_parameters();
        let direct = interconnect_mapped(
            battery.share_system(),
            core.system.share_storage(),
            &[],
            &[],
            &[],
            InterconnectionBudget {
                max_states: C + 3,
                max_total_ports: 3,
                max_dense_bytes: 1 << 19,
                max_coupling_multiply_terms: 256,
            },
        )?;
        let g = 1.0 / lead_resistance_ohm;
        let root = det::sqrt(g);
        let lead =
            ResistanceWitness::try_from_source_rows(2, vec![g, -g, -g, g], vec![root, -root])?;
        let electrical_system = direct.with_port_conductance_witness(&[0, 1], lead)?;
        let system = interconnect_mapped(
            electrical_system,
            DcMotor::new(motor_parameters)?.into_system(),
            &[2],
            &[0],
            &[1.0],
            InterconnectionBudget {
                max_states: N,
                max_total_ports: 6,
                max_dense_bytes: 1 << 20,
                max_coupling_multiply_terms: 512,
            },
        )?;
        let workspace = plan.workspace(&system)?;
        poll(&mut cancelled)?;
        Ok(Self {
            parameters,
            battery,
            core,
            lead_resistance_ohm,
            motor_parameters,
            system,
            workspace,
            plan,
        })
    }
    /// Explicit concatenation of original battery3, complete converter5/9,
    /// and original motor3. The gate lane includes its four gate charges before
    /// motor state. Initialize through the matching parameters().state method.
    pub fn state(
        &self,
        battery_state: &[f64; 3],
        converter_state: &[f64; C],
        motor_state: &[f64; 3],
    ) -> Result<[f64; N], PhsError> {
        self.core.energy(converter_state)?;
        let mut state = [0.0; N];
        state[..3].copy_from_slice(battery_state);
        state[3..3 + C].copy_from_slice(converter_state);
        state[3 + C..].copy_from_slice(motor_state);
        if state.iter().any(|v| !v.is_finite()) || !self.system.hamiltonian(&state).is_finite() {
            return Err(invalid(
                "finite supported original battery/converter/gate/motor coordinates required",
            ));
        }
        Ok(state)
    }
    pub fn battery_parameters(&self) -> BatteryParameters {
        self.battery.parameters()
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        self.battery.curve()
    }
    pub fn lead_resistance_ohm(&self) -> f64 {
        self.lead_resistance_ohm
    }
    pub fn motor_parameters(&self) -> DcMotorParameters {
        self.motor_parameters
    }
    pub fn intervals(&self) -> &[CommutationInterval] {
        &self.plan.intervals
    }
    pub fn period_s(&self) -> f64 {
        self.plan.period_s
    }
    pub fn integration_steps(&self) -> usize {
        self.plan.integration_steps
    }
    pub fn substeps_per_interval(&self) -> usize {
        self.plan.substeps
    }
    pub fn cached_interval_count(&self) -> usize {
        self.plan.caches.len()
    }
    pub fn step_controlled(
        &mut self,
        initial: &[f64; N],
        command: BatteryMotorSwitchCommand,
        next: &mut [f64; N],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveBatteryMotorCycleStep, PreparedStepError> {
        self.advance(initial, |_| command, next, cancelled)
    }
    pub fn step_interval_commands_controlled(
        &mut self,
        initial: &[f64; N],
        commands: &[BatteryMotorSwitchCommand],
        next: &mut [f64; N],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveBatteryMotorCycleStep, PreparedStepError> {
        if commands.len() != self.plan.intervals.len() {
            return Err(invalid(
                "one torque/auxiliary command required per original capacitive interval",
            )
            .into());
        }
        self.advance(initial, |k| commands[k], next, cancelled)
    }
    #[allow(clippy::too_many_arguments)]
    fn step_native(
        battery: &BatterySource,
        core: &ConverterCore<C>,
        lead: f64,
        m: DcMotorParameters,
        system: &PortHamiltonian,
        workspace: &mut StepWorkspace,
        gate: &GateStep,
        controls: PowerLoadControls,
        initial: &[f64; N],
        command: BatteryMotorSwitchCommand,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<([f64; N], CapacitiveBatteryMotorCycleStep), PreparedStepError> {
        poll(cancel)?;
        let initial_core = core_state(initial);
        let start = core.effort(&initial_core)?;
        core.point(gate, &start)?;
        load_point(command.auxiliary_load, start[2])?;
        let port = |_: &[f64], e: &[f64], out: &mut [f64]| {
            if e.len() != N || out.len() != N {
                return None;
            }
            let point = core.point(gate, &e[3..3 + C]).ok()?;
            let load = load_point(command.auxiliary_load, e[5]).ok()?;
            out.fill(0.0);
            let heat = point.flow(&mut out[3..3 + C])?;
            out[5] = number(point.legs[1].high_current + Dd::from_f64(load.current_a)).ok()?;
            Some(NonlinearPortPower {
                supplied_power_w: number(
                    point.supplied_power() - Dd::from_f64(load.delivered_power_w),
                )
                .ok()?,
                dissipated_power_w: heat,
            })
        };
        let tangent = |_: &[f64], e: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            if e.len() != N || de.len() != N || out.len() != N {
                return false;
            }
            let Ok(point) = core.point(gate, &e[3..3 + C]) else {
                return false;
            };
            let Ok(load) = load_point(command.auxiliary_load, e[5]) else {
                return false;
            };
            out.fill(0.0);
            if !point.tangent(&de[3..3 + C], &mut out[3..3 + C]) {
                return false;
            }
            out[5] = total([out[5], load.incremental_conductance_s * de[5]]);
            out.iter().all(|v| v.is_finite())
        };
        let dt = gate.duration_s;
        let (mut next, mut outputs) = ([0.0; N], [0.0; 4]);
        let native = workspace.step_into_supplied_controlled(
            system,
            initial,
            &[
                0.0,
                0.0,
                command.shaft_torque_nm[0],
                command.shaft_torque_nm[1],
            ],
            dt,
            &mut next,
            &mut outputs,
            &port,
            None,
            Some(&tangent),
            &mut *cancel,
        )?;
        let e = workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("capacitive drive requires its original accepted effort"))?;
        let ce: &[f64; C] = (&e[3..3 + C])
            .try_into()
            .map_err(|_| invalid("capacitive drive converter effort dimension"))?;
        let nodes = core.point(gate, ce)?;
        let load = load_point(command.auxiliary_load, e[5])?;
        let next_core = core_state(&next);
        let end = core.effort(&next_core)?;
        core.point(gate, &end)?;
        load_point(command.auxiliary_load, end[2])?;
        let lead_current = number((Dd::from_f64(e[2]) - Dd::from_f64(e[4])) / Dd::from_f64(lead))?;
        let b = battery.interval(
            &[initial[0], initial[1], initial[2]],
            &[next[0], next[1], next[2]],
            -lead_current,
            dt,
        )?;
        let output = number(-Dd::from_f64(e[3 + C]) - Dd::from_f64(load.current_a))?;
        let converter = core.account(
            &initial_core,
            &next_core,
            ce,
            nodes,
            [lead_current, output],
            dt,
            controls,
        )?;
        let relative_speed = total([e[4 + C], -e[5 + C]]);
        let step = CapacitiveBatteryMotorCycleStep {
            native,
            storage_change_j: [
                b.chemical_energy_change_j,
                b.polarization_energy_change_j,
                b.bus_energy_change_j,
                converter.storage_change_j,
                quadratic(initial[3 + C], next[3 + C], m.inductance_h),
                quadratic(initial[4 + C], next[4 + C], m.rotor_inertia_kg_m2),
                quadratic(initial[5 + C], next[5 + C], m.stator_inertia_kg_m2),
            ],
            converter,
            battery_terminal_work_j: b.terminal_work_j,
            motor_electrical_work_j: product(dt, e[5], e[3 + C]),
            electromagnetic_work_j: (Dd::from_f64(dt)
                * Dd::from_f64(m.coupling_nm_per_a)
                * Dd::from_f64(e[3 + C])
                * Dd::from_f64(relative_speed))
            .to_f64(),
            shaft_work_j: [
                product(dt, command.shaft_torque_nm[0], e[4 + C]),
                product(dt, command.shaft_torque_nm[1], e[5 + C]),
            ],
            load_work_j: (Dd::from_f64(dt) * Dd::from_f64(load.delivered_power_w)).to_f64(),
            curtailed_load_energy_j: (Dd::from_f64(dt) * Dd::from_f64(load.curtailed_power_w))
                .to_f64(),
            mean_bus_voltage_v: [e[2], e[4], e[5]],
            mean_lead_current_a: lead_current,
            mean_motor_current_a: e[3 + C],
            mean_shaft_speed_rad_s: [e[4 + C], e[5 + C]],
            last_auxiliary_load: load,
            angular_momentum_residual_nms: (Dd::from_f64(next[4 + C])
                - Dd::from_f64(initial[4 + C])
                + Dd::from_f64(next[5 + C])
                - Dd::from_f64(initial[5 + C])
                - Dd::from_f64(dt)
                    * (Dd::from_f64(command.shaft_torque_nm[0])
                        + Dd::from_f64(command.shaft_torque_nm[1])))
            .to_f64(),
            integration_steps: 1,
            battery_heat_j: [b.series_heat_j, b.polarization_heat_j],
            lead_heat_j: (Dd::from_f64(dt)
                * Dd::from_f64(lead)
                * Dd::from_f64(lead_current)
                * Dd::from_f64(lead_current))
            .to_f64(),
            copper_heat_j: (Dd::from_f64(dt)
                * Dd::from_f64(m.winding_resistance_ohm)
                * Dd::from_f64(e[3 + C])
                * Dd::from_f64(e[3 + C]))
            .to_f64(),
            bearing_heat_j: (Dd::from_f64(dt)
                * Dd::from_f64(m.bearing_drag_nm_s_per_rad)
                * Dd::from_f64(relative_speed)
                * Dd::from_f64(relative_speed))
            .to_f64(),
        };
        step.admit(
            step.storage_change_j.map(f64::abs),
            step.work().map(f64::abs),
            controls.absolute_energy_tolerance_j,
            controls.relative_energy_tolerance,
        )?;
        poll(cancel)?;
        Ok((next, step))
    }
    fn advance(
        &mut self,
        initial: &[f64; N],
        command: impl Fn(usize) -> BatteryMotorSwitchCommand,
        next: &mut [f64; N],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveBatteryMotorCycleStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let mut converter = ConverterAccumulator::new();
        let (mut changes, mut absolute_changes) = ([Dd::ZERO; 7], [Dd::ZERO; 7]);
        let (mut work, mut absolute_work) = ([Dd::ZERO; 12], [Dd::ZERO; 12]);
        let (mut means, mut heat, mut native) = ([Dd::ZERO; 7], [Dd::ZERO; 5], [Dd::ZERO; 3]);
        let (mut curtailed, mut impulse) = (Dd::ZERO, Dd::ZERO);
        let (mut iterations, mut residual) = (0usize, 0.0_f64);
        let mut last = None;
        for k in 0..self.plan.intervals.len() {
            poll(&mut cancelled)?;
            let command = command(k);
            for j in 0..self.plan.substeps {
                poll(&mut cancelled)?;
                let gate = self.plan.step(k, j);
                let (endpoint, step) = Self::step_native(
                    &self.battery,
                    &self.core,
                    self.lead_resistance_ohm,
                    self.motor_parameters,
                    &self.system,
                    &mut self.workspace,
                    gate,
                    self.plan.electrical,
                    &state,
                    command,
                    &mut cancelled,
                )?;
                let dt = Dd::from_f64(gate.duration_s);
                converter.push(&step.converter, gate.duration_s);
                for i in 0..7 {
                    changes[i] = changes[i] + Dd::from_f64(step.storage_change_j[i]);
                    absolute_changes[i] =
                        absolute_changes[i] + Dd::from_f64(step.storage_change_j[i].abs());
                }
                for (i, v) in step.work().into_iter().enumerate() {
                    work[i] = work[i] + Dd::from_f64(v);
                    absolute_work[i] = absolute_work[i] + Dd::from_f64(v.abs());
                }
                let sample = step.means();
                for i in 0..7 {
                    means[i] = means[i] + dt * Dd::from_f64(sample[i]);
                }
                for (s, v) in heat.iter_mut().zip([
                    step.battery_heat_j[0],
                    step.battery_heat_j[1],
                    step.lead_heat_j,
                    step.copper_heat_j,
                    step.bearing_heat_j,
                ]) {
                    *s = *s + Dd::from_f64(v);
                }
                for (s, v) in native.iter_mut().zip([
                    step.native.delta_h,
                    step.native.dissipated,
                    step.native.supplied,
                ]) {
                    *s = *s + Dd::from_f64(v);
                }
                curtailed = curtailed + Dd::from_f64(step.curtailed_load_energy_j);
                impulse = impulse
                    + dt * (Dd::from_f64(command.shaft_torque_nm[0])
                        + Dd::from_f64(command.shaft_torque_nm[1]));
                iterations = iterations
                    .checked_add(step.native.newton_iters)
                    .ok_or_else(|| invalid("capacitive drive Newton count overflow"))?;
                residual = residual.max(step.native.solver_residual);
                last = Some((sample, step.last_auxiliary_load));
                state = endpoint;
            }
        }
        let converter = converter.finish(
            &self.core,
            &core_state(initial),
            &core_state(&state),
            &self.plan,
        )?;
        let b = self.battery.parameters();
        let m = self.motor_parameters;
        let storage_change_j = [
            self.battery.curve().energy_change_j(initial[0], state[0])?,
            quadratic(initial[1], state[1], b.polarization_capacitance_f),
            quadratic(initial[2], state[2], b.bus_capacitance_f),
            converter.storage_change_j,
            quadratic(initial[3 + C], state[3 + C], m.inductance_h),
            quadratic(initial[4 + C], state[4 + C], m.rotor_inertia_kg_m2),
            quadratic(initial[5 + C], state[5 + C], m.stator_inertia_kg_m2),
        ];
        let delta_h = storage_energy_change(
            self.system.storage(),
            initial,
            &state,
            &mut [0.0; N],
            &mut [0.0; N],
        );
        let (last_mean, last_auxiliary_load) =
            last.ok_or_else(|| invalid("empty capacitive drive cycle"))?;
        let mean = if self.plan.integration_steps == 1 {
            last_mean
        } else {
            means.map(|v| (v / Dd::from_f64(self.plan.period_s)).to_f64())
        };
        let w = work.map(Dd::to_f64);
        let h = heat.map(Dd::to_f64);
        let step = CapacitiveBatteryMotorCycleStep {
            native: PreparedStepRecord {
                delta_h,
                dissipated: number(native[1])?,
                supplied: number(native[2])?,
                newton_iters: iterations,
                solver_residual: residual,
                solver_certificate: None,
            },
            converter,
            storage_change_j,
            battery_terminal_work_j: w[0],
            motor_electrical_work_j: w[3],
            electromagnetic_work_j: w[4],
            shaft_work_j: [w[5], w[6]],
            load_work_j: w[7],
            curtailed_load_energy_j: number(curtailed)?,
            mean_bus_voltage_v: [mean[0], mean[1], mean[2]],
            mean_lead_current_a: mean[3],
            mean_motor_current_a: mean[4],
            mean_shaft_speed_rad_s: [mean[5], mean[6]],
            last_auxiliary_load,
            angular_momentum_residual_nms: (Dd::from_f64(state[4 + C])
                - Dd::from_f64(initial[4 + C])
                + Dd::from_f64(state[5 + C])
                - Dd::from_f64(initial[5 + C])
                - impulse)
                .to_f64(),
            integration_steps: self.plan.integration_steps,
            battery_heat_j: [h[0], h[1]],
            lead_heat_j: h[2],
            copper_heat_j: h[3],
            bearing_heat_j: h[4],
        };
        let changes_abs = absolute_changes.map(Dd::to_f64);
        let work_abs = absolute_work.map(Dd::to_f64);
        let c = self.plan.controls;
        step.admit(
            changes_abs,
            work_abs,
            c.absolute_energy_tolerance_j,
            c.relative_energy_tolerance,
        )?;
        for i in 0..7 {
            admit(
                number(Dd::from_f64(storage_change_j[i]) - changes[i])?,
                changes_abs[i],
                c.absolute_energy_tolerance_j,
                c.relative_energy_tolerance,
                "capacitive drive endpoint and accumulated constituent storage disagree",
            )?;
        }
        let scale = total(
            changes_abs
                .into_iter()
                .chain(work_abs)
                .chain(step.all_heat_joules()),
        );
        admit(
            number(Dd::from_f64(delta_h) - native[0])?,
            scale,
            c.absolute_energy_tolerance_j,
            c.relative_energy_tolerance,
            "capacitive drive endpoint and accumulated native storage disagree",
        )?;
        // Converter work was independently accumulated by the shared account.
        for i in 0..2 {
            admit(
                total([step.converter.terminal_work_j[i], -w[1 + i]]),
                work_abs[1 + i],
                c.absolute_energy_tolerance_j,
                c.relative_energy_tolerance,
                "capacitive drive shared converter work accumulation disagrees",
            )?;
        }
        for k in 0..4 {
            admit(
                total([step.converter.gate_source_work_j()[k], -w[8 + k]]),
                work_abs[8 + k],
                c.absolute_energy_tolerance_j,
                c.relative_energy_tolerance,
                "finite-drive gate source work accumulation disagrees",
            )?;
        }
        poll(&mut cancelled)?;
        *next = state;
        Ok(step)
    }
}
