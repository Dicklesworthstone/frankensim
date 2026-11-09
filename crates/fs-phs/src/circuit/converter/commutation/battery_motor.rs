//! Finite battery energy through nonideal bridge commutation to both motor shafts.
//!
//! The original nine-state battery/lead/converter/motor interconnection retains
//! all storage and reciprocal motor reactions. Only its averaged bridge J and
//! channel R are replaced by the existing static HalfBridgeConduction laws.
//! A single native workspace resolves chemical depletion, polarization, every
//! gate's diode/channel currents, winding current and both momenta together.
//! No sampled ideal supply, sequential subsystem solve or efficiency map is used.

use super::{
    CommutationInterval, CommutationParameters, DeviceConductionPoint, DeviceGateState,
    HalfBridgeConduction, HalfBridgePoint, SwitchDevice, invalid, number, poll, positive, total,
};
use crate::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use crate::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use crate::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadPoint, RegulatedPowerLoad, load_point,
};
use crate::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use crate::circuit::converter::switching::power_load::MAX_POWER_SUBSTEPS_PER_INTERVAL;
use crate::circuit::converter::switching::{MAX_SWITCH_INTERVALS, SwitchCycleControls};
use crate::circuit::motor::{DcMotor, DcMotorParameters};
use crate::{
    NonlinearPortPower, PhsError, PreparedStepError, PreparedStepRecord, StepWorkspace,
    storage_energy_change,
};
use fs_math::dd::Dd;

mod cycle;

/// Original constituent balances on one COMPLETE authored independent-gate cycle.
/// All means are weighted accepted discrete efforts; last points are not means.
#[derive(Clone, Debug)]
pub struct CommutatedBatteryMotorCycleStep {
    /// Endpoint storage, sums of actual native work/heat, summed Newton count,
    /// maximum substep residual; no single solver certificate is fabricated.
    pub native: PreparedStepRecord,
    /// Chemical, polarization, battery bus, converter, motor magnetic, rotor,
    /// stator storage changes, evaluated independently at the original endpoints.
    pub storage_change_j: [f64; 7],
    /// Positive work enters the physical battery terminal, including regeneration.
    pub battery_terminal_work_j: f64,
    pub converter_terminal_work_j: [f64; 2],
    pub motor_electrical_work_j: f64,
    pub electromagnetic_work_j: f64,
    pub shaft_work_j: [f64; 2],
    pub load_work_j: f64,
    pub curtailed_load_energy_j: f64,
    /// Battery bus, converter input bus, converter output bus.
    pub mean_bus_voltage_v: [f64; 3],
    pub mean_lead_current_a: f64,
    pub mean_motor_current_a: f64,
    pub mean_shaft_speed_rad_s: [f64; 2],
    pub last_auxiliary_load: PowerLoadPoint,
    /// Input high/low then output high/low at the last accepted DISCRETE point.
    /// Every substep's initial, discrete and final device limits are checked.
    pub last_devices: [DeviceConductionPoint; 4],
    pub angular_momentum_residual_nms: f64,
    pub integration_steps: usize,
    heat_j: [f64; 14],
}

impl CommutatedBatteryMotorCycleStep {
    /// Converter winding/four channels, battery series/polarization, lead,
    /// motor copper, bearing, then four body diodes in the same device order.
    /// Neither auxiliary load delivery nor external shaft work is drive heat.
    pub fn heat_joules(&self) -> [f64; 14] {
        self.heat_j
    }
    pub fn total_storage_change_j(&self) -> f64 {
        total(self.storage_change_j)
    }
    pub fn balance_residual_j(&self) -> f64 {
        total(self.storage_change_j.into_iter().chain(self.heat_j).chain([
            self.load_work_j,
            -self.shaft_work_j[0],
            -self.shaft_work_j[1],
        ]))
    }
    fn work(&self) -> [f64; 8] {
        [
            self.battery_terminal_work_j,
            self.converter_terminal_work_j[0],
            self.converter_terminal_work_j[1],
            self.motor_electrical_work_j,
            self.electromagnetic_work_j,
            self.shaft_work_j[0],
            self.shaft_work_j[1],
            self.load_work_j,
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

    // Original constituent equations and interconnection exchanges, not just a
    // single final energy residual. Absolute transfers retain cycle cancellation.
    fn admit(
        &self,
        absolute_changes: [f64; 7],
        absolute_work: [f64; 8],
        absolute_allowance_j: f64,
        relative_allowance: f64,
    ) -> Result<(), PhsError> {
        let h = self.heat_j;
        let s = self.storage_change_j;
        let w = self.work();
        let bridge_heat = total(h[..5].iter().copied().chain(h[10..].iter().copied()));
        let heat = total(h);
        let battery_scale =
            total(
                absolute_changes[..3]
                    .iter()
                    .copied()
                    .chain([absolute_work[0], h[5], h[6]]),
            );
        let bridge_scale = total([
            absolute_changes[3],
            bridge_heat,
            absolute_work[1],
            absolute_work[2],
        ]);
        let lead_scale = total([absolute_work[0], absolute_work[1], h[7]]);
        let output_scale = total([absolute_work[2], absolute_work[3], absolute_work[7]]);
        let magnetic_scale = total([
            absolute_changes[4],
            h[8],
            absolute_work[3],
            absolute_work[4],
        ]);
        let mechanical_scale = total([
            absolute_changes[5],
            absolute_changes[6],
            h[9],
            absolute_work[4],
            absolute_work[5],
            absolute_work[6],
        ]);
        let scale = battery_scale
            .max(bridge_scale)
            .max(lead_scale)
            .max(output_scale)
            .max(magnetic_scale)
            .max(mechanical_scale)
            .max(heat);
        let balances = [
            (total([s[0], s[1], s[2], h[5], h[6], -w[0]]), battery_scale),
            (total([s[3], bridge_heat, -w[1], -w[2]]), bridge_scale),
            (total([-w[0], -w[1], -h[7]]), lead_scale),
            (total([w[2], w[3], w[7]]), output_scale),
            (total([s[4], h[8], w[4], -w[3]]), magnetic_scale),
            (
                total([s[5], s[6], h[9], -w[4], -w[5], -w[6]]),
                mechanical_scale,
            ),
            (total(s.into_iter().chain([-self.native.delta_h])), scale),
            (total([heat, -self.native.dissipated]), heat),
            (
                total([self.native.supplied, -w[5], -w[6], w[7]]),
                total([absolute_work[5], absolute_work[6], absolute_work[7]]),
            ),
            (self.balance_residual_j(), scale),
        ];
        if s.iter()
            .chain(&w)
            .chain(&self.means())
            .any(|v| !v.is_finite())
            || h.iter()
                .chain(&[self.load_work_j, self.curtailed_load_energy_j])
                .any(|v| !v.is_finite() || *v < 0.0)
            || !self.angular_momentum_residual_nms.is_finite()
            || !self.native.solver_residual.is_finite()
            || balances.iter().any(|&(residual, magnitude)| {
                let allowance = absolute_allowance_j + relative_allowance * magnitude;
                !magnitude.is_finite()
                    || !allowance.is_finite()
                    || !residual.is_finite()
                    || residual.abs() > allowance
            })
        {
            return Err(invalid(
                "finite battery, device, lead, winding and shaft accounts do not close on one trajectory",
            ));
        }
        Ok(())
    }
}

struct CachedInterval {
    gates: DeviceGateState,
    duration_s: f64,
    legs: [HalfBridgeConduction; 2],
}
#[derive(Clone, Copy)]
struct IntervalSteps {
    regular: usize,
    last: usize,
}

/// Frozen physical laws; all original accepted state remains in caller storage.
///
/// The SAME nine-state workspace serves every gate and duration. Its fixed
/// battery/lead/winding/motor interconnection does not change when a gate changes;
/// bounded static leg caches only provide the nonlinear commutation response.
/// Warm cycles allocate nothing. A failure or cancellation publishes no state.
pub struct PreparedCommutatedBatteryMotorCycle {
    model: BatteryMotorDrive,
    parameters: CommutationParameters,
    workspace: StepWorkspace,
    intervals: Vec<CommutationInterval>,
    interval_steps: Vec<IntervalSteps>,
    caches: Vec<CachedInterval>,
    substeps_per_interval: usize,
    integration_steps: usize,
    period_s: f64,
    electrical: PowerLoadControls,
    controls: SwitchCycleControls,
}

fn quadratic_change(before: f64, after: f64, coefficient: f64) -> f64 {
    let a = Dd::from_f64(before);
    let b = Dd::from_f64(after);
    ((b - a) * (Dd::from_f64(0.5) * a + Dd::from_f64(0.5) * b) / Dd::from_f64(coefficient)).to_f64()
}
fn product(a: f64, b: f64, c: f64) -> f64 {
    (Dd::from_f64(a) * Dd::from_f64(b) * Dd::from_f64(c)).to_f64()
}

struct OperatingPoint {
    legs: [HalfBridgePoint; 2],
    load: PowerLoadPoint,
}
fn point(
    legs: &[HalfBridgeConduction; 2],
    effort: &[f64],
    auxiliary: Option<RegulatedPowerLoad>,
) -> Result<OperatingPoint, PhsError> {
    if effort.len() != 9 || effort.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "commutated battery drive requires finite original nine-state efforts",
        ));
    }
    Ok(OperatingPoint {
        legs: [
            legs[0].evaluate(effort[4], effort[3])?,
            legs[1].evaluate(effort[5], -effort[3])?,
        ],
        load: load_point(auxiliary, effort[5])?,
    })
}

impl PreparedCommutatedBatteryMotorCycle {
    fn endpoint_support(
        parameters: CommutationParameters,
        legs: &[HalfBridgeConduction; 2],
        state: &[f64; 9],
        auxiliary: Option<RegulatedPowerLoad>,
    ) -> Result<(), PhsError> {
        if state.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "commutated battery drive requires finite original states",
            ));
        }
        let mut effort = [0.0; 9];
        effort[3] = state[3] / parameters.inductance_h;
        effort[4] = state[4] / parameters.bus_capacitance_f[0];
        effort[5] = state[5] / parameters.bus_capacitance_f[1];
        point(legs, &effort, auxiliary)?;
        Ok(())
    }

    // One shared original nine-state Newton solve. The interval receipt remains
    // private; no constituent subsystem receives a fictitious solver certificate.
    #[allow(clippy::too_many_arguments)]
    fn step_native(
        model: &BatteryMotorDrive,
        parameters: CommutationParameters,
        workspace: &mut StepWorkspace,
        cache: &CachedInterval,
        controls: PowerLoadControls,
        initial: &[f64; 9],
        command: BatteryMotorSwitchCommand,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; 9], CommutatedBatteryMotorCycleStep), PreparedStepError> {
        poll(cancelled)?;
        Self::endpoint_support(parameters, &cache.legs, initial, command.auxiliary_load)?;
        let legs = cache.legs;
        let auxiliary = command.auxiliary_load;
        let port = |_: &[f64], effort: &[f64], out: &mut [f64]| {
            if out.len() != 9 {
                return None;
            }
            let p = point(&legs, effort, auxiliary).ok()?;
            out.fill(0.0);
            out[3] = number(p.legs[1].node - p.legs[0].node).ok()?;
            out[4] = number(p.legs[0].high_current).ok()?;
            out[5] = number(p.legs[1].high_current + Dd::from_f64(p.load.current_a)).ok()?;
            Some(NonlinearPortPower {
                supplied_power_w: -p.load.delivered_power_w,
                dissipated_power_w: number(p.legs[0].loss + p.legs[1].loss).ok()?,
            })
        };
        let tangent = |_: &[f64], effort: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            if de.len() != 9 || out.len() != 9 {
                return false;
            }
            let Ok(p) = point(&legs, effort, auxiliary) else {
                return false;
            };
            out.fill(0.0);
            let mut nodes = [0.0; 2];
            for k in 0..2 {
                let draw = if k == 0 { de[3] } else { -de[3] };
                nodes[k] = total([
                    p.legs[k].d_node_voltage_d_bus_voltage * de[k + 4],
                    p.legs[k].d_node_voltage_d_draw_current_ohm * draw,
                ]);
                out[k + 4] = total([
                    p.legs[k].d_high_current_d_bus_voltage_s * de[k + 4],
                    p.legs[k].d_high_current_d_draw_current * draw,
                    if k == 1 {
                        p.load.incremental_conductance_s * de[5]
                    } else {
                        0.0
                    },
                ]);
            }
            out[3] = total([nodes[1], -nodes[0]]);
            out.iter().all(|v| v.is_finite())
        };
        let dt = cache.duration_s;
        let (mut candidate, mut outputs) = ([0.0; 9], [0.0; 4]);
        let native = workspace.step_into_supplied_controlled(
            model.system(),
            initial,
            &[
                0.0,
                0.0,
                command.shaft_torque_nm[0],
                command.shaft_torque_nm[1],
            ],
            dt,
            &mut candidate,
            &mut outputs,
            &port,
            None,
            Some(&tangent),
            &mut *cancelled,
        )?;
        let e = workspace
            .last_step_discrete_effort()
            .ok_or_else(|| invalid("commutated battery drive requires its actual native effort"))?;
        let operating = point(&legs, e, auxiliary)?;
        Self::endpoint_support(parameters, &legs, &candidate, auxiliary)?;
        let electrical = model.electrical().account_interval(
            &core::array::from_fn(|i| initial[i]),
            &core::array::from_fn(|i| candidate[i]),
            &core::array::from_fn(|i| e[i]),
            -e[6],
            auxiliary,
            dt,
        )?;
        let m = model.motor_parameters();
        let relative_speed = total([e[7], -e[8]]);
        let devices = [
            operating.legs[0].high,
            operating.legs[0].low,
            operating.legs[1].high,
            operating.legs[1].low,
        ];
        let mut heat_j = [0.0; 14];
        heat_j[0] = electrical.converter_conduction.joules[0];
        for i in 0..4 {
            heat_j[1 + i] = (Dd::from_f64(dt) * Dd::from_f64(devices[i].channel_power_w)).to_f64();
            heat_j[10 + i] = (Dd::from_f64(dt) * Dd::from_f64(devices[i].diode_power_w)).to_f64();
        }
        heat_j[5] = electrical.battery.series_heat_j;
        heat_j[6] = electrical.battery.polarization_heat_j;
        heat_j[7] = electrical.lead_heat_j;
        heat_j[8] = (Dd::from_f64(dt)
            * Dd::from_f64(m.winding_resistance_ohm)
            * Dd::from_f64(e[6])
            * Dd::from_f64(e[6]))
        .to_f64();
        heat_j[9] = (Dd::from_f64(dt)
            * Dd::from_f64(m.bearing_drag_nm_s_per_rad)
            * Dd::from_f64(relative_speed)
            * Dd::from_f64(relative_speed))
        .to_f64();
        let step = CommutatedBatteryMotorCycleStep {
            native,
            storage_change_j: [
                electrical.battery.chemical_energy_change_j,
                electrical.battery.polarization_energy_change_j,
                electrical.battery.bus_energy_change_j,
                electrical.converter_storage_change_j,
                quadratic_change(initial[6], candidate[6], m.inductance_h),
                quadratic_change(initial[7], candidate[7], m.rotor_inertia_kg_m2),
                quadratic_change(initial[8], candidate[8], m.stator_inertia_kg_m2),
            ],
            battery_terminal_work_j: electrical.battery.terminal_work_j,
            converter_terminal_work_j: electrical.converter_terminal_work_j,
            motor_electrical_work_j: product(dt, e[5], e[6]),
            electromagnetic_work_j: (Dd::from_f64(dt)
                * Dd::from_f64(m.coupling_nm_per_a)
                * Dd::from_f64(relative_speed)
                * Dd::from_f64(e[6]))
            .to_f64(),
            shaft_work_j: [
                product(dt, command.shaft_torque_nm[0], e[7]),
                product(dt, command.shaft_torque_nm[1], e[8]),
            ],
            load_work_j: electrical.load_work_j,
            curtailed_load_energy_j: electrical.curtailed_load_energy_j,
            mean_bus_voltage_v: electrical.discrete_bus_voltage_v,
            mean_lead_current_a: electrical.lead_current_a,
            mean_motor_current_a: e[6],
            mean_shaft_speed_rad_s: [e[7], e[8]],
            last_auxiliary_load: operating.load,
            last_devices: devices,
            angular_momentum_residual_nms: (Dd::from_f64(candidate[7]) - Dd::from_f64(initial[7])
                + Dd::from_f64(candidate[8])
                - Dd::from_f64(initial[8])
                - Dd::from_f64(dt)
                    * (Dd::from_f64(command.shaft_torque_nm[0])
                        + Dd::from_f64(command.shaft_torque_nm[1])))
            .to_f64(),
            integration_steps: 1,
            heat_j,
        };
        step.admit(
            step.storage_change_j.map(f64::abs),
            step.work().map(f64::abs),
            controls.absolute_energy_tolerance_j,
            controls.relative_energy_tolerance,
        )?;
        poll(cancelled)?;
        Ok((candidate, step))
    }
}
