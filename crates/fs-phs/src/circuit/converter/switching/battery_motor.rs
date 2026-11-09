//! Prescribed PWM from finite battery energy to both physical motor shafts.
//!
//! Each gate/subdivision uses the ORIGINAL nine-state battery/converter/motor
//! discrete-gradient owner. Chemical and polarization charge, all three bus
//! charges, both winding fluxes and both absolute momenta continue across every
//! edge. Signed armature current regenerates into the same finite battery.
//!
//! Work and ten heat channels retain native discrete-effort quadrature; explicit
//! subdivision is needed for temporal accuracy. This is ideal complementary
//! switching with finite on-resistance, not a dead-time, diode, transition-loss,
//! automatic control, battery-protection or semiconductor safe-area model.

use super::power_load::MAX_POWER_SUBSTEPS_PER_INTERVAL;
use super::{
    AveragedBuckBoost, BridgeModulation, BridgeSwitchState, BuckBoostParameters, Sum,
    SwitchCycleControls, SwitchInterval, admit_cycle, invalid, poll, positive,
};
use crate::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use crate::circuit::converter::battery_fed::motor::{BatteryMotorDrive, PreparedBatteryMotorDrive};
use crate::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadPoint, RegulatedPowerLoad,
};
use crate::circuit::motor::{DcMotor, DcMotorParameters};
use crate::{PreparedStepError, PreparedStepRecord, storage_energy_change};
use fs_math::dd::Dd;

/// Commands held on one ORIGINAL authored gate interval. Positive shaft torque
/// enters the named absolute rotor/stator coordinate; auxiliary demand exports
/// electrical work from the output bus and is not a drive heat channel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatteryMotorSwitchCommand {
    pub shaft_torque_nm: [f64; 2],
    pub auxiliary_load: Option<RegulatedPowerLoad>,
}

/// Complete-cycle physical accounts. Means are weighted native discrete efforts,
/// not endpoint observations or one constitutive discrete-gradient operating point.
#[derive(Clone, Debug)]
pub struct BatteryMotorSwitchCycleStep {
    /// Original endpoint storage and sums of accepted native work/heat. Iteration
    /// count is summed, residual is the maximum, and no single solve certificate
    /// represents the full cycle.
    pub native: PreparedStepRecord,
    /// Chemical, polarization, battery bus, converter, motor magnetic, rotor and
    /// stator storage changes [J], independently evaluated at the cycle endpoints.
    pub storage_change_j: [f64; 7],
    /// Positive enters the battery terminal; regeneration can make this positive.
    pub battery_terminal_work_j: f64,
    pub converter_terminal_work_j: [f64; 2],
    pub motor_electrical_work_j: f64,
    /// Electrical-to-mechanical conversion, signed, before bearing heat [J].
    pub electromagnetic_work_j: f64,
    pub shaft_work_j: [f64; 2],
    pub load_work_j: f64,
    pub curtailed_load_energy_j: f64,
    /// Battery bus, converter input, converter output [V].
    pub mean_bus_voltage_v: [f64; 3],
    pub mean_lead_current_a: f64,
    pub mean_motor_current_a: f64,
    pub mean_shaft_speed_rad_s: [f64; 2],
    /// The final accepted substep's load point; not a cycle-average branch.
    pub last_auxiliary_load: PowerLoadPoint,
    /// Independent complete-cycle angular impulse residual [N m s].
    pub angular_momentum_residual_nms: f64,
    pub integration_steps: usize,
    heat_j: [f64; 10],
}

fn total(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value))
        .to_f64()
}

impl BatteryMotorSwitchCycleStep {
    /// Converter five, battery series/polarization, lead, winding copper,
    /// bearing friction. Shaft and auxiliary-load exports never appear here.
    pub fn heat_joules(&self) -> [f64; 10] {
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
}

struct CachedStep {
    state: BridgeSwitchState,
    duration_s: f64,
    native: PreparedBatteryMotorDrive,
}

#[derive(Clone, Copy)]
struct IntervalSteps {
    regular: usize,
    last: usize,
}

/// Frozen physical laws with bounded, explicitly resolved gate intervals.
/// Identical binary connection/substep-duration pairs share native workspaces.
/// Immutable OCV data are shared by the original source owner. Successful warm
/// cycles allocate nothing and commit none of the caller's state on refusal.
pub struct PreparedBatteryMotorSwitchCycle {
    intervals: Vec<SwitchInterval>,
    interval_steps: Vec<IntervalSteps>,
    caches: Vec<CachedStep>,
    substeps_per_interval: usize,
    integration_steps: usize,
    period_s: f64,
    averaged_modulation: BridgeModulation,
    controls: SwitchCycleControls,
}

impl PreparedBatteryMotorSwitchCycle {
    /// Replace the model's averaged modulation with the explicitly authored
    /// binary schedule; all other original constitutive data remain unchanged.
    /// No physical state is stored in or reset by preparation.
    ///
    /// At most 64 physical intervals and 1..=256 subdivisions per interval are
    /// admitted. The native Newton/backtracking budget applies to every step.
    /// The last substep is a compensated remainder so every original gate edge
    /// survives exactly in the represented clock.
    pub fn new_controlled(
        model: BatteryMotorDrive,
        intervals: &[SwitchInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        let (period_s, averaged_modulation) = admit_cycle(intervals, controls, &mut cancelled)?;
        if !(1..=MAX_POWER_SUBSTEPS_PER_INTERVAL).contains(&substeps_per_interval)
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(invalid(
                "battery motor switching needs bounded subdivisions and relative allowance below one",
            ).into());
        }
        let integration_steps = intervals
            .len()
            .checked_mul(substeps_per_interval)
            .ok_or_else(|| invalid("battery motor switching step count is unrepresentable"))?;
        integration_steps
            .checked_mul(electrical.solve.max_iterations)
            .ok_or_else(|| {
                invalid("battery motor switching Newton work bound is unrepresentable")
            })?;
        let p = model.electrical().parameters();
        let b = model.electrical().battery_parameters();
        let curve = model.electrical().battery_curve();
        let lead = model.electrical().lead_resistance_ohm();
        let m = model.motor_parameters();
        let mut caches: Vec<CachedStep> = Vec::new();
        let mut interval_steps = Vec::with_capacity(intervals.len());
        for interval in intervals {
            poll(&mut cancelled)?;
            let regular_s = interval.duration_s / substeps_per_interval as f64;
            let last_s = (Dd::from_f64(interval.duration_s)
                - Dd::from_f64(regular_s) * Dd::from_f64((substeps_per_interval - 1) as f64))
            .to_f64();
            if !positive(regular_s) || !positive(last_s) {
                return Err(
                    invalid("battery motor subdivision duration is unrepresentable").into(),
                );
            }
            let mut represented = Sum::default();
            for j in 0..substeps_per_interval {
                poll(&mut cancelled)?;
                represented.add(if j + 1 == substeps_per_interval {
                    last_s
                } else {
                    regular_s
                });
            }
            if represented.total().to_bits() != interval.duration_s.to_bits() {
                return Err(
                    invalid("battery motor subdivisions must preserve each gate edge").into(),
                );
            }
            let mut indices = [0usize; 2];
            for (slot, duration_s) in [regular_s, last_s].into_iter().enumerate() {
                poll(&mut cancelled)?;
                let existing = caches.iter().position(|cache| {
                    cache.state == interval.state
                        && cache.duration_s.to_bits() == duration_s.to_bits()
                });
                indices[slot] = if let Some(index) = existing {
                    index
                } else {
                    let native = AveragedBuckBoost::new(p, interval.state.modulation())?
                        .with_battery_source(BatterySource::new(curve.clone(), b)?, lead)?
                        .with_dc_motor(DcMotor::new(m)?)?
                        .prepare_controlled(duration_s, electrical, &mut cancelled)?;
                    let index = caches.len();
                    caches.push(CachedStep {
                        state: interval.state,
                        duration_s,
                        native,
                    });
                    index
                };
            }
            interval_steps.push(IntervalSteps {
                regular: indices[0],
                last: indices[1],
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            intervals: intervals.to_vec(),
            interval_steps,
            caches,
            substeps_per_interval,
            integration_steps,
            period_s,
            averaged_modulation,
            controls,
        })
    }

    pub fn intervals(&self) -> &[SwitchInterval] {
        &self.intervals
    }
    pub const fn period_s(&self) -> f64 {
        self.period_s
    }
    pub const fn averaged_modulation(&self) -> BridgeModulation {
        self.averaged_modulation
    }
    pub const fn substeps_per_interval(&self) -> usize {
        self.substeps_per_interval
    }
    pub const fn integration_steps(&self) -> usize {
        self.integration_steps
    }
    pub fn cached_interval_count(&self) -> usize {
        self.caches.len()
    }
    pub fn parameters(&self) -> BuckBoostParameters {
        self.caches[0].native.electrical().parameters()
    }
    pub fn battery_parameters(&self) -> BatteryParameters {
        self.caches[0].native.electrical().battery_parameters()
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        self.caches[0].native.electrical().battery_curve()
    }
    pub fn lead_resistance_ohm(&self) -> f64 {
        self.caches[0].native.electrical().lead_resistance_ohm()
    }
    pub fn motor_parameters(&self) -> DcMotorParameters {
        self.caches[0].native.motor_parameters()
    }

    /// One command across the whole physical cycle. All nine caller coordinates
    /// publish together only after every native step and complete-cycle gate.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 9],
        command: BatteryMotorSwitchCommand,
        next: &mut [f64; 9],
        cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryMotorSwitchCycleStep, PreparedStepError> {
        self.advance(initial, |_| command, next, cancelled)
    }

    /// One command for each ORIGINAL authored interval, retaining gate/source
    /// correlations even when intervals share a native numerical cache.
    pub fn step_interval_commands_controlled(
        &mut self,
        initial: &[f64; 9],
        commands: &[BatteryMotorSwitchCommand],
        next: &mut [f64; 9],
        cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryMotorSwitchCycleStep, PreparedStepError> {
        if commands.len() != self.intervals.len() {
            return Err(invalid(
                "one battery motor command is required per original gate interval",
            )
            .into());
        }
        self.advance(initial, |k| commands[k], next, cancelled)
    }

    fn advance(
        &mut self,
        initial: &[f64; 9],
        command: impl Fn(usize) -> BatteryMotorSwitchCommand,
        next: &mut [f64; 9],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryMotorSwitchCycleStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let (mut heat, mut changes, mut absolute_changes) =
            ([Dd::ZERO; 10], [Dd::ZERO; 7], [Dd::ZERO; 7]);
        // Battery, converter input/output, motor electrical, electromagnetic,
        // rotor/stator shaft, auxiliary load. Signs remain physical.
        let (mut work, mut absolute_work) = ([Dd::ZERO; 8], [Dd::ZERO; 8]);
        let (mut means, mut native) = ([Dd::ZERO; 7], [Dd::ZERO; 3]);
        let (mut curtailed, mut impulse) = (Dd::ZERO, Dd::ZERO);
        let mut newton_iters = 0usize;
        let mut solver_residual = 0.0_f64;
        let mut last_means = [0.0; 7];
        let mut last_auxiliary_load = None;
        for k in 0..self.intervals.len() {
            poll(&mut cancelled)?;
            let command = command(k);
            let indices = self.interval_steps[k];
            for j in 0..self.substeps_per_interval {
                poll(&mut cancelled)?;
                let index = if j + 1 == self.substeps_per_interval {
                    indices.last
                } else {
                    indices.regular
                };
                let cache = &mut self.caches[index];
                let dt = Dd::from_f64(cache.duration_s);
                let mut endpoint = [0.0; 9];
                let step = cache.native.step_controlled(
                    &state,
                    &command.shaft_torque_nm,
                    command.auxiliary_load,
                    &mut endpoint,
                    &mut cancelled,
                )?;
                let e = &step.electrical;
                let delta = [
                    e.battery.chemical_energy_change_j,
                    e.battery.polarization_energy_change_j,
                    e.battery.bus_energy_change_j,
                    e.converter_storage_change_j,
                    step.magnetic_energy_change_j,
                    step.kinetic_energy_change_j[0],
                    step.kinetic_energy_change_j[1],
                ];
                for i in 0..7 {
                    changes[i] = changes[i] + Dd::from_f64(delta[i]);
                    absolute_changes[i] = absolute_changes[i] + Dd::from_f64(delta[i].abs());
                }
                for (sum, value) in heat.iter_mut().zip(step.heat_joules()) {
                    *sum = *sum + Dd::from_f64(value);
                }
                let em_work = (dt
                    * Dd::from_f64(step.discrete_back_emf_v)
                    * Dd::from_f64(step.discrete_motor_current_a))
                .to_f64();
                let transfers = [
                    e.battery.terminal_work_j,
                    e.converter_terminal_work_j[0],
                    e.converter_terminal_work_j[1],
                    step.motor_electrical_work_j,
                    em_work,
                    step.shaft_work_j[0],
                    step.shaft_work_j[1],
                    e.load_work_j,
                ];
                for i in 0..8 {
                    work[i] = work[i] + Dd::from_f64(transfers[i]);
                    absolute_work[i] = absolute_work[i] + Dd::from_f64(transfers[i].abs());
                }
                last_means = [
                    e.discrete_bus_voltage_v[0],
                    e.discrete_bus_voltage_v[1],
                    e.discrete_bus_voltage_v[2],
                    e.lead_current_a,
                    step.discrete_motor_current_a,
                    step.discrete_shaft_speed_rad_s[0],
                    step.discrete_shaft_speed_rad_s[1],
                ];
                for i in 0..7 {
                    means[i] = means[i] + dt * Dd::from_f64(last_means[i]);
                }
                curtailed = curtailed + Dd::from_f64(e.curtailed_load_energy_j);
                impulse = impulse
                    + dt * (Dd::from_f64(command.shaft_torque_nm[0])
                        + Dd::from_f64(command.shaft_torque_nm[1]));
                for (sum, value) in native.iter_mut().zip([
                    step.native.delta_h,
                    step.native.dissipated,
                    step.native.supplied,
                ]) {
                    *sum = *sum + Dd::from_f64(value);
                }
                newton_iters = newton_iters
                    .checked_add(step.native.newton_iters)
                    .ok_or_else(|| invalid("battery motor switching Newton count overflow"))?;
                solver_residual = solver_residual.max(step.native.solver_residual);
                last_auxiliary_load = Some(e.load);
                state = endpoint;
            }
        }
        let p = self.parameters();
        let b = self.battery_parameters();
        let m = self.motor_parameters();
        let quadratic = |i: usize, coefficient: f64| {
            let a = Dd::from_f64(initial[i]);
            let z = Dd::from_f64(state[i]);
            ((z - a) * (Dd::from_f64(0.5) * a + Dd::from_f64(0.5) * z) / Dd::from_f64(coefficient))
                .to_f64()
        };
        let storage_change_j = [
            self.battery_curve().energy_change_j(initial[0], state[0])?,
            quadratic(1, b.polarization_capacitance_f),
            quadratic(2, b.bus_capacitance_f),
            total([
                quadratic(3, p.inductance_h),
                quadratic(4, p.bus_capacitance_f[0]),
                quadratic(5, p.bus_capacitance_f[1]),
            ]),
            quadratic(6, m.inductance_h),
            quadratic(7, m.rotor_inertia_kg_m2),
            quadratic(8, m.stator_inertia_kg_m2),
        ];
        let delta_h = storage_energy_change(
            self.caches[0].native.system().storage(),
            initial,
            &state,
            &mut [0.0; 9],
            &mut [0.0; 9],
        );
        let heat_j = heat.map(Dd::to_f64);
        let work = work.map(Dd::to_f64);
        let absolute_work = absolute_work.map(Dd::to_f64);
        let absolute_changes = absolute_changes.map(Dd::to_f64);
        let means = if self.integration_steps == 1 {
            last_means
        } else {
            means.map(|v| (v / Dd::from_f64(self.period_s)).to_f64())
        };
        let native = native.map(Dd::to_f64);
        let step = BatteryMotorSwitchCycleStep {
            native: PreparedStepRecord {
                delta_h,
                dissipated: native[1],
                supplied: native[2],
                newton_iters,
                solver_residual,
                solver_certificate: None,
            },
            storage_change_j,
            battery_terminal_work_j: work[0],
            converter_terminal_work_j: [work[1], work[2]],
            motor_electrical_work_j: work[3],
            electromagnetic_work_j: work[4],
            shaft_work_j: [work[5], work[6]],
            load_work_j: work[7],
            curtailed_load_energy_j: curtailed.to_f64(),
            mean_bus_voltage_v: [means[0], means[1], means[2]],
            mean_lead_current_a: means[3],
            mean_motor_current_a: means[4],
            mean_shaft_speed_rad_s: [means[5], means[6]],
            last_auxiliary_load: last_auxiliary_load
                .ok_or_else(|| invalid("battery motor switching cycle has no accepted step"))?,
            angular_momentum_residual_nms: (Dd::from_f64(state[7]) - Dd::from_f64(initial[7])
                + Dd::from_f64(state[8])
                - Dd::from_f64(initial[8])
                - impulse)
                .to_f64(),
            integration_steps: self.integration_steps,
            heat_j,
        };
        let s = storage_change_j;
        let battery_scale = total(absolute_changes[..3].iter().copied().chain([
            heat_j[5],
            heat_j[6],
            absolute_work[0],
        ]));
        let lead_scale = total([absolute_work[0], absolute_work[1], heat_j[7]]);
        let converter_scale = total(heat_j[..5].iter().copied().chain([
            absolute_changes[3],
            absolute_work[1],
            absolute_work[2],
        ]));
        let magnetic_scale = total([
            absolute_changes[4],
            heat_j[8],
            absolute_work[3],
            absolute_work[4],
        ]);
        let shaft_scale = total([
            absolute_changes[5],
            absolute_changes[6],
            heat_j[9],
            absolute_work[4],
            absolute_work[5],
            absolute_work[6],
        ]);
        let scale = total(
            absolute_changes
                .into_iter()
                .chain(absolute_work)
                .chain(heat_j),
        );
        let controls = self.controls;
        let gate = |residual: f64, scale: f64| -> Result<(), PreparedStepError> {
            let tolerance =
                controls.absolute_energy_tolerance_j + controls.relative_energy_tolerance * scale;
            if !scale.is_finite()
                || !tolerance.is_finite()
                || !residual.is_finite()
                || residual.abs() > tolerance
            {
                return Err(
                    invalid("battery motor switching physical boundary does not close").into(),
                );
            }
            Ok(())
        };
        let residuals = [
            (
                total([s[0], s[1], s[2], heat_j[5], heat_j[6], -work[0]]),
                battery_scale,
            ),
            (total([work[0], work[1], heat_j[7]]), lead_scale),
            (
                total(
                    [s[3], -work[1], -work[2]]
                        .into_iter()
                        .chain(heat_j[..5].iter().copied()),
                ),
                converter_scale,
            ),
            (
                total([work[2], work[3], work[7]]),
                total([absolute_work[2], absolute_work[3], absolute_work[7]]),
            ),
            (total([s[4], heat_j[8], work[4], -work[3]]), magnetic_scale),
            (
                total([s[5], s[6], heat_j[9], -work[4], -work[5], -work[6]]),
                shaft_scale,
            ),
            (total(s.into_iter().chain([-delta_h])), scale),
            (total([delta_h, -native[0]]), scale),
            (total(heat_j.into_iter().chain([-native[1]])), total(heat_j)),
            (
                total([native[2], -work[5], -work[6], work[7]]),
                total([absolute_work[5], absolute_work[6], absolute_work[7]]),
            ),
            (step.balance_residual_j(), scale),
        ];
        for (residual, scale) in residuals {
            gate(residual, scale)?;
        }
        let component_scales = [
            battery_scale,
            battery_scale,
            battery_scale,
            converter_scale,
            magnetic_scale,
            shaft_scale,
            shaft_scale,
        ];
        for i in 0..7 {
            gate(
                (Dd::from_f64(s[i]) - changes[i]).to_f64(),
                component_scales[i],
            )?;
        }
        if state
            .iter()
            .chain(&heat_j)
            .chain(&work)
            .chain(&means)
            .chain(&s)
            .chain(&absolute_work)
            .chain(&absolute_changes)
            .any(|v| !v.is_finite())
            || heat_j.iter().any(|v| *v < 0.0)
            || !step.curtailed_load_energy_j.is_finite()
            || step.curtailed_load_energy_j < 0.0
            || !step.angular_momentum_residual_nms.is_finite()
            || !solver_residual.is_finite()
            || step.load_work_j < 0.0
        {
            return Err(
                invalid("battery motor switching physical account is unrepresentable").into(),
            );
        }
        poll(&mut cancelled)?;
        *next = state;
        Ok(step)
    }
}
