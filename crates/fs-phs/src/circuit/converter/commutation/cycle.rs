//! Bounded prescribed gate cycles on the shared commutated three-state owner.
use super::*;
use crate::circuit::converter::switching::power_load::{
    MAX_POWER_SUBSTEPS_PER_INTERVAL, PowerSwitchCommand,
};
use crate::circuit::converter::switching::{MAX_SWITCH_INTERVALS, SwitchCycleControls};
use crate::storage_energy_change;

struct CachedInterval {
    gates: DeviceGateState,
    duration_s: f64,
    native: PreparedInterval,
}
#[derive(Clone, Copy)]
struct IntervalSteps {
    regular: usize,
    last: usize,
}

/// Full original commutation cycles with explicit numerical refinement.
///
/// All device/source coefficients and limits are frozen over a cycle. Authored
/// gate edges remain in order, including dead time and deliberate overlap.
/// Every substep checks actual initial/final/discrete device support. Numerical
/// refinement is still required for unobserved peaks or stiff branch crossings.
pub struct PreparedCommutatedSwitchCycle {
    parameters: CommutationParameters,
    terminal_resistance_ohm: [f64; 2],
    intervals: Vec<CommutationInterval>,
    interval_steps: Vec<IntervalSteps>,
    caches: Vec<CachedInterval>,
    substeps_per_interval: usize,
    integration_steps: usize,
    period_s: f64,
    controls: SwitchCycleControls,
}

impl PreparedCommutatedSwitchCycle {
    #[allow(clippy::too_many_arguments)]
    pub fn new_controlled(
        parameters: CommutationParameters,
        terminal_resistance_ohm: [f64; 2],
        intervals: &[CommutationInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        poll(&mut cancelled)?;
        parameters.validate()?;
        VoltageSources::new([0.0; 2], terminal_resistance_ohm)?;
        if !(1..=MAX_SWITCH_INTERVALS).contains(&controls.max_intervals)
            || intervals.is_empty()
            || intervals.len() > controls.max_intervals
            || !(1..=MAX_POWER_SUBSTEPS_PER_INTERVAL).contains(&substeps_per_interval)
            || !controls.absolute_energy_tolerance_j.is_finite()
            || controls.absolute_energy_tolerance_j < 0.0
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
        {
            return Err(invalid("bounded gate intervals, subdivisions and finite commutation energy allowances required").into());
        }
        let integration_steps = intervals
            .len()
            .checked_mul(substeps_per_interval)
            .ok_or_else(|| invalid("commutation integration count overflow"))?;
        integration_steps
            .checked_mul(electrical.solve.max_iterations)
            .and_then(|count| count.checked_mul(usize::from(electrical.backtracking_halvings) + 1))
            .ok_or_else(|| invalid("commutation native nonlinear work bound overflow"))?;
        let mut period = Dd::ZERO;
        let mut caches: Vec<CachedInterval> = Vec::new();
        let mut interval_steps = Vec::with_capacity(intervals.len());
        for interval in intervals {
            poll(&mut cancelled)?;
            if !positive(interval.duration_s) {
                return Err(invalid("positive finite physical gate duration required").into());
            }
            period = period + Dd::from_f64(interval.duration_s);
            let regular = interval.duration_s / substeps_per_interval as f64;
            let last = number(
                Dd::from_f64(interval.duration_s)
                    - Dd::from_f64(regular) * Dd::from_f64((substeps_per_interval - 1) as f64),
            )?;
            if !positive(regular) || !positive(last) {
                return Err(invalid("commutation subdivision duration is unrepresentable").into());
            }
            let mut represented = Dd::ZERO;
            for j in 0..substeps_per_interval {
                poll(&mut cancelled)?;
                represented = represented
                    + Dd::from_f64(if j + 1 == substeps_per_interval {
                        last
                    } else {
                        regular
                    });
            }
            if number(represented)?.to_bits() != interval.duration_s.to_bits() {
                return Err(
                    invalid("commutation subdivisions must retain the original gate edge").into(),
                );
            }
            let mut indices = [0; 2];
            for (slot, duration_s) in [regular, last].into_iter().enumerate() {
                poll(&mut cancelled)?;
                let found = caches.iter().position(|c| {
                    c.gates == interval.gates && c.duration_s.to_bits() == duration_s.to_bits()
                });
                indices[slot] = if let Some(index) = found {
                    index
                } else {
                    let native = PreparedInterval::new(
                        parameters,
                        interval.gates,
                        duration_s,
                        electrical,
                        &mut cancelled,
                    )?;
                    let index = caches.len();
                    caches.push(CachedInterval {
                        gates: interval.gates,
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
        let period_s = number(period)?;
        if !positive(period_s) {
            return Err(invalid("commutation period is unrepresentable").into());
        }
        poll(&mut cancelled)?;
        Ok(Self {
            parameters,
            terminal_resistance_ohm,
            intervals: intervals.to_vec(),
            interval_steps,
            caches,
            substeps_per_interval,
            integration_steps,
            period_s,
            controls,
        })
    }
    pub const fn parameters(&self) -> CommutationParameters {
        self.parameters
    }
    pub const fn terminal_resistance_ohm(&self) -> [f64; 2] {
        self.terminal_resistance_ohm
    }
    pub fn intervals(&self) -> &[CommutationInterval] {
        &self.intervals
    }
    pub const fn period_s(&self) -> f64 {
        self.period_s
    }
    pub const fn integration_steps(&self) -> usize {
        self.integration_steps
    }
    pub const fn substeps_per_interval(&self) -> usize {
        self.substeps_per_interval
    }
    pub fn cached_interval_count(&self) -> usize {
        self.caches.len()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        source_voltage_v: &[f64; 2],
        loads: &[Option<RegulatedPowerLoad>; 2],
        next: &mut [f64; 3],
        discrete_mean_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedStep, PreparedStepError> {
        let command = PowerSwitchCommand {
            source_voltage_v: *source_voltage_v,
            loads: *loads,
        };
        self.advance(
            initial,
            |_| command,
            next,
            discrete_mean_bus_voltage_v,
            cancelled,
        )
    }

    /// One command per ORIGINAL gate interval. Work retains source/current and
    /// gate/load correlation; cache sharing never combines command histories.
    pub fn step_interval_commands_controlled(
        &mut self,
        initial: &[f64; 3],
        commands: &[PowerSwitchCommand],
        next: &mut [f64; 3],
        discrete_mean_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedStep, PreparedStepError> {
        if commands.len() != self.intervals.len() {
            return Err(invalid(
                "one source/load command is required per original commutation interval",
            )
            .into());
        }
        self.advance(
            initial,
            |k| commands[k],
            next,
            discrete_mean_bus_voltage_v,
            cancelled,
        )
    }

    fn advance(
        &mut self,
        initial: &[f64; 3],
        command: impl Fn(usize) -> PowerSwitchCommand,
        next: &mut [f64; 3],
        discrete_mean_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let (mut winding, mut change, mut absolute_change, mut supplied, mut native_heat) =
            (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
        let (mut channels, mut diodes) = ([Dd::ZERO; 4], [Dd::ZERO; 4]);
        let (mut source_work, mut terminal_work, mut source_heat, mut load_work, mut curtailed) = (
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
        );
        let (mut voltage_integral, mut source_charge, mut transfers) =
            ([Dd::ZERO; 2], [Dd::ZERO; 2], [Dd::ZERO; 2]);
        let (mut iterations, mut solver_residual) = (0usize, 0.0_f64);
        let (mut last_voltage, mut last_current, mut last_points) = ([0.0; 2], [0.0; 2], None);
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
                let (mut endpoint, mut voltage) = ([0.0; 3], [0.0; 2]);
                let step = cache.native.step(
                    &state,
                    &command.source_voltage_v,
                    &self.terminal_resistance_ohm,
                    &command.loads,
                    &mut endpoint,
                    &mut voltage,
                    &mut cancelled,
                )?;
                winding = winding + Dd::from_f64(step.winding_heat_j);
                for d in 0..4 {
                    channels[d] = channels[d] + Dd::from_f64(step.channel_heat_j[d]);
                    diodes[d] = diodes[d] + Dd::from_f64(step.diode_heat_j[d]);
                }
                for p in 0..2 {
                    source_work[p] = source_work[p] + Dd::from_f64(step.source_work_j[p]);
                    terminal_work[p] = terminal_work[p] + Dd::from_f64(step.terminal_work_j[p]);
                    source_heat[p] =
                        source_heat[p] + Dd::from_f64(step.terminal_resistor_heat_j[p]);
                    load_work[p] = load_work[p] + Dd::from_f64(step.load_work_j[p]);
                    curtailed[p] = curtailed[p] + Dd::from_f64(step.curtailed_load_energy_j[p]);
                    voltage_integral[p] = voltage_integral[p]
                        + Dd::from_f64(cache.duration_s) * Dd::from_f64(voltage[p]);
                    source_charge[p] = source_charge[p]
                        + Dd::from_f64(cache.duration_s)
                            * Dd::from_f64(step.mean_source_current_a[p]);
                    transfers[p] = transfers[p]
                        + Dd::from_f64(step.source_work_j[p].abs())
                        + Dd::from_f64(step.terminal_work_j[p].abs())
                        + Dd::from_f64(step.load_work_j[p]);
                }
                change = change + Dd::from_f64(step.native.delta_h);
                absolute_change = absolute_change + Dd::from_f64(step.native.delta_h.abs());
                native_heat = native_heat + Dd::from_f64(step.native.dissipated);
                supplied = supplied + Dd::from_f64(step.native.supplied);
                iterations = iterations
                    .checked_add(step.native.newton_iters)
                    .ok_or_else(|| invalid("commutation Newton count overflow"))?;
                if !step.native.solver_residual.is_finite() {
                    return Err(invalid("nonfinite native commutation residual").into());
                }
                solver_residual = solver_residual.max(step.native.solver_residual);
                state = endpoint;
                last_voltage = voltage;
                last_current = step.mean_source_current_a;
                last_points = Some((step.last_devices, step.last_loads));
            }
        }
        let delta_h = storage_energy_change(
            self.caches[0].native.model.system().storage(),
            initial,
            &state,
            &mut [0.0; 3],
            &mut [0.0; 3],
        );
        let winding_heat_j = number(winding)?;
        let channel_heat_j = channels.map(Dd::to_f64);
        let diode_heat_j = diodes.map(Dd::to_f64);
        let source_work_j = source_work.map(Dd::to_f64);
        let terminal_work_j = terminal_work.map(Dd::to_f64);
        let terminal_resistor_heat_j = source_heat.map(Dd::to_f64);
        let load_work_j = load_work.map(Dd::to_f64);
        let curtailed_load_energy_j = curtailed.map(Dd::to_f64);
        let physical_heat = channels
            .into_iter()
            .chain(diodes)
            .fold(winding, |s, h| s + h);
        let heat = number(physical_heat)?;
        let supplied = number(supplied)?;
        let dissipated = number(native_heat)?;
        let transfer_scale = number(transfers[0] + transfers[1])?;
        let scale = transfer_scale
            .max(number(absolute_change)?)
            .max(heat)
            .max(number(source_heat[0] + source_heat[1])?)
            .max(delta_h.abs());
        let allowance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        let residuals = [
            total([delta_h, heat, -supplied]),
            total([dissipated, -heat]),
            total(terminal_work_j.into_iter().chain([-supplied])),
            number(Dd::from_f64(delta_h) - change)?,
            total(
                [delta_h, heat]
                    .into_iter()
                    .chain(terminal_resistor_heat_j)
                    .chain(load_work_j)
                    .chain(source_work_j.map(|x| -x)),
            ),
        ];
        for p in 0..2 {
            let tolerance = self.controls.absolute_energy_tolerance_j
                + self.controls.relative_energy_tolerance * number(transfers[p] + source_heat[p])?;
            let residual =
                number(source_work[p] - source_heat[p] - load_work[p] - terminal_work[p])?;
            if !tolerance.is_finite() || residual.abs() > tolerance {
                return Err(
                    invalid("commutation cycle source/resistor/load work does not close").into(),
                );
            }
        }
        let (voltage, mean_source_current_a) = if self.integration_steps == 1 {
            (last_voltage, last_current)
        } else {
            (
                voltage_integral.map(|v| (v / Dd::from_f64(self.period_s)).to_f64()),
                source_charge.map(|q| (q / Dd::from_f64(self.period_s)).to_f64()),
            )
        };
        if !scale.is_finite()
            || !allowance.is_finite()
            || winding_heat_j < 0.0
            || state
                .iter()
                .chain(&voltage)
                .chain(&mean_source_current_a)
                .chain(&source_work_j)
                .chain(&terminal_work_j)
                .any(|x| !x.is_finite())
            || channel_heat_j
                .iter()
                .chain(&diode_heat_j)
                .chain(&terminal_resistor_heat_j)
                .chain(&load_work_j)
                .chain(&curtailed_load_energy_j)
                .any(|x| !x.is_finite() || *x < 0.0)
            || residuals
                .iter()
                .any(|r| !r.is_finite() || r.abs() > allowance)
        {
            return Err(
                invalid("complete commutation cycle physical work/heat balance refused").into(),
            );
        }
        let (last_devices, last_loads) =
            last_points.ok_or_else(|| invalid("empty accepted commutation cycle"))?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&state);
        discrete_mean_bus_voltage_v.copy_from_slice(&voltage);
        Ok(CommutatedStep {
            native: PreparedStepRecord {
                delta_h,
                dissipated,
                supplied,
                newton_iters: iterations,
                solver_residual,
                solver_certificate: None,
            },
            terminal_work_j,
            winding_heat_j,
            channel_heat_j,
            diode_heat_j,
            source_work_j,
            terminal_resistor_heat_j,
            load_work_j,
            curtailed_load_energy_j,
            mean_source_current_a,
            last_devices,
            last_loads,
            integration_steps: self.integration_steps,
        })
    }
}
