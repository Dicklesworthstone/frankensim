//! Prescribed PWM with nonlinear regulated demand and finite source impedance.
//!
//! Every authored gate interval uses the original three-state converter and
//! regulated-load Newton owner. Explicit subdivisions resolve its evolving
//! source current, load branch and winding ripple; all original flux/charge
//! coordinates continue across numerical steps and physical gate edges.
//!
//! Work and heat are sums of native discrete-effort quadrature, not exact
//! continuous integrals. Refine subdivisions to assess temporal accuracy. Gates
//! remain ideal complementary commands: no diode, dead-time, transition/recovery
//! heat, gate capacitance, controller, source capacity or device data is inferred.

use super::{
    AveragedBuckBoost, BridgeModulation, BridgeSwitchState, BuckBoostParameters,
    CONDUCTION_CHANNELS, ConductionHeat, ConverterStep, Sum, SwitchCycleControls, SwitchInterval,
    admit_cycle, invalid, poll, positive,
};
use crate::circuit::converter::power_load::{
    PowerLoadControls, PowerLoadPoint, PreparedPowerLoadedBuckBoost, RegulatedPowerLoad,
};
use crate::{PreparedStepError, PreparedStepRecord, storage_energy_change};
use fs_math::dd::Dd;

/// Explicit numerical work cap within EACH authored physical gate interval.
pub const MAX_POWER_SUBSTEPS_PER_INTERVAL: usize = 256;

/// A held voltage/demand command for one ORIGINAL physical gate interval.
/// Identical gate/duration caches do not merge or average these commands.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerSwitchCommand {
    pub source_voltage_v: [f64; 2],
    pub loads: [Option<RegulatedPowerLoad>; 2],
}

/// Complete-cycle native storage, separately owned work and physical heat.
#[derive(Clone, Debug)]
pub struct PowerLoadedSwitchCycleStep {
    /// Original endpoint storage change, NET bus work and five converter heats.
    /// Native Newton counts are summed and residual is the maximum accepted
    /// substep residual; no single algebraic certificate represents the cycle.
    pub converter: ConverterStep,
    /// Signed work leaving each held ideal voltage source [J].
    pub source_work_j: [f64; 2],
    /// Physical source-resistor heat, separate from converter heat [J].
    pub terminal_resistor_heat_j: [f64; 2],
    /// Actual delivered regulated-load energy, external by default [J].
    pub load_work_j: [f64; 2],
    pub curtailed_load_energy_j: [f64; 2],
    /// Duration-weighted native discrete source currents, positive into buses.
    pub mean_source_current_a: [f64; 2],
    /// Load points at the FINAL accepted discrete substep effort only.
    /// They are not cycle-average branch labels or endpoint observations.
    pub last_loads: [PowerLoadPoint; 2],
    pub integration_steps: usize,
}

struct CachedStep {
    state: BridgeSwitchState,
    duration_s: f64,
    native: PreparedPowerLoadedBuckBoost,
}

#[derive(Clone, Copy)]
struct IntervalSteps {
    regular: usize,
    last: usize,
}

/// Reusable explicitly subdivided physical gate cycle with nonlinear loads.
///
/// Parameters and terminal resistances are frozen over a complete cycle. Each
/// distinct binary connection/substep-duration pair shares one native workspace.
/// Held commands may change on every authored gate interval or between calls.
/// Preparation allocates; successful warm cycles allocate no new storage.
pub struct PreparedPowerLoadedSwitchCycle {
    parameters: BuckBoostParameters,
    terminal_resistance_ohm: [f64; 2],
    intervals: Vec<SwitchInterval>,
    interval_steps: Vec<IntervalSteps>,
    caches: Vec<CachedStep>,
    substeps_per_interval: usize,
    integration_steps: usize,
    period_s: f64,
    averaged_modulation: BridgeModulation,
    controls: SwitchCycleControls,
}

fn total(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |sum, value| sum + Dd::from_f64(value))
        .to_f64()
}

impl PreparedPowerLoadedSwitchCycle {
    /// Prepare the SAME implicit converter owner for every distinct substep.
    /// The existing at-most-64 physical-interval cap and an explicit 1..=256
    /// subdivision count bound work before native preparation. The caller's
    /// Newton/backtracking budget applies to each numerical integration step.
    /// Their complete-cycle product must be representable.
    ///
    /// The last substep is the compensated remainder of its authored interval.
    /// All substeps must be positive and their compensated duration must equal
    /// that interval exactly in the represented clock; otherwise preparation
    /// refuses. A numerical subdivision never creates an additional gate edge.
    #[allow(clippy::too_many_arguments)]
    pub fn new_controlled(
        parameters: BuckBoostParameters,
        terminal_resistance_ohm: [f64; 2],
        intervals: &[SwitchInterval],
        substeps_per_interval: usize,
        electrical: PowerLoadControls,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        let (period_s, averaged_modulation) = admit_cycle(intervals, controls, &mut cancelled)?;
        if !(1..=MAX_POWER_SUBSTEPS_PER_INTERVAL).contains(&substeps_per_interval)
            || !(0.0..1.0).contains(&controls.relative_energy_tolerance)
            || terminal_resistance_ohm
                .iter()
                .any(|r| !positive(*r) || !positive(1.0 / r))
        {
            return Err(invalid(
                "regulated switching requires bounded subdivisions, finite positive source resistance and relative allowance below one",
            ).into());
        }
        let integration_steps = intervals
            .len()
            .checked_mul(substeps_per_interval)
            .ok_or_else(|| invalid("regulated switching step count is unrepresentable"))?;
        integration_steps
            .checked_mul(electrical.solve.max_iterations)
            .ok_or_else(|| invalid("regulated switching Newton work bound is unrepresentable"))?;
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
                    invalid("regulated switch subdivision duration is unrepresentable").into(),
                );
            }
            let mut represented = Sum::default();
            for k in 0..substeps_per_interval {
                poll(&mut cancelled)?;
                represented.add(if k + 1 == substeps_per_interval {
                    last_s
                } else {
                    regular_s
                });
            }
            if represented.total().to_bits() != interval.duration_s.to_bits() {
                return Err(
                    invalid("regulated subdivisions must preserve the authored gate edge").into(),
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
                    let native = AveragedBuckBoost::new(parameters, interval.state.modulation())?
                        .prepare_power_loaded_controlled(
                        duration_s,
                        electrical,
                        &mut cancelled,
                    )?;
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
            parameters,
            terminal_resistance_ohm,
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

    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    pub const fn terminal_resistance_ohm(&self) -> [f64; 2] {
        self.terminal_resistance_ohm
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
    /// Number of native solves charged for one complete successful cycle.
    pub const fn integration_steps(&self) -> usize {
        self.integration_steps
    }
    pub fn cached_interval_count(&self) -> usize {
        self.caches.len()
    }

    /// Advance with the same voltage/demand command on all authored intervals.
    /// Output voltages are numerical duration-weighted discrete efforts. Both
    /// caller buffers remain untouched on any interval or whole-cycle refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        source_voltage_v: &[f64; 2],
        loads: &[Option<RegulatedPowerLoad>; 2],
        next: &mut [f64; 3],
        discrete_mean_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerLoadedSwitchCycleStep, PreparedStepError> {
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

    /// Apply one held command per ORIGINAL gate interval, independent of cache
    /// reuse and numerical subdivisions. Source work retains source/current
    /// correlation; delivered load work is never requested power times period
    /// when foldback or current limiting curtailed any accepted substep.
    pub fn step_interval_commands_controlled(
        &mut self,
        initial: &[f64; 3],
        commands: &[PowerSwitchCommand],
        next: &mut [f64; 3],
        discrete_mean_bus_voltage_v: &mut [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<PowerLoadedSwitchCycleStep, PreparedStepError> {
        if commands.len() != self.intervals.len() {
            return Err(
                invalid("one voltage/load command is required per original gate interval").into(),
            );
        }
        self.advance(
            initial,
            |interval| commands[interval],
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
    ) -> Result<PowerLoadedSwitchCycleStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let mut heat = [Sum::default(); CONDUCTION_CHANNELS];
        let mut terminal_heat = [Sum::default(); 2];
        let mut source_work = [Sum::default(); 2];
        let mut terminal_work = [Sum::default(); 2];
        let mut load_work = [Sum::default(); 2];
        let mut curtailed = [Sum::default(); 2];
        let mut voltage_integral = [Sum::default(); 2];
        let mut source_charge = [Sum::default(); 2];
        let mut port_transfers = [Sum::default(); 2];
        let mut native_change = Sum::default();
        let mut native_heat = Sum::default();
        let mut native_work = Sum::default();
        let mut absolute_change = Sum::default();
        let mut newton_iters = 0usize;
        let mut solver_residual = 0.0_f64;
        let (mut last_voltage, mut last_current, mut last_loads) = ([0.0; 2], [0.0; 2], None);
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
                let step = cache.native.step_voltage_driven_controlled(
                    &state,
                    &command.source_voltage_v,
                    &self.terminal_resistance_ohm,
                    &command.loads,
                    &mut endpoint,
                    &mut voltage,
                    &mut cancelled,
                )?;
                for channel in 0..CONDUCTION_CHANNELS {
                    heat[channel].add(step.converter.conduction.joules[channel]);
                }
                for port in 0..2 {
                    terminal_heat[port].add(step.terminal_resistor_heat_j[port]);
                    source_work[port].add(step.source_work_j[port]);
                    terminal_work[port].add(step.converter.terminal_work_j[port]);
                    load_work[port].add(step.load_work_j[port]);
                    curtailed[port].add(step.curtailed_load_energy_j[port]);
                    voltage_integral[port].add(cache.duration_s * voltage[port]);
                    source_charge[port].add(cache.duration_s * step.source_current_a[port]);
                    port_transfers[port].add(step.source_work_j[port].abs());
                    port_transfers[port].add(step.converter.terminal_work_j[port].abs());
                    port_transfers[port].add(step.load_work_j[port]);
                }
                native_change.add(step.converter.native.delta_h);
                absolute_change.add(step.converter.native.delta_h.abs());
                native_heat.add(step.converter.native.dissipated);
                native_work.add(step.converter.native.supplied);
                newton_iters = newton_iters
                    .checked_add(step.converter.native.newton_iters)
                    .ok_or_else(|| invalid("regulated switching Newton count overflow"))?;
                solver_residual = solver_residual.max(step.converter.native.solver_residual);
                state = endpoint;
                last_voltage = voltage;
                last_current = step.source_current_a;
                last_loads = Some(step.loads);
            }
        }
        let delta_h = storage_energy_change(
            self.caches[0].native.system().storage(),
            initial,
            &state,
            &mut [0.0; 3],
            &mut [0.0; 3],
        );
        let joules = heat.map(Sum::total);
        let terminal_resistor_heat_j = terminal_heat.map(Sum::total);
        let source_work_j = source_work.map(Sum::total);
        let terminal_work_j = terminal_work.map(Sum::total);
        let load_work_j = load_work.map(Sum::total);
        let curtailed_load_energy_j = curtailed.map(Sum::total);
        let transfers = port_transfers.map(Sum::total);
        let (voltage, mean_source_current_a) = if self.integration_steps == 1 {
            (last_voltage, last_current)
        } else {
            (
                voltage_integral.map(|value| value.total() / self.period_s),
                source_charge.map(|value| value.total() / self.period_s),
            )
        };
        let dissipated = native_heat.total();
        let supplied = native_work.total();
        let terminal_total = total(terminal_resistor_heat_j);
        let absolute_change = absolute_change.total();
        let scale = total(transfers)
            .max(absolute_change)
            .max(dissipated.abs())
            .max(terminal_total)
            .max(delta_h.abs());
        let tolerance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        let residuals = [
            total([delta_h, dissipated, -supplied]),
            total(
                [delta_h, dissipated]
                    .into_iter()
                    .chain(terminal_resistor_heat_j)
                    .chain(load_work_j)
                    .chain(source_work_j.map(|work| -work)),
            ),
            total([delta_h, -native_change.total()]),
            total(joules.into_iter().chain([-dissipated])),
            total(terminal_work_j.into_iter().chain([-supplied])),
        ];
        for port in 0..2 {
            let port_scale = transfers[port].max(terminal_resistor_heat_j[port]);
            let allowance = self.controls.absolute_energy_tolerance_j
                + self.controls.relative_energy_tolerance * port_scale;
            let residual = total([
                source_work_j[port],
                -terminal_resistor_heat_j[port],
                -load_work_j[port],
                -terminal_work_j[port],
            ]);
            if !allowance.is_finite() || !residual.is_finite() || residual.abs() > allowance {
                return Err(
                    invalid("regulated switch source/resistor/load work does not close").into(),
                );
            }
        }
        if !scale.is_finite()
            || !tolerance.is_finite()
            || !absolute_change.is_finite()
            || !solver_residual.is_finite()
            || state
                .iter()
                .chain(&voltage)
                .chain(&mean_source_current_a)
                .chain(&transfers)
                .chain(&joules)
                .chain(&terminal_resistor_heat_j)
                .chain(&source_work_j)
                .chain(&terminal_work_j)
                .chain(&load_work_j)
                .chain(&curtailed_load_energy_j)
                .any(|value| !value.is_finite())
            || joules
                .iter()
                .chain(&terminal_resistor_heat_j)
                .chain(&load_work_j)
                .chain(&curtailed_load_energy_j)
                .any(|value| *value < 0.0)
            || residuals
                .iter()
                .any(|value| !value.is_finite() || value.abs() > tolerance)
        {
            return Err(invalid(
                "complete regulated switching cycle physical work/heat balance refused",
            )
            .into());
        }
        let last_loads =
            last_loads.ok_or_else(|| invalid("regulated switch cycle has no accepted substep"))?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&state);
        discrete_mean_bus_voltage_v.copy_from_slice(&voltage);
        Ok(PowerLoadedSwitchCycleStep {
            converter: ConverterStep {
                native: PreparedStepRecord {
                    delta_h,
                    dissipated,
                    supplied,
                    newton_iters,
                    solver_residual,
                    solver_certificate: None,
                },
                terminal_work_j,
                conduction: ConductionHeat { joules },
            },
            source_work_j,
            terminal_resistor_heat_j,
            load_work_j,
            curtailed_load_energy_j,
            mean_source_current_a,
            last_loads,
            integration_steps: self.integration_steps,
        })
    }
}
