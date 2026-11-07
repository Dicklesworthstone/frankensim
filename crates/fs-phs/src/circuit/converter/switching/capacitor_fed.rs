//! PWM on a converter with finite external capacitor supplies.
//!
//! Every edge retains the original inductor flux, both bus charges and both
//! external capacitor charges. Terminal voltages evolve inside each interval,
//! so source sag, reverse charging and resistor ripple heat participate in the
//! same electrical trajectory. All gates and constitutive data remain the
//! prescribed ideal complementary, frozen-linear model of the parent owner.

use super::{
    AveragedBuckBoost, BridgeModulation, BuckBoostParameters, CONDUCTION_CHANNELS, ConductionHeat,
    ConverterStep, Sum, SwitchCycleControls, SwitchInterval, admit_cycle, invalid, poll,
};
use crate::circuit::converter::capacitor_terminals::{
    CapacitorFedStep, CapacitorTerminals, PreparedCapacitorFedBuckBoost,
};
use crate::linear_exact::{DenseExponentialBudget, LinearSplitBudget};
use crate::{PreparedStepError, PreparedStepRecord, storage_energy_change};
use fs_math::dd::Dd;

struct CachedInterval {
    interval: SwitchInterval,
    native: PreparedCapacitorFedBuckBoost,
}

/// Reusable gate cycle with five continuing storage states. Identical gate/
/// duration pairs share an exact propagator and its independent heat observers.
/// Warm stepping is O(intervals) fixed-size work, without new allocations,
/// factorization or exponentials. This is not a measured throughput claim.
pub struct PreparedCapacitorFedSwitchCycle {
    parameters: BuckBoostParameters,
    terminals: CapacitorTerminals,
    intervals: Vec<SwitchInterval>,
    cache_indices: Vec<usize>,
    caches: Vec<CachedInterval>,
    period_s: f64,
    averaged_modulation: BridgeModulation,
    controls: SwitchCycleControls,
}

impl PreparedCapacitorFedSwitchCycle {
    /// Admit all physical intervals before allocation. Budgets apply per unique
    /// cache; the existing interval cap bounds cache count and complete warm work.
    /// Capacitor values are fixed, so changing gate state performs no parameter
    /// work and requires no projection. Switching-transition heat is not modeled.
    pub fn new_controlled(
        parameters: BuckBoostParameters,
        terminals: CapacitorTerminals,
        intervals: &[SwitchInterval],
        budget: LinearSplitBudget,
        observation_budget: DenseExponentialBudget,
        controls: SwitchCycleControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PreparedStepError> {
        let (period_s, averaged_modulation) = admit_cycle(intervals, controls, &mut cancelled)?;
        let mut caches: Vec<CachedInterval> = Vec::new();
        let mut cache_indices = Vec::with_capacity(intervals.len());
        for &interval in intervals {
            poll(&mut cancelled)?;
            let cached = caches.iter().position(|cache| {
                cache.interval.state == interval.state
                    && cache.interval.duration_s.to_bits() == interval.duration_s.to_bits()
            });
            let index = match cached {
                Some(index) => index,
                None => {
                    let native = AveragedBuckBoost::new(parameters, interval.state.modulation())?
                        .with_capacitor_terminals(terminals)?
                        .prepare_controlled(
                            interval.duration_s,
                            budget,
                            observation_budget,
                            &mut cancelled,
                        )?;
                    let index = caches.len();
                    caches.push(CachedInterval { interval, native });
                    index
                }
            };
            cache_indices.push(index);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            parameters,
            terminals,
            intervals: intervals.to_vec(),
            cache_indices,
            caches,
            period_s,
            averaged_modulation,
            controls,
        })
    }

    pub const fn parameters(&self) -> BuckBoostParameters {
        self.parameters
    }
    pub const fn terminals(&self) -> CapacitorTerminals {
        self.terminals
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
    pub fn cached_interval_count(&self) -> usize {
        self.caches.len()
    }

    /// Held external charging currents for the complete cycle. The physical
    /// terminal currents are evolving state responses, not these source commands.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 5],
        external_current_a: &[f64; 2],
        next: &mut [f64; 5],
        average_voltage_v: &mut [f64; 4],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitorFedStep, PreparedStepError> {
        self.advance(
            initial,
            |_| *external_current_a,
            next,
            average_voltage_v,
            cancelled,
        )
    }

    /// Independently held charging currents at each authored gate interval.
    /// Their work retains the source/gate correlation instead of multiplying
    /// an averaged current by an independently averaged voltage.
    pub fn step_interval_sources_controlled(
        &mut self,
        initial: &[f64; 5],
        external_current_a: &[[f64; 2]],
        next: &mut [f64; 5],
        average_voltage_v: &mut [f64; 4],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitorFedStep, PreparedStepError> {
        if external_current_a.len() != self.intervals.len() {
            return Err(invalid(
                "one external current pair is required per capacitor-fed interval",
            )
            .into());
        }
        self.advance(
            initial,
            |k| external_current_a[k],
            next,
            average_voltage_v,
            cancelled,
        )
    }

    fn advance(
        &mut self,
        initial: &[f64; 5],
        source: impl Fn(usize) -> [f64; 2],
        next: &mut [f64; 5],
        average_voltage_v: &mut [f64; 4],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitorFedStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let mut state = *initial;
        let mut heat = [Sum::default(); CONDUCTION_CHANNELS];
        let mut terminal_heat = [Sum::default(); 2];
        let mut capacitor_work = [Sum::default(); 2];
        let mut capacitor_change = [Sum::default(); 2];
        let mut bus_work = [Sum::default(); 2];
        let mut external_work = [Sum::default(); 2];
        let mut voltage_integral = [Sum::default(); 4];
        let (mut original_change, mut complete_change, mut absolute_transfers) =
            (Sum::default(), Sum::default(), Sum::default());
        let mut last_voltage = [0.0; 4];
        for k in 0..self.intervals.len() {
            poll(&mut cancelled)?;
            let (mut endpoint, mut voltage) = ([0.0; 5], [0.0; 4]);
            let record = self.caches[self.cache_indices[k]].native.step_controlled(
                &state,
                &source(k),
                &mut endpoint,
                &mut voltage,
                &mut cancelled,
            )?;
            for channel in 0..CONDUCTION_CHANNELS {
                heat[channel].add(record.converter.conduction.joules[channel]);
            }
            for p in 0..2 {
                terminal_heat[p].add(record.terminal_resistor_heat_j[p]);
                capacitor_work[p].add(record.terminal_outflow_work_j[p]);
                capacitor_change[p].add(record.terminal_storage_change_j[p]);
                bus_work[p].add(record.converter.terminal_work_j[p]);
                external_work[p].add(record.external_work_j[p]);
                absolute_transfers.add(record.terminal_outflow_work_j[p].abs());
                absolute_transfers.add(record.converter.terminal_work_j[p].abs());
                absolute_transfers.add(record.external_work_j[p].abs());
                absolute_transfers.add(record.terminal_storage_change_j[p].abs());
            }
            for p in 0..4 {
                voltage_integral[p].add(self.intervals[k].duration_s * voltage[p]);
            }
            original_change.add(record.converter.native.delta_h);
            complete_change.add(record.total.delta_h);
            absolute_transfers.add(record.converter.native.delta_h.abs());
            absolute_transfers.add(record.total.dissipated);
            state = endpoint;
            last_voltage = voltage;
        }
        let system = &self.caches[0].native;
        let total_change = storage_energy_change(
            system.system().storage(),
            initial,
            &state,
            &mut [0.0; 5],
            &mut [0.0; 5],
        );
        let converter_change = storage_energy_change(
            system.converter_system().storage(),
            &initial[..3],
            &state[..3],
            &mut [0.0; 3],
            &mut [0.0; 3],
        );
        let terminal_storage_change_j = core::array::from_fn(|p| {
            let before = Dd::from_f64(initial[p + 3]);
            let after = Dd::from_f64(state[p + 3]);
            ((after - before) * (after * Dd::from_f64(0.5) + before * Dd::from_f64(0.5))
                / Dd::from_f64(self.terminals.capacitance_f[p]))
            .to_f64()
        });
        let joules = heat.map(Sum::total);
        let terminal_resistor_heat_j = terminal_heat.map(Sum::total);
        let terminal_outflow_work_j = capacitor_work.map(Sum::total);
        let terminal_work_j = bus_work.map(Sum::total);
        let external_work_j = external_work.map(Sum::total);
        let voltage = if self.intervals.len() == 1 {
            last_voltage
        } else {
            voltage_integral.map(|s| s.total() / self.period_s)
        };
        let sum = |values: &[f64]| {
            values
                .iter()
                .fold(Dd::ZERO, |s, &v| s + Dd::from_f64(v))
                .to_f64()
        };
        let converter_heat = sum(&joules);
        let total_heat = sum(&[converter_heat, sum(&terminal_resistor_heat_j)]);
        let supplied = sum(&external_work_j);
        let converter_supplied = sum(&terminal_work_j);
        let mut residuals = [0.0; 11];
        residuals[0] = sum(&[total_change, total_heat, -supplied]);
        residuals[1] = sum(&[converter_change, converter_heat, -converter_supplied]);
        residuals[2] = sum(&[total_change, -complete_change.total()]);
        residuals[3] = sum(&[converter_change, -original_change.total()]);
        residuals[4] = sum(&[
            converter_change,
            terminal_storage_change_j[0],
            terminal_storage_change_j[1],
            -total_change,
        ]);
        for p in 0..2 {
            residuals[5 + p] = sum(&[terminal_storage_change_j[p], -capacitor_change[p].total()]);
            residuals[7 + p] = sum(&[
                terminal_storage_change_j[p],
                terminal_outflow_work_j[p],
                -external_work_j[p],
            ]);
            residuals[9 + p] = sum(&[
                terminal_outflow_work_j[p],
                -terminal_work_j[p],
                -terminal_resistor_heat_j[p],
            ]);
        }
        let scale = absolute_transfers
            .total()
            .max(total_change.abs())
            .max(total_heat);
        let tolerance = self.controls.absolute_energy_tolerance_j
            + self.controls.relative_energy_tolerance * scale;
        if !scale.is_finite()
            || !tolerance.is_finite()
            || voltage
                .iter()
                .chain(&joules)
                .chain(&terminal_resistor_heat_j)
                .chain(&terminal_outflow_work_j)
                .chain(&terminal_work_j)
                .chain(&external_work_j)
                .chain(&terminal_storage_change_j)
                .any(|v| !v.is_finite())
            || joules
                .iter()
                .chain(&terminal_resistor_heat_j)
                .any(|v| *v < 0.0)
            || residuals
                .iter()
                .any(|v| !v.is_finite() || v.abs() > tolerance)
        {
            return Err(
                invalid("capacitor-fed cycle storage and physical work/heat do not close").into(),
            );
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&state);
        average_voltage_v.copy_from_slice(&voltage);
        let record = |delta_h, dissipated, supplied| PreparedStepRecord {
            delta_h,
            dissipated,
            supplied,
            newton_iters: 0,
            solver_residual: 0.0,
            solver_certificate: None,
        };
        Ok(CapacitorFedStep {
            converter: ConverterStep {
                native: record(converter_change, converter_heat, converter_supplied),
                terminal_work_j,
                conduction: ConductionHeat { joules },
            },
            total: record(total_change, total_heat, supplied),
            terminal_storage_change_j,
            terminal_outflow_work_j,
            external_work_j,
            terminal_resistor_heat_j,
        })
    }
}
