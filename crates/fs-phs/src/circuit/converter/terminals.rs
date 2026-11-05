//! Voltage sources and resistive loads on the original converter trajectory.
//!
//! Each bus is connected through a finite resistance to a held ideal voltage.
//! The physical terminal current is `(source_voltage - bus_voltage) / R`.
//! A zero source voltage is an ordinary resistor load. A bus above its source
//! regenerates into that source; no current or power sign is clipped.
//!
//! The existing simultaneous scattering closure propagates the original three
//! storage coordinates. Its reference impedance is this actual terminal
//! resistor, and its incident wave is half the source voltage. The converter
//! heat remains distinct from terminal heat, which independently integrates
//! `(source_voltage - bus_voltage(t))² / R` on the SAME loaded trajectory.
//! A Norton-equivalent shunt loss `bus_voltage² / R` is not this resistor heat.

use super::{
    AveragedBuckBoost, CONDUCTION_CHANNELS, ConductionHeat, ConverterStep, invalid, poll, positive,
};
use crate::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedMultiportScatteringZoh, PreparedSourceLoss,
    ScatteringPortRecord,
};
use crate::{PortHamiltonian, PreparedStepError, PreparedStepRecord};
use fs_math::det;

/// Physical work and heat for the complete converter/source/load boundary.
#[derive(Clone, Debug)]
pub struct VoltageDrivenStep {
    /// Original converter storage, its five conduction paths, and actual work
    /// entering its two bus terminals. This excludes the terminal resistors.
    pub converter: ConverterStep,
    /// Signed work leaving each ideal voltage source [J]. Negative means the
    /// source receives regenerated energy. A passive zero-volt load has zero.
    pub source_work_j: [f64; 2],
    /// Independent positive heat in each source/load resistor [J]. These are
    /// additional physical owners, not part of `converter.conduction`.
    pub terminal_resistor_heat_j: [f64; 2],
}

/// Exact held-voltage propagation with both finite terminal resistances.
/// Preparation owns all factors; warm stepping uses fixed-size scratch and
/// introduces no allocation, exponential, or accepted physical history.
pub struct PreparedVoltageDrivenBuckBoost {
    native: PreparedMultiportScatteringZoh,
    resistance_ohm: [f64; 2],
    fractions: [f64; CONDUCTION_CHANNELS],
    terminal_heat: [PreparedSourceLoss; 2],
}

impl AveragedBuckBoost {
    /// Connect each original bus to an ideal voltage through a finite resistor.
    /// Both voltages are supplied to each warm call, held for the entire
    /// interval. Zero voltage models a load; negative voltage and regeneration
    /// retain their physical signs. Original flux and bus charge are unchanged.
    ///
    /// The same native coupled generator advances both bus responses together.
    /// It retains the original internal heat. Positive source-square factors
    /// separately integrate each terminal resistor on that trajectory, including
    /// transient voltage ripple. Squared mean current never supplies that heat.
    ///
    /// The dense budget applies to each of the two fixed three-state heat
    /// observers. Source/device provenance, temperature dependence, ideal
    /// zero-resistance voltage clamps, controllers and switching transitions
    /// remain outside this held linear model.
    ///
    /// # Errors
    /// Nonpositive/unrepresentable resistance or duration, original native
    /// admission, heat-observer work/storage limits, or cancellation.
    pub fn prepare_voltage_driven_controlled(
        self,
        terminal_resistance_ohm: [f64; 2],
        duration_s: f64,
        budget: LinearSplitBudget,
        heat_budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedVoltageDrivenBuckBoost, PreparedStepError> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || terminal_resistance_ohm
                .iter()
                .any(|&r| !positive(r) || !positive(1.0 / r))
        {
            return Err(invalid("voltage-driven converter requires positive representable terminal resistance and duration").into());
        }
        let total: f64 = self.effective_resistance_ohm.iter().sum();
        let fractions = self
            .effective_resistance_ohm
            .map(|r| if total > 0.0 { r / total } else { 0.0 });
        let native = PreparedMultiportScatteringZoh::new_owned(
            self.system,
            &[
                (0, terminal_resistance_ohm[0]),
                (1, terminal_resistance_ohm[1]),
            ],
            duration_s,
            budget,
        )?;
        poll(&mut cancelled)?;

        // Reconstruct the exact SAME loaded generator used by the native
        // scattering owner, with input a = source_voltage/2. The original
        // storage and G columns supply every coordinate map; no independent
        // electrical ODE or approximate terminal solve is introduced.
        let system = native.system();
        let (j, r, g) = system.structure();
        let q = system.storage().quadratic_part().ok_or_else(|| {
            invalid("voltage-driven converter lost its original quadratic storage")
        })?;
        let mut loaded_r = [0.0; 9];
        loaded_r.copy_from_slice(r);
        let mut input = [0.0; 6];
        let root = terminal_resistance_ohm.map(|r| det::sqrt(1.0 / r));
        for port in 0..2 {
            let inverse = 1.0 / terminal_resistance_ohm[port];
            for i in 0..3 {
                let root_i = root[port] * g[i * 2 + port];
                for k in 0..3 {
                    loaded_r[i * 3 + k] += root_i * (root[port] * g[k * 2 + port]);
                }
                input[i * 2 + port] = 2.0 * inverse * g[i * 2 + port];
            }
        }
        let flow: [f64; 9] = core::array::from_fn(|slot| {
            let (i, k) = (slot / 3, slot % 3);
            (0..3)
                .map(|a| (j[i * 3 + a] - loaded_r[i * 3 + a]) * q.matrix[a * 3 + k])
                .sum()
        });
        let channel = |port: usize, cancelled: &mut _| {
            let source: [f64; 3] = core::array::from_fn(|k| {
                -root[port]
                    * (0..3)
                        .map(|i| g[i * 2 + port] * q.matrix[i * 3 + k])
                        .sum::<f64>()
            });
            let mut direct = [0.0; 2];
            direct[port] = 2.0 * root[port];
            PreparedSourceLoss::new_with_feedthrough_controlled(
                3,
                2,
                &flow,
                &input,
                &source,
                &direct,
                duration_s,
                heat_budget,
                cancelled,
            )
        };
        let terminal_heat = [channel(0, &mut cancelled)?, channel(1, &mut cancelled)?];
        poll(&mut cancelled)?;
        Ok(PreparedVoltageDrivenBuckBoost {
            native,
            resistance_ohm: terminal_resistance_ohm,
            fractions,
            terminal_heat,
        })
    }
}

impl PreparedVoltageDrivenBuckBoost {
    /// Original converter, with no source/load resistance added to its losses.
    pub fn system(&self) -> &PortHamiltonian {
        self.native.system()
    }
    pub fn duration_s(&self) -> f64 {
        self.native.duration()
    }
    pub const fn terminal_resistance_ohm(&self) -> [f64; 2] {
        self.resistance_ohm
    }

    /// Advance the complete original state under the two held source voltages.
    /// Returns actual mean bus voltages and independent physical work/heat.
    /// Every channel and both external source balances must close before any
    /// caller output publishes. Refusal/cancellation is safe to retry.
    pub fn step_controlled(
        &mut self,
        initial: &[f64; 3],
        source_voltage_v: &[f64; 2],
        next: &mut [f64; 3],
        average_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VoltageDrivenStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let waves = source_voltage_v.map(|v| 0.5 * v);
        if source_voltage_v.iter().any(|v| !v.is_finite())
            || (0..2).any(|p| 2.0 * waves[p] != source_voltage_v[p])
        {
            return Err(invalid(
                "converter source voltage is not representable in its native boundary chart",
            )
            .into());
        }
        let (mut candidate, mut voltage) = ([0.0; 3], [0.0; 2]);
        let mut ports = [ScatteringPortRecord::default(); 2];
        let account =
            self.native
                .step_into(initial, &waves, &mut candidate, &mut voltage, &mut ports)?;
        let mut terminal_resistor_heat_j = [0.0; 2];
        for (heat, observer) in terminal_resistor_heat_j.iter_mut().zip(&self.terminal_heat) {
            *heat = observer.evaluate_into_controlled(
                initial,
                &waves,
                &mut [0.0; 6],
                &mut cancelled,
            )?;
        }
        let duration = self.duration_s();
        let source_work_j: [f64; 2] = core::array::from_fn(|p| {
            let mean_current = (source_voltage_v[p] - voltage[p]) / self.resistance_ohm[p];
            (duration * source_voltage_v[p]) * mean_current
        });
        let terminal_work_j = ports.map(|p| p.boundary_work);
        let joules = self.fractions.map(|f| f * account.physical_heat);
        let scale = self
            .system()
            .hamiltonian(initial)
            .abs()
            .max(self.system().hamiltonian(&candidate).abs())
            .max(account.delta_h.abs())
            .max(account.physical_heat.abs())
            .max(source_work_j.iter().map(|w| w.abs()).sum::<f64>())
            .max(terminal_resistor_heat_j.iter().sum::<f64>())
            .max(terminal_work_j.iter().map(|w| w.abs()).sum::<f64>());
        let residuals = [
            source_work_j[0] - terminal_resistor_heat_j[0] - terminal_work_j[0],
            source_work_j[1] - terminal_resistor_heat_j[1] - terminal_work_j[1],
            joules.iter().sum::<f64>() - account.physical_heat,
            account.delta_h + account.physical_heat + terminal_resistor_heat_j.iter().sum::<f64>()
                - source_work_j.iter().sum::<f64>(),
        ];
        if !scale.is_finite()
            || source_work_j
                .iter()
                .chain(&terminal_resistor_heat_j)
                .chain(&joules)
                .any(|v| !v.is_finite())
            || residuals
                .iter()
                .any(|v| !v.is_finite() || v.abs() > 1e-12 * scale)
        {
            return Err(invalid(
                "voltage-driven converter physical source/load work and heat do not close",
            )
            .into());
        }
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        average_bus_voltage_v.copy_from_slice(&voltage);
        Ok(VoltageDrivenStep {
            converter: ConverterStep {
                native: PreparedStepRecord {
                    delta_h: account.delta_h,
                    dissipated: account.physical_heat,
                    supplied: account.boundary_work,
                    newton_iters: 0,
                    solver_residual: 0.0,
                    solver_certificate: None,
                },
                terminal_work_j,
                conduction: ConductionHeat { joules },
            },
            source_work_j,
            terminal_resistor_heat_j,
        })
    }
}
