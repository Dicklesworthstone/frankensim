//! Regulated loads and saturated magnetic storage in one original native step.
//!
//! The existing bipolar current-limit/foldback law supplies both its resisting
//! current and exact branch derivative. Its negative incremental conductance in
//! constant-power operation enters the same Newton equation as the source-owned
//! winding discrete-current derivative. No lagged load current or separate bus
//! integrator is introduced. Exported load work remains external work, not heat
//! in any of the converter's five conduction owners.

use super::{
    NonlinearConverterStep, PreparedNonlinearAveragedBuckBoost, invalid, poll, product, sum,
};
use crate::circuit::converter::power_load::{PowerLoadPoint, RegulatedPowerLoad, load_point};
use crate::{NonlinearPortPower, PreparedStepError};

/// All original converter accounts plus separately owned external load work.
#[derive(Clone, Debug)]
pub struct NonlinearPowerLoadedStep {
    /// Terminal work and native supply are NET source work minus exported work.
    /// Native dissipation and the five channels contain only converter heat.
    pub converter: NonlinearConverterStep,
    /// Authored load branch at each accepted discrete bus voltage.
    pub loads: [PowerLoadPoint; 2],
    /// Signed held-current source work entering input/output buses [J].
    pub source_work_j: [f64; 2],
    /// Positive energy exported to each external load [J], never converter heat.
    pub load_work_j: [f64; 2],
    /// Requested minus delivered load energy at native discrete quadrature [J].
    pub curtailed_load_energy_j: [f64; 2],
}

impl PreparedNonlinearAveragedBuckBoost {
    /// Solve saturated winding, finite buses and both regulated demands jointly.
    /// Sources are held currents positive INTO the buses; each optional load
    /// draws nonnegative power at either voltage polarity. Commands may change
    /// between calls while the original caller-owned flux and charges continue.
    ///
    /// Initial, final and discrete load voltages must lie in the authored domain;
    /// these checks do not certify unseen continuous-time overshoot. Work/heat
    /// retain native discrete-effort quadrature and all existing source limits.
    /// Nonlinear convergence does not establish stability or root uniqueness.
    /// No callback or account allocation is added. Every refusal/cancellation
    /// preserves BOTH output buffers, with no accepted history in this cache.
    #[allow(clippy::too_many_arguments)]
    pub fn step_power_loaded_controlled(
        &mut self,
        initial: &[f64; 3],
        terminal_current_a: &[f64; 2],
        loads: &[Option<RegulatedPowerLoad>; 2],
        next: &mut [f64; 3],
        discrete_bus_voltage_v: &mut [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearPowerLoadedStep, PreparedStepError> {
        poll(&mut cancelled)?;
        let inverse_c = self.converter.bus_capacitance_f.map(|c| 1.0 / c);
        for p in 0..2 {
            load_point(loads[p], initial[p + 1] * inverse_c[p])?;
        }
        let port = |_: &[f64], effort: &[f64], out: &mut [f64]| {
            if effort.len() != 3 || out.len() != 3 {
                return None;
            }
            out.fill(0.0);
            let mut power = [0.0; 2];
            for p in 0..2 {
                let point = load_point(loads[p], effort[p + 1]).ok()?;
                out[p + 1] = point.current_a;
                power[p] = point.delivered_power_w;
            }
            Some(NonlinearPortPower {
                supplied_power_w: -sum(power),
                dissipated_power_w: 0.0,
            })
        };
        let tangent = |_: &[f64], effort: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            if effort.len() != 3 || de.len() != 3 || out.len() != 3 {
                return false;
            }
            out.fill(0.0);
            for p in 0..2 {
                let Ok(point) = load_point(loads[p], effort[p + 1]) else {
                    return false;
                };
                out[p + 1] = point.incremental_conductance_s * de[p + 1];
            }
            out.iter().all(|v| v.is_finite())
        };
        let (mut candidate, mut voltage) = ([0.0; 3], [0.0; 2]);
        let native = self.workspace.step_into_supplied_controlled(
            &self.converter.system,
            initial,
            terminal_current_a,
            self.duration_s,
            &mut candidate,
            &mut voltage,
            &port,
            None, // Use the original winding's source-owned discrete derivative.
            Some(&tangent),
            &mut cancelled,
        )?;
        for p in 0..2 {
            load_point(loads[p], candidate[p + 1] * inverse_c[p])?;
        }
        let points = [
            load_point(loads[0], voltage[0])?,
            load_point(loads[1], voltage[1])?,
        ];
        let h = self.duration_s;
        let load_work_j = points.map(|p| product([h, p.delivered_power_w]));
        let curtailed_load_energy_j = points.map(|p| product([h, p.curtailed_power_w]));
        let source_work_j: [f64; 2] =
            core::array::from_fn(|p| product([h, terminal_current_a[p], voltage[p]]));
        let terminal_work_j = core::array::from_fn(|p| sum([source_work_j[p], -load_work_j[p]]));
        let work_scale_j = core::array::from_fn(|p| sum([source_work_j[p].abs(), load_work_j[p]]));
        if source_work_j.iter().any(|v| !v.is_finite())
            || load_work_j
                .iter()
                .chain(&curtailed_load_energy_j)
                .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "nonlinear converter source/load/curtailment work is unrepresentable",
            )
            .into());
        }
        let converter = self.account_interval(
            initial,
            &candidate,
            &voltage,
            terminal_work_j,
            work_scale_j,
            native,
        )?;
        poll(&mut cancelled)?;
        next.copy_from_slice(&candidate);
        discrete_bus_voltage_v.copy_from_slice(&voltage);
        Ok(NonlinearPowerLoadedStep {
            converter,
            loads: points,
            source_work_j,
            load_work_j,
            curtailed_load_energy_j,
        })
    }
}
