//! Physical-time thermal creep with the original Maxwell and Kelvin losses.

use super::{
    PreparedTensorKelvinStressControl, Result, TensorKelvinStressRampStep, TensorPronyHistory,
    Voigt, finite, invalid, mandel, poll, tensor,
};
use crate::visco::temperature::{TemperatureShift, TemperatureShiftState};
use crate::visco_tensor::temperature::{
    TemperatureRamp, TemperatureRampControls, TemperatureSchedule,
};
use fs_math::dd::Dd;
use std::sync::Arc;

/// Accepted prescribed-temperature interval with the original total viscous heat.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorKelvinStressTemperatureStep {
    /// Original constitutive response. The endpoint strain rate is per unit
    /// PHYSICAL time at the actual final temperature. Work and heat are sums of
    /// local increments, independent of the cumulative history's magnitude.
    pub response: TensorKelvinStressRampStep,
    /// Calibration evaluated at the actual initial temperature.
    pub start: TemperatureShiftState,
    /// Calibration evaluated at the actual final temperature.
    pub end: TemperatureShiftState,
    /// Exact constant-shift constitutive subintervals actually advanced.
    pub substeps: usize,
    /// Analytical subinterval ln(a_T) variation estimate, not a trajectory bound.
    pub maximum_log_shift_change: f64,
}

fn reduced_duration(physical: f64, scale: f64) -> Result<f64> {
    let reduced = physical / scale;
    if !reduced.is_finite() || (physical > 0.0 && reduced <= 0.0) {
        return Err(invalid("unrepresentable reduced-time Kelvin interval"));
    }
    Ok(reduced)
}

fn physical_endpoint_rate(response: &mut TensorKelvinStressRampStep, scale: f64) -> Result<()> {
    for rate in &mut response.strain_rate {
        let physical = *rate / scale;
        if !physical.is_finite() || (*rate != 0.0 && physical == 0.0) {
            return Err(invalid(
                "unrepresentable physical Kelvin endpoint strain rate",
            ));
        }
        *rate = physical;
    }
    Ok(())
}

impl PreparedTensorKelvinStressControl {
    /// Advance linear stress and temperature ramps in physical time.
    ///
    /// The reference law is declared at the calibration reference temperature.
    /// BOTH original Maxwell times and Kelvin viscosity scale by a_T. Each
    /// admitted midpoint-temperature subinterval therefore uses this SAME
    /// prepared controller for h/a_T: no factorization, refitting or memory
    /// translation occurs. Work and total dashpot heat need no rescaling.
    ///
    /// Stress interpolation remains linear in physical time, rather than in
    /// accumulated reduced time. The final strain rate uses the actual END
    /// shift, including when the final midpoint differs substantially from it.
    /// Constant reference temperature retains the original response bits.
    ///
    /// The original fixed moduli, branch tensors, strain bound and SPD Kelvin
    /// admission remain binding. This is a second-order coefficient-frozen
    /// approximation to a prescribed thermal path, with the existing explicit
    /// variation and work limits. Refine these limits and the loading history;
    /// no thermal expansion, caloric evolution or trajectory-error certificate
    /// is inferred. A temperature change requires positive physical duration.
    ///
    /// # Errors
    /// Original binding, source, strain and balance refusals, invalid thermal
    /// support or budgets, unrepresentable physical/reduced quantities, or
    /// cancellation. Every incoming history field remains unchanged on error,
    /// including endpoint-rate conversion and the final cancellation poll.
    #[allow(clippy::too_many_lines)] // One atomic physical-time constitutive transaction.
    pub fn advance_temperature_ramp(
        &self,
        history: &mut TensorPronyHistory,
        stress: Voigt,
        ramp: TemperatureRamp,
        calibration: &TemperatureShift,
        controls: TemperatureRampControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorKelvinStressTemperatureStep> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.model, &history.model) || !finite(&stress) {
            return Err(invalid(
                "thermal Kelvin control requires finite stress and original history",
            ));
        }
        let schedule =
            TemperatureSchedule::new(&self.model, ramp, calibration, controls, &mut cancelled)?;
        let mut candidate = history.clone();
        if ramp.start_k == ramp.end_k {
            let mut response = self.advance(
                &mut candidate,
                stress,
                reduced_duration(ramp.duration_s, schedule.start.factor())?,
                &mut cancelled,
            )?;
            physical_endpoint_rate(&mut response, schedule.end.factor())?;
            poll(&mut cancelled)?;
            *history = candidate;
            return Ok(TensorKelvinStressTemperatureStep {
                response,
                start: schedule.start,
                end: schedule.end,
                substeps: 1,
                maximum_log_shift_change: 0.0,
            });
        }
        let initial = history.stress;
        let target = mandel(&stress);
        let delta = core::array::from_fn::<_, 6, _>(|i| target[i] - initial[i]);
        if !finite(&target) || !finite(&delta) {
            return Err(invalid("unrepresentable physical Kelvin stress ramp"));
        }
        let mut work = Dd::ZERO;
        let mut heat = Dd::ZERO;
        let mut storage_change = Dd::ZERO;
        let mut increment = [Dd::ZERO; 6];
        let mut bounds = [history.strain(); 2];
        let mut elapsed = 0.0;
        let mut last = None;
        for step in 0..schedule.substeps {
            poll(&mut cancelled)?;
            let (endpoint, h, scale) = schedule.interval(step, elapsed, calibration)?;
            let next = if step + 1 == schedule.substeps {
                stress
            } else {
                tensor(&core::array::from_fn(|i| {
                    (endpoint / ramp.duration_s).mul_add(delta[i], initial[i])
                }))
            };
            let local = self.advance(
                &mut candidate,
                next,
                reduced_duration(h, scale)?,
                &mut cancelled,
            )?;
            work = work + Dd::from_f64(local.work);
            heat = heat + Dd::from_f64(local.dissipated);
            storage_change = storage_change + Dd::from_f64(local.stored_energy_change);
            for i in 0..6 {
                increment[i] = increment[i] + Dd::from_f64(local.strain_increment[i]);
                bounds[0][i] = bounds[0][i].min(local.strain_bounds[0][i]);
                bounds[1][i] = bounds[1][i].max(local.strain_bounds[1][i]);
            }
            last = Some(local);
            elapsed = endpoint;
        }
        let mut response = last.ok_or_else(|| invalid("missing thermal Kelvin interval"))?;
        response.strain_increment = increment.map(Dd::to_f64);
        response.stored_energy_change = storage_change.to_f64();
        response.work = work.to_f64();
        response.dissipated = heat.to_f64();
        response.balance_residual = (work - storage_change - heat).to_f64();
        response.strain_bounds = bounds;
        physical_endpoint_rate(&mut response, schedule.end.factor())?;
        let total_work = history.work + response.work;
        let total_heat = history.heat + response.dissipated;
        if !finite(&response.strain_increment)
            || response.dissipated < 0.0
            || ![
                response.stored_energy_change,
                response.work,
                response.dissipated,
                response.balance_residual,
                total_work,
                total_heat,
            ]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(invalid("unrepresentable thermal Kelvin work or heat"));
        }
        candidate.work = total_work;
        candidate.heat = total_heat;
        poll(&mut cancelled)?;
        *history = candidate;
        Ok(TensorKelvinStressTemperatureStep {
            response,
            start: schedule.start,
            end: schedule.end,
            substeps: schedule.substeps,
            maximum_log_shift_change: schedule.variation_per_step,
        })
    }
}
