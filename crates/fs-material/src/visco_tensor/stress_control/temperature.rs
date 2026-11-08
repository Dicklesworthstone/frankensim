//! Physical-time thermal creep and mixed control on the original prepared law.

use super::{
    PreparedTensorStressControl, Result, TensorPronyHistory, TensorStressRampStep, Voigt, finite,
    invalid, mandel, poll, tensor,
};
use crate::visco::temperature::{TemperatureShift, TemperatureShiftState};
use crate::visco_tensor::temperature::{
    TemperatureRamp, TemperatureRampControls, TemperatureSchedule,
};
use fs_math::dd::Dd;
use std::sync::Arc;

/// Accepted temperature-dependent stress or mixed-control interval.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorStressTemperatureStep {
    /// Full original-source response. Increments sum local reports and bounds
    /// enclose every subinterval; heat is the independent original dashpot integral.
    pub response: TensorStressRampStep,
    /// Calibration evaluated at the actual initial temperature.
    pub start: TemperatureShiftState,
    /// Calibration evaluated at the actual final temperature.
    pub end: TemperatureShiftState,
    /// Exact constitutive subintervals actually advanced.
    pub substeps: usize,
    /// Analytical subinterval ln(a_T) variation estimate, not a trajectory bound.
    pub maximum_log_shift_change: f64,
}

fn reduced_duration(physical: f64, scale: f64) -> Result<f64> {
    let reduced = if scale == 1.0 {
        physical
    } else {
        physical / scale
    };
    if !reduced.is_finite() || (physical > 0.0 && reduced <= 0.0) {
        return Err(invalid(
            "unrepresentable reduced-time constitutive interval",
        ));
    }
    Ok(reduced)
}

impl PreparedTensorStressControl {
    /// Advance prescribed stress and temperature, retaining the original memory.
    ///
    /// Loading and temperature are linear in PHYSICAL time. Each admitted
    /// midpoint-temperature subinterval reuses this prepared controller for
    /// duration h/a_T. Common relaxation-time scaling makes this the original
    /// reference law in reduced time; work and dashpot heat need no rescaling.
    /// No factorization, memory translation or relaxation-time grouping occurs.
    ///
    /// The original law is declared at the calibration reference temperature.
    /// Fixed moduli and the existing no-Kelvin admission remain binding. No
    /// thermal expansion, caloric feedback or trajectory-error certificate is
    /// inferred. Refine the coefficient limit and the physical loading history.
    /// Constant reference temperature retains the original response bits.
    ///
    /// # Errors
    /// Refuses invalid/binding-mismatched data, thermal applicability or budget
    /// exits, unrepresentable physical/reduced intervals, any original whole-path
    /// strain or source/energy refusal, overflow or cancellation. No history is
    /// published until the entire physical interval and final poll succeed.
    #[allow(clippy::too_many_lines)] // One atomic physical-time constitutive interval.
    pub fn advance_temperature_ramp(
        &self,
        history: &mut TensorPronyHistory,
        prescribed: Voigt,
        ramp: TemperatureRamp,
        calibration: &TemperatureShift,
        controls: TemperatureRampControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorStressTemperatureStep> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.model, &history.model) || !finite(&prescribed) {
            return Err(invalid(
                "thermal control requires finite controls and original history",
            ));
        }
        let schedule =
            TemperatureSchedule::new(&self.model, ramp, calibration, controls, &mut cancelled)?;
        if ramp.start_k == ramp.end_k {
            let duration = reduced_duration(ramp.duration_s, schedule.start.factor())?;
            let response =
                self.advance_controlled(history, prescribed, duration, &mut cancelled)?;
            return Ok(TensorStressTemperatureStep {
                response,
                start: schedule.start,
                end: schedule.end,
                substeps: 1,
                maximum_log_shift_change: 0.0,
            });
        }
        let initial = core::array::from_fn::<_, 6, _>(|i| {
            if self.stress_components[i] {
                history.stress[i]
            } else {
                history.strain[i]
            }
        });
        let target = mandel(&prescribed);
        let delta = core::array::from_fn::<_, 6, _>(|i| target[i] - initial[i]);
        if !finite(&target) || !finite(&delta) {
            return Err(invalid("unrepresentable physical control ramp"));
        }
        let mut candidate = history.clone();
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
                prescribed
            } else {
                tensor(&core::array::from_fn(|i| {
                    (endpoint / ramp.duration_s).mul_add(delta[i], initial[i])
                }))
            };
            let local = self.advance_controlled(
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
        let mut response = last.ok_or_else(|| invalid("missing thermal control interval"))?;
        response.strain_increment = increment.map(Dd::to_f64);
        response.stored_energy_change = storage_change.to_f64();
        response.work = work.to_f64();
        response.dissipated = heat.to_f64();
        response.balance_residual = (work - storage_change - heat).to_f64();
        response.strain_bounds = bounds;
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
            return Err(invalid("unrepresentable thermal control work or heat"));
        }
        candidate.work = total_work;
        candidate.heat = total_heat;
        poll(&mut cancelled)?;
        *history = candidate;
        Ok(TensorStressTemperatureStep {
            response,
            start: schedule.start,
            end: schedule.end,
            substeps: schedule.substeps,
            maximum_log_shift_change: schedule.variation_per_step,
        })
    }
}
