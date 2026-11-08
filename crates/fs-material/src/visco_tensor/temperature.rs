//! Physical-time temperature ramps through the original tensor Prony history.
//!
//! Only viscous kinetics vary: tau_j(T)=a_T tau_j(ref) and C_kv(T)=a_T C_kv(ref).
//! Equilibrium/branch tensors, their factors and every accepted memory coordinate
//! retain their original reference identity. A midpoint-frozen kinetic law on
//! bounded subintervals reuses the exact affine-strain/positive-heat transaction.
//! This is a second-order approximation to a prescribed thermal ramp, not an
//! error-certified coupled thermal/mechanical solve or a caloric constitutive law.

use super::{
    Result, Tangent6, TensorProny, TensorPronyError, TensorPronyHistory, TensorPronyStep, Voigt,
    apply, finite, invalid, mandel, poll, tensor, transpose_add,
};
use crate::visco::temperature::{
    TemperatureShift, TemperatureShiftState, scaled, shifted_tensor_frequency_band,
};
use crate::visco::{PronyRamp, ViscoError};
use fs_math::dd::Dd;

/// Prescribed linear absolute temperature as a function of physical time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureRamp {
    /// Absolute temperature at the interval start [K].
    pub start_k: f64,
    /// Absolute temperature at the interval end [K].
    pub end_k: f64,
    /// Original physical interval [s]. A temperature change requires positive time.
    pub duration_s: f64,
}

/// Explicit coefficient-resolution and work limits; neither is an error certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureRampControls {
    /// Maximum analytical variation estimate for ln(a_T) within a subinterval.
    /// Positive and finite. Refine this together with the caller's strain history.
    pub maximum_log_shift_change: f64,
    /// Maximum exact constitutive subintervals, checked before advancing memory.
    pub maximum_substeps: usize,
}

/// Accepted physical observations of the midpoint-frozen thermal approximation.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorPronyTemperatureStep {
    /// Endpoint stress uses the actual endpoint temperature. Mean stress is
    /// weighted in physical time; work and heat are summed from local increments,
    /// independently of any much larger cumulative history ledgers.
    pub response: TensorPronyStep,
    /// Exact calibration evaluation at the prescribed start temperature.
    pub start: TemperatureShiftState,
    /// Exact calibration evaluation at the prescribed endpoint temperature.
    pub end: TemperatureShiftState,
    /// Constitutive subintervals actually advanced.
    pub substeps: usize,
    /// Analytical derivative estimate for maximum subinterval ln(a_T) variation.
    /// Excludes floating-point effects and does not bound trajectory error.
    pub maximum_log_shift_change: f64,
}

fn shift_error(error: ViscoError) -> TensorPronyError {
    match error {
        ViscoError::Parameters { what } => invalid(what),
        _ => invalid("temperature shift leaves the original material applicability"),
    }
}

pub(super) fn scaled_parameter(value: f64, scale: f64) -> Result<f64> {
    scaled(value, scale).map_err(shift_error)
}

// Stack-only coefficient scaling: no material cloning, factorization or refitting.
// Scaling the existing Kelvin factor by sqrt(a_T) retains independent positive heat.
pub(super) fn scaled_kelvin(
    viscosity: &Tangent6,
    factor: &Tangent6,
    scale: f64,
) -> Result<(Tangent6, Tangent6)> {
    let mut viscosity = *viscosity;
    let mut factor = *factor;
    let root = fs_math::det::sqrt(scale);
    for row in &mut viscosity {
        for value in row {
            *value = scaled_parameter(*value, scale)?;
        }
    }
    for row in &mut factor {
        for value in row {
            *value = scaled_parameter(*value, root)?;
        }
    }
    Ok((viscosity, factor))
}

fn admit_scale(
    model: &TensorProny,
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    for branch in model.branches() {
        poll(cancelled)?;
        scaled_parameter(branch.tau_s, scale)?;
    }
    if let Some(viscosity) = model.kv_pa_s() {
        let factor = model
            .kv_factor
            .as_ref()
            .ok_or_else(|| invalid("missing admitted Kelvin viscosity factor"))?;
        scaled_kelvin(viscosity, factor, scale)?;
    }
    shifted_tensor_frequency_band(model.band_rad_s, scale).map_err(shift_error)?;
    Ok(())
}

// Observe the accepted original branch coordinates directly. Subtracting a large
// old Kelvin stress from a rounded total could erase much smaller elastic stress.
fn endpoint_stress(
    history: &TensorPronyHistory,
    rate: &Voigt,
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Voigt> {
    let equilibrium = &history.model.equilibrium_gram;
    let mut stress = [0.0; 6];
    transpose_add(
        equilibrium,
        &apply(equilibrium, &history.strain),
        &mut stress,
    );
    for ((q, branch), factor) in history
        .q
        .iter()
        .zip(history.model.branches())
        .zip(&history.model.factors)
    {
        poll(cancelled)?;
        transpose_add(&branch.resolved_gram_rows(factor), q, &mut stress);
    }
    if let Some(viscosity) = history.model.kv_pa_s() {
        let factor = history
            .model
            .kv_factor
            .as_ref()
            .ok_or_else(|| invalid("missing admitted Kelvin viscosity factor"))?;
        let (viscosity, _) = scaled_kelvin(viscosity, factor, scale)?;
        let viscous = apply(&viscosity, rate);
        for i in 0..6 {
            stress[i] += viscous[i];
        }
    }
    if !finite(&stress) {
        return Err(invalid("unrepresentable endpoint temperature stress"));
    }
    Ok(stress)
}

// The same physical-time partition serves strain, stress and mixed control.
pub(super) struct TemperatureSchedule {
    pub(super) start: TemperatureShiftState,
    pub(super) end: TemperatureShiftState,
    pub(super) substeps: usize,
    pub(super) variation_per_step: f64,
    ramp: TemperatureRamp,
    count: f64,
}

impl TemperatureSchedule {
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    pub(super) fn new(
        model: &TensorProny,
        ramp: TemperatureRamp,
        calibration: &TemperatureShift,
        controls: TemperatureRampControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if !ramp.duration_s.is_finite()
            || ramp.duration_s < 0.0
            || !controls.maximum_log_shift_change.is_finite()
            || controls.maximum_log_shift_change <= 0.0
            || controls.maximum_substeps == 0
        {
            return Err(invalid(
                "finite nonnegative duration and positive thermal work controls required",
            ));
        }
        let start = calibration.evaluate(ramp.start_k).map_err(shift_error)?;
        let end = calibration.evaluate(ramp.end_k).map_err(shift_error)?;
        if ramp.duration_s == 0.0 && ramp.start_k != ramp.end_k {
            return Err(invalid(
                "changing temperature requires positive physical duration",
            ));
        }
        admit_scale(model, start.factor(), cancelled)?;
        if ramp.start_k != ramp.end_k {
            admit_scale(model, end.factor(), cancelled)?;
        }
        // Monotone shift/derivative laws and convex temperature support admit the
        // whole path from these endpoint checks, without a sampling-band fiction.
        let variation = (ramp.end_k - ramp.start_k).abs()
            * start
                .log_derivative_per_k()
                .abs()
                .max(end.log_derivative_per_k().abs());
        let count = (variation / controls.maximum_log_shift_change)
            .ceil()
            .max(1.0);
        if !variation.is_finite()
            || !count.is_finite()
            || count > controls.maximum_substeps as f64
            || count >= usize::MAX as f64
            || count > 4_503_599_627_370_496.0
        {
            return Err(invalid("temperature ramp substep budget exhausted"));
        }
        let substeps = count as usize;
        let variation_per_step = variation / count;
        Ok(Self {
            start,
            end,
            substeps,
            variation_per_step,
            ramp,
            count,
        })
    }

    #[allow(clippy::cast_precision_loss)]
    pub(super) fn interval(
        &self,
        step: usize,
        elapsed: f64,
        calibration: &TemperatureShift,
    ) -> Result<(f64, f64, f64)> {
        let ramp = self.ramp;
        let endpoint = if step + 1 == self.substeps {
            ramp.duration_s
        } else {
            ((step + 1) as f64 / self.count) * ramp.duration_s
        };
        let h = endpoint - elapsed;
        let midpoint = (0.5 * h).mul_add(1.0, elapsed);
        if !h.is_finite() || h <= 0.0 || midpoint <= elapsed || midpoint >= endpoint {
            return Err(invalid(
                "thermal subinterval is below representable physical time",
            ));
        }
        let temperature =
            (midpoint / ramp.duration_s).mul_add(ramp.end_k - ramp.start_k, ramp.start_k);
        let shift = calibration.evaluate(temperature).map_err(shift_error)?;
        Ok((endpoint, h, shift.factor()))
    }
}

impl TensorPronyHistory {
    /// Advance linear physical strain and temperature ramps without resetting memory.
    ///
    /// The immutable history law is declared at the calibration's reference
    /// temperature. Every call scales that ORIGINAL law, never a previously
    /// shifted result. All moduli, factors, branch order and stored coordinates
    /// remain fixed; thermal expansion and modulus changes are not inferred.
    ///
    /// Uniform physical-time subintervals freeze kinetics at their midpoint.
    /// Their count uses |delta T| max |d ln(a_T)/dT| / the declared variation
    /// limit; both supported shift laws attain the derivative bound at an
    /// endpoint. The original exact affine-strain update supplies each local
    /// work/heat increment. Temperature/strain sampling still requires refinement:
    /// a coefficient-variation limit is not a certified physical error bound.
    ///
    /// Constant reference temperature delegates to the original advance, including
    /// its zero-time strain-jump behavior. Zero time with changing temperature
    /// refuses. Temperature changes at fixed strain retain branch storage and
    /// release only the existing dashpot heat; no caloric/entropy law is added.
    ///
    /// # Errors
    /// Refuses invalid ramp/controls, temperatures outside calibration, unrepresentable
    /// shifted coefficients or physical subintervals, source/strain-domain exits,
    /// work-budget exhaustion, overflow, or cancellation. The entire incoming
    /// history remains unchanged on every refusal, including the final poll.
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    pub fn advance_temperature_ramp(
        &mut self,
        strain_end: Voigt,
        ramp: TemperatureRamp,
        calibration: &TemperatureShift,
        controls: TemperatureRampControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorPronyTemperatureStep> {
        let schedule =
            TemperatureSchedule::new(&self.model, ramp, calibration, controls, &mut cancelled)?;
        let TemperatureSchedule {
            start,
            end,
            substeps,
            variation_per_step,
            ..
        } = schedule;
        if ramp.start_k == ramp.end_k {
            // This owner already publishes atomically after its final poll. There
            // are no fallible operations after it, retaining reference-step bits.
            let response = if start.factor() == 1.0 {
                self.advance(strain_end, ramp.duration_s, &mut cancelled)?
            } else {
                self.advance_with_kinetics(
                    strain_end,
                    ramp.duration_s,
                    start.factor(),
                    &mut cancelled,
                    |r, _| PronyRamp::new(r),
                )?
            };
            return Ok(TensorPronyTemperatureStep {
                response,
                start,
                end,
                substeps: 1,
                maximum_log_shift_change: 0.0,
            });
        }
        let target = mandel(&strain_end);
        self.model
            .validate_mandel_strain(&target)
            .map_err(|_| invalid("strain outside original tensor material domain"))?;
        let delta = core::array::from_fn::<_, 6, _>(|i| target[i] - self.strain[i]);
        let rate = delta.map(|value| value / ramp.duration_s);
        if !finite(&delta) || !finite(&rate) {
            return Err(invalid("unrepresentable physical strain ramp"));
        }
        let mut candidate = self.clone();
        let mut work = Dd::ZERO;
        let mut heat = Dd::ZERO;
        let mut mean = [Dd::ZERO; 6];
        let mut elapsed = 0.0;
        for step in 0..substeps {
            poll(&mut cancelled)?;
            let (endpoint, h, scale) = schedule.interval(step, elapsed, calibration)?;
            let next_strain = if step + 1 == substeps {
                strain_end
            } else {
                tensor(&core::array::from_fn(|i| {
                    (endpoint / ramp.duration_s).mul_add(delta[i], self.strain[i])
                }))
            };
            let local = candidate.advance_with_kinetics(
                next_strain,
                h,
                scale,
                &mut cancelled,
                |r, _| PronyRamp::new(r),
            )?;
            work = work + Dd::from_f64(local.work);
            heat = heat + Dd::from_f64(local.dissipated);
            for i in 0..6 {
                mean[i] = mean[i]
                    + Dd::from_f64(h / ramp.duration_s) * Dd::from_f64(local.mean_stress[i]);
            }
            elapsed = endpoint;
        }
        // Maxwell stress observes the actual carried endpoint memory; Kelvin
        // stress observes the actual END temperature and physical strain rate.
        candidate.stress = endpoint_stress(&candidate, &rate, end.factor(), &mut cancelled)?;
        let work = work.to_f64();
        let heat = heat.to_f64();
        let mean_stress = mean.map(Dd::to_f64);
        let residual = work - (candidate.stored - self.stored) - heat;
        let total_work = self.work + work;
        let total_heat = self.heat + heat;
        if !finite(&mean_stress)
            || heat < 0.0
            || ![work, heat, residual, total_work, total_heat]
                .iter()
                .all(|v| v.is_finite())
        {
            return Err(invalid("unrepresentable temperature-ramp work or heat"));
        }
        candidate.work = total_work;
        candidate.heat = total_heat;
        let response = TensorPronyStep {
            stress: tensor(&candidate.stress),
            mean_stress,
            stored_energy: candidate.stored,
            work,
            dissipated: heat,
            balance_residual: residual,
        };
        poll(&mut cancelled)?;
        *self = candidate;
        Ok(TensorPronyTemperatureStep {
            response,
            start,
            end,
            substeps,
            maximum_log_shift_change: variation_per_step,
        })
    }
}
