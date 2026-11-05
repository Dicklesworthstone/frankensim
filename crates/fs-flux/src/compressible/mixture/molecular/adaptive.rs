//! Physical local-accuracy control shared by the original Cartesian RKL2 owners.
//!
//! The finite-network owner defines the accuracy controls. Native field owners
//! retain all source evaluations, Legendre stages and signed face quadrature.
//! A probe uses one identical RKL2 scheme for the full step and both half steps;
//! the difference divided by three estimates the local error of the two halves.
//! Accepted fields are never extrapolated or clipped.

use super::interval::Stage;
use super::stabilized::{
    Action, IntervalFailure, Rkl2Scheme, StabilizedInterval, reserve_evaluation,
};
use super::{DissipationError, MixtureGas, allocated, poll};
use fs_qty::Temperature;

/// Original molecular local-accuracy policy. Cartesian owners apply its
/// velocity tolerance separately to every spatial component.
pub use fs_gas::multicomponent::network::dissipation::AdaptiveNetworkDissipationControl as AdaptiveDissipationControl;

mod planar;
pub use planar::AdaptiveMolecularReport;

/// Reused conservative anchor, full-probe physical coordinates and face quadratures.
#[derive(Default)]
pub(in crate::compressible) struct Workspace {
    pub anchor: Vec<f64>,
    /// Original species densities followed by admitted velocity components and
    /// temperature. Capturing the native endpoint cache avoids another EOS solve.
    pub coarse: Vec<f64>,
    pub coarse_boundary: Vec<f64>,
    pub first_boundary: Vec<f64>,
}
impl Workspace {
    pub fn required_bytes(
        state_len: usize,
        boundary_len: usize,
    ) -> Result<usize, DissipationError> {
        state_len
            .checked_add(boundary_len)
            .and_then(|n| n.checked_mul(2 * size_of::<f64>()))
            .ok_or(DissipationError::Invalid(
                "adaptive molecular workspace size overflow",
            ))
    }

    pub fn additional_bytes(
        &self,
        state_len: usize,
        boundary_len: usize,
    ) -> Result<usize, DissipationError> {
        let mut extra = 0usize;
        for (saved, count) in [
            (&self.anchor, state_len),
            (&self.coarse, state_len),
            (&self.coarse_boundary, boundary_len),
            (&self.first_boundary, boundary_len),
        ] {
            if saved.len() != count {
                extra = count
                    .checked_mul(size_of::<f64>())
                    .and_then(|bytes| extra.checked_add(bytes))
                    .ok_or(DissipationError::Invalid(
                        "adaptive molecular workspace size overflow",
                    ))?;
            }
        }
        Ok(extra)
    }

    /// The owning field admits the complete pending storage before this call.
    pub fn reserve(
        &mut self,
        state_len: usize,
        boundary_len: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        for (saved, count) in [
            (&mut self.anchor, state_len),
            (&mut self.coarse, state_len),
            (&mut self.coarse_boundary, boundary_len),
            (&mut self.first_boundary, boundary_len),
        ] {
            poll(cancelled)?;
            if saved.len() != count {
                *saved = allocated(count, 0.0)?;
            }
        }
        poll(cancelled)
    }

    pub fn workspace_bytes(&self) -> usize {
        (self.anchor.capacity()
            + self.coarse.capacity()
            + self.coarse_boundary.capacity()
            + self.first_boundary.capacity())
            * size_of::<f64>()
    }
}

#[derive(Clone, Copy, Default)]
pub(in crate::compressible) struct AdaptiveInterval {
    pub stabilized: StabilizedInterval,
    pub maximum_estimated_local_error: f64,
    pub temporal_rejections: usize,
}

/// The callback retains or restores its private anchor on every refused probe.
/// Only a trial with error <= 1 may accumulate its two-half face quadrature and
/// replace the private field. Its returned rate is the largest of the coarse
/// rate and HALF each fine rate, so multiplying by full dt retains actual probe
/// parabolic numbers. Source residual extrema include all three probes.
pub(in crate::compressible) fn advance<C: FnMut() -> bool>(
    duration: f64,
    controls: AdaptiveDissipationControl,
    evaluations: &mut usize,
    cancelled: &mut C,
    mut evaluate: impl FnMut(Action, f64, &mut usize, &mut C) -> Result<(Stage, f64), DissipationError>,
) -> Result<AdaptiveInterval, DissipationError> {
    poll(cancelled)?;
    controls.validate()?;
    if !duration.is_finite() || duration < 0.0 {
        return Err(DissipationError::Invalid(
            "nonnegative finite adaptive molecular duration required",
        ));
    }
    let stabilized = controls.stabilized;
    let base = stabilized.transport;
    let maximum_scheme = Rkl2Scheme::new(stabilized.maximum_stages)?;
    let initial_evaluations = *evaluations;
    let mut result = AdaptiveInterval::default();
    let mut next_step = base.maximum_step_s;
    while result.stabilized.interval.elapsed < duration {
        poll(cancelled)?;
        work_limit(&result, base.maximum_steps)?;
        let elapsed = result.stabilized.interval.elapsed;
        reserve_evaluation(
            evaluations,
            stabilized.maximum_operator_evaluations,
            elapsed,
        )?;
        let (first, _) = evaluate(Action::First, elapsed, evaluations, cancelled)?;
        if !first.rate.is_finite() || first.rate < 0.0 {
            return Err(DissipationError::Invalid(
                "adaptive molecular parabolic scale overflow",
            ));
        }
        let remaining = duration - elapsed;
        let mut dt = remaining.min(next_step);
        if first.rate > 0.0 {
            dt = dt.min(base.parabolic_cfl * maximum_scheme.stability_extent() / first.rate);
        }
        for retry in 0..=base.maximum_halvings {
            poll(cancelled)?;
            work_limit(&result, base.maximum_steps)?;
            if !dt.is_finite()
                || 0.5 * dt <= 0.0
                || elapsed + dt <= elapsed
                || (dt < base.minimum_step_s && dt != remaining)
            {
                return Err(DissipationError::Invalid(
                    "unrepresentable adaptive molecular substep",
                ));
            }
            let mut scheme = maximum_scheme;
            for stages in 2..=stabilized.maximum_stages {
                let candidate = Rkl2Scheme::new(stages)?;
                if dt * first.rate <= base.parabolic_cfl * candidate.stability_extent() {
                    scheme = candidate;
                    break;
                }
            }
            let attempt = evaluate(
                Action::Trial {
                    dt,
                    scheme,
                    first,
                    cached: retry == 0,
                },
                elapsed,
                evaluations,
                cancelled,
            )
            .and_then(|(stage, error)| {
                if error > 1.0 {
                    Err(DissipationError::EstimatedTemporalError {
                        observed: error,
                        maximum: 1.0,
                    })
                } else if !error.is_finite() || error < 0.0 {
                    Err(DissipationError::Invalid(
                        "non-finite adaptive molecular error estimate",
                    ))
                } else {
                    Ok((stage, error))
                }
            });
            match attempt {
                Ok((stage, error)) => {
                    let number = dt * stage.rate;
                    super::interval::parabolic(
                        number,
                        base.parabolic_cfl * scheme.stability_extent(),
                    )?;
                    let accepted = &mut result.stabilized;
                    accepted.interval.maximum_parabolic_number =
                        accepted.interval.maximum_parabolic_number.max(number);
                    accepted.maximum_scaled_parabolic_number = accepted
                        .maximum_scaled_parabolic_number
                        .max(number / scheme.stability_extent());
                    accepted.maximum_stages_used =
                        accepted.maximum_stages_used.max(scheme.stage_count());
                    accepted.interval.maximum_gradient_residual_per_m = accepted
                        .interval
                        .maximum_gradient_residual_per_m
                        .max(stage.gradient_residual);
                    accepted.interval.maximum_mass_frame_residual_kg_m2_s = accepted
                        .interval
                        .maximum_mass_frame_residual_kg_m2_s
                        .max(stage.mass_frame_residual);
                    accepted.interval.elapsed = if dt == remaining {
                        duration
                    } else {
                        elapsed + dt
                    };
                    accepted.interval.accepted_steps += 1;
                    result.maximum_estimated_local_error =
                        result.maximum_estimated_local_error.max(error);
                    next_step = (dt * factor(error, true)).min(base.maximum_step_s);
                    break;
                }
                Err(error) if error.is_terminal() => return Err(error),
                Err(error) => {
                    result.stabilized.interval.rejected_steps += 1;
                    let shrink =
                        if let DissipationError::EstimatedTemporalError { observed, .. } = error {
                            result.temporal_rejections += 1;
                            factor(observed, false)
                        } else {
                            0.5
                        };
                    if retry == base.maximum_halvings || dt * shrink < base.minimum_step_s {
                        return Err(error.step_failure(elapsed, dt));
                    }
                    dt *= shrink;
                }
            }
        }
    }
    poll(cancelled)?;
    result.stabilized.operator_evaluations = *evaluations - initial_evaluations;
    Ok(result)
}

fn work_limit(result: &AdaptiveInterval, maximum: usize) -> Result<(), DissipationError> {
    let interval = result.stabilized.interval;
    if interval.accepted_steps + interval.rejected_steps >= maximum {
        Err(DissipationError::WorkLimit {
            elapsed_time_s: interval.elapsed,
            accepted_steps: interval.accepted_steps,
            rejected_steps: interval.rejected_steps,
        })
    } else {
        Ok(())
    }
}

fn factor(error: f64, accepted: bool) -> f64 {
    if error == 0.0 {
        2.0
    } else {
        (0.9 / error.cbrt()).clamp(0.2, if accepted { 2.0 } else { 0.8 })
    }
}

fn ratio(coarse: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64, DissipationError> {
    let tolerance = absolute + relative * coarse.abs().max(fine.abs());
    if !tolerance.is_finite() || tolerance < 0.0 || !coarse.is_finite() || !fine.is_finite() {
        return Err(DissipationError::Invalid(
            "non-finite adaptive molecular error scale",
        ));
    }
    let difference = (fine - coarse).abs();
    Ok(if difference == 0.0 {
        0.0
    } else {
        (difference / tolerance) / 3.0
    })
}

/// Primitive values decoded by the original dimension-specific caloric owner.
#[derive(Clone, Copy)]
pub(in crate::compressible) struct PhysicalState<const D: usize> {
    pub density: f64,
    pub temperature: f64,
    pub velocity: [f64; D],
}

pub(in crate::compressible) fn state_error<const D: usize>(
    gas: &MixtureGas,
    coarse: &[f64],
    fine: &[f64],
    coarse_primitive: PhysicalState<D>,
    fine_primitive: PhysicalState<D>,
    controls: AdaptiveDissipationControl,
) -> Result<f64, DissipationError> {
    if !coarse_primitive.density.is_finite()
        || coarse_primitive.density <= 0.0
        || !fine_primitive.density.is_finite()
        || fine_primitive.density <= 0.0
    {
        return Err(DissipationError::Invalid(
            "adaptive molecular physical density",
        ));
    }
    let mut error = 0.0_f64;
    for k in 0..gas.species_count() {
        error = error.max(ratio(
            coarse[k],
            fine[k],
            controls.absolute_partial_density_tolerance_kg_m3,
            controls.relative_tolerance,
        )?);
    }
    error = error.max(ratio(
        coarse_primitive.temperature,
        fine_primitive.temperature,
        controls.absolute_temperature_tolerance_k,
        controls.relative_tolerance,
    )?);
    for axis in 0..D {
        error = error.max(ratio(
            coarse_primitive.velocity[axis],
            fine_primitive.velocity[axis],
            controls.absolute_velocity_tolerance_m_s,
            controls.relative_tolerance,
        )?);
    }
    Ok(error)
}

/// Inspect each signed face separately. Species/thermal/impulse absolute scales
/// use the incident cell's actual volume, heat capacity and mass. Stored energy
/// already uses the source's per-species Tmin sensible chart, so arbitrary
/// formation constants do not loosen the thermal quadrature tolerance.
#[allow(clippy::too_many_arguments)]
pub(in crate::compressible) fn face_error<const D: usize>(
    gas: &MixtureGas,
    coarse: &[f64],
    first: &[f64],
    second: &[f64],
    fine_densities: &[f64],
    fine_temperature: f64,
    cell_volume: f64,
    controls: AdaptiveDissipationControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, DissipationError> {
    let mut density = 0.0;
    let mut heat_capacity = 0.0;
    let mut correction = 0.0;
    for (partial, model) in fine_densities.iter().zip(gas.source().models()) {
        poll(cancelled)?;
        density += partial;
        if *partial == 0.0 {
            continue;
        }
        let properties = model
            .evaluate_properties(Temperature::new(fine_temperature))
            .map_err(fs_gas::multicomponent::Error::from)?;
        let cv = properties.cp().value() - fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
        let term = partial / model.molar_mass().value() * cv;
        let adjusted = term - correction;
        let next = heat_capacity + adjusted;
        correction = (next - heat_capacity) - adjusted;
        heat_capacity = next;
    }
    let mass = cell_volume * density;
    let capacity = cell_volume * heat_capacity;
    if !mass.is_finite() || mass <= 0.0 || !capacity.is_finite() || capacity <= 0.0 {
        return Err(DissipationError::Invalid(
            "adaptive molecular boundary storage scale",
        ));
    }
    let n = gas.species_count();
    let mut error = 0.0_f64;
    for k in 0..n + D + 1 {
        let absolute = if k < n {
            cell_volume * controls.absolute_partial_density_tolerance_kg_m3
        } else if k < n + D {
            mass * controls.absolute_velocity_tolerance_m_s
        } else {
            capacity * controls.absolute_temperature_tolerance_k
        };
        error = error.max(ratio(
            coarse[k],
            first[k] + second[k],
            absolute,
            controls.relative_tolerance,
        )?);
    }
    Ok(error)
}
