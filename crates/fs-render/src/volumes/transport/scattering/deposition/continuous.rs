//! Continuous absorption weighting for scalar incident-power heating.
//!
//! A flight samples scattering optical depth alone. Every traversed native
//! cell absorbs `P_in * -expm1(-rho * kappa_a * length)`; its surviving
//! radiation continues without an extra albedo debit at a scattering event.
//! This is the continuous-absorption estimator of Hayakawa et al.,
//! JOSA A 31 (2014), 301-311, DOI:10.1364/JOSAA.31.000301.
//!
//! Thin absorbers receive positive heat even in histories with no collision.
//! Variance and cost depend on the medium: this is an alternative to collision
//! capture, not a claim of universally better efficiency.

use super::*;
use crate::volumes::transport::cells::CellRay;
use crate::volumes::transport::deposition::BeamOpticalDepth;
use crate::volumes::transport::positive_product;
use crate::volumes::transport::scattering::{CellScatteringGrid, divide_product};

impl ScatteringGrid<'_> {
    /// Estimate scalar heat by analytically absorbing along scattering flights.
    ///
    /// The physical coefficients, phase laws, native cells and work limits are
    /// unchanged. No-scattering lanes use exact deterministic Beer transport.
    /// In sampled lanes, surviving power at the last admitted scattering event
    /// remains unresolved, as in `deposit_multiple_scattering_power`. A positive
    /// absorption contribution below representable watts likewise retains the
    /// entire incoming packet and marks its precision limitation.
    /// The report can be consumed directly by native `ScatteringHeating`.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_continuous_absorption_power<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
        self.deposit_continuous_observed(
            cx,
            origin,
            unit_direction,
            wavelengths_nm,
            incident_power_w,
            config,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn deposit_continuous_observed<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
        momentum: Option<&mut MomentumTallies<N>>,
    ) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, unit_direction, f64::INFINITY)?;
        if N == 0 {
            return Err(VolumeTransportError::Sample);
        }
        if self.absorption.absorption_palette().is_some() || self.scattering_palette().is_some() {
            let medium = PaletteOptics::try_new(self, *wavelengths_nm, &mut || checkpoint(cx))?;
            return deposit(
                &medium,
                &mut || checkpoint(cx),
                origin,
                unit_direction,
                incident_power_w,
                config,
                momentum,
            );
        }
        let (absorption, scattering, extinction) = self.coefficients(wavelengths_nm)?;
        let medium = UniformOptics {
            medium: self,
            wavelengths: *wavelengths_nm,
            absorption,
            scattering,
            extinction,
        };
        deposit(
            &medium,
            &mut || checkpoint(cx),
            origin,
            unit_direction,
            incident_power_w,
            config,
            momentum,
        )
    }
}

impl<const N: usize> CellScatteringGrid<'_, '_, N> {
    /// Continuous absorption with the original cell-local coefficients and phase.
    pub fn deposit_continuous_absorption_power(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
        self.deposit_continuous_absorption_power_with_checkpoint(
            &mut || checkpoint(cx),
            origin,
            unit_direction,
            incident_power_w,
            config,
        )
    }

    /// Hosting-solver variant preserving its exact cancellation/budget refusal.
    pub fn deposit_continuous_absorption_power_with_checkpoint(
        &self,
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
        origin: Point3,
        unit_direction: Vec3,
        incident_power_w: &[f64; N],
        config: ScatteringDepositionConfig,
    ) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
        deposit(
            self,
            check,
            origin,
            unit_direction,
            incident_power_w,
            config,
            None,
        )
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn deposit<const N: usize>(
    medium: &impl PreparedOptics<N>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    origin: Point3,
    unit_direction: Vec3,
    incident_power_w: &[f64; N],
    config: ScatteringDepositionConfig,
    momentum: Option<&mut MomentumTallies<N>>,
) -> Result<ScatteringPowerDeposition<N>, VolumeTransportError> {
    deposit_prepared_sampled(
        medium,
        check,
        origin,
        unit_direction,
        incident_power_w,
        config,
        momentum,
        &[ScalarPhase; N],
        AbsorptionEstimator::Continuous,
    )
}

/// One DDA per flight, no copied grid, temporary cell list or second traversal.
#[allow(clippy::too_many_arguments)]
pub(super) fn trace_history<const N: usize>(
    medium: &impl PreparedOptics<1>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    source: (Point3, Vec3, f64),
    limits: (u32, u64),
    traced_segments: &mut u64,
    random: &mut impl FnMut() -> f64,
    capture: &mut impl FnMut(usize, f64) -> Result<(), VolumeTransportError>,
    mut momentum: Option<&mut HistoryMomentum<'_, N>>,
) -> Result<HistoryOutcome, VolumeTransportError> {
    let (mut position, mut direction, source_w) = source;
    let mut power = source_w;
    let mut depth = BeamOpticalDepth::default();
    let mut captured = PowerSum::default();
    let mut precision_limited = false;
    let mut outcome = [0.0; 4];
    'flights: for collision in 0..limits.0 {
        check()?;
        count_segment(traced_segments, limits.1)?;
        // Carry one sampled optical distance across every crossed cell.
        // Absorbing-only cells attenuate power without consuming this distance.
        let mut remaining = PowerSum {
            sum: -det::ln(random().max(f64::MIN_POSITIVE)),
            correction: 0.0,
        };
        if let Some(mut ray) = CellRay::clip(
            &medium.geometry().grid,
            position,
            direction,
            f64::INFINITY,
            medium.geometry().max_steps,
        )? {
            while let Some(segment) = ray.next_with_checkpoint(check)? {
                let rho = medium.geometry().grid.data[segment.cell];
                if rho == 0.0 || segment.length_m == 0.0 {
                    continue;
                }
                let (a, s, _) = medium.coefficient(segment.cell, 0)?;
                let scattering_depth = if s == 0.0 {
                    0.0
                } else {
                    positive_product(rho, s, segment.length_m)
                };
                let scatters = s > 0.0 && remaining.value() < scattering_depth;
                let length = if scatters {
                    divide_product(remaining.value(), rho, s)
                } else {
                    segment.length_m
                };
                let parameter = segment.start_m + length;
                if scatters
                    && (!parameter.is_finite()
                        || parameter < segment.start_m
                        || parameter >= segment.start_m + segment.length_m
                        || (length > 0.0 && parameter == segment.start_m)
                        || (remaining.value() > 0.0 && length == 0.0))
                {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                if a > 0.0 && length > 0.0 {
                    let tau = positive_product(rho, a, length);
                    let heat = removed_power(power, -det::expm1(-tau), rho, length, a, a);
                    if !heat.is_finite() {
                        return Err(VolumeTransportError::Unrepresentable);
                    }
                    if heat == 0.0 {
                        // Retain the whole incoming packet before an absorption
                        // contribution falls below representable watts. Neither
                        // a rounded survivor nor discarded heat closes its account.
                        outcome[3] = power;
                        precision_limited = true;
                        if let Some(tally) = momentum.as_deref_mut() {
                            tally.leave(2, power, ray.position(segment.start_m), direction)?;
                        }
                        break 'flights;
                    }
                    if heat > 0.0 {
                        captured.add(heat)?;
                        capture(segment.cell, heat)?;
                        if let Some(tally) = momentum.as_deref_mut() {
                            // Along this straight segment, (r + t*d) x d is
                            // constant. Its positive distributed absorption
                            // has exactly the segment-start line of action.
                            tally.capture(
                                segment.cell,
                                ray.position(segment.start_m),
                                heat,
                                0.0,
                                direction,
                                direction,
                            )?;
                        }
                    }
                    // Compensated cumulative depth retains many tiny losses,
                    // including flights that revisit the same physical cell.
                    power = depth.transmit(source_w, tau);
                    if power == 0.0 {
                        break 'flights;
                    }
                }
                if scatters {
                    let next_position = ray.position(parameter);
                    let distance = ray.distance_from_origin(parameter);
                    if !distance.is_finite() || (distance > 0.0 && next_position == position) {
                        return Err(VolumeTransportError::Unrepresentable);
                    }
                    let capped = collision + 1 == limits.0;
                    if capped && momentum.is_none() {
                        outcome[3] = power;
                        break 'flights;
                    }
                    let outgoing = ScalarPhase.sample(
                        medium.phase(segment.cell),
                        medium.wavelengths()[0],
                        direction,
                        check,
                        random,
                    )?;
                    if let Some(tally) = momentum.as_deref_mut() {
                        // Absorption was already integrated along the flight.
                        // Only the elastic change of direction acts here.
                        tally.capture(
                            segment.cell,
                            next_position,
                            0.0,
                            power,
                            direction,
                            outgoing,
                        )?;
                        if capped {
                            tally.leave(2, power, next_position, outgoing)?;
                        }
                    }
                    if capped {
                        outcome[3] = power;
                        break 'flights;
                    }
                    direction = outgoing;
                    position = next_position;
                    continue 'flights;
                }
                remaining.add(-scattering_depth)?;
                if remaining.value() < 0.0 {
                    return Err(VolumeTransportError::Unrepresentable);
                }
            }
        }
        let escape = usize::from(collision != 0);
        outcome[escape] = power;
        if let Some(tally) = momentum.as_deref_mut() {
            tally.leave(escape, power, position, direction)?;
        }
        break;
    }
    outcome[2] = captured.value();
    check()?;
    Ok(HistoryOutcome {
        values: outcome,
        precision_limited,
    })
}

#[cfg(test)]
mod tests;
