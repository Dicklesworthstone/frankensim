//! Stored acoustic energy and wave momentum from the original coherent fields.
//!
//! Fields interfere before products. Storage changes use analytic field
//! increments and the polarized quadratic identity, never a subtraction of
//! two large endpoint energies. The nondispersive rest-medium law remains the
//! original owner's authority; these are not total fluid momentum or body loads.

use super::*;
use fs_fft::harmonic_window::phase_rad_s_dd_from_clock_dd;
use fs_math::det;

/// Work for one allocation-free endpoint or storage-change observation.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicVolumeFieldLimits {
    /// Original entries, including repeated receivers and zero weights.
    pub maximum_quadrature_points: usize,
    /// sum_b F_b*(1+4*Q): one phase pair and four field components per receiver
    /// and frequency. Each field term includes its bounded analytic increment.
    /// Fixed seven-channel products and the source-bank count remain bounded.
    pub maximum_field_terms: usize,
}

/// Instantaneous ensemble-mean storage in the authored physical volume.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicVolumeFieldStorage {
    /// Acoustic compression plus kinetic energy [J].
    pub energy_j: Dd,
    pub wave_momentum: HarmonicWaveMomentumStorage,
}

/// End minus start storage; energy changes may have either sign.
/// No enclosure, matched surface rule or mechanical reaction is inferred.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicVolumeFieldChange {
    pub energy_change_j: Dd,
    pub wave_momentum: HarmonicWaveMomentumTransfer,
}

impl PreparedHarmonicExposure<'_> {
    /// Original-field acoustic energy and linear/angular wave momentum at one
    /// precise source epoch. E=(p²/(rho*c²)+rho*u.u)/2 and g=p*u/c².
    /// Every original frequency owner, including silent lines, must authorize
    /// a homogeneous inviscid nondispersive rest medium for both storage laws.
    pub fn volume_field_storage_from_clock_dd(
        &self,
        quadrature: &[HarmonicVolumeQuadrature],
        moment_origin_m: [f64; 3],
        time_s: Dd,
        limits: HarmonicVolumeFieldLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicVolumeFieldStorage, TransferError> {
        let values = self.volume_fields(
            quadrature,
            moment_origin_m,
            time_s,
            None,
            limits,
            &mut cancelled,
        )?;
        Ok(HarmonicVolumeFieldStorage {
            energy_j: values[0],
            wave_momentum: HarmonicWaveMomentumStorage {
                linear_momentum_ns: [values[1], values[2], values[3]],
                angular_momentum_nms: [values[4], values[5], values[6]],
            },
        })
    }

    /// Analytic change over a finite nonnegative duration, independent of the
    /// DD start epoch. Zero duration gives exact zero after ordinary admission.
    /// Summed field increments retain short intervals and coherent nulls;
    /// constant storage cancels before products. Add outward surface transport
    /// only for an explicitly matched control-volume balance.
    #[allow(clippy::too_many_arguments)]
    pub fn volume_field_change_from_clock_dd(
        &self,
        quadrature: &[HarmonicVolumeQuadrature],
        moment_origin_m: [f64; 3],
        start_s: Dd,
        duration_s: f64,
        limits: HarmonicVolumeFieldLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicVolumeFieldChange, TransferError> {
        let values = self.volume_fields(
            quadrature,
            moment_origin_m,
            start_s,
            Some(duration_s),
            limits,
            &mut cancelled,
        )?;
        Ok(HarmonicVolumeFieldChange {
            energy_change_j: values[0],
            wave_momentum: HarmonicWaveMomentumTransfer {
                linear_momentum_ns: [values[1], values[2], values[3]],
                angular_momentum_nms: [values[4], values[5], values[6]],
            },
        })
    }

    fn volume_fields(
        &self,
        quadrature: &[HarmonicVolumeQuadrature],
        origin: [f64; 3],
        start: Dd,
        duration: Option<f64>,
        limits: HarmonicVolumeFieldLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 7], TransferError> {
        poll(cancelled)?;
        admit(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        admit(quadrature.len(), limits.maximum_quadrature_points)?;
        if !self.wave_momentum || origin.iter().any(|x| !x.is_finite()) {
            return Err(bad(
                "original acoustic owners do not authorize this volume storage query",
            ));
        }
        let start = Dd::from_pair(start.hi, start.lo);
        if !finite(start)
            || duration
                .is_some_and(|dt| !dt.is_finite() || dt < 0.0 || !finite(start + Dd::from_f64(dt)))
        {
            return Err(bad("acoustic volume source interval is invalid"));
        }
        for sample in quadrature {
            poll(cancelled)?;
            if sample.receiver >= self.reference.receivers.len()
                || !sample.volume_m3.is_finite()
                || sample.volume_m3 < 0.0
                || self.reference.receivers[sample.receiver]
                    .position_m
                    .iter()
                    .any(|x| !x.is_finite())
            {
                return Err(bad("acoustic volume field quadrature is invalid"));
            }
        }
        let banks = || core::iter::once(&self.mean).chain(self.factors.iter());
        let lines = banks().try_fold(0usize, |n, bank| {
            n.checked_add(bank.lines.len()).ok_or_else(overflow)
        })?;
        let work = quadrature
            .len()
            .checked_mul(4)
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(lines))
            .ok_or_else(overflow)?;
        admit(work, MAX_PRODUCTS.min(limits.maximum_field_terms))?;
        let mut sums = [Dd::ZERO; 7];
        if quadrature.is_empty() || duration == Some(0.0) {
            poll(cancelled)?;
            return Ok(sums);
        }
        // One cached phase/increment pair per original frequency serves every
        // receiver. The preparation's existing 4096-line cap bounds this fixed
        // 256 KiB scratch; no per-query allocation or receiver array is needed.
        let mut phases = [[[Dd::ZERO; 2]; 2]; MAX_EXPOSURE_LINES];
        for bank in banks() {
            poll(cancelled)?;
            let Some(first) = bank.lines.first() else {
                continue;
            };
            let carrier = Dd::from_f64(first.omega);
            let real = bank.statistics == HarmonicExposureFactorStatistics::Real;
            let rotation = if real {
                phase(carrier, start)?
            } else {
                [Dd::ONE, Dd::ZERO]
            };
            let carrier_delta = if real {
                increment(carrier, duration)?
            } else {
                [Dd::ZERO; 2]
            };
            let carrier_end = [Dd::ONE + carrier_delta[0], carrier_delta[1]];
            for (value, line) in phases.iter_mut().zip(&bank.lines) {
                poll(cancelled)?;
                // Proper quadratic moments are invariant under a common
                // carrier; real fields receive its exact rotation afterward.
                let relative = Dd::from_f64(line.omega) - carrier;
                *value = [phase(relative, start)?, increment(relative, duration)?];
            }
            for sample in quadrature {
                poll(cancelled)?;
                if sample.volume_m3 == 0.0 {
                    continue;
                }
                let mut fields = [[Dd::ZERO; 2]; 4];
                let mut deltas = [[Dd::ZERO; 2]; 4];
                for (line, rotations) in bank.lines.iter().zip(&phases) {
                    poll(cancelled)?;
                    for c in 0..4 {
                        let value =
                            multiply(component(&line.fields[sample.receiver], c), rotations[0]);
                        let delta = multiply(value, rotations[1]);
                        for part in 0..2 {
                            fields[c][part] = fields[c][part] + value[part];
                            deltas[c][part] = deltas[c][part] + delta[part];
                        }
                    }
                }
                if real {
                    for c in 0..4 {
                        let carrier_part = multiply(fields[c], carrier_delta);
                        let relative_part = multiply(deltas[c], carrier_end);
                        deltas[c] = [
                            multiply(
                                [
                                    carrier_part[0] + relative_part[0],
                                    carrier_part[1] + relative_part[1],
                                ],
                                rotation,
                            )[0],
                            Dd::ZERO,
                        ];
                        fields[c] = [multiply(fields[c], rotation)[0], Dd::ZERO];
                    }
                }
                if fields.iter().chain(&deltas).flatten().any(|x| !finite(*x)) {
                    return Err(bad("acoustic volume coherent field is unrepresentable"));
                }
                contract(
                    self.reference,
                    sample,
                    origin,
                    fields,
                    deltas,
                    !real,
                    duration.is_some(),
                    &mut sums,
                )?;
            }
        }
        poll(cancelled)?;
        if sums.iter().any(|x| !finite(*x)) || (duration.is_none() && sums[0].lt(Dd::ZERO)) {
            return Err(bad("acoustic volume field storage is unrepresentable"));
        }
        Ok(sums)
    }
}

fn multiply(a: [Dd; 2], b: [Dd; 2]) -> [Dd; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn phase(omega: Dd, time: Dd) -> Result<[Dd; 2], TransferError> {
    phase_rad_s_dd_from_clock_dd(omega, time).map_err(window_error)
}

fn increment(omega: Dd, duration: Option<f64>) -> Result<[Dd; 2], TransferError> {
    let Some(dt) = duration else {
        return Ok([Dd::ZERO; 2]);
    };
    let angle = omega * Dd::from_f64(dt);
    if !finite(angle) || (omega != Dd::ZERO && dt != 0.0 && angle == Dd::ZERO) {
        return Err(bad("acoustic volume phase increment is unrepresentable"));
    }
    let half = angle * Dd::from_f64(0.5);
    if half == Dd::ZERO {
        if angle != Dd::ZERO {
            return Err(bad("acoustic volume phase increment is unrepresentable"));
        }
        return Ok([Dd::ZERO; 2]);
    }
    let (sine, cosine) = det::sin_cos_dd(half);
    // exp(-i angle)-1: avoid cos(angle)-1 and preserve the small real term.
    let result = [
        -Dd::from_f64(2.0) * sine * sine,
        -Dd::from_f64(2.0) * sine * cosine,
    ];
    // Large fields can restore a finite change from an otherwise tiny phase
    // square. Refuse its lost real increment rather than inventing zero work.
    if result.iter().any(|x| !finite(*x)) || (sine != Dd::ZERO && result[0] == Dd::ZERO) {
        return Err(bad("acoustic volume phase increment is unrepresentable"));
    }
    Ok(result)
}

fn scaled_product(a: Dd, b: Dd, scale: Dd) -> Dd {
    // The physical positive weight precedes a potentially overflowing product.
    // Scale the larger field for weights <1 to protect a subnormal small field.
    let left = if scale.lt(Dd::ONE) {
        b.abs().lt(a.abs())
    } else {
        a.abs().lt(b.abs())
    };
    if left {
        (a * scale) * b
    } else {
        a * (b * scale)
    }
}

#[allow(clippy::too_many_arguments)]
fn contract(
    reference: &PreparedHarmonicObservation,
    sample: &HarmonicVolumeQuadrature,
    origin: [f64; 3],
    fields: [[Dd; 2]; 4],
    deltas: [[Dd; 2]; 4],
    proper: bool,
    change: bool,
    sums: &mut [Dd; 7],
) -> Result<(), TransferError> {
    let rho = Dd::from_f64(reference.medium.density);
    let c = Dd::from_f64(reference.medium.sound_speed);
    let volume = Dd::from_f64(sample.volume_m3);
    let statistical = if proper { Dd::from_f64(0.5) } else { Dd::ONE };
    let momentum_scale = ((volume / c) / c) * statistical;
    let pressure_scale = (momentum_scale / rho) * Dd::from_f64(0.5);
    let velocity_scale = (volume * (statistical * Dd::from_f64(0.5))) * rho;
    if [momentum_scale, pressure_scale, velocity_scale]
        .iter()
        .any(|&x| !finite(x) || !Dd::ZERO.lt(x))
    {
        return Err(bad(
            "acoustic volume physical storage scale is unrepresentable",
        ));
    }
    let moment = |i: usize, j: usize, scale: Dd| {
        let mut value = Dd::ZERO;
        for part in 0..2 {
            value = value
                + if change {
                    let mid_i = fields[i][part] + deltas[i][part] * Dd::from_f64(0.5);
                    let mid_j = fields[j][part] + deltas[j][part] * Dd::from_f64(0.5);
                    scaled_product(mid_i, deltas[j][part], scale)
                        + scaled_product(deltas[i][part], mid_j, scale)
                } else {
                    scaled_product(fields[i][part], fields[j][part], scale)
                };
        }
        value
    };
    sums[0] = sums[0] + moment(0, 0, pressure_scale);
    let mut momentum = [Dd::ZERO; 3];
    for i in 0..3 {
        sums[0] = sums[0] + moment(i + 1, i + 1, velocity_scale);
        momentum[i] = moment(0, i + 1, momentum_scale);
        sums[i + 1] = sums[i + 1] + momentum[i];
    }
    let position = reference.receivers[sample.receiver].position_m;
    let arm: [Dd; 3] =
        core::array::from_fn(|i| Dd::from_f64(position[i]) - Dd::from_f64(origin[i]));
    for i in 0..3 {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        sums[i + 4] = sums[i + 4] + arm[j] * momentum[k] - arm[k] * momentum[j];
    }
    if sums.iter().any(|x| !finite(*x)) {
        return Err(bad("acoustic volume field contraction is unrepresentable"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
