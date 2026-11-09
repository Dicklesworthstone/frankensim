//! Original-field surface power and wave-momentum transport before tone reduction.
//!
//! Short uniform and triangular windows share one factor per source bank. The original
//! pressure and three velocity components interfere before any product; signed
//! surface power and the full symmetric wave stress use those same projections.
//! Other banks retain their original real/proper pair law independently.

use super::*;
use fs_fft::harmonic_window::{
    HarmonicBandWorkspace, MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES,
};

/// Explicit per-query work admission. Queries allocate no heap memory and add
/// only fixed shared temporal-factor scratch and four 17-mode field vectors.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSurfaceFieldLimits {
    /// Authored entries, including repetitions and zero area vectors.
    pub maximum_quadrature_points: usize,
    /// Four field projections plus ten moment products per receiver and mode:
    /// sum_b M_b*(F_b + Q*(4*F_b + 10)) over admitted banks. Power-only queries
    /// use the same conservative count. Factor preparation is included.
    pub maximum_temporal_projection_terms: usize,
    /// Q*sum_b(F_b²) for banks using the original temporal pair kernels.
    /// Each pair has at most ten p/u products and seven signed output channels.
    pub maximum_quadrature_frequency_pairs: usize,
}

/// Expected finite-window power and wave-pseudomomentum flux on the authored
/// surface. All fields retain double-double low parts and physical signs.
/// These rates do not identify heat, a body's force, or total fluid momentum.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicSurfaceFieldMoments {
    pub mean_power_w: Dd,
    pub mean_linear_momentum_rate_n: [Dd; 3],
    pub mean_angular_momentum_rate_nm: [Dd; 3],
}

/// Signed transported energy and wave momentum over one uniform source interval.
/// No enclosure, matching volume storage, source term or body reaction is inferred.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicSurfaceFieldTransfer {
    pub energy_j: Dd,
    pub wave_momentum: HarmonicWaveMomentumTransfer,
}

impl PreparedHarmonicExposure<'_> {
    /// Signed original-field acoustic power [W] across the supplied oriented
    /// quadrature. This p*u observable also admits the original relaxing-medium
    /// owners; it does not apply the nondispersive wave-stress law to them.
    pub fn surface_field_power_window(
        &self,
        quadrature: &[HarmonicSurfaceQuadrature],
        window: HarmonicWindow,
        limits: HarmonicSurfaceFieldLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Dd, TransferError> {
        Ok(self.surface_fields(quadrature, None, window, Dd::ONE, limits, &mut cancelled)?[0])
    }

    /// Signed original-field acoustic energy [J], retaining a DD source epoch
    /// independently of the positive finite duration. No carrier-cycle average
    /// or normalized detector weight replaces the physical uniform interval.
    pub fn surface_field_energy_from_clock_dd(
        &self,
        quadrature: &[HarmonicSurfaceQuadrature],
        start_s: Dd,
        duration_s: f64,
        limits: HarmonicSurfaceFieldLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Dd, TransferError> {
        let window = HarmonicWindow::from_clock_dd(start_s, duration_s).map_err(window_error)?;
        Ok(self.surface_fields(
            quadrature,
            None,
            window,
            Dd::from_f64(duration_s),
            limits,
            &mut cancelled,
        )?[0])
    }

    /// Power [W] and signed linear/angular wave-momentum rates [N, N m].
    /// Every original frequency owner must authorize the homogeneous inviscid
    /// nondispersive rest-medium stress S=(p²/(rho*c²)-rho*u.u)I/2+rho*u*u.
    /// The exact receiver positions and supplied origin define moment arms.
    pub fn surface_field_window(
        &self,
        quadrature: &[HarmonicSurfaceQuadrature],
        moment_origin_m: [f64; 3],
        window: HarmonicWindow,
        limits: HarmonicSurfaceFieldLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicSurfaceFieldMoments, TransferError> {
        let values = self.surface_fields(
            quadrature,
            Some(moment_origin_m),
            window,
            Dd::ONE,
            limits,
            &mut cancelled,
        )?;
        Ok(HarmonicSurfaceFieldMoments {
            mean_power_w: values[0],
            mean_linear_momentum_rate_n: [values[1], values[2], values[3]],
            mean_angular_momentum_rate_nm: [values[4], values[5], values[6]],
        })
    }

    /// Signed energy [J] and wave momentum [N s, N m s] over the actual uniform
    /// source interval. The duration enters temporal products before reduction.
    /// Refusal/cancellation leaves the immutable original exposure reusable.
    #[allow(clippy::too_many_arguments)]
    pub fn surface_field_interval_from_clock_dd(
        &self,
        quadrature: &[HarmonicSurfaceQuadrature],
        moment_origin_m: [f64; 3],
        start_s: Dd,
        duration_s: f64,
        limits: HarmonicSurfaceFieldLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicSurfaceFieldTransfer, TransferError> {
        let window = HarmonicWindow::from_clock_dd(start_s, duration_s).map_err(window_error)?;
        let values = self.surface_fields(
            quadrature,
            Some(moment_origin_m),
            window,
            Dd::from_f64(duration_s),
            limits,
            &mut cancelled,
        )?;
        Ok(HarmonicSurfaceFieldTransfer {
            energy_j: values[0],
            wave_momentum: HarmonicWaveMomentumTransfer {
                linear_momentum_ns: [values[1], values[2], values[3]],
                angular_momentum_nms: [values[4], values[5], values[6]],
            },
        })
    }

    fn surface_fields(
        &self,
        quadrature: &[HarmonicSurfaceQuadrature],
        origin: Option<[f64; 3]>,
        window: HarmonicWindow,
        scale: Dd,
        limits: HarmonicSurfaceFieldLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 7], TransferError> {
        poll(cancelled)?;
        admit(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        admit(quadrature.len(), limits.maximum_quadrature_points)?;
        if origin.is_some_and(|o| !self.wave_momentum || o.iter().any(|x| !x.is_finite())) {
            return Err(bad(
                "original acoustic owners do not authorize this wave-stress query",
            ));
        }
        for sample in quadrature {
            poll(cancelled)?;
            if sample.receiver >= self.reference.receivers.len()
                || sample.oriented_area_m2.iter().any(|x| !x.is_finite())
                || (origin.is_some()
                    && self.reference.receivers[sample.receiver]
                        .position_m
                        .iter()
                        .any(|x| !x.is_finite()))
            {
                return Err(bad("original acoustic surface quadrature is invalid"));
            }
        }
        let banks = || core::iter::once(&self.mean).chain(self.factors.iter());
        let mut projection_terms = 0usize;
        let mut fallback_pairs = 0usize;
        if !quadrature.is_empty() {
            for bank in banks() {
                poll(cancelled)?;
                let f = bank.lines.len();
                if let Some(modes) = mode_count(bank, window)? {
                    projection_terms = f
                        .checked_mul(4)
                        .and_then(|n| n.checked_add(10))
                        .and_then(|n| n.checked_mul(quadrature.len()))
                        .and_then(|n| n.checked_add(f))
                        .and_then(|n| n.checked_mul(modes))
                        .and_then(|n| n.checked_add(projection_terms))
                        .ok_or_else(overflow)?;
                } else {
                    fallback_pairs = f
                        .checked_mul(f)
                        .and_then(|n| n.checked_mul(quadrature.len()))
                        .and_then(|n| n.checked_add(fallback_pairs))
                        .ok_or_else(overflow)?;
                }
            }
        }
        admit(projection_terms, limits.maximum_temporal_projection_terms)?;
        admit(fallback_pairs, MAX_PRODUCTS)?;
        admit(fallback_pairs, limits.maximum_quadrature_frequency_pairs)?;
        let mut result = [Dd::ZERO; 7];
        if !quadrature.is_empty() {
            for bank in banks() {
                poll(cancelled)?;
                if projected_bank(
                    bank,
                    self.reference,
                    quadrature,
                    origin,
                    window,
                    scale,
                    &mut result,
                    cancelled,
                )? {
                    continue;
                }
                for a in &bank.lines {
                    for b in &bank.lines {
                        poll(cancelled)?;
                        let kernel = bank.statistics.kernel(window, a.omega, b.omega)?;
                        for sample in quadrature {
                            poll(cancelled)?;
                            if zero_area(sample) {
                                continue;
                            }
                            let af = &a.fields[sample.receiver];
                            let bf = &b.fields[sample.receiver];
                            contract(self.reference, sample, origin, &mut result, |i, j| {
                                Ok(kernel.product(component(af, i), component(bf, j))? * scale)
                            })?;
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        if result.iter().any(|&x| !finite(x)) {
            return Err(bad("acoustic surface field exchange is unrepresentable"));
        }
        Ok(result)
    }
}

fn frequencies(bank: &Bank) -> [f64; MAX_UNIFORM_HARMONIC_BAND_TONES] {
    let mut result = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    for (value, line) in result.iter_mut().zip(&bank.lines) {
        *value = line.omega;
    }
    result
}
fn mode_count(bank: &Bank, window: HarmonicWindow) -> Result<Option<usize>, TransferError> {
    if !(3..=MAX_UNIFORM_HARMONIC_BAND_TONES).contains(&bank.lines.len()) {
        return Ok(None);
    }
    let values = frequencies(bank);
    match bank.statistics {
        HarmonicExposureFactorStatistics::Real => {
            window.real_band_mode_count_rad_s(&values[..bank.lines.len()])
        }
        HarmonicExposureFactorStatistics::ProperComplex => {
            window.band_mode_count_rad_s(&values[..bank.lines.len()])
        }
    }
    .map_err(window_error)
}
fn zero_area(sample: &HarmonicSurfaceQuadrature) -> bool {
    sample.oriented_area_m2.iter().all(|&x| x == 0.0)
}

#[allow(clippy::too_many_arguments)]
fn projected_bank(
    bank: &Bank,
    reference: &PreparedHarmonicObservation,
    quadrature: &[HarmonicSurfaceQuadrature],
    origin: Option<[f64; 3]>,
    window: HarmonicWindow,
    scale: Dd,
    sums: &mut [Dd; 7],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, TransferError> {
    if mode_count(bank, window)?.is_none() {
        return Ok(false);
    }
    let values = frequencies(bank);
    let count = bank.lines.len();
    let mut workspace = HarmonicBandWorkspace::new();
    let (factor, physical_scale) = match bank.statistics {
        HarmonicExposureFactorStatistics::Real => (
            window.prepare_real_band_rad_s_in(&values[..count], &mut workspace),
            scale,
        ),
        HarmonicExposureFactorStatistics::ProperComplex => (
            window.prepare_band_rad_s_in(&values[..count], &mut workspace),
            scale * Dd::from_f64(0.5),
        ),
    };
    // A positive interval must not silently disappear at the proper half.
    // Large fields could restore finite energy after that scale underflows.
    if physical_scale == Dd::ZERO && scale != Dd::ZERO {
        return Err(bad(
            "acoustic surface physical time scale is unrepresentable",
        ));
    }
    let factor = factor
        .map_err(window_error)?
        .ok_or_else(|| bad("admitted temporal factor absent"))?;
    let modes = factor.mode_count();
    let mut fields = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let mut projected = [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 4];
    for sample in quadrature {
        poll(cancelled)?;
        if zero_area(sample) {
            continue;
        }
        for c in 0..4 {
            for (field, line) in fields.iter_mut().zip(&bank.lines) {
                *field = component(&line.fields[sample.receiver], c);
            }
            factor
                .project_parts_dd(&fields[..count], &mut projected[c][..modes])
                .map_err(window_error)?;
            poll(cancelled)?;
        }
        for mode in 0..modes {
            poll(cancelled)?;
            contract(reference, sample, origin, sums, |i, j| {
                Ok(factor
                    .mode_cross_parts_dd_scaled(
                        mode,
                        projected[i][mode],
                        projected[j][mode],
                        physical_scale,
                    )
                    .map_err(window_error)?[0])
            })?;
        }
    }
    poll(cancelled)?;
    Ok(true)
}

fn contract(
    reference: &PreparedHarmonicObservation,
    sample: &HarmonicSurfaceQuadrature,
    origin: Option<[f64; 3]>,
    sums: &mut [Dd; 7],
    mut moment: impl FnMut(usize, usize) -> Result<Dd, TransferError>,
) -> Result<(), TransferError> {
    let area = sample.oriented_area_m2.map(Dd::from_f64);
    for axis in 0..3 {
        sums[0] = sums[0] + area[axis] * moment(0, axis + 1)?;
    }
    if let Some(origin) = origin {
        let rho = Dd::from_f64(reference.medium.density);
        let c = Dd::from_f64(reference.medium.sound_speed);
        let mut tensor = [[Dd::ZERO; 3]; 3];
        let mut isotropic = ((moment(0, 0)? / c) / c) / rho;
        for i in 0..3 {
            for j in i..3 {
                tensor[i][j] = rho * moment(i + 1, j + 1)?;
                tensor[j][i] = tensor[i][j];
            }
            isotropic = isotropic - tensor[i][i];
        }
        isotropic = isotropic * Dd::from_f64(0.5);
        let mut traction = [Dd::ZERO; 3];
        for i in 0..3 {
            traction[i] = isotropic * area[i];
            for j in 0..3 {
                traction[i] = traction[i] + tensor[i][j] * area[j];
            }
            sums[i + 1] = sums[i + 1] + traction[i];
        }
        let position = reference.receivers[sample.receiver].position_m;
        let arm: [Dd; 3] =
            core::array::from_fn(|i| Dd::from_f64(position[i]) - Dd::from_f64(origin[i]));
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            sums[i + 4] = sums[i + 4] + arm[j] * traction[k] - arm[k] * traction[j];
        }
    }
    if sums.iter().any(|&x| !finite(x)) {
        return Err(bad("acoustic surface field contraction is unrepresentable"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
