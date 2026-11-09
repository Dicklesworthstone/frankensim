//! Positive temporal field factors for unresolved, uniform-window bands.
//!
//! One factorization serves the pressure, every velocity component and their
//! signed cross moments. Real waveforms keep their sum-frequency terms; proper
//! complex factors use half the analytic-field outer product.

use super::{
    AeroacError, Dd, FactorKind, FieldMoments, HarmonicWindow, PreparedHarmonicField, poll,
    window_error,
};
use fs_fft::harmonic_window::{
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES, UniformHarmonicBandRef,
    UniformHarmonicBandWorkspace,
};

type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];

pub(super) fn accumulate<const D: usize>(
    field: &PreparedHarmonicField<D>,
    window: HarmonicWindow,
    kind: FactorKind,
    include_mean: bool,
    full_velocity: bool,
    sums: &mut FieldMoments<D>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, AeroacError> {
    let count = field.lines.len();
    if !(3..=MAX_UNIFORM_HARMONIC_BAND_TONES).contains(&count) {
        return Ok(false);
    }
    let angular = field.lines[0].frequency.is_angular();
    let mut frequencies = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    for (out, line) in frequencies.iter_mut().zip(&field.lines) {
        poll(cancelled)?;
        // Keep the original source-clock unit through the shared numerical
        // owner. A mixed internal bank retains its existing DD pair kernels.
        if line.frequency.is_angular() != angular {
            return Ok(false);
        }
        *out = line.frequency.native();
    }
    let frequencies = &frequencies[..count];
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let band = match (kind, angular) {
        (FactorKind::Deterministic, false) => {
            window.prepare_uniform_real_band_hz_in(frequencies, &mut workspace)
        }
        (FactorKind::Deterministic, true) => {
            window.prepare_uniform_real_band_rad_s_in(frequencies, &mut workspace)
        }
        (FactorKind::Proper, false) => {
            window.prepare_uniform_band_hz_in(frequencies, &mut workspace)
        }
        (FactorKind::Proper, true) => {
            window.prepare_uniform_band_rad_s_in(frequencies, &mut workspace)
        }
    }
    .map_err(window_error)?;
    poll(cancelled)?;
    let Some(band) = band else {
        return Ok(false);
    };
    let pressure = project(field, None, &band, cancelled)?;
    let mut velocity = [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; D];
    for (axis, modes) in velocity.iter_mut().enumerate() {
        *modes = project(field, Some(axis), &band, cancelled)?;
    }
    if include_mean {
        // Only deterministic banks include means. Their real-field plan keeps
        // the absolute carrier, and its constant mode is the actual time mean.
        sums.pressure_mean = sums.pressure_mean + pressure[0][0];
        for (sum, modes) in sums.velocity_mean.iter_mut().zip(&velocity) {
            *sum = *sum + modes[0][0];
        }
    }
    let normalization = match kind {
        FactorKind::Deterministic => Dd::ONE,
        FactorKind::Proper => Dd::from_f64(0.5),
    };
    for mode in 0..band.mode_count() {
        poll(cancelled)?;
        let product = |a, b| {
            band.mode_cross_parts_dd_scaled(mode, a, b, normalization)
                .map(|value| value[0])
                .map_err(window_error)
        };
        sums.pressure_square = sums.pressure_square + product(pressure[mode], pressure[mode])?;
        for (axis, sum) in sums.intensity.iter_mut().enumerate() {
            *sum = *sum + product(pressure[mode], velocity[axis][mode])?;
        }
        if full_velocity {
            for (row, values) in sums.velocity.iter_mut().enumerate() {
                poll(cancelled)?;
                for (column, sum) in values.iter_mut().enumerate().skip(row) {
                    *sum = *sum + product(velocity[row][mode], velocity[column][mode])?;
                }
            }
        }
    }
    poll(cancelled)?;
    Ok(true)
}

fn project<const D: usize>(
    field: &PreparedHarmonicField<D>,
    axis: Option<usize>,
    band: &UniformHarmonicBandRef<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Modes, AeroacError> {
    let mut components = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
    for (out, line) in components.iter_mut().zip(&field.lines) {
        poll(cancelled)?;
        *out = match axis {
            Some(axis) => line.velocity[axis].parts(),
            None => line.pressure.parts(),
        };
    }
    let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    band.project_parts_dd(
        &components[..field.lines.len()],
        &mut modes[..band.mode_count()],
    )
    .map_err(window_error)?;
    poll(cancelled)?;
    Ok(modes)
}
