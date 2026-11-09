//! Field-first temporal moments of the already calibrated/delayed channels.

use super::*;
use fs_fft::harmonic_window::{
    HARMONIC_BAND_PREPARATION_SCRATCH_BYTES, HarmonicBandWorkspace,
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES,
};

type Fields = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];

fn eligible(count: usize, statistics: HarmonicExposureFactorStatistics) -> bool {
    let maximum = match statistics {
        HarmonicExposureFactorStatistics::Real => MAX_UNIFORM_HARMONIC_BAND_TONES / 2,
        HarmonicExposureFactorStatistics::ProperComplex => MAX_UNIFORM_HARMONIC_BAND_TONES,
    };
    (2..=maximum).contains(&count)
}

/// Repeated windows are unknown at preparation. Admit the larger of the
/// original pair contraction and the bounded factor preparation/projection,
/// plus its fixed buffers. One factor serves every processed channel.
pub(super) fn admission(
    count: usize,
    channels: usize,
    statistics: HarmonicExposureFactorStatistics,
) -> Result<(usize, usize), TransferError> {
    let pairs = count
        .checked_mul(count)
        .and_then(|p| p.checked_mul(channels))
        .ok_or_else(overflow)?;
    if channels == 0 || !eligible(count, statistics) {
        return Ok((pairs, 0));
    }
    let columns = match statistics {
        HarmonicExposureFactorStatistics::Real => 2 * count,
        HarmonicExposureFactorStatistics::ProperComplex => count,
    };
    let projection = channels
        .checked_mul(count.checked_add(1).ok_or_else(overflow)?)
        .and_then(|p| p.checked_add(columns))
        .and_then(|p| p.checked_mul(MAX_UNIFORM_HARMONIC_BAND_MODES))
        .ok_or_else(overflow)?;
    let scratch = core::mem::size_of::<HarmonicBandWorkspace>()
        + HARMONIC_BAND_PREPARATION_SCRATCH_BYTES
        + core::mem::size_of::<[f64; MAX_UNIFORM_HARMONIC_BAND_TONES]>()
        + core::mem::size_of::<Fields>()
        + core::mem::size_of::<Modes>();
    Ok((pairs.max(projection), scratch))
}

/// False leaves both outputs unchanged for the exact original pair fallback.
/// Actual-real mode zero gives the profile-weighted field mean; proper mode
/// zero must not be used as an absolute-carrier mean or readout amplitude.
pub(super) fn accumulate(
    bank: &PressureBank,
    window: HarmonicWindow,
    sums: &mut [[Dd; 2]],
    integrated: &mut [[Dd; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, TransferError> {
    let count = bank.lines.len();
    if !eligible(count, bank.statistics) {
        return Ok(false);
    }
    poll(cancelled)?;
    let mut frequencies = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    for (out, line) in frequencies.iter_mut().zip(&bank.lines) {
        *out = line.omega;
    }
    let mut workspace = HarmonicBandWorkspace::new();
    let (factor, scale) = match bank.statistics {
        HarmonicExposureFactorStatistics::Real => (
            window.prepare_real_band_rad_s_in(&frequencies[..count], &mut workspace),
            Dd::ONE,
        ),
        HarmonicExposureFactorStatistics::ProperComplex => (
            window.prepare_band_rad_s_in(&frequencies[..count], &mut workspace),
            Dd::from_f64(0.5),
        ),
    };
    let Some(factor) = factor.map_err(window_error)? else {
        return Ok(false);
    };
    let modes = factor.mode_count();
    let mut fields: Fields = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let mut projected: Modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    for (channel, (sum, mean)) in sums.iter_mut().zip(integrated).enumerate() {
        poll(cancelled)?;
        for (field, line) in fields.iter_mut().zip(&bank.lines) {
            *field = line.amplitudes[channel];
        }
        factor
            .project_parts_dd(&fields[..count], &mut projected[..modes])
            .map_err(window_error)?;
        if bank.statistics == HarmonicExposureFactorStatistics::Real {
            *mean = [projected[0][0], Dd::ZERO];
        }
        for (mode, &value) in projected[..modes].iter().enumerate() {
            poll(cancelled)?;
            sum[1] = sum[1]
                + factor
                    .mode_cross_parts_dd_scaled(mode, value, value, scale)
                    .map_err(window_error)?[0];
        }
    }
    poll(cancelled)?;
    Ok(true)
}
