//! Shared temporal factors for selected, signed-lag channel products.

use super::*;
use fs_fft::harmonic_window::{
    HARMONIC_BAND_PREPARATION_SCRATCH_BYTES, HarmonicBandWorkspace,
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES,
};

type Fields = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
type Modes = [[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
type FirstModes = [Modes; MAX_PRESSURE_ARRAY_CHANNELS];

fn eligible(count: usize, statistics: HarmonicExposureFactorStatistics) -> bool {
    let maximum = match statistics {
        HarmonicExposureFactorStatistics::Real => MAX_UNIFORM_HARMONIC_BAND_TONES / 2,
        HarmonicExposureFactorStatistics::ProperComplex => MAX_UNIFORM_HARMONIC_BAND_TONES,
    };
    (2..=maximum).contains(&count)
}

pub(super) fn admission(
    count: usize,
    pairs: usize,
    first_channels: usize,
    statistics: HarmonicExposureFactorStatistics,
) -> Result<(usize, usize), TransferError> {
    let original = count
        .checked_mul(count)
        .and_then(|p| p.checked_mul(pairs))
        .ok_or_else(overflow)?;
    if pairs == 0 || !eligible(count, statistics) {
        return Ok((original, 0));
    }
    let columns = match statistics {
        HarmonicExposureFactorStatistics::Real => 2 * count,
        HarmonicExposureFactorStatistics::ProperComplex => count,
    };
    let projection = first_channels
        .checked_add(pairs)
        .and_then(|p| p.checked_mul(count))
        .and_then(|p| p.checked_add(pairs))
        .and_then(|p| p.checked_add(columns))
        .and_then(|p| p.checked_mul(MAX_UNIFORM_HARMONIC_BAND_MODES))
        .ok_or_else(overflow)?;
    let scratch = core::mem::size_of::<HarmonicBandWorkspace>()
        + HARMONIC_BAND_PREPARATION_SCRATCH_BYTES
        + core::mem::size_of::<[f64; MAX_UNIFORM_HARMONIC_BAND_TONES]>()
        + core::mem::size_of::<Fields>()
        + core::mem::size_of::<Modes>()
        + core::mem::size_of::<FirstModes>()
        + core::mem::size_of::<[bool; MAX_PRESSURE_ARRAY_CHANNELS]>();
    Ok((original.max(projection), scratch))
}

/// False preserves both outputs for the original exact pair fallback.
pub(super) fn accumulate(
    bank: &PressureBank,
    lagged: &CorrelationBank,
    pairs: &[HarmonicPressureCorrelationPair],
    window: HarmonicWindow,
    fluctuation: bool,
    sums: &mut [[Dd; 5]],
    integrated: &mut [[[Dd; 2]; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, TransferError> {
    let count = bank.lines.len();
    if pairs.is_empty() || !eligible(count, bank.statistics) {
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
    let mut first: FirstModes =
        [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; MAX_PRESSURE_ARRAY_CHANNELS];
    let mut cached = [false; MAX_PRESSURE_ARRAY_CHANNELS];
    let mut second: Modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    for (p, ((pair, sum), mean)) in pairs.iter().zip(sums).zip(integrated).enumerate() {
        poll(cancelled)?;
        if !cached[pair.first] {
            for (field, line) in fields.iter_mut().zip(&bank.lines) {
                *field = line.amplitudes[pair.first];
            }
            factor
                .project_parts_dd(&fields[..count], &mut first[pair.first][..modes])
                .map_err(window_error)?;
            cached[pair.first] = true;
        }
        for (line, field) in fields[..count].iter_mut().enumerate() {
            *field = lagged.shifted[line * pairs.len() + p];
        }
        factor
            .project_parts_dd(&fields[..count], &mut second[..modes])
            .map_err(window_error)?;
        if bank.statistics == HarmonicExposureFactorStatistics::Real {
            *mean = [
                [first[pair.first][0][0], Dd::ZERO],
                [second[0][0], Dd::ZERO],
            ];
        }
        let mut product = Dd::ZERO;
        for (mode, (&a, &b)) in first[pair.first][..modes]
            .iter()
            .zip(&second[..modes])
            .enumerate()
        {
            poll(cancelled)?;
            product = product
                + factor
                    .mode_cross_parts_dd_scaled(mode, a, b, scale)
                    .map_err(window_error)?[0];
        }
        sum[2] = sum[2] + product;
        if fluctuation {
            sum[3] = sum[3] + product;
        }
    }
    poll(cancelled)?;
    Ok(true)
}
