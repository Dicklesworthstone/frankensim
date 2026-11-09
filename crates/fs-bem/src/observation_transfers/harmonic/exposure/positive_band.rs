//! Shared positive temporal projection of original propagated BEM p/u fields.

use super::*;
use fs_fft::harmonic_window::{
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES, UniformHarmonicBandWorkspace,
};

/// Return false before changing sums when the original pair path is required.
/// Accepted factors share all source-frequency coefficients across receivers
/// and all four physical components. Independent banks are never combined.
pub(super) fn accumulate(
    bank: &Bank,
    window: HarmonicWindow,
    include_mean: bool,
    sums: &mut [[Dd; 14]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, TransferError> {
    let count = bank.lines.len();
    if !(3..=MAX_UNIFORM_HARMONIC_BAND_TONES).contains(&count) {
        return Ok(false);
    }
    let mut frequencies = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    for (frequency, line) in frequencies.iter_mut().zip(&bank.lines) {
        *frequency = line.omega;
    }
    poll(cancelled)?;
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let (prepared, scale) = match bank.statistics {
        HarmonicExposureFactorStatistics::Real => (
            window.prepare_uniform_real_band_rad_s_in(&frequencies[..count], &mut workspace),
            Dd::ONE,
        ),
        HarmonicExposureFactorStatistics::ProperComplex => (
            window.prepare_uniform_band_rad_s_in(&frequencies[..count], &mut workspace),
            Dd::from_f64(0.5),
        ),
    };
    let Some(prepared) = prepared.map_err(window_error)? else {
        return Ok(false);
    };
    poll(cancelled)?;
    let modes = prepared.mode_count();
    let mut fields = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let mut projected = [[[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES]; 4];
    for (receiver, sum) in sums.iter_mut().enumerate() {
        poll(cancelled)?;
        for c in 0..4 {
            for (value, line) in fields.iter_mut().zip(&bank.lines) {
                *value = component(&line.fields[receiver], c);
            }
            prepared
                .project_parts_dd(&fields[..count], &mut projected[c][..modes])
                .map_err(window_error)?;
            poll(cancelled)?;
            if include_mean {
                // Only the deterministic bank has a nonzero ensemble mean.
                // Its real factor's constant mode is the actual window mean.
                sum[c] = sum[c] + projected[c][0][0];
            }
        }
        for mode in 0..modes {
            poll(cancelled)?;
            let mut index = 4;
            for row in 0..4 {
                for col in row..4 {
                    let product = prepared
                        .mode_cross_parts_dd_scaled(
                            mode,
                            projected[row][mode],
                            projected[col][mode],
                            scale,
                        )
                        .map_err(window_error)?;
                    sum[index] = sum[index] + product[0];
                    index += 1;
                }
            }
        }
    }
    poll(cancelled)?;
    Ok(true)
}
