//! Shared positive temporal moments for native duct pressure/volume-flow banks.
//!
//! The original source-pair budgets remain structural input bounds. This
//! optional factor has its own hard cap of 64 tones and 17 temporal modes;
//! no temporal heap storage is added to the port or interior-field workspace.

use super::{Dd, DuctError, HarmonicWindow, poll, window_error};
use fs_fft::harmonic_window::{
    HarmonicBandWorkspace, MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES,
};

/// Return false when the original analytic pair kernels must be used.
/// Pressure, flow and signed work share one factor, preserving their Gram
/// relation. Proper factors use the physical real-peak half before products.
pub(in crate::network) fn accumulate(
    window: HarmonicWindow,
    frequencies: impl ExactSizeIterator<Item = f64>,
    fields: impl Fn(usize, usize) -> [Dd; 4],
    proper: bool,
    deterministic: bool,
    sums: &mut [[Dd; 5]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bool, DuctError> {
    let count = frequencies.len();
    let maximum = if proper {
        MAX_UNIFORM_HARMONIC_BAND_TONES
    } else {
        MAX_UNIFORM_HARMONIC_BAND_TONES / 2
    };
    if !(3..=maximum).contains(&count) {
        return Ok(false);
    }
    poll(cancelled, 0)?;
    let mut omega = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    for (slot, frequency) in omega.iter_mut().zip(frequencies) {
        *slot = frequency;
    }
    let mut workspace = HarmonicBandWorkspace::new();
    let band = if proper {
        window.prepare_band_rad_s_in(&omega[..count], &mut workspace)
    } else {
        window.prepare_real_band_rad_s_in(&omega[..count], &mut workspace)
    }
    .map_err(window_error)?;
    let Some(band) = band else { return Ok(false) };
    let modes = band.mode_count();
    let scale = if proper { Dd::from_f64(0.5) } else { Dd::ONE };
    let mut pressure = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let mut flow = pressure;
    let mut p = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    let mut u = p;
    for (probe, sum) in sums.iter_mut().enumerate() {
        poll(cancelled, probe)?;
        for i in 0..count {
            let f = fields(i, probe);
            pressure[i] = [f[0], f[1]];
            flow[i] = [f[2], f[3]];
        }
        band.project_parts_dd(&pressure[..count], &mut p[..modes])
            .map_err(window_error)?;
        poll(cancelled, probe)?;
        band.project_parts_dd(&flow[..count], &mut u[..modes])
            .map_err(window_error)?;
        if deterministic {
            // The real constant orthogonal coefficient is the weighted window
            // mean. No rounded carrier phases can swamp a cancelled mean.
            sum[0] = sum[0] + p[0][0];
            sum[1] = sum[1] + u[0][0];
        }
        for mode in 0..modes {
            poll(cancelled, probe)?;
            let pp = band
                .mode_cross_parts_dd_scaled(mode, p[mode], p[mode], scale)
                .map_err(window_error)?;
            let uu = band
                .mode_cross_parts_dd_scaled(mode, u[mode], u[mode], scale)
                .map_err(window_error)?;
            let pu = band
                .mode_cross_parts_dd_scaled(mode, p[mode], u[mode], scale)
                .map_err(window_error)?;
            sum[2] = sum[2] + pp[0];
            sum[3] = sum[3] + uu[0];
            sum[4] = sum[4] + pu[0];
        }
    }
    poll(cancelled, sums.len())?;
    Ok(true)
}

#[cfg(test)]
mod tests;
